from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, BinaryIO, Iterable

from . import __version__
from .baseline import _file_sha256
from .gff import parse_attributes


CONSENSUS_SCHEMA_VERSION = "ploidypatch.method_consensus.v1"
CHAIN_PRESERVING_SCHEMA_VERSION = "ploidypatch.method_candidate_pool.v2"
SAFE_METHOD = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
REDUNDANCY_INTERVAL_BIN_BP = 1_000_000
REDUNDANCY_POLICIES = frozenset(
    {"suppress_overlapping", "retain_distinct_chains"}
)


@dataclass(frozen=True)
class ConsensusInputModel:
    method: str
    upstream_model: str
    seqid: str
    strand: str
    cds: tuple[tuple[int, int, str], ...]

    @property
    def start(self) -> int:
        return min(start for start, _, _ in self.cds)

    @property
    def end(self) -> int:
        return max(end for _, end, _ in self.cds)

    @property
    def cds_intervals(self) -> tuple[tuple[int, int], ...]:
        return _merge_intervals((start, end) for start, end, _ in self.cds)

    @property
    def cds_bp(self) -> int:
        return sum(end - start + 1 for start, end in self.cds_intervals)

    @property
    def chain_key(self) -> tuple[str, str, tuple[tuple[int, int, str], ...]]:
        return self.seqid, self.strand, self.cds


def _merge_intervals(
    intervals: Iterable[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    merged: list[list[int]] = []
    for start, end in sorted(set(intervals)):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return tuple((start, end) for start, end in merged)


def _overlap_bp(
    left: tuple[tuple[int, int], ...],
    right: tuple[tuple[int, int], ...],
) -> int:
    left_index = 0
    right_index = 0
    overlap = 0
    while left_index < len(left) and right_index < len(right):
        left_start, left_end = left[left_index]
        right_start, right_end = right[right_index]
        overlap += max(0, min(left_end, right_end) - max(left_start, right_start) + 1)
        if left_end <= right_end:
            left_index += 1
        else:
            right_index += 1
    return overlap


def _consume_candidate_prefix(
    base_path: Path, candidate_handle: BinaryIO, candidate_path: Path
) -> None:
    last_byte = b""
    with base_path.open("rb") as base_handle:
        while chunk := base_handle.read(8 * 1024 * 1024):
            observed = candidate_handle.read(len(chunk))
            if observed != chunk:
                raise ValueError(
                    f"Candidate does not preserve the exact base GFF prefix: {candidate_path}"
                )
            last_byte = chunk[-1:]
    if last_byte not in {b"\n", b"\r"}:
        if candidate_handle.read(1) != b"\n":
            raise ValueError(f"Candidate lacks separator after base GFF: {candidate_path}")
    marker = candidate_handle.readline()
    if marker not in {b"###\n", b"###\r\n"}:
        raise ValueError(f"Candidate lacks adapter boundary marker: {candidate_path}")


def _read_candidate_models(
    method: str, candidate_path: Path, base_path: Path
) -> list[ConsensusInputModel]:
    transcripts: dict[str, tuple[str, str, int, int, dict[str, str]]] = {}
    children: dict[
        str, list[tuple[str, int, int, str, str, str]]
    ] = defaultdict(list)
    with candidate_path.open("rb") as raw_handle:
        _consume_candidate_prefix(base_path, raw_handle, candidate_path)
        text_handle = io.TextIOWrapper(
            raw_handle, encoding="utf-8", errors="strict", newline=None
        )
        try:
            for suffix_line, raw_line in enumerate(text_handle, start=1):
                stripped = raw_line.rstrip("\r\n")
                if not stripped or stripped.startswith("#"):
                    continue
                fields = stripped.split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"Malformed appended GFF line {suffix_line}: {candidate_path}"
                    )
                if fields[1] != "PloidyPatchBaseline":
                    raise ValueError(
                        f"Non-baseline feature after adapter marker: {candidate_path}"
                    )
                try:
                    start = int(fields[3])
                    end = int(fields[4])
                except ValueError as exc:
                    raise ValueError(
                        f"Non-integer appended coordinate: {candidate_path}"
                    ) from exc
                attributes, malformed = parse_attributes(fields[8])
                if malformed:
                    raise ValueError(f"Malformed appended attributes: {candidate_path}")
                feature_type = fields[2]
                if feature_type in {"mRNA", "transcript"}:
                    transcript_id = attributes.get("ID", "")
                    if not transcript_id or transcript_id in transcripts:
                        raise ValueError(
                            f"Empty or duplicate appended transcript ID: {candidate_path}"
                        )
                    transcripts[transcript_id] = (
                        fields[0],
                        fields[6],
                        start,
                        end,
                        attributes,
                    )
                elif feature_type in {"CDS", "exon"}:
                    parents = tuple(
                        value
                        for value in attributes.get("Parent", "").split(",")
                        if value
                    )
                    if len(parents) != 1:
                        raise ValueError(
                            f"Appended child must have one parent: {candidate_path}"
                        )
                    children[parents[0]].append(
                        (
                            feature_type,
                            start,
                            end,
                            fields[7],
                            fields[0],
                            fields[6],
                        )
                    )
                elif feature_type != "gene":
                    raise ValueError(
                        f"Unsupported appended feature type {feature_type}: {candidate_path}"
                    )
        finally:
            text_handle.detach()

    orphan_children = set(children) - set(transcripts)
    if orphan_children:
        raise ValueError(
            f"Appended children reference unknown transcripts: {candidate_path}"
        )
    models: list[ConsensusInputModel] = []
    observed_keys: set[tuple[str, str, tuple[tuple[int, int, str], ...]]] = set()
    for transcript_id, (seqid, strand, tx_start, tx_end, attributes) in transcripts.items():
        if strand not in {"+", "-"}:
            raise ValueError(f"Invalid appended transcript strand: {transcript_id}")
        cds = tuple(
            sorted(
                (start, end, phase)
                for child_type, start, end, phase, child_seqid, child_strand in children.get(
                    transcript_id, []
                )
                if child_type == "CDS"
                and child_seqid == seqid
                and child_strand == strand
            )
        )
        if not cds or any(phase not in {"0", "1", "2"} for _, _, phase in cds):
            raise ValueError(f"Appended transcript lacks a phased CDS: {transcript_id}")
        if any(start < tx_start or end > tx_end for start, end, _ in cds):
            raise ValueError(f"Appended CDS is outside transcript: {transcript_id}")
        raw_bp = sum(end - start + 1 for start, end, _ in cds)
        if raw_bp != sum(
            end - start + 1
            for start, end in _merge_intervals((start, end) for start, end, _ in cds)
        ):
            raise ValueError(f"Appended transcript has overlapping CDS segments: {transcript_id}")
        upstream_model = (
            attributes.get("upstream_model")
            or attributes.get("miniprot_model")
            or attributes.get("Derives_from")
            or transcript_id
        )
        model = ConsensusInputModel(
            method=method,
            upstream_model=upstream_model,
            seqid=seqid,
            strand=strand,
            cds=cds,
        )
        if model.chain_key in observed_keys:
            raise ValueError(
                f"Method {method} contains duplicate exact CDS chains: {transcript_id}"
            )
        observed_keys.add(model.chain_key)
        models.append(model)
    return models


def _chain_digest(
    key: tuple[str, str, tuple[tuple[int, int, str], ...]]
) -> str:
    payload = json.dumps(key, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _conflict_sets(
    eligible: list[
        tuple[
            tuple[str, str, tuple[tuple[int, int, str], ...]],
            list[ConsensusInputModel],
        ]
    ],
    *,
    max_redundancy_overlap: float,
) -> tuple[dict[str, tuple[str, int]], int]:
    """Return stable connected components of strongly overlapping CDS chains."""

    parents = list(range(len(eligible)))

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[max(left_root, right_root)] = min(left_root, right_root)

    interval_bins: dict[tuple[str, str, int], list[int]] = defaultdict(list)
    for index, (_, models) in enumerate(eligible):
        representative = models[0]
        first_bin = (representative.start - 1) // REDUNDANCY_INTERVAL_BIN_BP
        last_bin = (representative.end - 1) // REDUNDANCY_INTERVAL_BIN_BP
        nearby_indices: set[int] = set()
        for bin_index in range(first_bin, last_bin + 1):
            nearby_indices.update(
                interval_bins.get(
                    (representative.seqid, representative.strand, bin_index), ()
                )
            )
        for other_index in sorted(nearby_indices):
            other = eligible[other_index][1][0]
            overlap = _overlap_bp(
                representative.cds_intervals, other.cds_intervals
            )
            denominator = min(representative.cds_bp, other.cds_bp)
            if denominator and overlap / denominator > max_redundancy_overlap:
                union(index, other_index)
        for bin_index in range(first_bin, last_bin + 1):
            interval_bins[
                (representative.seqid, representative.strand, bin_index)
            ].append(index)

    members_by_root: dict[int, list[int]] = defaultdict(list)
    for index in range(len(eligible)):
        members_by_root[find(index)].append(index)
    annotations: dict[str, tuple[str, int]] = {}
    conflict_set_count = 0
    for members in members_by_root.values():
        if len(members) < 2:
            continue
        member_digests = sorted(_chain_digest(eligible[index][0]) for index in members)
        conflict_digest = hashlib.sha256(
            ("ploidypatch-conflict-v1:" + ",".join(member_digests)).encode("ascii")
        ).hexdigest()
        for digest in member_digests:
            annotations[digest] = (conflict_digest, len(member_digests))
        conflict_set_count += 1
    return annotations, conflict_set_count


def select_method_consensus(
    *,
    base_gff_path: str | Path,
    candidate_inputs: Iterable[tuple[str, str | Path]],
    output_gff_path: str | Path,
    decisions_tsv_path: str | Path,
    min_method_support: int = 2,
    max_redundancy_overlap: float = 0.5,
    redundancy_policy: str = "suppress_overlapping",
) -> dict[str, Any]:
    """Append exact phased-CDS chains supported by independent method families."""

    inputs = [(method, Path(path)) for method, path in candidate_inputs]
    labels = [method for method, _ in inputs]
    if len(inputs) < 2:
        raise ValueError("At least two candidate methods are required")
    if any(not SAFE_METHOD.fullmatch(label) for label in labels):
        raise ValueError("Method labels must be shell-safe identifiers")
    if len(set(labels)) != len(labels):
        raise ValueError("Method labels must be unique")
    if not 1 <= min_method_support <= len(inputs):
        raise ValueError("min_method_support must be within the method count")
    if not 0 <= max_redundancy_overlap <= 1:
        raise ValueError("max_redundancy_overlap must be between zero and one")
    if redundancy_policy not in REDUNDANCY_POLICIES:
        raise ValueError(
            "redundancy_policy must be suppress_overlapping or "
            "retain_distinct_chains"
        )
    base_path = Path(base_gff_path)
    output_path = Path(output_gff_path)
    decisions_path = Path(decisions_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    for required in (base_path, *(path for _, path in inputs)):
        if not required.is_file() or required.stat().st_size == 0:
            raise ValueError(f"Missing or empty consensus input: {required}")
        if required.suffix == ".gz":
            raise ValueError("Consensus selection requires uncompressed GFF3 inputs")
    collisions = [
        path for path in (output_path, decisions_path, manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite consensus artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    grouped: dict[
        tuple[str, str, tuple[tuple[int, int, str], ...]],
        list[ConsensusInputModel],
    ] = defaultdict(list)
    input_counts: dict[str, int] = {}
    for method, path in inputs:
        models = _read_candidate_models(method, path, base_path)
        input_counts[method] = len(models)
        for model in models:
            grouped[model.chain_key].append(model)

    eligible = [
        (key, models)
        for key, models in grouped.items()
        if len({model.method for model in models}) >= min_method_support
    ]
    eligible.sort(
        key=lambda item: (
            -len({model.method for model in item[1]}),
            -item[1][0].cds_bp,
            item[0][0],
            item[1][0].start,
            _chain_digest(item[0]),
        )
    )
    chain_preserving = redundancy_policy == "retain_distinct_chains"
    conflict_annotations: dict[str, tuple[str, int]] = {}
    conflict_set_count = 0
    if chain_preserving:
        conflict_annotations, conflict_set_count = _conflict_sets(
            eligible, max_redundancy_overlap=max_redundancy_overlap
        )
    accepted: list[
        tuple[
            tuple[str, str, tuple[tuple[int, int, str], ...]],
            list[ConsensusInputModel],
        ]
    ] = []
    accepted_interval_bins: dict[tuple[str, str, int], list[int]] = defaultdict(list)
    decisions: dict[str, tuple[str, str, str]] = {}
    for key, models in eligible:
        representative = models[0]
        digest = _chain_digest(key)
        if chain_preserving:
            accepted.append((key, models))
            decisions[digest] = ("accepted", "method_support_pass", "")
            continue
        redundant_with = ""
        first_bin = (representative.start - 1) // REDUNDANCY_INTERVAL_BIN_BP
        last_bin = (representative.end - 1) // REDUNDANCY_INTERVAL_BIN_BP
        nearby_indices: set[int] = set()
        for bin_index in range(first_bin, last_bin + 1):
            nearby_indices.update(
                accepted_interval_bins.get(
                    (representative.seqid, representative.strand, bin_index), ()
                )
            )
        for accepted_index in sorted(nearby_indices):
            accepted_key, accepted_models = accepted[accepted_index]
            other = accepted_models[0]
            overlap = _overlap_bp(
                representative.cds_intervals, other.cds_intervals
            )
            denominator = min(representative.cds_bp, other.cds_bp)
            if denominator and overlap / denominator > max_redundancy_overlap:
                redundant_with = _chain_digest(accepted_key)
                break
        if redundant_with:
            decisions[digest] = (
                "rejected",
                "redundant_consensus_candidate",
                redundant_with,
            )
        else:
            accepted_index = len(accepted)
            accepted.append((key, models))
            for bin_index in range(first_bin, last_bin + 1):
                accepted_interval_bins[
                    (representative.seqid, representative.strand, bin_index)
                ].append(accepted_index)
            decisions[digest] = ("accepted", "method_support_pass", "")
    for key, models in grouped.items():
        digest = _chain_digest(key)
        decisions.setdefault(
            digest,
            ("rejected", "method_support_below_threshold", ""),
        )

    with base_path.open("rb") as base_handle:
        while chunk := base_handle.read(8 * 1024 * 1024):
            if b"ID=PPCONS_" in chunk:
                raise ValueError("Base GFF already contains PloidyPatch consensus IDs")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with base_path.open("rb") as source_handle, output_path.open("xb") as output_handle:
        shutil.copyfileobj(source_handle, output_handle, length=8 * 1024 * 1024)
        source_handle.seek(0, 2)
        if source_handle.tell():
            source_handle.seek(-1, 2)
            if source_handle.read(1) not in {b"\n", b"\r"}:
                output_handle.write(b"\n")
        output_handle.write(b"###\n")
        for key, models in accepted:
            representative = models[0]
            digest = _chain_digest(key)
            short = digest[:20]
            gene_id = f"PPCONS_gene_{short}"
            transcript_id = f"PPCONS_tx_{short}"
            methods = sorted({model.method for model in models})
            attributes = (
                f"consensus_digest={digest};"
                f"support_method_count={len(methods)};"
                f"support_methods={','.join(methods)};coding_only=true"
            )
            if chain_preserving and digest in conflict_annotations:
                conflict_digest, conflict_member_count = conflict_annotations[digest]
                attributes += (
                    f";conflict_set_digest={conflict_digest};"
                    f"conflict_member_count={conflict_member_count};"
                    "mutually_exclusive_review=true"
                )
            prefix = f"{representative.seqid}\tPloidyPatchConsensus"
            lines = [
                f"{prefix}\tgene\t{representative.start}\t{representative.end}\t.\t"
                f"{representative.strand}\t.\tID={gene_id};{attributes}\n",
                f"{prefix}\tmRNA\t{representative.start}\t{representative.end}\t.\t"
                f"{representative.strand}\t.\tID={transcript_id};Parent={gene_id};"
                f"{attributes}\n",
            ]
            for exon_number, (start, end) in enumerate(
                representative.cds_intervals, start=1
            ):
                lines.append(
                    f"{prefix}\texon\t{start}\t{end}\t.\t{representative.strand}\t.\t"
                    f"ID=PPCONS_exon_{short}_{exon_number};Parent={transcript_id}\n"
                )
            for start, end, phase in representative.cds:
                lines.append(
                    f"{prefix}\tCDS\t{start}\t{end}\t.\t{representative.strand}\t"
                    f"{phase}\tParent={transcript_id}\n"
                )
            output_handle.write("".join(lines).encode("utf-8"))

    decisions_path.parent.mkdir(parents=True, exist_ok=True)
    reason_counts: Counter[str] = Counter()
    with decisions_path.open("x", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "consensus_digest",
            "seqid",
            "start",
            "end",
            "strand",
            "cds_segments",
            "cds_bp",
            "support_method_count",
            "support_methods",
            "upstream_models",
        ]
        if chain_preserving:
            fieldnames.extend(("conflict_set_digest", "conflict_member_count"))
        fieldnames.extend((
            "status",
            "reason",
            "redundant_with_digest",
        ))
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for key, models in sorted(grouped.items(), key=lambda item: _chain_digest(item[0])):
            representative = models[0]
            digest = _chain_digest(key)
            status, reason, redundant_with = decisions[digest]
            reason_counts[reason] += 1
            methods = sorted({model.method for model in models})
            row = {
                    "consensus_digest": digest,
                    "seqid": representative.seqid,
                    "start": representative.start,
                    "end": representative.end,
                    "strand": representative.strand,
                    "cds_segments": len(representative.cds),
                    "cds_bp": representative.cds_bp,
                    "support_method_count": len(methods),
                    "support_methods": ",".join(methods),
                    "upstream_models": ",".join(
                        sorted(f"{model.method}:{model.upstream_model}" for model in models)
                    ),
                    "status": status,
                    "reason": reason,
                    "redundant_with_digest": redundant_with,
                }
            if chain_preserving:
                conflict_digest, conflict_member_count = conflict_annotations.get(
                    digest, ("", 1)
                )
                row["conflict_set_digest"] = conflict_digest
                row["conflict_member_count"] = conflict_member_count
            writer.writerow(row)

    manifest: dict[str, Any] = {
        "schema_version": (
            CHAIN_PRESERVING_SCHEMA_VERSION
            if chain_preserving
            else CONSENSUS_SCHEMA_VERSION
        ),
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "base_gff": {
                "file_name": base_path.name,
                "sha256": _file_sha256(base_path),
            },
            "candidate_methods": [
                {
                    "method": method,
                    "file_name": path.name,
                    "sha256": _file_sha256(path),
                    "models": input_counts[method],
                }
                for method, path in inputs
            ],
        },
        "parameters": {
            "min_method_support": min_method_support,
            "max_redundancy_overlap": max_redundancy_overlap,
            "redundancy_interval_bin_bp": REDUNDANCY_INTERVAL_BIN_BP,
            "redundancy_index": "seqid_strand_fixed_interval_bins",
            "consensus_unit": "exact_seqid_strand_phased_cds_chain",
            "support_unit": "independent_method_family",
            "output_scope": "coding_model_only",
        },
        "counts": {
            "input_models": sum(input_counts.values()),
            "unique_cds_chains": len(grouped),
            "support_passing_chains": len(eligible),
            "accepted_models": len(accepted),
            "decision_counts": dict(sorted(reason_counts.items())),
        },
        "outputs": {
            "candidate_gff": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
            },
            "decisions": {
                "file_name": decisions_path.name,
                "sha256": _file_sha256(decisions_path),
                "rows": len(grouped),
            },
        },
    }
    if chain_preserving:
        manifest["parameters"].update(
            {
                "redundancy_policy": redundancy_policy,
                "conflict_definition": (
                    "connected_component_of_same_strand_chains_with_"
                    "cds_overlap_over_shorter_above_threshold"
                ),
                "conflict_action": "retain_all_for_ranking_and_review",
            }
        )
        manifest["counts"].update(
            {
                "conflict_sets": conflict_set_count,
                "conflicted_chains": len(conflict_annotations),
            }
        )
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
