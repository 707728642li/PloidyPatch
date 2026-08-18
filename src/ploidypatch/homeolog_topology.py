from __future__ import annotations

import csv
import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .baseline import _file_sha256
from .gff import TRANSCRIPT_TYPES, parse_attributes
from .io import open_text


HOMEOLOG_TOPOLOGY_SCHEMA_VERSION = "ploidypatch.homeolog_topology.v1"
HOMEOLOG_TOPOLOGY_COLUMNS = (
    "candidate_digest",
    "topology_available",
    "topology_reason",
    "candidate_gene_id",
    "partner_gene_id",
    "pair_id",
    "candidate_transcript_id",
    "partner_transcript_id",
    "candidate_cds_bp",
    "partner_cds_bp",
    "cds_bp_ratio",
    "candidate_cds_segments",
    "partner_cds_segments",
    "cds_segment_count_ratio",
    "phase_lcs_similarity",
    "junction_fraction_similarity",
    "candidate_coding_span_bp",
    "partner_coding_span_bp",
    "coding_span_ratio",
    "topology_coherence_score",
    "wgd_support_block_count",
    "wgd_longest_block_pairs",
)


@dataclass(frozen=True)
class TranscriptTopology:
    transcript_id: str
    seqid: str
    strand: str
    cds: tuple[tuple[int, int, str], ...]

    @property
    def cds_bp(self) -> int:
        return sum(end - start + 1 for start, end, _ in self.cds)

    @property
    def coding_span_bp(self) -> int:
        return max(end for _, end, _ in self.cds) - min(
            start for start, _, _ in self.cds
        ) + 1

    @property
    def phases(self) -> tuple[str, ...]:
        return tuple(phase for _, _, phase in self.cds)

    @property
    def junction_fractions(self) -> tuple[float, ...]:
        lengths = [end - start + 1 for start, end, _ in self.cds]
        total = sum(lengths)
        cumulative = 0
        fractions: list[float] = []
        for length in lengths[:-1]:
            cumulative += length
            fractions.append(cumulative / total)
        return tuple(fractions)


def _read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(f"{path} lacks required columns: {sorted(required)}")
        return list(reader)


def _read_gff_topologies(
    path: Path,
) -> tuple[dict[str, TranscriptTopology], dict[str, tuple[str, ...]]]:
    transcript_meta: dict[str, tuple[str, str]] = {}
    gene_to_transcripts: dict[str, set[str]] = defaultdict(set)
    cds_by_transcript: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GFF row in {path} line {line_number}")
            seqid, _, feature_type, start_text, end_text, _, strand, phase, raw = fields
            attributes, malformed = parse_attributes(raw)
            if malformed:
                raise ValueError(
                    f"Malformed GFF attributes in {path} line {line_number}"
                )
            if feature_type in TRANSCRIPT_TYPES:
                transcript_id = attributes.get("ID", "")
                if not transcript_id or transcript_id in transcript_meta:
                    raise ValueError(
                        f"Empty or duplicate transcript ID in {path} line {line_number}"
                    )
                transcript_meta[transcript_id] = (seqid, strand)
                for parent in attributes.get("Parent", "").split(","):
                    if parent:
                        gene_to_transcripts[parent].add(transcript_id)
            elif feature_type == "CDS":
                try:
                    start = int(start_text)
                    end = int(end_text)
                except ValueError as exc:
                    raise ValueError(
                        f"Non-numeric CDS coordinate in {path} line {line_number}"
                    ) from exc
                if start < 1 or end < start or phase not in {"0", "1", "2", "."}:
                    raise ValueError(f"Invalid CDS in {path} line {line_number}")
                parents = [
                    value
                    for value in attributes.get("Parent", "").split(",")
                    if value
                ]
                if not parents:
                    raise ValueError(f"CDS lacks Parent in {path} line {line_number}")
                for parent in parents:
                    cds_by_transcript[parent].append((start, end, phase))

    topologies: dict[str, TranscriptTopology] = {}
    for transcript_id, cds in cds_by_transcript.items():
        if transcript_id not in transcript_meta:
            continue
        seqid, strand = transcript_meta[transcript_id]
        ordered = tuple(sorted(cds, reverse=strand == "-"))
        topologies[transcript_id] = TranscriptTopology(
            transcript_id=transcript_id,
            seqid=seqid,
            strand=strand,
            cds=ordered,
        )
    indexed_genes = {
        gene_id: tuple(
            transcript_id
            for transcript_id in sorted(transcript_ids)
            if transcript_id in topologies
        )
        for gene_id, transcript_ids in gene_to_transcripts.items()
    }
    return topologies, indexed_genes


def _locus_transcripts(
    locus_id: str,
    topologies: dict[str, TranscriptTopology],
    gene_to_transcripts: dict[str, tuple[str, ...]],
) -> tuple[TranscriptTopology, ...]:
    if locus_id in topologies:
        return (topologies[locus_id],)
    return tuple(
        topologies[transcript_id]
        for transcript_id in gene_to_transcripts.get(locus_id, ())
    )


def _ratio(left: int, right: int) -> float:
    return min(left, right) / max(left, right)


def _lcs_similarity(left: tuple[str, ...], right: tuple[str, ...]) -> float:
    if not left and not right:
        return 1.0
    previous = [0] * (len(right) + 1)
    for left_value in left:
        current = [0]
        for index, right_value in enumerate(right, start=1):
            if left_value == right_value:
                current.append(previous[index - 1] + 1)
            else:
                current.append(max(current[-1], previous[index]))
        previous = current
    return 2.0 * previous[-1] / (len(left) + len(right))


def _junction_similarity(left: tuple[float, ...], right: tuple[float, ...]) -> float:
    if not left and not right:
        return 1.0
    if not left or not right:
        return 0.0

    def directed(first: tuple[float, ...], second: tuple[float, ...]) -> float:
        return sum(min(abs(value - other) for other in second) for value in first) / len(
            first
        )

    distance = (directed(left, right) + directed(right, left)) / 2.0
    return max(0.0, 1.0 - distance)


def _pair_metrics(
    candidate: TranscriptTopology, partner: TranscriptTopology
) -> dict[str, float | int | str]:
    cds_ratio = _ratio(candidate.cds_bp, partner.cds_bp)
    segment_ratio = _ratio(len(candidate.cds), len(partner.cds))
    phase_similarity = _lcs_similarity(candidate.phases, partner.phases)
    junction_similarity = _junction_similarity(
        candidate.junction_fractions, partner.junction_fractions
    )
    span_ratio = _ratio(candidate.coding_span_bp, partner.coding_span_bp)
    score = (
        cds_ratio
        + segment_ratio
        + phase_similarity
        + junction_similarity
        + span_ratio
    ) / 5.0
    return {
        "candidate_transcript_id": candidate.transcript_id,
        "partner_transcript_id": partner.transcript_id,
        "candidate_cds_bp": candidate.cds_bp,
        "partner_cds_bp": partner.cds_bp,
        "cds_bp_ratio": cds_ratio,
        "candidate_cds_segments": len(candidate.cds),
        "partner_cds_segments": len(partner.cds),
        "cds_segment_count_ratio": segment_ratio,
        "phase_lcs_similarity": phase_similarity,
        "junction_fraction_similarity": junction_similarity,
        "candidate_coding_span_bp": candidate.coding_span_bp,
        "partner_coding_span_bp": partner.coding_span_bp,
        "coding_span_ratio": span_ratio,
        "topology_coherence_score": score,
    }


def _best_pair(
    candidates: tuple[TranscriptTopology, ...],
    partners: tuple[TranscriptTopology, ...],
) -> dict[str, float | int | str]:
    best: dict[str, float | int | str] | None = None
    for candidate in sorted(candidates, key=lambda item: item.transcript_id):
        for partner in sorted(partners, key=lambda item: item.transcript_id):
            metrics = _pair_metrics(candidate, partner)
            if best is None or (
                metrics["topology_coherence_score"],
                metrics["phase_lcs_similarity"],
                metrics["cds_bp_ratio"],
                metrics["cds_segment_count_ratio"],
            ) > (
                best["topology_coherence_score"],
                best["phase_lcs_similarity"],
                best["cds_bp_ratio"],
                best["cds_segment_count_ratio"],
            ):
                best = metrics
    if best is None:
        raise ValueError("Cannot choose a topology pair from an empty locus")
    return best


def _format(value: float | int | str) -> str:
    if isinstance(value, float):
        return format(value, ".12g")
    return str(value)


def build_homeolog_topology_features(
    *,
    copy_feature_tsv_path: str | Path,
    wgd_selection_tsv_path: str | Path,
    candidate_gff_path: str | Path,
    base_gff_path: str | Path,
    output_tsv_path: str | Path,
) -> dict[str, Any]:
    """Compare candidate CDS topology with its truth-blind existing WGD partner."""

    feature_path = Path(copy_feature_tsv_path)
    selection_path = Path(wgd_selection_tsv_path)
    candidate_path = Path(candidate_gff_path)
    base_path = Path(base_gff_path)
    output_path = Path(output_tsv_path)
    manifest_path = Path(str(output_path) + ".manifest.json")
    collisions = [path for path in (output_path, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite homeolog-topology artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    feature_rows = _read_tsv(feature_path, {"candidate_digest"})
    selection_rows = _read_tsv(
        selection_path,
        {
            "gene_id",
            "consensus_digest",
            "pair_id",
            "partner_gene_id",
            "support_block_count",
            "longest_block_pairs",
            "status",
            "reason",
        },
    )
    selections: dict[str, dict[str, str]] = {}
    for row in selection_rows:
        digest = row["consensus_digest"]
        if not digest or digest in selections:
            raise ValueError("WGD selection has an empty or duplicate digest")
        selections[digest] = row

    candidate_topologies, candidate_genes = _read_gff_topologies(candidate_path)
    base_topologies, base_genes = _read_gff_topologies(base_path)
    output_rows: list[dict[str, str]] = []
    accepted = 0
    for feature in feature_rows:
        digest = feature["candidate_digest"]
        selection = selections.get(digest)
        if selection is None:
            raise ValueError(f"WGD selection lacks candidate digest {digest}")
        output = {column: "" for column in HOMEOLOG_TOPOLOGY_COLUMNS}
        output.update(
            {
                "candidate_digest": digest,
                "topology_available": "0",
                "topology_reason": selection["reason"],
                "candidate_gene_id": selection["gene_id"],
                "partner_gene_id": selection["partner_gene_id"],
                "pair_id": selection["pair_id"],
                "wgd_support_block_count": selection["support_block_count"],
                "wgd_longest_block_pairs": selection["longest_block_pairs"],
            }
        )
        if selection["status"] == "accepted":
            candidate_models = _locus_transcripts(
                selection["gene_id"], candidate_topologies, candidate_genes
            )
            partner_models = _locus_transcripts(
                selection["partner_gene_id"], base_topologies, base_genes
            )
            if not candidate_models:
                raise ValueError(
                    f"Accepted WGD candidate lacks coding topology: {digest}"
                )
            if not partner_models:
                raise ValueError(
                    f"Accepted WGD partner lacks coding topology: "
                    f"{selection['partner_gene_id']}"
                )
            output.update(
                {key: _format(value) for key, value in _best_pair(
                    candidate_models, partner_models
                ).items()}
            )
            output["topology_available"] = "1"
            output["topology_reason"] = "accepted_existing_wgd_partner"
            accepted += 1
        output_rows.append(output)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=HOMEOLOG_TOPOLOGY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(output_rows)

    manifest: dict[str, Any] = {
        "schema_version": HOMEOLOG_TOPOLOGY_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "truth_access": False,
        "selection_policy": (
            "maximize_mean_cds_length_segment_phase_junction_and_span_similarity_"
            "across_candidate_partner_isoforms"
        ),
        "inputs": {
            "copy_features": _file_sha256(feature_path),
            "wgd_selection": _file_sha256(selection_path),
            "candidate_gff": _file_sha256(candidate_path),
            "base_gff": _file_sha256(base_path),
        },
        "counts": {
            "candidates": len(output_rows),
            "topology_available": accepted,
            "topology_unavailable": len(output_rows) - accepted,
        },
        "outputs": {
            "features": {
                "file_name": output_path.name,
                "sha256": _file_sha256(output_path),
                "rows": len(output_rows),
            }
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
