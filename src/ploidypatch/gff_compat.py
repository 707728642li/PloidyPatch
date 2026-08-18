from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from . import __version__
from .io import open_text


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def _attributes(text: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in text.split(";"):
        key, separator, value = item.partition("=")
        if separator and key and key not in parsed:
            parsed[key] = value
    return parsed


def synthesize_root_transcript_genes(
    input_gff: str | Path,
    output_gff: str | Path,
    mapping_tsv: str | Path,
) -> dict[str, Any]:
    """Add stable gene parents when a GFF contains only root transcripts."""

    input_path = Path(input_gff)
    output_path = Path(output_gff)
    mapping_path = Path(mapping_tsv)
    manifest_path = Path(str(output_path) + ".manifest.json")
    collisions = [
        path for path in (output_path, mapping_path, manifest_path) if path.exists()
    ]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite root-transcript compatibility artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    with open_text(input_path) as handle:
        lines = handle.readlines()
    feature_ids: set[str] = set()
    gene_records = 0
    for raw_line in lines:
        if raw_line.startswith("#") or not raw_line.strip():
            continue
        fields = raw_line.rstrip("\r\n").split("\t")
        if len(fields) != 9:
            continue
        attrs = _attributes(fields[8])
        if attrs.get("ID"):
            feature_ids.add(attrs["ID"])
        if fields[2] == "gene":
            gene_records += 1
    if gene_records:
        raise ValueError(
            "Root-transcript compatibility mode requires an annotation with zero gene records"
        )

    transcript_to_gene: dict[str, str] = {}
    gene_groups: dict[str, dict[str, Any]] = {}
    generated_gene_ids: set[str] = set()
    for line_number, raw_line in enumerate(lines, start=1):
        if raw_line.startswith("#") or not raw_line.strip():
            continue
        fields = raw_line.rstrip("\r\n").split("\t")
        if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
            continue
        attrs = _attributes(fields[8])
        if attrs.get("Parent"):
            continue
        transcript_id = attrs.get("ID", "")
        if not transcript_id:
            raise ValueError(f"Root transcript lacks ID at line {line_number}")
        group_key = attrs.get("geneID") or transcript_id
        proposed_gene_id = group_key
        if proposed_gene_id in feature_ids or proposed_gene_id in generated_gene_ids:
            token = hashlib.sha256(group_key.encode("utf-8")).hexdigest()[:20]
            proposed_gene_id = f"PPROOTGENE-{token}"
        group = gene_groups.get(group_key)
        if group is None:
            if proposed_gene_id in feature_ids or proposed_gene_id in generated_gene_ids:
                raise ValueError(f"Synthesized gene ID collision: {proposed_gene_id}")
            generated_gene_ids.add(proposed_gene_id)
            group = {
                "gene_id": proposed_gene_id,
                "seqid": fields[0],
                "source": fields[1],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "name": attrs.get("Name", group_key),
            }
            gene_groups[group_key] = group
        elif fields[0] != group["seqid"] or fields[6] != group["strand"]:
            raise ValueError(
                f"Root transcript group spans incompatible loci: {group_key}"
            )
        else:
            group["start"] = min(group["start"], int(fields[3]))
            group["end"] = max(group["end"], int(fields[4]))
        transcript_to_gene[transcript_id] = str(group["gene_id"])

    output_lines: list[str] = []
    mappings: list[tuple[str, str]] = []
    emitted_genes: set[str] = set()
    for line_number, raw_line in enumerate(lines, start=1):
        if raw_line.startswith("#") or not raw_line.strip():
            output_lines.append(raw_line)
            continue
        ending = "\r\n" if raw_line.endswith("\r\n") else "\n" if raw_line.endswith("\n") else ""
        fields = raw_line.rstrip("\r\n").split("\t")
        if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
            output_lines.append(raw_line)
            continue
        attrs = _attributes(fields[8])
        if attrs.get("Parent"):
            output_lines.append(raw_line)
            continue
        transcript_id = attrs.get("ID", "")
        if not transcript_id:
            raise ValueError(f"Root transcript lacks ID at line {line_number}")
        gene_id = transcript_to_gene[transcript_id]
        if gene_id not in emitted_genes:
            group_key = attrs.get("geneID") or transcript_id
            group = gene_groups[group_key]
            gene_fields = fields.copy()
            gene_fields[2] = "gene"
            gene_fields[3] = str(group["start"])
            gene_fields[4] = str(group["end"])
            gene_fields[7] = "."
            gene_fields[8] = (
                f"ID={gene_id};Name={group['name']};"
                f"ploidypatch_source_gene_key={group_key}"
            )
            output_lines.append("\t".join(gene_fields) + ending)
            emitted_genes.add(gene_id)
        fields[8] = fields[8].rstrip(";") + f";Parent={gene_id};"
        output_lines.append("\t".join(fields) + ending)
        mappings.append((transcript_id, gene_id))
    if not mappings:
        raise ValueError("Annotation contains no parentless root transcripts")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        handle.write("".join(output_lines))
    with mapping_path.open("x", encoding="utf-8", newline="") as handle:
        handle.write("transcript_id\tsynthesized_gene_id\n")
        for transcript_id, gene_id in mappings:
            handle.write(f"{transcript_id}\t{gene_id}\n")
    report: dict[str, Any] = {
        "schema_version": "ploidypatch.root_transcript_gene_compat.v1",
        "generator": {"name": "PloidyPatch", "version": __version__},
        "policy": "only_when_input_contains_zero_gene_records",
        "counts": {
            "input_gene_records": 0,
            "root_transcripts": len(mappings),
            "synthesized_gene_records": len(emitted_genes),
        },
        "inputs": {"gff": {"file_name": input_path.name, "sha256": _sha256(input_path)}},
        "outputs": {
            "gff": {"file_name": output_path.name, "sha256": _sha256(output_path)},
            "mapping": {"file_name": mapping_path.name, "sha256": _sha256(mapping_path)},
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def synthesize_missing_transcript_exons(
    input_gff: str | Path,
    output_gff: str | Path,
    *,
    repair_parent_bounds: bool = False,
) -> dict[str, Any]:
    """Add deterministic exon rows only to transcripts that have none.

    Some provider GFF3 files represent transcript structure with CDS and UTR
    children but omit the redundant ``exon`` rows.  Homology tools such as
    LiftOn require exon children to extract transcript sequences.  This narrow
    compatibility adapter preserves every child-feature byte and appends exon
    rows derived solely from the union of that transcript's CDS/UTR intervals.
    The explicit ``repair_parent_bounds`` compatibility mode may expand gene
    and transcript start/end columns to contain same-seqid, same-strand child
    unions.  It never changes CDS/UTR/exon rows, phases, identifiers or
    parentage, and it never shrinks a parent interval.
    """

    input_path = Path(input_gff)
    output_path = Path(output_gff)
    manifest_path = Path(str(output_path) + ".manifest.json")
    collisions = [path for path in (output_path, manifest_path) if path.exists()]
    if collisions:
        raise FileExistsError(
            "Refusing to overwrite transcript-exon compatibility artifact(s): "
            + ", ".join(str(path) for path in collisions)
        )

    with open_text(input_path) as handle:
        lines = handle.readlines()

    transcript_types = {"mRNA", "transcript"}
    segment_types = {
        "CDS",
        "five_prime_UTR",
        "three_prime_UTR",
        "five_prime_utr",
        "three_prime_utr",
        "UTR",
        "utr",
        "start_codon",
        "stop_codon",
    }
    transcripts: dict[str, dict[str, Any]] = {}
    genes: dict[str, dict[str, Any]] = {}
    feature_records: dict[str, list[dict[str, Any]]] = {}
    feature_ids: set[str] = set()
    existing_exon_parents: set[str] = set()
    segments: dict[str, list[tuple[str, str, int, int, str]]] = {}
    structural_children: dict[str, list[tuple[str, int, int, str, str]]] = {}

    for line_number, raw_line in enumerate(lines, start=1):
        if raw_line.startswith("#") or not raw_line.strip():
            continue
        fields = raw_line.rstrip("\r\n").split("\t")
        if len(fields) != 9:
            raise ValueError(f"Malformed GFF3 row at line {line_number}")
        attrs = _attributes(fields[8])
        feature_id = attrs.get("ID", "")
        parents = tuple(
            item for item in attrs.get("Parent", "").split(",") if item
        )
        if feature_id:
            feature_ids.add(feature_id)
            feature_records.setdefault(feature_id, []).append(
                {
                    "feature_type": fields[2],
                    "seqid": fields[0],
                    "source": fields[1],
                    "start_text": fields[3],
                    "end_text": fields[4],
                    "strand": fields[6],
                    "parents": parents,
                    "line_index": line_number - 1,
                }
            )
        feature_type = fields[2]
        if repair_parent_bounds and feature_type == "gene":
            if not feature_id:
                raise ValueError(f"Gene lacks ID at line {line_number}")
            if feature_id in genes:
                raise ValueError(f"Duplicate gene feature ID: {feature_id}")
            start = int(fields[3])
            end = int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f"Invalid gene interval at line {line_number}")
            genes[feature_id] = {
                "seqid": fields[0],
                "start": start,
                "end": end,
                "strand": fields[6],
                "line_index": line_number - 1,
            }
        if feature_type in transcript_types:
            if not feature_id:
                raise ValueError(f"Transcript lacks ID at line {line_number}")
            if feature_id in transcripts:
                raise ValueError(f"Duplicate transcript feature ID: {feature_id}")
            transcript_start = int(fields[3])
            transcript_end = int(fields[4])
            if transcript_start < 1 or transcript_end < transcript_start:
                raise ValueError(f"Invalid transcript interval at line {line_number}")
            transcripts[feature_id] = {
                "seqid": fields[0],
                "source": fields[1],
                "start": transcript_start,
                "end": transcript_end,
                "strand": fields[6],
                "parents": parents,
                "line_index": line_number - 1,
            }
        if feature_type == "exon":
            existing_exon_parents.update(parents)
        elif feature_type in segment_types:
            start = int(fields[3])
            end = int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f"Invalid segment interval at line {line_number}")
            for parent in parents:
                segments.setdefault(parent, []).append(
                    (fields[0], fields[1], start, end, fields[6])
                )
        if repair_parent_bounds and (feature_type == "exon" or feature_type in segment_types):
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as error:
                raise ValueError(
                    f"Invalid structural-child coordinates at line {line_number}"
                ) from error
            if start < 1 or end < start:
                raise ValueError(f"Invalid structural-child interval at line {line_number}")
            if not parents:
                raise ValueError(f"Structural child lacks Parent at line {line_number}")
            for parent in parents:
                structural_children.setdefault(parent, []).append(
                    (fields[0], start, end, fields[6], feature_type)
                )

    parent_bound_repairs: list[dict[str, Any]] = []
    repaired_lines: dict[int, tuple[int, int]] = {}
    promoted_transcript_like_types: dict[str, int] = {}
    promoted_gene_like_root_types: dict[str, int] = {}
    preserved_exon_only_container_types: dict[str, int] = {}
    if repair_parent_bounds:
        unresolved_gene_parents = sorted(
            {
                parent
                for transcript in transcripts.values()
                for parent in transcript["parents"]
                if parent not in genes
            }
        )
        for parent in unresolved_gene_parents:
            records = feature_records.get(parent, [])
            if len(records) != 1:
                continue
            feature = records[0]
            feature_type = str(feature["feature_type"])
            if feature_type in transcript_types or tuple(feature["parents"]):
                continue
            try:
                feature_start = int(str(feature["start_text"]))
                feature_end = int(str(feature["end_text"]))
            except ValueError as error:
                raise ValueError(
                    f"Invalid gene-like root feature interval: {parent}"
                ) from error
            if feature_start < 1 or feature_end < feature_start:
                raise ValueError(f"Invalid gene-like root feature interval: {parent}")
            child_transcripts = [
                item
                for item in transcripts.values()
                if tuple(item["parents"]) == (parent,)
            ]
            if not child_transcripts:
                continue
            for transcript in child_transcripts:
                if transcript["seqid"] != feature["seqid"]:
                    raise ValueError(
                        f"Gene-like root transcript changes seqid: {parent}"
                    )
                if transcript["strand"] != feature["strand"]:
                    raise ValueError(
                        f"Gene-like root transcript changes strand: {parent}"
                    )
                if (
                    transcript["start"] < feature_start
                    or transcript["end"] > feature_end
                ):
                    raise ValueError(
                        f"Gene-like root transcript lies outside parent bounds: {parent}"
                    )
            genes[parent] = {
                "seqid": feature["seqid"],
                "start": feature_start,
                "end": feature_end,
                "strand": feature["strand"],
                "line_index": feature["line_index"],
            }
            promoted_gene_like_root_types[feature_type] = (
                promoted_gene_like_root_types.get(feature_type, 0) + 1
            )

        unknown_child_parents = sorted(set(structural_children) - set(transcripts))
        remaining_unknown: list[str] = []
        for parent in unknown_child_parents:
            records = feature_records.get(parent, [])
            if len(records) != 1:
                remaining_unknown.append(parent)
                continue
            feature = records[0]
            feature_type = str(feature["feature_type"])
            feature_parents = tuple(feature["parents"])
            children = structural_children[parent]
            if all(child[4] == "exon" for child in children):
                try:
                    feature_start = int(str(feature["start_text"]))
                    feature_end = int(str(feature["end_text"]))
                except ValueError as error:
                    raise ValueError(
                        f"Invalid exon-only feature interval: {parent}"
                    ) from error
                if feature_start < 1 or feature_end < feature_start:
                    raise ValueError(f"Invalid exon-only feature interval: {parent}")
                for seqid, start, end, strand, _feature_type in children:
                    if seqid != feature["seqid"]:
                        raise ValueError(
                            f"Exon-only feature child changes seqid: {parent}"
                        )
                    if strand != feature["strand"]:
                        raise ValueError(
                            f"Exon-only feature child changes strand: {parent}"
                        )
                    if start < feature_start or end > feature_end:
                        raise ValueError(
                            f"Exon-only feature child lies outside parent bounds: {parent}"
                        )
                preserved_exon_only_container_types[feature_type] = (
                    preserved_exon_only_container_types.get(feature_type, 0) + 1
                )
                continue
            if (
                feature_type == "gene"
                or len(feature_parents) != 1
                or feature_parents[0] not in genes
            ):
                remaining_unknown.append(parent)
                continue
            try:
                feature_start = int(str(feature["start_text"]))
                feature_end = int(str(feature["end_text"]))
            except ValueError as error:
                raise ValueError(
                    f"Invalid transcript-like feature interval: {parent}"
                ) from error
            if feature_start < 1 or feature_end < feature_start:
                raise ValueError(f"Invalid transcript-like feature interval: {parent}")
            gene = genes[feature_parents[0]]
            if feature["seqid"] != gene["seqid"]:
                raise ValueError(
                    f"Transcript-like feature changes gene seqid: {parent}"
                )
            if feature["strand"] != gene["strand"]:
                raise ValueError(
                    f"Transcript-like feature changes gene strand: {parent}"
                )
            transcripts[parent] = {
                "seqid": feature["seqid"],
                "source": feature["source"],
                "start": feature_start,
                "end": feature_end,
                "strand": feature["strand"],
                "parents": feature_parents,
                "line_index": feature["line_index"],
            }
            promoted_transcript_like_types[feature_type] = (
                promoted_transcript_like_types.get(feature_type, 0) + 1
            )
        unknown_child_parents = remaining_unknown
        if unknown_child_parents:
            raise ValueError(
                "Structural child Parent is not a transcript during bound repair: "
                + ", ".join(unknown_child_parents[:10])
            )
    if not transcripts:
        raise ValueError("Annotation contains no transcript records")

    if repair_parent_bounds:
        gene_transcripts: dict[str, list[str]] = {}
        for transcript_id, transcript in transcripts.items():
            parents = transcript["parents"]
            if len(parents) != 1:
                raise ValueError(
                    "Parent-bound repair requires exactly one gene Parent for transcript: "
                    + transcript_id
                )
            gene_id = str(parents[0])
            if gene_id not in genes:
                raise ValueError(
                    f"Transcript gene Parent is unresolved during bound repair: {transcript_id}"
                )
            gene_transcripts.setdefault(gene_id, []).append(transcript_id)
            children = structural_children.get(transcript_id, [])
            for seqid, _start, _end, strand, _feature_type in children:
                if seqid != transcript["seqid"]:
                    raise ValueError(
                        f"Transcript child changes seqid during parent-bound repair: {transcript_id}"
                    )
                if strand != transcript["strand"]:
                    raise ValueError(
                        f"Transcript child changes strand during parent-bound repair: {transcript_id}"
                    )
            if children:
                repaired_start = min(
                    int(transcript["start"]), *(start for _, start, _, _, _ in children)
                )
                repaired_end = max(
                    int(transcript["end"]), *(end for _, _, end, _, _ in children)
                )
                if (repaired_start, repaired_end) != (
                    transcript["start"],
                    transcript["end"],
                ):
                    parent_bound_repairs.append(
                        {
                            "feature_type": "transcript",
                            "feature_id": transcript_id,
                            "seqid": transcript["seqid"],
                            "strand": transcript["strand"],
                            "input_start": transcript["start"],
                            "input_end": transcript["end"],
                            "output_start": repaired_start,
                            "output_end": repaired_end,
                        }
                    )
                    transcript["start"] = repaired_start
                    transcript["end"] = repaired_end
                    repaired_lines[int(transcript["line_index"])] = (
                        repaired_start,
                        repaired_end,
                    )
        for gene_id, transcript_ids in gene_transcripts.items():
            gene = genes[gene_id]
            for transcript_id in transcript_ids:
                transcript = transcripts[transcript_id]
                if transcript["seqid"] != gene["seqid"]:
                    raise ValueError(
                        f"Gene transcript changes seqid during parent-bound repair: {gene_id}"
                    )
                if transcript["strand"] != gene["strand"]:
                    raise ValueError(
                        f"Gene transcript changes strand during parent-bound repair: {gene_id}"
                    )
            repaired_start = min(
                int(gene["start"]),
                *(int(transcripts[item]["start"]) for item in transcript_ids),
            )
            repaired_end = max(
                int(gene["end"]),
                *(int(transcripts[item]["end"]) for item in transcript_ids),
            )
            if (repaired_start, repaired_end) != (gene["start"], gene["end"]):
                parent_bound_repairs.append(
                    {
                        "feature_type": "gene",
                        "feature_id": gene_id,
                        "seqid": gene["seqid"],
                        "strand": gene["strand"],
                        "input_start": gene["start"],
                        "input_end": gene["end"],
                        "output_start": repaired_start,
                        "output_end": repaired_end,
                    }
                )
                repaired_lines[int(gene["line_index"])] = (repaired_start, repaired_end)

    synthesized: list[str] = []
    synthesized_transcripts = 0
    unresolved: list[str] = []
    for transcript_id in sorted(transcripts):
        if transcript_id in existing_exon_parents:
            continue
        child_segments = segments.get(transcript_id, [])
        if not child_segments:
            unresolved.append(transcript_id)
            continue
        transcript = transcripts[transcript_id]
        intervals: list[tuple[int, int]] = []
        for seqid, _source, start, end, strand in child_segments:
            if seqid != transcript["seqid"] or strand != transcript["strand"]:
                raise ValueError(
                    f"Transcript child changes seqid/strand: {transcript_id}"
                )
            if start < transcript["start"] or end > transcript["end"]:
                raise ValueError(
                    f"Transcript child lies outside parent bounds: {transcript_id}"
                )
            intervals.append((start, end))
        intervals.sort()
        merged: list[list[int]] = []
        for start, end in intervals:
            if not merged or start > merged[-1][1] + 1:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        for ordinal, (start, end) in enumerate(merged, start=1):
            token = hashlib.sha256(
                f"{transcript_id}\0{start}\0{end}".encode("utf-8")
            ).hexdigest()[:20]
            # LiftOn treats any trailing ``-N`` exon ID as a renumberable
            # family and, after chaining, may collapse a shared prefix to
            # ``prefix-1``, ``prefix-2`` independently in every transcript.
            # A delimiter-free hash keeps the source ID globally unique and,
            # if LiftOn must rebuild a transcript, forces its transcript-ID-
            # based fallback rather than a cross-transcript shared family.
            exon_id = f"PPX{token}"
            if exon_id in feature_ids:
                raise ValueError(f"Synthesized exon ID collision: {exon_id}")
            feature_ids.add(exon_id)
            synthesized.append(
                "\t".join(
                    (
                        str(transcript["seqid"]),
                        "PloidyPatch",
                        "exon",
                        str(start),
                        str(end),
                        ".",
                        str(transcript["strand"]),
                        ".",
                        (
                            f"ID={exon_id};Parent={transcript_id};"
                            f"ploidypatch_compat=synthesized_missing_exon;"
                            f"ploidypatch_exon_ordinal={ordinal}"
                        ),
                    )
                )
                + "\n"
            )
        synthesized_transcripts += 1

    output_lines = list(lines)
    for line_index, (start, end) in repaired_lines.items():
        raw_line = output_lines[line_index]
        ending = "\r\n" if raw_line.endswith("\r\n") else "\n" if raw_line.endswith("\n") else ""
        fields = raw_line.rstrip("\r\n").split("\t")
        fields[3] = str(start)
        fields[4] = str(end)
        output_lines[line_index] = "\t".join(fields) + ending

    protected_child_types = {"CDS", "exon", *segment_types}
    protected_child_records = 0
    for line_number, (input_line, output_line) in enumerate(
        zip(lines, output_lines, strict=True), start=1
    ):
        fields = input_line.rstrip("\r\n").split("\t")
        if len(fields) == 9 and fields[2] in protected_child_types:
            protected_child_records += 1
            if output_line != input_line:
                raise AssertionError(
                    "Parent-bound compatibility changed a protected CDS/UTR/exon row "
                    f"at line {line_number}"
                )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        handle.write("".join(output_lines))
        if output_lines and not output_lines[-1].endswith(("\n", "\r")):
            handle.write("\n")
        if synthesized:
            handle.write("# PloidyPatch deterministic missing-exon compatibility rows\n")
            handle.write("".join(synthesized))

    def cds_signature(path: Path) -> tuple[int, str]:
        digest = hashlib.sha256()
        count = 0
        with open_text(path) as handle:
            for raw_line in handle:
                fields = raw_line.rstrip("\r\n").split("\t")
                if len(fields) == 9 and fields[2] == "CDS":
                    count += 1
                    digest.update(raw_line.encode("utf-8"))
        return count, digest.hexdigest()

    input_cds_count, input_cds_sha = cds_signature(input_path)
    output_cds_count, output_cds_sha = cds_signature(output_path)
    if (input_cds_count, input_cds_sha) != (output_cds_count, output_cds_sha):
        raise AssertionError("Missing-exon compatibility adapter changed CDS rows")

    transcript_repairs = sum(
        item["feature_type"] == "transcript" for item in parent_bound_repairs
    )
    gene_repairs = sum(item["feature_type"] == "gene" for item in parent_bound_repairs)
    report: dict[str, Any] = {
        "schema_version": (
            "ploidypatch.missing_transcript_exon_compat.v2"
            if repair_parent_bounds
            else "ploidypatch.missing_transcript_exon_compat.v1"
        ),
        "generator": {"name": "PloidyPatch", "version": __version__},
        "policy": (
            "expand_gene_and_transcript_bounds_only_to_same_seqid_same_strand_child_union_then_"
            "append_union_of_CDS_UTR_intervals_only_when_transcript_has_zero_exons"
            if repair_parent_bounds
            else "append_union_of_CDS_UTR_intervals_only_when_transcript_has_zero_exons"
        ),
        "coordinate_or_cds_changes": bool(parent_bound_repairs),
        "child_coordinate_or_cds_changes": False,
        "repair_parent_bounds": repair_parent_bounds,
        "promoted_transcript_like_feature_types": dict(
            sorted(promoted_transcript_like_types.items())
        ),
        "promoted_gene_like_root_feature_types": dict(
            sorted(promoted_gene_like_root_types.items())
        ),
        "preserved_exon_only_feature_container_types": dict(
            sorted(preserved_exon_only_container_types.items())
        ),
        "counts": {
            "transcripts": len(transcripts),
            "promoted_transcript_like_features": sum(
                promoted_transcript_like_types.values()
            ),
            "promoted_gene_like_root_features": sum(
                promoted_gene_like_root_types.values()
            ),
            "preserved_exon_only_feature_containers": sum(
                preserved_exon_only_container_types.values()
            ),
            "transcripts_with_existing_exons": len(
                set(transcripts) & existing_exon_parents
            ),
            "transcripts_with_synthesized_exons": synthesized_transcripts,
            "synthesized_exon_records": len(synthesized),
            "unresolved_transcripts": len(unresolved),
            "input_cds_records": input_cds_count,
            "output_cds_records": output_cds_count,
            "transcript_parent_bounds_repaired": transcript_repairs,
            "gene_parent_bounds_repaired": gene_repairs,
            "parent_bounds_repaired": len(parent_bound_repairs),
            "protected_input_child_records_verified_unchanged": protected_child_records,
        },
        "cds_rows_sha256": {
            "input": input_cds_sha,
            "output": output_cds_sha,
        },
        "unresolved_transcript_ids": unresolved,
        "parent_bound_repairs": parent_bound_repairs,
        "inputs": {"gff": {"file_name": input_path.name, "sha256": _sha256(input_path)}},
        "outputs": {
            "gff": {"file_name": output_path.name, "sha256": _sha256(output_path)}
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report
