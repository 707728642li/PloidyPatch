"""Exact, auditable sequence-ID normalization for role-bound GFF3 inputs."""
from __future__ import annotations

from collections import Counter
import csv
import json
import os
from pathlib import Path
import tempfile
from typing import Any

from .artifact_manifest import sha256_file
from .io import open_text


SCHEMA_VERSION = "ploidypatch.exact_seqid_alias_rewrite.v1"


def _read_exact_table(path: str | Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    source = Path(path)
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked TSV: {source}")
    with source.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise ValueError(
                f"Exact TSV header differs for {source}: {reader.fieldnames!r}"
            )
        rows = list(reader)
    if not rows or any(not row[field] for row in rows for field in fields):
        raise ValueError(f"Exact TSV has empty data: {source}")
    return rows


def rewrite_gff3_seqids_exact(
    *,
    input_gff_path: str | Path,
    alias_table_path: str | Path,
    primary_seqid_table_path: str | Path,
    output_gff_path: str | Path,
    manifest_path: str | Path,
) -> dict[str, Any]:
    """Rewrite only explicitly mapped GFF3 seqids and drop all others.

    The alias target universe must equal the complete primary-seqid table.
    There is no prefix, substring, coordinate, order or fuzzy fallback.
    """

    source = Path(input_gff_path)
    output = Path(output_gff_path)
    manifest_output = Path(manifest_path)
    if not source.is_file() or source.is_symlink() or source.stat().st_size == 0:
        raise ValueError(f"Missing, empty or symlinked input GFF3: {source}")
    if output.resolve(strict=False) == source.resolve(strict=False):
        raise ValueError("Exact seqid rewrite cannot modify its input in place")
    if any(path.exists() or path.is_symlink() for path in (output, manifest_output)):
        raise FileExistsError("Refusing to overwrite exact seqid rewrite output")

    aliases = _read_exact_table(
        alias_table_path, ("gff_seqid", "genome_seqid", "chromosome_label")
    )
    primary = _read_exact_table(
        primary_seqid_table_path, ("seqid", "chromosome_label")
    )
    alias_sources = [row["gff_seqid"] for row in aliases]
    alias_targets = [row["genome_seqid"] for row in aliases]
    primary_seqids = [row["seqid"] for row in primary]
    if len(alias_sources) != len(set(alias_sources)):
        raise ValueError("Source GFF seqid alias is not unique")
    if len(alias_targets) != len(set(alias_targets)):
        raise ValueError("Target genome seqid alias is not unique")
    if len(primary_seqids) != len(set(primary_seqids)):
        raise ValueError("Primary genome seqid table is not unique")
    alias_labeled_targets = {
        (row["genome_seqid"], row["chromosome_label"]) for row in aliases
    }
    primary_labeled_targets = {
        (row["seqid"], row["chromosome_label"]) for row in primary
    }
    if alias_labeled_targets != primary_labeled_targets:
        raise ValueError("Alias target universe differs from primary seqid table")
    mapping = {row["gff_seqid"]: row["genome_seqid"] for row in aliases}

    output.parent.mkdir(parents=True, exist_ok=True)
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.working.", dir=output.parent
    )
    os.close(file_descriptor)
    temporary = Path(temporary_name)
    feature_counts: Counter[str] = Counter()
    source_feature_seqids: Counter[str] = Counter()
    output_feature_seqids: Counter[str] = Counter()
    comment_lines = 0
    skipped_nonprimary_features = 0
    skipped_nonprimary_regions = 0
    try:
        with open_text(source) as input_handle, temporary.open(
            "w", encoding="utf-8", newline=""
        ) as output_handle:
            for line_number, raw in enumerate(input_handle, start=1):
                line = raw.rstrip("\r\n")
                if not line:
                    output_handle.write("\n")
                    continue
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) != 4:
                        raise ValueError(
                            f"Malformed GFF3 sequence-region at line {line_number}"
                        )
                    mapped = mapping.get(parts[1])
                    if mapped is None:
                        skipped_nonprimary_regions += 1
                        continue
                    output_handle.write(
                        f"##sequence-region {mapped} {parts[2]} {parts[3]}\n"
                    )
                    comment_lines += 1
                    continue
                if line.startswith("#"):
                    output_handle.write(line + "\n")
                    comment_lines += 1
                    continue
                fields = line.split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"Malformed GFF3 line {line_number}: expected 9 fields"
                    )
                source_feature_seqids[fields[0]] += 1
                mapped = mapping.get(fields[0])
                if mapped is None:
                    skipped_nonprimary_features += 1
                    continue
                fields[0] = mapped
                feature_counts[fields[2]] += 1
                output_feature_seqids[mapped] += 1
                output_handle.write("\t".join(fields) + "\n")
            output_handle.flush()
            os.fsync(output_handle.fileno())

        missing_alias_sources = sorted(set(mapping) - set(source_feature_seqids))
        if missing_alias_sources:
            raise ValueError(
                "Alias source seqids have no GFF3 features: "
                + ",".join(missing_alias_sources)
            )
        if set(output_feature_seqids) != set(primary_seqids):
            raise ValueError("Normalized GFF3 does not cover the exact primary universe")
        os.replace(temporary, output)
        manifest = {
            "schema_version": SCHEMA_VERSION,
            "mapping_policy": "exact_bijective_alias_only_no_fuzzy_fallback",
            "primary_only": True,
            "counts": {
                "aliases": len(mapping),
                "output_primary_seqids": len(output_feature_seqids),
                "output_features": sum(feature_counts.values()),
                "skipped_nonprimary_features": skipped_nonprimary_features,
                "skipped_nonprimary_sequence_regions": skipped_nonprimary_regions,
                "preserved_comment_lines": comment_lines,
                "feature_types": dict(sorted(feature_counts.items())),
            },
            "inputs": {
                "gff3_sha256": sha256_file(source),
                "alias_table_sha256": sha256_file(alias_table_path),
                "primary_seqid_table_sha256": sha256_file(primary_seqid_table_path),
            },
            "outputs": {"gff3_sha256": sha256_file(output)},
            "coordinate_or_attribute_change": False,
            "truth_access": False,
        }
        manifest_output.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return manifest
    except BaseException:
        temporary.unlink(missing_ok=True)
        if output.exists() and not manifest_output.exists():
            output.unlink()
        raise
