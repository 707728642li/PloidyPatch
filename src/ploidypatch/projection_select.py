from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from . import __version__
from .baseline import _file_sha256
from .gff import parse_attributes
from .io import open_text


PROJECTION_SELECTION_SCHEMA_VERSION = "ploidypatch.projection_selection.v1"


def _read_source_groups(
    path: str | Path | None, observed_sources: set[str]
) -> tuple[dict[str, str], str, dict[str, str] | None]:
    if path is None:
        return (
            {source: source for source in observed_sources},
            "source",
            None,
        )
    mapping: dict[str, str] = {}
    source_path = Path(path)
    with source_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not {
            "source",
            "support_group",
        } <= set(reader.fieldnames):
            raise ValueError(
                "Source-group map requires source and support_group columns"
            )
        for line_number, row in enumerate(reader, start=2):
            source = row["source"]
            group = row["support_group"]
            if not source or not group or source in mapping:
                raise ValueError(
                    f"Empty or duplicate source-group mapping at line {line_number}"
                )
            mapping[source] = group
    missing = observed_sources - set(mapping)
    extra = set(mapping) - observed_sources
    if missing or extra:
        raise ValueError(
            "Source-group map must match projection sources exactly; "
            f"missing={sorted(missing)}, extra={sorted(extra)}"
        )
    return mapping, "explicit_group", {
        "file_name": source_path.name,
        "sha256": _file_sha256(source_path),
    }


def select_projection_support_models(
    *,
    candidate_gff_path: str | Path,
    projection_support_tsv_path: str | Path,
    output_gff_path: str | Path,
    selection_tsv_path: str | Path,
    min_support_group_count: int = 2,
    source_group_map_path: str | Path | None = None,
) -> dict[str, Any]:
    """Retain appended projection models supported by independent groups."""

    if min_support_group_count < 1:
        raise ValueError("min_support_group_count must be positive")
    output = Path(output_gff_path)
    selection = Path(selection_tsv_path)
    manifest_path = Path(str(output) + ".manifest.json")
    existing = [path for path in (output, selection, manifest_path) if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite projection-selection artifact(s): "
            + ", ".join(str(path) for path in existing)
        )

    support_rows: list[dict[str, str]] = []
    support_by_model: dict[str, tuple[str, ...]] = {}
    observed_sources: set[str] = set()
    with Path(projection_support_tsv_path).open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"model_id", "support_source_count", "support_sources"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise ValueError(
                "Projection support is missing model_id/support source columns"
            )
        for line_number, row in enumerate(reader, start=2):
            model_id = row["model_id"]
            sources = tuple(
                sorted(source for source in row["support_sources"].split(",") if source)
            )
            if not model_id or model_id in support_by_model:
                raise ValueError(
                    f"Empty or duplicate support model_id at line {line_number}"
                )
            try:
                declared_count = int(row["support_source_count"])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid support source count at line {line_number}"
                ) from exc
            if declared_count != len(sources) or len(set(sources)) != len(sources):
                raise ValueError(
                    f"Support source count mismatch at line {line_number}"
                )
            support_by_model[model_id] = sources
            observed_sources.update(sources)
            support_rows.append(row)
    if not support_rows:
        raise ValueError("Projection support table contains no models")

    source_groups, support_unit, group_manifest = _read_source_groups(
        source_group_map_path, observed_sources
    )
    selected_models = {
        model_id
        for model_id, sources in support_by_model.items()
        if len({source_groups[source] for source in sources})
        >= min_support_group_count
    }

    baseline_lines: dict[str, list[str]] = {}
    transcript_models: dict[str, str] = {}
    other_lines: list[str] = []
    with open_text(candidate_gff_path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            stripped = raw_line.rstrip("\r\n")
            if not stripped or stripped.startswith("#"):
                other_lines.append(raw_line)
                continue
            fields = stripped.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Malformed candidate GFF3 line {line_number}: expected 9 fields"
                )
            if fields[1] != "PloidyPatchBaseline":
                other_lines.append(raw_line)
                continue
            attributes, malformed = parse_attributes(fields[8])
            if malformed:
                raise ValueError(
                    f"Malformed baseline attributes at line {line_number}"
                )
            model_id = attributes.get("miniprot_model", "")
            if model_id:
                feature_id = attributes.get("ID", "")
                if fields[2] in {"mRNA", "transcript"} and feature_id:
                    transcript_models[feature_id] = model_id
            else:
                parents = [
                    value
                    for value in attributes.get("Parent", "").split(",")
                    if value
                ]
                parent_models = {
                    transcript_models[parent]
                    for parent in parents
                    if parent in transcript_models
                }
                if len(parent_models) != 1:
                    raise ValueError(
                        f"Cannot resolve one baseline model at line {line_number}"
                    )
                model_id = next(iter(parent_models))
            baseline_lines.setdefault(model_id, []).append(raw_line)

    missing_support = set(support_by_model) - set(baseline_lines)
    extra_baseline = set(baseline_lines) - set(support_by_model)
    if missing_support or extra_baseline:
        raise ValueError(
            "Candidate/support model mismatch; "
            f"missing_in_gff={sorted(missing_support)[:10]}, "
            f"missing_in_support={sorted(extra_baseline)[:10]}"
        )

    # The adapter writes every baseline model after the unchanged annotation.
    # Retain all non-baseline lines byte-for-byte and append selected complete
    # model hierarchies in their original model order.
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="") as handle:
        handle.writelines(other_lines)
        for model_id, lines in baseline_lines.items():
            if model_id in selected_models:
                handle.writelines(lines)

    selection.parent.mkdir(parents=True, exist_ok=True)
    selection_fields = (
        "model_id",
        "support_source_count",
        "support_sources",
        "support_group_count",
        "support_groups",
        "status",
        "reason",
    )
    selected_line_count = 0
    with selection.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=selection_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in support_rows:
            model_id = row["model_id"]
            sources = support_by_model[model_id]
            groups = tuple(sorted({source_groups[source] for source in sources}))
            retained = model_id in selected_models
            if retained:
                selected_line_count += len(baseline_lines[model_id])
            writer.writerow(
                {
                    "model_id": model_id,
                    "support_source_count": len(sources),
                    "support_sources": ",".join(sources),
                    "support_group_count": len(groups),
                    "support_groups": ",".join(groups),
                    "status": "accepted" if retained else "rejected",
                    "reason": (
                        "independent_support_pass"
                        if retained
                        else "support_below_threshold"
                    ),
                }
            )

    manifest: dict[str, Any] = {
        "schema_version": PROJECTION_SELECTION_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "inputs": {
            "candidate_gff": {
                "file_name": Path(candidate_gff_path).name,
                "sha256": _file_sha256(candidate_gff_path),
            },
            "projection_support": {
                "file_name": Path(projection_support_tsv_path).name,
                "sha256": _file_sha256(projection_support_tsv_path),
                "rows": len(support_rows),
            },
            "source_group_map": group_manifest,
        },
        "parameters": {
            "min_support_group_count": min_support_group_count,
            "support_unit": support_unit,
            "baseline_source": "PloidyPatchBaseline",
        },
        "counts": {
            "input_models": len(support_rows),
            "selected_models": len(selected_models),
            "rejected_models": len(support_rows) - len(selected_models),
            "selected_feature_lines": selected_line_count,
        },
        "outputs": {
            "candidate_gff": {
                "file_name": output.name,
                "sha256": _file_sha256(output),
            },
            "selection": {
                "file_name": selection.name,
                "sha256": _file_sha256(selection),
                "rows": len(support_rows),
            },
        },
    }
    with manifest_path.open("x", encoding="utf-8", newline="") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return manifest
