from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from . import __version__
from .io import open_text
from .perturb import _file_sha256, _text_sha256


PATCH_SCHEMA_VERSION = "ploidypatch.annotation_patch.v1"


def _json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def _write_text_exclusive(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(text)


def _line_sha(raw_line: str) -> str:
    return hashlib.sha256(raw_line.encode("utf-8")).hexdigest()


def _validate_replacement_line(raw_line: Any, operation_number: int) -> str:
    if not isinstance(raw_line, str):
        raise ValueError(
            f"Replacement in operation {operation_number} is not a string"
        )
    body = raw_line[:-2] if raw_line.endswith("\r\n") else raw_line.rstrip("\n")
    if "\n" in body or "\r" in body:
        raise ValueError(
            f"Replacement in operation {operation_number} contains multiple lines"
        )
    return raw_line


def _load_edit_spec(path: str | Path) -> tuple[list[str], list[dict[str, Any]]]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    event_ids = [str(value) for value in payload.get("event_ids", [])]
    if "operations" in payload:
        return event_ids, list(payload["operations"])
    if "events" in payload:
        operations: list[dict[str, Any]] = []
        for event in payload["events"]:
            if event.get("event_id"):
                event_ids.append(str(event["event_id"]))
            for edit in event.get("line_edits", []):
                operations.append(
                    {
                        "source_line_number": edit["source_line_number"],
                        "source_raw_line": edit.get("source_raw_line"),
                        "replacement_lines": [
                            replacement["raw_line"]
                            for replacement in edit.get("perturbed_lines", [])
                        ],
                    }
                )
        return event_ids, operations
    raise ValueError("Edit specification requires operations or events")


def create_annotation_patch(
    source_gff_path: str | Path,
    edits_json_path: str | Path,
    output_patch_path: str | Path,
) -> dict[str, Any]:
    """Freeze explicit line edits into an independently reversible patch."""

    output = Path(output_patch_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite annotation patch: {output}")
    with open_text(source_gff_path) as handle:
        source_lines = list(handle)
    event_ids, requested = _load_edit_spec(edits_json_path)
    if not requested:
        raise ValueError("Edit specification contains no operations")

    operations: list[dict[str, Any]] = []
    seen_lines: set[int] = set()
    for operation_number, operation in enumerate(requested, start=1):
        try:
            line_number = int(operation["source_line_number"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid source_line_number in operation {operation_number}"
            ) from exc
        if not 1 <= line_number <= len(source_lines):
            raise ValueError(
                f"Operation {operation_number} source line is outside the source"
            )
        if line_number in seen_lines:
            raise ValueError(f"Duplicate source line operation: {line_number}")
        seen_lines.add(line_number)
        source_raw_line = source_lines[line_number - 1]
        declared_source = operation.get("source_raw_line")
        if declared_source is not None and declared_source != source_raw_line:
            raise ValueError(
                f"Declared source bytes disagree at source line {line_number}"
            )
        replacements = tuple(
            _validate_replacement_line(raw_line, operation_number)
            for raw_line in operation.get("replacement_lines", [])
        )
        operations.append(
            {
                "source_line_number": line_number,
                "source_raw_line": source_raw_line,
                "source_line_sha256": _line_sha(source_raw_line),
                "replacement_lines": [
                    {"raw_line": raw_line, "line_sha256": _line_sha(raw_line)}
                    for raw_line in replacements
                ],
            }
        )
    operations.sort(key=lambda item: item["source_line_number"])
    replacements_by_line = {
        operation["source_line_number"]: tuple(
            replacement["raw_line"]
            for replacement in operation["replacement_lines"]
        )
        for operation in operations
    }
    patched_lines: list[str] = []
    for line_number, raw_line in enumerate(source_lines, start=1):
        patched_lines.extend(replacements_by_line.get(line_number, (raw_line,)))

    patch: dict[str, Any] = {
        "schema_version": PATCH_SCHEMA_VERSION,
        "generator": {"name": "PloidyPatch", "version": __version__},
        "source": {
            "file_name": Path(source_gff_path).name,
            "file_sha256": _file_sha256(source_gff_path),
            "text_sha256": _text_sha256(source_lines),
            "line_count": len(source_lines),
        },
        "edit_spec": {
            "file_name": Path(edits_json_path).name,
            "sha256": _file_sha256(edits_json_path),
        },
        "event_ids": sorted(set(event_ids)),
        "operations": operations,
        "patched": {
            "text_sha256": _text_sha256(patched_lines),
            "line_count": len(patched_lines),
        },
    }
    _write_text_exclusive(output, _json_text(patch))
    return patch


def _load_patch(path: str | Path) -> dict[str, Any]:
    patch = json.loads(Path(path).read_text(encoding="utf-8"))
    if patch.get("schema_version") != PATCH_SCHEMA_VERSION:
        raise ValueError("Unsupported or missing annotation patch schema")
    return patch


def _validated_operations(
    patch: dict[str, Any], source_line_count: int
) -> dict[int, tuple[str, tuple[str, ...]]]:
    operations: dict[int, tuple[str, tuple[str, ...]]] = {}
    for operation in patch.get("operations", []):
        line_number = int(operation["source_line_number"])
        if not 1 <= line_number <= source_line_count or line_number in operations:
            raise ValueError(f"Invalid or duplicate patch source line: {line_number}")
        source_raw_line = operation["source_raw_line"]
        if _line_sha(source_raw_line) != operation["source_line_sha256"]:
            raise ValueError(f"Patch source-line checksum failed at {line_number}")
        replacements: list[str] = []
        for replacement in operation.get("replacement_lines", []):
            raw_line = replacement["raw_line"]
            if _line_sha(raw_line) != replacement["line_sha256"]:
                raise ValueError(
                    f"Patch replacement checksum failed at source line {line_number}"
                )
            replacements.append(raw_line)
        operations[line_number] = (source_raw_line, tuple(replacements))
    if not operations:
        raise ValueError("Annotation patch contains no operations")
    return operations


def apply_annotation_patch(
    source_gff_path: str | Path,
    patch_path: str | Path,
    output_gff_path: str | Path,
) -> dict[str, Any]:
    """Apply a patch without modifying the source annotation."""

    output = Path(output_gff_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite patched GFF3: {output}")
    patch = _load_patch(patch_path)
    with open_text(source_gff_path) as handle:
        source_lines = list(handle)
    if len(source_lines) != int(patch["source"]["line_count"]):
        raise ValueError("Source line count does not match annotation patch")
    if _text_sha256(source_lines) != patch["source"]["text_sha256"]:
        raise ValueError("Source text checksum does not match annotation patch")
    operations = _validated_operations(patch, len(source_lines))
    patched_lines: list[str] = []
    for line_number, raw_line in enumerate(source_lines, start=1):
        operation = operations.get(line_number)
        if operation is None:
            patched_lines.append(raw_line)
            continue
        expected_source, replacements = operation
        if raw_line != expected_source:
            raise ValueError(f"Source bytes disagree at patch line {line_number}")
        patched_lines.extend(replacements)
    observed_sha = _text_sha256(patched_lines)
    if observed_sha != patch["patched"]["text_sha256"]:
        raise ValueError("Applied output checksum does not match annotation patch")
    _write_text_exclusive(output, "".join(patched_lines))
    return {
        "output_gff": str(output),
        "text_sha256": observed_sha,
        "operations": len(operations),
        "event_ids": patch.get("event_ids", []),
    }


def revert_annotation_patch(
    patched_gff_path: str | Path,
    patch_path: str | Path,
    output_gff_path: str | Path,
) -> dict[str, Any]:
    """Revert a patched annotation and require the exact original checksum."""

    output = Path(output_gff_path)
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite reverted GFF3: {output}")
    patch = _load_patch(patch_path)
    with open_text(patched_gff_path) as handle:
        patched_lines = list(handle)
    if len(patched_lines) != int(patch["patched"]["line_count"]):
        raise ValueError("Patched line count does not match annotation patch")
    if _text_sha256(patched_lines) != patch["patched"]["text_sha256"]:
        raise ValueError("Patched text checksum does not match annotation patch")
    source_line_count = int(patch["source"]["line_count"])
    operations = _validated_operations(patch, source_line_count)

    restored_lines: list[str] = []
    patched_index = 0
    for line_number in range(1, source_line_count + 1):
        operation = operations.get(line_number)
        if operation is None:
            if patched_index >= len(patched_lines):
                raise ValueError("Patched GFF3 ended before reversion completed")
            restored_lines.append(patched_lines[patched_index])
            patched_index += 1
            continue
        source_raw_line, replacements = operation
        observed = tuple(
            patched_lines[patched_index : patched_index + len(replacements)]
        )
        if observed != replacements:
            raise ValueError(
                f"Patched replacement bytes disagree at source line {line_number}"
            )
        patched_index += len(replacements)
        restored_lines.append(source_raw_line)
    if patched_index != len(patched_lines):
        raise ValueError("Patched GFF3 contains unexpected trailing lines")
    observed_sha = _text_sha256(restored_lines)
    if observed_sha != patch["source"]["text_sha256"]:
        raise ValueError("Reverted output checksum does not match original source")
    _write_text_exclusive(output, "".join(restored_lines))
    return {
        "output_gff": str(output),
        "text_sha256": observed_sha,
        "operations": len(operations),
        "event_ids": patch.get("event_ids", []),
    }
