"""Fail-closed verification for atomically published method outputs."""
from __future__ import annotations

import csv
from pathlib import Path, PurePosixPath
import re

from .artifact_manifest import sha256_file


SHA256 = re.compile(r"[0-9a-f]{64}")
BLIND_NAMESPACE_PROJECT_ROOT = PurePosixPath("/run/blind-run/project")
PUBLISHED_ROLES = {
    "gemoma": frozenset(
        {
            "final_annotation.gff",
            "unfiltered_predictions_from_species_0.gff",
            "reference_gene_table.tabular",
            "predicted_proteins.fasta",
            "protocol_GeMoMaPipeline.txt",
        }
    ),
    "lifton": frozenset({"lifton.gff3"}),
}


def verify_published_method_output(root: Path, expected_roles: set[str] | frozenset[str]) -> None:
    """Verify exact roles, bytes, hashes, and pre-rename path lineage."""
    root = root.resolve()
    manifest = root / "output_manifest.tsv"
    if (
        not root.is_dir()
        or root.is_symlink()
        or not manifest.is_file()
        or manifest.is_symlink()
    ):
        raise ValueError("Published method output or manifest is missing")
    with manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["role", "bytes", "sha256", "path"]:
            raise ValueError("Published method output manifest fields differ")
        rows = list(reader)
    roles = [row["role"] for row in rows]
    if set(roles) != set(expected_roles) or len(roles) != len(expected_roles):
        raise ValueError("Published method output roles differ")
    working = Path(str(root) + ".working")
    namespace_working: PurePosixPath | None = None
    project_indices = [index for index, part in enumerate(root.parts) if part == "project"]
    if project_indices:
        project_index = project_indices[-1]
        project_root = Path(*root.parts[: project_index + 1])
        relative_root = root.relative_to(project_root)
        namespace_working = (
            BLIND_NAMESPACE_PROJECT_ROOT
            / PurePosixPath(*relative_root.parent.parts)
            / f"{relative_root.name}.working"
        )
    for row in rows:
        role = row["role"]
        if not role or Path(role).name != role or role in {".", ".."}:
            raise ValueError(f"Unsafe published output role: {role}")
        recorded = row["path"]
        recorded_path = Path(recorded)
        host_lineage = (
            recorded_path.is_absolute()
            and recorded_path == working / "upstream" / role
        )
        recorded_posix = PurePosixPath(recorded)
        namespace_lineage = (
            namespace_working is not None
            and recorded_posix.is_absolute()
            and recorded_posix == namespace_working / "upstream" / role
        )
        if not host_lineage and not namespace_lineage:
            raise ValueError(f"Published output working-path lineage differs: {role}")
        path = root / "upstream" / role
        if (
            not path.is_file()
            or path.is_symlink()
            or not row["bytes"].isdigit()
            or path.stat().st_size != int(row["bytes"])
            or SHA256.fullmatch(row["sha256"]) is None
            or sha256_file(path) != row["sha256"]
        ):
            raise ValueError(f"Published method output differs: {role}")
