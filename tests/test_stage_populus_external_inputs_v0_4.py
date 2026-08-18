from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


MODULE_PATH = (
    Path(__file__).parents[1]
    / "scripts"
    / "stage_populus_external_inputs_v0.4.py"
)
SPEC = importlib.util.spec_from_file_location("populus_stage", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


@pytest.mark.parametrize(
    ("role", "artifact", "expected_prefix"),
    (
        ("target", "genome", "shared_target/Populus/"),
        ("target", "gff3", "evaluator_only/target_complete/Populus/"),
        ("target", "protein", "evaluator_only/target_complete/Populus/"),
        ("candidate_reference", "genome", "candidate_only/Populus/"),
        (
            "evaluator_reference",
            "genome",
            "evaluator_only/truth_references/Populus/",
        ),
    ),
)
def test_destination_enforces_role_separation(
    role: str, artifact: str, expected_prefix: str
) -> None:
    destination = MODULE.destination_for(
        {
            "role": role,
            "species_id": "Populus",
            "artifact": artifact,
            "source_path": "/source/file.gz",
        }
    )
    assert destination.as_posix().startswith(expected_prefix)


def test_destination_rejects_unknown_role() -> None:
    with pytest.raises(ValueError, match="Unknown"):
        MODULE.destination_for(
            {
                "role": "mixed_candidate_truth",
                "species_id": "x",
                "artifact": "genome",
                "source_path": "/source/file.gz",
            }
        )
