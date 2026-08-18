from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "run_helianthus_fixed_backbone_preflight_v0.6.sh"


def test_runner_is_target_only_frozen_and_candidate_free() -> None:
    text = SCRIPT.read_text(encoding="utf-8")
    assert "ploidypatch_v0.6_155f72a" in text
    assert "code_commit\\t155f72a" in text
    assert "runner_sha256" in text
    assert "wgdi_version\\t0.75" in text
    assert "candidate_access\\tfalse" in text
    assert "truth_or_label_access\\tfalse" in text
    assert "--max-target-seqs 20" in text
    assert "--min-block-pairs 20" in text
    assert "build-fixed-target-backbone" in text
    assert "candidate.gff" not in text
    assert "label" not in text.replace("truth_or_label_access", "")
    assert "evaluator" not in text
    assert "find . -type f ! -name SHA256SUMS -printf '%P\\0'" in text
