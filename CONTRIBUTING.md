# Contributing to PloidyPatch

Thank you for helping improve duplication-aware plant annotation repair.

## Development setup

```bash
conda create -n ploidypatch-dev -c conda-forge python=3.11
conda activate ploidypatch-dev
python -m pip install -e ".[test]"
python -m pytest -q
```

Before opening a pull request:

1. add a focused regression test for each behavior change;
2. keep candidate ranking separate from human approval;
3. preserve non-overwrite, provenance, and checksum invariants;
4. run the complete test suite and the five-minute example;
5. update user-facing documentation and `CHANGELOG.md` when appropriate; and
6. do not commit private data, unpublished access credentials, or large
   third-party reference files.

## Pull requests

Use a small, reviewable branch and explain the biological assumption, software
contract, failure behavior, and verification performed. Changes to scientific
thresholds must include a label-free rationale and must not reuse revealed
holdout labels for tuning.

## Reporting problems

Use the issue templates for reproducible bugs and feature proposals. Security
issues should follow [SECURITY.md](SECURITY.md). General questions can be sent
to [litaishan@caf.ac.cn](mailto:litaishan@caf.ac.cn).
