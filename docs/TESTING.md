# Testing PloidyPatch

PloidyPatch separates the reproducible public software suite from private
prepublication research material. Manuscripts, paper figures, biological result
trees, reviewer packets, and their internal archival tests remain in the local
scientific workspace and are not tracked by the software repository.

## Public suite

From a clean source checkout:

```bash
python -m pip install -e ".[test]"
python scripts/run_public_test_suite_v1.py
python examples/minimal_reviewed_patch/run_example.py \
  --output-dir work/minimal-reviewed-patch
```

The public runner executes every test shipped in the checkout: package, CLI,
normalization, graph, evidence, patching, release-build, report, and
scientific-core tests. It does not maintain a hidden exclusion list. Private
research tests are absent from the public tree together with the unpublished
artifacts they audit.

GitHub Actions runs this suite on Python 3.11, 3.12 and 3.13. The workflow also
checks generated CLI documentation under the canonical Python 3.11 maintenance
environment, builds distributions reproducibly, installs the wheel without the
source tree, and executes the checksum-bound example. Python 3.12 and 3.13 still
run the public suite; the documentation check is kept on 3.11 because argparse
help wrapping differs between Python minor versions without changing the CLI.

## Private archival suite

Maintainers run the private archival suite from the local scientific workspace,
where the checksum-bound manuscript and result trees are present:

```bash
python -m pytest -q
```

The distinction is deliberate: public CI validates the distributed software,
whereas the private suite validates unpublished evidence and submission
artifacts without placing those files in repository history.
