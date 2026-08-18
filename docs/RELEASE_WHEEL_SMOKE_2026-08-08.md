# Development wheel smoke test (2026-08-08)

Source commit: `b59f024af0a06ed68612fe6364f8a1a5bbfcf56a`.

This is a packaging smoke test for the development source tree.  It is not a
stable release, container, Bioconda package, or archived software artifact.

## Build

The wheel was built locally with the project development conda prefix using
PEP 517 isolation:

```powershell
conda run -p .\envs\ploidypatch-dev --no-capture-output `
  python -m pip wheel --no-deps `
  --wheel-dir .\work\release_smoke_b59f024\wheel .
```

Result:

- file: `ploidypatch-0.1.0a0-py3-none-any.whl`;
- bytes: 259,393;
- SHA-256: `9a32e8b0671026a6fe4413baf35d2b1fe2eb5d9884b56435b1e5296b242c3d72`.

## Clean installation

A new project-local conda prefix was created with Python 3.11.15 and pip.  The
wheel was installed with `--no-deps`, ensuring the core package did not rely on
an undeclared runtime package from the development environment.

The following checks passed:

```text
python -m pip check
No broken requirements found.

ploidypatch --version
0.1.0a0
```

`ploidypatch --help` successfully exposed all seven public command families:
`audit`, `benchmark`, `graph`, `patch`, `evidence`, `normalize`, and `baseline`.

## Boundary

This smoke test establishes that the pure-Python core can be built, installed,
imported, and launched from a clean prefix.  It does not bundle or validate
external miniprot, GeMoMa, LiftOn, SynGAP, WGDI, DIAMOND, or Pychopper
executables.  Those remain separately locked workflow dependencies.  A public
release still requires an explicit license decision, contributor/citation
metadata, archived artifacts, a small end-to-end example, container/Bioconda
packaging, and the frozen scaling result.

## Submission-tree revalidation (2026-08-10)

The post-Coffea manuscript tree at parent commit
`2d95a4b95a1d956fe26dbd2f7f48cf1dd66cb339` was rebuilt with the same isolated
PEP 517 wheel command. The resulting development wheel was:

- file: `ploidypatch-0.1.0a0-py3-none-any.whl`;
- bytes: 349,616;
- SHA-256: `b246f5a301365d428e2434c1b85dce19e4f30bce781eb037e496c129cf910bf0`.

It was force-installed with `--no-deps` into the project-local
`envs/ploidypatch-release-smoke` prefix. `pip check`, `python -m
ploidypatch.cli --version`, and the checksum-bound minimal reviewed-patch
example all passed; the example again emitted two explicitly reviewed
additions, `automatic_approval=false`, and byte-identical reversion. The full
source test suite passed 641 tests with nine platform-specific Windows skips.

The earlier boundary is now narrower: the minimal example and frozen scaling
study are complete. License choice, contributor/citation metadata, and public
archive/DOI creation remain release-administration actions that require the
project owners' final identities and repository destination. Container or
Bioconda packaging can follow publication and is not an evidence prerequisite.

## Public release-surface revalidation (2026-08-10)

The user-facing release candidate at source commit
`be06422217712ba3fcb0a64c3efe17e3f207f9da` was rebuilt with Python 3.11.15
and PyPA build 1.5.0. The build used PEP 517 isolation and produced both
distributions from the source tree:

- wheel: `ploidypatch-0.1.0a0-py3-none-any.whl`, 342,496 bytes,
  SHA-256 `28b8c8d1813acbde960422430194691209137a508f96ce9a8b53f52029b96418`;
- source distribution: `ploidypatch-0.1.0a0.tar.gz`, 1,232,403 bytes,
  SHA-256 `30c5efa53d6996198bf0ceab6a71a638f140ec4e158d026e3e8acabc083c9eb1`.

The source archive was inspected before installation and contained the
minimal reviewed-patch runner, user guide, reproducibility guide,
machine-generated CLI command inventory, Coffea holdout contract, and Coffea
blind pipeline entry point. This closes the earlier gap in which JSON CLI
documentation was tracked but not selected by `MANIFEST.in`; CI now checks the
same exact paths.

The wheel was installed with `--no-deps` into a new local Python 3.11 virtual
environment with no source-tree `PYTHONPATH`. The following all passed:

- `python -m pip check`;
- `python -m ploidypatch.cli --version` (`0.1.0a0`);
- top-level help containing all seven command families;
- the copied, installed-wheel minimal reviewed-patch example;
- two explicitly accepted additions, `automatic_approval=false`, and
  byte-identical source/reverted SHA-256; and
- the complete source suite: 659 passed and nine platform-specific Windows
  skips.

The public release surface now has a concise README, a task-oriented user
guide, a curated CLI reference, a reproducibility/custody guide, a research
status map, and an 86-leaf command inventory exported from the production
parser. Formal publication still requires owner-supplied authorship, license,
repository URL and archive DOI; none has been inferred.

## Manuscript-aligned release candidate (2026-08-10)

After freezing the v0.5 manuscript and v0.4 Genome Biology submission package,
the release surface was rebuilt from clean source commit
`34988bdfa15a5f3bd107cf11329f26e842042e66`. The fail-closed v0.2 release
builder used Python 3.11.15 and isolated PEP 517 construction, audited both
archive member universes, and produced:

- wheel: `ploidypatch-0.1.0a0-py3-none-any.whl`, 360,670 bytes,
  SHA-256 `a5287972f6bf4b60d044b466b565cab3ad0b30063189180dbbb607c6fb027f53`;
- source distribution: `ploidypatch-0.1.0a0.tar.gz`, 1,318,661 bytes,
  SHA-256 `3ad3e3a2ef65e5ed33c120113d56d27564790129c6df124c4292e702b4cfd916`.

The 78-member wheel has the expected console entry point and no mandatory
third-party runtime dependency. The 768-member source distribution contains
the user and reproducibility guides, generated CLI inventory, minimal reviewed
patch example, container definition and smoke runner. A new virtual environment
with no source-tree `PYTHONPATH` installed the wheel using `--no-deps`; `pip
check`, version query, seven-family help inventory and the copied example all
passed. The example accepted exactly two explicit review decisions, retained
`automatic_approval=false`, and reverted to the source GFF3 byte-for-byte.

Compact evidence is frozen under `docs/release_evidence/wheel_core_v0.2/`.
Its exact-universe `SHA256SUMS` has SHA-256
`8160101578dd8fa96863b33e13b626ebc2e89e40f523f4aa6923f41d3c3b1e97`.
After integrating that evidence and the current owner-handoff check, the
complete source suite passed 731 tests with thirteen platform-specific Windows
skips.
The distributions themselves remain in ignored local build storage; they are
development candidates, not a licensed public release. Author identities,
license choice, repository URL and archive DOI remain deliberately unresolved.

## Byte-reproducible distribution candidate (2026-08-10)

Release builder v0.3 closed the remaining archive-reproducibility gap at clean
source commit `e4d32e60a5ea16b09e565480cfb25a32cea6bde2`. It derives
`SOURCE_DATE_EPOCH` from that commit, performs two independent isolated PEP 517
builds, and canonicalizes sdist tar/gzip metadata without changing member
content. Both independent wheels were byte-identical and both independently
canonicalized sdists were byte-identical:

- wheel: 360,677 bytes, 78 members, SHA-256
  `b14592beea9735091e26eb7184c8a38a846616398ac329eb2b20da92eb66a246`;
- sdist: 1,291,167 bytes, 782 members, SHA-256
  `d4be89e48d7768a9449fe14aa0cf0ac0be32856dd5184dff7b26e35d32c0cf01`.

The wheel and sdist were then installed separately into clean virtual
environments with no source-tree `PYTHONPATH`. Both passed `pip check`, version
and seven-family CLI inspection, and the complete reviewed-patch example. The
five stable example outputs were byte-identical; path-normalized command logs
were semantically identical; both runs accepted exactly two explicit review
decisions, retained `automatic_approval=false`, and reverted byte-for-byte.

Compact evidence is frozen under `docs/release_evidence/wheel_core_v0.3/`.
Its exact-universe `SHA256SUMS` has SHA-256
`45d53a8d0fe6b803e52918430c2dfeacc187d8545fd81c8423eacde430b5e4f0`.
CI now executes the same two-build, two-install gate. This remains a development
candidate: reproducibility does not supply authorship, select a license, create
a public repository or mint an archive DOI.

## Core-container candidate (2026-08-10)

A digest-pinned, multi-stage, non-root core container was subsequently added at
commit `f75c35a343ad09418ac75948e9402f0d49a166e8`. The Linux amd64 image was
built on `bioserver` and passed the same reversible example with a read-only
root filesystem, no runtime network, read-only inputs and runtime UID/GID
`10001:10001`. The resulting image ID was
`sha256:d4457eb7c89e25d674bb53f26b78f483ad072029e26c1dc47431737ab811bbe7`
(133,761,776 bytes).

The container keeps the pure-Python core boundary: no external predictor,
aligner or synteny executable is bundled. Both retained failed attempts and
the final seven passing checks are documented in
`docs/CONTAINER_SMOKE_2026-08-10.md`; raw compact evidence is tracked under
`docs/release_evidence/container_core_v0.1/`. This is a verified local image
candidate, not a registry publication or Bioconda release.
