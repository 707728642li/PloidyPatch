# PloidyPatch core-container smoke audit (2026-08-10)

## Frozen implementation

The tested implementation is source commit
`f75c35a343ad09418ac75948e9402f0d49a166e8`.

- Base image: `python:3.11.15-slim-bookworm` at OCI index digest
  `sha256:d29f48a31a8b408ed19272ca1e7b10ebae13b240a27e862d3d4217c528e2e0c3`.
- `containers/Dockerfile`: SHA-256
  `147eb7873ae25ea094c5c9c1ee70d9a3a4ff3a17860f612298b1f8310a3a175c`.
- `scripts/smoke_container_v0.1.sh`: SHA-256
  `b8400509a622dccde6ed23c9fe7c31d94886eea64154231a466143ff88bbcd7a`.
- `.dockerignore`: SHA-256
  `1ed373ec8c5d2c397b5e44e838395161cbd540208a33a0a88b723f0caa5f7287`.

The Docker build was executed on the Linux amd64 `bioserver` host under the
project-only working directory
`/data/codexli/projects/PloidyPatch/work/container_smoke_20260810_v0.1`.
BuildKit was enabled. Both wheel construction and runtime wheel installation
ran with `--network=none` and `--no-deps`; wheel construction also used
`--no-build-isolation` against setuptools 79.0.1 already present in the pinned
base image.

## Retained failed attempts

The audit keeps two fail-closed engineering defects because both would have
invalidated a public container smoke claim.

1. Build attempt 1 successfully created the wheel but copied it to
   `/tmp/ploidypatch.whl`. Pip correctly rejected that filename as an invalid
   wheel filename. The preserved build log has SHA-256
   `5d284771e085d8d133eec422502138bf0f6eca39fe7f7a5da9718c849023c23d`.
2. Build attempt 2 succeeded after preserving the standards-compliant wheel
   filename, but smoke attempt 1 mounted the example at `/example`. The example
   runner deliberately checks for an adjacent source tree and requires a
   deeper path; it failed with `IndexError` before any PloidyPatch command ran.
   Its partial `.working` output remains on the server and was never promoted.

The final image precreates `/work/examples`, and the final smoke script mounts
the read-only example at `/work/examples/minimal_reviewed_patch`. At that
location the script can verify that `/work/src` is absent and therefore uses
the installed package rather than source checkout imports.

## Successful build and runtime result

Build attempt 3 succeeded. Its build log has SHA-256
`723dbd60d5dba70dd6b46a0ff59b4a419e3e5aa7eba76871ed709949a9f50fc7`.
The resulting local image was:

- image ID:
  `sha256:d4457eb7c89e25d674bb53f26b78f483ad072029e26c1dc47431737ab811bbe7`;
- image size: 133,761,776 bytes;
- platform: Linux amd64;
- runtime user and group: `10001:10001`;
- entry point: `ploidypatch`;
- reported version: `0.1.0a0`.

The smoke runner used a read-only root filesystem, `--network none`, a
read-only example mount, and a 16 MiB no-exec/nosuid tmpfs for ephemeral output.
All seven checks passed:

1. runtime user is `10001:10001`;
2. entry point is `ploidypatch`;
3. version is `0.1.0a0`;
4. the example contains exactly two accepted additions;
5. `automatic_approval=false`;
6. byte-identical reversion is true; and
7. source and reverted SHA-256 values are identical
   (`e2db45283602ab57a002c57f9b7627ef9fa79818e6cefd111584a3c84aa6131f`).

Top-level help also contained all seven command families. The final server
output contained exactly five raw files plus `SHA256SUMS`; every entry passed
`sha256sum -c`.

## Tracked evidence and boundary

The compact raw evidence is tracked under
`docs/release_evidence/container_core_v0.1/`. Its
`package_SHA256SUMS` has SHA-256
`240aa19b2f4a1deaad700600e5b22c17063dd93269bec141dd7d521bfb52e21a`.
The source distribution was rebuilt and confirmed to contain the Dockerfile,
`.dockerignore`, smoke runner, container guide and compact smoke evidence.

The source suite passed 666 tests with nine platform-specific Windows skips at
the implementation commit. This audit proves that a local core container can
be built and exercised under the declared constraints. It does not claim that
an image has been published to GHCR, Docker Hub or another registry. Registry
publication, license metadata and Bioconda packaging remain owner/release
actions.
