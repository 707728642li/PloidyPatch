# PloidyPatch core-container smoke evidence v0.1

This compact package is the direct output of
`scripts/smoke_container_v0.1.sh` on the Linux `bioserver` host. The image was
built from `containers/Dockerfile` and the allowlisted `.dockerignore` context.

Source implementation hashes used by the successful build:

- `containers/Dockerfile`: SHA-256
  `147eb7873ae25ea094c5c9c1ee70d9a3a4ff3a17860f612298b1f8310a3a175c`;
- `scripts/smoke_container_v0.1.sh`: SHA-256
  `b8400509a622dccde6ed23c9fe7c31d94886eea64154231a466143ff88bbcd7a`;
- `.dockerignore`: SHA-256
  `1ed373ec8c5d2c397b5e44e838395161cbd540208a33a0a88b723f0caa5f7287`.

The successful image is
`sha256:d4457eb7c89e25d674bb53f26b78f483ad072029e26c1dc47431737ab811bbe7`
(133,761,776 bytes, Linux amd64). It runs as `10001:10001`, has the
`ploidypatch` entry point, and reports development version `0.1.0a0`.

`upstream_SHA256SUMS` is the unmodified manifest written on the server. It
binds the five raw output files copied here. `package_SHA256SUMS` additionally
binds this README and the upstream manifest.

Two execution defects were retained during development rather than hidden:

1. Docker build attempt 1 generated the wheel but renamed it to an invalid
   wheel filename before installation; pip rejected it.
2. Container smoke attempt 1 mounted the example at a path too shallow for the
   example runner's source-tree detection and failed before invoking the CLI.

The final Dockerfile preserves the standards-compliant wheel filename and
precreates `/work/examples`; the smoke script mounts the example at
`/work/examples/minimal_reviewed_patch`. The final run used a read-only root
filesystem, disabled networking, mounted inputs read-only, and supplied only
an ephemeral `/tmp` for outputs.

`container_smoke.json` records seven passing checks, including two explicitly
accepted additions, `automatic_approval=false`, and matching source/reverted
SHA-256 values. This evidence validates a local core image; it is not evidence
of publication to a container registry.
