# PloidyPatch core container

## Scope

The PloidyPatch container packages the same pure-Python reviewed-patch core as
the wheel. It does not bundle miniprot, GeMoMa, LiftOn, DIAMOND, WGDI, SynGAP,
Pychopper, aligners, or gene predictors. Formal upstream workflows retain
their separate, tool-specific conda environments.

The container is intended for annotation-bundle audits, explicit review-ledger
compilation, immutable patch creation, patch application, and byte-identical
reversion. It is not an automatic annotation service.

## Security and reproducibility properties

- Both build stages use Python 3.11.15 slim-bookworm at the pinned OCI index
  digest
  `sha256:d29f48a31a8b408ed19272ca1e7b10ebae13b240a27e862d3d4217c528e2e0c3`.
- The wheel is built with the image's setuptools 79.0.1 under
  `RUN --network=none`, `--no-build-isolation`, and `--no-deps`.
- Runtime installation also uses `RUN --network=none` and `--no-deps`.
- The final image runs as numeric user and group `10001:10001`.
- The build context is allowlisted by `.dockerignore`; genomes, results,
  environments, manuscripts, Git history, tests, and server paths are not sent
  to the builder.
- OCI metadata identifies version `1.0.0`, the source repository, and the
  `BSD-3-Clause` license.

## Build

From the repository root:

```bash
docker build \
  --pull \
  --file containers/Dockerfile \
  --tag ploidypatch-core:1.0.0 \
  .
```

Inspect the runtime identity and version:

```bash
docker image inspect ploidypatch-core:1.0.0 \
  --format '{{.Config.User}} {{json .Config.Entrypoint}}'

docker run --rm --read-only --network none \
  ploidypatch-core:1.0.0 --version
```

The expected identity is `10001:10001`, the entry point is `ploidypatch`, and
the version is `1.0.0`.

## Run the reversible example in the container

The following command mounts the example read-only, gives the non-root process
an ephemeral writable `/tmp`, disables runtime networking, and executes the
installed package rather than the source tree:

```bash
docker run --rm \
  --read-only \
  --network none \
  --tmpfs /tmp:rw,noexec,nosuid,size=16m \
  --mount type=bind,src="$PWD/examples/minimal_reviewed_patch",dst=/work/examples/minimal_reviewed_patch,readonly \
  --entrypoint python \
  ploidypatch-core:1.0.0 \
  /work/examples/minimal_reviewed_patch/run_example.py \
  --input-dir /work/examples/minimal_reviewed_patch \
  --output-dir /tmp/reviewed-patch
```

The final JSON report must contain two accepted additions,
`automatic_approval=false`, `byte_identical_reversion=true`, and identical
source/reverted SHA-256 values.

## Use host input and output directories

Keep inputs read-only. On Linux, either make the output directory writable by
UID 10001 or run the container with the current non-root host identity:

```bash
mkdir -p output
docker run --rm \
  --user "$(id -u):$(id -g)" \
  --network none \
  --mount type=bind,src="$PWD/input",dst=/input,readonly \
  --mount type=bind,src="$PWD/output",dst=/output \
  ploidypatch-core:1.0.0 \
  patch apply \
  --source-gff /input/source.gff3 \
  --patch /input/reviewed.patch.json \
  --output-gff /output/annotation.patched.gff3
```

Do not mount evaluator-only truth into a blind candidate container. Container
packaging does not replace the role and custody controls in formal prospective
protocols.

## Release boundary

The core image packages the portable review-and-patch CLI. Bioconda packaging
and workflow-specific environments remain separate because upstream
bioinformatics engines have their own dependency and license constraints.
