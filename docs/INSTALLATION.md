# Installing PloidyPatch v1.0

PloidyPatch separates its dependency-light review/patch core from optional
scientific workflow programs. This keeps the command-line package portable and
prevents an installation from silently selecting third-party tools or versions.

## Supported platforms

- **Linux:** the formally supported platform for production annotation,
  projection, WGD, long-read and validation workflows.
- **macOS:** the portable review/patch core is covered by the Bioconda build;
  individual upstream bioinformatics programs may not be available.
- **Windows:** native Windows is not supported and there is no Windows-specific
  package. Use WSL2, a Linux server, or the container for real analyses.

The `py3-none-any` wheel describes a pure-Python core; it does not imply that
the complete external bioinformatics stack is supported on every operating
system.

## Recommended: Pip / PyPI

```bash
python -m pip install ploidypatch==1.0.0
ploidypatch --version
```

## Recommended: Conda

The tested noarch Conda package can be installed directly from the immutable
v1.0 release (registry version `1.0.0`):

```bash
conda create -n ploidypatch python=3.11
conda install -n ploidypatch \
  https://github.com/707728642li/PloidyPatch/releases/download/v1.0/ploidypatch-1.0.0-py_0.conda
conda run -n ploidypatch ploidypatch --version
```

After the Bioconda recipe is merged and propagated, the channel form is:

```bash
conda create -n ploidypatch -c conda-forge -c bioconda ploidypatch=1.0.0
conda activate ploidypatch
ploidypatch --version
ploidypatch --help
```

The tagged source can also be installed with pip inside a Conda environment:

```bash
conda create -n ploidypatch -c conda-forge python=3.11 pip
conda activate ploidypatch
python -m pip install \
  "ploidypatch @ git+https://github.com/707728642li/PloidyPatch.git@v1.0"
```

## Optional Python features

The reviewed-patch path has no mandatory third-party Python runtime
dependencies. Install optional groups only when needed:

```bash
# NumPy and scikit-learn for research evaluation commands
python -m pip install ".[evaluation]"

# Matplotlib for local figure generation
python -m pip install ".[visualization]"

# Complete developer and CI environment
python -m pip install -e ".[test]"
```

## Container

Build the digest-pinned, non-root core image:

```bash
docker build -f containers/Dockerfile -t ploidypatch:1.0.0 .
docker run --rm --read-only --network none ploidypatch:1.0.0 --version
```

The container does not bundle miniprot, GeMoMa, LiftOn, DIAMOND, WGDI,
SynGAP, aligners, or read preprocessors. It is deliberately scoped to the
portable core and the reviewed-patch workflow.

## External workflow programs

Research commands may consume outputs from the following separately managed
tools:

| Program | Role | Bundled? |
|---|---|---:|
| Miniprot | protein-to-genome projection | No |
| GeMoMa | homology-supported gene prediction | No |
| LiftOn | annotation lift-over | No |
| DIAMOND | protein similarity for synteny preparation | No |
| WGDI | collinearity and WGD analysis | No |
| minimap2 / long-read tools | transcript validation | No |

Formal manuscript workflows record exact versions and explicit environment
files. Do not substitute versions inside a frozen validation run.

## Verify an installation

```bash
ploidypatch --version
ploidypatch --help
python -c "import ploidypatch; print(ploidypatch.__version__)"
python examples/minimal_reviewed_patch/run_example.py \
  --output-dir work/install-smoke
```

Expected version: `1.0.0`. The example must report two accepted additions,
`automatic_approval=false`, and `byte_identical_reversion=true`.

## Troubleshooting

### `ModuleNotFoundError: numpy`

The requested research evaluation command uses the optional evaluation stack.
Install `.[evaluation]`; the core review and patch commands do not require it.

### Existing output path

PloidyPatch intentionally refuses to overwrite outputs. Choose a new path or
archive the previous result after checking its manifest.

### Upstream program missing

Install it in a workflow-specific Conda environment and record the explicit
environment. PloidyPatch will not install or invoke an undeclared substitute.
