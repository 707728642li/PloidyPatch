# SynGAP 1.2.5 environment audit

Date frozen: 2026-08-07  
Server prefix: `/data/codexli/projects/PloidyPatch/envs/ploidypatch-syngap`  
Status: primary `genblastg` environment validated

## Package identity

The environment was created through the server's configured conda-forge and
Bioconda channels, without a mirror or alternate source. The installed package
is `syngap 1.2.5 py312hdfd78af_0` from Bioconda. Its package archive is
515,837,080 bytes with SHA-256
`dc107d59bb6ecd3b76fe781b539fd60ba2521f303dae1d0ab4349a212f58408b`.
The installed entrypoint SHA-256 is
`7e43b773357220e75b53ecd3ed6dd641d35832e06d7ec7b0053035138390e428`.

The frozen tag and Bioconda metadata declare `CC-BY-NC-SA-4.0`. This overrides
the misleading GPL-3.0 classification shown for the repository's current
default branch. SynGAP remains a separately invoked, non-commercial external
comparator; its code is not incorporated into PloidyPatch.

The complete explicit specification, JSON/text package lists, conda metadata,
help text, hashes, and resource paths are under
`results/baselines/syngap_v1.2.5/environment` on the server. Two failed audit
probes are retained under the sibling `_incomplete` directory: one established
that EMBOSS wrappers require the activated PATH, and the other corrected the
actual bundled genBlast executable name. The final audit executes all probes
through `conda run -p` and passes.

## Selected dependency versions

| Component | Frozen version |
|---|---|
| Python | 3.11.15 |
| JCVI | 1.6.6 |
| bedtools | 2.31.1 |
| gffread | 0.12.9 |
| seqkit | 2.13.0 |
| DIAMOND | 2.2.2 |
| LAST | 1651 |
| EMBOSS | 6.6.0.0 |
| BioPerl | 1.7.8 |

Mamba resolved 496 packages with 684 MB of downloads. Environment creation
took 21:44.42 and peaked at 14,922,068 kB RSS; the installed prefix occupies
approximately 6.2 GB. The transaction reported stale legacy-package metadata
warnings and repaired an invalid cached R file, but exited zero and the
post-install executable/resource audit passed.

## Bundled prediction resources

The primary author's default backend and downstream filters depend on bundled
assets. Their installed hashes are:

| Resource | SHA-256 |
|---|---|
| CPC2 `CPC2.py` | `bc6bbfc37637ad907f5dade8b9e4d9139cee73edc87dad35385de3229218316d` |
| genBlast v1.38 executable | `50dcb1bedf3d810c8aef38f6f172b17506fe46c04f8fc3b5ec505dcdb68157a7` |
| bundled SwissProt archive | `17284b2caeaf5734e8878ce03775c874158491a88db105fc9353904bf2e377f1` |

The package does not declare or install miniprot. Therefore the frozen primary
`genblastg` comparison uses this environment unchanged; the secondary
documented miniprot mode, if run, must use a separately frozen environment and
cannot silently modify this prefix.

## Execution risk controls

Source inspection at tag `v1.2.5` shows that SynGAP builds many unquoted
`os.system` commands and often ignores their return values. The versioned
`scripts/run_syngap_dual.sh` wrapper consequently rejects unsafe paths and
requires output directories, clean/full annotations, anchors, completion
messages, absence of fatal log signatures, and exact target-GFF prefix
preservation. A top-level zero exit code alone is not accepted as a successful
baseline run.
