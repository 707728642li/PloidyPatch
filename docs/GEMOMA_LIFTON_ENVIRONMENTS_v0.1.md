# GeMoMa 1.9 and LiftOn 1.0.11 environment freeze

Date frozen: 2026-08-07  
Server project: `/data/codexli/projects/PloidyPatch`

## GeMoMa 1.9

GeMoMa is isolated at `envs/ploidypatch-gemoma`. It was installed from the
configured Bioconda channel as `gemoma-1.9-hdfd78af_0` under GPL-3.0. The
package archive SHA-256 is
`1c19229301adc2a3afc28e25b3d62e242ea92adf67f62fe7e0875316f39d581c`.

Frozen executable identities:

- launcher SHA-256:
  `e33205a7851f205f4a7017b27d57e1b98a833cfa1eae06770e619447fc5cab7b`;
- `GeMoMa-1.9.jar` SHA-256:
  `dd30c1a615539a45903aa5b162bdba23839fd45ebc8afbc5b61938b8bc1886a3`;
- BLAST 2.17.0+, MMseqs2 18.8cc5c, OpenJDK 11.0.30, and Python 3.14.5.

Creation took 22.95 seconds and peaked at 635,836 kB RSS. The immutable conda
specification, package list and metadata, tool help, dependency versions,
resource usage, per-file hashes, and audit manifest are stored at
`results/baselines/gemoma_v1.9/environment`. The audit is reproducible with
`scripts/audit_gemoma_environment.sh` and refuses to overwrite prior evidence.

## LiftOn 1.0.11

LiftOn is isolated at `envs/ploidypatch-lifton`. The official PyPI source
archive is retained at `third_party/downloads/lifton-1.0.11.tar.gz`, with
SHA-256
`c2125db9bede3640e13c6aa6a0672887aaf2611710442f9bf69017197d188f88`.
The installed release reports 1.0.11 from both its CLI and Python import and
passes `pip check`.

Frozen executable and aligner identities:

- LiftOn launcher SHA-256:
  `460f640530ca5a03f0cf4437df8f127c14febaaa3f9766da88fb4512893c65f8`;
- installed `lifton.py` SHA-256:
  `9ccbc82dc716bacdcb6a8d1bc113897fcbbc9964cfd84f0ba4a921ab8921d7e4`;
- minimap2 `2.30-r1287`, package archive SHA-256
  `a82a861a2c0e5a01bbe73701dcc6af8ba3f1da430a9fa6c94ed53b5c4a7a2b55`;
- miniprot `0.18-r281`, package archive SHA-256
  `a7785a0db01773819ba1c3b897d2987a5ec18f4c80165096c87d924c7e4e845b`;
- Python 3.11.15.

The conda prefix was created in 8.69 seconds with 230,436 kB peak RSS. The
successful source installation took 64.38 seconds with 221,908 kB peak RSS.

### Preserved installation failure

The first installation attempt from the same verified archive failed
transactionally while building `mappy`: compilation of `bseq.c` reported
`fatal error: zlib.h: No such file or directory`. Its stdout, stderr, and time
record remain in `logs/environments/lifton_1.0.11/pip_install.*`. We added
`zlib=1.3.2` to this project-only prefix, verified the header, and repeated the
same pip command against the same archive. The successful retry is retained as
`pip_install.retry1.*`; no alternate source or LiftOn revision was used.

The immutable conda and pip specifications, `pip inspect`, distribution
metadata, help flags for the 1.0.11 behavior switches, tool versions, hashes,
and resource logs are stored at
`results/baselines/lifton_v1.0.11/environment`. The audit is reproducible with
`scripts/audit_lifton_environment.sh` and refuses to overwrite prior evidence.

## Readiness decision

Both comparator environments pass the frozen executable and dependency
contracts. This only establishes reproducibility; it is not a biological
validation or a performance result. Raw prediction launchers must still freeze
their inputs, full commands, completion sentinels, output hashes, wall/CPU/RSS,
and disk use before the matched Glycine score is released.
