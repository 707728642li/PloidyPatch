# Bioconda submission

This directory contains the submission-ready Bioconda recipe. The local recipe
at `packaging/conda/meta.yaml` uses the current checkout for pre-release
validation; `packaging/bioconda/meta.yaml` is prepared for the immutable v1.0
GitHub release source distribution and its canonical SHA-256 digest. The
recipe is deliberately excluded from that source distribution, so updating the
recipe checksum cannot change the archive it verifies.

Submission procedure:

1. copy `meta.yaml` to `recipes/ploidypatch/meta.yaml` in a fork of
   `bioconda/bioconda-recipes`;
2. build and test with `bioconda-utils`, including the mulled container test;
3. submit the recipe to `bioconda/bioconda-recipes`.

The repository CI and `tests/test_v1_release_package.py` reject an official
recipe whose package version, license, entry point, or maintainer differs from
the release metadata.

This directory is submission-ready evidence, not proof that the recipe is
already merged into Bioconda. Availability may only be claimed after the
external Bioconda pull request and channel package are independently verified.
