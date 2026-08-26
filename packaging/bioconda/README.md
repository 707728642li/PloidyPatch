# Bioconda submission

PloidyPatch `1.0.1` is published and installable from the Bioconda channel.
This directory contains the published metadata-correction recipe, bound to the
immutable `v1.0.1` source release. The local recipe
at `packaging/conda/meta.yaml` uses the current checkout for pre-release
validation; `packaging/bioconda/meta.yaml` is prepared for the immutable v1.0
GitHub release source distribution and its canonical SHA-256 digest. The
recipe is deliberately excluded from that source distribution, so updating the
recipe checksum cannot change the archive it verifies.

Update procedure for a later patch release:

1. copy `meta.yaml` to `recipes/ploidypatch/meta.yaml` in a fork of
   `bioconda/bioconda-recipes`;
2. build and test with `bioconda-utils`, including the mulled container test;
3. submit the recipe to `bioconda/bioconda-recipes`.

The repository CI and `tests/test_v1_release_package.py` verify the published
recipe's version, license, entry point, maintainer and immutable source hash.
Availability claims are made only after the external channel exposes the
package; that verification has now been completed for `1.0.1`. The original
`1.0.0` channel package remains part of the historical release record.
