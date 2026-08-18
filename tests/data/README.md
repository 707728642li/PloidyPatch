# Test data

The public CI fixture is the checksum-bound dataset in
`examples/minimal_reviewed_patch/`. It is intentionally tiny, synthetic, and
redistributable. The fixture exercises candidate conflicts, explicit review,
patch creation/application, and byte-identical reversion without downloading
any reference genome or contacting the network.

Unit tests create additional synthetic GFF3/FASTA/TSV inputs inside temporary
directories. Large third-party genomes, annotations, reads, and alignments are
never required by GitHub Actions and are not stored in this directory.
