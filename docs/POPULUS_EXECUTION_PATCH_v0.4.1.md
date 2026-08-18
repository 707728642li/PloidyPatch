# Populus v0.4 execution patch 1

This is a pre-reveal, execution-only patch. It does not change the frozen
species, reference roles, target event, truth construction, candidate policy,
model, score, thresholds, endpoints, gates or statistics.

The first isolated blind attempt exited before candidate-pool construction,
self-WGD scoring, model scoring, custody sealing or truth reveal. Its two
deterministic execution failures were diagnosed only from candidate-safe logs
and provider annotations:

1. Ubuntu exposes `/bin` as a symlink to `/usr/bin`. The original bubblewrap
   launcher mounted `/usr` but failed to recreate `/bin`, while MMseqs generated
   scripts use `#!/bin/sh`. Patch 1 recreates the immutable `/bin -> usr/bin`
   alias inside the namespace.
2. The Salix purpurea provider GFF3 contains transcript CDS/UTR children but no
   explicit exon rows. LiftOn therefore generated an empty transcript FASTA.
   Patch 1 appends deterministic exon compatibility rows only for transcripts
   with zero exons, using the exact union of their CDS/UTR intervals. It asserts
   that all original CDS rows and hashes are unchanged and refuses unresolved
   transcripts.

The failed attempt is retained and hashed by the patch freeze. No method is
dropped, no input role is changed and no scientific rule is relaxed.
