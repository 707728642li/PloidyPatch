# Walnut Carya transcript-hierarchy patch v0.8

Walnut blind run 1 stopped before any candidate method executed. The sealed
candidate-input adapter rejected 35,226 exon rows whose immediate parent was
not typed `mRNA` or `transcript`. Read-only inspection showed that all 15,376
parent features are unique existing provider records, each has exactly one
resolved gene parent, and every feature, exon and gene agrees in sequence ID,
strand and containment. The provider feature types are noncoding RNA classes
or pseudogene records; their exon coordinates and parentage are not ambiguous.

This execution-only compatibility patch preserves an exon-only container
without treating it as an `mRNA` only when its ID occurs exactly once and all
of its exons agree with its sequence ID and strand and lie inside its existing
bounds. A non-`mRNA` container with CDS/UTR children can enter the internal
transcript-like set only when it also has exactly one resolved gene parent.
Any missing, duplicate, multi-parent, cross-sequence, cross-strand or
out-of-bounds case continues to fail closed. A unique parentless provider
feature that is the sole parent of an existing `mRNA`/`transcript` is treated
as the gene-like root only after the same sequence, strand and containment
checks; this covers 14 transcribed pseudogene roots in Carya. The adapter does
not add a
transcript, alter a child coordinate, phase, ID or Parent, or change any CDS
byte; valid provider-specific exon-only subtrees pass through byte-for-byte.

No candidate method, reference role, event, threshold, truth pair, sampling
rule, hypothesis or success gate changes. The failed run remains retained and
no candidate projection, pool, score or truth label existed when this patch
was defined.
