# Populus v0.4 execution patch 2

This is a second pre-reveal, execution-only compatibility patch. It changes no
species, reference role, CDS coordinate, candidate policy, model, threshold,
endpoint, gate or statistic.

Patch 1 successfully fixed the isolated `/bin/sh` failure and allowed LiftOn to
extract all Salix purpurea transcript/protein sequences. All four upstream
programs completed their computations. LiftOn then refused to publish the
S. purpurea result because its output validator detected 50 duplicate exon IDs.
No candidate pool, self-WGD candidate score, model score, custody artifact or
truth reveal existed.

The cause was diagnosed from the retained partial output and frozen LiftOn
source. LiftOn recognizes a trailing `-N` exon-ID family and may renumber that
family independently in each transcript after chaining. Patch 1's globally
unique `PPX-<hash>` source IDs therefore became repeated `PPX-1`, `PPX-2`, ...
families. Patch 2 uses the same deterministic coordinate hash without the
hyphen (`PPX<hash>`). It remains globally unique; if LiftOn rebuilds duplicated
exons, its own code now takes the transcript-ID-specific fallback. CDS rows,
coordinates, phases and hashes remain unchanged.

Output validation remains enabled. The failed patch-1 attempt is retained and
hashed by the patch-2 freeze; validation is not disabled or relaxed.
