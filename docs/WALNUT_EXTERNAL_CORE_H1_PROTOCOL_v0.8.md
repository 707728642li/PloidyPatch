# Walnut external core H1 protocol v0.8

Status: untouched pre-pair, pre-candidate and pre-label protocol skeleton,
2026-08-09. It is not an execution or result record.

## Confirmatory claim

The only scientific hypothesis is whether retaining distinct exact phased-CDS
chains recovers more deliberately hidden events than suppressing strongly
overlapping alternative chains. The estimand is the paired event-level
difference in exact phased-CDS-chain recall. Its uncertainty is a 20,000
replicate paired event bootstrap; success requires both a positive point
estimate and a two-sided 95% interval lower bound above zero.

There is no H2, candidate ranking, topology ranking, calibration, automatic
approval, or model threshold. Candidate-independent stable-reference ranker
v0.9 failed its development retention gates and is retired for Walnut. Its
only permitted consequence is to retain the chain-preserving H1 workflow.

## Biological event audit

Public comparative work places the juglandoid WGD in the ancestor of
Juglandaceae and finds the two surviving paleo-subgenomes in both *Juglans*
and *Carya*. Corylus (Betulaceae) and Castanea (Fagaceae) diverged outside this
event. Consequently, a target self-WGDI pair supported as the same exact 1:2
relationship by both evaluator families is a biologically defensible,
conservative subset of retained juglandoid homeologs.

The WGD is ancient and is commonly interpreted as allopolyploid. Biased
fractionation, paralog loss and deeper core-eudicot gamma homology make naive
multiplicity unsafe. Formal pairs therefore require all of the following:

1. a target self-WGDI block with at least 20 pairs, different primary
   chromosomes and pair-consistent reciprocal uniqueness;
2. closed `0.10 <= YN00 Ks <= 0.75`, fixed from published Juglandaceae ranges
   before any Walnut pair enumeration and separated from the reported older
   gamma peak;
3. one or more exact evaluator counterparts in each of Corylus and Castanea,
   with every accepted counterpart resolving to exactly the same two target
   genes; and
4. rejection, without truncation or best-hit rescue, of missing Ks,
   out-of-range Ks, more-than-two multiplicity, identifier ambiguity or
   discordant counterparts.

This rule does not assert that every true retained homeolog remains observable.
It prioritizes specificity and reports any event shortfall without relaxing
the rule. Supporting literature includes the
[walnut genome comparison](https://pmc.ncbi.nlm.nih.gov/articles/PMC4574618/),
[Juglandaceae synteny/subgenome analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC9899254/),
and [Juglandaceae genome evolution study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6497033/).

## Fixed data and role firewall

The five references and the historical `>=0.95` primary-chromosome base gate
are exactly those in the committed metadata record. Roles are one target, two
candidate-only references (*J. mandshurica*, *C. illinoinensis*) and two
evaluator-only references (*C. avellana*, *C. mollissima*). Candidate
execution receives only the target genome, candidate references and a sealed
truth-free blind benchmark. It cannot mount evaluator references, the complete
target annotation, truth, labels, NAS, or a network.

The requested sample is 800 nonoverlapping events; at least 500 events, 12 of
16 primary chromosomes (`ceil(0.75 * 16)`) and 20 events in each of the four
CDS-segment bins are required. Any shortfall is not evaluable without rule
relaxation.

## Controls and safety gates

Every arm must include no-op and oracle sentinels, byte-identical restoration,
identical complete/blind genome SHA-256, and zero loss of baseline transcript
structure. A failure invalidates the run rather than changing thresholds. The
complete control is evaluator-generated only after blind raw predictions and
pool custody hashes are frozen.

## Current implementation boundary

The contract, policy, event definition, metadata preflight, no-ranker sentinel,
generic role stager and protocol-freeze adapter are implementation skeletons.
They may validate and freeze bytes, but this commit does not enumerate pairs,
generate candidates, build truth, access labels, or execute the formal test.
