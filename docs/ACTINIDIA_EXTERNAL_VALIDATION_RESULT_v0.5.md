# Actinidia untouched secondary external replication result v0.5

Date: 2026-08-08

## Decision

The preregistered *Actinidia chinensis* Red5 evaluation is complete and
formally evaluable. Its fixed outcome is **formal negative external result**,
not a confirmatory pass. The only failed preregistered gate was
`topology_positive_coverage`. H1 and H2 were both positive in their fixed
sequence, and every other formal and safety gate passed.

This outcome is retained without changing the model, candidate policy,
reference roles, topology threshold, event definition, endpoint or formal
decision rule. The positive component results do not override the failed
composite gate.

## Immutable evidence binding

The manuscript source summary is a lossless selection of fields from the final
evaluator output and is checked against that source in the manuscript evidence
tests.

- final `evaluation.json` SHA256:
  `f3c187bd5e86cf172343b7ac2b4ad10f331b0d4ab8a2273c5540c333109a0c3e`;
- final evaluator output size: 23,412 bytes;
- protocol `SHA256SUMS` SHA256:
  `4cab1d6aa2eefcff4ce2a3d1d3fdecf8fd8d9dfc2fd6d506dd18cf1ec1a3bb8c`;
- composite model `SHA256SUMS` SHA256:
  `90a5d90b6056f883df32af5556d906849054b0cc1880537cd4c35bd07a3a5797`;
- blind custody manifest SHA256:
  `7ac3e1a30b325a1010264c96fa592f5b983340f8903474a95f3e4cdad7eb62c5`;
- tracked local evidence copy:
  `manuscript/source_data/actinidia_external_v0.5_evaluation.json`;
- frozen server path:
  `/data/codexli/projects/PloidyPatch/results/evaluator/actinidia/v0.5/reveal_patch1_run2/evaluation/evaluation.json`;
- manuscript/figure source:
  `manuscript/source_data/actinidia_external_v0.5_summary.json`.

## Evaluability

The evaluator obtained 800 events across all 29 target linkage groups. The
four CDS-complexity bins contained 192, 207, 206 and 195 events. No fixed rule
was relaxed. Formal evaluability, all four data gates and every sentinel
passed, including zero no-op recovery, 800/800 oracle recovery, byte-identical
restoration, identical blind/complete genome hashes and the pair-truth audit.

## Fixed-sequence results

| Test | Comparator | PloidyPatch v0.4 | Difference | 95% bootstrap interval | Gate |
| --- | ---: | ---: | ---: | ---: | --- |
| H1 exact phased-CDS event recall | legacy suppression: 333/800 (0.41625) | chain-preserving pool: 406/800 (0.50750) | +0.09125 | +0.07125 to +0.11125 | pass |
| H2 candidate average precision | frozen non-topology baseline: 0.383198 | v0.4 guard: 0.426015 | +0.042817 | +0.029529 to +0.056886 | pass |

H1 used 20,000 paired event-bootstrap replicates. H2 was tested after H1
passed and used 20,000/20,000 valid target-linkage-group bootstrap replicates.

## Sole failed gate and passed safety gates

Effective topology evidence was available for 252 of 406 exact-positive
candidates, a coverage fraction of 0.6206896552. This failed the preregistered
`topology_positive_coverage` gate and therefore fixed the formal-negative
outcome.

A checksum-bound post-reveal audit reproduced this numerator exactly. None of
the 406 positives was removed by the v0.4 conflict guard: 252 had usable
topology and 154 lacked it before guarding. The method-support decomposition
localized the limitation to sparse candidate support. Topology was available
for 186/210 positives supported by all three upstream methods (88.57%), but for
only 39/100 positives supported by exactly two methods and 27/96 supported by
one method. This audit is descriptive and does not change the frozen threshold,
scores or formal outcome.

Every other gate in `evaluation.json` was `true`:

- formal evaluability;
- H1 chain-preserving candidate ceiling;
- H2 numerical improvement and fixed-sequence testing;
- all requested bootstrap replicates valid;
- conflict-winner mapping identical to the baseline (`mismatch_count = 0`);
- top-1% review noninferiority;
- retention of at least 90% of the positive v0.3 AP gain;
- zero collateral loss in all seven arms;
- automatic approval absent (`automatic_approved = 0`).

The v0.4 guard retained 0.9985560569 of the positive v0.3 AP gain. The
evaluation contained 70,097 candidate chains and 406 exact-positive candidates.
These are review-ranking results on the controlled benchmark, not calibrated
probabilities and not authorization for automatic annotation changes.

## Interpretation

Actinidia independently supports both component findings: preserving distinct
phased CDS chains increased recoverable exact structures, and the frozen
ranker improved review prioritization over its non-topology baseline. It also
identifies a portability limit: usable topology evidence did not cover enough
positive candidates to satisfy the full preregistered policy. The publication
record must therefore report the positive H1 and H2 effects together with the
formal-negative composite outcome.
