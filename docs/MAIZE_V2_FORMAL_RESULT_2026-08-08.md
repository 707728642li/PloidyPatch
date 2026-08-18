# Maize v0.2 formal zero-retuning result (2026-08-08)

## Verdict and claim boundary

The predeclared maize external gate **failed**. This result is retained as the
formal v0.2 result and will not be replaced by maize-specific tuning. The
negative result shows that the homeolog-topology features contain biological
signal, but the jointly refitted v0.2 topology estimator does not provide a
portable ranking improvement in this distant monocot.

All values below were produced after the candidate, model, protocol, and
evaluation artifacts had been checksum-frozen. The model SHA-256 is
`e8d0227628bc0688504f1132cca3cd1e418b87a6dbdff39705d35c9c8f933958`.

## Formal method comparison

The endpoint is an exact phased CDS-chain match among 800 reversible maize
copy-collapse events. “Union” is the existing consensus selector with minimum
method support one; it is not a literal union because overlapping alternative
chains are suppressed.

| Method/policy | TP | FP | FN | Precision | Recall | F1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| miniprot | 459 | 279 | 341 | 0.6220 | 0.5738 | 0.5969 |
| GeMoMa | 609 | 171 | 191 | 0.7808 | 0.7612 | **0.7709** |
| LiftOn | 391 | 263 | 409 | 0.5979 | 0.4888 | 0.5378 |
| support >= 1 selector | 566 | 220 | 234 | 0.7201 | 0.7075 | 0.7137 |
| exact 2-of-3 consensus | 504 | 102 | 296 | 0.8317 | 0.6300 | 0.7169 |
| exact 3-of-3 consensus | 302 | 23 | 498 | 0.9292 | 0.3775 | 0.5369 |

GeMoMa exceeded the support-at-least-one selector by 0.0571 F1 with a paired
95% bootstrap interval of 0.0231 to 0.0917. The support-at-least-one and 2-of-3
F1 values were statistically indistinguishable; 2-of-3 offered the expected
higher-precision review tier.

## Frozen topology gate

There were 15,636 ranked candidates and 566 exact positives, so the frozen
candidate policy imposed a recall ceiling of 0.7075. Topology was available
for 423/566 positives (0.7473), passing the predeclared 0.70 coverage floor.

| Estimator | Average precision | ROC AUC |
| --- | ---: | ---: |
| paired baseline | **0.65285** | **0.96961** |
| homeolog topology | 0.65112 | 0.96768 |

Topology-minus-baseline AP was -0.00173. The predeclared 20,000-replicate
chromosome-grouped bootstrap 95% interval was -0.00887 to 0.00452. Both the
positive-delta requirement and the lower-bound-above-zero requirement failed.
At review budgets 79, 157, and 313 (0.5%, 1%, and 2%), topology also recovered
2, 4, and 6 fewer exact positives than baseline, respectively.

## Failure analysis (post hoc; diagnostic only)

No model was fitted or selected on these diagnostics. They do not alter the
formal verdict.

- Topology availability was itself strongly label-associated: 50.54% positive
  when available versus 0.97% when unavailable (risk ratio 52.3). Within the
  topology-available stratum, however, AP declined from 0.7581 to 0.7489.
- Individual coding-topology measures were directionally informative, but the
  joint refit nearly cancelled their contribution. For available positives,
  the topology add-on changed the linear score by +6.225 on average while
  refitting the original terms changed it by -6.133.
- The add-on helped sparse method-support strata, especially GeMoMa-only
  candidates (AP +0.155), but slightly harmed several robust multi-method
  strata. A global reorder is therefore too blunt; any next model should treat
  method-support regime explicitly.
- Of 800 hidden events, at least one individual method recovered 653 (oracle
  candidate recall 0.8163), whereas the current support-at-least-one selector
  retained only 566. It lost 87 exact structures present in an individual
  method, including 66 GeMoMa structures, because its >0.5 CDS-overlap rule
  suppresses distinct alternative chains. This loss occurs before ranking.

## Consequences for v0.3

The highest-priority correction is a chain-preserving candidate pool that
deduplicates identical phased CDS chains but retains overlapping alternative
chains as explicit conflicts for ranking and review. This is an information-
preserving change, not a claim that every retained model should be patched.

After that candidate contract is frozen using soybean and Brassica development
data, topology should be evaluated as a support-conditioned correction rather
than a global refit. Maize can be used only as a transparent post-hoc audit for
v0.3; a newly registered species is required for the next confirmatory external
test.

## Natural-genome firewall status

The unperturbed B73 NAM-5.0 discovery arm completed its candidate and top-K
rank freeze on 2026-08-08 before Iso-Seq sequence access. It contains 14,853
candidates ranked by each paired estimator, with fixed review budgets 25, 50,
100, and 200 and no automatic approval. Subsequent transcript and pan-genome
intersections are validation, not training.
