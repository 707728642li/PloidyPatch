# PloidyPatch candidate-pool ranker v0.3 result (2026-08-08)

## Verdict

The predeclared support-conditioned offset ranker passed all development
retention criteria and is retained for a new external freeze. This is a model-
development result, not an external confirmation and not an automatic gene-
model approval rule.

The primary estimator leaves the fold-specific baseline logit unchanged when
topology is unavailable. Where topology is available, it learns a correction
from five global topology components plus support-pattern-specific availability
and mean-coherence terms. The estimator, `C = 1`, review budgets, species roles,
and retention criteria were fixed before this evaluation.

## Audited repeated chromosome OOF result

Five randomized five-fold partitions were used per development species. Each
partition held out complete chromosomes, train and test chromosomes were
disjoint, every fold contained both classes, and all five partitions were
required to be distinct.

| Development species | Candidates | Exact positives | Baseline AP | Primary AP | Delta AP | 95% chromosome-bootstrap CI | Conflict top-1 delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Soybean | 20,504 | 630 | 0.8534 | 0.8766 | +0.0231 | 0.0106 to 0.0389 | +0.0465 |
| Brassica | 76,434 | 653 | 0.6136 | 0.6605 | +0.0469 | 0.0324 to 0.0597 | +0.0458 |
| Pooled, within-species percentile scores | 96,938 | 1,283 | 0.6017 | 0.6348 | +0.0331 | 0.0205 to 0.0448 | not pooled |

Every randomized partition gave a positive AP difference:

- soybean: +0.0207 to +0.0254;
- Brassica: +0.0454 to +0.0486.

At fixed review budgets, the primary estimator increased exact true positives
from 92 to 95 in the top 0.5% and from 187 to 193 in the top 1% for soybean. In
Brassica it increased them from 278 to 296 and from 461 to 483, respectively.
Conflict-set top-1 accuracy rose from 0.7674 to 0.8140 in soybean and from
0.7647 to 0.8105 in Brassica.

## Transfer and retrospective diagnostics

| Training -> evaluation | Evaluation role | Baseline AP | Primary AP | Delta AP |
| --- | --- | ---: | ---: | ---: |
| Brassica -> soybean | development transfer | 0.7966 | 0.8264 | +0.0298 |
| Soybean -> Brassica | development transfer | 0.5763 | 0.6425 | +0.0661 |
| Pooled development -> cotton | retrospective diagnostic | 0.7099 | 0.7475 | +0.0376 |
| Pooled development -> maize | retrospective diagnostic | 0.5950 | 0.5977 | +0.0027 |

Cotton conflict top-1 accuracy rose from 0.6540 to 0.6967. Maize conflict
top-1 was nearly unchanged (0.6374 to 0.6410), while primary top-budget counts
were slightly lower than baseline (for example, 78 versus 79 exact chains in
the top 0.5%). Thus v0.3 removes the clear development-species failure of the
old topology design and transfers strongly to cotton, but it does not turn
maize into positive confirmatory evidence.

Neither cotton nor maize was allowed to select the estimator. They remain
retrospective diagnostics because both had already been inspected during prior
method development.

## Splitter audit

The first execution used `StratifiedGroupKFold(shuffle=True)` and passed the
headline numerical criteria, but all five seeds produced the same chromosome
partition. That artifact is not used as the formal repeated-CV result. It is
preserved on the server under
`candidate_pool_ranker_v0.3_evaluation.invalid_repeated_partition_88da631`.

The evaluator was corrected without changing the model or criteria to use
randomized `GroupKFold`, with explicit distinct-partition, chromosome-leakage,
and label-class audits. The canonical rerun uses schema
`ploidypatch.candidate_pool_ranker_evaluation.v4` and code commit `d1b42c8`.

## Claim boundary and next gate

The result is sufficient to freeze v0.3 for a new preregistered external
species. It is not sufficient to claim broad portability, because soybean and
Brassica are development species, cotton and maize are no longer untouched,
and the maize review-budget result remains slightly unfavorable.

The next confirmatory gate must therefore use a species whose labels and
method outputs have not influenced candidate, feature, or model selection.
Natural maize full-chain cases also still require ORF, repeat, collision, and
raw-read audits before they can support annotation-correction claims.

## Canonical artifacts

- Server result:
  `/data/codexli/projects/PloidyPatch/results/copy_collapse/model_development/candidate_pool_ranker_v0.3_evaluation`
- Local compact mirror: `results/candidate_pool_ranker_v0.3_evaluation`.
- Code commit: `d1b42c8`
- All server files listed in `SHA256SUMS` pass `sha256sum -c`.
