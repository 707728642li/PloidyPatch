# PloidyPatch scaling benchmark protocol v0.1

Status: design freeze before the release scaling run.

## Purpose

Measure the computational behavior of the PloidyPatch-specific core separately
from mature upstream aligners and predictors.  The benchmark must support a
practical genome-scale software claim without presenting upstream miniprot,
GeMoMa, LiftOn, DIAMOND, or WGDI cost as a new PloidyPatch algorithmic result.

## Frozen workloads

Use complete, checksum-frozen real candidate universes rather than duplicated
synthetic rows.  The primary workloads are the current chain-preserving pools
for soybean, Brassica, cotton, maize, apple, and untouched Populus.  Smaller
TAIR10 and rice workloads may be included as low-end anchors.  For every
workload record:

- species, assembly and annotation release;
- genome size and primary chromosome count;
- base gene/transcript counts;
- method-family input and merged-candidate counts;
- conflict-set count and size distribution; and
- all input and policy SHA-256 hashes.

No species or replicate is removed because its timing is unfavorable.  Failed
runs remain in the raw table with the failure reason.

## Timed components

1. Chain-preserving exact method-family consensus and conflict construction.
2. Truth-blind candidate feature construction, including WGD/topology joins.
3. Frozen v0.3 portable rank scoring.
4. v0.4 conflict-winner guard and manifest verification.
5. End-to-end PloidyPatch core from already generated method candidates to the
   review-ranked, guarded table.

Upstream projection and synteny stages are reported in a separate descriptive
table using their frozen execution logs and exact thread counts.  They are not
combined with core timing when estimating algorithmic scaling.

## Execution contract

- Run on `bioserver`; record CPU model, physical/logical cores, RAM, kernel,
  filesystem, Python, environment lock, and full code commit.
- Use the production CLI and frozen model/composite artifact; no profiling-only
  alternate implementation.
- Preserve the production thread setting for each component.  Also run a
  thread-scaling arm at 1, 4, 16, 32, 64, and 128 threads only for components
  that actually expose parallelism.
- Do not drop operating-system caches or require privileged commands.  Label
  the first execution `cold_order_first`; run three subsequent repetitions and
  report their median and range as `warm_repeats`.
- Randomize species order for warm repetitions with literal seed `20261006`.
- Use `/usr/bin/time -v` plus a wrapper that records monotonic wall time, exit
  status, peak RSS, user/system CPU, filesystem inputs/outputs, and output
  bytes.  Capture stdout/stderr and refuse overwrite.
- Verify output SHA-256 and semantic equality against the production frozen
  result after every repetition.  A faster but different output is a failed
  replicate.
- Write into `.working` and atomically freeze a read-only result tree with
  `SHA256SUMS`.

## Reported endpoints

For every component/workload/repetition report wall seconds, CPU seconds, peak
RSS, output bytes, candidates per second, and conflict sets per second where
applicable.  Summaries report the first-run value and median/range of three warm
repetitions; no fastest-run-only reporting.

Fit descriptive log-log regressions of wall time and peak RSS against candidate
count, with all points shown.  Report slope, confidence interval, residuals,
and leave-one-species-out sensitivity.  These regressions describe empirical
scaling over the tested range and are not asymptotic complexity proofs.

The release claim is limited to the largest successfully verified real
workload.  A component is considered operationally acceptable only if it
completes that workload without swap, out-of-memory termination, output
divergence, or manual intervention.  There is no post hoc removal or input
simplification to make the gate pass.

## Figure and source-data outputs

- wall time versus candidate count, with cold and warm runs distinguished;
- peak RSS versus candidate count;
- component share of end-to-end core time;
- thread scaling and parallel efficiency for eligible components; and
- a machine-readable TSV containing every replicate and failure.

The plotting script must consume only the frozen timing TSV and environment
manifest, emit figure source data, and reproduce every value in the manuscript
resource statement.

## Frozen implementation

The tracked implementation is deliberately separate from the public CLI:

- `config/core_scaling_workloads_v0.1.json` fixes the six real workloads and
  role-separated inputs;
- `scripts/freeze_core_scaling_references_v0.1.py` reconstructs every workload
  once with the production CLI, proves byte identity to the existing candidate
  pools and feature matrices, and creates the v0.3/v0.4 reference outputs;
- `scripts/run_core_scaling_benchmark_v0.1.py` runs one first-order execution
  plus exactly three seeded warm repetitions for all five components; and
- `scripts/core_scaling_common_v0_1.py` contains the shared command builder,
  `/usr/bin/time -v` parser, exact-output verifier, and read-only freezer;
- `scripts/summarize_core_scaling_benchmark_v0.1.py` reports the first run and
  median/full range of all three warm repeats, component shares, log-log fits,
  residuals, and leave-one-species-out slopes; and
- `scripts/plot_core_scaling_benchmark_v0.1.py` renders a deterministic
  PNG/PDF/SVG figure only after the operational gate passes.

The six fixed workloads are soybean Wm82.a2, Brassica napus Da-Ae,
allotetraploid cotton PGCP v2.4.5, maize B73 NAM-5.0, apple GDDH13, and the untouched
Populus JGI v4 validation system.  Their development/diagnostic/confirmatory
roles remain explicit in the registry.  Runtime outcomes are never combined
with biological labels.

The reference freeze requires a full 40-character commit and its `git archive`.
The timed runner extracts and executes that frozen archive, refuses a changed
model, config, conda explicit lock, pip lock, runner, or reference checksum,
and retains every failed replicate in `replicates.tsv`.  Current core commands
do not expose a thread parameter, so the preregistered thread-scaling arm is
reported as not applicable; upstream aligners are not timed as PloidyPatch
core algorithms.
