#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then echo "usage: $0 PROJECT_ROOT" >&2; exit 2; fi
project_root=$(realpath "$1")
code_root=$project_root/code
python_bin=$project_root/envs/ploidypatch-dev/bin/python
validation=$project_root/results/natural/maize_v2/validation/isoseq_v0.1
rank_root=$project_root/results/natural/maize_v2/discovery/homeolog_ranker
evidence=$validation/evidence/candidate_isoseq_evidence.tsv
rankings=$rank_root/natural/review_rankings.tsv
result_root=$project_root/results/natural/maize_v2/statistics/isoseq_review_v0.1
working_root=${result_root}.working
for required in "$python_bin" "$validation/SHA256SUMS" "$rank_root/SHA256SUMS" \
                "$evidence" "$rankings"; do
    [[ -s $required ]] || { echo "missing maize Iso-Seq statistic input: $required" >&2; exit 1; }
done
(cd "$validation" && sha256sum -c SHA256SUMS >/dev/null)
(cd "$rank_root" && sha256sum -c SHA256SUMS >/dev/null)
if [[ -e $result_root || -e $working_root ]]; then
    echo "refusing to overwrite maize Iso-Seq statistics" >&2; exit 1
fi
mkdir -p "$working_root"
{
    printf 'field\tvalue\ncode_commit\t%s\n' \
        "${PLOIDYPATCH_CODE_COMMIT:-unavailable_server_mirror}"
    printf 'replicates\t20000\nseed\t20260808\nalpha\t0.05\n'
    printf 'primary_state\tfull_chain_supported\nprimary_budget\t100\n'
    printf 'resampling_unit\ttarget_chromosome\n'
    printf 'random_null\twithin_chromosome_observed_top_k_allocation\n'
} > "$working_root/run_contract.tsv"
cd "$code_root"
"$python_bin" -m ploidypatch.cli evidence bootstrap-isoseq-review-yield \
    --evidence "$evidence" --review-rankings "$rankings" \
    --review-budget 25 --review-budget 50 --review-budget 100 --review-budget 200 \
    --replicates 20000 --seed 20260808 --alpha 0.05 \
    --output-json "$working_root/bootstrap.json" \
    > "$working_root/stdout.json" 2> "$working_root/stderr.log"
[[ -s $working_root/bootstrap.json ]] || { echo "missing Iso-Seq bootstrap" >&2; exit 1; }
(
    cd "$working_root"
    find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS
    sha256sum -c SHA256SUMS >/dev/null
)
mv "$working_root" "$result_root"
printf 'maize Iso-Seq statistics frozen: %s\n' "$result_root"
