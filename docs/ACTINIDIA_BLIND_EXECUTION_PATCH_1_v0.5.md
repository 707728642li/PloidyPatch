# Actinidia blind execution patch 1 (v0.5)

Status: pre-score, pre-label and pre-reveal execution-only amendment

This amendment records a deterministic provider-GFF compatibility failure in
the first isolated Actinidia blind attempt. It does not change the frozen
species roles, source inputs, Red5 Ad-alpha truth construction, perturbation,
candidate methods, candidate policy, model, features, thresholds, endpoints or
safety gates.

## Superseded execution and retained attempt

- base code commit: `5d598cd394a8bb111bf31c8c4e87ed8d6a6136d0`
- base protocol `SHA256SUMS` SHA-256:
  `4cab1d6aa2eefcff4ce2a3d1d3fdecf8fd8d9dfc2fd6d506dd18cf1ec1a3bb8c`
- base execution `SHA256SUMS` SHA-256:
  `6c2cffc0d929b3d8128b75fa68a1e7b04ccf2fbc39af3123adf02b78c48aa82b`
- retained failed attempt:
  `results/blind_runs/actinidia_red5_v0.5_run1.working`
- exit status: `1`
- no candidate pool, score table, custody manifest, reveal authorization,
  label table or performance result was produced.

The failure occurred before homology projection. A read-only audit of all
285,988 direct child links in the normalized public A. eriantha GFF found
exactly three out-of-bounds rows, all CDS start underflows and all on the same
seqid and strand as their parents. `DTZ79_05g01990:CDS` at
`Chr05:2960679-2960853` precedes its provider gene/mRNA start `2960873` by 194
bp. Two `DTZ79_15g11550:CDS` rows precede their provider start `16224283` by
138 and 75 bp. The minimal repair changes four parent starts (two mRNA and two
genes). There are no cross-seqid, cross-strand, missing or ambiguous links.
A. rufa has 594,080 direct child links and zero bound violations, so the rule
is a no-op for that reference. The frozen compatibility adapter correctly
rejected the first invalid child.

## Exact repair

Patch 1 adds an explicit `--repair-parent-bounds` compatibility mode. For each
parent, it may only expand the start/end columns to the union of its original
same-seqid, same-strand structural children. It may never shrink a parent or
change any CDS, UTR or exon row, coordinate, phase, identifier or parentage.
Transcript bounds are repaired first; a gene bound may then expand only to the
union of its repaired transcript children. Cross-seqid, cross-strand,
multi-parent or unresolved-parent cases remain hard errors.

The adapter reports every repaired parent and proves byte-stable CDS rows with
input/output SHA-256 digests. The mode is enabled only in the Actinidia
candidate-reference LiftOn compatibility path; default behavior and the
already frozen Populus implementation remain unchanged.

## Patch-freeze conditions

The new execution freeze must bind this file, the immutable base protocol,
base execution, composite model, staged input checksums, retained failed
attempt and an exact changed-file whitelist. It must record
`post_evaluator_truth_failed_blind_pre_candidate_pre_score_pre_label_execution_patch`.
It must reject any failed attempt containing scores, labels, custody or reveal
artifacts and reject any change to scientific configuration, contract, source
roles, model or benchmark bytes.

The original freeze and failed attempt remain preserved. A second isolated
attempt is permitted only from the patch-frozen bytes, before any label or
performance reveal, as specified by the preregistered unrevealed execution
failure policy.
