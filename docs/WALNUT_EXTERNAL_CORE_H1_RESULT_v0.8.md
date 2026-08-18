# Walnut external core-H1 result v0.8

## Formal outcome

Walnut run 7 is a checksum-sealed **invalid external run**, not a biological
negative and not a not-evaluable dataset. The frozen evaluator was not
invoked, no H1 estimate or bootstrap interval was produced, and the Walnut
system contributes neither positive nor negative performance evidence.

The blind phase itself completed under the frozen isolation contract. All six
method/reference arms were atomically published, the retain-distinct pool had
73,382 accepted decisions, and the suppress-overlap pool had 52,093 accepted
and 21,289 rejected decisions. The pool manifest verified every listed file,
and the final blind-run tree verified 242/242 SHA256 entries with zero
failures. Custody recorded `truth_mounted=false`, `network_access=false`, the
`core_H1_only_no_ranker` profile, and zero model/ranker authorization.

The reveal builder stopped with the exact reason
`ValueError:Published output working-path lineage differs: final_annotation.gff`.
The GeMoMa output manifest correctly recorded its pre-rename path inside the
blind namespace, `/run/blind-run/project/.../candidate_carya.working/...`.
During host-side reveal, the verifier compared that immutable namespace path
to the host path below the finalized `run7` directory rather than comparing
their exact project-relative lineage. The same verification had correctly
passed inside the blind namespace. This is an execution-interface defect; it
does not concern candidate structure, truth construction, H1 scoring, genome
quality, or the biological applicability of the method.

## Truth custody and retry decision

The authorization artifact has `truth_reveal_authorized=true`. The builder
failed before `score_annotation_repair` read `hidden_truth.json`, but it had
already verified the checksum-sealed evaluator-only tree. We therefore use the
strict interpretation that the truth side was opened. Walnut v0.8 will not be
patched and retried, and no result from its completed blind pools will be
presented as a formal performance estimate.

The namespace-independent verifier correction is limited to future untouched
systems. It accepts only one of two exact path lineages: the native host
`<final>.working/upstream/<role>` path or the frozen
`/run/blind-run/project/<project-relative-final>.working/upstream/<role>`
path. Arbitrary prefixes, relative paths, wrong project-relative locations,
role drift, byte drift and SHA drift remain rejected.

## Immutable evidence

- Blind run root SHA256SUMS SHA256:
  `3d0afcbf1678b3572ab44cabdf22507f1b2640356bf08f16de44fc8fe9ca47a3`
- Blind custody manifest SHA256:
  `c32ccec00f1fea6dc84b90aaadd5353174bf86112b03cb1815ddf7e4bdb36526`
- Pool SHA256SUMS SHA256:
  `cd4558b54d627554b14a526a9cd235e2727ae8f2ecc69626a3b98b2a6030f062`
- Protocol SHA256SUMS SHA256:
  `2ba81a9bf24edf4054d310f23bdf3fac941651b70d7e5ab20fa43982b3d1559b`
- Execution patch-7 SHA256SUMS SHA256:
  `05d8484941f21751854adf10f953f3fbe940ed9af7b5626db28f750baeec2c6c`
- Reveal SHA256SUMS SHA256:
  `7d23135281ad626d5befec42cbc59e3ef613e05c434fdb8285c824c0cb03ae6e`
- Evidence copy:
  `manuscript/source_data/walnut_external_v0.8_invalid_run7/`

The copied reveal root and its nested reveal-input tree both pass their
original exact-universe SHA256SUMS manifests.
