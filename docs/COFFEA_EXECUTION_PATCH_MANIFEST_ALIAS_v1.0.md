# Coffea execution-only patch: duplicate role-manifest alias

The first isolated Coffea blind attempt exited with status 1 before candidate
input preparation or any projection method. The role root was mounted with the
same immutable `role_manifest.json` at two namespace paths. The strict
exact-universe checksum verifier correctly rejected the unlisted alias
`blind_role_manifest.json`.

This patch removes only that redundant namespace alias and records the canonical
`/holdout/role_manifest.json` in the namespace audit. It does not change the
protocol, input bytes, truth definition, candidate adapters, pool rules,
endpoints, thresholds, bootstrap settings, or failure policy. Evaluator truth
construction had completed before the patch; no blind candidate, formal score,
label, or performance result existed.
