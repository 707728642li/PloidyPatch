# Coffea reveal execution patch v1.2

The stable Coffea blind project completed custody before target labels were
interpreted. Reveal preflight exposed two packaging defects before the
evaluator was invoked: historical root checksum commands omitted nested
`SHA256SUMS` files, and the reveal orchestrator rejected the normal conda
`bin/python` symlink even when its resolved executable remained inside the
frozen environment.

The canonical frozen artifact-manifest writer replaces the historical shell
root writer. The reveal environment check now resolves `bin/python`, requires
an executable regular file, and requires the resolved path to remain within
the frozen conda prefix. The complete environment lock is still verified.

No candidate, reference, perturbation, H1 endpoint, bootstrap seed, threshold,
truth definition, or statistical rule changes. No builder, evaluator, ranker,
or performance calculation ran before this patch. The existing blind custody
and reveal authorization are checksum-bound; this is a reveal-execution repair
only, and the next evaluator result is retained in its resulting direction.
