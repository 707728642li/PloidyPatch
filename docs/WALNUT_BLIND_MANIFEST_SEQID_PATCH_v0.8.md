# Walnut blind manifest primary-seqid serialization patch v0.8

Walnut blind run 2 completed target and candidate-reference normalization and
all structural compatibility checks, then stopped before any candidate method
executed. The normalized-input manifest attempted to JSON-serialize the
`frozenset` returned by the generic primary-seqid table reader.

This execution patch returns the primary sequence IDs as a plain list in the
already frozen two-column seqid-table order. The same frozenset remains the
membership gate for FASTA extraction, and the extracted FASTA and FAI bytes
are unchanged. A regression test requires deterministic table order and
successful JSON serialization.

No reference, sequence, annotation, candidate method, event, threshold,
sampling rule, hypothesis or success gate changes. Run 2 is retained; no
projection, candidate pool, score or truth label existed at failure time.
