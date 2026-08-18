#!/usr/bin/env python3
"""Run Pychopper with a reproducible NumPy autotuning sample.

Pychopper 2.7.10 samples reads through ``numpy.random`` when its primer-score
cutoff is autotuned, but its command-line interface does not expose a seed.
This tiny adapter sets that seed and otherwise forwards the argument vector
unchanged to the upstream command.
"""

from __future__ import annotations

import sys
from collections.abc import Sequence


def main(arguments: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if arguments is None else arguments)
    if len(args) < 2:
        raise SystemExit(
            "usage: run_pychopper_seeded.py RANDOM_SEED PYCHOPPER_ARGUMENT ..."
        )
    try:
        random_seed = int(args[0])
    except ValueError as error:
        raise SystemExit("RANDOM_SEED must be an integer") from error
    if random_seed < 0 or random_seed > 2**32 - 1:
        raise SystemExit("RANDOM_SEED must be in [0, 2**32 - 1]")

    import numpy
    numpy.random.seed(random_seed)
    from pychopper.scripts.pychopper import main as pychopper_main

    sys.argv = [sys.argv[0], *args[1:]]
    result = pychopper_main()
    return 0 if result is None else int(result)


if __name__ == "__main__":
    raise SystemExit(main())
