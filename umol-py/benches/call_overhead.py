#!/usr/bin/env python3
"""Report empty-Python and trivial PyO3 call costs for binding comparisons."""

from __future__ import annotations

import argparse
import platform
import sys
import timeit

from umol import ValueAst


def empty_python_call() -> None:
    pass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=1_000_000)
    parser.add_argument("--repeat", type=int, default=5)
    args = parser.parse_args()

    value_as_lit = ValueAst.Lit(0).as_lit
    cases = [
        ("empty_python_call", empty_python_call),
        ("pyo3_value_as_lit", value_as_lit),
    ]

    print(f"python={sys.version.split()[0]} platform={platform.platform()}")
    for name, function in cases:
        samples = timeit.repeat(
            function,
            number=args.iterations,
            repeat=args.repeat,
        )
        nanoseconds = min(samples) / args.iterations * 1_000_000_000
        print(f"{name}: {nanoseconds:.1f} ns/call")


if __name__ == "__main__":
    main()
