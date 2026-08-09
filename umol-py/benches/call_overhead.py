#!/usr/bin/env python3
"""Report Python/PyO3 baselines and molecular-fingerprint call costs."""

from __future__ import annotations

import argparse
from functools import partial
import platform
import sys
import timeit

from umol import (
    HashedFingerprintConfig,
    MoleculeAst,
    RefinementRounds,
    StructuralFingerprintConfig,
    NumForm,
)


def empty_python_call() -> None:
    pass


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=1_000_000)
    parser.add_argument("--fingerprint-iterations", type=int, default=10_000)
    parser.add_argument("--repeat", type=int, default=5)
    args = parser.parse_args()

    value_as_lit = NumForm.Lit(0).as_lit
    molecule = MoleculeAst.from_smiles("CCO")
    hashed_wl = partial(
        molecule.hashed_fingerprint,
        config=HashedFingerprintConfig.Wl(
            rounds=RefinementRounds.Fixed(rounds=3),
        ),
    )
    hashed_ecfp = partial(
        molecule.hashed_fingerprint,
        config=HashedFingerprintConfig.Ecfp(radius=2),
    )
    hashed_morgan = partial(
        molecule.hashed_fingerprint,
        config=HashedFingerprintConfig.Morgan(),
    )
    counted_morgan = partial(
        molecule.counted_hashed_fingerprint,
        config=HashedFingerprintConfig.Morgan(),
    )
    structural = partial(
        molecule.structural_fingerprint,
        config=StructuralFingerprintConfig(max_bonds=2),
    )
    cases = [
        ("empty_python_call", empty_python_call, args.iterations),
        ("pyo3_value_as_lit", value_as_lit, args.iterations),
        ("pyo3_hashed_wl", hashed_wl, args.fingerprint_iterations),
        ("pyo3_hashed_ecfp", hashed_ecfp, args.fingerprint_iterations),
        ("pyo3_hashed_morgan", hashed_morgan, args.fingerprint_iterations),
        ("pyo3_counted_morgan", counted_morgan, args.fingerprint_iterations),
        (
            "pyo3_pattern",
            molecule.pattern_fingerprint,
            args.fingerprint_iterations,
        ),
        ("pyo3_structural", structural, args.fingerprint_iterations),
    ]

    print(f"python={sys.version.split()[0]} platform={platform.platform()}")
    for name, function, iterations in cases:
        samples = timeit.repeat(
            function,
            number=iterations,
            repeat=args.repeat,
        )
        nanoseconds = min(samples) / iterations * 1_000_000_000
        print(f"{name}: {nanoseconds:.1f} ns/call ({iterations} iterations)")


if __name__ == "__main__":
    main()
