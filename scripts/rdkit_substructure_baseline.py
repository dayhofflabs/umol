#!/usr/bin/env python
"""Actual-RDKit substructure-matching baseline for the basic_opensmiles corpus.

Non-permanent dev oracle (not a build dependency). Mirrors the Rust harness in
`umol-graph/benches/substructure.rs`: same corpus, same three patterns, and the
same matching semantics so the timings are directly comparable.

Semantics alignment with `MoleculeAst::substructure_matches`:
  - Atoms are element-only (`[#6]`/`[#7]`/`[#8]`), matching aromatic and aliphatic
    alike, like an element-only `AtomAst`.
  - All bonds are "any" (`~`), like a `BondAst(Undetermined)`.
  - All embeddings are enumerated (`uniquify=False`), like the Rust matcher, not a
    first-hit / symmetry-deduplicated count.

Run with the dedicated env:
    ~/.micromamba/envs/rdkit-ref/bin/python scripts/rdkit_substructure_baseline.py
"""

import os
import statistics
import time

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")  # silence per-molecule sanitization warnings

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CORPUS_DIR = os.path.join(
    REPO_ROOT, "umol-io", "tests", "smiles_parsing", "data", "basic_opensmiles"
)

# (name, SMARTS) — topology identical to the Rust patterns.
PATTERNS = [
    ("branched", "[#6](~[#6])~[#6](~[#6])~[#7]"),
    ("phenol", "[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#8]"),
    ("bicyclic", "[#6]~1~[#6]~2~[#6]~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~1"),
]

REPS = 5
MAX_MATCHES = 1_000_000


def load_corpus():
    """Parse the second line (line 1 is a `#` comment) of every .smiles file."""
    mols = []
    rejected = 0
    for root, _dirs, files in os.walk(CORPUS_DIR):
        for fname in files:
            if not fname.endswith(".smiles"):
                continue
            with open(os.path.join(root, fname), encoding="utf-8", errors="replace") as fh:
                lines = fh.read().splitlines()
            if len(lines) < 2 or not lines[1]:
                continue
            mol = Chem.MolFromSmiles(lines[1])
            if mol is None:
                rejected += 1
            else:
                mols.append(mol)
    return mols, rejected


def main():
    corpus, rejected = load_corpus()
    print(f"corpus: {len(corpus)} parsed, {rejected} rejected by RDKit")
    print(f"reps: {REPS}, uniquify=False, maxMatches={MAX_MATCHES}\n")

    queries = [(name, Chem.MolFromSmarts(smarts)) for name, smarts in PATTERNS]
    for name, query in queries:
        assert query is not None, f"bad SMARTS for {name}"

    print("| pattern  | hits | ms/pass (median) | ms/pass (min) |")
    print("|----------|------|------------------|---------------|")
    for name, query in queries:
        hits = sum(
            len(m.GetSubstructMatches(query, uniquify=False, maxMatches=MAX_MATCHES))
            for m in corpus
        )
        timings = []
        for _ in range(REPS):
            start = time.perf_counter()
            for m in corpus:
                m.GetSubstructMatches(query, uniquify=False, maxMatches=MAX_MATCHES)
            timings.append((time.perf_counter() - start) * 1000.0)
        print(
            f"| {name:<8} | {hits:>4} | {statistics.median(timings):>16.1f} "
            f"| {min(timings):>13.1f} |"
        )


if __name__ == "__main__":
    main()
