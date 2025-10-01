## Strategies for conformance testing

See 34-opensmiles-formal-grammar-2025-09-13.md for the previous iteration.

### General

- Differential testing: compare parse+IR/canonical SMILES against RDKit/Open Babel/ChemAxon on large corpora (PubChem/ChEMBL). Flag divergences; use as oracle.
- Grammar-based generation: generate valid-by-construction SMILES from the spec grammar (with constraints), and a companion “anti-grammar” for targeted invalids (e.g., unmatched rings, illegal bracket fields).
- Metamorphic testing: assert invariants under semantics-preserving edits:
  - Ring index renumbering/permutation → same IR
  - Branching parentheses reassociation (that doesn’t change connectivity) → same IR
  - Bracket field reordering → same IR
  - Inserting allowed whitespace/comments (lenient mode) → same IR
  - Pretty-printer idempotence: canon(parse(x)) == canon(parse(pretty(x)))
- Round-trip/canonicalization oracles:
  - parse → canonicalize → render → parse equivalence
  - parse → serialize (canonical) → parse → IR isomorphic
- Property checks on IR (parser-agnostic): connectivity, no dangling rings, bond endpoints in range, component counts, atom/bond count bounds, charge constraints (when valence rules added).
- Coverage-guided generation: grammar- or token-aware fuzzers (libFuzzer/AFL++ with custom mutator) plus coverage reports (llvm-cov/grcov) to drive inputs into rarely hit branches; distill corpora (afl-cmin).
- Mutation testing: mutate parser/linter rules and ensure tests fail; helps find missing assertions.
- Conformance corpus curation: curate minimal positive/negative exemplars for each spec clause; auto-generate parametric tables (e.g., ring sizes, bracket permutations, percent forms).
- Differential fuzzing: feed the same fuzzer stream to multiple engines (yours, RDKit, Open Babel); triage only disagreements.
- Snapshot (golden) diagnostics: stable expected lint/error codes + spans on a curated suite; changes require deliberate review.

### Use now (low lift, high value)
- Differential testing (CLI round-trips)
  - Parse/print with umol vs RDKit/OpenBabel/ChemAxon; snapshot deltas (success/failure, canonical SMILES, error spans).
- Snapshot testing with corpora
  - Insta suites over curated and external corpora (PubChem/ChEMBL slices). You already planned this; wire it up first.
- Metamorphic property tests
  - Ring renumbering/permutation; bracket field reordering; whitespace/comments (lenient flags); pretty-printer idempotence.
- IR invariants
  - One-pass assertions after parse: valid endpoints, no dangling rings, component integrity, simple charge/isotope bounds.
- Coverage-guided fuzz extension
  - Keep proptest; add a cargo-fuzz target with a token-aware mutator; gate on “no panics” + error-span in-bounds.

### Next (medium effort, strong assurance)
- Differential fuzzing
  - Feed identical fuzz streams to umol and RDKit/OpenBabel; snapshot only disagreements (oracle by consensus).
- Grammar-based generation
  - Valid-by-construction generator from a pared-down OpenSMILES grammar; small “anti-grammar” invalids for structure-negative cases.
- Canonicalization oracles
  - parse → canonicalize → render → parse equivalence; require graph isomorphism (node relabeling) rather than string equality.
- Mutation testing in CI
  - Use cargo-mutants (or mutagen) on parser/linter hotspots; require tests fail for mutants to catch gaps.
- Coverage reporting gates
  - llvm-cov/grcov in CI; track line/branch coverage; require coverage deltas non-negative on PRs.

### Later (heavier lift or diminishing returns)
- Full formal semantics/verification
  - Proof-oriented specs or model checking are likely overkill here.
- Symbolic execution/SMT or exhaustive combinatorics
  - Not practical given SMILES surface area.

### How to stage this in the repo
- Start by adding:
  - Differential harness binaries: “compare-rdkit”, “compare-obabel” that emit insta snapshots for a folder.
  - Metamorphic rstest/property modules for: ring renumbering, bracket field reorder, whitespace/comments, pretty printer idempotence.
  - IR invariant checks behind a feature flag you can call in tests.
  - cargo-fuzz target (retain proptest).
- Then add:
  - Differential fuzz job (nightly or CI optional job) with a triage artifact of mismatches.
  - Small grammar generator crate (valid + targeted invalid).
  - Coverage job and a simple quality gate.
