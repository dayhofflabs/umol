### Short answer
- Build a grammar-based test generator (bounded-depth) from the EBNF to produce many valid strings and assert the implementation accepts them.
- Construct negative tests by mutating valid strings to violate specific rules and assert rejection.
- Verify lexer with maximal‑munch equivalence tests (overlapping tokens) and token-by-token comparisons.
- Add invariant/property tests for semantic constraints (rings, directed bonds, bracket field rules, chirality ranges).
- Differential-test against a reference parser or another implementation for acceptance/rejection where semantics align.
- Fuzz with grammar-aware and token-level fuzzers; gate with coverage and CI.

### Practical process (sequenced)
1) Freeze a spec version
- Tag `opensmiles-umol.md` with a version (done).
- CI: any change to Tokens/Grammar/Semantics triggers the full suite.

2) Lexer conformance
- Deterministic tokenization tests for all ambiguous inputs (e.g., C vs Cl; : vs ::; %12 vs %1 2).
- Property tests: generate random token-like strings; check maximal‑munch and tie-breaking match the spec.
  - Curate a small, high-signal token corpus for maximal‑munch and tie-breaking: multi-char atoms (`Cl`, `Br`, `se`, `as`), percent-ring vs digits (`%12`, `%01`, `%1 2`), chirality (`@`, `@@`, `@TH1`, `@OH30`, out-of-range), charges (`+`, `-`, `++`, `--`, `+0`, `+07`, `+123`), bond runs (`:`, `::`, `//`, `\\`), bracket forms (`[C]`, `[CH3]`, order/multiplicity violations), whitespace (only trailing allowed). 
  - Write deterministic tokenization tests asserting exact token sequences (and spans) for these inputs.
  - Add a lightweight property test that generates short token-like strings from the token alphabet and asserts: no overlaps, full coverage of the string, and re-lexing of the concatenated token text yields the same token sequence (idempotence). Keep bounded length (<= 32).

3) Grammar acceptance
- Grammar-driven generator (bounded depth/width) from the EBNF:
  - Generate thousands of valid SMILES; assert `MoleculeParser` accepts all.
  - Keep seeds/artifacts for regressions.
- Negative generator:
  - Systematically delete/permute tokens to break specific productions (missing ‘]’, extra bond, inter-token whitespace) and assert rejection.

- Seeded valid generator (bounded depth): programmatically compose small `chain`s from a constrained set:
  - atoms: `C`, `N`, `c`, `*`, `[CH]`, `[C@H+]`, `[13C-:0]`
  - bonds: `-`, `=`, `:`, `/`, `\\`
  - rings: digits `1..3`, percent `%12`
  - branches: depth ≤ 2, width ≤ 2
  - connectors: optional `bond` or `.` between nodes
- Assert parse success and simple invariants (no open rings, bond count ≥ chain length-1, no self-loops).
- Negative generator (mutations of valid strings): single, targeted breaks per sample:
  - delete a ring closer; duplicate ring index to force self-loop; `12` cycle; conflicting `/` vs `\\` on ring closure; trailing bond; inter-token whitespace insertion; bracket field reorder/duplication; chirality out-of-range; `%0`/`%01`; negative class; `H10`.
- Start with a deterministic suite (dozens to low hundreds). Add a small property-based set after it’s stable.

4) Semantic constraints verification
- Map each normative rule to explicit tests (many already exist): ring closure rules, direction precedence/conflict, no self-loop/2-cycle, bracket field order/multiplicity, chirality ranges, percent-ring bounds, class/charge/hydrogen numeric behavior.
- Property tests with invariants post-parse:
  - No open rings; no self-loop bonds; two-member rings absent; directed-bond conflicts not present; bracket fields normalized.

5) Differential testing
- Reference parser for syntax-only (optional): transliterate EBNF to a PEG/`pest`/`pom` grammar and compare accept/reject with `MoleculeParser` on large corpora.
- Cross-impl sanity (optional): RDKit/Open Babel accept/reject for purely syntactic areas (not all semantics will align).

6) Fuzzing
- Token-level fuzzer (mutate token streams) and grammar-aware fuzzer (generate via EBNF).
- Run with `cargo-fuzz`; add dictionary entries for tokens (“Cl”, “Br”, “%12”, “[”, “]”, “@@”, etc.).
- Keep crashing/semi-ambiguous cases as regression tests.

7) Coverage and gating
- Measure coverage; ensure every production and semantic rule has hits.
- CI gates:
  - All generators/negative tests/differentials/fuzz seeds pass.
  - Coverage thresholds for parser and semantic checks.

8) Drift checks between spec and implementation
- Script to diff Tokens table in `opensmiles-umol.md` vs the lexer token set.
- Script to check grammar nonterminals vs LALRPOP production names/terminals.
- Fail CI on drift.

9) Release criteria for 1.0
- No red tests across the whole suite.
- Independent parser or tool passes the suite on the same corpus.
- No planned breaking changes over a stabilization window.

10) Error handling. Feasible and worth it. Plan for a small, stable error taxonomy with machine-readable codes, precise spans, and deterministic wording tied directly to the spec. A solid first release is ~2–3 weeks; a polished, “industry-quality” set with full coverage, docs, and a CLI linter is ~4–6 weeks for one engineer.

What to define
- Error taxonomy (stable codes):
  - LEX (lexer): invalid token, unterminated bracket, inter-token whitespace, bad percent form.
  - SYN (grammar): unexpected token, trailing bond, missing ‘]’, branch/dot misuse.
  - RING: unclosed ring, conflicting directions, self-loop, two-member ring, percent-range errors.
  - NUM: overflow, negative class, invalid H count, chirality range.
  - BRKT: field order, duplicate fields, adjacency issues (e.g., H10), empty class.
  - STEREO: inconsistent double-bond markers, ambiguous/insufficient markers where the spec requires unknown.
- Diagnostic schema (machine-usable):
  - code (string, stable), category, message (short), span {start,end}, primary token, expected/got, details (e.g., ring_index, atom_index), help (single short suggestion).
- Normative registry:
  - Put the canonical list in spec/errors.md (or a “Diagnostics” section in the spec) with: Code, Title, Condition (normative), Message template, Parameters, Example (valid/invalid), Help.
- Parser/lexer integration:
  - Replace ad hoc errors with structured diagnostics. Map every reject site to a specific code and fill parameters.
  - Keep a default “internal” code for unforeseen failures; avoid leaking it in user-visible paths.

Effort estimate
- Week 1:
  - Inventory all current rejects; map to draft codes; add spans and parameters in lexer and parser; wire through a Diagnostic type.
  - Write table-driven tests asserting code+span for each invalid case already in tests.
- Week 2:
  - Fill gaps per spec (every MUST reject must have a code). Add targeted invalid tests for each new code.
  - Add a simple CLI “smiles-lint” mode emitting JSON diagnostics and pretty text.
- Weeks 3–4 (quality pass):
  - Fuzzing hooks to ensure every reject path yields a structured code.
  - Documentation pages generated from the registry; examples verified by tests.
  - CI gates: 1:1 mapping between registry and implementation (no orphan codes, no unreachable codes).

Initial code set (indicative)
- LEX: LEX_INVALID_TOKEN, LEX_INTERTOKEN_WHITESPACE, LEX_BAD_PERCENT_LEADING_ZERO, LEX_BAD_PERCENT_RANGE
- SYN: SYN_UNEXPECTED_TOKEN, SYN_TRAILING_BOND, SYN_UNCLOSED_BRACKET, SYN_BAD_BRANCH_DOT
- RING: RING_UNCLOSED, RING_CONFLICT_DIR, RING_SELF_LOOP, RING_TWO_MEMBER
- NUM: NUM_OVERFLOW, NUM_CLASS_NEGATIVE, NUM_HCOUNT_BAD, NUM_CHIRAL_OUT_OF_RANGE
- BRKT: BRKT_ORDER, BRKT_DUP_FIELD, BRKT_HCOUNT_TWO_DIGITS, BRKT_EMPTY_CLASS
- STEREO: STEREO_DOUBLE_CONFLICT, STEREO_DOUBLE_INSUFFICIENT

Parsing tests:
- Relate mutations producing invalid strings tied to a normative rule so every MUST-reject is covered.

Placement and format
- Spec: `umol-models-graph/spec/errors.md` (normative registry).
- Implementation: `umol-models-graph/src/diagnostics.rs` (Diagnostic struct, codes enum with stable string repr), plus mappings in `lexer.rs`, `state.rs`, and parser actions.

Output formats
- Human: single-line message plus one help line.
- JSON: code, message, span, category, details. Deterministic wording is important for downstream linters/IDEs.

- Triage loop:
  - Failures imply either spec gaps or implementation gaps. For each, decide: fix code, refine spec, or adjust generator constraints.
  - Track a checklist mapping each normative rule in `spec/opensmiles-umol.md` to at least one passing test (valid and/or invalid).

- Optional quick win:
  - Add a debug mode in the parser to dump token spans and key decisions (e.g., ring open/close, chosen bond order). That will speed triage without committing to a full diagnostics framework yet.

Once this is green, you’ll have the “basic set of rejects” and acceptances to seed the error taxonomy (step 10).

