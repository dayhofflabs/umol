## OpenSMILES parser optimization: findings and plan (2025-09-21)

### Current findings (bench, macOS)
- lex_only/short ≈ 47 ns
- parse_minimal/short ≈ 1.57 µs (~33× lex)
- parse_only (Full)/short ≈ 3.62 µs (~77× lex)
- lint_only/short ≈ 16.9 µs
- lint_plus_parse_fast/short ≈ 18.8 µs (≈ lint_only)

Implications:
- LALRPOP grammar overhead (actions only) is already ~30× lex; this is the lower bound for this stack on our inputs.
- Full parse adds IR and late passes to reach ~80× lex.
- Lint dominates combined runs; however, parser competitiveness is primary, so we target parsing first.

### Hypotheses for LALRPOP costs
- LR driver overhead: table lookups, shift/reduce dispatch, stack ops.
- Action function call overhead and state method calls (even when minimal).
- Token-to-nonterminal mapping and lexer→parser handoff cost.

### Immediate measurements to run
- Flamegraph (Instruments on macOS or cargo-flamegraph) for parse_only and parse_minimal.
  - Identify hot symbols: state methods (bump_*, link_to, ring open/close), lalrpop runtime functions, Vec push/clone paths, HashMap.
  - Compare minimal vs full to isolate IR and late-pass contributions.

### Options forward
1) Optimize within LALRPOP
   - Remove or inline hot state calls; reduce logging (already off).
   - Replace HashMap in rings with smallvec path in Full (mirroring lint-fast) if safe.
   - Minimize allocations in actions; reuse small buffers.
   - Consider grammar refactors to reduce reductions/ambiguity.

2) Alternative parser infrastructures
   - Hand-rolled streaming recursive-descent/state machine with a compact stack.
   - Pratt/precedence-style where applicable (bonds/branches/rings are regular, not expressions; likely a custom FSM is best).
   - Other parser generators (e.g., LR with codegen focused on speed, or PEG with zero-copy; evaluate carefully).

3) Hybrid approach
   - Keep current lexer; build a manual SMILES FSM parser for chain/branch/ring handling, producing IR directly.
   - Reserve small fixed-size stacks for branches/components and a tiny ring index store.

### Success criteria
- Parser-only ≤ 2–4× lex on short/medium inputs.
- Identical surface behavior and diagnostics (modulo source switch from parser to linter where applicable).

### Next steps
1. Record flamegraphs for parse_only and parse_minimal; annotate hotspots.
2. Prototype a micro FSM for StartNode/Node/Branch/Ring over the existing lexer; benchmark vs parse_minimal.
3. If FSM approaches ≤ 4× lex, plan migration of grammar rules; if not, revisit LALRPOP tuning.

### Flamegraph analysis

- parse_minimal/short (~1.56 µs): dominated by the LALRPOP driver (shift/reduce and grammar action trampolines). Very little shows outside the parser loop; lexer cost is tiny. Our Minimal mode’s action counter (“tick”) appears but is a small slice.
- parse_only/short (~3.68 µs): all of the above plus IR/state work:
  - Frequent small calls: ParseState::bump_atom_idx / bump_bond_idx, stage_bond_*(), link_to().
  - IR pushes: push_atom(), push_resolved_bond() with Vec growth/field writes.
  - Ring bookkeeping: open_ring()/close_ring(); HashMap lookups/removes show up repeatedly.
  - Component/branch operations appear but are minor on short inputs.
- Takeaway:
  - LALRPOP baseline (driver + actions) already ~33× lex; that’s your hard floor on this stack.
  - IR/state doubles that for parse_only on short inputs.
  - Hitting 2–4× lex with LALRPOP is not realistic; a tight hand-rolled FSM (smallvec ring store, compact branch stack, no HashMap in hot paths, inlined defaults) is the viable path if the 2–4× target is firm.

## FSM parser plan

### Goals
- Single-pass streaming FSM over the existing lexer, targeting ≤ 2–4× lex latency.
- Minimal hot-path work: small branch stack, small ring store, pre-reserved atom/bond buffers.
- Defer semantics/stereo to optional post-passes.

### Milestones (M0→M6)
- M0: Chain-only, bare atoms, implicit single bonds. Bench vs lex-only.
- M1: Add branches ( and ).
- M2: Add explicit bonds (- = # : / \\) with staged bond state.
- M3: Add rings (digits/%), smallvec ring store, conflict/self/two-member checks.
- M4: Components (.).
- M5: Brackets: implement a tight bracket scanner for [ ... ].
- M6: Optional semantics (aromatic defaults, stereo post-passes).

### Broader benchmark (to add)
- Corpus:
  - Short typical: aromatic ring, one branch, one ring index.
  - With brackets: several bracket atoms, charges, hydrogens.
  - Rings-heavy: multiple closures including %xx.
  - Long linear chain (no branches) and long branched chain.
  - Edge cases: empty branches, dangling bonds, unclosed rings (still valid to parse until error).
- Metrics:
  - ns per atom, ns per token, total µs per string.
  - Compare lex_only vs FSM_M0..M5.
  - Include a long-run non-Criterion profile runner for flamegraphs.

### Lexer plan (bytestrings)
- Keep Logos for now but operate on &[u8] and ASCII-only paths per spec.
- Avoid String conversions; ensure tokens carry byte spans.
- Later, consider a lightweight bespoke scanner for the hot tokens.

### Notes
- After adding the benchmark corpus and adjusting the lexer interface to bytes, proceed with M0→M1 and re-measure.

## IR shape and manipulation notes

- Primary IR remains flat arrays: `atoms[]`, `bonds[]` with local indices; best for parser hot path and MOL unification.
- Do not embed parser-transient state (branch/ring stacks) into IR.
- For manipulation, layer views over arrays:
  - Adjacency/degree built lazily (CSR-like) and cached until invalidated.
  - Stereo/ring annotations as optional side arrays, populated by semantic passes.
- Lazy conversion: parse fills Atom/Bond fields; derived structures are optional and on-demand.
- Alignment: list-based IR suits both SMILES and MOL; graph conveniences should be layered, not baked in.

## Corpus (M0) adjustments

- Tiered to milestones. For M0 (chain-only, bare atoms):
  - Baseline: empty molecule.
  - Chain lengths: 1, 5, 10, 50, 100, 1000.
  - Two sets:
    - Same element: repeated `C`.
    - Organic-only mix (CHNOPSFClBrI) with target frequencies; omit bare `H` for M0.
  - Defer “all elements 1–118” to M5 (requires brackets).

## FSM parser implementation

### Milestones (M0→M6)
- M0: Chain-only, bare atoms, implicit single bonds. Bench vs lex-only.
- M1: Add branches ( and ).
- M2: Add explicit bonds (- = # : / \\) with staged bond state.
- M3: Add rings (digits/%), smallvec ring store, conflict/self/two-member checks.
  - M4: Components (.).
- M5: Brackets: implement a tight bracket scanner for [ ... ].
- M6: Optional semantics (aromatic defaults, stereo post-passes).


### M0 Plan

- Scope
  - Chain-only, bare organic atoms: B, C, N, O, P, S, F, I, and two-letter halogens Cl, Br
  - Implicit single bonds; single component; ASCII-only; no whitespace handling
  - Any other byte triggers “unsupported at M0” with position

- API
  - `parse_smiles_m0(input: &[u8]) -> Result<Molecule, M0Error>`
  - `M0Error { kind: UnsupportedToken, pos: usize }` (minimal set for M0)
  - Builds IR via `io::ir::builder::MoleculeBuilder` directly (no remap)

- FSM design (hot path)
  - Single pass over bytes with index `i`; track `last_atom_idx: Option<u32>`
  - Recognize atoms quickly:
    - 1-byte: B,C,N,O,P,S,F,I
    - 2-byte: Cl, Br (1-byte lookahead)
  - On atom:
    - If `last_atom_idx.is_some()`: `on_bond(last, curr, Single)`
    - Always `on_atom(element)` → set `last_atom_idx = Some(curr_idx)`
  - On any other byte: return error with `pos = i`
  - At end: `finish()` molecule

- Diagnostics
  - Only unsupported-token errors at a byte offset for out-of-scope features
  - Empty input returns empty molecule OK

- Benchmarks
  - Add `parse_m0_chain` group to `opensmiles_parsing.rs` mirroring chain corpora (1,5,10,50,100,1000)
  - Compare against `lex_only` and `parse_minimal`
  - Target: 2.5–5× lex on chain corpora

- Tests
  - Valid: `C`, `CC`, long chains, `CClC`, `CBrC`
  - Unsupported: `(`, `)`, `[`, `]`, digits/% rings, bond symbols `-=:#/\\`, `@`, `.`, lowercase, whitespace
  - Assert atom/bond counts and implicit single bonds

- Integration
  - New module: `io/smiles/fsm_m0.rs`; export via `io::smiles::fsm_m0` and re-export `parse_smiles_m0`
  - Keep legacy parser/linter on `lexer_old` during transition; remove after M2 parity

### M0 Tasks

1. Create `io/smiles/fsm_m0.rs` with byte FSM and `parse_smiles_m0(&[u8])`
2. Wire `MoleculeBuilder` calls: `on_atom`, `on_bond`, `finish`
3. Export module in `io/smiles.rs` and re-export the function
4. Add Criterion group `parse_m0_chain` to `benches/opensmiles_parsing.rs`
5. Add tests under `tests/opensmiles_parsing/` for valid chains and unsupported inputs
6. Benchmark vs `lex_only` and `parse_minimal`; record results in this doc
7. Decide go/no-go for M1 (trees) based on hitting ≤5× lex on chains

### M1 Tasks

1. Extend FSM to handle '(' and ')' with a small branch stack.
2. Connect branch bonds via MoleculeBuilder (push/pop attach points).
3. Add tests for single and nested branches; keep rings unsupported.
4. Add benches for branch inputs; compare vs current baseline.
5. Gate to M2 (bonds) if ≤5× lex on branch cases.

### M2 Tasks

1. Extend FSM with bond token handling ( - = # $ : / \ ).
2. Add rstest unit tests for bonds and errors (trailing/consecutive).
3. Add benches for M2 bonds vs M1 on chains/branches.
4. Update spec and errors: bond semantics, trailing-bond error.
5. Integrate style lints for explicit single and aromatic ':'.
6. Decide go/no-go for M3 (rings) based on ≤5× lex with bonds.

### M3 Plan

- Scope and goals
  - Add ring tokens and semantics on top of M2 (bonds), keeping M3 a strict superset of M2/M1.
  - Support single- and multi-digit ring indices, ring bonds, rings within branches/groups, fused and spiro junctions.
  - Emit precise, positionful errors for ring issues.

- Grammar additions (OpenSMILES-aligned)
  - Ring index tokens:
    - DIGIT ring indices: 1–9 following an atom.
    - Multi-digit ring indices: `%DD` (two digits). Decide on `%DDD` (100–999) support; default to two digits only unless needed.
  - Optional preceding bond symbol on either side of a ring index: `- = # $ : / \`
  - Multiple ring indices allowed on the same atom (e.g., fused `...12`).
  - Rings are allowed inside branches/groups.

- Error taxonomy (new M3 errors; keep all M2 errors unchanged)
  - SYN_RING_UNCLOSED { pos_open }: ring index seen once and not matched by EOI.
  - SYN_RING_SELF_LOOP { pos }: ring closure would connect an atom to itself.
  - SYN_RING_TWO_MEMBER { pos }: ring closure would connect immediately adjacent atoms to form a 2-member ring.
  - SYN_RING_BOND_CONFLICT { pos, open_pos }: both sides specified explicit bond types that disagree.
  - SYN_RING_DIR_CONFLICT { pos, open_pos }: directional bond markers across a ring that are inconsistent or lack a double-bond context.
  - LEX_RING_INDEX_INVALID { pos }: malformed ring token (`%` not followed by two digits, leading zeros rules if any).
  - Optional (decide): SYN_LEADING_RING { pos } if a ring token appears with no current atom.

- Implementation steps
  - M3-1: Create `io/smiles/fsm_m3.rs`; export in `io/smiles.rs` alongside M2.
    - Acceptance: module compiles; function signature mirrors M2 parse entry point; gated behind feature/tests only at first.
  - M3-2: Ring token lexing in FSM
    - Recognize DIGIT ring indices and `%DD`.
    - Capture optional preceding bond token and its position.
    - Acceptance: unit tests for tokenization-only surfaces pass (invalid `%`, out-of-place tokens).
  - M3-3: Ring open/close state table
    - Data: `Vec<Option<OpenRing>>` mapping ring index → open entry.
    - `OpenRing` fields: `atom_id`, `bond_opt`, `pos_open`, `dir_opt`, `aromatic_opt`.
    - Acceptance: opening stores entry; closing retrieves and clears entry.
  - M3-4: Close logic and bond precedence
    - Combine bond types: if both explicit and equal → use it; if both explicit and different → conflict error; if one explicit → use it; if none → infer default (aliphatic single; aromatic `:` when both atoms aromatic).
    - Directionality across ring: require consistent `/` or `\` on both sides only when part of an E/Z-capable double bond; otherwise error.
    - Acceptance: targeted unit tests for precedence and directionality pass.
  - M3-5: Structural validations
    - Self-loop detection: closing to same `atom_id` → error.
    - Two-member ring detection: closing to the immediately previous atom in the traversal → error.
    - Acceptance: unit tests covering `C1C1`, `C1CC1` negative/positive cases.
  - M3-6: Branch/group interaction
    - Ring indices can appear before, inside, or after branches/groups.
    - Ensure branch frames don’t interfere with ring table; `pos` of ring tokens is correct inside groups.
    - Acceptance: tests like `C1(C)CCC1`, `C1C(=O)NC1` pass.
  - M3-7: Unclosed rings at EOI (and component boundaries if/when supported)
    - On parse end, emit `SYN_RING_UNCLOSED` for any still-open ring entries.
    - Acceptance: tests for single unmatched ring and multiple unmatched rings (report last unmatched by spec policy) pass.
  - M3-8: Keep M3 as a superset of M2
    - Reuse M2’s atom/bond/branch error paths unchanged; add ring handling orthogonally.
    - Acceptance: run full M2 test suite against M3; zero regressions.
  - M3-9: Valid rings tests (rstest tables)
    - Aliphatic: `C1CCCCC1`, `C1CCC2CCC1CC2` (fused), `C1CCC2(CC1)CCC2` (spiro), multiple indices on one atom (`...12`).
    - Aromatic: `c1ccccc1`, `c1ccc2ccccc2c1`, and combinations with explicit `:` vs implicit aromatic bonding.
    - With bonds: `C1=CC=CC=C1`, `C1/C=C/C=C/C/1`.
    - With branches: `C1(C)CCC1`, `C1=CC(C)=CC=C1`.
    - `%DD`: `C%12CCCCC%12`.
    - Acceptance: expected IR matches builders; aromatic/default bond inference matches spec.
  - M3-10: Invalid rings tests
    - Unclosed: `C1CCC`.
    - Self-loop: `C1C1` where closure hits the same atom.
    - Two-member: adjacency closure creating a 2-cycle.
    - Bond conflict: `C=1...C#1`.
    - Direction conflict: `/` vs `\` on opposite sides or without a double-bond context.
    - Bad `%` forms: `C%1`, `C%1a`, `C%001` (if disallowed).
    - Acceptance: precise `pos` (and `open_pos` where applicable) and correct error variants.
  - M3-11: Builders for test expectations
    - Ensure `build_ring_c` and `build_two_rings_c` produce the same DFS order as the parser for aliphatic and aromatic cases, including fused and spiro junctions.
    - Acceptance: no IR ordering mismatches in snapshot tests.
  - M3-12: Benchmarks in `benches/opensmiles_parsing.rs`
    - Add groups: simple cycles (C6/C12), fused polycycles, spiro, multi-ring with bonds and branches, aromatic sets.
    - Compare M3 vs M2 on inputs without rings to catch regressions.
    - Acceptance: performance within ≤5× lex-with-bonds baseline; record timings.
  - M3-13: Spec and error registry
    - `spec/opensmiles-umol.md`: ring tokens, precedence, defaults, aromatic interplay, branch interactions, unclosed ring policy.
    - `spec/opensmiles-errors.md`: add ring error codes listed above with stable identifiers.
    - Acceptance: docs updated and consistent with tests.
  - M3-14: Linter (optional after core lands)
    - Style: conflicting dual explicit bond on same ring index (warn), normalizing ring index width policy (if desired), heuristics for redundant dual ring indices when a branch would be clearer.
    - Acceptance: compile and basic snapshots; no parser coupling.

- Out-of-scope for M3 (unless requested)
  - Extended ring index ranges beyond `%DD`.
  - SMARTS ring semantics and aromatic model changes.
  - Kekulization or aromaticity perception beyond what’s needed for bond defaulting at closure.

- Risks and mitigations
  - Aromatic defaulting: keep rules minimal—only default `:` when both atoms are aromatic and neither side specified a bond.
  - Directional bond validation: restrict to double-bond contexts to avoid false positives; add narrow tests first.
  - DFS ordering mismatches: lean on existing builders and snapshot diffs early.

- Revisions
  - Keep only %DD and lex %DDD as [%DD][D], which is consistent with the OpenSMILES spec. If we want indices > 100, we'll need to use another method. One possibility is %%DDD, %%%DDDD (n % symbols capture n + 1 following digits), another would be to introduce an optional stop symbol. For now, we should just note the current behavior and move on.
  - For the errors:
  RING_UNCLOSED, RING_SELF_LOOP, RING_TWO_MEMBER, RING_CONFLICT_DIR already exist in the RING category. Let's keep them
  Rename RING_CONFLICT_DIR  (exists) -> RING_DIR_CONFLICT, SYN_RING_BOND_CONFLICT -> RING_BOND_CONFLICT
  Rename LEX_BAD_PERCENT_FORM (exist) -> LEX_RING_INDEX_INVALID
  Add SYN_LEADING_RING -> LEX_LEADING_RING
  Leading zeros are allowed in percents, use STYLE_UNNECESSARY_PERCENT_RING_INDEX, STYLE_REUSED_RING_INDICES, STYLE_FIRST_RING_NOT_ONE, STYLE_NONCONSECUTIVE_RING_NUMBERING style lints.
  Not also that ring index 0 (also %00) is allowed.
  We'll need to review the categories and naming scheme a bit later. This should suffice for now.
  - We should try to remove hashing from the ring table, either using a Vec (probably easiest) or an HashMap with identity hash from SwissTable (hashbrown).

### M3 Tasks

1. Add module `io/smiles/fsm_m3.rs` and export the parser entrypoint.
2. Implement ring lexing: single digits 0–9 and `%DD`; tokenize `%DDD` as `%DD` + `D`.
3. Emit `LEX_RING_INDEX_INVALID` for bad `%` forms and `LEX_LEADING_RING` for ring tokens without a current atom.
4. Capture optional preceding bond/dir and precise positions for ring tokens.
5. Implement ring table as `Vec<Option<OpenRing>>` for indices 0..99 (including 0).
6. Record ring open state: `atom_id`, `bond_opt`, `dir_opt`, `pos_open`.
7. Close rings with bond precedence; default aromatic bond only if both atoms aromatic and no explicit bond.
8. Validate directional markers across ring closures; emit `RING_DIR_CONFLICT` on mismatch.
9. Detect self-loop and two-member rings; emit `RING_SELF_LOOP` and `RING_TWO_MEMBER`.
10. Emit `RING_UNCLOSED` at end-of-input for unmatched indices (using `pos_open`).
11. Ensure M3 is a strict superset of M2; run the full M2 test suite against M3.
12. Add valid ring tests: aliphatic/aromatic cycles, fused, spiro, with bonds/branches, `%00` and `%DD`.
13. Add invalid ring tests: unclosed, self-loop, two-member, bond/dir conflicts, bad `%` forms, leading ring.
14. Update specs: ring tokens, `%DDD` note, index 0 allowed, precedence/dir rules and interactions with groups/branches.
15. Add ring index edge tests (C1CC%01, C0CC%00, %99, invalid %0).
16. Map new M3 errors to diagnostics in the linter.
17. Add proptest fuzzing for random bytes/prefixes.
18. Optional micro-optimizations in ring closure.

### Fuzzing 

- 1) Timing
  - Do a minimal, crash-only fuzz pass now; defer deep fuzz until after M4/M5 stabilize. Cost now is low and can catch panics early; comprehensive fuzzing is more valuable once the grammar surface grows.

- 2) Grammar-generated fuzzers
  - Auto-generating from EBNF is nontrivial in Rust; there’s no turnkey “derive-fuzzer-from-EBNF.” You’d end up writing a generator (proptest strategies or an AST builder) from the grammar anyway. That’s effort-heavy and risks drift from the real parser semantics. I don’t recommend it now.

- 3) Without grammar: how to generate and get coverage
  - Use coverage-guided fuzz (cargo-fuzz/libFuzzer) on byte inputs to shake out panics, seeded with our existing SMILES corpus; add a small dictionary of tokens ('(', ')', '=', '#', '/', '\\', '%', digits, 'Cl', 'Br', 'c', 'n', etc.) to guide mutations.
  - Complement with a light token-stream proptest generator (handful of combinators for nodes, branches, ring closures) to exercise invariants deterministically. This avoids building a full grammar but still hits structured shapes.
  - Track coverage (llvm-cov or cargo-fuzz’s built-in stats), and iterate seeds from failing examples and tricky tests (rings, directions, percent indices, brackets) to push coverage up.

- 4) Library choice
  - Keep proptest for property tests (deterministic, CI-friendly).
  - Add cargo-fuzz (libFuzzer) for crash/UB hunting and coverage-guided exploration. No need to switch away from proptest; use both.

- 5) Effort estimate
  - Minimal setup (cargo-fuzz target that calls parse_smiles_m3 on byte slices, + token dictionary, + seed corpus from tests): ~1–2 hours.
  - Adding 2–3 simple proptest strategies for token streams and basic invariants (no panic; result or specific error): ~1–2 hours.
  - Full grammar-driven generator: 1–2 days. Not worth it before M4/M5.
  - So we can stay under half a day for a useful baseline and revisit after M4/M5.

If that sounds good, I’ll set up cargo-fuzz with:
- Target: parse_smiles_m3(bytes) with size cap, no panics allowed.
- Seeds: ring-heavy, stereo, percent, bracket samples from tests.
- Dictionary: core bond/direction/ring tokens.

### M4 Tasks

1. Create `fsm_m4.rs` with `parse_smiles_m4` as a strict superset of M3.
2. Implement '.' as component separator; reset adjacency while preserving open rings.
3. Allow ring closures across components (indices can open in one component and close in another).
4. Run full M3 test suite against `parse_smiles_m4`.
5. Add valid component tests: simple `CC.CC`, branched `C(C).C(C)`, mixed aromatic/aliphatic, multiple dots.
6. Add invalid tests: leading dot `.C`, trailing dot `C.`, consecutive dots `C..C`, dot before ring `C.%12` and `C.1` (invalid), direction/bond immediately after dot `C.-C` (invalid).
7. Add ring-across-components tests: `C1.CC1`, `C%12.CC%12`, stereo `C/1.CC/1`, mixed open/close orders where allowed by M3 rules.
8. Add ring stereo across components tests for up/down consistency.
9. Linter: map dot-related diagnostics (e.g., `LEX_DOT_BEFORE_RING`, `LEX_CONSECUTIVE_DOTS`, `LEX_LEADING_DOT`, `LEX_TRAILING_DOT`) and ensure pass-through.
10. Update `opensmiles-errors.md` with component diagnostics definitions and spans.
11. Update `opensmiles-umol.md` with M4 semantics: component separation, ring persistence across components, constraints.
12. Add Criterion benches targeting component-heavy inputs.
13. Proptest: extend ASCII alphabet to include '.', ensure no panics.
14. Fuzz: add '.' to dictionary and seed corpus.
