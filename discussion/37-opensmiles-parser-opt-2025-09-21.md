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

### M0 development plan (no Logos)

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

#### Tasks (M0)

1) Create `io/smiles/fsm_m0.rs` with byte FSM and `parse_smiles_m0(&[u8])`
2) Wire `MoleculeBuilder` calls: `on_atom`, `on_bond`, `finish`
3) Export module in `io/smiles.rs` and re-export the function
4) Add Criterion group `parse_m0_chain` to `benches/opensmiles_parsing.rs`
5) Add tests under `tests/opensmiles_parsing/` for valid chains and unsupported inputs
6) Benchmark vs `lex_only` and `parse_minimal`; record results in this doc
7) Decide go/no-go for M1 based on hitting ≤5× lex on chains



