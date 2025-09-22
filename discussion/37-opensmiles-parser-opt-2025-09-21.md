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