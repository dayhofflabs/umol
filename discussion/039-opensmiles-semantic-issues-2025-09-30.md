# Semantic issues for OpenSMILES implementation
## Post-parse semantic pass requires several model definitions
- Valence model
- Aromaticity model
- Chirality model?

## Valence
  Calculation of implicit Hs. Note:
  For example, for v = 4 carbon, C = [CH4], CC = [CH3][CH3], CCC = C[CH2]C, CC(C)C = C[CH](C)C, CC(C)(C)C = C[CH0](C)(C)C
  And similarly for other elements, charge states, etc. Problematic for S, N, P, As, Se, which allow for multiple stable valence states.
  OpenSMILES seems to indicate that one has to "round up" to the nearest higher valence for implicit Hs.

## Aromaticity 
### Goals and stance
- Faithful where the spec is precise; strict and deterministic where it isn’t.
- Separate parsing (bytes → graph + annotations) from validation (semantics, aromaticity).
- Make conflicts diagnosable (1:1 rules), not silently “interpreted.”

### Inputs/annotations we accept from SMILES
- Aromatic atom request: lowercase organic subset symbols (`b c n o p s`) and bracketed aromatic `se`, `as` (lowercase) per spec’s “allowed aromatic elements” set.
- Aromatic bond request: `:` between atoms.
- Explicit non‑aromatic bond orders: `- = # $` (including between aromatic atoms).
- All other fields (H, charge, isotope, class) from bracket atoms.

We record them as annotations; no aromaticity is “decided” during parsing.

### Formal model (post-parse validation domain)
- Graph G = (V, E) with:
  - Node attributes: element, charge, implicit/explicit H, bracket flags, aromatic_atom_req ∈ {true,false}.
  - Edge attributes: bond_order ∈ {single,double,triple,quad,unspecified}, aromatic_bond_req ∈ {true,false}.
- Cycles: computed over G ignoring bond order (simple undirected cycles).
  - For “in a ring” checks, presence in any simple cycle is sufficient.
  - For more structured checks, use a minimum cycle basis (not SSSR) to avoid ambiguity in fused systems.

### Strict interpretation of the three “aromatic” notions
1) Aromatic atoms (lowercase)
- Allowed elements: {B,C,N,O,P,S,As,Se} with per-element π‑electron contribution tables and allowed valence patterns.
- Structural constraints:
  - Must belong to at least one cycle. If not: error AROM_ATOM_NOT_IN_RING.
  - Must be able to support an sp2‑like configuration under valence rules given charge/H. If impossible: error AROM_ATOM_SP2_INFEASIBLE.
  - Lowercase outside those elements: error AROM_ATOM_ELEMENT_INVALID.

2) Aromatic bonds (`:`)
- Each `:` edge must lie on at least one cycle. If not: error AROM_BOND_NOT_IN_RING.
- A `:` between atoms that cannot be aromatic‑eligible (by element/valence) is an error AROM_BOND_ATOM_INELIGIBLE.
- An explicit non‑aromatic symbol (`- = # $`) is authoritative on that edge, even if both endpoints are aromatic atoms. Bonds can mix around a ring, but then see 3) below.

3) Aromatic rings (consistency of assignments)
- We do not try to “recognize chemistry by eye.” We require existence of at least one valence‑consistent Kekulé assignment for every aromatic subgraph consistent with the user’s requests.
- Formal constraint system:
  - Consider the subgraph H induced by edges that are either `:` or “unspecified between two aromatic‑requested atoms.” Edges explicitly `- = # $` are fixed.
  - Unknown edges in H can be assigned single/double subject to:
    - Atom valence capacity and charge/H constraints.
    - All aromatic_bond_req edges must be assigned “aromatic” (implemented as alternating participation in π in any valid Kekulé assignment; programmatically: allow an “aromatic” tag if the edge can be part of at least one alternating cycle).
  - If no assignment exists: error AROM_KEKULE_INCONSISTENT.
  - Optional (strict mode only): For isolated cycles, enforce Hückel 4n+2 on π count using a per‑atom π contribution table; if violated: AROM_HUCKEL_FAIL (error) or AROM_HUCKEL_WARN (warning), depending on mode.
  - For fused/polycyclic systems: do not force Hückel per‑cycle; require only a valid Kekulé assignment. (This avoids over‑constraining fused PAHs where “per‑cycle 4n+2” is not a valid criterion.)

Notes
- This makes lowercase a request (“try to make me aromatic under constraints”). If it cannot be satisfied, we produce a deterministic error rather than silently “kekulizing” in a surprising way.
- `:` is also a request on that edge. It must be satisfiable inside some alternating cycle.

### Conflict resolution (deterministic)
- Lowercase atom not on any cycle → AROM_ATOM_NOT_IN_RING (error).
- `:` edge not on any cycle → AROM_BOND_NOT_IN_RING (error).
- Lowercase atom whose valence/H/charge cannot permit sp2‑like pattern → AROM_ATOM_SP2_INFEASIBLE (error).
- Explicit `- = # $` contradicting a `:` request on the same edge → AROM_BOND_CONFLICT_EXPLICIT (error).
- Mixed explicit double/triple patterns around a requested aromatic cycle that preclude any Kekulé → AROM_KEKULE_INCONSISTENT (error).
- Optional strictness: isolated lowercase ring failing 4n+2 → AROM_HUCKEL_FAIL (error/warn per mode).

### Verification algorithm (post-parse)
1) Collect annotations on nodes/edges.
2) Find cycles; also build a minimum cycle basis for structural queries.
3) Validate local constraints quickly:
   - Aromatic atom element whitelist; degree/valence feasibility given H and charge.
   - `:` edges in cycles.
4) Build a constraint problem for the “aromatic subgraph” H:
   - Variables: for each undecided edge e in H, x_e ∈ {single,double}.
   - Constraints: per-atom valence count; fix explicit edges; ensure `:` edges can be part of an alternating cycle (implementation: tag edges on at least one alternating cycle of H; a practical approach is to run an alternating‑parity DFS/Kekulé feasibility check rather than a general SAT).
   - Solve via greedy Kekulé feasibility (standard algorithms used by RDKit/Open Babel), or formulate as bipartite matching on the line graph for each connected aromatic component.
5) Optional strict Hückel on isolated cycles:
   - If a component is a single simple cycle with no fusions and no external π contribution, count π per per‑element rules; enforce 4n+2.
6) Produce diagnostics with spans based on:
   - Parser side‑channel for ring indices and atom positions (already gated by a flag).
   - Otherwise, edge/node byte ranges.

### Lint policy (style) with PREFER/AVOID naming
- STYLE_PREFER_AROMATIC_FORM when a fully aromatic cycle is present but user provided a valid Kekulé input (suggest lowercase form).
- STYLE_AVOID_EXPLICIT_AROMATIC_BOND (prefer implicit between aromatic atoms).
- STYLE_PREFER_BOND_SYMBOL_AT_RING_OPEN (your new rule for where to place explicit bond symbol on a ring).
- STYLE_AVOID_AROMATIC_OUTSIDE_RINGS (lowercase on chains) – in strict mode this could be an error instead.

### Modes and flags
- STRICT_OPENSMILES (default): all errors above enforced; Hückel on isolated rings can be warn or error, configurable.
- LENIENT/COMPAT flags:
  - Allow lowercase not in rings to be auto‑kekulized (warn), if you decide to support legacy behavior behind a flag.
  - Disable Hückel checks entirely if interop needs it.

### Element π‑electron contribution (needed for optional Hückel)
- A small, explicit table for {B,C,N,O,P,S,As,Se} with charge/H‑dependent contributions (e.g., `n` in pyrrole contributes 2, `n` in pyridine contributes 1; bracketed charged variants adjust). Keep this compact and documented; deviations produce AROM_ATOM_SP2_INFEASIBLE or Hückel violations.

### Implementation notes
- Keep the parser fast and ignorant of aromatic semantics; it only annotates.
- The validator runs once post‑parse on the built graph; complexity is linear to near‑linear per component with the standard Kekulé feasibility algorithms.
- Fused systems: rely on Kekulé feasibility only; do not attempt multi‑ring Hückel counting.

This yields a coherent, testable definition:
- Lowercase/`:` are requests.
- Rings are required for aromaticity.
- A feasible Kekulé consistent with requests is required.
- Optional Hückel on the simple, isolated‑ring case makes the spec’s “chemistry” claims enforceable where they are actually meaningful.

## Stereochemistry resolution

- Do the stereochemistry resolution post-parse. It’s a local graph check around each non-aromatic double bond; no need to entangle the parser.
- Effort: small. Single pass O(E). Collect marked substituents per end, validate, and deduce.

Rule (post-parse):
- For each X=Y:
  - Gather all single bonds incident to X that carry / or \ (marked substituents on X’s side). Same for Y.
  - If either side has 0 marked substituents → STEREO_DOUBLE_INSUFFICIENT.
  - If a side has >1 and they don’t all share the same orientation (/ vs \) → STEREO_DOUBLE_CONFLICT.
  - Otherwise both sides have a single effective orientation sX, sY. The configuration (cis/trans) is implied by sX vs sY; both equal vs both different are both valid. No error.

This matches:
- C/C(/F)=C(\F)/C → conflict on the left end (/ vs \).
- All your “equivalent” examples are valid (each end has a consistent effective orientation).