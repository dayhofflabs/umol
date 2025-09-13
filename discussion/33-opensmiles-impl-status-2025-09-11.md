### 1) Parsers status and complexity

- Current complexity
  - The tree grammar handles components at top level and inside branches, with finalization centralized in `ParseState`. Most semantic weight (component boundaries, bookkeeping) is in state rather than productions, which keeps the grammar compact and LR-friendly.
  - With rings and stereo added, the grammar surface should remain small if we continue the same pattern: recognize tokens/structures in productions; do all pairing, checks, and IR writes in `ParseState`. This keeps LR conflicts low and productions readable.

- After rings and stereo
  - Rings: treat ring digits as tail forms that annotate the current atom (optionally with a preceding bond). All pairing and bond construction happens in `ParseState`. Grammar stays simple; state gets a “pending ring endpoints” map and validation rules.
  - Stereo: grammar only carries `/`, `\`, `@`, `@@`, bracket chiral classes through to `ParseState`. Computation/validation (e.g., double-bond E/Z resolution, atom-center parity) is post-parse or late-parse, not in productions.
  - Net: complexity remains tolerable if we avoid mixing ring/stereo branching into many productions. The coupling is semantic, not syntactic, so LR remains a good fit.

### 2) Acceptor grammars and feature plan

- Do we need to backport disconnected components to acceptors?
  - Not required. Acceptor grammars are most useful as minimal recognizers for targeted micro-tests and diagnostics, not feature mirrors. Value to backport only if:
    - You want micro-benchmarks that include dot-splitting at acceptor level.
    - You want a pedagogical “acceptance-only” spec mirror.
  - Otherwise, keep acceptors lean (they already prove the structural backbone). I’d only extend them for bracket atoms (to exercise the bracket micro-grammar) and maybe ring tokens for sanity checks. No strong value in fully mirroring components or stereo.

- Sequencing of features
  - Recommended order (optimize for implementation clarity and stable performance):
    1) Bracket atoms (non-organic elements, isotope, charge, H-count, atom class; capture chiral flags but defer final interpretation)
       - Independent of rings/aromaticity; ubiquitous. Good ROI early.
       - Parse bracket payload with a small, tight sub-parser over bytes (no allocation), then call `ParseState.add_atom_with_props(...)`. This aligns with your zero-copy preference [[memory:7124069]].
    2) Bond stereo markers `/` and `\` on double bonds (acyclic first)
       - Only record markers and adjacency; leave E/Z resolution for after the connectivity is fully known.
    3) Ring closures
       - Introduce ring digit handling and `%nnn`; maintain pending endpoints in `ParseState`. Enforce component boundaries (no cross-component links).
       - On closure, construct the bond with the correct order/aromatic flag and propagate slash markers if present.
    4) Stereo resolution
       - Now compute E/Z from slash patterns (including across rings) and interpret atom-centered stereo (`@`, `@@`, chiral classes) once neighbor sets are final.
       - Keep the interpretation and validation as a late pass inside `ParseState` (or immediately after parse completion).
    5) Aromaticity
       - Start by honoring input aromatic symbols (lowercase atoms) and aromatic bond symbols. Defer aromaticity perception.
       - Add perception later as a separate optional pass (kekulization or electron-count algorithms), keeping it strictly post-parse.

  - Why this order
    - Brackets unblock general chemistry first and are orthogonal.
    - Early partial stereo (markers captured, not resolved) is easy to add and test on acyclic chains.
    - Rings must precede full stereo resolution and aromaticity perception.
    - Aromaticity depends on chosen ring model and perceived/declared bond orders.

  - Alternative stepping
    - If you prefer minimal coupling during development, you can do i → iii → ii → iv. That slightly delays testing of stereo markers, but simplifies mental model: add edges (rings) before interpreting stereo. Either sequence is fine; I prefer i → ii (record-only) → iii → ii (resolve) → iv for earlier feedback on stereo token handling.

- Tooling choice for these steps
  - Keep LALRPOP + logos for the top-level grammar; push new semantics into `ParseState`. This fits your zero-copy, no-allocation parsing approach [[memory:7124069]] and preserves deterministic performance.
  - Use tiny helper parsers for bracket payloads (inside actions). Combinators (nom) or small handwritten routines are both viable for this micro-grammar; choose whichever is more direct to read. For the top-level, LR remains a better fit than PEG/combinators.
  - No need to change parser technology now; the current setup will scale through rings and stereo.

### 3) Background materials and how to leverage them

- RDKit/CDK/OpenBabel sources
  - Use them to confirm behavior for edge cases (e.g., ring closure across branches/components, interaction of bond markers and ring digits, bracket defaults like implicit hydrogens and charge).
  - Extract test motifs: small strings that probe single features (multi-digit rings `%nnn`, nested branches with dots, combined slash markers around double bonds, rare bracket chiral classes).
  - Compare error policies (dangling ring, duplicate ring index, invalid isotope ranges, inconsistent stereo).

- Daylight SMILES and OpenSMILES
  - Cross-check bracket atom rules (field order, defaults, allowed ranges).
  - Verify exact lexical forms for ring digits and `%nnn`, and the semantics for bond markers `/`, `\`, and `=` around closures.
  - Clarify what must be rejected vs normalized.

- Design notes worth capturing before coding
  - `ParseState` ring table: map ring index → {atom_id, optional bond order/aromatic flag, any pending slash info, location for diagnostics}.
  - Bracket parsing contract: inputs and outputs for a single helper function that returns a fully populated atom record (element, isotope, charge, H count, chiral flag, atom class).
  - Stereo capture vs resolution: explicit decision that capture is inline during parse; resolution is late (after graph is complete).
  - Aromaticity policy: “accept as-is first, perceive later” to keep the parser deterministic and fast.

If you want, I can draft the minimal additions needed for bracket atoms next (lexer tokens, a single `BracketAtom` nonterminal, and `ParseState` entry points), along with 3-5 focused tests per feature to keep the steps tight.

## 4) Update 2025-09-12

Here’s the OpenSMILES feature matrix and current status (source: OpenSMILES spec) [opensmiles.org](http://opensmiles.org/opensmiles.html).

Atoms and symbols
- Aliphatic organic subset (B,C,N,O,P,S,F,Cl,Br,I): Implemented
- Aromatic organic subset (b,c,n,o,p,s; se, as): Implemented (input-honoring only)
- Full bracket atoms: Implemented
  - Isotope [Number]: Implemented
  - Chirality (@, @@, @THn, @ALn, @SPn, @TBn, @OHn): Implemented (see Stereo below)
  - H-count (H, Hn): Implemented
  - Charge (+, -, ++, --, +n, -n): Implemented
  - Class (:Number): Implemented (parsed only)
  - Unknown atom (*): Implemented
- Radicals via lowercase non-aromatic (spec section 6.4): Not implemented

Bonds
- Orders: -, =, #, $: Implemented
- Aromatic bond “:”: Implemented; also default to aromatic when both endpoints are aromatic and bond is implicit
- Directional markers “/”, “\” (for E/Z): Implemented (resolution pass)
- Implicit single adjacency between atoms: Implemented

Branching and components
- Branching with parentheses (optional leading bond): Implemented
- Disconnected components with “.” (top-level and inside branches): Implemented

Rings
- Single-digit ring indices (0–9): Implemented
- Two-digit “%nn” (10–99 only): Implemented (lexer enforces two digits)
- Ring-spec with optional prior bond symbol (order/dir/“:”): Implemented
- Same-component rule (no ring closure across ‘.’): Needs explicit enforcement (ring table is not cleared on dot)

Aromaticity
- Honor input aromatic atoms (lowercase) and explicit “:” bonds: Implemented
- Perception/kekulization and aromatic valence checks: Not implemented (deferred to IR-level semantics)

Stereochemistry
- Double-bond E/Z from slash markers: Implemented (late-pass over IR bonds)
- Atom-based double-bond stereo using @/@@ on vinylic atoms (shown in spec examples): Not implemented (we use slash-based E/Z)
- Tetrahedral (@, @@, @THn): Implemented
  - @TH1/@TH2 → alias to @/@@: Implemented
  - Validation: requires exactly 4 substituents (explicit + bracket H, at most one implicit H): Implemented
  - TH arrangement index (n) semantics beyond viability: Not interpreted (kept as-is)
- Allene/axial (@ALn): Implemented
  - Alias @ → @AL1 and @@ → @AL2 when center has exactly two incident double bonds: Implemented
  - Viability check (axis present; each terminal has substituent or bracket H): Implemented
  - Longer odd-length cumulenes beyond simple allene: Not implemented
  - AL arrangement index semantics beyond viability: Not interpreted
- Square planar / trigonal bipyramidal / octahedral (@SPn/@TBn/@OHn): Implemented (syntax + viability only)
  - Validation: exact neighbor counts (4/5/6 respectively), no implicit H allowed: Implemented
  - Arrangement index semantics (positional patterns): Not interpreted

Syntax and error policy
- Invalid bracket forms, trailing bond, consecutive bond symbols, unclosed/unknown ring index, invalid %0/%09: Implemented (tests in place)
- Ring closure across components: Needs explicit invalidation on dot (not tested/blocked yet)
- Canonicalization/normalization: Not implemented (out of scope for parser)

Out of scope per spec
- Reaction SMILES: Not part of the spec (N/A)
- Relative/unknown stereo (@?): Not in spec (N/A)
- Twisted SMILES (conformational extension): Not implemented (extension)

If you want, I can:
- Enforce “no ring across ‘.’” by clearing/validating ring state on component boundaries.
- Add optional atom-based double-bond @/@@ handling (vinylic) to match the spec examples, alongside the existing slash-based E/Z.
- Track arrangement indexes for SP/TB/OH/AL (store and optionally validate positional patterns later).