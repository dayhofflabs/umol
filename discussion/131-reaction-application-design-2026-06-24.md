# 131 — Reaction application: prior-art review and design

Designing the end-user API for applying `ReactionAst` rules to molecules. The DPO
primitive (`MoleculeAst::apply_rule`) and the matcher (`substructure_matches`) exist
but are disconnected and uncalled (see doc 127). This doc records a review of how
existing codes do it, the resulting design decisions, and open directions.

Reaction is intrinsically multi-molecular: molecularity and product count are both
variable, and one rule application yields several reactions (one per match).

## Prior art

Reviewed under `materials/codes/`: RDKit, LillyMol, CDK, Indigo, MØD; plus
`/Users/dr/Dayhoff/molintern` (CGR on RDKit).

### Two rule representations

| | rule is | a change is | atom map |
| --- | --- | --- | --- |
| **Template + map-number** (RDKit, CDK, Indigo, LillyMol-SMIRKS) | LHS/RHS templates joined by integer `:n` | diff LHS↔RHS per map number (compiled to op-codes, or reconstructed at apply time) | a property to track/propagate |
| **Core graph + membership** (MØD) | one graph, each node/edge tagged `{L,K,R}` | a `K` element whose left label ≠ right label | free — `K` is the shared interface; the DPO comatch yields it |

`ReactionAst` today is template-plus-map (lhs/rhs `MoleculeAst` + `Option`-pair map,
K inferred) — correct as the SMIRKS-shaped surface. MØD's `lib::DPO::CombinedRule`
(one graph + `Membership{L,K,R}`, labels stored as parallel left/right value vectors,
`isChanged` = `K` with differing labels) is the form in which the DPO machinery is
clean: changes are explicit and the atom map is free.

### One core, many drivers

Every serious code is a single rewrite core with thin drivers — the "one core,
multiple APIs" plan is the norm:

- **RDKit**: `run_Reactants` (eager, combinatorial) is the core; `EnumerateLibrary`
  is a lazy index-space layer that just calls it; `run_Reactant` is an in-place
  single-molecule fast path.
- **LillyMol**: `_perform_reaction` is the core; `trxn` / `multiple_reactions` /
  `random_reactions` add scaffold-multiplicity, chaining, sampling on top.
- **CDK**: SMIRKS string → ordered **op-code program** (`TransformOp`) → generic
  `Transform` engine (parser/compiler cleanly separate from the rewrite engine;
  per-match apply is transactional with rollback). `SmirksOption` carries dialect
  presets (`RDKIT`, `DAYLIGHT`) because toolkits disagree on implicit-H / element-set
  / bond-overwrite semantics.
- **MØD**: no separate apply engine — a direct derivation **is** a rule composition
  (educt graph → bind-rule `∅←∅→G`, composed with the reaction rule). One mechanism
  serves single application and whole-network generation.

### Atom maps: lazy, not stored

- **LillyMol** strips atom maps by default (`reset_all_atom_map_numbers`); enumeration
  is embedding/index-driven, map-driven only on request.
- **MØD** stores only reaction topology (educt/product multisets + rule pointer) and
  **recomputes** the atom map on demand (`VertexMapper`) from (educts, rule, match).

So the atom map is not a mode to toggle on the product — it is re-derivable from the
match. The generation paths never materialize it; the mapped path re-derives it.

### The three application modes in the wild

- **Single mapped reaction** (our ii): RDKit `run_Reactant` in-place; CDK `Transform`
  exclusive apply; Indigo `ReactionTransformation` (loops over non-overlapping sites
  in place, returns the final AAM).
- **Combinatorial library** (our iii-a): LillyMol — scaffold + N reagent *lists*; a
  lazy mixed-radix **odometer** (`Reaction_Iterator`) over the lists; reagents
  pre-searched once at load; one product streamed per tuple; two-tier
  (cheap-SMILES → canonical) dedup. RDKit `EnumerateLibrary` is the same shape
  (strategy over an index space, checkpointable).
- **Reaction network** (our iii-b): MØD `DG::Builder` — bind educt graphs to the
  rule's L-components one graph per round (incremental composition), enumerate educt
  *combinations* via a multi-dimensional selector, dedup products by canonical-graph
  isomorphism, grow via `repeat(applyAll)` to a fixpoint with a subset/universe split.
  This is the 100k–1B-node case.
- **Fixpoint state enumeration** (our iv): RDKit tautomer + (de)protonation
  **bypass** the reaction engine — rule = SMARTS pattern + tiny payload (bond-order
  vector / charge deltas / acid-base pair); apply by plain match + in-place edit,
  iterate to a fixpoint, dedup by canonical SMILES, then a separate scoring pass picks
  the canonical form.

So iii is two distinct layers (library odometer vs network builder) over one core; iv
is a separate lightweight path, not the DPO engine.

### Input formats

- **SMIRKS** — universal interop, but dialect-dependent (CDK encodes the differences
  as option sets). Daylight semantics are the target (well-defined for the injective
  case `ReactionAst` encodes).
- **GML** (MØD) — DPO-native: `left`/`context`/`right` is L/K/R; a change is an edge
  in both `left` and `right` with different labels. Maps 1:1 onto a membership core.
- **Native EDN** — homoiconic with `ReactionAst`; none of the codes has this; umol's
  to own.

## Decisions

1. **Core rule form — decided: `ReactionAst = MoleculeAst (lhs) + Deltas`.** One full
   *state* (the `lhs` `MoleculeAst` — all 8 entity types, changed and unchanged) plus one
   *delta* (`Deltas` = `Vec<Delta>`, the resolved transformation). The membership view
   (`{L,K,R}`), the atom map, the R-side, and the CGR are all **derived**, not stored:
   - changes are the deltas (explicit);
   - atom map = the `lhs` id space carried through the deltas (`Add` = created with a fresh
     id, `Remove` = deleted), so K = surviving `lhs` atoms;
   - R-side = `lhs` with the deltas applied; CGR = `to_condensed()` (`lhs` superimposed with
     each delta's `(left, right)` labels; molintern grounds the CGR view).

   Why this rather than a *stored* symmetric combined graph (the earlier A+B hybrid):
   - **Uniform** — all 8 entity types live in `lhs`; all changes are `Delta`s. No
     atom/bond-vs-overlay asymmetry (the hybrid had it as an artifact of the
     "A = combined *graph*" framing). The graph/off-graph split is `MoleculeAst`'s concern,
     already solved.
   - **Full molecules** — `lhs` is complete, so unchanged entities (overlays included)
     persist into R. A pure edit list could not represent unchanged overlay *state* (edits
     are changes, not state); making `lhs` a full molecule closes that gap.
   - **Homoiconic, minimally** — a molecule is the empty-edits case; a rule is a *pattern*
     `lhs` + edits; applying a rule yields a concrete reaction of the same type. Exactly
     the base+delta form. One source of truth (`lhs + edits`), nothing to keep mutually
     consistent. Tradeoff: R / CGR are *computed* (an edit replay), not stored.
   - **Max reuse** — `MoleculeAst` payloads + the `transact`/`Undo` engine, reached by
     lowering `Delta → Vec<Edit>` at apply; the shared `*FieldChange` payloads already cover
     element/charge/order changes.

   Vocabulary: the `Delta` family (Add/Remove/SetField across all 8 entity types +
   constraints), plus the one addition for stereo — fold `StereoAtom/BondCorrespondence`
   (doc 127) into the `Delta` vocabulary: preserve→no-op, create→`AddStereo*`,
   destroy→`RemoveStereo*`, and add the
   **relative** rule-level ops `PermuteStereoAtom{perm}` (proper ligand-frame relabel),
   `MirrorStereoAtom` (the improper op / enantiomer — distinct from a permutation in
   general; coincides with a transposition only for tetrahedral via `S4≅Td`), and
   `SwapStereoAtom{a,b}` (sugar for `Permute` of a transposition); same for stereo bonds.
   Relative because a rule is a pattern that does not know the matched center's coset, so
   the concrete `SetStereoAtomField::Configuration{old,new}` cannot express it. Not
   remove+add: a reconfigured center is a preserved K element (Walden inversion keeps the
   carbon); `Permute` inverts (`perm⁻¹`) and composes (`p∘q`). At apply, each relative op
   *lowers* to `SetStereoAtomField::Configuration{matched, perm·matched}`, reusing
   `transact` + `Undo`.

   One further `Edit` variant pair — `TransformFrameStereoAtom` / `TransformFrameStereoBond`
   (verb-first, per the enum's op naming): the **coset counterpart of `IdRemapping`**.
   When reindexing changes a stereocenter's atom-index frame (composition's aliasing, or
   any renumbering), it keeps the *unchanged* configuration correct by applying the induced
   ligand-frame permutation's action to the coset (via the umol-perm coset algebra for the
   class). It carries the **before/after frames** (ordered ligands incl. virtual
   H / lone-pair positions, per the raise convention) — frames, not just atom indices, so
   it also covers virtual↔explicit swaps (an H-count change at the center); it is
   self-inverting (swap before/after). Distinct from the authored `PermuteStereo*`:
   `Permute` is a rule *intending* a frame change, `TransformFrame` is mechanical
   bookkeeping for a configuration that did *not* change. Genuine ligand-set changes
   (gain/loss of stereogenicity) are not a frame permutation — those go to membership
   (`Add`/`RemoveStereo*`).

   This makes **stereo-under-composition mechanical, not novel**: because the coset is
   defined in the atom-index frame, composition's reindexing yields the before/after frame,
   and `TransformFrame` applies the coset action. The action is a group action (associative)
   and frames are derived from the (associatively-composing) graph state, so the stereo
   layer composes associatively on top of the base — associativity is inherited, not
   re-earned. The algebra's novel surface thus reduces to the topology composition (the
   Behr papers) plus this coset bookkeeping.

   Naming — **it is an AST**: frozen (combined) topology, lattice-valued entities
   (`AtomAst`/`BondAst`), homoiconic (ground ⇒ concrete reaction, pattern ⇒ rule) — the
   `MoleculeAst` family. So the core is **`ReactionAst`**. Surface/input forms compile *to*
   it, mirroring `MoleculeDsl → MoleculeAst`: SMIRKS, GML, **`ReactionDsl`** (the EDN serde
   boundary), and an `ReactionAst::from_lhs_rhs_map(...)` constructor. The *current*
   `ReactionAst` struct (lhs/rhs + map) is renamed/demoted — it is a surface form, not the
   AST. Inherits AST treatment (canonical equality up to iso; `ReactionDsl` serde); a
   rule-level `Lattice` is not assumed (each L/R projection is a `MoleculeAst`, which is a
   lattice; whole-rule meet/join is not obviously meaningful).
2. **Swap boundary** = `compiled rule × host(s) → products`, with the rule-compile
   step *before* it (CDK op-codes / MØD composition) and the combination/enumeration
   layer *above* it (RDKit index space, LillyMol odometer, MØD builder). Atom map is
   re-derived from the match, never stored on the product. APIs depend on this
   contract so the generic core (today's `apply_rule` + matcher + a
   `MoleculeEmbedding→match` bridge) can be replaced by an op-program or composition
   engine later.
3. **iii is two layers over one core**: a combinatorial-library odometer (fixed roles,
   reagent lists, streamed, pre-searched reagents) and a network builder
   (open educt combinations, canonical-form dedup, repeat-to-fixpoint).
4. **iv is a separate saturation path**, not the DPO engine: a catalog of
   (pattern + payload) applied to a fixpoint, deduped by canonical form, with a cost
   function for canonical selection.

## Open directions (to digest)

### Rule composition — and why apply is a special case of it

DPO rules compose. Given `p1: L1←K1→R1` and `p2: L2←K2→R2` plus an overlap
(a partial match of `L2` into `R1`), the composite `p2∘p1` is built by a
pushout: anything `p1` creates that `p2` consumes becomes internal (in neither
composite L nor R); anything `p2` needs that `p1` didn't supply is added to the
composite L (pulled back through K1); anything `p1` produced that `p2` leaves alone
passes to the composite R; the composite K is the shared context. This is the DPO
concurrency theorem. MØD computes it directly (`lib::RC::CompositionHelper`).

**Apply = compose with a graph-as-rule.** A molecule G is the degenerate rule
`∅←∅→G` (pure R: it exists, consumes/creates nothing). Composing a reaction rule with
G's bind-rule = matching the rule's L into G and rewriting = applying the rule to G.
That is why MØD needs no separate apply engine. Composing reaction deltas = a single
combined delta (a reaction *path* condenses to one rule / one CGR).

**Algebraic structure — the rule algebra.** Composition is not unconditional
(`A*B=C`): it is parameterized by the overlap, so for a fixed pair it yields a *sum*
over all admissible overlaps, `A*B = C1 + C2 + …`. With `+` the formal sum (free module
over rules) and `*` this sum-over-overlaps, the result is bilinear (so left/right
**distributive**), **associative**, and unital. This is the **rule algebra** of Behr,
Danos & Garnier (2016) and Behr & Sobociński (2020) — PDFs in `materials/graph_rewriting/`.
Two traps to avoid: the overlap ranges over `R_A ↔ L_B` (including created/deleted
elements — e.g. A forms a bond B breaks), **not** `K_A ∩ K_B`; and composition is the
pushout case analysis above, **not** independent set-union/difference on the three
components (those drop the shared K interface).

The exact conditions, pinned from Behr-Sobociński 2020 (the categorical, DPO-with-monos
formulation — the one to follow):
- **Overlap class (Def 2.6):** an admissible match of `p2` into `p1` is a **span of
  monos** `I₂ ←m₂— M₂₁ —m₁→ O₁` — `M₂₁` is the overlap of rule-2's input (`lhs`) with
  rule-1's output (R-side); both legs are in the mono class `𝓜`. **Injective only** — no
  non-injective overlaps (Convention 2.4; 2016 paper: an overlap is "an injective partial
  morphism from the output of B to the input of A"). The set `M_{p2}(p1)` is **finite**
  (finitary 𝓜-adhesive category: finitely many 𝓜-subobjects), enumerated up to iso.
- **Admissibility (Def 2.3/2.6):** form the pushout `N₂₁` of the span; the overlap is
  admissible **iff both pushout complements exist** — one for `p1`'s context against
  `O₁→N₂₁`, one for `p2`'s context against `I₂→N₂₁`. POC-existence *is* the DPO gluing
  condition; in an adhesive category it subsumes the classical dangling + identification
  conditions. Missing POC ⇒ dangling/identification failure ⇒ that overlap drops from the
  sum. (The trivial overlap `M₂₁ = ∅` is always admissible.)
- **DPO vs SqPO is only a post-composition policy, not a different overlap class** (2016
  paper, fixing morphism `𝓕_T`): all variants enumerate the *same* overlaps; **DPO
  discards** dangling-edge composites (project to 0), while **SPO/SqPO auto-deletes** the
  dangling edges (`SPO_A` = SqPO for injective matches + linear rules). So one overlap
  enumeration serves both — the DPO/SqPO choice is a dangling-handling *policy* applied
  after composition, the same axis as `apply`'s dangling handling.
- **Algebra (Def 3.3/3.4, Thm 2.9/3.5):** `+` = free vector space over iso-classes of
  rules; `*` is **bilinear** (hence distributive); **associative** given 𝓜-effective
  unions; **unital** with the trivial empty rule `R_∅ = (∅ ← ∅)` given an 𝓜-initial
  object. Finite graphs are adhesive with both, so we inherit associativity + unit — do
  not re-derive.

**Composition on `(lhs, edits)`, and why the layout is unchanged.** `compose(A,B)` =
for each **overlap** of A's R-side with B's `lhs` → build the composite by **aliasing**
the overlap's B-atoms onto A's ids (id-reuse — "identification" is the map, not a
node-merge), **appending** B's non-overlap context and its remapped edits. Summed lazily
over admissible overlaps.

The overlap is a *partial common subgraph* `M` of (A's R-side, B's `lhs`), ranging from
`∅` up to the maximum — **not** subiso of the whole `lhs` (that is only the `M = lhs`
special case). So composition needs a **common-subgraph *enumeration* primitive**,
distinct from the apply path's subiso (apply needs the whole `lhs` in the host;
composition needs all partial overlaps). It is the all-solutions sibling of graph-core's
existing MCS (`maximum_common_*subgraph(s)`, McGregor, already returns the node
correspondence): enumerate *every* common subgraph with its embedding span, lazily and
bounded (≤ 2^(MCS size)), instead of only the maximum — a new mode on `mcs.rs`, not new
substrate. The empty overlap (always admissible) gives the parallel/independent
composite. Pin induced-vs-non-induced to the papers' `𝓜`-mono class (graph monos are
non-induced ⇒ the edge-subgraph / MCES family). This primitive is also the substrate for
critical-pair / confluence analysis, so it earns its place in graph-core.

Consequences for the type:
- **`+` is a collection, `*` is an operation** — `+` = `Vec<ReactionAst>` (lives in the
  saturation/DG layer, decision iii), `*` = `compose: ReactionAst × ReactionAst → impl
  Iterator<ReactionAst>` (lazy). Neither is a field; `ReactionAst` stays a single rule
  (a formal sum is neither a reaction nor a rule, so it must not enter the type). Unit =
  the empty-edits `ReactionAst`.
- **One new graph op**, not a layout change: `assemble`-along-a-partial-map (a pushout;
  generalizes disjoint `assemble` from the empty map). With mono overlaps it is clean
  1:1 aliasing — no "merge two existing nodes" case ever arises. `Edit` stays purely
  molecular (no `IdentifyAtoms`).
- **Composition runs on the derived membership view, not on the stored edits.** Only the
  naive `(L,K,R)` *set formula* (`L_C=(L_A\S)∪…`) is wrong; the proper `(L,K,R)`
  composition — MØD's `CompositionHelper`, a pushout with the per-element case analysis —
  is correct and **minimal by construction** (created-then-deleted elements become
  internal and vanish). For *that one operation* the membership form is the better
  carrier: `R1` is a first-class view (composition is `R1↔L2`), the case analysis yields
  the canonical composite directly, and the span is explicit. `(lhs, edits)` edit-concat
  is *operationally* correct (form-then-break = `[AddBond, RemoveBond]` applies to
  identity) but **non-minimal** — it keeps the cancelling pair and would need a simplify
  pass. Reconciliation: the condensed view is the same derived projection we
  need for the CGR and R-side, so **compose on the derived view** (project `lhs+edits` →
  `CondensedReactionAst` → MØD-style overlap+case-analysis → project the minimal result
  back to `lhs+edits`). `(lhs, edits)` stays the *store* (it wins on apply / uniformity /
  full molecules / homoiconicity, and apply ≫ compose for network building); the condensed
  form is the *working form* for compose and CGR. If a workflow ever became
  composition-heavy (rule-algebra analysis, pathway condensation), storing the membership
  form natively would start to win — revisit then.

### Primitive operations on `ReactionAst`

Store = `lhs: MoleculeAst` + `deltas: Deltas` (`Vec<Delta>`). Apply is the common case and is native on
the store; everything else (compose, CGR, SMIRKS, `(L,K,R)`/GML) flows through one
derivation hub plus its inverse.

**The hub derivation.**
- `to_condensed(&self) -> CondensedReactionAst` — replay `edits` on `lhs` to the
  **condensed reaction AST**, the superimposed union graph carrying a paired
  `(left, right)` value per element (from each edit's `old`/`new`) over the atom map (the
  K interface). This is a generalized, attributed CGR (= MØD's `CombinedRule`): the
  `{L,K,R}` membership is *derived* from each pair (`left` absent ⇒ created, `right`
  absent ⇒ deleted, `left==right` ⇒ context/K, `left!=right` ⇒ modified), not stored
  separately. Generalized past classic atoms-plus-dynamic-bonds CGR in three ways — values
  are lattice-valued `AtomAst`/`BondAst` (a *pattern* rule's condensed form is a query
  CGR), overlays carry their own `from→to`, and atoms may be created/deleted (atom balance
  is not assumed). Linear in `edits`. Convenience reads: `right() -> MoleculeAst`,
  `membership()`, `atom_map()`. Every item below except apply reads from this.
- It is the **symmetric pivot**: `ReactionAst` is `lhs`-anchored, but the condensed form is
  side-neutral, so `reverse()` routes through it (`to_condensed()` → swap `left`/`right` →
  project back to the rhs-anchored `ReactionAst`). The round-trip
  `(lhs, edits) ↔ CondensedReactionAst ↔ (rhs, reverse_edits)` is lossless on the **net
  transformation** (map and all values survive exactly); it canonicalizes the delta, so it
  is identity only up to edit-list normalization, not byte-identity on a hand-built
  journal.

**Apply (common case) — native on the store, does *not* need `combined()`.**
- `apply(&self, host, match) -> MoleculeAst` = remap `edits` onto the match, `transact`;
  match enumeration via subiso (`lhs` → host). [decision 2]

**Compose — via the hub.**
- `compose(&self, other) -> impl Iterator<ReactionAst>` = enumerate overlaps
  (common-subgraph enumeration of `self.right()` vs `other.lhs`), MØD-style
  case-analysis composite, project the minimal result back to `(lhs, edits)`. Lazy.

**Import — the inverse of the hub (produce `edits` by diffing).**
- `from_lhs_rhs_map(lhs, rhs, atom_map) -> ReactionAst` = diff L vs R per map number into
  an edit list (the CDK opcode-construction). Underlies SMIRKS / GML / EDN import and the
  lhs/rhs+map constructor.

**Export — NOT core methods; boundary types that *read* the core derivation.** Exactly as
`MoleculeAst` does not serialize SMILES/CTfile (that is `umol-io` + the DSL boundary
types), the reaction core exposes `to_condensed()` / `right()` / membership / atom-map and
the boundary layer renders from them:
- SMIRKS and GML (MØD-style `(L,K,R)`) → `umol-io` boundary types (external formats, like
  SMILES/CTfile), reading `to_condensed()`.
- CGR: `CondensedReactionAst` *is* the in-memory CGR (generalized); CGR *file* forms
  (CGR-SMILES, graphml, etc.) are `umol-io` boundary types serializing it.
- EDN `ReactionDsl` `FromEdn`/`ToEdn` → the `umol-ast` DSL boundary, like `MoleculeDsl`.

So the **core** primitive set is just: the hub `to_condensed()` (+ reads `right()`,
`membership()`, `atom_map()`), the diff constructor (import substrate), `apply`, and
`compose`. Everything format-facing is a boundary type reading the hub. Apply dominates
and stays native on `(lhs, edits)`, so only the rarer operations pay the (cheap)
projection; and because `CondensedReactionAst` is the MØD `CombinedRule`, the `(L,K,R)`
boundary type is nearly free and shares the hub with composition.

(`to_condensed`/`CondensedReactionAst` are settled; the remaining operation names above are
*descriptive placeholders* — final names TBD.)

### CGR — it *is* `CondensedReactionAst`

A CGR is the union graph with `(from→to)` labels per element over a shared atom map. That
is exactly `CondensedReactionAst`, so the CGR is not a separate type — it is the hub object
(`to_condensed()`), and CGR *files* are `umol-io` serializations of it.
`molintern/cgr_graph.py` builds the classic form from mapped reaction SMILES (nodes over
shared map indices with `from_charge/to_charge`; edges with `from_order/to_order`, order
`0` = bond absent on that side; `reduce` → the reaction center). Bond create/break is
`0→n` / `n→0`; the part a strict shared-atom CGR drops — atom create/destroy — is exactly
what the generalized form keeps (`left`/`right` absent), matching MØD's `CombinedRule`.

Decision implication: there is one object. It is the natural single-graph form for reaction
similarity / fingerprints / ML, the materialized reaction delta, the symmetric pivot for
forward/reverse, and the working form for compose and `(L,K,R)`/GML export. Import (mapped
reaction SMILES → reaction) lands here first, then projects to `(lhs, edits)`.

### `Deltas` as a `Canonicalize` AST component (not a `Lattice`)

The edit list is a secondary AST component of `ReactionAst`, exactly parallel to the flat
molecule-level `Constraints`:

```
MoleculeAst  = entities + constraints: Constraints   // Canonicalize only
ReactionAst  = lhs: MoleculeAst + deltas: Deltas     // Canonicalize only
```

Precedent confirmed in code: the flat `Constraints` (`constraint/molecule.rs:150`)
implements **`Canonicalize` only** (conjunctive flatten + sort + dedup); it does **not**
implement `Lattice`. `Lattice` lives on the *typed per-entity* containers
(`AtomConstraints`/`BondConstraints`/…), which are genuine value lattices. `Deltas` belongs
with the flat `Constraints`: wrap `Vec<Delta>` in a `Deltas` newtype that impls
`Canonicalize`; `ReactionAst::canonicalize` composes lhs's + deltas' canonical forms; lazy
equivalence (canonicalize-on-compare) is inherited — no bespoke reaction-equality code.
The compose "non-minimal cancelling pair" residue *is* `Edits::canonicalize`, not a
separate simplify pass.

**`Lattice` does not apply.** Edit lists form a non-commutative **monoid** under
concatenation (identity = empty edits = the homoiconic "reaction is just a molecule"
case). Canonicalization is reduction to a **monoid normal form**, not a meet/join: there
is no meaningful `join` of two edit *lists*. The value lattice still lives *inside* each
edit's payloads (a pattern edit's `new` is a lattice-valued `AtomAst`/`BondAst`), via the
entity ASTs — same split as `Constraints` (no `Lattice`) vs `AtomConstraints` (`Lattice`).
(Rule **composition** is a third, separate algebra — the Behr overlap product — neither
the `Deltas` monoid op nor a lattice.)

**The trait signature forces a clean two-scope split.** `Canonicalize::canonicalize(self)
-> Result<Self, Contradiction>` (`traits.rs:185`) takes **no context**, so
`Deltas::canonicalize` is necessarily lhs-free — it reads only the carried
`old`/`new`/`ast` values (chain connectivity), never lhs:

- **`Deltas::canonicalize`** (lhs-free, monoid normal form): the reduction system below.
- **lhs-relative minimization** (drop a delta whose `new` already equals lhs's value; spot
  an `Add` duplicating an lhs element) needs lhs → lives in `to_condensed`/`compose` on
  `ReactionAst`. `to_condensed` is the lhs-relative canonicalizer; bare
  `Deltas::canonicalize` is the lhs-free one. The signature *enforces* keeping `Deltas`
  self-contained.

**The reduction system (confirmed).** Ops interact iff they share a **target**: the entity
id for `Add`/`Remove`, the `(entity id, field)` for `SetField`. Fold each entity's
subsequence (in input order) to one normal form:

| Entity history | Normal form |
| --- | --- |
| preserved, untouched | (nothing) |
| preserved, field `f` changed `V0→…→Vn` | one `SetField{f, V0→Vn}` per changed field; dropped if `Vn==V0` |
| preserved, removed (with/without prior sets) | one `Remove{id, V0}` — prior sets erased by removal |
| created (`Add A0`), then modified `→An` | one `Add{id, An}` — sets absorbed into the `Add` ast |
| created, then removed | (nothing) — internal intermediate, fully cancelled |

So created entities never emit `SetField` (absorbed into `Add`); `Remove` subsumes prior
`SetField`s; the compose "form-then-break" residue (`AddBond`+`RemoveBond` on one bond)
cancels here — that *is* the non-minimality removal, no separate pass. `V0` for a fused
`Remove`/`SetField` comes from the first carried `old` in the chain (still lhs-free).
Constraint deltas reduce as a **set-diff** (not a chain): `Add C` ↔ `Remove C` cancel,
duplicate `Add C` dedups, net disjoint `{added}`/`{removed}` (leaning on `Constraints`'
own `Canonicalize` for C ordering).

`Err(Contradiction)` on: a disconnected chain (`Set{a→b}` then `Set{c→d}`, `b≠c`);
use-after-remove; double `Add` of one id; a `SetField` on a created entity whose `old` ≠
the `Add` ast's field value; and (cross-entity, after per-entity folding) a surviving
bond/overlay/stereo referencing a net-removed or never-created entity.

**Confluence (settled).** *Termination*: every fuse/cancel strictly drops the op count;
reordering is bounded. *Confluence (Newman)*: redexes are either same-entity (a
deterministic left-fold over the input-fixed intra-entity order — unaffected by reductions
on other entities) or different-entity (disjoint targets → commute); the cross-entity
dangling check is a property of the final set, order-independent. So local confluence holds
and, with termination, yields a unique normal form. Sequence order is **not** stored
(re-derived by the canonical topo sort at lowering).

**Open question — per-entity split.** Whether to split `Deltas` into typed
`AtomDelta`/`BondDelta`/… containers per the eight entity families, mirroring the
`*Constraints` split. The analogy is *weaker* than for constraints: constraints are
independent conjuncts (split is clean, each half a `Lattice`), but deltas are
**reference-coupled** — an atom removal cascades to its bonds / overlays / stereo, and
bond/stereo creations reference atoms, so one semantic delta spans families and the
canonical ordering is a cross-container dependency DAG. Grouping by `(family, id)` *inside*
a flat `Deltas` recovers the per-entity locality without fragmenting that DAG. Undecided.

**`Edit` enum suitability (it was designed for transactional editing, not reactions).**
Three properties fight canonicalization and need adaptation for the reaction store:

1. **Positional refs `New(usize)`** (an atom-ref to "the Nth creating edit in this batch")
   — the blocker: canonicalization reorders, which breaks ordinal references. Reactions
   need stable symbolic identity (lhs ids for preserved, fresh allocated ids for created;
   the K-interface map), never positional.
2. **Bulk variants** (`AddAtoms{Vec}`, `AddBonds{Vec}`, `RemoveTopology{atoms,bonds}`) give
   non-unique representations of the same change → normalize to one edit per target.
3. **Stack-discipline molecule-constraint edits** (`PushMoleculeConstraint`/`Pop…`) are
   order-coupled → for reactions use a constraint-set diff (the flat `Constraints` already
   `Canonicalize`s).

Already well-suited: `SetField` carrying `old`/`new` *is* the condensed `(left, right)`
pair; the per-`FieldChange` typing gives the `(target, field)` fusion key; and `inverse()`
on every `*FieldChange` is the `reverse()` machinery for free.

### Two edit types: deferred (`Edit`) vs resolved (`Delta`)

The reaction store and the transaction `Edit` are split into **two types, by type** (so
ops illegal in one context are unrepresentable in the other), not unified. The transaction
`Edit` stays untouched.

**Root axis — reference binding. Everything else is induced.** The single semantic
difference is how entity references are bound; it is *not* the value/lattice axis
(`Undetermined ↔ ground`), which lives in the shared payloads and is identical in both.

- **`Edit` — deferred / host-relative.** References are unbound (`New(n)` forward
  placeholders) or host-relative (`Id` to be matched), resolved against an environment
  *at apply time*. An open, operational recipe.
- **`Delta` — intrinsic / self-contained.** Every id (including created entities)
  is bound within the type's own frame (the lhs id space). A closed, order-free value.

The induction chain (binding decides the rest, they are not independent knobs):

> intrinsic ids → identity carried by id not position → order-free → canonical unordered
> set → per-element → persistent value
> unbound/host-relative refs → `New` is positional → order-significant → sequence → bulk
> for throughput → ephemeral recipe

**Conceptual raise/lower (the words stay reserved for `IntoAst`/`FromAst`; do not reuse
them in code).** `Delta` is the abstract/denotational end (a value); `Edit` is the
concrete/operational end (a recipe). So **apply lowers**: `apply(rule, host, match)` remaps
the rule's intrinsic ids onto the match (created entities become `New`), yielding a
`Vec<Edit>` fed to the existing `transact`. Reaction application reuses the transaction
engine — no second mutation path.

**`Edit`/`Undo` split is induced by deferral; `Delta` needs no `Undo`.** `Undo`
exists only because the deferred `Edit` cannot carry the information inversion needs: an
`Add` has no id until apply (so `Undo` carries `AddedAtom{id, ast}`), and a `Remove` is
destructive *and* triggers id remapping (so `Undo` carries the captured asts + the
`IdRemapping`). `Delta`'s vocabulary is instead **closed under inversion** — for every op there is
*another* op in the same vocabulary that inverts it (not that an op is its own inverse):
`Add{id, ast}` ↔ `Remove{id, ast}` (stable id, carried payload, no renumbering because
intrinsic ids are names not slots), `SetField{old,new}` → `SetField{new,old}`, constraint
`Add C` ↔ `Remove C` — so `reverse(): Delta → Delta` with no journal type. The carried `old`/`new` and
removed `ast` that make it self-inverting *are* the condensed `(left, right)` pair, and
cascades are explicit per-element ops over the dependency DAG (no implicit engine cascade
to capture).

**Constituent ops, factored from the variant.** The shared unit is the per-family op enum
— `AtomDelta`/`BondDelta`/… one per ABDAMNSS family, each `{Add, Remove, SetField}` over stable
ids + the shared payloads. Behavior (inversion, `(left,right)` projection, same-id fusion)
lives on these. The two container *shapes* are thin assemblies over the same ops: flat
tagged union (`enum { Atom(AtomDelta), … }`, order stored then canonicalized away) vs
per-entity struct-of-vecs (`{ atoms: Vec<AtomDelta>, … }`, order not stored, derived from
references). **The flat-vs-split question reduces to: is sequence order stored or
DAG-derived?** Ordering is a *lowering* concern (sequencing for `transact`/condensed), not
part of canonical equality — the canonical form is the unordered op set — which nudges
toward the split (struct-of-vecs) container.

**Sequencing = canonical topological sort.** Emitting a valid op sequence respects a
dependency DAG (atoms < bonds/overlays < stereo; removals = the same DAG reversed; acyclic
by that stratification). A plain topo sort is non-unique; for confluence use a
**deterministic** one — Kahn's algorithm with a key-ordered ready set — making the order a
function of the op set. This is a generic graph-core algorithm (joining subiso /
common-subgraph / `canonical_form`), called on the edit dependency graph.

**Naming (settled).** Keep `Edit`/`Undo` unchanged — they are exactly the deferred-world
pair and the `Undo` asymmetry is their defining trait; renaming `Edit` to `DeferredEdit`
would mislabel an immediately-applied transaction. The resolved form is named for **what it
is — the difference between two concrete graph states — not for its mechanism**: every word
for "references resolved" is already taken (`ground` = the value lattice, `Bound` =
`std::ops::Bound`, `anchor` = `SubPattern`), so a mechanism-name forces a collision or
jargon. "A difference between concrete states" already implies resolved refs (you can only
diff identified elements), and gives the honest low-jargon contrast `Edit` (a recipe you
apply) vs `Delta` (the difference itself, self-contained, self-inverting). So:

- resolved-edit op (sum): **`Delta`**
- per-family ops: **`AtomDelta`, `BondDelta`, `DativeBondDelta`, `AromaticSystemDelta`,
  `MulticenterBondDelta`, `NoncovalentBondDelta`, `StereoAtomDelta`, `StereoBondDelta`**

`Delta` is also the word this doc already uses for this object (reaction-as-delta,
base+delta), and the type is reusable beyond reactions (the base+delta molecule encoding,
MMP transforms below), which the non-use-case-scoped name preserves (`ReactionEdit`
rejected).

### Delta encoding (base + delta)

Two distinct ideas:

- **Reaction-as-delta** — store a derivation as (educts + rule + match), recompute
  products and atom map lazily. This is MØD's hyperedge + `VertexMapper` and is
  already implied by decision 2 (lazy map). Adopt. The CGR is the materialized form
  of this delta.
- **Molecule-as-base+delta** — store a product as parent + edit-list, materialize on
  demand. Value: memory at network scale, provenance, and composition (delta sum =
  rule composition). umol already has the substrate (`Edit` / `Transaction` /
  `MoleculeEmbedding::edits()`). Prior art is rich: CGR (reaction = delta), matched
  molecular pairs (transformation = delta), persistent/functional structures
  (structural sharing), VCS pack/delta. Caveat: at scale products are deduped by
  canonical form, so a node is reached by many paths (a DAG with multiple parents) —
  the delta layer rides *on top of* canonical-keyed node identity, not instead of it,
  and materialization for the next match is a memory↔CPU trade. So: a sound later
  optimization on the existing `Edit` infra, not the primary representation.

### Equality saturation for iv

The structure of iv (tautomers, (de)protonation, constitutional-isomer / resonance
generation by local rewiring) is exactly equality saturation: an equivalence class of
a molecule under a local rule set (the e-class), rewrites applied to a fixpoint
(saturation), then a cost function to extract the canonical/best form. RDKit's
`TautomerEnumerator` is already a hand-rolled version: canonical-SMILES-keyed set
(e-class; congruence = canonical form), fixpoint, `pickCanonical` by score
(extraction).

**Literature check (e-graphs for undirected graphs).** The e-graph *data structure*
does not generalize to undirected graphs, for a structural reason. Formally an e-graph
is a deterministic finite tree automaton — it represents a regular *tree* language, and
its engine is congruence closure over operators with **ordered children**; the
exponential compaction comes from subterm sharing via that congruence. Every published
"graph rewriting with equality saturation" result targets **rooted DAG IRs** (relational
query plans, tensor/ML computation graphs, logic-synthesis netlists, RTL), i.e. terms
with sharing — not undirected cyclic graphs. Unordered neighborhoods (the essence of an
undirected graph) are the commutative/AC case e-graphs handle worst (e-matching modulo
AC is a blow-up; matching modulo AC+idempotency is NP-complete; the mitigation is
canonical normalization). For a whole undirected graph "the same" means graph
isomorphism, whose canonical form is whole-graph (nauty) and **not compositional**, so
the e-graph's incremental congruence/sharing does not transfer — you get a *set of
canonical graphs*, not a congruence-compressed e-graph. The frontier (E-graphs Modulo
Theories) extends terms with theories/semantic e-ids but stays term-centric; no e-graph
variant whose congruence is graph isomorphism was found.

**Conclusion + the real primitive.** Adopt the equality-saturation *paradigm*, not the
e-graph structure. The paradigm realized over isomorphism classes is a
**canonical-form-keyed rewrite-closure / saturation** — non-destructive, keep-all,
saturate to a fixpoint, extract by cost — which is exactly a **derivation graph** (MØD's
`DG`) and, informally, RDKit's `TautomerEnumerator` (canonical-key set = the class,
fixpoint = saturation, `pickCanonical` = extraction). This *is* a good, graph-generic
primitive worth lifting into `umol-graph-core`, the same way subiso was — but named for
what it is (rewrite closure / saturation / derivation graph), **not** an e-graph. Shape:
inputs = seed graph(s), a `step: graph → set<graph>` rewrite, the existing nauty
`canonical_labeling` for dedup, an optional cost; output = the deduped closure (canonical
graphs + rewrite relation) to a fixpoint, plus cost-based extraction. It mirrors the
subiso split: `graph-core` owns the closure/dedup/fixpoint/extraction over an abstract
step; `umol-graph` supplies the DPO `step` and the chemistry rules. Payoff: this one
primitive underlies **both** iv (single seed, local rules, extract canonical) **and**
iii-b (the reaction network / DG) — they differ only in seed set and step. The `egg`
crate itself is out (term-vs-graph mismatch).

### Sketch: the saturation primitive

Two pieces. Only the first is genuinely graph-specific; the second is generic glue
co-located with it.

**(1) Congruence key — `graph-core`.** The isomorphism-class identity of a colored
graph as a hashable/orderable canonical form. graph-core already produces the canonical
permutation (`Graph::automorphisms(node_color, Nauty).canonical_labeling()`,
vertex-colored); the primitive wraps it into a comparable key that also folds in node
and edge labels:

```rust
// umol-graph-core
impl Graph {
    pub fn canonical_form(
        &self,
        node_color: impl Fn(NodeId) -> u32,
        edge_color: impl Fn(EdgeId) -> u32,
    ) -> CanonicalForm;
}
pub struct CanonicalForm(/* opaque; Clone + Eq + Hash + Ord */);
```

This is the analog of the e-graph's congruence — but whole-graph and non-incremental
(nauty), so it yields one representative per class rather than congruence-compressed
sharing. One extension needed: `canonical_labeling` is vertex-colored today; edge colors
fold in via the subdivision encoding already used in `fingerprint/substructure.rs`.

**(2) Saturation / derivation graph — generic, co-located in `graph-core`.** A
deduplicating fixpoint over an opaque state, parameterized by a rewrite system.
graph-core never sees labels; the domain supplies `key` (via `canonical_form`) and
`derive` (via DPO + rules):

```rust
pub trait RewriteSystem {
    type State;
    type Key: Clone + Eq + Hash;
    type Rule: Clone;
    fn key(&self, state: &Self::State) -> Self::Key;
    /// Direct derivations whose educts are drawn from `pool`, with at least one
    /// educt in `frontier` (the states added last round — incremental growth,
    /// avoids re-deriving). Molecularity = educts.len(); a mono rule (iv) returns
    /// one-educt derivations.
    fn derive(&self, frontier: &[usize], pool: &[Self::State])
        -> Vec<Derivation<Self::State, Self::Rule>>;
}
pub struct Derivation<S, R> { pub educts: Vec<usize>, pub rule: R, pub products: Vec<S> }
pub enum SaturationBound { ToFixpoint, MaxStates(u32), MaxRounds(u32) }

pub fn saturate<X: RewriteSystem>(
    seeds: impl IntoIterator<Item = X::State>,
    system: &X,
    bound: SaturationBound,
) -> DerivationGraph<X::State, X::Rule>;

pub struct DerivationGraph<S, R> { /* states deduped by Key; derivation hyperedges */ }
impl<S, R> DerivationGraph<S, R> {
    pub fn states(&self) -> &[S];                          // the saturated class (enumeration)
    pub fn derivations(&self) -> impl Iterator<Item = (&[usize], &R, &[usize])>;
    pub fn extract<C: Ord>(&self, cost: impl Fn(&S) -> C) -> Option<&S>; // canonical pick
}
```

Honesty: only `canonical_form` is graph-specific. `saturate` / `DerivationGraph` are a
generic deduplicating fixpoint — the graph-ness lives entirely in the domain's `key` +
`derive`. Co-locating them in graph-core is justified (the output is a graph; they pair
with the canonical key) but they could equally be a generic util. The
combination/molecularity enumeration stays in the domain `derive` — graph-core has no
generic multi-component matcher yet; if one is lifted later, that logic can move into
core.

**Mapping**

- **iv** (tautomers / (de)protonation / isomer rewiring): `State` = labeled molecule;
  `key` = `canonical_form`; `derive` applies the local rule catalog at every match to
  each frontier state → 1-educt derivations; `bound = ToFixpoint` (or `MaxStates`).
  `states()` = the enumeration; `extract(score)` = the canonical form. Mono ⇒ a plain
  derivation graph.
- **iii-b** (reaction network): `State` = molecule; `derive` binds the rule's
  L-components to educt combinations from `pool` (≥1 from `frontier`), DPO-rewrites,
  splits products → multi-educt derivations; `bound = MaxRounds`/`MaxStates`. The result
  is the reaction hypergraph; atom maps are recomputed lazily from `(educts, rule,
  match)`, not stored (decision 2).
- `frontier` / `pool` = MØD's subset/universe split; `ToFixpoint` = its
  repeat-to-fixpoint.

e-graph correspondence: e-class ↔ `Key` (iso class), but exactly one canonical
representative per class (no congruence compression) ⇒ this is a derivation graph, not an
e-graph — the literature conclusion above.

**The key is a strength parameter — and the fixed-atom regime is term-native.** The
"e-graphs don't fit graphs" conclusion is specifically about the *isomorphism* quotient.
Fix the atom set (distinguishable nuclei — the BO/PES regime, local rewirings only:
b2f2, homo/heterolytic make/break, `X=Y-ZH → XH-Y=Z`) and that quotient disappears: a
state is a fixed-arity record over the possible edges (`n choose 2` bond-order slots) and
the atom-property vector, with **syntactic** equality. There the rewrites are genuine
term rewrite rules and equality saturation applies natively; `key` collapses to a hash of
the vector (no nauty). So `canonical_form` / `key` is pluggable with three strengths,
same engine:

| regime | "same state" | key |
| --- | --- | --- |
| fixed distinguishable nuclei (PES, local rewirings) | syntactic vector equality | trivial hash — term-native |
| fixed framework modulo its symmetry | orbit under a fixed finite group | group / slotted canonicalization |
| variable labeling (general molecules/reactions) | graph isomorphism | nauty `canonical_form` |

Caveat on payoff: a *flat* record e-graph does not compress better than a hashset of
records — e-graph sharing needs *nested* subterms. Genuine compaction in the fixed-atom
regime would require encoding the state as a nested term over a **tree decomposition** of
the framework, so unchanged regions are shared subterms across configurations (high reuse
under local rewirings). That nested encoding is the open question for whether an e-graph
ever beats the canonical-key derivation graph here.

### SMARTS/SMIRKS → `ReactionAst`: the lift boundary (unscheduled; for understanding)

SMIRKS fixes only the **mapping** (total/bijective on the shared set), so it removes the
*correspondence* ambiguity. It does nothing for implicit H, valence, or aromaticity —
those are not mapping problems, they are the **chemical model** the whole SMARTS family
shares. The split is *correspondence* (SMIRKS solves) vs *resolution* (neither solves;
the lift must).

**Principle.** The lift is umol's existing raise + resolution applied to each side, then a
map-diff → `(lhs, deltas)`. Done right, **the DPO engine stays model-blind**: it sees
explicit atoms, explicit H, and explicit aromatic-system overlays — never "implicit H" or
"aromatic". Everything model-implicit collapses into explicit deltas *before* DPO. (RDKit's
pain is that this resolution is hidden inside the reaction engine as special cases;
relocating it into an inspectable resolved form is the win.)

**Alleviated by preprocessing:**

| Issue | Action | Class |
| --- | --- | --- |
| query ops / wildcards (`*`, `~`, `[C,N]`, recursive `$(...)`) | translate to lattice-valued `AtomAst`/`BondAst` + constraints | deterministic |
| atom map → K/L/R | diff shared map numbers → `(lhs, deltas)` | deterministic (SMIRKS); needs an unmapped-atom convention (SMARTS) |
| implicit H / valence / charge | valence solver per side → explicit H + charges as deltas | convention (model + preserve-vs-recompute) |
| aromaticity | re-perceive both sides under one model → explicit aromatic-system overlays + deltas | convention (which model) |
| dangling-on-template | detect/reject malformed rules at lift time | deterministic (validation) |

**Irreducible (preprocessing can only choose a convention or flag, not recover intent):**
unmapped-atom fate (SMARTS); stereo intent (inversion/retention/racemization often not
encodable — mapping the `@`/`/\` primitives is lossy); tautomer/resonance scope (the iv
**saturation** layer, not the lift); aromaticity-model disagreement (a model makes results
deterministic but not canonical). The lift's value here is that it makes each decision
**explicit** instead of silent engine magic.

**The distinction that makes it work:** the lift resolves **implicit-but-determinate**
(an omitted H on a fully-specified atom → concrete count) while **preserving
deliberate-wildcards** as lattice `Undetermined` (`*`, aromatic-or-not), carried into
matching via `matches()`. The resolved rule is concrete where the chemist was terse and
lattice-valued where they were agnostic — which is exactly why the AST-with-lattice is the
right lift target; neither raw SMARTS nor a fully-ground graph can hold both.

## Algorithm decisions and open points (inputs to the impl plan)

Two-increment structure: **increment 1** = σ-topology (atoms + bonds + atom/bond attribute
deltas; molecules keep overlays in `lhs`, unchanged) covering `apply` + **minimal
`compose`** + the primitives; **increment 2** = overlays + stereo + overlay composition.

### Settled

- **Application condition — DPO + injective matches.** When a delta deletes an atom the
  host bonds to outside the match, the match is **rejected** (DPO gluing/dangling
  condition); matches are injective (no identification). Rationale: DPO is complete and
  correct for the SMIRKS/bijective rule space (the actual usage), and models leaving groups
  correctly via *bond breaking* (the fragment leaves with substituents attached, as a
  separate product) — well-formed rules rarely create dangling at all (break bonds; delete
  only terminal atoms). SqPO is a strict extension (its final pullback complement coincides
  with DPO's pushout complement on well-formed matches — same construction and cost — and
  only *admits more*): cascade deletion, non-injective matching, cloning. For chemistry
  those are marginal — cascade silently breaks mass balance (useful only for
  abstraction/decomposition/excision rules), non-injective formalizes a rarely-intended
  SMARTS quirk, cloning is outside the `Delta` vocabulary. SqPO is a **future opt-in** for
  excision-style rules, not a prerequisite. Note: reaction-SMARTS messiness (implicit H,
  valence, unmapped atoms, aromaticity) is *upstream* in the SMARTS→rule lift (umol-io
  boundary), orthogonal to the pushout policy — neither DPO nor SqPO addresses it.
- **Canonical labeling (#6) — settled for σ-topology.** `canonical_form` = **nauty** (a
  dependency) on the **incidence graph already built for subiso** (the edge-color→
  vertex-color reduction: atom-nodes and bond-nodes in separate color classes, sub-colored
  by label; nauty's canonical labels pull back to the molecular graph). The only new
  canonicalization work is encoding overlays / multicenter / dative / stereo (hyperedges /
  relations) into that incidence graph faithfully — **increment 2**. WL alone is
  insufficient (incompleteness), hence nauty.
- **Delta canonicalization scope.** Tie-break for independent deltas is `lhs` ids (a fixed
  frame), so the delta normal form is **self-contained per frame** and does *not* depend on
  `canonical_form` (#6). `canonical_form` is needed for reaction equality *up to
  renumbering*, product dedup (iv), and composition output — not for bare apply + delta
  canonicalization.
- **#2 Delta reduction system — settled.** The per-entity normal form, contradiction
  conditions, and the confluence + termination argument are fixed; see "`Deltas` as a
  `Canonicalize` AST component" above for the full table.
- **#3 Created-entity numbering — subsumed by #6.** Cross-reaction equality is inherently
  up-to-renumbering, so created (R\K) entities are numbered by the condensed graph's
  canonical labeling (#6), not by a separate scheme. The per-frame delta normal form (#2)
  needs no canonical created-ids; #6 handles equality/dedup across frames.

### Open — gate the increment-1 plan

- **#4 Overlap enumeration** (minimal compose) — extend MCS (McGregor) to **all common
  subgraphs** (which overlaps; disconnected allowed; size range) + the **admissibility
  filter** (pushout-complement existence under DPO). (Behr–Sobociński.)
- **#5 σ-topology composite construction** — per-element case analysis over the condensed
  view → minimal `(lhs, deltas)`. "Implement the papers," construction to be written out.

### Open — increment 2 (design-then-plan, not blocking the first plan)

- **#7 Overlay composition** — Add/Remove/Modify overlays + attribute changes under an
  overlap; associativity (adhesive / attributed-DPO / SqPO / PBPO+).
- **#8 Stereo `TransformFrame`** — induced frame permutation from a reindexing + virtual-
  ligand (H / lone-pair) frames; coset action via umol-perm.
- **#9 Saturation for iv** — congruence key (= #6), fixpoint detection, cost/extraction
  function + rule set (domain choices).

## References

Code: `materials/codes/{rdkit,LillyMol,cdk,Indigo,mod}`; `molintern/cgr_graph.py`.
Prior docs: 050/051 (reaction research), 090 (DPO over the relational model), 127
(interim `ReactionAst`).

Graph-rewriting theory (`materials/graph_rewriting/`):

Rule algebra (composition `*`/`+`, the basis for §"Rule composition"):
- `Behr-Danos-Garnier-2016` — the original DPO rule algebra (relational formulation).
- `Behr-Sobocinski-2020` — rule algebra for `𝓜`-adhesive categories; the categorical
  DPO-with-monos formulation we follow (admissible-overlap class, associativity, unit).

Relations / overlays + (changing, lattice-valued) attributes — overlays are n-ary
relations carrying electron-count/charge/spin attributes that *change*:
- `Lack and Sobociński 2004 - Adhesive Categories` — DPO (and hence the rule algebra)
  generalizes to any adhesive category; typed/hyper/attributed graphs are adhesive, so
  overlays-as-relations are not special. The foundation under the rule-algebra papers.
- `Ehrig et al. - Fundamental theory for typed attributed graph transformation` — how
  attributes live in the adhesive structure so a change is a K element with `left≠right`
  labels (attributed DPO — umol's chosen form).
- `Orejas 2011 - Symbolic graphs for attributed graph constraints` and `Orejas and
  Lambers 2010 - Delaying constraint solving in symbolic graph transformation` —
  attributes as variables + constraints (symbolic), matching umol's *pattern/lattice-valued*
  attributes and constraint system; "delay solving" = match a rule without grounding
  attributes.

Relabeling-friendly alternatives (weigh against umol's DPO-with-explicit-removal choice,
doc 090):
- `Overbeek et al. 2021 - Graph rewriting and relabeling with PBPO+` — a **lattice on the
  label set** + order-preserving morphisms (= the `ValueAst` lattice directly); cleanest
  relabeling model, but no rule algebra yet (composition would be a gap).
- `Corradini et al. 2006 - Sesqui-Pushout Rewriting` + `Behr 2019 - Sesqui-pushout
  rewriting` — SqPO = "deletion in unknown context" = automatic overlay cleanup; Behr 2019
  gives its rule algebra. Take this route if auto-cleanup is wanted.
- `Behr and Krivine 2021 - Compositionality of rewriting rules with conditions` — rule
  algebra with application conditions (umol's constraints surviving composition).

Decision: `#1–#3` (adhesive + attributed/symbolic) are the necessary core for overlays —
umol's attributed DPO is covered by them plus the adhesive rule algebra already in hand;
SqPO / PBPO+ are alternatives, not prerequisites.
