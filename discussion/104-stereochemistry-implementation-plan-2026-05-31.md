# 104 — Stereochemistry implementation plan

Status: **Active / implementation plan.** Design record:
[103-stereochemistry-overlay-and-port-trajectory-2026-05-28.md](103-stereochemistry-overlay-and-port-trajectory-2026-05-28.md).
Step 1 (graph-core relation infrastructure) lands before the stereo phases; not yet implemented.

The full staged plan for the stereo deliverable. Design and decisions live in doc 103 (molecule-
level semantics, the `:stereo-atoms`/`:stereo-bonds` DSL, `#T`/`#C`, config algebra, contract D1/D2/D3, port trajectory,
superseded directions); this doc is the plan to build it, grounded against the real `relation.rs`
/ `graph.rs` / `remap.rs` / `builder.rs`.

## Scope

Organic central chirality (tetrahedral) + cis/trans stereobonds; parse SMILES/SMIRKS + MOL stereo,
assert/match, convert to 3D; **CIP not needed**; non-tetrahedral *perception*, ports, and reaction
stereo-transfer deferred (the coset algebra for all geometries is built in A′). The deliverable is **Step 1 +
Phases A, A′, B–E**; Phase F (3D) and Phase G follow.

## Layering (load-bearing)

- **Structural → relation (graph-core); predicate → AST.** A stereo element is a *birelation*
  `[site]-[ligands]`; the payloads `StereoAtomAst`/`StereoBondAst` hold `kind` + `configuration` + per-site `constraints`, no structural
  refs — mirroring `AromaticSystemAst`, which holds no atoms.
- Dependency is fixed: `graph-core (NodeId, EdgeId, RelationId)` → `chemistry (AtomId, BondId,
  Ligand, …)`. Nothing in graph-core names a chemistry type; the chemistry layer plugs concrete
  participant types into graph-core generics.
- **Coset algebra is pure and separate (`umol-perm`).** The Sₙ/R permutation-coset machinery (Phase A′) is a
  standalone crate — no chemistry, no geometry (*not* on `umol-msym`/libmsym, which is geometric point
  symmetry), no AST. umol-graph depends on it; the config's dense index is the OpenSMILES arrangement number it
  computes.

## Step 1 — graph-core relation infrastructure (prerequisite)

The generic relation/birelation family that stereo storage builds on. Load-bearing graph-core
change; lands before any stereo code.

### Why the full family (no simplification survives)

- **Keep set-like (sorted) relations** — a set stored unsorted costs an O(n log n) sort per
  comparison; the extra type is one-time. Aromatic systems / multicenter bonds are unordered atom
  sets, sorted on construction.
- **Cannot make everything a birelation** — aromatic systems / multicenter bonds are *one* atom
  set; a two-factor shape is meaningless for them.
- **Cannot drop birelations** — stereo relates two *different* participant types (site
  `NodeId`/`EdgeId` vs ligand `Ligand`/future `PortId`).

### Two orthogonal axes

- **Arity** — structural, different storage: `Fixed<N>` (a `[P; N]` array) vs `Var` (offset table).
- **Factor ordering** — behavioral, identical storage: `Unordered` (canonical = sorted) vs
  `Ordered` (input order is the datum). The sole difference is canonicalize-on-construction-and-
  remap: sort vs no-op. A type-parameter marker, not a struct axis.

The within-relation ordering is what stereo needs; the *outer* ordering (relations among
themselves — incidence CSR, `RelationId`) is orthogonal and unchanged. A birelation has two
**factors**; the elements within a factor are **participants**.

### The 5 structs

| kind | struct | storage |
| --- | --- | --- |
| uni | `FixedRelationSet<P, O, R, const N>` | `Vec<[P; N]>` |
| uni | `VarRelationSet<P, O, R>` | offsets + `Vec<P>` |
| bi | `FixedFixedBirelationSet<L1, O1, const N1, L2, O2, const N2, R>` | two arrays |
| bi | `FixedVarBirelationSet<L1, O1, const N1, L2, O2, R>` | array + offsets |
| bi | `VarVarBirelationSet<L1, O1, L2, O2, R>` | two offset tables |

Existing sets are `Unordered` instantiations (`VarRelationSet<NodeId, Unordered, AromaticSystemAst>`).
The stereo-atom table is `FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, StereoAtomAst>`
behind a `StereoAtom` alias (`Ordered`/`Unordered` coincide at N=1). Named `type` aliases keep call sites
readable.

### `FactorOrdering`

```rust
pub trait FactorOrdering { fn canonicalize<P: Ord>(participants: &mut [P]); }
pub struct Unordered;   // canonicalize = sort_unstable
pub struct Ordered;     // canonicalize = no-op  (monomorphizes away; zero cost in release)
```

Replaces the unconditional `p.sort_unstable()` (`relation.rs:64`) with `O::canonicalize`, in `new()`
and after a remap relabels participants (`Unordered` re-sorts; `Ordered` keeps positions so a stereo
config / bond direction stays valid).

### `RelationParticipant`

```rust
pub trait RelationParticipant: Copy + Ord + Hash {
    fn remap(self, r: &Remapping) -> Option<Self>;   // forward; None ⇒ removed
    fn unmap(self, r: &Remapping) -> Self;           // inverse; total (surviving ids)
    fn refs(self) -> ParticipantRefs;                // contained NodeId/EdgeId, for the incidence index
}

#[derive(Clone, Copy)]
pub struct ParticipantRefs { pub node: Option<NodeId>, pub edge: Option<EdgeId> }
// returned by value, not a callback: ≤1 ref per id-space (PortId spans node+edge), nothing to stream
```

graph-core impls it for `NodeId` (`map_node`/`unmap_node`) and `EdgeId` (`map_edge`/`unmap_edge`); the
chemistry layer for `Ligand` (and later `PortId`). Not `RelationId` — no relation's participants are
relations today, and a relation-ref could not remap under this `Remapping` anyway (see Forward-looking).
Routing is static per type — no runtime tag,
no dispatch; trait, not enum, because each factor is homogeneous at instantiation.

- **Removal is uniform drop-the-relation** — already the current behavior: `apply_to_var_relation_set`
  (`graph.rs:387`) collects `Option<Vec<_>>` and the `filter_map` drops the whole relation if any
  participant remaps to `None`. The trait preserves it; nothing extra is needed for removal.
- **`refs` feeds the incidence index** (design below).

### Incidence

Today one `NodeId → RelationId` inverted index (`relation.rs:19`); generalized to **one dumb union
index per id-space** — a node-index and an edge-index — built by routing each participant's `refs()`
to its space, both keyed to the same `RelationId`. `incident(node)` / `incident_edge(edge)` slice
their `(keys, rels)` pair.

- **No per-factor split.** The index is factor-agnostic; accessors filter the incident slice by role
  where needed (the existing `connecting_id` idiom, `aromatic_system.rs:87`). Per-factor indices give
  zero mapping benefit (incidence is rebuilt by `new()` on every edit, never touched by remap/unmap),
  partition rather than reduce entries, duplicate keys across indices, and force multi-index union
  queries — not worth saving an `== site` filter over a slice bounded by ≤1 aromatic system / a few
  stereo elements.
- **`bond.stereo()` needs no filter** — edges appear only as sites (ligands are nodes), so
  `incident_edge(bond)` is exactly the stereobond lookup. Node queries wanting site-only filter, since
  a node is a site in one element and a ligand in others.
- **Dedup `(key, rid)` per relation** in `build_incidence`: an `ImplicitHydrogen`/`LonePair` ligand
  holds the *site* node, so without dedup the site is emitted once per virtual ligand plus once for
  the site factor.

### graph-core changes

1. `Remapping`: rename the forward `node`/`edge` → `map_node`/`map_edge` (symmetry with the inverse;
   ripples into `IdRemapping::atom`/`bond` and `apply_to_*`), and add `unmap_node`/`unmap_edge` by lifting
   the inverse `unmap_dense` (`umol-ast/remap.rs:165`) down beside them (`removed_nodes`/`removed_edges`
   data is already on `Remapping`). **Done**
2. `FactorOrdering` + `Unordered`/`Ordered`. **Done**
3. `RelationParticipant` + `ParticipantRefs`; impls for `NodeId`/`EdgeId`. **Done**
4. Generalize `FixedRelationSet`/`VarRelationSet` over `<P: RelationParticipant, O: FactorOrdering>`
   (today `NodeId` + unconditional sort): `new()` uses `O::canonicalize`. The remap-rebuild moves *off*
   `Remapping` (the `apply_to_*_relation_set` methods, `graph.rs:364`,`382`) onto the set as
   `apply_remapping(self, &Remapping) -> Self` — keeping the existing builder name (`builder.rs:91`):
   `filter_map` relations, `p.remap(m)` per participant (drop the relation on any `None`), then
   `Self::new` (canonicalize + incidence). The inverse needs no set method — the builder's `restore_*`
   path already assembles entries from undo payloads + survivors, switching `undo.atom` → `p.unmap(m)`.
   This keeps `Remapping` a pure removal-data carrier and delegates id-space routing to the participant
   trait. **Done**
5. Add all three birelation sets now (not stubbed — context-retention; foreseen consumers: a ports
   table → `FixedFixed`, through-space donor–acceptor → `VarVar`). Factors use `1`/`2` suffixes
   throughout (`L1`/`L2`, `N1`/`N2`, `O1`/`O2`, `factor_1()`/`factor_2()`; the chemistry layer aliases
   `site`/`ligands`). Build the node+edge union incidence (dedup `(key, rid)` per relation) and the
   `incident_edge` query (see Incidence). **Done**

### Forward-looking: relations as participants

Excluding `RelationId` is a present-scope decision, not permanent. Foreseen: **fragments become relations** and serve as
stereo-element *sites* (an element sited on a fragment is a relation whose site participant is another
relation); a hapto bond between an aromatic system and a metal is a softer example (it can instead be
flattened to atoms on both sides). The use is **acyclic and layered** (relation → relation → … → atoms;
no cycles), so removal stays a terminating cascade — not runaway recursion.

The enabler is additive, and the framework already anticipates it: a generic relation-set id (`SetId`),
a relation-removal dimension on `Remapping` (today `Remapping` is *one-sided* — base graph only; relation
removal lives a layer up in `IdRemapping`, since `RelationId`s are per-set), a `relation` arm on
`ParticipantRefs`, and a relation key-space in the incidence index. `RelationParticipant::remap` takes
`&Remapping` and `refs()` returns a struct — both grow without breaking. The discipline to keep now: reach
`Remapping` only through its methods, never by destructuring its fields, so the relation dimension can be
added later. The genuine fork (decide then): name fragments as first-class relations and accept the
cascade, vs flatten to atoms — no recursion but no named unit. A third path — promote ports/fragments
to *primitives* — obviates relation-refs entirely (relations stay flat over primitives; the recursion
relocates to base-graph containment); see doc 103, maximal-port model.

## Step 2 — stereo chemistry layer

Built on Step 1. raise is mechanical (no chemistry); perception is the chemistry.

- **raise** (mechanical): TableIR `chirality`/`wedge` → `#T`/`#C` constraint projections in the
  neighbor-order frame (a reindexing); no element. The aromaticity template — `raise` → `#a`
  (`AromaticValence`) constraint (`raise.rs:87–99`).
- **perception** (chemistry): a perception core + drivers, mirroring the aromaticity ops
  (`ops/transformer/aromatizer.rs`, `ops/resolver/aromaticity`): `StereoResolver` (≈ `AromaticResolver`;
  partial `#T`/`#C` → elements — the deliverable's driver), `StereoInferrer` (≈ `Aromatizer`;
  marker-free → elements — later), stereo transformation ops (≈ `Kekulizer` — later).

Reuse the relation-set / builder / remap / DSL-round-trip / transformer patterns and `AtomAutomorphism`;
do not rebuild. Stereo test corpus: ~150 files under `umol-graph/tests/{mol,smiles,sdf}_parsing`.

Phases A–E are in scope; F (3D) and G follow.

- **Phase A — types** (`umol-ast/src/ast/stereo.rs`, new):
  - **A1** — `StereoKind { Tetrahedral, CisTrans, SquarePlanar, TrigonalBipyramidal, Octahedral }` (flat `Copy`
    enum). All five carry their `~` involution and `ClassKey` mapping — the umol-ast-side `~` table — wired to
    `umol-perm`; perception (B–E) still scopes to TH/CT. **Done**
  - **A2** — the index types: `StereoConfigurationAst { Undetermined, NotStereo, Stereo(StereoIndexAst) }`;
    `StereoIndexAst { Undetermined, Lit(u32), Expr(Box<Expr>) }`; recursive `Expr { Lit(u32), Var(String),
    SwapOp(Box<Expr>), ApplyOp(Box<Expr>, Permutation), LitSet(Vec<u32>), VarDomain(String, Vec<u32>) }` (stereo's
    own, ≠ `value::Expr`; `ApplyOp` carries the `umol-perm` `Permutation` parsed from the `^`-image, not a
    Lehmer code; `LitSet`/`VarDomain` parsed-but-inert until sets land). `Lit` duplicated top-level + `Expr::Lit`
    by design; `Undetermined` kept out of `Expr` ⇒ `~+` unrepresentable; `Var` scope-resolved (declare/use
    collapsed); `~`/`^k` recurse. Config value = dense coset index per stereo class (`u32`), equivariant;
    equality = same index up to a common frame. `StereoConfigurationAst::simplify` (AST method,
    never in the parser): lifts `Expr(Lit)`→`Lit` (folds the by-design duplication) and reduces closed
    operator-exprs (`~1`→the enantiomer coset, `~~e`→`e`) via `umol-perm` (so depends on A′); free-`Var`
    exprs are left as-is. **Done**
  - **A3** — payloads `StereoAtomAst` / `StereoBondAst`, each `{ kind: StereoKind, configuration:
    StereoConfigurationAst, constraints }` — predicate only, no site/ligands. **Done**
  - **A4** — per-site constraints `StereoAtomConstraint {}` / `StereoBondConstraint {}` (empty) + `…Constraints`
    collections + empty `…ConstraintDsl`, mirroring `NoncovalentBond` (split — projected structure differs per site).
    **Done**
  - **A5** — ligand participant `StereoLigand(NodeId, StereoLigandKind)` + `StereoLigandKind { Atom,
    ImplicitHydrogen, LonePair }`, with its `RelationParticipant` impl (route the inner `NodeId` via
    `map_node`/`unmap_node`, kind carried; `refs` = that node; stored as `NodeId`, views expose `AtomId`). **Done**
  - **A6** — local-frame constraints `AtomConstraint::TetrahedralStereo(StereoConfigurationAst)` /
    `BondConstraint::CisTransStereo(StereoConfigurationAst)` in `ast/constraint/{atom,bond}.rs` (`#T`/`#C`,
    uppercase derived-predicate namespace, `Th`→`#T`/`Ct`→`#C`; same arg as the element config, frame differs). **Done**
  - **A7** — ids via `define_id!` in `ast/ids.rs`. **Done**

- **Phase A′ — coset algebra** (`umol-perm`, new crate; pure permutation algebra — no chemistry, no geometry,
  no AST deps; a dependency of umol-graph). The dense coset index **is the OpenSMILES arrangement number**
  (`@TH/@AL/@SP/@TB/@OH`) — reproduced, never re-invented. **Built complete for all five geometries now** —
  it's a leaf crate, so finishing it here means later non-tetrahedral *perception* extends only the chemistry
  layers (C/D/E), never this one (B–E themselves still scope to TH/CT).
  - **A′1** — `Permutation` (one permutation, n ≤ 6, `Copy`, one-line notation): `identity`/`from_image`/
    `between(from, to)` (the relabel τ) / `apply`/`compose`/`inverse`/`sign`/`act`/`rank`/`unrank` (Lehmer —
    internal canonical ordering only); `Ord`+`Eq`+`Hash`. **Done**
  - **A′2** — `PermutationGroup` (subgroup of Sₙ, full enumeration): `generate(degree, gens)` (brute-force
    closure), `symmetric`/`alternating`/`cyclic`/`dihedral`, `order`/`contains`/`elements`. **Done**
  - **A′3** — `CosetSpace` = R (a `PermutationGroup`) + the coset **partition** (one canonical rep per coset):
    `count` (= n!/|R|), `coset_rep`. The algebraic layer — which orderings are equivalent, not yet a numbering.
    **Done**
  - **A′4** — `ClassKey` (`FromStr`/`Display` key, à la `SchoenfliesSymbol`): families `Sym/Alt/Cyc/Dih(u8)`
    (natural action) + geometry variants `Tetrahedral`/`CisTrans`/`SquarePlanar`/`TrigonalBipyramidal`/
    `Octahedral`. `static REGISTRY: LazyLock<Mutex<HashMap<ClassKey, &'static CosetSpace>>>`, `space(k)` =
    build-once + `Box::leak` (mirrors umol-msym `point_group.rs:71`); `Coset { space: &'static CosetSpace,
    index: u32 }` ties identity to the interned space (like `SymmetryOp`). Const generators per class:
    `Tetrahedral`→`alternating(4)` (A₄), `CisTrans`→double-swap (Z₂), `SquarePlanar`→`dihedral(4)` (D₄),
    `TrigonalBipyramidal`→D₃ on 5, `Octahedral`→O on 6 — all five built now. **Done**
  - **A′5** — the **index alignment** `CosetSpace::index(perm) -> u32` = the OpenSMILES decomposition for **all
    five** classes, transcribed from the spec and validated so its fibers equal R's coset partition (A′3):
    - **TH/CT** — the parity (the `sign` bit).
    - **SP** — the `U/4/Z` path shape (§3.8.5) → `@SP1–3`.
    - **TB** — apical-pair axis + 3-equatorial winding (§3.8.6) → `@TB1–20`.
    - **OH** — axis + 4-equatorial shape + winding (§3.8.7, recursive on the SP shapes) → `@OH1–30`.
    No ordering of our own — each is the spec's table, fiber-checked against R. Exhaustive round-trip tests
    (`index ∘ unindex = id` over all `n!` permutations) per class.
  - **A′6** — `reindex(k, input_index, τ) -> u32` (= `index(τ ∘ unindex(input_index))`) + the `~`/`^image`
    operator action on indices (over Sₙ/core(R); only `~` for the binary classes in scope, Phase G). raise
    (B3/B5) and perception (C) call `space(k).reindex(…)`. **Done**
  - **A′ as-built (verified against the spec, not assumed — load-bearing for Phase B):**
    1. **Coset side is the right coset `Rσ`.** A config `σ` is ligand→position; two configs are the same
       arrangement iff related by a rotation acting on positions (left mult `r∘σ`), so `coset_rep(σ) =
       min_{r∈R} r∘σ`. (Left coset `σR` is wrong — the §3.8.5 shape groups are exactly `Rσ`.)
    2. **`reindex` is a 2-arg method on the interned space:** `space(k).reindex(index, relabeling) =
       index(unindex(index) ∘ relabeling)`, with `relabeling = Permutation::between(order_in, order_out)`
       (`relabeling(i)` = position in `order_in` of `order_out`'s i-th neighbor). Supersedes the 3-arg sketch
       in A′6. B3/B5 pass `between(tableir_order, incidence_order)`.
    3. **Shape = diagonal-step count of the path** `σ(0)→σ(1)→σ(2)→σ(3)` on the square (diagonals `{0,2}`,
       `{1,3}`): **U = 0, Z = 1, 4 = 2**. With `@SP1/2/3 = U/4/Z` (§3.8.5), `@SP2`(4) has two diagonal steps and
       `@SP3`(Z) one — the opposite of the naïve reading.
    4. **`@`/`@@` are the two C₄ orbits of a shape**, *not* an ordering and its string-reverse (a reverse can be
       `ρ²` of the original — same orbit). OH equatorial reps: U `1234`/`4321`, Z `2314`/`2134`, 4 `2413`/`1324`.
    5. **Position frames match OpenSMILES:** TB axial = `0,4` / equatorial `1,2,3`; OH axis `0,5` / equatorial
       square `1,2,3,4` (cyclic, diagonals `1–3`,`2–4`); SP `dihedral(4)` already matched. TH/CT use the generic
       `CanonicalRank` numbering (= parity for a 2-coset space), not a bespoke parity decomposition.
    Validated: 7 equivalence-SMILES `reindex` cases (3 TB §3.8.6, 4 OH §3.8.7) + `fibers == cosets` +
    `index ∘ unindex` round-trip; 61 tests, clippy clean.

- **Phase B — raise → `#T`/`#C`** (mechanical; `table_ir/raise.rs`, mirrors the aromaticity→`#a` pass).
  Reads TableIR per-atom chirality / per-bond stereo, reindexes the input arrangement into umol's incidence
  frame, writes `#T`/`#C`; builds **no** element (Phase C).
  - **B1** — `incidence_ligand_order` on atom/bond views (= `stereo_ligands()`): real neighbors in adjacency
    order, then implicit-H, then lone-pairs.
  - **B2** — `tableir_ligand_order`: the arrangement in TableIR's neighbor frame, per the conventions below.
  - **B3** — `space(k).reindex(input_index, Permutation::between(input, inc))` (umol-perm, A′6); 2-coset ⇒
    the parity XOR. τ's only content is the virtual-ligand repositioning (reals keyed by identity).
  - **B4** — tetrahedral → `AtomConstraint::TetrahedralStereo` (pseudocode below).
  - **B5** — cis/trans → `BondConstraint::CisTransStereo`; τ splits per sp² carbon
    (`coset = input_coset ^ swap(C1) ^ swap(C2)`, incidence order = C1's neighbors then C2's).
  - **B6** — unspecified (wavy / `@?`) ⇒ `Stereo(Undetermined)` = `#T+`/`#C+`; absent ⇒ no constraint. raise
    never emits `*`/`!`/`Var`/operators (those are pattern-side).
  - **B7** — corpus tests over the conformance suite (`@`/`@@`, `/`,`\`, MOL wedges, CXSMILES).

  **Spec.** Config value = the **dense coset index per stereo class** — the SMILES arrangement number
  (`@TH1-2` … `@OH1-30` = n!/|R| cosets, `u32`; not a Lehmer rank). Verified input conventions
  (`materials/formats`):
  - **SMILES** `@`/`@@` = counterclockwise/clockwise viewed from the first-listed neighbor; order = from-atom
    then as-written; **implicit H immediately after the from-atom** (first if none). `/`,`\` = cis/trans of
    the flanking single bonds relative to the carbon (`F/C=C/F` trans, `F/C=C\F` cis).
  - **MOL V2000**: atom **parity field ignored** — stereo = the **wedge** (1 Up / 6 Down, narrow end at the
    center atom) + neighbors' **2D coords**, implicit H = the missing direction; double bond 0 = coords,
    3 = either. **CXSMILES**: `wU`/`wD`/`w` wedges, `@:`/`@@:` parity lists.

  ```text
  fn raise_tetrahedral(a: &table_ir::Atom) -> Option<AtomConstraint> {   // B4
      let cfg = match a.chirality? {                   // chirality: Option<Chirality>; None ⇒ no #T
          Chirality::Unspecified =>                    // @? / wavy — stereogenic, config undefined  (B6)
              Stereo(StereoIndexAst::Undetermined),    // #T+
          c @ (Chirality::CounterClockwise | Chirality::Clockwise | Chirality::Tetrahedral { .. }) => {
              let input_coset = tetra_index(c);        // @=TH1 / @@=TH2 shorthand, or explicit `arr` (the SMILES @TH number)
              let input = tableir_ligand_order(a);     // B2: from-atom, implicit-H after it, then as-written
              let inc   = incidence_ligand_order(a);   // B1: adjacency order; H then LP trailing
              let coset = space(ClassKey::Tetrahedral)                       // B3 — umol-perm (Phase A′)
                  .reindex(input_coset, Permutation::between(&input, &inc));
              Stereo(StereoIndexAst::Lit(coset))       // #T<coset>  (top-level Lit)
          }
          // Allenal / SquarePlanar / TrigonalBipyramidal / Octahedral — out of scope, parse-opaque
          _ => return None,
      };
      Some(AtomConstraint::TetrahedralStereo(cfg))
  }
  // CT is the bond analogue — read from Bond.stereo / wedge (not Chirality); space(ClassKey::CisTrans).reindex(…).
  ```
- **Phase C — perception → birelations** (chemistry; `umol-graph/src/ops/stereo.rs`). Lifts `#T`/`#C` +
  topology → the stereo birelation overlay. **Same-frame lift** — the element's `:ligands` order is
  `stereo_ligands()` is `#T`'s incidence frame (D2: equivariant, not canonical), so config copies through with
  no reindex; marker-free inference / meso / canonicalization stay in G.
  - **C1** — storage: `MoleculeAst` fields `stereo_atoms: StereoAtom` (= `FixedVarBirelationSet<NodeId,
    Ordered, 1, StereoLigand, Ordered, StereoAtomAst>`) + `stereo_bonds: StereoBond` (`EdgeId`-site analogue);
    ids `StereoAtomId`/`StereoBondId` (`define_id!`); extend `from_parts`/`Clone`/`PartialEq`.
  - **C2** — builder `add_stereo_atom` / `add_stereo_bond`, mirroring the four existing overlay builders.
  - **C3** — remap/restore: apply Step-1's `RelationParticipant` generalization to both tables —
    `apply_remapping` (`remap.node`→`p.remap`, `builder.rs:102`/`204`), `restore_participants`
    (`undo.atom`→`p.unmap`, `:270`/`277`), and `removed_stereo_*` in `IdRemapping`/`UndoRemapping` (hand-written
    per kind, no abstraction).
  - **C4** — `stereo_ligands()` on atom/bond views: real neighbors in adjacency order, then implicit-H
    (× H-count), then lone-pairs, each a `StereoLigand` (materializes the virtual ligands). Shared with B1;
    LP minimal/deferred for the organic deliverable.
  - **C5** — perception core (one element): from an atom/bond + its `#T`/`#C` + `stereo_ligands()`, build the
    `StereoAtomAst`/`StereoBondAst` element (focus, ligand set, kind, config copied through — same frame); plus
    the inverse projection element→`#T` (re-derivable, for write-back).
  - **C6** — `StereoResolver` (≈ `AromaticResolver`, `ops/resolver/aromaticity` template): scan a molecule's
    `#T`/`#C`, drive C5, populate both tables. Marker-driven and partial — only marked atoms/bonds become
    elements; marker-free `StereoInferrer` stays in G.
  - **C7** — tests: raise→resolve round-trip over the ~150-file corpus; assert focus/ligands/config of the
    perceived birelations for known R/S and E/Z cases.
  - **C8** - validators "Th3", "Ct3", ... are disallowed (tier-2 checks, doc 86)

- **Phase D — DSL round-trip** (`umol-ast/src/dsl/stereo.rs`, new; wired into `dsl/molecule.rs`, mirror
  `aromatic`). Faithful round-trip, no frame conversion (config stored relative to the written `:ligands`
  order; `#T`/`#C` relative to the local frame).
  - **D1** — config-string parser/writer: `class config` ↔ `(StereoKind, StereoConfigurationAst)`. Head
    `Th`/`Ct`; config `* | ! | + | <coset-term>`, coset-term recursive (`nat`→`Lit`, `?id`→`Expr(Var)`,
    `~e`→`Expr(SwapOp)`, `e^<image-number>`→`Expr(ApplyOp(perm))`); `Expr::LitSet`/`VarDomain` (`{…}`, `?o :: {…}`) reserved at the
    surface (deferred with non-tetrahedral). **One function** — D3's `:type` head and D5's `#T`/`#C` call it. **Done**
  - **D2 — AST + DSL stereo types** (the predicate/element split, doc 103): the `#T`/`#C` *predicate*
    is `StereoConfigurationAst` (`* | ! | Stereo(StereoCosetAst)` + full lattice); the stereo *element*
    carries only `StereoCosetAst` (always stereogenic — no `NotStereo`). The ligand EDN surface
    (`atom-ref | [:h ref] | [:lp ref]`, reserved `[:bond/:port/:fragment ref]`, unknown tags rejected),
    the `:stereo-atoms`/`:stereo-bonds` entries, and `MoleculeDsl` raise/lower wiring are **D3**, not here.
    - **D2a** — rename + field split. `StereoIndexAst → StereoCosetAst` (incl. inside
      `StereoConfigurationAst::Stereo(_)`); on `StereoAtomAst`/`StereoBondAst`, `configuration:
      StereoConfigurationAst` → `coset: StereoCosetAst`. `into_ground`/`into_zeroed` become no-ops
      (ground iff the coset is ground — no `NotStereo` to coerce). **Done**
    - **D2b** — `StereoCosetAst` (`Undetermined | Lit(u32) | Expr`): `AsLit` ✓, `Lattice` ✓,
      constructors ✓; **finish `matches_value()`**. **Done**
    - **D2c** — `StereoConfigurationAst` (`Undetermined | NotStereo | Stereo(StereoCosetAst)`): `AsLit`
      ✓, `Lattice` ✓, `From<u32>`/`From<Vec<u32>>` ✓ — rename ripple only. **Done**
    - **D2d** — `StereoAtomAst`/`StereoBondAst` (macro; fields `kind`, `coset`, `constraints`): add
      **`Lattice`** — per-kind committed-(A): `meet` cross-kind = `None`, `join`
      `debug_assert!`s equal kinds, `is_ground`/`is_undetermined`/`matches` over the three fields. **Done**
    - **D2e** — `StereoAtomConstraint`/`StereoBondConstraint` (`ast/constraint/stereo.rs`): uninhabited
      today; **trivial `Lattice`** (and the `StereoAtomConstraints`/`StereoBondConstraints` collections). **Done**
    - **D2f** — DSL `FromAst`/`IntoAst`. `StereoAtomDsl`/`StereoBondDsl` ↔ `StereoAtomAst`/`StereoBondAst`
      (`FromStr`/`Display`/`FromEdn`/`ToEdn` ✓; **add `FromAst`/`IntoAst`**, trivial now — no `NotStereo`
      default, `into_ground` no-op). `StereoAtomConstraintDsl`/`StereoBondConstraintDsl` ↔ the uninhabited
      AST constraints (`FromEdn`/`ToEdn` ✓; **add `FromAst`/`IntoAst`**). **Done**
  - **D3 — molecule overlays + EDN entry surface** (two sets: the `MoleculeAst` overlays in
    `ast/molecule.rs`, then the EDN entry surface in `dsl/molecule.rs`).
    - **D3a** — *AST overlays + `from_parts` construction.* Two `MoleculeAst` birelations:
      `stereo_atoms: FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, StereoAtomAst>`
      and `stereo_bonds: FixedVarBirelationSet<EdgeId, Ordered, 1, StereoLigand, Ordered, StereoBondAst>`
      (fixed arity-1 site → ordered ligands; payload the D2 element). Build them in `from_parts`
      (inputs `Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)>` / `Vec<(BondId, Vec<StereoLigand>,
      StereoBondAst)>`), plus `Clone` / `PartialEq` / derived `Default`. `from_arcs` defaults the two
      overlays to empty for now; the builder / `edit()` carry is split into D3h+. **Done**
    - **D3b** — *AST accessors + indexing.* Immutable and mutable accessors per relation;
      `Index<StereoAtomId>` / `Index<StereoBondId>` for `MoleculeAst`; include both relations in the
      subgraph definition. **Done**
    - **D3c** — *AST predicates + transforms.* `has_stereo_atoms()` / `has_stereo_bonds()` (folded into
      `has_overlays()`); extend `is_ground`, `simplify_values()`, `lift_constraints()`, and
      `inline_constraints()` over the two relations. **Done**
    - **D3d** — *DSL entry grammar + inputs.* Entry map `{ [:id <keyword>] :site <ref> :ligands
      [<ligand> …] :type <element-dsl> }` under the `:stereo-atoms` / `:stereo-bonds` keys. `:site` =
      atom-ref (atoms) or bond-ref (bonds); `:type` = the D1/D4 element string or keyword. Each
      `<ligand>` is a plain `<atom-ref>` (→ `StereoLigandKind::Atom`) or a keyword-headed vector
      `[:h <atom-ref>]` (→ `ImplicitHydrogen`) / `[:lp <atom-ref>]` (→ `LonePair`), the ref being the
      host atom. `StereoAtomEntryInput` / `StereoBondEntryInput`; add to `MoleculeInput` and `Metadata`.
      **Done**
    - **D3e** — *DSL parse + render.* A `read_stereo_ligand` / `parse_stereo_ligand` dispatching
      `<atom-ref>` vs `[:h <ref>]` / `[:lp <ref>]` → `StereoLigand`; `read_stereo_atom_dsl` (the
      `:type`) + `read_stereo_atom_entry` (the whole entry, streaming) and `parse_stereo_atom_entry`
      (tree), and the bond forms; wire into `parse_molecule_input`; emit both entry kinds (ligands
      back to ref / `[:h …]` / `[:lp …]`) in `render_molecule_edn`. **Done**
    - **D3f** — *`MoleculeDsl` raise/lower.* Thread the stereo entries through `FromAst`/`IntoAst` —1
      ref → id resolution in `into_ast`, the ligand list resolved like `:atoms`. Also add the
      **noncovalent bonds**, which are currently missing from `MoleculeDsl`'s `FromAst`/`IntoAst`.
    - **D3g** — *Tests.* EDN↔AST round-trip for both entry kinds (string and `:ccw`/`:z` payloads;
      Atom / `[:h]` / `[:lp]` ligands); ref resolution + unknown-ref errors; the new overlay accessors
      / predicates / `is_ground` / `simplify_values`. **Done**
    - **D3h** — *Builder carry + remap.* `from_arcs` + `MoleculeBuilder::from_parts` + `edit()` gain
      the two stereo Arcs; the builder stores them (shared storage) and round-trips them through
      `build`; `remove()` applies the node/edge remap so stereo refs stay valid after structural edits.
      Add StereoAtomRef and StereoBondRef to edit.rs. **Done**
    - **D3i** — *Builder stereo editing + undo* (full noncovalent-parity edit/undo surface; widened
      from the original "cascade-undo only" scope because undo is incomplete without the matching
      add/remove edits, and those need mutable storage). Pulls the `FixedVarSetStorage` enum forward
      from D3j — its first genuine consumer is exactly these mutators. Pieces:
      - `FixedVarSetStorage` (Shared/Mutable CoW) for the two stereo birelations; D3h's plain-Arc
        builder fields convert to it; `from_parts` / `build` route through it.
      - Builder mutators `add_stereo_atom` / `remove_stereo_atoms` / `add_stereo_bond` /
        `remove_stereo_bonds` (+ `remove_added_stereo_*` undo helpers).
      - `edit.rs`: `AddedStereoAtom` / `RemovedStereoAtom` (+ bond) — each carrying **site + ligands**
        (two factors), unlike the single-`atoms` overlays; `Edit::{AddStereoAtom, RemoveStereoAtom,
        AddStereoBond, RemoveStereoBond}`; `Undo::{RemoveAddedStereoAtom, RestoreRemovedStereoAtom,
        …Bond}`; `RemovedOverlays` += `stereo_atoms` / `stereo_bonds`.
      - `remap.rs`: `IdRemapping` / `UndoRemapping` gain `removed_stereo_atoms` / `removed_stereo_bonds`
        + accessors (ripples to `new` / `relations` / `empty` and their ~19 call sites).
      - `builder.rs`: `restore_stereo_atoms` / `restore_stereo_bonds` (+ singular wrappers); cascade
        capture so topology removal records dropped stereo into `RemovedOverlays` and `restore_topology`
        restores it.
      - `transact.rs`: dispatch the four new `Edit`s → apply + `Undo`; capture stereo in
        `capture_removed_topology`; undo dispatch for the four new `Undo`s.

      *Status: implementation landed, tree green.* Two carryovers to D3j: (1) the transactional
      `RemoveStereo*` apply + `capture_removed_topology` read the current element through a temporary
      `MoleculeBuilder::stereo_atom_entry` / `stereo_bond_entry` owned-snapshot stopgap (marked
      `TODO(D3j)`), because the builder read **views** were scoped to D3j; (2) D3i adds no tests of its
      own — they live with the view migration in D3j so they exercise the final (view-based) read path,
      not the stopgap.
    - **D3j** — *Builder views + field edits.* `BuilderView` / `BuilderViewMut` for stereo (read +
      field-mutation access), `SetStereoAtomField` / `SetStereoBondField` edits and their `Undo`s, plus
      remaining ergonomics. (The storage enum, add/remove mutators, and `Added*` / `Removed*` edit
      types moved to D3i.) Also folds in the D3i carryovers:
      - **Replace the D3i read stopgap.** Add `StereoAtomBuilderView` / `StereoBondBuilderView`
        (mirror `NoncovalentBondBuilderView`); switch the transact `RemoveStereo*` apply and
        `capture_removed_topology` to read through them (clone the `ast` off a borrow, like the sibling
        overlays); delete `MoleculeBuilder::stereo_atom_entry` / `stereo_bond_entry` and
        `FixedVarSetStorage::data` if then unused.
      - **Tests deferred from D3i** (write against the view-based path, not the stopgap):
        - Transactional `AddStereoAtom` / `AddStereoBond` + rollback restores the prior (stereo-free)
          state; round-trip via `transact` then `Transaction::rollback`.
        - Transactional `RemoveStereoAtom` / `RemoveStereoBond` + rollback restores the element
          (site / ligands / coset intact); `OldStateMismatch` on a stale recorded old-state.
        - Topology-removal cascade: a `transact` `RemoveTopology` that drops a stereo element (site or
          ligand atom, or site bond removed) and whose rollback restores it at its original id.
        - `IdRemapping` / `UndoRemapping` `stereo_atom` / `stereo_bond` accessor unit coverage
          (remap shift past removed indices; inverse `unmap`). **Done**
    - **D3k** - Verify the stereo entities are correctly included in the SubpatternAnchor and remapping
      data structs / APIs. **Done**
    - **D3l** — *`StereoAtomView` / `StereoBondView` relational + structural accessors.* Bring the stereo
      views to parity with the sibling element views. Conventions: bare noun → view, `_id`/`_ids` →
      id(s), `has_*` → bool — the lone exception is `ligands()`, which returns the inherent
      `&[StereoLigand]` (ligands have their own type; their atom projections are `ligand_atoms*`). Site
      lookups are unique (`Option`, precedent: `bonds().connecting`); all other lookups return
      iterators. Virtual ligands (`ImplicitHydrogen` / `LonePair`) carry their **bearing atom** as
      `.atom()` — the site for stereo atoms, the relevant end atom for stereo bonds. `#[inline]` on
      inherent-field getters. `coset_for` / `permutation_for` take `&[StereoLigand]` and return `Option`
      (None unless a permutation of `ligands()`); they live on the view because only it carries the
      ligand order **Done**
    - **D3m** — *Review relational constraints for stereo.* `RelationalConstraint`
      (`ast/constraint/relational.rs`) and its surface `RelationalConstraintDsl` (`dsl/relational.rs`)
      are DAMN-only (dative / aromatic / multicenter / noncovalent). Decide which stereo variants to
      add — site identity (`StereoAtomSite` → atom, `StereoBondSite` → bond), ligand-set membership /
      equality (`StereoAtom{Contains,Ligands}`, bond analogs), and ligand role predicates
      (`…AllLigands` / `…AnyLigand` delegating an `AtomConstraint`) — plus the `:<entity>-<role>` EDN
      keys in `RELATIONAL_KEYS`, `from_ast` / `into_ast`, and `simplify` / `remap`. **Done**
    - **D3n** Fix naming in remap.rs, index -> id, indices -> ids, review all field and method names **Done**
    - **D3o** Add coset constraint to StereoAtomConstraints and StereoBondConstraints **Rejected**
  - **D4** — sugar `:ccw`/`:cw`/`:e`/`:z` (each carries its class — `Th1`/`Th2`/`Ct1`/`Ct2`) ↔ the `:type` head.
    **Done**
  - **D5** — `#T`/`#C` atom/bond-string surface: the derived-predicate tokens in the existing atom/bond
    constraint-string parser (`dsl/constraint.rs`) — `#T<config>`/`#C<config>` (local-frame, the **same** D1
    `StereoConfigurationAst` parser) inside the atom-string (`C#h#T1`) / bond-string.
    Add `tetrahedral_stereo` method to atom view, `cis_trans_stereo` to bond view 
    (+ `tetrahedral_stereo_ligands` / `cis_trans_stereo_ligands` as a separate method or let `*_stereo` return
    a tuple).
    is_in_stereo_atom(), stereo_atoms(), stereo_atom_ids(), also separate query methods for site/ligands, names?
    is_in_stereo_bond(), stereo_bonds(), stereo_bond_ids(), also query methods incident atoms + ligands, names?
    add stereo atoms and stereo bonds to is_in_overlays.
    Add `#T/#C` to `derive_constraints()` method. **Done**
  - **D6** - add stereo atoms and bonds to SubPatternAnchorDsl **Done**
  - **D7** — round-trip tests: EDN↔AST for both surfaces (stereo elements *and* `#T`/`#C` strings) over the
    ~150-file corpus, under `--features conformance`. **Replaced by full conformance tests**
  - **D8** - macros stereo_atom!, stereo_atom_ground!, stereo_atom_zeroed!, stereo_bond!, stereo_bond_ground!,
     stereo_bond_zeroed! in `macros.rs`. **Done**
  - **D9** - update specifications in umol-dsl-spec.md **Done** (top-level keys fixed; grammar non-terminals
      aligned to key names; stereo elements + `#T`/`#C` constraints + relational + anchor + §7.14 subgrammar added)
  - **D10** - add to prop test and fuzzing **Done** (stereo elements + `#T`/`#C` + relational + anchor in `molecule_ast_strategy`/`constraint_leaf_strategy`/`sub_pattern_anchor_strategy`; stereo entity-string + keyword roundtrip tests; `parse_stereo_atom`/`parse_stereo_bond` in `fuzz_entity_strings`)

- **Phase E — matching** (the stereo ASTs' `AsLit` + `Lattice` impls — not a bespoke matcher; the existing
  substructure matcher is reused, and `umol-perm` enters exactly once, at the frame alignment).
  - **E1** — `AsLit` (trivial): `StereoConfigurationAst::as_lit` = the resolved config when `is_ground`
    (`Stereo(Lit)` / `NotStereo`), else `None`; likewise `StereoIndexAst`/`Expr`.
  - **E2** — `Lattice` on `StereoConfigurationAst` (+ `StereoIndexAst`/`Expr`; the `StereoAtomAst`/`StereoBondAst`
    payloads delegate — kind exact-match + config). The config Hasse: `*` top; `!` and `+` incomparable middles;
    `Stereo(Lit(k))` ground cosets. `is_undetermined`/`is_ground`/`meet`/`join`/`matches` — the **in-frame**
    partial order, `*` the wildcard (`join(Th1,Th2)=+`, `join(Th1,!)=*`, `meet(+,Th1)=Th1`, `meet(Th1,Th2)=None`).
    Free `Var`/`Expr` configs are **top** in the pure ops (full coset domain): `meet(?o,X)=X`, `join(?o,X)=*`;
    the variable's cross-occurrence correlation is NOT here — it is `unify`'s (E5). Closed operator-exprs reduce
    first (A2 `simplify`). **In-scope assumption:** every stereo variable is free (full domain) ⇒ `meet` returns
    the other operand; when `LitSet`/`VarDomain` land (non-TH), `meet` becomes a genuine domain intersection — the
    single leaf that grows. `meet`/`join` stay pure (never thread the env).
  - **E3** — frame alignment: the substructure matcher's pattern→target ligand correspondence → τ
    (`Permutation::between`) → reindex the pattern config to the target frame (`umol_perm::space(k).reindex`) →
    then E2's `matches`. The single point `umol-perm` enters E.
  - **E4** — `#T`/`#C` matching: `matches` on `AtomConstraint::TetrahedralStereo` / `BondConstraint::CisTransStereo`
    (the SMARTS-style atom/bond query), reusing E2 + E3 in the local-incidence frame.
  - **E5** — `unify` (relative stereo `?o`/`~?o`): the binding-aware match. Added to `Lattice` as
    `matches_capturing(target, &mut Env) -> bool` with a default body = `matches` (env ignored); stereo overrides
    it. It solves for the variable by inverting the operator chain against the frame-aligned target (`?o`→`o=c`;
    `~?o`→`o=~c`; `?o^k`→`o=c^{-k}`, via `umol-perm`), then unifies into a molecule-wide typed `Env` (bind if free,
    else equality-check), threaded as a post-correspondence pass over the matched stereo elements; the candidate
    matches iff every variable unifies. Deterministic (operators invertible) — the first live exercise of
    binding-correlation, and the pilot for the deferred cross-AST **Var/Bind sweep** (which flips `matches_capturing`
    to canonical with `matches` its empty-env wrapper, and grows `Env` from singletons to domains/relations for
    non-stereo variable types — `meet`/`join` stay pure throughout).
  - **E6** — tests: pattern↔target over the corpus — TH/CT absolute + relative-stereo binds; assert match/no-match.
- **Phase F — 3D (deferred).** umol-geometric (greenfield): config constrains local geometry at
  embedding (signed volume / cis-trans side). Substrate = doc 071.
- **Phase G — later.** Config operators `~`/`^` (relative stereo); deeper perception (stereogenicity on
  flat input); canonicalization / meso (derived `A_s` pass + para-stereo fixpoint).

## Naming (umol-ast vocabulary)

| concept | type |
| --- | --- |
| site tables | `StereoAtom` / `StereoBond` (aliases over `FixedVarBirelationSet`) |
| top-level DSL keys | `:stereo-atoms` / `:stereo-bonds` |
| payloads | `StereoAtomAst` / `StereoBondAst`, each `{ kind: StereoKind, configuration: StereoConfigurationAst, constraints }` |
| per-site constraints | `StereoAtomConstraint {}` / `StereoBondConstraint {}` (empty; + `…Constraints` collections, empty `…ConstraintDsl`) |
| shape | `StereoKind` (`Tetrahedral`, `CisTrans`, …) |
| configuration | `StereoConfigurationAst { Undetermined, NotStereo, Stereo(StereoIndexAst) }` |
| config value | `StereoIndexAst { Undetermined, Lit(u32), Expr(Box<Expr>) }` over a recursive stereo `Expr { Lit, Var, SwapOp, ApplyOp, LitSet, VarDomain }` (its own, ≠ `value::Expr`); `~`/`^image` recurse (`~1`, `0^2134` sayable; `ApplyOp` holds a `Permutation`); `Undetermined` out of `Expr` ⇒ `~+` unrepresentable; `Lit` duplicated by design; `LitSet`/`VarDomain` deferred; index = `u32` dense coset index per class (SMILES arrangement number, not Lehmer) |
| ligand | `StereoLigand` |
| ligand slot | `StereoLigandKind` { `Atom`, `ImplicitHydrogen`, `LonePair` } |
| local constraints | `AtomConstraint::TetrahedralStereo` / `BondConstraint::CisTransStereo` |
| DSL form | `StereoElementDsl` |

DSL key convention (whole molecule map): **`key` = kebab-plural of the type** — `:atoms`, `:bonds`,
`:dative-bonds`, `:aromatic-systems`, `:multicenter-bonds`, `:noncovalent-bonds`, `:stereo-atoms`,
`:stereo-bonds` (renamed from the bare-adjective `:dative`/`:aromatic`/… forms). Keys are per **site**,
not per kind — site is bounded/structural (one key per overlay, like the others); the kind stays in the
config head (`Th1`/`Ct2`).

## Notation (consistent within each layer)

- graph-core: `remap`/`unmap`; `Remapping::map_node`/`map_edge` (forward) + `unmap_node`/`unmap_edge` (inverse).
  Mechanical coordinate translation.
- umol-ast: `restore_participants` / `UndoRemapping` / `restore_*`, delegating down to `unmap_*`.
  Transaction-unwind semantics.

## Files

New: `umol-ast/src/ast/stereo.rs`, `umol-ast/src/dsl/stereo.rs`, `umol-graph/src/ops/stereo.rs`.
Modified — graph-core: `relation.rs`, `graph.rs`; umol-ast: `ast/{ids,molecule,remap}.rs`,
`ast/molecule/builder.rs`, `ast/views/{atom,bond}.rs`, `ast/constraint/{atom,bond}.rs`,
`dsl/molecule.rs`; umol-graph: `table_ir/raise.rs`; umol-geometric: `src/{coordinates,molecule}.rs`
(Phase F).

## Verification

Per-phase rstest table tests (relation/birelation generics + remap/unmap; types; raise `#T`/`#C` on
corpus; perception constraints→birelations; DSL round-trip; constraint↔element cross-check; matching).
End-to-end (A–E): `F[C@H](Cl)Br` / `F/C=C/F` → `#T`/`#C` → perceive → serialize `:stereo` → re-parse →
match. Run `mol/smiles/sdf_parsing` + resolution conformance (`--features conformance --test resolution`).

## Open

- Whether the `#T`/`#C` constraint is a stored `AtomConstraint`/`BondConstraint` kind (recommended,
  symmetric with `#a`) vs raise-transient.
