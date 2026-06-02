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
  - **A1** — `StereoKind { Tetrahedral, CisTrans }` (flat `Copy` enum, no algebra; further kinds added when built).
  - **A2** — the index types: `StereoConfigurationAst { Undetermined, NotStereo, Stereo(StereoIndexAst) }`;
    `StereoIndexAst { Undetermined, Lit(u32), Expr(Box<Expr>) }`; recursive `Expr { Lit(u32), Var(String),
    SwapOp(Box<Expr>), ApplyOp(Box<Expr>, u32), Set(Vec<u32>), VarDomain(String, Vec<u32>) }` (stereo's own,
    ≠ `value::Expr`; `Set`/`VarDomain` parsed-but-inert until sets land). `Lit` duplicated top-level + `Expr::Lit`
    by design; `Undetermined` kept out of `Expr` ⇒ `~+` unrepresentable; `Var` scope-resolved (declare/use
    collapsed); `~`/`^k` recurse. Config value = dense coset index per stereo class (`u32`), equivariant;
    equality = same index up to a common frame. `StereoConfigurationAst::simplify` (AST method,
    never in the parser): lifts `Expr(Lit)`→`Lit` (folds the by-design duplication) and reduces closed
    operator-exprs (`~1`→the enantiomer coset, `~~e`→`e`) via `umol-perm` (so depends on A′); free-`Var`
    exprs are left as-is.
  - **A3** — payloads `StereoAtomAst` / `StereoBondAst`, each `{ kind: StereoKind, configuration:
    StereoConfigurationAst, constraints }` — predicate only, no site/ligands.
  - **A4** — per-site constraints `StereoAtomConstraint {}` / `StereoBondConstraint {}` (empty) + `…Constraints`
    collections + empty `…ConstraintDsl`, mirroring `NoncovalentBond` (split — projected structure differs per site).
  - **A5** — ligand participant `StereoLigand(NodeId, StereoLigandKind)` + `StereoLigandKind { Atom,
    ImplicitHydrogen, LonePair }`, with its `RelationParticipant` impl (route the inner `NodeId` via
    `map_node`/`unmap_node`, kind carried; `refs` = that node; stored as `NodeId`, views expose `AtomId`).
  - **A6** — local-frame constraints `AtomConstraint::TetrahedralStereo(StereoConfigurationAst)` /
    `BondConstraint::CisTransStereo(StereoConfigurationAst)` in `ast/constraint/{atom,bond}.rs` (`#T`/`#C`,
    uppercase derived-predicate namespace, `Th`→`#T`/`Ct`→`#C`; same arg as the element config, frame differs).
  - **A7** — ids via `define_id!` in `ast/ids.rs`.
- **Phase A′ — coset algebra** (`umol-perm`, new crate; pure permutation algebra — no chemistry, no geometry,
  no AST deps; a dependency of umol-graph). The dense coset index **is the OpenSMILES arrangement number**
  (`@TH/@AL/@SP/@TB/@OH`) — reproduced, never re-invented. **Built complete for all five geometries now** —
  it's a leaf crate, so finishing it here means later non-tetrahedral *perception* extends only the chemistry
  layers (C/D/E), never this one (B–E themselves still scope to TH/CT).
  - **A′1** — `Permutation` (one permutation, n ≤ 6, `Copy`, one-line notation): `identity`/`from_image`/
    `between(from, to)` (the relabel τ) / `apply`/`compose`/`inverse`/`sign`/`act`/`rank`/`unrank` (Lehmer —
    internal canonical ordering only); `Ord`+`Eq`+`Hash`.
  - **A′2** — `PermutationGroup` (subgroup of Sₙ, full enumeration): `generate(degree, gens)` (brute-force
    closure), `symmetric`/`alternating`/`cyclic`/`dihedral`, `order`/`contains`/`elements`.
  - **A′3** — `CosetSpace` = R (a `PermutationGroup`) + the coset **partition** (one canonical rep per coset):
    `count` (= n!/|R|), `coset_rep`. The algebraic layer — which orderings are equivalent, not yet a numbering.
  - **A′4** — `ClassKey` (`FromStr`/`Display` key, à la `SchoenfliesSymbol`): families `Sym/Alt/Cyc/Dih(u8)`
    (natural action) + geometry variants `Tetrahedral`/`CisTrans`/`SquarePlanar`/`TrigonalBipyramidal`/
    `Octahedral`. `static REGISTRY: LazyLock<Mutex<HashMap<ClassKey, &'static CosetSpace>>>`, `space(k)` =
    build-once + `Box::leak` (mirrors umol-msym `point_group.rs:71`); `Coset { space: &'static CosetSpace,
    index: u32 }` ties identity to the interned space (like `SymmetryOp`). Const generators per class:
    `Tetrahedral`→`alternating(4)` (A₄), `CisTrans`→double-swap (Z₂), `SquarePlanar`→`dihedral(4)` (D₄),
    `TrigonalBipyramidal`→D₃ on 5, `Octahedral`→O on 6 — all five built now.
  - **A′5** — the **index alignment** `CosetSpace::index(perm) -> u32` = the OpenSMILES decomposition for **all
    five** classes, transcribed from the spec and validated so its fibers equal R's coset partition (A′3):
    - **TH/CT** — the parity (the `sign` bit).
    - **SP** — the `U/4/Z` path shape (§3.8.5) → `@SP1–3`.
    - **TB** — apical-pair axis + 3-equatorial winding (§3.8.6) → `@TB1–20`.
    - **OH** — axis + 4-equatorial shape + winding (§3.8.7, recursive on the SP shapes) → `@OH1–30`.
    No ordering of our own — each is the spec's table, fiber-checked against R. Exhaustive round-trip tests
    (`index ∘ unindex = id` over all `n!` permutations) per class.
  - **A′6** — `reindex(k, input_index, τ) -> u32` (= `index(τ ∘ unindex(input_index))`) + the `~`/`^k` operator
    action on indices (over Sₙ/core(R); only `~`/`^1` in scope, Phase G). raise (B3/B5) and perception (C)
    call `space(k).reindex(…)`.
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
- **Phase D — DSL round-trip** (`umol-ast/src/dsl/stereo.rs`, new; wired into `dsl/molecule.rs`, mirror
  `aromatic`). Faithful round-trip, no frame conversion (config stored relative to the written `:ligands`
  order; `#T`/`#C` relative to the local frame).
  - **D1** — config-string parser/writer: `class config` ↔ `(StereoKind, StereoConfigurationAst)`. Head
    `Th`/`Ct`; config `* | ! | + | <coset-term>`, coset-term recursive (`nat`→`Lit`, `?id`→`Expr(Var)`,
    `~e`→`Expr(SwapOp)`, `e^k`→`Expr(ApplyOp)`); `Expr::Set`/`VarDomain` (`{…}`, `?o :: {…}`) reserved at the
    surface (deferred with non-tetrahedral). **One function** — D3's `:type` head and D5's `#T`/`#C` call it.
  - **D2** — ligand surface: `atom-ref | [:h atom-ref] | [:lp atom-ref]` (kind-first; reserved
    `[:bond/:port/:fragment ref]`) ↔ the ordered `StereoLigand` list. Unknown tags rejected (no silent pass).
  - **D3** — `:stereo-atoms` / `:stereo-bonds` entry reader/writer: `{ :id? (keyword), :site ref,
    :ligands [ ligand+ ], :type config-string }` ↔ a `StereoAtom`/`StereoBond` (focus + D2 ligands + D1 config
    + kind).
  - **D4** — sugar `:ccw`/`:cw`/`:e`/`:z` (each carries its class — `Th1`/`Th2`/`Ct1`/`Ct2`) ↔ the `:type` head.
  - **D5** — `#T`/`#C` atom/bond-string surface: the derived-predicate tokens in the existing atom/bond
    constraint-string parser (`dsl/constraint.rs`) — `#T<config>`/`#C<config>` (local-frame, the **same** D1
    `StereoConfigurationAst` parser) inside the atom-string (`C#h#T1`) / bond-string.
  - **D6** — round-trip tests: EDN↔AST for both surfaces (elements *and* `#T`/`#C` strings) over the
    ~150-file corpus, under `--features conformance`.
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
    the other operand; when `Set`/`VarDomain` land (non-TH), `meet` becomes a genuine domain intersection — the
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
| config value | `StereoIndexAst { Undetermined, Lit(u32), Expr(Box<Expr>) }` over a recursive stereo `Expr { Lit, Var, SwapOp, ApplyOp, Set, VarDomain }` (its own, ≠ `value::Expr`); `~`/`^k` recurse (`~1`, `0^1^1` sayable); `Undetermined` out of `Expr` ⇒ `~+` unrepresentable; `Lit` duplicated by design; `Set`/`VarDomain` deferred; index = `u32` dense coset index per class (SMILES arrangement number, not Lehmer) |
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
