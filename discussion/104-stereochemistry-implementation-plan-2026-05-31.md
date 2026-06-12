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
- **stereo symmetry + drivers** (chemistry): the symmetry computation is umol-ast (`StereoSymmetry`, a
  graph-symmetry algorithm — **not** a umol-graph `…Perception` engine); the umol-graph drivers mirror the aromaticity ops
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

- **Phase B — raise → `#T`/`#C` assertions** (mechanical; `table_ir/raise.rs`). **DONE (2026-06-09)** — implemented
  and validated (SMILES `@`/`/`,`\`; MOL parity + wedges; SMILES/MOL cross-checks; `WedgeConflict`/`DanglingBondDirection`/
  `CisTransConflict`). MOL coordinate-derived cis/trans is out of B by design (perception). raise transcribes the source
  model's **atom/bond-based** stereo into the umol AST's per-atom / per-bond `#T`/`#C` constraints, reindexed
  into the atom's / bond's own umol neighbor frame. SMILES `@`/`@@` and MOL wedges are **per-atom** assertions
  (→ `#T` on the atom); `/`,`\` are **per-bond** assertions (→ `#C` on the bond) — the same atom/bond-based
  chirality model the source formats use. raise writes **only** these per-atom / per-bond assertions; it
  constructs **no** `StereoAtomAst`/`StereoBondAst` element and never touches the `stereo_atoms`/`stereo_bonds`
  tables. Discrete stereo elements appear **only downstream, in perception (Phase C)**, which consumes `#T`/`#C`
  as input. This is the **aromaticity pattern, exactly**: raise asserts the per-atom `#a` (`AromaticValence`)
  flag (`raise.rs:87–99`) and only later does `AromaticResolver` build aromatic-system elements from it; here
  raise asserts `#T`/`#C` and only later does the stereo resolver (C) build stereo elements from them. `#T`/`#C`
  are the resolver's **input**, never its output.

  **Current understanding — input conventions (research 2026-06-06).** Sources: `materials/codes/{cdk,rdkit}`
  (the reference impls), OpenSMILES §3.8, and a 95,170-record V2000 corpus scan under `materials/`.

  - **umol target frame (reshapes B1).** The coset index of a raised `#T`/`#C` is numbered against the atom's
    umol neighbor frame computed **directly in raise from TableIR**: real neighbors in **atom-list (index) order**
    (identical between TableIR and the AST), then the implicit H, then lone pairs. It is **not** `stereo_ligands()`
    on a view — at raise time the atom/bond views carry no stereo ligands, so the frame is a topology + field read,
    not an element reconstruction. The "shared with C4 / C4 = B1" claim in B1 does not hold for raise.
  - **source frame (reshapes B2).** The arrangement is read in the **source format's** neighbor order, converted
    to a source-frame coset index, then reindexed (B3) into the umol frame:
    - **SMILES (symbolic, coordinate-free).** Tetrahedral `@` = TH1 = anticlockwise, `@@` = TH2 = clockwise,
      viewed **from the first-listed neighbor** toward the center (OpenSMILES §3.8.2; CDK Beam `@TH1`=anticlockwise,
      `@TH2`=clockwise). Neighbor order = the preceding (from-)atom, then the **implicit H immediately after the
      from-atom** (first, if the chiral atom opens the SMILES), then the remaining neighbors as written;
      ring-closure neighbors take the position of their bond digit. Cis/trans `/`,`\` are the directional single
      bonds flanking the double bond, relative to the nearer sp² carbon: `F/C=C/F` trans, `F/C=C\F` cis (§3.8.3).
    - **MOL V2000 atom parity (symbolic, coordinate-free — the primary symbolic path).** The CTfile "ignored when
      read" note is itself ignored by CDK and RDKit; parity **is** read. CDK `MDLV2000Reader.createStereo0d`:
      carriers = neighbors in **connection (bond) order**; with 3 explicit neighbors the focus atom is appended
      last as the 4th carrier (placeholder for the implicit H / lone pair "at the back"); an explicit H is treated
      as "always at the back" — invert winding when it lands at carrier slot 0 or 2. Parity **1 = clockwise,
      2 = anticlockwise**, 3 = unspecified. Enabled by default ("Allow stereo created from parity value when no
      coordinates"). RDKit has the equivalent parity→chirality path.
    - **MOL V2000 wedge and double-bond E/Z (in CDK/RDKit, coordinate-based; umol's scope differs — below).** In **both** CDK and RDKit the
      wedge (1 Up / 6 Down, narrow end at the center) is resolved to a configuration only **with the neighbors'
      2D coordinates** (`StereoElementFactory.using2DCoordinates`; RDKit `assignChiralTypesFromBondDirs` over a
      conformer). Double-bond cis/trans likewise comes from the **2D coordinates** (`using2DCoordinates`;
      `MolOps::detectBondStereochemistry`); the molfile double-bond stereo field carries only `3` =
      either/crossed = unknown. The RDKit bond-stereo field map (`MolFileParser.cpp:1769–1786`) is a flat switch,
      **not** order-conditioned: `0`→none, `1`→BEGINWEDGE, `6`→BEGINDASH, `3`→EITHERDOUBLE+STEREOANY, `4`→UNKNOWN
      — so `1`/`6` are wedge directions, never cis/trans.
  - **empirical (95,170 V2000 records).** Double bonds: stereo `0` 518,674 (coords) · `3` 1,607 (either) ·
    `1`/`6` 3 each (noise). Single bonds: `1` 53,683 (up) · `6` 35,462 (down) · `4` 1,272 (either single). MOL
    double-bond E/Z is coordinate-derived in practice; `io::ctfile::parser::convert::convert_bond_stereo_direction_code`'s
    `1→Cis` / `6→Trans` is both reference-unbacked and effectively dead — those codes are single-bond wedge
    directions, and on double bonds they occur only as noise.
  - **`Chirality` is source-relative — verified, and the disambiguator must be explicit.** Parsing a SMILES and a
    MOL form of the same molecule (`(R)-CFClBrI`, L-alanine, (S)-butan-2-ol, (±)-menthol; ChemDraw + PubChem) shows
    each parser is individually faithful, but the same `Chirality` variant means different things:
    - SMILES sets `Atom.chirality` from `@`/`@@`. umol now maps **`@` → `CounterClockwise` (CCW), `@@` → `Clockwise`
      (CW)** to match the OpenSMILES winding (the prior mapping was inverted). `source_format = SMILES` for a
      non-empty parse; empty/whitespace input stays `UNKNOWN`.
    - MOL sets `Atom.chirality` from the **parity field** (`1 → Clockwise`, `2 → CounterClockwise`), in atom-number
      order with the highest/implicit-H neighbor behind. Common editors (ChemDraw) write **parity 0 + wedges + 2D
      coords**, so `Atom.chirality` is **absent** and the configuration lives in `Bond.wedge` + the depiction;
      PubChem/OEChem write parity, so it is present.
    - Thus `Chirality::Clockwise` denotes **different reference frames and opposite winding** across the two sources;
      the frame and the implicit-H placement are *implied by the convention*, not stored.
  - **decision — explicit `ChiralityFrame` flag (implemented 2026-06-06).** `table_ir::stereo::ChiralityFrame`
    `{ FirstNeighborToward, LastNeighborAway }`, stored as `Option<ChiralityFrame>` on `Molecule` /
    `ExtendedMolecule`. It names *what the rule does* — the viewing reference for reading a per-atom chirality
    descriptor — not who defined it; it governs **tetrahedral atom chirality only**, not E/Z bonds. `FirstNeighborToward`
    = SMILES (`@` viewed from the first neighbor); `LastNeighborAway` = CTAB parity (highest/implicit-H neighbor
    behind). raise reads **this** (not `source_format`) to map (`Chirality` token, atom-list neighbor order, the
    frame's H-placement + viewpoint + winding) → coset. Chosen over parser-level normalization (reorders neighbors,
    discards raw fidelity) and over `source_format`-implicit dispatch (undocumented coupling).
    - **Set with `source_format`, not with the scope flag.** A parser sets `chirality_frame` unconditionally for a
      parsed (non-empty) molecule — SMILES/SMARTS → `FirstNeighborToward` in `parse_smiles_inner` /
      `parse_extended_smiles_inner`; MOL → `LastNeighborAway` in `build_molecule` / `build_extended_molecule`. This is
      a per-format constant, unlike `ConfigurationScope` (the renamed `StereoInterpretation`, `{ Absolute, Relative }`),
      which is content-derived (MOL chiral flag, CXSMILES `a:`/`r`) and stays `None` absent such a marker.
    - **The parser parses; it does not interpret.** The parser only records the raw descriptors
      (`Atom.chirality`, `Bond.wedge`, `Bond.stereo`, `positions`); converting them is raise's job, not the parser's.

  **What raise consumes (in scope) vs what is external.** raise converts **every symbolic descriptor** into a
  `#T`/`#C` assertion, after which the descriptor is **gone** — the AST carries only `#T`/`#C`, never a wedge or a
  parity code:
  - atom → `#T`: `Atom.chirality` (MOL parity / SMILES `@`/`@@`) **and** wedges (MOL up/down when no parity, read
    against the 2D depiction). umol treats a wedge + its 2D depiction as a symbolic descriptor it consumes here —
    unlike CDK/RDKit, which resolve wedges through a coordinate pass.
  - bond → `#C`: the directional `/`,`\` (`Bond.direction`).

  **Out of scope — bare 3D coordinates.** A 3D conformer with no symbolic marker is **not** interpreted; raise
  leaves `positions` intact and passes it through, and E/Z (or chirality) is derived from it by conversion
  functions **outside umol-graph**. The MOL double-bond field `0` ("from coordinates") falls here. Perception
  (Phase C) consumes `#T`/`#C` only — it never reads wedges or coordinates.

  All of B lives in `umol-io/src/table_ir/raise.rs` (private fns; may move to a `raise` util module later) and uses
  `umol-perm` — a **direct `umol-io` dependency**. The wedge helpers resolve the up/down wedge + 2D depiction into a
  handedness via the `umol-geometric-core` leaf crate (doc 106): `wedge_winding` lifts the wedged neighbor off the
  depiction plane (z = ±1, sign from up/down), places the missing substituent at `complementary_direction`, and
  takes the sign of `signed_volume` of the four points. This is synthetic geometry over the depiction (the lift
  magnitude is arbitrary — only its sign matters), distinct from `umol-geometric`, the heavy 3D crate reserved for
  Phase F.

  RESOLVED (2026-06-07): the SMILES input-ordering problem (`C[C@]1(Cl)CC(C)CC1` vs `C[C@](Cl)1CC(C)CC1` — same atoms,
  ring vs branch bond written in a different order, opposite stereo) is fixed in the parser, not raise. The SMILES
  builder now stores each ring-closure bond at its **opening** position (reserve a `bond_table` entry at the open
  digit, fill it at close), so `mol.bonds` order = OpenSMILES write order and `input_neighbor_ordering` slots the
  ring-closure neighbor at its digit. The two cases now raise to **different** `#T` cosets (`Lit(0)` vs `Lit(1)`),
  verified in raise. (CXSMILES bond indices, which count ring bonds at their *closing* digit, are remapped at the CX
  boundary — `cx::BondIndexMap` — so they keep resolving correctly against the open-order list.)

  - **B1 — umol target ordering** (`tetrahedral_target_ordering` / `cis_trans_target_ordering`). The order the coset
    is numbered against, built from TableIR: real neighbors in **atom-index order**, then `virtual_ligands`
    (implicit H ×h, lone pairs ×l, each **hosted by the atom** — the center is never a `kind: Atom` ligand); for #C,
    `atom_1`'s side then `atom_2`'s. A topology + field read — **not** an AST view, **not** `stereo_ligands` (none
    exist at raise; perception rebuilds the same order in C4).
  - **B2 — source ordering, dispatched on `molecule.chirality_frame`.** Each constructor yields the source ordering
    + a source index; the target frame (B1) is atom-index order, so the frames are **not** symmetric:
    - `first_neighbor_toward_ordering` (SMILES/SMARTS): real neighbors in `input_neighbor_ordering` — the center's
      incident bonds in **parse order** (the OpenSMILES write order; basic `Atom` records no explicit order, so
      raise reads it off `mol.bonds`), implicit H at slot 0 if `atom_idx == 0` else slot 1. Differs from the target
      (H last; ring-closure slotting), so it **reorders**. `@` = CCW = TH1, `@@` = CW = TH2.
    - `last_neighbor_away_ordering` (MOL): atom-number order = atom-index order, H/lone pair **last** — *identical*
      to the target (`== tetrahedral_target_ordering`), so it passes through (`between` = identity); only the
      winding converts (parity `1` = CW, `2` = CW2, last ligand behind, → source index).
    - `tetrahedral_wedge_ordering` / `cis_trans_wedge_ordering` (no descriptor): read the up/down wedge (or `/`,`\`)
      against the 2D depiction. Geometry-free, so they return the **target ordering** plus the source index from
      `wedge_winding` (local 2D-winding × up/down sign).
    Every source ordering reorders the *same set* the target holds (shared `virtual_ligands`), so `between` never
    mismatches (it `assert_eq!`s the lengths).
  - **B3 — reindex.** `space(k).reindex(source_idx, Permutation::between(&source, &target))` → `target_idx`
    (umol-perm, A′6): the permutation carries the ordering / H-placement difference (identity for MOL), the B2
    winding the source index. 2-coset ⇒ a parity XOR.
  - **B4 — `raise_tetrahedral_stereo` → `#T`** (`AtomConstraint::TetrahedralStereo`). Whichever atom source is
    present, consumed into `Stereo(Lit(target_idx))`: the **descriptor** (MOL parity / SMILES `@`,`@@`, via B2/B3),
    or the **wedge** when `Atom.chirality` is absent (ChemDraw, parity 0). `Unspecified` (`@?` / wavy) → `#T+`.
    Non-tetrahedral `Chirality` (Allenal / SquarePlanar / TrigonalBipyramidal / Octahedral) is out of scope → no
    constraint.
  - **B5 — `raise_cis_trans_stereo` → `#C`** (`BondConstraint::CisTransStereo`). The flanking `/`,`\` (`Bond.direction`)
    → cis/trans → `Stereo(Lit(target_idx))`. MOL double bond: field `3` → `#C+`; field `0` ("from coordinates") →
    **not raised** (external).
  - **B6 — undetermined / absent.** The `#T+` / `#C+` cases above (`@?`, `Bond.stereo == Either`); absent ⇒ no
    constraint. raise never emits `*`/`!`/`Var`/operators (pattern-side only).
  - **B7 — wiring, sign-pinning, validation.** Wire into the molecule-level raise — extend the `raise.rs` atom loop
    (after the aromaticity block) and add a bond loop, mirroring that block: push the constraint when `Some`. Pin
    the unknowns with **inline `#[rstest]`** in `raise.rs`, not by guessing: (i) ordering tests assert the exact
    `StereoLigand` vectors and fix `input_neighbor_ordering`'s ring-closure behavior; (ii) raise tests assert the
    `#T`/`#C` coset and fix the sign constants (`@`/`@@`→idx, parity→idx, `wedge_winding`→idx, cis/trans→idx) from
    known configurations (CFClBrI, L-alanine, (±)-menthol, an E/Z, a ring stereocenter). Then a molecule-level test
    that **SMILES and MOL of the same molecule raise to the same physical `#T`/`#C`**; the conformance suites
    (MOL + SMILES; no element output — those are C7); and a cross-check of a few cosets against RDKit/CDK.
  - **B7a — umol-perm: coset spaces over an explicit parent group (prerequisite for `#C`).** `#C`'s realizable
    arrangements are not `Sₙ`: the two substituents on each sp² carbon are bonded to that carbon, and the bond may
    be written with either carbon first, so the parent is `S₂ ≀ S₂ = D₄` — within-side swaps `(0 1)`, `(2 3)` **and**
    the carbon swap `(0 2)(1 3)` (order 8) — not `S₄`. The cis/trans descriptor is `D₄ / V`, quotient by the Klein
    four `V = D₂ = {e, (0 1)(2 3), (0 2)(1 3), (0 3)(1 2)}` (face flip, carbon swap, both) → `8/4 = 2` cosets:
    `{e, (0 1)(2 3)}`-class = cis, `{(0 1),(2 3)}`-class = trans. The carbon swap living in `R = V` is what makes
    `#C` invariant to which carbon is written first (the source→target relabeling can itself be a carbon swap, so it
    must be inside the parent). `D₄` here is a *different* embedding than the existing `dihedral(4)` (the cyclic-square
    Sylow-2 used for `SquarePlanar`); it's the partition-respecting Sylow-2, built by `generate`. `CosetSpace`
    currently hardwires the parent to `Sₙ` (`coset.rs`: `count = n!/|R|`; `CanonicalRank` enumerates all `n!` via
    `Permutation::unrank`), so it cannot express a sub-`Sₙ` parent. Generalize `Sₙ/R → P/R`:
    - `CosetSpace` gains a `parent: PermutationGroup` field; `new(parent, group, decomposition)` asserts
      `group ⊆ parent` (`parent.contains`); `count = parent.order() / group.order()`.
    - `Decomposition::CanonicalRank` enumerates `parent.elements()` (was: `unrank` over `0..n!`), maps each to its
      `coset_rep` (min over `R∘σ`), sorts + dedups. The bespoke `SquarePlanar` / `TrigonalBipyramidal` / `Octahedral`
      decompositions are unchanged (explicit reps, all in `Sₙ`).
    - `coset_rep` / `index` / `unindex` / `reindex` logic unchanged; `reindex` requires the relabeling `∈ parent`
      (assert — `#C`'s within-side and carbon-swap relabelings are all in `D₄`).
    - `ClassKey::build()` arms return `(parent, R, decomposition)`. Existing arms set `parent =
      PermutationGroup::symmetric(degree)` (Tetrahedral → `symmetric(4)`/`alternating(4)`, SquarePlanar →
      `symmetric(4)`/`dihedral(4)` = `S₄/D₄` = 3, TB → `symmetric(5)`/…, OH → `symmetric(6)`/…) — behavior identical.
      `CisTrans` → `parent = generate(4, &[(0 1),(2 3),(0 2)(1 3)])` (`D₄`), `R = generate(4, &[(0 1)(2 3),(0 2)(1 3)])`
      (`V`), `CanonicalRank` → 2 cosets, **degree 4** (no longer `alternating(2)`/degree 2). The `space()` registry is
      unchanged.
    - Tests: `CisTrans` count = 2 and the eight `D₄` elements map to the two cosets (cis/trans, including the
      carbon-swap images); `reindex` over within-side **and** carbon-swap relabelings gives the same coset; the
      `group ⊆ parent` assertion fires on a bad pairing; existing counts unchanged (Tetrahedral 2, SquarePlanar 3,
      TB 20, OH 30). This is what lets B8's `#C` go through `space(ClassKey::CisTrans).reindex` on the 4-tuple
      `[atom_1's two ligands, atom_2's two ligands]` — consistent with `#T` — instead of a hand-rolled binary; B8's
      `#C` (below) uses it.
  - **B8 — under-determined ligands: no inference at raise (2026-06-07; supersedes the implicit-H / lone-pair
    inference in B1/B4/B5 and the degree-2 `#C` ordering).** A missing tetrahedral position or one-substituent
    `#C` side may be an implicit H or a lone pair; raise does **not** decide which and does **not** write
    `atom.lone_pairs` / `implicit_hydrogens` — that is resolution's job. `#T`/`#C` assert only the coset, never what
    the ligands are, which is the point: an opaque virtual ligand stands in, placed last in the target frame.
    - **`#T`** (keeps `space(ClassKey::Tetrahedral).reindex`). Count neighbors via `neighbor_count`: 4 → use them;
      3 → one `StereoLigand::Virtual` (last in `tetrahedral_ligand_ordering`; at `ligand_idx` 0/1 in
      `first_neighbor_toward_ordering`; the `complementary_direction` point in `coset_from_wedge_winding`); `< 3` or `> 4` →
      `RaiseError::TetrahedralLigandCount`. The virtual count is `4 − neighbors`, **not** an H/LP inference, so
      sulfoxides/sulfonium etc. raise to the same coset as before — only the (wrongly set) `#n` annotation is gone.
    - **`#C`** — via `space(ClassKey::CisTrans)` over parent `D₄` / quotient `V` (B7a), on the 4-tuple
      `[atom_1's two ligands, atom_2's two ligands]` (each side padded to two with a `StereoLigand::Virtual` when it
      has one substituent). Per side, resolve each substituent's **halfplane** (`StereoHalfplane {Top, Bottom}`) from
      its `/`,`\` (`Bond.direction`) — toward the substituent, with the geminal substituent's mark inverted
      (disagreement → `RaiseError::CisTransConflict`). Build `source` from the halfplanes and
      `target = [atom_1's ligands by index, atom_2's by index]`; then
      `coset = space(ClassKey::CisTrans).index(Permutation::between(&source, &target))` → `0` = cis, `1` = trans
      (`F/C=C\F → 0`, `F/C=C/F → 1`; the `BondDirection`→halfplane orientation is pinned by tests). The `V` quotient
      (face flip + carbon swap) makes the result invariant to which carbon is written first and to within-side order.
    - **Capability gate (`#C`), 2026-06-09 — supersedes `CisTransLigandCount` / `has_cis_trans_marker`.** A double bond
      is raised only when **cis/trans-capable** (`cis_trans_capable`: each end has a substituent besides the other);
      otherwise `raise_cis_trans_stereo` returns `Ok(None)` (terminal `=O`/`=CH2`, plain ethylene — no error, no
      marker needed). A directional `/`,`\` flanking **no** capable double bond is orphaned →
      `RaiseError::DanglingBondDirection`, checked per-bond by `validate_bond_direction` in the bond loop (the cis/trans
      analog of a chirality token on an under-coordinated atom). A capable side with no `/`,`\` → not raised
      (coordinates / external).
    - **Errors** (all raise-time): `#T` `< 3` or `> 4` neighbors (`TetrahedralLigandCount`); inconsistent MOL wedges at
      a `#T` center (`WedgeConflict`); contradictory geminal `/`,`\` on one carbon (`CisTransConflict`); orphaned
      `/`,`\` (`DanglingBondDirection`).
    - **`StereoLigand`** is raise-local (`raise/utils.rs`) `{ Atom(usize), Virtual(usize) }` (the `Virtual` carries its
      host atom so the two carbons' virtuals in a `#C` 4-tuple stay distinct — `Permutation::between<T: Eq>` matches by
      equality); distinct from `umol_ast::ast::StereoLigand`.
    - **As-built names / sign (2026-06-09).** `/`,`\` live on `Bond.direction: BondDirection {Rising, Falling}`, split
      from `Bond.wedge: BondWedge` (MOL tetrahedral wedges only) so the dangling-marker check sees only cis/trans
      markers (doc 106). Other renames: `coset_from_wedge_winding`, `tetrahedral_ligand_ordering`,
      `validate_bond_direction`, `MismatchedRingBondDirections`. Wedge sign (B7, “pinned by tests”):
      `signed_volume < 0 ⇒ coset 0` (matching SMILES `@` = CCW = 0; an inverted sign was found and fixed).

  **Coset index.** The config value is the **dense coset index per stereo class** — the OpenSMILES arrangement
  number (`@TH1-2` … `@OH1-30` = n!/|R| cosets, `u32`; not a Lehmer rank). Input conventions are in *Current
  understanding* above (which supersedes this doc's earlier wedge/parity notes).

  **One invariant, two pinned constants, one dependency.** (1) Every source ordering must be a reordering of the
  **same ligand set** as `tetrahedral_target_ordering` — `Permutation::between` `assert_eq!`s the lengths and
  `expect`s each element, so it panics on a set/length mismatch — hence both share `virtual_ligands`. (2) Which
  `wedge_winding` sign (#T) and which cis/trans relation (#C) map to coset `0` follow `umol-perm`'s `Tetrahedral` /
  `CisTrans` numbering, fixed by B7 (OpenSMILES §3.8.3). (3) the SMILES branch needs the **OpenSMILES neighbor
  order** (`input_neighbor_ordering`) — the bonds incident to `atom_idx` in write order, ring-closure bond at its
  digit slot; basic `Atom` records no explicit order (only `ExtendedAtom.ligand_order`, CXSMILES), so raise derives
  it from incident-bond parse order — **confirmed (2026-06-07): the builder preserves it** by storing ring-closure
  bonds at their opening position (see RESOLVED note above). `b.other` / `b.start_atom` are
  `Bond` methods; `neighbors` / `virtual_ligands` are local primitives, and `wedge_winding` calls
  `umol-geometric-core` (`signed_volume`, `complementary_direction`).
  `atom_idx` is the atom's `usize` index (0..n_atoms−1); `AtomId` appears only on the AST `StereoLigand`.

  ```rust
  // ---- ligand orderings ----------------------------------------------------------------------
  // A virtual ligand (ImplicitHydrogen / LonePair) is hosted by its atom — the center for #T, the sp²
  // carbon for #C; the center/carbon is never a `kind: Atom` ligand. Every source ordering is a
  // reordering of the same set the target ordering holds (so Permutation::between never mismatches).

  // neighbors of `atom_idx`, ascending index (= atom-list / MDL atom-number order).
  fn neighbors(mol: &TableMolecule, atom_idx: usize) -> Vec<usize> {
      let mut n: Vec<usize> = mol.bonds.iter().filter_map(|b| b.other(atom_idx as u32)).map(|x| x as usize).collect();
      n.sort();
      n
  }

  // `atom_idx`'s virtual ligands: implicit H ×h then lone pairs ×l, each hosted by `atom_idx`.
  fn virtual_ligands(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
      let (atom, host) = (&mol.atoms[atom_idx], AtomId(atom_idx as u32));
      (0..atom.implicit_hydrogens.unwrap_or(0)).map(|_| StereoLigand::new(host, ImplicitHydrogen))
          .chain((0..atom.lone_pairs.unwrap_or(0)).map(|_| StereoLigand::new(host, LonePair)))
          .collect()
  }

  // #T target ordering (umol DSL semantics): real neighbors ascending, then the virtual ligands.
  fn tetrahedral_target_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
      let mut o: Vec<_> = neighbors(mol, atom_idx).iter().map(|&n| StereoLigand::new(AtomId(n as u32), Atom)).collect();
      o.extend(virtual_ligands(mol, atom_idx));
      o
  }

  // #T source ordering, FirstNeighborToward (SMILES/SMARTS): neighbors in OpenSMILES write order (the
  // from-atom first, a ring-closure neighbor at its bond-digit slot, branches as written); the implicit
  // H goes to slot 0 if `atom_idx` opened the SMILES (`atom_idx == 0`) else slot 1; remaining virtuals appended.
  fn first_neighbor_toward_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
      let mut o: Vec<_> = input_neighbor_ordering(mol, atom_idx).map(|n| StereoLigand::new(AtomId(n as u32), Atom)).collect();
      let mut virt = virtual_ligands(mol, atom_idx).into_iter();
      if mol.atoms[atom_idx].implicit_hydrogens.unwrap_or(0) > 0 {
          o.insert(if atom_idx > 0 { 1 } else { 0 }, virt.next().unwrap());   // the implicit H, at its written slot
      }
      o.extend(virt);                                                          // any further H / lone pairs
      o
  }

  // #T source ordering, LastNeighborAway (MDL / MOL parity).
  // Semantics: number the neighbors by atom number (= ascending index); an implicit H or lone pair is
  // the highest-numbered, so it comes LAST. The parity token is read with that last ligand pointing
  // AWAY (behind the plane of the other three): parity 1 (odd) = the first three, in number order, run
  // clockwise; parity 2 (even) = counterclockwise. The center is never a ligand; the trailing ligand is
  // the H/LP hosted by `atom_idx`. This is exactly the target frame — umol numbers against atom-index
  // order, which equals MDL atom-number order, so the ordering passes through unchanged; only the
  // parity→source_idx winding differs, and the reindex below is the identity.
  fn last_neighbor_away_ordering(mol: &TableMolecule, atom_idx: usize) -> Vec<StereoLigand> {
      tetrahedral_target_ordering(mol, atom_idx)
  }

  // #T source ordering from the wedge: one up/down wedge at `atom_idx` + the 2D depiction. The result is
  // frame-free, so the ordering returned is the target frame and `source_idx` carries the config, from
  // `wedge_winding` — lift the wedged neighbor off the depiction plane (z = ±1 by up/down), place the missing
  // substituent at `complementary_direction`, sign of `signed_volume` of the four points. → source_idx (0/1, pinned by B7). None if no wedge.
  fn tetrahedral_wedge_ordering(mol: &TableMolecule, atom_idx: usize) -> Option<(Vec<StereoLigand>, usize)> {
      let pos = mol.positions.as_ref()?;
      let (w, up) = mol.bonds.iter().find_map(|b| match (b.other(atom_idx as u32), b.wedge) {
          (Some(w), Some(BondWedge::Up))   => Some((w as usize, true)),   // narrow end at `atom_idx`
          (Some(w), Some(BondWedge::Down)) => Some((w as usize, false)),
          _ => None,
      })?;
      let order = tetrahedral_target_ordering(mol, atom_idx);
      let source_idx = wedge_winding(pos, atom_idx, &order, w, up);
      Some((order, source_idx))
  }

  // #C target ordering (umol DSL): atom_1's substituents (asc index, then its virtuals) ++ atom_2's (same).
  fn cis_trans_target_ordering(mol: &TableMolecule, atom_1: usize, atom_2: usize) -> Vec<StereoLigand> {
      let side = |c: usize, other: usize| {
          let mut o: Vec<_> = neighbors(mol, c).into_iter().filter(|&n| n != other)
              .map(|n| StereoLigand::new(AtomId(n as u32), Atom)).collect();
          o.extend(virtual_ligands(mol, c));
          o
      };
      let mut o = side(atom_1, atom_2);
      o.extend(side(atom_2, atom_1));
      o
  }

  // SUPERSEDED by B8's `cis_trans_side` (see "#C cis/trans — verified design" below). This sketch reads
  // the *wedged* substituent's face, which is wrong when a side has two real substituents: for
  // F/C(C)=C(Cl)/C it compares F and the `/`-marked methyl (anti ⇒ 1), but the stored atom-index coset is
  // 0 (compare the first atom-index ligands F and Cl, both down ⇒ syn). The correct rule reads the FIRST
  // atom-index ligand's face (inferring from the geminal marker when that ligand is unmarked).
  //
  // #C source ordering from the directional `/`,`\` flanking bonds. For each sp² carbon read its
  // directional substituent's side ±1: `/`=Up, `\`=Down on the x–c bond, flipped if `c` is the bond's
  // start so the slope reads toward `c`. Geometry is frame-free, so the ordering returned is the target
  // frame and `source_idx` carries cis/trans: syn (d1 == d2) ⇒ 0 (cis), anti ⇒ 1 (trans) — pinned to
  // OpenSMILES §3.8.3 (F/C=C/F trans, F/C=C\F cis) by B7. None if either side lacks a directional bond
  // (⇒ MOL field 0 ⇒ coords ⇒ external).
  fn cis_trans_wedge_ordering(mol: &TableMolecule, atom_1: usize, atom_2: usize) -> Option<(Vec<StereoLigand>, usize)> {
      let side = |c: usize, other: usize| -> Option<i8> {
          mol.bonds.iter().find_map(|s| {
              if s.order != BondOrder::Single { return None; }
              let x = s.other(c as u32)? as usize;
              if x == other { return None; }
              let dir = match s.wedge? { BondWedge::Up => 1, BondWedge::Down => -1, _ => return None };
              Some(if s.start_atom() as usize == c { dir } else { -dir })
          })
      };
      let (d1, d2) = (side(atom_1, atom_2)?, side(atom_2, atom_1)?);
      let source_idx: usize = if d1 == d2 { 0 } else { 1 };
      Some((cis_trans_target_ordering(mol, atom_1, atom_2), source_idx))
  }

  // ---- raise: build a (source, target) ordering pair and reindex ------------------------------

  // B4 — tetrahedral atom → #T.
  fn raise_tetrahedral_stereo(mol: &TableMolecule, atom_idx: usize) -> Option<AtomConstraint> {
      if mol.atoms[atom_idx].chirality == Some(Chirality::Unspecified) {          // B6 — @? / wavy
          return Some(AtomConstraint::TetrahedralStereo(Stereo(StereoCosetAst::Undetermined)));
      }
      let (source, source_idx): (Vec<StereoLigand>, usize) = match (mol.atoms[atom_idx].chirality, mol.chirality_frame) {
          (Some(tok), Some(FirstNeighborToward)) => {
              let source_idx = match tok {                                         // @ = TH1, @@ = TH2
                  CounterClockwise | Tetrahedral { arr: 1 } => 0,
                  Clockwise        | Tetrahedral { arr: 2 } => 1,
                  _ => return None,                                                // AL/SP/TB/OH — out of scope
              };
              (first_neighbor_toward_ordering(mol, atom_idx), source_idx)
          }
          (Some(tok), Some(LastNeighborAway)) => {
              let source_idx = match tok { Clockwise => 0, CounterClockwise => 1, _ => return None };  // parity 1/2
              (last_neighbor_away_ordering(mol, atom_idx), source_idx)
          }
          (None, _)       => tetrahedral_wedge_ordering(mol, atom_idx)?,           // no descriptor ⇒ wedge, else no #T
          (Some(_), None) => return None,                                          // descriptor without a frame — unreachable
      };
      let target_idx = space(ClassKey::Tetrahedral)
          .reindex(source_idx, Permutation::between(&source, &tetrahedral_target_ordering(mol, atom_idx)));
      Some(AtomConstraint::TetrahedralStereo(Stereo(StereoCosetAst::Lit(target_idx))))
  }

  // B5 — double bond → #C.
  fn raise_cis_trans_stereo(mol: &TableMolecule, bond_idx: usize) -> Option<BondConstraint> {
      let bond = &mol.bonds[bond_idx];
      if bond.order != BondOrder::Double { return None; }
      if bond.stereo == Some(BondStereo::Either) {                                 // B6 — MOL "either"/crossed
          return Some(BondConstraint::CisTransStereo(Stereo(StereoCosetAst::Undetermined)));
      }
      let (atom_1, atom_2) = (bond.start_atom() as usize, bond.end_atom() as usize);
      let (source, source_idx) = cis_trans_wedge_ordering(mol, atom_1, atom_2)?;   // None ⇒ field 0 ⇒ external
      let target_idx = space(ClassKey::CisTrans)
          .reindex(source_idx, Permutation::between(&source, &cis_trans_target_ordering(mol, atom_1, atom_2)));
      Some(BondConstraint::CisTransStereo(Stereo(StereoCosetAst::Lit(target_idx))))
  }
  ```

  **B8 pseudocode** (revises the orderings + raise above; no H/LP inference, no `atom.lone_pairs` write).

  ```rust
  // Raise-local opaque ligand; Permutation::between is generic over Eq, so no StereoLigandKind change.
  enum TableLigand { Atom(usize), Virtual }

  // #T target ordering: neighbors ascending, then ONE Virtual when the center has 3 neighbors (last).
  fn tetrahedral_target_ordering(mol, atom_idx) -> Vec<TableLigand> {
      let mut order: Vec<_> =
          atom_ordering(mol, atom_idx).iter().map(|&n| TableLigand::Atom(n)).collect();
      if order.len() == 3 { order.push(TableLigand::Virtual(atom_idx)); }
      order
  }

  // #T source ordering (FirstNeighborToward): Virtual at ligand_idx (0 if atom_idx == 0 else 1).
  fn first_neighbor_toward_ordering(mol, atom_idx) -> Vec<TableLigand> {
      let mut order: Vec<_> =
          bond_neighbor_ordering(mol, atom_idx).iter().map(|&n| TableLigand::Atom(n)).collect();
      if order.len() == 3 {
          order.insert(if atom_idx > 0 { 1 } else { 0 }, TableLigand::Virtual);
      }
      order
  }
  // last_neighbor_away_ordering = tetrahedral_target_ordering (Virtual last);
  // wedge path uses tetrahedral_target_ordering, and wedge_winding maps TableLigand::Virtual to the
  // complementary_direction point.

  fn raise_tetrahedral_stereo(mol, atom_idx) -> Result<Option<AtomConstraint>, RaiseError> {
      // no chirality token and no wedge bond at atom_idx -> Ok(None)
      let neighbor_count = atom_ordering(mol, atom_idx).len();
      if neighbor_count < 3 || neighbor_count > 4 {
          return Err(RaiseError::TetrahedralLigandCount { atom: atom_idx, count: neighbor_count });
      }
      // Unspecified -> #T+; else build (source, source_idx) per chirality_frame / wedge as before:
      let target = tetrahedral_target_ordering(mol, atom_idx);
      let coset = space(ClassKey::Tetrahedral)
          .reindex(source_idx, Permutation::between(&source, &target));
      Ok(Some(/* TetrahedralStereo Stereo(Lit(coset)) */))
      // no (constraint, lone_pairs) return; no atom.lone_pairs / implicit_hydrogens write.
  }

  // #C via the coset machinery (B7a: space(CisTrans) over parent D₄, quotient V). Build each side as
  // [up-face, down-face] from the /,\, plus its atom-index pair, then take the coset of the relabeling.
  fn raise_cis_trans_stereo(mol, bond_idx) -> Result<Option<BondConstraint>, RaiseError> {
      let bond = &mol.bonds[bond_idx];
      if bond.order != BondOrder::Double { return Ok(None); }
      if bond.stereo == Some(BondStereo::Either) { return Ok(Some(/* #C Stereo(Undetermined) */)); }
      if !has_cis_trans_marker(mol, bond.start_atom(), bond.end_atom()) { return Ok(None); } // plain double bond
      let (Some(side_1), Some(side_2)) = (
          cis_trans_side(mol, bond.start_atom(), bond.end_atom())?,
          cis_trans_side(mol, bond.end_atom(), bond.start_atom())?,
      ) else {
          return Ok(None); // a side carries no /,\ -> coordinates / external
      };
      let source = [side_1.up_ligand, side_1.down_ligand, side_2.up_ligand, side_2.down_ligand];
      let target = [side_1.first_ligand, side_1.second_ligand, side_2.first_ligand, side_2.second_ligand]; // atom index
      let coset = space(ClassKey::CisTrans).index(Permutation::between(&source, &target));
      Ok(Some(/* CisTransStereo Stereo(Lit(coset)) */))
  }

  // One sp² atom's two #C ligands: by atom index (first_ligand, second_ligand|Virtual) and by face
  // (up_ligand, down_ligand). first = lowest-index substituent != other_atom_idx (atoms precede virtuals).
  // None if the side carries no /,\ (-> external).
  struct StereoBondAtom { first_ligand: StereoLigand, second_ligand: StereoLigand,
                          up_ligand: StereoLigand, down_ligand: StereoLigand }
  fn cis_trans_side(mol, atom_idx, other_atom_idx) -> Result<Option<StereoBondAtom>, RaiseError> {
      let subs: Vec<usize> =
          atom_ordering(mol, atom_idx).into_iter().filter(|&n| n != other_atom_idx).collect(); // ascending
      let first = *subs.first().ok_or(RaiseError::CisTransLigandCount { atom: atom_idx })?;
      let first_ligand = StereoLigand::Atom(first);
      let second_ligand = subs.get(1).map(|&o| StereoLigand::Atom(o)).unwrap_or(StereoLigand::Virtual(atom_idx));
      // face of `first`: from the bond toward it and the inverted bond toward the second; conflict -> Err
      let toward_first  = direction(mol, atom_idx, first);
      let toward_second = subs.get(1).and_then(|&o| direction(mol, atom_idx, o)).map(|up| !up);
      let up = match (toward_first, toward_second) {
          (Some(a), Some(b)) if a != b => return Err(RaiseError::CisTransConflict { atom: atom_idx }),
          (Some(face), _) | (_, Some(face)) => face,
          (None, None) => return Ok(None),
      };
      Ok(Some(if up { StereoBondAtom { first_ligand, second_ligand, up_ligand: first_ligand, down_ligand: second_ligand } }
              else   { StereoBondAtom { first_ligand, second_ligand, up_ligand: second_ligand, down_ligand: first_ligand } }))
  }
  // direction(mol, atom_idx, other): the /,\ (BondWedge) on the single bond atom_idx-other as a face
  // (true = up), oriented toward `atom_idx` (negate when `atom_idx` is the bond's stored end). Relies on the
  // builder's ring-closure flip (below) so a `\1`/`/1` written at the closing atom is already in the bond's
  // stored orientation. (#T MOL wedges: multiple wedge bonds at atom_idx must agree else WedgeConflict.)
  ```

  **#C cis/trans — verified design (2026-06-08).** Implemented in `umol-io/src/table_ir/raise.rs`;
  validation done by feeding reference SMILES and comparing `cargo run -p umol-graph --bin stereo`. The two
  templates (`X`,`Y` the sp² atoms; `0`,`1` on `X`, `2`,`3` on `Y`):

  ```
   0   2          0   3
    \ /            \ /
     X=Y   Z=0      X=Y   E=1
    / \            / \
   1   3          1   2
  ```
  "same side / opposite sides" is the `D₂`-invariant mnemonic for the two cosets (like CW/CCW for `#T`).
  The only freedom is how `0,1,2,3` map to substituents — CIP uses atomic-number rank, SMILES uses parse
  (bond) order, **umol uses atom-index order**. umol stores the coset in the **atom-index (target) frame**.

  - **Target ordering = atom index** (unchanged; `cis_trans_side` uses `atom_ordering`). Each side is its
    substituents ascending by index, then virtuals; an atom always gets a lower index than a virtual.
  - **Reduces to one bit.** `between([up,down,up,down], [first,second,first,second])`'s `D₄/V` coset is
    exactly: the **first** (atom-index) ligand on each side has the **same face → 0, opposite → 1**. The
    second ligand is the geminal opposite and never matters. The 4-field `StereoBondAtom` + `between` is
    kept only for symmetry with `#T`; `#C` could be a direct `(face_1 != face_2) as u32`.
  - **Face of `first`** = `direction(atom_idx, first)`, or `!direction(atom_idx, second)` when only the
    second is marked; the two disagreeing ⇒ `CisTransConflict`.
  - **Ring-closure wedge flip (parser `on_ring_bond`, both impls).** A ring bond is stored opening-atom
    first, but a `/`,`\` written at the *closing* atom is from that atom's perspective, so the builder
    **flips** it — `final_wedge = open.wedge` if the opening end set it, else `wedge_opt.flip()` —
    exactly as it already does for donation. Both ends set it ⇒ consistent only when the raw symbols are
    **opposite**; equal raw symbols ⇒ `MismatchedRingBondDirs` (mirrors the donation conflict rule).
    Without the flip `direction()` reads a ring substituent's face inverted. `BondWedge::flip` swaps
    Up↔Down, EitherUp↔EitherDown.

  **Worked example — `C/C=C1CO\1` ((E)-2-ethylideneoxirane); exercises both fixes.** Atoms
  `C0`(Me) `C1` `C2`(opens ring 1) `C3`(CH₂) `O4`(closes ring 1); double bond `C1=C2`, ring bond stored `[2 4]`.
  - In **SMILES (source) order**: left first = Me, right first = O (the ring digit precedes `C3`). Me down,
    O up (the `\1` flip) → opposite → **source coset 1**.
  - SMILES→umol reindex: indices follow string order, so O4 (idx 4) is *after* C3 (idx 3) though the bond to
    O is written first. Source `[O,CH₂]` vs target `[CH₂,O]` = one swap (odd) ⇒ coset flips ⇒ **stored
    atom-index `#C` = 0**. umol emits `2#C0`; test `ethylideneoxirane → Lit(0)`.
  - This is the only case where bond order ≠ atom index; for ring-free molecules they coincide, which is why
    the simpler cases didn't separate the two frames.

  **Validation methodology.** Reference cosets are given in SMILES/bond *source* order; the stored coset is
  that value reindexed by the parity of the source→target permutation (even = unchanged, odd = flipped).

- **Phase C — perception → stereo elements** (chemistry; `umol-graph/src/ops/stereo*`). Lifts `#T`/`#C` +
  topology into the molecule-level stereo overlay and provides the stereo query/operation surface — modeled on the
  aromaticity ops (`AromaticityModel` + `AromaticityPerception` + resolver/validator/transformer). The storage
  overlay is built: `stereo_atoms`/`stereo_bonds` (`FixedVarBirelationSet`), `StereoAtomId`/`StereoBondId`,
  builders, remap/restore, and `stereo_ligands()` (real neighbors in adjacency order, then implicit-H ×h, then
  lone pairs — materializes the virtual ligands). Same-frame lift: the element's ligand order is `#T`'s incidence
  frame (equivariant, not canonical), so config copies through with no reindex.

  - **Common core (the `AromaticityModel`/`AromaticityPerception` analog).**
    - `StereoModel` — policy (sibling of `AromaticityModel`, `ops/model.rs`): kinds in scope (TH/CT now),
      stereogenicity strictness (lean — emit nothing for a nonstereogenic site), para-stereocenter refinement.
    - `StereoSymmetry` — the shared engine that resolution and the queries consume, built on a stereo-aware
      color-refinement / automorphism pass: extend `umol-graph-core/algorithms/auto.rs` to fold the config
      descriptor into the color and re-refine to a fixpoint (the para-stereo case).
  - **Algorithms.**
    - **WL** (Weisfeiler–Lehman color refinement) — the refinement core: fast approximate hash + automorphism
      classes; the substrate for stereogenicity, canonical labeling, and isomorphism seeding.
    - **VF2 / VF2++** (long planned) — exact graph isomorphism **carrying the per-site stereo descriptor** in the
      match, seeded by the WL hash; VF2++ is the target node-ordering refinement.
  - **Resolution** (`StereoResolver`, ≈ `AromaticResolver`). Marker-driven: scan `#T`/`#C`, build the
    `StereoAtomAst`/`StereoBondAst` element (focus, ligand set, kind, config — same incidence frame, no reindex);
    plus the inverse projection element→`#T` for write-back. Lean: drops nonstereogenic sites (via the topicity
    test below). Tests: raise→resolve round-trip over the corpus, asserting focus/ligands/config for known R/S, E/Z.
  - **Validation** (≈ `validator/aromaticity`; pass/fail, never mutates). Three checks: (1) **coset range** — config
    index within the kind's coset count (TH/CT 0–1; `Th3`/`Ct3` over-range fails here); (2) **ligand count vs kind
    degree** — element ligand count matches the kind (TH/CT = 4, incl. virtuals); (3) **derived ↔ asserted
    projection** — the `#T`/`#C` on the atom/bond equals the projection re-derived from the resolved element. A
    stereo element on a **non-stereogenic** site is **not** a contradiction (the element is agnostic about ligand
    distinctness — it covers prochirality), so validation stays silent on it.
  - **Topicity / stereogenicity** — ligand distinctness is the **automorphism-orbit partition** (homotopic /
    enantiotopic / diastereotopic), computed by the WL/automorphism core (`StereoSymmetry`), **not** CIP. CIP is a
    separate derived module (priority total-order → R/S/E/Z) that may consume this partition — the dependency runs
    CIP → core, never core → CIP. Resolution's lean drop and the non-stereogenic handling below key off it.
  - **Handling non-stereogenic elements** (open). Validation is pass/fail and never alters structure, so a redundant
    stereo element on a non-stereogenic site is surfaced by one of: a **strip transformer** (mutates — removes them,
    a cleanup op) and/or a **linter** — a new op category that emits advisory lints, no fail, no mutate
    (clippy-style). Not mutually exclusive (lint to flag, transform to fix); split to decide. Prochirality ×
    reactions interactions are **out of scope** for now.
  - **Perception** (geometry → stereo). 3D conformer or 2D+wedges → `#T`/`#C`/elements — the coordinate path raise
    refuses (perception, not translation; MOL cis/trans and marker-free tetrahedral fall here). Needs the
    `MoleculeAst`+positions container; the inverse — config → wedges/2D/3D (depiction/embedding) — is the
    aromatizer/kekulizer shape. (Container + perceiver/depicter detailed in doc 103, Future plan.)
  - **Operations** (over `StereoSymmetry` + the `umol-perm` coset algebra).
    - **remove stereo** — strip elements / `#T`/`#C`.
    - **invert / manipulate** — `~` (enantiomer / other config) and `^` (group action), already in `umol-perm`.
    - **enumerate stereoisomers** — expand undetermined (`+`) configs over each site's coset space.
    - **canonical labeling** — stereo-aware canonical order (WL fixpoint; the storage frame stays equivariant).
    - **isomorphism incl. stereochemistry** — VF2(++) equality up-to-isomorphism with the descriptors.

  Minimum: `StereoModel` + the refinement engine + resolver + validator (graph-only). Perception (geometry) waits
  on the positions container; the heavier derived layers (CIP, symmetry numbers, meso) and reaction stereo are in
  doc 103, Future plan.
  - **C1/C2/C3 (done)** — storage, builders, remap/restore.
  - **C9** — ligand-count / ligand-order queries, separate from topology and from constraints. Needs design.

  **Section C design summary (2026-06-10).** Two layers, not to be conflated:
  - **`StereoModel` = perception *configuration*** (per-*kind* policy): which kinds are active; per
    kind, which **chemical elements + bond orders** participate in perception (KindScope); per kind, a
    **fluxionality on/off toggle**; the para-stereo fixpoint toggle. It holds no per-molecule data.
  - **Stereo-element *data*** (per element, stored, molecule-specific): the **arrangement record**
    (ordered ligands + coset) plus two per-element constraints — **`LigandSymmetry`** (static ligand
    distinctness) and **`Fluxionality`** (dynamic interchange).

  The deliverable scopes the kinds to TH/CT but the shape must admit non-tetrahedral centers, allenes,
  atropisomers, and planar/helical chirality (background: docs 049/101/102/103/107; StereoMolGraph).

  - **KindScope** — active kinds (each a `umol-perm` `ClassKey`) + assignment. Per-kind
    **element scope** (Th: C, N, S…; Sp: Cu, Pd, Pt…) **and bond-order / topology scope** (double
    bond → planar E/Z; cumulated → allene; hindered single pivot → axial). Perception assigns
    site → kind by scope + coordination/order match (mirrors `AromaticityModel` / `ElementScope`).
  - **Fluxionality (A2b — lazy quotient).** Each kind has a parent `P` and rigid proper group `R`;
    the stored config is the canonical `P/R` coset (unchanged by fluxionality). Adjoining the
    rearrangement moves gives `R ⊆ R′ ⊆ P`; observable configs are `P/R′`, with the merge table the
    surjection `P/R ↠ P/R′` (the `R′`-orbit partition of the base cosets). Rigid = `R′=R` (TH/CT);
    fully fluxional = `R′=P` (Berry pseudorotation → 1); partial = intermediate `R′`. The stored
    index never changes; the merge is applied on demand by equality / stereogenicity / enumeration.
    Fluxionality is written as a set of **additional generators only** (image notation, `⊆ P`)
    adjoined to `R`: `R′ = generate(R ∪ gens)`. Usually a singleton — `R` conjugates one move into
    the whole process (biaryl free rotation = one within-side swap `2134` → `R′=D₄` → 1 config;
    one Berry move + `D₃` → the full system). Named presets (FreeRotation / Berry / RingFlip) are
    optional sugar over the generators. The generators are **stored per stereo element** (a
    `Fluxionality` constraint — molecule-specific data, an external assertion); `StereoModel` carries
    only a per-*kind* on/off toggle (allow/exclude fluxionality during perception), **not** the
    generators and **not** a per-(chemical-element × kind) setting. Rejected **A1** (fluxionality as
    distinct kinds, à la StereoMolGraph `HinderedBond12/13/23/33`): conflates geometry with dynamics —
    a continuous change in hindrance or temperature would force a discrete kind swap.
  - **Topicity / prochirality.** Orbit relations under `Â` (proper automorphisms) vs `Â*`
    (proper + improper): homotopic / enantiotopic / diastereotopic; prochiral iff an enantiotopic
    ligand pair. **Derived** (not a model drop), recorded/asserted via the per-element `LigandSymmetry`
    constraint (below). The model owns only the para-stereo fixpoint toggle (opt-in).
  - **Kind anatomy.** kind = shared `CosetSpace` (proper part) + inversion/parity + site type +
    scope + fluxionality default. Double-bond E/Z, allene, and biaryl axial **share** the 2+2 `D₂`
    proper coset space and differ **only** in the inversion generator (StereoMolGraph `PlanarBond`
    inversion = None / achiral vs `AtropBond` inversion = pair-swap / chiral) — which is what `Â*`
    reads and what makes `~` a diastereomer-swap vs an enantiomer-swap.
  - **Engine** (settled). Optional WL fast first pass → nauty/VF2 exact, with the stereo descriptor
    folded into the node color → union-find the orbits (= StereoMolGraph
    `atom_automorphism_classes`). Para-stereo = an outer recolor-and-refine fixpoint, opt-in.
    Cordella VF2 is already in `umol-graph-core`; VF2++ node-ordering is the later optimization.
  - **Color derivation** (separate module, `umol-graph`, configurable). Selects which resolved
    atom/bond features fold into an `Ord + Copy` color (xxh3 `u64`) feeding `auto.rs` / `refine.rs`:
    element, isotope, charge, implicit-H, lone pairs, spin, aromatic / ring membership, localized
    valence; bond order / aromatic / dative. Not stereo-specific (reused by canonical labeling,
    isomorphism, symmetry numbers); the stereo overlay augments it with the coset descriptor.
  - **Prochirality / stereogenicity model (resolved 2026-06-10).**
    - **The element is an arrangement record**, not a stereogenicity claim: `:stereo-atoms` /
      `:stereo-bonds` = ordered ligands (index / atom-ref) + coset. Ligand distinctness is **not
      required and not asserted**. The element coset is `StereoCosetAst` = {`Th0`, `Th1`, `Th*`}
      (binary kinds), `Th*` = `Undetermined` (unknown *which* coset). There is **no `NotStereo` and
      no "could-be-NotStereo-or-Stereo" `Undetermined` at the element** — those are
      `StereoConfigurationAst`, specific to `#T`/`#C`, which need them because they are predicates on
      a *ligand-free* atom/bond where "is it even a stereocenter" is live. `#T`/`#C` have **implicit**
      ligands (topology-defined) ⇒ at most one virtual ligand; prochiral arrangements (multiple
      distinct virtual slots, e.g. 1,1-dichloroethene) are expressible only at the molecule-level
      element. (This is the existing D2 split in `stereo.rs`.)
    - **Labeled vs molecular level.** A concrete coset is a faithful, coordinate-testable fact about
      the *labeled* (Born–Oppenheimer-distinguishable) ligands — `Th0` on a CH₃ records that H₁,H₂,H₃
      are in a definite CW/CCW arrangement; it is **never vacuous**. **Stereogenicity** is "does a
      labeled-arrangement distinction survive the molecular automorphism quotient?" — a *derived,
      orthogonal* predicate, not a property of the coset. So the **lean/eager fork dissolves**: there
      is no policy about which sites store elements; you store an arrangement wherever you want to
      record/assert one, and stereogenicity is always derived.
    - **Static distinctness vs dynamic interchange — two axes, one machinery.** Ligand distinctness
      ("are the ligands identical?") and fluxionality ("can they interconvert by T-dependent
      dynamics?") both reduce cosets via `merge_under` (a permutation is a permutation), but are
      distinct chemistry concepts, asserted independently:
      - **`LigandSymmetry` (static, derived, per element).** The invariant ligand permutations Π,
        graded proper/improper — the primitive from which distinctness / stereogenicity / topicity all
        follow. **Derived** from the *global* molecular automorphism group (intrinsic; remote chirality
        can make isomorphic ligand subgraphs diastereotopic — the para-stereo case) and
        **cross-checkable** like valence / `#v`. `stereogenic(element)` ⟺ the stored coset is a
        singleton class of `merge_under(Π_proper)`; proper Π → homotopic + stereogenicity, improper Π →
        enantiotopic.
      - **`Fluxionality` (dynamic, external, per stereo element).** The rearrangement generators — an
        **external** assertion (barrier / temperature), **not derivable** and **not cross-checkable**,
        proper-only — **stored as a per-element `Fluxionality` constraint**. `StereoModel` carries only
        a per-*kind* on/off toggle; a per-(chemical-element × kind) setting is rejected as too
        fine-grained.
      They **compose**: same-molecule (static) = `merge_under(Π)`; same-observable-species (at T) =
      `merge_under(Π ∪ fluxional)`.
    - **Constraint shape — signed-permutation literals over the group Π (revised 2026-06-10).** The
      scalable primitive is the **group** Π (a subgroup of the signed / permutation-inversion group),
      **not** its orbit partition. A ground molecule carries the **full Π** (derived from
      automorphisms); **stereogenicity for any kind** (non-binary included) is the coset-merge of the
      full Π via `merge_under` — exact, no projection. The graded **equivalence relation (orbits) is a
      derived, legible view**, never the stored primitive — it is lossy (two Π with the same orbits
      merge SP/TB/OH cosets differently, e.g. `⟨(0 1)(2 3)(4 5)⟩` vs `⟨(0 1),(2 3),(4 5)⟩` on OH).
      The constraint is a `(B, N)` predicate of **signed-permutation literals**: positive `g ∈ Π`
      (symmetry present → ligands equivalent) and negative `g ∉ Π` (broken → distinct), over
      **arbitrary** elements (not just transpositions — `(0 1 2) ∈ Π ∧ (0 1) ∉ Π` pins C₃ vs S₃).
      Two-directional with **no complement**: a negative is an explicit literal tested against the
      concrete Π, not a (non-existent, non-Boolean) complement subgroup. The `∈`/`∉` asymmetry is
      structural — membership is closed (positives generate `B = ⟨…⟩`), non-membership is not
      (negatives stay a forbidden set `N`). Equivariant (relabeled with the element on remap).
      Supersedes the earlier positive-only `{ proper, improper }` generator-list shape.
    - **Lattice on the literal sets.** `meet` = **union** of literals (exact narrowing — the
      matching / resolution workhorse; `None` if a required element is forbidden,
      `⟨B₁,B₂⟩ ∩ (N₁∪N₂) ≠ ∅`, or unrealizable); `join` = **intersection** of literals (the LUB,
      over-approximating the concrete set-union — fine, disjunctive requirements don't arise);
      `matches(ground Π)` = every positive `∈ Π` and every negative `∉ Π`. CIP labels (pro-R/pro-S,
      Re/Si faces) deferred with CIP; Re/Si is a separate `trigonal` kind (doc 103), not this constraint.
    - **Constraint notation (settled 2026-06-10; `#o`/`#g` added 2026-06-11).** `#p` (LigandSymmetry)
      and `#f` (Fluxionality) are **non-unique** constraints — multiple per element, implicit conjunction (the `#R` convention);
      `{…}` stays reserved for disjunction and is **not** used here. Each entry is one signed literal,
      written as a cycle-notation permutation with two optional prefixes, `[!][']perm` (the cycle
      parens self-delimit — no wrapper):
      - `'` marks **improper** (the element's sign — binds the permutation); a leading **`!`** negates
        the literal (`∉ Π` — a propositional negation, outermost, matching `!`'s general negation role).
        Four atomic forms: `(…)` = `g∈Π` proper, `'(…)` = `g∈Π` improper, `!(…)` = `g∉Π` proper,
        `!'(…)` = `g∉Π` improper.
      - **`perm` = GAP-style product-of-cycles** — comma-separated, **0-indexed** (matching ligand
        positions; no DSL↔internal offset): `(0,1,2)(3,4)`, identity `()`. Prefixes bind the **whole**
        product (`'(0,1)(2,3)` is one improper element). `(0,1)(2,3)` = the double-transposition (e.g.
        CT's C₂), *not* the same as separate `(0,1)`, `(2,3)` (those generate the larger
        `⟨(0 1),(2 3)⟩`). (`()` = identity = vacuous; `!()` = always-false.)
      - `#f` uses the same form, proper + positive (`#f(…)`); `'` / `!` available but unusual.
      - **`^` migrates to the same GAP-style cycle notation** (replacing the one-line image `^2134`), so
        permutations read uniformly across `^`, `#p`, `#f`, and umol-perm. **0-indexed** throughout (the
        `^` migration drops its 1-indexed image for 0-indexed cycles). Ripple: `StereoExpr::ApplyOp`
        parse/print (D1); cycle **construct/decompose** on `Permutation` in umol-perm —
        `from_cycles(degree, cycles)`, `cycles() -> Vec<Vec<usize>>` (disjoint-cycle decomposition),
        `Display` (GAP cycle notation, 0-indexed) — the **string parse stays in the DSL** (it supplies
        the degree); and docs/tests using image notation.
    - **Topicity / stereogenicity as assertable lossy constraints (settled 2026-06-11).** Besides the
      precise `#p` group literals, the **lossy** topicity/stereogenicity facts are *also* assertable —
      the ring-count / ring-membership analog (lossy yet matchable; `#p` is the precise
      substructure-analog). The earlier "query-only" framing is withdrawn. Both ride a **unified
      same/different glyph language**, read at pair-level (`#o`) or site-level (`#g`):

      | glyph | `#o` topicity (pair) | `#g` stereogenicity (site) | meaning |
      | --- | --- | --- | --- |
      | `=` same      | homotopic      | symmetric   | equivalent under proper |
      | `'` mirror    | enantiotopic   | prochiral   | equivalent only under improper |
      | `/` different | diastereotopic | stereogenic | inequivalent under both |

      The mapping is **aligned, not inverted** — homotopic ligands ⟹ symmetric site, enantiotopic ⟹
      prochiral, diastereotopic ⟹ stereogenic — the same glyph at two granularities.
      - **`#o` topicity** — a relation between two ligands, reusing cycle notation **restricted to a
        transposition**: `#o=(i,j)`, `#o'(i,j)`, `#o/(i,j)` (explicit glyph, no bare form). A **distinct
        tag** from `#p` because `#p(i,j)` asserts the *transposition ∈ Π* (element-level) whereas topicity
        asserts the *orbit relation*, and the two differ for non-binary kinds (C₃: `0,1` homotopic yet
        `(0 1)∉Π`).
      - **`#g` stereogenicity** — a site-level flag: `#g=` / `#g'` / `#g/`. Derived: `#g/` ⟺ stored
        coset is a singleton class of `merge_under(Π_proper)`; `#g'` ⟺ not that but ∃ an enantiotopic
        ligand pair; `#g=` otherwise.
      - **`!` is a true negation** (not a value selector): because `/` carries "different" positively,
        `!` is freed for real negation, restoring the full **subset lattice** over `{=, ', /}` for both
        tags — `*` = any, singletons, `!X` = complement (`#g!/` = non-stereogenic = `{=, '}`; `#o!=` =
        not-homotopic = `{', /}`), sets = disjunction, `∅` = `None`. Value lattice: `meet`=∩, `join`=∪,
        `matches` = derived ground ∈ asserted set.
      - **Collection lattice.** `#g` is **unique per site** (one flag; `meet` = value-meet, `None` on
        conflict). `#o` is **keyed-unique per (unordered) pair** — a `{i,j} → value` map (the
        `JointDomain` keying), distinct pairs conjoined; collection `meet` = per-pair value-`meet`
        (`None` if a pair is double-asserted with disjoint sets, e.g. `#o=(0,1)` ∧ `#o/(0,1)`), union
        over pairs; `join` = per-pair value-`join`, keep pairs in both; `matches` = every asserted pair's
        value contains the target's derived topicity. CH₂Cl₂ = `{#o=(0,1), #o=(2,3), #o/(0,2)}` — three
        keys, conjoined.
      - **Stereo-bond carriers** use the **2-per-terminus, 4-position** ligand model (`[0,1]` terminus 1,
        `[2,3]` terminus 2); same-terminus pairs are the stereo-determining ones, and `#p`/`#o`/`#g`/`#f`
        apply unchanged over those positions (`#p(0,1)(2,3)` = the axial C₂; `#f(0,1)` = single-terminus
        swap merging cis/trans).
      - **cis/trans (achiral kinds).** `'` (enantiotopic / prochiral) is **meaningless** — no
        enantiomers — so `#o'`/`#g'` on a CT (or any achiral) kind is a **validator error**. The
        meaningful CT values are `=` and `/` taken w.r.t. the `~` involution (the 1,2-swap);
        1-chloroethene's "substitute at the CH₂ → diastereomers" is `#o/` on that terminus. The
        three-member nomenclature is kept; the validator gates `'` out for achiral kinds. EDN-reserved
        `'`/`/` are not a problem — the short-hand DSL is quoted.
      - **`~` sugar** — `~` = the class's canonical involution permutation (doc-103 `~`; binary classes
        → `(0,1)`, class-parameterized for TB/OH). Usable as the cycle operand everywhere: `#p~` ≡
        `#p(0,1)`, `#p'~` ≡ `#p'(0,1)`; `#o=~` ≡ `#o=(0,1)` (homotopic), `#o/~` ≡ `#o/(0,1)`
        (diastereotopic). The *same* involution the coset operator `~` uses — one applies it to a config
        value, the other names it as a literal. `~` resolves **eagerly** to a concrete `Permutation`: a
        stereo element's kind is always known and literal (a non-literal kind is an invalid element), so
        the class involution is well-defined at parse/lift — no deferred permutation expression is needed.
      - **EDN representation (structured; 2026-06-11).** The DSL uses **keywords, not EDN reader tags**
        (spec §2 / §7.1: `#` appears only *inside* atom/bond/stereo strings, never as a reader-dispatch
        tag). A **permutation** is a **vector of cycles**, each cycle a 0-indexed `[int…]` vector; identity
        `[]`: `[[0 1 2] [3 4]]` = `(0 1 2)(3 4)`. One vector-of-cycles is **one element** (vs two
        generators = two values); canonicalized least-element-first per cycle, cycles sorted, fixed points
        dropped, for deterministic roundtrip. `PermutationAst` hosts its `FromEdn`/`ToEdn`. This is the
        shared primitive reused by the molecule-level symmetry assertions ([[110-molecular-symmetry-structure]]).
        The constraint envelope inhabits the **already-provisioned** `StereoAtomConstraint` /
        `StereoBondConstraint` (empty today by state, not by prohibition — the type, `Lattice`, collection,
        and DSL plumbing are in place): a new `:stereo-atom` / `:stereo-bond` entity-constraint key carrying
        `{:ligand-symmetry [{:perm [[0 1]] :orientation :improper :member :not-in} …]}`,
        `{:fluxionality [[[0 1]] …]}`, `{:topicity {[0 1] :homotopic, [0 2] :diastereotopic}}`,
        `{:stereogenicity :stereogenic}` — `:orientation` (`:proper` default) / `:member` (`:in` default)
        omittable; relation values are the AST-enum keywords (`:homotopic`/`:enantiotopic`/`:diastereotopic`,
        `:symmetric`/`:prochiral`/`:stereogenic`), a **set** `#{…}` for `!`-complements, `:undetermined` for
        any. **Two serializations, mirroring atoms** (`dsl/atom.rs` + `dsl/constraint.rs`): the per-element
        constraints serialize **inline** in the §7.14 stereo-string (`#p`/`#f`/`#o`/`#g`), with
        `StereoAtomDsl` round-tripping via `Edn::Str(to_string())` — `dsl/stereo.rs`, the atom.rs pattern;
        the structured `:stereo-atom` / `:stereo-bond` entity-constraint key (above) is the molecule-scope
        peer in `dsl/constraint.rs` (parallel to `:atom` / `atom-constraint-form`), with
        `lift_constraints`/`inline_constraints` moving between them. This adds the `:stereo-atom`
        entity-constraint key + §7.14 stereo-string predicates to the spec (the planned evolution).
    - **No coherence machinery needed.** A concrete coset on a non-distinct site is legitimate
      labeled data, not an error to validate against. Molecular equality / stereoisomer enumeration
      uses the standard equivariant-storage → canonicalize-both-sides quotient umol already commits
      to for atom numbering (doc 103 D2) — not a stereo-special obligation. (StereoMolGraph's
      `parity ∈ {None,0,±1}` fuses which-coset + orientation-unknown + chiral/achiral that this model
      keeps orthogonal — `None` ≡ `Th*`, `0`-vs-`±1` ≡ the derived `improper`; it also bakes a
      partial molecular quotient into the descriptor, mixing the labeled and molecular levels.)
  - **Validator wiring** (settled). The stereo validator runs **last** in the composite validator
    chain, after aromaticity. A stereo linter is out of scope.
  - **Chirality & `umol-perm` (D).** Direct `umol-graph → umol-perm` dependency. Chirality is a
    derived property of `(CosetSpace, improper generator)`, **moved into `umol-perm`**: each
    `ClassKey` carries an **improper (orientation-reversing) generator** (StereoMolGraph's
    `inversion`); `is_chiral` ⇔ it maps some coset to a *different* one, and `enantiomer(coset)`,
    the self-enantiomeric / meso test, and `Â*` all derive from it. `StereoKind::is_chiral_class`
    delegates here instead of hardcoding. The improper generator is **distinct from `~`**: `~` is
    "the other config" (= the improper op for chiral kinds, but a deliberate non-trivial coset swap
    for achiral kinds, where the improper op is a no-op on cosets — intentional). The improper op is
    exposed in the coset AST as `StereoExpr::MirrorOp`, glyph `'` (prime) — an involution (`''=id`;
    a no-op for achiral kinds, where a config is its own mirror); `~` stays distinct, its chiral
    instances delegating to `improper`.
  - **Axial kind (naming resolved).** Add `ClassKey::Axial` + `StereoKind::Axial` (allene /
    cumulene / biaryl / spiro) — IUPAC-aligned ("axial chirality" / "chirality axis"; descriptors
    M/P, older Rₐ/Sₐ, deferred with CIP). It **shares CisTrans's proper coset space** (`D₄`-parent,
    `V`-quotient, 2 cosets, same `~`) and differs **only** in the improper generator (CisTrans
    trivial-on-cosets / achiral; Axial coset-swap / chiral). Factor the proper-space build, register
    two keys. `CisTrans` retained for double-bond E/Z; 1:1 `StereoKind ↔ ClassKey` preserved.
  - **Prerequisite API (settled 2026-06-10).**
    - `umol-perm` (group-theoretic names only): `CosetSpace::merge_under(&self, extra_gens) ->
      Vec<u32>` — the merge table (each base index → its canonical representative under
      `⟨R, extra_gens⟩`); **serves both** fluxionality (rearrangement gens) **and** stereogenicity (the
      ligand stabilizer Π) — a permutation is a permutation. Asserts `extra_gens ⊆ parent`; reuses
      `coset_rep` / `index` / `unindex`.
      `PermutationGroup::extend(&self, extra) -> Self` builds `⟨R, extra_gens⟩`. The `improper`
      generator on `CosetSpace` (4th arm of `ClassKey::build`, transcribed + test-pinned; identity
      for achiral CT/SP, coset-swap for Axial) with `is_chiral()` / `enantiomer(index)` derived via
      `reindex`. `ClassKey::Axial`. Fluxional generators are **proper** (`⊆ parent`); `improper` is
      the separate chirality datum — they never mix.
    - `umol-ast`: `StereoExpr::MirrorOp` (glyph `'`, involution; folds via `enantiomer`, no-op for
      achiral); `StereoKind::Axial`; `is_chiral_class` delegates to `umol-perm`; `~` retained, its
      chiral instances delegating to `improper`. Per-element constraints (filling the today-empty
      `StereoAtomConstraint` / `StereoBondConstraint`): `LigandSymmetry` (static distinctness) and
      `Fluxionality` (dynamic interchange, proper-only), each a set of **signed-permutation literals**
      (`g ∈ Π` / `g ∉ Π`) over the group, `meet` = literal-union, `join` = literal-intersection,
      `matches` = test against the concrete Π; equivariant.
    - `umol-graph`: `StereoModel` = perception config (active kinds; per-kind chemical-element +
      bond-order scope; per-kind fluxionality on/off; para-stereo toggle). The observable-coset lookup
      takes a stereo element's stored `Fluxionality` generators (gated by the kind's toggle) → calls
      `merge_under`; molecular-distinctness uses the element's `LigandSymmetry` Π → the same
      `merge_under`.
  - **Implementation plan — umol-ast topicity chunk (2026-06-10).** Build order: umol-perm primitive →
    color seam → perception loop → queries → constraints → threading.
    1. **umol-perm — `OrientedPermutation`** = `Permutation` + a Z₂ `Orientation { Proper, Improper }`
       (the **permutation-inversion** group `Sₙ×Z₂` — *not* a hyperoctahedral signed permutation; and
       `Orientation`, *not* `Sign`, which is the existing parity `Permutation::sign()`). `compose`
       multiplies orientations. `OrientedPermutationGroup` is Π; its *storage* (full group vs
       `proper: PermutationGroup` + `improper_rep`) is deferred.
    2. **`MoleculeColoring` seam (umol-ast)** — one method over the public `Entity` ref (`ast/molecule.rs`,
       generally useful), **not** views: `trait MoleculeColoring { fn color(&self, &MoleculeAst, Entity)
       -> u64; }`. **One impl now**, `ConstitutionColoring` (inherent-field scheme): the color is a hash of
       the entity's **inherent fields only**. This is *complete* for our consumer (a graph automorphism):
       every **derived** predicate is relational, hence a function of the structure and preserved for free;
       once each overlay is an **entity pseudonode** (loop, next), even overlay-derived predicates
       (aromaticity, dative/multicenter participation) are graph functions. Named `MoleculeColoring` (not
       `Invariants`/`Coloring`) to avoid conflation with graph-core's proper-coloring. nauty's refinement
       subsumes Morgan for perception; ECFP/FCFP/scaffold are **separate impls** owning their own (often
       derived) fields — their consumer is per-atom hashing, where "relational ⇒ free" does not hold —
       reusing this seam. Stereo folding is added by the perception loop (next), keeping the trait
       stereo-free and fingerprint-reusable. (Full fingerprint module is its own future effort; only this
       interface is built now.)
    3. **`stereo_symmetry` loop (umol-ast)** — `MoleculeAst::stereo_symmetry(&self,
       coloring: &impl MoleculeColoring, para_stereo: bool) -> StereoSymmetry`. Builds the graph (atoms + a
       pseudonode per relation/overlay, colored by `coloring.color(self, entity)`; dative direction
       gadget-encoded). Loop: `node color ⊕ canonical_stereo_descriptor(rel. to current partition)`;
       `auto = automorphisms(color)` (graded
       Â/Â*, extending `AtomAutomorphism`); return when the partition is stable **or** `!para`.
       `para=false` ⇒ one pass; `para=true` ⇒ InChI-style fixpoint — full structure, para a pure
       on/off. `StereoSymmetry` holds the converged molecule-wide graded automorphism. Config
       carrier: a small umol-ast `StereoSymmetryConfig { coloring, para_stereo }` (room for a max-iterations
       guard); `StereoModel` (umol-graph) composes it + the kind-scope and passes it down — umol-ast
       keeps **no** dependency on umol-graph.
    4. **Queries — read-only `StereoAtomView`/`StereoBondView`, taking `&StereoSymmetry`** (the consumer
       computes it and passes it to its own view queries; no interior mutability). The **validator** is the
       only pipeline consumer (computes it once, molecule-wide, on the resolved AST); the **resolver is
       structural** and computes none; query-time reuse caches on the immutable `Molecule` (doc 086). Nothing
       is shared op-to-op — and couldn't be, since the resolver mutates the molecule. `ligand_symmetry(&p) -> Π`
       (site-stabilizer, lift atom-perm → ligand-position perm with virtual-ligand blocks, grade
       proper/improper); `is_stereogenic(&p)` = stored coset is a singleton class of
       `merge_under(Π_proper)`; `is_chiral()` (kind-level, umol-perm, no perception); `topicity(a,b,&p)
       -> Topicity`; `is_homotopic`/`is_enantiotopic`/`is_diastereotopic(a,b,&p)`; `is_prochiral(&p)`.
    5. **Constraints — fill the empty `StereoAtomConstraint`/`StereoBondConstraint`.** Four variants, each
       wrapping a dedicated AST struct: `LigandSymmetry(LigandSymmetryAst)`, `Fluxionality(FluxionalityAst)`,
       `Topicity(TopicityAst)`, `Stereogenicity(StereogenicityAst)`, where
       `LigandSymmetryAst { perm: OrientedPermutationAst, mem: MemOp }` (one `±` literal over Π — `mem`
       `In`/`NotIn` = `g∈Π` / `g∉Π`, the `'` = `Orientation::Improper`); `FluxionalityAst { perm:
       PermutationAst }` (positive proper move); `TopicityAst { pair: LigandPairAst, rel:
       TopicityRelationAst }` (the `#o` pair relation; `LigandPairAst { first: StereoLigandId, second:
       StereoLigandId }` = an unordered pair of ligand positions, normalizing `first ≤ second` like `AtomPair`);
       `StereogenicityAst(StereogenicityRelationAst)` (the `#g` site flag — just the relation, no
       pair). The relation value (the **subset lattice over `{=, ', /}`**: `{ Undetermined, Lit, LitSet,
       NotSet }` — any / singleton / positive set / `!`-complement; `#o!=` lowers to `NotSet({Homotopic})`,
       **no** separate `Not`) is generated by a **declarative macro** over each ground enum
       (`Topicity { Homotopic, Enantiotopic, Diastereotopic }`, `Stereogenicity { Symmetric, Prochiral,
       Stereogenic }`) — **no** `StereoRelationKind` trait (not worth it for two; the macro bakes the domain
       for the full-set→`Undetermined` collapse). `meet`=∩, `join`=∪, `matches`=∈; the glyph↔variant map
       (`=`/`'`//') lives in the **DSL only**. `PermutationAst(Permutation)` and `OrientedPermutationAst { perm:
       PermutationAst, orientation: Orientation }` are **thin wrappers** over the umol-perm types — no added
       structure (a permutation literal is always concrete: no `Undetermined`, no per-permutation lattice),
       kept so the AST layer stays uniform and can host the EDN/DSL traits (orphan rule); plus the
       **local** `MemOp { In, NotIn }` (a duplicate — the `value.rs` `MemOp` is scoped to element/isotope
       sets; keep the stereo module decoupled).
       Lattice: `LigandSymmetry`/`Fluxionality` — `meet` = literal-union (`None` on a required-vs-forbidden
       clash or unrealizable), `join` = literal-intersection (over-approximating LUB), `matches(&p)` =
       derive Π and test each literal (frame-aligned via `permutation_for`). `Topicity` — keyed-unique per
       (unordered) pair: collection `meet` = per-pair value-`meet` (set ∩, `None` on disjoint) + union
       over pairs, `join` = per-pair value-`join` (set ∪) keeping shared pairs, `matches` = derived
       topicity ∈ asserted set per pair. `Stereogenicity` — unique per site, value-`meet`/`join` = set
       ∩/∪, `matches` = derived ∈ set. `remap` = no-op (position-indexed; `Ordered` preserves order). The
       validator gates `'` out for achiral kinds. Storage shape (flat literal `Vec` vs `(B, N)`) deferred.
    6. **Threading + wiring.** Representation-independent now: `from_parts` / view carry of the
       per-element constraints, umol-graph `StereoModel` (`StereoSymmetryConfig` + kind-scope) and
       `StereoValidator` last in the composite chain. Inhabiting the constraint enums (`match` arms,
       `Lattice`, DSL) follows the **step-5 representation/storage (TBD; notation settled)**.
  - **Artifacts (concrete surface — steps 1–4 and 6; step 5's representation/storage TBD, notation settled — see "Constraint notation").**
    - **1 · umol-perm `oriented.rs` (new):** `enum Orientation { Proper, Improper }` (`compose` / `flip`
      / `is_proper`); `struct OrientedPermutation { perm: Permutation, orientation: Orientation }`
      (`new` / `proper` / `improper` / `identity`; `perm`, `orientation`, `is_proper`, `degree`,
      `apply`, `compose`, `inverse`); `struct OrientedPermutationGroup` (`generate`, `order`,
      `contains`, `elements`, `proper -> PermutationGroup`, `improper_rep -> Option<OrientedPermutation>`,
      `proper_orbit_of`, `star_orbit_of`). Ext `permutation.rs`: `Permutation::from_cycles(degree,
      cycles)`, `cycles() -> Vec<Vec<usize>>` (disjoint-cycle decomposition), `Display` (GAP cycle
      notation, 0-indexed) — the construct/decompose algebra for `^` / `#p` / `#f`; the cycle-string
      parse lives in the umol-ast DSL (it supplies the degree).
    - **2 · umol-ast `ast/coloring.rs` (new) + `Entity` in `ast/molecule.rs`:** public `enum Entity` (a typed
      ref to any molecule entity: Atom/Bond/Dative/Aromatic/Multicenter/Noncovalent/StereoAtom/StereoBond);
      `trait MoleculeColoring { color(&MoleculeAst, Entity) -> u64 }` (one method, not views); `struct
      ConstitutionColoring { features }` (the one impl — hashes each entity kind's **inherent fields** + a
      kind tag; derived predicates excluded as automorphism-free); `struct ColorFeatures` (bitflags, inherent
      only — element / isotope / charge / implicit-h / lone-pairs / spin; bond order / charge / spin).
    - **3 · umol-ast `ast/incidence.rs` (new, generic) + `ast/stereo_symmetry.rs` (new):** generic
      `struct IncidenceGraph { graph: Graph, node_entity: Vec<Entity> }` (the molecule's incidence/Levi
      graph — relations lifted to CSR pseudonodes; dative direction in the gadget edges) + `MoleculeAst::
      incidence_graph()`. Then `struct StereoSymmetryConfig<C: MoleculeColoring> { coloring, para_stereo,
      max_iterations }`; `struct StereoSymmetry` (the *result*: `proper_orbit_of`, `star_orbit_of`,
      `same_proper_orbit`, `same_star_orbit`, `graded`, `ligand_symmetry(site, &[StereoLigand])`);
      `MoleculeAst::stereo_symmetry<C>(&self, &cfg) -> StereoSymmetry` (runs the automorphism over
      `incidence_graph` with stereo-folded colors). Ext `ast/automorphism.rs` (graded Â/Â* + expose
      generators); ext `umol-graph-core/algorithms/auto.rs` (`Automorphism::generators()`).
    - **4 · umol-ast `ast/stereo.rs` + `ast/views/stereo.rs` (ext):** `enum Topicity { Homotopic,
      Enantiotopic, Diastereotopic }` and `enum Stereogenicity { Symmetric, Prochiral, Stereogenic }`
      (the derived ground triples; the `=` / `'` / `/` glyph map is DSL-only). On `StereoAtomView` /
      `StereoBondView` —
      `ligand_symmetry(&p) -> OrientedPermutationGroup`, `stereogenicity(&p) -> Stereogenicity`,
      `is_stereogenic(&p)`, `is_chiral()`, `topicity(a: StereoLigandId, b: StereoLigandId, &p) -> Topicity`,
      `is_homotopic` / `is_enantiotopic` / `is_diastereotopic(a, b, &p)`, `is_prochiral(&p)`; helpers
      `ligand_position(AtomId) -> Option<StereoLigandId>`, `ligand_frame(&p) -> Vec<StereoLigand>`,
      `induced_ligand_permutation(atom_perm) -> OrientedPermutation`. New id `StereoLigandId(u8)` in
      `ast/ids.rs` — a frame-relative ligand position (Permutation-width, not `usize`); existing
      position-taking view methods migrate to it.
    - **5 · umol-ast `ast/constraint/stereo.rs` + `dsl/stereo.rs`:** `PermutationAst(Permutation)` (thin
      wrapper, cycle-parsed), `OrientedPermutationAst { perm: PermutationAst, orientation: Orientation }`,
      local `enum MemOp { In, NotIn }`, `struct LigandPairAst { first: StereoLigandId, second: StereoLigandId }`
      (normalizing `first ≤ second`, like `AtomPair`); a **declarative macro** `relation_ast! { TopicityRelationAst,
      Topicity }` / `relation_ast! { StereogenicityRelationAst, Stereogenicity }` generating the 4-variant
      subset lattice (`Undetermined`/`Lit`/`LitSet`/`NotSet`, `meet`/`join`/`matches`) — no trait; the
      constraint structs `LigandSymmetryAst { perm: OrientedPermutationAst, mem: MemOp }`,
      `FluxionalityAst { perm: PermutationAst }`, `TopicityAst { pair: LigandPairAst, rel:
      TopicityRelationAst }`, `StereogenicityAst(StereogenicityRelationAst)`; the glyph↔variant table lives
      in `dsl/stereo.rs`. **Trait coverage** (codebase convention): the macro gives `*RelationAst`
      `Lattice` + `AsLit` (`Lit` = the ground enum); `StereogenicityAst` delegates both; `TopicityAst`
      impls `Lattice` keyed per pair (à la `JointDomainAst`); the **concrete literals** (`PermutationAst`,
      `OrientedPermutationAst`, `LigandPairAst`, `LigandSymmetryAst`, `FluxionalityAst`) can't be lattices
      (no top), so they carry `Eq`/`Hash` **plus a standalone inherent `matches(&self, &Self)`** (the
      `StereoExpr::matches_value` precedent — matchable without `Lattice`; the collection's matching
      dispatches per kind). No `Matches` trait is extracted — `matches` is already a directly-written
      (non-`meet`-delegating) method by the `Lattice` doc-contract, so it stays co-located with `meet` for
      the lattice types. The `StereoAtomConstraint` enum impls neither `Lattice` nor `AsLit` (the collection
      dispatches), matching `AtomConstraint`. `StereoAtomConstraint::{ LigandSymmetry(LigandSymmetryAst),
      Fluxionality(FluxionalityAst), Topicity(TopicityAst), Stereogenicity(StereogenicityAst) }`
      (bond parallels). `#p`/`#f` non-unique (`is_unique() = false`); `#o` keyed-unique per pair; `#g`
      unique per site. `StereoAtomConstraints` mirrors the **full `AtomConstraints` surface** — the stored
      side of the derived↔stored pair: per-kind stored accessors (`ligand_symmetry()`/`fluxionality()` iter,
      `topicity(pair)`, `stereogenicity()`), `get`/`add`(insert policy)/`contains`/`get_all`/`remove`/`iter`/
      `simplify_each`/`remap` (no-op — positions are frame-relative), and `Lattice { is_ground, meet, join,
      matches }` (`matches` is plain AST↔AST, like `AtomConstraints`; `#p`/`#f` union+dedup, `#o` per-pair
      value-meet/join with clash-`None`, `#g` value-meet/join). **No `matches_derived`/`matches_group`** —
      the derived↔stored cross-check is the **C4 validator's** job (parallels the atom validator's
      topology-derived-vs-stored valence cross-check), comparing `view.{ligand_symmetry,topicity,
      stereogenicity}(&StereoSymmetry)` against the stored constraints. DSL tags `#p` `#f` `#o` `#g`,
      parse/render of `[!][']perm` (`#p`/`#f`) and `[!](=|'|/)[pair]` (`#o`/`#g`); `^` migrates to cycle
      notation; validator rejects `'` on achiral kinds. *(Storage shape — flat vs `(B,N)` — still deferred.)*
    - **6 · wiring (umol-graph, mirrors the aromaticity model/resolver/validator split):** `ops/model.rs`
      `StereoModel` (`KindScope` + fluxionality + para + `InconsistencyPolicy`, producing
      `StereoSymmetryConfig`); `ops/resolver/stereo.rs` `StereoResolver` (adds+resolves stereo elements
      from `#T`/`#C` — kind-by-scope, ligand frame, coset resolved — with **idempotency = site-membership**
      (`has_coincident`/`coincident`, **not** `is_ground` — this adds entities, it doesn't narrow atoms) and
      a **configurable coverage/inconsistency policy `{ Keep | Strip | Error }`**, never a silent drop);
      `ops/validator/stereo.rs`
      `StereoValidator` **last** after `AromaticityValidator`, checking coset-in-range, ligand-arity =
      kind degree, the achiral `'` gate, and the asserted↔derived cross-check (computing `StereoSymmetry`
      **once**, molecule-wide, on the resolved AST). The validator is the **sole** pipeline consumer; the
      resolver is **structural** (computes none) — so there's nothing to share op-to-op, and a shared
      result would be stale anyway (the resolver mutates the molecule). Representation-independent
      `from_parts` / view carry of constraints. *(Non-stereo follow-up — **C4d**: `AromaticityResolver` needs the same membership-based
      idempotency + shared `InconsistencyPolicy` — today it re-perceives on re-run (duplicating systems) and
      silently discards proposals.)*
  - **Deliverable unchanged:** TH/CT first — both rigid (`R′=R`) under any fluxionality choice.
  - **Open (impl-level, decide at build):** the `StereoSymmetryConfig` struct vs bare params (lean struct).
    (Resolved: only the **validator** computes `StereoSymmetry` (once, on the resolved AST); the resolver is
    **structural** and computes none — nothing is shared op-to-op, and couldn't be since the resolver mutates
    the molecule; query-time reuse caches on the immutable `Molecule`, doc 086.) Π
    storage = `OrientedPermutationGroup`; constraint-AST names settled. Literal storage is the **flat
    non-unique collection** (the `RingSize` pattern); `(B,N)` is the lattice *model*, not a stored form
    — `B` is computed on demand from the positive literals via Π's group-closure (transient in
    `matches`, never stored), `N` is the stored negative literals. The lone `simplify` detail: whether
    to canonicalize the positive literals to a minimal generator set for `Eq`/dedup. Section C /
    StereoModel **design + naming converged**.

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
  - **D10** - add to prop test and fuzzing **Done** (stereo elements + `#T`/`#C` + relational + anchor in
     `molecule_ast_strategy`/`constraint_leaf_strategy`/`sub_pattern_anchor_strategy`; stereo entity-string +
     keyword roundtrip tests; `parse_stereo_atom`/`parse_stereo_bond` in `fuzz_entity_strings`)

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
`dsl/molecule.rs`; umol-io: `table_ir/raise.rs` (Phase B, with `umol-geometric-core` for the wedge geometry);
umol-graph: `ops/stereo.rs` (Phase C); umol-geometric: `src/{coordinates,molecule}.rs` (Phase F).

## Verification

Per-phase rstest table tests (relation/birelation generics + remap/unmap; types; raise `#T`/`#C` on
corpus; perception constraints→birelations; DSL round-trip; constraint↔element cross-check; matching).
End-to-end (A–E): `F[C@H](Cl)Br` / `F/C=C/F` → `#T`/`#C` → perceive → serialize `:stereo` → re-parse →
match. Run `mol/smiles/sdf_parsing` + resolution conformance (`--features conformance --test resolution`).

## Open

- Whether the `#T`/`#C` constraint is a stored `AtomConstraint`/`BondConstraint` kind (recommended,
  symmetric with `#a`) vs raise-transient.
