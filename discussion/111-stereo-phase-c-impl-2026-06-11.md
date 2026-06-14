# Stereochemistry Phase C — detailed implementation plan

Status: **Active / implementation plan.** Design source:
[104-stereochemistry-implementation-plan-2026-05-31.md](104-stereochemistry-implementation-plan-2026-05-31.md)
§C (perception → stereo elements: ligand symmetry, topicity, stereogenicity, fluxionality). This doc is the
per-crate, per-module **API** breakdown; implementation bodies are settled at coding time. No code is written
until explicitly authorized.

Crates, in dependency order: **C1 umol-perm**, **C2 umol-graph-core**, **C3 umol-ast**, **C4 umol-graph**.
Each subphase is one module.

---

## C1 · umol-perm **Done**

### C1a · `oriented.rs` (new) **Done**

Orientation-graded permutation and its group (proper subgroup + one improper coset rep) — the
permutation-inversion element `π·E^k` of Sₙ×Z₂.

```rust
/// Orientation grade: proper rotation vs improper (mirror) op.
enum Orientation { Proper, Improper }

impl Orientation {
    fn compose(self, other: Orientation) -> Orientation;  // XOR: Improper∘Improper = Proper
    fn flip(self) -> Orientation;
    fn is_proper(self) -> bool;
}

#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct OrientedPermutation { perm: Permutation, orientation: Orientation }

impl OrientedPermutation {
    fn new(perm: Permutation, orientation: Orientation) -> Self;
    fn proper(perm: Permutation) -> Self;
    fn improper(perm: Permutation) -> Self;
    fn identity(degree: usize) -> Self;          // proper identity
    fn perm(self) -> Permutation;
    fn orientation(self) -> Orientation;
    fn is_proper(self) -> bool;
    fn degree(self) -> usize;
    fn apply(self, i: usize) -> usize;           // = perm.apply(i) (orientation is not a point map)
    fn compose(self, other: Self) -> Self;       // perm.compose, orientation.compose
    fn inverse(self) -> Self;                     // perm.inverse, same orientation
}

/// A subgroup of Sₙ×Z₂, stored as the proper subgroup plus (if any) one improper coset
/// representative — the improper elements are `{ improper_rep ∘ p : p ∈ proper }`.
struct OrientedPermutationGroup {
    proper: PermutationGroup,
    improper_rep: Option<OrientedPermutation>,
}

impl OrientedPermutationGroup {
    fn generate(degree: usize, generators: &[OrientedPermutation]) -> Self;
        // close under composition; split into proper subgroup + one improper rep.
    fn degree(&self) -> usize;
    fn order(&self) -> usize;                     // proper.order() * (1 + improper_rep.is_some())
    fn contains(&self, op: OrientedPermutation) -> bool;
    fn elements(&self) -> Vec<OrientedPermutation>;          // proper ∪ improper coset
    fn proper(&self) -> &PermutationGroup;
    fn improper_rep(&self) -> Option<OrientedPermutation>;
    fn proper_orbit_of(&self, point: usize) -> Vec<usize>;   // orbit under proper        (homotopic)
    fn star_orbit_of(&self, point: usize) -> Vec<usize>;     // orbit under proper∪improper (star)
}
```

### C1b · `permutation.rs` (extend `Permutation`) **Done**

Cycle construct/decompose + GAP-notation `Display` (0-indexed). Pure additions to the existing
degree-≤6 `Copy` type.

```rust
impl Permutation {
    /// Disjoint cycles `[c0,…,ck]` set σ(c0)=c1,…,σ(ck)=c0; unlisted points fixed. `[]` = identity.
    fn from_cycles(degree: usize, cycles: &[Vec<usize>]) -> Self;
    /// Disjoint-cycle decomposition, fixed points dropped; canonical (each cycle least-element-first,
    /// cycles sorted by least element). Identity → `[]`.
    fn cycles(self) -> Vec<Vec<usize>>;
}

impl Display for Permutation;   // "(0,1,2)(3,4)"; identity → "()"
```

### C1c · `class.rs` (extend `ClassKey`) **Done**

Add `Axial`; `build()` also yields the **improper (orientation-reversing) generator** per class (the
source of chirality). `Axial` shares `CisTrans`'s parent/group/decomposition but with a coset-swapping
improper generator (allenes are chiral).

```rust
enum ClassKey {
    Symmetric(u8), Alternating(u8), Cyclic(u8), Dihedral(u8),
    Tetrahedral, CisTrans, Axial,            // Axial: shares CisTrans coset space, chiral
    SquarePlanar, TrigonalBipyramidal, Octahedral,
}

impl ClassKey {
    fn build(self) -> CosetSpace;            // build() yields the improper generator and passes it to CosetSpace::new
    // improper generator (StereoMolGraph `inversion`):
    //   Tetrahedral → (0 1)  [odd → chiral];  CisTrans → identity [achiral];
    //   Axial → (0 1)  [swaps cosets → chiral];  SquarePlanar → identity [achiral];
    //   TrigonalBipyramidal, Octahedral → geometric reflection generator [chiral] (class-geometry data).
}
// Display/FromStr: add "AX" ↔ Axial.
```

The TB/OH improper generators are class-geometry data (same source as today's hardcoded
`is_chiral_class`), not invented in this plan.

### C1d · `coset.rs` (extend `CosetSpace`) **Done**

Carry the improper generator; add chirality, enantiomer, and the merge primitive shared by
stereogenicity and fluxionality.

```rust
struct CosetSpace {
    parent: PermutationGroup,
    group: PermutationGroup,         // R (proper)
    numbering: HashMap<Permutation, u32>,
    representatives: Vec<Permutation>,
    improper: Permutation,           // orientation-reversing generator
}

impl CosetSpace {
    fn new(parent, group, decomposition, improper: Permutation) -> Self;   // new param

    fn enantiomer(&self, index: u32) -> u32;     // index(unindex(index) ∘ improper)
    fn is_chiral(&self) -> bool;                  // any coset moved by the improper op

    /// Quotient the coset space by extra proper generators (right action on positions): cosets i,j
    /// merge iff related by ⟨generators⟩. Returns, per coset, its class's canonical (min) index.
    /// Shared by stereogenicity (Π_proper) and fluxionality (R′ moves).
    fn merge_under(&self, generators: &[Permutation]) -> Vec<u32>;

    /// The observable coset under a fluxional supergroup: the merged class id of `index`.
    fn observable_coset(&self, index: u32, fluxional: &[Permutation]) -> u32;
}
```

Stereogenicity reads `merge_under(Π_proper)` — the stored coset is stereogenic iff its class is a
singleton. The `~` / `MirrorOp` involution stays **umol-ast-side** (C3), delegating to `improper` for
chiral kinds.

---

## C2 · umol-graph-core **Done**

### C2a · `algorithms/auto.rs` (extend `Automorphism`) **Done**

Capture and expose nauty/Traces generators (currently discarded); no other graph-core change. Bond /
relation labels stay **caller-side** — C3c folds them into node colors via a colored incidence gadget —
so `Graph::automorphisms` keeps its node-color-only signature and graph-core stays label-agnostic.

```rust
struct Automorphism {
    orbits: Vec<NodeId>, canonical_lab: Vec<NodeId>, node_count: usize, orbit_count: usize,
    group_order: AutoGroupOrder,
    generators: Vec<Vec<NodeId>>,            // new: each a permutation image over 0..node_count
}
impl Automorphism {
    fn generators(&self) -> &[Vec<NodeId>];  // new accessor
}
// automorphisms_nauty: register a nauty `userautomproc` callback that pushes each emitted generator
// (a permutation image) into `generators`; nauty-Traces-sys already exposes the hook. `Graph::automorphisms`
// signature is unchanged. (Traces backend would wire the same callback when added.)
```

---

## C3 · umol-ast

### C3a · `ast/coloring.rs` (new) + `Entity` in `ast/entity.rs` (new) **Done**

The pluggable **round-0 color policy**: one `u64` per **graph node** — atom or pseudonode (a relation:
localized bond or an overlay) — that says what counts as distinguishable. One rule, applied to every node:
**hash the entity's inherent fields**. This is *complete* for our consumer (a graph automorphism): every
**derived** predicate is relational, so it is a function of the structure and is preserved by the
automorphism for free — once each overlay is an **entity pseudonode** (C3c), even overlay-derived predicates
(aromaticity, dative/multicenter participation) are graph functions. So the color carries only the
**irreducible inherent labels**; nothing derived enters it. Named `MoleculeColoring` to avoid conflation
with graph-core's proper-coloring `algorithms/coloring.rs`.

Because the rule is uniform across node kinds, the trait is a **single method over a public `Entity`** (a
typed ref to any molecule entity — lives in `ast/entity.rs`, generally useful beyond coloring; variant names
mirror the entity types exactly):

```rust
// ast/entity.rs  (public, general-purpose)
enum Entity {
    Atom(AtomId), Bond(BondId), DativeBond(DativeBondId), AromaticSystem(AromaticSystemId),
    MulticenterBond(MulticenterBondId), NoncovalentBond(NoncovalentBondId),
    StereoAtom(StereoAtomId), StereoBond(StereoBondId),
}

// ast/coloring.rs
trait MoleculeColoring {
    fn color(&self, mol: &MoleculeAst, entity: Entity) -> u64;   // stereo node → kind tag + geometric kind; coset folded by C3d
}

bitflags! struct ConstitutionFeatures: u32 {              // inherent fields only (§6.1) — derived excluded (free)
    ELEMENT, ISOTOPE, CHARGE, IMPLICIT_HYDROGENS, LONE_PAIRS, SPIN,   // atom
    BOND_ORDER, BOND_CHARGE, BOND_SPIN,                       // localized bond
    DATIVE_ORDER,                                             // dative bond
    AROMATIC_ELECTRONS, AROMATIC_CHARGE, AROMATIC_SPIN,       // aromatic system (electrons = order-indep π-count)
    MULTICENTER_ELECTRONS, MULTICENTER_CHARGE, MULTICENTER_SPIN, // multicenter bond
    NONCOVALENT_KIND,                                         // noncovalent bond
    STEREO_KIND,                                              // stereo atom/bond geometric kind (partition-free)
    // bit values are arbitrary (runtime config, never serialized); append freely.
}

struct ConstitutionColoring { features: ConstitutionFeatures }   // the one impl now (pre-stereo inherent fields)
impl ConstitutionColoring { fn new(features: ConstitutionFeatures) -> Self; fn full() -> Self; }
impl MoleculeColoring for ConstitutionColoring;
// color(): match entity kind → hash (kind tag, that kind's inherent fields, each gated by ConstitutionFeatures):
//   Atom → element/isotope/charge/#h/lone-pairs/spin;  Bond → order/charge/spin;
//   Aromatic → π-count/charge/spin;  Multicenter → e-count/charge/spin;
//   Dative → order (direction is gadget-encoded, C3c);  Noncovalent → kind.
//   StereoAtom/StereoBond → kind tag + geometric kind (STEREO_KIND); partition-dependent coset folded in C3d.
//   per-atom electron vectors are NOT node fields (order-dependent) → incidence-edge colors, C3c.
// the kind tag keeps kinds in disjoint color ranges (no pseudonode maps onto an atom or another kind).
```

`ConstitutionColoring` serves **our** distinctness / symmetry analysis. `MorganColoring`/`EcfpColoring`/`FcfpColoring`
are **separate** impls owning their own fields (ring, aromaticity, functional class) — they carry derived
props because their consumer is per-atom hashing, not an automorphism, so "relational ⇒ free" doesn't apply
to them. Stereo is **not** a coloring impl: its descriptor is partition-dependent, so the symmetry-refinement loop
folds it on top (C3d); the trait stays stereo-free (and fingerprint-reusable). Stereo atoms/bonds **are**
graph nodes (present in the full incidence selection); `ConstitutionColoring` gives them their kind tag +
partition-free geometric `kind` (gated by `STEREO_KIND`, like any field). Their partition-dependent observable
coset is the single term C3d folds on top — which is exactly why stereo is not itself a coloring impl.

### C3b · `ast/automorphism.rs` (extend) **Superseded — reverted by C3d (2026-06-13)**

Originally: expose generators on `AtomAutomorphism`; add the orientation-graded holder. The corrected C3d
grades **graph-core** generators over the incidence graph (`NodeId`) directly, so these atom-level wrappers are
unused — the additions below are reverted; the pre-existing `AtomAutomorphism` (atom-only molecule graph) stays.
Kept for the design record:

```rust
impl AtomAutomorphism {
    fn generators(&self) -> Vec<Vec<AtomId>>;        // graph-core generators as AtomId maps
}

struct GradedAutomorphism {
    proper: AtomAutomorphism,            // Â  (orientation-preserving)
    improper: Option<Vec<AtomId>>,       // Â* — one molecule→mirror rep; None ⇒ chiral (a coset has no generators)
}
impl GradedAutomorphism {
    fn new(proper: AtomAutomorphism, improper: Option<Vec<AtomId>>) -> Self;
    fn proper(&self) -> &AtomAutomorphism;
    fn improper_rep(&self) -> Option<&[AtomId]>;
}
```

### C3c · `ast/incidence.rs` (new) — generic **Done**

The molecule's **incidence graph** (Levi graph): each selected relation — localized bond and each overlay
(aromatic / multicenter / dative / noncovalent) — lifted to a pseudonode alongside the atoms, with incidence
edges to its participants. A **generic** `MoleculeAst` property (canonicalization, fingerprints, symmetry
numbers all reuse it), not stereo-specific; structure only — colors are applied at automorphism time.
Hyperedge overlays (>2 atoms) are why it's the incidence/Levi graph rather than an ordinary-edge subdivision;
it is **not** a port graph — attachment stays unordered (ligand order enters only in the per-site projection,
C3d). Localized bonds are lifted (not kept as edges) so their `BOND_*` labels can ride as node colors under
graph-core's node-color-only automorphism (C2a).

The node set is **configurable** over a forced atoms+bonds base: `IncidenceNodeSelection` toggles `OVERLAYS`
and `STEREO` (presets `topological` / `constitution` / `full`). Stereo atoms/bonds **are** graph nodes (just
not part of the constitution coloring); each attaches to **its site only** (the site atom, or the site bond's
pseudonode) — the ligand–site bonds already live in the topology, so re-linking ligands through the stereo
node would duplicate them; the new information a stereo node carries is its site plus (at color time) its
stereo label.

**No dative-direction gadget.** A directed bond between automorphism-equivalent atoms is a contradiction in
terms (dative-ness *is* the donor/acceptor asymmetry), so the coloring already separates the endpoints; the
direction itself is retained in the AST (`acceptor_slot`). Noncovalent relations are stored unordered and
need nothing. So no marker nodes.

Nodes are laid out in fixed kind-blocks (atoms, bonds, then each selected overlay/stereo kind in order), so
the node↔entity correspondence is fully determined by the **per-kind counts** — store those, not a per-node
vec. `entity(node)` / `node_of(entity)` are O(1) offset arithmetic (atoms are identity, the rest subtract a
block offset).

```rust
bitflags! struct IncidenceNodeSelection: u8 { OVERLAYS, STEREO }   // atoms + localized bonds always in
impl IncidenceNodeSelection { fn topological()/*empty*/; fn constitution()/*OVERLAYS*/; fn full()/*both*/; }

struct IncidenceGraph {
    graph: Graph,         // CSR (graph-core); nodes = atoms ++ selected relation-pseudonodes
    counts: [u32; 8],     // per-kind block sizes (Atom, Bond, Dative, Aromatic, Multicenter, Noncovalent, StereoAtom, StereoBond)
}
impl IncidenceGraph {
    fn graph(&self) -> &Graph;
    fn entity(&self, node: NodeId) -> Entity;     // O(1) from block offsets
    fn node_of(&self, entity: Entity) -> NodeId;  // inverse
}
impl MoleculeAst {
    fn incidence_graph(&self, selection: IncidenceNodeSelection) -> IncidenceGraph;
}
```

### C3d · `ast/symmetry.rs` (new) — two steps **Done**

The molecule's graph-automorphism symmetry under a coloring, **graded** into proper vs improper, and the
per-carrier projection derived from it. **Two results at two granularities:**

- **`GraphSymmetry`** — the general, *owned* artifact (RingSet rule: artifacts own, views borrow). The
  converged automorphism over the **full** incidence graph (C3c) under a coloring, graded. Self-contained: it
  stores the incidence graph + converged colors + proper orbits + star orbits + one improper rep, so every
  query (molecule-level orbits **and** per-carrier stabilizer re-runs) works from it alone — no `&MoleculeAst`
  ref, no generic-`C` leak (the coloring is applied during construction, then erased to the color vector).
  Reusable beyond stereo; cacheable later. Node ids never leave it — public queries speak `AtomId`; node↔entity
  is the compact block layout (C3c).
- **`StereoSymmetry`** — the compact per-carrier projection (`OrientedPermutationGroup` on ≤6 ligand positions +
  kind + stored coset). The assertion currency; produced from a `GraphSymmetry`, carries the predicates.

**Why grading, not a single partition (corrected 2026-06-13).** A scalar stereo-node color cannot forbid an
orientation-*reversing* ligand swap — the node is a degree-1 pendant carrying one number, and a swap that
inverts handedness leaves it unchanged. So nauty over the folded coloring returns the **full** automorphism
group, and its orbit array is the **star** partition (proper ∪ improper), *not* proper. Example: C(Cl,Cl,F,Br)
with a tetrahedral coset — the two Cl are enantiotopic, but the Cl↔Cl swap preserves the colored graph, so
nauty merges them; calling that "proper" reports them homotopic. To recover proper orbits we **grade the
generators**: a generator's orientation = its net action on stereocenter cosets (transport the stored coset by
the induced ligand permutation via `StereoCosetAst::apply_permutation` — equals the target center's coset ⇒
proper there, equals its enantiomer ⇒ improper; uniform across all centers ⇒ the generator's grade). This
**subsumes Route A**: an improper element exists iff some generator is improper (proper is index-2 normal), so a
separate molecule→mirror canonical run and `mirror_colors` are unnecessary.

**Engine + known limit.** nauty (fast) over the folded coloring + generator grading. Exact except for **false
automorphisms** at molecules with ≥2 independent prochiral centers: collapsed observable cosets can admit a
*mixed* generator (proper at one center, improper at another) that is not a real symmetry. Mixed generators are
**discarded** in grading (orbits become correctly finer); the residual gap (a true symmetry expressible only via
a mixed generator) is pathological and is the deferred **VF2-carrying-parity** exact path (doc 104). For the
common case (≤1 independent prochiral axis, resolved centers) this is exact.

**Step 1 — `GraphSymmetry` (build + molecule-level queries).**

```rust
struct GraphSymmetryConfig<C: MoleculeColoring> { coloring: C, iterate_to_fixpoint: bool, max_iterations: usize }

struct GraphSymmetry {
    incidence: IncidenceGraph,         // full selection
    colors: Vec<u64>,                  // converged node colors (for step-2 stabilizer re-runs)
    proper_orbits: Vec<NodeId>,        // orbit rep per node, under proper generators
    star_orbits: Vec<NodeId>,          // under proper ∪ improper generators
    chiral: bool,
}

impl MoleculeAst {
    fn graph_symmetry<C: MoleculeColoring>(&self, cfg: &GraphSymmetryConfig<C>) -> GraphSymmetry;
    // inc = self.incidence_graph(IncidenceNodeSelection::full());
    // base = cfg.coloring.color(self, inc.entity(node))                       // static, all kinds
    // fixpoint: recolor(node) = base ⊕ (stereo ? observable_coset(elem, orbits_{k-1}) : 0);
    //   orbits_k = inc.graph.automorphisms(recolor).orbit array (refinement signal, star granularity);
    //   stop when stable, else one pass if !iterate_to_fixpoint; max_iterations caps.
    // grade the converged run's generators by coset action → proper set / improper set (discard mixed);
    // proper_orbits = union-find(proper set); star_orbits = union-find(proper ∪ improper);
    // chiral = (some stereocenter is definite) && (no improper generator)  — an improper generator ⇒ achiral,
    //   but a molecule with no stereocenters is trivially achiral despite having none, hence the first clause.
}

impl GraphSymmetry {
    fn same_proper_orbit(&self, a: AtomId, b: AtomId) -> bool;     // proper_orbits[a] == proper_orbits[b]
    fn same_star_orbit(&self, a: AtomId, b: AtomId) -> bool;       // star_orbits[a] == star_orbits[b]
    fn proper_orbit_of(&self, a: AtomId) -> Vec<AtomId>;
    fn star_orbit_of(&self, a: AtomId) -> Vec<AtomId>;
    fn is_chiral(&self) -> bool;
    pub(crate) fn site_stabilizer(&self, site: NodeId) -> Vec<Vec<NodeId>>;  // raw stabilizer gens, carrier bumped (step 2)
}
```

The orientation grading lives **molecule-side** (it needs coset transport over the stereo elements) and is
shared by the build and the step-2 projection. `AtomAutomorphism::generators()` + `GradedAutomorphism` (C3b)
were scoped for an atom-level graded design; the corrected build grades **graph-core** generators over the
incidence graph (`NodeId`) directly, so those umol-ast wrappers are unused — **revert the C3b additions** (keep
the pre-existing `AtomAutomorphism`).

**Step 2 — `StereoSymmetry` (per-carrier projection + predicates).**

Entry on `MoleculeAst` (it owns the stereo elements / ligand order); `GraphSymmetry` is consulted and stays a
pure symmetry object.

```rust
enum Topicity { Homotopic, Enantiotopic, Diastereotopic }

struct StereoSymmetry { group: OrientedPermutationGroup, kind: StereoKind, coset: StereoCosetAst }

impl MoleculeAst {
    fn stereo_atom_symmetry(&self, gs: &GraphSymmetry, id: StereoAtomId) -> StereoSymmetry;
    fn stereo_bond_symmetry(&self, gs: &GraphSymmetry, id: StereoBondId) -> StereoSymmetry;
    // site = atom node (atom carrier) | bond pseudonode (bond carrier).
    // raw = gs.site_stabilizer(site); grade raw molecule-side → proper stabilizer gens + optional improper rep.
    // project each onto stored ligand positions (atom ligand p ↦ position holding σ(atom_p));
    //   add Sₖ per same-kind virtual-ligand block.
    // OrientedPermutationGroup::generate(ligand_count, proper ∪ improper rep).
}

impl StereoSymmetry {
    fn is_stereogenic(&self) -> bool;     // stored coset is a singleton class of merge_under(Π underlying perms,
                                          // proper ∪ improper) — see note
    fn topicity(&self, a: usize, b: usize) -> Topicity;  // same Π⁺ orbit → homotopic; same star, not proper → enantiotopic; else diastereotopic
}
```

**`is_stereogenic` uses the full local group, not just Π⁺ (corrected 2026-06-13).** Doc 104's
"singleton class of `merge_under(Π_proper)`" is necessary but not sufficient: it catches centers killed by a
*proper* ligand symmetry (homotopic ligands) but labels **prochiral** centers (enantiotopic ligands, e.g.
C(Cl,Cl,F,Br)) as stereogenic, since their identifying symmetry is *improper* (Π⁺ is trivial there). A center
is genuinely stereogenic iff its stored coset is not identified with any other coset by the local symmetry —
proper **or** improper. So merge over the underlying permutations of the whole oriented group `Π_s`; the coset's
class is a singleton iff stereogenic. (An improper transposition of two identical ligands maps the coset to its
enantiomer, collapsing the class — correctly non-stereogenic.)

### C3e · `ast/stereo.rs` (extend) **Done**

Add `Axial`; delegate chirality to umol-perm; add `MirrorOp` (`'`); add the two ground relation enums
(no trait — glyph mapping lives in the DSL, C3h).

Done notes: `is_chiral_class` and the `involution`'s chiral generator now delegate to umol-perm (added
`CosetSpace::improper()` accessor; achiral kinds keep their chosen swap). `Topicity`/`Stereogenicity` live
here (`Topicity` moved out of `symmetry.rs`). The DSL kind glyph for `Axial` is provisionally `Ax` and
`MirrorOp` renders/parses as `'` (in `dsl/stereo.rs`, render+parse symmetric for roundtrip) — C3h owns the
final glyph table and may revise.

```rust
enum StereoKind { Tetrahedral, CisTrans, Axial, SquarePlanar, TrigonalBipyramidal, Octahedral }  // +Axial

impl StereoKind {
    fn class_key(self) -> ClassKey;          // +Axial → ClassKey::Axial
    fn is_chiral_class(self) -> bool;        // now delegates to space(class_key).is_chiral()  (C1d)
    fn involution(self) -> Permutation;      // chiral kinds → improper generator (umol-perm); achiral → chosen swap
}

enum StereoExpr {                            // +MirrorOp (the `'` improper op; involution, folds via enantiomer)
    Lit(u32), Var(String), SwapOp(Box<StereoExpr>), ApplyOp(Box<StereoExpr>, Permutation),
    MirrorOp(Box<StereoExpr>), LitSet(Vec<u32>), VarDomain(String, Vec<u32>),
}

enum Topicity { Homotopic, Enantiotopic, Diastereotopic }      // derived ground (Lit payload)
enum Stereogenicity { Symmetric, Prochiral, Stereogenic }      // derived ground (Lit payload)
```

### C3f · `ast/views/stereo.rs` (extend); `StereoLigandId` in `ast/ids.rs` **Done**

Done notes: `p` is the **per-carrier** `StereoSymmetry` (C3d), pre-computed by the op and passed to the view's
pure-read queries — the methods delegate (`symmetry.topicity` etc.); `is_chiral`/`ligand_position`/
`ligand_frame` use only the view. The two views share one `stereo_view_queries!` macro. `StereoSymmetry` gained
`stereogenicity()` + `group()`/`kind()`/`coset()` accessors, and `topicity` migrated to `StereoLigandId`.
`StereoLigandId` is a hand-written `u8` newtype (not `define_id!`, which is `u32`); no EDN derives yet — add in
C3g if the constraint needs them. `ligand(idx)` migrated to `ligand(ligand_id: StereoLigandId)`.
**Deferred**: `induced_ligand_permutation` — it overlaps C3d's internal projection (`project_onto_ligands` +
grading) and has no consumer yet; expose it (with a shared helper) when one appears.

Ligand positions get a newtype — bare `usize` is too loose for the id-rich design. `StereoLigandId` is the
0-based position in a stereo element's **ordered ligand frame** (frame-relative, not a global id). It is `u8`,
matching umol-perm `Permutation`'s position width (the `[u8; 6]` image, bounded by the kind's degree) — *not*
`usize`; the boundary conversion is just `id.0` / `StereoLigandId(_)`. Defined in `ast/ids.rs` beside `AtomId`
etc.; existing position-taking view methods (`ligand(idx)`, `permutation_for`, `coset_for`) migrate to it too.

```rust
// ast/ids.rs
struct StereoLigandId(u8);   // position in a stereo element's ligand frame (≤ kind degree, Permutation-width)

// ast/views/stereo.rs — read-only queries taking &StereoSymmetry (the consuming op computes it and passes
//   it to its own view queries; ops do NOT share one symmetry result — resolver and validator are independent)
impl StereoAtomView<'_> {                    // StereoBondView parallel
    fn ligand_symmetry(&self, p: &StereoSymmetry) -> OrientedPermutationGroup;   // p.ligand_symmetry(site, frame)
    fn stereogenicity(&self, p: &StereoSymmetry) -> Stereogenicity;
    fn is_stereogenic(&self, p: &StereoSymmetry) -> bool;   // stored coset is a singleton merge_under(Π_proper) class
    fn is_chiral(&self) -> bool;                              // kind-level (umol-perm), no symmetry computation
    fn topicity(&self, a: StereoLigandId, b: StereoLigandId, p: &StereoSymmetry) -> Topicity;
    fn is_homotopic(&self, a: StereoLigandId, b: StereoLigandId, p: &StereoSymmetry) -> bool;
    fn is_enantiotopic(&self, a: StereoLigandId, b: StereoLigandId, p: &StereoSymmetry) -> bool;
    fn is_diastereotopic(&self, a: StereoLigandId, b: StereoLigandId, p: &StereoSymmetry) -> bool;
    fn is_prochiral(&self, p: &StereoSymmetry) -> bool;
    // projection helpers
    fn ligand_position(&self, atom: AtomId) -> Option<StereoLigandId>;
    fn ligand_frame(&self) -> Vec<StereoLigand>;                     // the ordered frame (atoms + virtual h/lp)
    fn induced_ligand_permutation(&self, atom_perm: &[AtomId]) -> OrientedPermutation;  // atom map → position perm
}
```

### C3g · `ast/constraint/stereo.rs` (extend — inhabit the empty enums)

**Staged plan (2026-06-13). Decisions locked:** (1) `relation_ast!` generates a **4-variant** subset lattice
`{ Undetermined /*full domain*/, Lit(T), LitSet(Vec<T>), NotSet(Vec<T>) }` (no special `Not`); meet=∩, join=∪,
matches=∈; domain baked so full-set ⇒ `Undetermined`. (2) Collections have the **same shape as
`AtomConstraints`** — flat storage, pos/neg encoded *in the constraint* (`LigandSymmetryAst.mem`), **no `(B,N)`
storage** (that's a matching-side concern, not storage). (3) **Macro-generate** the atom/bond parallel
constraint enum + collection now (split later if they diverge).

- **C3g.1 — leaf value types** (`PermutationAst`, `OrientedPermutationAst`, `MemOp`, `LigandPairAst`): concrete
  `Eq`/`Hash` + inherent `matches`; `LigandPairAst::new` normalizes (private fields + `first()`/`second()`).
  PermutationAst EDN host deferred to C3h (constraints round-trip via the inline stereo-string, not per-type EDN). **Done**
- **C3g.2 — `relation_ast!` macro** + `TopicityRelationAst`, `StereogenicityRelationAst` (Lattice + AsLit), the
  4-variant subset lattice above. **Done** — domain via `strum::VariantArray` on the domain enum (which also
  derives `Ord`); equality/hash/lattice ops are by the represented set (representation-independent), so
  `LitSet([H,E]) == NotSet([D])`; `from_set` canonicalizes (full⇒`Undetermined`, singleton⇒`Lit`, else smaller
  of subset/complement).
- **C3g.3 — constraint variant types**: `LigandSymmetryAst`, `FluxionalityAst` (concrete + `matches`),
  `TopicityAst` (keyed Lattice — `meet`/`join` debug-assert same pair, operate on the relation),
  `StereogenicityAst` (Lattice + AsLit, delegates to the inner relation). **Done**
- **C3g.4 — inhabit the unions**: `StereoAtomConstraint`/`StereoBondConstraint` (4 variants) +
  `…Kind` + `kind()`/`is_unique()`, replacing today's uninhabited `enum {}` (macro-generated atom/bond). **Done.**
  Inhabiting forced the `Constraint`-level exhaustive matches: `molecule.rs::inline_constraints` (now
  `stereo_{atom,bond}_mut(id).constraints.add(inner)`), `constraint/molecule.rs::simplify` (pass-through — no
  `ValueAst`) and `::remap` (`remap.stereo_{atom,bond}(id)?` + no-op constraint remap, positions are
  frame-relative). The collection `add` is a minimal push (cardinality is C3g.5). The `…Kind` re-exports are
  withheld until C3g.5 consumes them (would warn unused otherwise). The four DSL sites — `dsl/constraint.rs`
  (molecule-level) and `dsl/stereo.rs` (inline) `match c {}` over the now-inhabited enums — are **left red**:
  rendering/parsing stereo constraints is C3h's `ConstraintDsl::StereoAtom` + inline stereo-string work, done
  once there rather than as a throwaway `Err`. Tree is red on exactly those four sites until C3h.
- **C3g.5 — collections** `StereoAtomConstraints`/`StereoBondConstraints`: full surface mirroring
  `AtomConstraints` + `Lattice` (per-kind meet/join: `#p`/`#f` union+dedup, `#o` per-pair value-meet, `#g`
  value-meet). **Done.** Folded into `stereo_constraint!` (3rd param = collection name) since atom/bond are
  identical: kind-sorted `SmallVec<[_; 2]>`, `add` policy = `#g` unique-replace / `#o` keyed-replace per
  `LigandPairAst` / `#p`/`#f` append. `with_constraints` on the element (was a no-op stub) now wires
  `extend`. Verified by type-check (`cargo build --tests` clean apart from the four C3h DSL sites); tests
  cover add cardinality, kind-sort, accessors, meet (incl. `None` clash), join, matches, and
  is_undetermined/is_ground — they run once C3h closes the red.

The constraint AST. Trait coverage follows the codebase convention: **value-bearing AST types impl
`Lattice` (`is_undetermined`/`is_ground`/`meet`/`join`/`matches`) + `AsLit`; keyed values impl `Lattice`
(à la `JointDomainAst`); collections impl `Lattice`; the constraint enum impls neither** (the collection
dispatches). The **concrete literal wrappers can't be lattices** (no top → no `join`), but pattern matching
still needs `matches`, so each gets a **standalone inherent `matches(&self, target: &Self) -> bool`** (the
existing `StereoExpr::matches_value` precedent — matchable, no `Lattice`). `matches` stays a directly-written
method everywhere (the `Lattice` doc already mandates direct impls, not `meet`-delegation), so no `Matches`
trait is extracted and `matches`/`meet` stay co-located for the lattice types.

```rust
struct PermutationAst(Permutation);                       // concrete literal — Eq/Hash + inherent matches; hosts EDN; NOT Lattice
struct OrientedPermutationAst { perm: PermutationAst, orientation: Orientation }   // concrete — Eq/Hash + inherent matches
enum MemOp { In, NotIn }                                   // local (not value.rs MemOp); plain enum
struct LigandPairAst { first: StereoLigandId, second: StereoLigandId }   // key; Eq/Ord/Hash; new(a,b) normalizes; NOT Lattice

// the macro generates the 4-variant subset lattice over a 3-atom enum (domain baked for the
// full-set → Undetermined collapse) AND its `Lattice` (is_undetermined/is_ground/meet=∩/join=∪/matches=∈)
// + `AsLit { type Lit = <enum>; as_lit() = Some(t) when Lit(t) }`. No trait, no duplication; glyph map in DSL (C3h).
relation_ast! { TopicityRelationAst,       Topicity }       // impl Lattice + AsLit  (Lit = Topicity)
relation_ast! { StereogenicityRelationAst, Stereogenicity } // impl Lattice + AsLit  (Lit = Stereogenicity)

struct LigandSymmetryAst { perm: OrientedPermutationAst, mem: MemOp }   // ± literal over Π — concrete; Eq/Hash + inherent matches; NOT Lattice
struct FluxionalityAst { perm: PermutationAst }                        // proper move — concrete; Eq/Hash + inherent matches; NOT Lattice
struct TopicityAst { pair: LigandPairAst, rel: TopicityRelationAst }   // keyed → impl Lattice (JointDomainAst-style; rel meet/join per pair)
struct StereogenicityAst(StereogenicityRelationAst)                    // impl Lattice + AsLit (delegate to the inner)

enum StereoAtomConstraint {            // tagged union — impls neither Lattice nor AsLit (like AtomConstraint)
    LigandSymmetry(LigandSymmetryAst), Fluxionality(FluxionalityAst),
    Topicity(TopicityAst), Stereogenicity(StereogenicityAst),
}
impl StereoAtomConstraint { fn kind(&self) -> StereoAtomConstraintKind; fn is_unique(&self) -> bool; }
enum StereoBondConstraint { /* parallel */ }

// The collection mirrors AtomConstraints' surface — the constraint side of the derived↔stored pair.
impl StereoAtomConstraints {
    // stored-side accessors (the pair's constraint half; default Undetermined / empty if absent):
    fn ligand_symmetry(&self) -> impl Iterator<Item = &LigandSymmetryAst>;   // #p non-unique
    fn fluxionality(&self) -> impl Iterator<Item = &FluxionalityAst>;        // #f non-unique
    fn topicities(&self) -> impl Iterator<Item = &TopicityAst>;              // #o keyed per pair
    fn topicity(&self, pair: LigandPairAst) -> TopicityRelationAst;          // Undetermined if absent
    fn stereogenicity(&self) -> StereogenicityRelationAst;                   // #g unique; Undetermined if absent
    // management (parallel to AtomConstraints):
    fn new() -> Self; fn is_empty(&self) -> bool; fn len(&self) -> usize;
    fn contains(&self, kind: StereoAtomConstraintKind) -> bool;
    fn get(&self, kind) -> Option<&StereoAtomConstraint>; fn get_all(&self, kind) -> impl Iterator<…>;
    fn add(&mut self, c: StereoAtomConstraint) -> Option<StereoAtomConstraint>;   // per-kind/per-pair insert policy
    fn extend(…); fn retain(…); fn clear(&mut self); fn take(&mut self) -> impl Iterator<…>;
    fn remove(&mut self, kind) -> …; fn remove_all(&mut self, kind) -> Vec<…>; fn simplify_each(&mut self);
    fn iter(&self) -> …; fn iter_mut(&mut self) -> …; fn remap(self, _: &IdRemapping) -> Self;  // positions are frame-relative → no-op
}
impl Lattice for StereoAtomConstraints {
    fn is_ground(&self) -> bool;
    fn meet(&self, other) -> Option<Self>;   // #p/#f union+dedup; #o per-pair value-meet (None on clash, e.g. =(i,j)∧/(i,j)); #g value-meet
    fn join(&self, other) -> Self;           // #p/#f intersection; #o per-pair value-join; #g value-join
    fn matches(&self, target: &Self) -> bool; // AST↔AST, exactly like AtomConstraints (stored-vs-stored)
}
```

The **derived↔stored cross-check is the C4 validator's job**, not a constraint method (mirroring how the
atom validator cross-checks the topology-derived vs stored valence pairs): for each element it compares the
derived `view.{ligand_symmetry,topicity,stereogenicity}(&StereoSymmetry)` against the stored constraints,
erroring when both are ground and inconsistent. Query matching uses the standard `Lattice::matches`, with the
matcher feeding the target's derived stereo through the view's derived accessors. (No `matches_derived`.)

*(The validator (C4) also rejects `'` on achiral kinds. Storage shape — flat literal `Vec` vs `(B,N)` — deferred.)*

### C3h · `dsl/stereo.rs` + `dsl/constraint.rs` (extend — stereo-element constraint serialization)

Two surfaces: the **inline stereo-string** (`:type` payload, `dsl/stereo.rs`) and the **molecule-scope
structured EDN** (`:stereo-atom`/`:stereo-bond` keys, `dsl/constraint.rs`). The four currently-`match {}`
compile sites (`StereoAtomConstraintDsl`/`StereoBondConstraintDsl` `from_ast` at `dsl/stereo.rs:520`/`556`,
`ConstraintDsl::from_ast` at `dsl/constraint.rs:1904`/`1905`) belong to the structured surface.

Substeps:

- **C3h.1 — structured-EDN codecs** (`dsl/stereo.rs`): permutation ↔ vector-of-cycles EDN
  (`[[0 1 2] [3 4]]`, identity `[]`, via `Permutation::{from_cycles, cycles}`; degree not encoded — supplied
  by the reader from the kind; `perm_from_vov` validates range + disjointness itself since `from_cycles`
  panics on a non-bijection); relation ↔ keyword/set EDN (`:homotopic` / `#{:enantiotopic :diastereotopic}`
  / `:undetermined`) generated by a `relation_codec!` macro keyed on a per-type variant→keyword table. The
  AST carries no keywords; `relation_ast!`'s `to_set`/`from_set` are now `pub(crate)` for the codec. **Done**
  (type-checks; unused until C3h.2). The inline GAP-string + glyph codecs moved to **C3h.3** (they land with
  the inline grammar that consumes them), not here.
- **C3h.2 — molecule-scope structured form** (`ast/constraint/molecule.rs` + `dsl/stereo.rs` +
  `dsl/constraint.rs`). **Done.**

  **Where the kind lives (degree provenance).** The molecule-level `Constraint` carries the stereo **kind**
  beside the id — `Constraint::StereoAtom(StereoAtomId, StereoKind, StereoAtomConstraint)` (and `StereoBond`
  parallel). Rationale: stereo elements are subtyped by kind, and `StereoKind::degree` is many-to-one
  (Tetrahedral and SquarePlanar are both degree 4; the bond kinds collide too), so a permutation payload
  cannot recover its kind — a constraint *detached from its element* (the molecule-scope form) is not
  well-formed without the subtype tag. The kind lives **only** at the molecule level; the entity-level inline
  store stays bare `StereoAtomConstraint` (the element supplies the kind). `lift_constraints` copies the
  element's kind in (infallible); `inline_constraints` drops it and stays infallible — the **kind/degree
  consistency check is deferred to the C4 validator** (which already cross-checks stored-vs-derived and has the
  whole molecule). `is_vacuous` now delegates to the inner `is_undetermined` for the stereo leaves.

  **DSL boundary object.** `StereoAtomConstraintDsl(pub StereoKind, pub StereoAtomConstraint)` — a **tuple
  wrapper** that owns its full, self-contained EDN, so it impls the normal `FromEdn`/`ToEdn` and the generic
  2-field entity-leaf machinery (`parse_entity_leaf`/`entity_leaf_edn`/`read_entity_leaf`) applies unchanged.
  Mirror at the envelope: `ConstraintDsl::StereoAtom(StereoAtomRef, StereoAtomConstraintDsl)` (the kind rides
  inside the wrapper, not as a separate `ConstraintDsl` field); `from_ast`/`into_ast` construct/destructure the
  tuple inline (kind from the `Constraint` envelope). **No degree-free DSL intermediate, no `EntityCounts`
  change, no `MoleculeInput` change, no bespoke leaf parser** — the streaming reader bridges via
  `read_value_slice` + `read_string` + `FromEdn`.

  *Design history note:* a long exploration considered (a) putting the kind as a separate `ConstraintDsl`
  field / a 3-element leaf `[ref kind payload]`, which would break the generic leaf flow and force a bespoke
  parser + lose the payload `FromEdn`; and (b) the kind-inside-payload self-contained wrapper. The deciding
  insight was that EDN layout and struct layout are independent — `FromEdn` is free-form — so the wrapper keeps
  the kind *and* stays self-contained. Chosen EDN: flat payload map (not a nested `[kind payload]` vector).

  **EDN form** (the entry in the top-level `:constraints []`): `{:stereo-atom [<ref> {:kind <kind-kw>
  <constraint-kw> <value>}]}` — a normal `[ref payload]` leaf; the payload is one map carrying `:kind` plus the
  single constraint key. Kind keywords are the kebab `StereoKind` names (`:tetrahedral`, `:cis-trans`, `:axial`,
  `:square-planar`, `:trigonal-bipyramidal`, `:octahedral`). Examples:
  ```clojure
  {:stereo-atom [0 {:kind :tetrahedral :fluxionality [[0 1]]}]}
  {:stereo-atom [1 {:kind :tetrahedral :ligand-symmetry {:perm [[0 1]] :orientation :improper :member :not-in}}]}
  {:stereo-atom [0 {:kind :octahedral :topicity {:pair [0 1] :relation :enantiotopic}}]}
  {:stereo-atom [0 {:kind :tetrahedral :stereogenicity :stereogenic}]}
  ```
  (`:orientation` defaults `:proper`, `:member` defaults `:in` — omitted when default.) Verified by
  `test_constraint_dsl_stereo_atom_roundtrip` (5 cases: render-equals-expected-EDN **and** parse-back-equals-AST,
  across all four constraint kinds + the defaults). Full `umol-ast` suite green (3216 tests); workspace builds.
- **C3h.3 — inline stereo-string predicates** (`dsl/stereo.rs`). **Done.** Extends
  `parse_stereo_atom`/`fmt_stereo_atom` (and bond) so the `:type` payload is `class coset predicate*`. Note
  `~`/`'`/`^` are already coset-expression operators (swap/mirror/apply-1-indexed-image); the predicates get
  their own payload grammar that coexists, and the coset stays unchanged (predicates use 0-indexed
  disjoint-cycle notation — that is the doc's "`^` migrates to cycle notation").

  Grammar (degree = `kind.degree()` from the class):
  ```
  #p  ligand-symmetry   [!][']<cycles> | [!]~     ! → member NotIn (else In);  ' → Improper (else Proper)
  #f  fluxionality        <cycles> | ~            (no !/' — FluxionalityAst is a bare proper perm)
  #o  topicity           (* | [!]<glyph>)(i,j)    * → Undetermined; ! → NotSet([v]); else Lit(v); (i,j) = pair
  #g  stereogenicity     (* | [!]<glyph>)         * → Undetermined; ! → NotSet([v]); else Lit(v)
  <cycles>  product of disjoint cycles, 0-indexed, identity ()  e.g. (0,1,2)(3,4)  [Permutation::Display/from_cycles]
  <glyph>   =  '  /     #o: = homotopic ' enantiotopic / diastereotopic
                        #g: = symmetric ' prochiral    / stereogenic
  ```
  `(* | [!]<glyph>)` is exactly complete over the 3-element relation domain: `*` = `Undetermined`, glyph =
  `Lit`, `!`glyph = `NotSet([v])` (the 2-element complement) — all 7 non-empty subsets. Render: 1 → glyph,
  2 → `!`+complement-glyph. A full-domain (`Undetermined`) relation is **vacuous** — it follows the same
  convention as the atom `#a*`/`#T*` special forms: stored on parse (the collection keeps vacuous entries, per
  C3g.5) but **elided** from the canonical render (§7.1), so `#o*`/`#g*` parse but render as the bare element,
  equivalent to omitting the predicate. The element render loop skips `is_undetermined()` constraints.

  **`~` sugar = `kind.involution()`, eager, bidirectional.** It is the key binary-kind assertion: its
  Π-membership *is* the stereogenic-vs-symmetric bit. `#p~`/`#p!~` → `LigandSymmetryAst { perm: involution,
  orientation: Improper if chiral else Proper, mem }` (the orientation is fixed by the kind, so `~` subsumes
  `'`); `#f~` → `FluxionalityAst { perm: involution }` (no `!`; for chiral kinds the involution is improper →
  degenerate, C4 flags). Render emits `~` when a `#p`/`#f` perm equals the involution, else explicit cycles.

  Worked cases — **Th** (chiral, 4 ligands, 2 enantiomeric cosets): `#p~`/`!~` = prochiral vs stereocenter
  (higher symmetry → explicit proper cycles); `#f` rare (`~` degenerate); `#o`/`#g` all three glyphs. **Ct**
  (achiral, 2 diastereomeric configs): `#p~`/`!~` = symmetric vs stereogenic; `#f~` = fluxional E/Z
  interconversion (involution is proper here); `#o`/`#g` binary (`=`/`/`; `'` is chiral-only → C4).

  `'` on achiral kinds (and any glyph/kind/degree consistency) is accepted at parse and rejected by the C4
  validator, consistent with C3h.2. Implemented via a `stereo_predicate_parser!`/`stereo_constraint_fmt!`
  macro pair (atom+bond) over shared helpers (`cycle`/`perm_cycles`, `ligand_pair`, the glyph/relation
  codecs); `Permutation::involution` is now `pub(crate)`. Verified by `test_stereo_atom_inline_roundtrip`
  (10 cases incl. `~`/`!~`/`*`/`!glyph`/multiple), `test_stereo_atom_predicate` (+`_involution`) for parse
  semantics, and `test_stereo_bond_inline_roundtrip`.
- **C3h.4 — roundtrip tests** for both surfaces (structured EDN and inline string). **Done.** Full-molecule
  integration added to `dsl/molecule/tests.rs`: `test_molecule_dsl_stereo_edn_roundtrip` gains inline
  (`:type "Th1#f(0,1,2)#g/"`, `Ct1#g/`) and molecule-scope (`:constraints [{:stereo-atom …}]`,
  `{:stereo-bond …}`) cases for both atom and bond; `test_molecule_dsl_from_edn_str_matches_from_edn` gains the
  same (streaming-vs-tree agreement, exercising the bridge reader). This surfaced a real pre-existing bug:
  `stereo_atom_keyword_for`/`stereo_bond_keyword_for` rendered the canonical `:ccw`/`:cw`/`:z`/`:e` keyword
  shortcut from `(kind, coset)` alone, **silently dropping inline constraints** — fixed to require an empty
  constraint set before using the shortcut. Full suite 3245 green; workspace builds; clippy clean.

The per-element constraints serialize **inline in the stereo-string**, exactly as `dsl/atom.rs` does for the
atom-string: the `StereoAtomDsl` round-trips through `Edn::Str(to_string())`; the `:type` payload carries
`class coset` followed by the `#p`/`#f`/`#o`/`#g` predicates. Replicate the atom.rs scaffolding:

```rust
// mirror `AtomPredicate`:
enum StereoPredicate { Coset(StereoCosetAst), Constraint(StereoAtomConstraint) }

// parse (mirror parse_atom): class → coset → zero+ predicates → StereoAtomAst { coset, constraints }
fn parse_stereo_atom(input) -> StereoAtomDsl;       // extends today's class+coset parse with the predicates
// render (mirror AtomDsl::Display + fmt_constraint loop):
impl Display for StereoAtomDsl {                     // fmt_stereo_coset, then for c in constraints { fmt_stereo_constraint(f, c) }
    …
}
fn stereo_constraint_tag(kind: StereoAtomConstraintKind) -> &'static str;  // "#p" | "#f" | "#o" | "#g"
fn fmt_stereo_constraint(f, c: &StereoAtomConstraint) -> fmt::Result;      // mirror fmt_constraint, per variant
impl ToEdn for StereoAtomDsl { fn to_edn(&self) -> Edn { Edn::Str(self.to_string().into()) } }  // already the shape
// StereoBond* parallel.

// payload parse/render helpers:
//   permutation: GAP cycle string `(0,1,2)(3,4)` ↔ Permutation::{from_cycles(degree=kind.degree(), …), cycles(), Display};
//                prefixes `[!][']`; `~` → kind involution (eager); `^` migrates to cycle notation.
//   relation:    `[!](=|'|/)pair` (#o) / `[!](=|'|/)` (#g) via the DSL-only glyph table:
//                  Topicity:       = Homotopic   ' Enantiotopic / Diastereotopic
//                  Stereogenicity: = Symmetric   ' Prochiral    / Stereogenic
//                (one small per-type fn each way; the core AST carries no glyphs)
```

**Molecule-scope structured form** is the `dsl/constraint.rs` peer (parallel to `:atom` / `atom-constraint-form`,
not `dsl/stereo.rs`): a new `:stereo-atom` / `:stereo-bond` entity-constraint key. There the values are
structured EDN — `PermutationAst` = vector-of-cycles `[[0 1 2] [3 4]]` (identity `[]`), `OrientedPermutationAst`
= `{:perm <vov> :orientation :improper}` (default `:proper`), `LigandPairAst` = `[i j]`, `*RelationAst` =
`:homotopic | #{:enantiotopic :diastereotopic} | :undetermined`, `LigandSymmetryAst` = `{:perm <vov>
:member :not-in}` — paralleling how `coset-form` already admits `[int+]`. (This is the same vector-of-cycles
encoding the molecule-level symmetry assertions reuse — doc 110.) `lift_constraints`/`inline_constraints`
move per-element constraints between the inline and molecule-scope forms, exactly as for atoms.

### C3i · `umol-ast/spec/umol-dsl-spec.md` (update — normative surface) **Done**

Brought the spec in line with C3a–C3h. Edits made: **§7.14** — `stereo-string ::= class coset
stereo-predicate*` plus the predicate subgrammar (`#p`/`#f`/`#o`/`#g`; disjoint-cycle notation; `[!]`/`[']`;
`=`/`'`/`/` glyphs; `~` involution sugar; `*` = vacuous, **parse-admissible / render-elided** per §7.1, matching
the atom `#a*`/`#T*` special-forms convention; chiral-class restriction on `'` deferred to the validator).
**§7.9** — added `:stereo-atom`/`:stereo-bond` to `entity-constraint`; defined `stereo-{atom,bond}-constraint-form`
(`:kind` + one predicate key) with the structured encodings (`permutation-form` vector-of-cycles, `ligand-pair`,
`ligand-symmetry-form`, relation keyword/`#{…}`/`:undetermined`); revised the "stereo elements carry no entity
constraint" paragraph; updated inline-form coverage. **§6.1** — added the symmetry-derived stereo predicates
(ligand symmetry, topicity, stereogenicity) as derived predicates (filter, no grounding, validator
cross-check; `#f` is a stored dynamical assertion). The vacuous-`*` alignment also fixed the inline render to
elide `Undetermined` predicates (consistent with atoms); tests split into `_render` (elision) /
`_render_identity` (roundtrip) per the test-writing identity pattern. (Original section-edit plan below.)

Sections to edit:

- **§7.14 (stereo-string):** extend the `:type` grammar past `class coset` with the inline predicates
  `#p` (ligand-symmetry), `#f` (fluxionality), `#o` (topicity), `#g` (stereogenicity), and define their
  payload subgrammar: GAP cycle notation `(i,j,k)(l,m)` (0-indexed, identity `()`), the `[!]` (negation) and
  `[']` (improper) prefixes for `#p`/`#f`, the `=`/`'`//' relation glyphs for `#o`/`#g`, and `~` (the
  class involution sugar). Note `'` is invalid on achiral kinds.
- **§7.9 (constraint grammar):** **revise** the "stereo elements carry no entity constraint" statement —
  add `:stereo-atom` / `:stereo-bond` to `entity-constraint`, and define `stereo-atom-constraint-form` /
  `stereo-bond-constraint-form` (`{:ligand-symmetry …} | {:fluxionality …} | {:topicity …} |
  {:stereogenicity …}`) with the structured value encodings (permutation = vector-of-cycles `[[i j]…]`,
  oriented-permutation map, ligand-pair `[i j]`, relation keyword / `#{…}` set / `:undetermined`). Add
  `:stereo-atom`/`:stereo-bond` to the **inline-form coverage** and **lift/inline** lists.
- **§6.1 (inherent vs derived):** add the **symmetry-derived** stereo predicates (ligand symmetry,
  topicity, stereogenicity) as derived predicates — filter matches, don't affect grounding, cross-checked
  against the stored constraints by the validator (like the topology-derived/stored field pairs).
- **§7.2 / ids:** note `StereoLigandId` (the 0-based ligand-frame position, `u8`, ≤ kind degree).
- **shared primitive:** document the **vector-of-cycles** permutation EDN once (reused by the molecule-level
  symmetry surface, doc 110).

(Spec edits are part of the implementation — recorded here, not made now.)

### C3j · add to property testing and fuzz corpora **Done**

**Prereq (MemOp unification):** the stereo ligand-symmetry `MemOp` duplicated `value::MemOp` (byte-identical
`{In, NotIn}`) and the stereo copy wasn't publicly reachable (clash at the crate root), blocking randomized
`#p`. Resolved by extracting the shared operator enums (`ArithOp`, `RelOp`, `MemOp`) into a new foundation
module `ast/operators.rs` (re-exported at `umol_ast::ast::{ArithOp, MemOp, RelOp}`); `value.rs` and the stereo
constraints both import from it. One `MemOp` now, publicly reachable.

**Property tests** (`tests/property.rs`): added `permutation_strategy(degree)` (shuffled one-line image),
`ligand_pair_strategy`, `topicity`/`stereogenicity` relation strategies (non-vacuous only — `Undetermined`
elides on render), `orientation`/`mem_op` strategies, and a `stereo_constraint_strategy!` macro generating
`StereoAtomConstraint`/`StereoBondConstraint` (all four kinds, perm degree = `kind.degree()`).
`stereo_atom_ast_strategy`/`stereo_bond_ast_strategy` now attach inline constraints, so the existing
Display↔FromStr and molecule EDN roundtrip props cover the **inline shorthand**. Two new props
(`test_stereo_{atom,bond}_constraint_dsl_to_edn_from_edn_roundtrip`) cover the **EDN-shaped** molecule-scope
form. 42 props green.

**Fuzz corpora** (`fuzz/seeds/`): added `fuzz_entity_strings/stereo_{atom,bond}_predicates` (inline
`#p`/`#f`/`#o`/`#g`), `fuzz_constraints/stereo_atom_{ligand_symmetry,fluxionality,topicity}` +
`stereo_bond_stereogenicity` (EDN entity-constraint forms), and `fuzz_molecule/stereo_atom_{inline_predicates,
entity_constraint}` (both surfaces in a full molecule).

---

## C4 · umol-graph

The **symmetry computation** (`StereoSymmetry`) lives in umol-ast (C3d) — a graph-symmetry algorithm, **not**
a chemistry-perception engine, so it is deliberately *not* a umol-graph `…Perception` struct paralleling
`AromaticityPerception`. umol-graph's stereo **ops** (model/resolver/validator) do mirror the aromaticity ops
(`model.rs` + `resolver/aromaticity.rs` + `validator/aromaticity.rs`), but they **consume** umol-ast's
`StereoSymmetry` rather than housing a perception engine. **Resolver and validator are independent ops** —
and there is **no `StereoSymmetry` to share between them**:

- the **resolver is structural** (C4b) — kind/coset/ligands from scope + topology — and computes **no**
  `StereoSymmetry`;
- the **validator** is the **only** pipeline consumer — it computes `MoleculeAst::graph_symmetry` **once**,
  molecule-wide, on the *resolved* AST, then projects per element (cheap, `image_under`, doc 110).

Sharing would be unsound anyway: the resolver **mutates** the molecule, so any symmetry computed before/during
resolution is stale for the validator (which needs the final ground state). The reuse that *does* matter —
many `topicity`/`stereogenicity`/`ligand_symmetry` queries on one resolved molecule — is served by a cache on
the **immutable `Molecule`** (doc 086, keyed by config), not by op-to-op passing.

### C4a · `ops/model.rs` (extend) — `StereoModel` **Done**

Folded into `ChemistryModel` beside `valence`/`aromaticity`; `ElementScope` (existing) is reused.

```rust
struct ChemistryModel { valence: ValenceModel, aromaticity: AromaticityModel, stereo: StereoModel }

struct StereoModel {
    kind_models: [Option<StereoKindModel>; StereoKind::COUNT],  // slotmap by kind discriminant; None = off
    para_stereo: bool,
    max_iterations: usize,
    inconsistency: InconsistencyPolicy,
}
struct StereoKindModel { scope: ElementScope, fluxionality: bool }   // per-kind: which elements, fluxional?
enum InconsistencyPolicy { Keep, Strip, Error }                     // resolver coverage policy (C4b)
```

Design notes vs the original sketch: `kind_scopes`/`KindScope`/`FluxionalitySetting` collapsed into the
per-kind-slot array `[Option<StereoKindModel>; StereoKind::COUNT]` (`StereoKind` gained `strum::EnumCount`;
indexed by `kind as usize`), with `fluxionality` a per-kind `bool` (not a global setting). There is no
`StereoSymmetryConfig` in umol-ast — `StereoModel::graph_symmetry_config()` builds the existing
`GraphSymmetryConfig<ConstitutionColoring> { coloring: full(), iterate_to_fixpoint: para_stereo, max_iterations }`
(the doc's `para_stereo` *is* `iterate_to_fixpoint` — the fixpoint pass resolves para-stereocenters).
**Default**: tetrahedral + cis-trans on (`ElementScope::Any`), the higher geometries off, `para_stereo: false`,
`max_iterations: 16`, `inconsistency: Error`. `strum` added to umol-graph (trait only). Tests:
`test_stereo_model_default`, `_kind_model` (per-kind table), `_graph_symmetry_config`.

### C4b · `ops/resolver/stereo.rs` (new) — `StereoResolver`

Structural resolver, mirroring `AromaticityResolver` (thin, `Solution`-typed): from each atom `#T` / bond `#C`
it adds a `:stereo-atom` / `:stereo-bond` element with the canonical ligand frame and the coset copied
**verbatim**. Computes no `StereoSymmetry`.

**Changes:**
- New `StereoResolver { model: StereoModel }` with `new(&StereoModel)` and
  `resolve(&self, &mut MoleculeAst) -> Result<Solution<(), StereoContradiction>, StereoError>`.
- Wire into `ops/resolver.rs` **directly after `aromaticity`**: the `Resolver` field, `new()`, the `resolve()`
  call (after aromaticity, before bonds), and a `Stereo` arm in `ResolverContradiction` / `ResolverError`.
- **No coset conversion.** raise (`tetrahedral_ligand_ordering` / `raise_cis_trans_stereo`) already stored the
  coset in the canonical ligand frame; the resolver reproduces that frame and copies the coset — no
  permutation.

**Deferred (later pass):** `InconsistencyPolicy` (Keep/Strip/Error) and the `StereoContradiction` variants it
produces. For now a `#T`/`#C` that can't be realized is simply **not added** (`None`); the resolver never
mutates atom constraints. When that policy lands it is shared with `AromaticityResolver`, not stereo-specific.

```rust
fn resolve_atom(&self, ast: &MoleculeAst, id: AtomId) -> Option<Edit> {
    if ast.stereo_atoms().has_coincident(id) { return None; }   // already a stereo center here → skip
    let atom = ast.atom(id);
    if atom.is_in_aromatic_system() { return None; }            // aromatic atom → not tetrahedral

    let kind = StereoKind::Tetrahedral;
    let StereoConfigurationAst::Stereo(coset) = atom.ast.constraints.tetrahedral_stereo()
        else { return None; };                                  // Undetermined / NotStereo (#T!) → skip
    let model = self.model.kind_model(kind)?;                    // kind off → skip
    if !model.scope.contains(atom.element()?) { return None; }  // out of scope → skip

    let mut ligands: Vec<(AtomRef, StereoLigandKind)> =
        atom.neighbors().map(|n| (AtomRef::Id(n.id()), StereoLigandKind::Atom)).collect();
    if ligands.len() + 1 == kind.degree() {                     // 3 real → one virtual, appended last
        let ValueAst::Lit(h)  = *atom.implicit_hydrogens() else { return None; };
        let ValueAst::Lit(lp) = *atom.lone_pairs()         else { return None; };
        ligands.push((AtomRef::Id(id),
            if h >= 1 { ImplicitHydrogen } else if lp >= 1 { LonePair } else { return None }));
    }
    if ligands.len() != kind.degree() { return None; }          // arity off → skip

    Some(Edit::AddStereoAtom { site: AtomRef::Id(id), ligands, ast: StereoAtomAst::new(kind, coset) })
}
```

`resolve` collects the `Some`s over all atom ids, then all bond ids, `transact`s them, returns
`Solution::Determined(())`. `resolve_bond` mirrors `resolve_atom`: `cis_trans_stereo()`, `kind = CisTrans`,
`site: BondRef`, and the 4-ligand frame reproducing `raise_cis_trans_stereo`'s order
(`[side_1.first, side_1.second, side_2.first, side_2.second]`).

**Substeps:**
- **C4b.1** **Done** — `StereoResolver` skeleton + `resolve` loop + wiring into `Resolver` (field / `new` /
  `resolve` / enums, after aromaticity). `resolve_atom`/`resolve_bond` stubbed to `None`; `StereoContradiction`/
  `StereoError` are empty (mirror `BondsContradiction`/`BondsError`). `model` field unread until C4b.2.
- **C4b.2** **Done** — `resolve_atom` (`#T` → `Tetrahedral`): guards (coincident, aromatic, `#T` config,
  kind/scope via `as_lit`), canonical frame (neighbors + one virtual last), arity check, coset verbatim.
- **C4b.3** — `resolve_bond` (`#C` → `CisTrans`): mirror with the cis/trans 4-ligand frame.
- **C4b.4** — tests: atom + bond add, idempotency skip, aromatic skip, out-of-scope / arity skip.

### C4c · `ops/validator/stereo.rs` (new) + `ops/validator.rs` (wire) — `StereoValidator` (mirror `ops/validator/aromaticity.rs`)

Carries the `StereoModel` and computes the `StereoSymmetry` **once** (molecule-wide, on the resolved AST),
projecting per element — like `AromaticityValidator` carries its model; runs **last**, after
`AromaticityValidator`. (It is the **only** pipeline consumer of `StereoSymmetry`; the resolver is structural.) Pass/fail, never mutates (doc 092: validators are tier-2
+ constraint↔data agreement; tier-1 structural is enforced at `MoleculeAst::new`; always-on). Invariants:

1. **Coset in range** — `0 ≤ index < kind.count()` for the resolved kind.
2. **Ligand-frame arity = kind degree** (port count).
3. **Achiral-kind gate** — improper `'` relations/literals (`#o'`/`#g'`) are invalid on achiral kinds (no
   enantiomers).
4. **Asserted-constraint ↔ derived agreement** (tier-2 constraint↔data, docs 087/092): every ground asserted
   `#p`/`#o`/`#f` / `#g` is verified against the validator's `StereoSymmetry`; mismatch = contradiction. This
   is the derived-vs-stored cross-check — the stereo analog of the topology-derived/stored field pairs
   (doc 086) that the constraint validator already cross-checks.

**Not** validator errors: a concrete coset on a non-stereogenic site (legitimate labeled data, doc 104 — the
resolver's `InconsistencyPolicy` governs redundancy, not the validator); CIP/chirality-label invention
(doc 080 — deferred, out of scope). No transformer at this stage.

### C4d · `ops/resolver/aromaticity.rs` (fix — idempotency)

`AromaticityResolver::resolve` calls `find_systems` + `add_systems` unconditionally, so re-running an
already-aromatized molecule **duplicates** systems. Skip atoms already in an aromatic system (the membership
query), making re-runs a stable no-op — the analog of the stereo resolver's `has_coincident` guard (C4b). The
inconsistency policy that previously also lived here is deferred together with C4b's.

### C4e · fuzz corpora + proptests

- **Fuzz** (`umol-ast/fuzz`): extend the DSL/EDN round-trip target to cover stereo-strings carrying
  `#p`/`#f`/`#o`/`#g` and the `:stereo-atom`/`:stereo-bond` molecule-scope forms; seed the `corpus` with
  stereo examples (each kind; the four constraints; `~`/`'`/`!` variants; cycle products; vector-of-cycles
  EDN). (`umol-edn`/`umol-io` fuzzers unaffected.)
- **Proptests** (per crate, where the type lives):
  - **umol-perm:** `Permutation` cycle round-trip (`from_cycles(cycles()) == self`; `Display` parse-back),
    `compose`/`inverse`/`sign` laws, `OrientedPermutationGroup` closure, and `CosetSpace::{merge_under,
    observable_coset, enantiomer, is_chiral}` invariants.
  - **umol-ast:** `StereoAtomDsl`/`StereoBondDsl` parse↔render and EDN round-trip; `StereoRelationAst` and
    `StereoAtomConstraints` `Lattice` laws (meet commutative/associative/absorptive, `matches` ⇔
    `meet == Some(target)`); `LigandPairAst` normalization; the concrete-literal inherent `matches`.

### C4f · conformance suite (large chunk)

Extend the **feature-gated** conformance suite (`umol-graph`, `--features conformance`, the resolution test)
with stereo — the bulk of Phase C's verification:

- **OpenSMILES arrangements** (docs 038/047): the `@TH`/`@SP`/`@TB`/`@OH` numbering cases (the `class.rs`
  reindex pairs, scaled to the full corpus) through parse → resolve → coset.
- **Constraints:** `#p`/`#o`/`#g`/`#f` parse + resolve + validate on reference molecules; topicity /
  stereogenicity on known cases (CH₂Cl₂, ethanol's prochiral CH₂, 1,2- vs 1,1-dichloroethene, allene);
  the achiral `'`-gate rejection; resolver idempotency (re-resolve is a no-op); `InconsistencyPolicy`
  outcomes.
- **Highly-symmetric integration** (doc 110): cubane, o-carborane `B10C2H12`, Cu-phthalocyanine, `C70` —
  orbit/topicity/symmetry-number correctness.
- **End-to-end round-trip:** parse (`F[C@H](Cl)Br`, `F/C=C/F`, …) → resolve → serialize `:stereo-*` →
  reparse → equal.
