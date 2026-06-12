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

## C2 · umol-graph-core

### C2a · `algorithms/auto.rs` (extend `Automorphism`)

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

### C3a · `ast/coloring.rs` (new) + `Entity` in `ast/molecule.rs` (extend)

The pluggable **round-0 color policy**: one `u64` per **graph node** — atom or pseudonode (a relation:
localized bond or an overlay) — that says what counts as distinguishable. One rule, applied to every node:
**hash the entity's inherent fields**. This is *complete* for our consumer (a graph automorphism): every
**derived** predicate is relational, so it is a function of the structure and is preserved by the
automorphism for free — once each overlay is an **entity pseudonode** (C3c), even overlay-derived predicates
(aromaticity, dative/multicenter participation) are graph functions. So the color carries only the
**irreducible inherent labels**; nothing derived enters it. Named `MoleculeColoring` to avoid conflation
with graph-core's proper-coloring `algorithms/coloring.rs`.

Because the rule is uniform across node kinds, the trait is a **single method over a public `Entity`** (a
typed ref to any molecule entity — lives in `ast/molecule.rs`, generally useful beyond coloring):

```rust
// ast/molecule.rs  (public, general-purpose)
enum Entity {
    Atom(AtomId), Bond(BondId), Dative(DativeBondId), Aromatic(AromaticSystemId),
    Multicenter(MulticenterBondId), Noncovalent(NoncovalentBondId),
    StereoAtom(StereoAtomId), StereoBond(StereoBondId),
}

// ast/coloring.rs
trait MoleculeColoring {
    fn color(&self, mol: &MoleculeAst, entity: Entity) -> u64;   // invoked only for graph-participating kinds
}

bitflags! struct ColorFeatures: u16 {              // inherent fields only (§6.1) — derived excluded (free)
    ELEMENT, ISOTOPE, CHARGE, IMPLICIT_H, LONE_PAIRS, SPIN,   // atom
    BOND_ORDER, BOND_CHARGE, BOND_SPIN,                       // localized bond
}

struct ConstitutionColoring { features: ColorFeatures }   // the one impl now (pre-stereo inherent fields)
impl ConstitutionColoring { fn new(features: ColorFeatures) -> Self; fn full() -> Self; }
impl MoleculeColoring for ConstitutionColoring;
// color(): match entity kind → hash (kind tag, that kind's inherent fields):
//   Atom → ColorFeatures-selected atom fields;  Bond → order/charge/spin;
//   Aromatic → charge/spin/π-count;  Multicenter → charge/spin/e-count;
//   Dative → order (direction is gadget-encoded, C3c);  Noncovalent → kind.
// the kind tag keeps kinds in disjoint color ranges (no pseudonode maps onto an atom or another kind).
```

`ConstitutionColoring` serves **our** distinctness / symmetry analysis. `MorganColoring`/`EcfpColoring`/`FcfpColoring`
are **separate** impls owning their own fields (ring, aromaticity, functional class) — they carry derived
props because their consumer is per-atom hashing, not an automorphism, so "relational ⇒ free" doesn't apply
to them. Stereo is **not** a coloring impl: its descriptor is partition-dependent, so the symmetry-refinement loop
folds it on top (C3d); the trait stays stereo-free (and fingerprint-reusable). `Entity`'s `StereoAtom`/
`StereoBond` variants exist for the general enum but are never passed to `color` (stereo isn't a graph node).

### C3b · `ast/automorphism.rs` (extend)

Expose generators on `AtomAutomorphism`; add the orientation-graded holder (Â proper + Â* mirror).

```rust
impl AtomAutomorphism {
    fn generators(&self) -> Vec<Vec<AtomId>>;        // graph-core generators as AtomId maps
}

struct GradedAutomorphism {
    proper: AtomAutomorphism,                         // Â  (orientation-preserving)
    improper: Option<Vec<Vec<AtomId>>>,              // Â* (molecule→mirror mapping generators)
}
impl GradedAutomorphism {
    fn proper(&self) -> &AtomAutomorphism;
    fn improper_generators(&self) -> Option<&[Vec<AtomId>]>;
}
```

### C3c · `ast/incidence.rs` (new) — generic

The molecule's **incidence graph** (Levi graph): every relation — localized bond and each overlay
(aromatic / multicenter / dative / noncovalent) — lifted to a node alongside the atoms, with incidence edges
to its participants. A **generic** `MoleculeAst` property (canonicalization, fingerprints, symmetry numbers
all reuse it), not stereo-specific; structure only — colors are applied at automorphism time. Hyperedge
overlays (>2 atoms) are why it's the incidence/Levi graph rather than an ordinary-edge subdivision; it is
**not** a port graph — attachment stays unordered (ligand order enters only in the per-site projection, C3d).

```rust
struct IncidenceGraph {
    graph: Graph,             // CSR (graph-core); nodes = atoms ++ relation-pseudonodes
    node_entity: Vec<Entity>, // node → the molecule Entity it represents
}
impl MoleculeAst {
    fn incidence_graph(&self) -> IncidenceGraph;   // dative donor→acceptor direction = distinguished gadget edges
}
```

### C3d · `ast/stereo_symmetry.rs` (new)

`stereo_symmetry` runs the automorphism over the incidence graph (C3c) with colors folded each round;
**`StereoSymmetry`** is the *result* (the converged graded automorphism + partition + queries) — the noun
names the result, not the process. Molecule-level orbit queries + the per-carrier projection
(the `image_under` bridge to a small local group).

```rust
struct StereoSymmetryConfig<C: MoleculeColoring> { coloring: C, para_stereo: bool, max_iterations: usize }

struct StereoSymmetry { automorphism: GradedAutomorphism /* + converged partition */ }

impl MoleculeAst {
    fn stereo_symmetry<C: MoleculeColoring>(&self, cfg: &StereoSymmetryConfig<C>) -> StereoSymmetry;
    // inc = self.incidence_graph();
    // loop: color(node) = cfg.coloring.color(self, inc.node_entity[node]) ⊕ stereo_descriptor(partition);
    //   auto = inc.graph.automorphisms(color); converge when partition stable or !para (one pass).
    //   Â* from a molecule↔mirror isomorphism.
}

impl StereoSymmetry {
    // molecule-level
    fn same_proper_orbit(&self, a: AtomId, b: AtomId) -> bool;
    fn same_star_orbit(&self, a: AtomId, b: AtomId) -> bool;
    fn proper_orbit_of(&self, a: AtomId) -> Vec<AtomId>;
    fn star_orbit_of(&self, a: AtomId) -> Vec<AtomId>;
    fn graded(&self) -> &GradedAutomorphism;
    // per-carrier projection: site stabilizer (targeted re-run, now) → local ligand-position group.
    // Later a BSGS stabilizer feeds the same projection — the API is unchanged (doc 110).
    fn ligand_symmetry(&self, site: AtomId, ligands: &[StereoLigand]) -> OrientedPermutationGroup;
}
```

### C3e · `ast/stereo.rs` (extend)

Add `Axial`; delegate chirality to umol-perm; add `MirrorOp` (`'`); add the two ground relation enums
(no trait — glyph mapping lives in the DSL, C3h).

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

### C3f · `ast/views/stereo.rs` (extend); `StereoLigandId` in `ast/ids.rs`

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

### C3h · `dsl/stereo.rs` (extend — inline stereo-string constraints, mirroring `dsl/atom.rs`)

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

### C3i · `umol-ast/spec/umol-dsl-spec.md` (update — normative surface)

Bring the spec in line with C3a–C3h (the planned evolution: stereo elements gain entity constraints).
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
- the **validator** is the **only** pipeline consumer — it computes `MoleculeAst::stereo_symmetry` **once**,
  molecule-wide, on the *resolved* AST, then projects per element (cheap, `image_under`, doc 110).

Sharing would be unsound anyway: the resolver **mutates** the molecule, so any symmetry computed before/during
resolution is stale for the validator (which needs the final ground state). The reuse that *does* matter —
many `topicity`/`stereogenicity`/`ligand_symmetry` queries on one resolved molecule — is served by a cache on
the **immutable `Molecule`** (doc 086, keyed by config), not by op-to-op passing.

### C4a · `ops/model.rs` (extend) — `StereoModel` + `KindScope`

Folded into `ChemistryModel` beside `valence`/`aromaticity`; mirrors `AromaticityModel` / `ElementScope`.

```rust
struct ChemistryModel { valence: ValenceModel, aromaticity: AromaticityModel, stereo: StereoModel }

struct StereoModel { kind_scopes: KindScopes, fluxionality: FluxionalitySetting, para_stereo: bool, inconsistency: InconsistencyPolicy }
enum KindScope { /* per-kind: element + coordination/bond-order scope, mirror ElementScope */ }
// StereoModel produces the umol-ast `StereoSymmetryConfig { coloring, para_stereo, max_iterations }`.
```

### C4b · `ops/resolver/stereo.rs` (new) — `StereoResolver` (mirror `ops/resolver/aromaticity.rs`)

Unlike valence (which narrows atoms in place), this **adds and resolves new entities**: from each atom `#T` /
bond `#C` it builds a `:stereo-atoms` / `:stereo-bonds` element — kind by `KindScope` + coordination/order,
ligand frame from neighbors + virtuals, and the **coset resolved** against the kind's `CosetSpace`. Stereo
elements are **not** ground-by-default (unlike aromatic systems) — they need this resolution.

- **Idempotency is site-membership, not `is_ground`.** The valence `is_ground` skip is the wrong test here
  (valence narrows atoms; this adds entities). Re-running must not re-add — skip sites that already bear a
  stereo element, via the existing site query (`StereoAtomViews::has_coincident(atom)` / `coincident`;
  `StereoBondViews` parallel). (Aromaticity's analog: skip atoms already in an aromatic system — that
  membership query exists.)
- **Coverage + inconsistency policy (configurable, never silent).** When a `#T`/`#C` assertion can't be
  (fully) realized, a policy `{ Keep | Strip | Error }` decides — never a silent drop. It covers both
  (i) **partial coverage** (some asserting sites get an element, others can't) and (ii) an **unsatisfiable
  assertion** (a site that can't bear a stereo element of the scoped kind). All of this is **structural**
  (scope / coordination / topology), so the **resolver computes no `StereoSymmetry`** — non-stereogenicity is
  *not* an inconsistency here (a coset on a non-stereogenic site is legitimate labeled data, doc 104).
- The **same two gaps affect `AromaticityResolver`** — these are general entity-adding-resolver concerns,
  not stereo, so they're pulled out into **C4d**.

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

### C4d · `ops/resolver/aromaticity.rs` (fix — non-stereo, surfaced here)

`AromaticityResolver` has the same two entity-adding-resolver gaps the stereo resolver must avoid (C4b).
They're **general, not stereo-specific** — fixed here so both resolvers behave consistently:

- **Idempotency = membership, not re-perception.** Today `resolve` calls `find_systems` + `add_systems`
  unconditionally, so re-running an already-aromatized molecule **duplicates** systems. Skip atoms already in
  an aromatic system (the membership query), making re-runs a stable no-op.
- **Coverage + inconsistency policy** (`{ Keep | Strip | Error }`, never silent) — replaces today's silent
  discard of failing proposals (doc 080). Handles **partial coverage** (some but not all `#a` atoms land in
  systems) and an **unsatisfiable assertion** (a lone `#a+` atom that can't form a ring system).

The policy type is shared with the stereo resolver (one `InconsistencyPolicy`, model-layer); this is the
principled home for it rather than duplicating per resolver.

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
