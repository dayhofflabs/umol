# Molecular permutation analysis and scalable permutation groups

Status: **Decision / staged roadmap.** Genesis: surfaced while designing the stereo ligand-symmetry
constraints (doc 104 §C). Builds on the permutation/group infrastructure roadmap in
[109-permutation-infrastructure-2026-06-09.md](109-permutation-infrastructure-2026-06-09.md).

## Decision summary

- Build the **stereo ligand-symmetry, topicity, and stereogenicity queries first**, without first
  implementing molecule-scale BSGS machinery.
- The stereo implementation must reason over the **whole molecule**, because remote stereochemistry
  can distinguish otherwise-isomorphic ligands. Its result for one stereo element is nevertheless a
  **small local image group** acting on the stereo kind's bounded ligand set — at most six positions
  for the current classes — which should remain fully enumerated in `umol-perm`.
- Put the chemistry-aware query engine in **`umol-ast`**, as an explicit
  `MolecularAutomorphismAnalysis` result plus stereo-view projection methods. Put representation-
  independent permutation-group operations in **`umol-perm`**. Keep graph-automorphism discovery in
  **`umol-graph-core`**.
- Initially obtain site stabilizers using targeted exact graph-automorphism runs with the site
  distinguished. Later, a molecule-scale generated group + BSGS can provide the same stabilizers
  more efficiently. This is an internal backend change: the stereo element and query APIs do not
  change.
- Add molecule-scale generated groups and BSGS when molecule-level consumers need repeated
  stabilizers, membership, or symmetry-aware search. Reaction mapping, canonical-image search,
  symmetry-pruned enumeration, and symmetry-aware RMSD are likely early consumers.
- Keep the split between **enumerated small groups** and **generated large groups**. The two
  representations have different semantics and operations, not merely different performance.

The sequencing decision is therefore doc-104's option **(a)**, but with an explicit reusable
molecule-analysis boundary. The progression to a top-level permutation-analysis facility is
incremental rather than a later rewrite.

## Scope and terminology

The core object available from the molecular graph is the automorphism group of a **chosen colored
graph representation**. The coloring policy determines which atom, bond, relation, and stereo facts
must be preserved. It is therefore better called a **molecular automorphism/permutation analysis**
than an unqualified molecular-symmetry group.

Three related objects must not be conflated:

1. **Colored constitutional automorphisms** preserve the selected atom/bond/relation features.
2. **Orientation-preserving molecular mappings** preserve the relevant oriented stereochemical
   frames under the selected stereo/perception model.
3. **Orientation-reversing / molecule-to-mirror mappings** reverse those frames consistently. These
   supply the mappings needed for enantiotopicity and molecular achirality tests.

A plain graph automorphism has no intrinsic proper/improper grade. That grade must come from how the
mapping acts on oriented stereochemical frames, or from an explicit molecule-to-mirror isomorphism.
When no stored stereo element supplies the needed frame — for example, a queried prochiral candidate
site — the analysis must materialize the candidate frame transiently or use an equivalent
replacement/mirror test. Defining and validating this grading is a chemistry-aware analysis task,
not a BSGS operation.

This graph-derived analysis is useful for atom/bond equivalence, topicity, regiochemical site
reduction, canonicalization, and symmetry-pruned enumeration. It is **not by itself** the
Longuet-Higgins feasible permutation-inversion group: connectivity alone does not establish which
internal motions are dynamically feasible. It also does not alone establish NMR magnetic
equivalence or thermochemical symmetry numbers; those require additional physical/geometric input.
The geometric point group remains `umol-msym`'s domain.

## Central decomposition

The reusable design has four layers.

### 1. Graph automorphism discovery — `umol-graph-core`

`umol-graph-core` discovers automorphisms of a supplied graph representation. It should remain
chemistry-agnostic and expose ordinary index mappings:

```rust
pub struct Automorphism {
    orbits: Vec<NodeId>,
    canonical_lab: Vec<NodeId>,
    generators: Vec<Vec<NodeId>>,
    group_order: AutoGroupOrder,
}
```

The immediate missing pieces are:

- expose nauty/Traces automorphism generators instead of discarding them;
- support the caller's edge/relation labels, normally by constructing an appropriate colored
  incidence/gadget graph before invoking nauty;
- support distinguished vertices, edges, and relation gadgets for stabilizer queries;
- keep automorphism discovery pluggable independently of the later `umol-perm` group backend.

The bundled `nauty-Traces-sys` already exposes generator callbacks and ships nauty's random Schreier
implementation. The latter is useful internally to nauty, but is probabilistic and too limited to be
the generic exact group API of `umol-perm`.

### 2. Chemistry-aware molecular analysis — `umol-ast`

`umol-ast` owns the query engine because it understands molecular fields, overlays, stereo
descriptors, virtual ligands, and structural carriers. The initial API should return an explicit
analysis result:

```rust
pub struct MolecularAutomorphismConfig<C> {
    pub coloring: C,
    pub para_stereo: bool,
    pub max_iters: usize,
}

pub struct MolecularAutomorphismAnalysis<'a> {
    molecule: &'a MoleculeAst,
    // Orientation-preserving and molecule-to-mirror mapping results.
}

impl MoleculeAst {
    pub fn automorphism_analysis<C: AtomColoring>(
        &self,
        config: &MolecularAutomorphismConfig<C>,
    ) -> MolecularAutomorphismAnalysis<'_>;
}
```

The analysis exposes molecule-level queries such as:

```rust
impl MolecularAutomorphismAnalysis<'_> {
    pub fn same_constitutional_orbit(&self, a: AtomId, b: AtomId) -> bool;
    pub fn same_proper_orbit(&self, a: AtomId, b: AtomId) -> bool;
    pub fn same_star_orbit(&self, a: AtomId, b: AtomId) -> bool;

    pub fn constitutional_atom_orbit(&self, atom: AtomId) -> Vec<AtomId>;
    pub fn proper_atom_orbit(&self, atom: AtomId) -> Vec<AtomId>;
    pub fn star_atom_orbit(&self, atom: AtomId) -> Vec<AtomId>;
}
```

Do **not** initially store this as a ring-style `OnceLock` field on `MoleculeAst`. Unlike topology-
only rings, the result depends on atom/bond attributes and on a caller-selected coloring/stereo
policy. `MoleculeAst` currently permits those attributes to mutate in place, so such a cache could
become stale. A cache can be added later in a service/context keyed by molecule revision and complete
analysis configuration.

There is also no single policy-independent "true symmetry" cache entry. A chemistry-default coloring
can be convenient, but it is still one explicitly named policy.

This refines doc 104's proposed `StereoPerception` artifact. Prefer naming the reusable result
`MolecularAutomorphismAnalysis` now, with the stereo perception/model producing its configuration.
If `StereoPerception` remains as a public name, it should be only a thin stereo-configured wrapper
over the general analysis. The stereo-view API still takes the analysis by reference either way.

### 3. Site projection and stereo queries — `umol-ast`

For a stereo carrier `s`, the needed object is not the full molecule-scale group. It is the image of
the carrier stabilizer on the carrier's ordered ligand positions:

```text
Π_s^+ = image(Stab_{A+}(s) -> S_ligands(s))
Π_s^- = image({molecule-to-mirror mappings stabilizing s} -> S_ligands(s))
```

`A+` is the orientation-preserving mapping subgroup under the selected stereo/perception model.
`Π_s^+` is the proper local subgroup; `Π_s^-`, when non-empty, is its improper coset. Their oriented
union is the local ligand symmetry `Π_s`.

The stereo-facing API remains local and stable:

```rust
pub struct LigandSymmetry {
    group: OrientedPermutationGroup,
}

impl StereoAtomView<'_> {
    pub fn ligand_symmetry(
        &self,
        analysis: &MolecularAutomorphismAnalysis<'_>,
    ) -> LigandSymmetry;

    pub fn topicity(
        &self,
        a: usize,
        b: usize,
        analysis: &MolecularAutomorphismAnalysis<'_>,
    ) -> Topicity;

    pub fn is_stereogenic(
        &self,
        analysis: &MolecularAutomorphismAnalysis<'_>,
    ) -> bool;
}
```

`StereoBondView` has the parallel API. Extended carriers later define their own stabilizer semantics
(pointwise, setwise, ordered path, directed axis, and so on) while returning the same
`LigandSymmetry` type.

Assertions should evaluate against one derived local result rather than repeatedly querying the
molecule:

```rust
let symmetry = stereo.ligand_symmetry(&analysis);

constraint.matches(&symmetry);
symmetry.topicity(0, 1);
symmetry.is_stereogenic(stereo.kind(), stereo.coset());
```

`LigandSymmetry` is a chemistry-semantic wrapper owned by `umol-ast`; its underlying oriented
permutation group and generic group operations are owned by `umol-perm`.

### 4. Generic permutation algebra — `umol-perm`

`umol-perm` owns no graph or chemistry discovery. It consumes permutations and supplies:

- small and dynamic permutation types with the same action/composition convention;
- enumerated and generated groups;
- proper/star orbit and membership operations;
- stabilizers and BSGS when the generated-group layer lands;
- local coset merging and the existing stereo arrangement algebra.

For the stereo queries, the local group degree is bounded by the stereo kind — at most six for the
current classes — and should be enumerated:

```rust
pub struct OrientedPermutationGroup {
    proper: EnumeratedGroup,
    improper_rep: Option<OrientedPermutation>,
}

impl OrientedPermutationGroup {
    pub fn contains(&self, op: OrientedPermutation) -> bool;
    pub fn proper_orbit_of(&self, point: usize) -> Vec<usize>;
    pub fn star_orbit_of(&self, point: usize) -> Vec<usize>;
}
```

`CosetSpace::merge_under` consumes the proper permutations from this group. Thin semantic methods
such as `LigandSymmetry::topicity` and `LigandSymmetry::is_stereogenic` stay in `umol-ast`. The
important boundary is that `umol-perm` knows group actions, while `umol-ast` knows what ligands,
topicity, and stereogenicity mean chemically.

## Computing stereo queries without molecule-scale BSGS

BSGS is not required for doc-104's stereo-specific assertions.

An exact initial implementation can compute `Π_s` as follows:

1. Build the chemistry- and stereo-colored graph representation for the molecule.
2. Distinguish the carrier `s` and rerun exact graph automorphism discovery. For an atom carrier,
   uniquely color its vertex. For a bond or extended carrier, distinguish an appropriate edge or
   relation gadget with the required pointwise/setwise/ordered semantics.
3. Extract generators of the resulting carrier stabilizer.
4. Project each atom-scale generator onto the carrier's ligand-position action, including the
   defined handling of virtual-ligand blocks.
5. Enumerate the generated image group on the stereo kind's bounded position set.
6. Obtain orientation-reversing local mappings from a constrained molecule-to-mirror isomorphism.
   If one improper representative exists, composing it with the proper subgroup gives the improper
   coset.

The local result directly answers:

- positive `g in Π_s` assertions;
- negative `g not in Π_s` assertions;
- proper and star orbits;
- homotopic/enantiotopic/diastereotopic classification;
- prochirality;
- stereogenicity through `CosetSpace::merge_under(Π_s^+)`.

Negative assertions require no complement-group construction: they are ordinary failed membership
tests against the concrete local enumerated group.

This approach may rerun nauty for each queried carrier. That is acceptable for the bounded stereo
deliverable and gives a correctness reference for the later BSGS path. A BSGS-backed analysis can
replace the targeted reruns when profiling or molecule-level consumers justify it.

## Incremental progression

The intended progression is:

1. Build chemistry-correct colored graph representations, generator extraction, and mirror-aware
   analysis.
2. Implement `MolecularAutomorphismAnalysis` and project targeted carrier stabilizers into enumerated
   `LigandSymmetry`.
3. Implement stereo ligand-symmetry literals, topicity, prochirality, and stereogenicity against
   `LigandSymmetry`.
4. Expose molecule-level atom/bond orbit and equivalence queries on the same analysis object.
5. Add arbitrary-degree permutations and generated groups in `umol-perm`.
6. Add BSGS and use it internally for repeated stabilizer/membership queries.
7. Add reaction mapping, canonical-image, double-coset, and other molecule-scale search consumers.

Provided stereo views consume `MolecularAutomorphismAnalysis` and return `LigandSymmetry`, steps 5–7
do not change the stereo-element API.

The doc-104 signed-permutation literal model is also reusable, but the molecule-level pattern
surface is not necessarily identical to the stereo surface: stereo literals address ligand
positions, while molecule-level constraints address atoms, bonds, or other structural entities.
They should share predicate/evaluation machinery without forcing one AST notation onto both.

## Scope of the scalable `umol-perm` extension

The core algorithms are well understood. For the target of roughly 100 atoms, one exact
deterministic incremental Schreier-Sims implementation should be sufficient; a runtime-pluggable
permutation-group backend is not initially justified.

The work should be staged:

| Stage | Capability | Scope |
| --- | --- | --- |
| 1 | `DynPermutation`, generated group, generator validation, point orbits, bounded enumeration | Small/moderate |
| 2 | Exact deterministic BSGS construction, sifting, exact order, membership, pointwise stabilizers | Moderate |
| 3 | Base changes, generator filtering, randomized-and-verified construction, memory/performance tuning | Moderate |
| 4 | Setwise stabilizers, transporters, subgroup intersections, canonical images, double-coset search | Large; separate search algorithms |

BSGS makes membership, exact order, point orbits, and pointwise stabilizers efficient. It does **not**
make every advanced group operation automatic. Setwise stabilizers, intersections, transporters,
canonical images, and double-coset representatives require additional backtracking/search algorithms
and practical heuristics.

The named molecular targets — cubane, `B10C2H12`, copper phthalocyanine, and `C70` — are reasonable
integration tests, but are unlikely to be the hardest BSGS cases. Chemistry-correct graph encoding
and automorphism discovery may dominate their cost. Stress tests should also include explicit
identical hydrogens, many repeated disconnected components, and wreath-product-like generated groups,
whose element counts are enormous despite compact generators.

Keep algorithm selection pluggable where it already matters:

- graph automorphism/isomorphism discovery in `umol-graph-core`;
- later canonical-image/transporter search strategies in `umol-perm`.

Do not introduce a general runtime backend abstraction around the basic BSGS representation until
there is evidence that two implementations with meaningfully different tradeoffs are needed.

## Reference implementations and validation

- **[GAP permutation groups](https://docs.gap-system.org/doc/ref/chap43.html)** are the production
  reference and should be the differential-test oracle for generated groups, orders, membership,
  stabilizers, and advanced operations.
- **[SymPy permutation groups](https://docs.sympy.org/latest/modules/combinatorics/perm_groups.html)**
  have a readable exact incremental and randomized Schreier-Sims implementation suitable for
  studying the algorithm and constructing independent test cases.
- **[nauty/Traces](https://pallini.di.uniroma1.it/)**, already bundled through
  `nauty-Traces-sys`, supplies automorphism generators and contains a random Schreier
  implementation. Its documented results are probabilistic and its API is focused on nauty's
  search, so it should not define the public `umol-perm` semantics.
- A crates.io search found no mature general-purpose Rust permutation-group/BSGS library suitable as
  a dependency. The [`butler-portugal`](https://crates.io/crates/butler-portugal) crate contains a
  small tensor-focused implementation, but it is not a reliable reference implementation for this
  work.

The test strategy should compare generated groups against:

- the existing enumerated groups for every group small enough to enumerate;
- GAP/SymPy results for standard and randomized larger groups;
- graph automorphism group orders and stabilizer orders reported by nauty;
- projection invariants: mapped stabilizer generators must generate exactly the local
  `LigandSymmetry` obtained by the targeted-run reference path.

## Enumerated versus generated groups

Keep the split proposed in doc 109:

```rust
pub struct SmallPermutation { /* degree <= 6, Copy */ }
pub struct DynPermutation { /* arbitrary degree */ }

pub struct EnumeratedGroup { /* every element stored */ }
pub struct GeneratedGroup { /* generators + optional stabilizer chain */ }
```

The existing `CosetSpace` fundamentally depends on enumeration:

- every arrangement has a stable OpenSMILES-compatible index;
- representatives are stored explicitly;
- canonical representatives are selected by testing the small rotation group;
- stereo `merge_under` acts on a small finite arrangement set.

A BSGS scales down computationally, but it does not remove the need for explicit enumeration where
enumeration is part of the API's meaning. Conversely, molecule-scale groups cannot generally expose
`elements()`.

Useful bridges are:

```rust
impl GeneratedGroup {
    pub fn try_enumerate(&self, limit: usize) -> Result<EnumeratedGroup, EnumerationLimit>;

    pub fn image_under<H: PermutationHomomorphism>(
        &self,
        homomorphism: &H,
        limit: usize,
    ) -> Result<EnumeratedGroup, EnumerationLimit>;
}
```

The stereo path is exactly such a bridge: take a potentially huge molecule-scale carrier stabilizer,
map its generators through the carrier-action homomorphism into a tiny ligand-position action, and
enumerate the image.

## `MoleculeAst::canonical_eq` as a generated-group consumer

Full semantic equality for independently numbered `MoleculeAst` values is a concrete consumer of the
scalable generated-group and canonical-image work. The nearer-term comparison API is deliberately
narrower: [doc 156](156-ast-comparison-and-property-suite-2026-07-20.md) defines `MoleculeAst::equiv`
for the current ID frame and `MoleculeAst::equiv_under` for a supplied total
`MoleculeCorrespondence`. Neither operation searches for a correspondence or requires canonical
labeling.

Future `MoleculeAst::canonical_eq` should compare full semantic canonical forms independently of atom,
bond, and overlay IDs. For a graph-plus-overlays implementation, the operation decomposes as follows:

1. canonicalize every entity AST and constraint without changing the molecule ID frame;
2. obtain a canonical labeling and automorphism generators for the ordinary atom-bond topology;
3. compute the canonical image of the complete molecule under that automorphism action, including
   bond ASTs, all overlay/stereo families, participant frames, and molecule constraints;
4. derive the complete `MoleculeCorrespondence` and `IdRemapping` selected by the canonical image;
5. rebuild the molecule in that canonical ID frame.

Step 3 is the generated-group dependency. Topology-only canonical labeling does not determine a
unique complete molecule when bond ASTs, overlays, stereo, or constraints break a graph
automorphism. Canonical-image search selects one complete molecule image without enumerating the
automorphism group. A faithful full-Levi canonical labeling avoids this separate search by encoding
ordinary bonds as vertices, but pays the cost of splitting every bond. The graph-plus-overlays and
full-Levi paths therefore need comparative prototypes and benchmarks before selecting the production
implementation.

This consumer requires the scalable `umol-perm` roadmap through more than BSGS construction:

- `DynPermutation` and `GeneratedGroup`;
- exact stabilizer-chain construction and sifting;
- a canonical-image action and transporter contract;
- symmetry-pruned canonical-image search;
- differential validation against enumerated small groups and GAP/SymPy;
- benchmarks for asymmetric molecules, overlay-broken symmetry, highly symmetric cages, and repeated
  disconnected components.

The exact generated-group representation, canonical-image API, search algorithm, and greenfield versus
external implementation decision remain open. They must be settled in the scalable `umol-perm` design
before this path becomes an implementation plan. `MoleculeAst::canonical_eq` must not make these
choices implicitly.

Once that dependency is available, the AST-side work consists of a molecule-wide remapping operation,
`Canonicalize for MoleculeAst`, canonicalization laws (idempotence and numbering invariance), and
comparison tests distinguishing `==`, `equiv`, `equiv_under`, and `canonical_eq`. No internal
canonical-state flag or cache is part of this design; callers that need a proof-bearing canonical value
use the existing `Canonical<MoleculeAst>` wrapper.

## Placement and materialization decision

The top-level facility is an **analysis service over `MoleculeAst`**, not initially a mandatory stored
property of `MoleculeAst`:

```text
MoleculeAst
    + MolecularAutomorphismConfig
        -> MolecularAutomorphismAnalysis
            -> molecule-level orbit/equivalence queries
            -> per-carrier LigandSymmetry projections
```

This placement preserves the general nature of the facility while respecting that:

- the result is coloring- and stereo-policy-relative;
- attribute mutation currently prevents a simple internal cache from remaining valid;
- callers often need one analysis reused across many stereo queries;
- future BSGS storage is an implementation detail of the analysis result.

`umol-msym` remains complementary: it owns geometric point-group analysis, while this facility owns
colored graph automorphism and permutation-action analysis.
