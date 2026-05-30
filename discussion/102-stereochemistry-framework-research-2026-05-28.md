# Stereochemistry framework research

Status: **Research / design proposal**. Created independently from
`101-stereochemistry-framework-2026-05-28.md`.

## Summary

umol should model stereochemistry as a molecule-level configuration layer, not
as atom and bond fields. The layer should contain stereochemical elements that
refer to structural carriers and ligand sites, plus a separate assertion layer
that says whether the configuration is absolute, relative, racemic, unknown, or
a constrained set of possibilities.

The most useful conceptual core is:

```
Molecule = ConstitutionGraph + StereoLayer

StereoLayer =
  elements: [StereoElement]
  assertions: StereoAssertionExpression
  provenance: parse/perception/source annotations
```

Each `StereoElement` is an equivalence class of ordered ligand-site tuples under
a declared local symmetry group. This generalizes tetrahedral centers and E/Z
bonds to square-planar, trigonal-bipyramidal, octahedral, allene, axial,
planar, spiro, atropisomeric, helical, and ring-conformational elements.

The representation should distinguish three things that current toolkits often
collapse:

1. **Geometric fact**: a local or extended arrangement of sites in 3D/topological
   stereospace.
2. **Nomenclature descriptor**: R/S, r/s, E/Z, M/P, Delta/Lambda, etc., derived
   from CIP or another ranking system.
3. **Sample assertion**: absolute single stereoisomer, relative single
   stereoisomer, racemate, mixture, unknown, undefined, or a more precise
   logical constraint over assignments.

## Precedents and useful ideas

### IUPAC/CIP

CIP remains the best available basis for atom-index-independent ordering. The
2018 Hanson/Mayfield/Yerin/Redkin/Musacchio paper is the practical starting
point for implementation: it is explicitly a guide for machine implementation
of CIP rules and is the basis of modern RDKit CIP labelling. RDKit's
`rdCIPLabeler` documents that it is a port of John Mayfield's `centres`
implementation and assigns CIP labels using an accurate CIP implementation:
<https://www.rdkit.org/docs/source/rdkit.Chem.rdCIPLabeler.html>.

`centres` is worth studying directly. It is a Java library for perception and
labelling of stereogenic centres using CIP priority rules, supports R/S/r/s and
E/Z, and includes CDK/OPSIN/JChem integration points:
<https://github.com/SiMolecule/centres>.

Important conclusion: CIP is a ranking and descriptor system, not a complete
internal stereochemical data model. Store CIP-derived orderings/descriptors as
derived annotations or canonicalized views over explicit stereo elements. Do
not make R/S or E/Z the only representation of the stereochemical state.

### InChI stereo layers

InChI is a strong precedent for separating whole-molecule stereo interpretation
from local parities. It stores tetrahedral configuration as local parities in
the `/t` layer, then uses `/m` and `/s` to express mirror choice and whether the
stereo is absolute, relative, or racemic. The InChI FAQ explicitly says that
enantiomers can have identical `/t` layers but differ by `/m0` or `/m1`, and
that `/s1`, `/s2`, `/s3` encode absolute, relative, and racemic stereo:
<https://www.inchi-trust.org/technical-faq/>.

Useful idea for umol: do not encode "relative" by weakening each center. Encode
it as a whole-molecule assertion over a vector of local states.

### V3000 / enhanced stereochemistry / RDKit StereoGroups

V3000 enhanced stereo and RDKit StereoGroups are practical precedents for
sample-level assertions. RDKit stores enhanced stereochemistry as molecule-level
`StereoGroup` objects of type ABS/AND/OR. RDKit describes AND as a mixture, OR
as an unknown single substance, and ABS as a single substance, following V3000
molfile conventions:
<https://www.rdkit.org/docs/RDKit_Book.html#support-for-enhanced-stereochemistry>.

Limit: ABS/AND/OR groups are useful but not a full logical language. Depth-First
notes that enhanced representation is flexible but cannot represent all
mixtures:
<https://depth-first.com/articles/2022/02/09/v3000-molfile-enhanced-stereochemistry-representation/>.

Useful idea for umol: support ABS/AND/OR import/export, but internally consider
a more general `StereoAssertionExpression`.

### OpenSMILES

OpenSMILES confirms that the notation already has more than tetrahedral atoms
and double bonds: `TH`, `AL`, `SP`, `TB`, and `OH` cover tetrahedral,
allene-like, square planar, trigonal bipyramidal, and octahedral classes. It
also notes that few SMILES systems implement SP/TB/OH:
<https://opensmiles.org/opensmiles.html>.

This supports prioritizing a generic shape/permutation framework instead of
hard-coding tetrahedral and cis/trans fields.

### Open Babel

Open Babel stores stereo as molecule-attached stereo data, not directly as atom
or bond properties. Its docs say stereo records use stable atom/bond Ids rather
than mutable indices and include a special implicit reference for hydrogens or
lone pairs. It stores data such as "looking from atom 2, atoms 4, 5 and 6 are
arranged clockwise around atom 3":
<https://openbabel.org/docs/Stereochemistry/stereo.html>.

Open Babel also has enum support for cis/trans, extended cis/trans, square
planar, tetrahedral, extended tetrahedral, trigonal bipyramidal, and octahedral
stereo:
<https://openbabel.org/api/3.0/structOpenBabel_1_1OBStereo.shtml>.

Useful idea: stable identity references matter. Avoid vector indices as the
semantic identity of stereo carriers.

### Chemical graph transformation with stereo-information

Andersen, Flamm, Merkle, and Stadler propose a graph transformation framework
where local geometry is represented by ordered incident-edge lists and
permutation groups define equivalence classes of orderings that correspond to
the same spatial embedding. They explicitly support partially specified
stereoinformation and suggest generalization beyond tetrahedral geometry:
<https://publica.fraunhofer.de/entities/publication/203a47e4-2594-4818-9275-8694c16ac638>.

This is directly aligned with umol's needs for transformations and reactions.

### Molassembler stereopermutators

SCINE Molassembler is one of the closest open-source precedents to the desired
model. It splits the molecule into an undirected graph plus stereopermutators.
Stereopermutators are not named stereocenters because they also represent
nonstereogenic shape/ranking cases. Atom-centered stereopermutators handle
local shapes, and bond-centered stereopermutators handle configurational
permutations due to rotational barriers:
<https://pubs.acs.org/doi/10.1021/acs.jcim.0c00503>.

The Python docs say Molassembler ranks substituents, converts rankings into an
abstract case such as octahedral `(A-A)BBCD`, and symbolically computes which
permutations are not superimposable by rotations:
<https://scine.ethz.ch/static/download/documentation/molassembler/v3.0.0/py/stereopermutators.html>.

This is the strongest precedent for using explicit local spatial symmetry and
permutation spaces.

### StereoMolGraph

StereoMolGraph is a new open-source Python library focused on graph
representations of molecules and reactions with stereochemistry. Its docs say
it uses local stereodescriptors derived from group-theoretical principles and
implements graph algorithms for stereochemical equivalence, symmetry, and
transformations:
<https://stereomolgraph.readthedocs.io/dev/>.

Its descriptor docs define stereo as an ordered atom tuple plus parity, with a
`PERMUTATION_GROUP` giving allowed permutations. It includes tetrahedral,
square planar, trigonal bipyramidal, octahedral, planar-bond, and atropisomeric
bond descriptors:
<https://stereomolgraph.readthedocs.io/dev/reference/stereodescriptors.html>.

The package is MIT licensed, beta status, and interoperates with RDKit:
<https://pypi.org/project/StereoMolGraph/>.

This is a good reference implementation for permutation-invariant descriptor
equality and graph isomorphism with stereo.

### MCDL

MCDL uses separate stereochemical modules and canonical ordering to make
descriptors independent of input drawings. It explores equivalent numbering
schemes in quasi-symmetric structures and selects canonical stereo descriptors.
It is limited mostly to atom and double-bond stereo, but the modular approach is
useful:
<https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-3-5>.

### Stereoisomer enumeration and automorphism groups

The `enu` isomer enumerator uses graph automorphism groups to enumerate
stereoisomers of a constitutional isomer and is open source as part of CombiFF:
<https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00677-6>.

For umol, automorphism groups are not optional. They are needed for:

- meso detection;
- equality/canonicalization of stereo layers;
- deciding whether an apparent center is stereogenic;
- determining enantiomeric/diastereomeric relationships;
- generating nonredundant stereoisomers;
- deciding whether transformations preserve, invert, erase, or create stereo.

General-purpose tools worth knowing: bliss computes automorphism groups and
canonical labels for colored graphs (<https://www.tcs.hut.fi/Software/bliss/>);
Rust has `nauty-pet` for petgraph graphs (<https://docs.rs/nauty-pet>) and
`graph-canon` for canonical labeling/automorphisms
(<https://docs.rs/graph-canon>).

## Proposed umol representation

### Core types

Use names like these conceptually. Exact Rust names can follow local style.

```rust
struct StereoLayer {
    elements: Vec<StereoElement>,
    assertions: StereoAssertionExpr,
}

struct StereoElement {
    id: StereoElementId,
    kind: StereoKind,
    carrier: StereoCarrier,
    sites: Vec<StereoSite>,
    shape: StereoShape,
    state: StereoState,
    ranking: Option<RankingSnapshot>,
    source: StereoSource,
}
```

`StereoElement` is molecule-level data. It refers to atoms, bonds, axes,
planes, paths, fragments, and implicit or virtual sites, but it does not live
inside those referenced structures.

### Carriers

A carrier is the structural locus whose embedding is being constrained.

Suggested carrier variants:

| Carrier | Examples |
| --- | --- |
| `Atom(atom)` | tetrahedral carbon, pyramidal nitrogen, sulfoxide sulfur, metal center |
| `Bond(bond)` | alkene E/Z, atropisomeric single bond, hindered amide |
| `AtomPath([atoms])` | allene/cumulene axis, spiro center, ring trans-cycloalkene axis |
| `BondPath([bonds])` | cumulene, conjugated chain, axis along multiple bonds |
| `Axis { near, far }` | atropisomerism, allenes, biphenyls, helicity projections |
| `Plane { frame, out_of_plane_sites }` | planar chirality, metallocenes, paracyclophanes |
| `Helix { path, sense_sites }` | helicenes, conformational helicity |
| `ConformationLock { scope }` | trans-cycloalkene, restricted ring conformations |

The key point: a carrier can be extended and need not be an atom or a bond.

### Sites

For central chirality and stereobond-like elements, a site is not merely a
neighbor atom. It is a ligand position in a local embedding.

Suggested `StereoSite` variants:

| Site | Purpose |
| --- | --- |
| `Atom(atom)` | direct neighbor when sufficient |
| `Bond(bond)` | distinguishes attachment bond, useful for reaction mapping |
| `HalfEdge { carrier, neighbor, bond }` | best default for atom-centered stereo |
| `FragmentRoot { root, through_bond }` | ligand identity begins at an attachment |
| `FragmentOrbit { root, through_bond, orbit }` | symmetric/equivalent ligand set |
| `ImplicitHydrogen(carrier)` | implicit H without inventing an atom |
| `LonePair(carrier, slot)` | pyramidal N, sulfoxide-like or nonbonding site |
| `Vacancy(carrier, slot)` | square planar/open coordination sites if needed |
| `VirtualDuplicate(node)` | CIP duplicate nodes for multiple bonds/rings |
| `PseudoSite(label)` | import/export placeholder for unknown attachment |

Recommendation: use `HalfEdge` or `FragmentRoot` as the normal site identity,
not raw neighboring atoms alone. Neighbor atoms are convenient but lose
information under graph edits, reactions, haptic/coordinate bonds, duplicate
ligands, and implicit-site handling.

The full substituent fragment should not be stored inside every stereo element.
Instead, store a root site plus an optional cached ranking/orbit snapshot. The
fragment is derived from the molecule graph when needed. This avoids duplicated
structure while still allowing CIP and symmetry to use full substituent
fragments.

### Shape and permutation group

Every stereo element should declare:

```
shape.vertices = n
shape.rotation_group = subgroup of S_n preserving handedness/configuration
shape.full_group = rotation_group plus reflections/inversions if relevant
site_order = ordered tuple of n sites
state = equivalence class of site_order under rotation_group
```

Examples:

| Shape | Vertices | State model |
| --- | ---: | --- |
| Tetrahedral | 4 | two enantiomeric classes if ABCD |
| Trigonal pyramidal | 4 including lone pair | usually tetrahedral-like, but inversion barrier policy matters |
| Square planar | 4 | usually achiral local classes, multiple ligand arrangements |
| Trigonal bipyramidal | 5 | apical/equatorial-sensitive permutations |
| Octahedral | 6 | fac/mer, Delta/Lambda, cis/trans relations |
| Planar double bond | 4 plus central bond | E/Z or cis/trans after ranking |
| Allene/cumulene | terminal sites around an axis | axial chirality class |
| Atropisomeric bond | four ortho/flanking sites around a hindered axis | axial chirality class |
| Planar chirality | plane frame plus out-of-plane site(s) | pR/pS or equivalent |
| Helical | ordered path plus sense | P/M or Delta/Lambda-like |

This mirrors Molassembler and StereoMolGraph while keeping the data independent
of current toolkit atom/bond field conventions.

### State

Represent state as a small algebra, not as strings:

```rust
enum StereoState {
    Specified(StereoClass),      // an equivalence class in the shape's state space
    Unknown,                    // explicitly unknown, e.g. wavy bond
    Undefined,                  // stereogenic but input did not specify
    NotStereogenic,             // shape present but only one arrangement possible
    Invalid,                    // input specified impossible stereo
}
```

`StereoClass` can be a canonical representative permutation plus parity, or a
compact ordinal into the shape/ranking-specific state space. For external
formats, derive R/S/E/Z/TH1/etc. as views.

### Configurational stability

Some local arrangements are geometrically chiral but not configurationally
stable under ordinary interpretation. Amines with rapidly inverting lone pairs,
conformational helicity, and rotatable atropisomeric candidates should not be
treated the same way as stable sulfoxides, quaternary ammonium centers, or
locked biaryls.

Add a policy field or derived classification:

```rust
enum StereoStability {
    Configurational,       // stable enough to identify a stereoisomer
    Conformational,        // conformer-level, not constitution-level by default
    Fluxional,             // rapidly interconverting
    Unknown,
}
```

This lets umol represent lone-pair and conformational stereo without pretending
that every geometric descriptor creates a separable stereoisomer.

### Assertion layer

Enhanced stereo groups and InChI show that local states are not enough.

Represent sample-level meaning separately:

```rust
enum StereoAssertionExpr {
    Any,                         // no sample-level assertion
    Absolute(Vec<ElementState>),  // known absolute stereoisomer
    Relative(Vec<ElementState>),  // one member of an enantiomeric pair, absolute unknown
    Racemic(Vec<ElementState>),   // equal mixture of enantiomeric pair
    Mixture(Vec<StereoAssertionExpr>),
    OneOf(Vec<StereoAssertionExpr>),
    AllOf(Vec<StereoAssertionExpr>),
    Not(Box<StereoAssertionExpr>),
}
```

For a first implementation, this can be limited to:

- `Absolute(group)`;
- `Relative(group)`;
- `Racemic(group)`;
- `Mixture(groups)`;
- `OneOf(groups)`;
- `Unspecified`.

But the data model should not make ABS/AND/OR the ceiling. A logical expression
over element states can represent cases that V3000 enhanced stereo cannot.

### Meso forms

Do not store "meso" as a local stereo tag. Meso is a property of the molecule
under the stereo layer and graph automorphism/mirror operation:

```
is_meso = has_stereogenic_elements
          && mirror(stereo_layer) is automorphic to stereo_layer
          && not all local elements are individually achiral
```

This requires graph automorphism and stereo-aware isomorphism. Store a cached
classification if useful, but derive it from the graph plus stereo layer.

## CIP and canonical ordering strategy

The user's desired property, atom-index-independent internal notation, should
be achieved in two layers:

1. `StereoElement.sites` store stable structural references.
2. `RankingSnapshot` stores a derived canonical order of those sites using CIP
   or a chosen ranking policy.

```rust
struct RankingSnapshot {
    policy: RankingPolicy,       // CIP2013, CIP2018, umol-canonical, input-order
    ranks: Vec<RankClass>,       // ties allowed
    ordered_sites: Vec<StereoSite>,
    descriptor: Option<StereoDescriptor>,
    dependencies: RankingDeps,   // atoms/bonds/isotopes/stereo elements used
}
```

Why not store only CIP labels?

- CIP descriptors can change when substituent priorities change even if the
  spatial arrangement is retained.
- Reactions often alter the ranking context without altering the local
  geometry.
- CIP has pseudoasymmetric and ligand-reflection cases where the surrounding
  stereochemistry affects ranking.
- Non-CIP shapes and partially specified stereo need representation before a
  descriptor is available.

So: store the geometric/permutation state as truth, and store CIP as a
canonical view. When a graph edit affects `RankingDeps`, invalidate and
recompute descriptors.

## Central and stereobond substituent representation

For chiral-center/stereobond-style elements, the best default is:

```
StereoSite = HalfEdge { carrier, neighbor, bond }
```

with optional fragment/ranking caches.

Comparison of alternatives:

| Representation | Pros | Cons | Recommendation |
| --- | --- | --- | --- |
| Neighbor atoms | compact, matches many toolkits, easy for SMILES | ambiguous with implicit H/lone pairs, weak under edits, not enough for haptic/duplicate sites | use only as shorthand |
| Connecting bonds | stable through atom substitution, reaction-friendly | insufficient without endpoint/site role | include in `HalfEdge` |
| Full substituent fragments | chemically complete, direct ranking | duplicates graph, hard to update, expensive | derive, do not store |
| Half-edge/rooted ligand site | compact, explicit attachment, derives full fragment | requires site abstraction | **best default** |
| Fragment orbit | handles symmetry and equivalent ligands | needs automorphism computation | use for canonicalization/perception |

For double bonds and axial elements, use role-labelled sites:

```
PlanarBond {
  carrier: Bond(c1-c2),
  left_sites: [site at c1, site at c1],
  right_sites: [site at c2, site at c2],
}

Axis {
  near_sites: [site_a, site_b],
  far_sites: [site_c, site_d],
}
```

This avoids encoding E/Z as a property of the central bond. The central bond is
only the carrier; the stereo fact is an arrangement of sites around it.

## Perception pipeline

A robust implementation can be staged:

1. **Import preservation**
   - Parse SMILES `@`, `@@`, `/`, `\`, `@TH`, `@AL`, `@SP`, `@TB`, `@OH` into
     raw `StereoElement`s using input-local ordering provenance.
   - Parse CTAB wedge/hash and V3000 enhanced stereo into elements plus
     assertions.
   - Preserve invalid/impossible stereo as diagnostics rather than silently
     dropping it.

2. **Candidate generation**
   - Find possible local elements by valence/geometry/topology:
     tetrahedral/trigonal pyramidal, planar double bonds, allenes/cumulenes,
     square planar, TB, OH.
   - Add extended candidates later: atropisomeric bonds, spiro, planar,
     helical, trans-cycloalkene.

3. **Site construction**
   - Build `HalfEdge` sites, implicit H, lone pair, vacancy, and virtual sites.
   - Attach source provenance: SMILES order, 2D wedge, 3D coordinates, inferred.

4. **Symmetry and stereogenicity**
   - Compute graph automorphism or local ligand orbits.
   - Determine whether a candidate has more than one non-superimposable state.
   - Mark nonstereogenic or invalid specified elements explicitly.

5. **State assignment**
   - From 2D: wedge/hash, crossed/wavy, bond directions.
   - From 3D: signed volumes, dihedrals, local shape fitting.
   - From 0D: input stereo tokens only; otherwise `Undefined`.

6. **Canonical ranking and descriptors**
   - Assign CIP ranks/descriptors with a CIP module.
   - Cache descriptors separately from geometric state.

7. **Stereo-aware canonicalization**
   - Canonicalize constitution graph and stereo layer together.
   - Use automorphism groups to avoid duplicate stereoisomers and detect meso.

## Algorithm buckets and implementation references

1. **CIP ranking and descriptors**
   - Study Hanson/Mayfield's algorithmic CIP analysis and the `centres` Java
     code: <https://github.com/SiMolecule/centres>.
   - RDKit's current CIP labeler is a port of `centres`:
     <https://www.rdkit.org/docs/source/rdkit.Chem.rdCIPLabeler.html>.
   - Critical implementation details: hierarchical digraphs, duplicate nodes
     for multiple bonds/rings, isotope ordering, pseudoasymmetric r/s,
     Rule 5 enantiomorphic/diastereomorphic ligand comparison, and recursive
     termination in symmetric graphs.

2. **Format stereo preservation**
   - OpenSMILES `TH`, `AL`, `SP`, `TB`, `OH`:
     <https://opensmiles.org/opensmiles.html>.
   - Open Babel molecule-level stereo data and stable Id references:
     <https://openbabel.org/docs/Stereochemistry/stereo.html>.
   - RDKit enhanced stereo groups:
     <https://www.rdkit.org/docs/RDKit_Book.html#support-for-enhanced-stereochemistry>.

3. **Shape/permutation stereochemistry**
   - Molassembler stereopermutators and abstract shape cases:
     <https://scine.ethz.ch/static/download/documentation/molassembler/v3.0.0/py/stereopermutators.html>.
   - StereoMolGraph permutation groups and local descriptors:
     <https://stereomolgraph.readthedocs.io/dev/reference/stereodescriptors.html>.

4. **Graph transformations with stereo**
   - Ordered-list plus permutation-group method from Andersen/Flamm/Merkle/Stadler:
     <https://publica.fraunhofer.de/entities/publication/203a47e4-2594-4818-9275-8694c16ac638>.

5. **Automorphism/canonicalization**
   - bliss for automorphism groups and canonical forms:
     <https://www.tcs.hut.fi/Software/bliss/>.
   - Rust options: `nauty-pet` (<https://docs.rs/nauty-pet>) and
     `graph-canon` (<https://docs.rs/graph-canon>).

6. **Stereoisomer enumeration**
   - RDKit `EnumerateStereoisomers`:
     <https://rdkit.org/docs/source/rdkit.Chem.EnumerateStereoisomers.html>.
   - `enu`/CombiFF uses automorphism groups for stereoisomer enumeration:
     <https://jcheminf.biomedcentral.com/articles/10.1186/s13321-022-00677-6>.

## Implementation roadmap

### Phase 1: data model and import preservation

Implement a molecule-level `StereoLayer` with:

- tetrahedral elements;
- planar double-bond elements;
- extended tetrahedral/allene elements;
- `StereoSite::HalfEdge`, `ImplicitHydrogen`, `LonePair`;
- `StereoState::{Specified, Unknown, Undefined, Invalid}`;
- assertion groups that can round-trip ABS/AND/OR.

Do not attempt full CIP first. Preserve input stereo and make it inspectable.

### Phase 2: CIP ranking and descriptors

Implement or port a CIP ranking engine inspired by Hanson/Mayfield and
`centres`. Start with:

- tetrahedral R/S and pseudoasymmetric r/s;
- double-bond E/Z;
- isotopes;
- duplicate nodes for multiple bonds and rings;
- Rule 5 handling for enantiomorphic/diastereomorphic ligands.

Use the CIP validation suite referenced by `centres`:
<https://cipvalidationsuite.github.io/ValidationSuite/>.

### Phase 3: symmetry and canonical stereo

Add stereo-aware automorphism/canonicalization:

- constitution graph automorphisms;
- ligand orbit detection;
- mirror/inversion operation on the stereo layer;
- enantiomer/diastereomer/identical classification;
- meso detection.

Consider Rust graph tooling (`graph-canon`, `nauty-pet`) or a local
individualization/refinement implementation if external dependencies are too
heavy.

### Phase 4: non-tetrahedral and extended carriers

Add shape/permutation descriptors for:

- square planar;
- trigonal bipyramidal;
- octahedral;
- allene/cumulene;
- atropisomeric axis;
- spiro;
- planar chirality;
- helicity;
- trans-cycloalkene and other conformationally locked elements.

Molassembler and StereoMolGraph are the best guides here.

### Phase 5: reactions and transformations

Stereo changes in reactions should be operations over elements:

- retain element if carrier/sites map and orientation is preserved;
- invert if the local embedding operation changes parity;
- erase if carrier becomes planar/free-rotating or mapping is ambiguous;
- create `Undefined` if a new stereogenic element is formed without stereo
  control;
- create `Relative`/`Racemic` assertions when the reaction specification says
  so.

This is much easier if stereo is a top-level relation over structural elements
instead of a property hidden on atoms and bonds.

## Open-source code to study

| Project | Language | What to inspect |
| --- | --- | --- |
| RDKit | C++/Python | `CIPLabeler`, stereo perception, `StereoGroup`, stereo enumeration |
| Centres | Java | CIP implementation, validation suite, Rule 5 handling |
| CDK | Java | stereo element classes, CIP APIs, non-tetrahedral classes |
| Open Babel | C++ | molecule-level `StereoData`, stable refs, perception flow |
| Molassembler | C++/Python | stereopermutators, shape classification, ranking propagation |
| StereoMolGraph | Python | permutation groups, stereo-aware isomorphism, reaction graphs |
| Indigo | C++ | enhanced stereo and production cheminformatics behavior |
| OPSIN | Java | parsing IUPAC stereochemical descriptors into structures |
| enu / CombiFF | C++ | stereoisomer enumeration via automorphism groups |

## Design decisions recommended for umol

1. Store stereo only at molecule scope.
2. Represent stereo as shape-specific ordered site tuples modulo explicit
   permutation groups.
3. Use `HalfEdge`/rooted ligand sites as the default substituent representation.
4. Store full substituent fragments only as derived ranking context, never as
   duplicated stereo payload.
5. Separate geometric state, CIP descriptor, and sample assertion.
6. Treat meso as derived from mirror plus automorphism, not as a stored local
   descriptor.
7. Make unknown and undefined distinct.
8. Preserve invalid input stereo with diagnostics.
9. Make ABS/AND/OR import/export a compatibility layer over a more general
   assertion expression.
10. Prioritize tetrahedral, double bond, allene/extended tetrahedral, and
    non-tetrahedral SMILES classes before planar/helical/conformational cases.

## Main risk

The hardest part is not the data structure; it is canonicalization under
symmetry while rankings may depend on other stereo elements. That argues for an
incremental architecture:

- first preserve and expose stereo elements;
- then add CIP as a derived pass;
- then add automorphism-aware canonicalization;
- then broaden shapes.

This keeps the framework flexible without requiring all stereochemistry to be
solved before umol can round-trip and inspect useful stereo data.
