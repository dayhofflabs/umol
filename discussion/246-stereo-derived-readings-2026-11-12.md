# 246 — Derived stereo readings

Status: In Progress
Date: 2026-11-12
Relates: [245](245-canonical-eq-cost-2026-11-12.md)

## Purpose

A stereo atom or stereo bond has three derived readings: its stereogenicity, the topicity of a
pair of its ligands, and whether a permutation of its ligands is a symmetry. The molecule has two:
whether it is chiral, and which atoms are equivalent. The stereo transformers select sites by a
further property, which this document separates from stereogenicity: whether a site is a
stereogenic unit of the compound.

Issue 43 reports a CH₂ carbon beside an assigned stereocenter that is classified as `Stereogenic`.
The review of that report found that the present derivation gives wrong results for molecules with
two or more assigned stereo entities, in all five readings.

This document records

- the classes of `Stereogenicity` and `Topicity` and their definitions,
- where each reading is asked,
- the one definition of the colors that the readings and canonicalization share, and
- how the readings are computed.

S4f, S4g and S5a-S5c are implemented. The reviewed S4/S5 semantics, interfaces and names are
settled. Main has the selected Constitution search, flat auxiliary representation, narrower
generator filter, FieldColor names and the three derived stereo-view readings. S5d's stereogenic
units are next; the workspace-wide S5 gate remains ahead.

## Starting state

Before S5b, the readings were derived in `umol-graph-ir/src/ir/molecule/symmetry.rs` through two
public helper values. S1 moved that implementation without changing these computations. This
section describes the implementation used for the original baselines.

1. `Molecule::graph_symmetry(config)` returns `GraphSymmetry`. It runs the automorphism search on
   the incidence graph. The color of a stereo node includes its stored coset index, reduced by
   exchanges of ligands in one orbit, and with `iterate_to_fixpoint` the colors are recomputed
   until the orbits are stable. Each generator is then graded by its action on the stereo
   entities with a literal coset: proper if it preserves every configuration, improper if it
   carries every configuration to its mirror image, and discarded otherwise. `GraphSymmetry`
   answers `is_chiral`, `same_proper_orbit`, `same_star_orbit`, `proper_orbit_of` and
   `star_orbit_of`.
2. `Molecule::stereo_atom_symmetry(&graph_symmetry, id)` and `stereo_bond_symmetry` return
   `StereoSymmetry`, the local group of one stereo entity. It is built from the generators of the
   stabilizer of the site, graded one at a time and projected onto the ligand positions. For a
   tetrahedral stereo atom the parity of a generator's permutation of its own ligands counts as
   its orientation there. `StereoSymmetry` answers `stereogenicity`, `is_stereogenic`, `topicity`
   and `group`.

`is_stereogenic` is true when the coset is alone in its orbit under the permutations of the local
group. `stereogenicity` returns `Stereogenic` in that case, `Prostereogenic` when some pair of
ligand positions is enantiotopic, and `Nonstereogenic` otherwise. `topicity` returns `Homotopic`
for two positions in one proper orbit of the local group, `Enantiotopic` for two positions that
only an improper element relates, and `Diastereotopic` for every other pair, including two
ligands of different constitution.

The readings are consumed as follows.

- `StereoPerceiver::classify`, the conformance validator and the stereo enumerator build a
  `GraphSymmetry`, derive a `StereoSymmetry` per stereo entity, and call `stereogenicity`,
  `topicity` and `group().contains(..)` on it. The validator and the enumerator compare the result
  with the assertion themselves.
- The stereo views have methods that take a `&StereoSymmetry` and forward to it: `stereogenicity`,
  `is_stereogenic`, `is_prostereogenic`, `topicity`, `is_homotopic`, `is_enantiotopic`,
  `is_diastereotopic` and `ligand_symmetry`. They do not read the view, and nothing ties the
  passed value to the view's entity.
- The keyed reading `constraints().derived(key)` of a stereo atom or stereo bond returns the
  undetermined form for every key. `evidenced`, `satisfies` and `is_compatible` therefore read the
  assertion alone.
- `StereoModel::para_stereo` sets `iterate_to_fixpoint`. The stereo operation configs carry the
  limit `para_stereo_max_iterations`.
- No operation calls `is_chiral` or the orbit queries.

### Measurements

A temporary test built three compounds with explicit hydrogen atoms on the stereocenters and
literal cosets, and read the present functions with the fixpoint on.

| Molecule and reading | Expected | Returned |
| --- | --- | --- |
| 1,2-dichloro-1,2-difluoroethane, meso form: `is_chiral` | false | true under one atom numbering, false under another |
| the same compound, chiral form: the two carbons in one proper orbit | yes | yes under one numbering, no under the other |
| 2,3,4-trichloropentane, both meso forms: class of C3 | `Stereogenic` | `Prochiral` |
| 2,3,4-trichloropentane, chiral form: class of C3 | not stereogenic | `Stereogenic` |
| 1,3-dichlorocyclobutane, cis and trans: class of C1 and C3 | `Stereogenic` | `Prochiral` |
| 1,3-dichlorocyclobutane, cis and trans: C1 and C3 in one proper orbit | yes | yes in one isomer, no in the other |

The two numberings differ in the ids of H and F on the second carbon. At C3 of
2,3,4-trichloropentane the topicity of the two CHCl branches was returned as enantiotopic in the
meso forms and diastereotopic in the chiral form, which is correct.

The results have three causes.

- **Stored coset index in the color.** The index is relative to the stored ligand order. Under
  some numberings two equivalent stereo entities get different colors, and the automorphism that
  relates them is not found.
- **Generators graded one at a time.** The proper orbits are built from the proper generators
  alone, so a proper operation that is the product of two improper generators is lost.
- **Stereogenicity read from symmetry operations.** A remapping that exchanges two equal ligands
  and leaves every other configuration in place shows that the two configurations of the center
  are the same molecule. It inverts the center's own configuration, so it is not a symmetry
  operation of the molecule, and it is discarded. This is the case of issue 43 and of C3 in the
  chiral form above. In the meso forms and in the ring the exchange of the two equal ligands
  inverts the other centers as well. It is kept as an improper operation, and the center is
  reported as not stereogenic.

The original reading tests used molecules with a single stereo entity. S3 adds remapping-check
tests with multiple stereo entities; the new reading tests are part of the remaining plan.

## Sources

- K. Mislow and J. Siegel, "Stereoisomerism and Local Chirality", J. Am. Chem. Soc. 1984, 106,
  3319–3328. In `materials/stereochemistry`.
- IUPAC, Rules for the Nomenclature of Organic Chemistry, Section E: Stereochemistry
  (Recommendations 1974), Pure Appl. Chem. 1976, 45, 11–30. In `materials/stereochemistry`.
- IUPAC, Basic Terminology of Stereochemistry (Recommendations 1996), Pure Appl. Chem. 1996, 68,
  2193–2222.
- M. Razinger, K. Balasubramanian, M. Perdih, M. E. Munk, "Stereoisomer Generation in
  Computer-Enhanced Structure Elucidation", J. Chem. Inf. Comput. Sci. 1993, 33, 812–825. In
  `materials/stereochemistry`.
- M. Perdih, M. Razinger, "Stereochemistry and Sequence Rules: A Proposal for Modification of
  Cahn-Ingold-Prelog System", Tetrahedron: Asymmetry 1994, 5, 835–861. In
  `materials/stereochemistry`.

Mislow and Siegel separate two properties of a center. Stereogenicity is a property of
permutations: a center is stereogenic when the exchange of two ligands gives a stereoisomer.
Topic relations are statements about symmetry: two ligands are homotopic when a proper operation
relates them, enantiotopic when only an improper one does, diastereotopic when they have the same
constitution and no operation relates them, and constitutionally heterotopic when their
constitution differs (their Figure 1). Their compounds 13 and 14 have a center that is not
stereogenic and whose two equal ligands are diastereotopic, so one property cannot be read from
the other.

| Center | Class | Equal ligands | Source |
| --- | --- | --- | --- |
| C of CHBrClF | `Stereogenic` | none | Mislow–Siegel; Rule E-4.3 |
| C-3 of meso-2,3,4-trihydroxyglutaric acid (pseudoasymmetric) | `Stereogenic` | none; two mirror-image groups | Mislow–Siegel compounds 11, 12; Rule E-4.8 |
| C-1 of ethanol | `Prostereogenic` | enantiotopic H | Rule E-4.12, example A |
| C-3 of glyceraldehyde, C-3 of butan-2-ol | `Prostereogenic` | diastereotopic H | Rule E-4.12, example B; IUPAC 1996 |
| C-3 of chiral 2,3,4-trihydroxyglutaric acid | `Prostereogenic` | diastereotopic CHOHCOOH groups | Mislow–Siegel compounds 13, 14 |
| methyl C of acetic acid | `Nonstereogenic` | homotopic H | Mislow–Siegel p. 3326 |

Razinger et al. generate the stereoisomers of a constitution. Every potential stereocenter is
classified as a non-stereocenter, a true-stereocenter or a para-stereocenter. A true-stereocenter
has constitutionally different substituents. A para-stereocenter is "a structural feature which,
although not a true-stereocenter, does influence the determination of the number and kinds of
stereoisomers"; three rules that count stereocenters inside the identical substituents decide it,
the first applied recursively to rings. Every stereocenter, true or para, gets a parity, and the
stereoisomers are the classes of parity vectors under the automorphisms of the graph, each of
which flips the parity of a para-stereocenter whose equivalent neighbors it exchanges.

Perdih and Razinger define the same by automorphisms. A potential stereogenic center "really is a
stereogenic centre if there exist no automorphic permutation causing an inversion of the
configuration of only the potential stereogenic centre under consideration". If one exists, it
"can be a stereogenic centre if the number of topologically equivalent neighbours (ligands) ... is
less than or equal to the number of different configurations of these ligands", and it is an
actual one only when those ligands differ in configuration. They name the two kinds
constitutional and configurational stereogenic centers.

## Design

### Classes

`Stereogenicity` is a property of one stereo atom or stereo bond.

| Variant | Before | Meaning |
| --- | --- | --- |
| `Stereogenic` | `Stereogenic` | Another coset of its kind gives a different molecule. |
| `Prostereogenic` | `Prochiral` | Not stereogenic; replacing one ligand by a new one makes it stereogenic. |
| `Nonstereogenic` | `Symmetric` | Not stereogenic, and no single replacement makes it stereogenic. |

`Stereogenic` follows Mislow and Siegel (footnote 24): the center can have two non-equivalent
configurations. The definition is the same for every kind.

| Kind | Cosets |
| --- | ---: |
| tetrahedral, cis/trans, axial | 2 |
| square planar | 3 |
| trigonal bipyramidal | 20 |
| octahedral | 30 |

With two cosets there is one other coset, and the definition reads "the other configuration is a
different molecule". With more cosets the present rule, "the coset is alone in its orbit", is a
different statement and is wrong. Square-planar Ma₂b₂ has three labeled arrangements, one trans
and two cis that are the same molecule. The cis coset is not alone in its orbit, and the center is
stereogenic because trans is a different molecule. Trigonal-bipyramidal Pa₄b, the example of
Mislow and Siegel, has twenty arrangements in two classes, and none is alone in its orbit. The
class of a stereo entity does not depend on its own coset.

`Topicity` is a relation between two ligands of one stereo atom or stereo bond. The symmetry
operations are those that keep the stereo entity in place.

| Variant | Before | Meaning |
| --- | --- | --- |
| `Homotopic` | `Homotopic` | The same ligands, related by a proper symmetry operation. |
| `Enantiotopic` | `Enantiotopic` | The same ligands, related only by an improper symmetry operation. |
| `Diastereotopic` | `Diastereotopic` | The same ligands, not related by a symmetry operation. |
| `Distinct` | part of `Diastereotopic` | The ligands themselves differ. |

The names follow the sources where the sources have a usable word.

- `Prostereogenic` is the word of Mislow and Siegel (footnote 74) for the class that IUPAC calls
  prochiral. `Prochiral` describes a permutation class by a chirality word. The two terms cover
  the same centers except a CH₂ between two mirror-image groups: Rule E-4.12 excludes it, and
  the replacement definition includes it.
- `Nonstereogenic` names the remainder. The sources count replacement steps and have a word only
  for the two-step case, "proprochiral" (IUPAC 1996). `Symmetric` describes symmetry, which does
  not select this class: a stereogenic center can lie on a twofold axis, and a center without
  any symmetry can be not stereogenic.
- `Distinct` replaces the source term "constitutionally heterotopic". In IUPAC 1996 the bare word
  "heterotopic" refers to "stereoheterotopic", which means enantiotopic or diastereotopic.

### Colors

The search starts from the colors of the entities and of the incidences. The code has two
definitions of the color of an entity of a molecule.

- `ir::coloring` has the trait `MoleculeColoring`, its implementation `ConstitutionColoring`, and
  the switches `MoleculeColoringFeatures`. The color is a hash of the selected inherent fields.
  Its only user is the present `Molecule::graph_symmetry`.
- `AsColor` with `initial_colors` gives the colors of canonicalization and of `canonical_eq`. The
  color of a form is its inherent fields and its entity constraints, as far as they do not
  depend on a ligand frame. Each incidence has a color of its own.

| | `ConstitutionColoring` | `AsColor` with `initial_colors` |
| --- | --- | --- |
| inherent fields | selected by 18 switches | all |
| entity constraints | never part of a color | part of the color |
| entity kinds | all eight | chosen by the description level |
| incidences | no color | the role and, for an atom of an aromatic system or a multicenter bond, its electron count |
| value | a hash | integers ranked from an ordered comparison |
| present caller's input | the molecule as stored | canonicalization supplies normalized input |

One concept has one definition. The readings and canonicalization both read `AsColor`, and the
three types of `ir::coloring` are removed. No reading gets its meaning by changing the molecule
it is given. An earlier state of this design cleared the constraints on a copy of the molecule,
to reach the meaning of `ConstitutionColoring` through `AsColor`; that is rejected.

The shared color calculation reads the supplied forms. Normalization is the responsibility of
the operation that requires it. Canonicalization already normalizes its inputs before calling
`initial_colors`; graph-symmetry operations read the molecule as supplied. A literal and a
singleton set can therefore have different colors in an unnormalized input. Sharing `AsColor`
does not require the two callers to use the same equality relation.

**Incidences.** The present orbits are too coarse where an incidence carries information, and
they change with the one definition. A temporary test read both on two molecules:

| Molecule | Atoms 0 and 1 in one orbit, present `GraphSymmetry` | The same with `initial_colors` |
| --- | --- | --- |
| a dative bond between two boron atoms, donor 0 and acceptor 1 | yes | no |
| a five-carbon aromatic ring with the electron counts 2, 1, 1, 1, 1 | yes | no |

The exchange of donor and acceptor is no symmetry of the molecule as described. Multicenter
bonds were not run; the mechanism is the same.

**Entity constraints.** The two definitions part in meaning at one point: whether an entity
constraint makes an entity different from an otherwise equal entity without it. The present
readings say no, and canonicalization says yes. Two examples, worked from the rules and not run:

| Molecule | Constraints not compared | Constraints compared |
| --- | --- | --- |
| CH₂ClF with explicit hydrogens and an entity constraint on one H | the two H are equivalent; the center is `Prostereogenic` and the pair `Enantiotopic` | the two H differ; the center is `Stereogenic` and the pair `Distinct` |
| a meso compound with a stereogenicity assertion on one of its two centers | achiral | the two centers differ, no remapping carries the configurations to their mirror images, and the molecule is chiral |

For not comparing them:

- An assertion does not change the molecule. If it is true of one atom it is true of every
  equivalent atom, stated or not, so annotating part of a molecule leaves the readings as they
  are.
- The readings are what the validator, `evidenced` and `satisfies` compare assertions with. If
  constraints are compared, an assertion can change the value it is compared with, as in the
  second example.
- It is what the readings do today.

For comparing them:

- Constraints participate in both symmetry and canonicalization. On normalized inputs the
  colors agree; graph-symmetry readings still compare the supplied forms after transport,
  while canonicalization retains its normalized comparisons.
- In a pattern a constraint restricts what an atom matches. Two pattern atoms with different
  constraints are different, and orbits that ignore this must not be used to identify pattern
  atoms or matches.
- In a description that is not concrete a constraint can be the only statement of a difference,
  for example a degree or a valence where the hydrogen count is undetermined.

Both answers are justified, so neither is prescribed. The caller chooses.

**The option.** `include_constraints` is a `bool` field of `GraphSymmetryConfig`. It says whether
the constraints of two entities are compared along with their inherent fields when the entities
are compared. Its default is `false`, which is what the readings do today. It is separate from
the description level: the level says which entity kinds take part, and the option says whether
the constraints of those entities are compared. It acts in two places.

- In the colors it decides whether the constraints are part of the color.
- In the stereo checks it decides what is compared: the configurations only, or the whole
  stereo forms with their frame-relative constraints, such as the topicity of a pair.

The readings never use the `Full` level, so molecule-level constraints are compared under
neither value. Canonicalization always compares constraints, as the guides require of every
equality, and what it computes does not change.

**`AsColor`.** The trait remains AsColor. Its full and field readings have separate associated
color types. The full reading retains its original representation and comparisons. The field
reading contains the frame-independent inherent fields without constraints. The names below
are the approved names; S4g applies them to the implementation-time InherentColor spellings.

```rust
pub(crate) trait AsColor {
    type Color<'a>: Ord + Hash
    where
        Self: 'a;

    type FieldColor<'a>: Ord + Hash
    where
        Self: 'a;

    /// The inherent fields and the constraints.
    fn as_color(&self) -> Self::Color<'_>;

    /// The inherent fields without constraints.
    fn as_field_color(&self) -> Self::FieldColor<'_>;
}
```

- The full colors keep their existing field order and comparison implementations. The atom,
  bond, dative-bond and noncovalent-bond colors remain single-reference wrappers. A color carries
  no runtime constraints flag.
- The field colors contain the selected inherent fields. Their names are AtomFieldColor,
  BondFieldColor, DativeBondFieldColor, AromaticSystemFieldColor, MulticenterBondFieldColor,
  NoncovalentBondFieldColor, StereoAtomFieldColor and StereoBondFieldColor. For stereo atoms
  and stereo bonds that is the kind; the stereogenicity assertion belongs to the full color. Each span reading
  composes its selected color type and collapses a modification whose two projected sides agree.
- `molecule::symmetry::initial_colors` and `reaction_span::symmetry::initial_colors` take
  `include_constraints` and select `as_color` or `as_field_color` directly at each entity-kind
  ranking site. Both readings use the existing generic ranker, which sorts the selected values
  and groups them with equality. The repeated branches are deliberate; no selector trait,
  marker types or additional ranking helper are needed. S4f records this selected replacement
  for the implemented selector machinery. The graph-symmetry checks retain their constraints
  option.
- The molecule and reaction-span canonicalization callers pass `true`. Their direct `as_color`
  reads, full-color order, upfront normalization and normalized stereo comparisons are unchanged.

**Rejected: constraints never part of a color.** A color would be the inherent fields for every
user, and the check would compare the constraints where they count. Canonicalization would then
start from coarser cells for every molecule with entity constraints and search more, and its
keys would have to carry the constraints separately.

**Field switches.** The switches of `MoleculeColoringFeatures` have no replacement. No operation
selects fewer fields than all of them. The orbit test of `tests/canonicalization.rs` reads
orbits without isotopes, and one unit test selects no field. A selection of fields could be
added to the one definition later, by dividing the inherent part of a color.

### `GraphSymmetry`

`GraphSymmetry` is auxiliary data computed once for a molecule. All readings use one
Constitution incidence graph, its colored search adapter, refined partition and node-to-cell
index. The value also holds its config and the stereogenic-unit cache described below. It holds
no molecule, borrowed or copied. A reading runs a search when asked, unless an exit by cell
decides it. Stereo records are read from the supplied molecule when checking a remapping.

`GraphSymmetryView` follows `RingViews`. `Molecule::graph_symmetry` borrows the molecule and
returns a view owning the computed data. `from_graph_symmetry` borrows existing data and the
molecule without recomputing or copying either. `into_graph_symmetry` consumes the view: it
moves owned data or clones borrowed data. The owning view constructor is crate-private.

The five molecule-wide queries have the molecule as their receiver and take the auxiliary
data. The view exposes the same query names and forwards to those molecule operations. Every
search check reads that borrowed molecule. The stereo views likewise already borrow their
molecule; attaching the auxiliary data gives their readings both inputs.

```rust
impl Molecule {
    pub fn graph_symmetry(&self, config: &GraphSymmetryConfig) -> GraphSymmetryView<'_>;

    pub fn is_chiral(&self, symmetry: &GraphSymmetry) -> bool;
    pub fn same_proper_orbit(&self, symmetry: &GraphSymmetry, a: AtomId, b: AtomId) -> bool;
    pub fn same_star_orbit(&self, symmetry: &GraphSymmetry, a: AtomId, b: AtomId) -> bool;
    pub fn proper_orbit_of(&self, symmetry: &GraphSymmetry, a: AtomId) -> Vec<AtomId>;
    pub fn star_orbit_of(&self, symmetry: &GraphSymmetry, a: AtomId) -> Vec<AtomId>;
}

impl<'a> GraphSymmetryView<'a> {
    pub fn from_graph_symmetry(symmetry: &'a GraphSymmetry, molecule: &'a Molecule) -> Self;
    pub fn into_graph_symmetry(self) -> GraphSymmetry;

    pub fn is_chiral(&self) -> bool;
    pub fn same_proper_orbit(&self, a: AtomId, b: AtomId) -> bool;
    pub fn same_star_orbit(&self, a: AtomId, b: AtomId) -> bool;
    pub fn proper_orbit_of(&self, a: AtomId) -> Vec<AtomId>;
    pub fn star_orbit_of(&self, a: AtomId) -> Vec<AtomId>;
}

pub struct GraphSymmetryConfig {
    pub automorphism_algorithm: AutomorphismAlgorithm,
    pub max_assignments: u32,
    pub include_constraints: bool,
}
```

- **Config.** `GraphSymmetryConfig` holds the automorphism algorithm, the limit
  `max_assignments` of the stereogenic units, and `include_constraints`. It has a `Default`: the
  algorithm the stereo operation configs default to, 4,096 assignments, and `false`. The caller
  selects all three; they are not part of a chemistry model. The coloring, the fixpoint switch
  and the iteration limit of the present config are removed.
- **Constraints.** `include_constraints` decides whether the constraints of the entities are
  compared, as described under Colors. With `false` a derived reading does not depend on an
  assertion it is compared with. Molecule-level constraints are never compared.
- **Input.** The builder and the readings do not normalize the molecule or its entity forms.
  An excluded constraint is not normalized, transported or validated by the symmetry checks.
  An included constraint participates in its stored representation after the required frame
  transport; a contradictory assertion is not by itself a normalization failure here.
- **Malformed molecules.** `Molecule::graph_symmetry` remains infallible. Preparation checks
  literal electron counts against aromatic or multicenter participants and checks the selected
  stereo attributes for kind, frame and coset-action compatibility. The stereo checks still
  run even though the Constitution graph omits stereo nodes. Compatibility checks of entity
  constraints apply only when `include_constraints` is `true`. For an input that cannot be
  encoded, every entity in the search graph gets a number of its own in place of a color. No two
  atoms then share a cell, and the exits decide the readings: every atom is alone in its orbit, each
  stereo record with a kind is `Stereogenic` and `Constitutional`, each pair of ligand atoms is
  `Distinct`, and `is_chiral` is decided by the identity alone. These are the readings of a
  molecule without symmetry. They are not correct for the malformed molecule, and the symmetry
  code does not report the malformation. This fallback is for encoding and transport
  incompatibility. It does not add canonicalization's normalization-failure rule to
  graph-symmetry operations.
- **Validity.** The value describes the molecule as it was when the value was built. After a
  change to the molecule a new value is built. Nothing ties a value to a molecule, as for
  `RingViews::from_ring_set`. Pairing data with a molecule does not validate that relationship;
  mismatched context has no specified semantic result and must not panic.
- **Caches.** Every query shares the Constitution search prepared by the factory. The first
  constitutional topicity check uses that same data. The stereogenic units are the one answer
  kept, because they are computed for all records together in S5d. Other query answers are not
  cached.

**Selected layout: one Constitution search.** The search graph omits stereo entities. At a
leaf, the image of a stereo record is found from its mapped atom or bond site; its kind,
configuration and selected constraints are then checked against that image record. A proper,
mirror or coset-change query adds its stereo predicate to incidence preservation. A Constitution
query checks incidence preservation without the stereo predicate. All these queries use the
same prepared graph and the existing search algorithm.

The trial retained SymmetrySearch in both candidates to isolate the graph-layout comparison.
The subsequent decision removes that private bundle: incidence, adapter, partition and
cell_indices belong directly to GraphSymmetry, followed by config. Preparation stays in
Molecule::graph_symmetry and the existing private cell method moves to GraphSymmetry.
GraphSymmetryConfig remains the last type in the module. Inlining the bundle adds no semantic
distinction or additional heap allocation. The rejected second search and its lazy cache have
no place in this representation.

The cell_indices field stores one integer per adapter node for constant-time cell comparisons.
It is the reverse lookup of OrderedPartition::cells, which stores the nodes of each cell.
OrderedPartition does not already store that reverse lookup. Retain the index, at four bytes
per adapter node plus vector storage, and construct it with the existing
OrderedPartition::cell_indices method, made pub(crate), instead of repeating its calculation.

| Layout | Preparation and checking | Disposition |
| --- | --- | --- |
| One Constitution search | Prepare one graph without stereo entities. Check stereo through mapped sites at the leaves. Reuse that graph for every reading. | Selected. Preparation was 35–48% faster and fixed build-and-read 37–60% faster on the four measured fixtures. |
| Structure search plus lazy Constitution search | Prepare the graph including stereo entities for stereo queries, and build a second graph for the constitutional fallback. | Rejected for this implementation. It examined the same number of final candidates on these fixtures while preparing larger backend graphs and retaining more storage. |

The Structure graph can offer extra pruning when stereo kinds or selected constraints
distinguish otherwise equivalent sites. That possible benefit was not observed in this focused
comparison. Stereo nodes are unnecessary for the query semantics, and no public layout selector
or persistent second implementation is added.

The separate-worktree comparison in S5b preserves the exact workloads, all medians, paired
changes, cell and candidate-check counts, allocation measurements and source identities. Both
candidates passed the same 254 molecule-symmetry cases, and the eight fixed fixture/config
outputs matched exactly. These results describe four fixtures rather than performance across
all molecules. The selected layout and its correctness fixes remain to be integrated into
main; selection does not complete S5b or close the other review findings.

Molecule chirality and atom equivalence have no assertion in the constraint language. They are
derived readings on the molecule and its symmetry view.

### Readings

Three constraint keys of a stereo entity are derivable: stereogenicity, the topicity of a pair,
and the ligand symmetry of a permutation. Each has an asserted, a derived and an evidenced
reading. The receiver states which: the entity view gives the derived value, and its constraints
view gives the assertion and the keyed readings.

| Key | Asserted | Derived | Evidenced |
| --- | --- | --- | --- |
| stereogenicity | `constraints().stereogenicity()`, a `StereogenicityForm` | `stereogenicity()`, a `Stereogenicity` | `constraints().evidenced(key)` |
| topicity of a pair | `constraints().topicity(pair)`, a `TopicityRelationForm` | `topicity(a, b)`, a `Topicity` | `constraints().evidenced(key)` |
| ligand symmetry of a permutation | `constraints().ligand_symmetry(permutation)`, a `BooleanForm` | `ligand_symmetry(permutation)`, a `bool` | `constraints().evidenced(key)` |

The derived and the evidenced readings need the `GraphSymmetry`. It is attached to the stereo view
once, the way a ring set is attached to the atom constraints view with `with_rings`, and the
readings take no symmetry argument:

```rust
let symmetry = mol.graph_symmetry(&config).into_graph_symmetry();
let stereo = mol.stereo_atom(id).with_symmetry(&symmetry);

stereo.stereogenicity();
stereo.topicity(a, b);
stereo.ligand_symmetry(permutation);
stereo.stereogenic_unit();
stereo.constraints().evidenced(key);
```

On `StereoAtomView` and `StereoBondView`:

```rust
pub fn with_symmetry(self, symmetry: &'a GraphSymmetry) -> Self;
pub fn stereogenicity(&self) -> Stereogenicity;
pub fn topicity(&self, a: StereoLigandPosition, b: StereoLigandPosition) -> Topicity;
pub fn ligand_symmetry(&self, permutation: OrientedLigandPermutation) -> bool;
pub fn stereogenic_unit(&self) -> StereogenicUnit;
```

- A derived reading on a view without the value attached is a caller error and panics, as a ring
  key does without a ring set. The readings also panic for an entity without a kind, as
  `stereo_atom_symmetry` does today.
- The constraints view reached from a stereo view carries the attached value.
  `StereoAtomConstraintsView` and `StereoBondConstraintsView` also have `with_symmetry`.
- `derived(key)` returns the derived value as a literal form for the stereogenicity key, a
  topicity key and a ligand-symmetry key. For the fluxionality key and for an entity without a
  kind it stays vacuous.
- `evidenced`, `satisfies` and `is_compatible` include the derivation through the keyed reading.

The validator and the enumerator read the derived values from the stereo views and keep their own
comparison with the assertion, which reports the derived value in a contradiction.

### Stereogenic units

`StereogenicUnit` is a property of a stereo record on the constitution, decided among the stereo
records the molecule has. It does not depend on the configurations stored in them. The adder, the
remover and the enumerator select by it. The adder supplies every candidate before it asks. A
single record on one center of 1,3-dichlorocyclobutane, with none on the other, is not a unit.

| Variant | Meaning | Published as |
| --- | --- | --- |
| `Constitutional` | A stereogenic unit whose ligands are constitutionally different. | "true-stereocenter" (Razinger et al.); "constitutional stereogenic centre" (Perdih and Razinger) |
| `Configurational` | A stereogenic unit with constitutionally equal ligands: its configuration changes the stereoisomer for some assignment of the other units. | "para-stereocenter"; "configurational stereogenic centre" |
| `Prostereogenic` | Not a stereogenic unit, and one after a single ligand is replaced by a new one. | "prostereogenic" (Mislow and Siegel) |
| `Nonstereogenic` | None of these. | the non-stereocenters of Razinger et al. without the prostereogenic ones |

The two kinds are told apart by the ligands, as in the papers. For the two-coset kinds a
`Constitutional` unit is one that does not depend on other units. A square-planar Ma₂b₂ center
has equal ligands and is `Configurational`, although it depends on no other unit.

"Stereogenic unit" is the term of IUPAC 1996 and of the CIP papers for a focus of
stereoisomerism. The rustdoc of `StereogenicUnit`, `Stereogenicity` and `Topicity` and of their
variants cites the papers their definitions are taken from.

`Stereogenicity` and `StereogenicUnit` are two types because they differ on the same site.

| Site | `Stereogenicity` | `StereogenicUnit` |
| --- | --- | --- |
| C1 of 1,3-dichlorocyclobutane, the other center undetermined | `Prostereogenic` | `Configurational` |
| one of the two like centers of cis,trans-1,3,5-trimethylcyclohexane | `Prostereogenic` | `Configurational` |
| C3 of 2,3,4-trichloropentane, chiral form | `Prostereogenic` | `Configurational` |
| the CH₂ carbon of 1-chloro-1-fluoropropane with explicit hydrogens | `Prostereogenic` | `Prostereogenic` |

In cis,trans-1,3,5-trimethylcyclohexane inverting one of the two like centers gives the same
molecule, so the center is not stereogenic, and all three records are needed to tell the isomer
from the all-cis one.

The reading is `stereogenic_unit()` on the stereo views. It has no asserted side. Whether one site
is a unit depends on which of the others are, so the value of every stereo record of the molecule
is computed together, on first use, and kept in the `GraphSymmetry`.

The procedure decides the two conditions of Perdih and Razinger with the search.

1. A record whose ligands lie in pairwise different refined cells is `Constitutional`.
2. A record is dropped when every other coset of its kind is reached by remappings that change
   its coset in place and leave every other kept record in place with its coset unchanged. For a
   kind with two cosets one such remapping suffices. This is repeated until no record is dropped.
3. A kept record whose coset some remapping changes alone, moving other kept records without
   changing their cosets in place, is confirmed: it stays if, for some assignment of the kept
   records that such remappings can move, its `Stereogenicity` is `Stereogenic`. These are the
   kept records that do not sit alone in their cell when the record under test is held fixed.
   Every other record is taken with the first coset of its kind. If no assignment exists the
   record is dropped, and step 2 is repeated. The kept records that are not `Constitutional` are
   `Configurational`.
4. A dropped record is `Prostereogenic` if it is kept when one of its explicit ligands is made
   distinct from every other atom, and `Nonstereogenic` otherwise.

Steps 2 and 3 read the frame action of a remapping on a record and apply it to the cosets of
the kind. They do not read the stored cosets, which is why the units do not depend on them.

Step 3 tries assignments one by one and is exponential in the number of moved records when none of
them works. It tries at most `max_assignments` assignments for one record. When the limit is
reached the confirmation stops and the record stays a unit, as the iteration limit of the present
fixpoint ends its work with the result so far. `stereogenic_unit()` does not fail. A record can
then be kept whose configuration never changes the stereoisomer.

### Transformers

`StereoSelection` selects by `StereogenicUnit`: the field `stereogenic_units` replaces
`stereogenicities`. The default of the adder is `Constitutional` and `Configurational`.
`StereoPerceiver::classify` returns the `StereogenicUnit` of each requested site.

The adder adds every recognized site to a copy of the molecule, builds its `GraphSymmetry`, and
keeps the records whose unit is selected.

- The sequence of issue 43, add, enumerate, add, no longer adds records at the CH₂ carbons. They
  are `Prostereogenic`, and the units do not depend on the configuration that the enumeration
  assigned.
- Ring cis/trans centers and pseudoasymmetric centers are `Configurational` and are added by the
  default adder from the start. One add and one enumerate reach every stereoisomer.
- A caller who wants only the units with constitutionally different ligands selects
  `Constitutional` alone. For tetrahedral atoms and cis/trans bonds these are the units that do
  not depend on others.
- The remover and the enumerator use the same selection. A remover asked for the records that are
  not stereogenic units leaves the three records of cis,trans-1,3,5-trimethylcyclohexane in
  place.
- The enumerator evaluates asserted stereo constraints on each candidate with the readings of
  the stereo views, and builds one `GraphSymmetry` per candidate for that.

### Computation

Every reading is decided by the search of [245](245-canonical-eq-cost-2026-11-12.md) for a
remapping of the molecule onto itself that passes a check. Source and target use the same
Constitution graph, whose colors come from S4f's shared calculation with the selected constraints
setting. Stereo entities have no nodes in this graph, so their kinds, configurations and
constraints do not refine its cells. Those selected stereo forms are checked at the leaves,
using the image record found at the mapped site. The refined cells do not depend on stored
configurations.

| Reading | Remapping searched for | Result |
| --- | --- | --- |
| stereogenicity of X | preserves every configuration except that of X, which it carries to another coset | found for every other coset: X is not `Stereogenic` |
| topicity of a and b at X | fixes X, sends a to b, and preserves every configuration | found: `Homotopic` |
| | fixes X, sends a to b, and carries every configuration to its mirror image | found: `Enantiotopic` |
| | fixes X and sends a to b, configurations not compared | found: `Diastereotopic`; none: `Distinct` |
| ligand symmetry of a permutation at X | fixes X, sends each ligand to its image, and preserves every configuration (proper permutation) or carries every configuration to its mirror image (improper permutation) | found: true |
| `is_chiral` | carries every configuration to its mirror image | found: the molecule is achiral |
| same proper orbit of a and b | sends a to b and preserves every configuration | found: true |
| same star orbit of a and b | sends a to b and either preserves every configuration or carries every configuration to its mirror image | found: true |

- **Mirror image.** The coset space of each kind supplies the mirror image of a coset. For a
  tetrahedral or axial configuration it is the other coset; for a cis/trans configuration it is
  the same coset. The graph check composes the new mirror action with the remapping's frame
  action before applying it once to the source coset. Check that the kind and frame degrees
  agree before composition, which otherwise panics. Constraints receive only the frame action.
- **Application.** `coset_apply_permutation` evaluates the new action at the outer
  `StereoCoset` level. It transports literals and set members, preserves `Undetermined` and
  the set variant, and rejects an incompatible action or an out-of-range literal. For a `Term`
  an action that is the identity on the kind's cosets returns the supplied term unchanged;
  otherwise it constructs `Apply` around that whole term. It does not enter existing operators,
  normalize a variable domain or call `canon_coset` on the supplied expression. Combining newly
  requested actions before application preserves the raw operand when their product is the
  identity. Two separately constructed symbolic applications may have a different expression
  tree with the same value; graph checks do not normalize that tree for comparison.
- **Constraints.** With `include_constraints` unset, each check compares the configurations of
  the stereo entities and not their constraints. With it set, each check compares the whole
  transported stereo forms. These graph-operation comparisons are structural. The existing
  canonicalization checks retain `normalized_eq` and their normalized-input caller contract;
  sharing colors and transport does not replace that relation.
- **Stereogenicity.** The check takes the coset that X is carried to. The search runs once per
  other coset of the kind. X is `Stereogenic` at the first coset for which no remapping is
  found, and not stereogenic if one is found for every coset. A kind with two cosets needs one
  search. For a coset of X that is not a literal the cosets are compared with the first one.
- **Generator pruning.** The target automorphisms kept by the search must preserve the
  particular query predicate. Preserving generators alone do not guarantee this for a
  coset-change query whose focus can move or whose kind has more than two cosets. For these
  queries, retain a preserving generator only if it also fixes the focus record and preserves
  the requested target coset under the induced ligand-frame action. Reject a generator if its
  focus record, kind or frame cannot be read or transported. This is the selected correction
  from the generator-filter comparison under S5b. Clearing every generator is a correct
  alternative with the same query semantics, but loses this pruning. Bypassing the filter
  would instead keep all generators. Proper and mirror queries retain their preserving filter.
- **Prostereogenic.** For a stereo entity that is not stereogenic, the same decision is repeated
  once per explicit ligand with that ligand held fixed. If X comes out stereogenic for some
  ligand, it is `Prostereogenic`. Otherwise it is `Nonstereogenic`.
- **Topicity.** The three searches run in the order of the table, and the first that finds a
  remapping decides. All three use the prepared Constitution graph. The first two check the
  selected stereo forms at the leaves; the third omits that predicate. It prepares no second
  graph.
- **Fixed and sent atoms.** An atom is kept fixed, or sent to another atom, by partitions in
  which it is a single cell, one partition for the source and one for the target. The existing
  `OrderedPartition::individualize` makes an atom a single cell, once per atom.
- **Virtual ligands.** A remapping keeps every virtual ligand at its position, which is the rule
  of the data types guide when the stereo record is fixed. Holding a virtual ligand fixed
  therefore changes nothing. Two equal virtual ligands at one site are excluded by molecule
  integrity. A virtual ligand and a ligand of another kind are `Distinct`.
- **Undetermined cosets elsewhere.** An entity with an undetermined coset is carried onto one
  with an undetermined coset. The class is then a statement about the description as stored: it
  says whether another coset of X changes what the description covers. In
  1,3-dichlorocyclobutane with one center undetermined, the other center is not stereogenic.
- **No assigned configuration.** A molecule without one is achiral.
- **Records that cannot be carried.** A check fails for a stereo record that cannot be carried
  onto its image, among them a record whose coset is outside the range of its kind, as
  `remapping_preserves_stereo` does.

The exits come before any search.

1. If every refined cell holds a single node, the identity is the only remapping. Every stereo
   entity with a kind is `Stereogenic`, and no two atoms share an orbit.
2. If the ligands of X lie in pairwise different cells, X is `Stereogenic`. A remapping that
   carries X to another coset has a power that fixes X and permutes its ligands, which needs two
   ligands in one cell. The argument does not use the number of cosets.
3. Two atoms in different refined cells are not related by a remapping. Two ligands in different
   cells are neither `Homotopic` nor `Enantiotopic`, and two ligands in different cells of the
   Constitution level are `Distinct`.

The new code is private and carries the root of the functions beside it. The search is not
changed.

| Name | Decides or gives |
| --- | --- |
| `remapping_preserves_stereo`, existing | the normalized stereo forms are preserved for canonicalization |
| `remapping_mirrors_stereo` | every configuration goes to its mirror image |
| `remapping_changes_coset` | one stereo entity goes to a given coset, and the others are preserved |
| `remapping_frame_action` | the permutation of a record's ligands under a remapping that sends the record to itself |
| `molecule::symmetry::remapping_exists_at` | a remapping at a description level and of an `Orientation` sends given atoms to given atoms |
| `molecule::symmetry::remapping_exists_for` | a remapping that keeps given atoms fixed changes the coset of a stereo entity to a given one |
| molecule::symmetry::molecule_remapping_exists | shared atom validation, partition preparation and search for the two entry points |
| molecule::symmetry::remapping_matches_stereo | the shared stereo walk, selected by Orientation and an optional entity/coset substitution |
| molecule::symmetry::entity_image | the checked image of an entity under the incidence-node mapping |

The two search entry points take the molecule and its auxiliary symmetry data and delegate
to the private molecule::symmetry::molecule_remapping_exists. That function calls the existing
generic ir::symmetry::remapping_exists; their names differ without relying on module paths.
The stereo checks share the private molecule::symmetry::remapping_matches_stereo walk with its
existing Orientation and optional entity/coset arguments. No new selector type is introduced.
Keep molecule::symmetry::entity_image as a private function with its existing signature. Its
single consumer uses it at several mapping sites; relocating it into a closure would merely
move the same helper. Their graph checks compare transported representations directly. The
canonicalization consumer continues to use its normalized comparisons. The molecular operations
and the stereo views supply their borrowed molecule; the auxiliary data does not supply one.

### Modules

S1 placed each type's symmetry and canonicalization beside one another. Symmetry belongs to
the molecule, and canonicalization is one of its users. The adapter belongs with the views and
forwards to molecule operations.

| Module | Holds | Before S1 |
| --- | --- | --- |
| `ir::entity` | `DescriptionLevel`, beside `EntityKind` | `ir/canonicalize.rs` |
| `ir::symmetry` | the generic search: the colors' ranking, the adapter, the partition and its refinement, the search, `remapping_exists` | `ir/canonicalize.rs` |
| `ir::canonicalize` | what is independent of the type: the `Canonicalize` trait, the context, the choice of level, the errors | `ir/canonicalize.rs` |
| `ir::molecule::symmetry` | the molecule's colors and remapping checks, `GraphSymmetry` and its readings | `ir/canonicalize/molecule.rs`, `ir/symmetry.rs` |
| `ir::molecule::canonicalize` | the molecule's leaf keys, canonical form and `Canonicalize` impl | `ir/canonicalize/molecule.rs` |
| `ir::reaction_span::symmetry` | the reaction span's colors and remapping checks | `ir/canonicalize/reaction_span.rs` |
| `ir::reaction_span::canonicalize` | the reaction span's leaf keys, canonical form and impls, with the `Reaction` impl that calls them | `ir/canonicalize/reaction_span.rs` |
| `ir::view::symmetry` | `GraphSymmetryView`, borrowing a molecule and owning or borrowing `GraphSymmetry` | added by S5b |

A description level is a cumulative group of entity kinds, and the readings and canonicalization
are two of its users. `DescriptionLevel` therefore sits with the entity kinds and belongs to
neither. The symmetry modules import nothing from the canonicalization modules.

The colors and remapping checks of a type are `initial_colors`, `incidence_color`,
`remapping_preserves_stereo` and `retain_stereo_preserving_generators`. S1 moved them unchanged
with their tests. The remaining color work changes `AsColor` and the shared color construction;
the canonicalization checks retain their normalized comparisons. The graph-operation checks
read the molecule supplied by their callers and use structural comparisons after shallow
application. Search machinery remains independent of both consumers.

### Removed

| Item | Replaced by |
| --- | --- |
| `StereoSymmetry`, `Molecule::stereo_atom_symmetry`, `Molecule::stereo_bond_symmetry` | the three derived readings on the stereo views |
| the five molecule-wide query methods on `GraphSymmetry` | the same query names on `Molecule` and `GraphSymmetryView` |
| the present content of `GraphSymmetry` and its derivation | the value described above |
| the view methods that take a `&StereoSymmetry`, among them `ligand_symmetry` returning the local group | the three derived readings |
| the view conveniences `is_stereogenic`, `is_prostereogenic`, `is_homotopic`, `is_enantiotopic`, `is_diastereotopic` | comparison with the variant |
| `MoleculeColoring`, `ConstitutionColoring`, `MoleculeColoringFeatures` | `AsColor`, the one definition of a color, with `include_constraints`; the field switches have no replacement |
| `StereoModel::para_stereo` | nothing; selecting `Constitutional` alone leaves the units that depend on others out |
| `para_stereo_max_iterations` of the stereo operation configs and perceiver methods | `max_assignments` |
| `StereoSelection::stereogenicities` | `StereoSelection::stereogenic_units` |

The search has no iteration to switch on. A pseudoasymmetric center and a ring cis/trans center
are `Configurational` units, and in a molecule whose other configurations make them so they have
the class `Stereogenic`.

### Consequences

- **`umol-graph`.** `StereoPerceiver::classify`, the conformance validator and the enumerator
  build one `GraphSymmetry` per molecule and read the stereo views. `StereoValidateConfig`,
  `AddStereoConfig`, `RemoveStereoConfig` and `EnumerateStereoConfig` carry `max_assignments` in
  place of the iteration limit, and `StereoModel` loses `para_stereo`. In a molecule whose other
  configurations make them so, pseudoasymmetric and ring cis/trans centers have the class
  `Stereogenic`. The operations get no error from `graph_symmetry`. In a malformed molecule
  every recognized site reads `Constitutional`, and the default adder keeps all of them.
- **Stereo transformers.** The adder, the remover and the enumerator select by `StereogenicUnit`,
  as described above. Results of the default adder change for compounds with configurational
  units.
- **`umol-io`.** The table reader names the two renamed variants.
- **DSL.** The EDN keywords are `:` and the kebab case of the variant: `:prostereogenic`,
  `:nonstereogenic` and `:distinct`. In the compact attribute strings the glyphs `=`, `'` and `/`
  stand for `Nonstereogenic`, `Prostereogenic` and `Stereogenic`, and `Distinct` has the fourth
  topicity glyph `%` beside `=`, `'` and `/`.
- **Python.** The `Stereogenicity` and `Topicity` enums, `StereoSelection`, the stereo operation
  configs and `StereoModel` follow the Rust changes, and `StereogenicUnit` is added. The Python
  stereo constraints views have the keyed readings, so Python gets `GraphSymmetry`,
  `GraphSymmetryConfig`, `GraphSymmetryView`, `Molecule.graph_symmetry`, `with_symmetry`, and
  the four derived readings on the stereo views. A derived reading without the value raises
  `ValueError`, as a ring key does.
- **Guides.** `nomenclature.md`, `data-types.md`, `operations.md`, `python-api.md` and
  `docs/chemistry/representation.md` name the classes `Symmetric` and `Prochiral`, describe the
  para-stereo setting, and describe the coloring features.
- **Tests.** The tests of `symmetry.rs`, among them the fifteen cases of
  `test_stereo_symmetry_is_stereogenic`, state the rule of the local group, and two of those
  cases state it for a kind with more than two cosets. The property tests under
  `tests/property/stereo`, `tests/canonicalization.rs` and `tests/integrity.rs` use
  `GraphSymmetry`; `tests/canonicalization.rs` reads orbits with isotopes ignored, and
  `tests/integrity.rs` states that symmetry cannot fail. Several tests assert `Diastereotopic`
  for two different ligands. The compounds of the measurements are the first fixtures with two
  or more stereo entities.

### Cost

One `GraphSymmetry` per molecule carries one Constitution incidence graph, its colors and
refinement. It does not pay for molecule copying or normalization. A view borrows the molecule and either owns or borrows
that auxiliary data. Extracting owned data moves it; extracting borrowed data clones the auxiliary
data. A reading costs an exit by cell or one or more searches. The stereogenic units are computed
together and cached.

The timing trials below used a setup that retained a normalized molecule. Their measurements
describe those trials, not the cost of the revised setup and checks. The correctness trial of
shallow application below did not measure time. Required benchmarks precede the changes in S4a,
and comparisons follow the shared-color work and the consumer migration. Setup-only timings are
reported separately: the old setup precomputes chirality and orbits, while the new setup prepares
queries. The comparison also measures a fixed build-and-read workload.

The later search-layout comparison in S5b measured the selected setup directly. One Constitution
search was faster on every measured workload and retained 43–59% less Rust-owned storage
initially, and 63–70% less after the constitutional fallback. The final-candidate counts were
identical on those fixtures. The record distinguishes the primitive constitutional query from
the complete topicity method, and Rust-owned storage from native allocations.

The separate generator-filter comparison under S5b measures one coset-change query with prepared
data. Retaining the generators that also fix the focus and requested coset costs 5.2–7.1% in
six selected queries where it removes no candidates, and is about ten times faster in two
queries with symmetry away from the focus. This is a defensible selected tradeoff, not evidence
that it is the best choice across workloads. The full operation's slowdown is expected to be
smaller because the extra filter work is only part of that operation; it has not been measured.

Canonicalization keeps its input normalization and normalized comparisons. Its color ordering
and results must remain unchanged, and its benchmarks are checked before proceeding from the
shared foundations to the symmetry rewrite. No cost claim about the encoding-failure fallback is
made without those measurements.

## Trial

The trial ran in a scratch worktree and left the search unchanged. It added to the
canonicalization module

- a check that compares as `remapping_preserves_stereo` does, after replacing the literal coset of
  one stereo entity by the other coset or of every stereo entity by its mirror image,
- partitions in which chosen atoms are single cells, which keep an atom fixed or send one atom to
  another, and
- readings built on these: the class of a stereo atom and of a stereo bond, the topicity of two
  ligand atoms, `is_chiral`, and the two orbit queries for a pair of atoms.

The readings took literal cosets and explicit ligand atoms. The normalized molecule, its incidence
graph, its colors and its refined cells were held in one value per molecule, called the setup
below.

### Results

Six tests pass. They cover

- 1,2-dichloro-1,2-difluoroethane, meso and chiral, under both atom numberings: `is_chiral`, the
  orbits of the two carbons, and the classes;
- 2,3,4-trichloropentane, both meso forms and the chiral form: C3 is `Stereogenic` in the meso
  forms and `Prostereogenic` in the chiral form, and its two branches are enantiotopic and
  diastereotopic respectively;
- 1,3-dichlorocyclobutane, cis and trans: both centers are `Stereogenic`, the molecule is achiral,
  and the two centers share a proper orbit;
- 1-chloro-1-fluoropropane, the molecule of issue 43: the CH₂ carbon is `Prostereogenic` with
  diastereotopic hydrogens beside the assigned center and with enantiotopic hydrogens without it;
- single centers with four, three and two different ligands, and stereo bonds with the three
  classes.

### Timing

Medians in microseconds. "Today" is `graph_symmetry`, the class of every stereo entity, and
`is_chiral`. The topicity columns read every pair of explicit ligands of every stereo atom.

| Molecule | Today | Today, fixpoint | Setup | Classes, one setup | Classes, setup per call | `is_chiral` | Topicity today | Topicity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `tetrahedral_stereo` | 4.7 | 5.0 | 4.1 | 4.1 | 4.1 | 7.3 | 4.7 | 6.0 |
| `cis_trans_stereo_bond` | 5.0 | 5.0 | 4.5 | 4.5 | 4.5 | 8.2 | 3.4 | 4.7 |
| `mixed_atom_and_bond_stereo` | 9.3 | 9.2 | 8.2 | 8.3 | 16.7 | 16.1 | 7.4 | 12.5 |
| `meso_dichlorobutane` | 7.9 | 9.3 | 5.4 | 5.4 | 10.8 | 23.0 | 8.0 | 7.7 |
| `para_stereo_trichloropentane` | 11.4 | 13.2 | 8.4 | 45.5 | 61.9 | 34.3 | 11.6 | 32.0 |
| `para_stereo_cascade` | 70.6 | 78.0 | 25.9 | 581.6 | 806.7 | 299.5 | 71.1 | 861.5 |
| `symmetric_stereo_triisopropanolamine` | 22.7 | 22.7 | 9.3 | 9.3 | 27.6 | 50.5 | 22.7 | 14.0 |
| `symmetric_stereo_butanol` | 24.8 | 24.9 | 9.5 | 9.6 | 9.5 | 79.7 | 24.9 | 16.6 |
| `symmetric_stereo_dimethylbutanol` | 68.9 | 68.8 | 11.9 | 12.2 | 12.0 | 196.4 | 68.5 | 21.5 |
| chain, ten centers, substituents F, Cl, Br in turn | 53.5 | 54.0 | 31.4 | 32.8 | 314.7 | 71.4 | 53.9 | 44.6 |
| chain, ten centers, substituent Cl | 54.8 | 60.7 | 28.7 | 29.2 | 288.1 | 130.2 | 56.0 | 43.4 |

The named rows are stereo fixtures of the canonicalization benchmark. The two chains are a carbon
chain of twelve atoms with one substituent and one implicit hydrogen on each inner carbon.
`para_stereo_cascade` is a synthetic fixture with ten stereo entities on four shared ligand atoms,
among them axial and square-planar kinds, for which the trial's replacement of a coset is not
defined; its row shows the order of magnitude only.

### Conclusions

- **The search carries the readings.** The checks and the fixed or sent atoms need no change to
  `remapping_exists` or to the search below it, and the results do not depend on atom numbering.
- **Classes with one setup per molecule cost less than today where an exit applies**, between
  0.2 and 0.9 of today's time. Where a search runs, as in `para_stereo_trichloropentane`, they
  cost four times today's time. One search costs about as much as the setup.
- **The setup per call grows with the number of stereo entities.** At ten centers the classes
  cost between 5.3 and 5.9 times today's time. An operation that reads every stereo entity of a
  molecule needs the setup once per molecule.
- **`is_chiral` costs the setup and one search.** For a chiral molecule the search runs to the
  end. It costs between 1.3 and 4.2 times today's whole computation.
- **Topicity needs exits by cell.** Two ligands in different refined cells are not related by a
  remapping, and two ligands in different cells of the Constitution level are `Distinct`. Without
  these exits every pair ran three searches and cost 10 to 60 times today's time; with them it
  costs between 0.3 and 2.8 times.

The trial did not run the operations, the ligand symmetry reading, virtual ligands, or
undetermined cosets.

### Second trial: stereogenic units

The second trial ran in the same worktree. It tested the selection of the adder against the two
papers.

- **Reference.** Every assignment of the cosets of a compound was enumerated and the assignments
  were grouped by canonical form. A site counts as a stereocenter when its coset changes the
  group for some assignment of the others.
- **Compounds.** Thirteen compounds of Table I of Razinger et al. with their counts, and the
  structures a, c to k and n of their Figure 1 with the classification stated in the text. Each
  carbon with four ligands, at most one of them an implicit hydrogen, and each C=C bond with two
  ligands at each end was a candidate.

| Check | Outcome |
| --- | --- |
| counts of stereoisomers of the reference against Table I, total, chiral and achiral | all 13 compounds agree |
| stereocenters of the reference against the paper | all 21 structures agree |
| steps 1 and 2 of the procedure against the paper | 20 of 21; structure k keeps its central carbon |
| steps 1 to 3 against the paper | all 21 agree |
| the four variants on seven compounds | as expected on all seven |

Among the counts are inositol with 9 stereoisomers, 2 chiral and 7 achiral, compound I of the
paper with 10, 6 and 4, and its double-bond variant with 7, 4 and 3. Step 3 was needed in 9 of the
21 structures and tried at most 11 assignments; structure k tries all 8 and drops the central
carbon.

Step 2 as written runs one search per record per round, without step 1 as an exit. Medians in
microseconds:

| Molecule | Classes today | Setup | Step 2 |
| --- | ---: | ---: | ---: |
| `tetrahedral_stereo` | 4.7 | 4.1 | 7.2 |
| `cis_trans_stereo_bond` | 5.0 | 4.5 | 8.0 |
| `mixed_atom_and_bond_stereo` | 9.3 | 8.2 | 22.9 |
| `meso_dichlorobutane` | 7.9 | 5.4 | 38.6 |
| `para_stereo_trichloropentane` | 11.4 | 8.4 | 83.4 |
| `para_stereo_cascade` | 70.6 | 25.9 | 2698.6 |
| `symmetric_stereo_triisopropanolamine` | 22.7 | 9.3 | 232.4 |
| `symmetric_stereo_butanol` | 24.8 | 9.5 | 78.1 |
| `symmetric_stereo_dimethylbutanol` | 68.9 | 11.9 | 195.6 |
| chain, ten centers, substituents F, Cl, Br in turn | 53.5 | 31.4 | 439.4 |
| chain, ten centers, substituent Cl | 54.8 | 28.7 | 1158.1 |

Conclusions:

- **The representation and canonical equality reproduce the published stereoisomer counts.** The
  reference is a sound yardstick for the adder followed by the enumerator.
- **The published definition is decided by the search and the class reading.** No analysis of
  rings by rule is needed.
- **Step 1 is needed as an exit.** Without it step 2 costs between 1.5 and 21 times today's
  classification on the chemical fixtures. With it a record with constitutionally different
  ligands needs no search.
- **Step 3 is exponential in the worst case** and small in the structures of the paper.

Step 3 as stated in the design takes the assignments over fewer records than the trial first
did: over the kept records that do not sit alone in their cell when the record under test is
held fixed, and not over every other kept record. This was run afterwards on the same 21
structures. All agree with the paper, and the numbers of assignments are the same, so these
structures do not show the saving.

The second trial did not run a molecule in which the two sets of records differ, records with
undetermined cosets, cumulated double bonds, kinds with more than two cosets, the
structures l and m of Figure 1, the other 56 compounds of Table I, or the adder and the enumerator
themselves.

### Third trial: shallow application and comparison

A standalone correctness trial used the workspace crates at `b98caac18227` as dependencies,
without a repository copy or production edits. All 39 cases passed. It enumerated complete small
candidate remapping sets directly, without search pruning; its non-stereo entity domain was
atoms and localized bonds.

The trial evaluated literals, sets and undetermined cosets at their outer form. For a symbolic
term it preserved the whole operand and recorded only the new action; a coset-identity action
preserved the operand exactly. Graph checks composed a new mirror and frame action before
application and compared representations directly. Canonical comparison controls kept
`normalized_eq` and normalized inputs.

| Check | Evidence and result |
| --- | --- |
| chirality and complete proper/star orbit vectors | all 12 molecule cases passed, including a meso/chiral pair with two centers and both cis/trans configurations |
| literal, undetermined and selected finite-set actions, all six kinds | 27,768 inverse checks and 222,144 bounded composition checks passed |
| canonical comparison after transport, with and without mirroring | all 115,712 decisions agreed with the current normalized-comparison path |
| ignored and included raw constraints | excluded contradictions and invalid constraint positions were not read; included payloads were compared structurally |
| accepted symbolic remappings | proper and star group checks passed on 1,000 bounded three-record assignments, including 2,600 group products |

The 12 molecule cases covered literal, undetermined, finite-set, singleton-set and symbolic
chlorofluoromethane configurations; contradictory atom constraints under both option settings;
a stored singleton element versus a literal; meso and chiral dichlorodifluoroethane; and both
cis/trans stereo-bond configurations. They checked complete orbit vectors and retained the
source molecules structurally. The raw element example deliberately gave different symmetry
readings before and after explicit normalization, while canonical forms agreed.

For composition, every permitted first action was paired with at most the first eight second
actions in permutation order. Canonical comparison used eight normalized input rows per kind,
including symbolic expressions, under every permitted frame action, with and without mirroring.
The 115,712 agreements establish local comparison parity; they do not represent a rebuilt
canonicalizer using the prototype. The group checks used ten term shapes, all six permutations
of the three records and the specified compatible frames. These are bounded checks rather than
exhaustive checks of every expression or pair of actions.

An initial interpretation contracted an operator already inside the supplied term. It lost raw
syntax in `Mirror(Mirror(x))` across an action and its inverse. That interpretation was rejected.
The final checks assert the complete newly constructed expression trees and separately check
their value equivalence; exact operand preservation is checked when the new actions compose to
the identity. In particular, separately applying an action and its inverse to
`Mirror(Mirror(x))` constructs
`Apply(Apply(Mirror(Mirror(x)), action), inverse)`. It preserves the original operand as a
subexpression and is value-equivalent to it. Composing those two actions before application
preserves the original expression exactly. This does not impose normalization on graph inputs
or alter the public `FrameTransport` contract.

The trial checked the application and comparison decisions using direct enumeration. It did
not implement the new `GraphSymmetry`, its view or the shared `AsColor` path, and it measured no
timings. The implementation tests and measurements below precede those changes.

## Review dispositions: S4 and unfinished S5

The supplied review is evidence for discussion, not a blanket approval of its findings or
proposed corrections. The identifiers below refer to that review. This record distinguishes
changes verified in main, behavior checked only in the trial, remaining implementation work,
and design questions that still need a decision. S4f completes the selector replacement; it
does not close every S4 finding. Selecting one Constitution search likewise does not complete
S5b or approve all of its private machinery. The table records the findings and work assigned
before implementation: the requirements assigned to S4g and S5b are now complete, as recorded
in the completion checkpoint under S5b. Coverage and consumers assigned to S5c and later
subitems remain pending.

| Review items | Current finding and disposition | Remaining work |
| --- | --- | --- |
| P1, B1, B3, N1–N3 | S4f removed ColorProjection, its marker types and constant, and both generic dispatch bodies. Full colors remain reference-sized and the ranker remains generic. | Closed for the selector mechanism. No selector or replacement ranking abstraction is added by the S5 port. |
| P2, P3 | The operative design and S4e approval wording have been corrected. S4c now explicitly labels the selector as an initial implementation superseded by S4f. The earlier mechanism was not approved by accepting the reference-sized color outcome. | Keep completion claims scoped to their verified subitems; the remaining S4/S5 items below are open. |
| B2 | Both molecule and reaction-span constraints flags were retained by the accepted S4f design. Current reaction-span production callers pass true. | No new span consumer is claimed or added. This explicit design decision closes the flag question for the current work. |
| B4 | The proposed stereo-bond frame-admissibility helper and the analysis offered to justify it were rejected. The user then dropped this review item. | Dropped. Do not introduce the proposed helper or treat this review finding as an approved consolidation task. |
| P4, B5, N5 | The selected representation inlines incidence, adapter, partition and cell_indices into GraphSymmetry, with config last. Retain the node-to-cell index and reuse OrderedPartition::cell_indices with pub(crate) visibility. | Remove SymmetrySearch and the second-search cache. Prepare the fields in Molecule::graph_symmetry and move the private cell method onto GraphSymmetry. These are settled choices. |
| B6, N6 | Shared search setup and the shared stereo walk are retained. Their approved private names are molecule::symmetry::molecule_remapping_exists and molecule::symmetry::remapping_matches_stereo. Keep the existing Orientation/optional-substitution interface and the private entity_image function, including its name and signature. | Apply the two renames and retain the shared implementation. Do not add selector types or replace entity_image with a closure. The generic ir::symmetry::remapping_exists keeps its existing name. |
| B7 | The trial reproduced the degree-mismatch panic under both constraint settings, then passed the cases with a degree check before mirror composition. That correction is accepted. Main still needs the guard. | Put the two exact encoding-failure cases in the owning main tests before porting the correction. The public operation remains infallible and preserves its input. |
| B8 | Explicit mappings now demonstrate both failures of predicate invariance: a preserving generator can move the focus, or fix its stored coset while moving the requested coset. The two safe filters have identical query semantics and each passes all 261 trial molecule-symmetry cases. The narrower filter is selected on the qualified cost/complexity grounds recorded under S5b; the comparison does not establish that it is the best choice overall. | Put the explicit invariance cases and exact coset answers in the owning tests before porting the narrower filter. Extend production reading coverage in S5c. Main still has the old filter; no wrong public answer from that backend has been reproduced. |
| B9 | The single-search decision already eliminates the provisional second-search cache. Its initialization, contents and pointer-identity assertions describe the rejected implementation. | Remove those assertions and the duplicate call while retaining the Constitution-query assertion when integrating the selected layout. This is a consequence of the settled decision, not a new approval question. |
| B10 | All nine form/span tests compute expected_color by calling as_inherent_color on another form. The test-writing skill requires independent expected values; this is a repository requirement, not an unresolved design choice. | In S4g, use explicit expected color fields and field-distinguishing cases while preserving the existing exclusion, hashing and input-preservation requirements. The exact-rank cases remain. |
| B11 | Module grouping and the test-writing skill determine the ownership and order of these tests. The four color-size controls currently live in the generic search tests, and GraphSymmetryView tests live in the molecule tests. | Move size controls beside their color types in S4g. Move view tests to the view module in S5b and order tests with their subjects. Do not add a production seam merely to share fixtures. No separate design decision is needed. |
| N4 | Retain AsColor as the trait. The approved associated type and method are FieldColor and as_field_color; the eight concrete names are <Kind>FieldColor. These supersede the briefly accepted <Kind>InherentColor names. | Apply the settled names consistently to forms, EntitySpan, aggregate color consumers and their tests in S4g. Full Color and as_color retain their names and semantics. |
| Correction 7 | Existing symbol names and the nomenclature guide supply symmetry for GraphSymmetry locals and config for operational choices. The module description must cover the colors and checks it owns. | Apply these narrow consistency and documentation corrections in S5b. They introduce no new concept or API and do not require another naming decision. |
| P5 | The supplied review identifies six tests written after the code. The recorded initial compile failure has 25 errors, including authoring errors involving push, TopicityForm::single and StereoLigandPair. S4f subsequently ran the whole current main library suite: 12,859 passed, three ignored. | Preserve the actual chronology. The earlier claim that some current library cases had no recorded run is outdated; passing them later does not make them test-first. The new degree cases did fail on the panic before their trial fix. |
| P6 | The old retirement list contained eight names without their assertions and omitted the changed stereo-view test. The ledger under S5b below now records the requirements and remaining production coverage. | Complete the pending S5c replacements and make the mixed stereo-view test's name, inputs and assertions agree. A ledger entry is not a claim that its replacement already exists. |
| B12 | The exhaustive class/topicity calculations establish values but do not yet compare those readings with production. The molecule symmetry cases already compare the calculation and the molecule/view queries with explicit expectations. | Preserve the accepted exhaustive calculations as ongoing correctness checks. S5c must compare its production classes and topicities with the established values; testing a calculation alone is not production coverage. |
| D: gates and consumers | Main's canonicalization/integrity tests and stereo property consumers still use removed APIs; Clippy and formatting failures remain recorded. No direct removed symmetry API use was found in the graph-IR benchmark sources during this pass. | Finish the named S5b consumers and its cleanup; migrate stereo properties with S5c. Build graph-IR benchmark targets to verify them rather than assuming every target needs source migration. The higher-level stereo benchmarks still migrate with S5e. |
| D: later consumers | remapping_exists_for and remapping_changes_coset are currently test-only; max_assignments is for the unit computation. | Enable the production consumers in S5c and S5d. Rust operations and Python follow in S5e/S5f; their absence here is planned sequencing, not a reason to add speculative paths. |

No further semantic, interface or naming decision has been identified for this reviewed S4 and
unfinished S5b scope. The explicit decisions above settle the representation and helper names;
the remaining test independence, ownership, ordering and terminology corrections follow the
repository rules. They are implementation requirements, not additional approval questions.
The stereo-bond proposal was dropped and is not pending. This does not claim that implementation,
coverage or the S5 gates are complete.

Test corrections follow the repository's current requirements: independent expected values,
owner placement and ongoing semantic checks. Temporary drivers, allocation counters and
unselected implementations stay in scratch. No historical-role label or duplicate algorithm is
retained in production or tests solely because it was used during this development.

## Staged implementation plan

S0-S3, S4a-S4f and S5a have implementation records. S4g applies the settled naming and test
corrections before the unfinished S5b port of the selected single Constitution search.
The remaining work uses one color definition through AsColor and
shared search-color construction. Each subitem carries its own tests under the test-writing
skill. Required tests, independent expected values and benchmark harnesses are
written before the implementation they exercise. An intended initial test failure
belongs inside that subitem; it does not justify leaving a nonbreaking subitem red. Where an
existing assertion changes or is removed, its old and new statements are shown before editing it.
Fixtures stay beside their tests and are added when approved.

A subitem is checked with the tests of the crates it touches. A subitem that moves code without
changing it gets one build with those tests. A subitem that changes an API a benchmark uses also
builds the benchmarks of that crate, which the tests do not compile. Python is rebuilt and
tested after every subitem that touches `umol-py`.

Formatting is `cargo +nightly fmt --all`. Every command that compiles or loads `umol-py`, including
a workspace Cargo command, activates `umol-py/.venv` in that same shell; Python 3.13 is checked
before the first such build. Baselines use the shared Cargo target directory. Source files are
not edited during a benchmark run. Existing `pre-246` measurements remain the baseline for their
unchanged workloads; newly added cases get their own baseline before the relevant code changes.

### S0 — Fixtures and baseline

#### S0a — Fixtures **Done**

- **Module:** the test modules of S5 and the benchmark of S0b. The fixtures have no file of
  their own: each molecule is a `#[fixture]` function or a case row beside the tests that use it.
- **Adds:** the compounds with two or more stereo entities and the two compounds with a kind of
  more than two cosets, each with explicit ligand atoms:
  - 1,2-dichloro-1,2-difluoroethane, meso and chiral, each under two atom numberings;
  - 2,3,4-trichloropentane, both meso forms and the chiral form;
  - 1,3-dichlorocyclobutane, cis and trans;
  - 1-chloro-1-fluoropropane with and without the configuration of C1;
  - cis- and trans-diamminedichloroplatinum, square planar;
  - chlorotetrafluorophosphorane, trigonal bipyramidal, chlorine axial and equatorial;
  - 1,3,5-trimethylcyclohexane, all-cis and cis,trans;
  - the thirteen compounds of Table I and the eleven structures of Figure 1 of Razinger et al.
    that the second trial used, with the published counts and classifications.
- **Kind:** no change to the tree. The molecules are checked in the scratch worktree and enter
  the tree with their tests and benchmark cases.

The 38 molecules were written and checked in the scratch worktree. Every one passes the integrity
check. The first group gives the readings of the trial. Canonical equality puts the three
square-planar cosets and the twenty trigonal-bipyramidal cosets into the two classes each that
the fixtures name. The brute-force reference gives the counts of Table I for the 13 compounds and
the stereocenters of the paper for the 21 structures.

#### S0b — Benchmarks of the operations and baseline **Done**

- **Module:** `umol-graph/benches/transform.rs`, `umol-graph/benches/transform/stereo.rs`,
  `umol-graph-ir/benches/canonicalize.rs`.
- **Adds:** cases for `StereoPerceiver::classify`, the conformance validator and the enumerator
  with five compounds of S0a, each chosen for the piece of the new code it exercises:

  | Compound | Exercises |
  | --- | --- |
  | 1-chloro-1-fluoropropane, with an asserted class on the CH₂ record | a class search that finds a remapping, the prostereogenic step, a record dropped by step 2, the validator's comparison, the enumerator's `GraphSymmetry` per candidate |
  | meso-2,3,4-trichloropentane | a class search among equal ligands that finds no remapping |
  | inositol | six equal centers decided by step 2 alone; an enumeration over 64 assignments |
  | compound I | step 3 confirming records, seven assignments in all |
  | structure k | step 3 trying all eight assignments and dropping the record |

  The validator has no benchmark today and gets one function per case beside `classify`. The
  timings of these cases and of the canonicalization groups of
  [245](245-canonical-eq-cost-2026-11-12.md) are saved as a criterion baseline, the record that
  the final comparison reads.
- **Kind:** additive (green). `[dep: S0a]`

The five compounds are cases of the group `stereo_recognition` with the functions `classify` and
`validate`. Two are cases of the group `stereo_enumerator`: 1-chloro-1-fluoropropane with an
asserted class, and inositol with the selection of the class `Stereogenic`. Clippy passes for all
targets of `umol-graph`.

Two results of the present code shaped the enumerator cases.

- With the class asserted on the CH₂ record of 1-chloro-1-fluoropropane the present enumerator
  returns no result, which is issue 43 in the enumerator. The enumerator case therefore asserts
  `Stereogenic` on C1 and gives four results. The case of `classify` and `validate` keeps the
  assertion on the CH₂ record.
- For inositol the present classification selects no site, and the enumerator returns the input
  as its one result.

The runs use `--warm-up-time 0.5 --measurement-time 1.5`, as in
[245](245-canonical-eq-cost-2026-11-12.md):

```text
cargo bench -p umol-graph --bench transform -- \
    --warm-up-time 0.5 --measurement-time 1.5 --save-baseline pre-246 \
    'stereo_recognition/(classify|validate)/|stereo_enumerator/.*/(chlorofluoropropane|inositol)'
cargo bench -p umol-graph-ir --bench canonicalize -- \
    --warm-up-time 0.5 --measurement-time 1.5 --save-baseline pre-246 'canonical_eq'
```

Criterion keeps both runs as the baseline `pre-246` in the cargo target directory: 51 benchmarks
of the operations and 215 of canonicalization. The first run includes `classify` on the cases
the group already had, which cover the readings that an exit decides.

The tables below hold the whole record: the median of each benchmark of the operations, and the
medians of the canonicalization benchmarks by group. All values are in µs and are read from the
stored record.

The five compounds:

| Compound | Case | `classify` | `validate` |
| --- | --- | ---: | ---: |
| 1-chloro-1-fluoropropane | `chlorofluoropropane` | 11.6 | 12.0 |
| meso-2,3,4-trichloropentane | `trichloropentane_meso` | 18.0 | 17.8 |
| inositol | `inositol` | 38.7 | 14.6 |
| compound I | `bis_fluoroethyl_difluorohexane` | 43.0 | 16.8 |
| structure k | `tris_hydroxyethyl_methane` | 26.0 | 10.6 |

The two enumerator cases:

| `stereo_enumerator` | `chlorofluoropropane` | `inositol` |
| --- | ---: | ---: |
| `transform_shared` | 10.6 | 42.6 |
| `transform_into_shared` | 10.8 | 42.1 |
| `transform_unique` | 10.5 | 40.0 |
| `transform_into_unique` | 10.4 | 40.0 |
| `iter_setup` | 0.50 | 39.9 |
| `iter_first` | 10.0 | 0.07 |
| `iter_prefix_4` | 39.1 | 0.12 |
| `iter_all` | 40.1 | 0.12 |

`classify` on the cases the group already had:

| Case | `classify` | Case | `classify` |
| --- | ---: | --- | ---: |
| `atoms_5` | 4.5 | `fixture_tetrahedral` | 4.7 |
| `atoms_12` | 7.7 | `fixture_tetrahedral_literal` | 4.5 |
| `atoms_36` | 20.9 | `fixture_cis_trans` | 4.9 |
| `atoms_132` | 101.6 | `fixture_cis_trans_literal` | 4.8 |
| `ring_4` | 0.005 | `fixture_implicit_hydrogen` | 4.1 |
| `ring_6` | 0.005 | `fixture_lone_pair` | 4.2 |
| `ring_8` | 8.2 | `fixture_imine` | 4.0 |
| `ring_12` | 10.6 | `fixture_small_ring` | 0.005 |
| `sites_1` | 244.0 | `fixture_evidenced_stereo` | 7.5 |
| `sites_2` | 307.0 | `fixture_remote_stereo` | 15.8 |
| `sites_4` | 450.9 | `fixture_prostereogenic` | 7.2 |
| `sites_8` | 762.4 | `fixture_nonstereogenic` | 22.7 |
| `sites_16` | 1611.7 | | |

The three cases at 0.005 µs have no stereo site to classify.

Canonicalization, lowest, median and highest of the medians of each group:

| Group | Benchmarks | Lowest | Median | Highest |
| --- | ---: | ---: | ---: | ---: |
| `canonicalize/complete/canonical_eq_equal` | 20 | 7.6 | 15.2 | 314.6 |
| `canonicalize/complete/canonical_eq_unequal_colors` | 20 | 3.6 | 5.8 | 20.4 |
| `canonicalize/complete/canonical_eq_unequal_constitution` | 3 | 5.6 | 8.2 | 18.6 |
| `canonicalize/complete/canonical_eq_unequal_counts` | 20 | 0.010 | 0.011 | 0.012 |
| `canonicalize/complete/canonical_eq_unequal_stereo` | 10 | 9.5 | 36.1 | 312.7 |
| `canonicalize/inline_constraint/molecule/canonical_eq_equal` | 3 | 7.9 | 12.5 | 13.1 |
| `canonicalize/inline_constraint/reaction_span/canonical_eq_equal` | 3 | 8.7 | 13.8 | 14.6 |
| `canonicalize/reaction/canonical_eq_equal` | 20 | 16.8 | 31.5 | 362.6 |
| `canonicalize/reaction/canonical_eq_unequal` | 20 | 12.2 | 19.6 | 60.1 |
| `canonicalize/reaction_span/canonical_eq_equal` | 20 | 8.4 | 16.3 | 327.0 |
| `canonicalize/reaction_span_modified/canonical_eq_equal` | 20 | 7.2 | 15.0 | 202.2 |
| `canonicalize/reaction_span_modified/canonical_eq_unequal_colors` | 20 | 3.4 | 6.5 | 24.2 |
| `canonicalize/reaction_span_modified/canonical_eq_unequal_counts` | 20 | 0.009 | 0.010 | 0.011 |
| `canonicalize/reaction_span_modified/canonical_eq_unequal_stereo` | 10 | 10.5 | 30.7 | 204.2 |
| `canonicalize/scaling/canonical_eq_equal` | 3 | 7.5 | 12.0 | 12.7 |
| `canonicalize/scaling/canonical_eq_unequal` | 3 | 0.011 | 0.011 | 0.011 |

The measured code is unchanged since these runs: S1 moved code, S2 renamed variants and added
one, and S3 added two functions compiled for tests only, tests, and an enum. The record can
therefore be taken again with the two commands above, on the tree as it is after S3c.

### S1 — Modules

Every subitem of this stage moves code without changing it.

#### S1a — Present symmetry code under the molecule **Done**

- **Module:** `ir/symmetry.rs` → `ir/molecule/symmetry.rs`.
- **Changes:** the file moves with its tests. The public paths `ir::GraphSymmetry`,
  `ir::GraphSymmetryConfig` and `ir::StereoSymmetry` stay.
- **Kind:** move (green).

Only import paths changed: those of the moved file and those of its one user,
`ir/view/stereo.rs`. `ir::molecule` re-exports the three types, as it does for its other
submodules. The tests of `umol-graph-ir` pass with the default features.

#### S1b — Generic search **Done**

- **Module:** `ir/canonicalize.rs` → `ir/entity.rs` and `ir/symmetry.rs`.
- **Changes:** `DescriptionLevel` moves to `ir/entity.rs`, beside `EntityKind`. `InitialColors`,
  `rank_colors`, `AutomorphismColor`, `AutomorphismAdapter`, `OrderedPartition`, `same_cells`,
  `remapping_preserves_incidences`, the canonical search and `remapping_exists` move to
  `ir/symmetry.rs` with their tests. `ir::canonicalize` keeps the trait, the context, the choice
  of level, the leaf-key helpers and the errors, and imports the search.
- **Kind:** move (green). `[dep: S1a]`

The listed items moved together with what belongs to them: `AutomorphismAdapterOutput`,
`source_orbits_from_generators`, and the types of the canonical search, `CanonicalSearchStats`,
`CanonicalCandidate` and `CanonicalSearchResult`. `ir::canonicalize` keeps `DeferredKey`,
`TopologySearchInput`, `color_representatives`, `molecule_remapping` and `hash_value`.

The moved code is unchanged. The items and fields that `ir::canonicalize` and the tests use are
visible within the crate, and rustfmt rewrapped three signatures that the visibility made longer
than the line. `ir/canonicalize.rs` differs from before only in its imports.

The tests that moved are `test_description_level_order`, to `ir/entity.rs`, and, to
`ir/symmetry/tests.rs`, `test_source_orbits_from_generators`, the five tests of
`OrderedPartition`, `test_canonical_search_stats` with its helper, and `test_rank_colors`. The
tests of the search that take molecules as input are in the molecule's test file and move with
it. The tests of `umol-graph-ir` pass with the default features, the same 12,443 unit tests as
before.

#### S1c — Colors and remapping checks of the two types **Done**

- **Module:** `ir/canonicalize/molecule.rs` → `ir/molecule/symmetry.rs`;
  `ir/canonicalize/reaction_span.rs` → `ir/reaction_span/symmetry.rs`.
- **Changes:** `initial_colors`, `incidence_color`, `remapping_preserves_stereo` and
  `retain_stereo_preserving_generators` of each type move with their tests.
- **Kind:** move (green). `[dep: S1b]`

The four functions of each type are unchanged apart from their visibility within the crate, and
the two canonicalization files differ only in their imports. `ir::molecule::symmetry` and the new
`ir::reaction_span::symmetry` are visible within the crate.

Six tests of each type moved unchanged, to `ir/molecule/symmetry/tests.rs` and
`ir/reaction_span/symmetry/tests.rs`. The fixtures `initial_color_molecule` and
`initial_color_span` moved with the tests of `initial_colors`; the canonicalization tests that
use them import them from there, and `ir/canonicalize/tests.rs` passes
`initial_color_molecule` on to the two test files below it until they move. The test module of
`ir/molecule/symmetry.rs`, which was inline, is now that file.

The tests of the search that take molecules as input stay with the canonicalization tests of the
molecule, because four of the seven use canonicalization helpers. The tests of `umol-graph-ir`
pass with the default features, the same 12,443 unit tests.

#### S1d — Canonicalization of the two types **Done**

- **Module:** `ir/canonicalize/molecule.rs` → `ir/molecule/canonicalize.rs`;
  `ir/canonicalize/reaction_span.rs` → `ir/reaction_span/canonicalize.rs`.
- **Changes:** the two files move with their tests. The tests of `ir/canonicalize/tests.rs`
  follow the functions they test.
- **Kind:** move (green). `[dep: S1c]`

The two files differ from before only in their imports. In `ir/canonicalize.rs` the items the two
files use are visible within the crate: the three functions that choose the level, `hash_value`,
`color_representatives`, `molecule_remapping`, `DeferredKey` and `TopologySearchInput`.

Of the tests of `ir/canonicalize/tests.rs`, the ten of the molecule's `Canonicalize` impl went
to `ir/molecule/canonicalize/tests.rs` with the four molecule fixtures and `reverse_remapping`,
and the eleven of the reaction span and the reaction went to
`ir/reaction_span/canonicalize/tests.rs` with their fixture. The file keeps the tests of what
`ir::canonicalize` holds, the test that compares the three types together, and what both types'
tests share: `canonicalize_context`, `permutations` and `explicitly_dense_equivalent`. No test
changed. The tests of `umol-graph-ir` pass with the default features, the same 12,443 unit
tests.

#### S1e — Comparison **Done**

- **Check:** the canonicalization groups against the baseline of S0b. No change is expected.
  `[dep: S1d]`

The 215 benchmarks were run with the settings of S0b and `--baseline pre-246`. The ratio of the
time after to the time before, from the medians, by group:

| Group | Benchmarks | Lowest | Median | Highest |
| --- | ---: | ---: | ---: | ---: |
| `canonicalize/complete/canonical_eq_equal` | 20 | 0.97 | 0.99 | 1.00 |
| `canonicalize/complete/canonical_eq_unequal_colors` | 20 | 0.94 | 0.98 | 1.00 |
| `canonicalize/complete/canonical_eq_unequal_constitution` | 3 | 0.98 | 0.98 | 0.99 |
| `canonicalize/complete/canonical_eq_unequal_counts` | 20 | 0.90 | 0.97 | 1.11 |
| `canonicalize/complete/canonical_eq_unequal_stereo` | 10 | 0.98 | 1.00 | 1.01 |
| `canonicalize/inline_constraint/molecule/canonical_eq_equal` | 3 | 0.98 | 0.98 | 0.98 |
| `canonicalize/inline_constraint/reaction_span/canonical_eq_equal` | 3 | 0.98 | 0.98 | 0.99 |
| `canonicalize/reaction/canonical_eq_equal` | 20 | 0.99 | 1.01 | 1.03 |
| `canonicalize/reaction/canonical_eq_unequal` | 20 | 0.96 | 1.00 | 1.03 |
| `canonicalize/reaction_span/canonical_eq_equal` | 20 | 0.98 | 0.99 | 1.01 |
| `canonicalize/reaction_span_modified/canonical_eq_equal` | 20 | 0.98 | 0.99 | 1.00 |
| `canonicalize/reaction_span_modified/canonical_eq_unequal_colors` | 20 | 0.99 | 1.00 | 1.03 |
| `canonicalize/reaction_span_modified/canonical_eq_unequal_counts` | 20 | 0.93 | 1.03 | 1.12 |
| `canonicalize/reaction_span_modified/canonical_eq_unequal_stereo` | 10 | 0.98 | 0.99 | 1.05 |
| `canonicalize/scaling/canonical_eq_equal` | 3 | 1.00 | 1.01 | 1.01 |
| `canonicalize/scaling/canonical_eq_unequal` | 3 | 0.98 | 1.01 | 1.02 |

- **Above one microsecond.** The 172 benchmarks have a median ratio of 0.99 and lie between 0.94
  and 1.05. The two highest, 1.05 and 1.03, are in
  `reaction_span_modified/canonical_eq_unequal_stereo`; a second run of that group gave 0.98 and
  0.99 for them.
- **About ten nanoseconds.** The 43 benchmarks that the entity counts decide lie between 0.90 and
  1.12, which is one nanosecond.

The moves changed no timing beyond the noise of a run.

The stage ends green. `umol-graph-ir` passes with all features, 12,443 unit tests and 559
property tests, and the crates that depend on it pass, 10,458 tests.

### S2 — Class names

#### S2a — `Stereogenicity` variants **Done**

- **Module:** `ir::stereo`, then every user in `umol-graph-ir`, `umol-graph`, `umol-io` and
  `umol-py`.
- **Changes:** `Prochiral` becomes `Prostereogenic` and `Symmetric` becomes `Nonstereogenic`.
  The keywords become `:prostereogenic` and `:nonstereogenic`, and the glyphs `=`, `'` and `/`
  stay with the renamed variants. The Python enum follows, and so does
  `umol-graph-ir/spec/umol-dsl-spec.md`, which states the keywords and the glyphs.
- **Kind:** breaking (red→green within the subitem).

The two variants are renamed in the enum, in the Python enum and at every use, 160 in Rust and
34 in the Python tests. The EDN keywords are `:nonstereogenic` and `:prostereogenic`; six test
strings spell them. The specification states the new keywords and the new words for the glyphs.

Wherever the words "prochiral" and "symmetric" named the two classes, they follow: 73 test case
labels, fixture names, benchmark case names, comments, the view method `is_prostereogenic`, and
one test constant in `umol-io`. The word "symmetric" stays where it means something else: a
relation, the symmetric group, or a molecule with many automorphisms, as in the
`symmetric_stereo_*` cases of the canonicalization benchmark. No test body changed. The guides
and `README.md` follow in S7a.

Two cases of the `classify` baseline changed their names, to `fixture_prostereogenic` and
`fixture_nonstereogenic`. Their record in `pre-246` was saved again under the new names, on code
whose `classify` is unchanged since the first run.

The tests pass: the Rust tests of the four crates, with all features for `umol-graph-ir` and
`umol-graph`, and 3,420 Python tests after a rebuild.

#### S2b — `Topicity::Distinct` **Done**

- **Module:** `ir::stereo`, `dsl::stereo`, the topicity forms, `umol-py`.
- **Adds:** the variant, its keyword `:distinct`, its glyph in the compact attribute strings, its
  Python member, and its entry in `umol-graph-ir/spec/umol-dsl-spec.md`. The glyph is chosen at
  this subitem.
- **Changes:** expected values of the lattice and normalization tests of `TopicityRelationForm`
  that enumerate all variants.
- **Kind:** breaking (red→green within the subitem). `[dep: S2a]`

The glyph is `%`. The variant is the last of the enum, its keyword is `:distinct`, and the
Python enum has the member `Distinct`. Nothing derives the value before S5c.

With four values a set of two has a complement of two, and normalization keeps the positive
side. The tests that wrote the whole domain or relied on a set of two becoming a complement
changed as follows, with H, E, D and X for the four variants.

| Test and case | Before | After |
| --- | --- | --- |
| `normalize`, `litset_polarity_to_notset` | `LitSet{H,E}` gives `NotSet{D}` | `LitSet{H,E,D}` gives `NotSet{X}` |
| `normalize`, `full_litset_to_undetermined` | `LitSet{H,E,D}` gives `Undetermined` | `LitSet{H,E,D,X}` gives `Undetermined` |
| `normalize`, `notset_complement_to_lit` | `NotSet{H,E}` gives `Lit(D)` | `NotSet{H,E,D}` gives `Lit(X)` |
| `normalize`, `full_notset_err` | `NotSet{H,E,D}` is a contradiction | `NotSet{H,E,D,X}` is a contradiction |
| `normalized_eq_default`, the full set among the values | `LitSet{H,E,D}` | `LitSet{H,E,D,X}` |
| `meet`, `notset_notset_to_lit` | `NotSet{D}` and `NotSet{E}` give `Lit(H)` | `NotSet{D,X}` and `NotSet{E,X}` give `Lit(H)` |
| `join`, `lit_union_to_notset` | `Lit(H)` and `Lit(E)` give `NotSet{D}` | renamed `lit_union_to_litset`: they give `LitSet{H,E}` |
| join of two topicity constraints, `topicity_widens` and `shared_topicity_widens` | `Lit(H)` and `Lit(E)` give `NotSet{D}` | they give `LitSet{H,E}` |
| meet of two topicity constraints, `shared_topicity_value_meet`, stereo atom and stereo bond | `NotSet{D}` and `NotSet{H}` give `Lit(E)` | they give `LitSet{E,X}` |
| `is_compatible`, `full_complement` | `Undetermined` with `NotSet{H,E,D}` is incompatible | with `NotSet{H,E,D,X}` |
| `is_compatible`, `complement_union` | `NotSet{H}` with `NotSet{E,D}` is incompatible | with `NotSet{E,D,X}` |

New rows state the tie and the new variant: `NotSet{H,E}` normalizes to `LitSet{D,X}`,
`LitSet{H,E}` is its own normal form, `LitSet{H,E}` joined with `Lit(D)` gives `NotSet{X}`, the
string `Th1#o(0,1)%` parses to `Lit(Distinct)` and renders unchanged, `:distinct` reads and
writes in EDN, the Python conversion round-trips, and the property strategies generate the
variant.

The specification states the glyph, the keyword, the two domain sizes and the tie. The tests
pass: 27,955 Rust tests of the four crates with all features for `umol-graph-ir` and
`umol-graph`, and 3,421 Python tests after a rebuild.

### S3 — Check and partitions

#### S3a — Checks with given cosets **Done**

- **Module:** `ir::molecule::symmetry`.
- **Adds:** the private checks `remapping_mirrors_stereo` and `remapping_changes_coset` beside
  `remapping_preserves_stereo`. They compare as it does, the first with the mirror image of
  every stored coset, the second with the coset of one stereo entity replaced by a given one.
  Virtual ligands keep their positions. A record whose coset is outside the range of its kind
  fails either check.
- **Tests:** each check on remappings of the S0a fixtures, among them a stereo bond, a center
  with an implicit hydrogen, a stereo entity with an undetermined coset, and a record with a
  coset outside the range of its kind.
- **Kind:** additive (green). `[dep: S1c, S0a]`

Both checks are in `ir/molecule/symmetry.rs`, right after `remapping_preserves_stereo`, which is
unchanged. Each repeats its loop over the stereo atoms and the stereo bonds.

- `remapping_mirrors_stereo` takes the mirror image of each form with the existing `mirror` and
  compares as `remapping_preserves_stereo` does.
- `remapping_changes_coset` gives the stereo entity its own coset as a source and the given
  coset as an image. Its own coset is the stored one, or the first coset of its kind when the
  stored one is not a literal. A stereo entity without a kind has no coset to change, and the
  check fails.

The two functions are compiled for tests only until the S5b-S5c readings call them.

`test_remapping_mirrors_stereo` has seven cases and `test_remapping_changes_coset` nine, on
1,2-dichloro-1,2-difluoroethane, 1-chloro-1-fluoropropane, dichlorofluoromethane with its
implicit hydrogen, but-2-ene, 1-bromo-1-chloro-2,2-difluoroethene and
bromochlorofluoroiodomethane. They cover an exchange of two equal ligands, an exchange of two
equivalent centers, the identity, a second stereo entity that must stay, an undetermined coset,
a coset outside the range of its kind and an entity without a kind. The tests of
`umol-graph-ir` pass with the default features, and clippy passes for all its targets.

#### S3b — Partitions with single cells **Done**

- **Module:** `ir::molecule::symmetry`, tests only.
- **Adds:** no function. A given atom is made a single cell with the existing
  `OrderedPartition::individualize`, once per atom, in the partition of the source and in that
  of the target of `remapping_exists`.
- **Tests:** an atom kept fixed and an atom sent to another, with `remapping_exists`, on the
  S0a fixtures; two atoms of different colors give no remapping.
- **Kind:** additive (green). `[dep: S1b, S0a]`

`test_remapping_exists_individualized` is in `ir/molecule/symmetry/tests.rs`, where the colors
and the checks of a molecule are. Each case gives a molecule, the atoms made single cells in the
partition of the source, and those in the partition of the target, in the same order. The call
is the one `canonical_eq` makes at the Structure level, with the molecule as source and as
target.

| Case | Molecule | Source | Target | Remapping |
| --- | --- | --- | --- | --- |
| `equivalent_atom` | 1,2-dichloro-1,2-difluoroethane, chiral form | C0 | C1 | yes |
| `mirror_image_atom` | 1,2-dichloro-1,2-difluoroethane, meso form | C0 | C1 | no |
| `atom_fixed_equivalent_atom` | trans-1,3-dichlorocyclobutane | C1, C0 | C1, C2 | yes |
| `atom_fixed_mirror_image_atom` | cis-1,3-dichlorocyclobutane | C1, C0 | C1, C2 | no |
| `different_colors` | 1,2-dichloro-1,2-difluoroethane, chiral form | C0 | F3 | no |

The tests of `umol-graph-ir` pass with the default features, and clippy passes for all its
targets.

#### S3c — `StereogenicUnit` **Done**

- **Module:** `ir::stereo`.
- **Adds:** the enum with its four variants. The rustdoc of `StereogenicUnit`, `Stereogenicity`
  and `Topicity` and of their variants cites the papers of the Sources section.
- **Kind:** additive (green). `[dep: S2a, S2b]`

`StereogenicUnit` is in `ir/stereo.rs` after `Stereogenicity` and is exported from `ir`. Its
variants are in the order `Constitutional`, `Configurational`, `Prostereogenic`,
`Nonstereogenic`. It derives what `Stereogenicity` derives except `VariantArray`, which nothing
needs for it.

The rustdoc of the three types links each paper by its DOI, in the form the workspace uses.

| Type | Cited |
| --- | --- |
| `Topicity` | Mislow and Siegel, Figure 1; "constitutionally heterotopic" at `Distinct` |
| `Stereogenicity` | Mislow and Siegel, footnote 24, p. 3326 and footnote 74; IUPAC 1974, Rule E-4.12, and IUPAC 1996 for "prochiral" |
| `StereogenicUnit` | IUPAC 1996 for "stereogenic unit"; Razinger et al. and Perdih and Razinger for the two kinds, with their terms at the variants; Mislow and Siegel, footnote 74 |

The five DOIs match the Crossref records of the papers. The unit tests of `umol-graph-ir` pass,
the documentation builds without warnings, and clippy passes for all targets of the crate.

### S4 — Shared colors and shallow coset application

This stage establishes the shared calculation before replacing its graph-symmetry consumer.
Every subitem ends green. Canonicalization retains its normalization and normalized stereo
comparisons; only its calls of the shared color constructor gain the explicit constraints flag.

#### S4a — Color controls and baselines **Done**

- **Module:** the form and `EntitySpan` test modules, `ir/molecule/symmetry/tests.rs`,
  `ir/reaction_span/symmetry/tests.rs`, the canonicalization tests,
  `umol-graph-ir/benches/canonicalize.rs` and `umol-graph/benches/transform/stereo.rs`.
- **Tests and benches first:** establish exact color-order and canonical-result controls before
  changing color types. Cover all eight forms and spans, entity constraints, nonuniform aromatic
  and multicenter electron counts, and normalizable raw expressions. Add cases to the existing
  canonicalization benchmark for these inputs before changing their implementation. Use exact
  expected values in the existing test modules; the current graph-symmetry answers are not the
  oracle for a known defect. Extend the existing `stereo` benchmark group with fixed workloads
  that build the data and query the molecule or its stereo records on the approved S0a
  molecules, including raw and constrained forms, before changing the shared coset helper.
- **Baselines:** retain S0b's 215 canonicalization and 51 operation cases under `pre-246`. Save
  the added canonicalization cases under `pre-246-colors` and the fixed reading workloads and
  existing setup/local-reading cases under `pre-246-readings`, using S0b's settings. Record the
  filters, input manifests and medians before S4b-S4d. Public operations provide the measured
  surface; no private colors or second color implementation are exposed. Setup-only results
  stay separately labelled because the old setup precomputes answers. Old incorrect answers
  are measured as old behavior, not accepted as expected results.
- **Kind:** additive (green). `[dep: S0a, S0b, S1d, S2b, S3c]`

S4a adds 67 exact cases to the existing form, span, search-color and canonicalization test
modules. They fix the full-color order, including field precedence and entity constraints;
distinguish raw literals, singleton sets and expressions; and cover the positional payloads that
are excluded from entity colors. The span cases cover all four variants, projection of equal
colors, and a modification confined to constraints. Canonical-result controls give explicit
expected molecule and span values after caller-owned normalization.

The following runs use `b98caac18` with only the S4a test and benchmark additions. Production
color, transport and symmetry implementations are unchanged. The compiler was rustc 1.99.1
(58cf8f9c3) on arm64. Criterion 0.8.1 used Nauty,
100 samples, 0.5 s warm-up and 1.5 s measurement time, matching S0b. All times below are medians in
microseconds, read from the named baseline's `estimates.json`.

```text
cargo bench --locked -p umol-graph-ir --bench canonicalize -- \
    --warm-up-time 0.5 --measurement-time 1.5 --save-baseline pre-246-colors \
    '^canonicalize/colors/'
cargo bench --locked -p umol-graph --bench transform -- \
    --warm-up-time 0.5 --measurement-time 1.5 --save-baseline pre-246-readings \
    '^stereo/(graph_symmetry|local_symmetry|classification|build_molecule_readings|build_stereo_readings)/'
```

`pre-246-colors` holds 24 cases. The six inputs are:

| Case | Input |
| --- | --- |
| entity_constraints | The eight-atom overlay input, with a constraint on each of the eight entity kinds, aromatic contributions 2,1,1,1 and multicenter contributions 2,1,1. |
| raw_forms | A copy of that input with singleton element, order, charge and coset forms, explicit zero charges, and a cis/trans coset term over literal 0. |
| aromatic_nonuniform | The five-carbon ring with contributions 2,1,1,1,1. |
| multicenter_nonuniform | Three carbon atoms with contributions 2,1,1. |
| dative_roles | Two boron atoms, one donor and one acceptor. |
| stereo_frame_constraints | The tetrahedral input with entity- and molecule-level topicity constraints from the existing corpus. |

The input manifest recorded counts in the order atom, localized bond, dative bond, aromatic
system, multicenter bond, noncovalent bond, stereo atom, stereo bond:

| Case | Molecule counts | Span counts |
| --- | --- | --- |
| entity_constraints | `[8, 8, 1, 1, 1, 1, 1, 1]` | `[9, 8, 1, 1, 1, 1, 1, 1]` |
| raw_forms | `[8, 8, 1, 1, 1, 1, 1, 1]` | `[9, 8, 1, 1, 1, 1, 1, 1]` |
| aromatic_nonuniform | `[5, 5, 0, 1, 0, 0, 0, 0]` | `[6, 5, 0, 1, 0, 0, 0, 0]` |
| multicenter_nonuniform | `[3, 0, 0, 0, 1, 0, 0, 0]` | `[4, 0, 0, 0, 1, 0, 0, 0]` |
| dative_roles | `[2, 0, 1, 0, 0, 0, 0, 0]` | `[3, 0, 1, 0, 0, 0, 0, 0]` |
| stereo_frame_constraints | `[5, 4, 0, 0, 0, 0, 1, 0]` | `[6, 4, 0, 0, 0, 0, 1, 0]` |

Each span adds one nitrogen atom to its molecule. Equality compares the input with a reversed
atom/entity numbering. Construction, renumbering and cloning for the consuming calls are outside
timing. The four columns are `canonicalize` and `canonical_eq` on the molecule and its span.

| Case | Molecule C | Molecule eq | Span C | Span eq |
| --- | ---: | ---: | ---: | ---: |
| entity_constraints | 19.361 | 38.567 | 22.057 | 41.911 |
| raw_forms | 19.961 | 39.442 | 22.478 | 42.864 |
| aromatic_nonuniform | 12.105 | 12.770 | 12.954 | 13.723 |
| multicenter_nonuniform | 7.327 | 6.233 | 7.722 | 6.652 |
| dative_roles | 4.664 | 4.431 | 4.829 | 4.540 |
| stereo_frame_constraints | 8.621 | 12.997 | 9.741 | 13.993 |

`pre-246-readings` holds 116 cases: the 84 existing setup, local-projection and class-reading
cases, plus four workloads on eight additional inputs. The new inputs are the approved meso and
chiral dichlorodifluoroethane fixtures, chlorofluoropropane, and two sampled meso
trichloropentane forms. Three more inputs copy the meso dichlorodifluoroethane fixture:
`readings_raw_term` stores `Term(Lit(0))` at its first center; `readings_raw_singleton` stores a
singleton carbon element and singleton coset there; `readings_constraints` stores singleton
valence 1 on H2 and singleton stereogenicity `Stereogenic` on the first stereo record.

The historical benchmark ID `readings_trichloropentane_chiral` inverts the middle record at
C3. S5a's exhaustive check identifies it as the second meso form. The approved chiral fixture
instead inverts a terminal record. The benchmark ID and input remain unchanged to preserve the
fixed workload and baseline; the new query cases use the correctly named chiral fixture.

The reading config is `ConstitutionColoring::full()`, fixpoint disabled, `max_iterations = 16`
and Nauty. These are the current controls, before the config replacement. The setup and local
columns below measure the existing `graph_symmetry` and `local_symmetry` workloads. The two
fixed workloads include setup once per iteration:

- `build_molecule_readings` asks chirality, both same-orbit predicates for atoms 0 and 1, and the
  complete proper and star orbit vector for every atom.
- `build_stereo_readings` asks the class, every unordered ligand pair's topicity, and membership
  of the proper ligand permutation [1,0,2,3] at each stereo record. Local projection is included.

| Case | Setup | Local | Build + molecule readings | Build + stereo readings |
| --- | ---: | ---: | ---: | ---: |
| readings_meso | 5.382 | 2.970 | 5.730 | 10.060 |
| readings_chiral | 6.489 | 3.018 | 7.008 | 11.242 |
| readings_prostereogenic | 7.074 | 4.586 | 6.920 | 12.957 |
| readings_trichloropentane_meso | 10.417 | 7.445 | 10.031 | 20.299 |
| readings_trichloropentane_chiral | 7.482 | 5.797 | 7.918 | 15.587 |
| readings_raw_term | 5.247 | 2.996 | 5.574 | 9.954 |
| readings_raw_singleton | 4.499 | 3.045 | 4.912 | 9.329 |
| readings_constraints | 5.444 | 3.000 | 5.732 | 10.120 |

The old molecule-query outputs were also recorded. Here `singletons` means that each atom
forms its own orbit; the other entries list complete orbit partitions. The raw-term,
raw-singleton and constrained variants all had the same eight-atom output as `readings_meso`.
These are observations of the old implementation; S5a establishes the expected replacement readings
independently. In particular, the old code incorrectly reported `readings_meso` as chiral.

| Case | Atoms | Old chiral | Old proper orbits | Old star orbits |
| --- | ---: | --- | --- | --- |
| readings_meso, readings_raw_term, readings_raw_singleton, readings_constraints | 8 | true | singletons | singletons |
| readings_chiral | 8 | true | `{0,1}, {2,5}, {3,6}, {4,7}` | `{0,1}, {2,5}, {3,6}, {4,7}` |
| readings_prostereogenic | 8 | true | singletons | singletons |
| readings_trichloropentane_meso | 11 | false | singletons | `{0,4}, {1,3}, {2}, {5,7}, {6}, {8,10}, {9}` |
| readings_trichloropentane_chiral | 11 | true | singletons | singletons |

The existing 84 cases retain their inputs and workload definitions:

| Case | Setup | Local projection | Class reading |
| --- | ---: | ---: | ---: |
| atoms_12 | 5.406 | 1.826 | 0.299 |
| atoms_132 | 76.063 | 27.544 | 0.299 |
| atoms_36 | 15.859 | 5.072 | 0.296 |
| atoms_5 | 2.992 | 1.190 | 0.299 |
| fixture_cis_trans | 3.239 | 1.248 | 0.180 |
| fixture_cis_trans_literal | 3.397 | 1.254 | 0.185 |
| fixture_evidenced_stereo | 3.950 | 2.875 | 0.483 |
| fixture_hydrogen_and_lone_pair | 2.389 | 1.003 | 0.295 |
| fixture_imine | 2.630 | 1.044 | 0.180 |
| fixture_implicit_hydrogen | 2.669 | 1.090 | 0.295 |
| fixture_lone_pair | 2.660 | 1.096 | 0.294 |
| fixture_nonstereogenic | 14.482 | 10.532 | 6.201 |
| fixture_partial_stereo | 4.085 | 3.010 | 0.675 |
| fixture_prostereogenic | 5.562 | 1.772 | 1.112 |
| fixture_remote_stereo | 7.042 | 5.845 | 0.892 |
| fixture_small_ring | 4.947 | 1.995 | 0.181 |
| fixture_square_planar | 3.161 | 1.172 | 0.286 |
| fixture_tetrahedral | 2.976 | 1.167 | 0.295 |
| fixture_tetrahedral_literal | 3.121 | 1.187 | 0.298 |
| ring_12 | 7.683 | 2.844 | 0.179 |
| ring_4 | 4.007 | 1.710 | 0.179 |
| ring_6 | 4.918 | 1.998 | 0.179 |
| ring_8 | 5.799 | 2.202 | 0.179 |
| sites_1 | 170.845 | 76.780 | 0.294 |
| sites_16 | 214.262 | 1411.280 | 4.718 |
| sites_2 | 171.040 | 137.126 | 0.594 |
| sites_4 | 172.686 | 286.103 | 1.180 |
| sites_8 | 183.406 | 596.208 | 2.351 |

These readings measure the current implementation, including its known incorrect results; they
are not correctness oracles for the replacement. Setup-only and projection-only costs do not
measure the same work as a fixed build-and-read request.

All 215 canonicalization and 51 active operation estimates in `pre-246` are unchanged, checked by
hash. The two retired-name records from S2a are also preserved. For a later `pre-246` comparison,
select the original 215 IDs with
`^canonicalize/(complete|inline_constraint|reaction|reaction_span|reaction_span_modified|scaling)/.*canonical_eq`.
The broader S0b `canonical_eq` filter now also selects 12 added color cases, whose baseline is
`pre-246-colors`.

The tables above record all 24 color-control and all 116 reading medians. Their input
constructors and workload definitions are committed in `umol-graph-ir/benches/canonicalize.rs`
and `umol-graph/benches/transform/stereo.rs` at `16f06bae5542`, the completed S4a revision.
Source hashes matched before and after both runs. The graph-IR package tests pass
(12,538 library tests, three ignored, plus integration and documentation tests). Both benchmark
targets pass the preflight checks, and clippy with `--all-targets -- -D warnings` passes for
`umol-graph-ir` and `umol-graph`. Formatting and the diff check pass.

#### S4b — `AsColor::as_inherent_color` **Done**

- **Module:** `ir::traits`, the eight entity-form modules and `ir::span`.
- **Tests first:** add exact cases for both methods on each form. `as_color` retains the current
  order and constraint distinctions; `as_inherent_color` omits only constraints and has a separate
  associated color type. Cover empty and nonempty stores, raw singleton/set expressions, each
  span variant, and a modified span whose sides differ only in excluded constraints. Check raw
  input retention and equal-color hash consistency.
- **Changes:** add the associated `InherentColor` type, the second trait method and their
  implementations. Inherent colors omit the constraints, including the stereogenicity assertion
  for stereo colors. Preserve the representation and declaration/comparison order of the existing
  full colors. Neither method normalizes a form. The trait and color types keep their existing
  crate visibility.
- **Kind:** additive (green). `[dep: S4a]`

The method has implementations for all eight forms and `EntitySpan`. Each inherent color
borrows the selected fields, or holds the stereo kind, without a constraints container. Full
colors retain their original representations and comparisons. The four whole-form wrappers are
single references, and the stereo and positional projections keep their full-color definitions.
S4e restores these wrappers and separates the inherent types after measuring the initial wider
representation.

The 42 new cases cover empty and nonempty stores, contradictory and raw constraints, raw scalar
and coset forms, all four span variants, and collapse of a modification confined to constraints.
Both color readings have equal-color hash checks. The cases assert that the supplied forms stay
structurally unchanged. Direct inherent-color assertions use the inherent projection of the
expected forms because the two readings have different Rust types. The exact shared-rank and
canonicalization expectations are unchanged. S4c enables the method for production use.

The initial graph-IR package gate passed: 12,580 library tests, three ignored, plus integration
and documentation tests. Clippy passed for all targets without warnings. Formatting and the
diff check passed.

Before S4e's representation correction, the 24 `pre-246-colors` cases were compared with the
same inputs, Nauty, 100 samples, 0.5 s warm-up and 1.5 s measurement time. The original baseline
estimates and benchmark source match their saved hashes. Median changes in the first comparison
are:

| Case | Molecule C | Molecule eq | Span C | Span eq |
| --- | ---: | ---: | ---: | ---: |
| entity_constraints | -1.53% | -1.77% | -1.67% | +0.49% |
| raw_forms | -2.59% | -1.70% | -1.50% | +0.66% |
| aromatic_nonuniform | -1.36% | -1.34% | -0.22% | +0.69% |
| multicenter_nonuniform | -1.56% | -1.05% | -0.25% | +1.02% |
| dative_roles | -1.07% | -0.50% | +0.13% | +2.13% |
| stereo_frame_constraints | -0.89% | -1.17% | -0.25% | +1.53% |

Criterion initially flagged the raw and dative span-equality cases. A focused rerun of those
cases and their molecule counterparts did not repeat the flags. The raw-span median was
42.777 µs against 42.864 µs; the dative-span median was 4.551 µs against 4.540 µs, within
Criterion's noise threshold. This check found no reproducible regression on these workloads.

```text
cargo bench --locked -p umol-graph-ir --bench canonicalize -- \
    --warm-up-time 0.5 --measurement-time 1.5 --baseline pre-246-colors \
    '^canonicalize/colors/'
cargo bench --locked -p umol-graph-ir --bench canonicalize -- \
    --warm-up-time 0.5 --measurement-time 1.5 --baseline pre-246-colors \
    '^canonicalize/colors/(molecule|reaction_span)/canonical_eq_equal/(raw_forms|dative_roles)$'
```

The baseline medians in S4a, the first-run changes in the table above and the two focused
recheck medians preserve the evidence used for this subitem's conclusion.

#### S4c — Shared search colors and canonicalization callers **Done**

- **Module:** `ir::symmetry`, `ir::molecule::symmetry`, `ir::reaction_span::symmetry`, and the two
  owning canonicalization modules.
- **Tests first:** extend `initial_colors` cases for both values of `include_constraints`, for
  all entity kinds and spans. Check incidence roles and per-participant electron counts, raw
  literals versus singleton sets, and constraint distinctions. Excluded constraints must not
  be normalized or checked for frame compatibility. Included constraints are read as stored.
  Keep the existing incompatible electron-count, kind, frame and literal-coset cases. Establish
  unchanged canonical representatives and comparison results before changing the call sites.
- **Changes:** add `include_constraints` to the shared color constructors. Select the full or
  inherent projection once and instantiate the same generic calculation for its associated color
  type. Keep the existing ranker unchanged. Enable `as_inherent_color` outside tests as these
  production callers land. Limit stereo compatibility checks to the data selected by the caller.
  Canonicalization passes `true` at every molecule and reaction-span call and retains its own
  upfront normalization. Its direct `as_color` reads, leaf ordering, normalized stereo checks
  and generator filters retain their behavior.
- **Kind:** breaking internally, red-to-green within the subitem as all call sites migrate.
  `[dep: S4b]`

The following describes the initial S4c implementation, superseded by S4f. Its selector
mechanism was not approved as part of the reference-sized color decision.

Both shared color constructors took `include_constraints`. A dispatcher selected the full or
inherent `ColorProjection` once; each aggregate's generic calculation then ranked the selected
associated color type with the existing ranker. The projection supplies a constant for the stereo
compatibility checks. Excluded stereo constraints are neither transported nor checked.
Configuration compatibility remains checked in both modes, on each present side of a span.
The electron-count, site-kind, frame and outer literal-coset checks remain. Neither constructor
normalizes the supplied forms.

Canonicalization changes only at its nine call sites, which pass `true`. Its normalization,
direct color reads, leaf ordering, normalized comparisons and generator filters are unchanged.
The trait method is enabled outside tests.

The 224 color cases cover both modes, all entity kinds, span variants, raw forms, incidence
roles, electron counts and incompatible selected data. Excluded incompatible constraints
succeed; incompatible configurations fail in both modes. Raw contradictory constraints are
read as stored when included. The existing canonicalization expectations are unchanged.

At `61f7918b6`, the all-features graph-IR test run passes: 12,732 library tests, three ignored,
plus 559 property tests, one ignored, and the integration and documentation tests. Clippy passes
with all features and targets and warnings denied. Formatting and the diff check pass. The
benchmark comparison remains the S4e gate, after S4d.

#### S4d — Shallow coset application **Done**

- **Module:** `ir::stereo` and its existing tests.
- **Tests first:** state exact results for literals, undetermined values, finite sets and terms,
  compatible and incompatible actions, and out-of-range outer literals. Preserve singleton/set
  variants and explicit variable domains. Cover identity actions on nested raw terms without
  entering their operators. Check ground inverse/composition laws and canonical value parity
  with normalized comparisons over all six kinds. For separately constructed symbolic actions,
  assert the exact resulting expression trees and state the value-level comparison separately.
- **Changed assertions:** show the old and new expectations before editing them. A term over a
  literal acted on by a nonidentity permutation stays an acted-on term instead of becoming a
  literal. An identity action on a term with an empty variable domain preserves that raw term
  instead of rejecting it through recursive normalization. Outer literal range and action
  compatibility failures remain failures.
- **Changes:** give `coset_apply_permutation` the shallow behavior in Computation. Remove its
  recursive `canon_coset` call on supplied terms. Recognize a new action that acts identically
  on the kind's cosets without rewriting its operand. Keep canonicalization's normalized
  comparisons and the existing public expression constructors.
- **Kind:** internal behavior change, green at the subitem boundary. `[dep: S4a, S3a]`

`coset_apply_permutation` checks action compatibility before dispatching on the outer form.
Literals and set members are transported directly; empty and singleton sets keep their variants.
An action that fixes every coset of the kind preserves a supplied term exactly, including its
operators and explicit variable domain. Other actions construct `Apply` around that whole term.
The helper no longer calls `canon_coset`. Incompatible actions fail for every outer variant,
including `Undetermined` and empty sets. Outer literal and set-member range checks remain.

The 71 application cases cover exact results, raw terms, identity actions on cosets, and action
and range failures. Across all six kinds, the algebra cases check inverses under every allowed
action and compositions with up to eight partner actions per action. Symbolic cases assert the
exact trees produced by separately applied actions and check their normalized value equivalence
independently. They also check exact operand preservation when the new actions are composed first.
Comparison controls use normalized inputs and retain `normalized_eq`, with and without a
composed mirror action.
The public expression constructors, frame transport and canonicalization sources are unchanged.

At `73698c3fe`, the all-features graph-IR test run passes: 12,791 library tests, three ignored,
plus 559 property tests, one ignored, and the integration and documentation tests. Clippy passes
with all features and targets and warnings denied. Formatting and the diff check pass. Reproduce
the S4c and S4d package gates with:

```text
cargo test --locked -p umol-graph-ir --all-features --quiet
cargo clippy --locked -p umol-graph-ir --all-features --all-targets -- -D warnings
```

S4e supplies the benchmark comparison and dependent-crate gate.

#### S4e — Shared-foundation gate **Done**

- **Module:** the S4b-S4d modules, their tests and the canonicalization benchmarks.
- **Check:** run the affected form, span, color and canonicalization tests and the explicitly
  enabled graph-IR property suite. Compare the canonicalization groups against `pre-246` and
  the added groups against `pre-246-colors`. Report medians and affected cases. Resolve measured
  regressions or settle the accepted cost before S5 begins. Run the dependent-crate stage gate
  after the internal signature migration. Corrections established by this gate are completed
  before the graph-symmetry consumer migration.
- **Kind:** verification and its shared-color correction (green). `[dep: S4c, S4d]`

The initial comparison at `73698c3fe` covered all 215 original canonicalization cases and all
24 added color cases. It used Nauty, 100 samples, a 0.5 s warm-up and a 1.5 s measurement time,
with the S0b and S4a baselines. These measurements precede the representation correction and
are not timings of the integrated generic version:

```text
cargo bench --locked -p umol-graph-ir --bench canonicalize -- \
    --warm-up-time 0.5 --measurement-time 1.5 --sample-size 100 --baseline pre-246 \
    '^canonicalize/(complete|inline_constraint|reaction|reaction_span|reaction_span_modified|scaling)/.*canonical_eq'
cargo bench --locked -p umol-graph-ir --bench canonicalize -- \
    --warm-up-time 0.5 --measurement-time 1.5 --sample-size 100 --baseline pre-246-colors \
    '^canonicalize/colors/'
```

The first original-case comparison produced 86 Criterion regression flags. A focused recheck
selected ten flagged cases above 1 µs, six constraint controls and two extreme count checks.
Four of the ten microsecond-scale flags repeated:

| Case | Baseline µs | First µs | Recheck µs | Recheck change |
| --- | ---: | ---: | ---: | ---: |
| scaling, equal, symmetry-heavy radicals | 7.5146 | 7.6649 | 7.6445 | +1.73% |
| molecule, equal, overlay-heavy | 35.8111 | 36.8263 | 36.1987 | +1.08% |
| molecule, unequal stereo, symmetric triisopropanolamine | 69.8670 | 71.1654 | 71.1130 | +1.78% |
| modified span, unequal colors, symmetric charge sum | 4.7768 | 5.0303 | 4.9947 | +4.56% |

The initially flagged count checks and newly flagged control varied between runs. The old
baseline comparison alone did not attribute the changes to this work. None of the 24 added
color-control cases was flagged in the initial S4e comparison; their median changes were:

| Input | Molecule C | Molecule eq | Span C | Span eq |
| --- | ---: | ---: | ---: | ---: |
| entity_constraints | -3.40% | -2.99% | -2.86% | -1.84% |
| raw_forms | -4.21% | -2.84% | -2.85% | -2.18% |
| aromatic_nonuniform | -2.84% | -2.64% | -1.30% | -1.97% |
| multicenter_nonuniform | -0.25% | -2.38% | -1.57% | -1.41% |
| dative_roles | +0.90% | -1.61% | -0.16% | -1.47% |
| stereo_frame_constraints | -2.40% | -2.59% | -0.81% | -1.13% |

A fresh comparison against `16f06bae5`, before S4b and with the identical benchmark harness,
localized an added cost to S4b's wider full-color values. The correction retains the original
full-color representations and derived comparisons, adds separate inherent-color types, and
selects the generic calculation once before the entity loops. The existing ranker and both
canonicalization modules are unchanged by this correction. A Boolean is not stored in a color.

The stage-isolation comparison used reference-S4b-S4c-S4c-S4b-reference order, with the same
100-sample settings. S4b was `2beacd9c9`; S4c was `61f7918b6`. Six cases were checked; the two
below drove the final focused comparison. The geometric means of the paired median estimates
give:

| Case | Reference µs | S4b µs | S4c µs | S4b/reference | S4c/S4b |
| --- | ---: | ---: | ---: | ---: | ---: |
| molecule, unequal colors, symmetric butan-2-ol | 8.66962 | 8.77744 | 8.82573 | +1.24% | +0.55% |
| modified span, unequal colors, symmetric charge sum | 4.84031 | 5.06407 | 5.04906 | +4.62% | -0.30% |

The full atom-color ranking entry had grown from 16 to 64 bytes; the corresponding span entry
had grown from 32 to 128 bytes. An AtomColor-only control added a Boolean to the form reference,
with the original ranker, callers and other color types unchanged. It measured +1.44% on the
molecule case and +0.18% on the span case in its focused run. Its atom entry was 24 bytes and
retained runtime flag checks. These observations motivated separate projection types rather
than carrying the reading choice in every value.

The four size controls were written and failed before the correction, then passed with the
single-reference atom, bond, dative-bond and noncovalent-bond wrappers. The scratch trial passed
467 focused color cases and 520 canonicalization cases; those filtered counts may overlap.
All eight full-color definitions and methods were checked against the reference implementation.

The accepted implementation is `93c1fa476` (`Add InherentColor`); the reference is
`16f06bae55423acab07c9d3fdefb8e5fab281ba9`. Both used the same committed
`umol-graph-ir/benches/canonicalize.rs` input definitions and workloads. The benchmark source's
Git blob is `3f7ae0fc4fab34998b007f68b7024b3300506823`. The compiler was rustc 1.99.1
(`58cf8f9c3`, 2026-10-13), on arm64, with Criterion 0.8.1 and Nauty.

The focused reference-candidate-candidate-reference comparison used 80 samples, a 0.2 s warm-up
and a 0.8 s measurement time per case. It took 11.8 s. The exact benchmark IDs were
`canonicalize/complete/canonical_eq_unequal_colors/symmetric_stereo_butanol` and
`canonicalize/reaction_span_modified/canonical_eq_unequal_colors/symmetric_charge_sum`.
The table uses the geometric mean of the two median estimates for each binary:

| Case | Reference µs | Generic µs | Change |
| --- | ---: | ---: | ---: |
| molecule, unequal colors, symmetric butan-2-ol | 8.61063 | 8.61913 | +0.10% |
| modified span, unequal colors, symmetric charge sum | 4.79817 | 4.82826 | +0.63% |

The individual median estimates, in run order and in microseconds, were:

| Case | Reference 1 | Generic 1 | Generic 2 | Reference 2 |
| --- | ---: | ---: | ---: | ---: |
| molecule | 8.58210 | 8.62674 | 8.61154 | 8.63925 |
| modified span | 4.80492 | 4.82534 | 4.83119 | 4.79143 |

The paired molecule changes were +0.52% and -0.32%; the paired span changes were +0.43% and
+0.83%. To repeat this comparison, build the benchmark once at each revision and save a copy of
each executable identified by Cargo's JSON `compiler-artifact` record before building the other
revision. Run those saved executables in the order above, using a distinct baseline label per
run. The executable and run-label arguments below stand for those paths and labels:

```text
cargo bench --locked -p umol-graph-ir --bench canonicalize --no-run --message-format=json
<benchmark-executable> --bench --warm-up-time 0.2 --measurement-time 0.8 --sample-size 80 \
    --save-baseline <run-label> \
    '^(?:canonicalize/complete/canonical_eq_unequal_colors/symmetric_stereo_butanol|canonicalize/reaction_span_modified/canonical_eq_unequal_colors/symmetric_charge_sum)$'
```

The focused comparison supported retaining reference-sized full colors and separate inherent
types. The selector mechanism itself was not approved; S4f below records its selected
replacement. These are preliminary results for the two measured workloads; they do not
establish performance parity across the complete suite.

The integrated Rust sources match the reviewed trial. The workspace test gate passes with
Python 3.13.18 active and `umol-graph-ir/proptest`, `umol-graph/proptest` and
`umol-graph/conformance` enabled: 32,287 passed and ten ignored. This includes the affected form,
span, color and canonicalization cases and the required property and conformance suites.
Workspace clippy passes with those features, all targets and warnings denied. The integrated
gate used:

```text
source umol-py/.venv/bin/activate && python --version
source umol-py/.venv/bin/activate && cargo test --locked --workspace \
    --features umol-graph-ir/proptest,umol-graph/proptest,umol-graph/conformance --quiet
source umol-py/.venv/bin/activate && cargo clippy --locked --workspace \
    --features umol-graph-ir/proptest,umol-graph/proptest,umol-graph/conformance \
    --all-targets -- -D warnings
cargo +nightly fmt --all -- --check
git diff --check
```

Formatting, the scope review and the diff check pass. These checks describe the selector
implementation; S4f below records its replacement.

#### S4f — Direct selection of full and inherent colors **Done**

The 2026-11-14 decision selects direct branching at each entity kind. The corrected trial is
integrated into the main worktree's existing color function bodies. The integration and
standalone stage-gate results are recorded below.
This subitem replaces the selector machinery recorded in S4c and S4e and precedes resumption of
unfinished S5b. It does not repeat S5a or settle the other S4/S5 review questions.

- **Module:** `ir::molecule::symmetry`, `ir::reaction_span::symmetry` and `ir::symmetry`.
- **Tests and measurements first:** reuse the existing S4a-S4c full/inherent color, span,
  encoding-failure and canonicalization cases with their exact expectations. Run the focused
  color and canonicalization tests on the integration base before changing production code.
  The trial and the three-version comparison below already precede main-worktree integration.
  This refactor requires no new correctness algorithm, fixture or changed expectation.
- **Changes:** retain the existing signatures of `molecule::symmetry::initial_colors` and
  `reaction_span::symmetry::initial_colors`, including `include_constraints`. At each entity-kind
  ranking site, branch on that argument, convert forms with `as_color` or `as_inherent_color`,
  and call `symmetry::rank_colors` with the corresponding associated color type. Both aggregate
  flags remain. Define each producer iterator and its encoding checks once; execute only the
  selected branch. For aromatic, multicenter and stereo entities, retain explicit vector
  capacity and push loops rather than introducing fallible collection. Preserve kind order,
  incidence colors and integer ranks.
- **Removal and boundaries:** remove `ColorProjection`, `FullColorProjection`,
  `InherentColorProjection`, `INCLUDE_CONSTRAINTS`, `molecule::symmetry::initial_colors_with`
  and `reaction_span::symmetry::initial_colors_with`. Introduce no replacement trait, marker,
  helper or color flag. Retain the two AsColor associated types and the reference-sized full
  colors. Preserve selected-data compatibility checks and failure behavior: excluded constraints
  are not transported or checked; included
  constraints are read as stored; configuration compatibility is checked in both modes. Add
  no normalization. Canonicalization continues to pass `true` and owns its normalization and
  normalized comparisons. Port only this selection change into the current function bodies;
  do not replace the containing files with older trial copies.
- **Gate:** rerun the focused tests, then the all-features graph-IR tests, including the
  explicit property suite, all-target Clippy with warnings denied, formatting and diff checks.
  Complete the existing S4 dependent-crate gate on a green integration base before resuming
  S5b. Review the full diff for extra abstractions, names, visibility changes and altered test
  coverage. Reuse the recorded timings when the port matches the reviewed trial; further edits
  to the measured path require a bounded comparison of the affected cases, not a tuning loop.
- **Kind:** internal refactor, green at completion; existing signatures are unchanged.
  `[dep: S4c, S4e]`

**Measurements before integration.** The three versions are original AsColor at `16f06bae5`,
the implemented selector at `93c1fa476`, and the corrected direct-branching trial based on
`93c1fa476`. Only the three source modules above differ in the measured trial. Its source-patch
SHA-256 is
`e31a888b02bbf8b5784d605aceded2db3eee157cbc22c65adffa7ae5edecae29`.
The original executable's hash matches its build record, and all 143 recorded source-file
hashes match its revision. The benchmark is unchanged, blob
`3f7ae0fc4fab34998b007f68b7024b3300506823`. All versions use rustc 1.99.1 (`58cf8f9c3`,
2026-10-13), Nauty and matching profiles: optimization level 3, debuginfo 2, debug assertions
and overflow checks disabled.

The exact workload IDs are
`canonicalize/complete/canonical_eq_unequal_colors/symmetric_stereo_butanol` and
`canonicalize/reaction_span_modified/canonical_eq_unequal_colors/symmetric_charge_sum`.
One original/selector/direct/direct/selector/original comparison used 80 samples, a 0.2 s
warm-up and a 0.8 s measurement time per case; the six runs took 17.80 s. The table gives the
geometric mean of the two median estimates per implementation, in microseconds. Negative
changes mean faster execution.

| Workload | Original AsColor µs | Selector µs | Direct µs | Direct vs original | Direct vs selector |
| --- | ---: | ---: | ---: | ---: | ---: |
| molecule, unequal colors, symmetric butan-2-ol | 8.66205 | 8.64902 | 8.50111 | −1.86% | −1.71% |
| modified span, unequal colors, symmetric charge sum | 4.84906 | 4.83831 | 4.82556 | −0.48% | −0.26% |

The individual medians, in run order and microseconds, were:

| Workload | Original 1 | Selector 1 | Direct 1 | Direct 2 | Selector 2 | Original 2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| molecule | 8.70996 | 8.69755 | 8.49195 | 8.51029 | 8.60077 | 8.61440 |
| modified span | 4.86579 | 4.84756 | 4.82362 | 4.82749 | 4.82908 | 4.83239 |

Paired direct-versus-original changes were −2.50% and −1.21% for the molecule, and −0.87% and
−0.10% for the span. Against the selector they were −2.36% and −1.05%, and −0.49% and −0.03%,
respectively. No penalty was observed on these two coloring-sensitive workloads. The small
span differences are treated as roughly equal performance; these results do not establish a
persistent speedup or describe the complete suite.

An earlier two-version trial also changed reserved vectors and push loops to fallible
collection. It measured +1.11% for the molecule and +0.72% for the span against the selector,
and omitted original AsColor. It does not establish the cost of direct branching. The corrected
comparison above preserves the allocation strategy and supplies the original baseline.

The corrected trial passes 467 focused color cases and 520 canonicalization cases; the filters
may overlap. All-features, all-target Clippy, formatting, diff checks and both workload
preflights for each executable pass. The trial's commands were:

```text
cargo test --locked -p umol-graph-ir --lib color --quiet
cargo test --locked -p umol-graph-ir --lib canonical --quiet
cargo clippy --locked -p umol-graph-ir --all-features --all-targets -- -D warnings
cargo +nightly fmt -p umol-graph-ir -- --check
git diff --check
cargo bench --locked -p umol-graph-ir --bench canonicalize --no-run --message-format=json
<saved-executable> --test <exact-workload-filter>
<saved-executable> --bench --warm-up-time 0.2 --measurement-time 0.8 --sample-size 80 \
    --save-baseline <unique-run-label> <exact-workload-filter>
```

The filter selects exactly the two IDs above. Separate Cargo and Criterion directories preserve
other builds and timing baselines. These inline measurements and source identities remain the
record when temporary trial directories are removed.

**Alternatives not implemented as the replacement.** A generic form-ranking helper could select
the reading and reuse one body per entity kind while retaining the two associated color types.
Its proposed iterator, lifetime, associated-type and fallibility requirements complicate the
signature to avoid repeated direct branches. It was not implemented or measured; no relative
performance claim is made. The repetition in the selected design is acceptable, and a later
abstraction remains possible if a concrete need justifies it.

A common output type for `as_color` and `as_inherent_color`, carrying the reading choice, remains
the fallback rather than the selected replacement. Earlier wider-color and flag experiments
motivated concern about approximately 1–4% costs on affected cases; that is not a measurement
of this unimplemented final design. This option trades simpler type machinery for a possible
persistent cost and would be reconsidered only if the selected direct design proved inadequate.

The selector trait and marker types were implemented at `93c1fa476`; they were not an
unimplemented alternative. They retained one generic aggregate body and chose the reading once,
but added an abstraction, names, visibility and placement decisions that were not approved.
S4f removes them; their measured timings appear in the comparison above.

**Integration checks.** The two integrated color function bodies are byte-for-byte identical
to the measured corrected trial. The port changes only these bodies, their imports and the
selector declarations in `ir::symmetry`. It preserves the unfinished S5b code around them,
the form/color types, the ranker, both canonicalization modules, all tests and the benchmark
harness. No additional timing run was needed.

Before the port, the main worktree passed 454 focused color cases and 520 canonicalization
cases. Both counts and results are unchanged after the port. The main all-features library
suite passes: 12,859 passed, three ignored. The lower color count than the trial reflects the
existing S5b changes; this port changes no tests or expectations.

The standalone S4f gate uses the green S4 base at `93c1fa476` with exactly the measured source
patch. With the repository's Python 3.13.18 environment active, the workspace tests pass:
32,287 passed, ten ignored. This includes every graph-IR feature and the explicitly enabled
property and conformance suites. Workspace Clippy with all targets and warnings denied,
workspace formatting and the diff check pass. The gate commands were:

```text
cargo test --locked --workspace \
    --features umol-graph-ir/proptest,umol-graph/proptest,umol-graph/conformance --quiet
cargo clippy --locked --workspace \
    --features umol-graph-ir/proptest,umol-graph/proptest,umol-graph/conformance \
    --all-targets -- -D warnings
cargo +nightly fmt --all -- --check
git diff --check
```

The following describes the S4f integration checkpoint. The main worktree remains in the
unfinished S5b migration. Its all-features package compile has
the same 22 error diagnostics before and after the port, from the retired coloring/config
surface and stereo methods in `tests/canonicalization.rs`, `tests/integrity.rs` and the stereo
property tests. Main Clippy reports the existing unused `StereoCoset` import and two absolute
paths in the S5b cache test. Main workspace formatting still reports the existing S5b edits;
the molecule module has exactly the same 27 formatter edits before and after S4f. These main
checks remain outstanding for the S5b/S5c migration. The standalone S4f gate above passes;
it does not assert that the main workspace migration is complete.

#### S4g — Field-color names and test corrections **Done**

- **Module:** ir::traits; the eight form/color owners; ir::span; molecule and reaction-span
  color consumers; their owning tests and the misplaced size controls in ir::symmetry::tests.
- **Tests first:** before the production rename, replace the nine method-derived expected
  values with explicit expected color fields and include cases distinguishing retained fields.
  Preserve the existing constraint-exclusion, equal-color hash and input-preservation checks.
  Move each of the four reference-size checks beside its color type, preserving the exact
  requirement. Use the existing exact-rank and canonicalization cases as correctness controls.
- **Changes:** keep AsColor and its full Color/as_color reading. Rename InherentColor to
  FieldColor, as_inherent_color to as_field_color, and all eight Inherent<Kind>Color types to
  <Kind>FieldColor. Update EntitySpan and both aggregate consumers in the same change. This
  changes names and test independence; the selected fields, layouts, comparisons and ranking
  algorithm are unchanged. Add no new abstraction or public surface.
- **Gate:** the 43 corrected color cases pass before the rename and as part of the all-features
  library suite afterward. The library, canonicalization/integrity integration suites, Clippy,
  benchmark linking and formatting pass. This correction ran after S5b within the ongoing S5
  migration; its workspace gate remains S5g. No new timing comparison was required.
- **Kind:** additive (green); internal names and test corrections only. `[dep: S4f]`

The completed S4b-S4f records retain their implementation-time inherent-color names. The
operative design above uses the approved field-color names; this subitem applies them to code.

### S5 — Symmetry data, molecule-backed queries and consumer migration

This is the breaking stage. S5a establishes evidence before the replacement. S5b retires the old
public surface; the workspace can remain red through S5f while its consumers migrate. The
`umol-graph-ir` library tests are green after S5b; its integration and property suites follow
after S5c supplies the stereo-view readings. `umol-graph` is green after S5e; S5f restores the
Rust/Python workspace. Each implementation subitem writes its required tests before its
production changes. S5g is the green stage boundary.

#### S5a — Exhaustive correctness checks and workload controls **Done**

- **Module:** `ir/molecule/symmetry/tests.rs`, the stereo property/exhaustive modules,
  `umol-graph/benches/transform/stereo.rs` and `umol-graph-ir/benches/canonicalize.rs`.
- **Tests and benches first:** bring in the needed approved S0a fixtures beside their tests.
  Establish expected remappings, chirality, complete orbit vectors, classes and published unit
  counts through exact cases and bounded independent enumeration. Include the dative donor and
  acceptor, the five-carbon aromatic system with electron counts 2, 1, 1, 1, 1, and a nonuniform
  multicenter case. The query assertions are written first in their owning implementation
  subitems; their expected values and fixture inputs are established here.
- **Benchmark controls:** check that S4a's `pre-246-readings` cases still measure the fixed
  build-and-read workloads before the API migration. Preserve their inputs and scope as their
  calls migrate. Do not replace that baseline after the shared transport change: the old
  symmetry implementation also calls `coset_apply_permutation`. Keep the 51 operation cases of
  `pre-246` and report setup-only and fixed-workload comparisons separately.
- **Kind:** additive (green). `[dep: S4e, S0a, S0b, S3b]`

All 38 approved S0a fixtures now live beside their tests in the permanent
[molecule symmetry tests](../umol-graph-ir/src/ir/molecule/symmetry/tests.rs) and
[exhaustive stereo checks](../umol-graph-ir/tests/property/stereo/exhaustive.rs).

The exhaustive symmetry calculation enumerates atom bijections on dense, simple constraint-free
graphs with literal binary configurations and literal incidence weights. It compares atom
and bond forms directly, preserves dative roles and weighted overlay participants, and checks
tetrahedral parity or cis/trans endpoint blocks without production colors, transport or
remapping search. It establishes 39 exact cases:

- 13 chirality and complete proper/star orbit vectors, including renumberings and virtual
  hydrogen ligands;
- 13 stereogenicity vectors on the approved multirecord fixtures;
- three complete remapping lists and orbit vectors for dative roles, aromatic counts
  `[2,1,1,1,1]` and multicenter counts `[2,1,1]`;
- four complete proper/improper remapping lists for the chiral and meso dichlorodifluoroethanes;
- six topicity cases covering all four variants.

An orbit vector records each atom's minimum numbered representative. These calculations establish
the expected values for the replacement queries, while existing production assertions and
property laws remain unchanged. Raw symbolic forms and constraint handling receive their query
assertions in S5b and S5c before those implementations change.

The literature checks exhaust all 420 binary assignments across 21 approved inputs, with
at most seven records and 128 assignments in one input. It groups complete canonical molecules
by structural equality using explicitly selected Nauty. Mirrored assignments establish the
chiral and achiral counts. Changing each record's bit over all assignments establishes the
unit set independently of the old stereo-symmetry readings. The separate published
Constitutional and Configurational sets remain explicit in the case rows for S5d; S5a validates
their union, rather than claiming to derive that subdivision.

The 13 Table I count controls are listed below; all 21 cases check complete unit sets. The
remaining eight inputs are Figure 1 structures a, c, f, g, i, j, k and n. Source figure numbers
are recorded beside every permanent fixture.

| Figure 7 no. | Input | Total | Chiral | Achiral |
| --- | --- | ---: | ---: | ---: |
| 1 | bromochloroethanol | 2 | 2 | 0 |
| 13 | bis(fluoroethyl)difluorohexane | 10 | 6 | 4 |
| 14 | dimethylcyclopentane | 3 | 2 | 1 |
| 16 | 1,2-dimethylcyclohexane | 3 | 2 | 1 |
| 17 | 1,3-dimethylcyclohexane | 3 | 2 | 1 |
| 18 | 1,4-dimethylcyclohexane | 2 | 0 | 2 |
| 19 | 1,2,3-trimethylcyclohexane | 4 | 2 | 2 |
| 22 | tetramethylcyclobutane | 4 | 0 | 4 |
| 28 | cyclohexanetriol | 2 | 0 | 2 |
| 30 | butene | 2 | 0 | 2 |
| 38 | bis(fluoroethyl)difluorohexene | 7 | 4 | 3 |
| 48 | decalin | 2 | 0 | 2 |
| 69 | inositol | 9 | 2 | 7 |

Four further cases enumerate every square-planar and trigonal-bipyramidal coset.
The platinum cis/trans fixtures have class vector `[0,1,0]`; the phosphorane axial/equatorial
fixtures have one class for cosets 0–7 and another for 8–19. Each case also checks the class of
its named stored configuration. Together these are 25 new literature/coset tests.

The three exhaustive calculations have continuing checking roles. Atom-remapping enumeration
supplies the independent calculation for molecule and per-record symmetry comparisons in S5b
and S5c. Assignment enumeration checks the canonicalizer's literature counts and supplies the
unit-membership comparison in S5d. Coset enumeration checks canonical configuration equivalence
and supplies the nonbinary stereo comparisons. Their private StereoRecord and StereoSite types
carry the data those calculations need. Names describe these computations and data, rather than
their role during development.

S5a's symmetry assertions check the exhaustive calculation against explicit expected values;
they do not establish correctness of the replacement GraphSymmetry queries. S5b and S5c add
those production comparisons before changing the query implementations. One-off derivation
drivers, diagnostic output and performance probes remain in scratch; none is imported into
these test modules. The expected cases and necessary evidence remain in permanent tests and
this record.

Both benchmark source files are byte-identical to S4a revision `16f06bae5542`. All 116 reading
workloads pass their Criterion preflight in test mode, including the unchanged historical
trichloropentane input described in S4a. All 408 saved baseline estimate files have unchanged
hashes before and after that check. This subitem records no new timing measurements and
replaces no baseline.

The all-features graph-IR gate passes with 13,472 tests and seven ignored, including 12,834
library tests and 584 property tests. The four stored nonbinary configuration assertions were
also checked after their addition. Clippy passes with all features, all targets and warnings
denied. The commands are:

```text
cargo test --locked -p umol-graph-ir --all-features --quiet
cargo test --locked -p umol-graph-ir --features proptest --test property test_exhaustive_coset_classes --quiet
cargo clippy --locked -p umol-graph-ir --all-features --all-targets -- -D warnings
cargo bench --locked -p umol-graph --bench transform --no-run --message-format=json
<transform-benchmark-executable> --test '^stereo/(graph_symmetry|local_symmetry|classification|build_molecule_readings|build_stereo_readings)/'
cargo +nightly fmt --all -- --check
git diff --check
```

Production behavior and public interfaces are unchanged. S5a ends green.

#### S5b — `GraphSymmetryConfig`, auxiliary data and molecule queries **Done**

- **Module:** `ir::molecule::symmetry`, `ir::view::symmetry`, `ir::view`, `ir::coloring` and
  the affected exports and graph-IR tests/benchmarks.
- **Readiness:** the single Constitution layout, flat auxiliary representation, field-color
  names, helper names and narrower generator filter are settled. S4g is independent. The remaining
  review corrections below follow the repository rules and do not need separate approval.
- **Tests first:** retain the trial's degree-mismatch cases as permanent cases before changing
  the main checks. Before replacing the coset-change generator filter, add the explicit mappings
  for a moving focus and a fixed square-planar focus whose other cosets are permuted, and the
  exact coset answers from the generator-filter comparison below. Reuse the existing raw-form,
  constraints, orbit and exact coset cases and the recorded comparison. Replace the provisional
  cache/pointer assertions with the query assertion; the selected layout has no second-search
  cache. Cover config defaults and both
  constraints settings; the five molecule queries and the same five view queries on S0a fixtures
  and their renumberings; owned and borrowed views; extraction, source/data independence and compatible reuse. Test the raw
  literal/singleton distinction, excluded contradictory or frame-incompatible constraints,
  included raw constraints, molecule-level constraints and encoding failures. Test the exact
  shallow symbolic operands and composed mirror/frame actions from the correctness trial.
  Check that mismatched context does not panic without promising a semantic result.
- **Changes:** complete the auxiliary-data implementation with one Constitution incidence
  graph, shared colors, adapter, refined partition and cell_indices. Inline incidence, adapter,
  partition and cell_indices into GraphSymmetry, remove SymmetrySearch, and construct the index
  with OrderedPartition::cell_indices, made pub(crate). Prepare the fields in
  Molecule::graph_symmetry and move the private cell method onto GraphSymmetry. Remove the
  provisional Structure input, its unused preparation branch and the lazy Constitution-search
  cache. The value holds
  no molecule, and config stays last. Preserve the nongeneric config with `automorphism_algorithm`, `max_assignments` and
  `include_constraints`. Build colors from supplied forms through S4f and check selected stereo
  encoding compatibility even though the graph omits stereo nodes. Complete the removal of
  `ir::coloring`, its three types, exports and feature-selector callers.
- **Interface:** `Molecule::graph_symmetry` returns `GraphSymmetryView`, owning the auxiliary
  data and borrowing the molecule. The view's owning constructor is crate-private;
  `from_graph_symmetry` borrows data and molecule, and `into_graph_symmetry` moves owned data
  or clones borrowed data. The auxiliary data supports that explicit cloning and has no
  arbitrary-parts constructor. The five query names move to `Molecule`, taking the data, and
  the view forwards them. Export the adapter with the other views and from `ir`.
- **Checks:** find image stereo records through their mapped sites in the supplied molecule.
  Keep the shared search setup as molecule::symmetry::molecule_remapping_exists and the shared
  stereo walk as molecule::symmetry::remapping_matches_stereo, with its existing Orientation
  and optional entity/coset arguments. Retain the private entity_image function and its
  signature. The generic ir::symmetry::remapping_exists keeps its name.
  Compare selected transported forms structurally, check degrees before composing mirror and
  frame actions, apply the composition shallowly, and skip excluded constraint transport.
  For coset-change queries, retain a preserving generator only if it fixes the focus record
  and preserves the requested target coset under its induced ligand-frame action. Use the
  selected direct implementation from the generator-filter comparison below, with no new
  production types or helpers. Proper and mirror queries retain their preserving filter.
  Preserve the normalized-comparison path for canonicalization, the generic search and the
  established cell exits. Every query uses the same prepared graph; the Constitution query
  omits the stereo predicate. Encoding failure uses distinct colors; lack of a normal form is
  not a reason for that fallback.
- **Migration:** retire `StereoSymmetry`, the two molecule methods returning it, the forwarding
  view methods and their conveniences with the old derivation. Remove their graph-IR tests
  only after recording which definitions the replacement cases state. Migrate the molecule
  query consumers in `tests/canonicalization.rs` and the config use in `tests/integrity.rs`;
  remove the isotope-ignoring orbit selection. Audit and build graph-IR benchmark targets,
  adapting actual affected callers. The stereo property consumers migrate in S5c, when their
  replacement readings exist.
- **Review corrections:** move GraphSymmetryView tests to its owning view module and keep
  tests in source order. Correct the mixed stereo-view test so its name, values and assertions
  describe the current operations. Use symmetry for GraphSymmetry locals, config in its
  rustdoc, and an accurate module description of its colors and checks. Use the retirement
  ledger below to account for every
  removed assertion, recording the actual test/code ordering. Finish the graph-IR integration
  consumers, formatting and Clippy work; build affected benchmark targets before claiming
  benchmark migration is needed. Do not copy temporary counters or drivers into main.
- **Kind:** breaking; graph-IR library tests green, integration/property consumers can remain
  red until S5c and the workspace until S5f.
  `[dep: S5a, S4f, S4d]`

**Completion checkpoint, 2026-11-14.** S5b was completed first, followed by the independent
S4g correction. Before the search port, the focused suite had 259 passing and six failing cases:
two kind/frame-degree panics and four Constitution-remapping checks. After the port all 265
passed. Existing late-written tests retain their earlier chronology; this does not claim that
they were test-first.

The completed code has one Constitution graph, flat GraphSymmetry fields with config last,
the existing cell-index calculation, the approved private helper names, the degree guard and
the selected narrower generator filter. Adapter tests now live in the view module. The mixed
stereo-view test is separated into kind chirality and exact ligand-position cases; its irrelevant
own-coset variation is removed. The isotope-ignoring corpus selector is retired as planned.

S4g gives the eight atomic color tests explicit expected fields or kinds. The span test constructs
its expected variants explicitly, using the separately tested atomic color producer for its
opaque atomic values; the expected carrier never calls the span projection or collapse. All
43 field-color cases passed before the production rename. The four reference-size assertions
now live beside their types. Full and field-color implementations retain their fields, layout,
comparison and hash behavior; only the agreed names change.

Final verification in main passes: 12,874 all-features library cases, three ignored; nine
canonicalization and 19 integrity integration cases; Clippy with warnings denied for those
targets and all benchmarks; all ten benchmark executables linked in the test profile; workspace
formatting and the diff check. The commands were:

```text
cargo test --locked -p umol-graph-ir --all-features --lib --test canonicalization --test integrity
cargo clippy --locked -p umol-graph-ir --all-features --lib --test canonicalization --test integrity --benches -- -D warnings
cargo test --locked -p umol-graph-ir --all-features --benches --no-run
cargo +nightly fmt --all -- --check
git diff --check
```

The scope audit confirms that both canonicalization modules are unchanged, the generic search
changes only the approved cell_indices visibility, and both existing canonical stereo checks
and full-color bodies are unchanged apart from formatting. Field-color bodies change only
names and formatting. No temporary benchmark driver, counter, selector trait or extra production
helper was retained. The unrelated umol-graph stereo edit is preserved. No new timing run or
mutating git command was used. The stereo property consumers remain assigned to S5c, and the
higher Rust/Python consumers to S5e/S5f; this is not a workspace-wide completion claim.

**Retired and changed test requirements.** This ledger covers the current S5b working-tree
changes. A pending replacement remains work for the named subitem.

| Test or group | Requirement represented | Current disposition and replacement |
| --- | --- | --- |
| test_graph_symmetry_site_stabilizer, test_molecule_observable_descriptor, test_project_onto_ligands | Details of the removed stabilizer/descriptor/local-projection implementation. | Those internals have no requirement in the replacement. Their helper-specific assertions are retired. |
| test_project_onto_ligands_virtual | Implicit-hydrogen, lone-pair and virtual ligands keep their positions under a remapping. | The rule remains. Implicit hydrogen is covered by the remapping-exists and S3a cases; lone-pair and virtual-ligand production coverage remains required in S5c. |
| test_molecule_stereo_atom_symmetry_stereogenic, test_molecule_stereo_atom_symmetry_prostereogenic, test_molecule_stereo_bond_symmetry | Classes of the elementary atom/bond examples and topicities of selected ligand pairs. | Preserve the class requirements in S5c. Pairs of constitutionally different ligands now read Distinct instead of Diastereotopic under the settled definition; show the changed exact values in those cases. |
| test_stereo_symmetry_is_stereogenic, 15 cases | The removed local-group rule and independence from the own coset for binary kinds. | The local-group rule is replaced by the remapping definition. Binary own-coset independence remains required in S5c. The existing multicoset remapping cases are partial coverage, not a replacement for the derived reading. |
| test_stereo_atom_view_stereo_queries, now test_stereo_atom_view_is_chiral | Forwarded stereo readings, plus kind chirality and ligand lookup over four own-coset forms. | S5b separates kind chirality and ligand-position tests and retains the existing exact ligand/count checks. The irrelevant own-coset variation is removed. Restore derived-reading coverage with attached symmetry data in S5c. |
| test_graph_symmetry_is_chiral_no_stereo, now test_graph_symmetry_view_is_chiral | Benzene is achiral. | The query assertion remains; place it with the owning view tests. |
| test_exhaustive_symmetry, now test_molecule_graph_symmetry | Explicit chirality and full orbit expectations established by exhaustive atom-bijection calculations. | Kept and extended to compare the molecule queries and owned/borrowed views with the same exact values. |
| rdkit_isotope_excluded in tests/canonicalization.rs | The old feature selector could ignore isotope fields. | Feature selection is retired. The isotope-including exact partition and renumbering case remain; there is no isotope-ignoring query in the new interface. |

**Search-layout comparison and decision, 2026-11-14.** The single Constitution search is
selected for S5b. Both designs were implemented in a separate worktree, based on the current
S5b source over `a02f535b6`. Main still contains the provisional two-search implementation;
porting the selected layout and its correctness fixes is pending. The search-bundle and field
names were accepted separately, and config ordering is already applied. The temporary drivers
and the rejected candidate are not part of the production implementation.

The two-search candidate prepares Structure data and lazily prepares Constitution data. The
one-search candidate prepares Constitution data once and uses it for all queries; stereo images
are identified through their mapped atom or bond sites. Both retain SymmetrySearch for this
comparison, use the same ranker and search kernel, read the supplied molecule, and preserve raw
forms without normalization. Both check selected stereo encoding compatibility when preparing
Constitution data, so omitting stereo nodes does not bypass the encoding-failure fallback.

Before timing, a kind/frame-degree mismatch reproduced the mirror-composition panic with both
constraint settings. Both candidates check the degree before composing. Coset-change searches
clear the backend generators to disable generator pruning for that predicate; ordinary preserving
and mirror searches retain their filter. Bypassing the filter would retain all backend generators,
which is not the same operation. The existing exact coset cases pass with pruning disabled.
This describes the candidates used for the layout timings. The later generator-filter
comparison below selects a narrower preserving filter for integration; it does not change
these recorded layout results.
The provisional cache/pointer assertions were removed from the trial's Constitution-query test;
its query assertion remains. These correctness changes are in the trial, pending integration
with the selected layout.

Each candidate passes the same 254 molecule-symmetry test cases, including raw forms, selected
and excluded constraints, malformed encodings, the new degree cases and the existing coset
cases. The eight outputs for the four fixed fixtures under both constraint settings agree
exactly. This is a focused correctness check, not the completed S5b workspace gate.

The fixtures are the accepted meso and chiral dichlorodifluoroethanes and chlorofluoropropane
from the permanent molecule-symmetry tests, and the inositol input from the stereo operation
benchmarks: six `C#h` ring atoms, six `O#h` substituents, and six tetrahedral records with
undetermined cosets. Timings use the default config, Nauty and excluded constraints. The fixed
workloads, established before adapting the candidates, are:

| Workload ID component | Timed operation |
| --- | --- |
| prepare | Build and drop the owning symmetry view. |
| prepared_readings | On prepared data, read chirality, proper/star equivalence of atoms 0 and 1, and both complete orbits of atom 0. |
| first_constitution | Query constitutional equivalence with the site fixed, on fresh prepared data. Preparation and input destruction are outside timing through Criterion PerIteration batching by reference. |
| build_and_read | Prepare data, run the five prepared readings, then the Constitution query; includes destruction. |
| repeated_readings | After the Constitution query has initialized any cache, repeat the five readings and Constitution query three times. |

The Constitution query fixes atom 0 and sends atom 1 to atom 2 for both dichlorodifluoroethanes;
fixes atom 1 and sends atom 6 to atom 7 for chlorofluoropropane; and fixes atom 0 and sends atom 1
to atom 6 for inositol. It measures the primitive needed by topicity's constitutional fallback,
not the as-yet-unimplemented complete topicity method. All but the chlorofluoropropane query
exit by different cells; the latter executes one candidate check. The repeated workload is a
fixed comparison, with no claim about real query frequencies.

The run used rustc 1.99.1 (`58cf8f9c3`, 2026-10-13), matching optimized benchmark profiles
(optimization level 3, debuginfo 2, debug assertions and overflow checks disabled), 80 samples,
0.1 s warm-up and 0.4 s measurement per workload. The fixed order was two/one/one/two, with
separate Criterion output directories and no concurrent builds. The four runs took 12.75,
12.61, 12.66 and 12.70 s: 50.72 s total. Allocation and check-count instrumentation is enabled
only in the debug binaries; it is absent from the timed binaries.

Every full benchmark ID is `layout/<workload>/<fixture>`, with the five workload components
above and fixtures `meso`, `chiral`, `chlorofluoropropane` and `inositol`. The table gives the
geometric mean of each candidate's two median estimates, in microseconds. Negative changes
mean that one search is faster; paired changes compare corresponding repeats.

| Workload | Fixture | Two µs | One µs | One vs two | Pair 1 | Pair 2 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| build_and_read | chiral | 45.27982 | 23.91083 | -47.19% | -47.90% | -46.48% |
| build_and_read | chlorofluoropropane | 26.79639 | 16.87100 | -37.04% | -37.97% | -36.09% |
| build_and_read | inositol | 368.93548 | 148.07002 | -59.87% | -60.06% | -59.67% |
| build_and_read | meso | 60.57229 | 32.86440 | -45.74% | -46.56% | -44.92% |
| first_constitution | chiral | 2.58497 | 0.02180 | -99.16% | -99.16% | -99.16% |
| first_constitution | chlorofluoropropane | 6.85427 | 3.52282 | -48.60% | -49.56% | -47.63% |
| first_constitution | inositol | 3.58724 | 0.02188 | -99.39% | -99.39% | -99.39% |
| first_constitution | meso | 2.60629 | 0.02181 | -99.16% | -99.17% | -99.15% |
| prepare | chiral | 4.05133 | 2.62747 | -35.15% | -35.84% | -34.44% |
| prepare | chlorofluoropropane | 5.31540 | 3.37707 | -36.47% | -37.49% | -35.42% |
| prepare | inositol | 6.77280 | 3.55159 | -47.56% | -48.19% | -46.93% |
| prepare | meso | 4.06545 | 2.65204 | -34.77% | -35.47% | -34.05% |
| prepared_readings | chiral | 37.93832 | 21.18839 | -44.15% | -44.62% | -43.68% |
| prepared_readings | chlorofluoropropane | 14.49105 | 9.77477 | -32.55% | -33.51% | -31.57% |
| prepared_readings | inositol | 356.14304 | 145.00263 | -59.29% | -59.99% | -58.57% |
| prepared_readings | meso | 53.54934 | 29.90544 | -44.15% | -45.22% | -43.07% |
| repeated_readings | chiral | 113.76004 | 63.48991 | -44.19% | -44.84% | -43.54% |
| repeated_readings | chlorofluoropropane | 53.85235 | 40.15540 | -25.43% | -26.76% | -24.09% |
| repeated_readings | inositol | 1069.08068 | 433.81367 | -59.42% | -59.93% | -58.90% |
| repeated_readings | meso | 160.83263 | 89.53043 | -44.33% | -45.06% | -43.59% |

The individual median estimates, in run order and microseconds, are:

| Workload | Fixture | Two 1 | One 1 | One 2 | Two 2 |
| --- | --- | ---: | ---: | ---: | ---: |
| build_and_read | chiral | 45.44777 | 23.67895 | 24.14499 | 45.11250 |
| build_and_read | chlorofluoropropane | 26.97102 | 16.72947 | 17.01372 | 26.62289 |
| build_and_read | inositol | 370.15325 | 147.82668 | 148.31377 | 367.72171 |
| build_and_read | meso | 61.01548 | 32.60815 | 33.12266 | 60.13233 |
| first_constitution | chiral | 2.59685 | 0.02191 | 0.02169 | 2.57314 |
| first_constitution | chlorofluoropropane | 6.91326 | 3.48680 | 3.55922 | 6.79578 |
| first_constitution | inositol | 3.60744 | 0.02190 | 0.02185 | 3.56716 |
| first_constitution | meso | 2.62760 | 0.02177 | 0.02185 | 2.58515 |
| prepare | chiral | 4.07278 | 2.61296 | 2.64206 | 4.03000 |
| prepare | chlorofluoropropane | 5.34543 | 3.34124 | 3.41329 | 5.28553 |
| prepare | inositol | 6.82578 | 3.53674 | 3.56650 | 6.72023 |
| prepare | meso | 4.07928 | 2.63233 | 2.67189 | 4.05168 |
| prepared_readings | chiral | 38.15405 | 21.13138 | 21.24556 | 37.72382 |
| prepared_readings | chlorofluoropropane | 14.58917 | 9.70032 | 9.84979 | 14.39359 |
| prepared_readings | inositol | 360.72768 | 144.32376 | 145.68470 | 351.61667 |
| prepared_readings | meso | 54.13189 | 29.65351 | 30.15951 | 52.97305 |
| repeated_readings | chiral | 114.30966 | 63.05780 | 63.92499 | 113.21306 |
| repeated_readings | chlorofluoropropane | 54.35749 | 39.81334 | 40.50039 | 53.35191 |
| repeated_readings | inositol | 1080.94580 | 433.10417 | 434.52433 | 1057.34580 |
| repeated_readings | meso | 162.05672 | 89.02683 | 90.03689 | 159.61779 |

The untimed debug-build counts explain the work and storage difference. Storage below is the inline value
plus live requested Rust heap bytes retained by it; it excludes allocator overhead and native
C allocations. Allocation counts cover Rust allocation/reallocation calls during preparation.
The prepared-read column counts final candidate predicates / stereo-check calls, including
stereo checks used by the generator filter. These counts are identical between candidates.

| Fixture | Backend vertices two/one | Initial cells two/one | Prepared bytes two/one | After fallback bytes two/one | Preparation allocations two/one | Prepared checks, both |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| meso | 27/15 | 14/8 | 3368/1838 | 5022/1838 | 60/47 | 9/10 |
| chiral | 27/15 | 14/8 | 3368/1838 | 5022/1838 | 60/47 | 6/7 |
| chlorofluoropropane | 27/15 | 24/13 | 3836/2178 | 5830/2178 | 80/58 | 3/4 |
| inositol | 60/24 | 9/4 | 6126/2504 | 8446/2504 | 68/48 | 13/13 |

The inline GraphSymmetry values are 368 bytes for two searches and 184 for one. The initial
refined cell sizes are:

| Fixture | Two searches | One search |
| --- | --- | --- |
| meso | `[2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2]` | `[2, 2, 2, 2, 2, 1, 2, 2]` |
| chiral | `[2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2]` | `[2, 2, 2, 2, 2, 1, 2, 2]` |
| chlorofluoropropane | `[2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1]` | `[2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1]` |
| inositol | `[6, 6, 6, 6, 6, 6, 12, 6, 6]` | `[6, 6, 6, 6]` |

On these fixtures, one search examines the same number of final candidates and stereo checks.
It prepares and queries smaller backend graphs, clones smaller partitions, and retains one
input instead of two after the fallback. Identifying stereo records through mapped sites adds
work at the leaves, but the measured totals favor the smaller input. This is not a claim that
the number of candidates stays equal for every molecule: the Structure graph can provide
additional pruning when stereo kinds or selected constraints distinguish otherwise equivalent
sites.

Preparation is 35–48% faster, prepared readings 33–59% faster, build-and-read 37–60% faster,
and repeated readings 25–59% faster with one search. Both paired comparisons agree in direction
for every workload. The approximately 22 ns first-Constitution results for three fixtures are
cell exits after preparation, not complete topicity timings. The one-search design also retains
43–59% less Rust storage initially and 63–70% less after the fallback on this subset.
The single Constitution search is selected on this evidence. Its implementation in main and
the remaining S5b corrections are still pending. The comparison supplies the layout decision;
no additional timing or tuning run is a prerequisite for the port.

The source identities below preserve the evidence independently of temporary directories.
The main S5b source snapshot before the ordering change was:
`9575d751646db8197fda700e58be1354f017364bd921ebdf635cb56f2659ff46`.

| Candidate | Molecule-symmetry source SHA-256 | Optimized executable SHA-256 |
| --- | --- | --- |
| two | `635678775ce628b689c1abdf0813fd1c77f2e5d72d781b5987495f800938a2ab` | `24485cc07a3d81bb181f7e386ccdad676909e45875ce30677361a50a9b5d2878` |
| one | `eef5912e1c408730c7bd143c4cead2df019d5110c64d0b50de56840c22b36d42` | `4ac5b999a4a752c4b2384e3fb46c8484c8e4053298029113a8b31ca5d77fefa9` |

The identical owning test module has SHA-256
`db7a994d99ac956e1aa3636b6df5a7b6b2d225560952975ec53f1ce6e790f725`.
The temporary workload-driver hashes are `390bed9c89a4c4a9bec8e4534c34afee6372f1001f9864a4bb06073a143e6bb7`
for two searches and `31902a5581a47d0add66890793060aa18e2f02e2f6cb58376ec0346874bb3a80`
for one. They differ only in the field used to report the initial partition.

Build each fixed source version and save the executable identified by Cargo's JSON artifact
before building the other version. The commands are:

```text
cargo test --locked -p umol-graph-ir --lib --no-run --message-format=json
<debug-executable> ir::molecule::symmetry::tests --test-threads=1
UMOL_LAYOUT_MODE=counts <debug-executable> \
    ir::molecule::symmetry::tests::test_molecule_graph_symmetry_layout \
    --exact --nocapture --test-threads=1
cargo test --locked -p umol-graph-ir --lib --no-run --message-format=json --profile bench
UMOL_LAYOUT_MODE=time UMOL_LAYOUT_CRITERION=<distinct-output-directory> <saved-executable> \
    ir::molecule::symmetry::tests::test_molecule_graph_symmetry_layout \
    --exact --nocapture --test-threads=1
```

The candidate sources and temporary driver remain in the separate checkout for review. The
main code changes in this comparison are limited to placing config last, both as a field and
as the module's final type. The main graph-IR library compiles with its existing unused-import
warning; the unfinished S5b migration and its previously recorded gate failures remain open.

**Generator-filter comparison and decision, 2026-11-14.** Retaining only preserving generators
that also fix the focus record and requested target coset is selected for the coset-change
search. Clearing all generators is the alternative. They have the same query semantics and
differ in pruning cost and implementation complexity. Both still compute the backend
generators; clearing discards them before they can prune the search.

This is a defensible accepted tradeoff, not an established best choice. Trading a gain in a
particular case against slowdowns in common cases is ordinarily disfavored. The reasons for
accepting it here are the size of the roughly tenfold gain and the expectation that the
slowdown will be smaller within the entire operation, of which this filter is only one part.
The latter is an expectation, not an end-to-end measurement. The selected inputs do not
establish workload frequencies or an overall stereogenicity speedup.

The filter question arises in molecule::symmetry::remapping_exists_for, currently test-only,
which S5c will use when asking whether another coset changes a stereo entity's description.
The search requires its final check to be invariant under every retained target generator.
Preserving the stored molecule does not suffice when the target changes one record's coset.
For example, exchanging two equivalent records with coset 0 preserves the stored molecule but
exchanges which record has the requested coset 1. For a square-planar record, a permutation
can fix coset 0 while exchanging cosets 1 and 2. In either case a generator can relate a
passing candidate to a failing candidate, making its use for pruning unjustified.

The selected filter first identifies the focus record's kind, ligands and atom or bond site.
It retains a generator only if the site maps to itself, the stored stereo forms are preserved,
and the induced permutation of the focus's ligands carries the requested literal coset to
itself. The unique record at that site is then fixed. An unavailable record or kind clears
the list; an unavailable image or frame rejects the generator. Coset transport uses the
existing shallow operation. No normalization or change to the query's final predicate is
introduced.

**Correctness before the filter change.** The new fixture assertions and timing driver were
installed before the narrower filter. The initial focused run passed 33 generator-related
cases, including the following explicit mappings. Each is checked with constraints included
and excluded. The generator preserves the stored forms; the candidate reaches target coset
1; composing that generator after the candidate fails the same coset-change check.

- For the two-center input below, the candidate swaps atoms 4 and 5. The generator exchanges
  the two ends with atom images `[1, 0, 3, 2, 6, 7, 4, 5]`; incident bonds follow those atom
  images. It moves the record carrying the change.
- For `{:atoms ["Pt#c-2" "F" "F" "F" "F"] :bonds [[0 1 "1"] [0 2 "1"]
  [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]}`,
  the candidate has atom images `[0, 1, 3, 2, 4]` and the generator `[0, 3, 2, 1, 4]`.
  The generator fixes the record and coset 0 but moves the candidate's target coset 1.

These demonstrate failure of the old filter's invariance requirement. They are not a claim
that the current backend has produced a wrong public answer. Each safe candidate then passes
all 261 molecule-symmetry cases. The five inputs below supply 30 explicit coset answers, checked
with both constraint settings for 60 expected-answer checks per candidate, together with input
preservation. Temporary drivers and counters are excluded from the main implementation.

**Fixed inputs.** Every focus is stereo atom 0, every stored coset is 0, and no atoms are held
fixed. The complete inputs are retained here so the evidence survives removal of the trial
directories.

Tetrahedral, with expected answers `[true, true]` for cosets 0 and 1:

```edn
{:atoms ["C#h" "F" "Cl" "Cl"]
 :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
 :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th0"}]}
```

Two centers, with expected answers `[true, true]`:

```edn
{:atoms ["C" "C" "Cl" "Cl" "H" "H" "H" "H"]
 :bonds [[0 1 "1"] [0 2 "1"] [1 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
 :stereo-atoms [{:site 0 :ligands [1 2 4 5] :attrs "Th0"}
                {:site 1 :ligands [0 3 6 7] :attrs "Th0"}]}
```

Square planar, with expected answers `[true, false, true]`:

```edn
{:atoms ["Pt#c-2" "N#c+#h3" "N#c+#h3" "Cl" "Cl"]
 :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
 :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]}
```

Trigonal bipyramidal, with true for cosets 0–7 and false for cosets 8–19:

```edn
{:atoms ["P" "Cl" "F" "F" "F" "F"]
 :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"]]
 :stereo-atoms [{:site 0 :ligands [1 2 3 4 5] :attrs "Tb0"}]}
```

Square planar with explicit remote ammonia hydrogens, with expected answers
`[true, false, true]`:

```edn
{:atoms ["Pt#c-2" "N#c+" "N#c+" "Cl" "Cl" "H" "H" "H" "H" "H" "H"]
 :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]
         [1 5 "1"] [1 6 "1"] [1 7 "1"] [2 8 "1"] [2 9 "1"] [2 10 "1"]]
 :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]}
```

**Timing boundary and method.** Each iteration makes one complete call to
molecule::symmetry::remapping_exists_for for one focus and one requested coset. It includes
partition cloning and refinement, backend automorphism computation, generator filtering,
candidate search, stereo comparisons, and temporary allocation and deallocation. The molecule
and GraphSymmetry are prepared once outside the timed loop and reused. Their construction,
incidence graphs, colors and initial refined partition are excluded. This does not time the
full stereogenicity calculation, which can take early exits or try several cosets.

The eight workloads were fixed before implementing the narrower filter. Each uses Nauty,
include_constraints=false, 40 Criterion samples, 100 ms warm-up and 400 ms measurement.
One clear/narrow/narrow/clear sequence was run, taking 5.27, 5.04, 4.97 and 5.20 seconds for
the respective eight-workload timing runs. No further timing or tuning iterations followed.
The builds use rustc 1.99.1 (58cf8f9c3, 2026-10-13), aarch64-apple-darwin, LLVM 23.1.1,
opt-level 3, debuginfo 2, and disabled debug assertions and overflow checks. Debug counters and
the allocation instrumentation are compiled out of the timed binaries.

Values below are geometric means of the two run medians, in microseconds per query. A positive
change means the narrower filter is slower. Every paired comparison agrees on the direction.

| Input / target coset | Clear | Narrow | Narrow change | Paired changes |
| --- | ---: | ---: | ---: | --- |
| Tetrahedral / 1 | 5.64877 | 5.94425 | +5.23% | +5.27%, +5.19% |
| Two centers / 1 | 35.71322 | 38.25153 | +7.11% | +7.08%, +7.13% |
| Square planar / 1 | 15.16459 | 16.20217 | +6.84% | +6.68%, +7.00% |
| Square planar / 2 | 15.20295 | 16.23231 | +6.77% | +6.56%, +6.98% |
| Trigonal bipyramidal / 7 | 88.96847 | 94.24794 | +5.93% | +6.98%, +4.90% |
| Trigonal bipyramidal / 8 | 88.70805 | 94.39726 | +6.41% | +7.46%, +5.38% |
| Remote ammonia / 1 | 799.05416 | 78.12945 | -90.22% | -90.05%, -90.39% |
| Remote ammonia / 2 | 799.96942 | 77.76531 | -90.28% | -90.17%, -90.38% |

The raw run medians, also in microseconds, are:

| Input / target coset | Clear 1 | Narrow 1 | Narrow 2 | Clear 2 |
| --- | ---: | ---: | ---: | ---: |
| Tetrahedral / 1 | 5.64005 | 5.93749 | 5.95103 | 5.65751 |
| Two centers / 1 | 35.74306 | 38.27429 | 38.22879 | 35.68341 |
| Square planar / 1 | 15.20078 | 16.21616 | 16.18818 | 15.12848 |
| Square planar / 2 | 15.22287 | 16.22192 | 16.24271 | 15.18306 |
| Trigonal bipyramidal / 7 | 88.77812 | 94.97465 | 93.52678 | 89.15922 |
| Trigonal bipyramidal / 8 | 88.90993 | 95.54127 | 93.26696 | 88.50664 |
| Remote ammonia / 1 | 798.00450 | 79.37326 | 76.90512 | 800.10521 |
| Remote ammonia / 2 | 800.91522 | 78.70946 | 76.83248 | 799.02473 |

**Search work and complexity.** These counts come from separate debug runs. Allocation counts
are Rust allocation/reallocation calls during the query, excluding preparation and native C
allocation. Generator totals accumulate across filter invocations rather than counting
independent generators of one group. Stereo-check counts include the original preservation
checks made by the filter; the additional direct focus-frame checks are reflected in allocations
and time.

| Input / target coset | Clear/narrow candidates | Clear/narrow stereo checks | Clear/narrow allocations | Narrow generators seen/kept |
| --- | ---: | ---: | ---: | ---: |
| Tetrahedral / 1 | 3/3 | 3/4 | 137/139 | 1/0 |
| Two centers / 1 | 9/9 | 9/15 | 640/658 | 7/0 |
| Square planar / 1 | 5/5 | 5/9 | 316/324 | 4/0 |
| Square planar / 2 | 5/5 | 5/9 | 316/324 | 4/0 |
| Trigonal bipyramidal / 7 | 25/25 | 25/48 | 1599/1668 | 23/0 |
| Trigonal bipyramidal / 8 | 25/25 | 25/48 | 1599/1668 | 23/0 |
| Remote ammonia / 1 | 145/5 | 145/31 | 12779/1033 | 26/16 |
| Remote ammonia / 2 | 145/5 | 145/31 | 12779/1033 | 26/16 |

The narrower filter retains no generators in the first six queries, so their 5.2–7.1% cost
buys no search reduction. It reduces the remote queries from 145 to five candidate mappings
and is 10.23–10.29 times faster. The direct filter bodies are five versus 32 lines in the trial,
including both atom and bond records. The extra 27 lines repeat focus-record lookup, ligand
mapping and frame derivation from the common check. A generator reaching the frame test gets
one temporary ligand vector. There are no new production types, traits, helpers, flags or
public interfaces. Integration remains pending in main.

**Source identities.** Both candidates use the selected single Constitution layout in a
separate checkout at `93c1fa476f9eae32fa1ec0852c9d17f8cee077c1`, with the current S5b overlay
from main HEAD `a02f535b6dab6a2570c76ab0a442d3861c61963d`. They differ only in the generator
filter for coset-change queries. The source and binary SHA-256 values are retained here:

| Candidate | Molecule-symmetry source | Timed binary |
| --- | --- | --- |
| Clear | `0c58cbfe8b925be869660ba7fca413cff29f38afa1ab18fcf93ea71c3268ed78` | `cc71a280adfe0e0b84b9885de573d4f09c6b580bc0784f1c111804f718b46239` |
| Narrow | `aa8c09d45e5837777ea556e60823db817e4cf4d0ddf1cb29375e7a752f814a54` | `d443e9443452553ed895664494937532c2b5e7a914c7d2bd6cda35c488d54bb2` |

The control build's final metadata snapshot originally picked up the next variant's source
hash because its copy started before metadata recording finished. That bookkeeping error was
corrected against the frozen control input before timing. Binary inspection found the narrowed
filter's StereoLigand mapping-allocation specialization only in the narrower binary; the
unrelated shared specialization was present in both. The build driver now captures source
hashes before compilation and verifies that they remain unchanged. Main code and documents
were unchanged by the trial; no mutating git operation was performed.

#### S5c — Derived stereo-view readings **Done**

- **Module:** `ir::view::stereo` and `ir::molecule::symmetry`.
- **Tests and benches first:** add the exact three-class and four-topicity cases on the approved
  fixtures, proper/improper ligand permutations, the kinds with more than two cosets, virtual
  ligands, undetermined cosets, kindless records and encoding failures. Test missing attached
  symmetry data and equivalence under renumbering. Write the per-record reading workloads
  before implementing these methods, using S5a's established inputs. Map the old local-group
  assertions in `tests/property/stereo/exhaustive.rs` and `semantics.rs` to exhaustive
  calculations and new readings, preserving each stated law. Show changed topicity expectations
  before editing them; any necessary change to a law is an explicit design decision.
- **Search layout:** use S5b's selected single Constitution search and verify the derived
  readings against it. Extend S5b's primitive coset-change cases for a moved focus and for a
  kind with more than two cosets to the production readings, using the selected narrower
  generator filter; compare production answers with the explicit expected values. The
  recorded comparison already includes preparation, the primitive constitutional fallback and
  repeated reads, with untimed search and allocation counts.
  Benchmark the new derived methods on the established inputs; repeat the layout comparison
  only if subsequent changes materially alter its costs. Keep experiment drivers and the
  unselected candidate in scratch, with no public selector.
- **Changes:** add `with_symmetry` to the stereo views and the three derived readings without
  a symmetry argument. The caller extracts `GraphSymmetry` from the owning view to attach it.
  Searches read the stereo view's molecule and the shared auxiliary data; they implement the coset,
  prostereogenic and topicity decisions in the design. Retain the kindless and missing-context
  contracts. Migrate those stereo property consumers after their replacements exist, and
  keep stored-coset assertions as tests of the supplied forms. No answer is recovered from a
  copied or normalized molecule.
- **Kind:** additive after the breaking removal; graph IR green, workspace still red.
  `[dep: S5b, S2a, S2b]`

**Completion, 2026-11-14.** StereoAtomView and StereoBondView have with_symmetry,
stereogenicity, topicity and ligand_symmetry. Each view holds an optional borrowed
GraphSymmetry beside its molecule reference and id; its existing producers leave the attachment
empty. The readings use the single prepared Constitution graph and the existing
molecule::symmetry::remapping_exists_at and molecule::symmetry::remapping_exists_for functions.
Those functions and GraphSymmetry::cell are crate-visible for the view consumers. The
coset-change check now has its production caller. The auxiliary fields, generator filter,
canonicalization checks, color calculation and generic search algorithm are unchanged.
Constraint-view propagation remains S6a; stereogenic_unit remains S5d.

The first 80 view cases and all 15 benchmark workloads were written before production changes.
Their initial library build failed only on the missing with_symmetry methods. They check the
approved classes and topicities, raw own cosets, proper/improper ligand actions, virtual
positions, missing attachments, kindless records, incompatible constraints, moved focus records,
renumbered fixtures, and square-planar, trigonal-bipyramidal and octahedral cosets. Nine
supplementary cases for mismatched context, unavailable topicity positions, degree/coset
encoding failures and virtual bond ligands were added after the first implementation; these
nine are not test-first evidence. All 89 pass. The tests exercise production readings against
explicit expected values; the S5a independent calculations remain separate checks.

The property consumers now enumerate label-preserving ligand permutations and apply
stereo::exhaustive::binary_coset_action independently of the production search. They retain
binary own-coset independence, the generated element/coset distribution, all 256 four-label
patterns for each kind and own-coset form, and the stored-coset assertions on the supplied
views. The retired iteration flag adds no domain variation and is removed. The cis/trans
expectation is now Prostereogenic when exactly one endpoint has equal ligands, as stated before
editing the cases: replacing one ligand then makes the site stereogenic. Equal pairs at both
ends remain Nonstereogenic. Constitutionally distinct ligand pairs use S5a's approved Distinct
expectation. Literature counts and canonicalization laws are unchanged.

The all-features graph-IR gate passes: 12,963 library tests (three ignored), 52 integration and
compile-fail tests, 576 property tests (one ignored), and one doctest (three ignored). Clippy
passes with all features and all targets under -D warnings. Its fixes include a Range import in
the new tests and two boolean assertions in the GraphSymmetry view tests; the five affected
view cases pass afterward. Formatting and git diff --check pass. Workspace restoration remains
S5f's boundary.

**Per-record timing.** One bounded Criterion pass on arm64 macOS 15.7.3 with rustc 1.99.1,
default graph-symmetry config and the bench profile measured these medians in microseconds:

| Fixture | stereogenicity | topicity | ligand_symmetry |
| --- | ---: | ---: | ---: |
| Chlorofluoropropane, CH₂ record | 33.2364 | 10.8579 | 5.39852 |
| Meso dichlorodifluoroethane, first record | 0.0127228 | 0.00459169 | 4.67822 |
| Square-planar Ma₂b₂ | 16.122 | 8.01221 | 3.13394 |
| Trigonal-bipyramidal Pa₄b | 744.753 | 25.6256 | 3.87548 |
| Cis/trans FClC=CFCl | 6.68062 | 5.8237 | 8.11575 |

Each sample measures one derived reading on an already attached view. Molecule construction,
GraphSymmetry preparation and attachment are outside the timed region. Searches, temporary
allocations and checks needed by the reading are inside it. The existing counting allocator
remains installed with counting disabled. The meso class and topicity workloads finish at the
cell exits. Topicity positions are [2, 3], [0, 1], [0, 1], [1, 2] and [0, 2] in table order.
Ligand symmetry requests proper identity for atoms and proper endpoint-block exchange for the
bond. These costs establish the new workloads; they do not compare the removed local-group
computation or the complete operation pipeline. The previous layout and generator-filter
comparisons were not repeated.

The command used 40 samples, 0.1-second warmup, a 0.4-second measurement target and 1,000
bootstrap resamples. Criterion extended the trigonal-bipyramidal class collection to about
0.615 seconds to obtain 40 samples. The entire command, including the release build, took
50.1 seconds:

```text
CARGO_TARGET_DIR=/Users/dr/.cargo-target cargo bench --locked -p umol-graph-ir --bench constraint_readings -- '^stereo_derived_readings/' --warm-up-time 0.1 --measurement-time 0.4 --sample-size 40 --nresamples 1000 --noplot --save-baseline s5c-readings-2026-11-14
```

The source was the S5c working-tree overlay on afffc23f7b499e9f9eb30dafdade4be235b3e411.
Source hashes were captured before compilation and checked unchanged afterward. SHA-256 values:

| Artifact | SHA-256 |
| --- | --- |
| umol-graph-ir/src/ir/view/stereo.rs | `7f3a9b6511b99ad6f6764a86de0dcf76b020815c9d043781cb5ff779425d4920` |
| umol-graph-ir/src/ir/molecule/symmetry.rs | `0df6aeef62cdccd7dacc6ede265911d1a39edf9844658e06ada5c311a30ea5e5` |
| umol-graph-ir/benches/constraint_readings.rs | `fa8ff2c8b0adb335a737bc1dba2e82dc5f8b33e1510c7098fd3647f54d21e27a` |
| Timed constraint_readings binary | `d7c8bcd261a074acb55c9cdd98b3b2ce75328a3fc4af8ae47c78286434772011` |

#### S5d — Stereogenic units

- **Module:** `ir::molecule::symmetry` and `ir::view::stereo`.
- **Tests and benches first:** establish the unit of every approved Razinger candidate against
  the published classification and exhaustive assignment enumeration. Cover all four
  variants, literal/undetermined/mixed stored cosets, a single record on one cyclobutane center,
  square-planar and trigonal-bipyramidal fixtures, encoding failure and the assignment limit
  below structure k's eight assignments. Write the unit-reading benchmark on the prepared
  S0a inputs before the new computation.
- **Changes:** add the frame-action reading used by the unit procedure, implement its four
  steps and `stereogenic_unit` on the stereo views. Compute all records together on first use
  and cache the result in the auxiliary data, using the supplied molecule for the computation.
  The stored configurations do not decide the units. The assignment limit retains the result
  so far as specified; it does not add an error surface.
- **Kind:** additive; graph IR green, workspace still red. `[dep: S5c, S3c]`

#### S5e — Rust operation consumers

- **Module:** `umol-graph` stereo perception, validation, transformation and enumeration,
  their configs/models, tests, properties and `benches/transform/stereo.rs`.
- **Tests and benches first:** state the operation outcomes before changing consumers: issue
  43's add/enumerate/add sequence, the default selection of configurational units, removal on
  cis,trans-1,3,5-trimethylcyclohexane, all four config defaults, intersecting selectors with
  unrestricted and empty cases, and the published stereoisomer counts grouped by canonical
  equality. Show old/new expected values where behavior changes.
  Keep the S0b and S5a workloads and inputs ready before migration; adapt their calls with the
  API change so each still measures its stated workload.
- **Changes:** build one auxiliary value per input or enumeration candidate, extracting it
  from the owning adapter for stereo views. Use the derived readings while keeping the
  operation's assertion comparison and diagnostics. Select by `stereogenic_units`, return
  `StereogenicUnit` from `StereoPerceiver::classify`, and migrate the selection's constructor
  and `matches` argument with its field. Use the agreed adder defaults and replace
  `para_stereo_max_iterations` by `max_assignments` in the operation configs and in
  `StereoPerceiver::recognize` and `classify`. Remove `StereoModel::para_stereo` and the
  enumerator's per-site cache. Keep the existing recognition kinds and other operation contracts.
- **Kind:** breaking; graph IR and graph green, workspace red until S5f. `[dep: S5d]`

#### S5f — Python follows the Rust removals

- **Module:** `umol-py` model, stereo perception, selection and transformation bindings and tests.
- **Tests first:** write the Python counterparts of the changed model/config/selector and
  classification contracts before their bindings. Preserve the corresponding Rust ownership
  and error behavior; name each required Python/Rust difference explicitly.
- **Changes:** expose `StereogenicUnit`, change `StereoSelection` to `stereogenic_units`, remove
  `para_stereo`, replace the iteration limits by `max_assignments`, and migrate the existing
  operation bindings and tests. The new graph-symmetry/view surface is added in S6b.
- **Gate:** rebuild the extension in `umol-py/.venv`, then run the affected Python tests and
  Rust binding tests. Restore the workspace build and tests before S5g.
- **Kind:** breaking, red-to-green for the workspace. `[dep: S5e, S3c]`

#### S5g — Query and operation gate

- **Module:** the graph-IR, graph and Python consumers and their tests and benchmarks.
- **Check:** run the previously written graph-IR and graph acceptance suites, their explicitly
  enabled properties, and the graph conformance suite. Rebuild and run Python after S5f. Build
  the affected benchmark targets and compare the fixed reading workloads with
  `pre-246-readings` and the operation workloads with `pre-246`. Recheck canonicalization
  against S4e if any later change reached its shared code. Report changed-workload costs
  separately from unchanged canonicalization costs. Resolve discrepancies before S6.
- **Kind:** verification (green). `[dep: S5f]`

### S6 — Keyed readings and the Python symmetry surface

S6a changes the Rust keyed behavior and can leave its Python consumers red until S6b supplies
the attached-data surface and migrates their assertions. Both Rust and Python are green at the
stage boundary.

#### S6a — Stereo constraint views and keyed readings

- **Module:** `ir::view::constraints`, `ir::view::stereo`, their tests and
  `umol-graph-ir/benches/constraint_readings.rs`.
- **Tests and benches first:** state exact results for each derived key with attached data and
  the missing-context failure before changing the keyed methods. Keep fluxionality and
  kindless cases vacuous. Cover `derived_complete`, `evidenced`, `evidenced_complete`,
  `satisfies` and `is_compatible`, and both constraints settings. Show the old vacuous derived
  results beside the new literal values before editing assertions. Add the keyed workloads
  and record their old timings as `pre-246-keyed` before changing their behavior; the old and
  new work performed by each reading is stated with the results.
- **Changes:** add `with_symmetry` to the two stereo constraint views and propagate attached
  data when a stereo view yields its constraints view. Derive the three keys through the
  molecule-backed readings. Preserve the existing lattice semantics of evidence and assertion
  comparison; those operations do not require changing the raw symmetry predicates.
- **Gate:** run the affected tests and property cases, build the benchmark and measure the new
  keyed workloads. Recheck graph consumers if the shared view plumbing changes.
- **Kind:** breaking; Rust green, Python keyed consumers migrate in S6b. `[dep: S5g]`

#### S6b — Python symmetry data, view and readings

- **Module:** `umol-py` molecule, stereo and constraint views, the corresponding tests and the
  existing ring bindings as the ownership/lifecycle reference.
- **Tests first:** cover all config fields and defaults, including `include_constraints`; the
  factory's adapter return; owned/borrowed pairing, extraction and compatible data reuse; all
  five molecule-wide queries and their view equivalents; stereo `with_symmetry`, the four
  derived readings and the keyed readings. Cover missing context, consumed views and molecule
  mutation using the existing error vocabulary. Map the old vacuous expectations in
  `tests/test_constraint_readings.py` to attached-data readings and missing-context failures,
  preserving the kindless cases. Tests are written before the bindings.
- **Changes:** bind `GraphSymmetry`, `GraphSymmetryConfig` and `GraphSymmetryView` with the
  settled Rust surface. `Molecule.graph_symmetry` returns the view. Follow `RingViews` for
  molecule ownership, owned/borrowed auxiliary storage, view counters, extraction costs and
  consumption. Attach the extracted data to stereo/constraint views. A derived stereo reading
  without attached data raises `ValueError`; kindless keyed readings retain their vacuous
  behavior. Use `from_rust` and `to_rust` where Rust/Python conversions are needed.
- **Gate:** rebuild with Python 3.13 and run Rust binding tests and the Python counterparts of
  S5b-S5d and S6a. Compare the public symbols, names, semantics and lifecycle against Rust.
- **Kind:** additive surface plus keyed behavior change, green at completion.
  `[dep: S6a, S5f]`

### S7 — Guides, final comparisons and closeout

#### S7a — Guides and interface reconciliation

- **Module:** `docs/development/nomenclature.md`, `data-types.md`, `operations.md`,
  `python-api.md`, `property-tests.md`, `docs/chemistry/representation.md`, `README.md` and
  `umol-graph-ir/spec/umol-dsl-spec.md` where already implemented enum terminology is stated.
- **Changes:** describe the one color definition, caller-owned normalization, constraint
  selection, shallow application, the molecule/data/view relationship and the RingViews-style
  lifecycle. Carry the class/unit/transformer decisions into the guides. Reconcile the
  description-level entries: levels select entity kinds, canonicalization includes their
  constraints, and graph readings choose whether to compare entity constraints. Replace the
  obsolete coloring-features examples. Check the eight form and span implementations, Rust
  and Python exports, consumers, examples and DSL names against the agreed surface. Public
  rustdoc states the operational domain and comparison relation of each relevant law.
- **Kind:** documentation (green). `[dep: S6b]`

#### S7b — Final performance record

- **Module:** this record and the affected graph-IR and graph benchmarks.
- **Check:** collect the comparisons from S4e, S4f, S5g and S6a, distinguishing canonicalization,
  setup, fixed query workloads and operation workloads. Reuse valid measurements; rerun only
  groups affected by subsequent implementation changes. Record medians and conclusions here.
  Record source revisions, benchmark IDs, commands, relevant medians and conclusions inline;
  temporary logs are not persistent evidence references. New correct stereo behavior does not
  establish a runtime improvement by itself.
- **Kind:** verification (green). `[dep: S7a, S4e, S4f, S5g, S6a]`

#### S7c — Closeout

- **Module:** the complete affected workspace surface and this discussion's status record.
- **Check:** review the complete implementation diff against the settled public-symbol and
  ownership contracts, confirm that `ir::coloring` and stale copy/normalize paths are gone,
  and check the discussion's status language. Run the final gates below. Record the results
  and set the document and index status only when the agreed implementation scope is complete.
- **Kind:** verification and documentation (green). `[dep: S7a, S7b]`

### Stage gates

S4f, S4g, S5b and S5c are complete. The sustained breaking migration continues with S5d, and S5f
restores the workspace; S5g is its acceptance
and performance gate. S6 ends green with keyed and Python parity complete. No later stage is
used to repair an earlier stage's compilation.

During implementation, run the tests of the affected crate/subitem, then the stage gate; do not
repeat unrelated full suites for a documentation-only subitem. Explicit feature-gated targets
include:

```text
cargo test --locked -p umol-graph-ir --features proptest --test property
cargo test --locked -p umol-graph --features proptest --test property
cargo test --locked -p umol-graph --features conformance --test resolution
```

At each completed implementation stage, the workspace tests pass with the Python environment
active. At closeout, use the all-features workspace gate, rebuild the extension and run Python:

```text
cargo +nightly fmt --all
source umol-py/.venv/bin/activate && python --version
source umol-py/.venv/bin/activate && cargo test --locked --workspace --all-features
source umol-py/.venv/bin/activate && cargo clippy --locked --workspace --all-features --all-targets -- -D warnings
source umol-py/.venv/bin/activate && maturin develop --locked --manifest-path umol-py/Cargo.toml
source umol-py/.venv/bin/activate && pytest -q umol-py/tests
```

### Critical path

S0-S3, S4a-S4e and S5a have implementation records. S4f is integrated and its standalone gate
passes. The remaining path is:

S5d -> S5e -> S5f -> S5g -> S6a -> S6b -> S7a -> S7b -> S7c.

S4f replaces the selector after S4e's measurements; S4g applies the later approved names and
the repository-required color-test corrections. The single Constitution search and narrower
generator filter from the recorded comparisons and the other settled review corrections are
integrated by S5b. S5c supplies the implemented derived readings and property migration. The reviewed SIN
questions are settled as recorded in Review dispositions.
The completed trial fixtures and measurements precede the port. New required tests are written
before their production changes; the existing late-written cases retain their actual history.
No remaining stage is deferred: the shared colors, migrated consumers, keyed readings and Python
surface are all part of this work.
