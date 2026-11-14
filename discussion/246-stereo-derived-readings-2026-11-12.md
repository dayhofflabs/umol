# 246 — Derived stereo readings

Status: Proposed
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

## Current state

The readings are derived in `umol-graph-ir/src/ir/molecule/symmetry.rs` through two public helper
values. S1 moved the implementation without changing these computations.

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

**`AsColor`.** A color consists of the inherent fields of the form and of its constraints. The
trait gets a second method, which returns the same color type with an empty constraints
container:

```rust
pub(crate) trait AsColor {
    type Color<'a>: Ord + Hash
    where
        Self: 'a;

    /// The inherent fields and the constraints.
    fn as_color(&self) -> Self::Color<'_>;

    /// The inherent fields, with an empty constraints container.
    fn as_inherent_color(&self) -> Self::Color<'_>;
}
```

- A constraints form is a container and can be empty, so "not compared" has no representation
  of its own. The color of a stereo atom or stereo bond holds the kind and the stereogenicity
  assertion, and an absent assertion is `StereogenicityForm::Undetermined` already.
- There is one color type per form, which states what a color consists of. Its fields are
  compared in the order of the form, so `as_color` gives canonicalization the order of colors it
  has today.
- `initial_colors` and the graph-symmetry checks take `include_constraints` as a parameter.
  `initial_colors` picks the method and ranks the result on one code path, because the type is
  the same for both. The molecule and reaction-span canonicalization callers pass `true`.
  Their stereo checks continue to compare the whole forms with `normalized_eq`.
- Canonicalization reads colors directly in 30 places, and none of them changes. A parameter on
  `as_color` in place of the second method would put an argument at each of them.

**Rejected: constraints never part of a color.** A color would be the inherent fields for every
user, and the check would compare the constraints where they count. Canonicalization would then
start from coarser cells for every molecule with entity constraints and search more, and its
keys would have to carry the constraints separately.

**Field switches.** The switches of `MoleculeColoringFeatures` have no replacement. No operation
selects fewer fields than all of them. The orbit test of `tests/canonicalization.rs` reads
orbits without isotopes, and one unit test selects no field. A selection of fields could be
added to the one definition later, by dividing the inherent part of a color.

### `GraphSymmetry`

`GraphSymmetry` is auxiliary data computed once for a molecule. It holds the incidence graph,
the ranked entity and incidence colors, the search adapter and refined cells, the operation
config and the caches described below. It holds no molecule, borrowed or copied. A reading
runs a search when it is asked, unless an exit by cell decides it.

`GraphSymmetryView` follows `RingViews`. `Molecule::graph_symmetry` borrows the molecule and
returns a view owning the computed data. `from_graph_symmetry` borrows existing data and the
molecule without recomputing or copying either. `into_graph_symmetry` consumes the view: it
moves owned data or clones borrowed data. The owning view constructor is crate-private.

The five molecule-wide queries have the molecule as their receiver and take the auxiliary
data. The view exposes the same query names and forwards to those molecule operations. Every
search check reads that borrowed molecule. The stereo views likewise already borrow their
molecule; attaching the auxiliary data gives their readings both inputs.

```rust
pub struct GraphSymmetryConfig {
    pub automorphism_algorithm: AutomorphismAlgorithm,
    pub max_assignments: u32,
    pub include_constraints: bool,
}

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
- **Malformed molecules.** `graph_symmetry` remains infallible. `initial_colors` can fail when
  literal electron counts do not match their aromatic or multicenter participants, or when a
  stereo kind, frame or literal coset cannot support the required action. Compatibility checks
  of entity constraints apply only when `include_constraints` is `true`. For an input that
  cannot be encoded, every entity gets a number of its own in place of a color. No two atoms
  then share a cell, and the exits decide the readings: every atom is alone in its orbit, each
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
- **Caches.** The incidence graph and the cells of the Constitution level, which the `Distinct`
  test needs, are built on first use and kept in the value. The stereogenic units are the one
  answer that is kept, because they are computed for all records together. No other answer is
  kept.

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
remapping of the molecule onto itself that passes a check. Source and target of the search are
the molecule itself. The color of a stereo entity holds its kind, with `include_constraints`
also its stereogenicity assertion, and not its coset, so the refined cells do not depend on the
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
  action before applying it once to the source coset. Constraints receive only the frame action.
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
- **Prostereogenic.** For a stereo entity that is not stereogenic, the same decision is repeated
  once per explicit ligand with that ligand held fixed. If X comes out stereogenic for some
  ligand, it is `Prostereogenic`. Otherwise it is `Nonstereogenic`.
- **Topicity.** The three searches run in the order of the table, and the first that finds a
  remapping decides. The third runs on the incidence graph of the Constitution level, which has
  no stereo entities.
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

The two search entry points take the molecule and its auxiliary symmetry data and call
`remapping_exists`. Their graph checks compare transported representations directly. The
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

One `GraphSymmetry` per molecule carries the incidence data, colors and refinement. It does not
pay for molecule copying or normalization. A view borrows the molecule and either owns or borrows
that auxiliary data. Extracting owned data moves it; extracting borrowed data clones the auxiliary
data. A reading costs an exit by cell or one or more searches. The stereogenic units are computed
together and cached.

The timing trials below used a setup that retained a normalized molecule. Their measurements
describe those trials, not the cost of the revised setup and checks. The correctness trial of
shallow application below did not measure time. Required benchmarks precede the changes in S4a,
and comparisons follow the shared-color work and the consumer migration. Setup-only timings are
reported separately: the old setup precomputes chirality and orbits, while the new setup prepares
queries. The comparison also measures a fixed build-and-read workload.

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

A standalone crate in `scratch/246-normalization-trial` used the current workspace crates as
dependencies, without a repository copy or production edits. The final run passed all 39 cases.
The [report](../scratch/246-normalization-trial/report.txt) records its scope and command.

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

An initial interpretation contracted an operator already inside the supplied term. It lost raw
syntax in `Mirror(Mirror(x))` across an action and its inverse. That interpretation was rejected.
The final checks assert the complete newly constructed expression trees and separately check
their value equivalence; exact operand preservation is checked when the new actions compose to
the identity. This does not impose normalization on graph inputs or alter the public
`FrameTransport` contract.

The trial checked the application and comparison decisions using direct enumeration. It did
not implement the new `GraphSymmetry`, its view or the shared `AsColor` path, and it measured no
timings. The implementation tests and measurements below precede those changes.

## Staged implementation plan

The remaining work starts at S4. Its foundation is one color definition through `AsColor` and
shared search-color construction. S0 to S3 retain their completed scope. Each subitem carries its
own tests under the test-writing skill. Required tests, independent expected values and benchmark
harnesses are written before the implementation they exercise. An intended initial test failure
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
chiral dichlorodifluoroethane fixtures, chlorofluoropropane, and the meso and chiral
trichloropentane fixtures. Three more inputs copy the meso dichlorodifluoroethane fixture:
`readings_raw_term` stores `Term(Lit(0))` at its first center; `readings_raw_singleton` stores a
singleton carbon element and singleton coset there; `readings_constraints` stores singleton
valence 1 on H2 and singleton stereogenicity `Stereogenic` on the first stereo record.

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

The [run record](../scratch/246-s4a/context.json), [input manifest](../scratch/246-s4a/input-manifest.txt)
and [full medians](../scratch/246-s4a/baseline-results.json) retain the commands and machine-readable
results. Source hashes match before and after both runs. The graph-IR package tests pass
(12,538 library tests, three ignored, plus integration and documentation tests). Both benchmark
targets pass the preflight checks, and clippy with `--all-targets -- -D warnings` passes for
`umol-graph-ir` and `umol-graph`. Formatting and the diff check pass.

#### S4b — `AsColor::as_inherent_color` **Done**

- **Module:** `ir::traits`, the eight entity-form modules and `ir::span`.
- **Tests first:** add exact cases for both methods on each form. `as_color` retains the current
  order and constraint distinctions; `as_inherent_color` omits only constraints and has the same
  color type. Cover empty and nonempty stores, raw singleton/set expressions, each span variant,
  and a modified span whose sides differ only in excluded constraints. Check raw input retention
  and equal-color hash consistency.
- **Changes:** add the second trait method and its implementations. Represent the excluded
  constraints by the empty container, or the existing undetermined stereogenicity assertion for
  stereo colors. Preserve the declaration/comparison order of the existing full colors. Neither
  method normalizes a form. The trait and color types keep their existing crate visibility.
- **Kind:** additive (green). `[dep: S4a]`

The method has implementations for all eight forms and `EntitySpan`. Each result borrows the
inherent fields and uses an empty constraint store of the same color type. The empty stores are
shared per type; stereo colors use `StereogenicityForm::Undetermined`. Full colors borrow the
stored constraints. The four former whole-form wrappers now name their borrowed fields in the
same declaration order as their forms. The stereo and positional projections keep their existing
full-color definitions.

The 42 new cases cover empty and nonempty stores, contradictory and raw constraints, raw scalar
and coset forms, all four span variants, and collapse of a modification confined to constraints.
Both color readings have equal-color hash checks. The cases assert that the supplied forms stay
structurally unchanged. S4c enables the new method and its empty-store imports for production use.

The graph-IR package tests pass: 12,580 library tests, three ignored, plus integration and
documentation tests. Clippy passes for all targets without warnings. Formatting and the diff
check pass.

The 24 `pre-246-colors` cases were compared with the same inputs, Nauty, 100 samples, 0.5 s
warm-up and 1.5 s measurement time. The original baseline estimates and benchmark source match
their saved hashes. Median changes in the first comparison are:

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

The [comparison record](../scratch/246-s4b/comparison-results.json) retains both runs' medians.

#### S4c — Shared search colors and canonicalization callers **Done**

- **Module:** `ir::molecule::symmetry`, `ir::reaction_span::symmetry`, and the two owning
  canonicalization modules.
- **Tests first:** extend `initial_colors` cases for both values of `include_constraints`, for
  all entity kinds and spans. Check incidence roles and per-participant electron counts, raw
  literals versus singleton sets, and constraint distinctions. Excluded constraints must not
  be normalized or checked for frame compatibility. Included constraints are read as stored.
  Keep the existing incompatible electron-count, kind, frame and literal-coset cases. Establish
  unchanged canonical representatives and comparison results before changing the call sites.
- **Changes:** add `include_constraints` to the shared color constructors and choose between the
  two `AsColor` methods on the same ranking path. Enable `as_inherent_color` and its empty-store
  imports outside tests as these production callers land. Limit stereo compatibility checks to
  the data selected by the caller. Canonicalization passes `true` at every molecule and reaction-span
  call and retains its own upfront normalization. Its direct `as_color` reads, leaf ordering,
  normalized stereo checks and generator filters retain their behavior.
- **Kind:** breaking internally, red-to-green within the subitem as all call sites migrate.
  `[dep: S4b]`

Both shared color constructors now take `include_constraints` and select the corresponding
`AsColor` reading on the existing ranking path. Excluded stereo constraints are neither
transported nor checked. Configuration compatibility remains checked in both modes, on each
present side of a span. The electron-count, site-kind, frame and outer literal-coset checks
remain. Neither constructor normalizes the supplied forms.

Canonicalization changes only at its nine call sites, which pass `true`. Its normalization,
direct color reads, leaf ordering, normalized comparisons and generator filters are unchanged.
The trait method and empty-store imports are enabled outside tests.

The 224 color cases cover both modes, all entity kinds, span variants, raw forms, incidence
roles, electron counts and incompatible selected data. Excluded incompatible constraints
succeed; incompatible configurations fail in both modes. Raw contradictory constraints are
read as stored when included. The existing canonicalization expectations are unchanged.

The [all-features graph-IR test run](../scratch/246-s4c/package-tests.txt) passes: 12,732 library
tests, three ignored, plus 559 property tests, one ignored, and the integration and documentation
tests. [Clippy](../scratch/246-s4c/clippy.txt) passes with all features and targets and warnings
denied. Formatting and the diff check pass. The benchmark comparison remains the S4e gate,
after S4d.

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

The [all-features graph-IR test run](../scratch/246-s4d/package-tests.txt) passes: 12,791 library
tests, three ignored, plus 559 property tests, one ignored, and the integration and documentation
tests. [Clippy](../scratch/246-s4d/clippy.txt) passes with all features and targets and warnings
denied. Formatting and the diff check pass. S4e supplies the benchmark comparison and
dependent-crate gate.

#### S4e — Shared-foundation gate

- **Module:** the S4b-S4d modules, their tests and the canonicalization benchmarks.
- **Check:** run the affected form, span, color and canonicalization tests and the explicitly
  enabled graph-IR property suite. Compare the canonicalization groups against `pre-246` and
  the added groups against `pre-246-colors`. Report medians and affected cases; a reproducible
  regression in unchanged canonicalization behavior is resolved before S5 begins. Run the
  dependent-crate stage gate after the internal signature migration. No new implementation is
  added in this subitem.
- **Kind:** verification (green). `[dep: S4c, S4d]`

### S5 — Symmetry data, molecule-backed queries and consumer migration

This is the breaking stage. S5a establishes evidence before the replacement. S5b retires the old
public surface; the workspace can remain red through S5f while its consumers migrate. The
`umol-graph-ir` library tests are green after S5b; its integration and property suites follow
after S5c supplies the stereo-view readings. `umol-graph` is green after S5e; S5f restores the
Rust/Python workspace. Each implementation subitem writes its required tests before its
production changes. S5g is the green stage boundary.

#### S5a — Query oracles and workload controls

- **Module:** `ir/molecule/symmetry/tests.rs`, the existing stereo property/reference modules,
  `umol-graph/benches/transform/stereo.rs` and `umol-graph-ir/benches/canonicalize.rs`.
- **Tests and benches first:** bring in the needed approved S0a fixtures beside their tests.
  Establish expected remappings, chirality, complete orbit vectors, classes and published unit
  counts through exact cases and bounded independent enumeration. Include the dative donor and
  acceptor, the five-carbon aromatic system with electron counts 2, 1, 1, 1, 1, and a nonuniform
  multicenter case. The query assertions are written first in their owning implementation
  subitems; their expected values and reference inputs are established here.
- **Benchmark controls:** check that S4a's `pre-246-readings` cases still measure the fixed
  build-and-read workloads before the API migration. Preserve their inputs and scope as their
  calls migrate. Do not replace that baseline after the shared transport change: the old
  symmetry implementation also calls `coset_apply_permutation`. Keep the 51 operation cases of
  `pre-246` and report setup-only and fixed-workload comparisons separately.
- **Kind:** additive (green). `[dep: S4e, S0a, S0b, S3b]`

#### S5b — `GraphSymmetryConfig`, auxiliary data and molecule queries

- **Module:** `ir::molecule::symmetry`, `ir::view::symmetry`, `ir::view`, `ir::coloring` and
  the affected exports and graph-IR tests/benchmarks.
- **Tests first:** cover config defaults and both constraints settings; the five molecule
  queries and the same five view queries on S0a fixtures and their renumberings; owned and
  borrowed views; extraction, source/data independence and compatible reuse. Test the raw
  literal/singleton distinction, excluded contradictory or frame-incompatible constraints,
  included raw constraints, molecule-level constraints and encoding failures. Test the exact
  shallow symbolic operands and composed mirror/frame actions from the correctness trial.
  Check that mismatched context does not panic without promising a semantic result.
- **Changes:** replace the old prepared value with auxiliary incidence data, shared colors,
  adapter, cells, config and caches. It holds no molecule. Make the config nongeneric, with
  `automorphism_algorithm`, `max_assignments` and `include_constraints`. Build colors from the
  supplied forms through S4c, and remove `ir::coloring` and its three types when this consumer
  switches. Remove their exports and feature-selector callers.
- **Interface:** `Molecule::graph_symmetry` returns `GraphSymmetryView`, owning the auxiliary
  data and borrowing the molecule. The view's owning constructor is crate-private;
  `from_graph_symmetry` borrows data and molecule, and `into_graph_symmetry` moves owned data
  or clones borrowed data. The auxiliary data supports that explicit cloning and has no
  arbitrary-parts constructor. The five query names move to `Molecule`, taking the data, and
  the view forwards them. Export the adapter with the other views and from `ir`.
- **Checks:** implement the graph-operation remapping checks and search entry points with the
  molecule supplied by the query. Compare selected transported forms structurally, compose the
  new mirror and frame actions before shallow application, and skip excluded constraint
  transport. Preserve the existing normalized-comparison path for canonicalization. Use the
  existing generic search and the established cell exits. Build Constitution-level data and
  cells on first use. Encoding failure uses distinct colors; lack of a normal form is not a
  reason for that fallback.
- **Migration:** retire `StereoSymmetry`, the two molecule methods returning it, the forwarding
  view methods and their conveniences with the old derivation. Remove their graph-IR tests
  only after recording which definitions the replacement cases state. Migrate the molecule
  query consumers in `tests/canonicalization.rs`, the config use in `tests/integrity.rs` and
  affected graph-IR benchmarks; remove the isotope-ignoring orbit selection. The stereo
  property consumers migrate in S5c, when their replacement readings exist.
- **Kind:** breaking; graph-IR library tests green, integration/property consumers can remain
  red until S5c and the workspace until S5f.
  `[dep: S5a, S4c, S4d]`

#### S5c — Derived stereo-view readings

- **Module:** `ir::view::stereo` and `ir::molecule::symmetry`.
- **Tests and benches first:** add the exact three-class and four-topicity cases on the approved
  fixtures, proper/improper ligand permutations, the kinds with more than two cosets, virtual
  ligands, undetermined cosets, kindless records and encoding failures. Test missing attached
  symmetry data and equivalence under renumbering. Write the per-record reading workloads
  before implementing these methods, using S5a's established inputs. Map the old local-group
  assertions in `tests/property/stereo/reference.rs` and `semantics.rs` to the independent
  reference and new readings, preserving each stated law. Show changed topicity expectations
  before editing them; any necessary change to a law is an explicit design decision.
- **Changes:** add `with_symmetry` to the stereo views and the three derived readings without
  a symmetry argument. The caller extracts `GraphSymmetry` from the owning view to attach it.
  Searches read the stereo view's molecule and the shared auxiliary data; they implement the coset,
  prostereogenic and topicity decisions in the design. Retain the kindless and missing-context
  contracts. Migrate those stereo property consumers after their replacements exist, and
  keep stored-coset assertions as tests of the supplied forms. No answer is recovered from a
  copied or normalized molecule.
- **Kind:** additive after the breaking removal; graph IR green, workspace still red.
  `[dep: S5b, S2a, S2b]`

#### S5d — Stereogenic units

- **Module:** `ir::molecule::symmetry` and `ir::view::stereo`.
- **Tests and benches first:** establish the unit of every approved Razinger candidate against
  the published classification and the independent assignment reference. Cover all four
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
- **Check:** collect the comparisons from S4e, S5g and S6a, distinguishing canonicalization,
  setup, fixed query workloads and operation workloads. Reuse valid measurements; rerun only
  groups affected by subsequent implementation changes. Record medians and conclusions here.
  New correct stereo behavior does not establish a runtime improvement by itself.
- **Kind:** verification (green). `[dep: S7a, S4e, S5g, S6a]`

#### S7c — Closeout

- **Module:** the complete affected workspace surface and this discussion's status record.
- **Check:** review the complete implementation diff against the settled public-symbol and
  ownership contracts, confirm that `ir::coloring` and stale copy/normalize paths are gone,
  and check the discussion's status language. Run the final gates below. Record the results
  and set the document and index status only when the agreed implementation scope is complete.
- **Kind:** verification and documentation (green). `[dep: S7a, S7b]`

### Stage gates

S4 ends with the shared-color work green. S5b begins the sustained breaking migration and S5f
restores the workspace; S5g is its acceptance and performance gate. S6 ends green with keyed and
Python parity complete. No later stage is used to repair an earlier stage's compilation.

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

S0-S3 and S4a-S4d are complete. The remaining path is:

S4e -> S5a -> S5b -> S5c -> S5d -> S5e -> S5f -> S5g -> S6a -> S6b -> S7a -> S7b -> S7c.

S4e checks the shared color calculation and shallow application before their graph-symmetry
consumer is replaced. Baselines precede each affected implementation, and every required test
is written before the production change in its subitem.
No remaining stage is deferred: the shared colors, migrated consumers, keyed readings and Python
surface are all part of this work.
