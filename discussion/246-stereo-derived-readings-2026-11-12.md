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
- where each reading is asked, and
- how the readings are computed.

## Current state

The readings are derived in `umol-graph-ir/src/ir/symmetry.rs` through two public helper values.

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
group. `stereogenicity` returns `Stereogenic` in that case, `Prochiral` when some pair of ligand
positions is enantiotopic, and `Symmetric` otherwise. `topicity` returns `Homotopic` for two
positions in one proper orbit of the local group, `Enantiotopic` for two positions that only an
improper element relates, and `Diastereotopic` for every other pair, including two ligands of
different constitution.

The readings are consumed as follows.

- `StereoPerceiver::classify`, the conformance validator and the stereo enumerator build a
  `GraphSymmetry`, derive a `StereoSymmetry` per stereo entity, and call `stereogenicity`,
  `topicity` and `group().contains(..)` on it. The validator and the enumerator compare the result
  with the assertion themselves.
- The stereo views have methods that take a `&StereoSymmetry` and forward to it: `stereogenicity`,
  `is_stereogenic`, `is_prochiral`, `topicity`, `is_homotopic`, `is_enantiotopic`,
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

Every test in `symmetry.rs` uses a molecule with a single stereo entity.

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
operations are those that keep the stereo entity in place, as today.

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

### `GraphSymmetry`

`GraphSymmetry` is an owned value built once per molecule, as a `RingSet` is. It can be stored and
asked without the molecule. It holds the prepared input of the search: a normalized copy of the
molecule without its constraints, the incidence graph, the colors and the refined cells. A
reading runs a search when it is asked, unless an exit by cell decides it.

```rust
pub struct GraphSymmetryConfig {
    pub automorphism_algorithm: AutomorphismAlgorithm,
    pub max_assignments: u32,
}

impl Molecule {
    pub fn graph_symmetry(&self, config: &GraphSymmetryConfig) -> GraphSymmetry;
}

impl GraphSymmetry {
    pub fn is_chiral(&self) -> bool;
    pub fn same_proper_orbit(&self, a: AtomId, b: AtomId) -> bool;
    pub fn same_star_orbit(&self, a: AtomId, b: AtomId) -> bool;
    pub fn proper_orbit_of(&self, a: AtomId) -> Vec<AtomId>;
    pub fn star_orbit_of(&self, a: AtomId) -> Vec<AtomId>;
}
```

- **Config.** `GraphSymmetryConfig` holds the automorphism algorithm and the limit
  `max_assignments` of the stereogenic units, and has a `Default`: the algorithm the stereo
  operation configs default to, and 4,096 assignments. The caller selects both; they are not part
  of a chemistry model. The coloring, the fixpoint switch and the iteration limit of the present
  config are removed.
- **Constraints.** The copy has neither the inline constraints of the entities nor the
  molecule-level constraints. A derived reading then does not depend on an assertion it is
  compared with. The present coloring uses inherent fields only, so the meaning is the present
  one. The constraints are removed with the existing `clear` of the constraint stores.
- **Malformed molecules.** `graph_symmetry` does not fail, as today. `normalize` and
  `initial_colors` return a contradiction for a stored form that is contradictory, for an
  aromatic system or a multicenter bond whose literal electron counts do not match its atoms,
  and for a stereo record whose kind does not fit its site or whose coset is outside the range
  of its kind. For such a molecule the value is built from the molecule as stored, and every
  entity gets a number of its own in place of a color. No two atoms then share a cell, and the
  exits decide the readings: every atom is alone in its orbit, each stereo record with a kind is
  `Stereogenic` and `Constitutional`, each pair of ligand atoms is `Distinct`, and `is_chiral`
  is decided by the identity alone. These are the readings of a molecule without symmetry. They
  are not correct for the malformed molecule, and the symmetry code does not report the
  malformation. `canonical_eq` treats two such molecules in the same way: it compares them as
  stored.
- **Validity.** The value describes the molecule as it was when the value was built. After a
  change to the molecule a new value is built. Nothing ties a value to a molecule, as for
  `RingViews::from_ring_set`.
- **Caches.** The incidence graph and the cells of the Constitution level, which the `Distinct`
  test needs, are built on first use and kept in the value. The stereogenic units are the one
  answer that is kept, because they are computed for all records together. No other answer is
  kept.

Molecule chirality and atom equivalence have no assertion in the constraint language. They are
derived readings only and are asked of the value itself.

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
readings take no further argument:

```rust
let symmetry = mol.graph_symmetry(&config);
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

The adder works as it does today with the class. It adds every recognized site to a copy of the
molecule, builds the `GraphSymmetry` of the copy, and keeps the records whose unit is selected.

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
the copy held by the `GraphSymmetry`. The color of a stereo entity holds its kind and not its
coset, so the refined cells do not depend on the configurations.

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
  the same coset. The check uses the existing `mirror` of the stereo forms, which also covers a
  coset given as a set.
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
  of the data types guide. Holding a virtual ligand fixed therefore changes nothing, and two
  virtual ligands of one kind at one site are `Diastereotopic`, as today. A virtual ligand and
  a ligand of another kind are `Distinct`.
- **Undetermined cosets elsewhere.** Stereo entities are compared as `remapping_preserves_stereo`
  compares them: an entity with an undetermined coset is carried onto one with an undetermined
  coset. The class is then a statement about the description as stored: it says whether another
  coset of X changes what the description covers. In 1,3-dichlorocyclobutane with one center
  undetermined, the other center is not stereogenic, as today.
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
| `remapping_preserves_stereo`, existing | every configuration is preserved |
| `remapping_mirrors_stereo` | every configuration goes to its mirror image |
| `remapping_changes_coset` | one stereo entity goes to a given coset, and the others are preserved |
| `remapping_frame_action` | the permutation of a record's ligands under a remapping that sends the record to itself |
| `GraphSymmetry::remapping_exists_at` | a remapping at a description level and of an `Orientation` sends given atoms to given atoms |
| `GraphSymmetry::remapping_exists_for` | a remapping that keeps given atoms fixed changes the coset of a stereo entity to a given one |

The two methods call `remapping_exists`. The readings that the stereo views call on the
`GraphSymmetry` have the names of the view methods and are visible within the crate.

### Modules

Canonicalization is organized by operation today, with `canonicalize/molecule.rs` and
`canonicalize/reaction_span.rs`, while the rest of `ir` is organized by type. Symmetry belongs to
the molecule, and canonicalization is one of its users. Each type gets its symmetry and its
canonicalization side by side, and canonicalization depends on symmetry.

| Module | Holds | Today |
| --- | --- | --- |
| `ir::entity` | `DescriptionLevel`, beside `EntityKind` | `ir/canonicalize.rs` |
| `ir::symmetry` | the generic search: the colors' ranking, the adapter, the partition and its refinement, the search, `remapping_exists` | `ir/canonicalize.rs` |
| `ir::canonicalize` | what is independent of the type: the `Canonicalize` trait, the context, the choice of level, the errors | `ir/canonicalize.rs` |
| `ir::molecule::symmetry` | the molecule's colors and remapping checks, `GraphSymmetry` and its readings | `ir/canonicalize/molecule.rs`, `ir/symmetry.rs` |
| `ir::molecule::canonicalize` | the molecule's leaf keys, canonical form and `Canonicalize` impl | `ir/canonicalize/molecule.rs` |
| `ir::reaction_span::symmetry` | the reaction span's colors and remapping checks | `ir/canonicalize/reaction_span.rs` |
| `ir::reaction_span::canonicalize` | the reaction span's leaf keys, canonical form and impls, with the `Reaction` impl that calls them | `ir/canonicalize/reaction_span.rs` |

A description level is a cumulative group of entity kinds, and the readings and canonicalization
are two of its users. `DescriptionLevel` therefore sits with the entity kinds and belongs to
neither. The symmetry modules import nothing from the canonicalization modules.

The colors and remapping checks of a type are `initial_colors`, `incidence_color`,
`remapping_preserves_stereo` and `retain_stereo_preserving_generators`. Everything moves
unchanged. Items used from another module become visible within the crate, and nothing becomes
public. The test files move with the functions they test.

### Removed

| Item | Replaced by |
| --- | --- |
| `StereoSymmetry`, `Molecule::stereo_atom_symmetry`, `Molecule::stereo_bond_symmetry` | the three derived readings on the stereo views |
| the present content of `GraphSymmetry` and its derivation | the value described above |
| the view methods that take a `&StereoSymmetry`, among them `ligand_symmetry` returning the local group | the three derived readings |
| the view conveniences `is_stereogenic`, `is_prochiral`, `is_homotopic`, `is_enantiotopic`, `is_diastereotopic` | comparison with the variant |
| `MoleculeColoring`, `ConstitutionColoring`, `MoleculeColoringFeatures` | nothing; every operation used the full coloring |
| `StereoModel::para_stereo` | nothing; selecting `Constitutional` alone leaves the units that depend on others out |
| `para_stereo_max_iterations` of the stereo operation configs | `max_assignments` |
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
  stand for `Nonstereogenic`, `Prostereogenic` and `Stereogenic`, and `Distinct` gets a fourth
  topicity glyph beside `=`, `'` and `/`.
- **Python.** The `Stereogenicity` and `Topicity` enums, `StereoSelection`, the stereo operation
  configs and `StereoModel` follow the Rust changes, and `StereogenicUnit` is added. The Python
  stereo constraints views have the keyed readings, so Python gets `GraphSymmetry`,
  `GraphSymmetryConfig`, `Molecule.graph_symmetry`, `with_symmetry`, and the four derived
  readings on the stereo views. A derived reading without the value raises `ValueError`, as a
  ring key does.
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

One `GraphSymmetry` per molecule carries the normalization, the colors and the refinement, which
the trial below calls the setup. A reading then costs nothing more when an exit decides it and
one or more searches otherwise. Today one `GraphSymmetry` per molecule and one `StereoSymmetry`
per stereo entity answer every question about that entity. The trial measures the readings; the
plan measures the operations.

The handling of a malformed molecule adds nothing for a well-formed one. `normalize` and
`initial_colors` return a `Result` in either case, and only their error case differs. The
additions are in `ir::molecule::symmetry`, about twenty lines:

| Place | Addition | Runs |
| --- | --- | --- |
| building the value | the error case numbers the entities of the stored molecule | malformed molecules only |
| the Constitution-level cells, built on first use | the same error case | malformed molecules only |
| the check with given cosets | a coset outside the range of its kind fails the check | every check that takes a mirror image |

`initial_colors`, the search, `remapping_exists` and canonicalization are not touched. This path
was not run in the trial, and that a well-formed molecule costs the same was not measured.

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

## Staged implementation plan

The plan sequences the design above. Each subitem carries its tests, written by the test-writing
skill. Where a subitem changes the expected value of an existing test or removes a test, the old
and the new statement are shown before the edit. A fixture is added to the tree when it has been
approved. Formatting is `cargo +nightly fmt --all`.

A subitem is checked with the tests of the crates it touches. A subitem that moves code without
changing it gets one build with those tests. A subitem that changes an API a benchmark uses also
builds the benchmarks of that crate, which the tests do not compile. Python is rebuilt and
tested after every subitem that touches `umol-py`.

### S0 — Fixtures and baseline

#### S0a — Fixtures **Done**

- **Module:** the test modules of S4 and the benchmark of S0b. The fixtures have no file of
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

Baseline of the new cases, medians in µs:

| Compound | `classify` | `validate` |
| --- | ---: | ---: |
| 1-chloro-1-fluoropropane | 11.6 | 12.0 |
| meso-2,3,4-trichloropentane | 18.0 | 17.8 |
| inositol | 38.7 | 14.5 |
| compound I | 43.2 | 16.9 |
| structure k | 26.0 | 10.7 |

| `stereo_enumerator` | 1-chloro-1-fluoropropane | inositol |
| --- | ---: | ---: |
| `transform_shared` | 10.7 | 42.6 |
| `transform_into_shared` | 10.8 | 44.0 |
| `transform_unique` | 10.8 | 40.7 |
| `transform_into_unique` | 10.5 | 40.2 |
| `iter_setup` | 0.50 | 40.1 |
| `iter_first` | 10.1 | 0.07 |
| `iter_prefix_4` | 39.8 | 0.12 |
| `iter_all` | 40.1 | 0.12 |

### S1 — Modules

Every subitem of this stage moves code without changing it.

#### S1a — Present symmetry code under the molecule

- **Module:** `ir/symmetry.rs` → `ir/molecule/symmetry.rs`.
- **Changes:** the file moves with its tests. The public paths `ir::GraphSymmetry`,
  `ir::GraphSymmetryConfig` and `ir::StereoSymmetry` stay.
- **Kind:** move (green).

#### S1b — Generic search

- **Module:** `ir/canonicalize.rs` → `ir/entity.rs` and `ir/symmetry.rs`.
- **Changes:** `DescriptionLevel` moves to `ir/entity.rs`, beside `EntityKind`. `InitialColors`,
  `rank_colors`, `AutomorphismColor`, `AutomorphismAdapter`, `OrderedPartition`, `same_cells`,
  `remapping_preserves_incidences`, the canonical search and `remapping_exists` move to
  `ir/symmetry.rs` with their tests. `ir::canonicalize` keeps the trait, the context, the choice
  of level, the leaf-key helpers and the errors, and imports the search.
- **Kind:** move (green). `[dep: S1a]`

#### S1c — Colors and remapping checks of the two types

- **Module:** `ir/canonicalize/molecule.rs` → `ir/molecule/symmetry.rs`;
  `ir/canonicalize/reaction_span.rs` → `ir/reaction_span/symmetry.rs`.
- **Changes:** `initial_colors`, `incidence_color`, `remapping_preserves_stereo` and
  `retain_stereo_preserving_generators` of each type move with their tests.
- **Kind:** move (green). `[dep: S1b]`

#### S1d — Canonicalization of the two types

- **Module:** `ir/canonicalize/molecule.rs` → `ir/molecule/canonicalize.rs`;
  `ir/canonicalize/reaction_span.rs` → `ir/reaction_span/canonicalize.rs`.
- **Changes:** the two files move with their tests. The tests of `ir/canonicalize/tests.rs`
  follow the functions they test.
- **Kind:** move (green). `[dep: S1c]`

#### S1e — Comparison

- **Check:** the canonicalization groups against the baseline of S0b. No change is expected.
  `[dep: S1d]`

### S2 — Class names

#### S2a — `Stereogenicity` variants

- **Module:** `ir::stereo`, then every user in `umol-graph-ir`, `umol-graph`, `umol-io` and
  `umol-py`.
- **Changes:** `Prochiral` becomes `Prostereogenic` and `Symmetric` becomes `Nonstereogenic`.
  The keywords become `:prostereogenic` and `:nonstereogenic`, and the glyphs `=`, `'` and `/`
  stay with the renamed variants. The Python enum follows, and so does
  `umol-graph-ir/spec/umol-dsl-spec.md`, which states the keywords and the glyphs.
- **Kind:** breaking (red→green within the subitem).

#### S2b — `Topicity::Distinct`

- **Module:** `ir::stereo`, `dsl::stereo`, the topicity forms, `umol-py`.
- **Adds:** the variant, its keyword `:distinct`, its glyph in the compact attribute strings, its
  Python member, and its entry in `umol-graph-ir/spec/umol-dsl-spec.md`. The glyph is chosen at
  this subitem.
- **Changes:** expected values of the lattice and normalization tests of `TopicityRelationForm`
  that enumerate all variants.
- **Kind:** breaking (red→green within the subitem). `[dep: S2a]`

### S3 — Check and partitions

#### S3a — Checks with given cosets

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

#### S3b — Partitions with single cells

- **Module:** `ir::symmetry`.
- **Adds:** no function. A given atom is made a single cell with the existing
  `OrderedPartition::individualize`, once per atom, in the partition of the source and in that
  of the target of `remapping_exists`.
- **Tests:** an atom kept fixed and an atom sent to another, with `remapping_exists`, on the
  S0a fixtures; two atoms of different colors give no remapping.
- **Kind:** additive (green). `[dep: S1b, S0a]`

#### S3c — `StereogenicUnit`

- **Module:** `ir::stereo`.
- **Adds:** the enum with its four variants. The rustdoc of `StereogenicUnit`, `Stereogenicity`
  and `Topicity` and of their variants cites the papers of the Sources section.
- **Kind:** additive (green). `[dep: S2a, S2b]`

### S4 — `GraphSymmetry` and the derived readings

The workspace is red from S4a until S4f. The unit tests of `umol-graph-ir` pass after each of
S4a, S4b and S4c, because S4a removes what reads the present derivation together with it.

#### S4a — `GraphSymmetryConfig`, `GraphSymmetry`, `Molecule::graph_symmetry`

- **Module:** `ir::molecule::symmetry`, `ir::coloring`, `ir::view::stereo`.
- **Changes:** the config holds the algorithm and `max_assignments` and has a `Default`.
  `graph_symmetry` builds the value from the molecule without constraints and does not fail. For
  a molecule that does not normalize or has no colors it builds the value from the molecule as
  stored, with a number of its own for every entity, and the Constitution-level cells follow.
  `is_chiral` and the four orbit queries run the search through the private
  `remapping_exists_at`, after the exits by cell. The caches are `OnceLock`s. The present
  derivation and the three coloring types are removed, and with them what reads the derivation:
  `StereoSymmetry`, `stereo_atom_symmetry`, `stereo_bond_symmetry`, the forwarding view methods
  and the five conveniences, and the unit tests that state the rule of the local group.
- **Tests:** `is_chiral` and the orbits on the S0a fixtures, under both numberings; a molecule
  with an inline constraint on one of two equivalent atoms; a molecule that does not normalize,
  one with electron counts of the wrong length and one with a coset outside the range of its
  kind, each with every atom alone in its orbit.
- **Kind:** breaking (red for the workspace). `[dep: S3a, S3b]`

#### S4b — Derived readings on the stereo views

- **Module:** `ir::view::stereo`, `ir::molecule::symmetry`.
- **Changes:** `StereoAtomView` and `StereoBondView` get `with_symmetry`. `stereogenicity`,
  `topicity` and `ligand_symmetry` take no symmetry argument and run the search of the design,
  with the loop over the other cosets through the private `remapping_exists_for`, the
  prostereogenic step, the exits by cell, and the Constitution-level cells built on first use.
  They panic without the value attached.
- **Tests:** each reading on the S0a fixtures: the three classes for stereo atoms and stereo
  bonds, the four topicity variants, proper and improper ligand permutations, the two kinds with
  more than two cosets, a stereo entity with an undetermined coset, an entity without a kind, a
  malformed molecule.
- **Kind:** additive within the red stage. `[dep: S4a]`

#### S4c — Stereogenic units

- **Module:** `ir::molecule::symmetry`, `ir::view::stereo`.
- **Adds:** `remapping_frame_action`, with which a remapping is checked to keep a stereo record
  at its site with its coset unchanged or changed, the four steps of the procedure, the result
  kept in the `GraphSymmetry` on first use, and `stereogenic_unit` on the stereo views.
- **Tests:** the unit of every candidate of the Razinger fixtures against the published
  classification; the four variants on the compounds of the second trial; a molecule whose
  records have assigned, undetermined and mixed cosets gives the same units; a single record on
  one center of 1,3-dichlorocyclobutane; the square-planar and trigonal-bipyramidal fixtures;
  structure k with a limit below its eight assignments keeps its central carbon; a malformed
  molecule, whose records are all `Constitutional`.
- **Kind:** additive within the red stage. `[dep: S4b, S3c]`

#### S4d — Tests of `umol-graph-ir`

- **Module:** `tests/canonicalization.rs`, `tests/integrity.rs`, `tests/property/stereo`.
- **Changes:** the orbit tests of `tests/canonicalization.rs` drop the reading with isotopes
  ignored. `tests/integrity.rs` and the test of an out-of-range coset keep stating that symmetry
  does not fail. Tests that assert `Diastereotopic` for two different ligands assert `Distinct`.
- **Kind:** breaking (red→green for `umol-graph-ir`). `[dep: S4c]`

#### S4e — Operations of `umol-graph`

- **Module:** `ops::stereo` with `ops::stereo::recognize`, `ops::validate` with
  `ops::validate::stereo`, `ops::transform::stereo`, `ops::model`, the property tests under
  `tests/property`, the benchmarks.
- **Changes:** the perceiver, the validator and the enumerator build one `GraphSymmetry` per
  molecule and read the stereo views; the enumerator's cache per site is removed.
  `StereoSelection` selects by `stereogenic_units`, `StereoPerceiver::classify` returns the
  `StereogenicUnit` of each site and takes `max_assignments` in place of
  `para_stereo_max_iterations`, and the default of the adder is `Constitutional` and
  `Configurational`.
  `StereoModel` loses `para_stereo`, and `StereoValidateConfig`, `AddStereoConfig`,
  `RemoveStereoConfig` and `EnumerateStereoConfig` carry `max_assignments` in place of
  `para_stereo_max_iterations`. The benchmarks follow the changed API. The functions of
  `benches/transform/stereo.rs` that call the removed API are replaced by functions that measure
  building the `GraphSymmetry` and reading the stereo records.
- **Tests:** the operations on the S0a fixtures; the sequence add, enumerate, add of issue 43
  leaves the records unchanged; the results of add followed by enumerate, grouped by canonical
  equality, give the published counts of stereoisomers for the Razinger fixtures; the remover
  leaves the records of
  cis,trans-1,3,5-trimethylcyclohexane; the tests of the removed setting are removed; expected
  values that change with the definitions are shown old beside new.
- **Kind:** breaking (red→green for `umol-graph`). `[dep: S4c]`

#### S4f — `umol-py` follows the removals

- **Module:** `umol-py/src/model.rs`, `umol-py/src/model/stereo.rs`,
  `umol-py/src/stereo/perception.rs`, `umol-py/src/transform.rs`, their tests.
- **Changes:** the Python `StereoModel` loses `para_stereo`, the Python operation configs carry
  `max_assignments` in place of `para_stereo_max_iterations`, `StereogenicUnit` is added, and the
  Python `StereoSelection` selects by `stereogenic_units`.
- **Kind:** breaking (red→green for the workspace). `[dep: S4e]`

### S5 — Keyed readings and Python

#### S5a — `with_symmetry` and the keyed derived reading

- **Module:** `ir::view::constraints`, `ir::view::stereo`.
- **Adds:** `with_symmetry` on `StereoAtomConstraintsView` and `StereoBondConstraintsView`, and
  the constraints view reached from a stereo view carries its attached value.
- **Changes:** `derived` and `derived_complete` return the derived value as a literal form for
  the three keys and panic without the value attached; `evidenced`, `satisfies` and
  `is_compatible` follow.
- **Tests:** each keyed reading per key with and without the value; the fluxionality key and an
  entity without a kind stay vacuous; existing tests that state the vacuous derived side are
  shown old beside new.
- **Kind:** breaking (red→green within the subitem). `[dep: S4b]`

#### S5b — Python

- **Module:** `umol-py`.
- **Adds:** `GraphSymmetry`, `GraphSymmetryConfig`, `Molecule.graph_symmetry`, `with_symmetry` on
  the stereo views and the stereo constraints views, and the four derived readings on the stereo
  views. A derived reading without the value raises `ValueError`.
- **Tests:** the Python counterparts of S4a, S4b, S4c and S5a.
- **Kind:** additive, with the change of the keyed readings (green). `[dep: S5a, S4f]`

### S6 — Guides, final comparison, closeout

#### S6a — Guides

- **Module:** `docs/development/nomenclature.md`, `data-types.md`, `operations.md`,
  `python-api.md`, `property-tests.md`, `docs/chemistry/representation.md`, `README.md`.
- **Changes:** the classes and their definitions, the stereogenic units, `GraphSymmetry` and the
  three kinds of reading, the selection of the transformers, the modules, and the removal of the
  para-stereo setting and the coloring features. The entry of `nomenclature.md` on the
  description level calls the hierarchy private canonicalization machinery; it names
  `ir::entity` and both users, the readings and canonicalization.
  `[dep: S5b]`

#### S6b — Final comparison

- **Check:** the benchmarks of S0b against the baseline, with conclusions on the operations and
  on canonicalization. `[dep: S5b]`

#### S6c — Closeout

- **Check:** the workspace tests with all features, the Python tests, and clippy on the
  workspace. The results and the tables go into this document, and the status is set.
  `[dep: S6a, S6b]`

### Stage gate

A stage ends green: the tests of every crate pass, and after S2, S4 and S5 the Python tests pass.
S4 is the one stage that is red between its subitems.

### Critical path

S0 → S1 → S2 → S3 → S4 → S5 → S6. S1 and S2 do not depend on each other and can be done in
either order. No stage is deferrable: S5 completes the keyed readings that Python already
exposes.
