# Stereochemistry as molecule-level structure — framework and design space — 2026-05-28

Grounding note for adding stereo perception to umol-graph. This is analysis, not a
plan: it reviews the reference implementations, surveys the conceptual/literature
landscape (including non-standard, symmetry-explicit treatments), and lays out the
design space with trade-offs. No implementation choices are finalized here; the
genuine forks are collected in §6.

Supersedes nothing; extends [049-stereochemistry-encoding-2026-01-29.md](049-stereochemistry-encoding-2026-01-29.md),
which already argued for ports + a group-theoretic coset model + a molecule-level
artifact. Doc 049 was written against the older "GraphIR/TableIR/`sir_to_gir`"
naming; the semantic layer it calls *GraphIR* is today's `MoleculeAst` (umol-ast),
and the conversion site is `table_ir::raise` (umol-graph).

## 0. Current state in umol (May 2026)

- **`MoleculeAst` carries no stereo at all.** `AtomAst` has element / isotope /
  charge / implicit_hydrogens / lone_pairs / spin / constraints; `BondAst` has
  order / charge / spin / constraints. Stereo in the semantic layer is greenfield.
- **`MoleculeAst` already stores several molecule-level relations** —
  `aromatic_systems`, `multicenter_bonds`, `dative_bonds` (`VarRelationSet<T>`),
  `noncovalent_bonds` (`FixedRelationSet<T, 2>`), each keyed by a `RelationId`-backed
  id (`AromaticSystemId`, …) and referencing a set of `AtomId`s. Aromaticity in
  particular is a molecule-level relation over an atom set, **not** a per-atom flag.
  This is the existing precedent for the kind of top-level, element-referencing
  description the stereo positions below call for.
- **TableIR (parse container) already captures local stereo faithfully:**
  - `Atom.chirality: Option<Chirality>` — `Clockwise`/`CounterClockwise`/`Unspecified`
    plus arrangement-indexed `Tetrahedral{arr}`, `Allenal{arr}`, `SquarePlanar{arr}`,
    `TrigonalBipyramidal{arr}`, `Octahedral{arr}`.
  - `Bond.stereo: Option<BondStereo>` (`Cis`/`Trans`/`Either`), `Bond.wedge: Option<BondWedge>`.
  - `BicycloStereo` (CXSMILES THB/TLB/TEB).
  - `ExtendedAtom.attachment_order` / `ligand_order: Option<Vec<(u32, u8)>>` — the
    input-defined **port order** the §2.1 model needs is already preserved for CTFile.
  - molecule-level `StereoInterpretation` (`Absolute`/`Relative`) and
    `cx_data.stereo_groups: HashMap<u32, StereoSet { atoms, mode∈{Correlated,Independent} }>`.

So the open question is purely about the **semantic layer**: what stereo object lives
in `MoleculeAst`, and how `raise` builds it from the TableIR signals above.

### The four positions, restated as constraints

1. Stereo is a **molecule-level description referring to structural elements**, not a
   property of atoms (central chirality) or bonds (E/Z). Rationale: avoids
   context-dependence under renumbering / substructure extraction / reaction rewriting,
   and does not privilege central chirality over axial/planar/helical.
2. Internal configuration is **index-order-independent**, expressed against a canonical
   ordering of the inducing elements rather than SMILES `@`/`@@` traversal order.
3. The framework must leave room for non-tetrahedral centers, lone-pair centers,
   allene/cumulene, spiro, *trans*-cycloalkene, and non-central chirality — without
   model changes. SMILES-expressible cases (allene, cumulene, SP/TB/OH) are near-term.
4. Absolute / relative / racemic / meso must all be representable.

## 1. Where the field actually stores stereo (reference review)

`materials/codes` review. mx, octet, aefw_v1_0, PPP.jl are electronic-structure codes
with no stereo handling. The cheminformatics codes:

| Toolkit | Storage locus | Focus | Carriers | Configuration encoding | Enhanced stereo |
|---|---|---|---|---|---|
| **CDK** | molecule-level `IStereoElement<F,C>` list; **no per-atom chiral field** | atom *or* bond (`F`) | atoms *or* bonds (`C[]`), order significant | packed int: low byte order (`LEFT`/`RIGHT`=1/2 or larger), next byte class (`TH/CT/AL/CU/AT/SP/SPY/TBPY/OC/PBPY/HBPY8/9`, coordination number in 4th nibble) | high bytes of same int: `GRP_ABS/RAC/REL` + group number |
| **RDKit** | per-atom `ChiralType`, per-bond `BondStereo`+`stereoAtoms`; `StereoGroup` overlay; newer unified `Chirality::StereoInfo` | atom or bond (`centeredOn`) | ordered `controllingAtoms` | `StereoType` + `StereoDescriptor` + `permutation`; `specified∈{Unspecified,Specified,Unknown}` | `StereoGroup{type∈ABS/OR/AND, atoms, bonds, id}` |
| **OpenBabel** | molecule-level `OBStereoBase` generic data | center atom / bond | ordered `Refs` (stable atom **ids**, not indices); `ImplicitRef` sentinel for implicit H / lone pair | `Config` = winding (CW/ACW) + view (from/towards); equality via `NumInversions` parity | per-config `specified` flag |
| **Indigo** | molecule-level `MoleculeStereocenters` + `MoleculeAlleneStereo` + cis/trans side-tables | atom / bond | ordered "pyramid" neighbor list | type (`ABS/AND/OR/ANY`) + neighbor order | the type field |
| **LillyMol** | list of `Chiral_Centre` on the molecule | center atom | 4 ordered neighbors (incl. implicit-H / lone-pair slots) | clockwise/anticlockwise pyramid | — |
| **StereoMolGraph** (2025) | graph-level `dict[AtomId,AtomStereo]` / `dict[Bond,BondStereo]` | atom (`central_atom`) or bond | ordered `atoms` tuple (id), `None` slot = implicit/lone-pair | `parity∈{None,0,±1}` + explicit `PERMUTATION_GROUP` (rotation group as index perms) | (planned) |

**Finding.** Every mature general-chemistry toolkit *except RDKit* stores stereo as a
molecule-level collection keyed by a focus, with **ordered carriers** and a
configuration that is a permutation class. RDKit's per-atom chiral tag is a historical
artifact of SMILES round-tripping and is itself now overlaid by a unified descriptor
(`StereoInfo`) and a molecule-level group object (`StereoGroup`). Position 1 is not
contrarian — it is where the field's best designs already sit, and CDK has shipped it
as the *primary* store for years.

The convergent descriptor across all of them is the tuple

> **(focus, ordered carriers, configuration, specified-ness, group membership)**

and the only real disagreements are (a) what a "carrier" is (§4), (b) what the
configuration integer means (parity vs coset index, §2.1), and (c) what the ordering
reference is (input order vs CIP rank vs canonical rank, §2.2 / §6).

CDK's two extended classes show the focus is genuinely polymorphic:
`ExtendedTetrahedral` (allene) uses focus = central cumulated atom + 4 peripheral
carriers, collapsing focus+terminals into one pseudo-node so the **tetrahedral**
permutation group is reused; `Atropisomeric` uses focus = the **σ-bond** + 4 ortho
carriers. Same shape, different focus kind and symmetry group — the basis for
extensibility (position 3).

## 2. (i) Conceptual framework and literature precedents

### 2.1 The permutation / coset core (rigorous)

For a stereogenic element of coordination number *n*:

- Let `P = {p_1,…,p_n}` be the **ports** (attachment half-edges) on the carrier, and
  `Π = {1,…,n}` the idealized **positions** of the reference polytope.
- A **labeling** is a bijection `f: Π → P`; fixing any reference order makes it an
  element of the symmetric group `S_n`.
- Let `R ≤ S_n` be the **proper rotation group** of the idealized geometry acting on `Π`.

Two labelings represent the same configuration iff they differ by a proper rotation,
so a **configuration is a right coset in `S_n / R`**, and for distinct carriers the
count of configurations is `n! / |R|`. This reproduces the OpenSMILES arrangement
ranges exactly:

| Geometry | n | R | \|R\| | n!/\|R\| | OpenSMILES |
|---|---:|---|---:|---:|---|
| Tetrahedral | 4 | A₄ | 12 | 2 | `@`/`@@`, `@TH1..2` |
| Square planar | 4 | D₄ | 8 | 3 | `@SP1..3` |
| Trigonal bipyramidal | 5 | D₃ | 6 | 20 | `@TB1..20` |
| Octahedral | 6 | O | 24 | 30 | `@OH1..30` |

**Chirality of the local descriptor** is then a separate question from how many
configurations there are: the descriptor is chiral iff no odd permutation of positions
realizing the mirror is already in `R` (i.e. the full idealized point group, restricted
to position permutations, strictly contains `R`). Tetrahedral: mirror swaps two carriers,
not in A₄ → two distinct (enantiomeric) configurations → chiral. Square planar: the
mirror image is superimposable → the three SP arrangements are all achiral (they are
diastereomeric, not enantiomeric).

**StereoMolGraph is the concrete realization of exactly this**, published Dec 2025
(independent of doc 049). Its `Stereo` protocol stores an ordered `atoms` tuple, a
`parity ∈ {None, 0, ±1}`, and an explicit `PERMUTATION_GROUP` given as a list of index
permutations:

- `None` — orientation undefined (unspecified / "either").
- `0` — orientation defined but **achiral** (square planar, E/Z double bond): the
  geometry is its own mirror, so there is one orbit, not a ± pair.
- `±1` — chiral handedness (the two mirror cosets).

`canonical_form()` = the lexicographic minimum over the permutation-group orbit (and,
for chiral, the inversion orbit); `__eq__`/`__hash__` quotient by the group. That is
position 2 made operational: **index-order-independence is the orbit-min under the
geometry's symmetry group.** The concrete classes — `Tetrahedral` (inversion = swap two
carriers, group = A₄), `SquarePlanar` (parity 0), `TrigonalBipyramidal` (inversion =
swap axials, group = D₃), `Octahedral` (group = O), `PlanarBond` (E/Z, parity 0),
`AtropBond` (parity ±1), `NonRotatableBond` — are exactly the user's near-term + future
list, each an instance of one shape differing only in `(focus kind, symmetry group)`.

The three-valued-plus-achiral `parity` is the cleanest answer to position 4 I found in
any code: a meso center is a `±1` local descriptor that cancels under a **molecular
automorphism**, so molecular achirality is *derived* (automorphisms acting on the set of
local descriptors), never a per-center flag — see §2.2 for why this separation is the
right one.

### 2.2 Non-standard / symmetry-explicit precedents

The user asked specifically for treatments that make permutational or spatial symmetry
explicit, not necessarily implemented in current toolkits. These are the relevant lines:

- **Mislow & Siegel, "Stereoisomerism and Local Chirality," JACS 1984, 106, 3319.**
  The foundational separation umol should adopt vocabulary from: a site is **chirotopic**
  (local chirality — a spatial-symmetry property of the point) vs **stereogenic**
  (permuting its ligands yields a distinct stereoisomer — a permutation-group property
  of the constitution). The two are independent. This is precisely the split between
  "does the local descriptor have a handedness" (chirotopic, §2.1 chirality test) and
  "is the molecule chiral / are two molecules different" (needs molecular automorphisms).
  A clean stereo model keeps these in separate layers; conflating them is the root of
  the per-atom-flag pathology.

- **Ruch's algebraic theory of chirality / chirality functions** (Ruch & Schönhofer;
  Ruch, *The Permutation Group in Physics and Chemistry*, Lecture Notes in Chemistry).
  Classifies ligand permutations at a center via the symmetric group and Young tableaux
  / "chirality order," giving an algebraic measure of chirality rather than a binary
  R/S. Of interest for the user's tensor-spectroscopy angle (ROA, hyperpolarizabilities):
  chirality functions are the algebraic objects whose lowest non-vanishing terms control
  pseudoscalar observables. Likely overkill for a representation, but it is the rigorous
  source for "chirality is a representation-theoretic property of `S_n` acting on ligands."

- **Fujita's stereoisogram / RS-stereoisomeric group framework** (Fujita, many papers;
  e.g. "Stereogenicity/Astereogenicity as Global/Local Permutation-Group Symmetry,"
  *J. Math. Chem.*; the stereoisogram series in *J. Theor./Comput. chemistry*). The most
  complete symmetry treatment available: it integrates **three** groups in one structure
  — the **point group** (spatial symmetry), the **RS-permutation group** (ligand
  permutations, stereogenicity), and the **ligand-reflection group** (chirality) — and
  derives enantiomer / RS-diastereomer / holantimer relations as orbits. This is the
  literature realization of "use explicit permutational *and* spatial symmetry," and it
  maps onto umol cleanly: §2.1's `R` is the proper-rotation part, the chirality test is
  the ligand-reflection part, and molecular automorphisms are the global counterpart.
  Worth reading before fixing the configuration encoding, because it tells you which
  quotients are meaningful.

- **InChI's layered stereo** (IUPAC; InChI Technical Manual). A working,
  canonical, *molecule-level* stereo description that is **not** CIP-based: tetrahedral
  parity (`/t`) is the sign of the oriented tetrahedron volume w.r.t. InChI's own
  **canonical atom numbering**; double bonds in `/b`; `/m0`,`/m1` flags the enantiomer;
  `/s1`,`/s2`,`/s3` flags absolute / relative / racemic. Two lessons: (a)
  index-independence does **not** require CIP — a graph-canonical numbering is sufficient
  and far cheaper (relevant to position 2 and the §6 ordering fork); (b)
  absolute/relative/racemic as a small molecule-level enum (`/s`) plus a global parity
  flip (`/m`) is a proven minimal encoding for position 4.

- **OpenSMILES non-tetrahedral arrangement indices** are themselves a coset enumeration
  (`@SPn`, `@TBn`, `@OHn`) — the spec already speaks the §2.1 language; a coset-index
  encoding round-trips to/from SMILES with no translation layer.

### 2.3 Connection to the umol-msym symmetry stack

umol already has a point-group / symmetry module (umol-msym; see the symmetry roadmap in
doc 074 and the SALC work in doc 072). The idealized-geometry groups `R` of §2.1 are
small fixed point groups (A₄, D₄, D₃, O, …); they could be hard-coded per stereo class
(as StereoMolGraph and CDK do) **or** sourced from the same machinery that produces
point groups for the spectroscopy work, which would make "add a new coordination
geometry" a data change rather than a code change. This is a genuine architectural
synergy given the user's research priorities (chirality is the symmetry property ROA /
optical activity probe), and it is one of the forks in §6.

## 3. (ii) Algorithms and open-source implementations (web search)

### Configuration assignment / CIP labeling

- **Hanson, Musacchio, Mayfield, Vainio, Yerin, Redkin, "Algorithmic Analysis of
  Cahn–Ingold–Prelog Rules of Stereochemistry," *J. Chem. Inf. Model.* 2018, 58,
  1755–1765.** The reference machine-implementable specification of CIP, including
  corrections to Rules 1b and 2. Any rigorous R/S layer should follow this, not the
  textbook prose.
- **`centres`** (John Mayfield) — the standalone reference CIP library the paper
  describes; **RDKit `rdCIPLabeler`** is a C++ port of it; **Jmol** has an independent
  rigorous CIP implementation by Hanson. These are the three trustworthy CIP codes;
  legacy RDKit `AssignStereochemistry` and CDK's older `CIPTool` are approximate.
- CIP is **derived presentation**, not storage: it is convention-heavy, expensive, and
  orthogonal to the stored state (doc 049 §8.3). It layers on top of whatever
  ordering-independent state umol stores.

### Perception (which elements are stereogenic; from 0D/2D/3D)

- **RDKit** `Chirality::findPotentialStereo` → `std::vector<StereoInfo>`, with the
  `detail::isAtomPotential*` / `isBondPotentialStereoBond` predicates; toggle flags
  `useLegacyStereoPerception`, `setAllowNontetrahedralChirality`.
- **CDK** `StereoElementFactory` (from 2D coords + wedges, or 3D) and `Stereocenters`
  (true/para/potential classification).
- **OpenBabel** `StereoFrom0D` / `StereoFrom2D` / `StereoFrom3D` + `OBStereoFacade`.
- **LillyMol** `is_actually_chiral` — automorphism-based test of whether a candidate
  center is genuinely stereogenic (distinguishes from constitutionally symmetric).
- Perception is the hard part (doc 049 §8.2); detection scope is a §6 fork.

### Canonicalization, automorphism, meso / diastereomer detection

- **Schneider, Sayle, Landrum, "Get Your Atoms in Order," *J. Chem. Inf. Model.* 2015,
  55, 2111** — canonical ranking including stereo (the basis of RDKit canonical SMILES).
- **InChI** canonical labeling + stereo normalization (above).
- **nauty / Traces** (McKay) — graph automorphism groups; the principled route to meso
  detection and to collapsing equivalent configurations (doc 049 §3.3). umol already has
  an `automorphism` module in umol-ast — the molecular-symmetry half of §2.1 is in part
  available.

### Enhanced (mixture / relative) stereo

- **V3000 MOLfile** `MDLV30/STEABS/STERAC/STEREL` and **CXSMILES** `a`/`&n`/`o n` groups
  → RDKit `StereoGroup{ABSOLUTE/AND/OR}`, CDK `GRP_ABS/RAC/REL` bits. umol's TableIR
  `cx_data.stereo_groups` (`Correlated`/`Independent`) and `StereoInterpretation` already
  capture the inputs; the semantic layer needs the matching per-element group concept.

### Group-theoretic graph representation (the closest precedent)

- **StereoMolGraph** — Papusha et al., ChemRxiv 2025 (DOI `10.26434/chemrxiv-2025-0g4wn`),
  GitHub `maxim-papusha/StereoMolGraph`. Permutation-invariant local stereodescriptors
  grounded in group theory; supports non-tetrahedral centers and *changing*
  stereochemistry in reactions / transition states (condensed reaction graphs). Closest
  existing realization of §2.1 and of doc 049; the `stereodescriptors.py` design is worth
  reading in full before fixing umol's encoding. Relevant to the user's reaction-network
  work (memory: combinatorial reaction networks) because it handles stereo *change*.

## 4. (iii) Substituent / carrier representation

What identifies the elements that induce the stereochemistry? Four options:

| Option | Used by | Resolves multi-edge? | Lone pair / implicit H | Survives edits | Cost |
|---|---|---|---|---|---|
| **Neighbor atoms** (atom ids) | RDKit `controllingAtoms`, OB `Refs`, SMG `atoms`, Indigo, LillyMol | **No** — two bonds to the same neighbor are indistinguishable | sentinel (`ImplicitRef` / `None` / placeholder terminal atom) | atom ids stable under attribute edits | minimal |
| **Connecting bonds** (edge ids) | CDK `IStereoElement<F, IBond>` for CT / atropisomer | Yes | a bond to nothing can't exist → still needs a sentinel for lone pair / implicit H | edge ids stable | minimal |
| **Ports / half-edges** `(focus, incidence)` | doc 049; not first-class in any reviewed toolkit | **Yes** — each endpoint is a distinct port | a lone pair or implicit H is just a port with no atom on the far side; no sentinel | needs port identity preserved across edits | one indirection |
| **Full substituent fragments** | nobody (storage); what CIP/canonical ranking *computes over* | Yes | n/a | duplicates the graph; breaks under any edit | heavy |

**Why this matters specifically for umol.** umol's scope is general chemistry, not
drug-like (memory: multicenter / dative bonds, organometallics, mixed-valence are core).
In that scope "the neighbor atom" is genuinely ambiguous or absent:

- multi-hapto / multicenter bonds (already a `VarRelationSet` in `MoleculeAst`): a metal
  center's "ligand" may be an η²/η⁵ system, not one atom;
- dative bonds and multiple bonds between the same pair (multi-edges);
- lone-pair stereocenters (position 3) — the inducing element is **not an atom**;
- aromatic-system carriers.

The neighbor-atom option papers over all of these with sentinels (the `ImplicitRef` /
`None` slots every atom-based code carries). The **port** option is the one that makes
them first-class: a carrier is a half-edge at the focus, an atom / bond / fragment is a
*derived view* of where that half-edge lands, and a lone pair is simply a port whose far
side is empty. This is the carrier identity consistent with constraints 1 + 3; its price
is that the port table must keep stable identity across structural edits (which already
go through `MoleculeBuilder`, so there is a single place to maintain it).

Full **fragments** are never the stored identity — they are needed only transiently, by
CIP and by canonical ranking, to order the carriers, and are derivable from topology on
demand. Storing them would duplicate the graph and invalidate on every edit.

This is the one place where I think the constraints select an answer rather than leaving
a balanced choice: atoms vs ports is a real trade (simplicity vs the multicenter / lone-
pair cases), but fragments-as-identity is dominated, and bonds-as-identity is a strict
subset of ports. The §6 fork is therefore narrowed to **atom-ref vs port**.

## 5. A candidate umol model (for discussion, not decided)

Synthesizing §1–§4 against the existing `MoleculeAst` architecture, the shape that fits
is a **molecule-level stereo relation, parallel to `aromatic_systems`**:

- A `StereoElement` per stereogenic feature, holding:
  - **focus** — an atom, a bond, or (later) a higher carrier; polymorphic as in CDK.
  - **ordered carriers** — ports (§4 fork), as many as the coordination number.
  - **geometry / symmetry group** — a small tag (`Tetrahedral`, `SquarePlanar`,
    `TrigonalBipyramidal`, `Octahedral`, `Allene`, `Cumulene`, `Atropisomer`, `LonePair…`)
    carrying its proper-rotation group `R` and its chirality test.
  - **configuration** — an orbit / coset in `S_n/R`, encoded either as `parity∈{None,0,±1}`
    + group (StereoMolGraph) or as a coset index (OpenSMILES/CDK); fork in §6.
  - **specified-ness** — three-valued (`Specified` / `Unspecified` / `Unknown`), as RDKit
    `StereoSpecified` and OB; distinct from "achiral" (`parity 0`).
  - **enhanced-stereo group membership** — `Absolute` / `Racemic(n)` / `Relative(n)`,
    fed by TableIR `StereoInterpretation` + `cx_data.stereo_groups`.
- Stored in a new relation set keyed by a `StereoElementId` (`RelationId`-backed), exactly
  like the other relations — which gives remap / equivariance under renumbering for free,
  because the relation sets already participate in `remap`.
- **Absolute / relative / racemic / meso (position 4):** per-element group membership for
  absolute/relative/racemic; **meso and molecular chirality are derived** from the
  molecule's automorphism group acting on the set of local descriptors (Mislow–Siegel
  separation, §2.2), not stored as a flag.
- **Ordering-independence (position 2):** carriers are ordered by a canonical reference;
  the stored configuration is the orbit-min under `R` w.r.t. that reference. The reference
  is either CIP rank or a graph-canonical rank — §6 fork; InChI shows canonical rank
  suffices and is cheaper, with CIP R/S as a derived presentation layer.
- **Conversion (`table_ir::raise`):** for each TableIR `chirality{arr}` / `BondStereo` /
  `wedge`, build the port set from the graph, read the input port order from
  `attachment_order` / `ligand_order` / parse spans, compute the permutation to canonical
  order, combine with the parsed arrangement to select the coset, and emit a
  `StereoElement`; molecule-level `StereoInterpretation` / `stereo_groups` populate the
  group memberships. This is doc 049 §6 in current terms.
- **Extensibility (position 3):** allene/cumulene = focus-as-bond(s) collapsing the
  cumulated axis to one pseudo-node reusing the tetrahedral/CT group (CDK
  `ExtendedTetrahedral`); atropisomer / *trans*-cycloalkene = focus-as-bond with the
  atrop group; spiro / planar / helical = new geometry tags with their own `R`; lone-pair
  center = a port with empty far side. None require touching the `StereoElement` shape.

## 6. Open decisions to weigh

Genuine forks; each changes downstream code materially. No pick here.

1. **Carrier identity — atom-ref vs port (§4).** Atom-ref: minimal, matches every
   toolkit, needs implicit/lone-pair sentinels and cannot disambiguate multi-edges.
   Port: first-class implicit/lone-pair/multicenter, costs a port table with stable
   identity across builder edits. (Fragments and bonds-only are dominated — see §4.)
2. **Configuration encoding — parity+explicit group (StereoMolGraph) vs packed coset
   index (OpenSMILES/CDK).** Parity+group: uniform across geometries, `canonical_form`
   and equality fall out of the group, reaction/TS-friendly. Coset index: zero-translation
   round-trip to SMILES `@TBn`/`@OHn`, compact, but per-geometry index tables and
   bespoke equality.
3. **Ordering reference — CIP rank vs graph-canonical rank.** CIP: chemistry-meaningful,
   R/S immediate, but heavy and convention-laden (needs a `centres`-class implementation).
   Canonical rank (InChI-style): cheap, self-consistent, index-independent; R/S becomes a
   separate derived layer. These are not exclusive — canonical rank for storage, CIP for
   presentation — but which is *authoritative* is a real choice.
4. **Where the geometry symmetry groups live — hard-coded per class vs sourced from
   umol-msym (§2.3).** Hard-coded: simple, matches StereoMolGraph/CDK. Sourced: "new
   geometry = data, not code," and reuses the spectroscopy point-group machinery.
5. **Detection scope, now vs later.** Which elements `raise`/perception emits initially
   (tetrahedral + E/Z + the SMILES-native SP/TB/OH/allene), and which are
   representable-but-not-yet-perceived (lone pair, atropisomer, spiro, helical).
6. **Storage — a dedicated `StereoElement` relation set vs reusing molecule-level
   `Constraints`.** A relation set matches `aromatic_systems` and gets remap for free; the
   `Constraints` `Vec<Constraint>` is for predicate/combinator constraints and is a poorer
   fit, but avoids a new entity type.

## Sources

- [Algorithmic Analysis of CIP Rules (Hanson et al., JCIM 2018)](https://pubs.acs.org/doi/abs/10.1021/acs.jcim.8b00324) · [RDKit rdCIPLabeler docs](https://www.rdkit.org/docs/source/rdkit.Chem.rdCIPLabeler.html) · [RDKit CIPLabeler PR #3234](https://github.com/rdkit/rdkit/pull/3234)
- [Mislow & Siegel, "Stereoisomerism and Local Chirality," JACS 1984 (PDF)](https://moodle2.units.it/pluginfile.php/344286/mod_resource/content/1/JACS_1984_Mislow_Siegel.pdf)
- [Fujita, Stereogenicity/Astereogenicity as Global/Local Permutation-Group Symmetry (J. Math. Chem.)](https://link.springer.com/article/10.1023/A:1023251932146) · [Stereoisograms for reorganizing stereochemistry (ScienceDirect)](https://www.sciencedirect.com/science/article/abs/pii/S0957416614003334)
- [Applications of the Permutation Group in Dynamic Stereochemistry (Springer)](https://link.springer.com/chapter/10.1007/978-3-642-93124-6_3)
- [InChI, the IUPAC International Chemical Identifier (J. Cheminform. 2015)](https://link.springer.com/article/10.1186/s13321-015-0068-4) · [InChI Technical Manual (PDF)](https://www.inchi-trust.org/download/104/InChI_TechMan.pdf)
- [StereoMolGraph (ChemRxiv 2025, PDF)](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2025-0g4wn) · [StereoMolGraph (GitHub)](https://github.com/maxim-papusha/StereoMolGraph)
- [Stereochemistry-aware string-based molecular generation (ChemRxiv 2024)](https://chemrxiv.org/doi/pdf/10.26434/chemrxiv-2024-tkjr1)
- [Open Babel Stereochemistry docs](https://open-babel.readthedocs.io/en/latest/Stereochemistry/stereo.html)
