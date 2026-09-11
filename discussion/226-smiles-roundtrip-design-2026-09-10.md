# 226 — SMILES roundtrip design

Status: In Progress
Date: 2026-09-10
Relates: [155](155-smiles-io-and-resolve-configuration-2026-07-19.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[166](166-molecule-ops-2026-07-27.md),
[170](170-reaction-smiles-python-2026-07-28.md),
[224](224-smiles-ring-closure-frame-2026-09-08.md),
[data type contracts](../docs/development/data-types.md)

Correction (2026-09-11): the projection contract and S4 plan below replace the runtime
reconstruction requirement. The S4a and S4b corrections are complete. Their earlier
implementation and timing records are retained as historical evidence. S4a0a–S4a0c remain complete.

## Scope

Design the reverse molecular-format path: resolver-owned projection within graph IR, followed by
Convey conversion to TableIR and construction of Smiles or ReactionSmiles for rendering.
TableIR is the shared SMILES/CTfile
boundary representation. Decide what it must carry for aromatic systems and stereo, and whether
shared storage suffices or explicit format separation is needed. The settled design and staged
implementation plan are recorded below; S0 inventory, semantic fixtures, and baselines are complete.

Responsibilities, preservation guarantees, boundary types, and core algorithms are settled.
The stages below sequence implementation and verification within that agreed scope.
CTfile output requirements inform the shared representation; implementing a CTfile writer is a
separate scope. CXSMILES and broader parser work remain with doc 153.

## Existing foundations

- Doc 155 §§1.1–1.3 defines ordered, syntax-normal parse/render roundtrips, distinguishes them
  from graph canonicalization, and places representability checks at graph-to-SMILES conversion.
  Its completed configuration work did not implement a renderer.
- Smiles and ReactionSmiles privately contain TableIR values. Their current public surface parses
  text/bytes and exposes borrowed or consuming TableIR access. Neither has rendering or public
  construction from TableIR: see [molecule.rs](../umol-io/src/smiles/molecule.rs) and
  [reaction.rs](../umol-io/src/smiles/reaction.rs).
- TableIR Molecule and Reaction are open field carriers. Atom has aromatic and chirality fields;
  Bond has order, direction, stereo, and wedge fields. Molecule additionally has explicit
  tetrahedral stereo_atoms, positions, and configuration_scope. There is no parallel explicit
  stereo-bond frame collection: see [table_ir](../umol-io/src/table_ir.rs).
- Doc 224 introduced explicit tetrahedral ligand frames while retaining opening-order bond storage.
  TableIR raise converts these frames directly to graph-IR stereo atoms and skips the wedge-derived
  tetrahedral constraint for a site with an explicit frame. Atom.chirality is not the operative
  tetrahedral frame. Cis/trans raising still interprets directions, annotations, and geometry:
  see [raise.rs](../umol-io/src/table_ir/raise.rs).
- Resolver and the aromaticity/stereo resolvers currently have no project operation.
  Both resolve and the proposed project operate exclusively on graph IR; resolvers know nothing
  about TableIR. Convey owns conversion from projected graph IR to TableIR:
  see [resolve.rs](../umol-graph/src/ops/resolve.rs).
- Graph ingestion composes boundary interpretation with resolution. ReactionSmiles interpretation
  rejects nonempty agents and repeated members of a map class on either side before resolving the
  sides: see [ingest.rs](../umol-graph/src/ingest.rs). Boundary parsing retains those representations.

## Two operations and their boundaries

```text
graph-IR Molecule
    -- resolver projection --> graph-IR Molecule
    -- Convey: GraphIR-to-TableIR conversion + boundary construction --> Smiles
    -- rendering under IO configuration --> text

text -- parse --> Smiles
    -- Interpret: TableIR raise --> graph-IR Molecule
    -- resolver resolution --> graph-IR Molecule
```

Boundary construction establishes that the projected table can be expressed as SMILES. It is a
contract boundary within the second operation, not a proposed extra processing framework.

Resolver projection transforms graph IR into graph IR with resolve's mutation semantics.
Convey clones the caller's graph IR for this working transformation, then expresses the projected
atom/bond information and supported overlays
in TableIR terms, with explicit frame transport. Resolver projection has no knowledge of TableIR
or the downstream output format and cannot choose format fields. Aromaticity and stereo work
within graph IR remains with the owning resolver machinery; boundary encoding belongs to Convey.
Projection does not rerun resolution implicitly, select a surviving interpretation, or use recovery
policies to discard inconvenient information. The exact model dependence must be stated: a model
may be needed to establish that notation reinterprets to the supplied state, rather than merely to
copy fields.

Formatting owns traversal, components, branches, ring labels, bracket/implicit-H spelling, and
stereo markers relative to the emitted traversal. The writer consumes boundary information without
aromaticity perception or chemical resolution. SmilesIoConfig retains doc 155's bidirectional role;
chemistry model choices remain separate. Algorithm selectors, if exposed at the IO layer, are explicit.

## General export representability rule

Export must fail for information that the supported output notation cannot
represent. Nonzero charge or non-singlet spin carried by bonds or aromatic
systems cannot be silently transferred to atoms or discarded; explicit
localization must precede export. Already atom-localized attributes remain
subject to the output format's representability requirements. Undetermined
values do not authorize assuming zero charge or singlet spin.

Unsupported bond types likewise cause export failure. Do not replace them with
supported bond types, omit them, or otherwise approximate the input. Any
required localization or other transformation is an explicit operation outside
projection, boundary conversion, and formatting. This is a general rule, not
an aromaticity-specific exception.

## Shared TableIR and apparently duplicated information

Two distinct questions must not be conflated: whether TableIR can carry both formats' input
notations, and whether projection must populate both output notations simultaneously.

### Aromatic systems

Atom.aromatic and BondOrder::Aromatic already express different incidences. Aromatic endpoint atoms
do not by themselves determine whether their connecting bond is aromatic. The projected table must
distinguish an aromatic-system bond from a localized connection between aromatic atoms. The writer
then chooses lowercase symbols and explicit or omitted bond tokens from those values.

The Convey boundary-encoding rule is settled: mark participating atoms aromatic and encode
the existing systems' bonds as aromatic bonds. Preserve localized connections
between separate systems. Do not perceive new systems or choose a Kekulé form
as an implicit output transformation.

An aromatic system with nonzero system-level charge or non-singlet system-level
spin causes export failure. These quantities must be explicitly localized before
export; projection does not distribute them over atoms. This restriction concerns
the aromatic-system attributes, not already atom-localized charges. It does not
authorize treating undetermined values as zero or singlet.

This needs coordinated atom and bond information, but not necessarily separate SMILES and CTfile
copies of an aromatic system. A CTfile writer can consume aromatic bond orders; output requiring a
localized Kekulé representation introduces a separate selection question outside this writer's
scope. Projection uses the existing systems and their electron contributions directly; it does
not search for a reconstructible assignment.

### Stereo atoms and bonds

The proposed direction is to retain explicit configurations in stated ligand frames as the common
output information. For tetrahedral sites this builds on StereoAtom from doc 224. A SMILES writer
transports that frame into its actual emitted neighbor order before choosing a marker. It cannot
copy the old chirality token after changing traversal. CTfile boundary conversion must determine
which supported notation expresses the configuration with the coordinates actually present, or
without coordinates when they are absent. Neither projection, boundary conversion, nor formatting
generates coordinates.

Stereo bonds use the explicit StereoBond/BondConfiguration/BondRelation design below. Raw
Bond.direction and Bond.stereo are consumed at format reading and retired from the shared table;
wedges remain separate. The alternatives considered were coordinated duplicate encodings and
separate format representations. Neither is needed for the settled ordinary stereo scope.
A format-independent resolver cannot choose downstream format fields; Convey owns that encoding.

### Retained stereo-atom frames

Retain StereoAtom and StereoLigand::LonePair from doc 224. The proposed removal was withdrawn
following the stereo/valence experiments: bracket H omission means H0, so a three-incidence marked
frame does not authorize inferred hydrogen. Its lone-pair participant is an assertion whose
chemical consistency is checked after valence resolution. Preserve source encounter order and
root/incoming-neighbor placement of virtual participants. The separate specification staging list
below records the clarification; the carbon shell-capacity defect remains deferred to doc 166.
Do not replace the current four-participant frame with a new placeholder or three-ligand type.

### Explicit stereo-bond frames

The complete GraphIR frame motivating the smaller TableIR record is:

```text
site bond: u--v
endpoint u: [a, b]
endpoint v: [c, d]
configuration: a and c on the same side, or on opposite sides
```

The two ordered ligand pairs are the frame. The configuration relates their first members; their
second members occupy the complementary positions. These are frame-relative descriptors, not
CIP E/Z labels or stored SMILES slash directions. No coordinates, drawing orientation, traversal,
or priority ranking are needed to interpret the record. The smaller TableIR layout below is settled.

This follows the existing graph-IR stereo-bond contract: a bond site and two consecutive
two-ligand endpoint blocks. The frame actions are:

| Change to the frame | Change to its same/opposite descriptor |
| --- | --- |
| Exchange a and b only | Flip |
| Exchange c and d only | Flip |
| Exchange both pairs internally | Preserve |
| Exchange the complete endpoint blocks and their endpoints | Preserve |
| Move one ligand to the other endpoint | Not a frame change of the same site incidence |

For example, endpoint blocks `[F, H]; [Cl, H]` with F and Cl opposite describe the same configuration
as `[H, F]; [Cl, H]` with the first H and Cl on the same side. Actual configuration conversion should use
the existing umol-perm frame action; the descriptor's correspondence to numeric cosets must be
pinned explicitly rather than inferred from enum or integer ordering.

#### Site, ligand identity, and integrity

The record needs one site identity, not independently supplied bond and endpoint identities that
can disagree. A table bond index can provide it, paralleling StereoAtom's atom index. TableIR's
AtomPair orders the endpoints; define the first ligand block against its first endpoint and the
second against its second endpoint. Projection must transport the source frame into that endpoint
order. Atom renumbering that reverses the ordered pair exchanges the complete blocks; reordering
the bond table transports the site index. Neither operation reconstructs the frame from neighbors.

A complete graph-IR frame distinguishes actual atoms, implicit hydrogens, and lone pairs, with
virtual ligands anchored at their endpoint. This does not establish that the external formats
represent implicit lone-pair ligands. The OpenSMILES double-bond convention relates actual bonds;
CTfile lists LP as an atom symbol. Explicit LP pseudoatoms are a separate topic and are outside
this design. Doc 224's tetrahedral LonePair slot is an umol boundary convention, not evidence for
an implicit lone-pair participant in the double-bond notation.

This shape covers an ordinary bond site. Extended cumulene stereo relates terminal ligand frames
across an axis, which the current local bond-site contract does not express. Do not pretend that
adding this record implements that coverage or silently attach remote ligands to a local bond.

#### A smaller explicit boundary frame

The complete graph-IR frame explains the frame action, but need not dictate the TableIR record.
For an ordinary double bond with an actual substituent at each endpoint, the boundary can instead
state:

```text
reference atom a -- site endpoint u == site endpoint v -- reference atom c
configuration: a and c on the same side, or on opposite sides
```

This is a four-atom reference frame, with two actual reference substituents. It does not assert
what occupies either complementary position. A site bond index plus one reference atom for each
of its ordered endpoints supplies the frame without duplicate endpoint storage. Both references
must exist and have the required incidence, and each site has one operative record. These are
representation requirements, not a stereogenicity test or a chemistry-model judgment.

Swapping the complete endpoint/reference pairs preserves the descriptor. Selecting the other
actual substituent at one endpoint flips it. Selecting both alternatives preserves it. No CIP
ranking, coordinates, or SMILES traversal is needed to interpret this relative configuration.
The boundary record and reference-selection policy below are settled design;
implementation remains future work.

The reference-selection policy for definite configurations is settled: at each endpoint select the actual
substituent with the lowest TableIR atom index, excluding the opposite site
endpoint. Apply this rule in SMILES reading, CTfile interpretation, and graph-IR
projection, independently of which source bonds carry direction markers.
Interpret and check all supplied directional markers before expressing the
configuration in the selected frame. Both substituent bonds at an endpoint may
be marked; their consistency is judged relative to that endpoint, accounting
for traversal direction, not by comparing the slash characters literally.
Switching the reference at one endpoint flips the descriptor; switching both
preserves it. This is deterministic within the table's numbering, not chemical
canonicalization. A writer may choose other bonds after transporting the
descriptor to its output frame.

The agreed Rust shape preserves both framed and site-only assertions:

```rust
pub struct StereoBond {
    pub bond: u32,
    pub configuration: BondConfiguration,
}

pub enum BondConfiguration {
    Either,
    Framed {
        references: [u32; 2],
        relation: BondRelation,
    },
}

pub enum BondRelation {
    SameSide,
    OppositeSide,
}
```

Here `references[0]` belongs to the referenced bond's AtomPair::first endpoint,
and `references[1]` to its second endpoint. All indices belong to the owning
TableIR molecule. The record is an open carrier, like StereoAtom; no constructor
or independently supplied endpoint pair is proposed. Checked consumers that
require a meaningful frame must establish index validity, substituent incidence,
and consistency with the site and other operative assertions before using it.
Malformed references or conflicting assertions are conversion failures, not
indexing panics or invitations to repair the input. Chemical validity remains
with molecular interpretation. The minimum-index rule is a producer policy,
not a reason to reject an otherwise meaningful independently supplied frame.

The agreed additional public names are StereoBond, BondConfiguration, and BondRelation;
collection placement, exact conversion/error APIs, and Python exposure remain
to be settled. BondConfiguration::Either preserves a site-only assertion without
references, using the CTfile name. Framed carries the selected references with
a definite SameSide or OppositeSide relation. No record means no stereo
assertion; it is distinct from Either. No source changes are authorized
by this design agreement.

Projection from a graph-IR stereo bond chooses actual reference ligands from the endpoint blocks
and expresses their relative configuration. The writer can use different neighboring bonds by
transporting that descriptor. A definite source frame lacking an actual reference on either side
does not fit the Framed contract; its output representability requires a separate decision rather
than an invented reference atom or a downgrade to Either. Projection must still account for
source distinctions the boundary cannot carry.

The reverse path retains a useful separation. Format reading determines the reference relation;
molecular interpretation completes any graph-IR ligand identities. In current code,
[raise/utils.rs](../umol-io/src/table_ir/raise/utils.rs) uses an endpoint-anchored Virtual placeholder,
and StereoPerception::bond_side_ligands in [stereo.rs](../umol-graph/src/ops/stereo.rs) later selects
H or lone pair from molecular values. That is umol interpretation machinery, not evidence that
SMILES or CTfile supplies an unidentified implicit ligand. The proposed boundary frame does not
need a new unidentified-ligand variant merely to duplicate that machinery.

The selected smaller frame avoids copying all four graph-IR ligands. It carries
explicit stereo independently of raw slash marks or drawing coordinates, while retaining molecular
ligand completion in raise/resolution. Unlike a fully identified four-ligand frame, it is not in
general a direct graph-IR stereo-entity constructor. The symmetry with StereoAtom should be shared
frame ownership and explicit configuration, not necessarily identical participant storage.

#### Settled raise and convey contract

TableIR's two-reference frame need not identify the complementary molecular
ligands before resolution. Raise converts Framed to the existing bond `#C`
assertion, transporting SameSide/OppositeSide into its canonical endpoint and
neighbor frame. Orientation is established without assigning chemical identity
to a complementary virtual position. Either becomes `#C` with an undetermined
configuration; absence of a record supplies no assertion. Resolution retains
its existing responsibility to derive the complete molecular ligand frame from
the resolved atom states and apply the assertion.

Convey conversion of definite graph-IR stereo selects the lowest-indexed actual
substituent at each endpoint in the output TableIR numbering and transports the
source configuration to those references. An undetermined configuration becomes
Either without references. A definite configuration that cannot be expressed
in this boundary frame is a representability failure, never a downgrade to Either.

No semantic choice remains open here for ordinary cis/trans configurations.
The exact frame permutation and coset correspondence must be derived and
verified against the existing stereo algebra before implementation.

#### Settled replacement of bond-stereo fields

Consume and check Bond.direction markers during format reading, derive the
stereo-bond records under the agreed equivalence rules, and remove the stored
direction field. Replace Bond.stereo with StereoBond records, translating its
supported source annotations, including Either and the existing CXSMILES
Cis/Trans producers, into the explicit boundary meaning. There is one operative
representation of interpreted double-bond configuration.

Wedges are separate from this work. Axial and other extended stereo, including
their atom-chirality annotations, require no redesign for this task and remain
outside its scope. No broader atom-stereo representation decision is pending
under this item. Point 2 is settled; exact producer migration is implementation
work, not a reason to reopen the representation choice.

#### Absence, unknown configuration, and partial notation

Three meanings must remain distinct:

- no stereo assertion at the site;
- an asserted cis/trans site whose configuration is unknown;
- a definite configuration in a stated frame.

The current reader already distinguishes these: absent or nondetermining geometry produces no
cis/trans constraint, while CTfile Either produces a cis/trans assertion with an undetermined
coset. The latter does not require reference substituents in the current raise path. The CTfile
bond record likewise does not require them, so an assertion cannot be discarded merely because
references are unavailable. Preserve this assertion as BondConfiguration::Either,
whether or not reference atoms could be selected. Projection of an undetermined
graph-IR stereo configuration likewise produces Either without references.

A wavy single bond adjacent to a double bond can also assert an unknown
double-bond configuration. Its interpreted assertion is Either, not a framed
relation. Which source bond was drawn wavy is a separate annotation concern;
minimum-index references would not preserve that placement. There is no
BondRelation::Either variant.

No frame therefore must not mean an explicit NotStereo assertion. Likewise, omitting output slash
marks cannot automatically stand for every explicitly unknown source assertion. The destination
boundary must establish whether those meanings coincide in its supported contract.

Partial source directions also need classification before removal from TableIR. A marker on only
one side need not determine a configuration. Multiple markers can share a neighboring single bond.
Retain their source meaning where required by the boundary roundtrip; do not promote an incomplete
reading to a definite frame. This work must reconcile the repository specification with current
raise behavior rather than changing parser acceptance incidentally.

#### Direction-marker examples

The following examples separate marker interpretation from the proposed stored
frames. A scratch probe ran all 29 inputs in
`scratch/stereo-valence-scan/directions.smi` through Smiles::parse and TableIR
raise on 2026-09-10; exact results are in `directions.jsonl`. This probes the
current reader, not an implemented frame producer or writer. No production
code or permanent tests changed.

In these tables, S and O mean SameSide and OppositeSide relative to the agreed
minimum-index references. A dash means no configuration assertion, not Either.
For chains with multiple double bonds, entries follow their order in the text.
The probe's numeric cosets were interpreted against the current sorted-neighbor
frame; the proposed enum's numeric representation is not thereby specified.

| Input | Interpreted relation | Point to preserve |
| --- | --- | --- |
| `FC=CF` | — | No markers, no assertion. |
| `F/C=CF` | — | One marked endpoint does not determine a relation. |
| `F\C=CF` | — | Reversing that lone marker still determines no relation. |
| `FC=C/F` | — | The same applies at the other endpoint. |
| `F/C=C/F` | O | A complete pair determines a frame. |
| `F\C=C\F` | O | Reversing both markers preserves the configuration. |
| `F/C=C\F` | S | Reversing one marker changes the configuration. |
| `C(/F)=C/F` | S | A branch marker is read from its own endpoint viewpoint. |
| `F/C(Cl)=C(/Br)I` | O | References are F and Br; complementary ligands need no markers. |
| `F/C(/Cl)=C(/Br)\I` | O | Marking all four ligands adds no configuration information. |
| `C(/F)(\Cl)=C(/Br)\I` | S | Both slash glyphs at one endpoint are permitted when consistent. |
| `F/C(\Cl)=C/Br` | Raise: CisTransConflict at atom 1 | Conflicting endpoint markers must not be hidden by selecting one reference. |
| `F/C(\Cl)=CBr` | Raise: CisTransConflict at atom 1 | Even partial notation must be checked before discarding markers. |
| `C/C` | Raise: DanglingBondDirection at bond 0 | Current raise rejects a directional bond without an eligible neighboring double bond. |

The endpoint-viewpoint rule and redundant complete markings are described in
[OpenSMILES §3.8.3](https://opensmiles.org/opensmiles.html). The table records
umol's current behavior for partial and rejected inputs; it does not settle a
new parser acceptance policy.

| Input | Relations | Point to preserve |
| --- | --- | --- |
| `C/C=C/C=C/C` | O, O | The central single-bond marker participates in both configurations. |
| `C/C=C/C=C\C` | O, S | The two configurations differ while sharing a marker. |
| `C/C=C/C=CC` | O, — | The shared marker completes the first site but leaves the second unspecified. |
| `CC=C/C=CC` | —, — | A shared marker alone completes neither site. |
| `C/C=CC=C/C` | —, — | Separated outer markers do not establish a relation across the intervening unmarked single bond. |
| `F/C=C/CC=CC` | O, — | A saturated spacer separates the two sites; only the first is specified. |
| `C/C=C/C=C/C=C/C` | O, O, O | Marking both outer sites through the backbone also specifies the middle site. |
| `C/C=C/C=CC=C/C` | O, —, — | Removing one internal marker also loses the third site's configuration. |
| `C/C=C/C=CC(/[H])=C/C` | O, —, O | A marker to an explicit H completes the third site without completing the middle site. |
| `C/C=C/C=CC(\[H])=C/C` | O, —, S | Reversing that H marker changes only the third site. |

The triene examples are a writer constraint, not a proposal to unfold hydrogens
implicitly. With only the backbone reference bonds available, independently
emitting the two outer configurations can invent a middle configuration. The
explicit-H examples demonstrate another encoding when that atom is present.
Boundary conversion must preserve the implicit/explicit hydrogen distinction;
it cannot introduce an explicit H to obtain that encoding. No coordinate
generation is involved.
OpenSMILES's discussion of shared directional bonds likewise identifies the
need for coordinated marker assignment. The proposed per-site frames do not
remove that output problem.

| Input | Current interpretation | Point to preserve |
| --- | --- | --- |
| `C/C=C/1CO1` | S | Direction supplied at ring opening. |
| `C/C=C1CO\1` | S | Equivalent direction supplied at ring closure. |
| `C/C=C/1CO\1` | S | Both ends consistently specify the same ring bond. |
| `C/C=C/1CO/1` | Parse: MismatchedRingBondDirections | Equal glyphs at opposite endpoints conflict. |

These examples support one configuration authority for complete ordinary
cis/trans input. Under the settled equivalence rules, partial spellings that
establish no configuration are normalized away after checking consistency.
None of the successful partial examples produces Either in current raise. Do not use Either
as a generic container for incomplete slash notation. Cumulene axes and other
extended stereo remain outside this local bond-frame example set.

The probe is reproducible with:

```sh
cargo run --offline --manifest-path scratch/stereo-valence-scan/Cargo.toml \
  --target-dir /Users/dr/.cargo-target --bin directions \
  < scratch/stereo-valence-scan/directions.smi \
  > scratch/stereo-valence-scan/directions.jsonl
```

#### Settled output expressibility contract

Assign direction markers jointly across each connected group of double-bond
sites coupled through eligible marker bonds. The assignment must express every
definite frame with the correct relation and leave sites without assertions
unspecified. It may select different existing substituent bonds from the stored
references, with the corresponding configuration transport, but must preserve
the implicit/explicit hydrogen distinction.

If no compatible assignment exists in the supported output notation, report a
representability failure. Failure of one greedy assignment does not establish
that no encoding exists. An explicit Either assertion also needs supported
notation: plain SMILES omission of markers supplies no assertion and cannot
silently substitute for Either.

Point 3's preservation and failure policy is settled. Marker assignment and
establishing its completeness are implementation work. First identify the groups
of sites coupled by shared eligible marker bonds, including unasserted sites
whose accidental specification must be prevented, then solve each group jointly.
This search may warrant a reusable graph primitive. Inspect existing graph views
and traversal operations before adding one; its required relation is marker-sharing
connectivity, not necessarily the chemistry model's broader notion of conjugation.
No new primitive or public API is prescribed by this implementation consideration.

#### Downstream consequences

| Consumer or producer | Consequence of explicit frames |
| --- | --- |
| SMILES reading | Interpret directional markers in their endpoint viewpoints and publish their configuration/frame where determined, preserving partial information separately as required. |
| CTfile reading | Interpret the supported codes and supplied geometry into the same frame meaning; geometry absence or degeneracy does not authorize completion. |
| Graph-IR raise | Complete identified ligand frames can map to stereo-bond entries. The smaller reference frame expresses a relative assertion for molecular interpretation; site-only unknown assertions remain distinguishable. |
| SMILES boundary conversion and writing | Find directional marks expressing the supplied frames collectively. A single bond shared by two sites must receive one compatible direction. Preserve unspecified sites as well as definite ones. |
| CTfile boundary conversion and writing | Determine whether the chosen format and supplied coordinates can express the stored configuration. A complete frame preserves the information but does not make every destination capable of writing it. |

The SMILES writer's choices are therefore coordinated across adjacent sites, even though TableIR
stores configurations locally. It must encode the given configurations without introducing a
configuration at another site. This is notation assignment, not chemical resolution. OpenSMILES
§§3.8.3 and 3.8.8 describe directional and partially specified double-bond stereo in the
[external specification](https://opensmiles.org/opensmiles.html).

Once frames carry the interpreted configuration, retained directions, codes, and coordinates cannot
remain competing authorities that consumers silently select between. The format-reading boundary
must establish their relationship to the frame. Adding a frame and merely skipping all old
information during raise would leave contradictions and partial information unresolved.

The smaller actual-atom frame for definite configurations and site-only Either assertions are settled above.
Partial-source equivalence and output preservation are also settled above.
Store stereo-bond records alongside stereo_atoms in the applicable TableIR molecule carriers.
The conversion API direction below governs boundary construction; implementation follows the plan.

### Coordinates and cross-format conversion

For ordinary MOL double-bond stereo, coordinates carry the definite cis/trans configuration.
V2000's double-bond stereo field has code 0 (derive cis/trans from atom coordinates) and code 3
(either). It has no separate definite cis and definite trans codes. The V3000 bond CFG field is
also not an E/Z field: its values describe none/up/either/down. See the bond-block tables in the
[CTfile specification](../materials/formats/mol/ctfile.pdf), printed pages 12 and 83.
[RDKit's MOL stereo documentation](https://rdkit.org/docs/RDKit_Book.html#from-mol) confirms
coordinate-based double-bond interpretation for both 2D and 3D and the V3000 CFG=2 unknown marker.
This concerns ordinary molecular stereo, not an encoding in custom properties or query constraints.

Consequently, without usable coordinates ordinary MOL cannot serialize the distinction between
the two definite configurations. All-zero coordinates likewise provide no orientation. Atom parity,
the chiral flag, and stereo-care flags do not substitute for a definite double-bond descriptor.
An explicit either mark can survive the absence of geometry; a definite frame cannot be replaced
by that mark without losing information.

The smaller TableIR reference frame makes this asymmetry explicit:

- SMILES can supply a definite reference frame without coordinates.
- MOL without usable coordinates can supply no definite configuration through its ordinary
  double-bond stereo fields; an explicit either assertion remains distinguishable from no assertion.
- MOL with usable coordinates can be interpreted into a definite reference frame. Once extracted,
  the frame can support SMILES output without retaining the coordinates.
- A definite TableIR frame without usable coordinates cannot be faithfully converted back to
  ordinary MOL. If Ctfile promises a serializable format value, its construction must report this
  representability failure. Keeping hidden configuration data in the wrapper would not make its
  serialized output faithful.

TableIR therefore carries the combined expressive capabilities of the supported formats, not only
their intersection. Adding duplicate stereo fields cannot make MOL express a geometry-dependent
configuration without geometry. Reading into TableIR and then removing coordinates also differs
from discarding coordinates before stereo interpretation: only the former can retain an already
extracted definite frame. Coordinate removal itself is an explicit transformation.

Coordinates are supplied data. The proposed Ctfile boundary either has them or does not. Missing
coordinates must remain missing; coordinate generation, including any future distance-geometry
operation, is an explicit caller step outside projection, conversion, and output. Degenerate or
contradictory supplied coordinates likewise do not authorize implicit generation or repair.

The design must cover the following paths before choosing the TableIR stereo representation.
These are required conversion cases, not claims of implemented CTfile output support:

| Path | Information available and required decision |
| --- | --- |
| Smiles → TableIR → Ctfile, no coordinates supplied | A definite ordinary double-bond configuration is not expressible in MOL without usable geometry; report the conversion failure and retain the TableIR frame. |
| Smiles → TableIR → Ctfile, coordinates explicitly supplied | Establish correspondence between coordinates and atoms and the compatibility required by the chosen stereo encoding. Preserve the supplied geometry. |
| Ctfile with coordinates → TableIR → Smiles | Interpret supported source stereo using its declared conventions and geometry where required, then carry explicit configurations into SMILES traversal. Ordinary SMILES does not carry the coordinates. |
| Ctfile without coordinates → TableIR → Smiles | Interpret only what the remaining notation actually determines. Preserve expressible unspecified state; do not infer a definite configuration from missing geometry or discard an unrepresentable assertion. |
| Ctfile → TableIR → Ctfile, with or without coordinates | Preserve coordinate presence and supported stereo meaning. Distinguish retaining source notation from reconstructing it after changes to atom order or configuration. |

Coordinate presence alone is insufficient: the CTfile variant, coordinate dimensionality, and
specific stereo encoding determine what can be interpreted and emitted. The supported cases need
an explicit representability inventory. Storing two encodings neither proves nor disproves that
inventory, and missing coordinates are not by themselves an argument against shared storage.

The proposed Ctfile wrapper and its ownership contract are not implemented by the current parser;
doc 153 T2 tracks the missing MOL/SDF boundaries. TableIR currently has optional positions. Their
role in a future Ctfile boundary must be settled without adding coordinates to graph-IR Molecule
or expecting projection from that molecule to recover them. A path through graph IR alone cannot
promise preservation of source coordinates; a boundary-level roundtrip must retain them separately.

## Smiles and ReactionSmiles as output boundaries

Keep their private payloads. Convey returns Smiles or ReactionSmiles, using
TableIR as its internal conversion intermediate. Checked from_table_ir constructors
on those boundary types accept a table and IO configuration and establish renderability.
This supersedes the earlier rejection of public TableIR-to-boundary conversion. Existing
as_table_ir and into_table_ir accessors remain unchanged.

The boundary guarantee is settled: successful convey or from_table_ir construction guarantees
rendering under the same SmilesIoConfig, following the DSL/defaults precedent. Construction establishes
representability, including feasibility of coordinated stereo-marker assignment;
rendering under that configuration cannot subsequently report an unsupported
representation. Writing to an external sink may still fail with an I/O error.
Rendering under another configuration may fail if that configuration cannot
express the value. Neither operation performs implicit chemical transformations.

Rendering recomputes traversal and marker assignment for now. The guarantee does
not itself prescribe cached data, a stored configuration, or another public type.
It applies to exported values and
does not silently strengthen the existing parser's acceptance contract. Exact
export and rendering names and error types remain to be specified.

ReactionSmiles composes the molecular writer over the three ordered sections.
Rendering a parsed boundary value writes Atom.class, preserving input labels,
repeated classes, and agent annotations even where interpretation as a graph-IR
Reaction is unavailable.

The label/index relationship is settled by the existing parser: Atom.class
stores the label, and Reaction.atom_mapping indexes those labels as class to
reactant and product atom-index lists. Builder::on_atom populates both together.
Agent atoms retain their classes but do not contribute to atom_mapping.
Interpret consumes that index to construct graph-IR correspondence and rejects
repeated classes that prevent an unambiguous correspondence; parsing retains them.

Graph-IR convey assigns label i + 1 to both atoms of matched correspondence pair i, using
the existing matched-pair order (sorted by left ID). It stores those labels in Atom.class
and derives atom_mapping from them; unmatched atoms receive no correspondence label. The export producer
establishes consistency just as parsing does. There is no independent authority
choice between these fields and no need to introduce general validation machinery
for independently assembled TableIR on this path. Point 6's mapping policy is settled.

Graph-IR reaction export materializes its sides and uses their correspondence to assign map labels.
It must not infer a new atom map. A Reaction may fail to materialize a ReactionSpan; that failure
precedes molecular projection. Numeric source map labels and agents are not recoverable from the
graph-IR Reaction alone. Export can produce an empty agent section and deterministic fresh map
labels, but cannot claim a boundary-exact roundtrip of discarded input metadata or arbitrary
reaction constraints.

Existing Python from_smiles and from_reaction_smiles ingestion paths must migrate with the Rust
resolver-argument change. Additional Python output or direct boundary-type APIs require a separate
public-surface decision; the plan does not invent Python Resolver or boundary wrappers.

## Roundtrip guarantees and evidence

There are two laws, with different domains and equality relations:

1. Boundary parse/render follows doc 155: preserve the ordered boundary meaning, accept the result
   under the same IO configuration, and make syntax normalization idempotent. Source spans and
   arbitrary original spelling are not reconstructed. Raw TableIR equality is not this law.
2. For supported, successfully ingested foreign input, exporting the ingested result preserves
   the input's molecular meaning under the agreed equivalences, allowing entity renumbering and
   stereo-frame transport. At the resolver layer this is project(resolve(input)) ≡ input.
   Boundary metadata absent from graph IR is outside this law. A general inverse for arbitrary
   constructed GraphIR is not required, nor is runtime re-resolution to verify projection.

A parser and writer can agree on the same mistake. Verification therefore needs independent known
configurations and semantic comparisons, including changed traversals at ring-closing stereo sites,
virtual ligands, and coupled directional-bond choices. Doc 224 supplies concrete motivating failures.
These are evidence requirements, not new test edits or an implementation sequence. Existing parser
and resolver benchmarks are starting points; projection and formatting need their own measurements
when algorithm work begins.

Property coverage must enforce the construction boundary guarantee for both Smiles
and ReactionSmiles: every successful convey or from_table_ir under an IO configuration renders
successfully under that same configuration. Include coupled stereo-marker cases.

### Independent parser-validation follow-up from doc 224

The user reports completing the at-least-60-minute SMILES parser fuzzing campaign. The run's
terminal log, execution count, and loaded-corpus/configuration details have not been independently
inspected here. Do not invalidate that run because of the inner catch_unwind: the 2026-09-10
follow-up below verified that libfuzzer-sys aborts in its panic hook before unwinding. A run under
that normal runtime still exercises crash detection, though not chemical/stereo equivalence.
The target and seed follow-up is independent of the roundtrip implementation stages.

Local formatting references include RDKit's SmilesWrite::FragmentSmilesConstruct in
[SmilesWrite.cpp](../materials/codes/rdkit/Code/GraphMol/SmilesParse/SmilesWrite.cpp), which consumes a
traversal stack and emits atom/bond/ring tokens, and Indigo's traversal-relative tetrahedral mapping
in [smiles_saver.cpp](../materials/codes/Indigo/core/indigo-core/molecule/src/smiles_saver.cpp).
These establish useful writer structures to study. They do not establish umol projection semantics
or require their chemistry models or canonicalization policies.

## Reference implementation findings

Reviewed local `materials/codes` sources on 2026-09-10 for the remaining writer
algorithms and API vocabulary. This is source inspection, not an implementation
plan, benchmark, or claim of equivalence to umol. No external writer was executed.
The inspected checkouts had no tracked modifications:

| Source | Commit |
| --- | --- |
| RDKit | `dfbce43dd0c2a78a36c48be25903b78b8aea166a` |
| Indigo | `8642d73cca2b827a5351d8a8eb21816492381497` |
| CDK | `09edb823cb2b3ea15c929634a521ba2be46ba575` |
| Open Babel | `1afe81019b781c857eecb86251e106e2f31da928` |

### Traversal, output order, and atom-frame transport

RDKit's SmilesWrite::FragmentSmilesConstruct in
[SmilesWrite.cpp](../materials/codes/rdkit/Code/GraphMol/SmilesParse/SmilesWrite.cpp)
(line 366) first obtains a MolStack from Canon::canonicalizeFragment, then emits
atom, bond, ring, and branch entries. It records atom and bond output order.
MolToSmiles supplies canonical ranks or sequential atom-index ranks (line 632),
so the traversal machinery is not intrinsically tied to chemical canonicalization.
Ring labels are allocated from currently available numbers and released after
closure; the delayed release prevents reuse at the same atom. The numeric labels
are output bookkeeping, separate from identifying closure edges.

Canon::canonicalizeFragment in
[Canon.cpp](../materials/codes/rdkit/Code/GraphMol/Canon.cpp) (line 1406) collects
the actual traversal bond order and computes tetrahedral permutation parity
against it, with root/ring treatment, before text emission. Its hydrogen-suppression
adjustments are not transferable to umol's explicit/implicit-H contract.

Indigo's SmilesSaver in
[smiles_saver.cpp](../materials/codes/Indigo/core/indigo-core/molecule/src/smiles_saver.cpp)
uses DfsWalk with optional vertex ranks (line 187). It then builds per-atom
written-neighbor order (line 254), reserving ring-opening slots and filling them
when the closing edge is encountered. The tetrahedral pass (line 321) orders the
incoming ligand, virtual participant, and remaining neighbors before calculating
pyramid permutation parity. This is directly relevant to doc 224: output stereo
must use planned encounter order, not merely the bond table's order. The code
also records written atoms/bonds and manages reusable ring numbers separately.

Open Babel's OBMol2Cansmi::BuildCanonTree in
[smilesformat.cpp](../materials/codes/openbabel/src/formats/smilesformat.cpp)
(line 3130) similarly separates a traversal tree from output. It orders neighbors
by supplied ranks, with a multiple-bond preference to avoid putting those bonds
on closures; explicit user ordering disables that preference. Such preferences
are algorithm choices, not necessary molecular semantics.

The common useful structure is a traversal description sufficient to recover
parent edges, closure encounters, branch order, and output atom/bond order before
formatting stereo. This does not imply a new public plan type or a canonical
SMILES promise.

### Shared marker assignment and preservation limits

RDKit's Canon::getReferenceDirection (Canon.cpp, line 114) transports a stored
two-reference cis/trans relation to chosen controlling bonds, accounting for
reference substitutions and endpoint direction. Canon::canonicalizeDoubleBonds
(line 720) builds adjacency between definite stereo bonds through eligible
neighbor bonds. It prioritizes sites by neighboring stereo-site count and output
visit order, then processes connected sites with a queue. This is concrete
precedent for marker-sharing connectivity without a general aromatic/conjugation
perception operation.

The subsequent Canon::removeUnwantedBondDirSpecs (line 1223) checks double bonds
without definite stereo for accidental specification. It attempts to remove a
marker where an adjacent definite site has a redundant marker. Its source comment
explicitly describes the repair as incomplete, and a branch leaves the unwanted
assignment when removal is unavailable. Canon::canonicalizeDoubleBond also has
conflict paths that log warnings rather than fail (for example line 556).
Canon::removeRedundantBondDirSpecs is a further cleanup pass, not a completeness
proof. These routines offer useful ordering and direction-transport mechanics,
but do not establish umol's required preservation/failure guarantee.

Indigo's SmilesSaver::_banSlashes (line 985) explicitly examines unasserted sites
that would receive markers on both ends and bans implicated bonds. _markCisTrans
(line 1130 area) classifies definite sites with no allowed marker at an endpoint
as complicated. _updateSideBonds (line 1231) propagates parity through assigned
side bonds and throws on incompatible assignments; _calcBondDirection (line 1330)
repeatedly propagates assignments and seeds a previously unassigned direction.
The bans are conservative: inspecting this procedure does not establish that it
finds every possible marker selection.

Complicated cases are routed through writeRingCisTrans (line 2152), which writes
CXSMILES cis/trans annotations when extension output is enabled. In the ordinary
molecular saver this call is gated by write_extra_info and chemaxon (line 702).
This is not proof of faithful plain-SMILES output or a fallback umol can adopt
without the corresponding supported notation and boundary checks.

Open Babel's OBMol2Cansmi::GetCisTransBondSymbol (line 2566) assigns directions
lazily, remembers already assigned bond directions, and adjusts a site's signs
to shared assignments. Its CreateCisTrans filters out small-ring sites, and
GetSmilesElement can suppress explicit H atoms by default. Those policies do not
meet umol's preservation contract. This review found useful local direction
bookkeeping, not evidence of a complete solver preserving every unasserted site.

For umol, grouping, selection of which bonds carry markers, and orientation
assignment are distinct algorithmic concerns. Connected-component search only
solves the first. Once marker choices are fixed, the same/opposite relations
supply parity constraints; preserving unasserted sites also constrains marker
presence. That is an inference from the inspected mechanisms and our settled
contract, not a claim that any reviewed writer supplies the complete algorithm.

### Boundary conversion, reaction composition, and naming

CDK separates molecular-model conversion from its Beam writer. In
[CDKToBeam.java](../materials/codes/cdk/storage/smiles/src/main/java/org/openscience/cdk/smiles/CDKToBeam.java),
addGeometricConfiguration (line 324) passes the bond endpoints and one reference
at each end to GraphBuilder.geometric(...).together/opposite. Tetrahedral
configuration is passed as an ordered ligand frame. This closely matches the
boundary information we have selected, without establishing umol's construction
guarantee. Beam's implementation is an external dependency and was not present
among the inspected sources, so its solver behavior is not inferred here.

CDK's SmilesGenerator.create in
[SmilesGenerator.java](../materials/codes/cdk/storage/smiles/src/main/java/org/openscience/cdk/smiles/SmilesGenerator.java)
(line 442) converts to Beam, optionally permutes canonical labels and normalizes
directional labels, then calls toSmiles with an output-order array. Its optional
resonance and atom-map-renumbering operations must not be imported as implicit
umol output steps. The reaction overload (line 576) joins the three molecular
sections and adjusts their output-index offsets.

RDKit's
[ReactionWriter.cpp](../materials/codes/rdkit/Code/GraphMol/ChemReactions/ReactionWriter.cpp)
composes molecular writers with dots and reaction separators, optionally sorting
component strings for canonical output. Atom map numbers are written by the
molecular atom formatter (SmilesWrite.cpp, line 235). Indigo's
[rsmiles_saver.cpp](../materials/codes/Indigo/core/indigo-core/reaction/src/rsmiles_saver.cpp)
likewise delegates each molecule to SmilesSaver and collects written atom/bond
orders for reaction-wide extensions. These support the settled composition and
label-ownership design; they do not choose umol's fresh correspondence-label order.

The names MolToSmiles, SmilesGenerator.create, and SmilesSaver.saveMolecule reflect
different library boundaries. None is a reason to copy a public arbitrary-TableIR
conversion or a new intermediate wrapper into umol. Public export/render names
must follow the existing umol operation vocabulary. The sources do not settle
whether an export result retains a traversal/assignment or recomputes it later;
the inspected APIs predominantly construct text within the same operation.

### Existing umol primitives and discussion consequences

[traversal.rs](../umol-graph-core/src/algorithms/traversal.rs) currently exposes
Graph::neighborhood, a distance-limited BFS, not a general DFS event/tree API.
[components.rs](../umol-graph-core/src/algorithms/connectivity/components.rs)
provides Graph::enumerate_connected_components with an explicit algorithm selector.
If the IO operation constructs the marker-sharing graph, its connected components
can use that existing operation. A new chemistry-wide conjugation primitive is
not established as necessary by this review. Whether to materialize that graph
or traverse its incidence directly remains an algorithm/storage decision.

The review informs the settled designs recorded below:

- shared graph-core DFS/BFS traversal with the settled callback contract and API names;
- apply the marker component, selection, and orientation design recorded below;
- rendering recomputes traversal and marker assignment; do not retain them in the boundary for now;
- correspondence-pair indices plus one for fresh reaction labels, and the convey/export/render APIs.

The settled semantic rules remain unchanged. The foreign writers' recovery,
hydrogen suppression, canonicalization, and extension fallbacks are not adopted.
No implementation sequence or new public symbols are approved by these findings.

## Conversion API direction

Use Convey and convey as the counterparts to Interpret and interpret. Both traits belong to
umol-graph, where the chemistry-dependent conversion is implemented for the umol-io boundary
types. Convey is implemented on the output boundary type and has an associated input type;
Interpret remains on the input boundary type with its associated output type. Do not introduce
a generic Export target on Molecule or Reaction.

The intended call shapes are:

```text
smiles.interpret(&resolver)                    -> Molecule
Smiles::convey(&molecule, &resolver, io_config) -> Smiles

reaction_smiles.interpret(&resolver)                       -> Reaction
ReactionSmiles::convey(&reaction, &resolver, io_config)    -> ReactionSmiles
```

These sketches omit Result and error types. Interpret is to receive a supplied Resolver rather
than separate chemistry-model and resolve-config arguments. Ingest and export likewise use a
supplied Resolver. Convey invokes graph-IR projection, converts the resulting graph IR to TableIR,
and constructs the checked boundary. Resolver::project operates only on graph IR, as does resolve;
neither resolver operation knows about TableIR. Project follows resolve's mutation semantics:
work on an intermediate result and publish to the supplied molecule only on determined success.
Convey clones the caller's graph IR before projection. Parse and render
remain boundary operations in umol-io. This separates the conversion name convey from the
convenience-function name export_smiles.

The agreed trait shape is:

```rust
pub trait Convey: Sized {
    type Input;
    type Config;
    type Error;

    fn convey(
        input: &Self::Input,
        resolver: &Resolver<'_>,
        config: &Self::Config,
    ) -> Result<Self, Self::Error>;
}
```

Smiles uses Molecule as Input; ReactionSmiles uses Reaction. Both use SmilesIoConfig as Config.
Interpret retains associated Output and Error and takes `&self` and `&Resolver`.
Convey uses checked from_table_ir constructors owned by umol-io; the payload fields stay private.
These constructors establish representability under the supplied IO configuration, including
marker-assignment feasibility. They provide the construction counterpart to as_table_ir and
into_table_ir without exposing unchecked fields or direct TableIR formatting.

export_smiles and export_reaction_smiles return text, composing convey and render, opposite to
ingestion from text. Ingest/export take a supplied Resolver; their default forms use the default
OpenSMILES IO configuration, and `_with` forms additionally accept an explicit IO configuration.
Parse/render and from_table_ir require no resolver.

Error responsibilities follow the operations: ProjectError covers operational projection failure,
with the result structure following resolve; SmilesConstructionError covers boundary
unrepresentability; ConveyError distinguishes projection from construction failure;
SmilesRenderError covers rendering under the requested configuration; SmilesOutputError composes
convey and render failures. Reaction errors add side context where appropriate. Construction and
rendering share representability diagnostics. Successful construction precludes a later
representability failure under the same configuration. Exact variants and reaction error names
remain implementation-design details. No production implementation is authorized yet.

## Resolver projection: roundtrip target and initial scope

The required law, under the same resolver, is:

```text
project(resolve(input)) ≡ input
```

Its domain is supported, successfully ingested SMILES/MOL input. At the GraphIR layer, input means
the raised representation; equivalence allows the already agreed defaults, redundant notation,
ordering, and stereo-frame transport. It does not require raw GraphIR or source-text equality.
Convey/render/parse tests establish the corresponding boundary law. Implicit H and explicit H
atoms remain distinct. CTfile inputs inform projection coverage; a CTfile writer remains separate.

The opposite law, resolve(project(molecule)) ≡ molecule for arbitrary constructed GraphIR, is
not a required deliverable. Experiments in that direction on ingestion-produced values may help
find defects, but do not expand the contract. No provenance tracking or runtime test of ingestion
origin is needed.

Projection translates the existing state directly. It does not re-admit valence candidates,
apply tie-breaks, rediscover rings, select aromatic assignments, or resolve a copy to verify its
answer. Uncertainty about reconstruction, a hypothetical different partition, or a solver search
limit is not an export failure category.

The lowering rule is simple: lower a projected value if TableIR can carry it; otherwise return
an error. The SMILES boundary rejects TableIR features that its supported notation/configuration
cannot encode. These checks belong to the first operation requiring the corresponding property.
Resolver projection remains GraphIR-only; Convey owns GraphIR-to-TableIR conversion. A failure
identifies a concrete unsupported property or value, not a failed reconstruction proof. Start
with a useful writable subset and expand it through concrete cases; minimum encodings and general
valence inversion are not prerequisites.

Use property tests of the required input-domain law and boundary roundtrips, backed by independent
expected outputs and required-success cases. Runtime projection does not perform these tests.

### Ordinary atom valence

Preserve element, isotope subject to its settled default policy, charge, H count, and bond orders.
Bracket H counts remain implicit hydrogen participants, not explicit hydrogen atoms. Apply the
TableIR-capacity rule to atom electron fields and the supported-notation rule at the SMILES
boundary. Do not clear lone pairs, unpaired electrons, or multiplicity merely to test whether a
valence resolver would derive them again. Do not infer unrestricted SMILES support from a numeric
storage field: radical/spin encodings have a limited set of supported states, with extensions
depending on configuration. A concrete unencodable value fails at its owning conversion.

S4a removed the implemented candidate reconstruction and comparison machinery, including clearing
done solely for that check. With atom fields retained, there is no ordinary valence projection
operation: ValenceResolver::project and ValenceProjectError are removed rather than replaced by
an empty method. Forward valence resolution, registry invariants, and exact post-meet duplicate
elimination retain their own contracts. Neither a registry pruning API nor a reverse counts
algorithm is needed for this work.

### Isotope resolution and registry raising (settled during S4a)

The isotope rejection exposed during S4a is a resolution inconsistency to repair, not an
unsupported projection category. Before S4a0a, registry raise inherited Natural from concrete atom
defaults, while counts admission independently filled an undetermined isotope with Natural.
Consequently an explicit isotope could conflict with an otherwise applicable default registry row.

IsotopeResolver is a separate phase, after assertion placement and before valence admission.
Its configuration is independent of valence's Strict/MostSaturated policy. The isotope policies
are Strict, which leaves Undetermined unchanged, and Natural, which completes Undetermined with
Natural. Both preserve supplied Natural and explicit isotope values; neither selects a mass
from an isotope set or overrides a variable. The isotope policy therefore takes effect before
custom registry patterns are considered. Format-owned interpretation already present in GraphIR
is preserved. Strict is the general default isotope policy. The SMILES defaults select Natural,
independently of their MostSaturated valence policy. Explicitly supplied resolver configuration
remains authoritative.

IsotopeResolver::project uses the same policy as resolve, analogous to omitting defaults during
lowering. Under Natural, project replaces Natural with Undetermined because resolution restores
Natural. Under Strict, project retains Natural because resolution would not restore an omitted
value. Both policies retain explicit isotope masses. No separate projection policy or provenance
distinguishing supplied Natural from defaulted Natural is needed: the agreed default equivalence
does not preserve that distinction. This operation stays within GraphIR; boundary
conversion and rendering own the output notation.

Valence resolution must stop supplying isotope defaults. It may complete its own fields while
isotope information remains undetermined; the composite resolver retains its final completeness
check and atomic publication semantics. Aromaticity remains directly after valence.

Registry raise uses concrete defaults for inherent fields other than isotope, together with its
existing valence-related constraint defaults. Override isotope with the existing Required policy,
which leaves an omitted isotope Undetermined. Every stored AtomTypeRegistry entry must have
isotope Undetermined. Enforce this invariant at every construction and insertion boundary,
including programmatic AtomForms, TOML loading, the registry macro, and Python construction.
Reject supplied Natural, masses, sets, variables, and every other non-Undetermined isotope form;
do not strip information to make an entry admissible. Validate before insertion so rejection
leaves the registry's contents and hash unchanged. Private storage and read-only accessors preserve
the invariant after construction; lookup and resolution need no repeated checks.

Programmatic construction and insertion have checked try_from_atoms and try_add methods returning
ConfigError, with from_atoms and add as asserted counterparts using the same entry check. TOML
loading and Python construction use the checked gate; Python reports ValueError. The gate also
owns literal element/charge requirements and the i8 charge-key range check already used by Python.

The broader RequireUndetermined DSL policy and fallible raise/lower migration are dropped.
IntoIr, FromIr, and general molecule/reaction defaults keep their existing contracts. The
restriction belongs to the registry that requires it. ValenceTable stores numeric states rather
than AtomForms, so it needs no corresponding isotope validation after counts admission's
programmatic isotope completion is removed.

Exact post-meet duplicate elimination is also settled: AtomTypingValence admission retains one
copy of each exactly equal AtomForm, including constraints, in existing first-occurrence order.
This is independent of isotope policy and valence tie-breaking. Distinct completed states remain
distinct. Custom overlapping registry patterns supplied the motivating reproduction. Its isotope
variants become inadmissible under the registry invariant; regressions must exercise admissible
non-isotope overlaps instead, preserving the exact-equality law.

S4a0a–S4a0c remain complete. Their registry and isotope corrections are independent of the
S4a projection correction. Natural and explicit masses remain supported under both
candidate sources without a valence-specific isotope policy.

### Direct aromatic projection

Projection fails for nonzero charge or `#u > 0` on any bond or aromatic system; it does not
discard or implicitly localize them. This projection restriction complements the general export
representability rule above.

Read each existing system's members and their electron contributions. Write concrete atom
`#a<n>` using each member's contribution n, mark that system's bonds `#a`, and remove the system
entity. Preserve localized connections between systems. This phase preserves ordinary atom
lone-pair, unpaired-electron, and multiplicity fields. No valence source, tie-break, ring discovery,
or joint selection is involved; a molecule without aromatic systems needs no aromatic edit.

Convey marks participating TableIR atoms aromatic and encodes the systems' bonds as aromatic
bonds. SMILES/MOL cannot explicitly encode contribution n; their inputs did not supply that
quantity directly. This does not justify a runtime reconstruction check or a blanket export
failure for aromatic molecules.

The information accounting is:

| Source information | Projection or boundary treatment |
| --- | --- |
| System members and member electron contributions | Concrete atom `#a<n>`; aromatic marking at the boundary |
| Bonds belonging to the system | Bond #a; aromatic bond order at the boundary |
| Localized links and atom-localized charge | Preserve |
| Nonzero system or bond charge, #u > 0, or non-singlet system spin | Error; no implicit localization or redistribution |
| Other projected values and bond kinds | Lower if TableIR supports them; otherwise error. The SMILES boundary separately enforces its supported notation. |

Do not retain blanket SystemConstraints/MoleculeConstraints rejection merely to enforce an
arbitrary-GraphIR inverse. Forward resolution discharges the system ElectronCount assertion when
its derived value is ground. General preservation of independently constructed assertions about
deleted system IDs is not an additional deliverable for the supported ingestion domain.

### Stereo projection

For ordinary tetrahedral and cis/trans stereo, reuse the assertion frames constructed by
StereoPerception::derive_stereo_atom and derive_stereo_bond in
[stereo.rs](../umol-graph/src/ops/stereo.rs):

- The `#T` frame has actual neighbors in GraphIR neighbor order, followed by an implicit-H or
  lone-pair participant when needed to complete four ligands.
- The `#C` frame follows the bond's ordered endpoints, each with its actual substituent neighbors
  excluding the opposite endpoint and a virtual participant when needed to complete that side.

Obtain the frame the same stereo model derives, use coset_for to transport the stored
entity configuration into it, write the corresponding `#T` or `#C` assertion, and remove the
entity. Failure to construct a frame, incompatible ligand identities, or unrepresentable
configuration causes projection failure.
Do not relabel lone pairs as hydrogen, fold or expand hydrogen atoms, or silently drop stereo
under permissive failure policies. Use configured perception where needed to derive the frame;
do not run forward stereo resolution as an export check. TableIR frame construction remains
entirely with Convey.

## Marker-assignment design

### Candidate components and eligibility

Construct marker-assignment components before choosing which bonds receive markers. Candidate
bonds are vertices of an auxiliary graph; each relevant double-bond site groups its candidate
side bonds, and shared candidate bonds merge those groups. A spanning star suffices to connect
each site's group. Preserve the site's candidates separately by endpoint alongside component
membership; a flat component list alone loses the relationships needed for assignment.

Eligibility is structural, without element lists, CIP ranking, chemical equivalence filtering,
or ring-size restrictions:

- A site is a localized ordinary double bond in the settled two-endpoint stereo scope. Each
  endpoint has one or two actual substituent neighbors excluding the opposite endpoint.
  Extended cumulene and axial stereo remain separate work.
- A candidate is an actual single bond from a site endpoint to a substituent that can carry
  a SMILES direction marker. Branch and ring-closure bonds are eligible. Virtual participants
  provide no candidate bonds; aromatic, double, dative, and noncovalent bonds are not candidates.
- Definite, Either, and unasserted sites participate. A site with candidates at both endpoints
  connects all those candidates. Include unspecified sites because marker choices can couple
  through them even though they require no definite configuration.
- Retain components containing a definite assertion. A definite site without candidates at
  both endpoints fails export rather than disappearing from consideration. An unspecified
  site missing candidates at an endpoint cannot acquire a complete directional assertion and
  need not connect groups across that site.

Branches retain every eligible continuation; cycles are handled by ordinary component traversal.
No molecular path or ring opening is selected during this construction. Molecular cycles and
cycles in candidate connectivity are distinct. The auxiliary graph is marker-sharing connectivity,
not a general chemical conjugation classification. The current cis_trans_capable helper is not
the full eligibility predicate and must not define it by precedent.

### Marker selection

Use one Boolean per candidate bond for marked/unmarked. An endpoint is covered when at least one
of its candidates is marked. Definite sites require both endpoints covered. Either and unasserted
sites forbid simultaneous coverage of both endpoints. A shared candidate has one Boolean across
all its occurrences.

Use deterministic constraint search with propagation:

- A required endpoint with one remaining candidate forces that candidate marked.
- Coverage at one endpoint of an unspecified site forces all candidates at its other endpoint
  unmarked.
- A required endpoint with no marked or remaining candidates is a conflict.
- After propagation, branch on an undecided candidate, trying unmarked first.

For A–B=C–D=E–F=G–H with definite B=C and F=G but unasserted D=E, the outer assertions force
C–D and E–F marked, which would specify D=E. That candidate arrangement has no solution;
additional branches can provide alternative candidates. A failed greedy choice is not proof
of unrepresentability: exhaust the search before reporting no faithful assignment.

There is no minimum-marker objective. Accept the first complete solution, including consistent
orientation. A selection satisfying coverage can still fail orientation around a cycle; reject
that selection and continue searching. Initially check orientation after complete selections.
Earlier parity propagation is a possible optimization only if later benchmarks justify it.

### Direction propagation and deterministic ordering

For a fixed selection, transport the stereo frames into same/opposite parity constraints between
selected marker bonds, accounting for reference substitutions and bond-endpoint orientation.
Seed one direction per connected constraint group and propagate the rest. Revisiting an assigned
bond checks consistency; disagreement rejects the selection. Selection may split the candidate
component into several constraint groups, each needing its own seed. Flipping all directions in
one group preserves the represented stereo, so seed direction requires no search.

Use TableIR bond order throughout:

- Process candidate components by their lowest bond index.
- Branch on the lowest-index undecided candidate, trying unmarked first.
- Seed each constraint group at its lowest-index selected bond with `/` in the chosen output
  traversal orientation, then propagate directions.

No canonical ranking is introduced. Rendering recomputes traversal and marker assignments;
do not retain them in the boundary for now. The resulting markers must imply exactly the intended
definite sites and configurations, without introducing definite assertions at unspecified sites.
Explicit Either remains distinct from absence: unmarked bonds alone do not encode it, and its
output remains subject to the supported-notation boundary rules.

## Shared traversal design

The agreed direction is a graph-core traversal primitive with ControlFlow event callbacks.
Connectivity access is shared by DFS and BFS; their traversal state and event contracts remain
algorithm-specific. The signatures and names below are sketches, not implementation approval.

### Table ownership and temporary connectivity

TableIR remains atom and bond tables. It must not own a Graph or persist a second authoritative
connectivity representation. Operations may derive a temporary neighbor index over unchanged
tables. Graph-core owns the domain-independent traversal; IO owns the index over bond rows and
the interpretation of traversal events as SMILES structure.

Currently, raise allocates AtomNeighbors unconditionally before constructing graph IR. That
index covers every table bond and preserves table bond indices. Molecule::try_from_entries then
constructs Graph adjacency for localized bonds; dative and noncovalent table bonds become
relations instead. The final Graph therefore cannot generally substitute for the table index.
Explicit StereoAtom frames are translated directly; source wedge and double-bond interpretation
use neighbor lookup. The design direction is to construct temporary lookup only where needed,
and reassess those needs as explicit stereo-bond frames replace source interpretation in raise.
This ownership decision does not depend on benchmarks. A separate specialized IO traversal
would be a last resort requiring benchmark evidence against the shared implementation.

### Connectivity access

Use an iterator-producing callback rather than introducing a graph trait hierarchy or another
public storage type:

```rust
N: Fn(NodeId) -> I,
I: Iterator<Item = Neighbor>,
```

Neighbor is the existing graph-core node/edge pair. Graph supplies
`|node| graph.neighbors(node).iter().copied()`. A temporary TableIR index maps its atom/bond
indices to NodeId/EdgeId while iterating. Both are statically dispatched, borrowing iterators;
neither requires boxing, dynamic dispatch, per-call collection, or a generic associated type.
The callback may also supply filtered or differently ordered incidence, subject to the input
contract. The traversal does not impose chemical ranking or sort neighbors.

Node IDs index traversal state. DFS also needs an edge bound for edge visitation state; filtering
must not require renumbering retained edges. Each undirected edge has a stable identity, including
parallel edges. Input connectivity remains fixed during traversal. Incident relations do not
automatically become traversal edges: the caller chooses the connectivity being traversed.

### DFS events and ordering

```rust
pub enum DepthFirstEvent {
    Discover { node: NodeId, parent: Option<Neighbor> },
    NonTreeEdge { from: NodeId, to: NodeId, edge: EdgeId },
    Finish { node: NodeId },
    FinishTree { root: NodeId },
}
```

The visitor has the shape `FnMut(DepthFirstEvent) -> ControlFlow<B>`. A parent entry names the parent
node and connecting edge; None identifies a component root. Discover also reports tree edges,
so no separate tree-edge event is needed.

- Consider candidate roots in supplied order, skipping already reached nodes. Supplying all
  node IDs visits the whole graph; each first unreached candidate starts another component.
- Examine neighbors in iterator order. Retain a suspended neighbor iterator for each active
  DFS frame, using an explicit stack rather than recursion or collecting adjacency anew.
- On normal completion, emit Discover and Finish once per reached node with nested DFS
  lifetimes. Report each non-tree edge once, using edge identity to distinguish parallel edges
  and report a self-loop once even if represented by two incidences.
- Break returns immediately with the visitor's value. No synthetic Finish events follow it.

The writer collects the traversal structure needed for branches, ring labels, and stereo-frame
transport. Events need not emit text immediately: a ring edge discovered later can affect an
earlier atom's output. Component roots supply component boundaries without a separate component
enumeration pass.

### BFS events and existing consumers

```rust
pub enum BreadthFirstEvent {
    Discover { node: NodeId, parent: Option<Neighbor>, depth: usize },
    Finish { node: NodeId, depth: usize },
    FinishTree { root: NodeId },
}
```

The visitor has the shape `FnMut(BreadthFirstEvent) -> ControlFlow<B>`. BFS uses the same connectivity
callback as DFS, with a FIFO queue and node visitation state; it does not need DFS edge tracking.

- Traverse each candidate root's component before considering the next candidate, skipping
  nodes already reached. This is sequential component traversal, not multi-source BFS.
- Mark nodes when enqueued and emit Discover then. Parent records the first reaching edge;
  depth is shortest distance from that traversal root in the connectivity excluding nodes
  reached by earlier roots.
- Emit Finish after examining a node's neighbors. An optional maximum depth suppresses
  expansion at the limit, but nodes there still receive Discover and Finish.
- Preserve supplied neighbor order and stop immediately on Break.

Graph::neighborhood supplies one root and its depth limit, collecting `(node, depth)` from
Discover. Connected-component enumeration supplies all nodes without a depth limit, starts
a component at each parentless Discover, collects its nodes, and sorts each completed component.
Sorting is an output operation, not a change to BFS queue order. Component order follows the
first unreached candidate roots. With a depth limit, later roots can reach previously unvisited
parts of the same connected component; the component interpretation requires unrestricted depth.

### Replaced LIFO flood fill

Before S1c, despite ConnectedComponentsAlgorithm::Bfs and its documentation, the implementation in
[components.rs](../umol-graph-core/src/algorithms/connectivity/components.rs) used Vec::pop and
marked all unseen neighbors when pushing them. It was a LIFO reachability flood fill. It found
components correctly, and sorting their members hid visitation order in the returned result,
but it supplied neither BFS distances nor the proposed nested DFS event semantics.

For edges A–B, A–C, and B–C, that flood fill discovers B and C from A before exploring either.
Recursive-style DFS instead explores one and discovers the other through it. A stack alone does
not establish the DFS frame contract.

S1c replaces this component flood-fill loop with collection over the shared BFS implementation,
retaining sorted component results and the explicit Bfs selection. No separate
LIFO algorithm selector is needed for the initial design. DFS and BFS remain separate traversal
implementations using the common connectivity interface; a universal event engine is unnecessary.

### Callback contract

Retain the simple callback interface. Consistent undirected connectivity gives the specified
traversal; inconsistent input has no correctness guarantee, but the algorithm's own handling must
remain panic-free. There is no defensive integrity-validation pass, TraversalError, or Result
return. Panic-free handling of supplied indices does not establish or certify input integrity.
Callback execution remains the caller's responsibility: a callback that panics or hangs can cause
the operation to panic or hang. No isolation or recovery mechanism is required.

Graph construction establishes its connectivity integrity. A temporary index producer establishes
the index's internal integrity; subsequent mutation of the original tables carries no guarantee
that the index still corresponds to them. Do not introduce correspondence validation or wrappers
to preserve such a guarantee. The first operation requiring a property owns its check where that
property has not already been established; traversal does not defensively revalidate its producers.

### Traversal API and collectors

Use visit_depth_first and visit_breadth_first as graph-core free functions, with algorithm choice
explicit in the names. Their proposed signatures are:

```rust
pub fn visit_depth_first<R, N, I, V, B>(
    node_bound: usize,
    edge_bound: usize,
    roots: R,
    neighbors: N,
    visitor: V,
) -> ControlFlow<B>
where
    R: IntoIterator<Item = NodeId>,
    N: Fn(NodeId) -> I,
    I: Iterator<Item = Neighbor>,
    V: FnMut(DepthFirstEvent) -> ControlFlow<B>;

pub fn visit_breadth_first<R, N, I, V, B>(
    node_bound: usize,
    roots: R,
    max_depth: Option<usize>,
    neighbors: N,
    visitor: V,
) -> ControlFlow<B>
where
    R: IntoIterator<Item = NodeId>,
    N: Fn(NodeId) -> I,
    I: Iterator<Item = Neighbor>,
    V: FnMut(BreadthFirstEvent) -> ControlFlow<B>;
```

Provide inherent Graph::visit_depth_first and Graph::visit_breadth_first methods that supply
the bounds and CSR neighbor access to the free functions. Graph::enumerate_depth_first_events
and Graph::enumerate_breadth_first_events collect the corresponding event vectors through those
visitor methods. The `_events` suffix is settled and distinguishes these results from node lists.

Both event enums include FinishTree { root: NodeId }, emitted when the active DFS stack or BFS
queue for that root becomes empty, before considering another root. Finish remains per-node;
in BFS, finishing the root does not finish its tree. A depth-limited tree need not cover an entire
connected component, so the event is named FinishTree rather than FinishComponent. Break suppresses
all subsequent events, including FinishTree.

Add Graph::visit_connected_components with the following visitor contract:

```rust
pub fn visit_connected_components<B, V>(
    &self,
    algorithm: ConnectedComponentsAlgorithm,
    visitor: V,
) -> ControlFlow<B>
where
    V: FnMut(&[NodeId]) -> ControlFlow<B>;
```

Its Bfs implementation uses Graph::visit_breadth_first, buffering one component and sorting it
at FinishTree before visiting the borrowed slice. Break prevents traversal of subsequent
components. Graph::enumerate_connected_components collects owned vectors through that component
visitor, rather than maintaining another traversal implementation.

Keep Graph::neighborhood and ConnectedComponentsAlgorithm. Rename the neighborhood-specific
TraversalAlgorithm to NeighborhoodAlgorithm, retaining Bfs. Use usize consistently for traversal
depth, maximum-depth arguments, and returned distances. Update Graph::neighborhood and related
methods, callers, and bindings accordingly, including its current u32 depth argument and result.
This follows graph-core's size-argument convention; NodeId and EdgeId remain u32-backed.

### Semantic naming

The redesign includes revising existing and proposed names to match their actual semantics,
including algorithm selectors, operations, events, and connectivity terminology. Existing names
are not compatibility constraints. In particular, BFS must mean FIFO breadth-first traversal,
and DFS must mean the nested traversal contract rather than arbitrary LIFO flood fill. Reconcile
names with the final contracts before implementation. The API names settled above are part of
that redesign; no production implementation is authorized by this design record.

## Implementation handoff

The initial semantic scope, responsibility boundaries, traversal API, and marker-assignment
algorithm are settled above. Exact diagnostic variants and cohesive helper signatures must be
reconciled with the public contract before their implementing subitem; this is not permission to
introduce new wrappers, strategies, or chemistry transformations. Implement direct projection
and the supported boundary encodings with explicit failures for unrepresentable values.
The plan below sequences the work; S0–S3 are complete.

## Staged implementation plan

S0–S3, including allocation follow-ups S3d1–S3d7, are complete.
S4a0a–S4a0c and the S4a–S4b corrections are complete. S4c–S4d remain pending.
S5 and later stages are pending.
Every stage ends with a green tree; additive subitems remain green individually. Breaking
subitems include the consumer migration needed to restore green within that subitem or stage.
Each subitem is a reviewable unit, not an instruction to commit. No commits are authorized.

### Scope and module inventory

The consumer requirements determine these work groups; the stages order them foundation-first.

| Consumer/module | Required work |
| --- | --- |
| umol-graph ingest/export | Interpret with a supplied Resolver; Convey on boundary types; text-returning export conveniences and reaction composition |
| umol-io smiles boundaries | Checked from_table_ir construction, render/render_with, shared representability diagnostics, same-config guarantee |
| umol-io smiles writer | Deterministic traversal, atom/ring/branch spelling, stereo-frame transport, component-wide marker selection and parity |
| umol-graph resolver | Direct GraphIR-only project, isotope default elision, aromatic/stereo assertion recovery, atomic publication |
| umol-io TableIR and readers | Explicit stereo-bond frames; transient source markers; SMILES/CX/CTfile producers and raise migration; temporary neighbor access |
| umol-graph-core | Iterator-producing connectivity callback, DFS/BFS visitors and collectors, component visitor, semantic renames and usize depth |
| Tests/specification/consumers | Independent correctness cases and properties, existing Python ingestion migration, separate spec updates and bounded performance evidence |

No CTfile writer, coordinate generation, wedge-selection policy from 225, extended stereo,
canonical SMILES, general conjugation API, minimum-marker optimization, or carbon shell-capacity
fix from 166 is included. Existing supported input information must not be silently erased while
adapting a producer. Unrepresentable output must fail explicitly.

### Verification rules for all stages

Every subitem includes exact fixtures and relevant properties with the changed implementation;
unit tests follow test-writing and properties follow docs/development/property-tests.md. State
laws independently of the generated distributions. Assert successful cases as well as errors so
an implementation that rejects everything cannot satisfy roundtrip coverage. Do not weaken laws
or remove difficult fixtures to obtain a green stage.

Use narrow crate checks while editing, then the stage's affected-crate tests and feature suites.
Run cargo fmt --all and git diff --check, and review the full stage diff for unauthorized public
symbols, renames, wrappers, dropped data, and changes in what tests exercise. For workspace-wide
commands or any PyO3 compilation, apply python-build and activate umol-py/.venv with Python 3.13.

Keep benchmark fixtures from the start of algorithm implementation. Record one bounded comparison
at a meaningful stage boundary, including allocations where directly measurable. No broad tuning
campaign or speculative scale assumptions. A specialized IO traversal requires evidence against
the shared graph-core implementation and a separate decision.

### S0 — Executable contracts and source inventory

- **S0a — Representation and API inventory (completed 2026-09-10).** Modules: graph-core traversal/connectivity,
  TableIR stereo/raise, smiles boundaries, resolver, ingest. **Additive (green).** [dep: none]
  Enumerate changed public symbols and their settled signatures, construction guarantees, and
  error ownership. Pin current producers/consumers with workspace searches, including reexports,
  graph views, parser targets, CX annotations, CTfile formats, benchmarks, examples, and Python
  adapters. Check current doc 224 implementation rather than assuming an earlier snapshot.
  Translate every implementation-facing statement in this document into this inventory; do not
  implement a discarded alternative. Verify existing behavior with focused current tests.
- **S0b — Independent semantic fixtures and comparisons (completed 2026-09-10).** Modules: existing IO/graph property
  targets and test support. **Additive (green).** [dep: S0a]
  Preserve the direction-marker examples, stereo/valence scratch findings, and ring-encounter
  regressions as independent expected configurations where relevant. Define boundary equivalence
  excluding spans/spelling but preserving explicit versus implicit H, definite versus Either versus
  absence, labels, and stereo-frame action. Use existing molecular equivalence/remapping machinery
  for independent GraphIR expectations; compare atom states, system members/contributions, and
  stereo explicitly rather than formula or entity counts. Add named successful ordinary, charged, radical, aromatic,
  tetrahedral, and cis/trans cases without claiming unimplemented output support.
- **S0c — Evidence fixtures and independent fuzz follow-up (completed 2026-09-10).** Modules: existing benches and
  umol-io/fuzz. **Additive (green).** [dep: S0a]
  Prepare a bounded fixture set for traversal, raise, projection, marker assignment, and rendering;
  measure existing operations only and add new operations with their implementing stages. Preserve
  parser/raise separation. Check whether the user-run doc 224 fuzz campaign has completed; record
  its evidence or leave it explicitly outstanding. The required campaign is at least 60 minutes
  of fuzz execution using the existing fuzz_parse_opensmiles target, not setup/build time. Do not
  duplicate an active campaign or treat it as a prerequisite for independent implementation work.

**Gate:** Current affected tests remain green; the contract and fixture inventory distinguish
implemented behavior, future assertions, and independent outstanding validation.

### S1 — Shared graph-core traversal, additive APIs

- **S1a — DFS events, callback kernel, and Graph visitor (completed 2026-09-10).** Module:
  umol-graph-core/src/algorithms/traversal.rs. **Additive (green).** [dep: S0a, S0b]
  Add DepthFirstEvent and visit_depth_first with the settled callback contract, bounds, ordered
  candidate roots, suspended neighbor iterators, edge identity tracking, and FinishTree. Add the
  inherent Graph visitor supplying CSR access and the agreed public exports. No TraversalError or
  defensive validation pass. Verify nested events against a small independent DFS, disconnected
  roots, supplied order, loops, parallel edges, isolated nodes, and immediate Break at every event
  kind. Malformed finite callbacks exercise panic-free internal handling without asserting useful
  output. Callback panics/hangs remain caller behavior. Include a long-chain case for no recursion.
- **S1b — BFS events, callback kernel, and Graph visitor (completed 2026-09-10).** Same module.
  **Additive (green).** [dep: S0a, S0b]
  Add BreadthFirstEvent and visit_breadth_first with FIFO discovery, sequential candidate roots,
  usize depth/limit, per-node Finish, and FinishTree before the next root. Verify shortest distances
  independently, maximum-depth boundaries, duplicate roots, disconnected components, exact order,
  and Break with no further callbacks. Document that depth-limited multi-root trees are not
  necessarily full components and are not multi-source BFS.
- **S1c — Event collectors and component visitor (completed 2026-09-10).** Modules: traversal.rs and
  algorithms/connectivity/components.rs. **Additive (green).** [dep: S1a, S1b]
  Add Graph::enumerate_depth_first_events and enumerate_breadth_first_events through their visitors.
  Add visit_connected_components using unrestricted BFS, collecting/sorting one component at
  FinishTree. Make existing enumerate_connected_components collect through this visitor. Remove
  its mislabeled LIFO loop. Check exact visitor/collector agreement, sorted members, component order,
  and early termination before a later component is explored. Compare component membership with
  definition-level reachability, not only with the old implementation.

**Gate:** graph-core unit and proptest targets pass. Add traversal benchmarks on the S0 fixtures;
no specialization or new graph storage is introduced.

### S2 — Neighborhood and size API migration

- **S2a — Neighborhood naming and usize contract (completed 2026-09-10).** Modules: graph-core traversal, refinement,
  reexports, and every affected graph view/consumer. **Breaking (red→green).** [dep: S1b, S1c]
  Rename TraversalAlgorithm to NeighborhoodAlgorithm and retain Bfs. Reimplement neighborhood
  through the shared BFS visitor, returning usize distances and accepting a usize depth limit.
  Update related depth/radius size arguments where they propagate this contract, all callers,
  examples, benchmarks, properties, and exposed bindings in the same stage. NodeId/EdgeId storage
  does not change. Do not retain a compatibility alias or cast through u32. Tests preserve shortest
  distance, ordering, zero limit, disconnected behavior, and refinement results.

**Gate:** graph-core and dependent graph-IR/graph targets compile and pass, including properties;
workspace search finds no stale TraversalAlgorithm or affected u32 depth contract. Reconcile every
changed public signature against S0a, without sweeping unrelated numeric fields into the migration.

### S3 — Explicit TableIR stereo-bond frames and producer migration

- **S3a — Stereo-bond vocabulary and frame algebra (completed 2026-09-10).** Module: umol-io/src/table_ir/stereo.rs
  and its exports. **Additive (green).** [dep: S0a, S0b]
  Add StereoBond { bond, configuration }, BondConfiguration::{Either, Framed { references,
  relation }}, and BondRelation::{SameSide, OppositeSide}. Pin endpoint/reference ordering and
  correspondence with existing permutation/coset operations. Verify one-reference swap, both swaps,
  endpoint exchange, and transport through atom/bond remapping with independent expected frames.
  Preserve StereoAtom, Winding, and LonePair; no generic placeholder or extra Either relation.
- **S3b — Source direction normalization (completed 2026-09-10).** Modules: smiles parser builder and existing annotation
  processing. **Additive (green).** [dep: S3a]
  Implement the cohesive derivation functions with direct unit tests before switching producers.
  Consume all source markers in endpoint viewpoint, including ring-opening/closing spelling and
  shared bonds. Select lowest-index actual references, check same-side redundancy and conflicts,
  and produce no assertion for consistent partial notation. Keep lexical directions in parser
  working state; do not create another persistent authority. Verify all recorded equivalence and
  rejection cases, including conflicts when the other endpoint is unmarked and the doc 224 order.
- **S3c — CTfile and annotation frame derivation (completed 2026-09-10).** Modules: CTfile readers and CX bond annotation
  producers. **Additive (green).** [dep: S3a]
  Reuse supported geometry interpretation to derive definite frames, and preserve site-only Either
  without inventing reference atoms. Treat adjacent wavy-bond effects under the existing convention;
  wedges remain separate. Verify usable 2D/3D, absent/all-zero/degenerate geometry, explicit Either,
  and supported Cis/Trans annotations with exact frames. No coordinate generation or CTfile output.
- **S3d — Publish frames and retire duplicate fields (implemented 2026-09-10; allocation closeout completed 2026-09-11).** Modules: TableIR Molecule/ExtendedMolecule,
  Bond/ExtendedBond, parser targets/conversions, raise, all affected consumers/tests.
  **Breaking (red→green).** [dep: S3a, S3b, S3c]
  Add the stereo_bonds collection consistently to applicable carriers, wire every producer, and
  remove persistent Bond.direction and Bond.stereo plus superseded equivalents. Retain lexical
  marker types only where they still serve parsing. Raise Framed through the settled #C frame
  transport and Either through an undetermined configuration assertion. Preserve out-of-range,
  incidence, and conflicting-annotation failures at their owning boundary. Move source interpretation
  out of raise once frames own its result. Update struct literals, fixtures, conversions, and
  benchmarks together; no consumer may silently prefer one of two operative stereo encodings.
- **S3d1 — Move storage through consuming table conversion (completed 2026-09-11).** Module: table_ir/molecule.rs.
  **Refactor (green).** [dep: S3d]
  Consume ExtendedMolecule atom/bond vectors in the existing fallible conversion; move positions,
  stereo frames, multicenter bonds, comments, and properties into the result. Preserve the existing
  public signature, feature rejection, and complete basic/extended preservation law. Do not add a
  borrowed conversion or compatibility wrapper. Verify success and rejection plus preservation of
  nonempty frame/position collections; inspect ownership and allocations rather than treating value
  equality alone as proof that cloning was removed.
- **S3d2 — Direct CTfile bond construction and borrowed derivation input (completed 2026-09-11).** Modules: CTfile bond
  readers/builders, table_ir/stereo/derive.rs, and CX derivation callers.
  **Refactor (green).** [dep: S3d]
  Have each CTfile bond block accumulate final Bond/ExtendedBond storage and a separate list of
  bond-indexed stereo assertions. Move the bond vector into its table; apply properties before
  deriving configurations. Remove the records vectors from CTfile and CX callers and read the
  original bonds directly. Separate basic/extended derivation functions are acceptable; a shared
  implementation must not require copied tables, dynamic dispatch, or new public adapter types.
  Use bond_stereo_assertions for the code-bearing list. Migrate all affected private callers in
  this subitem. Verify basic/extended V2000, property processing, CX indices, and exact frames.
- **S3d3 — SMILES finalization without full-bond scratch copies (completed 2026-09-11).** Modules: smiles/parser/builder.rs
  and smiles/parser/stereo.rs. **Refactor (green).** [dep: S3d]
  Remove the completed intermediate and lexical copied bond table. Read the pending/final bond
  storage directly during derivation and perform only the necessary ownership transfer to the
  published table. Keep marker participation with transient marker bookkeeping rather than in a
  separately allocated vector covering every bond. Preserve marker values until every sharing
  site has been checked; do not consume a direction on first use. Retain ring-opening table order,
  endpoint viewpoint, atom-stereo incidence order, and CX completion-index remapping. No redesign
  of unrelated parser factoring is included. Verify shared/branched/cyclic markers, partial and
  conflicting notation, ring spelling, and the existing atom/bond remapping laws.
- **S3d4 — Sparse stereo assertion processing and CX frame updates (completed 2026-09-11).** Modules:
  table_ir/stereo/derive.rs and smiles/parser/cx.rs. **Refactor (green).** [dep: S3d2]
  Replace the bond-sized codes array with processing proportional to explicit assertions; group
  assertions by bond using their owned storage and check repeated evidence before consolidation.
  Retain wavy evidence across later wedge overwrites. Merge/check incoming evidence against
  existing frames without expanding every frame back into a code list and rebuilding unchanged
  output. Preserve the reference frame when comparing configurations. Skip stereo work when CX
  changes nothing relevant to it; topology-affecting entries, wedges, and supplied positions must
  be considered explicitly. Geometry-derived frames still cover applicable unannotated double
  bonds, and Either still suppresses geometry at its site. Verify all combinations of existing
  direction frames, repeated codes, Either, geometry, wavy evidence, and reaction-section remapping;
  include labels-only CX input with existing frames. No reinterpretation or acceptance expansion.
- **S3d5 — Bounded local stereo scratch (completed 2026-09-11).** Modules: both stereo derivation kernels and
  table_ir/raise.rs. **Refactor (green).** [dep: S3d3, S3d4]
  Replace per-endpoint substituent vectors and raise's two-element blocks vector with bounded
  local storage for the supported two-substituent-per-endpoint domain. Check site relevance before
  gathering ligands. Detect excess distinct incidences without allocating to accommodate them;
  preserve duplicate-incidence treatment, canonical reference selection, unsupported-site checks,
  and frame transport. Verify zero/one/two/excess substituents, shared ligands, cumulated axes,
  and single/both reference swaps. Persistent StereoAtom ligand vectors and wedge-algorithm
  redesign are not part of this scratch-storage change.
- **S3d6 — Apply bond frames without a constraint HashMap (completed 2026-09-11).** Module: table_ir/raise.rs.
  **Refactor (green).** [dep: S3d5]
  Order borrowed frame references by their table-bond index and consume them alongside the bond
  construction loop. Construct each constraint at its destination instead of retaining all
  constraints in a HashMap. Do not sort/mutate the input collection or substitute a dense
  all-bond mapping. Preserve arbitrary input frame order, duplicate-site and contextual failures,
  and the localized/dative/noncovalent partition. Verify unsorted frames, duplicate and invalid
  sites, and relation bonds preceding/between framed localized bonds. Allocation scales with the
  frame ordering actually required, not every table bond.
- **S3d7 — Initial SMILES table capacity (completed 2026-09-11).** Modules: both SMILES inner parsers and builder
  initialization. **Refactor (green).** [dep: S3d3]
  Retain one-pass parsing and grow atom/bond tables as entries are produced; do not add a
  preliminary counting pass. Initial reservation is a tunable implementation estimate, not a
  requirement to start at zero capacity. Reduce the current overestimation, including reservations
  inflated by CX text and later reaction sections. Preserve parsing and diagnostic semantics;
  check long annotations, bracket-heavy input, and uneven reaction sections using inline examples.
  Capacity estimates must not become acceptance limits. Further tuning does not block the other
  allocation corrections or require an open-ended benchmark campaign.
- **S3e — Operation-local neighbor lookup (completed 2026-09-11).** Modules: TableIR raise and parser normalization.
  **Additive/refactor (green).** [dep: S3d3, S3d4, S3d6]
  Remove unconditional adjacency construction where frame-based raise needs none; derive temporary
  lookup only for remaining operations that require it. Preserve all-table-bond indices and the
  partition into localized edges versus relations. TableIR remains tables. Use iterator adaptation
  to graph-core Neighbor for traversal rather than adding Graph storage or an adjacency trait.
  Verify ordinary/stereo/mixed-bond raising and malformed-index behavior; use the bounded raise
  fixtures to report allocation changes without opening a benchmark-driven redesign.
  This remains the owner of the audit's neighbor-allocation work: skip lookup where unused and
  share one temporary lookup across the phases that need the same unchanged table. Preserve the
  public AtomNeighbors/Neighbor surface and table order; contiguous private storage may replace
  the nested vectors without introducing a foundational adjacency type or storing a graph in
  TableIR. Final Graph adjacency cannot generally replace the all-table incidence index.

The allocation follow-ups preserve the existing representation and source semantics; each closes
with its affected tests green. S3d1 and S3d2 are independent of S3d3. The dependency paths converge
at S3d5, then S3d6 and S3e; S3d7 depends only on S3d3. All are complete.
The allocation policies are settled; private factoring and capacity tuning remain implementation
choices within the recorded constraints.

After each allocation subitem, assess whether its measurements support the rationale, what
concrete costs remain for downstream improvements, and whether the change should be retained.
Separate measured allocation savings from runtime evidence and future expectations; mixed timings
are not grounds for an open-ended local tuning campaign.

For each affected path, review the resulting allocation/ownership flow and reuse the existing
bounded inline fixtures and allocation probe. Account for in-place reuse by owning Vec collections;
counting collect expressions is not an allocation audit. Verify removal of identified copies and
unnecessary allocations, retain the established semantic/property suites, and run the S3 gate.
Measurements provide evidence of the final change, not a search for the design or an invitation to
extend the benchmark campaign. No fixture files are introduced into benchmarks.

**Gate:** IO unit/property/CTfile suites and dependent graph interpretation tests pass. Every former
raw-direction/stereo-field producer and consumer is accounted for. Parser acceptance and error-layer
changes are reviewed for semantic preservation, not accepted solely because tests were rewritten.

### S4 — Direct GraphIR-only projection

S4a0 grouped three isotope-repair prerequisites, completed before S4a. The draft's default-registry
isotope rejection was repaired rather than retained as a support boundary.
S4a0a's registry change exposed its dependency on isotope resolution and was not an independent
green milestone. S4a0c has restored the conformance gate through the isotope phase with Natural
resolver policy; the test harness does not fill isotope fields, and the corpus inputs and
snapshots remain unchanged.
The phase order, Strict/Natural behaviors, Strict as the general default, Natural for SMILES, and
the registry's isotope-Undetermined invariant are settled. No general DSL conversion change is
required.

- **S4a0a — Secure the atom-type registry's isotope invariant.** Modules: ops/valence/registry,
  its Rust callers/macros, and umol-py/model/valence.
  **Breaking (green milestone with S4a0c); completed 2026-09-11.** [dep: S0a]
  Require isotope Undetermined in every stored row. Set registry raise's isotope default to
  Required while retaining its other concrete and valence-constraint defaults. Cover from_atoms,
  add, TOML string/file loading, the registry macro, built-in/default construction, and Python
  producers. TOML currently builds storage directly, so securing add alone is insufficient.
  Establish one authoritative entry check at construction/insertion; preserve the invariant through
  cloning and read-only access without checking again in admission or lookup. Add try_from_atoms
  and try_add returning ConfigError; from_atoms and add assert the same checks. Fallible loaders
  and Python report invalid input as errors. Reject before changing contents or content_hash.
  Test all producer routes, acceptance of Undetermined, rejection of Natural, masses, sets,
  variables and other non-Undetermined forms, failed-insertion atomicity, and built-in row/hash
  consistency. Migrate callers and fixtures that previously inherited or explicitly asserted
  Natural, preserving what their tests exercise. Rewrite duplicate-admission regressions using
  admissible non-isotope overlaps. Leave general DSL defaults and IntoIr/FromIr unchanged.
- **S4a0b — Standalone isotope resolver and policies.** Module: ops/resolve/isotope.rs and
  the phase's configuration/error vocabulary. **Additive (green); completed 2026-09-11.** [dep: S0a]
  Implement IsotopeResolver with independent Strict and Natural policies, defaulting to Strict.
  Strict preserves an
  unresolved isotope; Natural fills only Undetermined. Preserve supplied Natural, masses, sets,
  and variables. Report unresolved isotope information without requiring valence completion or
  modifying unrelated atom fields. Follow the existing phase planning and atomic application
  contracts, with the exact public symbols reconciled before implementation.
  Implement project under that same policy: Natural becomes Undetermined only under Natural;
  Strict retains Natural, and both retain explicit masses. No separate projection setting.
  Test exact edits/results, explicit-value preservation, policy differences, idempotence, and
  failure/underdetermination behavior. Add inverse-law properties for resolve(project(source))
  under each policy over resolved isotope values, preserving unrelated fields. Do not wire it
  into the composite pipeline yet.
- **S4a0c — Integrate isotope resolution and remove valence defaults.** Modules: composite
  resolve/config, valence counts/registry/admission consumers, and umol-py resolve/model bindings.
  **Breaking/refactor (red→green within subitem); completed 2026-09-11.** [dep: S4a0a, S4a0b]
  Run isotope resolution after assertion placement and before valence admission; aromaticity
  remains directly after valence. Expose the independent isotope setting through ResolveConfig
  and its existing Python configuration surface. Keep Strict as the general default and select
  Natural in the molecular and reaction SMILES default entry points; preserve supplied settings.
  Remove counts admission's isotope completion. ValenceTable has no AtomForms and needs no
  isotope validation. Atom typing consumes the registry invariant established in S4a0a without
  revalidating entries; molecule AtomForms may still carry explicit isotope information.
  Audit phase progress so an unresolved isotope does not prevent determination of valence-owned
  fields; retain final composite completeness and atomic publication. Test natural and explicit
  isotope input with both valence sources, Strict/Natural crossed with both valence tie-breaks,
  isotope preservation through registry meets and counts completion, stereo ordering, and
  composite nonpublication when isotope remains unresolved. Update Rust/Python docs and tests
  together. Apply the staged isotope-composition spec correction in this subitem, alongside the
  SMILES Natural preset; the remaining specification updates retain their later stage. Use a
  bounded inline benchmark comparison for default resolve/ingest paths; no fixture
  files or tuning campaign. Revise the S4a draft tests that encoded isotope rejection as a limitation.

S4a0a and S4a0b are independent foundations; both precede S4a0c, which precedes S4a. None is
deferrable for S4a completion. Each subitem includes its own tests and caller migration; run the
affected graph-IR/IO/graph suites and Python checks where bindings change. The S4a0 integration gate
includes graph conformance and property suites for both valence sources and both isotope policies.

S4a0a implementation checkpoint (2026-09-11): try_from_atoms and try_add enforce the registry
invariant, with asserted counterparts and one shared entry check for programmatic and TOML paths.
Python delegates to checked Rust construction. Registry raise leaves isotope Undetermined; failed
insertion preserves contents and hash. The shared check also retains Python's charge-range
rejection for all producers. No isotope checks were added to lookup or admission, and no general
DSL defaults, conversion traits, or counts behavior changed.

Focused verification: 1,186 graph library tests and the whitepaper test pass; all seven graph
property tests pass, including projection with carbon-13 under both valence sources. Rust binding
tests pass (1,642 passed, two ignored); after rebuilding with Python 3.13, all 145 Python model
tests pass. Graph/Python all-target clippy passes with warnings denied. The complete Python suite
has 1,523 passed, two skipped, and four failures in stereo-error expectations for F/C=C and
conflicting directional markers. Those parser and Python error paths are unchanged from HEAD;
their error-layer expectations are separate from the registry work.

At the S4a0a checkpoint, the full graph conformance run had 23 passed and 660 snapshot failures
after removal of the registry's isotope default. The corpus and snapshots remained unchanged;
the gate awaited the explicit isotope phase and Natural resolver policy in S4a0c.

S4a0b completion (2026-09-11): ops::resolve::isotope defines IsotopePolicy and IsotopeResolver,
with new(policy), Default (Strict), plan, resolve, and project. IsotopeContradiction is empty;
IsotopeError carries transaction failures, and IsotopeProjectError additionally identifies the
first atom with a non-ground isotope. The resolve facade re-exports the resolver, policy, and
resolution errors; the projection error remains module-qualified, as for valence projection.

Planning retains proposed Natural defaults even when another atom's isotope is unresolved.
Direct resolve publishes only determined plans. Project requires ground isotope fields but does
not require completed valence or other fields. It elides Natural only under Natural and retains
explicit masses under both policies. Sets and variables are preserved without normalization or
selection. Both operations apply edits atomically and skip the editor's molecule clone when no
edits are needed. No composite configuration, valence behavior, or format defaults changed.

Verification: all 52 isotope unit cases pass within the general module, including exact plans,
stale-plan transaction rollback, partial-result nonpublication, and projection failures. All
1,237 graph library tests and all nine graph properties pass (256 cases per property); the two new
properties cover exact projection/recovery under both policies and resolution idempotence with
plan/application agreement. Graph all-target clippy with conformance and proptest enabled passes
with warnings denied; formatting and diff checks pass.

S4a0c completion (2026-09-11): ResolveConfig.isotope and Resolver.isotope integrate the standalone
phase before valence admission. Placement edits affect assertions only, so its editor also applies
the isotope plan against the unchanged inherent fields; no extra molecule clone is introduced.
Both determined and underdetermined isotope plans continue through the remaining phases. The final
completeness check and atomic publication remain intact. ResolveError::Isotope reports edit
application failure. Counts admission no longer supplies Natural, and both candidate sources
preserve unresolved isotope forms while determining valence fields. The completion documentation
now explicitly distinguishes phase-owned fields from still-unresolved fields.

The Python surface includes IsotopePolicy and ResolveConfig's isotope keyword, getter, equality,
and representation. General defaults remain Strict; Rust string/byte and Python molecule/reaction
SMILES convenience entry points select Natural when configuration is omitted. Explicit configs
remain authoritative. TableIR raise already supplies Natural for an omitted SMILES isotope, so
Strict preserves that format-owned value too; it does not undo input interpretation. The isotope
composition correction is applied to the SMILES specification. No general DSL defaults or
conversions changed.

Verification: the combined graph-IR/IO/graph/Python Rust gate passes, including 1,295 graph,
6,816 graph-IR, 3,916 IO, and 1,646 binding unit tests. All 683 resolution conformance cases pass
with unchanged snapshots; all ten graph properties pass with 256 cases each. New cases cover
both isotope policies, both valence sources and tie-breaks, exact isotope preservation, partial
nonpublication, and reaching stereo contradictions despite unresolved isotopes. The composite
property checks exact results and idempotence across those configurations. Graph/Python all-target
clippy passes with warnings denied. The rebuilt Python 3.13 suite has 1,541 passed, two skipped,
and the same four previously recorded stereo-error expectation failures; those expectations and
parser/error paths are unchanged by S4a0c.

The bounded Criterion comparison against 24d9b51a8 uses the existing inline cases, ten samples,
0.5 s warm-up and 1 s measurement per case. All seven ingest and sixteen isolated resolve
measurements completed before and after integration. Ingest point estimates increased 0.5–4.1%; isolated resolve changes
ranged from −1.5% to +6.1%. These results cover SMILES-raised inputs, whose isotope fields are
already Natural or explicit; no tuning campaign was undertaken. Selected point estimates in
microseconds:

| Case | Before | After | Change |
| --- | ---: | ---: | ---: |
| Ingest methane | 3.291 | 3.425 | +4.1% |
| Ingest octane | 15.226 | 15.295 | +0.5% |
| Ingest benzene | 22.328 | 22.685 | +1.6% |
| Ingest purine | 59.263 | 59.583 | +0.5% |
| Resolve chain of 64 carbons | 98.781 | 100.970 | +2.2% |
| Resolve branched molecule | 11.344 | 12.036 | +6.1% |
| Resolve shared-marker cycle | 33.340 | 33.500 | +0.5% |

Both timing processes subsequently encountered the existing atom-typing projection benchmark's
unresolved-source setup assertion. Its setup now obtains a resolved source through ingest_smiles
before either candidate-source projection; preparation remains outside timing. The complete
benchmark target passes in Criterion test mode. Logs and the full comparison are in
scratch/s4a0c-before.log, scratch/s4a0c-after.log, and scratch/s4a0c-benchmark.csv.

- **S4a — Correct ordinary atom projection.** Module: resolve/valence.rs and its projection
  callers/tests/benchmarks. **Breaking/refactor (red→green); correction completed 2026-09-11.** [dep: S0a, S0b, S4a0c]
  Remove candidate re-admission, tie-breaking, state comparison, and field clearing performed only
  for reconstruction. Retain direct GraphIR field handling under the settled capacity rule;
  actual TableIR narrowing belongs to Convey and notation support to the SMILES boundary. No
  ordinary atom projection operation remains; remove ValenceResolver::project and its error type
  rather than add an identity stub. Do not replace reconstruction with unchecked clearing.
  Remove projection-only dependencies and errors
  that exist solely for admission ambiguity or mismatch; migrate affected callers in this subitem.
  This code knows no TableIR or SmilesIoConfig. Retain forward atom-typing/counts behavior, the
  registry isotope invariant, exact post-meet duplicate elimination, and isotope default elision.
  Preserve useful ordinary atom/charge/radical cases as independent forward-resolution and
  supported-input fixtures. Test exact preservation in the aromatic caller's no-system path.
  Replace arbitrary-state inverse rejection expectations with the corrected contract; field
  lowering and its errors are tested with the actual conversion in S7a.
- **S4b — Aromatic assertion recovery.** Modules: resolve/aromaticity.rs and graph-core induced edges.
  **Breaking/refactor (red→green); correction completed 2026-09-11.** [dep: S4a]
  Replace reconstruction with direct atom `#a<n>` and bond #a assertions from existing systems,
  then remove those systems. Preserve atom electron fields and localized links. Remove the
  ValenceResolver/ValenceTieBreak dependency, candidate
  admission, joint selection/ring discovery, and recovered-state comparison. Remove errors and
  blanket assertion checks justified only by the arbitrary-GraphIR inverse. Migrate the public
  signature, callers, tests, and benchmarks together. Reject nonzero system charge, #u > 0, and
  non-singlet spin; do not localize them or treat undetermined required values as defaults.
  Test exact `#a<n>`/#a output and retained fields for benzene, heteroaromatics, localized links,
  fused systems, and atom-localized charge. Retain independent input expectations and frame-order
  variations. Remove required rejection cases based solely on alternative partitions, candidate
  ambiguity, scope/selection differences, or search bounds.
- **S4c — Stereo assertion recovery.** Module: resolve/stereo.rs, reusing ops/stereo.rs and
  graph-IR coset transport. **Additive (green).** [dep: S4a]
  Recover #T/#C assertions in the frames that the same model derives, remove entities after
  transporting their configurations, and preserve typed virtual participants. Test reordered
  actual ligands, implicit versus explicit H, LP frames, endpoint flips, Either where in domain,
  required frame-derivation failures, and preservation under permissive resolver policies.
  Unsupported extended kinds fail rather than partially projecting. No runtime re-resolution
  or comparison of recovered configurations is added.
- **S4d — Atomic Resolver::project composition.** Module: resolve.rs and related
  error/report vocabulary. **Additive (green).** [dep: S4a, S4b, S4c]
  Compose the projection plans with resolve's Result/Solution and commit semantics. Reject nonzero
  charge or #u > 0 on bonds and aromatic systems. Preserve non-elidable information on other
  unsupported structures for explicit failure. Only successful projection publishes the candidate;
  every other outcome leaves the caller unchanged. Do not resolve a copy or compare a recovered
  molecule. Test project(resolve(input)) under the agreed equivalence for supported raised foreign
  inputs, direct field/frame expectations, and atomic failure independently. Reverse-direction
  experiments are optional evidence, not a gate for arbitrary GraphIR.

#### S4a correction completion — 2026-09-11

ValenceResolver::project and ValenceProjectError are removed. This eliminates the private molecule
clone, whole-atom field clearing, candidate admission, tie-breaking, and exact recovery comparison
from ordinary valence projection. No replacement phase or empty method is introduced. Ordinary
atom values remain for Convey's later TableIR conversion; isotope default elision remains with
IsotopeResolver. Forward resolution and the registry invariant/deduplication implementation are
unchanged.

The only production caller was AromaticityResolver::project. Its no-system path now returns
Determined with an empty report and preserves the entire molecule, including non-ground fields.
The six diagnostics still used by its existing system reconstruction moved from the removed
ValenceProjectError wrapper directly into AromaticityProjectError: NonConcreteAtom,
NonLiteralBondOrder, DativeBonds, MulticenterBonds, IncompleteAtom, and AtomMismatch. This is caller
migration, not completion of S4b: aromatic candidate admission, ring discovery/selection, and
recovery comparisons still await removal there. No other public consumers or Python bindings used
the removed valence projection surface.

Fourteen independent atom expectations and the custom-registry overlap regression now exercise
Resolver::resolve directly. Seven molecular fixtures now start from SMILES with fixed bracket H
counts, preserving the original determined-input domain under both valence sources and tie-breaks.
Their complete expected molecules cover branches, a ring, nitrile, amide, sulfoxide, zwitterion,
and disconnected charged/isotopic components. Carbon/radical/isotope and mixed-element generated
cases retain their independent expected states as forward-completion properties. The obsolete
projection mismatch, candidate ambiguity, and pairing-rejection tests are removed. Existing
isotope, registry, electron-pairing, and aromatic projection coverage remains.

Verification: PROPTEST_CASES=256 cargo test -p umol-graph --features conformance,proptest --offline
passes, including 1,411 unit cases, 683 conformance cases, 11 property tests, and integration tests.
All-target graph Clippy with the same features and -D warnings passes. Logs:
scratch/s4a-correction-gate.log and scratch/s4a-correction-clippy.log. Formatting was run, with
unrelated formatter-only changes removed from the diff. The remaining resolve benchmark target
passes in Criterion test mode (scratch/s4a-correction-bench.log); git diff --check passes.

The valence projection and projection-plus-resolve benchmark groups are removed with the operation.
There is no replacement operation to time or new speedup estimate to report. Retain the historical
costs below; the next bounded timing comparison belongs to S4b's actual aromatic transformation.

#### S4b correction completion — 2026-09-11

AromaticityResolver::project now takes only the mutable molecule and returns
Result<Solution<(), AromaticityContradiction>, AromaticityProjectError>. It reads each existing
system's member-aligned contributions into atom `#a<n>` assertions, marks the system's induced
bonds #a, and removes the system. Atom electron fields and links between systems are preserved.
Edits are published together; errors preserve the caller's molecule. System removal uses the
editor's existing constraint-compaction semantics.

Projection no longer performs valence admission, candidate selection, ring discovery, or a recovery
comparison. It does not depend on valence policy or aromatic perception scope. The remaining
diagnostics are NonConcreteSystem, ChargedSystem, SystemSpin, AtomAssertion, and BondAssertion:
required contributions/charge/spin must be concrete, system charge must be zero, system spin must
be closed-shell singlet, and existing atom/bond assertions must admit the projected assertions.
The temporary reconstruction diagnostics described in the S4a completion record are removed.
Forward aromatic resolution methods and helpers are unchanged from the S4a closeout.

Exact fixtures cover contributions 0, 1, and 2, preserved lone-pair fields, localized charge,
heteroaromatics, fused rings, and links between separate systems. Identity cases include non-ground
ordinary atoms. Failure cases check exact errors and atomic publication. The property suite checks
direct mapping, idempotence, rotated/reversed member frames with transported contributions, and
projection of independently expected SMILES inputs under both valence sources. It does not enforce
arbitrary GraphIR reconstruction.

The frame-order property exposed the sorted-input restriction in Graph::induced_edges: the valid
benzene frame [1,2,3,4,5,0] omitted two bonds. The corrected method implements endpoint membership
for arbitrary node order and repetitions, yielding each edge ID once, including self-loops and
parallel edges. Its signature is unchanged. It builds a subset-sized HashSet and scans stored
edges once; expected work is O(|nodes| + |E|), with no graph-sized membership allocation. This
also corrects aromatic-system bond views without a local workaround. Exact core regressions and
a generated multigraph property compare against the endpoint definition. Kekulizer plan tests
now compare exact matched/unmatched bond sets independently of unspecified iteration order,
while retaining repeatability checks; no kekulization implementation changes were needed.

Inline core benchmarks include reversed node lists. With 10 samples, 0.5 s warmup, and 1 s
measurement, full rings of 6, 64, and 1,024 nodes take 0.2020, 1.6827, and 26.660 microseconds.
Selecting six adjacent nodes from a 4,096-node ring takes 31.837 microseconds: the whole-edge
scan remains visible for small subsets of large graphs. These are current implementation costs,
not before/after comparisons. Log: scratch/induced-edges-bench.log.

Verification: 1,013 graph-core unit cases, the generated induced-edge property and related graph
properties at PROPTEST_CASES=256, and 6,816 graph-IR unit cases pass (three existing ignored cases).
The graph gate with features conformance,proptest and PROPTEST_CASES=256 passes: 1,355 unit cases,
683 conformance cases, 11 property tests, and integration tests. All-target Clippy for graph-core
and graph with those features and -D warnings passes. Formatting and git diff --check pass.
Logs: scratch/induced-edges-{core-gate,property-gate,ir-gate}.log,
scratch/s4b-correction-gate.log, and scratch/s4b-correction-clippy.log.

The inline aromatic projection benchmark retains the four prior inputs and excludes parsing,
source construction, and caller cloning. Projection now has one policy-independent measurement;
the projection-plus-resolve benchmark is removed. Central estimates in microseconds, using
10 samples, 0.5 s warmup, and 1 s measurement:

| Input | Previous counts project | Previous atom-typing project | Direct project |
| --- | ---: | ---: | ---: |
| Benzene | 12.487 | 35.008 | 1.6204 |
| Pyrrole | 11.873 | 28.182 | 1.4085 |
| Naphthalene | 23.632 | 61.803 | 2.3157 |
| Biphenyl | 26.002 | 72.129 | 2.9099 |

These bounded results support the expected reduction from removing reconstruction; they do not
establish scaling for large collections of systems. Full output:
scratch/s4b-correction-benchmark.log. S4c is next.

#### Historical S4a/S4b implementation and measurements — 2026-09-11

The following records describe the implementation before the contract correction. Its runtime
reconstruction checks and the tests enforcing them are superseded by the corrected subitems above;
passing those checks does not complete the corrected work. Retain the measurements for a bounded
before/after comparison. Projection-plus-resolve measurements are historical diagnostics, not a
required output-pipeline benchmark or deliverable.

S4a prior implementation (2026-09-11): ValenceResolver::project retains element, isotope, charge, fixed
implicit H, bonds, entities, and assertions. It opens lone pairs and both unpaired-electron fields
on a private molecule, admits that reduced input, applies the supplied existing tie-break, and
publishes only when the selected inherent atom fields exactly match the source. Candidate
constraints remain admission evidence; projection does not copy them into the output. Exact
post-meet duplicate elimination remains in the atom-typing producer. No registry pruning API or
additional production change was needed after S4a0.

The phase requires concrete atoms, literal localized bond orders, and no aromatic evidence,
dative bonds, or multicenter bonds. Non-concrete atoms, unsupported structures, incomplete
admissions, and selected-state mismatches have distinct ValenceProjectError variants. Admission
contradictions and unresolved candidate plurality retain Solution semantics. Every unsuccessful
outcome leaves the caller unchanged. This is the ordinary valence phase: other entities are
preserved without interpretation, and the full molecule inverse and bond/system export restrictions
remain S4b–S4d.

Exact fixtures cover ordinary main-group atoms, charges, radicals, Natural and explicit isotope
masses, implicit versus explicit H, branches, rings, nitriles, amides, sulfoxides, zwitterions, and
disconnected components. Custom registries exercise overlapping rows, strict ambiguity,
MostSaturated selection, unreconstructible pairing/multiplicity, and incomplete candidates.
Generated inverse tests retain the carbon/radical/isotope domain and add C/N/O/S chains with
localized charges, lone pairs, single/double bonds, and disconnected components, checking both
candidate sources and both tie-breaks. Expected source states are constructed independently of
resolution. Tests remain in the general valence test module and existing project property module.

Verification: graph unit, integration, conformance, and property suites pass with
PROPTEST_CASES=256 and features conformance,proptest, including 1,326 unit cases, 683 resolution
conformance cases, and 11 property tests. Graph all-target Clippy with those features and
-D warnings passes. Logs: scratch/s4a-gate.log and scratch/s4a-clippy.log.

The existing inline projection benchmarks ran with 10 samples, 0.5 s warmup, and 1 s measurement.
Times below are central estimates in microseconds; these establish current costs, not an
improvement claim. Projection includes its private copy and candidate admission; the second
measurement also includes deliberate full re-resolution. Input construction and the caller's
input clone are outside timing. Full results: scratch/s4a-benchmark.log.

| Input | Counts project | Counts project + resolve | Atom typing project | Atom typing project + resolve |
| --- | ---: | ---: | ---: | ---: |
| Methane | 0.393 | 1.991 | 3.740 | 8.804 |
| Chain, 64 atoms | 22.841 | 74.987 | 244.420 | 524.240 |
| Methyl radical | 0.431 | 2.018 | 3.794 | 9.427 |
| Ammonium | 0.401 | 2.030 | 2.121 | 5.911 |

S4b prior implementation (2026-09-11): AromaticityResolver::project takes the mutable molecule, a borrowed
ValenceResolver, and the existing ValenceTieBreak. Its result is
Result<Solution<ResolveReport, ResolveContradiction>, AromaticityProjectError>. The new error
type stays in resolve::aromaticity; no facade re-export, registry API, or boundary type is added.

On a private candidate, project sets atom lone-pair count, unpaired-electron count, and
multiplicity to Undetermined, meets generic aromatic assertions onto system atoms and their
induced bonds, and removes the systems. Element, isotope, charge, fixed implicit H, bond orders,
and unrelated entities remain. Clearing the electron fields tests whether retained information
reconstructs the source atom state; retaining them would supply part of the answer to admission.
Compatible existing aromatic assertions retain their more specific information. An incompatible
assertion causes an error. Marking only each system's induced bonds preserves localized links
between distinct systems. A molecule without systems delegates to ordinary valence projection.

The supplied valence source admits the reduced atoms, then the existing joint selector runs
without stored systems constraining its answer. The remaining valence tie-break follows the
forward resolver's rule. Success requires exact inherent atom recovery and equality of system
membership, atom-aligned electron contributions, charge, and spin, ignoring system ids and
participant order. Every other outcome preserves the caller. Candidate constraints remain solver
evidence. Atom changes share the editor's mutable atom storage rather than rebuilding the full
atom array after assertion insertion.

Unsupported input includes non-concrete atoms/systems, non-literal bond orders, dative or
multicenter bonds, nonzero system charge, and non-singlet system spin. System-local assertions
without a projected representation are rejected, as are molecule assertions that system deletion
would remove. Incomplete candidate fields and changed reconstructed states are errors; ambiguity,
the existing assignment bound, and chemistry contradictions retain Solution semantics. This is
the constitution-phase inverse check. Stereo recovery, full resolver equivalence, and the general
bond charge/spin export restrictions remain S4c–S4d.

Verification covers benzene, pyridine, pyrrole, furan, atom-localized anionic/cationic rings,
naphthalene, and linked benzene rings under both valence sources, both valence tie-breaks, and both
aromaticity tie-breaks. Exact output verifies generic atom/bond assertions and removal of system
entities. A connected 18-atom case rejects a source partition into three six-atom systems when
selection reconstructs one 18-atom system; the whole-system source succeeds. Failure tests cover
charge, spin, contribution changes, assertion loss/conflicts, scope exclusion even with Keep,
custom-registry ambiguity, and the assignment limit. Generated independent ring states verify
the inverse with rotated/reversed member frames and nonuniform electron contributions.

Graph unit, integration, conformance, and property suites pass with PROPTEST_CASES=256 and
features conformance,proptest: 1,408 unit cases, 683 resolution conformance cases, and 12 property
tests. Graph all-target Clippy with these features and -D warnings passes. Logs are in
scratch/s4b-gate.log and scratch/s4b-clippy.log.

The new inline aromaticity projection benchmark group measures candidate construction, valence
admission, joint selection, and recovery comparison. Its paired measurement adds full
re-resolution; parsing, source construction, and caller input cloning are outside timing.
Central estimates in microseconds from 10 samples, 0.5 s warmup, and 1 s measurement:

| Input | Counts project | Counts project + resolve | Atom typing project | Atom typing project + resolve |
| --- | ---: | ---: | ---: | ---: |
| Benzene | 12.487 | 31.546 | 35.008 | 77.069 |
| Pyrrole | 11.873 | 29.123 | 28.182 | 63.196 |
| Naphthalene | 23.632 | 59.901 | 61.803 | 135.560 |
| Biphenyl | 26.002 | 64.160 | 72.129 | 155.510 |

These establish initial costs, not a before/after performance claim. All benchmark input checks
and measured cases pass. Full output: scratch/s4b-benchmark.log.

**Corrected S4 gate:** graph unit, conformance, and property suites pass for both valence strategies.
Record the supported input domain and concrete unsupported values. Verify direct projection and
atomic publication without runtime reconstruction. Run a bounded projection benchmark comparison
against the historical measurements; no re-resolution measurement or tuning campaign is required.

### S5 — IO traversal and direction assignment kernels

- **S5a — Output traversal.** Module: smiles writer implementation within umol-io.
  **Additive (green).** [dep: S1a, S1c, S3d, S3e]
  Adapt temporary bond-row incidence to graph-core DFS. Derive ordered roots, parent/child edges,
  non-tree edges, ring encounters, and emitted atom/bond order without storing a Graph in TableIR.
  Pin deterministic TableIR ordering and ring-label allocation with exact fixtures. Retain the
  traversal only within the operation. Verify branches, disconnected/empty molecules, rings,
  multiple closures at one atom, and reference transport under changed traversal.
- **S5b — Marker-assignment components.** Module: smiles writer stereo functions.
  **Additive (green).** [dep: S3a, S3d, S1c]
  Implement structural site/candidate eligibility and auxiliary connectivity, retaining per-endpoint
  candidate membership and table bond IDs. Include unspecified coupling sites, branches and cycles;
  retain components requiring definite output. Definite sites with no candidate on an endpoint fail.
  Check membership against a direct connectivity definition, including separate systems, shared
  candidates, localized links, rings, Either, and rejected unsupported side-bond types.
- **S5c — Boolean selection and parity assignment.** Same module.
  **Additive (green).** [dep: S5a, S5b]
  Implement required endpoint coverage, forbidden double coverage at unspecified sites, forced
  propagation, lowest-index branching with unmarked first, and complete backtracking. Check parity
  after a complete selection; seed each selected constraint group with `/` at its lowest-index
  marker in output orientation. On conflict continue search; only exhaustive failure means no
  assignment. Use a small exhaustive marker/sign oracle independent of the production propagation
  to test soundness and completeness, plus chain/branch/cycle conflicts and alternative solutions.
  Do not minimize marker count or adopt the reviewed incomplete cleanup/conservative-ban policies.
- **S5d — Token and atom-stereo formatting.** Module: smiles writer.
  **Additive (green).** [dep: S5a, S5c]
  Emit atom/bracket fields, preserved implicit-H counts, isotope/charge/class, bonds, aromatic
  tokens, components, branches, and ring labels. Transport StereoAtom into emitted encounter order
  before choosing @/@@; use assigned directions with endpoint/viewpoint correction at closures.
  Explicit H nodes remain explicit. Ordinary non-aromatic GraphIR convey may use bracket H counts;
  parsed boundary normalization must preserve its separate H-ownership semantics. Test exact
  spellings where policy determines them and independent configurations where spelling is flexible.

**Gate:** kernel tests and exhaustive bounded assignment comparison pass; emitted fixtures parse
under the intended IO configuration and preserve independent stereo expectations. Formatting and
marker-search benchmarks use S0 fixtures without making optimization a gate to semantic progress.

### S6 — Checked SMILES boundaries and rendering

- **S6a — Smiles construction/rendering contract.** Modules: smiles/molecule.rs, config.rs,
  error.rs, writer. **Additive (green).** [dep: S5d]
  Add checked from_table_ir(table, config), render, and render_with. Keep private payloads and
  existing read/consume accessors. Validate the supplied open table's required integrity and
  establish actual format representability, including assignment feasibility, at construction.
  Share the operative checks/formatting logic with rendering; do not cache traversal or markers.
  Return distinct construction/render diagnostics backed by shared representability reasons.
  Test successful construction implies same-config rendering, narrower-config failures, and
  preservation of parsed boundary semantics. Unsupported Either/annotations must fail in a
  notation that cannot express them, never become mere absence. Do not add Display with weaker
  guarantees or an arbitrary-parts shortcut.
- **S6b — ReactionSmiles construction/rendering.** Module: smiles/reaction.rs and shared writer.
  **Additive (green).** [dep: S6a]
  Add the corresponding checked constructor and render methods. Compose three ordered sections
  and molecular components. Preserve parsed Atom.class values, repeated classes, and agents;
  establish consistency of the supplied table's derived mapping index at the constructor boundary.
  Test empty sections, multiple components, repeated labels, agent-only labels, molecular errors
  with section context, and the same-config construction guarantee.

**Gate:** molecular/reaction boundary roundtrip properties and normalization idempotence pass under
supported IO configurations. No public direct TableIR formatter or stored marker assignment exists.
A parser/writer agreement alone is not enough: independent frame fixtures from S0 remain required.

### S7 — Convey and text export

- **S7a — Convey trait and molecular conversion.** Module: new umol-graph export module and
  its owning conversion functions. **Additive (green).** [dep: S4d, S6a]
  Add Convey with associated Input/Config/Error and implement it for Smiles. Clone source GraphIR,
  invoke Resolver::project, and convert the projected GraphIR fields/assertions into TableIR.
  Build StereoAtom/StereoBond frames here, not in the resolver. Lower projected values that TableIR
  can carry and reject concrete unsupported values/ranges. The checked boundary constructor owns
  SMILES notation/configuration support. Do not add reconstruction checks at either boundary.
  Test the supported foreign-input roundtrip through convey/render/parse under the agreed
  equivalences, with required-success fixtures. Include original-input immutability on success
  and failure and exact diagnostics for unsupported encodings.
- **S7b — Reaction convey.** Same module, ReactionSmiles implementation.
  **Additive (green).** [dep: S7a, S6b]
  Materialize reaction sides through existing ReactionSpan/correspondence operations, then convey
  both sides. Assign paired atoms label i + 1 in existing matched-pair order; derive atom_mapping
  from Atom.class and leave the agent section empty for GraphIR export. Do not infer correspondence
  or recover discarded source labels. Preserve the correspondence through any entity remapping.
  Test creation/deletion, nontrivial pair ordering, disconnected sides, materialization failure,
  side-specific convey errors, and no accidental preservation claim for source-only metadata.
- **S7c — Export conveniences and composed errors.** Same module and exports.
  **Additive (green).** [dep: S7a, S7b]
  Add export_smiles/export_reaction_smiles and explicit-IO `_with` variants taking a supplied
  Resolver and returning text by convey then render. Compose projection, construction, rendering,
  and reaction-side errors with source chains. Verify equality with the explicit composition,
  same-config guarantees, deterministic repeated output, and failure without mutation. Do not add
  redundant byte-output or sink APIs without a concrete approved contract.

**Gate:** a useful molecular and reaction text roundtrip works through the public output API,
including ordinary atoms, supported aromatic systems, tetrahedral sites, and coupled cis/trans
cases. Assert successful coverage in each category; document actual limitations precisely.

### S8 — Interpret/ingest resolver migration and existing consumers

- **S8a — Supplied-resolver interpretation and ingestion.** Modules: ingest.rs and all Rust
  callers/tests/examples/benchmarks. **Breaking (red→green).** [dep: S7c]
  Change Interpret to accept &Resolver and have molecule/reaction interpretation reuse it after
  raise. Change ingest_smiles/ingest_reaction_smiles and byte/default/`_with` families consistently:
  chemistry comes from the supplied resolver, default variants select only OpenSMILES IO config.
  Preserve existing raise, contradiction, underdetermination, and execution diagnostics. Migrate
  every caller in this subitem; do not silently select a different chemistry preset while replacing
  model/config arguments. Compare old fixture outcomes and test export/ingest composition directly.
- **S8b — Existing Python adapters.** Modules: umol-py molecule/reaction ingestion and associated
  tests. **Breaking caller migration (red→green with S8a).** [dep: S8a]
  Construct and pass the appropriate Rust Resolver from the existing Python model/config inputs,
  preserving their semantic defaults and lifetimes. Apply python-build, use Python 3.13, rebuild
  the native extension, and run affected tests. This required migration adds no new Python boundary
  or resolver wrapper. New Python export surface is deferred pending its explicit naming/ownership
  decision; do not pretend existing Python ingestion exposes the new output operations.

**Gate:** the whole workspace builds and tests with no stale Interpret/ingest signatures. Python
uses a freshly rebuilt extension. S8a and S8b form one green stage; no intermediate broken Python
consumer may be treated as completion.

### S9 — Integration evidence, specification, and closeout

- **S9a — Full property/conformance coverage.** Modules: IO and graph property/conformance targets.
  **Additive (green).** [dep: S2a, S3e, S7c, S8b]
  Run project(resolve(input)) for supported raised foreign inputs, boundary normalization,
  same-config construction/rendering, and ingestion-to-export preservation under the agreed
  equivalences. No arbitrary-GraphIR resolve(project(molecule)) law is required. Include both
  valence strategies, supported model policies, source-order variations, explicit/implicit H
  distinctions, aromatic grouping, marker-selection impossibility and alternatives, reactions
  and correspondence. Add retained
  regression cases for discovered defects. A success-conditional property must be paired with
  required-success fixtures; no shrinking the feature domain merely to pass tests.
- **S9b — Specification and public documentation.** Modules: umol-io/spec/opensmiles-spec.md,
  affected public rustdoc and usage examples. **Documentation (green).** [dep: S9a]
  Apply only the separate Staged specification updates list below to the spec: interpretation and
  equivalence, not pipeline/roundtrip promises. State constructor, traversal, projection, and
  rendering laws in their owning public APIs without citing this discussion. Examples show the
  parse/interpret and convey/render pairs and ingest/export conveniences with supplied resolvers.
  Describe the supported domain and explicit failure cases, including Either notation limitations.
- **S9c — Final gates and evidence reconciliation.** Modules: workspace and this record/index.
  **Verification/documentation (green).** [dep: S9a, S9b]
  Run the final commands below, inspect the entire change against the settled symbol inventory,
  and report bounded parser/raise/project/render benchmark evidence. Reconcile the independent
  doc 224 fuzz campaign with actual execution evidence; it must not be silently claimed complete.
  Record implemented scope, limitations, and deferred decisions here. Update status only when
  required work is complete; do not close outstanding scope by relabeling it a future enhancement.

### Final verification commands

Activate umol-py/.venv and confirm Python 3.13 before the workspace/Python gates, following
python-build. Commands are intended gates for implementation, not commands run during planning:

```text
cargo fmt --all
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo test -p umol-graph-core --features proptest --test property
cargo test -p umol-graph-ir --features proptest --test property
cargo test -p umol-io --features proptest --test smiles_property
cargo test -p umol-graph --features proptest --test property
cargo test -p umol-io --features conformance
cargo test -p umol-graph --features conformance --test resolution
maturin develop --manifest-path umol-py/Cargo.toml
pytest -q umol-py/tests
git diff --check
```

Run explicitly any additional feature-gated target introduced by a stage; default workspace tests
are not proof that those suites ran. Fuzz execution is a separate recorded campaign, not a cargo
test substitute. Benchmarks use the inspected existing targets plus focused new cases/targets,
with commands and results recorded when actually run.

### Critical path and deferred work

The writer path is S0 → S1 → S3 → S5 → S6 → S7. Projection S4 can proceed independently of
TableIR/writer implementation after S0 and joins at S7. S2 is the required traversal-consumer
migration before final integration. S8 restores the complete ingest/export symmetry and existing
Python consumers; S9 supplies the acceptance evidence. No stage may leave public invariants
unfinished behind temporary visibility or compatibility layers.

Deferred and not prerequisites for the Rust roundtrip: minimum-information projection, early
parity pruning, minimum-marker optimization, canonical output, cached traversal/markers,
new Python boundary/output surface, CTfile writing,
extended stereo, coordinate generation, and depiction wedge policy. Expand writable coverage
through concrete supported encodings and tests. A general counts inverse or a proof of arbitrary
GraphIR reconstruction is not a prerequisite or planned deliverable. No separate optimization
implementation stages are approved by this plan.

## S0a inventory and verification — 2026-09-10

Completed against source commit `0d8122af465a2347ce0ccdac984f36eb3271e905`, with a clean worktree
at the start of this subitem. This section inventories implementation work; it adds no production
API, tests, or algorithm. Tests below establish current behavior only. Later subitems remain pending.

### Public surface inventory

The signatures already specified above are authoritative. The following inventory identifies
changed and retained symbols, their consumers, and unresolved mechanical details rather than
inventing extra public seams.

| Owner and symbols | Planned change and contract | Consumers/migration |
| --- | --- | --- |
| graph-core traversal: DepthFirstEvent, BreadthFirstEvent | New public event enums with Discover/Finish/FinishTree; DFS also NonTreeEdge. Values are traversal output, not graph handles or persistent correspondence witnesses. | New visitors and event collectors; IO writer consumes node/edge identities. |
| traversal::visit_depth_first, traversal::visit_breadth_first | New public callback kernels with the bound/root/iterator/ControlFlow signatures above. No validation Result; correctness assumes consistent connectivity, internal handling remains panic-free. | Graph methods supply CSR access; IO supplies temporary table incidence. |
| Graph::visit_depth_first, Graph::visit_breadth_first | New inherent methods omit the free functions' bounds/neighbor arguments; borrow Graph and retain explicit roots, BFS depth limit, and visitor. | Event collectors, component visitor, neighborhood. |
| Graph::enumerate_depth_first_events, enumerate_breadth_first_events | New inherent collectors return the corresponding Vec<Event>, retaining root order and BFS depth limit. | Tests, callers needing retained traversal descriptions. No new persistent traversal type. |
| Graph::visit_connected_components | New ControlFlow visitor over one sorted borrowed node slice, with ConnectedComponentsAlgorithm. | Existing enumerate_connected_components becomes its collector; FinishTree owns completion/early-exit boundary. |
| Graph::enumerate_connected_components, ConnectedComponentsAlgorithm::Bfs | Retain names, selector, sorted members and component order; replace mislabeled LIFO flood fill with shared FIFO BFS. | Cycle-basis implementation; graph-IR GraphView wrapper; graph constraint validators and aromatic HMO; algorithms benchmark. |
| Graph::neighborhood; TraversalAlgorithm | Keep operation name, rename selector to NeighborhoodAlgorithm; depth limit and returned distances become usize. | graph-core refinement.rs, traversal tests, lib.rs reexport. No current GraphView neighborhood wrapper was found. |
| CircularRefinementAlgorithm::Ec { radius }, Graph::circular_refine | Radius feeds neighborhood through circular_refine_ec/remove_duplicate_environments; enumerate the transitive usize migration before S2 edits. | CircularRefinementHash::combine round argument; graph Morgan/ECFP featurizers, hash schemes, Python fingerprint configs; details below. |
| table_ir::StereoBond, BondConfiguration, BondRelation | New open records/enums with the settled site-only Either or two-reference Framed layout. No arbitrary complete GraphIR ligand frame copied into this type. | Parser finalization, CTfile/CX conversion, raise, checked boundaries and writer. |
| table_ir::Molecule::stereo_bonds, ExtendedMolecule::stereo_bonds | New collections alongside stereo_atoms; preserve through empty constructors and basic/extended conversion. | All molecule struct literals, reaction component carriers, SMILES targets, CTfile construction, tests. |
| Bond/ExtendedBond::direction and ::stereo | Retire persistent fields and update constructors, updates/conversions and literals. Retain wedge/donation/ring/source fields for their own meanings. | Parser builder and utils, CX application, CTfile conversion, TableIR raise and tests. |
| BondDirection, BondStereo and lexical conversion helpers | Inventory lexical use separately from table-field retirement. BondDirection is still needed to interpret source tokens; BondStereo currently represents CT/CX codes. Do not delete lexical vocabulary just because its stored field disappears. | Basic/extended parser targets, ring conflict handling, CT conversion, CX entries. Exact residual public exposure is reconciled in S3. |
| StereoAtom, StereoLigand, Winding, ConfigurationScope | Retain public layout and the implicit-H/LP participant convention. | Doc 224 parser finalizer, raise, future writer. No LP-removal migration. |
| AtomNeighbors and table_ir::Neighbor; both molecule atom_neighbors methods | Remain operation-local lookup; no new foundational graph field or adjacency trait. Eliminate unnecessary allocation at callers without assuming final Graph connectivity matches the table. | Raise helpers, table tests; future parser/writer incidence access. Retirement/replacement of the public helper itself is not approved. |
| TryIntoIr<Molecule> for &table_ir::Molecule; RaiseError | Raise explicit frames instead of interpreting retained directions; preserve model-independent representation errors at the owning conversion. | Interpret, parser/raise fixtures, fuzz target. Source marker conflicts move with their producer; do not erase diagnostics. |
| Resolver::project and per-resolver project functions | Direct GraphIR-only transformation, same mutation/publication semantics as resolve; no TableIR or IO-config parameters. Ordinary valence has no projection phase; S4a removed ValenceResolver::project and ValenceProjectError. S4b reconciled aromatic projection's signature and five diagnostics in its completion record above. | Convey; supported-input roundtrip properties and projection benchmarks. |
| Resolver::resolve, ResolveConfig, ResolveState and existing reports/errors | Retain forward semantics, including completed isotope corrections. Project does not invoke resolution for validation. | Existing ingestion and resolution tests; supported-input projection properties. |
| AtomTypeRegistry lookup/admission operations | Retain forward admission and completed registry invariants/deduplication; no projection candidate filtering or new pruning API. | Forward atom typing. |
| Smiles/ReactionSmiles::from_table_ir | New checked owned-table constructors with &SmilesIoConfig; Result establishes renderability. Private table fields remain private. | Convey across the crate boundary; tests of independent open tables. |
| Smiles/ReactionSmiles::render, render_with | New String-returning Result methods; OpenSMILES default or explicit IO config. Recompute traversal/assignment. | export conveniences; boundary normalization properties. No Display or unchecked formatter is approved. |
| Smiles/ReactionSmiles::parse, parse_bytes, parse_with, parse_bytes_with, FromStr | Retain parse surface and syntax-config defaults; frame normalization changes their produced internal tables. | Existing ingestion, direct boundary users, properties and fuzzing. Parsing does not become chemistry resolution. |
| Smiles/ReactionSmiles::as_table_ir, into_table_ir | Retain signatures; borrow/consume private payload, no mutable accessor. | Interpret is the actual production consumer. |
| Interpret::{Output, Error, interpret} | Change method to interpret(&self, &Resolver); keep associated output/error types. | Both boundary implementations in ingest.rs and every ingestion caller. |
| Convey::{Input, Config, Error, convey} | New local trait in umol-graph, implemented on output boundary types with the settled signature; Config = SmilesIoConfig here. | Molecular/reaction export. No Export trait or generic format target on Molecule. |
| ingest_smiles and ingest_reaction_smiles families | All eight existing text/bytes/default/explicit-IO functions receive a supplied Resolver; remove separate model/resolve-config arguments. | Full caller list below, including benches and Python adapters. Default variants retain only default IO selection. |
| export_smiles/export_reaction_smiles and their _with variants | New four functions return text through convey + render, taking a Resolver and default/explicit IO config. | Rust output callers. No additional sink/byte family or Python wrapper is inferred. |
| ProjectError, construction/convey/render/output diagnostics | Operation boundaries and proposed names are settled; exact variants, report types, reaction names and shared-reason placement remain for their cohesive implementation subitems. | Error source chains, section context and exact failure fixtures. No catch-all silent fallback. |

Argument ordering not explicitly shown in the approved signatures, error enum layout, and
per-resolver helper signatures must be reconciled at their owning subitems; the inventory does
not label invented spellings as approved APIs. New symbols should be exported only through the
owning module and the repository's established public reexport surface, not from every crate root.

### Producer and consumer paths

- **Current doc 224 implementation:** smiles/parser/builder.rs has PendingStereo, opening-order
  bond slots completed by complete_bond, and finish that builds explicit StereoAtom records.
  Bracket H count and the three-incidence/H0 LP convention are applied there; the ring completion
  map remains separate for CX indexing. Both Basic and Extended implementations of Target use it.
  The completed document contains historical descriptions of the old bug; those are not current
  implementation requirements.
- **SMILES and CX:** smiles/parser.rs parses basic/extended molecules and reaction sections,
  then remaps CX completion-order bond indices through BondIndexMap before applying annotations.
  smiles/parser/cx.rs has separate basic and extended branches setting Cis, Trans and Either.
  The frame producer must run with final interpreted annotations and detect conflicting evidence;
  early finalization cannot silently bypass later CX assertions.
- **CTfile:** ctfile/parser/bond.rs and parser/convert.rs interpret the bond codes and wedge
  orientation, while parser.rs constructs basic/extended tables with supplied positions. Current
  definite double-bond geometry is interpreted later in table_ir/raise.rs and raise/utils.rs.
  S3 moves the appropriate interpretation to frame production while preserving wedge responsibility.
  V2000 readers, accumulator overrides, and basic/extended conversion must all reach the same
  result. A future V3000 reader must follow the same contract; no V3000 parser exists currently.
- **Table conversions:** table_ir/molecule.rs implements From<Molecule> for ExtendedMolecule and
  TryFrom<ExtendedMolecule> for Molecule. ExtendedBond has both fields being retired, not merely
  Bond. Reaction and ExtendedReaction contain three corresponding molecule carriers and derived
  mapping indexes. Their from_molecules/from_extended_molecules constructors currently initialize
  that index empty; checked boundary construction must not assume an independently assembled
  table carries a valid index.
- **Raise:** table_ir/raise.rs unconditionally builds AtomNeighbors, constructs atoms/bonds,
  partitions donating/noncovalent records into relations, copies explicit atom frames, and invokes
  GraphIR Molecule::try_from_entries. That constructor validates references, builds Graph from
  localized bonds, and checks integrity. Its graph-accepting try_from_arcs is private. No new
  public arbitrary-parts GraphIR constructor is needed for this task.
- **Resolvers:** resolve.rs coordinates valence admission/aromatic selection, atom/system edits,
  stereo planning, localized/multicenter defaults, discharge, and a final concrete-only publication.
  resolve/{valence,aromaticity,stereo,bonds,multicenter}.rs own those domains. ops/stereo.rs supplies
  the existing #T/#C frame derivation and coset transport. Projection must preserve these layer
  boundaries and must not call IO from a resolver to test expressibility.
- **Reaction conversion:** ingest.rs reads Atom.class-derived atom_mapping, rejects agents and
  repeated classes for GraphIR interpretation, raises/resolves both sides and constructs their
  correspondence. Output must use existing side materialization and correspondence, not infer a
  new map. Correspondence::matched_pairs is sorted by left ID; pair i supplies label i + 1.
- **Ingestion callers:** ingest.rs tests; graph benches resolve.rs, fingerprint.rs, substructure.rs;
  fingerprint/{morgan,reaction,pattern}.rs tests; graph tests fingerprint.rs, whitepaper.rs and
  property/publication.rs; Python molecule.rs, reaction.rs and error.rs, including embedded tests.
  Workspace search found no separate graph/examples directory; executable examples live in the
  declared bin/test targets and must be searched again when signatures change.
- **Python:** from_smiles/from_reaction_smiles currently accept optional IO/model/resolve configs,
  construct the SMILES valence preset by default, and call Rust ingestion. S8 preserves that
  behavior by constructing a Rust Resolver with the required lifetime. There is no existing
  Python Resolver or Smiles boundary wrapper to assume. New Python output exposure is deferred.
- **Evidence targets:** graph-core benches/algorithms.rs and feature-gated tests/property.rs;
  IO benches/smiles_parsing.rs, tests/smiles_property.rs, parser unit tests and conformance targets;
  graph benches/resolve.rs, tests/property.rs and conformance resolution target; graph-IR property
  tests and GraphView component fixtures. The 60-minute parser campaign is independent evidence.

### Construction, mutation, and failure ownership

| Value/operation | Integrity and contextual contract | Failure and preservation rule |
| --- | --- | --- |
| Graph and traversal callback | Graph establishes its stored connectivity. An arbitrary callback carries the consistency precondition; traversal does not certify it. | No defensive validation Result; inconsistent finite values have no correctness promise, internal processing must avoid panics. Callback execution is caller-owned. |
| Open TableIR and temporary index | Tables remain authoritative. Index construction establishes its own storage; later table edits promise no index/table correspondence. | No stale-index repair or general validation wrapper. Operations needing table references first establish their required properties; no chemistry inference in parse. |
| Explicit stereo frames | Site and references belong to the table; endpoint ordering and frame action are fixed. Either has no references. | Source interpretation checks marker consistency; raise/checked construction enforce required references. Do not silently pad, retag virtual ligands, or downgrade configuration. |
| Checked SMILES boundary | Private payload plus successful construction under IO config establishes representability; parse retains its existing acceptance contract. | Mutable fields are not exposed. Same-config rendering succeeds; other-config rendering may fail. No invented coordinates or lost annotations. |
| Resolver::project | Caller supplies GraphIR; work is GraphIR throughout, on an intermediate candidate. | Result/Solution follows resolve; successful projection publishes atomically, without re-resolution. Failure preserves caller input. |
| Convey and text export | Convey clones, projects, lowers values supported by TableIR, and calls checked construction for notation support. Export composes convey/render. | Preserve source graph and supported input semantics; reject concrete unrepresentable values. No runtime roundtrip verification. |
| Reaction boundary construction | A supplied table may have independently assembled labels/index; the constructor owns the required consistency check. | Preserve label semantics including parsed repeats/agents where representable; GraphIR convey produces its index from correspondence labels. |

The laws to verify in later subitems are traversal visitor/collector equivalence and completion,
project(resolve(input)) for supported foreign input plus atomic failure, boundary normalization
and same-config renderability, and molecular/reaction ingestion-to-export preservation with
explicit/implicit H and stereo-frame action.
These are distinct operational domains; no structural count comparison substitutes for them.

### Migration findings to retain in later stages

1. **Size migration reaches fingerprints.** neighborhood is called by
   refinement.rs::remove_duplicate_environments with the circular radius. CircularRefinementAlgorithm
   stores radius: u32, CircularRefinementHash::combine takes round: u32, graph MorganFeaturizer and
   EcfpFeaturizer expose u32 radii, and Python HashedFingerprintConfig has Morgan/Ecfp radii.
   Include these callers in S2's size review, including graph/src/hash.rs implementations and
   fingerprint golden fixtures. usize sizes do not by themselves authorize changing a named hash
   recipe's serialization or output bits; reconcile that detail before changing its combine input.
2. **Stored versus lexical directions.** BondDirection remains source-token data in parser
   state. ParseError::MismatchedRingBondDirections stores source positions, not a BondDirection
   value. Keep the mismatch diagnostic when removing Bond.direction. Likewise, CT/CX decoding
   still needs to distinguish source codes after Bond.stereo is removed.
3. **Panic-reporting correction after S0a.** The initial inventory incorrectly concluded that
   the target's discarded catch_unwind result hides panics under cargo-fuzz. libfuzzer-sys 0.4.13
   installs an aborting panic hook before target execution; the hook runs before unwinding.
   A controlled caught-panic reproducer confirms libFuzzer reports a crash. The redundant inner
   catch is now removed. See the follow-up below; the user's completed run is not invalidated by
   this wrapper.
4. **Checks are baseline evidence only.** Current IO tests exercise persistent directions and the
   old raise path. Their passing result neither validates future frame normalization nor proves
   the inverse projection/writer. S0b supplies independent expectations; later stages implement
   and test the new laws without rewriting them to match behavior.

### S0a verification results

All commands ran with --offline on the pinned source; no production or test code changed.

| Command | Result |
| --- | --- |
| cargo test -p umol-graph-core --lib algorithms::traversal --offline | 5 passed |
| cargo test -p umol-graph-core --lib algorithms::connectivity::components --offline -- --quiet | 5 passed |
| cargo test -p umol-io --lib --offline | 3,589 passed; includes parser, CTfile and TableIR/raise unit coverage |
| cargo test -p umol-graph --lib ingest:: --offline -- --quiet | 161 passed |
| cargo test -p umol-graph --lib ops::resolve:: --offline -- --quiet | 141 passed |

Total: 3,901 passing tests, no failures. cargo fmt --all -- --check also passed, with stable-rustfmt
warnings about the repository's nightly-only formatting options; git diff --check passed.
Feature-gated property/conformance suites, Python builds,
benchmarks, and fuzz execution were not run in this inventory subitem. Those remain in their
own planned work. S0b was the next subitem at this inventory's closeout.

## Parser fuzzing follow-up — 2026-09-10

This work was explicitly authorized separately from S0b: remove the target's panic catch, audit
what it exercises, and improve stereo/ring seed representation. It preceded the S0b work below.

The target executes Smiles::parse_bytes on arbitrary bytes with the OpenSMILES configuration,
then TryIntoIr<Molecule> on successful parses. Returned parse/raise errors are ordinary outcomes.
It performs no resolution, chemical conformance test, or stereo-equivalence comparison. Its
useful evidence is crash/timeout/sanitizer robustness over the executed inputs, not stereochemical
correctness. Panic-free but incorrect stereo requires the independent unit/property expectations.

The earlier claim that the 60-minute run lost panic detection was wrong. The fuzz workspace locks
libfuzzer-sys 0.4.13; its initialize function installs a hook that aborts before unwinding, before
an inner catch_unwind can return. The same behavior is present in locally available 0.4.12.
A scratch cargo-fuzz target containing a deliberate panic inside a discarded catch_unwind
reported `SUMMARY: libFuzzer: deadly signal` and target exit status 77. The source wrapper was
still removed so the target clearly lets failures propagate without redundant recovery code.

The user reports the previous run completed for 60 minutes. That is meaningful robustness work
under the normal cargo-fuzz runtime; no exact historical execution/coverage count is claimed here.
The live corpus currently contains 13,627 inputs totaling 5,991,687 bytes. A byte-level scan found
2,113 with @, 3,374 with slash/backslash, and 2,029 with @ plus a digit. These are lexical counts,
not evidence that those inputs parse, reach a stereo frame, or came exclusively from that run.

### Seed audit and changes

The checked-in seed set previously had 13 inputs. Its stereo-marked entries exercised rejection
paths; none successfully raised with a tetrahedral frame or a cis/trans assertion. The old
ring_dir seed itself is a mismatched-direction parse rejection. Retain these regression seeds,
but do not mistake them for successful stereo-path coverage.

Add 38 named seeds, for 51 total. They cover four actual tetrahedral ligands, bracket and explicit
H, LP assertions, opening/closing ring digits, percent labels, multiple closures, fused/bridged/
spiro rings, adjacent centers, later components, ring-label reuse and closure across a dot;
also definite/partial/redundant cis/trans, four substituents, shared chains, branching and cycles,
and ring markers supplied at opening/closing/both ends. Add explicit malformed stereo/ring cases
and aromatic fused systems, bracket H, and localized links. No existing seed was removed.

An actual parse/raise census of the final seeds found 41 successful raises, 3 parse rejections,
and 7 raise rejections. Twelve successful inputs contain tetrahedral frames, and eleven produce
cis/trans assertions, including branched and cyclic shared-marker cases. This verifies that
seeds reach the intended broad paths; it does not validate each configuration's chemical meaning.
The scratch audit source/results are in scratch/stereo-valence-scan/src/bin/seed_audit.rs and
scratch/fuzz-seed-audit.jsonl. The tracked fuzz/README.md documents the target, seed categories,
and commands that explicitly load both seeds and the dictionary. Neither is automatically
loaded by the default cargo-fuzz command. New corpus inputs should go to the first, mutable corpus
directory, preserving the checked-in seed directory.

### Verification

- Controlled caught-panic cargo-fuzz reproducer: crash detected, target exit status 77, as expected;
  source/log under scratch/fuzz-panic-proof. No deliberate panic was inserted in production code.
- Updated real target built with nightly cargo-fuzz and default AddressSanitizer/debug assertions.
  Replay of the initial expanded 47 seeds with the existing 98-entry dictionary succeeded.
- Final 51-seed set was copied to scratch/fuzz-smiles-check-corpus and used for a bounded 30-second
  mutation run with the dictionary: 547,910 executions in 31 seconds, no crash reported. libFuzzer
  counters increased from cov 2,764 / ft 5,293 to cov 4,761 / ft 17,113. These are instrumentation
  counters, not percentages of code or semantic coverage. The log is scratch/fuzz-smiles-check.log.
- The live accumulated corpus was not modified, and no second 60-minute campaign was started.
  Broader campaign evidence remains separate from S0b's semantic fixtures.

## S0b semantic fixtures and comparisons — 2026-09-10

This subitem adds tests only. Parser, raise, resolver, and boundary APIs are unchanged.
Existing independent DSL expectations already cover ordinary molecules, localized charges,
radicals, aromatic system membership/contributions, and ring-opening/closing stereo frames.
Retain those expectations and the explicit-remapping frame tests; do not replace them with
parse/render agreement when output becomes available.

### Comparison contract

- **Boundary law:** compare interpreted boundary information under an atom/bond correspondence,
  ignoring source spans and equivalent spelling (ring numbers, traversal, branch order, and the
  agreed redundant direction markers). Preserve each atom/bond property, explicit atoms versus
  bracket-H counts, atom classes/labels, and supported extension data. Compare stereo in the
  corresponding ligand frame, including the winding/relation action. Definite, Either, and absent
  assertions are distinct; a consistent partial directional marking contributes no assertion.
  Check marker consistency before discarding redundant notation. For reactions, also preserve
  side/agent membership and mapping classes. GraphIR equality alone is insufficient for this law:
  raise does not retain all boundary data. A general boundary comparator follows the explicit
  stereo-bond representation in S3; S0b records concrete expected semantic observations without
  introducing a helper that erases fields from today's mixed source representation.
- **GraphIR comparison tools:** the corrected required law is project(resolve(input)) ≡ input
  over supported raised foreign inputs, under the agreed representation equivalences. It is not
  raw equality of the raised and projected molecules. For comparisons of complete molecular states
  in independent ingestion fixtures or optional reverse-direction experiments, use
  Molecule::framed_eq when entity IDs are retained,
  framed_eq_under when an explicit remapping is supplied, and aggregate canonical_eq when text
  traversal changes IDs. Compare the complete molecule, including atom states, localized bond
  states, aromatic system members and electron contributions, and stereo entities/configurations.
  No hydrogen folding, charge relocation, or overlay removal is part of comparison. Independent
  expected molecules in ingestion tests use stored equality where their exact frame/order is
  specified. Existing mirror/alkene controls and the new explicit-versus-implicit-H control
  prevent equivalence checks from erasing those distinctions.

### Executable evidence

- umol-graph/src/ingest.rs extends independent complete-molecule DSL expectations with four-atom
  tetrahedral carbon, bracket H, explicit H, a carbanion lone pair, both sulfoxide bond/charge
  representations, and both alkene configurations. Additional counts/atom-typing cases check
  exact central atom states and ligand identities for C/N/P/S, including N/P cations with four
  actual ligands. Neutral trigonal C and C+ reject at stereo resolution; too few/many ligands
  and duplicate bracket H reject at raise with exact errors. Atom typing rejects the recorded
  overvalent charged C/S cases. Counts-model acceptance of the known problematic states is not
  promoted into an expected-success contract; that independent issue remains with doc 166.
- umol-io/src/table_ir/raise.rs retains ring-frame and explicit-remapping regressions and adds
  independent coset/absence expectations for redundant four-substituent notation, endpoint
  viewpoint, partial conjugated systems, trienes, explicit-H marker choices, and two-ended ring
  markers. A coordinate-free CX case distinguishes Either from absence and definite configuration;
  its atom class is checked before raise.
- umol-io/tests/smiles_property.rs generates chains of one through eight double bonds with each
  adjacent single bond unmarked, forward, or backward. Expected assertions are derived directly
  from the chosen neighboring glyphs, independently of parser/raise results. A second property
  varies ring numbers, atom classes, and both tetrahedral windings against a literal encounter
  frame. These test input semantics; projection and render laws become executable with their
  implementing stages.

### S0b verification

- cargo test -p umol-io --lib --features proptest --test smiles_property --offline -- --quiet:
  3,602 unit tests and eight property tests passed.
- cargo test -p umol-graph --lib --test property --features proptest --offline -- --quiet:
  1,061 unit tests and five property tests passed.
- cargo clippy -p umol-io -p umol-graph --lib --tests --features proptest --offline -- -D warnings:
  passed; log under scratch/s0b-clippy.log.
- cargo fmt --all and git diff --check passed. Test logs are under scratch/s0b-io-tests.log and
  scratch/s0b-graph-tests.log. No production APIs, algorithms, resolver acceptance rules, or
  existing assertions were changed. Projection, rendering, and their future property gates are
  still unimplemented; this subitem establishes their independent input-side expectations.

S0c followed with the bounded baselines below. The separately completed fuzz follow-up supplies
its fuzz status; do not duplicate the user's campaign.

## S0c inline benchmark examples and baselines — 2026-09-10

All new benchmark examples are inline in the existing benchmark sources. The 16 molecular cases
cover a 64-atom chain, branching, components, localized charge/radical states, fused aromaticity,
aromatic lone-pair donation, tetrahedral actual/bracket-H/explicit-H/LP frames, and definite,
partial, branched, and cyclic directional systems. The six graph examples cover paths, a binary
tree, a cycle, disconnected cycles, and a graph with a loop, parallel edges, and an isolated node.
These are bounded synthetic workloads, not a representative corpus or scale claim.

### Timing boundaries and reuse

- umol-io/benches/smiles_parsing.rs adds smiles_roundtrip/parse and smiles_roundtrip/raise.
  Raise borrows a pre-parsed TableIR; parsing is not inside its timing. Both use Criterion::iter,
  including destruction of their returned boundary/GraphIR values.
- umol-graph/benches/resolve.rs adds smiles_roundtrip/resolve using ValenceModel::smiles()
  (counts, MostSaturated) with the other chemistry/resolve defaults. Each input must resolve to
  Determined during setup. Criterion::iter_batched_ref clones the raised input outside timing;
  parsing, raising, resolver construction, and destruction of the mutated input are also outside.
  This isolates resolution from the existing end-to-end ingest benchmark.
- umol-graph-core/benches/algorithms.rs adds traversal_baseline/components and
  traversal_baseline/neighborhood over prebuilt graphs. Components uses the current Bfs selector
  whose implementation is still a LIFO flood fill. Neighborhood uses Bfs from node 0 without a
  practical depth bound. On disconnected graphs it visits only node 0's component; components
  enumerates all components. These timings are not measurements of identical operations.
- S1 adds visitor/event measurements on these graph shapes. S3 reuses the isolated raise cases.
  S4 prepares resolved inputs from the molecular examples outside projection timing; S5/S6 add
  marker-assignment and rendering measurements on the corresponding normalized tables. Those
  operations do not exist yet and have no measurements or placeholder implementations here.

### Bounded timing run

Source: e64f787f99d4704bda27e27bf8fa1d0ab8f214c0 plus the S0c benchmark additions;
rustc 1.96.0, macOS 15.7.3, arm64, normal workspace bench profile (optimized with debug info).
The three benchmark runs were sequential, with compilation outside measurement and the normal
allocator. Each case used 20 samples, 0.1-second warm-up, 0.2-second measurement, and 1,000
resamples. No tuning or repeat campaign was performed.

Run these commands from the workspace root, appending
`--sample-size 20 --warm-up-time 0.1 --measurement-time 0.2 --nresamples 1000 --noplot --save-baseline s0c`
to each:

```sh
cargo bench -p umol-io --bench smiles_parsing --offline -- smiles_roundtrip
cargo bench -p umol-graph --bench resolve --offline -- smiles_roundtrip
cargo bench -p umol-graph-core --bench algorithms --offline -- traversal_baseline
```

Criterion time point estimates in microseconds per operation:

| Example | Parse | Raise | Resolve |
| --- | ---: | ---: | ---: |
| chain_64 | 0.983 | 12.232 | 104.240 |
| branched | 0.274 | 1.882 | 11.551 |
| components | 0.194 | 1.012 | 5.965 |
| aromatic_fused | 0.323 | 2.869 | 37.258 |
| aromatic_lone_pair | 0.244 | 1.659 | 18.483 |
| charged | 0.151 | 0.576 | 1.766 |
| radical | 0.148 | 0.571 | 1.778 |
| tetra_four | 0.297 | 2.225 | 10.345 |
| tetra_ring | 0.325 | 2.577 | 20.556 |
| tetra_explicit_h | 0.371 | 2.907 | 15.617 |
| tetra_lone_pair | 0.285 | 2.231 | 10.756 |
| alkene_four | 0.244 | 2.015 | 11.873 |
| shared_chain | 0.216 | 2.236 | 14.390 |
| partial_triene | 0.246 | 2.576 | 14.993 |
| shared_branch | 0.272 | 2.553 | 16.096 |
| shared_cycle | 0.276 | 2.760 | 33.894 |

| Graph | Components | Neighborhood |
| --- | ---: | ---: |
| path_64 | 0.537 | 0.556 |
| path_1024 | 6.565 | 6.119 |
| binary_tree_255 | 2.752 | 2.060 |
| cycle_64 | 0.744 | 0.556 |
| four_hexagons | 0.508 | 0.158 |
| loops_parallel_isolated | 0.114 | 0.086 |

These short runs establish initial baselines, not fine performance rankings. For example,
tetra_ring resolve has a wide reported interval, 16.122–26.991 microseconds. Preserve the intervals
in scratch/s0c-timings.csv and original scratch/s0c-{io,graph,core}-bench.log files when comparing
later changes. Criterion also saved the s0c baseline under /Users/dr/.cargo-target/criterion.
Do not sum the columns as an end-to-end time: their setup/destruction boundaries differ.

### Separate allocation snapshot

The scratch/stereo-valence-scan/src/bin/s0c_allocations.rs probe uses the same inline molecular
examples and a counting System allocator in a release build. It warms parse/raise, then measures
one parse, neighbor lookup, and raise per example separately. Input creation is outside each
measurement. The allocator is not installed in the timing benchmarks or production crates.

Counts below include successful alloc, alloc_zeroed, and realloc calls. Requested bytes sum full
requested sizes, including the full new size on realloc; they are not net retained memory.
Peak bytes mean added live requested Rust allocation sizes above the pre-operation baseline,
with the result retained. They exclude allocator overhead, native allocations, and any internal
old/new-buffer overlap during realloc; they are not process RSS.

| Example | Parse calls | Neighbor calls | Raise calls (including realloc) | Raise requested bytes | Raise peak added bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| chain_64 | 2 | 65 | 102 (9) | 51,828 | 27,320 |
| aromatic_fused | 3 | 11 | 52 (5) | 13,636 | 8,520 |
| tetra_ring | 5 | 8 | 53 (5) | 7,284 | 4,688 |
| shared_cycle | 3 | 9 | 55 (3) | 7,140 | 4,824 |

The chain's neighbor lookup alone requests 3,584 bytes in 65 allocations; full raise includes
that lookup. This is a concrete S3 allocation baseline, not an estimate of time attributable to
the lookup. The complete 48-row snapshot is scratch/s0c-allocations.csv. Reproduce with:

```sh
cargo run --release --offline --manifest-path scratch/stereo-valence-scan/Cargo.toml \
  --target-dir /Users/dr/.cargo-target --bin s0c_allocations
```

### Fuzz status and verification

The user-reported completed 60-minute parser/raise campaign and the subsequent target/seed audit
are recorded above. The original hour-long run log was not available in the inspected scratch
outputs; no historical execution or coverage count is inferred. The corrected target, 51 named
seeds, and documented command are already in the repository. S0c did not start another campaign.

All 60 timing benchmarks and the 48 allocation observations completed successfully. The two
molecular benchmark groups preflight parse/raise success, and all 16 resolver inputs preflight
Determined outcomes. These checks establish successful workloads, not independent chemistry
correctness; S0b's exact expectations and properties supply that evidence.

Validation passed: 3,602 IO unit tests and eight SMILES properties; 1,061 graph unit tests and five
graph properties; five graph-core traversal and five connected-component tests. Strict Clippy
passed for all three changed benchmark targets; cargo fmt --all and git diff --check passed.
Logs are scratch/s0c-tests.log, scratch/s0c-core-{traversal,components}-tests.log, and
scratch/s0c-clippy.log. Only benchmark code and this record/index changed; S1a is next.

## S1a DFS contract and verification — 2026-09-10

The public additions are DepthFirstEvent, traversal::visit_depth_first (also exported at the
crate root), and Graph::visit_depth_first. DepthFirstEvent is an open descriptive enum: its
NodeId and EdgeId fields do not certify membership in a particular graph. There are no new
constructors, conversions, validators, stored graph views, or Python bindings.

For fixed, consistent undirected connectivity, the callback kernel preserves candidate-root
and neighbor order, emits nested discovery/finish events, reports each non-tree edge once by
identity, and emits FinishTree after each completed root. Discover's parent names the parent
node and connecting edge. Graph supplies its existing CSR. Bounds are usize; ids retain their
existing types. Suspended iterators live on an explicit stack.

The failure boundary is ControlFlow: Break stops immediately, including further root/neighbor
iteration and callback calls. Index access must remain panic-free for inconsistent finite
callback inputs, without a validation pass or a correctness guarantee for those inputs.
Callback and iterator panics or nontermination remain caller-owned. No separate error type or
integrity certificate is introduced. Tests cover ordered event traces, an independent recursive
reference, loops and parallel edges, sparse ids, early termination, malformed inputs, and a
long path. Benchmarks extend the existing inline traversal examples.

Implemented against source commit 71bf6a8f66c98902de62c25aa0e7b7fddb2a64be. Public-surface
reconciliation matches the three additions above, including enum fields, Fn neighbor callbacks,
usize bounds, and ControlFlow return values. No other public symbols changed. The kernel uses
node/edge visitation arrays and a stack retaining each active neighbor iterator; Graph supplies
borrowed CSR iterators directly.

Twelve new unit cases cover exact events and early termination, including a 100,000-node path.
Three new properties compare the kernel and Graph adapter against recursive DFS with set-based
visitation, check early-break prefixes, and exercise raw finite callback tables with invalid and
reused ids. Valid generated multigraphs have at most eight nodes and twenty edges; roots can be
subsets or duplicates, incidence order can reverse, and edge ids can remain sparse. These are
bounded generative checks, not an exhaustive proof. The immediate-break unit cases also reject
any subsequent root iteration, neighbor callback, neighbor iteration, or visitor call.

Validation passed: 979 graph-core unit tests, 54 integration tests, and 129 properties with the
proptest feature enabled. Strict Clippy passed for all graph-core targets with that feature;
rustdoc built without warnings. Commands and logs:

```text
cargo test -p umol-graph-core --features proptest --offline
    scratch/s1a-tests.log
cargo clippy -p umol-graph-core --all-targets --features proptest --offline -- -D warnings
    scratch/s1a-clippy.log
cargo doc -p umol-graph-core --no-deps --offline
    scratch/s1a-doc.log
```

All six new traversal/depth_first benchmarks use the existing inline graph examples and visit
every component with a black-box event consumer. Criterion's reported time estimates are:

| Inline example | DFS visitor |
| --- | ---: |
| path_64 | 959.25 ns |
| path_1024 | 12.635 µs |
| binary_tree_255 | 2.4593 µs |
| cycle_64 | 962.82 ns |
| four_hexagons | 296.18 ns |
| loops_parallel_isolated | 84.738 ns |

This establishes the new visitor's baseline; it is not a comparison with the different work
performed by neighborhood or component enumeration. The bounded run used 20 samples, 100 ms
warm-up and 200 ms measurement time per example, with 1,000 resamples:

```text
cargo bench -p umol-graph-core --bench algorithms --offline -- traversal/depth_first
    --sample-size 20 --warm-up-time 0.1 --measurement-time 0.2
    --nresamples 1000 --noplot --save-baseline s1a
```

The run log is scratch/s1a-core-bench.log. cargo fmt --all and git diff --check passed. Final
diff review found only the S1a API, tests, inline benchmark additions, and this record/index.
S1b is next; existing neighborhood and component implementations remain for their planned stages.

## S1b BFS contract and verification — 2026-09-10

The public additions are BreadthFirstEvent, traversal::visit_breadth_first (also exported at the
crate root), and Graph::visit_breadth_first. The event enum is an open descriptive carrier of
node ids, an optional parent node/edge pair, and usize depths. Its construction certifies no
relationship to a particular graph. There are no new constructors, conversions, validators,
transformations, or Python bindings.

For fixed, consistent undirected connectivity, the kernel preserves candidate-root and neighbor
order, discovers each node once when enqueued, finishes each node after its permitted expansion,
and emits FinishTree when that root's queue empties. Roots are sequential, not simultaneous
sources. Previously visited nodes remain visited across roots; shortest distances for a later
depth-limited tree are therefore measured in the remaining unvisited graph. At the depth limit,
nodes receive Discover and Finish without a neighbor callback. A limited tree need not be a
whole connected component. Graph supplies its existing CSR directly.

The kernel accepts usize node bounds and Option<usize> maximum depth, uses Fn neighbor callbacks,
and returns ControlFlow. It needs no edge visitation state or edge bound. Consistency is a caller
precondition for correctness, not a validation pass. Internal indexing remains panic-free for
inconsistent finite inputs; callback and iterator panics or nontermination remain caller-owned.
Break immediately stops all event, callback, and iterator execution. Ordered traces, independently
computed shortest distances, depth boundaries, duplicate roots, malformed inputs, and early
termination provide the tests. Benchmark inputs remain inline.

Implemented against source commit a47159f267ff2daeb25ab712ecb04b0b92ad2cfe. The public surface
matches the three additions above and their settled fields/signatures. A VecDeque and node
visitation array supply FIFO traversal; no common traversal engine or connectivity storage was
introduced. Existing DFS, neighborhood, and component implementations are unchanged.

Seventeen new unit cases cover exact event order, root order and duplicates, loops and parallel
edges, sparse edge ids, empty inputs, depth zero/exact/unrestricted/usize::MAX, and limited
successive trees. Depth-limit tests observe neighbor callback calls directly; Break tests reject
all later callback and iterator execution, including when other nodes remain queued. A limited
successive-tree case preserves previously visited nodes as a barrier to a later root's shortcut.

Three new properties cover the kernel, Graph adapter, and malformed finite callback tables.
Valid generated multigraphs have at most eight nodes and twenty edges. A layer-by-layer reference
checks exact events and early-break prefixes; independent repeated edge relaxation checks each
tree's reached set and shortest distances after excluding earlier trees. This is bounded
generative evidence, not an exhaustive proof. The malformed-input property requires normal
completion without asserting a useful event stream.

Validation passed: 996 graph-core unit tests, 54 integration tests, and 132 properties. Strict
Clippy passed for all graph-core targets with the proptest feature; rustdoc built without warnings.
Formatting and diff checks passed. Commands and logs:

```text
cargo test -p umol-graph-core --features proptest --offline
    scratch/s1b-tests.log
cargo clippy -p umol-graph-core --all-targets --features proptest --offline -- -D warnings
    scratch/s1b-clippy.log
cargo doc -p umol-graph-core --no-deps --offline
    scratch/s1b-doc.log
cargo fmt --all -- --check
    scratch/s1b-fmt-check.log
```

The six traversal/breadth_first benchmarks reuse the existing inline graphs, visit every
component without a depth limit, and consume each event through black_box. Criterion's reported
time estimates establish a new BFS baseline:

| Inline example | BFS visitor |
| --- | ---: |
| path_64 | 443.99 ns |
| path_1024 | 6.3370 µs |
| binary_tree_255 | 1.9734 µs |
| cycle_64 | 410.53 ns |
| four_hexagons | 201.12 ns |
| loops_parallel_isolated | 62.672 ns |

```text
cargo bench -p umol-graph-core --bench algorithms --offline -- traversal/breadth_first
    --sample-size 20 --warm-up-time 0.1 --measurement-time 0.2
    --nresamples 1000 --noplot --save-baseline s1b
```

All six measurements completed; the run log is scratch/s1b-core-bench.log. Final diff review
matches the S1b surface, tests, inline benchmarks, and this record/index. S1c is next: event
collectors and the component visitor/migration remain pending.

## S1c collectors and components — 2026-09-10

The additions are Graph::enumerate_depth_first_events, Graph::enumerate_breadth_first_events,
and Graph::visit_connected_components. The collectors retain their visitor's ordered roots
and, for BFS, Option<usize> depth limit, returning the exact complete event sequence as Vec.
These open collections add no graph-membership certificate, new constructors, validators,
transformations, or Python surface. Their contextual guarantees are those of the existing
Graph traversal visitors; invalid roots retain the visitors' panic-free handling.

The component visitor takes ConnectedComponentsAlgorithm and FnMut(&[NodeId]) -> ControlFlow<B>.
Graph supplies valid undirected connectivity. With Bfs, unrestricted traversal collects one
component, sorts its node ids, and visits the borrowed slice at FinishTree. Every node belongs
to exactly one emitted component; each component is a maximal connected node set. Members are
ascending and components are ordered by their least node id. Break returns immediately before
another component is explored. Callback panics/nontermination remain caller-owned. The slice
is borrowed only for the callback; callers retain a component by copying it.

Graph::enumerate_connected_components retains its signature, selector, ordering, and result
type, collecting through this visitor. Its old LIFO flood-fill implementation is removed.
No integrity pass, alternate traversal engine, or additional public helper is introduced.
Tests compare collectors with visitor sequences and independent references, and components
with definition-level reachability. Existing BFS FinishTree break tests establish immediate
termination of the underlying traversal; component tests check propagation and emitted prefixes.

Implemented against source commit 79694704162697e7096384475fd848ee092896da. The public-surface
audit matches the three additions and the retained enumerate_connected_components method above.
The event collectors contain only collection through their existing visitors. Component visitation
uses the existing BFS with one reusable component buffer and sorts it at FinishTree. The old
connected_components_bfs helper is removed. There are no new reexports, boundary types, or changes
to the neighborhood/depth APIs scheduled for S2.

Seven new unit cases check empty/isolated graphs, interleaved component membership, loops,
parallel edges, and Break at the first/middle/last component. Three new properties compare both
event collectors with their visitors and definition-level references, and components with repeated
closure under edge incidence. Generated multigraphs are bounded at eight nodes and twenty edges;
the original DFS/BFS properties retain their laws and generators. The component property also
checks collector agreement and early-break prefixes. Source review confirms Break is propagated
from the component callback through the BFS FinishTree callback; S1b's instrumented callback tests
check that BFS then makes no further root/neighbor/visitor calls.

Validation passed: 1,003 graph-core unit tests, 54 integration tests, and 135 properties, plus all
three graph-IR component consumer cases. Strict graph-core Clippy with all targets and proptest
passed, and rustdoc built without warnings. Commands and logs:

```text
cargo test -p umol-graph-core --features proptest --offline
    scratch/s1c-tests.log
cargo test -p umol-graph-ir enumerate_connected_components --offline
    scratch/s1c-graph-ir-tests.log
cargo clippy -p umol-graph-core --all-targets --features proptest --offline -- -D warnings
    scratch/s1c-clippy.log
cargo doc -p umol-graph-core --no-deps --offline
    scratch/s1c-doc.log
```

All 24 benchmark measurements completed on the existing inline examples: six each for DFS event
collection, BFS event collection, component visitation, and the migrated component enumerator.
Criterion time point estimates below are in microseconds per operation:

| Example | DFS events | BFS events | Component visitor | Components |
| --- | ---: | ---: | ---: | ---: |
| path_64 | 1.5719 | 0.72706 | 0.55926 | 0.58502 |
| path_1024 | 15.270 | 7.3811 | 6.8372 | 6.6054 |
| binary_tree_255 | 2.7983 | 2.5289 | 2.2349 | 2.3918 |
| cycle_64 | 1.3698 | 0.70065 | 0.69514 | 0.72852 |
| four_hexagons | 0.45190 | 0.40137 | 0.24097 | 0.34119 |
| loops_parallel_isolated | 0.15954 | 0.16455 | 0.077054 | 0.14111 |

Compared with S0c, component collection has mixed timing changes: the 1,024-node path is
6.605 µs versus 6.565 µs, the binary tree is 2.392 µs versus 2.752 µs, and the small
loops/parallel/isolated graph is 0.141 µs versus 0.114 µs. These bounded runs provide baseline
evidence, not fine rankings; the DFS path_64 interval is notably wider at 1.391–1.894 µs.
No tuning or repeat campaign was performed. Reproduce with:

```text
cargo bench -p umol-graph-core --bench algorithms --offline --
    '(traversal/(depth_first_events|breadth_first_events|component_visitor)|traversal_baseline/components)'
    --sample-size 20 --warm-up-time 0.1 --measurement-time 0.2
    --nresamples 1000 --noplot --save-baseline s1c
```

The run log is scratch/s1c-core-bench.log. cargo fmt --all, its check mode, and git diff --check
passed. Final diff review covers the agreed public surface, direct visitor composition, unchanged
existing test laws, and inline benchmark additions. S1 is complete. S2a's neighborhood naming and
usize migration are next.

## S2a neighborhood and size migration — 2026-09-10

The public selector becomes NeighborhoodAlgorithm::Bfs, without a TraversalAlgorithm alias.
Graph::neighborhood takes max_depth: usize and returns Vec<(NodeId, usize)> by collecting Discover
events from one-root BFS with that limit. It preserves CSR neighbor order within equal-distance
shells, shortest distances, the source at depth zero, and exclusion of disconnected nodes.
Invalid source ids inherit the traversal kernel's panic-free handling. No new validator or error
type is introduced; Graph already owns consistent connectivity.

The related public size changes are CircularRefinementAlgorithm::Ec::radius,
CircularRefinementHash::combine's round argument, MorganFeaturizer::radius and ::new,
EcfpFeaturizer::radius and ::new, and Python HashedFingerprintConfig::{Morgan,Ecfp}::radius.
These remain open size/configuration values; internal round indices and duplicate-removal keys
also become usize. Graph node/edge ids, chemistry invariant components, bond labels, hash words,
and the independent WL RefinementRounds type retain their existing widths. GraphView exposes
no neighborhood method requiring migration.

Hash recipes retain their frozen encodings: Morgan's zero-based layer salt is a 32-bit hash word
(low 32 bits), while RogersHahn serializes the round as an eight-byte little-endian word. The
conversion at that hash-word boundary does not narrow the radius or iteration counter. Existing
fingerprint golden outputs must remain unchanged. Python radius constructors accept platform-size
nonnegative integers and preserve them through conversion to Rust; out-of-range integers are
rejected by the binding's usize conversion. Constructing a large radius does not execute refinement.

Tests retain existing fingerprint goldens, add ordered neighborhood and independent-distance
checks, exercise refinement beyond radius zero, and check hash encoding and Python size boundaries.
No new binding wrapper or chemistry interpretation is added.

Implemented against source commit 9c86723d26724be31f047d15921c682df60b21f4. Public-surface review
matches the selector, neighborhood signature, circular-refinement size fields/trait argument,
featurizer fields/constructors, and Python radius fields listed above. The circular fingerprint
benchmark's CIRCULAR_RADIUS constant also becomes usize. Source search finds no TraversalAlgorithm
or remaining affected u32 radius/depth/round signature; historical discussion records retain their
original names. The WL round type and subgraph-isomorphism depth counter are independent paths.

Neighborhood fixtures now assert exact returned order instead of sorting before comparison, and
include depth zero and usize::MAX. A generated property checks every source in bounded multigraphs
against both BFS event order and independently relaxed distances. Refinement fixtures retain
radius-zero behavior and check one/two-round duplicate removal. All existing Morgan/ECFP golden
identifiers are unchanged. Direct hash checks cover the first rounds and rounds above u32::MAX;
the ECFP expectation uses an explicit 16-byte round/current buffer. Python tests cover zero and
usize::MAX construction, negative and overflowing integers, and unchanged fingerprint results.

Validation passed:

| Crate | Library unit tests passed | Property tests passed | Additional validation |
| --- | ---: | ---: | --- |
| graph-core | 1,007 | 136 | 54 integration tests |
| graph | 1,066 | 5 | 54 integration tests and 21 binary tests |
| graph-IR | 6,816 | 373 | 39 integration tests, including the compile-fail suite |
| Python | 50 selected Rust config tests | — | 94 Python fingerprint tests against the rebuilt extension |

Graph-IR retains three ignored unit tests and one ignored property. The full three-crate test
gate passed; graph-core and graph were rerun after the final round-iteration lint correction.
Strict Clippy passed across graph-core, graph-IR, graph, and Python with all targets and the three
Rust property features enabled. Rustdoc for graph-core and graph built without warnings.

```text
cargo test -p umol-graph-core -p umol-graph-ir -p umol-graph
    --features umol-graph-core/proptest,umol-graph-ir/proptest,umol-graph/proptest --offline
    scratch/s2a-tests.log
cargo test -p umol-graph-core -p umol-graph
    --features umol-graph-core/proptest,umol-graph/proptest --offline
    scratch/s2a-final-core-graph-tests.log
cargo test -p umol-py --lib fingerprint::config --offline
    scratch/s2a-python-rust-tests.log
cargo clippy -p umol-graph-core -p umol-graph-ir -p umol-graph -p umol-py --all-targets
    --features umol-graph-core/proptest,umol-graph-ir/proptest,umol-graph/proptest --offline -- -D warnings
    scratch/s2a-clippy.log
cargo doc -p umol-graph-core -p umol-graph --no-deps --offline
    scratch/s2a-doc.log
maturin develop --manifest-path umol-py/Cargo.toml --offline
    scratch/s2a-maturin-final.log
pytest -q umol-py/tests/test_fingerprint.py
    scratch/s2a-python-tests.log
```

Every Python build/test command ran with umol-py/.venv activated (Python 3.13.15). The native
rebuild used UV_CACHE_DIR=/private/tmp/umol-s2a-uv-cache because the default uv cache was not writable.

All six existing inline neighborhood benchmarks completed with the usize limit and shared BFS.
Criterion time point estimates are in microseconds per operation:

| Example | S0c neighborhood | S2a neighborhood |
| --- | ---: | ---: |
| path_64 | 0.556 | 0.53127 |
| path_1024 | 6.119 | 5.7548 |
| binary_tree_255 | 2.060 | 2.2877 |
| cycle_64 | 0.556 | 0.51933 |
| four_hexagons | 0.158 | 0.11249 |
| loops_parallel_isolated | 0.086 | 0.070969 |

These bounded measurements retain the S0c workloads and establish the migrated baseline without
tuning. The full intervals are in scratch/s2a-core-bench.log. Reproduce with:

```text
cargo bench -p umol-graph-core --bench algorithms --offline -- traversal_baseline/neighborhood
    --sample-size 20 --warm-up-time 0.1 --measurement-time 0.2
    --nresamples 1000 --noplot --save-baseline s2a
```

cargo fmt --all, its check mode, and git diff --check passed. Final diff review confirms the
size migration is limited to the agreed path, with no compatibility alias, radius truncation,
fingerprint golden changes, or modifications to the independent WL size contract. S2 is complete;
S3a's explicit TableIR stereo-bond vocabulary and frame algebra are next.

## S3a stereo-bond vocabulary and frame algebra — 2026-09-10

Contract for this additive subitem:

- **Types and role:** StereoBond, BondConfiguration, and BondRelation are open TableIR carriers.
  Public fields and enum variants have exactly the settled shape above; the existing table_ir
  reexport exposes them. No constructor, validator, conversion, transformation method, or Python
  binding is added.
- **Intrinsic representation:** Either has no references; Framed has exactly two actual atom ids
  and a definite relation. No record means absence of an assertion, distinct from Either.
- **Context:** the site id indexes the owning molecule's bond table. Reference positions follow
  that bond's ordered AtomPair endpoints. Index validity, substituent incidence, and compatibility
  of operative assertions belong to the first consumer requiring those properties, introduced in
  later subitems. Bare construction does not inspect a molecule or certify chemistry.
- **Preservation and failures:** this subitem performs no conversion or normalization and adds no
  fallible operation. Renumbering must preserve reference identity, exchanging the reference slots
  if the site's ordered endpoints reverse. Bond-table renumbering transports the site id. Selecting
  minimum-index references is a producer policy, not a representation invariant.
- **Algebra:** with complete endpoint blocks `[a, b]; [c, d]` and references `[a, c]`, SameSide
  corresponds to ClassKey::CisTrans coset 0 and OppositeSide to coset 1. Changing one reference
  exchanges the cosets; changing both or exchanging complete endpoint blocks preserves them.
  These correspondences are explicit, independent of enum discriminants. Tests use the existing
  permutation/coset operations and independently specified target frames; they do not introduce
  a second production frame-transport implementation.

StereoAtom, Winding, and LonePair retain their existing definitions. Parser production, collection
storage, source-marker retirement, and raise migration remain S3b–S3e work.

Implemented against source commit c59421b1c82932b3989ef855d22353bed3e69015. The public-surface
review matches the three records/enums, their fields/variants, and the existing table_ir reexport;
there is no new conversion or remapping API. Rustdoc states endpoint order, site-id transport,
the distinction between absence and Either, contextual checks, and the explicit coset convention.

The 38 contract cases comprise all eight admissible endpoint-block actions on both configurations,
six atom/bond renumberings for Either and both definite configurations, and two inadmissible
cross-endpoint actions on both cosets. Exact target references and complete frames are independent
case data. Permutation/coset operations and GraphIR StereoBondForm frame transport are checked
against those expectations, with independent signed side labels for the reference relation.
Renumbering cases include reversed endpoints, changed site ids, and retained references that are
no longer minimum-indexed. These are vocabulary/algebra fixtures; integrated producer and consumer
transport is exercised when those paths adopt the records. No executable algorithm is added in
S3a, so no new benchmark is needed for this subitem.

Validation passed:

```text
cargo test -p umol-io --lib table_ir::stereo --offline
cargo test -p umol-io -p umol-graph
    --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline
cargo clippy -p umol-io -p umol-graph --all-targets
    --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline -- -D warnings
RUSTDOCFLAGS='-D warnings' cargo doc -p umol-io --no-deps --offline
cargo fmt --all -- --check
git diff --check
```

The final focused run passes all 38 frame-contract cases. The broader gate includes IO unit,
SMILES property, SMILES/MOL/SDF conformance, and dependent graph unit/integration/property tests.
Full diff review confirms only the agreed vocabulary, documentation, and contract tests changed;
the existing atom-stereo types and parser/raise paths are unchanged. S3b source direction
normalization is next.

## S3b source direction normalization — 2026-09-10

The derivation kernel stays private to the SMILES parser. Its input is a borrowed sequence of
completed source bonds (ordered AtomPair, order, optional lexical direction) and an atom count;
its output is StereoBond records in table-bond order. Directions are viewed from AtomPair::first,
as established by the existing builder's ring reconciliation and bond construction. Input markers
are not changed. No new persistent table field, public constructor, or public conversion is added.

Temporary incidence lookup belongs to this operation. It selects the minimum actual reference
at each endpoint, consumes every incident single-bond marker, and checks endpoint agreement before
deciding whether both endpoints determine a frame. Consistent partial input supplies no assertion.
Each marker must belong to an eligible local double-bond site; shared markers can participate in
several sites. Failures retain the responsible atom or bond id in private derivation diagnostics.
Malformed indices are checked before accessing endpoint neighborhoods. A marked site outside the
local two-ligand-per-endpoint domain fails explicitly rather than discarding evidence or selecting
two ligands from a larger set. Marked cumulated axes likewise fail rather than being treated as
redundant local partial notation. This is frame representability, not chemical valence resolution.

S3b tests the kernel before publication. Parser acceptance, CX updates, and raise remain unchanged;
S3c handles annotation/geometry derivation and S3d wires producers and maps failures into boundary
diagnostics. The temporary dead-code expectation is confined to this unconnected parser module
and is removed when it is wired. Benchmarks use inline examples in scratch to call this private
kernel without adding a public benchmark-only API.

Implemented against source commit d9950f36e51c36f12013f8489365c0765df8902a in
smiles/parser/stereo.rs, with separate unit and property modules. The builder's opening-order
bond slots, ring direction reconciliation, basic/extended conversion, CX application, and existing
raise functions are unchanged. The kernel borrows transient source triples; it adds no second
persistent direction authority and does not assign H, lone pairs, coordinates, or molecular states.

The final focused run passes 92 tests. Exact-frame fixtures cover the recorded direction examples
through both parser targets, complementary references, shared chains/branches/cycles, partial
trienes, explicit-H markers, and opening/closing/both-end ring spellings. Raw source fixtures check
opening versus completion order and reversed bond storage, invalid endpoints, self-incidence,
and conflicting parallel markers. Unsupported marked cumulenes fail explicitly. An independent
half-plane oracle exhausts all 81 absent/rising/falling assignments on four substituent bonds,
including contradictory and partial assignments. A generated property checks atom/bond
renumbering, new minimum references, endpoint reversal, global marker reversal, and partiality
against physical side assignments.

Validation passed:

```text
cargo test -p umol-io -p umol-graph
    --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline
cargo test -p umol-io --lib smiles::parser::stereo --features proptest --offline
cargo clippy -p umol-io -p umol-graph --all-targets
    --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline -- -D warnings
cargo fmt --all -- --check
git diff --check
```

The broad gate includes 10,223 SMILES, 2,253 MOL, and 407 SDF conformance cases, the existing
SMILES properties, and dependent graph unit/integration/property tests. The final focused and
Clippy runs also include the cumulene checks added during diff review. No existing assertion or
property law was weakened, and no production acceptance/error-layer change is published here.

The standalone scratch harness borrows the private kernel by source path and contains eight
inline SMILES examples. Parsing and source-triple adaptation occur outside the timed loop;
temporary incidence construction, frame derivation, and result destruction are measured.
These are kernel measurements, not end-to-end parser timings or an optimization gate. Median
microseconds per call over seven batches of 10,000 calls:

| Example | Time (µs) | Frames |
| --- | ---: | ---: |
| chain_64 | 0.0231 | 0 |
| alkene_four | 0.2472 | 1 |
| shared_chain | 0.2770 | 2 |
| partial_triene | 0.3726 | 1 |
| shared_branch | 0.3562 | 2 |
| shared_cycle | 0.3357 | 2 |
| ring_closure | 0.2139 | 1 |
| redundant | 0.2426 | 1 |

```text
cargo run --release --manifest-path scratch/s3b-direction-bench/Cargo.toml --offline
```

Final surface/diff review confirms one private derivation function and its private diagnostic
enum, without a public helper or changes to the atom-stereo vocabulary. S3c CTfile and annotation
frame derivation is next; publication and retirement of duplicate fields remain S3d.

## S3c CTfile and annotation frame derivation — 2026-09-10

Contract: a private shared IO kernel consumes completed bond endpoints/orders/wedges, optional
borrowed positions, and the ordered list of explicit bond-code annotations. It produces the
existing StereoBond vocabulary in table-bond order. No public constructor, conversion, error,
or additional persistent field is introduced. The producer remains unconnected until S3d.

CTfile code 0 supplies geometry evidence; code 3 (V2000) or CFG=2 (V3000) supplies site-only Either.
The latter is a format contract, not implemented reader coverage: current counts parsing accepts
V2000 only. S3c tests the existing V2000 readers and decoded annotation semantics; it does not
introduce a V3000 parser.
Either requires no references or positions and suppresses coordinate interpretation at that site.
The existing narrow-end wavy-bond convention supplies Either only when that endpoint has exactly
one double-bond partner. Definite Up/Down wedges remain atom-stereo data. Missing or unusable
geometry supplies no assertion. Local definite frames use minimum-index actual references and
the existing same_side_of_axis predicate, including its 3D and relative-tolerance convention.

CX Cis/Trans become SameSide/OppositeSide in the same selected reference frame; Either remains
site-only. The local RDKit CXSmilesOps.cpp parse_doublebond_stereo path requests CX ordering;
cxsmiles_test.cpp's regression cases around lines 1547–1597 pin lowest-index references for both
ring and acyclic examples. Current umol stores Cis/Trans but raise only interprets Either and
geometry/directions. S3c implements their explicit frame meaning without changing current parse
acceptance. CX completion-order bond ids must be remapped before invoking the kernel.

Annotations are consumed as a list so conflicting repeated codes cannot disappear through field
overwrite. Repeated equal evidence is redundant; conflicting explicit evidence fails with the
site id. Explicit definite evidence requires a representable local frame. Usable geometry must
agree with a definite annotation; absent/degenerate geometry does not invalidate that independent
annotation. Contextual index/position checks occur where their values are required. No coordinates
are generated or changed, and no chemistry, CIP ranking, or stereogenicity judgment is added.

Implemented in table_ir/stereo/derive.rs, with exact-frame tests and a generated similarity law
beside the kernel. Tests cover two, three, and four actual substituents; supplied 2D/3D positions;
absent, zero, degenerate, non-finite, and very large coordinates; annotation conflicts; narrow-end
wavy evidence; and contextual index errors. Temporary coordinate scaling bounds intermediate
products before invoking the existing geometric predicate. Basic and extended V2000/CX readers
are exercised through their decoded tables, including CX ring-completion bond-index remapping.

Verification on base 30e1ce58e2e1a571d48ad5b375d68a08f72c5e53:

- All 67 focused derivation tests pass, including the generated law for signed coordinate
  permutations, translation, positive scaling, and atom/bond reordering.
- The IO/graph gate passes with properties and IO conformance enabled: 3,799 IO unit tests,
  1,066 graph unit tests, 10,223 SMILES conformance cases, 2,253 MOL cases, 407 SDF cases,
  and the remaining integration/property suites. One existing graph doctest remains ignored.
- Strict Clippy for both crates and all targets with the same features passes, as do formatting
  and diff checks. The public-symbol review confirms no added public API or carrier field.

The release-mode kernel baseline uses eight inline inputs in scratch/s3c-stereo-bench, with seven
batches of 10,000 calls per input. Median times include temporary lookup allocation and result
destruction, excluding input construction; these are kernel measurements, not parser throughput.

| Input | Frames | Median µs/call |
| --- | --- | --- |
| Absent positions | 0 | 0.1409 |
| All-zero positions | 0 | 0.2840 |
| Same-side 2D | 1 | 0.3456 |
| Opposite-side 2D | 1 | 0.3510 |
| Opposite-side 3D | 1 | 0.3289 |
| Explicit Either | 1 | 0.2210 |
| Explicit Cis | 1 | 0.2832 |
| Adjacent wavy bond | 1 | 0.2041 |

Commands and logs are retained under scratch/s3c-*. The principal gate is:

```sh
cargo test -p umol-io -p umol-graph --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline
cargo clippy -p umol-io -p umol-graph --all-targets --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline -- -D warnings
cargo run --release --manifest-path scratch/s3c-stereo-bench/Cargo.toml --offline
```

S3d is next: publish the frame collections, wire producers, preserve repeated annotation evidence
before current field overwrites, and migrate raise/consumers while retiring duplicate fields.

## S3d implementation record — published bond frames (2026-09-10)

### Contract and public surface

Molecule and ExtendedMolecule now carry a stereo_bonds collection of StereoBond records. They remain open tables;
empty construction initializes an empty collection, basic-to-extended conversion preserves it,
and extended-to-basic conversion copies it without reinterpreting coordinates. No arbitrary-parts
constructor for Smiles or ReactionSmiles was added. The S3a StereoBond, BondConfiguration, and
BondRelation vocabulary is unchanged.

Bond and ExtendedBond no longer have direction or stereo fields. BondDirection and BondStereo
remain source vocabulary used by parsers; they are not a second persistent configuration authority.
Wedges, StereoAtom, Winding, and LonePair retain their separate existing roles. TableIR stores no
Graph or adjacency structure.

SMILES ParseError now owns DanglingBondDirection and CisTransConflict, previously raise errors,
and adds UnsupportedStereoBond, ConflictingBondConfiguration, and MissingPosition. CTfile
ParseError adds the latter three variants. Existing boundary index errors carry derivation index
failures. RaiseError adds StereoBondIndexOutOfBounds, UnsupportedStereoBond,
InvalidStereoBondReference, and DuplicateStereoBond for independently supplied open frames.
These checks occur when source evidence is interpreted or a frame is consumed, respectively;
construction and mutation do not validate the entire table.

### Producer and consumer migration

- The SMILES builder holds directional markers beside pending bond records until ring slots are
  complete, normalizes endpoint viewpoints, and publishes the S3b-derived frames. Consistent
  partial markers leave no assertion. Marker spellings are then discarded.
- The V2000 readers retain bond codes beside decoded records, apply property records, and derive
  frames from the completed tables and supplied coordinates. SDF follows the same molecule
  construction. This does not add V3000 parsing or coordinate generation.
- CX processing preserves repeated Cis/Trans/Either annotations as evidence rather than
  overwriting a bond field. Parser-produced direction frames use the same minimum-index
  references as CX codes, allowing the two sources to be checked together with geometry.
  Wavy-bond evidence is retained even if a later entry overwrites that wedge. Reaction CX splitting
  remaps indices before each section is normalized.
- Raise consumes the frame collection. Framed transports the selected references into the sorted
  endpoint blocks with the existing cis/trans permutation action; Either produces an undetermined
  #C assertion. It does not interpret coordinates or direction markers for bond stereo. Frame sites
  retain all-table-bond indices while localized graph bonds and relations are partitioned.
  Complementary virtual slots are not assigned a hydrogen or lone-pair identity by this step.
- Parser targets, empty tables, conversions, ingest tests, table-remapping tests, and benchmark
  inputs now use frames. Expected-table test builders describe bond orders and assign expected
  frames explicitly; they no longer accept a direction that is silently ignored.

### Reviewed acceptance changes

The existing 10,223-case SMILES corpus has eleven changed classifications. Inputs and source
comments are preserved byte-for-byte; only outcome categories and their reviewed snapshots change.
The corpus category invalid denotes rejection by the current parser, including unsupported stereo,
not a claim that every rejected string is universally invalid SMILES.

| Source cases | New outcome | Reason |
| --- | --- | --- |
| NextMove 0007, 0008 | chemaxon_invalid | CX Either conflicts with definite direction evidence at bond 3. OpenSMILES parsing still succeeds when CX is not selected. |
| OpenBabel 0105 | invalid | Direction at bond 1 has no supported local double-bond frame. |
| Indigo 2163 | invalid | Contradictory endpoint markers at atom 9. |
| Indigo 0173, 1442; RDKit 0528, 0611, 0637, 0654, 0704 | invalid | Marked cumulated local axes are outside the S3b-supported domain. |

The cumulated-axis group includes shared markers that also touch an ordinary alkene. This follows
the explicit S3b rejection rule; the unsupported interpretation is not silently dropped. The exact
inputs and diagnostics are retained in scratch/s3d-corpus-changes.tsv and the corpus snapshots.

### Verification and measurements

Exact-frame tests cover basic/extended SMILES, V2000 and CX producers, reaction sections,
ring-completion remapping, repeated and contradictory annotations, wavy evidence, and source
errors. Raise tests cover reference swaps, endpoint/atom/bond remapping, site-only Either,
malformed indices and incidences, duplicates, and all-table indices with a preceding dative bond.
Generated laws exercise reference-swap transport, exact basic/extended preservation, and coordinate
edits leaving published bond configurations unchanged. Arbitrary open bond frames are checked for
absence of panics. Existing atom-stereo coverage remains in place.

On base 4d3366bed05c71953b60ad5ad687291c8ea2210a plus S3d:

- The IO/graph gate passes with proptest and IO conformance enabled: 1,066 graph unit tests,
  10,223 SMILES cases, 2,253 MOL cases, 407 SDF cases, and the remaining integration/property
  suites. The final IO unit run has 3,800 passing tests. One existing graph doctest is ignored.
- Strict Clippy passes for both crates, all targets, and those features. Formatting and diff checks
  pass. The fuzz_parse_opensmiles target builds; no additional timed fuzz campaign was run here.
- Every smiles_parsing benchmark passes Criterion's test mode. Three old inputs with dangling
  directions are retained as explicit basic/extended rejection benchmarks; success groups use
  supported inline replacements. No fixture files were introduced into benchmarks.

```sh
cargo test -p umol-io -p umol-graph --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline
cargo clippy -p umol-io -p umol-graph --all-targets --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline -- -D warnings
cargo check --manifest-path umol-io/fuzz/Cargo.toml --bin fuzz_parse_opensmiles --offline
cargo bench -p umol-io --bench smiles_parsing --offline -- --test
cargo bench -p umol-io --bench smiles_parsing --offline -- smiles_roundtrip --sample-size 20 --warm-up-time 0.1 --measurement-time 0.2 --nresamples 1000 --noplot
```

The final bounded run uses the same sixteen inline S0c inputs and timing boundaries, with
20 samples, 0.1-second warm-up, 0.2-second measurement, and 1,000 resamples. Compilation and
other verification completed before measurement. Criterion point estimates are shown below;
intervals are retained in scratch/s3d-timings.csv and scratch/s3d-bench-final.log.

| Example | Parse µs | Raise µs |
| --- | ---: | ---: |
| chain_64 | 1.029 | 12.484 |
| branched | 0.296 | 1.909 |
| components | 0.207 | 1.018 |
| aromatic_fused | 0.346 | 2.889 |
| aromatic_lone_pair | 0.266 | 1.677 |
| charged | 0.169 | 0.584 |
| radical | 0.166 | 0.570 |
| tetra_four | 0.327 | 2.200 |
| tetra_ring | 0.450 | 2.576 |
| tetra_explicit_h | 0.411 | 2.886 |
| tetra_lone_pair | 0.288 | 2.204 |
| alkene_four | 0.545 | 1.938 |
| shared_chain | 0.533 | 2.136 |
| partial_triene | 0.688 | 2.206 |
| shared_branch | 0.689 | 2.404 |
| shared_cycle | 0.679 | 2.645 |

Source normalization now runs during parsing; raise transports the published frames. These short
measurements describe that division of work and do not establish fine performance rankings.

### Allocation review and revised closeout (2026-09-11)

The S3d correctness gate passed, but its completion claim was premature: the implementation added
full-bond copies to fit shared signatures and retained avoidable allocation in the connected paths.
S3d1–S3d7 now follow S3d in the plan; S3e continues to own conditional temporary neighbor lookup.
The earlier test and timing results describe the implementation reviewed, not acceptance of its
allocation design. No allocation corrections have been implemented by this plan update.

The settled constraints are direct borrowed access or separate basic/extended methods, acceptance
of deliberate duplication, movement of consumed storage, and no persistent graph in TableIR.
Choosing how to factor private functions is an implementation responsibility within those
constraints; it is not an open question that requires another abstract API design. Existing public
boundary types, methods, stereo semantics, and failure ownership remain unchanged.

| Audit finding | Disposition/owner |
| --- | --- |
| Clones in consuming extended-to-basic conversion, including stereo_bonds | S3d1: move fields and consume elements. |
| CTfile tuple staging and CTfile/CX records copies | S3d2: accumulate final bonds plus stereo assertions; derive through borrowed access. |
| SMILES completed/lexical intermediates and bond-sized participation flags | S3d3: remove copies and consolidate transient marker bookkeeping. |
| Dense code array and expansion/reconstruction of existing CX frames | S3d4: group sparse assertions and update/check existing frames directly. |
| Repeated endpoint and blocks heap vectors | S3d5: bounded local scratch for the supported frame domain. |
| Raise's temporary constraint HashMap | S3d6: ordered borrowed frames consumed at bond construction. |
| Reservations proportional to unrelated trailing input | S3d7: one-pass growth with a tunable initial reservation. |
| Nested/repeated/unconditional AtomNeighbors allocation | S3e: conditional operation-local lookup, retaining table incidence semantics. |
| Final table atoms/bonds, supplied positions, published frames, and GraphIR result | Retain required output storage; eliminate redundant copies, not the representation. |

### Capacity policy for S3d7 (settled 2026-09-11)

Use one-pass parsing and let the owned atom/bond tables grow as entries are produced. Initial
reservation remains a tunable implementation estimate; zero initial capacity is not required.
The current reservation overestimates entry counts by using the remaining input's byte length,
including bracket syntax, CX annotations, and subsequent reaction sections. Improve that estimate
without making exact capacity prediction a prerequisite for parsing or the other allocation fixes.

A preliminary atom/bond counting pass is not selected. Growth reallocations are an accepted
tradeoff, and tuning can refine the initial reservation later. This policy introduces no input-size
limit, algorithm selector, or public configuration. Specification edits remain separately staged
below.

### S3d1 pre-change measurements (2026-09-11)

The baseline measures the existing consuming ExtendedMolecule-to-Molecule conversion, before
S3d1 code changes. The standalone probe is scratch/s3d1-conversion-bench; its four examples are
inline and use the public conversion. stereo_positions has six supplied positions, one atom
frame, and one bond frame. stereo_metadata adds atom labels/values, a comment, and a property to
that case. multicenter contains two multicenter bonds. Preparation verifies complete equality
against the basic input; every measured result is checked outside the measurement boundary.

Allocation measurement reuses the S0c counting System allocator, enabled only by the scratch
crate's allocations feature. After warm-up, input cloning/preparation occurs before counters reset;
conversion consumes that input, and the result remains alive while counters are read. Calls include
successful alloc, alloc_zeroed, and realloc; requested bytes include the full new size for every
reallocation. Peak added live bytes are above the operation-entry baseline, which already includes
the input. They exclude allocator overhead and are not process RSS.

Timing is a separate build with the normal allocator, using eleven batches of 1,024 conversions
per case after 64 warm-up conversions. Input batches and output vector capacity are prepared
outside timing; the timed loop converts inputs and stores outputs. Input destruction intrinsic
to conversion is included; result destruction, result verification, and batch-buffer destruction
are excluded. The table reports the median of batch means and their observed range, not a
statistical confidence interval. No timing repeats or tuning campaign were performed.

| Inline case | Atoms / bonds | Allocation calls (reallocations) | Requested bytes | Peak added live bytes | Median µs / conversion | Batch-mean range µs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| chain_64 | 64 / 63 | 10 (8) | 17,856 | 9,216 | 4.742 | 4.583–5.228 |
| stereo_positions | 6 / 5 | 8 (2) | 1,952 | 1,376 | 0.669 | 0.660–0.734 |
| stereo_metadata | 6 / 5 | 26 (2) | 2,330 | 1,754 | 0.998 | 0.988–1.040 |
| multicenter | 11 / 11 | 7 (4) | 4,448 | 2,720 | 1.094 | 1.088–1.119 |

The chain's ten calls consist of two initial output-vector allocations and eight growth
reallocations. The additional calls for positions, stereo, and metadata expose the clone work
that S3d1 is intended to remove. These numbers are the baseline, not an after-change claim.

Environment: rustc 1.96.0 (ac68faa20), arm64 macOS 15.7.3, standalone Cargo release profile.
Base revision: 4d3366bed05c71953b60ad5ad687291c8ea2210a plus the uncommitted S3d work.
Measured table_ir/molecule.rs SHA-256:
6d79906050d9cd410e2edbec7d267e147b3003c09fb24e1924bd2f653c9fb5d6.
Raw output is scratch/s3d1-before-allocations.csv and scratch/s3d1-before-timings.csv; corresponding
build logs use the same prefixes. Reuse this probe unchanged for the S3d1 after measurement.

```sh
cargo run --release --manifest-path scratch/s3d1-conversion-bench/Cargo.toml --offline --features allocations
cargo run --release --manifest-path scratch/s3d1-conversion-bench/Cargo.toml --offline
```

### S3d1 implementation and after measurements (2026-09-11)

Molecule::try_from(ExtendedMolecule) now consumes the atom and bond vectors through their existing
fallible element conversions. Positions, multicenter bonds, both stereo collections, comments,
and properties move into the result. The signature, ConversionError, atom/bond rejection rules,
open-table semantics, and basic/extended preservation law are unchanged. No constructor, borrowed
conversion, helper type, or other public symbol was added. CTfile/CX-only fields retain their
existing conversion treatment.

The exact conversion tests retain atom-feature rejection and now cover bond-feature rejection.
The roundtrip case additionally carries nonempty positions, a bond frame, multicenter data,
comments, and properties. Existing atom-frame cases and generated frame/coordinate preservation
laws remain in place. The code review confirms that the conversion contains no field or element
cloning; allocation measurements below verify the ownership change independently of value equality.

The same four inline inputs and measurement boundaries were used for one after allocation run
and one normal-allocator timing run. Preparation and verification remain outside measurement.

| Case | Allocation calls before → after | Requested bytes before → after | Peak added live bytes before → after | Median µs before → after | After batch-mean range µs |
| --- | ---: | ---: | ---: | ---: | ---: |
| chain_64 | 10 → 6 | 17,856 → 23,368 | 9,216 → 2,536 | 4.742 → 1.831 | 1.679–2.108 |
| stereo_positions | 8 → 3 | 1,952 → 2,144 | 1,376 → 256 | 0.669 → 0.267 | 0.262–0.288 |
| stereo_metadata | 26 → 3 | 2,330 → 2,144 | 1,754 → 256 | 0.998 → 0.267 | 0.262–0.280 |
| multicenter | 7 → 4 | 4,448 → 4,240 | 2,720 → 592 | 1.094 → 0.448 | 0.418–0.541 |

After-run reallocations are 5, 2, 2, and 3 respectively; the remaining one call per case is a
fresh allocation. The metadata case now has the same allocation cost as the coordinate/stereo
case. The short timing run shows lower conversion time for all four inputs, without establishing
performance across arbitrary molecule sizes or metadata distributions.

There is a capacity tradeoff: the owning atom conversion reuses the larger ExtendedAtom buffer.
A single diagnostic trace of the 64-atom case shows input atom capacity 64 at 288 bytes per element
and output atom capacity 177 at 104 bytes per element. The allocator receives a shrink request
from 18,432 to 18,408 bytes. The bond result grows through 160, 320, 640, 1,280, and 2,560 bytes.
This explains the larger cumulative requested-byte total despite fewer calls and a smaller
additional live-memory peak. It also means the result retains more spare atom capacity; the
smaller peak-above-input metric must not be read as a smaller retained result. No shrink-to-fit
policy or alternative collection strategy is introduced in S3d1. The trace is recorded in
scratch/s3d1-allocation-trace.log and its separate allocation_trace binary; the comparison probe's
measurement code is unchanged.

Verification passes: 32 focused molecule tests; the IO/graph gate with proptest and IO conformance
(3,801 IO unit tests, 1,066 graph unit tests, 10,223 SMILES cases, 2,253 MOL cases, 407 SDF cases,
and the remaining integration/property suites); and strict Clippy for both crates/all targets with
those features. One existing graph doctest remains ignored. Formatting and diff checks pass.
The incremental source/test diff was reviewed against S3d1; other S3d work was preserved.

After raw results are scratch/s3d1-after-allocations.csv and scratch/s3d1-after-timings.csv.
The measured table_ir/molecule.rs SHA-256 is
bb26ebc70a7bdd6b7037d7ef01bb87d99a52df30fb2ba8f87d4ccd9819866bb9.
S3d2 is next; no other allocation subitem was implemented here.


### S3d2 implementation and measurements (2026-09-11)

Both V2000 bond blocks accumulate final Bond/ExtendedBond vectors plus separate
bond_stereo_assertions lists. Each assertion carries its original table-bond index and source
stereo code; absent codes add no entry. The builders move those vectors into their tables,
apply properties, then derive stereo from the resulting bond storage. Neither block stages a
full vector of bond/code tuples, and neither builder recollects the bonds.

CTfile and CX basic/extended callers now pass borrowed bond slices into the existing derivation
kernel. A statically dispatched field callback reads each bond's endpoints, order, and wedge;
there is no copied records vector, dynamic dispatch, or adapter type. The four production
callbacks read those three fields directly. CX's existing frame-to-code reconciliation remains
unchanged, and its code-bearing lists also use the name bond_stereo_assertions. The dense codes
array, local substituent allocations, and temporary neighbor storage remain for their owning
follow-ups; this subitem does not claim to eliminate all derivation allocations.

Public-symbol reconciliation: no public type, constructor, conversion, signature, visibility,
error variant, or Python surface changed. TableIR remains an open table carrier. Parser ordering,
property precedence, wedge viewpoint, CX completion-index remapping, frame reference selection,
and derivation failure behavior are preserved. Only private bond-block/builder and derivation
signatures changed. Existing property assertions and generated domains are unchanged.

Exact tests now cover sparse assertions at bond indices 1 and 3, an empty bond block, Either
codes, and full basic/extended bond values. Additional end-to-end V2000 cases show that ZBO
promotion/demotion changes geometry-derived frames, unchanged double-bond order preserves Either,
and demotion with an explicit Either code reports UnsupportedStereoBond. Those cases also check
charge/isotope property preservation. Existing CX ring/completion-index, wavy, geometry, conflict,
and coordinate-similarity tests continue to exercise the borrowed path.

One before and one after run use scratch/s3d2-parser-bench with four inline examples, each through
basic and extended parsing: a 64-atom V2000 chain; a four-atom V2000 alkene with Either, CHG, and
ISO properties; a 64-atom CX chain with an atom label; and a CX ring with an explicit cis code.
V2000 text is generated inline before measurement; no fixture files are loaded. Each run uses
64 warm-ups. Allocation counting covers one complete parse while retaining the result. Timings
use the normal allocator and 11 batches of 512 parses, retaining results in preallocated output
vectors until after timing. Input/configuration construction, verification, result destruction,
and output-vector allocation are excluded. This measures parsing into TableIR, not raise/resolve.
The probe is unchanged between the two runs, which ran separately from tests and compilation of
other targets. The source baseline includes completed S3d1 and the existing uncommitted S3d work;
its affected files are retained under scratch/s3d2-before.

Allocation calls include reallocations; requested bytes sum complete new allocation sizes,
including realloc requests. Peak added live bytes are tracked relative to the pre-parse baseline,
not RSS or retained-result size. Ranges below are batch-mean ranges, not confidence intervals.

| Case | Allocation calls before → after | Requested bytes before → after | Peak added live bytes before → after | Median µs before → after | Before / after batch-mean ranges µs |
| --- | ---: | ---: | ---: | ---: | --- |
| mol_chain_64_basic | 72 → 70 | 18,199 → 14,431 | 15,427 → 14,431 | 10.565 → 11.281 | 10.502–11.978 / 10.854–12.068 |
| mol_chain_64_extended | 72 → 70 | 37,523 → 29,735 | 30,971 → 29,735 | 13.433 → 13.583 | 12.269–14.624 / 12.511–14.161 |
| mol_stereo_properties_basic | 18 → 16 | 3,531 → 3,363 | 3,072 → 3,072 | 1.187 → 1.145 | 1.147–1.329 / 1.119–1.181 |
| mol_stereo_properties_extended | 18 → 16 | 4,615 → 4,267 | 3,976 → 3,976 | 1.243 → 1.231 | 1.233–1.401 / 1.194–1.359 |
| cx_chain_64_basic | 74 → 73 | 18,709 → 17,953 | 14,869 → 14,242 | 2.993 → 2.773 | 2.828–3.392 / 2.692–3.392 |
| cx_chain_64_extended | 75 → 74 | 43,173 → 42,417 | 32,053 → 31,426 | 4.894 → 4.948 | 4.465–5.540 / 4.177–5.768 |
| cx_ring_basic | 31 → 30 | 10,254 → 10,086 | 7,656 → 7,656 | 1.240 → 1.263 | 1.211–1.430 / 1.198–1.623 |
| cx_ring_extended | 31 → 30 | 18,350 → 18,182 | 14,192 → 14,192 | 1.453 → 1.393 | 1.422–1.934 / 1.378–1.530 |

CTfile removes two allocation calls in each case: the tuple-to-bond collection's realloc and
the copied derivation table allocation. CX removes one copied-table allocation. Requested bytes
fall in every case; peaks fall for the chains and are unchanged for the smaller examples.
Runtime results are mixed: the basic 64-atom MOL median increases about 6.8%, while several
small/CX cases improve. The short runs and overlapping ranges do not establish a general speedup
or a stable regression. No further tuning or benchmark campaign was undertaken.

Verification passes: 3,815 IO unit tests and 1,066 graph unit tests with proptest enabled;
10,223 SMILES, 2,253 MOL, and 407 SDF conformance cases; the remaining integration/property suites;
and strict IO/graph Clippy for all targets with those features. One existing graph doctest is
ignored. The final bond-test assertion review was followed by a focused run of all 115 CTfile
bond tests. Formatting and diff checks pass. The incremental diff was reviewed against S3d2;
pre-existing S3d/S3d1 work was preserved.

Commands:

```sh
cargo test -p umol-io -p umol-graph --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline
cargo test -p umol-io --lib ctfile::parser::bond --offline
cargo clippy -p umol-io -p umol-graph --all-targets --features umol-io/proptest,umol-io/conformance,umol-graph/proptest --offline -- -D warnings
cargo run --release --manifest-path scratch/s3d2-parser-bench/Cargo.toml --offline --features allocations
cargo run --release --manifest-path scratch/s3d2-parser-bench/Cargo.toml --offline
```

Raw results are scratch/s3d2-{before,after}-{allocations,timings}.csv, with corresponding build
logs. Verification logs use scratch/s3d2-{gate,bond-tests,clippy,fmt}.log. S3d3 is next.


### S3d3 implementation and measurements (2026-09-11)

SMILES finalization derives bond frames directly from the builder's bond table after the existing
unclosed-ring/branch checks. One consuming collection then moves the final bonds into table order.
The completed intermediate and copied lexical bond table are removed. The statically dispatched
field callback reads original bonds and borrows their transient markers; it creates no adapter
collection or dynamic-dispatch boundary.

DirectionMarker is a private parser value containing a direction and a Cell<bool> participation
flag. Both ordinary bond insertion and ring completion create fresh markers. Derivation reads
the direction unchanged at every sharing site, marks participation, then rejects dangling markers
after the full scan. Cell permits participation updates while topology and marker values remain
borrowed; it requires no separate bond-sized flag vector. Markers are discarded with the temporary
builder storage, and the derivation contract requires fresh markers. No public type, constructor,
conversion, signature, error, or Python surface changed. Existing ring-slot completion establishes
the internal expect precondition; no new public checks or acceptance rules were introduced.

Atom-stereo incidence order, normalized endpoint viewpoints, ring-opening bond-table order,
CX completion-index remapping, partial/conflicting-marker handling, and frame reference selection
are preserved. Existing exact, exhaustive 81-assignment, and generated remapping tests retain
their assertions and domains. A new basic/extended regression checks a later conflict after a
shared marker has already participated at an earlier double bond. Neighbor construction and
per-endpoint substituent vectors remain for the later subitems.

The standalone scratch/s3d3-parser-bench uses the S3d2 measurement method: separate release builds
for allocation counting and normal-allocator timing, 64 warm-ups, then one allocation observation
or 11 batches of 512 parses. Parsing results remain alive until measurement ends; input/config
construction, output-vector allocation, verification, and result destruction are outside timing.
One before and one after run use the same six inline inputs through basic and extended parsing:
a bare 64-atom chain, a 42-atom chain with 20 marked double bonds, a triene, a shared branch,
a ring with atom and bond stereo, and a CX cis ring without lexical direction markers. The longer
directional string is built inline before timing. No fixture files are read. Tests and other
compilation were excluded from the measurement runs. The baseline includes S3d2; affected source
files are retained under scratch/s3d3-before.

Allocation calls include reallocations. Requested bytes sum full allocation/reallocation requests;
peak added live bytes are allocator-tracked above the pre-parse baseline, not RSS. Timing ranges
are batch-mean ranges, not confidence intervals.

| Case | Allocation calls before → after | Requested bytes before → after | Peak added live bytes before → after | Median µs before → after | Before / after batch-mean ranges µs |
| --- | ---: | ---: | ---: | ---: | --- |
| chain_64_basic | 3 → 3 | 12,188 → 12,188 | 9,428 → 9,428 | 0.998 → 0.994 | 0.942–1.683 / 0.979–1.828 |
| chain_64_extended | 3 → 3 | 31,512 → 31,512 | 24,984 → 24,984 | 3.182 → 2.935 | 2.877–4.194 / 2.549–3.579 |
| directional_42_basic | 92 → 90 | 20,325 → 19,792 | 15,669 → 15,136 | 2.915 → 3.046 | 2.826–3.675 / 2.862–3.849 |
| directional_42_extended | 92 → 90 | 45,365 → 44,832 | 35,861 → 35,328 | 5.972 → 5.244 | 4.722–6.582 / 5.111–5.895 |
| triene_basic | 21 → 19 | 3,475 → 3,384 | 2,811 → 2,720 | 0.674 → 0.642 | 0.660–0.825 / 0.634–0.804 |
| triene_extended | 21 → 19 | 7,915 → 7,824 | 6,411 → 6,320 | 0.810 → 0.763 | 0.804–1.166 / 0.752–1.300 |
| shared_branch_basic | 21 → 19 | 3,839 → 3,748 | 3,151 → 3,060 | 0.665 → 0.696 | 0.651–0.875 / 0.631–0.816 |
| shared_branch_extended | 21 → 19 | 8,579 → 8,488 | 6,995 → 6,904 | 0.886 → 0.732 | 0.789–1.323 / 0.728–0.811 |
| stereo_ring_basic | 19 → 17 | 3,959 → 3,868 | 3,287 → 3,196 | 0.613 → 0.610 | 0.604–0.774 / 0.583–0.757 |
| stereo_ring_extended | 19 → 17 | 8,699 → 8,608 | 7,131 → 7,040 | 0.760 → 0.734 | 0.748–0.836 / 0.686–1.074 |
| cx_ring_basic | 30 → 30 | 10,086 → 10,086 | 7,656 → 7,656 | 1.294 → 1.352 | 1.217–1.468 / 1.208–1.587 |
| cx_ring_extended | 30 → 30 | 18,182 → 18,182 | 14,192 → 14,192 | 1.412 → 1.352 | 1.379–1.878 / 1.330–1.918 |

The directional cases remove two allocations: the lexical table and participation vector. Their
requested and peak bytes fall by 13 bytes per bond on this build: 533 bytes for the 41-bond chain
and 91 bytes for the seven-bond cases. The completed collection had already reused its input
buffer; removing it does not eliminate an allocation call. Unmarked controls have identical
allocation counts and requested/peak bytes. Reallocation counts are unchanged in all cases.

Assessment against the three continuing questions:

- **Does the evidence support the case?** It supports removing redundant work, but the memory
  reduction here is modest and there is no uniform runtime improvement. For example, the longer
  basic directional case changes from 2.915 to 3.046 µs, while its extended counterpart changes
  from 5.972 to 5.244 µs. The before/after ranges overlap. These observations do not establish a
  general parser speedup, and fewer intermediate collections must not be counted as independent
  allocation savings when they reuse storage.
- **What can later subitems improve?** S3d4 can skip entire stereo derivations for irrelevant CX
  changes. S3d5 still targets two endpoint-vector allocations per considered double bond: forty
  such vectors in this 20-double-bond chain, compared with the two allocations removed here.
  S3e targets unnecessary neighbor construction on other paths, including raise, which this probe
  does not measure. These are concrete remaining costs, not demonstrated future speedups.
- **Should S3d3 be retained?** Yes. It removes a whole copied table, a separate participation
  allocation, and an intermediate collection pass. The replacement is limited to private marker
  bookkeeping and borrowed field access; it preserves the ordinary path's measured allocation
  footprint. The benefit is bounded, and the mixed timings do not justify further local tuning.

Verification: the focused parser suite passes 1,963 tests; the IO/graph gate passes 3,817 IO
and 1,066 graph unit tests with proptest enabled, all 10,223 SMILES / 2,253 MOL / 407 SDF
conformance cases, and the remaining integration/property suites. One existing graph doctest
remains ignored. Strict IO/graph Clippy passes for all targets with these features. Formatting
and diff checks pass. The complete incremental source/test diff was reviewed against S3d3;
pre-existing changes and later subitem scope were preserved.

Raw results are scratch/s3d3-{before,after}-{allocations,timings}.csv, with corresponding build
logs; the probe commands are the S3d2 commands with s3d3-parser-bench as the manifest directory.
Verification logs are scratch/s3d3-{focused-tests,gate,clippy,fmt}.log. The incremental review is
saved in scratch/s3d3-incremental.diff. S3d4 is next; no later subitem was implemented here.


### S3d4 implementation and measurements (2026-09-11)

The CTfile/CX derivation kernel now consumes the existing bond_stereo_assertions vector, validates
its sites, adds wavy evidence, sorts by bond index, checks conflicting repeats, and deduplicates
in place. An ordered iterator replaces the dense codes array. CX retains the earlier wavy
assertions when later entries overwrite the wedge field, so consolidation cannot discard that
evidence. Grouping may change which of several independent invalid assertions is diagnosed first;
the error variants and accepted source semantics are unchanged.

The kernel consumes and reuses the existing frame vector. It compares new codes and supplied
geometry against a retained frame after accounting for complementary reference choices at either
endpoint, without changing the retained references or relation. Existing Either suppresses geometry
and conflicts with definite codes; incoming Either conflicts with an existing definite frame.
New frames append to the reused vector, followed by an in-place table-order sort only when an
existing collection gained entries. Existing frames are unique and ordered by their private
parser producers. No frame-to-code expansion or rebuilt unchanged frame collection remains.

Basic and extended CX updates skip derivation when only unrelated metadata changes. Coordinates,
wedges, coordinate/hydrogen bonds, and multicenter entries mark stereo input as changed; explicit
bond stereo assertions also trigger derivation. The trigger is conservative for bond-related
entries. Labels, values, radicals, atom properties, groups, and other unchanged bond-frame inputs
do not cause an extra stereo pass. Relevant changes still check site order, references, missing
positions, and contradictory evidence at derivation. The redundant finish_stereo_bonds helper is
removed. CTfile builders transfer their owned assertion vectors into the same kernel.

No externally public Rust/Python type, constructor, conversion, signature, visibility, or error
variant changed. TableIR remains an open table carrier. The source interpretation uses supplied
coordinates only and performs no chemistry transformation. Local substituent vectors and neighbor
representation remain for their owning later subitems.

Exact tests cover all four reference choices, preserved vector storage when frames are unchanged,
Either/definite disagreement, invalid references, newly added frames before and after retained
frames, unsorted duplicates and conflicts, labels on existing directional frames, overwritten
wavy evidence, geometry agreement/conflict, short coordinate lists, and a CX hydrogen-bond order
change at an existing stereo site. Existing basic/extended parity, reaction-section remapping,
ring-completion indices, and CTfile property-order tests remain green. A new generated property
checks that any tested number of repeated equivalent codes preserves an existing frame across
both endpoint-complement choices; the coordinate-similarity property is unchanged.

The scratch/s3d4-parser-bench probe uses the established release measurement method: 64 warm-ups,
one allocation observation with result retention, and normal-allocator timing over 11 batches of
512 parses with output-vector allocation and result destruction outside timing. The eight inline
inputs each use basic and extended parsing: a bare 64-atom chain; that chain with a label; a
42-atom, 20-double-bond directional chain with a label; repeated cis codes; an existing directional
frame with agreeing trans code; terminal Either; a wavy bond; and supplied geometry. No fixture
files are read. Construction and verification stay outside measurement, and measurements run
separately from the test/lint gates. The baseline includes S3d3; affected source is retained under
scratch/s3d4-before. These are parser measurements, not raise/resolve timings or a corpus-wide
performance claim.

Allocation calls include reallocations; requested bytes sum full allocation/reallocation requests.
Peak added live bytes are allocator-tracked above the pre-parse baseline, not RSS or retained-result
size. Timing ranges are batch-mean ranges, not confidence intervals.

| Case | Allocation calls before → after | Requested bytes before → after | Peak added live bytes before → after | Median µs before → after | Before / after batch-mean ranges µs |
| --- | ---: | ---: | ---: | ---: | --- |
| bare_chain_64_basic | 3 → 3 | 12,188 → 12,188 | 9,428 → 9,428 | 1.003 → 0.990 | 0.977–1.686 / 0.971–1.770 |
| bare_chain_64_extended | 3 → 3 | 31,512 → 31,512 | 24,984 → 24,984 | 2.528 → 2.473 | 2.385–3.503 / 2.396–3.808 |
| labels_chain_64_basic | 74 → 8 | 21,600 → 17,953 | 14,677 → 14,677 | 2.731 → 1.262 | 2.678–3.232 / 1.209–1.469 |
| labels_chain_64_extended | 74 → 8 | 43,980 → 40,333 | 32,677 → 32,677 | 4.387 → 2.931 | 4.201–4.822 / 2.729–3.354 |
| labels_directional_42_basic | 184 → 95 | 29,710 → 25,557 | 18,005 → 18,005 | 5.257 → 3.159 | 5.180–6.517 / 3.077–3.904 |
| labels_directional_42_extended | 184 → 95 | 57,806 → 53,653 | 40,581 → 40,581 | 7.995 → 5.529 | 7.301–9.077 / 5.248–6.143 |
| repeated_code_basic | 18 → 17 | 7,579 → 7,576 | 6,432 → 6,432 | 0.850 → 0.829 | 0.817–0.898 / 0.809–0.891 |
| repeated_code_extended | 18 → 17 | 13,219 → 13,216 | 10,992 → 10,992 | 1.047 → 1.015 | 0.961–1.677 / 0.934–1.263 |
| existing_frame_basic | 24 → 22 | 6,739 → 6,672 | 5,600 → 5,600 | 0.863 → 0.818 | 0.852–0.940 / 0.787–1.026 |
| existing_frame_extended | 23 → 21 | 9,331 → 9,264 | 8,720 → 8,720 | 0.879 → 0.839 | 0.867–0.992 / 0.834–1.497 |
| either_basic | 11 → 10 | 5,457 → 5,456 | 5,248 → 5,248 | 0.498 → 0.487 | 0.494–0.522 / 0.483–0.748 |
| either_extended | 12 → 11 | 9,041 → 9,040 | 7,792 → 7,792 | 0.560 → 0.758 | 0.559–0.646 / 0.689–1.075 |
| wavy_basic | 15 → 14 | 6,415 → 6,412 | 5,532 → 5,532 | 0.654 → 0.643 | 0.636–0.691 / 0.622–0.792 |
| wavy_extended | 14 → 13 | 9,007 → 9,004 | 8,652 → 8,652 | 0.692 → 0.660 | 0.665–0.721 / 0.654–0.823 |
| geometry_basic | 19 → 18 | 9,771 → 9,768 | 7,960 → 7,960 | 0.839 → 0.994 | 0.828–0.911 / 0.865–1.089 |
| geometry_extended | 19 → 18 | 18,467 → 18,464 | 14,976 → 14,976 | 0.916 → 0.888 | 0.884–1.271 / 0.878–0.975 |

The first after run showed increases for extended Either and basic geometry. One unchanged-code
confirmation timing run after the gate checked repeatability; no allocation rerun or tuning was
performed. Both timing runs are retained rather than selecting only the faster result.

| Case | Before median µs | First after median µs | Confirmation median µs | Confirmation batch-mean range µs |
| --- | ---: | ---: | ---: | --- |
| bare_chain_64_basic | 1.003 | 0.990 | 1.001 | 0.975–1.777 |
| bare_chain_64_extended | 2.528 | 2.473 | 2.638 | 2.358–3.243 |
| labels_chain_64_basic | 2.731 | 1.262 | 1.340 | 1.208–2.058 |
| labels_chain_64_extended | 4.387 | 2.931 | 2.844 | 2.756–3.548 |
| labels_directional_42_basic | 5.257 | 3.159 | 3.162 | 3.080–3.282 |
| labels_directional_42_extended | 7.995 | 5.529 | 5.898 | 5.210–6.771 |
| repeated_code_basic | 0.850 | 0.829 | 0.838 | 0.814–0.865 |
| repeated_code_extended | 1.047 | 1.015 | 0.942 | 0.931–1.242 |
| existing_frame_basic | 0.863 | 0.818 | 0.796 | 0.791–0.817 |
| existing_frame_extended | 0.879 | 0.839 | 0.852 | 0.840–0.882 |
| either_basic | 0.498 | 0.487 | 0.490 | 0.487–0.648 |
| either_extended | 0.560 | 0.758 | 0.571 | 0.556–0.802 |
| wavy_basic | 0.654 | 0.643 | 0.650 | 0.627–0.677 |
| wavy_extended | 0.692 | 0.660 | 0.689 | 0.674–0.768 |
| geometry_basic | 0.839 | 0.994 | 0.837 | 0.826–0.859 |
| geometry_extended | 0.916 | 0.888 | 0.939 | 0.877–1.628 |

The two initial small-case increases did not repeat. The metadata-path improvements did repeat:
label-only chains lose 66 allocation calls per parse, and the labeled directional chains lose 89.
Requested bytes fall by 3,647 and 4,153 respectively. Peak added live bytes are unchanged in every
case: eliminating these later allocations does not lower the parser's earlier peak in these
inputs. Bare-chain allocation controls are unchanged. Relevant stereo cases remove one dense-array
allocation, or two allocations where unchanged frame rebuilding also disappears; their runtime
differences are much smaller and do not establish a broad speedup.

Assessment against the continuing questions:

- **Does the evidence support the case?** Yes, strongly for skipping unnecessary derivation. The
  metadata cases improve roughly 26–54% in median parse time across the two after runs, with large
  allocation-count reductions. This establishes the benefit of removing that whole
  pass on these inputs, not that allocator overhead alone explains its runtime or that all SMILES
  parsing becomes faster. Sparse storage and frame reuse alone show smaller benefits.
- **What remains promising?** S3d5 still removes the two endpoint-vector allocations per considered
  double bond, including forty in the 20-double-bond example's required first derivation. S3e
  addresses other unnecessary neighbor construction, including raise. Those costs remain concrete,
  but their future timing gains must be measured separately; they cannot be added arithmetically
  to this stage's gains.
- **Should S3d4 be retained?** Yes. It has a substantial measured benefit where CX metadata formerly
  repeated stereo work, reuses existing frame storage, and preserves the reference frame during
  evidence comparison. The extra merge logic serves the settled semantics, and the unchanged-code
  confirmation does not support pursuing the initial small-case timing increases further.

Verification passes: 3,844 IO unit tests and 1,066 graph unit tests with proptest; 10,223 SMILES,
2,253 MOL, and 407 SDF conformance cases; remaining integration/property suites; and strict
IO/graph Clippy for all targets with those features. One existing graph doctest remains ignored.
Formatting and diff checks pass. The complete incremental diff was reviewed against S3d4, with
pre-existing work and later allocation subitems preserved.

Raw results are scratch/s3d4-{before,after}-{allocations,timings}.csv and
scratch/s3d4-after-confirmation-timings.csv, with corresponding build logs. The probe uses the
same cargo run commands as S3d3 with s3d4-parser-bench as its manifest directory. Verification logs
are scratch/s3d4-{focused-tests,gate,clippy,fmt}.log; the reviewed diff is
scratch/s3d4-incremental.diff. S3d5 is next.


### S3d5 implementation and measurements (2026-09-11)

Both bond-stereo derivation kernels now gather at most two distinct substituents per endpoint
in inline SmallVec storage, using the existing dependency. Every push is guarded by the two-entry
limit; excess incidences never spill onto the heap. Duplicate incidences do not count toward that
limit. The retained entries are sorted to preserve the minimum table-index reference. SMILES
checks marker relevance before gathering endpoint ligands.

Raise uses a fixed array of two inline ligand blocks. It scans all qualifying incidences even after
finding excess ligands, so a reference appearing later still produces UnsupportedStereoBond and
an absent reference still produces InvalidStereoBondReference. Donation/noncovalent exclusions,
cumulated-axis rejection, shared-ligand rejection, and endpoint reference transport are unchanged.
SMILES still skips an empty-ended site before rejecting excess on the other endpoint; table
geometry still skips unsupported sites unless an explicit assertion requires a frame.

Public-symbol reconciliation: no public symbols, constructors, visibility, fields, or error
variants change. TableIR remains an open table carrier; raise remains the contextual consumer of
its explicit frames. There is no new normalization or chemistry inference. Persistent stereo-atom
ligand vectors, wedge interpretation, neighbor storage, and raise's constraint HashMap are unchanged.

The inline probe in scratch/s3d5-stereo-bench measures basic/extended parsing and borrowed TableIR
raise separately. Parsing, configuration, and table setup for raise occur outside its measurement;
raise does not clone its input. The probe checks complete output equality against its warmup result.
Timing uses the normal allocator, 64 warmups and eleven batches of 512 operations, retaining outputs
and checking/dropping them outside the timed region. Allocation measurements use a separate build
with the counting allocator. No fixture files are loaded. Before and after runs use the same probe;
no tests or other builds ran concurrently with timing. These are bounded examples and batch medians,
not confidence intervals or a workload-wide performance claim.

| Input | Parse allocation calls, basic and extended | Raise allocation calls | Basic parse median, ns | Extended parse median, ns | Raise median, ns |
| --- | --- | --- | --- | --- | --- |
| Bare chain, 64 atoms | 3 → 3 | 102 → 102 | 1053.6 → 1023.6 (-2.9%) | 2857.6 → 2946.5 (+3.1%) | 12630.5 → 12710.1 (+0.6%) |
| Directional chain, 20 double bonds | 90 → 50 | 160 → 100 | 2964.4 → 2502.7 (-15.6%) | 5080.6 → 4951.3 (-2.5%) | 12987.3 → 12303.4 (-5.3%) |
| Triene | 19 → 13 | 49 → 40 | 643.6 → 571.7 (-11.2%) | 754.1 → 707.0 (-6.2%) | 2441.3 → 2296.2 (-5.9%) |
| Branched conjugated sites | 19 → 14 | 45 → 39 | 636.7 → 566.8 (-11.0%) | 758.6 → 702.2 (-7.4%) | 2239.7 → 2131.8 (-4.8%) |
| Ring with atom and bond stereo | 17 → 15 | 58 → 55 | 588.6 → 557.9 (-5.2%) | 693.8 → 663.0 (-4.4%) | 2644.7 → 2599.3 (-1.7%) |
| CX coordinates | 18 → 16 | 38 → 35 | 953.4 → 802.1 (-15.9%) | 899.7 → 900.6 (+0.1%) | 1376.6 → 1270.3 (-7.7%) |
| Four substituents | 14 → 12 | 39 → 36 | 495.3 → 488.1 (-1.4%) | 576.5 → 582.0 (+1.0%) | 1724.9 → 1765.4 (+2.3%) |

The directional chain removes forty allocation calls and 640 requested bytes from either parse
path, and sixty calls and 1,600 requested bytes from raise. Its parser peak added live memory falls
by only 32 bytes; its raise peak is unchanged. Across the probe, parser peaks fall by 32 bytes on
lexically marked inputs, with unchanged peaks for the bare chain and CX geometry. All raise peaks
are unchanged. These savings remove small, short-lived allocations, not large simultaneously live
buffers. The branched parser loses five calls: it also avoids gathering the one nonempty endpoint
of an irrelevant terminal double bond. Bare-chain allocation counts and requested bytes are unchanged.

Assessment against the continuing questions:

- **Does the evidence support the case?** Yes. The expected per-site allocations disappear exactly.
  Basic parsing improves about 16% and raise about 5% on the directional chain; its extended parse
  improves about 2.5%. Triene and branched examples also improve. Small cases remain mixed, including
  a 2.3% raise increase for four substituents, so this does not establish a universal speedup.
- **What can later subitems improve?** S3d6 still targets the constraint HashMap and construction
  of constraints away from their destination. S3e targets the remaining neighbor allocation and
  unnecessary lookup construction; the directional parser still makes fifty allocation calls and
  raise one hundred. S3d7 addresses initial table reservations. Those costs remain, but this run
  does not predict their timing gains.
- **Should S3d5 be retained?** Yes. It removes all identified local ligand/block heap allocations,
  gives measurable gains on the multiple-stereo-site examples, and keeps the implementation local
  without new types or API layers. No further tuning is needed to close this subitem.

Verification passes: 3,880 IO and 1,066 graph unit tests with properties; 10,223 SMILES, 2,253 MOL,
and 407 SDF conformance cases; remaining integration/property suites; and strict all-target
IO/graph Clippy with those features. One existing graph doctest remains ignored. Final review
separated cumulation from excess ligands in a new raise fixture; its focused rerun also passes.
Formatting and diff checks pass. The full diff was reviewed against S3d5; no later subitem was
implemented. New exact cases cover zero, one,
two, duplicate, and excess incidences; descending reference order; first/second/both reference swaps;
shared ligands; cumulated axes; missing versus late references at overfull sites; and donation and
noncovalent exclusions. Existing stereo properties and conformance suites remain unchanged.

Raw allocation/timing results and build logs are scratch/s3d5-{before,after}-{allocations,timings}.*;
the complete measurement table is scratch/s3d5-comparison.md. Verification logs use
scratch/s3d5-{focused-tests,raise-tests,gate,clippy,fmt}.log. The reviewed diff is
scratch/s3d5-incremental.diff. S3d6 is next.


### S3d6 implementation and measurements (2026-09-11)

Raise now collects borrowed StereoBond references into one vector, sorts them in place by table-bond
index, and consumes them alongside bond construction. Each constraint is derived for its destination
and moved directly into that bond's constraints. The input frames remain untouched. The temporary
vector contains one reference per frame, not one entry per table bond; empty frame input allocates
nothing for this ordering. Sorting uses sort_unstable_by_key and adds no sorting buffer.

Frame derivation occurs before the localized/dative/noncovalent dispatch so a frame on a relation
bond cannot disappear silently, including the Shared donation fallback. Repeated sites produce
DuplicateStereoBond; a frame left after the bond table is exhausted produces
StereoBondIndexOutOfBounds. Existing invalid-reference and unsupported-site errors remain at frame
consumption. Atoms are now raised before bond frames, and frame errors are encountered in table-bond
order rather than input-frame order. Inputs containing multiple faults may therefore report a
different fault first; the eager all-frame validation/construction pass is not retained.

Public-symbol reconciliation: no constructor, conversion signature, field, visibility, error variant,
or Python surface changes. TableIR remains open, and raise checks frame context when it consumes the
frame. No new validation pass, dense bond-index mapping, or persistent adjacency is introduced.
Neighbor construction and the later allocation subitems remain unchanged.

The inline scratch/s3d6-raise-bench probe reuses S3d5's borrowed-raise measurement method: parsing,
configuration, and table preparation outside measurement; 64 warmups; eleven batches of 512 calls;
outputs retained and checked/dropped outside timing. Normal-allocator timing and counting-allocator
builds are separate. Reversed frame order is also measured on examples with multiple frames, with
reversal and cloning outside measurement. No fixture files are loaded and no tests/builds ran
concurrently with timing. The table reports batch medians, not confidence intervals.

| Case | Calls before → after | Requested bytes before → after | Peak live bytes before → after | Median ns before → after | Change |
| --- | --- | --- | --- | --- | --- |
| bare_chain_64_raise | 102 → 102 | 51828 → 51828 | 27320 → 27320 | 12827.5 → 11834.0 | -7.7% |
| directional_42_raise | 100 → 100 | 52348 → 51188 | 28312 → 27152 | 12891.2 → 11209.6 | -13.0% |
| directional_42_reversed_raise | 100 → 100 | 52348 → 51188 | 28312 → 27152 | 13136.3 → 11715.3 | -10.8% |
| triene_raise | 40 → 40 | 6968 → 6820 | 4984 → 4836 | 2432.2 → 2301.3 | -5.4% |
| triene_reversed_raise | 40 → 40 | 6968 → 6820 | 4984 → 4836 | 2514.7 → 2400.7 | -4.5% |
| branched_raise | 39 → 39 | 6840 → 6684 | 4856 → 4700 | 2175.9 → 2128.5 | -2.2% |
| branched_reversed_raise | 39 → 39 | 6840 → 6684 | 4856 → 4700 | 2337.7 → 2072.4 | -11.3% |
| ring_stereo_raise | 55 → 55 | 7584 → 7420 | 4988 → 4824 | 2745.6 → 2790.8 | +1.6% |
| geometry_raise | 35 → 35 | 3600 → 3436 | 3160 → 2996 | 1452.4 → 1235.6 | -14.9% |
| four_substituents_raise | 36 → 36 | 6512 → 6348 | 4544 → 4380 | 1871.1 → 1698.2 | -9.2% |

Allocation calls do not decrease: the frame-reference vector replaces the HashMap's one allocation.
For twenty frames, that replaces 1,320 requested bytes with 160, reducing both total requested bytes
and peak added live memory by 1,160 bytes. Three frames save 148 bytes, two save 156, and one saves
164. The bare-chain control's allocation metrics are unchanged. A constraint is no longer retained
in temporary map storage while waiting for all atoms and earlier bonds to be raised.

Assessment against the continuing questions:

- **Does the evidence support the case?** Yes for memory and storage ownership. The measured peak
  reduction matches the smaller temporary container. Timing is encouraging: about 13% lower for
  the ordered directional chain and 11% for reversed frames. However, the no-stereo control also
  improves about 8%, and ring stereo increases about 1.6%; these short runs do not isolate a
  stereo-specific speedup or support assigning all observed improvement to HashMap removal.
- **What remains for downstream improvement?** S3d7 still addresses initial SMILES reservations.
  S3e still addresses unnecessary neighbor construction and the nested neighbor allocations.
  The directional raise still makes one hundred allocation calls, and the bare-chain raise still
  makes 102; this stage does not remove those costs or establish their future timing benefit.
- **Should S3d6 be retained?** Yes. It directly constructs constraints where they are used, removes
  hash insertion/removal and temporary constraint storage, and gives a measurable peak-memory
  reduction with a small local implementation. Further timing experiments are not needed to close
  this subitem.

Verification passes: 3,889 IO and 1,066 graph unit tests with properties; 10,223 SMILES, 2,253 MOL,
and 407 SDF conformance cases; remaining integration/property suites; and strict IO/graph Clippy
for all targets with those features. One existing graph doctest remains ignored. Formatting and
diff checks pass, and the full diff was reviewed against S3d6. No later subitem was implemented.
Added exact tests cover full output and
input preservation for sorted/reversed frames with dative and noncovalent records before/between
framed localized bonds; nonadjacent duplicate frames; invalid and maximum-u32 sites; and rejection
of frames on donating, accepting, Shared, and noncovalent bonds. Existing reference transport and
stereo property tests are unchanged.

Raw results/build logs are scratch/s3d6-{before,after}-{allocations,timings}.*;
the full measurement table is scratch/s3d6-comparison.md. Verification logs use
scratch/s3d6-{focused-tests,gate,clippy,fmt}.log. The reviewed diff is
scratch/s3d6-incremental.diff. S3d7 is next.


### S3d7 implementation and measurements (2026-09-11)

Both SMILES inner parsers now cap their initial atom and bond reservations at 64 entries.
For short input, atom reservation remains the byte length and bond reservation remains its
saturating predecessor. Empty inner-parser input reserves zero. The existing builder receives
these bounded estimates; normal Vec growth accommodates arbitrarily larger tables. No counting
pass, grammar change, new configuration, or builder abstraction is introduced.

This bounds over-reservation caused by bracket syntax, CX text, and subsequent reaction sections.
It does not predict their exact entry counts or eliminate every surplus slot: trailing input can
still inflate a short section's reservation up to the cap. Geometric growth can also leave more
spare capacity than the former byte-length reservation on dense input just beyond a growth boundary.
Those are allocation tradeoffs, not acceptance limits. Public symbols, constructors, conversions,
visibility, diagnostics, source spans, and stereo semantics are unchanged.

A first measurement with a 32-entry cap reduced sparse-input reservations but introduced growth
for the 64-atom and 42-atom controls. One adjustment to 64 removes those growth reallocations while
retaining the major over-reservation reductions. No threshold sweep was performed; 64 is a bounded
initial policy, not a claim of an optimal or representative molecular size.

The inline scratch/s3d7-parser-bench probe uses the previous method: 64 warmups, eleven batches of
512 operations, results retained and checked/dropped outside the timed region, and separate normal-
and counting-allocator builds. Inputs and configuration are prepared outside measurement. It covers
basic/extended molecules and reactions, 64/256-atom chains, a 64-atom [13CH2] chain, twenty directional
double bonds, ring/atom stereo, a 1,280-character label on CC or C>>C, and reactions with a 256-atom
section preceding/following small or empty sections. No fixture files are loaded; tests/builds do
not run concurrently with timings. These are batch medians, not confidence intervals.

| Case | Calls before → after | Requested bytes before → after | Peak live bytes before → after | Median ns before → after | Change |
| --- | --- | --- | --- | --- | --- |
| empty_basic | 0 → 0 | 0 → 0 | 0 → 0 | 23.1 → 23.0 | -0.4% |
| empty_extended | 0 → 0 | 0 → 0 | 0 → 0 | 26.0 → 25.9 | -0.6% |
| small_basic | 4 → 4 | 1424 → 1424 | 1184 → 1184 | 244.5 → 250.2 | +2.3% |
| small_extended | 4 → 4 | 3408 → 3408 | 2832 → 2832 | 288.2 → 301.3 | +4.5% |
| bare_chain_64_basic | 3 → 3 | 12188 → 12188 | 9428 → 9428 | 995.0 → 1051.9 | +5.7% |
| bare_chain_64_extended | 3 → 3 | 31512 → 31512 | 24984 → 24984 | 2664.3 → 3335.5 | +25.2% |
| bare_chain_256_basic | 3 → 7 | 49044 → 77544 | 37844 → 37888 | 5450.2 → 7277.9 | +33.5% |
| bare_chain_256_extended | 3 → 7 | 126744 → 202208 | 100248 → 100352 | 14258.0 → 15182.5 | +6.5% |
| brackets_64_basic | 3 → 3 | 85900 → 12272 | 66260 → 9472 | 4257.6 → 1932.5 | -54.6% |
| brackets_64_extended | 3 → 3 | 221976 → 31712 | 175512 → 25088 | 6277.1 → 3968.3 | -36.8% |
| directional_42_basic | 50 → 50 | 19152 → 15584 | 15104 → 12336 | 2588.5 → 2670.3 | +3.2% |
| directional_42_extended | 50 → 50 | 44192 → 35024 | 35296 → 27952 | 4900.6 → 4436.2 | -9.5% |
| ring_stereo_basic | 15 → 15 | 3836 → 3836 | 3164 → 3164 | 594.4 → 622.2 | +4.7% |
| ring_stereo_extended | 15 → 15 | 8576 → 8576 | 7008 → 7008 | 708.7 → 714.0 | +0.7% |
| long_label_basic | 8 → 8 | 252112 → 17392 | 195464 → 14512 | 4150.0 → 2497.9 | -39.8% |
| long_label_extended | 8 → 8 | 643248 → 36832 | 509440 → 30112 | 4692.1 → 3518.5 | -25.0% |
| empty_reaction_basic | 4 → 3 | 396 → 292 | 352 → 252 | 238.6 → 261.7 | +9.7% |
| empty_reaction_extended | 4 → 3 | 1064 → 776 | 960 → 680 | 334.6 → 317.0 | -5.3% |
| long_product_basic | 6 → 10 | 98652 → 89816 | 76100 → 47344 | 7669.6 → 7117.9 | -7.2% |
| long_product_extended | 6 → 10 | 254952 → 233920 | 201624 → 125408 | 16525.6 → 17838.5 | +7.9% |
| long_agent_basic | 7 → 11 | 99320 → 89920 | 76616 → 47424 | 7335.7 → 7751.4 | +5.7% |
| long_agent_extended | 7 → 11 | 256704 → 234208 | 203040 → 125664 | 16163.3 → 17192.7 | +6.4% |
| long_reactant_basic | 4 → 8 | 49712 → 77648 | 38360 → 37968 | 5411.3 → 6770.9 | +25.1% |
| long_reactant_extended | 4 → 8 | 128496 → 202496 | 101664 → 100608 | 13903.7 → 14887.1 | +7.1% |
| reaction_label_basic | 13 → 13 | 501244 → 31616 | 386016 → 23968 | 4539.6 → 3630.7 | -20.0% |
| reaction_label_extended | 13 → 13 | 1283816 → 70496 | 1014208 → 55168 | 4940.8 → 5495.2 | +11.2% |

The long-label molecule requests 252,112 → 17,392 bytes in the basic parser and
643,248 → 36,832 in the extended parser (93–94% lower). Label-heavy reactions request
501,244 → 31,616 and 1,283,816 → 70,496 bytes, also about 94% lower. Peak added live memory falls
by a similar proportion. The bracket-heavy chain's requested bytes and peak fall about 86% with
no extra allocation calls. Uneven reactions with the large section later reduce peaks about 38%.

Conversely, the plain 256-atom chain gains four growth reallocations; requested bytes rise roughly
58–60% and its measured peak remains nearly unchanged. The large-reactant-first reaction has the
same growth cost with little peak reduction. These costs are real: bounding the reservation trades
away the previous exact estimate for dense single-character chains. Small controls and the plain
64-atom chain retain their allocation metrics. The empty reaction loses one allocation.

One unchanged-code timing confirmation was run after verification. The bracket-heavy chain
improves 55–57% basic and 37–45% extended across the two after runs; the long-label molecule
improves 40–48% basic and 25–37% extended. The 256-atom plain chain instead slows 17–34% basic
and 3–6% extended; the large-reactant-first reaction slows 18–25% basic and about 7% extended.
The label-heavy extended reaction remains slower by 11–13%, despite its large memory saving.
The unchanged-allocation 64-atom extended control varies from +25% to +14%, so timing shifts
cannot be assigned wholly to reallocations. These are observed ranges across two runs, not
confidence bounds. This closes the bounded measurement round without further tuning.

Assessment against the continuing questions:

- **Does the evidence support the case?** Yes for bounding surplus storage. The savings are large
  for long labels and bracket-heavy input, and uneven reactions with large later sections retain
  much less memory. Timing benefits are workload-dependent; larger plain chains gain reallocations
  and can become slower. The unchanged-code confirmation below retains both improvements and
  regressions; it does not support a universal speedup.
- **What remains for downstream improvement?** S3e still owns unnecessary/repeated neighbor
  construction and the nested neighbor allocations. This stage changes reservations, not those
  costs. Capacity tuning can be revisited separately; it is not required before S3e.
- **Should S3d7 be retained?** The bounded reservation is worth retaining for its large reductions
  in surplus storage, with the dense-input growth cost reported explicitly. It is not a universal
  speed optimization, and no further threshold tuning is part of this subitem.

Verification passes: 3,909 IO and 1,066 graph unit tests with properties; 10,223 SMILES, 2,253 MOL,
and 407 SDF conformance cases; remaining integration/property suites; and strict IO/graph Clippy
for all targets with those features. One existing graph doctest remains ignored. Formatting and
diff checks pass. The full diff was reviewed against S3d7, including the unchanged builder contract
and both parser call sites. No later subitem was implemented. Twenty new exact cases compare complete
basic and extended outputs, including spans and labels, for organic/bracketed chains of 1, 33, and
257 atoms and uneven reaction sections, with and without long labels. They verify growth without
pinning Vec capacities. Existing stereo properties, diagnostic tests, and conformance suites are
unchanged.

Raw results/build logs are scratch/s3d7-{before,after}-{allocations,timings}.*;
the first candidate uses scratch/s3d7-cap32-* and the final table is scratch/s3d7-comparison.md.
The unchanged-code confirmation is scratch/s3d7-after-confirmation-timings.csv with its build log;
its table is scratch/s3d7-confirmation-comparison.md. Verification logs use
scratch/s3d7-{focused-tests,gate,clippy,fmt}.log. The reviewed diff is scratch/s3d7-incremental.diff.
S3e is next.


### S3e implementation and measurements (2026-09-11)

AtomNeighbors now stores contiguous neighbor entries and row offsets privately. Its public
constructor, Neighbor entries, neighbors slice accessor, and distinct degree accessor are unchanged.
Construction collects compact endpoint pairs so the existing single-pass iterator input can be
replayed, counts incidences, and fills rows backwards in original bond-table order. The offset
array doubles as the fill cursor; no separate cursor allocation or adjacency sorting is needed.
For a nonempty borrowed table iterator, construction uses three allocations: temporary endpoint
pairs and the two retained vectors. No full Bond records are copied. Duplicate incidences,
single self-incidence, invalid-pair omission, and original bond-index gaps are preserved.

Raise uses one lazy operation-local lookup. Ordinary raising and explicit atom-frame raising
without wedges no longer construct it; Either bond frames also need none. Wedge interpretation
and framed bond transport share the lookup when needed. The wedge presence test is deliberately
coarse; this does not add a separate per-atom wedge index. Parser derivation similarly waits for
an incidence-dependent operation, avoiding lookup for absent evidence, reference-free Either,
and coordinates without relevant bond sites.

SMILES lexical normalization and subsequent CX updates share one lazy lookup, with independent
lookups for reactants, agents, and products. Reuse requires unchanged atom count, bond-table
indices, and endpoints. The current CX updates preserve those properties; they may change bond
attributes, which derivation reads from the current tables. CTfile normalization uses the same
lazy derivation path. These caches are local variables, not TableIR fields. No public symbols,
visibility, construction boundaries, or conversion guarantees changed. TableIR remains tables;
final Graph adjacency does not replace the all-table incidence lookup. No unused graph-core
adapter or new adjacency abstraction was introduced; traversal adaptation belongs at its consumer.

New tests compare contiguous lookup against an independent direct incidence scan, including
duplicates, self pairs, invalid endpoints, and bond-index gaps. Exact cases verify skipped lookup
for absent evidence and Either, and independent reaction sections with uneven indices. Fresh and
reused lookup paths produce identical complete results or errors for CX codes, compatible and
contradictory geometry, and bond-order updates. Existing raise, mixed-relation, malformed-index,
and stereo property laws retain their assertions.

The bounded scratch/s3e-neighbors-bench probe uses inline examples, 64 warmups, and eleven
512-operation samples. Timings use the normal allocator; a separate build counts allocations.
Raise input setup and output destruction stay outside timing, and results are checked against
the warmup output. Before and after timing runs had no concurrent tests or builds. The following
representative results are medians; requested bytes count cumulative allocation requests, while
peak bytes measure maximum live allocation during the operation.

| Operation | Allocation calls before → after | Requested bytes before → after | Peak bytes before → after | Time ns before → after |
| --- | --- | --- | --- | --- |
| Plain 64-atom MOL parse, basic | 69 → 4 | 14,368 → 10,784 | 14,368 → 10,784 | 10,844.8 → 9,482.7 |
| Plain 64-atom MOL parse, extended | 69 → 4 | 29,672 → 26,088 | 29,672 → 26,088 | 12,664.1 → 11,140.1 |
| Plain 64-atom raise | 102 → 37 | 51,828 → 48,244 | 27,320 → 23,736 | 11,774.3 → 10,321.3 |
| Directional 42-atom parse, basic | 50 → 10 | 15,584 → 14,560 | 12,336 → 10,984 | 2,491.1 → 1,861.8 |
| Directional 42-atom parse, extended | 50 → 10 | 35,024 → 34,000 | 27,952 → 26,600 | 3,921.2 → 3,251.5 |
| Directional 42-atom raise | 100 → 60 | 51,188 → 50,164 | 27,152 → 25,800 | 11,506.3 → 10,489.6 |
| Atom-frame raise | 47 → 42 | 4,216 → 3,992 | 3,224 → 3,000 | 1,648.2 → 1,556.8 |
| Either-frame raise | 31 → 26 | 3,412 → 3,188 | 2,996 → 2,772 | 1,031.3 → 860.4 |
| Direction plus CX code parse, basic | 18 → 11 | 6,608 → 6,272 | 5,600 → 5,688 | 754.6 → 632.2 |
| Direction plus CX code parse, extended | 17 → 10 | 9,200 → 8,864 | 8,720 → 8,808 | 1,075.6 → 851.0 |
| Branched stereo raise | 39 → 33 | 6,684 → 6,476 | 4,700 → 4,436 | 2,010.8 → 2,034.7 |
| Plain 64-atom SMILES parse, extended control | 3 → 3 | 31,512 → 31,512 | 24,984 → 24,984 | 2,697.2 → 2,924.3 |

Assessment against the continuing questions:

- **Does the evidence support the case?** Yes. Skipping the plain 64-atom lookup removes 65
  allocation calls; contiguous construction removes 40 calls for the 42-atom directional chain.
  MOL parsing and ordinary raise improve about 12%; directional parsing improves 17–25%.
  Sharing lexical/CX lookup saves seven calls in the small example, but retaining it through CX
  raises peak live storage by 88 bytes while reducing total requested storage. Timings are not
  universally better: branched stereo raise is 1.2% slower and the unchanged-allocation extended
  plain-chain control is 8.4% slower. This is bounded before/after evidence, not a universal speed claim.
- **What remains for downstream improvement?** Required lookups still replay compact endpoint
  pairs. Final table/Graph construction, wedge-specific scratch collections, and the previously
  recorded capacity tradeoffs remain. These results establish no need for another allocation
  subitem or further tuning; S4a moves to ordinary valence reconstruction.
- **Should S3e be retained?** Yes. It removes unused whole lookup construction and per-atom
  allocations, with substantial measured reductions and unchanged public semantics. The small
  CX peak-lifetime tradeoff is explicit and does not outweigh those savings.

Verification passes: 3,920 IO and 1,066 graph unit tests with properties; 10,223 SMILES, 2,253 MOL,
and 407 SDF conformance cases; remaining integration/property suites; and strict IO/graph Clippy
for all targets with those features. One existing graph doctest remains ignored. Formatting and
diff checks pass. Review covered the complete diff, including lookup producers, cache lifetimes,
CX mutation boundaries, public contract preservation, and unchanged property assertions.

Raw results are scratch/s3e-{before,after}-{allocations,timings}.csv with corresponding build logs;
the full measurement table is scratch/s3e-comparison.md. Verification logs are
scratch/s3e-{focused-tests,gate,clippy,fmt-check}.log. The reviewed diff is
scratch/s3e-incremental.diff. S3 is complete; S4a is next.

## Staged specification updates

This is the separate staging list for changes to `umol-io/spec/opensmiles-spec.md`.
It records specification content, not implementation tasks or roundtrip guarantees.
The specification defines interpretation and equivalence; the roundtrip design consumes
that equivalence relation. The isotope correction below is applied; the remaining changes are staged.

### Isotope composition

- Applied in S4a0c (2026-09-11) to the Numeric and Bracket Fields and Organic Subset sections:
  an omitted isotope specification denotes naturally occurring isotopic composition, not an
  undetermined isotope composition. Natural composition does not select one mass number. This follows OpenSMILES
  section 3.1.4.

### Implicit hydrogen section

- State explicitly that an implicit hydrogen, including a bracket H-count participant,
  is not equivalent to an explicit hydrogen atom connected by a bond. Preserve this
  representational distinction even when both descriptions denote the same chemistry.
- Retain the existing bracket rule: omitted H means H0, and bracket atoms undergo no
  implicit hydrogen calculation. A stereo marker does not change that rule.

### Stereo-atom interpretation

- A bracket hydrogen count contributes the corresponding implicit-hydrogen participants
  to the tetrahedral frame. An explicit `[H]` neighbor is an actual atom participant.
- With exactly three actual incident ligands and bracket H0, the fourth participant
  of the supported tetrahedral frame is a lone pair. This is the participant asserted
  by the descriptor, not a determination that the atom has a chemically valid lone pair.
  It does not set molecular lone-pair or unpaired-electron counts. Molecular interpretation
  checks whether the atom's resolved state supports the asserted frame.
- Specify participant ordering: insert bracket-H participants, or the lone-pair participant
  in the three-neighbor H0 case, before the first actual ligand at a traversal root and
  immediately after the incoming ligand otherwise. Preserve the written incidence order,
  including ring-closure encounters. Identify the lone-pair placement as the umol convention;
  do not claim that OpenSMILES fully specifies it.
- Add no other padding to an incomplete frame. A descriptor does not authorize inventing
  an atom or hydrogen, or distinct occurrences of an identical virtual participant.
  Unsupported frame arity or repeated participants are not repaired by normalization.
- These are interpretation clarifications. They introduce no additional stereo-atom
  equivalences and no equivalence between implicit and explicit hydrogen representations.

### Stereo-bond equivalence and redundant notation

- Valid partial direction markers that establish no double-bond configuration are redundant
  after all markers have been checked for consistency. For example, `F/C=CF`, `F\C=CF`,
  and `FC=CF` are equivalent with respect to double-bond configuration.
- Additional consistent markers that restate an established configuration are redundant.
  For example, `F/C(Cl)=C(/Br)I` and `F/C(/Cl)=C(/Br)\I` encode the same configuration.
- Reversing all markers in a connected directional assignment preserves its configurations;
  for example, `F/C=C/F` and `F\C=C\F` are equivalent.
- Selecting different substituent bonds or ring-closure endpoints to encode the same
  configuration is equivalent when the direction changes account for endpoint viewpoint
  and reference changes. This does not permit creating or removing explicit hydrogen atoms.
- Assess equivalence over the complete directional system. It preserves exactly which
  double bonds have definite configurations and the configurations asserted at those sites.
  Shared markers in conjugated systems must not be discarded or generated independently
  when doing so changes another site's assertion.
- Preserve an explicit Either assertion as distinct from absence wherever the boundary
  supports that assertion. Partial slash notation is not an explicit Either assertion.
  Correct the current Double-Bond Stereochemistry paragraph that conflates incomplete
  direction information with unknown/either, including its ring-closure wording.
- Conflicting markers remain errors, even if another endpoint is unmarked. Equivalence
  and normalization do not erase conflicts. The Direction-marker examples above supply
  concrete equivalence and rejection cases; they are evidence for this separate list.
