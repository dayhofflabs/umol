# 226 — SMILES roundtrip design

Status: In Progress
Date: 2026-09-10
Relates: [155](155-smiles-io-and-resolve-configuration-2026-07-19.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[166](166-molecule-ops-2026-07-27.md),
[170](170-reaction-smiles-python-2026-07-28.md),
[224](224-smiles-ring-closure-frame-2026-09-08.md),
[data type contracts](../docs/development/data-types.md)

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
localized Kekulé representation introduces a separate selection question. We must not assume that
an arbitrary localized assignment preserves every graph-IR aromatic electron-count, charge, or
spin distinction. The projection contract must identify which states the notation can reconstruct
under the selected model and reject unsupported distinctions rather than silently erase them.

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
2. Convey/render/parse/interpret preserves the supported graph-IR molecular meaning under the
   same chemistry model and interpretation settings, allowing entity renumbering and stereo-frame
   transport. This requires a declared representable domain, including retained constraints,
   electron accounting, hydrogen conventions, and stereo state. Being concrete alone is not a
   sufficient representability condition. Boundary metadata absent from graph IR is outside this law.

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

## Resolver projection: inverse target and initial scope

The desired law, under the same resolver, is:

```text
resolve(project(resolved_molecule)) ≡ resolved_molecule
```

Equivalence accounts for representation ordering and stereo-frame transport. This is recovery
of the resolved state, not recovery of the original unresolved input or its discarded assertions.
The full roundtrip is the target, but initial support may cover a small, explicit feature domain
and expand incrementally. Project must return an error where faithful reconstruction cannot be
established; it must not guess a different state or silently discard information.

Finding a minimum set of retained properties that selects the original valence state is not an
initial requirement. Retain enough information to recover that state. For atom typing, investigate
narrowing the candidate list against retained properties, potentially through a registry pruning
operation; no registry method or signature is settled yet. Counts-based valence reconstruction
requires separate analysis of the reverse algorithm and need not have a general solution in the
initial scope.

For aromaticity, the starting approach is to set atom and bond `#a` assertions in projected graph
IR. Establish the domain in which those assertions, together with retained properties, reconstruct
the original resolved systems and atom states under the same resolver. This is a proposed starting
encoding, not a claim that aromatic flags suffice for every system. If the required properties
cannot be set unambiguously for faithful reconstruction, fail rather than approximate.

Use property tests of the inverse law over each supported domain, alongside the full
convey/render/parse/interpret roundtrip properties. Broaden support as reconstruction rules become
established; neither general inversion nor minimal encodings block an initially useful subset.

### Ordinary atom valence

Start reconstruction analysis with non-aromatic atoms, localized bonds, and explicit bracket H
counts. Preserve element, isotope, charge, H count, and bond orders. Bracket H counts remain
implicit hydrogen participants, not explicit hydrogen atoms. Avoid arbitrary element restrictions:
the initial domain should be useful for ordinary molecules, charges, and reconstructible radicals.
Compact unbracketed H-inference spelling can follow once its reconstruction rule is established.

TableIR can retain lone-pair, unpaired-electron, and multiplicity fields that ordinary SMILES
cannot explicitly carry. Their presence in TableIR alone is not a roundtrip proof. Candidate
reconstruction must use properties that actually survive the selected notation.

For non-aromatic counts resolution with fixed H, the current algorithm derives nonbonding
electrons from element valence electrons minus charge, bond valence, and H count. With lone pairs
and unpaired electrons unspecified, it pairs the remainder except for its parity. Compare the
resulting state, including multiplicity, against the source; differing pairing or spin fails
unless another supported encoding preserves it. This handles a useful initial domain without
solving the general counts inverse.

Atom typing already filters registry rows by element, charge, retained properties, and
asserted/derived incidence constraints. A conservative initial criterion is one completed state
matching the source; multiple rows yielding the same state need not imply ambiguity. Evaluate
the proposed retained properties, not the fully ground source atom, which admission may skip.
These are initial reconstruction criteria to verify, not completed roundtrip coverage.

### Aromatic reconstruction and the inverse check

Projection fails for nonzero charge or `#u > 0` on any bond or aromatic system; it does not
discard or implicitly localize them. This projection restriction complements the general export
representability rule above.

Use the existing joint aromaticity selection machinery, not a second selection algorithm.
Individual atoms may have multiple candidates while the system has one selected assignment.
Removing stored aromatic systems also removes the selector's requirement to reproduce their
member sets; `#a` assertions and retained properties must reconstruct the original grouping,
atom states, and electron contributions rather than merely some accepted aromatic assignment.

For the initial implementation, construct a proposed projected GraphIR value, resolve a copy
under the same resolver, and require determined success and equivalence to the source. Changed
assignments, changed grouping, unresolved ambiguity, or an existing search limit cause projection
failure. This reconstruction check establishes the inverse guarantee rather than defensively
revalidating representation integrity. Convey must separately establish that the required
information survives the boundary notation; a successful GraphIR check does not prove that.

### Stereo projection

For ordinary tetrahedral and cis/trans stereo, reuse the assertion frames constructed by
StereoPerception::derive_stereo_atom and derive_stereo_bond in
[stereo.rs](../umol-graph/src/ops/stereo.rs):

- The `#T` frame has actual neighbors in GraphIR neighbor order, followed by an implicit-H or
  lone-pair participant when needed to complete four ligands.
- The `#C` frame follows the bond's ordered endpoints, each with its actual substituent neighbors
  excluding the opposite endpoint and a virtual participant when needed to complete that side.

Obtain the frame the same stereo model would reconstruct, use coset_for to transport the stored
entity configuration into it, write the corresponding `#T` or `#C` assertion, and remove the
entity. The resolve-and-compare check must recover the original stereo. Failure to construct a
frame, incompatible ligand identities, or unrepresentable configuration causes projection failure.
Do not relabel lone pairs as hydrogen, fold or expand hydrogen atoms, or silently drop stereo
under permissive failure policies. Configured scope and perception exclusions are also subject
to the reconstruction check. TableIR frame construction remains entirely with Convey.

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
  depth is shortest distance from that traversal root.
- Emit Finish after examining a node's neighbors. An optional maximum depth suppresses
  expansion at the limit, but nodes there still receive Discover and Finish.
- Preserve supplied neighbor order and stop immediately on Break.

Graph::neighborhood supplies one root and its depth limit, collecting `(node, depth)` from
Discover. Connected-component enumeration supplies all nodes without a depth limit, starts
a component at each parentless Discover, collects its nodes, and sorts each completed component.
Sorting is an output operation, not a change to BFS queue order. Component order follows the
first unreached candidate roots. With a depth limit, later roots can reach previously unvisited
parts of the same connected component; the component interpretation requires unrestricted depth.

### The current LIFO flood fill

Despite ConnectedComponentsAlgorithm::Bfs and its documentation, the current implementation in
[components.rs](../umol-graph-core/src/algorithms/connectivity/components.rs) uses Vec::pop and
marks all unseen neighbors when pushing them. It is a LIFO reachability flood fill. It finds
components correctly, and sorting their members hides visitation order in the returned result,
but it supplies neither BFS distances nor the proposed nested DFS event semantics.

For edges A–B, A–C, and B–C, that flood fill discovers B and C from A before exploring either.
Recursive-style DFS instead explores one and discovers the other through it. A stack alone does
not establish the DFS frame contract.

The agreed direction is to replace this component flood-fill loop with collection over the shared
BFS implementation, retaining sorted component results and the explicit Bfs selection. No separate
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
introduce new wrappers, strategies, or chemistry transformations. Implement the supported inverse
domain with explicit failures, not a promise of general valence inversion or minimum encodings.
The plan below sequences the work; S0 is complete.

## Staged implementation plan

S0a–S0c and S1a are complete; S1b and later subitems are pending.
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
| umol-graph resolver | GraphIR-only project, reconstructible atom-state reductions, aromatic/stereo assertion recovery, atomic inverse verification |
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
  for the GraphIR law; compare atom states, system members/contributions, and stereo explicitly
  rather than formula or entity counts. Add named successful ordinary, charged, radical, aromatic,
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
- **S1b — BFS events, callback kernel, and Graph visitor.** Same module.
  **Additive (green).** [dep: S0a, S0b]
  Add BreadthFirstEvent and visit_breadth_first with FIFO discovery, sequential candidate roots,
  usize depth/limit, per-node Finish, and FinishTree before the next root. Verify shortest distances
  independently, maximum-depth boundaries, duplicate roots, disconnected components, exact order,
  and Break with no further callbacks. Document that depth-limited multi-root trees are not
  necessarily full components and are not multi-source BFS.
- **S1c — Event collectors and component visitor.** Modules: traversal.rs and
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

- **S2a — Neighborhood naming and usize contract.** Modules: graph-core traversal, refinement,
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

- **S3a — Stereo-bond vocabulary and frame algebra.** Module: umol-io/src/table_ir/stereo.rs
  and its exports. **Additive (green).** [dep: S0a, S0b]
  Add StereoBond { bond, configuration }, BondConfiguration::{Either, Framed { references,
  relation }}, and BondRelation::{SameSide, OppositeSide}. Pin endpoint/reference ordering and
  correspondence with existing permutation/coset operations. Verify one-reference swap, both swaps,
  endpoint exchange, and transport through atom/bond remapping with independent expected frames.
  Preserve StereoAtom, Winding, and LonePair; no generic placeholder or extra Either relation.
- **S3b — Source direction normalization.** Modules: smiles parser builder and existing annotation
  processing. **Additive (green).** [dep: S3a]
  Implement the cohesive derivation functions with direct unit tests before switching producers.
  Consume all source markers in endpoint viewpoint, including ring-opening/closing spelling and
  shared bonds. Select lowest-index actual references, check same-side redundancy and conflicts,
  and produce no assertion for consistent partial notation. Keep lexical directions in parser
  working state; do not create another persistent authority. Verify all recorded equivalence and
  rejection cases, including conflicts when the other endpoint is unmarked and the doc 224 order.
- **S3c — CTfile and annotation frame derivation.** Modules: CTfile readers and CX bond annotation
  producers. **Additive (green).** [dep: S3a]
  Reuse supported geometry interpretation to derive definite frames, and preserve site-only Either
  without inventing reference atoms. Treat adjacent wavy-bond effects under the existing convention;
  wedges remain separate. Verify usable 2D/3D, absent/all-zero/degenerate geometry, explicit Either,
  and supported Cis/Trans annotations with exact frames. No coordinate generation or CTfile output.
- **S3d — Publish frames and retire duplicate fields.** Modules: TableIR Molecule/ExtendedMolecule,
  Bond/ExtendedBond, parser targets/conversions, raise, all affected consumers/tests.
  **Breaking (red→green).** [dep: S3a, S3b, S3c]
  Add the stereo_bonds collection consistently to applicable carriers, wire every producer, and
  remove persistent Bond.direction and Bond.stereo plus superseded equivalents. Retain lexical
  marker types only where they still serve parsing. Raise Framed through the settled #C frame
  transport and Either through an undetermined configuration assertion. Preserve out-of-range,
  incidence, and conflicting-annotation failures at their owning boundary. Move source interpretation
  out of raise once frames own its result. Update struct literals, fixtures, conversions, and
  benchmarks together; no consumer may silently prefer one of two operative stereo encodings.
- **S3e — Operation-local neighbor lookup.** Modules: TableIR raise and parser normalization.
  **Additive/refactor (green).** [dep: S3d]
  Remove unconditional adjacency construction where frame-based raise needs none; derive temporary
  lookup only for remaining operations that require it. Preserve all-table-bond indices and the
  partition into localized edges versus relations. TableIR remains tables. Use iterator adaptation
  to graph-core Neighbor for traversal rather than adding Graph storage or an adjacency trait.
  Verify ordinary/stereo/mixed-bond raising and malformed-index behavior; use the bounded raise
  fixtures to report allocation changes without opening a benchmark-driven redesign.

**Gate:** IO unit/property/CTfile suites and dependent graph interpretation tests pass. Every former
raw-direction/stereo-field producer and consumer is accounted for. Parser acceptance and error-layer
changes are reviewed for semantic preservation, not accepted solely because tests were rewritten.

### S4 — GraphIR-only inverse projection

- **S4a — Ordinary valence reconstruction.** Modules: resolve/valence.rs and valence
  atom_typing/counts/registry as needed. **Additive (green).** [dep: S0a, S0b]
  Implement reconstruction from retained GraphIR properties for the initial non-aromatic domain.
  Reuse candidate admission and compare completed states, deduplicating equivalent completions;
  do not check the already-ground source through the admission fast path. For counts, start with
  fixed H and the established electron pairing/multiplicity reconstruction. Retain enough evidence
  rather than search for a minimum encoding. No registry public pruning API is added without its
  contract review; reuse existing capabilities where sufficient. Cover useful ordinary atoms,
  charges and radicals across both candidate sources, model/tie-break differences, and exact
  rejection of unreconstructible pairing/spin. This code knows no TableIR or SmilesIoConfig.
- **S4b — Aromatic assertion recovery.** Module: resolve/aromaticity.rs.
  **Additive (green).** [dep: S4a]
  Produce atom/bond #a assertions while removing the corresponding system entities and retaining
  reconstruction evidence. Use existing joint selection to recover atom contributions and system
  grouping. Test benzene, heteroaromatics, localized links between aromatic systems, fused systems
  within the supported domain, atom-localized charge, and differing admissible partitions. Reject
  nonzero system charge or #u > 0 and undetermined unsupported attributes; do not localize them.
- **S4c — Stereo assertion recovery.** Module: resolve/stereo.rs, reusing ops/stereo.rs and
  graph-IR coset transport. **Additive (green).** [dep: S4a]
  Recover #T/#C assertions in the frames that the same model derives, remove entities after
  transporting their configurations, and preserve typed virtual participants. Test reordered
  actual ligands, implicit versus explicit H, LP frames, endpoint flips, Either where in domain,
  model scope exclusions, and permissive resolver policies that would otherwise lose stereo.
  Unsupported extended kinds fail rather than partially projecting.
- **S4d — Atomic Resolver::project and inverse verification.** Module: resolve.rs and related
  error/report vocabulary. **Additive (green).** [dep: S4a, S4b, S4c]
  Compose the projection plans with resolve's Result/Solution and commit semantics. Reject nonzero
  charge or #u > 0 on bonds and aromatic systems. Preserve non-elidable information on other
  unsupported structures for explicit failure, never erase it to make the inverse pass. Resolve
  a copy of the candidate with the same resolver and compare the recovered molecular meaning,
  including aromatic partitions and transported stereo. Only determined matching success publishes
  the projected value; all other outcomes leave the caller unchanged. Test the inverse law and
  atomic failure independently, including candidate ambiguity and configured search limits.

**Gate:** graph unit, conformance, and property suites pass for both valence strategies. Record the
supported domain and exact unsupported categories. Projection benchmarks include candidate building
and the deliberate re-resolution cost; do not replace the latter with an unproven shortcut.

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
  Build StereoAtom/StereoBond frames here, not in the resolver. Check narrowing ranges and every
  nondefault property/constraint that the boundary cannot preserve; unsupported information errors
  rather than disappearing. Use the checked boundary constructor. Verify the graph inverse and
  the stronger conveyed/rendered/parsed/interpreted law separately so retained internal electron
  fields cannot hide loss in text. Tests include original-input immutability on success/failure.
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
  Run the separate laws for GraphIR inverse projection, boundary normalization, same-config
  construction/rendering, and graph export/ingest. Include both valence strategies, supported
  model policies, source-order variations, explicit/implicit H distinctions, aromatic grouping,
  marker-selection impossibility and alternatives, reactions and correspondence. Add retained
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

Deferred and not prerequisites for the Rust roundtrip: minimum-information projection, general
counts inversion beyond established support, early parity pruning, minimum-marker optimization,
canonical output, cached traversal/markers, new Python boundary/output surface, CTfile writing,
extended stereo, coordinate generation, and depiction wedge policy. Expanding the chemistry domain
requires demonstrated reconstruction rules, not automatic acceptance. No separate optimization
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
| Resolver::project and per-resolver project functions | New GraphIR-only inverse transformation, same mutation/publication semantics as resolve; no TableIR or IO-config parameters. | Convey; inverse-law properties and projection benchmarks. Concrete report/contradiction types are not yet named by the design. |
| Resolver::resolve, ResolveConfig, ResolveState and existing reports/errors | Retain existing semantics. Project's validation uses the same resolver; do not change resolution or its chemistry defaults to force roundtrips. | Existing ingestion and resolution tests; new projection check. |
| AtomTypeRegistry lookup/admission operations | Reuse candidate filtering initially; no approved new registry pruning signature. | Inverse valence analysis. Any new invariant-bearing public operation needs its own contract within S4a. |
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
  V2000/V3000, accumulator overrides, and basic/extended conversion must all reach the same result.
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
| Resolver::project | Caller supplies GraphIR; work is GraphIR throughout, on an intermediate candidate. | Result/Solution follows resolve; only matching determined reconstruction publishes. Failure preserves caller input. |
| Convey and text export | Convey clones, projects, translates, and calls checked construction. Export composes convey/render. | Preserve source graph and supported semantics; reject unrepresentable fields/constraints. A GraphIR inverse check alone does not prove text preservation. |
| Reaction boundary construction | A supplied table may have independently assembled labels/index; the constructor owns the required consistency check. | Preserve label semantics including parsed repeats/agents where representable; GraphIR convey produces its index from correspondence labels. |

The laws to verify in later subitems are traversal visitor/collector equivalence and completion,
projection inverse plus atomic failure, boundary normalization and same-config renderability,
and molecular/reaction export-ingest equivalence with explicit/implicit H and stereo-frame action.
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
- **GraphIR inverse law:** use Molecule::framed_eq when entity IDs are retained,
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

## Staged specification updates

This is the separate staging list for changes to `umol-io/spec/opensmiles-spec.md`.
It records specification content, not implementation tasks or roundtrip guarantees.
The specification defines interpretation and equivalence; the roundtrip design consumes
that equivalence relation. The specification itself has not yet been edited.

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
