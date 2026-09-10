# 226 — SMILES roundtrip design

Status: Proposed
Date: 2026-09-10
Relates: [155](155-smiles-io-and-resolve-configuration-2026-07-19.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[170](170-reaction-smiles-python-2026-07-28.md),
[224](224-smiles-ring-closure-frame-2026-09-08.md),
[data type contracts](../docs/development/data-types.md)

## Scope

Design the reverse molecular-format path: graph-IR Molecule through resolver-owned projection
to TableIR, then formatting through Smiles or ReactionSmiles. TableIR is the shared SMILES/CTfile
boundary representation. Decide what it must carry for aromatic systems and stereo without
splitting it by output format. This document records the initial design discussion; the choices
below are proposals, not an approved API or implementation plan.

The immediate outcome is agreement on responsibilities, preservation guarantees, and the boundary
types. Sequencing, algorithms, exact signatures, and implementation follow only after that agreement.
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
- Resolver and the aromaticity/stereo resolvers currently have no TableIR project operation.
  Resolution is a transformation; the proposed reverse projection must have a separate contract:
  see [resolve.rs](../umol-graph/src/ops/resolve.rs). umol-graph already depends on umol-io, so
  resolver-owned TableIR production fits the existing dependency direction.
- Graph ingestion composes boundary interpretation with resolution. ReactionSmiles interpretation
  rejects nonempty agents and repeated members of a map class on either side before resolving the
  sides: see [ingest.rs](../umol-graph/src/ingest.rs). Boundary parsing retains those representations.

## Two operations and their boundaries

```text
graph-IR Molecule
    -- resolver projection --> TableIR Molecule
    -- SMILES boundary construction --> Smiles
    -- rendering under IO configuration --> text

text -- parse --> Smiles -- TableIR raise + explicit resolution --> graph-IR Molecule
```

Boundary construction establishes that the projected table can be expressed as SMILES. It is a
contract boundary within the second operation, not a proposed extra processing framework.

Projection reads the supplied molecule without modifying it. Its umol-specific responsibility is
to express the stored atom/bond information and supported overlays in TableIR terms, with explicit
frame transport. Aromaticity and stereo work should remain with their owning resolver machinery.
Projection does not rerun resolution implicitly, select a surviving interpretation, or use recovery
policies to discard inconvenient information. The exact model dependence must be stated: a model
may be needed to establish that notation reinterprets to the supplied state, rather than merely to
copy fields.

Formatting owns traversal, components, branches, ring labels, bracket/implicit-H spelling, and
stereo markers relative to the emitted traversal. The writer consumes boundary information without
aromaticity perception or chemical resolution. SmilesIoConfig retains doc 155's bidirectional role;
chemistry model choices remain separate. Algorithm selectors, if exposed at the IO layer, are explicit.

## Shared TableIR and apparently duplicated information

Two distinct questions must not be conflated: whether TableIR can carry both formats' input
notations, and whether projection must populate both output notations simultaneously.

### Aromatic systems

Atom.aromatic and BondOrder::Aromatic already express different incidences. Aromatic endpoint atoms
do not by themselves determine whether their connecting bond is aromatic. The projected table must
distinguish an aromatic-system bond from a localized connection between aromatic atoms. The writer
then chooses lowercase symbols and explicit or omitted bond tokens from those values.

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
copy the old chirality token after changing traversal. A CTfile writer derives its own supported
stereo notation in relation to its output coordinates and numbering; projection need not invent
wedges or geometry in advance.

Stereo bonds expose the remaining representation choice. Bond.direction describes a single-bond
direction; BondStereo alone does not supply an explicit pair of endpoint ligand frames. Possible
designs are:

| Choice | Consequence |
| --- | --- |
| Populate format-specific fields for the intended output | Keeps the current field family, but makes projection output-dependent and requires the SMILES writer to preserve or transport the reference convention. |
| Add explicit stereo-bond frame information alongside StereoAtom | Gives both writers a common configuration to translate; requires coordinated TableIR construction and raise semantics. |
| Populate both format-specific encodings | Requires a defined agreement rule, and CTfile geometric encodings need output geometry before they can be produced. |

The second is the initial preference, not a settled type change. It does not require splitting
TableIR or copying the full graph-IR overlay model into it. Generalizing StereoAtom beyond
tetrahedral configurations is another explicit scope decision; parsed chirality variants alone do
not establish end-to-end support for every stereo kind.

If explicit frames and source-format annotations coexist, their roles must be stated. Retained
source notation can be evidence without being a second independent configuration. Independently
supplied operative encodings must agree; a consumer must not silently choose between contradictory
configurations. The current tetrahedral raise precedence is implemented behavior, not proof that a
future general multi-encoding contract is already enforced.

## Smiles and ReactionSmiles as output boundaries

Keep their private payloads. Broaden their role from parsed values to SMILES-representable boundary
values, constructible by parsing or checked conversion from TableIR. Public arbitrary-parts
construction is not implied. Exact conversion names and errors remain undecided.

TableIR is an open carrier: the first public conversion that promises SMILES output must check
references, usable stereo frames, mutually consistent operative fields, and representability under
the requested IO configuration. A trusted projection producer can establish some of these by
construction. Smiles construction does not certify chemical validity under a chemistry model.
Rendering with a different, narrower syntax configuration may still fail; construction under one
configuration cannot promise output under every other configuration. No stored configuration or
additional witness type is proposed merely to avoid this distinction.

ReactionSmiles should compose the molecular writer over the three ordered sections and preserve
map-class information. Its boundary roundtrip must retain agents and repeated classes even where
conversion to one graph-IR Reaction is unavailable. Atom.class and Reaction.atom_mapping already
offer two related representations; their authority and consistency need an explicit rule for output.

Graph-IR reaction export materializes its sides and uses their correspondence to assign map labels.
It must not infer a new atom map. A Reaction may fail to materialize a ReactionSpan; that failure
precedes molecular projection. Numeric source map labels and agents are not recoverable from the
graph-IR Reaction alone. Export can produce an empty agent section and deterministic fresh map
labels, but cannot claim a boundary-exact roundtrip of discarded input metadata or arbitrary
reaction constraints.

Rust and Python should expose the same semantic operations. The existing Python from_smiles and
from_reaction_smiles ingestion paths provide the composition precedent; names and the extent of
direct Python boundary-type exposure remain for API discussion.

## Roundtrip guarantees and evidence

There are two laws, with different domains and equality relations:

1. Boundary parse/render follows doc 155: preserve the ordered boundary meaning, accept the result
   under the same IO configuration, and make syntax normalization idempotent. Source spans and
   arbitrary original spelling are not reconstructed. Raw TableIR equality is not this law.
2. Project/render/parse/raise/resolve preserves the supported graph-IR molecular meaning under the
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

Local formatting references include RDKit's SmilesWrite::FragmentSmilesConstruct in
[SmilesWrite.cpp](../materials/codes/rdkit/Code/GraphMol/SmilesParse/SmilesWrite.cpp), which consumes a
traversal stack and emits atom/bond/ring tokens, and Indigo's traversal-relative tetrahedral mapping
in [smiles_saver.cpp](../materials/codes/Indigo/core/indigo-core/molecule/src/smiles_saver.cpp).
These establish useful writer structures to study. They do not establish umol projection semantics
or require their chemistry models or canonicalization policies.

## Decisions for the next discussion

- Confirm resolver-owned, nonmutating projection and the split between semantic projection and
  traversal-relative output notation.
- Choose common explicit stereo configurations versus output-specific population, particularly
  the missing stereo-bond frame representation and the role of retained source annotations.
- Define the supported molecular roundtrip domain and the model dependence of aromatic and
  valence reconstruction, without silently stripping unrepresentable state.
- Settle checked boundary construction, IO-config-dependent rendering failure, and reaction
  map-field authority. Then settle public names and Rust/Python exposure before sequencing work.
