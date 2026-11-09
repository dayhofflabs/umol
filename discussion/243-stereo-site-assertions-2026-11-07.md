# 243 — Frame-independent atom and bond stereo assertions

Status: In Progress
Date: 2026-11-07
Relates: [214](214-aggregate-frame-semantics-2026-08-28.md),
[226](226-smiles-roundtrip-design-2026-09-10.md),
[240](240-validation-constraint-ops-2026-10-18.md),
[242](242-stereo-transformers-2026-10-29.md),
[data-type guide](../docs/development/data-types.md),
[operations guide](../docs/development/operations.md)

## Scope

Restrict stored atom- and bond-level stereo assertions to frame-independent stereo constraints.
Configuration data is stored only on explicit stereo atoms and stereo bonds, together
with their ligand frames. A generic virtual stereo ligand lets those entities hold a
frame whose virtual ligand is not yet identified. SMILES/MOL raising populates the
entities; export conveys them directly to the format boundary. The DSL uses
explicit stereo records for configuration data. Configuration-bearing #T/#C
syntax is rejected in every DSL context.

This records the selected design direction and its consequences. S0–S3 and
S4a–S4b are complete, including the combined S3f verification gate. S4c's operation
documentation inventory is next; S4d remains pending.
The semantics, value interface, operation interfaces, and names below are the
working design. S2 records the implemented perception/resolution interfaces.
The shared resolver conflict/discharge policy and transformer contracts are
recorded below. The remaining subitems implement those contracts; missing plan
details do not reopen the semantics or authorize additional public interfaces.
Remaining subitems require explicit instruction. The staged implementation plan
follows the change inventory.

Molecule-level constraints and RelationalConstraint are outside scope. Reaction-span
canonicalization is separate work. This change does not remove kindless stereo
entities or redesign the existing transport of explicit stereo entities.
Extending the unfinished mol! / frag! syntax is outside scope.

## Motivation

A coset index is not an intrinsic scalar property of an atom or bond. Its integer
value denotes a configuration only relative to an ordered ligand frame. Changing
that order can change the index without changing the represented configuration.
Valence and aromatic-valence assertions do not have this dependence on participant
ordering.

The current #T/#C configuration payloads use implicit local frames. Consequently,
atom or bond attributes can change meaning when atom ids are reordered, even
though their stored values remain identical. This breaks the assumption that
attributes can be compared or combined independently of the graph identification.
Hydrogen folding is the only operation that transports these cosets. Remapping,
pushout, composition, edits, and reaction application treat them as ordinary values.

Current behavior on C#T0 with neighbors Cl, F, Br, I, and on the corresponding
molecule with the first two neighbors exchanged:

- Canonicalization relabels the neighbors without changing the stored coset, so its
  result perceives as the enantiomer of its input.
- Canonical equality reports the enantiomer pair equal and the same configuration in
  another numbering unequal.
- Substructure matching compares the stored values directly: the enantiomer pattern
  matches, and the pattern with the same configuration does not. Doc 103 specified
  matching through the pattern-to-target incidence mapping; that was not implemented.

No scalar attribute can repair this. Doc 103 records that a scalar stereo descriptor
is either local and index-dependent, or index-independent and global. A global
descriptor is still not pointwise. In 1,4-dimethylcyclohexane, an automorphism fixes
both ring carbons and exchanges the two ring branches at each, so every order-based
parity flips at both sites. One molecule then has two value assignments, and a
pointwise comparison cannot identify them because that automorphism moves no atom.

These payloads also serve as a channel for importing stereo data from external
formats. That engineering convenience does not justify making ordinary atom and
bond attributes frame-relative. The graph model already has explicit stereo
entities that carry both configuration and the frame needed to interpret it.

## Reduced assertions

A payload belongs on an atom or bond when its meaning is complete without knowing
the neighbors or their order. #v4 qualifies because it is a scalar with the
neighbors summed inside it. #T! qualifies because it states that no stereo of that
kind sits at the site. #T0 does not, because it selects one of two arrangements of
ligands that it does not identify. Among AtomConstraintForm and BondConstraintForm
keys, TetrahedralStereo and CisTransStereo are the only ones with frame-relative
payloads, so no other key changes.

Every atom and bond attribute is then pointwise. Framed payloads exist only on
ordered relations. Existing frame transport remains necessary for explicit stereo
entities and other positional payloads, including electron-count vectors. The
boundary follows from the carrier type: no key inside an atom or
bond constraint store behaves differently from the others.

For each supported stereo kind, an atom or bond assertion can be unconstrained,
negative, or positive. A positive assertion states stereo of that kind with the
configuration unspecified. For the existing constraints, #T concerns tetrahedral
stereo at an atom and #C concerns cis/trans stereo at a bond. A negative assertion
excludes that particular kind; it does not exclude every possible stereo kind.

The positive stereo constraint currently written #T+ or #C+ illustrates the intended meaning.
Coset literals, finite coset sets, variables, and configuration operators cease to
be stored payloads of these assertions. StereoCoset and StereoTerm remain available
on explicit stereo entities. Configuration-valued conversions into atom/bond stereo
constraints cease to be valid.

This reduces the current TetrahedralStereoForm and CisTransStereoForm.
StereoAtomForm, StereoBondForm, and StereoConfigurationForm retain their
configuration data.

These stereo constraints do not classify a site as stereogenic, prochiral, or symmetric.
An explicit stereo entity can describe an arrangement without asserting
stereogenicity. Stereo constraints and derived stereogenicity retain their
distinct meanings. Stereo descriptions can be placed on every carbon, including
CH3 groups, provided their ligand frames satisfy the existing integrity rules.
Neither a stereo atom nor a positive stereo constraint implies stereogenicity.

Ordinary normalized equality and lattice operations suffice for the reduced
assertions. For a fixed kind, positive and negative assertions have an incompatible
meet and an unconstrained join. Removing coset-specific branches requires no new
equality, matching, or lattice semantics. Atom and bond relabeling cannot change
the stereo constraint's value merely by changing a local ordering.

### Value interface

TetrahedralStereoForm and CisTransStereoForm remain named enums with the unit
variants Undetermined, NotStereo, and Stereo. The stored forms are not replaced
by an option. The compact spellings #T*, #T!, and #T+ correspond to the EDN
values :undetermined, :not-stereo, and :stereo. Configuration maps are not accepted
as standalone stereo-constraint values.

The ground enums TetrahedralStereo and CisTransStereo are removed. Both forms
implement AsLit with Lit = bool:

| Form | as_lit result |
| --- | --- |
| Undetermined | None |
| NotStereo | Some(false) |
| Stereo | Some(true) |

Both forms accept boolean construction through From<bool>: false produces
NotStereo and true produces Stereo. Coset-specific methods such as
matches_value(u32) are removed; ordinary lattice matching compares these constraints.

TetrahedralConfiguration (Ccw/Cw) and CisTransConfiguration (Z/E) remain
configuration types for explicit stereo entities. Their conversions into the
atom/bond constraint forms are removed. A configuration is not silently converted
to a positive stereo constraint.

Both named configuration enums gain From conversions into StereoConfigurationForm,
carrying their kind and the existing 0/1 coset mapping. Python's existing
configuration shorthand on explicit entities delegates to those Rust conversions.

Python preserves the same forms with payload-free variants. Constraint inputs
accept bool or the corresponding form, and as_lit returns bool or None. The ground
enum bindings and configuration-valued constraint coercions are removed.

### Uses of stereo constraints

Stereo constraints behave like the derived predicates #D and #R: computed from structure,
assertable in patterns, and pointwise.

- As pattern predicates, C#T+ matches an atom carrying a tetrahedral stereo entity,
  and C#T! matches an atom without one, including as single-atom patterns. This
  includes symmetric and prochiral entities. Bond stereo constraints behave the
  same way. A pattern entity needs a frame and has no negative form.
- A negative stereo constraint excludes a site from recognition and remains on the molecule.
  StereoSelection provides a comparable exclusion only per call.
- A positive stereo constraint requests an entity from StereoResolver and bypasses the
  eligibility check in recognition.

### Uses of the former payload

| Use | Home after this change |
|---|---|
| Staging imported stereo | Explicit entity, with a Virtual ligand where needed |
| Reading a coset in a chosen frame | The coset_for reading, computed on read and never stored |
| Query predicate on configuration | Pattern stereo entity, aligned through the match |
| Stereo kind at a site | Stereo constraints |
| Relative stereo through variables and operators | StereoTerm on explicit entities |
| Local literal in a molecule | Explicit stereo entity with its stated ligand frame |
| Assertion stored beside an entity and cross-checked against it | None; the entity is the only statement of its configuration |
| Coset on a standalone atom or bond value | None; such a value has no frame |

## Generic virtual stereo ligands

External formats can supply a partial frame. A SMILES slash pair names one reference
substituent at each endpoint, and a MOL wedge at a three-neighbor atom leaves the
fourth ligand unidentified. Which virtual ligand completes the frame depends on atom
state that is not known before resolution: F/C=C/F needs implicit hydrogen, C/N=N/C
needs a lone pair, and C/C=N/O needs one of each. The guess matters. An azo stereo
bond with implicit-hydrogen ligands resolves to a contradiction, and the same bond
with lone-pair ligands resolves. Raise has no chemistry model and does no chemical
interpretation.

StereoLigandKind gains a Virtual variant alongside Atom, ImplicitHydrogen, and
LonePair. Its DSL spelling is :virt. It denotes a generic virtual substituent
whose specific kind has not yet been identified. Undetermined remains reserved for
lattice states; StereoLigandKind has no Lattice implementation.

Virtual is a distinct ligand kind, not a wildcard. Ligand equality and frame
alignment remain exact in matching, gluing, canonicalization, and old-state
checks. No ligand lattice or separate compatibility operation is introduced.

A Virtual ligand is borne by its site atom or bond endpoint, like the other
virtual ligands. Frames remain complete, bounded, and pairwise-distinct, so a site or
endpoint holds at most one Virtual ligand. Remapping and frame reordering preserve the ligand kind.
Integrity accepts the new kind where it accepts virtual ligands.

StereoPerceiver::derive considers H/LP replacements compatible with the atom's
counts, fixed ligands, and the supported frame. One admissible completion identifies
the kinds; several leave Virtual in an explicitly Underdetermined derivation;
none produces the existing entity-failure diagnostic. Insufficient evidence also
returns Underdetermined. StereoResolver propagates that outcome independently of
the failure policy and leaves an in-place input unchanged.
Only impossible completion uses Error/Keep/Remove. Partial counts alone establish
neither ambiguity nor impossibility.

For a unique completion, StereoResolver writes the existing stereo-ligand replacement
edit. Positions, configurations, and entity constraints are preserved. Identification
does not change the atom's counts. This is a model-dependent transformation, not a
matching relation.

A pattern containing Virtual must be resolved before it can match a host whose
corresponding ligand has a specific virtual kind. Structural operations do not
perform that identification implicitly.

## DSL representation

#T and #C accept only the three states: unconstrained, negative, and
positive. Literals such as #T0 and #C1, coset sets, variables, and configuration
operators are rejected in every context, including molecule maps, reaction LHS
maps, standalone values, aliases, edits, deltas, and reaction-span sides.
Configuration-valued EDN constraint maps are rejected as well.

Configuration data is written on explicit :stereo-atoms / :stereo-bonds records,
with ordered ligands and attributes such as Th0 or Ct0. Their full configuration
expression grammar remains available. AtomDsl and BondDsl continue to wrap their
forms transparently, and molecule input entries keep their existing attribute
representation. No configuration payload is extracted from an atom or bond
string to create another entity.

Rendering uses the existing explicit stereo records. There is no implicit-frame
inference, inline/explicit combination rule, extra stereo id allocation rule,
conditional inline rendering, or stereo-sugar configuration. Existing entity
keywords, atom aliases, and defaults retain their current roles.

### Decision and migration cost

The compact spelling was considered and excluded on 2026-11-08. It would make an
attribute token construct an ordered relation using surrounding topology, then
require additional rules to recover the spelling during rendering. Preserving
the direct DSL-to-IR mapping and avoiding that machinery outweighs the shorter
input. A future compact representation can be considered as separate work.

The source audit found no affected whitepaper listings. Repository tests contain
59 EDN fixtures and 179 inline or parameterized molecule examples using the old
spelling, including repeated input/expected values. Some cases already require
rewriting under the representation change. Migrate configuration descriptions
to explicit records; use the reduced stereo constraints when a test asserts no
configuration.

## Raising SMILES and MOL data

Raise must construct stereo atoms/bonds with their ligand frames and configuration
data, instead of placing configuration indices in #T/#C for later materialization.
This explicitly reverses the earlier instruction to use those constraints as an
import channel.

The [TableIR raise](../umol-io/src/table_ir/raise.rs) maps explicit tetrahedral table
frames, wedge-derived tetrahedral configurations, and cis/trans bond frames
directly into graph-IR stereo entities.

Source ordering, endpoint reference choices, and actual/virtual ligand identity
must be accounted for while constructing the graph-IR frame. Raise stays mechanical:
where the source does not identify a virtual ligand, the frame holds a Virtual
ligand and resolution identifies it. Existing requirements for complete,
bounded, pairwise-distinct ligand frames remain in force. A failure to construct a
valid frame is handled at the conversion boundary; it is not a reason to retain a
frame-dependent atom attribute or silently discard configuration.

Absent stereo information and explicitly indefinite stereo remain distinct. Imported
stereo with an unspecified configuration must retain that assertion, including the
distinction represented by TableIR's Either bond configuration. No concrete coset
is chosen merely to make conversion succeed.

## Shared conflict and discharge policy

Valence, aromaticity, and stereo resolution use the same policy. Constraints
restrict admissible resolutions: unresolved data may be refined, and candidates
must satisfy the constraints. Resolution must not overwrite explicit data or
delete a constraint merely to hide a disagreement. A definite conflict between
independently meaningful descriptions is contradictory under every failure policy.

Discharge removes a constraint only when the resulting structure independently
entails it. For a normalized constraint C and its independent structural reading D:

- An incompatible meet is a conflict.
- If D meets C to D, the constraint adds no information and is discharged.
- If the readings are compatible but C adds information, retain the constraint.

Use the open derived reading of the planned result for this test. An evidenced
reading includes the assertion itself and cannot prove its redundancy. Closing
an absent relation to a negative reading cannot justify discharge either. The
rule applies equally to newly constructed and pre-existing entities; creation
history does not determine whether a constraint is retained.

| Constraint | Independent structural information | Result |
| --- | --- | --- |
| #v4 | Valence 4 | Discharge |
| #v4 | Valence 3 | Conflict |
| Atom #a1 | Aromatic-system contribution 1 | Discharge |
| Atom #a1 | Aromatic-system contribution 2 | Conflict |
| #T+ or #C+ | Entity of the asserted kind, with any coset | Discharge |
| #T! or #C! | Entity of the excluded kind | Conflict |
| #T+ or #C+ | No entity | Construct a compatible entity or retain the constraint; absence alone is not a conflict |
| #T! or #C! | No entity | Retain; the open derived reading is undetermined |

Realizability failures are separate from conflicts. Where an existing failure
policy offers Keep or Remove, it may retain or explicitly remove the failed
source. It must not choose a winner between independently conflicting
descriptions. Ambiguity is underdetermination, not a realizability failure.

S3c, S3d, and S3e implement this policy for stereo, valence, and aromaticity.
Earlier S2 completion notes record the behavior implemented at that stage.

## Resolution and stereo transformers

StereoResolver remains useful for evidence-driven resolution of stereo constraints. A positive
stereo constraint can cause it to construct a stereo entity of the asserted kind with an
undetermined coset. Imported configurations already exist on entities after raise,
and StereoResolver replaces their Virtual ligands with the identified kinds.
Shared perception and consistency checks must read the reduced stereo constraints
without extracting configurations from them.

The stereo transformers operate on explicit stereo entities, including before resolution.

In StereoPerceiver::recognize and StereoAdder, read positive and negative
#T/#C as kind constraints. Positive assertions supply kind evidence and
negatives exclude sites. New entities have undetermined cosets; existing
entities and all input assertions remain intact. Preserve rejection of
contradictory selected evidence and the existing required-input errors.

StereoRemover removes only selected entities and clears the corresponding
#T/#C assertion only where an entity is removed. Keep assertions at sites
without a removed entity and keep unrelated assertions. Preserve All and
Undetermined mode semantics and the shared selection rules.

StereoEnumerator obtains domains only from existing entities' configurations.
Remove intersections with configuration domains formerly stored in #T/#C;
retain checks of the reduced scalar constraints and entity-local constraints.
Preserve scalar assertions, frame positions, kinds, unselected data, assignment
order, multiplicity, iterator ownership, and the existing rejection of selected
unsupported symbolic configurations. An assertion alone supplies no enumeration
target.

HydrogenFolder and HydrogenUnfolder leave #T/#C unchanged. Delete scalar
coset transport. Convert actual hydrogen ligands and ImplicitHydrogen ligands
at the corresponding explicit frame positions, preserving configurations and
entity-local constraints. Retain explicit frame transport where positions
change. Virtual is never identified as hydrogen by these operations.

Apply these same contracts to inputs before resolution, including assertions
without entities and partial entity data. Retain the existing operation
interfaces, selectors, configs, and error variants.

StereoResolveConfig merges its constraint and entity failure policies, leaving
stereo_atom_failure and stereo_bond_failure. Each covers both a failed stereo constraint and a
failed entity at its site kind: Keep retains the failed value, and Remove clears the
stereo constraint or removes the entity. StereoInconsistency keeps its separate variants. At
ingest every entity comes from the import, so a caller distinguishes import failures
from stored failures by the config passed to each call.

The shared conflict/discharge policy above governs stereo resolution as well as
valence and aromaticity resolution.

Fixtures that convey configurations through #T0/#T1 or #C0/#C1 must be
rewritten using explicit stereo entities. Cases without configuration data use
the reduced stereo constraints.

## Standalone projection and export

Keep a standalone stereo-to-constraint projection as a deliberately lossy operation.
For a known tetrahedral or cis/trans kind, it replaces the entity by the corresponding
positive stereo constraint and discards its ligand frame, configuration, and
entity-local constraints. Every coset is accepted, including undetermined, finite-set,
and symbolic configurations, because the configuration is discarded. Different
configurations can therefore have the same projection. No reference-frame derivation
is needed, and recovery of the original configuration is not promised.

Unsupported or undetermined stereo kinds remain rejected. Unrelated atom and bond
data are preserved. Projection meets the positive constraint with any constraint
already at the site; a negative constraint is a conflict. A Virtual ligand does
not prevent this projection, since the whole frame is discarded.

That operation cannot remain a stage of the overall projection used by Convey.
Reducing a configuration to a stereo constraint there would lose information before the output
boundary has a chance to encode it. Other graph-IR projection stages must retain
the stereo entities needed by export.

The configuration path changes from

```text
graph-IR stereo entities -> #T/#C cosets -> TableIR stereo records -> SMILES
```

to

```text
graph-IR stereo entities -> TableIR stereo records -> SMILES
```

[Convey](../umol-graph/src/export.rs) reads explicit entities directly. Tetrahedral
conversion carries their ligand frames and winding into TableIR. Cis/trans conversion
chooses the boundary's actual reference substituents and transports the configuration
to those references. TableIR already has carriers for the supported output.

When encoding a stereo entity, Convey rejects the Virtual stereo ligand kind
rather than identifying it during export. This concerns ligand kind, not an
undetermined coset; the existing output support rules for cosets still apply.
TableIR is not extended for the new ligand kind in this change. If a later use case
requires it at that boundary, TableIR can gain the variant and format writers can
reject it when SMILES or MOL cannot express it.

Convey also rejects entity-local constraints that TableIR cannot preserve. Unlike
the standalone lossy projection, it cannot discard them merely because it can
encode the configuration. A positive #T/#C constraint can be discharged by an
encoded entity of that kind; it does not cause entity construction during export.
Remaining constraints without a TableIR representation are rejected, including a
positive constraint with no entity and the negative stereo constraints.

Tetrahedral boundary output requires a supported determined winding.
Undetermined bond configuration uses TableIR's Either; the format writer retains
its own narrower output limits. Configurations not reducible to supported
boundary values remain rejected.

Hydrogen-omission and representability checks in Convey must inspect explicit stereo
entities instead of relying on configuration-bearing #T/#C. Molecular and reaction
SMILES export share this boundary conversion change. The projection/export design
in doc 226 remains historical evidence for the current implementation; its stereo
intermediate is the part this proposal replaces.

## Ring-size model and operation configuration

Ring-size limits are chemistry choices owned by StereoModel. They define which
stereo descriptions the selected model admits. RingConfig supplies ring-detection
algorithms; para_stereo_max_iterations bounds symmetry refinement when the model
enables para stereo. These algorithmic settings belong to operation configs.

| Owner | Field | Default | Scope |
| --- | --- | --- | --- |
| StereoModel | stereo_bond_min_ring_size: u32 | 8 | Minimum ring size admitting cis/trans stereo on a double bond; implemented in S2g. |
| StereoModel | stereo_bond_e_min_ring_size: u32 | 8 | Minimum ring size admitting a ring-relative E configuration; checked by conformance validation. |
| AddStereoConfig, RemoveStereoConfig, EnumerateStereoConfig, StereoValidateConfig | para_stereo_max_iterations: u32 | 16 | Maximum para-stereo symmetry-refinement iterations. |
| AddStereoConfig, StereoResolveConfig, StereoValidateConfig | ring_config: RingConfig | RingConfig::default() | Ring-detection algorithms. |

A zero minimum ring size imposes no ring-size restriction; acyclic bonds are
unaffected. Ring searches cover the model's threshold even when it exceeds the
default search bound. Nitrogen ring recognition retains its own search and is
independent of the stereo-bond threshold.

StereoPerceiver captures the selected StereoModel at construction. Its queries
receive algorithm settings, with no numeric model-parameter override:

```rust
StereoPerceiver::derive(
    &self, molecule: &Molecule, ring_config: RingConfig,
) -> Solution<StereoDerivation, StereoContradiction>

StereoPerceiver::recognize(
    &self, molecule: &Molecule, selection: &StereoSelection,
    automorphism_algorithm: AutomorphismAlgorithm, ring_config: RingConfig,
    para_stereo_max_iterations: u32,
) -> Result<StereoDerivation, StereoPerceptionError>

StereoPerceiver::classify(
    &self, molecule: &Molecule, sites: &BTreeSet<StereoSite>,
    automorphism_algorithm: AutomorphismAlgorithm,
    para_stereo_max_iterations: u32,
) -> Result<BTreeMap<StereoSite, Stereogenicity>, StereoPerceptionError>
```

Adder and resolver apply the model's minimum when creating missing entities.
Existing stereo bonds retain their frames, configurations, and constraints;
Virtual-kind identification and structural checks still apply. A filtered #C
request is retained. Filtering supplies neither an inconsistency nor a reason
to invoke StereoFailurePolicy. Remover and enumerator operate on supplied
entities rather than applying a creation cutoff.

Addition and removal pass their configured iteration limit through recognition
or classification. Enumeration passes it to both classification and assignment
checks. Validation uses its configured limit for graph symmetry. In graph IR,
GraphSymmetryConfig::max_iterations and the refinement counter are u32. Runtime
stereo calls take their limits from configs; 16 is the operation default.

Conformance validation checks both model ring-size limits without altering
supplied descriptions. RingTooSmall and RingTooSmallForE carry the topology
BondId, ring_size: u32, and minimum_ring_size: u32. The general minimum also
applies to positive #C evidence whose entity creation was filtered.

The E check considers every simple ring below its threshold, using that ring's
continuation ligands and existing frame transport. A raw coset number alone does
not identify a ring-relative E configuration; one bond can be Z relative to one
ring and E relative to another. A required unresolved configuration contributes
Underdetermined; established contradictions take precedence. Both findings use
StereoConformanceContradiction and retain the existing error channel.

Parsing, raising, Interpret, and Convey receive no separate ring filter.
Wrappers invoking resolution use its ordinary config and the selected chemistry
model. Export preserves supplied configurations, including small-ring stereo.

## Operation interfaces

These choices keep the operations on their existing owners and avoid adding a
framed-operations trait, a frame wrapper, or a second family of stereo engines.

### Ligand identity and resolution

StereoLigand retains its existing public operations. Matching compares complete
ligands by equality after id mapping and uses the existing frame transport for
stereo attributes. No StereoLigand::matches method or Python equivalent is added.
Identification of Virtual belongs to StereoPerceiver and StereoResolver, using
resolved atom-state evidence and the existing ligand-replacement edits.

### Perception and resolution

Remove configuration extraction and implicit-frame coset transport from
StereoPerceiver::derive. Construct missing requested entities with an
undetermined coset; identify Virtual ligand kinds on retained entities while
preserving frame positions, configurations, and entity-local constraints.
Compare scalar constraints by kind independently of coset interpretation or
virtual-ligand ambiguity. Coset validation remains with the validator.

In StereoResolver::plan_resolve, discharge positive constraints entailed by
either new or retained entities of the corresponding kind. Retain negatives
at sites without an entity. Reject definite kind conflicts under every
failure policy; use those policies only for realizability failures. Preserve
the existing determined/underdetermined planning and application contracts.

Rename StereoPerception to StereoPerceiver and retain StereoDerivation. The
existing stereo_atoms and stereo_bonds rows can carry both entities requested by
positive constraints and completed versions of supported stored entities. For a
stored entity, completion preserves its ligand positions, configuration, and
entity-local constraints; it
replaces only the Virtual ligand kinds. StereoResolver compares those rows with
the input and emits the existing ReplaceStereoAtomLigands or
ReplaceStereoBondLigands edits. No public completion-result type is needed.

derive returns Solution<StereoDerivation, StereoContradiction> directly. The
non-contradictory outcomes carry the derivation, including established findings;
an unresolved stored row preserves its original ligand kinds and attributes.
The resolver checks those findings
before propagating underdetermination. Retain StereoInconsistency and rename its
existing contradiction wrapper to StereoContradiction. Remove the empty
StereoResolveError and the outer Result from stereo derivation and resolution;
this design introduces no operational failure for them. S2 records the complete
working contract and dependent changes.

StereoPerceiver::recognize continues to return missing entities only, so
StereoAdder retains its existing-entity preservation contract. The public
derive_stereo_atom and derive_stereo_bond methods still accept a coset for
construction in their documented frames. They remain useful independently of
stereo constraints. Recognition and derivation use the captured model's ring-size
limit and the operation's RingConfig. Recognition and classification receive
para_stereo_max_iterations. StereoPerceptionError is unchanged. The coset_for readings retain exact
frame transport. Remove the Python perception class and StereoDerivation;
supported Python operations retain their configuration and diagnostic values.

Remove StereoPerceiver::skipped_stereo_bonds and
StereoDerivation::skipped_stereo_bonds in S2g. Intentional exclusion of a new
candidate supplies no addition; it must not be translated into a deletion of
stored assertions or entities. Share the ring eligibility decision inside the
perceiver. Other ring-dependent recognition rules and RingConfig remain.

### Conformance validation

Apply the tier-3 contract in the operations guide through
StereoConformanceValidator::validate. Consume the kind agreement and
realizability findings from StereoPerceiver::derive. For supported known
kinds, reject an explicit entity's unrealizability whether or not its site
carries a scalar assertion. Report entity failures through the existing
StereoInconsistency::StereoAtomFailure and StereoBondFailure variants, wrapped
by StereoConformanceContradiction::Inconsistency. Keep the existing scalar
constraint failure and mismatch variants for their respective findings.

In validate, establish each explicit entity's kind/frame/configuration
prerequisites before calling validate_rings or the symmetry checks in
validate_stereo_atom and validate_stereo_bond. Check admissible coset indices
directly and report invalid configurations through the corresponding entity
failure diagnostic. These prerequisites remain model-independent; ring and
symmetry conformance remain separate checks. Configuration uncertainty leaves
dependent checks underdetermined without suppressing a known kind conflict.

Retain public entry points, configs, diagnostic variants, and the existing
ring and symmetry semantics. ValenceConformanceValidator and
AromaticityConformanceValidator retain their semantics.
IncidenceConstraintInvariantsValidator retains DerivedKind::Derived and
DerivedComplete and their distinct absence readings. Established contradictions
take precedence over unrelated uncertainty; validation leaves every input unchanged.

Verify entity failures with and without site assertions, kind agreement with
unspecified or symbolic cosets, invalid literal and finite-set coset indices,
virtual-ligand ambiguity, and established conflicts on partial inputs. Retain
ring-limit, ring-relative E, stereogenicity, topicity, and ligand-symmetry
cases. Check exact outcomes through standalone and composite validation,
unchanged inputs on every outcome, and open/complete incidence validation
of a positive assertion without an entity.

### Projection

Keep projection on StereoResolver, but make its entry points associated
functions with no resolver instance or model:

```rust
StereoResolver::plan_project(&Molecule) -> Result<Edits, StereoProjectError>
StereoResolver::project(Molecule) -> Result<Molecule, StereoProjectError>
StereoResolver::project_into(&mut Molecule) -> Result<(), StereoProjectError>
```

Projection has a deliberately lossy, model-independent contract. Keeping the
existing owner avoids adding a StereoProjector solely to house three methods.
The direct Result removes a Solution wrapper that cannot report a distinct
resolution outcome. Planning remains read-only, and project_into leaves the
input unchanged on error.

StereoProjectError retains UnsupportedStereoAtom, UnsupportedStereoBond,
AtomAssertion, and BondAssertion. StereoAtomFrame and StereoBondFrame disappear
because this operation no longer derives or transports frames. Conflict with a
negative stereo constraint uses the existing assertion error.

Remove ProjectFlags::STEREO and ProjectError::Stereo from aggregate projection.
ProjectContradiction::Stereo also disappears; no remaining projection stage produces it.
Resolver::project and its other projection stages keep their interfaces and
preserve stereo entities. Convey keeps its existing public arguments because
the other projection stages and boundary conversion still require them.

## Other affected surfaces

| Surface | Implication |
| --- | --- |
| Rust stereo site forms | Keep the named three-state forms, use bool literals and construction, and remove the ground enums and coset-specific methods/conversions. Keep configuration types for explicit stereo entities. |
| Stereo ligands | Add the Virtual variant and :virt spelling. Update exhaustive kind dispatch, canonical keys, and Python enum conversion. Keep ligand identity and matching exact. |
| Atom/bond constraint readings | Project the stereo kind without reconstructing a neighbor-relative coset. A known stereo kind can give a positive stereo constraint even when its coset is undetermined. Existing open/complete reading distinctions still apply. A caller that needs a coset uses coset_for with a stated frame. |
| Validation | Check stereo constraint realizability and kind agreement with explicit entities. Validate those entities independently of accompanying scalar assertions, including the frame/configuration prerequisites needed by model checks. Preserve valence and aromaticity validator semantics. |
| Hydrogen folding/unfolding | Remove special transport of inline #T/#C cosets. Continue maintaining actual/virtual ligands on explicit stereo entities. Hydrogen folding currently has separate code for both representations. |
| Canonicalization | Remove coset data from atom/bond constraint keys. Those attributes then have their ordinary relabeling semantics; explicit stereo canonicalization still transports its frames and configurations. |
| Matching and gluing | Ordinary attribute operations compare reduced stereo constraints. Configuration-sensitive comparisons use explicit stereo entities and existing exact frame alignment and transport. |
| DSL parser | Reject configuration-bearing #T/#C in every context. Retain the current AtomDsl/BondDsl wrappers and molecule input entries; configuration syntax remains available in the explicit stereo-entity DSLs. |
| DSL renderer | Retain explicit stereo-record rendering. No implicit-frame inference, conditional inline output, or new serialization configuration. |
| Python bindings | Use bool or the corresponding form for constraint input, and bool or None for literal extraction. Remove ground enum bindings and configuration-valued constraint coercions. Mirror the merged failure policies and uniform DSL rejection. |
| Fixtures and documentation | Rewrite configuration-bearing atom/bond strings using explicit stereo records. Cases without configuration data use the reduced constraints. Update examples, properties, benchmarks, and API descriptions that assumed a stored coset. |
| Depiction | The current depiction code already reads explicit graph-IR stereo entities. It needs no corresponding architectural change. |

Atom/bond constraint edits and deltas, including those carried by reactions, inherit
the reduced payloads without new frame-dependent edit semantics.

## Change inventory

Inspected against code at 9715d5176 on 2026-11-07. This is a symbol and consumer
inventory, not an implementation sequence. Source paths below are relative to the
repository root. A retained interface can still require changes to its
implementation, documentation, fixtures, or exhaustive matches.

### Public value and construction surface

| Owner | Symbol or family | Disposition |
| --- | --- | --- |
| umol-graph-ir/src/ir/stereo.rs | TetrahedralStereoForm, CisTransStereoForm | Retain names and Undetermined / NotStereo / Stereo cases; make Stereo payload-free. |
| Same | undetermined(), not_stereo(), stereo(...), is_stereo() | Retain the constructors and query; stereo becomes stereo(). |
| Same | Normalize and Lattice implementations on the two forms | Retain the trait surface with three-state semantics; remove coset-specific normalization and domain operations from these implementations. |
| Same | AsLit on the two forms | Change Lit to bool; positive and negative are ground. Add `From<bool>`. |
| Same | matches_value(u32) | Remove from both forms. Do not replace it with another inherent comparison method; ordinary lattice matching compares these boolean constraints. |
| Same | TetrahedralStereo, CisTransStereo | Remove both ground enums and all conversions to or from them. Remove their re-exports in umol-graph-ir/src/ir.rs. |
| Same | `From<TetrahedralConfiguration>` for TetrahedralStereoForm; `From<CisTransConfiguration>` for CisTransStereoForm | Remove. Configuration input must not implicitly become a stereo constraint. |
| Same | TetrahedralConfiguration, CisTransConfiguration | Retain the named configurations; add `From<TetrahedralConfiguration>` and `From<CisTransConfiguration>` for StereoConfigurationForm, preserving kind and the existing 0/1 mappings. |
| Same | StereoConfigurationForm, StereoConfiguration, StereoCoset, StereoTerm; StereoAtomForm, StereoBondForm and updates | Retain the framed configuration and constraint interfaces. The boolean literal change does not extend to these types. Correct rustdoc that equates an arrangement with stereogenicity. |
| umol-graph-ir/src/ir/constraint/atom.rs and bond.rs | AtomConstraintForm::TetrahedralStereo, BondConstraintForm::CisTransStereo and corresponding keys | Retain these variants and keys with reduced forms. Their tetrahedral_stereo / cis_trans_stereo constructors keep `impl Into<Form>` and thereby accept bool. |
| Same | AtomConstraintsForm / BondConstraintsForm getters, setters, asserted readings, normalization and lattice methods | Retain signatures; propagate the reduced values and remove configuration-valued fixtures. |
| umol-graph-ir/src/ir/ligand.rs | StereoLigandKind::Virtual | Add. Preserve existing Atom / ImplicitHydrogen / LonePair meanings. |
| Same | StereoLigand::new, map, try_map; RelationParticipant implementation | Retain construction and reference-transport interfaces. Transport preserves Virtual as a distinct kind. A standalone ligand remains an open value; owning-molecule integrity checks its incidence. |
| umol-graph-ir/src/ir/molecule/spec.rs | StereoLigandArg | Add Virtual alongside the other virtual arguments. Existing stereo_atom / stereo_bond spec constructors fill the bearer as they already do for H and LP. |
| umol-graph-ir/src/ir/view/ligand.rs and view/stereo.rs | StereoLigandView, StereoAtomView, StereoBondView; ligand and coset_for readings | Retain signatures. Expose the new kind through the existing kind getter; coset_for remains exact. |
| umol-graph-ir/src/ir/edit.rs and delta.rs | Stereo additions/removals, ligand replacement edits, stereo deltas, form-valued atom/bond changes | Retain carriers and entry points. They accept the extended ligand enum and reduced constraint forms; old-state comparisons stay exact. |

### Operation and boundary surface

| Owner | Symbol or family | Disposition |
| --- | --- | --- |
| umol-graph/src/ops/resolve/stereo.rs | StereoResolveConfig | Retain the merged failure fields, both defaulting to Error. Add ring_config; the ring-size threshold comes from StereoModel. |
| Same | StereoFailurePolicy; StereoResolveContradiction; StereoResolveError | Retain the policy. Rename the contradiction wrapper to shared StereoContradiction in ops/stereo.rs; retain its StereoInconsistency payload. Remove the empty error type. |
| Same | StereoResolver::new, with_config, plan_resolve, resolve, resolve_into | Retain construction; resolution returns Solution directly, without the outer Result. Add missing requested entities and replace Virtual kinds on retained entities. Discharge entailed positive constraints beside both new and retained entities. |
| Same | StereoResolver::plan_project, project, project_into | Keep names; remove the receiver and return the direct Result shapes recorded above. |
| Same | StereoProjectError | Retain UnsupportedStereoAtom, UnsupportedStereoBond, AtomAssertion and BondAssertion; remove StereoAtomFrame and StereoBondFrame. |
| umol-graph/src/ops/valence.rs and ops/resolve/valence.rs | ValenceContradiction | Move the definition to ops::valence, matching AromaticityContradiction and StereoContradiction in their domain modules. Preserve variants, conversions, and messages; migrate imports and remove the former resolver re-export. |
| umol-graph/src/ops/resolve/valence.rs and ops/resolve.rs | ValenceResolver::admit; Resolver::plan_constitution | Retain candidate admission. Apply independent entailment to discharge of #v, #d, and #t in the aggregate plan. |
| umol-graph/src/ops/resolve/aromaticity.rs, ops/aromaticity.rs, and ops/resolve.rs | AromaticityResolver::plan_resolve, plan_discharge; AromaticityPerceiver::plan_aromatic_system; Resolver::plan_constitution | Discharge atom and bond aromatic assertions entailed by new or retained systems. Retain negative assertions whose only support is absence. Keep standalone and aggregate resolution aligned. |
| umol-graph/src/ops/resolve.rs | ProjectFlags, ProjectError, ProjectContradiction | Remove STEREO, Stereo and Stereo respectively. Retain non-stereo flags and variants; aggregate projection preserves stereo entities. Remove the uninhabited ResolveError::Stereo; retain ResolveContradiction::Stereo with the renamed payload. |
| umol-graph/src/ops/stereo.rs | StereoPerception; derive; StereoDerivation::stereo_atoms / stereo_bonds | Rename the actor to StereoPerceiver. derive returns Solution explicitly; retain row types for completed stored frames and requested additions. |
| Same | StereoPerceiver::skipped_stereo_bonds; StereoDerivation::skipped_stereo_bonds | Remove both surfaces in S2g. |
| Same and ops/stereo/recognize.rs | StereoPerceiver::derive, recognize, new, derive_stereo_atom, derive_stereo_bond, classify; StereoPerceptionError and StereoInconsistency | Pass ring_config to derive and recognize. Pass para_stereo_max_iterations to recognize and classify. Read chemistry thresholds from the captured model. Preserve structural support boundaries and undetermined cosets on newly recognized entities. |
| umol-graph/src/ops/transform/stereo.rs and ops/transform.rs | StereoAdder, StereoRemover, StereoEnumerator, their configs/errors/iterator; Transform convenience methods | Expose para_stereo_max_iterations on all three configs and forward it to every symmetry path. Adder uses the model threshold and its RingConfig. Preserve supplied descriptions. Remove source-constraint coset handling in S3. |
| umol-graph/src/ops/transform/hydrogens.rs | HydrogenFolder, HydrogenUnfolder and their errors/configs | Retain interfaces; remove inline-constraint coset transport and maintain explicit frames through existing edits. Virtual is not implicitly an H ligand. |
| umol-graph/src/ops/validate/stereo.rs | StereoConformanceValidator, StereoValidateConfig, StereoConformanceContradiction | Use ring_config and para_stereo_max_iterations from StereoValidateConfig. Validate both model ring-size limits through RingTooSmall and RingTooSmallForE. Compare scalar constraints by kind and validate explicit descriptions independently of site assertions, with direct configuration prerequisites. Retain public entry points and diagnostic families. |
| umol-graph/src/ops/model.rs | StereoModel | Own stereo_bond_min_ring_size and stereo_bond_e_min_ring_size, both defaulting to 8. Preserve the other model fields. |
| umol-graph-ir/src/ir/symmetry.rs | GraphSymmetryConfig::max_iterations; refinement counter | Use u32 for the limit and counter. Preserve the refinement algorithm. |
| umol-io/src/table_ir/raise.rs | `TryIntoIr<Molecule>` for &table_ir::Molecule; RaiseError | Retain public conversion shape and existing error vocabulary. Wedge and bond-frame paths produce explicit entities and report unrepresentable frames at this boundary. |
| umol-graph/src/ingest.rs | Interpret for Smiles / ReactionSmiles; ingest_smiles and ingest_reaction_smiles families | Retain entry points and errors. Resolved routes use the ordinary resolver config; no import-specific ring filter is added. Unresolved routes retain the raised representation. |
| umol-io/src/ctfile/parser.rs | parse_mol_to_ir, parse_mol_bytes_to_ir; TableIR parser families | Retain signatures. These are model-free parse/raise paths, so they retain small-ring stereo. There is currently no resolved MOL Interpret implementation; no new resolved MOL API is required by this inventory. |
| umol-graph/src/export.rs | Convey for Smiles / ReactionSmiles; export_smiles and export_reaction_smiles families; ConveyError | Retain public shapes. Read explicit frames and preserve their configurations, including small-ring stereo. Reject unrepresentable data through existing errors; add no ring-output filter. |

### DSL ownership and consumers

| Files or public surface | Required migration |
| --- | --- |
| umol-graph-ir/src/dsl/atom.rs, bond.rs | AtomDsl / BondDsl, update wrappers, predicate enums, and parse_atom / parse_bond / update entry points remain form-valued. Their standalone grammar accepts only the three stereo constraint states. Public AtomPredicate / BondPredicate need no configuration-payload variant. |
| umol-graph-ir/src/dsl/stereo.rs | Reduce the internal TetrahedralStereoDsl / CisTransStereoDsl wrappers and their parsers/formatters. Retain the full StereoAtomDsl / StereoBondDsl and configuration-expression grammar. The two site wrappers are not re-exported from dsl.rs. |
| umol-graph-ir/src/dsl/constraint.rs | AtomConstraintDsl / BondConstraintDsl and general constraint readers use reduced stereo values. Configuration expressions remain on explicit stereo entities. Molecule-level constraint semantics remain outside this work. |
| umol-graph-ir/src/dsl/molecule.rs | Retain the existing atom/bond input entry types, entity registration, and explicit stereo-record rendering. The shared atom/bond grammar rejects configuration-bearing #T/#C here as it does elsewhere. |
| umol-graph-ir/src/dsl/refs.rs and namespace.rs | Read and render :virt in tree and streaming paths. Preserve its bearer and kind in existing structural references and participant indexes, without new stereo id allocation rules. |
| umol-graph-ir/src/dsl/reaction.rs | LHS maps use the same reduced attribute grammar as molecule inputs. Delta values reject configuration-bearing #T/#C. Existing explicit stereo records and their reference handling remain. |
| umol-graph-ir/src/dsl/reaction_span.rs and edit.rs | Atom/bond sides and edit values reject configuration-bearing #T/#C. Explicit stereo records, participant references, and ligand-replacement payloads accept Virtual. |
| umol-graph-ir/src/dsl/config.rs; MoleculeDsl FromIr / IntoIr | Retain StereoDefault, AtomDefaults, BondDefaults, overrides, and aggregate defaults. There is no additional stereo-sugar configuration or graph-aware atom/bond attribute conversion. |
| umol-graph-ir/src/ir/molecule/spec.rs | Standalone atom/bond specs reject configuration syntax; explicit stereo spec terms remain available, including the new virtual argument. |
| umol-graph-ir/src/macros.rs | Existing DSL macros inherit molecule acceptance or standalone rejection through their owning parsers. |
| umol-graph-ir/spec/umol-dsl-spec.md | Specify the three values and :virt; remove configuration-bearing #T/#C from every context. Keep full configuration syntax on explicit stereo records and preserve tree/streaming parity. |

The same attribute grammar applies in every context. Full stereo configurations
roundtrip through the existing explicit entity representation; no molecule-only
interpretation of an atom or bond attribute is needed.

### Graph algorithms and existing framed infrastructure

| Files | Disposition |
| --- | --- |
| umol-graph-ir/src/ir/view/atom.rs, view/bond.rs, view/constraints.rs | Derive a positive form from the supported known kind without sorting ligands or transporting cosets. Preserve open/complete absence and unknown-kind behavior. Matching and consistency consume these revised readings. |
| umol-graph-ir/src/ir/canonicalize.rs | Remove cosets from site-constraint keys; add the distinct Virtual kind key. Preserve existing explicit-entity frame actions and configuration canonicalization. |
| umol-graph-ir/src/ir/substructure.rs | Keep exact ligand alignment and existing stereo-form transport/matching. Reduced atom/bond predicates use the revised readings; no new ligand compatibility rule is introduced. |
| umol-graph-ir/src/ir/molecule/integrity.rs and ir/stereo/integrity.rs | Admit Virtual at the bearing site/endpoint, retaining incidence, arity, bounded degree, and pairwise-distinctness checks. Existing construction and edit publication boundaries still enforce integrity. |
| umol-graph-ir/src/ir/symmetry.rs, coloring.rs, incidence.rs | Keep Virtual distinct in graph incidence and coloring; graph-derived ligand actions fix it as they fix other virtual kinds. Existing symmetry projection is representation-based and does not identify virtual kinds. |
| umol-graph-ir/src/ir/molecule/build.rs, editor.rs, editor/transact.rs, apply.rs, transact.rs; edit.rs, delta.rs | Carry the new enum value through existing construction, mutation, rollback, removal, and restoration paths. No new edit or transaction family. Test old/new ligand lists and integrity failures. |
| umol-graph-ir/src/ir/molecule/pushout.rs; compose.rs, correspondence.rs, reaction.rs, reaction_span.rs | Keep exact overlap and old-state frame alignment. Reduced attributes propagate through existing operations. Existing ligand equality and frame transport continue to apply. |
| umol-graph/src/ops/stereo.rs, ops/stereo/recognize.rs, ops/resolve/stereo.rs | Remove constraint-coset sourcing and general ring suppression. Preserve stored frame order during virtual-kind identification and preserve all entity attributes. |
| umol-graph/src/ops/transform/stereo/enumeration.rs | Remove coset-domain intersections and configuration checks sourced from #T/#C. Preserve explicit-entity configuration constraints, assignment traversal, and unrelated constraint handling. |
| umol-graph/src/ops/transform/hydrogens.rs | Delete the inline #T/#C transport branches and their now-unused local machinery; retain explicit ligand conversion and configuration transport when actual frame positions change. |
| umol-io/src/layout/stereo.rs, layout/coordgen.rs, depict/molecule.rs and depict/verify.rs | Preserve explicit-entity depiction and verification. Extend exhaustive virtual-kind dispatch without identifying Virtual as H or LP. Named CisTransConfiguration remains used here. |

umol-perm, domain-independent graph-core algorithms, geometric models, and TableIR's
public stereo carriers require no new abstraction. The SMILES/CTfile parsing and
SMILES rendering architectures keep their boundary model; conversion expectations
and end-to-end cases change. Generic EDN parsing does not interpret #T/#C strings.

### Python surface and parity

| Files | Required migration |
| --- | --- |
| umol-py/src/stereo.rs | Make the two site forms payload-free; return bool or None from as_lit; remove ground enum classes and conversions. Add Virtual. Retain named configuration enums and explicit-entity inputs; update repr, hashing, equality, lattice methods, and variant construction consistently. |
| umol-py/src/constraint/atom.rs and constraint/bond.rs; atom.rs and bond.rs | TetrahedralStereoLike / CisTransStereoLike accept bool or the corresponding form. Accept True as positive; remove named-configuration coercions. Migrate container/view setters and constraint constructors together. |
| umol-py/src/resolve.rs and transform.rs | Expose ring_config on StereoResolveConfig and para_stereo_max_iterations on all three transformer configs. Keep constructors, getters, Rust conversions, repr, and defaults aligned; enforce the u32 range. Retain merged failure-policy fields. |
| umol-py/src/model/stereo.rs | Expose stereo_bond_min_ring_size and stereo_bond_e_min_ring_size on StereoModel construction, getters, conversion, equality, and repr; enforce the u32 range and reject the former stereo_bond_minimum_ring_size keyword. |
| umol-py/src/stereo/perception.rs | Remove StereoPerception and StereoDerivation with their constructors and queries. Retain StereoSite, StereoSelection, and StereoInconsistency used by supported operations. |
| umol-py/src/lib.rs and python/umol/__init__.py | Remove the perception and derivation classes in S2a and both ground enums in S3 from registration, imports, and __all__. Keep configuration enums and existing operation bindings. |
| umol-py/src/molecule.rs and reaction.rs | Existing from_smiles / to_smiles and stereo transformer methods inherit the revised Rust behavior; keep signatures, ownership, copying, and defaults. Molecule DSL parsing enforces the same reduced attribute grammar as other contexts. |
| umol-py/src/edit.rs, delta.rs, reaction_span.rs; stereo.rs participant extraction | Existing enum conversion carries Virtual through additions, replacement edits, deltas, and span construction. Do not add Python-only stereo interpretation. |
| umol-py/src/error.rs | Reuse existing exception mappings for conversion, perception, and contradiction errors. Check removed variants and new conversion failures against the existing vocabulary. |

StereoResolver, ProjectFlags, and standalone projection currently have no Python
binding. Remove the existing Python perception/query surface without replacing it;
Python exposes the supported molecular operations and their required values.
Named-configuration shorthand on explicit entities is already present in Python;
the conversion clarification
in the value interface gives it a Rust-owned implementation.

### Verification, fixtures, and documentation

| Existing coverage or artifact | Migration and required coverage |
| --- | --- |
| umol-graph-ir/src/ir/stereo.rs tests; tests/property/lattice.rs, strategies.rs | Replace reduced-form coset generators with all three states. Preserve equality/lattice/ground-literal laws. Keep full configuration generators and independent coset references on explicit entities. |
| ir/ligand.rs, ir/stereo/tests, ir/view/constraints/tests/stereo_atom.rs and stereo_bond.rs; tests/integrity.rs | Cover Virtual identity, bearer checks, duplicate ligands, exact matching, and known-kind readings with undetermined cosets. |
| dsl/atom.rs, bond.rs, stereo.rs, molecule/tests.rs, reaction.rs, reaction_span.rs, edit.rs; tests/property/stereo/serialization.rs and molecule/serialization.rs | Reject configuration-bearing #T/#C consistently in every context. Test explicit stereo records with literal, set, variable, and operator configurations, defaults, aliases, and exact tree/streaming/render roundtrips. Preserve unrelated kindless-payload equivalence contracts. |
| umol-graph-ir/src/macros.rs tests | Check uniform rejection of configuration-bearing #T/#C through the existing DSL macros. |
| tests/property/molecule/{canonicalize,remap,reframe,meet_pushout}.rs; property/substructure.rs; property/reaction and property/edit families | Retain relabeling, frame-transport, pushout, reaction, edit, and rollback laws. Check same-configuration and enantiomer distinctions, including Virtual ligands, without weakening exact gluing or old-state checks. |
| umol-graph/src/ops/resolve/valence.rs, resolve/aromaticity.rs, resolve.rs; tests/property/resolve.rs | Cover the common conflict/discharge policy: independently entailed constraints, compatible retained information, definite conflicts, new and existing systems, and negative assertions without a supporting relation. Preserve candidate-selection and resolution laws. |
| umol-graph/src/ops/stereo.rs, recognize/tests.rs, resolve/stereo.rs and resolve.rs; ops/validate/stereo.rs | Cover direct imported entities, virtual-kind replacement with unchanged coset, merged policies, unconditional mismatch contradictions, removal of general ring suppression, and lossy standalone projection. |
| umol-graph/src/ops/validate/stereo.rs, valence.rs, aromaticity.rs, constraint/incidence.rs, and validate.rs | Cover kind agreement independently of coset uncertainty, explicit entity failures without site assertions, direct coset-range rejection, unchanged inputs, and composite validation. Preserve valence/aromaticity outcomes and the tier-2 open/complete distinction. |
| ops/transform/stereo/tests, stereo/enumeration/tests.rs, hydrogens/tests; tests/property/transform/stereo and hydrogens.rs | Preserve selection, enumeration, reconstruction, and existing-entity laws. Move configuration descriptions to explicit entities. Isolated atom examples that only test the stereo constraint use #T+. |
| umol-io/src/table_ir/raise.rs; CTfile and SMILES parser/render tests; umol-graph/src/ingest.rs and export.rs; tests/property/project.rs and resolve.rs | Update expected imported carriers. Cover alkene H, azo LP, mixed imine ligands, Either/unspecified stereo, output rejection, ring-policy placement, explicit frame transport, and molecule/reaction roundtrips. Aggregate projection assertions must now retain entities. |
| umol-graph/tests/resolution/data and tests/kekulization/data | Migrate the 59 fixture files using the old syntax: 21 in stereo_tetrahedral, 19 in stereo_cis_trans, 17 elsewhere in resolution data, and 2 in kekulization data. Preserve their intended chemistry and update carrier and ring-policy expectations. |
| umol-py/tests/test_stereo.py, test_stereo_perception.py, test_constraint.py, test_constraint_readings.py, test_structural_pattern_matching.py, test_class_constructor_signatures.py | Remove perception/derivation query and construction cases; preserve operation-required value coverage. Cover unit variants, bool input/output, removed imports and coercions, Virtual, merged keyword signatures, readings, and exact exception/roundtrip behavior. |
| umol-py/tests/test_molecule_stereo.py, test_molecule_stereo_enumeration.py, test_import.py, test_smiles.py, test_smiles_export.py, test_resolve.py, test_workflow.py; edit/delta/transaction tests | Exercise the existing high-level entry points and participant carriers after rebuilding the extension. Other Atom/Bond/Molecule view and constructor fixtures that construct old forms migrate with them. |
| umol-graph/benches/resolve.rs, transform.rs and transform/stereo.rs | Update instance projection calls, Solution assumptions, aggregate output checks, and configuration-bearing input strings. Use explicit stereo records while preserving the measured operation and iterator lifetime controls. |
| umol-graph-ir/benches/canonicalize.rs, constraint_readings.rs, editor.rs, integrity.rs; IO parsing/layout/SVG benchmarks | Add or adapt representative Virtual and explicit stereo fixtures where those operations are measured. Preserve the measured operation; parser or entity-only fixtures are not interchangeable workloads. |
| umol-graph-ir/fuzz/fuzz_targets/fuzz_entity_strings.rs, fuzz_molecule.rs, fuzz_reaction.rs, fuzz_reaction_span.rs and fuzz_edits.rs | Preserve parser and tree/streaming parity laws. Configuration-bearing #T/#C inputs become rejection cases in every context; add valid explicit stereo-record and Virtual seeds. Keep generated corpora as historical inputs, not bulk migration targets. |
| umol-graph-ir/spec/umol-dsl-spec.md; docs/chemistry/representation.md; docs/development/data-types.md, operations.md, and python-api.md; relevant nomenclature entries | Synchronize carrier, conversion, frame, and binding descriptions. Record the shared resolver conflict/discharge policy in the operations guide. Remove wording that equates a stereo entity or NotStereo with stereogenicity. |
| .agents/skills/mol-notation/SKILL.md and its .claude redirect | Update canonical notation guidance to reject configuration-bearing #T/#C and document explicit stereo records and :virt. Preserve shared-skill ownership; update other examples only where their actual claims become false. |

The gated suites are graph-ir's property target with proptest; graph's property
target with proptest and resolution target with conformance; and IO's
smiles_property, table_ir_property, and mol_property targets with proptest.
IO's MOL/SDF/SMILES parser conformance modules require conformance as well.
CoordGen/layout and SVG paths additionally require coordgen/depiction. Existing
crate and workspace checks still apply. Python compilation and tests use the
repository Python 3.13 environment and a rebuilt extension. These are inventory
entries for later verification; no tests or benchmarks were run for this audit.
No additional runtime dependency is identified.

## Staged implementation plan

Each subitem identifies its implementation files and owns the corresponding tests
and documentation. Record implementation results beneath that subitem, as in
242. Completion and verification are recorded beneath each finished subitem.

Every stage ends with a passing build and tests. S3 contains the coordinated
removal of configuration-bearing site constraints: its callers must be migrated
in the same stage. The detailed change inventory above is the checklist for
those callers.

### S0 — Establish correctness references and measurements

- **S0a — Add independent stereo reference cases** — additive (green).
  Modules: umol-graph-ir/src/ir/canonicalize/tests.rs,
  umol-graph-ir/tests/property/stereo/reference.rs,
  umol-graph-ir/tests/property/substructure.rs. [dep: none]

  Define small tetrahedral and cis/trans examples using explicit stereo entities
  in two different ligand orders. Check that transporting the configuration
  preserves equality and matching, while changing the configuration distinguishes
  the stereoisomers. Historical reproduction details belong in Motivation;
  tests assert the current grammar and explicit-entity equality and matching.

  **Completed — 2026-11-08.** Added 40 canonicalization cases covering independent
  ligand orders, an atom-id swap, both para_stereo settings, and opposite
  configurations. The framed-equality reference checks all 128 valid combinations
  of kind, ligand permutation, source coset, and comparison coset in both
  directions, using the independent parity/endpoint-block oracle. The matching
  table covers seven arrangements, with and without atom relabeling, for equal
  and opposite configurations. It checks exact atom mappings in both directions
  under both matching strategies and all six graph algorithms.

  The inline-canonicalization defect is recorded in Motivation above. Production
  code and APIs are unchanged. Validation passed: all 72 new parameterized cases
  and the existing stereo-filtered library/property tests, using
  `cargo test -p umol-graph-ir --features proptest --lib --test property stereo`;
  `cargo clippy -p umol-graph-ir --features proptest --lib --tests -- -D warnings`;
  `cargo +nightly fmt --all`; and `git diff --check`.

- **S0b — Measure the existing operations** — additive (green).
  Modules: umol-graph-ir/benches/parsing.rs,
  umol-graph-ir/benches/rendering.rs,
  umol-graph-ir/benches/fixtures/stereo.rs,
  umol-graph-ir/benches/canonicalize.rs,
  umol-graph-ir/benches/constraint_readings.rs,
  umol-graph/benches/resolve.rs,
  umol-graph/benches/transform.rs. [dep: S0a]

  Select existing cases covering molecule parsing, rendering, canonicalization,
  constraint readings, import, projection, and stereo transformation. Use explicit
  stereo records for comparable baseline inputs and record old attribute-based
  workloads separately. Run correctness checks and record the comparison baseline
  and relevant measurement conditions here. Keep fixed explicit-entity workloads
  separate from import/parser workloads whose stored
  representation will change. Preserve the existing iterator lifetime controls.

  **Completed — 2026-11-08.** Added tetrahedral and cis/trans parsing/rendering,
  explicit projection/removal, and SMILES ingest cases; reused the existing
  canonicalization, readings, addition, and enumeration loops. All 55 selected
  paths passed benchmark test mode. Timing-loop boundaries are unchanged.

  Two runs on verified AC power replace the battery measurements. They use
  production revision eefd1daa64b8a597abd822b0f6ee0a04330d497a, rustc 1.99.1,
  macOS arm64, default features, and Low Power Mode off: 50 samples, 300 ms
  warm-up, and one-second measurement. Of the 55 paths, 54 agree within 2.71%.
  SmallInput iter_first varies by 12.27%; a longer repeat reduces this to 3.20%,
  versus 0.86% for PerIteration. Retain both controls when assessing enumeration
  changes; their setup and destruction lifetimes differ.

  Compare unchanged explicit-entity workloads directly. Treat old #T/#C inputs
  as historical once rejected, and report changed import/projection semantics
  alongside their timings. Library tests, stereo properties, benchmark Clippy,
  formatting, and diff checks passed.

**S0 gate:** the explicit-entity reference tests pass and the benchmark results
are reproducible. Rejection tests for the obsolete attribute syntax land in S3.
**Passed — 2026-11-08.**

### S1 — Add generic virtual ligands and configuration conversions

- **S1a — Give named configurations an explicit-entity conversion** — additive (green).
  Modules: umol-graph-ir/src/ir/stereo.rs,
  umol-py/src/stereo.rs,
  umol-py/tests/test_stereo.py. [dep: S0b]

  Add From conversions from TetrahedralConfiguration and CisTransConfiguration
  into StereoConfigurationForm. Preserve the current Ccw/Cw and Z/E mappings to
  their kind and coset index. Replace Python's duplicate mapping with these Rust
  conversions. Test all four named values through Rust and Python construction.
  Removal of their old constraint conversions belongs to S3a.

  **Completed — 2026-11-08.** Added both From implementations into
  StereoConfigurationForm. Ccw and Z map to coset 0; Cw and E map to coset 1,
  with the corresponding stereo kind. StereoConfigurationLike::to_rust now
  delegates to those Rust conversions. The constructor documentation names the
  accepted configuration enums. The two From implementations are the only public
  API additions.

  Exact tables cover all four Rust conversions, binding coercions, and Python
  StereoAtomForm/StereoBondForm construction. The existing full-form passthrough
  assertion is retained. Validation passed: 153 graph-IR configuration tests,
  12 Rust binding configuration tests, and 131 Python stereo
  tests after rebuilding the extension with maturin develop under Python 3.13.16.
  Clippy passed for both crates' libraries/tests with warnings denied; formatting
  and `git diff --check` passed.

- **S1b — Make Virtual a supported stored ligand kind** — breaking (red→green).
  Modules: umol-graph-ir/src/ir/ligand.rs,
  umol-graph-ir/src/ir/molecule/spec.rs,
  umol-graph-ir/src/ir/molecule/integrity.rs,
  umol-graph-ir/src/ir/canonicalize.rs,
  umol-graph-ir/src/ir/reaction.rs,
  umol-graph-ir/src/ir/symmetry.rs,
  umol-graph-ir/src/dsl/refs.rs,
  umol-graph-ir/src/dsl/molecule.rs,
  umol-graph-ir/src/dsl/edit.rs,
  umol-io/src/layout/stereo.rs,
  umol-io/src/depict/molecule.rs,
  umol-py/src/stereo.rs. [dep: S0b]

  Add StereoLigandKind::Virtual, StereoLigandArg::Virtual, and the :virt
  spelling in explicit stereo records and edit values. Fill the bearer through
  the existing construction paths. Update exhaustive enum conversions and the
  canonical key, integrity, layout, and depiction consumers together.

  Test correct and incorrect bearers, duplicate ligands, ordered frames,
  relabeling, compaction/restoration, and exact equality with the other ligand
  kinds. Verify that substructure matching still requires equal ligand kinds.
  Check tree/streaming serialization and Python enum conversion. Exercise
  the layout and depiction feature paths. Existing graph-derived ligand actions
  keep virtual ligands fixed; they do not identify their chemical kind.

  **Completed — 2026-11-08.** Virtual is supported in Rust/Python construction,
  molecule/edit DSL values (`[:virt bearer]`), integrity, canonicalization,
  layout, and depiction. Equality and matching retain exact ligand kinds.
  The DSL specification and development guides are updated.

  Tests cover invalid bearers and duplicates, frame transport and canonical
  relabeling, exact matching, tree/streaming serialization, native layout,
  depiction, and bindings. Property generators include Virtual with their
  existing laws unchanged.

  Validation passed: graph-IR and IO library tests, graph-IR properties,
  native layout tests, Rust binding tests, and the full Python suite after
  rebuilding with Python 3.13.16. Workspace Clippy passed for all targets with
  proptest/depiction enabled and warnings denied; formatting and diff checks
  passed. Identifying Virtual as hydrogen or a lone pair remains S2c–S2d.

**S1 gate:** the new kind can be constructed, serialized, transported, and depicted
through the supported interfaces. The affected Rust and Python tests
pass together. Configuration-bearing site forms still exist at this point.
**Passed — 2026-11-08.**

### S2 — Import stereo directly and resolve virtual ligand kinds

**Implemented — 2026-11-08.** S2a–S2i are complete and the stage gate passed.
Site-value reduction and the projection/export migration remain S3.

Rename the Rust actor to StereoPerceiver. Its policy-free derive operation
checks supported stored frames and requests missing entities from positive
constraints. Recognize continues to return missing entities only; classify
continues to assess existing arrangements. Their Rust query contracts and
StereoPerceptionError remain separate from the changed derive interface.

The shared derivation/resolution interfaces are:

```rust
StereoPerceiver::derive(&Molecule, ring_config: RingConfig)
    -> Solution<StereoDerivation, StereoContradiction>
StereoResolver::plan_resolve(&Molecule)
    -> Solution<Edits, StereoContradiction>
StereoResolver::resolve(Molecule)
    -> Solution<Molecule, StereoContradiction, ()>
StereoResolver::resolve_into(&mut Molecule)
    -> Solution<(), StereoContradiction>
```

StereoContradiction retains the existing Inconsistency(StereoInconsistency)
wrapper. Remove the empty StereoResolveError and the outer Result; no
StereoError replaces it. Remove the resulting uninhabited ResolveError::Stereo
branch from the aggregate resolver and its adapters.

StereoDerivation keeps its entity row types and finding fields; the old
skipped-site field is retired in S2g:

```rust
stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>
stereo_bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>
inconsistencies: Vec<StereoInconsistency>
```

Retain all six StereoInconsistency variants: TetrahedralStereoFailure { atom },
StereoAtomFailure { stereo_atom }, TetrahedralStereoMismatch { atom, stereo_atom },
CisTransStereoFailure { bond }, StereoBondFailure { stereo_bond }, and
CisTransStereoMismatch { bond, stereo_bond }. Derivation collects findings;
resolver policy decides recovery or rejection.

Identify Virtual only by H/LP replacements consistent with available counts,
fixed ligands, topology, and the supported complete frame. One established
completion identifies the kinds. Several completions, or insufficient evidence
to establish the answer, produce Underdetermined. Ruling out every completion
produces the existing entity-failure finding. Partial counts alone imply
neither ambiguity nor impossibility. No ligand lattice or compatibility API
is added.

Both non-contradictory outcomes carry StereoDerivation. An underdetermined
derivation retains established rows, unresolved stored frames in their original
form (including Virtual where present), and
established inconsistencies; it does not enumerate alternative molecules.
Stored rows preserve ligand order, configuration, and all entity constraints.
Identification changes only ligand kinds, never atom counts or cosets. An
undetermined coset is a valid stored state, not by itself an underdetermined
identification result. Stereo descriptions do not imply stereogenicity.

Resolver policy applies only to impossible completion. Error rejects, Keep
retains the failed description, and Remove clears the failed constraint or
removes the failed entity. Mismatches remain contradictory under every policy.
Check established findings before propagating underdetermination, so unrelated
ambiguity cannot conceal a known rejection. Otherwise return an empty
underdetermined edit plan, without publishing partial changes. A consuming
unresolved outcome returns unit; an in-place unresolved outcome leaves the
input unchanged. Planning establishes applicability and integrity preservation;
application makes no additional chemical decisions.

Remove the Python perception class, its query methods, and StereoDerivation.
Keep supported molecular operations and their required models, configurations,
selectors, and diagnostics. Do not add a replacement perceiver binding or a
Python derivation-result tuple. Existing constraint payload handling remains
until S3; model-owned ring-size filtering is implemented in S2g.

- **S2a — Remove the Python perception and derivation interface** — breaking (red→green).
  Modules: umol-py/src/stereo/perception.rs, umol-py/src/stereo.rs,
  umol-py/src/lib.rs, umol-py/python/umol/__init__.py,
  umol-py/tests/test_stereo_perception.py, umol-py/tests/test_import.py,
  and docs/development/python-api.md. [dep: S1b]

  Remove StereoPerception and StereoDerivation construction, queries,
  registration, and exports. Retain StereoSite, StereoSelection, and
  StereoInconsistency, including their operation-facing conversions and error
  payloads. Update the Python API inventory without introducing another query
  surface. Preserve tests of retained values; remove cases whose sole subject
  is a deleted query or result carrier. Check imports and the existing
  addition, removal, enumeration, and resolution bindings.

  **Completed — 2026-11-08.** Removed the Python StereoPerception and
  StereoDerivation wrappers, registration, exports, and query/construction tests.
  Retained StereoSite, StereoSelection, StereoInconsistency, their value tests,
  and the error adapters used by operations. The export inventory verifies that
  both retired classes are absent from the package and native module. Updated
  the Python API guide. The Rust perception/resolution implementation is unchanged
  by this subitem.

  On Python 3.13.16, `cargo test -p umol-py --lib --quiet` passed 1,773 cases
  with two ignored. After `maturin develop --locked --manifest-path
  umol-py/Cargo.toml`, `pytest -q umol-py/tests` passed 3,138 cases with two
  skipped. `cargo clippy -p umol-py --all-targets --no-deps -- -D warnings`,
  formatting, and `git diff --check` passed. Dependency-inclusive Clippy still
  fails on type_complexity in the unfinished virtual-ligand identification code
  in umol-graph/src/ops/stereo.rs; that code is outside this subitem and unchanged.

- **S2b — Rename the Rust actor and align domain contradiction types** — breaking (red→green).
  Modules: umol-graph/src/ops/stereo.rs,
  umol-graph/src/ops/stereo/recognize.rs, umol-graph/src/ops/resolve/stereo.rs,
  umol-graph/src/ops/valence.rs, umol-graph/src/ops/resolve/valence.rs,
  umol-graph/src/ops/resolve.rs, umol-graph/src/ops/validate/stereo.rs,
  umol-graph/src/ops/transform/stereo.rs, umol-graph/src/ingest.rs,
  umol-graph/benches/transform/stereo.rs,
  umol-graph/tests/property/stereo.rs,
  umol-graph/tests/property/transform/stereo.rs,
  docs/development/data-types.md, docs/development/nomenclature.md,
  and the owning unit tests.
  [dep: S2a]

  Rename StereoPerception to StereoPerceiver throughout Rust callers and
  rustdoc. Rename StereoResolveContradiction to StereoContradiction and own it
  in ops/stereo.rs alongside StereoInconsistency. Retain the wrapper and all
  six finding variants. Preserve recognition/classification contracts and
  diagnostics; add no selectors, query methods, or compatibility aliases.
  Run the affected perception, transformer, and resolver cases with the
  renamed symbols.

  Use one placement convention for the parallel domain contradiction types:
  ValenceContradiction in ops::valence, AromaticityContradiction in
  ops::aromaticity, and StereoContradiction in ops::stereo. Move
  ValenceContradiction from the resolver and migrate its consumers, preserving
  its variants and conversions. Remove its former resolver re-export.

  **Completed — 2026-11-08.** Renamed the actor to StereoPerceiver throughout
  Rust code, rustdoc, tests, properties, and the benchmark caller. Moved the
  renamed StereoContradiction to ops/stereo.rs and ValenceContradiction to
  ops/valence.rs, matching AromaticityContradiction's existing placement.
  Migrated resolver and ingest imports and documented the convention in the
  nomenclature guide. Retained the wrappers, their messages and variants, and
  the existing query/operation result shapes. StereoPerceptionError is
  unchanged; no compatibility aliases or re-exports were added.

  Validation passed: graph library, property, and resolution suites with
  `cargo test -p umol-graph --features proptest,conformance --lib --test property
  --test resolution --quiet`; all targets, including benchmarks, with
  `cargo check -p umol-graph --all-targets --features proptest,conformance`;
  and 1,773 Rust binding tests with two ignored on Python 3.13.16. Formatting
  and `git diff --check` passed. Review against the pre-subitem source confirmed
  that changes are limited to the two renames, type relocations, imports,
  documentation, and formatting; method bodies and test assertions are preserved.
  The type_complexity lint recorded under S2a remains in the unchanged
  virtual-ligand identification implementation.

- **S2c — Derive identified frames and explicit underdetermination** — breaking (red until S2d).
  Modules: umol-graph/src/ops/stereo.rs and its tests. [dep: S2b]

  Implement H/LP identification and the direct Solution return from derive.
  Use the existing row types for requested additions and supported stored
  frames. Preserve stored positions, configurations, constraints, and atom
  counts. Retain unresolved rows and known findings in the underdetermined
  payload; do not infer the operation outcome from a Virtual ligand alone.

  Cover hydrogen, lone pair, mixed cis/trans endpoints, fixed ligands that
  force a completion, compatible partial counts, ambiguous evidence, and
  impossible completion. Check exact frames and payloads, attribute
  preservation, and that an open coset alone does not imply ambiguity.
  Migrate direct consumers in S2d before requiring a green tree; introduce
  no temporary public interface to bridge the signature change.

  **Completed — 2026-11-08.** derive now returns the direct Solution, with explicit
  determinacy from H/LP identification. Unique completions identify the kinds;
  ambiguity or missing evidence retains the original row in Underdetermined;
  impossible completions produce the existing finding. The partial payload
  retains completed rows and independently established failures/mismatches.
  Unresolved frames do not supply configuration-mismatch findings. Stored frame
  order, cosets, entity constraints, and atom counts are preserved. No public
  completion type, error channel, or compatibility interface was added.

  Added 49 parameterized cases for finite/ranged/symbolic counts, fixed ligands,
  mixed cis/trans endpoints, open cosets, unresolved typed frames without Virtual,
  and partial payloads containing established findings. The 269 perceiver tests
  passed in an isolated copy with temporary consumer adaptations solely for
  compilation; the checkout's consumers remain unchanged. Library Clippy with
  warnings denied passed in that copy, including removal of the earlier
  type_complexity lint. Formatting and `git diff --check` passed.

  At S2c completion the checkout intentionally remained red: `cargo check -p
  umol-graph --lib` reported eight field-access errors in the unchanged stereo
  resolver and validator. S2d restores those consumers without a temporary bridge.

- **S2d — Consume derivation outcomes in resolution and validation** — breaking (red→green).
  Modules: umol-graph/src/ops/resolve/stereo.rs, umol-graph/src/ops/resolve.rs,
  umol-graph/src/ops/validate/stereo.rs, umol-graph/src/ingest.rs,
  umol-graph/benches/resolve.rs, umol-graph/benches/transform.rs,
  umol-graph/benches/transform/stereo.rs,
  umol-py/src/error.rs, and the owning tests. [dep: S2c]

  Apply the existing ligand-replacement edits only from a determined plan.
  Check findings before propagating underdetermination, retain unconditional
  mismatch rejection, and use failure policy only for impossible descriptions.
  Remove StereoResolveError and the outer Result from the stereo resolver
  methods; migrate the aggregate resolver, direct callers, and exception
  adapters. Validation reports unresolved identification as Underdetermined,
  rather than as an unrealizable frame, without mutating the molecule.

  Test all recovery policies, mixed failed and unresolved sites, empty
  underdetermined plans, no partial publication, consuming/in-place agreement,
  idempotence, and composite-resolver rollback after earlier phases. Keep
  existing stereogenicity/topicity checks and transformer preservation laws.

  **Completed — 2026-11-08.** Resolution and validation consume the explicit
  derivation outcome. Established rejections take precedence over ambiguity;
  failure policies apply only to impossible descriptions. An underdetermined
  plan is empty and publishes no partial replacements or removals. Removed
  StereoResolveError, its re-export, ResolveError::Stereo, and the resolution
  methods' outer Result. Migrated aggregate resolution and property callers.
  The existing Python exception adapter required no change. Benchmark preflight
  now checks the direct outcome and preservation of stored open cosets.

  Added 23 regression cases covering unresolved typed and Virtual frames,
  mixed failures/mismatches under all policies, unchanged rejected/unresolved
  inputs, and composite rollback after ligand replacement. Existing
  stereogenicity/topicity, idempotence, and transformer properties passed.
  Verification passed: graph library (3,711), property (223), and resolution
  conformance (4,219) tests; graph all-target Clippy with proptest/conformance
  and warnings denied; and 28 perception benchmark cases in test mode.
  On Python 3.13.16, 1,773 Rust binding tests passed with two ignored; after
  rebuilding the extension, 3,138 Python tests passed with two skipped.
  Formatting and `git diff --check` passed.

- **S2e — Merge constraint and entity failure configuration** — breaking (red→green).
  Modules: umol-graph/src/ops/resolve/stereo.rs, umol-graph/src/ops/resolve.rs,
  umol-py/src/resolve.rs, umol-py/tests/test_class_constructor_signatures.py,
  umol-py/tests/test_resolve.py, umol-py/tests/test_molecule.py,
  umol-py/tests/test_workflow.py,
  and docs/development/python-api.md. [dep: S2d]

  Leave only stereo_atom_failure and stereo_bond_failure, both defaulting to
  Error. Each covers constraint and entity realizability at its site kind.
  Migrate native callers and Python construction, getters, conversion, and
  repr together. Test the merged fields against both finding variants, reject
  obsolete Python keywords, remove redundant default updates from native config
  fixtures, and verify ambiguity is independent of policy.

  **Completed — 2026-11-08.** The two-field surface was already present in Rust
  and Python. Audited defaults, native callers, Python construction, getters,
  from_rust/to_rust, and repr against the merged contract. Both fields default
  to Error and govern constraint and entity failures at their own site kind.
  Documented that ambiguity stays underdetermined and mismatches always reject.

  Native constraint cases now distinguish atom and bond policies; unresolved
  resolution is checked under all nine policy combinations. Added 72 Python
  operation cases covering the four failure variants, all nine combinations,
  and both resolve and resolve_into, with exact diagnostics and molecules.
  Added rejection cases for both obsolete constructor keywords. Existing
  signature, getter, conversion, repr, and workflow coverage passed.

  Verification passed: graph library (3,717), property (223), and resolution
  conformance (4,219) tests; 1,773 Rust binding tests with two ignored; and
  all-target graph/binding Clippy with warnings denied. After rebuilding on
  Python 3.13.16, all 3,212 Python tests passed with two skipped. Formatting
  and `git diff --check` passed.

- **S2f — Raise SMILES and MOL stereo into explicit entities** — breaking (red→green).
  Modules: umol-io/src/table_ir/raise.rs,
  umol-io/src/table_ir/raise/utils.rs,
  umol-io/src/ctfile/parser.rs,
  umol-graph/src/ingest.rs,
  umol-io/tests/smiles_property.rs,
  umol-py/tests/test_molecule.py and umol-py/tests/test_import.py. [dep: S2e]

  Replace the wedge-derived tetrahedral and cis/trans constraint outputs with
  complete stereo entities. Transport source order and reference choices into
  the stored frame. Retain source-identified hydrogen/lone-pair ligands and use
  Virtual where the source does not identify the kind. Keep raise mechanical;
  invalid or impossible frames fail at the conversion boundary.

  Check exact imported entities for tetrahedral examples, alkenes, azo compounds,
  and imines, including source permutations and Either/unspecified stereo.
  Update import expectations in the owning tests, preserving the distinction
  between absent stereo and an entity with an unspecified configuration.

  **Completed — 2026-11-08.** TableIR raise now stores wedge-derived tetrahedral
  stereo and bond configurations on explicit entities. Stored frames account for
  source winding, endpoint references, and localized-bond reindexing when dative
  or noncovalent rows are separated. Explicit TableIR hydrogen/lone-pair ligands
  retain their kinds; unnamed positions in wedges and bond-reference frames use
  Virtual. Raise makes no valence-based identification. Either retains an entity
  with an undetermined coset; absent stereo produces no entity.

  Both definite and indefinite bond frames undergo the complete-frame checks.
  Unrepresentable frames and missing/incomplete wedge geometry report existing
  RaiseError variants. Updated MOL string/byte entry-point documentation and
  expected import carriers. Retained the source-permutation and reference-swap
  laws; the property generators and independent parity expectations are unchanged.
  Native and Python cases check Virtual before resolution and H/LP afterward for
  alkenes, azo compounds, and both imine orientations.

  Verification passed: 301 raise tests; all I/O library/test targets with
  proptest/conformance; 4,322 I/O library tests with depiction; graph library
  (3,721), property (223), and resolution conformance (4,219) tests; 1,773 Rust
  binding tests with two ignored; and all-target I/O, graph, and binding Clippy
  with warnings denied. After rebuilding on Python 3.13.16, 3,228 Python tests
  passed with two skipped. Nine ingest benchmark cases passed in test mode.
  Formatting and `git diff --check` passed. Revised ring-size operation alignment
  is covered by S2g–S2i.

- **S2g — Align model parameters, algorithm settings, and the perceiver** — breaking (red→green).
  Modules: umol-graph/src/ops/model.rs, umol-graph/src/ops/stereo.rs,
  umol-graph/src/ops/stereo/recognize.rs, umol-graph/src/ops/resolve/stereo.rs,
  umol-graph/src/ops/transform/stereo.rs, umol-graph/src/ops/validate/stereo.rs,
  umol-graph-ir/src/ir/symmetry.rs, umol-py/src/model/stereo.rs,
  umol-py/src/resolve.rs, umol-py/src/transform.rs, and their tests and benchmark
  callers. [dep: S2e, S2f]

  Keep stereo_bond_min_ring_size on StereoModel. Pass RingConfig explicitly from
  adder, resolver, and validator. Expose para_stereo_max_iterations: u32 on the
  adder, remover, enumerator, and validator configs; pass it through recognition,
  classification, enumeration checks, and validation. Use u32 for the graph-IR
  limit and counter. Migrate Python construction, properties, conversions, repr,
  and signature tests together.

  Share the model-based creation rule between derive and recognize. Preserve
  existing stereo, retained #C requests, structural checks, Virtual-kind
  identification, mismatch handling, and failure policies. Retire the skipped-site
  query, result field, and resolver deletion path. Test threshold boundaries,
  zero, acyclic sites, large rings, model snapshots, nondefault iteration limits,
  and u32 boundary handling. Preserve property laws and input domains.

  **Completed — 2026-11-08.** stereo_bond_min_ring_size remains on StereoModel;
  the perceiver reads its captured model. RingConfig is explicit in adder,
  resolver, and validator queries. All stereo symmetry paths forward
  para_stereo_max_iterations: u32, including both enumeration paths. The graph-IR
  limit and counter are u32. Rust/Python construction, getters, conversion, and
  repr are aligned. Existing descriptions and filtered assertions are preserved;
  the skipped-site deletion channel is removed. Property generators and laws
  are unchanged. The two conformance snapshots retain filtered source assertions.

  Verification passed: graph-IR library (11,825) and property (545) tests;
  graph library (3,943), property (223), and resolution conformance (4,219)
  tests; Rust binding tests (1,791); and Python tests (3,250) after rebuilding
  on Python 3.13.16. Existing ignored/skipped tests remain excluded. All 435
  transform benchmark correctness checks, all-target Clippy for graph IR,
  graph, and bindings with warnings denied, formatting, and `git diff --check`
  passed. E-specific model parameters and conformance checks remain S2h.

- **S2h — Validate model ring-size limits** — breaking (red→green).
  Modules: umol-graph/src/ops/model.rs, umol-graph/src/ops/validate/stereo.rs,
  umol-py/src/model/stereo.rs, and their tests. [dep: S2g]

  Add stereo_bond_e_min_ring_size, defaulting to 8, to StereoModel and its Python
  surface. Use RingTooSmall and RingTooSmallForE, each carrying the topology
  BondId, ring_size: u32, and minimum_ring_size: u32.
  Validate supplied descriptions against both model limits using the configured
  ring algorithms. Determine E relative to each ring's ligand frame. Retain
  Solution outcomes, StereoConformanceContradiction, and the existing error
  channel; introduce no separate inconsistency type.

  Test boundaries for both model limits, unrestricted checking, E/Z,
  positive #C evidence, open cosets, atom/frame permutations, and multiple rings.
  A determinate violation takes precedence over unrelated ambiguity. Assert input
  preservation on every outcome. The validator config remains Rust-only; add no
  Python validator API.

  **Completed — 2026-11-08.** Added the model field, both diagnostics, and Python
  model parity. Added 75 validation cases, including all endpoint-block frame
  actions with atom renumbering, a bond that is Z in its smaller ring and E in
  its larger ring, thresholds beyond the default ring-search bound, and known
  violations alongside unresolved stereo. Every new case checks input preservation.

  Verification passed: graph library 4,019; graph properties 223; resolution
  conformance 4,219; Rust bindings 1,791; Python 3,260 on Python 3.13.16.
  Existing ignored/skipped cases remain excluded. All-target Clippy for graph
  and bindings with property/conformance features and warnings denied,
  formatting, and `git diff --check` passed.

- **S2i — Verify operation and representation alignment through callers** — verification (green).
  Modules: umol-graph/src/ops/transform/stereo/tests,
  umol-graph/tests/property/transform/stereo.rs,
  umol-graph/tests/property/resolve.rs, umol-graph/tests/resolution;
  tests in umol-graph/src/ingest.rs, umol-graph/src/export.rs,
  and umol-graph/src/parse.rs; umol-py/tests/test_model.py,
  umol-py/tests/test_resolve.py, umol-py/tests/test_molecule_stereo.py,
  umol-py/tests/test_molecule_stereo_enumeration.py,
  umol-py/tests/test_molecule.py, umol-py/tests/test_smiles_export.py,
  umol-py/tests/test_class_constructor_signatures.py;
  affected benchmarks and docs/development/python-api.md. [dep: S2g, S2h]

  Compare missing stereo, #C+, and explicit stereo under the same operation
  config. Recognition and derivation agree on creation eligibility; adder and
  resolver preserve their own assertion-discharge contracts. Existing stereo
  stays intact through addition/resolution, can be enumerated, and is removed
  only by an explicitly requested removal operation or existing failure policy,
  never by the creation cutoff. Preserve the property laws and test domains.

  Verify parsing/raising fidelity and that wrappers using resolution forward
  the ordinary resolver config. Check molecule/reaction import and export with
  explicit small-ring stereo and unchanged borrowed inputs. Rebuild Python,
  run native and binding suites plus property/conformance tests, and run the
  affected benchmark correctness checks. Update current docs to the implemented
  contracts; record actual outcomes under these subitems.

  **Completed — 2026-11-08.** Added 71 Rust and 12 Python cases. Combined
  operation checks compare missing stereo, #C+, and supplied configurations under
  the same model and algorithm settings. SMILES and MOL callers preserve E/Z
  descriptions through raising and resolution under all failure policies;
  molecule and mapped-reaction export roundtrips retain small-ring stereo and
  borrowed inputs. Python covers model thresholds and small-ring enumeration.
  Existing property laws and generated domains remain intact.

  Verification passed: workspace 26,105; graph-IR properties 545; graph properties
  223; resolution conformance 4,219; I/O with property, conformance, and depiction
  features 17,260; Python 3,272 after rebuilding the extension on Python 3.13.16.
  Existing ignored/skipped cases remain excluded. All 561 resolve/transform
  benchmark correctness checks, workspace and feature-enabled I/O Clippy with
  warnings denied, formatting, and `git diff --check` passed.

**S2 gate:** external imports carry configurations on explicit entities, Virtual
identification preserves those configurations, and the perceiver consistently uses
model ring-size limits and supplied algorithm configs. Defined stereo is preserved;
small-ring E validation reports without mutation. Rust/Python config
surfaces and all affected consumers pass their checks. The site-value type change
remains S3.
**Passed — 2026-11-08.**

### S3 — Remove configuration-bearing site constraints

S3a and S3b changed the value types and shared DSL grammar. S3c completed their
coupled consumer migration and restored a passing Rust/Python test baseline,
including the required tests, fixtures, generators, and callers.
S3d and S3e completed valence and aromaticity discharge independently, with
their own passing checks. S3f completed the combined workflow verification.

S3c closed the compilation break. Every subsequent subitem starts from a passing
baseline and must pass its own checks before completion. Molecule
entry storage, frame construction, entity id allocation, and explicit-record
rendering retain their existing design.

- **S3a — Reduce site values and their structural readings** — breaking (red→green).
  Modules: umol-graph-ir/src/ir/stereo.rs,
  umol-graph-ir/src/ir/constraint/atom.rs,
  umol-graph-ir/src/ir/constraint/bond.rs,
  umol-graph-ir/src/ir/view/atom.rs,
  umol-graph-ir/src/ir/view/bond.rs,
  umol-graph-ir/src/ir/view/constraints.rs,
  umol-graph-ir/src/ir/canonicalize.rs,
  umol-graph-ir/src/ir.rs,
  umol-py/src/stereo.rs,
  umol-py/src/constraint/atom.rs,
  umol-py/src/constraint/bond.rs,
  umol-py/src/lib.rs,
  umol-py/python/umol/__init__.py. [dep: S1a, S2i]

  Make Stereo a unit variant in the two site forms, use bool as their literal
  type, add From<bool>, and remove the ground enums, configuration conversions,
  and coset-specific comparisons. Update the Python variants, bool coercions,
  literal extraction, setters, registration, and exports in the same change.
  Derive the stereo constraint from a known kind without deriving a coset, and
  remove configuration data from site-constraint canonical keys.

  Test all three states and their lattice laws, bool extraction, open/complete
  readings, unknown kinds, undetermined cosets, per-kind negatives, and relabeling.
  Reduce the site-value generators while retaining the full entity-configuration
  generators and their independent references.

  **Implemented — 2026-11-08; verification pending S3.** Both site forms now
  have unit Undetermined, NotStereo, and Stereo variants, with bool literals and
  From<bool>. The ground enums and configuration-to-constraint conversions are
  removed. Structural readings use only the entity's known kind; canonical keys
  carry no coset. Python variants, coercions, getters/setters, and exports follow
  the same contract. Explicit-entity configuration types remain unchanged.

  Added finite-state Rust/Python tables and structural-reading cases for known
  and unknown kinds, open cosets, virtual ligands, and atom relabeling. Site
  generators use the reduced domain; full entity-configuration generators and
  reference properties are retained. Formatting, Python syntax checks, and
  git diff --check pass. Graph-IR runtime and property tests pass with S3b on
  2026-11-09. Binding verification and the Python rebuild await the remaining
  S3 migrations.

- **S3b — Use the same reduced stereo grammar in every DSL context** — breaking (red→green).
  Modules: umol-graph-ir/src/dsl/stereo.rs,
  umol-graph-ir/src/dsl/atom.rs,
  umol-graph-ir/src/dsl/bond.rs,
  umol-graph-ir/src/dsl/constraint.rs,
  umol-graph-ir/src/dsl/molecule.rs,
  umol-graph-ir/src/dsl/reaction.rs,
  umol-graph-ir/src/dsl/reaction_span.rs,
  umol-graph-ir/src/dsl/edit.rs,
  umol-graph-ir/src/ir/molecule/spec.rs,
  umol-graph-ir/src/macros.rs. [dep: S3a]

  Restrict #T/#C and their EDN values to the three states, including
  :stereo. Reject all configuration-bearing forms, including in molecule maps,
  reaction LHS maps, aliases, edits, deltas, and reaction-span sides. Retain the
  full configuration grammar on explicit stereo entities and the existing
  molecule parser, registration, and rendering paths.

  Test uniform rejection of literals, sets, variables, and operators in attribute
  positions. Check exact tree/streaming/render roundtrips of explicit stereo
  records, including symbolic configurations, defaults, aliases, entity keywords,
  and ids referenced by reaction deltas. Retain the explicit-entity reference
  checks that distinguish configurations.

  **Completed — 2026-11-09.** The shared compact grammar accepts only *, !,
  and +; EDN values are :undetermined, :not-stereo, and :stereo. Tree and streaming
  readers reject configuration payloads. The attribute-only coset DSL and its
  readers are removed; explicit stereo entities retain their full grammar.

  Tests cover the current grammar across standalone values, aliases, molecule
  maps, reaction LHS values and deltas, edits, and both reaction-span sides.
  Exact roundtrips preserve explicit configurations, symbolic terms, virtual
  ligands, defaults, keywords, aliases, and stereo ids referenced by deltas.
  Parser seeds use the current constraint values. The unused coset comparison
  helper is removed; its configuration-matching cases remain on
  StereoConfigurationForm. The last removed-ground-enum fixture and unused
  property imports are updated.

  Validation passed: cargo test -p umol-graph-ir --features proptest, including
  11,958 library tests and 545 property tests, plus the integration, macro, and
  documentation tests. Clippy passed for graph-IR with proptest, all targets,
  and warnings denied; formatting and git diff --check passed. The coordinated
  consumer migration and its first passing Rust/Python gate are S3c.
  The umol-graph check stops in those consumers. S3f is the later combined gate.

- **S3c — Finish #T/#C consumer migration and restore passing test suites** — breaking (red→green).
  Modules: umol-graph/src/ops/stereo.rs,
  umol-graph/src/ops/stereo/recognize.rs,
  umol-graph/src/ops/resolve/stereo.rs,
  umol-graph/src/ops/validate/stereo.rs,
  umol-graph/src/ops/transform/stereo.rs,
  umol-graph/src/ops/transform/stereo/enumeration.rs,
  umol-graph/src/ops/transform/hydrogens.rs,
  umol-graph/src/ops/resolve.rs, umol-graph/src/export.rs, and their tests;
  umol-graph-ir/tests/property/strategies.rs,
  umol-graph/tests/property, umol-graph/tests/resolution/data,
  umol-graph/tests/kekulization/data, umol-py/src/molecule.rs,
  umol-py/src/reaction.rs, umol-py/tests, and the Rust/Python, benchmark,
  and fuzz consumers listed in the change inventory. [dep: S3a, S3b, S2i]

  Complete the coupled migration of the removed configuration-bearing attribute
  representation. This single breaking subitem owns stereo derivation and
  resolution, conformance validation, recognition and transformations, standalone
  projection, aggregate projection, and direct export. Their required migrations
  must finish together because every graph test compiles the whole graph library.

  Implement the recorded [derivation and resolution](#perception-and-resolution),
  [transformer contracts](#resolution-and-stereo-transformers),
  [conformance validation](#conformance-validation), and
  [projection and export](#standalone-projection-and-export) semantics. Use only
  the public interface changes specified under [Operation interfaces](#operation-interfaces).
  Valence and aromaticity discharge changes remain separate in S3d and S3e.

  Migrate every Rust test body, fixture, generator, property reference, benchmark
  caller, fuzz input, and Python consumer needed by these changes in this subitem.
  Move configuration descriptions to explicit stereo records while preserving
  intended chemistry, property laws, and measured workloads. Each changed
  operation owns its exact success, failure, unchanged-input, and partial-input
  cases. Cover frame/configuration preservation, virtual-ligand identification,
  scalar constraint handling, selection, enumeration order and ownership,
  hydrogen reconstruction, projection loss, and molecule/reaction export.
  Check standalone projection's idempotence and negative-constraint conflicts,
  and Convey's exact frame-preserving roundtrips and unchanged borrowed inputs
  on rejection.

  Restore a passing library build, Rust unit and integration suites, affected
  feature-gated property and conformance suites, and benchmark compilation.
  Rebuild the extension with the repository Python 3.13 environment and pass the
  Python suites. Use the verification commands below. Record results here;
  compilation-only success is insufficient. This subitem closes only once the
  representation migration has a passing test baseline, without compatibility
  shims or disabled coverage.

  **Completed — 2026-11-09.** Stereo derivation and recognition use kind-only
  constraints; new entities have undetermined cosets. Resolution preserves stored
  configurations and discharges independently entailed constraints beside both new
  and retained entities. Validation checks explicit frame/configuration prerequisites
  without requiring a scalar marker. Enumeration uses entity configuration domains,
  and hydrogen conversion preserves scalar constraints while replacing explicit
  hydrogen ligands in their stored slots.

  Standalone stereo projection is model-independent and deliberately discards
  configurations, frames, and entity-local constraints. Aggregate projection retains
  explicit stereo; its stereo flag and phase errors are removed. Convey lowers those
  entities directly, preserving configurations under their stored frames and rejecting
  unsupported data. Rust/Python tests, generators, fixtures, and benchmarks use the
  current representation. Two conformance snapshots now retain specified small-ring
  stereo as explicit entities; the remaining resolution outputs are unchanged.

  Verification passed: 26,214 workspace tests (nine ignored), 545 graph-IR
  properties, 223 graph properties, 4,219 resolution cases, and the I/O
  property/conformance/depiction suites. The rebuilt Python 3.13 extension passed
  3,307 tests (two skipped). Resolve/transform benchmarks passed 561 correctness
  checks; all 255 curated graph-IR parser seeds replayed successfully. Workspace
  and feature-enabled graph/I/O Clippy passed with warnings denied; formatting
  and git diff --check passed. Valence and aromaticity discharge remain S3d/S3e.

- **S3d — Discharge entailed #v, #d, and #t constraints** — behavior change (green).
  Modules: umol-graph/src/ops/resolve/valence.rs,
  umol-graph/src/ops/resolve.rs, and their tests;
  umol-graph/tests/property/resolve.rs, umol-py/tests/test_resolve.py,
  umol-graph/benches/resolve.rs. [dep: S3c]

  Retain candidate admission and completion selection. In
  Resolver::plan_constitution, discharge #v, #d, and #t only when their
  independent open derived readings entail them. This extends the existing
  #v discharge to redundant donated- and accepted-pair constraints, rather
  than retaining those keys unconditionally. Keep compatible constraints that
  still add information and preserve rejection of definite conflicts.

  Verify both valence strategies, entailed #v with and without field changes,
  entailed #d/#t on explicit dative incidence, and retained zero-pair assertions
  when no dative relation supplies evidence. Assert complete output molecules,
  exact conflict outcomes, idempotence, and unchanged in-place inputs on
  rejection. Preserve the existing candidate-selection laws.

  Update affected fixture, Python, and benchmark expectations here, alongside
  the behavior change. Pass the owning Rust unit cases and affected resolution
  properties/conformance cases before closing this subitem; rebuild the extension
  before any affected Python checks. These checks run against the passing S3c
  baseline and do not wait for S3f.

  **Completed — 2026-11-09.** Resolver::plan_constitution uses the independent
  open derived reading to discharge #v, #d, and #t on all atoms, including
  already concrete atoms without a completion. Compatible constraints that
  still add information remain; #d0/#t0 are retained without dative incidence.
  ValenceResolver::admit and completion selection are unchanged.

  Exact Rust/Python outputs cover both valence sources, field completion and
  concrete inputs, literal and set-valued assertions, dative donor/acceptor
  readings, idempotence, and unchanged in-place inputs on conflicts. The
  aromaticity property retains its standalone #v expectation and separately
  checks aggregate discharge. Conformance snapshots and benchmark expectations
  required no changes.

  Verification passed: 26,262 workspace tests (nine ignored), 224 graph
  properties, 4,219 resolution conformance cases, and 126 resolve benchmark
  correctness checks. The rebuilt Python 3.13 extension passed 3,387 tests
  (two skipped). Workspace and feature-enabled graph Clippy passed with warnings
  denied; formatting and git diff --check passed.

- **S3e — Discharge aromaticity constraints beside new and retained systems** — behavior change (green).
  Modules: umol-graph/src/ops/resolve/aromaticity.rs,
  umol-graph/src/ops/aromaticity.rs,
  umol-graph/src/ops/resolve.rs, and their tests;
  umol-graph/src/ops/transform/aromaticity.rs;
  umol-graph/tests/property/resolve.rs, umol-py/tests/test_resolve.py,
  umol-graph/benches/resolve.rs. [dep: S3c]

  Apply independent entailment to atom aromatic-valence and bond-aromaticity
  assertions against the planned aromatic systems. Discharge redundant
  assertions beside retained systems as well as new systems. Stop treating
  absence of a selected system as sufficient to discharge a negative assertion.
  Apply the same rule in standalone aromaticity planning and aggregate
  Resolver::plan_constitution. Preserve selection and mismatch rejection under
  both failure policies. Preserve Aromatizer's insertion and assertion-removal
  contract when adjusting the shared AromaticityPerceiver::plan_aromatic_system
  method.

  Verify atom contributions and bond assertions with new and existing systems,
  retained negatives outside systems, positive requests retained under Keep
  when unrealizable, and exact mismatch contradictions under Error and Keep.
  Check full outputs, idempotence, standalone/aggregate agreement, and unchanged
  inputs on rejected in-place resolution. Preserve the aromaticity search laws.

  Update affected fixture, Python, and benchmark expectations here, including
  coverage of the shared Aromatizer caller. Pass the owning Rust unit cases and
  affected resolution properties/conformance cases before closing this subitem;
  rebuild the extension before any affected Python checks. This subitem is
  independently testable after S3c and does not depend on S3d's behavior change.

  **Completed — 2026-11-09.** Standalone and aggregate aromaticity resolution
  discharge independently entailed atom and bond assertions beside new and
  retained systems. AromaticityResolver::plan_discharge also covers retained
  atom assertions; AromaticityPerceiver::plan_aromatic_system checks the new
  system's contributions before discharging its member assertions. More specific
  assertions remain when stored contributions are undetermined. Negative
  assertions outside systems remain. Selection and mismatch rejection under
  Error and Keep are unchanged, as is Aromatizer's insertion and assertion removal.

  Exact tests cover stored participant order, literal and set-valued assertions,
  undetermined stored contributions, failed requests under Keep, idempotence,
  standalone/aggregate agreement, and unchanged rejected inputs. Mixed-system
  properties retain their generated domains and check the updated discharge
  contract. Six conformance snapshots now retain #a! on bonds outside systems;
  the ingestion expectation for explicitly linked biphenyl changes likewise.
  Benchmark expectations required no changes.

  Verification passed: 26,286 workspace tests (nine ignored), 224 graph
  properties, 4,219 resolution conformance cases, and 561 resolve/transform
  benchmark correctness checks. The rebuilt Python 3.13 extension passed 3,427
  tests (two skipped). Workspace and feature-enabled graph Clippy passed with
  warnings denied; formatting and git diff --check passed.

- **S3f — Verify combined resolution, stereo transformations, and export** — verification (green).
  Modules: umol-graph/tests/property/resolve.rs,
  umol-graph/tests/property/project.rs,
  umol-graph/tests/property/transform, umol-graph/tests/resolution;
  tests in umol-graph/src/ingest.rs and umol-graph/src/export.rs;
  umol-io property/conformance suites, umol-py/tests, and the affected benchmark
  targets. [dep: S3c, S3d, S3e]

  Verify the combined public workflows after the independently tested changes:
  resolution followed by validation, stereo addition/removal/enumeration, hydrogen
  reconstruction, standalone projection, and molecule/reaction import and export.
  Check exact stored results, unchanged rejected inputs, and the supported
  Rust/Python interfaces.

  Run the workspace and affected feature suites, rebuild and test Python, and
  check the affected benchmark targets. Reconcile the public-symbol and consumer
  inventory with completed S3c–S3e work. This is an integration gate; required
  code, fixture, and consumer migrations belong to their owning implementation
  subitem and are not deferred here.

  **Completed — 2026-11-09.** The combined checks cover resolution, validation,
  stereo addition/removal/enumeration, hydrogen conversion, standalone and
  aggregate projection, and molecule/reaction import and export. The existing
  generated mixed aromatic/stereo SMILES roundtrip now also calls Validator
  before export and after re-import, using the same chemistry model. Its input
  domain and roundtrip assertions are unchanged.

  The public-symbol and consumer inventory agrees with the completed migration:
  scalar stereo forms use bool literals; configurations remain on explicit
  entities; Virtual retains exact ligand identity; all DSL contexts use the
  reduced attribute grammar. Standalone stereo projection has its model-free
  associated functions; aggregate projection and export retain explicit stereo.
  Model ring thresholds, operation configurations, shared contradictions, and
  the supported Python exports match the recorded interfaces. No further
  production migration was needed for this gate.

  Verification passed: 26,286 workspace tests (nine ignored), 545 graph-IR
  properties (one ignored molecule-scope matching case), 224 graph properties,
  4,219 resolution conformance cases, and 17,260 feature-enabled I/O tests,
  including property, conformance, layout, and depiction coverage. Resolve and
  transform benchmarks passed 561 correctness checks; I/O benchmarks passed
  247. The rebuilt Python 3.13 extension passed 3,427 tests (two skipped).
  Workspace and feature-enabled graph/I/O Clippy passed with warnings denied;
  formatting and git diff --check passed. No snapshots changed in this subitem.

**S3 completed — 2026-11-09.**

**S3 gate:** site constraints contain only frame-independent information, all DSL
contexts reject configuration-bearing #T/#C, and configuration descriptions use
explicit records. Chemistry operations and export obtain configuration from those
entities. All affected callers compile and pass their tests; ligand identity and
matching remain exact. Valence, aromaticity, and stereo resolution follow the
shared conflict/discharge policy.
Conformance validation checks explicit stereo descriptions independently of
site assertions and preserves the tier-2 open/complete reading choice.

### S4 — Documentation and final verification

- **S4a — Document reduced #T/#C constraints and explicit stereo entities** — additive (green).
  Modules: umol-graph-ir/spec/umol-dsl-spec.md,
  docs/chemistry/representation.md,
  docs/development/data-types.md,
  docs/development/python-api.md,
  docs/development/nomenclature.md,
  .agents/skills/mol-notation/SKILL.md. [dep: S3f]

  Document the settled representation, construction, syntax, conversion, and
  binding contracts and verify implementation compliance. Record an implementation
  discrepancy under its owning S3 subitem and correct it there. Correct affected
  examples and wording that equates a stereo description with stereogenicity.
  Preserve shared-skill ownership and review the complete source/documentation
  diff against the inventory. Source comments and rustdoc describe the implemented
  contract without referring to discussion records.

  **Completed — 2026-11-09.** The DSL specification, representation and development
  guides, and canonical notation skill now distinguish scalar kind constraints
  from explicit framed configurations. They document :virt, exact ligand identity,
  construction and conversion boundaries, bool literals and Python inputs, and
  stereo descriptions independently of stereogenicity. The notation skill uses
  the current DSL macros; its .claude redirect retains matching metadata.

  Corrected stale rustdoc in graph-ir's stereo value and DSL modules and graph's
  stereo transformers, including claims that #T/#C carried frames or cosets.
  The inventory review found no implementation discrepancy requiring an S3 change.
  No executable code or tests changed.

  Verification passed: cargo doc for umol-graph-ir and umol-graph, Python 3.13
  checks of both new examples and the documented syntax/input distinctions,
  notation-skill validation and relative links, formatting, and git diff --check.

- **S4b — Document resolver conflict and discharge semantics** — additive (green).
  Modules: docs/development/operations.md,
  docs/development/README.md, docs/development/data-types.md,
  docs/development/nomenclature.md, and AGENTS.md. [dep: S4a]

  Extend the operations guide with the conflict/discharge policy for valence,
  aromaticity, and stereo, including refinement, independent structural
  entailment, open readings, and the separation of conflicts from
  realizability-failure recovery.

  Consolidate resolver conflict, discharge, and recovery-policy explanations
  from nomenclature.md's Discharge, Failure, Inconsistency, Mismatch, Policy,
  Recovery action, and Resolution entries. Retain the term definitions there
  and link to the operations guide for the rules. Review data-types.md's
  Validation and Resolution, projection, and transformations sections for
  overlapping operation-policy prose; replace duplication with links while
  preserving its construction, ownership, and publication contracts.

  Keep integrity.md separate. Preserve the recorded tier-3 conformance contract,
  keep guide links synchronized, and check implementation against the recorded
  contracts. Broader operation-documentation work is inventoried in S4c.

  **Completed — 2026-11-09.** The operations guide defines refinement, unconditional
  conflict rejection, independent structural entailment, open readings for discharge,
  and realizability-failure recovery across valence, aromaticity, and stereo.
  It identifies each resolver's owned keys and the current failure-policy actions.
  Nomenclature retains the definitions and links to those rules, replacing the old
  creation-dependent discharge descriptions and removed mismatch-recovery actions.

  Data type contracts links to the shared policy and retains construction,
  ownership, planning, and publication guarantees. The development-guide index is
  updated; AGENTS.md already registers the operations guide. The existing tier-3
  conformance section is unchanged, and integrity remains separate. Implementation
  review found no discrepancy requiring a code change.

  Verification passed: 1,159 resolver unit tests with
  cargo test -p umol-graph --lib ops::resolve, guide-link checks, formatting, and
  git diff --check. No executable code or tests changed.

- **S4c — Inventory operations needing explicit semantic documentation** — additive (green).
  Module: discussion/243-stereo-site-assertions-2026-11-07.md. [dep: S4b]

  Review public operation families against data-types.md, nomenclature.md,
  algorithms.md, and the operations guide. Cover perception and derivation,
  resolution, projection and boundary conversion, validation, transformers,
  normalization, canonicalization, frame transport, matching, and rewriting.
  Enumerate the operations whose intended contracts need explicit documentation
  or consolidation.

  Record a table here with the owning API, a summary of its intended contract,
  existing references, missing documentation, and any confirmed implementation
  discrepancies. Include effects on attributes, entities, and constraints,
  partial-input behavior, failure outcomes, and preservation obligations where
  relevant. Identify the documentation or implementation follow-up for each gap.
  This subitem ends with the inventory; writing the remaining guide sections
  and implementing fixes are separately tracked work. Existing semantics remain
  the basis for those follow-ups.

- **S4d — Run final checks and compare measurements** — verification (green).
  Modules: the test and benchmark targets listed below, and
  discussion/243-stereo-site-assertions-2026-11-07.md,
  discussion/000-status.md. [dep: S4c]

  Run the completion commands and replay the curated parser inputs. Repeat the
  selected S0 benchmarks with the same inputs and timing settings. Compare fixed
  explicit-entity workloads directly and report changed parser/import workloads
  separately. Investigate repeatable regressions outside baseline variation.

  Record commands, outcomes, and measurements beneath the owning subitems.
  Reconcile the public-symbol inventory, confirm the scoped implementation and
  supported bindings agree, and mark the document complete only after all
  required work and checks are finished.

### Verification commands

Use the owning modules' tests during each subitem. Every stage ends with the
workspace suite and its affected property, feature, and Python checks passing.
The final gate runs these commands from the repository root in one activated
shell; separate invocations repeat activation before compiling or loading Python:

```sh
source umol-py/.venv/bin/activate
python --version
cargo +nightly fmt --all -- --check
cargo test --workspace
cargo test -p umol-graph-ir --features proptest --test property
cargo test -p umol-graph --features proptest,conformance --test property --test resolution
cargo test -p umol-io --features proptest,conformance,depiction
cargo clippy --workspace --all-targets -- -D warnings
cargo clippy -p umol-io --features proptest,conformance,depiction --all-targets -- -D warnings
maturin develop --locked --manifest-path umol-py/Cargo.toml --features depiction
pytest -q umol-py/tests
git diff --check
```

Confirm that Python is the repository environment's Python 3.13. The feature
commands cover IO conformance, property tests, and native layout/depiction; the
all-targets checks compile benchmark callers. Replay the existing graph-IR entity,
molecule, reaction, span, and edit fuzz targets on their curated seed inputs.

S0–S3 are complete. S3c restored the passing baseline; S3d and S3e passed their
independent gates, followed by S3f's combined workflow checks. S4a–S4b are complete;
S4c is next.
S4d reuses applicable recorded S3f results while source and test inputs remain
unchanged, and runs the remaining final checks and
measurement comparisons.
There are no deferred stages in this deliverable. The exclusions in Scope remain
unchanged. S4c–S4d remain pending.
