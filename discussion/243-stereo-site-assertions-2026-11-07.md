# 243 — Frame-independent atom and bond stereo assertions

Status: In Progress
Date: 2026-11-07
Relates: [214](214-aggregate-frame-semantics-2026-08-28.md),
[226](226-smiles-roundtrip-design-2026-09-10.md),
[240](240-validation-constraint-ops-2026-10-18.md),
[242](242-stereo-transformers-2026-10-29.md),
[data-type guide](../docs/development/data-types.md)

## Scope

Restrict stored atom- and bond-level stereo assertions to frame-independent stereo constraints.
Configuration data is stored only on explicit stereo atoms and stereo bonds, together
with their ligand frames. A generic virtual stereo ligand lets those entities hold a
frame whose virtual ligand is not yet identified. SMILES/MOL raising populates the
entities; export conveys them directly to the format boundary. The DSL uses
explicit stereo records for configuration data. Configuration-bearing #T/#C
syntax is rejected in every DSL context.

This records the selected design direction and its consequences. S0, S1, and S2a–S2e
are complete; the remaining S2–S4 work is pending. The semantics, value interface, operation
interfaces, and names below are the working design. S2 records the revised
perception/resolution proposal; it is provisional, not final design approval.
Implementation is paused pending explicit instruction. The staged implementation
plan follows the change inventory.

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
stereogenicity. Presence/kind assertions and derived stereogenicity retain their
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
matches_value(u32) are removed; ordinary lattice matching compares the presence
assertions.

TetrahedralConfiguration (Ccw/Cw) and CisTransConfiguration (Z/E) remain
configuration types for explicit stereo entities. Their conversions into the
atom/bond constraint forms are removed. A configuration is not silently converted
to a positive presence assertion.

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
| Presence or kind | Stereo constraints |
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

#T and #C accept only the three presence states: unconstrained, negative, and
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
to explicit records; use presence markers where presence is what a test asserts.

## Raising SMILES and MOL data

Raise must construct stereo atoms/bonds with their ligand frames and configuration
data, instead of placing configuration indices in #T/#C for later materialization.
This explicitly reverses the earlier instruction to use those constraints as an
import channel.

The current [TableIR raise](../umol-io/src/table_ir/raise.rs) already maps explicit
tetrahedral table frames directly into graph-IR stereo atoms. Wedge-derived
tetrahedral configurations and cis/trans bond frames still become #T/#C payloads;
those paths need the same explicit-entity outcome.

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

## Resolution and stereo transformers

StereoResolver remains useful for evidence-driven resolution of stereo constraints. A positive
stereo constraint can cause it to construct a stereo entity of the asserted kind with an
undetermined coset. Imported configurations already exist on entities after raise,
and StereoResolver replaces their Virtual ligands with the identified kinds.
Shared perception and consistency checks must read stereo constraints as
presence/kind assertions and cease extracting configurations from them.

The core contracts of the transformers in doc 242 survive:

- StereoAdder recognizes missing selected entities. New entities have undetermined
  cosets; existing entities, including imported configurations, remain intact.
  Positive stereo constraints supply kind evidence and negative stereo constraints still exclude sites.
- StereoRemover removes selected entities and clears corresponding stereo constraints at their
  sites. Its Undetermined mode continues to inspect the explicit entity's coset,
  so that mode retains its meaning.
- StereoEnumerator obtains its configuration domains from explicit entities.
  Its assignment traversal does not need redesign. Checks of configuration domains
  stored in #T/#C disappear; stereo constraint checks and constraints on explicit stereo
  entities remain applicable.

StereoResolveConfig merges its constraint and entity failure policies, leaving
stereo_atom_failure and stereo_bond_failure. Each covers both a failed stereo constraint and a
failed entity at its site kind: Keep retains the failed value, and Remove clears the
stereo constraint or removes the entity. StereoInconsistency keeps its separate variants. At
ingest every entity comes from the import, so a caller distinguishes import failures
from stored failures by the config passed to each call.

Mismatches between independently realizable constraints and entities remain
contradictory under every failure policy. Keep and Remove apply to realizability
failures, not to choosing which of two contradictory statements to believe.

Fixtures that convey configurations through #T0/#T1 or #C0/#C1 must be
rewritten using explicit stereo entities. Presence-only cases use the reduced
stereo constraints.

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

Hydrogen-omission and representability checks in Convey must inspect explicit stereo
entities instead of relying on configuration-bearing #T/#C. Molecular and reaction
SMILES export share this boundary conversion change. The projection/export design
in doc 226 remains historical evidence for the current implementation; its stereo
intermediate is the part this proposal replaces.

## Minimum ring size

stereo_bond_minimum_ring_size applies at external-format I/O only. Its purpose is
output: without it, every ring double bond would be written with stereo notation.
StereoResolver and StereoPerceiver no longer apply it.

Resolved external-format ingest applies the rule after raise and before resolution, where the
chemistry model is available. Stereo bonds and positive #C constraints in rings below the
minimum are dropped, so such a bond cannot fail during resolution. Reaction SMILES
ingest shares this path. Convey omits the same sites when writing TableIR. This is an
exception to the Convey contract that a definite configuration is never downgraded.

The umol DSL is not an external format. An explicit stereo-bond record in a
ring below the minimum keeps its entity through resolution. Unresolved ingest has no
chemistry model, so it retains the stereo bond, and later resolution retains it too.
Convey still omits it.

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
stereo constraints. Recognition and classification retain their existing Rust
query contracts and StereoPerceptionError. The coset_for readings retain exact
frame transport. Remove the Python perception class and StereoDerivation;
supported Python operations retain their configuration and diagnostic values.

Remove StereoPerceiver::skipped_stereo_bonds and
StereoDerivation::skipped_stereo_bonds in S2g. They describe the ring-size
suppression that now belongs exclusively to I/O. Other
ring-dependent recognition rules and their RingConfig argument remain.

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
negative presence constraint uses the existing assertion error.

Remove ProjectFlags::STEREO and ProjectError::Stereo from aggregate projection.
ProjectContradiction::Stereo also disappears; no remaining projection stage produces it.
Resolver::project and its other projection stages keep their interfaces and
preserve stereo entities. Convey keeps its existing public arguments because
the other projection stages and the I/O ring policy still require them.

## Other affected surfaces

| Surface | Implication |
| --- | --- |
| Rust stereo site forms | Keep the named three-state forms, use bool literals and construction, and remove the ground enums and coset-specific methods/conversions. Keep configuration types for explicit stereo entities. |
| Stereo ligands | Add the Virtual variant and :virt spelling. Update exhaustive kind dispatch, canonical keys, and Python enum conversion. Keep ligand identity and matching exact. |
| Atom/bond constraint readings | Report presence/kind without reconstructing a neighbor-relative coset. A known stereo kind can give a positive stereo constraint even when its coset is undetermined. Existing open/complete reading distinctions still apply. A caller that needs a coset uses coset_for with a stated frame. |
| Validation | Check stereo constraint realizability and agreement with entity presence/kind. Configuration and ligand-frame validation continue on explicit entities. |
| Hydrogen folding/unfolding | Remove special transport of inline #T/#C cosets. Continue maintaining actual/virtual ligands on explicit stereo entities. Hydrogen folding currently has separate code for both representations. |
| Canonicalization | Remove coset data from atom/bond constraint keys. Those attributes then have their ordinary relabeling semantics; explicit stereo canonicalization still transports its frames and configurations. |
| Matching and gluing | Ordinary attribute operations compare reduced stereo constraints. Configuration-sensitive comparisons use explicit stereo entities and existing exact frame alignment and transport. |
| DSL parser | Reject configuration-bearing #T/#C in every context. Retain the current AtomDsl/BondDsl wrappers and molecule input entries; configuration syntax remains available in the explicit stereo-entity DSLs. |
| DSL renderer | Retain explicit stereo-record rendering. No implicit-frame inference, conditional inline output, or new serialization configuration. |
| Python bindings | Use bool or the corresponding form for constraint input, and bool or None for literal extraction. Remove ground enum bindings and configuration-valued constraint coercions. Mirror the merged failure policies and uniform DSL rejection. |
| Fixtures and documentation | Rewrite configuration-bearing atom/bond strings using explicit stereo records. Presence-only cases use the reduced constraints. Update examples, properties, benchmarks, and API descriptions that assumed a stored coset. |
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
| Same | matches_value(u32) | Remove from both forms. Do not replace it with another inherent comparison method; ordinary lattice matching handles boolean presence. |
| Same | TetrahedralStereo, CisTransStereo | Remove both ground enums and all conversions to or from them. Remove their re-exports in umol-graph-ir/src/ir.rs. |
| Same | `From<TetrahedralConfiguration>` for TetrahedralStereoForm; `From<CisTransConfiguration>` for CisTransStereoForm | Remove. Configuration input must not implicitly become a presence constraint. |
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
| umol-graph/src/ops/resolve/stereo.rs | StereoResolveConfig | Remove tetrahedral_stereo_failure and cis_trans_stereo_failure; retain stereo_atom_failure and stereo_bond_failure with the merged meanings. Keep Default = Error. |
| Same | StereoFailurePolicy; StereoResolveContradiction; StereoResolveError | Retain the policy. Rename the contradiction wrapper to shared StereoContradiction in ops/stereo.rs; retain its StereoInconsistency payload. Remove the empty error type. |
| Same | StereoResolver::new, with_config, plan_resolve, resolve, resolve_into | Retain construction; resolution returns Solution directly, without the outer Result. Add missing requested entities and replace Virtual kinds on retained entities. |
| Same | StereoResolver::plan_project, project, project_into | Keep names; remove the receiver and return the direct Result shapes recorded above. |
| Same | StereoProjectError | Retain UnsupportedStereoAtom, UnsupportedStereoBond, AtomAssertion and BondAssertion; remove StereoAtomFrame and StereoBondFrame. |
| umol-graph/src/ops/valence.rs and ops/resolve/valence.rs | ValenceContradiction | Move the definition to ops::valence, matching AromaticityContradiction and StereoContradiction in their domain modules. Preserve variants, conversions, and messages; migrate imports and remove the former resolver re-export. |
| umol-graph/src/ops/resolve.rs | ProjectFlags, ProjectError, ProjectContradiction | Remove STEREO, Stereo and Stereo respectively. Retain non-stereo flags and variants; aggregate projection preserves stereo entities. Remove the uninhabited ResolveError::Stereo; retain ResolveContradiction::Stereo with the renamed payload. |
| umol-graph/src/ops/stereo.rs | StereoPerception; derive; StereoDerivation::stereo_atoms / stereo_bonds | Rename the actor to StereoPerceiver. derive returns Solution explicitly; retain row types for completed stored frames and requested additions. |
| Same | StereoPerceiver::skipped_stereo_bonds; StereoDerivation::skipped_stereo_bonds | Remove both surfaces in S2g. |
| Same and ops/stereo/recognize.rs | StereoPerceiver::new, derive_stereo_atom, derive_stereo_bond, recognize, classify; StereoPerceptionError and StereoInconsistency | Retain signatures and existing support boundaries. New recognized entities have undetermined cosets; classify still reads stored arrangements. |
| umol-graph/src/ops/transform/stereo.rs and ops/transform.rs | StereoAdder, StereoRemover, StereoEnumerator, their configs/errors/iterator; Transform convenience methods | Retain public interfaces and shared StereoSelection. Remove source-constraint coset handling from the consumers. |
| umol-graph/src/ops/transform/hydrogens.rs | HydrogenFolder, HydrogenUnfolder and their errors/configs | Retain interfaces; remove inline-constraint coset transport and maintain explicit frames through existing edits. Virtual is not implicitly an H ligand. |
| umol-graph/src/ops/validate/stereo.rs | StereoConformanceValidator and result/error surface | Retain signatures; consume the revised derivation and presence constraints. |
| umol-graph/src/ops/model.rs | StereoModel::stereo_bond_minimum_ring_size | Retain field and default; update documentation and consumers to the I/O-only scope. Other stereo-model fields do not change. |
| umol-io/src/table_ir/raise.rs | `TryIntoIr<Molecule>` for &table_ir::Molecule; RaiseError | Retain public conversion shape and existing error vocabulary. Wedge and bond-frame paths produce explicit entities and report unrepresentable frames at this boundary. |
| umol-graph/src/ingest.rs | Interpret for Smiles / ReactionSmiles; ingest_smiles and ingest_reaction_smiles families | Retain entry points and errors, including byte and unresolved variants. Resolved routes apply the ring filter between raise and resolution; unresolved routes retain entities. |
| umol-io/src/ctfile/parser.rs | parse_mol_to_ir, parse_mol_bytes_to_ir; TableIR parser families | Retain signatures. These are model-free parse/raise paths, so they retain small-ring stereo. There is currently no resolved MOL Interpret implementation; no new resolved MOL API is required by this inventory. |
| umol-graph/src/export.rs | Convey for Smiles / ReactionSmiles; export_smiles and export_reaction_smiles families; ConveyError | Retain public shapes. Read explicit frames, apply the ring-output policy, and reject unrepresentable ligand kinds, cosets, and entity-local constraints through existing errors. |

### DSL ownership and consumers

| Files or public surface | Required migration |
| --- | --- |
| umol-graph-ir/src/dsl/atom.rs, bond.rs | AtomDsl / BondDsl, update wrappers, predicate enums, and parse_atom / parse_bond / update entry points remain form-valued. Their standalone grammar accepts only presence states. Public AtomPredicate / BondPredicate need no configuration-payload variant. |
| umol-graph-ir/src/dsl/stereo.rs | Reduce the internal TetrahedralStereoDsl / CisTransStereoDsl wrappers and their parsers/formatters. Retain the full StereoAtomDsl / StereoBondDsl and configuration-expression grammar. The two site wrappers are not re-exported from dsl.rs. |
| umol-graph-ir/src/dsl/constraint.rs | AtomConstraintDsl / BondConstraintDsl and general constraint readers use reduced stereo values. Configuration expressions remain on explicit stereo entities. Molecule-level constraint semantics remain outside this work. |
| umol-graph-ir/src/dsl/molecule.rs | Retain the existing atom/bond input entry types, entity registration, and explicit stereo-record rendering. The shared atom/bond grammar rejects configuration-bearing #T/#C here as it does elsewhere. |
| umol-graph-ir/src/dsl/refs.rs and namespace.rs | Read and render :virt in tree and streaming paths. Preserve its bearer and kind in existing structural references and participant indexes, without new stereo id allocation rules. |
| umol-graph-ir/src/dsl/reaction.rs | LHS maps use the same reduced attribute grammar as molecule inputs. Delta values reject configuration-bearing #T/#C. Existing explicit stereo records and their reference handling remain. |
| umol-graph-ir/src/dsl/reaction_span.rs and edit.rs | Atom/bond sides and edit values reject configuration-bearing #T/#C. Explicit stereo records, participant references, and ligand-replacement payloads accept Virtual. |
| umol-graph-ir/src/dsl/config.rs; MoleculeDsl FromIr / IntoIr | Retain StereoDefault, AtomDefaults, BondDefaults, overrides, and aggregate defaults. There is no additional stereo-sugar configuration or graph-aware atom/bond attribute conversion. |
| umol-graph-ir/src/ir/molecule/spec.rs | Standalone atom/bond specs reject configuration syntax; explicit stereo spec terms remain available, including the new virtual argument. |
| umol-graph-ir/src/macros.rs | Existing DSL macros inherit molecule acceptance or standalone rejection through their owning parsers. |
| umol-graph-ir/spec/umol-dsl-spec.md | Specify the three presence values and :virt; remove configuration-bearing #T/#C from every context. Keep full configuration syntax on explicit stereo records and preserve tree/streaming parity. |

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
| umol-py/src/resolve.rs | Constructor keywords, getters, Rust conversion, and repr expose only stereo_atom_failure / stereo_bond_failure. |
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
| umol-graph-ir/src/ir/stereo.rs tests; tests/property/lattice.rs, strategies.rs | Replace reduced-form coset generators with all three presence states. Preserve equality/lattice/ground-literal laws. Keep full configuration generators and independent coset references on explicit entities. |
| ir/ligand.rs, ir/stereo/tests, ir/view/constraints/tests/stereo_atom.rs and stereo_bond.rs; tests/integrity.rs | Cover Virtual identity, bearer checks, duplicate ligands, exact matching, and known-kind readings with undetermined cosets. |
| dsl/atom.rs, bond.rs, stereo.rs, molecule/tests.rs, reaction.rs, reaction_span.rs, edit.rs; tests/property/stereo/serialization.rs and molecule/serialization.rs | Reject configuration-bearing #T/#C consistently in every context. Test explicit stereo records with literal, set, variable, and operator configurations, defaults, aliases, and exact tree/streaming/render roundtrips. Preserve unrelated kindless-payload equivalence contracts. |
| umol-graph-ir/src/macros.rs tests | Check uniform rejection of configuration-bearing #T/#C through the existing DSL macros. |
| tests/property/molecule/{canonicalize,remap,reframe,meet_pushout}.rs; property/substructure.rs; property/reaction and property/edit families | Retain relabeling, frame-transport, pushout, reaction, edit, and rollback laws. Add same-configuration versus enantiomer regression cases and Virtual cases without weakening exact gluing or old-state checks. |
| umol-graph/src/ops/stereo.rs, recognize/tests.rs, resolve/stereo.rs and resolve.rs; ops/validate/stereo.rs | Cover direct imported entities, virtual-kind replacement with unchanged coset, merged policies, unconditional mismatch contradictions, removal of general ring suppression, and lossy standalone projection. |
| ops/transform/stereo/tests, stereo/enumeration/tests.rs, hydrogens/tests; tests/property/transform/stereo and hydrogens.rs | Preserve selection, enumeration, reconstruction, and existing-entity laws. Move configuration descriptions to explicit entities. Isolated atom examples that only test presence use #T+. |
| umol-io/src/table_ir/raise.rs; CTfile and SMILES parser/render tests; umol-graph/src/ingest.rs and export.rs; tests/property/project.rs and resolve.rs | Update expected imported carriers. Cover alkene H, azo LP, mixed imine ligands, Either/unspecified stereo, output rejection, ring-policy placement, explicit frame transport, and molecule/reaction roundtrips. Aggregate projection assertions must now retain entities. |
| umol-graph/tests/resolution/data and tests/kekulization/data | Migrate the 59 fixture files using the old syntax: 21 in stereo_tetrahedral, 19 in stereo_cis_trans, 17 elsewhere in resolution data, and 2 in kekulization data. Preserve their intended chemistry and update carrier and ring-policy expectations. |
| umol-py/tests/test_stereo.py, test_stereo_perception.py, test_constraint.py, test_constraint_readings.py, test_structural_pattern_matching.py, test_class_constructor_signatures.py | Remove perception/derivation query and construction cases; preserve operation-required value coverage. Cover unit variants, bool input/output, removed imports and coercions, Virtual, merged keyword signatures, readings, and exact exception/roundtrip behavior. |
| umol-py/tests/test_molecule_stereo.py, test_molecule_stereo_enumeration.py, test_import.py, test_smiles.py, test_smiles_export.py, test_resolve.py, test_workflow.py; edit/delta/transaction tests | Exercise the existing high-level entry points and participant carriers after rebuilding the extension. Other Atom/Bond/Molecule view and constructor fixtures that construct old forms migrate with them. |
| umol-graph/benches/resolve.rs, transform.rs and transform/stereo.rs | Update instance projection calls, Solution assumptions, aggregate output checks, and configuration-bearing input strings. Use explicit stereo records while preserving the measured operation and iterator lifetime controls. |
| umol-graph-ir/benches/canonicalize.rs, constraint_readings.rs, editor.rs, integrity.rs; IO parsing/layout/SVG benchmarks | Add or adapt representative Virtual and explicit stereo fixtures where those operations are measured. Preserve the measured operation; parser or entity-only fixtures are not interchangeable workloads. |
| umol-graph-ir/fuzz/fuzz_targets/fuzz_entity_strings.rs, fuzz_molecule.rs, fuzz_reaction.rs, fuzz_reaction_span.rs and fuzz_edits.rs | Preserve parser and tree/streaming parity laws. Configuration-bearing #T/#C inputs become rejection cases in every context; add valid explicit stereo-record and Virtual seeds. Keep generated corpora as historical inputs, not bulk migration targets. |
| umol-graph-ir/spec/umol-dsl-spec.md; docs/chemistry/representation.md; docs/development/data-types.md and python-api.md; relevant nomenclature entries | Synchronize carrier, conversion, frame, and binding descriptions. Remove wording that equates a stereo entity or NotStereo with stereogenicity. |
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
  the stereoisomers. Keep the existing failing inline-canonicalization reproduction
  as historical evidence. S3 verifies that its obsolete input syntax is rejected
  and that explicit stereo records retain the intended behavior.

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

**Working proposal — 2026-11-08.** Recorded for continuation, not final design
approval. S2a–S2e are complete; implementation of the remaining subitems awaits
explicit instruction. Record completion and verification under each subitem.

Rename the Rust actor to StereoPerceiver. Its policy-free derive operation
checks supported stored frames and requests missing entities from positive
constraints. Recognize continues to return missing entities only; classify
continues to assess existing arrangements. Their Rust query contracts and
StereoPerceptionError remain separate from the changed derive interface.

The shared derivation/resolution interfaces are:

```rust
StereoPerceiver::derive(&Molecule)
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

StereoDerivation keeps its existing fields and row types:

```rust
stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>
stereo_bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>
inconsistencies: Vec<StereoInconsistency>
skipped_stereo_bonds: Vec<BondId> // retained until S2g
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
until S3; ring suppression remains until S2g.

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
  payload; do not infer the operation outcome from the presence of Virtual.

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
  umol-py/tests/test_import.py. [dep: S2e]

  Replace the wedge-derived tetrahedral and cis/trans constraint outputs with
  complete stereo entities. Transport source order and reference choices into
  the stored frame. Retain source-identified hydrogen/lone-pair ligands and use
  Virtual where the source does not identify the kind. Keep raise mechanical;
  invalid or impossible frames fail at the conversion boundary.

  Check exact imported entities for tetrahedral examples, alkenes, azo compounds,
  and imines, including source permutations and Either/unspecified stereo.
  Update import expectations in the owning tests, preserving the distinction
  between absent stereo and an entity with an unspecified configuration.

- **S2g — Move small-ring suppression to external-format I/O** — breaking (red→green).
  Modules: umol-graph/src/ops/stereo.rs,
  umol-graph/src/ops/stereo/recognize.rs,
  umol-graph/src/ingest.rs,
  umol-graph/src/export.rs,
  umol-graph/src/ops/model.rs. [dep: S2e, S2f]

  Remove ring-size suppression from general perception and resolution. Remove
  StereoPerceiver::skipped_stereo_bonds and the corresponding StereoDerivation
  field. The Python perception surface is removed in S2a. Retain the model
  setting and apply its policy after raise on resolved import. On export, filter
  the operation's private input before Resolver::project: that projection still
  removes stereo entities until S3d replaces the export path. This filter remains
  in place after S3d. The
  borrowed export input is preserved.

  Test resolved and unresolved molecule/reaction imports, later general
  resolution, explicit DSL stereo, and export at both sides of the threshold.
  Direct MOL raising remains unresolved and retains the imported entity. Update
  recognition, addition, and enumeration tests affected by removal of the general
  suppression rule.

**S2 gate:** external imports now carry their configurations on explicit entities,
virtual kinds resolve without changing configuration, and ring suppression has
only its agreed I/O scope. Rust consumers of the derivation outcomes and Rust/Python
operation consumers of the changed policies pass their tests. The site-value type
change is still S3.

### S3 — Remove configuration-bearing site constraints

This is the coordinated breaking stage. The value types, shared DSL grammar,
and their Rust/Python consumers must change together. The five subitems divide
that change by responsibility; the stage ends only when their combined build and
tests pass. Molecule entry storage, frame construction, entity id allocation,
and explicit-record rendering retain their existing design.

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
  umol-py/python/umol/__init__.py. [dep: S1a, S2g]

  Make Stereo a unit variant in the two site forms, use bool as their literal
  type, add From<bool>, and remove the ground enums, configuration conversions,
  and coset-specific comparisons. Update the Python variants, bool coercions,
  literal extraction, setters, registration, and exports in the same change.
  Derive presence from a known stereo kind without deriving a coset, and remove
  configuration data from site-constraint canonical keys.

  Test all three states and their lattice laws, bool extraction, open/complete
  readings, unknown kinds, undetermined cosets, per-kind negatives, and relabeling.
  Reduce the site-value generators while retaining the full entity-configuration
  generators and their independent references.

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

  Restrict #T/#C and their EDN values to the three presence states, including
  :stereo. Reject all configuration-bearing forms, including in molecule maps,
  reaction LHS maps, aliases, edits, deltas, and reaction-span sides. Retain the
  full configuration grammar on explicit stereo entities and the existing
  molecule parser, registration, and rendering paths.

  Test uniform rejection of literals, sets, variables, and operators in attribute
  positions. Check exact tree/streaming/render roundtrips of explicit stereo
  records, including symbolic configurations, defaults, aliases, entity keywords,
  and ids referenced by reaction deltas. The old canonicalization reproduction
  must reject at its attribute syntax; independently constructed explicit-entity
  reference cases continue to distinguish the configurations.

- **S3c — Remove site-coset handling from chemistry operations** — breaking (red→green).
  Modules: umol-graph/src/ops/stereo.rs,
  umol-graph/src/ops/stereo/recognize.rs,
  umol-graph/src/ops/resolve/stereo.rs,
  umol-graph/src/ops/validate/stereo.rs,
  umol-graph/src/ops/transform/stereo.rs,
  umol-graph/src/ops/transform/stereo/enumeration.rs,
  umol-graph/src/ops/transform/hydrogens.rs. [dep: S3a, S2d]

  Remove every remaining source of configuration data from site constraints.
  Resolve positive markers into missing entities with undetermined cosets.
  Preserve the addition/removal selectors and existing-entity behavior; obtain
  enumeration domains only from explicit configurations and retain entity-local
  constraint checks. Delete hydrogen-folding transport of site-constraint cosets
  while retaining transport and ligand conversion on explicit entities.

  Verify each operation's full outputs and failure behavior. Include retained
  imported configurations, unspecified and symbolic cosets, enumeration order
  and iterator ownership, selected/unselected sites, hydrogen reconstruction,
  and exact frame/configuration preservation. A Virtual ligand must not be
  treated as an implicit hydrogen during folding or unfolding.

- **S3d — Separate lossy stereo projection from export** — breaking (red→green).
  Modules: umol-graph/src/ops/resolve/stereo.rs,
  umol-graph/src/ops/resolve.rs,
  umol-graph/src/export.rs,
  umol-graph/tests/property/project.rs,
  umol-graph/benches/resolve.rs. [dep: S3a, S2g]

  Make standalone stereo projection model-free, using the agreed associated
  functions and direct Result returns. Replace supported entities with positive
  markers and remove obsolete frame errors. Remove the stereo stage and variants
  from aggregate projection so that it retains explicit stereo entities.

  Make Convey encode those entities directly, transporting cis/trans reference
  choices and reading their frames for hydrogen omission. Reject Virtual kinds,
  unsupported cosets, unpreservable entity constraints, and remaining site
  constraints that cannot be encoded. Preserve the I/O ring policy from S2g.

  Test standalone projection's deliberate loss and idempotence, negative-marker
  conflicts, aggregate projection's retention of entities, and exact molecule
  and reaction export roundtrips. Verify unchanged borrowed inputs on rejection
  and update the affected projection benchmark calls and correctness checks.

- **S3e — Complete the remaining consumer and fixture migration** — breaking (red→green).
  Modules: umol-graph-ir/tests/property/strategies.rs,
  umol-graph/tests/property/strategies.rs,
  umol-graph/tests/resolution/data,
  umol-graph/tests/kekulization/data,
  umol-py/src/molecule.rs,
  umol-py/src/reaction.rs,
  umol-py/tests, and the remaining test, fuzz, and benchmark files identified in
  the change inventory. [dep: S3b, S3c, S3d]

  Reconcile every inventoried consumer with its owning subitem. Preserve the
  intended chemistry while rewriting configuration descriptions as explicit
  stereo records and updating expected stored carriers. Replace isolated-atom
  configuration syntax with a presence marker only where the test
  concerns presence. Preserve the laws exercised by the property suites and the
  behavior exercised by parser seeds. Update remaining benchmark expectations
  to match the measured operation.

  Rebuild Python and test the existing high-level import, export, transformer,
  edit, delta, and transaction interfaces. Run the Rust suites, including the
  affected feature suites, and compile all benchmark targets. This subitem closes
  only when the coordinated change builds and its tests pass together.

**S3 gate:** site constraints contain only presence/kind information, all DSL
contexts reject configuration-bearing #T/#C, and configuration descriptions use
explicit records. Chemistry operations and export obtain configuration from those
entities. All affected callers compile and pass their tests; ligand identity and
matching remain exact.

### S4 — Documentation and final verification

- **S4a — Update the permanent descriptions and examples** — additive (green).
  Modules: umol-graph-ir/spec/umol-dsl-spec.md,
  docs/chemistry/representation.md,
  docs/development/data-types.md,
  docs/development/python-api.md,
  docs/development/nomenclature.md,
  .agents/skills/mol-notation/SKILL.md. [dep: S3e]

  Update the representation, construction, syntax, conversion, and binding
  contracts to match the implementation. Correct affected examples and wording
  that equates a stereo description with stereogenicity. Preserve shared-skill
  ownership and review the complete source/documentation diff against the
  inventory. Source comments and rustdoc describe current behavior without
  referring to discussion records.

- **S4b — Run final checks and compare measurements** — verification (green).
  Modules: the test and benchmark targets listed below, and
  discussion/243-stereo-site-assertions-2026-11-07.md,
  discussion/000-status.md. [dep: S4a]

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

The critical sequence is S0 → S1 → S2 → S3 → S4. Within S3, value reduction
feeds the shared DSL grammar, chemistry operations, and projection/export
consumers. Their complete caller and fixture migration restores the passing build.
There are no deferred stages in this deliverable. The exclusions in Scope remain
unchanged. S0 and S1 are complete; S2–S4 remain pending.
