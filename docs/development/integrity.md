# Representation integrity

## Purpose

This is the normative inventory and justification for eager representation-integrity checks. The
general construction and validation tiers are defined in [Data type contracts](data-types.md). This
guide records the narrower question: which properties must hold for the low-level aggregate types
to work coherently, and what concrete failure each check prevents.

The default preference is for open data types and first-requiring-operation validation. Low-level
types should impose only the constraints required to interpret and operate on their representation.
Do not build large defensive API moats, eagerly normalize inputs, or validate semantic properties
merely because a constructor has enough information to do so. Repeating defensive checks in every
method obscures the actual preconditions, adds cost, and makes the happy path harder to understand.

`Molecule`, `Reaction`, and `ReactionSpan` are the deliberate exception. Their entity families,
typed ids, participant frames, constraints, and two-sided projections interact across too many
operations for every operation to rediscover whether the stored representation is coherent. They
therefore establish a small tier-1 contract when a value is published. This closure is reluctant:
it is permission to enforce the minimum common representation contract, not permission to move
chemistry, satisfiability, normalization, or operation-specific validation into construction.

## Admission test for an integrity check

A property belongs to representation integrity only when omitting it would do at least one of the
following during ordinary operations on the type:

- expose an out-of-bounds lookup, degree assertion, impossible variant, or other internal panic;
- leave two stored fields without one coherent positional or referential interpretation;
- make an entity identity or lookup key non-unique where the API promises a single entity;
- allow frame transport, remapping, projection, or matching to return a silently incorrect result;
- force most operations on the type to repeat the same prerequisite check before they can proceed.

The check must reject the malformed representation at publication more clearly and cheaply than
rechecking it throughout the type. If a coherent stored value can violate the property and only
some operations care, the first operation that requires the property verifies it instead.

The failure descriptions below identify the primary concrete problem prevented by each current
check. They are not invitations to add checks for every hypothetical misuse. When a check no longer
prevents the named failure because the representation or operation changes, remove or reclassify it.

## Enforcement model

Independent assembly uses open carriers. Publication turns an accepted carrier into a closed
aggregate:

| Aggregate | Open input | Crate-private authoritative check | Checked publication | Asserted publication |
| --- | --- | --- | --- | --- |
| `Molecule` | `MoleculeEntries`; transient `MoleculeEditor` state | `Molecule::check_integrity` | `Molecule::try_from_entries`, `MoleculeEditor::snapshot`, `MoleculeEditor::try_build`, and checked integrity-sensitive mutations | `Molecule::from_entries`, `MoleculeEditor::build`, and trusted internal publishers |
| `Reaction` | a closed lhs `Molecule` plus independently assembled `Deltas` | `Reaction::check_integrity` | `Reaction::try_new` | `Reaction::new` and trusted internal publishers |
| `ReactionSpan` | `ReactionSpanEntries` | `ReactionSpan::check_integrity` | `ReactionSpan::try_from_entries` | `ReactionSpan::from_entries` and trusted internal publishers |

A checked route returns the aggregate's typed `*IntegrityError`. Its asserted sibling runs the same
check and panics when its producer contract is broken. Boundary adapters translate the checked
error into their own parse, conversion, or binding error; they do not reproduce the checks.

Once published, an operation accepting only a closed aggregate relies on the contract. Do not call
`check_integrity` defensively in every method. A trusted transformation must preserve the complete
contract by construction and test that preservation. A new public raw constructor or mutable escape
hatch is not harmless convenience: it reopens the container and would require defensive checks
throughout its operations.

Integrity-sensitive live mutation follows the same rule. Raw whole-form aromatic-system and
multicenter-bond mutation is restricted to graph IR. The public singular and family-wide
`try_modify_aromatic_system*` and `try_modify_multicenter_bond*` operations modify a private
candidate, publish it only after the authoritative check succeeds, return the exact
`MoleculeIntegrityError` on rejection, and leave the source unchanged. Python live views use those
checked operations and translate rejection to `ValueError`; they do not publish a partial change.

## `Molecule` integrity inventory

### References and parallel storage

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `InvalidReference` | A bond endpoint, relation participant or site, stereo-ligand anchor, or constraint refers to an entity outside the owning molecule. | Internal entity, relation, constraint, remapping, and projection code uses dense ids to index aggregate storage. A dangling stored id would otherwise become an out-of-bounds panic or be remapped as the wrong entity. |
| `ElectronCountLengthMismatch` | A literal aromatic or multicenter electron-count vector does not have one value per participant. | Counts are transported position by position with the participant frame. A mismatch currently turns reframing into `None` or leaves `permute` unchanged; accepting it would silently detach counts from atoms and later surface as an unrelated contradiction. |

### Fixed entity identity

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `DuplicateParticipant` | A participant frame or undistinguished factor repeats an actual atom: a bond or noncovalent self-loop, a repeated dative donor, a repeated aromatic or multicenter member, a repeated actual stereo ligand, or a stereo-atom site repeated as an actual ligand. A dative acceptor may also occur once as a donor because the two roles are distinguished factors. | Relation coincidence, incidence, and frame operations assume each individual participant frame is simple. Repetition within one frame would make identity and participant actions ambiguous or make one occurrence masquerade as two positions; occurrence across distinguished dative factors remains unambiguous. |
| `BondsParallel` | Two covalent bonds have the same unordered endpoint pair. | The graph-IR molecule gives a covalent bond identity by its endpoints. Single-edge lookup, correspondence induction, and bond matching would otherwise have multiple answers. |
| `DativeBondsIdentical` | Two dative bonds have the same acceptor and donor multiset, including when their stored donor orders differ. Shared acceptors or donors are permitted when the complete keys differ. | The complete `(acceptor, donor multiset)` is the dative identity and singular coincidence key. Duplicate complete keys would make lookup, correspondence, and delta targeting non-unique. |
| `NoncovalentBondsParallel` | Two noncovalent bonds have the same unordered endpoint pair, even if their kinds differ. | Noncovalent bond identity is the endpoint pair. Multiple entries would make `coincident_id`, matching, and delta targeting ambiguous. Combined interaction kinds must be represented in one form instead. |
| `AromaticSystemsOverlap` | One atom belongs to more than one aromatic system. | Aromatic membership names a unique owning system. Algorithms that recover the system from a member atom would otherwise choose an arbitrary incident relation. |
| `MulticenterBondsIdentical` | Two multicenter bonds have the same participant set. | The participant set is the multicenter bond's uniqueness key. Duplicate sets would make coincidence lookup, correspondence, and delta targeting non-unique. |
| `StereoAtomSitesDuplicate` | More than one stereo-atom entity is borne by the same atom. | The site identifies the stereo entity. Site-based constraints, lookup, perception, and reaction edits require one answer. |
| `StereoBondSitesDuplicate` | More than one stereo-bond entity is borne by the same bond. | The site bond identifies the stereo entity. Site-based constraints, lookup, perception, and reaction edits require one answer. |

### Stereo frames and domains

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `DuplicateStereoLigand` | A stereo frame repeats the same `StereoLigand`, including an identical implicit hydrogen or lone pair anchored at the same atom. | Equal frame positions do not determine a unique permutation action. Accepting them would require orbit search in every reframe, comparison, matching, pushout, and canonicalization path and could silently transport configurations or constraints by different actions. |
| `StereoFrameDegreeTooLarge` | A stereo frame has more than `umol_perm::MAX_DEGREE` ligands, whether or not a kind is asserted. | `Permutation` is a bounded representation whose constructors and actions assert the maximum degree. Rejecting at publication prevents degree assertions and fixed-array indexing failures in later frame operations. |
| `StereoLigandIncidenceMismatch` | A stereo-atom ligand is not borne by or adjacent to its site as required, or a stereo-bond frame is not two consecutive endpoint blocks with each ligand borne by or adjacent to the corresponding endpoint. | Stereo frames are site-relative, and bond frame actions preserve or swap whole endpoint blocks. Invalid incidence would attach a ligand value to the wrong site or endpoint and make reframing, matching, and application return incorrect stereochemistry. |
| `StereoKindSiteMismatch` | A kind is asserted on a site type that cannot carry its geometry. | The kind selects the action group and frame interpretation. Treating, for example, an atom as a `CisTrans` bond site would apply the wrong group despite the same numerical degree and produce incorrect transport. |
| `StereoLigandArity` | A kinded configuration or stereo constraint is paired with a frame whose length differs from the kind's degree. | Coset and permutation actions are defined for one degree. A mismatch would turn malformed input into `None` or `Contradiction`, or index a frame under the wrong action domain. |
| `StereoCosetOutOfRange` | A literal coset, literal set, or variable domain contains an index outside the kind's dense coset range. | The value does not name a configuration in the selected kind. Without the check, group action and normalization would report an unrelated failed action or contradiction instead of malformed representation. |
| `StereoPermutationDegree` | A ligand-symmetry, fluxionality, or expression permutation has a degree different from the frame or kind it acts on. | Applying the permutation would use positions from a different action domain, causing a failed action, an indexing panic in positional code, or incorrect constraint transport. |
| `StereoLigandPositionOutOfRange` | A topicity pair names a position outside its stereo frame. | Topicity evaluation and frame transport use these positions to address ligands or inverse actions. The check prevents out-of-bounds access and prevents a missing position from being misreported as no match or contradiction. |

These checks define only whether the stereo representation can be interpreted. They do not decide
whether the site is stereogenic, physically realizable, or accepted by a stereo model.

## `Reaction` integrity inventory

`Reaction` closes the representation formed by a closed lhs molecule and an open delta collection.
It does not require that the deltas can already materialize a consistent reaction span.

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `InvalidReference` | A delta or nested constraint refers to neither an lhs entity nor a uniquely added entity; an added id collides with the lhs or another addition. | Delta execution and remapping index entities by id, while removal integrity indexes source-frame maps after reference validation. Missing ids would panic; colliding created ids would give later deltas and correspondences two incompatible meanings. |
| `StereoIntegrityError` | A stereo addition has an invalid local frame, configuration, or inline constraint, or a stereo constraint wrapper asserts a kind inadmissible for its site type. | These values enter reaction operations without first becoming a molecule. Reusing the molecule-local checks prevents bounded-permutation panics, invalid position domains, and wrong site-group interpretation before span materialization. |
| `IncidenceMismatch` | An overlay removal names a site or structured participant incidence different from the lhs entity or same-reaction addition it removes. Factor-local reordering and complete stereo-bond endpoint-block exchange preserve incidence; moving individual ligands between blocks does not. | A removal id and its recorded incidence would describe different entities. Span conversion and application could then delete one entity while matching, transporting, or reporting another. |
| `StereoKindModified` | One stereo entity's configuration change replaces one determined kind with another. | The old and new configurations would require different action groups, so no single entity frame can transport the change coherently. A kind change is represented by removing the old entity and adding a new one. |

An overlay removal may record compatible incidence in a participant order different from its source.
That sequence is an explicit local frame, not malformed representation. Because complete participant
values cannot repeat, the entity kind determines one local-to-source action. Matching transports the
recorded payload through that action before comparing it with the source. Aggregate reaction
reframing composes the local alignment with the one owning action for the lhs or `Add` entity.

Reaction integrity intentionally does not establish delta normal form, old/new continuity,
constraint satisfiability, two-sided span materializability, DPO gluing conditions, host
applicability, or chemistry. The operation that first requires one of those properties checks it.

## `ReactionSpan` integrity inventory

`ReactionSpan` stores one union namespace and projects it to two closed molecules. The projection
checks are representation integrity because `lhs()` and `rhs()` promise infallible `Molecule`
results.

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `InvalidReference` | A union-frame participant, site, ligand, or constraint refers outside the union namespace. | Projection uses dense union-to-side maps and indexes them by stored ids. A missing union id would panic during map indexing before either molecule projection could report its own integrity error. |
| `Lhs` | The lhs projection fails any `Molecule` integrity check. | `ReactionSpan::lhs` uses the asserted `Molecule::from_entries` path. Establishing the projection at span publication prevents that infallible accessor and every lhs-consuming operation from panicking or receiving an incoherent molecule. |
| `Rhs` | The rhs projection fails any `Molecule` integrity check. | `ReactionSpan::rhs` uses the asserted `Molecule::from_entries` path. Establishing the projection at span publication prevents that infallible accessor and every rhs-consuming operation from panicking or receiving an incoherent molecule. |
| `StereoKindModified` | The two determined sides of one stereo entity assert different kinds against their shared participant frame. | Reframing the span requires one action for the complete entity. Different kinds have different admissible groups, so selecting or applying one action would misinterpret at least one side. A kind change uses removal plus addition. |

Reaction-span integrity does not establish a DPO dangling condition, reaction applicability,
chemistry, satisfiability, or canonical form.

## Properties deliberately kept lazy

The closed-container exception does not alter first-requiring-operation validation for other
properties. In particular, none of the following belongs to integrity merely because it can be
checked eagerly:

- normalization, deduplication, canonical ordering, repair, or loss of representable distinctions;
- groundness or satisfiability of forms and constraints;
- model-independent physical invariants or chemistry-model conformance;
- whether reaction deltas are mutually consistent or materialize a two-sided span;
- DPO gluing conditions, host-dependent applicability, matching, or product existence;
- canonical form, canonical equality, resolution, perception, or source-format interpretation.

An equivalent `Modified` reaction-span entry is not an integrity failure. Checked and asserted span
construction preserve that raw tag. The current internal reaction-span normalization path collapses
it to `Unchanged`, canonicalization invokes that path, and `ReactionSpan::superimpose` may emit the
standardized form directly under its separate deriving contract. Doc 214 will expose the same
normal-form behavior through the public normalization and reframing pipeline. A semantically invalid
but representation-coherent value remains representable until the named operation that needs the
stronger property is invoked.

## Maintenance rule

Any change to `MoleculeIntegrityError`, `ReactionIntegrityError`, or
`ReactionSpanIntegrityError` updates this guide in the same work item. Each added check requires:

1. a concrete malformed representation;
2. the specific panic, ambiguous interpretation, or incorrect result prevented;
3. one authoritative implementation at the owning aggregate boundary;
4. a focused regression for the checked boundary and, where separately meaningful, its asserted
   sibling;
5. confirmation that the rejected property is not chemistry, normalization, or an
   operation-specific precondition.

If the second item cannot be stated precisely, the proposed check does not belong in integrity.
When a public raw construction or mutation route is closed, remove defensive rechecks that existed
only to compensate for that route rather than retaining both the moat and its guards.
