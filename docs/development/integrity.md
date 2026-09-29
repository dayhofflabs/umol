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

`Molecule`, `Reaction`, and `ReactionSpan` are the deliberate exception. Their entity sets,
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
| `Molecule` | `MoleculeEntries`; transient editor or transaction state | `Molecule::check_integrity` | Molecule::try_from_entries, MoleculeEditor::probe/finish, Transaction::probe/commit/tracked_commit | Molecule::from_entries, MoleculeBuilder::build, and trusted internal publishers |
| `Reaction` | a closed lhs `Molecule` plus independently assembled `Deltas` | `Reaction::check_integrity` | `Reaction::try_new` | `Reaction::new` and trusted internal publishers |
| `ReactionSpan` | `ReactionSpanEntries` | `ReactionSpan::check_integrity` | `ReactionSpan::try_from_entries` | `ReactionSpan::from_entries` and trusted internal publishers |

A checked route reports the aggregate's typed `*IntegrityError`, directly or within its operation
error. Its asserted sibling runs the same check and panics when its producer contract is broken.
Boundary adapters translate the checked error into their own parse, conversion, or binding error;
they do not reproduce the checks.

MoleculeEditor::apply checks each edit's preconditions and returns a transient editor. Its probe
checks integrity before lending a molecule; failure retains the editor. Its finish checks integrity
before returning the owned molecule; failure drops the editor. Molecule::apply and tracked_apply
include finish. Transaction::apply also permits intermediate states that violate integrity;
commit checks integrity and restores transaction-entry state on failure. Molecule::transact and
tracked_transact include commit. A failed transaction probe leaves the transaction active.

Once published, an operation accepting only a closed aggregate relies on the contract. Do not call
`check_integrity` defensively in every method. A trusted transformation must preserve the complete
contract by construction and test that preservation. A new public raw constructor or mutable escape
hatch is not harmless convenience: it reopens the container and would require defensive checks
throughout its operations.

Entity attributes and entity-level constraints do not establish frame agreement. Constructors
preserve supplied electron counts, configurations, and constraints without padding, truncation,
normalization, or repair. Length, coset, kind/site, action-degree, and local constraint-position
checks belong to the operations that need them. Molecule-level constraints retain their integrity checks.

## `Molecule` integrity inventory

### References and parallel storage

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `InvalidReference` | A bond endpoint, relation participant or site, stereo-ligand anchor, or constraint refers to an entity outside the owning molecule. | Internal entity, relation, constraint, remapping, and projection code uses dense ids to index aggregate storage. A dangling stored id would otherwise become an out-of-bounds panic or be remapped as the wrong entity. |

### Fixed entity identity

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `DuplicateAtom` | An atom occurs twice in one entity's atom references: a bond or noncovalent self-loop, a dative atom repeated among the donors or between a donor and the acceptor, a repeated aromatic or multicenter member, or a stereo-atom site repeated as an actual atom ligand. | Relation coincidence, incidence, and frame operations assume that every actual-atom occurrence identifies a distinct participant. A dative donor equal to its acceptor creates parallel, differently labelled incidences between the same atom and dative entity; ordinary incidence matching can then select one edge for both roles and return an incorrect result. Other repetitions make identity and participant actions ambiguous or make one occurrence masquerade as two positions. |
| `ParallelBonds` | Two localized bonds have the same unordered endpoint pair. | The graph-IR molecule gives a localized bond identity by its endpoints. Single-edge lookup, correspondence induction, and bond matching would otherwise have multiple answers. |
| `IdenticalDativeBonds` | Two dative bonds have the same acceptor and donor multiset, including when their stored donor orders differ. Shared acceptors or donors are permitted when the complete keys differ. | The complete `(acceptor, donor multiset)` is the dative identity and singular coincidence key. Duplicate complete keys would make lookup, correspondence, and delta targeting non-unique. |
| `ParallelNoncovalentBonds` | Two noncovalent bonds have the same unordered endpoint pair, even if their kinds differ. | Noncovalent bond identity is the endpoint pair. Multiple entries would make `coincident_id`, matching, and delta targeting ambiguous. Combined interaction kinds must be represented in one form instead. |
| `AromaticSystemsOverlap` | One atom belongs to more than one aromatic system. | Aromatic membership names a unique owning system. Algorithms that recover the system from a member atom would otherwise choose an arbitrary incident relation. |
| `IdenticalMulticenterBonds` | Two multicenter bonds have the same participant set. | The participant set is the multicenter bond's uniqueness key. Duplicate sets would make coincidence lookup, correspondence, and delta targeting non-unique. |
| `DuplicateStereoAtomSites` | More than one stereo-atom entity is borne by the same atom. | The site identifies the stereo entity. Site-based constraints, lookup, perception, and reaction edits require one answer. |
| `DuplicateStereoBondSites` | More than one stereo-bond entity is borne by the same bond. | The site bond identifies the stereo entity. Site-based constraints, lookup, perception, and reaction edits require one answer. |

### Stereo frames and domains

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `DuplicateStereoLigand` | A stereo frame repeats the same `StereoLigand`, including an identical implicit hydrogen or lone pair anchored at the same atom. | Equal frame positions do not determine a unique permutation action. Accepting them would require orbit search in every reframe, comparison, matching, pushout, and canonicalization path and could silently transport configurations or constraints by different actions. |
| `StereoFrameDegreeTooLarge` | A stereo frame has more than `umol_perm::MAX_DEGREE` ligands, whether or not a kind is asserted. | `Permutation` is a bounded representation whose constructors and actions assert the maximum degree. Rejecting at publication prevents degree assertions and fixed-array indexing failures in later frame operations. |
| `StereoLigandIncidenceMismatch` | A stereo-atom ligand is not borne by or adjacent to its site as required, or a stereo-bond frame is not two consecutive endpoint blocks with each ligand borne by or adjacent to the corresponding endpoint. | Stereo frames are site-relative, and bond frame actions preserve or swap whole endpoint blocks. Invalid incidence would attach a ligand value to the wrong site or endpoint and make reframing, matching, and application return incorrect stereochemistry. |
| `StereoKindSiteMismatch` | A top-level stereo constraint asserts a kind inadmissible for its atom or bond site. | Top-level constraint interpretation requires the named site action group. |
| `StereoLigandArity` | A top-level stereo constraint declares a kind whose degree differs from its referenced frame. | Positional interpretation requires agreement with the referenced frame. |
| `StereoPermutationDegree` | A top-level stereo constraint has a permutation of the wrong degree. | Constraint transport composes that permutation with the referenced frame action. |
| `StereoLigandPositionOutOfRange` | A top-level topicity constraint names a position outside its referenced frame. | Constraint evaluation and transport index those positions. |

These checks establish structural frames and top-level constraint integrity. They do not decide
whether the site is stereogenic, physically realizable, or accepted by a stereo model.

## `Reaction` integrity inventory

`Reaction` closes the representation formed by a closed lhs molecule and an open delta collection.
It does not require that the deltas can already materialize a consistent reaction span.

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `InvalidReference` | A delta or nested constraint refers to neither an lhs entity nor a uniquely added entity. | Delta execution and remapping index entities by id, while removal integrity reads the source from the lhs or its addition after reference validation. A missing id would panic or select no source. |
| `DuplicateReference` | An `Add` uses an entity ID already present in the lhs or used by an earlier `Add`. | The same entity reference would name two different entities, giving later deltas and correspondences incompatible meanings. |
| `DuplicateAtom` | A stereo Add or Remove repeats an actual atom ligand, or a stereo-atom entry uses its site atom as an actual ligand. | Distinct actual atom occurrences are required for unambiguous incidence and frame actions. Virtual ligands anchored at the site remain allowed. |
| `DuplicateStereoLigand` | A stereo Add or Remove repeats the same complete ligand value. | Equal frame positions do not determine a unique permutation action for reaction transport. |
| `StereoFrameDegreeTooLarge` | A stereo Add or Remove has more ligands than the bounded permutation representation supports. | Frame-action construction would otherwise reach a degree assertion. |
| `StereoKindSiteMismatch` | A top-level stereo constraint asserts a kind inadmissible for its atom or bond site. | Top-level constraint interpretation requires the named site action group. |
| `StereoLigandArity` | A top-level stereo constraint declares a kind whose degree differs from its referenced frame. | Positional interpretation requires agreement with the referenced frame. |
| `StereoPermutationDegree` | A top-level stereo constraint has a permutation of the wrong degree. | Constraint transport composes that permutation with the referenced frame action. |
| `StereoLigandPositionOutOfRange` | A top-level topicity constraint names a position outside its referenced frame. | Constraint evaluation and transport index those positions. |
| `IncidenceMismatch` | A bond or overlay removal records endpoints, a site, or structured participant incidence different from the lhs entity or same-reaction addition it removes. Factor-local reordering and complete stereo-bond endpoint-block exchange preserve incidence; moving individual ligands between blocks does not. | A removal id and its recorded incidence would describe different entities. Span conversion and application could then delete one entity while matching, transporting, or reporting another. |

The local stereo failures are direct Reaction variants with the same fields
as their Molecule counterparts. Molecule and Reaction share the local validation
rules, but each aggregate owns its reference, incidence, and public error
contract. Reaction does not wrap `MoleculeIntegrityError` for a delta payload;
ReactionSpan's `Lhs` and `Rhs` variants still wrap actual failed Molecule
projections.

An overlay removal may record compatible incidence in a participant order different from its source.
That sequence is an explicit local frame, not malformed representation. Because complete participant
values cannot repeat, the entity kind determines one local-to-source action. Matching transports the
recorded payload through that action before comparing it with the source. Aggregate reaction
transport conjugates the one owning action for the lhs or `Add` entity by the local alignment, so
the removal retains the same relation to its owner. Normalization instead uses the local-to-source
action directly to align the removal with that owner before reframing.

Reaction integrity does not establish delta normal form, old/new continuity, constraint
satisfiability, two-sided span materializability, DPO gluing conditions, host applicability, or
chemistry. The operation that first requires each deferred property checks it. A ModifyField may
carry old and new stereo configurations of different kinds or configurations incompatible with
their frame. Application checks
whether that change can execute, and span conversion checks whether it can form one preserved
stereo entity.

## `ReactionSpan` integrity inventory

`ReactionSpan` stores one union namespace and projects it to two closed molecules. The projection
checks are representation integrity because `lhs()` and `rhs()` promise infallible `Molecule`
results.

| Error | Rejected representation | Concrete failure prevented |
| --- | --- | --- |
| `InvalidReference` | A union-frame participant, site, ligand, or constraint refers outside the union namespace. | Projection uses dense union-to-side maps and indexes them by stored ids. A missing union id would panic during map indexing before either molecule projection could report its own integrity error. |
| `Lhs` | The lhs projection fails any `Molecule` integrity check. | `ReactionSpan::lhs` uses the asserted `Molecule::from_entries` path. Establishing the projection at span publication prevents that infallible accessor and every lhs-consuming operation from panicking or receiving an incoherent molecule. |
| `Rhs` | The rhs projection fails any `Molecule` integrity check. | `ReactionSpan::rhs` uses the asserted `Molecule::from_entries` path. Establishing the projection at span publication prevents that infallible accessor and every rhs-consuming operation from panicking or receiving an incoherent molecule. |
| `StereoKindModified` | The two determined sides of one stereo entity assert different kinds against their shared participant frame. | One preserved stereo entity has one kind on both sides of a span. Equal frame degrees or compatible actions do not make distinct kinds the same entity; a kind change uses removal plus addition. |

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
construction preserve that raw tag. ReactionSpan::normalize collapses it to `Unchanged`;
canonicalization invokes normalization, and ReactionSpan::superimpose may emit that form directly.
A semantically invalid but representation-coherent value remains representable until the named
operation that needs the stronger property is invoked.

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
