# Data type contracts

## Purpose

This is a normative developer guide for deciding which properties belong in constructors,
converters, integrity checks, validators, and transformations. It applies primarily to aggregate
model types such as `Molecule`, `Reaction`, and `ReactionSpan`. Small value types may
enforce stronger invariants when those invariants define the value represented by the type.

The central rule is that an aggregate constructor establishes that its input can be represented; it
does not establish every useful property of the resulting value. Semantic properties are lazy and
must be requested through named operations.

An invariant is therefore enforced by the first operation that requires it, not by the earliest
operation capable of checking it. A caller may request an explicit validator earlier. If it does
not, a later conversion or operation reports the failure when that invariant becomes a precondition
of producing its result.

This distinction is already visible in the reaction representations: `Reaction` construction
does not require its deltas to materialize an actual two-sided reaction span. An operation that
requires such a span reports the failure at that boundary. Host-dependent DPO application
conditions are a separate contextual concern. This guide generalizes that decision and makes it
available independently of the reaction implementation history.

## Operation taxonomy

| Operation | Establishes | Does not implicitly do |
| --- | --- | --- |
| Construction | Representation and referential integrity | Semantic validation, canonicalization, resolution, repair |
| Conversion | Faithful representation in the target type | Semantic validation, normalization, silent loss or repair |
| Integrity check | A tier-1 representation contract | Semantic validation, mutation, or repair |
| Validation | A tier-2 invariant or tier-3 conformance property | Mutation or repair |
| Transformation | An explicitly named change of representation or state | Pretend that the input was already valid or canonical |

## Integrity checks and validation tiers

The three tiers are ordered by the information needed to state the failed property. The numbering
does not imply one API family: tier 1 is enforced by graph-IR-owned integrity checks, while tiers 2 and 3
are semantic validation in the chemistry layer.

| tier | category | required context | consequence of failure |
| --- | --- | --- | --- |
| 1 | integrity | the stored representation alone | the representation is malformed and cannot satisfy the type's construction contract |
| 2 | invariants | the representation and model-independent physical or mathematical rules | the representation is coherent but violates a universal invariant |
| 3 | conformance | the representation and a selected chemistry model | the representation is coherent but is not accepted by that model |

Boundary reference resolution is not a validation tier. Names, aliases, and other boundary-only
references that are erased during conversion must resolve before the target graph IR can be assembled;
failure belongs to that conversion.

Tier 1 is representation integrity. It is an invariant of the data type and is established by every
public constructor, conversion, and mutation that produces the type. It includes resolvable stored
references, required collection shapes, participant-frame arity, and values lying in the
kind-dependent domain needed to interpret their representation. A checked constructor returns an
error for a violation; an asserted constructor may panic when its documented producer contract is
broken.

For stereo, a ligand-frame length different from the declared kind's degree and a concrete coset
outside that kind's coset space are representation-integrity failures. The same applies to an
explicit coset set or variable domain containing an out-of-range member and to an explicit frame
permutation of the wrong degree. These values do not denote malformed-but-chemically-questionable
stereo configurations; they fail to denote configurations at all. Whether a structurally valid
configuration is stereogenic, physically realizable, or admitted by a selected stereo model belongs
to the later tiers.

Tier 2 contains universal conditions over an otherwise coherent representation, such as electron
and angular-momentum consistency and agreement between independently meaningful stored constraints
and their model-independent derived values. These states remain representable so that patterns,
rules, diagnostics, and explicit validation can operate on them.
Tier 3 contains choices such as valence, aromaticity, stereo, and connectivity models. A different
model may legitimately produce a different conformance result for the same structurally valid IR.

The composite order is integrity check, invariants validation, then conformance validation. Within
conformance, the cheap topology-only connectivity check runs before valence conformance. A later
tier may assume the earlier tiers only when its public contract says so or its producer establishes
them. It must not turn an earlier-tier failure into a chemistry contradiction.

### Invoking an integrity check

Representation integrity has one authoritative implementation in the crate that owns the graph IR. It
is exposed as `check_integrity` with a corresponding `*IntegrityError`; there is no `*Checker`
object. An integrity check returns `Result<(), *IntegrityError>`, never `Solution`,
`Underdetermined`, or `Contradictory`. It includes stored entity and constraint references, parallel
collection shapes such as participant and electron-count lengths, and kind-dependent data needed to
interpret a value such as stereo frame arity, coset domains, and permutation degree.

Every path that publishes an aggregate IR value uses that same implementation:

1. A boundary parser or format converter first resolves names, aliases, and other boundary-only
   references, assembles entries, then invokes the checked IR constructor. It translates the IR
   construction error into its boundary-specific parse or conversion error without reproducing the
   checks.
2. A public constructor from independently supplied entries invokes the checked implementation and
   returns its error. An asserted constructor invokes the same implementation and may panic only
   when its documented producer contract is broken.
3. A builder, editor, transaction, or multi-step transformation may hold incomplete intermediate
   state. The integrity gate is at the operation that publishes a `Molecule`, not at every
   primitive mutation. A trusted internal producer may take the asserted route only when its
   construction establishes the complete invariant.
4. A public operation that requires interpretable representation, including aggregate
   canonicalization, invokes the same check defensively if unchecked or compromised values can
   reach it. It returns a typed integrity error and never turns the failure into a chemistry
   contradiction or indexing panic.

For `Molecule`, this places the same gate behind direct checked entry construction, molecule-DSL
raise, TableIR raise, Python construction from entries, and finalization of builders and editors.
Reaction code uses it whenever it materializes a molecule side or projection; it does not grow a
reaction-specific copy of molecule integrity. A route that already establishes the contract may use
the asserted constructor, but the asserted and checked routes must share the implementation and the
set of accepted values.

Do not run a semantic validator at these boundaries. DSL and external-format raise may preserve an
atom `#T` or bond `#C` assertion before a stereo entity has been perceived or resolved. Validating
constraint satisfaction during construction would reject a representable input and silently impose
an operation order. Source-format checks such as wedge, directional-bond, or chirality-frame
interpretation likewise stay in the source-format layer; only the resulting graph IR's representation
is checked by the shared integrity gate.

Canonicalization is a transformation of a representation-integrity-valid value. It does not repair,
preserve as opaque data, or select a canonical form for malformed representation state. If an
unchecked or compromised value reaches a public canonicalization operation, the operation reports a
typed integrity error and does not panic. Canonicalization otherwise preserves tier-2 and tier-3
invalid states; an intrinsic lattice `Contradiction` remains distinct from malformed storage.

## Representation ownership and crate layering

An operation belongs with the crate that owns the representation when its correctness requires
coordinated access to that representation's complete internal shape. For `Molecule` and the
reaction representations, this includes operations that must keep all entity families, participant frames,
typed ids, and references inside constraints synchronized while rebuilding or remapping the value.
Such operations belong in `umol-graph-ir`; moving them to `umol-graph` because a higher-level model
currently contains one of their inputs would invite the higher layer to reconstruct graph-IR internals.

The layers divide responsibility as follows:

- `umol-graph-core` supplies chemistry-independent graph and relation primitives and takes explicit
  algorithm selectors;
- `umol-graph-ir` owns graph-IR construction, complete incidence encoding, frame and id remapping, and
  transformations that must preserve the whole representation; and
- `umol-graph` supplies chemistry models and model-aware operations over the public graph-IR API. It may
  construct a public graph-IR operation context from a model and config, but it must not reach into IR
  storage or maintain a parallel representation of the same operation.

When a graph-IR-owned operation needs only one semantic consequence of a higher-level model, expose that
requirement directly in the graph-IR operation context. Do not move the operation upward or pass the
entire higher-level model merely because that model is the current source of the value. Conversely,
do not duplicate a reduced model type when the operation needs only a primitive parameter.

Aggregate canonicalization is the reference case. `umol-graph-core` supplies canonical labels and
the automorphism algorithm selector. `umol-graph-ir` owns typed incidence, stereo-frame actions,
complete remapping, constraint transport, and canonical representative construction.
`umol-graph` may project `StereoModel::para_stereo` and its operation config into the public graph-IR
context, but canonicalization itself remains in `umol-graph-ir`.

### Construction

Construction enforces the minimum invariants needed for the value to be a coherent instance of its
type. For aggregate graph-IR values, these are representation invariants:

- stored references name entities in the owning namespace;
- participant, site, ligand, anchor, and constraint references are resolvable;
- correspondence pairs are in range and form a partial bijection;
- parallel collections have the shape required by the representation;
- a representation variant contains the data required to interpret that variant.

Construction does not establish model-independent semantics merely because they can be checked
without a chemistry model. In particular, construction does not imply:

- that a reaction's deltas can materialize an actual two-sided span;
- host-dependent application conditions such as the DPO gluing conditions;
- chemical validity or conformance to a model;
- agreement between independently meaningful entities and constraints;
- groundness;
- canonical form;
- satisfaction of an operation's preconditions.

An asserted constructor such as `from_entries` and a checked constructor such as
`try_from_entries` establish the same invariant. They differ only in how a violated construction
contract is reported. The asserted form is for producers that establish the invariant by
construction; the checked form is for untrusted or independently assembled input.

### Conversion

A conversion preserves every source state representable by its target, including semantically
invalid states. It fails only when the target representation cannot structurally encode the source
or when the source contains contradictory information that cannot be represented as one target
value.

A conversion must not silently make its input acceptable. Dropping dangling entities, removing
constraints, selecting a semantic interpretation, canonicalizing, or resolving undetermined state
are transformations rather than conversion mechanics. If such behavior is useful, expose it as a
separately named operation.

The error belongs to the failed boundary:

- failure to assemble the target's entries uses the target entry error;
- contradictory source information uses the relevant contradiction or a conversion error that
  preserves that cause;
- failure of an optional semantic precondition belongs to the operation requesting that property,
  not to the target constructor.

### Raw representation and semantic normal form

Derived structural equality may observe storage details that are deliberately absent from the
semantic model. Their presence in a raw IR value does not require every operation or intermediate
representation to preserve them. In particular, do not introduce ordering or multiset semantics
solely to make an exact structural roundtrip retain redundant input.

Molecule-level `Constraints` is the reference case. Its raw `Vec` representation permits ordering
and duplicate entries, and `==` observes both. Its semantic value is an implicit conjunction:
conjunct order is irrelevant and repeated equal constraints are idempotent. Normalization
therefore normalizes the individual constraints, flattens conjunctions, sorts them, and removes
duplicates. Reaction deltas, reaction spans, matching, and aggregate canonicalization must use this
set-like semantic value; they must not treat duplicate occurrences as distinct chemical assertions.
Context-free constraint-delta normalization deduplicates repeated additions or removals and cancels
an equal addition/removal pair. Materializing the normalized delta against an LHS then enforces
continuity: adding an already-present canonical constraint or removing an absent one is a
`Contradiction`, not an idempotent no-op. This matches old/new continuity for entity changes and
keeps addition and removal as meaningful inverses.

This does not license silent normalization in an operation that promises faithful representation
conversion. Instead, classify the operation accurately:

- a boundary parse/render roundtrip or a conversion documented as representation-preserving retains
  every representable raw distinction;
- an operation that derives a semantic representation or delta normal form may discard redundant
  ordering and multiplicity, but must document the normal form and state its roundtrip under the
  relevant semantic equality rather than raw `==`; and
- a transaction whose contract is exact rollback may retain positions and occurrences in its undo
  data without making those details part of the semantic constraint model.

When exact raw preservation conflicts with a smaller, closed, or algebraically coherent semantic
API, preserve the semantic API and weaken the raw roundtrip to the documented canonical relation.
Do not add variants, counters, or special cases merely to preserve redundant syntax or collection
layout.

### Validation

Validation checks one named tier-2 invariant or tier-3 conformance property without mutation. The
validator concept belongs in `umol-graph`, where `*InvariantsValidator` names model-independent
validation and `*ConformanceValidator` names validation under a selected chemistry model. Both
return `Result<Solution<_, _>, _>` because a coherent non-ground representation may leave a semantic
question underdetermined.

Validation remains separate even when a property is model-independent. Whether reaction deltas can
materialize a two-sided span is a precondition of operations that require that representation, not
a DPO condition. DPO dangling and identification conditions depend on a host and match and belong
to reaction application. Canonicality is similarly useful but not required for representation.

Do not define `*Validator` types in `umol-graph-ir`. A former IR validator is either a tier-1 integrity
check that stays with its owning IR type and returns `Result`, or a tier-2/tier-3 validator that
belongs in `umol-graph`. The composite graph `Validator` may group invariants and conformance passes,
but it does not expose `validate_integrity`: successfully constructed IR inputs already satisfy
their integrity contract.

### Transformation

Canonicalization, resolution, repair, stripping, cascading removal, and closure under a rewriting
semantics are transformations. Their names and return types must expose the change. A constructor or
ordinary conversion must not perform one as an incidental implementation step.

## Identity and constraints

Derived `==` compares the stored IR structure exactly, including constraints, ids, ordering, and
non-normal value encodings. `equiv` compares normalized forms in the current frame;
`equiv_under` performs the same comparison after an explicit entity or participant-frame mapping.
Aggregate `canonical_eq` compares complete canonical IR values after selecting the canonical frame.
All three semantic comparisons include constraints. This distinction matters for patterns, where
constraints are not redundant with the structural description.

Entity and molecular structural identity is established from inherent fields and structural
incidence. Constraints restrict the states admitted by an entity or molecule but do not establish
that structural identity. They therefore do not distinguish the initial structural automorphism
orbits. Aggregate canonicalization may use normalized constraints after structural labeling to
select among the remaining structurally equivalent frames. Remapping transports every reference
they contain. This post-hoc participation gives the complete IR assertion a unique canonical form
without turning constraints into structural identity features.

## Aggregate canonicalization

> **TODO (2026-08-07):** This section records the approved target contract from discussion doc 186.
> Fixed-frame normalization now uses `Normalize`; the context-bearing aggregate trait remains to be
> implemented. Remove this marker when doc 186 is implemented.

Aggregate canonicalization selects an entity-id and participant frame for a complete indexed graph IR.
It is distinct from `Normalize`, which puts values into normal form without changing their id or
participant frame. `Normalize` is context-free and remains the supertrait of `Lattice`;
`Canonicalize` takes an explicit canonicalization context and is implemented by complete indexed
aggregates rather than their leaf values. Do not add a context parameter to `Normalize` merely to
share a trait between the two operations.

The aggregate operation preserves the represented molecule: it remaps every entity and reference,
transports position-sensitive relation data, carries stereo cosets through stereo-frame
permutations, and normalizes every carried value in the selected frame. It does not perceive,
resolve, strip, repair, or validate chemistry.

The equality operations form three levels:

- `==` compares the exact stored representation;
- `equiv` compares normalized values in the current frame, while `equiv_under` applies an explicitly
  supplied correspondence or participant order before that comparison; and
- `canonical_eq` compares complete aggregate canonical forms under a shared context, selecting the
  frame rather than receiving it from the caller.

For inputs in the aggregate operation domain, `canonical_eq` holds exactly when canonicalization
produces the same complete IR value. Equivalently, an admissible remapping exists under which
`equiv_under` holds. `equiv_under` is therefore the explicit-map member of the `equiv` family, not a
fourth equality relation.

The context-bearing trait has the semantic shape

```rust
pub trait Canonicalize: Sized {
    type Error;

    fn canonicalize(
        self,
        context: &CanonicalizationContext,
    ) -> Result<Self, Self::Error>;

    fn canonicalize_by(
        self,
        level: CanonicalizationLevel,
        context: &CanonicalizationContext,
    ) -> Result<Self, Self::Error>;

    fn canonical_eq(
        &self,
        other: &Self,
        context: &CanonicalizationContext,
    ) -> bool;

    fn canonical_eq_by(
        &self,
        other: &Self,
        level: CanonicalizationLevel,
        context: &CanonicalizationContext,
    ) -> bool;
}
```

Use one concrete context for `Molecule`, `ReactionSpan`, and `Reaction` unless an
implementation establishes a real need for distinct context types. Canonical-form construction is
fallible; equality is total. Structurally equal inputs compare equal immediately. Two successfully
canonicalized inputs compare by structural equality, and two intrinsic contradictions compare
equal because both denote the empty semantic value. One contradiction and one successful form do
not compare equal. Integrity failures never make distinct inputs equal; callers that need the
diagnostic invoke `canonicalize` directly.

Stored stereo entities participate whether or not para-stereo refinement is enabled. Without
para-stereo refinement, perform one stereo-sensitive refinement from the constitution-level
partition. With it enabled, feed each stereo-sensitive partition back into refinement until the
partition stabilizes. Do not expose an iteration cutoff that can change the canonical form.

A parameterized operation may select a frame using only a coarser structural layer while returning
the complete original molecule in that frame. Its guarantee is deliberately limited:

- the selected structural layer is in canonical form;
- the complete result is a remapping of the input and retains excluded features semantically,
  while still normalizing their carried form values; and
- the ordering of excluded features within an automorphism class of the selected layer is not
  determined.

An excluded feature must not break such a tie. Complete outputs from differently numbered inputs
may therefore differ by an automorphism of the selected layer. Complete `canonical_eq` does not use
this coarser operation; it compares canonical representatives formed from every available entity
kind and constraint.

Because `canonicalize_by` returns the complete normalized aggregate, an intrinsic contradiction in
excluded data still makes that transformation fail. `canonical_eq_by` compares only the selected
layer and must not be implemented by comparing complete `canonicalize_by` results; contradictions
outside the selected layer do not affect that reduced relation.

Backend canonical labels are search inputs, not the stable numbering contract. The canonical frame
is the minimum under the library's typed comparison order. Automorphism generators and orbits may
prune equivalent branches, but changing the selected graph algorithm must not change the resulting
canonical representative.

The typed comparison schema is a compatibility contract, but that contract must admit additive
entity-model extensions. Existing entity-kind blocks, field components, constraint variants, and
their order have explicit stable schema positions. New entity kinds and constraint variants occupy
append-only extension positions and do not renumber or reinterpret existing positions. Schema
positions must not be inferred from Rust enum declaration order.

Adding an extension must leave the canonical numbering and canonical form of every molecule that
contains none of the new entity kinds or constraint variants unchanged. Molecules that use the new
extension acquire a coherent order that is frozen from the version introducing the extension
onward. No comparison exists with versions that could not represent that extension. More generally,
if one schema is an append-only extension of another, every molecule expressible in the earlier
schema has the same canonical representative under both. This cumulative promise concerns the
canonical IR, not an internal comparison-key encoding.

### Canonical comparison schema

The following schema is the normative typed order for aggregate canonicalization. Numeric positions
are local to the table in which they occur. A position is permanent once published: declarations may
move in Rust, but an existing position must not be reused or assigned another meaning. New entries
are appended after the highest assigned position in their table or structural domain.

The entity model has three ordered structural domains. Topology is AB, non-stereo is DAMN, and
stereo is SS. Constitution is topology plus non-stereo. Overlays are non-stereo plus stereo. Thus the
public cumulative canonicalization levels are `Topology`, `Constitution`, and `Full`, where `Full`
is topology plus overlays. `NonStereo` names the middle entity domain; it is not another cumulative
level.

| Domain position | Structural domain | Entity slots |
| ---: | --- | --- |
| 0 | Topology | Atom = 0, Bond = 1 |
| 1 | NonStereo | Dative bond = 0, Aromatic system = 1, Multicenter bond = 2, Noncovalent bond = 3 |
| 2 | Stereo | Stereo atom = 0, Stereo bond = 1 |

An entity-block position is the composite `(domain position, entity slot)`. Entity blocks compare by
this position and then by their dense row sequence. This hierarchy makes `Topology`, `Constitution`,
and `Full` exact domain prefixes while allowing a future entity kind to be appended within the
domain to which it belongs. An absent block contributes no entry. Constraints form a separate
terminal section, so extending an entity domain does not move the constraint section or alter keys
for molecules that lack the new kind.

Rows compare their components in the following local field order. The dense row index is implicit
in the row sequence. Inline constraints are excluded here and enter through the constraint section.

| Entity row | Components from position 0 onward |
| --- | --- |
| Atom | element, isotope mass, charge, implicit hydrogens, lone pairs, unpaired electrons |
| Bond | endpoint pair, order, charge, unpaired electrons |
| Dative bond | donors, acceptor, order |
| Aromatic system | participants, participant electron counts, charge, unpaired electrons |
| Multicenter bond | participants, participant electron counts, charge, unpaired electrons |
| Noncovalent bond | endpoint pair, kind |
| Stereo atom | site, ligand frame, configuration |
| Stereo bond | site, ligand frame, configuration |

Canonical-search initial classes retain these published field positions even when participant data
is represented by incidence occurrences. Omitted participant-bearing fields are not renumbered: for
example, a bond node uses positions 1 through 3 for order, charge, and unpaired electrons, while its
endpoint pair is represented by two incidences. Entity-node classes contain only normalized,
constraint-free, frame-independent values. Aromatic and multicenter participant electron counts and
stereo ligand kinds occur on their corresponding incidences; raw stereo configurations do not enter
the initial node classes.

Typed incidences use the following frozen order:

| Position | Incidence |
| ---: | --- |
| 0 | Bond endpoint |
| 1 | Dative donor |
| 2 | Dative acceptor |
| 3 | Aromatic participant, followed by its normalized electron-count value |
| 4 | Multicenter participant, followed by its normalized electron-count value |
| 5 | Noncovalent endpoint |
| 6 | Stereo site |
| 7 | Stereo ligand, followed by its ligand kind |

The public `Incidence` total order follows this table and agrees with the typed canonical key for
normalized incidence values. Entity-node and incidence classes occupy disjoint key domains.

Endpoint pairs and unordered participant sets are in their normalized participant order. Dative
donors are ordered independently of the acceptor. A stereo ligand is the product `(atom id, ligand
kind)`, with ligand kinds `Atom = 0`, `ImplicitHydrogen = 1`, and `LonePair = 2`. Unpaired electrons
are the product `(count = 0, multiplicity = 1)`. A stereo configuration is the product of its kind
and coset after the complete participant-frame action described above.

Products compare components in position order. Variants compare their explicit tag first and then
their payload components in written order. Sequences compare lexicographically, with a proper prefix
sorting before the longer sequence. Sets compare as sorted, duplicate-free member sequences. Options
use `None = 0` and `Some = 1`, followed by the value. Booleans use `false < true`; signed and unsigned
integers use numeric order; strings use lexicographic Unicode scalar-value order; entity ids use
their numeric indices. Elements use atomic-number order. A permutation compares as
`(degree, one-line image)`, and orientation uses `Proper = 0`, `Improper = 1`.

The normalized form variants have these tag positions:

| Type | Variant positions from 0 onward |
| --- | --- |
| `BooleanForm` | `Undetermined`, `Lit` |
| `NumForm` | `Undetermined`, `Lit`, `LitSet`, `RangeFrom`, `RangeTo`, `ArithExpr`, `PredExpr` |
| `ArithExpr` | `Lit`, `Var`, `Neg`, `Sum`, `Product`, `Div`, `Rem` |
| `PredExpr` | `Rel`, `Mem`, `Not`, `And`, `Or` |
| `RelOp` | `Le`, `Ge`, `Eq`, `Lt`, `Gt`, `Ne` |
| `MemOp` | `In`, `NotIn` |
| `ElementForm` | `Undetermined`, `Lit`, `LitSet`, `NotSet`, `Var` |
| `IsotopeMassForm` | `Undetermined`, `Natural`, `Lit`, `LitSet`, `Var` |
| `ElectronCountsForm` | `Undetermined`, `Lit` |
| `AromaticValenceForm` | `Undetermined`, `NotAromatic`, `Aromatic` |
| `MulticenterValenceForm` | `Undetermined`, `NotMulticenter`, `Multicenter` |
| `TetrahedralStereoForm`, `CisTransStereoForm` | `Undetermined`, `NotStereo`, `Stereo` |
| `NoncovalentBondKindForm` | `Undetermined`, `Lit` |
| `NoncovalentBondKind` | `HydrogenBond`, `HalogenBond`, `ChalcogenBond`, `Ionic`, `VanDerWaals` |
| `StereoConfigurationForm` | `Undetermined`, `Kinded` |
| `StereoKind` | `Tetrahedral`, `CisTrans`, `Axial`, `SquarePlanar`, `TrigonalBipyramidal`, `Octahedral` |
| `StereoCoset` | `Undetermined`, `Lit`, `LitSet`, `Term` |
| `StereoTerm` | `Var`, `Lit`, `LitSet`, `Swap`, `Mirror`, `Apply` |
| `RingScope` | `All`, `Size` |
| `TopicityRelationForm`, `StereogenicityForm` | `Undetermined`, `Lit`, `LitSet`, `NotSet` |
| `Topicity` | `Homotopic`, `Enantiotopic`, `Diastereotopic` |
| `Stereogenicity` | `Symmetric`, `Prochiral`, `Stereogenic` |

Inline constraint blocks use the same composite entity-block positions and therefore the same
topology, non-stereo, and stereo domain hierarchy. The molecule-level constraint tree is a terminal
block after every inline block. An inline constraint row is `(entity id, constraint tag, payload)`.
The currently assigned inline constraint tags are:

| Constraint block | Variant positions from 0 onward |
| --- | --- |
| Atom | `Valence`, `DonatedPairs`, `AcceptedPairs`, `AromaticValence`, `MulticenterValence`, `TetrahedralStereo`, `Degree`, `TotalDegree`, `TotalValence`, `RingDegree`, `RingValence`, `TotalHydrogens`, `RingMembership` |
| Bond | `Aromatic`, `CisTransStereo`, `RingMembership` |
| Dative bond | `Aromatic`, `RingMembership` |
| Aromatic system | `ElectronCount` |
| Multicenter bond | `ElectronCount` |
| Noncovalent bond | `Intramolecular` |
| Stereo atom, stereo bond | `LigandSymmetry`, `Fluxionality`, `Topicity`, `Stereogenicity` |

The molecule-level `Constraint` comparison key uses `EntityLeaf = 0`, `Relational = 1`,
`Molecule = 2`, `And = 3`, `Or = 4`, and `Not = 5`. Every public entity-leaf variant maps to
`EntityLeaf`, followed by its composite entity-block position, referenced entity id, stereo kind
where present, and inline constraint key. This comparison-key grouping does not change the public
`Constraint` enum. `And` and `Or` compare their normalized set-like child sequences.

`MoleculeConstraint` uses `ChargeSum = 0`, `UnpairedElectronCoupling = 1`, `BondOrderSum = 2`,
and `Connected = 3`. Payload fields compare in their declaration's semantic order: the optional
entity subset first, then the asserted value where one exists.

`RelationalConstraint` uses the composite position `(owning entity-block position, local slot)`.
The assigned local slots are:

| Owning entity | Local slots and variants in order |
| --- | --- |
| Dative bond | 0–7: `DativeBondDonors`, `DativeBondDonor`, `DativeBondContainsAllDonors`, `DativeBondAllDonors`, `DativeBondAnyDonor`, `DativeBondAcceptor`, `DativeBondAcceptorSatisfies`, `DativeBondParallels` |
| Aromatic system | 0–4: `AromaticSystemAtoms`, `AromaticSystemContains`, `AromaticSystemContainsAll`, `AromaticSystemAllAtoms`, `AromaticSystemAnyAtom` |
| Multicenter bond | 0–4: `MulticenterBondAtoms`, `MulticenterBondContains`, `MulticenterBondContainsAll`, `MulticenterBondAllAtoms`, `MulticenterBondAnyAtom` |
| Noncovalent bond | 0–2: `NoncovalentBondEnds`, `NoncovalentBondContains`, `NoncovalentBondEndsSatisfy` |
| Stereo atom | 0–4: `StereoAtomSite`, `StereoAtomContains`, `StereoAtomLigands`, `StereoAtomAllLigands`, `StereoAtomAnyLigand` |
| Stereo bond | 0–4: `StereoBondSite`, `StereoBondContains`, `StereoBondLigands`, `StereoBondAllLigands`, `StereoBondAnyLigand` |

Each relational payload compares its named fields from left to right as shown by the public variant.
Atom collections with set semantics are normalized before comparison; ordered predicate pairs retain
their order.

Reaction spans reuse the complete molecule schema and lift each entity row through the span tags
`Unchanged = 0`, `Added = 1`, `Removed = 2`, and `Modified = 3`. `Modified` compares its lhs value
before its rhs value. Constraint spans use the first three tags. These positions are deliberately
independent of the Rust declaration order.

An absent extension contributes no positioned entry. Therefore appending a later field, variant,
constraint, or entity kind within its assigned domain leaves the key of every earlier-schema value
byte-for-byte equivalent at the typed-key level. A genuinely new structural category may be
appended as a domain after stereo. Moving an existing entity kind between domains or inserting a
domain is schema-breaking because it changes both comparison order and cumulative-level semantics.
The concrete Rust key storage remains private and may change, but exact ordering tests must
instantiate these published positions and verify this append-only property.

## Provenance and contextual validity

Whether a consumer must re-establish a contextual property depends on what kind of value it
receives.

An **open data carrier** can be constructed independently of the objects with which it will later be
used. Its intrinsic constructor checks only its own representation, but a contextual consumer checks
the properties required to combine it with those objects. `Correspondence` is such a carrier: atom
pairs may come from SMIRKS or another mapped external format, and Rust and Python both permit direct
construction. A valid partial bijection is not necessarily a correspondence over the particular two
molecules supplied to a later operation. `MoleculeCorrespondence` is likewise not bound to molecule
instances by its ids.

An **operation-issued value** may instead be provenance-bound. A `Transaction` is issued by applying
edits and records how to undo that particular successful application. It has no public constructor
for independently asserting that provenance. Replacing it with another transaction or mutating the
object independently violates the operation contract; rollback must not panic, but it does not owe
correct restoration for the compromised pairing.

The resulting rules are:

- do not add validation merely to defend against swapping an opaque, provenance-bound result into a
  different operation history;
- do validate an open carrier when a public operation legitimately accepts it alongside independent
  objects and needs contextual agreement to produce a correct result;
- an internal producer-consumer path may assert the contextual property established by its producer;
- simple ids and open correspondences remain deliberately unbound to object identity; this avoids
  heavier identity infrastructure but makes contextual checks the responsibility of operations that
  promise a relationship between independently supplied objects.

For molecule correspondences, an atom correspondence read from an external format is a normal input
to induction over a supplied molecule pair. Count agreement and any structural uniqueness required
to derive the remaining entity families are therefore operation preconditions, not defenses against
deliberate tampering. A full correspondence produced for the same molecule pair by a conforming
operation satisfies those checks by provenance. Reusing that result with the same unchanged pair
cannot newly produce a contextual mismatch. Supplying an atom correspondence and molecule pair
independently is different: this is the ordinary bridge from mapped formats such as SMIRKS, and
induction may report incompatible carrier sizes or non-unique entity incidence. Supplying a full
correspondence for another molecule pair is likewise an ordinary public-input case whenever the
consuming API accepts the correspondence and pair independently; it is not treated as tampering
with an opaque operation result.

### Containing fallibility

Fallibility is not propagated merely because an implementation calls another fallible operation.
It belongs on the first public operation whose promised result cannot be produced without the
property:

- use `Option` when there is one ordinary absence condition and the cause carries no useful detail;
- use `Result` when the caller can distinguish or act on failure causes;
- preserve an existing operation-specific error surface when an internal implementation route
  changes;
- assert only on internal paths whose producer establishes the required property.

For example, converting a partial correspondence to a total-on-source remapping may return `None`
because no total-left mapping exists. This can occur for a correspondence correctly produced for
its molecule pair: partiality is part of the correspondence model, whereas a remapping is total on
its source.
It does not make correspondence construction, composition, reversal, or unrelated consumers
fallible. Exact error taxonomy remains subject to the repository-wide error review; the
construction/validation boundary does not require introducing a new error type for each method.

## Remapping

Remapping is an explicit transformation between id spaces. It transports represented values and
incidence along a total function; it does not validate chemistry, normalize attributes, repair
references, or remove entities. The image vectors passed to a remapping constructor define its
source domains, and every id in those domains has an image. Construction is therefore infallible.
A general remapping need not be injective or surjective; injection or bijection is a contextual
precondition of operations that require distinct or dense target entities. Removal uses compaction
instead because a removed id has no image.

The facility has two coordinated levels:

- `umol_graph_core::Remapping` maps node and edge ids used by graphs and relation participants. A
  relation-set remapping must relabel each factor, canonicalize that factor according to its
  `Ordered` or `Unordered` marker, and apply the induced position permutation to its `RelationData`.
  Positional payloads therefore remain aligned with their participants. An ordered stereo ligand
  frame retains its positions under a pure id remapping, so its induced position permutation is the
  identity and its coset is unchanged.
- `IdRemapping` maps all eight molecule entity-id families. It is used for graph-IR values that contain
  entity references, including constraints and deltas; it does not duplicate graph-core participant
  canonicalization.

A higher-level operation that moves molecular data into another namespace derives both mappings
from the same correspondence or construction result: the graph-core mapping transports topology
and relation participants, while the IR mapping transports references to owned entity rows. Do not
manually sort remapped relation participants or permute their payloads at individual call sites;
that behavior belongs to relation-set remapping through `RelationData`.

### Participant frames and payload equivalence

`RelationData::on_permutation` and `BiRelationData::on_permutation` are the sole primitive actions
of participant-position permutations on relation payloads. `Equiv` and `BiEquiv` derive comparison
from that action rather than defining a second remapping protocol:

- `equiv` compares two payloads expressed in the same participant frame;
- `equiv_under` first expresses `self` in the other payload's frame and then performs the same
  comparison;
- when `is_permutation_invariant` is true, the frame change is observationally irrelevant and
  `equiv_under` reduces to `equiv`.

This frame action is part of normalized equivalence rather than a separate comparison relation. It
records the dominant semantic case: most
molecular relation payloads do not assign values to participant positions. Dative-bond order,
noncovalent-bond kind, stereo configuration carried by an ordered ligand factor, and their
constraints are position-independent. In the current model, the only position-sensitive payload
fields are the per-participant electron counts of aromatic systems and multicenter bonds. Those two
implementations permute their electron-count vectors; an undetermined electron-count value is
itself permutation-invariant.

The graph-core traits nevertheless use a conservative default of `false` for
`is_permutation_invariant`, and every payload implementation must supply its position action. This
keeps a newly introduced position-sensitive payload from silently inheriting a no-op action. A
separate marker-trait hierarchy or a default no-op action would encode the common case with less
boilerplate but weaker review pressure and more public machinery.

`Equiv` and `BiEquiv` therefore remain derived, blanket-implemented operations. Payload types do not
override their comparison independently of `on_permutation` and normalized equality. Relation-set
remapping applies the position action to stored data, while read-only comparisons and matching may
apply it to a temporary value before using their own comparison relation. Graph core does not need
a second payload-equivalence API.

Applying a remapping to an independently supplied graph or relation set introduces a contextual
coverage condition: every participant must lie in the remapping's declared source domain. Public
relation-set APIs provide paired routes. `remap` asserts coverage for producer paths that
establish it and documents that a mismatch panics; `try_remap` checks coverage and returns
`None` for an independently supplied mismatch. Both use the same transport implementation. The
checked route belongs at application because construction cannot know which carrier will later be
supplied; `map_node`, `map_edge`, and participant-level remapping remain direct indexing operations.

A total mapping may target a sparse or larger ambient namespace, as when the rhs of a reaction is
embedded into an lhs-anchored union. Such a mapping can transport relation entries and referenced
values but cannot by itself produce a standalone `Molecule`, whose eight entity tables use dense
ids. End-to-end molecule remapping is defined only when every entity-family mapping is a bijection
onto a dense target id space. An embedding into a union and a remapping of a standalone molecule are
therefore related operations with different codomains.

### Dense molecule remapping

A public end-to-end remapping operation on `Molecule` accepts a `MoleculeCorrespondence` that
describes the complete old and new id spaces. The correspondence source counts must equal the
molecule counts, and every component correspondence must be total on both sides. The operation
returns `None` when these structural conditions do not hold or when the source molecule fails its
representation-integrity contract. The asserted route panics under either condition.

On success, it transports topology, every relation participant, position-sensitive relation data,
stereo frames, entity forms, and every typed reference in constraints. It does not validate
chemistry, resolve values, normalize attributes, repair references, compact tables, or remove
entities. Identity remapping is exact; applying a remapping and its inverse recovers the original;
and sequential remapping agrees with correspondence composition.

Dense molecule remapping is semantics-preserving alpha-renaming. Its primary semantic law is
`source.equiv_under(&remapped, &correspondence)`. Property tests must state this law directly in
addition to testing identity, inverse, composition, and referential integrity. Generated cases must
include crossing permutations, all entity families, position-sensitive relation data and stereo
frames, and constraints containing typed entity references; testing only reordered atom and bond
tables is insufficient.

Canonicalization is a consumer of this operation, not an alternative implementation of it. It
derives a complete correspondence from a canonical labeling and applies the ordinary molecule
remapping operation. Canonicalization code must not introduce a second path for transporting entity
tables or referenced values.

## Reaction-span construction

`ReactionSpan` stores the union-frame encoding of an actual span of two molecules. The union
namespace is necessary to share entity identity across the sides, but union-frame integrity alone is
not sufficient: the entries selected on the left and right must each form a referentially intact
`Molecule`.

For example, an `Unchanged` bond incident to a `Removed` atom is representable as an arbitrary
annotated union graph, but its right projection contains an edge without both endpoints. It is
therefore not a reaction span. Accepting it and dropping the bond during projection would impose a
cascading or SqPO-style transformation rather than faithfully represent the supplied entries.

This is a representation invariant of `ReactionSpan`, not chemical validation. A type named as a
span must contain both objects of the span. Consequently:

- `ReactionSpan::try_from_entries` checks the union namespace and the referential integrity of
  each projected side;
- `ReactionSpan::from_entries` asserts the same invariant for trusted producers;
- parsing accepts only entries that construct an actual two-sided span;
- `ReactionSpan::lhs`, `rhs`, `to_reaction`, and `correspondence` are infallible;
- projection retains every entity and constraint present on the selected side and remaps its
  references without repair or silent loss.

The union frame is lhs-anchored. Preserved entities retain their lhs ids and rhs-only entities are
appended. Consequently, `superimpose(lhs, rhs, correspondence).lhs() == lhs`, while the rhs
projection is equivalent to the supplied `rhs` under the induced total correspondence rather than
necessarily structurally equal to it. A crossing source correspondence changes rhs entity order;
neither a reaction's deltas nor the span membership columns encode that otherwise semantically
irrelevant permutation. `ReactionSpan::correspondence` relates the normalized projections and
does not retain a source correspondence whose rhs frame differs.

The two reaction representations deliberately have different construction boundaries.
`Reaction` is a permissive lhs-plus-deltas carrier and may contain deltas that cannot materialize
an actual two-sided span. `Reaction::to_reaction_span` is therefore fallible: it rejects deltas
whose projected product cannot form the second molecule of a span. The conversion retains its
existing `Contradiction` surface pending the repository-wide error review. Conversely, every
constructed `ReactionSpan` has a valid lhs and rhs, so `ReactionSpan::to_reaction` is
infallible.

### Span integrity and DPO application

There is no separate DPO-validity predicate on a reaction span. A removed atom with a surviving bond
or overlay would violate the right-side construction invariant: the purported `R` is not a molecule
and the value is not a span. The constructor's symmetric check also rejects an added atom required
by a surviving left-side entity and covers stereo and constraint references.

Remove the reaction-span validator entry point rather than retaining a validator that can only
confirm a type invariant. A check over permissive `Reaction` values that predicts whether their
deltas can materialize a span must be named and documented as reaction consistency or
materializability rather than DPO validation. The public name remains to be settled. DPO dangling
and identification conditions remain part of reaction application because they concern the
supplied host and match, not the reaction or reaction span alone.

This keeps fallibility at the representation boundary. `Reaction::reverse`, reaction
composition, and reaction fingerprints retain their operation-specific result surfaces; an
internal route through `Reaction::to_reaction_span` does not require side projection itself to
become fallible. Python checked construction reports invalid span entries as `ValueError`, while
`lhs`, `rhs`, and `to_reaction` return their values directly.

## Review questions

When adding or revising an API, ask in order:

1. Can the type store this value without unresolved references or malformed representation state?
2. If yes, is the rejected property instead a semantic predicate that should have a named
   validator?
3. Does the conversion preserve all information the target can represent?
4. Is any dropping, normalization, repair, or resolution being hidden inside construction or
   conversion?
5. Does each failure report the boundary that actually failed rather than a broader semantic label?

If the answer to question 1 is yes, a semantic failure alone is not a reason to reject construction.
