# Data type contracts

## Purpose

This is a normative developer guide for deciding which properties belong in constructors,
converters, integrity checks, validators, and transformations. It applies primarily to aggregate
model types such as `Molecule`, `Reaction`, and `ReactionSpan`. Small value types may
enforce stronger invariants when those invariants define the value represented by the type.

The central rule is that an aggregate constructor establishes that its input can be represented; it
does not establish every useful property of the resulting value. Semantic properties are lazy and
must be requested through named operations.

### Closed containers and the minimum eager contract

`Molecule`, `Reaction`, and `ReactionSpan` are closed containers. Every value obtainable through
their public construction, conversion, mutation, and transformation surfaces satisfies the type's
tier-1 representation-integrity contract. Their independently assembled entry and delta types are
open carriers: they may be incomplete or inconsistent until a checked or asserted aggregate
constructor accepts them.

The closed-container contract is deliberately the smallest eager validation surface that makes the
stored representation coherent and mutually interpretable. An operation on a closed value may rely
on that contract instead of repeatedly checking references, collection shapes, participant frames,
and other representation invariants inside the type's implementation. Every publisher must either
invoke the authoritative integrity gate or establish preservation of the complete contract by
construction.

Closed does not mean fully validated. A coherent value may still violate model-independent physical
invariants, fail a selected chemistry model, be unsatisfiable, fail to match a host, or violate an
operation-specific precondition. Those questions remain lazy and are checked by the named operation
that first requires them. Do not enlarge a closed container's eager contract merely to save a later
operation-specific check, and do not defer a representation invariant when doing so would force most
operations on the type to rediscover whether its stored fields have a coherent interpretation.

Use this tradeoff when defining other aggregate types: eagerly establish exactly the invariants
needed for all operations to interpret the value consistently, then leave properties needed by only
some operations to their explicit validation or execution boundaries. The current check-by-check
inventory and the concrete failure prevented by each check are maintained in
[Representation integrity](integrity.md).

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

`# Semantic properties` states the laws of an operation algebra — composition, roundtrips,
delegation, normal forms — while `# Assumes` and `# Establishes` state data properties:
predicates on values that downstream operations rely on. A property established by every public
producer of a type is stated once in the type's documentation rather than repeated in each
constructor's `# Establishes`.

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

The fixed relation semantics of every entity kind are part of molecule integrity. Each entity's
actual-atom participant relation is simple, even when storage divides it into distinguished
factors. This excludes localized and noncovalent self-loops, repeated dative participants across
the donors and acceptor, repeated aromatic and multicenter participants, repeated actual-atom
stereo ligands, and a stereo-atom site reused as an actual-atom ligand. Virtual ligands are not
actual-atom participants, but an identical virtual ligand occurrence may not repeat within one
stereo frame. An implicit hydrogen and a lone pair anchored at the same atom remain distinct
ligands.

Stereo ligand frames also carry site-relative incidence. An actual ligand of a stereo atom must be
a covalent graph neighbor of the site, while an implicit hydrogen or lone pair is borne by the site
itself. A stereo-bond frame consists of two consecutive two-ligand endpoint blocks. Each block's
actual ligands must be adjacent to its endpoint and cannot be the opposite site endpoint; its
virtual ligands must be borne by that endpoint. Exchanging the two complete endpoint blocks is
valid, but moving one ligand across the endpoint boundary is not.

The same contract supplies the cross-entity uniqueness needed to interpret each relation-valued
entity kind:
localized bonds have unique unordered endpoint pairs; dative bonds have unique complete
`(acceptor, donor multiset)` keys while distinct keys may share incidences; aromatic systems are
atom-disjoint; multicenter participant sets are unique; noncovalent bonds have unique unordered
endpoint pairs regardless of interaction kind; and stereo-atom and stereo-bond sites are unique
within their kinds. These are not deferred semantic judgments. They define the stored relation
represented by each entity kind and are established whenever a `Molecule` is published.

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

Representation integrity has one authoritative implementation in the crate that owns the graph IR.
It is the crate-private `check_integrity` operation with a corresponding public `*IntegrityError`;
there is no public validator or `*Checker` object. An integrity check returns
`Result<(), *IntegrityError>`, never `Solution`,
`Underdetermined`, or `Contradictory`. It includes stored entity and constraint references, parallel
collection shapes such as participant and electron-count lengths, fixed entity-relation semantics,
and kind-dependent data needed to interpret a value such as stereo frame arity, coset domains, and
permutation degree.

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
4. An operation accepting a published aggregate relies on the closed-container contract and does
   not invoke the check defensively. Trusted remapping, canonicalization, projection, reversal, and
   application publishers establish preservation by construction and property tests.

For `Molecule`, this places the same gate behind direct checked entry construction, molecule-DSL
raise, TableIR raise, Python construction from entries, and finalization of builders and editors.
The public aromatic-system and multicenter-bond mutation operations use a private candidate and
commit only after the same gate succeeds; their raw whole-form mutation kernels remain crate-private.
Reaction code uses the molecule constructor whenever it materializes a molecule side or projection;
it does not grow a reaction-specific copy of molecule integrity. A route that already establishes
the contract may use the asserted constructor, but the asserted and checked routes must share the
implementation and the set of accepted values.

Do not run a semantic validator at these boundaries. DSL and external-format raise may preserve an
atom `#T` or bond `#C` assertion before a stereo entity has been perceived or resolved. Validating
constraint satisfaction during construction would reject a representable input and silently impose
an operation order. Source-format checks such as wedge, directional-bond, or chirality-frame
interpretation likewise stay in the source-format layer; only the resulting graph IR's representation
is checked by the shared integrity gate.

Those assertions must eventually be realized against an integrity-valid stereo frame. They do not
create distinguishable occurrences of an equal virtual ligand. If perception could satisfy `#T` or
`#C` only by repeating an implicit hydrogen or lone pair, it perceives no such stereo entity and
follows the operation's existing absence policy. Any boundary that supplies an explicit repeated or
oversized frame instead rejects it during raise; it must not normalize, deduplicate, or choose an
arbitrary coset for malformed input.

Canonicalization is a transformation of a closed, representation-integrity-valid value. It does not
repair or revalidate malformed representation state, and its error types contain no unreachable
integrity arm. Canonicalization preserves tier-2 and tier-3 invalid states; intrinsic normalization
or reaction-span materialization may still report `Contradiction`.

## Representation ownership and crate layering

An operation belongs with the crate that owns the representation when its correctness requires
coordinated access to that representation's complete internal shape. For `Molecule` and the
reaction representations, this includes operations that must keep all entity sets, participant frames,
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
- entity relations satisfy their fixed participant and cross-entity uniqueness semantics;
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

Flat aggregate entries preserve each entity kind's relation structure. A unipartite relation has
one participant collection paired with its attributes: variable-degree aromatic systems and
multicenter bonds use `Vec<AtomId>`, while fixed-degree noncovalent bonds use `[AtomId; 2]`.
Bipartite or site-bearing entities keep their distinguished factors separate, as in dative bonds
and stereo entities. Do not flatten a unipartite fixed relation into separate participant tuple
fields merely because its degree is two.

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

An equivalent `Modified` reaction-span entry is another raw/normal distinction. The tag carries no
assertion beyond its two side values. If those values are equal under `normalized_eq`, the entry is
semantically a no-op and normalizes to `Unchanged`. Checked and asserted span construction preserve
the explicitly supplied raw tag. `ReactionSpan::normalize` collapses it, and reframing and
canonicalization inherit the collapse because normalization is their first pipeline step.

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
non-normal value encodings. `normalized_eq` compares normal forms in the current participant and
entity-id frame. `framed_eq` additionally selects participant frames, while
`framed_eq_under` first applies an explicitly supplied entity-id remapping and then performs
framed equality. Aggregate `canonical_eq` compares complete canonical IR values after selecting
participant frames and entity ids. Every semantic comparison includes constraints. This
distinction matters for patterns, where constraints are not redundant with the structural
description.

Entity and molecular structural identity is established from inherent fields and structural
incidence. Constraints restrict the states admitted by an entity or molecule but do not establish
that structural identity. They therefore do not distinguish the initial structural automorphism
orbits. After structural labeling establishes the minimum structure key, aggregate canonicalization
normalizes entity-level and molecule-level constraints and minimizes the typed canonical constraint
key over every admissible entity remapping and participant-frame action that attains that structure
minimum. Equivalently, complete canonicalization minimizes `(structure key, constraint key)`
lexicographically: constraints distinguish only tied structural frames and never make a larger
structure key win. Remapping transports every reference they contain. Graph-only automorphism
orbits may discard a tied frame only when the orbit action carries the complete covariant
participant-frame and constraint action; otherwise those frames are searched exactly. This
post-hoc participation gives the complete IR assertion a unique canonical form without turning
constraints into structural identity features. Complete canonicalization must therefore retain or
reconstruct every action attaining the minimum structure key until the constraint minimum is
selected; it cannot commit the single `Structure` representative and inspect only that frame.

## Aggregate canonicalization

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

The equality operations form the same nested quotient pipeline:

- `==` compares the exact stored representation;
- `normalized_eq` compares normal forms in the current participant and entity-id frame;
- `framed_eq` compares after selecting participant frames, while `framed_eq_under` first applies an
  explicitly supplied entity-id remapping and then performs framed equality; and
- `canonical_eq` compares complete aggregate canonical forms under a shared context, selecting the
  participant frames and entity ids rather than receiving an id witness from the caller.

For integrity-valid inputs whose complete canonicalization succeeds, `canonical_eq` holds exactly
when an admissible total dense correspondence exists under which `framed_eq_under` holds. Two
intrinsic contradictions still compare equal under canonical equality's failure-totalization rule,
but that convention does not require a correspondence witness. `framed_eq_under` is the explicit
entity-id-witness form of `framed_eq`, not another quotient level.

The context-bearing trait has the semantic shape

```rust
pub trait Canonicalize: Reframe {
    type Error;

    fn canonicalize(
        self,
        context: &CanonicalizeContext,
    ) -> Result<Self, Self::Error>;

    fn tracked_canonicalize(
        self,
        context: &CanonicalizeContext,
    ) -> Result<(Self, MoleculeRemapping), Self::Error>;

    fn canonical_hash(
        self,
        context: &CanonicalizeContext,
    ) -> Result<u64, Self::Error>;

    fn canonical_eq(
        &self,
        other: &Self,
        context: &CanonicalizeContext,
    ) -> bool;
}
```

`tracked_canonicalize` returns the entity-id-renumbering witness, not the participant
frame action. Applying that remapping and then `reframe` reconstructs the canonical value;
remapping alone need not do so.

Use one concrete context for `Molecule`, `ReactionSpan`, and `Reaction` unless an
implementation establishes a real need for distinct context types. Canonical-form construction is
fallible; equality is total. Structurally equal inputs compare equal immediately. Two successfully
canonicalized inputs compare by structural equality, and two intrinsic contradictions compare
equal because both denote the empty semantic value. One contradiction and one successful form do
not compare equal. Integrity failures never make distinct inputs equal; callers that need the
diagnostic invoke `canonicalize` directly.

The public operation is complete-only. Callers cannot select topology, constitution, or structure
as a reduced comparison surface, and a molecule does not publish a description-level query.
Canonicalization may use a private effective structural prefix to avoid constructing empty search
domains, but every present entity, inline constraint, and molecule-level constraint participates in
the result and in `canonical_eq`.

Stored stereo entities participate whether or not para-stereo refinement is enabled. Without
para-stereo refinement, perform one stereo-sensitive refinement from the constitution-level
partition. With it enabled, feed each stereo-sensitive partition back into refinement until the
partition stabilizes. The refinement state is the ordered partition of the complete structure-level
incidence graph. Each round recomputes stereo descriptors from the current classes of the stereo
entity's site and ligand occurrences, refines every incidence node by its previous class together
with its applicable stereo descriptor, and then applies ordinary exact graph refinement. Including
the previous class makes refinement split-only. Stability means equality of the ordered partition
cells, not merely an unchanged cell count or unchanged backend color integers. Every non-stable
round strictly increases the number of cells, so the finite incidence-node set supplies the
termination bound. Do not expose an iteration cutoff that can change the canonical form.

Backend canonical labels are search inputs, not the numbering authority. The canonical frame
is the minimum under the library's typed comparison order. Automorphism generators and orbits may
prune equivalent branches, but changing the selected graph algorithm must not change the resulting
canonical representative.

For a fixed umol release, canonicalization is deterministic under a fixed context. The
returned canonical form is an ordinary IR value: it carries no canonicalization-schema version or
producer provenance and is not a persistent identifier. During the 0.x series, the typed comparison
schema and resulting canonical numbering may change between releases as the entity model and
canonicalization rules are corrected or extended. Persisted canonical forms must therefore record
their producing umol version externally and must not be compared across releases as stable ids.

The explicit positions below define the current implementation and keep its order independent of
Rust enum declaration order. An append-only extension that leaves every earlier value unchanged is
desirable when the model permits it, but it is not a 0.x compatibility guarantee. New entity kinds,
fields, and constraint variants—and corrections to existing ones—must choose their positions
deliberately and update the normative table and exact ordering tests together. A future durable
canonicalization profile requires its own explicit version, compatibility rules, and conformance
fixtures; the unversioned aggregate API does not imply that contract.

### Canonical comparison schema

The following schema is the normative typed order for the current aggregate-canonicalization
implementation. Numeric positions are local to the table in which they occur. Code must use these
explicit positions instead of inferring order from Rust declarations. During 0.x, an intentional
schema revision may change them together with this table and the corresponding exact tests.

The entity model has three ordered structural domains. Topology is AB, non-stereo is DAMN, and
stereo is SS. Constitution is topology plus non-stereo. Overlays are non-stereo plus stereo, and
structure is topology plus overlays. The private canonicalization implementation uses those
prefixes to select the least search carrier containing the complete input; this is an optimization,
not a public description-level operation. Complete canonicalization appends normalized entity-level
and molecule-level constraints after the complete structure key.

| Domain position | Structural domain | Entity indices |
| ---: | --- | --- |
| 0 | Topology | Atom = 0, Bond = 1 |
| 1 | NonStereo | Dative bond = 0, Aromatic system = 1, Multicenter bond = 2, Noncovalent bond = 3 |
| 2 | Stereo | Stereo atom = 0, Stereo bond = 1 |

An entity-block position is the composite `(domain position, entity index)`. Entity blocks compare by
this position and then by their dense row sequence. This hierarchy makes topology, constitution,
and structure exact internal domain prefixes while allowing a future entity kind to be appended
within the domain to which it belongs. An absent block contributes no entry. The separate terminal
constraint section follows every structural domain. Extending an entity domain does not move that
section or alter keys for molecules that lack the new kind.

Rows compare their components in the following local field order. Every present entity kind
contributes all of its inherent fields; participant topology,
participant-indexed values, and frame-dependent values occupy the listed structural components
rather than being omitted. The dense row index is implicit in the row sequence. Inline constraints
are excluded here and enter through the constraint section.

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

Canonical-search initial colors retain these listed field positions even when participant data
is represented by incidence occurrences. Omitted participant-bearing fields are not renumbered: for
example, a bond node uses positions 1 through 3 for order, charge, and unpaired electrons, while its
endpoint pair is represented by two incidences. Entity-node colors contain only normalized,
constraint-free, frame-independent values. Aromatic and multicenter participant electron counts and
stereo ligand kinds occur on their corresponding incidences; raw stereo configurations do not enter
the initial node colors.

Typed incidences use the following order:

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
normalized incidence values. Entity-node and incidence colors occupy disjoint key domains.

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

`RelationalConstraint` uses the composite position `(owning entity-block position, local index)`.
The assigned local indices are:

| Owning entity | Local indices and variants in order |
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
constraint, or entity kind within its assigned domain can leave the key of every earlier value
byte-for-byte equivalent at the typed-key level. This is a useful extension property, not a
cross-release promise. A genuinely new structural category may instead require another domain;
moving an existing entity kind between domains or inserting a domain changes both comparison order
and the private structural-prefix semantics. The concrete Rust key storage remains private and may
change.
Exact ordering tests instantiate the current positions to detect accidental drift; an intentional
schema revision updates those tests and this table together.

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

Reaction iterators are operation-issued values with a different lifecycle. `Reaction::apply`
checks reaction-wide preconditions and then issues a `ReactionApplicationIter<T>` that owns snapshots
of the reaction and host, normalized application state, and an eagerly enumerated correspondence
set. `React::react` issues a `ReactionProductsIter` over the same application lifecycle. Neither
iterator has a public constructor. The output type is selected by `apply` (molecule),
`tracked_apply` (molecule and correspondence), `apply_to_reaction`, or `apply_to_reaction_span`.
Outputs and product-component lists are realized lazily in
match order, so failures at that stage remain iterator items; a fatal item error is yielded once and
terminates the iterator. Later mutation of the source reaction or molecules cannot change the
issued operation.

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
to derive the remaining entity kinds are therefore operation preconditions, not defenses against
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

For example, converting a partial correspondence to a semantic remapping may return `None` because
the correspondence is not total on both sides. This can occur for a correspondence correctly
produced for its molecule pair: partiality is part of the correspondence model, whereas a
remapping is a total bijection. Such narrowing is not currently exposed; reference transport
consumes correspondence directly.
It does not make correspondence construction, composition, reversal, or unrelated consumers
fallible. Exact error taxonomy remains subject to the repository-wide error review; the
construction/validation boundary does not require introducing a new error type for each method.

### Supplied-match reaction application

The four supplied-match reaction methods return `Result<Option<T>, ApplyError>`:
`apply_at` produces a molecule, `tracked_apply_at` produces a molecule and its
host-to-product correspondence, and `apply_at_to_reaction` /
`apply_at_to_reaction_span` produce the realized primary objects.
`Ok(None)` is ordinary non-applicability caused by dangling incidence or a structural
conflict; it is not an execution error. Iteration skips that outcome directly.
Actual errors are yielded once and terminate iteration; error variants are not
classified to recover continuation policy.

All four methods preserve the same application and applicability behavior. The supplied
correspondence maps the rule into the host; the produced witness maps host to product.
Other entity pairings are induced from atom provenance and compatible incidence. Equal
stereo kinds remain matchable, while different determined kinds identify different entities
and produce removal plus addition. An undetermined kind does not assert a conflicting geometry.
Realized reactions and spans use the existing lhs-anchored representation and constraint
normal form; their product projections preserve the produced molecule's semantics rather
than promising its exact row ordering.

### Consuming correspondence updates

`Correspondence::identity(count)` initializes a declared identity without temporary images.
`extend_right` consumes the correspondence and adds unmatched right-domain ids by increasing
its right count. `compact_right` discards affected pairs and compacts surviving right ids;
`uncompact_right` expands those ids through the inverse compaction, leaving restored positions
unmatched. Both check the applicable intermediate count using the existing composition errors.
They preserve the left count and reuse the pair-vector allocation, including when all pairs
are removed. Compaction followed by expansion does not recreate discarded pairings.

Graph and molecule correspondences expose component-wise versions. Graph extension takes node
and edge increments; molecule extension selects an `EntityKind`. Molecule compaction adapts only
the removed atom/bond id lists from graph ids, not the vectors of matched pairs. These operations
do not clone molecular payloads or require intermediate molecules. Ordinary `compose` remains
borrowed and returns a separate correspondence.

### Editor session correspondence

`MoleculeEditor` accumulates an initial-to-current correspondence over all eight entity kinds.
Tracking retains only id pairs and counts. Additions extend the right domains, removals compact
them, and undo restoration applies inverse compaction without recreating discarded pairs.
Attribute changes preserve pairings. A failed transaction whose rollback succeeds restores the
pre-transaction session correspondence, alongside the molecular state.

`tracked_snapshot` returns the same integrity-checked molecule as `snapshot`, plus an independent
copy of the session correspondence. `try_tracked_build` and `tracked_build` transfer the accumulated
vectors into the result and share the ordinary checked/asserted publication boundary. Plain
publication returns only the molecule; plain snapshots do not copy the correspondence.

## Pushout results

Graph and relation-set `pushout` methods return only the resulting object.
Their `tracked_pushout` companions return that object and its two input-to-result
mappings: `(Graph, GraphPushoutCorrespondence)` and
`Option<(Self, RelationPushoutCorrespondence)>`, respectively. Each pair has identical
outputs and failure behavior; mapping carriers are optional API output.
The carriers contain public `left` and `right` fields, not the result object; no constructors
are added.
Each operation produces components with equal target counts covering their respective
inputs. Independently assembled fields do not establish agreement with a result object.
Relation payload-combination failure still returns `None`.

Molecule `meet_pushout` returns `Option<Molecule>`; `tracked_meet_pushout` returns
`Option<(Molecule, MoleculePushoutCorrespondence)>`. Its carrier likewise contains only public
`left` and `right` correspondences, with no constructor. Operation-produced components map each
input into the result and share its target counts for all eight entity kinds. Both forms preserve
the same attribute meets, constraint transport, and inadmissibility behavior.

## Pullback and pushout-complement results

Graph and relation-set `pullback` methods return the object; `tracked_pullback`
also returns `PullbackCorrespondence` or `RelationPullbackCorrespondence`.
Their public `left` and `right` components map the result into each input and
have equal source counts when produced by the operation. Graph pullback retains
its count-mismatch `Result`; relation-set pullback retains `Option` for rejected
payload combinations.

Graph `pushout_complement` returns the context graph; `tracked_pushout_complement`
also returns `PushoutComplementCorrespondence`. Its public `context` maps context
to host, and `interface` maps interface to context. Both forms retain `None` for
dangling or count disagreement. Each plain/tracked pair has identical results
and failure behavior. These carriers contain no result object or new constructors;
agreement of independently supplied fields with their graphs remains contextual.

## Compaction

`Compaction<Id>` declares a finite source count and stores sorted, distinct removed ids.
Construction rejects out-of-range removals with `CompactionError`; the result count is the
source count minus the number removed. Identity requires an explicit source count, not `Default`.
`empty()` on compaction or remapping types means zero source and result counts in every component;
it is not a domain-independent identity. `Remapping::identity(count)` and
`GraphRemapping::identity(node_count, edge_count)` construct declared-domain identities. Existing remapping defaults remain zero-domain values.
Graph and molecule aggregates accept validated component compactions, including count-bearing
identities for untouched entity kinds. Producers capture counts before removal; `UndoCompaction`
preserves them.

`compact` returns `None` for removed or out-of-source-range ids. `uncompact` asserts result-domain
membership; `try_uncompact` returns `None` outside that domain. `compact_vec` requires exactly the
source count of values, and `try_compact_vec` returns `None` on a length mismatch. Survivor order
is preserved. Borrowed conversions to single-space, graph, and molecule correspondences preserve
both counts and every survivor pairing in each entity kind.

Graph `remove_cascading` returns `()`, and `try_remove` returns `Option<()>`;
`tracked_remove_cascading` and `try_tracked_remove` return the graph compaction instead.
The checked removal leaves the graph unchanged when the dangling condition fails.
Single-node/edge cascading-removal conveniences return no witness. Relation-set
`compact` returns the resulting set; `tracked_compact` returns that set
and its relation-id compaction. Each pair preserves the same output and failure behavior.

Molecule editor `remove` and the six bulk relation-removal methods return `()`; their
`tracked_` companions return a full `MoleculeCompaction`, including unchanged-family counts.
Molecule `extract` returns only the molecule; `tracked_extract` pairs it with the actual
host-to-result compaction. Its sub-to-host input describes selection, not result numbering:
extraction preserves host order. Both forms retain the same cascades and constraint transport.

Rollback checks result counts against the current editor before applying an inverse compaction.
A mismatch remains `TransactionError::RollbackStateMismatch`, not an indexing panic.

## Remapping

Remapping is an explicit total bijection between dense id spaces. It transports represented values
and incidence without adding, dropping, or identifying entities; it does not validate chemistry,
normalize attributes, repair references, or remove entities. Removal uses compaction instead.

`Remapping::new` checks that the image vector is a permutation of `0..images.len()` and returns
`RemappingError` for the first out-of-range or repeated image. It does not modify the images.
The private vector has no public mutable access. Graph and molecule aggregate constructors take
already-valid component remappings and do not repeat validation. Borrowed `From` conversions
widen remappings to correspondences without changing pairings or source and target counts.
Agreement with an independently supplied object remains contextual.

The private storage is `Vec<Id>`. Construction and lookup require `Copy + Into<usize>`;
lookup accepts and returns the same `Id` type, with integer extraction confined to storage access.
Operations that enumerate source ids additionally require `From<usize>`. These conversions do not
bind an id to a particular object or establish its validity in that object's namespace.

The facility has three coordinated levels:

- `umol_graph_core::Remapping<Id>` stores the dense image vector for one typed id space;
- `GraphRemapping` aggregates node and edge remappings used by graphs and relation participants; and
- `MoleculeRemapping` aggregates remappings for all eight molecule entity kinds.

A relation-set remapping relabels each factor and leaves both the participant sequence and the
payload as supplied. Graph core never reorders a frame or reads a payload; a positional payload
stays aligned because nothing moved.

Constraint reference transport consumes `MoleculeCorrespondence`. `try_map` requires an image for
every referenced entity and returns `None` if any is absent; `map` asserts that same coverage.
Unused correspondence entries need not be matched. Predicates and frame positions are preserved.
Edit handle resolution is separate: it uses only the entities referenced by the edited constraint,
not a molecule-wide correspondence.

When a higher-level operation moves molecular data with coordinated graph-core and graph-IR maps,
both must derive from the same operation witness: the graph-core mapping transports topology and
relation participants, while the IR mapping transports references to owned entity rows. Do not
manually sort remapped relation participants or permute their payloads at individual call sites;
that belongs to graph IR, which owns frame selection and transport.

### Participant frames and payload equivalence

A relation set stores the participant sequence it is given. The participant multiset is the
relation's identity; the stored sequence is the coordinate frame its payload is expressed in. Graph
core never interprets that frame — it relabels ids and preserves sequence — so all frame semantics
live in graph IR.

Published entity frames contain pairwise-distinct complete participant values. In particular, a
stereo frame cannot repeat an equal `StereoLigand`, including an implicit hydrogen or lone pair with
the same anchor and kind, and its length cannot exceed `umol_perm::MAX_DEGREE`. These integrity
rules make the action between two compatible stored frames unique and keep every stereo action in
the bounded `Permutation` representation.

Comparing two entries therefore has three independent parts, and a site chooses each:

- **identity** — do the two hold the same structured participants under the entity kind's factor
  semantics. Each ordinary factor is compared as a multiset; stereo-bond endpoint blocks may be
  exchanged only as complete blocks. Graph-core coincidence supplies ordinary factor comparison
  through `coincident`, `coincident_edge`, `is_coincident`, and `participants_match`; graph IR adds
  the entity-kind structure that storage alone cannot express.
- **frame transport** — restating one side's payload in the other's frame. Every entity frame has
  distinct complete participant values, so equal structured incidence determines one action in the
  entity kind's group. `DynPermutation::between` or `Permutation::between` derives that action, and
  `FrameTransport::reframe_by` transports the payload.
- **the value relation** — `normalized_eq`, `matches`, or `meet`, the caller's own.

`FrameTransport::reframe_by` returning `Some` does not establish identity: a frame-invariant payload
reads neither frame, and an undetermined electron-count vector has nothing to reorder. A site that
needs both must ask for both.

The six overlay kinds select and transport frames as follows:

| Entity kind | Representative frame | Local action | Position-sensitive values |
| --- | --- | --- | --- |
| Dative bond | donors sorted by atom id; acceptor fixed | `DynPermutation` on donors | donor sequence; the current form and constraints are invariant |
| Aromatic system | participants sorted by atom id | `DynPermutation` | participant electron counts |
| Multicenter bond | participants sorted by atom id | `DynPermutation` | participant electron counts |
| Noncovalent bond | endpoints sorted by atom id | degree-two `DynPermutation` | ordered predicate pairs in `NoncovalentBondEndsSatisfy`; the entity form is invariant |
| Stereo atom | complete ligands sorted | bounded `Permutation` in the full ligand symmetric group | configuration, ligand-symmetry and fluxionality permutations, topicity positions, and stereo constraints |
| Stereo bond | sort each two-ligand endpoint block, then order the blocks | bounded `Permutation` in `S_2 wr S_2` | configuration and every frame-relative stereo-bond constraint |

The stereo-bond group may permute within each endpoint block and exchange the two complete blocks;
it cannot move one ligand across the block boundary. Dative-bond order, noncovalent-bond kind,
aromatic and multicenter charge and spin, and their inline constraints are currently
frame-invariant. Their `FrameTransport` implementations destructure the complete form, so a new
field must be classified explicitly rather than being silently left unframed.

`FrameTransport` is the transport-only operation for a receiver and an independently supplied
action. Entity forms, `EntitySpan<form>`, individual overlay delta payloads, and constraint values
implement it because they can consume a compatible action but do not own enough participant data
to select one. `Reframe` is implemented only by a frame-owning carrier. It extends `Normalize`,
derives a representative action, and therefore represents the second prefix in the
normalize–reframe–canonicalize pipeline.

The associated action is complete for its receiver:

| Carrier | `FrameTransport::Action` |
| --- | --- |
| One ordinary overlay form or form span | one `DynPermutation` |
| One stereo form or form span | one bounded `Permutation` |
| One entity-kind aggregate or `*Spans` peer | one typed local action per entity id |
| `Molecule` or `ReactionSpan` | one six-component `OverlaysFrameAction` |
| `Reaction` | one `OverlaysFrameAction` covering every lhs and `Add`-owned overlay id |
| `Constraint`, `ConstraintSpan`, or `ConstraintDelta` | a covering `OverlaysFrameAction` |
| `Deltas` without its reaction | none; removals may use local frames whose owners live on the reaction |

Each entity-kind aggregate action and `OverlaysFrameAction` is an operation-issued witness with
private construction. Identity and inverse preserve its exact typed id-and-degree domain, and
composition is defined only for equal domains and degrees. Consumer compatibility is weaker and
receiver-relative: a receiver requires coverage for every frame-relative value it contains, ignores
irrelevant entries, and returns `None` for a missing degree, inadmissible subgroup action, or other
observable incompatibility. A witness may therefore be reused with another compatible carrier; it
is not bound to one object identity.

`representative_action` is derived from the input's frame owners before normalization and remains
total when later normalization finds an intrinsic contradiction or erases an entry.
`tracked_reframe` returns that input-domain witness with the result. Plain aggregate `reframe`
derives and immediately consumes local actions as it visits entries; it does not pre-emptively
allocate the complete aggregate witness merely to discard it. Sparse action storage used for
frame-relative constraints allocates no backing map for an empty domain.

Each constraint form classifies its own participant-frame use through an exhaustive match. Each
entity delta likewise classifies its fields exhaustively and delegates constraint payloads to that
form-level decision. Aggregate domain collection and application delegate to those classifications
rather than repeating variant lists, so adding a field, constraint form, or delta variant cannot be
silently omitted from frame transport.

A reaction removal may carry the source entity's structured incidence in another participant
ordering. Its repeated sequence is an explicit local frame for the recorded payload. Construction
preserves that frame; a consumer derives the unique local-to-source action, transports the payload,
and then applies its value relation. For generic aggregate frame transport, let `q` be that
local-to-owner action and `a` the action on the owning lhs or `Add` frame. The removal consumes the
conjugate `q.compose(a).compose(q.inverse())`; this preserves its local-to-owner relation and the
identity, inverse, and composition laws. Reaction normalization instead applies `q` directly to
align the removal with its owner. Reframing normalizes first, so `q` is identity when it applies the
selected owning action. Application-specific rule-to-host alignment may likewise move an already
normalized removal directly into the host frame; that is not generic frame transport of a raw
locally framed reaction. An incompatible distinguished factor or stereo-bond endpoint-block
assignment changes structured incidence rather than producing a second frame of the same entity.
Reaction integrity reports that case as `IncidenceMismatch`.

Reaction application derives the unique action from each mapped rule frame to the host frame and
transports every frame-relative field and constraint delta before matching. After the pattern
`old` value matches, lowering records the realized host value as the transaction's `old`; it does
not retain the rule pattern as if it were the concrete value being replaced.

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
ids. End-to-end molecule remapping is defined only when every entity-kind mapping is a bijection
onto a dense target id space. An embedding into a union and a remapping of a standalone molecule are
therefore related operations with different codomains.

### Dense molecule remapping

A public end-to-end remapping operation on `Molecule` accepts a `MoleculeRemapping`.
Each component length must equal the corresponding molecule entity count; bijectivity is already
established by construction. The checked operation returns `None` on a count mismatch, and the
asserted route panics under the same condition. Component `remap_vec` operations move table values;
constraint transport widens the remapping to a correspondence and uses the existing traversal.

On success, it transports topology, every relation participant, position-sensitive relation data,
stereo frames, entity forms, and every typed reference in constraints. It does not validate
chemistry, resolve values, normalize attributes, repair references, compact tables, or remove
entities. Identity remapping is exact; applying a remapping and its inverse recovers the original;
and sequential remapping agrees with correspondence composition.

Dense molecule remapping is semantics-preserving alpha-renaming. Its primary semantic law is
`source.framed_eq_under(&remapped, &remapping)`. Property tests must state this law directly in
addition to testing identity, inverse, composition, and referential integrity. Generated cases must
include crossing permutations, all entity kinds, position-sensitive relation data and stereo
frames, and constraints containing typed entity references; testing only reordered atom and bond
tables is insufficient.

Canonicalization is a consumer of this operation, not an alternative implementation of it. It
constructs a complete remapping directly from its selected entity order and applies the ordinary molecule
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
- both entry constructors preserve an explicitly supplied `Modified` tag, including one whose
  sides are `normalized_eq`;
- parsing accepts only entries that construct an actual two-sided span;
- `ReactionSpan::lhs`, `rhs`, `to_reaction`, and `correspondence` are infallible;
- projection retains every entity and constraint present on the selected side and remaps its
  references without repair or silent loss.

Construction permits any dense union-frame ordering whose two projections satisfy the integrity
contract. LHS anchoring is a normal-form property, not an additional construction invariant. In the
LHS-anchored form, entries present on the LHS form each entity kind's dense prefix and right-only
entries are appended. `ReactionSpan::to_reaction` reanchors an arbitrary valid union order into this
form by assigning compact projected LHS ids first and fresh ids to additions. Aggregate
canonicalization likewise always emits the LHS-anchored form because it carries no less information
and corresponds directly to the LHS-plus-deltas representation of `Reaction`.

`ReactionSpan::superimpose` already produces the LHS-anchored form: preserved entities retain their
lhs ids and rhs-only entities are appended. Consequently,
`superimpose(lhs, rhs, correspondence).lhs() == lhs`, while the rhs projection is equivalent to the
supplied `rhs` under the induced total correspondence rather than necessarily structurally equal to
it. A crossing source correspondence changes rhs entity order; neither a reaction's deltas nor the
span membership columns encode that otherwise semantically irrelevant permutation.
For a paired entity whose values are `normalized_eq`, `superimpose` emits `Unchanged` carrying the
lhs value. This is permitted because `superimpose` derives a span rather than faithfully converting
open entries: the choice retains exact lhs projection and the same semantic rhs-projection
guarantee without eagerly normalizing the selected payload.
`ReactionSpan::correspondence` relates the normalized projections and does not retain a source
correspondence whose rhs frame differs. A span-to-reaction-to-span roundtrip is therefore exact only
for an LHS-anchored span in the documented constraint normal form; for another valid union order it
returns the equivalent LHS-anchored normal form.

Dense span remapping follows the molecule remapping API. `ReactionSpan::remap` and
`ReactionSpan::try_remap` accept a `MoleculeRemapping` over the eight union tables. Component
lengths must equal the span's table sizes; each component is already a dense bijection. The operation transports both values of every `EntitySpan`, relation participants,
position-sensitive relation data, stereo frames, and all constraint references. It does not
normalize, repair, project, or reanchor the span implicitly. The asserted route panics when its
documented producer contract is violated; the checked route returns `None` only for a component-count mismatch. The published source and remapped result are closed by construction and are not
rechecked.

Canonicalization constructs an LHS-anchored remapping and applies this ordinary remapping
operation. It does not maintain a second canonicalization-only transport path. For
`S(s) = s.to_reaction().to_reaction_span()?` and span canonicalization `C`, the normal-form laws are
`S(C(s)) == C(s)`, `S(S(s)) == S(s)`, and `C(S(s)) == C(s)`. Only the first requires exact
roundtripping of the canonicalized input; arbitrary valid union ordering is normalized rather than
preserved.

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
