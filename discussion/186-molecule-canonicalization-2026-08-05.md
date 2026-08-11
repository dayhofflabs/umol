# 186 - Molecule canonicalization

Status: In Progress
Date: 2026-08-05
Relates: [156](156-ast-comparison-and-property-suite-2026-07-20.md),
[113](113-ast-canonical-equality-and-lattice-2026-06-14.md),
[117](117-entity-model-extensibility-2026-06-20.md),
[185](185-python-reaction-span-2026-08-04.md),
[188](188-inchi-migration-2026-08-09.md)

`Molecule` does not provide a canonical representation with respect to entity renumbering.
`Graph::canonical_key` and canonical labels are available in
`umol-graph-core::algorithms::automorphism`; the molecule-level incidence (Levi) graph is available
in `umol-graph-ir::ir::incidence`. Canonicalization must derive a canonical frame and apply it
through the ordinary end-to-end molecule-remapping operation. Consult doc 156 for the relation
between `equiv`, frame changes, and canonical equality.

## Normalization, canonicalization, and equality

The existing `Canonicalize` trait actually performs intrinsic, fixed-frame form normalization. It
folds lattice expressions, fields, constraints, entity forms, and `Deltas` without changing entity
ids or participant frames. It is context-free and returns `Contradiction` when the represented value
is unsatisfiable. This operation will be renamed to `Normalize`; its by-value operation becomes
`normalize`, its borrowed normal-form projection becomes `normalized`, and `Canonical<T>` becomes
`Normalized<T>`. `Lattice` depends on `Normalize`, not on aggregate canonicalization.

Canonical representative selection for `Molecule`, `Reaction`, and `ReactionSpan` is the
operation that retains the conventional name `Canonicalize`. It selects an entity-id and participant
frame modulo structural isomorphism, uses a graph canonical-labeling algorithm, and depends on an
explicit canonicalization context. These aggregate types do not implement `Lattice` and do not
acquire a lattice or a context-free normalization requirement merely to share this trait. Aggregate
canonicalization normalizes every carried value after transporting it into the selected frame.

The equality API follows the same separation:

- `==` compares the exact stored representation;
- `equiv` compares normalized values in the current id and participant frame;
- `equiv_under` first transports the receiver through an explicitly supplied participant order or
  entity correspondence and then applies `equiv`; and
- `canonical_eq` compares complete aggregate canonical forms after the implementation has selected
  the canonical frame.

Thus `equiv` and `equiv_under` are one semantic-equivalence family, not two additional equality
relations. Aggregate `canonical_eq` is the search-based extension of `equiv_under`: for inputs in
the operation domain, it holds exactly when canonicalization produces the same complete IR form,
and equivalently when an admissible remapping exists under which `equiv_under` holds. The current
fixed-frame `Canonicalize::canonical_eq` becomes `Equiv::equiv`; no fourth equality relation is
introduced.

The aggregate trait has the following intended shape:

```rust
pub trait Canonicalize: Sized {
    type Error;

    fn canonicalize(
        self,
        context: &CanonicalizationContext,
    ) -> Result<Self, Self::Error>;

    fn canonical_eq(
        &self,
        other: &Self,
        context: &CanonicalizationContext,
    ) -> bool;
}
```

The context is concrete rather than an associated type unless the implementations establish a
real need for different context types. `Molecule`, `ReactionSpan`, and `Reaction` share the
same canonicalization semantics and context. Canonical-form construction is fallible; equality is
not. A caller that needs a diagnostic invokes `canonicalize` directly rather than using equality as
an error-reporting operation.

Normalized and canonical equality totalize failures according to their meaning. Structurally equal
inputs compare equal immediately. Two successful normal or canonical forms compare structurally.
Two intrinsic `Contradiction` results compare equal because both denote the empty semantic value;
one contradiction and one successful form do not. Integrity or operational failures do not make
distinct inputs equal. Each aggregate implements this distinction directly rather than relying on a
generic comparison of its error values.

The aggregate error types are operation- and carrier-specific:
`MoleculeCanonicalizationError`, `ReactionSpanCanonicalizationError`, and
`ReactionCanonicalizationError`. The molecule error contains molecule integrity and intrinsic
contradiction failures. The reaction and span errors preserve their own integrity and conversion
boundaries rather than exposing an accidental internal call path. Their exact variants are fixed
after the integrity-error migration establishes those underlying error types.

`Reaction` currently implements the old `Canonicalize` by normalizing only its deltas while
preserving the LHS frame. That partial aggregate implementation conflates the two domains and must
be removed during the rename. `Deltas` continues to implement `Normalize`; callers that require
fixed-frame delta normalization invoke it explicitly.

As an adjacent nomenclature migration, rename the associated type `Ctx` to `Context` on `FromIr`,
`IntoIr`, `TryFromIr`, and `TryIntoIr`, and use `context` in their public trait signatures. This
does not belong to the normalization/canonicalization semantics, but it shares the terminology and
should land with the trait-surface cleanup. Public identifiers and associated types use complete
words rather than clipped abbreviations; the permanent rule is recorded in the nomenclature guide.

Making aromatic-system or multicenter-bond participants ordered is not an alternative to this
separation. Their participant order is a representational frame, not chemical identity; the
position-sensitive electron-count data must travel through the corresponding frame action.

The final implementation stage must update the permanent development guides and property-test
descriptions to the implemented names, then remove the dated doc-186 TODO markers from
`docs/development/nomenclature.md`, `docs/development/data-types.md`, and
`docs/development/property-tests.md`. Until then, those guides deliberately state the target
semantics so that the implementation can be reviewed against them.

## Crate boundary

Aggregate canonicalization belongs in `umol-graph-ir`, beside the representation whose complete
internal shape it must inspect and remap. The crate responsibilities are:

- `umol-graph-core` supplies automorphism groups, canonical labels, and the explicit algorithm
  selector;
- `umol-graph-ir` owns the aggregate trait, typed incidence encoding, stereo coset and frame
  handling, complete entity remapping, and post-hoc canonical placement of constraints; and
- `umol-graph` may translate higher-level model and operation configuration into the public graph-IR
  context, but does not inspect, reorder, or reconstruct graph-IR internals.

The aggregate context does not require the complete `StereoModel`. Its `kind_models` and element
scopes govern which stereo entities perception attempts to derive, while fluxionality belongs to
the chemical interpretation of that perception. Canonicalization neither perceives nor validates
stereo: every stereo entity already stored in the IR participates regardless of whether a
perception model would have generated it. Of the current model parameters, only `para_stereo`
affects canonicalization by enabling iterative stereo-sensitive symmetry refinement.

The graph-IR context therefore needs only the following inputs at present:

```rust
pub struct CanonicalizationContext {
    pub para_stereo: bool,
    pub automorphism_algorithm: AutomorphismAlgorithm,
}
```

`CanonicalizationContext` is the explicit low-level carrier needed by `umol-graph-ir`; it is not
itself a model or config. At the `umol-graph` layer, the semantic and operational sources remain
separate: `StereoModel::para_stereo` supplies the semantic choice and `CanonicalizationConfig`
contains the `automorphism_algorithm`. The graph layer constructs the graph-IR context from those
two inputs. Do not add a one-field `CanonicalizationModel` that merely duplicates the projection of
`StereoModel`.

Stored stereo entities participate under either `para_stereo` setting. When it is false,
canonicalization performs one stereo-sensitive refinement from the constitution-level partition
and does not feed the result back into another stereo refinement. When it is true, the
stereo-sensitive result is fed back until the partition stabilizes. The context does not carry a
configurable maximum iteration count: an operational cutoff must not determine the canonical form.
The implementation must make the fixpoint refinement monotonic and establish its finite termination
condition.

### Domain and failure semantics

Aggregate canonicalization supports non-ground form values; groundness is not a precondition.
Intrinsic fixed-frame normalization is applied to every carried form value, and an intrinsic
`Contradiction` remains an error of the aggregate operation rather than being ignored or repaired.

Aggregate canonicalization is likewise partial for `Reaction`. After representation integrity
has been established, it reports `Contradiction` when intrinsic normalization of an LHS or delta
value fails, when the deltas cannot be normalized into a coherent before/after transition against
the LHS, or when the normalized transition cannot materialize two referentially intact sides of a
reaction span. The last case is reaction consistency or span materializability, not a DPO
condition. Canonicalization does not test chemistry invariants or conformance and does not test the
host- and match-dependent DPO dangling or identification conditions. A constructed
`ReactionSpan` has already established side integrity and can fail aggregate canonicalization
only through its defensive integrity check or intrinsic normalization of a carried value.

Valid stereo ligand frames may be reordered as part of canonical frame selection only when the
stereo coset is transported through the corresponding frame permutation. Undetermined stereo
configurations are invariant under frame changes, so their ligand frames may be reordered without
introducing a determined assignment.

Malformed stereo data are outside the canonicalization domain. A ligand count inconsistent with the
declared stereo kind, an out-of-range concrete coset or coset-domain member, and an explicit
permutation with the wrong degree are tier-1 representation-integrity failures: they do not denote a
stereo configuration that canonicalization could preserve. Checked molecule construction and every
public mutation or conversion producing a molecule must reject them. Aggregate canonicalization
returns a typed integrity error, rather than a canonical representative or a panic, if a compromised
or unchecked value nevertheless reaches it. This error is distinct from the intrinsic lattice
`Contradiction` returned for a validly represented but unsatisfiable form value.

The stereo checks divide across the validation tiers as follows:

- tier-1 representation integrity covers the ligand-frame length for a declared kind, every
  explicit coset index or variable-domain member, the degree of explicit permutations, ligand
  positions named by topicity constraints, the kinds permitted on stereo atoms and stereo bonds,
  and the referential and incidence rules for ordinary and virtual ligands;
- tier-2 invariants validation covers agreement between atom `#T` or bond `#C` constraints and a
  separately stored stereo entity, including agreement between inline and molecule-level forms;
- tier-3 conformance covers the kinds and sites admitted by `StereoModel`, whether a structurally
  valid entity can be perceived under that model, and the graph-symmetry-derived `#g`, `#o`, and
  `#p` assertions. The stored `#f` assertion is dynamical and is not statically derived; and
- no additional tier-2 stereo property beyond model-independent constraint satisfaction is
  currently justified. Malformed encoding is tier 1, while stereogenicity and realizability require
  the selected stereo model and are tier 3. A concrete label on a structurally valid site that the
  selected model does not regard as stereogenic remains representable and is not a universal
  physical contradiction.

`ImproperOnAchiral` does not belong in another tier and should be removed. Whether a stereo class's
coset algebra can encode handedness does not imply that the represented site cannot be prochiral or
carry an enantiotopic assertion; those are established by the symmetry-derived conformance checks.

Non-ground stereo values remain in the canonicalization domain. A fully
`StereoConfigurationAst::Undetermined` value has no declared kind and therefore no kind-relative
coset or frame-arity check. `Kinded(kind, StereoCoset::Undetermined)` is likewise valid, but its
ligand frame must already have `kind.degree()` positions. Literal sets and symbolic terms are valid
provided every explicit coset index in the value or variable domain is in range and every explicit
permutation has the declared kind's degree. A later semantic validator may return
`Underdetermined`; neither integrity validation nor canonicalization may assume a literal coset or
panic while extracting one.

The same representation-shape rule is already relevant outside stereo. Aromatic systems and
multicenter bonds store electron-count vectors positionally alongside their participant lists. The
reference traversal used by `Molecule::try_from_entries` verifies that every participant atom
exists, but it does not compare the two lengths. `EntityStructureValidator` does perform both length
checks for literal electron-count vectors. These are representation-integrity checks, not chemistry
contradictions: an undetermined vector is valid, while a concrete vector of a different length
cannot assign one contribution to each participant. They must move into the same shared
representation-integrity gate as the stereo shape checks.

The current implementation does not yet establish this contract. `Molecule::try_from_entries`
checks entity and constraint references only. `EntityStructureValidator` separately mixes the two
electron-vector shape checks with semantic entity-structure checks such as graph simplicity and
relation uniqueness, so invoking that validator wholesale is not the construction fix. Molecule
construction, molecule-DSL raise, and editor publication need one graph-IR-owned
representation-integrity implementation containing the reference, parallel-collection, and stereo
shape checks. Model-independent constraint satisfaction remains an explicit tier-2 validation pass,
and source-format-specific stereo interpretation remains in the source-format layer.

The existing entity validator also does not check stereo ligand arity or coset domains. Those two
checks currently live in `StereoConformanceValidator` as `LigandArity` and `CosetOutOfRange`
contradictions, which incorrectly classifies them as tier-3 chemistry conformance. At the leaf level,
`canon_coset` deliberately omits the range check while `StereoKind::act` assumes a valid index and
can panic. Construction and mutation integrity must therefore be brought into conformance before
aggregate canonicalization relies on valid stereo frame actions; canonicalization must call the
same integrity implementation defensively rather than add a second partial implementation.

### Required integrity-check and validator restructuring

Tier 1 and semantic validation must use different operation families. `umol-graph-ir` owns integrity
checks, named `check_integrity`, with `*IntegrityError` failures. These checks return
`Result<(), *IntegrityError>` and never return `Solution`, `Underdetermined`, or `Contradictory`.
There is no `*Checker` object. `umol-graph` owns validators: `*InvariantsValidator` for tier 2 and
`*ConformanceValidator` for tier 3. Both validator families retain
`Result<Solution<_, _>, _>` because a coherent non-ground value may leave a semantic question
underdetermined.

The graph-IR restructuring is:

- add `Molecule::check_integrity`, `Reaction::check_integrity`, and the corresponding operation
  for `ReactionSpan`, with `MoleculeIntegrityError`, `ReactionIntegrityError`, and the analogous
  span error;
- make checked entry construction, molecule and reaction DSL raise, TableIR raise, Python entry
  construction, builder/editor finalization, reaction-side materialization, and aggregate
  canonicalization share those checks. Boundary layers translate the integrity error but do not
  reproduce the checks;
- keep asserted construction only for producers that establish the same contract by construction;
  the checked and asserted routes differ in failure reporting, not accepted values;
- move aromatic and multicenter electron-vector length checks out of
  `EntityStructureValidator` and into molecule integrity, together with the stereo representation
  checks listed above; and
- replace `ReactionIntegrityValidator` with `Reaction::check_integrity`. Its invalid-reference
  and incidence-mismatch cases become integrity errors rather than semantic contradictions.

After this migration, `umol-graph-ir` must not define `*Validator` types. Existing graph-IR
validators split by tier rather than moving as one block:

- the remaining model-independent entity-structure checks move to an invariants validator in
  `umol-graph`;
- remove the span-level DPO-validator entry point because it only confirms span integrity. Review
  the reaction-level entry point separately: a check that predicts whether deltas materialize a
  span is not DPO validation, while actual DPO dangling and identification conditions belong to
  application with a supplied host and match;
- `ConstraintValidator` and its incidence, ring, relational, and molecule-scope components move to
  tier-2 invariants validation in `umol-graph`; and
- `ConnectivityValidator` and `ConnectivityModel` move together to `umol-graph`, with the validator
  renamed `ConnectivityConformanceValidator` because its result depends on the selected model.
  Add `ConnectivityModel` to `ChemistryModel` and run connectivity first in the conformance pass,
  immediately before valence conformance.

The existing graph validators then need the following cleanup:

- remove `validate_integrity`, `EntityStructureValidator`, and `ConstraintValidator` from the
  composite `Validator`. The composite runs `validate_invariants` followed by
  `validate_conformance`; successful graph-IR construction is its integrity precondition.
  Connectivity runs first within `validate_conformance`, immediately before valence conformance;
- remove `Validator::validate_atom`. It is only a partial combination of the valence and spin
  invariants for a free-standing `AtomForm`, yet constructing the model-bearing composite provides
  no input to either operation. Keep the focused `validate_atom` methods on
  `ValenceInvariantsValidator` and `SpinInvariantsValidator` public; a top-level `AtomValidator` may
  be designed separately if a complete atom-validation operation is later required;
- remove `LigandArity`, `CosetOutOfRange`, `ImproperOnAchiral`, `MissingStereoAtom`, and
  `MissingStereoBond` from `StereoConformanceValidator`. The first two become integrity errors, the
  third disappears, and the last two are covered by tier-2 constraint satisfaction;
- leave `StereoConformanceValidator` with model/perception failures and the symmetry-derived `#g`,
  `#o`, and `#p` checks, handling a non-ground configuration as `Underdetermined` rather than
  extracting a literal with an assertion; and
- align contradiction and error names with their validator types. In particular,
  `AromaticityValidatorContradiction` becomes `AromaticityConformanceContradiction`, while
  `StereoValidatorContradiction` and `StereoValidatorError` become
  `StereoConformanceContradiction` and `StereoConformanceError`.

The domain stem of every migrated subordinate invariants validator still requires review; the
settled rule here is the `*InvariantsValidator` or `*ConformanceValidator` suffix and graph-layer
ownership. All validators remain public, including focused validators used independently by
resolvers, transformers, or callers that want to check only one semantic property. The migration
must select and approve the public names rather than preserving misleading names or reducing
visibility to avoid the decision.

The integrity and validator migration proceeds in four green stages:

1. Add the three graph-IR-owned `check_integrity` operations and `*IntegrityError` types. Molecule
   integrity consolidates reference checks, aromatic and multicenter electron-vector lengths,
   stereo ligand arity, explicit coset domains, and explicit permutation degree. Route checked and
   asserted construction, DSL and external-format raise, Python construction, and builder/editor
   publication through the shared checks.
2. Replace `ReactionIntegrityValidator` with `Reaction::check_integrity`, and unify
   `ReactionSpan` checked construction with its integrity operation. Preserve the distinction
   between a permissive lhs-plus-deltas reaction and a span whose two projections are molecules.
3. Move every semantic validator to `umol-graph`: the non-integrity portion of
   `EntityStructureValidator`, the aggregate and focused constraint validators, and connectivity.
   Rebalance stereo validation and apply the approved invariants/conformance names. All resulting
   validators are public.
4. Rewire the composite `Validator`, reaction application, resolvers, transformers, substructure
   operations, Rust and Python callers, tests, property tests, specifications, and rustdoc. Remove
   the old graph-IR validator modules and exports only after every consumer has moved and the
   workspace is green.

The phrase *non-simple input* is too narrow for the current reaction-application precondition.
`Reaction::apply` presently runs the whole `EntityStructureValidator` over both the reaction LHS
and the host. In addition to localized-bond self-loops and parallel bonds, that validator rejects
duplicate or role-conflicting dative participants, parallel dative and noncovalent relations,
duplicate or overlapping aromatic participants, duplicate or identical multicenter participant
sets, and repeated stereo sites. Its electron-vector length cases move to integrity; the remaining
cases are model-independent structural invariants. Before the migration rewires `apply`, decide
which of those invariants are actual application preconditions. The decision must be stated per
property rather than summarized as graph simplicity.

Nauty handles disconnected graphs directly. Splitting a molecule into connected components is
therefore not required for correctness and would complicate constraints or other data that span
components. Component-wise processing may be considered later as a measured optimization, not as
part of the canonicalization semantics.

Canonicalization does not identify or add stereo entities. The context's `para_stereo` setting
determines whether stereo-sensitive graph refinement iterates to a fixpoint. Meso compounds and
other cases in which constitution and stereo symmetry interact must be included in the correctness
corpus.

The construction of a canonical Kekule form and aromaticity or stereo perception are outside this
work. Canonicalization preserves and renumbers the structure represented by the input IR; it does
not perceive, resolve, or replace that structure.

Previous benchmarks found the full incidence graph comparatively slow next to the graph-and-
overlays representation. Those measurements were made for a different operation and did not
separate construction cost from algorithm cost, but they make representation cost an early design
question here. The likely cost is subdivision of every localized bond: for `V` atoms and `E` bonds,
the topological incidence graph has `V + E` nodes and `2E` edges, while molecular overlays are
usually sparse and add far fewer nodes. In a connected molecular graph, `E` is ordinarily at least
of the same order as `V` and is often greater.

This explanation remains a hypothesis until measured. The first canonicalization benchmark must
compare the compact graph-and-overlays representation with materialized full incidence, record
construction and canonical-labeling time separately, and include node and edge counts. The corpus
must cover ordinary molecules as well as deliberately overlay-heavy cases. The comparison must be
made separately for topology, non-stereo constitution, and stereo-sensitive canonicalization,
because their refinement costs differ.

Nauty consumes vertex colors rather than edge colors, so some exact encodings may still require
subdivision at its adapter boundary. That does not establish that every localized bond should be a
node in the common molecular representation. The semantic coverage of the full incidence model is
required; its current fully materialized Levi-graph shape is not presumed to be the final
canonicalization carrier. If the compact representation wins materially, the common incidence
facility should be revised to retain localized bonds as graph edges and lift only the entities or
typed incidences that require it. A second canonicalization-only molecule graph must not be added.

## Canonicalization graph encoding

There is no separate "exact canonicalization graph." The existing incidence-graph facility must be
extended so that the encoding supplied to the canonical-labeling algorithm preserves every
structural distinction selected by the requested canonicalization level. Constraints are excluded
from canonical labeling at every level. Within the selected structural features, isomorphism of the
colored encodings must agree with equivalence under entity renumbering and the supplied
`para_stereo` semantics where stereo is included.

The current graph representations are:

- the raw molecular graph, with atoms as nodes and localized bonds as edges;
- `IncidenceGraph` with `IncidenceNodeSelection::topological()`, which represents atoms and
  localized bonds as nodes;
- `IncidenceGraph` with `IncidenceNodeSelection::constitution()`, which additionally represents
  dative bonds, aromatic systems, multicenter bonds, and noncovalent bonds;
- `IncidenceGraph` with `IncidenceNodeSelection::full()`, which additionally represents stereo
  atoms and stereo bonds and is the representation currently used by `graph_symmetry`; and
- `SubdividedGraph` in `umol-graph-core`, the lower-level graph obtained by subdividing each edge
  once. It is not a molecule incidence graph; it is a mechanism for presenting edge distinctions to
  a vertex-colored algorithm.

The semantic coverage of `IncidenceGraph::full()` is the starting point for canonicalization
because it includes all eight entity kinds, but neither its current information nor its fully
materialized representation is assumed to be sufficient for exact, efficient canonicalization:

- dative donor and acceptor roles are not distinguished structurally;
- stereo nodes attach only to their sites, while their ligand frames and cosets are interpreted by
  separate symmetry code;
- per-participant aromatic and multicenter electron counts must remain associated with the correct
  participant positions;
- `ConstitutionColoring` produces `u64` hashes, whereas canonicalization must construct
  collision-free equality classes from the represented values.

These are deficiencies of the common incidence representation, not reasons to introduce a parallel
canonicalization-only system. They also matter to other consumers that require exact molecular
automorphisms. Incidence-based substructure matching currently tolerates an over-approximation
because it post-verifies overlay roles and values; canonical labeling has no corresponding filter.

`IncidenceGraph` must therefore be extended to retain typed incidence information. Dative
donor/acceptor roles, participant-associated aromatic and multicenter electron counts, and stereo
site and ligand-frame roles belong to the incidences between entity nodes. The underlying `Graph`
has uncolored edges, so the nauty adapter may subdivide those typed incidences into colored vertices;
that is an algorithm-input transformation of the same `IncidenceGraph`, not another molecular graph
model. Node values must likewise be ranked into collision-free equality classes rather than reduced
to `u64` hashes.

Entity-level and molecule-level constraints do not participate in canonical labeling. They are not
part of constitution and do not establish entity or molecular structural identity. They therefore
neither distinguish incidence nodes nor break structural automorphism orbits. After structural
canonical labeling establishes the canonical image and its remaining automorphisms, normalized
entity-level and molecule-level constraints select a canonical placement among those structurally
equivalent frames. The resulting remapping transports every constraint reference together with the
entities it names.

Derived `==` is exact equality of the stored IR and therefore includes constraints, ids, ordering,
and non-normal value encodings. `canonical_eq` compares complete canonical IR forms and also
includes constraints. This is required for patterns, where constraints are not redundant with the
structural description. The constraints' post-hoc participation in selecting a complete canonical
form does not make them structural colors or allow them to alter the underlying structural orbits.

## Canonical comparison order

For a fixed aggregate-canonicalization context, consider every complete `Molecule` obtained by a
valid dense remapping of all eight entity families and every corresponding participant-frame action.
Before comparison, every carried entity and constraint value undergoes intrinsic, fixed-frame
normalization. Unordered aromatic-system and multicenter-bond participant reorders
permute their participant-indexed electron counts. Stereo ligand frames may likewise be reordered,
provided that the stereo coset is transported through the same frame permutation. The canonical
form is the minimum resulting complete IR under a specified typed total order.

The typed order, rather than a rendered DSL string, is normative. A canonical DSL rendering is a
derived representation of the canonical IR; parsing, rendering, whitespace, shorthand, or other
surface-syntax changes must not alter the selected entity numbering. The implementation may use a
collision-free comparison key or another optimized representation, but it must agree with the typed
order and must not replace it with a hash.

The comparison key is a private typed data object with explicit ordering semantics, not a rendered
encoding and not a tuple-of-vectors public contract. A straightforward implementation may contain a
vector block for each of the eight entity families followed by constraints, with typed row keys for
participants, frames, and normalized form values. Its storage shape is private; its comparison order
is the contract. The choice among valid total orders is mathematically arbitrary, like choosing a
fixed seed, but becomes a public compatibility contract once published. For the same input and
semantic context, canonical numbering must remain stable across library releases, platforms, and
supported algorithm implementations.

The comparison schema therefore assigns explicit, platform-independent positions to every entity
block, row component, and variant tag. Its `Ord` implementation is manual or delegates only to
types whose order is explicitly frozen. It must not depend on hash iteration, compiler layout,
implicit Rust enum discriminants, locale, or the output of a particular canonical-labeling backend.
An implementation may optimize the search freely but must reproduce the scheme's exact result.
Correct algorithms selected through the operational algorithm enum cannot choose different
canonical representatives.

Protocol schemas are useful governance precedents, not encodings to reuse. [Protocol
Buffers](https://protobuf.dev/programming-guides/proto3/#assigning-field-numbers) gives fields
permanent numeric identities and forbids reusing reserved numbers; [Cap'n
Proto](https://capnproto.org/language.html#evolving-your-protocol) requires consecutive ordinals in
addition order while permitting declarations to be rearranged without changing those ordinals;
[Thrift](https://thrift.apache.org/docs/idl#field-id) provides explicit field ids but does not
require them in the grammar. The applicable lesson is explicit stable positions plus append-only
extension. Protocol Buffers explicitly [does not promise canonical
serialization](https://protobuf.dev/programming-guides/serialization-not-canonical/), so protocol
wire bytes are not a comparison key and the canonical molecule schema is not implemented by
serializing through one of these formats.

The comparison scheme is fixed by the library contract rather than selected through the operation
context. If an incompatible replacement ever becomes unavoidable, it requires a separately
versioned API and cannot silently change the published scheme or its existing results.

Freezing the schema does not freeze the entity model. The compatibility unit is the canonical
representative of a value expressible in the schema version that first covered it. The existing
entity-kind blocks, field components, constraint variants, and their comparison order receive
explicit stable schema positions. A future entity kind or constraint variant is added at an
append-only extension position; it must not renumber, reorder, or reinterpret an existing position.
The implementation must not derive these positions from a Rust enum's declaration order.

Consequently, extending the schema must leave the canonical numbering and complete canonical form
of every molecule containing none of the new entity kinds or constraint variants unchanged. For a
molecule that does use an extension, the extended schema defines a coherent total order, but its
structure and canonical ordering carry no compatibility guarantee with library versions predating
that extension because those versions could not represent the molecule. From the version that
introduces the extension onward, the new schema positions are themselves frozen: every later
append-only extension must preserve the molecule's canonical numbering and complete canonical form
whenever the molecule uses none of the still-later additions. Equivalently, if one schema is an
append-only extension of another, every molecule expressible in the earlier schema has the same
canonical representative under both. This cumulative promise permits the entity-model growth
described in doc 117 while keeping already representable molecules stable. It applies to the public
canonical representative; an internal comparison-key encoding may change provided it reproduces
that representative.

The architectural requirement is that `Molecule::incidence_graph` remains the single
molecule-to-incidence facility used by symmetry, substructure, and canonicalization. The facility
may need to expose a compact graph-and-overlays form and explicit adapters, rather than requiring
all consumers to start from a fully subdivided graph.

Stereo-sensitive colors are refined according to the context's `para_stereo` setting; when enabled,
refinement proceeds to a stable partition before the canonical labeling is accepted. The existing
stereo symmetry operations provide the frame and coset actions needed for this step. They are reused
for canonicalization rather than invoking stereo perception.

## External precedents and points of comparison

The source trees under `materials/codes` are references for particular parts of the design, not
drop-in specifications. Their chemical models and canonical outputs differ from `Molecule`, so
literal canonical numbers need not agree. The useful comparisons are invariance under renumbering,
symmetry partitions, treatment of selected features, stereo behavior, and representation cost.

### RDKit

RDKit's current canonical-ranking implementation is in
`materials/codes/rdkit/Code/GraphMol/new_canon.h` and `new_canon.cpp`. It is useful for two API and
algorithmic distinctions:

- ranking without tie breaking yields symmetry classes, while tie breaking produces a canonical
  order; and
- chirality, isotopes, atom maps, and related features can participate in ranking independently.

The first distinction parallels the separation between an orbit/equivalence partition used during
refinement and the final canonical frame. The second is a precedent for the parameterized ordering
operation, while `Molecule::canonical_eq` remains the unparameterized complete comparison.
RDKit's stereo regression corpus is also useful for meso and symmetry-dependent stereo cases. Its
canonical numbers are implementation-specific and are not an oracle for umol's numbering.

### CDK

CDK's individualization-refinement implementation under
`materials/codes/cdk/tool/group` is a comparatively readable reference for partition refinement,
automorphism handling, and optional inclusion of elements and bond orders. In particular,
`AbstractDiscretePartitionRefiner.java` and `AtomRefinable.java` are useful implementation
comparisons. CDK's separate `base/standard/.../Canon.java` explicitly documents known examples for
which its symmetry classes are incorrect, so it must not be treated as an exact oracle. It remains
useful for test cases and for comparing feature-selection APIs.

### InChI

InChI's `CanonGraph` in
`materials/codes/InChI/INCHI-1-SRC/INCHI_BASE/src/ichican2.c` is an adapted implementation of the
McKay 1981 algorithm. It exposes both symmetry ranks and canonical ranks and supports progressive
additional coloring: when a later layer can only reduce the earlier automorphism group, the earlier
canonical numbering and orbits constrain the subsequent search. This is a useful precedent for
reusing topology and constitution refinement when adding isotope, tautomer, or stereo layers.

The graph passed to `CanonGraph` does **not** use general edge colors. `Graph` is an adjacency-list
type, and the ordinary calls to `CreateNeighList` pass `0` for its optional double-bond-neighbor
duplication mode. Element, connection count, hydrogen, isotope, and tautomer information instead
enters through initial atom invariants and additional canonicalization layers. InChI therefore
demonstrates that a McKay-style chemical canonicalizer need not subdivide every bond. It does not,
however, show that umol can simply omit bond values: InChI canonicalizes its own normalized,
layered representation, whereas `Molecule` must preserve explicit localized bonds and all other
selected structural distinctions.

This precedent reinforces the need to benchmark the canonicalization carrier before committing to
full edge subdivision. It also suggests measuring whether progressive refinement over a compact
base graph can avoid materializing bond nodes, without importing InChI's chemical normalization
semantics.

Doc 188 records a full review of the InChI 1.07 implementation. Four of its findings bear on this
plan directly; none changes the plan's structure or ordering.

- The internal search strategy of S6–S8 has a proven implementation shape in InChI's chained
  canonicalization: up to eight `CanonGraph` passes, each constrained to reproduce the previous
  layer's certificate (`zb_rho_fix`) while minimizing only the new layer, with per-layer
  certificate memoization, a guard that skips a pass when the coarser orbits already separate the
  new colors, and an auxiliary pre-ranking that steers branch order. When the S4b comparison
  schema orders entity blocks coarse to fine, this chaining is a correct search strategy for the
  typed-order minimum. It constrains rather than commits, so it is compatible with the rule that
  intermediate numberings are never committed as canonical frames; all of it is internal
  optimization under S9d's fixture-identity rule.
- InChI discards the automorphism generators after `CanonGraph` and consequently re-derives the
  group in a second backtracking search when minimizing stereo descriptors (`map_stereo_bonds4` in
  `ichimap4.c`). The nauty adapter already returns generators and orbits; S8a's coset actions and
  S9b's constraint placement should consume them directly, avoiding that second search by design.
- InChI's stereo canonicalization loop is not monotonic: discovering a non-stereocenter mid-search
  forces a total restart, because InChI perceives and removes stereo entities during
  canonicalization. This plan preserves every stereo entity and performs no perception, so the
  restart pathology cannot arise; S9a's monotonic-refinement and termination argument is easier
  than this closest precedent.
- The no-edge-color observation above has a sharper form: bond orders are not part of InChI
  identity at all, which is the entire reason its carrier needs no bond information. The S5d
  compact-carrier alternative for umol therefore remains typed incidences presented as adapter
  colors, never omission of bond values.

Doc 188 also describes a dev-only vendored build of the C source with internal symbols exposed
(its verification step 0). Its `CanonGraph` entry is an independent implementation of
vertex-colored-graph canonical labeling and orbit partitions, usable offline to generate
checked-in orbit fixtures (S0b) and as an algorithm-layer cross-check while nauty is the sole
automorphism backend (S5d, S6b), keeping external programs out of the test dependencies. It is an
oracle for the graph layer only, never for umol's molecule numbering: InChI's graph carries no
bond values, so agreement holds only where the colored encodings coincide. This build is not a
prerequisite for any stage of this plan: S0b is satisfiable with hand-derived precedent cases and
bounded internal cases, and its fixtures remain extensible later. The earliest point where the
oracle adds strength is S5d/S6b, so whether to build it can be decided at S5d without affecting
the critical path.

### Internal verification

For ordinary molecular graphs, RDKit, CDK, and InChI can supply comparison cases and expected
symmetry behavior. None covers the full DAMNSS entity model. Exact verification of dative,
aromatic, multicenter, noncovalent, and stereo overlay handling therefore requires bounded
exhaustive internal cases in addition to renumbering-invariance and idempotence properties.

## Three meanings of progressive canonicalization

Three related but distinct progressions must not be conflated.

### Internal canonicalization strategy

Canonical labeling may refine structural equivalence in stages: topology, non-stereo constitution,
then stereo. The useful intermediate result is an orbit partition or a progressively refined
labeling problem. An intermediate numbering must not be committed as the canonical molecule frame,
because a higher level may distinguish entities that were symmetric at a lower level. The final
canonical frame is applied to the original molecule only after the requested refinement has
finished. Stereo refinement uses the constitution-level classes and, when
the context's `para_stereo` setting is enabled, iterates to a fixpoint.

### Canonicalization excluding higher-level features

A caller may intentionally request a coarser canonicalization relation. Topology-only
canonicalization ignores non-topological structure; non-stereo canonicalization ignores stereo and
therefore treats stereo isomers as equal for the purpose of choosing an order. This is analogous to
supplying a comparator to `sort_by`: values that compare equal are retained, not collapsed or
removed.

The closed selector is:

```rust
pub enum CanonicalizationLevel {
    Topology,
    Constitution,
    Full,
}
```

`Topology` contains atoms and localized bonds. `Constitution` adds dative bonds, aromatic systems,
multicenter bonds, and noncovalent bonds. `Full` adds stereo atoms and stereo bonds. Constraints are
excluded from the structural levels. Para-stereo is not a fourth level: both one-pass stereo
refinement and the para-stereo fixpoint operate at `Full`, with the behavior selected by
`CanonicalizationContext::para_stereo`.

The parameterized operation returns the complete original `Molecule` in the frame selected under
the requested relation. Every entity, field, stereo assignment, and constraint is retained and
transported through the resulting remapping. The operation performs no projection, stripping,
resolution, or normalization. When excluded features distinguish entities tied by the selected
relation, the returned complete molecule need not be a unique canonical representative with respect
to those excluded features; that is inherent in requesting the coarser ordering.

More precisely, the selected structural layer of the result is in canonical form and the complete
result is a remapping of the original molecule. The ordering of excluded features within an
automorphism class of that selected layer is not determined. Two differently numbered inputs may
therefore produce complete outputs that differ by an automorphism of the selected layer even though
their selected layers have the same canonical form. Breaking such a tie with an excluded feature
would make that feature part of the ordering and is deliberately not done.

Constraints are excluded from all of these structural levels.

The public operations are `canonicalize_by(level, context)` and
`canonical_eq_by(other, level, context)`. The first returns the complete input transported into the
frame selected at `level`; the second compares only the selected structural layers. It must not be
implemented by comparing the complete outputs of `canonicalize_by`, because excluded features may
remain differently ordered within an automorphism class of the selected level.

These operations are distinct from `Molecule::canonical_eq`, whose name and semantics are not
parameterized by a structural level. `canonical_eq` compares the complete set of available entity
kinds, including stereo entities. It must not adopt the reduced equality relation supplied to the
ordering operation or use the coarser result as its canonical representative.

Feature and level selectors follow one repository rule. A `*Features` type is a bitflag set of
independently combinable switches; a `*Level` type is a closed enum of nested named layers.
Accordingly, rename `ConstitutionFeatures` to `MoleculeColoringFeatures`, because it also contains
`STEREO_KIND`, and retain its bitflag semantics. Replace `IncidenceNodeSelection` with the
`IncidenceLevel` enum using `Topology`, `Constitution`, and `Full`. Canonicalization uses the parallel
`CanonicalizationLevel`. The level enums remove the misleading `OVERLAYS` flag, which currently
excludes the two stereo overlay kinds.

### Progressive implementation

Implementation proceeds in four stages: complete and verify topology-level canonicalization; extend
it to all non-stereo structural entities and typed incidences; add one-pass stereo refinement with
`para_stereo == false`; then add the full para-stereo fixpoint. The latter two stages both implement
the `Full` structural level. Each implementation stage must provide a complete, tested contract for
its supported level and reuse the same incidence and remapping facilities; it is not a temporary
parallel implementation that is discarded at the next stage. The unqualified complete
`canonicalize` and `canonical_eq` surface lands only with the final stage.

Before reaction canonicalization is implemented, the core stereo delta vocabulary must become
absolute-only. Remove `Apply`, `Swap`, and `Mirror` from `StereoAtomDelta` and `StereoBondDelta`.
Every retained modification then states an explicit old and new value, like the modifications for
the other entity kinds. The reaction DSL and Python bindings must expose the same core vocabulary;
they must not preserve the removed variants as a second route into host-relative delta semantics.

This is a semantic simplification rather than a canonicalization workaround. Absolute deltas form a
uniform declarative core: they are closed under inversion, denote explicit before/after states, and
have a faithful reaction-span image. Relative operations such as "invert this stereocenter" remain
reasonable higher-level operations, but they produce one absolute reaction when the input
configuration determines the result or a finite collection of absolute reactions when the rule is
configuration-generic. The collection and its surface API are separate design work; they are not
additional core delta variants. Removing the variants does not remove the stereo permutation and
frame-action machinery used by remapping and canonicalization.

This changes the implementation order and removes the need for a parallel reaction-canonicalization
path:

- migrate the Rust delta enums, reaction DSL, Python bindings, specification, examples, and
  generated inputs to the absolute-only vocabulary before implementing reaction canonicalization;
- simplify delta application, inversion, fixed-frame normalization, remapping, and composition
  by removing their host-relative stereo branches;
- establish that reaction/span conversion preserves the match domain and application result of the
  normalized reaction, in addition to reproducing the same materialized span; and
- canonicalize `Reaction` through `ReactionSpan` after the molecule and span operations are
  complete. A direct action-level canonicalizer for relative delta variants is no longer needed.

`Reaction` currently implements only the partial fixed-frame normalization described above, and
`ReactionSpan` has no canonicalization operation. Numbering-invariant reaction canonicalization
cannot in general be implemented by
canonicalizing the bare LHS and then remapping the deltas, because a transformation can distinguish
entities that are symmetric in the LHS and added entities do not occur in the LHS at all. The
reaction change structure must therefore participate in canonical labeling.

`ReactionSpan` is the primary aggregate to canonicalize. It has the same eight entity families
and typed incidence structure as `Molecule`, expressed in the union frame of the two sides, while
each entity value is lifted into an `EntitySpan<T>`. Its collision-free structural comparison key
therefore extends the corresponding molecule entity key as follows:

- `Unchanged(value)` is the explicit `Unchanged` schema tag followed by the value key;
- `Added(value)` is the explicit `Added` schema tag followed by the value key;
- `Removed(value)` is the explicit `Removed` schema tag followed by the value key; and
- `Modified { lhs, rhs }` is the explicit `Modified` schema tag followed by the ordered lhs and rhs
  value keys.

The tags have fixed comparison-schema positions and must not be derived from Rust enum
discriminants or replaced by hashes. `Added` and `Removed`, and the two positions within `Modified`,
remain distinct so canonicalization never identifies a reaction with its reverse. A modified entity
uses one incidence frame on both sides; a change of participants is represented structurally as a
removal and an addition. Participant permutations act on both values in `Modified`, including both
sets of participant-indexed electron counts and both stereo configurations. Constraint spans remain
outside the structural coloring, like molecule constraints, and participate post-hoc in selecting
the complete canonical span among structurally equivalent frames.

Molecule constraints have set-like conjunction semantics in this process, not multiset semantics.
Normalize each constraint value, discard duplicate normal-form values, then classify each value by
side membership: present on both sides is `Unchanged`, left only is `Removed`, and right only is
`Added`. Constraint order and duplicate occurrences in a raw `Constraints` store do not contribute
to the reaction or span canonical form. `ConstraintDelta` normalization must use the same set
difference: repeated equal additions or removals collapse, and an equal addition/removal pair
cancels. Materialization against the LHS then enforces continuity: adding a canonical constraint
already present on the LHS or removing one absent from it is a `Contradiction`, not an idempotent
no-op. Occurrence counts must not survive merely to reproduce redundant raw constraint storage.

The numbering-invariant canonical form of `Reaction` is obtained through the span:

1. `Reaction::to_reaction_span` normalizes and materializes the before/after state, retaining its
   existing `Contradiction` result for intrinsically contradictory values, incoherent delta
   transitions, or a product that cannot form a reaction span;
2. aggregate canonicalization selects the canonical union frame of the `ReactionSpan`; and
3. `ReactionSpan::to_reaction` infallibly derives the LHS and deltas in delta normal form.

This is an explicit partial canonicalization contract. Failure to materialize the span is not
silently repaired and is not reported as DPO invalidity. A valid span may still fail a
host-dependent DPO application condition later. Non-ground values and tier-2 or tier-3 chemistry
invalidity remain in the canonicalization domain.

Because the core delta vocabulary is absolute-only, this round trip preserves the complete
declarative reaction semantics rather than merely its result on the stored LHS. It may still replace
an input sequence by its delta normal form, but every normalized modification states the same
before/after transition. The roundtrip properties must therefore establish both that materializing
the recovered reaction reproduces the same span and that the original and recovered reactions have
the same match domain and produce the same result for every accepted match.

#### Reaction/span conversion and property contract

Removing relative stereo deltas does not change the conversion signatures:

```rust
Reaction::to_reaction_span(&self) -> Result<ReactionSpan, Contradiction>
ReactionSpan::to_reaction(&self) -> Reaction
```

`to_reaction_span` remains fallible because an absolute delta collection may still be intrinsically
contradictory, inconsistent with its LHS, or unable to materialize a referentially intact right
side. `to_reaction` remains infallible because construction of `ReactionSpan` establishes both
side projections as part of representation integrity. Removing host-relative operations changes
neither boundary.

For a reaction `r` for which span materialization succeeds, define its delta normalization as

```text
N(r) = r.to_reaction_span()?.to_reaction()
```

The reaction/span APIs and aggregate canonicalization assert the following properties:

- `N(r).lhs == r.lhs`;
- `N(r).to_reaction_span() == r.to_reaction_span()`;
- `N(N(r)) == N(r)` by exact structural equality;
- for every host `h` and explicit molecule correspondence `c` satisfying the correspondence
  preconditions, `r.apply_at(h, c)` and `N(r).apply_at(h, c)` have the same failure or produce
  exactly equal `ReactionDerivation` values.

For a constructed span `s`, define its reaction-derived normal form as

```text
S(s) = s.to_reaction().to_reaction_span().expect("constructed span materializes")
```

Then `S(S(s)) == S(s)`. If `s` already uses canonical constraint values in the set-difference normal
form, `S(s) == s`; otherwise the two spans have the same structural entities and canonically equal
side constraints, but need not be structurally equal. This is deliberate normalization, not loss of
reaction semantics.

The application property uses `apply_at` so that it tests reaction semantics independently of
substructure enumeration and its algorithm selectors. It must be exercised on generated hosts and
explicit correspondences, not only by applying a reaction to the identity occurrence of its own
LHS. Materializable reactions and constructed spans are generated by construction for these success
properties. Separately generated non-materializable reactions introduce one named defect and assert
the exact `to_reaction_span` failure; success properties must not silently narrow an arbitrary
reaction strategy with `if let Ok(...)`.

Let `C_S` be aggregate span canonicalization and define reaction canonicalization by

```text
C_R(r) = C_S(r.to_reaction_span()?).to_reaction()
```

For values on which the displayed partial operations succeed, the canonicalization properties
include:

- `C_S(C_S(s)) == C_S(s)`;
- `C_S(remap(s, p)) == C_S(s)` for every valid dense renumbering `p` of the complete span;
- the left and right projections of `C_S(s)` are respectively equivalent to the projections of `s`
  under the induced side remappings;
- `C_R(C_R(r)) == C_R(r)`;
- `C_R(N(r)) == C_R(r)`;
- `C_R(reverse(C_R(r))) == C_R(reverse(r))`; and
- `C_R` is invariant under every valid dense renumbering of the complete reaction, while preserving
  the represented transformation.

The last property compares exact canonical representatives after transporting the complete LHS,
deltas, participant frames, position-sensitive relation data, stereo configurations, and constraint
references. It does not canonicalize the LHS independently of the reaction change structure.

Canonicalizing only the LHS is a possible optimization, not the baseline semantics. With no added
entities, it is sufficient when every changed or removed LHS entity lies in a singleton LHS
automorphism orbit; changes applied invariantly across a complete orbit are also harmless. A change
that distinguishes one member of an LHS orbit must participate in canonical labeling. Added
entities require additional work even when the LHS is rigid, because the added substructure may
have its own automorphisms. Constraint-only changes may likewise select among structurally
equivalent frames during the post-hoc step.

The general shortcut condition is invariance of the complete transformation under the remaining
LHS automorphism group, including a compatible action on added entities. Establishing that condition
can approach the cost of canonicalizing the span. The initial implementation therefore uses the
span route uniformly, apart from a possible trivial empty-deltas path. A later measured optimization
may reuse the LHS refinement or orbit partition as the initial partition for span canonicalization;
it must not commit to an LHS numbering before reaction changes have refined the remaining symmetry.

### Remapping nomenclature prerequisite

Canonicalization must use one remapping operation vocabulary across the two remapping carriers:

- `umol_graph_core::Remapping` maps graph `NodeId` and `EdgeId` values; and
- `IdRemapping` maps the eight molecule entity-id families.

The carrier exposes `map_*` methods for looking up the image of one id. A represented value or
container exposes `remap` for transporting itself through the carrier. The receiver determines the
work required by that transport: relation-set remapping additionally restores canonical participant
order and permutes position-sensitive relation data, but this is the relation set's implementation
of remapping rather than a distinct operation.

Before canonicalization is implemented, graph-core relation sets use `remap` and `try_remap` across
all five forms. The ordinary `remap` route asserts that the remapping covers the receiver;
`try_remap` checks coverage for an independently supplied receiver/remapping pair and returns `None`
on a mismatch. The analogous removal-driven operation is `compact`: it drops relations containing
removed participants and relabels every survivor into the compacted id space. These names must not
change transport, participant-order, payload-permutation, dropping, or failure semantics.

Graph-IR values transported through `IdRemapping` retain the existing `remap` spelling. The
end-to-end `Molecule` operation follows the same pair: `remap` is the asserted route for a
producer-known dense correspondence, while `try_remap` checks an independently supplied
correspondence. The argument type identifies whether graph ids, typed molecule ids, or the complete
molecule namespace are being transported.

## Required molecule-remapping operation

Canonicalization requires a general end-to-end `Molecule` remapping operation. This is a public
graph-IR transformation rather than a canonicalization-specific helper: canonicalization derives a
canonical labeling, represents it as a `MoleculeCorrespondence`, and applies the same operation that
other callers can use to transport a standalone molecule between dense id spaces. Its asserted and
checked routes are `remap` and `try_remap`, respectively.

The operation has the following contract:

- the correspondence source counts equal the molecule's counts for all eight entity families;
- every component correspondence is total on both sides and therefore defines a bijection onto a
  dense target id space;
- topology, relation participants, position-sensitive relation data, stereo frames, entity forms,
  and all references in constraints are transported together;
- it performs no chemistry validation, resolution, attribute canonicalization, repair, compaction,
  or entity removal; and
- failure of the source molecule's representation integrity or of the correspondence to describe
  such a dense renumbering is reported with `Option`.

The implementation coordinates `umol_graph_core::Remapping`, which owns topology and relation
participant transport, with `IdRemapping`, which owns typed references across all eight entity
families. It must not reimplement relation participant sorting or payload permutation at the
`Molecule` layer. The graph-core relation-remapping correction and its immediate reaction-span
consumer remain in doc 185, S3c; this work consumes that corrected facility.

Required properties are:

- semantic preservation, stated directly as
  `source.equiv_under(&remapped, &correspondence)`;
- exact identity remapping;
- inverse round-tripping;
- agreement between sequential remapping and correspondence composition; and
- preservation of referential integrity for every entity family and constraint reference.

The generated cases must exercise crossing permutations, all eight entity families,
position-sensitive relation data and stereo frames, and constraints that reference remapped
entities. These are properties of molecule remapping itself, independent of its use by
canonicalization. Canonicalization properties additionally verify that applying the remapping
derived from the canonical labeling produces the canonical representative.

This operation does not replace embedding into an ambient union namespace. A reaction-side mapping
may target a sparse or larger union id space and can transport entries into that space, but it cannot
produce a standalone dense `Molecule` without a separate dense reindexing.

## Staged implementation plan

Doc 176's aggregate-type and crate renames are complete, so the public canonicalization surface
below uses `Molecule`, `ReactionSpan`, and `Reaction` in `umol-graph-ir`. Doc 185 is complete; its
corrected reaction-span and relation-remapping semantics are prerequisites rather than work
repeated here. The broader trait renames excluded by doc 176 remain outside this plan.

Every subitem ends green unless it is explicitly marked breaking; a breaking subitem includes its
complete caller migration and ends green. New and changed tests follow the test-writing conventions.
Public operation rustdoc states representation and contextual preconditions, exact failure behavior,
and the semantic properties validated by the corresponding property tests.

### S0 — Baselines and selector vocabulary

- **S0a — Canonicalization benchmark baseline.** Add a dedicated graph-IR canonicalization benchmark
  target with a fixed corpus covering ordinary, disconnected, overlay-heavy, stereo, meso, and
  para-stereo structures. Measure raw atom/bond topology canonical labeling, then separately measure
  incidence construction and canonical labeling at topology, constitution, and full levels. Record
  raw and incidence node and edge counts in the benchmark ids. The raw and incidence paths are not
  semantic parity beyond topology: the exact compact overlay-aware path does not yet exist, so its
  comparison with incidence remains S5d. This is additive and establishes the measurable baseline
  before a carrier is changed. [dep: doc 176] **Done.**
- **S0b — Canonicalization correctness corpus.** Add checked-in exact cases derived from the RDKit,
  CDK, and InChI precedents plus bounded internal DAMNSS cases. Record expected orbit partitions and
  renumbering invariance rather than requiring another library's canonical numbers. Keep external
  programs out of the test dependencies. This is additive. [dep: S0a] **Done.**
- **S0c — Coloring-feature terminology.** Rename `ConstitutionFeatures` to
  `MoleculeColoringFeatures` throughout `umol-graph-ir` and its callers, preserving the independent
  bitflag semantics and the existing coloring behavior. Update imports, rustdoc, unit tests, and
  benchmarks in the same breaking red-to-green subitem. [dep: S0a] **Done.**
- **S0d — Structural-level terminology.** Replace the bitflag-shaped `IncidenceNodeSelection` with
  `IncidenceLevel::{Topology, Constitution, Full}` and migrate every incidence, symmetry,
  substructure, test, and benchmark caller. Verify the exact entity-kind membership of each level;
  do not retain a second bitflag route. This is breaking red-to-green. [dep: S0c] **Done.**
- **S0e — Conversion-context terminology.** Rename the `Ctx` associated type to `Context` on
  `FromIr`, `IntoIr`, `TryFromIr`, and `TryIntoIr`, migrate every implementation and caller, and
  update public rustdoc. This is a mechanical breaking red-to-green migration and does not rename the
  traits or crate. [dep: doc 176] **Done.**

### S1 — Fixed-frame normalization vocabulary

- **S1a — `Normalize`, `Equiv`, and `Normalized<T>`.** Replace the current context-free
  `Canonicalize` trait with `Normalize`; rename `canonicalize` to `normalize`, `canonical` or the
  borrowed normal-form projection to `normalized`, `Canonical<T>` to `Normalized<T>`, and the current
  fixed-frame `canonical_eq` operation to `Equiv::equiv`. Keep `Lattice: Normalize`. Migrate every
  leaf, entity, constraint, update, delta, and container implementation together with exact and
  property tests for idempotence, contradictions, and equality on normal forms. This is breaking
  red-to-green. [dep: S0e] **Done.**
- **S1b — Frame-aware equivalence and callers.** Rebase `RelationEquiv::equiv_under`,
  `BiRelationEquiv::equiv_under`, relation-data permutation actions, molecule comparisons, edit
  continuity checks, and transaction checks on `Equiv` without introducing another comparison
  protocol. Migrate every Rust caller and preserve the current participant-order actions for
  aromatic and multicenter electron counts. This is breaking red-to-green. [dep: S1a] **Done.**
- **S1c — Boundary and documentation migration.** Update the existing Python lattice operations,
  DSL and EDN tests, specifications, examples, fuzz targets, and current rustdoc to the normalization
  vocabulary. Remove the old fixed-frame `Canonicalize` implementation on `Reaction`; `Deltas`
  implements `Normalize`. Do not introduce aggregate canonicalization until S4. This is breaking
  red-to-green. [dep: S1b] **Done.**

### S2 — Representation integrity and public validators

- **S2a — Molecule integrity contract.** Add `MoleculeIntegrityError` and the authoritative
  `Molecule::check_integrity`. Consolidate reference resolution, aromatic and multicenter
  participant/electron-vector shape, stereo ligand arity, explicit coset domains, topicity ligand
  positions, permitted stereo entity kinds, and explicit permutation degree. Route checked and
  asserted entry construction through the same implementation and add exact success and error cases
  plus properties that generated checked molecules pass their integrity check. This is additive
  until the constructors are rewired, then breaking red-to-green within the subitem. [dep: S1c]
  **Done.**
- **S2b — Reaction and span integrity contracts.** Add `ReactionIntegrityError` and
  `ReactionSpanIntegrityError` with public `check_integrity` operations. Replace
  `ReactionIntegrityValidator`; make checked construction, DSL raise, side materialization, and span
  projection use the shared checks. Preserve permissive lhs-plus-deltas construction for reactions
  and the stronger two-sided construction invariant for spans. Verify malformed stored references
  and side projections as typed integrity errors rather than contradictions or panics. This is
  breaking red-to-green. [dep: S2a] **Done.**
- **S2c — Public invariants validators.** Move the non-integrity portion of
  `EntityStructureValidator` and the aggregate, incidence, ring, relational, and molecule-scope
  constraint validators from `umol-graph-ir` to `umol-graph`. Apply the approved public
  `*InvariantsValidator` names and keep every aggregate and focused validator public. Preserve
  `Result<Solution<_, _>, _>` and non-ground `Underdetermined` behavior. This is breaking
  red-to-green; the exact public domain stems must be approved before this subitem starts. [dep: S2a]
  **Done.**
- **S2d — Public conformance validators.** Move `ConnectivityValidator` and `ConnectivityModel`
  together to `umol-graph`, rename the validator `ConnectivityConformanceValidator`, add
  `ConnectivityModel` to `ChemistryModel`, and run connectivity first in the conformance pass,
  immediately before valence conformance. Rebalance stereo conformance. Remove `LigandArity`,
  `CosetOutOfRange`, `ImproperOnAchiral`,
  `MissingStereoAtom`, and `MissingStereoBond` from stereo conformance; retain perception/model and
  symmetry-derived checks, with non-ground configuration producing `Underdetermined`. Align the
  aromaticity and stereo contradiction/error names. All resulting validators are public. This is
  breaking red-to-green. [dep: S2b, S2c] **Done.**
- **S2e — Composite validation and operation callers.** Remove `validate_integrity` and the moved
  graph-IR validators from the composite `Validator`; run invariants and then conformance, with
  connectivity first in the conformance pass. Remove the partial composite
  `Validator::validate_atom` operation; keep the focused valence- and spin-invariants atom
  operations public without adding a top-level `AtomValidator`. Remove the duplicate public
  `ValenceInvariants::check` and `check_atom` surface, consolidate that semantic operation behind
  `ValenceInvariantsValidator::validate` and `validate_atom`, and do not retain `check` aliases.
  Rewire resolvers, transformers, substructure operations, Rust and Python callers, specifications,
  and rustdoc. This is breaking red-to-green. [dep: S2c, S2d] **Done.**
- **S2f — Reaction-application preconditions.** Replace `Reaction::apply`'s wholesale use of
  `EntityStructureValidator` with integrity checks and the explicitly approved operation-specific
  preconditions. Before match enumeration, check the reaction LHS and host for localized-bond
  simplicity, dative role and uniqueness, noncovalent uniqueness, aromatic participant uniqueness
  and overlap, multicenter participant and participant-set uniqueness, and repeated stereo sites.
  These remain tier-2 invariants generally, but application requires them because matching,
  incidence-induced correspondence, participant-position transport, and singular aromatic/stereo
  views rely on the corresponding uniqueness contracts. Do not encode them under a single
  “non-simple” condition.

  Rename `Reaction::validate_application` to `Reaction::check_preconditions`; it is an
  operation-precondition probe, not a tier-2 or tier-3 validator. `Reaction::apply` and the Python
  eager-correspondence path run this gate once before enumerating matches. The lower-level
  `Reaction::apply_at` does not repeat the global gate: its supplied `MoleculeCorrespondence`
  discharges the relevant matching ambiguity, while its lowering remains panic-free and checks the
  generated product. Keep that per-match product conflict gate separately, because a structurally
  valid LHS and host can still conflict under one embedding while another embedding succeeds.

  Graph IR uses the existing per-family conflict predicates for the operation precondition; the
  graph-layer `EntityStructureInvariantsValidator` remains responsible for detailed semantic
  diagnosis through `EntityStructureInvariantsContradiction` and
  `EntityStructureInvariantsError`. Report the coarse operation failures as
  `ApplyPreconditionError::ReactionStructureInvariant { kind: EntityKind }` and
  `ApplyPreconditionError::HostStructureInvariant { kind: EntityKind }`. Retain the existing
  `ApplyError::StructuralConflict` for the match-local product gate unchanged; its name belongs to
  the later application-error hierarchy review, and this work introduces no additional
  `*Conflict` name. Exercise every retained precondition, the distinction between global
  precondition and match-local product rejection, and independently supplied malformed inputs
  without panics. This is breaking red-to-green. [dep: S2b, S2c, S2e] **Done.**

### S3 — Complete remapping

- **S3a — Relation remapping and compaction names.** Rename public graph-core relation-set
  `apply_remapping` to `remap`, `try_apply_remapping` to `try_remap`, and `apply_compaction` to
  `compact` across all five relation-set forms and every caller. Preserve asserted coverage,
  checked `None`, removal-driven relation dropping, participant canonical ordering, and
  `RelationData`/`BiRelationData` permutation behavior. This is breaking red-to-green.
  [dep: doc 176] **Done.**
- **S3b — Dense molecule remapping.** Add public `Molecule::remap` and `try_remap` over a complete
  `MoleculeCorrespondence`. The checked route returns `None` when source counts differ, an
  entity-family mapping is not a bijection onto a dense target, or the source molecule fails its
  representation-integrity contract. Rebuild topology and every entity table through graph-core
  remapping plus `IdRemapping`; transport all constraint references and reuse relation payload
  permutation rather than reimplementing it. This is additive. [dep: S2a, S3a] **Done.**
- **S3c — Remapping semantic properties.** Add generated crossing permutations over all eight entity
  families, position-sensitive aromatic and multicenter data, stereo frames, and reference-bearing
  constraints. Validate `equiv_under`, exact identity, inverse roundtrip, composition agreement, and
  post-remap integrity. Keep asserted-producer tests separate from independently supplied coverage
  failures. This is additive. [dep: S3b] **Done.**

### S4 — Aggregate canonicalization contract and comparison schema

- **S4a — Public aggregate types.** Add `CanonicalizationContext`,
  `CanonicalizationLevel::{Topology, Constitution, Full}`, the aggregate `Canonicalize` trait, and
  the three carrier-specific `*CanonicalizationError` types. `canonicalize` is fallible;
  `canonical_eq` is total and implements the structural-equality, successful-form, contradiction,
  and integrity/operational-failure cases specified above. Add API-contract unit cases without
  implementing a fake canonicalizer. This is additive because S1 freed the name. [dep: S1c, S2b]
- **S4b — Stable typed comparison schema.** Define private typed comparison-key components with
  explicit entity-block, field, variant, span-tag, and constraint-extension positions. Implement
  ordering manually or through explicitly frozen component orders; do not use Rust discriminants,
  hashes, rendered DSL, or protocol bytes. Add exact ordering and append-only compatibility cases,
  including values with future-extension positions absent. This is additive. [dep: S1a]
- **S4c — Graph-layer operation inputs.** Add `CanonicalizationConfig` in `umol-graph` containing the
  automorphism algorithm and construct the graph-IR `CanonicalizationContext` from that config plus
  `StereoModel::para_stereo`. Do not add a duplicate one-field canonicalization model or move
  graph-IR reconstruction into `umol-graph`. This is additive. [dep: S4a]

### S5 — Exact common incidence encoding

- **S5a — Typed incidences.** Extend the existing `IncidenceGraph` to retain dative donor/acceptor
  roles, aromatic and multicenter participant positions, and stereo site and ligand-frame roles at
  every selected level. Keep one common incidence facility; do not add a canonicalization-only
  molecule graph. Add exact entity/node and typed-incidence recovery cases. This is additive to the
  carrier and may require a breaking constructor migration that ends green. [dep: S0d, S2a]
- **S5b — Collision-free initial classes.** Add collision-free equality-class ranking from normalized
  inherent values and typed incidences. Keep `MoleculeColoringFeatures` for consumers that genuinely
  need selectable hashed colors, but canonicalization never treats a `u64` hash as identity. Add
  cases proving equal represented values share a class and every selected distinction separates
  classes. This is additive. [dep: S1a, S4b, S5a]
- **S5c — Canonical-labeling adapter.** Translate the exact incidence carrier into the vertex-colored
  graph required by the selected automorphism backend, subdividing typed edges only at that adapter
  boundary when necessary. Return canonical labels and automorphisms in terms of the common
  incidence entities. This is additive. [dep: S5a, S5b]
- **S5d — Encoding verification and benchmark decision.** On bounded generated molecules, compare
  colored-encoding isomorphism with explicit dense remapping equivalence at each implemented level.
  Re-run S0's construction/labeling benchmark with the exact carrier and select the common
  representation on measured evidence. If the compact carrier wins, revise `IncidenceGraph` itself;
  do not retain two molecular encodings. This is additive or a measured breaking red-to-green
  carrier replacement. [dep: S0b, S5c]

### S6 — Topology canonicalization

- **S6a — Topology frame selection.** Use exact topology colors and canonical labels to derive a
  complete dense entity correspondence, choose the minimum typed topology key, and apply the public
  molecule remapping operation. Implement `canonicalize_by(Topology, context)` and
  `canonical_eq_by(..., Topology, context)`; the equality compares only the selected layer while the
  transformation returns the complete remapped molecule. This is additive. [dep: S3c, S4b, S5d]
- **S6b — Topology properties.** Validate exact idempotence of the selected topology layer,
  invariance under dense atom and bond renumbering, preservation of the complete input under the
  induced remapping, integrity, disconnected structures, and agreement across supported
  automorphism algorithms. Freeze the first topology canonical-number fixtures. This is additive.
  [dep: S6a]

### S7 — Constitution canonicalization

- **S7a — Non-stereo overlay frame selection.** Extend the canonical search and typed key to dative,
  aromatic, multicenter, and noncovalent entities, including participant roles and positional
  electron counts. Implement `canonicalize_by(Constitution, context)` and its selected-layer
  equality without using stereo or constraints to break constitution-level ties. This is additive.
  [dep: S6b]
- **S7b — Constitution properties.** Add renumbering, participant-permutation, inverse/composition,
  algorithm-agreement, and exact idempotence properties over all six non-stereo entity kinds.
  Include undetermined and non-literal inherent values and freeze constitution fixtures without
  changing the topology fixtures. This is additive. [dep: S7a]

### S8 — Full canonicalization without para-stereo refinement

- **S8a — Stereo frame action.** Extend the exact encoding and typed comparison key to stereo atoms
  and bonds. Reorder valid ligand frames only with the corresponding coset action; treat undetermined
  configurations as frame-invariant and rely on S2 integrity for arity and coset domains. This is
  additive. [dep: S2b, S7b]
- **S8b — One-pass full level.** Implement `canonicalize_by(Full, context)` for
  `para_stereo == false` using one stereo-sensitive refinement from the constitution partition.
  Preserve every stereo entity rather than perceiving, resolving, or validating it. This is
  additive. [dep: S5c, S8a]
- **S8c — Stereo properties.** Validate frame/coset covariance, exact idempotence of the selected
  full layer, renumbering invariance, meso cases, undetermined configurations, and algorithm
  agreement. Malformed stereo reaches the exact integrity-error cases and never a panic. This is
  additive. [dep: S8b]

### S9 — Para-stereo fixpoint and complete molecule API

- **S9a — Monotonic para-stereo refinement.** Implement the `para_stereo == true` fixpoint with no
  caller-selected iteration cutoff. Establish monotonic partition refinement and finite termination,
  then validate cases where constitution and stereo symmetry interact. This is additive. [dep: S8c]
- **S9b — Constraint placement and complete canonical form.** After full structural labeling,
  normalize entity and molecule constraints and use their frozen typed order only to select among
  remaining structurally equivalent frames. Constraints do not alter structural orbits. Preserve
  set-like conjunction semantics and transport every reference. This is additive. [dep: S4b, S9a]
- **S9c — `Canonicalize` for molecules.** Implement unqualified `canonicalize` and total
  `canonical_eq` for `Molecule`, using `Full` plus the context's para-stereo behavior and S9b's
  complete-key selection. Add contradiction and integrity totalization cases, exact canonical-number
  fixtures, full idempotence, renumbering invariance, `canonical_eq` equivalence laws, and
  cross-algorithm identity. This is additive. [dep: S4a, S9b]
- **S9d — Molecule benchmark gate.** Re-run the S0 corpus for topology, constitution, one-pass full,
  and para-stereo full operations. Record construction, refinement, labeling, key comparison, and
  remapping costs separately; optimization may change internals only if canonical fixtures remain
  identical. This is additive. [dep: S9c]

### S10 — Absolute stereo delta vocabulary

- **S10a — Rust delta core.** Remove `Apply`, `Swap`, and `Mirror` from `StereoAtomDelta` and
  `StereoBondDelta`; retain explicit before/after modifications only. Simplify delta application,
  inversion, normalization, remapping, composition, and reaction-span conversion. Migrate generated
  strategies and verify closure under inversion and faithful absolute before/after semantics. This is
  breaking red-to-green. [dep: S1c]
- **S10b — DSL, specification, and Python.** Remove the relative variants from reaction DSL parsing
  and rendering, the DSL specification, examples, fixtures, fuzz seeds, and Python bindings. Do not
  retain compatibility variants or a second core route. Higher-level generic “invert” operations
  remain separate future work. This is breaking red-to-green. [dep: S10a]
- **S10c — Reaction/span normalization properties.** Strengthen the roundtrip suite so
  `N(r) = r.to_reaction_span()?.to_reaction()` is idempotent, reproduces the span, and has the same
  `apply_at` failure or exact derivation on generated hosts and explicit correspondences. Generate
  materializable successes by construction and introduce one named defect for exact failure cases.
  This is additive. [dep: S2b, S10b]

### S11 — Reaction-span canonicalization

- **S11a — Span comparison schema.** Extend the private typed schema with fixed `Unchanged`, `Added`,
  `Removed`, and `Modified { lhs, rhs }` tags and ordered lhs/rhs components. Apply participant-frame
  actions to both modified values. Normalize constraint spans as set difference rather than a
  multiset and align that semantics with `ConstraintDelta` normalization. This is additive.
  [dep: S4b, S10c]
- **S11b — Span remapping and canonical frame.** Reuse the molecule incidence, exact class ranking,
  and remapping facilities for `EntitySpan<T>` values in the union namespace. Implement aggregate
  `Canonicalize` and the three level-specific operations for `ReactionSpan`; do not canonicalize
  either side independently. This is additive. [dep: S9c, S11a]
- **S11c — Span properties.** Validate exact canonical idempotence, invariance under every valid
  dense union-frame renumbering, preservation of lhs/rhs projections under induced side remappings,
  reversal distinction, integrity, and algorithm agreement. Freeze canonical-span fixtures with
  additions, removals, modifications, all entity kinds, and constraint-only changes. This is
  additive. [dep: S11b]

### S12 — Reaction canonicalization

- **S12a — Canonicalize through the span.** Implement `Reaction::canonicalize` as
  `to_reaction_span`, span canonicalization, then infallible `to_reaction`. Preserve the existing
  materialization contradiction boundary and wrap integrity/conversion causes in
  `ReactionCanonicalizationError`; do not add a direct LHS-only canonicalizer. This is additive.
  [dep: S11c]
- **S12b — Total reaction equality.** Implement total `Reaction::canonical_eq` with the agreed
  structural-equality, successful-form, contradiction, and operational/integrity-failure semantics.
  Ensure it compares complete reaction canonical forms rather than canonicalized LHS values. This is
  additive. [dep: S12a]
- **S12c — Reaction properties.** Validate exact idempotence, invariance under complete reaction
  renumbering, `C_R(N(r)) == C_R(r)`, the weakened reversal law
  `C_R(reverse(C_R(r))) == C_R(reverse(r))`, preservation of match domain and application result,
  and algorithm agreement. Include non-materializable and intrinsically contradictory cases with
  exact errors and no panics. This is additive. [dep: S10c, S12b]

### S13 — Surface audit and closeout

- **S13a — Public API and rustdoc audit.** Review exports for the context, levels, configs, errors,
  traits, integrity checks, and every public validator. State the semantic properties and failure
  domains in rustdoc; do not cite this dated discussion document from code. Audit Python only for
  surfaces changed by the migrations above; a new Python canonicalization API requires its own
  explicit binding design if it has not already been approved. This is additive or a final
  red-to-green export migration. [dep: S2f, S9c, S12c]
- **S13b — Repository-wide verification.** Run formatting, clippy, workspace tests, the full
  graph-IR and graph property targets at the agreed larger case count, conformance targets, affected
  fuzz builds, and the canonicalization benchmarks. Confirm that no old `Canonicalize`
  normalization names, relative stereo delta variants, `IncidenceNodeSelection`,
  `ConstitutionFeatures`, graph-IR validator modules, or
  `apply_remapping`/`try_apply_remapping`/`apply_compaction` spellings remain. This is additive.
  [dep: S13a]
- **S13c — Permanent documentation and status.** Update the DSL specification, current examples,
  nomenclature, data-type, and property-test guides to the implemented API; remove the dated doc-186
  TODO markers; record benchmark results and exact compatibility promises; then mark this document
  completed and update `000-status.md`. The whitepaper remains author-managed. [dep: S13b]

The critical path is S0 → S1 → S2a/S2b → S3 → S4 → S5 → S6 → S7 → S8 → S9 → S11 → S12 →
S13. S10 is independent after S1 and must complete before S11. S2c–S2f can proceed alongside the
remapping and benchmark work after the integrity contracts land, but S2f remains a semantic blocker
until the reaction-application preconditions are approved. No stage in the core path is deferrable.
The LHS-only reaction shortcut, connected-component optimization, alternate canonical-labeling
backend, higher-level relative stereo operations, and new Python canonicalization bindings are
explicitly deferred; none may change the frozen canonical comparison schema or canonical numbering.
