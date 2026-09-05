# Nomenclature guide

## Purpose

This is the living guide to terms coined or assigned a repository-specific meaning in umol. It is
normative for new public names and explanatory documentation once a term is recorded as settled. It
does not rename existing APIs by itself; a disagreement between the guide and existing code is a
separate migration task.

Prefer the established domain noun over a newly generalized synonym. Public names need not be
artificially parallel when the underlying constraints, entities, or operations have different names.
A shared name is useful only when the semantics and available operations are genuinely shared.

Use complete words in public identifiers and associated-type names; do not introduce clipped
abbreviations. In particular, a conversion trait's associated type is `Context`, not `Ctx`.
Established repository terms recorded by this guide, including `Ir`, `Dsl`, `Id`, and `Config`, do
not license new abbreviations by analogy.

## How to use this guide

- **Before naming something new**, check *Retired and discouraged* first. It is indexed by the wrong
  word, not the right one.
- **Every glossary entry is self-contained.** Nothing depends on the surrounding section, so landing
  mid-file by search is safe.
- **`Not:` lines name the confusable neighbours verbatim**, so searching for a near-miss term reaches
  the entry that corrects it.
- **`In code:` gives the identifiers**, which is the bridge from this prose to the API surface being
  edited.
- Headings name the **concept**. If a type is renamed, the `In code` line changes and the heading
  does not, so links and citations to this guide stay valid.

## Boundaries

Orientation, read once. The contrasts below are the ones most often collapsed; each participant has
its own glossary entry.

**Operations.** Three disjoint kinds: *resolution* fills undetermined state using a chemistry model;
*transformation* rewrites one resolved representation into another; *validation* checks integrity,
model-independent invariants, and model-dependent conformance without mutation. Kekulization,
aromatization, and charge delocalization are transformations, not resolver behaviour, because they
alter determined representation.

**Integrity and validation tiers.** *Integrity* is the graph-IR-owned construction contract and is
enforced by an integrity check. *Invariants* are model-independent semantic properties and
*conformance* is acceptance by a selected chemistry model; both are checked by validators in
`umol-graph`. Run integrity, invariants, and conformance in that order; within conformance, run the
cheap topology-only connectivity check before valence conformance. Do not give integrity and
semantic validation one API family merely because the tiers are ordered.

**Derivation and policy.** *Perception* produces a policy-free *derivation*, which may carry
*inconsistencies*. A *policy* maps a classified inconsistency to a *recovery action*. Policies belong
to resolvers only: perception, chemistry models, and validators have none.

**Choices.** A *model* decides which result is chemically accepted. A *config* decides how an
operation is performed. A *policy* decides what happens after an inconsistency is established. An
*algorithm* selects one implementation of one algorithmic problem. These are four different things
and none is a synonym for another.

**Determinacy.** *Undetermined* is stored lattice state on a value. *Underdetermined* is an operation
outcome. The words differ by two letters and by kind; the glossary keeps them apart deliberately.

**Id-space relations.** A *correspondence* pairs some entities in two declared carriers and is a
partial bijection. A *remapping* transports every source id into another id space and is therefore
total on its source, although its image may occupy only part of a larger target. A *compaction* is
the partial old-to-new mapping produced by removal: surviving ids have dense images and removed ids
have none. Partiality from pairing and partiality from removal are different semantics.

**Position-space actions.** A *permutation* reorders positions within one participant frame; it does
not pair, relabel, or remove entities. `Permutation` is the fixed-maximum, `Copy` carrier used by
bounded stereo algebra. `DynPermutation` is the arbitrary-degree carrier used by ordinary overlay
participant frames. Both act on integer positions and are distinct from every id-space relation
above.

## Suffixes

The generative core of this guide. Each family gains members as work lands, so a new type's suffix
should be chosen here rather than by analogy with whichever neighbour was read last. Counts are
approximate.

| Suffix | Means | Count | Crates |
| --- | --- | --- | --- |
| `*Algorithm` | selects one implementation of one algorithmic problem | 36 | graph-core (primitives), graph-ir, graph |
| `*Config` | composite operational configuration | 30 | graph-ir (ops), graph, io, py — never graph-core |
| `*Model` | semantic choices deciding chemical acceptance | 13 | graph |
| `*Policy` | maps a classified inconsistency to a recovery action | 11 | edn, graph, py |
| `*Kind` | unit-variant enum discriminating a closed set of alternatives | 11 | graph-ir, geometric, graph-core, msym, py |
| `*Features` | bitflag set of independently combinable switches | 1 | graph-ir |
| `*Level` | closed enum selecting one of several nested named layers | 2 | graph-ir |
| `*Constraint` | one assertable predicate over an entity | 6 | graph-ir, py |
| `*Constraints` | the container holding an entity's constraints | 9 | graph-ir, py — as `*ConstraintsForm`, because the container is lattice-shaped |
| `*Key` | identifies a constraint entry within a container | 13 | graph-ir, perm, py |
| `*View` / `*Views` / `*ViewMut` | borrowed accessor into a molecule; plural is the collection; `Mut` is the editing form | 37 / 15 / 16 | graph-ir, io, py |
| `*Delta` | an encoded change belonging to a reaction side | 18 | graph-ir, py |
| `*Update` | a field-level change applied to one entity | 10 | graph-ir |
| `*Defaults` | values used where an input states nothing | 13 | graph-ir, py |
| `*Overrides` | values that replace what an input stated | 10 | graph-ir |
| `*Entry` | one row of a table-shaped registry or format | 41 | several |
| `*Contradiction` | a semantic rejection, the `Solution::Contradictory` payload | 18 | graph-ir, graph |
| `*Mismatch` | two independently meaningful things that disagree | 3 | graph |
| `*Error` | operational or setup failure, the `Err` side | 58 | all |
| `*Validator` | validates a tier-2 invariant or tier-3 conformance property | 12 | graph after the integrity-check migration |
| `*Resolver` | fills undetermined state under a chemistry model | 5 | graph |
| `*Form` / `*Dsl` | lattice-shaped graph-IR value / boundary surface | 69 / 38 | graph-ir, py |

Families with no member yet, named in settled vocabulary and expected to gain one:
`*Failure` and transformers (see *Transformer naming* below).

### Stacking order

Suffixes compose, and 93 type names carry more than one. The order is fixed and each position answers
a different question:

```
<entity kind> <role> <plurality> <representation> <mutability>
    Atom      Constraint    s          Form
    Atom      Constraint               Key
    Atom      Editor                   View            Mut
    AromaticSystem          s          View
```

- **entity kind** — which of the eight, omitted when the type is not per-kind;
- **role** — `Constraint`, `Update`, `Delta`, `Editor`, `Defaults`, `Overrides`;
- **plurality** — a trailing `s` on the role marks the container, not "several of these";
- **representation** — `Form` (implements `Lattice`), `Dsl` (boundary surface), `View` (borrowed
  accessor), `Key` (identifies an entry);
- **mutability** — `Mut`, last, only on views.

`AtomConstraintsForm` is therefore correct in every position: one atom's constraint *container*, which
is itself lattice-shaped and so carries `Form` by the same rule every other lattice type does. A new
type should be assembled in this order rather than by analogy with the nearest neighbour.

The entity-kind aggregate types are the plurality marker with no role and no representation:
`AromaticSystems`, `MulticenterBonds`, `NoncovalentBonds`, `DativeBonds`, `StereoAtoms`, and
`StereoBonds` each name the container of one entity kind. They own what their storage shape cannot
state — which factor bears the participant frame, which is a site, and the entity kind's uniqueness key —
so they take a bare plural rather than a fourth position, exactly as the rule below requires.

**Two known tensions, neither to be "fixed".** The plural is the least visible marker in the stack —
a single trailing `s` against three CamelCase words. It is not to be replaced with `Store`, `Set` or
another suffix: the names are already long, and a fourth CamelCase position costs more legibility
than the `s` costs. What makes it learnable instead is that the `s` is positional and invariant — it
attaches to the role, never elsewhere, and always means the container, never "several of these".

The length itself comes from the entity kind appearing in both the module path and the type name
(`ir::constraint::atom::AtomConstraintForm`). Dropping it from the name is the idiomatic Rust fix and
would shorten roughly ninety types, but it makes bare imports collide across kinds — `use
atom::ConstraintsForm` and `use bond::ConstraintsForm` in one scope — which is why most Rust codebases
keep the prefix and silence the lint. Recorded as a considered trade, not an oversight.

### Constraint singular and plural

`*Constraint` is one assertable predicate. `*Constraints` is the container holding an entity's
constraints. The plural is never a synonym for "several constraints" in a signature; it names the
store.

**Not:** using the plural for a `Vec<Constraint>` parameter.
**In code:** `AtomConstraintForm` against `AtomConstraintsForm`; likewise for every entity kind.

### View and Views

A **`*View`** is a record bundling an index with the underlying data, so a consumer never assembles
an `(id, data, participants)` tuple by hand. A **`*Views`** is a *namespace*, not a collection: it
groups the per-relation accessors — `count`, `ids`, `iter`, `get`, `Index` — so they do not have to
be buried on `Molecule` itself. Adding an entity kind therefore adds a `*Views` namespace rather
than five more methods on the molecule.

`*ViewMut` is the mutable form; `*EditorView` and `*EditorViewMut` are the editor-scope bundles used
inside an edit session. Views also exist over derived things, not only entities: `GraphView`,
`RingView`, `RingViews`, `NeighborView`, `StereoLigandView`.

**Views are receivers, never arguments.** A function takes ids and the molecule, or takes the owned
`*Form`; it does not take a view. A view borrows its molecule and exists to be called *on*, so passing
one propagates a borrow through a signature that did not need it.

**Views are presentation facades.** A view never takes another view as an argument nor holds one
as state, and its implementation lives in functions of the molecule and entity ids beneath it — a
facade's implementation never builds another facade. Accessor chains that *return* narrower views
(`atoms().get`, `rings().atom`, `atom().constraints()`) are the namespacing mechanism: scope lives
in the receiver, not in name prefixes.

**Not:** the owned representation, which is the `*Form` type. A view does not survive its molecule.
**In code:** `AtomView`, `AtomViews`, `AtomViewMut`, `AtomEditorView`, `AtomEditorViewMut`.

### Delta and Update

A `*Delta` is an encoded change belonging to a reaction side. A `*Update` is a field-level change
applied to one entity. They are not interchangeable and exist per entity kind in both families.

**Not:** each other. `AtomDelta` is reaction-side change encoding; `AtomUpdate` is an entity edit.
**In code:** `AtomDelta`, `Deltas`, `ConstraintDelta` against `AtomUpdate`, `AtomFieldChange`.

### Defaults and Overrides

`*Defaults` supplies values where an input states nothing. `*Overrides` replaces values the input did
state.

**Not:** each other. The difference is whether the input spoke.
**In code:** `MoleculeDefaults`, `AtomDefaults` against `MoleculeOverrides`, `AtomOverrides`.

### Transformer naming — unsettled

Three patterns are in use for the same family:

- agent noun: `Kekulizer`, `Aromatizer`, each having grown `*Error` and `*Config`;
- verb phrase: `DelocalizeCharge`;
- target phrase: `ToExplicitHydrogens`, `ToImplicitHydrogens` (planned).

The house pattern for a thing that performs an operation is the agent noun — `Resolver`, `Validator`,
`Kekulizer`. The counter-argument is that `Kekulizer` and `Aromatizer` are configured engines with
their own error and config types, while `DelocalizeCharge` and the hydrogen transforms are
parameterless operations, which may justify two conventions rather than one.

The engine-adjacent side is settled by the *Operation names* entry: the agent noun names engines
only, and every run artifact (`*Error`, `*Config`, state) takes the verb stem. Open here remains
only the choice between agent-noun engines and verb-phrase parameterless transforms.

## Module names

A module exporting a result-object noun is named that noun (`coloring`, `embedding`,
`matching`, `value`); an operation module whose only defining type is an agent is named the
verb (`resolve`, `validate`, `transform`, `parse`, `edit`), with the agent noun inside
(`resolve::Resolver`). Do not mix agent-noun and verb module names within one operation
family.

A data module may take a concise adjective when that is the most direct description of its
content (`oriented`, `dynamic`). Related modules are named recognizably alike.

## Retired and discouraged

Indexed by the word not to use. Add a row whenever an entry's `Not:` line forbids a specific
spelling.

| Do not write | Write instead | Because |
| --- | --- | --- |
| `projection` for a stored constraint or the incidence category | `constraint`, `incidence constraint` | `projection` names an actual mapping between representations |
| `predicate` or `representation` for the stored object | `constraint` | the repository term for a possibly non-ground assertion |
| `undetermined` for an operation outcome | `underdetermined` | stored state versus outcome |
| `underdetermined` for stored lattice state | `undetermined` | as above |
| `config` for semantic acceptance choices | `model` | config is operational |
| `policy` for chemical acceptance | `model` | policy acts after acceptance is established |
| `validate_integrity` | `check_integrity` | tier-1 integrity is a graph-IR construction check, not semantic validation |
| `*Validator` in `umol-graph-ir` | `check_integrity` for tier 1, or move the validator to `umol-graph` for tiers 2 and 3 | validators are chemistry-layer semantic operations |
| `*Selection` for nested structural layers | `*Level` | selection does not state that the alternatives form an ordered, nested hierarchy |
| `*Features` for mutually exclusive nested presets | `*Level` | features are independently combinable switches |
| `try_*` solely because an operation returns `Option` or `Result` | the ordinary operation name | fallibility is part of the signature; `try_*` distinguishes a checked counterpart or restricted dispatch |
| agent-stem composites for run artifacts (`ResolverError`, `ValidatorError`, `KekulizerError`) | verb stem (`ResolveError`, `ValidateError`, `KekulizeError`) | errors, configs, and state belong to the run, not the engine |
| operation-noun composites for run artifacts (`KekulizationConfig`) | verb stem (`KekulizeConfig`) | the operation noun names a completed act, not a run's parameters |
| `CanonicalizeLevel`, `CanonicalizationLevel` | private `DescriptionLevel` | the internal hierarchy describes represented prefixes; public canonicalization is complete-only |
| `ground` for a molecule's or entity's chemistry-level determinedness | `concrete` | ground is the lattice term (constraint coordinates included) and risks the ground-state reading |

## Open issues

Findings raised while compiling the guide, each also stated in its entry. This index exists so they
can be reviewed without reading the whole file. Remove a row when the issue is decided; the entry
text changes at the same time.

| # | Issue | Where | Kind |
| --- | --- | --- | --- |
| 1 | Transformer naming has three patterns in use — agent noun (`Kekulizer`, `Aromatizer`), verb phrase (`DelocalizeCharge`), target phrase (`ToExplicitHydrogens`). Needs deciding before the planned members land. | *Transformer naming* | decision |
| 3 | One patch law, two spellings: `apply`/`diff` on `EntityPatch`, `update`/`difference_to` on the entity update surface. | *Patch algebra* | naming split |
| 4 | `umol-perm/src/coset.rs` line 1 says the coset space is `R\P`; the `CosetSpace` doc on line 25 says `P/R`. The prose establishes right cosets `Rσ`, so line 25 is wrong. | *Coset* | doc error |
| 6 | `*Failure` and the transformer family are named in settled vocabulary but have no members yet, so the conventions are untested. | *Suffixes* | latent |
| 7 | Six of eleven `Solution` accessors have no call sites outside `umol-utils`: `into_determined`, `into_data`, `into_contradiction`, `is_determined`, `is_underdetermined`, `contradiction`. API surface rather than nomenclature, but cheap to trim now. | *Solution* | unused surface |
| 8 | `Derivation`, `Failure` and `Projection` have no `In code` names. | those entries | incomplete |

## Glossary

Alphabetical by concept. Entry shape: definition, optional detail, then `Not:` and `In code:`.

### Algorithm

An **algorithm** enum selects a concrete implementation of one algorithmic problem. `*Algorithm` is
the suffix for algorithm selectors: `umol-graph-core` defines the graph-algorithm primitives, and
higher layers follow the same suffix.

**Not:** config (which may *contain* an algorithm selection), model, policy.
**In code:** `*Algorithm` enums in `umol-graph-core`, e.g. `CommonSubgraphEnumerationAlgorithm`,
`RelevantCycleEnumerationAlgorithm`; also higher-layer selectors such as `SubstructureMatchAlgorithm`.

### Application

**Application** realizes an operation against a concrete host. Edit application executes a complete
edit plan and publishes the result only when edit execution and the result's publication checks
succeed. Reaction
application matches a reaction rule against a host and emits one `ReactionDerivation` per successful
match.

**Not:** plan (which is derived without mutating), `react` (which intentionally discards the
derivation and emits only product components).
**In code:** `apply`, `apply_at`, `ReactionApplicationIter`.

### Attributes

An entity's **attributes** are its complete entity form excluding identity and participants. Rust and
Python fields, constructor arguments, and pattern-matching attributes use `attributes`. The DSL map
uses `:attrs` for the same payload. Use `form` for an arbitrary lattice value and `ir` for a converted
graph-IR value; do not extend `attributes` to values that are not complete entity payloads.

**Not:** `ast`, `type`, or the entity's participants.
**In code:** `attributes` on entity views, edits, undos, and delta variants; `:attrs` in the DSL.

### Canonical and canonicalize

**Canonicalize** selects a canonical entity-id and participant frame for a complete indexed graph IR
modulo the admissible remappings. It uses canonical labeling, transports every entity and reference
through the selected frame, applies the corresponding participant actions, and normalizes the
carried forms. The operation takes an explicit canonicalization context. The library's typed order
defines the selected frame; backend canonical labels may guide and prune the search but do not
define that order.

**Canonical equality** compares the complete canonical forms produced under the same context. It is
the search-based counterpart of `framed_eq_under`: the caller does not supply a remapping because
canonicalization selects the frame.

The public operation is complete-only: topology, constitution, and structure are private search
prefixes, not caller-selectable canonicalization levels. For a fixed umol release and context,
canonicalization is deterministic. During the 0.x
series, its typed order and canonical representatives may change between releases. A canonical form
is therefore not a persistent identifier unless a future API supplies an explicitly versioned
canonicalization profile.

Automorphism **colors** are opaque equality labels supplied to the graph algorithm. Canonicalization
uses `Color` in code for these labels and the normalized keys ranked into them. They are distinct
from the cells or equivalence classes of a refined partition and from a stereo class.

**Not:** *normalize*, which operates within an existing id and participant frame. Not *canonical
labeling* either: canonical labeling is the graph-algorithm component used to select the frame,
whereas aggregate canonicalization constructs the complete remapped graph IR.
**In code:** `Canonicalize`, `canonicalize`, `canonicalize_with_remapping`, `canonical_eq`.

### Class

A **class** is the stereo family a configuration belongs to: either a named permutation-group family
with its degree, or a coordination geometry whose proper-rotation group fixes the cosets. The class
determines the coset space, and therefore what a configuration index means.

Named families: `Symmetric(n)`, `Alternating(n)`, `Cyclic(n)`, `Dihedral(n)`. Geometries:
`Tetrahedral`, `CisTrans`, `Axial`, `SquarePlanar`, `TrigonalBipyramidal`, `Octahedral`.

**Not:** *kind*, which selects a parameterless semantic alternative. A stereo atom's kind is
`StereoKind`; its class is `ClassKey`.
**In code:** `ClassKey`, and the `class` field of the stereo `:attrs` payload.

### Combine

**Combine** forms the disjoint union of two structures, yielding one structure with two components.

Named `join` originally; renamed to avoid collision with the lattice join. Implemented and unlikely
to gain analogues, so it is recorded rather than generative.

**Not:** *join*, which is the lattice least upper bound.
**In code:** `combine`, `combine_all`, `combine_from`.

### Compaction

A **compaction** is the partial old-to-new id mapping produced by removal. A surviving id maps to
its position in the closed-up post-removal table; a removed id has no image.

Each component declares its finite source count and derives its result count from the removed ids.
Identity requires a declared source count. `compact_vec` applies that operation to a source-sized
value column while preserving survivor order.

It is layered on both axes, matching correspondence: `Compaction<Id>` over one id space,
`GraphCompaction` pairing the node and edge spaces, and `MoleculeCompaction` carrying a typed
compaction per entity kind so every stale reference can be updated or discarded consistently.

**`UndoCompaction`** is the inverse view of a `MoleculeCompaction` for rollback: it carries surviving
post-removal ids back into the pre-removal coordinate system, and removed entities are restored from
the explicit `Undo` payloads rather than from the mapping, because a compaction has no image for
them.

**Not:** a correspondence, whose unmatched ids remain members of their respective carriers; not a
remapping, which gives every source id an image and never expresses removal.
**In code:** `Compaction<Id>`, `GraphCompaction`, `MoleculeCompaction`, `UndoCompaction`;
`compact`, `compact_*`, `uncompact_*`.

### Completion

A **completion** is one admissible ground assignment for an entity's underdetermined attributes,
produced by a resolution phase from the chemistry model. It is represented by the entity's form in
its ground state — completion is a role a ground `AtomForm` plays, as pattern is a role a
`Molecule` plays — never by a dedicated primitive-valued type, which would reintroduce the retired
ground/pattern type split. The disjunction of completions for one atom is extensional, a vector of
forms; a single non-ground form would denote the Cartesian product of its fields and lose the
coupling between them. A phase emits the set of completions that survive narrowing; a later phase
selects among them, and an underdetermined verdict reports the survivors.

**Not:** a stored assertion — completions live in solver state and are never written to the
constraint channel; not the resolution result, which is the committed outcome; not a struct of
primitive fields.
**In code:** ground `AtomForm` members of `AtomCompletions`, the keyed carrier.

### Concrete

An entity is **concrete** when its chemical attributes — the inherent fields — are determined;
a molecule is concrete when every entity form is. The constraint channel does not bear on
concreteness: assertions are not attributes, so a concrete molecule may carry non-ground
assertions until discharge removes the determined-redundant ones. Concreteness is the
chemistry-level determinedness of a structure; lattice groundness is the algebraic notion over
the full product including the constraint coordinates, and the two decompose per form as
`is_ground() == is_concrete() && constraints.is_ground()`. The pattern/concrete asymmetry names
this axis: patterns are permanently non-concrete descriptions.

**Not:** *ground* — the lattice term, which includes the constraint coordinates and stays on
`Lattice::is_ground`; *complete* — the closure/completion family; the electronic *ground state*.
**In code:** `Molecule::is_concrete`, `is_concrete` on each entity form (exhaustive destructure
with `constraints: _`); `FingerprintError::NotConcrete`.

### Config

A **config** contains operational choices controlling how an operation is performed, including
algorithm selection and iteration limits. `*Config` is the suffix for composite configuration: in
`umol-graph-ir` it belongs to ops rather than to individual methods, and it is used throughout
`umol-graph`, `umol-io`, and `umol-py`. `umol-graph-core` defines no configs; its operations take
algorithm selectors directly.

**Not:** model, policy.
**In code:** `*Config`, e.g. `ResolveConfig`, `ValidateConfig`, `SubstructureSearchConfig`,
`ConstraintValidateConfig`, `SmilesIoConfig`.

### Configuration

A **configuration** is a ligand-to-position map. Two configurations describe the same arrangement
when a proper rotation relates them, so the observable descriptor is the coset rather than the map.

**Not:** *coset*, which is the equivalence class the configuration falls into; not *conformation*,
which umol does not represent.
**In code:** `StereoConfigurationForm`; the `configuration` attribute of stereo atoms and bonds.

### Conformance

**Conformance** asks whether a structure is accepted by a selected chemistry model. Conformance
validators may reuse policy-free derivations, but they convert every determined inconsistency into a
contradiction.

**Not:** integrity or an invariant. Conformance is tier-3 semantic validation in `umol-graph`.
**In code:** `validate_conformance`, `*ConformanceValidator`.

### Constitution

**Constitution** is the structural level containing topology (AB) plus the non-stereo domain (DAMN):
dative bonds, aromatic systems, multicenter bonds, and noncovalent bonds. It distinguishes
structural composition and connectivity without adding stereo configuration.

**Not:** all overlays, because stereo atoms and stereo bonds belong only to the structure level. Not
constraints, which do not contribute to structural identity.
**In code:** private `DescriptionLevel::Constitution`; public `IncidenceLevel::Constitution`.

### Constraint

A **constraint** is a possibly non-ground assertion represented by the constraint forms. This is
the public repository term. A constraint restricts the states admitted by an entity or molecule but
does not contribute to its identity. It therefore does not distinguish structural automorphism
orbits or select a canonical entity frame. It remains part of the complete IR assertion and is
transported and, where required, compared after a structural correspondence has been established.

**Not:** `projection`, `predicate`, or `representation`, when naming the stored object or an
operation over it. Not an inherent field, which does contribute to identity.
**In code:** `AtomConstraintForm`, `AtomConstraintsForm`, and the per-entity equivalents.

### Constraints view

A **constraints view** is the constraint reading of one entity, reached by accessor chaining
(`molecule.atom(id).constraints()`): the stored container's read API with its meanings intact —
the typed getters, `iter`, and `is_empty` read the asserted side — plus the keyed accessors
`asserted`, `derived`, and `derived_complete` and the comparisons `satisfies` and
`is_compatible`. Scope lives in the receiver, not in name prefixes: `atom(i).valence()` is the
derived quantity, `atom(i).constraints().valence()` the asserted payload. Mutation never routes
through a constraints view; it belongs to the stored container.

**Not:** the container (`*ConstraintsForm`), which is the storage the view reads; like every view
it is a receiver, never an argument.
**In code:** `AtomConstraintsView`, `BondConstraintsView`, and the views for every entity kind, from
`AtomView::constraints` and its peers.

### Contradiction

A **contradiction** is a semantic rejection represented by `Solution::Contradictory`. Validators have
no recovery policy: every determined failure or mismatch in their scope becomes a contradiction.

**Not:** error (operational, outside `Solution`), failure or inconsistency (policy-free
classifications that a resolver may still act on).
**In code:** `Solution::Contradictory`, `*Contradiction`.

### Correspondence

A **correspondence** is a partial bijection between two declared carriers' id spaces. It records
which entities are paired and the full size of each carrier; unmatched ids remain entities of their
own carrier rather than being interpreted as removed. A correspondence may relate a pattern and a
host, two reaction sides, or two id frames of the same semantic structure.

`MoleculeCorrespondence` holds one correspondence for each of the eight entity kinds. Its atom
component is a `Correspondence<AtomId>` aligned with the molecular graph, so the bond component can
be induced from the atom pairing and the two topologies. The remaining entity kinds carry their
own pairings.

A correspondence is **valueless** — it records pairing and nothing else. Adding values and a
direction is what lifts it to a reaction span.

Correspondences compose and reverse, which is what lets a chain of operations be followed end to
end. Reference transport consumes a correspondence directly; no correspondence-to-remapping
conversion is exposed.

**Not:** a compaction, because being unmatched does not mean that an entity was removed; not a
remapping, because a correspondence may be partial and records a relation rather than performing
transport.
**In code:** `MoleculeCorrespondence`, `Correspondence<T>`, `GraphCorrespondence`, `induce`,
`compose`, `reverse`, `left_of`, `right_of`, `is_total`, `map`, `try_map`.

### Coset

A **coset** is the observable stereo descriptor: the equivalence class of ligand-to-position maps
under the proper-rotation group of the class. Two configurations `σ` and `r∘σ` are the same
arrangement for any proper rotation `r`, so the class is the right coset `Rσ`, and the canonical
representative is its minimum-rank element.

A `CosetSpace` is `R\P` for the proper-rotation group `R` inside a parent group `P`, where `P` is the
group of realizable arrangements — `Sₙ` for the geometry classes, a partition subgroup for cis/trans,
where substituents are bonded to fixed sp² carbons. A *decomposition* numbers the cosets:
`CanonicalRank` is the generic Lehmer-minimum ordering, which reduces to the parity bit for a
two-coset space; the geometry decompositions reproduce the OpenSMILES arrangement numbers.

This is the standard group-theoretic term used in its ordinary sense; the entry exists because the
coset, not the configuration, is what the notation stores and what `#T` and `#C` carry.

**Not:** *configuration*, which is one representative map. Not *orbit* — the ligand permutation group
acts, but the stored descriptor is a coset of the rotation subgroup.
**In code:** `CosetSpace`, `coset_rep`, `observable_coset`, `orbit_reps`; the coset index in
`#T<n>` / `#C<n>`.

**Notation slip to fix:** `umol-perm/src/coset.rs` line 1 says the coset space is `R\P`, and the
`CosetSpace` doc comment on line 25 says `P/R`. The prose establishes right cosets `Rσ`, so `R\P` is
correct and the struct comment is wrong.

### Covalence

**Covalence** is the number of electrons an atom gains by sharing:

```text
covalence = valence + implicit_hydrogens + aromatic_covalence
aromatic_covalence = 1 if aromatic_valence == 1 else 0
```

The aromatic case is where the word earns its separate existence. **Aromatic valence** is what the
atom *contributes* to the π system; **aromatic covalence** is what it *gains back*. An atom
contributing one electron pairs with a neighbour's and gains one; an atom contributing a whole lone
pair — pyrrole-type nitrogen, furan oxygen — shares electrons it already owned and gains nothing.
Hence `{0, 2} → 0` and `1 → 1`.

Close to Irving Langmuir's sense, who introduced the term in 1919 as "the number of pairs of
electrons that a given atom shares with its neighbors". The two readings agree for localized bonds,
since each shared pair supplies exactly one electron from the partner, and diverge exactly where a
lone pair is donated into a delocalized system: Langmuir would still count the pair as shared, umol
counts nothing gained.

**Not:** *valence*, which counts localized bond orders and says nothing about what the atom receives.
Not *total valence* (`#V`), which sums contributions rather than gains. See *Valence* for the three
side by side.
**In code:** `AtomView::covalence`, `AtomView::aromatic_covalence`, `aromatic_covalence`;
`target_covalences` in the valence table.

### Deltas

**`Deltas`** is the resolved-delta collection carried by a reaction — the ordered transformation
applied to the left-hand side. The plural names the container, following the same convention as
`*Constraints`.

**Not:** a `Vec<Delta>` parameter. The plural type is the reaction's transformation, not any
sequence of deltas.
**In code:** `Deltas`, `Delta`.

### Derivation

**Derivation** carries two senses, and unlike *matching* this one is worth watching, because the
coined sense collides with the established one.

- **Reaction derivation** — one firing of a reaction: the two concrete molecule sides plus the
  correspondence between them. This is the standard sense in algebraic graph transformation, where a
  derivation is the result of applying a production, and it is the sense the whitepaper glossary
  uses. `ReactionDerivation` is `apply`'s codomain; a rule is to a derivation as a function is to one
  evaluation.
- **Perception derivation** — the policy-free result of perception, including candidates and exact
  inconsistencies.

**Not:** resolution (which applies policy), validation. The two senses are not related; the second
would be the one to rename if either were.
**In code:** `ReactionDerivation`, `apply`, `to_reaction`.

### Derived and asserted

A derivable constraint has an **asserted side** — the stored value, open-world by definition:
absence is the vacuous constraint — and a **derived side**, obtained by projection from the
relations the entity takes part in. The constraints views name the readings: `asserted(key)`
reads storage; `asserted_complete(key)` reads storage under resolution's closed-world claim —
absence of a constraint is actual absence, so the absence cell of an entity-creating overlay
key closes to its definite negative, and aromatic evidence merges both dialect placements
(the atom's own assertion or an incident bond's mark) without ever reading relations;
`derived(key)` reads present relations only and is vacuous on absence;
`derived_complete(key)` adds the closure — under the caller's claim that the relation set is
complete, absence of a resolution-written overlay closes to its definite negative. Topology keys
(valence, degree, the totals, the ring keys) read identically under both derived accessors and
have no absence cell on the asserted side; only the overlay-incidence keys have one. The two
closures are justified by different claims, so the four readings stay four flat accessors —
each consumer picks exactly one statically: matching reads `asserted` (open-world — a pattern
constrains only what it mentions), resolution's admission reads `asserted_complete`.

**Not:** `projected` for the reading — *projection* names the operation yielding the derived
side; not a mode parameter or stored molecule state — the closure choice is per call; no
parameterized accessor unifying the readings.
**In code:** `asserted`, `asserted_complete`, `derived`, `derived_complete` on the
constraints views.

### Description level

A **description level** is one member of the cumulative graph-IR description hierarchy:
`Topology`, `Constitution`, `Structure`, and `Full`. Each level includes every lower level.

The hierarchy is private canonicalization machinery. It selects the least search prefix containing
the complete input so empty higher domains need not enlarge the incidence graph. It is not a public
projection, a molecule property, or a reduced canonicalization surface: public canonicalization
always includes every present entity and constraint.

`Topology` contains the topology domain. `Constitution = Topology + NonStereo`, and
`Structure = Constitution + Stereo`. `Full` adds inline and molecule-level constraints. The domains
are disjoint groups of entity kinds; the levels are their cumulative prefixes.

**Not:** *model*, which decides semantic acceptance; *features*, which are independently combinable
switches; *domain*, which names a compositional entity group; *scope*, which does not express the
cumulative order; *representation level*, which compares distinct models; a public
*canonicalization level*, because the public operation is complete-only.
**In code:** private `DescriptionLevel` in aggregate canonicalization.

### Determined

**Determined** is an operation outcome: the pass produced a fully resolved result.

**Not:** *ground*, which is a property of stored state rather than of an outcome. An operation may
return `Determined` with a payload that is not ground if the operation's own contract is satisfied.
**In code:** `Solution::Determined`, `is_determined`, `into_determined`.

### Discharge

**Discharge** is the closing pass of the resolve pipeline that removes stored assertions the
structure has come to determine: per key, a ground derived value under the closure that refines
the assertion makes it redundant and removes it; a meet of `⊥` is a contradiction; a non-ground
derived value keeps the assertion. Discharge is the only removal concept for the constraint
channel. An operation that realizes an assertion may remove it early in its own transaction as an
implementation liberty; a transform that deletes a relation removes the assertions its output no
longer satisfies, which is the separate emit-compliance policy.

**Not:** cache invalidation — nothing derivable from structure is ever stored; not *consume* or a
per-operation contract.
**In code:** — (the planned closing `Resolver` stage).

### Donor and acceptor

The two roles of a dative bond: one or more **donors** contribute the electron pair, one
**acceptor** receives it. The asymmetry is the reason a dative bond is an overlay rather than an
edge — an undirected edge cannot carry a direction.

The order attribute counts donated pairs, not per-atom contributions, which is what distinguishes a
dative bond from a multicenter bond over the same atoms.

Donors form one participant frame and the acceptor is a distinguished second factor. All atoms in
the complete dative participant sequence are pairwise distinct: donors may not repeat, and the
acceptor may not also occur as a donor. A dative bond's identity is its complete
`(acceptor, donor multiset)` key; distinct dative bonds may share individual donors or their
acceptor.

**Not:** interchangeable, and not *participants* used flatly; a diagnostic should name the role.
**In code:** `:donors`, `:acceptor`; `DativeBondForm`.

### Edit and undo

An **`Edit`** is caller-facing mutation data: the symbolic vocabulary handed to
`MoleculeEditor::apply` or `MoleculeEditor::transact`. An **`Undo`** is realized rollback data
produced by the transaction path. The two are not two views of one thing — an edit says what was
asked for, an undo says what actually has to be reversed.

The symbolic references (`AtomHandle`, `BondHandle`, …) appear only inside `Edit`. `Id(_)` names an
existing entity, `New(n)` names the entity created by the nth edit earlier in the same batch, which
is why a handle is meaningless outside its batch.

**Not:** *delta*, which is the reaction-side change encoding, nor `*Update`, which is a field-level
change to one entity. Edits are molecule-level application data; only `transact` realizes undos.
**In code:** `Edit`, `Undo`, `UndoCompaction`, `MoleculeEditor::apply`,
`MoleculeEditor::transact`.

### Electron counts

**Electron counts** are the per-participant vector recording what each atom contributes to a
delocalized entity: π electrons for an aromatic system, shared electrons for a multicenter bond. The
vector is positional — entry *i* belongs to participant *i* — so reordering the participants must
permute it.

One leaf type serves both entity kinds rather than being duplicated per kind.

**Not:** the entity's total, which is a separate asserted constraint (`#e`) cross-checked against the
vector's sum. Not the entity's charge, which is carried on the entity and is why the five carbons of
a cyclopentadienyl ring can hold equal contributions.
**In code:** `ElectronCountsForm`, the leading `electron-counts` of an aromatic or multicenter string.

### Embedding kind

**Embedding kind** selects how strictly a correspondence must preserve non-edges, and is a required
argument to every common-subgraph and substructure operation.

- **Induced** — a pair of mapped nodes is adjacent in one graph exactly when it is adjacent in the
  other. An edge present on one side and absent on the other is inadmissible.
- **Monomorphism** — edges must map to edges, but the host may carry additional edges between mapped
  nodes.

The distinction is what Raymond and Willett call maximum common *induced* subgraph against maximum
common *edge* subgraph, and it decides whether the "exactly one source edge exists" case is
admissible in the modular product.

**Not:** an algorithm. `EmbeddingKind` says which answer is wanted; `*Algorithm` says how to compute
it.
**In code:** `EmbeddingKind::Induced`, `EmbeddingKind::Monomorphism`; also `McisAlgorithm` and
`McesAlgorithm`.

### Entity

An **entity** is one of the eight kinds represented by `Entity`: atom, localized bond, dative bond,
aromatic system, multicenter bond, noncovalent bond, stereo atom, or stereo bond. Use the concrete
entity name in diagnostics and action fields when it matters which kind is affected.

**Not:** *kind*, which is the discriminant naming an entity's category rather than an instance. Not
*overlay* either: entity covers all eight kinds, overlay covers the six that are not
topology.
**In code:** `Entity`, `EntityKind`.

### Entity constraint

An **entity constraint** is the broad category of constraints stored on or addressed to an entity. It
includes ring constraints.

**Not:** a name for the non-ring subset — use *incidence constraint* for that.
**In code:** `EntityConstraint` must not be narrowed to the non-ring subset.

### Equality ladder

The equality ladder follows the normalization, reframing, and canonicalization pipeline.
`framed_eq_under` is the explicit entity-id-witness form of `framed_eq`, not another quotient
level.

- **`==`** — derived structural equality of the stored IR. Same structure, constraints, ids, and order.
  Deliberately *not* chemical identity, so it stays cheap on the hot path.
- **`normalized_eq`** — equality of normal forms in the current entity-id and participant frame.
- **`framed_eq`** — equality after normalization and participant-frame selection.
  **`Molecule::framed_eq_under`** first remaps entity ids through an explicitly supplied
  `MoleculeRemapping`, then performs framed equality.
- **`canonical_eq`** — equality of complete aggregate canonical forms under a shared context. The
  implementation selects participant frames and entity ids rather than receiving an id witness
  from the caller.

For integrity-valid inputs whose complete canonicalization succeeds, `canonical_eq` holds exactly
when a remapping exists under which `framed_eq_under` holds. Equality
totalization for two intrinsic contradictions does not require such a witness.

Structural canonical labeling initially establishes automorphism orbits from inherent fields and
incidence without constraints. Complete aggregate canonicalization then uses normalized
constraints to select among structurally equivalent frames. This distinction is especially
important for patterns: constraints do not define the underlying structural orbits, but they remain
meaningful parts of the canonical IR assertion.

**Not:** each other. Reaching for `==` when `normalized_eq` or `framed_eq` is meant is the common
error, because `==` exists on everything and silently answers a different question.
**In code:** `PartialEq`, `Normalize::normalized_eq`, `Reframe::framed_eq`,
`Molecule::framed_eq_under`, and `Canonicalize::canonical_eq`.

### Error

An **error** is an operational or setup failure outside the semantic `Solution`, such as a failed
transaction or unavailable model parameters.

**Not:** contradiction, failure, inconsistency — all of which are semantic.
**In code:** the `Err` side of `Result<Solution<_, _>, _>`. Every module error type implements
`umol_utils::UmolError`, which supplies `as_any` for downcasting; `Box<dyn UmolError>` is the
cross-module boundary form and `?` promotes into it.

### Failure

A **failure** is a policy-free derivation classification. A constraint failure means that a
non-vacuous constraint cannot produce a valid entity under the selected topology and model. An entity
failure means that a structurally readable entity is not realizable under them.

**Not:** error (operational, not semantic); contradiction (a failure is not yet one — a resolver may
be configured to retain or remove the affected input).
**In code:** —

### Features

A **features** type is a bitflag set of independently combinable switches. Callers may select any
meaningful subset; the values do not represent steps in a single hierarchy.

Use a domain-qualified name when the switches describe one operation's inputs. The molecule-coloring
flags are therefore `MoleculeColoringFeatures`, not `ConstitutionFeatures`: the set also contains a
stereo-kind switch and is not limited to constitution.

**Not:** *level*, whose alternatives are mutually exclusive nested layers.
**In code:** `MoleculeColoringFeatures`.

### Form

A **`*Form`** is a graph-IR type that implements `Lattice`. The suffix marks exactly that property:
every Rust `*Form` implements `Lattice`, and a non-lattice type must not carry the suffix. Python
`*Form` classes bind the corresponding Rust forms and expose the same operations. Forms range from
leaf values such as `NumForm` to complete entity forms and constraint containers.

The non-lattice aggregate roots are `Molecule`, `Reaction`, and `ReactionSpan`. They use bare names
because they are complete graph-model objects, not lattice values. Role-bearing types such as
`MoleculeDsl`, `MoleculeDefaults`, `AtomUpdate`, and `AtomDelta` retain their role suffixes.

**Not:** an abstract syntax tree, a DSL boundary type, or a suffix for every graph-IR value.
**In code:** `AtomForm`, `AtomConstraintsForm`, `NumForm`; `Molecule`, `Reaction`, `ReactionSpan`.

### Full

**Full** is the private terminal description level: the complete molecular structure plus inline
and molecule-level constraints. It is not a public canonicalization selector.

`IncidenceLevel::Full` has a narrower carrier-specific use: it means every structural entity is
present in the incidence graph. Constraints are not incidence entities, so this carrier corresponds
to `DescriptionLevel::Structure`, not `DescriptionLevel::Full`.

**Not:** a chemistry model or conformance pass. The full description preserves the complete
represented assertion; it does not validate it.
**In code:** private `DescriptionLevel::Full`; public `IncidenceLevel::Full` in its
carrier-specific sense.

### Ground term

A **ground term** is a structure in which every inherent field holds a definite value. Groundness is
structural: it says the lattice is resolved to a bottom element, not that the structure satisfies
chemistry invariants or that its entities are mutually consistent.

`AsLit` is the exact projection out of a ground value, bound by the totality law
`value.is_ground() == value.as_lit().is_some()`. It does not normalize, apply defaults, validate,
or merge ground states that happen to have the same downstream numerical effect.

**Not:** valid, or chemically admissible. Structural groundness is separate from both. Not the chemistry-level determinedness of a molecule or entity, which is *concrete*.
**In code:** `Lattice::is_ground`, `AsLit::as_lit`, `Ground<T>` (planned).

### Graph IR

The **graph IR** is the explicit representation of the graph model. Its Rust package and public
module path are `umol-graph-ir` and `umol_graph_ir::ir`. Boundary conversions use `FromIr`,
`IntoIr`, `TryFromIr`, and `TryIntoIr`; their associated type is `Context`.

**Not:** the surface DSL, the format-level table IR, or a module retained under the renamed
crate.
**In code:** `umol_graph_ir::ir`; the `*Ir` conversion traits.

### Id, handle, and argument

Three ways to refer to an entity, at three stages.

- **`*Id`** — a resolved index into the entity table, a transparent `u32`. `AtomId`, `BondId`,
  `AromaticSystemId`.
- **`*Handle`** — a reference *within an edit batch*: either an existing `*Id`, or the Nth
  entity-creating edit earlier in the same batch. It refers to something that may not exist yet.
- **`*Arg`** — a builder argument that may create or refer: a spec creates, an integer selects by
  creation position, a name selects by name.

**Not:** each other. A `*Handle` outside a batch and an `*Arg` outside a builder are both meaningless.
**In code:** `AtomId`, `AtomHandle`, `AtomArg`; likewise per entity kind.

### Incidence

**Incidence** is the structural relationship between an entity and a relation it participates in: an
atom is incident with the bonds and overlays that name it. The word names the relationship, never the
value derived from it, and never the value's type — an incidence-derived value may be a Boolean, a
count, or a weighted sum.

The notion appears at three levels and they are the same idea:

- the `incident*` accessors on an entity view;
- the **incidence index** carried by a relation set, routed from every participant's `refs()`, which
  is what makes a relation reachable from any of its participants; in a birelation set the index
  spans both factors, so reachability does not depend on which id-space a participant belongs to;
- the **incidence graph**, where the relationship is made explicit as edges.

**Not:** *projection*, which is the operation that reads a value across an incidence. Incidence is
the relationship; projection is what you do with it.
**In code:** `incident*` methods, `IncidenceGraph`.

### Incidence constraint

An **incidence constraint** is an entity-local constraint whose value is derived from the entity's
fields and directly incident localized bonds or overlay relations, without running a separate graph
algorithm. The category includes:

- degree and valence aggregates over incident localized bonds;
- donated- and accepted-pair aggregates over incident dative bonds;
- aromatic-valence and aromatic-bond values from incident aromatic systems;
- multicenter-valence values from incident multicenter bonds;
- the corresponding directly incident stereo relations.

`Incidence` names the structural relationship, not the value type. An incidence-derived value may be
Boolean, a count, or a weighted sum.

**Not:** *relational constraint* (which has a narrower existing meaning); *projection*; a synonym for
the whole entity-constraint category.
**In code:** the incidence-constraint component of tier-2 invariants validation; the established
`incident*` methods expose the same relationship.

### Incidence graph

The **incidence graph**, also called the Levi graph, makes relations explicit as nodes: atoms occupy
the node space, every selected relation becomes a **pseudonode**, and each pseudonode is wired to its
participants. Everything defined over graphs then applies to the whole entity model unchanged —
canonical labeling, automorphisms, subgraph isomorphism — with no bespoke hypergraph algorithm.

Stereo elements attach to their site only, because the ligand topology is already present through the
bonds; the only new information a stereo node carries is its site and, at colouring time, its label.
Bond direction is not encoded structurally either: the colouring separates the endpoints of a
directed bond, so a dative donor and acceptor are never automorphism-equivalent.

`IncidenceLevel` names three carrier levels, which land on the chemist's own hierarchy: `Topology`
is atoms and localized bonds (AB), `Constitution` adds the non-stereo domain (DAMN), and `Full` adds
the stereo domain (SS). The final name means the complete incidence carrier; constraints have no
incidence nodes. The parallel description level for this carrier is `Structure`.

**Cost, and why it is not the default.** Because nauty does not accept edge colours, every localized
bond must also become a pseudonode, and a molecule has far more bonds than overlays. So
`SubstructureMatchAlgorithm::GraphAndOverlays` — match the topology, then verify overlays against the
correspondence — is the default, and `Incidence` is chosen when the connectivity a pattern turns on is
carried by an overlay rather than by bonds. Both return the same matches; the difference is pruning,
not semantics.

**Not:** the molecular topology, which is the atom-and-bond graph itself. The incidence graph is a
construction over a structure, not a representation of one.
**In code:** `IncidenceGraph`, `incidence_graph`, `IncidenceLevel`,
`SubstructureMatchAlgorithm::Incidence`.

### Inconsistency

An **inconsistency** is the policy-free umbrella classification for a constraint failure, entity
failure, or entity/constraint mismatch found during derivation. An inconsistency does not by itself
select an authority or recovery action.

**Not:** contradiction. A resolver policy may retain the inconsistency, remove one or both inputs
where sound, replace the entity, or report a contradiction.
**In code:** `AromaticityInconsistency`, `StereoInconsistency`, which identify the exact constraint
site and entity involved.

### Inherent field

An **inherent field** identifies an entity and is always present on it, though its value may be
undetermined.

**Not:** constraint, which restricts admitted states without contributing to identity and may be
absent.
**In code:** the non-constraint fields of each entity `*Form`.

### Integrity

**Integrity** is the tier-1 construction contract of a graph-IR representation: well-formed storage,
resolvable stored references, required parallel-collection shapes, fixed entity-relation semantics,
and kind-dependent values needed to interpret the representation. Constraint satisfaction and other
model-independent semantic conditions are invariants rather than integrity.

**Not:** an invariant or conformance. Integrity is established by construction and checked in
`umol-graph-ir`; it is not a `Solution` verdict.
**In code:** crate-private `check_integrity`, public `*IntegrityError`.

### Integrity check

An **integrity check** is the graph-IR-owned, error-valued operation that enforces tier 1. It returns
`Result<(), *IntegrityError>` and is shared by checked constructors, boundary conversions, and every
path that publishes the IR type. Trusted asserted constructors use the same implementation and
change only the failure reporting. The operation is crate-private: callers publish through the
checked or asserted boundary rather than validating an already published aggregate. There is no
public validator or `*Checker` object.

**Not:** a validator. Validators return semantic `Solution` values and belong in `umol-graph`.
**In code:** crate-private `Molecule::check_integrity`, `Reaction::check_integrity`, and
`ReactionSpan::check_integrity`; the public corresponding `*IntegrityError` types.

### Invariant

An **invariant** is a model-independent physical or mathematical condition over an otherwise
well-formed structure. Agreement between an independently meaningful stored constraint and its
model-independent derived value is also an invariant.

**Not:** integrity or conformance. Invariants are tier-2 semantic validation in `umol-graph`.
**In code:** `validate_invariants`, `*InvariantsValidator`, including
`ValenceInvariantsValidator` and `SpinInvariantsValidator`.

### Key and kind

A **kind** is a unit-variant discriminant enum naming an alternative without parameters; a
**key** is a discriminant that carries parameters. `AtomFieldKind` enumerates atom field keys;
`AtomConstraintKey` is a key because `RingMembership(RingScope)` carries its scope.
Peer entity forms gain `*FieldKind` enums as consumers arrive.

**Not:** interchangeable — a parameterless discriminant enum is a kind, not a key.
**In code:** `AtomFieldKind`; `AtomConstraintKey` and the per-entity-kind constraint keys.

### Lattice

**Lattice** is the internal term for the partial order on attribute values under which matching,
meet, and join are defined. The whitepaper calls it the *attribute lattice*; "attribute" is a
clarification added for chemist readers and is not repository terminology.

**Not:** *attribute lattice* in code or internal documentation.
**In code:** `Lattice`, `meet`, `join`, `matches`, `is_compatible`.

### Leaf type

A **leaf type** is a lattice-valued attribute that bottoms out in a concrete domain rather than in
other forms. Every leaf follows one shape: an `Undetermined` variant as the top of the lattice, a
`Lit` variant carrying a definite value, and whatever enrichments the domain admits.

`BooleanForm` is the minimal case, `Undetermined | Lit(bool)`. `NumForm` is the richest, adding
`LitSet`, `RangeFrom`, `RangeTo`, and the two expression arms. `ElectronCountsForm` is shared by
aromatic systems and multicenter bonds rather than duplicated per kind.

A new leaf type should follow the same shape, and should implement `AsLit` so that
`is_ground() == as_lit().is_some()` holds.

**Not:** an entity type, which is a record of leaves plus a constraint store.
**In code:** `BooleanForm`, `NumForm`, `ElectronCountsForm`, `UnpairedElectronsForm`, `ElementForm`,
`IsotopeMassForm`.

### Level

A **level** is one member of a closed enum of nested named layers. Selecting a level includes every
lower layer; the alternatives cannot be combined independently.

Private `DescriptionLevel` classifies cumulative parts of the graph-IR description for complete
canonicalization. Future structural entity kinds enter the first applicable level without changing
the meanings of earlier levels; future constraint variants enter `Full`.

`IncidenceLevel` stops at `Full`, using that name for the complete structural carrier because
constraints do not occur in an incidence graph. Thus `IncidenceLevel::Full` supplies the carrier for
`DescriptionLevel::Structure`.

The corresponding disjoint entity domains are topology (AB), non-stereo (DAMN), and stereo (SS).
They form the cumulative levels `Constitution = Topology + NonStereo` and
`Structure = Constitution + Stereo`; overlays are `NonStereo + Stereo`. A domain groups the entity
kinds introduced together, whereas a level includes that domain and every preceding domain.

**Not:** *features*, which are independently combinable bitflags; *selection*, which does not express
the nested relation.
**In code:** private `DescriptionLevel`; public `IncidenceLevel`.

### Ligand and site

A stereo entity has two role-distinguished participant groups.

- The **site** is the atom or bond bearing the configuration. A stereo atom's site is an atom; a
  stereo bond's site is a bond.
- The **ligands** are the ordered frame against which the configuration index is read. Order is
  significant: the coset is meaningless without it.

A ligand need not be an atom of the graph. A **virtual ligand** occupies a coordination position
without a corresponding node — an implicit hydrogen or a lone pair, each borne by a named atom. This
is how a sulfoxide or a pyramidal amine is expressed without materializing the hydrogen or inventing
an atom for the lone pair.

A published stereo frame contains pairwise-distinct complete ligand values and at most
`MAX_DEGREE` ligands. Two equal virtual ligands with the same anchor and kind are therefore
prohibited; an implicit hydrogen and a lone pair on the same anchor remain distinct. Explicit
hydrogens are atoms and remain distinct by atom id.

**Not:** *participant* used flatly. Site and ligands are both participants, but they are not
interchangeable, and diagnostics should say which.
**In code:** `StereoLigand`, `StereoLigandKind::{Atom, ImplicitHydrogen, LonePair}`; `:site` and
`:ligands` in the notation.

### Localized bond

A **localized bond** is a σ bond between exactly two atoms, carried as an edge of the molecular
topology. Unqualified, *bond* means localized bond.

Dative, multicenter and noncovalent bonds are overlay entities and are never called *bond* without
their qualifier.

**Not:** any overlay entity. The distinction is the topology/overlay split: a localized bond is an
edge, the others are relations.
**In code:** `BondForm`, `BondId`, `:bonds`.

### Matching

**Matching** carries two senses, both standard in their own literature and neither available to be
renamed. Qualify it where context does not.

- **Graph matching** — a set of edges with no shared endpoints. This is the combinatorial object:
  maximum, perfect and bipartite matchings, used by the kekulizer to assign double bonds.
- **Pattern matching** — finding embeddings of a pattern in a host, and separately the lattice
  relation `matches` where a pattern admits a value.

`BondMatching` is the first sense — a typed wrapper over `umol_graph_core::Matching` exposing matched
bonds and matched-atom membership — while `substructure_matches` and `SubstructureMatchAlgorithm` are
the second. The two sit close together in the module tree, so a reader arriving from substructure
search will read `BondMatching` as a match result unless told otherwise.

**`satisfies`** is the receiver-inverted reading of the lattice relation:
`target.satisfies(pattern)` iff `pattern.matches(target)`. It exists because a view must be a
receiver and stands on the target side; the `Lattice` default is the definition and is never
overridden.

**Not:** interchangeable, and neither is *correspondence*, which is the *result* of pattern matching
rather than the operation.
**In code:** first sense — `Matching`, `BondMatching`, `PerfectMatchingAlgorithm`,
`MatchingEnumerationAlgorithm`; second sense — `substructure_matches`, `SubstructureMatchAlgorithm`,
`Lattice::matches`.

### Mismatch

A **mismatch** means that a constraint and entity, or a constraint and its derived value, are each
independently meaningful but disagree.

**Not:** failure (where one side is not realizable at all).
**In code:** concrete names in diagnostic variants and config fields, such as
`AromaticValenceMismatch`, `CisTransStereoMismatch`.

### Model

A **model** contains semantic choices defining which result is chemically accepted. Named preset
constructors freeze a format's reading — `ValenceModel::smiles()` (the umol-owned SMILES table
with the `MostSaturated` tie-break), `ValenceModel::mdl()` (the frozen MDL table, likewise),
`AromaticityModel::daylight()` and `::mdl()`. Preset tables are frozen releases whose revisions
are new names; the living defaults promise only candidate-set monotonicity under additions.

**Not:** config, policy, algorithm.
**In code:** `ChemistryModel`, `ValenceModel`, `AromaticityModel`, `StereoModel`, `RingModel`.

### Molecular structure

A **molecular structure** is a molecular topology together with its overlay entities.

**Not:** molecular topology, which excludes them.
**In code:** `Molecule::is_concrete` reads the structure's determinedness.

### Molecular topology

The **molecular topology** is the attributed undirected simple graph of atoms and localized bonds,
carrying no aromatic, stereo, or coordination information.

**Not:** molecular structure, which includes the overlay entities.
**In code:** —

### Molecule atom

A **molecule atom** is an atom addressed within its molecule, passed as `&Molecule` plus `AtomId`
per the view rule. Operation names distinguish it from the owned form: `*_molecule_atom` takes the
molecule and the id, while bare `*_atom` takes an `AtomForm` outside any molecule context. The
pairs `check_atom`/`check_molecule_atom`, `for_atom`/`for_molecule_atom`, and
`resolve_atom`/`resolve_molecule_atom` follow this convention, which extends to every entity kind.

**Not:** "standalone atom" or "free-standing atom" for the form-level case — the bare operation
name with its `AtomForm` parameter already carries it.
**In code:** `*_molecule_atom` operations; bare `*_atom` operations on `AtomForm`.

### Narrow and widen

**Narrowing** moves a value down the attribute lattice and **widening** moves it up: the in-place
meet and the in-place join. Narrowing is what resolution does, and it returns whether the value
actually changed; widening returns `Err(NoJoin)` where no join exists.

These are the plain English words for descending and ascending an order of admitted-value sets, and
they are the right ones. `meet_with` and `join_with` would name the operation rather than its effect;
`refine_with` would collide with refinement, which is the order relation itself; `restrict` and
`relax` are vaguer; `intersect` and `union` are wrong for lattices that are not sets.

**Note for readers arriving from abstract interpretation.** That field uses both words for specific
operators that are *not* the lattice operations: widening (∇) deliberately over-approximates the join
so that iteration over an infinite-height lattice terminates, and narrowing (Δ) is the bounded descent
that recovers precision afterwards. Neither applies here — there is no fixpoint iteration and no
termination problem — so `widen_with` really is the join and `narrow_from` really is the meet.

**Not:** the abstract-interpretation operators of the same names. Not *refinement*, which is the order
relation rather than an operation.
**In code:** `Lattice::narrow_from`, `Lattice::widen_with`, `meet`, `join`.

### Noncovalent kind

The **kind** of a noncovalent bond is its interaction type, and it is the entity's only inherent
attribute beyond charge and spin: `HydrogenBond`, `HalogenBond`, `ChalcogenBond`, `Ionic`,
`VanDerWaals`.

**Not:** an order or a strength. Noncovalent bonds carry no order, which is one reason they are
overlays despite their binary shape.
**In code:** `NoncovalentBondKind`, `NoncovalentBondKindForm`; the notation literals `Hbd`, `Xbd`,
`Ybd`, `Ion`, `Vdw`.

### Non-stereo

**Non-stereo** is the DAMN entity domain: dative bonds, aromatic systems, multicenter bonds, and
noncovalent bonds. It excludes topology (AB) and stereo entities (SS). Constitution is topology plus
non-stereo; overlays are non-stereo plus stereo.

**Not:** constitution, which also includes topology; overlays, which also include stereo.
**In code:** the `NonStereo` position in the canonical comparison schema.

### Normalize

**Normalize** puts a form into a deterministic normal form without changing entity ids or
participant frames. It folds value expressions, normalizes set representations, flattens and
deduplicates logical constraints, normalizes entity fields and constraints, and normalizes
fixed-frame transformation values such as `Deltas`. It is context-free, idempotent on satisfiable
values, and returns `Err(Contradiction)` for an unsatisfiable represented value.

**Normalized equality** is `normalized_eq`: two values are equal in the current frame when their
normal forms are structurally equal. `Normalized<T>` carries the guarantee that normalization has
already run, so its derived `Eq`, `Hash`, and `Ord` operate on normal forms and can be used for
semantic deduplication.

**Not:** aggregate canonicalization, which selects an entity and participant frame and requires an
explicit context. Not chemical standardization, resolution, validation, or repair.
**In code:** `Normalize`, `normalize`, `normalized`, `Normalized<T>`,
`Normalize::normalized_eq`.

### Normalizer

The allowed permutation carrying a ligand frame to its least presentation. In code:
`CosetSpace::normalizer`, with `CosetSpace::allows` testing membership in the parent group. The
noun names the returned action, not an engine; the frame selection it performs is reframe-side
machinery despite the shared stem with *normalize*. **Not:** the group-theoretic normalizer
subgroup `N_P(R)`, which does not occur in code.

### Operation names

An operation family has up to three stems — the verb (*resolve*), the agent noun (*resolver*),
and the operation noun (*resolution*) — and composites choose by referent:

- The **agent noun** names the engine, its qualified engines, and agent classifications, and
  nothing else: `Resolver`, `ValenceResolver`, `ConstraintValidator`, `ParserType`.
- The **verb stem** names everything a run consumes, threads, or emits — configs, state,
  reports, errors, contradictions: `ResolveConfig`, `ResolveState`, `ResolveReport`,
  `ResolveError`, `ParseError`, `ValidateConfig`, `KekulizeConfig`. This is the default for
  operation-adjacent types.
- The **operation noun** is reserved for prose and for result objects that name a completed act
  as data (`AromaticityDerivation`); when the result object has a more specific name
  (`ResolveReport`, a coloring, a matching), prefer it.

Errors, configs, and state belong to the run, not to the engine that performs it — engines are
built *from* configs and *produce* reports, so `ResolverError` misattributes.

One level below the engine facade, method-qualified domain engines are named
`<Method><Domain>` with no agent suffix: `CountsValence`, `AtomTypingValence`,
`HueckelAromaticity`, `HmoAromaticity`, `ClarAromaticity`. The method qualifier does the
naming work; appending the agent noun to a dispatch variant's payload adds length and
stutters against the facade that selects it.

**Not:** agent-stem composites for run artifacts (`ResolverError`, `ValidatorError`,
`KekulizerError`); operation-noun composites for run artifacts (`KekulizationConfig`); agent
suffixes on method-qualified engines
(`HueckelAromaticityPerceiver`).
**In code:** `Resolver`/`ResolveConfig`/`ResolveState`/`ResolveReport`, `ValidateError`,
`KekulizeConfig`, `AromaticityPerceiver`, `CountsValence`.

### Overlay

An **overlay** is one of the six entity kinds that are not molecular topology: the non-stereo DAMN
entities plus the stereo SS entities. Atoms and localized bonds are the topology (AB) and are not
overlays.

**Not:** *relation* or *hyperedge*, which are whitepaper framings for the same thing and appear in
source comments descriptively; overlay is the repository term. Not *entity*, which is the umbrella
over all eight kinds.
**In code:** `GraphAndOverlays`, `verify_overlays`, `RemovedOverlays`.

### Participant

A **participant** is an entity referenced by an overlay entity — the atoms of an aromatic system, the
donors and acceptor of a dative bond, the site and ligands of a stereo element.

The word carries the same idea at two layers. In `umol-graph-core` a participant is a typed value
occupying a relation factor, routed through the incidence index; in the entity model it is the atom
or bond an overlay refers to. The lower layer is the mechanism for the upper one.

**Not:** member, constituent, or argument.
**In code:** `ParticipantPosition`, `RelationParticipant`.

### Patch algebra

A **delta** is the morphism between two entity states, and the pair of operations over it is a patch
algebra: `apply` carries a state forward by a delta, `diff` factors two states back into the deltas
between them.

The law is `apply(lhs, diff(lhs, rhs)) == rhs`. The entity update API states the same law as
`x.update(x.difference_to(y)) == y`. Both are read under `normalized_eq`, not `==`: a delta's payload is
compared up to normal form, so the law holds for values that are equivalent rather than identical.
Applying a delta checks its `old` against the stored state the same way.

**Naming split to resolve.** One law, two spellings — `apply`/`diff` on `EntityPatch`,
`update`/`difference_to` on the entity update surface. Pick one pair before either grows further; a
reader who learns the law under one name will not find it under the other.

**Not:** edit application, which executes a whole edit plan and publishes only on success. A patch
is an entity-level morphism; a transaction is a rollback-capable molecule-level lifecycle.
**In code:** `EntityPatch::apply`, `EntityPatch::diff`, `update`, `difference_to`.

### Perception

**Perception** applies a selected chemistry model and algorithms to identify or assess structural
entities. It produces a policy-free derivation.

**Not:** resolution (which applies policy and edits), validation.
**In code:** `AromaticityPerception`.

### Permutation and dynamic permutation

A **permutation** is a bijective action on the integer positions of one frame. Both carriers use a
one-line image with the direction `new[i] = old[action[i]]`; `between(from, to)` derives the unique
action carrying `from` to `to` and returns `None` when no unique action exists.

`Permutation` stores up to `MAX_DEGREE` positions in a fixed array and is `Copy`. It supports the
bounded permutation-group, coset, and stereo algebra. `DynPermutation` owns a `Vec<usize>`, has no
fixed maximum degree, and is not `Copy`; it carries actions for ordinary relation frames such as an
aromatic system with arbitrarily many participants. The two names distinguish representation and
degree, not two action conventions.

The image contains participant positions, never graph or entity ids. `ParticipantPosition` remains
the graph-core storage-facing position type and does not enter either action carrier.

**Not:** a remapping, which relabels every source id into another id space; a correspondence, which
pairs ids in two declared carriers and may be partial; or a compaction, which maps surviving ids and
has no image for removals.
**In code:** `Permutation`, `DynPermutation`, `MAX_DEGREE`, `FrameTransport::Action`.

### Plan

A **plan** is the complete edit sequence derived without mutating the source object.

**Not:** application (which executes a plan).
**In code:** `plan`.

### Policy

A **policy** is operational configuration mapping a classified inconsistency to a recovery action.
Within the chemistry layer, policies belong to resolvers, not perception, derivation, chemistry
models, or validators. The suffix is also used outside that layer for the same shape of decision —
`DuplicateKeyPolicy` in `umol-edn` maps a duplicate map key to `Error` or `LastWins` — so the
resolver restriction scopes where policies may live in `umol-graph`, not what the suffix means. A policy
enum is shared when the available action set is identical; the config field supplies the concrete
constraint or entity context.

**Not:** model. A model determines chemical acceptance; a policy determines what an operation does
after acceptance or inconsistency has been established.
**In code:** `AromaticityInconsistencyPolicy`, `StereoInconsistencyPolicy`.

### Projection

Use **projection** for an actual mapping from one representation or indexed relation to another, as
in pullback projections, reaction-side projections, or cycle projection from a subdivision graph. It
may be used descriptively when explaining that an entity relation induces participant values.

**Not:** the public name of stored constraints, or of the incidence-constraint category.
**In code:** —

### React

**React** is the product-oriented convenience operation over reaction application. For one molecule
it applies the reaction and splits each successful right-hand side. For several molecules it first
combines them by disjoint union in input order, then performs the same apply-and-split operation. It
emits one product-component collection per successful match and discards the derivation and split
correspondences.

**Not:** `apply`, which returns complete `ReactionDerivation` values; reaction construction or
composition.
**In code:** `React`, `react`, `react_all`, `ReactionProductsIter`.

### Reaction

A **reaction** is a left-hand-side molecule together with a resolved transformation. Everything else
— the atom map, the right-hand side, the condensed form, the reverse reaction — is *derived* from
`(lhs, deltas)` rather than stored.

The type is homoiconic in the sense the whitepaper claims: a molecule is the empty-deltas case, a
rule is a pattern left-hand side plus deltas, and applying a rule yields a concrete reaction of the
same type.

**Not:** a derivation, which is one firing of a reaction against a concrete host. A reaction is the
rule; a derivation is an evaluation of it.
**In code:** `Reaction`, `lhs`, `Deltas`.

### Reaction span

A **reaction span** is the superimposed `L ∪_K R` graph encoding a reaction's double-pushout rule
span. It carries, per atom and bond, both the before and after state plus a membership tag, and the
span `L ←K→ R` is read off those tags: `K = Unchanged ∪ Modified`, `L = K ∪ Removed`,
`R = K ∪ Added`.

`Modified` — a preserved entity relabeled across the reaction — is the relabeling-DPO reading: the
entity persists in `K` and its label is resolved per side. The tag asserts nothing beyond those two
side values. If they are `normalized_eq`, the entry is semantically a no-op whose normal form is
`Unchanged`. Raw span construction preserves an explicitly supplied tag; `normalize`, `reframe`,
and `canonicalize` collapse it, while `superimpose` may emit `Unchanged` directly because it derives
a standardized span.

A correspondence with values and a direction added is what lifts it to a span.

The materialized union is lhs-anchored: preserved entities retain lhs ids and rhs-only entities are
appended. Its lhs projection is structurally identical to the source lhs. Its rhs projection is the
source rhs reindexed into that reaction frame and is compared under the induced total
correspondence, not by structural equality when the source correspondence crosses entity order.

**Not:** a correspondence, which is valueless pairing; not a reaction, which is the rule itself.
**In code:** `ReactionSpan`, `lhs()`, `rhs()`.

### Recovery action

A **recovery action** is a policy variant naming what an operation does about a classified
inconsistency:

- `Keep` retains the affected inputs;
- `Remove` removes the target identified by the config field;
- `RemoveConstraint` retains the independently valid entity;
- `ReplaceEntity` atomically replaces the entity with the valid constraint-derived result;
- `RemoveBoth` removes both members of the pair identified by the mismatch policy;
- `Error` converts the classified inconsistency into a contradiction.

**Not:** separate enums created solely to repeat identical variants with longer target-specific
spellings. Policy enums may be shared when they admit exactly the same action set; the config field
supplies the context.
**In code:** the variants above.

### Refinement

**Refinement** names three unrelated operations. Always qualify it.

- **Lattice refinement** — the order relation on attribute values; `b` refines `a` when `a ∧ b = b`.
- **Colour refinement** — the automorphism and canonical-labelling procedure. `RefinementAlgorithm`,
  `umol-graph-core/src/algorithms/refinement.rs`.
- **Circular refinement** — the fingerprint iteration over atom environments.
  `CircularRefinementAlgorithm`, `CircularRefinementHash`.

Both graph refinements take **rounds**: `ToFixpoint` runs until the colouring stabilizes, `Fixed(n)`
runs exactly `n`. A fixed count is what makes a fingerprint reproducible; a fixpoint is what makes a
canonical form.

**Not:** any of the three used unqualified where another could be meant.
**In code:** `RefinementAlgorithm`, `CircularRefinementAlgorithm`, `RefinementRounds`,
`Lattice::matches`.

### Reframe

A **reframe** restates a frame-relative value in a different participant frame. It is the middle of
three nested quotients on a value: `normalize` reduces it, `reframe` reduces and then selects a
frame, `canonicalize` reduces, reframes and then selects ids. Their equalities nest the same way:
`==` refines `normalized_eq`, which refines `framed_eq`, which refines `canonical_eq`.

Pairwise alignment uses `DynPermutation::between` or `Permutation::between` to derive the unique
participant-frame action between two supplied frames, then applies it through
`FrameTransport::reframe_by`.
Integrity prohibits repeated complete participant values, so equal structured incidence determines
one action. `FrameTransport` is the transport-only operation for forms, form spans, delta payloads,
and constraints that can consume an action but do not select a frame. `Reframe` extends `Normalize`
for a frame-owning aggregate, derives its representative action, and selects the frame. A reaction
removal may carry another explicit local ordering of its source incidence: reaction transport
conjugates the owning action by the derived local-to-owner alignment, preserving the relation
between the removal frame and its owner. Normalization instead aligns the removal directly with the
owner before reframing selects another frame.

An action over an entity aggregate uses the aggregate's plural name followed by singular
`FrameAction`: `DativeBondsFrameAction`, `AromaticSystemsFrameAction`,
`MulticenterBondsFrameAction`, `NoncovalentBondsFrameAction`, `StereoAtomsFrameAction`, and
`StereoBondsFrameAction`. `OverlaysFrameAction` is the complete six-component action. The plurality
belongs to the carrier being acted on; `FrameAction` remains singular because the value is one
complete action.

Aggregate actions are operation-issued witnesses with private construction. Their identity,
inverse, and composition operations preserve an exact typed id-and-degree domain, with composition
defined only for equal domains. Consumption is receiver-relative: a frame-relative receiver needs
coverage only for the values it contains and ignores irrelevant action entries. A missing action,
wrong degree, or inadmissible entity-kind subgroup returns `None`.

`representative_action` and `reframe_with_action` materialize the complete input-domain witness for
a downstream consumer. Plain aggregate `reframe` fuses local action derivation and transport and
does not allocate that complete witness merely to discard it. `OverlaysFrameAction` covers
`Molecule`, `Reaction`, and `ReactionSpan` roots and their frame-relative constraints; a bare
`Deltas` value has no complete action because its owning frames live on the reaction.

**Not:** a remapping, which relabels ids across id spaces and does not touch order; not
canonicalization, which also selects ids.
**In code:** `Reframe`, `FrameTransport`, `representative_action`, `reframe_with_action`, `reframe`,
`reframe_by`, `framed_eq`, `OverlaysFrameAction`.

### Relation set

A **relation set** stores n-ary relations over typed participants — `NodeId`, `EdgeId`, or an
external type implementing `RelationParticipant` — with a shared incidence index so that a relation
is reachable from any of its participants.

Three axes parameterize a set, and the vocabulary is worth keeping straight:

- **arity** — `FixedRelationSet` has compile-time-known arity; `VarRelationSet` is variable-arity;
- **factor** — a birelation set relates *two* factors, each with its own participant type and arity.
  `FixedFixedBirelationSet`, `FixedVarBirelationSet`, `VarVarBirelationSet`.

A set stores the participant sequence it is given. The multiset is the relation's identity; the
stored sequence is the coordinate frame its payload is expressed in, and only graph IR interprets
that frame.

The incidence index of a birelation set spans both factors, so a relation is reachable from a
participant regardless of which id-space it belongs to.

**Not:** a relational-database table, though the whitepaper uses that reading to explain overlays.
Not the overlay itself: an overlay entity is *stored in* a relation set.
**In code:** `FixedRelationSet`, `VarRelationSet`, `FixedFixedBirelationSet`, `RelationParticipant`.

### Relational constraint

A **relational constraint** is a molecule-scope, reference-bearing constraint relating an overlay
entity to atoms, bonds, roles, or predicates.

**Not:** a synonym for *incidence constraint*, merely because the derived value depends on a
relation.
**In code:** `RelationalConstraint`.

### Remapping

A **remapping** is a total bijective old-to-new relabeling with a dense source id space. Every
source id has exactly one image, every target id has exactly one preimage, and no entity is added or
dropped. It is semantics-preserving alpha-renaming rather than a structural edit.

`umol_graph_core::Remapping<Id>` checks that its image vector is a permutation of its dense
source space. `GraphRemapping` aggregates validated node and edge remappings, and
`MoleculeRemapping` aggregates validated remappings for all eight molecule entity kinds.
Borrowed conversions to correspondences preserve all pairings and both counts. Agreement with
an independently supplied object's id-space sizes remains a contextual consumer requirement.

A correspondence describes a semantic remapping only when it is total on both sides and both id
spaces are dense. Relabeling relation participants preserves their stored sequence, so a positional
payload stays aligned without being touched.

**Not:** a correspondence, which may be partial and only records pairing; not a compaction, which
expresses removal by leaving removed source ids without images.
**In code:** `Remapping`, `GraphRemapping`, `MoleculeRemapping`, `remap`.

### Reset

**Reset** clears a constraint by setting it to its undetermined form.

**Not:** a general synonym for removing an entity or replacing determined structural information.
**In code:** `reset_aromatic_valence`, `reset_stereo_constraints`.

### Resolution

**Resolution** applies configured recovery policy, constructs edit plans against private drafts,
and publishes the result only when it is determined. An underdetermined, contradictory, or failed
resolution leaves its source unchanged.

**Not:** transformation (which rewrites determined representation), validation (which does not
mutate), perception (which is policy-free).
**In code:** `Resolver`, `resolve`, `ResolveConfig`.

### Result delivery

**Result delivery** is how a multi-result operation hands results to the caller, and it is
independent of algorithm selection — the `*Algorithm` enum says how to search, the delivery prefix
says how results arrive.

- **`visit_*`** — callback delivery; the visitor returns `ControlFlow`, and `Break` terminates the
  search.
- **`enumerate_*`** — eager collection of every result. Every eager operation returning a
  collection of results carries this prefix, whether or not a visitor form exists.
- **`iter_*`** — a resumable iterator with an explicit search cursor, added only when suspension
  across calls has a consumer.

An operation returning a single value — one output struct, one set, one coloring, a count, an
`Option` — takes none of these prefixes and has no visitor form; plural-sounding names for
single-value operations, such as `automorphisms`, are not delivery prefixes. In this family,
`try_*` marks restricted input-domain dispatch (simple graph against subdivision fallback), while
`*_fallback` marks the corresponding fallback route. Both compose with the delivery prefix. The
general meaning of `try_*` is defined under *Try prefix*. The execution contract behind the
prefixes — streamability, visitor payload, and emission order — is defined in the algorithm
execution guide.

**Not:** *iterative* for visitor delivery — an implementation may remain recursive while visiting
results. Not a bare plural name for an eager collection-returning operation.
**In code:** `visit_simple_cycles`, `try_visit_relevant_cycles`, `visit_perfect_matchings`,
`enumerate_relevant_cycles`; per the settled migration, `enumerate_subgraph_isomorphisms` and
`visit_subgraph_isomorphisms`.

### Ring constraint

A **ring constraint** is a constraint whose value requires ring enumeration. Ring membership, ring
degree, and ring valence belong here even though they are stored in atom, bond, or dative-bond
constraint containers.

Ring constraints are separate from incidence constraints because their evaluation has an algorithmic
dependency. Their fixed molecular semantics are the Relevant ring set through size 22; the selected
relevant-cycle enumeration algorithm is operational configuration.

**Not:** incidence constraint.
**In code:** the ring-constraint component of tier-2 invariants validation in `umol-graph`.

### Solution

**`Solution<T, C>`** is the three-valued outcome of an engine pass: `Determined`, `Underdetermined`,
or `Contradictory(C)` with a typed diagnostic payload. Setup and parameter failures travel separately
in `Result` and never collapse into it, so every umol engine returns `Result<Solution<_, _>, _>`.

An operation must decide which outcomes it treats as success, and the two reasonable answers differ
precisely on `Underdetermined`: a validator accepts it (only `Contradictory` fails), a transformer
does not (it cannot rewrite a representation on partial information). State the choice when adding an
operation; it is the outcome most easily overlooked. `into_observation` and `into_decisive` encode
the two, but with 4 and 1 call sites they are conveniences rather than established vocabulary.

**Not:** `Result`. `Solution` carries the semantic verdict; `Result` carries operational success.
Both appear in one signature and mean different things.
**In code:** `Solution`, `umol_utils::solution`.

### Snapshot

A **snapshot** is a non-consuming materialization of transient working state whose ordinary
finalization would consume the working object. Taking a snapshot preserves that object so subsequent
operations can continue from the same transient state. The implementation may clone, but the name
describes the lifecycle operation rather than the copying mechanism.

**Not:** ordinary finalization, which consumes the working object; transaction rollback, which
replays a realized undo journal; a general synonym for `clone`.
**In code:** planned `MoleculeEditor::snapshot`.

### Split

**Split** decomposes a structure into its connected components.

Implemented and unlikely to gain analogues; recorded rather than generative.

**Not:** a partition by any other criterion.
**In code:** `split`.

### Structural domain

A **structural domain** is one of three disjoint groups of entity kinds: topology (AB), non-stereo
(DAMN), or stereo (SS). The domains state which entity kinds enter the structural hierarchy
together. They are not cumulative: an atom belongs only to topology, a dative bond only to
non-stereo, and a stereo atom only to stereo.

The structure levels are cumulative prefixes of these domains:

```text
Topology     = topology
Constitution = topology + non-stereo
Structure    = topology + non-stereo + stereo
Full         = structure + constraints
Overlays     = non-stereo + stereo
```

**Not:** a structure level. `NonStereo` excludes topology even though `Constitution` includes both;
`Stereo` excludes constitution even though `Structure` includes both.
**In code:** `StructuralDomainPosition` and the `domain` field of `EntityBlockPosition` in the
canonical comparison schema.

### Structure

**Structure** is the complete structural description level: topology plus every overlay,
AB + DAMN + SS. It includes all inherent fields and typed incidences of every entity kind while
excluding constraints.

**Not:** `Full`, which additionally includes constraints; molecular topology, which excludes
overlays.
**In code:** private `DescriptionLevel::Structure`.

### Tie-break

A **tie-break** is the disposal policy for plural valence candidate survivors: a named
lexicographic key over atom fields — candidates ordered by each (field, direction) pair in
sequence, greatest selected — evaluated only after every model criterion has voted, aromaticity
selection included. `Strict` selects nothing and leaves survivors in the resolve report;
`MostSaturated` prefers max implicit hydrogens, then max lone pairs, then min unpaired
electrons. Policies are named variants only; keys are not openly constructible.

**Not:** a model criterion — it never preempts aromaticity; not a hidden preference — its use is
recorded per atom in the resolve report.
**In code:** `ValenceTieBreak`, `ValenceTieBreak::key`, `compare_by_key`; `SortingDirection` in
`umol_graph::utils`.

### Topology

**Topology** is the AB entity domain and the lowest structural level: atoms and localized bonds,
including the inherent values carried by those entities. It does not include overlay entities.

**Not:** constitution, which adds the non-stereo overlays; the incidence graph, which is an
algorithmic representation constructed from selected structure; constraints.
**In code:** private `DescriptionLevel::Topology`; public `IncidenceLevel::Topology`.

### Transaction

A **transaction** is the journal of realized undos for one batch of edits. `transact` applies a batch
and returns it; `rollback` consumes it to restore the prior state; `append` concatenates two so a
multi-stage operation can be reversed as a unit.

`transact` borrows an editor mutably, restores it on application failure, and returns a journal on
success. Editor `apply` consumes the editor and returns its modified state without producing a
journal; on failure, the consumed partial state is dropped. Molecule `apply` additionally publishes
the modified state as a checked `Molecule` while leaving its source unchanged, and distinguishes
transaction failure from failure of the molecule-integrity publication gate.

**Not:** *application*, which executes a plan and publishes only on success; a transaction is the
rollback mechanism for callers that retain an editor. Not *patch algebra*, which is the entity-level
`apply`/`diff` pair.
**In code:** `Transaction`, `transact`, `rollback`, `append`, `TransactionError`;
`MoleculeEditor::apply`, `Molecule::apply`.

### Transformation

A **transformation** explicitly rewrites one valid representation into another. Kekulization,
aromatization, and charge delocalization are transformations because they alter determined
representation rather than fill undetermined state.

**Not:** a resolver policy; not resolution.
**In code:** `umol-graph/src/ops/transform`.

### Try prefix

The **`try_*` prefix** distinguishes a checked operation from an unprefixed asserted or panicking
counterpart, or marks an attempt at a restricted input-domain dispatch for which the caller may use
a fallback. It does not merely announce that a method returns `Option` or `Result`. An operation
whose ordinary contract is fallible keeps its ordinary verb and expresses failure in its return
type, as `apply`, `parse`, `rollback`, `reverse`, and `to_reaction_span` do. Standard trait names
such as `TryFrom` and `TryInto` are unaffected.

Use checked/asserted pairs only when both public routes are useful. An internal producer that knows
a checked operation cannot fail handles that result explicitly rather than requiring a second
public method solely to avoid it.

**Not:** a general marker for fallibility.
**In code:** `map`/`try_map`, `remap`/`try_remap`, `from_entries`/`try_from_entries`, and
`try_visit_relevant_cycles` for restricted dispatch.

### Underdetermined

**Underdetermined** is an operation outcome: the available input does not justify a determined answer
or a contradiction.

**Not:** *undetermined*, which is stored state. Do not use the two interchangeably.
**In code:** `Solution::Underdetermined`.

### Undetermined

**Undetermined** is stored lattice state: a form asserts no concrete value at that position. An
absent or undetermined constraint is vacuous and does not assert that an entity is missing.

**Not:** *underdetermined*, which is an operation outcome. Do not use the two interchangeably.
**In code:** `NumForm::Undetermined` and the per-type equivalents.

### Valence

Three related counts, distinguished by what they measure.

- **Valence** — the sum of the orders of an atom's localized bonds to non-hydrogen neighbours.
  Excludes implicit hydrogens, dative, aromatic, multicenter and noncovalent contributions, each of
  which is a separate field. Tag `#v`.
- **Total valence** — valence plus implicit hydrogens plus aromatic valence plus multicenter
  valence. What the atom *contributes* in total. Tag `#V`, derived.
- **Covalence** — valence plus implicit hydrogens plus aromatic *covalence*. What the atom *gains*
  by sharing. No tag; it is a computed quantity rather than a constraint.

Total valence and covalence differ only in the aromatic term, and only for an atom contributing a
lone pair: pyrrole-type nitrogen has aromatic valence 2 and aromatic covalence 0.

**Not:** each other. The lowercase and uppercase tags are distinct predicates, and *covalence* is
neither. See *Covalence*.
**In code:** `#v`, `#V`, `AtomView::valence`, `AtomView::total_valence`, `AtomView::covalence`.

### Validation

**Validation** reports tier-2 invariant or tier-3 conformance outcomes without repairing or
selecting an authoritative representation. Validators live in `umol-graph`, use the
`*InvariantsValidator` or `*ConformanceValidator` suffix, and return
`Result<Solution<_, _>, _>` so non-ground semantic questions may be underdetermined.

**Not:** an integrity check, resolution, or transformation. Integrity is ordered before validation
but is enforced by graph-IR construction rather than by a validator.
**In code:** `Validator`, `validate`.

## Maintaining this guide

Add a term when it has a durable meaning that future code or documentation must preserve. Every
entry carries three fields:

- the definition, in one sentence where possible;
- `Not:` — the nearby terms it must not be used to mean, spelled out so a search for the wrong word
  reaches this entry;
- `In code:` — the public type or method names implied by the decision, or `—` if none yet.

Whenever a `Not:` line forbids a specific spelling, add the corresponding row to *Retired and
discouraged*.

Do not turn the guide into an inventory of every exported type. Ordinary Rust or chemistry terms need
entries only when umol narrows, distinguishes, or coordinates their meanings.

Keep the file loadable in one piece — under roughly a thousand lines. If it outgrows that, split by
domain deliberately rather than letting it sprawl past the point where a reader or an agent sees all
of it.

Record findings in *Open issues* as well as in the entry, so they can be reviewed without reading the
whole file. An entry states the issue where a reader will meet it; the index exists so nothing has to
be rediscovered. Remove the row and amend the entry together when an issue is decided.
