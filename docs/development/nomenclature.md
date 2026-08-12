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

## Suffixes

The generative core of this guide. Each family gains members as work lands, so a new type's suffix
should be chosen here rather than by analogy with whichever neighbour was read last. Counts are as of
2026-08-10.

| Suffix | Means | Count | Crates |
| --- | --- | --- | --- |
| `*Algorithm` | selects one implementation of one algorithmic problem | 36 | graph-core (primitives), graph-ir, graph |
| `*Config` | composite operational configuration | 30 | graph-ir (ops), graph, io, py — never graph-core |
| `*Model` | semantic choices deciding chemical acceptance | 13 | graph |
| `*Policy` | maps a classified inconsistency to a recovery action | 11 | edn, graph, py |
| `*Kind` | unit-variant enum discriminating a family | 11 | graph-ir, geometric, graph-core, msym, py |
| `*Features` | bitflag set of independently combinable switches | 1 | graph-ir |
| `*Level` | closed enum selecting one of several nested named layers | 2 | graph-ir |
| `*Constraint` | one assertable predicate over an entity | 6 | graph-ir, py |
| `*Constraints` | the container holding an entity's constraints | 9 | graph-ir, py — as `*ConstraintsForm`, because the container is lattice-shaped |
| `*Key` | identifies a constraint slot within a container | 13 | graph-ir, perm, py |
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
  accessor), `Key` (identifies a slot);
- **mutability** — `Mut`, last, only on views.

`AtomConstraintsForm` is therefore correct in every position: one atom's constraint *container*, which
is itself lattice-shaped and so carries `Form` by the same rule every other lattice type does. A new
type should be assembled in this order rather than by analogy with the nearest neighbour.

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

This needs a decision before the planned transformers land, not after.

## Retired and discouraged

Indexed by the word not to use. Add a row whenever an entry's `Not:` line forbids a specific
spelling.

| Do not write | Write instead | Because |
| --- | --- | --- |
| `projection` for a stored constraint or the incidence category | `constraint`, `incidence constraint` | `projection` names an actual mapping between representations |
| `predicate` or `representation` for the stored object | `constraint` | the repository term for a possibly non-ground assertion |
| `EntityConstraint` for the non-ring subset | `incidence constraint` | entity constraints include ring constraints |
| `relational constraint` for an incidence constraint | `incidence constraint` | `RelationalConstraint` already means a molecule-scope, reference-bearing constraint |
| `undetermined` for an operation outcome | `underdetermined` | stored state versus outcome |
| `underdetermined` for stored lattice state | `undetermined` | as above |
| `reset` for removing or replacing an entity | `remove`, `replace` | reset clears a constraint to its undetermined form |
| `config` for semantic acceptance choices | `model` | config is operational |
| `policy` for chemical acceptance | `model` | policy acts after acceptance is established |
| `validate_integrity` | `check_integrity` | tier-1 integrity is a graph-IR construction check, not semantic validation |
| `*Validator` in `umol-graph-ir` | `check_integrity` for tier 1, or move the validator to `umol-graph` for tiers 2 and 3 | validators are chemistry-layer semantic operations |
| `*Selection` for nested structural layers | `*Level` | selection does not state that the alternatives form an ordered, nested hierarchy |
| `*Features` for mutually exclusive nested presets | `*Level` | features are independently combinable switches |
| `iterative` for visitor delivery | visitor delivery, `visit_*` | the implementation may remain recursive while visiting results |
| a bare plural method name for an eager collection-returning operation | `enumerate_*` | delivery must be legible at the call site; the bare plural is reserved for single-value operations |
| `*Ast` for graph-IR values | `*Form` for a `Lattice` type; bare aggregate name otherwise | the suffix records lattice semantics, not syntax-tree shape |
| `ast` as the graph-IR module or entity payload member | `ir` for the module; `attributes` for the payload | the graph model is an IR, and the payload is the entity's complete attribute form |
| `:type` for an entity payload in the DSL | `:attrs` | the payload contains attributes, not the entity kind |
| `Ctx` in a public identifier | `Context` | public identifiers use complete words |
| `apply_remapping`, `try_apply_remapping` | `remap`, `try_remap` | the receiver is transported through the supplied remapping |
| `apply_compaction` | `compact` | the receiver is transported through the supplied compaction |

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

**Application** executes a complete edit plan transactionally and publishes the result only when
every edit succeeds.

**Not:** plan (which is derived without mutating), transformation.
**In code:** `apply`, `apply_at`.

### Attributes

An entity's **attributes** are its complete entity form excluding identity and participants. Rust and
Python fields, constructor arguments, and pattern-matching attributes use `attributes`. The DSL map
uses `:attrs` for the same payload. Use `form` for an arbitrary lattice value and `ir` for a converted
graph-IR value; do not extend `attributes` to values that are not complete entity payloads.

**Not:** `ast`, `type`, or the entity's participants.
**In code:** `attributes` on entity views, edits, undos, and delta variants; `:attrs` in the DSL.

### Canonical and canonicalize

> **TODO (2026-08-07):** The definitions below are the approved target semantics from discussion
> doc 186. Fixed-frame form normalization now uses `Normalize` and `Equiv`; the aggregate
> canonicalization API remains to be implemented. Remove this marker when that API lands.

**Canonicalize** selects a canonical entity-id and participant frame for a complete indexed graph IR
modulo the admissible remappings. It uses canonical labeling, transports every entity and reference
through the selected frame, applies the corresponding participant actions, and normalizes the
carried forms. The operation takes an explicit canonicalization context. The library's typed order
defines the selected frame; backend canonical labels may guide and prune the search but do not
define that order.

**Canonical equality** compares the complete canonical forms produced under the same context. It is
the search-based counterpart of `equiv_under`: the caller does not supply a correspondence because
canonicalization selects the frame.

**Not:** *normalize*, which operates within an existing id and participant frame. Not *canonical
labeling* either: canonical labeling is the graph-algorithm component used to select the frame,
whereas aggregate canonicalization constructs the complete remapped graph IR.
**In code:** `Canonicalize`, `canonicalize`, `canonical_eq`.

### Class

A **class** is the stereo family a configuration belongs to: either a named permutation-group family
with its degree, or a coordination geometry whose proper-rotation group fixes the cosets. The class
determines the coset space, and therefore what a configuration index means.

Named families: `Symmetric(n)`, `Alternating(n)`, `Cyclic(n)`, `Dihedral(n)`. Geometries:
`Tetrahedral`, `CisTrans`, `Axial`, `SquarePlanar`, `TrigonalBipyramidal`, `Octahedral`.

**Not:** *kind*, which discriminates entity families. A stereo atom's kind is `StereoKind`; its class
is `ClassKey`.
**In code:** `ClassKey`, and the `class` field of the stereo `:attrs` payload.

### Combine

**Combine** forms the disjoint union of two structures, yielding one structure with two components.

Named `join` originally; renamed to avoid collision with the lattice join. Implemented and unlikely
to gain analogues, so it is recorded rather than generative.

**Not:** *join*, which is the lattice least upper bound.
**In code:** `combine`, `combine_all`, `combine_from`.

### Compaction

A **compaction** is the partial old-to-new id mapping produced by removal. A surviving id maps to
its position in the closed-up post-removal table; a removed id has no image. `IdCompaction` wraps
`umol_graph_core::Compaction` for atoms and bonds and carries the removed ids for the six relation
kinds so every stale reference can be updated or discarded consistently.

**`UndoCompaction`** is the inverse view of an `IdCompaction` for rollback: it carries surviving
post-removal ids back into the pre-removal coordinate system, and removed entities are restored from
the explicit `Undo` payloads rather than from the mapping, because a compaction has no image for
them.

**Not:** a correspondence, whose unmatched ids remain members of their respective carriers; not a
remapping, which gives every source id an image and never expresses removal.
**In code:** `IdCompaction`, `UndoCompaction`, `umol_graph_core::Compaction`; `compact`,
`compact_*`, `uncompact_*`.

### Completion

A **completion** is one admissible ground assignment for an entity's underdetermined attributes,
produced by a resolution phase from the chemistry model. `AtomCompletion` pairs the inherent-field
completion (`AtomFields`: implicit hydrogens, lone pairs, unpaired electrons, spin) with the model
values that selection votes on (valence, donated and accepted pairs, aromatic and multicenter
valence). A phase emits the set of completions that survive narrowing; a later phase selects among
them, and an underdetermined verdict reports the survivors.

**Not:** a stored assertion — completions live in solver state and are never written to the
constraint channel; not the resolution result, which is the committed outcome.
**In code:** `AtomCompletion`, `AtomFields`.

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
**In code:** the `Constitution` variants of `IncidenceLevel` and `CanonicalizationLevel`.

### Constraint

A **constraint** is a possibly non-ground assertion represented by the constraint forms. This is
the public repository term. A constraint restricts the states admitted by an entity or molecule but
does not contribute to its identity. It therefore does not distinguish structural automorphism
orbits or select a canonical entity frame. It remains part of the complete IR assertion and is
transported and, where required, compared after a structural correspondence has been established.

**Not:** `projection`, `predicate`, or `representation`, when naming the stored object or an
operation over it. Not an inherent field, which does contribute to identity.
**In code:** `AtomConstraintForm`, `AtomConstraintsForm`, and the per-entity equivalents.

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

`MoleculeCorrespondence` holds one correspondence for each of the eight entity families. Its atom
component is a `Correspondence<AtomId>` aligned with the molecular graph, so the bond component can
be induced from the atom pairing and the two topologies. The remaining entity families carry their
own pairings.

A correspondence is **valueless** — it records pairing and nothing else. Adding values and a
direction is what lifts it to a reaction span.

Correspondences compose and reverse, which is what lets a chain of operations be followed end to end,
and `to_remapping` converts one into a total-on-source remapping when it is total on the left. The
result may map into a larger target id space. End-to-end remapping of a standalone `Molecule`
requires the stronger condition that every entity-family correspondence is total on both sides, so
the target tables are dense and contain exactly the mapped entities.

**Not:** a compaction, because being unmatched does not mean that an entity was removed; not a
remapping, because a correspondence may be partial and records a relation rather than performing
transport.
**In code:** `MoleculeCorrespondence`, `Correspondence<T>`, `GraphCorrespondence`, `induce`,
`compose`, `reverse`, `left_of`, `right_of`, `is_total`, `to_remapping`.

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

### Determined

**Determined** is an operation outcome: the pass produced a fully resolved result.

**Not:** *ground*, which is a property of stored state rather than of an outcome. An operation may
return `Determined` with a payload that is not ground if the operation's own contract is satisfied.
**In code:** `Solution::Determined`, `is_determined`, `into_determined`.

### Donor and acceptor

The two roles of a dative bond: one or more **donors** contribute the electron pair, one
**acceptor** receives it. The asymmetry is the reason a dative bond is an overlay rather than an
edge — an undirected edge cannot carry a direction.

The order attribute counts donated pairs, not per-atom contributions, which is what distinguishes a
dative bond from a multicenter bond over the same atoms.

**Not:** interchangeable, and not *participants* used flatly; a diagnostic should name the role.
**In code:** `:donors`, `:acceptor`; `DativeBondForm`.

### Edit and undo

An **`Edit`** is caller-facing mutation data: the symbolic vocabulary handed to
`MoleculeEditor::transact`. An **`Undo`** is realized rollback data produced by the checked
transaction path. The two are not two views of one thing — an edit says what was asked for, an undo
says what actually has to be reversed.

The symbolic references (`AtomHandle`, `BondHandle`, …) appear only inside `Edit`. `Id(_)` names an
existing entity, `New(n)` names the entity created by the nth edit earlier in the same batch, which
is why a handle is meaningless outside its batch.

**Not:** *delta*, which is the reaction-side change encoding, nor `*Update`, which is a field-level
change to one entity. Edits are molecule-level and transactional.
**In code:** `Edit`, `Undo`, `UndoCompaction`, `MoleculeEditor::transact`.

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

**Not:** *kind*, which names what distinguishes one entity family from another rather than an
instance. Not *overlay* either: entity covers all eight kinds, overlay covers the six that are not
topology.
**In code:** `Entity`, `EntityKind`.

### Entity constraint

An **entity constraint** is the broad category of constraints stored on or addressed to an entity. It
includes ring constraints.

**Not:** a name for the non-ring subset — use *incidence constraint* for that.
**In code:** `EntityConstraint` must not be narrowed to the non-ring subset.

### Equality ladder

Three levels of equality exist on forms and molecules. `equiv_under` is the mapped form of
`equiv`, not a fourth relation.

- **`==`** — derived structural equality of the stored IR. Same structure, constraints, ids, and order.
  Deliberately *not* chemical identity, so it stays cheap on the hot path.
- **`equiv`** — equality of normalized forms in the current id and participant frame.
  **`equiv_under`** is the same after reindexing the receiver into the other's frame via an explicit
  correspondence or participant order. The work is skipped when the payload is
  permutation-invariant.
- **`canonical_eq`** — equality of complete aggregate canonical forms under a shared context. The
  implementation selects the canonical frame rather than receiving one from the caller.

For inputs in the canonicalization domain, `canonical_eq` holds exactly when an admissible remapping
exists under which `equiv_under` holds. The two-factor frame-aware operation for birelation payloads
reindexes each factor before comparing.

Structural canonical labeling initially establishes automorphism orbits from inherent fields and
incidence without constraints. Complete aggregate canonicalization then uses normalized
constraints to select among structurally equivalent frames. This distinction is especially
important for patterns: constraints do not define the underlying structural orbits, but they remain
meaningful parts of the canonical IR assertion.

**Not:** each other. Reaching for `==` when `equiv` is meant is the common error, because `==` exists
on everything and silently answers a different question.
**In code:** `PartialEq`, `Equiv::equiv`, the frame-aware `equiv_under` traits, and
`Canonicalize::canonical_eq`.

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

**Full** is the terminal aggregate canonicalization level: the complete molecular structure plus
normalized entity-level and molecule-level constraints. Constraints enter only after the minimum
structure key has been established and distinguish only tied structural frames. `Full` is exactly
equivalent to the unqualified canonicalization operation, so callers may opt into the complete
operation either directly or through the level selector.

`IncidenceLevel::Full` has a narrower carrier-specific use: it means every structural entity is
present in the incidence graph. Constraints are not incidence entities, so this carrier corresponds
to `CanonicalizationLevel::Structure`, not to its `Full` constraint-selection step.

**Not:** a chemistry model or conformance pass. Full canonicalization preserves and orders the
complete represented assertion; it does not validate it.
**In code:** `CanonicalizationLevel::Full`; `IncidenceLevel::Full` in its carrier-specific sense.

### Ground term

A **ground term** is a structure in which every inherent field holds a definite value. Groundness is
structural: it says the lattice is resolved to a bottom element, not that the structure satisfies
chemistry invariants or that its entities are mutually consistent.

`AsLit` is the exact projection out of a ground value, bound by the totality law
`value.is_ground() == value.as_lit().is_some()`. It does not normalize, apply defaults, validate,
or merge ground states that happen to have the same downstream numerical effect.

**Not:** valid, or chemically admissible. Structural groundness is separate from both.
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
**In code:** `incident*` methods, `ParticipantAnchor`, `IncidenceGraph`.

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
incidence nodes. The parallel canonicalization level for this carrier is `Structure`.

**Cost, and why it is not the default.** Because nauty does not accept edge colours, every localized
bond must also become a pseudonode, and a molecule has far more bonds than overlays. So
`SubstructureMatchAlgorithm::GraphAndOverlays` — match the skeleton, then verify overlays against the
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
**In code:** `check_integrity`, `*IntegrityError`.

### Integrity check

An **integrity check** is the graph-IR-owned, error-valued operation that enforces tier 1. It returns
`Result<(), *IntegrityError>` and is shared by checked constructors, boundary conversions, and every
path that publishes the IR type. Trusted asserted constructors use the same implementation and
change only the failure reporting. There is no `*Checker` object.

**Not:** a validator. Validators return semantic `Solution` values and belong in `umol-graph`.
**In code:** `Molecule::check_integrity`, `Reaction::check_integrity`,
`ReactionSpan::check_integrity`; the corresponding `*IntegrityError` types.

### Invariant

An **invariant** is a model-independent physical or mathematical condition over an otherwise
well-formed structure. Agreement between an independently meaningful stored constraint and its
model-independent derived value is also an invariant.

**Not:** integrity or conformance. Invariants are tier-2 semantic validation in `umol-graph`.
**In code:** `validate_invariants`, `*InvariantsValidator`, including
`ValenceInvariantsValidator` and `SpinInvariantsValidator`.

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

The aggregate canonicalization levels are `Topology`, `Constitution`, `Structure`, and `Full`.
`Topology` contains atoms and localized bonds. `Constitution` adds dative bonds, aromatic systems,
multicenter bonds, and noncovalent bonds. `Structure` adds stereo atoms and stereo bonds while
excluding constraints. `Full` adds normalized constraints through post-hoc selection among tied
structural frames and is identical to unqualified canonicalization. Para-stereo is context-dependent
refinement within the structure pass shared by `Structure` and `Full`, not another level. Future
structural entity kinds extend the first applicable structural level without changing the meanings
of the earlier ones; future constraint variants extend `Full` append-only.

`IncidenceLevel` stops at `Full`, using that name for the complete structural carrier because
constraints do not occur in an incidence graph. Thus `IncidenceLevel::Full` supplies the carrier for
`CanonicalizationLevel::Structure` and for the structural phase of `CanonicalizationLevel::Full`.

The corresponding entity domains are topology (AB), non-stereo (DAMN), and stereo (SS).
Constitution is topology plus non-stereo; overlays are non-stereo plus stereo. Domains are
compositional groups, whereas levels are the nested public selectors.

**Not:** *features*, which are independently combinable bitflags; *selection*, which does not express
the nested relation.
**In code:** `IncidenceLevel`, `CanonicalizationLevel`.

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

A **model** contains semantic choices defining which result is chemically accepted.

**Not:** config, policy, algorithm.
**In code:** `ChemistryModel`, `ValenceModel`, `AromaticityModel`, `StereoModel`, `RingModel`.

### Molecular structure

A **molecular structure** is a molecular topology together with its overlay entities.

**Not:** molecular topology, which excludes them.
**In code:** —

### Molecular topology

The **molecular topology** is the attributed undirected simple graph of atoms and localized bonds,
carrying no aromatic, stereo, or coordination information.

**Not:** molecular structure, which includes the overlay entities.
**In code:** —

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

**Normalized equivalence** is `equiv`: two values are equivalent in the current frame when their
normal forms are structurally equal. `Normalized<T>` carries the guarantee that normalization has
already run, so its derived `Eq`, `Hash`, and `Ord` operate on normal forms and can be used for
semantic deduplication.

**Not:** aggregate canonicalization, which selects an entity and participant frame and requires an
explicit context. Not chemical standardization, resolution, validation, or repair.
**In code:** `Normalize`, `normalize`, `normalized`, `Normalized<T>`, `Equiv::equiv`.

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
**In code:** `ParticipantPosition`, `RelationParticipant`, `ParticipantAnchor`.

### Patch algebra

A **delta** is the morphism between two entity states, and the pair of operations over it is a patch
algebra: `apply` carries a state forward by a delta, `diff` factors two states back into the deltas
between them.

The law is `apply(lhs, diff(lhs, rhs)) == rhs`. The entity update API states the same law as
`x.update(x.difference_to(y)) == y`.

**Naming split to resolve.** One law, two spellings — `apply`/`diff` on `EntityPatch`,
`update`/`difference_to` on the entity update surface. Pick one pair before either grows further; a
reader who learns the law under one name will not find it under the other.

**Not:** transaction application, which executes a whole edit plan and publishes only on success. A
patch is an entity-level morphism; a transaction is a molecule-level lifecycle.
**In code:** `EntityPatch::apply`, `EntityPatch::diff`, `update`, `difference_to`.

### Perception

**Perception** applies a selected chemistry model and algorithms to identify or assess structural
entities. It produces a policy-free derivation.

**Not:** resolution (which applies policy and edits), validation.
**In code:** `AromaticityPerception`.

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
entity persists in `K` and its label is resolved per side.

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

### Relation set

A **relation set** stores n-ary relations over typed participants — `NodeId`, `EdgeId`, or an
external type implementing `RelationParticipant` — with a shared incidence index so that a relation
is reachable from any of its participants.

Three axes parameterize a set, and the vocabulary is worth keeping straight:

- **arity** — `FixedRelationSet` has compile-time-known arity; `VarRelationSet` is variable-arity;
- **ordering** — `Unordered` or `Ordered`, which is what controls canonicalization;
- **factor** — a birelation set relates *two* factors, each with its own participant type, ordering,
  and arity. `FixedFixedBirelationSet`, `FixedVarBirelationSet`, `VarVarBirelationSet`.

The incidence index of a birelation set spans both factors, so a relation is reachable from a
participant regardless of which id-space it belongs to.

**Not:** a relational-database table, though the whitepaper uses that reading to explain overlays.
Not the overlay itself: an overlay entity is *stored in* a relation set.
**In code:** `FixedRelationSet`, `VarRelationSet`, `FixedFixedBirelationSet`, `ParticipantAnchor`,
`RelationParticipant`.

### Relational constraint

A **relational constraint** is a molecule-scope, reference-bearing constraint relating an overlay
entity to atoms, bonds, roles, or predicates.

**Not:** a synonym for *incidence constraint*, merely because the derived value depends on a
relation.
**In code:** `RelationalConstraint`.

### Remapping

A **remapping** is a total old-to-new relabeling: every source id has an image and no entity is
dropped. Totality is directional. A remapping may inject a source into a larger ambient id space,
such as an lhs-anchored reaction union, so it is not necessarily surjective, bijective, or
reversible.

`umol_graph_core::Remapping` transports graph nodes, edges, and relation participants.
`IdRemapping` transports typed references across all eight molecule entity families. When relation
participant order changes during transport, the relation-set operation also permutes positional
relation data so that values remain aligned with their participants.

A correspondence that is total on the left can produce a remapping. Applying a remapping to a
complete standalone molecule requires a bijection onto dense target tables, which corresponds to
totality on both sides. In that case remapping is semantics-preserving alpha-renaming rather than a
structural edit.

**Not:** a correspondence, which may be partial and only records pairing; not a compaction, which
expresses removal by leaving removed source ids without images.
**In code:** `umol_graph_core::Remapping`, `IdRemapping`, `remap`, `to_remapping`.

### Reset

**Reset** clears a constraint by setting it to its undetermined form.

**Not:** a general synonym for removing an entity or replacing determined structural information.
**In code:** `reset_aromatic_valence`, `reset_stereo_constraints`.

### Resolution

**Resolution** applies configured recovery policy, constructs one atomic edit plan, and may apply
that plan transactionally.

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
single-value operations, such as `automorphisms`, are not delivery prefixes. `try_*` and
`*_fallback` mark input-domain dispatch (simple graph against subdivision fallback) and compose
with the delivery prefix. The execution contract behind the prefixes — streamability, visitor
payload, and emission order — is defined in the algorithm execution guide.

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

### Structure

**Structure** is the complete structural canonicalization level: topology plus every overlay,
AB + DAMN + SS. It includes all inherent fields and typed incidences of every entity kind while
excluding constraints. Para-stereo changes refinement within this level rather than defining
another level.

**Not:** `Full`, which additionally uses normalized constraints to select among tied structure
frames; molecular topology, which excludes overlays.
**In code:** `CanonicalizationLevel::Structure`.

### Topology

**Topology** is the AB entity domain and the lowest structural level: atoms and localized bonds,
including the inherent values carried by those entities. It does not include overlay entities.

**Not:** constitution, which adds the non-stereo overlays; the incidence graph, which is an
algorithmic representation constructed from selected structure; constraints.
**In code:** the `Topology` variants of `IncidenceLevel` and `CanonicalizationLevel`.

### Transaction

A **transaction** is the journal of realized undos for one batch of edits. `transact` applies a batch
and returns it; `rollback` consumes it to restore the prior state; `append` concatenates two so a
multi-stage operation can be reversed as a unit.

`transact` is the checked path. `transact_unchecked` applies without producing a journal and cannot
be rolled back.

**Not:** *application*, which executes a plan and publishes only on success; a transaction is the
mechanism application uses. Not *patch algebra*, which is the entity-level `apply`/`diff` pair.
**In code:** `Transaction`, `transact`, `transact_unchecked`, `rollback`, `append`,
`TransactionError`.

### Transformation

A **transformation** explicitly rewrites one valid representation into another. Kekulization,
aromatization, and charge delocalization are transformations because they alter determined
representation rather than fill undetermined state.

**Not:** a resolver policy; not resolution.
**In code:** `umol-graph/src/ops/transform`.

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
entry carries three slots:

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
