# 194 — Constraint channel: assertions only, projections read on demand

Status: In Progress
Date: 2026-08-10
Relates: [053](053-molecule-validation-scheme-2026-02-17.md),
[125](125-constraints-as-projections-2026-06-22.md),
[128](128-substructure-derived-predicates-2026-06-23.md),
[138](138-constraint-container-api-2026-07-06.md),
[149](149-molecule-ring-cache-and-hashing-2026-07-13.md),
[174](174-aromatic-hydrogen-resolution-2026-07-31.md),
[193](193-subpattern-constraints-2026-08-09.md)

## Problem

Stored constraints go stale under mutation: editing topology leaves derived constraints
(`#v`, `#a`, …) on entities untouched, so a mutated molecule can carry constraint values that
contradict its own relations. Doc 125 diagnosed the cause — projections of the relational structure
are denormalized into the entity rows — and this document settles the remedy. The whitepaper's
Mutation section currently has no honest account of what happens to constraints under edits, which
blocks publication.

Doc 125 identified three roles muxed into one channel: projections of relations, query predicates,
and un-normalized staging input. The nomenclature guide already fixes the semantics the channel is
supposed to have: a constraint is *a possibly non-ground assertion*, and *projection* names a
mapping, not a stored object. Two code paths violate that definition by writing projections into
the assertion channel:

- atom-typing valence resolution extends the stored container with derived constraints before
  candidate selection and writes the narrowed result back
  (`umol-graph/src/ops/valence/atom_typing.rs`);
- substructure matching materializes per-atom clones of the host with derived constraints extended
  in (`host_match_targets` in `umol-graph-ir/src/ir/substructure.rs`).

Everything else is compensation: the `zeroed()` omission threshold hides materialized redundancy on
lowering, hand-rolled constraint removals in the kekulizer, aromatizer, and aromaticity resolver
undo materialization locally, and `ConstraintValidator` detects — but cannot prevent — conflicts
after mutation.

## Model

In database terms:

| umol concept                        | database analogue               | staleness behavior                          |
| ----------------------------------- | ------------------------------- | ------------------------------------------- |
| topology, overlays, inherent fields | base relations                  | primary data; mutation operates here        |
| projection (valence, degree, rings) | virtual generated column (view) | computed on read; cannot go stale           |
| stored constraint                   | CHECK-style assertion           | may become unsatisfiable; a validity finding |
| pattern constraint                  | WHERE predicate                 | queries are not updated by data changes     |
| `ConstraintValidator`               | deferred CHECK evaluation       | unchanged                                   |

The rule: **the stored constraint channel holds only assertions; projections are never stored.**

The pattern/concrete asymmetry dissolves under this rule. The distinction is per key,
determined versus undetermined: an assertion is redundant exactly when the relational structure
determines that key and the derived value refines the assertion. A pattern's topology never determines
its constraint keys, so pattern assertions are primary data — the same rule, not an exception. A
molecule whose stored assertion conflicts with its structure is carrying an unsatisfiable
assertion; validation reports it, consistent with edits enforcing no chemistry.

Constraint placement (settled 2026-08-12): an entity's assertion home is its inline container;
the molecule-scope list carries `Relational` and `Molecule` leaves, combinators, and — as input —
bare entity leaves, which remain a valid DSL spelling (the spec's lift/inline section stands).
Invalidating the bare spelling was considered and dropped: a singleton `:or`/`:and` over an
entity leaf denotes exactly the bare assertion, so a validity rule distinguishing them would
legislate spelling, not content. Placement is instead normalized by resolution: the pipeline
opens with a placement stage that applies the context-free constraint normalization, reduces
trivial wrappers (a singleton `:or`/`:and` is its element), and moves bare entity leaves inline —
`inline_constraints` gains its production caller, and the bridges stay. After the stage,
`asserted(key)` is complete for entity-addressed content; assertions inside irreducible
disjunctions remain molecule-scope logic, evaluated where the list is evaluated (validation,
discharge, doc 195 matching). Patterns are never resolved, so both spellings persist there and
pattern evaluation handles them (doc 195).

Canonicalization is unaffected by placement (settled 2026-08-12). Placement normalization is
fallible — colliding assertions combine by meet, and `⊥` is `Contradictory` — while
canonicalization preserves the complete represented assertion and cannot fail; placement therefore
belongs to resolution and to no canonicalization level. Full remains full, the default level
stands, and differently-placed descriptions stay canonically distinct, exactly as a Kekulé input
and its perceived form do; the equivalence that merges placements is canonical equality after
resolution. Trivial-wrapper reduction (a singleton `:or`/`:and` is its element) is by contrast
infallible and content-preserving and belongs in the context-free `Normalize`, completing the
conjunction-flattening family; Full canonical equality thereby strengthens for those degenerate
spellings alone, changing their canonical keys across versions.

## Three flows, three fates

| flow                        | example                                             | fate                                                            |
| --------------------------- | --------------------------------------------------- | --------------------------------------------------------------- |
| projection materialization  | `derive_constraints` extended into storage          | eliminated; replaced by transient per-key reads                 |
| model-supplied narrowing    | registry rows surviving for one atom                | candidate sets in solver state; never stored; commit lands primary data |
| input assertions            | raised SMILES `#a` negatives, user-supplied `#v4`   | stored; removed by discharge once the structure determines them |

A single registry row meeting into an atom would be an assertion *by the chemistry model*, but the
valence phase does not produce a single row: it produces the set of rows that survive narrowing,
and doc 174 shows that collapsing the set locally is the defect — `compare_valence_preference`
cannot know which member lets the ring satisfy 4n+2; aromaticity never votes. The constraint
channel cannot carry the set: per-key assertions are conjunctive, so a container denotes the
*product* of its per-key value sets, while a candidate set is a disjunction of correlated tuples.
Pyrrole-type nitrogen pairs `#h1` with `#a2` and pyridine-type pairs `#h0` with `#a1`; the per-key
join `#h{0,1} #a{1,2}` admits the non-candidate `#h1 #a1`. Any single stored value is therefore
premature collapse or correlation loss. The candidate set is solver state by type — the doc 053
phase shape `Set<CandidateState> → Set<CandidateState>`: the valence phase emits the surviving
set, the aromaticity phase selects with the criterion it alone has, and only the committed outcome
reaches the molecule, as primary data. The stored channel thus holds only input assertions, and
discharge — not cache invalidation — removes them once the structure determines them; nothing
derivable-from-structure is ever stored.

## Read-side restructure

The `Lattice` trait gains one provided method and is otherwise untouched: `satisfies`, the
receiver-inverted reading of `matches` — `fn satisfies(&self, pattern: &Self) -> bool
{ pattern.matches(self) }`, a universal default never overridden. `matches` keeps its
SMARTS-anchored meaning and pattern-as-receiver direction; `satisfies` exists because a view must
be a receiver and the view stands on the target side.

What restructures are the container-level comparison entry points (`matches`, `is_compatible`,
match-target construction), which currently demand two materialized containers. They move onto
per-entity constraint views reached by accessor chaining — `atom(i).constraints()` — exposing the
two sides of every key:

- `asserted(key)` — the stored side; open-world by definition, absence is the vacuous constraint;
- `derived(key)` — from present relations only; vacuous on absence;
- `derived_complete(key)` — the closure reading: under the caller's claim that the relation set
  is complete, absence of a resolution-written table closes to the definite negatives.

The view is always the receiver of a comparison, never a function argument; ring context attaches
as data (`with_rings(&RingSet)`, owned data per the views rule). Comparisons fix the reading
their semantics implies: `satisfies(&pattern)` evaluates the host under `derived_complete` (query
against host), `is_compatible(&other)` under `derived` (narrowing admissibility). There is no
`effective(key)` accessor: no consumer needs the meet of the sides as a value — comparisons meet
internally, and discharge and validation read the sides and apply value-level `meet` directly,
where `None` already means `⊥` and folds to `false` or `Contradictory`.

- Evaluation is driven by the pattern's keys, so only requested derivations are computed. The
  per-key selector of doc 125 falls out, and the monolithic `derive_constraints` retires.
- A bare entity form needs no view: its asserted side is its container, read directly (the
  electron-count invariant check already does), and it has no relations to derive from.
- Ring context is built once per run by the consumer that scans keys, as `ConstraintEvaluation`
  already memoizes for validation.

This restructure is unconditional, not an alternative to consumption: matching takes the host by
shared reference and cannot write-then-remove, and candidate selection needs the combined reading
regardless of where narrowing is stored.

Consumers that change:

- `host_match_targets`: per-key `satisfies` evaluation replaces the `Cow` clone-and-extend path
  (also the doc 122 direction: do not materialize the host).
- atom-typing candidate admission: `is_compatible` against the view; the phase emits the
  surviving candidate set instead of writing narrowed constraints back — derived values
  participate in admission but are never written.
- the electron-count invariant check: already reads the container — the asserted side of a bare
  form — directly; no change.
- `ConstraintValidator`: already per-key derived-versus-asserted; may reuse the keyed accessors.

## Discharge

Only input assertions are ever stored, so removal concerns them alone, and it has one home: the
resolve pipeline ends with a *discharge* pass, after which no stored assertion is
determined-redundant. Determination needs no new predicate — it is `derived_complete(key)`
yielding a ground value (`is_ground`), existing vocabulary composed. A ground derived value that
refines the assertion ⇒ redundant ⇒ removed; the meet is `⊥` ⇒ contradiction; a non-ground
derived value with a strictly narrower assertion ⇒ kept. An operation that
realizes an assertion may remove it early in its own transaction — the aromaticity resolver does
this for `#a` today — but that is an implementation liberty inside the pipeline, not a second
concept; the guarantee is the closing pass. Distributed keys are why the pass must exist: `#v` is
realized by no single operation, only by all incident bond orders and the implicit hydrogen count
becoming literal, and nothing creates rings, so `#R` is determined by topology alone.

The kekulizer and aromatizer removals are a different contract: a transform that deletes a
relation must delete the assertions its output no longer satisfies, which is the existing
emit-compliance policy and stays per-operation.

When `derived_complete` is ground, per key (to be confirmed row by row during implementation
planning):

| key                          | determined when                                                    |
| ---------------------------- | ------------------------------------------------------------------ |
| `#v` valence                 | all incident bond orders literal                                   |
| `#D` degree                  | always (topology)                                                  |
| `#X`, `#H`, `#V` totals      | implicit hydrogens literal (+ explicit H neighbors' elements literal) |
| `#d`, `#t` dative pairs      | dative overlay incidence under the closure                         |
| `#a` aromatic valence        | when the aromatic system is created (removable early)              |
| `#m` multicenter valence     | when the multicenter bond is created (removable early)             |
| `#T`, `#C` stereo            | when the stereo overlay is created (design 103 flow)               |
| `#R`, `#x`, `#y` ring keys   | always (pure function of topology); checked by discharge, gated on presence |
| overlay-entity constraints   | at the resolution stage of the owning overlay                      |

## Resolution carrier

Settled 2026-08-12. The valence phase's output for an atom is a set of *completions* — one per
registry row surviving `is_compatible` admission against the atom's constraints view:

- A completion is a ground `AtomForm` — the unified form in its ground state, per the doc 080
  commitment that ground is a state of the one type, not a separate type. No primitive-valued
  completion struct exists: it would reintroduce the retired ground/pattern split and decouple
  the coupled spin fields. Element and charge are fields like any other — today every disjunct
  carries the same literals because admission is keyed on them; a strategy enumerating charge
  states varies them with no type change.
- The disjunction is extensional: `SmallVec<[AtomForm; 1]>` per atom, nearly always a singleton.
  Each emitted disjunct is ground in every key it carries — points, never rectangles; a single
  non-ground form would denote the Cartesian product of its fields and lose coupling. This is the
  emit contract of the phase; the type stays the permissive unified form.
- The carrier maps `AtomId` to its disjunction (`AtomCompletions`). "Completion" is a role a
  ground `AtomForm` plays, as "pattern" is a role a `Molecule` plays. Admission emits disjuncts
  by meet; commit applies the chosen form via `difference_to` — existing machinery. A candidate
  set renders as atom strings for inspection.
- Openness (recorded 2026-08-12): extensional members and atom-to-atom decoupling are the *emit
  contract of the current valence phase*, not commitments of the carrier types, and must not be
  written into their rustdoc. A future relaxation takes the shape of constraint equations over
  the existing value-expression variables — intensional members are already type-inhabitants
  (non-ground forms with variables), and inter-atom coupling is an additive molecule-scope
  component behind the carrier's private storage. The equations build on the variable facility,
  not beside it — the parallel-mechanism collision is what killed the earlier joint-domain
  constraint. Consumers read members through the `AsLit` extraction discipline (explicit
  non-literal handling), so intensional members degrade to defined behavior, not panics.

Phase composition keeps early commit: the valence phase applies edits for singleton atoms as today
and forwards only plural sets; the aromaticity phase selects jointly per candidate aromatic system
— enumerating assignments over the product of the member atoms' candidate sets and keeping those
whose π-electron sum satisfies the model — and emits its edits including the chosen completions.
Survivor count zero is `Contradictory`; one is `Determined`; more than one falls to
`compare_valence_preference` as a tie-break of last resort whose use is visible in the verdict,
else `Underdetermined`. Single commit was considered and rejected: every failure path — error,
contradiction, mid-pipeline underdetermined — already rolls back the whole journal, and `resolve`
holds the molecule exclusively, so intermediate states are unobservable either way.

The `Underdetermined` payload of `Resolver::resolve` becomes the carrier itself — per-atom
surviving completions, not a cardinality summary — so a caller can inspect and refine: assert what
disambiguates through ordinary edits and re-resolve. The normal path subsumes a select-and-commit
interface, so none is added; `Solution` itself is unchanged. The SMILES reader stops pinning
`#h0` on bare aromatic heteroatoms (the doc 174 prerequisite). Acceptance: the doc 174 regression
triple and the pyrrolyl DSL case; the corrected F420 structure joins once the zero-contributor
perception and registry-coverage items of doc 174 land.

## The closure is per call

Whether absence of an overlay reads as a definite negative (`NotAromatic`) or as no evidence is a
property of the consuming operation, not of the molecule:

- a stored molecule-level status flag is itself a materialized projection of "is everything
  determined?" and would go stale under mutation — the same defect one level up;
- "perception ran" is history, not structure; it is not recoverable from the data, but every caller
  knows its own context (matching presents hosts as fully perceived; resolution knows it is
  staging);
- a two-type split (ground versus staging molecule) is the type-level version but is heavy and
  Python-visible, and the boundary is per key, not per molecule.

The choice therefore lives in the accessor pair, not in a mode parameter: `derived` reads present
relations only, `derived_complete` adds the closure. The closure's license tracks exactly the
tables resolution writes: phases create overlays and never atoms or bonds, so the topology keys
(`#v`, `#D`, `#X`, `#H`, `#V`) read identically under both accessors, and the overlay-incidence
keys (`#a`, `#m`, `#d`/`#t`, `#T`) differ in the absence cell alone. The un-closed reading is not
a resolution internal: patterns are permanently partial public descriptions, and validating one is
coherent only under `derived` — the closure would turn every overlay assertion into `⊥` against
overlays a pattern never stores. An earlier design carried a `RelationCoverage`
(`Complete`/`Partial`) parameter on view construction; it was replaced because it encoded a 2 × 2
scheme with an uninhabitable cell (asserted constraints are partial by definition) and a mode that
every caller fixes statically — the three-accessor surface states the actual structure.

## Rejected alternatives

- **Refresh-on-write** (stored generated column; triggers in `transact`): the mutation layer would
  need every derivation rule including ring perception — inverted layering — and would have to be
  suppressed for patterns, dragging a world-assumption flag into `transact`.
- **Provenance tag** (asserted/derived bit per constraint): constraint containers are open data
  carriers; an operation-issued marker inside one must be maintained by every operation and
  distrusted by every consumer. It systematizes the trigger-maintained column instead of
  eliminating it.
- **Effective-value reads over retained storage** (no removal): fixes readers but keeps redundant
  storage, keeps the `zeroed()` elision, and makes equality of resolved molecules depend on
  materialization history.
- **Transient cache**: addresses cost, and this is a semantics problem. Caching remains available
  later as a pure optimization over immutable snapshots (the pattern doc 149 removed), fully
  decoupled from constraint semantics.
- **Stored staging narrowing** (write model-supplied narrowing into the constraint channel during
  resolution, remove when determined): a per-key container can only hold the collapse of a
  candidate set, so the doc 174 premature-collapse defect is structural under this option, not a
  preference-tuning problem. It also leaves mid-pipeline molecules carrying chemistry-model
  opinions in storage and produces create-then-remove churn in transactions. A candidate-set
  carrier between phases is unavoidable for doc 174 regardless; storing narrowing as well would
  add a third channel next to it.

## Consequences

- Ground molecules obtained through resolution store no constraints; staleness under mutation is
  impossible by construction. Equality of resolved molecules no longer depends on materialization
  history.
- Lowering: nothing to elide on resolved molecules; the constraint-elision role of `zeroed()` is
  retired. Conformance snapshots change wholesale (the 652-snapshot experiment of doc 125, in
  reverse: outputs carry only primary data).
- Raising is unchanged: input dialect assertions land as staging and are discharged by resolution.
- Whitepaper: the Mutation section states that edits operate on relations and inherent fields;
  projections are computed, so they cannot go stale; stored assertions may become unsatisfiable,
  which validation reports; resolution discharges assertions the structure comes to determine.

## Staged implementation plan

Modules, bottom-up: `umol-graph-ir` (constraint views, substructure), `umol-graph`
(valence, aromaticity, resolve, invariant), `umol-io` (SMILES ingest), `umol-py`, conformance.

### S0 — constraints-view foundation (`umol-graph-ir`); green throughout

- S0a: `Lattice::satisfies` as a provided method (universal default: `pattern.matches(self)`);
  the dual recorded in the nomenclature Matching entry. Additive. **Done 2026-08-12:** duality law
  added to `assert_lattice_laws`, covering every `Lattice` impl in the property sweep.
- S0b: extract the per-quantity derivation functions from the `AtomView`/`RingAtomView` method
  bodies (typed accessors become delegates; signatures and their tests unchanged); retype
  `AtomView::constraints()` to `AtomConstraintsView` — the container's read API inherited with
  meanings intact, plus `asserted`/`derived`/`derived_complete` and `with_rings`. Near-additive:
  the accessor was read-only, so the retype is almost source-compatible. **Done 2026-08-13:** two
  call-site edits workspace-wide; suite and clippy green.
- S0c: `satisfies`/`is_compatible` on `AtomConstraintsView`. Additive. [dep: S0b] **Done
  2026-08-13:** pattern-key-driven, sides met internally, `⊥` folds to `false`; conflicted-host
  and absent-overlay-skip semantics pinned by cases.
- S0d: the same retype and keyed surface for `BondView::constraints()` (`BondConstraintsView`
  with `BondConstraintKey`). Additive. [dep: S0b] **Done 2026-08-13:** bond `#a` becomes derivable
  for the first time (both-endpoints-in-system incidence, closure negative `Aromatic(false)`);
  `derive_constraints` unchanged in behavior (still emits only `#C`), reimplemented over the
  dispatch; two validator call sites migrated.
- S0e: the same for the remaining six entity views (uniform surface). Additive. [dep: S0b] —
  gates only S6a. **Done 2026-08-13,** mirroring the validator's established derivations: dative
  aromatic incidence is binary-only (doc 117 stub; ring key unprojected); system/multicenter
  electron counts are self-projections with no absence cell; noncovalent intramolecularity is
  shared localized-bond component membership; the stereo constraint kinds have no projection —
  their derived side is vacuous under both modes. The six views are one `constraints_view!`
  macro plus per-family typed getters.
- S0f: nomenclature entries: view family, the `derived`/`derived_complete` readings, discharge,
  the `satisfies` dual in the Matching entry; reword the two `invariant.rs` "standalone" doc
  comments to the molecule-atom convention. Additive. **Done 2026-08-13:** glossary gains
  "Constraints view", "Derived and asserted", and "Discharge" (the `satisfies` dual landed with
  S0a); the two comments now say "a bare `AtomForm`".
- S0g: facade-rule enforcement from the doc 165 ring-view worklist: delete the view-as-argument
  relations `RingView::shared_atoms`/`shared_bonds` (`view/ring.rs:101-105`) — the id-keyed
  `RingSet::shared_atoms`/`shared_bonds` already exist and are the only form retained; the two
  `RingView` unit tests retarget to the `RingSet` forms. No production callers; green. The
  remaining doc 165 ring-view items (id/view accessor conventions) stay in doc 165. **Done
  2026-08-13:** the id-keyed `RingSet` forms already had their own tests, so the two `RingView`
  tests were redundant and deleted rather than retargeted; `intersection` import dropped.
- S0h: extend the context-free `Normalize` for constraint trees with trivial-wrapper reduction
  (a singleton `:or`/`:and` is its element), completing the conjunction-flattening family. Full
  canonical equality strengthens for those degenerate spellings alone (their canonical keys
  change); placement remains outside canonicalization by the fallibility criterion (Model
  section). Normal-form change; green. [dep: none] **Done 2026-08-13:** two existing empty-child
  unit expectations updated per this spec, six new cases (including dedup-to-singleton and
  cross-combinator reduction); new property module `property/constraint.rs` states the tree's
  normalize laws directly — idempotence and wrapper-free normal form — over a raw tree domain
  whose combinator arity 1..=3 generates the degenerate spellings; zero conformance churn in the
  full-suite stage-end run.

The S0 surface (final form, settled 2026-08-13 after two revisions: a standalone constraints
view was built, dissolved into `AtomView` when its accessor proved unnameable, then reinstated as
a *chained* facade when the dissolution proved to trade the naming problem for a scope problem —
an atom-scoped receiver with constraint-scoped comparisons):

```rust
molecule.atom(i).valence()            // AtomView typed quantities — unchanged
molecule.atom(i).constraints()        // -> AtomConstraintsView<'a> (retyped accessor)
    .valence(), .iter(), .is_empty()  // inherited container read API — asserted, meanings intact
    .asserted(key)                    // keyed core: the stored side
    .derived(key)                     // present relations only; vacuous on absence
    .derived_complete(key)            // the closure reading
    .with_rings(&ring_set)            // ring context; ring keys panic without it
    .satisfies(&pattern)              // S0c — derived_complete reading
    .is_compatible(&other)            // S0c — derived reading
```

- Scope lives in the receiver, not in name prefixes: `atom(i).valence()` is the derived quantity
  and `atom(i).constraints().valence()` the asserted payload — the pre-existing pattern,
  preserved verbatim; the view slots under the established names. Accessor chains returning
  narrower views are the namespacing mechanism (`RingViews::atom → RingAtomView` is the
  precedent); a facade's *implementation* still never builds another facade — the logic lives in
  the derivation functions beneath (`view/atom.rs`, `view/ring.rs`), with `derived_constraint` as
  the keyed dispatch and `derive_constraints` reimplemented over it until S4g deletes it.
- `AtomConstraintsView { molecule, atom, rings }` — no context enum, no form constructor, no
  `Molecule` accessor. The bare-form case needs no view (its asserted side is its container,
  read directly — stage S2 dissolved), and the only path is the chain.
- No typed derived accessors (`derived_valence()`, …): the keyed core is the complete surface;
  typed sugar returns if usage shows the need. `.valence(coverage)` on `AtomView` was examined
  and rejected — the parameter is a no-op for topology-derived quantities and lossy for overlay
  quantities (`NumForm` cannot carry `NotAromatic` versus `Aromatic(0)`), the same exit that
  killed `RelationCoverage`.
- No iteration on the derived side — `RingMembership(RingScope)` is open-ended, so evaluation is
  driven by the consumer's key set; asserted iteration is the inherited `iter()`.
- The retype was near-source-compatible: `constraints()` always returned an immutable reference,
  so no write site existed, and read sites compile through the inherited API. Two call sites
  needed edits (a validator `get` → `asserted`; ring validation taking the container as a value
  via `attributes.constraints`).
- Inherent `satisfies`/`is_compatible` must never land on a type that impls `Lattice` (forms,
  containers): inherent methods silently shadow trait methods. `AtomConstraintsView` impls no
  `Lattice`.
- No `effective` accessor and no `determined` predicate: consumers meet the sides with
  value-level `meet` (`None` = `⊥`, folding to `false` in comparisons and `Contradictory` in
  discharge), and determination is `derived_complete(key)` ground (`is_ground`) — existing
  vocabulary composed.
- Naming record: `effective(key)` designed then deleted (no consumer needs the meet as a value);
  `meet(key)` rejected (uniformly binary over like values); `derived_partial` rejected (the
  closure is the marked step, not the partiality); `projected()`/`projected_complete()`
  considered and dropped for the paper's side vocabulary (`asserted`/`derived`); the standalone
  accessor names (`atom_constraints_view`, `atom_constraints`) rejected (convention break versus
  container-word collision); the full dissolution rejected in turn (hidden scope change); the
  chained facade is the resolution. "Skeleton keys" retired for *topology keys*, aligning with
  the molecular-topology glossary term.

### S1 — matching on views (`umol-graph-ir::ir::substructure`)

- S1a: `host_match_targets` and the predicate closures evaluate per pattern key via the views'
  `satisfies` (the `derived_complete` reading); ring context built once per run as today; the
  `Cow` target materialization is deleted. The pattern scan gains a gate: a non-empty
  molecule-scope `Constraints` list on the pattern is an error naming the construct — matching
  becomes fallible for exactly that input class, sound-but-incomplete until doc 195 replaces the
  gate with evaluation (the current silent ignore admits false positives). Breaking, green at
  stage end: substructure, fingerprint, and reaction-matching suites unchanged. A host carrying
  unconsumed input assertions now meets them against the closure; any test relying on such a
  staging host is surfaced, not silently rewritten. [dep: S0c, S0d] **Done 2026-08-13:**
  `SubstructureMatchError::MoleculeScopeConstraints`; both public entry points return `Result`;
  propagation through reaction application adds `ApplyPreconditionError::Match(#[from] ...)`;
  `host_match_targets` replaced by `host_ring_context` + per-family satisfies tables; field
  matchers destructure exhaustively; matcher table gains rows per derived key (`#v #D #X #H #V`,
  `#a` positive/negative/Kekulé-flag, `#m`, `#d #t`, `#T`, bond `#a` in-system and
  asserted-vs-closure `⊥`, bond `#C`). The gate surfaced six `MissingEntry`-via-matcher tests
  (constraint-remove LHS carries the removed molecule-scope constraint); they are `#[ignore]`d
  pending doc 195 evaluation, which lists them as an explicit unskip item (2026-08-13).

### S2 — dissolved

- S2a is no work: a bare form's asserted side is its container, which the electron-count
  invariant check already reads directly. The stage label is retained so later references stay
  stable.

### S3 — completion carrier (`umol-graph::ops::valence`); additive, green

- S3a: `AtomCompletions`, the keyed carrier. Additive. **Done 2026-08-13:**
  `ops/valence/completion.rs`, re-exported from `ops::valence`; the settled surface exactly
  (private `BTreeMap` storage, `insert` asserting the non-empty invariant, slice-typed reads);
  rustdoc carries the invariant and deterministic order only, per the openness note.
- S3b: `ResolveReport`, the resolver verdict payload. Additive. [dep: S3a] **Done 2026-08-13:**
  same module as the carrier, re-exported from `ops::valence`; public fields per the settled
  surface, derives only — no methods and no tests to carry, since a descriptive record with no
  invariant has nothing but derived behavior.

The S3 surface (settled 2026-08-12):

```rust
pub struct AtomCompletions {
    entries: BTreeMap<AtomId, SmallVec<[AtomForm; 1]>>,   // deterministic order
}
// insert (asserts entry non-empty — the one representation invariant),
// get, remove, iter, len, is_empty

pub struct ResolveReport {
    pub unresolved: AtomCompletions,   // plural survivors; empty under Determined
    pub tie_breaks: Vec<AtomId>,       // preference-selected atoms; sorted, deduplicated
}
```

- No completion value type: a completion is a ground `AtomForm` (Resolution carrier section).
  `AtomFields`, `AtomCompletion`, and ground-literal valence enums were considered and rejected —
  the retired ground/pattern split, decoupled spin fields.
- `AtomCompletions` holds the one representation invariant (entries non-empty — an empty set is
  `Contradictory`, not an entry), hence private storage with accessors; `ResolveReport` is a
  descriptive record with no invariant, hence public fields. Both are operation-issued but open
  construction is harmless and undefended.
- Tie-breaks record use, not an audit trail: the selected completions are committed into the
  molecule; a richer per-atom record waits for a consumer.
- Python (S4f): both bound read-only per the doc 192 type roles.

### S4 — pipeline rework; the one red stage, green only at its end

- S4a: atom-typing: admission via the view's `is_compatible`; no stored-constraint
  extension, no narrowed write-back; singleton atoms produce edits, plural atoms produce
  completions. Breaking. [dep: S0c, S3b] **Done 2026-08-13:** admission = field compatibility
  plus the view's `is_compatible` over row-constraint keys; disjuncts are full meets (row
  constraints ride in solver state only — a committed singleton restores the atom's own
  constraint container, so neither derived nor registry constraints are written back);
  `plan → Solution<(Edits, AtomCompletions), _>`, `resolve → Solution<AtomCompletions, _>`
  with singleton edits committed under `Underdetermined` (early commit).
  `classify_molecule_atom` composes the closure reading from the view's keyed core
  (meet of asserted, `derived_complete`, and the row entry per key; ring keys asserted-only) —
  sharpened from the old derived-replaces-asserted pattern, and topology keys outside the old
  six-key derivation now participate; S6a's reading selector absorbs this composed site.
  `ValenceResolver::plan` interim-drops the completions payload until S4d threads the carrier.
  `derive_constraints` has no callers left workspace-wide. Red baseline captured: 517
  failures in 4 targets — umol-graph lib 81 (ingest/morgan/pattern/parse), conformance
  resolution 402, fingerprint featurizer 10, umol-py lib 24 — in two classes:
  plural-admission `Underdetermined` (closed-shell selection returns as the S4c/S4d
  tie-break) and conformance snapshots whose outputs are already correct but record the
  retired `#v` write-back (e.g. `atomic_ions_hg+2_dimer`) — expectation deltas awaiting the
  S5d regeneration, per the S4h scoping. atom-typing unit tests rewritten to the new
  surface (28/28 green).
- S4b: counts: plural `#h` candidates emitted through the same carrier; `CountsInput` readings
  aligned with the view. Breaking. [dep: S3b]
- S4c: aromaticity resolve: joint selection per candidate system over the carrier;
  `compare_valence_preference` demoted to a visible last-resort tie-break. Breaking. [dep: S4a]
- S4d: `Resolver::resolve`: thread the carrier valence → aromaticity; `Solution<(), _>` →
  `Solution<ResolveReport, _>`; rollback paths preserved. Breaking. [dep: S4c]
- S4e: SMILES ingest stops pinning `#h0` on bare aromatic heteroatoms (`umol-io`). Breaking.
  [dep: S4d]
- S4f: `umol-py`: bind the report for inspection; update the resolution verdict mapping.
  Breaking. [dep: S4d]
- S4g: retire `derive_constraints` and `include_missing` (last callers died in S1a/S4a).
  Breaking. [dep: S4a]
- S4h: acceptance: the doc 174 regression triple and the pyrrolyl DSL case; full suite green
  (`--all-features --tests`, clippy) **except conformance snapshots**, whose remaining
  failures are audited to be exactly the planned expectation-delta class — outputs correct,
  snapshots recording constraints that are no longer written back (atom-typing `#v` since S4a,
  counts since S4b, `#h0` pinning since S4e). There is nothing to fix in that class; the
  expectations changed as planned, and the one regeneration stays at S5d, after discharge
  (S5b) and elision retirement (S5c) change the outputs once more. [dep: S4a–S4g]
- S4i: audit the stereo phase for the same premature collapse (doc 174's remaining open item) —
  whether a local preference selects before a later criterion can vote. Read-and-report while the
  pipeline is open; any fix is its own proposal, not S4 scope. [dep: S4d]

### S5 — discharge and output cleanup; one conformance regeneration

- S5a: per-key determination checks on the view layer; the determination table is confirmed row
  by row here. Additive. [dep: S0b]
- S5b: the placement and discharge stages bookending `Resolver::resolve`, in the same journal.
  Opening placement stage: apply the context-free constraint normalization (including the S0h
  trivial-wrapper reduction) and move bare top-level entity leaves inline, colliding assertions
  combining by meet with `⊥` ⇒ `Contradictory` — `inline_constraints` gains its production
  caller with its collision policy corrected from last-wins to meet, and the DSL spec's
  lift/inline collision note updated accordingly. Closing discharge pass: `⊥` ⇒
  `Contradictory`; also evaluates the remaining molecule-scope list with the validator's
  machinery: decided true with ground inputs ⇒ implied ⇒ removed; decided false ⇒
  `Contradictory`; undecided ⇒ kept (patterns are never resolved, so pattern constraints
  persist). Breaking — resolved outputs lose stored assertions. [dep: S4d, S5a]
- S5c: lowering: retire the `zeroed()` elision-only paths; raise-side dialect filling stays.
  Breaking. [dep: S5b]
- S5d: regenerate conformance snapshots once; final green: `--all-features --tests`, clippy.
  [dep: S5b, S5c]

### S6 — deferrable

- S6a: `ConstraintValidator` internals on the entity views. After S5 the constraint pass is
  vacuous on resolved molecules, so the validator becomes a staging/pattern tool; the per-kind
  hand-built comparisons collapse to "violation ⇔ the per-key meet of `asserted` and the selected
  derived reading is `⊥`". `ConstraintValidateConfig` selects the reading (`derived` versus
  `derived_complete`; selector representation decided here) with current per-kind absence
  semantics read off `incidence.rs` and preserved as the default. [dep: S0e]
- S6b: registry row `"C #v4 #a0"` (grounds as `#h0 #n0 #v4 #a0`; the exocyclic-carbonyl aromatic
  carbon of 2-pyridone, uracil, 4-pyranone, tropone, and tropolone) plus their snapshots. Pure
  data; independent — may land as early as S4h.
- S6c: correct the doc 174 "zero-contributor" diagnosis and pin it with tests. The perception
  refusal is the model's element scope, not zero handling: `AromaticityModel::daylight()` (the
  `ChemistryModel` default) excludes boron, so `is_atom_eligible` fails on scope
  (`hueckel_rule.rs:117`), `filter_ring`'s all-members rule drops the ring, and `derive`'s
  membership sweep (`aromaticity.rs:240-248`) reports `AromaticValenceFailure` for every atom
  asserting `Aromatic(_)` — doc 174's witnesses (borepin, borazine, 1,2-azaborine) are all boron
  cases. Zero contributions are handled correctly throughout: `Aromatic(Lit(0))` → `Some(0)`
  eligibility, plain summation, membership counted. Work: regression tests — the tropone family
  under `daylight()` (needs only S6b) and the borepin family under a boron-including scope
  (`permissive()`); a dated correction note in doc 174. No production logic change expected;
  scope membership stays a model decision that `ElementScope` already expresses.
- S6d: the F420 acceptance case (unique completion, C29H36N5O18P); its C/N ring is within the
  Daylight scope, so S6c is not a dependency. [dep: S4, S6b]

Critical path: S0 → S1 → S3 → S4 → S5. S1 must precede S5 because discharge strips assertions
from resolved hosts, after which matching must project every pattern key. S2 is dissolved; S3
is parallel to S1. The core deliverable — staleness-free semantics and the whitepaper Mutation
story — completes at S5; S6 is not required for it.

## Open items

- F420 enablement via the doc 174 zero-contributor and registry-coverage items (S6b).
- Matching of molecule-level pattern constraints: doc 195 (out of this document's critical path).
