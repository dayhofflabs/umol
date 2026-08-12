# 194 — Constraint channel: assertions only, projections read on demand

Status: Proposed
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
determines that key and the projection refines the assertion. A pattern's topology never determines
its constraint keys, so pattern assertions are primary data — the same rule, not an exception. A
molecule whose stored assertion conflicts with its structure is carrying an unsatisfiable
assertion; validation reports it, consistent with edits enforcing no chemistry.

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

The `Lattice` trait is untouched: `meet`, `join`, and `matches` on values are already
source-agnostic. What restructures are the container-level comparison entry points (`matches`,
`is_compatible`, match-target construction), which currently demand two materialized containers.
They move onto a per-entity constraints view — `AtomConstraintsView` and the per-entity family —
constructed per call from the assertions container, the entity's relational context, and a reading
mode, exposing the per-key *effective value*:

```
effective(key) = asserted(key) ∧ projected(key)    // meet; either side may be absent
```

The view is always the receiver of a comparison, never a function argument.

- Evaluation is driven by the pattern's keys, so only requested projections are computed. The
  per-key selector of doc 125 falls out, and the monolithic `derive_constraints` retires.
- A `⊥` effective value needs no special path: `matches`/`is_compatible` return false.
- A standalone entity form (no molecule context, e.g. the electron-count invariant check) is the
  same abstraction with an assertion-only view.
- Ring keys obtain their `RingViews` context inside the view, as `ConstraintEvaluation` already
  memoizes for validation.

This restructure is unconditional, not an alternative to consumption: matching takes the host by
shared reference and cannot write-then-remove, and candidate selection needs the combined reading
regardless of where narrowing is stored.

Consumers that change:

- `host_match_targets`: per-key source evaluation replaces the `Cow` clone-and-extend path (also
  the doc 122 direction: do not materialize the host).
- atom-typing candidate admission: compatibility against effective values; the phase emits the
  surviving candidate set instead of writing narrowed constraints back — derived values
  participate in admission but are never written.
- the electron-count invariant check: assertion-only source, semantics unchanged.
- `ConstraintValidator`: already per-key derived-versus-asserted; may reuse the view.

## Discharge

Only input assertions are ever stored, so removal concerns them alone, and it has one home: the
resolve pipeline ends with a *discharge* pass, after which no stored assertion is
determined-redundant. Per key, discharge asks whether the structure now determines the key: the
projection refines the assertion ⇒ redundant ⇒ removed; the meet is `⊥` ⇒ contradiction; the
assertion is strictly narrower than a not-yet-determined projection ⇒ kept. An operation that
realizes an assertion may remove it early in its own transaction — the aromaticity resolver does
this for `#a` today — but that is an implementation liberty inside the pipeline, not a second
concept; the guarantee is the closing pass. Distributed keys are why the pass must exist: `#v` is
realized by no single operation, only by all incident bond orders and the implicit hydrogen count
becoming literal, and nothing creates rings, so `#R` is determined by topology alone.

The kekulizer and aromatizer removals are a different contract: a transform that deletes a
relation must delete the assertions its output no longer satisfies, which is the existing
emit-compliance policy and stays per-operation.

Determination criteria per key (to be confirmed row by row during implementation planning):

| key                          | determined when                                                    |
| ---------------------------- | ------------------------------------------------------------------ |
| `#v` valence                 | all incident bond orders literal                                   |
| `#D` degree                  | always (topology)                                                  |
| `#X`, `#H`, `#V` totals      | implicit hydrogens literal (+ explicit H neighbors' elements literal) |
| `#d`, `#t` dative pairs      | dative overlay incidence, closed-world reading                     |
| `#a` aromatic valence        | when the aromatic system is created (removable early)              |
| `#m` multicenter valence     | when the multicenter bond is created (removable early)             |
| `#T`, `#C` stereo            | when the stereo overlay is created (design 103 flow)               |
| `#R`, `#x`, `#y` ring keys   | always (pure function of topology); checked by discharge, gated on presence |
| overlay-entity constraints   | at the resolution stage of the owning overlay                      |

## Resolution carrier

Settled 2026-08-12. The valence phase's output for an atom is a set of *completions* — one per
registry row surviving admission against the atom's effective values:

- `AtomFields` is the inherent-field completion (implicit hydrogens, lone pairs, unpaired
  electrons, spin); element and charge are fixed inputs of admission.
- `AtomCompletion` pairs `AtomFields` with the model values selection votes on: valence, donated
  and accepted pairs, aromatic valence, multicenter valence.
- The carrier maps `AtomId` to `SmallVec<[AtomCompletion; 1]>`; nearly all atoms are singletons.

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

## World assumption is per-call

Whether absence of an overlay reads as a definite negative (`NotAromatic`) or as unknown is a
property of the consuming operation, not of the molecule:

- a stored molecule-level status flag is itself a materialized projection of "is everything
  determined?" and would go stale under mutation — the same defect one level up;
- "perception ran" is history, not structure; it is not recoverable from the data, but every caller
  knows its own context (matching presents hosts as fully perceived; resolution knows it is
  staging);
- a two-type split (ground versus staging molecule) is the type-level version but is heavy and
  Python-visible, and the boundary is per key, not per molecule.

The existing `include_missing: bool` therefore becomes the reading mode supplied when
constructing the constraints view: `Complete` (the deriving relation set is complete; absence is a
definite negative) versus `Partial` (the relation set is partial; absence is no evidence).

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

### S0 — constraints-view foundation (`umol-graph-ir`); all additive, green throughout

- S0a: reading-mode enum with `Complete`/`Partial` in the view module (type name proposal:
  `RelationReading`, to confirm at review). Additive.
- S0b: `AtomConstraintsView`: construction from assertions container, atom relational context,
  reading mode, optional ring context; `effective(key)` over every `AtomConstraintKey`. Additive.
  [dep: S0a]
- S0c: comparison methods on the view — pattern-driven `matches`/`is_compatible` against a pattern
  container; the view is the receiver. Additive. [dep: S0b]
- S0d: `BondConstraintsView`, same shape. Additive. [dep: S0a]
- S0e: views for the remaining six entity families (uniform surface). Additive. [dep: S0a] —
  gates only S6a.
- S0f: nomenclature entries: view family, reading mode, effective value, discharge. Additive.

### S1 — matching on views (`umol-graph-ir::ir::substructure`)

- S1a: `host_match_targets` and the predicate closures evaluate per pattern key on the views
  (`Complete` reading); ring context built once per run as today; the `Cow` target
  materialization is deleted. Breaking, green at stage end: substructure, fingerprint, and
  reaction-matching suites unchanged. A host carrying unconsumed input assertions now meets them
  against `Complete` projections; any test relying on such a staging host is surfaced, not
  silently rewritten. [dep: S0c, S0d]

### S2 — standalone reads (`umol-graph::ops::invariant`)

- S2a: the electron-count invariant check reads through an assertion-only `AtomConstraintsView`;
  semantics unchanged. Green. [dep: S0b]

### S3 — completion vocabulary (`umol-graph::ops::valence`); additive, green

- S3a: `AtomFields`, `AtomCompletion`. Additive.
- S3b: the carrier map (`AtomId → SmallVec<[AtomCompletion; 1]>`) and the resolver report payload.
  Additive. [dep: S3a]

### S4 — pipeline rework; the one red stage, green only at its end

- S4a: atom-typing: admission against effective values via the view; no stored-constraint
  extension, no narrowed write-back; singleton atoms produce edits, plural atoms produce
  completions. Breaking. [dep: S0c, S3b]
- S4b: counts: plural `#h` candidates emitted through the same carrier; `CountsInput` readings
  aligned with the view. Breaking. [dep: S3b]
- S4c: aromaticity resolve: joint selection per candidate system over the carrier;
  `compare_valence_preference` demoted to a visible last-resort tie-break. Breaking. [dep: S4a]
- S4d: `Resolver::resolve`: thread the carrier valence → aromaticity; `Solution<(), _>` →
  `Solution<Report, _>`; rollback paths preserved. Breaking. [dep: S4c]
- S4e: SMILES ingest stops pinning `#h0` on bare aromatic heteroatoms (`umol-io`). Breaking.
  [dep: S4d]
- S4f: `umol-py`: bind the report for inspection; update the resolution verdict mapping.
  Breaking. [dep: S4d]
- S4g: retire `derive_constraints` and `include_missing` (last callers died in S1a/S4a).
  Breaking. [dep: S4a]
- S4h: acceptance: the doc 174 regression triple and the pyrrolyl DSL case; full suite green
  (`--all-features --tests`, clippy). [dep: S4a–S4g]

### S5 — discharge and output cleanup; one conformance regeneration

- S5a: per-key determination checks on the view layer; the determination table is confirmed row
  by row here. Additive. [dep: S0b]
- S5b: the discharge pass as the closing `Resolver` stage in the same journal; `⊥` ⇒
  `Contradictory`. Breaking — resolved outputs lose stored assertions. [dep: S4d, S5a]
- S5c: lowering: retire the `zeroed()` elision-only paths; raise-side dialect filling stays.
  Breaking. [dep: S5b]
- S5d: regenerate conformance snapshots once; final green: `--all-features --tests`, clippy.
  [dep: S5b, S5c]

### S6 — deferrable

- S6a: `ConstraintValidator` internals on the entity views; semantics unchanged. [dep: S0e]
- S6b: the F420 acceptance case; requires the doc 174 zero-contributor perception fix and the
  `#a0`-carbon registry coverage, tracked there.

Critical path: S0 → S1 → S3 → S4 → S5. S1 must precede S5 because discharge strips assertions
from resolved hosts, after which matching must project every pattern key. S2 floats after S0; S3
is parallel to S1. The core deliverable — staleness-free semantics and the whitepaper Mutation
story — completes at S5; S6 is not required for it.

## Open items

- Reading-mode enum type name (S0a) at review.
- F420 enablement via the doc 174 zero-contributor and registry-coverage items (S6b).
