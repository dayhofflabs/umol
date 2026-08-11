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
| input assertions            | raised SMILES `#a` negatives, user-supplied `#v4`   | stored; consumed at the commit that realizes them (aromaticity does this now) |

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
reaches the molecule, as primary data, with the input assertions it realizes consumed in the same
transaction. The stored channel thus holds only input assertions, with consumption semantics, not
cache-invalidation semantics; nothing derivable-from-structure is ever stored.

## Read-side restructure

The `Lattice` trait is untouched: `meet`, `join`, and `matches` on values are already
source-agnostic. What restructures are the container-level comparison entry points (`matches`,
`is_compatible`, match-target construction), which currently demand two materialized containers.
They become generic over a per-key value source:

```
effective(key) = asserted(key) ∧ projected(key)    // meet; either side may be absent
```

- Evaluation is driven by the pattern's keys, so only requested projections are computed. The
  per-key selector of doc 125 falls out, and the monolithic `derive_constraints` retires.
- A `⊥` effective value needs no special path: `matches`/`is_compatible` return false.
- A standalone entity form (no molecule context, e.g. the electron-count invariant check) is the
  same abstraction with an assertion-only source.
- Ring keys obtain their `RingViews` context inside the source, as `ConstraintEvaluation` already
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
- `ConstraintValidator`: already per-key derived-versus-asserted; may reuse the source.

## Consumption at commit

Only input assertions are ever stored, so consumption concerns them alone. The commit that
realizes an assertion into relations removes it in the same transaction. The aromaticity resolver
already does exactly this for `#a`; the kekulizer and aromatizer removals are the same contract.
The final commit of a resolution consumes the input assertions its chosen completion realizes
(`#v4` once implicit hydrogens and all incident bond orders are literal). Keys no single commit
owns are verified and removed by a final normalization stage, gated on presence.

Determination criteria per key (to be confirmed row by row during implementation planning):

| key                          | determined when                                                    |
| ---------------------------- | ------------------------------------------------------------------ |
| `#v` valence                 | all incident bond orders literal                                   |
| `#D` degree                  | always (topology)                                                  |
| `#X`, `#H`, `#V` totals      | implicit hydrogens literal (+ explicit H neighbors' elements literal) |
| `#d`, `#t` dative pairs      | dative overlay incidence, closed-world reading                     |
| `#a` aromatic valence        | consumed when the aromatic system is created                       |
| `#m` multicenter valence     | consumed when the multicenter bond is created                      |
| `#T`, `#C` stereo            | consumed when the stereo overlay is created (design 103 flow)      |
| `#R`, `#x`, `#y` ring keys   | always (pure function of topology); verified and removed by final normalization, gated on presence |
| overlay-entity constraints   | at the resolution stage of the owning overlay                      |

Direction of the check: the projection `matches` the assertion ⇒ redundant ⇒ remove; meet is `⊥` ⇒
resolution error; assertion strictly narrower than a not-yet-determined projection ⇒ keep.

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

The existing `include_missing: bool` therefore becomes a named reading mode supplied when
constructing the value source.

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
- Raising is unchanged: input dialect assertions land as staging and are consumed by resolution.
- Whitepaper: the Mutation section states that edits operate on relations and inherent fields;
  projections are computed, so they cannot go stale; stored assertions may become unsatisfiable,
  which validation reports; resolution consumes assertions the structure comes to determine.

## Open items

- Naming, subject to the nomenclature process: the combined per-key read (candidate: *effective
  value*); the per-entity value source abstraction (candidates: `*ConstraintSource`,
  `*ValueSource`); the realization-time removal (candidate: *consume*); the final pipeline
  normalization stage (candidate: *discharge*); the reading mode (candidates: perceived/staging,
  closed/open world).
- Confirmation of the per-key determination table.
- The candidate-set carrier between resolution phases: its types (doc 053 shape; the existing
  `Solution`/`Progress` vocabulary), how surviving multiplicity surfaces in verdicts, and the doc
  174 regression triple as its acceptance test.
- Staged implementation plan, to be added here once the above are fixed.
