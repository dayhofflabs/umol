# 125 — What constraints represent: the projection / view model

Status: Active
Date: 2026-06-22

## Draft API

- `AtomView::derive(&[AtomConstraintKey]) -> AtomConstraints`: derive for given set of constraints
- `AtomView::derive_for(&AtomConstraints) -> AtomConstraints`: derive for constraints in container
- `AtomView::matches(&AtomConstraints) -> bool`: derive and match constraints in container
- `AtomView::is_compatible(&AtomConstraints) -> bool`: derive and check for compatibility for constraints in container
- `AtomConstraints::append(&AtomConstraints) -> Option<AtomConstraints>`: meet-add constraints in container (add semantics)

## Why this doc

Resolving stereo through atom-typing exposed a structural problem that isn't
specific to stereo: the `AtomConstraints` channel is doing several unrelated jobs
with one representation, and the resolver assumed only one of them. The immediate
case is fixed (see "Forcing example"), but the fix only works because we stopped
the registry from asserting a constraint it had no business asserting. That points
at a question we should answer deliberately rather than case by case: **what is a
constraint, and where does derived information belong?**

The useful lens is database modeling. The derived constraints (valence, degree,
ring membership, aromatic-valence, donated/accepted pairs) are **projections of
the relational structure** — the graph and the overlays. Attaching them to the
atom entity is denormalization, and the symptoms we have been hitting are the
classic symptoms of denormalization.

## Forcing example (context, already resolved)

Atom-typing valence resolution matched each atom against a registry of valence
templates and narrowed toward the match. The registry was parsed with
`AtomDefaults::zeroed()`, which sets `tetrahedral_stereo: NotStereo`. A real
stereocenter carries `#T1`. `meet(Stereo(coset), NotStereo) = None`, so:

- a stereocenter matched no registry entry (`no atom-typing match`), or
- where it did, `narrow_from` forced `NotStereo` onto it and the stereo was lost
  before the stereo stage could promote it to an overlay.

counts resolution was correct throughout because it never touches stereo: it
computes only valence fields and combines via `meet`, which *preserves* what it
does not constrain.

Resolution that converges with counts:

- the registry leaves tetrahedral stereo unconstrained (`StereoDefault::Required`
  → emit nothing), so an entry never asserts `NotStereo`;
- candidate selection uses `is_compatible` (`meet(..).is_some()` — "no conflict"),
  not `matches` ("target refines pattern"), because the atom carries the specific
  `#T1` and we only need consistency, not refinement;
- `narrow_from` (already `meet`) preserves the atom's `#T1` for the stereo stage.

This works, but the reason it had to work this way is the subject below.

## The three roles "constraint" currently muxes

| role | example | database analogue |
|---|---|---|
| projection of relations | valence = Σ bond orders; degree; ring membership; aromatic-valence; donated/accepted pairs | a **view** (or materialized view) |
| query predicate | a substructure pattern's `#v4`, `#T1` | a **WHERE clause** |
| un-normalized input | SMILES `@` → `#T` before any stereo overlay exists | **staging** data, pre-normalization |

These are genuinely different. Treating them as one ("whatever is useful to
attach to an atom") is the denormalization. The stereo bug is an **update
anomaly**: a projection (`NotStereo`) stored on the entity contradicted the fact
living in the relation.

## Principle: do not denormalize projections into entities

How efficient databases handle on-the-fly projections of base relations:

- a **view** — recomputed on read. Right for cheap projections (valence, degree).
- a **materialized view** — precomputed, stored *separately* from the base tables,
  with an invalidation/refresh contract. Right for expensive projections (rings).

In neither case is the projection written into the base entity row. That is
exactly to avoid the anomaly we hit. umol already has both primitives:

- the view accessors (`AtomView::valence`, `degree`, …) are views.
- `rings_cache` is a materialized view of a pure function over the immutable AST.

## Re-reading the two "warts"

The two cases of interior mutability in `MoleculeAst` are opposites, not two
instances of the same smell:

- `rings_cache` is the **correct** pattern: a transparent cache of a pure
  projection over an immutable base. The interior mutability is lazy-init of a
  cache; it changes no semantics, and it is the model for how an expensive
  projection should be handled — materialized *beside* the entity, never inside it.
- in-place constraint materialization (extending an atom with derived constraints
  during resolution) is the **anti-pattern**: it mutates *asserted content* — what
  the molecule claims — by writing a projection into the entity.

The line is not "interior mutability good/bad." It is **transparent cache of a
pure projection (fine) vs. mutation of asserted content (not fine).**

## Corroboration: `zeroed()` as an omission threshold

A later experiment made `AtomDefaults`/`BondDefaults::zeroed()` equal to `ground()`
(stop zeroing the constraints) and ran the suite. Resolution *behavior* was
unchanged; the only fallout was that **all 652 conformance snapshots changed** —
the previously-omitted derived constraints (`#a!`, valence, …) now rendered in the
lowered output (AST → DSL). So the zeroing is not superfluous, but only because it
is compensating for something.

That something is eager materialization. In a fully resolved molecule the derived
constraints are redundant with the topology and overlays — they are projections,
written onto the entity eagerly. `zeroed()`'s job on the lowering path is solely to
*hide that redundancy*, omitting any constraint equal to its zero value. If those
predicates were views (computed on read), there would be nothing to omit: lowered
output would naturally carry only the primary, non-derivable data. The omission
threshold is a band-aid over materialization, not a feature — the same root cause
as the resolver difficulties, seen from the lowering (output) side instead of the
raising/matching side. It also shows the cost is currently only aesthetic (the
materialized predicates are redundant-but-harmless), which is why it went unnoticed
until stereo made a materialized predicate *conflict* rather than merely repeat.

## Lifecycle this implies

```
input (staging: un-normalized constraints)
  → resolution (normalization: constraints become relations — bonds, overlays)
    → ground molecule (relations + views)
```

On a fully resolved molecule, a projection-constraint stored on the atom is
redundant with the relations; it should be a view, not a stored field. Stored
constraints are only legitimate as (a) query predicates or (b) staging not yet
normalized.

## What the model validates

- **Resolution needs only a selected few projections** — valence, aromatic
  valence, donated pairs, accepted pairs, multicenter valence — computed on
  demand, not a monolithic `derive_constraints` that materializes all of them
  (rings, stereo, noncovalent included). The others cannot participate in valence
  resolution; including them costs correctness (the stereo bug) and performance.
  counts already approximates this with `derive_fields` + `meet`.
- **A per-`*ConstraintKind` selector** ("project these columns") is the right
  shape for requesting exactly the projections a consumer needs, instead of
  bespoke `derive_constraints` / `zeroed` that hard-code the full set.
- **TableIR is a staging format**, one input dialect — not the schema. Resolution
  being shaped by the atom/bond-list view is staging defining the model, which the
  projection lens says it should not.

## The real work: demultiplex the roles

Deleting `derive_constraints` is not the hard part. The hard part is giving each
role its own channel so they stop sharing one representation:

- base relations — the graph and overlays (already exist).
- views — computed projections (accessors exist); materialized where expensive
  (`rings_cache` exists).
- a predicate representation for queries — what a substructure pattern carries.
- a staging representation for un-normalized input.

Substructure matching then becomes its clean database form: a **query (predicates)
evaluated against views (projections)** of the host. This is also the direction
the doc-122 host-borrow / don't-materialize-the-host work was already pulling.

"Constraint" is load-bearing in all three roles today, so the migration is wide,
not deep — which is the argument for staging it.

## First steps (low regret, do not commit to the full demux)

1. Selective, per-kind derivation: resolution requests only the valence
   projections it needs (valence / aromatic-valence / donated / accepted /
   multicenter), via `*ConstraintKind`, computed not materialized. Retire the
   monolithic `derive_constraints` / `zeroed` full-set assumptions in favor of
   generic per-constraint selection over the existing defaults/overrides facility.
2. Stop materializing projections into the entity during resolution (the in-place
   extend). Compute what selection needs on the fly (as counts does) and write
   only resolved primary data plus the normalized relations.

## Open questions

- Should a ground molecule store projection-constraints at all, or only expose
  them as views? If views only, what is the storage of record for, e.g., an
  asserted valence on a partial molecule (staging) versus a resolved one
  (projection)?
- Predicate vs. value: is a query predicate a different type from a stored value,
  or the same type read in a different context? The lattice `matches` /
  `is_compatible` distinction (refinement vs. consistency) is the operational
  shadow of this question.
- Where does staging live — does the AST carry un-normalized input constraints, or
  is there a pre-resolution representation distinct from the resolved AST?
- Cost model: which projections are cheap enough to be pure views and which need
  materialization (the `rings_cache` treatment) under reaction-network scale.
- Relationship to doc 123 (allocation survey) and doc 122 (read-path references):
  both are consistent with "projections are views, computed not stored," and would
  fold into this model rather than being separate optimizations.
