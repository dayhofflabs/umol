# 239 — Potential removal of kindless stereo

Status: Proposed
Date: 2026-10-17
Relates: [238](238-non-cascading-removal-2026-10-16.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md),
[Python API guide](../docs/development/python-api.md),
[DSL specification](../umol-graph-ir/spec/umol-dsl-spec.md)

## Scope and disposition

Record the findings of the isolated stereo experiment for a possible future
removal of kindless stereo forms and their Default implementations. The experiment
supports some simplifications but does not establish an adoptable implementation.
Production code is unchanged. There is no implementation plan or commitment to
complete this redesign.

The discussion distinguished three changes:

1. Preserve actual entity removal/addition identity in reaction application and
   composition by using the complete edit-produced correspondence.
2. Require reaction modifications to preserve stereo kind, with kind changes
   expressed as removal/addition. General Edits and direct mutation remain broader.
3. Remove kindless stereo forms and the associated Default and wildcard handling.

The first belongs to [238](238-non-cascading-removal-2026-10-16.md). On 2026-10-17,
the decision there returned to the smallest tracking correction and directly
related consistency fixes, retaining the uniform entity interfaces. The second
change is not being pursued in that work. This document concerns the third change.
Reaction kind preservation does not require removing kindless forms, and removing
those forms does not itself prohibit a reaction from changing one known kind to
another.

The experiment combined all three changes. Its reaction and timing results must
not be presented as evidence that the third is necessary. The shared delta/update
interface is retained in 238; this record does not authorize its redesign.

## Motivation and capability trade-off

Stereo kinds describe different kinds of stereogenic entity within the common
stereo-atom and stereo-bond entity families. A definite kind determines the
interpretation of the coset and entity constraints. Requiring it would remove
missing-kind branches in readings, normalization, frame handling, and the DSL.

The removed state would be StereoConfigurationForm::Undetermined, which currently
asserts a stereo entity without selecting its geometry. This is a real loss of
pattern expressivity: omitting the entity does not assert that it exists with an
unknown kind. The proposal also removes that cross-kind wildcard/generalization.
Known-kind unknown cosets, such as Th* and Ct*, remain meaningful.

Default is not a semantic reason to retain a kindless form. Normalize and
FrameTransport do not require Default, and Lattice already accommodates types
without a common top. The concrete ownership dependency is in collection and
span code using mem::take. Parsing need not manufacture a kindless placeholder.

The implemented kindless constraint-vacuity treatment in
[181](181-python-boundary-ownership-2026-08-03.md#kindless-stereo-constraints--decision-2026-10-14)
remains current unless a later redesign is adopted. Its historical decisions are
not rewritten by this proposal.

## Candidate exercised by the experiment

The following are trial choices, not settled public interfaces:

- StereoConfigurationForm becomes a value with public kind and coset fields and
  a constructor requiring StereoKind. StereoAtomForm and StereoBondForm retain
  their configuration and constraints. These forms have no Default.
- Unknown coset remains top within a known kind. Different satisfiable kinds
  meet to bottom and have no join through the existing NoJoinError channel.
- StereoConfigurationUpdate carries a required kind and optional coset. An
  omitted coset preserves the current coset for the same kind; a general update
  that changes kind without supplying a coset starts with an unknown coset.
  StereoAtomUpdate and StereoBondUpdate require that configuration, including for constraint-only
  updates. Omission remains distinct from explicitly assigning an unknown coset.
- Stereo constraint changes carry definite kind context. Individual constraints
  do not gain their own intrinsic kind. EntityPatch gains a Context associated
  type, unit for ordinary entities and StereoKind for stereo; that context adds
  no runtime allocation.
- Direct mutation and general Edits still permit changes between known kinds,
  including incompatible raw cosets and constraints. Attribute validity does
  not become a mutation failure boundary or a new molecule-integrity requirement.
- The combined trial rejects kind-changing modifications at the existing
  Reaction::try_new boundary. Raw delta values retain complete old/new
  configurations and can still express such changes. The experiment therefore
  does not establish a delta representation incapable of changing kind.

Empty collections retain Default. Generic EntitySpan need not gain an empty
variant to replace the lost Default on stereo attributes. A degree-zero Any or
Null kind introduced solely as a default is not part of the proposal.

The fixed-kind TetrahedralStereoForm and CisTransStereoForm, including their
Undetermined and NotStereo alternatives and the #T/#C notation, remain. Missing
overlays, missing constraints, absent span sides, and bottom readings retain their
separate meanings.

## Findings

### Moving attributes without Default

The trial replaces mem::take with a temporary value containing the original
kind, an unknown coset, and empty constraints. The complete original attributes
are moved out and normalized or transported; the temporary does not determine
the frame action or participate in the operation. Both sides of a Modified span
receive the same action.

Source inspection found no payload clone or heap allocation introduced by
constructing this temporary. Existing Arc detachment and normalization work
remain. These operations consume their owning collection; failure drops that
working value, and surviving clones retain their original attributes and frames.
No graph-core storage redesign or additional rollback mechanism was needed.

This resolves the ownership obstacle in the exercised paths. It does not prove
every normalization or frame law, and it is not a measured speedup claim.

### Trial join rule fails on raw contradictions

The trial treats one contradictory operand as bottom and returns the other
normalized value. It returns NoJoinError when both operands are contradictory.
That extension is not associative on raw contradictory inputs.

Let b be a tetrahedral configuration with an empty literal coset set, and v be
Th0. Normalizing b returns Contradiction. Then:

```text
b.join(b).and_then(|x| x.join(&v)) = Err(NoJoinError)
b.join(&v).and_then(|x| b.join(&x)) = Ok(v)
```

The same issue arises with entity-constraint contradictions. Passing laws over
normalizable values do not settle this case. The raw-input domain and return
contract require revision before adoption; the trial rule is not acceptable as
an associative extension. This finding does not imply that the whole workspace
needs a new lattice API.

The baseline uses a global-undetermined fallback when a configuration operand
cannot normalize. Removing that state requires an explicit replacement contract;
substituting a per-kind unknown value mechanically does not establish the laws.
This question is independent of the temporary replacement used for ownership.

### Existing symbolic-expression/frame defect

A new property comparing normalization before and after frame transport fails
for SquarePlanar, a free coset variable s, and successive permutations
[0,3,2,1] and [1,3,0,2]. Normalizing after the first action and then applying the
second yields Apply(s, [1,3,2,0]); applying both actions before normalization
yields Apply(s, [3,1,0,2]).

An independently built baseline reproduces exactly the same result. This defect
predates the candidate representation and is not caused by removing Default.
The failing property, regression seed, and minimized standalone probe were kept;
no symbolic-algebra repair was included. Its disposition is separate from the
decision to remove kindless forms.

### Reaction results belong to the combined trial

Focused checks established that same-kind modification retains identity, while
explicit removal/addition remains unmatched even when participants and kind
coincide. These cases covered stereo atoms and bonds, realized reactions and
spans, and stereo-atom composition. General known-kind mutation remained
available. Application reused Molecule::tracked_apply's complete correspondence;
composition used complete correspondences instead of reconstructing identity
from atom incidence.

Those results support the tracking route and demonstrate the tested reaction
restriction. They do not establish a need to remove kindless forms or justify
adopting that restriction at the cost of the uniform entity interfaces retained
in 238.

## Verification evidence and limits

The experiment used independent source snapshots of commit
1f782a5f1147e7c70165bd4e84ba1fbae998d33c. The baseline library source was unchanged;
only common benchmark and diagnostic instrumentation was added there.

| Check | Recorded result |
| --- | --- |
| Baseline graph-IR unit suite | 11,293 passed; 3 ignored |
| Baseline stereo-filtered properties | 79 passed |
| Existing trial property suite, before the added diagnostic | 482 passed; 1 existing ignored test |
| Trial canonicalization, integrity, and span integration targets | 54 passed |
| Final stereo-filtered properties | 86 passed; the new expression/frame property failed |
| Independent baseline expression/frame probe | Reproduced that failure |
| Trial library Clippy with warnings denied | Passed |
| Trial full unit suite | Did not compile: old APIs remained in fixtures and expectations |

The full property run preceded six additional passing focused cases and the
failing expression/frame property. Python, downstream-crate, doctest, and full
workspace gates were not run. The prototype is not a green workspace migration.

The changed test domain was explicit: kindless generator branches were removed;
known kinds, unknown cosets, literal/set ranges, raw constraints, and general
kind-changing edits remained. Kindless accessor/reference fixtures used explicit
kinds with unknown cosets while retaining their original assertions. Kindless
export-elision cases became required-kind roundtrips. Reaction kind-mismatch
rejection moved to construction; the separate discontinuous-field application
failure property remained.

An accidental narrowing of the normalized configuration generator's literal/set
range was caught and restored before the final focused run. No tests were
disabled to produce the reported results. Existing tests were not all unchanged,
and the unfinished unit migration must not be represented as passing coverage.

## Bounded performance evidence

The first comparison was invalid: sharing one Cargo target directory allowed a
baseline rebuild to reuse the trial library. All figures from that comparison,
including an apparent composition regression, were discarded.

The valid runs used separate baseline/trial target directories, Rust 1.99.1 on
Darwin 24.6.0 arm64, and the optimized Cargo bench profile. Benchmark additions
were identical. Each case used 20 Criterion samples, one second of warm-up, and
one second of measurement. There was one valid baseline/trial pair per case.

| Case | Baseline sample median | Trial sample median | Change |
| --- | ---: | ---: | ---: |
| tracked_apply_at, 5 atoms | 3.208 µs | 2.449 µs | -23.7% |
| tracked_apply_at, 128 atoms | 21.790 µs | 5.739 µs | -73.7% |
| Realized span, 5 atoms | 12.173 µs | 11.502 µs | -5.5% |
| Realized span, 128 atoms | 135.412 µs | 119.662 µs | -11.6% |
| Composition, 5 atoms | 737.081 µs | 703.452 µs | -4.6% |

Application used one tetrahedral stereo atom; the 128-atom host added a chain
outside the five-atom match. Composition used a same-kind change and its inverse
on five distinct elements over all enumerated overlaps. These workloads do not
establish production frequencies or scaling beyond those cases.

All five comparisons favor the combined trial. Eliminating incidence-based
correspondence reconstruction is a plausible source of the application gains,
but the measurements do not isolate it from the representation and reaction
changes. No compiler attributes or layout tuning were introduced. These results
cannot justify removing kindless forms on performance grounds, or substitute for
measuring the eventual narrower tracking implementation.

## Migration surface if reconsidered

This is a scope inventory, not an implementation sequence. Reconcile it with
the actual outcome of 238 before reopening the design.

| Surface | Consequences to address |
| --- | --- |
| Forms and updates in graph-IR stereo | Constructors, variants, Default, kind/coset access, update/difference laws, conversions, equality/hash/order, literal extraction, normalization, lattice operations, and frame transport |
| Stereo collections and spans | Replace Default-dependent moves without changing entity identity, participant frames, or failure/clone behavior |
| Entity constraints, views, and projections | Remove missing-kind bypasses and vacuity branches; retain owning-kind interpretation, Option-for-bottom readings, and existing projection-kind restrictions |
| Molecule-level stereo constraints | Retain their independently supplied definite kind and existing integrity contract, including logical constraints |
| Edits, updates, deltas, and EntityPatch | Preserve general mutation and journal restoration; reconcile field-change/update shapes, diff laws, and kind context with the reaction design settled in 238 |
| Matching, correspondence, canonicalization, and symmetry | Remove missing-kind branches without broadening reaction identity rules to arbitrary correspondences; preserve open-attribute failure behavior and frame-aware distinctions |
| DSL forms, updates, edits, reactions, and spans | Migrate compact syntax, tree and streaming parsing, rendering, and lowering together; reject removed syntax rather than inventing a kind or silently dropping payloads |
| Python | Mirror the final Rust construction and observation contracts, conversions, stubs, and retained immutable payloads; keep live mutation and ordinary Python 3.9 access |
| Chemistry and format consumers | Migrate stereo perception/resolution/validation, raising, layout, and depiction consumers; keep unknown cosets and external-format policies |
| Public exports and normative documentation | Reconcile data-type, integrity, nomenclature, Python API, and DSL specification contracts with the actual adopted interface |

The Python audit was source-only. The candidate turns the existing kind and coset
properties into total readings and removes the Undetermined/Kinded configuration
variants. StereoAtomForm and StereoBondForm constructors already require a
configuration; their update constructors would newly require one. Definite kind
context affects four constraint Edit variants, two Delta variants, and two batch
builder methods. These are consequences of the trial shape, not approved binding
changes or reasons to add new Python-only names.

The experiment required no graph-core or umol-perm representation changes.
Additional geometries, stronger molecule attribute-integrity checks, and manuscript
edits are outside this proposal.

## Contracts and evidence still needed for adoption

- Resolve joins on contradictory raw operands and preserve the stated laws over
  their explicit domain. Settle whether is_undetermined means top within a kind;
  it must not erase kind or entity existence.
- Settle final constructors, fields/accessors, updates, and constraint kind
  context. Constraint-only and no-op updates must serialize faithfully. Removing
  kindless entity forms alone does not supply context to independent changes.
- Preserve infallible raw attribute mutation and the distinction between omitted
  updates and explicit unknown cosets. Do not add a global top, a dummy kind, or
  fallibility merely to accommodate an implementation helper.
- Complete the unit and downstream migration. Preserve normalizing, lattice,
  frame action/composition, update/difference, DSL roundtrip, and Python parity
  laws. Exercise all kind domains and incompatible raw attributes permitted by
  the open-form contract. Account explicitly for the existing frame failure.
- Measure any identified additional allocation or traversal in the chosen
  implementation. The temporary replacement needs no separate benchmark
  campaign; a claimed performance benefit must be attributable to this change.

## Experiment artifacts

The durable findings are recorded above. Temporary source snapshots, the patch,
raw samples, logs, and diagnostic probes remain under
[scratch/238-stereo-experiment](../scratch/238-stereo-experiment/experiment.md).
The independent build targets were /Users/dr/.cargo-target/238-stereo-baseline
and /Users/dr/.cargo-target/238-stereo-trial. Artifacts marked discarded-* or
discarded-shared-target belong to the invalid initial comparison.
