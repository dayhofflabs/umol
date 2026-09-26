# 231 — View arguments in operations

Status: Proposed
Date: 2026-09-25
Relates: [213](213-editor-overlay-storage-2026-08-27.md),
[149](149-molecule-ring-cache-and-hashing-2026-07-13.md),
[086](086-molecule-ast-api-2026-04-16.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Purpose

Narrow review of the workspace for functions that take a view as an argument and for structs that
hold one, against the nomenclature guide's rule for views. The rule has held since views were
introduced, yet such signatures keep reappearing. This document records every occurrence as
one work unit so that none is fixed piecemeal or overlooked again. Doc 213 S2f is cancelled;
its two listed functions are no longer scheduled for correction there. This document records
facts and the required outcome; it contains no implementation plan. Line references are at
commit `d19648bae`, at which none of the cited files had uncommitted changes.

## Rule

[Nomenclature guide](../docs/development/nomenclature.md) §"View and Views":

> **Views are receivers, never arguments.** A function takes ids and the molecule, or takes the
> owned `*Form`; it does not take a view. A view borrows its molecule and exists to be called *on*,
> so passing one propagates a borrow through a signature that did not need it.
>
> **Views are presentation facades.** A view never takes another view as an argument nor holds one
> as state, and its implementation lives in functions of the molecule and entity ids beneath it.

## Method

Every `fn` item in `umol-graph-ir`, `umol-graph`, `umol-io`, `umol-py`, `umol-graph-core`,
`umol-geometric`, and `umol-geometric-graph` (sources and tests) was scanned with a signature
parser that handles multi-line parameter lists, generic parameters, and `where` clauses. A
parameter or bound counts when its type names an identifier ending in `View`, `Views`, or
`ViewMut`; `self` receivers are excluded, since a view's own methods are the intended use. Two
line-based greps independently reproduced the same 14 functions; a field grep over struct
definitions found the one stored view. The working tree was read only.

## Findings

### umol-io — depiction and layout

All five sit in the depiction and layout path and share one shape: the molecule is passed
alongside a view obtained from it.

| Function | Signature (view parameter) | Callers | Scheduled |
| --- | --- | --- | --- |
| `depict/molecule.rs:224` `tetrahedral_candidates` | `stereo: StereoAtomView<'_>` beside `molecule: &Molecule`, `layout: &MoleculeLayout` | `:186` | no (S2f cancelled) |
| `layout/coordgen.rs:60` `cis_trans_bond` | `stereo: StereoBondView<'_>` beside `molecule: &Molecule` | `:43`, `:234` | no (S2f cancelled) |
| `depict/molecule.rs:350` `atom_label` | `atom: AtomView<'_>` (sole parameter) | `:91` from a view iteration; test `:1083` | no |
| `depict/molecule.rs:405` `aromatic_contour` | `system: AromaticSystemView<'_>` beside `molecule`, `layout` | `:102` | no |
| `depict/molecule.rs:484` `aromatic_annotation` | `system: AromaticSystemView<'_>` beside `molecule`, `layout`, `contour` | `:105` | no |

### umol-graph — constraint validation

| Site | Shape | Callers |
| --- | --- | --- |
| `ops/validate/constraint/ring.rs:177` `validate_atom_constraints` | `view: &RingAtomView<'_>` beside `atom_id`, `constraints` | `:40`, `:112`; both build the view at the call site as `&rings.atom(id)` |
| `ops/validate/constraint/ring.rs:202` `validate_bond_constraints` | `view: &RingBondView<'_>` beside `bond_id`, `constraints` | `:53`, `:144`; same pattern with `&rings.bond(id)` |
| `ops/validate/constraint/relational.rs:344` `evaluate_atom`, `:391` `all_atoms`, `:404` `any_atom` | `rings: Option<&RingViews<'_>>` beside `molecule`, ids, `predicate` | the relational entry computes `molecule.rings(..)` once and threads `rings.as_ref()` through these three helpers at 20 call sites |
| `ops/validate/constraint.rs:209` `ConstraintEvaluation` | field `rings: Option<RingViews<'a>>` beside `molecule`, `config`, `component_by_atom` | the lazily perceived ring context of one validation pass, held as a view |

The `ConstraintEvaluation` field is the "holds one as state" clause of the rule rather than an
argument; it is the same category of misuse and the source of the threaded `Option<&RingViews>`.

These views are not molecule views. `RingAtomView` and `RingBondView` hold a `RingSet` and the
molecule, and every quantity the validators read from them (`ring_degree`, `ring_valence`,
`ring_membership(scope)`) already delegates to a derivation layer of free functions in
`ir/view/ring.rs`: `atom_ring_membership(&RingSet, AtomId, RingScope)`,
`atom_ring_degree(&Molecule, &RingSet, AtomId)`, `atom_ring_valence(&Molecule, &RingSet, AtomId)`,
and `bond_ring_membership(&RingSet, BondId, RingScope)`. The layering the guide asks for therefore
exists; the layer is `pub(crate)`, used by the ring views and by `AtomConstraintsView` and
`BondConstraintsView`, and unreachable from umol-graph. `RingSet::enumerate` is `pub(super)` and
takes graph-core's `Graph`, so the only public route to an owned `RingSet` is
`molecule.rings(model, config).into_ring_set()`.

Assessment: the ring-view arguments are a workaround, not a design. The validators reached the
ring derivations through the only public path, the views, instead of exposing the derivation
layer, and the threaded `Option<&RingViews>` and the stored `RingViews` followed from that
choice.

### umol-graph-ir — `GraphView` isomorphism family

`ir/view/graph.rs:249` `visit_subgraph_isomorphisms`, `:278` `enumerate_subgraph_isomorphisms`,
`:296` `visit_subgraph_isomorphisms_at`, and `:327` `enumerate_subgraph_isomorphisms_at` each take
`query: &GraphView<'_>`: a view method taking another view. The production matcher does not use
them. `ir/substructure.rs:247` and `:312` call graph-core's `Graph::visit_subgraph_isomorphisms`
directly on `Molecule::raw_graph()` and on the Levi graph. The only callers of the four methods
are the unit tests at `ir/molecule/tests.rs:3562–3655`.

### Compliant

- `umol-py`: the four Rust-side view mentions (`dative.rs:270`, `aromatic.rs:319`,
  `multicenter.rs:328`, `noncovalent.rs:391`) are return types of accessor constructors, not
  parameters.
- `umol-graph-core`, `umol-geometric`, `umol-geometric-graph`: no view parameters.
- `umol-io` TableIR `SgroupsView` and `RgroupsView`: never taken as arguments.
- No test helper in any crate takes a view; no closure bound such as `Fn(AtomView)` appears in
  any signature; no view type nests inside another view in `ir/view`.
- `ConstraintsViewMut` (doc 213 S2a) is a view handed out as a mutable object rather than taken
  as an argument; S2m already removes it.

## Required change

Every function above loses its view parameter and receives what the view wraps, obtaining any
view it needs inside. `ConstraintEvaluation` stops holding a view as state. The inventory
contains fourteen functions and one field; no pair is scheduled under doc 213 following
S2f's cancellation. The cases fall into three groups:

- **Molecule and id (the five umol-io functions).** The view parameter becomes the entity id;
  the four that already take the molecule change nothing else, and `atom_label` gains the
  molecule.
- **Molecule, `RingSet`, and id (the umol-graph validators and the field).** The validators
  receive `(molecule, &RingSet, id)` and call the derivation layer directly; the relational
  helpers thread `Option<&RingSet>` beside the molecule they already carry; `ConstraintEvaluation`
  holds the owned `RingSet`. This requires the derivation layer to become public, which is the
  one interface decision below.
- **Remove or redefine (the `GraphView` isomorphism family).** With no production caller, the
  decision is whether the four methods exist, not how they take their query.

## Open questions

1. **Scheduling of the depiction/layout functions.** No changes to these five functions are
   scheduled under doc 213 after S2f's cancellation. Their scheduling remains open here;
   cancellation does not approve a replacement implementation item.
2. **The public home of the ring derivations.** The four `pub(crate)` functions in
   `ir/view/ring.rs` must become reachable from umol-graph with signatures over the molecule,
   the `RingSet`, and the id. The alternatives are public free functions exported from `ir`;
   methods on `RingSet`, two of which take the molecule as an argument; or methods on `Molecule`
   taking the `RingSet`. The ring views and the constraints views keep delegating to whichever
   form is chosen. `RingViews` remains the constructor of the owned set through `into_ring_set`
   unless `RingSet::enumerate` also changes, which this document does not propose.
3. **The `GraphView` isomorphism family.** Remove the four methods, moving their tests to the
   graph-core surface `substructure.rs` already uses, or keep them with the query molecule as
   the argument. Removal shrinks `GraphView` to the operations doc 086 listed for it;
   redefinition keeps a typed `AtomId` isomorphism entry point on the view.
