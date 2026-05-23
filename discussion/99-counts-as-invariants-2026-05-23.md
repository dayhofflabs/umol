# Counts model on top of invariants

Status: Settled — ready for implementation.

## Decisions

| #   | Decision                                                                                              |
| --- | ----------------------------------------------------------------------------------------------------- |
| 1   | Counts restated as invariant enumeration + min-u + max-n (see below).                                 |
| 2   | `ValenceTable` field renames only: `allowed_valences` → `valence_set`, `allowed_aromatic_valences` → `aromatic_valence_set`. No schema change. |
| 3   | `NormalValenceTable` removed. Doc 96 step 7 folds into this redesign.                                 |
| 4   | Counts and AtomTyping coexist. AtomTyping covers multicenter / dative / specific Lewis modes Counts cannot express. |
| 5   | Min-u and max-n apply per atom (matches existing per-atom resolve loop).                              |
| 6   | `ValenceInvariants::solve_atom` returns `Vec<AtomAst>` (ground candidates). Counts iterates for min-u + max-n. JointDomain return revisitable once doc 97/98 operations settle. |

## Context

Doc 96 step 4 calls for `ValenceModel::Counts` to expose `resolve` / `validate` methods, with `Counts::resolve` delegating to `ValenceInvariants::solve`. Inspecting the current Counts shows three semantic mechanisms layered into the implementation:

1. `ValenceTable.allowed_valences` (RDKit-style per-element valence list) used for non-aromatic σ-sum picking, via `compute_implicit_hydrogens` (`umol-graph/src/ops/valence/table.rs:127`).
2. `ValenceTable.allowed_aromatic_valences` used as the enumeration list for aromatic-valence in the aromatic branch.
3. A separate per-element aromatic-h table (currently `NormalValenceTable.aromatic_normal_valence_for`) defaults implicit_h for aromatic atoms when the input leaves h `Undetermined`.

The non-aromatic branch picks `v_total = first allowed ≥ topology_v`. This is equivalent to **max lone_pairs**: with everything else pinned, the conservation equation `ve − q − v_total = h + 2·lp + u` has `h = v_total − topology_v` already determined per trial, so leftover for (lp, u) shrinks as v_total grows. Picking min v_total maximises leftover, and with u minimised by parity, max lp.

The aromatic branch picks via the separate aromatic-h table. With h fixed, the aromatic_valence is enumerated and the (lp, u) split follows from parity. Empirically (resolution conformance suite, pyrrole + pyridine), the same outcome falls out of **min unpaired then max lone_pairs** over the joint (h, av) enumeration without the second table. See trace below.

The mechanism that *actually* disambiguates is therefore two-rule (min-u then max-n) applied over the invariant-filtered enumeration, with element-level ranges optionally narrowed by per-element enumeration sets. The three-table layout is a precomputed shortcut for two of those rules.

## Proposed restatement

`Counts` resolves a single atom by:

1. Enumerate `(charge, implicit_h, lone_pairs, unpaired, valence, aromatic_valence, donated_pairs, accepted_pairs, multicenter_valence)` subject to:
   - Pinned (`Lit`) fields use the pin.
   - `Undetermined` fields range over the element-derived bound (`umol-shared::element` helpers), optionally narrowed by `ValenceTable` entries (`valence_set`, `aromatic_valence_set`).
   - Counts model fixes `donated_pairs = accepted_pairs = multicenter_valence = 0` when unpinned (model identity).
2. Filter by orbital == electron (`ValenceInvariants::solve_atom`'s conservation equation).
3. Apply min-unpaired (rule i): drop candidates whose unpaired is non-minimal among the filtered set.
4. Apply max-lone-pairs (rule ii): drop candidates whose lone_pairs is non-maximal among the survivors.
5. Aggregate: 0 → `NoMatch`, 1 → narrow, >1 → `Ambiguous`.

Rules i and ii apply only to `Undetermined` fields. A user pin on `#u` or `#n` overrides the rule (the pin enters step 1 as the value used for enumeration).

## Pyridine trace under the restated rules

Input: N at topology_v=2, q=0, h `Undetermined`. Equation: `h + 2·lp + u = 5 − 0 − 2 − av = 3 − av`.

Enumeration over (av ∈ {1, 2}, h, lp, u) satisfying invariant:

| av | h | lp | u |
| --- | --- | --- | --- |
| 1 | 0 | 1 | 0 |
| 1 | 1 | 0 | 1 |
| 1 | 2 | 0 | 0 |
| 2 | 0 | 0 | 1 |
| 2 | 1 | 0 | 0 |

Min-unpaired drops u=1 rows → three survivors: `(1,0,1,0)`, `(1,2,0,0)`, `(2,1,0,0)`.
Max-lone-pairs picks `(av=1, h=0, lp=1, u=0)` → pyridine ✓.

Pyrrole input pins h=1: rows with h ≠ 1 drop out; min-u leaves `(2,1,0,0)` → pyrrole ✓.

Both resolve via min-u + max-n without a separate aromatic-h table.

## What the table still does

After consolidating, `ValenceTable` entries serve only to **narrow enumeration ranges**, not to encode disambiguation. Per `(element, charge)`:

```rust
pub struct ValenceEntry {
    pub valence_set: Vec<u8>,              // narrowed v_total options; empty = element default
    pub aromatic_valence_set: Vec<u8>,     // narrowed aromatic_v options; empty = element default
}
```

Naming follows `ValueAst::Set` — the table entries are per-element narrowings expressible as `ValueAst::Set` predicates on the corresponding atom fields. Empty vec means "no narrowing — fall back to the element-level range from `umol-shared::element`."

Without the table at all, Counts still works (enumerates the full element-default range). The table is a performance / chemistry-knowledge layer, not a semantic requirement.

## Pyrrole nuance: input pin vs. model rule

The current Counts resolves pyrrole correctly because the input pins `#h` explicitly (Lit(1)). The model never has to disambiguate h=1 (pyrrole) from h=0 (pyridine) — the input tells it. Without the pin, max-n picks pyridine-mode (lp=1) over pyrrole-mode (lp=0). This is **not a Counts bug** — it is the model expressing a chemistry-justified preference (pyridine-N is more lp-rich than pyrrole-N) that the user can override with a pin.

Implication: Counts is a **default-providing model**. The user supplies pins for cases where the default chemistry-preference is wrong. AtomTyping is a **mode-fixing model** where the registry entries themselves encode the chemistry choice. The two address the same disambiguation problem at different levels.

## Doc 96 sequencing impact

- Step 4 (`ValenceModel` API methods, two-variant dispatch): implement `Counts::resolve` per the restatement above (invariant enumeration + min-u + max-n); `Counts::validate` calls `ValenceInvariants::check`. `AtomTyping::resolve` / `validate` per doc 96.
- Step 5 (collapse resolvers): unchanged.
- Step 6 (collapse validators): unchanged.
- Step 7 (`NormalValenceTable` removal): folds into step 4 of this redesign.

## Critical files

- `umol-graph/src/ops/valence/table.rs` — schema change (`valence_set` / `aromatic_valence_set` field renames; `compute_implicit_hydrogens` removal)
- `umol-graph/src/ops/valence/normal_valence.rs` — removal
- `umol-graph/src/ops/valence/counts.rs` — replace with invariant-based loop + min-u + max-n
- `umol-graph/src/ops/valence/invariants.rs` — possibly extend to accept narrowing sets per field (or Counts wraps externally)
- `umol-graph/src/ops/config.rs` — `ValenceModel` enum shape decision
- `umol-graph/config/default-valence-table.toml` — rename TOML keys; add per-charge entries to replace charge-shift fallback
- `umol-graph/config/default-normal-valence-table.toml` — removal
- `discussion/96-valence-resolution-plan-2026-05-21.md` — step 4–7 sequencing update
