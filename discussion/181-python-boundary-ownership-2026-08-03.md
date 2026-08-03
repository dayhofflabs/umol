# 181 — Python boundary ownership and copy costs

Status: Proposed
Date: 2026-08-03
Relates: [179](179-python-editing-and-transactions-2026-08-02.md),
[151](151-python-molecule-workflows-2026-07-13.md)

The Rust/Python boundary should be reviewed for unnecessary cloning and reconstruction after the
editing and transaction work in doc 179 is complete. Python is expected to be the primary interface,
so ownership choices made for convenient Rust bindings must not impose avoidable recurring costs on
Python operations. This document records preliminary findings only. It is not yet an implementation
plan.

## Correction to the preliminary inventory

Counting explicit `.clone()` calls does not measure copying across the boundary. For example,
`IsotopeMassAst::from_rust` directly clones owned sets and strings, while
`ElementAst::from_rust` reconstructs sets through iteration and `collect` and clones a variable name.
Both produce detached owned Python values from borrowed Rust values. The apparent difference between
the two types is an implementation detail, not a meaningful ownership distinction.

The review must therefore inspect boundary operations rather than count methods or syntax. For every
Rust-to-Python and Python-to-Rust path, establish:

- whether the source is owned or borrowed;
- whether the result is detached, parent-backed, or shared;
- which data is copied, reconstructed, or moved;
- whether the same data is copied again on the reverse conversion;
- the expected size and call frequency of the value;
- whether the Rust operation consumes an argument only to move its contents efficiently.

Rust-to-Python and Python-to-Rust costs must be reported separately. The latter may be more important
for reusable Python configuration and edit values.

## Preliminary ownership classes

### Detached values obtained from borrowed parents

Field values returned from molecule views, table entries returned from a valence table, and individual
entries obtained from `Edits` or `Deltas` are borrowed on the Rust side but represented as detached
Python values. Those values must acquire ownership somewhere. Copying is expected under the current
semantics.

Small leaf values such as `ElementAst`, `IsotopeMassAst`, and `ValenceEntry` are unlikely to justify
parent-backed Python views. Such views would complicate lifetime, mutation, equality, and detachment
semantics. The review should still measure these paths and identify accidental duplicate conversions,
but it should not assume that eliminating every copy is desirable.

Edit and delta materialization similarly copies the entity AST payloads stored in borrowed container
entries. Returning views instead would change the semantics of indexing and iteration, especially for
containers that can subsequently be extended. Preserve detached results unless measurements justify a
different public contract.

### Owned results passed through borrowed conversion helpers

Some callers may own a Rust result but pass it to a conversion helper taking `&T`, which then clones or
reconstructs it. These copies should be removable by consistently distinguishing owned conversion from
borrowed projection. This is expected to be an internal binding cleanup with no visible Python API
change and usually no public Rust API change.

### Large immutable model data

`AtomTypeRegistry` and `ValenceTable` require focused review. Built-in instances are ordinarily borrowed
static data, but custom instances are owned and may be large. The current Python `ChemistryModel`
contains Python-owned model values and converts them back into a `ChemistryModel` whose custom registry
or table is `Cow::Owned`. Repeated operations with one Python model can therefore deep-copy reusable
model data.

This is the strongest candidate for shared or borrowed boundary representation. Possible solutions
include shared ownership for immutable model data, an operation input that borrows the Python-owned
model for the duration of the call, or a binding-specific conversion path that does not first construct
another owned aggregate. The review must examine the current `Cow<'static, _>` representation rather
than assume that cloning is unavoidable.

This area may require a targeted public Rust API change while leaving the Python call shape unchanged.
The default static model and custom owned models must both retain straightforward construction and use.

### Reusable Python values passed to consuming Rust operations

Rust transaction application consumes `Edits`, allowing its payloads to move into the editor. Python
normally treats an `Edits` argument as reusable, so the binding clones the full container before calling
the consuming Rust operation. Changing this has real semantic trade-offs:

- consuming the Python wrapper would make an ordinary argument one-shot;
- making Rust application borrow the batch could move copying into the application of each payload;
- shared or copy-on-write storage would add representation complexity.

The same analysis applies to other Rust operations that consume values supplied by reusable Python
objects. Preserve the current behavior unless measured costs justify changing either API. In
particular, do not weaken efficient Rust move semantics merely to remove a clone visible in the binding
source.

## Expected API impact

The preliminary expectation is:

- **Python API:** low likelihood of broad changes. Owned conversion helpers, shared internals, and
  borrowed operation bridges can preserve current call shapes. Parent-backed return values or
  consuming Python arguments would be visible changes and require specific justification.
- **Rust API:** moderate likelihood of targeted ownership changes, particularly for chemistry-model
  data and operation inputs. These changes may be additive or may replace exposed `Cow<'static, _>`
  fields if that representation obstructs efficient Python reuse.
- **Binding internals:** high likelihood of cleanup. Owned and borrowed conversions should be
  explicit, duplicate reconstruction should be removed, and detached versus parent-backed behavior
  should be documented consistently.

## Review requirements

The eventual review should produce an inventory grouped by user operation, not a raw list of types or
`.clone()` sites. It should include at least:

- construction and repeated use of a custom `ChemistryModel` with a large registry or table;
- molecule and reaction operations repeatedly using the same model;
- small and large `Edits` and `Deltas` crossings;
- indexing and iteration over edit and delta containers;
- large molecule, reaction, correspondence, and fingerprint results;
- field access through molecule-backed views;
- owned Rust results whose conversion currently starts from a borrow.

Benchmarks belong at the beginning of the work. They must separate Python object-allocation cost,
Rust data-copy cost, and the cost of the underlying operation where practical. The number of cloned
types is not a useful optimization target by itself.

## Sequencing

Finish the implementation work tracked by doc 179 before starting this review. The boundary audit may
use the completed editing and transaction API as one of its representative workloads, but it must not
expand or interrupt that work list.
