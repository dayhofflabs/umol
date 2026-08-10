# 192 — Python API type roles

Status: Proposed
Date: 2026-08-09
Relates: [139](139-mutability-hashability-equality-2026-07-09.md),
[176](176-ast-naming-2026-07-31.md), [178](178-python-lattice-ops-2026-08-01.md),
[179](179-python-editing-and-transactions-2026-08-02.md),
[181](181-python-boundary-ownership-2026-08-03.md)

## Purpose

Review every public `umol-py` type, assign one observable Python ownership and mutation role, and
split types whose current behavior combines incompatible roles. The normative rules live in
[`docs/development/python-api.md`](../docs/development/python-api.md); this document holds the
inventory, design decisions, open questions, and eventual implementation plan.

The public Rust API is the baseline. The Python binding should preserve its concepts and operation
semantics where reasonable, omit surfaces that are not useful in Python, and deviate only for a
specific Python usability or correctness reason. This is not a greenfield redesign of the Rust type
hierarchy. In particular, a property must not return a disconnected mutable copy: assignment either
writes through to the object from which the value was obtained or is rejected.

## Current public surface

`umol-py/python/umol/__init__.py` currently lists 215 public names:

- 206 classes registered with `module.add_class`;
- 7 exception classes registered with `module.add`;
- the element namespace `E`; and
- `__version__`.

The public surface therefore contains 213 class-like types. The 232 `#[pyclass]` declarations in
`umol-py/src` are not the public count because they include private iterators and adapters.

The following suffix counts overlap but show the scale and repeated structure of the review:

| family | count |
| --- | ---: |
| `*Form` | 42 |
| `*Update` | 10 |
| `*View` | 17 |
| `*Views` | 8 |
| `*Config` | 13 |
| `*Model` | 5 |
| `*Delta` | 10 |
| `Deltas` | 1 |
| `*FieldChange` | 8 |
| `*Algorithm` | 9 |
| exception classes | 7 |

The inventory must use the package's `__all__` as the public boundary. Counting Rust wrappers alone
would include implementation classes and omit exception classes.

## Whitepaper yardstick

The Python examples in the whitepaper are the current target workflows. They do not completely
specify the API, but a role assignment that makes them awkward is suspect.

### Persistent molecule and reaction values

The primer parses molecules from SMILES and the umol notation, renders them, preserves metadata
through an explicit parse/render pair, and builds reactions from reaction SMILES. Operations are
value-producing:

- `combine` returns a combined molecule and a correspondence;
- `split` returns molecule/correspondence pairs;
- reaction application returns derivations whose `rhs` values are molecules; and
- reaction composition returns composite reactions.

No example requires direct mutation of a molecule or reaction. The ordinary mutation path applies
an edit batch and returns a new molecule. This establishes the important value-producing operations,
but absence from the examples does not prove that the complete `Molecule` or `Reaction` surface must
be frozen. Their Rust APIs and other intended Python workflows remain part of the decision.

### Declarative forms and lattice values

Entity forms, updates, and lattice members are parsed or constructed as values. `meet`, `join`,
`matches`, and `is_compatible` do not mutate their operands. The whitepaper requires these values to
be easy to construct and pass to operations; it does not require field assignment after
construction. Whether the concrete Python forms also expose mutation should follow the corresponding
Rust API unless doing so creates incompatible nested-value behavior.

### Models and configurations

Chemistry models, resolution configuration, matching configuration, algorithm selectors, and
fingerprint configuration are constructed and passed to operations. The examples do not mutate them
after construction. Their natural initial role is an immutable owned value unless a stateful model
facility establishes a different need.

### Results

Correspondences, derivations, metadata, feature sets, and fingerprints are queried, iterated,
composed, folded, or passed to later operations. The examples do not mutate them. Correspondence
properties such as `matched_pairs` must therefore be safe value observations rather than mutable
aliases into correspondence storage.

### Edits

The mutation section deliberately uses `Edits` as a mutable ordered builder:

```python
edits = Edits()
carbon = edits.add_atom(AtomForm.parse("C#h3"))
edits.add_bond(0, carbon, BondForm.parse("1"))
edits.update_atom(0, AtomForm.parse("N#h3"), AtomUpdate.parse("#h2"))

methylamine = ammonia.apply(edits)
```

The exact type names in the whitepaper still reflect the aggregate and entity naming migration, but
the lifecycle is the point: `Edits` mutates in place, additions return stable handles, and applying
the reified batch returns a new molecule. `Edits` must also remain inspectable, parseable, renderable,
and replayable; it is not merely an opaque editor.

### Views and editors

Entity views and `MoleculeEditor` do not appear in the ordinary whitepaper workflow. They may support
advanced mutation without changing the primary value API, but their existence must not force
value-like forms or molecules to acquire ambiguous mutation behavior. The whitepaper's ordinary path
does not require `snapshot` or direct editor use.

## Initial role assignments

These assignments identify likely directions to verify against the corresponding Rust APIs. They
are not a substitute for reading each complete Rust and Python surface.

| role | initial families | required behavior |
| --- | --- | --- |
| immutable owned value | individual `Edit` and `Delta` variants; other Rust value types without a public mutation surface, to be established by the inventory | no assignment or mutating methods; nested umol values are also immutable |
| mutable owned object | `Edits`, `Deltas`, and other types when ordinary mutation is established | mutation changes that object; nested mutable access writes through |
| live view | `*View` and `*Views` families where they access an owner | reads track the owner; supported writes update the owner; owner lifetime is retained |
| editor | `MoleculeEditor` | staged mutation is explicit and finalization produces an owned value |

`Transaction` is an operation-issued mutable object with a consuming lifecycle rather than an
editor. Stateful and one-shot lifecycles are review properties within these roles, not additional
naming families. Value equality and copying are usually inappropriate for such resources.

## Known mixed-role areas

### Initial Rust findings

`AtomForm` is not an immutable Rust value hidden behind a mutable Python wrapper. Its fields are
public, its constraint container is directly mutable, and its `with_*` methods provide an additional
owned construction style. The existing Python field setters therefore parallel the Rust surface.
The problem is narrower: access to an `AtomForm` nested in another object must say whether it is a
live mutable value, a read-only observation, or an explicit copy.

Rust likewise has one `Molecule` type. It exposes mutable per-entity and molecule-constraint access
through `atom_mut`, the other `*_mut` entity accessors, and `constraints_mut`, while structural
changes go through `MoleculeEditor`. Python retains the same division: one `Molecule`, direct
attribute and constraint mutation through the existing live entity views, and `MoleculeEditor` for
topology and overlay mutation. Attribute mutation must not require an editor. A parallel
`Molecule`/`MutableMolecule` pair is not needed for this boundary.

Individual `Edit` and `Delta` variants are immutable in Python. Their containers are independently
mutable: `Edits` and `Deltas` are constructible from iterables and support `append` and `extend`,
preserving order and duplicates. This is list-like accumulation, not a commitment to item
assignment, insertion, deletion, sorting, or reordering. Extension appends raw entries without
rebasing identifiers or `New(n)` handles.

### Entity forms

Entity forms currently expose setters while also serving as nested values inside deltas and other
declarative objects. Returning a cloned mutable form from a delta property produces misleading code:

```python
atom = delta.attributes
atom.charge = 1
```

If the assignment does not modify the delta, it must be rejected. The review must decide whether
the existing form can preserve the Rust semantics consistently in every Python context or whether a
split is unavoidable. It must not assume a parallel mutable/frozen hierarchy in advance. If a split
is needed, `T` / `TMut` is an available pattern when mutation capability is exactly the distinction
and both wrappers have clear consumers. The decision is local; it does not imply generating a
mutable/frozen pair for every public type.

### Constraint containers

The existing container design is accepted and its distinction is behavioral:

- value-backed `*ConstraintsForm` classes are owned mutable containers;
- an entity view exposes constraints whose writes change its owning molecule; and
- explicit copying may produce an independent mutable container when the operation says so.

The remaining question belongs to a frozen parent that contains constraints: its property must not
return a mutable disconnected copy that appears to update the parent. Resolve that access boundary
without replacing the general constraint-container design.

### Updates and deltas

Individual delta wrappers are immutable, but several nested form or constraint properties remain
mutable. PyO3's frozen marker prevents mutable Rust borrowing of the wrapper; it does not make the
Python-visible object graph immutable. The remaining delta problem is therefore deep immutable
access to payloads, not whether the variant itself can be assigned to.

`Deltas` already provides the selected constructor, `append`, and `extend` shape. `Edits` provides
iterable construction and `append` but still needs `extend`. Extension has ordinary append-only list
semantics: it preserves entry order and every handle exactly as written. It does not rebase an
independently constructed sequence. A `New(0)` in the appended entries remains `New(0)` and is
interpreted in the resulting sequence.

This matches the whitepaper construction example. `add_atom` issues a `New` handle and the later
`add_bond` entry stores and reuses that same handle. Handles are part of the immutable edit entry;
container operations do not reinterpret or rewrite them.

### Molecules and entity views

The molecule boundary is settled. Functional operations such as `combine`, `split`, and `apply`
coexist with direct attribute and constraint mutation through the existing live entity views.
Topology and overlay mutation goes through `MoleculeEditor`. The views are not candidates for
replacement merely because the whitepaper does not exercise direct attribute assignment.

## Review matrix

Every public class-like name receives one row with the following fields:

| field | question |
| --- | --- |
| public name | What name appears in `umol.__all__` after the doc 176 migrations? |
| family | Is it an entity form, update, delta, view, config, model, result, container, error, or selector? |
| Rust surface | Which public Rust type, constructors, accessors, transformations, and mutating operations define the baseline? |
| current construction | Which constructors, variant constructors, parsers, and conversions create it? |
| current mutation | Which setters, deletion operations, mutating methods, or one-shot transitions exist? |
| nested returns | Which properties return umol objects, and are they copies, immutable values, or live views? |
| current aliasing | Can a returned object modify its owner, or can a retained input modify the constructed object? |
| whitepaper use | Is it named, returned, or implied by a target workflow? |
| target role | Immutable value, mutable owner, live view, or editor? |
| lifecycle | Is it reusable, append-only, staged, consuming, or operation-issued? |
| split required | Does one current wrapper cover more than one target role? |
| Python deviation | Which concrete Python constraint justifies behavior different from the Rust API? |
| equality and hash | Are equality and hashing stable under the target mutation behavior? |
| conversion cost | Does enforcing the role require a clone, move, immutable projection, or back-reference? |
| verification | Which construction-isolation, mutation-rejection, write-through, equality, and hash tests express the contract? |

The matrix records public semantics first. Private wrapper consolidation and copy optimization follow
from the selected roles; they do not decide them.

## Review order

1. Classify the whitepaper-critical path: molecule, reaction, forms, updates, `Edits`,
   correspondences, derivations, configs/models, metadata, and fingerprint values.
2. Apply the settled entity-form decision uniformly across all eight entity families.
3. Verify the accepted constraint-container semantics across entity-local and molecule constraints,
   concentrating changes on frozen-parent access boundaries.
4. Review deltas, field changes, individual edits, and their containers recursively.
5. Review live views, view namespaces, `MoleculeEditor`, and `Transaction` against their owners and
   lifecycle guarantees.
6. Review remaining selectors, models, configs, registries, tables, errors, iterators, and utility
   values.
7. Convert the settled matrix into a staged implementation plan that keeps each family internally
   consistent at stage boundaries.

## Open questions

- Can each entity form preserve its Rust mutation semantics in Python without a second public form
  type, especially when nested inside immutable deltas and updates?
- Which frozen-parent properties need a read-only form wrapper, potentially using a local `T` /
  `TMut` pair, to expose nested forms and constraints without misleading copy semantics?
- Should metadata remain immutable values with functional `remap`, or is there a supported metadata
  editing workflow?
- Which model and registry types are values, and which represent intentionally stateful reusable
  facilities?
- Are any current hash implementations invalidated by setters or live nested values?
- Which Python properties currently return mutable disconnected clones, beyond the known delta and
  update cases?

## External precedents

Rust backing does not determine Python mutability. [Polars](https://docs.pola.rs/) uses predominantly
value-producing frame operations, while [rustworkx](https://www.rustworkx.org/install.html) exposes
directly mutable graphs. [RDKit](https://www.rdkit.org/docs/GettingStartedInPython.html), a
domain-relevant C++ precedent, uses a separate `RWMol` for structural editing. Python's concrete
mutable built-ins are simply `list`, `dict`, and `set`; the `Mutable*` spelling in the standard
[container ABCs](https://docs.python.org/3/library/collections.abc.html) names protocols rather than
establishing a general concrete-type convention. These precedents support preserving the underlying
domain API and making only necessary Python adaptations, not manufacturing a parallel mutable and
frozen hierarchy.
