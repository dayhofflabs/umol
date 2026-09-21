# Python API

## Purpose

This guide defines the public ownership and mutation contract of the Python API. The Rust storage
strategy and PyO3 implementation are not themselves part of that contract. A Python user should be
able to tell from a type and operation whether a value is immutable, independently mutable, or a
live accessor into another object.

The public Rust API defines the operations and contracts to preserve. Every difference in Python
requires a specific justification, without exceptions: exposed types and methods, names, semantics,
construction, lifecycle, ownership, aliasing, failure behavior, and copying. The binding may omit
Rust APIs or adapt call syntax for a concrete usability reason; it must not invent domain
operations. Existing code, tests, documentation, and prior acceptance do not justify a deviation.

For each deviation, identify the public Rust counterpart, the precise difference, the concrete
Python requirement, and the evidence supporting the chosen adaptation. Claims that PyO3 or the
Python ABI requires a design need evidence for the relevant version and code path. An unexplained
or unreviewed difference remains unresolved; familiarity and passing tests do not establish parity.

Cloning and snapshots require the same justification. State what is copied and why, distinguishing
shared ownership from copying owned data and subsequent copy-on-write costs. Examine whether
borrowing, ownership transfer with an explicit consumed state, or a supported live view can preserve
the Rust contract before selecting a copy. These are alternatives to evaluate, not mandatory new
wrappers. Converting a consuming Rust operation into a reusable copy-producing Python operation is
a lifecycle change. Neither defensive cloning nor a universal ban on copying decides its contract.
Do not change the Rust API solely to legitimize a wrapper convenience.

Preserve consumption when the owning Rust operation consumes its input. Possible reuse does not
justify an unconditional copy: callers explicitly copy when they need to retain the input, as in
Rust. A departure requires a concrete justification beyond possible reuse or a general claim that
consumption is unsuitable for Python. Failure ordering and accessor implementation are interface
questions, not reasons to assume copying is necessary.

Consumption leaves all aliases of the Python owner in a consumed state and invalidates its
outstanding views and iterators, including nested accessors. Access raises one of two shared
exceptions, both subclasses of `RuntimeError`:

- `ConsumedError`: access to an owned object whose Rust contents have been consumed.
- `InvalidatedViewError`: access through a view or iterator invalidated by owner consumption or
  structural replacement.

Messages identify the consumed type or the invalidated accessor and its owning type. Explicit
independent copies made before consumption remain usable. Do not copy automatically on invalid
access or block consumption solely because an accessor remains alive.

The central rule is that apparent mutation must have one observable meaning. An assignment either
changes the object through which the value was obtained or is rejected. A property must not return
a disconnected mutable copy on which assignments succeed but have no effect on the parent.
This includes nested containers and mutating methods. Mutable access updates backing storage;
immutable observations reject mutation. An explicitly requested independent copy may be mutable,
with mutations affecting its own storage. Calling a property result detached does not justify
silent, ineffective writes.

Keep bindings thin and cheap. Prefer read-only access to backing Rust storage over copied
observations where possible. Preserve iterator laziness: creating an iterator must not materialize
its entries, and requesting one entry must not convert the unrequested remainder. Prefer a thin
read-only Python accessor for each yielded entry, including nested access, over copying its payload.
Owner lifetime and interaction with mutation or consumption need explicit contracts; they do not
justify defaulting to snapshots. Python being the primary application interface does not justify
separate semantics or additional copying.

## Graph-IR naming at the boundary

The Python surface follows the graph IR's public concepts without retaining obsolete Rust-layer
names. Lattice types use the `*Form` suffix, while the non-lattice aggregate roots are `Molecule`,
`Reaction`, and `ReactionSpan`. No `*Ast` compatibility classes are exported.

Entity forms are the `attributes` payload of entity views, edits, and deltas. Constructors,
properties, annotations, representations, and structural pattern-matching fields use that name.
The corresponding DSL map key is `:attrs`. Rust-side binding conversions use `from_rust`, `to_rust`,
and `to_rust_mut`; graph-IR boundary conversions use the `*Ir` trait family.

Recursive subpattern constraints are currently absent from both Rust and Python. They must not be
reintroduced on only one side of the binding boundary.

## Public object roles

Every exported class has an explicit ownership/access role. Combining owned and accessor
representations in one class requires a type-specific justification and an explicit mutation
contract.

### Reusable owner-backed access pattern

Where Rust exposes borrowed access to stored data, a Python accessor can retain its owner and a
location within that owner. Each operation briefly borrows the owner, checks availability, locates
the data, and performs the requested access. No Rust reference into relocatable storage survives
between Python calls. Nested accessors preserve the same ownership and invalidation contract.

Read-only access rejects mutation; mutable access updates backing storage. Explicit copying
produces an independent owned value. Consumption invalidates dependent accessors as specified
above. The location and its behavior under structural mutation must be defined for the owning type.

For the Edit/Delta nested-access pattern, structural replacement invalidates all existing child
accessors, including container accessors and their descendants. Fresh accessors expose the new
state. Ordinary mutation within an existing value remains visible through valid accessors.
Replacement does not preserve selected children, retarget them to replacement values, or copy
their former contents. Independent explicit copies remain usable.

This is a reusable pattern, not a requirement to give every type an owned-or-accessor representation.
Owned-only values and accessor-only views remain appropriate. Where one Python class legitimately
supports both roles, an internal owned-or-owner-and-location representation is possible; it is not
ordinary Rust `Cow` and must not clone automatically on mutation. Public role distinctions are
settled per type. Introduce shared implementation machinery only when concrete uses support it;
this pattern does not authorize a universal wrapper framework or a blanket binding rewrite.

For Edit/Delta entries, owned values and read-only accessors use the same Python variant classes
and field names. This preserves one variant-inspection and pattern-matching interface instead of
duplicating the variant hierarchy. A read-only `readonly` property reports mutation permissions;
it is not a general ownership indicator. Accessors obtained from batches are read-only, including
nested access. `copy()` returns an independent owned entry of the same variant. The backing-storage
distinction remains internal; mutation never silently detaches an accessor by copying.

### Immutable owned value

An immutable owned value represents data rather than a mutable identity. Its public fields cannot
be assigned or deleted. Operations that change it return another value or an explicitly mutable
counterpart.

Immutability applies transitively to the umol-owned objects exposed through its properties. An
immutable `AtomForm` must not expose a mutable constraints container, for example. Returning an
immutable value copy is acceptable because the returned object cannot imply write-through
behavior.

Small forms, enum variants, field changes, deltas, and other declarative values are candidates for
this role, but being data-like does not establish immutability. Preserve a Rust type's deliberate
public fields and mutating operations unless a Python-specific conflict requires adaptation. Do not
create mutable counterparts for immutable leaf values merely for symmetry.

### Mutable owned object

A mutable owned object has independent identity and state. Assignments and mutating methods change
that object. A value retrieved from one of its properties either writes through to the object or is
itself immutable.

Use this role when mutation is an ordinary, long-lived way to work with the object. Do not make a
type mutable solely because its Rust representation is mutable internally.

### Live view

A live view is an accessor into an owning object. It retains the owner for at least as long as the
view can be used. Reads observe the owner's current state and supported writes update the owner.
The `*View` suffix is reserved for this role; it does not mean a copied value.

An `AtomView` obtained from a molecule is therefore distinct from an owned `AtomForm`. The view may
write through to the molecule, while the form has its own ownership and mutation contract.

### Editor

An editor represents a staged or transactional mutation lifecycle. It is appropriate when users
begin with an owned value, perform several related structural changes, and then finalize or commit
a new value. The `*Editor` suffix communicates that lifecycle and is preferable to presenting the
editor as a generally mutable form of the original type.

For molecules, the editor is the structural mutation boundary: adding or removing topology and
overlays goes through `MoleculeEditor`. Ordinary attribute and constraint mutation does not require
an editor; it writes through the molecule's entity and constraint views.

### Operation-issued iterator

An operation-issued iterator captures one configured execution and cannot be constructed
independently. Reaction application is the reference case. `Reaction.apply` performs the
reaction-wide precondition check and returns a one-shot iterator over product molecules;
`tracked_apply` yields `(product, correspondence)` tuples, while `apply_to_reaction` and
`apply_to_reaction_span` yield realized reactions and spans.
`Molecule.react` and `Molecule.react_all` return a one-shot iterator over product-component lists.
These iterator classes are not exported or directly constructible.

The iterator owns snapshots of every input needed by the operation, so mutating the source reaction
or molecules after the call does not change the pending results. Matching is completed when the
iterator is issued; products, realized reactions/spans, and product splitting are realized lazily
in match order.
An eager reaction-wide failure is raised by the operation call. An execution failure while
realizing a selected match is raised by `next` once and terminates the iterator. This placement
parallels the outer `Result` and fallible iterator item on the Rust side.

Use `tracked_apply` when product provenance is required; use `react` when only disconnected
product molecules are needed:

```python
for product, correspondence in reaction.tracked_apply(host):
    print(correspondence)

for products in Molecule.react_all([first, second], reaction):
    print(products)
```

## Choosing the Python mutation model

Begin with the public behavior of the corresponding Rust type. An owned transformation normally
remains value-producing in Python; a Rust operation that deliberately mutates through `&mut self`
normally mutates the same Python object; a Rust view remains a view; and an editor remains an
editor. Internal Rust mutability used only to implement an operation does not make the Python object
mutable.

Rust attaches mutability to bindings and references, so one Rust type may be observed through both
`&T` and `&mut T`. Python normally exposes one method and property surface on a concrete object.
This mismatch may require a live view or another explicit role at a particular access boundary, but
it does not justify duplicating the entire type family into mutable and frozen classes.

Expose one concrete Python type per Rust concept by default. Do not generate a complete mutable and
frozen pair for every value type. Python's concrete mutable built-ins use unqualified names such as
`list`, `dict`, and `set`; `Mutable*` is primarily the spelling of container ABCs, not a general rule
for concrete classes. Likewise, do not export a Rust-style `*Mut` suffix merely to encode how a
reference is borrowed.

Add a qualified counterpart only when there are genuinely two public concepts. Preserve an
established Rust role name such as `*Editor` or `*View` when it describes the Python behavior. If a
long-lived mutable and immutable pair is unavoidable, settle its public names from the domain and
Python usage rather than applying a repository-wide suffix rule.

`T` / `TMut` remains available for a local pair whose only meaningful distinction is mutation
capability, especially when it directly parallels a carefully designed Rust pair. It is neither
mandatory nor prohibited. Use it only where both types have a clear public consumer and one wrapper
cannot express both contracts without ambiguous writes.

The short, unqualified name belongs to the ordinary domain type. The existence of a mutating method
does not by itself require renaming that type: Python users ordinarily expect a concrete object to
advertise mutation through its methods and properties. A split is justified when one name would
otherwise cover observably incompatible aliasing or lifecycle behavior.

An immutable entry type may belong to a mutable container; container mutation does not imply entry
mutation. Python collection idioms must preserve the owning Rust container's reference semantics.
A spelling such as `append` may delegate to Rust `push`; it does not authorize batch concatenation.
Entries in an `Edits` sequence share one handle namespace. Combining independent batches requires
a separately designed Rust operation on `Edits`; raw list extension cannot establish that contract.
Any supported bulk operation on `Deltas` likewise belongs to its Rust container. Do not implement
collection composition in the bindings merely because Python lists support `extend`.

## Properties and nested objects

A property contract is determined by the returned object's role:

- an immutable property may return the same immutable object or an immutable copy;
- a property of a mutable owner may return a live mutable view when writes are meant to propagate;
- a property must not return a mutable disconnected copy;
- an explicit `copy`, `to_mutable`, or equivalent operation may return an independently mutable
  object when that distinction is useful; and
- returning a live view solely to avoid a copy is not sufficient reason to introduce aliasing.

These rules apply recursively. Marking the outer PyO3 class as frozen is insufficient when a nested
constraints container or other umol object remains mutable.

The existing constraint-container distinction is retained. A value-backed `*ConstraintsForm` is an
owned mutable container; a molecule-backed `*ConstraintsView` is a live container whose writes
update the molecule. Entity views use the latter and do not require an editor. If a frozen parent
needs to expose constraints, resolve that parent's access contract without replacing this general
container design or returning a misleading mutable disconnected copy.

## Construction and conversion

Construction of an immutable owner must not retain a mutable Python alias that can later change the
owner. Convert mutable inputs to the owner's immutable representation at the construction boundary.
This is an ownership rule, not semantic validation; validation follows the data-type contracts.

Entity forms are ordinarily writable owned objects. Access through a read-only entry is read-only,
including nested fields and constraints. Such access rejects mutation with `TypeError`.
`copy`, `normalize`, `meet`, and `join` produce ordinary
writable forms. `Deltas` provides append-only container mutation; its entry accessors do not permit
mutation of stored entries.

Use `*Like` argument adapters for accepted alternate input representations. A `*Like` type is an
argument boundary and is not exposed as the stored or returned type.

Rust/Python boundary conversions use `from_rust`, `to_rust`, and `to_rust_mut` where mutable Rust
access is genuinely part of the wrapper's role. Private PyO3 adapter types are named for the
representation behavior they implement, not for one consumer such as a delta. They must not create
a second public data model.

## Aggregate canonicalization

`Molecule`, `Reaction`, and `ReactionSpan` expose complete-only `canonicalize`,
`tracked_canonicalize`, and `canonical_eq`. Each accepts the same optional
`StereoModel` and `CanonicalizeConfig`; canonicalization of an intrinsically contradictory value
raises `ContradictionError`. Python does not expose `DescriptionLevel`, a molecule
`description_level` query, or level-parameterized `*_by` operations. Topology, constitution, and
structure remain private search prefixes rather than reduced public comparison surfaces.

`tracked_canonicalize` returns exactly the same canonical aggregate as `canonicalize`
plus a `MoleculeRemapping` containing the source-to-canonical entity-id permutations. The
participant-frame action is the separate witness consumed internally by Rust reframing. Python does not expose frame-action classes
until a supported Python operation has an independent action consumer; the remapping does not
implicitly promise participant-frame transport by itself.

Canonical representatives may change between umol 0.x releases and are not persistent identifiers.
Persist the molecular assertion, not a canonical entity numbering, unless a future API supplies an
explicitly versioned canonicalization profile.

## Equality and hashing

Immutable values may be hashable when their Rust semantics define stable equality and hashing.
Mutable objects are not value-hashable. Views use the equality semantics of the value they expose
only when those semantics remain stable and unsurprising while the owner changes; otherwise they
must not imply immutable-value behavior.

Equality across immutable and mutable representations is a deliberate API decision. Do not obtain
it incidentally from conversion or wrapper internals.

## PyO3 implementation notes

`#[pyclass(frozen)]` controls mutable Rust borrowing through PyO3. It does not by itself make the
Python-visible object graph immutable. The public contract also depends on setters, mutating
methods, nested return values, aliases, and views.

Cloning a Rust value may be a correct implementation of an immutable return. It is not a semantic
justification for returning a mutable Python object that looks connected to its source. Copy and
alias behavior must be chosen first; cloning follows from that choice.

## Contract tests

Tests for an exported aggregate should cover the behavior that its role promises:

- assignment to immutable fields is rejected;
- nested umol values obtained from an immutable owner are also immutable;
- mutating an independently mutable object changes that object;
- mutating a live view changes its owner;
- constructing an immutable owner from a mutable input does not retain a mutation alias;
- an explicit mutable copy does not change its source; and
- hashing is unavailable for mutable objects and stable for hashable immutable values.

Tests should assert these positive contracts. They do not need to enumerate every internal helper
or unexported adapter.
