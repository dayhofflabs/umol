# Python API

## Purpose

This guide defines the public ownership and mutation contract of the Python API. The Rust storage
strategy and PyO3 implementation are not themselves part of that contract. A Python user should be
able to tell from a type and operation whether a value is immutable, independently mutable, or a
live accessor into another object.

The public Rust API is the starting point. Its types, constructors, value-producing operations,
mutating operations, views, editors, and lifecycle boundaries carry deliberate semantics that the
Python binding preserves where Python can express them reasonably. The binding may omit Rust APIs
and may adapt call syntax, but it does not redesign the type hierarchy merely because Python has a
different ownership model. A Python-specific deviation needs a concrete usability or correctness
reason.

The central rule is that apparent mutation must have one observable meaning. An assignment either
changes the object through which the value was obtained or is rejected. A property must not return
a disconnected mutable copy on which assignments succeed but have no effect on the parent.

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

Every exported class has one of the following roles.

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

An immutable entry type may belong to a mutable container. `Edit` and `Delta` variants are immutable
values, while `Edits` and `Deltas` provide iterable construction, `append`, and `extend` in insertion
order. Container mutation does not imply entry mutation. These operations establish list-like
append-only accumulation, not unrestricted list replacement, insertion, deletion, or reordering.
Extension appends entries exactly as written and does not rebase or otherwise rewrite identifiers or
`New(n)` handles.

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

Entity forms are ordinarily writable owned objects. When an entity form is stored in an immutable
delta, construction clones its value into a form instance whose `readonly` property is true. The
state cannot be changed in place, getters retain that same read-only instance, and field and
constraint mutation raise `TypeError`. `copy`, `normalize`, `meet`, and `join` produce ordinary
writable forms. Individual delta variants are immutable; `Deltas` provides only append-only
container mutation.

Use `*Like` argument adapters for accepted alternate input representations. A `*Like` type is an
argument boundary and is not exposed as the stored or returned type.

Rust/Python boundary conversions use `from_rust`, `to_rust`, and `to_rust_mut` where mutable Rust
access is genuinely part of the wrapper's role. Private PyO3 adapter types are named for the
representation behavior they implement, not for one consumer such as a delta. They must not create
a second public data model.

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
