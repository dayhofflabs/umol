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

Python aliases do not follow Rust's move rules. Molecule remains reusable: edit, apply,
tracked_apply, resolve, resolve_with_report, canonicalize, tracked_canonicalize, combine, and the
base transformation methods below clone the Rust receiver and pass that clone to the consuming
Rust operation. The Python source, its aliases, and its views remain unchanged on every outcome.
The initial clone shares topology, entity tables, and constraints; shared payloads are copied
when the operation modifies or takes ownership of them. This is the same ownership choice as a
Rust caller explicitly cloning its input to keep it available. Python exposes named transformation
methods on Molecule without Transformer engine classes. The project operation remains Rust-only.

Edits retains consumption on transfer: implicitly cloning a complete batch would copy its entries
and payloads, so callers request Edits.copy() when they need reuse. MoleculeEditor.finish consumes
the editor on success or failure. An into_… conversion retains its consuming lifecycle;
RingViews.into_ring_set consumes the views while leaving the molecule usable. Deltas has no
consuming Python operation; its explicit copy and value-producing operations preserve the source.

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

Named Rust variant fields retain their names in Python properties, constructor keywords,
representations, and `__match_args__`. Genuine tuple fields use `_0`, `_1`, and so on.
All payloads are accessible through properties on Python 3.9. Structural pattern matching is
an optional convenience on Python 3.10 and later, never the only way to inspect a value.

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

For molecule access, each view retains the owner and its accessor counter. Starting transact,
tracked_transact, resolve_into, resolve_into_with_report, combine_from, or any transformation's
_into method invalidates all existing views and iterators, including nested accessors. The counter
advances immediately before execution, including for a no-op, rejection, or rollback. No list of
views is maintained. Ordinary attribute and entity-constraint assignments leave accessors valid
and their changes visible. Fresh accessors read the current molecule; independent copies remain usable.
Argument extraction or borrowing failure before execution leaves the counter unchanged.

RingViews retains its molecule and any borrowed RingSet. into_ring_set moves owned ring data or
copies borrowed ring data into an independent RingSet without repeating perception. It consumes
the RingViews, whose later access raises ConsumedError; its dependent ring, atom, bond, and
iterator accessors raise InvalidatedViewError. The molecule and any source RingSet remain usable.

This is a reusable pattern, not a requirement to give every type an owned-or-accessor representation.
Owned-only values and accessor-only views remain appropriate. Where one Python class legitimately
supports both roles, an internal owned-or-owner-and-location representation is possible; it is not
ordinary Rust `Cow` and must not clone automatically on mutation. Public role distinctions are
settled per type. Introduce shared implementation machinery only when concrete uses support it;
this pattern does not authorize a universal wrapper framework or a blanket binding rewrite.

Edit and Delta are immutable owned values. Indexing Edits or Deltas returns an independent entry;
iteration copies one entry per step and excludes entries appended after iterator creation.
Entity forms retained within an entry are read-only, including nested constraints; their readonly
property reports that restriction, and copy() returns an ordinary writable form. These copies
isolate immutable mutation data from later changes to the supplied mutable forms without requiring
consumed state on every form. Consuming Edits invalidates its iterators; already returned entries
remain usable.

Edits.copy() explicitly clones the complete Rust batch into an independent, unconsumed owner.
It preserves entry order, duplicates, handles, payloads, and creation ordinals without checking
applicability against a molecule. Appending to or consuming either batch leaves the other usable;
copying leaves existing source iterators valid. A consumed source raises ConsumedError.

Deltas.copy() clones the complete Rust collection into an independent owner, preserving order,
duplicates, ids, and payloads without normalization or applicability checks. Both containers and
existing source iterators remain usable; appending to either leaves the other unchanged.

Delta.inverse() and each entity delta's inverse() return an independent immutable inverse,
leaving the source usable. Retained entity forms remain read-only. Inversion preserves the Rust
operation: it does not normalize payloads or check applicability against a molecule.

### Immutable owned value

An immutable owned value represents data rather than a mutable identity. Its public fields cannot
be assigned or deleted. Operations that change it return another value or an explicitly mutable
counterpart.

Collection payload properties use tuples for sequences and frozensets for sets, including
nested restriction sets and participant rows. StereoModel.kind_models returns a mapping proxy
over an independent dictionary of the existing model values. Constructors accept their ordinary
list, set, and dictionary inputs as well as the corresponding immutable getter results.
These observation containers preserve the stored contents, order, duplicates, and absence;
they do not normalize the underlying values.

Immutability applies transitively to the umol-owned objects exposed through its properties. An
immutable `AtomForm` must not expose a mutable constraints container, for example. Returning an
immutable value copy is acceptable because the returned object cannot imply write-through
behavior.

Small forms, enum variants, field changes, deltas, and other declarative values are candidates for
this role, but being data-like does not establish immutability. Preserve a Rust type's deliberate
public fields and mutating operations unless a Python-specific conflict requires adaptation. Do not
create mutable counterparts for immutable leaf values merely for symmetry.

Each Python ChemistryModel owns a complete Rust ChemistryModel. Construction converts the supplied
connectivity, valence, aromaticity, and stereo components once; custom registries and tables may be
copied at that boundary. ChemistryModel.default() preserves the Rust default's borrowed static
registry. Operations borrow the retained model for the synchronous call, without consuming the
Python object or rebuilding its registry or table. Omitting the model selects the operation's Rust
default, which can differ from ChemistryModel.default().

The four component properties return independent values of their existing Python types. Reading a
component may copy its registry, table, or other owned data. Reusing the ChemistryModel directly
avoids that extraction and reconstruction cost; rebuilding it from its properties does not.

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

An editor owns transient working state. It permits direct changes and edit batches before finish
checks integrity and returns a molecule. It provides no rollback.

Structural changes use MoleculeEditor or prepared transaction batches. Ordinary attributes and
entity-level constraints are assigned through the molecule's entity views. Molecule-level
constraints are read-only on the molecule; changing them requires editing and publication.

| Python operation | Result | Receiver and batch ownership |
| --- | --- | --- |
| Molecule.edit | MoleculeEditor | Edits an independent copy; source and views remain usable |
| Molecule.apply | Molecule | Consumes Edits and publishes a changed copy; source and views remain unchanged |
| Molecule.tracked_apply | (Molecule, MoleculeCorrespondence) | Same ownership as apply |
| MoleculeEditor.apply | None | Mutates the editor and consumes Edits |
| MoleculeEditor.tracked_apply | MoleculeCorrespondence | Same ownership as editor apply; maps only this batch's input to its result |
| MoleculeEditor.finish | Molecule | Consumes the editor on success or failure; checks representation integrity |

Editor application checks edit preconditions but does not check final molecule integrity or record
undo. Failure raises TransactionError and retains the completed prefix for further editing; it
returns no correspondence. finish raises InvalidStructureError if publication fails, leaving every
alias of the editor consumed. These are the Rust editor's mutation and publication contracts.
Molecule.apply and tracked_apply include finish; either kind of failure discards the working copy
and leaves the source unchanged.

Python omits Rust Molecule::mutate, the direct batch primitive for producers that establish
successful application in advance. Python's general batch operations accept independently
assembled Edits and report ordinary exceptions; mutate instead panics and empties the receiver
during unwinding. That failure contract is not exposed as a second general Python batch API.
High-level operations may use mutate inside their Rust implementation while preserving their
own checked Python outcomes. This omission is a usability choice, not a PyO3 limitation.

MoleculeEditor.remove_topology(atoms, bonds) removes exactly the supplied integer IDs and returns
None. It raises TransactionError before changing storage if a stored participant or site is
unavailable (InvalidReference), or a retained entity would reference a removed entity
(DanglingReference). The editor remains usable after rejection. ID bounds errors remain IndexError.
MoleculeEditor.remove_topology_cascading(atoms, bonds) explicitly includes incident bonds and
overlays with unavailable stored references, including references already absent in an editor
draft. It returns None and compacts molecule constraints with the completed removal.

Each of remove_topology, remove_topology_cascading, remove_dative_bonds,
remove_aromatic_systems, remove_multicenter_bonds, remove_noncovalent_bonds,
remove_stereo_atoms, and remove_stereo_bonds has a tracked_ counterpart returning
MoleculeCompaction. Tracking preserves the corresponding mutation and failure behavior. Every
component records its source count immediately before the call, including identities for unchanged
entity kinds. The binding moves the Rust result into the existing immutable wrapper.
All direct removals mutate the editor without consuming it. They accept current integer IDs in
any order, including repetitions; unlike Edit application, direct calls deduplicate removals.
The existing Python ID bounds check translates Rust's asserted-domain prerequisite into IndexError
before mutation. A consumed editor raises ConsumedError; bounds and reference errors leave an
available editor usable.

Compaction exposes source_count and removed as immutable properties. It omits Rust's into_parts:
that method transfers the owned Rust vector, while a Python tuple still requires constructing
Python values. Consuming the Python compaction would invalidate aliases without preserving that
allocation-reuse benefit.

Edit.RemoveTopology(atoms=..., bonds=...) has the same exact semantics, with batch handles instead
of current editor IDs. Constructing an Edit does not check it against a molecule; application does.
Edit.RemoveTopologyCascading(atoms=..., bonds=...) explicitly includes incident bonds and overlays
with unavailable stored references. Both variants expose immutable atoms and bonds tuples and
support ordinary construction and property access on Python 3.9. Molecule constraint compaction,
completed-prefix behavior, batch transfer, and transaction recovery follow the existing contracts.
Edits.remove_topology and Edits.remove_topology_cascading append the corresponding variant and
return None. Both accept initial integer IDs and New handles; contextual checks occur on application.

A batch stays consumed after transfer, even when application or publication fails. An unavailable
batch leaves the molecule or editor unchanged. A consumed editor rejects application before taking
the batch. Use Edits.copy() to apply a prepared batch more than once:

```python
from umol import AtomForm, Edits, Element, Molecule

source = Molecule.from_entries([AtomForm(Element("C"))])
edits = Edits()
edits.add_atom(AtomForm(Element("N")))
product = source.apply(edits.copy())

editor = source.edit()
editor.apply(edits)
assert editor.finish() == product
assert source == Molecule.from_entries([AtomForm(Element("C"))])
```

Molecule.transact and tracked_transact accept an iterable of Edits and borrow the receiver.
transact returns None; tracked_transact returns the whole transaction's correspondence. On
application or integrity failure they restore the receiver under Rust's Molecule::normalized_eq,
preserving participant order. Python exposes these prepared-batch operations without scoped
Transaction or TrackedTransaction handles, editor probe, or edit_into callbacks, avoiding Rust
borrows retained across Python calls.

Both methods accept keyword-only rollback_method, defaulting to RollbackMethod.Journal;
RollbackMethod.Clone remains selectable. Rust requires the selector explicitly. Journal records
undo for the applied changes, while Clone edits one working copy and installs it only on success.
Journal is the prepared-batch default because its recovery storage follows the changes rather
than a whole working molecule. The choice changes storage and execution cost, not recovery or
correspondence guarantees.

The iterable is fully extracted before any batch is consumed. Batches are then consumed in order,
before molecule execution. If a later batch cannot be transferred, earlier batches remain consumed
and the molecule and its views remain unchanged. Thus `[edits, edits]` consumes the first occurrence
and fails on the second. Once execution starts, rollback restores the molecule, not its old views
or consumed Edits. Edit failures raise TransactionError; publication failures raise
InvalidStructureError.

### SMILES ingestion

Molecule.from_smiles uses resolve=True by default and delegates to Rust's
ingest::ingest_smiles_with. With resolve=False it delegates to
ingest::ingest_smiles_unresolved_with, applying format-owned defaults and converting
to graph IR without resolution. Undetermined fields and raised constraints are preserved.
io_config applies in both modes; chemistry_model and resolve_config are ignored when
resolve=False. Parsing and conversion retain the existing ParseError and
ModelConversionError mappings.

Reaction.from_reaction_smiles likewise defaults to resolve=True and delegates to
ingest::ingest_reaction_smiles_with. With resolve=False it delegates to
ingest::ingest_reaction_smiles_unresolved_with, raising both sides without resolution
and constructing the reaction from their atom-map correspondence. Raised constraints
and overlays are preserved; the product is reindexed into the lhs-anchored reaction
frame. io_config applies in both modes; chemistry_model and resolve_config are ignored
when resolve=False. Both modes reject nonempty agents and atom-map labels repeated
on either side; labels present on only one side remain permitted. Side-conversion
failures retain the existing ModelConversionError mapping and identify the failing side.

The keyword selects between Rust operations while retaining one Python constructor
per input format. It introduces no additional interpretation in the binding.

### Resolution

Molecule.resolve resolves a copy of its receiver and returns Solution. Only Determined contains a
molecule; Underdetermined is payload-free and Contradictory contains its diagnostic. Planning errors
raise RuntimeError. Every outcome preserves the source and its existing views.

Molecule.resolve_into mutates the receiver only on Determined, whose molecule field is None.
Underdetermined, Contradictory, and planning errors restore the receiver under Rust's
Molecule::normalized_eq, preserving participant order. Existing views are invalidated when
execution starts, even if the molecule is restored.
Its keyword-only rollback_method selects RollbackMethod.Clone (the default) or
RollbackMethod.Journal with the same recovery guarantee. The selector covers all phases in one
transaction. Rust requires the selector explicitly; Python defaults to Clone here because one
working copy avoids recording undo across the resolution pipeline. Neither default promises the
cheaper choice for every workload; callers can select either method.

Molecule.resolve_with_report and resolve_into_with_report return `(solution, report)`. The report
is present on Determined and Underdetermined, and None on Contradictory. Solution itself has no
report field. Returning a separate report preserves one Python Solution variant hierarchy for
the Rust methods' different payload types. Ordinary resolution does not construct a report.
resolve_with_report preserves the source and its views; resolve_into_with_report follows
resolve_into's mutation, invalidation, and rollback selection contract.

### Transformations

Rust's Transform extension trait supplies named methods on Molecule, delegating to engines that
implement Transformer. Transformer retains transform, transform_into, and transform_iter for
configured engine use. Python exposes the named Molecule methods and their configuration values:

| Copy-producing method | In-place method | Arguments after self |
| --- | --- | --- |
| aromatize | aromatize_into | model, *, config=None |
| kekulize | kekulize_into | node_order, *, config=None |
| delocalize_charge | delocalize_charge_into | none |
| fold_hydrogens | fold_hydrogens_into | none |
| unfold_hydrogens | unfold_hydrogens_into | none |

model is an explicitly supplied AromaticityModel. node_order is an explicit sequence of atom
indices determining the localization order. The methods do not choose a model or compute an order.
Omitting config, or passing None, selects the corresponding Rust Default: AromaticityConfig for
aromatization and KekulizeConfig for kekulization. Rust requires these configuration arguments
explicitly; Python's keyword-only defaults follow the existing configuration calling convention.

KekulizeConfig accepts the keyword-only algorithm argument and exposes it as a read-only property.
The default is MaximumMatchingAlgorithm.Edmonds(); MaximumMatchingAlgorithm.HopcroftKarp() is also
available. Configuration and algorithm values are immutable and compare by value.

Rust's base methods consume the molecule. Python's base methods clone the Rust receiver once,
sharing storage initially, and return an independently mutable Molecule. Every outcome preserves
the source, its aliases, and its views. The _into methods return None and change the receiver on
success. Each operation plans before mutation, so rejection leaves the molecule exactly unchanged
without a transaction or rollback. Starting in-place execution invalidates existing views and
iterators even on rejection or a no-op; argument extraction or borrowing failure before execution
preserves them.

Compose transformations through successive calls. Each publishes independently; failure of a
later _into call leaves earlier successful changes applied. To retain a source through a sequence,
use a copying method first, then _into methods on its result:

```python
from umol import Molecule, MoleculeDefaults

source = Molecule.parse(
    '{:atoms ["C#h3" "H"] :bonds [[0 1 "1"]]}',
    defaults=MoleculeDefaults.concrete(),
)
product = source.fold_hydrogens()
assert len(source.atoms) == 2
assert len(product.atoms) == 1

product.unfold_hydrogens_into()
assert len(product.atoms) == 5
assert len(source.atoms) == 2
```

Both methods in a pair use the same exceptions and retain the complete Rust diagnostic:

| Operation failure | Python exception |
| --- | --- |
| Aromatization setup failure | RuntimeError |
| Missing data required by aromatization or kekulization; nonliteral count backing a virtual stereo hydrogen during unfolding | UnderdeterminedError |
| Aromatization contradiction; contradicted assertions or valence invariants in aromatization or kekulization | ContradictionError |
| Missing or duplicate atoms in the supplied Kekulé order | ValueError |
| Kekulization electron-count length mismatch; negative affected hydrogen count; unfolding count/stereo-frame mismatch | InvalidStructureError |
| Unsupported kekulization state or unavailable matching; folding incompatible with a stereo frame | ModelConversionError |

Charge delocalization has no operation error. Ordinary Python argument-extraction errors retain
their existing behavior.

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
Entries in an Edits sequence share one handle namespace. Independent batches are supplied
separately to transact or tracked_transact, each with its own namespace; raw list extension
cannot establish that contract.
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

A value-backed entity `*ConstraintsForm` is an owned mutable container; an entity
`*ConstraintsView` writes through to the stored entity constraints. These writes do not require
an editor. A view obtained from an entity Update retains that immutable owner and borrows its
constraints. All writes through that view, including indexed ring counts, raise TypeError.
Molecule-level `ConstraintsView` is a read-only sequence. Its iterator yields one
stored constraint per step, including an entire `And` or `Or` as one entry. Accessing a
composition's children produces another read-only sequence; iteration never descends implicitly.

Molecule-backed `Constraint` entries retain their owner and location. Variant fields read the
stored constraint; leaf payloads are converted to the existing immutable value forms when
requested. `copy()` produces an independent constraint value. Structural mutation invalidates
the collection, entries, child sequences, and iterators through the molecule's accessor counter.
Ordinary entity assignment leaves these constraint accessors valid. Standalone `Constraints`
remains mutable and returns independent entry values.

## Construction and conversion

Construction of an immutable owner must not retain a mutable Python alias that can later change the
owner. Convert mutable inputs to the owner's immutable representation at the construction boundary.
This is an ownership rule, not semantic validation; validation follows the data-type contracts.

Entity forms are ordinarily writable owned objects. Forms retained in Edit or Delta are read-only,
including nested fields and constraints. Such access rejects mutation with `TypeError`.
`copy`, `normalize`, `meet`, and `join` produce ordinary
writable forms. Edits and Deltas provide append-only container mutation; retrieved entries do not
permit mutation of stored entries.

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
