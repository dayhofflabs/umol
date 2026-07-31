# 175 — Ground AST API

Status: **Proposed**
Date: 2026-07-31
Relates: [172](172-ast-literal-extraction-2026-07-30.md),
[173](173-ground-literals-and-spin-2026-07-30.md)

## Scope

This document defines the concrete useful API built around `Ground<T>` across the AST layer. Doc
173 settles the generic witness shape and corrects the exact literal types on which it depends;
this document owns:

- recursive groundness across leaf, entity, relation, constraint, molecule, reaction, and update
  AST families;
- the supported concrete `Ground<T>` forms;
- ground-preserving navigation through molecule and entity views;
- infallible exact access to present ground values while preserving optional constraints;
- the boundary between stored values, topology-derived values, structural integrity, and chemistry
  validation;
- migration of operations that genuinely require fully ground inputs.

The work includes the DSL specification where groundness is described, unit and Python tests,
property strategies and properties, conformance suites, fuzz targets and named seed corpora,
fixtures, examples, benchmarks, and snapshots affected by the final contract. It does not redefine
physical spin validity or the exact literal carriers settled in doc 173.

No staged implementation plan should be written until the concrete API below has been completed.

## Settled generic wrapper

Doc 173 approves a new top-level `umol-ast/src/ast/ground.rs` module and this generic surface:

```rust
pub struct Ground<T>(T);

Ground::new(value) -> Option<Self>
ground.as_ref() -> &T
ground.into_inner() -> T
```

The field is private. `Ground<T>` owns the supplied value or handle, implements `AsRef<T>`, and
does not implement `Deref`. Construction checks structural groundness. There is no public
`Groundable` trait.

`Ground<T>` is evidence of structural groundness only. It does not silently add entity-integrity,
chemistry-model, or physical-validity guarantees.

## Groundness of constraints

Constraints are part of the groundness guarantee. Optionality and underdetermination remain
distinct:

- an absent optional constraint is compatible with a ground entity;
- every present constraint value must be ground;
- a ground constraint accessor returns `Option<Lit>`, where `None` means absent and `Some` is
  necessarily literal;
- there is no present-but-nonliteral result on the ground surface.

Molecule-level constraint groundness is recursive:

- entity and molecule leaves delegate to their contained AST values;
- `And`, `Or`, and `Not` recurse into their children;
- structural predicates without lattice-valued data are ground by construction;
- `SubPattern` delegates to its nested `MoleculeAst`.

`MoleculeAst::is_ground()` must include its molecule-level `Constraints` in addition to all entity
ASTs. The current implementation checks the entity families but omits the molecule constraint
tree; that omission is part of this work unit.

## Known API distinction

Local entity groundness and molecule-context groundness prove different facts.
`AtomView::is_ground()`, for example, currently checks only its stored `AtomAst`; it says nothing
about neighboring bonds or overlays used by topology-derived methods. The API must not construct
the same ground view type from proofs of different strength and then expose methods justified only
by the stronger proof.

The design must distinguish:

- a ground standalone entity AST such as `Ground<&AtomAst>`;
- an entity view obtained from a ground molecule, such as `Ground<AtomView<'a>>`;
- any additional integrity precondition required by a topology-derived projection.

One candidate is to make molecule entity views constructible only through navigation from
`Ground<&MoleculeAst>`, while locally checked values use `Ground<&AtomAst>` and the corresponding
entity AST forms. This is not yet a settled API.

## Concrete design work

The following must be settled before implementation planning:

1. Enumerate every AST type for which `Ground::new` is public and every ground wrapper produced
   only by navigation from another ground value.
2. Define molecule and collection navigation without introducing a parallel hierarchy of wrapper
   types merely to reproduce `AtomViews`, `BondViews`, and the relation-view namespaces.
3. Define the exact stored-field and optional-constraint accessors for every entity family.
4. Classify each topology-derived atom, bond, and overlay method as:
   - total under structural molecule groundness;
   - total only after entity-integrity validation; or
   - inherently fallible/optional even for a valid ground molecule.
5. Decide how the stronger validated boundary is represented where structural groundness alone is
   insufficient, without folding validators into the meaning of `Ground<T>`.
6. Determine whether reaction, delta, edit/update, and partial DSL ASTs have a useful ground wrapper
   or only participate in recursive groundness checks.
7. Define the first consumer migration. Current fingerprint implementations are the clearest
   ground-only consumers, but their required surface must not dictate misleading general AST
   semantics.
8. Specify the property laws, conformance cases, fuzz inputs, named seeds, and benchmark baselines
   for each approved wrapper and navigation operation.

## Initial consumer evidence

The current ECFP, Morgan, WL, pattern, and structural fingerprint implementations validate
`MoleculeAst::is_ground()` and then repeatedly extract:

- atom element, isotope mass, charge, implicit/total hydrogens, heavy-atom degree, and heavy-atom
  valence;
- localized bond order and aromatic-system membership.

These are useful evidence for the first checked surface, not an automatic specification of it.
Stored fields can follow directly from exact leaf projection. Derived degree, valence, hydrogen,
and aromatic-membership methods must be justified against the precise groundness and integrity
guarantees before becoming infallible.

## Testing requirements

At minimum, the eventual plan must cover these property families:

- checked construction succeeds exactly when the corresponding recursive `is_ground` contract
  holds;
- every present value reachable through a ground accessor has an exact literal projection;
- optional constraints preserve absence as `None` and never hide a present non-ground value;
- ground-preserving navigation cannot produce a non-ground child;
- `as_ref` and `into_inner` preserve identity and do not clone the underlying molecule;
- structural groundness alone never claims entity integrity or chemistry validity;
- migrated consumers retain their fixed reference outputs and existing public error behavior.

The property suite is part of the public communication of these guarantees. Unit cases should
pin individual accessor shapes, while property tests state the cross-family laws. Fuzz targets and
named seeds must cover recursive constraint trees, nested subpatterns, partial entity fields, and
every supported ground-construction boundary.
