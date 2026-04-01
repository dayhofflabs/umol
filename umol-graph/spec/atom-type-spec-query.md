# Atom Type Spec and Query Notation

This file defines the compact notation used by `graph_ir::AtomTypeSpec` and
`graph_ir::AtomTypeQuery`, including TOML serialization and macro usage.

## Common field syntax

Both specs and queries use the same token syntax after the element symbol:

- `+n` / `-n`: formal charge (`+`/`-` without number means `1`)
- `/n`: lone pairs
- `>n`: donated pairs
- `<n`: accepted pairs
- `^n`: unpaired electrons
- `xn`: multiplicity (defaults to `unpaired_electrons + 1`)
- `Hn`: attached/implicit hydrogens
- `vn`: sigma valence (bond order sum)
- `an`: aromatic valence contribution
- `mn`: multicenter valence contribution

## AtomTypeSpec — `{El...}`

A fully determined atom type. All fields have concrete values; omitted tokens default to `0`
(except multiplicity, which defaults to `unpaired_electrons + 1`).

Examples:

- `{C+0v4}` neutral tetra-valent carbon
- `{N+0/1v3a1}` aromatic N with one lone pair
- `{O-1/3v1}` anionic O with three lone pairs
- `{C-1/1^2x1H1v2a1m2}` full form with explicit spin/aromatic/multicenter fields

## AtomTypeQuery — `?{El...}`

A partial query with optional constraints. Omitted fields are unconstrained (`None`);
explicitly specified fields (including zero values like `+0`) are match constraints.

Examples:

- `?{H}` hydrogen, all fields unconstrained
- `?{H+0}` hydrogen with charge explicitly 0
- `?{C+0v4}` carbon with charge 0 and valence 4, other fields unconstrained
- `?{N+1}` positively charged nitrogen, valence unconstrained

The distinction between `?{H}` (charge undetermined) and `?{H+0}` (charge 0) is the
primary motivation for using queries in test inputs.

## TOML format

`AtomTypeSpec` is serialized as a string:

```toml
atom_types = [
  "{C+0v4a0m0}",
  "{N+0/1v3a1m0}",
  "{O-1/3v1a0m0}"
]
```

`AtomTypeQuery` is serialized as a string:

```toml
atoms = "?{H+0} ?{C+0}"
```

Loading API:

- `AtomTypeRegistry::from_toml_str(...)`
- `AtomTypeRegistry::from_toml_file(...)`
- `AtomTypeOverrides::from_toml_str(...)`
- `AtomTypeOverrides::from_toml_file(...)`

## Macros

Public macros:

- `spec!("{...}")` -> `AtomTypeSpec`
- `query!("?{...}")` -> `AtomTypeQuery`
- `registry!["{...}", "{...}", ...]` -> `AtomTypeRegistry`

The default atom typing registry is also defined from this string notation
using a lazy-initialized macro-expanded registry.
