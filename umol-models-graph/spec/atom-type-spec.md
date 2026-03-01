# Atom Type Spec Notation

This file defines the compact notation used by `graph_ir::valence::AtomTypeSpec`,
including TOML serialization format and macro usage.

## String notation

An atom type spec is a bracketed token list:

`[El ...tokens...]`

- `El`: element symbol (`H`, `C`, `N`, `Fe`, ...)
- Tokens are optional and may appear in any order:
  - `+n` / `-n`: formal charge (`+`/`-` without number means `1`)
  - `/n`: lone pairs
  - `>n`: donated pairs
  - `<n`: accepted pairs
  - `^n`: unpaired electrons
  - `*n`: multiplicity (defaults to `unpaired_electrons + 1`)
  - `Hn`: attached/implicit hydrogens
  - `vn`: sigma valence (bond order sum)
  - `an`: aromatic valence contribution
  - `mn`: multicenter valence contribution

If a token is omitted, the default is `0` (except multiplicity default above).

## Examples

- `[C+0v4]` neutral tetra-valent carbon
- `[N+0/1v3a1]` aromatic N with one lone pair
- `[O-1/3v1]` anionic O with three lone pairs
- `[C-1/1^2*1H1v2a1m2]` full form with explicit spin/aromatic/multicenter fields

## TOML format

`AtomTypeSpec` is serialized as a string, so the registry TOML format is compact:

```toml
atom_types = [
  "[C+0v4a0m0]",
  "[N+0/1v3a1m0]",
  "[O-1/3v1a0m0]"
]
```

Loading API:

- `AtomTypeRegistry::from_toml_str(...)`
- `AtomTypeRegistry::from_toml_file(...)`
- `AtomTypeOverrides::from_toml_str(...)`
- `AtomTypeOverrides::from_toml_file(...)`

## Macros

Public macros:

- `spec!("...")` -> `AtomTypeSpec`
- `registry!("[...]", "[...]", ...)` -> `AtomTypeRegistry`

The default atom typing registry is also defined from this string notation
using a lazy-initialized macro-expanded registry.
