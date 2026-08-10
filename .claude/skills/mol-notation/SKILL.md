---
description: How to denote molecules, atoms, and bonds in the umol DSL — the EDN molecule map (`{:atoms [...] :bonds [...]}`), the atom-string and bond-string `#`-predicate grammars (`#h #c #a #R …`), element/value wildcards and query patterns (`*`, `#R+`, sets, bool-exprs), and the construction macros (`mol!`, `mol_ground!`, `dsl!`, `atom!`, `bond!`, `atom_ground!`, `bond_ground!`). Use whenever writing a molecule/atom/bond from a string — test fixtures, examples, inline data, or substructure query patterns — or when unsure which macro or token to use. Normative grammar: `umol-graph-ir/spec/umol-dsl-spec.md`.
---

# umol molecule notation

The normative grammar is `umol-graph-ir/spec/umol-dsl-spec.md`. Read it for edge cases
(value-expr precedence, stereo strings, constraints, dative/noncovalent/multicenter).
This is the practical reference for the common 90%.

## Construction macros (umol-graph-ir)

| Macro | Input | Produces | Use for |
|---|---|---|---|
| `mol!("{...}")` | molecule EDN | `MoleculeAst`, **un-grounded** (wildcards preserved, metadata dropped) | **query/substructure patterns**; molecules with wildcards |
| `mol_ground!("{...}")` | molecule EDN | `MoleculeAst` with `MoleculeDefaults::ground()` applied | **concrete test molecules** |
| `dsl!("{...}")` | molecule EDN | `MoleculeDsl` (keeps metadata: ids, aliases) | when ids/aliases matter |
| `atom!("C#h3")` | atom-string | `AtomAst`, un-grounded | a single pattern/atom |
| `atom_ground!("C#h3")` | atom-string | grounded `AtomAst` | a single concrete atom |
| `bond!("2")` / `bond_ground!("2")` | bond-string | `BondAst` | a single bond |

`mol!` vs `mol_ground!` is the key choice: **patterns → `mol!`** (Undetermined stays
a match wildcard), **concrete molecules → `mol_ground!`** (defaults fill the slots).

## Molecule EDN map

```clojure
{:atoms ["C#h3" "C#h2" "O#h1"] :bonds [[0 1 "1"] [1 2 "1"]]}   ; ethanol
```

- `:atoms` — vector of atom-strings. An entry may also be `[:id "C#h3"]` (inline id
  keyword) or a bare keyword referencing `:atom-aliases`.
- `:bonds` — vector of `[a b "bondstr"]` (endpoints are 0-based indices or id keywords),
  or `{:a 0 :b 1 :attrs "1"}`. `:attrs` may be a keyword shorthand: `:single` `:double`
  `:triple` `:quadruple`.
- Optional keys: `:dative-bonds :aromatic-systems :multicenter-bonds :noncovalent-bonds
  :stereo-atoms :stereo-bonds :atom-aliases :constraints :guards`. Each structural
  entry needs `:attrs`.
- **Aromaticity is not a bond order** — never write order 1.5/"aromatic". Use an
  `:aromatic-systems` entry plus ordinary (Kekulé) `:bonds`; the atom `#a` / bond `#a`
  flags mark membership.

## Atom-string: `element` then `#tag payload`

Element first, then zero+ predicates in any order (canonical order is `#i #c #h #n #u
#s #v #d #t #a #m #D #X #x #H #R #T`). At most one predicate per tag.

**Element**: `C`, `Cl` (IUPAC casing) · `*` any · `{C,N,O}` set · `!H` not · `?e ::
{C,N}` bind · `?e` ref. (`*`, sets, binds, refs are query/rule only — invalid in ground.)

**Inherent-field tags** (identify the atom; participate field-wise in matching):

| Tag | Meaning | Notes |
|---|---|---|
| `#i` | isotope mass | `#i=` natural · `#i13` mass 13 · `#i*` any · bare `#i`=1 |
| `#c` | formal charge | **sign required**: `#c+` (=+1) `#c-` `#c+2` `#c-2` `#c0`; bare `#c` invalid |
| `#h` | implicit H | `#h3` · bare `#h`=1 · `#h*` wildcard |
| `#n` | lone pairs | bare=1 |
| `#u` | unpaired electrons | bare=1 |
| `#s` | spin multiplicity (2S+1) | bare=1 |
| `#a` | aromatic π contribution | `#a`/`#a1`=1 · `#a0` (in system, 0 π) · `#a2` · `#a+` member, π unspecified · `#a*` no constraint · `#a!` **not** in any aromatic system |
| `#m` | multicenter valence | same shape as `#a` (`#m+ #m* #m!`) |
| `#v` | localized valence (σ; heavy-neighbor bond-order sum, excludes H/dative/aromatic) | |
| `#d` / `#t` | dative donated / accepted pairs | |

**Derived predicates** (topology queries against the target; never affect grounding):

| Tag | SMARTS | Meaning |
|---|---|---|
| `#D` | D | degree (neighbor count) |
| `#X` | X | connectivity (degree + implicit H) |
| `#x` | x | ring connectivity (ring-bond count) |
| `#H` | H | total H (implicit + explicit) |
| `#R` | R | ring membership — `#Rn` total count n · `#R(6)n` rings of size 6 · bare `#R`/`#R(6)`=1 · `#R+` ≥1 ring · `#R!` acyclic · `#R*` no constraint |
| `#T` | — | tetrahedral stereo — `#Tn` coset · `#T+` stereocenter, coset unspecified · `#T!` not a stereocenter · `#T*` no constraint |

## Bond-string: `order` then `#tag payload`

```
"1"  "2"  "3"  "4"      ; discrete orders
"*"                     ; any order (query wildcard)  ← use for ~ patterns
"{1,2}"                 ; order set
```

Bond predicates (own namespace): `#c` charge (sign rules as atoms), `#u`/`#s` spin,
`#a` aromatic membership (no payload), `#R` ring membership (as atoms), `#C` cis/trans
stereo (`#Cn`/`#C+`/`#C!`/`#C*`). Keyword shorthands `:single :double :triple
:quadruple` stand in for `"1".."4"`.

## Value-payload conventions (atom + bond)

| Payload | Means |
|---|---|
| bare (e.g. `#h`) | 1 (decimal-only count tags); `#i`=mass 1; `#c` bare is **invalid** |
| `n` | exact literal |
| `+` / `-` on `#c` | +1 / −1 (charge only) |
| `*` | wildcard / Undetermined (query) |
| `{a,b,c}` | finite set (matches if target is one of) |
| `!x`, `!{…}` | negation / negative set (query) |
| `?id`, `?id :: {…}`, `?h+1`, `?h>=1` | bind / arithmetic / bool-expr (query/rule) |
| special `+` `!` `*` on `#a` `#m` `#R` `#T` | symbolic states (see tables) — not arithmetic |

## Ground vs query (matching)

- **Ground** (`mol_ground!`, `*_ground!`): every inherent slot concrete; no `*`, sets,
  `?id`, bool-exprs, `element-bind/ref`.
- **Query/pattern** (`mol!`, `atom!`, `bond!`): wildcards, sets, `#R+`, bool-exprs
  allowed. A pattern matches a target iff, per slot, the pattern's solution-set ⊇ the
  target's (the pattern admits everything the target does). `*` matches anything; an
  empty/`default()` atom matches any atom (`AtomAst::default()` ≡ `"*"`).

## Common patterns

```rust
mol_ground!(r#"{:atoms ["C#h3" "C#h3"] :bonds [[0 1 "1"]]}"#)        // ethane (concrete)
mol!(r#"{:atoms ["*" "*"] :bonds [[0 1 "*"]]}"#)                      // [*]~[*] (any 2 atoms, any bond)
mol!(r#"{:atoms ["*#R+" "*#R+" "*#R+"] :bonds [[0 1 "*"] [1 2 "*"] [2 0 "*"]]}"#)  // 3-ring of ring atoms
atom_ground!("N#c+#h4")                                               // ammonium N
```

## Gotchas

- `#c` charge **must** carry sign/zero; bare `#c` is a parse error.
- Bare count payload = 1, not 0 (`#h` = 1 H). Use `#h0` / `#R!` for zero.
- Aromatic ≠ a bond order; declare `:aromatic-systems` + Kekulé `:bonds`.
- Wildcard bond order is `"*"` (e.g. `[0 1 "*"]`), wildcard atom is `"*"`.
- `#R`/`#D`/`#X`/`#x`/`#H`/`#T` are **derived** — they filter matches against the
  target graph and never change whether an atom is "ground".
- Use `mol!` for anything with wildcards; `mol_ground!` will try to fill them.
