# DSL ↔ AST I/O ergonomics

## Motivation

Adding test cases for the aromatic-system equalization rule (doc 93 §3) surfaced a broader friction: writing molecule inputs against `MoleculeAst::new(vec![…], vec![…], vec![], vec![], vec![], vec![], Constraints::default())` is the same friction the EDN DSL was supposed to remove. The DSL is in place but the I/O surface has gaps and asymmetries — the AST types lack the `FromStr`/`Display`/`FromEdn`/`ToEdn` that the DSL types have, there are no construction macros, and the `FromAst`/`IntoAst` traits declare a fallibility that the implementation does not produce.

The aim of this doc is to capture the redesign before any code lands.

## Current state

| Type | `FromEdn` | `ToEdn` | `FromStr` | `Display` | `FromAst` | `IntoAst` |
|---|---|---|---|---|---|---|
| `MoleculeDsl` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `MoleculeAst` | ✓ | ✓ | — | — | n/a | n/a |
| `AtomDsl` | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| `AtomAst` | — | — | — | — | n/a | n/a |
| `BondDsl` / `BondAst` | DSL only | DSL only | DSL only | DSL only | ✓ | ✓ |
| `AromaticSystemDsl` / `AromaticSystemAst` | DSL only | DSL only | DSL only | DSL only | ✓ | ✓ |
| `DativeBondDsl` / `MulticenterBondDsl` / `NoncovalentBondDsl` and their `*Ast` partners | DSL only | DSL only | DSL only | DSL only | ✓ | ✓ |

DSL-only state on the DSL side: `Metadata` carries `atom_ids`, `atom_aliases`, `bond_ids`, `dative_bond_ids`, `aromatic_system_ids`, `multicenter_bond_ids`, `noncovalent_bond_ids`. The AST is positional only.

The lower/raise direction is currently expressed by two traits in `umol_ast::ast::traits`:

```rust
trait FromAst<A>      { type Ctx; type Error; fn from_ast(ast: &A, cfg: &Self::Ctx) -> Result<Self, Self::Error>; }
trait IntoAst<A>      { type Ctx; type Error; fn into_ast(self, cfg: &Self::Ctx) -> Result<A,    Self::Error>; }
```

Every entity-level impl declares `Error = ParseError` but the underlying `lower_*` / `raise_*` functions return `()`; the impls wrap unconditionally with `Ok(...)`. `ConstraintsDsl::from_ast` is already invoked with `.expect("ConstraintsDsl::from_ast is infallible for a well-formed AST")` at `umol-ast/src/dsl/molecule.rs:673`. The fallibility is fictional.

`MoleculeAst::new` takes seven positional `vec!`-shaped arguments; the only common-case shortcut is going through the DSL.

No construction macros exist (`mol!`, `dsl!`, `atom!`, `bond!` are absent).

## Decisions

### 1. Symmetric I/O surface on the AST

Every AST entity type implements `FromStr`, `Display`, `FromEdn`, `ToEdn` — the same surface the DSL types have. Default context is the corresponding `*Defaults::default()`. Concretely:

```rust
impl FromStr for AtomAst { … parse via AtomDsl, lower with AtomDefaults::default() … }
impl Display  for AtomAst { … }
impl<'de> FromEdn<'de> for AtomAst { … }
impl       ToEdn        for AtomAst { … }
```

Same for `BondAst`, `AromaticSystemAst`, `DativeBondAst`, `MulticenterBondAst`, `NoncovalentBondAst`. `MoleculeAst` already has `FromEdn`/`ToEdn`; gains `FromStr` (delegating to `from_edn_str`) and `Display` (delegating to `to_edn().to_string()`).

The default context is the empty / no-substitution form: any value not literally present in the input stays `Undetermined` in the AST, no values filled in. Round-tripping `Display ∘ FromStr` is the identity. The DSL path (with explicit opinionated context, with metadata) is unchanged.

### 1a. `*Defaults::new()` and `Default` impl

The current `*Defaults::verbatim()` constructor (`Required` for every field, no substitution) is renamed to `*Defaults::new()` and a `Default` impl is added pointing to `new()`. Convention pull: `T::new() == T::default() == empty` is what `Vec`, `HashMap`, `String` etc. do; for a config-of-defaults type, "empty" is the no-substitution case. `*Defaults::zeroed()` (the opinionated charge=0/lone_pairs=0/etc. preset) keeps its name.

Applies to: `AtomDefaults`, `BondDefaults`, `AromaticSystemDefaults`, `DativeBondDefaults`, `MulticenterBondDefaults`, `NoncovalentBondDefaults`, `MoleculeDefaults`.

No deprecation path. `verbatim()` goes away in the same change.

### 2. Construction macros — `macro_rules!`

Four runtime-parsing macros, panic on bad input. Definitions live alongside the types they construct (`umol_ast::macros` re-exported at the crate root).

```rust
mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#)   // -> MoleculeAst (metadata dropped)
dsl!(r#"{:atom-aliases {…} :atoms […] :bonds […]}"#)  // -> MoleculeDsl (metadata preserved)
atom!("C#h=#a+")                                    // -> AtomAst
bond!("1#a")                                        // -> BondAst
```

Implementation skeleton (representative):

```rust
#[macro_export]
macro_rules! mol {
    ($s:expr) => {{ <$crate::ast::MoleculeAst as ::core::str::FromStr>::from_str($s).unwrap() }};
}
```

No proc-macro crate. A future `mol_const!` proc-macro for compile-time validation can be added later without breaking call sites — the input shape is the same.

### 3. `MoleculeAst` constructor surface

| Method | Signature | Purpose |
|---|---|---|
| `new()` | `() -> Self` | empty molecule (zero atoms, zero bonds) |
| `from_atoms_and_bonds(atoms, bonds)` | `(Vec<AtomAst>, Vec<(AtomIdx, AtomIdx, BondAst)>) -> Self` | the common case |
| `from_parts(atoms, bonds, dative, aromatic, multicenter, noncovalent, constraints)` | full positional | escape hatch for round-trips and tests covering every entity type |
| `builder()` | `() -> MoleculeBuilder` | fluent / programmatic construction |

Breaking change: the current `MoleculeAst::new(7-args)` is renamed to `from_parts`. Every existing call site migrates by name change only.

### 4. Macro return shape

`mol!` returns `MoleculeAst`. `dsl!` returns `MoleculeDsl`. They are not interchangeable: code that needs metadata uses `dsl!`, code that doesn't uses `mol!`. There is no auto-coercing macro that returns one or the other based on inference.

### 5. Trait split — `FromAst` / `TryFromAst`

```rust
trait FromAst<A>: Sized {
    type Ctx;
    fn from_ast(ast: &A, cfg: &Self::Ctx) -> Self;
}

trait IntoAst<A>: Sized {
    type Ctx;
    fn into_ast(self, cfg: &Self::Ctx) -> A;
}

trait TryFromAst<A>: Sized {
    type Ctx;
    type Error;
    fn try_from_ast(ast: &A, cfg: &Self::Ctx) -> Result<Self, Self::Error>;
}

trait TryIntoAst<A>: Sized {
    type Ctx;
    type Error;
    fn try_into_ast(self, cfg: &Self::Ctx) -> Result<A, Self::Error>;
}
```

Rationale:

- The current Dsl ↔ Ast pair is lossless and infallible. It uses `FromAst` / `IntoAst`.
- `TryFromAst` / `TryIntoAst` is reserved for `TableIR → MoleculeAst`, where raising can fail (the AST is the tree; raising goes from a concrete table-shaped representation up to it). `ExtendedMolecule` Sgroups, for instance, have no faithful AST representation.
- Mirrors `From` / `TryFrom` in `std`. No blanket impl needed.

The `.expect("infallible")` in `MoleculeDsl::to_edn` collapses to direct destructuring.

### 6. Metadata access from `mol!`

`mol!` drops `Metadata`. Tests that need IDs or aliases use `dsl!` and call `.into_parts()` if they also want the AST. The two-macro design is intentional: keeping the most common case (no metadata) as the simplest call.

### 7. Out of scope here

- Schema / typestate replacement of `GroundMolecule` wrappers — interesting direction, distinct work, separate doc when concrete need arises.
- Free functions at the `umol_ast` crate root (`umol_ast::mol(s)`, etc.) — deferred. To be revisited together with SMILES/SMARTS/MOL/SDF/TableIR conversions.
- Compile-time validation of `mol!` / `atom!` / `bond!` literals — deferred until there's evidence the runtime panic is insufficient.
- Source-position tracking in the DSL parser — separate concern.

## Migration plan

The work is strictly additive plus two renames. Nothing in the resolver / validator / transformer surfaces changes.

### Phase 1 — Trait split

1. Define `FromAst` / `IntoAst` (infallible) and `TryFromAst` / `TryIntoAst` (fallible) in `umol_ast::ast::traits`.
2. Rewrite all current Dsl ↔ Ast impls (`AtomDsl ↔ AtomAst`, `BondDsl ↔ BondAst`, `AromaticSystemDsl ↔ AromaticSystemAst`, `DativeBondDsl ↔ DativeBondAst`, `MulticenterBondDsl ↔ MulticenterBondAst`, `NoncovalentBondDsl ↔ NoncovalentBondAst`, `MoleculeDsl ↔ MoleculeAst`) under the infallible traits.
3. Remove `Result` / `?` from composite calls; remove the `.expect("…infallible…")` line in `dsl/molecule.rs`.

### Phase 2 — Constructor surface

1. Rename `MoleculeAst::new(7-args)` → `from_parts(7-args)`. Migrate call sites.
2. Add empty `MoleculeAst::new() -> Self`.
3. Add `MoleculeAst::from_atoms_and_bonds(atoms, bonds) -> Self`.
4. Expose `MoleculeAst::builder() -> MoleculeBuilder` directly (it already exists via `.edit()` for editing — make a clean empty-start variant).

### Phase 3 — `*Defaults` rename + symmetric AST I/O surface

1. Rename `*Defaults::verbatim()` → `*Defaults::new()` across all six `*Defaults` types and the composite `MoleculeDefaults`. Add `impl Default` pointing to `new()`. Migrate every call site of `verbatim()`. Drop the `verbatim` name entirely.
2. For each entity AST type (`AtomAst`, `BondAst`, `AromaticSystemAst`, `DativeBondAst`, `MulticenterBondAst`, `NoncovalentBondAst`):
   - `impl FromStr` — parse via the existing DSL parser, raise with `*Defaults::default()`.
   - `impl Display` — lower with `*Defaults::default()`, format via existing DSL renderer.
   - `impl FromEdn<'de>` — same path.
   - `impl ToEdn` — same path.
3. `MoleculeAst` adds `FromStr` and `Display` (it already has `FromEdn` / `ToEdn`).

### Phase 4 — Macros

Define `mol!`, `dsl!`, `atom!`, `bond!` as `macro_rules!` in `umol_ast::macros`. Re-export at crate root. Each calls the corresponding `FromStr::from_str(s).unwrap()` (or `MoleculeDsl::from_edn_str(s).unwrap()` for `dsl!`).

### Phase 5 — Equalization tests

Now write the equalization tests planned in doc 93 §3, using `mol!` for inputs and direct AST accessors for assertions. This is the original motivation; it lands last because the prior phases unblock it.

## Open

- Connectivity-shortcut name is `from_atoms_and_bonds` for now. If the call site count grows past a small number and the verbosity becomes a real friction, revisit.
- Whether `FromAst` and `TryFromAst` should share a marker / blanket-impl story (e.g., `T: FromAst<A>` automatically gives `T: TryFromAst<A, Error = Infallible>`). Not required for the immediate work.
