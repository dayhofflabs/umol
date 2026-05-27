# Fix DSL parsing

## Context

The current molecule DSL implementation in `umol-graph/src/dsl/molecule.rs` has accumulated multiple layers of responsibility in one place:

- EDN map parsing
- fast-path `from_edn_str` parsing
- alias and id resolution
- semantic lowering into `MoleculeAst`
- compact DSL rendering
- metadata management for ids and aliases

This produced a design with:

- `MoleculeAstWrapper`
- `RawMoleculeAst`
- sidecar `Metadata`
- mixed tree-walking and streaming parsing logic
- `edn.to_string()` reparsing in `FromEdn`
- `AtomPattern` / `BondPattern` leaking into the molecule representation layer

The code works in many cases, but the structure is no longer readable enough to evolve safely.

## Goals

The cleanup should preserve the actual requirements:

1. Keep both EDN entry paths:
   - `from_edn`
   - `from_edn_str`
2. Keep the fast path. `from_edn_str` has repeatedly shown real performance wins and stays.
3. Allow `MoleculeAst` to serialize without metadata, even if the output is less compact.
4. Keep compact DSL rendering as an explicit capability rather than a hidden requirement of the semantic AST.
5. Make the architecture legible:
   - EDN shape parsing
   - semantic lowering
   - compact rendering
   - configuration
   should each have one obvious home.

## Core design

The key split is between:

- `*Dsl` types: surface representation
- `*Ast` types: semantic representation

Recommended names:

- `AtomDsl`
- `BondDsl`
- `MoleculeDsl`

These are distinct from:

- `AtomAst`
- `BondAst`
- `MoleculeAst`

### `*Dsl` layer

The `*Dsl` layer owns syntax-level concerns:

- EDN map shape
- ids
- aliases
- symbolic refs
- compact forms
- nested subpatterns in DSL form
- the choice between positional and symbolic references

This is the correct place for:

- `FromEdn`
- `ToEdn`
- `FromEdnMap`
- `ToEdnMap`
- the fused `from_edn_str` parser

### `*Ast` layer

The `*Ast` layer owns semantic concerns:

- resolved indices
- semantic constraints
- no alias table
- no syntax metadata
- no requirement to preserve compactness

This is the correct place for:

- canonical serialization
- matcher / solver / transformation semantics
- semantic invariants

## Trait split

Two trait families should exist because they solve different problems.

### 1. EDN shape traits

`FromEdnMap` / `ToEdnMap` are small convenience traits for map-shaped EDN documents.

They are not lowering traits. They only express:

- "construct this type from an EDN map"
- "render this type as an EDN map"

Suggested shape:

```rust
pub trait FromEdnMap<'de>: Sized {
    fn from_edn_map(map: &EdnMap<'de>) -> Result<Self, DeError>;
}

pub trait ToEdnMap {
    fn to_edn_map(&self) -> EdnMap<'static>;
}
```

These belong primarily on `*Dsl` types.

Benefits:

- reduces repeated `Edn::Map` matching
- centralizes field extraction and unknown-field checks
- gives a principled tree-based path for `from_edn`
- keeps the EDN-facing code small and explicit

### 2. Raising / lowering traits

`FromAst` / `ToAst` are representation-conversion traits.

They should be reintroduced as the main architectural seam and should support configuration.

Suggested shape:

```rust
pub trait ToAst<T> {
    type Error;
    type Config;

    fn to_ast(&self, config: &Self::Config) -> Result<T, Self::Error>;
}

pub trait FromAst<T>: Sized {
    type Error;
    type Config;

    fn from_ast(value: &T, config: &Self::Config) -> Result<Self, Self::Error>;
}
```

These traits belong across representation layers:

- `AtomDsl <-> AtomAst`
- `BondDsl <-> BondAst`
- `MoleculeDsl <-> MoleculeAst`
- `TableIR <-> *Ast`

This is also the right place to retire the current `coerce` / `release` framework. Representation conversion becomes explicit instead of being expressed indirectly through constraint machinery.

## Parsing architecture

Both EDN entry paths should target the same `*Dsl` type.

### Tree path

`from_edn` should decode the EDN tree into `MoleculeDsl`.

Shape:

```rust
impl<'de> FromEdn<'de> for MoleculeDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        // tree-based structural decode
    }
}
```

This path should not stringify and reparse the tree.

### Fast path

`from_edn_str` should remain the fused parser, but it should also target `MoleculeDsl`.

Shape:

```rust
impl<'de> FromEdn<'de> for MoleculeDsl {
    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        // fused fast path
    }
}
```

### Shared semantic lowering

After either parse path succeeds:

```rust
MoleculeDsl -> MoleculeAst
```

via `ToAst<MoleculeAst>`.

That gives one semantic lowering implementation instead of duplicating logic across tree and fused parsing code.

## Serialization architecture

Two serialization targets should exist.

### 1. Canonical semantic serialization

`MoleculeAst` should be serializable directly.

This output is allowed to be less compact:

- no aliases required
- no metadata required
- explicit structure preferred

This path is the semantic one and should not depend on sidecar metadata.

### 2. Compact DSL serialization

`MoleculeDsl` should own compact, human-oriented surface rendering:

- ids
- aliases
- compact endpoint refs
- nested subpatterns in DSL form

Compactness is a rendering policy, not a property of `MoleculeAst`.

## Auto-compaction

Auto-compaction is a legitimate feature, but it belongs in raising, not in the semantic AST.

That means:

- `MoleculeAst` stays clean
- `MoleculeDsl::from_ast(value, config)` may choose:
  - canonical rendering
  - preserve ids if available
  - generate ids automatically
  - introduce aliases
  - prefer positional refs

This is exactly why `FromAst` / `ToAst` should accept config.

## Data ownership

### `MoleculeDsl`

`MoleculeDsl` should contain surface-only state such as:

- atom entries
- bond entries
- dative / aromatic / multicenter / noncovalent entries
- atom aliases
- surface constraint forms
- nested subpattern DSL values

Representative helper types:

```rust
enum AtomRefDsl {
    Index(usize),
    Id(String),
}

enum AtomEntryDsl {
    Atom(AtomDsl),
    Alias(String),
    WithId(String, Box<AtomEntryDsl>),
}
```

The exact names may vary, but the principle is fixed: surface refs stay in the DSL layer.

### `MoleculeAst`

`MoleculeAst` should contain only semantic data:

- atoms
- bonds
- dative bonds
- aromatic systems
- multicenter bonds
- noncovalent bonds
- semantic constraints

No ids, aliases, or compactness metadata.

## Constraints

The same split should apply to constraints.

### DSL constraints

`MoleculeConstraintDsl` stores surface references:

- `int | keyword`
- nested `MoleculeDsl` in subpatterns
- any surface forms that still require resolution

### AST constraints

`MoleculeConstraint` stores resolved references:

- `AtomIdx`
- `BondIdx`
- `AromaticSystemIdx`
- `MulticenterBondIdx`
- nested `MoleculeAst` in subpatterns

This fixes the current problem where nested pattern metadata is needed transiently and then lost in ways that make round-tripping awkward.

## Error split

The cleanup should also simplify the error model.

### EDN parsing / shape errors

These belong to:

- `EdnError`
- `DeError`
- `FromEdn`
- `FromEdnMap`

### DSL semantic lowering errors

These belong to:

- alias resolution
- id resolution
- namespace disjointness
- invalid symbolic references
- structural validation during lowering

This should be one molecule-layer lowering error type rather than a maze of remapping.

### Atom / bond subgrammar errors

These remain local subgrammar errors and are embedded where needed.

## Migration plan

### Phase 1: introduce `*Dsl` types

Create:

- `AtomDsl`
- `BondDsl`
- `MoleculeDsl`

Start by moving surface-only representation there without changing behavior.

### Phase 2: move EDN impls to `*Dsl`

Implement:

- `FromEdn`
- `ToEdn`
- `FromEdnMap`
- `ToEdnMap`

for `MoleculeDsl`.

Move the fused parser into:

- `MoleculeDsl::from_edn_str`

The tree path and fast path must both construct the same DSL object model.

### Phase 3: add raising / lowering traits

Reintroduce:

- `ToAst`
- `FromAst`

with config.

Implement:

- `AtomDsl <-> AtomAst`
- `BondDsl <-> BondAst`
- `MoleculeDsl <-> MoleculeAst`

### Phase 4: move resolution logic out of parsing

Alias resolution, id resolution, endpoint resolution, and surface constraint lowering should move from parse-time ad hoc logic into:

- `MoleculeDsl::to_ast(&config)`

This is the main step that removes the current knot in `molecule.rs`.

### Phase 5: add canonical `MoleculeAst` serialization

Implement direct `ToEdn` / `FromEdn` for `MoleculeAst` in canonical, explicit form.

This path should require no metadata.

### Phase 6: move compact rendering into `FromAst`

Implement:

- `MoleculeDsl::from_ast(&ast, &config)`

This is where:

- auto-generated ids
- alias introduction
- compaction strategy
- stable rendering policy

should live.

### Phase 7: retire old framework

Delete:

- `MoleculeAstWrapper`
- `RawMoleculeAst`
- sidecar `Metadata`
- wrapper-specific APIs
- `edn.to_string()` reparsing in `FromEdn`

Replace the old `coerce` / `release` usage with `FromAst` / `ToAst`.

## Immediate code target

The first concrete target is `umol-graph/src/dsl/molecule.rs`.

The likely end state is:

- `dsl/molecule.rs` defines `MoleculeDsl` and its EDN parsing/rendering
- lowering lives in explicit `ToAst` impls
- semantic `MoleculeAst` serialization no longer depends on wrapper state

It may be worth splitting this file once the new architecture is in place, but file decomposition is secondary. The main problem is responsibility separation, not line count by itself.

## Non-goals

This cleanup does not require simplifying the DSL itself.

Dynamic resolution, symbolic ids, and aliases are still acceptable features. The problem is not that these features exist. The problem is that the current implementation mixes:

- surface representation
- semantic representation
- EDN parsing
- optimization
- rendering policy

into one structure.

## Success criteria

The refactor is successful if:

1. `MoleculeDsl::from_edn` and `MoleculeDsl::from_edn_str` both exist and produce the same DSL object model.
2. `MoleculeAst` can serialize without metadata.
3. Compact rendering is possible, but not required for semantic serialization.
4. Nested subpatterns are represented cleanly in both DSL and AST layers.
5. `coerce` / `release` are no longer needed for representation conversion.
6. The molecule DSL code is organized by responsibility rather than by historical accumulation.
