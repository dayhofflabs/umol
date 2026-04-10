# EDN Subgrammar Design

## Problem

EDN strings and tagged literals frequently embed domain-specific subgrammars.
Examples in umol-graph:

- **Untagged**: `"C#h4#c+"` inside `:atoms` vector — an atom DSL string parsed by nom
- **Untagged**: `"2#c+#u"` inside `:bonds` — a bond DSL string
- **Tagged**: `#umol/atom "C#h4"` — hypothetical tagged variant (not yet used)
- **EDN-native**: `{:atoms [...] :bonds [...]}` — MoleculeAst, parsed by EDN itself

Each subgrammar has its own error type (`dsl::ParseError` for atom/bond, `EdnError` for
molecule). Currently, all subgrammar errors are flattened to `DeError::Custom(String)`,
losing type information and making programmatic error handling impossible.

## Two Paths

### Tree path (`from_edn`)

1. EDN text → `Edn` tree (reader)
2. `Edn` tree → Rust value (`FromEdn::from_edn`)

On this path, strings are opaque `Edn::Str` nodes. The `from_edn` impl must parse the
string content via the subgrammar parser. The subgrammar error must travel through
`DeError`.

### Fused path (`from_edn_str`)

1. EDN text → Rust value directly (streaming reader + inline subgrammar parsing)

On this path, the deserializer has full positional context. It reads an EDN string token,
then feeds the string content to the subgrammar parser. The subgrammar error can be
richer (byte offset into the outer EDN source).

## Design: `DeError::Subgrammar`

### Variant shape

```rust
#[derive(Clone, Debug, Error)]
pub enum DeError {
    // ... existing variants ...

    /// A domain-specific subgrammar embedded in a string or tagged literal
    /// failed to parse or validate.
    #[error("{grammar}: {message}")]
    Subgrammar {
        /// Identifies the subgrammar (e.g. "atom", "bond", "molecule").
        grammar: &'static str,
        /// Human-readable error message from the subgrammar parser.
        message: String,
        /// Path within the EDN structure where the error occurred.
        path: Vec<String>,
    },
}
```

### Why not `Box<dyn Error>`?

`DeError` currently derives `Clone, PartialEq, Eq`. `Box<dyn Error>` breaks all three.

Options considered:
- `Arc<dyn Error + Send + Sync>` — preserves `Clone` but breaks `PartialEq/Eq`
- Drop `PartialEq/Eq` from `DeError` — only `matches!()` used in tests, no `==` checks
- Type-erased but structured — keep `grammar: &'static str` + `message: String`

The structured approach (`grammar` + `message`) is chosen because:
- Preserves all existing derives
- The `grammar` field enables programmatic dispatch without downcasting
- Downstream consumers needing the original error type can re-parse the message or
  keep the source error alongside the `DeError`
- No trait object machinery needed

### Derive impact

`PartialEq` and `Eq` are preserved. `Clone` is preserved. No changes to existing code
that pattern-matches on `DeError`.

### Migration

All current `DeError::Custom(e.to_string())` calls in subgrammar contexts become:

```rust
DeError::Subgrammar {
    grammar: "atom",
    message: e.to_string(),
    path: Vec::new(),
}
```

A convenience constructor avoids boilerplate:

```rust
impl DeError {
    pub fn subgrammar(grammar: &'static str, err: impl std::fmt::Display) -> Self {
        Self::Subgrammar {
            grammar,
            message: err.to_string(),
            path: Vec::new(),
        }
    }
}
```

## Design: Fallible Tag Readers

### Current state

`TagFn` signature:

```rust
pub type TagFn = Arc<dyn for<'a> Fn(Edn<'a>) -> Result<Edn<'a>, ParseError> + Send + Sync>;
```

Tag readers receive the parsed inner `Edn` value and can transform or validate it.
They return `ParseError` on failure, which is a syntactic error type — semantically
wrong for domain validation.

### Proposed change

Tag readers should be allowed to fail with subgrammar errors. Two options:

**Option A: Tag readers return `EdnError`**

```rust
pub type TagFn = Arc<dyn for<'a> Fn(Edn<'a>) -> Result<Edn<'a>, EdnError> + Send + Sync>;
```

Pro: maximally flexible. Con: forces `EdnError` into the parser, which currently only
deals with `ParseError`.

**Option B: Tag readers return `ParseError` (status quo + subgrammar as ParseError variant)**

Add a `ParseError::TagReaderFailed` variant that wraps a subgrammar error. This keeps
the parser's error type uniform.

```rust
ParseError::TagReaderFailed {
    offset: usize,
    tag: String,
    message: String,
}
```

**Option C: Tag readers return a dedicated `TagError`**

A new error type that can represent either parse or domain failures:

```rust
pub enum TagError {
    Parse(ParseError),
    Domain { grammar: &'static str, message: String },
}
```

### Recommendation

**Option B** is simplest. The parser already handles `ParseError`. A tag reader that
needs to validate domain content (e.g. `#umol/atom "C#h4"`) would parse the string,
run the subgrammar, and wrap failures in `ParseError::TagReaderFailed`. The error
message carries the subgrammar detail. On the `FromEdn` side, the same error becomes
`DeError::Subgrammar`.

This avoids introducing a new error type and keeps the parser monomorphic.

## Design: Untagged Subgrammar Readers (Fused Path)

### Problem

On the fused path, `EdnStreamDeserializer` reads tokens directly. When it encounters a
string that is an untagged subgrammar, the caller currently does:

```rust
let s = de.read_string()?;
s.parse::<AtomAst>().map_err(|e| DeError::Custom(e.to_string()).into())
```

This loses the byte offset of the string within the outer EDN source.

### Proposed: `read_subgrammar`

```rust
impl<'de> EdnStreamDeserializer<'de> {
    /// Read a string token and parse its contents as a subgrammar.
    ///
    /// On success, returns the parsed value. On failure, wraps the
    /// subgrammar error in `DeError::Subgrammar` with the byte offset
    /// of the string token in the outer EDN source.
    pub fn read_subgrammar<T: FromStr>(
        &mut self,
        grammar: &'static str,
    ) -> Result<T, EdnError>
    where
        T::Err: std::fmt::Display,
    {
        let offset = self.position();
        let s = self.read_string()?;
        s.parse::<T>().map_err(|e| {
            DeError::Subgrammar {
                grammar,
                message: e.to_string(),
                path: vec![format!("@{offset}")],
            }
            .into()
        })
    }
}
```

For subgrammars that accept keywords too (bond aliases):

```rust
    pub fn read_subgrammar_or_keyword<T: FromStr>(
        &mut self,
        grammar: &'static str,
        aliases: &bimap::BiMap<String, T>,
    ) -> Result<T, EdnError>
    where
        T::Err: std::fmt::Display,
        T: Clone,
    {
        let offset = self.position();
        let s = self.read_string_or_keyword()?;
        if let Some(v) = aliases.get_by_left(s.as_ref()) {
            return Ok(v.clone());
        }
        s.parse::<T>().map_err(|e| {
            DeError::Subgrammar {
                grammar,
                message: e.to_string(),
                path: vec![format!("@{offset}")],
            }
            .into()
        })
    }
```

These are convenience methods — not traits. They avoid repetitive error mapping at each
call site.

## Design: Tree Postprocessing (`FromEdnPostprocess`)

### Problem

Some types deserialize from EDN in two stages:
1. Parse EDN structure into an intermediate (e.g. `MoleculeInput`)
2. Postprocess: resolve aliases, validate cross-references → final type (e.g. `MoleculeAst`)

The intermediate implements `FromEdn`. The final type does not directly correspond to
any single EDN shape — it requires semantic postprocessing.

### Current approach

`MoleculeAst` has a manual `FromEdn` impl that delegates through `MoleculeInput`:

```rust
impl<'de> FromEdn<'de> for MoleculeAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        MoleculeInput::from_edn(edn)?.into_ast().map_err(...)
    }
}
```

This works. The question is whether to formalize the pattern.

### Proposed: No trait

The `MoleculeInput` → `MoleculeAst` pattern is clean enough as a manual `FromEdn` impl.
A `FromEdnPostprocess` or `FromEdnAdapter` trait would add:
- A trait definition
- A blanket `FromEdn` impl (or a macro)
- Associated type for the intermediate

For one or two types, this is more machinery than value. The manual impl is 5 lines and
fully transparent. If more types need this pattern in the future, revisit.

**Decision: no trait. Manual `FromEdn` impls for postprocessed types.**

## Design: Serialization Side

### Subgrammar → EDN string

`ToEdn` is infallible. Subgrammars that serialize to strings just produce `Edn::Str`:

```rust
impl ToEdn for AtomAst {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}
```

No design changes needed — `Display` → string is straightforward.

### Keyword aliases on serialization

Bond types check aliases and emit keywords when possible:

```rust
impl ToEdn for BondAst {
    fn to_edn(&self) -> Edn<'static> {
        if let Some(name) = aliases.get_by_right(self) {
            Edn::Keyword(EdnKeyword::owned(name.clone()))
        } else {
            Edn::Str(Cow::Owned(self.to_string()))
        }
    }
}
```

This is already implemented and correct. No changes.

### Tagged subgrammar serialization

If tagged subgrammars are introduced (e.g. `#umol/atom "C#h4"`), the `ToEdn` impl
would produce `Edn::Tagged`:

```rust
Edn::Tagged(Cow::Borrowed("umol/atom"), Box::new(Edn::Str(Cow::Owned(self.to_string()))))
```

Not needed yet. Untagged strings work for positional contexts. Tagged forms would be
useful when atom/bond strings appear in non-positional contexts where the reader
cannot infer the subgrammar from position.

### `ToEdnStr` / fused serialization

`ToEdn` produces a tree. Printing requires a second pass. For hot types, a fused
`to_edn_string` could skip the tree:

```rust
pub trait ToEdnStr {
    fn to_edn_string(&self) -> String;
}
```

Not needed yet. `to_edn().to_string()` is sufficient — the tree is tiny for atom/bond
types, and molecule serialization is dominated by the vector/map construction, not the
final printing.

**Decision: no `ToEdnStr` trait. Revisit if profiling shows serialization overhead.**

## Summary of Changes

### Implemented

1. `DeError::Subgrammar { grammar, message, path }` variant — added
2. `DeError::subgrammar()` convenience constructor — added
3. `ParseError::TagReaderFailed { offset, tag, message }` variant — added
4. `EdnStreamDeserializer::read_subgrammar()` convenience method — added
5. `AtomDslError` and `BondDslError` — independent error types for atom/bond subgrammars
6. `AtomAst::from_str` returns `AtomDslError`, `BondAst::from_str` returns `BondDslError`
7. All subgrammar `DeError::Custom(e.to_string())` → `DeError::subgrammar(...)` in umol-graph
8. `read_subgrammar_or_keyword` NOT added to umol-edn (bimap is not a dep; callers handle keyword aliases)

### Deferred

- `FromEdnPostprocess` trait — not needed, manual impls suffice
- `ToEdnStr` trait — not needed, `to_edn().to_string()` is fast enough
- Tagged subgrammar forms (`#umol/atom`) — not needed for positional contexts
