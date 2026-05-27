# clojure-reader 0.5.1: issues for EDN serde integration

Crate: `clojure-reader = { version = "0.5.1", features = ["derive"] }`

We use `clojure-reader` for reading and writing EDN in the molecule DSL.
This documents the issues encountered and the workarounds in `dsl/edn_serde.rs`.

**Resolution:** We dropped EDN tagged literals (`#atom`, `#bond`) entirely. Atom and
bond specs are plain strings — context (`:atoms` map values, `:bond` fields, `:aliases`
values) is sufficient for disambiguation. This sidesteps issues 1, 2, 3, and 4 and
is EDN-spec compliant (bare `#tag` names are reserved by the spec for built-in use).

## 1. Tagged literals crash `deserialize_any`

The built-in `Edn<'de>` deserializer's `deserialize_any` does not handle `Edn::Tagged`.
It falls through to a catch-all that returns an error:

```rust
// de.rs line 78
_ => Err(de::Error::custom(format!("Don't know how to convert {self:?} into any")))
```

This makes `#atom "C"` and `#bond "2"` unusable through `from_str::<T>()`.

**Workaround:** `EdnDeserializer` wrapper that unwraps `Tagged(_, inner)` recursively
before dispatching to the visitor.

**Ideal fix in crate:** `deserialize_any` should either unwrap the tag transparently or
delegate to `deserialize_newtype_struct` with the tag name, letting the visitor decide.

## 2. Serializer ignores the tag name in `serialize_newtype_struct`

The built-in serializer's `serialize_newtype_struct` discards the struct name and just
serializes the inner value:

```rust
// ser.rs line 159
fn serialize_newtype_struct<T>(self, _name: &'static str, value: &T) -> Result<()> {
    value.serialize(self)  // tag name is ignored
}
```

So `serializer.serialize_newtype_struct("atom", &"C")` emits `"C"` instead of `#atom "C"`.

**Workaround:** `EdnSerializer` wrapper that prepends `#name ` before the inner value.

**Ideal fix in crate:** emit `#name <value>` when `serialize_newtype_struct` is called,
since this is exactly what EDN tagged literals are for.

## 3. Serializer uses `name/variant` for all enum/variant forms

The built-in serializer emits `#Name/Variant` for unit variants, newtype variants, tuple
variants, and struct variants. This is the Clojure namespaced-keyword convention but does
not match our use case where keywords are simple (`:single`, `:double`).

Example: `serialize_unit_variant("BondSpec", 0, "single")` emits `#BondSpec/single nil`
instead of `:single`.

**Workaround:** custom `Serialize` impls that emit `:keyword` strings directly via the
`EdnKeyword` helper, bypassing the enum serialization path.

## 4. `deserialize_enum` requires namespaced tags (`#ns/variant`)

The built-in deserializer's `deserialize_enum` splits the tag on `/` and rejects bare
tags:

```rust
// de.rs line 250
let mut split = tag.split('/');
let (Some(tag_first), Some(tag_second)) = (split.next(), split.next()) else {
    return Err(de::Error::custom(format!("Expected namespace in {tag} for Tagged for enum")));
};
```

So `#atom "C"` fails with `"Expected namespace in atom for Tagged for enum"`. Only
`#domain/atom "C"` is accepted, which matches Clojure convention but not our DSL.

Combined with issue 3 (serializer emits `#Name/Variant`), the crate's serde layer is
designed exclusively around namespaced tags. Bare tags — which are valid EDN — are
unsupported on both sides.

**Workaround:** `EdnDeserializer` unwraps tagged values before they reach `deserialize_enum`,
so the tag name is discarded and the inner value is deserialized directly. Tag validation
happens in the domain parsers (`parse_atom_dsl`, `parse_bond_dsl`).

**Ideal fix in crate:** support bare tags in `deserialize_enum`, or at least in
`deserialize_newtype_struct` which is the natural serde mapping for `#tag value`.

## 5. `Error` type lacks `Clone` and `PartialEq`

```rust
// error.rs
pub struct Error {
    pub code: Code,        // derives Debug, Eq, PartialEq
    pub line: Option<usize>,
    pub column: Option<usize>,
    pub ptr: Option<usize>,
}
```

`Error` implements `Debug` and `Display` only. `Code` derives `Debug + Eq + PartialEq`
but not `Clone`. Since our `ParseError` requires `Clone + PartialEq` (for `assert_eq!`
in tests), we cannot wrap `Error` or `Code` directly.

**Workaround:** `recode_edn_error()` maps known `Code` variants to `ParseError` variants
(`UnexpectedEOF` -> `Incomplete`, `Serde(msg)` -> `EdnParse(msg)`) and falls back to
`EdnParse(format!("{code:?}"))` for the rest.

**Ideal fix in crate:** derive `Clone` on both `Code` and `Error`.

## 5. `Display` delegates to `Debug`

```rust
// error.rs line 51
impl Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}
```

So `.to_string()` produces `EdnError { code: Serde("..."), line: None, column: None, ptr: None }`
instead of a human-readable message. This made error messages unreadable before we added
the recoding step.

**Ideal fix in crate:** human-readable `Display` impl, e.g. `"unexpected EOF"` or
`"invalid character at line 3, column 7"`.

## 6. `deserialize_any` treats empty maps as `visit_unit`

```rust
// de.rs line 65
Edn::Map(mut map) => {
    if map == BTreeMap::new() {
        visitor.visit_unit()
    } else {
        visitor.visit_map(MapEdn::new(&mut map))
    }
}
```

An empty EDN map `{}` deserializes as unit (`nil`) rather than an empty map. This is
surprising — `{}` is a valid empty map, not nil.

We replicated this behavior in our `EdnDeserializer` for compatibility but it is arguably
a bug.

## 8. Sequence deserialization reverses element order

```rust
// de.rs line 61
Edn::Vector(mut list) | Edn::List(mut list) => {
    list.reverse();
    Ok(visitor.visit_seq(SeqEdn::new(list))?)
}
```

The built-in deserializer reverses vectors and lists before visiting. The `SeqEdn`
implementation then pops from the end. This works but means the `Edn` representation
stores elements in reverse order internally after the reverse call. Our `EdnDeserializer`
avoids this by iterating forward with `into_iter()`.

## 9. `from_str` vs `read` API gap

- `edn::read(s)` returns `(Edn<'_>, &str)` — the remaining unparsed input.
- `de::from_str::<T>(s)` calls `edn::read_string(s)` which silently discards trailing
  content.

There is no `from_str` variant that uses `read` semantics (error on trailing content).
We call `edn::read` ourselves, check for trailing content, then feed the `Edn` value
to our deserializer.

## 10. No `serialize_str` escaping

```rust
// ser.rs line 106
fn serialize_str(self, v: &str) -> Result<()> {
    self.output += "\"";
    self.output += v;
    self.output += "\"";
    Ok(())
```

String serialization does not escape special characters (`"`, `\`, newlines). A string
containing a double quote will produce malformed EDN. Our serializer has the same
limitation since we copied the pattern; we rely on the fact that atom/bond DSL strings
never contain characters that need escaping.

## Summary of workarounds in `dsl/edn_serde.rs`

| Issue | Workaround |
|---|---|
| Tagged in `deserialize_any` | dropped tags — atom/bond specs are plain strings |
| No tag in `serialize_newtype_struct` | dropped tags — `serialize_str` instead |
| `name/variant` enum format | custom `Serialize` with `EdnKeyword` |
| Bare tags rejected in `deserialize_enum` | dropped tags entirely |
| `Error` lacks `Clone`/`PartialEq` | `recode_edn_error()` maps to `ParseError` variants |
| `Display` = `Debug` | recoding extracts the inner message |
| Empty map = unit | replicated for compatibility |
| Sequence reversal | iterate forward instead |
| No trailing-content check in `from_str` | use `edn::read` + manual check |
| No string escaping | acceptable for our DSL strings |
