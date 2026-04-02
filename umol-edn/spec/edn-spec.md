# umol-edn Specification

This document clarifies the
[EDN specification](https://github.com/edn-format/edn) for the purposes of the
umol-edn implementation, resolving points where the original specification is
ambiguous.

## 0. Encoding

EDN input is UTF-8. In Rust this is enforced by the `&str` type — all public
parsing entry points accept `&str`, which is guaranteed valid UTF-8 by the
language. There is no encoding detection or BOM handling.

## 1. Whitespace and commas

Whitespace characters: space (0x20), tab (0x09), newline (0x0A), carriage return
(0x0D). Comma (0x2C) is also treated as whitespace. Form feed (0x0C) is **not**
whitespace.

Whitespace separates tokens and is otherwise insignificant. Leading and trailing
whitespace around values is ignored.

## 2. Comments

Line comments begin with `;` and extend to the end of the line (or end of
input). They are equivalent to whitespace.

## 3. Discard

`#_` followed by a value discards that value. The discarded value must be
syntactically valid. `#_ #_ a b` discards both `a` and `b`. `#_` applied to a
tagged literal discards the entire tagged form (e.g., `#_ #inst "2024"` discards
the whole `#inst` form).

## 4. Nil, booleans

- `nil` — the nil value
- `true`, `false` — boolean values

These are case-sensitive. `Nil`, `TRUE`, etc. are symbols, not literals.

## 5. Integers

Grammar: `[+-]? digit+`

Parsed as signed 64-bit integer (`i64`). Overflow is an error (without `bignum`
feature).

The `N` suffix (`42N`) requests arbitrary-precision. Without the `bignum`
feature this is an error. With `bignum`, parsed as `BigInt`.

No integer other than 0 may begin with 0. Leading zeros (`007`) are rejected.

**Ambiguity resolution:** Without `bignum`, integer overflow is always an error.
With `bignum`, overflow promotes to `BigInt`.

## 6. Floating-point numbers

Grammar: `[+-]? digit+ ('.' digit+)? ([eE] [+-]? digit+)?`

A number containing `.` or `e`/`E` is a float. Parsed as `f64`.

The `M` suffix (`3.14M`) requests arbitrary-precision decimal. Without the
`bignum` feature this is an error.

`-0.0` is valid and distinct from `0.0` at the bit level.

## 7. Strings

Delimited by `"`. Supports escape sequences:

| Escape   | Character                         |
| -------- | --------------------------------- |
| `\\`     | backslash                         |
| `\"`     | double quote                      |
| `\n`     | newline (0x0A)                    |
| `\r`     | carriage return (0x0D)            |
| `\t`     | tab (0x09)                        |
| `\uNNNN` | unicode code point (4 hex digits) |

Unterminated strings are an error. Unknown escape sequences are an error.

**Ambiguity resolution:** The EDN spec lists "standard C/Java escape characters
`\t, \r, \n, \\ and \"`" — exactly five. `\uNNNN` is not explicitly listed but
is universally expected. umol-edn accepts `\uNNNN` because rejecting it would
break real-world EDN files.

## 8. Characters

Preceded by `\`. Forms:

| Form              | Example                |
| ----------------- | ---------------------- |
| Single character  | `\a`, `\Z`, `\!`       |
| Named: `\newline` | newline (0x0A)         |
| Named: `\return`  | carriage return (0x0D) |
| Named: `\space`   | space (0x20)           |
| Named: `\tab`     | tab (0x09)             |
| Unicode: `\uNNNN` | 4 hex digits           |

A character literal must be followed by a non-symbol character, whitespace, or
end of input. `\abc` is an error (multi-character sequence that isn't a named
character).

`\ ` (backslash followed by literal space) is invalid.

## 9. Keywords

Grammar: `:` followed by a symbol name.

The name follows symbol rules: must start with a letter or permitted
punctuation. Digit-starting keywords (`:0`, `:1`) are not valid.

Keywords may be namespace-qualified: `:ns/name`. The `/` separates the namespace
from the name. A bare `:` with no following name is an error.

## 10. Symbols

Grammar: starts with a letter or one of `. * ! _ ? $ % & = < > /`, followed by
any of those characters plus digits, `+`, `-`, `#`, `:`, `'`.

If `-`, `+` or `.` are the first character, the second character (if any) must
be non-numeric.

Reserved symbols: `nil`, `true`, `false` are parsed as their respective
literals, not as symbols.

`/` alone is a valid symbol (the division symbol in Clojure).

Symbols may be namespace-qualified: `ns/name`. `/` can be used once only in the
middle of a symbol. Neither the prefix nor the name part can be empty when `/`
is present. The name after `/` must follow first-character restrictions for
symbols.

## 11. Lists

Delimited by `(` `)`. Ordered sequence of zero or more values. `()` is an empty
list.

## 12. Vectors

Delimited by `[` `]`. Ordered sequence of zero or more values. `[]` is an empty
vector.

## 13. Maps

Delimited by `{` `}`. Sequence of key-value pairs. Keys and values are arbitrary
EDN values. `{}` is an empty map.

An odd number of elements (key without value) is an error.

**Duplicate key handling:** Configurable via `DuplicateKeyPolicy`:

- `Error` (default) — duplicate keys are rejected
- `LastWins` — last value for a given key is kept

**Ambiguity resolution:** The EDN spec is silent on duplicate keys. Clojure's
reader allows them (last wins). umol-edn defaults to error because silent
overwrites mask data bugs.

**Ordering:** Maps use `BTreeMap` internally. Iteration order is deterministic
(sorted by key) regardless of insertion order.

## 14. Sets

Introduced by `#{`, closed by `}`. Unordered collection of unique values. `#{}`
is an empty set.

Duplicate elements in a set are not currently detected at parse time. Sets use
`BTreeSet` internally for deterministic iteration order.

## 15. Tagged literals

Grammar: `#` followed by a symbol (the tag), then a value.

The tag symbol must follow symbol rules. Unqualified tags are reserved for
built-in use; user-defined tags must be namespace-qualified
(`#myapp/Person {...}`). Unqualified tags that are not registered in
`TagReaders` are rejected.

Built-in tag support (feature-gated):

| Tag     | Feature  | Rust type               |
| ------- | -------- | ----------------------- |
| `#inst` | `chrono` | `chrono::DateTime<Utc>` |
| `#uuid` | `uuid`   | `uuid::Uuid`            |

Without the corresponding feature, tagged values parse as
`Tagged(tag_string, Box<Edn>)` and are available for application-level dispatch.

## 16. Ambiguity resolutions

Summary of decisions on underspecified areas of the EDN spec:

| Topic                                 | Decision                                      | Section |
| ------------------------------------- | --------------------------------------------- | ------- |
| Integer overflow                      | Error without `bignum`; promote with `bignum` | 5       |
| String escapes beyond the listed five | `\uNNNN` accepted                             | 7       |
| Duplicate map keys                    | Configurable; error by default                | 13      |
| `#_` nesting (`#_ #_ a b`)            | Each `#_` discards next form                  | 3       |
| Whitespace definition                 | Space, tab, newline, CR, comma. Not form feed | 1       |
| `#_` + tagged literal                 | Discards entire tagged form                   | 3       |
| Leading zeros (`007`)                 | Rejected                                      | 5       |
| `-0` / `-0.0`                         | Valid                                         | 5, 6    |
| `/` as symbol                         | Valid                                         | 10      |

## Round-tripping

Semantic round-trip is guaranteed: parse → display → parse produces an equal
`Edn` value. Lossless round-trip (preserving whitespace, comments, formatting)
is not a goal; that requires a CST representation.
