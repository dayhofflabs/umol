# umol-edn Specification

This document attempts to clarify the
[EDN specification](https://github.com/edn-format/edn) for the purposes of this
implementation where the original contains some ambiguity. umol-edn supports two
dialects: **Edn** (strict) and **Clojure** (default). The Edn dialect tracks the
spec text with one extension (leading `+` on numbers). The Clojure dialect adds
features present in Clojure's reader.

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

**Dialect:** Clojure only. In Edn mode, `#_` is parsed as a tagged literal with
tag `_`.

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

Leading zeros (`007`) are parsed as decimal, not octal. The parser does not
reject leading zeros.

**Extension:** Leading `+` is accepted in both dialects. The EDN spec only
mentions `-` as a sign prefix. umol-edn accepts `+` for symmetry.

**Ambiguity resolution (D2b):** Without `bignum`, integer overflow is always an
error. With `bignum`, overflow promotes to `BigInt`.

## 6. Floating-point numbers

Grammar: `[+-]? digit+ ('.' digit+)? ([eE] [+-]? digit+)?`

A number containing `.` or `e`/`E` is a float. Parsed as `f64`.

The `M` suffix (`3.14M`) requests arbitrary-precision decimal. Without the
`bignum` feature this is an error.

`-0.0` is valid and distinct from `0.0` at the bit level.

### Special float literals

`##NaN`, `##Inf`, `##-Inf` produce `f64` NaN, positive infinity, and negative
infinity respectively.

**Dialect:** Clojure only. In Edn mode, `##` is an error.

## 7. Strings

Delimited by `"`. Supports escape sequences:

| Escape      | Character                         | Dialect      |
| ----------- | --------------------------------- | ------------ |
| `\\`        | backslash                         | both         |
| `\"`        | double quote                      | both         |
| `\n`        | newline (0x0A)                    | both         |
| `\r`        | carriage return (0x0D)            | both         |
| `\t`        | tab (0x09)                        | both         |
| `\b`        | backspace (0x08)                  | Clojure only |
| `\f`        | form feed (0x0C)                  | Clojure only |
| `\uNNNN`    | unicode code point (4 hex digits) | both         |
| `\0`–`\377` | octal byte value                  | Clojure only |

Unterminated strings are an error. Unknown escape sequences are an error.

**Ambiguity resolution (D11):** The EDN spec lists "standard C/Java escape
characters `\t, \r, \n, \\ and \"`" — exactly five. `\uNNNN` is not explicitly
listed but is universally expected. umol-edn accepts `\uNNNN` in both dialects
because rejecting it would break real-world EDN files.

## 8. Characters

Preceded by `\`. Forms:

| Form                | Example                | Dialect      |
| ------------------- | ---------------------- | ------------ |
| Single character    | `\a`, `\Z`, `\!`       | both         |
| Named: `\newline`   | newline (0x0A)         | both         |
| Named: `\return`    | carriage return (0x0D) | both         |
| Named: `\space`     | space (0x20)           | both         |
| Named: `\tab`       | tab (0x09)             | both         |
| Named: `\formfeed`  | form feed (0x0C)       | Clojure only |
| Named: `\backspace` | backspace (0x08)       | Clojure only |
| Unicode: `\uNNNN`   | 4 hex digits           | both         |

A character literal must be followed by a non-symbol character, whitespace, or
end of input. `\abc` is an error (multi-character sequence that isn't a named
character).

`\ ` (backslash followed by literal space) is invalid.

## 9. Keywords

Grammar: `:` followed by a symbol name.

In Edn mode, the name follows symbol rules (must start with a letter or
permitted punctuation). In Clojure mode, the name may also start with a digit
(`:0`, `:1`, `:123abc`).

Keywords may be namespace-qualified: `:ns/name`. The `/` separates the namespace
from the name. A bare `:` with no following name is an error.

## 10. Symbols

Grammar: starts with a letter or one of `. * ! _ ? $ % & = < > /`, followed by
any of those characters plus digits, `+`, `-`, `#`, `:`, `'`.

Reserved symbols: `nil`, `true`, `false` are parsed as their respective
literals, not as symbols.

`/` alone is a valid symbol (the division symbol in Clojure).

Symbols may be namespace-qualified: `ns/name`.

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

The tag symbol must follow symbol rules. Tags are typically namespace-qualified
for user-defined types (`#myapp/Person {...}`). Unqualified tags are reserved
for built-in use.

Built-in tag support (feature-gated):

| Tag     | Feature  | Rust type               |
| ------- | -------- | ----------------------- |
| `#inst` | `chrono` | `chrono::DateTime<Utc>` |
| `#uuid` | `uuid`   | `uuid::Uuid`            |

Without the corresponding feature, tagged values parse as
`Tagged(tag_string, Box<Edn>)` and are available for application-level dispatch.

## 16. Clojure extensions not supported

These Clojure reader features are not part of EDN and are not supported in
either dialect:

| Feature                                   | Reason                                              |
| ----------------------------------------- | --------------------------------------------------- |
| `#'var`                                   | Clojure-specific                                    |
| `@deref`                                  | Clojure-specific                                    |
| `#()` anonymous functions                 | Clojure-specific                                    |
| `'quote`, `` `syntax-quote ``, `~unquote` | Clojure-specific                                    |
| `#?` reader conditionals                  | Clojure-specific                                    |
| `#=` read-eval                            | Arbitrary expression execution                      |
| `#:ns{...}` namespaced maps               | Clojure extension, not in EDN spec                  |
| Rationals (`3/4`)                         | Not in spec; behind `bignum` feature as future work |

## Dialect feature matrix

| Feature                                 | Edn   | Clojure | Section |
| --------------------------------------- | ----- | ------- | ------- |
| `nil`, `true`, `false`                  | yes   | yes     | 4       |
| Integers with `+`/`-` sign              | yes   | yes     | 5       |
| Floats                                  | yes   | yes     | 6       |
| `##NaN`, `##Inf`, `##-Inf`              | no    | yes     | 6       |
| Strings                                 | yes   | yes     | 7       |
| `\t`, `\r`, `\n`, `\\`, `\"` in strings | yes   | yes     | 7       |
| `\uNNNN` in strings                     | yes   | yes     | 7       |
| `\b`, `\f` in strings                   | no    | yes     | 7       |
| Octal escapes in strings                | no    | yes     | 7       |
| Characters (single, named, `\uNNNN`)    | yes   | yes     | 8       |
| `\formfeed`, `\backspace` characters    | no    | yes     | 8       |
| Keywords (`:foo`, `:ns/name`)           | yes   | yes     | 9       |
| Digit-starting keywords (`:0`, `:1`)    | no    | yes     | 9       |
| Symbols                                 | yes   | yes     | 10      |
| Lists, vectors, maps, sets              | yes   | yes     | 11–14   |
| Tagged literals                         | yes   | yes     | 15      |
| `#_` discard                            | no    | yes     | 3       |
| `;` line comments                       | yes   | yes     | 2       |
| `,` as whitespace                       | yes   | yes     | 1       |
| Bignum suffixes (`N`, `M`)              | error | error   | 5, 6    |

## Ambiguity resolutions

Summary of decisions on underspecified areas of the EDN spec (D15):

| Topic                                   | Decision                                                     | Section |
| --------------------------------------- | ------------------------------------------------------------ | ------- |
| Integer overflow                        | Error without `bignum`; promote with `bignum`                | 5       |
| String escapes beyond the listed five   | `\uNNNN` in both dialects; `\b`, `\f`, octal in Clojure only | 7       |
| `\formfeed`, `\backspace` char literals | Clojure only                                                 | 8       |
| Duplicate map keys                      | Configurable; error by default                               | 13      |
| `#_` nesting (`#_ #_ a b`)              | Each `#_` discards next form                                 | 3       |
| Whitespace definition                   | Space, tab, newline, CR, comma. Not form feed                | 1       |
| `#_` + tagged literal                   | Discards entire tagged form                                  | 3       |
| Leading `+` on numbers                  | Accepted (umol-edn extension)                                | 5       |
| Leading zeros (`007`)                   | Parsed as decimal                                            | 5       |
| `-0` / `-0.0`                           | Valid                                                        | 5, 6    |
| `/` as symbol                           | Valid                                                        | 10      |
| `#_` in Edn dialect                     | Parsed as tagged literal with tag `_`                        | 3       |

## Round-tripping

Semantic round-trip is guaranteed: parse → display → parse produces an equal
`Edn` value. Lossless round-trip (preserving whitespace, comments, formatting)
is not a goal; that requires a CST representation.
