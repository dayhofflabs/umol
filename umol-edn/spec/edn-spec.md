# umol-edn Dialect Specification

umol-edn supports two dialects: **Edn** (strict) and **Clojure** (default).
The Edn dialect follows the [EDN specification](https://github.com/edn-format/edn)
with the single extension of allowing `+` as a numeric sign prefix (see below).
The Clojure dialect adds several features present in the Clojure reader but
absent from the EDN spec.

## Dialect feature matrix

| Feature | Edn | Clojure | Notes |
|---|---|---|---|
| `nil`, `true`, `false` | yes | yes | |
| Integers (`-1`, `0`, `+5`) | yes | yes | Leading `+` is a umol-edn extension to EDN |
| Floats (`1.5`, `1e10`) | yes | yes | |
| Strings (`"hello"`) | yes | yes | |
| `\t`, `\r`, `\n`, `\\`, `\"` in strings | yes | yes | |
| `\uNNNN` unicode escapes in strings | yes | yes | |
| `\b` (backspace) in strings | no | yes | |
| `\f` (formfeed) in strings | no | yes | |
| Octal escapes (`\0`–`\377`) in strings | no | yes | |
| Characters (`\a`, `\newline`, `\space`, `\tab`, `\return`) | yes | yes | |
| `\uNNNN` unicode character literals | yes | yes | |
| `\formfeed` character literal | no | yes | |
| `\backspace` character literal | no | yes | |
| Keywords (`:foo`, `:ns/name`) | yes | yes | |
| Digit-starting keywords (`:0`, `:1`) | no | yes | |
| Symbols (`foo`, `ns/name`) | yes | yes | |
| Lists (`(1 2 3)`) | yes | yes | |
| Vectors (`[1 2 3]`) | yes | yes | |
| Maps (`{:a 1}`) | yes | yes | |
| Sets (`#{1 2 3}`) | yes | yes | |
| Tagged literals (`#tag value`) | yes | yes | |
| `##NaN`, `##Inf`, `##-Inf` | no | yes | |
| `#_` discard | no | yes | In Edn mode, `#_` parses as tagged literal with tag `_` |
| `;` line comments | yes | yes | |
| `,` as whitespace | yes | yes | |
| Bignum suffixes (`N`, `M`) | error | error | Requires `bignum` feature (not yet implemented) |

## Extension: leading `+` on numbers

The EDN spec defines integer as `/ -? digit+ /` and float similarly. umol-edn
extends both dialects to accept an optional leading `+` sign (`+5`, `+3.14`).
This is a deliberate deviation: forbidding `+` while allowing `-` would be an
arbitrary asymmetry with no practical benefit.

## Duplicate key handling

Both dialects support configurable duplicate key policy via `ParseConfig`:

- `DuplicateKeyPolicy::Error` (default) — reject maps with duplicate keys
- `DuplicateKeyPolicy::LastWins` — last value for a given key wins
