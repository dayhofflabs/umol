# EDN spec conformance audit

2026-04-01

Systematic walk-through of the [EDN spec](https://github.com/edn-format/edn).
For each statement we check:

1. What the spec says (verbatim)
2. What Clojure actually does (`clj -e '...'`)
3. What umol-edn does in Edn dialect (strict)
4. What umol-edn does in Clojure dialect (default)
5. What `umol-edn/spec/edn-spec.md` says
6. Whether conformance tests exist

Key: E = Edn dialect, C = Clojure dialect.

---

## S1. Encoding

> edn elements, streams and files should be encoded using UTF-8.

- **Clojure:** UTF-8 (JVM default for tools.reader).
- **umol-edn E/C:** Enforced by `&str` (always UTF-8).
- **Spec doc:** Section 0. Covered.
- **Tests:** Implicit (all inputs are `&str`).

**Status:** Done.

---

## S2. Whitespace

> Whitespace, other than within strings, is not otherwise significant, nor need
> redundant whitespace be preserved during transmissions. Commas `,` are also
> considered whitespace, other than within strings.

- **Clojure:** space, tab, newline, CR, comma. Form feed is whitespace in Clojure.
- **umol-edn E/C:** space, tab, newline, CR, comma. Form feed is NOT whitespace.
- **Spec doc:** Section 1. Covered.
- **Tests:** `test_conformance_whitespace`, `test_conformance_comma_in_collections`.

**Status:** Done. Dialect difference: Clojure treats form feed as whitespace; we do not. Spec is silent on form feed.

---

## S3. Delimiters

> The delimiters `{ } ( ) [ ]` need not be separated from adjacent elements by
> whitespace.

- **Clojure:** `[1[2]]` => `[1 [2]]`, `{:a{:b 1}}` => `{:a {:b 1}}`, `(1(2))` => `(1 (2))`.
- **umol-edn E/C:** Works. `[1]2` parses as two values via `read_all`.
- **Tests:** `test_s3_delimiters_no_whitespace`.

**Status:** Done.

---

## S4. Dispatch character

> Tokens beginning with `#` are reserved. The character following `#` determines
> the behavior. The dispatches `#{` (sets), `#_` (discard), #alphabetic-char
> (tag) are defined below. `#` is not a delimiter.

Two sub-statements:

### S4a. Unknown `#` dispatch is reserved (should error)

- **Clojure:** `#!foo` => "EOF while reading" (error).
- **umol-edn E/C:** Error. Covered.
- **Tests:** `test_read_string_error_unexpected_token` covers `##xyz`.

**Status:** Done.

### S4b. `#` is not a delimiter

> `#` does not terminate tokens. `foo#bar` is a single symbol.

- **Clojure:** `(read-string "foo#bar")` => symbol `foo#bar`.
- **umol-edn E/C:** TODO verify
- **Tests:** `test_conformance_hash_is_not_delimiter`.

---

## S5. Nil

> `nil` represents nil, null or nothing.

- **Clojure:** `(read-string "nil")` => `nil`.
- **umol-edn E/C:** `Edn::Nil`.
- **Tests:** Covered.

**Status:** Done.

---

## S6. Booleans

> `true` and `false` should be mapped to booleans.

- **Clojure:** `(read-string "true")` => `true`.
- **umol-edn E/C:** `Edn::Bool(true)` / `Edn::Bool(false)`.
- **Tests:** Covered.

> If a platform has canonic values for true and false, it is a further semantic
> of booleans that all instances of `true` yield that (identical) value.

- **umol-edn:** `Edn::Bool(true) == Edn::Bool(true)`. Satisfied.

**Status:** Done.

---

## S7. Strings

> Strings are enclosed in `"double quotes"`. May span multiple lines. Standard
> C/Java escape characters `\t, \r, \n, \\ and \"` are supported.

- **Clojure:** Also supports `\b`, `\f`, `\uNNNN`, octal escapes.
- **umol-edn E:** Only `\t \r \n \\ \"` per spec. Rejects `\b`, `\f`, `\uNNNN`, octal.
- **umol-edn C:** Also `\b`, `\f`, `\uNNNN`, octal.
- **Spec doc:** Section 7. Covered.
- **Tests:** `test_conformance_string_escapes`, `test_s7_string_escapes_both_dialects`, `test_s7_string_clojure_escapes_accepted`, `test_s7_string_clojure_escapes_rejected_edn`.

**Status:** Done.

---

## S8. Characters

> Characters are preceded by a backslash: `\c`, `\newline`, `\return`, `\space`
> and `\tab` yield the corresponding characters. Unicode characters are
> represented with `\uNNNN` as in Java. Backslash cannot be followed by
> whitespace.

- **Clojure:** Also `\formfeed`, `\backspace`.
- **umol-edn E:** `\c`, `\newline`, `\return`, `\space`, `\tab`, `\uNNNN`. Rejects `\formfeed`, `\backspace`.
- **umol-edn C:** Also `\formfeed`, `\backspace`.
- **Tests:** `test_conformance_characters`, `test_s8_char_both_dialects`, `test_s8_char_clojure_named_accepted`, `test_s8_char_clojure_named_rejected_edn`.

**Status:** Done.

---

## S9. Symbols

### S9a. Start characters

> Symbols begin with a non-numeric character and can contain alphanumeric
> characters and `. * + ! - _ ? $ % & = < >`.

- **Clojure:** TODO verify full set
- **umol-edn E/C:** `is_symbol_start` covers these.
- **Tests:** `test_conformance_symbol_start_chars`.

### S9b. Sign/dot restriction

> If `-`, `+` or `.` are the first character, the second character (if any)
> must be non-numeric.

- **Clojure:** `+1` => Long 1. `+a` => symbol. `.a` => symbol. `.1` => symbol (not a number!).
- **umol-edn E/C:** `+a`, `-a`, `.a` are symbols. `+1`, `-1` are integers. `.1` is a symbol.
- **Tests:** `test_s9b_sign_dot_first_char`.

**Status:** Done.

### S9c. Interior characters

> Additionally, `: #` are allowed as constituent characters in symbols other
> than as the first character.

- **Clojure:** `(read-string "foo#bar")` => `foo#bar`. `(read-string "foo:bar")` => `foo:bar`.
- **umol-edn E/C:** Covered by `is_symbol_char`.
- **Tests:** `test_conformance_symbol_with_all_char_types`.

### S9d. Slash — single use

> `/` has special meaning in symbols. It can be used once only in the middle of
> a symbol to separate the prefix (often a namespace) from the name, e.g.
> `my-namespace/foo`. `/` by itself is a legal symbol, but otherwise neither
> the prefix nor the name part can be empty when the symbol contains `/`.

Sub-checks:
- `/` alone is valid: covered (`test_conformance_slash_symbol`).
- `ns/name` valid: covered.
- `ns/` invalid (empty name): **Clojure rejects.** umol-edn rejects. Test: `test_s9d_empty_name_after_slash`.
- `/name` invalid (empty prefix): **Clojure rejects.** umol-edn rejects. Test: `test_s9d_empty_prefix_before_slash`.
- `a/b/c` (multiple slashes): **Clojure accepts** (diverges from spec). E: rejects, C: accepts. Test: `test_s9d_multiple_slashes`.

**Status:** Done.

### S9e. Post-slash first-character restriction

> If a symbol has a prefix and `/`, the following name component should follow
> the first-character restrictions for symbols as a whole.

- `foo/bar` valid, `foo/1bar` invalid, `foo/_bar` valid.
- **Clojure:** `foo/1bar` => error. `foo/#bar`, `foo/:bar` => ok (lenient).
- **umol-edn C:** Digits rejected. `#`/`:` after slash accepted (matches Clojure).
- **umol-edn E:** Digits, `#`, `:` after slash all rejected (spec: must be symbol-start).
- **Tests:** `test_s9e_post_slash_digit_rejected`, `test_s9e_keyword_post_slash_digit_rejected`, `test_s9e_post_slash_valid_start`.

**Status:** Done.

---

## S10. Keywords

> Keywords follow the rules of symbols, except they can (and must) begin with
> `:`, e.g. `:fred` or `:my/fred`.

"Follow the rules of symbols" is load-bearing. It means:

### S10a. Bare keywords are valid

`:foo` is valid — no namespace required (unlike `/` which cannot stand alone
for keywords since `:` consumes the prefix slot).

- **Tests:** `test_s10_bare_keyword_valid`.

### S10b. :/ and :/anything are not legal keywords

- **Clojure:** `:/` is accepted (diverges from spec). `:/foo` => error.
- **umol-edn E/C:** Rejects both. Follows spec.
- **Tests:** `test_conformance_keyword_invalid_slash`.

### S10c. The namespace part must be a valid symbol

In `:ns/name`, `ns` must follow symbol first-character rules.

- **Clojure:** `:0/foo` accepted (lenient).
- **umol-edn E:** Rejects — `0` is not a symbol-start char.
- **umol-edn C:** Accepts (matches Clojure).
- **Tests:** `test_s10_keyword_namespace_is_symbol`.

### S10d. The keyword part after : follows symbol first-character rules

`:0` is invalid under strict EDN — `0` cannot start a symbol.

- **Clojure:** `:0`, `:0foo` accepted (lenient).
- **umol-edn E:** Rejects.
- **umol-edn C:** Accepts (matches Clojure).
- **Tests:** `test_s10_keyword_first_char_restriction`.

Note: `:#foo` is rejected in both dialects — `#` is an interior-only char (S9c), not a valid symbol start.

### S10e. Post-slash first-character rules apply to keyword names

Same as S9e but for keywords. The name after `/` must start with a symbol-start
char. This means `#` and `:` are rejected as first char after slash (they are
interior-only per S9c), not just digits.

- **Clojure:** `:foo/#bar` and `:foo/:bar` accepted (lenient). `:foo/0bar` rejected.
- **umol-edn C:** Matches Clojure — digits rejected, `#`/`:` accepted.
- **umol-edn E:** All three rejected (spec: must be symbol-start).
- `#` and `:` as interior chars within the name part are fine in both dialects:
  `:foo/bar#baz`, `:foo/bar:baz`.
- **Tests:** `test_s10_keyword_post_slash_symbol_start`.

### S10f. Special symbol-start chars are valid after :

`.`, `+`, `-` are symbol-start chars, so `:.foo`, `:+foo`, `:-foo` are valid.

- **Tests:** `test_s10_keyword_special_start_chars`.

### S10g. :: auto-resolve keywords

- **Clojure:** `::foo` auto-resolves to `:current-ns/foo`. `::alias/name` resolves
  the alias via the namespace alias map.
- **umol-edn E:** Rejects `::` (not part of EDN spec).
- **umol-edn C:** Supported when `auto_resolve` config is provided. `::foo` →
  `:current_ns/foo`, `::alias/name` → `:resolved_ns/name`. Without config,
  `MissingAutoResolve` error. Unknown alias → `UnknownAlias` error.
- **Config:** `ParseConfig::auto_resolve: Option<AutoResolve>` with `current_ns: String`
  and `aliases: HashMap<String, String>`.
- **Tests:** `test_s10g_auto_resolve` (parser), `test_s10g_auto_resolve_streaming`,
  `test_s10g_unknown_alias`, `test_s10g_missing_config`, `test_s10g_in_map`,
  `test_s10b_double_colon_rejected` (Edn).

**Status:** Done.

---

## S11. Integers

### S11a. Optional sign

> Integers consist of the digits `0` - `9`, optionally prefixed by `-` to
> indicate a negative number, or (redundantly) by `+`.

- **Clojure:** `(read-string "+5")` => `5`.
- **umol-edn E/C:** Supported.
- **Tests:** `test_conformance_leading_plus_both_dialects`.

### S11b. No leading zeros

> No integer other than 0 may begin with 0.

- **Clojure:** `007` => 7, `00` => 0 (lenient, accepts leading zeros).
- **umol-edn E:** Rejects `007` and `00`. Follows spec.
- **umol-edn C:** Accepts (matches Clojure behavior).
- **Tests:** `test_s11b_leading_zeros_rejected_edn`, `test_s11b_leading_zeros_accepted_clojure`.

**Status:** Done.

### S11c. -0

> -0 is a valid integer not distinct from 0.

- **umol-edn E/C:** `Edn::Int(0)`.
- **Tests:** `test_conformance_negative_zero_int`.

### S11d. N suffix

> An integer can have the suffix `N` to indicate that arbitrary precision is
> desired.

- **umol-edn:** Error without `bignum` feature. Correct.
- **Tests:** `test_conformance_bigint_suffix_error`.

---

## S12. Floating point

> 64-bit (double) precision is expected.

- **umol-edn:** `f64`. Correct.

### S12a. M suffix

> a floating-point number may have the suffix `M` to indicate that exact
> precision is desired.

- **umol-edn:** Error without `bignum` feature.
- **Tests:** `test_conformance_bigdec_suffix_error`.

---

## S13. Lists

> A list is a sequence of values. Lists are represented by zero or more
> elements enclosed in parentheses `()`. Note that lists can be heterogeneous.

- **umol-edn E/C:** `Edn::List(Vec<Edn>)`.
- **Tests:** `test_conformance_seqs`.

**Status:** Done.

---

## S14. Vectors

> A vector is a sequence of values that supports random access. Vectors are
> represented by zero or more elements enclosed in square brackets `[]`. Note
> that vectors can be heterogeneous.

- **umol-edn E/C:** `Edn::Vector(Vec<Edn>)`.
- **Tests:** `test_conformance_seqs`.

**Status:** Done.

---

## S15. Maps

> A map is a collection of associations between keys and values. Maps are
> represented by zero or more key and value pairs enclosed in curly braces `{}`.
> Each key should appear at most once. No semantics should be associated with
> the order in which the pairs appear.

- **umol-edn E/C:** `Edn::Map(EdnMap)`. Duplicate keys: error by default, last-wins opt-in.
- **Tests:** `test_conformance_map_duplicate_key_error`, `test_conformance_map_duplicate_key_last_wins`, `test_conformance_map_contains_all_entries`.

> Note that keys and values can be elements of any type. The use of commas above
> is optional, as they are parsed as whitespace.

- **Tests:** `test_conformance_map_complex_keys`, `test_conformance_comma_in_collections`.

**Status:** Done.

---

## S16. Sets

> A set is a collection of unique values. Sets are represented by zero or more
> elements enclosed in curly braces preceded by `#` `#{}`. No semantics should
> be associated with the order in which the elements appear. Note that sets can
> be heterogeneous.

- **umol-edn E/C:** `Edn::Set(EdnSet)`.
- **Tests:** `test_conformance_set_contains_all_elements`.

**Status:** Done.

---

## S17. Tagged elements

> `#` followed immediately by a symbol starting with an alphabetic character
> indicates that symbol is a tag. A tag indicates the semantic interpretation
> of the following element.

- **umol-edn E/C:** `Edn::Tagged(tag, Box<Edn>)`.
- **Tests:** `test_conformance_tagged_basic`, `test_conformance_tagged_unqualified`.

### S17a. Tag without element is error

> Tags themselves are not elements. It is an error to have a tag without a
> corresponding tagged element.

- **umol-edn E/C:** Error. Correct.
- **Tests:** `test_conformance_tag_without_value_error`, `test_conformance_tag_at_eof_error`.

**Status:** Done.

### S17b. Unknown tag handling

> If a reader encounters a tag for which no handler is registered, the
> implementation can either report an error, call a designated 'unknown element'
> handler, or create a well-known generic representation that contains both the
> tag and the tagged element, as it sees fit.

- **umol-edn E/C:** If a `TagFn` is registered for the tag, it is called with
  the parsed value. Otherwise, creates `Tagged(tag, value)`. Valid per spec.
- **Tag reader registry:** `TagReaders` type in `config.rs`. Prepopulated with
  `#inst` (chrono) and `#uuid` (uuid) when features enabled. Custom readers via
  `TagReaders::insert()`.

### S17c. Reserved tags

> Tag symbols without a prefix are reserved by edn for built-ins. User tags
> must contain a prefix component.

- **Clojure:** No enforcement. Both `clj` and edamame ignore this rule.
- **umol-edn E:** Rejects bare (unqualified) tags unless registered in
  `TagReaders`. Built-in `inst`/`uuid` are pre-registered. Error: `InvalidTag`.
- **umol-edn C:** Accepts all tags (matches Clojure behavior).
- **Tests:** `test_s17c_bare_tag_rejected_edn`, `test_s17c_bare_tag_accepted_clojure`,
  `test_s17c_qualified_tag_accepted_both`, `test_s17c_inst_accepted_edn`,
  `test_s17c_uuid_accepted_edn`.

**Status:** Done.

---

## S18. Built-in tagged elements

> `#inst "rfc-3339-format"` — an instant in time.
> `#uuid "f81d4fae-..."` — a UUID.

- **umol-edn:** Feature-gated (`chrono`, `uuid`). When enabled, `TagReaders`
  pre-registers validation handlers. `read_inst` validates RFC 3339, `read_uuid`
  validates UUID format. Invalid values produce `InvalidInst`/`InvalidUuid`
  errors. Without features, stored as `Tagged` (no validation).
- **Conversion helpers:** `inst_to_edn()`, `uuid_to_edn()` in `tags.rs`.
- **Tests:** `test_s17c_inst_accepted_edn`, `test_s17c_uuid_accepted_edn`,
  `test_s17c_inst_invalid_rejected`, `test_s17c_uuid_invalid_rejected`,
  `test_s17c_inst_non_string_rejected`.

**Status:** Done.

---

## S19. Comments

> If a `;` character is encountered outside of a string, that character and all
> subsequent characters to the next newline should be ignored.

- **umol-edn E/C:** Supported.
- **Tests:** `test_conformance_comments`, `test_conformance_comment_at_eof`, `test_conformance_comment_inside_vector`.

**Status:** Done.

---

## S20. Discard

> `#` followed immediately by `_` is the discard sequence, indicating that the
> next element (whether separated from `#_` by whitespace or not) should be
> read and discarded. Note that the next element must still be a readable
> element. A reader should not call user-supplied tag handlers during the
> processing of the element to be discarded.

- **umol-edn E/C:** Supported in both dialects.
- **Tests:** `test_conformance_discard_nested`, `test_conformance_discard_in_map`,
  `test_conformance_discard_at_eof_error`, `test_conformance_discard_only_content_error`,
  `test_s20_discard_both_dialects`.

> The discard sequence is not an element. It is an error to have a discard
> sequence without a following element.

- **Tests:** `test_conformance_discard_at_eof_error`.

**Status:** Done.

---

## S21. Equality

> nil, booleans, strings, characters, and symbols are equal to values of the
> same type with the same edn representation.

- **umol-edn:** `PartialEq + Eq` derived/implemented for `Edn`. Correct.

> integers and floating point numbers should be considered equal to values only
> of the same magnitude, type, and precision.

- **umol-edn:** `Edn::Int` and `Edn::Float` are distinct variants, never equal to each other.

> sequences (lists and vectors) are equal to other sequences whose count of
> elements is the same, and for which each corresponding pair of elements (by
> ordinal) is equal.

- **umol-edn:** `Vec` equality. Correct.

> sets are equal if they have the same count of elements and, for every element
> in one set, an equal element is in the other.

- **umol-edn:** `HashSet` equality. Correct.

> maps are equal if they have the same number of entries, and for every
> key/value entry in one map an equal key is present and mapped to an equal
> value in the other.

- **umol-edn:** `HashMap` equality. Correct.

**Status:** Done. Tests implicit via `PartialEq` usage throughout.

---

---

## Ambiguity resolutions

These are decisions on underspecified areas of the EDN spec, documented in
`umol-edn/spec/edn-spec.md` section "Ambiguity resolutions".

### D1. Integer overflow

> Error without `bignum`; promote with `bignum`.

- **umol-edn E/C:** `i64::MAX + 1` produces `InvalidNumber`. No `bignum` feature yet.
- **Tests (E):** `test_s11_overflow`.
- **Tests (C):** `test_s11_overflow`.

**Status:** Done.

---

### D2. String escapes beyond the listed five

> `\uNNNN` in both dialects; `\b`, `\f`, octal in Clojure only.

- **umol-edn E:** Rejects `\b`, `\f`, `\uNNNN`, octal.
- **umol-edn C:** Accepts all.
- **Note:** The spec doc claims `\uNNNN` is accepted in both dialects, but the
  code currently rejects `\uNNNN` in Edn mode. This is a known spec-code discrepancy
  (see also S7).
- **Tests (E):** `test_s7_string_clojure_escapes_rejected` (rejects `\b \f \u \octal`).
- **Tests (C):** `test_s7_string_backspace_formfeed`, `test_s7_string_unicode_escape`,
  `test_s7_string_octal`.

**Status:** Tests cover current behavior. Spec doc discrepancy on `\uNNNN` to be resolved.

---

### D3. \formfeed, \backspace character literals

> Clojure only.

- **Tests (E):** `test_s8_clojure_named_rejected`.
- **Tests (C):** `test_s8_formfeed_backspace`.

**Status:** Done.

---

### D4. Duplicate map keys

> Configurable; error by default.

- **umol-edn E/C:** `DuplicateKeyPolicy::Error` (default) or `LastWins`.
- **Tests (E):** `test_s15_map_duplicate_key_error`, `test_s15_map_duplicate_key_last_wins`.
- **Tests (C):** `test_s15_map_duplicate_key_error`, `test_s15_map_duplicate_key_last_wins`.

**Status:** Done.

---

### D5. #_ nesting

> Each `#_` discards next form. `#_ #_ a b` discards both.

- **umol-edn E/C:** Supported in both dialects.
- **Tests (E):** `test_s20_discard_nested`, `test_s20_discard_nested_streaming`.
- **Tests (C):** `test_s20_discard_nested`, `test_s20_discard_streaming`.

**Status:** Done.

---

### D6. Whitespace definition

> Space, tab, newline, CR, comma. Not form feed.

- **Tests (E):** `test_s2_formfeed_not_whitespace`.
- **Tests (C):** `test_s2_formfeed_not_whitespace`.

**Status:** Done.

---

### D7. #_ + tagged literal

> Discards entire tagged form (e.g., `#_ #inst "2024"` discards the whole
> `#inst` form).

- **Tests (E):** `test_s20_discard_tagged_literal`.
- **Tests (C):** `test_s20_discard_tagged_literal`.

**Status:** Done.

---

### D8. Leading zeros

> Parsed as decimal, not octal.

Already covered: S11b. E: rejects, C: accepts (parsed as decimal).

**Status:** Done.

---

### D9. -0 / -0.0

> Valid.

- `-0` → `Edn::Int(0)`. `-0.0` → `Edn::Float(-0.0)` (sign-negative).
- **Tests (E):** `test_s11_negative_zero`, `test_s12_negative_zero_float`.
- **Tests (C):** `test_s11_negative_zero`, `test_s12_negative_zero_float`.

**Status:** Done.

---

### D10. / as symbol

> Valid.

Already covered: S9d. `/` alone is a legal symbol.

- **Tests (E):** `test_s9_slash_alone`.
- **Tests (C):** `test_s9_slash_alone`.

**Status:** Done.

---

### D11. #_ in Edn dialect

> Spec doc says: "In Edn mode, `#_` is parsed as a tagged literal with tag `_`."

**Note:** The code was updated (2026-04-01) to support `#_` as discard in both
dialects. The spec doc is outdated on this point. Current behavior: `#_` is
discard in both E and C.

- **Tests (E):** `test_s20_discard`, `test_s20_discard_streaming`.
- **Tests (C):** `test_s20_discard`, `test_s20_discard_nested`, etc.

**Status:** Done. Spec doc updated (section 3 dialect note).

---

## Clojure vs spec divergences

Where Clojure's reader diverges from the spec text, and what umol-edn does:

| Spec statement | Clojure behavior | umol-edn decision |
|---|---|---|
| `/` once only in symbols | Accepts `a/b/c` | E: rejects, C: accepts |
| `:/` not legal keyword | Accepts `:/` | Rejects (follows spec) |
| No leading zeros | Accepts `007` | E: rejects, C: accepts |
| `::` not legal keyword | Auto-resolves `::foo` to `:ns/foo` | E: rejects, C: auto-resolves (with config) |
| `:0`, `:0foo` digit-start kw | Accepts | E: rejects, C: accepts |
| `:0/foo` digit-start namespace | Accepts | E: rejects, C: accepts |
| `:foo/#bar`, `:foo/:bar` | Accepts | E: rejects, C: accepts |
| `foo/#bar`, `foo/:bar` symbols | Accepts | E: rejects, C: accepts |
| `\b`, `\f`, `\uNNNN` in strings | Accepts | E: rejects, C: accepts |
| `\formfeed`, `\backspace` chars | Accepts | E: rejects, C: accepts |

## Summary of resolved items

| ID | Issue | Status |
|---|---|---|
| S3 | Delimiters without whitespace | Done — tested |
| S4a | Unknown # dispatch errors | Done — tested |
| S4b | # is not a delimiter | Done — tested |
| S5–S6 | nil, booleans | Done |
| S9a | Symbol start chars | Done — tested |
| S9b | Sign/dot first char restriction | Done — fixed + tested |
| S9c | Interior chars (: # allowed) | Done — tested |
| S9d | Slash rules | Done — fixed + tested |
| S9e | Post-slash first-char restriction | Done — fixed + tested |
| S10a | Bare keywords valid | Done — tested |
| S10b | :/ and :/foo rejected | Done — tested |
| S10c | Namespace must be valid symbol | Done — tested |
| S10d | Keyword after : follows symbol-start rules | Done — tested |
| S10e | Post-slash symbol-start rules for keywords | Done — fixed + tested |
| S10f | Special start chars valid after : | Done — tested |
| S10g | :: auto-resolve keywords | Done — E rejects, C auto-resolves with config |
| S11a | Optional + sign | Done — tested |
| S11b | Leading zeros | Done — fixed + tested |
| S11c | -0 equals 0 | Done — tested |
| S13–S16 | Collections | Done — tested |
| S17a | Tag without element errors | Done — tested |
| S17b | Unknown tag handling | Done — TagReaders dispatch + Tagged fallback |
| S17c | Reserved tags | Done — E rejects bare unless registered; C accepts all |
| S18 | #inst and #uuid | Done — feature-gated validation handlers |
| S19 | Comments | Done — tested |
| S7 | String escape strict mode | Done — fixed + tested |
| S8 | Character literal strict mode | Done — tested |
| S20 | Discard | Done — both dialects |
| S21 | Equality | Done |

| D1 | Integer overflow | Done — both dialects |
| D2 | String escapes (\uNNNN discrepancy) | Spec doc discrepancy — \uNNNN rejected in E, spec says both |
| D3 | \formfeed, \backspace chars | Done |
| D4 | Duplicate map keys | Done — both dialects |
| D5 | #_ nesting | Done — both dialects |
| D6 | Whitespace: not form feed | Done — both dialects |
| D7 | #_ + tagged literal discard | Done — both dialects |
| D8 | Leading zeros | Done (see S11b) |
| D9 | -0 / -0.0 | Done — both dialects |
| D10 | / as symbol | Done — both dialects |
| D11 | #_ in Edn dialect | Done — spec doc updated |

## Summary of open items

| Item | What's needed |
|------|---------------|
| D2 | Resolve `\uNNNN` spec-code discrepancy (spec says both dialects, code rejects in Edn) |
