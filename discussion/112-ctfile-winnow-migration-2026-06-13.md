# 112 — CTfile parser migration from nom to winnow

Status: Active · 2026-06-13

## Scope

Migrate the MOL/SDF CTfile parser in `umol-io/src/ctfile` from nom 8 to winnow 1.0.
This is the last production nom dependency after the molecule DSL parser moved to
winnow.

The migration should not be a mechanical replacement of combinator names. CTfile is a
line-oriented, fixed-column record format, and several current parsers are written in a
nom-specific style:

- parser factories return `impl Parser` closures;
- closures explicitly return `(remaining, output)`;
- block parsers manually preserve and replace immutable remaining slices;
- low-level parsers return `NomError<&[u8]>`, then block parsers erase those errors into
  line-level `ParseError` variants;
- nearly all parse failures are recoverable `nom::Err::Error`, including failures that
  occur after a record or field has already been identified.

The target should use winnow's mutable-input model directly, preserve the current
line-oriented block structure, and make recoverable versus committed failures explicit.

## Target parser API

Use `PResult` and `.parse_next()` throughout:

```rust
use winnow::error::ErrMode;
use winnow::{LocatingSlice, Parser};

type Input<'i> = LocatingSlice<&'i [u8]>;
type PResult<T> = Result<T, ErrMode<ParseError>>;

fn atom_input(input: &mut Input<'_>, flags: CtabParseFlags) -> PResult<(Atom, Point3D)> {
    // ...
}
```

Parser functions should normally accept `&mut Input` and return only their output.
Use a closure only where a combinator needs to capture configuration:

```rust
repeat(atom_count as usize, |input: &mut Input<'_>| {
    atom_input(input, flags)
})
.parse_next(input)
```

Do not preserve the nom pattern of returning an `impl Parser` closure solely to capture
flags or line offsets:

```rust
// Avoid as the default migration shape.
fn atom_input(flags: CtabParseFlags) -> impl Parser<...> {
    move |input| {
        // returns (remaining, output)
    }
}
```

At public API boundaries, construct the input stream and unwrap `ErrMode`:

```rust
fn parse_mol_bytes_to_table_ir_with(
    bytes: &[u8],
    config: &CtfileIoConfig,
) -> Result<Molecule, ParseError> {
    let mut input = Input::new(bytes);
    mol_input(&mut input, config).map_err(unwrap_err)
}
```

## `LocatingSlice`

`umol-edn` demonstrates the basic winnow location pattern. Its input is
`LocatingSlice<&str>`, and its `ParserError` implementation obtains an absolute offset
with `current_token_start()`. The molecule DSL implements `ParserError`, but does not use
`LocatingSlice`.

For CTfile, the useful location unit is normally one physical line. Construct a fresh
`LocatingSlice<&[u8]>` for each line:

```rust
fn parse_atom_line(
    line: &[u8],
    line_number: u32,
    flags: CtabParseFlags,
) -> Result<(Atom, Point3D), ParseError> {
    let mut input = Input::new(line);
    atom_input(&mut input, flags)
        .map_err(|error| promote_line_error(error, line_number, LineKind::Atom))
}
```

With this arrangement:

- `current_token_start()` is the column within the line;
- block parsers already know the physical line number;
- line-level errors can be promoted without pointer subtraction;
- the parser does not need a custom stream or a global byte-offset-to-line index.

A whole-file `LocatingSlice` would only provide an absolute byte offset. It is still
useful at top-level boundaries, but line and column reporting would require an additional
line index. The existing parser already identifies physical lines, so per-line locating
slices are the simpler fit.

One caveat: parsing a fixed-width field with `take(width)` produces the taken `&[u8]`.
Parsing that detached slice directly loses the parent `LocatingSlice` position. Fixed
width helpers should therefore remember the field-start column and promote child errors
at that position, or parse through a nested located input and add the field-start offset.

## Error model

### Current state

The CTfile `ParseError` currently contains nom-specific infrastructure:

- `Incomplete`;
- `NomError(NomErrorKind)`;
- `impl nom::error::ParseError`;
- `impl From<nom::Err<ParseError>>`;
- `from_nom` and the block-specific `*_from_nom` conversion helpers.

Low-level parsers generally return `NomError<&[u8]>`. Atom, bond, counts, property, and
other block parsers convert those errors into broad line-level variants, often discarding
the underlying reason.

### Recommended state

The `ParseError` API is in flux, so remove the nom-specific compatibility layer rather
than reproduce it with winnow names. Keep syntax location and semantic/build failures in
the domain error:

```rust
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum ParseError {
    #[error("{kind} syntax error at byte {offset}: {reason}")]
    Syntax {
        kind: LineKind,
        offset: usize,
        reason: SyntaxErrorKind,
    },

    #[error("invalid {kind} line at line {line}, col {column}: {reason}")]
    InvalidLine {
        kind: LineKind,
        line: u32,
        column: u32,
        reason: SyntaxErrorKind,
    },

    // Existing semantic/build errors remain domain-specific.
    PropertyMismatch(String),
    InconsistentSgroups(String),
    // ...
}
```

`ParserError<Input>` can construct a generic located syntax error:

```rust
impl<'i> ParserError<Input<'i>> for ParseError {
    type Inner = Self;

    fn from_input(input: &Input<'i>) -> Self {
        ParseError::Syntax {
            kind: LineKind::Unknown,
            offset: input.current_token_start(),
            reason: SyntaxErrorKind::UnexpectedInput,
        }
    }

    fn into_inner(self) -> Result<Self::Inner, Self> {
        Ok(self)
    }
}
```

Block boundaries should promote low-level syntax errors into line-aware errors. More
specific variants such as `MissingMEndTag`, `MissingDelimiter`, `UnexpectedEof`, and
semantic accumulator errors remain useful because they communicate format concepts
rather than parser-library concepts.

### `FromExternalError`

The DSL `ParseError` implements `ParserError<I>` but does **not** implement
`FromExternalError`. It manually creates `ErrMode::Backtrack` and `ErrMode::Cut` values.
`umol-edn` follows the same general approach.

Do not add a generic `FromExternalError` implementation solely to enable `.try_map()`.
That would tend to erase useful conversion failures into a generic syntax variant.
Prefer:

- `.verify_map()` for `Option`-returning conversions;
- explicit `.map_err(...)` for domain conversions;
- direct `ErrMode::Cut(ParseError::...)` or `ErrMode::Backtrack(ParseError::...)` where
  commitment is part of the grammar.

A targeted `FromExternalError<Input, SpecificDomainError>` may be worthwhile if a
specific conversion occurs repeatedly and retains all useful detail.

## Recoverable errors and `Cut`

The current CTfile parser creates no explicit `nom::Err::Failure` values. The only
`Failure` handling found is propagation/conversion code. Consequently, effectively every
CTfile syntax failure is currently recoverable.

In winnow:

- `ErrMode::Backtrack` means another grammar branch may own the same input;
- `ErrMode::Cut` means this parser has recognized the construct and failure is final;
- `alt`, `opt`, and repetition restore checkpoints only for `Backtrack`;
- a top-level `Backtrack` versus `Cut` may produce the same rejection, but the
  distinction matters as soon as the parser is nested.

Important CTfile commit points:

1. **Fixed-width fields.** Once a nonblank field is present, malformed contents are
   committed. `fixed_width_partial` already documents this as a fatal error, but the
   implementation currently returns recoverable errors.
2. **Counted blocks.** Once the counts line says an atom, bond, legacy atom-list, or
   property record is required, a malformed required record cannot be skipped.
3. **Recognized property tags.** After consuming `M  CHG`, `M  STY`, `V  `, and similar
   tags, failures in the corresponding body are committed.
4. **Length-prefixed entries.** Once an entry count is parsed, every required entry must
   parse. Use `length_repeat(count, cut_err(entry))` or an explicit loop that returns
   `Cut`.
5. **Optional fields with a separator.** If the separator is absent, the field is absent.
   If the separator is present, malformed contents are committed:

   ```rust
   opt(preceded(separator, cut_err(field))).parse_next(input)
   ```

6. **SDF records.** After recognizing `>`, a malformed field header or value is committed.
   After recognizing `$$$$`, a malformed delimiter line is committed.
7. **SGroup and RGroup expressions.** After recognizing an operator or range marker,
   failure to parse the right-hand side is committed.
8. **Semantic construction.** Molecule construction, property accumulation, duplicate
   detection, and consistency validation failures are committed.

Do not apply `Cut` before a branch has been identified. Unknown property tags, for
example, may remain recoverable while choosing among record kinds, but the body of a
known property tag must not backtrack into another record parser.

The DSL is useful precedent for explicit `ErrMode`, but it also underuses `Cut`. After a
known predicate prefix such as `#h` is consumed, malformed predicate bodies can remain
`Backtrack`, allowing enclosing repetition to stop and report trailing input instead.
The CTfile migration should establish commit points deliberately rather than copy that
behavior.

## Nom-to-winnow mapping

| Current nom combinator or pattern | Winnow mapping / recommendation |
| --- | --- |
| `nom::Parser::parse(input)` | `parser.parse_next(&mut input)` inside parsers; `.parse(...)` only at a complete-input boundary when appropriate |
| `IResult<&[u8], T, E>` | `PResult<T>` with mutable `Input` |
| closure returning `(remaining, output)` | direct `fn(&mut Input, ...) -> PResult<Output>` |
| `nom::Err::Error` | `ErrMode::Backtrack` |
| `nom::Err::Failure` | `ErrMode::Cut` |
| `nom::Err::Incomplete` | normally remove; CTfile uses complete input |
| `tag(...)` | `winnow::token::literal(...)` or byte/string literals implementing `Parser` |
| `tag_no_case(...)` | `winnow::ascii::Caseless(...)` with `literal`, or explicit normalized comparison |
| `take(n)` | `winnow::token::take(n)` |
| `take_while_m_n(min, max, p)` | `winnow::token::take_while(min..=max, p)` |
| `is_not(chars)` | `take_while(1.., \|c\| !chars.contains(c))` |
| `take_until1(tag)` | `winnow::token::take_until(1.., tag)` |
| `space0`, `multispace0`, integer parsers | corresponding `winnow::ascii` parsers; use `dec_int`/`dec_uint` or field-specific parsing |
| `map(parser, f)` | `parser.map(f)` |
| `value(value, parser)` | `parser.value(value)` |
| `map_opt(parser, f)` | `parser.verify_map(f)` |
| `map_res(parser, f)` | explicit `.map_err(...)`; use `.try_map(...)` only with a deliberate external-error conversion |
| `verify(parser, predicate)` | `parser.verify(predicate)` |
| `opt(parser)` | `opt(parser)`; add `cut_err` after an optional prefix has committed |
| `cond(condition, parser)` | `cond(condition, parser)` or clearer direct control flow |
| `alt((...))` | `alt((...))`; only use for genuine alternatives |
| `preceded`, `terminated`, `delimited`, `separated_pair` | same winnow combinators or parser tuple methods |
| parser tuples | parser tuples remain supported |
| `count(parser, n)` | `repeat(n, parser)` collected into the requested output |
| `length_count(count, parser)` | `length_repeat(count, parser)` |
| `separated_list1(separator, parser)` | `separated(1.., parser, separator)` |
| `success(value)` | `empty.value(value)` |
| `rest` | `winnow::token::rest` or read `input.as_ref()` where consumption is explicit |
| `all_consuming(parser)` | complete-input entry point or `terminated(parser, eof)` |
| `map_parser(outer, inner)` | `outer.and_then(inner)`; take care to preserve/promote location |
| manual property tag `match` | direct line-prefix dispatch or `dispatch!` after consuming the tag |
| manual error remapping by pointer subtraction | per-line `LocatingSlice` plus block-level promotion |

## Better winnow-native structure

### Keep block loops imperative

Atom, bond, property, and SDF blocks carry line numbers, counts, feature flags, and
cross-line state. Direct loops are clearer than forcing all block behavior into
`repeat_till` or nested combinators. The loop should mutate the input stream and call
line parsers with `.parse_next()`.

### Treat a physical line as the normal parsing unit

Counts, atom, bond, and most property records are fixed-column lines. Parse the line
boundary once, then parse a located line input. This makes end-of-line validation,
column reporting, and block-level error promotion explicit.

The optimized atom, bond, and counts parsers already use direct indexing and field
helpers rather than deeply nested combinators. That remains appropriate. Their migration
should replace nom error construction with located `PResult` failures, not rewrite them
into large declarative combinator expressions.

### Use deterministic record dispatch

Property records have unique prefixes. A parser written from scratch should inspect or
consume the prefix once and dispatch directly:

```rust
fn property_input(input: &mut Input<'_>, flags: CtabParseFlags) -> PResult<PropertyEntries> {
    let tag = take(6usize).parse_next(input)?;

    match tag {
        b"M  CHG" => cut_err(charge_entries).parse_next(input),
        b"M  STY" if flags.contains(CtabParseFlags::SGROUPS) => {
            cut_err(sgroup_type_entries).parse_next(input)
        }
        // ...
        _ => Err(ErrMode::Backtrack(ParseError::UnknownPropertyTag { /* ... */ })),
    }
}
```

`dispatch!` is also a good fit where its feature-flag guards remain readable. Either
form is preferable to repeatedly reparsing prefixes. Once a known tag has been consumed,
its body error must be `Cut`.

### Reserve `alt` for real grammar alternatives

`alt` remains appropriate for:

- SGroup multiplier expression forms;
- RGroup occurrence forms;
- SDF field versus delimiter at the item boundary;
- other cases where multiple parsers genuinely may own the same starting input.

It is less useful for deterministic fixed-column records whose prefix or current block
already identifies the parser.

## SDF data block

### Why a choice parser is needed

After the CTAB block and `M  END`, `sdf_data_block` is always parsed. The block contains
zero or more data fields followed by the required record delimiter:

```text
(data-field)* "$$$$"
```

The SDF property fields are optional/repeated; the entire SDF data block and its
delimiter are not optional. At each iteration, the next record is one of:

1. a data field beginning with `>`;
2. the terminal delimiter beginning with `$$$$`.

The current nom implementation first tries `sdf_data_field`, then tries
`sdf_delimiter`. This works because nom parsers receive immutable input and return a new
remaining slice. A failed field parser cannot mutate the original `remaining` variable.

With winnow, naively calling both with `.parse_next(input)` is incorrect: a failed field
parser may have advanced the mutable input. `alt` checkpoints the input, restores it on
`Backtrack`, and stops immediately on `Cut`.

`parse_peek()` would also preserve the caller's input by parsing a copy and returning
the remaining input, but it is mainly useful at migration/testing or API boundaries.
Inside a parser, `alt(...).parse_next(input)` expresses the grammar more directly.

### `alt` example

```rust
enum SdfItem {
    Field { name: String, value: String },
    Delimiter,
}

fn sdf_item(input: &mut Input<'_>, line: u32) -> PResult<SdfItem> {
    alt((
        preceded(
            literal(b">"),
            cut_err(|input: &mut Input<'_>| sdf_field_after_marker(input, line)),
        )
        .map(|(name, value)| SdfItem::Field { name, value }),

        preceded(literal(b"$$$$"), cut_err(sdf_delimiter_tail))
            .value(SdfItem::Delimiter),
    ))
    .parse_next(input)
}
```

The first `>` commits to an SDF field. A malformed field header must not reset and then
be reported as a missing delimiter. The delimiter branch remains available only while
the field branch has not claimed the input.

The delimiter parser should validate the complete delimiter record. Once its prefix is
recognized, malformed trailing content is also committed.

### Preferred record-oriented implementation

Because SDF records are line-oriented and their first byte identifies their role, a
direct line-dispatch loop may be clearer than speculative parsing:

```rust
loop {
    let line = next_line(input)?;

    match line.as_ref() {
        bytes if bytes.starts_with(b">") => {
            let (name, value) = sdf_field_after_header(line, input)?;
            data.insert(name, value);
        }
        b"$$$$" => return Ok(data),
        _ => {
            return Err(ErrMode::Cut(ParseError::InvalidSdfDataHeader {
                line: line_number,
            }));
        }
    }
}
```

This avoids speculative consumption entirely. `alt` is still the correct winnow
translation of the current grammar, but deterministic line dispatch is likely the
better from-scratch design.

## Detailed implementation plan

This migration should be implemented as one coordinated change set. The parser modules,
their shared helpers, `ParseError`, tests, and benchmarks are tightly coupled through
nom's input and error types. Preserving a compiling dual-parser state at every
intermediate step would require temporary adapters that add work and obscure the target
design.

The sequence below is an implementation order within that one change set, not a series
of independently compiling commits. Intermediate states may be broken. The required
compilation boundary is the completed migration.

### 0. Record the baseline and migration boundary

Before changing code:

1. Run the existing `umol-io` unit and integration tests and record the baseline.
2. Record the existing MOL/SDF benchmark results used to detect a material regression.
3. Inventory every nom reference under `umol-io/src/ctfile`, its tests, and the exported
   parsers used by benchmarks.
4. Treat changes outside CTfile parser plumbing, CTfile errors, tests, benchmarks, and
   `umol-io/Cargo.toml` as out of scope unless compilation requires them.
5. Preserve accepted input behavior by default. Intentional changes should be limited to
   improved error classification, location reporting, and committed-failure behavior.

### 1. Establish shared winnow parser infrastructure

Add the final shared types before migrating individual modules:

```rust
type Input<'i> = LocatingSlice<&'i [u8]>;
type PResult<T> = Result<T, ErrMode<ParseError>>;
```

Then add:

1. An `unwrap_err(ErrMode<ParseError>) -> ParseError` public-boundary helper.
2. Small constructors/helpers for domain `Backtrack` and `Cut` errors where they improve
   readability.
3. A line reader that advances a whole-file `Input` and returns the physical line without
   its terminator while preserving LF/CRLF byte consumption.
4. A helper that creates a fresh line-local `Input` from a physical line.
5. A line parser boundary that:
   - calls the line parser with `.parse_next()`;
   - accepts only permitted trailing whitespace;
   - promotes the line-local offset to a line-and-column `ParseError`;
   - preserves whether the child failure was `Backtrack` or `Cut`.
6. A consistent line-number convention. Retain the current external convention unless
   intentionally changing all errors and tests together.

Both whole-file and line-local streams can use `Input`; a whole-file stream advances
across records, while a newly constructed line-local stream resets location so
`current_token_start()` is a column.

### 2. Replace the nom-specific `ParseError` layer

Change `umol-io/src/ctfile/error.rs` before rewriting parser failures:

1. Add the final `LineKind` and `SyntaxErrorKind` types, or equivalent structured
   variants.
2. Implement `ParserError<Input<'_>> for ParseError`.
3. Add block-level promotion from line-local syntax errors to line-and-column errors.
4. Preserve domain-specific variants such as `UnexpectedEof`, `MissingMEndTag`,
   `MissingDelimiter`, property consistency errors, and SGroup errors.
5. Remove `Incomplete`, `NomError`, `NomErrorKind`, `NomParseError`, `From<nom::Err<_>>`,
   `from_nom`, and all block-specific `*_from_nom` helpers.
6. Do not implement generic `FromExternalError`.
7. Decide which errors are branch-selection failures and which are committed errors at
   their construction sites; do not defer all commitment decisions to block wrappers.

At the end of this step, many existing parsers will be broken because they still produce
nom errors. That is expected.

### 3. Rebuild primitive and fixed-width helpers

Migrate `parser/utils.rs` first because every line parser depends on it:

1. Replace `IntParser::nom_parser` with a winnow-aware integer-field abstraction or
   explicit typed field parsers.
2. Replace `parse_int_opt` with a parser/helper that:
   - recognizes an all-whitespace field as `None`;
   - rejects trailing garbage;
   - reports the field-start column;
   - returns `Cut` once nonblank field contents are present but malformed.
3. Replace `parse_float_f10_4` and `fixed_width_float_f10_4`, retaining Fortran `Fw.d`
   semantics and `fast_float2`.
4. Replace `fixed_width_partial`, `fixed_width_opt`, integer-range, minus-one, string,
   element, and unused-field helpers.
5. Make truncated-field behavior explicit and preserve the existing `partial_ok`
   semantics.
6. Ensure nested parsing of a taken field preserves or reconstructs the field-start
   column.
7. Migrate `rgroup_occurrence` and `rgroup_occurrences`; cut after `-`, `>`, or `<` has
   identified an occurrence form.
8. Retain `LinesWithOffset` only if the new advancing line reader still needs it;
   otherwise remove it rather than maintaining two line-consumption models.

Test the helpers conceptually while writing them, but migrate their test module later in
the same change set after the final error API is stable.

### 4. Migrate fixed-column line parsers

Migrate the parsers whose record type is already known from block context:

1. `counts_input`
2. `atom_input` and `extended_atom_input`
3. `bond_input` and `extended_bond_input`
4. `legacy_atom_list_input`

For each parser:

1. Change the signature to `fn(&mut Input, flags...) -> PResult<Output>`.
2. Keep direct indexing where it remains clearer and faster than combinators.
3. Read `input.as_ref()` for fixed-column validation, then advance the mutable input by
   the exact consumed width rather than returning `&input[offset..]`.
4. Replace `NomErrorKind::{Eof, Digit, Verify, Char, MapRes, Tag}` with structured domain
   syntax reasons.
5. Return `Cut` for malformed required fields because block context already establishes
   the record type.
6. Leave optional/truncated tail fields optional only where the CTfile rules and flags
   allow them.
7. Preserve the distinction between basic and extended atom/bond parsers.

The goal is not to force these optimized fixed-column parsers into nested combinators.
The goal is to make their input advancement, locations, and failure modes winnow-native.

### 5. Migrate counts, atom, bond, and legacy block parsers

Rewrite the corresponding block parsers around the advancing line helper:

1. Consume the counts line and parse it through a line-local `Input`.
2. For atom, bond, and legacy atom-list blocks, iterate exactly the number of records
   declared by the counts line.
3. Return committed `UnexpectedEof` when a counted record is missing.
4. Promote line parser errors with the correct physical line and local column.
5. Preserve position collection, `IGNORE_POSITIONS`, all-zero position handling, and
   basic/extended output types.
6. Remove immutable `remaining`/`byte_offset` bookkeeping where mutable input advancement
   now provides it.

After the counts line has been accepted, malformed required records must never be
interpreted as the end of a block or as another record kind.

### 6. Migrate RGroup and SGroup auxiliary grammars

Migrate `parser/rgroup.rs` and `parser/sgroup.rs` before property entries:

1. Convert all parsers to direct `PResult` functions.
2. Replace `map_res` enum/tag conversions with explicit domain errors or
   `.verify_map()` where lossless.
3. Use `alt` only for genuine expression alternatives.
4. In SGroup multiplier expressions, cut after an arithmetic operator is consumed.
5. In RGroup ranges and comparisons, cut after the range/comparison marker is consumed.
6. Preserve default values only when the complete optional construct is absent, not when
   its recognized prefix has malformed contents.

### 7. Migrate property entry-body parsers

Migrate the individual entry parsers in `parser/properties.rs` before the property
dispatcher:

1. Convert simple mappings to parser `.map(...)` or explicit direct functions.
2. Replace `length_count` with `length_repeat` or an explicit counted loop.
3. Cut every required entry after its count has been parsed.
4. Replace `count` with `repeat` where collection through a combinator is clearer.
5. For optional values introduced by a space or other separator, use
   `opt(preceded(separator, cut_err(value)))`.
6. Replace `map_parser` with `.and_then(...)` only where nested parsing remains clearer;
   otherwise parse fixed-width fields directly.
7. Preserve feature-flag gates and produce a domain error when a recognized but disabled
   property is encountered.
8. Preserve the special two-line `A  ` and `G  ` property forms, but convert their first
   line parsers to the common line-local error model.

This is the largest mechanical portion of the migration. Complete it before changing
property dispatch so the dispatcher can target final direct functions.

### 8. Replace property dispatch and property block parsing

Rewrite `property_input`, `extended_property_input`, `properties_block`, and
`extended_properties_block`:

1. Inspect/consume the `V  ` or six-byte `M  ...` tag once.
2. Dispatch with a direct `match` or `dispatch!`; do not reproduce nested
   `.parse(remaining).map(|(i, o)| ...)` calls.
3. Keep an unknown/unrecognized record tag as `Backtrack` only while selecting a record
   kind.
4. Wrap every recognized property body in `cut_err`.
5. Treat a known but disabled extension tag as a committed property error rather than
   allowing another parser to reinterpret it.
6. Keep `M  END`, `A  `, and `G  ` handling at block level.
7. For two-line properties, commit after the first-line prefix and return
   `UnexpectedEof` if the second line is missing.
8. Preserve `NO_V2000_END_TAGS`; otherwise require and commit to `M  END`.
9. Promote all body errors to the current physical property line and local column.

### 9. Rewrite SDF data parsing as record-oriented parsing

Prefer deterministic line dispatch over a literal translation of the current
try-field-then-delimiter loop:

1. At the SDF item boundary, inspect the next physical line without losing the ability to
   report its location.
2. If it begins with `>`, commit to a data field and parse the complete header.
3. Parse value lines until the required blank separator, preserving the current
   multi-line value joining behavior.
4. If the line is exactly a valid `$$$$` delimiter record, finish the molecule.
5. If input ends before `$$$$`, return committed `MissingDelimiter`/`UnexpectedEof`.
6. If a line begins like a delimiter but is malformed, return a committed delimiter
   error.
7. Otherwise return a committed invalid SDF data-record/header error.
8. Preserve zero-field SDF records: an immediate valid `$$$$` is accepted.
9. Preserve line offsets across multiple SDF molecules and LF/CRLF inputs.

If `alt` is retained at the SDF item boundary instead, the field branch must become
`Cut` immediately after `>`, and the delimiter branch must become `Cut` immediately
after its identifying prefix. Do not use `parse_peek()` as routine internal control
flow.

### 10. Migrate CTAB composition and public MOL/SDF entry points

Rewrite `parser.rs` after all child parsers have final signatures:

1. Convert `ctab_block` and `extended_ctab_block` to direct `PResult` functions.
2. Call counts, atom, bond, legacy-property, and property blocks with `.parse_next()`.
3. Convert unsupported-feature and semantic molecule-build failures to `Cut`.
4. Convert single-MOL entry points to construct a whole-file `Input`, parse one record,
   consume permitted trailing whitespace, and reject other trailing input.
5. Convert SDF entry points to loop over the mutable input until permitted trailing
   whitespace/end-of-file.
6. Preserve Unicode whitespace normalization before constructing the input.
7. Preserve comments, source format, chirality frame, positions, properties, and
   basic/extended output behavior.
8. Unwrap `ErrMode` only at public `Result<_, ParseError>` boundaries.

### 11. Migrate tests and benchmarks in the same change set

Update all CTfile parser tests after the final parser and error APIs exist:

1. Replace nom imports, `IResult`, `Finish`, `.parse(...)`, and tuple
   `(remaining, output)` assertions.
2. Add a small test helper that constructs `Input`, calls `.parse_next()`, and exposes the
   remaining slice when a test needs to assert consumption.
3. Replace assertions on `NomErrorKind` and `nom::Err` with assertions on domain
   `ParseError`, `ErrMode::Backtrack`, or `ErrMode::Cut`.
4. Keep existing accepted/rejected fixture coverage unchanged unless the new commitment
   behavior intentionally improves the reported error.
5. Add focused tests for every commit point listed in this document.
6. Update benchmark-facing exported parsers to accept/construct winnow input without
   adding nom-compatibility wrappers.
7. Compare parser benchmark results with the recorded baseline and investigate material
   regressions, especially in atom/bond fixed-column parsing.

### 12. Remove nom and compatibility remnants

Once production code, tests, and benchmarks use winnow:

1. Remove nom from `umol-io/Cargo.toml`.
2. Remove all nom imports and aliases.
3. Remove obsolete conversion helpers, compatibility wrappers, and comments describing
   nom behavior.
4. Search all of `umol-io` for `nom`, `IResult`, `NomError`, `NomErrorKind`, `.parse(`,
   `Finish`, and tuple-return parser signatures.
5. Confirm that any remaining `.parse(...)` calls are intentional complete-input winnow
   boundaries rather than missed migrations.

### 13. Final verification

The completed one-block migration is ready only after:

1. Formatting succeeds.
2. `cargo check -p umol-io` succeeds.
3. All `umol-io` unit and integration tests pass.
4. Workspace checks/tests covering `umol-io` consumers pass.
5. Clippy reports no new warnings in the changed crates.
6. No nom references or dependency remain in `umol-io`.
7. MOL/SDF conformance counts match the baseline.
8. Targeted error tests demonstrate the intended `Backtrack`/`Cut` behavior and correct
   line/column reporting.
9. Benchmark comparison shows no unexplained material regression.

## Verification focus

The existing MOL/SDF conformance suites should remain the primary behavioral guard.
Add targeted error tests for the newly explicit commit behavior:

- nonblank malformed fixed-width fields do not become absent/default fields;
- malformed bodies after known property tags report that property line;
- counted blocks cannot stop early or reinterpret malformed required records;
- optional fields backtrack only when their prefix/separator is absent;
- malformed SDF headers beginning with `>` do not become missing-delimiter errors;
- zero-field SDF records accept an immediate valid `$$$$`;
- malformed delimiters are committed delimiter errors;
- line and column reporting remains correct for LF, CRLF, and truncated lines.

## Decision summary

- Use `PResult<T>` and `.parse_next()` as the normal internal API.
- Use per-line `LocatingSlice<&[u8]>` so locations naturally report columns.
- Remove nom-specific error variants and conversions instead of reproducing them.
- Do not add a generic `FromExternalError`; preserve domain conversion details manually.
- Keep stateful block parsing imperative and line-oriented.
- Use deterministic prefix dispatch for CTfile records.
- Use `alt` for genuine alternatives, with `Cut` immediately after a branch is
  recognized.
- In SDF, data fields are optional/repeated, but the terminal delimiter is required.
