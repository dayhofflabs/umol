# Diagnostics for OpenSMILES strings

### Scope and goals
- Define a stable diagnostics taxonomy with severity (Error, Warning), category, and code.
- Enable linting of SMILES with both errors (blocking) and warnings (style/normalization).
- Instrument lexer, parser, and semantic passes to emit structured diagnostics with spans.

### Deliverables (initial)
- Spec file: `umol-models-graph/spec/opensmiles-errors.md` — canonical registry of codes (normative).
- Core types: `umol-models-graph/src/diagnostics.rs` — Diagnostic, Severity, Category, Code, Span, Details.
- Mappings:
  - Lexer: map tokenization failures and percent forms to diagnostics.
  - Parser: map rejects (unexpected token, trailing bond, bracket errors) to diagnostics.
  - Semantics (ring/state): map ring rules (unclosed, self-loop, 2-member, conflict) to diagnostics.
- Lint API: `lint_smiles(input)` returning a report with errors and warnings.
- Tests: table-driven asserts of code+span; a few cross-layer cases.

### Taxonomy shape
- Severity: Error | Warning
- Category: LEX (lexing), SYN (syntax/grammar), RING (ring rules), NUM (numeric constraints), BRKT (bracket fields), STEREO (stereochemistry), STYLE (style/normalization), INTERNAL (fallback).
- Code: stable UPPER_SNAKE strings tied to the spec, e.g., LEX_INVALID_TOKEN. No numeric prefixes initially; add later only if needed.
- Message: deterministic, short; parameters carried in Details (e.g., ring_index, expected, got).
- Span: byte start/end (inclusive, exclusive).

### Initial error codes (Errors)
- LEX: LEX_INVALID_TOKEN, LEX_INTERTOKEN_WHITESPACE, LEX_BAD_PERCENT_FORM, LEX_BAD_PERCENT_RANGE
- SYN: SYN_UNEXPECTED_TOKEN, SYN_TRAILING_BOND, SYN_UNCLOSED_BRACKET, SYN_DOT_BEFORE_RING
- RING: RING_UNCLOSED, RING_CONFLICT_DIR, RING_SELF_LOOP, RING_TWO_MEMBER
- NUM: NUM_OVERFLOW, NUM_CLASS_NEGATIVE, NUM_HCOUNT_BAD, NUM_CHIRAL_OUT_OF_RANGE
- BRKT: BRKT_DUP_FIELD, BRKT_HCOUNT_TWO_DIGITS, BRKT_EMPTY_CLASS
- STEREO: STEREO_DOUBLE_CONFLICT, STEREO_DOUBLE_INSUFFICIENT
- INTERNAL: INTERNAL_PARSER_STATE

Notes:
- R50(iv) is acceptance: bracket field order is no longer an error. Any order is accepted; we’ll lint as STYLE_BRKT_ORDER (warning), not BRKT_ORDER (error).

### Initial warning codes (Warnings; “Standard Form” 50)
- STYLE_BRKT_ORDER: Suggest “[chirality][H][charge][class]” ordering when brackets are used.
- STYLE_CHARGE_SIGN_SIMPLE: Prefer “[X+]” over “[X+1]” and “[X-]” over “[X-1]”.
- STYLE_HCOUNT_ONE_SIMPLE: Prefer “H” over “H1” in bracket H-count.
- STYLE_BARE_ORGANIC: Prefer bare atom for organic subset when equivalent to bracketed form.
- STYLE_SINGLE_DIGIT_RING: Prefer single-digit ring numbers for 1..9 instead of “%01..%09”.
- STYLE_AVOID_AROMATIC_BOND_SYMBOL: Avoid explicit “:” when aromatic default applies.
- STYLE_NO_REUSE_RING_DIGITS: Discourage reusing the same ring digit within a connected component.

Deferrals for later (require deeper semantics): avoid ring closures on multiple bonds, pick main chain, aromatic preference normalization.

### Mapping plan (incremental)
- Lexer (`io/smiles/lexer.rs`):
  - Convert Logos errors and bad percent lexemes into LEX_* diagnostics with spans from token slices.
- Parser (`io/smiles/parser/grammar.lalrpop`):
  - Wrap lalrpop `ParseError` with SYN_* codes; attach expected/got in Details.
  - Emit BRKT_* on duplicate fields, bad H count digits, empty class.
- Semantics (`io/smiles/state.rs`):
  - On finish/commit: RING_* diagnostics with ring_index, atom indices, and involved bond span if available.
- Lint pass (`io/smiles/linter.rs` or `io/smiles/mod.rs::linter::lint_smiles`):
  - Run parse (collect errors). If parse succeeds, run style checks for STYLE_*.
  - Return DiagnosticsReport { diagnostics, counts, has_errors }.

### Testing
- Deterministic table tests in existing suites:
  - Lexing: invalid token, bad percent; assert LEX_* code and exact span.
  - Parsing: trailing bond, unmatched bracket; assert SYN_* with span.
  - Rings: self-loop, 2-member, conflict; assert RING_* with ring_index.
  - Style: bracket order, charge ±1, single-digit ring preference; assert STYLE_* and suggested replacement in Details.
- Keep messages stable; assert against code and span; message formatting is secondary but deterministic.

### Spec integration
- `spec/opensmiles-errors.md`:
  - For each code: Code, Title, Category, Severity, Condition (normative), Parameters, Example (valid/invalid), Help.
- Cross-reference from `spec/opensmiles-umol.md` “Error Policy” to the registry file.

### Acceptance criteria for “initial”
- All listed error codes are emitted from at least one site; the corresponding tests pass.
- STYLE_BRKT_ORDER, STYLE_CHARGE_SIGN_SIMPLE, STYLE_HCOUNT_ONE_SIMPLE, STYLE_SINGLE_DIGIT_RING, STYLE_AVOID_AROMATIC_BOND_SYMBOL, STYLE_BARE_ORGANIC, STYLE_NO_REUSE_RING_DIGITS implemented and tested.
- Lint API returns stable codes and spans.

### Interaction with umol Error hierarchy

Before we proceed, I wanted to ask about the interaction between the error codes and the Error types in @error.rs . Is it necessary / useful to add new error types or adapt existing ones to simplify linting?

- Keep `umol/src/error.rs` as-is. Linting should use a separate, non-exceptional diagnostics layer (Diagnostic, Severity, Category, Code, Span, Details) in `umol-models-graph`.
- Provide two APIs:
  - parse_strict: returns Result<… , ParseError> (convert first Error-severity diagnostic to a `ParseError` for strict flows).
  - lint_smiles: returns DiagnosticsReport with errors and warnings (no exceptions).
- Optional bridge (only if you want tighter integration): add a helper like `DiagnosticsReport::to_parse_error()` or a new `ParseError::Diagnostics(Vec<Diagnostic>)`. I’d defer this per YAGNI and start with a simple adapter in the parse_strict wrapper.

This keeps error handling for operations separate from lint diagnostics, while allowing strict callers to translate diagnostics into existing `ParseError` when needed.

If this looks good, I’ll proceed to draft `spec/errors.md` with the initial registry and then define the core types and LEX/SYN mappings.

- Added a diagnostics TODO list and set “Create diagnostics registry spec” in progress. 
- Proposed a focused plan: taxonomy, initial error/warn codes, mapping points in `lexer.rs`, `grammar.lalrpop`, `state.rs`, a `lint_smiles` API, and table tests.