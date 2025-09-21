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
- LEX: LEX_INVALID_TOKEN, LEX_INTERTOKEN_WHITESPACE, LEX_BAD_PERCENT_FORM
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
- STYLE_BRACKET_ORGANIC: Prefer bare atom for organic subset when equivalent to bracketed form.
- STYLE_UNNECESSARY_PERCENT_RING_INDEX: Prefer single-digit ring numbers for 1..9 instead of “%01..%09”.
- STYLE_EXPLICIT_AROMATIC_BOND: Avoid explicit “:” when aromatic default applies.
- STYLE_REUSED_RING_INDICES: Discourage reusing the same ring digit within a connected component.

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

### Additional tasks

I am going over the list of rules in @35-opensmiles-missing-specs-2025-09-13.md . Please verify that each of these requirements has an associated error / warning and list them. Implement the TODOs included below. By "verify" I mean adding a test for this behavior, if it's not already available, and verifying that it passes. For each of the error / warning behaviors, make sure that there is at least one test emitting them. Ignore the points that are not listed below. Make a complete list of todos and execute them systematically. For implementations that require significant effort or design decisions, add to todo list and defer implementation. No need to provide intermediate reports until this list is completed.
1. Invalid atom symbol - error
4. Hydrogen has a hydrogen count - error
5. Hydrogen count > 9 or < 0 - error  TODO: verify that H0 is valid
H count exceed max number of implicit hydrogens - warning (see @element.rs , ELEMENT_DATA)
6. Charge < -15 or > +15 - error. TODO: consider if the limits make sense
charge > valence electrons - warning, charge outside of [min_charge, max_charge] range - warning (see @element.rs ELEMENT_DATA).
7, 9. Isotope mass number < 0 or > 999 - error. TODO: consider if the upper limit makes sense.
not a stable / metastable isotope - warning (see @isotope.rs , is_catalogued() function can be used).
8. Isotope 0 is valid  TODO: consider if we should make it invalid. How can it be usefully interpreted?
10. Invalid atom symbol (not organic subset or *) outside of brackets - error
12, 13. TODO: verify that [*] can have charge, chirality, H count, class.
15. Atom class < 0  - error
16. TODO: verify that class 0 is valid
17. Inconsistent bond symbols - error
18. Unclosed ring - error
19. Repeated ring indices - warning
20. TODO: verify that ring index 0 is valid
21. Invalid percent forms: %[^0-9] or %[0-9][^0-9] - error. TODO: verify that %123 is valid but parsed as "%12" + "3".
22. TODO: verify that ring closure C1CCC%01 is valid.
23. TODO: verify that multiple ring closures per atom are allowed, C1CC12CC2 as well as C1CC1%10CC%10, C%10CC%101CC1, C%10CC%10%11CC%11 .
24. (i) Self-bond - error, example: C11C 
(ii) Multiple bonds between atoms (chain + ring) - error, example: C1C1C
(iii) Multiple bonds between atoms (rings) - error, example: C12CCCCC12 
25. TODO: verify that c1ccc1 is valid
26. ":" in aromatic rings - warning
27. TODO: verify that c1ccccc1-c2ccccc2 is valid
28. (i) TODO: verify that [H+] is valid
(ii) TODO: verify that [H][H] is valid
(iii) TODO: verify that bridging H atoms are valid: [BH2]1[H][BH2][H]1
(iv) TODO: verify that [2H] and [3H] is valid
29. TODO: verify that [H][CH3] is valid
30. Dot before ring closure - error, example:  C.1CCCCC.1. TODO: consider if this rules makes a lot of sense.
32. TODO: verify that  C1.C1 is valid and equivalent to CC. 
(i) Leading dot - error, (ii) trailing dot - error, (iii) multiple dots - error
41. String too long TODO: consider if that is a reasonable limitation?
42. Class < 0 or > 9999 - error . TODO: consider if the upper limit makes sense?
43. TODO: verify that molecule with 100 rings can be parsed. Sorry ... :-)
44. TODO: verify that molecule with 100 nested branches can be parsed.
45. TODO: verify that atom can have 10 bonds.
46. TODO: verify that SMILES with 100 disconnected atoms can be parsed.
50. Bunch of warnings here
(i) Unnecessary bracket notation: [CH4] == C - warning. Brackets necessary if charge /= 0, atom class /= 0, chirality is given TODO: requires valence model, add to list but defer implementation.
(ii) Unnecessary charge index "+1" == "+", "-1" == "-" - warning
(iii) Unnecessary H count H1 = H - warning
(iv) Incorrect property ordering in brackets - warning
(v) Unnecessary explicit hydrogens [H][CH3] == C, [H][CH2+] == [CH3+] - warning
(vi) Explicit single bond "-" - warning (except between aromatic rings). Explicit aromatic bond ":" - warning
(va) Reuse of ring indices - c1ccccc1C1CCCC1 == c1ccccc1C2CCCC2 - warning
(via) Non-consecutive ring numbering C1CCC1CCC3CCC3 = C1CCC1CCC2CCC2 - warning,
First ring number not 1  C0CCCC0 == C1CCCC1  - warning.
(vii) Unnecessary ring closure on multiple bond C1CCC=1 == C1C=CC1 - perhaps warning. TODO: consider if this is too much effort to check
(viii) Multiple ring closures on the same atom (except for spiro atoms) - perhaps warning. TODO: consider if this is worth the effort. How hard is it to identify spiro atoms?
(ix) %0[0-9] ring index C%01CC%01 = C1CC1 - warning
(xiii) Unnecessary dot C1.C1C - perhaps warning. TODO: consider if this is worth the effort.
(xv) Unnecessary chiral markers - perhaps warning. TODO: add to todo list but defer.
 