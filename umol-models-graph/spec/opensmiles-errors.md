### Diagnostics registry (initial)

This file defines the initial diagnostics taxonomy and stable codes for OpenSMILES parsing and linting in `umol-models-graph`.

#### Taxonomy

- Severity: Error | Warning
- Category: LEX (lexing), SYN (syntax/grammar), RING (ring rules), NUM (numeric constraints), BRKT (bracket fields), STEREO (stereochemistry), STYLE (style/normalization), INTERNAL (fallback)
- Code: stable UPPER_SNAKE identifiers (no numeric prefixes); messages are short and deterministic; parameters are carried as named fields
- Span: byte range [start, end) over the original input

#### Errors

- LEX_INVALID_TOKEN (LEX, Error): input slice cannot be tokenized per lexical rules
- LEX_BAD_PERCENT_FORM (LEX, Error): percent form not followed by two digits
- LEX_BAD_PERCENT_RANGE (LEX, Error): percent form out of accepted range

- SYN_UNEXPECTED_TOKEN (SYN, Error): parser encountered a token not valid in the current production
- SYN_TRAILING_BOND (SYN, Error): bond symbol found at end of input or before terminator
- SYN_UNCLOSED_BRACKET (SYN, Error): bracket atom not properly closed
- SYN_DOT_BEFORE_RING (SYN, Error): dot placed before a ring index

- RING_UNCLOSED (RING, Error): ring index opened but not closed by end of component/molecule
- RING_CONFLICT_DIR (RING, Error): conflicting up/down directions on the same ring closure
- RING_SELF_LOOP (RING, Error): ring closure creates a self-loop
- RING_TWO_MEMBER (RING, Error): ring closure creates a two‑member ring

Emitted by parser/state during ring processing; codes may appear in `lint_smiles_parse` reports even when lexing succeeds.

- NUM_OVERFLOW (NUM, Error): numeric literal exceeds supported bounds
- NUM_CLASS_NEGATIVE (NUM, Error): negative atom class is not permitted
- NUM_HCOUNT_BAD (NUM, Error): hydrogen count invalid (e.g., >9)
- NUM_CHIRAL_OUT_OF_RANGE (NUM, Error): chirality parameter out of accepted range

- BRKT_DUP_FIELD (BRKT, Error): duplicate bracket field of the same kind
- BRKT_HCOUNT_TWO_DIGITS (BRKT, Error): H-count with two or more digits
- BRKT_EMPTY_CLASS (BRKT, Error): class field without a numeric value

- STEREO_DOUBLE_CONFLICT (STEREO, Error): conflicting cis/trans specifications
- STEREO_DOUBLE_INSUFFICIENT (STEREO, Error): insufficient markers to define double-bond stereo

- INTERNAL_PARSER_STATE (INTERNAL, Error): unexpected internal parser state

Notes:
- R50(iv) order is accepted by grammar; misordering is a warning (STYLE_BRKT_ORDER), not an error

#### Warnings (style)

- STYLE_BRKT_ORDER (STYLE, Warning): prefer [chirality][H][charge][class] ordering in bracket atoms
- STYLE_CHARGE_SIGN_SIMPLE (STYLE, Warning): prefer [+]/[-] over [+1]/[-1]
- STYLE_HCOUNT_ONE_SIMPLE (STYLE, Warning): prefer H over H1 in bracket H-count
- STYLE_BARE_ORGANIC (STYLE, Warning): prefer bare organic subset atom over equivalent bracketed form
- STYLE_SINGLE_DIGIT_RING (STYLE, Warning): prefer single‑digit ring numbers for 1..9 instead of %01..%09
- STYLE_AVOID_AROMATIC_BOND_SYMBOL (STYLE, Warning): avoid explicit ':' when aromatic default applies
- STYLE_NO_REUSE_RING_DIGITS (STYLE, Warning): avoid reusing the same ring digit in a connected component


