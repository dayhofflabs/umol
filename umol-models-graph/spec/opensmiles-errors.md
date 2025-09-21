### Diagnostics registry

This file defines the diagnostics taxonomy and stable codes for OpenSMILES parsing and linting in `umol-models-graph`.

#### Taxonomy

- Severity: Error | Warning
- Category: LEX (lexing), BRANCH (branch structure), RING (ring rules), NUM (numeric constraints), BRKT (bracket fields), STEREO (stereochemistry), STYLE (style/normalization), SYN (parser-only fallback), INTERNAL (fallback)
- Code: stable UPPER_SNAKE identifiers (no numeric prefixes); messages are short and deterministic; parameters are carried as named fields
- Span: byte range [start, end) over the original input

#### Errors

- LEX_INVALID_TOKEN (LEX, Error): input slice cannot be tokenized per lexical rules
- LEX_BAD_PERCENT_FORM (LEX, Error): percent form not followed by two digits
- LEX_TRAILING_BOND (LEX, Error): bond symbol found at end of input or before terminator
- LEX_DOT_BEFORE_RING (LEX, Error): dot placed before a ring index
- LEX_LEADING_DOT (LEX, Error): dot at start of input or component
- LEX_TRAILING_DOT (LEX, Error): dot at end of input or component
- LEX_MULTIPLE_DOTS (LEX, Error): consecutive dots without an intervening component

- BRKT_UNCLOSED (BRKT, Error): bracket atom not properly closed
- BRKT_UNEXPECTED_CLOSE (BRKT, Error): ']' outside of a bracket atom
- BRKT_FIELD_OUTSIDE (BRKT, Error): bracket-only field (e.g., @, +, ++) outside of a bracket atom
- BRKT_DUP_FIELD (BRKT, Error): duplicate bracket field of the same kind
- BRKT_HCOUNT_TWO_DIGITS (BRKT, Error): H-count with two or more digits
- BRKT_EMPTY_CLASS (BRKT, Error): class field without a numeric value
- BRKT_H_ON_H (BRKT, Error): hydrogen element carries an H-count field

- BRCH_UNEXPECTED_CLOSE (BRANCH, Error): ')' without a matching '('
- BRCH_UNCLOSED (BRANCH, Error): open '(' not closed before component/end
- BRCH_DANGLING_BOND (BRANCH, Error): bond before ')' or component end
- BRCH_EMPTY_BRANCH (BRANCH, Error): empty branch (e.g., '()')

- RING_UNCLOSED (RING, Error): ring index opened but not closed by end of component/molecule
- RING_CONFLICT_DIR (RING, Error): conflicting up/down directions on the same ring closure
- RING_SELF_LOOP (RING, Error): ring closure creates a self-loop
- RING_TWO_MEMBER (RING, Error): ring closure creates a two‑member ring

- SYN_UNEXPECTED_TOKEN (SYN, Error): parser encountered a token not valid in the current production

Emitted by parser/state during ring processing; codes may appear in `lint_smiles_parse` reports even when lexing succeeds.

- NUM_OVERFLOW (NUM, Error): numeric literal exceeds supported bounds
- NUM_CLASS_NEGATIVE (NUM, Error): negative atom class is not permitted
- NUM_CLASS_TOO_LARGE (NUM, Error): atom class exceeds 4 digits (max 9999)
- NUM_HCOUNT_OUT_OF_RANGE (NUM, Error): hydrogen count invalid (e.g., >9)
- NUM_CHARGE_OUT_OF_RANGE (NUM, Error): absolute charge exceeds hard limit (|q| > 15)
- NUM_ISOTOPE_TOO_LARGE (NUM, Error): isotope mass number exceeds 999
- NUM_CHIRAL_OUT_OF_RANGE (NUM, Error): chirality parameter out of accepted range

- STEREO_DOUBLE_CONFLICT (STEREO, Error): conflicting cis/trans specifications
- STEREO_DOUBLE_INSUFFICIENT (STEREO, Error): insufficient markers to define double-bond stereo

- INTERNAL_PARSER_STATE (INTERNAL, Error): unexpected internal parser state

#### Warnings (style)

- STYLE_BRKT_ORDER (STYLE, Warning): prefer [chirality][H][charge][class] ordering in bracket atoms
- STYLE_CHARGE_SIGN_SIMPLE (STYLE, Warning): prefer [+]/[-] over [+1]/[-1]
- STYLE_HCOUNT_ONE_SIMPLE (STYLE, Warning): prefer H over H1 in bracket H-count
- STYLE_BRKT_ORGANIC (STYLE, Warning): prefer bare organic subset atom over equivalent bracketed form
- STYLE_UNNECESSARY_PERCENT_RING_INDEX (STYLE, Warning): prefer single‑digit ring numbers for 1..9 instead of %01..%09
- STYLE_EXPLICIT_AROMATIC_BOND (STYLE, Warning): avoid explicit ':' when aromatic default applies
- STYLE_REUSED_RING_INDICES (STYLE, Warning): avoid reusing the same ring digit in a connected component
- STYLE_EXPLICIT_SINGLE_BOND (STYLE, Warning): avoid explicit '-' when default applies
- STYLE_UNNECESSARY_EXPLICIT_H (STYLE, Warning): avoid explicit hydrogen when implicit is preferred
- STYLE_FIRST_RING_NOT_ONE (STYLE, Warning): prefer starting ring numbering at 1
- STYLE_NONCONSECUTIVE_RING_NUMBERING (STYLE, Warning): prefer consecutive ring numbers where possible
- NUM_HCOUNT_EXCEEDS_MAX_IMPLICIT (NUM, Warning): H-count exceeds element's max implicit hydrogens
- NUM_CHARGE_OUTSIDE_ELEMENT_RANGE (NUM, Warning): charge outside element-supported bounds
- NUM_CHARGE_EXCEEDS_VALENCE_ELECTRONS (NUM, Warning): charge exceeds valence electron count
- NUM_ISOTOPE_UNCATALOGUED (NUM, Warning): isotope is not catalogued (unstable or too short-lived)
