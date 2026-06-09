### Diagnostics registry

This file defines the diagnostics taxonomy and stable codes for OpenSMILES parsing and linting in `umol-models-graph`.

#### Taxonomy

- Severity: Error | Warning
- Category: LEX (lexing), SYN (parsing), TOPO (graph topology), VAL (valence checks), AROM (aromaticity),  STEREO (stereochemistry), NUM (numeric constraints), STYLE (style/normalization), INTERNAL (fallback)
- Code: stable UPPER_SNAKE identifiers (no numeric prefixes); messages are short and deterministic; parameters are carried as named fields
- Span: byte range [start, end) over the original input

#### Errors

##### Lexical / syntactic (during parse)

- LEX_LEADING_WHITESPACE (LEX, Error): whitespace at start of input before SMILES string
- LEX_INVALID_WHITESPACE (LEX, Error): inter-token whitespace encountered
- LEX_INVALID_ELEMENT (LEX, Error): invalid element symbol
- LEX_INVALID_TOKEN (LEX, Error): input slice cannot be tokenized per lexical rules

- SYN_UNBALANCED_OPEN_PAREN (SYN, Error): open `(` not closed before component/end
- SYN_UNBALANCED_CLOSE_PAREN (SYN, Error): `)` without a matching `(`
- SYN_EMPTY_BRANCH (SYN, Error): empty branch (e.g., `()`)
- SYN_EMPTY_GROUP (SYN, Error): empty group (e.g., `()`)
- SYN_NONFINAL_GROUP (SYN, Error): group closing `)` is not last before component/end

- SYN_LEADING_BOND (SYN, Error): bond token appears with no preceding atom (start of input or after a top-level group)
- SYN_TRAILING_BOND (SYN, Error): bond token found at end of input or component
- SYN_CONSECUTIVE_BONDS (SYN, Error): consecutive bond tokens without an intervening atom or ring index

- SYN_LEADING_RING (SYN, Error): ring index appears with no preceding atom (start of input or after a top-level group)
- SYN_UNBALANCED_RING_INDEX (SYN, Error): ring index opened but not closed by end of component/molecule
- SYN_INVALID_RING_INDEX (SYN, Error): percent ring form invalid (must be exactly two digits)
- SYN_MISMATCHED_RING_BOND_DIRS (SYN, Error): conflicting up/down bond directions on the same ring closure
- SYN_MISMATCHED_RING_BOND_ORDERS (SYN, Error): explicit bond orders on the two endpoints of a ring closure differ

- SYN_LEADING_DOT (SYN, Error): dot at start of input or component
- SYN_TRAILING_DOT (SYN, Error): dot at end of input or component
- SYN_CONSECUTIVE_DOTS (SYN, Error): consecutive dots without an intervening component
- SYN_DOT_BEFORE_RING (SYN, Error): dot placed before a ring index

- SYN_UNBALANCED_OPEN_BRACKET (SYN, Error): `[` without a matching closing `]` (reported by parser)
- SYN_UNBALANCED_CLOSE_BRACKET (SYN, Error): `]` without a matching `[` (reported by parser)
- SYN_STRAY_BRACKET_FIELD (SYN, Error): bracket-only field (e.g., `@`, `+`, `++`) outside of a bracket atom
- SYN_DUPLICATE_BRACKET_FIELD (SYN, Error): duplicate bracket field of the same kind
- SYN_MISSING_CLASS_INDEX (SYN, Error): class field without numeric value
- SYN_MISSING_CHIRALITY_INDEX (SYN, Error): chirality marker (`@TH`, `@AL`, `@SP`, `@TB`, `@OH`) without numeric value
- SYN_EMPTY_BRACKET (SYN, Error): no bracket fields
- SYN_BRACKET_H_WITH_HCOUNT (SYN, Error): hydrogen element carries an H-count field
- SYN_INVALID_BRACKET (SYN, Error): bracket errors (e.g., no element in [1])

##### Semantic (post-parse)

- TOPO_SELF_LOOP_RING (TOPO, Error): self-loop edge present
- TOPO_PARALLEL_EDGES (TOPO, Error): more than one edge between the same atom pair

- VAL_OUT_OF_ELEMENT_RANGE (VALENCE, Error): specified valence outside of the element range
- VAL_HCOUNT_OUT_OF_ELEMENT_RANGE (VALENCE, Error): specified H-count outside of the element range
- VAL_CHARGE_OUT_OF_ELEMENT_RANGE (VALENCE, Error): specified charge outside of the element range
- VAL_HCOUNT_MISMATCH (VALENCE, Error): explicit bracket H-count conflicts with valence-based implicit hydrogen count
- VAL_NO_MATCH (VALENCE, Warning): no valence pattern matched the atom (element, bond sum, charge, etc.)
- VAL_AMBIGUOUS_MATCH (VALENCE, Warning): multiple valence patterns matched; selected the most specific

- AROM_ATOM_NOT_IN_RING (AROM, Error): aromatic atom not in ring
- AROM_BOND_NOT_IN_RING (AROM, Error): aromatic bond symbol not in ring
- AROM_NO_MATCHING_AROMATIC_ATOM_CONFIG (AROM, Error): atom must support aromatic valences
- AROM_INVALID_AROMATIC_ATOM (AROM, Error): atom ineligible for aromatic systems
- AROM_INVALID_AROMATIC_BOND_ATOM (AROM, Error): aromatic bond between atoms that are ineligible for aromatic systems
- AROM_BOND_ORDER_MISMATCH (AROM, Error): Explicit `- = # $` contradicting a `:` request on the same edge 
- AROM_KEKULE_INCONSISTENT (AROM: Error): no valid Kekule assignment is possible
- AROM_HUCKEL_FAIL (AROM: Error): Does not follow (4n + 2) rule

- STEREO_DOUBLE_CONFLICT (STEREO, Error): conflicting cis/trans specifications
- STEREO_DOUBLE_INSUFFICIENT (STEREO, Error): insufficient markers to define double-bond stereo

- NUM_OVERFLOW (NUM, Error): numeric literal exceeds supported bounds
- NUM_CLASS_OUT_OF_RANGE (NUM, Error): atom class exceeds 4 digits (max 9999)
- NUM_HCOUNT_OUT_OF_RANGE (NUM, Error): hydrogen count invalid (e.g., >9)
- NUM_CHARGE_OUT_OF_RANGE (NUM, Error): absolute charge exceeds hard limit (|q| > 15)
- NUM_ISOTOPE_TOO_LARGE (NUM, Error): isotope mass number exceeds 999
- NUM_CHIRALITY_OUT_OF_RANGE (NUM, Error): chirality parameter out of accepted range

- INTERNAL_ERROR (INTERNAL, Error): unexpected internal parser state

#### Warnings (style)

- STYLE_PREFER_BARE_ORGANIC_ATOM (STYLE, Warning): prefer bare organic subset atom over equivalent bracketed form 
- STYLE_PREFER_IMPLICIT_H (STYLE, Warning): prefer C over `[CH4]`. Exceptions: charge `[H+]`, `[H-]`, H-H bonds `[H][H]`, bridging H `[BH2]1([H])[BH2][H]1`, isotopes `[2H]`, `[3H]`
- STYLE_PREFER_BRACKET_FIELD_ORDER (STYLE, Warning): prefer {chirality}{H}{charge}{class} ordering in bracket atoms
- STYLE_PREFER_SIMPLE_CHARGE_SIGN (STYLE, Warning): prefer `+`/`-` over `+1`/`-1`
- STYLE_AVOID_DOUBLE_CHARGE_SIGN (STYLE, Warning): avoid `++`/`--`, use `+2`/`-2` instead
- STYLE_PREFER_SIMPLE_HCOUNT (STYLE, Warning): prefer H over H1 in bracket H-count
- STYLE_AVOID_EXPLICIT_SINGLE_BOND (STYLE, Warning): avoid explicit `-` when default applies
- STYLE_AVOID_EXPLICIT_AROMATIC_BOND (STYLE, Warning): avoid explicit `:` when aromatic default applies
- STYLE_AVOID_UNNECESSARY_GROUP (STYLE, Warning): avoid top-level grouping parentheses that do not affect connectivity
- STYLE_AVOID_REDUNDANT_NESTED_PARENS (STYLE, Warning): avoid redundant nested grouping parentheses such as `((...))`
- STYLE_PREFER_BRANCHES_BEFORE_RING_BONDS (STYLE, Warning): prefer branches, then ring bonds after atom
- STYLE_PREFER_FIRST_RING_ONE (STYLE, Warning): prefer starting ring numbering at 1
- STYLE_PREFER_CONSECUTIVE_RING_NUMBERING (STYLE, Warning): prefer consecutive ring numbering in the parsing sequence
- STYLE_AVOID_REUSED_RING_INDICES (STYLE, Warning): avoid reusing the same ring digit in a connected component
- STYLE_PREFER_SINGLE_DIGIT_RING_INDEX (STYLE, Warning): prefer single‑digit ring numbers for 1..9 instead of %01..%09
- STYLE_PREFER_SINGLE_RING_CLOSURE (STYLE, Warning): prefer single bond for ring-closure digits; avoid `=` or `#` at ring closure
- STYLE_AVOID_ADJACENT_RING_CLOSURES (STYLE, Warning): avoid starting a ring system on an atom that will carry two ring-closure bonds
- STYLE_PREFER_BOND_SYMBOL_AT_RING_OPEN (STYLE, Warning): prefer placing an explicit bond symbol at the ring opening rather than the closing index
- STYLE_AVOID_RING_CLOSURE_ACROSS_DOT (STYLE, Warning): avoid making ring closures across dot bonds, replace with simple bonds
- STYLE_PREFER_AROMATIC_FORM (STYLE, Warning): fully aromatic cycle is present but user provided a valid Kekulé input

- AROM_AVOID_MIXED_AROMATICITY (AROM, Warning): avoid mixing aromatic and non-aromatic atoms and/or bonds in the same ring
- AROM_AVOID_INCONSISTENT_AROMATICITY (AROM, Warning): Non-aromatic bonds (single, double, triple, quadruple) between aromatic atoms or aromatic bonds between non-aromatic atoms
- AROM_HUCKEL_INCONSISTENT (AROM, Warning): aromatic tokens used but HMO shows no significant delocalization

- STEREO_AVOID_UNNECESSARY_STEREO_DESCRIPTOR (STEREO, Warning): avoid stereo descriptors (bonds or chiral centers) at non-stereogenic elements

- NUM_ISOTOPE_UNCATALOGUED (NUM, Warning): isotope is not catalogued (unstable or too short-lived)
