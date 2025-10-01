### Diagnostics registry

This file defines the diagnostics taxonomy and stable codes for OpenSMILES parsing and linting in `umol-models-graph`.

#### Taxonomy

- Severity: Error | Warning
- Category: LEX (lexing), BRKT (bracket) BRCH (branch structure), GRP (groups), RING (ring rules), TOPO (graph topology), NUM (numeric constraints), BRKT (bracket fields), STEREO (stereochemistry), STYLE (style/normalization), SYN (parser-only fallback), INTERNAL (fallback)
- Code: stable UPPER_SNAKE identifiers (no numeric prefixes); messages are short and deterministic; parameters are carried as named fields
- Span: byte range [start, end) over the original input

#### Errors

##### Lexical / syntactic (during parse)

- LEX_INVALID_TOKEN (LEX, Error): input slice cannot be tokenized per lexical rules
- LEX_INVALID_ELEMENT (LEX, Error): invalid element symbol
- LEX_RING_INDEX_INVALID (LEX, Error): percent ring form invalid (must be exactly two digits)
- LEX_TRAILING_BOND (LEX, Error): bond symbol found at end of input or before terminator
- LEX_INTERTOKEN_WHITESPACE (LEX, Error): inter-token whitespace encountered
- LEX_COMMENT (LEX, Error): C-style comment encountered
- LEX_UNTERMINATED_BLOCK_COMMENT (LEX, Error): block comment not terminated before end of input
- LEX_DOT_BEFORE_RING (LEX, Error): dot placed before a ring index
- LEX_LEADING_DOT (LEX, Error): dot at start of input or component
- LEX_TRAILING_DOT (LEX, Error): dot at end of input or component
- LEX_MULTIPLE_DOTS (LEX, Error): consecutive dots without an intervening component
- LEX_LEADING_BOND (LEX, Error): bond token appears with no preceding atom (start of input or after a top-level group)
- LEX_CONSECUTIVE_BONDS (LEX, Error): consecutive bond tokens without an intervening atom or ring index
- LEX_TRAILING_BOND (LEX, Error): bond token found at end of input or component
- LEX_LEADING_RING (LEX, Error): ring index appears with no preceding atom (start of input or after a top-level group)

- BRKT_UNBALANCED_OPEN (BRKT, Error): '[' without a matching closing ']' (reported by parser)
- BRKT_UNBALANCED_CLOSE (BRKT, Error): ']' without a matching '[' (reported by parser)
- BRKT_UNEXPECTED_CLOSE (BRKT, Error): ']' outside of a bracket atom
- BRKT_FIELD_OUTSIDE (BRKT, Error): bracket-only field (e.g., @, +, ++) outside of a bracket atom
- BRKT_DUP_FIELD (BRKT, Error): duplicate bracket field of the same kind
- BRKT_EMPTY_CLASS (BRKT, Error): class field without a numeric value
- BRKT_H_ON_H (BRKT, Error): hydrogen element carries an H-count field

- BRCH_UNEXPECTED_CLOSE (BRANCH, Error): ')' without a matching '('
- BRCH_UNCLOSED (BRANCH, Error): open '(' not closed before component/end
- BRCH_DANGLING_BOND (BRANCH, Error): bond before ')' or component end
- BRCH_EMPTY_BRANCH (BRANCH, Error): empty branch (e.g., '()')
- GRP_LEADING_DOT (BRANCH, Error): group begins with a dot; groups must start with an atom or '('
- GRP_LEADING_BOND (BRANCH, Error): group begins with a bond; groups must start with an atom or '('

- RING_UNCLOSED (RING, Error): ring index opened but not closed by end of component/molecule
- RING_BOND_DIR_CONFLICT (RING, Error): conflicting up/down bond directions on the same ring closure
- RING_BOND_ORDER_CONFLICT (RING, Error): explicit bond orders on the two endpoints of a ring closure differ

##### Semantic (post-parse)

- TOPO_SELF_LOOP (TOPOLOGY, Error): self-loop edge present
- TOPO_PARALLEL_EDGES (TOPOLOGY, Error): more than one edge between the same atom pair

// TODO: Syntactic error = STYLE_AVOID_AROMATICITY_OUTSIDE_RINGS 
- AROM_ATOM_NOT_IN_RING (AROM, Error): aromatic atom not in ring
// TODO: Requires a defined valence model
- AROM_ATOM_SP2_INFEASIBLE (AROM, Error): atom must support at least one sp2-like configuration under valence model
// TODO: Syntactic error
- AROM_ATOM_ELEMENT_INVALID (AROM, Error): atom ineligible for aromatic systems
// TODO: Syntactic error = STYLE_AVOID_AROMATICITY_OUTSIDE_RINGS
- AROM_BOND_NOT_IN_RING (AROM, Error): aromatic bond symbol not in ring
// TODO: Syntactic error = RING_BOND_ORDER_CONFLICT?
- AROM_BOND_CONFLICT_EXPLICIT (AROM, Error): Explicit `- = # $` contradicting a `:` request on the same edge 
- AROM_BOND_ATOM_INELIGIBLE (AROM, Error): aromatic bond between atoms that are ineligible for aromatic systems
// TODO: How expensive is this check?
- AROM_KEKULE_INCONSISTENT (AROM: Error): no valid Kekule assignment is possible
// TODO: Mode-dependent error/warning
- AROM_HUCKEL_FAIL (AROM: Error): Does not follow (4n + 2) rule
- AROM_HUCKEL_WARN (AROM: Warning): Does not follow (4n + 2) rule

- NUM_OVERFLOW (NUM, Error): numeric literal exceeds supported bounds
- NUM_CLASS_OUT_OF_RANGE (NUM, Error): atom class exceeds 4 digits (max 9999)
- NUM_HCOUNT_OUT_OF_RANGE (NUM, Error): hydrogen count invalid (e.g., >9)
- NUM_CHARGE_OUT_OF_RANGE (NUM, Error): absolute charge exceeds hard limit (|q| > 15)
- NUM_ISOTOPE_TOO_LARGE (NUM, Error): isotope mass number exceeds 999
- NUM_CHIRAL_OUT_OF_RANGE (NUM, Error): chirality parameter out of accepted range
 - BRKT_CHIRAL_OUT_OF_RANGE (BRKT, Error): bracket chirality parameter out of accepted range

- STEREO_DOUBLE_CONFLICT (STEREO, Error): conflicting cis/trans specifications
- STEREO_DOUBLE_INSUFFICIENT (STEREO, Error): insufficient markers to define double-bond stereo

- INTERNAL_PARSER_STATE (INTERNAL, Error): unexpected internal parser state

#### Warnings (style)

- STYLE_PREFER_BRKT_FIELD_ORDER (STYLE, Warning): prefer [chirality][H][charge][class] ordering in bracket atoms
- STYLE_PREFER_SIMPLE_CHARGE_SIGN (STYLE, Warning): prefer [+]/[-] over [+1]/[-1]
- STYLE_PREFER_IMPLICIT_H (STYLE, Warning): prefer C over [CH4]. Exceptions: charge [H+], [H-], H-H bonds [H][H], bridging H [BH2]1([H])[BH2][H]1, isotopes [2H], [3H]
- STYLE_PREFER_H_OVER_H1 (STYLE, Warning): prefer H over H1 in bracket H-count
- STYLE_AVOID_UNNECESSARY_TOPLEVEL_PARENS (STYLE, Warning): avoid top-level grouping parentheses that do not affect connectivity
- STYLE_AVOID_REDUNDANT_NESTED_PARENS (STYLE, Warning): avoid redundant nested grouping parentheses such as '((...))'
- STYLE_PREFER_BRANCHES_BEFORE_RING_BONDS (STYLE, Warning): prefer branches, then ring bonds after atom
// TODO: Equivalence depends on valence model. See 39-opensmiles-semantic-issues-2025-09-30.md
- STYLE_PREFER_BARE_ORGANIC_ATOM (STYLE, Warning): prefer bare organic subset atom over equivalent bracketed form 
- STYLE_PREFER_AROMATIC_FORM (STYLE, Warning): fully aromatic cycle is present but user provided a valid Kekulé input
- STYLE_PREFER_SINGLE_DIGIT_RING_INDEX (STYLE, Warning): prefer single‑digit ring numbers for 1..9 instead of %01..%09
- STYLE_AVOID_EXPLICIT_AROMATIC_BOND (STYLE, Warning): avoid explicit ':' when aromatic default applies
- STYLE_AVOID_REUSED_RING_INDICES (STYLE, Warning): avoid reusing the same ring digit in a connected component
- STYLE_AVOID_EXPLICIT_SINGLE_BOND (STYLE, Warning): avoid explicit '-' when default applies
- STYLE_AVOID_EXPLICIT_H (STYLE, Warning): avoid explicit hydrogen when implicit is preferred
- STYLE_PREFER_FIRST_RING_ONE (STYLE, Warning): prefer starting ring numbering at 1
- STYLE_PREFER_CONSECUTIVE_RING_NUMBERING (STYLE, Warning): prefer consecutive ring numbering in the parsing sequence
- STYLE_AVOID_UNNECESSARY_STEREO_DESCRIPTOR (STEREO, Warning): avoid stereo descriptors (bonds or chiral centers) at non-stereogenic elements
// TODO: Same as AROM_ATOM_NOT_IN_RING?
- STYLE_AVOID_AROMATICITY_OUTSIDE_RINGS (AROM, Warning): avoid using aromatic atom symbols (bare or bracketed) and aromatic bonds outside of rings
- STYLE_AVOID_MIXED_AROMATICITY (AROM, Warning): avoid mixing aromatic and non-aromatic atoms and/or bonds in the same ring
- STYLE_PREFER_BOND_SYMBOL_AT_RING_OPEN (STYLE, Warning): prefer placing an explicit bond symbol at the ring opening rather than the closing index
- STYLE_PREFER_SINGLE_RING_CLOSURE (STYLE, Warning): prefer single bond for ring-closure digits; avoid '=' or '#' at ring closure
- STYLE_AVOID_ADJACENT_RING_CLOSURES (STYLE, Warning): avoid starting a ring system on an atom that will carry two ring-closure bonds
- STYLE_PREFER_DIRECT_BOND_OVER_DOT_RING_CLOSURE (STYLE, Warning): avoid making ring closures across dot bonds, replace with simple bonds
// TODO: What about the biphenyl example?
- STYLE_AVOID_INCONSISTENT_AROMATICITY (AROM, Warning): Non-aromatic bonds (single, double, triple, quadruple) between aromatic atoms or aromatic bonds between non-aromatic atoms

- NUM_HCOUNT_EXCEEDS_MAX_IMPLICIT (NUM, Warning): H-count exceeds element's max implicit hydrogens
- NUM_CHARGE_OUTSIDE_ELEMENT_RANGE (NUM, Warning): charge outside element-supported bounds
- NUM_CHARGE_EXCEEDS_VALENCE_ELECTRONS (NUM, Warning): charge exceeds valence electron count
- NUM_ISOTOPE_UNCATALOGUED (NUM, Warning): isotope is not catalogued (unstable or too short-lived)

#### Sources (post-parse modality)

- Parser-mapped (reported directly from FSM parser `ParseError`):
  - BRKT_UNBALANCED_OPEN, BRKT_UNBALANCED_CLOSE
  - BRCH_UNCLOSED, BRCH_UNEXPECTED_CLOSE, BRCH_EMPTY_BRANCH
  - LEX_LEADING_BOND, LEX_CONSECUTIVE_BONDS, LEX_TRAILING_BOND
  - LEX_LEADING_DOT, LEX_TRAILING_DOT, LEX_MULTIPLE_DOTS
  - LEX_LEADING_RING, LEX_RING_INDEX_INVALID
  - LEX_INTERTOKEN_WHITESPACE, LEX_COMMENT, LEX_UNTERMINATED_BLOCK_COMMENT
  - RING_UNCLOSED, RING_SELF_LOOP, RING_TWO_MEMBER, RING_MULTIPLE_RINGS, RING_BOND_DIR_CONFLICT, RING_BOND_ORDER_CONFLICT

- Post-parse derived (simple context around `ParseError`):
  - LEX_DOT_BEFORE_RING (from LeadingRing with preceding '.')
  - BRCH_DANGLING_BOND (from TrailingBond before ')'/component end)
  - STYLE_FIRST_RING_NOT_ONE, STYLE_NONCONSECUTIVE_RING_NUMBERING (sequence over parsed ring indices when parse succeeds)
  - STYLE_UNNECESSARY_PERCENT_RING_INDEX (scan for `%0n` where n in 1..9)

- Bracket-utils-based (parsed via `parser::utils::parse_bracket`):
  - STYLE_BRKT_ORDER, STYLE_BRKT_ORGANIC, STYLE_CHARGE_SIGN_SIMPLE, STYLE_HCOUNT_ONE_SIMPLE
  - NUM_ISOTOPE_TOO_LARGE, NUM_CHIRAL_OUT_OF_RANGE
