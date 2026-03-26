## OpenSMILES — UMOL Formalization

Version 1.1 (2026-02-12)

### Note

Review the BALSA specification [ChemRxiv](https://doi.org/10.26434/chemrxiv-2022-01ltp)

### Preface

This document is not intended as an extension or alternative to the official OpenSMILES specification (http://opensmiles.org/opensmiles.html). Rather, it formalizes aspects of the official spec where the original is ambiguous or lacks formal structure, providing precise lexical rules, grammar, and semantic constraints suitable for implementation.

The EBNF grammar and lexical rules in this document are normative for the UMOL implementation.

### Versioning

This specification is considered stable. Future revisions are limited to bugfixes and clarifications; no new features will be added. Any extensions (e.g., CXSMILES, SMARTS) will be documented as separate addenda rather than modifications to this base specification.

### Conformance

This specification uses the terminology of RFC 2119 and RFC 8174. The key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD", "SHOULD NOT", "RECOMMENDED", "MAY", and "OPTIONAL" in this document are to be interpreted as described in those RFCs.

A conforming parser MUST accept inputs matching the grammar and token rules herein and MUST reject inputs that violate the grammar. Semantic constraints SHALL be validated in a separate pass after parsing. Specific error messages are not required.

### Document Structure

This specification is organized into three parts:

1. **Lexical** — character classes, tokens, and lexical conventions
2. **Syntax** — grammar rules and syntactic constraints
3. **Semantics** — validation rules applied after parsing

---

## Part 1: Lexical

This section defines the character classes, tokens, and lexical conventions for SMILES strings.

### Lexical Conventions

Input MUST be ASCII and is case-sensitive. Tokens SHALL be recognized using a maximal-munch policy: at each source position the longest possible token is chosen. If multiple tokens of equal maximal length can begin at a position, the one appearing earlier in the token definition order below SHALL be selected. Whitespace (SPACE, TAB, LINEFEED, CARRIAGE_RETURN) MAY appear as the only content (empty molecule) or as trailing terminator at end of input; inter-token whitespace MUST NOT appear.

### Tokens (normative summary)

| Name          | Definition                                                             | Notes                                                       |
| ------------- | ---------------------------------------------------------------------- | ----------------------------------------------------------- |
| DIGIT         | 0-9                                                                    |                                                             |
| NUMBER        | DIGIT+                                                                 | Base-10, unsigned                                           |
| PERCENT_RING  | '%' DIGIT DIGIT                                                        | Exactly two digits; 00-99; leading zeros allowed            |
| DOT           | '.'                                                                    | Component separator                                         |
| BOND          | '-' '=' '#' '$' ':' '/' '\'                                            | Bond symbols                                                |
| BRACKET_OPEN  | '['                                                                    |                                                             |
| BRACKET_CLOSE | ']'                                                                    |                                                             |
| CHIRALITY     | '@' '@@' '@TH' DIGIT '@AL' DIGIT '@SP' DIGIT '@TB' NUMBER '@OH' NUMBER | Ranges defined below                                        |
| CHARGE        | '+' '-' '++' '--' '+' DIGIT [DIGIT]? '-' DIGIT [DIGIT]?                | +/- optionally with one or two digits; '++'/'--' equal +/-2 |
| WHITESPACE    | SPACE TAB LINEFEED CARRIAGE_RETURN                                     | Trailing only                                               |
| END_OF_STRING | (end sentinel)                                                         | Not a character                                             |

---

## Part 2: Syntax

This section defines the grammar rules and syntactic constraints for well-formed SMILES strings. The parser SHALL accept any input conforming to this syntax; semantic validity is checked separately.

### Grammar (normative EBNF)

```ebnf
DIGIT ::= '0' | '1' | '2' | '3' | '4' | '5' | '6' | '7' | '8' | '9'
NUMBER ::= DIGIT+
PERCENT_RING ::= '%' DIGIT DIGIT

two_digits ::= DIGIT | DIGIT DIGIT

atom ::= bracket_atom | aliphatic_organic | aromatic_organic | '*'

aliphatic_organic ::= 'B' | 'C' | 'N' | 'O' | 'S' | 'P' | 'F' | 'Cl' | 'Br' | 'I'
aromatic_organic ::= 'b' | 'c' | 'n' | 'o' | 's' | 'p'

bracket_atom ::= '[' isotope? symbol bracket_field* ']'
bracket_field ::= chiral | hcount | charge | class
symbol ::= element_symbols | aromatic_symbols | '*'
isotope ::= NUMBER

element_symbols ::= 'H' | 'He' | 'Li' | 'Be' | 'B' | 'C' | 'N' | 'O' | 'F' | 'Ne' | 'Na' | 'Mg' | 'Al' | 'Si' | 'P' | 'S' | 'Cl' | 'Ar' | 'K' | 'Ca' | 'Sc' | 'Ti' | 'V' | 'Cr' | 'Mn' | 'Fe' | 'Co' | 'Ni' | 'Cu' | 'Zn' | 'Ga' | 'Ge' | 'As' | 'Se' | 'Br' | 'Kr' | 'Rb' | 'Sr' | 'Y' | 'Zr' | 'Nb' | 'Mo' | 'Tc' | 'Ru' | 'Rh' | 'Pd' | 'Ag' | 'Cd' | 'In' | 'Sn' | 'Sb' | 'Te' | 'I' | 'Xe' | 'Cs' | 'Ba' | 'La' | 'Ce' | 'Pr' | 'Nd' | 'Pm' | 'Sm' | 'Eu' | 'Gd' | 'Tb' | 'Dy' | 'Ho' | 'Er' | 'Tm' | 'Yb' | 'Lu' | 'Hf' | 'Ta' | 'W' | 'Re' | 'Os' | 'Ir' | 'Pt' | 'Au' | 'Hg' | 'Tl' | 'Pb' | 'Bi' | 'Po' | 'At' | 'Rn' | 'Fr' | 'Ra' |  'Ac' | 'Th' | 'Pa' | 'U' | 'Np' | 'Pu' | 'Am' | 'Cm' | 'Bk' | 'Cf' | 'Es' | 'Fm' | 'Md' | 'No' | 'Lr' | 'Rf' | 'Db' | 'Sg' | 'Bh' | 'Hs' | 'Mt' | 'Ds' | 'Rg' | 'Cn' | 'Nh' | 'Fl' | 'Mc' | 'Lv' | 'Ts' | 'Og'

aromatic_symbols ::= 'b' | 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'

chiral ::= '@' | '@@' | '@TH' DIGIT | '@AL' DIGIT | '@SP' DIGIT | '@TB' NUMBER | '@OH' NUMBER

hcount ::= 'H' | 'H' DIGIT

charge ::= '+' | '-' | '++' | '--' | '+' two_digits | '-' two_digits

class ::= ':' NUMBER

bond ::= '-' | '=' | '#' | '$' | ':' | '/' | '\\'

ringbond ::= bond? DIGIT | bond? PERCENT_RING

node ::= atom ( ringbond | branch )*
branch ::= '(' connector? chain ')'   // branch may start with bond/dot (connector)
connector ::= bond | dot
group ::= '(' chain ')'               // group cannot start with connector
chain ::= (node | group) ( connector? (node | group) )*

dot ::= '.'

smiles ::= ws END_OF_STRING | chain ws END_OF_STRING
ws ::= ( ' ' | '\t' | '\n' | '\r' )*
```

### Implicit Bond Default

When two atoms are adjacent in the chain without an explicit bond token (i.e., the `connector` in the grammar is absent), the bond order SHALL default as follows:

- If both adjacent atoms are aromatic, the bond order SHALL default to aromatic.
- Otherwise, the bond order SHALL default to single.

This rule applies during chain growth and ring closure.

### Syntactic Constraints

Parentheses have two roles:

- Branch: '(' following an atom opens a branch attached to that atom; ')' closes the branch and restores the attach point. A branch MUST contain at least one atom; empty branches '()' are invalid. An atom MAY be followed by more than one branch definition.
- Group: '(' at top level (i.e., when no branch attach point is pending) is grouping only and does not create or remove bonds. Grouping parentheses MAY nest. Redundant grouping such as '(CC)' or '((CC))' is valid and connectivity-preserving; implementations MAY warn as a style issue (STYLE_AVOID_UNNECESSARY_GROUP, STYLE_AVOID_REDUNDANT_NESTED_PARENS). An empty top-level group '()' MUST be rejected (SYN_EMPTY_GROUP). A group MUST span the entire component or chain at its level; atoms or groups following a closed group at the same level MUST be rejected (SYN_NONFINAL_GROUP). For example, '(CC)' is valid, but '(CC)C' and '(CC)(CC)' are both invalid.

All parentheses MUST be paired. A ')' without a matching '(' MUST be rejected (SYN_UNBALANCED_CLOSE_PAREN). At end of input, any unclosed '(' MUST be rejected (SYN_UNBALANCED_OPEN_PAREN). Grouping does not alter dot/component semantics; for example, '(CC.CC)' is equivalent to 'CC.CC'. The empty molecule SHALL be represented only by the empty input (or whitespace-only input).

In branches, an initial connector (bond or dot) after '(' is permitted and applies to the first edge in the branch. In groups, the first token after '(' MUST be an atom or another group; a leading connector (bond or dot) inside a group MUST be rejected. A connector (bond) at the very start of input or immediately following a top-level group MUST be rejected as there is no attach point (SYN_LEADING_BOND). A ring index at the very start of input or immediately following a top-level group MUST be rejected as there is no attach point (SYN_LEADING_RING).

Ring closures and branches MAY appear in any order after an atom and MAY interleave arbitrarily; the grammar `node ::= atom ( ringbond | branch )*` permits both orders and mixtures.

### Ring Closures

A ring index is either a single digit in the range 0–9 or a percent form in the range 00–99. The maximum ring index is 99. Percent forms MUST consist of exactly two digits; leading zeros are permitted. Ring indices SHALL compare by numeric value (leading zeros ignored), e.g., '1' matches '%01' and '0' matches '%00'. '%0' (percent with a single digit) MUST be rejected. Open rings MAY be closed across component separators; at the end of input, any unclosed ring MUST cause rejection. Self-loops (e.g., 'C11') and two-member cycles (e.g., 'C1C1') are syntactically valid but SHALL be rejected during semantic validation (see Part 3).

Ring indices MAY be reused after closure. Each occurrence of a ring index opens or closes a ring bond: the first occurrence opens the ring, the second closes it, the third reopens it, and so on. Consequently, every ring index MUST occur an even number of times in the input; an odd count for any ring index MUST be rejected as an unclosed ring. When the same ring index appears on two adjacent atoms (e.g., 'C1C1'), each pair forms a separate ring bond. Implementations MAY emit a style warning when a ring index is reused.

When both sides of a ring closure specify a bond direction and they differ, the input MUST be rejected (SYN_MISMATCHED_RING_BOND_DIRS). When directions do not conflict, the closing-side specification takes precedence; otherwise the opening-side specification applies. Consecutive bond tokens without an intervening atom or ring index MUST be rejected (SYN_CONSECUTIVE_BONDS), and a trailing bond token at end-of-input MUST be rejected (SYN_TRAILING_BOND).

Ring bond order semantics:

- If both endpoints of a ring closure specify a bond order and they differ, the input MUST be rejected (SYN_MISMATCHED_RING_BOND_ORDERS).
- If exactly one endpoint specifies a non-single order, the input is valid; the specified order applies to the ring bond. If both endpoints specify orders, they MUST be the same or the input MUST be rejected (SYN_MISMATCHED_RING_BOND_ORDERS).
- A ring bond with a non-single order MUST NOT carry a direction marker; such combinations MUST be rejected (SYN_MISMATCHED_RING_BOND_ORDERS).
- Error position convention: implementations SHOULD report errors at the closing ring index token (for single-digit rings, the digit; for percent rings, the '%').

### Aromaticity

Lowercase aromatic atom tokens designate aromatic atoms. The '*' token MAY appear adjacent to or within aromatic systems; its presence does not by itself imply aromaticity. An explicit ':' SHALL always produce an aromatic bond without altering atom aromaticity. Implicit bond defaulting between aromatic atoms is defined in the Implicit Bond Default section above.

### Double-Bond Stereochemistry

Directional markers '/' and '\' attached to single bonds adjacent to a double bond encode relative geometry. The parser SHALL collect such markers and classify double-bond geometry when it is unambiguous; otherwise the geometry SHALL remain unknown ("either"). In ring-directed closures where only one side supplies a determinative marker, the geometry SHALL be recorded as unknown/either.

Cumulenes (consecutive double bonds) follow the same stereochemistry rules extended across the chain:

- **Odd number of double bonds** (1, 3, 5, ...): Extended cis/trans (E/Z) geometry applies. Directional markers '/' and '\' on bonds to the terminal atoms establish the relative configuration across the entire cumulene axis. Examples: `F/C=C=C=C/F` (trans-difluorobutatriene), `F/C=C=C=C\F` (cis-difluorobutatriene).
- **Even number of double bonds** (2, 4, 6, ...): Extended tetrahedral (axial) chirality applies. The '@AL1' and '@AL2' markers on the central atom of the cumulene system specify the configuration; the "neighbor" atoms to which chirality refers are at the ends of the allene system. For simple allenes (2 double bonds), '@' and '@@' are aliased to '@AL1' and '@AL2' when the center has exactly two incident double bonds. Example: `NC(Br)=[C@]=C(O)C`.

Validation of cumulene stereochemistry SHALL be performed in the semantic pass.

### Top-Level Structure

A SMILES string represents zero or more molecular components. The dot separator '.' delimits components within the string. Each component is a connected molecular graph; the complete SMILES string MAY represent a disconnected system of multiple molecules. An empty SMILES string (or whitespace-only input) SHALL represent the empty molecule (zero atoms, zero components).

### Dot Constraints

A dot '.' separates components within a SMILES string. The following constraints apply:

- A dot at the start of input or at the start of a component MUST be rejected (SYN_LEADING_DOT).
- A dot at the end of input or at the end of a component MUST be rejected (SYN_TRAILING_DOT).
- Consecutive dots without an intervening atom MUST be rejected (SYN_CONSECUTIVE_DOTS).
- A dot immediately before a ring index MUST be rejected (SYN_DOT_BEFORE_RING).

Every component delimited by dots MUST contain at least one atom.

### Components and Molecule Finalization

A dot separates components. The first atom of a component SHALL NOT connect to the previous component unless a subsequent ring closure connects them, in which case the components merge into a single connected graph. If no rings are open when a component ends, the molecule MAY be finalized immediately; if rings are open at a component boundary, finalization SHALL be deferred until those rings are closed. At end of input, open rings MUST cause rejection.

### Branch-Local Components

Parentheses introduce a branch attached at the current node. Inside a branch, dots separate branch-local components. These components belong to the enclosing molecule. Ring indices inside a branch MAY connect to atoms outside the branch, and such connections are permitted. Upon leaving a branch, branch-local components that did not merge by bonding or ring closure SHALL remain as separate components within the same molecule unless later merged by ring closure.

### Numeric and Bracket Fields

In a bracket atom, the isotope field is a non-negative decimal integer; zero is permitted and SHALL NOT alter the element identity. The parser SHALL NOT validate chemical plausibility of isotopes. If the isotope field is absent, the isotope composition SHALL be considered undetermined (no specific isotope is implied).

**Deviation from OpenSMILES:** The chiral, hydrogen-count, charge, and class fields MAY appear in any order, at most once each. The OpenSMILES specification mandates a fixed field order; this specification relaxes that constraint.

The charge field is one of '+', '-', '++', '--', or a '+' or '-' followed by one or two decimal digits; '+0' and '-0' SHALL be accepted as zero charge. If the charge field is absent, the charge value SHALL default to 0.

In a bracket atom, 'H' without a digit SHALL set the hydrogen count to 1; 'H' followed by a single digit in 0–9 SHALL set it to that value; 'H0' is accepted and yields zero. If the 'H' field is omitted inside brackets, the hydrogen count SHALL default to 0. Bracket atoms SHALL NOT undergo implicit hydrogen calculation; the hydrogen count is exactly as specified (or 0 if omitted).

A bracket atom whose element is hydrogen ('H') MUST NOT include an 'H' count field; such forms (e.g., '[HH]', '[HH1]') MUST be rejected.

The class field consists of ':' followed by a non-negative integer; leading zeros are permitted; zero is allowed. **Deviation from OpenSMILES:** The absence of the class field SHALL NOT assign the atom to class 0. A bracket atom without a ':' field has no class assignment. This matches the treatment of organic subset atoms, which also have no class assignment.

### Organic Subset

Only atoms B, C, N, O, P, S, F, Cl, Br, and I MAY appear without brackets. For these atoms:

- Isotope composition SHALL be considered undetermined (no specific isotope is implied).
- Charge SHALL default to 0.
- No chiral specification SHALL be assumed.
- No class SHALL be assigned.
- The implicit hydrogen count SHALL be determined by summing the bond orders of the bonds connected to the atom. If that sum is equal to a known valence for the element or is greater than any known valence, then the implicit hydrogen count SHALL be 0. Otherwise the implicit hydrogen count SHALL be the difference between that sum and the next highest known valence.

This implicit hydrogen calculation applies exclusively to organic subset atoms. Bracket atoms SHALL NOT undergo implicit hydrogen calculation.

| Element   | Normal valences |
| --------- | --------------- |
| B         | 3               |
| C         | 4               |
| N         | 3, 5            |
| O         | 2               |
| P         | 3, 5            |
| S         | 2, 4, 6         |
| F,Cl,Br,I | 1               |

### Unknown Atom

The '*' token denotes an unknown atom without element assignment and does not imply aromaticity. It is permitted anywhere a normal atom may appear, including in aromatic rings. No implicit hydrogens SHALL be inferred for '*' outside brackets. Inside a bracket atom, '*' MAY carry fields (isotope, chirality, hydrogen count, charge, class) subject to the same syntax as elements.

**Note:** The '*' production is included in the grammar for completeness. The basic SMILES parser (`parse_smiles`) SHALL reject wildcard atoms; they are accepted only by the extended parser (`parse_extended_smiles`). This separation exists because wildcard atoms do not denote a concrete element and cannot be resolved to a valid molecular graph; the extended parser produces a distinct representation that accommodates this.

### Numeric Limits and Overflow

Numeric literals are base-10 and unsigned. Implementations MUST reject any numeric literal that overflows their supported domain (NUM_OVERFLOW). In UMOL, the following domain-specific limits apply:

- **Ring indices**: single-digit indices SHALL be in the range 0–9; percent-form indices SHALL be in the range 00–99 (exactly two digits), yielding a maximum ring index of 99.
- **Hydrogen count**: the digit after 'H' SHALL be in the range 0–9 (NUM_HCOUNT_OUT_OF_RANGE).
- **Charge**: the absolute charge magnitude MUST NOT exceed 15; values with |q| > 15 MUST be rejected (NUM_CHARGE_OUT_OF_RANGE).
- **Isotope**: the isotope mass number MUST NOT exceed 999 (NUM_ISOTOPE_TOO_LARGE).
- **Atom class**: the class number MUST NOT exceed 9999 (NUM_CLASS_OUT_OF_RANGE).
- **Chirality**: chirality parameters for @TB (1–20) and @OH (1–30) MUST be within their respective ranges (NUM_CHIRALITY_OUT_OF_RANGE).

---

## Part 3: Semantics

This section defines validation rules applied after parsing to ensure chemical and topological validity. These constraints SHALL be checked on the parsed molecule structure, not during lexical analysis or grammar parsing.

### Topology Constraints

The following topological structures are syntactically valid but semantically invalid:

- **Self-loops**: A ring closure connecting an atom to itself (e.g., 'C11') MUST be rejected.
- **Two-member cycles**: A ring closure creating a cycle of length two (e.g., 'C1C1') MUST be rejected. Note that multiple bonds between atoms via separate ring closures (e.g., 'C12C12') are valid and create parallel edges.

### Valence and Hydrogen Validation

_To be specified. Validation of implicit hydrogen counts, valence limits, and charge compatibility._

### Stereochemistry Validation

_To be specified. Validation of chiral center neighbor counts, allene geometry, and cis/trans bond consistency._

### Aromaticity Validation

_To be specified. Validation of aromatic ring membership and Kekulization._

---

## Appendix

### Error Policy and Diagnostics

Any violation of the grammar or a semantic constraint SHALL result in rejection. Implementations SHOULD provide structured diagnostics with a stable code, severity, category, span, and message. The canonical diagnostics registry is maintained in `spec/opensmiles-errors.md`.

### Conformance Criteria

A conforming implementation MUST accept all valid examples in the accompanying test suite and MUST reject all invalid examples. It SHALL enforce the ring, directed-ring, bond, bracket, chirality, aromaticity, hydrogen, and numeric constraints described above.

### Test Suite Reference

The repository includes a machine-readable test suite organized into valid, edge, and invalid categories under `umol-models-graph/src/io/smiles/parser/tests/`. Implementations SHOULD use these examples to verify conformance.

### Implementation Notes

This implementation accepts '0' as a single-digit ring index. Percent ring indices MUST be exactly two digits; leading-zero forms are accepted as numeric equivalents (e.g., '1' == '%01', '0' == '%00'). When both endpoints are aromatic and a bond is implicit, the bond defaults to aromatic. Dots inside branches start branch-local components which may later merge via rings. The tokens '++' and '--' are recognized as charge values equal to +/-2.

### Invalid Forms

**Syntactic errors** (rejected during parsing):

- `%0` — SYN_INVALID_RING_INDEX (requires exactly two digits)
- `C1CC` — SYN_UNBALANCED_RING_INDEX
- `C/1CC\1` — SYN_MISMATCHED_RING_BOND_DIRS
- `C-` — SYN_TRAILING_BOND
- `C==C` — SYN_CONSECUTIVE_BONDS
- `.C` — SYN_LEADING_DOT
- `C.` — SYN_TRAILING_DOT
- `C..C` — SYN_CONSECUTIVE_DOTS
- `C.1CC1` — SYN_DOT_BEFORE_RING
- `1CC1` — SYN_LEADING_RING
- `()` — SYN_EMPTY_GROUP
- `C()` — SYN_EMPTY_BRANCH
- `(CC)C` — SYN_NONFINAL_GROUP
- `(CC)(CC)` — SYN_NONFINAL_GROUP

**Semantic errors** (rejected during validation):

- `C11` — TOPO_SELF_LOOP_RING
- `C1C1` — TOPO_SELF_LOOP_RING (two-member cycle; but `C12C12` is valid: parallel edges)
