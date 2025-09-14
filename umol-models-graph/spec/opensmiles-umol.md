## OpenSMILES (UMOL) — Grammar and Normative Rules
Version 0.1 (2025-09-13)

This document specifies the UMOL dialect of the SMILES language. The EBNF grammar and lexical rules in this document are normative. The implementation grammar (e.g., LALRPOP) is non-normative; where any discrepancy exists, this document prevails.

### Versioning and Compatibility

This specification is versioned as part of the UMOL project. Minor revisions may clarify wording and add non‑breaking rules that do not invalidate previously valid inputs nor validate previously invalid inputs. Major revisions may introduce breaking changes to tokens, grammar, or semantics. Implementations claiming conformance must indicate the spec version they target. When behavior differs across versions, the earlier version’s behavior remains valid only for that version; new implementations should target the latest version.

### Conformance

This specification uses the terminology of RFC 2119 and RFC 8174. A conforming parser must accept inputs matching the grammar and token rules herein and must reject inputs that violate the grammar or any semantic constraint defined in this document. Specific error messages are not required.

### Lexical Conventions

Input is ASCII and case-sensitive. Tokens are recognized using a maximal‑munch policy: at each source position the longest possible token is chosen. If multiple tokens of equal maximal length can begin at a position, the one appearing earlier in the token definition order below is selected. Whitespace (SPACE, TAB, LINEFEED, CARRIAGE_RETURN) may appear only as trailing terminator at end of input; inter-token whitespace is not permitted. A valid SMILES is followed only by optional trailing whitespace and end-of-string.

### Tokens (normative summary)

| Name | Definition | Notes |
|---|---|---|
| DIGIT | 0–9 | |
| NUMBER | DIGIT+ | Base‑10, unsigned |
| PERCENT_RING | '%' DIGIT DIGIT | Exactly two digits; 00–99; leading zeros allowed |
| DOT | '.' | Component separator |
| BOND | '-' '=' '#' '$' ':' '/' '\\' | Bond symbols |
| BRACKET_OPEN | '[' | |
| BRACKET_CLOSE | ']' | |
| CHIRALITY | '@' '@@' '@TH' DIGIT '@AL' DIGIT '@SP' DIGIT '@TB' NUMBER '@OH' NUMBER | Ranges defined below |
| CHARGE | '+' '-' '++' '--' '+' DIGIT [DIGIT]? '-' DIGIT [DIGIT]? | +/− optionally with one or two digits; '++'/'--' equal ±2 |
| WHITESPACE | SPACE TAB LINEFEED CARRIAGE_RETURN | Trailing only |
| END_OF_STRING | (end sentinel) | Not a character |

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

element_symbols ::= 'H' | 'He' | 'Li' | 'Be' | 'B' | 'C' | 'N' | 'O' | 'F' | 'Ne' | 'Na' | 'Mg' | 'Al' | 'Si' | 'P' | 'S' | 'Cl' | 'Ar' | 'K' | 'Ca' | 'Sc' | 'Ti' | 'V' | 'Cr' | 'Mn' | 'Fe' | 'Co' | 'Ni' | 'Cu' | 'Zn' | 'Ga' | 'Ge' | 'As' | 'Se' | 'Br' | 'Kr' | 'Rb' | 'Sr' | 'Y' | 'Zr' | 'Nb' | 'Mo' | 'Tc' | 'Ru' | 'Rh' | 'Pd' | 'Ag' | 'Cd' | 'In' | 'Sn' | 'Sb' | 'Te' | 'I' | 'Xe' | 'Cs' | 'Ba' | 'Hf' | 'Ta' | 'W' | 'Re' | 'Os' | 'Ir' | 'Pt' | 'Au' | 'Hg' | 'Tl' | 'Pb' | 'Bi' | 'Po' | 'At' | 'Rn' | 'Fr' | 'Ra' | 'Rf' | 'Db' | 'Sg' | 'Bh' | 'Hs' | 'Mt' | 'Ds' | 'Rg' | 'Cn' | 'Fl' | 'Lv' | 'La' | 'Ce' | 'Pr' | 'Nd' | 'Pm' | 'Sm' | 'Eu' | 'Gd' | 'Tb' | 'Dy' | 'Ho' | 'Er' | 'Tm' | 'Yb' | 'Lu' | 'Ac' | 'Th' | 'Pa' | 'U' | 'Np' | 'Pu' | 'Am' | 'Cm' | 'Bk' | 'Cf' | 'Es' | 'Fm' | 'Md' | 'No' | 'Lr'

aromatic_symbols ::= 'b' | 'c' | 'n' | 'o' | 'p' | 's' | 'se' | 'as'

chiral ::= '@' | '@@' | '@TH' DIGIT | '@AL' DIGIT | '@SP' DIGIT | '@TB' NUMBER | '@OH' NUMBER

hcount ::= 'H' | 'H' DIGIT

charge ::= '+' | '-' | '++' | '--' | '+' two_digits | '-' two_digits

class ::= ':' NUMBER

bond ::= '-' | '=' | '#' | '$' | ':' | '/' | '\\'

ringbond ::= bond? DIGIT | bond? PERCENT_RING

node ::= atom ringbond* branch*
branch ::= '(' connector? chain ')'
connector ::= bond | dot
chain ::= node ( connector? node )*

dot ::= '.'

smiles ::= ws END_OF_STRING | chain ws END_OF_STRING
ws ::= ( ' ' | '\t' | '\n' | '\r' )*
```

### Semantic Constraints

A ring index is either a single digit in the range 0–9 or a percent form in the range 00–99. Percent forms are exactly two digits; leading zeros are permitted. Ring indices compare by numeric value (leading zeros ignored), e.g., '1' matches '%01' and '0' matches '%00'. '%0' (percent with a single digit) is invalid. Open rings may be closed across component separators; at the end of input, any unclosed ring causes the input to be rejected. A ring closure must not create a self-loop or a two‑member cycle.

When both sides of a ring closure specify a bond direction and they differ, the input is rejected. When directions do not conflict, the closing-side specification takes precedence; otherwise the opening-side specification applies. Consecutive bond tokens without an intervening atom or ring index are invalid, and a trailing bond token at end‑of‑input is invalid.

### Aromaticity

Lowercase aromatic atom tokens designate aromatic atoms. The '*' token may appear adjacent to or within aromatic systems; its presence does not by itself imply aromaticity. An explicit ':' always produces an aromatic bond without altering atom aromaticity. When a bond is implicit and both adjacent atoms are aromatic, the bond order defaults to aromatic; otherwise the implicit default is single. This defaulting applies during chain growth and ring closure.

### Implicit Hydrogens

Inside a bracket atom, 'H' without a digit sets the hydrogen count to one; 'H' followed by a single digit sets the count to that value. A bracket atom whose element is hydrogen ('H') MUST NOT include an 'H' count field; such forms are invalid (e.g., '[HH]', '[HH1]'). Outside brackets, the parser does not infer hydrogen counts from valence; an unspecified hydrogen count is treated as zero. Stereochemical validations may consult the bracket hydrogen count when checking neighbor totals.

### Double‑Bond Stereochemistry

Directional markers '/' and '\\' attached to single bonds adjacent to a double bond encode relative geometry. The parser collects such markers and classifies double‑bond geometry when it is unambiguous; otherwise the geometry remains unknown ("either"). In ring‑directed closures where only one side supplies a determinative marker, the geometry is recorded as unknown/either.

For cumulenes, UMOL does not propagate endpoint up/down markers across chains; only the local alkene E/Z rule is applied. Endpoint stereomarkers (up/down) on cumulenes with more than two consecutive double bonds are not interpreted and do not establish geometry. Even‑length systems such as allenes may be characterized by allenic stereochemistry via '@ALn' on the central atom; '@' and '@@' are aliased to '@AL1' and '@AL2' when the center has exactly two incident double bonds. Allenic stereomarkers are not available for even‑count cumulenes with more than two consecutive double bonds (n > 2). Validation of '@ALn' assignments is performed in the semantic pass.

### Components and Molecule Finalization

A dot separates components. The first atom of a component does not connect to the previous component unless a subsequent ring closure connects them, in which case the components merge into a single molecule. If no rings are open when a component ends, the molecule may be finalized immediately; if rings are open at a component boundary, finalization is deferred until those rings are closed. At end of input, open rings cause rejection.

### Branch‑Local Components

Parentheses introduce a branch attached at the current node. Inside a branch, dots separate branch‑local components. These components belong to the enclosing molecule. Ring indices inside a branch may connect to atoms outside the branch, and such connections are allowed. Upon leaving a branch, branch‑local components that did not merge by bonding or ring closure remain as separate components within the same molecule unless later merged by ring closure.

### Numeric and Bracket Fields

In a bracket atom, the isotope field is a non‑negative decimal integer; zero is permitted and does not alter the element identity. The parser does not validate chemical plausibility of isotopes. The chiral, hydrogen‑count, charge, and class fields may appear in any order, at most once each. The charge field is one of '+', '-', '++', '--', or a '+' or '-' followed by one or two decimal digits; '+0' and '-0' are accepted as zero. In a bracket atom, 'H' without a digit sets hydrogen_count to one; 'H' followed by a single digit in 0–9 sets it to that value; 'H0' is accepted and yields zero. The class field consists of ':' followed by a non‑negative integer; leading zeros are permitted; zero is allowed. If the class field is omitted, the class value defaults to 0.

### Unknown Atom

The '*' token denotes an unknown atom without element assignment and does not imply aromaticity. It is permitted anywhere a normal atom may appear, including in aromatic rings. No implicit hydrogens are inferred for '*' outside brackets. Inside a bracket atom, '*' may carry fields (isotope, chirality, hydrogen count, charge, class) subject to the same syntax as elements; unspecified class defaults to 0.

### Numeric Limits and Overflow

Numeric literals are base‑10 and unsigned. Implementations must reject any numeric literal that overflows their supported domain. In UMOL, numeric fields are interpreted as unsigned 32‑bit values; values outside the range 0 to 4,294,967,295 are rejected. Ring index ranges are further constrained to 0–9 for single‑digit and 00–99 for percent forms (exactly two digits). The hydrogen count after 'H' is a single digit in 0–9.

### Error Policy and Diagnostics

Any violation of the grammar or a semantic constraint results in rejection. Implementations SHOULD provide structured diagnostics with a stable code, severity, category, span, and message. The canonical diagnostics registry is maintained in `spec/opensmiles-errors.md`.

### Conformance Criteria

A conforming implementation accepts all valid examples in the accompanying test suite and rejects all invalid examples. It enforces the ring, directed‑ring, bond, bracket, chirality, aromaticity, hydrogen, and numeric constraints described above.

### Test Suite Reference

The repository includes a machine‑readable test suite organized into valid, edge, and invalid categories under `umol-models-graph/src/io/smiles/parser/tests/`. Implementations should use these examples to verify conformance.

### Differences from OpenSMILES (UMOL) — informational

UMOL accepts '0' as a single‑digit ring index. UMOL requires percent ring indices to be exactly two digits and accepts leading‑zero forms as numeric equivalents (e.g., '1' == '%01', '0' == '%00'). When both endpoints are aromatic and a bond is implicit, UMOL defaults that bond to aromatic. UMOL allows dots inside branches to start branch‑local components which may later merge via rings. UMOL recognizes '++' and '--' as charge tokens equal to ±2.

### Invalid Forms — informational

The following examples illustrate inputs that must be rejected: "%0" (invalid percent ring index), "C1CC" (unclosed ring index), "C11" (self‑loop ring closure), "C12C21" (two‑member ring), "C/1CC\\1" (conflicting directions on a ring bond), "C-" (trailing bond), and "C==C" (consecutive bond tokens without an intervening atom or ring index).

### Non‑Normative Appendix: Implementation Grammar

The implementation grammar (e.g., LALRPOP) may be included as an appendix for reference. When discrepancies arise, the EBNF grammar above is authoritative.
