# SMILES Parser Configuration Design

## Background

The current parser has premature configuration features (UMOL dialect, comments, extended
whitespace) that were added speculatively. This document outlines a principled approach to
parser configuration based on actual dialect requirements.

## Approach

1. Remove premature parse flags and corresponding code
2. Enumerate features from existing SMILES implementations
3. Reverse-engineer Daylight and CXSMILES dialects
4. Group features logically with feature gates
5. Consider SMARTS/SMIRKS and other extensions
6. Decide on presets vs dialects
7. Create implementation plan

## Priorities

a. Feature parity with open source parsers (RDKit, OpenBabel, CDK, Indigo)
b. SMARTS parsing support
c. Daylight SMILES / CXSMILES extensions (from docs and examples)
d. Additional formats, inputs, extensions

## Design Decisions

### Separation of Concerns

Parsing (syntax) is strictly separated from semantic processing:
- **Parser**: Faithful translation of input string to table IR; handles only syntactic issues
- **Semantic validation**: Happens during table IR → graph IR conversion (separate step)
- **Sanitization/normalization**: Post-parsing transformations, not part of the parser

This avoids the entanglement of syntactic/semantic processing that causes issues in RDKit.

### Four-Parser Architecture

For efficiency and clarity, we maintain four separate parsers sharing common infrastructure:

1. **`parse_smiles`** - Molecule SMILES
2. **`parse_reaction_smiles`** - Reaction SMILES (split at `>>`, parse components as SMILES)
3. **`parse_smarts`** - Substructure SMARTS (query patterns)
4. **`parse_reaction_smarts`** - Reaction SMARTS (split, parse components as SMARTS)

**Rationale**: Including SMARTS syntax in the SMILES parser would slow down the common case.
This mirrors the basic/extended parser split in MOL parsing.

### Type Hierarchy

```
Molecule              - parsed from SMILES
ExtendedMolecule      - parsed from SMARTS (query molecule)
Reaction              - container for Molecule instances + atom mapping
ExtendedReaction      - container for ExtendedMolecule instances + atom mapping
```

### Shared Infrastructure

All four parsers share:
- Tokenization (atoms, bonds, ring digits, branches)
- Bracket atom parsing (element, isotope, charge, chirality, class)
- Ring closure tracking
- Branch stack management

Parser-specific logic:
- **SMILES**: strict atom/bond validation, reject query syntax
- **SMARTS**: query atoms (`[#6]`, `[C,N]`, `[!#1]`), query bonds (`~`, `@`, `!-`), 
  recursive SMARTS (`$(...)`), logical operators (`&`, `,`, `;`, `!`)
- **Reaction SMILES**: split on `>>` or `>agent>`, validate atom mapping consistency
- **Reaction SMARTS**: same splitting, parse components as SMARTS

Previous decision: No SMIRKS support since it's been absorbed into reaction SMARTS. 
Should reconsider and include more explicit support for SMIRKS-like semantics (1:1 atom mapping).
See discussion/50-reaction-design-research-2026-02-08.md and
discussion/51-reaction-design-research-claude-2026-02-08.md for more details.

### Feature Grouping Philosophy

Feature flags are grouped logically rather than mimicking tool-specific quirks.
Tool-specific presets combine logical groups as needed for compatibility.

---

## 1. Current State

**Cleanup completed 2026-01-23.** Removed all speculative flags:
- `EXTENDED_WS`, `ALLOWS_COMMENTS`, `EXPLICIT_EOI`, `UMOL_DIALECT`, `LENIENT` - removed
- `NO_METADATA` - removed (spans always included)
- `STRICT_OPENSMILES` - retained as baseline (value 0)

Parser is now a clean OpenSMILES baseline. All 11,483 tests pass.

---

## 2. Feature Enumeration from Existing Implementations

### 2.1 RDKit SMILES Features

Source: RDKit docs, source code, test suite

**Atomic:**
- [x] Extended aromatic atoms: `se`, `as` (selenium, arsenic) - always enabled
- [x] Extended aromatic atoms: `si`, `te` (silicon, tellurium) - EXTENDED_AROMATICS flag
- [ ] Tri-character atoms: `Uu<x>` (<x> = n, b, y, q, p, h, s, o)
- [x] Wildcard atom: `*` (extended parser)
- [x] Atom class: `[C:1]` (stored as atom class; reaction semantics are separate)
- [x] Isotope: `[13C]`
- [x] Charge: `[NH4+]`, `[O-]`
- [x] Hydrogen count: `[CH4]`
- [x] Chirality: `@`, `@@`, `@TH1`, etc.
- [ ] Atom map numbers (reactions): `[C:1]` (requires reaction parser + mapping semantics)

**Bonding:**
- [x] Aromatic bonds (implicit in lowercase)
- [x] Single/double/triple: `-`, `=`, `#`
- [x] Quadruple: `$` (rare)
- [x] Aromatic: `:`
- [x] Up/down stereo: `/`, `\`
- [x] Disconnected components: `.`

**Ring closures:**
- [x] Single digit: `C1CCCCC1`
- [x] Two-digit: `C%10...%10`
- [x] Bond type on closure: `C1=CC=CC=C1` vs `c1ccccc1`

**Stereo:**
- [x] Tetrahedral: `@`, `@@`
- [x] Allene-like: `@AL1`, `@AL2`
- [x] Square planar: `@SP1`, `@SP2`, `@SP3`
- [x] Trigonal bipyramidal: `@TB1`-`@TB20`
- [x] Octahedral: `@OH1`-`@OH30`
- [x] Double bond: `/`, `\`

**Extensions (RDKit-specific):**
- [x] CXSMILES support: basic/extended CX parsing and application (Phase 6)
- [x] Dative bonds: `->`, `<-` (behind `EXTENDED_BONDS`)
- [x] Any bond: `~` (behind `EXTENDED_BONDS`, non-OpenSMILES)
- [ ] Hypervalent atoms accepted (semantic layer; parser does not do valence checking)
- [ ] Sanitization can be disabled (not a parser concern)

### 2.2 CDK/Beam SMILES Features

Source: CDK docs, Beam source

- [ ] OpenSMILES-adjacent parsing
- [ ] CXSMILES extensions (atom values, labels, coordinates)
- [ ] Preserves explicit hydrogens as-is
- [ ] Kekulization options

### 2.3 OpenBabel Features

Source: OpenBabel docs, "Universal SMILES"

- [ ] Radicals: `[CH2.]`
- [ ] InChI-based canonicalization
- [ ] Tautomer normalization
- [ ] Nitro normalization

### 2.4 Indigo Features

Source: Indigo docs, test suite

- [ ] Own dialect variations
- [ ] Reaction support
- [ ] Query features

---

## 3. Daylight SMILES Specification

The original Daylight spec (proprietary, partially documented):

**Core syntax:**
- Atoms: organic subset (B, C, N, O, P, S, F, Cl, Br, I) + brackets
- Bonds: implicit, `-`, `=`, `#`, `:`
- Branches: `()`
- Ring closures: digits 1-9, `%nn`
- Disconnection: `.`
- Aromaticity: lowercase atoms

**Stereo:**
- Tetrahedral: `@`, `@@`
- Double bond: `/`, `\`

**Extensions (later Daylight):**
- Atom class: `[atom:n]`
- Isotopes: `[nX]` where n is mass number

---

## 4. CXSMILES Specification

ChemAxon extension format: `SMILES<ws>|...|`

- The SMILES part is parsed up to the first ASCII whitespace character.
- If the remaining input (after trimming leading ASCII whitespace) starts with `|`, it is treated as
  a CX annotation block.
- Inside the `|...|` block, entries are **comma-separated**.

**Common entry types (examples):**
- Coordinates: `(x,y,z;...)`
- Atom labels: `$name;name;...$`
- Atom values: `$_AV:value;value;...$`
- Radicals: `^n:idx,idx,...`
- Lone pairs: `LP:...`, `lp:...`
- Wiggly bonds: `w:`, `wU:`, `wD:`
- Bond stereo markers: `c:`, `t:`, `ctu:`
- Local parity: `@:...`, `@@:...`
- Local bicyclo stereo: `THB:...`, `TLB:...`, `TEB:...`
- Fragment grouping: `f:...`
- Enhanced stereo: `a:...`, `o<n>:...`, `&<n>:...`
- Relative stereo: `r` (flag); `r:...` (fragment indices in reaction/multicomponent cases)
- Atom properties: `atomProp:...`
- Polymer/data/Markush: `Sg:...`, `SgD:...`, `SgH:...`, `m:...`, `RG:...`, `LOG:...`, `LO:...`
- Query-only features (CXSMARTS / query SMILES): `rb:...`, `s:...`, `u:...`

---

## 5. SMARTS Considerations

SMARTS extends SMILES with query features. Will be implemented as separate parser
(`parse_smarts`) sharing infrastructure with SMILES parser.

**Atomic primitives:**
- `*` - any atom
- `a` - aromatic
- `A` - aliphatic
- `D<n>` - degree
- `H<n>` - total H count
- `h<n>` - implicit H count
- `R<n>` - ring membership count
- `r<n>` - smallest ring size
- `v<n>` - valence
- `X<n>` - total connectivity
- `x<n>` - ring connectivity
- `#n` - atomic number
- `+<n>`, `-<n>` - charge
- `@` - chirality (in atom context)

**Logical operators:**
- `&` - and (high precedence)
- `,` - or
- `;` - and (low precedence)
- `!` - not

**Bonds:**
- `~` - any bond
- `@` - ring bond
- `!-` - not single bond (and other negated bonds)

**Recursive SMARTS:**
- `$(...)` - recursive SMARTS pattern

**Output type:** `ExtendedMolecule` (distinct from `Molecule`)

---

## 6. Feature Grouping Proposal

### Group 1: Core OpenSMILES (baseline)
- Organic subset atoms
- Standard bonds
- Branches, ring closures (1-9, %nn)
- Basic stereo (@, @@, /, \)
- Atom class, isotope, charge, H count

### Group 2: Extended Atoms ✓ COMPLETE
- Aromatic `se`, `as` (always enabled)
- Aromatic `si`, `te` (EXTENDED_AROMATICS flag)
- Wildcard `*` (extended parser)
- Tri-character `Uu<x>` (<x> = n, b, y, q, p, h, s, o) - pending, low priority

### Group 3: Extended Stereo ✓ COMPLETE (parsing)
- Allene (@AL), square planar (@SP)
- Trigonal bipyramidal (@TB), octahedral (@OH)
- Semantic validation pending (separate pass)
- Relaxed stereo validation

### Group 4: Extended Bonds ✓  COMPLETE
- Dative bonds (`->`, `<-`)
- Any bond (`~`)
- Quadruple bonds (`$`)

### Group 5: CXSMILES
- Extension block parsing `|...|`
- Coordinates, labels, radicals
- Enhanced stereo groups
- S-groups

### Group 6: Reactions (separate parser: `parse_reaction_smiles`)
- Reaction arrows `>>`
- Agent notation `>`
- Atom mapping `:n` in reactions
- Output: `Reaction` type containing `Molecule` instances

### Group 7: SMARTS (separate parser: `parse_smarts`)
- Query atom primitives (`*`, `a`, `A`, `D<n>`, `H<n>`, `#n`, `R<n>`, `r<n>`, `v<n>`, `X<n>`, `x<n>`)
- Query bonds (`~`, `@`, `!-`, `!:`, etc.)
- Logical operators (`&`, `,`, `;`, `!`)
- Recursive patterns `$(...)`
- Output: `ExtendedMolecule` type (query molecule)

### Group 8: Reaction SMARTS (separate parser: `parse_reaction_smarts`)
- Same splitting logic as reaction SMILES
- Components parsed as SMARTS
- Output: `ExtendedReaction` type containing `ExtendedMolecule` instances

---

## 7. Implementation Plan

### Phase 1: Cleanup ✓ COMPLETE
- [x] Remove UMOL_DIALECT, ALLOWS_COMMENTS, EXTENDED_WS, EXPLICIT_EOI, NO_METADATA
- [x] Remove corresponding parser code paths
- [x] Establish clean OpenSMILES baseline

### Phase 2: Feature Audit ✓ COMPLETE

Analyzed 473 failing SMILES. Breakdown by category:

**By error type:**
| Error Type | Count |
|------------|-------|
| InvalidElement | 351 |
| ConsecutiveBonds | 52 |
| InvalidBracket | 43 |
| InvalidToken | 23 |
| MismatchedRingBondDirs | 2 |
| LeadingBond | 1 |
| InvalidRingIndex | 1 |

**By root cause:**

| Category | Count | Action |
|----------|-------|--------|
| Reaction SMILES (`>>`) | 62 | Out of scope - `parse_reaction_smiles` |
| Garbage text (not SMILES) | 61 | Data quality issue |
| SMARTS syntax (`,`, `;`, `!`) | ~40 | Out of scope - `parse_smarts` |
| Wildcard `*` atom | 32 | Group 2: Extended Atoms |
| Escaped backslashes (`\\`) | ~52 | Data quality issue (double-escaped) |
| Extended ring `%(nnn)` | 1 | Low priority extension |

**Key findings:**
1. Most "InvalidElement" errors are from `*` (wildcard) - 30+ cases
2. "ConsecutiveBonds" errors are from double-escaped backslashes (`\\` in files)
3. Many failures are reaction SMILES/SMARTS, not molecule SMILES
4. Data quality issues (garbage text, escaped chars) account for ~110 cases

**Priority for SMILES parser:**
1. **Wildcard `*` support** - fixes 32 cases, needed for superatom notation
2. **Extended ring closures `%(nnn)`** - fixes 1 case, rare but valid

**Data quality fixes needed in conformance suite:**
1. Remove garbage text entries (ChemDraw:, Error., CHEMBL IDs without SMILES)
2. Fix double-escaped backslashes (`\\` should be `\`)
3. Separate reaction SMILES/SMARTS into dedicated test sets

### Phase 3: Extended SMILES Parser + Wildcards + Aromatics ✓ COMPLETE

- [x] `ExtendedMolecule` type exists and is compatible
- [x] Implement `parse_extended_smiles_bytes_with` and all wrapper variants
- [x] Add wildcard `*` support (standalone and `[*:n]` with atom class)
- [x] `From<Molecule> for ExtendedMolecule` works
- [x] Update conformance suite categories: basic_opensmiles, opensmiles, basic_chemaxon, chemaxon, chemaxon_invalid, invalid, bug
- [x] Parser hierarchy enforcement:
  - basic_opensmiles ⊆ opensmiles
  - basic_chemaxon ⊆ chemaxon
  - opensmiles ⊆ chemaxon only when there is no CX block (otherwise CX failures are classified as `chemaxon_invalid`)
- [x] Sum formula for `ExtendedMolecule` now includes wildcards (appended as `*`, `*2`, etc.)
- [x] Add aromatic `si`, `te` (silicon, tellurium) support (`se`, `as` already present)
- [x] Feature gate: `EXTENDED_AROMATICS` (applies to both parsers)

**Result:** The conformance suite and classifier use the categories listed above.

### Phase 4: Extended Stereo ✓ COMPLETE

Extended stereo types already parsed per OpenSMILES spec:
- [x] `@AL1`, `@AL2` (allenal/axial chirality)
- [x] `@SP1`, `@SP2`, `@SP3` (square planar)
- [x] `@TB1`-`@TB20` (trigonal bipyramidal)
- [x] `@OH1`-`@OH30` (octahedral)
- [x] Range validation and error handling

**Pending:** Semantic validation of stereo centers (valence, neighbor count). This is part of the semantic pass, not the parser.

### Phase 5: Extended Bonds (Group 4) ✓ COMPLETE

- [x] Dative bonds (`->`, `<-`) with `BondDonation` enum
- [x] Any bond (`~`) with `BondOrder::Any`
- [x] Quadruple bonds (`$`) - already OpenSMILES, no gating needed
- [x] Feature gate: `EXTENDED_BONDS` (gates `->`, `<-`, `~`)
- [x] `Bond::new_dative()` constructor handles `AtomPair` normalization correctly
- [x] Applied to both basic and extended parsers

### Phase 6: CXSMILES (Group 5) ✓ COMPLETE

#### 6.1 Integration point: `SMILES<ws>|...|`

The SMILES parsers stop at the first ASCII whitespace character and return the remaining input.
If the remainder (after trimming leading ASCII whitespace) starts with `|` and
`SmilesParseFlags::CHEMAXON_EXTENSIONS` is enabled, the CX block is parsed and applied to the IR.
Otherwise the remainder is ignored (per the “data after whitespace is ignored” convention).

This yields four “dialects” that are used consistently across classification and conformance:

- `basic_opensmiles`: SMILES only, no wildcards, ignore CX suffix
- `opensmiles`: SMILES + wildcards, ignore CX suffix
- `basic_chemaxon`: SMILES + basic CX parsing
- `chemaxon`: SMILES + wildcards + extended CX parsing

The classifier/conformance suite also uses:
- `chemaxon_invalid`: SMILES part parses, but CX block is invalid/unhandled in strict CX mode

#### 6.2 CX block parsing: `CxEntry` + basic/extended split

`io/smiles/parser/cx.rs` defines a `CxEntry` enum and two block parsers:

- `parse_cx_annotations(input: &[u8], flags: SmilesParseFlags) -> Result<Vec<CxEntry>, ParseError>`
  - Accepts the basic subset.
  - Returns `InvalidCxTag` on any CX tag that is not part of the basic subset (including extended-only tags),
    unless `flags` contains `SKIP_UNKNOWN_CHEMAXON_TAGS` (then the tag is ignored).
- `parse_extended_cx_annotations(input: &[u8], flags: SmilesParseFlags) -> Result<Vec<CxEntry>, ParseError>`
  - Accepts the full subset implemented so far.

Both parsers:
- split entries on commas inside `|...|`
- if `flags` contains `SKIP_UNKNOWN_CHEMAXON_TAGS`: skip unknown/unrecognized entries (consume up to the next
  entry boundary)
- do **not** skip malformed *known* tags, even with `SKIP_UNKNOWN_CHEMAXON_TAGS` enabled
- support HTML entity escapes for labels/values

Coordinate parsing notes:
- `()` means “no atoms have coordinates” (empty vector)
- missing components are allowed (`(x,,)`, `(,y,)`, `(,,z)`, `(,,)` etc.)
- 4D coordinates are rejected as a hard failure

#### 6.3 Applying CX entries to TableIR

Two application functions exist:

- `update_molecule(&mut Molecule, entries: Vec<CxEntry>)`
  - applies positions/labels/values/radicals/bond markers/bond annotations
  - ignores extended-only entries
- `update_extended_molecule(&mut ExtendedMolecule, entries: Vec<CxEntry>)`
  - applies the same per-atom/per-bond data using extended atom/bond types
  - sets `ExtendedMolecule.stereo_interpretation: Option<StereoInterpretation>` directly from:
    - `a:` → `Absolute`
    - `r` → `Relative`
  - populates `ExtendedMolecule.cx_data: Option<CxAnnotationData>` when any of:
    - `stereo_groups: HashMap<u32, StereoSet>` (from `o<n>:` and `&<n>:`)
    - `components: Option<Vec<Vec<u32>>>` (from `f:`)

Out-of-range atom/bond indices are **hard errors** (no skipping):
- `AtomIndexOutOfBounds`
- `BondIndexOutOfBounds`
- `MismatchedAtomBondIndices` (bond-indexed tags where the referenced atom is not incident on the bond)

#### 6.4 CXSMILES completeness + correctness work plan

The implementation above covers the “common” CX blocks seen in typical datasets.
ChemAxon’s CXSMILES spec includes additional tags, and a few tags we parse today need a
spec-correct interpretation.

**Correctness fixes (implemented):**

- `C:` / `H:`:
  - **Spec**: list of `<first_atom_idx>.<bond_idx>` pairs referring to an existing bond in the SMILES part.
  - **Implemented**: apply by mutating the referenced bond (set `donation` for `C:`, set
    `noncovalent=Hydrogen` and `order=Zero` for `H:`).
- `w:` / `wU:` / `wD:`:
  - **Spec**: list of `<atom_idx>.<bond_idx>` pairs.
  - **Implemented**: apply wedge to the referenced bond. (TableIR stores wedges relative to the first
    atom of the normalized `AtomPair`; we may be unable to preserve the “wedge endpoint” if the
    CX atom index is not that first atom.)
- `ctu:`:
  - **Spec**: bond indices with UNSPEC ring double-bond stereo.
  - **Implemented**: apply as `BondStereo::Either` for the listed bonds.
- `r` / `r:...`:
  - **Spec**: `r` marks relative stereoconfiguration; `r:...` lists fragment indices with relative configuration
    in reaction/multicomponent cases.
  - **Implemented**: accept standalone `r` for molecule parsing; reject `r:...` as `InvalidCxTag` (not meaningful
    for `ExtendedMolecule` parsing).

**Added CXSMILES tags (implemented):**

- `Sg:` / `SgD:` / `SgH:`:
  - **MOL-equivalent**: yes (maps naturally to `ExtendedMolecule.ctfile_data.sgroups: BTreeMap<u32, SGroup>`).
  - **Implemented**: parse S-groups in CX order, assign sequential S-group indices, store into `ctfile_data`.
  - **Stored**: Sg: type, subtype, atoms, subscript, connectivity, bond_indices (head+tail), bracket_orientation, bracket_style, bracket_coords (2 or 4 points), connectivity_flip. SgD: atoms, field_name, data_content, query_operator, field_units, query_identifier, bracket_coords (when numeric). SgH: hierarchy_parent on child S-groups.
  - **For discussion**:
    - SgD coords `(-1)` — atom-attached: we do not set `bracket_coords` and do not store that it was `(-1)` vs absent. Round-trip or display may need an explicit flag.
    - `SGroupBracketStyle::TypeR` and `TypeS` (CX bracket types `r`/`s`): exact semantics TBD; stored as-is.
    - Sg type keywords for Superatom (SUP) and MultipleGroup (MUL): not mapped; unknown keywords rejected.
    - Connectivity flip format: parsed as `ht,1` / `ht,flip` / `ht,0`; ChemAxon format not confirmed.
    - SgD multi-line `data_content`: spec allows it; we use `Some(vec![s])`; multi-line would need different delimiter/escaping.
    - `take_until_entry_boundary`: now requires `:` after tag name (e.g. `,Sg:`) so bracket coords `s,b,1,2,3,4` are not split at `,b`; any CX tag with comma+letter without `:` could be affected.
- `LN:` link nodes:
  - **Mostly MOL-equivalent**: partially (closest match is `ExtendedAtom.link_atom: Option<LinkAtom>`).
  - **Implemented**: parse min/max repetition + outer atoms; store in `ExtendedAtom.link_atom`.
- `LO:` ligand order:
  - **MOL-equivalent**: partially (port/ligand ordering metadata).
  - **Implemented**: store ordered neighbor list in `ExtendedAtom.ligand_order`.
- `rb:` / `s:` / `u:` query features:
  - **MOL-equivalent**: yes (maps to `ExtendedAtom.ring_bond_count`, `substitution_count`, `unsaturated`).
  - **Implemented**: apply to `ExtendedAtom.ring_bond_count`, `substitution_count`, `unsaturated`.
- `LP:` / `lp:` lone electron pairs:
  - **MOL-equivalent**: partially (`AtomSymbol::LonePair` exists, but the CX tags also encode counts).
  - **Implemented**: store count in `ExtendedAtom.lone_pairs`.
- `m:` multicenter bonds / variable attachment:
  - **Implemented**: parse and store in `ExtendedMolecule.multicenter_bonds` as `MulticenterBond` (center + ligand sets).
- `THB:` / `TLB:` / `TEB:` (bicyclo stereo):
  - **Implemented**: parse and store in `ExtendedAtom.bicyclo_stereo` as `BicycloStereoData`.

**Remaining CXSMILES tags (todo):**

- `RG:` / `LOG:` Markush:
  - **LOG** is MOL-equivalent (maps to `cx_data.rgroups` metadata).
      The `RGroup` struct has: `label`, `dependent_label`, `rgroup_or_h`, `occurrence`
  - **RG** requires additional representation for member structures (embedded `{...}` CXSMILES blocks).
      CXSMILES embeds member structures as `{...}` blocks (SMILES strings). Store as raw strings for roundtripping.
  - **Status**: Removed from parser (code quality insufficient); pending reimplementation.
- `@:` / `@@:` (Local parity):
  - CXSMILES spec is ambiguous, no examples found.
  - **Status**: Implementation postponed until format can be established

#### 6.5 Conformance + classification details

Both the conformance driver and the `classify_smiles_strings` tool treat “has CX annotations”
as an input property, not a parser outcome. The detection uses this regex:

`^\S+\s+\|.*\|`

This distinguishes “plain SMILES” from “SMILES + CX suffix” even when OpenSMILES parsers ignore the suffix.

Conformance tests are gated behind the `conformance` Cargo feature and do not run under plain `cargo test`.
Run with `cargo test -p umol-models-graph --features conformance --test smiles_parsing` (and similarly for `mol_parsing` / `sdf_parsing`).

### Phase 7: Presets ✓ COMPLETE

Parsing presets are exposed via `SmilesIoConfig` constructors and map to `SmilesParseFlags`:

- `SmilesIoConfig::basic_opensmiles()` → `BASIC_OPENSMILES`
- `SmilesIoConfig::opensmiles()` → `OPENSMILES`
- `SmilesIoConfig::basic_chemaxon()` → `BASIC_CHEMAXON`
- `SmilesIoConfig::chemaxon()` → `CHEMAXON`

Additional presets exist for testing/debugging (`basic_max`, `extended_max`, `strict`, `extended`, `lenient`).

### Phase 8: Reaction SMILES Parser (Group 6)

Design research: [50-reaction-design-research-2026-01-29.md](50-reaction-design-research-2026-01-29.md)

- [x] Create `Reaction` type (container for `Molecule` + atom mapping)
- [x] Implement `parse_reaction_smiles`
- [x] Split on `>>` (reactants >> products) and `>` (agents)
- [x] Parse components using `parse_smiles`
- [ ] Validate atom mapping consistency

### Phase 9: SMARTS Parser
- [ ] Extend `ExtendedMolecule` for SMARTS query primitives (if needed)
- [ ] Implement `parse_smarts_bytes_with` and wrappers
- [ ] Implement query atom primitives (`a`, `A`, `D<n>`, `H<n>`, `#n`, `R<n>`, etc.)
- [ ] Implement query bonds (`~`, `@`, `!-`)
- [ ] Implement logical operators (`&`, `,`, `;`, `!`)
- [ ] Implement recursive SMARTS (`$(...)`)

Note: Builds on `ExtendedMolecule` from Phase 3. Wildcard `*` already supported.

### Phase 10: Reaction SMARTS Parser
- [ ] Create `ExtendedReaction` type (container for `ExtendedMolecule` + atom mapping)
- [ ] Implement `parse_reaction_smarts`
- [ ] Split on `>>` and `>`
- [ ] Parse components using `parse_smarts`

---

## 8. Resolved Questions

1. **Should SMARTS use the same lexer as SMILES?**
   → Separate parsers sharing common infrastructure. Including SMARTS syntax would slow 
   down SMILES parsing. Same approach as basic/extended MOL parsers.

2. **How to handle toolkit-specific quirks (e.g., RDKit sanitization)?**
   → Sanitization and semantic processing are out of scope for the parser. Parser produces
   faithful table IR; semantic validation happens during table IR → graph IR conversion.

3. **Should we support SMIRKS?**
   → No. SMIRKS and reaction SMARTS have converged. All major toolkits use reaction SMARTS.
   No significant data sources provide SMIRKS exclusively.

4. **How granular should feature gates be?**
   → Logical grouping (extended atoms, extended stereo, etc.) combined into presets.
   Tool-specific presets combine logical groups as needed, not the other way around.

## 9. Remaining Questions

1. **ExtendedMolecule structure for SMARTS** - Postponed until after SMILES implementation.
   Recursive SMARTS (`$(...)`) will require careful design.

## 10. Preliminary Designs

### Wildcard and Extended Atom Handling

**Decision:** SMILES with wildcards (`*`, `[*]`, `[*:n]`) parse into `ExtendedMolecule`, not `Molecule`.

**Rationale:**
- Type safety: `Molecule` implies fully specified atoms; wildcards break downstream calculations
- Alignment: Matches MOL parsing where wildcards are extended features
- Semantic clarity: Both "unknown" and "query" atoms require identical handling

**API pattern (aligned with MOL/SDF parsers):**
```rust
// Basic SMILES → Molecule (fails on wildcards)
parse_smiles_bytes_with(input: &[u8], config: &SmilesIoConfig) -> Result<Molecule, ParseError>
parse_smiles_bytes(input: &[u8]) -> Result<Molecule, ParseError>
parse_smiles_with(input: &str, config: &SmilesIoConfig) -> Result<Molecule, ParseError>
parse_smiles(input: &str) -> Result<Molecule, ParseError>

// Extended SMILES → ExtendedMolecule (handles wildcards, query atoms)
parse_extended_smiles_bytes_with(input: &[u8], config: &SmilesIoConfig) -> Result<ExtendedMolecule, ParseError>
parse_extended_smiles_bytes(input: &[u8]) -> Result<ExtendedMolecule, ParseError>
parse_extended_smiles_with(input: &str, config: &SmilesIoConfig) -> Result<ExtendedMolecule, ParseError>
parse_extended_smiles(input: &str) -> Result<ExtendedMolecule, ParseError>
```

The `_bytes_with` variants are the low-level implementations; others layer on top.

**Goal:** `ExtendedMolecule` should be a superset of `Molecule` (`From<Molecule> for ExtendedMolecule`).
This allows code that works on extended molecules to accept plain molecules.
May require careful design for recursive SMARTS features.

### Atom Mapping for Reactions

Maps atom class (mapping number) to atom indices in reactant/product atom tables:

```rust
// Simple case: at most one occurrence per side
type AtomMapping = BTreeMap<u32, (Option<u32>, Option<u32>)>;
// key: atom class, value: (reactant_atom_idx, product_atom_idx)

// If multiple occurrences per side are needed:
type AtomMapping = BTreeMap<u32, (Vec<u32>, Vec<u32>)>;
// key: atom class, value: (reactant_atom_indices, product_atom_indices)
```

Bond mapping (useful for reaction templates) can follow the same pattern.

Indices refer to positions in the atom table IR (u32 or usize).

---

## 11. SMILES Extension Analysis

This section catalogs SMILES extensions from various sources, classified by priority and parser category.

### 11.1 Classification Criteria

**Priority levels:**
- **P1 (High)**: Required for interoperability with major tools/databases
- **P2 (Medium)**: Useful for advanced use cases, common in research workflows
- **P3 (Low)**: Rare, niche, or can be approximated by other means

**Parser category:**
- **Basic**: Extends `parse_smiles` → `Molecule`
- **Extended**: Extends `parse_extended_smiles` → `ExtendedMolecule`
- **Reaction**: For `parse_reaction_smiles` / `parse_reaction_smarts`
- **SMARTS**: For `parse_smarts` → `ExtendedMolecule` (query)
- **CXSMILES**: Requires extension block parsing

### 11.2 Open-Source Implementation Features

#### RDKit

| Feature | Status | Priority | Category | Notes |
|---------|--------|----------|----------|-------|
| Wildcard `*` | ✓ Done | P1 | Extended | Standalone and `[*:n]` |
| Aromatic `se`, `as` | ✓ Done | P2 | Basic | Selenium, arsenic (always enabled) |
| Aromatic `si`, `te` | ✓ Done | P2 | Basic | Silicon, tellurium (EXTENDED_AROMATICS flag) |
| Tri-character `Uun` | Pending | P2 | Basic | Obsolete symbols for Ds--Og |
| Dative bonds `->` `<-` | ✓ Done | P2 | Basic | Behind `EXTENDED_BONDS` |
| Quadruple bond `$` | Done | P3 | Basic | Very rare (Mo-Mo) |
| Any bond `~` | ✓ Done | P3 | Basic | Behind `EXTENDED_BONDS` (non-OpenSMILES) |
| Extended stereo `@AL`, `@SP`, `@TB`, `@OH` | ✓ Parsed | P2 | Basic | Semantic validation pending |
| CXSMILES extension block | ✓ Done | P1 | CXSMILES | Parsed + applied when enabled |
| Enhanced stereo groups | ✓ Parsed | P2 | CXSMILES | `a:`, `o<n>:`, `&<n>:` |
| Radicals (CXSMILES) | ✓ Parsed | P2 | CXSMILES | `^1:` through `^7:` |
| Atom map numbers `:n` in reactions | Pending | P1 | Reaction | |
| Hypervalent atoms | N/A | - | - | Semantic layer, not parser |
| Sanitization disable | N/A | - | - | Semantic layer |

#### OpenBabel

| Feature | Status | Priority | Category | Notes |
|---------|--------|----------|----------|-------|
| External bond `&` | Skip | P3 | Extended | Fragment attachment points (`CC&1.C&1`) |
| Radicals `[CH2.]` | Pending | P2 | Basic | Dot notation for unpaired e- |
| Extended ring `%(nnn)` | Pending | P3 | Basic | 3+ digit ring indices |
| Atom typing | N/A | - | - | Semantic layer |

**External bond `&`:** OpenBabel-specific syntax for combinatorial chemistry fragments.
`CC&1.C&1` connects two fragments at `&1` positions. Creates dummy atoms for unmatched
attachment points. Low priority - no other toolkit supports this syntax.

#### CDK/Beam

| Feature | Status | Priority | Category | Notes |
|---------|--------|----------|----------|-------|
| OpenSMILES-compliant | ✓ Baseline | P1 | Basic | Via Beam library |
| CXSMILES extensions | Pending | P1 | CXSMILES | Full support (see SmiFlavor) |
| Extended stereo | Pending | P2 | Basic | SP, TB, OH in SmiFlavor |
| Explicit H preservation | N/A | - | - | Already preserved in IR |

**SmiFlavor.java** (output flags, implies parsing support):
- `StereoSquarePlanar`, `StereoTrigonalBipyramidal`, `StereoOctahedral` - extended stereo
- `CxAtomLabel`, `CxAtomValue`, `CxRadical`, `CxMulticenter`, `CxPolymer` - CXSMILES
- `CxEnhancedStereo`, `CxLigandOrder`, `CxDataSgroups` - advanced CXSMILES
- `InChILabelling` - Universal SMILES (with caveats about canonicalization)

#### Indigo

| Feature | Status | Priority | Category | Notes |
|---------|--------|----------|----------|-------|
| CXSMILES parsing | Pending | P1 | CXSMILES | Full `\|...\|` block support |
| Query SMILES/SMARTS | Pending | P2 | SMARTS | Separate `loadSMARTS()` method |
| Reaction support | Pending | P1 | Reaction | `loadMolecule()` handles reactions |
| R-group notation | Pending | P2 | Extended | `[R1]`, `[R2]` |
| Polymer notation | Pending | P2 | CXSMILES | `{...}n` curly brace syntax |

**Note:** `smiles_loader_parsers.cpp` primarily parses CXSMILES extension blocks, not
Indigo-specific extensions. Features include: `w:` (wiggly stereo), `a:` (absolute),
`o1:/&1:` (OR/AND groups), `^n:` (radicals), `$...$` (pseudoatoms), `c:` (coords),
`Sg:` (S-groups).

### 11.3 Daylight/ChemAxon Documentation

#### Daylight (Original Spec)

| Feature | Status | Priority | Category | Notes |
|---------|--------|----------|----------|-------|
| Core SMILES | ✓ Done | P1 | Basic | OpenSMILES baseline |
| Vector binding `&&` | Skip | P3 | - | Obsolete |
| Unique SMILES | N/A | - | - | Canonicalization, not parsing |

#### ChemAxon (CXSMILES)

| Feature | Status | Priority | Category | Notes |
|---------|--------|----------|----------|-------|
| Extension block `\|...\|` | ✓ Parsed | P1 | CXSMILES | Primary container |
| 2D/3D coordinates `(x,y,z;...)` | ✓ Parsed | P1 | CXSMILES | Essential for structure |
| Atom labels `$name$` | ✓ Parsed | P2 | CXSMILES | Display names |
| Atom values `$_AV:...$` | ✓ Parsed | P2 | CXSMILES | Atom values |
| Radicals `^1:` through `^7:` | ✓ Parsed | P2 | CXSMILES | With spin multiplicity |
| Enhanced stereo `a:`, `o<n>:`, `&<n>:` | ✓ Parsed | P2 | CXSMILES | Absolute/OR/AND groups |
| Wiggly bonds `w:`, `wU:`, `wD:` | ✓ Parsed + applied | P2 | CXSMILES | Indexed by bond in spec |
| CIS/TRANS `c:`, `t:`, `ctu:` | ✓ Parsed + applied | P2 | CXSMILES | Ring double bond stereo |
| Fragment groups `f:` | ✓ Parsed | P2 | CXSMILES | Fragment grouping |
| Atom properties `atomProp:` | ✓ Parsed | P3 | CXSMILES | Generic properties |
| Coordinate bonds `C:` | ✓ Parsed + applied | P3 | CXSMILES | Indexed by bond in spec |
| Hydrogen bonds `H:` | ✓ Parsed + applied | P3 | CXSMILES | Indexed by bond in spec |
| Relative stereo `r` | ✓ Parsed | P3 | CXSMILES | Flag-only (`r`); `r:...` (fragment indices) rejected in molecule parsing |
| Lone pairs `LP:`, `lp:` | ✓ Parsed | P3 | CXSMILES | Electron pairs / explicit counts |
| Local parity `@:`, `@@:` | Pending | P3 | CXSMILES | Additional stereo metadata |
| Local bicyclo stereo `THB:`, `TLB:`, `TEB:` | ✓ Parsed | P3 | CXSMILES | Bicyclo stereo metadata |
| S-groups `Sg:` | ✓ Parsed | P2 | CXSMILES | Polymers, abbreviations |
| Data S-groups `SgD:` | ✓ Parsed | P3 | CXSMILES | Embedded data |
| S-group hierarchy `SgH:` | ✓ Parsed | P3 | CXSMILES | Parent-child relations (CX order) |
| Ligand order `LO:` | ✓ Parsed | P3 | CXSMILES | Ligand/port ordering metadata |
| Link nodes `LN:` | ✓ Parsed | P3 | CXSMILES | Link nodes (repeat ranges) |
| Ring bond count `rb:` | ✓ Parsed | P3 | CXSMILES | Query feature (MDL query) |
| Substitution count `s:` | ✓ Parsed | P3 | CXSMILES | Query feature (MDL query) |
| Unsaturation `u:` | ✓ Parsed | P3 | CXSMILES | Query feature (MDL query) |
| Multicenter bonds `m:` | ✓ Parsed | P3 | CXSMILES | Variable attachment / organometallics |
| R-groups `RG:` | Pending | P2 | CXSMILES | Markush definitions (removed; pending reimplementation) |
| R-logic `LOG:` | Pending | P2 | CXSMILES | Markush logic (removed; pending reimplementation) |

### 11.4 Uncommon but Interesting Variants

#### DeepSMILES

| Feature | Assessment | Priority | Notes |
|---------|------------|----------|-------|
| Ring closure pairs | Interesting | P3 | Removes ambiguity in ring closures |
| Branch depth encoding | Interesting | P3 | `))` instead of `))` with matching `(` |
| ML-friendly design | - | - | Designed for sequence models |

**Verdict:** Potentially useful as alternative *output* format for ML applications. Not priority for input parsing since data sources don't use it.

#### SELFIES

| Feature | Assessment | Priority | Notes |
|---------|------------|----------|-------|
| 100% syntactically valid | Claimed | - | But semantically invalid molecules common |
| Token-based design | - | - | For RNN/Transformer input |

**Verdict:** Skip. The "100% valid" claim is misleading - produces syntactically valid but chemically nonsensical structures. No significant data sources use it.

#### Jmol Extensions

| Feature | Assessment | Priority | Notes |
|---------|------------|----------|-------|
| 3D coordinate embedding | Interesting | P3 | Inline coordinates |
| Animation directives | Skip | - | Visualization-specific |

**Verdict:** Low priority. Jmol's SMILES extensions are mostly for visualization, not structure interchange.

#### TwistSMILES / CanonSMILES

| Feature | Assessment | Priority | Notes |
|---------|------------|----------|-------|
| Canonical forms | - | - | Canonicalization, not parsing |

**Verdict:** N/A for parser. These are about canonical output, not input parsing.

#### Extended Connectivity Fingerprints (ECFP) in SMILES

Some tools embed ECFP-like atom environments in SMILES-like notation. Not standard; skip.

### 11.5 Commercial Tool Outputs

#### ChemDraw

| Feature | Observed | Priority | Notes |
|---------|----------|----------|-------|
| Standard SMILES | Yes | P1 | Usually OpenSMILES-compatible |
| Abbreviated groups | Sometimes | P2 | `Ph`, `Me`, `Et` as pseudoatoms |
| Stereochemistry | Yes | P1 | Standard notation |

#### ChemDoodle

| Feature | Observed | Priority | Notes |
|---------|----------|----------|-------|
| Standard SMILES | Yes | P1 | OpenSMILES-compatible |
| 3D export | Optional | P2 | Via CXSMILES or mol2 |

#### Marvin/ChemAxon Tools

| Feature | Observed | Priority | Notes |
|---------|----------|----------|-------|
| CXSMILES | Default | P1 | Full CXSMILES support |

### 11.6 Prioritized Feature List

**P1 (High Priority) - Required for interoperability:**

1. ✓ Wildcard atoms `*` (Done)
2. ✓ CXSMILES extension block parsing `|...|` (Parsed)
3. ✓ CXSMILES coordinates `(x,y,z;...)` (Parsed)
4. Reaction SMILES `>>` splitting
5. Atom mapping in reactions `[C:1]`

**P2 (Medium Priority) - Common advanced features:**

1. ✓ Aromatic `se`, `te` (Done, `b` already basic)
2. Dative bonds `->` `<-`
3. ✓ Radicals CXSMILES `^n:` (Parsed); OpenBabel dot notation pending
4. ✓ CXSMILES atom labels `$name$` (Parsed)
5. ✓ CXSMILES enhanced stereo groups (Parsed)
6. ✓ CXSMILES S-groups (basic: abbreviations, polymers) (Parsed)
7. SMARTS query primitives
8. R-group notation `[R1]`

**P3 (Low Priority) - Rare or niche:**

1. Quadruple bonds `$`
2. Extended ring closures `%(nnn)`
3. ✓ CXSMILES link nodes, data S-groups (Parsed)
4. DeepSMILES (as output format only)
5. ✓ Multicenter bonds (Parsed)

**Skip (Not implementing):**

1. SELFIES (misleading validity claims)
2. Jmol animation directives
3. Daylight vector binding `&&`
4. OpenBabel external bonds `&` (no other toolkit supports)
5. Any bond `~` in SMILES (SMARTS only; CXSMILES `|Z:|` is different)

---

## 12. References

- OpenSMILES spec: http://opensmiles.org/
- Daylight theory manual (archived): http://www.daylight.com/dayhtml/doc/theory/
- ChemAxon CXSMILES docs: https://docs.chemaxon.com/display/docs/chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.md
- RDKit SMILES parsing docs: https://www.rdkit.org/docs/RDKit_Book.html
- OpenBabel SMILES: http://openbabel.org/wiki/SMILES
- CDK/Beam: https://github.com/cdk/cdk (Beam integrated)
- Indigo: https://lifescience.opensource.epam.com/indigo/
- DeepSMILES paper: https://doi.org/10.26434/chemrxiv.7097960
- Depth-First blog (ChemCore comparisons)
