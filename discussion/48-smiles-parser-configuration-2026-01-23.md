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

### No SMIRKS Support

SMIRKS and reaction SMARTS have converged. Both use atom mapping, and reaction SMARTS 
covers the same ground with clearer semantics. RDKit, CDK, and OpenBabel all use 
reaction SMARTS. No major database provides SMIRKS exclusively.

### Feature Grouping Philosophy

Feature flags are grouped logically rather than mimicking tool-specific quirks.
Tool-specific presets combine logical groups as needed for compatibility.

---

## 1. Current State (COMPLETED)

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
- [ ] Wildcard atom: `*`
- [ ] Atom class: `[C:1]`
- [ ] Isotope: `[13C]`
- [ ] Charge: `[NH4+]`, `[O-]`
- [ ] Hydrogen count: `[CH4]`
- [ ] Chirality: `@`, `@@`, `@TH1`, etc.
- [ ] Atom map numbers (reactions): `[C:1]`

**Bonding:**
- [ ] Aromatic bonds (implicit in lowercase)
- [ ] Single/double/triple: `-`, `=`, `#`
- [ ] Quadruple: `$` (rare)
- [ ] Aromatic: `:`
- [ ] Up/down stereo: `/`, `\`
- [ ] Disconnected components: `.`

**Ring closures:**
- [ ] Single digit: `C1CCCCC1`
- [ ] Two-digit: `C%10...%10`
- [ ] Bond type on closure: `C1=CC=CC=C1` vs `c1ccccc1`

**Stereo:**
- [ ] Tetrahedral: `@`, `@@`
- [ ] Allene-like: `@AL1`, `@AL2`
- [ ] Square planar: `@SP1`, `@SP2`, `@SP3`
- [ ] Trigonal bipyramidal: `@TB1`-`@TB20`
- [ ] Octahedral: `@OH1`-`@OH30`
- [ ] Double bond: `/`, `\`

**Extensions (RDKit-specific):**
- [ ] CXSMILES support (partial): coordinates, atom labels, radicals
- [ ] Dative bonds: `->`, `<-`
- [ ] Query bonds: `~`
- [ ] Hypervalent atoms accepted
- [ ] Sanitization can be disabled

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

ChemAxon extension format: `SMILES |ext1;ext2;...|`

**Extension types:**
- `c:coords` - 2D/3D coordinates
- `atomProp:` - atom properties
- `$name$` - atom labels
- `r:n,radical` - radicals
- `^n:atoms` - enhanced stereo (absolute, or, and groups)
- `Sg:` - S-groups (polymers, etc.)
- `rb:` - ring bond count
- `s:` - substitution count
- `u:` - unsaturation
- `o:` - atom ordering (for canonical output)

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

### Group 4: Extended Bonds
- Dative bonds (`->`, `<-`)
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
- [x] Update conformance suite with three categories: basic_opensmiles, opensmiles, invalid
- [x] Parser hierarchy enforcement: extended parser must be superset of basic parser
- [x] Classifier and test suite both verify hierarchy invariant
- [x] Sum formula for `ExtendedMolecule` now includes wildcards (appended as `*`, `*2`, etc.)
- [x] Add aromatic `si`, `te` (silicon, tellurium) support (`se`, `as` already present)
- [x] Feature gate: `EXTENDED_AROMATICS` (applies to both parsers)

**Result:** 32 SMILES now parse as `opensmiles` category (wildcards); 8979 as `basic_opensmiles`

### Phase 4: Extended Stereo ✓ COMPLETE (parsing)

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

### Phase 6: CXSMILES (Group 5)

#### 6.1 Architecture

The CX annotation block (`|...|`) uses an **accumulator pattern** similar to CTab M-lines.
A `CxAccumulator` collects properties, then `update_molecule` / `update_extended_molecule`
applies them to the IR.

Two separate parsers (like CTab `properties_block` / `extended_properties_block`):
- `parse_cx_annotations` - basic properties only, errors on extended
- `parse_extended_cx_annotations` - all properties

**Basic annotations (properties of basic molecules):**
- Coordinates `(x,y,z;...)`
- Radicals `^1:`, `^2:`, etc.
- Atom labels `$...;...$`
- Atom values `$_AV:...$`
- Wiggly bonds `w:`, `wU:`, `wD:` (undefined stereo indicator)
- CIS/TRANS markers `c:`, `t:`, `ctu:` (for double bonds in rings)

**Extended annotations (require ExtendedMolecule):**
- All basic annotations plus:
- Fragment grouping `f:`
- Enhanced stereo `a:`, `o<n>:`, `&<n>:`
- Atom properties `atomProp:`
- Coordinate/dative bonds `C:`
- Hydrogen bonds `H:`
- Relative stereo flag `r`

**Deferred (CTab legacy):**
- Pseudo/special atoms (`*_p`, `Q_e`, etc.)
- S-groups, R-groups, link nodes
- MDL query features (`rb:`, `s:`, `u:`)
- Ligand order, bicycloalkane stereo, lone pairs

#### 6.2 IR Changes

- [ ] Add `stereo_groups: Vec<StereoGroup>` to `ExtendedMolecule` only
- [ ] Define `StereoGroup { group_type: StereoGroupType, atoms: Vec<u32> }`
- [ ] Define `StereoGroupType { Absolute, Or(u32), And(u32) }`

#### 6.3 Parser Implementation

**Step 1: Create module and accumulator**
- [ ] Create `io/smiles/parser/cx.rs` module
- [ ] Define `CxAccumulator` struct:
```rust
#[derive(Debug, Default)]
pub(super) struct CxAccumulator {
    // Per-atom (indexed by atom position in SMILES)
    pub atom_labels: BTreeMap<u32, String>,
    pub atom_values: BTreeMap<u32, String>,
    pub atom_radicals: BTreeMap<u32, u8>,
    pub atom_properties: BTreeMap<u32, HashMap<String, String>>,
    
    // Per-bond (indexed by bond position or atom pair)
    pub wiggly_bonds: BTreeMap<u32, BondWedge>,
    pub coordinate_bonds: Vec<(u32, u32)>,
    pub hydrogen_bonds: Vec<(u32, u32)>,
    pub cis_bonds: Vec<u32>,
    pub trans_bonds: Vec<u32>,
    
    // Molecule-level
    pub coordinates: Option<Vec<Point3D>>,
    pub stereo_groups: Vec<StereoGroup>,
    pub relative_stereo: bool,
    pub fragment_groups: Vec<Vec<u32>>,
}
```

**Step 2: Tokenizer**
- [ ] Tokenize on `,` respecting nested `()`, `{}`, `$...$`
- [ ] Handle `&#code;` unescaping

**Step 3: Two parser functions**
- [ ] `parse_cx_annotations(input: &[u8]) -> Result<CxAccumulator, ParseError>`
  - Handles basic tags only
  - Returns error on extended-only tags
- [ ] `parse_extended_cx_annotations(input: &[u8]) -> Result<CxAccumulator, ParseError>`
  - Handles all tags

**Step 4: Tag parsers**
Basic tags:
- [ ] Coordinates `(x,y,z;...)`
- [ ] Radicals `^1:` through `^7:`
- [ ] Atom labels `$...;...$`
- [ ] Atom values `$_AV:...$`
- [ ] Wiggly bonds `w:`, `wU:`, `wD:`
- [ ] CIS/TRANS `c:`, `t:`, `ctu:`

Extended-only tags:
- [ ] Fragment grouping `f:`
- [ ] Enhanced stereo `a:`, `o<n>:`, `&<n>:`
- [ ] Atom properties `atomProp:`
- [ ] Coordinate bonds `C:`
- [ ] Hydrogen bonds `H:`
- [ ] Relative stereo `r` or `r:idx,...`

**Step 5: Apply to IR**
- [ ] `update_molecule(&mut Molecule)` - basic features only
- [ ] `update_extended_molecule(&mut ExtendedMolecule)` - all features
- [ ] Validate atom/bond indices against molecule size

**Step 6: Integration with SMILES parser**
- [ ] Add `CXSMILES` flag to `SmilesParseFlags`
- [ ] Detect `|...|` suffix after SMILES string
- [ ] Basic parser calls `parse_cx_annotations`
- [ ] Extended parser calls `parse_extended_cx_annotations`

#### 6.4 Testing
- [ ] Unit tests for each tag parser
- [ ] Integration tests with real CXSMILES from RDKit/CDK test suites
- [ ] Conformance suite: add `cxsmiles` category

### Phase 7: Presets
- [ ] `basic_opensmiles` - Group 1 only (current strict baseline)
- [ ] `opensmiles` - Groups 1-4 (with wildcards, extended aromatics)
- [ ] `cxsmiles` - Groups 1-5
- [ ] Tool-specific presets as needed (rdkit_compat, etc.)

### Phase 8: Reaction SMILES Parser
- [ ] Create `Reaction` type (container for `Molecule` + atom mapping)
- [ ] Implement `parse_reaction_smiles`
- [ ] Split on `>>` (reactants >> products) and `>` (agents)
- [ ] Parse components using `parse_smiles`
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
| Dative bonds `->` `<-` | Pending | P2 | Basic | Metal complexes |
| Quadruple bond `$` | Done | P3 | Basic | Very rare (Mo-Mo) |
| Any bond `~` | N/A | - | SMARTS | Query bond only; in SMILES only via CXSMILES `\|Z:\|` |
| Extended stereo `@AL`, `@SP`, `@TB`, `@OH` | ✓ Parsed | P2 | Basic | Semantic validation pending |
| CXSMILES extension block | Pending | P1 | CXSMILES | Coordinates, labels, radicals |
| Enhanced stereo groups | Pending | P2 | CXSMILES | `^1:`, `^2:` |
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
| Extension block `\|...\|` | Pending | P1 | CXSMILES | Primary container |
| 2D/3D coordinates `c:` | Pending | P1 | CXSMILES | Essential for structure |
| Atom labels `$name$` | Pending | P2 | CXSMILES | Display names |
| Radicals `^1:`, `^2:` | Pending | P2 | CXSMILES | Radical notation |
| Enhanced stereo `&1:`, `o1:` | Pending | P2 | CXSMILES | AND/OR groups |
| S-groups `Sg:` | Pending | P2 | CXSMILES | Polymers, abbreviations |
| Atom properties `atomProp:` | Pending | P3 | CXSMILES | Generic properties |
| Ring bond count `rb:` | Skip | P3 | CXSMILES | Query feature |
| Substitution count `s:` | Skip | P3 | CXSMILES | Query feature |
| Unsaturation `u:` | Skip | P3 | CXSMILES | Query feature |
| Atom ordering `o:` | Pending | P3 | CXSMILES | Canonical output |
| Link nodes `LN:` | Pending | P3 | CXSMILES | Variable attachments |
| Data S-groups `SgD:` | Pending | P3 | CXSMILES | Embedded data |
| Hydrogen bonding `H:` | Skip | P3 | CXSMILES | Query feature |
| Wedge bond info `w:` | Pending | P3 | CXSMILES | Stereo display |
| Multicenter bonds `m:` | Pending | P3 | CXSMILES | Organometallics |

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
2. CXSMILES extension block parsing `|...|`
3. CXSMILES coordinates `c:x,y,z;...`
4. Reaction SMILES `>>` splitting
5. Atom mapping in reactions `[C:1]`

**P2 (Medium Priority) - Common advanced features:**

1. Aromatic `se`, `te`, `b`
2. Dative bonds `->` `<-`
3. Radicals (OpenBabel dot notation and CXSMILES `^n:`)
4. CXSMILES atom labels `$name$`
5. CXSMILES enhanced stereo groups
6. CXSMILES S-groups (basic: abbreviations, polymers)
7. SMARTS query primitives
8. R-group notation `[R1]`

**P3 (Low Priority) - Rare or niche:**

1. Quadruple bonds `$`
2. Extended ring closures `%(nnn)`
3. CXSMILES link nodes, data S-groups
4. DeepSMILES (as output format only)
5. Multicenter bonds

**Skip (Not implementing):**

1. SELFIES (misleading validity claims)
2. Jmol animation directives
3. CXSMILES query-only features (`rb:`, `s:`, `u:`)
4. Daylight vector binding `&&`
5. OpenBabel external bonds `&` (no other toolkit supports)
6. Any bond `~` in SMILES (SMARTS only; CXSMILES `|Z:|` is different)

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
