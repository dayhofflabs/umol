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
- [ ] Extended aromatic atoms: `te`, `se` (tellurium, selenium)
- [ ] Aromatic boron: `b`
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

### Group 2: Extended Atoms
- Aromatic `se`, `te`, `b`
- Wildcard `*`

### Group 3: Extended Stereo
- Allene (@AL), square planar (@SP)
- Trigonal bipyramidal (@TB), octahedral (@OH)
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

### Phase 3: Extended SMILES Parser + Wildcards
- [ ] Create `ExtendedMolecule` type (or verify existing structure suffices)
- [ ] Implement `parse_extended_smiles_bytes_with` and wrappers
- [ ] Add wildcard `*` support (standalone and `[*:n]` with atom class)
- [ ] Ensure `From<Molecule> for ExtendedMolecule`
- [ ] Update conformance suite to test extended parser

**Impact:** Fixes 32+ conformance failures (wildcard cases)

### Phase 3b: Extended Aromatic Atoms
- [ ] Add aromatic `se`, `te` (selenium, tellurium) support
- [ ] Add aromatic `b` (boron) support
- [ ] Feature gate: `EXTENDED_AROMATICS` (applies to both parsers)

### Phase 4: Extended Stereo (Group 3)
- [ ] Implement @AL, @SP, @TB, @OH parsing
- [ ] Relaxed stereo validation mode
- [ ] Feature gate: `EXTENDED_STEREO`

### Phase 5: Extended Bonds (Group 4)
- [ ] Dative bonds (`->`, `<-`)
- [ ] Quadruple bonds (`$`)
- [ ] Feature gate: `EXTENDED_BONDS`

### Phase 6: CXSMILES (Group 5)
- [ ] Parse `|...|` extension block
- [ ] Extract coordinates, labels, radicals
- [ ] Enhanced stereo groups
- [ ] Feature gate: `CXSMILES`

### Phase 7: Presets
- [ ] `opensmiles_strict` - Group 1 only (current)
- [ ] `opensmiles_extended` - Groups 1-4
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

## 11. References

- OpenSMILES spec: http://opensmiles.org/
- Daylight theory manual (archived)
- ChemAxon CXSMILES docs
- RDKit SMILES parsing docs
- Depth-First blog (ChemCore comparisons)
