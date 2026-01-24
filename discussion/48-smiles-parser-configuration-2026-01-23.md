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

---

## 1. Current State (to be removed)

```rust
SmilesParseFlags {
    EXTENDED_WS,        // speculative - remove
    ALLOWS_COMMENTS,    // speculative - remove
    EXPLICIT_EOI,       // speculative - remove
    NO_METADATA,        // keep? review
    STRICT_OPENSMILES,  // keep as baseline
    UMOL_DIALECT,       // remove
    LENIENT,            // remove
}
```

**Action:** Strip back to minimal baseline parser without these extensions.

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

SMARTS extends SMILES with query features:

**Atomic primitives:**
- `*` - any atom
- `a` - aromatic
- `A` - aliphatic
- `D<n>` - degree
- `H<n>` - H count
- `h<n>` - implicit H count
- `R<n>` - ring membership
- `r<n>` - ring size
- `v<n>` - valence
- `X<n>` - connectivity
- `x<n>` - ring connectivity
- `#n` - atomic number

**Logical operators:**
- `&` - and (high precedence)
- `,` - or
- `;` - and (low precedence)
- `!` - not

**Bonds:**
- `~` - any bond
- `@` - ring bond

**Recursive SMARTS:**
- `$()` - recursive SMARTS pattern

**Design question:** Separate parser or shared infrastructure with SMILES?

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

### Group 6: Reactions
- Reaction arrows `>>`
- Agent notation `>`
- Atom mapping for reactions

### Group 7: SMARTS (separate parser?)
- Query primitives
- Logical operators
- Recursive patterns

---

## 7. Implementation Plan

### Phase 1: Cleanup
- [ ] Remove UMOL_DIALECT, ALLOWS_COMMENTS, EXTENDED_WS
- [ ] Remove corresponding parser code paths
- [ ] Establish clean OpenSMILES baseline

### Phase 2: Feature Audit
- [ ] Analyze 473 failing SMILES from conformance suite
- [ ] Categorize by failure reason
- [ ] Map to feature groups above

### Phase 3: Extended Atoms (Group 2)
- [ ] Add aromatic `se`, `te`, `b` support
- [ ] Add wildcard `*` support
- [ ] Feature gate: `EXTENDED_AROMATICS`

### Phase 4: Extended Stereo (Group 3)
- [ ] Implement @AL, @SP, @TB, @OH parsing
- [ ] Relaxed stereo validation mode
- [ ] Feature gate: `EXTENDED_STEREO`

### Phase 5: CXSMILES (Group 5)
- [ ] Parse `|...|` extension block
- [ ] Extract coordinates, labels, radicals
- [ ] Feature gate: `CXSMILES`

### Phase 6: Presets
- [ ] `opensmiles_strict` - Group 1 only
- [ ] `rdkit_compat` - Groups 1-4
- [ ] `cxsmiles` - Groups 1-5

### Phase 7: SMARTS
- [ ] Decide: separate parser or shared tokenizer
- [ ] Implement SMARTS primitives
- [ ] Implement logical operators

---

## 8. Open Questions

1. Should SMARTS use the same lexer as SMILES with different parser rules?
2. How to handle toolkit-specific quirks (e.g., RDKit sanitization bypass)?
3. Should we support SMIRKS (reaction transforms) as part of SMARTS?
4. How granular should feature gates be? Individual flags vs grouped presets?

---

## 9. References

- OpenSMILES spec: http://opensmiles.org/
- Daylight theory manual (archived)
- ChemAxon CXSMILES docs
- RDKit SMILES parsing docs
- Depth-First blog (ChemCore comparisons)
