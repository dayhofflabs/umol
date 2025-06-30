# Geometric Molecular Models Design Discussion

## Overview

Discussion about implementing a molecular model that considers only nuclear (atomic) positions as points in 3D space without bond definitions or explicit electronic structure. The goal is to design clean semantic models and APIs for different spatial representations.

## Key Decisions

### 1. Crate Naming
- **Chosen name:** `umol-models-geometric` 
- **Rationale:** Emphasizes purely geometric nature, distinguishes from graph-based models
- **Rejected alternatives:** 
  - `umol-models-bo` (Born-Oppenheimer implies electronic structure)
  - `umol-models-classical` (too vague)
  - `umol-models-spatial` (less domain-specific)

### 2. Semantic Model Architecture

**Two distinct semantic model categories identified:**
1. **Core Geometric Models:** Different mathematical representations of atomic positions
2. **I/O Format Models:** Models that capture exactly what each file format can represent

**Key insight:** Different coordinate systems (Cartesian, internal, fractional) are distinct representations, not just views of the same data, due to:
- Different degrees of freedom (Cartesian has 6 degenerate SE(3) DOF)
- Non-linear, iterative transformations between them
- Different optimization landscapes and numerical properties

### 3. Module Structure

```
umol-models-geometric/
├── cartesian/           # Cartesian coordinate representation
│   ├── Molecule, Atom
├── valence/             # Internal/valence coordinate representation  
│   ├── Molecule, Atom, Coordinate
├── redundant/           # Redundant internal coordinates (future)
│   ├── Molecule, Atom, Coordinate
├── io/                  # Format-specific semantic models
│   ├── xyz/
│   │   ├── Molecule, Atom
│   ├── zmatrix/
│   │   ├── Molecule, Atom, Distance, Angle, Dihedral
│   └── pdb/             # (future - coordinates only)
└── conversions/         # Inter-representation conversions
```

**Naming Philosophy:**
- Short, simple names over verbose specificity
- `cartesian::Molecule` not `CartesianMolecule`
- Emphasizes these are different *representations* of molecules, not different *types*

### 4. API Design Principles

**Core Model APIs:**

**Cartesian Model:**
- Pure positional data, no connectivity except when needed for internal coordinate definition
- SE(3) operations (translation, rotation, alignment)
- Geometric queries (distances, angles, dihedrals)

**Valence (Internal) Model:**
- Tree-structured internal coordinates
- Forward kinematics to Cartesian
- Optimization-friendly gradient transformations

**I/O Format Models:**
- Exact representation of what each format can store
- No hidden assumptions or implicit conversions
- Explicit conversion to/from geometric models

### 5. Conversion Strategy

**Multiple Internal Coordinate Problem:**
Any `cartesian::Molecule` can have many valid `valence::Molecule` representations due to:
- Connectivity choice
- Reference frame choice  
- Coordinate tree structure
- Coordinate ordering

**Solution: Hybrid Approach (Default + Builder)**
```rust
impl cartesian::Molecule {
    // Convenient default for most users
    fn to_valence(&self) -> Result<valence::Molecule>
    
    // Builder pattern for explicit control
    fn valence_builder(&self) -> ValenceCoordinateBuilder
}
```

**Error Handling:**
- All conversions return `Result` types from umol crate
- Rich error types for common pathologies (singular coordinates, non-convergence, etc.)
- Incorporate quantum chemistry best practices for coordinate choice

### 6. Conversion Complexities

**Forward (Cartesian → Valence):**
- Non-unique: multiple valid representations possible
- Requires connectivity determination
- Can fail at singularities (linear molecules, overlapping atoms)

**Reverse (Valence → Cartesian):**
- Forward kinematics, generally more reliable
- Can fail with overcomplete/redundant coordinates
- May require iterative solvers with convergence issues
- Need different strategies for different molecular types

**I/O Conversions:**
- XYZ → Cartesian: straightforward parsing
- Z-matrix → Valence: direct mapping
- Complex formats may compose through multiple models

## Open Questions

1. **Valence naming conflict:** Potential confusion with `umol-models-valence` crate (valence graph model). May need to revisit naming.

2. **Default strategies:** Need to incorporate quantum chemistry best practices for:
   - Connectivity determination (distance-based with element-specific cutoffs)
   - Root selection strategies
   - Tree building algorithms that avoid singularities
   - Coordinate ordering preferences

3. **Performance:** Caching strategies for expensive conversions, especially for repeated Cartesian ↔ Valence transformations.

4. **Redundant coordinates:** Future support for overcomplete internal coordinate sets, which are more robust for optimization but create their own conversion challenges.

## Implementation Priority

1. Basic `cartesian::Molecule` with geometric operations
2. `io::xyz::Molecule` with conversion to/from Cartesian  
3. `valence::Molecule` with default conversion strategy
4. Builder pattern for controlled valence coordinate generation
5. Rich error types for common failure modes
6. `io::zmatrix::Molecule` support

## Notes

- Heavy influence from quantum chemistry best practices
- Emphasis on explicit semantics over hidden assumptions  
- Performance important but semantics come first
- Learn from OpenBabel's over-abstraction mistakes 