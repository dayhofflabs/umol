# MOL Parsing IR Implementation Plan

* Implementation Plan
* Date: 2025-08-24

## Overview

This document outlines the implementation plan for refactoring the MOL parser to use an intermediate representation (IR) with a generic `ParseTarget` trait system. This refactoring addresses several issues:

1. **Unicode whitespace handling** - Immediate issue with non-breaking spaces in MOL files
2. **Type naming clarity** - "Standard" vs "General" is confusing terminology  
3. **Parser duplication** - Separate standard/general parsers with duplicated logic
4. **Configuration flexibility** - Need fine-grained control over parsing behavior
5. **Future IR compatibility** - Step toward full IR architecture from `28-intermediate-representation-2025-08-18.md`

## Goals

- **Single unified parser** with generic return types
- **Clear type naming** that matches chemist expectations
- **Configurable parsing behavior** for different use cases
- **Intermediate representation** as stepping stone to full IR architecture
- **Immediate Unicode whitespace fix** for compliance test issues
- **No backward compatibility constraints** - clean break acceptable

## Type System Redesign

### Current → New Type Names

| Current Type | New Type | Rationale |
|--------------|----------|-----------|
| `AtomStandard` | `Atom` | What chemists expect "Atom" to mean |
| `Atom` | `AtomLike` | Can represent queries, wildcards, etc. |
| `BondStandard` | `Bond` | Concrete bond with defined properties |
| `Bond` | `BondLike` | Can represent query bonds, variable order |
| `MoleculeStandard` | `Molecule` | Single, concrete molecule |
| `Molecule` | `MoleculeLike` | Can represent queries, polymers, etc. |

### Type Characteristics

```rust
// Concrete types - represent real chemical entities
pub struct Atom {
    pub element: Element,
    pub position: (f64, f64, f64),
    pub formal_charge: i8,
    pub isotope: Option<u16>,
    // Only concrete properties, no queries
}

pub struct Bond {
    pub order: BondOrder,  // Single, Double, Triple, Aromatic
    pub stereo: Option<BondStereo>,
    // No query bond types
}

pub struct Molecule {
    pub graph: petgraph::Graph<Atom, Bond>,
    // No S-groups, R-groups, or query features
}

// General types - can represent queries and complex structures  
pub struct AtomLike {
    pub element_or_query: ElementOrQuery,  // Element, Any, AtomList, etc.
    pub position: Option<(f64, f64, f64)>,
    pub formal_charge: i8,
    pub query_properties: Option<AtomQueryProperties>,
    // All possible atom features
}

pub struct BondLike {
    pub order: BondOrderOrQuery,  // Including query bond types
    pub stereo: Option<BondStereo>,
    pub query_properties: Option<BondQueryProperties>,
    // All possible bond features
}

pub struct MoleculeLike {
    pub graph: petgraph::Graph<AtomLike, BondLike>,
    pub sgroups: Vec<SGroup>,
    pub rgroups: Vec<RGroup>,
    // All possible molecular features
}
```

## Intermediate Representation (IR)

### RawMolecule Structure

Located in `io::ir::molecule.rs`:

```rust
/// Raw intermediate representation from MOL file parsing
/// Contains all information from the file without validation or interpretation
pub struct RawMolecule {
    pub header: Header,
    pub atoms: Vec<RawAtom>,
    pub bonds: Vec<RawBond>, 
    pub properties: Vec<PropertyEntries>,
    pub sgroups: Vec<SGroup>,
}

/// Union of all possible atom data from any supported format
pub struct RawAtom {
    // Core properties (always present)
    pub element: Element,
    pub position: Option<(f64, f64, f64)>,
    pub formal_charge: i8,
    pub isotope: Option<u16>,
    
    // Query-specific properties (MOL queries, future SMARTS)
    pub query_type: Option<AtomQueryType>,
    pub atom_list: Option<Vec<Element>>,
    pub attachment_point: Option<u8>,
    pub ring_bond_count: Option<u8>,
    pub substitution_count: Option<u8>,
    pub unsaturated: Option<bool>,
    
    // Format-specific metadata
    pub source_format: SourceFormat,
    pub original_text: Option<String>,  // For debugging
}

pub struct RawBond {
    pub atom_indices: (usize, usize),
    pub order: BondOrderOrQuery,
    pub stereo: Option<BondStereo>,
    pub topology: Option<BondTopology>,
    pub reacting_center: Option<ReactingCenter>,
    
    // Format-specific metadata  
    pub source_format: SourceFormat,
    pub original_text: Option<String>,
}

pub enum SourceFormat {
    MOL,
    SMILES,    // For future use
    SMARTS,    // For future use
}

impl RawMolecule {
    pub fn has_query_features(&self) -> bool {
        self.atoms.iter().any(|a| a.has_query_features()) ||
        self.bonds.iter().any(|b| b.has_query_features()) ||
        !self.sgroups.is_empty()
    }
    
    pub fn has_3d_coordinates(&self) -> bool {
        self.atoms.iter().any(|a| a.position.is_some())
    }
}
```

### ParseTarget Trait

Located in `io::ir::traits.rs`:

```rust
/// Trait for types that can be constructed from parsed molecular data
pub trait ParseTarget: Sized {
    /// Whether this target type accepts query features (atom lists, wildcards, etc.)
    fn allows_query_features() -> bool;
    
    /// Whether this target type accepts S-groups
    fn allows_sgroups() -> bool;
    
    /// Whether this target type accepts R-groups  
    fn allows_rgroups() -> bool;
    
    /// Convert parsed data to target type with validation
    fn from_parsed_data(
        parsed: RawMolecule,
        config: &MolParsingConfig,
    ) -> Result<Self, ParseError>;
}

impl ParseTarget for Molecule {
    fn allows_query_features() -> bool { false }
    fn allows_sgroups() -> bool { false }
    fn allows_rgroups() -> bool { false }
    
    fn from_parsed_data(
        parsed: RawMolecule,
        config: &MolParsingConfig,
    ) -> Result<Self, ParseError> {
        // Strict validation
        if parsed.has_query_features() && !config.force_allow_queries {
            return Err(ParseError::QueryFeaturesNotAllowed {
                target_type: "Molecule",
                found_features: parsed.describe_query_features(),
            });
        }
        
        if !parsed.sgroups.is_empty() && !config.force_allow_sgroups {
            return Err(ParseError::SGroupsNotAllowed {
                target_type: "Molecule",
                sgroup_count: parsed.sgroups.len(),
            });
        }
        
        // Convert to concrete molecule
        Ok(Self {
            graph: build_concrete_graph(&parsed.atoms, &parsed.bonds)?,
            properties: convert_concrete_properties(&parsed.properties)?,
        })
    }
}

impl ParseTarget for MoleculeLike {
    fn allows_query_features() -> bool { true }
    fn allows_sgroups() -> bool { true }
    fn allows_rgroups() -> bool { true }
    
    fn from_parsed_data(
        parsed: RawMolecule,
        _config: &MolParsingConfig,
    ) -> Result<Self, ParseError> {
        // Accept everything - no validation restrictions
        Ok(Self {
            graph: build_general_graph(&parsed.atoms, &parsed.bonds)?,
            properties: convert_all_properties(&parsed.properties)?,
            sgroups: parsed.sgroups,
        })
    }
}
```

## Configuration System

### MolParsingConfig Structure

Located in `io::config.rs`:

```rust
#[derive(Debug, Clone)]
pub struct MolParsingConfig {
    // === Immediate Unicode Whitespace Fix ===
    /// Accept Unicode whitespace characters (U+00A0, U+2000-U+200A, etc.) as equivalent to ASCII spaces
    pub accept_unicode_whitespace: bool,
    
    // === Field Parsing Behavior ===
    /// Allow truncated lines (current behavior)
    pub allow_truncated_lines: bool,
    /// Allow missing optional fields at end of lines
    pub allow_missing_trailing_fields: bool,
    /// Strict validation of numeric field formats
    pub strict_numeric_fields: bool,
    /// Normalize line endings (\r\n, \r to \n)
    pub normalize_line_endings: bool,
    
    // === Format Compliance ===
    /// Enforce 80-character line length limit
    pub enforce_v2000_line_length: bool,
    /// Validate proper element symbol capitalization
    pub validate_atom_symbol_case: bool,
    /// Strict counts line format validation
    pub require_counts_line_format: bool,
    
    // === Target Type Override Controls ===
    /// Allow query features even when target type normally rejects them
    pub force_allow_queries: bool,
    /// Allow S-groups even when target type normally rejects them  
    pub force_allow_sgroups: bool,
    /// Allow R-groups even when target type normally rejects them
    pub force_allow_rgroups: bool,
    
    // === Error Handling & Diagnostics ===
    /// Collect non-fatal parsing warnings
    pub collect_warnings: bool,
    /// Preserve original text for debugging
    pub preserve_original_text: bool,
    /// Include detailed line/column info in errors
    pub detailed_error_positions: bool,
    /// Stop on first error vs. collect all errors
    pub fail_on_first_error: bool,
    
    // === Performance Optimization ===
    /// Skip query-specific parsing when target type doesn't allow queries
    pub optimize_for_concrete_types: bool,
    /// Fail fast if query features detected for strict target types
    pub early_query_detection: bool,
}

impl MolParsingConfig {
    /// Configuration optimized for real-world files with common quirks
    pub fn lenient() -> Self {
        Self {
            accept_unicode_whitespace: true,
            allow_truncated_lines: true,
            allow_missing_trailing_fields: true,
            strict_numeric_fields: false,
            normalize_line_endings: true,
            enforce_v2000_line_length: false,
            validate_atom_symbol_case: false,
            require_counts_line_format: false,
            force_allow_queries: false,
            force_allow_sgroups: false, 
            force_allow_rgroups: false,
            collect_warnings: true,
            preserve_original_text: true,
            detailed_error_positions: true,
            fail_on_first_error: false,
            optimize_for_concrete_types: true,
            early_query_detection: false,
        }
    }
    
    /// Configuration for strict format compliance testing
    pub fn strict() -> Self {
        Self {
            accept_unicode_whitespace: false,
            allow_truncated_lines: false,
            strict_numeric_fields: true,
            enforce_v2000_line_length: true,
            validate_atom_symbol_case: true,
            require_counts_line_format: true,
            fail_on_first_error: true,
            early_query_detection: true,
            ..Self::lenient()
        }
    }
    
    /// Configuration automatically tuned for target type
    pub fn for_target<T: ParseTarget>() -> Self {
        Self {
            optimize_for_concrete_types: !T::allows_query_features(),
            early_query_detection: !T::allows_query_features(),
            ..Self::default()
        }
    }
}

impl Default for MolParsingConfig {
    fn default() -> Self {
        Self {
            // Solve immediate Unicode whitespace issue
            accept_unicode_whitespace: true,
            // Keep current behavior for compatibility
            allow_truncated_lines: true,
            // Reasonable defaults balancing usability and correctness
            strict_numeric_fields: false,
            collect_warnings: true,
            detailed_error_positions: true,
            fail_on_first_error: false,
            // Default other fields...
            normalize_line_endings: true,
            enforce_v2000_line_length: false,
            allow_missing_trailing_fields: true,
            validate_atom_symbol_case: false,
            require_counts_line_format: false,
            force_allow_queries: false,
            force_allow_sgroups: false,
            force_allow_rgroups: false,
            preserve_original_text: false,
            optimize_for_concrete_types: false,
            early_query_detection: false,
        }
    }
}
```

## Generic Parser Implementation

### Core Parser Function

Located in `io::mol::parser.rs`:

```rust
/// Generic MOL parser that returns any type implementing ParseTarget
pub fn parse_mol<T: ParseTarget>(
    input: &[u8],
    config: &MolParsingConfig,
) -> Result<T, ParseError> {
    // Phase 1: Parse to intermediate representation
    let parsed = parse_to_ir(input, config)?;
    
    // Phase 2: Early validation based on target type capabilities
    if !T::allows_query_features() && parsed.has_query_features() && config.early_query_detection {
        return Err(ParseError::QueryFeaturesNotAllowed {
            target_type: std::any::type_name::<T>(),
            found_features: parsed.describe_query_features(),
        });
    }
    
    // Phase 3: Convert to target type with validation
    T::from_parsed_data(parsed, config)
}

/// Parse MOL file including header information
pub fn parse_mol_file<T: ParseTarget>(
    input: &[u8], 
    config: &MolParsingConfig
) -> Result<MolFile<T>, ParseError> {
    let molecule = parse_mol::<T>(input, config)?;
    let header = extract_header_from_parsed(&molecule, input)?;
    Ok(MolFile::new(header, molecule))
}

/// Internal function: Parse MOL text to intermediate representation
fn parse_to_ir(
    input: &[u8],
    config: &MolParsingConfig,
) -> Result<RawMolecule, ParseError> {
    let mut warnings = if config.collect_warnings { Some(Vec::new()) } else { None };
    
    // Preprocess input based on config
    let processed_input = preprocess_input(input, config)?;
    
    // Parse header (3 lines)
    let (remaining, header) = parse_header(&processed_input)?;
    
    // Parse counts line
    let (remaining, counts) = parse_counts_line(remaining, config)?;
    
    // Parse atoms
    let (remaining, atoms) = parse_atoms_to_ir(
        remaining, 
        counts.atoms() as usize, 
        config,
        &mut warnings
    )?;
    
    // Parse bonds  
    let (remaining, bonds) = parse_bonds_to_ir(
        remaining,
        counts.bonds() as usize,
        config,
        &mut warnings
    )?;
    
    // Parse properties
    let (remaining, properties, sgroups) = parse_properties_to_ir(
        remaining,
        config,
        &mut warnings
    )?;
    
    Ok(RawMolecule {
        header,
        atoms,
        bonds,
        properties,
        sgroups,
    })
}

/// Preprocess input according to configuration
fn preprocess_input(input: &[u8], config: &MolParsingConfig) -> Result<Vec<u8>, ParseError> {
    let mut processed = input.to_vec();
    
    if config.normalize_line_endings {
        // Convert \r\n and \r to \n
        processed = processed.replace(b"\r\n", b"\n").replace(b"\r", b"\n");
    }
    
    if config.accept_unicode_whitespace {
        // Replace Unicode whitespace with ASCII spaces
        processed = replace_unicode_whitespace(processed)?;
    }
    
    Ok(processed)
}

/// Replace Unicode whitespace characters with ASCII spaces
fn replace_unicode_whitespace(input: Vec<u8>) -> Result<Vec<u8>, ParseError> {
    let text = String::from_utf8(input)
        .map_err(|e| ParseError::InvalidUtf8(e))?;
    
    // Replace common Unicode whitespace with ASCII space
    let normalized = text
        .replace('\u{00A0}', " ")  // Non-breaking space
        .replace('\u{2000}', " ")  // En quad  
        .replace('\u{2001}', " ")  // Em quad
        .replace('\u{2002}', " ")  // En space
        .replace('\u{2003}', " ")  // Em space
        .replace('\u{2004}', " ")  // Three-per-em space
        .replace('\u{2005}', " ")  // Four-per-em space
        .replace('\u{2006}', " ")  // Six-per-em space
        .replace('\u{2007}', " ")  // Figure space
        .replace('\u{2008}', " ")  // Punctuation space
        .replace('\u{2009}', " ")  // Thin space
        .replace('\u{200A}', " ")  // Hair space
        .replace('\u{3000}', " "); // Ideographic space
    
    Ok(normalized.into_bytes())
}
```

### Convenience APIs

```rust
// Type-specific convenience constructors
impl Molecule {
    pub fn from_mol_str(s: &str) -> Result<Self, ParseError> {
        parse_mol(s.as_bytes(), &MolParsingConfig::for_target::<Self>())
    }
    
    pub fn from_mol_bytes(bytes: &[u8]) -> Result<Self, ParseError> {
        parse_mol(bytes, &MolParsingConfig::for_target::<Self>())
    }
    
    pub fn from_mol_bytes_with_config(
        bytes: &[u8], 
        config: &MolParsingConfig
    ) -> Result<Self, ParseError> {
        parse_mol(bytes, config)
    }
}

impl MoleculeLike {
    pub fn from_mol_str(s: &str) -> Result<Self, ParseError> {
        parse_mol(s.as_bytes(), &MolParsingConfig::for_target::<Self>())
    }
    
    pub fn from_mol_bytes(bytes: &[u8]) -> Result<Self, ParseError> {
        parse_mol(bytes, &MolParsingConfig::for_target::<Self>())
    }
    
    pub fn from_mol_bytes_with_config(
        bytes: &[u8], 
        config: &MolParsingConfig
    ) -> Result<Self, ParseError> {
        parse_mol(bytes, config)
    }
}
```

## Module Structure

```
umol-models-graph/src/io/
├── config.rs                 # MolParsingConfig
├── ir.rs
├── ir/                       # Intermediate representation
│   ├── molecule.rs           # RawMolecule, RawAtom, RawBond  
│   ├── traits.rs             # ParseTarget trait
│   └── convert.rs            # Conversion utilities
├── mol.rs
├── mol/
│   ├── parser.rs             # Main generic parser functions
│   ├── parser/
│   │   ├── mod.rs
│   │   ├── header.rs         # Header parsing (unchanged)
│   │   ├── preprocess.rs     # Input preprocessing
│   │   └── tests.rs          # Updated tests
├── ctab.rs
├── ctab/                     # CTAB parsing logic (mostly unchanged)
│   ├── parser.rs             # Updated to produce RawAtom/RawBond
│   └── ...
└── ...
```

## Implementation Phases

### Phase 1: Foundation (Week 1)
1. **Create module structure** - Set up `io::ir` and `io::config` modules
2. **Define new types** - Rename existing types (breaking change)
3. **Implement RawMolecule** - Basic IR structure
4. **Create MolParsingConfig** - Configuration system
5. **Update imports** throughout codebase

### Phase 2: ParseTarget Trait (Week 1-2)
1. **Define ParseTarget trait** in `io::ir::traits`
2. **Implement for Molecule** with strict validation
3. **Implement for MoleculeLike** with permissive validation
4. **Create conversion utilities** in `io::ir::convert`

### Phase 3: Generic Parser Core (Week 2)
1. **Implement parse_to_ir()** - Convert existing parsing to produce IR
2. **Update CTAB parsers** to produce RawAtom/RawBond
3. **Implement Unicode whitespace handling** - Immediate issue fix
4. **Add input preprocessing** - Line ending normalization, etc.

### Phase 4: Generic API (Week 2-3)
1. **Implement parse_mol<T>()** - Main generic function
2. **Add parse_mol_file<T>()** - File variant
3. **Create convenience methods** - Type-specific constructors
4. **Add error types** - Detailed error messages

### Phase 5: Testing & Migration (Week 3-4)
1. **Update existing tests** - Convert to new API
2. **Add configuration tests** - Test all config options
3. **Add Unicode whitespace tests** - Test the immediate fix
4. **Performance testing** - Ensure optimization flags work
5. **Create migration examples** - Document API changes

### Phase 6: Documentation & Cleanup (Week 4)
1. **Update API documentation** - Document new generic approach
2. **Create migration guide** - Help users transition
3. **Add examples** - Show common usage patterns
4. **Remove deprecated code** - Clean up old APIs

## Migration Strategy

### Immediate Breaking Changes
Since backward compatibility is not required, we can make clean breaks:

```rust
// Old API (remove immediately)
pub fn parse_mol(input: &[u8]) -> Result<Molecule>           // Remove
pub fn parse_mol_standard(input: &[u8]) -> Result<MoleculeStandard> // Remove

// New API
pub fn parse_mol<T: ParseTarget>(input: &[u8], config: &MolParsingConfig) -> Result<T>
```

### Type Migration
```rust
// Users need to update their code:
// Old:
let mol: MoleculeStandard = parse_mol_standard(input)?;

// New:
let mol: Molecule = parse_mol(input, &MolParsingConfig::default())?;
// Or more concisely:
let mol: Molecule = Molecule::from_mol_bytes(input)?;
```

### Gradual Feature Addition
- Start with basic configuration options
- Add more sophisticated options as needed
- Unicode whitespace fix is highest priority

## Testing Strategy

### Test Categories

1. **Type System Tests**
   - Verify ParseTarget implementations work correctly
   - Test conversion from RawMolecule to concrete types
   - Validate error handling for incompatible features

2. **Configuration Tests**  
   - Test each configuration option individually
   - Test configuration combinations
   - Verify performance optimization flags

3. **Unicode Whitespace Tests**
   - Test various Unicode whitespace characters
   - Verify correct parsing of problematic files like `hydrogen_isotopes.mol`
   - Test configuration toggle works correctly

4. **Parsing Correctness Tests**
   - Convert existing compliance tests to new API
   - Ensure identical results between old and new parsers
   - Test edge cases and error conditions

5. **Performance Tests**
   - Benchmark generic vs. specialized paths
   - Verify optimization flags improve performance
   - Test memory usage with large files

### Test Files

Update existing test structure:
```
tests/
├── compliance.rs             # Updated for new API
├── mol_parsing/
│   ├── mod.rs                # Updated test infrastructure  
│   ├── unicode_tests.rs      # New: Unicode whitespace tests
│   ├── config_tests.rs       # New: Configuration tests
│   └── data/                 # Existing test files
└── ...
```

## Error Handling

### ParseError Enum

```rust
#[derive(Debug, thiserror::Error)]
pub enum ParseError {
    #[error("Query features not allowed in target type {target_type}: {found_features:?}")]
    QueryFeaturesNotAllowed {
        target_type: &'static str,
        found_features: Vec<String>,
    },
    
    #[error("S-groups not allowed in target type {target_type}: found {sgroup_count} S-groups")]
    SGroupsNotAllowed {
        target_type: &'static str,
        sgroup_count: usize,
    },
    
    #[error("Invalid UTF-8 in MOL file: {0}")]
    InvalidUtf8(std::string::FromUtf8Error),
    
    #[error("Invalid numeric field at line {line}, column {column}: {value:?}")]
    InvalidNumericField {
        line: usize,
        column: usize,
        value: String,
    },
    
    #[error("Format compliance error: {message}")]
    FormatCompliance { message: String },
    
    // ... other error variants
}
```

## Success Metrics

1. **Immediate Issue Resolution** - Unicode whitespace in `hydrogen_isotopes.mol` parses correctly
2. **API Clarity** - Users can intuitively choose between `Molecule` and `MoleculeLike`
3. **Performance** - No regression in parsing speed for optimized cases
4. **Flexibility** - Configuration system handles diverse real-world files
5. **Future Compatibility** - Clear path to full IR architecture implementation

## Future Integration with Full IR

This implementation provides a natural stepping stone to the full IR architecture described in `28-intermediate-representation-2025-08-18.md`:

- `RawMolecule` becomes the "Raw IR"
- `ParseTarget::from_parsed_data()` becomes the "Validation" phase  
- Final types (`Molecule`, `MoleculeLike`) become the "Validated Molecule" types
- SMILES parser can target same `RawMolecule` IR
- Validation logic can be extracted into separate validator components

The modular design ensures this refactoring provides immediate value while aligning with long-term architectural goals.
