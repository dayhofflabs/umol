# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Development Commands

### Build Commands
- `cargo build` - Build all workspace crates
- `cargo build --package umol-models-graph` - Build specific package
- `cargo build --release` - Build optimized release version

### Testing Commands
- `cargo test` - Run all tests across the workspace
- `cargo test --package umol` - Test specific package
- `cargo test --package umol-models-graph` - Test the graph models (includes MOL parsing tests)
- `cargo test <test_name>` - Run specific test by name
- `cargo test -- --nocapture` - Run tests with output visible
- `cargo test --lib` - Test only library code (no integration tests)

### Linting and Formatting
- `cargo clippy` - Run Clippy linter
- `cargo clippy -- -D warnings` - Run Clippy treating warnings as errors
- `cargo fmt` - Format code using rustfmt (configured in rustfmt.toml)
- `cargo check` - Fast compilation check without producing binaries

### Specialized Tools
- `cargo run --bin mol_classifier -- <mol_file>` - Classify MOL files as molecule/moleculelike/invalid
- `cargo run --bin test_mol_file -- <mol_file>` - Test individual MOL files with both parsers
- `cargo bench` - Run benchmarks (parsing performance tests)
- `typos` - Check for typos (configured in typos.toml, excludes chemical element symbols)

### Documentation
- `cargo doc --open` - Generate and open documentation
- `cargo doc --no-deps` - Generate docs only for workspace crates

## Project Architecture

### Workspace Structure
This is a multi-crate Rust workspace for cheminformatics with specialized responsibilities:

- **`umol`** - Core library defining fundamental traits and abstractions:
  - Entity/Relation system for semantic chemical objects
  - Model trait with capability system for molecular representations
  - Property computation framework
  - Conversion and operation abstractions

- **`umol-data`** - Chemical element and isotope data:
  - Periodic table information from authoritative sources
  - NUBASE2020 isotope data
  - Atomic mass evaluation data

- **`umol-macros`** - Procedural macros for code generation

- **`umol-models`** - Basic molecular models and stoichiometry

- **`umol-models-valence`** - Valence-based molecular modeling:
  - Atom and bond specifications with validation
  - Matcher systems for structural queries
  - Registry-based atom/bond type management

- **`umol-models-graph`** - Graph-based molecular representations:
  - Comprehensive MOL file parser using nom combinators
  - Two-tier parsing: basic molecules and extended "moleculelike" structures
  - SMILES parsing foundation with Logos lexer and LALRPOP grammar
  - Intermediate representation (IR) for format-agnostic handling

- **`umol-models-geometric`** - Geometric molecular models and spatial representations

### Key Design Patterns

#### Entity-Capability Architecture
The codebase uses an entity-capability pattern where:
- `Entity` trait represents semantic objects (molecules, conformers, reactions) with namespaced IDs
- `Model` trait represents computational molecular models with declared capabilities
- `Capability` system allows querying what operations models support
- `Property` trait enables polymorphic computation over different model types

#### Parsing Strategy
- **MOL files**: nom combinator-based parsing for fixed-width format compliance
- **SMILES**: Generated parser using LALRPOP + Logos for recursive grammar
- **Intermediate Representation**: Format-agnostic IR for unified handling across file types
- **Two-tier validation**: Basic structural parsing + extended feature parsing

#### Performance Optimizations
- Zero-copy parsing with byte slice operations
- SmallVec for typical molecular sizes
- Build-time code generation (LALRPOP, custom build scripts)
- Compile-time capability checking where possible

### Testing Philosophy

#### Automated Test Generation
- Build scripts dynamically generate test functions based on directory structure
- MOL files automatically classified as molecule/moleculelike/invalid
- Snapshot testing with `insta` for regression prevention

#### Compliance Testing
The project emphasizes chemical format compliance:
- Real-world MOL file testing against specification
- OpenSMILES specification adherence for SMILES parsing
- Extensive edge case coverage for chemical file format quirks

### Toolchain Requirements
- Uses nightly Rust (pinned to `nightly-2025-08-27` in rust-toolchain.toml)
- Required for advanced procedural macro features and parsing optimizations
- LALRPOP for parser generation
- Build scripts for dynamic test generation

### Chemical Domain Knowledge
This codebase handles sophisticated cheminformatics concepts:
- Molecular graph representations with aromaticity
- Stereochemistry and chirality validation
- Query atoms and substructure matching
- S-groups and R-groups for extended molecular features
- Ring closure tracking for SMILES parsing
- Valence state validation and chemical correctness

## Development Context

### File Format Expertise
When working with chemical file formats, understand that:
- MOL format has fixed-width fields requiring exact spacing
- Unicode whitespace preprocessing is critical for non-ASCII spaces
- Coordinate precision varies (handles `.0000` without leading zeros)
- Extended features (S-groups, R-groups) require separate parsing tier

### Parser Development
- Use `build.rs` for code generation that needs file system access
- Leverage nom for fixed-format parsing, LALRPOP for recursive grammars
- State tracking is essential for SMILES (ring closures, branches, chirality)
- Error recovery and detailed error reporting are prioritized

### Performance Considerations
- Chemical datasets can be large - optimize for batch processing
- Memory efficiency matters for molecular graph storage
- Build-time optimizations reduce runtime overhead
- Zero-allocation parsing where possible
