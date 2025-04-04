# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Test Commands
- Build: `cargo build`
- Run all tests: `cargo test`
- Run a specific test: `cargo test test_name`
- Run tests in a specific module: `cargo test module::tests`
- Run property tests: `cargo test proptest`
- Format code: `cargo fmt`
- Lint code: `cargo clippy`

## Code Style Guidelines
- Use Rust 2021 edition
- Follow standard Rust naming conventions (snake_case for functions/variables, CamelCase for types)
- Organize imports by external crates first, then internal modules
- Use thiserror for error handling with descriptive error messages
- Implement Display trait for user-facing types
- Use builder pattern for complex object construction
- Write comprehensive tests including:
  - Unit tests with rstest for parameterized testing
  - Property-based tests with proptest
- Document public APIs with clear comments
- Use strong typing with enums for domain concepts
- Validate input values and return Result<T, Error> instead of panicking

## Project Architecture
- Four-domain semantic model:
  1. Structure Domain: Chemical structures and transformations between them
  2. Model Domain: Representation models with different capabilities
  3. Instance Domain: Structure-model pairs and operations on them
  4. Property Domain: Properties calculated from instances
- GraphMolecule is the primary data structure for molecular representation
- Builder pattern APIs for molecular construction
- Strong type safety with validation of chemical rules (valence, etc.)
- Favor composition over inheritance, using traits for shared behaviors
- Use petgraph for the underlying graph data structure
- Support for various molecular representations (2D graph, 3D coordinates)