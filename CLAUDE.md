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
- Don't use mod.rs modules, instead use the named modules one level up
  (core.rs instead of core/mod.rs).
- Use thiserror for error handling with descriptive error messages
- Implement Display trait for user-facing types
- Use builder pattern for complex object construction
- For tests use rstest for parametrized testing.
- Use criterion for benchmarking.
- Use strong typing with enums for domain concepts
- Validate input values and return Result<T> instead of panicking.
  Use error and result types from the umol crate.

## Project Architecture
- Four-domain semantic model:
  1. Structure Domain: Chemical structures and transformations between them
  2. Model Domain: Representation models with different capabilities
  3. Instance Domain: Structure-model pairs and operations on them
  4. Property Domain: Properties calculated from instances
- Models are grouped in umol-models-* crates. Each of them represents a specific semantic model
  of molecular structure, e.g., chemical graph, point cloud in 3D space, Born-Oppenheimer model
  (classical point cloud for nuclei, quantum mechanical model for electrons).
- I/O formats are represented by dedicated semantic models. All model conversions should be
  explicit and parametrizable.