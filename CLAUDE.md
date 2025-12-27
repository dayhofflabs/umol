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
- Put all imports at the top of file. Organize imports by external crates first, then internal modules.
  Exception: tests
- Don't use mod.rs modules, instead use the named modules one level up
  (core.rs instead of core/mod.rs).
- Use thiserror for error handling with descriptive error messages
- Implement Display trait for user-facing types
- Use builder pattern for complex object construction
- For tests use rstest for parametrized testing.
- Use criterion for benchmarking.
- Use strong typing with enums for domain concepts
- Validate input values and return Result<T> instead of panicking.
- Use the following import Conventions
   1. Internal symbols (same project, any crate): should be imported and used unqualified, rename with crate/module prefix only if ambiguous.
   2. Types/Structs/Enums from external crates: import directly, use unqualified, rename with the crate prefix in case of naming clashes
   3. Functions/Constants from external crates: import parent module, call qualified (`io::stdin()`)
   4. Enum variants: import enum, do not import variants directly. Exceptions: In large match blocks or tests.
- Write comments for an external reader, do not refer to aspects of editing or refactoring history
- Refrain from comment ornamentation. Do not introduce separator lines or headings like ========== in comments.
  If the file is so long that it warrants headings, it should be split into multiple files instead.
- Avoid aliasing of types defined in the codebase. Aliasing of core language or library types is ok if it serves
  semantic clarity.

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