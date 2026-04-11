//! Structural AST: a generalization that can express both ground molecules and
//! pattern queries. AST types are independent of any particular text format;
//! the [`crate::dsl`] module provides one EDN-based serialization.

pub mod atom;
pub mod bond;
pub mod config;
pub mod constraint;
pub mod error;
pub mod molecule;

use error::LoweringError;

/// AST marker trait carrying the lowering/raising config type.
pub trait Ast {
    type Config;
}

/// Lowering targets that can be built from an AST.
pub trait FromAst<A: Ast>: Sized {
    fn from_ast(ast: &A, cfg: &A::Config) -> Result<Self, LoweringError>;
}

/// Raise a ground or pattern type to its AST representation.
pub trait ToAst<A: Ast> {
    fn to_ast(&self, cfg: &A::Config) -> A;
}
