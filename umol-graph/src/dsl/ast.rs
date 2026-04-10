//! AST lowering and raising traits.

use super::error::LoweringError;

/// DSL AST marker trait; carries the config type shared by lowering and raising.
pub trait DslAst {
    type Config;
}

/// AST lowering targets.
pub trait FromAst<A: DslAst>: Sized {
    fn from_ast(ast: &A, cfg: &A::Config) -> Result<Self, LoweringError>;
}

/// Raise a ground or pattern type to its AST representation, for formatting.
///
/// Config-aware: uses the same [`DslAst::Config`] as lowering to determine which
/// fields are implied defaults and can be omitted from the AST.
pub trait ToAst<A: DslAst> {
    fn to_ast(&self, cfg: &A::Config) -> A;
}
