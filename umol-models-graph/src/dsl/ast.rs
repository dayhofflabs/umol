//! AST lowering and raising traits.

use super::error::LoweringError;

/// Lowering trait for DSL AST types; carries the config type shared by all lowering targets.
pub trait LowerAst {
    type Config: Default;

    fn lower<T: FromAst<Self>>(&self) -> Result<T, LoweringError>
    where
        Self: Sized + Clone,
    {
        T::from_ast(self.clone(), &Self::Config::default())
    }

    fn lower_with<T: FromAst<Self>>(&self, cfg: &Self::Config) -> Result<T, LoweringError>
    where
        Self: Sized + Clone,
    {
        T::from_ast(self.clone(), cfg)
    }

    fn lower_into<T: FromAst<Self>>(self) -> Result<T, LoweringError>
    where
        Self: Sized,
    {
        T::from_ast(self, &Self::Config::default())
    }

    fn lower_into_with<T: FromAst<Self>>(self, cfg: &Self::Config) -> Result<T, LoweringError>
    where
        Self: Sized,
    {
        T::from_ast(self, cfg)
    }
}

/// AST lowering targets.
pub trait FromAst<A: LowerAst>: Sized {
    fn from_ast(ast: A, cfg: &A::Config) -> Result<Self, LoweringError>;
}

/// Raise a ground or pattern type to its AST representation, for formatting.
///
/// Unlike [`FromAst`], raising is infallible and requires no configuration —
/// every value maps unambiguously to an AST node.
pub trait ToAst<A> {
    fn to_ast(&self) -> A;
}
