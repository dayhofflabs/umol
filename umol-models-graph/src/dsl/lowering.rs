//! AST lowering into IR targets

use super::error::LoweringError;

/// Lowering trait for DSL AST types, config shared by all targets
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

/// AST lowering targets
pub trait FromAst<A: LowerAst>: Sized {
    fn from_ast(ast: A, cfg: &A::Config) -> Result<Self, LoweringError>;
}
