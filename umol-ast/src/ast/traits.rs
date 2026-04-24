//! AST conversion traits.
//!
//! `FromAst` / `IntoAst` are conversion traits between the AST and external representations.

/// Build `Self` from a borrowed AST of type `A` plus a configuration context.
/// Conventionally the AST → DSL (rendering / lowering) direction.
pub trait FromAst<A>: Sized {
    /// Configuration consumed during conversion (e.g. entity defaults).
    type Ctx;
    type Error;
    fn from_ast(ast: &A, ctx: &Self::Ctx) -> Result<Self, Self::Error>;
}

/// Consume `self` to produce an AST of type `A` plus a configuration context.
/// Conventionally the DSL → AST (parsing / raising) direction.
pub trait IntoAst<A>: Sized {
    /// Configuration consumed during conversion (e.g. entity defaults).
    type Ctx;
    type Error;
    fn into_ast(self, ctx: &Self::Ctx) -> Result<A, Self::Error>;
}
