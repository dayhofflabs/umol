//! AST conversion traits.
//!
//! `FromAst` / `IntoAst` are conversion traits between the AST and external representations.

pub trait FromAst<A>: Sized {
    type Ctx;
    type Error;
    fn from_ast(ast: &A, ctx: &Self::Ctx) -> Result<Self, Self::Error>;
}

pub trait IntoAst<A>: Sized {
    type Ctx;
    type Error;
    fn into_ast(self, ctx: &Self::Ctx) -> Result<A, Self::Error>;
}
