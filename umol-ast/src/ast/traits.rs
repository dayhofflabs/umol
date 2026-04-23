//! AST conversion traits.
//!
//! `FromAst` / `IntoAst` are conversion traits between the AST and external representations.

pub trait FromAst<A>: Sized {
    type Ctx<'a>;
    type Error;
    fn from_ast<'a>(ast: &A, ctx: &Self::Ctx<'a>) -> Result<Self, Self::Error>;
}

pub trait IntoAst<A>: Sized {
    type Ctx<'a>;
    type Error;
    fn into_ast<'a>(self, ctx: &Self::Ctx<'a>) -> Result<A, Self::Error>;
}
