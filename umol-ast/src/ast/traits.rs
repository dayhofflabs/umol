//! AST conversion traits.
//!
//! `FromAst` / `IntoAst` are the infallible pair: lossless conversions between
//! AST and an external representation, parameterized by a configuration
//! context. Used for the DSL ↔ AST pair, where lower/raise is a structural
//! transformation that cannot fail.
//!
//! `TryFromAst` / `TryIntoAst` are the fallible pair: conversions that can
//! reject the input on chemistry or representation grounds. Used for
//! `TableIR → MoleculeAst` raising (e.g., `ExtendedMolecule` Sgroups have no
//! faithful AST representation).
//!
//! The split mirrors `From` / `TryFrom` in `std`. There is no blanket impl
//! between the pairs.

/// Build `Self` from a borrowed AST of type `A` plus a configuration context.
/// AST → DSL direction. Infallible.
pub trait FromAst<A>: Sized {
    /// Configuration consumed during conversion (e.g. entity defaults).
    type Ctx;
    fn from_ast(ast: &A, ctx: &Self::Ctx) -> Self;
}

/// Consume `self` to produce an AST of type `A` plus a configuration context.
/// DSL → AST direction. Infallible.
pub trait IntoAst<A>: Sized {
    /// Configuration consumed during conversion (e.g. entity defaults).
    type Ctx;
    fn into_ast(self, ctx: &Self::Ctx) -> A;
}

/// Build `Self` from a borrowed AST of type `A` plus a configuration context.
/// Fallible variant: the conversion can reject the input.
pub trait TryFromAst<A>: Sized {
    type Ctx;
    type Error;
    fn try_from_ast(ast: &A, ctx: &Self::Ctx) -> Result<Self, Self::Error>;
}

/// Consume `self` to produce an AST of type `A` plus a configuration context.
/// Fallible variant: used by raising paths whose source carries information
/// without a faithful AST representation (e.g. TableIR Sgroups).
pub trait TryIntoAst<A>: Sized {
    type Ctx;
    type Error;
    fn try_into_ast(self, ctx: &Self::Ctx) -> Result<A, Self::Error>;
}
