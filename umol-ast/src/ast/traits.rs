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

/// Literal extraction for AST types whose value space includes an
/// `Undetermined` / non-literal branch. `as_lit` returns `Some(lit)` only
/// when the AST is fully resolved to a concrete literal; the derived methods
/// are pure pass-throughs to the corresponding `Option` combinators.
pub trait AsLit {
    /// Concrete literal type (e.g. `i64` for `ValueAst`, `Element` for
    /// `ElementAst`, `SpinState` for `SpinStateAst`).
    type Lit;

    /// `Some(lit)` when the AST resolves to a concrete literal; `None` for
    /// `Undetermined`, expression patterns, sets, sentinels, or physics-invalid
    /// composites.
    fn as_lit(&self) -> Option<Self::Lit>;

    #[inline]
    fn as_lit_ok_or<E>(&self, err: E) -> Result<Self::Lit, E> {
        self.as_lit().ok_or(err)
    }

    #[inline]
    fn as_lit_ok_or_else<E, F: FnOnce() -> E>(&self, err: F) -> Result<Self::Lit, E> {
        self.as_lit().ok_or_else(err)
    }

    #[inline]
    fn as_lit_or(&self, default: Self::Lit) -> Self::Lit {
        self.as_lit().unwrap_or(default)
    }

    #[inline]
    fn as_lit_or_else<F: FnOnce() -> Self::Lit>(&self, default: F) -> Self::Lit {
        self.as_lit().unwrap_or_else(default)
    }

    #[inline]
    fn as_lit_expect(&self, msg: &str) -> Self::Lit {
        self.as_lit().expect(msg)
    }
}
