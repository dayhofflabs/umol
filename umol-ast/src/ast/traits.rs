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

/// Refinement lattice on AST value types.
///
/// `Undetermined` is the top (most general / "any value"); fully ground
/// concrete values are bottom (most specific). `meet` is the greatest lower
/// bound (most specific common refinement; `None` when incompatible);
/// `join` is the least upper bound (most general common generalization).
///
/// `narrow_from` and `widen_with` are the in-place counterparts of `meet`
/// and `join`; both return `true` iff `self` actually changed.
///
/// `matches` is the partial-order check: `pattern.matches(target)` is true
/// iff `target` refines `pattern`, i.e. `pattern.meet(target) == Some(target)`
/// up to canonicalization of set-valued representations. Each impl provides
/// its own direct implementation rather than delegating to `meet`.
pub trait Lattice: Sized + Clone + PartialEq {
    /// Top of the lattice — `self` carries no value information.
    fn is_undetermined(&self) -> bool;

    /// Bottom of the lattice — `self` resolves to a single concrete value.
    fn is_ground(&self) -> bool;

    /// Greatest lower bound. `None` when `self` and `other` are mutually
    /// incompatible (no value can satisfy both).
    fn meet(&self, other: &Self) -> Option<Self>;

    /// Least upper bound. Total — incompatible pairs widen toward the
    /// nearest common generalization (typically `Undetermined`).
    fn join(&self, other: &Self) -> Self;

    /// Partial-order check: `self` (pattern) is true on `target` iff every
    /// value `target` admits is also admitted by `self`. Equivalent to
    /// `self.meet(target) == Some(target)` modulo canonicalization.
    fn matches(&self, target: &Self) -> bool;

    /// Symmetric compatibility: there exists a ground value that refines
    /// both `self` and `other`. Equivalent to `self.meet(other).is_some()`;
    /// override when a direct check is cheaper than constructing the meet.
    fn is_compatible(&self, other: &Self) -> bool {
        self.meet(other).is_some()
    }

    /// In-place `meet`. Returns `true` iff `self` actually changed. When
    /// `self` and `other` are incompatible, leaves `self` unchanged and
    /// returns `false`.
    fn narrow_from(&mut self, other: &Self) -> bool {
        match self.meet(other) {
            Some(new) if new != *self => {
                *self = new;
                true
            }
            _ => false,
        }
    }

    /// In-place `join`. Returns `true` iff `self` actually changed.
    fn widen_with(&mut self, other: &Self) -> bool {
        let new = self.join(other);
        if new != *self {
            *self = new;
            true
        } else {
            false
        }
    }

    /// Cross-field constraint propagation. Default no-op for types without
    /// relational constraints. Types with relational constraints (e.g.,
    /// `AtomAst` carrying `JointDomain` entries) override this to project
    /// constraints across fields. `Err(Contradiction)` signals that no
    /// admissible value assignment remains; the derived `meet` converts this
    /// to `None` so it propagates through the standard `Option<Self>`
    /// contract.
    fn saturate(&mut self) -> Result<(), super::error::Contradiction> {
        Ok(())
    }
}
