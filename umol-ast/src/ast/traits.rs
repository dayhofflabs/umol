//! AST-related traits.
//!
//! `FromAst` / `IntoAst` are infallible conversions parameterized by context.
//! `TryFromAst` / `TryIntoAst` are fallible, reject invalid boundary inputs.
//! `AsLit` extracts a literal value from an AST type.
//! `Lattice` is a refinement lattice on AST value types.
//! `Canonicalize` is a canonical form of an AST value.
//! `Canonical` is a value carrying the guarantee that it is canonical.

use std::borrow::Cow;
use std::hash::Hash;

use super::error::Contradiction;

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

    /// `true` iff `self` resolves to exactly the literal `value`. Non-canonicalizing,
    /// mirroring `as_lit`: a non-literal that would fold to `value` is not matched.
    /// The clean replacement for the old ad-hoc `matches_value`.
    #[inline]
    fn as_lit_matches(&self, value: Self::Lit) -> bool
    where
        Self::Lit: PartialEq,
    {
        self.as_lit() == Some(value)
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
/// iff `target` refines `pattern`, i.e. `pattern.meet(target) == canonical(target)`.
/// It has a `meet`-derived default; impls override it only as a cheaper shortcut.
pub trait Lattice: Canonicalize {
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
    /// value `target` admits is also admitted by `self` — i.e. the meet refines
    /// to `target`'s canonical form. Default is `meet`-derived; override only as
    /// a cheaper shortcut that yields the same result.
    fn matches(&self, target: &Self) -> bool {
        match (self.meet(target), target.canonical()) {
            (Some(meet), Ok(target)) => meet == *target,
            _ => false,
        }
    }

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
}

/// Normal (canonical) form of an AST value. The per-type `canonicalize` puts a
/// value into its one canonical form — sorted/deduped sets, singleton collapse,
/// folded decidable expressions — returning `Err(Contradiction)` for an
/// unsatisfiable value (e.g. an empty set).
///
/// Equality is **lazy**: `==`/`Hash`/`Ord` stay derived-structural ("same
/// tree"); semantic equality is `equiv`, comparing canonical forms. The hot
/// path is cheap — `canonical` borrows values that are already canonical.
pub trait Canonicalize: Sized + Clone + PartialEq {
    /// By-value canonical form (the folding lives here). Idempotent.
    fn canonicalize(self) -> Result<Self, Contradiction>;

    /// By-reference canonical form: `Cow::Borrowed` when already canonical (the
    /// fast path), else clone + `canonicalize`. `Err` = unsatisfiable. The
    /// default always clones; override to borrow the trivial variants.
    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        Ok(Cow::Owned(self.clone().canonicalize()?))
    }

    /// Semantic equality: equal canonical forms (two unsatisfiable values count
    /// as equal). Structural short-circuit first, then compare canonical forms
    /// with the derived structural `==` (on `Result<Cow<_>, _>`) — no recursion.
    fn equiv(&self, other: &Self) -> bool {
        self == other || self.canonical() == other.canonical()
    }
}

/// A value carrying the guarantee that it is canonical. Built via `new` (which
/// canonicalizes once); its derived structural `Eq`/`Hash`/`Ord` are therefore
/// *semantic*, so it can key a `HashMap` / `BiBTreeMap` for semantic dedup.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Canonical<T>(T);

impl<T: Canonicalize> Canonical<T> {
    /// Canonicalize `value` once; `Err` if it is unsatisfiable.
    pub fn new(value: T) -> Result<Self, Contradiction> {
        Ok(Self(value.canonicalize()?))
    }

    pub fn get(&self) -> &T {
        &self.0
    }

    pub fn into_inner(self) -> T {
        self.0
    }
}

/// The patch algebra over one entity's fields and constraints. A delta is the morphism between two
/// states (`Ast`s): `apply` is its action — carrying a state forward by a `ModifyField` /
/// `ModifyConstraint` delta — and `diff` is the inverse, factoring two states back into the deltas
/// between them, with `apply(lhs, diff(lhs, rhs)) == rhs`.
pub trait EntityPatch: Sized {
    type Id: Copy + Eq + Hash + From<usize>;
    type Ast: Clone;
    type FieldChange;
    type Constraint: Clone + PartialEq;

    fn modify_field(id: Self::Id, change: Self::FieldChange) -> Self;
    fn modify_constraint(
        id: Self::Id,
        old: Option<Self::Constraint>,
        new: Option<Self::Constraint>,
    ) -> Self;
    fn apply_field(ast: &mut Self::Ast, change: Self::FieldChange) -> Result<(), Contradiction>;
    fn diff_field(lhs: &Self::Ast, rhs: &Self::Ast) -> Vec<Self::FieldChange>;
    fn apply_constraint(
        ast: &mut Self::Ast,
        old: Option<Self::Constraint>,
        new: Option<Self::Constraint>,
    ) -> Result<(), Contradiction>;
    #[allow(clippy::type_complexity)]
    fn diff_constraints(
        lhs: &Self::Ast,
        rhs: &Self::Ast,
    ) -> Vec<(Option<Self::Constraint>, Option<Self::Constraint>)>;

    /// The `ModifyField` / `ModifyConstraint` deltas carrying `lhs` to `rhs` for one entity — the
    /// inverse of `apply_*_change`, recovering a `Modified` entity's deltas.
    fn diff(id: Self::Id, lhs: &Self::Ast, rhs: &Self::Ast) -> Vec<Self> {
        let mut out: Vec<Self> = Self::diff_field(lhs, rhs)
            .into_iter()
            .map(|change| Self::modify_field(id, change))
            .collect();
        out.extend(
            Self::diff_constraints(lhs, rhs)
                .into_iter()
                .map(|(old, new)| Self::modify_constraint(id, old, new)),
        );
        out
    }
}
