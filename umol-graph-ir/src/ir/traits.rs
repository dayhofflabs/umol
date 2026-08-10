//! IR-related traits.
//!
//! `FromIr` / `IntoIr` are infallible conversions parameterized by context.
//! `TryFromIr` / `TryIntoIr` are fallible, reject invalid boundary inputs.
//! `AsLit` extracts a literal value from an AST type.
//! `Lattice` is a refinement lattice on AST value types.
//! `Canonicalize` is a canonical form of an AST value.
//! `Canonical` is a value carrying the guarantee that it is canonical.

use std::borrow::Cow;
use std::hash::Hash;

use umol_graph_core::{BiRelationData, ParticipantPosition, RelationData};

use super::error::{Contradiction, NoJoin};

/// Build `Self` from a borrowed IR value of type `A` plus a configuration context.
/// IR → DSL direction. Infallible.
pub trait FromIr<A>: Sized {
    /// Configuration consumed during conversion (e.g. entity defaults).
    type Ctx;
    fn from_ir(ir: &A, ctx: &Self::Ctx) -> Self;
}

/// Consume `self` to produce an IR value of type `A` plus a configuration context.
/// DSL → IR direction. Infallible.
pub trait IntoIr<A>: Sized {
    /// Configuration consumed during conversion (e.g. entity defaults).
    type Ctx;
    fn into_ir(self, ctx: &Self::Ctx) -> A;
}

/// Build `Self` from a borrowed IR value of type `A` plus a configuration context.
/// Fallible variant: the conversion can reject the input.
pub trait TryFromIr<A>: Sized {
    type Ctx;
    type Error;
    fn try_from_ir(ir: &A, ctx: &Self::Ctx) -> Result<Self, Self::Error>;
}

/// Consume `self` to produce an IR value of type `A` plus a configuration context.
/// Fallible variant: used by raising paths whose source carries information
/// without a faithful IR representation (e.g. TableIR Sgroups).
pub trait TryIntoIr<A>: Sized {
    type Ctx;
    type Error;
    fn try_into_ir(self, ctx: &Self::Ctx) -> Result<A, Self::Error>;
}

/// Exact literal projection for AST types whose value space includes an
/// `Undetermined` or otherwise non-ground branch.
///
/// For a type implementing both `Lattice` and `AsLit`, projection is total
/// exactly on structurally ground values:
/// `value.is_ground() == value.as_lit().is_some()`. It does not canonicalize,
/// apply defaults, validate domain invariants, or identify distinct canonical
/// ground values that happen to have the same downstream numerical effect.
pub trait AsLit {
    /// Exact non-lattice representation of a ground value.
    type Lit;

    /// `Some(lit)` exactly when the AST value is structurally ground.
    fn as_lit(&self) -> Option<Self::Lit>;
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

    /// Least upper bound. `Err(NoJoin)` when `self` and `other` have no common
    /// generalization — a top-less (meet-semilattice) type whose operands lie in
    /// different fibers (e.g. two `AtomConstraintForm`s of different kind). Bounded
    /// lattices always return `Ok`.
    fn join(&self, other: &Self) -> Result<Self, NoJoin>;

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

    /// In-place `join`. Returns `Ok(true)` iff `self` actually changed;
    /// `Err(NoJoin)` (leaving `self` unchanged) when the join does not exist.
    fn widen_with(&mut self, other: &Self) -> Result<bool, NoJoin> {
        let new = self.join(other)?;
        Ok(if new != *self {
            *self = new;
            true
        } else {
            false
        })
    }
}

/// Normal (canonical) form of an AST value. The per-type `canonicalize` puts a
/// value into its one canonical form — sorted/deduped sets, singleton collapse,
/// folded decidable expressions — returning `Err(Contradiction)` for an
/// unsatisfiable value (e.g. an empty set).
///
/// Equality is **lazy**: `==`/`Hash`/`Ord` stay derived-structural ("same
/// tree"); semantic equality is `canonical_eq`, comparing canonical forms. The hot
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
    fn canonical_eq(&self, other: &Self) -> bool {
        self == other || self.canonical() == other.canonical()
    }
}

/// Framed value equivalence for a single-factor relation payload: the value axis (`canonical_eq`)
/// composed with the position axis (`on_permutation`). `equiv` is the frameless case; `equiv_under`
/// reindexes `self` into `other`'s frame first, skipping the work when the payload is permutation-invariant.
pub trait Equiv: RelationData + Canonicalize {
    fn equiv(&self, other: &Self) -> bool {
        self.canonical_eq(other)
    }

    fn equiv_under(&self, other: &Self, order: &[ParticipantPosition]) -> bool {
        if self.is_permutation_invariant() {
            self.canonical_eq(other)
        } else {
            let mut probe = self.clone();
            probe.on_permutation(order);
            probe.canonical_eq(other)
        }
    }
}

impl<T: RelationData + Canonicalize> Equiv for T {}

/// Two-factor analog of [`Equiv`] for a birelation payload: `equiv_under` reindexes `self` per factor
/// before comparing canonical forms.
pub trait BiEquiv: BiRelationData + Canonicalize {
    fn equiv(&self, other: &Self) -> bool {
        self.canonical_eq(other)
    }

    fn equiv_under(
        &self,
        other: &Self,
        order_1: &[ParticipantPosition],
        order_2: &[ParticipantPosition],
    ) -> bool {
        if self.is_permutation_invariant() {
            self.canonical_eq(other)
        } else {
            let mut probe = self.clone();
            probe.on_permutation(order_1, order_2);
            probe.canonical_eq(other)
        }
    }
}

impl<T: BiRelationData + Canonicalize> BiEquiv for T {}

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
/// states (`Attributes`): `apply` is its action — carrying a state forward by a `ModifyField` /
/// `ModifyConstraint` delta — and `diff` is the inverse, factoring two states back into the deltas
/// between them, with `apply(lhs, diff(lhs, rhs)) == rhs`.
pub trait EntityPatch: Sized {
    type Id: Copy + Eq + Hash + From<usize>;
    type Attributes: Clone;
    type FieldChange;
    type Constraint: Canonicalize;

    fn modify_field(id: Self::Id, change: Self::FieldChange) -> Self;
    fn modify_constraint(
        id: Self::Id,
        old: Option<Self::Constraint>,
        new: Option<Self::Constraint>,
    ) -> Self;
    fn apply_field(
        attributes: &mut Self::Attributes,
        change: Self::FieldChange,
    ) -> Result<(), Contradiction>;
    fn diff_field(lhs: &Self::Attributes, rhs: &Self::Attributes) -> Vec<Self::FieldChange>;
    fn apply_constraint(
        attributes: &mut Self::Attributes,
        old: Option<Self::Constraint>,
        new: Option<Self::Constraint>,
    ) -> Result<(), Contradiction>;
    #[allow(clippy::type_complexity)]
    fn diff_constraints(
        lhs: &Self::Attributes,
        rhs: &Self::Attributes,
    ) -> Vec<(Option<Self::Constraint>, Option<Self::Constraint>)>;

    /// The `ModifyField` / `ModifyConstraint` deltas carrying `lhs` to `rhs` for one entity — the
    /// inverse of `apply_*_change`, recovering a `Modified` entity's deltas.
    fn diff(id: Self::Id, lhs: &Self::Attributes, rhs: &Self::Attributes) -> Vec<Self> {
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
