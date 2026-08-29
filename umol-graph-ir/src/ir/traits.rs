//! IR-related traits.
//!
//! `FromIr` / `IntoIr` are infallible conversions parameterized by context.
//! `TryFromIr` / `TryIntoIr` are fallible, reject invalid boundary inputs.
//! `AsLit` extracts a literal value from a form.
//! `Lattice` defines the refinement lattice on forms.
//! `Normalize` puts a form into normal fixed-frame representation.
//! `Normalized` is a value carrying the guarantee that it is normalized.

use std::borrow::Cow;
use std::hash::Hash;

use super::error::{Contradiction, NoJoin};

/// Build `Self` from a borrowed IR value of type `A` plus a configuration context.
/// IR → DSL direction. Infallible.
pub trait FromIr<A>: Sized {
    /// Configuration consumed during conversion (e.g. entity defaults).
    type Context;
    fn from_ir(ir: &A, context: &Self::Context) -> Self;
}

/// Consume `self` to produce an IR value of type `A` plus a configuration context.
/// DSL → IR direction. Infallible.
pub trait IntoIr<A>: Sized {
    /// Configuration consumed during conversion (e.g. entity defaults).
    type Context;
    fn into_ir(self, context: &Self::Context) -> A;
}

/// Build `Self` from a borrowed IR value of type `A` plus a configuration context.
/// Fallible variant: the conversion can reject the input.
pub trait TryFromIr<A>: Sized {
    /// Configuration consumed during conversion.
    type Context;
    /// Error returned when the source cannot be represented as `Self`.
    type Error;
    fn try_from_ir(ir: &A, context: &Self::Context) -> Result<Self, Self::Error>;
}

/// Consume `self` to produce an IR value of type `A` plus a configuration context.
/// Fallible variant: used by raising paths whose source carries information
/// without a faithful IR representation (e.g. TableIR Sgroups).
pub trait TryIntoIr<A>: Sized {
    /// Configuration consumed during conversion.
    type Context;
    /// Error returned when `self` cannot be represented as `A`.
    type Error;
    fn try_into_ir(self, context: &Self::Context) -> Result<A, Self::Error>;
}

/// Exact literal projection for graph-IR forms whose value space includes an
/// `Undetermined` or otherwise non-ground branch.
///
/// For a type implementing both `Lattice` and `AsLit`, projection is total
/// exactly on structurally ground values:
/// `value.is_ground() == value.as_lit().is_some()`. It does not normalize,
/// apply defaults, validate domain invariants, or identify distinct normalized
/// ground values that happen to have the same downstream numerical effect.
pub trait AsLit {
    /// Exact non-lattice representation of a ground value.
    type Lit;

    /// `Some(lit)` exactly when the form is structurally ground.
    fn as_lit(&self) -> Option<Self::Lit>;
}

/// Refinement lattice on graph-IR forms.
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
/// iff `target` refines `pattern`, i.e. `pattern.meet(target) == normalized(target)`.
/// It has a `meet`-derived default; impls override it only as a cheaper shortcut.
/// `satisfies` is its receiver-inverted reading for target-side receivers.
pub trait Lattice: Normalize {
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
    /// to `target`'s normal form. Default is `meet`-derived; override only as
    /// a cheaper shortcut that yields the same result.
    fn matches(&self, target: &Self) -> bool {
        match (self.meet(target), target.normalized()) {
            (Some(meet), Ok(target)) => meet == *target,
            _ => false,
        }
    }

    /// Receiver-inverted [`Lattice::matches`]: `target.satisfies(pattern)` is
    /// true iff `pattern.matches(target)`. For receivers standing on the
    /// target side; the default is the definition and is never overridden.
    fn satisfies(&self, pattern: &Self) -> bool {
        pattern.matches(self)
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

/// Normal form of a graph-IR value. The per-type `normalize` puts a
/// value into its one normal form — sorted/deduped sets, singleton collapse,
/// folded decidable expressions — returning `Err(Contradiction)` for an
/// unsatisfiable value (e.g. an empty set).
///
/// Equality is **lazy**: `==`/`Hash`/`Ord` stay derived-structural ("same
/// tree"); semantic equality is [`Equiv::equiv`], comparing normal forms. The hot
/// path is cheap — `normalized` borrows values that are already normalized.
pub trait Normalize: Sized + Clone + PartialEq {
    /// By-value normal form (the folding lives here). Idempotent.
    fn normalize(self) -> Result<Self, Contradiction>;

    /// By-reference normal form: `Cow::Borrowed` when already normalized (the
    /// fast path), else clone + `normalize`. `Err` = unsatisfiable. The
    /// default always clones; override to borrow the trivial variants.
    fn normalized(&self) -> Result<Cow<'_, Self>, Contradiction> {
        Ok(Cow::Owned(self.clone().normalize()?))
    }
}

/// Semantic equality of graph-IR values in their current id and participant frame.
pub trait Equiv: Normalize {
    /// Equal normal forms. Two unsatisfiable values count as equal.
    /// Structural equality short-circuits normalization; otherwise the normal forms
    /// with the derived structural `==` (on `Result<Cow<_>, _>`) — no recursion.
    fn equiv(&self, other: &Self) -> bool {
        self == other || self.normalized() == other.normalized()
    }
}

/// Transport a frame-relative value through an independently supplied compatible action.
///
/// This operation neither normalizes the value nor selects a participant frame. Compatibility is
/// receiver-relative: an implementation checks every action-domain, degree, positional-length, and
/// subgroup condition represented by `self`, and returns `None` when any such condition fails.
/// Information available only from an owning aggregate is checked by that aggregate.
///
/// # Semantic properties
///
/// For every compatible action family, identity leaves the value unchanged, applying an action and
/// its inverse recovers the value, and sequential application agrees with action composition.
pub trait FrameTransport: Sized {
    /// The complete frame action on `Self`.
    type Action;

    /// Restate `self` under `action`, or return `None` when the receiver exposes an incompatibility.
    fn reframe_by(self, action: &Self::Action) -> Option<Self>;
}

/// The frame quotient over an entity family: select a determinate participant frame and restate the
/// frame-relative payload accordingly.
///
/// The family is the carrier that knows both which factor bears a frame and what the payload means,
/// so it owns the quotient rather than the storage shape or the form. One member is required and
/// the other two are laws over it: `reframe` is the selection alone, and `framed_eq` is equality of
/// selected values.
pub trait Reframe: Sized + PartialEq {
    /// The frame action selected for one entry, keyed by the family's own id type. Four families
    /// carry a position order; the two stereo families carry a `Permutation`, since a stereo frame
    /// is bounded by its kind's degree.
    type Action;

    /// Reduce every entry, then present each in its selected frame, returning the action per entry.
    fn reframe_with_action(&self) -> Result<(Self, Vec<Self::Action>), Contradiction>;

    /// Reduce every entry, then present each in its selected frame.
    fn reframe(&self) -> Result<Self, Contradiction> {
        Ok(self.reframe_with_action()?.0)
    }

    /// Equality modulo the stored frame: the reframed values agree. Two unsatisfiable families
    /// count as equal, as they do under [`Equiv`].
    fn framed_eq(&self, other: &Self) -> bool {
        match (self.reframe(), other.reframe()) {
            (Ok(left), Ok(right)) => left == right,
            (Err(_), Err(_)) => true,
            _ => false,
        }
    }
}

impl<T: Normalize> Equiv for T {}

/// A value carrying the guarantee that it is normalized. Built via `new` (which
/// normalizes once); its derived structural `Eq`/`Hash`/`Ord` are therefore
/// *semantic*, so it can key a `HashMap` / `BiBTreeMap` for semantic dedup.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Normalized<T>(T);

impl<T: Normalize> Normalized<T> {
    /// Normalize `value` once; `Err` if it is unsatisfiable.
    pub fn new(value: T) -> Result<Self, Contradiction> {
        Ok(Self(value.normalize()?))
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
    type Constraint: Normalize;

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
