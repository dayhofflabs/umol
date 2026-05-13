//! Atom-level AST fragments shared across crates.

use std::mem;
use std::ops::{Add, Div, Mul, Sub};

use umol_shared::element::Element;

use super::constraint::{
    AromaticValenceAst, AtomConstraint, AtomConstraintKind, AtomConstraints, MulticenterValenceAst,
};
use super::spin::SpinStateAst;
use super::value::{Expr, ValueAst};

/// Atom AST: structural representation of an atom plus the atom-level
/// constraints (valence, degree, ring membership, etc.) that pattern
/// against the surrounding topology.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomAst {
    pub element: ElementAst,
    pub isotope_mass: IsotopeAst,
    pub charge: ValueAst,
    pub implicit_hydrogens: ImplicitHydrogensAst,
    pub lone_pairs: ValueAst,
    pub spin: SpinStateAst,
    pub constraints: AtomConstraints,
}

impl AtomAst {
    pub fn new(element: ElementAst) -> Self {
        Self {
            element,
            ..Default::default()
        }
    }

    pub fn from_element(element: Element) -> Self {
        Self::new(ElementAst::Lit(element))
    }

    pub fn with_element(mut self, element: impl Into<ElementAst>) -> Self {
        self.element = element.into();
        self
    }

    pub fn with_isotope_mass(mut self, mass: impl Into<IsotopeAst>) -> Self {
        self.isotope_mass = mass.into();
        self
    }

    pub fn with_charge(mut self, charge: impl Into<ValueAst>) -> Self {
        self.charge = charge.into();
        self
    }

    pub fn with_implicit_hydrogens(mut self, hydrogens: impl Into<ImplicitHydrogensAst>) -> Self {
        self.implicit_hydrogens = hydrogens.into();
        self
    }

    pub fn with_lone_pairs(mut self, lone_pairs: impl Into<ValueAst>) -> Self {
        self.lone_pairs = lone_pairs.into();
        self
    }

    pub fn with_spin(mut self, spin: impl Into<SpinStateAst>) -> Self {
        self.spin = spin.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// kind (last-wins per `AtomConstraints::add`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<AtomConstraint>) -> Self {
        self.constraints.add(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `AtomConstraints::add`). Does not
    /// clear existing constraints; use `atom.constraints.clear()` or direct
    /// field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<AtomConstraint>,
    {
        for c in constraints {
            self.constraints.add(c.into());
        }
        self
    }

    /// Fill `Undetermined` value-bearing struct fields with zero defaults
    /// (isotope→Natural, charge / implicit hydrogens / lone pairs → 0, spin →
    /// closed-shell singlet). Existing literal or expression values and all
    /// constraints are preserved. The result is ground iff `element` is
    /// already ground.
    pub fn into_ground(mut self) -> Self {
        if self.isotope_mass.is_undetermined() {
            self.isotope_mass = IsotopeAst::Natural;
        }
        if self.charge.is_undetermined() {
            self.charge = ValueAst::Lit(0);
        }
        if self.implicit_hydrogens.is_undetermined() {
            self.implicit_hydrogens = ImplicitHydrogensAst::Lit(0);
        }
        if self.lone_pairs.is_undetermined() {
            self.lone_pairs = ValueAst::Lit(0);
        }
        if self.spin.is_undetermined() {
            self.spin = SpinStateAst::from((0_u8, 1_u8));
        }
        self
    }

    /// `into_ground()` plus chemistry-default constraints for an isolated,
    /// non-bonded, non-aromatic, non-multicenter atom: `Valence(0)`,
    /// `DonatedPairs(0)`, `AcceptedPairs(0)`, `MulticenterValence(NotMulticenter)`,
    /// `AromaticValence(NotAromatic)`. Each is added only if the corresponding
    /// constraint kind is not already present; existing entries are preserved.
    /// Matches the `atom_zeroed!` / `mol_zeroed!` macro semantics.
    pub fn into_zeroed(mut self) -> Self {
        self = self.into_ground();
        if !self.constraints.contains(AtomConstraintKind::Valence) {
            self.constraints
                .add(AtomConstraint::Valence(ValueAst::Lit(0)));
        }
        if !self.constraints.contains(AtomConstraintKind::DonatedPairs) {
            self.constraints
                .add(AtomConstraint::DonatedPairs(ValueAst::Lit(0)));
        }
        if !self.constraints.contains(AtomConstraintKind::AcceptedPairs) {
            self.constraints
                .add(AtomConstraint::AcceptedPairs(ValueAst::Lit(0)));
        }
        if !self
            .constraints
            .contains(AtomConstraintKind::MulticenterValence)
        {
            self.constraints.add(AtomConstraint::MulticenterValence(
                MulticenterValenceAst::NotMulticenter,
            ));
        }
        if !self
            .constraints
            .contains(AtomConstraintKind::AromaticValence)
        {
            self.constraints.add(AtomConstraint::AromaticValence(
                AromaticValenceAst::NotAromatic,
            ));
        }
        self
    }

    pub fn is_ground(&self) -> bool {
        self.element.is_ground()
            && self.isotope_mass.is_ground()
            && self.charge.is_ground()
            && self.implicit_hydrogens.is_ground()
            && self.lone_pairs.is_ground()
            && self.spin.is_ground()
    }

    /// `self` (pattern) matches `target` iff every admissible assignment
    /// of `target` is also admissible by `self`, checked field-wise.
    /// See per-field `matches` for the scalar rules.
    pub fn matches(&self, target: &AtomAst) -> bool {
        self.element.matches(&target.element)
            && self.isotope_mass.matches(&target.isotope_mass)
            && self.charge.matches(&target.charge)
            && self.implicit_hydrogens.matches(&target.implicit_hydrogens)
            && self.lone_pairs.matches(&target.lone_pairs)
            && self.spin.matches(&target.spin)
    }

    /// Simplify every value-bearing field in place: `isotope_mass`,
    /// `charge`, `implicit_hydrogens`, `lone_pairs`, both `spin` slots,
    /// and each constraint. `element` has no value to simplify.
    pub fn simplify_values(&mut self) {
        self.isotope_mass = mem::take(&mut self.isotope_mass).simplify();
        self.charge = mem::take(&mut self.charge).simplify();
        self.implicit_hydrogens = mem::take(&mut self.implicit_hydrogens).simplify();
        self.lone_pairs = mem::take(&mut self.lone_pairs).simplify();
        self.spin.simplify_values();
        self.constraints.simplify_each();
    }
}

/// Element expressions
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ElementAst {
    Lit(Element),
    #[default]
    Undetermined,
    Set(Vec<Element>),
    Bind {
        id: String,
        set: Vec<Element>,
    },
    Ref(String),
}

impl From<Element> for ElementAst {
    fn from(element: Element) -> Self {
        Self::Lit(element)
    }
}

impl ElementAst {
    pub fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    /// The element this value denotes when ground; `None` otherwise. Aligned
    /// with [`Self::is_ground`].
    #[inline]
    pub fn literal(&self) -> Option<Element> {
        match self {
            Self::Lit(e) => Some(*e),
            _ => None,
        }
    }

    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    /// Pattern matches target iff every element the target admits is also
    /// admitted by the pattern (superset semantics).
    pub fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Ref(_), _) | (_, Self::Ref(_)) => false,
            (Self::Lit(p), Self::Lit(t)) => p == t,
            (Self::Lit(p), Self::Set(ts) | Self::Bind { set: ts, .. }) => ts.iter().all(|t| t == p),
            (Self::Set(ps) | Self::Bind { set: ps, .. }, Self::Lit(t)) => ps.contains(t),
            (
                Self::Set(ps) | Self::Bind { set: ps, .. },
                Self::Set(ts) | Self::Bind { set: ts, .. },
            ) => ts.iter().all(|t| ps.contains(t)),
        }
    }
}

/// Isotope-mass expressions. `Natural` denotes the naturally most abundant
/// isotope (`#i=`); numeric variants mirror `ValueAst` and are flattened here
/// to keep `Undetermined` as a single top-level state.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IsotopeAst {
    #[default]
    Undetermined,
    Natural,
    Lit(i64),
    LitSet(Vec<i64>),
    Expr(Expr),
}

impl IsotopeAst {
    /// Semantic ground: `Natural` and `Lit(_)` are the primary ground forms;
    /// `LitSet`/`Expr` delegate through the shared helpers so constant-valued
    /// expressions and singleton sets are also ground. Same fast-path /
    /// cold-slow-path split as [`ValueAst::is_ground`].
    #[inline]
    pub fn is_ground(&self) -> bool {
        match self {
            Self::Natural | Self::Lit(_) => true,
            Self::Undetermined => false,
            _ => self.is_ground_slow(),
        }
    }

    #[inline(never)]
    #[cold]
    fn is_ground_slow(&self) -> bool {
        match self {
            Self::LitSet(s) => super::value::litset_is_ground(s),
            Self::Expr(e) => e.is_ground(),
            Self::Natural | Self::Lit(_) | Self::Undetermined => unreachable!(),
        }
    }

    /// The mass number this value denotes when ground; `None` otherwise.
    /// `Natural` returns `Some(0)` as the sentinel for "natural isotopic
    /// abundance — no specific mass committed". Aligned with
    /// [`Self::is_ground`].
    #[inline]
    pub fn literal(&self) -> Option<u32> {
        match self {
            Self::Natural => Some(0),
            Self::Lit(n) => u32::try_from(*n).ok(),
            Self::Undetermined => None,
            _ => self.literal_slow(),
        }
    }

    #[inline(never)]
    #[cold]
    fn literal_slow(&self) -> Option<u32> {
        match self {
            Self::LitSet(s) => super::value::litset_is_ground(s)
                .then(|| u32::try_from(s[0]).ok())
                .flatten(),
            Self::Expr(e) => e
                .evaluate_checked(&super::value::Bindings::new())
                .and_then(|n| u32::try_from(n).ok()),
            Self::Natural | Self::Lit(_) | Self::Undetermined => unreachable!(),
        }
    }

    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    pub fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Natural, Self::Natural) => true,
            (Self::Natural, _) | (_, Self::Natural) => false,
            (p, t) => p.as_value().matches(&t.as_value()),
        }
    }

    fn as_value(&self) -> ValueAst {
        match self {
            Self::Undetermined => ValueAst::Undetermined,
            Self::Lit(n) => ValueAst::Lit(*n),
            Self::LitSet(s) => ValueAst::LitSet(s.clone()),
            Self::Expr(e) => ValueAst::Expr(e.clone()),
            Self::Natural => unreachable!("Natural handled before as_value"),
        }
    }

    /// Simplify the inner `Expr` of `Expr(_)` and lift `Expr(Lit(n))` /
    /// `Expr(Neg(Lit(n)))` to `Lit(n)` / `Lit(-n)` (mirrors
    /// `ValueAst::simplify`). Other variants are unchanged.
    pub fn simplify(self) -> Self {
        match self {
            Self::Expr(e) => Self::from(ValueAst::Expr(e).simplify()),
            other => other,
        }
    }
}

impl From<ValueAst> for IsotopeAst {
    fn from(v: ValueAst) -> Self {
        match v {
            ValueAst::Undetermined => Self::Undetermined,
            ValueAst::Lit(n) => Self::Lit(n),
            ValueAst::LitSet(s) => Self::LitSet(s),
            ValueAst::Expr(e) => Self::Expr(e),
        }
    }
}

impl From<i64> for IsotopeAst {
    fn from(value: i64) -> Self {
        Self::Lit(value)
    }
}

/// Implicit hydrogen expressions. `Normal` denotes the valence-model default
/// (`#h=`); numeric variants mirror `ValueAst` and are flattened here to keep
/// `Undetermined` as a single top-level state.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ImplicitHydrogensAst {
    #[default]
    Undetermined,
    Normal,
    Lit(i64),
    LitSet(Vec<i64>),
    Expr(Expr),
}

impl ImplicitHydrogensAst {
    /// Semantic ground: only `Lit(_)`, ground singleton `LitSet`, and ground
    /// `Expr` count. `Normal` is **not** ground — it's a placeholder for
    /// "compute via valence model"; the resolver lowers it to `Lit(n)`.
    #[inline]
    pub fn is_ground(&self) -> bool {
        match self {
            Self::Lit(_) => true,
            Self::Normal | Self::Undetermined => false,
            _ => self.is_ground_slow(),
        }
    }

    #[inline(never)]
    #[cold]
    fn is_ground_slow(&self) -> bool {
        match self {
            Self::LitSet(s) => super::value::litset_is_ground(s),
            Self::Expr(e) => e.is_ground(),
            Self::Lit(_) | Self::Normal | Self::Undetermined => unreachable!(),
        }
    }

    /// The single integer this value denotes when ground; `None` otherwise
    /// (including `Normal`). Aligned with [`Self::is_ground`].
    #[inline]
    pub fn literal(&self) -> Option<i64> {
        match self {
            Self::Lit(n) => Some(*n),
            Self::Normal | Self::Undetermined => None,
            _ => self.literal_slow(),
        }
    }

    #[inline(never)]
    #[cold]
    fn literal_slow(&self) -> Option<i64> {
        match self {
            Self::LitSet(s) => super::value::litset_is_ground(s).then(|| s[0]),
            Self::Expr(e) => e.evaluate_checked(&super::value::Bindings::new()),
            Self::Lit(_) | Self::Normal | Self::Undetermined => unreachable!(),
        }
    }

    pub fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    pub fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Normal, Self::Normal) => true,
            (Self::Normal, _) | (_, Self::Normal) => false,
            (p, t) => p.as_value().matches(&t.as_value()),
        }
    }

    fn as_value(&self) -> ValueAst {
        match self {
            Self::Undetermined => ValueAst::Undetermined,
            Self::Lit(n) => ValueAst::Lit(*n),
            Self::LitSet(s) => ValueAst::LitSet(s.clone()),
            Self::Expr(e) => ValueAst::Expr(e.clone()),
            Self::Normal => unreachable!("Normal handled before as_value"),
        }
    }

    /// Simplify the inner `Expr` of `Expr(_)` and lift `Expr(Lit(n))` /
    /// `Expr(Neg(Lit(n)))` to `Lit(n)` / `Lit(-n)`. Other variants are
    /// unchanged.
    pub fn simplify(self) -> Self {
        match self {
            Self::Expr(e) => Self::from(ValueAst::Expr(e).simplify()),
            other => other,
        }
    }
}

impl From<ValueAst> for ImplicitHydrogensAst {
    fn from(v: ValueAst) -> Self {
        match v {
            ValueAst::Undetermined => Self::Undetermined,
            ValueAst::Lit(n) => Self::Lit(n),
            ValueAst::LitSet(s) => Self::LitSet(s),
            ValueAst::Expr(e) => Self::Expr(e),
        }
    }
}

impl From<i64> for ImplicitHydrogensAst {
    fn from(value: i64) -> Self {
        Self::Lit(value)
    }
}

impl Add for ImplicitHydrogensAst {
    type Output = ImplicitHydrogensAst;
    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Lit(a), Self::Lit(b)) => Self::Lit(a + b),
            _ => Self::Undetermined,
        }
    }
}

impl Sub for ImplicitHydrogensAst {
    type Output = ImplicitHydrogensAst;
    fn sub(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Lit(a), Self::Lit(b)) => Self::Lit(a - b),
            _ => Self::Undetermined,
        }
    }
}

impl Mul for ImplicitHydrogensAst {
    type Output = ImplicitHydrogensAst;
    fn mul(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Lit(a), Self::Lit(b)) => Self::Lit(a * b),
            _ => Self::Undetermined,
        }
    }
}

impl Div for ImplicitHydrogensAst {
    type Output = ImplicitHydrogensAst;
    fn div(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Lit(a), Self::Lit(b)) => Self::Lit(a / b),
            _ => Self::Undetermined,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::{AtomConstraint, AtomConstraintKind};
    use crate::atom_zeroed;

    #[rstest]
    fn test_atom_ast_from_element() {
        assert_eq!(
            AtomAst::from_element(Element::C),
            AtomAst {
                element: ElementAst::Lit(Element::C),
                ..Default::default()
            },
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_element_ast(AtomAst::default().with_element(ElementAst::Lit(Element::C)), AtomAst { element: ElementAst::Lit(Element::C), ..Default::default() })]
    #[case::with_element_primitive(AtomAst::default().with_element(Element::N), AtomAst { element: ElementAst::Lit(Element::N), ..Default::default() })]
    #[case::with_isotope_mass(AtomAst::default().with_isotope_mass(12_i64), AtomAst { isotope_mass: IsotopeAst::Lit(12), ..Default::default() })]
    #[case::with_charge(AtomAst::default().with_charge(1_i64), AtomAst { charge: ValueAst::Lit(1), ..Default::default() })]
    #[case::with_implicit_hydrogens(AtomAst::default().with_implicit_hydrogens(3_i64), AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Lit(3), ..Default::default() })]
    #[case::with_lone_pairs(AtomAst::default().with_lone_pairs(2_i64), AtomAst { lone_pairs: ValueAst::Lit(2), ..Default::default() })]
    #[case::with_spin_tuple(AtomAst::default().with_spin((0_u8, 1_u8)), AtomAst { spin: SpinStateAst::from((0_u8, 1_u8)), ..Default::default() })]
    #[case::with_constraint(AtomAst::default().with_constraint(AtomConstraint::valence(4_i64)),
        AtomAst { constraints: AtomConstraints::from(AtomConstraint::valence(4)),..Default::default() })]
    #[case::with_constraints_extends(
        AtomAst::default()
            .with_constraint(AtomConstraint::valence(4_i64))
            .with_constraints([AtomConstraint::donated_pairs(1_i64), AtomConstraint::ring_size(6_i64)]),
        AtomAst { constraints: AtomConstraints::from_iter([
            AtomConstraint::valence(4),
            AtomConstraint::donated_pairs(1),
            AtomConstraint::ring_size(6),
        ]), ..Default::default() },
    )]
    #[case::with_constraint_replaces_same_kind(
        AtomAst::default()
            .with_constraint(AtomConstraint::valence(3_i64))
            .with_constraint(AtomConstraint::valence(4_i64)),
        AtomAst { constraints: AtomConstraints::from(AtomConstraint::valence(4)), ..Default::default() },
    )]
    fn test_atom_ast_with_methods(#[case] actual: AtomAst, #[case] expected: AtomAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_element(AtomAst::from_element(Element::C).into_ground(),
        AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeAst::Natural, charge: ValueAst::Lit(0), implicit_hydrogens: ImplicitHydrogensAst::Lit(0),
        lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AtomConstraints::new() })]
    #[case::with_charge(AtomAst::from_element(Element::C).with_charge(1_i64).into_ground(),
        AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeAst::Natural, charge: ValueAst::Lit(1), implicit_hydrogens: ImplicitHydrogensAst::Lit(0),
        lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AtomConstraints::new() })]
    #[case::constraint(AtomAst::from_element(Element::C).with_constraint(AtomConstraint::valence(4_i64)).into_ground(),
        AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeAst::Natural, charge: ValueAst::Lit(0), implicit_hydrogens: ImplicitHydrogensAst::Lit(0),
        lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AtomConstraints::from(AtomConstraint::valence(4)) })]
    fn test_atom_ast_into_ground(#[case] actual: AtomAst, #[case] expected: AtomAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::element(AtomAst::from_element(Element::C).into_zeroed(), atom_zeroed!("C"))]
    #[case::element_charge(AtomAst::from_element(Element::C).with_charge(1_i64).into_zeroed(), atom_zeroed!("C #i= #c+"))]
    #[case::constraint(AtomAst::from_element(Element::C).with_constraint(AtomConstraint::valence(3_i64)).into_zeroed(),
        atom_zeroed!("C #v3"))]
    fn test_atom_ast_into_zeroed(#[case] actual: AtomAst, #[case] expected: AtomAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(AtomAst::default(), false)]
    #[case::all_ground(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ImplicitHydrogensAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, true)]
    #[case::element_undetermined(AtomAst { element: ElementAst::Undetermined, isotope_mass: IsotopeAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ImplicitHydrogensAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::isotope_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeAst::Undetermined, charge: ValueAst::Lit(0),
        implicit_hydrogens: ImplicitHydrogensAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::charge_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeAst::Lit(12), charge: ValueAst::Undetermined,
        implicit_hydrogens: ImplicitHydrogensAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::hydrogens_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ImplicitHydrogensAst::Undetermined, lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::lone_pairs_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ImplicitHydrogensAst::Lit(4), lone_pairs: ValueAst::Undetermined, spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::spin_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ImplicitHydrogensAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::default(),
        constraints: AtomConstraints::new() }, false)]
    fn test_atom_ast_is_ground(#[case] ast: AtomAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard_vs_ground(AtomAst::default(), AtomAst::from_element(Element::C), true)]
    #[case::same_element(AtomAst::from_element(Element::C), AtomAst::from_element(Element::C), true)]
    #[case::element_mismatch(AtomAst::from_element(Element::C), AtomAst::from_element(Element::N), false)]
    #[case::pattern_more_specific_than_target(AtomAst::from_element(Element::C), AtomAst::default(), false)]
    #[case::charge_mismatch(AtomAst::from_element(Element::C).with_charge(1_i64), AtomAst::from_element(Element::C).with_charge(0_i64), false)]
    #[case::charge_wildcard_pattern(AtomAst::from_element(Element::C), AtomAst::from_element(Element::C).with_charge(1_i64), true)]
    #[case::isotope_mismatch(AtomAst::from_element(Element::C).with_isotope_mass(12_i64), AtomAst::from_element(Element::C).with_isotope_mass(13_i64), false)]
    #[case::hydrogens_mismatch(AtomAst::from_element(Element::C).with_implicit_hydrogens(3_i64), AtomAst::from_element(Element::C).with_implicit_hydrogens(4_i64), false)]
    #[case::lone_pairs_mismatch(AtomAst::from_element(Element::C).with_lone_pairs(1_i64), AtomAst::from_element(Element::C).with_lone_pairs(2_i64), false)]
    #[case::spin_mismatch(AtomAst::from_element(Element::C).with_spin((2_u8, 3_u8)), AtomAst::from_element(Element::C).with_spin((0_u8, 1_u8)), false)]
    fn test_atom_ast_matches(
        #[case] pattern: AtomAst,
        #[case] target: AtomAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_atom_ast_simplify_values() {
        let mut atom = AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Expr(Expr::Lit(12)),
            charge: ValueAst::Expr(Expr::Lit(1)),
            implicit_hydrogens: ImplicitHydrogensAst::Expr(Expr::Lit(3)),
            lone_pairs: ValueAst::Expr(Expr::Neg(Box::new(Expr::Lit(2)))),
            spin: SpinStateAst {
                unpaired: ValueAst::Expr(Expr::Lit(0)),
                multiplicity: ValueAst::Expr(Expr::Lit(1)),
            },
            constraints: AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Expr(
                Expr::Lit(4),
            ))]),
        };
        atom.simplify_values();
        assert_eq!(atom.isotope_mass, IsotopeAst::Lit(12));
        assert_eq!(atom.charge, ValueAst::Lit(1));
        assert_eq!(atom.implicit_hydrogens, ImplicitHydrogensAst::Lit(3));
        assert_eq!(atom.lone_pairs, ValueAst::Lit(-2));
        assert_eq!(atom.spin.unpaired, ValueAst::Lit(0));
        assert_eq!(atom.spin.multiplicity, ValueAst::Lit(1));
        assert_eq!(
            atom.constraints.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::valence(4)),
        );
    }

    #[rstest]
    #[case::carbon(Element::C, ElementAst::Lit(Element::C))]
    #[case::nitrogen(Element::N, ElementAst::Lit(Element::N))]
    fn test_element_ast_from(#[case] element: Element, #[case] expected: ElementAst) {
        assert_eq!(ElementAst::from(element), expected);
    }

    #[rstest]
    #[case::lit_carbon(ElementAst::Lit(Element::C), Some(Element::C))]
    #[case::lit_nitrogen(ElementAst::Lit(Element::N), Some(Element::N))]
    #[case::wildcard(ElementAst::Undetermined, None)]
    #[case::set(ElementAst::Set(vec![Element::C, Element::N]), None)]
    #[case::bind(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, None)]
    #[case::reference(ElementAst::Ref("e".into()), None)]
    fn test_element_ast_literal_and_is_ground(
        #[case] ast: ElementAst,
        #[case] expected: Option<Element>,
    ) {
        assert_eq!(ast.literal(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rstest]
    #[case::lit(ElementAst::Lit(Element::C), false)]
    #[case::wildcard(ElementAst::Undetermined, true)]
    #[case::set(ElementAst::Set(vec![Element::C, Element::N]), false)]
    #[case::bind(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, false)]
    #[case::reference(ElementAst::Ref("e".into()), false)]
    fn test_element_ast_is_undetermined(#[case] ast: ElementAst, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_lit(ElementAst::Undetermined, ElementAst::Lit(Element::C), true)]
    #[case::undetermined_undetermined(ElementAst::Undetermined, ElementAst::Undetermined, true)]
    #[case::undetermined_set(ElementAst::Undetermined, ElementAst::Set(vec![Element::C, Element::N]), true)]
    #[case::lit_undetermined(ElementAst::Lit(Element::C), ElementAst::Undetermined, false)]
    #[case::set_undetermined(ElementAst::Set(vec![Element::C]), ElementAst::Undetermined, false)]
    #[case::lit_lit_match(ElementAst::Lit(Element::C), ElementAst::Lit(Element::C), true)]
    #[case::lit_lit_mismatch(ElementAst::Lit(Element::C), ElementAst::Lit(Element::N), false)]
    #[case::lit_singleton_set(ElementAst::Lit(Element::C), ElementAst::Set(vec![Element::C]), true)]
    #[case::lit_multi_set(ElementAst::Lit(Element::C), ElementAst::Set(vec![Element::C, Element::N]), false)]
    #[case::set_lit_in(ElementAst::Set(vec![Element::C, Element::N]), ElementAst::Lit(Element::N), true)]
    #[case::set_lit_out(ElementAst::Set(vec![Element::C, Element::N]), ElementAst::Lit(Element::O), false)]
    #[case::set_set_subset(ElementAst::Set(vec![Element::C, Element::N, Element::O]), ElementAst::Set(vec![Element::C, Element::N]), true)]
    #[case::set_set_equal(ElementAst::Set(vec![Element::C, Element::N]), ElementAst::Set(vec![Element::C, Element::N]), true)]
    #[case::set_set_superset(ElementAst::Set(vec![Element::C]), ElementAst::Set(vec![Element::C, Element::N]), false)]
    #[case::bind_lit_match(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, ElementAst::Lit(Element::C), true)]
    #[case::bind_lit_mismatch(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, ElementAst::Lit(Element::N), false)]
    #[case::bind_set_subset(ElementAst::Bind { id: "e".into(), set: vec![Element::C, Element::N] }, ElementAst::Set(vec![Element::C]), true)]
    #[case::set_bind_subset(ElementAst::Set(vec![Element::C, Element::N]), ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, true)]
    #[case::bind_bind_subset(ElementAst::Bind { id: "p".into(), set: vec![Element::C, Element::N] }, ElementAst::Bind { id: "t".into(), set: vec![Element::N] }, true)]
    #[case::bind_bind_superset(ElementAst::Bind { id: "p".into(), set: vec![Element::C] }, ElementAst::Bind { id: "t".into(), set: vec![Element::C, Element::N] }, false)]
    #[case::undetermined_bind(ElementAst::Undetermined, ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, true)]
    #[case::bind_undetermined(ElementAst::Bind { id: "e".into(), set: vec![Element::C] }, ElementAst::Undetermined, false)]
    #[case::ref_lit(ElementAst::Ref("e".into()), ElementAst::Lit(Element::C), false)]
    #[case::lit_ref(ElementAst::Lit(Element::C), ElementAst::Ref("e".into()), false)]
    #[case::ref_set(ElementAst::Ref("e".into()), ElementAst::Set(vec![Element::C]), false)]
    #[case::set_ref(ElementAst::Set(vec![Element::C]), ElementAst::Ref("e".into()), false)]
    #[case::ref_bind(ElementAst::Ref("e".into()), ElementAst::Bind { id: "f".into(), set: vec![Element::C] }, false)]
    #[case::ref_ref(ElementAst::Ref("e".into()), ElementAst::Ref("f".into()), false)]
    #[case::ref_undetermined(ElementAst::Ref("e".into()), ElementAst::Undetermined, false)]
    fn test_element_ast_matches(
        #[case] pattern: ElementAst,
        #[case] target: ElementAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::from_lit(IsotopeAst::from(ValueAst::Lit(13)), IsotopeAst::Lit(13))]
    #[case::from_undetermined(IsotopeAst::from(ValueAst::Undetermined), IsotopeAst::Undetermined)]
    fn test_isotope_ast_from_value(#[case] actual: IsotopeAst, #[case] expected: IsotopeAst) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::positive(IsotopeAst::from(13_i64), IsotopeAst::Lit(13))]
    #[case::zero(IsotopeAst::from(0_i64), IsotopeAst::Lit(0))]
    fn test_isotope_ast_from_i64(#[case] actual: IsotopeAst, #[case] expected: IsotopeAst) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::natural(IsotopeAst::Natural, false)]
    #[case::lit(IsotopeAst::Lit(12), false)]
    #[case::undetermined(IsotopeAst::Undetermined, true)]
    #[case::lit_set(IsotopeAst::LitSet(vec![12, 13]), false)]
    #[case::expr(IsotopeAst::Expr(Expr::Lit(12)), false)]
    fn test_isotope_ast_is_undetermined(#[case] ast: IsotopeAst, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rstest]
    #[case::natural(IsotopeAst::Natural, Some(0))]
    #[case::lit(IsotopeAst::Lit(12), Some(12))]
    #[case::lit_zero(IsotopeAst::Lit(0), Some(0))]
    #[case::wildcard(IsotopeAst::Undetermined, None)]
    #[case::set_singleton(IsotopeAst::LitSet(vec![14]), Some(14))]
    #[case::set_multi(IsotopeAst::LitSet(vec![12, 13]), None)]
    #[case::expr_lit(IsotopeAst::Expr(Expr::Lit(15)), Some(15))]
    #[case::expr_var(IsotopeAst::Expr(Expr::Var("x".into())), None)]
    fn test_isotope_ast_literal_and_is_ground(
        #[case] ast: IsotopeAst,
        #[case] expected: Option<u32>,
    ) {
        assert_eq!(ast.literal(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_natural(IsotopeAst::Undetermined, IsotopeAst::Natural, true)]
    #[case::undetermined_value(IsotopeAst::Undetermined, IsotopeAst::Lit(12), true)]
    #[case::undetermined_undetermined(IsotopeAst::Undetermined, IsotopeAst::Undetermined, true)]
    #[case::natural_undetermined(IsotopeAst::Natural, IsotopeAst::Undetermined, false)]
    #[case::value_undetermined(IsotopeAst::Lit(12), IsotopeAst::Undetermined, false)]
    #[case::natural_natural(IsotopeAst::Natural, IsotopeAst::Natural, true)]
    #[case::natural_value(IsotopeAst::Natural, IsotopeAst::Lit(12), false)]
    #[case::value_natural(IsotopeAst::Lit(12), IsotopeAst::Natural, false)]
    #[case::value_lit_match(IsotopeAst::Lit(12), IsotopeAst::Lit(12), true)]
    #[case::value_lit_mismatch(IsotopeAst::Lit(12), IsotopeAst::Lit(13), false)]
    #[case::value_wildcard_lit(IsotopeAst::Undetermined, IsotopeAst::Lit(12), true)]
    #[case::value_set_lit_in(IsotopeAst::LitSet(vec![12, 13]), IsotopeAst::Lit(13), true)]
    #[case::value_set_lit_out(IsotopeAst::LitSet(vec![12, 13]), IsotopeAst::Lit(14), false)]
    #[case::value_set_set_subset(IsotopeAst::LitSet(vec![12, 13, 14]), IsotopeAst::LitSet(vec![12, 13]), true)]
    #[case::value_set_set_superset(IsotopeAst::LitSet(vec![12]), IsotopeAst::LitSet(vec![12, 13]), false)]
    #[case::value_lit_wildcard(IsotopeAst::Lit(12), IsotopeAst::Undetermined, false)]
    fn test_isotope_ast_matches(
        #[case] pattern: IsotopeAst,
        #[case] target: IsotopeAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::expr_lit(IsotopeAst::Expr(Expr::Lit(13)), IsotopeAst::Lit(13))]
    #[case::expr_neg_lit(IsotopeAst::Expr(Expr::Neg(Box::new(Expr::Lit(12)))), IsotopeAst::Lit(-12))]
    fn test_isotope_ast_simplify(#[case] input: IsotopeAst, #[case] expected: IsotopeAst) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::natural(IsotopeAst::Natural)]
    #[case::lit(IsotopeAst::Lit(12))]
    #[case::undetermined(IsotopeAst::Undetermined)]
    #[case::lit_set(IsotopeAst::LitSet(vec![12, 13]))]
    #[case::expr_var(IsotopeAst::Expr(Expr::Var("x".into())))]
    fn test_isotope_ast_simplify_identity(#[case] input: IsotopeAst) {
        assert_eq!(input.clone().simplify(), input);
    }

    #[rstest]
    #[case::from_lit_set(ImplicitHydrogensAst::from(ValueAst::LitSet(vec![0, 1])), ImplicitHydrogensAst::LitSet(vec![0, 1]))]
    fn test_implicit_hydrogens_ast_from_value(
        #[case] actual: ImplicitHydrogensAst,
        #[case] expected: ImplicitHydrogensAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::positive(ImplicitHydrogensAst::from(3_i64), ImplicitHydrogensAst::Lit(3))]
    #[case::zero(ImplicitHydrogensAst::from(0_i64), ImplicitHydrogensAst::Lit(0))]
    fn test_implicit_hydrogens_ast_from_i64(
        #[case] actual: ImplicitHydrogensAst,
        #[case] expected: ImplicitHydrogensAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::normal(ImplicitHydrogensAst::Normal, false)]
    #[case::lit(ImplicitHydrogensAst::Lit(2), false)]
    #[case::undetermined(ImplicitHydrogensAst::Undetermined, true)]
    #[case::lit_set(ImplicitHydrogensAst::LitSet(vec![1, 2]), false)]
    #[case::expr(ImplicitHydrogensAst::Expr(Expr::Lit(2)), false)]
    fn test_implicit_hydrogens_ast_is_undetermined(
        #[case] ast: ImplicitHydrogensAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rstest]
    #[case::normal(ImplicitHydrogensAst::Normal, None)]
    #[case::lit(ImplicitHydrogensAst::Lit(2), Some(2))]
    #[case::lit_zero(ImplicitHydrogensAst::Lit(0), Some(0))]
    #[case::wildcard(ImplicitHydrogensAst::Undetermined, None)]
    #[case::set_singleton(ImplicitHydrogensAst::LitSet(vec![3]), Some(3))]
    #[case::set_multi(ImplicitHydrogensAst::LitSet(vec![1, 2]), None)]
    #[case::expr_lit(ImplicitHydrogensAst::Expr(Expr::Lit(4)), Some(4))]
    #[case::expr_var(ImplicitHydrogensAst::Expr(Expr::Var("x".into())), None)]
    fn test_implicit_hydrogens_ast_literal_and_is_ground(
        #[case] ast: ImplicitHydrogensAst,
        #[case] expected: Option<i64>,
    ) {
        assert_eq!(ast.literal(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_normal(ImplicitHydrogensAst::Undetermined, ImplicitHydrogensAst::Normal, true)]
    #[case::undetermined_value(ImplicitHydrogensAst::Undetermined, ImplicitHydrogensAst::Lit(3), true)]
    #[case::normal_undetermined(ImplicitHydrogensAst::Normal, ImplicitHydrogensAst::Undetermined, false)]
    #[case::normal_normal(ImplicitHydrogensAst::Normal, ImplicitHydrogensAst::Normal, true)]
    #[case::normal_value(ImplicitHydrogensAst::Normal, ImplicitHydrogensAst::Lit(0), false)]
    #[case::value_normal(ImplicitHydrogensAst::Lit(0), ImplicitHydrogensAst::Normal, false)]
    #[case::value_lit_match(ImplicitHydrogensAst::Lit(2), ImplicitHydrogensAst::Lit(2), true)]
    #[case::value_lit_mismatch(ImplicitHydrogensAst::Lit(2), ImplicitHydrogensAst::Lit(3), false)]
    #[case::value_wildcard(ImplicitHydrogensAst::Undetermined, ImplicitHydrogensAst::Lit(2), true)]
    #[case::value_set_subset(ImplicitHydrogensAst::LitSet(vec![1, 2]), ImplicitHydrogensAst::LitSet(vec![1]), true)]
    fn test_implicit_hydrogens_ast_matches(
        #[case] pattern: ImplicitHydrogensAst,
        #[case] target: ImplicitHydrogensAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::expr_lit(ImplicitHydrogensAst::Expr(Expr::Lit(4)), ImplicitHydrogensAst::Lit(4))]
    #[case::expr_neg_lit(ImplicitHydrogensAst::Expr(Expr::Neg(Box::new(Expr::Lit(3)))), ImplicitHydrogensAst::Lit(-3))]
    fn test_implicit_hydrogens_ast_simplify(
        #[case] input: ImplicitHydrogensAst,
        #[case] expected: ImplicitHydrogensAst,
    ) {
        assert_eq!(input.simplify(), expected);
    }

    #[rstest]
    #[case::normal(ImplicitHydrogensAst::Normal)]
    #[case::lit(ImplicitHydrogensAst::Lit(2))]
    #[case::undetermined(ImplicitHydrogensAst::Undetermined)]
    #[case::lit_set(ImplicitHydrogensAst::LitSet(vec![1, 2]))]
    #[case::expr_var(ImplicitHydrogensAst::Expr(Expr::Var("h".into())))]
    fn test_implicit_hydrogens_ast_simplify_identity(#[case] input: ImplicitHydrogensAst) {
        assert_eq!(input.clone().simplify(), input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_lit(
        ImplicitHydrogensAst::Lit(2), ImplicitHydrogensAst::Lit(1),
        ImplicitHydrogensAst::Lit(3),
    )]
    #[case::normal_collapses(
        ImplicitHydrogensAst::Normal, ImplicitHydrogensAst::Lit(1),
        ImplicitHydrogensAst::Undetermined,
    )]
    #[case::undetermined_collapses(
        ImplicitHydrogensAst::Lit(2), ImplicitHydrogensAst::Undetermined,
        ImplicitHydrogensAst::Undetermined,
    )]
    fn test_implicit_hydrogens_ast_add(
        #[case] lhs: ImplicitHydrogensAst,
        #[case] rhs: ImplicitHydrogensAst,
        #[case] expected: ImplicitHydrogensAst,
    ) {
        assert_eq!(lhs + rhs, expected);
    }

    #[rstest]
    #[case::lit_lit(
        ImplicitHydrogensAst::Lit(5), ImplicitHydrogensAst::Lit(2),
        ImplicitHydrogensAst::Lit(3),
    )]
    #[case::normal_collapses(
        ImplicitHydrogensAst::Normal, ImplicitHydrogensAst::Lit(1),
        ImplicitHydrogensAst::Undetermined,
    )]
    fn test_implicit_hydrogens_ast_sub(
        #[case] lhs: ImplicitHydrogensAst,
        #[case] rhs: ImplicitHydrogensAst,
        #[case] expected: ImplicitHydrogensAst,
    ) {
        assert_eq!(lhs - rhs, expected);
    }

    #[rstest]
    #[case::lit_lit(
        ImplicitHydrogensAst::Lit(3), ImplicitHydrogensAst::Lit(4),
        ImplicitHydrogensAst::Lit(12),
    )]
    #[case::normal_collapses(
        ImplicitHydrogensAst::Lit(3), ImplicitHydrogensAst::Normal,
        ImplicitHydrogensAst::Undetermined,
    )]
    fn test_implicit_hydrogens_ast_mul(
        #[case] lhs: ImplicitHydrogensAst,
        #[case] rhs: ImplicitHydrogensAst,
        #[case] expected: ImplicitHydrogensAst,
    ) {
        assert_eq!(lhs * rhs, expected);
    }

    #[rstest]
    #[case::lit_lit(
        ImplicitHydrogensAst::Lit(10), ImplicitHydrogensAst::Lit(2),
        ImplicitHydrogensAst::Lit(5),
    )]
    #[case::normal_collapses(
        ImplicitHydrogensAst::Normal, ImplicitHydrogensAst::Lit(2),
        ImplicitHydrogensAst::Undetermined,
    )]
    fn test_implicit_hydrogens_ast_div(
        #[case] lhs: ImplicitHydrogensAst,
        #[case] rhs: ImplicitHydrogensAst,
        #[case] expected: ImplicitHydrogensAst,
    ) {
        assert_eq!(lhs / rhs, expected);
    }

    #[rstest]
    #[should_panic]
    fn test_implicit_hydrogens_ast_div_by_zero_panics() {
        let _ = ImplicitHydrogensAst::Lit(5) / ImplicitHydrogensAst::Lit(0);
    }
}
