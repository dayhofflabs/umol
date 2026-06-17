//! Atom-level AST fragments shared across crates.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::mem;

use umol_ast_macros::Lattice;
use umol_shared::element::Element;

use super::constraint::{
    AromaticValenceAst, AtomConstraint, AtomConstraintKind, AtomConstraints, MulticenterValenceAst,
};
use super::error::Contradiction;
use super::operators::MemOp;
use super::spin::SpinStateAst;
use super::stereo::StereoConfigurationAst;
use super::traits::{AsLit, Canonicalize, Lattice};
use super::value::ValueAst;

/// Atom AST: structural representation of an atom plus the atom-level
/// constraints (valence, degree, ring membership, etc.) that pattern
/// against the surrounding topology.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Lattice)]
pub struct AtomAst {
    pub element: ElementAst,
    pub isotope_mass: IsotopeMassAst,
    pub charge: ValueAst,
    pub implicit_hydrogens: ValueAst,
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

    pub fn with_isotope_mass(mut self, mass: impl Into<IsotopeMassAst>) -> Self {
        self.isotope_mass = mass.into();
        self
    }

    pub fn with_charge(mut self, charge: impl Into<ValueAst>) -> Self {
        self.charge = charge.into();
        self
    }

    pub fn with_implicit_hydrogens(mut self, hydrogens: impl Into<ValueAst>) -> Self {
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

    /// Fill `Undetermined` value-bearing struct fields with defaults: isotope→
    /// Natural; charge / implicit hydrogens / lone pairs → 0; spin → 0 unpaired
    /// and, for the (possibly already-set) unpaired count, the maximal
    /// multiplicity `unpaired + 1` (so a fully unset spin becomes the
    /// closed-shell singlet). Existing literal or expression values and all
    /// constraints are preserved. The result is ground iff `element` is
    /// already ground.
    pub fn into_ground(mut self) -> Self {
        if self.isotope_mass.is_undetermined() {
            self.isotope_mass = IsotopeMassAst::Natural;
        }
        if self.charge.is_undetermined() {
            self.charge = ValueAst::Lit(0);
        }
        if self.implicit_hydrogens.is_undetermined() {
            self.implicit_hydrogens = ValueAst::Lit(0);
        }
        if self.lone_pairs.is_undetermined() {
            self.lone_pairs = ValueAst::Lit(0);
        }
        if self.spin.unpaired.is_undetermined() {
            self.spin.unpaired = ValueAst::Lit(0);
        }
        if self.spin.multiplicity.is_undetermined() {
            let unpaired = self.spin.unpaired.as_lit().unwrap_or(0);
            self.spin.multiplicity = ValueAst::Lit(unpaired + 1);
        }
        self
    }

    /// `into_ground()` plus chemistry-default constraints for an isolated,
    /// non-bonded, non-aromatic, non-multicenter atom: `Valence(0)`,
    /// `DonatedPairs(0)`, `AcceptedPairs(0)`, `MulticenterValence(NotMulticenter)`,
    /// `AromaticValence(NotAromatic)`, `TetrahedralStereo(NotStereo)`. Each is
    /// added only if the corresponding constraint kind is not already present;
    /// existing entries are preserved. Matches the `atom_zeroed!` macro semantics.
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
        if !self
            .constraints
            .contains(AtomConstraintKind::TetrahedralStereo)
        {
            self.constraints.add(AtomConstraint::TetrahedralStereo(
                StereoConfigurationAst::NotStereo,
            ));
        }
        self
    }

    /// Simplify every value-bearing field in place: `isotope_mass`,
    /// `charge`, `implicit_hydrogens`, `lone_pairs`, both `spin` slots,
    /// and each constraint. `element` has no value to simplify.
    pub fn simplify_values(&mut self) {
        self.charge = mem::take(&mut self.charge).simplify();
        self.implicit_hydrogens = mem::take(&mut self.implicit_hydrogens).simplify();
        self.lone_pairs = mem::take(&mut self.lone_pairs).simplify();
        self.spin.simplify_values();
        self.constraints.simplify_each();
    }
}

/// Element expression: undetermined, a single element, a finite element set, a
/// complement set (`!{…}`), or a variable (free `?x`, or membership-restricted
/// `?x :: {…}` / `?x :: !{…}`). Sets are cardinality-canonical and
/// universe-relative: a semantic set of more than `⌊118/2⌋` elements is stored
/// as its complement `NotSet(U∖S)`.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ElementAst {
    #[default]
    Undetermined,
    Lit(Element),
    LitSet(Box<BTreeSet<Element>>),
    NotSet(Box<BTreeSet<Element>>),
    #[allow(clippy::type_complexity)]
    Var(Box<(String, Option<(MemOp, BTreeSet<Element>)>)>),
}

impl ElementAst {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn lit(element: Element) -> Self {
        Self::Lit(element)
    }

    pub fn lit_set(elements: impl IntoIterator<Item = Element>) -> Self {
        Self::LitSet(Box::new(elements.into_iter().collect()))
    }

    pub fn not(element: Element) -> Self {
        Self::NotSet(Box::new(BTreeSet::from([element])))
    }

    pub fn not_set(elements: impl IntoIterator<Item = Element>) -> Self {
        Self::NotSet(Box::new(elements.into_iter().collect()))
    }

    /// A free element variable `?name`.
    pub fn var(name: impl Into<String>) -> Self {
        Self::Var(Box::new((name.into(), None)))
    }

    /// A variable restricted to membership `?name :: {…}`.
    pub fn var_in(name: impl Into<String>, elements: impl IntoIterator<Item = Element>) -> Self {
        Self::Var(Box::new((
            name.into(),
            Some((MemOp::In, elements.into_iter().collect())),
        )))
    }

    /// A variable restricted to non-membership `?name :: !{…}`.
    pub fn var_not_in(name: impl Into<String>, elements: impl IntoIterator<Item = Element>) -> Self {
        Self::Var(Box::new((
            name.into(),
            Some((MemOp::NotIn, elements.into_iter().collect())),
        )))
    }
}

impl From<Element> for ElementAst {
    fn from(element: Element) -> Self {
        Self::Lit(element)
    }
}

impl Canonicalize for ElementAst {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            ElementAst::Undetermined => ElementAst::Undetermined,
            ElementAst::Lit(e) => ElementAst::Lit(e),
            ElementAst::LitSet(s) => canon_set(*s, false)?,
            ElementAst::NotSet(s) => canon_set(*s, true)?,
            ElementAst::Var(v) => {
                let (name, domain) = *v;
                let domain = match domain {
                    None => None,
                    Some((op, set)) => canon_var_domain(op, set)?,
                };
                ElementAst::Var(Box::new((name, domain)))
            }
        })
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            ElementAst::Undetermined | ElementAst::Lit(_) => Ok(Cow::Borrowed(self)),
            _ => Ok(Cow::Owned(self.clone().canonicalize()?)),
        }
    }
}

fn universe_size() -> usize {
    Element::all().len()
}

fn complement(s: &BTreeSet<Element>) -> BTreeSet<Element> {
    Element::all()
        .iter()
        .copied()
        .filter(|e| !s.contains(e))
        .collect()
}

/// Canonicalize a set given its polarity (`negated` = the set is a complement,
/// denoting `U∖s`). The semantic set is stored on its smaller side: `≤ ⌊|U|/2⌋`
/// → positive (`Lit`/`LitSet`), else `NotSet` of its complement (tiebreak
/// positive). Empty semantic set → `Err`; full → `Undetermined`.
fn canon_set(s: BTreeSet<Element>, negated: bool) -> Result<ElementAst, Contradiction> {
    let u = universe_size();
    let semantic_len = if negated { u - s.len() } else { s.len() };
    if semantic_len == 0 {
        return Err(Contradiction);
    }
    if semantic_len == u {
        return Ok(ElementAst::Undetermined);
    }
    if semantic_len <= u / 2 {
        let positive = if negated { complement(&s) } else { s };
        Ok(if positive.len() == 1 {
            ElementAst::Lit(*positive.iter().next().unwrap())
        } else {
            ElementAst::LitSet(Box::new(positive))
        })
    } else {
        let excluded = if negated { s } else { complement(&s) };
        Ok(ElementAst::NotSet(Box::new(excluded)))
    }
}

/// Canonicalize a `Var` membership domain by the same cardinality polarity.
/// Vacuous (`In U` / `NotIn ∅`) → `None` (free); empty admissible set → `Err`.
fn canon_var_domain(
    op: MemOp,
    set: BTreeSet<Element>,
) -> Result<Option<(MemOp, BTreeSet<Element>)>, Contradiction> {
    let u = universe_size();
    let admissible_len = match op {
        MemOp::In => set.len(),
        MemOp::NotIn => u - set.len(),
    };
    if admissible_len == 0 {
        return Err(Contradiction);
    }
    if admissible_len == u {
        return Ok(None);
    }
    Ok(Some(if admissible_len <= u / 2 {
        let admissible = match op {
            MemOp::In => set,
            MemOp::NotIn => complement(&set),
        };
        (MemOp::In, admissible)
    } else {
        let excluded = match op {
            MemOp::In => complement(&set),
            MemOp::NotIn => set,
        };
        (MemOp::NotIn, excluded)
    }))
}

impl AsLit for ElementAst {
    type Lit = Element;

    /// The single element this denotes, only when it is a literal.
    /// Non-canonicalizing (mirrors `ValueAst::as_lit`).
    #[inline]
    fn as_lit(&self) -> Option<Element> {
        match self {
            Self::Lit(e) => Some(*e),
            _ => None,
        }
    }
}

/// The concrete forms as `(set, negated)`: `Lit`/`LitSet` positive, `NotSet` the
/// complement. `Undetermined`/`Var` have no finite-set view.
fn element_set_view(e: &ElementAst) -> Option<(BTreeSet<Element>, bool)> {
    match e {
        ElementAst::Lit(x) => Some((BTreeSet::from([*x]), false)),
        ElementAst::LitSet(s) => Some(((**s).clone(), false)),
        ElementAst::NotSet(s) => Some(((**s).clone(), true)),
        ElementAst::Undetermined | ElementAst::Var(_) => None,
    }
}

fn intersect(s: &BTreeSet<Element>, t: &BTreeSet<Element>) -> BTreeSet<Element> {
    s.intersection(t).copied().collect()
}

fn union(s: &BTreeSet<Element>, t: &BTreeSet<Element>) -> BTreeSet<Element> {
    s.union(t).copied().collect()
}

fn difference(s: &BTreeSet<Element>, t: &BTreeSet<Element>) -> BTreeSet<Element> {
    s.difference(t).copied().collect()
}

impl Lattice for ElementAst {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    /// Bottom — resolves to a single element. Literal-only, non-canonicalizing
    /// (aligned with `as_lit`).
    #[inline]
    fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    /// Greatest lower bound (set intersection), canonicalizing operands and
    /// result. `Var` meets only an equal `Var`; `Var` vs concrete → `None`.
    fn meet(&self, other: &Self) -> Option<Self> {
        let a = self.canonical().ok()?;
        let b = other.canonical().ok()?;
        use ElementAst::*;
        match (a.as_ref(), b.as_ref()) {
            (Undetermined, _) => Some(b.as_ref().clone()),
            (_, Undetermined) => Some(a.as_ref().clone()),
            (Var(x), Var(y)) => (x == y).then(|| a.as_ref().clone()),
            (Var(_), _) | (_, Var(_)) => None,
            _ => {
                let (sa, na) = element_set_view(a.as_ref()).unwrap();
                let (sb, nb) = element_set_view(b.as_ref()).unwrap();
                let (set, negated) = match (na, nb) {
                    (false, false) => (intersect(&sa, &sb), false),
                    (false, true) => (difference(&sa, &sb), false),
                    (true, false) => (difference(&sb, &sa), false),
                    (true, true) => (union(&sa, &sb), true),
                };
                canon_set(set, negated).ok()
            }
        }
    }

    /// Least upper bound (set union), canonicalizing operands and result.
    fn join(&self, other: &Self) -> Self {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        use ElementAst::*;
        match (a.as_ref(), b.as_ref()) {
            (Undetermined, _) | (_, Undetermined) => Undetermined,
            (Var(x), Var(y)) if x == y => a.as_ref().clone(),
            (Var(_), _) | (_, Var(_)) => Undetermined,
            _ => {
                let (sa, na) = element_set_view(a.as_ref()).unwrap();
                let (sb, nb) = element_set_view(b.as_ref()).unwrap();
                let (set, negated) = match (na, nb) {
                    (false, false) => (union(&sa, &sb), false),
                    (false, true) => (difference(&sb, &sa), true),
                    (true, false) => (difference(&sa, &sb), true),
                    (true, true) => (intersect(&sa, &sb), true),
                };
                canon_set(set, negated).unwrap_or(Undetermined)
            }
        }
    }

    /// `target` refines `self`: `self.meet(target) == canonical(target)`.
    fn matches(&self, target: &Self) -> bool {
        match (self.meet(target), target.canonical()) {
            (Some(meet), Ok(target)) => meet == *target,
            _ => false,
        }
    }
}

/// Isotope-mass expression: undetermined, the natural isotopic mixture
/// (`#i=`), a single mass number, a finite mass set, or a variable (free
/// `?x`, or membership-restricted `?x :: {…}`). Positive-only — no negation
/// and no complement (the mass domain is open). `Natural` is a distinct
/// ground, disjoint from every specific mass.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IsotopeMassAst {
    #[default]
    Undetermined,
    Natural,
    Lit(u32),
    LitSet(Box<BTreeSet<u32>>),
    Var(Box<(String, Option<BTreeSet<u32>>)>),
}

impl IsotopeMassAst {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn natural() -> Self {
        Self::Natural
    }

    pub fn lit(mass: u32) -> Self {
        Self::Lit(mass)
    }

    pub fn lit_set(masses: impl IntoIterator<Item = u32>) -> Self {
        Self::LitSet(Box::new(masses.into_iter().collect()))
    }

    /// A free isotope variable `?name`.
    pub fn var(name: impl Into<String>) -> Self {
        Self::Var(Box::new((name.into(), None)))
    }

    /// A variable restricted to a finite mass set `?name :: {…}`.
    pub fn var_in(name: impl Into<String>, masses: impl IntoIterator<Item = u32>) -> Self {
        Self::Var(Box::new((name.into(), Some(masses.into_iter().collect()))))
    }
}

impl From<u32> for IsotopeMassAst {
    fn from(mass: u32) -> Self {
        Self::Lit(mass)
    }
}

impl Canonicalize for IsotopeMassAst {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            IsotopeMassAst::Undetermined => IsotopeMassAst::Undetermined,
            IsotopeMassAst::Natural => IsotopeMassAst::Natural,
            IsotopeMassAst::Lit(n) => IsotopeMassAst::Lit(n),
            IsotopeMassAst::LitSet(s) => canon_mass_set(*s)?,
            IsotopeMassAst::Var(v) => {
                let (name, domain) = *v;
                let domain = match domain {
                    None => None,
                    Some(set) => {
                        if set.is_empty() {
                            return Err(Contradiction);
                        }
                        Some(set)
                    }
                };
                IsotopeMassAst::Var(Box::new((name, domain)))
            }
        })
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            IsotopeMassAst::Undetermined | IsotopeMassAst::Natural | IsotopeMassAst::Lit(_) => {
                Ok(Cow::Borrowed(self))
            }
            _ => Ok(Cow::Owned(self.clone().canonicalize()?)),
        }
    }
}

/// Canonicalize a mass set: empty → `Err`; singleton → `Lit`; else `LitSet`.
fn canon_mass_set(s: BTreeSet<u32>) -> Result<IsotopeMassAst, Contradiction> {
    match s.len() {
        0 => Err(Contradiction),
        1 => Ok(IsotopeMassAst::Lit(*s.iter().next().unwrap())),
        _ => Ok(IsotopeMassAst::LitSet(Box::new(s))),
    }
}

impl AsLit for IsotopeMassAst {
    type Lit = u32;

    /// The single mass number this denotes, only when it is a literal.
    /// `Natural` yields `None` (it commits to no specific mass).
    /// Non-canonicalizing (mirrors `ElementAst::as_lit`).
    #[inline]
    fn as_lit(&self) -> Option<u32> {
        match self {
            Self::Lit(n) => Some(*n),
            _ => None,
        }
    }
}

/// The positive mass set a concrete form denotes. `Undetermined`, `Natural`,
/// and `Var` have no mass-set view.
fn mass_set_view(iso: &IsotopeMassAst) -> Option<BTreeSet<u32>> {
    match iso {
        IsotopeMassAst::Lit(n) => Some(BTreeSet::from([*n])),
        IsotopeMassAst::LitSet(s) => Some((**s).clone()),
        IsotopeMassAst::Undetermined | IsotopeMassAst::Natural | IsotopeMassAst::Var(_) => None,
    }
}

impl Lattice for IsotopeMassAst {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    /// `Natural` and a single mass are ground; sets, `Var`, and
    /// `Undetermined` are not.
    #[inline]
    fn is_ground(&self) -> bool {
        matches!(self, Self::Natural | Self::Lit(_))
    }

    /// Greatest lower bound. `Undetermined` is top; `Natural` is an isolated
    /// ground (meets only itself); mass sets meet by intersection (empty →
    /// `None`). `Var` meets only an equal `Var`; `Var` vs concrete → `None`.
    fn meet(&self, other: &Self) -> Option<Self> {
        let a = self.canonical().ok()?;
        let b = other.canonical().ok()?;
        use IsotopeMassAst::*;
        match (a.as_ref(), b.as_ref()) {
            (Undetermined, _) => Some(b.as_ref().clone()),
            (_, Undetermined) => Some(a.as_ref().clone()),
            (Natural, Natural) => Some(Natural),
            (Natural, _) | (_, Natural) => None,
            (Var(x), Var(y)) => (x == y).then(|| a.as_ref().clone()),
            (Var(_), _) | (_, Var(_)) => None,
            _ => {
                let sa = mass_set_view(a.as_ref()).unwrap();
                let sb = mass_set_view(b.as_ref()).unwrap();
                let intersection: BTreeSet<u32> = sa.intersection(&sb).copied().collect();
                canon_mass_set(intersection).ok()
            }
        }
    }

    /// Least upper bound. `Undetermined` absorbs; `Natural` joins only
    /// itself (else `Undetermined`); mass sets join by union.
    fn join(&self, other: &Self) -> Self {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        use IsotopeMassAst::*;
        match (a.as_ref(), b.as_ref()) {
            (Undetermined, _) | (_, Undetermined) => Undetermined,
            (Natural, Natural) => Natural,
            (Natural, _) | (_, Natural) => Undetermined,
            (Var(x), Var(y)) if x == y => a.as_ref().clone(),
            (Var(_), _) | (_, Var(_)) => Undetermined,
            _ => {
                let sa = mass_set_view(a.as_ref()).unwrap();
                let sb = mass_set_view(b.as_ref()).unwrap();
                let union: BTreeSet<u32> = sa.union(&sb).copied().collect();
                canon_mass_set(union).unwrap_or(Undetermined)
            }
        }
    }

    /// `target` refines `self`: `self.meet(target) == canonical(target)`.
    fn matches(&self, target: &Self) -> bool {
        match (self.meet(target), target.canonical()) {
            (Some(meet), Ok(target)) => meet == *target,
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::{AtomConstraint, AtomConstraintKind};
    use crate::ast::value::ValueTerm;
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
    #[case::with_isotope_mass(AtomAst::default().with_isotope_mass(12_u32), AtomAst { isotope_mass: IsotopeMassAst::Lit(12), ..Default::default() })]
    #[case::with_charge(AtomAst::default().with_charge(1_i64), AtomAst { charge: ValueAst::Lit(1), ..Default::default() })]
    #[case::with_implicit_hydrogens(AtomAst::default().with_implicit_hydrogens(3_i64), AtomAst { implicit_hydrogens: ValueAst::Lit(3), ..Default::default() })]
    #[case::with_lone_pairs(AtomAst::default().with_lone_pairs(2_i64), AtomAst { lone_pairs: ValueAst::Lit(2), ..Default::default() })]
    #[case::with_spin_tuple(AtomAst::default().with_spin((0_u8, 1_u8)), AtomAst { spin: SpinStateAst::from((0_u8, 1_u8)), ..Default::default() })]
    #[case::with_constraint(AtomAst::default().with_constraint(AtomConstraint::valence(4_i64)),
        AtomAst { constraints: AtomConstraints::from(AtomConstraint::valence(4)),..Default::default() })]
    #[case::with_constraints_extends(AtomAst::default().with_constraint(AtomConstraint::valence(4_i64)).with_constraints([AtomConstraint::donated_pairs(1_i64), AtomConstraint::ring_size(6_i64)]),
        AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::valence(4), AtomConstraint::donated_pairs(1), AtomConstraint::ring_size(6)]), ..Default::default() })]
    #[case::with_constraint_replaces_same_kind(AtomAst::default().with_constraint(AtomConstraint::valence(3_i64)).with_constraint(AtomConstraint::valence(4_i64)),
        AtomAst { constraints: AtomConstraints::from(AtomConstraint::valence(4)), ..Default::default() })]
    fn test_atom_ast_with_methods(#[case] actual: AtomAst, #[case] expected: AtomAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_element(AtomAst::from_element(Element::C).into_ground(),
        AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeMassAst::Natural, charge: ValueAst::Lit(0), implicit_hydrogens: ValueAst::Lit(0),
        lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AtomConstraints::new() })]
    #[case::with_charge(AtomAst::from_element(Element::C).with_charge(1_i64).into_ground(),
        AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeMassAst::Natural, charge: ValueAst::Lit(1), implicit_hydrogens: ValueAst::Lit(0),
        lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AtomConstraints::new() })]
    #[case::constraint(AtomAst::from_element(Element::C).with_constraint(AtomConstraint::valence(4_i64)).into_ground(),
        AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeMassAst::Natural, charge: ValueAst::Lit(0), implicit_hydrogens: ValueAst::Lit(0),
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
    #[case::all_ground(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeMassAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ValueAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, true)]
    #[case::element_undetermined(AtomAst { element: ElementAst::Undetermined, isotope_mass: IsotopeMassAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ValueAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::isotope_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeMassAst::Undetermined, charge: ValueAst::Lit(0),
        implicit_hydrogens: ValueAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::charge_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeMassAst::Lit(12), charge: ValueAst::Undetermined,
        implicit_hydrogens: ValueAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::hydrogens_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeMassAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ValueAst::Undetermined, lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::lone_pairs_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeMassAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ValueAst::Lit(4), lone_pairs: ValueAst::Undetermined, spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AtomConstraints::new() }, false)]
    #[case::spin_undetermined(AtomAst { element: ElementAst::Lit(Element::C), isotope_mass: IsotopeMassAst::Lit(12), charge: ValueAst::Lit(0),
        implicit_hydrogens: ValueAst::Lit(4), lone_pairs: ValueAst::Lit(0), spin: SpinStateAst::default(),
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
    #[case::isotope_mismatch(AtomAst::from_element(Element::C).with_isotope_mass(12_u32), AtomAst::from_element(Element::C).with_isotope_mass(13_u32), false)]
    #[case::hydrogens_mismatch(AtomAst::from_element(Element::C).with_implicit_hydrogens(3_i64), AtomAst::from_element(Element::C).with_implicit_hydrogens(4_i64), false)]
    #[case::lone_pairs_mismatch(AtomAst::from_element(Element::C).with_lone_pairs(1_i64), AtomAst::from_element(Element::C).with_lone_pairs(2_i64), false)]
    #[case::spin_mismatch(AtomAst::from_element(Element::C).with_spin((2_u8, 3_u8)), AtomAst::from_element(Element::C).with_spin((0_u8, 1_u8)), false)]
    #[case::constraint_required_present(
        AtomAst::from_element(Element::C).with_constraint(AtomConstraint::valence(4)),
        AtomAst::from_element(Element::C).with_constraint(AtomConstraint::valence(4)),
        true)]
    #[case::constraint_required_absent(
        AtomAst::from_element(Element::C).with_constraint(AtomConstraint::valence(4)),
        AtomAst::from_element(Element::C),
        false)]
    #[case::constraint_value_mismatch(
        AtomAst::from_element(Element::C).with_constraint(AtomConstraint::valence(4)),
        AtomAst::from_element(Element::C).with_constraint(AtomConstraint::valence(3)),
        false)]
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
            isotope_mass: IsotopeMassAst::Lit(12),
            charge: ValueAst::term(ValueTerm::Lit(1)),
            implicit_hydrogens: ValueAst::term(ValueTerm::Lit(3)),
            lone_pairs: ValueAst::term(ValueTerm::Neg(Box::new(ValueTerm::Lit(2)))),
            spin: SpinStateAst {
                unpaired: ValueAst::term(ValueTerm::Lit(0)),
                multiplicity: ValueAst::term(ValueTerm::Lit(1)),
            },
            constraints: AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::term(
                ValueTerm::Lit(4),
            ))]),
        };
        atom.simplify_values();
        assert_eq!(atom.isotope_mass, IsotopeMassAst::Lit(12));
        assert_eq!(atom.charge, ValueAst::Lit(1));
        assert_eq!(atom.implicit_hydrogens, ValueAst::Lit(3));
        assert_eq!(atom.lone_pairs, ValueAst::Lit(-2));
        assert_eq!(atom.spin.unpaired, ValueAst::Lit(0));
        assert_eq!(atom.spin.multiplicity, ValueAst::Lit(1));
        assert_eq!(
            atom.constraints.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::valence(4)),
        );
    }

    #[rstest]
    #[case::both_default(AtomAst::default(), AtomAst::default(), Some(AtomAst::default()))]
    #[case::element_mismatch(
        AtomAst::from_element(Element::C),
        AtomAst::from_element(Element::N),
        None
    )]
    #[case::narrows_charge(AtomAst::from_element(Element::C), AtomAst::from_element(Element::C).with_charge(1),
        Some(AtomAst::from_element(Element::C).with_charge(1)))]
    fn test_atom_ast_meet(
        #[case] a: AtomAst,
        #[case] b: AtomAst,
        #[case] expected: Option<AtomAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::element_mismatch_widens(AtomAst::from_element(Element::C), AtomAst::from_element(Element::N),
        ElementAst::lit_set(vec![Element::C, Element::N]))]
    fn test_atom_ast_join_element(
        #[case] a: AtomAst,
        #[case] b: AtomAst,
        #[case] expected: ElementAst,
    ) {
        assert_eq!(a.join(&b).element, expected);
    }

    #[rstest]
    #[case::charge_change(AtomAst::from_element(Element::C), AtomAst::from_element(Element::C).with_charge(1), true,
        AtomAst::from_element(Element::C).with_charge(1))]
    #[case::no_change(
        AtomAst::from_element(Element::C),
        AtomAst::from_element(Element::C),
        false,
        AtomAst::from_element(Element::C)
    )]
    fn test_atom_ast_narrow_from(
        #[case] mut target: AtomAst,
        #[case] source: AtomAst,
        #[case] expected_changed: bool,
        #[case] expected_after: AtomAst,
    ) {
        let changed = target.narrow_from(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }


    #[rustfmt::skip]
    #[rstest]
    #[case::lit_set(ElementAst::lit_set([Element::C, Element::N]), ElementAst::LitSet(Box::new(BTreeSet::from([Element::C, Element::N]))))]
    #[case::not(ElementAst::not(Element::H), ElementAst::NotSet(Box::new(BTreeSet::from([Element::H]))))]
    #[case::not_set(ElementAst::not_set([Element::F, Element::Cl]), ElementAst::NotSet(Box::new(BTreeSet::from([Element::F, Element::Cl]))))]
    #[case::var(ElementAst::var("x"), ElementAst::Var(Box::new(("x".to_string(), None))))]
    #[case::var_in(ElementAst::var_in("x", [Element::C]), ElementAst::Var(Box::new(("x".to_string(), Some((MemOp::In, BTreeSet::from([Element::C])))))))]
    #[case::var_not_in(ElementAst::var_not_in("x", [Element::C]), ElementAst::Var(Box::new(("x".to_string(), Some((MemOp::NotIn, BTreeSet::from([Element::C])))))))]
    fn test_element_ast_constructors(#[case] actual: ElementAst, #[case] expected: ElementAst) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::carbon(Element::C, ElementAst::Lit(Element::C))]
    #[case::nitrogen(Element::N, ElementAst::Lit(Element::N))]
    fn test_element_ast_from(#[case] element: Element, #[case] expected: ElementAst) {
        assert_eq!(ElementAst::from(element), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::litset_singleton(ElementAst::lit_set([Element::C]), Ok(ElementAst::Lit(Element::C)))]
    #[case::litset_empty(ElementAst::lit_set([]), Err(Contradiction))]
    #[case::notset_empty(ElementAst::not_set([]), Ok(ElementAst::Undetermined))]
    #[case::var_in_empty(ElementAst::var_in("x", []), Err(Contradiction))]
    #[case::var_not_in_vacuous(ElementAst::var_not_in("x", []), Ok(ElementAst::var("x")))]
    fn test_element_ast_canonicalize(
        #[case] input: ElementAst,
        #[case] expected: Result<ElementAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ElementAst::Undetermined)]
    #[case::lit(ElementAst::Lit(Element::C))]
    #[case::litset(ElementAst::lit_set([Element::C, Element::N]))]
    #[case::notset(ElementAst::not(Element::H))]
    #[case::var_free(ElementAst::var("x"))]
    #[case::var_in(ElementAst::var_in("x", [Element::C]))]
    fn test_element_ast_canonicalize_identity(#[case] input: ElementAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    /// Cardinality polarity + universe boundaries (need the 118-element universe,
    /// so expected sets are computed from `Element::all()`, not hardcoded).
    #[rstest]
    fn test_element_ast_canonicalize_cardinality() {
        let take = |n: usize| -> BTreeSet<Element> { Element::all().iter().take(n).copied().collect() };
        let skip = |n: usize| -> BTreeSet<Element> { Element::all().iter().skip(n).copied().collect() };

        // Tiebreak: 59 stays positive; 60 flips to the complement.
        assert_eq!(ElementAst::LitSet(Box::new(take(59))).canonicalize(), Ok(ElementAst::LitSet(Box::new(take(59)))));
        assert_eq!(ElementAst::LitSet(Box::new(take(60))).canonicalize(), Ok(ElementAst::NotSet(Box::new(skip(60)))));
        // Full positive set → Undetermined; NotSet of the full set → Err.
        assert_eq!(ElementAst::LitSet(Box::new(take(118))).canonicalize(), Ok(ElementAst::Undetermined));
        assert_eq!(ElementAst::NotSet(Box::new(take(118))).canonicalize(), Err(Contradiction));
        // Large NotSet flips to a positive LitSet of its (small) complement.
        assert_eq!(ElementAst::NotSet(Box::new(take(60))).canonicalize(), Ok(ElementAst::LitSet(Box::new(skip(60)))));
        // Var In over a large domain flips to NotIn of the complement; full domain → free.
        assert_eq!(ElementAst::var_in("x", take(60)).canonicalize(), Ok(ElementAst::var_not_in("x", skip(60))));
        assert_eq!(ElementAst::var_in("x", take(118)).canonicalize(), Ok(ElementAst::var("x")));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_carbon(ElementAst::Lit(Element::C), Some(Element::C))]
    #[case::lit_nitrogen(ElementAst::Lit(Element::N), Some(Element::N))]
    #[case::undetermined(ElementAst::Undetermined, None)]
    #[case::litset(ElementAst::lit_set([Element::C, Element::N]), None)]
    #[case::notset(ElementAst::not(Element::H), None)]
    #[case::var_in(ElementAst::var_in("e", [Element::C]), None)]
    #[case::var(ElementAst::var("e"), None)]
    fn test_element_ast_as_lit(#[case] ast: ElementAst, #[case] expected: Option<Element>) {
        assert_eq!(ast.as_lit(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_match(ElementAst::Lit(Element::C), Element::C, true)]
    #[case::lit_mismatch(ElementAst::Lit(Element::C), Element::N, false)]
    #[case::undetermined(ElementAst::Undetermined, Element::C, false)]
    #[case::litset(ElementAst::lit_set([Element::C, Element::N]), Element::C, false)]
    fn test_element_ast_as_lit_matches(
        #[case] ast: ElementAst,
        #[case] value: Element,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.as_lit_matches(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ElementAst::Lit(Element::C), false)]
    #[case::undetermined(ElementAst::Undetermined, true)]
    #[case::litset(ElementAst::lit_set([Element::C, Element::N]), false)]
    #[case::notset(ElementAst::not(Element::H), false)]
    #[case::var_in(ElementAst::var_in("e", [Element::C]), false)]
    #[case::var(ElementAst::var("e"), false)]
    fn test_element_ast_is_undetermined(#[case] ast: ElementAst, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ElementAst::Undetermined, ElementAst::Lit(Element::C), Some(ElementAst::Lit(Element::C)))]
    #[case::lit_und(ElementAst::Lit(Element::C), ElementAst::Undetermined, Some(ElementAst::Lit(Element::C)))]
    #[case::lit_lit_eq(ElementAst::Lit(Element::C), ElementAst::Lit(Element::C), Some(ElementAst::Lit(Element::C)))]
    #[case::lit_lit_neq(ElementAst::Lit(Element::C), ElementAst::Lit(Element::N), None)]
    #[case::lit_set_in(ElementAst::Lit(Element::C), ElementAst::lit_set([Element::C, Element::N]), Some(ElementAst::Lit(Element::C)))]
    #[case::lit_set_out(ElementAst::Lit(Element::O), ElementAst::lit_set([Element::C, Element::N]), None)]
    #[case::set_set_singleton(ElementAst::lit_set([Element::C, Element::N]), ElementAst::lit_set([Element::N, Element::O]), Some(ElementAst::Lit(Element::N)))]
    #[case::set_set_multi(ElementAst::lit_set([Element::C, Element::N, Element::O]), ElementAst::lit_set([Element::N, Element::O, Element::F]), Some(ElementAst::lit_set([Element::N, Element::O])))]
    #[case::set_set_disjoint(ElementAst::lit_set([Element::C, Element::N]), ElementAst::lit_set([Element::O, Element::F]), None)]
    #[case::set_notset(ElementAst::lit_set([Element::C, Element::N]), ElementAst::not(Element::N), Some(ElementAst::Lit(Element::C)))]
    #[case::notset_notset(ElementAst::not(Element::C), ElementAst::not(Element::N), Some(ElementAst::not_set([Element::C, Element::N])))]
    #[case::lit_notset_in(ElementAst::Lit(Element::C), ElementAst::not(Element::N), Some(ElementAst::Lit(Element::C)))]
    #[case::lit_notset_out(ElementAst::Lit(Element::C), ElementAst::not(Element::C), None)]
    #[case::var_var_eq(ElementAst::var("e"), ElementAst::var("e"), Some(ElementAst::var("e")))]
    #[case::var_var_neq(ElementAst::var("e"), ElementAst::var("f"), None)]
    #[case::var_lit(ElementAst::var("e"), ElementAst::Lit(Element::C), None)]
    #[case::lit_var(ElementAst::Lit(Element::C), ElementAst::var("e"), None)]
    fn test_element_ast_meet(
        #[case] a: ElementAst,
        #[case] b: ElementAst,
        #[case] expected: Option<ElementAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ElementAst::Undetermined, ElementAst::Lit(Element::C), ElementAst::Undetermined)]
    #[case::lit_lit_eq(ElementAst::Lit(Element::C), ElementAst::Lit(Element::C), ElementAst::Lit(Element::C))]
    #[case::lit_lit_neq(ElementAst::Lit(Element::C), ElementAst::Lit(Element::N), ElementAst::lit_set([Element::C, Element::N]))]
    #[case::lit_set(ElementAst::Lit(Element::O), ElementAst::lit_set([Element::C, Element::N]), ElementAst::lit_set([Element::C, Element::N, Element::O]))]
    #[case::set_set(ElementAst::lit_set([Element::C, Element::N]), ElementAst::lit_set([Element::N, Element::O]), ElementAst::lit_set([Element::C, Element::N, Element::O]))]
    #[case::lit_notset_out(ElementAst::Lit(Element::C), ElementAst::not(Element::N), ElementAst::not(Element::N))]
    #[case::lit_notset_in(ElementAst::Lit(Element::N), ElementAst::not(Element::N), ElementAst::Undetermined)]
    #[case::notset_notset_disjoint(ElementAst::not(Element::C), ElementAst::not(Element::N), ElementAst::Undetermined)]
    #[case::notset_notset_overlap(ElementAst::not_set([Element::C, Element::N]), ElementAst::not_set([Element::N, Element::O]), ElementAst::not(Element::N))]
    #[case::var_var_eq(ElementAst::var("e"), ElementAst::var("e"), ElementAst::var("e"))]
    #[case::var_var_neq(ElementAst::var("e"), ElementAst::var("f"), ElementAst::Undetermined)]
    #[case::var_lit(ElementAst::var("e"), ElementAst::Lit(Element::C), ElementAst::Undetermined)]
    fn test_element_ast_join(
        #[case] a: ElementAst,
        #[case] b: ElementAst,
        #[case] expected: ElementAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ElementAst::Undetermined, ElementAst::Lit(Element::C), true)]
    #[case::und_und(ElementAst::Undetermined, ElementAst::Undetermined, true)]
    #[case::und_set(ElementAst::Undetermined, ElementAst::lit_set([Element::C, Element::N]), true)]
    #[case::und_var(ElementAst::Undetermined, ElementAst::var("e"), true)]
    #[case::lit_und(ElementAst::Lit(Element::C), ElementAst::Undetermined, false)]
    #[case::set_und(ElementAst::lit_set([Element::C, Element::N]), ElementAst::Undetermined, false)]
    #[case::var_und(ElementAst::var("e"), ElementAst::Undetermined, false)]
    #[case::lit_lit_match(ElementAst::Lit(Element::C), ElementAst::Lit(Element::C), true)]
    #[case::lit_lit_mismatch(ElementAst::Lit(Element::C), ElementAst::Lit(Element::N), false)]
    #[case::lit_singleton_set(ElementAst::Lit(Element::C), ElementAst::lit_set([Element::C]), true)]
    #[case::lit_multi_set(ElementAst::Lit(Element::C), ElementAst::lit_set([Element::C, Element::N]), false)]
    #[case::set_lit_in(ElementAst::lit_set([Element::C, Element::N]), ElementAst::Lit(Element::N), true)]
    #[case::set_lit_out(ElementAst::lit_set([Element::C, Element::N]), ElementAst::Lit(Element::O), false)]
    #[case::set_set_subset(ElementAst::lit_set([Element::C, Element::N, Element::O]), ElementAst::lit_set([Element::C, Element::N]), true)]
    #[case::set_set_equal(ElementAst::lit_set([Element::C, Element::N]), ElementAst::lit_set([Element::C, Element::N]), true)]
    #[case::set_set_superset(ElementAst::lit_set([Element::C]), ElementAst::lit_set([Element::C, Element::N]), false)]
    #[case::notset_lit_admitted(ElementAst::not(Element::N), ElementAst::Lit(Element::C), true)]
    #[case::notset_lit_excluded(ElementAst::not(Element::C), ElementAst::Lit(Element::C), false)]
    #[case::var_var_equal(ElementAst::var("e"), ElementAst::var("e"), true)]
    #[case::var_var_distinct(ElementAst::var("e"), ElementAst::var("f"), false)]
    #[case::var_lit(ElementAst::var("e"), ElementAst::Lit(Element::C), false)]
    #[case::lit_var(ElementAst::Lit(Element::C), ElementAst::var("e"), false)]
    #[case::set_var(ElementAst::lit_set([Element::C]), ElementAst::var("e"), false)]
    fn test_element_ast_matches(
        #[case] pattern: ElementAst,
        #[case] target: ElementAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_set(IsotopeMassAst::lit_set([12, 13]), IsotopeMassAst::LitSet(Box::new(BTreeSet::from([12, 13]))))]
    #[case::var(IsotopeMassAst::var("m"), IsotopeMassAst::Var(Box::new(("m".to_string(), None))))]
    #[case::var_in(IsotopeMassAst::var_in("m", [12, 13]), IsotopeMassAst::Var(Box::new(("m".to_string(), Some(BTreeSet::from([12, 13]))))))]
    fn test_isotope_mass_ast_constructors(#[case] actual: IsotopeMassAst, #[case] expected: IsotopeMassAst) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::positive(13, IsotopeMassAst::Lit(13))]
    #[case::zero(0, IsotopeMassAst::Lit(0))]
    fn test_isotope_mass_ast_from(#[case] mass: u32, #[case] expected: IsotopeMassAst) {
        assert_eq!(IsotopeMassAst::from(mass), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::litset_singleton(IsotopeMassAst::lit_set([12]), Ok(IsotopeMassAst::Lit(12)))]
    #[case::litset_empty(IsotopeMassAst::lit_set([]), Err(Contradiction))]
    #[case::var_in_empty(IsotopeMassAst::var_in("m", []), Err(Contradiction))]
    fn test_isotope_mass_ast_canonicalize(
        #[case] input: IsotopeMassAst,
        #[case] expected: Result<IsotopeMassAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(IsotopeMassAst::Undetermined)]
    #[case::natural(IsotopeMassAst::Natural)]
    #[case::lit(IsotopeMassAst::Lit(12))]
    #[case::litset(IsotopeMassAst::lit_set([12, 13]))]
    #[case::var_free(IsotopeMassAst::var("m"))]
    #[case::var_in(IsotopeMassAst::var_in("m", [12, 13]))]
    #[case::var_in_singleton(IsotopeMassAst::var_in("m", [12]))]
    fn test_isotope_mass_ast_canonicalize_identity(#[case] input: IsotopeMassAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(IsotopeMassAst::Lit(12), Some(12))]
    #[case::lit_zero(IsotopeMassAst::Lit(0), Some(0))]
    #[case::natural(IsotopeMassAst::Natural, None)]
    #[case::undetermined(IsotopeMassAst::Undetermined, None)]
    #[case::litset(IsotopeMassAst::lit_set([12, 13]), None)]
    #[case::var(IsotopeMassAst::var("m"), None)]
    #[case::var_in(IsotopeMassAst::var_in("m", [12]), None)]
    fn test_isotope_mass_ast_as_lit(#[case] ast: IsotopeMassAst, #[case] expected: Option<u32>) {
        assert_eq!(ast.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_match(IsotopeMassAst::Lit(12), 12, true)]
    #[case::lit_mismatch(IsotopeMassAst::Lit(12), 13, false)]
    #[case::natural(IsotopeMassAst::Natural, 0, false)]
    #[case::undetermined(IsotopeMassAst::Undetermined, 12, false)]
    #[case::litset(IsotopeMassAst::lit_set([12, 13]), 12, false)]
    fn test_isotope_mass_ast_as_lit_matches(
        #[case] ast: IsotopeMassAst,
        #[case] value: u32,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.as_lit_matches(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(IsotopeMassAst::Undetermined, true)]
    #[case::natural(IsotopeMassAst::Natural, false)]
    #[case::lit(IsotopeMassAst::Lit(12), false)]
    #[case::litset(IsotopeMassAst::lit_set([12, 13]), false)]
    #[case::var(IsotopeMassAst::var("m"), false)]
    #[case::var_in(IsotopeMassAst::var_in("m", [12]), false)]
    fn test_isotope_mass_ast_is_undetermined(#[case] ast: IsotopeMassAst, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::natural(IsotopeMassAst::Natural, true)]
    #[case::lit(IsotopeMassAst::Lit(12), true)]
    #[case::undetermined(IsotopeMassAst::Undetermined, false)]
    #[case::litset(IsotopeMassAst::lit_set([12, 13]), false)]
    #[case::var(IsotopeMassAst::var("m"), false)]
    #[case::var_in(IsotopeMassAst::var_in("m", [12]), false)]
    fn test_isotope_mass_ast_is_ground(#[case] ast: IsotopeMassAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(IsotopeMassAst::Undetermined, IsotopeMassAst::Lit(12), Some(IsotopeMassAst::Lit(12)))]
    #[case::lit_und(IsotopeMassAst::Lit(12), IsotopeMassAst::Undetermined, Some(IsotopeMassAst::Lit(12)))]
    #[case::und_natural(IsotopeMassAst::Undetermined, IsotopeMassAst::Natural, Some(IsotopeMassAst::Natural))]
    #[case::natural_natural(IsotopeMassAst::Natural, IsotopeMassAst::Natural, Some(IsotopeMassAst::Natural))]
    #[case::natural_lit(IsotopeMassAst::Natural, IsotopeMassAst::Lit(12), None)]
    #[case::lit_natural(IsotopeMassAst::Lit(12), IsotopeMassAst::Natural, None)]
    #[case::lit_lit_eq(IsotopeMassAst::Lit(12), IsotopeMassAst::Lit(12), Some(IsotopeMassAst::Lit(12)))]
    #[case::lit_lit_neq(IsotopeMassAst::Lit(12), IsotopeMassAst::Lit(13), None)]
    #[case::lit_set_in(IsotopeMassAst::Lit(12), IsotopeMassAst::lit_set([12, 13]), Some(IsotopeMassAst::Lit(12)))]
    #[case::lit_set_out(IsotopeMassAst::Lit(14), IsotopeMassAst::lit_set([12, 13]), None)]
    #[case::set_set_singleton(IsotopeMassAst::lit_set([12, 13]), IsotopeMassAst::lit_set([13, 14]), Some(IsotopeMassAst::Lit(13)))]
    #[case::set_set_multi(IsotopeMassAst::lit_set([12, 13, 14]), IsotopeMassAst::lit_set([13, 14, 15]), Some(IsotopeMassAst::lit_set([13, 14])))]
    #[case::set_set_disjoint(IsotopeMassAst::lit_set([12, 13]), IsotopeMassAst::lit_set([14, 15]), None)]
    #[case::var_var_eq(IsotopeMassAst::var("m"), IsotopeMassAst::var("m"), Some(IsotopeMassAst::var("m")))]
    #[case::var_var_neq(IsotopeMassAst::var("m"), IsotopeMassAst::var("n"), None)]
    #[case::var_lit(IsotopeMassAst::var("m"), IsotopeMassAst::Lit(12), None)]
    fn test_isotope_mass_ast_meet(
        #[case] a: IsotopeMassAst,
        #[case] b: IsotopeMassAst,
        #[case] expected: Option<IsotopeMassAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(IsotopeMassAst::Undetermined, IsotopeMassAst::Lit(12), IsotopeMassAst::Undetermined)]
    #[case::natural_natural(IsotopeMassAst::Natural, IsotopeMassAst::Natural, IsotopeMassAst::Natural)]
    #[case::natural_lit(IsotopeMassAst::Natural, IsotopeMassAst::Lit(12), IsotopeMassAst::Undetermined)]
    #[case::lit_lit_eq(IsotopeMassAst::Lit(12), IsotopeMassAst::Lit(12), IsotopeMassAst::Lit(12))]
    #[case::lit_lit_neq(IsotopeMassAst::Lit(12), IsotopeMassAst::Lit(13), IsotopeMassAst::lit_set([12, 13]))]
    #[case::lit_set(IsotopeMassAst::Lit(14), IsotopeMassAst::lit_set([12, 13]), IsotopeMassAst::lit_set([12, 13, 14]))]
    #[case::set_set(IsotopeMassAst::lit_set([12, 13]), IsotopeMassAst::lit_set([13, 14]), IsotopeMassAst::lit_set([12, 13, 14]))]
    #[case::var_var_eq(IsotopeMassAst::var("m"), IsotopeMassAst::var("m"), IsotopeMassAst::var("m"))]
    #[case::var_var_neq(IsotopeMassAst::var("m"), IsotopeMassAst::var("n"), IsotopeMassAst::Undetermined)]
    #[case::var_lit(IsotopeMassAst::var("m"), IsotopeMassAst::Lit(12), IsotopeMassAst::Undetermined)]
    fn test_isotope_mass_ast_join(
        #[case] a: IsotopeMassAst,
        #[case] b: IsotopeMassAst,
        #[case] expected: IsotopeMassAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(IsotopeMassAst::Undetermined, IsotopeMassAst::Lit(12), true)]
    #[case::und_natural(IsotopeMassAst::Undetermined, IsotopeMassAst::Natural, true)]
    #[case::und_und(IsotopeMassAst::Undetermined, IsotopeMassAst::Undetermined, true)]
    #[case::und_set(IsotopeMassAst::Undetermined, IsotopeMassAst::lit_set([12, 13]), true)]
    #[case::und_var(IsotopeMassAst::Undetermined, IsotopeMassAst::var("m"), true)]
    #[case::lit_und(IsotopeMassAst::Lit(12), IsotopeMassAst::Undetermined, false)]
    #[case::natural_und(IsotopeMassAst::Natural, IsotopeMassAst::Undetermined, false)]
    #[case::natural_natural(IsotopeMassAst::Natural, IsotopeMassAst::Natural, true)]
    #[case::natural_lit(IsotopeMassAst::Natural, IsotopeMassAst::Lit(12), false)]
    #[case::lit_natural(IsotopeMassAst::Lit(12), IsotopeMassAst::Natural, false)]
    #[case::lit_lit_match(IsotopeMassAst::Lit(12), IsotopeMassAst::Lit(12), true)]
    #[case::lit_lit_mismatch(IsotopeMassAst::Lit(12), IsotopeMassAst::Lit(13), false)]
    #[case::set_lit_in(IsotopeMassAst::lit_set([12, 13]), IsotopeMassAst::Lit(13), true)]
    #[case::set_lit_out(IsotopeMassAst::lit_set([12, 13]), IsotopeMassAst::Lit(14), false)]
    #[case::set_set_subset(IsotopeMassAst::lit_set([12, 13, 14]), IsotopeMassAst::lit_set([12, 13]), true)]
    #[case::set_set_superset(IsotopeMassAst::lit_set([12]), IsotopeMassAst::lit_set([12, 13]), false)]
    #[case::var_var_equal(IsotopeMassAst::var("m"), IsotopeMassAst::var("m"), true)]
    #[case::var_var_distinct(IsotopeMassAst::var("m"), IsotopeMassAst::var("n"), false)]
    #[case::var_lit(IsotopeMassAst::var("m"), IsotopeMassAst::Lit(12), false)]
    #[case::lit_var(IsotopeMassAst::Lit(12), IsotopeMassAst::var("m"), false)]
    fn test_isotope_mass_ast_matches(
        #[case] pattern: IsotopeMassAst,
        #[case] target: IsotopeMassAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
