//! Atom-level AST fragments shared across crates.

use std::borrow::Cow;
use std::collections::BTreeSet;

use umol_chem::element::{Element, MAX_ATOMIC_NUMBER};
use umol_graph_ir_macros::{Canonicalize, Lattice};

use super::constraint::{AtomConstraintForm, AtomConstraintsForm};
use super::error::{Contradiction, NoJoin};
use super::num::NumForm;
use super::operators::MemOp;
use super::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
use super::traits::{AsLit, Canonicalize, Lattice};

/// Atom AST: structural representation of an atom plus the atom-level
/// constraints (valence, degree, ring membership, etc.) that pattern
/// against the surrounding topology.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
pub struct AtomForm {
    pub element: ElementForm,
    pub isotope_mass: IsotopeMassForm,
    pub charge: NumForm,
    pub implicit_hydrogens: NumForm,
    pub lone_pairs: NumForm,
    pub unpaired_electrons: UnpairedElectronsForm,
    pub constraints: AtomConstraintsForm,
}

/// Attribute update for an atom. `None` leaves an ordinary scalar field unchanged;
/// `Some(value)` sets it exactly, including to an undetermined value. Unpaired-electron components
/// are updated independently. Constraints are keyed updates, with undetermined entries removing
/// their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomUpdate {
    pub element: Option<ElementForm>,
    pub isotope_mass: Option<IsotopeMassForm>,
    pub charge: Option<NumForm>,
    pub implicit_hydrogens: Option<NumForm>,
    pub lone_pairs: Option<NumForm>,
    pub unpaired_electrons: UnpairedElectronsUpdate,
    pub constraints: AtomConstraintsForm,
}

impl AtomForm {
    pub fn new(element: ElementForm) -> Self {
        Self {
            element,
            ..Default::default()
        }
    }

    pub fn from_element(element: Element) -> Self {
        Self::new(ElementForm::Lit(element))
    }

    pub fn with_element(mut self, element: impl Into<ElementForm>) -> Self {
        self.element = element.into();
        self
    }

    pub fn with_isotope_mass(mut self, mass: impl Into<IsotopeMassForm>) -> Self {
        self.isotope_mass = mass.into();
        self
    }

    pub fn with_charge(mut self, charge: impl Into<NumForm>) -> Self {
        self.charge = charge.into();
        self
    }

    pub fn with_implicit_hydrogens(mut self, hydrogens: impl Into<NumForm>) -> Self {
        self.implicit_hydrogens = hydrogens.into();
        self
    }

    pub fn with_lone_pairs(mut self, lone_pairs: impl Into<NumForm>) -> Self {
        self.lone_pairs = lone_pairs.into();
        self
    }

    pub fn with_unpaired_electrons(
        mut self,
        unpaired_electrons: impl Into<UnpairedElectronsForm>,
    ) -> Self {
        self.unpaired_electrons = unpaired_electrons.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// kind (last-wins per `AtomConstraintsForm::add`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<AtomConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `AtomConstraintsForm::add`). Does not
    /// clear existing constraints; use `atom.constraints.clear()` or direct
    /// field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<AtomConstraintForm>,
    {
        self.constraints
            .extend(constraints.into_iter().map(Into::into));
        self
    }

    /// Apply an attribute update, leaving omitted fields and constraint keys unchanged.
    pub fn update(&self, update: &AtomUpdate) -> AtomForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        AtomForm {
            element: update
                .element
                .clone()
                .unwrap_or_else(|| self.element.clone()),
            isotope_mass: update
                .isotope_mass
                .clone()
                .unwrap_or_else(|| self.isotope_mass.clone()),
            charge: update.charge.clone().unwrap_or_else(|| self.charge.clone()),
            implicit_hydrogens: update
                .implicit_hydrogens
                .clone()
                .unwrap_or_else(|| self.implicit_hydrogens.clone()),
            lone_pairs: update
                .lone_pairs
                .clone()
                .unwrap_or_else(|| self.lone_pairs.clone()),
            unpaired_electrons: self.unpaired_electrons.update(&update.unpaired_electrons),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> AtomUpdate {
        let mut constraints = AtomConstraintsForm::new();
        for new in other.constraints.iter() {
            if self
                .constraints
                .get(new.key())
                .is_none_or(|old| !old.canonical_eq(new))
            {
                constraints.set(new.clone());
            }
        }
        for old in self.constraints.iter() {
            if other.constraints.get(old.key()).is_none() {
                constraints.set(old.as_undetermined());
            }
        }
        AtomUpdate {
            element: (!self.element.canonical_eq(&other.element)).then(|| other.element.clone()),
            isotope_mass: (!self.isotope_mass.canonical_eq(&other.isotope_mass))
                .then(|| other.isotope_mass.clone()),
            charge: (!self.charge.canonical_eq(&other.charge)).then(|| other.charge.clone()),
            implicit_hydrogens: (!self
                .implicit_hydrogens
                .canonical_eq(&other.implicit_hydrogens))
            .then(|| other.implicit_hydrogens.clone()),
            lone_pairs: (!self.lone_pairs.canonical_eq(&other.lone_pairs))
                .then(|| other.lone_pairs.clone()),
            unpaired_electrons: self
                .unpaired_electrons
                .difference_to(&other.unpaired_electrons),
            constraints,
        }
    }

    /// Fill `Undetermined` value-bearing struct fields with defaults: isotope→
    /// Natural; charge / implicit hydrogens / lone pairs → 0; unpaired-electron count → 0
    /// and, for the (possibly already-set) count, the maximal
    /// multiplicity `count + 1` (so a fully unset pair becomes the
    /// closed-shell singlet). Existing literal or expression values and all
    /// constraints are preserved. The result is ground iff `element` is
    /// already ground.
    pub fn into_ground(mut self) -> Self {
        if self.isotope_mass.is_undetermined() {
            self.isotope_mass = IsotopeMassForm::Natural;
        }
        if self.charge.is_undetermined() {
            self.charge = NumForm::Lit(0);
        }
        if self.implicit_hydrogens.is_undetermined() {
            self.implicit_hydrogens = NumForm::Lit(0);
        }
        if self.lone_pairs.is_undetermined() {
            self.lone_pairs = NumForm::Lit(0);
        }
        if self.unpaired_electrons.count.is_undetermined() {
            self.unpaired_electrons.count = NumForm::Lit(0);
        }
        if self.unpaired_electrons.multiplicity.is_undetermined() {
            let count = self.unpaired_electrons.count.as_lit().unwrap_or(0);
            self.unpaired_electrons.multiplicity = NumForm::Lit(count + 1);
        }
        self
    }
}

impl From<Element> for AtomForm {
    fn from(element: Element) -> Self {
        Self::from_element(element)
    }
}

/// Construction sugar for `b.atom("C#h3")`: parse a compact atom-string, panicking on
/// invalid input — a bad literal is a programmer error, like the `atom_dsl!` macro.
impl From<&str> for AtomForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid atom string")
    }
}

/// Element expression: undetermined, a single element, a finite element set, a
/// complement set (`!{…}`), or a variable (free `?x`, or membership-restricted
/// `?x :: {…}` / `?x :: !{…}`). Sets are cardinality-canonical and
/// universe-relative: a semantic set of more than `⌊118/2⌋` elements is stored
/// as its complement `NotSet(U∖S)`.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ElementForm {
    #[default]
    Undetermined,
    Lit(Element),
    LitSet(Box<BTreeSet<Element>>),
    NotSet(Box<BTreeSet<Element>>),
    #[allow(clippy::type_complexity)]
    Var(Box<(String, Option<(MemOp, BTreeSet<Element>)>)>),
}

impl ElementForm {
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
    pub fn var_not_in(
        name: impl Into<String>,
        elements: impl IntoIterator<Item = Element>,
    ) -> Self {
        Self::Var(Box::new((
            name.into(),
            Some((MemOp::NotIn, elements.into_iter().collect())),
        )))
    }
}

impl From<Element> for ElementForm {
    fn from(element: Element) -> Self {
        Self::Lit(element)
    }
}

impl Canonicalize for ElementForm {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            ElementForm::Undetermined => ElementForm::Undetermined,
            ElementForm::Lit(e) => ElementForm::Lit(e),
            ElementForm::LitSet(s) => canon_set(*s, false)?,
            ElementForm::NotSet(s) => canon_set(*s, true)?,
            ElementForm::Var(v) => {
                let (name, domain) = *v;
                let domain = match domain {
                    None => None,
                    Some((op, set)) => canon_var_domain(op, set)?,
                };
                ElementForm::Var(Box::new((name, domain)))
            }
        })
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            ElementForm::Undetermined | ElementForm::Lit(_) => Ok(Cow::Borrowed(self)),
            _ => Ok(Cow::Owned(self.clone().canonicalize()?)),
        }
    }
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
fn canon_set(s: BTreeSet<Element>, negated: bool) -> Result<ElementForm, Contradiction> {
    let semantic_len = if negated {
        MAX_ATOMIC_NUMBER as usize - s.len()
    } else {
        s.len()
    };
    if semantic_len == 0 {
        return Err(Contradiction);
    }
    if semantic_len == MAX_ATOMIC_NUMBER as usize {
        return Ok(ElementForm::Undetermined);
    }
    if semantic_len <= MAX_ATOMIC_NUMBER as usize / 2 {
        let positive = if negated { complement(&s) } else { s };
        Ok(if positive.len() == 1 {
            ElementForm::Lit(*positive.iter().next().unwrap())
        } else {
            ElementForm::LitSet(Box::new(positive))
        })
    } else {
        let excluded = if negated { s } else { complement(&s) };
        Ok(ElementForm::NotSet(Box::new(excluded)))
    }
}

/// Canonicalize a `Var` membership domain by the same cardinality polarity.
/// Vacuous (`In U` / `NotIn ∅`) → `None` (free); empty admissible set → `Err`.
fn canon_var_domain(
    op: MemOp,
    set: BTreeSet<Element>,
) -> Result<Option<(MemOp, BTreeSet<Element>)>, Contradiction> {
    let admissible_len = match op {
        MemOp::In => set.len(),
        MemOp::NotIn => MAX_ATOMIC_NUMBER as usize - set.len(),
    };
    if admissible_len == 0 {
        return Err(Contradiction);
    }
    if admissible_len == MAX_ATOMIC_NUMBER as usize {
        return Ok(None);
    }
    Ok(Some(if admissible_len <= MAX_ATOMIC_NUMBER as usize / 2 {
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

impl AsLit for ElementForm {
    type Lit = Element;

    /// The single element this denotes, only when it is a literal.
    /// Non-canonicalizing (mirrors `NumForm::as_lit`).
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
fn element_set_view(e: &ElementForm) -> Option<(BTreeSet<Element>, bool)> {
    match e {
        ElementForm::Lit(x) => Some((BTreeSet::from([*x]), false)),
        ElementForm::LitSet(s) => Some(((**s).clone(), false)),
        ElementForm::NotSet(s) => Some(((**s).clone(), true)),
        ElementForm::Undetermined | ElementForm::Var(_) => None,
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

impl Lattice for ElementForm {
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
        use ElementForm::*;
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
                let raw = if negated {
                    NotSet(Box::new(set))
                } else {
                    LitSet(Box::new(set))
                };
                raw.canonicalize().ok()
            }
        }
    }

    /// Least upper bound (set union), canonicalizing operands and result.
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        use ElementForm::*;
        Ok(match (a.as_ref(), b.as_ref()) {
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
                let raw = if negated {
                    NotSet(Box::new(set))
                } else {
                    LitSet(Box::new(set))
                };
                raw.canonicalize().unwrap_or(Undetermined)
            }
        })
    }

    /// Partial-order check `target ⊑ self`, allocation-free for the literal
    /// cases. `NotSet` (complement) and `Var` on either side fall back to the
    /// canonicalizing `meet`-derived default, which this must equal.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::NotSet(_) | Self::Var(_), _) | (_, Self::NotSet(_) | Self::Var(_)) => {
                match (self.meet(target), target.canonical()) {
                    (Some(meet), Ok(target)) => meet == *target,
                    _ => false,
                }
            }
            (Self::Undetermined, Self::Undetermined | Self::Lit(_)) => true,
            (Self::Lit(_), Self::Undetermined) => false,
            (Self::Lit(p), Self::Lit(t)) => p == t,
            (Self::Undetermined, Self::LitSet(t)) => !t.is_empty(),
            (Self::LitSet(_), Self::Undetermined) => false,
            (Self::Lit(p), Self::LitSet(t)) => t.len() == 1 && t.contains(p),
            (Self::LitSet(p), Self::Lit(t)) => p.contains(t),
            (Self::LitSet(p), Self::LitSet(t)) => !t.is_empty() && t.iter().all(|x| p.contains(x)),
        }
    }
}

/// Isotope-mass expression: undetermined, the natural isotopic mixture
/// (`#i=`), a single mass number, a finite mass set, or a variable (free
/// `?x`, or membership-restricted `?x :: {…}`). Positive-only — no negation
/// and no complement (the mass domain is open). `Natural` is a distinct
/// ground, disjoint from every specific mass.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IsotopeMass {
    Natural,
    MassNumber(u32),
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IsotopeMassForm {
    #[default]
    Undetermined,
    Natural,
    Lit(u32),
    LitSet(Box<BTreeSet<u32>>),
    Var(Box<(String, Option<BTreeSet<u32>>)>),
}

impl IsotopeMassForm {
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

impl From<u32> for IsotopeMassForm {
    fn from(mass: u32) -> Self {
        Self::Lit(mass)
    }
}

impl From<IsotopeMass> for IsotopeMassForm {
    fn from(mass: IsotopeMass) -> Self {
        match mass {
            IsotopeMass::Natural => Self::Natural,
            IsotopeMass::MassNumber(mass) => Self::Lit(mass),
        }
    }
}

impl Canonicalize for IsotopeMassForm {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            IsotopeMassForm::Undetermined => IsotopeMassForm::Undetermined,
            IsotopeMassForm::Natural => IsotopeMassForm::Natural,
            IsotopeMassForm::Lit(n) => IsotopeMassForm::Lit(n),
            IsotopeMassForm::LitSet(s) => canon_mass_set(*s)?,
            IsotopeMassForm::Var(v) => {
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
                IsotopeMassForm::Var(Box::new((name, domain)))
            }
        })
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            IsotopeMassForm::Undetermined | IsotopeMassForm::Natural | IsotopeMassForm::Lit(_) => {
                Ok(Cow::Borrowed(self))
            }
            _ => Ok(Cow::Owned(self.clone().canonicalize()?)),
        }
    }
}

/// Canonicalize a mass set: empty → `Err`; singleton → `Lit`; else `LitSet`.
fn canon_mass_set(s: BTreeSet<u32>) -> Result<IsotopeMassForm, Contradiction> {
    match s.len() {
        0 => Err(Contradiction),
        1 => Ok(IsotopeMassForm::Lit(*s.iter().next().unwrap())),
        _ => Ok(IsotopeMassForm::LitSet(Box::new(s))),
    }
}

impl AsLit for IsotopeMassForm {
    type Lit = IsotopeMass;

    /// The exact natural-composition or mass-number value when ground.
    #[inline]
    fn as_lit(&self) -> Option<IsotopeMass> {
        match self {
            Self::Natural => Some(IsotopeMass::Natural),
            Self::Lit(n) => Some(IsotopeMass::MassNumber(*n)),
            _ => None,
        }
    }
}

/// The positive mass set a concrete form denotes. `Undetermined`, `Natural`,
/// and `Var` have no mass-set view.
fn mass_set_view(iso: &IsotopeMassForm) -> Option<BTreeSet<u32>> {
    match iso {
        IsotopeMassForm::Lit(n) => Some(BTreeSet::from([*n])),
        IsotopeMassForm::LitSet(s) => Some((**s).clone()),
        IsotopeMassForm::Undetermined | IsotopeMassForm::Natural | IsotopeMassForm::Var(_) => None,
    }
}

impl Lattice for IsotopeMassForm {
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
        use IsotopeMassForm::*;
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
                LitSet(Box::new(intersection)).canonicalize().ok()
            }
        }
    }

    /// Least upper bound. `Undetermined` absorbs; `Natural` joins only
    /// itself (else `Undetermined`); mass sets join by union.
    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        use IsotopeMassForm::*;
        Ok(match (a.as_ref(), b.as_ref()) {
            (Undetermined, _) | (_, Undetermined) => Undetermined,
            (Natural, Natural) => Natural,
            (Natural, _) | (_, Natural) => Undetermined,
            (Var(x), Var(y)) if x == y => a.as_ref().clone(),
            (Var(_), _) | (_, Var(_)) => Undetermined,
            _ => {
                let sa = mass_set_view(a.as_ref()).unwrap();
                let sb = mass_set_view(b.as_ref()).unwrap();
                let union: BTreeSet<u32> = sa.union(&sb).copied().collect();
                LitSet(Box::new(union))
                    .canonicalize()
                    .unwrap_or(Undetermined)
            }
        })
    }

    /// Partial-order check `target ⊑ self`, allocation-free for the literal
    /// cases. `Natural` is an isolated ground (matches only `Natural`); only
    /// `Var` on either side falls back to the `meet`-derived default.
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Var(_), _) | (_, Self::Var(_)) => {
                match (self.meet(target), target.canonical()) {
                    (Some(meet), Ok(target)) => meet == *target,
                    _ => false,
                }
            }
            (Self::Undetermined, Self::Undetermined | Self::Natural | Self::Lit(_)) => true,
            (Self::Natural, Self::Natural) => true,
            (Self::Natural, _) | (_, Self::Natural) => false,
            (Self::Lit(_), Self::Undetermined) => false,
            (Self::Lit(p), Self::Lit(t)) => p == t,
            (Self::Undetermined, Self::LitSet(t)) => !t.is_empty(),
            (Self::LitSet(_), Self::Undetermined) => false,
            (Self::Lit(p), Self::LitSet(t)) => t.len() == 1 && t.contains(p),
            (Self::LitSet(p), Self::Lit(t)) => p.contains(t),
            (Self::LitSet(p), Self::LitSet(t)) => !t.is_empty() && t.iter().all(|x| p.contains(x)),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::constraint::{AtomConstraintForm, RingScope};

    #[rstest]
    fn test_atom_form_from_element() {
        assert_eq!(
            AtomForm::from_element(Element::C),
            AtomForm {
                element: ElementForm::Lit(Element::C),
                ..Default::default()
            },
        );
    }

    #[rstest]
    fn test_atom_form_from() {
        let expected = AtomForm {
            element: ElementForm::Lit(Element::C),
            ..Default::default()
        };
        assert_eq!(AtomForm::from(Element::C), expected);
        assert_eq!(AtomForm::from("C"), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_element_form(AtomForm::default().with_element(ElementForm::Lit(Element::C)), AtomForm { element: ElementForm::Lit(Element::C), ..Default::default() })]
    #[case::with_element_primitive(AtomForm::default().with_element(Element::N), AtomForm { element: ElementForm::Lit(Element::N), ..Default::default() })]
    #[case::with_isotope_mass(AtomForm::default().with_isotope_mass(12_u32), AtomForm { isotope_mass: IsotopeMassForm::Lit(12), ..Default::default() })]
    #[case::with_charge(AtomForm::default().with_charge(1_i64), AtomForm { charge: NumForm::Lit(1), ..Default::default() })]
    #[case::with_implicit_hydrogens(AtomForm::default().with_implicit_hydrogens(3_i64), AtomForm { implicit_hydrogens: NumForm::Lit(3), ..Default::default() })]
    #[case::with_lone_pairs(AtomForm::default().with_lone_pairs(2_i64), AtomForm { lone_pairs: NumForm::Lit(2), ..Default::default() })]
    #[case::with_unpaired_electrons_tuple(AtomForm::default().with_unpaired_electrons((0_u8, 1_u8)), AtomForm { unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), ..Default::default() })]
    #[case::with_constraint(AtomForm::default().with_constraint(AtomConstraintForm::valence(4_i64)),
        AtomForm { constraints: AtomConstraintsForm::from(AtomConstraintForm::valence(4)),..Default::default() })]
    #[case::with_constraints_extends(AtomForm::default().with_constraint(AtomConstraintForm::valence(4_i64)).with_constraints([AtomConstraintForm::donated_pairs(1_i64), AtomConstraintForm::ring_membership(RingScope::Size(6), 1)]),
        AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::valence(4), AtomConstraintForm::donated_pairs(1), AtomConstraintForm::ring_membership(RingScope::Size(6), 1)]), ..Default::default() })]
    #[case::with_constraint_replaces_same_kind(AtomForm::default().with_constraint(AtomConstraintForm::valence(3_i64)).with_constraint(AtomConstraintForm::valence(4_i64)),
        AtomForm { constraints: AtomConstraintsForm::from(AtomConstraintForm::valence(4)), ..Default::default() })]
    fn test_atom_form_with_methods(#[case] actual: AtomForm, #[case] expected: AtomForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::element(AtomForm::from_element(Element::C), AtomUpdate { element: Some(ElementForm::Lit(Element::N)), ..Default::default() }, AtomForm::from_element(Element::N))]
    #[case::element_undetermined(AtomForm::from_element(Element::C), AtomUpdate { element: Some(ElementForm::Undetermined), ..Default::default() }, AtomForm::default())]
    #[case::isotope_mass(AtomForm::from_element(Element::C).with_isotope_mass(12_u32), AtomUpdate { isotope_mass: Some(IsotopeMassForm::Lit(13)), ..Default::default() }, AtomForm::from_element(Element::C).with_isotope_mass(13_u32))]
    #[case::isotope_mass_undetermined(AtomForm::from_element(Element::C).with_isotope_mass(12_u32), AtomUpdate { isotope_mass: Some(IsotopeMassForm::Undetermined), ..Default::default() }, AtomForm::from_element(Element::C))]
    #[case::charge(AtomForm::from_element(Element::C).with_charge(0_i64), AtomUpdate { charge: Some(NumForm::Lit(1)), ..Default::default() }, AtomForm::from_element(Element::C).with_charge(1_i64))]
    #[case::charge_undetermined(AtomForm::from_element(Element::C).with_charge(1_i64), AtomUpdate { charge: Some(NumForm::Undetermined), ..Default::default() }, AtomForm::from_element(Element::C))]
    #[case::implicit_hydrogens(AtomForm::from_element(Element::C).with_implicit_hydrogens(4_i64), AtomUpdate { implicit_hydrogens: Some(NumForm::Lit(3)), ..Default::default() }, AtomForm::from_element(Element::C).with_implicit_hydrogens(3_i64))]
    #[case::implicit_hydrogens_undetermined(AtomForm::from_element(Element::C).with_implicit_hydrogens(4_i64), AtomUpdate { implicit_hydrogens: Some(NumForm::Undetermined), ..Default::default() }, AtomForm::from_element(Element::C))]
    #[case::lone_pairs(AtomForm::from_element(Element::N).with_lone_pairs(1_i64), AtomUpdate { lone_pairs: Some(NumForm::Lit(2)), ..Default::default() }, AtomForm::from_element(Element::N).with_lone_pairs(2_i64))]
    #[case::lone_pairs_undetermined(AtomForm::from_element(Element::N).with_lone_pairs(1_i64), AtomUpdate { lone_pairs: Some(NumForm::Undetermined), ..Default::default() }, AtomForm::from_element(Element::N))]
    #[case::unpaired_electrons(AtomForm::from_element(Element::C).with_unpaired_electrons((2_u8, 3_u8)), AtomUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }, AtomForm::from_element(Element::C).with_unpaired_electrons((0_u8, 1_u8)))]
    #[case::unpaired_electrons_count(AtomForm::from_element(Element::C).with_unpaired_electrons((2_u8, 3_u8)), AtomUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: None }, ..Default::default() }, AtomForm::from_element(Element::C).with_unpaired_electrons((0_u8, 3_u8)))]
    #[case::unpaired_electrons_multiplicity(AtomForm::from_element(Element::C).with_unpaired_electrons((2_u8, 3_u8)), AtomUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }, AtomForm::from_element(Element::C).with_unpaired_electrons((2_u8, 1_u8)))]
    #[case::unpaired_electrons_count_undetermined(AtomForm::from_element(Element::C).with_unpaired_electrons((2_u8, 3_u8)), AtomUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: None }, ..Default::default() }, AtomForm::from_element(Element::C).with_unpaired_electrons(UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) }))]
    #[case::unpaired_electrons_multiplicity_undetermined(AtomForm::from_element(Element::C).with_unpaired_electrons((2_u8, 3_u8)), AtomUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Undetermined) }, ..Default::default() }, AtomForm::from_element(Element::C).with_unpaired_electrons(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined }))]
    #[case::constraint_set(AtomForm::from_element(Element::C), AtomUpdate { constraints: AtomConstraintsForm::from(AtomConstraintForm::valence(4_i64)), ..Default::default() }, AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4_i64)))]
    #[case::constraint_replace(AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(3_i64)), AtomUpdate { constraints: AtomConstraintsForm::from(AtomConstraintForm::valence(4_i64)), ..Default::default() }, AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4_i64)))]
    #[case::constraint_remove(AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4_i64)), AtomUpdate { constraints: AtomConstraintsForm::from(AtomConstraintForm::valence(NumForm::Undetermined)), ..Default::default() }, AtomForm::from_element(Element::C))]
    fn test_atom_form_update(#[case] atom: AtomForm, #[case] update: AtomUpdate, #[case] expected: AtomForm) {
        assert_eq!(atom.update(&update), expected);
    }

    #[rstest]
    #[case::empty(AtomForm::from_element(Element::C).with_charge(1_i64).with_constraint(AtomConstraintForm::valence(4_i64)))]
    fn test_atom_form_update_identity(#[case] atom: AtomForm) {
        assert_eq!(atom.update(&AtomUpdate::default()), atom);
    }

    #[rstest]
    fn test_atom_form_difference_to() {
        let atom = AtomForm::from_element(Element::C)
            .with_isotope_mass(12_u32)
            .with_charge(0_i64)
            .with_implicit_hydrogens(4_i64)
            .with_lone_pairs(0_i64)
            .with_unpaired_electrons((0_u8, 1_u8))
            .with_constraints([
                AtomConstraintForm::valence(4_i64),
                AtomConstraintForm::donated_pairs(1_i64),
            ]);
        let other = AtomForm::from_element(Element::N)
            .with_isotope_mass(13_u32)
            .with_implicit_hydrogens(3_i64)
            .with_lone_pairs(1_i64)
            .with_constraints([
                AtomConstraintForm::valence(3_i64),
                AtomConstraintForm::degree(2_i64),
            ]);
        assert_eq!(
            atom.difference_to(&other),
            AtomUpdate {
                element: Some(ElementForm::Lit(Element::N)),
                isotope_mass: Some(IsotopeMassForm::Lit(13)),
                charge: Some(NumForm::Undetermined),
                implicit_hydrogens: Some(NumForm::Lit(3)),
                lone_pairs: Some(NumForm::Lit(1)),
                unpaired_electrons: UnpairedElectronsUpdate {
                    count: Some(NumForm::Undetermined),
                    multiplicity: Some(NumForm::Undetermined),
                },
                constraints: AtomConstraintsForm::from_iter([
                    AtomConstraintForm::valence(3_i64),
                    AtomConstraintForm::donated_pairs(NumForm::Undetermined),
                    AtomConstraintForm::degree(2_i64),
                ]),
            }
        );
    }

    #[rstest]
    #[case::canonical(AtomForm::from_element(Element::C).with_charge(1_i64), AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([1])))]
    fn test_atom_form_difference_to_identity(#[case] atom: AtomForm, #[case] other: AtomForm) {
        assert_eq!(atom.difference_to(&other), AtomUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_element(AtomForm::from_element(Element::C).into_ground(),
        AtomForm { element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Natural, charge: NumForm::Lit(0), implicit_hydrogens: NumForm::Lit(0),
        lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), constraints: AtomConstraintsForm::new() })]
    #[case::with_charge(AtomForm::from_element(Element::C).with_charge(1_i64).into_ground(),
        AtomForm { element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Natural, charge: NumForm::Lit(1), implicit_hydrogens: NumForm::Lit(0),
        lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), constraints: AtomConstraintsForm::new() })]
    #[case::constraint(AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4_i64)).into_ground(),
        AtomForm { element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Natural, charge: NumForm::Lit(0), implicit_hydrogens: NumForm::Lit(0),
        lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), constraints: AtomConstraintsForm::from(AtomConstraintForm::valence(4)) })]
    fn test_atom_form_into_ground(#[case] actual: AtomForm, #[case] expected: AtomForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(AtomForm::default(), false)]
    #[case::all_ground(AtomForm { element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Lit(12), charge: NumForm::Lit(0),
        implicit_hydrogens: NumForm::Lit(4), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AtomConstraintsForm::new() }, true)]
    #[case::element_undetermined(AtomForm { element: ElementForm::Undetermined, isotope_mass: IsotopeMassForm::Lit(12), charge: NumForm::Lit(0),
        implicit_hydrogens: NumForm::Lit(4), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AtomConstraintsForm::new() }, false)]
    #[case::isotope_undetermined(AtomForm { element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Undetermined, charge: NumForm::Lit(0),
        implicit_hydrogens: NumForm::Lit(4), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AtomConstraintsForm::new() }, false)]
    #[case::charge_undetermined(AtomForm { element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Lit(12), charge: NumForm::Undetermined,
        implicit_hydrogens: NumForm::Lit(4), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AtomConstraintsForm::new() }, false)]
    #[case::hydrogens_undetermined(AtomForm { element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Lit(12), charge: NumForm::Lit(0),
        implicit_hydrogens: NumForm::Undetermined, lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AtomConstraintsForm::new() }, false)]
    #[case::lone_pairs_undetermined(AtomForm { element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Lit(12), charge: NumForm::Lit(0),
        implicit_hydrogens: NumForm::Lit(4), lone_pairs: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
        constraints: AtomConstraintsForm::new() }, false)]
    #[case::unpaired_electrons_undetermined(AtomForm { element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Lit(12), charge: NumForm::Lit(0),
        implicit_hydrogens: NumForm::Lit(4), lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::default(),
        constraints: AtomConstraintsForm::new() }, false)]
    fn test_atom_form_is_ground(#[case] form: AtomForm, #[case] expected: bool) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_charge(
        AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([4])),
        Ok(AtomForm::from_element(Element::C).with_charge(4)),
    )]
    #[case::charge_empty_litset_contradiction(
        AtomForm::from_element(Element::C).with_charge(NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_atom_form_canonicalize(
        #[case] input: AtomForm,
        #[case] expected: Result<AtomForm, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard_vs_ground(AtomForm::default(), AtomForm::from_element(Element::C), true)]
    #[case::same_element(AtomForm::from_element(Element::C), AtomForm::from_element(Element::C), true)]
    #[case::element_mismatch(AtomForm::from_element(Element::C), AtomForm::from_element(Element::N), false)]
    #[case::pattern_more_specific_than_target(AtomForm::from_element(Element::C), AtomForm::default(), false)]
    #[case::charge_mismatch(AtomForm::from_element(Element::C).with_charge(1_i64), AtomForm::from_element(Element::C).with_charge(0_i64), false)]
    #[case::charge_wildcard_pattern(AtomForm::from_element(Element::C), AtomForm::from_element(Element::C).with_charge(1_i64), true)]
    #[case::isotope_mismatch(AtomForm::from_element(Element::C).with_isotope_mass(12_u32), AtomForm::from_element(Element::C).with_isotope_mass(13_u32), false)]
    #[case::hydrogens_mismatch(AtomForm::from_element(Element::C).with_implicit_hydrogens(3_i64), AtomForm::from_element(Element::C).with_implicit_hydrogens(4_i64), false)]
    #[case::lone_pairs_mismatch(AtomForm::from_element(Element::C).with_lone_pairs(1_i64), AtomForm::from_element(Element::C).with_lone_pairs(2_i64), false)]
    #[case::unpaired_electrons_mismatch(AtomForm::from_element(Element::C).with_unpaired_electrons((2_u8, 3_u8)), AtomForm::from_element(Element::C).with_unpaired_electrons((0_u8, 1_u8)), false)]
    #[case::constraint_required_present(
        AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
        AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
        true)]
    #[case::constraint_required_absent(
        AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
        AtomForm::from_element(Element::C),
        false)]
    #[case::constraint_value_mismatch(
        AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
        AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(3)),
        false)]
    fn test_atom_form_matches(
        #[case] pattern: AtomForm,
        #[case] target: AtomForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(AtomForm::default(), AtomForm::default(), Some(AtomForm::default()))]
    #[case::element_mismatch(
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::N),
        None
    )]
    #[case::narrows_charge(AtomForm::from_element(Element::C), AtomForm::from_element(Element::C).with_charge(1),
        Some(AtomForm::from_element(Element::C).with_charge(1)))]
    fn test_atom_form_meet(
        #[case] a: AtomForm,
        #[case] b: AtomForm,
        #[case] expected: Option<AtomForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::element_mismatch_widens(AtomForm::from_element(Element::C), AtomForm::from_element(Element::N),
        ElementForm::lit_set(vec![Element::C, Element::N]))]
    fn test_atom_form_join_element(
        #[case] a: AtomForm,
        #[case] b: AtomForm,
        #[case] expected: ElementForm,
    ) {
        assert_eq!(a.join(&b).unwrap().element, expected);
    }

    #[rstest]
    #[case::charge_change(AtomForm::from_element(Element::C), AtomForm::from_element(Element::C).with_charge(1), true,
        AtomForm::from_element(Element::C).with_charge(1))]
    #[case::no_change(
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::C),
        false,
        AtomForm::from_element(Element::C)
    )]
    fn test_atom_form_narrow_from(
        #[case] mut target: AtomForm,
        #[case] source: AtomForm,
        #[case] expected_changed: bool,
        #[case] expected_after: AtomForm,
    ) {
        let changed = target.narrow_from(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_set(ElementForm::lit_set([Element::C, Element::N]), ElementForm::LitSet(Box::new(BTreeSet::from([Element::C, Element::N]))))]
    #[case::not(ElementForm::not(Element::H), ElementForm::NotSet(Box::new(BTreeSet::from([Element::H]))))]
    #[case::not_set(ElementForm::not_set([Element::F, Element::Cl]), ElementForm::NotSet(Box::new(BTreeSet::from([Element::F, Element::Cl]))))]
    #[case::var(ElementForm::var("x"), ElementForm::Var(Box::new(("x".to_string(), None))))]
    #[case::var_in(ElementForm::var_in("x", [Element::C]), ElementForm::Var(Box::new(("x".to_string(), Some((MemOp::In, BTreeSet::from([Element::C])))))))]
    #[case::var_not_in(ElementForm::var_not_in("x", [Element::C]), ElementForm::Var(Box::new(("x".to_string(), Some((MemOp::NotIn, BTreeSet::from([Element::C])))))))]
    fn test_element_form_constructors(#[case] actual: ElementForm, #[case] expected: ElementForm) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::carbon(Element::C, ElementForm::Lit(Element::C))]
    #[case::nitrogen(Element::N, ElementForm::Lit(Element::N))]
    fn test_element_form_from(#[case] element: Element, #[case] expected: ElementForm) {
        assert_eq!(ElementForm::from(element), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::litset_singleton(ElementForm::lit_set([Element::C]), Ok(ElementForm::Lit(Element::C)))]
    #[case::litset_empty(ElementForm::lit_set([]), Err(Contradiction))]
    #[case::notset_empty(ElementForm::not_set([]), Ok(ElementForm::Undetermined))]
    #[case::var_in_empty(ElementForm::var_in("x", []), Err(Contradiction))]
    #[case::var_not_in_vacuous(ElementForm::var_not_in("x", []), Ok(ElementForm::var("x")))]
    fn test_element_form_canonicalize(
        #[case] input: ElementForm,
        #[case] expected: Result<ElementForm, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ElementForm::Undetermined)]
    #[case::lit(ElementForm::Lit(Element::C))]
    #[case::litset(ElementForm::lit_set([Element::C, Element::N]))]
    #[case::notset(ElementForm::not(Element::H))]
    #[case::var_free(ElementForm::var("x"))]
    #[case::var_in(ElementForm::var_in("x", [Element::C]))]
    fn test_element_form_canonicalize_identity(#[case] input: ElementForm) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    /// Cardinality polarity + universe boundaries (need the 118-element universe,
    /// so expected sets are computed from `Element::all()`, not hardcoded).
    #[rstest]
    fn test_element_form_canonicalize_cardinality() {
        let take =
            |n: usize| -> BTreeSet<Element> { Element::all().iter().take(n).copied().collect() };
        let skip =
            |n: usize| -> BTreeSet<Element> { Element::all().iter().skip(n).copied().collect() };

        // Tiebreak: 59 stays positive; 60 flips to the complement.
        assert_eq!(
            ElementForm::LitSet(Box::new(take(59))).canonicalize(),
            Ok(ElementForm::LitSet(Box::new(take(59))))
        );
        assert_eq!(
            ElementForm::LitSet(Box::new(take(60))).canonicalize(),
            Ok(ElementForm::NotSet(Box::new(skip(60))))
        );
        // Full positive set → Undetermined; NotSet of the full set → Err.
        assert_eq!(
            ElementForm::LitSet(Box::new(take(118))).canonicalize(),
            Ok(ElementForm::Undetermined)
        );
        assert_eq!(
            ElementForm::NotSet(Box::new(take(118))).canonicalize(),
            Err(Contradiction)
        );
        // Large NotSet flips to a positive LitSet of its (small) complement.
        assert_eq!(
            ElementForm::NotSet(Box::new(take(60))).canonicalize(),
            Ok(ElementForm::LitSet(Box::new(skip(60))))
        );
        // Var In over a large domain flips to NotIn of the complement; full domain → free.
        assert_eq!(
            ElementForm::var_in("x", take(60)).canonicalize(),
            Ok(ElementForm::var_not_in("x", skip(60)))
        );
        assert_eq!(
            ElementForm::var_in("x", take(118)).canonicalize(),
            Ok(ElementForm::var("x"))
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_carbon(ElementForm::Lit(Element::C), Some(Element::C))]
    #[case::lit_nitrogen(ElementForm::Lit(Element::N), Some(Element::N))]
    #[case::undetermined(ElementForm::Undetermined, None)]
    #[case::litset(ElementForm::lit_set([Element::C, Element::N]), None)]
    #[case::notset(ElementForm::not(Element::H), None)]
    #[case::var_in(ElementForm::var_in("e", [Element::C]), None)]
    #[case::var(ElementForm::var("e"), None)]
    fn test_element_form_as_lit(#[case] form: ElementForm, #[case] expected: Option<Element>) {
        assert_eq!(form.as_lit(), expected);
        assert_eq!(form.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ElementForm::Lit(Element::C), false)]
    #[case::undetermined(ElementForm::Undetermined, true)]
    #[case::litset(ElementForm::lit_set([Element::C, Element::N]), false)]
    #[case::notset(ElementForm::not(Element::H), false)]
    #[case::var_in(ElementForm::var_in("e", [Element::C]), false)]
    #[case::var(ElementForm::var("e"), false)]
    fn test_element_form_is_undetermined(#[case] form: ElementForm, #[case] expected: bool) {
        assert_eq!(form.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ElementForm::Undetermined, ElementForm::Lit(Element::C), Some(ElementForm::Lit(Element::C)))]
    #[case::lit_und(ElementForm::Lit(Element::C), ElementForm::Undetermined, Some(ElementForm::Lit(Element::C)))]
    #[case::lit_lit_eq(ElementForm::Lit(Element::C), ElementForm::Lit(Element::C), Some(ElementForm::Lit(Element::C)))]
    #[case::lit_lit_neq(ElementForm::Lit(Element::C), ElementForm::Lit(Element::N), None)]
    #[case::lit_set_in(ElementForm::Lit(Element::C), ElementForm::lit_set([Element::C, Element::N]), Some(ElementForm::Lit(Element::C)))]
    #[case::lit_set_out(ElementForm::Lit(Element::O), ElementForm::lit_set([Element::C, Element::N]), None)]
    #[case::set_set_singleton(ElementForm::lit_set([Element::C, Element::N]), ElementForm::lit_set([Element::N, Element::O]), Some(ElementForm::Lit(Element::N)))]
    #[case::set_set_multi(ElementForm::lit_set([Element::C, Element::N, Element::O]), ElementForm::lit_set([Element::N, Element::O, Element::F]), Some(ElementForm::lit_set([Element::N, Element::O])))]
    #[case::set_set_disjoint(ElementForm::lit_set([Element::C, Element::N]), ElementForm::lit_set([Element::O, Element::F]), None)]
    #[case::set_notset(ElementForm::lit_set([Element::C, Element::N]), ElementForm::not(Element::N), Some(ElementForm::Lit(Element::C)))]
    #[case::notset_notset(ElementForm::not(Element::C), ElementForm::not(Element::N), Some(ElementForm::not_set([Element::C, Element::N])))]
    #[case::lit_notset_in(ElementForm::Lit(Element::C), ElementForm::not(Element::N), Some(ElementForm::Lit(Element::C)))]
    #[case::lit_notset_out(ElementForm::Lit(Element::C), ElementForm::not(Element::C), None)]
    #[case::var_var_eq(ElementForm::var("e"), ElementForm::var("e"), Some(ElementForm::var("e")))]
    #[case::var_var_neq(ElementForm::var("e"), ElementForm::var("f"), None)]
    #[case::var_lit(ElementForm::var("e"), ElementForm::Lit(Element::C), None)]
    #[case::lit_var(ElementForm::Lit(Element::C), ElementForm::var("e"), None)]
    fn test_element_form_meet(
        #[case] a: ElementForm,
        #[case] b: ElementForm,
        #[case] expected: Option<ElementForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ElementForm::Undetermined, ElementForm::Lit(Element::C), ElementForm::Undetermined)]
    #[case::lit_lit_eq(ElementForm::Lit(Element::C), ElementForm::Lit(Element::C), ElementForm::Lit(Element::C))]
    #[case::lit_lit_neq(ElementForm::Lit(Element::C), ElementForm::Lit(Element::N), ElementForm::lit_set([Element::C, Element::N]))]
    #[case::lit_set(ElementForm::Lit(Element::O), ElementForm::lit_set([Element::C, Element::N]), ElementForm::lit_set([Element::C, Element::N, Element::O]))]
    #[case::set_set(ElementForm::lit_set([Element::C, Element::N]), ElementForm::lit_set([Element::N, Element::O]), ElementForm::lit_set([Element::C, Element::N, Element::O]))]
    #[case::lit_notset_out(ElementForm::Lit(Element::C), ElementForm::not(Element::N), ElementForm::not(Element::N))]
    #[case::lit_notset_in(ElementForm::Lit(Element::N), ElementForm::not(Element::N), ElementForm::Undetermined)]
    #[case::notset_notset_disjoint(ElementForm::not(Element::C), ElementForm::not(Element::N), ElementForm::Undetermined)]
    #[case::notset_notset_overlap(ElementForm::not_set([Element::C, Element::N]), ElementForm::not_set([Element::N, Element::O]), ElementForm::not(Element::N))]
    #[case::var_var_eq(ElementForm::var("e"), ElementForm::var("e"), ElementForm::var("e"))]
    #[case::var_var_neq(ElementForm::var("e"), ElementForm::var("f"), ElementForm::Undetermined)]
    #[case::var_lit(ElementForm::var("e"), ElementForm::Lit(Element::C), ElementForm::Undetermined)]
    fn test_element_form_join(
        #[case] a: ElementForm,
        #[case] b: ElementForm,
        #[case] expected: ElementForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(ElementForm::Undetermined, ElementForm::Lit(Element::C), true)]
    #[case::und_und(ElementForm::Undetermined, ElementForm::Undetermined, true)]
    #[case::und_set(ElementForm::Undetermined, ElementForm::lit_set([Element::C, Element::N]), true)]
    #[case::und_var(ElementForm::Undetermined, ElementForm::var("e"), true)]
    #[case::lit_und(ElementForm::Lit(Element::C), ElementForm::Undetermined, false)]
    #[case::set_und(ElementForm::lit_set([Element::C, Element::N]), ElementForm::Undetermined, false)]
    #[case::var_und(ElementForm::var("e"), ElementForm::Undetermined, false)]
    #[case::lit_lit_match(ElementForm::Lit(Element::C), ElementForm::Lit(Element::C), true)]
    #[case::lit_lit_mismatch(ElementForm::Lit(Element::C), ElementForm::Lit(Element::N), false)]
    #[case::lit_singleton_set(ElementForm::Lit(Element::C), ElementForm::lit_set([Element::C]), true)]
    #[case::lit_multi_set(ElementForm::Lit(Element::C), ElementForm::lit_set([Element::C, Element::N]), false)]
    #[case::set_lit_in(ElementForm::lit_set([Element::C, Element::N]), ElementForm::Lit(Element::N), true)]
    #[case::set_lit_out(ElementForm::lit_set([Element::C, Element::N]), ElementForm::Lit(Element::O), false)]
    #[case::set_set_subset(ElementForm::lit_set([Element::C, Element::N, Element::O]), ElementForm::lit_set([Element::C, Element::N]), true)]
    #[case::set_set_equal(ElementForm::lit_set([Element::C, Element::N]), ElementForm::lit_set([Element::C, Element::N]), true)]
    #[case::set_set_superset(ElementForm::lit_set([Element::C]), ElementForm::lit_set([Element::C, Element::N]), false)]
    #[case::notset_lit_admitted(ElementForm::not(Element::N), ElementForm::Lit(Element::C), true)]
    #[case::notset_lit_excluded(ElementForm::not(Element::C), ElementForm::Lit(Element::C), false)]
    #[case::var_var_equal(ElementForm::var("e"), ElementForm::var("e"), true)]
    #[case::var_var_distinct(ElementForm::var("e"), ElementForm::var("f"), false)]
    #[case::var_lit(ElementForm::var("e"), ElementForm::Lit(Element::C), false)]
    #[case::lit_var(ElementForm::Lit(Element::C), ElementForm::var("e"), false)]
    #[case::set_var(ElementForm::lit_set([Element::C]), ElementForm::var("e"), false)]
    fn test_element_form_matches(
        #[case] pattern: ElementForm,
        #[case] target: ElementForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_set(IsotopeMassForm::lit_set([12, 13]), IsotopeMassForm::LitSet(Box::new(BTreeSet::from([12, 13]))))]
    #[case::var(IsotopeMassForm::var("m"), IsotopeMassForm::Var(Box::new(("m".to_string(), None))))]
    #[case::var_in(IsotopeMassForm::var_in("m", [12, 13]), IsotopeMassForm::Var(Box::new(("m".to_string(), Some(BTreeSet::from([12, 13]))))))]
    fn test_isotope_mass_form_constructors(#[case] actual: IsotopeMassForm, #[case] expected: IsotopeMassForm) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::positive(13, IsotopeMassForm::Lit(13))]
    #[case::zero(0, IsotopeMassForm::Lit(0))]
    fn test_isotope_mass_form_from(#[case] mass: u32, #[case] expected: IsotopeMassForm) {
        assert_eq!(IsotopeMassForm::from(mass), expected);
    }

    #[rstest]
    #[case::natural(IsotopeMass::Natural, IsotopeMassForm::Natural)]
    #[case::mass_number(IsotopeMass::MassNumber(13), IsotopeMassForm::Lit(13))]
    fn test_isotope_mass_form_from_isotope_mass(
        #[case] mass: IsotopeMass,
        #[case] expected: IsotopeMassForm,
    ) {
        assert_eq!(IsotopeMassForm::from(mass), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::litset_singleton(IsotopeMassForm::lit_set([12]), Ok(IsotopeMassForm::Lit(12)))]
    #[case::litset_empty(IsotopeMassForm::lit_set([]), Err(Contradiction))]
    #[case::var_in_empty(IsotopeMassForm::var_in("m", []), Err(Contradiction))]
    fn test_isotope_mass_form_canonicalize(
        #[case] input: IsotopeMassForm,
        #[case] expected: Result<IsotopeMassForm, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(IsotopeMassForm::Undetermined)]
    #[case::natural(IsotopeMassForm::Natural)]
    #[case::lit(IsotopeMassForm::Lit(12))]
    #[case::litset(IsotopeMassForm::lit_set([12, 13]))]
    #[case::var_free(IsotopeMassForm::var("m"))]
    #[case::var_in(IsotopeMassForm::var_in("m", [12, 13]))]
    #[case::var_in_singleton(IsotopeMassForm::var_in("m", [12]))]
    fn test_isotope_mass_form_canonicalize_identity(#[case] input: IsotopeMassForm) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(IsotopeMassForm::Lit(12), Some(IsotopeMass::MassNumber(12)))]
    #[case::lit_zero(IsotopeMassForm::Lit(0), Some(IsotopeMass::MassNumber(0)))]
    #[case::natural(IsotopeMassForm::Natural, Some(IsotopeMass::Natural))]
    #[case::undetermined(IsotopeMassForm::Undetermined, None)]
    #[case::litset(IsotopeMassForm::lit_set([12, 13]), None)]
    #[case::var(IsotopeMassForm::var("m"), None)]
    #[case::var_in(IsotopeMassForm::var_in("m", [12]), None)]
    fn test_isotope_mass_form_as_lit(
        #[case] form: IsotopeMassForm,
        #[case] expected: Option<IsotopeMass>,
    ) {
        assert_eq!(form.as_lit(), expected);
        assert_eq!(form.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(IsotopeMassForm::Undetermined, true)]
    #[case::natural(IsotopeMassForm::Natural, false)]
    #[case::lit(IsotopeMassForm::Lit(12), false)]
    #[case::litset(IsotopeMassForm::lit_set([12, 13]), false)]
    #[case::var(IsotopeMassForm::var("m"), false)]
    #[case::var_in(IsotopeMassForm::var_in("m", [12]), false)]
    fn test_isotope_mass_form_is_undetermined(#[case] form: IsotopeMassForm, #[case] expected: bool) {
        assert_eq!(form.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::natural(IsotopeMassForm::Natural, true)]
    #[case::lit(IsotopeMassForm::Lit(12), true)]
    #[case::undetermined(IsotopeMassForm::Undetermined, false)]
    #[case::litset(IsotopeMassForm::lit_set([12, 13]), false)]
    #[case::var(IsotopeMassForm::var("m"), false)]
    #[case::var_in(IsotopeMassForm::var_in("m", [12]), false)]
    fn test_isotope_mass_form_is_ground(#[case] form: IsotopeMassForm, #[case] expected: bool) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(IsotopeMassForm::Undetermined, IsotopeMassForm::Lit(12), Some(IsotopeMassForm::Lit(12)))]
    #[case::lit_und(IsotopeMassForm::Lit(12), IsotopeMassForm::Undetermined, Some(IsotopeMassForm::Lit(12)))]
    #[case::und_natural(IsotopeMassForm::Undetermined, IsotopeMassForm::Natural, Some(IsotopeMassForm::Natural))]
    #[case::natural_natural(IsotopeMassForm::Natural, IsotopeMassForm::Natural, Some(IsotopeMassForm::Natural))]
    #[case::natural_lit(IsotopeMassForm::Natural, IsotopeMassForm::Lit(12), None)]
    #[case::lit_natural(IsotopeMassForm::Lit(12), IsotopeMassForm::Natural, None)]
    #[case::lit_lit_eq(IsotopeMassForm::Lit(12), IsotopeMassForm::Lit(12), Some(IsotopeMassForm::Lit(12)))]
    #[case::lit_lit_neq(IsotopeMassForm::Lit(12), IsotopeMassForm::Lit(13), None)]
    #[case::lit_set_in(IsotopeMassForm::Lit(12), IsotopeMassForm::lit_set([12, 13]), Some(IsotopeMassForm::Lit(12)))]
    #[case::lit_set_out(IsotopeMassForm::Lit(14), IsotopeMassForm::lit_set([12, 13]), None)]
    #[case::set_set_singleton(IsotopeMassForm::lit_set([12, 13]), IsotopeMassForm::lit_set([13, 14]), Some(IsotopeMassForm::Lit(13)))]
    #[case::set_set_multi(IsotopeMassForm::lit_set([12, 13, 14]), IsotopeMassForm::lit_set([13, 14, 15]), Some(IsotopeMassForm::lit_set([13, 14])))]
    #[case::set_set_disjoint(IsotopeMassForm::lit_set([12, 13]), IsotopeMassForm::lit_set([14, 15]), None)]
    #[case::var_var_eq(IsotopeMassForm::var("m"), IsotopeMassForm::var("m"), Some(IsotopeMassForm::var("m")))]
    #[case::var_var_neq(IsotopeMassForm::var("m"), IsotopeMassForm::var("n"), None)]
    #[case::var_lit(IsotopeMassForm::var("m"), IsotopeMassForm::Lit(12), None)]
    fn test_isotope_mass_form_meet(
        #[case] a: IsotopeMassForm,
        #[case] b: IsotopeMassForm,
        #[case] expected: Option<IsotopeMassForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(IsotopeMassForm::Undetermined, IsotopeMassForm::Lit(12), IsotopeMassForm::Undetermined)]
    #[case::natural_natural(IsotopeMassForm::Natural, IsotopeMassForm::Natural, IsotopeMassForm::Natural)]
    #[case::natural_lit(IsotopeMassForm::Natural, IsotopeMassForm::Lit(12), IsotopeMassForm::Undetermined)]
    #[case::lit_lit_eq(IsotopeMassForm::Lit(12), IsotopeMassForm::Lit(12), IsotopeMassForm::Lit(12))]
    #[case::lit_lit_neq(IsotopeMassForm::Lit(12), IsotopeMassForm::Lit(13), IsotopeMassForm::lit_set([12, 13]))]
    #[case::lit_set(IsotopeMassForm::Lit(14), IsotopeMassForm::lit_set([12, 13]), IsotopeMassForm::lit_set([12, 13, 14]))]
    #[case::set_set(IsotopeMassForm::lit_set([12, 13]), IsotopeMassForm::lit_set([13, 14]), IsotopeMassForm::lit_set([12, 13, 14]))]
    #[case::var_var_eq(IsotopeMassForm::var("m"), IsotopeMassForm::var("m"), IsotopeMassForm::var("m"))]
    #[case::var_var_neq(IsotopeMassForm::var("m"), IsotopeMassForm::var("n"), IsotopeMassForm::Undetermined)]
    #[case::var_lit(IsotopeMassForm::var("m"), IsotopeMassForm::Lit(12), IsotopeMassForm::Undetermined)]
    fn test_isotope_mass_form_join(
        #[case] a: IsotopeMassForm,
        #[case] b: IsotopeMassForm,
        #[case] expected: IsotopeMassForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_lit(IsotopeMassForm::Undetermined, IsotopeMassForm::Lit(12), true)]
    #[case::und_natural(IsotopeMassForm::Undetermined, IsotopeMassForm::Natural, true)]
    #[case::und_und(IsotopeMassForm::Undetermined, IsotopeMassForm::Undetermined, true)]
    #[case::und_set(IsotopeMassForm::Undetermined, IsotopeMassForm::lit_set([12, 13]), true)]
    #[case::und_var(IsotopeMassForm::Undetermined, IsotopeMassForm::var("m"), true)]
    #[case::lit_und(IsotopeMassForm::Lit(12), IsotopeMassForm::Undetermined, false)]
    #[case::natural_und(IsotopeMassForm::Natural, IsotopeMassForm::Undetermined, false)]
    #[case::natural_natural(IsotopeMassForm::Natural, IsotopeMassForm::Natural, true)]
    #[case::natural_lit(IsotopeMassForm::Natural, IsotopeMassForm::Lit(12), false)]
    #[case::lit_natural(IsotopeMassForm::Lit(12), IsotopeMassForm::Natural, false)]
    #[case::lit_lit_match(IsotopeMassForm::Lit(12), IsotopeMassForm::Lit(12), true)]
    #[case::lit_lit_mismatch(IsotopeMassForm::Lit(12), IsotopeMassForm::Lit(13), false)]
    #[case::set_lit_in(IsotopeMassForm::lit_set([12, 13]), IsotopeMassForm::Lit(13), true)]
    #[case::set_lit_out(IsotopeMassForm::lit_set([12, 13]), IsotopeMassForm::Lit(14), false)]
    #[case::set_set_subset(IsotopeMassForm::lit_set([12, 13, 14]), IsotopeMassForm::lit_set([12, 13]), true)]
    #[case::set_set_superset(IsotopeMassForm::lit_set([12]), IsotopeMassForm::lit_set([12, 13]), false)]
    #[case::var_var_equal(IsotopeMassForm::var("m"), IsotopeMassForm::var("m"), true)]
    #[case::var_var_distinct(IsotopeMassForm::var("m"), IsotopeMassForm::var("n"), false)]
    #[case::var_lit(IsotopeMassForm::var("m"), IsotopeMassForm::Lit(12), false)]
    #[case::lit_var(IsotopeMassForm::Lit(12), IsotopeMassForm::var("m"), false)]
    fn test_isotope_mass_form_matches(
        #[case] pattern: IsotopeMassForm,
        #[case] target: IsotopeMassForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
