//! Atom-level AST fragments shared across crates.

use std::{mem, slice};

use umol_ast_macros::Lattice;
use umol_shared::element::Element;

use super::constraint::joint_domain::{JointDomainAst, JointValue, JointVar};
use super::constraint::{
    AromaticValenceAst, AtomConstraint, AtomConstraintKind, AtomConstraints, MulticenterValenceAst,
};
use super::error::Contradiction;
use super::operators::MemOp;
use super::spin::SpinStateAst;
use super::stereo::StereoConfigurationAst;
use super::traits::{AsLit, Lattice};
use super::value::{set_is_ground, ValueAst};

/// Atom AST: structural representation of an atom plus the atom-level
/// constraints (valence, degree, ring membership, etc.) that pattern
/// against the surrounding topology.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Lattice)]
#[lattice(saturate = "saturate_atom")]
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

/// Cross-field saturation for `AtomAst`. Walks every `JointDomain`
/// constraint, projects its tuples against current field values (forward
/// propagation only), and either keeps the pruned JD, narrows fields to
/// the surviving tuple's literals and drops the JD (single tuple remains),
/// or reports `Err(Contradiction)` (no tuples remain). Narrowing is
/// inline — subsequent JDs in the same pass observe the just-narrowed
/// fields. Loops to a fixpoint so that cascading narrowings (one JD
/// narrowing a field that another JD reads) are resolved.
///
/// Wired into `AtomAst`'s derived `Lattice::meet` via
/// `#[lattice(saturate = "saturate_atom")]`; every `meet` call runs this
/// after the field-wise meet completes.
pub fn saturate_atom(atom: &mut AtomAst) -> Result<(), Contradiction> {
    loop {
        let jds: Vec<JointDomainAst> = atom.constraints.joint_domains().cloned().collect();
        atom.constraints.remove_all(AtomConstraintKind::JointDomain);
        let mut changed = false;
        for jd in jds {
            let pruned = jd.project(|var| field_value_for_joint_var(atom, var))?;
            let Some(tuples) = pruned.tuples() else {
                continue;
            };
            if tuples.len() == 1 {
                let vars = pruned.vars().expect("Domain has vars");
                let tuple = tuples[0].clone();
                for (var, value) in vars.iter().zip(tuple.iter()) {
                    narrow_joint_var_to_lit(atom, *var, value);
                }
                changed = true;
            } else {
                if pruned != jd {
                    changed = true;
                }
                atom.constraints.add(AtomConstraint::JointDomain(pruned));
            }
        }
        if !changed {
            return Ok(());
        }
    }
}

/// Read the atom field corresponding to a `JointVar`. Struct fields
/// (`charge`, etc.) read directly; constraint-backed vars (`Valence`,
/// `DonatedPairs`, `AcceptedPairs`) read through the constraint
/// container's typed accessors, which return `Undetermined` when absent.
fn field_value_for_joint_var(atom: &AtomAst, var: JointVar) -> ValueAst {
    match var {
        JointVar::Charge => atom.charge.clone(),
        JointVar::ImplicitHydrogens => atom.implicit_hydrogens.clone(),
        JointVar::LonePairs => atom.lone_pairs.clone(),
        JointVar::UnpairedElectrons => atom.spin.unpaired.clone(),
        JointVar::Multiplicity => atom.spin.multiplicity.clone(),
        JointVar::Valence => atom.constraints.valence(),
        JointVar::DonatedPairs => atom.constraints.donated_pairs(),
        JointVar::AcceptedPairs => atom.constraints.accepted_pairs(),
    }
}

/// Narrow the atom field corresponding to `var` to `ValueAst::Lit(value)`.
/// For struct fields, calls `narrow_from`. For constraint-backed vars,
/// `add`s a `Lit(value)` constraint (last-wins per `AtomConstraints::add`
/// for unique-kind variants).
fn narrow_joint_var_to_lit(atom: &mut AtomAst, var: JointVar, value: &JointValue) {
    let JointValue::Int(n) = value;
    let lit = ValueAst::Lit(*n);
    match var {
        JointVar::Charge => {
            atom.charge.narrow_from(&lit);
        }
        JointVar::ImplicitHydrogens => {
            atom.implicit_hydrogens.narrow_from(&lit);
        }
        JointVar::LonePairs => {
            atom.lone_pairs.narrow_from(&lit);
        }
        JointVar::UnpairedElectrons => {
            atom.spin.unpaired.narrow_from(&lit);
        }
        JointVar::Multiplicity => {
            atom.spin.multiplicity.narrow_from(&lit);
        }
        JointVar::Valence => {
            atom.constraints.add(AtomConstraint::Valence(lit));
        }
        JointVar::DonatedPairs => {
            atom.constraints.add(AtomConstraint::DonatedPairs(lit));
        }
        JointVar::AcceptedPairs => {
            atom.constraints.add(AtomConstraint::AcceptedPairs(lit));
        }
    }
}

/// Element expressions
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ElementAst {
    #[default]
    Undetermined,
    Lit(Element),
    Set(Vec<Element>),
    Not(Element),
    NotSet(Vec<Element>),
    Bind(Box<(String, MemOp, Vec<Element>)>),
    Ref(String),
}

impl ElementAst {
    pub fn bind(id: impl Into<String>, set: Vec<Element>, op: MemOp) -> Self {
        Self::Bind(Box::new((id.into(), op, set)))
    }
}

impl From<Element> for ElementAst {
    fn from(element: Element) -> Self {
        Self::Lit(element)
    }
}

impl AsLit for ElementAst {
    type Lit = Element;

    #[inline]
    fn as_lit(&self) -> Option<Element> {
        match self {
            Self::Lit(e) => Some(*e),
            Self::Undetermined | Self::Ref(_) | Self::Bind(_) | Self::Not(_) | Self::NotSet(_) => {
                None
            }
            Self::Set(s) => element_set_is_ground(s).then(|| s[0]),
        }
    }
}

impl Lattice for ElementAst {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        match self {
            Self::Lit(_) => true,
            Self::Undetermined | Self::Ref(_) | Self::Bind(_) | Self::Not(_) | Self::NotSet(_) => {
                false
            }
            Self::Set(s) => element_set_is_ground(s),
        }
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::Lit(a), Self::Lit(b)) => (a == b).then_some(Self::Lit(*a)),
            (Self::Lit(a), Self::Set(s)) | (Self::Set(s), Self::Lit(a)) => {
                s.contains(a).then_some(Self::Lit(*a))
            }
            (Self::Lit(a), Self::Not(b)) | (Self::Not(b), Self::Lit(a)) => {
                (a != b).then_some(Self::Lit(*a))
            }
            (Self::Lit(a), Self::NotSet(s)) | (Self::NotSet(s), Self::Lit(a)) => {
                (!s.contains(a)).then_some(Self::Lit(*a))
            }
            (Self::Set(s), Self::Set(t)) => {
                let intersection: Vec<Element> =
                    s.iter().filter(|x| t.contains(x)).copied().collect();
                canonicalize_set(intersection)
            }
            (Self::Set(s), Self::Not(b)) | (Self::Not(b), Self::Set(s)) => {
                let filtered: Vec<Element> = s.iter().filter(|x| x != &b).copied().collect();
                canonicalize_set(filtered)
            }
            (Self::Set(s), Self::NotSet(t)) | (Self::NotSet(t), Self::Set(s)) => {
                let filtered: Vec<Element> = s.iter().filter(|x| !t.contains(x)).copied().collect();
                canonicalize_set(filtered)
            }
            (Self::Not(a), Self::Not(b)) => {
                if a == b {
                    Some(Self::Not(*a))
                } else {
                    let mut v = vec![*a];
                    if !v.contains(b) {
                        v.push(*b);
                    }
                    Some(Self::NotSet(v))
                }
            }
            (Self::Not(a), Self::NotSet(s)) | (Self::NotSet(s), Self::Not(a)) => {
                let mut v: Vec<Element> = s.clone();
                if !v.contains(a) {
                    v.push(*a);
                }
                Some(Self::NotSet(v))
            }
            (Self::NotSet(s), Self::NotSet(t)) => {
                let mut v: Vec<Element> = s.clone();
                for &x in t.iter() {
                    if !v.contains(&x) {
                        v.push(x);
                    }
                }
                Some(canonicalize_not_set(v))
            }
            (Self::Ref(a), Self::Ref(b)) if a == b => Some(Self::Ref(a.clone())),
            (Self::Bind(a), Self::Bind(b)) if a == b => Some(self.clone()),
            _ => None,
        }
    }

    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Lit(a), Self::Lit(b)) => {
                if a == b {
                    Self::Lit(*a)
                } else {
                    Self::Set(vec![*a, *b])
                }
            }
            (Self::Lit(a), Self::Set(s)) | (Self::Set(s), Self::Lit(a)) => {
                let mut v: Vec<Element> = s.clone();
                if !v.contains(a) {
                    v.push(*a);
                }
                canonicalize_set(v).unwrap_or(Self::Undetermined)
            }
            (Self::Set(s), Self::Set(t)) => {
                let mut v: Vec<Element> = s.clone();
                for &x in t.iter() {
                    if !v.contains(&x) {
                        v.push(x);
                    }
                }
                canonicalize_set(v).unwrap_or(Self::Undetermined)
            }
            (Self::Lit(a), Self::Not(b)) | (Self::Not(b), Self::Lit(a)) => {
                if a == b {
                    Self::Undetermined
                } else {
                    Self::Not(*b)
                }
            }
            (Self::Lit(a), Self::NotSet(s)) | (Self::NotSet(s), Self::Lit(a)) => {
                let remaining: Vec<Element> = s.iter().filter(|x| x != &a).copied().collect();
                canonicalize_not_set(remaining)
            }
            (Self::Set(s), Self::Not(b)) | (Self::Not(b), Self::Set(s)) => {
                if s.contains(b) {
                    Self::Undetermined
                } else {
                    Self::Not(*b)
                }
            }
            (Self::Set(s), Self::NotSet(t)) | (Self::NotSet(t), Self::Set(s)) => {
                let remaining: Vec<Element> =
                    t.iter().filter(|x| !s.contains(x)).copied().collect();
                canonicalize_not_set(remaining)
            }
            (Self::Not(a), Self::Not(b)) => {
                if a == b {
                    Self::Not(*a)
                } else {
                    Self::Undetermined
                }
            }
            (Self::Not(a), Self::NotSet(s)) | (Self::NotSet(s), Self::Not(a)) => {
                if s.contains(a) {
                    Self::Not(*a)
                } else {
                    Self::Undetermined
                }
            }
            (Self::NotSet(s), Self::NotSet(t)) => {
                let intersection: Vec<Element> =
                    s.iter().filter(|x| t.contains(x)).copied().collect();
                canonicalize_not_set(intersection)
            }
            (Self::Ref(a), Self::Ref(b)) if a == b => Self::Ref(a.clone()),
            (Self::Bind(a), Self::Bind(b)) if a == b => self.clone(),
            _ => Self::Undetermined,
        }
    }

    /// Pattern matches target iff every element the target admits is also
    /// admitted by the pattern (superset semantics).
    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Ref(_), _) | (_, Self::Ref(_)) => false,
            _ => {
                let (Some((ps, pp)), Some((ts, tp))) =
                    (element_set_view(self), element_set_view(target))
                else {
                    return false;
                };
                match (pp, tp) {
                    (MemOp::In, MemOp::In) => ts.iter().all(|t| ps.contains(t)),
                    (MemOp::In, MemOp::NotIn) => false,
                    (MemOp::NotIn, MemOp::In) => ts.iter().all(|t| !ps.contains(t)),
                    (MemOp::NotIn, MemOp::NotIn) => ps.iter().all(|p| ts.contains(p)),
                }
            }
        }
    }
}

fn element_set_is_ground(s: &[Element]) -> bool {
    match s {
        [] => false,
        [first, rest @ ..] => rest.iter().all(|x| x == first),
    }
}

/// Canonicalize a literal-set intersection result: empty → None (contradiction),
/// singleton → Lit, otherwise → Set. Used by meet operations.
fn canonicalize_set(v: Vec<Element>) -> Option<ElementAst> {
    match v.len() {
        0 => None,
        1 => Some(ElementAst::Lit(v[0])),
        _ => Some(ElementAst::Set(v)),
    }
}

/// Canonicalize a not-set result: empty → Undetermined, singleton → Not,
/// otherwise → NotSet.
fn canonicalize_not_set(v: Vec<Element>) -> ElementAst {
    match v.len() {
        0 => ElementAst::Undetermined,
        1 => ElementAst::Not(v[0]),
        _ => ElementAst::NotSet(v),
    }
}

/// View any concrete `ElementAst` (Lit / Set / Not / NotSet / Bind) as a
/// `(set, polarity)` pair describing the admissible domain. Returns None for
/// `Undetermined` and `Ref` which don't have a finite set encoding.
fn element_set_view(ast: &ElementAst) -> Option<(&[Element], MemOp)> {
    match ast {
        ElementAst::Lit(x) => Some((slice::from_ref(x), MemOp::In)),
        ElementAst::Set(s) => Some((s, MemOp::In)),
        ElementAst::Not(x) => Some((slice::from_ref(x), MemOp::NotIn)),
        ElementAst::NotSet(s) => Some((s, MemOp::NotIn)),
        ElementAst::Bind(b) => Some((&b.2, b.1)),
        ElementAst::Undetermined | ElementAst::Ref(_) => None,
    }
}

/// Isotope-mass expressions. `Natural` denotes the naturally most abundant
/// isotope (`#i=`) and is its own channel — does not meet with numeric
/// literals. `Not`/`NotSet` are cofinite exclusion forms (`#i!12`,
/// `#i!{12,13}`). `Bind`/`Ref` are named-bind / named-reference for
/// joint-domain constraints, mirroring `ElementAst`.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum IsotopeMassAst {
    #[default]
    Undetermined,
    Natural,
    Lit(i64),
    Set(Box<Vec<i64>>),
    Not(i64),
    NotSet(Box<Vec<i64>>),
    Bind(Box<(String, MemOp, Vec<i64>)>),
    Ref(String),
}

impl IsotopeMassAst {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn natural() -> Self {
        Self::Natural
    }

    pub fn lit(n: i64) -> Self {
        Self::Lit(n)
    }

    pub fn set(values: Vec<i64>) -> Self {
        Self::Set(Box::new(values))
    }

    pub fn not(n: i64) -> Self {
        Self::Not(n)
    }

    pub fn not_set(values: Vec<i64>) -> Self {
        Self::NotSet(Box::new(values))
    }

    pub fn bind(id: impl Into<String>, set: Vec<i64>, op: MemOp) -> Self {
        Self::Bind(Box::new((id.into(), op, set)))
    }

    pub fn reference(id: impl Into<String>) -> Self {
        Self::Ref(id.into())
    }
}

impl From<i64> for IsotopeMassAst {
    fn from(value: i64) -> Self {
        Self::Lit(value)
    }
}

impl AsLit for IsotopeMassAst {
    type Lit = u32;

    /// Mass number when ground; `None` otherwise. `Natural` returns
    /// `Some(0)` as the sentinel for "natural isotopic abundance — no
    /// specific mass committed".
    #[inline]
    fn as_lit(&self) -> Option<u32> {
        match self {
            Self::Natural => Some(0),
            Self::Lit(n) => u32::try_from(*n).ok(),
            Self::Set(s) => set_is_ground(s).then(|| u32::try_from(s[0]).ok()).flatten(),
            Self::Undetermined | Self::Not(_) | Self::NotSet(_) | Self::Bind(_) | Self::Ref(_) => {
                None
            }
        }
    }
}

impl Lattice for IsotopeMassAst {
    #[inline]
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    /// `Natural` and singleton `Lit`/`Set` are ground; cofinite forms,
    /// `Bind`, `Ref`, and `Undetermined` are not.
    #[inline]
    fn is_ground(&self) -> bool {
        match self {
            Self::Natural | Self::Lit(_) => true,
            Self::Set(s) => set_is_ground(s),
            Self::Undetermined | Self::Not(_) | Self::NotSet(_) | Self::Bind(_) | Self::Ref(_) => {
                false
            }
        }
    }

    /// `Natural` is its own channel — meet with any numeric form yields
    /// `None`. `Undetermined` is the top. Other combinations follow the
    /// finite-Boolean algebra extended with bind/ref equality.
    fn meet(&self, other: &Self) -> Option<Self> {
        match (self, other) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::Natural, Self::Natural) => Some(Self::Natural),
            (Self::Natural, _) | (_, Self::Natural) => None,
            (Self::Lit(a), Self::Lit(b)) => (a == b).then_some(Self::Lit(*a)),
            (Self::Lit(a), Self::Set(s)) | (Self::Set(s), Self::Lit(a)) => {
                s.contains(a).then_some(Self::Lit(*a))
            }
            (Self::Lit(a), Self::Not(b)) | (Self::Not(b), Self::Lit(a)) => {
                (a != b).then_some(Self::Lit(*a))
            }
            (Self::Lit(a), Self::NotSet(s)) | (Self::NotSet(s), Self::Lit(a)) => {
                (!s.contains(a)).then_some(Self::Lit(*a))
            }
            (Self::Set(s), Self::Set(t)) => {
                let v: Vec<i64> = s.iter().filter(|x| t.contains(x)).copied().collect();
                canonicalize_isotope_set(v)
            }
            (Self::Set(s), Self::Not(b)) | (Self::Not(b), Self::Set(s)) => {
                let v: Vec<i64> = s.iter().filter(|x| x != &b).copied().collect();
                canonicalize_isotope_set(v)
            }
            (Self::Set(s), Self::NotSet(t)) | (Self::NotSet(t), Self::Set(s)) => {
                let v: Vec<i64> = s.iter().filter(|x| !t.contains(x)).copied().collect();
                canonicalize_isotope_set(v)
            }
            (Self::Not(a), Self::Not(b)) => {
                if a == b {
                    Some(Self::Not(*a))
                } else {
                    let mut v = vec![*a];
                    if !v.contains(b) {
                        v.push(*b);
                    }
                    Some(Self::NotSet(Box::new(v)))
                }
            }
            (Self::Not(a), Self::NotSet(s)) | (Self::NotSet(s), Self::Not(a)) => {
                let mut v: Vec<i64> = (**s).clone();
                if !v.contains(a) {
                    v.push(*a);
                }
                Some(Self::NotSet(Box::new(v)))
            }
            (Self::NotSet(s), Self::NotSet(t)) => {
                let mut v: Vec<i64> = (**s).clone();
                for &x in t.iter() {
                    if !v.contains(&x) {
                        v.push(x);
                    }
                }
                Some(canonicalize_isotope_not_set(v))
            }
            (Self::Ref(a), Self::Ref(b)) if a == b => Some(Self::Ref(a.clone())),
            (Self::Bind(a), Self::Bind(b)) if a == b => Some(self.clone()),
            _ => None,
        }
    }

    fn join(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Natural, Self::Natural) => Self::Natural,
            (Self::Natural, _) | (_, Self::Natural) => Self::Undetermined,
            (Self::Lit(a), Self::Lit(b)) => {
                if a == b {
                    Self::Lit(*a)
                } else {
                    Self::Set(Box::new(vec![*a, *b]))
                }
            }
            (Self::Lit(a), Self::Set(s)) | (Self::Set(s), Self::Lit(a)) => {
                let mut v: Vec<i64> = (**s).clone();
                if !v.contains(a) {
                    v.push(*a);
                }
                canonicalize_isotope_set(v).unwrap_or(Self::Undetermined)
            }
            (Self::Set(s), Self::Set(t)) => {
                let mut v: Vec<i64> = (**s).clone();
                for &x in t.iter() {
                    if !v.contains(&x) {
                        v.push(x);
                    }
                }
                canonicalize_isotope_set(v).unwrap_or(Self::Undetermined)
            }
            (Self::Lit(a), Self::Not(b)) | (Self::Not(b), Self::Lit(a)) => {
                if a == b {
                    Self::Undetermined
                } else {
                    Self::Not(*b)
                }
            }
            (Self::Lit(a), Self::NotSet(s)) | (Self::NotSet(s), Self::Lit(a)) => {
                let remaining: Vec<i64> = s.iter().filter(|x| x != &a).copied().collect();
                canonicalize_isotope_not_set(remaining)
            }
            (Self::Set(s), Self::Not(b)) | (Self::Not(b), Self::Set(s)) => {
                if s.contains(b) {
                    Self::Undetermined
                } else {
                    Self::Not(*b)
                }
            }
            (Self::Set(s), Self::NotSet(t)) | (Self::NotSet(t), Self::Set(s)) => {
                let remaining: Vec<i64> = t.iter().filter(|x| !s.contains(x)).copied().collect();
                canonicalize_isotope_not_set(remaining)
            }
            (Self::Not(a), Self::Not(b)) => {
                if a == b {
                    Self::Not(*a)
                } else {
                    Self::Undetermined
                }
            }
            (Self::Not(a), Self::NotSet(s)) | (Self::NotSet(s), Self::Not(a)) => {
                if s.contains(a) {
                    Self::Not(*a)
                } else {
                    Self::Undetermined
                }
            }
            (Self::NotSet(s), Self::NotSet(t)) => {
                let intersection: Vec<i64> = s.iter().filter(|x| t.contains(x)).copied().collect();
                canonicalize_isotope_not_set(intersection)
            }
            (Self::Ref(a), Self::Ref(b)) if a == b => Self::Ref(a.clone()),
            (Self::Bind(a), Self::Bind(b)) if a == b => self.clone(),
            _ => Self::Undetermined,
        }
    }

    fn matches(&self, target: &Self) -> bool {
        match (self, target) {
            (Self::Undetermined, _) => true,
            (_, Self::Undetermined) => false,
            (Self::Natural, Self::Natural) => true,
            (Self::Natural, _) | (_, Self::Natural) => false,
            (Self::Ref(_), _) | (_, Self::Ref(_)) => false,
            _ => {
                let (Some((ps, pp)), Some((ts, tp))) =
                    (isotope_set_view(self), isotope_set_view(target))
                else {
                    return false;
                };
                match (pp, tp) {
                    (MemOp::In, MemOp::In) => ts.iter().all(|t| ps.contains(t)),
                    (MemOp::In, MemOp::NotIn) => false,
                    (MemOp::NotIn, MemOp::In) => ts.iter().all(|t| !ps.contains(t)),
                    (MemOp::NotIn, MemOp::NotIn) => ps.iter().all(|p| ts.contains(p)),
                }
            }
        }
    }
}

/// Canonicalize an isotope-set intersection result.
fn canonicalize_isotope_set(v: Vec<i64>) -> Option<IsotopeMassAst> {
    match v.len() {
        0 => None,
        1 => Some(IsotopeMassAst::Lit(v[0])),
        _ => Some(IsotopeMassAst::Set(Box::new(v))),
    }
}

/// Canonicalize an isotope not-set result.
fn canonicalize_isotope_not_set(v: Vec<i64>) -> IsotopeMassAst {
    match v.len() {
        0 => IsotopeMassAst::Undetermined,
        1 => IsotopeMassAst::Not(v[0]),
        _ => IsotopeMassAst::NotSet(Box::new(v)),
    }
}

/// View any concrete `IsotopeAst` (Lit / Set / Not / NotSet / Bind) as a
/// `(set, polarity)` pair describing the admissible domain. Returns None
/// for `Undetermined`, `Natural`, and `Ref` which don't have a finite set
/// encoding.
fn isotope_set_view(ast: &IsotopeMassAst) -> Option<(&[i64], MemOp)> {
    match ast {
        IsotopeMassAst::Lit(x) => Some((slice::from_ref(x), MemOp::In)),
        IsotopeMassAst::Set(s) => Some((s, MemOp::In)),
        IsotopeMassAst::Not(x) => Some((slice::from_ref(x), MemOp::NotIn)),
        IsotopeMassAst::NotSet(s) => Some((s, MemOp::NotIn)),
        IsotopeMassAst::Bind(b) => Some((&b.2, b.1)),
        IsotopeMassAst::Undetermined | IsotopeMassAst::Natural | IsotopeMassAst::Ref(_) => None,
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::{AtomConstraint, AtomConstraintKind};
    use crate::ast::value::ValueExpr;
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
    #[case::with_isotope_mass(AtomAst::default().with_isotope_mass(12_i64), AtomAst { isotope_mass: IsotopeMassAst::Lit(12), ..Default::default() })]
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
    #[case::isotope_mismatch(AtomAst::from_element(Element::C).with_isotope_mass(12_i64), AtomAst::from_element(Element::C).with_isotope_mass(13_i64), false)]
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
            charge: ValueAst::Expr(Box::new(ValueExpr::Lit(1))),
            implicit_hydrogens: ValueAst::Expr(Box::new(ValueExpr::Lit(3))),
            lone_pairs: ValueAst::Expr(Box::new(ValueExpr::Neg(Box::new(ValueExpr::Lit(2))))),
            spin: SpinStateAst {
                unpaired: ValueAst::Expr(Box::new(ValueExpr::Lit(0))),
                multiplicity: ValueAst::Expr(Box::new(ValueExpr::Lit(1))),
            },
            constraints: AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Expr(
                Box::new(ValueExpr::Lit(4)),
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
    #[case::saturate_contradiction(
        AtomAst::from_element(Element::C).with_constraint(AtomConstraint::JointDomain(
            JointDomainAst::from_ints(
                vec![JointVar::Charge, JointVar::ImplicitHydrogens],
                vec![vec![0, 1], vec![1, 2]],
            )
            .unwrap(),
        )),
        AtomAst::from_element(Element::C).with_charge(5_i64),
        None)]
    fn test_atom_ast_meet(
        #[case] a: AtomAst,
        #[case] b: AtomAst,
        #[case] expected: Option<AtomAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::element_mismatch_widens(AtomAst::from_element(Element::C), AtomAst::from_element(Element::N),
        ElementAst::Set(vec![Element::C, Element::N]))]
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

    #[rstest]
    #[case::prunes_to_pruned_jd(
        AtomAst::from_element(Element::C)
            .with_charge(ValueAst::Set(Box::new(vec![0, 1])))
            .with_constraint(AtomConstraint::JointDomain(
                JointDomainAst::from_ints(
                    vec![JointVar::Charge, JointVar::ImplicitHydrogens],
                    vec![vec![0, 1], vec![1, 2], vec![2, 3]],
                )
                .unwrap(),
            )),
        AtomAst::from_element(Element::C)
            .with_charge(ValueAst::Set(Box::new(vec![0, 1])))
            .with_constraint(AtomConstraint::JointDomain(
                JointDomainAst::from_ints(
                    vec![JointVar::Charge, JointVar::ImplicitHydrogens],
                    vec![vec![0, 1], vec![1, 2]],
                )
                .unwrap(),
            )),
    )]
    #[case::collapses_to_lits(
        AtomAst::from_element(Element::C)
            .with_charge(0_i64)
            .with_constraint(AtomConstraint::JointDomain(
                JointDomainAst::from_ints(
                    vec![JointVar::Charge, JointVar::ImplicitHydrogens],
                    vec![vec![0, 1], vec![1, 2]],
                )
                .unwrap(),
            )),
        AtomAst::from_element(Element::C)
            .with_charge(0_i64)
            .with_implicit_hydrogens(1_i64),
    )]
    #[case::cascades_across_domains(
        AtomAst::from_element(Element::C)
            .with_charge(0_i64)
            .with_constraint(AtomConstraint::JointDomain(
                JointDomainAst::from_ints(
                    vec![JointVar::Charge, JointVar::ImplicitHydrogens],
                    vec![vec![0, 1], vec![1, 2]],
                )
                .unwrap(),
            ))
            .with_constraint(AtomConstraint::JointDomain(
                JointDomainAst::from_ints(
                    vec![JointVar::ImplicitHydrogens, JointVar::LonePairs],
                    vec![vec![1, 3], vec![2, 4]],
                )
                .unwrap(),
            )),
        AtomAst::from_element(Element::C)
            .with_charge(0_i64)
            .with_implicit_hydrogens(1_i64)
            .with_lone_pairs(3_i64),
    )]
    fn test_saturate_atom(#[case] mut atom: AtomAst, #[case] expected: AtomAst) {
        saturate_atom(&mut atom).unwrap();
        assert_eq!(atom, expected);
    }

    #[rstest]
    #[case::no_admissible_tuple(
        AtomAst::from_element(Element::C)
            .with_charge(5_i64)
            .with_constraint(AtomConstraint::JointDomain(
                JointDomainAst::from_ints(
                    vec![JointVar::Charge, JointVar::ImplicitHydrogens],
                    vec![vec![0, 1], vec![1, 2]],
                )
                .unwrap(),
            )),
    )]
    fn test_saturate_atom_error(#[case] mut atom: AtomAst) {
        assert_eq!(saturate_atom(&mut atom).unwrap_err(), Contradiction);
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
    #[case::bind(ElementAst::bind("e", vec![Element::C], MemOp::In), None)]
    #[case::reference(ElementAst::Ref("e".into()), None)]
    fn test_element_ast_literal_and_is_ground(
        #[case] ast: ElementAst,
        #[case] expected: Option<Element>,
    ) {
        assert_eq!(ast.as_lit(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rstest]
    #[case::lit(ElementAst::Lit(Element::C), false)]
    #[case::wildcard(ElementAst::Undetermined, true)]
    #[case::set(ElementAst::Set(vec![Element::C, Element::N]), false)]
    #[case::bind(ElementAst::bind("e", vec![Element::C], MemOp::In), false)]
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
    #[case::bind_lit_match(ElementAst::bind("e", vec![Element::C], MemOp::In), ElementAst::Lit(Element::C), true)]
    #[case::bind_lit_mismatch(ElementAst::bind("e", vec![Element::C], MemOp::In), ElementAst::Lit(Element::N), false)]
    #[case::bind_set_subset(ElementAst::bind("e", vec![Element::C, Element::N], MemOp::In), ElementAst::Set(vec![Element::C]), true)]
    #[case::set_bind_subset(ElementAst::Set(vec![Element::C, Element::N]), ElementAst::bind("e", vec![Element::C], MemOp::In), true)]
    #[case::bind_bind_subset(ElementAst::bind("p", vec![Element::C, Element::N], MemOp::In), ElementAst::bind("t", vec![Element::N], MemOp::In), true)]
    #[case::bind_bind_superset(ElementAst::bind("p", vec![Element::C], MemOp::In), ElementAst::bind("t", vec![Element::C, Element::N], MemOp::In), false)]
    #[case::undetermined_bind(ElementAst::Undetermined, ElementAst::bind("e", vec![Element::C], MemOp::In), true)]
    #[case::bind_undetermined(ElementAst::bind("e", vec![Element::C], MemOp::In), ElementAst::Undetermined, false)]
    #[case::ref_lit(ElementAst::Ref("e".into()), ElementAst::Lit(Element::C), false)]
    #[case::lit_ref(ElementAst::Lit(Element::C), ElementAst::Ref("e".into()), false)]
    #[case::ref_set(ElementAst::Ref("e".into()), ElementAst::Set(vec![Element::C]), false)]
    #[case::set_ref(ElementAst::Set(vec![Element::C]), ElementAst::Ref("e".into()), false)]
    #[case::ref_bind(ElementAst::Ref("e".into()), ElementAst::bind("f", vec![Element::C], MemOp::In), false)]
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
    #[case::positive(IsotopeMassAst::from(13_i64), IsotopeMassAst::Lit(13))]
    #[case::zero(IsotopeMassAst::from(0_i64), IsotopeMassAst::Lit(0))]
    fn test_isotope_ast_from_i64(#[case] actual: IsotopeMassAst, #[case] expected: IsotopeMassAst) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::natural(IsotopeMassAst::Natural, false)]
    #[case::lit(IsotopeMassAst::Lit(12), false)]
    #[case::undetermined(IsotopeMassAst::Undetermined, true)]
    #[case::set(IsotopeMassAst::Set(Box::new(vec![12, 13])), false)]
    #[case::not(IsotopeMassAst::Not(14), false)]
    #[case::not_set(IsotopeMassAst::NotSet(Box::new(vec![12, 13])), false)]
    #[case::bind(IsotopeMassAst::bind("m", vec![12], MemOp::In), false)]
    #[case::reference(IsotopeMassAst::Ref("m".into()), false)]
    fn test_isotope_ast_is_undetermined(#[case] ast: IsotopeMassAst, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rstest]
    #[case::natural(IsotopeMassAst::Natural, Some(0))]
    #[case::lit(IsotopeMassAst::Lit(12), Some(12))]
    #[case::lit_zero(IsotopeMassAst::Lit(0), Some(0))]
    #[case::wildcard(IsotopeMassAst::Undetermined, None)]
    #[case::set_singleton(IsotopeMassAst::Set(Box::new(vec![14])), Some(14))]
    #[case::set_multi(IsotopeMassAst::Set(Box::new(vec![12, 13])), None)]
    #[case::not(IsotopeMassAst::Not(14), None)]
    #[case::not_set(IsotopeMassAst::NotSet(Box::new(vec![12, 13])), None)]
    #[case::bind(IsotopeMassAst::bind("m", vec![12], MemOp::In), None)]
    #[case::reference(IsotopeMassAst::Ref("m".into()), None)]
    fn test_isotope_ast_literal_and_is_ground(
        #[case] ast: IsotopeMassAst,
        #[case] expected: Option<u32>,
    ) {
        assert_eq!(ast.as_lit(), expected);
        assert_eq!(ast.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_natural(IsotopeMassAst::Undetermined, IsotopeMassAst::Natural, true)]
    #[case::undetermined_value(IsotopeMassAst::Undetermined, IsotopeMassAst::Lit(12), true)]
    #[case::undetermined_undetermined(IsotopeMassAst::Undetermined, IsotopeMassAst::Undetermined, true)]
    #[case::natural_undetermined(IsotopeMassAst::Natural, IsotopeMassAst::Undetermined, false)]
    #[case::value_undetermined(IsotopeMassAst::Lit(12), IsotopeMassAst::Undetermined, false)]
    #[case::natural_natural(IsotopeMassAst::Natural, IsotopeMassAst::Natural, true)]
    #[case::natural_value(IsotopeMassAst::Natural, IsotopeMassAst::Lit(12), false)]
    #[case::value_natural(IsotopeMassAst::Lit(12), IsotopeMassAst::Natural, false)]
    #[case::value_lit_match(IsotopeMassAst::Lit(12), IsotopeMassAst::Lit(12), true)]
    #[case::value_lit_mismatch(IsotopeMassAst::Lit(12), IsotopeMassAst::Lit(13), false)]
    #[case::value_wildcard_lit(IsotopeMassAst::Undetermined, IsotopeMassAst::Lit(12), true)]
    #[case::value_set_lit_in(IsotopeMassAst::Set(Box::new(vec![12, 13])), IsotopeMassAst::Lit(13), true)]
    #[case::value_set_lit_out(IsotopeMassAst::Set(Box::new(vec![12, 13])), IsotopeMassAst::Lit(14), false)]
    #[case::value_set_set_subset(IsotopeMassAst::Set(Box::new(vec![12, 13, 14])), IsotopeMassAst::Set(Box::new(vec![12, 13])), true)]
    #[case::value_set_set_superset(IsotopeMassAst::Set(Box::new(vec![12])), IsotopeMassAst::Set(Box::new(vec![12, 13])), false)]
    #[case::value_lit_wildcard(IsotopeMassAst::Lit(12), IsotopeMassAst::Undetermined, false)]
    fn test_isotope_ast_matches(
        #[case] pattern: IsotopeMassAst,
        #[case] target: IsotopeMassAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
