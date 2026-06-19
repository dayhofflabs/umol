//! Stereochemistry AST: the configuration value and the operator-expression
//! tree over it.
//!
//! A configuration value is a dense coset index per stereo kind, corresponds to OpenSMILES
//! numbering for SP, TB, and OH.
//! `~` and `^` are group actions on the index; the owning configuration's
//! `canonicalize` folds closed operator-expressions against the coset algebra.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::mem;

use strum::VariantArray;
use umol_ast_macros::Canonicalize;
use umol_perm::{space, ClassKey, Permutation};

use super::constraint::{
    StereoAtomConstraint, StereoAtomConstraints, StereoBondConstraint, StereoBondConstraints,
};
use super::error::Contradiction;
use super::traits::{AsLit, Canonicalize, Lattice};

/// Stereo kind: the atom-centered coordination geometries and the bond-centered cis/trans kind.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, strum::EnumCount)]
pub enum StereoKind {
    Tetrahedral,
    CisTrans,
    Axial,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
}

impl StereoKind {
    /// The `umol-perm` class key for this stereo kind.
    pub fn class_key(self) -> ClassKey {
        match self {
            StereoKind::Tetrahedral => ClassKey::Tetrahedral,
            StereoKind::CisTrans => ClassKey::CisTrans,
            StereoKind::Axial => ClassKey::Axial,
            StereoKind::SquarePlanar => ClassKey::SquarePlanar,
            StereoKind::TrigonalBipyramidal => ClassKey::TrigonalBipyramidal,
            StereoKind::Octahedral => ClassKey::Octahedral,
        }
    }

    /// Number of ligand positions in this stereo kind.
    pub fn degree(self) -> usize {
        space(self.class_key()).degree()
    }

    /// Number of cosets/configurations in this stereo kind.
    pub fn count(self) -> usize {
        space(self.class_key()).count()
    }

    /// Whether this stereo kind can encode local handedness.
    pub fn is_chiral_class(self) -> bool {
        space(self.class_key()).is_chiral()
    }

    /// Kind-specific `~` involution. Chiral kinds borrow the orientation-reversing
    /// generator from umol-perm; achiral kinds use a chosen ligand swap (no improper
    /// generator to borrow — theirs is the identity):
    /// - cis/trans: swap the two configurations
    /// - square-planar: swap the diagonal ligand pair
    pub fn involution(self) -> Permutation {
        let coset_space = space(self.class_key());
        if coset_space.is_chiral() {
            coset_space.improper()
        } else {
            match self {
                StereoKind::CisTrans => Permutation::from_image(4, &[1, 0, 2, 3]),
                StereoKind::SquarePlanar => Permutation::from_image(4, &[2, 1, 0, 3]),
                _ => unreachable!("only achiral kinds reach the chosen-swap branch"),
            }
        }
    }

    /// Act on coset index `index` by `permutation`, through the class's coset algebra.
    #[allow(unused)]
    fn act(self, index: u32, permutation: Permutation) -> u32 {
        space(self.class_key()).reindex(index, permutation)
    }
}

/// Topicity of two ligand positions of a stereo carrier (derived ground value).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, VariantArray)]
pub enum Topicity {
    Homotopic,
    Enantiotopic,
    Diastereotopic,
}

/// Stereogenicity classification of a stereo carrier (derived ground value).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, VariantArray)]
pub enum Stereogenicity {
    Symmetric,
    Prochiral,
    Stereogenic,
}

/// Element-side stereo configuration: either undetermined (geometry not yet
/// known, so no coset) or `Kinded` — a concrete geometry bound to a coset that
/// may still be open. `*` (`Undetermined`) and `Th*` (`Kinded(Tetrahedral,
/// Undetermined)`) are distinct. `canonicalize` folds the coset under the kind;
/// no physical range-check (tier-2; the validator does it).
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoConfigurationAst {
    #[default]
    Undetermined,
    Kinded(StereoKind, StereoCosetAst),
}

impl StereoConfigurationAst {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn kinded(kind: StereoKind, coset: impl Into<StereoCosetAst>) -> Self {
        Self::Kinded(kind, coset.into())
    }

    /// The coordination-geometry kind, or `None` when undetermined.
    pub fn kind(&self) -> Option<StereoKind> {
        match self {
            Self::Kinded(kind, _) => Some(*kind),
            Self::Undetermined => None,
        }
    }

    /// The coset, or `None` when undetermined.
    pub fn coset(&self) -> Option<&StereoCosetAst> {
        match self {
            Self::Kinded(_, coset) => Some(coset),
            Self::Undetermined => None,
        }
    }

    /// Mutable access to the coset, or `None` when undetermined.
    pub fn coset_mut(&mut self) -> Option<&mut StereoCosetAst> {
        match self {
            Self::Kinded(_, coset) => Some(coset),
            Self::Undetermined => None,
        }
    }
}

impl From<(StereoKind, u32)> for StereoConfigurationAst {
    fn from((kind, coset): (StereoKind, u32)) -> Self {
        Self::Kinded(kind, StereoCosetAst::Lit(coset))
    }
}

impl Canonicalize for StereoConfigurationAst {
    fn canonicalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Kinded(kind, coset) => Self::Kinded(kind, canon_coset(coset, kind)?),
            Self::Undetermined => Self::Undetermined,
        })
    }

    fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
        match self {
            Self::Kinded(..) => Ok(Cow::Owned(self.clone().canonicalize()?)),
            Self::Undetermined => Ok(Cow::Borrowed(self)),
        }
    }
}

impl AsLit for StereoConfigurationAst {
    type Lit = StereoConfiguration;

    fn as_lit(&self) -> Option<StereoConfiguration> {
        match self {
            Self::Kinded(kind, coset) => coset
                .as_lit()
                .map(|coset| StereoConfiguration { kind: *kind, coset }),
            Self::Undetermined => None,
        }
    }
}

impl Lattice for StereoConfigurationAst {
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        matches!(self, Self::Kinded(_, StereoCosetAst::Lit(_)))
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        let a = self.canonical().ok()?;
        let b = other.canonical().ok()?;
        match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
            (Self::Kinded(k1, ca), Self::Kinded(k2, cb)) => {
                if k1 != k2 {
                    return None;
                }
                Some(Self::Kinded(*k1, coset_meet(ca, cb, *k1)?))
            }
        }
    }

    fn join(&self, other: &Self) -> Self {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Kinded(k1, ca), Self::Kinded(k2, cb)) => {
                if k1 != k2 {
                    Self::Undetermined
                } else {
                    Self::Kinded(*k1, coset_join(ca, cb, *k1))
                }
            }
        }
    }

    fn matches(&self, target: &Self) -> bool {
        match (self.meet(target), target.canonical()) {
            (Some(meet), Ok(target)) => meet == *target,
            _ => false,
        }
    }
}

/// Generates a constraint-side stereo state for a fixed geometry (`#T`/`#C`):
/// undetermined, explicitly not-stereo, or a stereo center with a coset. The
/// geometry is the type's identity (`$kind`), so the coset folds/meets under that
/// constant kind — no kind field.
macro_rules! stereo_site {
    ($name:ident, $kind:expr) => {
        #[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum $name {
            #[default]
            Undetermined,
            NotStereo,
            Stereo(StereoCosetAst),
        }

        impl $name {
            pub fn undetermined() -> Self {
                Self::Undetermined
            }

            pub fn not_stereo() -> Self {
                Self::NotStereo
            }

            pub fn stereo(coset: impl Into<StereoCosetAst>) -> Self {
                Self::Stereo(coset.into())
            }

            pub fn is_stereo(&self) -> bool {
                matches!(self, Self::Stereo(_))
            }

            /// Matches literal coset index `value` under the type's kind.
            pub fn matches_value(&self, value: u32) -> bool {
                match self {
                    Self::Stereo(coset) => coset_matches(coset, &StereoCosetAst::Lit(value), $kind),
                    Self::NotStereo => false,
                    Self::Undetermined => true,
                }
            }
        }

        impl Canonicalize for $name {
            fn canonicalize(self) -> Result<Self, Contradiction> {
                Ok(match self {
                    Self::Stereo(coset) => Self::Stereo(canon_coset(coset, $kind)?),
                    other => other,
                })
            }

            fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction> {
                match self {
                    Self::Stereo(_) => Ok(Cow::Owned(self.clone().canonicalize()?)),
                    _ => Ok(Cow::Borrowed(self)),
                }
            }
        }

        impl AsLit for $name {
            type Lit = StereoConfiguration;

            /// The ground configuration (kind = the type's constant) when the coset
            /// is a literal; `NotStereo`/`Undetermined`/non-literal → `None`.
            fn as_lit(&self) -> Option<StereoConfiguration> {
                match self {
                    Self::Stereo(StereoCosetAst::Lit(coset)) => Some(StereoConfiguration {
                        kind: $kind,
                        coset: *coset,
                    }),
                    _ => None,
                }
            }
        }

        impl Lattice for $name {
            fn is_undetermined(&self) -> bool {
                matches!(self, Self::Undetermined)
            }

            fn is_ground(&self) -> bool {
                match self {
                    Self::NotStereo => true,
                    Self::Stereo(coset) => matches!(coset, StereoCosetAst::Lit(_)),
                    Self::Undetermined => false,
                }
            }

            fn meet(&self, other: &Self) -> Option<Self> {
                let a = self.canonical().ok()?;
                let b = other.canonical().ok()?;
                match (a.as_ref(), b.as_ref()) {
                    (Self::Undetermined, x) | (x, Self::Undetermined) => Some(x.clone()),
                    (Self::NotStereo, Self::NotStereo) => Some(Self::NotStereo),
                    (Self::NotStereo, Self::Stereo(_)) | (Self::Stereo(_), Self::NotStereo) => None,
                    (Self::Stereo(ca), Self::Stereo(cb)) => {
                        Some(Self::Stereo(coset_meet(ca, cb, $kind)?))
                    }
                }
            }

            fn join(&self, other: &Self) -> Self {
                let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
                let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
                match (a.as_ref(), b.as_ref()) {
                    (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
                    (Self::NotStereo, Self::NotStereo) => Self::NotStereo,
                    (Self::NotStereo, Self::Stereo(_)) | (Self::Stereo(_), Self::NotStereo) => {
                        Self::Undetermined
                    }
                    (Self::Stereo(ca), Self::Stereo(cb)) => Self::Stereo(coset_join(ca, cb, $kind)),
                }
            }

            fn matches(&self, target: &Self) -> bool {
                match (self.meet(target), target.canonical()) {
                    (Some(meet), Ok(target)) => meet == *target,
                    _ => false,
                }
            }
        }
    };
}

stereo_site! { TetrahedralStereoAst, StereoKind::Tetrahedral }
stereo_site! { CisTransStereoAst, StereoKind::CisTrans }

/// Operator-expression term: a `Var` (with optional finite domain), a literal
/// `Lit`/`LitSet` base, or one of these under the permutation-action operators
/// `~` (swap), `'` (mirror), `^` (apply). Kind-relative — **no
/// `Lattice`/`Canonicalize`** (structural `Eq` only); the owning configuration
/// normalizes it under its concrete kind. Canonicalization composes the operator
/// word into one net permutation: over a literal base it folds to a concrete
/// coset; over a `Var` it leaves at most one operator layer.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoTerm {
    Var(Box<(String, Option<BTreeSet<u32>>)>),
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Swap(Box<StereoTerm>),
    Mirror(Box<StereoTerm>),
    Apply(Box<StereoTerm>, Permutation),
}

impl StereoTerm {
    /// A free coset variable `?name`.
    pub fn var(name: impl Into<String>) -> Self {
        Self::Var(Box::new((name.into(), None)))
    }

    /// A variable restricted to a finite coset domain `?name :: {…}`.
    pub fn var_in(name: impl Into<String>, domain: impl IntoIterator<Item = u32>) -> Self {
        Self::Var(Box::new((name.into(), Some(domain.into_iter().collect()))))
    }

    /// A finite literal-set base `{…}` (folds under the owner's kind).
    pub fn lit_set(values: impl IntoIterator<Item = u32>) -> Self {
        Self::LitSet(values.into_iter().collect())
    }

    /// `~inner` — the kind involution applied to `inner`.
    pub fn swap(inner: Self) -> Self {
        Self::Swap(Box::new(inner))
    }

    /// `'inner` — the enantiomer (mirror) of `inner`.
    pub fn mirror(inner: Self) -> Self {
        Self::Mirror(Box::new(inner))
    }

    /// `inner ^ permutation` — the group action of `permutation` on `inner`.
    pub fn apply(inner: Self, permutation: Permutation) -> Self {
        Self::Apply(Box::new(inner), permutation)
    }
}

/// Dense coset-index expression (0-indexed per stereo kind): undetermined, a
/// single index, a finite set, a complement set, or an operator `Term` over a
/// variable. Kind-relative — **no `Lattice`/`Canonicalize`/`AsLit`**; the owning
/// the owning configuration or site normalizes it under its kind.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoCosetAst {
    #[default]
    Undetermined,
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Term(Box<StereoTerm>),
}

impl StereoCosetAst {
    pub fn lit_set(values: impl IntoIterator<Item = u32>) -> Self {
        Self::LitSet(values.into_iter().collect())
    }

    pub fn term(term: StereoTerm) -> Self {
        Self::Term(Box::new(term))
    }
}

impl From<u32> for StereoCosetAst {
    fn from(index: u32) -> Self {
        Self::Lit(index)
    }
}

impl AsLit for StereoCosetAst {
    type Lit = u32;

    /// The single coset index, only when literal. Kind-independent — so it lives on
    /// the bare coset, unlike the kind-aware folding ops.
    #[inline]
    fn as_lit(&self) -> Option<u32> {
        match self {
            Self::Lit(i) => Some(*i),
            _ => None,
        }
    }
}

/// A ground stereo configuration: a concrete geometry plus its coset index. The
/// `AsLit` target of `StereoConfigurationAst` and the per-kind site types.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct StereoConfiguration {
    pub kind: StereoKind,
    pub coset: u32,
}

impl StereoKind {
    /// The mirror (improper, μ) generator as a permutation: chiral kinds use the
    /// orientation-reversing generator; achiral kinds act trivially on cosets.
    pub(crate) fn mirror_permutation(self) -> Permutation {
        if self.is_chiral_class() {
            space(self.class_key()).improper()
        } else {
            Permutation::identity(self.degree())
        }
    }

    /// Whether `g` and `h` induce the same coset relabeling for this kind.
    fn coset_action_eq(self, g: Permutation, h: Permutation) -> bool {
        let s = space(self.class_key());
        (0..s.count() as u32).all(|i| s.reindex(i, g) == s.reindex(i, h))
    }
}

/// The literal coset-index set a positive coset denotes; `None` for the wildcard
/// `Undetermined` and the symbolic `Term`. Used by `coset_meet`/`coset_join`
/// after those two cases are handled.
fn coset_to_set(coset: &StereoCosetAst) -> Option<BTreeSet<u32>> {
    match coset {
        StereoCosetAst::Lit(i) => Some(BTreeSet::from([*i])),
        StereoCosetAst::LitSet(s) => Some(s.clone()),
        StereoCosetAst::Undetermined | StereoCosetAst::Term(_) => None,
    }
}

/// Walk a term's operator word into one net coset permutation (composed inner →
/// outer), returning the base leaf (`Var`/`Lit`/`LitSet`) and that permutation.
fn compose_term(term: &StereoTerm, kind: StereoKind) -> (&StereoTerm, Permutation) {
    match term {
        StereoTerm::Swap(inner) => {
            let (base, g) = compose_term(inner, kind);
            (base, g.compose(kind.involution()))
        }
        StereoTerm::Mirror(inner) => {
            let (base, g) = compose_term(inner, kind);
            (base, g.compose(kind.mirror_permutation()))
        }
        StereoTerm::Apply(inner, p) => {
            let (base, g) = compose_term(inner, kind);
            (base, g.compose(*p))
        }
        base => (base, Permutation::identity(kind.degree())),
    }
}

/// Canonicalize a coset under `kind`. A `Term` over a `Var` renders by priority
/// Mirror > Swap > Apply (canonicalizing the domain); every other form reduces to
/// a literal index set that folds: ∅ → `Err` (the bottom `meet` uses to signal
/// incompatible cosets), singleton → `Lit`, else `LitSet`. No universe folding
/// (`full → Undetermined`) and no range-check — both are tier-2 (the validator).
pub(crate) fn canon_coset(
    coset: StereoCosetAst,
    kind: StereoKind,
) -> Result<StereoCosetAst, Contradiction> {
    let s = space(kind.class_key());
    let set: BTreeSet<u32> = match &coset {
        StereoCosetAst::Undetermined => return Ok(StereoCosetAst::Undetermined),
        StereoCosetAst::Lit(i) => BTreeSet::from([*i]),
        StereoCosetAst::LitSet(values) => values.clone(),
        StereoCosetAst::Term(t) => {
            let (base, g) = compose_term(t, kind);
            match base {
                StereoTerm::Var(v) => {
                    let n = kind.count() as u32;
                    let domain = match &v.1 {
                        Some(set) if set.is_empty() => return Err(Contradiction),
                        Some(set) if set.len() as u32 == n => None,
                        Some(set) => Some(set.clone()),
                        None => None,
                    };
                    let var = StereoTerm::Var(Box::new((v.0.clone(), domain)));
                    let id = Permutation::identity(kind.degree());
                    let term = if kind.coset_action_eq(g, id) {
                        var
                    } else if kind.is_chiral_class()
                        && kind.coset_action_eq(g, kind.mirror_permutation())
                    {
                        StereoTerm::Mirror(Box::new(var))
                    } else if kind.coset_action_eq(g, kind.involution()) {
                        StereoTerm::Swap(Box::new(var))
                    } else {
                        StereoTerm::Apply(Box::new(var), g)
                    };
                    return Ok(StereoCosetAst::term(term));
                }
                StereoTerm::Lit(i) => BTreeSet::from([s.reindex(*i, g)]),
                StereoTerm::LitSet(values) => values.iter().map(|i| s.reindex(*i, g)).collect(),
                StereoTerm::Swap(_) | StereoTerm::Mirror(_) | StereoTerm::Apply(..) => {
                    unreachable!("compose_term returns a base leaf")
                }
            }
        }
    };
    if set.is_empty() {
        Err(Contradiction)
    } else if set.len() == 1 {
        Ok(StereoCosetAst::Lit(set.into_iter().next().unwrap()))
    } else {
        Ok(StereoCosetAst::LitSet(set))
    }
}

/// Greatest lower bound of two cosets under `kind` (canonicalizing operands);
/// `Term` meets only an equal canonical `Term`.
pub(crate) fn coset_meet(
    a: &StereoCosetAst,
    b: &StereoCosetAst,
    kind: StereoKind,
) -> Option<StereoCosetAst> {
    let ca = canon_coset(a.clone(), kind).ok()?;
    let cb = canon_coset(b.clone(), kind).ok()?;
    use StereoCosetAst::{Term, Undetermined};
    match (&ca, &cb) {
        (Undetermined, _) => Some(cb),
        (_, Undetermined) => Some(ca),
        (Term(_), Term(_)) => (ca == cb).then_some(ca),
        (Term(_), _) | (_, Term(_)) => None,
        _ => {
            let sa = coset_to_set(&ca).unwrap();
            let sb = coset_to_set(&cb).unwrap();
            canon_coset(
                StereoCosetAst::LitSet(sa.intersection(&sb).copied().collect()),
                kind,
            )
            .ok()
        }
    }
}

/// Least upper bound of two cosets under `kind`.
pub(crate) fn coset_join(
    a: &StereoCosetAst,
    b: &StereoCosetAst,
    kind: StereoKind,
) -> StereoCosetAst {
    let ca = canon_coset(a.clone(), kind).unwrap_or(StereoCosetAst::Undetermined);
    let cb = canon_coset(b.clone(), kind).unwrap_or(StereoCosetAst::Undetermined);
    use StereoCosetAst::{Term, Undetermined};
    match (&ca, &cb) {
        (Undetermined, _) | (_, Undetermined) => StereoCosetAst::Undetermined,
        (Term(_), Term(_)) if ca == cb => ca,
        (Term(_), _) | (_, Term(_)) => StereoCosetAst::Undetermined,
        _ => {
            let sa = coset_to_set(&ca).unwrap();
            let sb = coset_to_set(&cb).unwrap();
            canon_coset(
                StereoCosetAst::LitSet(sa.union(&sb).copied().collect()),
                kind,
            )
            .unwrap_or(StereoCosetAst::Undetermined)
        }
    }
}

/// `target` refines `pattern` under `kind` (meet-derived).
pub(crate) fn coset_matches(
    pattern: &StereoCosetAst,
    target: &StereoCosetAst,
    kind: StereoKind,
) -> bool {
    match (
        coset_meet(pattern, target, kind),
        canon_coset(target.clone(), kind),
    ) {
        (Some(m), Ok(ct)) => m == ct,
        _ => false,
    }
}

/// Apply a ligand-order permutation to a coset under `kind`.
pub(crate) fn coset_apply_permutation(
    coset: &StereoCosetAst,
    permutation: Permutation,
    kind: StereoKind,
) -> StereoCosetAst {
    let s = space(kind.class_key());
    match coset {
        StereoCosetAst::Undetermined => StereoCosetAst::Undetermined,
        StereoCosetAst::Lit(i) => StereoCosetAst::Lit(s.reindex(*i, permutation)),
        StereoCosetAst::LitSet(set) => {
            StereoCosetAst::LitSet(set.iter().map(|i| s.reindex(*i, permutation)).collect())
        }
        StereoCosetAst::Term(t) => canon_coset(
            StereoCosetAst::term(StereoTerm::apply((**t).clone(), permutation)),
            kind,
        )
        .unwrap_or(StereoCosetAst::Undetermined),
    }
}

/// StereoAtomAst and StereoBondAst generator
macro_rules! stereo_element {
    (
        $(#[doc = $doc:literal])+
        $name:ident, $constraints:ident, $constraint:ident
    ) => {
        $(#[doc = $doc])+
        #[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize)]
        pub struct $name {
            pub configuration: StereoConfigurationAst,
            pub constraints: $constraints,
        }

        impl $name {
            pub fn new(kind: StereoKind, coset: impl Into<StereoCosetAst>) -> Self {
                Self {
                    configuration: StereoConfigurationAst::kinded(kind, coset),
                    constraints: $constraints::new(),
                }
            }

            /// Add a single constraint.
            pub fn with_constraint(self, _constraint: impl Into<$constraint>) -> Self {
                self
            }

            /// Add each constraint from iterator.
            pub fn with_constraints<I>(mut self, constraints: I) -> Self
            where
                I: IntoIterator,
                I::Item: Into<$constraint>,
            {
                self.constraints.extend(constraints.into_iter().map(Into::into));
                self
            }

            /// No-op. A stereo element is always stereogenic, so its coset has no
            /// zero default; it is ground iff its coset is ground.
            pub fn into_ground(self) -> Self {
                self
            }

            /// Equivalent to `into_ground()`; there are no constraint defaults.
            pub fn into_zeroed(self) -> Self {
                self.into_ground()
            }

            /// Canonicalize the configuration (fold the coset under its kind) and
            /// simplify each constraint in place.
            pub fn simplify_values(&mut self) {
                let cfg = mem::take(&mut self.configuration);
                self.configuration = cfg.clone().canonicalize().unwrap_or(cfg);
                self.constraints.simplify_each();
            }
        }

        impl Lattice for $name {
            fn is_undetermined(&self) -> bool {
                self.configuration.is_undetermined() && self.constraints.is_undetermined()
            }

            fn is_ground(&self) -> bool {
                self.configuration.is_ground() && self.constraints.is_ground()
            }

            fn meet(&self, other: &Self) -> Option<Self> {
                Some(Self {
                    configuration: self.configuration.meet(&other.configuration)?,
                    constraints: self.constraints.meet(&other.constraints)?,
                })
            }

            fn join(&self, other: &Self) -> Self {
                Self {
                    configuration: self.configuration.join(&other.configuration),
                    constraints: self.constraints.join(&other.constraints),
                }
            }

            fn matches(&self, target: &Self) -> bool {
                self.configuration.matches(&target.configuration)
                    && self.constraints.matches(&target.constraints)
            }
        }
    };
}

stereo_element! {
    /// StereoAtomAst with geometry class, configuration, and per-site constraints.
    StereoAtomAst, StereoAtomConstraints, StereoAtomConstraint
}

stereo_element! {
    /// StereoBondAst with cis/trans configuration and per-site constraints.
    StereoBondAst, StereoBondConstraints, StereoBondConstraint
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, ClassKey::Tetrahedral)]
    #[case::cis_trans(StereoKind::CisTrans, ClassKey::CisTrans)]
    #[case::axial(StereoKind::Axial, ClassKey::Axial)]
    #[case::square_planar(StereoKind::SquarePlanar, ClassKey::SquarePlanar)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, ClassKey::TrigonalBipyramidal)]
    #[case::octahedral(StereoKind::Octahedral, ClassKey::Octahedral)]
    fn test_stereo_kind_class_key(#[case] kind: StereoKind, #[case] expected: ClassKey) {
        assert_eq!(kind.class_key(), expected);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, 4)]
    #[case::cis_trans(StereoKind::CisTrans, 4)]
    #[case::axial(StereoKind::Axial, 4)]
    #[case::square_planar(StereoKind::SquarePlanar, 4)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, 5)]
    #[case::octahedral(StereoKind::Octahedral, 6)]
    fn test_stereo_kind_degree(#[case] kind: StereoKind, #[case] expected: usize) {
        assert_eq!(kind.degree(), expected);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, 2)]
    #[case::cis_trans(StereoKind::CisTrans, 2)]
    #[case::axial(StereoKind::Axial, 2)]
    #[case::square_planar(StereoKind::SquarePlanar, 3)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, 20)]
    #[case::octahedral(StereoKind::Octahedral, 30)]
    fn test_stereo_kind_count(#[case] kind: StereoKind, #[case] expected: usize) {
        assert_eq!(kind.count(), expected);
    }

    #[rstest]
    #[case::tetrahedral(StereoKind::Tetrahedral, true)]
    #[case::cis_trans(StereoKind::CisTrans, false)]
    #[case::axial(StereoKind::Axial, true)]
    #[case::square_planar(StereoKind::SquarePlanar, false)]
    #[case::trigonal_bipyramidal(StereoKind::TrigonalBipyramidal, true)]
    #[case::octahedral(StereoKind::Octahedral, true)]
    fn test_stereo_kind_is_chiral_class(#[case] kind: StereoKind, #[case] expected: bool) {
        assert_eq!(kind.is_chiral_class(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::tetrahedral((StereoKind::Tetrahedral, 1), StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)))]
    #[case::octahedral((StereoKind::Octahedral, 5), StereoConfigurationAst::Kinded(StereoKind::Octahedral, StereoCosetAst::Lit(5)))]
    fn test_stereo_configuration_ast_from(#[case] input: (StereoKind, u32), #[case] expected: StereoConfigurationAst) {
        assert_eq!(StereoConfigurationAst::from(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::term_swap_folds_to_lit(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::term(StereoTerm::swap(StereoTerm::Lit(0)))), StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)))]
    fn test_stereo_configuration_ast_canonicalize(#[case] input: StereoConfigurationAst, #[case] expected: StereoConfigurationAst) {
        assert_eq!(input.canonicalize(), Ok(expected));
    }

    #[rstest]
    #[case::undetermined(StereoConfigurationAst::Undetermined)]
    #[case::kind_lit(StereoConfigurationAst::Kinded(
        StereoKind::Tetrahedral,
        StereoCosetAst::Lit(0)
    ))]
    #[case::kind_open(StereoConfigurationAst::Kinded(
        StereoKind::Tetrahedral,
        StereoCosetAst::Undetermined
    ))]
    // Multi-element / full coset sets are preserved (no complement or full→Undetermined fold).
    #[case::multi_element_set(StereoConfigurationAst::Kinded(StereoKind::SquarePlanar, StereoCosetAst::lit_set([0, 1])))]
    #[case::full_set(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::lit_set([0, 1])))]
    fn test_stereo_configuration_ast_canonicalize_identity(#[case] input: StereoConfigurationAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rstest]
    #[case::empty_set(StereoConfigurationAst::Kinded(
        StereoKind::SquarePlanar,
        StereoCosetAst::LitSet(BTreeSet::new())
    ))]
    fn test_stereo_configuration_ast_canonicalize_error(#[case] input: StereoConfigurationAst) {
        assert_eq!(input.canonicalize(), Err(Contradiction));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_lit(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)), Some(StereoConfiguration { kind: StereoKind::Tetrahedral, coset: 1 }))]
    #[case::undetermined(StereoConfigurationAst::Undetermined, None)]
    #[case::kind_open(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), None)]
    fn test_stereo_configuration_ast_as_lit(#[case] config: StereoConfigurationAst, #[case] expected: Option<StereoConfiguration>) {
        assert_eq!(config.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationAst::Undetermined, true)]
    #[case::kind_open(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), false)]
    #[case::kind_lit(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)), false)]
    fn test_stereo_configuration_ast_is_undetermined(#[case] config: StereoConfigurationAst, #[case] expected: bool) {
        assert_eq!(config.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_lit(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)), true)]
    #[case::undetermined(StereoConfigurationAst::Undetermined, false)]
    #[case::kind_open(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), false)]
    fn test_stereo_configuration_ast_is_ground(#[case] config: StereoConfigurationAst, #[case] expected: bool) {
        assert_eq!(config.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_narrows(StereoConfigurationAst::Undetermined, StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), Some(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0))))]
    #[case::coset_same(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), Some(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0))))]
    #[case::open_narrows(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), Some(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0))))]
    #[case::coset_conflict(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::Tetrahedral, 1)), None)]
    #[case::kind_conflict(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::CisTrans, 0)), None)]
    fn test_stereo_configuration_ast_meet(#[case] a: StereoConfigurationAst, #[case] b: StereoConfigurationAst, #[case] expected: Option<StereoConfigurationAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_absorbs(StereoConfigurationAst::Undetermined, StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::Undetermined)]
    #[case::coset_same(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)))]
    #[case::coset_widens(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::Tetrahedral, 1)), StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::lit_set([0, 1])))]
    #[case::kind_conflict(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::CisTrans, 0)), StereoConfigurationAst::Undetermined)]
    fn test_stereo_configuration_ast_join(#[case] a: StereoConfigurationAst, #[case] b: StereoConfigurationAst, #[case] expected: StereoConfigurationAst) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_matches_any(StereoConfigurationAst::Undetermined, StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), true)]
    #[case::open_matches_lit(StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), true)]
    #[case::specific_vs_undetermined(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::Undetermined, false)]
    #[case::coset_match(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), true)]
    #[case::coset_mismatch(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::Tetrahedral, 1)), false)]
    #[case::kind_mismatch(StereoConfigurationAst::from((StereoKind::Tetrahedral, 0)), StereoConfigurationAst::from((StereoKind::CisTrans, 0)), false)]
    fn test_stereo_configuration_ast_matches(#[case] pattern: StereoConfigurationAst, #[case] target: StereoConfigurationAst, #[case] expected: bool) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::term_swap_folds(TetrahedralStereoAst::Stereo(StereoCosetAst::term(StereoTerm::swap(StereoTerm::Lit(0)))), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1)))]
    fn test_tetrahedral_stereo_ast_canonicalize(#[case] input: TetrahedralStereoAst, #[case] expected: TetrahedralStereoAst) {
        assert_eq!(input.canonicalize(), Ok(expected));
    }

    #[rstest]
    #[case::undetermined(TetrahedralStereoAst::Undetermined)]
    #[case::not_stereo(TetrahedralStereoAst::NotStereo)]
    #[case::stereo_lit(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)))]
    #[case::stereo_open(TetrahedralStereoAst::Stereo(StereoCosetAst::Undetermined))]
    fn test_tetrahedral_stereo_ast_canonicalize_identity(#[case] input: TetrahedralStereoAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereo_lit(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1)), Some(StereoConfiguration { kind: StereoKind::Tetrahedral, coset: 1 }))]
    #[case::not_stereo(TetrahedralStereoAst::NotStereo, None)]
    #[case::undetermined(TetrahedralStereoAst::Undetermined, None)]
    #[case::stereo_open(TetrahedralStereoAst::Stereo(StereoCosetAst::Undetermined), None)]
    fn test_tetrahedral_stereo_ast_as_lit(#[case] site: TetrahedralStereoAst, #[case] expected: Option<StereoConfiguration>) {
        assert_eq!(site.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(TetrahedralStereoAst::Undetermined, true)]
    #[case::not_stereo(TetrahedralStereoAst::NotStereo, false)]
    #[case::stereo(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), false)]
    fn test_tetrahedral_stereo_ast_is_undetermined(#[case] site: TetrahedralStereoAst, #[case] expected: bool) {
        assert_eq!(site.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::not_stereo(TetrahedralStereoAst::NotStereo, true)]
    #[case::stereo_lit(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), true)]
    #[case::undetermined(TetrahedralStereoAst::Undetermined, false)]
    #[case::stereo_open(TetrahedralStereoAst::Stereo(StereoCosetAst::Undetermined), false)]
    fn test_tetrahedral_stereo_ast_is_ground(#[case] site: TetrahedralStereoAst, #[case] expected: bool) {
        assert_eq!(site.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoAst::Undetermined, TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), Some(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0))))]
    #[case::not_stereo_same(TetrahedralStereoAst::NotStereo, TetrahedralStereoAst::NotStereo, Some(TetrahedralStereoAst::NotStereo))]
    #[case::not_stereo_vs_stereo(TetrahedralStereoAst::NotStereo, TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), None)]
    #[case::stereo_same(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), Some(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0))))]
    #[case::stereo_disjoint(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1)), None)]
    #[case::open_narrows(TetrahedralStereoAst::Stereo(StereoCosetAst::Undetermined), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), Some(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0))))]
    fn test_tetrahedral_stereo_ast_meet(#[case] a: TetrahedralStereoAst, #[case] b: TetrahedralStereoAst, #[case] expected: Option<TetrahedralStereoAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoAst::Undetermined, TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Undetermined)]
    #[case::not_stereo_same(TetrahedralStereoAst::NotStereo, TetrahedralStereoAst::NotStereo, TetrahedralStereoAst::NotStereo)]
    #[case::not_stereo_vs_stereo(TetrahedralStereoAst::NotStereo, TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Undetermined)]
    #[case::stereo_same(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)))]
    #[case::stereo_widens(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1)), TetrahedralStereoAst::Stereo(StereoCosetAst::lit_set([0, 1])))]
    fn test_tetrahedral_stereo_ast_join(#[case] a: TetrahedralStereoAst, #[case] b: TetrahedralStereoAst, #[case] expected: TetrahedralStereoAst) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoAst::Undetermined, TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), true)]
    #[case::open_matches_lit(TetrahedralStereoAst::Stereo(StereoCosetAst::Undetermined), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), true)]
    #[case::specific_vs_undetermined(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Undetermined, false)]
    #[case::lit_match(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), true)]
    #[case::lit_mismatch(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1)), false)]
    #[case::not_stereo_match(TetrahedralStereoAst::NotStereo, TetrahedralStereoAst::NotStereo, true)]
    #[case::not_stereo_vs_stereo(TetrahedralStereoAst::NotStereo, TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), false)]
    fn test_tetrahedral_stereo_ast_matches(#[case] pattern: TetrahedralStereoAst, #[case] target: TetrahedralStereoAst, #[case] expected: bool) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoAst::Undetermined, 0, true)]
    #[case::not_stereo(TetrahedralStereoAst::NotStereo, 0, false)]
    #[case::stereo_match(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), 0, true)]
    #[case::stereo_miss(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)), 1, false)]
    fn test_tetrahedral_stereo_ast_matches_value(#[case] site: TetrahedralStereoAst, #[case] value: u32, #[case] expected: bool) {
        assert_eq!(site.matches_value(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::term_swap_folds(CisTransStereoAst::Stereo(StereoCosetAst::term(StereoTerm::swap(StereoTerm::Lit(0)))), CisTransStereoAst::Stereo(StereoCosetAst::Lit(1)))]
    fn test_cis_trans_stereo_ast_canonicalize(#[case] input: CisTransStereoAst, #[case] expected: CisTransStereoAst) {
        assert_eq!(input.canonicalize(), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereo_lit(CisTransStereoAst::Stereo(StereoCosetAst::Lit(0)), Some(StereoConfiguration { kind: StereoKind::CisTrans, coset: 0 }))]
    #[case::not_stereo(CisTransStereoAst::NotStereo, None)]
    fn test_cis_trans_stereo_ast_as_lit(#[case] site: CisTransStereoAst, #[case] expected: Option<StereoConfiguration>) {
        assert_eq!(site.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(CisTransStereoAst::Undetermined, CisTransStereoAst::Stereo(StereoCosetAst::Lit(0)), Some(CisTransStereoAst::Stereo(StereoCosetAst::Lit(0))))]
    #[case::not_stereo_vs_stereo(CisTransStereoAst::NotStereo, CisTransStereoAst::Stereo(StereoCosetAst::Lit(0)), None)]
    #[case::stereo_disjoint(CisTransStereoAst::Stereo(StereoCosetAst::Lit(0)), CisTransStereoAst::Stereo(StereoCosetAst::Lit(1)), None)]
    fn test_cis_trans_stereo_ast_meet(#[case] a: CisTransStereoAst, #[case] b: CisTransStereoAst, #[case] expected: Option<CisTransStereoAst>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCosetAst::Lit(2), Some(2))]
    #[case::undetermined(StereoCosetAst::Undetermined, None)]
    #[case::lit_set(StereoCosetAst::lit_set([1, 3]), None)]
    #[case::term(StereoCosetAst::term(StereoTerm::var("o")), None)]
    fn test_stereo_coset_ast_as_lit(#[case] coset: StereoCosetAst, #[case] expected: Option<u32>) {
        assert_eq!(coset.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_identity(StereoCosetAst::Lit(1), StereoKind::Tetrahedral, StereoCosetAst::Lit(1))]
    #[case::swap_lit_even(StereoCosetAst::term(StereoTerm::swap(StereoTerm::Lit(0))), StereoKind::Tetrahedral, StereoCosetAst::Lit(1))]
    #[case::swap_lit_odd(StereoCosetAst::term(StereoTerm::swap(StereoTerm::Lit(1))), StereoKind::Tetrahedral, StereoCosetAst::Lit(0))]
    #[case::mirror_chiral(StereoCosetAst::term(StereoTerm::mirror(StereoTerm::Lit(0))), StereoKind::Tetrahedral, StereoCosetAst::Lit(1))]
    #[case::mirror_achiral_noop(StereoCosetAst::term(StereoTerm::mirror(StereoTerm::Lit(0))), StereoKind::CisTrans, StereoCosetAst::Lit(0))]
    #[case::apply_lit(StereoCosetAst::term(StereoTerm::apply(StereoTerm::Lit(0), Permutation::from_image(4, &[1, 0, 2, 3]))), StereoKind::Tetrahedral, StereoCosetAst::Lit(1))]
    #[case::sp_swap_four(StereoCosetAst::term(StereoTerm::swap(StereoTerm::Lit(1))), StereoKind::SquarePlanar, StereoCosetAst::Lit(2))]
    #[case::swap_var_chiral_to_mirror(StereoCosetAst::term(StereoTerm::swap(StereoTerm::var("o"))), StereoKind::Tetrahedral, StereoCosetAst::term(StereoTerm::mirror(StereoTerm::var("o"))))]
    #[case::swap_var_achiral_stays(StereoCosetAst::term(StereoTerm::swap(StereoTerm::var("o"))), StereoKind::CisTrans, StereoCosetAst::term(StereoTerm::swap(StereoTerm::var("o"))))]
    #[case::multi_element_set_preserved(StereoCosetAst::lit_set([0, 1]), StereoKind::SquarePlanar, StereoCosetAst::lit_set([0, 1]))]
    #[case::singleton_set_to_lit(StereoCosetAst::lit_set([1]), StereoKind::Octahedral, StereoCosetAst::Lit(1))]
    fn test_canon_coset(#[case] coset: StereoCosetAst, #[case] kind: StereoKind, #[case] expected: StereoCosetAst) {
        assert_eq!(canon_coset(coset, kind), Ok(expected));
    }

    #[rstest]
    #[case::empty_set(StereoCosetAst::LitSet(BTreeSet::new()), StereoKind::SquarePlanar)]
    fn test_canon_coset_error(#[case] coset: StereoCosetAst, #[case] kind: StereoKind) {
        assert_eq!(canon_coset(coset, kind), Err(Contradiction));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCosetAst::Undetermined, StereoCosetAst::Lit(1), StereoKind::Tetrahedral, Some(StereoCosetAst::Lit(1)))]
    #[case::lit_same(StereoCosetAst::Lit(1), StereoCosetAst::Lit(1), StereoKind::Tetrahedral, Some(StereoCosetAst::Lit(1)))]
    #[case::lit_disjoint(StereoCosetAst::Lit(0), StereoCosetAst::Lit(1), StereoKind::Tetrahedral, None)]
    #[case::set_intersect(StereoCosetAst::lit_set([1, 3]), StereoCosetAst::lit_set([3, 5]), StereoKind::Octahedral, Some(StereoCosetAst::Lit(3)))]
    #[case::term_equal(StereoCosetAst::term(StereoTerm::var("o")), StereoCosetAst::term(StereoTerm::var("o")), StereoKind::Tetrahedral, Some(StereoCosetAst::term(StereoTerm::var("o"))))]
    #[case::term_distinct(StereoCosetAst::term(StereoTerm::var("o")), StereoCosetAst::term(StereoTerm::var("p")), StereoKind::Tetrahedral, None)]
    fn test_coset_meet(#[case] a: StereoCosetAst, #[case] b: StereoCosetAst, #[case] kind: StereoKind, #[case] expected: Option<StereoCosetAst>) {
        assert_eq!(coset_meet(&a, &b, kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCosetAst::Undetermined, StereoCosetAst::Lit(1), StereoKind::Tetrahedral, StereoCosetAst::Undetermined)]
    #[case::lit_same(StereoCosetAst::Lit(1), StereoCosetAst::Lit(1), StereoKind::Tetrahedral, StereoCosetAst::Lit(1))]
    #[case::lit_union(StereoCosetAst::Lit(0), StereoCosetAst::Lit(1), StereoKind::Tetrahedral, StereoCosetAst::lit_set([0, 1]))]
    #[case::set_union(StereoCosetAst::lit_set([1, 3]), StereoCosetAst::lit_set([3, 5]), StereoKind::Octahedral, StereoCosetAst::lit_set([1, 3, 5]))]
    fn test_coset_join(#[case] a: StereoCosetAst, #[case] b: StereoCosetAst, #[case] kind: StereoKind, #[case] expected: StereoCosetAst) {
        assert_eq!(coset_join(&a, &b, kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCosetAst::Undetermined, StereoCosetAst::Lit(1), StereoKind::Tetrahedral, true)]
    #[case::lit_match(StereoCosetAst::Lit(1), StereoCosetAst::Lit(1), StereoKind::Tetrahedral, true)]
    #[case::lit_miss(StereoCosetAst::Lit(0), StereoCosetAst::Lit(1), StereoKind::Tetrahedral, false)]
    #[case::set_member(StereoCosetAst::lit_set([1, 3]), StereoCosetAst::Lit(3), StereoKind::Octahedral, true)]
    #[case::set_nonmember(StereoCosetAst::lit_set([1, 3]), StereoCosetAst::Lit(2), StereoKind::Octahedral, false)]
    #[case::specific_vs_wildcard(StereoCosetAst::Lit(0), StereoCosetAst::Undetermined, StereoKind::Tetrahedral, false)]
    fn test_coset_matches(#[case] pattern: StereoCosetAst, #[case] target: StereoCosetAst, #[case] kind: StereoKind, #[case] expected: bool) {
        assert_eq!(coset_matches(&pattern, &target, kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCosetAst::Lit(0), Permutation::from_image(4, &[1, 0, 2, 3]), StereoKind::Tetrahedral, StereoCosetAst::Lit(1))]
    #[case::undetermined(StereoCosetAst::Undetermined, Permutation::from_image(4, &[1, 0, 2, 3]), StereoKind::Tetrahedral, StereoCosetAst::Undetermined)]
    fn test_coset_apply_permutation(#[case] coset: StereoCosetAst, #[case] permutation: Permutation, #[case] kind: StereoKind, #[case] expected: StereoCosetAst) {
        assert_eq!(coset_apply_permutation(&coset, permutation, kind), expected);
    }

    #[rstest]
    fn test_stereo_atom_ast_new() {
        let stereo_atom = StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined);
        assert_eq!(
            stereo_atom.configuration,
            StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCosetAst::Undetermined)
        );
        assert_eq!(stereo_atom.constraints, StereoAtomConstraints::new());
    }

    #[rstest]
    fn test_stereo_atom_ast_simplify_values() {
        let mut atom = StereoAtomAst::new(
            StereoKind::Tetrahedral,
            StereoCosetAst::term(StereoTerm::swap(StereoTerm::Lit(0))),
        );
        atom.simplify_values();
        assert_eq!(
            atom.configuration,
            StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(1))
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), false)]
    #[case::ground(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), true)]
    fn test_stereo_atom_ast_is_ground(#[case] atom: StereoAtomAst, #[case] expected: bool) {
        assert_eq!(atom.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::open_coset(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined))]
    #[case::ground(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32))]
    fn test_stereo_atom_ast_into_ground(#[case] atom: StereoAtomAst) {
        assert_eq!(atom.clone().into_ground(), atom);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32),
        Some(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32)))]
    #[case::different_kind(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined),
        StereoAtomAst::new(StereoKind::SquarePlanar, StereoCosetAst::Undetermined), None)]
    #[case::config_conflict(StereoAtomAst::new(StereoKind::Tetrahedral, 0u32), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), None)]
    fn test_stereo_atom_ast_meet(
        #[case] a: StereoAtomAst,
        #[case] b: StereoAtomAst,
        #[case] expected: Option<StereoAtomAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_coset(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32))]
    #[case::distinct_cosets_widen(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), StereoAtomAst::new(StereoKind::Tetrahedral, 2u32), StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::lit_set([1, 2])))]
    fn test_stereo_atom_ast_join(#[case] a: StereoAtomAst, #[case] b: StereoAtomAst, #[case] expected: StereoAtomAst) {
        assert_eq!(a.join(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_match(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined), StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), true)]
    #[case::different_kind(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32), StereoAtomAst::new(StereoKind::SquarePlanar, 1u32), false)]
    fn test_stereo_atom_ast_matches(
        #[case] pattern: StereoAtomAst,
        #[case] target: StereoAtomAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_coset(
        StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::lit_set([1])),
        Ok(StereoAtomAst::new(StereoKind::Tetrahedral, 1u32)),
    )]
    #[case::empty_coset_litset_contradiction(
        StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::lit_set(Vec::<u32>::new())),
        Err(Contradiction),
    )]
    fn test_stereo_atom_ast_canonicalize(
        #[case] input: StereoAtomAst,
        #[case] expected: Result<StereoAtomAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    fn test_stereo_bond_ast_new() {
        let stereo_bond = StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Undetermined);
        assert_eq!(
            stereo_bond.configuration,
            StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCosetAst::Undetermined)
        );
        assert_eq!(stereo_bond.constraints, StereoBondConstraints::new())
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Undetermined), StereoBondAst::new(StereoKind::CisTrans, 1u32),
        Some(StereoBondAst::new(StereoKind::CisTrans, 1u32)))]
    #[case::config_conflict(StereoBondAst::new(StereoKind::CisTrans, 0u32), StereoBondAst::new(StereoKind::CisTrans, 1u32), None)]
    fn test_stereo_bond_ast_meet(
        #[case] a: StereoBondAst,
        #[case] b: StereoBondAst,
        #[case] expected: Option<StereoBondAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_coset(
        StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::lit_set([1])),
        Ok(StereoBondAst::new(StereoKind::CisTrans, 1u32)),
    )]
    #[case::empty_coset_litset_contradiction(
        StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::lit_set(Vec::<u32>::new())),
        Err(Contradiction),
    )]
    fn test_stereo_bond_ast_canonicalize(
        #[case] input: StereoBondAst,
        #[case] expected: Result<StereoBondAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }
}
