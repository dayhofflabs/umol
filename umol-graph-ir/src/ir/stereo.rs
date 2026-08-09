//! Stereochemistry forms: the configuration value and the operator-expression
//! tree over it.
//!
//! A configuration value is a dense coset index per stereo kind, corresponds to OpenSMILES
//! numbering for SP, TB, and OH.
//! `~` and `^` are group actions on the index; the owning configuration's
//! `canonicalize` folds closed operator-expressions against the coset algebra.

use std::borrow::Cow;
use std::collections::BTreeSet;

use strum::VariantArray;
use umol_graph_core::{BiRelationData, ParticipantPosition};
use umol_graph_ir_macros::{Canonicalize, Lattice};
use umol_perm::{ClassKey, Permutation};

use super::constraint::{
    StereoAtomConstraintForm, StereoAtomConstraintsForm, StereoBondConstraintForm,
    StereoBondConstraintsForm,
};
use super::error::{Contradiction, NoJoin};
use super::ligand::StereoLigand;
use super::traits::{AsLit, Canonicalize, Lattice};

/// Defines the stereo entity forms.
macro_rules! stereo_element {
    (
        $(#[doc = $doc:literal])+
        $name:ident, $constraints:ident, $constraint:ident
    ) => {
        $(#[doc = $doc])+
        #[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
        pub struct $name {
            pub configuration: StereoConfigurationForm,
            pub constraints: $constraints,
        }

        impl BiRelationData for $name {
            /// The ligands are an `Ordered` factor, so `order_2` from `canonicalize_positions` is
            /// always the identity — the frame-relative configuration needs no reindex here.
            fn on_permutation(
                &mut self,
                _order_1: &[ParticipantPosition],
                _order_2: &[ParticipantPosition],
            ) {
            }

            fn is_permutation_invariant(&self) -> bool {
                true
            }
        }

        impl From<&str> for $name {
            fn from(s: &str) -> Self {
                s.parse()
                    .expect(concat!("invalid ", stringify!($name), " string"))
            }
        }

        impl $name {
            pub fn new(kind: StereoKind, coset: impl Into<StereoCoset>) -> Self {
                Self {
                    configuration: StereoConfigurationForm::kinded(kind, coset),
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

            /// Relabel the ligand positions (`^`); constraints are positionless and unchanged.
            pub fn apply(&self, permutation: Permutation) -> Self {
                Self {
                    configuration: self.configuration.apply(permutation),
                    constraints: self.constraints.clone(),
                }
            }

            /// The kind involution (`~`).
            pub fn swap(&self) -> Self {
                Self {
                    configuration: self.configuration.swap(),
                    constraints: self.constraints.clone(),
                }
            }

            /// The enantiomer / mirror (`'`).
            pub fn mirror(&self) -> Self {
                Self {
                    configuration: self.configuration.mirror(),
                    constraints: self.constraints.clone(),
                }
            }

            /// Restate the coset in the `after` ligand frame given it was stated in `before` — the
            /// coset action of the induced frame permutation. Mechanical bookkeeping for a relabel
            /// that leaves the physical configuration unchanged; self-inverting. `before`/`after`
            /// must be the same ligand multiset reordered (a genuine ligand-set change is membership,
            /// not a frame permutation).
            pub fn transform_frame(
                &self,
                before: &[StereoLigand],
                after: &[StereoLigand],
            ) -> Option<Self> {
                Permutation::between(before, after).map(|permutation| self.apply(permutation))
            }
        }

    };
}

stereo_element! {
    /// Stereo atom form with geometry class, configuration, and per-site constraints.
    StereoAtomForm, StereoAtomConstraintsForm, StereoAtomConstraintForm
}

stereo_element! {
    /// Stereo bond form with cis/trans configuration and per-site constraints.
    StereoBondForm, StereoBondConstraintsForm, StereoBondConstraintForm
}

/// Configuration portion of a stereo-element update.
///
/// `Unchanged` omits the field, `Undetermined` explicitly clears it, and
/// `Kinded` carries either an absolute coset or a kind-only relative update.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoConfigurationUpdate {
    #[default]
    Unchanged,
    Undetermined,
    Kinded {
        kind: StereoKind,
        coset: Option<StereoCoset>,
    },
}

impl StereoConfigurationUpdate {
    fn apply_to(&self, current: &StereoConfigurationForm) -> StereoConfigurationForm {
        match self {
            Self::Unchanged => current.clone(),
            Self::Undetermined => StereoConfigurationForm::Undetermined,
            Self::Kinded {
                kind,
                coset: Some(coset),
            } => StereoConfigurationForm::kinded(*kind, coset.clone()),
            Self::Kinded { kind, coset: None } => match current {
                StereoConfigurationForm::Kinded(current_kind, current_coset)
                    if current_kind == kind =>
                {
                    StereoConfigurationForm::kinded(*kind, current_coset.clone())
                }
                _ => StereoConfigurationForm::kinded(*kind, StereoCoset::Undetermined),
            },
        }
    }

    pub(crate) fn kind(&self) -> Option<StereoKind> {
        match self {
            Self::Kinded { kind, .. } => Some(*kind),
            Self::Unchanged | Self::Undetermined => None,
        }
    }
}

/// Attribute update for a stereo atom.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoAtomUpdate {
    pub configuration: StereoConfigurationUpdate,
    pub constraints: StereoAtomConstraintsForm,
}

/// Attribute update for a stereo bond.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoBondUpdate {
    pub configuration: StereoConfigurationUpdate,
    pub constraints: StereoBondConstraintsForm,
}

impl StereoAtomForm {
    /// Apply an attribute update.
    pub fn update(&self, update: &StereoAtomUpdate) -> Self {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        Self {
            configuration: update.configuration.apply_to(&self.configuration),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> StereoAtomUpdate {
        let mut constraints = StereoAtomConstraintsForm::new();
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
        let configuration = if self.configuration.canonical_eq(&other.configuration) {
            match self.configuration.kind() {
                Some(kind) if !constraints.is_empty() => {
                    StereoConfigurationUpdate::Kinded { kind, coset: None }
                }
                _ => StereoConfigurationUpdate::Unchanged,
            }
        } else {
            match &other.configuration {
                StereoConfigurationForm::Undetermined => StereoConfigurationUpdate::Undetermined,
                StereoConfigurationForm::Kinded(kind, coset) => StereoConfigurationUpdate::Kinded {
                    kind: *kind,
                    coset: Some(coset.clone()),
                },
            }
        };
        StereoAtomUpdate {
            configuration,
            constraints,
        }
    }
}

impl StereoBondForm {
    /// Apply an attribute update.
    pub fn update(&self, update: &StereoBondUpdate) -> Self {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        Self {
            configuration: update.configuration.apply_to(&self.configuration),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> StereoBondUpdate {
        let mut constraints = StereoBondConstraintsForm::new();
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
        let configuration = if self.configuration.canonical_eq(&other.configuration) {
            match self.configuration.kind() {
                Some(kind) if !constraints.is_empty() => {
                    StereoConfigurationUpdate::Kinded { kind, coset: None }
                }
                _ => StereoConfigurationUpdate::Unchanged,
            }
        } else {
            match &other.configuration {
                StereoConfigurationForm::Undetermined => StereoConfigurationUpdate::Undetermined,
                StereoConfigurationForm::Kinded(kind, coset) => StereoConfigurationUpdate::Kinded {
                    kind: *kind,
                    coset: Some(coset.clone()),
                },
            }
        };
        StereoBondUpdate {
            configuration,
            constraints,
        }
    }
}

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
        self.class_key().space().degree()
    }

    /// Number of cosets/configurations in this stereo kind.
    pub fn count(self) -> usize {
        self.class_key().space().count()
    }

    /// Whether this stereo kind can encode local handedness.
    pub fn is_chiral_class(self) -> bool {
        self.class_key().space().is_chiral()
    }

    /// Kind-specific `~` involution. Chiral kinds borrow the orientation-reversing
    /// generator from umol-perm; achiral kinds use a chosen ligand swap (no improper
    /// generator to borrow — theirs is the identity):
    /// - cis/trans: swap the two configurations
    /// - square-planar: swap the diagonal ligand pair
    pub fn involution(self) -> Permutation {
        let coset_space = self.class_key().space();
        if coset_space.is_chiral() {
            coset_space.improper()
        } else {
            match self {
                StereoKind::CisTrans => Permutation::from_image(&[1, 0, 2, 3]),
                StereoKind::SquarePlanar => Permutation::from_image(&[2, 1, 0, 3]),
                _ => unreachable!("only achiral kinds reach the chosen-swap branch"),
            }
        }
    }

    /// Act on coset index `index` by `permutation`, through the class's coset algebra.
    pub fn act(self, index: u32, permutation: Permutation) -> u32 {
        self.class_key()
            .space()
            .reindex(index, permutation)
            .expect("act: valid coset index and permutation")
    }

    /// The mirror (improper, μ) generator as a permutation: chiral kinds use the
    /// orientation-reversing generator; achiral kinds act trivially on cosets.
    pub fn mirror_permutation(self) -> Permutation {
        if self.is_chiral_class() {
            self.class_key().space().improper()
        } else {
            Permutation::identity(self.degree())
        }
    }

    /// Whether `g` and `h` induce the same coset permutation for this kind.
    fn coset_action_eq(self, g: Permutation, h: Permutation) -> bool {
        let s = self.class_key().space();
        (0..s.count() as u32).all(|i| s.reindex(i, g) == s.reindex(i, h))
    }

    /// Canonicalize coset permutation, priority `Mirror > Swap > Apply`; `None`
    /// when it acts as the identity on cosets.
    pub fn canonicalize_permutation(self, g: Permutation) -> Option<CosetOp> {
        if self.coset_action_eq(g, Permutation::identity(self.degree())) {
            None
        } else if self.is_chiral_class() && self.coset_action_eq(g, self.mirror_permutation()) {
            Some(CosetOp::Mirror)
        } else if self.coset_action_eq(g, self.involution()) {
            Some(CosetOp::Swap)
        } else {
            Some(CosetOp::Apply(g))
        }
    }
}

/// Permutation in canonical priority form, `Mirror` > `Swap` > `Apply`, kind-dependent.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CosetOp {
    Swap,
    Mirror,
    Apply(Permutation),
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
pub enum StereoConfigurationForm {
    #[default]
    Undetermined,
    Kinded(StereoKind, StereoCoset),
}

impl StereoConfigurationForm {
    pub fn undetermined() -> Self {
        Self::Undetermined
    }

    pub fn kinded(kind: StereoKind, coset: impl Into<StereoCoset>) -> Self {
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
    pub fn coset(&self) -> Option<&StereoCoset> {
        match self {
            Self::Kinded(_, coset) => Some(coset),
            Self::Undetermined => None,
        }
    }

    /// Mutable access to the coset, or `None` when undetermined.
    pub fn coset_mut(&mut self) -> Option<&mut StereoCoset> {
        match self {
            Self::Kinded(_, coset) => Some(coset),
            Self::Undetermined => None,
        }
    }

    /// Relabel the ligand positions (`^`); `Undetermined` is fixed.
    pub fn apply(&self, permutation: Permutation) -> Self {
        self.map_kinded(|kind, coset| coset.apply(kind, permutation))
    }

    /// The kind involution (`~`).
    pub fn swap(&self) -> Self {
        self.map_kinded(|kind, coset| coset.swap(kind))
    }

    /// The enantiomer / mirror (`'`).
    pub fn mirror(&self) -> Self {
        self.map_kinded(|kind, coset| coset.mirror(kind))
    }

    fn map_kinded(&self, f: impl FnOnce(StereoKind, &StereoCoset) -> StereoCoset) -> Self {
        match self {
            Self::Undetermined => Self::Undetermined,
            Self::Kinded(kind, coset) => Self::Kinded(*kind, f(*kind, coset)),
        }
    }

    /// Overwrite with `other`, field-wise: an `Undetermined` `other` keeps `self`; a same-kind
    /// `other` with an `Undetermined` coset keeps `self`'s coset (a partial that fixes only the
    /// kind); a differing kind or a determined coset overrides wholesale.
    pub fn update(&self, other: &Self) -> Self {
        match (self, other) {
            (_, Self::Undetermined) => self.clone(),
            (Self::Kinded(ks, cs), Self::Kinded(ko, StereoCoset::Undetermined)) if ks == ko => {
                Self::Kinded(*ks, cs.clone())
            }
            (_, Self::Kinded(..)) => other.clone(),
        }
    }
}

impl From<(StereoKind, u32)> for StereoConfigurationForm {
    fn from((kind, coset): (StereoKind, u32)) -> Self {
        Self::Kinded(kind, StereoCoset::Lit(coset))
    }
}

impl Canonicalize for StereoConfigurationForm {
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

impl AsLit for StereoConfigurationForm {
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

impl Lattice for StereoConfigurationForm {
    fn is_undetermined(&self) -> bool {
        matches!(self, Self::Undetermined)
    }

    fn is_ground(&self) -> bool {
        matches!(self, Self::Kinded(_, StereoCoset::Lit(_)))
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

    fn join(&self, other: &Self) -> Result<Self, NoJoin> {
        let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
        Ok(match (a.as_ref(), b.as_ref()) {
            (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
            (Self::Kinded(k1, ca), Self::Kinded(k2, cb)) => {
                if k1 != k2 {
                    Self::Undetermined
                } else {
                    Self::Kinded(*k1, coset_join(ca, cb, *k1))
                }
            }
        })
    }
}

/// Generates a constraint-side stereo state for a fixed geometry (`#T`/`#C`):
/// undetermined, explicitly not-stereo, or a stereo center with a coset. The
/// geometry is the type's identity (`$kind`), so the coset folds/meets under that
/// constant kind — no kind field.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TetrahedralStereo {
    NotStereo,
    Stereo(u32),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CisTransStereo {
    NotStereo,
    Stereo(u32),
}

/// Named tetrahedral configurations and their canonical coset indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TetrahedralConfiguration {
    Ccw,
    Cw,
}

/// Named cis/trans configurations and their canonical coset indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CisTransConfiguration {
    Z,
    E,
}

macro_rules! stereo_site {
    ($name:ident, $lit:ident, $kind:expr) => {
        #[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub enum $name {
            #[default]
            Undetermined,
            NotStereo,
            Stereo(StereoCoset),
        }

        impl $name {
            pub fn undetermined() -> Self {
                Self::Undetermined
            }

            pub fn not_stereo() -> Self {
                Self::NotStereo
            }

            pub fn stereo(coset: impl Into<StereoCoset>) -> Self {
                Self::Stereo(coset.into())
            }

            pub fn is_stereo(&self) -> bool {
                matches!(self, Self::Stereo(_))
            }

            /// Matches literal coset index `value` under the type's kind.
            pub fn matches_value(&self, value: u32) -> bool {
                match self {
                    Self::Stereo(coset) => coset_matches(coset, &StereoCoset::Lit(value), $kind),
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
            type Lit = $lit;

            /// The exact absence or stereo-coset value when ground.
            fn as_lit(&self) -> Option<$lit> {
                match self {
                    Self::NotStereo => Some($lit::NotStereo),
                    Self::Stereo(StereoCoset::Lit(coset)) => Some($lit::Stereo(*coset)),
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
                    Self::Stereo(coset) => matches!(coset, StereoCoset::Lit(_)),
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

            fn join(&self, other: &Self) -> Result<Self, NoJoin> {
                let a = self.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
                let b = other.canonical().unwrap_or(Cow::Owned(Self::Undetermined));
                Ok(match (a.as_ref(), b.as_ref()) {
                    (Self::Undetermined, _) | (_, Self::Undetermined) => Self::Undetermined,
                    (Self::NotStereo, Self::NotStereo) => Self::NotStereo,
                    (Self::NotStereo, Self::Stereo(_)) | (Self::Stereo(_), Self::NotStereo) => {
                        Self::Undetermined
                    }
                    (Self::Stereo(ca), Self::Stereo(cb)) => Self::Stereo(coset_join(ca, cb, $kind)),
                })
            }
        }
    };
}

stereo_site! { TetrahedralStereoForm, TetrahedralStereo, StereoKind::Tetrahedral }
stereo_site! { CisTransStereoForm, CisTransStereo, StereoKind::CisTrans }

impl From<TetrahedralStereo> for TetrahedralStereoForm {
    fn from(stereo: TetrahedralStereo) -> Self {
        match stereo {
            TetrahedralStereo::NotStereo => Self::NotStereo,
            TetrahedralStereo::Stereo(coset) => Self::Stereo(StereoCoset::Lit(coset)),
        }
    }
}

impl From<CisTransStereo> for CisTransStereoForm {
    fn from(stereo: CisTransStereo) -> Self {
        match stereo {
            CisTransStereo::NotStereo => Self::NotStereo,
            CisTransStereo::Stereo(coset) => Self::Stereo(StereoCoset::Lit(coset)),
        }
    }
}

impl From<TetrahedralConfiguration> for TetrahedralStereo {
    fn from(configuration: TetrahedralConfiguration) -> Self {
        match configuration {
            TetrahedralConfiguration::Ccw => Self::Stereo(0),
            TetrahedralConfiguration::Cw => Self::Stereo(1),
        }
    }
}

impl From<CisTransConfiguration> for CisTransStereo {
    fn from(configuration: CisTransConfiguration) -> Self {
        match configuration {
            CisTransConfiguration::Z => Self::Stereo(0),
            CisTransConfiguration::E => Self::Stereo(1),
        }
    }
}

impl From<TetrahedralConfiguration> for TetrahedralStereoForm {
    fn from(configuration: TetrahedralConfiguration) -> Self {
        TetrahedralStereo::from(configuration).into()
    }
}

impl From<CisTransConfiguration> for CisTransStereoForm {
    fn from(configuration: CisTransConfiguration) -> Self {
        CisTransStereo::from(configuration).into()
    }
}

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
/// single index, a finite set, or an operator `Term` over a variable.
/// Kind-relative — no `Lattice` or `Canonicalize`; the owning configuration or
/// site normalizes it under its kind.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoCoset {
    #[default]
    Undetermined,
    Lit(u32),
    LitSet(BTreeSet<u32>),
    Term(Box<StereoTerm>),
}

impl StereoCoset {
    pub fn lit_set(values: impl IntoIterator<Item = u32>) -> Self {
        Self::LitSet(values.into_iter().collect())
    }

    pub fn term(term: StereoTerm) -> Self {
        Self::Term(Box::new(term))
    }

    /// Relabel the ligand positions (the `^` op): move each literal coset index through the kind's
    /// coset algebra, eager on `Lit`/`LitSet`; an open `Term` keeps the operator layer.
    fn apply(&self, kind: StereoKind, permutation: Permutation) -> Self {
        self.map_index(
            |c| kind.act(c, permutation),
            |t| StereoTerm::apply(t, permutation),
        )
    }

    /// The kind involution (the `~` op).
    fn swap(&self, kind: StereoKind) -> Self {
        self.map_index(|c| kind.act(c, kind.involution()), StereoTerm::swap)
    }

    /// The enantiomer / mirror (the `'` op).
    fn mirror(&self, kind: StereoKind) -> Self {
        self.map_index(
            |c| kind.act(c, kind.mirror_permutation()),
            StereoTerm::mirror,
        )
    }

    /// Map each literal index by `lit`; an open `Term` is wrapped by `term` (the only case that keeps
    /// an operator layer — a bare variable cannot be evaluated). `Undetermined` is fixed.
    fn map_index(
        &self,
        lit: impl Fn(u32) -> u32,
        term: impl FnOnce(StereoTerm) -> StereoTerm,
    ) -> Self {
        match self {
            Self::Undetermined => Self::Undetermined,
            Self::Lit(c) => Self::Lit(lit(*c)),
            Self::LitSet(s) => Self::LitSet(s.iter().map(|&c| lit(c)).collect()),
            Self::Term(t) => Self::term(term((**t).clone())),
        }
    }
}

impl From<u32> for StereoCoset {
    fn from(index: u32) -> Self {
        Self::Lit(index)
    }
}

impl AsLit for StereoCoset {
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
/// `AsLit` target of `StereoConfigurationForm` and the per-kind site types.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct StereoConfiguration {
    pub kind: StereoKind,
    pub coset: u32,
}

/// The literal coset-index set a positive coset denotes; `None` for the wildcard
/// `Undetermined` and the symbolic `Term`. Used by `coset_meet`/`coset_join`
/// after those two cases are handled.
fn coset_to_set(coset: &StereoCoset) -> Option<BTreeSet<u32>> {
    match coset {
        StereoCoset::Lit(i) => Some(BTreeSet::from([*i])),
        StereoCoset::LitSet(s) => Some(s.clone()),
        StereoCoset::Undetermined | StereoCoset::Term(_) => None,
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
    coset: StereoCoset,
    kind: StereoKind,
) -> Result<StereoCoset, Contradiction> {
    let s = kind.class_key().space();
    let set: BTreeSet<u32> = match &coset {
        StereoCoset::Undetermined => return Ok(StereoCoset::Undetermined),
        StereoCoset::Lit(i) => BTreeSet::from([*i]),
        StereoCoset::LitSet(values) => values.clone(),
        StereoCoset::Term(t) => {
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
                    let term = match kind.canonicalize_permutation(g) {
                        None => var,
                        Some(CosetOp::Mirror) => StereoTerm::Mirror(Box::new(var)),
                        Some(CosetOp::Swap) => StereoTerm::Swap(Box::new(var)),
                        Some(CosetOp::Apply(g)) => StereoTerm::Apply(Box::new(var), g),
                    };
                    return Ok(StereoCoset::term(term));
                }
                StereoTerm::Lit(i) => BTreeSet::from([s.reindex(*i, g).ok_or(Contradiction)?]),
                StereoTerm::LitSet(values) => values
                    .iter()
                    .map(|i| s.reindex(*i, g).ok_or(Contradiction))
                    .collect::<Result<_, _>>()?,
                StereoTerm::Swap(_) | StereoTerm::Mirror(_) | StereoTerm::Apply(..) => {
                    unreachable!("compose_term returns a base leaf")
                }
            }
        }
    };
    if set.is_empty() {
        Err(Contradiction)
    } else if set.len() == 1 {
        Ok(StereoCoset::Lit(set.into_iter().next().unwrap()))
    } else {
        Ok(StereoCoset::LitSet(set))
    }
}

/// Greatest lower bound of two cosets under `kind` (canonicalizing operands);
/// `Term` meets only an equal canonical `Term`.
pub(crate) fn coset_meet(
    a: &StereoCoset,
    b: &StereoCoset,
    kind: StereoKind,
) -> Option<StereoCoset> {
    let ca = canon_coset(a.clone(), kind).ok()?;
    let cb = canon_coset(b.clone(), kind).ok()?;
    use StereoCoset::{Term, Undetermined};
    match (&ca, &cb) {
        (Undetermined, _) => Some(cb),
        (_, Undetermined) => Some(ca),
        (Term(_), Term(_)) => (ca == cb).then_some(ca),
        (Term(_), _) | (_, Term(_)) => None,
        _ => {
            let sa = coset_to_set(&ca).unwrap();
            let sb = coset_to_set(&cb).unwrap();
            canon_coset(
                StereoCoset::LitSet(sa.intersection(&sb).copied().collect()),
                kind,
            )
            .ok()
        }
    }
}

/// Least upper bound of two cosets under `kind`.
pub(crate) fn coset_join(a: &StereoCoset, b: &StereoCoset, kind: StereoKind) -> StereoCoset {
    let ca = canon_coset(a.clone(), kind).unwrap_or(StereoCoset::Undetermined);
    let cb = canon_coset(b.clone(), kind).unwrap_or(StereoCoset::Undetermined);
    use StereoCoset::{Term, Undetermined};
    match (&ca, &cb) {
        (Undetermined, _) | (_, Undetermined) => StereoCoset::Undetermined,
        (Term(_), Term(_)) if ca == cb => ca,
        (Term(_), _) | (_, Term(_)) => StereoCoset::Undetermined,
        _ => {
            let sa = coset_to_set(&ca).unwrap();
            let sb = coset_to_set(&cb).unwrap();
            canon_coset(StereoCoset::LitSet(sa.union(&sb).copied().collect()), kind)
                .unwrap_or(StereoCoset::Undetermined)
        }
    }
}

/// `target` refines `pattern` under `kind` (meet-derived).
pub(crate) fn coset_matches(pattern: &StereoCoset, target: &StereoCoset, kind: StereoKind) -> bool {
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
    coset: &StereoCoset,
    permutation: Permutation,
    kind: StereoKind,
) -> Option<StereoCoset> {
    let s = kind.class_key().space();
    match coset {
        StereoCoset::Undetermined => Some(StereoCoset::Undetermined),
        StereoCoset::Lit(i) => Some(StereoCoset::Lit(s.reindex(*i, permutation)?)),
        StereoCoset::LitSet(set) => Some(StereoCoset::LitSet(
            set.iter()
                .map(|i| s.reindex(*i, permutation))
                .collect::<Option<_>>()?,
        )),
        StereoCoset::Term(t) => Some(
            canon_coset(
                StereoCoset::term(StereoTerm::apply((**t).clone(), permutation)),
                kind,
            )
            .unwrap_or(StereoCoset::Undetermined),
        ),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::super::constraint::StereogenicityForm;
    use super::super::id::AtomId;
    use super::super::ligand::StereoLigandKind;
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
    #[case::tetrahedral((StereoKind::Tetrahedral, 1), StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)))]
    #[case::octahedral((StereoKind::Octahedral, 5), StereoConfigurationForm::Kinded(StereoKind::Octahedral, StereoCoset::Lit(5)))]
    fn test_stereo_configuration_form_from(#[case] input: (StereoKind, u32), #[case] expected: StereoConfigurationForm) {
        assert_eq!(StereoConfigurationForm::from(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::term_swap_folds_to_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0)))), StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)))]
    fn test_stereo_configuration_form_canonicalize(#[case] input: StereoConfigurationForm, #[case] expected: StereoConfigurationForm) {
        assert_eq!(input.canonicalize(), Ok(expected));
    }

    #[rstest]
    #[case::undetermined(StereoConfigurationForm::Undetermined)]
    #[case::kind_lit(StereoConfigurationForm::Kinded(
        StereoKind::Tetrahedral,
        StereoCoset::Lit(0)
    ))]
    #[case::kind_open(StereoConfigurationForm::Kinded(
        StereoKind::Tetrahedral,
        StereoCoset::Undetermined
    ))]
    // Multi-element / full coset sets are preserved (no complement or full→Undetermined fold).
    #[case::multi_element_set(StereoConfigurationForm::Kinded(StereoKind::SquarePlanar, StereoCoset::lit_set([0, 1])))]
    #[case::full_set(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::lit_set([0, 1])))]
    fn test_stereo_configuration_form_canonicalize_identity(
        #[case] input: StereoConfigurationForm,
    ) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rstest]
    #[case::empty_set(StereoConfigurationForm::Kinded(
        StereoKind::SquarePlanar,
        StereoCoset::LitSet(BTreeSet::new())
    ))]
    fn test_stereo_configuration_form_canonicalize_error(#[case] input: StereoConfigurationForm) {
        assert_eq!(input.canonicalize(), Err(Contradiction));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)), Some(StereoConfiguration { kind: StereoKind::Tetrahedral, coset: 1 }))]
    #[case::undetermined(StereoConfigurationForm::Undetermined, None)]
    #[case::kind_open(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), None)]
    fn test_stereo_configuration_form_as_lit(#[case] config: StereoConfigurationForm, #[case] expected: Option<StereoConfiguration>) {
        assert_eq!(config.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationForm::Undetermined, true)]
    #[case::kind_open(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), false)]
    #[case::kind_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)), false)]
    fn test_stereo_configuration_form_is_undetermined(#[case] config: StereoConfigurationForm, #[case] expected: bool) {
        assert_eq!(config.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)), true)]
    #[case::undetermined(StereoConfigurationForm::Undetermined, false)]
    #[case::kind_open(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), false)]
    fn test_stereo_configuration_form_is_ground(#[case] config: StereoConfigurationForm, #[case] expected: bool) {
        assert_eq!(config.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_narrows(StereoConfigurationForm::Undetermined, StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), Some(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0))))]
    #[case::coset_same(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), Some(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0))))]
    #[case::open_narrows(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), Some(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0))))]
    #[case::coset_conflict(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 1)), None)]
    #[case::kind_conflict(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::CisTrans, 0)), None)]
    fn test_stereo_configuration_form_meet(#[case] a: StereoConfigurationForm, #[case] b: StereoConfigurationForm, #[case] expected: Option<StereoConfigurationForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_absorbs(StereoConfigurationForm::Undetermined, StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::Undetermined)]
    #[case::coset_same(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)))]
    #[case::coset_widens(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 1)), StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::lit_set([0, 1])))]
    #[case::kind_conflict(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::CisTrans, 0)), StereoConfigurationForm::Undetermined)]
    fn test_stereo_configuration_form_join(#[case] a: StereoConfigurationForm, #[case] b: StereoConfigurationForm, #[case] expected: StereoConfigurationForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined_matches_any(StereoConfigurationForm::Undetermined, StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), true)]
    #[case::open_matches_lit(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), true)]
    #[case::specific_vs_undetermined(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::Undetermined, false)]
    #[case::coset_match(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), true)]
    #[case::coset_mismatch(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::Tetrahedral, 1)), false)]
    #[case::kind_mismatch(StereoConfigurationForm::from((StereoKind::Tetrahedral, 0)), StereoConfigurationForm::from((StereoKind::CisTrans, 0)), false)]
    fn test_stereo_configuration_form_matches(#[case] pattern: StereoConfigurationForm, #[case] target: StereoConfigurationForm, #[case] expected: bool) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::not_stereo(TetrahedralStereo::NotStereo, TetrahedralStereoForm::NotStereo)]
    #[case::stereo(
        TetrahedralStereo::Stereo(1),
        TetrahedralStereoForm::Stereo(StereoCoset::Lit(1))
    )]
    fn test_tetrahedral_stereo_form_from(
        #[case] stereo: TetrahedralStereo,
        #[case] expected: TetrahedralStereoForm,
    ) {
        assert_eq!(TetrahedralStereoForm::from(stereo), expected);
    }

    #[rstest]
    #[case::ccw(
        TetrahedralConfiguration::Ccw,
        TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))
    )]
    #[case::cw(
        TetrahedralConfiguration::Cw,
        TetrahedralStereoForm::Stereo(StereoCoset::Lit(1))
    )]
    fn test_tetrahedral_stereo_form_from_configuration(
        #[case] configuration: TetrahedralConfiguration,
        #[case] expected: TetrahedralStereoForm,
    ) {
        assert_eq!(TetrahedralStereoForm::from(configuration), expected);
    }

    #[rstest]
    #[case::not_stereo(CisTransStereo::NotStereo, CisTransStereoForm::NotStereo)]
    #[case::stereo(
        CisTransStereo::Stereo(1),
        CisTransStereoForm::Stereo(StereoCoset::Lit(1))
    )]
    fn test_cis_trans_stereo_form_from(
        #[case] stereo: CisTransStereo,
        #[case] expected: CisTransStereoForm,
    ) {
        assert_eq!(CisTransStereoForm::from(stereo), expected);
    }

    #[rstest]
    #[case::z(
        CisTransConfiguration::Z,
        CisTransStereoForm::Stereo(StereoCoset::Lit(0))
    )]
    #[case::e(
        CisTransConfiguration::E,
        CisTransStereoForm::Stereo(StereoCoset::Lit(1))
    )]
    fn test_cis_trans_stereo_form_from_configuration(
        #[case] configuration: CisTransConfiguration,
        #[case] expected: CisTransStereoForm,
    ) {
        assert_eq!(CisTransStereoForm::from(configuration), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::term_swap_folds(TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0)))), TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)))]
    fn test_tetrahedral_stereo_form_canonicalize(#[case] input: TetrahedralStereoForm, #[case] expected: TetrahedralStereoForm) {
        assert_eq!(input.canonicalize(), Ok(expected));
    }

    #[rstest]
    #[case::undetermined(TetrahedralStereoForm::Undetermined)]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo)]
    #[case::stereo_lit(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)))]
    #[case::stereo_open(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined))]
    fn test_tetrahedral_stereo_form_canonicalize_identity(#[case] input: TetrahedralStereoForm) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereo_lit(TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)), Some(TetrahedralStereo::Stereo(1)))]
    #[case::stereo_zero(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), Some(TetrahedralStereo::Stereo(0)))]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo, Some(TetrahedralStereo::NotStereo))]
    #[case::undetermined(TetrahedralStereoForm::Undetermined, None)]
    #[case::stereo_open(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined), None)]
    fn test_tetrahedral_stereo_form_as_lit(#[case] site: TetrahedralStereoForm, #[case] expected: Option<TetrahedralStereo>) {
        assert_eq!(site.as_lit(), expected);
        assert_eq!(site.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(TetrahedralStereoForm::Undetermined, true)]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo, false)]
    #[case::stereo(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), false)]
    fn test_tetrahedral_stereo_form_is_undetermined(#[case] site: TetrahedralStereoForm, #[case] expected: bool) {
        assert_eq!(site.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo, true)]
    #[case::stereo_lit(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), true)]
    #[case::undetermined(TetrahedralStereoForm::Undetermined, false)]
    #[case::stereo_open(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined), false)]
    fn test_tetrahedral_stereo_form_is_ground(#[case] site: TetrahedralStereoForm, #[case] expected: bool) {
        assert_eq!(site.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoForm::Undetermined, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), Some(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))))]
    #[case::not_stereo_same(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::NotStereo, Some(TetrahedralStereoForm::NotStereo))]
    #[case::not_stereo_vs_stereo(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), None)]
    #[case::stereo_same(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), Some(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))))]
    #[case::stereo_disjoint(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)), None)]
    #[case::open_narrows(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), Some(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))))]
    fn test_tetrahedral_stereo_form_meet(#[case] a: TetrahedralStereoForm, #[case] b: TetrahedralStereoForm, #[case] expected: Option<TetrahedralStereoForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoForm::Undetermined, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Undetermined)]
    #[case::not_stereo_same(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::NotStereo)]
    #[case::not_stereo_vs_stereo(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Undetermined)]
    #[case::stereo_same(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)))]
    #[case::stereo_widens(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)), TetrahedralStereoForm::Stereo(StereoCoset::lit_set([0, 1])))]
    fn test_tetrahedral_stereo_form_join(#[case] a: TetrahedralStereoForm, #[case] b: TetrahedralStereoForm, #[case] expected: TetrahedralStereoForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoForm::Undetermined, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), true)]
    #[case::open_matches_lit(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), true)]
    #[case::specific_vs_undetermined(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Undetermined, false)]
    #[case::lit_match(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), true)]
    #[case::lit_mismatch(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)), false)]
    #[case::not_stereo_match(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::NotStereo, true)]
    #[case::not_stereo_vs_stereo(TetrahedralStereoForm::NotStereo, TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), false)]
    fn test_tetrahedral_stereo_form_matches(#[case] pattern: TetrahedralStereoForm, #[case] target: TetrahedralStereoForm, #[case] expected: bool) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(TetrahedralStereoForm::Undetermined, 0, true)]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo, 0, false)]
    #[case::stereo_match(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), 0, true)]
    #[case::stereo_miss(TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)), 1, false)]
    fn test_tetrahedral_stereo_form_matches_value(#[case] site: TetrahedralStereoForm, #[case] value: u32, #[case] expected: bool) {
        assert_eq!(site.matches_value(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::term_swap_folds(CisTransStereoForm::Stereo(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0)))), CisTransStereoForm::Stereo(StereoCoset::Lit(1)))]
    fn test_cis_trans_stereo_form_canonicalize(#[case] input: CisTransStereoForm, #[case] expected: CisTransStereoForm) {
        assert_eq!(input.canonicalize(), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stereo_zero(CisTransStereoForm::Stereo(StereoCoset::Lit(0)), Some(CisTransStereo::Stereo(0)))]
    #[case::stereo_lit(CisTransStereoForm::Stereo(StereoCoset::Lit(1)), Some(CisTransStereo::Stereo(1)))]
    #[case::not_stereo(CisTransStereoForm::NotStereo, Some(CisTransStereo::NotStereo))]
    #[case::undetermined(CisTransStereoForm::Undetermined, None)]
    #[case::stereo_open(CisTransStereoForm::Stereo(StereoCoset::Undetermined), None)]
    fn test_cis_trans_stereo_form_as_lit(#[case] site: CisTransStereoForm, #[case] expected: Option<CisTransStereo>) {
        assert_eq!(site.as_lit(), expected);
        assert_eq!(site.is_ground(), expected.is_some());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(CisTransStereoForm::Undetermined, CisTransStereoForm::Stereo(StereoCoset::Lit(0)), Some(CisTransStereoForm::Stereo(StereoCoset::Lit(0))))]
    #[case::not_stereo_vs_stereo(CisTransStereoForm::NotStereo, CisTransStereoForm::Stereo(StereoCoset::Lit(0)), None)]
    #[case::stereo_disjoint(CisTransStereoForm::Stereo(StereoCoset::Lit(0)), CisTransStereoForm::Stereo(StereoCoset::Lit(1)), None)]
    fn test_cis_trans_stereo_form_meet(#[case] a: CisTransStereoForm, #[case] b: CisTransStereoForm, #[case] expected: Option<CisTransStereoForm>) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCoset::Lit(2), Some(2))]
    #[case::undetermined(StereoCoset::Undetermined, None)]
    #[case::lit_set(StereoCoset::lit_set([1, 3]), None)]
    #[case::term(StereoCoset::term(StereoTerm::var("o")), None)]
    fn test_stereo_coset_as_lit(#[case] coset: StereoCoset, #[case] expected: Option<u32>) {
        assert_eq!(coset.as_lit(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_identity(StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::swap_lit_even(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(0))), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::swap_lit_odd(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))), StereoKind::Tetrahedral, StereoCoset::Lit(0))]
    #[case::mirror_chiral(StereoCoset::term(StereoTerm::mirror(StereoTerm::Lit(0))), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::mirror_achiral_noop(StereoCoset::term(StereoTerm::mirror(StereoTerm::Lit(0))), StereoKind::CisTrans, StereoCoset::Lit(0))]
    #[case::apply_lit(StereoCoset::term(StereoTerm::apply(StereoTerm::Lit(0), Permutation::from_image(&[1, 0, 2, 3]))), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::sp_swap_four(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))), StereoKind::SquarePlanar, StereoCoset::Lit(2))]
    #[case::swap_var_chiral_to_mirror(StereoCoset::term(StereoTerm::swap(StereoTerm::var("o"))), StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::mirror(StereoTerm::var("o"))))]
    #[case::swap_var_achiral_stays(StereoCoset::term(StereoTerm::swap(StereoTerm::var("o"))), StereoKind::CisTrans, StereoCoset::term(StereoTerm::swap(StereoTerm::var("o"))))]
    #[case::multi_element_set_preserved(StereoCoset::lit_set([0, 1]), StereoKind::SquarePlanar, StereoCoset::lit_set([0, 1]))]
    #[case::singleton_set_to_lit(StereoCoset::lit_set([1]), StereoKind::Octahedral, StereoCoset::Lit(1))]
    fn test_canon_coset(#[case] coset: StereoCoset, #[case] kind: StereoKind, #[case] expected: StereoCoset) {
        assert_eq!(canon_coset(coset, kind), Ok(expected));
    }

    #[rstest]
    #[case::empty_set(StereoCoset::LitSet(BTreeSet::new()), StereoKind::SquarePlanar)]
    fn test_canon_coset_error(#[case] coset: StereoCoset, #[case] kind: StereoKind) {
        assert_eq!(canon_coset(coset, kind), Err(Contradiction));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCoset::Undetermined, StereoCoset::Lit(1), StereoKind::Tetrahedral, Some(StereoCoset::Lit(1)))]
    #[case::lit_same(StereoCoset::Lit(1), StereoCoset::Lit(1), StereoKind::Tetrahedral, Some(StereoCoset::Lit(1)))]
    #[case::lit_disjoint(StereoCoset::Lit(0), StereoCoset::Lit(1), StereoKind::Tetrahedral, None)]
    #[case::set_intersect(StereoCoset::lit_set([1, 3]), StereoCoset::lit_set([3, 5]), StereoKind::Octahedral, Some(StereoCoset::Lit(3)))]
    #[case::term_equal(StereoCoset::term(StereoTerm::var("o")), StereoCoset::term(StereoTerm::var("o")), StereoKind::Tetrahedral, Some(StereoCoset::term(StereoTerm::var("o"))))]
    #[case::term_distinct(StereoCoset::term(StereoTerm::var("o")), StereoCoset::term(StereoTerm::var("p")), StereoKind::Tetrahedral, None)]
    fn test_coset_meet(#[case] a: StereoCoset, #[case] b: StereoCoset, #[case] kind: StereoKind, #[case] expected: Option<StereoCoset>) {
        assert_eq!(coset_meet(&a, &b, kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCoset::Undetermined, StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::Undetermined)]
    #[case::lit_same(StereoCoset::Lit(1), StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::lit_union(StereoCoset::Lit(0), StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::lit_set([0, 1]))]
    #[case::set_union(StereoCoset::lit_set([1, 3]), StereoCoset::lit_set([3, 5]), StereoKind::Octahedral, StereoCoset::lit_set([1, 3, 5]))]
    fn test_coset_join(#[case] a: StereoCoset, #[case] b: StereoCoset, #[case] kind: StereoKind, #[case] expected: StereoCoset) {
        assert_eq!(coset_join(&a, &b, kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(StereoCoset::Undetermined, StereoCoset::Lit(1), StereoKind::Tetrahedral, true)]
    #[case::lit_match(StereoCoset::Lit(1), StereoCoset::Lit(1), StereoKind::Tetrahedral, true)]
    #[case::lit_miss(StereoCoset::Lit(0), StereoCoset::Lit(1), StereoKind::Tetrahedral, false)]
    #[case::set_member(StereoCoset::lit_set([1, 3]), StereoCoset::Lit(3), StereoKind::Octahedral, true)]
    #[case::set_nonmember(StereoCoset::lit_set([1, 3]), StereoCoset::Lit(2), StereoKind::Octahedral, false)]
    #[case::specific_vs_wildcard(StereoCoset::Lit(0), StereoCoset::Undetermined, StereoKind::Tetrahedral, false)]
    fn test_coset_matches(#[case] pattern: StereoCoset, #[case] target: StereoCoset, #[case] kind: StereoKind, #[case] expected: bool) {
        assert_eq!(coset_matches(&pattern, &target, kind), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCoset::Lit(0), Permutation::from_image(&[1, 0, 2, 3]), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::undetermined(StereoCoset::Undetermined, Permutation::from_image(&[1, 0, 2, 3]), StereoKind::Tetrahedral, StereoCoset::Undetermined)]
    fn test_coset_apply_permutation(#[case] coset: StereoCoset, #[case] permutation: Permutation, #[case] kind: StereoKind, #[case] expected: StereoCoset) {
        assert_eq!(coset_apply_permutation(&coset, permutation, kind), Some(expected));
    }

    #[rstest]
    fn test_stereo_atom_form_new() {
        let stereo_atom = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined);
        assert_eq!(
            stereo_atom.configuration,
            StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)
        );
        assert_eq!(stereo_atom.constraints, StereoAtomConstraintsForm::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined), false)]
    #[case::ground(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), true)]
    fn test_stereo_atom_form_is_ground(#[case] atom: StereoAtomForm, #[case] expected: bool) {
        assert_eq!(atom.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::open_coset(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined))]
    #[case::ground(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32))]
    fn test_stereo_atom_form_into_ground(#[case] atom: StereoAtomForm) {
        assert_eq!(atom.clone().into_ground(), atom);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absolute(
        StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() },
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    )]
    #[case::relative(
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() },
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    )]
    #[case::undetermined(
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() },
        StereoAtomForm::default(),
    )]
    #[case::explicit_open(
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Undetermined) }, ..Default::default() },
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
    )]
    #[case::constraint_set(
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), ..Default::default() },
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    )]
    #[case::constraint_remove(
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() },
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    )]
    fn test_stereo_atom_form_update(
        #[case] atom: StereoAtomForm,
        #[case] update: StereoAtomUpdate,
        #[case] expected: StereoAtomForm,
    ) {
        assert_eq!(atom.update(&update), expected);
    }

    #[rstest]
    #[case::empty(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32))]
    fn test_stereo_atom_form_update_identity(#[case] atom: StereoAtomForm) {
        assert_eq!(atom.update(&StereoAtomUpdate::default()), atom);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::configuration_and_constraint(
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoAtomForm::default(),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    )]
    #[case::constraint_context(
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    )]
    fn test_stereo_atom_form_difference_to(
        #[case] atom: StereoAtomForm,
        #[case] other: StereoAtomForm,
        #[case] expected: StereoAtomUpdate,
    ) {
        assert_eq!(atom.difference_to(&other), expected);
    }

    #[rstest]
    #[case::same(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32))]
    fn test_stereo_atom_form_difference_to_identity(#[case] atom: StereoAtomForm) {
        assert_eq!(atom.difference_to(&atom), StereoAtomUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::absolute(
        StereoBondForm::new(StereoKind::CisTrans, 0_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() },
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    )]
    #[case::relative(
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() },
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    )]
    #[case::undetermined(
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() },
        StereoBondForm::default(),
    )]
    #[case::explicit_open(
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Undetermined) }, ..Default::default() },
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
    )]
    #[case::constraint_set(
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))), ..Default::default() },
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
    )]
    #[case::constraint_remove(
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() },
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
    )]
    fn test_stereo_bond_form_update(
        #[case] bond: StereoBondForm,
        #[case] update: StereoBondUpdate,
        #[case] expected: StereoBondForm,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(StereoBondForm::new(StereoKind::CisTrans, 1_u32))]
    fn test_stereo_bond_form_update_identity(#[case] bond: StereoBondForm) {
        assert_eq!(bond.update(&StereoBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::configuration_and_constraint(
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoBondForm::default(),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    )]
    #[case::constraint_context(
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoBondForm::new(StereoKind::CisTrans, 1_u32),
        StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)) },
    )]
    fn test_stereo_bond_form_difference_to(
        #[case] bond: StereoBondForm,
        #[case] other: StereoBondForm,
        #[case] expected: StereoBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::same(StereoBondForm::new(StereoKind::CisTrans, 1_u32))]
    fn test_stereo_bond_form_difference_to_identity(#[case] bond: StereoBondForm) {
        assert_eq!(bond.difference_to(&bond), StereoBondUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        Some(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)))]
    #[case::different_kind(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Undetermined), None)]
    #[case::config_conflict(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), None)]
    fn test_stereo_atom_form_meet(
        #[case] a: StereoAtomForm,
        #[case] b: StereoAtomForm,
        #[case] expected: Option<StereoAtomForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_coset(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32))]
    #[case::distinct_cosets_widen(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::Tetrahedral, 2u32), StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([1, 2])))]
    fn test_stereo_atom_form_join(#[case] a: StereoAtomForm, #[case] b: StereoAtomForm, #[case] expected: StereoAtomForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_match(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), true)]
    #[case::different_kind(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32), StereoAtomForm::new(StereoKind::SquarePlanar, 1u32), false)]
    fn test_stereo_atom_form_matches(
        #[case] pattern: StereoAtomForm,
        #[case] target: StereoAtomForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_coset(
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([1])),
        Ok(StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)),
    )]
    #[case::empty_coset_litset_contradiction(
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set(Vec::<u32>::new())),
        Err(Contradiction),
    )]
    fn test_stereo_atom_form_canonicalize(
        #[case] input: StereoAtomForm,
        #[case] expected: Result<StereoAtomForm, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    fn test_stereo_bond_form_new() {
        let stereo_bond = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined);
        assert_eq!(
            stereo_bond.configuration,
            StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined)
        );
        assert_eq!(stereo_bond.constraints, StereoBondConstraintsForm::new())
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::same_kind_narrows(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined), StereoBondForm::new(StereoKind::CisTrans, 1u32),
        Some(StereoBondForm::new(StereoKind::CisTrans, 1u32)))]
    #[case::config_conflict(StereoBondForm::new(StereoKind::CisTrans, 0u32), StereoBondForm::new(StereoKind::CisTrans, 1u32), None)]
    fn test_stereo_bond_form_meet(
        #[case] a: StereoBondForm,
        #[case] b: StereoBondForm,
        #[case] expected: Option<StereoBondForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_coset(
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set([1])),
        Ok(StereoBondForm::new(StereoKind::CisTrans, 1u32)),
    )]
    #[case::empty_coset_litset_contradiction(
        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::lit_set(Vec::<u32>::new())),
        Err(Contradiction),
    )]
    fn test_stereo_bond_form_canonicalize(
        #[case] input: StereoBondForm,
        #[case] expected: Result<StereoBondForm, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rstest]
    #[case::identity(StereoKind::Tetrahedral, 0, Permutation::identity(4), 0)]
    #[case::involution(StereoKind::Tetrahedral, 0, StereoKind::Tetrahedral.involution(), 1)]
    #[case::involution_back(StereoKind::Tetrahedral, 1, StereoKind::Tetrahedral.involution(), 0)]
    fn test_stereo_kind_act(
        #[case] kind: StereoKind,
        #[case] index: u32,
        #[case] permutation: Permutation,
        #[case] expected: u32,
    ) {
        assert_eq!(kind.act(index, permutation), expected);
    }

    #[rstest]
    #[case::undetermined(
        StereoCoset::Undetermined,
        StereoKind::Tetrahedral,
        Permutation::identity(4),
        StereoCoset::Undetermined
    )]
    #[case::lit_identity(
        StereoCoset::Lit(0),
        StereoKind::Tetrahedral,
        Permutation::identity(4),
        StereoCoset::Lit(0)
    )]
    #[case::lit_involution(StereoCoset::Lit(0), StereoKind::Tetrahedral, StereoKind::Tetrahedral.involution(), StereoCoset::Lit(1))]
    #[case::lit_set(StereoCoset::lit_set([0]), StereoKind::Tetrahedral, StereoKind::Tetrahedral.involution(), StereoCoset::lit_set([1]))]
    #[case::term_layers(
        StereoCoset::term(StereoTerm::var("x")),
        StereoKind::Tetrahedral,
        Permutation::identity(4),
        StereoCoset::term(StereoTerm::apply(StereoTerm::var("x"), Permutation::identity(4)))
    )]
    fn test_stereo_coset_apply(
        #[case] coset: StereoCoset,
        #[case] kind: StereoKind,
        #[case] permutation: Permutation,
        #[case] expected: StereoCoset,
    ) {
        assert_eq!(coset.apply(kind, permutation), expected);
    }

    #[rstest]
    #[case::undetermined(
        StereoCoset::Undetermined,
        StereoKind::Tetrahedral,
        StereoCoset::Undetermined
    )]
    #[case::tetrahedral_0(StereoCoset::Lit(0), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::tetrahedral_1(StereoCoset::Lit(1), StereoKind::Tetrahedral, StereoCoset::Lit(0))]
    #[case::cis_trans(StereoCoset::Lit(0), StereoKind::CisTrans, StereoCoset::Lit(1))]
    #[case::term_layers(
        StereoCoset::term(StereoTerm::var("x")),
        StereoKind::Tetrahedral,
        StereoCoset::term(StereoTerm::swap(StereoTerm::var("x")))
    )]
    fn test_stereo_coset_swap(
        #[case] coset: StereoCoset,
        #[case] kind: StereoKind,
        #[case] expected: StereoCoset,
    ) {
        assert_eq!(coset.swap(kind), expected);
    }

    #[rstest]
    #[case::undetermined(
        StereoCoset::Undetermined,
        StereoKind::Tetrahedral,
        StereoCoset::Undetermined
    )]
    #[case::chiral(StereoCoset::Lit(0), StereoKind::Tetrahedral, StereoCoset::Lit(1))]
    #[case::achiral_noop(StereoCoset::Lit(0), StereoKind::CisTrans, StereoCoset::Lit(0))]
    #[case::term_layers(
        StereoCoset::term(StereoTerm::var("x")),
        StereoKind::Tetrahedral,
        StereoCoset::term(StereoTerm::mirror(StereoTerm::var("x")))
    )]
    fn test_stereo_coset_mirror(
        #[case] coset: StereoCoset,
        #[case] kind: StereoKind,
        #[case] expected: StereoCoset,
    ) {
        assert_eq!(coset.mirror(kind), expected);
    }

    #[rstest]
    #[case::undetermined(
        StereoConfigurationForm::Undetermined,
        Permutation::identity(4),
        StereoConfigurationForm::Undetermined
    )]
    #[case::kinded(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)), StereoKind::Tetrahedral.involution(), StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)))]
    fn test_stereo_configuration_form_apply(
        #[case] config: StereoConfigurationForm,
        #[case] permutation: Permutation,
        #[case] expected: StereoConfigurationForm,
    ) {
        assert_eq!(config.apply(permutation), expected);
    }

    #[rstest]
    #[case::undetermined(
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Undetermined
    )]
    #[case::kinded(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1))
    )]
    fn test_stereo_configuration_form_swap(
        #[case] config: StereoConfigurationForm,
        #[case] expected: StereoConfigurationForm,
    ) {
        assert_eq!(config.swap(), expected);
    }

    #[rstest]
    #[case::undetermined(
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Undetermined
    )]
    #[case::chiral(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1))
    )]
    #[case::achiral_noop(
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0))
    )]
    fn test_stereo_configuration_form_mirror(
        #[case] config: StereoConfigurationForm,
        #[case] expected: StereoConfigurationForm,
    ) {
        assert_eq!(config.mirror(), expected);
    }

    #[rstest]
    #[case::other_undetermined_keeps_self(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1))
    )]
    #[case::same_kind_undetermined_coset_keeps_coset(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1))
    )]
    #[case::same_kind_determined_coset_overrides(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0))
    )]
    #[case::different_kind_overrides(
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Undetermined),
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Undetermined)
    )]
    #[case::self_undetermined_takes_other(
        StereoConfigurationForm::Undetermined,
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)
    )]
    fn test_stereo_configuration_form_update(
        #[case] base: StereoConfigurationForm,
        #[case] other: StereoConfigurationForm,
        #[case] expected: StereoConfigurationForm,
    ) {
        assert_eq!(base.update(&other), expected);
    }

    #[rstest]
    #[case::apply(StereoAtomForm::new(StereoKind::Tetrahedral, 0u32), StereoKind::Tetrahedral.involution(), StereoAtomForm::new(StereoKind::Tetrahedral, 1u32))]
    fn test_stereo_atom_form_apply(
        #[case] input: StereoAtomForm,
        #[case] permutation: Permutation,
        #[case] expected: StereoAtomForm,
    ) {
        assert_eq!(input.apply(permutation), expected);
    }

    #[rstest]
    #[case::tetrahedral(
        StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)
    )]
    fn test_stereo_atom_form_swap(#[case] input: StereoAtomForm, #[case] expected: StereoAtomForm) {
        assert_eq!(input.swap(), expected);
    }

    #[rstest]
    #[case::chiral(
        StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)
    )]
    fn test_stereo_atom_form_mirror(
        #[case] input: StereoAtomForm,
        #[case] expected: StereoAtomForm,
    ) {
        assert_eq!(input.mirror(), expected);
    }

    #[rstest]
    #[case::cis_trans(
        StereoBondForm::new(StereoKind::CisTrans, 0u32),
        StereoBondForm::new(StereoKind::CisTrans, 1u32)
    )]
    fn test_stereo_bond_form_swap(#[case] input: StereoBondForm, #[case] expected: StereoBondForm) {
        assert_eq!(input.swap(), expected);
    }

    #[rstest]
    #[case::identity(
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        0,
    )]
    #[case::transposition(
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        [StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        1,
    )]
    #[case::even_cycle(
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        [StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        0,
    )]
    #[case::virtual_explicit_swap(
        [StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen)],
        [StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
        1,
    )]
    fn test_stereo_atom_form_transform_frame(
        #[case] before: [StereoLigand; 4],
        #[case] after: [StereoLigand; 4],
        #[case] expected_coset: u32,
    ) {
        let atom = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
        assert_eq!(
            atom.transform_frame(&before, &after),
            Some(StereoAtomForm::new(StereoKind::Tetrahedral, expected_coset,)),
        );
    }

    #[rstest]
    #[case::length(
        &[StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        &[StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom)],
    )]
    #[case::repetition(
        &[StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        &[StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
    )]
    #[case::membership(
        &[StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)],
        &[StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(5), StereoLigandKind::Atom)],
    )]
    fn test_stereo_atom_form_transform_frame_error(
        #[case] before: &[StereoLigand],
        #[case] after: &[StereoLigand],
    ) {
        assert_eq!(
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32).transform_frame(before, after),
            None,
        );
    }

    #[rstest]
    fn test_stereo_atom_form_transform_frame_self_inverse() {
        let before = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ];
        let after = [
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ];
        let atom = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
        assert_eq!(
            atom.transform_frame(&before, &after)
                .and_then(|transformed| transformed.transform_frame(&after, &before)),
            Some(atom),
        );
    }

    #[rstest]
    #[case::reordered(
        &[StereoLigand::new(AtomId(0), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen)],
        &[StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(0), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen)],
        Some(StereoBondForm::new(StereoKind::CisTrans, 1u32)),
    )]
    #[case::membership(
        &[StereoLigand::new(AtomId(0), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen)],
        &[StereoLigand::new(AtomId(0), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(1), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::LonePair)],
        None,
    )]
    fn test_stereo_bond_form_transform_frame(
        #[case] before: &[StereoLigand],
        #[case] after: &[StereoLigand],
        #[case] expected: Option<StereoBondForm>,
    ) {
        assert_eq!(
            StereoBondForm::new(StereoKind::CisTrans, 0u32).transform_frame(before, after),
            expected,
        );
    }
}
