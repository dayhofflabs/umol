//! Resolved edit vocabulary: the `Delta` counterpart of the deferred `Edit`.
//!
//! A `Delta` is one resolved edit over a `MoleculeAst`, referencing entities by stable
//! ids in the molecule's own id space (no positional `New`). The vocabulary is closed
//! under inversion — every delta's inverse is another delta.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::Hash;
use std::iter;
use std::mem::{discriminant, Discriminant};
use std::slice::{Iter, IterMut};

use umol_graph_core::{
    BiRelationData, FactorOrdering, ParticipantPosition, RelationData, Unordered,
};
use umol_perm::Permutation;

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::{
    AromaticSystemConstraintAst, AromaticSystemConstraintKey, AtomConstraintAst, AtomConstraintKey,
    BondConstraintAst, BondConstraintKey, Constraint, DativeBondConstraintAst,
    DativeBondConstraintKey, MulticenterBondConstraintAst, MulticenterBondConstraintKey,
    NoncovalentBondConstraintAst, NoncovalentBondConstraintKey, StereoAtomConstraintAst,
    StereoAtomConstraintKey, StereoBondConstraintAst, StereoBondConstraintKey,
};
use super::dative::DativeBondAst;
use super::edit::{
    AromaticSystemFieldChange, AtomFieldChange, BondFieldChange, DativeBondFieldChange,
    MulticenterBondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange,
    StereoBondFieldChange,
};
use super::error::Contradiction;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::remap::IdRemapping;
use super::stereo::{CosetOp, StereoAtomAst, StereoBondAst, StereoConfigurationAst, StereoKind};
use super::traits::{Canonicalize, EntityPatch};

/// A resolved edit to a single atom.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum AtomDelta {
    Add {
        id: AtomId,
        ast: AtomAst,
    },
    Remove {
        id: AtomId,
        ast: AtomAst,
    },
    ModifyField {
        id: AtomId,
        change: AtomFieldChange,
    },
    ModifyConstraint {
        id: AtomId,
        old: Option<AtomConstraintAst>,
        new: Option<AtomConstraintAst>,
    },
}

impl AtomDelta {
    /// The inverse delta: `Add`↔`Remove`; `ModifyField` / `ModifyConstraint` swap old/new.
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, ast } => Self::Remove { id, ast },
            Self::Remove { id, ast } => Self::Add { id, ast },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: new,
                new: old,
            },
        }
    }
}

/// A resolved edit to a single bond.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum BondDelta {
    Add {
        id: BondId,
        atoms: [AtomId; 2],
        ast: BondAst,
    },
    Remove {
        id: BondId,
        atoms: [AtomId; 2],
        ast: BondAst,
    },
    ModifyField {
        id: BondId,
        change: BondFieldChange,
    },
    ModifyConstraint {
        id: BondId,
        old: Option<BondConstraintAst>,
        new: Option<BondConstraintAst>,
    },
}

impl BondDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, atoms, ast } => Self::Remove { id, atoms, ast },
            Self::Remove { id, atoms, ast } => Self::Add { id, atoms, ast },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: new,
                new: old,
            },
        }
    }
}

/// A resolved edit to a single dative bond. `donors`/`acceptor` are the directed
/// participants (structural payload, like `BondDelta::atoms`); identity is the id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum DativeBondDelta {
    Add {
        id: DativeBondId,
        donors: Vec<AtomId>,
        acceptor: AtomId,
        ast: DativeBondAst,
    },
    Remove {
        id: DativeBondId,
        donors: Vec<AtomId>,
        acceptor: AtomId,
        ast: DativeBondAst,
    },
    ModifyField {
        id: DativeBondId,
        change: DativeBondFieldChange,
    },
    ModifyConstraint {
        id: DativeBondId,
        old: Option<DativeBondConstraintAst>,
        new: Option<DativeBondConstraintAst>,
    },
}

impl DativeBondDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                donors,
                acceptor,
                ast,
            } => Self::Remove {
                id,
                donors,
                acceptor,
                ast,
            },
            Self::Remove {
                id,
                donors,
                acceptor,
                ast,
            } => Self::Add {
                id,
                donors,
                acceptor,
                ast,
            },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: new,
                new: old,
            },
        }
    }
}

/// A resolved edit to a single aromatic system. `atoms` are the member atoms
/// (structural payload); identity is the id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum AromaticSystemDelta {
    Add {
        id: AromaticSystemId,
        atoms: Vec<AtomId>,
        ast: AromaticSystemAst,
    },
    Remove {
        id: AromaticSystemId,
        atoms: Vec<AtomId>,
        ast: AromaticSystemAst,
    },
    ModifyField {
        id: AromaticSystemId,
        change: AromaticSystemFieldChange,
    },
    ModifyConstraint {
        id: AromaticSystemId,
        old: Option<AromaticSystemConstraintAst>,
        new: Option<AromaticSystemConstraintAst>,
    },
}

impl AromaticSystemDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, atoms, ast } => Self::Remove { id, atoms, ast },
            Self::Remove { id, atoms, ast } => Self::Add { id, atoms, ast },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: new,
                new: old,
            },
        }
    }
}

/// A resolved edit to a single multicenter bond. `atoms` are the member atoms
/// (structural payload); identity is the id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum MulticenterBondDelta {
    Add {
        id: MulticenterBondId,
        atoms: Vec<AtomId>,
        ast: MulticenterBondAst,
    },
    Remove {
        id: MulticenterBondId,
        atoms: Vec<AtomId>,
        ast: MulticenterBondAst,
    },
    ModifyField {
        id: MulticenterBondId,
        change: MulticenterBondFieldChange,
    },
    ModifyConstraint {
        id: MulticenterBondId,
        old: Option<MulticenterBondConstraintAst>,
        new: Option<MulticenterBondConstraintAst>,
    },
}

impl MulticenterBondDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, atoms, ast } => Self::Remove { id, atoms, ast },
            Self::Remove { id, atoms, ast } => Self::Add { id, atoms, ast },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: new,
                new: old,
            },
        }
    }
}

/// A resolved edit to a single noncovalent bond. `atoms` are its two participants
/// (structural payload); identity is the id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum NoncovalentBondDelta {
    Add {
        id: NoncovalentBondId,
        atoms: [AtomId; 2],
        ast: NoncovalentBondAst,
    },
    Remove {
        id: NoncovalentBondId,
        atoms: [AtomId; 2],
        ast: NoncovalentBondAst,
    },
    ModifyField {
        id: NoncovalentBondId,
        change: NoncovalentBondFieldChange,
    },
    ModifyConstraint {
        id: NoncovalentBondId,
        old: Option<NoncovalentBondConstraintAst>,
        new: Option<NoncovalentBondConstraintAst>,
    },
}

impl NoncovalentBondDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, atoms, ast } => Self::Remove { id, atoms, ast },
            Self::Remove { id, atoms, ast } => Self::Add { id, atoms, ast },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: new,
                new: old,
            },
        }
    }
}

/// A resolved edit to a single stereo atom. `site` + `ligands` are the structural payload; identity
/// is the id. `ModifyField`/`ModifyConstraint` are the absolute set-ops (as for DAMN); `Apply`,
/// `Swap`, `Mirror` are the frame-relative coset ops — they carry no pre-state, resolving against
/// the matched host configuration at apply.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum StereoAtomDelta {
    Add {
        id: StereoAtomId,
        site: AtomId,
        ligands: Vec<StereoLigand>,
        ast: StereoAtomAst,
    },
    Remove {
        id: StereoAtomId,
        site: AtomId,
        ligands: Vec<StereoLigand>,
        ast: StereoAtomAst,
    },
    ModifyField {
        id: StereoAtomId,
        change: StereoAtomFieldChange,
    },
    ModifyConstraint {
        id: StereoAtomId,
        /// Serialization context: the geometry kind the constraint renders/parses against (its
        /// permutation degree, `~` shortcut). `None` for a kind-free constraint on an
        /// `Undetermined`-geometry center. Not read by apply/canonicalize/diff.
        kind: Option<StereoKind>,
        old: Option<StereoAtomConstraintAst>,
        new: Option<StereoAtomConstraintAst>,
    },
    Apply {
        id: StereoAtomId,
        kind: StereoKind,
        permutation: Permutation,
    },
    Swap {
        id: StereoAtomId,
        kind: StereoKind,
    },
    Mirror {
        id: StereoAtomId,
        kind: StereoKind,
    },
}

impl StereoAtomDelta {
    pub fn id(&self) -> StereoAtomId {
        match self {
            Self::Add { id, .. }
            | Self::Remove { id, .. }
            | Self::ModifyField { id, .. }
            | Self::ModifyConstraint { id, .. }
            | Self::Apply { id, .. }
            | Self::Swap { id, .. }
            | Self::Mirror { id, .. } => *id,
        }
    }

    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                site,
                ligands,
                ast,
            } => Self::Remove {
                id,
                site,
                ligands,
                ast,
            },
            Self::Remove {
                id,
                site,
                ligands,
                ast,
            } => Self::Add {
                id,
                site,
                ligands,
                ast,
            },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, kind, old, new } => Self::ModifyConstraint {
                id,
                kind,
                old: new,
                new: old,
            },
            Self::Apply {
                id,
                kind,
                permutation,
            } => Self::Apply {
                id,
                kind,
                permutation: permutation.inverse(),
            },
            Self::Swap { id, kind } => Self::Swap { id, kind },
            Self::Mirror { id, kind } => Self::Mirror { id, kind },
        }
    }
}

/// A resolved edit to a single stereo bond. `site` (a bond) + `ligands` are the structural payload;
/// identity is the id. Same op vocabulary as `StereoAtomDelta`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum StereoBondDelta {
    Add {
        id: StereoBondId,
        site: BondId,
        ligands: Vec<StereoLigand>,
        ast: StereoBondAst,
    },
    Remove {
        id: StereoBondId,
        site: BondId,
        ligands: Vec<StereoLigand>,
        ast: StereoBondAst,
    },
    ModifyField {
        id: StereoBondId,
        change: StereoBondFieldChange,
    },
    ModifyConstraint {
        id: StereoBondId,
        /// Serialization context — see `StereoAtomDelta::ModifyConstraint`.
        kind: Option<StereoKind>,
        old: Option<StereoBondConstraintAst>,
        new: Option<StereoBondConstraintAst>,
    },
    Apply {
        id: StereoBondId,
        kind: StereoKind,
        permutation: Permutation,
    },
    Swap {
        id: StereoBondId,
        kind: StereoKind,
    },
    Mirror {
        id: StereoBondId,
        kind: StereoKind,
    },
}

impl StereoBondDelta {
    pub fn id(&self) -> StereoBondId {
        match self {
            Self::Add { id, .. }
            | Self::Remove { id, .. }
            | Self::ModifyField { id, .. }
            | Self::ModifyConstraint { id, .. }
            | Self::Apply { id, .. }
            | Self::Swap { id, .. }
            | Self::Mirror { id, .. } => *id,
        }
    }

    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                site,
                ligands,
                ast,
            } => Self::Remove {
                id,
                site,
                ligands,
                ast,
            },
            Self::Remove {
                id,
                site,
                ligands,
                ast,
            } => Self::Add {
                id,
                site,
                ligands,
                ast,
            },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, kind, old, new } => Self::ModifyConstraint {
                id,
                kind,
                old: new,
                new: old,
            },
            Self::Apply {
                id,
                kind,
                permutation,
            } => Self::Apply {
                id,
                kind,
                permutation: permutation.inverse(),
            },
            Self::Swap { id, kind } => Self::Swap { id, kind },
            Self::Mirror { id, kind } => Self::Mirror { id, kind },
        }
    }
}

/// A resolved change to the molecule-level constraint set, as a set-diff.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ConstraintDelta {
    Add(Constraint),
    Remove(Constraint),
}

impl ConstraintDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add(constraint) => Self::Remove(constraint),
            Self::Remove(constraint) => Self::Add(constraint),
        }
    }
}

/// One resolved edit across the localized-topology families.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Delta {
    Atom(AtomDelta),
    Bond(BondDelta),
    DativeBond(DativeBondDelta),
    AromaticSystem(AromaticSystemDelta),
    MulticenterBond(MulticenterBondDelta),
    NoncovalentBond(NoncovalentBondDelta),
    StereoAtom(StereoAtomDelta),
    StereoBond(StereoBondDelta),
    Constraint(ConstraintDelta),
}

impl Delta {
    /// The inverse delta.
    pub fn inverse(self) -> Self {
        match self {
            Self::Atom(delta) => Self::Atom(delta.inverse()),
            Self::Bond(delta) => Self::Bond(delta.inverse()),
            Self::DativeBond(delta) => Self::DativeBond(delta.inverse()),
            Self::AromaticSystem(delta) => Self::AromaticSystem(delta.inverse()),
            Self::MulticenterBond(delta) => Self::MulticenterBond(delta.inverse()),
            Self::NoncovalentBond(delta) => Self::NoncovalentBond(delta.inverse()),
            Self::StereoAtom(delta) => Self::StereoAtom(delta.inverse()),
            Self::StereoBond(delta) => Self::StereoBond(delta.inverse()),
            Self::Constraint(delta) => Self::Constraint(delta.inverse()),
        }
    }
}

/// Per-variant diff/apply ops for the `EntityPatch` impl, from the `(variant => ast field)` map:
/// `apply_field`, `diff_field`, `diff_constraints`.
macro_rules! diff_field_ops {
    ($change:ident, $ast:ident, $constraint:ident, { $($variant:ident => $field:ident),+ $(,)? }) => {
        fn apply_field(ast: &mut $ast, change: $change) -> Result<(), Contradiction> {
            match change {
                $(
                    $change::$variant { old, new } => {
                        if !ast.$field.canonical_eq(&old) {
                            return Err(Contradiction);
                        }
                        ast.$field = new;
                    }
                )+
            }
            Ok(())
        }

        fn diff_field(lhs: &$ast, rhs: &$ast) -> Vec<$change> {
            let mut out = Vec::new();
            $(
                if !lhs.$field.canonical_eq(&rhs.$field) {
                    out.push($change::$variant {
                        old: lhs.$field.clone(),
                        new: rhs.$field.clone(),
                    });
                }
            )+
            out
        }

        #[allow(clippy::type_complexity)]
        fn diff_constraints(
            lhs: &$ast,
            rhs: &$ast,
        ) -> Vec<(Option<$constraint>, Option<$constraint>)> {
            let mut lhs_by_key: HashMap<_, $constraint> = HashMap::new();
            for constraint in lhs.constraints.iter() {
                lhs_by_key.insert(constraint.key(), constraint.clone());
            }
            let mut rhs_by_key: HashMap<_, $constraint> = HashMap::new();
            for constraint in rhs.constraints.iter() {
                rhs_by_key.insert(constraint.key(), constraint.clone());
            }
            let mut keys: HashSet<_> = lhs_by_key.keys().cloned().collect();
            keys.extend(rhs_by_key.keys().cloned());
            let mut out = Vec::new();
            for key in keys {
                let l = lhs_by_key.get(&key).cloned();
                let r = rhs_by_key.get(&key).cloned();
                if !options_canonical_eq(&l, &r) {
                    out.push((l, r));
                }
            }
            out
        }
    };
}

/// Canonical equivalence over optional payloads: both absent is equal, both present compares by
/// `canonical_eq`, presence mismatch is unequal.
fn options_canonical_eq<T: Canonicalize>(l: &Option<T>, r: &Option<T>) -> bool {
    match (l, r) {
        (None, None) => true,
        (Some(a), Some(b)) => a.canonical_eq(b),
        _ => false,
    }
}

/// Per-variant fold ops for the crate-private `EntityFold` impl: `fuse_field`,
/// `field_is_identity`.
macro_rules! fold_field_ops {
    ($change:ident, { $($variant:ident),+ $(,)? }) => {
        fn fuse_field(prev: $change, next: $change) -> Option<$change> {
            match (prev, next) {
                $(
                    (
                        $change::$variant { old, new: prev_new },
                        $change::$variant { old: next_old, new },
                    ) if prev_new.canonical_eq(&next_old) => Some($change::$variant { old, new }),
                )+
                #[allow(unreachable_patterns)]
                _ => None,
            }
        }

        fn field_is_identity(change: &$change) -> bool {
            match change {
                $( $change::$variant { old, new } => old.canonical_eq(new), )+
            }
        }
    };
}

/// One entity's span across a reaction — its slice of the superimposed `L`∪`K`∪`R`. A *state*, not
/// an operation (unlike `Edit` / `Delta`). `lhs()` / `rhs()` read the side values.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EntitySpan<T> {
    /// In the interface `K` — present and identical on both sides.
    Unchanged(T),
    /// In the interface `K` — present on both sides but relabeled (a dynamic entity).
    Modified { lhs: T, rhs: T },
    /// In `R` only — created.
    Added(T),
    /// In `L` only — deleted.
    Removed(T),
}

impl<T> EntitySpan<T> {
    /// The lhs (`L`) value, or `None` if the entity is created.
    pub fn lhs(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Removed(value) | Self::Modified { lhs: value, .. } => {
                Some(value)
            }
            Self::Added(_) => None,
        }
    }

    /// The rhs (`R`) value, or `None` if the entity is deleted.
    pub fn rhs(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Added(value) | Self::Modified { rhs: value, .. } => {
                Some(value)
            }
            Self::Removed(_) => None,
        }
    }
}

/// A span stored as relation data (a `ReactionSpanAst` overlay) reindexes each present side's payload
/// and compares side-wise, delegating to the underlying payload's [`RelationData`].
impl<U: RelationData> RelationData for EntitySpan<U> {
    fn on_permutation(&mut self, order: &[ParticipantPosition]) {
        match self {
            Self::Unchanged(value) | Self::Added(value) | Self::Removed(value) => {
                value.on_permutation(order)
            }
            Self::Modified { lhs, rhs } => {
                lhs.on_permutation(order);
                rhs.on_permutation(order);
            }
        }
    }

    fn is_permutation_invariant(&self) -> bool {
        self.lhs().is_none_or(U::is_permutation_invariant)
            && self.rhs().is_none_or(U::is_permutation_invariant)
    }
}

impl<U: BiRelationData> BiRelationData for EntitySpan<U> {
    fn on_permutation(&mut self, order_1: &[ParticipantPosition], order_2: &[ParticipantPosition]) {
        match self {
            Self::Unchanged(value) | Self::Added(value) | Self::Removed(value) => {
                value.on_permutation(order_1, order_2)
            }
            Self::Modified { lhs, rhs } => {
                lhs.on_permutation(order_1, order_2);
                rhs.on_permutation(order_1, order_2);
            }
        }
    }

    fn is_permutation_invariant(&self) -> bool {
        self.lhs().is_none_or(U::is_permutation_invariant)
            && self.rhs().is_none_or(U::is_permutation_invariant)
    }
}

impl<T: Canonicalize> EntitySpan<T> {
    /// Superimpose an entity's optional lhs and rhs values into a span — the per-entity kernel of
    /// `ReactionSpanAst::superimpose`: present-both maps to `Unchanged` (equal) or `Modified`,
    /// lhs-only to `Removed`, rhs-only to `Added`, neither to `None`.
    pub fn superimpose(lhs: Option<T>, rhs: Option<T>) -> Option<Self> {
        match (lhs, rhs) {
            (Some(lhs), Some(rhs)) if lhs.canonical_eq(&rhs) => Some(Self::Unchanged(lhs)),
            (Some(lhs), Some(rhs)) => Some(Self::Modified { lhs, rhs }),
            (Some(lhs), None) => Some(Self::Removed(lhs)),
            (None, Some(rhs)) => Some(Self::Added(rhs)),
            (None, None) => None,
        }
    }
}

/// A molecule-level constraint's span across a reaction — its slice of the superimposed `L`∪`K`∪`R`.
/// A *state*, not an operation (unlike `ConstraintDelta`). `lhs()` / `rhs()` read the side values.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum ConstraintSpan {
    /// In the interface `K` — present and identical on both sides.
    Unchanged(Constraint),
    /// In `R` only — created.
    Added(Constraint),
    /// In `L` only — deleted.
    Removed(Constraint),
}

impl ConstraintSpan {
    /// The lhs (`L`) value, or `None` if the constraint is created.
    pub fn lhs(&self) -> Option<&Constraint> {
        match self {
            Self::Unchanged(value) | Self::Removed(value) => Some(value),
            Self::Added(_) => None,
        }
    }

    /// The rhs (`R`) value, or `None` if the constraint is deleted.
    pub fn rhs(&self) -> Option<&Constraint> {
        match self {
            Self::Unchanged(value) | Self::Added(value) => Some(value),
            Self::Removed(_) => None,
        }
    }
}

/// The per-entity op the fold operates on, abstracting `AtomDelta`/`BondDelta`. `atoms` carries
/// the entity's participant atoms in `Add`/`Remove` (`()` for an atom, its two ids for a bond).
pub(crate) enum EntityOp<F: EntityFold> {
    Add {
        atoms: F::Atoms,
        ast: F::Ast,
    },
    Remove {
        atoms: F::Atoms,
        ast: F::Ast,
    },
    ModifyField(F::FieldChange),
    ModifyConstraint {
        old: Option<F::Constraint>,
        new: Option<F::Constraint>,
    },
}

/// The canonicalize-fold extension of `EntityPatch` — the `EntityOp` deconstruction
/// (`split`/`rebuild`), per-field/constraint fold helpers, and span→deltas recovery.
/// Crate-internal: only the fold and span lowering use it.
pub(crate) trait EntityFold: EntityPatch {
    type ConstraintKey: Clone + Eq + Hash;
    type Atoms;

    fn id(&self) -> Self::Id;
    fn split(self) -> EntityOp<Self>;
    fn rebuild(id: Self::Id, op: EntityOp<Self>) -> Self;
    fn into_delta(self) -> Delta;

    fn fuse_field(prev: Self::FieldChange, next: Self::FieldChange) -> Option<Self::FieldChange>;
    fn field_is_identity(change: &Self::FieldChange) -> bool;
    fn field_inverse(change: Self::FieldChange) -> Self::FieldChange;
    fn constraint_key(constraint: &Self::Constraint) -> Self::ConstraintKey;

    /// Recover this kind's deltas from its lhs/rhs state column: `Added`/`Removed` become
    /// structural `Add`/`Remove` (their `atoms` from `atoms(index)`), `Modified` becomes the
    /// field/constraint `diff`. The id of entity `i` is `i` (the column is id-indexed).
    fn deltas_from_states(
        states: &[EntitySpan<Self::Ast>],
        atoms: impl Fn(usize) -> Self::Atoms,
    ) -> Vec<Delta> {
        let mut out = Vec::new();
        for (index, state) in states.iter().enumerate() {
            let id = Self::Id::from(index);
            match state {
                EntitySpan::Unchanged(_) => {}
                EntitySpan::Added(ast) => out.push(
                    Self::rebuild(
                        id,
                        EntityOp::Add {
                            atoms: atoms(index),
                            ast: ast.clone(),
                        },
                    )
                    .into_delta(),
                ),
                EntitySpan::Removed(ast) => out.push(
                    Self::rebuild(
                        id,
                        EntityOp::Remove {
                            atoms: atoms(index),
                            ast: ast.clone(),
                        },
                    )
                    .into_delta(),
                ),
                EntitySpan::Modified { lhs, rhs } => {
                    out.extend(Self::diff(id, lhs, rhs).into_iter().map(Self::into_delta));
                }
            }
        }
        out
    }
}

/// Fold one entity's ops (input order) to its normal form, branching on the `created`
/// (an `Add` is present) vs `preserved` (no `Add`) path.
fn fold_group<F: EntityFold>(id: F::Id, group: Vec<F>) -> Result<Vec<F>, Contradiction> {
    let ops: Vec<EntityOp<F>> = group.into_iter().map(F::split).collect();
    let created = ops.iter().any(|op| matches!(op, EntityOp::Add { .. }));
    let folded = if created {
        fold_created(ops)?
    } else {
        fold_preserved(ops)?
    };
    Ok(folded.into_iter().map(|op| F::rebuild(id, op)).collect())
}

/// Created entity: seed `ast` from `Add`, absorb subsequent field/constraint changes; an
/// `Add`+`Remove` cancels. Yields one `Add` with the final ast, or nothing.
fn fold_created<F: EntityFold>(ops: Vec<EntityOp<F>>) -> Result<Vec<EntityOp<F>>, Contradiction> {
    let mut state: Option<(F::Atoms, F::Ast)> = None;
    let mut removed = false;
    for op in ops {
        if removed {
            return Err(Contradiction);
        }
        match op {
            EntityOp::Add { atoms, ast } => {
                if state.is_some() {
                    return Err(Contradiction);
                }
                state = Some((atoms, ast));
            }
            EntityOp::ModifyField(change) => {
                let (_, ast) = state.as_mut().ok_or(Contradiction)?;
                F::apply_field(ast, change)?;
            }
            EntityOp::ModifyConstraint { old, new } => {
                let (_, ast) = state.as_mut().ok_or(Contradiction)?;
                F::apply_constraint(ast, old, new)?;
            }
            EntityOp::Remove { .. } => {
                if state.is_none() {
                    return Err(Contradiction);
                }
                state = None;
                removed = true;
            }
        }
    }
    Ok(match state {
        Some((atoms, ast)) => vec![EntityOp::Add { atoms, ast }],
        None => Vec::new(),
    })
}

/// Preserved entity: fuse `ModifyField` chains per field and `ModifyConstraint` chains per key. A
/// `Remove` subsumes the prior changes and carries the *original* value (the changes are
/// reverted on the removed ast).
#[allow(clippy::type_complexity)]
fn fold_preserved<F: EntityFold>(ops: Vec<EntityOp<F>>) -> Result<Vec<EntityOp<F>>, Contradiction> {
    let mut fields: HashMap<Discriminant<F::FieldChange>, F::FieldChange> = HashMap::new();
    let mut constraints: HashMap<F::ConstraintKey, (Option<F::Constraint>, Option<F::Constraint>)> =
        HashMap::new();
    let mut removed: Option<(F::Atoms, F::Ast)> = None;
    for op in ops {
        if removed.is_some() {
            return Err(Contradiction);
        }
        match op {
            EntityOp::Add { .. } => return Err(Contradiction),
            EntityOp::ModifyField(change) => {
                let slot = discriminant(&change);
                let fused = match fields.remove(&slot) {
                    Some(prev) => F::fuse_field(prev, change).ok_or(Contradiction)?,
                    None => change,
                };
                fields.insert(slot, fused);
            }
            EntityOp::ModifyConstraint { old, new } => {
                let key = match old.as_ref().or(new.as_ref()) {
                    Some(constraint) => F::constraint_key(constraint),
                    None => continue,
                };
                match constraints.remove(&key) {
                    Some((first_old, prev_new)) => {
                        if !options_canonical_eq(&prev_new, &old) {
                            return Err(Contradiction);
                        }
                        constraints.insert(key, (first_old, new));
                    }
                    None => {
                        constraints.insert(key, (old, new));
                    }
                }
            }
            EntityOp::Remove { atoms, ast } => {
                removed = Some((atoms, ast));
            }
        }
    }
    if let Some((atoms, mut ast)) = removed {
        for (_slot, change) in fields {
            F::apply_field(&mut ast, F::field_inverse(change))?;
        }
        for (_key, (old, new)) in constraints {
            F::apply_constraint(&mut ast, new, old)?;
        }
        return Ok(vec![EntityOp::Remove { atoms, ast }]);
    }
    let mut out = Vec::new();
    for (_slot, change) in fields {
        if !F::field_is_identity(&change) {
            out.push(EntityOp::ModifyField(change));
        }
    }
    for (_key, (old, new)) in constraints {
        if !options_canonical_eq(&old, &new) {
            out.push(EntityOp::ModifyConstraint { old, new });
        }
    }
    Ok(out)
}

impl EntityPatch for AtomDelta {
    type Id = AtomId;
    type Ast = AtomAst;
    type FieldChange = AtomFieldChange;
    type Constraint = AtomConstraintAst;

    fn modify_field(id: AtomId, change: AtomFieldChange) -> Self {
        AtomDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: AtomId,
        old: Option<AtomConstraintAst>,
        new: Option<AtomConstraintAst>,
    ) -> Self {
        AtomDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(AtomFieldChange, AtomAst, AtomConstraintAst, {
        Element => element,
        IsotopeMass => isotope_mass,
        Charge => charge,
        ImplicitHydrogens => implicit_hydrogens,
        LonePairs => lone_pairs,
        Spin => spin,
    });

    fn apply_constraint(
        ast: &mut AtomAst,
        old: Option<AtomConstraintAst>,
        new: Option<AtomConstraintAst>,
    ) -> Result<(), Contradiction> {
        ast.constraints.compare_and_set(old, new)
    }
}

impl EntityFold for AtomDelta {
    type ConstraintKey = AtomConstraintKey;
    type Atoms = ();

    fn id(&self) -> AtomId {
        match self {
            AtomDelta::Add { id, .. }
            | AtomDelta::Remove { id, .. }
            | AtomDelta::ModifyField { id, .. }
            | AtomDelta::ModifyConstraint { id, .. } => *id,
        }
    }

    fn split(self) -> EntityOp<Self> {
        match self {
            AtomDelta::Add { ast, .. } => EntityOp::Add { atoms: (), ast },
            AtomDelta::Remove { ast, .. } => EntityOp::Remove { atoms: (), ast },
            AtomDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            AtomDelta::ModifyConstraint { old, new, .. } => EntityOp::ModifyConstraint { old, new },
        }
    }

    fn rebuild(id: AtomId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { ast, .. } => AtomDelta::Add { id, ast },
            EntityOp::Remove { ast, .. } => AtomDelta::Remove { id, ast },
            EntityOp::ModifyField(change) => Self::modify_field(id, change),
            EntityOp::ModifyConstraint { old, new } => Self::modify_constraint(id, old, new),
        }
    }

    fn into_delta(self) -> Delta {
        Delta::Atom(self)
    }

    fn field_inverse(change: AtomFieldChange) -> AtomFieldChange {
        change.inverse()
    }

    fn constraint_key(constraint: &AtomConstraintAst) -> AtomConstraintKey {
        constraint.key()
    }

    fold_field_ops!(AtomFieldChange, {
        Element,
        IsotopeMass,
        Charge,
        ImplicitHydrogens,
        LonePairs,
        Spin,
    });
}

impl EntityPatch for BondDelta {
    type Id = BondId;
    type Ast = BondAst;
    type FieldChange = BondFieldChange;
    type Constraint = BondConstraintAst;

    fn modify_field(id: BondId, change: BondFieldChange) -> Self {
        BondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: BondId,
        old: Option<BondConstraintAst>,
        new: Option<BondConstraintAst>,
    ) -> Self {
        BondDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(BondFieldChange, BondAst, BondConstraintAst, {
        Order => order,
        Charge => charge,
        Spin => spin,
    });

    fn apply_constraint(
        ast: &mut BondAst,
        old: Option<BondConstraintAst>,
        new: Option<BondConstraintAst>,
    ) -> Result<(), Contradiction> {
        ast.constraints.compare_and_set(old, new)
    }
}

impl EntityFold for BondDelta {
    type ConstraintKey = BondConstraintKey;
    type Atoms = [AtomId; 2];

    fn id(&self) -> BondId {
        match self {
            BondDelta::Add { id, .. }
            | BondDelta::Remove { id, .. }
            | BondDelta::ModifyField { id, .. }
            | BondDelta::ModifyConstraint { id, .. } => *id,
        }
    }

    fn split(self) -> EntityOp<Self> {
        match self {
            BondDelta::Add { atoms, ast, .. } => EntityOp::Add { atoms, ast },
            BondDelta::Remove { atoms, ast, .. } => EntityOp::Remove { atoms, ast },
            BondDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            BondDelta::ModifyConstraint { old, new, .. } => EntityOp::ModifyConstraint { old, new },
        }
    }

    fn rebuild(id: BondId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { atoms, ast } => BondDelta::Add { id, atoms, ast },
            EntityOp::Remove { atoms, ast } => BondDelta::Remove { id, atoms, ast },
            EntityOp::ModifyField(change) => Self::modify_field(id, change),
            EntityOp::ModifyConstraint { old, new } => Self::modify_constraint(id, old, new),
        }
    }

    fn into_delta(self) -> Delta {
        Delta::Bond(self)
    }

    fn field_inverse(change: BondFieldChange) -> BondFieldChange {
        change.inverse()
    }

    fn constraint_key(constraint: &BondConstraintAst) -> BondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(BondFieldChange, {
        Order,
        Charge,
        Spin,
    });
}

impl EntityPatch for DativeBondDelta {
    type Id = DativeBondId;
    type Ast = DativeBondAst;
    type FieldChange = DativeBondFieldChange;
    type Constraint = DativeBondConstraintAst;

    fn modify_field(id: DativeBondId, change: DativeBondFieldChange) -> Self {
        DativeBondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: DativeBondId,
        old: Option<DativeBondConstraintAst>,
        new: Option<DativeBondConstraintAst>,
    ) -> Self {
        DativeBondDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(DativeBondFieldChange, DativeBondAst, DativeBondConstraintAst, {
        Order => order,
    });

    fn apply_constraint(
        ast: &mut DativeBondAst,
        old: Option<DativeBondConstraintAst>,
        new: Option<DativeBondConstraintAst>,
    ) -> Result<(), Contradiction> {
        ast.constraints.compare_and_set(old, new)
    }
}

impl EntityFold for DativeBondDelta {
    type ConstraintKey = DativeBondConstraintKey;
    type Atoms = (Vec<AtomId>, AtomId);

    fn id(&self) -> DativeBondId {
        match self {
            DativeBondDelta::Add { id, .. }
            | DativeBondDelta::Remove { id, .. }
            | DativeBondDelta::ModifyField { id, .. }
            | DativeBondDelta::ModifyConstraint { id, .. } => *id,
        }
    }

    fn split(self) -> EntityOp<Self> {
        match self {
            DativeBondDelta::Add {
                donors,
                acceptor,
                ast,
                ..
            } => EntityOp::Add {
                atoms: (donors, acceptor),
                ast,
            },
            DativeBondDelta::Remove {
                donors,
                acceptor,
                ast,
                ..
            } => EntityOp::Remove {
                atoms: (donors, acceptor),
                ast,
            },
            DativeBondDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            DativeBondDelta::ModifyConstraint { old, new, .. } => {
                EntityOp::ModifyConstraint { old, new }
            }
        }
    }

    fn rebuild(id: DativeBondId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add {
                atoms: (donors, acceptor),
                ast,
            } => DativeBondDelta::Add {
                id,
                donors,
                acceptor,
                ast,
            },
            EntityOp::Remove {
                atoms: (donors, acceptor),
                ast,
            } => DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                ast,
            },
            EntityOp::ModifyField(change) => Self::modify_field(id, change),
            EntityOp::ModifyConstraint { old, new } => Self::modify_constraint(id, old, new),
        }
    }

    fn into_delta(self) -> Delta {
        Delta::DativeBond(self)
    }

    fn field_inverse(change: DativeBondFieldChange) -> DativeBondFieldChange {
        change.inverse()
    }

    fn constraint_key(constraint: &DativeBondConstraintAst) -> DativeBondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(DativeBondFieldChange, { Order });
}

impl EntityPatch for AromaticSystemDelta {
    type Id = AromaticSystemId;
    type Ast = AromaticSystemAst;
    type FieldChange = AromaticSystemFieldChange;
    type Constraint = AromaticSystemConstraintAst;

    fn modify_field(id: AromaticSystemId, change: AromaticSystemFieldChange) -> Self {
        AromaticSystemDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: AromaticSystemId,
        old: Option<AromaticSystemConstraintAst>,
        new: Option<AromaticSystemConstraintAst>,
    ) -> Self {
        AromaticSystemDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(
        AromaticSystemFieldChange,
        AromaticSystemAst,
        AromaticSystemConstraintAst,
        {
            Electrons => electrons,
            Charge => charge,
            Spin => spin,
        }
    );

    fn apply_constraint(
        ast: &mut AromaticSystemAst,
        old: Option<AromaticSystemConstraintAst>,
        new: Option<AromaticSystemConstraintAst>,
    ) -> Result<(), Contradiction> {
        ast.constraints.compare_and_set(old, new)
    }
}

impl EntityFold for AromaticSystemDelta {
    type ConstraintKey = AromaticSystemConstraintKey;
    type Atoms = Vec<AtomId>;

    fn id(&self) -> AromaticSystemId {
        match self {
            AromaticSystemDelta::Add { id, .. }
            | AromaticSystemDelta::Remove { id, .. }
            | AromaticSystemDelta::ModifyField { id, .. }
            | AromaticSystemDelta::ModifyConstraint { id, .. } => *id,
        }
    }

    fn split(self) -> EntityOp<Self> {
        match self {
            AromaticSystemDelta::Add { atoms, ast, .. } => EntityOp::Add { atoms, ast },
            AromaticSystemDelta::Remove { atoms, ast, .. } => EntityOp::Remove { atoms, ast },
            AromaticSystemDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            AromaticSystemDelta::ModifyConstraint { old, new, .. } => {
                EntityOp::ModifyConstraint { old, new }
            }
        }
    }

    fn rebuild(id: AromaticSystemId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { atoms, ast } => AromaticSystemDelta::Add { id, atoms, ast },
            EntityOp::Remove { atoms, ast } => AromaticSystemDelta::Remove { id, atoms, ast },
            EntityOp::ModifyField(change) => Self::modify_field(id, change),
            EntityOp::ModifyConstraint { old, new } => Self::modify_constraint(id, old, new),
        }
    }

    fn into_delta(self) -> Delta {
        Delta::AromaticSystem(self)
    }

    fn field_inverse(change: AromaticSystemFieldChange) -> AromaticSystemFieldChange {
        change.inverse()
    }

    fn constraint_key(constraint: &AromaticSystemConstraintAst) -> AromaticSystemConstraintKey {
        constraint.key()
    }

    fold_field_ops!(AromaticSystemFieldChange, { Electrons, Charge, Spin });
}

impl EntityPatch for MulticenterBondDelta {
    type Id = MulticenterBondId;
    type Ast = MulticenterBondAst;
    type FieldChange = MulticenterBondFieldChange;
    type Constraint = MulticenterBondConstraintAst;

    fn modify_field(id: MulticenterBondId, change: MulticenterBondFieldChange) -> Self {
        MulticenterBondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: MulticenterBondId,
        old: Option<MulticenterBondConstraintAst>,
        new: Option<MulticenterBondConstraintAst>,
    ) -> Self {
        MulticenterBondDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(
        MulticenterBondFieldChange,
        MulticenterBondAst,
        MulticenterBondConstraintAst,
        {
            Electrons => electrons,
            Charge => charge,
            Spin => spin,
        }
    );

    fn apply_constraint(
        ast: &mut MulticenterBondAst,
        old: Option<MulticenterBondConstraintAst>,
        new: Option<MulticenterBondConstraintAst>,
    ) -> Result<(), Contradiction> {
        ast.constraints.compare_and_set(old, new)
    }
}

impl EntityFold for MulticenterBondDelta {
    type ConstraintKey = MulticenterBondConstraintKey;
    type Atoms = Vec<AtomId>;

    fn id(&self) -> MulticenterBondId {
        match self {
            MulticenterBondDelta::Add { id, .. }
            | MulticenterBondDelta::Remove { id, .. }
            | MulticenterBondDelta::ModifyField { id, .. }
            | MulticenterBondDelta::ModifyConstraint { id, .. } => *id,
        }
    }

    fn split(self) -> EntityOp<Self> {
        match self {
            MulticenterBondDelta::Add { atoms, ast, .. } => EntityOp::Add { atoms, ast },
            MulticenterBondDelta::Remove { atoms, ast, .. } => EntityOp::Remove { atoms, ast },
            MulticenterBondDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            MulticenterBondDelta::ModifyConstraint { old, new, .. } => {
                EntityOp::ModifyConstraint { old, new }
            }
        }
    }

    fn rebuild(id: MulticenterBondId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { atoms, ast } => MulticenterBondDelta::Add { id, atoms, ast },
            EntityOp::Remove { atoms, ast } => MulticenterBondDelta::Remove { id, atoms, ast },
            EntityOp::ModifyField(change) => Self::modify_field(id, change),
            EntityOp::ModifyConstraint { old, new } => Self::modify_constraint(id, old, new),
        }
    }

    fn into_delta(self) -> Delta {
        Delta::MulticenterBond(self)
    }

    fn field_inverse(change: MulticenterBondFieldChange) -> MulticenterBondFieldChange {
        change.inverse()
    }

    fn constraint_key(constraint: &MulticenterBondConstraintAst) -> MulticenterBondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(MulticenterBondFieldChange, { Electrons, Charge, Spin });
}

impl EntityPatch for NoncovalentBondDelta {
    type Id = NoncovalentBondId;
    type Ast = NoncovalentBondAst;
    type FieldChange = NoncovalentBondFieldChange;
    type Constraint = NoncovalentBondConstraintAst;

    fn modify_field(id: NoncovalentBondId, change: NoncovalentBondFieldChange) -> Self {
        NoncovalentBondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: NoncovalentBondId,
        old: Option<NoncovalentBondConstraintAst>,
        new: Option<NoncovalentBondConstraintAst>,
    ) -> Self {
        NoncovalentBondDelta::ModifyConstraint { id, old, new }
    }

    // Hand-written (not `diff_field_ops!`): `NoncovalentBondConstraintAst` is uninhabited,
    // so the macro's constraint loop would be unreachable code.
    fn apply_field(
        ast: &mut NoncovalentBondAst,
        change: NoncovalentBondFieldChange,
    ) -> Result<(), Contradiction> {
        match change {
            NoncovalentBondFieldChange::Kind { old, new } => {
                if !ast.kind.canonical_eq(&old) {
                    return Err(Contradiction);
                }
                ast.kind = new;
            }
        }
        Ok(())
    }

    fn diff_field(
        lhs: &NoncovalentBondAst,
        rhs: &NoncovalentBondAst,
    ) -> Vec<NoncovalentBondFieldChange> {
        if !lhs.kind.canonical_eq(&rhs.kind) {
            vec![NoncovalentBondFieldChange::Kind {
                old: lhs.kind.clone(),
                new: rhs.kind.clone(),
            }]
        } else {
            Vec::new()
        }
    }

    fn diff_constraints(
        _lhs: &NoncovalentBondAst,
        _rhs: &NoncovalentBondAst,
    ) -> Vec<(
        Option<NoncovalentBondConstraintAst>,
        Option<NoncovalentBondConstraintAst>,
    )> {
        Vec::new()
    }

    /// `NoncovalentBondConstraintAst` is uninhabited, so `old`/`new` are always `None`.
    fn apply_constraint(
        _ast: &mut NoncovalentBondAst,
        old: Option<NoncovalentBondConstraintAst>,
        new: Option<NoncovalentBondConstraintAst>,
    ) -> Result<(), Contradiction> {
        debug_assert!(
            old.is_none() && new.is_none(),
            "noncovalent constraints are uninhabited"
        );
        Ok(())
    }
}

impl EntityFold for NoncovalentBondDelta {
    type ConstraintKey = NoncovalentBondConstraintKey;
    type Atoms = [AtomId; 2];

    fn id(&self) -> NoncovalentBondId {
        match self {
            NoncovalentBondDelta::Add { id, .. }
            | NoncovalentBondDelta::Remove { id, .. }
            | NoncovalentBondDelta::ModifyField { id, .. }
            | NoncovalentBondDelta::ModifyConstraint { id, .. } => *id,
        }
    }

    fn split(self) -> EntityOp<Self> {
        match self {
            NoncovalentBondDelta::Add { atoms, ast, .. } => EntityOp::Add { atoms, ast },
            NoncovalentBondDelta::Remove { atoms, ast, .. } => EntityOp::Remove { atoms, ast },
            NoncovalentBondDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            NoncovalentBondDelta::ModifyConstraint { old, new, .. } => {
                EntityOp::ModifyConstraint { old, new }
            }
        }
    }

    fn rebuild(id: NoncovalentBondId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { atoms, ast } => NoncovalentBondDelta::Add { id, atoms, ast },
            EntityOp::Remove { atoms, ast } => NoncovalentBondDelta::Remove { id, atoms, ast },
            EntityOp::ModifyField(change) => Self::modify_field(id, change),
            EntityOp::ModifyConstraint { old, new } => Self::modify_constraint(id, old, new),
        }
    }

    fn into_delta(self) -> Delta {
        Delta::NoncovalentBond(self)
    }

    fn field_inverse(change: NoncovalentBondFieldChange) -> NoncovalentBondFieldChange {
        change.inverse()
    }

    fn constraint_key(constraint: &NoncovalentBondConstraintAst) -> NoncovalentBondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(NoncovalentBondFieldChange, { Kind });
}

// Stereo entities impl only `EntityPatch` (diff / apply of the four DAMN arms), not `EntityFold`:
// the relative ops `Apply`/`Swap`/`Mirror` have no `EntityOp` image, so `canonicalize` folds stereo
// on a bespoke path (the four arms still route through these `diff`/`apply` methods).
impl EntityPatch for StereoAtomDelta {
    type Id = StereoAtomId;
    type Ast = StereoAtomAst;
    type FieldChange = StereoAtomFieldChange;
    type Constraint = StereoAtomConstraintAst;

    fn modify_field(id: StereoAtomId, change: StereoAtomFieldChange) -> Self {
        StereoAtomDelta::ModifyField { id, change }
    }

    /// Kind-less fallback (the trait signature has no kind); the real producer is the overridden
    /// `diff`, which stamps `kind` from the entity's config. Stereo's flow never uses this arm.
    fn modify_constraint(
        id: StereoAtomId,
        old: Option<StereoAtomConstraintAst>,
        new: Option<StereoAtomConstraintAst>,
    ) -> Self {
        StereoAtomDelta::ModifyConstraint {
            id,
            kind: None,
            old,
            new,
        }
    }

    diff_field_ops!(StereoAtomFieldChange, StereoAtomAst, StereoAtomConstraintAst, {
        Configuration => configuration,
    });

    /// Stamp each `ModifyConstraint` with the config's kind (the serialization context the
    /// constraint needs) — the default `diff` routes through `modify_constraint`, which can't.
    fn diff(id: StereoAtomId, lhs: &StereoAtomAst, rhs: &StereoAtomAst) -> Vec<Self> {
        let kind = lhs
            .configuration
            .kind()
            .or_else(|| rhs.configuration.kind());
        let mut out: Vec<Self> = Self::diff_field(lhs, rhs)
            .into_iter()
            .map(|change| StereoAtomDelta::ModifyField { id, change })
            .collect();
        out.extend(
            Self::diff_constraints(lhs, rhs)
                .into_iter()
                .map(|(old, new)| StereoAtomDelta::ModifyConstraint { id, kind, old, new }),
        );
        out
    }

    fn apply_constraint(
        ast: &mut StereoAtomAst,
        old: Option<StereoAtomConstraintAst>,
        new: Option<StereoAtomConstraintAst>,
    ) -> Result<(), Contradiction> {
        ast.constraints.compare_and_set(old, new)
    }
}

impl EntityPatch for StereoBondDelta {
    type Id = StereoBondId;
    type Ast = StereoBondAst;
    type FieldChange = StereoBondFieldChange;
    type Constraint = StereoBondConstraintAst;

    fn modify_field(id: StereoBondId, change: StereoBondFieldChange) -> Self {
        StereoBondDelta::ModifyField { id, change }
    }

    /// Kind-less fallback — see `StereoAtomDelta::modify_constraint`.
    fn modify_constraint(
        id: StereoBondId,
        old: Option<StereoBondConstraintAst>,
        new: Option<StereoBondConstraintAst>,
    ) -> Self {
        StereoBondDelta::ModifyConstraint {
            id,
            kind: None,
            old,
            new,
        }
    }

    diff_field_ops!(StereoBondFieldChange, StereoBondAst, StereoBondConstraintAst, {
        Configuration => configuration,
    });

    /// Stamp each `ModifyConstraint` with the config's kind — see `StereoAtomDelta::diff`.
    fn diff(id: StereoBondId, lhs: &StereoBondAst, rhs: &StereoBondAst) -> Vec<Self> {
        let kind = lhs
            .configuration
            .kind()
            .or_else(|| rhs.configuration.kind());
        let mut out: Vec<Self> = Self::diff_field(lhs, rhs)
            .into_iter()
            .map(|change| StereoBondDelta::ModifyField { id, change })
            .collect();
        out.extend(
            Self::diff_constraints(lhs, rhs)
                .into_iter()
                .map(|(old, new)| StereoBondDelta::ModifyConstraint { id, kind, old, new }),
        );
        out
    }

    fn apply_constraint(
        ast: &mut StereoBondAst,
        old: Option<StereoBondConstraintAst>,
        new: Option<StereoBondConstraintAst>,
    ) -> Result<(), Contradiction> {
        ast.constraints.compare_and_set(old, new)
    }
}

/// Apply a resolved per-entity change to a value AST, reusing the `EntityPatch` apply that
/// `canonicalize` uses. `ModifyField` / `ModifyConstraint` mutate the ast; `Add` / `Remove` are
/// no-ops (they carry a whole ast, not a change). Materializes the rhs-hand value of a
/// preserved entity for a `ReactionSpanAst`.
pub(crate) fn apply_atom_change(ast: &mut AtomAst, delta: &AtomDelta) -> Result<(), Contradiction> {
    match delta {
        AtomDelta::ModifyField { change, .. } => AtomDelta::apply_field(ast, change.clone()),
        AtomDelta::ModifyConstraint { old, new, .. } => {
            AtomDelta::apply_constraint(ast, old.clone(), new.clone())
        }
        AtomDelta::Add { .. } | AtomDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_bond_change(ast: &mut BondAst, delta: &BondDelta) -> Result<(), Contradiction> {
    match delta {
        BondDelta::ModifyField { change, .. } => BondDelta::apply_field(ast, change.clone()),
        BondDelta::ModifyConstraint { old, new, .. } => {
            BondDelta::apply_constraint(ast, old.clone(), new.clone())
        }
        BondDelta::Add { .. } | BondDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_dative_change(
    ast: &mut DativeBondAst,
    delta: &DativeBondDelta,
) -> Result<(), Contradiction> {
    match delta {
        DativeBondDelta::ModifyField { change, .. } => {
            DativeBondDelta::apply_field(ast, change.clone())
        }
        DativeBondDelta::ModifyConstraint { old, new, .. } => {
            DativeBondDelta::apply_constraint(ast, old.clone(), new.clone())
        }
        DativeBondDelta::Add { .. } | DativeBondDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_aromatic_change(
    ast: &mut AromaticSystemAst,
    delta: &AromaticSystemDelta,
) -> Result<(), Contradiction> {
    match delta {
        AromaticSystemDelta::ModifyField { change, .. } => {
            AromaticSystemDelta::apply_field(ast, change.clone())
        }
        AromaticSystemDelta::ModifyConstraint { old, new, .. } => {
            AromaticSystemDelta::apply_constraint(ast, old.clone(), new.clone())
        }
        AromaticSystemDelta::Add { .. } | AromaticSystemDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_multicenter_change(
    ast: &mut MulticenterBondAst,
    delta: &MulticenterBondDelta,
) -> Result<(), Contradiction> {
    match delta {
        MulticenterBondDelta::ModifyField { change, .. } => {
            MulticenterBondDelta::apply_field(ast, change.clone())
        }
        MulticenterBondDelta::ModifyConstraint { old, new, .. } => {
            MulticenterBondDelta::apply_constraint(ast, old.clone(), new.clone())
        }
        MulticenterBondDelta::Add { .. } | MulticenterBondDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_noncovalent_change(
    ast: &mut NoncovalentBondAst,
    delta: &NoncovalentBondDelta,
) -> Result<(), Contradiction> {
    match delta {
        NoncovalentBondDelta::ModifyField { change, .. } => {
            NoncovalentBondDelta::apply_field(ast, change.clone())
        }
        NoncovalentBondDelta::ModifyConstraint { old, new, .. } => {
            NoncovalentBondDelta::apply_constraint(ast, old.clone(), new.clone())
        }
        NoncovalentBondDelta::Add { .. } | NoncovalentBondDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_stereo_atom_change(
    ast: &mut StereoAtomAst,
    delta: &StereoAtomDelta,
) -> Result<(), Contradiction> {
    match delta {
        StereoAtomDelta::ModifyField { change, .. } => {
            StereoAtomDelta::apply_field(ast, change.clone())
        }
        StereoAtomDelta::ModifyConstraint { old, new, .. } => {
            StereoAtomDelta::apply_constraint(ast, old.clone(), new.clone())
        }
        StereoAtomDelta::Apply { permutation, .. } => {
            *ast = ast.apply(*permutation);
            Ok(())
        }
        StereoAtomDelta::Swap { .. } => {
            *ast = ast.swap();
            Ok(())
        }
        StereoAtomDelta::Mirror { .. } => {
            *ast = ast.mirror();
            Ok(())
        }
        StereoAtomDelta::Add { .. } | StereoAtomDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_stereo_bond_change(
    ast: &mut StereoBondAst,
    delta: &StereoBondDelta,
) -> Result<(), Contradiction> {
    match delta {
        StereoBondDelta::ModifyField { change, .. } => {
            StereoBondDelta::apply_field(ast, change.clone())
        }
        StereoBondDelta::ModifyConstraint { old, new, .. } => {
            StereoBondDelta::apply_constraint(ast, old.clone(), new.clone())
        }
        StereoBondDelta::Apply { permutation, .. } => {
            *ast = ast.apply(*permutation);
            Ok(())
        }
        StereoBondDelta::Swap { .. } => {
            *ast = ast.swap();
            Ok(())
        }
        StereoBondDelta::Mirror { .. } => {
            *ast = ast.mirror();
            Ok(())
        }
        StereoBondDelta::Add { .. } | StereoBondDelta::Remove { .. } => Ok(()),
    }
}

/// One config-affecting op of a preserved stereo entity: an absolute set, or a relative coset op
/// (its net permutation). Constraint / membership ops fold separately.
enum StereoConfigOp {
    Set {
        old: StereoConfigurationAst,
        new: StereoConfigurationAst,
    },
    Relative(Permutation),
}

/// The net config action after folding a preserved entity's config ops (rules ii–vi).
enum StereoConfigFold {
    Identity,
    Relative(Permutation),
    Set {
        old: StereoConfigurationAst,
        new: StereoConfigurationAst,
    },
}

/// Fold a preserved entity's config ops, in order, to one net action. Relatives compose by
/// permutation (ii/iii); a set absorbs a leading relative by pulling `old` back through its inverse
/// (vi) and transforms `new` through a trailing relative (v); two sets fuse (`new₁ == old₂`).
fn fold_stereo_config(
    kind: StereoKind,
    ops: Vec<StereoConfigOp>,
) -> Result<StereoConfigFold, Contradiction> {
    enum State {
        Relative(Permutation),
        Set {
            old: StereoConfigurationAst,
            new: StereoConfigurationAst,
        },
    }
    let mut state = State::Relative(Permutation::identity(kind.degree()));
    for op in ops {
        state = match (state, op) {
            (State::Relative(sigma), StereoConfigOp::Relative(p)) => {
                State::Relative(sigma.compose(p))
            }
            (State::Relative(sigma), StereoConfigOp::Set { old, new }) => State::Set {
                old: old.apply(sigma.inverse()),
                new,
            },
            (State::Set { old, new }, StereoConfigOp::Relative(p)) => State::Set {
                old,
                new: new.apply(p),
            },
            (State::Set { old, new }, StereoConfigOp::Set { old: o2, new: n2 }) => {
                if new.clone().canonicalize()? != o2.clone().canonicalize()? {
                    return Err(Contradiction);
                }
                State::Set { old, new: n2 }
            }
        };
    }
    Ok(match state {
        State::Relative(sigma) => match kind.canonicalize_permutation(sigma) {
            None => StereoConfigFold::Identity,
            Some(_) => StereoConfigFold::Relative(sigma),
        },
        State::Set { old, new } => {
            if old.clone().canonicalize()? == new.clone().canonicalize()? {
                StereoConfigFold::Identity
            } else {
                StereoConfigFold::Set { old, new }
            }
        }
    })
}

/// Fold one stereo atom's deltas to normal form (input order). Created: seed from `Add`, apply each
/// op to the ast, `Add`+`Remove` cancels. Preserved: fold config ops (`fold_stereo_config`) +
/// constraints (by key); a `Remove` reverts both onto the removed (original) ast.
fn fold_stereo_atom_group(
    id: StereoAtomId,
    group: Vec<StereoAtomDelta>,
) -> Result<Vec<StereoAtomDelta>, Contradiction> {
    if group
        .iter()
        .any(|d| matches!(d, StereoAtomDelta::Add { .. }))
    {
        let mut state: Option<(AtomId, Vec<StereoLigand>, StereoAtomAst)> = None;
        let mut removed = false;
        for delta in group {
            if removed {
                return Err(Contradiction);
            }
            match delta {
                StereoAtomDelta::Add {
                    site, ligands, ast, ..
                } => {
                    if state.is_some() {
                        return Err(Contradiction);
                    }
                    state = Some((site, ligands, ast));
                }
                StereoAtomDelta::Remove { .. } => {
                    if state.is_none() {
                        return Err(Contradiction);
                    }
                    state = None;
                    removed = true;
                }
                other => {
                    let (_, _, ast) = state.as_mut().ok_or(Contradiction)?;
                    apply_stereo_atom_change(ast, &other)?;
                }
            }
        }
        return Ok(match state {
            Some((site, ligands, ast)) => vec![StereoAtomDelta::Add {
                id,
                site,
                ligands,
                ast,
            }],
            None => Vec::new(),
        });
    }
    let mut kind: Option<StereoKind> = None;
    let mut config_ops: Vec<StereoConfigOp> = Vec::new();
    let mut constraints: HashMap<
        StereoAtomConstraintKey,
        (
            Option<StereoAtomConstraintAst>,
            Option<StereoAtomConstraintAst>,
        ),
    > = HashMap::new();
    let mut removed: Option<(AtomId, Vec<StereoLigand>, StereoAtomAst)> = None;
    for delta in group {
        if removed.is_some() {
            return Err(Contradiction);
        }
        match delta {
            StereoAtomDelta::Add { .. } => return Err(Contradiction),
            StereoAtomDelta::ModifyField {
                change: StereoAtomFieldChange::Configuration { old, new },
                ..
            } => {
                kind = kind.or_else(|| old.kind());
                config_ops.push(StereoConfigOp::Set { old, new });
            }
            StereoAtomDelta::Apply {
                kind: k,
                permutation,
                ..
            } => {
                kind = Some(k);
                config_ops.push(StereoConfigOp::Relative(permutation));
            }
            StereoAtomDelta::Swap { kind: k, .. } => {
                kind = Some(k);
                config_ops.push(StereoConfigOp::Relative(k.involution()));
            }
            StereoAtomDelta::Mirror { kind: k, .. } => {
                kind = Some(k);
                config_ops.push(StereoConfigOp::Relative(k.mirror_permutation()));
            }
            StereoAtomDelta::ModifyConstraint {
                kind: constraint_kind,
                old,
                new,
                ..
            } => {
                kind = kind.or(constraint_kind);
                let key = match old.as_ref().or(new.as_ref()) {
                    Some(c) => c.key(),
                    None => continue,
                };
                match constraints.remove(&key) {
                    Some((first_old, prev_new)) => {
                        if !options_canonical_eq(&prev_new, &old) {
                            return Err(Contradiction);
                        }
                        constraints.insert(key, (first_old, new));
                    }
                    None => {
                        constraints.insert(key, (old, new));
                    }
                }
            }
            StereoAtomDelta::Remove {
                site, ligands, ast, ..
            } => {
                removed = Some((site, ligands, ast));
            }
        }
    }
    let config = match kind {
        Some(k) => fold_stereo_config(k, config_ops)?,
        None => StereoConfigFold::Identity,
    };
    if let Some((site, ligands, mut ast)) = removed {
        match config {
            StereoConfigFold::Identity => {}
            StereoConfigFold::Relative(sigma) => ast = ast.apply(sigma.inverse()),
            StereoConfigFold::Set { old, new } => {
                if ast.configuration.clone().canonicalize()? != new.clone().canonicalize()? {
                    return Err(Contradiction);
                }
                ast.configuration = old;
            }
        }
        for (_key, (old, new)) in constraints {
            StereoAtomDelta::apply_constraint(&mut ast, new, old)?;
        }
        return Ok(vec![StereoAtomDelta::Remove {
            id,
            site,
            ligands,
            ast,
        }]);
    }
    let mut out = Vec::new();
    match config {
        StereoConfigFold::Identity => {}
        StereoConfigFold::Relative(sigma) => {
            let k = kind.expect("a relative fold implies a kind");
            match k.canonicalize_permutation(sigma) {
                None => {}
                Some(CosetOp::Swap) => out.push(StereoAtomDelta::Swap { id, kind: k }),
                Some(CosetOp::Mirror) => out.push(StereoAtomDelta::Mirror { id, kind: k }),
                Some(CosetOp::Apply(g)) => out.push(StereoAtomDelta::Apply {
                    id,
                    kind: k,
                    permutation: g,
                }),
            }
        }
        StereoConfigFold::Set { old, new } => out.push(StereoAtomDelta::ModifyField {
            id,
            change: StereoAtomFieldChange::Configuration { old, new },
        }),
    }
    for (_key, (old, new)) in constraints {
        if !options_canonical_eq(&old, &new) {
            out.push(StereoAtomDelta::ModifyConstraint { id, kind, old, new });
        }
    }
    Ok(out)
}

/// Fold one stereo bond's deltas to normal form — the `fold_stereo_atom_group` twin (bond ids/ast).
fn fold_stereo_bond_group(
    id: StereoBondId,
    group: Vec<StereoBondDelta>,
) -> Result<Vec<StereoBondDelta>, Contradiction> {
    if group
        .iter()
        .any(|d| matches!(d, StereoBondDelta::Add { .. }))
    {
        let mut state: Option<(BondId, Vec<StereoLigand>, StereoBondAst)> = None;
        let mut removed = false;
        for delta in group {
            if removed {
                return Err(Contradiction);
            }
            match delta {
                StereoBondDelta::Add {
                    site, ligands, ast, ..
                } => {
                    if state.is_some() {
                        return Err(Contradiction);
                    }
                    state = Some((site, ligands, ast));
                }
                StereoBondDelta::Remove { .. } => {
                    if state.is_none() {
                        return Err(Contradiction);
                    }
                    state = None;
                    removed = true;
                }
                other => {
                    let (_, _, ast) = state.as_mut().ok_or(Contradiction)?;
                    apply_stereo_bond_change(ast, &other)?;
                }
            }
        }
        return Ok(match state {
            Some((site, ligands, ast)) => vec![StereoBondDelta::Add {
                id,
                site,
                ligands,
                ast,
            }],
            None => Vec::new(),
        });
    }
    let mut kind: Option<StereoKind> = None;
    let mut config_ops: Vec<StereoConfigOp> = Vec::new();
    let mut constraints: HashMap<
        StereoBondConstraintKey,
        (
            Option<StereoBondConstraintAst>,
            Option<StereoBondConstraintAst>,
        ),
    > = HashMap::new();
    let mut removed: Option<(BondId, Vec<StereoLigand>, StereoBondAst)> = None;
    for delta in group {
        if removed.is_some() {
            return Err(Contradiction);
        }
        match delta {
            StereoBondDelta::Add { .. } => return Err(Contradiction),
            StereoBondDelta::ModifyField {
                change: StereoBondFieldChange::Configuration { old, new },
                ..
            } => {
                kind = kind.or_else(|| old.kind());
                config_ops.push(StereoConfigOp::Set { old, new });
            }
            StereoBondDelta::Apply {
                kind: k,
                permutation,
                ..
            } => {
                kind = Some(k);
                config_ops.push(StereoConfigOp::Relative(permutation));
            }
            StereoBondDelta::Swap { kind: k, .. } => {
                kind = Some(k);
                config_ops.push(StereoConfigOp::Relative(k.involution()));
            }
            StereoBondDelta::Mirror { kind: k, .. } => {
                kind = Some(k);
                config_ops.push(StereoConfigOp::Relative(k.mirror_permutation()));
            }
            StereoBondDelta::ModifyConstraint {
                kind: constraint_kind,
                old,
                new,
                ..
            } => {
                kind = kind.or(constraint_kind);
                let key = match old.as_ref().or(new.as_ref()) {
                    Some(c) => c.key(),
                    None => continue,
                };
                match constraints.remove(&key) {
                    Some((first_old, prev_new)) => {
                        if !options_canonical_eq(&prev_new, &old) {
                            return Err(Contradiction);
                        }
                        constraints.insert(key, (first_old, new));
                    }
                    None => {
                        constraints.insert(key, (old, new));
                    }
                }
            }
            StereoBondDelta::Remove {
                site, ligands, ast, ..
            } => {
                removed = Some((site, ligands, ast));
            }
        }
    }
    let config = match kind {
        Some(k) => fold_stereo_config(k, config_ops)?,
        None => StereoConfigFold::Identity,
    };
    if let Some((site, ligands, mut ast)) = removed {
        match config {
            StereoConfigFold::Identity => {}
            StereoConfigFold::Relative(sigma) => ast = ast.apply(sigma.inverse()),
            StereoConfigFold::Set { old, new } => {
                if ast.configuration.clone().canonicalize()? != new.clone().canonicalize()? {
                    return Err(Contradiction);
                }
                ast.configuration = old;
            }
        }
        for (_key, (old, new)) in constraints {
            StereoBondDelta::apply_constraint(&mut ast, new, old)?;
        }
        return Ok(vec![StereoBondDelta::Remove {
            id,
            site,
            ligands,
            ast,
        }]);
    }
    let mut out = Vec::new();
    match config {
        StereoConfigFold::Identity => {}
        StereoConfigFold::Relative(sigma) => {
            let k = kind.expect("a relative fold implies a kind");
            match k.canonicalize_permutation(sigma) {
                None => {}
                Some(CosetOp::Swap) => out.push(StereoBondDelta::Swap { id, kind: k }),
                Some(CosetOp::Mirror) => out.push(StereoBondDelta::Mirror { id, kind: k }),
                Some(CosetOp::Apply(g)) => out.push(StereoBondDelta::Apply {
                    id,
                    kind: k,
                    permutation: g,
                }),
            }
        }
        StereoConfigFold::Set { old, new } => out.push(StereoBondDelta::ModifyField {
            id,
            change: StereoBondFieldChange::Configuration { old, new },
        }),
    }
    for (_key, (old, new)) in constraints {
        if !options_canonical_eq(&old, &new) {
            out.push(StereoBondDelta::ModifyConstraint { id, kind, old, new });
        }
    }
    Ok(out)
}

/// Re-anchor a delta's ids and participant atoms through a total id relabeling. Used to move
/// deltas between id spaces (reverse re-anchoring, composition). The relabeling must cover every id
/// the delta references. Overlay participants on `Unordered` factors are re-sorted to canonical
/// order; aromatic/multicenter electrons are permuted to stay aligned with their atoms.
pub fn remap_delta(delta: Delta, map: &IdRemapping) -> Delta {
    match delta {
        Delta::Atom(a) => Delta::Atom(match a {
            AtomDelta::Add { id, ast } => AtomDelta::Add {
                id: map.map_atom(id),
                ast,
            },
            AtomDelta::Remove { id, ast } => AtomDelta::Remove {
                id: map.map_atom(id),
                ast,
            },
            AtomDelta::ModifyField { id, change } => AtomDelta::ModifyField {
                id: map.map_atom(id),
                change,
            },
            AtomDelta::ModifyConstraint { id, old, new } => AtomDelta::ModifyConstraint {
                id: map.map_atom(id),
                old,
                new,
            },
        }),
        Delta::Bond(b) => Delta::Bond(match b {
            BondDelta::Add { id, atoms, ast } => BondDelta::Add {
                id: map.map_bond(id),
                atoms: [map.map_atom(atoms[0]), map.map_atom(atoms[1])],
                ast,
            },
            BondDelta::Remove { id, atoms, ast } => BondDelta::Remove {
                id: map.map_bond(id),
                atoms: [map.map_atom(atoms[0]), map.map_atom(atoms[1])],
                ast,
            },
            BondDelta::ModifyField { id, change } => BondDelta::ModifyField {
                id: map.map_bond(id),
                change,
            },
            BondDelta::ModifyConstraint { id, old, new } => BondDelta::ModifyConstraint {
                id: map.map_bond(id),
                old,
                new,
            },
        }),
        Delta::DativeBond(d) => Delta::DativeBond(match d {
            // Donors are the unordered factor with no per-participant ast data, so canonicalize the
            // order after remap (acceptor is the single ordered factor). No permutation to track.
            DativeBondDelta::Add {
                id,
                donors,
                acceptor,
                ast,
            } => {
                let mut donors: Vec<AtomId> = donors.iter().map(|a| map.map_atom(*a)).collect();
                Unordered::canonicalize(&mut donors);
                DativeBondDelta::Add {
                    id: map.map_dative(id),
                    donors,
                    acceptor: map.map_atom(acceptor),
                    ast,
                }
            }
            DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                ast,
            } => {
                let mut donors: Vec<AtomId> = donors.iter().map(|a| map.map_atom(*a)).collect();
                Unordered::canonicalize(&mut donors);
                DativeBondDelta::Remove {
                    id: map.map_dative(id),
                    donors,
                    acceptor: map.map_atom(acceptor),
                    ast,
                }
            }
            DativeBondDelta::ModifyField { id, change } => DativeBondDelta::ModifyField {
                id: map.map_dative(id),
                change,
            },
            DativeBondDelta::ModifyConstraint { id, old, new } => {
                DativeBondDelta::ModifyConstraint {
                    id: map.map_dative(id),
                    old,
                    new,
                }
            }
        }),
        Delta::AromaticSystem(a) => Delta::AromaticSystem(match a {
            AromaticSystemDelta::Add { id, atoms, mut ast } => {
                let mut atoms: Vec<AtomId> = atoms.iter().map(|a| map.map_atom(*a)).collect();
                let order = Unordered::canonicalize_positions(&mut atoms);
                ast.permute(&order);
                AromaticSystemDelta::Add {
                    id: map.map_aromatic(id),
                    atoms,
                    ast,
                }
            }
            AromaticSystemDelta::Remove { id, atoms, mut ast } => {
                let mut atoms: Vec<AtomId> = atoms.iter().map(|a| map.map_atom(*a)).collect();
                let order = Unordered::canonicalize_positions(&mut atoms);
                ast.permute(&order);
                AromaticSystemDelta::Remove {
                    id: map.map_aromatic(id),
                    atoms,
                    ast,
                }
            }
            AromaticSystemDelta::ModifyField { id, change } => AromaticSystemDelta::ModifyField {
                id: map.map_aromatic(id),
                change,
            },
            AromaticSystemDelta::ModifyConstraint { id, old, new } => {
                AromaticSystemDelta::ModifyConstraint {
                    id: map.map_aromatic(id),
                    old,
                    new,
                }
            }
        }),
        Delta::MulticenterBond(m) => Delta::MulticenterBond(match m {
            MulticenterBondDelta::Add { id, atoms, mut ast } => {
                let mut atoms: Vec<AtomId> = atoms.iter().map(|a| map.map_atom(*a)).collect();
                let order = Unordered::canonicalize_positions(&mut atoms);
                ast.permute(&order);
                MulticenterBondDelta::Add {
                    id: map.map_multicenter(id),
                    atoms,
                    ast,
                }
            }
            MulticenterBondDelta::Remove { id, atoms, mut ast } => {
                let mut atoms: Vec<AtomId> = atoms.iter().map(|a| map.map_atom(*a)).collect();
                let order = Unordered::canonicalize_positions(&mut atoms);
                ast.permute(&order);
                MulticenterBondDelta::Remove {
                    id: map.map_multicenter(id),
                    atoms,
                    ast,
                }
            }
            MulticenterBondDelta::ModifyField { id, change } => MulticenterBondDelta::ModifyField {
                id: map.map_multicenter(id),
                change,
            },
            MulticenterBondDelta::ModifyConstraint { id, old, new } => {
                MulticenterBondDelta::ModifyConstraint {
                    id: map.map_multicenter(id),
                    old,
                    new,
                }
            }
        }),
        Delta::NoncovalentBond(n) => Delta::NoncovalentBond(match n {
            // Both participants are the unordered factor with no per-participant ast data, so
            // canonicalize the order after remap. No permutation to track.
            NoncovalentBondDelta::Add { id, atoms, ast } => {
                let mut atoms = [map.map_atom(atoms[0]), map.map_atom(atoms[1])];
                Unordered::canonicalize(&mut atoms);
                NoncovalentBondDelta::Add {
                    id: map.map_noncovalent(id),
                    atoms,
                    ast,
                }
            }
            NoncovalentBondDelta::Remove { id, atoms, ast } => {
                let mut atoms = [map.map_atom(atoms[0]), map.map_atom(atoms[1])];
                Unordered::canonicalize(&mut atoms);
                NoncovalentBondDelta::Remove {
                    id: map.map_noncovalent(id),
                    atoms,
                    ast,
                }
            }
            NoncovalentBondDelta::ModifyField { id, change } => NoncovalentBondDelta::ModifyField {
                id: map.map_noncovalent(id),
                change,
            },
            NoncovalentBondDelta::ModifyConstraint { id, old, new } => {
                NoncovalentBondDelta::ModifyConstraint {
                    id: map.map_noncovalent(id),
                    old,
                    new,
                }
            }
        }),
        // Stereo: ids relabel and the ligand-frame atom ids relabel in place; the frame is `Ordered`
        // (not re-sorted on remap), so the coset stays valid — no `transform_frame`. `Apply`'s
        // permutation is position-space, and constraints reference ligand positions — both untouched.
        Delta::StereoAtom(s) => Delta::StereoAtom(match s {
            StereoAtomDelta::Add {
                id,
                site,
                ligands,
                ast,
            } => StereoAtomDelta::Add {
                id: map.map_stereo_atom(id),
                site: map.map_atom(site),
                ligands: ligands
                    .into_iter()
                    .map(|l| StereoLigand::new(map.map_atom(l.atom_id), l.kind))
                    .collect(),
                ast,
            },
            StereoAtomDelta::Remove {
                id,
                site,
                ligands,
                ast,
            } => StereoAtomDelta::Remove {
                id: map.map_stereo_atom(id),
                site: map.map_atom(site),
                ligands: ligands
                    .into_iter()
                    .map(|l| StereoLigand::new(map.map_atom(l.atom_id), l.kind))
                    .collect(),
                ast,
            },
            StereoAtomDelta::ModifyField { id, change } => StereoAtomDelta::ModifyField {
                id: map.map_stereo_atom(id),
                change,
            },
            StereoAtomDelta::ModifyConstraint { id, kind, old, new } => {
                StereoAtomDelta::ModifyConstraint {
                    id: map.map_stereo_atom(id),
                    kind,
                    old,
                    new,
                }
            }
            StereoAtomDelta::Apply {
                id,
                kind,
                permutation,
            } => StereoAtomDelta::Apply {
                id: map.map_stereo_atom(id),
                kind,
                permutation,
            },
            StereoAtomDelta::Swap { id, kind } => StereoAtomDelta::Swap {
                id: map.map_stereo_atom(id),
                kind,
            },
            StereoAtomDelta::Mirror { id, kind } => StereoAtomDelta::Mirror {
                id: map.map_stereo_atom(id),
                kind,
            },
        }),
        Delta::StereoBond(s) => Delta::StereoBond(match s {
            StereoBondDelta::Add {
                id,
                site,
                ligands,
                ast,
            } => StereoBondDelta::Add {
                id: map.map_stereo_bond(id),
                site: map.map_bond(site),
                ligands: ligands
                    .into_iter()
                    .map(|l| StereoLigand::new(map.map_atom(l.atom_id), l.kind))
                    .collect(),
                ast,
            },
            StereoBondDelta::Remove {
                id,
                site,
                ligands,
                ast,
            } => StereoBondDelta::Remove {
                id: map.map_stereo_bond(id),
                site: map.map_bond(site),
                ligands: ligands
                    .into_iter()
                    .map(|l| StereoLigand::new(map.map_atom(l.atom_id), l.kind))
                    .collect(),
                ast,
            },
            StereoBondDelta::ModifyField { id, change } => StereoBondDelta::ModifyField {
                id: map.map_stereo_bond(id),
                change,
            },
            StereoBondDelta::ModifyConstraint { id, kind, old, new } => {
                StereoBondDelta::ModifyConstraint {
                    id: map.map_stereo_bond(id),
                    kind,
                    old,
                    new,
                }
            }
            StereoBondDelta::Apply {
                id,
                kind,
                permutation,
            } => StereoBondDelta::Apply {
                id: map.map_stereo_bond(id),
                kind,
                permutation,
            },
            StereoBondDelta::Swap { id, kind } => StereoBondDelta::Swap {
                id: map.map_stereo_bond(id),
                kind,
            },
            StereoBondDelta::Mirror { id, kind } => StereoBondDelta::Mirror {
                id: map.map_stereo_bond(id),
                kind,
            },
        }),
        Delta::Constraint(c) => Delta::Constraint(c),
    }
}

/// The resolved-delta collection.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct Deltas(Vec<Delta>);

impl Deltas {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn as_slice(&self) -> &[Delta] {
        &self.0
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn iter(&self) -> Iter<'_, Delta> {
        self.0.iter()
    }

    pub fn iter_mut(&mut self) -> IterMut<'_, Delta> {
        self.0.iter_mut()
    }

    pub fn push(&mut self, delta: Delta) {
        self.0.push(delta);
    }
}

impl FromIterator<Delta> for Deltas {
    fn from_iter<I: IntoIterator<Item = Delta>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl Canonicalize for Deltas {
    /// Per-entity fold to the normal form, then a stable sort. Different entities are
    /// independent and each entity's fold is deterministic over input order, so the result is
    /// a unique normal form; sequence order is not stored. `Err(Contradiction)` on an
    /// inconsistent set.
    #[allow(clippy::mutable_key_type)]
    fn canonicalize(self) -> Result<Self, Contradiction> {
        let mut atoms: HashMap<AtomId, Vec<AtomDelta>> = HashMap::new();
        let mut bonds: HashMap<BondId, Vec<BondDelta>> = HashMap::new();
        let mut dative: HashMap<DativeBondId, Vec<DativeBondDelta>> = HashMap::new();
        let mut aromatic: HashMap<AromaticSystemId, Vec<AromaticSystemDelta>> = HashMap::new();
        let mut multicenter: HashMap<MulticenterBondId, Vec<MulticenterBondDelta>> = HashMap::new();
        let mut noncovalent: HashMap<NoncovalentBondId, Vec<NoncovalentBondDelta>> = HashMap::new();
        let mut stereo_atom: HashMap<StereoAtomId, Vec<StereoAtomDelta>> = HashMap::new();
        let mut stereo_bond: HashMap<StereoBondId, Vec<StereoBondDelta>> = HashMap::new();
        let mut constraints: Vec<ConstraintDelta> = Vec::new();
        for delta in self.0 {
            match delta {
                Delta::Atom(d) => atoms.entry(d.id()).or_default().push(d),
                Delta::Bond(d) => bonds.entry(d.id()).or_default().push(d),
                Delta::DativeBond(d) => dative.entry(d.id()).or_default().push(d),
                Delta::AromaticSystem(d) => aromatic.entry(d.id()).or_default().push(d),
                Delta::MulticenterBond(d) => multicenter.entry(d.id()).or_default().push(d),
                Delta::NoncovalentBond(d) => noncovalent.entry(d.id()).or_default().push(d),
                Delta::StereoAtom(d) => stereo_atom.entry(d.id()).or_default().push(d),
                Delta::StereoBond(d) => stereo_bond.entry(d.id()).or_default().push(d),
                Delta::Constraint(d) => constraints.push(d),
            }
        }

        let mut out: Vec<Delta> = Vec::new();
        let mut removed_atoms: HashSet<AtomId> = HashSet::new();
        for (id, group) in atoms {
            let folded = fold_group::<AtomDelta>(id, group)?;
            if folded.iter().any(|d| matches!(d, AtomDelta::Remove { .. })) {
                removed_atoms.insert(id);
            }
            out.extend(folded.into_iter().map(Delta::Atom));
        }
        for (id, group) in bonds {
            let folded = fold_group::<BondDelta>(id, group)?;
            for delta in &folded {
                if let BondDelta::Add { atoms, .. } = delta {
                    if atoms.iter().any(|atom| removed_atoms.contains(atom)) {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::Bond));
        }
        // Overlay families: same fold; a created overlay must not reference a net-removed atom.
        for (id, group) in dative {
            let folded = fold_group::<DativeBondDelta>(id, group)?;
            for delta in &folded {
                if let DativeBondDelta::Add {
                    donors, acceptor, ..
                } = delta
                {
                    if donors
                        .iter()
                        .chain(iter::once(acceptor))
                        .any(|atom| removed_atoms.contains(atom))
                    {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::DativeBond));
        }
        for (id, group) in aromatic {
            let folded = fold_group::<AromaticSystemDelta>(id, group)?;
            for delta in &folded {
                if let AromaticSystemDelta::Add { atoms, .. } = delta {
                    if atoms.iter().any(|atom| removed_atoms.contains(atom)) {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::AromaticSystem));
        }
        for (id, group) in multicenter {
            let folded = fold_group::<MulticenterBondDelta>(id, group)?;
            for delta in &folded {
                if let MulticenterBondDelta::Add { atoms, .. } = delta {
                    if atoms.iter().any(|atom| removed_atoms.contains(atom)) {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::MulticenterBond));
        }
        for (id, group) in noncovalent {
            let folded = fold_group::<NoncovalentBondDelta>(id, group)?;
            for delta in &folded {
                if let NoncovalentBondDelta::Add { atoms, .. } = delta {
                    if atoms.iter().any(|atom| removed_atoms.contains(atom)) {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::NoncovalentBond));
        }
        // Stereo: bespoke fold (coset ops), same created-references-removed-atom guard (site +
        // ligands for a stereo atom; ligands for a stereo bond, whose site is a bond).
        for (id, group) in stereo_atom {
            let folded = fold_stereo_atom_group(id, group)?;
            for delta in &folded {
                if let StereoAtomDelta::Add { site, ligands, .. } = delta {
                    if iter::once(site)
                        .chain(ligands.iter().map(|l| &l.atom_id))
                        .any(|atom| removed_atoms.contains(atom))
                    {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::StereoAtom));
        }
        for (id, group) in stereo_bond {
            let folded = fold_stereo_bond_group(id, group)?;
            for delta in &folded {
                if let StereoBondDelta::Add { ligands, .. } = delta {
                    if ligands.iter().any(|l| removed_atoms.contains(&l.atom_id)) {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::StereoBond));
        }
        // Molecule-level constraints are a multiset: net multiplicity per constraint
        // (`Add`/`Remove` cancel one-for-one; duplicates are kept, not deduped).
        let mut net: BTreeMap<Constraint, i64> = BTreeMap::new();
        for delta in constraints {
            match delta {
                ConstraintDelta::Add(constraint) => *net.entry(constraint).or_insert(0) += 1,
                ConstraintDelta::Remove(constraint) => *net.entry(constraint).or_insert(0) -= 1,
            }
        }
        for (constraint, count) in net {
            if count > 0 {
                for _ in 0..count {
                    out.push(Delta::Constraint(ConstraintDelta::Add(constraint.clone())));
                }
            } else if count < 0 {
                for _ in 0..(-count) {
                    out.push(Delta::Constraint(ConstraintDelta::Remove(
                        constraint.clone(),
                    )));
                }
            }
        }

        out.sort();
        Ok(Self(out))
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::constraint::MoleculeConstraint;
    use super::super::noncovalent::NoncovalentBondKind;
    use super::super::value::ValueAst;
    use super::*;
    use crate::ast::{
        BooleanAst, StereoConfigurationAst, StereoCosetAst, StereoKind, StereoLigandKind,
        StereogenicityAst,
    };

    #[rstest]
    #[case::add_remove(
        AtomDelta::Add { id: AtomId(0), ast: AtomAst::from_element(Element::C) },
        AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::C) }
    )]
    #[case::set_field(
        AtomDelta::ModifyField {
            id: AtomId(1),
            change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(1) },
        },
        AtomDelta::ModifyField {
            id: AtomId(1),
            change: AtomFieldChange::Charge { old: ValueAst::Lit(1), new: ValueAst::Lit(0) },
        }
    )]
    #[case::set_constraint(
        AtomDelta::ModifyConstraint {
            id: AtomId(2),
            old: Some(AtomConstraintAst::Valence(ValueAst::Lit(4))),
            new: Some(AtomConstraintAst::Valence(ValueAst::Lit(3))),
        },
        AtomDelta::ModifyConstraint {
            id: AtomId(2),
            old: Some(AtomConstraintAst::Valence(ValueAst::Lit(3))),
            new: Some(AtomConstraintAst::Valence(ValueAst::Lit(4))),
        }
    )]
    fn test_atom_delta_inverse(#[case] input: AtomDelta, #[case] expected: AtomDelta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    #[case::singleton_set(ValueAst::Lit(1), ValueAst::lit_set([1]))]
    fn test_atom_delta_diff_canonical(#[case] lhs: ValueAst, #[case] rhs: ValueAst) {
        // Canonically-equal charges that are structurally distinct → `diff` emits nothing.
        let lhs = AtomAst::from_element(Element::C).with_charge(lhs);
        let rhs = AtomAst::from_element(Element::C).with_charge(rhs);
        assert_eq!(AtomDelta::diff(AtomId(0), &lhs, &rhs), Vec::new());
    }

    #[rstest]
    #[case::add_remove(
        BondDelta::Add {
            id: BondId(0),
            atoms: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        },
        BondDelta::Remove {
            id: BondId(0),
            atoms: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        }
    )]
    #[case::set_field(
        BondDelta::ModifyField {
            id: BondId(2),
            change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
        },
        BondDelta::ModifyField {
            id: BondId(2),
            change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(1) },
        }
    )]
    #[case::set_constraint(
        BondDelta::ModifyConstraint {
            id: BondId(3),
            old: None,
            new: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        },
        BondDelta::ModifyConstraint {
            id: BondId(3),
            old: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
            new: None,
        }
    )]
    fn test_bond_delta_inverse(#[case] input: BondDelta, #[case] expected: BondDelta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    #[case::add_remove(
        StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            ast: StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
        },
        StereoAtomDelta::Remove {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            ast: StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
        }
    )]
    #[case::set_field(
        StereoAtomDelta::ModifyField {
            id: StereoAtomId(1),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)),
                new: StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
            },
        },
        StereoAtomDelta::ModifyField {
            id: StereoAtomId(1),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
                new: StereoConfigurationAst::Kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)),
            },
        }
    )]
    #[case::set_constraint(
        StereoAtomDelta::ModifyConstraint {
            id: StereoAtomId(2),
            kind: Some(StereoKind::Tetrahedral),
            old: None,
            new: Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)),
        },
        StereoAtomDelta::ModifyConstraint {
            id: StereoAtomId(2),
            kind: Some(StereoKind::Tetrahedral),
            old: Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)),
            new: None,
        }
    )]
    #[case::apply(
        StereoAtomDelta::Apply { id: StereoAtomId(3), kind: StereoKind::Tetrahedral, permutation: Permutation::from_image(4, &[1, 2, 0, 3]) },
        StereoAtomDelta::Apply { id: StereoAtomId(3), kind: StereoKind::Tetrahedral, permutation: Permutation::from_image(4, &[2, 0, 1, 3]) }
    )]
    #[case::swap(
        StereoAtomDelta::Swap { id: StereoAtomId(4), kind: StereoKind::Tetrahedral },
        StereoAtomDelta::Swap { id: StereoAtomId(4), kind: StereoKind::Tetrahedral }
    )]
    #[case::mirror(
        StereoAtomDelta::Mirror { id: StereoAtomId(5), kind: StereoKind::Tetrahedral },
        StereoAtomDelta::Mirror { id: StereoAtomId(5), kind: StereoKind::Tetrahedral }
    )]
    fn test_stereo_atom_delta_inverse(
        #[case] input: StereoAtomDelta,
        #[case] expected: StereoAtomDelta,
    ) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    #[case::add_remove(
        StereoBondDelta::Add {
            id: StereoBondId(0),
            site: BondId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            ast: StereoBondAst::new(StereoKind::CisTrans, 0u32),
        },
        StereoBondDelta::Remove {
            id: StereoBondId(0),
            site: BondId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            ast: StereoBondAst::new(StereoKind::CisTrans, 0u32),
        }
    )]
    #[case::apply(
        StereoBondDelta::Apply { id: StereoBondId(1), kind: StereoKind::CisTrans, permutation: Permutation::from_image(4, &[1, 2, 0, 3]) },
        StereoBondDelta::Apply { id: StereoBondId(1), kind: StereoKind::CisTrans, permutation: Permutation::from_image(4, &[2, 0, 1, 3]) }
    )]
    #[case::swap(
        StereoBondDelta::Swap { id: StereoBondId(2), kind: StereoKind::CisTrans },
        StereoBondDelta::Swap { id: StereoBondId(2), kind: StereoKind::CisTrans }
    )]
    #[case::mirror(
        StereoBondDelta::Mirror { id: StereoBondId(3), kind: StereoKind::CisTrans },
        StereoBondDelta::Mirror { id: StereoBondId(3), kind: StereoKind::CisTrans }
    )]
    fn test_stereo_bond_delta_inverse(
        #[case] input: StereoBondDelta,
        #[case] expected: StereoBondDelta,
    ) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    fn test_stereo_atom_delta_diff() {
        assert_eq!(
            StereoAtomDelta::diff(
                StereoAtomId(0),
                &StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
                &StereoAtomAst::new(StereoKind::Tetrahedral, 1u32),
            ),
            vec![StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCosetAst::Lit(0)
                    ),
                    new: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCosetAst::Lit(1)
                    ),
                },
            }],
        );
    }

    #[rstest]
    fn test_stereo_atom_delta_apply_field() {
        let mut ast = StereoAtomAst::new(StereoKind::Tetrahedral, 0u32);
        StereoAtomDelta::apply_field(
            &mut ast,
            StereoAtomFieldChange::Configuration {
                old: StereoConfigurationAst::Kinded(
                    StereoKind::Tetrahedral,
                    StereoCosetAst::Lit(0),
                ),
                new: StereoConfigurationAst::Kinded(
                    StereoKind::Tetrahedral,
                    StereoCosetAst::Lit(1),
                ),
            },
        )
        .unwrap();
        assert_eq!(ast, StereoAtomAst::new(StereoKind::Tetrahedral, 1u32));
    }

    #[rstest]
    fn test_stereo_atom_delta_apply_field_error() {
        let mut ast = StereoAtomAst::new(StereoKind::Tetrahedral, 1u32);
        assert_eq!(
            StereoAtomDelta::apply_field(
                &mut ast,
                StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCosetAst::Lit(0)
                    ),
                    new: StereoConfigurationAst::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCosetAst::Lit(1)
                    ),
                },
            ),
            Err(Contradiction),
        );
    }

    #[rstest]
    // ii: swap/mirror are involutions — two of the same cancel.
    #[case::swap_swap(
        vec![
            StereoAtomDelta::Swap { id: StereoAtomId(0), kind: StereoKind::Tetrahedral },
            StereoAtomDelta::Swap { id: StereoAtomId(0), kind: StereoKind::Tetrahedral },
        ],
        vec![],
    )]
    #[case::mirror_mirror(
        vec![
            StereoAtomDelta::Mirror { id: StereoAtomId(0), kind: StereoKind::Tetrahedral },
            StereoAtomDelta::Mirror { id: StereoAtomId(0), kind: StereoKind::Tetrahedral },
        ],
        vec![],
    )]
    // iii: apply composes by permutation — a transposition twice is the identity.
    #[case::apply_apply_identity(
        vec![
            StereoAtomDelta::Apply { id: StereoAtomId(0), kind: StereoKind::Tetrahedral, permutation: Permutation::from_image(4, &[1, 0, 2, 3]) },
            StereoAtomDelta::Apply { id: StereoAtomId(0), kind: StereoKind::Tetrahedral, permutation: Permutation::from_image(4, &[1, 0, 2, 3]) },
        ],
        vec![],
    )]
    // iv: swap∘mirror — both are the improper generator for a chiral kind — composes to identity.
    #[case::swap_mirror(
        vec![
            StereoAtomDelta::Swap { id: StereoAtomId(0), kind: StereoKind::Tetrahedral },
            StereoAtomDelta::Mirror { id: StereoAtomId(0), kind: StereoKind::Tetrahedral },
        ],
        vec![],
    )]
    // Priority Mirror > Swap: for a chiral kind swap and mirror are the same improper generator,
    // so a lone `Swap` normalizes to `Mirror`.
    #[case::swap_normalizes_to_mirror(
        vec![StereoAtomDelta::Swap { id: StereoAtomId(0), kind: StereoKind::Tetrahedral }],
        vec![StereoAtomDelta::Mirror { id: StereoAtomId(0), kind: StereoKind::Tetrahedral }],
    )]
    // Created + coset op: the Add absorbs the transform (config Th0 → Th1 under swap).
    #[case::add_then_swap(
        vec![
            StereoAtomDelta::Add {
                id: StereoAtomId(0),
                site: AtomId(0),
                ligands: vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                ast: StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
            },
            StereoAtomDelta::Swap { id: StereoAtomId(0), kind: StereoKind::Tetrahedral },
        ],
        vec![StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            ast: StereoAtomAst::new(StereoKind::Tetrahedral, 1u32),
        }],
    )]
    // coset op then Remove: reverts onto the removed (original) ast — recorded Th1 → Th0.
    #[case::swap_then_remove(
        vec![
            StereoAtomDelta::Swap { id: StereoAtomId(0), kind: StereoKind::Tetrahedral },
            StereoAtomDelta::Remove {
                id: StereoAtomId(0),
                site: AtomId(0),
                ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                ast: StereoAtomAst::new(StereoKind::Tetrahedral, 1u32),
            },
        ],
        vec![StereoAtomDelta::Remove {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            ast: StereoAtomAst::new(StereoKind::Tetrahedral, 0u32),
        }],
    )]
    fn test_deltas_canonicalize_stereo_atom(
        #[case] input: Vec<StereoAtomDelta>,
        #[case] expected: Vec<StereoAtomDelta>,
    ) {
        let canon = Deltas::from_iter(input.into_iter().map(Delta::StereoAtom))
            .canonicalize()
            .unwrap();
        assert_eq!(
            canon,
            Deltas::from_iter(expected.into_iter().map(Delta::StereoAtom)),
        );
    }

    #[rstest]
    fn test_constraint_delta_inverse() {
        let constraint = Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: ValueAst::Lit(0),
        });
        assert_eq!(
            ConstraintDelta::Add(constraint.clone()).inverse(),
            ConstraintDelta::Remove(constraint.clone()),
        );
        assert_eq!(
            ConstraintDelta::Add(constraint.clone()).inverse().inverse(),
            ConstraintDelta::Add(constraint),
        );
    }

    #[rstest]
    #[case::atom(
        Delta::Atom(AtomDelta::Add { id: AtomId(0), ast: AtomAst::from_element(Element::C) }),
        Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::C) })
    )]
    #[case::bond(
        Delta::Bond(BondDelta::Add {
            id: BondId(0),
            atoms: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        }),
        Delta::Bond(BondDelta::Remove {
            id: BondId(0),
            atoms: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        })
    )]
    fn test_delta_inverse(#[case] input: Delta, #[case] expected: Delta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(5), Some(&5))]
    #[case::modified(EntitySpan::Modified { lhs: 1, rhs: 2 }, Some(&1))]
    #[case::removed(EntitySpan::Removed(7), Some(&7))]
    #[case::added(EntitySpan::Added(9), None)]
    fn test_entity_span_lhs(#[case] state: EntitySpan<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(state.lhs(), expected);
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(5), Some(&5))]
    #[case::modified(EntitySpan::Modified { lhs: 1, rhs: 2 }, Some(&2))]
    #[case::added(EntitySpan::Added(9), Some(&9))]
    #[case::removed(EntitySpan::Removed(7), None)]
    fn test_entity_span_rhs(#[case] state: EntitySpan<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(state.rhs(), expected);
    }

    #[rstest]
    #[case::singleton_set(ValueAst::Lit(1), ValueAst::lit_set([1]))]
    fn test_entity_span_superimpose_canonical(#[case] lhs: ValueAst, #[case] rhs: ValueAst) {
        // Canonically-equal sides that are structurally distinct superimpose to `Unchanged`,
        // not `Modified`.
        let lhs = AtomAst::from_element(Element::C).with_charge(lhs);
        let rhs = AtomAst::from_element(Element::C).with_charge(rhs);
        assert_eq!(
            EntitySpan::superimpose(Some(lhs.clone()), Some(rhs)),
            Some(EntitySpan::Unchanged(lhs))
        );
    }

    #[fixture]
    fn remapping() -> IdRemapping {
        IdRemapping::new(
            HashMap::from([
                (AtomId(0), AtomId(2)),
                (AtomId(1), AtomId(0)),
                (AtomId(2), AtomId(1)),
            ]),
            HashMap::from([(BondId(0), BondId(1)), (BondId(1), BondId(0))]),
            HashMap::from([(DativeBondId(0), DativeBondId(1))]),
            HashMap::from([(AromaticSystemId(0), AromaticSystemId(1))]),
            HashMap::from([(MulticenterBondId(0), MulticenterBondId(1))]),
            HashMap::from([(NoncovalentBondId(0), NoncovalentBondId(1))]),
            HashMap::new(),
            HashMap::new(),
        )
    }

    #[rstest]
    #[case::atom(
        Delta::Atom(AtomDelta::Add { id: AtomId(1), ast: AtomAst::from_element(Element::C) }),
        Delta::Atom(AtomDelta::Add { id: AtomId(0), ast: AtomAst::from_element(Element::C) })
    )]
    #[case::bond(
        Delta::Bond(BondDelta::Add {
            id: BondId(0),
            atoms: [AtomId(2), AtomId(1)],
            ast: BondAst::default(),
        }),
        Delta::Bond(BondDelta::Add {
            id: BondId(1),
            atoms: [AtomId(1), AtomId(0)],
            ast: BondAst::default(),
        })
    )]
    #[case::dative_resort(
        Delta::DativeBond(DativeBondDelta::Add {
            id: DativeBondId(0),
            donors: vec![AtomId(0), AtomId(2)],
            acceptor: AtomId(1),
            ast: DativeBondAst::from_order(1),
        }),
        Delta::DativeBond(DativeBondDelta::Add {
            id: DativeBondId(1),
            donors: vec![AtomId(1), AtomId(2)],
            acceptor: AtomId(0),
            ast: DativeBondAst::from_order(1),
        })
    )]
    #[case::aromatic_resort_permute(
        Delta::AromaticSystem(AromaticSystemDelta::Add {
            id: AromaticSystemId(0),
            atoms: vec![AtomId(0), AtomId(1)],
            ast: AromaticSystemAst::from_electrons(vec![1, 2]),
        }),
        Delta::AromaticSystem(AromaticSystemDelta::Add {
            id: AromaticSystemId(1),
            atoms: vec![AtomId(0), AtomId(2)],
            ast: AromaticSystemAst::from_electrons(vec![2, 1]),
        })
    )]
    #[case::aromatic_remove(
        Delta::AromaticSystem(AromaticSystemDelta::Remove {
            id: AromaticSystemId(0),
            atoms: vec![AtomId(0), AtomId(1)],
            ast: AromaticSystemAst::from_electrons(vec![1, 2]),
        }),
        Delta::AromaticSystem(AromaticSystemDelta::Remove {
            id: AromaticSystemId(1),
            atoms: vec![AtomId(0), AtomId(2)],
            ast: AromaticSystemAst::from_electrons(vec![2, 1]),
        })
    )]
    #[case::multicenter_resort_permute(
        Delta::MulticenterBond(MulticenterBondDelta::Add {
            id: MulticenterBondId(0),
            atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
            ast: MulticenterBondAst::from_electrons(vec![3, 5, 7]),
        }),
        Delta::MulticenterBond(MulticenterBondDelta::Add {
            id: MulticenterBondId(1),
            atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
            ast: MulticenterBondAst::from_electrons(vec![5, 7, 3]),
        })
    )]
    #[case::noncovalent_resort(
        Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(0),
            atoms: [AtomId(2), AtomId(1)],
            ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        }),
        Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(1),
            atoms: [AtomId(0), AtomId(1)],
            ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        })
    )]
    #[case::overlay_modify_field(
        Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(0),
            change: AromaticSystemFieldChange::Charge {
                old: ValueAst::Lit(0),
                new: ValueAst::Lit(1),
            },
        }),
        Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(1),
            change: AromaticSystemFieldChange::Charge {
                old: ValueAst::Lit(0),
                new: ValueAst::Lit(1),
            },
        })
    )]
    #[case::constraint(
        Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(
            MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) },
        ))),
        Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(
            MoleculeConstraint::ChargeSum { atoms: None, sum: ValueAst::Lit(0) },
        )))
    )]
    fn test_remap_delta(remapping: IdRemapping, #[case] input: Delta, #[case] expected: Delta) {
        assert_eq!(remap_delta(input, &remapping), expected);
    }

    fn charge_set(id: u32, old: i64, new: i64) -> Delta {
        Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(id),
            change: AtomFieldChange::Charge {
                old: ValueAst::Lit(old),
                new: ValueAst::Lit(new),
            },
        })
    }

    #[rstest]
    fn test_deltas_canonicalize_field_fusion() {
        let deltas = Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 1, 2)]);
        assert_eq!(
            deltas.canonicalize().unwrap(),
            Deltas::from_iter([charge_set(0, 0, 2)]),
        );
    }

    #[rstest]
    fn test_deltas_canonicalize_field_noop_dropped() {
        let deltas = Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 1, 0)]);
        assert_eq!(deltas.canonicalize().unwrap(), Deltas::new());
    }

    #[rstest]
    fn test_deltas_canonicalize_created_absorbs_field() {
        let deltas = Deltas::from_iter([
            Delta::Atom(AtomDelta::Add {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C).with_charge(ValueAst::Lit(0)),
            }),
            charge_set(0, 0, 1),
        ]);
        assert_eq!(
            deltas.canonicalize().unwrap(),
            Deltas::from_iter([Delta::Atom(AtomDelta::Add {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C).with_charge(ValueAst::Lit(1)),
            })]),
        );
    }

    #[rstest]
    fn test_deltas_canonicalize_created_then_removed_cancels() {
        let deltas = Deltas::from_iter([
            Delta::Atom(AtomDelta::Add {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            }),
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            }),
        ]);
        assert_eq!(deltas.canonicalize().unwrap(), Deltas::new());
    }

    #[rstest]
    fn test_deltas_canonicalize_remove_subsumes_field() {
        // ModifyField then Remove must canonicalize to a Remove carrying the original value.
        let deltas = Deltas::from_iter([
            charge_set(0, 0, 1),
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C).with_charge(ValueAst::Lit(1)),
            }),
        ]);
        assert_eq!(
            deltas.canonicalize().unwrap(),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C).with_charge(ValueAst::Lit(0)),
            })]),
        );
    }

    #[rstest]
    fn test_deltas_canonicalize_constraint_chain() {
        let deltas = Deltas::from_iter([
            Delta::Atom(AtomDelta::ModifyConstraint {
                id: AtomId(0),
                old: None,
                new: Some(AtomConstraintAst::Valence(ValueAst::Lit(4))),
            }),
            Delta::Atom(AtomDelta::ModifyConstraint {
                id: AtomId(0),
                old: Some(AtomConstraintAst::Valence(ValueAst::Lit(4))),
                new: Some(AtomConstraintAst::Valence(ValueAst::Lit(3))),
            }),
        ]);
        assert_eq!(
            deltas.canonicalize().unwrap(),
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyConstraint {
                id: AtomId(0),
                old: None,
                new: Some(AtomConstraintAst::Valence(ValueAst::Lit(3))),
            })]),
        );
    }

    #[rstest]
    fn test_deltas_canonicalize_order_independent() {
        let order_set = Delta::Bond(BondDelta::ModifyField {
            id: BondId(0),
            change: BondFieldChange::Order {
                old: ValueAst::Lit(1),
                new: ValueAst::Lit(2),
            },
        });
        let forward = Deltas::from_iter([charge_set(0, 0, 1), order_set.clone()]);
        let reverse = Deltas::from_iter([order_set, charge_set(0, 0, 1)]);
        assert_eq!(
            forward.canonicalize().unwrap(),
            reverse.canonicalize().unwrap()
        );
    }

    #[rstest]
    fn test_deltas_canonicalize_idempotent() {
        let once = Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 1, 2)])
            .canonicalize()
            .unwrap();
        assert_eq!(once.clone().canonicalize().unwrap(), once);
    }

    #[rstest]
    fn test_deltas_canonicalize_dangling_bond_error() {
        let deltas = Deltas::from_iter([
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            }),
            Delta::Bond(BondDelta::Add {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: BondAst::default(),
            }),
        ]);
        assert!(matches!(deltas.canonicalize(), Err(Contradiction)));
    }

    #[rstest]
    fn test_deltas_canonicalize_discontinuous_chain_error() {
        let deltas = Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 2, 3)]);
        assert!(matches!(deltas.canonicalize(), Err(Contradiction)));
    }

    fn charge_sum(sum: i64) -> Constraint {
        Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: ValueAst::Lit(sum),
        })
    }

    #[rstest]
    fn test_deltas_canonicalize_molecule_constraint_cancels() {
        let deltas = Deltas::from_iter([
            Delta::Constraint(ConstraintDelta::Add(charge_sum(0))),
            Delta::Constraint(ConstraintDelta::Remove(charge_sum(0))),
        ]);
        assert_eq!(deltas.canonicalize().unwrap(), Deltas::new());
    }

    #[rstest]
    fn test_deltas_canonicalize_molecule_constraint_multiplicity() {
        // Two adds and one remove net to one add — multiset, not dedup.
        let deltas = Deltas::from_iter([
            Delta::Constraint(ConstraintDelta::Add(charge_sum(0))),
            Delta::Constraint(ConstraintDelta::Add(charge_sum(0))),
            Delta::Constraint(ConstraintDelta::Remove(charge_sum(0))),
        ]);
        assert_eq!(
            deltas.canonicalize().unwrap(),
            Deltas::from_iter([Delta::Constraint(ConstraintDelta::Add(charge_sum(0)))]),
        );
    }
}
