//! Resolved edit vocabulary: the `Delta` counterpart of the deferred `Edit`.
//!
//! A `Delta` is one resolved edit over a `MoleculeAst`, referencing entities by stable
//! ids in the molecule's own id space (no positional `New`). The vocabulary is closed
//! under inversion — every delta's inverse is another delta.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::Hash;
use std::mem::{discriminant, Discriminant};
use std::slice::{Iter, IterMut};

use umol_graph_core::{FactorOrdering, Unordered};

use super::aromatic::AromaticSystemAst;
use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::{
    AromaticSystemConstraint, AromaticSystemConstraintKey, AtomConstraint, AtomConstraintKey,
    BondConstraint, BondConstraintKey, Constraint, DativeBondConstraint, DativeBondConstraintKey,
    MulticenterBondConstraint, MulticenterBondConstraintKey, NoncovalentBondConstraint,
    NoncovalentBondConstraintKey,
};
use super::dative::DativeBondAst;
use super::edit::{
    AromaticSystemFieldChange, AtomFieldChange, BondFieldChange, DativeBondFieldChange,
    MulticenterBondFieldChange, NoncovalentBondFieldChange,
};
use super::error::Contradiction;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::NoncovalentBondAst;
use super::remap::IdRemapping;
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
        old: Option<AtomConstraint>,
        new: Option<AtomConstraint>,
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
        old: Option<BondConstraint>,
        new: Option<BondConstraint>,
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
        old: Option<DativeBondConstraint>,
        new: Option<DativeBondConstraint>,
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
        old: Option<AromaticSystemConstraint>,
        new: Option<AromaticSystemConstraint>,
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
        old: Option<MulticenterBondConstraint>,
        new: Option<MulticenterBondConstraint>,
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
        old: Option<NoncovalentBondConstraint>,
        new: Option<NoncovalentBondConstraint>,
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
                        if ast.$field != old {
                            return Err(Contradiction);
                        }
                        ast.$field = new;
                    }
                )+
            }
            Ok(())
        }

        fn diff_field(left: &$ast, right: &$ast) -> Vec<$change> {
            let mut out = Vec::new();
            $(
                if left.$field != right.$field {
                    out.push($change::$variant {
                        old: left.$field.clone(),
                        new: right.$field.clone(),
                    });
                }
            )+
            out
        }

        #[allow(clippy::type_complexity)]
        fn diff_constraints(
            left: &$ast,
            right: &$ast,
        ) -> Vec<(Option<$constraint>, Option<$constraint>)> {
            let mut left_by_key: HashMap<_, $constraint> = HashMap::new();
            for constraint in left.constraints.iter() {
                left_by_key.insert(constraint.key(), constraint.clone());
            }
            let mut right_by_key: HashMap<_, $constraint> = HashMap::new();
            for constraint in right.constraints.iter() {
                right_by_key.insert(constraint.key(), constraint.clone());
            }
            let mut keys: HashSet<_> = left_by_key.keys().cloned().collect();
            keys.extend(right_by_key.keys().cloned());
            let mut out = Vec::new();
            for key in keys {
                let l = left_by_key.get(&key).cloned();
                let r = right_by_key.get(&key).cloned();
                if l != r {
                    out.push((l, r));
                }
            }
            out
        }
    };
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
                    ) if prev_new == next_old => Some($change::$variant { old, new }),
                )+
                #[allow(unreachable_patterns)]
                _ => None,
            }
        }

        fn field_is_identity(change: &$change) -> bool {
            match change {
                $( $change::$variant { old, new } => old == new, )+
            }
        }
    };
}

/// One entity's span across a reaction — its slice of the superimposed `L`∪`K`∪`R`. A *state*, not
/// an operation (unlike `Edit` / `Delta`). `left()` / `right()` read the side values.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EntitySpan<T> {
    /// In the interface `K` — present and identical on both sides.
    Unchanged(T),
    /// In the interface `K` — present on both sides but relabeled (a dynamic entity).
    Modified { left: T, right: T },
    /// In `R` only — created.
    Added(T),
    /// In `L` only — deleted.
    Removed(T),
}

impl<T> EntitySpan<T> {
    /// The left-side (`L`) value, or `None` if the entity is created.
    pub fn left(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Removed(value) | Self::Modified { left: value, .. } => {
                Some(value)
            }
            Self::Added(_) => None,
        }
    }

    /// The right-side (`R`) value, or `None` if the entity is deleted.
    pub fn right(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Added(value) | Self::Modified { right: value, .. } => {
                Some(value)
            }
            Self::Removed(_) => None,
        }
    }
}

/// A molecule-level constraint's span across a reaction — its slice of the superimposed `L`∪`K`∪`R`.
/// A *state*, not an operation (unlike `ConstraintDelta`). `left()` / `right()` read the side values.
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
    /// The left-side (`L`) value, or `None` if the constraint is created.
    pub fn left(&self) -> Option<&Constraint> {
        match self {
            Self::Unchanged(value) | Self::Removed(value) => Some(value),
            Self::Added(_) => None,
        }
    }

    /// The right-side (`R`) value, or `None` if the constraint is deleted.
    pub fn right(&self) -> Option<&Constraint> {
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

    /// Recover this kind's deltas from its left/right state column: `Added`/`Removed` become
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
                EntitySpan::Modified { left, right } => {
                    out.extend(
                        Self::diff(id, left, right)
                            .into_iter()
                            .map(Self::into_delta),
                    );
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
                        if prev_new != old {
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
        if old != new {
            out.push(EntityOp::ModifyConstraint { old, new });
        }
    }
    Ok(out)
}

impl EntityPatch for AtomDelta {
    type Id = AtomId;
    type Ast = AtomAst;
    type FieldChange = AtomFieldChange;
    type Constraint = AtomConstraint;

    fn modify_field(id: AtomId, change: AtomFieldChange) -> Self {
        AtomDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: AtomId,
        old: Option<AtomConstraint>,
        new: Option<AtomConstraint>,
    ) -> Self {
        AtomDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(AtomFieldChange, AtomAst, AtomConstraint, {
        Element => element,
        IsotopeMass => isotope_mass,
        Charge => charge,
        ImplicitHydrogens => implicit_hydrogens,
        LonePairs => lone_pairs,
        Spin => spin,
    });

    fn apply_constraint(
        ast: &mut AtomAst,
        old: Option<AtomConstraint>,
        new: Option<AtomConstraint>,
    ) -> Result<(), Contradiction> {
        if let Some(old) = old {
            if ast.constraints.remove_entry(&old).is_none() {
                return Err(Contradiction);
            }
        }
        if let Some(new) = new {
            ast.constraints.add(new);
        }
        Ok(())
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

    fn constraint_key(constraint: &AtomConstraint) -> AtomConstraintKey {
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
    type Constraint = BondConstraint;

    fn modify_field(id: BondId, change: BondFieldChange) -> Self {
        BondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: BondId,
        old: Option<BondConstraint>,
        new: Option<BondConstraint>,
    ) -> Self {
        BondDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(BondFieldChange, BondAst, BondConstraint, {
        Order => order,
        Charge => charge,
        Spin => spin,
    });

    fn apply_constraint(
        ast: &mut BondAst,
        old: Option<BondConstraint>,
        new: Option<BondConstraint>,
    ) -> Result<(), Contradiction> {
        if let Some(old) = old {
            if ast.constraints.remove_entry(&old).is_none() {
                return Err(Contradiction);
            }
        }
        if let Some(new) = new {
            ast.constraints.add(new);
        }
        Ok(())
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

    fn constraint_key(constraint: &BondConstraint) -> BondConstraintKey {
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
    type Constraint = DativeBondConstraint;

    fn modify_field(id: DativeBondId, change: DativeBondFieldChange) -> Self {
        DativeBondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: DativeBondId,
        old: Option<DativeBondConstraint>,
        new: Option<DativeBondConstraint>,
    ) -> Self {
        DativeBondDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(DativeBondFieldChange, DativeBondAst, DativeBondConstraint, {
        Order => order,
    });

    fn apply_constraint(
        ast: &mut DativeBondAst,
        old: Option<DativeBondConstraint>,
        new: Option<DativeBondConstraint>,
    ) -> Result<(), Contradiction> {
        if let Some(old) = old {
            if ast.constraints.remove_entry(&old).is_none() {
                return Err(Contradiction);
            }
        }
        if let Some(new) = new {
            ast.constraints.add(new);
        }
        Ok(())
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

    fn constraint_key(constraint: &DativeBondConstraint) -> DativeBondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(DativeBondFieldChange, { Order });
}

impl EntityPatch for AromaticSystemDelta {
    type Id = AromaticSystemId;
    type Ast = AromaticSystemAst;
    type FieldChange = AromaticSystemFieldChange;
    type Constraint = AromaticSystemConstraint;

    fn modify_field(id: AromaticSystemId, change: AromaticSystemFieldChange) -> Self {
        AromaticSystemDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: AromaticSystemId,
        old: Option<AromaticSystemConstraint>,
        new: Option<AromaticSystemConstraint>,
    ) -> Self {
        AromaticSystemDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(
        AromaticSystemFieldChange,
        AromaticSystemAst,
        AromaticSystemConstraint,
        {
            Electrons => electrons,
            Charge => charge,
            Spin => spin,
        }
    );

    fn apply_constraint(
        ast: &mut AromaticSystemAst,
        old: Option<AromaticSystemConstraint>,
        new: Option<AromaticSystemConstraint>,
    ) -> Result<(), Contradiction> {
        if let Some(old) = old {
            if ast.constraints.remove_by_key(old.key()).is_none() {
                return Err(Contradiction);
            }
        }
        if let Some(new) = new {
            ast.constraints.add(new);
        }
        Ok(())
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

    fn constraint_key(constraint: &AromaticSystemConstraint) -> AromaticSystemConstraintKey {
        constraint.key()
    }

    fold_field_ops!(AromaticSystemFieldChange, { Electrons, Charge, Spin });
}

impl EntityPatch for MulticenterBondDelta {
    type Id = MulticenterBondId;
    type Ast = MulticenterBondAst;
    type FieldChange = MulticenterBondFieldChange;
    type Constraint = MulticenterBondConstraint;

    fn modify_field(id: MulticenterBondId, change: MulticenterBondFieldChange) -> Self {
        MulticenterBondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: MulticenterBondId,
        old: Option<MulticenterBondConstraint>,
        new: Option<MulticenterBondConstraint>,
    ) -> Self {
        MulticenterBondDelta::ModifyConstraint { id, old, new }
    }

    diff_field_ops!(
        MulticenterBondFieldChange,
        MulticenterBondAst,
        MulticenterBondConstraint,
        {
            Electrons => electrons,
            Charge => charge,
            Spin => spin,
        }
    );

    fn apply_constraint(
        ast: &mut MulticenterBondAst,
        old: Option<MulticenterBondConstraint>,
        new: Option<MulticenterBondConstraint>,
    ) -> Result<(), Contradiction> {
        if let Some(old) = old {
            if ast.constraints.remove_by_key(old.key()).is_none() {
                return Err(Contradiction);
            }
        }
        if let Some(new) = new {
            ast.constraints.add(new);
        }
        Ok(())
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

    fn constraint_key(constraint: &MulticenterBondConstraint) -> MulticenterBondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(MulticenterBondFieldChange, { Electrons, Charge, Spin });
}

impl EntityPatch for NoncovalentBondDelta {
    type Id = NoncovalentBondId;
    type Ast = NoncovalentBondAst;
    type FieldChange = NoncovalentBondFieldChange;
    type Constraint = NoncovalentBondConstraint;

    fn modify_field(id: NoncovalentBondId, change: NoncovalentBondFieldChange) -> Self {
        NoncovalentBondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: NoncovalentBondId,
        old: Option<NoncovalentBondConstraint>,
        new: Option<NoncovalentBondConstraint>,
    ) -> Self {
        NoncovalentBondDelta::ModifyConstraint { id, old, new }
    }

    // Hand-written (not `diff_field_ops!`): `NoncovalentBondConstraint` is uninhabited,
    // so the macro's constraint loop would be unreachable code.
    fn apply_field(
        ast: &mut NoncovalentBondAst,
        change: NoncovalentBondFieldChange,
    ) -> Result<(), Contradiction> {
        match change {
            NoncovalentBondFieldChange::Kind { old, new } => {
                if ast.kind != old {
                    return Err(Contradiction);
                }
                ast.kind = new;
            }
        }
        Ok(())
    }

    fn diff_field(
        left: &NoncovalentBondAst,
        right: &NoncovalentBondAst,
    ) -> Vec<NoncovalentBondFieldChange> {
        if left.kind != right.kind {
            vec![NoncovalentBondFieldChange::Kind {
                old: left.kind.clone(),
                new: right.kind.clone(),
            }]
        } else {
            Vec::new()
        }
    }

    fn diff_constraints(
        _left: &NoncovalentBondAst,
        _right: &NoncovalentBondAst,
    ) -> Vec<(
        Option<NoncovalentBondConstraint>,
        Option<NoncovalentBondConstraint>,
    )> {
        Vec::new()
    }

    /// `NoncovalentBondConstraint` is uninhabited, so `old`/`new` are always `None`.
    fn apply_constraint(
        _ast: &mut NoncovalentBondAst,
        old: Option<NoncovalentBondConstraint>,
        new: Option<NoncovalentBondConstraint>,
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

    fn constraint_key(constraint: &NoncovalentBondConstraint) -> NoncovalentBondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(NoncovalentBondFieldChange, { Kind });
}

/// Apply a resolved per-entity change to a value AST, reusing the `EntityPatch` apply that
/// `canonicalize` uses. `ModifyField` / `ModifyConstraint` mutate the ast; `Add` / `Remove` are
/// no-ops (they carry a whole ast, not a change). Materializes the right-hand value of a
/// preserved entity for a `ReactionSpanAst`.
pub(crate) fn apply_atom_change(ast: &mut AtomAst, delta: &AtomDelta) -> Result<(), Contradiction> {
    match delta {
        AtomDelta::ModifyField { change, .. } => {
            <AtomDelta as EntityPatch>::apply_field(ast, change.clone())
        }
        AtomDelta::ModifyConstraint { old, new, .. } => {
            <AtomDelta as EntityPatch>::apply_constraint(ast, old.clone(), new.clone())
        }
        AtomDelta::Add { .. } | AtomDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_bond_change(ast: &mut BondAst, delta: &BondDelta) -> Result<(), Contradiction> {
    match delta {
        BondDelta::ModifyField { change, .. } => {
            <BondDelta as EntityPatch>::apply_field(ast, change.clone())
        }
        BondDelta::ModifyConstraint { old, new, .. } => {
            <BondDelta as EntityPatch>::apply_constraint(ast, old.clone(), new.clone())
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
            <DativeBondDelta as EntityPatch>::apply_field(ast, change.clone())
        }
        DativeBondDelta::ModifyConstraint { old, new, .. } => {
            <DativeBondDelta as EntityPatch>::apply_constraint(ast, old.clone(), new.clone())
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
            <AromaticSystemDelta as EntityPatch>::apply_field(ast, change.clone())
        }
        AromaticSystemDelta::ModifyConstraint { old, new, .. } => {
            <AromaticSystemDelta as EntityPatch>::apply_constraint(ast, old.clone(), new.clone())
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
            <MulticenterBondDelta as EntityPatch>::apply_field(ast, change.clone())
        }
        MulticenterBondDelta::ModifyConstraint { old, new, .. } => {
            <MulticenterBondDelta as EntityPatch>::apply_constraint(ast, old.clone(), new.clone())
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
            <NoncovalentBondDelta as EntityPatch>::apply_field(ast, change.clone())
        }
        NoncovalentBondDelta::ModifyConstraint { old, new, .. } => {
            <NoncovalentBondDelta as EntityPatch>::apply_constraint(ast, old.clone(), new.clone())
        }
        NoncovalentBondDelta::Add { .. } | NoncovalentBondDelta::Remove { .. } => Ok(()),
    }
}

/// Re-anchor a delta's ids and participant atoms through a total id relabeling. Used to move
/// deltas between id spaces (reverse re-anchoring, composition). The relabeling must cover every id
/// the delta references. Overlay participants on `Unordered` factors are re-sorted to canonical
/// order; aromatic/multicenter electrons are permuted to stay aligned with their atoms.
pub(crate) fn remap_delta(delta: Delta, map: &IdRemapping) -> Delta {
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
            DativeBondDelta::Add {
                id,
                donors,
                acceptor,
                ast,
            } => DativeBondDelta::Add {
                id: map.map_dative(id),
                donors: donors.iter().map(|a| map.map_atom(*a)).collect(),
                acceptor: map.map_atom(acceptor),
                ast,
            },
            DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                ast,
            } => DativeBondDelta::Remove {
                id: map.map_dative(id),
                donors: donors.iter().map(|a| map.map_atom(*a)).collect(),
                acceptor: map.map_atom(acceptor),
                ast,
            },
            DativeBondDelta::ModifyField { id, change } => DativeBondDelta::ModifyField {
                id: map.map_dative(id),
                change,
            },
            DativeBondDelta::ModifyConstraint { id, old, new } => DativeBondDelta::ModifyConstraint {
                id: map.map_dative(id),
                old,
                new,
            },
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
            NoncovalentBondDelta::Add { id, atoms, ast } => NoncovalentBondDelta::Add {
                id: map.map_noncovalent(id),
                atoms: [map.map_atom(atoms[0]), map.map_atom(atoms[1])],
                ast,
            },
            NoncovalentBondDelta::Remove { id, atoms, ast } => NoncovalentBondDelta::Remove {
                id: map.map_noncovalent(id),
                atoms: [map.map_atom(atoms[0]), map.map_atom(atoms[1])],
                ast,
            },
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
        let mut constraints: Vec<ConstraintDelta> = Vec::new();
        for delta in self.0 {
            match delta {
                Delta::Atom(d) => atoms.entry(d.id()).or_default().push(d),
                Delta::Bond(d) => bonds.entry(d.id()).or_default().push(d),
                Delta::DativeBond(d) => dative.entry(d.id()).or_default().push(d),
                Delta::AromaticSystem(d) => aromatic.entry(d.id()).or_default().push(d),
                Delta::MulticenterBond(d) => multicenter.entry(d.id()).or_default().push(d),
                Delta::NoncovalentBond(d) => noncovalent.entry(d.id()).or_default().push(d),
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
                        .chain(std::iter::once(acceptor))
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
    use crate::ast::BooleanAst;

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
            old: Some(AtomConstraint::Valence(ValueAst::Lit(4))),
            new: Some(AtomConstraint::Valence(ValueAst::Lit(3))),
        },
        AtomDelta::ModifyConstraint {
            id: AtomId(2),
            old: Some(AtomConstraint::Valence(ValueAst::Lit(3))),
            new: Some(AtomConstraint::Valence(ValueAst::Lit(4))),
        }
    )]
    fn test_atom_delta_inverse(#[case] input: AtomDelta, #[case] expected: AtomDelta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
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
            new: Some(BondConstraint::Aromatic(BooleanAst::Lit(true))),
        },
        BondDelta::ModifyConstraint {
            id: BondId(3),
            old: Some(BondConstraint::Aromatic(BooleanAst::Lit(true))),
            new: None,
        }
    )]
    fn test_bond_delta_inverse(#[case] input: BondDelta, #[case] expected: BondDelta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
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
    #[case::modified(EntitySpan::Modified { left: 1, right: 2 }, Some(&1))]
    #[case::removed(EntitySpan::Removed(7), Some(&7))]
    #[case::added(EntitySpan::Added(9), None)]
    fn test_entity_span_left(#[case] state: EntitySpan<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(state.left(), expected);
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(5), Some(&5))]
    #[case::modified(EntitySpan::Modified { left: 1, right: 2 }, Some(&2))]
    #[case::added(EntitySpan::Added(9), Some(&9))]
    #[case::removed(EntitySpan::Removed(7), None)]
    fn test_entity_span_right(#[case] state: EntitySpan<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(state.right(), expected);
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
    #[case::dative(
        Delta::DativeBond(DativeBondDelta::Add {
            id: DativeBondId(0),
            donors: vec![AtomId(0), AtomId(2)],
            acceptor: AtomId(1),
            ast: DativeBondAst::from_order(1),
        }),
        Delta::DativeBond(DativeBondDelta::Add {
            id: DativeBondId(1),
            donors: vec![AtomId(2), AtomId(1)],
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
    #[case::noncovalent(
        Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(0),
            atoms: [AtomId(2), AtomId(1)],
            ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        }),
        Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(1),
            atoms: [AtomId(1), AtomId(0)],
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
                new: Some(AtomConstraint::Valence(ValueAst::Lit(4))),
            }),
            Delta::Atom(AtomDelta::ModifyConstraint {
                id: AtomId(0),
                old: Some(AtomConstraint::Valence(ValueAst::Lit(4))),
                new: Some(AtomConstraint::Valence(ValueAst::Lit(3))),
            }),
        ]);
        assert_eq!(
            deltas.canonicalize().unwrap(),
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyConstraint {
                id: AtomId(0),
                old: None,
                new: Some(AtomConstraint::Valence(ValueAst::Lit(3))),
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
