//! Resolved edit vocabulary: the `Delta` counterpart of the deferred `Edit`.
//!
//! A `Delta` is one resolved edit over a `MoleculeAst`, referencing entities by stable
//! ids in the molecule's own frame (no positional `New`). The vocabulary is closed
//! under inversion — every delta's inverse is another delta.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::Hash;
use std::mem::{discriminant, Discriminant};
use std::slice::Iter;

use super::atom::AtomAst;
use super::bond::BondAst;
use super::constraint::{
    AtomConstraint, AtomConstraintKey, BondConstraint, BondConstraintKey, Constraint,
};
use super::edit::{AtomFieldChange, BondFieldChange};
use super::error::Contradiction;
use super::id::{AtomId, BondId};
use super::traits::Canonicalize;

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
    SetField {
        id: AtomId,
        change: AtomFieldChange,
    },
    SetConstraint {
        id: AtomId,
        old: Option<AtomConstraint>,
        new: Option<AtomConstraint>,
    },
}

impl AtomDelta {
    /// The inverse delta: `Add`↔`Remove`; `SetField` / `SetConstraint` swap old/new.
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, ast } => Self::Remove { id, ast },
            Self::Remove { id, ast } => Self::Add { id, ast },
            Self::SetField { id, change } => Self::SetField {
                id,
                change: change.inverse(),
            },
            Self::SetConstraint { id, old, new } => Self::SetConstraint {
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
        endpoints: [AtomId; 2],
        ast: BondAst,
    },
    Remove {
        id: BondId,
        endpoints: [AtomId; 2],
        ast: BondAst,
    },
    SetField {
        id: BondId,
        change: BondFieldChange,
    },
    SetConstraint {
        id: BondId,
        old: Option<BondConstraint>,
        new: Option<BondConstraint>,
    },
}

impl BondDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, endpoints, ast } => Self::Remove { id, endpoints, ast },
            Self::Remove { id, endpoints, ast } => Self::Add { id, endpoints, ast },
            Self::SetField { id, change } => Self::SetField {
                id,
                change: change.inverse(),
            },
            Self::SetConstraint { id, old, new } => Self::SetConstraint {
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
    Constraint(ConstraintDelta),
}

impl Delta {
    /// The inverse delta.
    pub fn inverse(self) -> Self {
        match self {
            Self::Atom(delta) => Self::Atom(delta.inverse()),
            Self::Bond(delta) => Self::Bond(delta.inverse()),
            Self::Constraint(delta) => Self::Constraint(delta.inverse()),
        }
    }
}

/// Generates the per-variant field-change operations (`fuse_field`, `field_is_identity`,
/// `apply_field`) of a `DeltaFamily` impl from the `(variant => ast field)` map.
macro_rules! field_ops {
    ($change:ident, $ast:ident, { $($variant:ident => $field:ident),+ $(,)? }) => {
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
    };
}

/// The per-entity op the fold operates on, abstracting `AtomDelta`/`BondDelta`. `payload`
/// is the structural data of `Add`/`Remove` (`()` for atoms, endpoints for bonds).
enum EntityOp<F: DeltaFamily> {
    Add {
        payload: F::Payload,
        ast: F::Ast,
    },
    Remove {
        payload: F::Payload,
        ast: F::Ast,
    },
    SetField(F::FieldChange),
    SetConstraint {
        old: Option<F::Constraint>,
        new: Option<F::Constraint>,
    },
}

/// A per-entity delta family the canonicalize fold is generic over. Atoms and bonds (and,
/// later, the overlay families) supply the structural deconstruction (`id`/`split`/
/// `rebuild`) and per-variant field/constraint operations; the fold itself is written once.
trait DeltaFamily: Sized {
    type Id: Copy + Eq + Hash;
    type Ast;
    type FieldChange;
    type Constraint: PartialEq;
    type ConstraintKey: Eq + Hash;
    type Payload;

    fn id(&self) -> Self::Id;
    fn split(self) -> EntityOp<Self>;
    fn rebuild(id: Self::Id, op: EntityOp<Self>) -> Self;

    fn fuse_field(prev: Self::FieldChange, next: Self::FieldChange) -> Option<Self::FieldChange>;
    fn field_is_identity(change: &Self::FieldChange) -> bool;
    fn field_inverse(change: Self::FieldChange) -> Self::FieldChange;
    fn apply_field(ast: &mut Self::Ast, change: Self::FieldChange) -> Result<(), Contradiction>;

    fn constraint_key(constraint: &Self::Constraint) -> Self::ConstraintKey;
    fn apply_constraint(
        ast: &mut Self::Ast,
        old: Option<Self::Constraint>,
        new: Option<Self::Constraint>,
    ) -> Result<(), Contradiction>;
}

/// Fold one entity's ops (input order) to its normal form. `created` (an `Add` is present)
/// vs `preserved` paths per doc 131.
fn fold_group<F: DeltaFamily>(id: F::Id, group: Vec<F>) -> Result<Vec<F>, Contradiction> {
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
fn fold_created<F: DeltaFamily>(ops: Vec<EntityOp<F>>) -> Result<Vec<EntityOp<F>>, Contradiction> {
    let mut state: Option<(F::Payload, F::Ast)> = None;
    let mut removed = false;
    for op in ops {
        if removed {
            return Err(Contradiction);
        }
        match op {
            EntityOp::Add { payload, ast } => {
                if state.is_some() {
                    return Err(Contradiction);
                }
                state = Some((payload, ast));
            }
            EntityOp::SetField(change) => {
                let (_, ast) = state.as_mut().ok_or(Contradiction)?;
                F::apply_field(ast, change)?;
            }
            EntityOp::SetConstraint { old, new } => {
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
        Some((payload, ast)) => vec![EntityOp::Add { payload, ast }],
        None => Vec::new(),
    })
}

/// Preserved entity: fuse `SetField` chains per field and `SetConstraint` chains per key. A
/// `Remove` subsumes the prior changes and carries the *original* value (the changes are
/// reverted on the removed ast).
#[allow(clippy::type_complexity)]
fn fold_preserved<F: DeltaFamily>(
    ops: Vec<EntityOp<F>>,
) -> Result<Vec<EntityOp<F>>, Contradiction> {
    let mut fields: HashMap<Discriminant<F::FieldChange>, F::FieldChange> = HashMap::new();
    let mut constraints: HashMap<F::ConstraintKey, (Option<F::Constraint>, Option<F::Constraint>)> =
        HashMap::new();
    let mut removed: Option<(F::Payload, F::Ast)> = None;
    for op in ops {
        if removed.is_some() {
            return Err(Contradiction);
        }
        match op {
            EntityOp::Add { .. } => return Err(Contradiction),
            EntityOp::SetField(change) => {
                let slot = discriminant(&change);
                let fused = match fields.remove(&slot) {
                    Some(prev) => F::fuse_field(prev, change).ok_or(Contradiction)?,
                    None => change,
                };
                fields.insert(slot, fused);
            }
            EntityOp::SetConstraint { old, new } => {
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
            EntityOp::Remove { payload, ast } => {
                removed = Some((payload, ast));
            }
        }
    }
    if let Some((payload, mut ast)) = removed {
        for (_slot, change) in fields {
            F::apply_field(&mut ast, F::field_inverse(change))?;
        }
        for (_key, (old, new)) in constraints {
            F::apply_constraint(&mut ast, new, old)?;
        }
        return Ok(vec![EntityOp::Remove { payload, ast }]);
    }
    let mut out = Vec::new();
    for (_slot, change) in fields {
        if !F::field_is_identity(&change) {
            out.push(EntityOp::SetField(change));
        }
    }
    for (_key, (old, new)) in constraints {
        if old != new {
            out.push(EntityOp::SetConstraint { old, new });
        }
    }
    Ok(out)
}

impl DeltaFamily for AtomDelta {
    type Id = AtomId;
    type Ast = AtomAst;
    type FieldChange = AtomFieldChange;
    type Constraint = AtomConstraint;
    type ConstraintKey = AtomConstraintKey;
    type Payload = ();

    fn id(&self) -> AtomId {
        match self {
            AtomDelta::Add { id, .. }
            | AtomDelta::Remove { id, .. }
            | AtomDelta::SetField { id, .. }
            | AtomDelta::SetConstraint { id, .. } => *id,
        }
    }

    fn split(self) -> EntityOp<Self> {
        match self {
            AtomDelta::Add { ast, .. } => EntityOp::Add { payload: (), ast },
            AtomDelta::Remove { ast, .. } => EntityOp::Remove { payload: (), ast },
            AtomDelta::SetField { change, .. } => EntityOp::SetField(change),
            AtomDelta::SetConstraint { old, new, .. } => EntityOp::SetConstraint { old, new },
        }
    }

    fn rebuild(id: AtomId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { ast, .. } => AtomDelta::Add { id, ast },
            EntityOp::Remove { ast, .. } => AtomDelta::Remove { id, ast },
            EntityOp::SetField(change) => AtomDelta::SetField { id, change },
            EntityOp::SetConstraint { old, new } => AtomDelta::SetConstraint { id, old, new },
        }
    }

    field_ops!(AtomFieldChange, AtomAst, {
        Element => element,
        IsotopeMass => isotope_mass,
        Charge => charge,
        ImplicitHydrogens => implicit_hydrogens,
        LonePairs => lone_pairs,
        Spin => spin,
    });

    fn field_inverse(change: AtomFieldChange) -> AtomFieldChange {
        change.inverse()
    }

    fn constraint_key(constraint: &AtomConstraint) -> AtomConstraintKey {
        constraint.key()
    }

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

impl DeltaFamily for BondDelta {
    type Id = BondId;
    type Ast = BondAst;
    type FieldChange = BondFieldChange;
    type Constraint = BondConstraint;
    type ConstraintKey = BondConstraintKey;
    type Payload = [AtomId; 2];

    fn id(&self) -> BondId {
        match self {
            BondDelta::Add { id, .. }
            | BondDelta::Remove { id, .. }
            | BondDelta::SetField { id, .. }
            | BondDelta::SetConstraint { id, .. } => *id,
        }
    }

    fn split(self) -> EntityOp<Self> {
        match self {
            BondDelta::Add { endpoints, ast, .. } => EntityOp::Add {
                payload: endpoints,
                ast,
            },
            BondDelta::Remove { endpoints, ast, .. } => EntityOp::Remove {
                payload: endpoints,
                ast,
            },
            BondDelta::SetField { change, .. } => EntityOp::SetField(change),
            BondDelta::SetConstraint { old, new, .. } => EntityOp::SetConstraint { old, new },
        }
    }

    fn rebuild(id: BondId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { payload, ast } => BondDelta::Add {
                id,
                endpoints: payload,
                ast,
            },
            EntityOp::Remove { payload, ast } => BondDelta::Remove {
                id,
                endpoints: payload,
                ast,
            },
            EntityOp::SetField(change) => BondDelta::SetField { id, change },
            EntityOp::SetConstraint { old, new } => BondDelta::SetConstraint { id, old, new },
        }
    }

    field_ops!(BondFieldChange, BondAst, {
        Order => order,
        Charge => charge,
        Spin => spin,
    });

    fn field_inverse(change: BondFieldChange) -> BondFieldChange {
        change.inverse()
    }

    fn constraint_key(constraint: &BondConstraint) -> BondConstraintKey {
        constraint.key()
    }

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

/// Apply a resolved per-entity change to a value AST, reusing the `DeltaFamily` apply that
/// `canonicalize` uses. `SetField` / `SetConstraint` mutate the ast; `Add` / `Remove` are
/// no-ops (they carry a whole ast, not a change). Materializes the right-hand value of a
/// preserved entity for a `ReactionSpanAst`.
pub(crate) fn apply_atom_change(ast: &mut AtomAst, delta: &AtomDelta) -> Result<(), Contradiction> {
    match delta {
        AtomDelta::SetField { change, .. } => {
            <AtomDelta as DeltaFamily>::apply_field(ast, change.clone())
        }
        AtomDelta::SetConstraint { old, new, .. } => {
            <AtomDelta as DeltaFamily>::apply_constraint(ast, old.clone(), new.clone())
        }
        AtomDelta::Add { .. } | AtomDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_bond_change(ast: &mut BondAst, delta: &BondDelta) -> Result<(), Contradiction> {
    match delta {
        BondDelta::SetField { change, .. } => {
            <BondDelta as DeltaFamily>::apply_field(ast, change.clone())
        }
        BondDelta::SetConstraint { old, new, .. } => {
            <BondDelta as DeltaFamily>::apply_constraint(ast, old.clone(), new.clone())
        }
        BondDelta::Add { .. } | BondDelta::Remove { .. } => Ok(()),
    }
}

/// Re-anchor a delta's ids and bond endpoints through total atom/bond id maps. Used to move
/// deltas between frames (reverse re-anchoring, composition). The maps must cover every id the
/// delta references.
pub(crate) fn remap_delta(
    delta: Delta,
    atom: &HashMap<AtomId, AtomId>,
    bond: &HashMap<BondId, BondId>,
) -> Delta {
    match delta {
        Delta::Atom(a) => Delta::Atom(match a {
            AtomDelta::Add { id, ast } => AtomDelta::Add {
                id: atom[&id],
                ast,
            },
            AtomDelta::Remove { id, ast } => AtomDelta::Remove {
                id: atom[&id],
                ast,
            },
            AtomDelta::SetField { id, change } => AtomDelta::SetField {
                id: atom[&id],
                change,
            },
            AtomDelta::SetConstraint { id, old, new } => AtomDelta::SetConstraint {
                id: atom[&id],
                old,
                new,
            },
        }),
        Delta::Bond(b) => Delta::Bond(match b {
            BondDelta::Add { id, endpoints, ast } => BondDelta::Add {
                id: bond[&id],
                endpoints: [atom[&endpoints[0]], atom[&endpoints[1]]],
                ast,
            },
            BondDelta::Remove { id, endpoints, ast } => BondDelta::Remove {
                id: bond[&id],
                endpoints: [atom[&endpoints[0]], atom[&endpoints[1]]],
                ast,
            },
            BondDelta::SetField { id, change } => BondDelta::SetField {
                id: bond[&id],
                change,
            },
            BondDelta::SetConstraint { id, old, new } => BondDelta::SetConstraint {
                id: bond[&id],
                old,
                new,
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
        let mut constraints: Vec<ConstraintDelta> = Vec::new();
        for delta in self.0 {
            match delta {
                Delta::Atom(d) => atoms.entry(d.id()).or_default().push(d),
                Delta::Bond(d) => bonds.entry(d.id()).or_default().push(d),
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
                if let BondDelta::Add { endpoints, .. } = delta {
                    if endpoints.iter().any(|atom| removed_atoms.contains(atom)) {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::Bond));
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
    use super::super::value::ValueAst;
    use super::*;

    #[rstest]
    #[case::add_remove(
        AtomDelta::Add { id: AtomId(0), ast: AtomAst::from_element(Element::C) },
        AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::C) }
    )]
    #[case::set_field(
        AtomDelta::SetField {
            id: AtomId(1),
            change: AtomFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Lit(1) },
        },
        AtomDelta::SetField {
            id: AtomId(1),
            change: AtomFieldChange::Charge { old: ValueAst::Lit(1), new: ValueAst::Lit(0) },
        }
    )]
    #[case::set_constraint(
        AtomDelta::SetConstraint {
            id: AtomId(2),
            old: Some(AtomConstraint::Valence(ValueAst::Lit(4))),
            new: Some(AtomConstraint::Valence(ValueAst::Lit(3))),
        },
        AtomDelta::SetConstraint {
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
            endpoints: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        },
        BondDelta::Remove {
            id: BondId(0),
            endpoints: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        }
    )]
    #[case::set_field(
        BondDelta::SetField {
            id: BondId(2),
            change: BondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
        },
        BondDelta::SetField {
            id: BondId(2),
            change: BondFieldChange::Order { old: ValueAst::Lit(2), new: ValueAst::Lit(1) },
        }
    )]
    #[case::set_constraint(
        BondDelta::SetConstraint {
            id: BondId(3),
            old: None,
            new: Some(BondConstraint::Aromatic),
        },
        BondDelta::SetConstraint {
            id: BondId(3),
            old: Some(BondConstraint::Aromatic),
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
            endpoints: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        }),
        Delta::Bond(BondDelta::Remove {
            id: BondId(0),
            endpoints: [AtomId(0), AtomId(1)],
            ast: BondAst::default(),
        })
    )]
    fn test_delta_inverse(#[case] input: Delta, #[case] expected: Delta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    fn charge_set(id: u32, old: i64, new: i64) -> Delta {
        Delta::Atom(AtomDelta::SetField {
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
        // SetField then Remove must canonicalize to a Remove carrying the original value.
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
            Delta::Atom(AtomDelta::SetConstraint {
                id: AtomId(0),
                old: None,
                new: Some(AtomConstraint::Valence(ValueAst::Lit(4))),
            }),
            Delta::Atom(AtomDelta::SetConstraint {
                id: AtomId(0),
                old: Some(AtomConstraint::Valence(ValueAst::Lit(4))),
                new: Some(AtomConstraint::Valence(ValueAst::Lit(3))),
            }),
        ]);
        assert_eq!(
            deltas.canonicalize().unwrap(),
            Deltas::from_iter([Delta::Atom(AtomDelta::SetConstraint {
                id: AtomId(0),
                old: None,
                new: Some(AtomConstraint::Valence(ValueAst::Lit(3))),
            })]),
        );
    }

    #[rstest]
    fn test_deltas_canonicalize_order_independent() {
        let order_set = Delta::Bond(BondDelta::SetField {
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
                endpoints: [AtomId(0), AtomId(1)],
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
