//! Molecule-scope constraints, the `Constraint` combinator tree, and the
//! molecule-level `Constraints` store.

use std::collections::HashSet;
use std::mem;
use std::slice::Iter;
use std::vec::IntoIter;

use umol_perm::{DynPermutation, Permutation};

use super::super::edit::{CascadedConstraints, ModifiedConstraint, RemovedConstraint};
use super::super::entity::EntityKind;
use super::super::error::Contradiction;
use super::super::frame::OverlaysFrameAction;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::num::NumForm;
use super::super::remap::{IdRemapping, MoleculeCompaction};
use super::super::spin::UnpairedElectronsForm;
use super::super::stereo::StereoKind;
use super::super::traits::{FrameTransport, Lattice, Normalize};
use super::aromatic::AromaticSystemConstraintForm;
use super::atom::AtomConstraintForm;
use super::bond::BondConstraintForm;
use super::dative::DativeBondConstraintForm;
use super::multicenter::MulticenterBondConstraintForm;
use super::noncovalent::NoncovalentBondConstraintForm;
use super::relational::RelationalConstraint;
use super::stereo::{StereoAtomConstraintForm, StereoBondConstraintForm};

pub(crate) trait ConstraintFrameActions {
    fn dative_bond_action(&self, id: DativeBondId) -> Option<&DynPermutation>;
    fn aromatic_system_action(&self, id: AromaticSystemId) -> Option<&DynPermutation>;
    fn multicenter_bond_action(&self, id: MulticenterBondId) -> Option<&DynPermutation>;
    fn noncovalent_bond_action(&self, id: NoncovalentBondId) -> Option<&DynPermutation>;
    fn stereo_atom_action(&self, id: StereoAtomId) -> Option<&Permutation>;
    fn stereo_bond_action(&self, id: StereoBondId) -> Option<&Permutation>;
}

impl ConstraintFrameActions for OverlaysFrameAction {
    fn dative_bond_action(&self, id: DativeBondId) -> Option<&DynPermutation> {
        self.dative_bonds().action(id)
    }

    fn aromatic_system_action(&self, id: AromaticSystemId) -> Option<&DynPermutation> {
        self.aromatic_systems().action(id)
    }

    fn multicenter_bond_action(&self, id: MulticenterBondId) -> Option<&DynPermutation> {
        self.multicenter_bonds().action(id)
    }

    fn noncovalent_bond_action(&self, id: NoncovalentBondId) -> Option<&DynPermutation> {
        self.noncovalent_bonds().action(id)
    }

    fn stereo_atom_action(&self, id: StereoAtomId) -> Option<&Permutation> {
        self.stereo_atoms().action(id)
    }

    fn stereo_bond_action(&self, id: StereoBondId) -> Option<&Permutation> {
        self.stereo_bonds().action(id)
    }
}

#[derive(Default)]
pub(crate) struct ConstraintFrameActionDomain {
    dative_bonds: HashSet<DativeBondId>,
    aromatic_systems: HashSet<AromaticSystemId>,
    multicenter_bonds: HashSet<MulticenterBondId>,
    noncovalent_bonds: HashSet<NoncovalentBondId>,
    stereo_atoms: HashSet<StereoAtomId>,
    stereo_bonds: HashSet<StereoBondId>,
}

impl ConstraintFrameActionDomain {
    pub(super) fn insert_dative_bond(&mut self, id: DativeBondId) {
        self.dative_bonds.insert(id);
    }

    pub(super) fn insert_aromatic_system(&mut self, id: AromaticSystemId) {
        self.aromatic_systems.insert(id);
    }

    pub(super) fn insert_multicenter_bond(&mut self, id: MulticenterBondId) {
        self.multicenter_bonds.insert(id);
    }

    pub(super) fn insert_noncovalent_bond(&mut self, id: NoncovalentBondId) {
        self.noncovalent_bonds.insert(id);
    }

    pub(super) fn insert_stereo_atom(&mut self, id: StereoAtomId) {
        self.stereo_atoms.insert(id);
    }

    pub(super) fn insert_stereo_bond(&mut self, id: StereoBondId) {
        self.stereo_bonds.insert(id);
    }

    pub(crate) fn contains_dative_bond(&self, id: DativeBondId) -> bool {
        self.dative_bonds.contains(&id)
    }

    pub(crate) fn contains_aromatic_system(&self, id: AromaticSystemId) -> bool {
        self.aromatic_systems.contains(&id)
    }

    pub(crate) fn contains_multicenter_bond(&self, id: MulticenterBondId) -> bool {
        self.multicenter_bonds.contains(&id)
    }

    pub(crate) fn contains_noncovalent_bond(&self, id: NoncovalentBondId) -> bool {
        self.noncovalent_bonds.contains(&id)
    }

    pub(crate) fn contains_stereo_atom(&self, id: StereoAtomId) -> bool {
        self.stereo_atoms.contains(&id)
    }

    pub(crate) fn contains_stereo_bond(&self, id: StereoBondId) -> bool {
        self.stereo_bonds.contains(&id)
    }

    pub(crate) fn count(&self, entity_kind: EntityKind) -> usize {
        match entity_kind {
            EntityKind::Atom | EntityKind::Bond => 0,
            EntityKind::DativeBond => self.dative_bonds.len(),
            EntityKind::AromaticSystem => self.aromatic_systems.len(),
            EntityKind::MulticenterBond => self.multicenter_bonds.len(),
            EntityKind::NoncovalentBond => self.noncovalent_bonds.len(),
            EntityKind::StereoAtom => self.stereo_atoms.len(),
            EntityKind::StereoBond => self.stereo_bonds.len(),
        }
    }
}

/// Tree node type: per-entity leaf, molecule-scope leaf, relational leaf, or
/// combinator. The bare entity-leaf forms appear only inside a combinator
/// (e.g. `And(Atom(..), Bond(..))`) or a molecule-scope predicate;
/// unconditional per-entity value-only constraints live inline on the entity
/// form and are lifted there at DSL → IR conversion time. Cross-entity
/// ref-bearing constraints (e.g. a dative-bond donor identity, aromatic
/// system membership, noncovalent endpoints) live only at molecule scope
/// via `Relational`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Constraint {
    Atom(AtomId, AtomConstraintForm),
    Bond(BondId, BondConstraintForm),
    DativeBond(DativeBondId, DativeBondConstraintForm),
    AromaticSystem(AromaticSystemId, AromaticSystemConstraintForm),
    MulticenterBond(MulticenterBondId, MulticenterBondConstraintForm),
    NoncovalentBond(NoncovalentBondId, NoncovalentBondConstraintForm),
    StereoAtom(StereoAtomId, StereoKind, StereoAtomConstraintForm),
    StereoBond(StereoBondId, StereoKind, StereoBondConstraintForm),
    Relational(RelationalConstraint),
    Molecule(MoleculeConstraint),
    And(Vec<Constraint>),
    Or(Vec<Constraint>),
    Not(Box<Constraint>),
}

impl Constraint {
    /// A constraint is vacuous when it asserts nothing: the entity-leaf
    /// variants delegate to their inner `is_undetermined`; molecule-scope
    /// leaves delegate to `MoleculeConstraint::is_vacuous`. Combinators
    /// (`And`/`Or`) are vacuous only when empty; `Not(c)` is treated as
    /// non-vacuous (negating a vacuous claim is a meaningful unsat claim).
    /// `Relational` is always non-vacuous (no `Undetermined` payload to
    /// elide); the stereo leaves delegate to the inner `is_undetermined` (a
    /// `#o`/`#g` with an `Undetermined` relation constrains nothing).
    pub fn is_vacuous(&self) -> bool {
        match self {
            Self::Atom(_, c) => c.is_undetermined(),
            Self::Bond(_, c) => c.is_undetermined(),
            Self::DativeBond(_, c) => c.is_undetermined(),
            Self::AromaticSystem(_, c) => c.is_undetermined(),
            Self::MulticenterBond(_, c) => c.is_undetermined(),
            Self::NoncovalentBond(_, c) => c.is_undetermined(),
            Self::StereoAtom(_, _, c) => c.is_undetermined(),
            Self::StereoBond(_, _, c) => c.is_undetermined(),
            Self::Relational(_) => false,
            Self::Molecule(c) => c.is_vacuous(),
            Self::And(xs) | Self::Or(xs) => xs.is_empty(),
            Self::Not(_) => false,
        }
    }

    pub fn compact(self, compaction: &MoleculeCompaction) -> Option<Self> {
        match self {
            Constraint::Atom(id, c) => {
                let i = compaction.compact_atom(id)?;
                c.compact(compaction).map(|c| Constraint::Atom(i, c))
            }
            Constraint::Bond(id, c) => {
                let i = compaction.compact_bond(id)?;
                c.compact(compaction).map(|c| Constraint::Bond(i, c))
            }
            Constraint::DativeBond(id, c) => {
                let i = compaction.compact_dative_bond(id)?;
                c.compact(compaction).map(|c| Constraint::DativeBond(i, c))
            }
            Constraint::AromaticSystem(id, c) => {
                let i = compaction.compact_aromatic_system(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::AromaticSystem(i, c))
            }
            Constraint::MulticenterBond(id, c) => {
                let i = compaction.compact_multicenter_bond(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::MulticenterBond(i, c))
            }
            Constraint::NoncovalentBond(id, c) => {
                let i = compaction.compact_noncovalent_bond(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::NoncovalentBond(i, c))
            }
            Constraint::StereoAtom(id, kind, c) => {
                let i = compaction.compact_stereo_atom(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::StereoAtom(i, kind, c))
            }
            Constraint::StereoBond(id, kind, c) => {
                let i = compaction.compact_stereo_bond(id)?;
                c.compact(compaction)
                    .map(|c| Constraint::StereoBond(i, kind, c))
            }
            Constraint::Relational(r) => r.compact(compaction).map(Constraint::Relational),
            Constraint::Molecule(m) => m.compact(compaction).map(Constraint::Molecule),
            Constraint::And(xs) => xs
                .into_iter()
                .map(|c| c.compact(compaction))
                .collect::<Option<Vec<_>>>()
                .map(Constraint::And),
            Constraint::Or(xs) => xs
                .into_iter()
                .map(|c| c.compact(compaction))
                .collect::<Option<Vec<_>>>()
                .map(Constraint::Or),
            Constraint::Not(x) => x.compact(compaction).map(|c| Constraint::Not(Box::new(c))),
        }
    }

    /// Re-anchor every entity ref through a total id remapping (match-based: lhs → host,
    /// created → appended). Total — never drops (the parallel of `compact`, which compacts after
    /// removal). Stays a separate flow from `compact`; the two are not bridged.
    pub fn remap(self, map: &IdRemapping) -> Self {
        match self {
            Constraint::Atom(id, c) => {
                let i = map.map_atom(id);
                Constraint::Atom(i, c.remap(map))
            }
            Constraint::Bond(id, c) => {
                let i = map.map_bond(id);
                Constraint::Bond(i, c.remap(map))
            }
            Constraint::DativeBond(id, c) => {
                Constraint::DativeBond(map.map_dative(id), c.remap(map))
            }
            Constraint::AromaticSystem(id, c) => {
                Constraint::AromaticSystem(map.map_aromatic(id), c.remap(map))
            }
            Constraint::MulticenterBond(id, c) => {
                Constraint::MulticenterBond(map.map_multicenter(id), c.remap(map))
            }
            Constraint::NoncovalentBond(id, c) => {
                Constraint::NoncovalentBond(map.map_noncovalent(id), c.remap(map))
            }
            Constraint::StereoAtom(id, kind, c) => {
                Constraint::StereoAtom(map.map_stereo_atom(id), kind, c.remap(map))
            }
            Constraint::StereoBond(id, kind, c) => {
                Constraint::StereoBond(map.map_stereo_bond(id), kind, c.remap(map))
            }
            Constraint::Relational(r) => Constraint::Relational(r.remap(map)),
            Constraint::Molecule(m) => Constraint::Molecule(m.remap(map)),
            Constraint::And(xs) => Constraint::And(xs.into_iter().map(|c| c.remap(map)).collect()),
            Constraint::Or(xs) => Constraint::Or(xs.into_iter().map(|c| c.remap(map)).collect()),
            Constraint::Not(x) => Constraint::Not(Box::new(x.remap(map))),
        }
    }

    fn collect_frame_action_domain(&self, domain: &mut ConstraintFrameActionDomain) {
        match self {
            Self::DativeBond(id, constraint) if constraint.uses_participant_frame() => {
                domain.insert_dative_bond(*id)
            }
            Self::AromaticSystem(id, constraint) if constraint.uses_participant_frame() => {
                domain.insert_aromatic_system(*id)
            }
            Self::MulticenterBond(id, constraint) if constraint.uses_participant_frame() => {
                domain.insert_multicenter_bond(*id)
            }
            Self::NoncovalentBond(id, constraint) if constraint.uses_participant_frame() => {
                domain.insert_noncovalent_bond(*id)
            }
            Self::StereoAtom(id, _, constraint) if constraint.uses_participant_frame() => {
                domain.insert_stereo_atom(*id)
            }
            Self::StereoBond(id, _, constraint) if constraint.uses_participant_frame() => {
                domain.insert_stereo_bond(*id)
            }
            Self::Relational(constraint) => constraint.collect_frame_action_domain(domain),
            Self::And(constraints) | Self::Or(constraints) => {
                for constraint in constraints {
                    constraint.collect_frame_action_domain(domain);
                }
            }
            Self::Not(constraint) => constraint.collect_frame_action_domain(domain),
            Self::Atom(..)
            | Self::Bond(..)
            | Self::DativeBond(..)
            | Self::AromaticSystem(..)
            | Self::MulticenterBond(..)
            | Self::NoncovalentBond(..)
            | Self::StereoAtom(..)
            | Self::StereoBond(..)
            | Self::Molecule(..) => {}
        }
    }

    pub(crate) fn reframe_by_actions(self, actions: &impl ConstraintFrameActions) -> Option<Self> {
        Some(match self {
            Self::DativeBond(id, constraint) => {
                if constraint.uses_participant_frame() {
                    let action = actions.dative_bond_action(id)?;
                    Self::DativeBond(id, constraint.reframe_by(action)?)
                } else {
                    Self::DativeBond(id, constraint)
                }
            }
            Self::AromaticSystem(id, constraint) => {
                if constraint.uses_participant_frame() {
                    let action = actions.aromatic_system_action(id)?;
                    Self::AromaticSystem(id, constraint.reframe_by(action)?)
                } else {
                    Self::AromaticSystem(id, constraint)
                }
            }
            Self::MulticenterBond(id, constraint) => {
                if constraint.uses_participant_frame() {
                    let action = actions.multicenter_bond_action(id)?;
                    Self::MulticenterBond(id, constraint.reframe_by(action)?)
                } else {
                    Self::MulticenterBond(id, constraint)
                }
            }
            Self::NoncovalentBond(id, constraint) => {
                if constraint.uses_participant_frame() {
                    let action = actions.noncovalent_bond_action(id)?;
                    Self::NoncovalentBond(id, constraint.reframe_by(action)?)
                } else {
                    Self::NoncovalentBond(id, constraint)
                }
            }
            Self::StereoAtom(id, kind, constraint) => {
                if !constraint.uses_participant_frame() {
                    Self::StereoAtom(id, kind, constraint)
                } else {
                    let action = actions.stereo_atom_action(id)?;
                    if !kind.class_key().space().allows(*action) {
                        return None;
                    }
                    Self::StereoAtom(id, kind, constraint.reframe_by(action)?)
                }
            }
            Self::StereoBond(id, kind, constraint) => {
                if !constraint.uses_participant_frame() {
                    Self::StereoBond(id, kind, constraint)
                } else {
                    let action = actions.stereo_bond_action(id)?;
                    if !kind.class_key().space().allows(*action) {
                        return None;
                    }
                    Self::StereoBond(id, kind, constraint.reframe_by(action)?)
                }
            }
            Self::Relational(constraint) => {
                Self::Relational(constraint.reframe_by_actions(actions)?)
            }
            Self::And(constraints) => Self::And(
                constraints
                    .into_iter()
                    .map(|constraint| constraint.reframe_by_actions(actions))
                    .collect::<Option<Vec<_>>>()?,
            ),
            Self::Or(constraints) => Self::Or(
                constraints
                    .into_iter()
                    .map(|constraint| constraint.reframe_by_actions(actions))
                    .collect::<Option<Vec<_>>>()?,
            ),
            Self::Not(constraint) => Self::Not(Box::new(constraint.reframe_by_actions(actions)?)),
            invariant @ (Self::Atom(..) | Self::Bond(..) | Self::Molecule(..)) => invariant,
        })
    }
}

impl FrameTransport for Constraint {
    type Action = OverlaysFrameAction;

    fn reframe_by(self, actions: &Self::Action) -> Option<Self> {
        self.reframe_by_actions(actions)
    }
}

impl Normalize for Constraint {
    /// Normalize the inner predicate of each leaf; for `And`/`Or`, recurse,
    /// flatten the same combinator, drop empty `And`/`Or`, sort + dedup
    /// children by the `Constraint` order, then reduce a trivial wrapper — a
    /// singleton `And`/`Or` is its element. `Not` canonicalizes its inner
    /// node.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Atom(id, c) => Self::Atom(id, c.normalize()?),
            Self::Bond(id, c) => Self::Bond(id, c.normalize()?),
            Self::DativeBond(id, c) => Self::DativeBond(id, c.normalize()?),
            Self::AromaticSystem(id, c) => Self::AromaticSystem(id, c.normalize()?),
            Self::MulticenterBond(id, c) => Self::MulticenterBond(id, c.normalize()?),
            Self::NoncovalentBond(id, c) => Self::NoncovalentBond(id, c.normalize()?),
            Self::StereoAtom(id, kind, c) => Self::StereoAtom(id, kind, c.normalize()?),
            Self::StereoBond(id, kind, c) => Self::StereoBond(id, kind, c.normalize()?),
            Self::Relational(r) => Self::Relational(r.normalize()?),
            Self::Molecule(m) => Self::Molecule(m.normalize()?),
            Self::And(xs) => {
                let mut children = normalize_logical_constraints(xs, true)?;
                if children.len() == 1 {
                    children.remove(0)
                } else {
                    Self::And(children)
                }
            }
            Self::Or(xs) => {
                let mut children = normalize_logical_constraints(xs, false)?;
                if children.len() == 1 {
                    children.remove(0)
                } else {
                    Self::Or(children)
                }
            }
            Self::Not(c) => Self::Not(Box::new((*c).normalize()?)),
        })
    }
}

impl Normalize for Constraints {
    /// The store is an implicit conjunction, so it canonicalizes like an `And`:
    /// flatten top-level `And` entries, drop empty `And`/`Or`, sort + dedup.
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(Self(normalize_logical_constraints(self.0, true)?))
    }
}

impl FrameTransport for Constraints {
    type Action = OverlaysFrameAction;

    fn reframe_by(self, actions: &Self::Action) -> Option<Self> {
        self.reframe_by_actions(actions)
    }
}

impl Constraints {
    pub(crate) fn frame_action_domain(&self) -> ConstraintFrameActionDomain {
        let mut domain = ConstraintFrameActionDomain::default();
        for constraint in &self.0 {
            constraint.collect_frame_action_domain(&mut domain);
        }
        domain
    }

    pub(crate) fn reframe_by_actions(self, actions: &impl ConstraintFrameActions) -> Option<Self> {
        self.into_iter()
            .map(|constraint| constraint.reframe_by_actions(actions))
            .collect::<Option<Vec<_>>>()
            .map(Self)
    }
}

/// Normalize each child, splice same-combinator children (flatten), drop
/// empty `And`/`Or`, then sort + dedup. `is_and` selects which combinator is the
/// parent: `true` flattens nested `And` (and the conjunctive top-level store),
/// `false` flattens nested `Or`.
fn normalize_logical_constraints(
    constraints: Vec<Constraint>,
    is_and: bool,
) -> Result<Vec<Constraint>, Contradiction> {
    let mut out = Vec::new();
    for child in constraints {
        match child.normalize()? {
            Constraint::And(inner) if is_and => out.extend(inner),
            Constraint::Or(inner) if !is_and => out.extend(inner),
            Constraint::And(inner) | Constraint::Or(inner) if inner.is_empty() => {}
            other => out.push(other),
        }
    }
    out.sort();
    out.dedup();
    Ok(out)
}

/// Molecule-level constraint store: a flat list of `Constraint` tree nodes
/// (molecule-scope predicates, combinators, and entity-leaves that appear
/// inside combinators). Unconditional per-entity constraints live on the
/// entity form's own `constraints` field; the DSL parser lifts them there.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Constraints(Vec<Constraint>);

impl Constraints {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn as_slice(&self) -> &[Constraint] {
        &self.0
    }

    pub fn iter(&self) -> Iter<'_, Constraint> {
        self.0.iter()
    }

    pub fn push(&mut self, c: Constraint) {
        self.0.push(c);
    }

    pub fn remove_at(&mut self, position: usize) -> Constraint {
        self.0.remove(position)
    }

    pub fn retain(&mut self, mut f: impl FnMut(&Constraint) -> bool) {
        self.0.retain(|c| f(c));
    }

    pub fn take(&mut self) -> Vec<Constraint> {
        mem::take(&mut self.0)
    }

    pub fn clear(&mut self) {
        self.0.clear();
    }

    /// Remap entity indices. Entries that reference a removed entity (directly
    /// or via a combinator subtree) are dropped.
    pub fn compact(&mut self, compaction: &MoleculeCompaction) {
        self.0 = mem::take(&mut self.0)
            .into_iter()
            .filter_map(|c| c.compact(compaction))
            .collect();
    }

    /// Remap entity indices and return the patch needed to restore or inspect
    /// constraints that were dropped or rewritten by the compaction.
    pub fn compact_with_update(&mut self, compaction: &MoleculeCompaction) -> CascadedConstraints {
        let mut update = CascadedConstraints::default();
        let mut next = Vec::new();
        for (position, constraint) in mem::take(&mut self.0).into_iter().enumerate() {
            match constraint.clone().compact(compaction) {
                Some(mapped) => {
                    if mapped != constraint {
                        update.modified.push(ModifiedConstraint {
                            position,
                            old: constraint,
                            new: mapped.clone(),
                        });
                    }
                    next.push(mapped);
                }
                None => update.removed.push(RemovedConstraint {
                    position,
                    constraint,
                }),
            }
        }
        self.0 = next;
        update
    }
}

impl FromIterator<Constraint> for Constraints {
    fn from_iter<I: IntoIterator<Item = Constraint>>(iter: I) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl IntoIterator for Constraints {
    type Item = Constraint;
    type IntoIter = IntoIter<Constraint>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl From<Constraint> for Constraints {
    fn from(c: Constraint) -> Self {
        Self(vec![c])
    }
}

impl From<Vec<Constraint>> for Constraints {
    fn from(cs: Vec<Constraint>) -> Self {
        Self(cs)
    }
}

/// Molecule-scope predicates: non-logical, unanchored assertions whose scope
/// is the molecule as a whole or a declared subset of entities.
///
/// For `ChargeSum` / `UnpairedElectronCoupling` / `BondOrderSum` / `Connected`,
/// an `atoms` (or `bonds`) value of `None` denotes the entire molecule's atoms
/// (or bonds), making the predicate stable across structural growth.
/// `Some(vec)` denotes a fixed subset.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MoleculeConstraint {
    ChargeSum {
        atoms: Option<Vec<AtomId>>,
        sum: NumForm,
    },
    UnpairedElectronCoupling {
        atoms: Option<Vec<AtomId>>,
        unpaired_electrons: UnpairedElectronsForm,
    },
    BondOrderSum {
        bonds: Option<Vec<BondId>>,
        sum: NumForm,
    },
    Connected {
        atoms: Option<Vec<AtomId>>,
    },
}

impl MoleculeConstraint {
    /// A constraint is vacuous when its value-bearing payload is
    /// `Undetermined`: `ChargeSum`/`BondOrderSum` with `Undetermined` sum,
    /// `UnpairedElectronCoupling` with both unpaired-electron fields
    /// `Undetermined`. `Connected` is structural and never vacuous in this
    /// sense.
    pub fn is_vacuous(&self) -> bool {
        match self {
            Self::ChargeSum { sum, .. } => sum.is_undetermined(),
            Self::BondOrderSum { sum, .. } => sum.is_undetermined(),
            Self::UnpairedElectronCoupling {
                unpaired_electrons, ..
            } => unpaired_electrons.is_undetermined(),
            Self::Connected { .. } => false,
        }
    }

    pub fn compact(self, compaction: &MoleculeCompaction) -> Option<Self> {
        match self {
            MoleculeConstraint::ChargeSum { atoms, sum } => {
                let atoms = compact_atom_subset(atoms, compaction)?;
                Some(MoleculeConstraint::ChargeSum { atoms, sum })
            }
            MoleculeConstraint::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => {
                let atoms = compact_atom_subset(atoms, compaction)?;
                Some(MoleculeConstraint::UnpairedElectronCoupling {
                    atoms,
                    unpaired_electrons,
                })
            }
            MoleculeConstraint::BondOrderSum { bonds, sum } => {
                let bonds = compact_bond_subset(bonds, compaction)?;
                Some(MoleculeConstraint::BondOrderSum { bonds, sum })
            }
            MoleculeConstraint::Connected { atoms } => {
                let atoms = compact_atom_subset(atoms, compaction)?;
                Some(MoleculeConstraint::Connected { atoms })
            }
        }
    }

    pub fn remap(self, map: &IdRemapping) -> Self {
        match self {
            MoleculeConstraint::ChargeSum { atoms, sum } => MoleculeConstraint::ChargeSum {
                atoms: remap_atom_subset(atoms, map),
                sum,
            },
            MoleculeConstraint::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => MoleculeConstraint::UnpairedElectronCoupling {
                atoms: remap_atom_subset(atoms, map),
                unpaired_electrons,
            },
            MoleculeConstraint::BondOrderSum { bonds, sum } => MoleculeConstraint::BondOrderSum {
                bonds: remap_bond_subset(bonds, map),
                sum,
            },
            MoleculeConstraint::Connected { atoms } => MoleculeConstraint::Connected {
                atoms: remap_atom_subset(atoms, map),
            },
        }
    }
}

impl Normalize for MoleculeConstraint {
    /// Normalize the value payload (`sum` / `unpaired_electrons`) and sort and deduplicate
    /// each atom/bond subset; refs are otherwise unchanged.
    fn normalize(self) -> Result<Self, Contradiction> {
        fn normalize_subset<T: Ord>(mut values: Vec<T>) -> Vec<T> {
            values.sort_unstable();
            values.dedup();
            values
        }
        Ok(match self {
            Self::ChargeSum { atoms, sum } => Self::ChargeSum {
                atoms: atoms.map(normalize_subset),
                sum: sum.normalize()?,
            },
            Self::UnpairedElectronCoupling {
                atoms,
                unpaired_electrons,
            } => Self::UnpairedElectronCoupling {
                atoms: atoms.map(normalize_subset),
                unpaired_electrons: unpaired_electrons.normalize()?,
            },
            Self::BondOrderSum { bonds, sum } => Self::BondOrderSum {
                bonds: bonds.map(normalize_subset),
                sum: sum.normalize()?,
            },
            Self::Connected { atoms } => Self::Connected {
                atoms: atoms.map(normalize_subset),
            },
        })
    }
}

/// Remap an `Option<Vec<AtomId>>`. `None` (all atoms) passes through.
/// `Some(vec)` compacts each element; if any atom was removed the whole
/// constraint is dropped (returns outer `None`).
fn compact_atom_subset(
    atoms: Option<Vec<AtomId>>,
    compaction: &MoleculeCompaction,
) -> Option<Option<Vec<AtomId>>> {
    match atoms {
        None => Some(None),
        Some(vec) => vec
            .into_iter()
            .map(|a| compaction.compact_atom(a))
            .collect::<Option<Vec<_>>>()
            .map(Some),
    }
}

/// Remap an `Option<Vec<BondId>>`. Same semantics as `compact_atom_subset`.
fn compact_bond_subset(
    bonds: Option<Vec<BondId>>,
    compaction: &MoleculeCompaction,
) -> Option<Option<Vec<BondId>>> {
    match bonds {
        None => Some(None),
        Some(vec) => vec
            .into_iter()
            .map(|b| compaction.compact_bond(b))
            .collect::<Option<Vec<_>>>()
            .map(Some),
    }
}

/// Re-anchor an `Option<Vec<AtomId>>` through a total id remapping (the parallel of
/// `compact_atom_subset`). `None` (all atoms) passes through; total — never drops.
fn remap_atom_subset(atoms: Option<Vec<AtomId>>, map: &IdRemapping) -> Option<Vec<AtomId>> {
    atoms.map(|vec| vec.into_iter().map(|a| map.map_atom(a)).collect())
}

/// Re-anchor an `Option<Vec<BondId>>`. Same semantics as `remap_atom_subset`.
fn remap_bond_subset(bonds: Option<Vec<BondId>>, map: &IdRemapping) -> Option<Vec<BondId>> {
    bonds.map(|vec| vec.into_iter().map(|b| map.map_bond(b)).collect())
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::{EdgeId, GraphCompaction, NodeId};
    use umol_perm::{DynPermutation, Permutation};

    use super::*;
    use crate::ir::constraint::RingScope;
    use crate::ir::frame::{
        AromaticSystemsFrameAction, DativeBondsFrameAction, MulticenterBondsFrameAction,
        NoncovalentBondsFrameAction, StereoAtomsFrameAction, StereoBondsFrameAction,
    };
    use crate::ir::id::{
        AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
        StereoLigandPosition,
    };
    use crate::ir::num::{ArithExpr, NumForm};
    use crate::ir::spin::UnpairedElectronsForm;
    use crate::ir::{
        BooleanForm, StereoLigandPair, Stereogenicity, StereogenicityForm, Topicity, TopicityForm,
        TopicityRelationForm,
    };

    #[fixture]
    fn overlays_frame_action() -> OverlaysFrameAction {
        OverlaysFrameAction::new(
            DativeBondsFrameAction::from_vec(vec![]).expect("actions are admissible"),
            AromaticSystemsFrameAction::from_vec(vec![]).expect("actions are admissible"),
            MulticenterBondsFrameAction::from_vec(vec![]).expect("actions are admissible"),
            NoncovalentBondsFrameAction::from_vec(vec![
                DynPermutation::try_from(vec![1, 0]).expect("action is a permutation")
            ])
            .expect("action is admissible"),
            StereoAtomsFrameAction::from_vec(vec![Permutation::from_image(&[1, 0, 2, 3])])
                .expect("action is admissible"),
            StereoBondsFrameAction::from_vec(vec![]).expect("actions are admissible"),
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::frame_invariant(
        Constraints::from(vec![Constraint::And(vec![
            Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(true)),
            Constraint::AromaticSystem(AromaticSystemId(0), AromaticSystemConstraintForm::electron_count(6)),
            Constraint::MulticenterBond(MulticenterBondId(0), MulticenterBondConstraintForm::electron_count(2)),
            Constraint::NoncovalentBond(NoncovalentBondId(0), NoncovalentBondConstraintForm::intramolecular(true)),
            Constraint::StereoAtom(
                StereoAtomId(0),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
            ),
            Constraint::StereoBond(
                StereoBondId(0),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
            ),
        ])]),
        [false, false, false, false, false, false],
    )]
    #[case::frame_relative(
        Constraints::from(vec![Constraint::Not(Box::new(Constraint::And(vec![
            Constraint::Relational(RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::valence(4)),
                    Box::new(AtomConstraintForm::degree(2)),
                ],
            }),
            Constraint::StereoAtom(
                StereoAtomId(0),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
                    relation: TopicityRelationForm::Lit(Topicity::Homotopic),
                }),
            ),
            Constraint::StereoBond(
                StereoBondId(0),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
                    relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                }),
            ),
        ])))]),
        [false, false, false, true, true, true],
    )]
    fn test_constraints_frame_action_domain(
        #[case] constraints: Constraints,
        #[case] expected: [bool; 6],
    ) {
        let domain = constraints.frame_action_domain();

        assert_eq!(domain.contains_dative_bond(DativeBondId(0)), expected[0]);
        assert_eq!(domain.contains_aromatic_system(AromaticSystemId(0)), expected[1]);
        assert_eq!(domain.contains_multicenter_bond(MulticenterBondId(0)), expected[2]);
        assert_eq!(domain.contains_noncovalent_bond(NoncovalentBondId(0)), expected[3]);
        assert_eq!(domain.contains_stereo_atom(StereoAtomId(0)), expected[4]);
        assert_eq!(domain.contains_stereo_bond(StereoBondId(0)), expected[5]);
        assert_eq!(domain.count(EntityKind::DativeBond), usize::from(expected[0]));
        assert_eq!(domain.count(EntityKind::AromaticSystem), usize::from(expected[1]));
        assert_eq!(domain.count(EntityKind::MulticenterBond), usize::from(expected[2]));
        assert_eq!(domain.count(EntityKind::NoncovalentBond), usize::from(expected[3]));
        assert_eq!(domain.count(EntityKind::StereoAtom), usize::from(expected[4]));
        assert_eq!(domain.count(EntityKind::StereoBond), usize::from(expected[5]));
    }

    fn id_compaction(removed_nodes: Vec<u32>, removed_edges: Vec<u32>) -> MoleculeCompaction {
        MoleculeCompaction::new(
            GraphCompaction::new(
                removed_nodes.into_iter().map(NodeId).collect(),
                removed_edges.into_iter().map(EdgeId).collect(),
            ),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    }

    fn relation_compaction(
        removed_dative: Vec<u32>,
        removed_aromatic: Vec<u32>,
        removed_multicenter: Vec<u32>,
        removed_noncovalent: Vec<u32>,
        removed_stereo_atoms: Vec<u32>,
        removed_stereo_bonds: Vec<u32>,
    ) -> MoleculeCompaction {
        MoleculeCompaction::new(
            GraphCompaction::new(Vec::new(), Vec::new()),
            removed_dative.into_iter().map(DativeBondId).collect(),
            removed_aromatic.into_iter().map(AromaticSystemId).collect(),
            removed_multicenter
                .into_iter()
                .map(MulticenterBondId)
                .collect(),
            removed_noncovalent
                .into_iter()
                .map(NoncovalentBondId)
                .collect(),
            removed_stereo_atoms.into_iter().map(StereoAtomId).collect(),
            removed_stereo_bonds.into_iter().map(StereoBondId).collect(),
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_lit(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)), false)]
    #[case::atom_undetermined(Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::Undetermined)), true)]
    #[case::bond_lit(Constraint::Bond(BondId(0), BondConstraintForm::ring_membership(RingScope::Size(6), 1)), false)]
    #[case::bond_undetermined(Constraint::Bond(BondId(0), BondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined)), true)]
    #[case::bond_aromatic_flag(Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))), false)]
    #[case::dative_undetermined(Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined)), true)]
    #[case::aromatic_system_undetermined(Constraint::AromaticSystem(AromaticSystemId(0),
        AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)), true)]
    #[case::multicenter_undetermined(Constraint::MulticenterBond(MulticenterBondId(0),
        MulticenterBondConstraintForm::ElectronCount(NumForm::Undetermined)), true)]
    #[case::relational(Constraint::Relational(RelationalConstraint::DativeBondDonor {
        bond: DativeBondId(0), atom: AtomId(0) }), false)]
    #[case::molecule_undetermined(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: None, sum: NumForm::Undetermined }), true)]
    #[case::molecule_lit(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: None, sum: NumForm::Lit(0) }), false)]
    #[case::and_empty(Constraint::And(vec![]), true)]
    #[case::and_nonempty(Constraint::And(vec![Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))]), false)]
    #[case::or_empty(Constraint::Or(vec![]), true)]
    #[case::or_nonempty(Constraint::Or(vec![Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true)))]), false)]
    #[case::not(Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)))), false)]
    fn test_constraint_is_vacuous(#[case] c: Constraint, #[case] expected: bool) {
        assert_eq!(c.is_vacuous(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::leaf_folds(
        Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::arith_expr(ArithExpr::Lit(4)))),
        Ok(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))),
    )]
    #[case::and_flattens_nested(
        Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::And(vec![Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true)))]),
        ]),
        Ok(Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        ])),
    )]
    #[case::and_drops_empty_or_child_and_reduces(
        Constraint::And(vec![Constraint::Or(vec![]), Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))]),
        Ok(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))),
    )]
    #[case::and_sorts_and_dedups(
        Constraint::And(vec![
            Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
        ]),
        Ok(Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        ])),
    )]
    #[case::or_flattens_nested(
        Constraint::Or(vec![
            Constraint::Or(vec![Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))]),
            Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        ]),
        Ok(Constraint::Or(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        ])),
    )]
    #[case::or_drops_empty_and_child_and_reduces(
        Constraint::Or(vec![Constraint::And(vec![]), Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true)))]),
        Ok(Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true)))),
    )]
    #[case::and_singleton_reduces(
        Constraint::And(vec![Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))]),
        Ok(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))),
    )]
    #[case::or_singleton_reduces(
        Constraint::Or(vec![Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))]),
        Ok(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))),
    )]
    #[case::and_dedup_reduces_to_element(
        Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
        ]),
        Ok(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))),
    )]
    #[case::cross_combinator_singleton_reduces(
        Constraint::And(vec![Constraint::Or(vec![Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))])]),
        Ok(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))),
    )]
    #[case::not_singleton_child_reduces(
        Constraint::Not(Box::new(Constraint::Or(vec![Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))]))),
        Ok(Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))))),
    )]
    #[case::not_folds_child(
        Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::arith_expr(ArithExpr::Lit(4)))))),
        Ok(Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4))))),
    )]
    #[case::inner_contradiction_propagates(
        Constraint::And(vec![Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::lit_set(Vec::<i64>::new())))]),
        Err(Contradiction),
    )]
    fn test_constraint_normalize(
        #[case] input: Constraint,
        #[case] expected: Result<Constraint, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_bond(
        Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(true)),
        Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(true)),
    )]
    #[case::aromatic_system(
        Constraint::AromaticSystem(AromaticSystemId(0), AromaticSystemConstraintForm::electron_count(6)),
        Constraint::AromaticSystem(AromaticSystemId(0), AromaticSystemConstraintForm::electron_count(6)),
    )]
    #[case::multicenter_bond(
        Constraint::MulticenterBond(MulticenterBondId(0), MulticenterBondConstraintForm::electron_count(2)),
        Constraint::MulticenterBond(MulticenterBondId(0), MulticenterBondConstraintForm::electron_count(2)),
    )]
    #[case::noncovalent_bond(
        Constraint::NoncovalentBond(NoncovalentBondId(1), NoncovalentBondConstraintForm::intramolecular(true)),
        Constraint::NoncovalentBond(NoncovalentBondId(1), NoncovalentBondConstraintForm::intramolecular(true)),
    )]
    #[case::stereo_atom_stereogenicity(
        Constraint::StereoAtom(StereoAtomId(1), StereoKind::Tetrahedral, StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        Constraint::StereoAtom(StereoAtomId(1), StereoKind::Tetrahedral, StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
    )]
    #[case::stereo_bond_stereogenicity(
        Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
        Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
    )]
    #[case::recursive_relational(
        Constraint::Not(Box::new(Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                ],
            },
        ))),
        Constraint::Not(Box::new(Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                ],
            },
        ))),
    )]
    #[case::stereo_atom_topicity(
        Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(2)),
                relation: TopicityRelationForm::Undetermined,
            }),
        ),
        Constraint::StereoAtom(
            StereoAtomId(0),
            StereoKind::Tetrahedral,
            StereoAtomConstraintForm::Topicity(TopicityForm {
                pair: StereoLigandPair::new(StereoLigandPosition(1), StereoLigandPosition(2)),
                relation: TopicityRelationForm::Undetermined,
            }),
        ),
    )]
    fn test_constraint_reframe_by(
        #[case] input: Constraint,
        #[case] expected: Constraint,
        overlays_frame_action: OverlaysFrameAction,
    ) {
        assert_eq!(input.reframe_by(&overlays_frame_action), Some(expected));
    }

    #[rstest]
    fn test_constraints_reframe_by(overlays_frame_action: OverlaysFrameAction) {
        let input = Constraints::from(vec![Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                ],
            },
        )]);
        let expected = Constraints::from(vec![Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                ],
            },
        )]);

        assert_eq!(input.reframe_by(&overlays_frame_action), Some(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_shifts(
        Constraint::Atom(AtomId(2), AtomConstraintForm::valence(4)),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Atom(AtomId(1), AtomConstraintForm::valence(4))),
    )]
    #[case::atom_dropped(
        Constraint::Atom(AtomId(1), AtomConstraintForm::valence(4)),
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::bond_shifts(
        Constraint::Bond(BondId(3), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        id_compaction(vec![], vec![1]),
        Some(Constraint::Bond(BondId(2), BondConstraintForm::Aromatic(BooleanForm::Lit(true)))),
    )]
    #[case::bond_dropped(
        Constraint::Bond(BondId(1), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        id_compaction(vec![], vec![1]),
        None,
    )]
    #[case::dative_shifts(
        Constraint::DativeBond(DativeBondId(2), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        relation_compaction(vec![0], vec![], vec![], vec![], vec![], vec![]),
        Some(Constraint::DativeBond(DativeBondId(1), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)))),
    )]
    #[case::dative_dropped(
        Constraint::DativeBond(DativeBondId(1), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        relation_compaction(vec![1], vec![], vec![], vec![], vec![], vec![]),
        None,
    )]
    #[case::aromatic_system_shifts(
        Constraint::AromaticSystem(AromaticSystemId(2), AromaticSystemConstraintForm::electron_count(6)),
        relation_compaction(vec![], vec![0], vec![], vec![], vec![], vec![]),
        Some(Constraint::AromaticSystem(AromaticSystemId(1), AromaticSystemConstraintForm::electron_count(6))),
    )]
    #[case::aromatic_system_dropped(
        Constraint::AromaticSystem(AromaticSystemId(1), AromaticSystemConstraintForm::electron_count(6)),
        relation_compaction(vec![], vec![1], vec![], vec![], vec![], vec![]),
        None,
    )]
    #[case::multicenter_shifts(
        Constraint::MulticenterBond(MulticenterBondId(2), MulticenterBondConstraintForm::electron_count(2)),
        relation_compaction(vec![], vec![], vec![0], vec![], vec![], vec![]),
        Some(Constraint::MulticenterBond(MulticenterBondId(1), MulticenterBondConstraintForm::electron_count(2))),
    )]
    #[case::multicenter_dropped(
        Constraint::MulticenterBond(MulticenterBondId(1), MulticenterBondConstraintForm::electron_count(2)),
        relation_compaction(vec![], vec![], vec![1], vec![], vec![], vec![]),
        None,
    )]
    #[case::relational_dative_donor_shifts_atom(
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(2), atom: AtomId(3) }),
        id_compaction(vec![0], vec![]),
        Some(Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(2), atom: AtomId(2) })),
    )]
    #[case::relational_dative_donor_dropped_when_bond_removed(
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(1), atom: AtomId(0) }),
        relation_compaction(vec![1], vec![], vec![], vec![], vec![], vec![]),
        None,
    )]
    #[case::relational_aromatic_system_contains_shifts_atom(
        Constraint::Relational(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(1), atom: AtomId(2) }),
        id_compaction(vec![0], vec![]),
        Some(Constraint::Relational(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(1), atom: AtomId(1) })),
    )]
    #[case::relational_multicenter_contains_shifts_atom(
        Constraint::Relational(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondId(0), atom: AtomId(2) }),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Relational(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondId(0), atom: AtomId(1) })),
    )]
    #[case::relational_noncovalent_contains_shifts_atom(
        Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(0), atom: AtomId(3) }),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Relational(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(0), atom: AtomId(2) })),
    )]
    #[case::molecule_charge_sum_shifts(
        Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(2)]), sum: NumForm::Lit(1) }),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(1)]), sum: NumForm::Lit(1) })),
    )]
    #[case::and_all_survive(
        Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintForm::valence(2)),
        ]),
        id_compaction(vec![1], vec![]),
        Some(Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(1), AtomConstraintForm::valence(2)),
        ])),
    )]
    #[case::and_drops_if_any_leaf_drops(
        Constraint::And(vec![
            Constraint::Atom(AtomId(1), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintForm::valence(2)),
        ]),
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::or_all_survive(
        Constraint::Or(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintForm::valence(2)),
        ]),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Or(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(1), AtomConstraintForm::valence(2)),
        ])),
    )]
    #[case::or_drops_if_any_leaf_drops(
        Constraint::Or(vec![
            Constraint::Atom(AtomId(1), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintForm::valence(2)),
        ]),
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::not_wraps_child(
        Constraint::Not(Box::new(Constraint::Atom(AtomId(2), AtomConstraintForm::valence(4)))),
        id_compaction(vec![1], vec![]),
        Some(Constraint::Not(Box::new(Constraint::Atom(AtomId(1), AtomConstraintForm::valence(4))))),
    )]
    #[case::not_drops_child(
        Constraint::Not(Box::new(Constraint::Atom(AtomId(1), AtomConstraintForm::valence(4)))),
        id_compaction(vec![1], vec![]),
        None,
    )]
    fn test_constraint_compact(
        #[case] c: Constraint,
        #[case] compaction: MoleculeCompaction,
        #[case] expected: Option<Constraint>,
    ) {
        assert_eq!(c.compact(&compaction), expected);
    }

    fn id_remapping(
        atom: &[(u32, u32)],
        bond: &[(u32, u32)],
        dative: &[(u32, u32)],
    ) -> IdRemapping {
        IdRemapping::new(
            atom.iter().map(|&(a, b)| (AtomId(a), AtomId(b))).collect(),
            bond.iter().map(|&(a, b)| (BondId(a), BondId(b))).collect(),
            dative
                .iter()
                .map(|&(a, b)| (DativeBondId(a), DativeBondId(b)))
                .collect(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
            Default::default(),
        )
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom(
        Constraint::Atom(AtomId(2), AtomConstraintForm::valence(4)),
        id_remapping(&[(2, 5)], &[], &[]),
        Constraint::Atom(AtomId(5), AtomConstraintForm::valence(4)),
    )]
    #[case::bond(
        Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        id_remapping(&[], &[(0, 3)], &[]),
        Constraint::Bond(BondId(3), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
    )]
    #[case::dative_leaf(
        Constraint::DativeBond(DativeBondId(1), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        id_remapping(&[], &[], &[(1, 0)]),
        Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
    )]
    #[case::molecule_charge_sum(
        Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(2)]), sum: NumForm::Lit(1) }),
        id_remapping(&[(0, 3), (2, 4)], &[], &[]),
        Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(3), AtomId(4)]), sum: NumForm::Lit(1) }),
    )]
    #[case::molecule_unpaired_electron_coupling_subset(
        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0), AtomId(2)]), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) }),
        id_remapping(&[(0, 3), (2, 4)], &[], &[]),
        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(3), AtomId(4)]), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) }),
    )]
    #[case::molecule_unpaired_electron_coupling_all_atoms(
        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) }),
        id_remapping(&[(0, 3), (2, 4)], &[], &[]),
        Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) }),
    )]
    #[case::relational_dative_donor(
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(1), atom: AtomId(2) }),
        id_remapping(&[(2, 5)], &[], &[(1, 0)]),
        Constraint::Relational(RelationalConstraint::DativeBondDonor { bond: DativeBondId(0), atom: AtomId(5) }),
    )]
    #[case::and(
        Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintForm::valence(2)),
        ]),
        id_remapping(&[(0, 1), (2, 3)], &[], &[]),
        Constraint::And(vec![
            Constraint::Atom(AtomId(1), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(3), AtomConstraintForm::valence(2)),
        ]),
    )]
    #[case::not(
        Constraint::Not(Box::new(Constraint::Atom(AtomId(2), AtomConstraintForm::valence(4)))),
        id_remapping(&[(2, 0)], &[], &[]),
        Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)))),
    )]
    fn test_constraint_remap(
        #[case] c: Constraint,
        #[case] map: IdRemapping,
        #[case] expected: Constraint,
    ) {
        assert_eq!(c.remap(&map), expected);
    }

    #[rstest]
    fn test_constraints_new() {
        let cs = Constraints::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
        assert_eq!(cs.as_slice(), &[] as &[Constraint]);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(vec![], 0)]
    #[case::molecule_leaves(vec![
            Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(1)]), sum: NumForm::Lit(0) }),
            Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0)]), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) }),
        ], 2)]
    #[case::combinator(vec![Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        ])], 1)]
    fn test_constraints_push(
        #[case] items: Vec<Constraint>,
        #[case] expected_len: usize,
    ) {
        let mut cs = Constraints::new();
        for c in items {
            cs.push(c);
        }
        assert_eq!(cs.len(), expected_len);
        assert_eq!(cs.is_empty(), expected_len == 0);
    }

    #[rstest]
    fn test_constraints_retain() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: NumForm::Lit(0),
        }));
        cs.push(Constraint::And(vec![]));

        cs.retain(|c| matches!(c, Constraint::Molecule(_)));
        assert_eq!(cs.len(), 1);
    }

    #[rstest]
    fn test_constraints_take() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
        }));
        cs.push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: NumForm::Lit(0),
        }));

        let taken = cs.take();
        assert_eq!(taken.len(), 2);
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_constraints_clear() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0)]),
        }));
        cs.clear();
        assert!(cs.is_empty());
    }

    #[rstest]
    fn test_constraints_iter() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)));
        cs.push(Constraint::Bond(
            BondId(0),
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
        ));
        let collected: Vec<_> = cs.iter().cloned().collect();
        assert_eq!(
            collected,
            vec![
                Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
                Constraint::Bond(
                    BondId(0),
                    BondConstraintForm::Aromatic(BooleanForm::Lit(true))
                ),
            ],
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::drops_entity_leaf_on_removed_atom(
        vec![
            Constraint::Atom(AtomId(1), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintForm::valence(3)),
        ],
        id_compaction(vec![1], vec![]),
        vec![Constraint::Atom(AtomId(1), AtomConstraintForm::valence(3))],
    )]
    #[case::shifts_remaining_leaves(
        vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(2), AtomConstraintForm::degree(3)),
        ],
        id_compaction(vec![1], vec![]),
        vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(1), AtomConstraintForm::degree(3)),
        ],
    )]
    #[case::drops_combinator_if_any_leaf_dropped(
        vec![Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(1), AtomConstraintForm::degree(3)),
        ])],
        id_compaction(vec![1], vec![]),
        vec![],
    )]
    fn test_constraints_compact(
        #[case] items: Vec<Constraint>,
        #[case] compaction: MoleculeCompaction,
        #[case] expected: Vec<Constraint>,
    ) {
        let mut cs = Constraints::new();
        for c in items {
            cs.push(c);
        }
        cs.compact(&compaction);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_constraints_compact_with_update() {
        let mut cs = Constraints::new();
        cs.push(Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)));
        cs.push(Constraint::Atom(AtomId(1), AtomConstraintForm::degree(3)));
        cs.push(Constraint::Bond(
            BondId(2),
            BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        ));
        cs.push(Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(2)]),
        }));

        let update = cs.compact_with_update(&id_compaction(vec![1], vec![1]));

        assert_eq!(
            cs.as_slice(),
            &[
                Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
                Constraint::Bond(
                    BondId(1),
                    BondConstraintForm::ring_membership(RingScope::Size(6), 1)
                ),
                Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(0), AtomId(1)]),
                }),
            ],
        );
        assert_eq!(
            update,
            CascadedConstraints {
                removed: vec![RemovedConstraint {
                    position: 1,
                    constraint: Constraint::Atom(AtomId(1), AtomConstraintForm::degree(3)),
                }],
                modified: vec![
                    ModifiedConstraint {
                        position: 2,
                        old: Constraint::Bond(
                            BondId(2),
                            BondConstraintForm::ring_membership(RingScope::Size(6), 1)
                        ),
                        new: Constraint::Bond(
                            BondId(1),
                            BondConstraintForm::ring_membership(RingScope::Size(6), 1)
                        ),
                    },
                    ModifiedConstraint {
                        position: 3,
                        old: Constraint::Molecule(MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(0), AtomId(2)]),
                        }),
                        new: Constraint::Molecule(MoleculeConstraint::Connected {
                            atoms: Some(vec![AtomId(0), AtomId(1)]),
                        }),
                    },
                ],
            },
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::flattens_top_level_and_then_sorts(
        Constraints::from(vec![
            Constraint::And(vec![Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true)))]),
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
        ]),
        Ok(Constraints::from(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        ])),
    )]
    #[case::drops_empty_or_and_dedups(
        Constraints::from(vec![
            Constraint::Or(vec![]),
            Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        ]),
        Ok(Constraints::from(vec![Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true)))])),
    )]
    #[case::inner_contradiction_propagates(
        Constraints::from(vec![Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::lit_set(Vec::<i64>::new())))]),
        Err(Contradiction),
    )]
    fn test_constraints_normalize(
        #[case] input: Constraints,
        #[case] expected: Result<Constraints, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_sum_lit(MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Lit(0) }, false)]
    #[case::charge_sum_undetermined(MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Undetermined }, true)]
    #[case::unpaired_electron_coupling_ground(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) }, false)]
    #[case::unpaired_electron_coupling_undetermined(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsForm::default() }, true)]
    #[case::bond_order_sum_lit(MoleculeConstraint::BondOrderSum { bonds: None, sum: NumForm::Lit(4) }, false)]
    #[case::bond_order_sum_undetermined(MoleculeConstraint::BondOrderSum { bonds: None, sum: NumForm::Undetermined }, true)]
    #[case::connected(MoleculeConstraint::Connected { atoms: None }, false)]
    fn test_molecule_constraint_is_vacuous(
        #[case] c: MoleculeConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.is_vacuous(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_sum_sorts_and_folds(
        MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(2), AtomId(0), AtomId(2)]), sum: NumForm::arith_expr(ArithExpr::Lit(1)) },
        Ok(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(2)]), sum: NumForm::Lit(1) }),
    )]
    #[case::unpaired_electron_coupling_sorts_and_folds(
        MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(2), AtomId(0)]),
            unpaired_electrons: UnpairedElectronsForm { count: NumForm::arith_expr(ArithExpr::Lit(0)), multiplicity: NumForm::arith_expr(ArithExpr::Lit(1)) } },
        Ok(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0), AtomId(2)]), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) }),
    )]
    #[case::bond_order_sum_sorts_and_folds(
        MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(2), BondId(0), BondId(2)]), sum: NumForm::arith_expr(ArithExpr::Lit(4)) },
        Ok(MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(0), BondId(2)]), sum: NumForm::Lit(4) }),
    )]
    #[case::connected_sorts(
        MoleculeConstraint::Connected { atoms: Some(vec![AtomId(3), AtomId(1), AtomId(2)]) },
        Ok(MoleculeConstraint::Connected { atoms: Some(vec![AtomId(1), AtomId(2), AtomId(3)]) }),
    )]
    #[case::charge_sum_empty_litset_contradiction(
        MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::lit_set(Vec::<i64>::new()) },
        Err(Contradiction),
    )]
    fn test_molecule_constraint_normalize(
        #[case] input: MoleculeConstraint,
        #[case] expected: Result<MoleculeConstraint, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rstest]
    #[case::connected_none(MoleculeConstraint::Connected { atoms: None })]
    fn test_molecule_constraint_normalize_identity(#[case] input: MoleculeConstraint) {
        assert_eq!(input.clone().normalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_before_connected(
        MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Lit(0) },
        MoleculeConstraint::Connected { atoms: None },
        Ordering::Less,
    )]
    fn test_molecule_constraint_cmp(
        #[case] a: MoleculeConstraint,
        #[case] b: MoleculeConstraint,
        #[case] expected: Ordering,
    ) {
        assert_eq!(a.cmp(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_sum_shifts(
        MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(2)]), sum: NumForm::Lit(1) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(1)]), sum: NumForm::Lit(1) }),
    )]
    #[case::charge_sum_drops_when_atom_removed(
        MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(1), AtomId(2)]), sum: NumForm::Lit(0) },
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::charge_sum_all_atoms_passes_through(
        MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Lit(0) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Lit(0) }),
    )]
    #[case::unpaired_electron_coupling_shifts(
        MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0), AtomId(2)]), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(0), AtomId(1)]), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) }),
    )]
    #[case::unpaired_electron_coupling_drops_when_atom_removed(
        MoleculeConstraint::UnpairedElectronCoupling { atoms: Some(vec![AtomId(1), AtomId(2)]), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) },
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::unpaired_electron_coupling_all_atoms_passes_through(
        MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::UnpairedElectronCoupling { atoms: None, unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)) }),
    )]
    #[case::bond_order_sum_shifts(
        MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(0), BondId(2)]), sum: NumForm::Lit(4) },
        id_compaction(vec![], vec![1]),
        Some(MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(0), BondId(1)]), sum: NumForm::Lit(4) }),
    )]
    #[case::bond_order_sum_drops_when_bond_removed(
        MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(1)]), sum: NumForm::Lit(2) },
        id_compaction(vec![], vec![1]),
        None,
    )]
    #[case::bond_order_sum_all_bonds_passes_through(
        MoleculeConstraint::BondOrderSum { bonds: None, sum: NumForm::Lit(0) },
        id_compaction(vec![], vec![1]),
        Some(MoleculeConstraint::BondOrderSum { bonds: None, sum: NumForm::Lit(0) }),
    )]
    #[case::connected_shifts(
        MoleculeConstraint::Connected { atoms: Some(vec![AtomId(0), AtomId(2), AtomId(3)]) },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::Connected { atoms: Some(vec![AtomId(0), AtomId(1), AtomId(2)]) }),
    )]
    #[case::connected_drops_when_atom_removed(
        MoleculeConstraint::Connected { atoms: Some(vec![AtomId(1)]) },
        id_compaction(vec![1], vec![]),
        None,
    )]
    #[case::connected_all_atoms_passes_through(
        MoleculeConstraint::Connected { atoms: None },
        id_compaction(vec![1], vec![]),
        Some(MoleculeConstraint::Connected { atoms: None }),
    )]
    fn test_molecule_constraint_compact(
        #[case] c: MoleculeConstraint,
        #[case] compaction: MoleculeCompaction,
        #[case] expected: Option<MoleculeConstraint>,
    ) {
        assert_eq!(c.compact(&compaction), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(
        vec![Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true)))],
        vec![Constraint::Bond(BondId(0), BondConstraintForm::Aromatic(BooleanForm::Lit(true)))],
    )]
    #[case::preserves_order_and_duplicates(
        vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3)),
        ],
        vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3)),
        ],
    )]
    #[case::empty(vec![], vec![])]
    fn test_constraints_from_iter(
        #[case] input: Vec<Constraint>,
        #[case] expected: Vec<Constraint>,
    ) {
        let cs = Constraints::from_iter(input);
        assert_eq!(cs.as_slice(), expected.as_slice());
    }

    #[rstest]
    fn test_constraints_into_iter() {
        let cs = Constraints::from_iter([
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Bond(
                BondId(0),
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            ),
        ]);
        let collected: Vec<_> = cs.into_iter().collect();
        assert_eq!(
            collected,
            vec![
                Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
                Constraint::Bond(
                    BondId(0),
                    BondConstraintForm::Aromatic(BooleanForm::Lit(true))
                ),
            ],
        );
    }

    #[rstest]
    fn test_constraints_from_constraint() {
        let cs: Constraints = Constraint::Bond(
            BondId(0),
            BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
        )
        .into();
        assert_eq!(
            cs.as_slice(),
            &[Constraint::Bond(
                BondId(0),
                BondConstraintForm::Aromatic(BooleanForm::Lit(true))
            )],
        );
    }

    #[rstest]
    fn test_constraints_from_vec() {
        let cs: Constraints = vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
            Constraint::Bond(
                BondId(0),
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            ),
        ]
        .into();
        assert_eq!(
            cs.as_slice(),
            &[
                Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4)),
                Constraint::Bond(
                    BondId(0),
                    BondConstraintForm::Aromatic(BooleanForm::Lit(true))
                ),
            ],
        );
    }
}
