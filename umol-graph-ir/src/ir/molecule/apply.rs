//! Single-edit execution and undo replay on molecule storage.

use std::borrow::Cow;
use std::collections::HashSet;
use std::hash::Hash;

use umol_graph_core::{Compaction, GraphCompaction};
use umol_perm::{ClassKey, Permutation, MAX_DEGREE};

use super::{Molecule, TransactionError};
use crate::ir::aromatic::AromaticSystemEntry;
use crate::ir::compact::MoleculeCompaction;
use crate::ir::constraint::Constraint;
use crate::ir::correspondence::MoleculeCorrespondence;
use crate::ir::dative::DativeBondEntry;
use crate::ir::edit::{
    AddAromaticSystem, AddBond, AddMulticenterBond, AddStereoAtom, AddStereoBond,
    AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, AddedStereoAtom, AddedStereoBond, AromaticSystemFieldChange,
    AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange, BondHandle,
    CascadedConstraints, ConstraintEdit, DativeBondFieldChange, DativeBondHandle, Edit,
    MulticenterBondFieldChange, MulticenterBondHandle, NoncovalentBondFieldChange,
    NoncovalentBondHandle, RemovedAromaticSystem, RemovedAromaticSystems, RemovedAtom, RemovedBond,
    RemovedConstraint, RemovedDativeBond, RemovedDativeBonds, RemovedMulticenterBond,
    RemovedMulticenterBonds, RemovedNoncovalentBond, RemovedNoncovalentBonds, RemovedOverlays,
    RemovedStereoAtom, RemovedStereoAtoms, RemovedStereoBond, RemovedStereoBonds, RemovedTopology,
    StereoAtomFieldChange, StereoAtomHandle, StereoBondFieldChange, StereoBondHandle, Undo,
};
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::ligand::{StereoLigand, StereoLigandKind};
use crate::ir::multicenter::MulticenterBondEntry;
use crate::ir::noncovalent::NoncovalentBondEntry;
use crate::ir::stereo::{StereoAtomEntry, StereoBondEntry};
use crate::ir::traits::{Normalize, Reframe};

#[derive(Default)]
struct HandleTable<I> {
    initial_count: usize,
    initial: Option<Vec<Option<I>>>,
    created: Vec<Option<I>>,
}

impl<I: Copy + From<usize>> HandleTable<I> {
    fn new(initial_count: usize) -> Self {
        Self {
            initial_count,
            initial: None,
            created: Vec::new(),
        }
    }

    fn initial(&self, kind: EntityKind, index: usize) -> Result<I, TransactionError> {
        if index >= self.initial_count {
            return Err(TransactionError::HandleOutOfRange {
                kind,
                index,
                count: self.initial_count,
            });
        }
        match &self.initial {
            Some(initial) => Self::resolve(initial, kind, index),
            None => Ok(I::from(index)),
        }
    }

    fn created(&self, kind: EntityKind, index: usize) -> Result<I, TransactionError> {
        Self::resolve(&self.created, kind, index)
    }

    fn removed(&self) -> Vec<I> {
        self.initial
            .iter()
            .flat_map(|initial| initial.iter().enumerate())
            .filter(|(_, current)| current.is_none())
            .map(|(index, _)| I::from(index))
            .collect()
    }

    fn push(&mut self, id: I) {
        self.created.push(Some(id));
    }

    fn compact(&mut self, mut compact: impl FnMut(I) -> Option<I>) {
        match &mut self.initial {
            Some(initial) => {
                for id in initial {
                    if let Some(current) = *id {
                        *id = compact(current);
                    }
                }
            }
            None => {
                self.initial = Some(
                    (0..self.initial_count)
                        .map(|index| compact(I::from(index)))
                        .collect(),
                );
            }
        }
        for id in &mut self.created {
            if let Some(current) = *id {
                *id = compact(current);
            }
        }
    }

    fn resolve(
        entries: &[Option<I>],
        kind: EntityKind,
        index: usize,
    ) -> Result<I, TransactionError> {
        entries
            .get(index)
            .ok_or(TransactionError::HandleOutOfRange {
                kind,
                index,
                count: entries.len(),
            })?
            .ok_or(TransactionError::HandleRemoved { kind, index })
    }
}

/// Resolves initial and newly created handles within one edit batch.
///
/// Initial ids remain implicit until removal compacts them. Removed entries retain
/// their handle positions in both namespaces.
pub(crate) struct ApplicationState {
    atoms: HandleTable<AtomId>,
    bonds: HandleTable<BondId>,
    dative_bonds: HandleTable<DativeBondId>,
    aromatic_systems: HandleTable<AromaticSystemId>,
    multicenter_bonds: HandleTable<MulticenterBondId>,
    noncovalent_bonds: HandleTable<NoncovalentBondId>,
    stereo_atoms: HandleTable<StereoAtomId>,
    stereo_bonds: HandleTable<StereoBondId>,
}

impl ApplicationState {
    pub(crate) fn new(molecule: &Molecule) -> Self {
        Self {
            atoms: HandleTable::new(molecule.atoms().count()),
            bonds: HandleTable::new(molecule.bonds().count()),
            dative_bonds: HandleTable::new(molecule.dative_bonds().count()),
            aromatic_systems: HandleTable::new(molecule.aromatic_systems().count()),
            multicenter_bonds: HandleTable::new(molecule.multicenter_bonds().count()),
            noncovalent_bonds: HandleTable::new(molecule.noncovalent_bonds().count()),
            stereo_atoms: HandleTable::new(molecule.stereo_atoms().count()),
            stereo_bonds: HandleTable::new(molecule.stereo_bonds().count()),
        }
    }

    fn atom(&self, handle: AtomHandle) -> Result<AtomId, TransactionError> {
        match handle {
            AtomHandle::Id(id) => self.atoms.initial(EntityKind::Atom, id.index()),
            AtomHandle::New(index) => self.atoms.created(EntityKind::Atom, index),
        }
    }

    fn bond(&self, handle: BondHandle) -> Result<BondId, TransactionError> {
        match handle {
            BondHandle::Id(id) => self.bonds.initial(EntityKind::Bond, id.index()),
            BondHandle::New(index) => self.bonds.created(EntityKind::Bond, index),
        }
    }

    fn dative_bond(&self, handle: DativeBondHandle) -> Result<DativeBondId, TransactionError> {
        match handle {
            DativeBondHandle::Id(id) => self
                .dative_bonds
                .initial(EntityKind::DativeBond, id.index()),
            DativeBondHandle::New(index) => {
                self.dative_bonds.created(EntityKind::DativeBond, index)
            }
        }
    }

    fn aromatic_system(
        &self,
        handle: AromaticSystemHandle,
    ) -> Result<AromaticSystemId, TransactionError> {
        match handle {
            AromaticSystemHandle::Id(id) => self
                .aromatic_systems
                .initial(EntityKind::AromaticSystem, id.index()),
            AromaticSystemHandle::New(index) => self
                .aromatic_systems
                .created(EntityKind::AromaticSystem, index),
        }
    }

    fn multicenter_bond(
        &self,
        handle: MulticenterBondHandle,
    ) -> Result<MulticenterBondId, TransactionError> {
        match handle {
            MulticenterBondHandle::Id(id) => self
                .multicenter_bonds
                .initial(EntityKind::MulticenterBond, id.index()),
            MulticenterBondHandle::New(index) => self
                .multicenter_bonds
                .created(EntityKind::MulticenterBond, index),
        }
    }

    fn noncovalent_bond(
        &self,
        handle: NoncovalentBondHandle,
    ) -> Result<NoncovalentBondId, TransactionError> {
        match handle {
            NoncovalentBondHandle::Id(id) => self
                .noncovalent_bonds
                .initial(EntityKind::NoncovalentBond, id.index()),
            NoncovalentBondHandle::New(index) => self
                .noncovalent_bonds
                .created(EntityKind::NoncovalentBond, index),
        }
    }

    fn stereo_atom(&self, handle: StereoAtomHandle) -> Result<StereoAtomId, TransactionError> {
        match handle {
            StereoAtomHandle::Id(id) => self
                .stereo_atoms
                .initial(EntityKind::StereoAtom, id.index()),
            StereoAtomHandle::New(index) => {
                self.stereo_atoms.created(EntityKind::StereoAtom, index)
            }
        }
    }

    fn stereo_bond(&self, handle: StereoBondHandle) -> Result<StereoBondId, TransactionError> {
        match handle {
            StereoBondHandle::Id(id) => self
                .stereo_bonds
                .initial(EntityKind::StereoBond, id.index()),
            StereoBondHandle::New(index) => {
                self.stereo_bonds.created(EntityKind::StereoBond, index)
            }
        }
    }

    fn push_atom(&mut self, id: AtomId) {
        self.atoms.push(id);
    }

    fn push_bond(&mut self, id: BondId) {
        self.bonds.push(id);
    }

    fn push_dative_bond(&mut self, id: DativeBondId) {
        self.dative_bonds.push(id);
    }

    fn push_aromatic_system(&mut self, id: AromaticSystemId) {
        self.aromatic_systems.push(id);
    }

    fn push_multicenter_bond(&mut self, id: MulticenterBondId) {
        self.multicenter_bonds.push(id);
    }

    fn push_noncovalent_bond(&mut self, id: NoncovalentBondId) {
        self.noncovalent_bonds.push(id);
    }

    fn push_stereo_atom(&mut self, id: StereoAtomId) {
        self.stereo_atoms.push(id);
    }

    fn push_stereo_bond(&mut self, id: StereoBondId) {
        self.stereo_bonds.push(id);
    }

    fn compact(&mut self, compaction: &MoleculeCompaction) {
        self.atoms.compact(|id| compaction.compact_atom(id));
        self.bonds.compact(|id| compaction.compact_bond(id));
        self.dative_bonds
            .compact(|id| compaction.compact_dative_bond(id));
        self.aromatic_systems
            .compact(|id| compaction.compact_aromatic_system(id));
        self.multicenter_bonds
            .compact(|id| compaction.compact_multicenter_bond(id));
        self.noncovalent_bonds
            .compact(|id| compaction.compact_noncovalent_bond(id));
        self.stereo_atoms
            .compact(|id| compaction.compact_stereo_atom(id));
        self.stereo_bonds
            .compact(|id| compaction.compact_stereo_bond(id));
    }

    /// Update an input correspondence with this batch's removals and surviving additions.
    pub(crate) fn update_correspondence(&self, correspondence: &mut MoleculeCorrespondence) {
        if self.atoms.initial.is_none()
            && self.atoms.created.is_empty()
            && self.bonds.initial.is_none()
            && self.bonds.created.is_empty()
            && self.dative_bonds.initial.is_none()
            && self.dative_bonds.created.is_empty()
            && self.aromatic_systems.initial.is_none()
            && self.aromatic_systems.created.is_empty()
            && self.multicenter_bonds.initial.is_none()
            && self.multicenter_bonds.created.is_empty()
            && self.noncovalent_bonds.initial.is_none()
            && self.noncovalent_bonds.created.is_empty()
            && self.stereo_atoms.initial.is_none()
            && self.stereo_atoms.created.is_empty()
            && self.stereo_bonds.initial.is_none()
            && self.stereo_bonds.created.is_empty()
        {
            return;
        }
        self.update_changed_correspondence(correspondence);
    }

    fn update_changed_correspondence(&self, correspondence: &mut MoleculeCorrespondence) {
        if self.atoms.initial.is_some()
            || self.bonds.initial.is_some()
            || self.dative_bonds.initial.is_some()
            || self.aromatic_systems.initial.is_some()
            || self.multicenter_bonds.initial.is_some()
            || self.noncovalent_bonds.initial.is_some()
            || self.stereo_atoms.initial.is_some()
            || self.stereo_bonds.initial.is_some()
        {
            let compaction = MoleculeCompaction::new(
                GraphCompaction::new(
                    Compaction::new(
                        self.atoms.initial_count,
                        self.atoms.removed().into_iter().map(Into::into).collect(),
                    )
                    .expect("removed handles belong to the batch input"),
                    Compaction::new(
                        self.bonds.initial_count,
                        self.bonds.removed().into_iter().map(Into::into).collect(),
                    )
                    .expect("removed handles belong to the batch input"),
                ),
                Compaction::new(self.dative_bonds.initial_count, self.dative_bonds.removed())
                    .expect("removed handles belong to the batch input"),
                Compaction::new(
                    self.aromatic_systems.initial_count,
                    self.aromatic_systems.removed(),
                )
                .expect("removed handles belong to the batch input"),
                Compaction::new(
                    self.multicenter_bonds.initial_count,
                    self.multicenter_bonds.removed(),
                )
                .expect("removed handles belong to the batch input"),
                Compaction::new(
                    self.noncovalent_bonds.initial_count,
                    self.noncovalent_bonds.removed(),
                )
                .expect("removed handles belong to the batch input"),
                Compaction::new(self.stereo_atoms.initial_count, self.stereo_atoms.removed())
                    .expect("removed handles belong to the batch input"),
                Compaction::new(self.stereo_bonds.initial_count, self.stereo_bonds.removed())
                    .expect("removed handles belong to the batch input"),
            );
            correspondence
                .compact_right(&compaction)
                .expect("correspondence ends in the batch input id spaces");
        }
        for (kind, count) in [
            (
                EntityKind::Atom,
                self.atoms.created.iter().filter(|id| id.is_some()).count(),
            ),
            (
                EntityKind::Bond,
                self.bonds.created.iter().filter(|id| id.is_some()).count(),
            ),
            (
                EntityKind::DativeBond,
                self.dative_bonds
                    .created
                    .iter()
                    .filter(|id| id.is_some())
                    .count(),
            ),
            (
                EntityKind::AromaticSystem,
                self.aromatic_systems
                    .created
                    .iter()
                    .filter(|id| id.is_some())
                    .count(),
            ),
            (
                EntityKind::MulticenterBond,
                self.multicenter_bonds
                    .created
                    .iter()
                    .filter(|id| id.is_some())
                    .count(),
            ),
            (
                EntityKind::NoncovalentBond,
                self.noncovalent_bonds
                    .created
                    .iter()
                    .filter(|id| id.is_some())
                    .count(),
            ),
            (
                EntityKind::StereoAtom,
                self.stereo_atoms
                    .created
                    .iter()
                    .filter(|id| id.is_some())
                    .count(),
            ),
            (
                EntityKind::StereoBond,
                self.stereo_bonds
                    .created
                    .iter()
                    .filter(|id| id.is_some())
                    .count(),
            ),
        ] {
            correspondence.extend_right(kind, count);
        }
    }

    fn stereo_ligands(
        &self,
        ligands: Vec<(AtomHandle, StereoLigandKind)>,
    ) -> Result<Vec<StereoLigand>, TransactionError> {
        ligands
            .into_iter()
            .map(|(atom, kind)| Ok(StereoLigand::new(self.atom(atom)?, kind)))
            .collect()
    }

    fn resolve_constraint(&self, edit: ConstraintEdit) -> Result<Constraint, TransactionError> {
        edit.resolve(
            |handle| self.atom(handle),
            |handle| self.bond(handle),
            |handle| self.dative_bond(handle),
            |handle| self.aromatic_system(handle),
            |handle| self.multicenter_bond(handle),
            |handle| self.noncovalent_bond(handle),
            |handle| self.stereo_atom(handle),
            |handle| self.stereo_bond(handle),
        )
    }
}

impl Molecule {
    pub(crate) fn apply_edit(
        &mut self,
        edit: Edit,
        state: &mut ApplicationState,
    ) -> Result<(), TransactionError> {
        match edit {
            Edit::AddAtoms { atoms } => {
                for id in self.add_atoms(atoms) {
                    state.push_atom(id);
                }
                Ok(())
            }
            Edit::AddBonds { bonds } => {
                let bonds: Vec<_> = bonds
                    .into_iter()
                    .map(
                        |AddBond {
                             endpoints: [first, second],
                             attributes,
                         }| {
                            Ok(([state.atom(first)?, state.atom(second)?], attributes))
                        },
                    )
                    .collect::<Result<_, TransactionError>>()?;
                for id in self.add_bonds(bonds) {
                    state.push_bond(id);
                }
                Ok(())
            }
            Edit::RemoveTopology { atoms, bonds } => {
                let atoms: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|id| {
                        let id = state.atom(id)?;
                        Ok(id)
                    })
                    .collect::<Result<_, _>>()?;
                let bonds: Vec<BondId> = bonds
                    .into_iter()
                    .map(|id| {
                        let id = state.bond(id)?;
                        Ok(id)
                    })
                    .collect::<Result<_, _>>()?;
                ensure_unique(&atoms, EntityKind::Atom)?;
                ensure_unique(&bonds, EntityKind::Bond)?;
                let graph = self.tracked_remove_topology(&atoms, &bonds);
                let dative_bonds = self.compact_dative_bonds(&graph);
                let aromatic_systems = self.compact_aromatic_systems(&graph);
                let multicenter_bonds = self.compact_multicenter_bonds(&graph);
                let noncovalent_bonds = self.compact_noncovalent_bonds(&graph);
                let stereo_atoms = self.compact_stereo_atoms(&graph);
                let stereo_bonds = self.compact_stereo_bonds(&graph);
                let compaction = MoleculeCompaction::new(
                    graph,
                    dative_bonds,
                    aromatic_systems,
                    multicenter_bonds,
                    noncovalent_bonds,
                    stereo_atoms,
                    stereo_bonds,
                );
                self.compact_constraints(&compaction);
                state.compact(&compaction);
                Ok(())
            }
            Edit::ModifyAtomElements { changes } => {
                for (handle, old, new) in changes {
                    let id = state.atom(handle)?;
                    let mut view = self.atom_mut(id);
                    let attributes = view.attributes_mut();
                    if !attributes.element.normalized_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    attributes.element = new;
                }
                Ok(())
            }
            Edit::ModifyAtomIsotopeMasses { changes } => {
                for (handle, old, new) in changes {
                    let id = state.atom(handle)?;
                    let mut view = self.atom_mut(id);
                    let attributes = view.attributes_mut();
                    if !attributes.isotope_mass.normalized_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    attributes.isotope_mass = new;
                }
                Ok(())
            }
            Edit::ModifyAtomCharges { changes } => {
                for (handle, old, new) in changes {
                    let id = state.atom(handle)?;
                    let mut view = self.atom_mut(id);
                    let attributes = view.attributes_mut();
                    if !attributes.charge.normalized_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    attributes.charge = new;
                }
                Ok(())
            }
            Edit::ModifyAtomImplicitHydrogens { changes } => {
                for (handle, old, new) in changes {
                    let id = state.atom(handle)?;
                    let mut view = self.atom_mut(id);
                    let attributes = view.attributes_mut();
                    if !attributes.implicit_hydrogens.normalized_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    attributes.implicit_hydrogens = new;
                }
                Ok(())
            }
            Edit::ModifyAtomLonePairs { changes } => {
                for (handle, old, new) in changes {
                    let id = state.atom(handle)?;
                    let mut view = self.atom_mut(id);
                    let attributes = view.attributes_mut();
                    if !attributes.lone_pairs.normalized_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    attributes.lone_pairs = new;
                }
                Ok(())
            }
            Edit::ModifyAtomUnpairedElectrons { changes } => {
                for (handle, old, new) in changes {
                    let id = state.atom(handle)?;
                    let mut view = self.atom_mut(id);
                    let attributes = view.attributes_mut();
                    if !attributes.unpaired_electrons.normalized_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    attributes.unpaired_electrons = new;
                }
                Ok(())
            }
            Edit::ModifyAtomField { id, change } => {
                let id = state.atom(id)?;
                let mut view = self.atom_mut(id);
                let attributes = view.attributes_mut();
                match change {
                    AtomFieldChange::Element { old, new } => {
                        if !attributes.element.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.element = new;
                    }
                    AtomFieldChange::IsotopeMass { old, new } => {
                        if !attributes.isotope_mass.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.isotope_mass = new;
                    }
                    AtomFieldChange::Charge { old, new } => {
                        if !attributes.charge.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.charge = new;
                    }
                    AtomFieldChange::ImplicitHydrogens { old, new } => {
                        if !attributes.implicit_hydrogens.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.implicit_hydrogens = new;
                    }
                    AtomFieldChange::LonePairs { old, new } => {
                        if !attributes.lone_pairs.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.lone_pairs = new;
                    }
                    AtomFieldChange::UnpairedElectrons { old, new } => {
                        if !attributes.unpaired_electrons.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.unpaired_electrons = new;
                    }
                }
                Ok(())
            }
            Edit::ModifyBondField { id, change } => {
                let id = state.bond(id)?;
                let mut view = self.bond_mut(id);
                let attributes = view.attributes_mut();
                match change {
                    BondFieldChange::Order { old, new } => {
                        if !attributes.order.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.order = new;
                    }
                    BondFieldChange::Charge { old, new } => {
                        if !attributes.charge.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.charge = new;
                    }
                    BondFieldChange::UnpairedElectrons { old, new } => {
                        if !attributes.unpaired_electrons.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.unpaired_electrons = new;
                    }
                }
                Ok(())
            }
            Edit::AddDativeBond {
                donors,
                acceptor,
                attributes,
            } => {
                let donors: Vec<AtomId> = donors
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let acceptor = state.atom(acceptor)?;
                let id = self.add_dative_bond(&donors, acceptor, attributes);
                state.push_dative_bond(id);
                Ok(())
            }
            Edit::RemoveDativeBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, donors, acceptor, attributes) in removes {
                    let id = state.dative_bond(id)?;
                    let donors: Vec<AtomId> = donors
                        .into_iter()
                        .map(|r| state.atom(r))
                        .collect::<Result<_, _>>()?;
                    let acceptor = state.atom(acceptor)?;
                    let set = self.raw_dative_bonds();
                    let old = DativeBondEntry {
                        donors: Cow::Borrowed(&donors),
                        acceptor,
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if !set.is_coincident(id, acceptor, &donors) || !set.entry(id).framed_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::DativeBond)?;
                let dative_bonds = self.tracked_remove_dative_bonds(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    dative_bonds,
                    Compaction::identity(self.aromatic_systems().count()),
                    Compaction::identity(self.multicenter_bonds().count()),
                    Compaction::identity(self.noncovalent_bonds().count()),
                    Compaction::identity(self.stereo_atoms().count()),
                    Compaction::identity(self.stereo_bonds().count()),
                );
                self.compact_constraints(&forward);
                state.compact(&forward);
                Ok(())
            }
            Edit::ReplaceDativeBondDonors { id, old, new } => {
                let id = state.dative_bond(id)?;
                let old = old
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let new = new
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut view = self.dative_bond_view_mut(id);
                if !view.donor_ids().eq(old.iter().copied()) {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_donors(&new);
                Ok(())
            }
            Edit::ReplaceDativeBondAcceptor { id, old, new } => {
                let id = state.dative_bond(id)?;
                let old = state.atom(old)?;
                let new = state.atom(new)?;
                let mut view = self.dative_bond_view_mut(id);
                if view.acceptor_id() != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_acceptor(new);
                Ok(())
            }
            Edit::ModifyDativeBondField { id, change } => {
                let id = state.dative_bond(id)?;
                let mut view = self.dative_bond_mut(id);
                let attributes = view.attributes_mut();
                match change {
                    DativeBondFieldChange::Order { old, new } => {
                        if !attributes.order.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.order = new;
                    }
                }
                Ok(())
            }
            Edit::AddAromaticSystem(addition) => {
                let AddAromaticSystem { atoms, attributes } = *addition;
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_aromatic_system(&resolved, attributes);
                state.push_aromatic_system(id);
                Ok(())
            }
            Edit::RemoveAromaticSystems { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, atoms, attributes) in removes {
                    let id = state.aromatic_system(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    let set = self.raw_aromatic_systems();
                    let old = AromaticSystemEntry {
                        atoms: Cow::Borrowed(&saved_atoms),
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if !set.is_coincident(id, &saved_atoms) || !set.entry(id).framed_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::AromaticSystem)?;
                let aromatic_systems = self.tracked_remove_aromatic_systems(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    aromatic_systems,
                    Compaction::identity(self.multicenter_bonds().count()),
                    Compaction::identity(self.noncovalent_bonds().count()),
                    Compaction::identity(self.stereo_atoms().count()),
                    Compaction::identity(self.stereo_bonds().count()),
                );
                self.compact_constraints(&forward);
                state.compact(&forward);
                Ok(())
            }
            Edit::ReplaceAromaticSystemAtoms { id, old, new } => {
                let id = state.aromatic_system(id)?;
                let old = old
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let new = new
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut view = self.aromatic_system_view_mut(id);
                if !view.atom_ids().eq(old.iter().copied()) {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_atoms(&new);
                Ok(())
            }
            Edit::ModifyAromaticSystemField { id, change } => {
                let id = state.aromatic_system(id)?;
                let mut view = self.aromatic_system_mut(id);
                let attributes = view.attributes_mut();
                match change {
                    AromaticSystemFieldChange::Electrons { old, new } => {
                        if !attributes.electrons.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.electrons = new;
                    }
                    AromaticSystemFieldChange::Charge { old, new } => {
                        if !attributes.charge.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.charge = new;
                    }
                    AromaticSystemFieldChange::UnpairedElectrons { old, new } => {
                        if !attributes.unpaired_electrons.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.unpaired_electrons = new;
                    }
                }
                Ok(())
            }
            Edit::AddMulticenterBond(addition) => {
                let AddMulticenterBond { atoms, attributes } = *addition;
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_multicenter_bond(&resolved, attributes);
                state.push_multicenter_bond(id);
                Ok(())
            }
            Edit::RemoveMulticenterBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, atoms, attributes) in removes {
                    let id = state.multicenter_bond(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    let set = self.raw_multicenter_bonds();
                    let old = MulticenterBondEntry {
                        atoms: Cow::Borrowed(&saved_atoms),
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if !set.is_coincident(id, &saved_atoms) || !set.entry(id).framed_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::MulticenterBond)?;
                let multicenter_bonds = self.tracked_remove_multicenter_bonds(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    Compaction::identity(self.aromatic_systems().count()),
                    multicenter_bonds,
                    Compaction::identity(self.noncovalent_bonds().count()),
                    Compaction::identity(self.stereo_atoms().count()),
                    Compaction::identity(self.stereo_bonds().count()),
                );
                self.compact_constraints(&forward);
                state.compact(&forward);
                Ok(())
            }
            Edit::ReplaceMulticenterBondAtoms { id, old, new } => {
                let id = state.multicenter_bond(id)?;
                let old = old
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let new = new
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut view = self.multicenter_bond_view_mut(id);
                if !view.atom_ids().eq(old.iter().copied()) {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_atoms(&new);
                Ok(())
            }
            Edit::ModifyMulticenterBondField { id, change } => {
                let id = state.multicenter_bond(id)?;
                let mut view = self.multicenter_bond_mut(id);
                let attributes = view.attributes_mut();
                match change {
                    MulticenterBondFieldChange::Electrons { old, new } => {
                        if !attributes.electrons.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.electrons = new;
                    }
                    MulticenterBondFieldChange::Charge { old, new } => {
                        if !attributes.charge.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.charge = new;
                    }
                    MulticenterBondFieldChange::UnpairedElectrons { old, new } => {
                        if !attributes.unpaired_electrons.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.unpaired_electrons = new;
                    }
                }
                Ok(())
            }
            Edit::AddNoncovalentBond { atoms, attributes } => {
                let a = state.atom(atoms[0].clone())?;
                let b = state.atom(atoms[1].clone())?;
                let id = self.add_noncovalent_bond([a, b], attributes);
                state.push_noncovalent_bond(id);
                Ok(())
            }
            Edit::RemoveNoncovalentBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, atoms, attributes) in removes {
                    let id = state.noncovalent_bond(id)?;
                    let saved_atoms =
                        [state.atom(atoms[0].clone())?, state.atom(atoms[1].clone())?];
                    let set = self.raw_noncovalent_bonds();
                    let old = NoncovalentBondEntry {
                        atoms: saved_atoms,
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if !set.is_coincident(id, saved_atoms[0], saved_atoms[1])
                        || !set.entry(id).framed_eq(&old)
                    {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::NoncovalentBond)?;
                let noncovalent_bonds = self.tracked_remove_noncovalent_bonds(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    Compaction::identity(self.aromatic_systems().count()),
                    Compaction::identity(self.multicenter_bonds().count()),
                    noncovalent_bonds,
                    Compaction::identity(self.stereo_atoms().count()),
                    Compaction::identity(self.stereo_bonds().count()),
                );
                self.compact_constraints(&forward);
                state.compact(&forward);
                Ok(())
            }
            Edit::ReplaceNoncovalentBondAtoms { id, old, new } => {
                let id = state.noncovalent_bond(id)?;
                let old = [state.atom(old[0].clone())?, state.atom(old[1].clone())?];
                let new = [state.atom(new[0].clone())?, state.atom(new[1].clone())?];
                let mut view = self.noncovalent_bond_view_mut(id);
                if view.atom_ids() != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_atoms(new);
                Ok(())
            }
            Edit::ModifyNoncovalentBondField { id, change } => {
                let id = state.noncovalent_bond(id)?;
                let mut view = self.noncovalent_bond_mut(id);
                let attributes = view.attributes_mut();
                match change {
                    NoncovalentBondFieldChange::Kind { old, new } => {
                        if !attributes.kind.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.kind = new;
                    }
                }
                Ok(())
            }
            Edit::AddStereoAtom(addition) => {
                let AddStereoAtom {
                    site,
                    ligands,
                    attributes,
                } = *addition;
                let site = state.atom(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_atom(site, &ligands, attributes);
                state.push_stereo_atom(id);
                Ok(())
            }
            Edit::RemoveStereoAtoms { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, site, ligands, attributes) in removes {
                    let id = state.stereo_atom(id)?;
                    let site = state.atom(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    let set = self.raw_stereo_atoms();
                    let old = StereoAtomEntry {
                        site,
                        ligands: Cow::Borrowed(&ligands),
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if set.site(id) != site
                        || ligands.len() > MAX_DEGREE
                        || Permutation::between(&ligands, set.ligands(id)).is_none()
                        || !set.entry(id).framed_eq(&old)
                    {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::StereoAtom)?;
                let stereo_atoms = self.tracked_remove_stereo_atoms(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    Compaction::identity(self.aromatic_systems().count()),
                    Compaction::identity(self.multicenter_bonds().count()),
                    Compaction::identity(self.noncovalent_bonds().count()),
                    stereo_atoms,
                    Compaction::identity(self.stereo_bonds().count()),
                );
                self.compact_constraints(&forward);
                state.compact(&forward);
                Ok(())
            }
            Edit::ReplaceStereoAtomSite { id, old, new } => {
                let id = state.stereo_atom(id)?;
                let old = state.atom(old)?;
                let new = state.atom(new)?;
                let mut view = self.stereo_atom_view_mut(id);
                if view.site_id() != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_site(new);
                Ok(())
            }
            Edit::ReplaceStereoAtomLigands { id, old, new } => {
                let id = state.stereo_atom(id)?;
                let old = state.stereo_ligands(old)?;
                let new = state.stereo_ligands(new)?;
                let mut view = self.stereo_atom_view_mut(id);
                if view.ligand_ids() != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_ligands(&new);
                Ok(())
            }
            Edit::ModifyStereoAtomField { id, change } => {
                let id = state.stereo_atom(id)?;
                let mut view = self.stereo_atom_mut(id);
                let attributes = view.attributes_mut();
                match change {
                    StereoAtomFieldChange::Configuration { old, new } => {
                        if !attributes.configuration.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.configuration = new;
                    }
                }
                Ok(())
            }
            Edit::AddStereoBond(addition) => {
                let AddStereoBond {
                    site,
                    ligands,
                    attributes,
                } = *addition;
                let site = state.bond(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_bond(site, &ligands, attributes);
                state.push_stereo_bond(id);
                Ok(())
            }
            Edit::RemoveStereoBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, site, ligands, attributes) in removes {
                    let id = state.stereo_bond(id)?;
                    let site = state.bond(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    let set = self.raw_stereo_bonds();
                    let old = StereoBondEntry {
                        site,
                        ligands: Cow::Borrowed(&ligands),
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if set.site(id) != site
                        || ligands.len() != 4
                        || !Permutation::between(&ligands, set.ligands(id))
                            .is_some_and(|action| ClassKey::CisTrans.space().allows(action))
                        || !set.entry(id).framed_eq(&old)
                    {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::StereoBond)?;
                let stereo_bonds = self.tracked_remove_stereo_bonds(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    Compaction::identity(self.aromatic_systems().count()),
                    Compaction::identity(self.multicenter_bonds().count()),
                    Compaction::identity(self.noncovalent_bonds().count()),
                    Compaction::identity(self.stereo_atoms().count()),
                    stereo_bonds,
                );
                self.compact_constraints(&forward);
                state.compact(&forward);
                Ok(())
            }
            Edit::ReplaceStereoBondSite { id, old, new } => {
                let id = state.stereo_bond(id)?;
                let old = state.bond(old)?;
                let new = state.bond(new)?;
                let mut view = self.stereo_bond_view_mut(id);
                if view.site_id() != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_site(new);
                Ok(())
            }
            Edit::ReplaceStereoBondLigands { id, old, new } => {
                let id = state.stereo_bond(id)?;
                let old = state.stereo_ligands(old)?;
                let new = state.stereo_ligands(new)?;
                let mut view = self.stereo_bond_view_mut(id);
                if view.ligand_ids() != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_ligands(&new);
                Ok(())
            }
            Edit::ModifyStereoBondField { id, change } => {
                let id = state.stereo_bond(id)?;
                let mut view = self.stereo_bond_mut(id);
                let attributes = view.attributes_mut();
                match change {
                    StereoBondFieldChange::Configuration { old, new } => {
                        if !attributes.configuration.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.configuration = new;
                    }
                }
                Ok(())
            }
            Edit::ModifyAtomConstraints { key, changes } => {
                for (handle, old, new) in changes {
                    let id = state.atom(handle)?;
                    if old.as_ref().is_some_and(|value| value.key() != key)
                        || new.as_ref().is_some_and(|value| value.key() != key)
                    {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    if old.is_none() && new.is_none() {
                        continue;
                    }
                    let mut view = self.atom_mut(id);
                    let constraints = &mut view.attributes_mut().constraints;
                    match (constraints.get(key), old.as_ref()) {
                        (None, None) => {}
                        (Some(current), Some(old)) if current.normalized_eq(old) => {}
                        _ => return Err(TransactionError::OldStateMismatch),
                    }
                    match new {
                        Some(value) => constraints.set(value),
                        None => {
                            constraints.remove(key);
                        }
                    }
                }
                Ok(())
            }
            Edit::ModifyAtomConstraint { id, old, new } => {
                let id = state.atom(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(()),
                };
                let mut view = self.atom_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(())
            }
            Edit::ModifyBondConstraint { id, old, new } => {
                let id = state.bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(()),
                };
                let mut view = self.bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(())
            }
            Edit::ModifyDativeBondConstraint { id, old, new } => {
                let id = state.dative_bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(()),
                };
                let mut view = self.dative_bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(())
            }
            Edit::ModifyAromaticSystemConstraint { id, old, new } => {
                let id = state.aromatic_system(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(()),
                };
                let mut view = self.aromatic_system_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(())
            }
            Edit::ModifyMulticenterBondConstraint { id, old, new } => {
                let id = state.multicenter_bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(()),
                };
                let mut view = self.multicenter_bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(())
            }
            Edit::ModifyNoncovalentBondConstraint { id, old, new } => {
                let id = state.noncovalent_bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(()),
                };
                let mut view = self.noncovalent_bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(())
            }
            Edit::ModifyStereoAtomConstraint {
                id,
                kind: _,
                old,
                new,
            } => {
                let id = state.stereo_atom(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(()),
                };
                let mut view = self.stereo_atom_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(())
            }
            Edit::ModifyStereoBondConstraint {
                id,
                kind: _,
                old,
                new,
            } => {
                let id = state.stereo_bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(()),
                };
                let mut view = self.stereo_bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(())
            }
            Edit::AddMoleculeConstraint { constraint } => {
                let constraint = state.resolve_constraint(*constraint)?;
                self.push_constraint(constraint);
                Ok(())
            }
            Edit::RemoveMoleculeConstraint { constraint } => {
                let constraint = state.resolve_constraint(*constraint)?;
                let positions = self
                    .constraints()
                    .as_slice()
                    .iter()
                    .enumerate()
                    .filter(|(_, stored)| stored.normalized_eq(&constraint))
                    .map(|(position, _)| position)
                    .collect::<Vec<_>>();
                if positions.is_empty() {
                    return Err(TransactionError::MissingEntry);
                }
                for &position in positions.iter().rev() {
                    self.remove_constraint_at(position);
                }
                Ok(())
            }
        }
    }

    pub(crate) fn apply_edit_with_undo(
        &mut self,
        edit: Edit,
        state: &mut ApplicationState,
    ) -> Result<Option<Undo>, TransactionError> {
        match edit {
            Edit::AddAtoms { atoms } => {
                let added = self
                    .add_atoms(atoms)
                    .map(|id| {
                        state.push_atom(id);
                        AddedAtom {
                            id,
                            attributes: self.atom(id).attributes().clone(),
                        }
                    })
                    .collect();
                Ok(Undo::RemoveAddedTopology {
                    atoms: added,
                    bonds: Vec::new(),
                })
            }
            Edit::AddBonds { bonds } => {
                let bonds: Vec<_> = bonds
                    .into_iter()
                    .map(
                        |AddBond {
                             endpoints: [first, second],
                             attributes,
                         }| {
                            Ok(([state.atom(first)?, state.atom(second)?], attributes))
                        },
                    )
                    .collect::<Result<_, TransactionError>>()?;
                let added = self
                    .add_bonds(bonds)
                    .map(|id| {
                        state.push_bond(id);
                        let view = self.bond(id);
                        AddedBond {
                            id,
                            endpoints: view.atom_ids(),
                            attributes: view.attributes().clone(),
                        }
                    })
                    .collect();
                Ok(Undo::RemoveAddedTopology {
                    atoms: Vec::new(),
                    bonds: added,
                })
            }
            Edit::RemoveTopology { atoms, bonds } => {
                let atoms: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|id| {
                        let id = state.atom(id)?;
                        Ok(id)
                    })
                    .collect::<Result<_, _>>()?;
                let bonds: Vec<BondId> = bonds
                    .into_iter()
                    .map(|id| {
                        let id = state.bond(id)?;
                        Ok(id)
                    })
                    .collect::<Result<_, _>>()?;
                ensure_unique(&atoms, EntityKind::Atom)?;
                ensure_unique(&bonds, EntityKind::Bond)?;
                let (removed_atoms, removed_bonds, overlays) =
                    self.capture_removed_topology(&atoms, &bonds);
                let graph = self.tracked_remove_topology(&atoms, &bonds);
                let dative_bonds = self.compact_dative_bonds(&graph);
                let aromatic_systems = self.compact_aromatic_systems(&graph);
                let multicenter_bonds = self.compact_multicenter_bonds(&graph);
                let noncovalent_bonds = self.compact_noncovalent_bonds(&graph);
                let stereo_atoms = self.compact_stereo_atoms(&graph);
                let stereo_bonds = self.compact_stereo_bonds(&graph);
                let compaction = MoleculeCompaction::new(
                    graph,
                    dative_bonds,
                    aromatic_systems,
                    multicenter_bonds,
                    noncovalent_bonds,
                    stereo_atoms,
                    stereo_bonds,
                );
                let cascade = self.tracked_compact_constraints(&compaction);
                state.compact(&compaction);
                let undo_compaction = compaction.undo_compaction();
                Ok(Undo::RestoreRemovedTopology(Box::new(RemovedTopology {
                    atoms: removed_atoms,
                    bonds: removed_bonds,
                    overlays,
                    compaction,
                    undo_compaction,
                    cascade,
                })))
            }
            Edit::ModifyAtomElements { changes } => {
                let mut inverse = Vec::with_capacity(changes.len());
                for (handle, old, new) in changes {
                    let result = state.atom(handle).and_then(|id| {
                        let mut view = self.atom_mut(id);
                        let attributes = view.attributes_mut();
                        if !attributes.element.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.element = new.clone();
                        Ok(id)
                    });
                    match result {
                        Ok(id) => inverse.push((id, new, old)),
                        Err(error) => {
                            for (id, _, old) in inverse.into_iter().rev() {
                                self.atom_mut(id).attributes_mut().element = old;
                            }
                            return Err(error);
                        }
                    }
                }
                Ok(Undo::ModifyAtomElements { changes: inverse })
            }
            Edit::ModifyAtomIsotopeMasses { changes } => {
                let mut inverse = Vec::with_capacity(changes.len());
                for (handle, old, new) in changes {
                    let result = state.atom(handle).and_then(|id| {
                        let mut view = self.atom_mut(id);
                        let attributes = view.attributes_mut();
                        if !attributes.isotope_mass.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.isotope_mass = new.clone();
                        Ok(id)
                    });
                    match result {
                        Ok(id) => inverse.push((id, new, old)),
                        Err(error) => {
                            for (id, _, old) in inverse.into_iter().rev() {
                                self.atom_mut(id).attributes_mut().isotope_mass = old;
                            }
                            return Err(error);
                        }
                    }
                }
                Ok(Undo::ModifyAtomIsotopeMasses { changes: inverse })
            }
            Edit::ModifyAtomCharges { changes } => {
                let mut inverse = Vec::with_capacity(changes.len());
                for (handle, old, new) in changes {
                    let result = state.atom(handle).and_then(|id| {
                        let mut view = self.atom_mut(id);
                        let attributes = view.attributes_mut();
                        if !attributes.charge.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.charge = new.clone();
                        Ok(id)
                    });
                    match result {
                        Ok(id) => inverse.push((id, new, old)),
                        Err(error) => {
                            for (id, _, old) in inverse.into_iter().rev() {
                                self.atom_mut(id).attributes_mut().charge = old;
                            }
                            return Err(error);
                        }
                    }
                }
                Ok(Undo::ModifyAtomCharges { changes: inverse })
            }
            Edit::ModifyAtomImplicitHydrogens { changes } => {
                let mut inverse = Vec::with_capacity(changes.len());
                for (handle, old, new) in changes {
                    let result = state.atom(handle).and_then(|id| {
                        let mut view = self.atom_mut(id);
                        let attributes = view.attributes_mut();
                        if !attributes.implicit_hydrogens.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.implicit_hydrogens = new.clone();
                        Ok(id)
                    });
                    match result {
                        Ok(id) => inverse.push((id, new, old)),
                        Err(error) => {
                            for (id, _, old) in inverse.into_iter().rev() {
                                self.atom_mut(id).attributes_mut().implicit_hydrogens = old;
                            }
                            return Err(error);
                        }
                    }
                }
                Ok(Undo::ModifyAtomImplicitHydrogens { changes: inverse })
            }
            Edit::ModifyAtomLonePairs { changes } => {
                let mut inverse = Vec::with_capacity(changes.len());
                for (handle, old, new) in changes {
                    let result = state.atom(handle).and_then(|id| {
                        let mut view = self.atom_mut(id);
                        let attributes = view.attributes_mut();
                        if !attributes.lone_pairs.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.lone_pairs = new.clone();
                        Ok(id)
                    });
                    match result {
                        Ok(id) => inverse.push((id, new, old)),
                        Err(error) => {
                            for (id, _, old) in inverse.into_iter().rev() {
                                self.atom_mut(id).attributes_mut().lone_pairs = old;
                            }
                            return Err(error);
                        }
                    }
                }
                Ok(Undo::ModifyAtomLonePairs { changes: inverse })
            }
            Edit::ModifyAtomUnpairedElectrons { changes } => {
                let mut inverse = Vec::with_capacity(changes.len());
                for (handle, old, new) in changes {
                    let result = state.atom(handle).and_then(|id| {
                        let mut view = self.atom_mut(id);
                        let attributes = view.attributes_mut();
                        if !attributes.unpaired_electrons.normalized_eq(&old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.unpaired_electrons = new.clone();
                        Ok(id)
                    });
                    match result {
                        Ok(id) => inverse.push((id, new, old)),
                        Err(error) => {
                            for (id, _, old) in inverse.into_iter().rev() {
                                self.atom_mut(id).attributes_mut().unpaired_electrons = old;
                            }
                            return Err(error);
                        }
                    }
                }
                Ok(Undo::ModifyAtomUnpairedElectrons { changes: inverse })
            }
            Edit::ModifyAtomField { id, change } => {
                let id = state.atom(id)?;
                let mut view = self.atom_mut(id);
                let attributes = view.attributes_mut();
                match &change {
                    AtomFieldChange::Element { old, new } => {
                        if !attributes.element.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.element = new.clone();
                    }
                    AtomFieldChange::IsotopeMass { old, new } => {
                        if !attributes.isotope_mass.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.isotope_mass = new.clone();
                    }
                    AtomFieldChange::Charge { old, new } => {
                        if !attributes.charge.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.charge = new.clone();
                    }
                    AtomFieldChange::ImplicitHydrogens { old, new } => {
                        if !attributes.implicit_hydrogens.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.implicit_hydrogens = new.clone();
                    }
                    AtomFieldChange::LonePairs { old, new } => {
                        if !attributes.lone_pairs.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.lone_pairs = new.clone();
                    }
                    AtomFieldChange::UnpairedElectrons { old, new } => {
                        if !attributes.unpaired_electrons.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.unpaired_electrons = new.clone();
                    }
                }
                Ok(Undo::ModifyAtomField {
                    id,
                    change: change.inverse(),
                })
            }
            Edit::ModifyBondField { id, change } => {
                let id = state.bond(id)?;
                let mut view = self.bond_mut(id);
                let attributes = view.attributes_mut();
                match &change {
                    BondFieldChange::Order { old, new } => {
                        if !attributes.order.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.order = new.clone();
                    }
                    BondFieldChange::Charge { old, new } => {
                        if !attributes.charge.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.charge = new.clone();
                    }
                    BondFieldChange::UnpairedElectrons { old, new } => {
                        if !attributes.unpaired_electrons.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.unpaired_electrons = new.clone();
                    }
                }
                Ok(Undo::ModifyBondField {
                    id,
                    change: change.inverse(),
                })
            }
            Edit::AddDativeBond {
                donors,
                acceptor,
                attributes,
            } => {
                let donors: Vec<AtomId> = donors
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let acceptor = state.atom(acceptor)?;
                let id = self.add_dative_bond(&donors, acceptor, attributes);
                state.push_dative_bond(id);
                let view = self.dative_bond(id);
                Ok(Undo::RemoveAddedDativeBond(AddedDativeBond {
                    id,
                    donors: view.donor_ids().collect(),
                    acceptor: view.acceptor_id(),
                    attributes: view.attributes().clone(),
                }))
            }
            Edit::RemoveDativeBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, donors, acceptor, attributes) in removes {
                    let id = state.dative_bond(id)?;
                    let donors: Vec<AtomId> = donors
                        .into_iter()
                        .map(|r| state.atom(r))
                        .collect::<Result<_, _>>()?;
                    let acceptor = state.atom(acceptor)?;
                    let set = self.raw_dative_bonds();
                    let old = DativeBondEntry {
                        donors: Cow::Borrowed(&donors),
                        acceptor,
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if !set.is_coincident(id, acceptor, &donors) || !set.entry(id).framed_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.dative_bond(id);
                    removed.push(RemovedDativeBond {
                        id,
                        donors: view.donor_ids().collect(),
                        acceptor: view.acceptor_id(),
                        attributes: view.attributes().clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::DativeBond)?;
                let dative_bonds = self.tracked_remove_dative_bonds(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    dative_bonds,
                    Compaction::identity(self.aromatic_systems().count()),
                    Compaction::identity(self.multicenter_bonds().count()),
                    Compaction::identity(self.noncovalent_bonds().count()),
                    Compaction::identity(self.stereo_atoms().count()),
                    Compaction::identity(self.stereo_bonds().count()),
                );
                let cascade = self.tracked_compact_constraints(&forward);
                state.compact(&forward);
                Ok(Undo::RestoreRemovedDativeBonds(Box::new(
                    RemovedDativeBonds {
                        removed,
                        undo_compaction: forward.undo_compaction(),
                        cascade,
                    },
                )))
            }
            Edit::ReplaceDativeBondDonors { id, old, new } => {
                let id = state.dative_bond(id)?;
                let old = old
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let new = new
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut view = self.dative_bond_view_mut(id);
                let donors = view.donor_ids().collect::<Vec<_>>();
                if donors != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_donors(&new);
                Ok(Undo::RestoreDativeBondDonors { id, donors })
            }
            Edit::ReplaceDativeBondAcceptor { id, old, new } => {
                let id = state.dative_bond(id)?;
                let old = state.atom(old)?;
                let new = state.atom(new)?;
                let mut view = self.dative_bond_view_mut(id);
                let acceptor = view.acceptor_id();
                if acceptor != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_acceptor(new);
                Ok(Undo::RestoreDativeBondAcceptor { id, acceptor })
            }
            Edit::ModifyDativeBondField { id, change } => {
                let id = state.dative_bond(id)?;
                let mut view = self.dative_bond_mut(id);
                let attributes = view.attributes_mut();
                match &change {
                    DativeBondFieldChange::Order { old, new } => {
                        if !attributes.order.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.order = new.clone();
                    }
                }
                Ok(Undo::ModifyDativeBondField {
                    id,
                    change: change.inverse(),
                })
            }
            Edit::AddAromaticSystem(addition) => {
                let AddAromaticSystem { atoms, attributes } = *addition;
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_aromatic_system(&resolved, attributes);
                state.push_aromatic_system(id);
                let view = self.aromatic_system(id);
                Ok(Undo::RemoveAddedAromaticSystem(Box::new(
                    AddedAromaticSystem {
                        id,
                        atoms: view.atom_ids().collect(),
                        attributes: view.attributes().clone(),
                    },
                )))
            }
            Edit::RemoveAromaticSystems { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, atoms, attributes) in removes {
                    let id = state.aromatic_system(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    let set = self.raw_aromatic_systems();
                    let old = AromaticSystemEntry {
                        atoms: Cow::Borrowed(&saved_atoms),
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if !set.is_coincident(id, &saved_atoms) || !set.entry(id).framed_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.aromatic_system(id);
                    let current_atoms: Vec<AtomId> = view.atom_ids().collect();
                    removed.push(RemovedAromaticSystem {
                        id,
                        atoms: current_atoms,
                        attributes: view.attributes().clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::AromaticSystem)?;
                let aromatic_systems = self.tracked_remove_aromatic_systems(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    aromatic_systems,
                    Compaction::identity(self.multicenter_bonds().count()),
                    Compaction::identity(self.noncovalent_bonds().count()),
                    Compaction::identity(self.stereo_atoms().count()),
                    Compaction::identity(self.stereo_bonds().count()),
                );
                let cascade = self.tracked_compact_constraints(&forward);
                state.compact(&forward);
                Ok(Undo::RestoreRemovedAromaticSystems(Box::new(
                    RemovedAromaticSystems {
                        removed,
                        undo_compaction: forward.undo_compaction(),
                        cascade,
                    },
                )))
            }
            Edit::ReplaceAromaticSystemAtoms { id, old, new } => {
                let id = state.aromatic_system(id)?;
                let old = old
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let new = new
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut view = self.aromatic_system_view_mut(id);
                let atoms = view.atom_ids().collect::<Vec<_>>();
                if atoms != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_atoms(&new);
                Ok(Undo::RestoreAromaticSystemAtoms { id, atoms })
            }
            Edit::ModifyAromaticSystemField { id, change } => {
                let id = state.aromatic_system(id)?;
                let mut view = self.aromatic_system_mut(id);
                let attributes = view.attributes_mut();
                match &change {
                    AromaticSystemFieldChange::Electrons { old, new } => {
                        if !attributes.electrons.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.electrons = new.clone();
                    }
                    AromaticSystemFieldChange::Charge { old, new } => {
                        if !attributes.charge.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.charge = new.clone();
                    }
                    AromaticSystemFieldChange::UnpairedElectrons { old, new } => {
                        if !attributes.unpaired_electrons.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.unpaired_electrons = new.clone();
                    }
                }
                Ok(Undo::ModifyAromaticSystemField {
                    id,
                    change: change.inverse(),
                })
            }
            Edit::AddMulticenterBond(addition) => {
                let AddMulticenterBond { atoms, attributes } = *addition;
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_multicenter_bond(&resolved, attributes);
                state.push_multicenter_bond(id);
                let view = self.multicenter_bond(id);
                Ok(Undo::RemoveAddedMulticenterBond(Box::new(
                    AddedMulticenterBond {
                        id,
                        atoms: view.atom_ids().collect(),
                        attributes: view.attributes().clone(),
                    },
                )))
            }
            Edit::RemoveMulticenterBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, atoms, attributes) in removes {
                    let id = state.multicenter_bond(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    let set = self.raw_multicenter_bonds();
                    let old = MulticenterBondEntry {
                        atoms: Cow::Borrowed(&saved_atoms),
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if !set.is_coincident(id, &saved_atoms) || !set.entry(id).framed_eq(&old) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.multicenter_bond(id);
                    let current_atoms: Vec<AtomId> = view.atom_ids().collect();
                    removed.push(RemovedMulticenterBond {
                        id,
                        atoms: current_atoms,
                        attributes: view.attributes().clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::MulticenterBond)?;
                let multicenter_bonds = self.tracked_remove_multicenter_bonds(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    Compaction::identity(self.aromatic_systems().count()),
                    multicenter_bonds,
                    Compaction::identity(self.noncovalent_bonds().count()),
                    Compaction::identity(self.stereo_atoms().count()),
                    Compaction::identity(self.stereo_bonds().count()),
                );
                let cascade = self.tracked_compact_constraints(&forward);
                state.compact(&forward);
                Ok(Undo::RestoreRemovedMulticenterBonds(Box::new(
                    RemovedMulticenterBonds {
                        removed,
                        undo_compaction: forward.undo_compaction(),
                        cascade,
                    },
                )))
            }
            Edit::ReplaceMulticenterBondAtoms { id, old, new } => {
                let id = state.multicenter_bond(id)?;
                let old = old
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let new = new
                    .into_iter()
                    .map(|atom| state.atom(atom))
                    .collect::<Result<Vec<_>, _>>()?;
                let mut view = self.multicenter_bond_view_mut(id);
                let atoms = view.atom_ids().collect::<Vec<_>>();
                if atoms != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_atoms(&new);
                Ok(Undo::RestoreMulticenterBondAtoms { id, atoms })
            }
            Edit::ModifyMulticenterBondField { id, change } => {
                let id = state.multicenter_bond(id)?;
                let mut view = self.multicenter_bond_mut(id);
                let attributes = view.attributes_mut();
                match &change {
                    MulticenterBondFieldChange::Electrons { old, new } => {
                        if !attributes.electrons.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.electrons = new.clone();
                    }
                    MulticenterBondFieldChange::Charge { old, new } => {
                        if !attributes.charge.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.charge = new.clone();
                    }
                    MulticenterBondFieldChange::UnpairedElectrons { old, new } => {
                        if !attributes.unpaired_electrons.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.unpaired_electrons = new.clone();
                    }
                }
                Ok(Undo::ModifyMulticenterBondField {
                    id,
                    change: change.inverse(),
                })
            }
            Edit::AddNoncovalentBond { atoms, attributes } => {
                let a = state.atom(atoms[0].clone())?;
                let b = state.atom(atoms[1].clone())?;
                let id = self.add_noncovalent_bond([a, b], attributes);
                state.push_noncovalent_bond(id);
                let view = self.noncovalent_bond(id);
                Ok(Undo::RemoveAddedNoncovalentBond(AddedNoncovalentBond {
                    id,
                    atoms: view.atom_ids(),
                    attributes: view.attributes().clone(),
                }))
            }
            Edit::RemoveNoncovalentBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, atoms, attributes) in removes {
                    let id = state.noncovalent_bond(id)?;
                    let saved_atoms =
                        [state.atom(atoms[0].clone())?, state.atom(atoms[1].clone())?];
                    let set = self.raw_noncovalent_bonds();
                    let old = NoncovalentBondEntry {
                        atoms: saved_atoms,
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if !set.is_coincident(id, saved_atoms[0], saved_atoms[1])
                        || !set.entry(id).framed_eq(&old)
                    {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.noncovalent_bond(id);
                    removed.push(RemovedNoncovalentBond {
                        id,
                        atoms: view.atom_ids(),
                        attributes: view.attributes().clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::NoncovalentBond)?;
                let noncovalent_bonds = self.tracked_remove_noncovalent_bonds(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    Compaction::identity(self.aromatic_systems().count()),
                    Compaction::identity(self.multicenter_bonds().count()),
                    noncovalent_bonds,
                    Compaction::identity(self.stereo_atoms().count()),
                    Compaction::identity(self.stereo_bonds().count()),
                );
                let cascade = self.tracked_compact_constraints(&forward);
                state.compact(&forward);
                Ok(Undo::RestoreRemovedNoncovalentBonds(Box::new(
                    RemovedNoncovalentBonds {
                        removed,
                        undo_compaction: forward.undo_compaction(),
                        cascade,
                    },
                )))
            }
            Edit::ReplaceNoncovalentBondAtoms { id, old, new } => {
                let id = state.noncovalent_bond(id)?;
                let old = [state.atom(old[0].clone())?, state.atom(old[1].clone())?];
                let new = [state.atom(new[0].clone())?, state.atom(new[1].clone())?];
                let mut view = self.noncovalent_bond_view_mut(id);
                let atoms = view.atom_ids();
                if atoms != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_atoms(new);
                Ok(Undo::RestoreNoncovalentBondAtoms { id, atoms })
            }
            Edit::ModifyNoncovalentBondField { id, change } => {
                let id = state.noncovalent_bond(id)?;
                let mut view = self.noncovalent_bond_mut(id);
                let attributes = view.attributes_mut();
                match &change {
                    NoncovalentBondFieldChange::Kind { old, new } => {
                        if !attributes.kind.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.kind = new.clone();
                    }
                }
                Ok(Undo::ModifyNoncovalentBondField {
                    id,
                    change: change.inverse(),
                })
            }
            Edit::AddStereoAtom(addition) => {
                let AddStereoAtom {
                    site,
                    ligands,
                    attributes,
                } = *addition;
                let site = state.atom(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_atom(site, &ligands, attributes);
                state.push_stereo_atom(id);
                let view = self.stereo_atom(id);
                Ok(Undo::RemoveAddedStereoAtom(Box::new(AddedStereoAtom {
                    id,
                    site,
                    ligands,
                    attributes: view.attributes().clone(),
                })))
            }
            Edit::RemoveStereoAtoms { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, site, ligands, attributes) in removes {
                    let id = state.stereo_atom(id)?;
                    let site = state.atom(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    let set = self.raw_stereo_atoms();
                    let old = StereoAtomEntry {
                        site,
                        ligands: Cow::Borrowed(&ligands),
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if set.site(id) != site
                        || ligands.len() > MAX_DEGREE
                        || Permutation::between(&ligands, set.ligands(id)).is_none()
                        || !set.entry(id).framed_eq(&old)
                    {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.stereo_atom(id);
                    removed.push(RemovedStereoAtom {
                        id,
                        site: view.site_id(),
                        ligands: view.ligand_ids().to_vec(),
                        attributes: view.attributes().clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::StereoAtom)?;
                let stereo_atoms = self.tracked_remove_stereo_atoms(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    Compaction::identity(self.aromatic_systems().count()),
                    Compaction::identity(self.multicenter_bonds().count()),
                    Compaction::identity(self.noncovalent_bonds().count()),
                    stereo_atoms,
                    Compaction::identity(self.stereo_bonds().count()),
                );
                let cascade = self.tracked_compact_constraints(&forward);
                state.compact(&forward);
                Ok(Undo::RestoreRemovedStereoAtoms(Box::new(
                    RemovedStereoAtoms {
                        removed,
                        undo_compaction: forward.undo_compaction(),
                        cascade,
                    },
                )))
            }
            Edit::ReplaceStereoAtomSite { id, old, new } => {
                let id = state.stereo_atom(id)?;
                let old = state.atom(old)?;
                let new = state.atom(new)?;
                let mut view = self.stereo_atom_view_mut(id);
                let site = view.site_id();
                if site != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_site(new);
                Ok(Undo::RestoreStereoAtomSite { id, site })
            }
            Edit::ReplaceStereoAtomLigands { id, old, new } => {
                let id = state.stereo_atom(id)?;
                let old = state.stereo_ligands(old)?;
                let new = state.stereo_ligands(new)?;
                let mut view = self.stereo_atom_view_mut(id);
                let ligands = view.ligand_ids().to_vec();
                if ligands != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_ligands(&new);
                Ok(Undo::RestoreStereoAtomLigands { id, ligands })
            }
            Edit::ModifyStereoAtomField { id, change } => {
                let id = state.stereo_atom(id)?;
                let mut view = self.stereo_atom_mut(id);
                let attributes = view.attributes_mut();
                match &change {
                    StereoAtomFieldChange::Configuration { old, new } => {
                        if !attributes.configuration.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.configuration = new.clone();
                    }
                }
                Ok(Undo::ModifyStereoAtomField {
                    id,
                    change: change.inverse(),
                })
            }
            Edit::AddStereoBond(addition) => {
                let AddStereoBond {
                    site,
                    ligands,
                    attributes,
                } = *addition;
                let site = state.bond(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_bond(site, &ligands, attributes);
                state.push_stereo_bond(id);
                let view = self.stereo_bond(id);
                Ok(Undo::RemoveAddedStereoBond(Box::new(AddedStereoBond {
                    id,
                    site,
                    ligands,
                    attributes: view.attributes().clone(),
                })))
            }
            Edit::RemoveStereoBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, site, ligands, attributes) in removes {
                    let id = state.stereo_bond(id)?;
                    let site = state.bond(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    let set = self.raw_stereo_bonds();
                    let old = StereoBondEntry {
                        site,
                        ligands: Cow::Borrowed(&ligands),
                        attributes: Cow::Borrowed(&attributes),
                    };
                    if set.site(id) != site
                        || ligands.len() != 4
                        || !Permutation::between(&ligands, set.ligands(id))
                            .is_some_and(|action| ClassKey::CisTrans.space().allows(action))
                        || !set.entry(id).framed_eq(&old)
                    {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.stereo_bond(id);
                    removed.push(RemovedStereoBond {
                        id,
                        site: view.site_id(),
                        ligands: view.ligand_ids().to_vec(),
                        attributes: view.attributes().clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::StereoBond)?;
                let stereo_bonds = self.tracked_remove_stereo_bonds(&ids);
                let forward = MoleculeCompaction::new(
                    GraphCompaction::new(
                        Compaction::identity(self.atoms().count()),
                        Compaction::identity(self.bonds().count()),
                    ),
                    Compaction::identity(self.dative_bonds().count()),
                    Compaction::identity(self.aromatic_systems().count()),
                    Compaction::identity(self.multicenter_bonds().count()),
                    Compaction::identity(self.noncovalent_bonds().count()),
                    Compaction::identity(self.stereo_atoms().count()),
                    stereo_bonds,
                );
                let cascade = self.tracked_compact_constraints(&forward);
                state.compact(&forward);
                Ok(Undo::RestoreRemovedStereoBonds(Box::new(
                    RemovedStereoBonds {
                        removed,
                        undo_compaction: forward.undo_compaction(),
                        cascade,
                    },
                )))
            }
            Edit::ReplaceStereoBondSite { id, old, new } => {
                let id = state.stereo_bond(id)?;
                let old = state.bond(old)?;
                let new = state.bond(new)?;
                let mut view = self.stereo_bond_view_mut(id);
                let site = view.site_id();
                if site != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_site(new);
                Ok(Undo::RestoreStereoBondSite { id, site })
            }
            Edit::ReplaceStereoBondLigands { id, old, new } => {
                let id = state.stereo_bond(id)?;
                let old = state.stereo_ligands(old)?;
                let new = state.stereo_ligands(new)?;
                let mut view = self.stereo_bond_view_mut(id);
                let ligands = view.ligand_ids().to_vec();
                if ligands != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                view.replace_ligands(&new);
                Ok(Undo::RestoreStereoBondLigands { id, ligands })
            }
            Edit::ModifyStereoBondField { id, change } => {
                let id = state.stereo_bond(id)?;
                let mut view = self.stereo_bond_mut(id);
                let attributes = view.attributes_mut();
                match &change {
                    StereoBondFieldChange::Configuration { old, new } => {
                        if !attributes.configuration.normalized_eq(old) {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        attributes.configuration = new.clone();
                    }
                }
                Ok(Undo::ModifyStereoBondField {
                    id,
                    change: change.inverse(),
                })
            }
            Edit::ModifyAtomConstraints { key, changes } => {
                let mut inverse = Vec::with_capacity(changes.len());
                for (handle, old, new) in changes {
                    let result = (|| {
                        let id = state.atom(handle)?;
                        if old.as_ref().is_some_and(|value| value.key() != key)
                            || new.as_ref().is_some_and(|value| value.key() != key)
                        {
                            return Err(TransactionError::OldStateMismatch);
                        }
                        if old.is_none() && new.is_none() {
                            return Ok(None);
                        }
                        let mut view = self.atom_mut(id);
                        let constraints = &mut view.attributes_mut().constraints;
                        match (constraints.get(key), old.as_ref()) {
                            (None, None) => {}
                            (Some(current), Some(old)) if current.normalized_eq(old) => {}
                            _ => return Err(TransactionError::OldStateMismatch),
                        }
                        match new {
                            Some(value) => constraints.set(value),
                            None => {
                                constraints.remove(key);
                            }
                        }
                        Ok(Some(id))
                    })();
                    match result {
                        Ok(Some(id)) => inverse.push((id, old)),
                        Ok(None) => {}
                        Err(error) => {
                            for (id, old) in inverse.into_iter().rev() {
                                let mut view = self.atom_mut(id);
                                let constraints = &mut view.attributes_mut().constraints;
                                match old {
                                    Some(value) => constraints.set(value),
                                    None => {
                                        constraints.remove(key);
                                    }
                                }
                            }
                            return Err(error);
                        }
                    }
                }
                Ok(Undo::RestoreAtomConstraints {
                    key,
                    constraints: inverse,
                })
            }
            Edit::ModifyAtomConstraint { id, old, new } => {
                let id = state.atom(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(None),
                };
                let mut view = self.atom_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(Undo::RestoreAtomConstraint {
                    id,
                    key,
                    constraint: old,
                })
            }
            Edit::ModifyBondConstraint { id, old, new } => {
                let id = state.bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(None),
                };
                let mut view = self.bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(Undo::RestoreBondConstraint {
                    id,
                    key,
                    constraint: old,
                })
            }
            Edit::ModifyDativeBondConstraint { id, old, new } => {
                let id = state.dative_bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(None),
                };
                let mut view = self.dative_bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(Undo::RestoreDativeBondConstraint {
                    id,
                    key,
                    constraint: old,
                })
            }
            Edit::ModifyAromaticSystemConstraint { id, old, new } => {
                let id = state.aromatic_system(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(None),
                };
                let mut view = self.aromatic_system_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(Undo::RestoreAromaticSystemConstraint {
                    id,
                    key,
                    constraint: old,
                })
            }
            Edit::ModifyMulticenterBondConstraint { id, old, new } => {
                let id = state.multicenter_bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(None),
                };
                let mut view = self.multicenter_bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(Undo::RestoreMulticenterBondConstraint {
                    id,
                    key,
                    constraint: old,
                })
            }
            Edit::ModifyNoncovalentBondConstraint { id, old, new } => {
                let id = state.noncovalent_bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(None),
                };
                let mut view = self.noncovalent_bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(Undo::RestoreNoncovalentBondConstraint {
                    id,
                    key,
                    constraint: old,
                })
            }
            Edit::ModifyStereoAtomConstraint {
                id,
                kind: _,
                old,
                new,
            } => {
                let id = state.stereo_atom(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(None),
                };
                let mut view = self.stereo_atom_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(Undo::RestoreStereoAtomConstraint {
                    id,
                    key,
                    constraint: old,
                })
            }
            Edit::ModifyStereoBondConstraint {
                id,
                kind: _,
                old,
                new,
            } => {
                let id = state.stereo_bond(id)?;
                let key = match (&old, &new) {
                    (Some(old), Some(new)) if old.key() != new.key() => {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    (Some(constraint), _) | (_, Some(constraint)) => constraint.key(),
                    (None, None) => return Ok(None),
                };
                let mut view = self.stereo_bond_mut(id);
                let constraints = &mut view.attributes_mut().constraints;
                match (constraints.get(key), old.as_ref()) {
                    (None, None) => {}
                    (Some(current), Some(old)) if current.normalized_eq(old) => {}
                    _ => return Err(TransactionError::OldStateMismatch),
                }
                match new {
                    Some(constraint) => constraints.set(constraint),
                    None => {
                        constraints.remove(key);
                    }
                }
                Ok(Undo::RestoreStereoBondConstraint {
                    id,
                    key,
                    constraint: old,
                })
            }
            Edit::AddMoleculeConstraint { constraint } => {
                let constraint = state.resolve_constraint(*constraint)?;
                let position = self.constraints().len();
                self.push_constraint(constraint);
                Ok(Undo::RemoveAddedMoleculeConstraint { position })
            }
            Edit::RemoveMoleculeConstraint { constraint } => {
                let constraint = state.resolve_constraint(*constraint)?;
                let positions = self
                    .constraints()
                    .as_slice()
                    .iter()
                    .enumerate()
                    .filter(|(_, stored)| stored.normalized_eq(&constraint))
                    .map(|(position, _)| position)
                    .collect::<Vec<_>>();
                if positions.is_empty() {
                    return Err(TransactionError::MissingEntry);
                }
                let mut removed = positions
                    .into_iter()
                    .rev()
                    .map(|position| RemovedConstraint {
                        position,
                        constraint: self.remove_constraint_at(position),
                    })
                    .collect::<Vec<_>>();
                removed.reverse();
                Ok(Undo::RestoreMoleculeConstraints(CascadedConstraints {
                    removed,
                    modified: Vec::new(),
                }))
            }
        }
        .map(Some)
    }

    pub(crate) fn apply_undo(&mut self, undo: Undo) {
        match undo {
            Undo::RemoveAddedTopology { atoms, bonds } => {
                let atoms = atoms
                    .into_iter()
                    .map(|entry| entry.id)
                    .filter(|id| id.index() < self.atoms().count())
                    .collect::<Vec<_>>();
                let bonds = bonds
                    .into_iter()
                    .map(|entry| entry.id)
                    .filter(|id| id.index() < self.bonds().count())
                    .collect::<Vec<_>>();
                self.remove_topology(&atoms, &bonds);
            }
            Undo::RestoreRemovedTopology(removed) => {
                let RemovedTopology {
                    atoms,
                    bonds,
                    overlays,
                    undo_compaction,
                    cascade,
                    ..
                } = *removed;
                let compaction = undo_compaction.forward();
                self.restore_topology(compaction.graph(), atoms, bonds);
                self.restore_dative_bond_topology_ids(compaction.graph());
                self.restore_dative_bonds(
                    compaction.dative_bonds(),
                    overlays
                        .dative_bonds
                        .into_iter()
                        .map(|entry| (entry.id, entry.donors, entry.acceptor, entry.attributes))
                        .collect(),
                );
                self.restore_aromatic_system_topology_ids(compaction.graph());
                self.restore_aromatic_systems(
                    compaction.aromatic_systems(),
                    overlays
                        .aromatic_systems
                        .into_iter()
                        .map(|entry| (entry.id, entry.atoms, entry.attributes))
                        .collect(),
                );
                self.restore_multicenter_bond_topology_ids(compaction.graph());
                self.restore_multicenter_bonds(
                    compaction.multicenter_bonds(),
                    overlays
                        .multicenter_bonds
                        .into_iter()
                        .map(|entry| (entry.id, entry.atoms, entry.attributes))
                        .collect(),
                );
                self.restore_noncovalent_bond_topology_ids(compaction.graph());
                self.restore_noncovalent_bonds(
                    compaction.noncovalent_bonds(),
                    overlays
                        .noncovalent_bonds
                        .into_iter()
                        .map(|entry| (entry.id, entry.atoms, entry.attributes))
                        .collect(),
                );
                self.restore_stereo_atom_topology_ids(compaction.graph());
                self.restore_stereo_atoms(
                    compaction.stereo_atoms(),
                    overlays
                        .stereo_atoms
                        .into_iter()
                        .map(|entry| (entry.id, entry.site, entry.ligands, entry.attributes))
                        .collect(),
                );
                self.restore_stereo_bond_topology_ids(compaction.graph());
                self.restore_stereo_bonds(
                    compaction.stereo_bonds(),
                    overlays
                        .stereo_bonds
                        .into_iter()
                        .map(|entry| (entry.id, entry.site, entry.ligands, entry.attributes))
                        .collect(),
                );
                self.restore_constraints(&cascade);
            }
            Undo::RemoveAddedDativeBond(added) => {
                if added.id.index() < self.dative_bonds().count() {
                    self.remove_dative_bonds(&[added.id]);
                }
            }
            Undo::RestoreRemovedDativeBonds(removed) => {
                let RemovedDativeBonds {
                    removed,
                    undo_compaction,
                    cascade,
                } = *removed;
                self.restore_dative_bonds(
                    undo_compaction.forward().dative_bonds(),
                    removed
                        .into_iter()
                        .map(|entry| (entry.id, entry.donors, entry.acceptor, entry.attributes))
                        .collect(),
                );
                self.restore_constraints(&cascade);
            }
            Undo::RestoreDativeBondDonors { id, donors } => {
                if id.index() < self.dative_bonds().count() {
                    self.dative_bond_view_mut(id).replace_donors(&donors);
                }
            }
            Undo::RestoreDativeBondAcceptor { id, acceptor } => {
                if id.index() < self.dative_bonds().count() {
                    self.dative_bond_view_mut(id).replace_acceptor(acceptor);
                }
            }
            Undo::RemoveAddedAromaticSystem(added) => {
                if added.id.index() < self.aromatic_systems().count() {
                    self.remove_aromatic_systems(&[added.id]);
                }
            }
            Undo::RestoreRemovedAromaticSystems(removed) => {
                let RemovedAromaticSystems {
                    removed,
                    undo_compaction,
                    cascade,
                } = *removed;
                self.restore_aromatic_systems(
                    undo_compaction.forward().aromatic_systems(),
                    removed
                        .into_iter()
                        .map(|entry| (entry.id, entry.atoms, entry.attributes))
                        .collect(),
                );
                self.restore_constraints(&cascade);
            }
            Undo::RestoreAromaticSystemAtoms { id, atoms } => {
                if id.index() < self.aromatic_systems().count() {
                    self.aromatic_system_view_mut(id).replace_atoms(&atoms);
                }
            }
            Undo::RemoveAddedMulticenterBond(added) => {
                if added.id.index() < self.multicenter_bonds().count() {
                    self.remove_multicenter_bonds(&[added.id]);
                }
            }
            Undo::RestoreRemovedMulticenterBonds(removed) => {
                let RemovedMulticenterBonds {
                    removed,
                    undo_compaction,
                    cascade,
                } = *removed;
                self.restore_multicenter_bonds(
                    undo_compaction.forward().multicenter_bonds(),
                    removed
                        .into_iter()
                        .map(|entry| (entry.id, entry.atoms, entry.attributes))
                        .collect(),
                );
                self.restore_constraints(&cascade);
            }
            Undo::RestoreMulticenterBondAtoms { id, atoms } => {
                if id.index() < self.multicenter_bonds().count() {
                    self.multicenter_bond_view_mut(id).replace_atoms(&atoms);
                }
            }
            Undo::RemoveAddedNoncovalentBond(added) => {
                if added.id.index() < self.noncovalent_bonds().count() {
                    self.remove_noncovalent_bonds(&[added.id]);
                }
            }
            Undo::RestoreRemovedNoncovalentBonds(removed) => {
                let RemovedNoncovalentBonds {
                    removed,
                    undo_compaction,
                    cascade,
                } = *removed;
                self.restore_noncovalent_bonds(
                    undo_compaction.forward().noncovalent_bonds(),
                    removed
                        .into_iter()
                        .map(|entry| (entry.id, entry.atoms, entry.attributes))
                        .collect(),
                );
                self.restore_constraints(&cascade);
            }
            Undo::RestoreNoncovalentBondAtoms { id, atoms } => {
                if id.index() < self.noncovalent_bonds().count() {
                    self.noncovalent_bond_view_mut(id).replace_atoms(atoms);
                }
            }
            Undo::RemoveAddedStereoAtom(added) => {
                if added.id.index() < self.stereo_atoms().count() {
                    self.remove_stereo_atoms(&[added.id]);
                }
            }
            Undo::RestoreRemovedStereoAtoms(removed) => {
                let RemovedStereoAtoms {
                    removed,
                    undo_compaction,
                    cascade,
                } = *removed;
                self.restore_stereo_atoms(
                    undo_compaction.forward().stereo_atoms(),
                    removed
                        .into_iter()
                        .map(|entry| (entry.id, entry.site, entry.ligands, entry.attributes))
                        .collect(),
                );
                self.restore_constraints(&cascade);
            }
            Undo::RestoreStereoAtomSite { id, site } => {
                if id.index() < self.stereo_atoms().count() {
                    self.stereo_atom_view_mut(id).replace_site(site);
                }
            }
            Undo::RestoreStereoAtomLigands { id, ligands } => {
                if id.index() < self.stereo_atoms().count() {
                    self.stereo_atom_view_mut(id).replace_ligands(&ligands);
                }
            }
            Undo::RemoveAddedStereoBond(added) => {
                if added.id.index() < self.stereo_bonds().count() {
                    self.remove_stereo_bonds(&[added.id]);
                }
            }
            Undo::RestoreRemovedStereoBonds(removed) => {
                let RemovedStereoBonds {
                    removed,
                    undo_compaction,
                    cascade,
                } = *removed;
                self.restore_stereo_bonds(
                    undo_compaction.forward().stereo_bonds(),
                    removed
                        .into_iter()
                        .map(|entry| (entry.id, entry.site, entry.ligands, entry.attributes))
                        .collect(),
                );
                self.restore_constraints(&cascade);
            }
            Undo::RestoreStereoBondSite { id, site } => {
                if id.index() < self.stereo_bonds().count() {
                    self.stereo_bond_view_mut(id).replace_site(site);
                }
            }
            Undo::RestoreStereoBondLigands { id, ligands } => {
                if id.index() < self.stereo_bonds().count() {
                    self.stereo_bond_view_mut(id).replace_ligands(&ligands);
                }
            }
            Undo::ModifyAtomElements { changes } => {
                for (id, _, new) in changes.into_iter().rev() {
                    if id.index() < self.atoms().count() {
                        self.atom_mut(id).attributes_mut().element = new;
                    }
                }
            }
            Undo::ModifyAtomIsotopeMasses { changes } => {
                for (id, _, new) in changes.into_iter().rev() {
                    if id.index() < self.atoms().count() {
                        self.atom_mut(id).attributes_mut().isotope_mass = new;
                    }
                }
            }
            Undo::ModifyAtomCharges { changes } => {
                for (id, _, new) in changes.into_iter().rev() {
                    if id.index() < self.atoms().count() {
                        self.atom_mut(id).attributes_mut().charge = new;
                    }
                }
            }
            Undo::ModifyAtomImplicitHydrogens { changes } => {
                for (id, _, new) in changes.into_iter().rev() {
                    if id.index() < self.atoms().count() {
                        self.atom_mut(id).attributes_mut().implicit_hydrogens = new;
                    }
                }
            }
            Undo::ModifyAtomLonePairs { changes } => {
                for (id, _, new) in changes.into_iter().rev() {
                    if id.index() < self.atoms().count() {
                        self.atom_mut(id).attributes_mut().lone_pairs = new;
                    }
                }
            }
            Undo::ModifyAtomUnpairedElectrons { changes } => {
                for (id, _, new) in changes.into_iter().rev() {
                    if id.index() < self.atoms().count() {
                        self.atom_mut(id).attributes_mut().unpaired_electrons = new;
                    }
                }
            }
            Undo::ModifyAtomField { id, change } => {
                if id.index() < self.atoms().count() {
                    let mut view = self.atom_mut(id);
                    let attributes = view.attributes_mut();
                    match change {
                        AtomFieldChange::Element { new, .. } => {
                            attributes.element = new;
                        }
                        AtomFieldChange::IsotopeMass { new, .. } => {
                            attributes.isotope_mass = new;
                        }
                        AtomFieldChange::Charge { new, .. } => {
                            attributes.charge = new;
                        }
                        AtomFieldChange::ImplicitHydrogens { new, .. } => {
                            attributes.implicit_hydrogens = new;
                        }
                        AtomFieldChange::LonePairs { new, .. } => {
                            attributes.lone_pairs = new;
                        }
                        AtomFieldChange::UnpairedElectrons { new, .. } => {
                            attributes.unpaired_electrons = new;
                        }
                    }
                }
            }
            Undo::ModifyBondField { id, change } => {
                if id.index() < self.bonds().count() {
                    let mut view = self.bond_mut(id);
                    let attributes = view.attributes_mut();
                    match change {
                        BondFieldChange::Order { new, .. } => {
                            attributes.order = new;
                        }
                        BondFieldChange::Charge { new, .. } => {
                            attributes.charge = new;
                        }
                        BondFieldChange::UnpairedElectrons { new, .. } => {
                            attributes.unpaired_electrons = new;
                        }
                    }
                }
            }
            Undo::ModifyDativeBondField { id, change } => {
                if id.index() < self.dative_bonds().count() {
                    let mut view = self.dative_bond_mut(id);
                    let attributes = view.attributes_mut();
                    match change {
                        DativeBondFieldChange::Order { new, .. } => {
                            attributes.order = new;
                        }
                    }
                }
            }
            Undo::ModifyAromaticSystemField { id, change } => {
                if id.index() < self.aromatic_systems().count() {
                    let mut view = self.aromatic_system_mut(id);
                    let attributes = view.attributes_mut();
                    match change {
                        AromaticSystemFieldChange::Electrons { new, .. } => {
                            attributes.electrons = new;
                        }
                        AromaticSystemFieldChange::Charge { new, .. } => {
                            attributes.charge = new;
                        }
                        AromaticSystemFieldChange::UnpairedElectrons { new, .. } => {
                            attributes.unpaired_electrons = new;
                        }
                    }
                }
            }
            Undo::ModifyMulticenterBondField { id, change } => {
                if id.index() < self.multicenter_bonds().count() {
                    let mut view = self.multicenter_bond_mut(id);
                    let attributes = view.attributes_mut();
                    match change {
                        MulticenterBondFieldChange::Electrons { new, .. } => {
                            attributes.electrons = new;
                        }
                        MulticenterBondFieldChange::Charge { new, .. } => {
                            attributes.charge = new;
                        }
                        MulticenterBondFieldChange::UnpairedElectrons { new, .. } => {
                            attributes.unpaired_electrons = new;
                        }
                    }
                }
            }
            Undo::ModifyNoncovalentBondField { id, change } => {
                if id.index() < self.noncovalent_bonds().count() {
                    let mut view = self.noncovalent_bond_mut(id);
                    let attributes = view.attributes_mut();
                    match change {
                        NoncovalentBondFieldChange::Kind { new, .. } => {
                            attributes.kind = new;
                        }
                    }
                }
            }
            Undo::ModifyStereoAtomField { id, change } => {
                if id.index() < self.stereo_atoms().count() {
                    let mut view = self.stereo_atom_mut(id);
                    let attributes = view.attributes_mut();
                    match change {
                        StereoAtomFieldChange::Configuration { new, .. } => {
                            attributes.configuration = new;
                        }
                    }
                }
            }
            Undo::ModifyStereoBondField { id, change } => {
                if id.index() < self.stereo_bonds().count() {
                    let mut view = self.stereo_bond_mut(id);
                    let attributes = view.attributes_mut();
                    match change {
                        StereoBondFieldChange::Configuration { new, .. } => {
                            attributes.configuration = new;
                        }
                    }
                }
            }
            Undo::RestoreAtomConstraints { key, constraints } => {
                for (id, constraint) in constraints.into_iter().rev() {
                    if id.index() < self.atoms().count() {
                        let mut view = self.atom_mut(id);
                        let constraints = &mut view.attributes_mut().constraints;
                        match constraint {
                            Some(value) => constraints.set(value),
                            None => {
                                constraints.remove(key);
                            }
                        }
                    }
                }
            }
            Undo::RestoreAtomConstraint {
                id,
                key,
                constraint,
            } => {
                if id.index() < self.atoms().count() {
                    let mut view = self.atom_mut(id);
                    let constraints = &mut view.attributes_mut().constraints;
                    match constraint {
                        Some(constraint) => constraints.set(constraint),
                        None => {
                            constraints.remove(key);
                        }
                    }
                }
            }
            Undo::RestoreBondConstraint {
                id,
                key,
                constraint,
            } => {
                if id.index() < self.bonds().count() {
                    let mut view = self.bond_mut(id);
                    let constraints = &mut view.attributes_mut().constraints;
                    match constraint {
                        Some(constraint) => constraints.set(constraint),
                        None => {
                            constraints.remove(key);
                        }
                    }
                }
            }
            Undo::RestoreDativeBondConstraint {
                id,
                key,
                constraint,
            } => {
                if id.index() < self.dative_bonds().count() {
                    let mut view = self.dative_bond_mut(id);
                    let constraints = &mut view.attributes_mut().constraints;
                    match constraint {
                        Some(constraint) => constraints.set(constraint),
                        None => {
                            constraints.remove(key);
                        }
                    }
                }
            }
            Undo::RestoreAromaticSystemConstraint {
                id,
                key,
                constraint,
            } => {
                if id.index() < self.aromatic_systems().count() {
                    let mut view = self.aromatic_system_mut(id);
                    let constraints = &mut view.attributes_mut().constraints;
                    match constraint {
                        Some(constraint) => constraints.set(constraint),
                        None => {
                            constraints.remove(key);
                        }
                    }
                }
            }
            Undo::RestoreMulticenterBondConstraint {
                id,
                key,
                constraint,
            } => {
                if id.index() < self.multicenter_bonds().count() {
                    let mut view = self.multicenter_bond_mut(id);
                    let constraints = &mut view.attributes_mut().constraints;
                    match constraint {
                        Some(constraint) => constraints.set(constraint),
                        None => {
                            constraints.remove(key);
                        }
                    }
                }
            }
            Undo::RestoreNoncovalentBondConstraint {
                id,
                key,
                constraint,
            } => {
                if id.index() < self.noncovalent_bonds().count() {
                    let mut view = self.noncovalent_bond_mut(id);
                    let constraints = &mut view.attributes_mut().constraints;
                    match constraint {
                        Some(constraint) => constraints.set(constraint),
                        None => {
                            constraints.remove(key);
                        }
                    }
                }
            }
            Undo::RestoreStereoAtomConstraint {
                id,
                key,
                constraint,
            } => {
                if id.index() < self.stereo_atoms().count() {
                    let mut view = self.stereo_atom_mut(id);
                    let constraints = &mut view.attributes_mut().constraints;
                    match constraint {
                        Some(constraint) => constraints.set(constraint),
                        None => {
                            constraints.remove(key);
                        }
                    }
                }
            }
            Undo::RestoreStereoBondConstraint {
                id,
                key,
                constraint,
            } => {
                if id.index() < self.stereo_bonds().count() {
                    let mut view = self.stereo_bond_mut(id);
                    let constraints = &mut view.attributes_mut().constraints;
                    match constraint {
                        Some(constraint) => constraints.set(constraint),
                        None => {
                            constraints.remove(key);
                        }
                    }
                }
            }
            Undo::RemoveAddedMoleculeConstraint { position } => {
                if position < self.constraints().len() {
                    self.remove_constraint_at(position);
                }
            }
            Undo::RestoreMoleculeConstraints(update) => {
                self.restore_constraints(&update);
            }
        }
    }

    fn capture_removed_topology(
        &self,
        atoms: &[AtomId],
        bonds: &[BondId],
    ) -> (Vec<RemovedAtom>, Vec<RemovedBond>, RemovedOverlays) {
        let atom_set: HashSet<AtomId> = atoms.iter().copied().collect();
        let bond_set: HashSet<BondId> = bonds.iter().copied().collect();
        let removed_atoms = atoms
            .iter()
            .map(|&id| RemovedAtom {
                id,
                attributes: self.atom(id).attributes().clone(),
            })
            .collect();
        let removed_bonds = (0..self.bonds().count())
            .map(BondId::from)
            .filter(|&id| {
                let view = self.bond(id);
                bond_set.contains(&id) || view.atom_ids().iter().any(|atom| atom_set.contains(atom))
            })
            .map(|id| {
                let view = self.bond(id);
                RemovedBond {
                    id,
                    endpoints: view.atom_ids(),
                    attributes: view.attributes().clone(),
                }
            })
            .collect();

        let dative_bonds = (0..self.dative_bonds().count())
            .map(DativeBondId::from)
            .filter_map(|id| {
                let view = self.dative_bond(id);
                let atoms: Vec<AtomId> = view.atom_ids().collect();
                atoms
                    .iter()
                    .any(|a| atom_set.contains(a))
                    .then(|| RemovedDativeBond {
                        id,
                        donors: view.donor_ids().collect(),
                        acceptor: view.acceptor_id(),
                        attributes: view.attributes().clone(),
                    })
            })
            .collect();
        let aromatic_systems = (0..self.aromatic_systems().count())
            .map(AromaticSystemId::from)
            .filter_map(|id| {
                let view = self.aromatic_system(id);
                let atoms: Vec<AtomId> = view.atom_ids().collect();
                atoms
                    .iter()
                    .any(|a| atom_set.contains(a))
                    .then(|| RemovedAromaticSystem {
                        id,
                        atoms,
                        attributes: view.attributes().clone(),
                    })
            })
            .collect();
        let multicenter_bonds = (0..self.multicenter_bonds().count())
            .map(MulticenterBondId::from)
            .filter_map(|id| {
                let view = self.multicenter_bond(id);
                let atoms: Vec<AtomId> = view.atom_ids().collect();
                atoms
                    .iter()
                    .any(|a| atom_set.contains(a))
                    .then(|| RemovedMulticenterBond {
                        id,
                        atoms,
                        attributes: view.attributes().clone(),
                    })
            })
            .collect();
        let noncovalent_bonds = (0..self.noncovalent_bonds().count())
            .map(NoncovalentBondId::from)
            .filter_map(|id| {
                let view = self.noncovalent_bond(id);
                view.atom_ids()
                    .iter()
                    .any(|a| atom_set.contains(a))
                    .then(|| RemovedNoncovalentBond {
                        id,
                        atoms: view.atom_ids(),
                        attributes: view.attributes().clone(),
                    })
            })
            .collect();

        // A stereo atom drops when its site atom or any ligand atom is removed;
        // a stereo bond drops when its site bond (directly or via a removed
        // endpoint) or any ligand atom is removed. Mirrors `birelation_removed`.
        let stereo_atoms = (0..self.stereo_atoms().count())
            .map(StereoAtomId::from)
            .filter_map(|id| {
                let view = self.stereo_atom(id);
                let dropped = atom_set.contains(&view.site_id())
                    || view
                        .ligand_ids()
                        .iter()
                        .any(|l| atom_set.contains(&l.atom_id));
                dropped.then(|| RemovedStereoAtom {
                    id,
                    site: view.site_id(),
                    ligands: view.ligand_ids().to_vec(),
                    attributes: view.attributes().clone(),
                })
            })
            .collect();
        let stereo_bonds = (0..self.stereo_bonds().count())
            .map(StereoBondId::from)
            .filter_map(|id| {
                let view = self.stereo_bond(id);
                let site = view.site_id();
                let site_dropped = bond_set.contains(&site)
                    || self
                        .bond(site)
                        .atom_ids()
                        .iter()
                        .any(|a| atom_set.contains(a));
                let ligand_dropped = view
                    .ligand_ids()
                    .iter()
                    .any(|l| atom_set.contains(&l.atom_id));
                (site_dropped || ligand_dropped).then(|| RemovedStereoBond {
                    id,
                    site,
                    ligands: view.ligand_ids().to_vec(),
                    attributes: view.attributes().clone(),
                })
            })
            .collect();

        (
            removed_atoms,
            removed_bonds,
            RemovedOverlays {
                dative_bonds,
                aromatic_systems,
                multicenter_bonds,
                noncovalent_bonds,
                stereo_atoms,
                stereo_bonds,
            },
        )
    }
}

fn ensure_unique<I>(ids: &[I], kind: EntityKind) -> Result<(), TransactionError>
where
    I: Copy + Eq + Hash,
{
    let mut seen = HashSet::with_capacity(ids.len());
    if ids.iter().copied().all(|id| seen.insert(id)) {
        Ok(())
    } else {
        Err(TransactionError::DuplicateRemoval { kind })
    }
}

#[cfg(test)]
mod tests;
