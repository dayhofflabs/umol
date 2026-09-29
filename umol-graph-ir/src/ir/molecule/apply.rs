//! Single-edit execution and undo replay on molecule storage.

use std::collections::HashSet;
use std::hash::Hash;

use umol_graph_core::{Compaction, GraphCompaction};
use umol_perm::{DynPermutation, Permutation};

use super::{Molecule, TransactionError};
use crate::ir::aromatic::AromaticSystemForm;
use crate::ir::compact::{MoleculeCompaction, UndoCompaction};
use crate::ir::constraint::{
    AromaticSystemConstraintForm, AtomConstraintForm, BondConstraintForm, Constraint,
    DativeBondConstraintForm, MulticenterBondConstraintForm, NoncovalentBondConstraintForm,
    StereoAtomConstraintForm, StereoBondConstraintForm,
};
use crate::ir::dative::DativeBondForm;
use crate::ir::edit::{
    AddBond, AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, AddedStereoAtom, AddedStereoBond, AromaticSystemFieldChange,
    AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange, BondHandle,
    CascadedConstraints, ConstraintEdit, DativeBondFieldChange, DativeBondHandle, Edit,
    MulticenterBondFieldChange, MulticenterBondHandle, NoncovalentBondFieldChange,
    NoncovalentBondHandle, RemovedAromaticSystem, RemovedAtom, RemovedBond, RemovedConstraint,
    RemovedDativeBond, RemovedMulticenterBond, RemovedNoncovalentBond, RemovedOverlays,
    RemovedStereoAtom, RemovedStereoBond, StereoAtomFieldChange, StereoAtomHandle,
    StereoBondFieldChange, StereoBondHandle, Undo,
};
use crate::ir::entity::EntityKind;
use crate::ir::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use crate::ir::ligand::{StereoLigand, StereoLigandKind};
use crate::ir::multicenter::MulticenterBondForm;
use crate::ir::noncovalent::NoncovalentBondForm;
use crate::ir::stereo::{StereoAtomForm, StereoBondForm};
use crate::ir::traits::{FrameTransport, Normalize};

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
                let dative_bonds = self.tracked_compact_dative_bonds(&graph);
                let aromatic_systems = self.tracked_compact_aromatic_systems(&graph);
                let multicenter_bonds = self.tracked_compact_multicenter_bonds(&graph);
                let noncovalent_bonds = self.tracked_compact_noncovalent_bonds(&graph);
                let stereo_atoms = self.tracked_compact_stereo_atoms(&graph);
                let stereo_bonds = self.tracked_compact_stereo_bonds(&graph);
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
                    if !self.dative_bond_equiv(id, acceptor, &donors, &attributes) {
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
            Edit::AddAromaticSystem { atoms, attributes } => {
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
                    if !self.aromatic_system_equiv(id, &saved_atoms, &attributes) {
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
            Edit::AddMulticenterBond { atoms, attributes } => {
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
                    if !self.multicenter_bond_equiv(id, &saved_atoms, &attributes) {
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
                    if !self.noncovalent_bond_equiv(id, saved_atoms, &attributes) {
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
            Edit::AddStereoAtom {
                site,
                ligands,
                attributes,
            } => {
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
                    if !self.stereo_atom_equiv(id, site, &ligands, &attributes) {
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
            Edit::AddStereoBond {
                site,
                ligands,
                attributes,
            } => {
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
                    if !self.stereo_bond_equiv(id, site, &ligands, &attributes) {
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
            Edit::ModifyAtomConstraint { id, old, new } => {
                let id = state.atom(id)?;
                self.apply_modify_atom_constraint(id, old, new)
            }
            Edit::ModifyBondConstraint { id, old, new } => {
                let id = state.bond(id)?;
                self.apply_modify_bond_constraint(id, old, new)
            }
            Edit::ModifyDativeBondConstraint { id, old, new } => {
                let id = state.dative_bond(id)?;
                self.apply_modify_dative_bond_constraint(id, old, new)
            }
            Edit::ModifyAromaticSystemConstraint { id, old, new } => {
                let id = state.aromatic_system(id)?;
                self.apply_modify_aromatic_system_constraint(id, old, new)
            }
            Edit::ModifyMulticenterBondConstraint { id, old, new } => {
                let id = state.multicenter_bond(id)?;
                self.apply_modify_multicenter_bond_constraint(id, old, new)
            }
            Edit::ModifyNoncovalentBondConstraint { id, old, new } => {
                let id = state.noncovalent_bond(id)?;
                self.apply_modify_noncovalent_bond_constraint(id, old, new)
            }
            Edit::ModifyStereoAtomConstraint {
                id,
                kind: _,
                old,
                new,
            } => {
                let id = state.stereo_atom(id)?;
                self.apply_modify_stereo_atom_constraint(id, old, new)
            }
            Edit::ModifyStereoBondConstraint {
                id,
                kind: _,
                old,
                new,
            } => {
                let id = state.stereo_bond(id)?;
                self.apply_modify_stereo_bond_constraint(id, old, new)
            }
            Edit::AddMoleculeConstraint { constraint } => {
                let constraint = state.resolve_constraint(constraint)?;
                self.push_constraint(constraint);
                Ok(())
            }
            Edit::RemoveMoleculeConstraint { constraint } => {
                let constraint = state.resolve_constraint(constraint)?;
                let list = self.constraints_mut();
                let position = list
                    .as_slice()
                    .iter()
                    .rposition(|c| *c == constraint)
                    .ok_or(TransactionError::MissingEntry)?;
                list.remove_at(position);
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
                let dative_bonds = self.tracked_compact_dative_bonds(&graph);
                let aromatic_systems = self.tracked_compact_aromatic_systems(&graph);
                let multicenter_bonds = self.tracked_compact_multicenter_bonds(&graph);
                let noncovalent_bonds = self.tracked_compact_noncovalent_bonds(&graph);
                let stereo_atoms = self.tracked_compact_stereo_atoms(&graph);
                let stereo_bonds = self.tracked_compact_stereo_bonds(&graph);
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
                Ok(Undo::RestoreRemovedTopology {
                    atoms: removed_atoms,
                    bonds: removed_bonds,
                    overlays,
                    compaction,
                    undo_compaction,
                    cascade,
                })
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
                    if !self.dative_bond_equiv(id, acceptor, &donors, &attributes) {
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
                Ok(Undo::RestoreRemovedDativeBonds {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
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
            Edit::AddAromaticSystem { atoms, attributes } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_aromatic_system(&resolved, attributes);
                state.push_aromatic_system(id);
                let view = self.aromatic_system(id);
                Ok(Undo::RemoveAddedAromaticSystem(AddedAromaticSystem {
                    id,
                    atoms: view.atom_ids().collect(),
                    attributes: view.attributes().clone(),
                }))
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
                    if !self.aromatic_system_equiv(id, &saved_atoms, &attributes) {
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
                Ok(Undo::RestoreRemovedAromaticSystems {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
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
            Edit::AddMulticenterBond { atoms, attributes } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_multicenter_bond(&resolved, attributes);
                state.push_multicenter_bond(id);
                let view = self.multicenter_bond(id);
                Ok(Undo::RemoveAddedMulticenterBond(AddedMulticenterBond {
                    id,
                    atoms: view.atom_ids().collect(),
                    attributes: view.attributes().clone(),
                }))
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
                    if !self.multicenter_bond_equiv(id, &saved_atoms, &attributes) {
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
                Ok(Undo::RestoreRemovedMulticenterBonds {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
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
                    if !self.noncovalent_bond_equiv(id, saved_atoms, &attributes) {
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
                Ok(Undo::RestoreRemovedNoncovalentBonds {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
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
            Edit::AddStereoAtom {
                site,
                ligands,
                attributes,
            } => {
                let site = state.atom(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_atom(site, &ligands, attributes);
                state.push_stereo_atom(id);
                let view = self.stereo_atom(id);
                Ok(Undo::RemoveAddedStereoAtom(AddedStereoAtom {
                    id,
                    site,
                    ligands,
                    attributes: view.attributes().clone(),
                }))
            }
            Edit::RemoveStereoAtoms { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, site, ligands, attributes) in removes {
                    let id = state.stereo_atom(id)?;
                    let site = state.atom(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    if !self.stereo_atom_equiv(id, site, &ligands, &attributes) {
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
                Ok(Undo::RestoreRemovedStereoAtoms {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
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
            Edit::AddStereoBond {
                site,
                ligands,
                attributes,
            } => {
                let site = state.bond(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_bond(site, &ligands, attributes);
                state.push_stereo_bond(id);
                let view = self.stereo_bond(id);
                Ok(Undo::RemoveAddedStereoBond(AddedStereoBond {
                    id,
                    site,
                    ligands,
                    attributes: view.attributes().clone(),
                }))
            }
            Edit::RemoveStereoBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, site, ligands, attributes) in removes {
                    let id = state.stereo_bond(id)?;
                    let site = state.bond(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    if !self.stereo_bond_equiv(id, site, &ligands, &attributes) {
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
                Ok(Undo::RestoreRemovedStereoBonds {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
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
            Edit::ModifyAtomConstraint { id, old, new } => {
                let id = state.atom(id.clone())?;
                if old.is_none() && new.is_none() {
                    return Ok(None);
                }
                let undo = Undo::ApplyEdit(Box::new(Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(id),
                    old: new.clone(),
                    new: old.clone(),
                }));
                self.apply_modify_atom_constraint(id, old, new)?;
                Ok(undo)
            }
            Edit::ModifyBondConstraint { id, old, new } => {
                let id = state.bond(id.clone())?;
                if old.is_none() && new.is_none() {
                    return Ok(None);
                }
                let undo = Undo::ApplyEdit(Box::new(Edit::ModifyBondConstraint {
                    id: BondHandle::Id(id),
                    old: new.clone(),
                    new: old.clone(),
                }));
                self.apply_modify_bond_constraint(id, old, new)?;
                Ok(undo)
            }
            Edit::ModifyDativeBondConstraint { id, old, new } => {
                let id = state.dative_bond(id.clone())?;
                if old.is_none() && new.is_none() {
                    return Ok(None);
                }
                let undo = Undo::ApplyEdit(Box::new(Edit::ModifyDativeBondConstraint {
                    id: DativeBondHandle::Id(id),
                    old: new.clone(),
                    new: old.clone(),
                }));
                self.apply_modify_dative_bond_constraint(id, old, new)?;
                Ok(undo)
            }
            Edit::ModifyAromaticSystemConstraint { id, old, new } => {
                let id = state.aromatic_system(id.clone())?;
                if old.is_none() && new.is_none() {
                    return Ok(None);
                }
                let undo = Undo::ApplyEdit(Box::new(Edit::ModifyAromaticSystemConstraint {
                    id: AromaticSystemHandle::Id(id),
                    old: new.clone(),
                    new: old.clone(),
                }));
                self.apply_modify_aromatic_system_constraint(id, old, new)?;
                Ok(undo)
            }
            Edit::ModifyMulticenterBondConstraint { id, old, new } => {
                let id = state.multicenter_bond(id.clone())?;
                if old.is_none() && new.is_none() {
                    return Ok(None);
                }
                let undo = Undo::ApplyEdit(Box::new(Edit::ModifyMulticenterBondConstraint {
                    id: MulticenterBondHandle::Id(id),
                    old: new.clone(),
                    new: old.clone(),
                }));
                self.apply_modify_multicenter_bond_constraint(id, old, new)?;
                Ok(undo)
            }
            Edit::ModifyNoncovalentBondConstraint { id, old, new } => {
                let id = state.noncovalent_bond(id.clone())?;
                if old.is_none() && new.is_none() {
                    return Ok(None);
                }
                let undo = Undo::ApplyEdit(Box::new(Edit::ModifyNoncovalentBondConstraint {
                    id: NoncovalentBondHandle::Id(id),
                    old: new.clone(),
                    new: old.clone(),
                }));
                self.apply_modify_noncovalent_bond_constraint(id, old, new)?;
                Ok(undo)
            }
            Edit::ModifyStereoAtomConstraint { id, kind, old, new } => {
                let id = state.stereo_atom(id.clone())?;
                if old.is_none() && new.is_none() {
                    return Ok(None);
                }
                let undo = Undo::ApplyEdit(Box::new(Edit::ModifyStereoAtomConstraint {
                    id: StereoAtomHandle::Id(id),
                    kind,
                    old: new.clone(),
                    new: old.clone(),
                }));
                self.apply_modify_stereo_atom_constraint(id, old, new)?;
                Ok(undo)
            }
            Edit::ModifyStereoBondConstraint { id, kind, old, new } => {
                let id = state.stereo_bond(id.clone())?;
                if old.is_none() && new.is_none() {
                    return Ok(None);
                }
                let undo = Undo::ApplyEdit(Box::new(Edit::ModifyStereoBondConstraint {
                    id: StereoBondHandle::Id(id),
                    kind,
                    old: new.clone(),
                    new: old.clone(),
                }));
                self.apply_modify_stereo_bond_constraint(id, old, new)?;
                Ok(undo)
            }
            Edit::AddMoleculeConstraint { constraint } => {
                let constraint = state.resolve_constraint(constraint)?;
                self.push_constraint(constraint.clone());
                Ok(Undo::ApplyEdit(Box::new(Edit::RemoveMoleculeConstraint {
                    constraint: constraint.into(),
                })))
            }
            Edit::RemoveMoleculeConstraint { constraint } => {
                let constraint = state.resolve_constraint(constraint)?;
                let list = self.constraints_mut();
                let position = list
                    .as_slice()
                    .iter()
                    .rposition(|c| *c == constraint)
                    .ok_or(TransactionError::MissingEntry)?;
                let constraint = list.remove_at(position);
                Ok(Undo::ApplyCascadedConstraints(CascadedConstraints {
                    removed: vec![RemovedConstraint {
                        position,
                        constraint,
                    }],
                    modified: Vec::new(),
                }))
            }
        }
        .map(Some)
    }

    pub(crate) fn apply_undo(&mut self, undo: Undo) {
        self.validate_undo(&undo)?;
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
            Undo::RestoreRemovedTopology {
                atoms,
                bonds,
                overlays,
                undo_compaction,
                cascade,
                ..
            } => {
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
            Undo::RestoreRemovedDativeBonds {
                removed,
                undo_compaction,
                cascade,
            } => {
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
                self.dative_bond_view_mut(id).replace_donors(&donors);
            }
            Undo::RestoreDativeBondAcceptor { id, acceptor } => {
                self.dative_bond_view_mut(id).replace_acceptor(acceptor);
            }
            Undo::RemoveAddedAromaticSystem(added) => {
                if added.id.index() < self.aromatic_systems().count() {
                    self.remove_aromatic_systems(&[added.id]);
                }
            }
            Undo::RestoreRemovedAromaticSystems {
                removed,
                undo_compaction,
                cascade,
            } => {
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
                self.aromatic_system_view_mut(id).replace_atoms(&atoms);
            }
            Undo::RemoveAddedMulticenterBond(added) => {
                if added.id.index() < self.multicenter_bonds().count() {
                    self.remove_multicenter_bonds(&[added.id]);
                }
            }
            Undo::RestoreRemovedMulticenterBonds {
                removed,
                undo_compaction,
                cascade,
            } => {
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
                self.multicenter_bond_view_mut(id).replace_atoms(&atoms);
            }
            Undo::RemoveAddedNoncovalentBond(added) => {
                if added.id.index() < self.noncovalent_bonds().count() {
                    self.remove_noncovalent_bonds(&[added.id]);
                }
            }
            Undo::RestoreRemovedNoncovalentBonds {
                removed,
                undo_compaction,
                cascade,
            } => {
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
                self.noncovalent_bond_view_mut(id).replace_atoms(atoms);
            }
            Undo::RemoveAddedStereoAtom(added) => {
                if added.id.index() < self.stereo_atoms().count() {
                    self.remove_stereo_atoms(&[added.id]);
                }
            }
            Undo::RestoreRemovedStereoAtoms {
                removed,
                undo_compaction,
                cascade,
            } => {
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
                self.stereo_atom_view_mut(id).replace_site(site);
            }
            Undo::RestoreStereoAtomLigands { id, ligands } => {
                self.stereo_atom_view_mut(id).replace_ligands(&ligands);
            }
            Undo::RemoveAddedStereoBond(added) => {
                if added.id.index() < self.stereo_bonds().count() {
                    self.remove_stereo_bonds(&[added.id]);
                }
            }
            Undo::RestoreRemovedStereoBonds {
                removed,
                undo_compaction,
                cascade,
            } => {
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
                self.stereo_bond_view_mut(id).replace_site(site);
            }
            Undo::RestoreStereoBondLigands { id, ligands } => {
                self.stereo_bond_view_mut(id).replace_ligands(&ligands);
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
            Undo::ApplyCascadedConstraints(update) => {
                self.restore_constraints(&update);
            }
            Undo::ApplyEdit(edit) => {
                let mut state = ApplicationState::new(self);
                self.apply_edit(*edit, &mut state)
                    .map_err(|_| rollback_mismatch())?;
            }
        }
    }

    /// Compare an acceptor, unordered donors, and attributes in the stored donor frame.
    fn dative_bond_equiv(
        &self,
        id: DativeBondId,
        acceptor: AtomId,
        donors: &[AtomId],
        attributes: &DativeBondForm,
    ) -> bool {
        let set = self.raw_dative_bonds();
        let stored: Vec<AtomId> = set.donors(id).collect();
        set.is_coincident(id, acceptor, donors)
            && DynPermutation::between(donors, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// `true` iff aromatic system `id` structurally equals `(atoms, attributes)`.
    fn aromatic_system_equiv(
        &self,
        id: AromaticSystemId,
        atoms: &[AtomId],
        attributes: &AromaticSystemForm,
    ) -> bool {
        let set = self.raw_aromatic_systems();
        let stored: Vec<AtomId> = set.atoms(id).collect();
        set.is_coincident(id, atoms)
            && DynPermutation::between(atoms, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// `true` iff multicenter bond `id` structurally equals `(atoms, attributes)`.
    fn multicenter_bond_equiv(
        &self,
        id: MulticenterBondId,
        atoms: &[AtomId],
        attributes: &MulticenterBondForm,
    ) -> bool {
        let set = self.raw_multicenter_bonds();
        let stored: Vec<AtomId> = set.atoms(id).collect();
        set.is_coincident(id, atoms)
            && DynPermutation::between(atoms, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// Compare an unordered atom pair and its attributes in the stored frame.
    fn noncovalent_bond_equiv(
        &self,
        id: NoncovalentBondId,
        atoms: [AtomId; 2],
        attributes: &NoncovalentBondForm,
    ) -> bool {
        let set = self.raw_noncovalent_bonds();
        let stored = set.atoms(id);
        set.is_coincident(id, atoms[0], atoms[1])
            && DynPermutation::between(&atoms, &stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// `true` iff stereo atom `id` structurally equals `(site, ligands, attributes)`.
    fn stereo_atom_equiv(
        &self,
        id: StereoAtomId,
        site: AtomId,
        ligands: &[StereoLigand],
        attributes: &StereoAtomForm,
    ) -> bool {
        let set = self.raw_stereo_atoms();
        let stored = set.ligands(id);
        set.site(id) == site
            && Permutation::between(ligands, stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
    }

    /// `true` iff stereo bond `id` structurally equals `(site, ligands, attributes)`.
    fn stereo_bond_equiv(
        &self,
        id: StereoBondId,
        site: BondId,
        ligands: &[StereoLigand],
        attributes: &StereoBondForm,
    ) -> bool {
        let set = self.raw_stereo_bonds();
        let stored = set.ligands(id);
        set.site(id) == site
            && Permutation::between(ligands, stored)
                .and_then(|action| attributes.clone().reframe_by(&action))
                .is_some_and(|restated| restated.normalized_eq(set.attributes(id)))
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

    fn apply_modify_atom_constraint(
        &mut self,
        id: AtomId,
        old: Option<AtomConstraintForm>,
        new: Option<AtomConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.atom_mut(id)
            .attributes_mut()
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_bond_constraint(
        &mut self,
        id: BondId,
        old: Option<BondConstraintForm>,
        new: Option<BondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.bond_mut(id)
            .attributes_mut()
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_dative_bond_constraint(
        &mut self,
        id: DativeBondId,
        old: Option<DativeBondConstraintForm>,
        new: Option<DativeBondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.dative_bond_mut(id)
            .attributes_mut()
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_aromatic_system_constraint(
        &mut self,
        id: AromaticSystemId,
        old: Option<AromaticSystemConstraintForm>,
        new: Option<AromaticSystemConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.aromatic_system_mut(id)
            .attributes_mut()
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_multicenter_bond_constraint(
        &mut self,
        id: MulticenterBondId,
        old: Option<MulticenterBondConstraintForm>,
        new: Option<MulticenterBondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.multicenter_bond_mut(id)
            .attributes_mut()
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_noncovalent_bond_constraint(
        &mut self,
        id: NoncovalentBondId,
        old: Option<NoncovalentBondConstraintForm>,
        new: Option<NoncovalentBondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.noncovalent_bond_mut(id)
            .attributes_mut()
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_stereo_atom_constraint(
        &mut self,
        id: StereoAtomId,
        old: Option<StereoAtomConstraintForm>,
        new: Option<StereoAtomConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.stereo_atom_mut(id)
            .attributes_mut()
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_stereo_bond_constraint(
        &mut self,
        id: StereoBondId,
        old: Option<StereoBondConstraintForm>,
        new: Option<StereoBondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.stereo_bond_mut(id)
            .attributes_mut()
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }
    fn validate_undo(&self, undo: &Undo) -> Result<(), TransactionError> {
        let compaction = match undo {
            Undo::RestoreRemovedTopology {
                undo_compaction, ..
            }
            | Undo::RestoreRemovedDativeBonds {
                undo_compaction, ..
            }
            | Undo::RestoreRemovedAromaticSystems {
                undo_compaction, ..
            }
            | Undo::RestoreRemovedMulticenterBonds {
                undo_compaction, ..
            }
            | Undo::RestoreRemovedNoncovalentBonds {
                undo_compaction, ..
            }
            | Undo::RestoreRemovedStereoAtoms {
                undo_compaction, ..
            }
            | Undo::RestoreRemovedStereoBonds {
                undo_compaction, ..
            } => Some(undo_compaction.forward()),
            _ => None,
        };
        if let Some(compaction) = compaction {
            let counts_match = [
                (
                    self.atoms().count(),
                    compaction.graph().nodes().result_count(),
                ),
                (
                    self.bonds().count(),
                    compaction.graph().edges().result_count(),
                ),
                (
                    self.dative_bonds().count(),
                    compaction.dative_bonds().result_count(),
                ),
                (
                    self.aromatic_systems().count(),
                    compaction.aromatic_systems().result_count(),
                ),
                (
                    self.multicenter_bonds().count(),
                    compaction.multicenter_bonds().result_count(),
                ),
                (
                    self.noncovalent_bonds().count(),
                    compaction.noncovalent_bonds().result_count(),
                ),
                (
                    self.stereo_atoms().count(),
                    compaction.stereo_atoms().result_count(),
                ),
                (
                    self.stereo_bonds().count(),
                    compaction.stereo_bonds().result_count(),
                ),
            ]
            .into_iter()
            .all(|(actual, expected)| actual == expected);
            if !counts_match {
                return Err(rollback_mismatch());
            }
        }
        let fits = match undo {
            Undo::RemoveAddedTopology { atoms, bonds } => {
                ids_fit(
                    atoms.iter().map(|entry| entry.id.index()),
                    self.atoms().count(),
                ) && ids_fit(
                    bonds.iter().map(|entry| entry.id.index()),
                    self.bonds().count(),
                )
            }
            Undo::RestoreRemovedTopology {
                atoms,
                bonds,
                overlays,
                compaction,
                undo_compaction,
                cascade: _,
            } => {
                let atom_count = self.atoms().count() + atoms.len();
                let bond_count = self.bonds().count() + bonds.len();
                undo_compaction.forward() == compaction
                    && reconstruction_fits(
                        self.atoms().count(),
                        atoms.iter().map(|entry| entry.id.index()),
                        |id| undo_compaction.uncompact_atom(AtomId::from(id)).index(),
                    )
                    && reconstruction_fits(
                        self.bonds().count(),
                        bonds.iter().map(|entry| entry.id.index()),
                        |id| undo_compaction.uncompact_bond(BondId::from(id)).index(),
                    )
                    && bonds
                        .iter()
                        .all(|entry| entry.endpoints.iter().all(|id| id.index() < atom_count))
                    && removed_overlays_fit(self, overlays, undo_compaction, atom_count, bond_count)
            }
            Undo::RemoveAddedDativeBond(entry) => entry.id.index() < self.dative_bonds().count(),
            Undo::RestoreRemovedDativeBonds {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.dative_bonds().count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_dative_bond(DativeBondId::from(id))
                            .index()
                    },
                ) && removed_dative_bonds_fit(removed, self.atoms().count())
            }
            Undo::RestoreDativeBondDonors { id, .. }
            | Undo::RestoreDativeBondAcceptor { id, .. } => {
                id.index() < self.dative_bonds().count()
            }
            Undo::RemoveAddedAromaticSystem(entry) => {
                entry.id.index() < self.aromatic_systems().count()
            }
            Undo::RestoreRemovedAromaticSystems {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.aromatic_systems().count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_aromatic_system(AromaticSystemId::from(id))
                            .index()
                    },
                ) && removed_aromatic_systems_fit(removed, self.atoms().count())
            }
            Undo::RestoreAromaticSystemAtoms { id, .. } => {
                id.index() < self.aromatic_systems().count()
            }
            Undo::RemoveAddedMulticenterBond(entry) => {
                entry.id.index() < self.multicenter_bonds().count()
            }
            Undo::RestoreRemovedMulticenterBonds {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.multicenter_bonds().count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_multicenter_bond(MulticenterBondId::from(id))
                            .index()
                    },
                ) && removed_multicenter_bonds_fit(removed, self.atoms().count())
            }
            Undo::RestoreMulticenterBondAtoms { id, .. } => {
                id.index() < self.multicenter_bonds().count()
            }
            Undo::RemoveAddedNoncovalentBond(entry) => {
                entry.id.index() < self.noncovalent_bonds().count()
            }
            Undo::RestoreRemovedNoncovalentBonds {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.noncovalent_bonds().count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_noncovalent_bond(NoncovalentBondId::from(id))
                            .index()
                    },
                ) && removed_noncovalent_bonds_fit(removed, self.atoms().count())
            }
            Undo::RestoreNoncovalentBondAtoms { id, .. } => {
                id.index() < self.noncovalent_bonds().count()
            }
            Undo::RemoveAddedStereoAtom(entry) => entry.id.index() < self.stereo_atoms().count(),
            Undo::RestoreRemovedStereoAtoms {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.stereo_atoms().count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_stereo_atom(StereoAtomId::from(id))
                            .index()
                    },
                ) && removed_stereo_atoms_fit(removed, self.atoms().count())
            }
            Undo::RestoreStereoAtomSite { id, .. } | Undo::RestoreStereoAtomLigands { id, .. } => {
                id.index() < self.stereo_atoms().count()
            }
            Undo::RemoveAddedStereoBond(entry) => entry.id.index() < self.stereo_bonds().count(),
            Undo::RestoreRemovedStereoBonds {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.stereo_bonds().count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_stereo_bond(StereoBondId::from(id))
                            .index()
                    },
                ) && removed_stereo_bonds_fit(removed, self.atoms().count(), self.bonds().count())
            }
            Undo::RestoreStereoBondSite { id, .. } | Undo::RestoreStereoBondLigands { id, .. } => {
                id.index() < self.stereo_bonds().count()
            }
            Undo::ModifyAtomField { id, .. } => id.index() < self.atoms().count(),
            Undo::ModifyBondField { id, .. } => id.index() < self.bonds().count(),
            Undo::ModifyDativeBondField { id, .. } => id.index() < self.dative_bonds().count(),
            Undo::ModifyAromaticSystemField { id, .. } => {
                id.index() < self.aromatic_systems().count()
            }
            Undo::ModifyMulticenterBondField { id, .. } => {
                id.index() < self.multicenter_bonds().count()
            }
            Undo::ModifyNoncovalentBondField { id, .. } => {
                id.index() < self.noncovalent_bonds().count()
            }
            Undo::ModifyStereoAtomField { id, .. } => id.index() < self.stereo_atoms().count(),
            Undo::ModifyStereoBondField { id, .. } => id.index() < self.stereo_bonds().count(),
            Undo::ApplyCascadedConstraints(_) => true,
            Undo::ApplyEdit(_) => true,
        };
        fits.then_some(()).ok_or_else(rollback_mismatch)
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

fn ids_fit(ids: impl IntoIterator<Item = usize>, count: usize) -> bool {
    let mut seen = HashSet::new();
    ids.into_iter().all(|id| id < count && seen.insert(id))
}

fn reconstruction_fits(
    current_count: usize,
    removed: impl IntoIterator<Item = usize>,
    mut uncompact: impl FnMut(usize) -> usize,
) -> bool {
    let removed: Vec<_> = removed.into_iter().collect();
    let restored_count = current_count + removed.len();
    let mut occupied = vec![false; restored_count];
    for id in removed {
        if id >= restored_count || occupied[id] {
            return false;
        }
        occupied[id] = true;
    }
    for id in 0..current_count {
        let restored = uncompact(id);
        if restored >= restored_count || occupied[restored] {
            return false;
        }
        occupied[restored] = true;
    }
    occupied.into_iter().all(|is_occupied| is_occupied)
}

fn removed_dative_bonds_fit(removed: &[RemovedDativeBond], atom_count: usize) -> bool {
    removed.iter().all(|entry| {
        entry.acceptor.index() < atom_count && entry.donors.iter().all(|id| id.index() < atom_count)
    })
}

fn removed_aromatic_systems_fit(removed: &[RemovedAromaticSystem], atom_count: usize) -> bool {
    removed
        .iter()
        .all(|entry| entry.atoms.iter().all(|id| id.index() < atom_count))
}

fn removed_multicenter_bonds_fit(removed: &[RemovedMulticenterBond], atom_count: usize) -> bool {
    removed
        .iter()
        .all(|entry| entry.atoms.iter().all(|id| id.index() < atom_count))
}

fn removed_noncovalent_bonds_fit(removed: &[RemovedNoncovalentBond], atom_count: usize) -> bool {
    removed
        .iter()
        .all(|entry| entry.atoms.iter().all(|id| id.index() < atom_count))
}

fn removed_stereo_atoms_fit(removed: &[RemovedStereoAtom], atom_count: usize) -> bool {
    removed.iter().all(|entry| {
        entry.site.index() < atom_count
            && entry
                .ligands
                .iter()
                .all(|ligand| ligand.atom_id.index() < atom_count)
    })
}

fn removed_stereo_bonds_fit(
    removed: &[RemovedStereoBond],
    atom_count: usize,
    bond_count: usize,
) -> bool {
    removed.iter().all(|entry| {
        entry.site.index() < bond_count
            && entry
                .ligands
                .iter()
                .all(|ligand| ligand.atom_id.index() < atom_count)
    })
}

fn removed_overlays_fit(
    molecule: &Molecule,
    removed: &RemovedOverlays,
    undo_compaction: &UndoCompaction,
    atom_count: usize,
    bond_count: usize,
) -> bool {
    reconstruction_fits(
        molecule.dative_bonds().count(),
        removed.dative_bonds.iter().map(|entry| entry.id.index()),
        |id| {
            undo_compaction
                .uncompact_dative_bond(DativeBondId::from(id))
                .index()
        },
    ) && reconstruction_fits(
        molecule.aromatic_systems().count(),
        removed
            .aromatic_systems
            .iter()
            .map(|entry| entry.id.index()),
        |id| {
            undo_compaction
                .uncompact_aromatic_system(AromaticSystemId::from(id))
                .index()
        },
    ) && reconstruction_fits(
        molecule.multicenter_bonds().count(),
        removed
            .multicenter_bonds
            .iter()
            .map(|entry| entry.id.index()),
        |id| {
            undo_compaction
                .uncompact_multicenter_bond(MulticenterBondId::from(id))
                .index()
        },
    ) && reconstruction_fits(
        molecule.noncovalent_bonds().count(),
        removed
            .noncovalent_bonds
            .iter()
            .map(|entry| entry.id.index()),
        |id| {
            undo_compaction
                .uncompact_noncovalent_bond(NoncovalentBondId::from(id))
                .index()
        },
    ) && reconstruction_fits(
        molecule.stereo_atoms().count(),
        removed.stereo_atoms.iter().map(|entry| entry.id.index()),
        |id| {
            undo_compaction
                .uncompact_stereo_atom(StereoAtomId::from(id))
                .index()
        },
    ) && reconstruction_fits(
        molecule.stereo_bonds().count(),
        removed.stereo_bonds.iter().map(|entry| entry.id.index()),
        |id| {
            undo_compaction
                .uncompact_stereo_bond(StereoBondId::from(id))
                .index()
        },
    ) && removed_dative_bonds_fit(&removed.dative_bonds, atom_count)
        && removed_aromatic_systems_fit(&removed.aromatic_systems, atom_count)
        && removed_multicenter_bonds_fit(&removed.multicenter_bonds, atom_count)
        && removed_noncovalent_bonds_fit(&removed.noncovalent_bonds, atom_count)
        && removed_stereo_atoms_fit(&removed.stereo_atoms, atom_count)
        && removed_stereo_bonds_fit(&removed.stereo_bonds, atom_count, bond_count)
}

fn rollback_mismatch() -> TransactionError {
    TransactionError::RollbackStateMismatch
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::ir::atom::{AtomForm, ElementForm, IsotopeMassForm};
    use crate::ir::bond::BondForm;
    use crate::ir::constraint::{RelationalConstraint, StereogenicityForm};
    use crate::ir::edit::ModifiedConstraint;
    use crate::ir::electrons::ElectronCountsForm;
    use crate::ir::molecule::MoleculeEntries;
    use crate::ir::noncovalent::{NoncovalentBondKind, NoncovalentBondKindForm};
    use crate::ir::num::NumForm;
    use crate::ir::spin::UnpairedElectronsForm;
    use crate::ir::stereo::{StereoConfigurationForm, StereoCoset, StereoKind};

    #[rstest]
    #[case::first(3, 0, Ok(AtomId(0)))]
    #[case::last(3, 2, Ok(AtomId(2)))]
    #[case::outside(3, 3, Err(TransactionError::HandleOutOfRange {
        kind: EntityKind::Atom, index: 3, count: 3,
    }))]
    #[case::empty(0, 0, Err(TransactionError::HandleOutOfRange {
        kind: EntityKind::Atom, index: 0, count: 0,
    }))]
    fn test_handle_table_initial(
        #[case] count: usize,
        #[case] index: usize,
        #[case] expected: Result<AtomId, TransactionError>,
    ) {
        let table = HandleTable::new(count);

        assert_eq!(table.initial(EntityKind::Atom, index), expected);
        assert_eq!(table.initial, None);
    }

    #[rstest]
    #[case::first(0, Ok(AtomId(3)))]
    #[case::second(1, Ok(AtomId(4)))]
    #[case::outside(2, Err(TransactionError::HandleOutOfRange {
        kind: EntityKind::Atom, index: 2, count: 2,
    }))]
    fn test_handle_table_created(
        #[case] index: usize,
        #[case] expected: Result<AtomId, TransactionError>,
    ) {
        let mut table = HandleTable::new(3);
        table.push(AtomId(3));
        table.push(AtomId(4));

        assert_eq!(table.created(EntityKind::Atom, index), expected);
        assert_eq!(table.initial, None);
    }

    #[rstest]
    #[case::identity(vec![], vec![Some(AtomId(0)), Some(AtomId(1)), Some(AtomId(2))], vec![Some(AtomId(3)), Some(AtomId(4))])]
    #[case::initial_and_created(vec![AtomId(1), AtomId(3)], vec![Some(AtomId(0)), None, Some(AtomId(1))], vec![None, Some(AtomId(2))])]
    #[case::all(vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)], vec![None, None, None], vec![None, None])]
    fn test_handle_table_compact(
        #[case] removed: Vec<AtomId>,
        #[case] initial: Vec<Option<AtomId>>,
        #[case] created: Vec<Option<AtomId>>,
    ) {
        let mut table = HandleTable::new(3);
        table.push(AtomId(3));
        table.push(AtomId(4));
        let compaction = Compaction::new(5, removed).unwrap();

        table.compact(|id| compaction.compact(id));

        assert_eq!(table.initial_count, 3);
        assert_eq!(table.initial, Some(initial));
        assert_eq!(table.created, created);
    }

    #[rstest]
    fn test_handle_table_compact_sequence() {
        let mut table = HandleTable::new(3);
        table.push(AtomId(3));
        table.push(AtomId(4));
        let first = Compaction::new(5, vec![AtomId(1), AtomId(3)]).unwrap();
        table.compact(|id| first.compact(id));
        table.push(AtomId(3));
        let second = Compaction::new(4, vec![AtomId(0), AtomId(2)]).unwrap();
        table.compact(|id| second.compact(id));

        assert_eq!(table.initial_count, 3);
        assert_eq!(table.initial, Some(vec![None, None, Some(AtomId(0))]));
        assert_eq!(table.created, vec![None, None, Some(AtomId(1))]);
        assert_eq!(table.initial(EntityKind::Atom, 2), Ok(AtomId(0)));
        assert_eq!(table.created(EntityKind::Atom, 2), Ok(AtomId(1)));
        assert_eq!(
            table.initial(EntityKind::Atom, 1),
            Err(TransactionError::HandleRemoved {
                kind: EntityKind::Atom,
                index: 1
            }),
        );
        assert_eq!(
            table.created(EntityKind::Atom, 1),
            Err(TransactionError::HandleRemoved {
                kind: EntityKind::Atom,
                index: 1
            }),
        );
        assert_eq!(
            table.initial(EntityKind::Atom, 3),
            Err(TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 3,
                count: 3
            }),
        );
        assert_eq!(
            table.created(EntityKind::Atom, 3),
            Err(TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 3,
                count: 3
            }),
        );
    }

    #[rstest]
    fn test_molecule_apply_edit_additions(#[values(false, true)] journaled: bool) {
        let initial = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
            ],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        });
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let edits = [
            Edit::AddAtoms {
                atoms: vec![
                    AtomForm::from_element(Element::O),
                    AtomForm::from_element(Element::F),
                ],
            },
            Edit::AddBonds {
                bonds: vec![
                    AddBond {
                        endpoints: [AtomHandle::Id(AtomId(1)), AtomHandle::New(0)],
                        attributes: BondForm::from_order(2),
                    },
                    AddBond {
                        endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::New(1)],
                        attributes: BondForm::from_order(1),
                    },
                    AddBond {
                        endpoints: [AtomHandle::New(0), AtomHandle::New(1)],
                        attributes: BondForm::from_order(1),
                    },
                ],
            },
            Edit::AddDativeBond {
                donors: vec![AtomHandle::New(1), AtomHandle::Id(AtomId(0))],
                acceptor: AtomHandle::New(0),
                attributes: DativeBondForm::from_order(1),
            },
            Edit::AddAromaticSystem {
                atoms: vec![
                    AtomHandle::New(0),
                    AtomHandle::Id(AtomId(1)),
                    AtomHandle::Id(AtomId(0)),
                ],
                attributes: AromaticSystemForm::from_electrons(vec![2, 1, 1]),
            },
            Edit::AddMulticenterBond {
                atoms: vec![AtomHandle::New(1), AtomHandle::New(0)],
                attributes: MulticenterBondForm::from_electrons(vec![1, 2]),
            },
            Edit::AddNoncovalentBond {
                atoms: [AtomHandle::New(1), AtomHandle::Id(AtomId(0))],
                attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            },
            Edit::AddStereoAtom {
                site: AtomHandle::Id(AtomId(0)),
                ligands: vec![
                    (AtomHandle::New(1), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(0)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                ],
                attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
            },
            Edit::AddStereoBond {
                site: BondHandle::New(0),
                ligands: vec![
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                    (
                        AtomHandle::Id(AtomId(1)),
                        StereoLigandKind::ImplicitHydrogen,
                    ),
                    (AtomHandle::New(1), StereoLigandKind::Atom),
                    (AtomHandle::New(0), StereoLigandKind::ImplicitHydrogen),
                ],
                attributes: StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            },
        ];
        let mut undos = Vec::new();
        for edit in edits {
            if journaled {
                undos.push(
                    molecule
                        .apply_edit_with_undo(edit, &mut state)
                        .unwrap()
                        .unwrap(),
                );
            } else {
                molecule.apply_edit(edit, &mut state).unwrap();
            }
        }

        let stereo_atom_ligands = vec![
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
        ];
        let stereo_bond_ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
        ];
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
                AtomForm::from_element(Element::F),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![(
                vec![AtomId(3), AtomId(0)],
                AtomId(2),
                DativeBondForm::from_order(1),
            )],
            aromatic: vec![(
                vec![AtomId(2), AtomId(1), AtomId(0)],
                AromaticSystemForm::from_electrons(vec![2, 1, 1]),
            )],
            multicenter: vec![(
                vec![AtomId(3), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            )],
            noncovalent: vec![(
                [AtomId(3), AtomId(0)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            stereo_atoms: vec![(
                AtomId(0),
                stereo_atom_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
            )],
            stereo_bonds: vec![(
                BondId(1),
                stereo_bond_ligands.clone(),
                StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            )],
            ..Default::default()
        });
        assert_eq!(molecule, expected);
        assert_eq!(state.atom(AtomHandle::Id(AtomId(1))), Ok(AtomId(1)));
        assert_eq!(state.atom(AtomHandle::New(0)), Ok(AtomId(2)));
        assert_eq!(state.atom(AtomHandle::New(1)), Ok(AtomId(3)));
        assert_eq!(state.bond(BondHandle::Id(BondId(0))), Ok(BondId(0)));
        assert_eq!(state.bond(BondHandle::New(0)), Ok(BondId(1)));
        assert_eq!(state.bond(BondHandle::New(2)), Ok(BondId(3)));
        assert_eq!(
            state.dative_bond(DativeBondHandle::New(0)),
            Ok(DativeBondId(0))
        );
        assert_eq!(
            state.aromatic_system(AromaticSystemHandle::New(0)),
            Ok(AromaticSystemId(0))
        );
        assert_eq!(
            state.multicenter_bond(MulticenterBondHandle::New(0)),
            Ok(MulticenterBondId(0))
        );
        assert_eq!(
            state.noncovalent_bond(NoncovalentBondHandle::New(0)),
            Ok(NoncovalentBondId(0))
        );
        assert_eq!(
            state.stereo_atom(StereoAtomHandle::New(0)),
            Ok(StereoAtomId(0))
        );
        assert_eq!(
            state.stereo_bond(StereoBondHandle::New(0)),
            Ok(StereoBondId(0))
        );

        if journaled {
            assert_eq!(
                undos,
                vec![
                    Undo::RemoveAddedTopology {
                        atoms: vec![
                            AddedAtom {
                                id: AtomId(2),
                                attributes: AtomForm::from_element(Element::O)
                            },
                            AddedAtom {
                                id: AtomId(3),
                                attributes: AtomForm::from_element(Element::F)
                            },
                        ],
                        bonds: vec![],
                    },
                    Undo::RemoveAddedTopology {
                        atoms: vec![],
                        bonds: vec![
                            AddedBond {
                                id: BondId(1),
                                endpoints: [AtomId(1), AtomId(2)],
                                attributes: BondForm::from_order(2)
                            },
                            AddedBond {
                                id: BondId(2),
                                endpoints: [AtomId(0), AtomId(3)],
                                attributes: BondForm::from_order(1)
                            },
                            AddedBond {
                                id: BondId(3),
                                endpoints: [AtomId(2), AtomId(3)],
                                attributes: BondForm::from_order(1)
                            },
                        ],
                    },
                    Undo::RemoveAddedDativeBond(AddedDativeBond {
                        id: DativeBondId(0),
                        donors: vec![AtomId(3), AtomId(0)],
                        acceptor: AtomId(2),
                        attributes: DativeBondForm::from_order(1)
                    }),
                    Undo::RemoveAddedAromaticSystem(AddedAromaticSystem {
                        id: AromaticSystemId(0),
                        atoms: vec![AtomId(2), AtomId(1), AtomId(0)],
                        attributes: AromaticSystemForm::from_electrons(vec![2, 1, 1])
                    }),
                    Undo::RemoveAddedMulticenterBond(AddedMulticenterBond {
                        id: MulticenterBondId(0),
                        atoms: vec![AtomId(3), AtomId(2)],
                        attributes: MulticenterBondForm::from_electrons(vec![1, 2])
                    }),
                    Undo::RemoveAddedNoncovalentBond(AddedNoncovalentBond {
                        id: NoncovalentBondId(0),
                        atoms: [AtomId(3), AtomId(0)],
                        attributes: NoncovalentBondForm::from_kind(
                            NoncovalentBondKind::HydrogenBond
                        )
                    }),
                    Undo::RemoveAddedStereoAtom(AddedStereoAtom {
                        id: StereoAtomId(0),
                        site: AtomId(0),
                        ligands: stereo_atom_ligands,
                        attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32)
                    }),
                    Undo::RemoveAddedStereoBond(AddedStereoBond {
                        id: StereoBondId(0),
                        site: BondId(1),
                        ligands: stereo_bond_ligands,
                        attributes: StereoBondForm::new(StereoKind::CisTrans, 1_u32)
                    }),
                ]
            );
            for undo in undos.into_iter().rev() {
                molecule.apply_undo(undo);
            }
            assert_eq!(molecule, initial);
        }
    }

    #[rstest]
    #[case::first(0)]
    #[case::middle(1)]
    #[case::last(2)]
    fn test_molecule_apply_edit_add_bonds_error(
        #[case] invalid_position: usize,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            ..Default::default()
        });
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let mut bonds = vec![
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                attributes: BondForm::from_order(1),
            },
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
                attributes: BondForm::from_order(2),
            },
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(2)), AtomHandle::Id(AtomId(0))],
                attributes: BondForm::from_order(1),
            },
        ];
        bonds[invalid_position].endpoints[1] = AtomHandle::Id(AtomId(3));
        let edit = Edit::AddBonds { bonds };
        let error = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).unwrap_err()
        } else {
            molecule.apply_edit(edit, &mut state).unwrap_err()
        };

        assert_eq!(
            error,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 3,
                count: 3
            }
        );
        assert_eq!(molecule, initial);
        assert_eq!(
            state.bond(BondHandle::New(0)),
            Err(TransactionError::HandleOutOfRange {
                kind: EntityKind::Bond,
                index: 0,
                count: 0
            })
        );
    }

    #[fixture]
    fn removal_entries() -> MoleculeEntries {
        MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
            ],
            dative: vec![
                (vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1)),
                (vec![AtomId(2)], AtomId(3), DativeBondForm::from_order(1)),
            ],
            aromatic: vec![
                (
                    vec![AtomId(0), AtomId(1)],
                    AromaticSystemForm::from_electrons(vec![1, 2]),
                ),
                (
                    vec![AtomId(2), AtomId(3)],
                    AromaticSystemForm::from_electrons(vec![2, 1]),
                ),
            ],
            multicenter: vec![
                (
                    vec![AtomId(0), AtomId(1)],
                    MulticenterBondForm::from_electrons(vec![1, 2]),
                ),
                (
                    vec![AtomId(2), AtomId(3)],
                    MulticenterBondForm::from_electrons(vec![2, 1]),
                ),
            ],
            noncovalent: vec![
                (
                    [AtomId(0), AtomId(1)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                ),
                (
                    [AtomId(2), AtomId(3)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                ),
            ],
            stereo_atoms: vec![
                (
                    AtomId(0),
                    vec![
                        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    ],
                    StereoAtomForm::default(),
                ),
                (
                    AtomId(2),
                    vec![
                        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                    ],
                    StereoAtomForm::default(),
                ),
            ],
            stereo_bonds: vec![
                (
                    BondId(0),
                    vec![
                        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
                ),
                (
                    BondId(1),
                    vec![
                        StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                        StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                        StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, 1_u32),
                ),
            ],
            constraints: Default::default(),
        }
    }

    #[rstest]
    fn test_molecule_apply_edit_remove_topology(
        mut removal_entries: MoleculeEntries,
        #[values(false, true)] journaled: bool,
    ) {
        let removed_constraint = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
        let old_constraint = Constraint::StereoBond(
            StereoBondId(1),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        );
        let new_constraint = Constraint::StereoBond(
            StereoBondId(0),
            StereoKind::CisTrans,
            StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
        );
        removal_entries.constraints = vec![
            removed_constraint.clone(),
            old_constraint.clone(),
            removed_constraint.clone(),
        ]
        .into();
        let initial = Molecule::from_entries(removal_entries);
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        molecule
            .apply_edit(
                Edit::AddAtoms {
                    atoms: vec![AtomForm::from_element(Element::F)],
                },
                &mut state,
            )
            .unwrap();
        let before_removal = molecule.clone();
        let edit = Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: vec![],
        };
        let undo = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).unwrap()
        } else {
            molecule.apply_edit(edit, &mut state).unwrap();
            None
        };
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
            ],
            bonds: vec![(AtomId(1), AtomId(2), BondForm::from_order(1))],
            dative: vec![(vec![AtomId(1)], AtomId(2), DativeBondForm::from_order(1))],
            aromatic: vec![(
                vec![AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            )],
            multicenter: vec![(
                vec![AtomId(1), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            )],
            noncovalent: vec![(
                [AtomId(1), AtomId(2)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
            stereo_atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::default(),
            )],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, 1_u32),
            )],
            constraints: vec![new_constraint.clone()].into(),
        });
        assert_eq!(molecule, expected);
        assert_eq!(
            state.atom(AtomHandle::Id(AtomId(0))),
            Err(TransactionError::HandleRemoved {
                kind: EntityKind::Atom,
                index: 0,
            })
        );
        assert_eq!(state.atom(AtomHandle::Id(AtomId(3))), Ok(AtomId(2)));
        assert_eq!(state.atom(AtomHandle::New(0)), Ok(AtomId(3)));
        assert_eq!(state.bond(BondHandle::Id(BondId(1))), Ok(BondId(0)));
        assert_eq!(
            state.dative_bond(DativeBondHandle::Id(DativeBondId(1))),
            Ok(DativeBondId(0))
        );
        assert_eq!(
            state.aromatic_system(AromaticSystemHandle::Id(AromaticSystemId(1))),
            Ok(AromaticSystemId(0))
        );
        assert_eq!(
            state.multicenter_bond(MulticenterBondHandle::Id(MulticenterBondId(1))),
            Ok(MulticenterBondId(0))
        );
        assert_eq!(
            state.noncovalent_bond(NoncovalentBondHandle::Id(NoncovalentBondId(1))),
            Ok(NoncovalentBondId(0))
        );
        assert_eq!(
            state.stereo_atom(StereoAtomHandle::Id(StereoAtomId(1))),
            Ok(StereoAtomId(0))
        );
        assert_eq!(
            state.stereo_bond(StereoBondHandle::Id(StereoBondId(1))),
            Ok(StereoBondId(0))
        );
        if journaled {
            let undo = undo.unwrap();
            let Undo::RestoreRemovedTopology { cascade, .. } = &undo else {
                panic!("expected topology restoration");
            };
            assert_eq!(
                cascade,
                &CascadedConstraints {
                    removed: vec![
                        RemovedConstraint {
                            position: 0,
                            constraint: removed_constraint.clone()
                        },
                        RemovedConstraint {
                            position: 2,
                            constraint: removed_constraint
                        },
                    ],
                    modified: vec![ModifiedConstraint {
                        position: 1,
                        old: old_constraint,
                        new: new_constraint
                    }],
                }
            );
            molecule.apply_undo(undo);
            assert_eq!(molecule, before_removal);
        }
    }

    #[rstest]
    #[case::dative(EntityKind::DativeBond)]
    #[case::aromatic(EntityKind::AromaticSystem)]
    #[case::multicenter(EntityKind::MulticenterBond)]
    #[case::noncovalent(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_molecule_apply_edit_remove_overlays(
        mut removal_entries: MoleculeEntries,
        #[case] kind: EntityKind,
        #[values(false, true)] journaled: bool,
    ) {
        let (edit, removed_constraint, old_constraint) = match kind {
            EntityKind::DativeBond => (
                Edit::RemoveDativeBonds {
                    removes: vec![(
                        DativeBondHandle::Id(DativeBondId(0)),
                        vec![AtomHandle::Id(AtomId(0))],
                        AtomHandle::Id(AtomId(1)),
                        DativeBondForm::from_order(1),
                    )],
                },
                Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(false)),
                Constraint::DativeBond(DativeBondId(1), DativeBondConstraintForm::aromatic(false)),
            ),
            EntityKind::AromaticSystem => (
                Edit::RemoveAromaticSystems {
                    removes: vec![(
                        AromaticSystemHandle::Id(AromaticSystemId(0)),
                        vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
                        AromaticSystemForm::from_electrons(vec![2, 1]),
                    )],
                },
                Constraint::AromaticSystem(
                    AromaticSystemId(0),
                    AromaticSystemConstraintForm::electron_count(3),
                ),
                Constraint::AromaticSystem(
                    AromaticSystemId(1),
                    AromaticSystemConstraintForm::electron_count(3),
                ),
            ),
            EntityKind::MulticenterBond => (
                Edit::RemoveMulticenterBonds {
                    removes: vec![(
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
                        MulticenterBondForm::from_electrons(vec![2, 1]),
                    )],
                },
                Constraint::MulticenterBond(
                    MulticenterBondId(0),
                    MulticenterBondConstraintForm::ElectronCount(3.into()),
                ),
                Constraint::MulticenterBond(
                    MulticenterBondId(1),
                    MulticenterBondConstraintForm::ElectronCount(3.into()),
                ),
            ),
            EntityKind::NoncovalentBond => (
                Edit::RemoveNoncovalentBonds {
                    removes: vec![(
                        NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                        [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    )],
                },
                Constraint::NoncovalentBond(
                    NoncovalentBondId(0),
                    NoncovalentBondConstraintForm::intramolecular(true),
                ),
                Constraint::NoncovalentBond(
                    NoncovalentBondId(1),
                    NoncovalentBondConstraintForm::intramolecular(true),
                ),
            ),
            EntityKind::StereoAtom => (
                Edit::RemoveStereoAtoms {
                    removes: vec![(
                        StereoAtomHandle::Id(StereoAtomId(0)),
                        AtomHandle::Id(AtomId(0)),
                        vec![
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                            (
                                AtomHandle::Id(AtomId(0)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                        ],
                        StereoAtomForm::default(),
                    )],
                },
                Constraint::Relational(RelationalConstraint::StereoAtomAllLigands {
                    stereo_atom: StereoAtomId(0),
                    predicate: Box::new(AtomConstraintForm::degree(1)),
                }),
                Constraint::Relational(RelationalConstraint::StereoAtomAllLigands {
                    stereo_atom: StereoAtomId(1),
                    predicate: Box::new(AtomConstraintForm::degree(1)),
                }),
            ),
            EntityKind::StereoBond => (
                Edit::RemoveStereoBonds {
                    removes: vec![(
                        StereoBondHandle::Id(StereoBondId(0)),
                        BondHandle::Id(BondId(0)),
                        vec![
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::LonePair),
                            (
                                AtomHandle::Id(AtomId(0)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (
                                AtomHandle::Id(AtomId(1)),
                                StereoLigandKind::ImplicitHydrogen,
                            ),
                            (AtomHandle::Id(AtomId(1)), StereoLigandKind::LonePair),
                        ],
                        StereoBondForm::new(StereoKind::CisTrans, 0_u32),
                    )],
                },
                Constraint::StereoBond(
                    StereoBondId(0),
                    StereoKind::CisTrans,
                    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
                ),
                Constraint::StereoBond(
                    StereoBondId(1),
                    StereoKind::CisTrans,
                    StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
                ),
            ),
            EntityKind::Atom | EntityKind::Bond => unreachable!(),
        };
        let unchanged = Constraint::Atom(AtomId(3), AtomConstraintForm::degree(1));
        removal_entries.constraints = vec![
            removed_constraint.clone(),
            unchanged.clone(),
            old_constraint.clone(),
            removed_constraint.clone(),
        ]
        .into();
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let undo = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).unwrap()
        } else {
            molecule.apply_edit(edit, &mut state).unwrap();
            None
        };
        match kind {
            EntityKind::DativeBond => {
                removal_entries.dative.remove(0);
            }
            EntityKind::AromaticSystem => {
                removal_entries.aromatic.remove(0);
            }
            EntityKind::MulticenterBond => {
                removal_entries.multicenter.remove(0);
            }
            EntityKind::NoncovalentBond => {
                removal_entries.noncovalent.remove(0);
            }
            EntityKind::StereoAtom => {
                removal_entries.stereo_atoms.remove(0);
            }
            EntityKind::StereoBond => {
                removal_entries.stereo_bonds.remove(0);
            }
            EntityKind::Atom | EntityKind::Bond => unreachable!(),
        }
        removal_entries.constraints = vec![unchanged, removed_constraint.clone()].into();
        assert_eq!(molecule, Molecule::from_entries(removal_entries));
        assert_eq!(state.atom(AtomHandle::Id(AtomId(3))), Ok(AtomId(3)));
        assert_eq!(state.bond(BondHandle::Id(BondId(1))), Ok(BondId(1)));
        for (current_kind, removed, survivor) in [
            (
                EntityKind::DativeBond,
                state
                    .dative_bond(DativeBondHandle::Id(DativeBondId(0)))
                    .map(|id| id.index()),
                state
                    .dative_bond(DativeBondHandle::Id(DativeBondId(1)))
                    .map(|id| id.index()),
            ),
            (
                EntityKind::AromaticSystem,
                state
                    .aromatic_system(AromaticSystemHandle::Id(AromaticSystemId(0)))
                    .map(|id| id.index()),
                state
                    .aromatic_system(AromaticSystemHandle::Id(AromaticSystemId(1)))
                    .map(|id| id.index()),
            ),
            (
                EntityKind::MulticenterBond,
                state
                    .multicenter_bond(MulticenterBondHandle::Id(MulticenterBondId(0)))
                    .map(|id| id.index()),
                state
                    .multicenter_bond(MulticenterBondHandle::Id(MulticenterBondId(1)))
                    .map(|id| id.index()),
            ),
            (
                EntityKind::NoncovalentBond,
                state
                    .noncovalent_bond(NoncovalentBondHandle::Id(NoncovalentBondId(0)))
                    .map(|id| id.index()),
                state
                    .noncovalent_bond(NoncovalentBondHandle::Id(NoncovalentBondId(1)))
                    .map(|id| id.index()),
            ),
            (
                EntityKind::StereoAtom,
                state
                    .stereo_atom(StereoAtomHandle::Id(StereoAtomId(0)))
                    .map(|id| id.index()),
                state
                    .stereo_atom(StereoAtomHandle::Id(StereoAtomId(1)))
                    .map(|id| id.index()),
            ),
            (
                EntityKind::StereoBond,
                state
                    .stereo_bond(StereoBondHandle::Id(StereoBondId(0)))
                    .map(|id| id.index()),
                state
                    .stereo_bond(StereoBondHandle::Id(StereoBondId(1)))
                    .map(|id| id.index()),
            ),
        ] {
            if current_kind == kind {
                assert_eq!(
                    removed,
                    Err(TransactionError::HandleRemoved { kind, index: 0 })
                );
                assert_eq!(survivor, Ok(0));
            } else {
                assert_eq!(removed, Ok(0));
                assert_eq!(survivor, Ok(1));
            }
        }
        if journaled {
            let undo = undo.unwrap();
            let cascade = match &undo {
                Undo::RestoreRemovedDativeBonds { cascade, .. }
                | Undo::RestoreRemovedAromaticSystems { cascade, .. }
                | Undo::RestoreRemovedMulticenterBonds { cascade, .. }
                | Undo::RestoreRemovedNoncovalentBonds { cascade, .. }
                | Undo::RestoreRemovedStereoAtoms { cascade, .. }
                | Undo::RestoreRemovedStereoBonds { cascade, .. } => cascade,
                _ => panic!("expected overlay restoration"),
            };
            assert_eq!(
                cascade,
                &CascadedConstraints {
                    removed: vec![
                        RemovedConstraint {
                            position: 0,
                            constraint: removed_constraint.clone()
                        },
                        RemovedConstraint {
                            position: 3,
                            constraint: removed_constraint.clone()
                        },
                    ],
                    modified: vec![ModifiedConstraint {
                        position: 2,
                        old: old_constraint,
                        new: removed_constraint
                    }],
                }
            );
            molecule.apply_undo(undo);
            assert_eq!(molecule, initial);
        }
    }

    #[rstest]
    #[case::duplicate_atom(Edit::RemoveTopology {
        atoms: vec![AtomHandle::Id(AtomId(0)); 2], bonds: vec![],
    }, TransactionError::DuplicateRemoval { kind: EntityKind::Atom })]
    #[case::duplicate_bond(Edit::RemoveTopology {
        atoms: vec![], bonds: vec![BondHandle::Id(BondId(0)); 2],
    }, TransactionError::DuplicateRemoval { kind: EntityKind::Bond })]
    #[case::missing_atom(Edit::RemoveTopology {
        atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(4))], bonds: vec![],
    }, TransactionError::HandleOutOfRange { kind: EntityKind::Atom, index: 4, count: 4 })]
    #[case::missing_bond(Edit::RemoveTopology {
        atoms: vec![AtomHandle::Id(AtomId(0))], bonds: vec![BondHandle::Id(BondId(2))],
    }, TransactionError::HandleOutOfRange { kind: EntityKind::Bond, index: 2, count: 2 })]
    #[case::duplicate_dative_bond(Edit::RemoveDativeBonds { removes: vec![(
        DativeBondHandle::Id(DativeBondId(0)), vec![AtomHandle::Id(AtomId(0))],
        AtomHandle::Id(AtomId(1)), DativeBondForm::from_order(1),
    ); 2] }, TransactionError::DuplicateRemoval { kind: EntityKind::DativeBond })]
    #[case::dative_bond_old_state(Edit::RemoveDativeBonds { removes: vec![
        (DativeBondHandle::Id(DativeBondId(0)), vec![AtomHandle::Id(AtomId(0))],
         AtomHandle::Id(AtomId(1)), DativeBondForm::from_order(1)),
        (DativeBondHandle::Id(DativeBondId(1)), vec![AtomHandle::Id(AtomId(2))],
         AtomHandle::Id(AtomId(3)), DativeBondForm::from_order(2)),
    ] }, TransactionError::OldStateMismatch)]
    fn test_molecule_apply_edit_remove_error(
        removal_entries: MoleculeEntries,
        #[case] edit: Edit,
        #[case] expected: TransactionError,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(removal_entries);
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let result = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(result, Err(expected));
        assert_eq!(molecule, initial);
        assert_eq!(state.atoms.initial, None);
        assert_eq!(state.bonds.initial, None);
        assert_eq!(state.dative_bonds.initial, None);
        assert_eq!(state.aromatic_systems.initial, None);
        assert_eq!(state.multicenter_bonds.initial, None);
        assert_eq!(state.noncovalent_bonds.initial, None);
        assert_eq!(state.stereo_atoms.initial, None);
        assert_eq!(state.stereo_bonds.initial, None);
    }

    #[rstest]
    #[case::identity(vec![])]
    #[case::site(vec![BondHandle::Id(BondId(0))])]
    fn test_molecule_apply_edit_remove_bonds(
        mut removal_entries: MoleculeEntries,
        #[case] bonds: Vec<BondHandle>,
        #[values(false, true)] journaled: bool,
    ) {
        removal_entries.stereo_atoms.clear();
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        if !bonds.is_empty() {
            removal_entries.bonds.remove(0);
            removal_entries.stereo_bonds.remove(0);
            removal_entries.stereo_bonds[0].0 = BondId(0);
        }
        let expected = Molecule::from_entries(removal_entries);
        let edit = Edit::RemoveTopology {
            atoms: vec![],
            bonds,
        };
        let undo = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).unwrap()
        } else {
            molecule.apply_edit(edit, &mut state).unwrap();
            None
        };
        assert_eq!(molecule, expected);
        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert_eq!(molecule, initial);
        }
    }

    #[rstest]
    #[case::element(
        AtomFieldChange::Element { old: ElementForm::Lit(Element::C), new: ElementForm::Lit(Element::N) },
        AtomForm::from_element(Element::N),
    )]
    #[case::element_equivalent(
        AtomFieldChange::Element { old: ElementForm::lit_set([Element::C]), new: ElementForm::Lit(Element::N) },
        AtomForm::from_element(Element::N),
    )]
    #[case::isotope_mass(
        AtomFieldChange::IsotopeMass { old: IsotopeMassForm::Undetermined, new: IsotopeMassForm::Lit(13) },
        AtomForm { isotope_mass: IsotopeMassForm::Lit(13), ..AtomForm::from_element(Element::C) },
    )]
    #[case::charge(
        AtomFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(-1) },
        AtomForm { charge: NumForm::Lit(-1), ..AtomForm::from_element(Element::C) },
    )]
    #[case::implicit_hydrogens(
        AtomFieldChange::ImplicitHydrogens { old: NumForm::Undetermined, new: NumForm::Lit(2) },
        AtomForm { implicit_hydrogens: NumForm::Lit(2), ..AtomForm::from_element(Element::C) },
    )]
    #[case::lone_pairs(
        AtomFieldChange::LonePairs { old: NumForm::Undetermined, new: NumForm::Lit(1) },
        AtomForm { lone_pairs: NumForm::Lit(1), ..AtomForm::from_element(Element::C) },
    )]
    #[case::unpaired_electrons(
        AtomFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
        AtomForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..AtomForm::from_element(Element::C) },
    )]
    fn test_molecule_apply_edit_modify_atom_field(
        mut removal_entries: MoleculeEntries,
        #[case] change: AtomFieldChange,
        #[case] expected: AtomForm,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let edit = Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change,
        };
        let undo = if journaled {
            molecule
                .apply_edit_with_undo(edit.clone(), &mut state)
                .unwrap()
        } else {
            molecule.apply_edit(edit.clone(), &mut state).unwrap();
            None
        };
        removal_entries.atoms[0] = expected;
        let expected = Molecule::from_entries(removal_entries);
        assert_eq!(molecule, expected);

        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);

        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert!(molecule.normalized_eq(&initial));
        }
    }

    #[rstest]
    #[case::order(
        BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
        BondForm::from_order(2),
    )]
    #[case::order_equivalent(
        BondFieldChange::Order { old: NumForm::lit_set([1]), new: NumForm::Lit(2) },
        BondForm::from_order(2),
    )]
    #[case::charge(
        BondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) },
        BondForm { charge: NumForm::Lit(1), ..BondForm::from_order(1) },
    )]
    #[case::unpaired_electrons(
        BondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
        BondForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..BondForm::from_order(1) },
    )]
    fn test_molecule_apply_edit_modify_bond_field(
        mut removal_entries: MoleculeEntries,
        #[case] change: BondFieldChange,
        #[case] expected: BondForm,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let edit = Edit::ModifyBondField {
            id: BondHandle::Id(BondId(0)),
            change,
        };
        let undo = if journaled {
            molecule
                .apply_edit_with_undo(edit.clone(), &mut state)
                .unwrap()
        } else {
            molecule.apply_edit(edit.clone(), &mut state).unwrap();
            None
        };
        removal_entries.bonds[0].2 = expected;
        let expected = Molecule::from_entries(removal_entries);
        assert_eq!(molecule, expected);

        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);

        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert!(molecule.normalized_eq(&initial));
        }
    }

    #[rstest]
    #[case::order(
        DativeBondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
        DativeBondForm::from_order(2),
    )]
    #[case::order_equivalent(
        DativeBondFieldChange::Order { old: NumForm::lit_set([1]), new: NumForm::Lit(2) },
        DativeBondForm::from_order(2),
    )]
    fn test_molecule_apply_edit_modify_dative_bond_field(
        mut removal_entries: MoleculeEntries,
        #[case] change: DativeBondFieldChange,
        #[case] expected: DativeBondForm,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let edit = Edit::ModifyDativeBondField {
            id: DativeBondHandle::Id(DativeBondId(0)),
            change,
        };
        let undo = if journaled {
            molecule
                .apply_edit_with_undo(edit.clone(), &mut state)
                .unwrap()
        } else {
            molecule.apply_edit(edit.clone(), &mut state).unwrap();
            None
        };
        removal_entries.dative[0].2 = expected;
        let expected = Molecule::from_entries(removal_entries);
        assert_eq!(molecule, expected);

        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);

        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert!(molecule.normalized_eq(&initial));
        }
    }

    #[rstest]
    #[case::electrons(
        AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2, 2]) },
        AromaticSystemForm::from_electrons(vec![2, 2]),
    )]
    #[case::electron_count_length(
        AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2]) },
        AromaticSystemForm::from_electrons(vec![2]),
    )]
    #[case::charge(
        AromaticSystemFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(1) },
        AromaticSystemForm { charge: NumForm::Lit(1), ..AromaticSystemForm::from_electrons(vec![1, 2]) },
    )]
    #[case::unpaired_electrons(
        AromaticSystemFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
        AromaticSystemForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..AromaticSystemForm::from_electrons(vec![1, 2]) },
    )]
    fn test_molecule_apply_edit_modify_aromatic_system_field(
        mut removal_entries: MoleculeEntries,
        #[case] change: AromaticSystemFieldChange,
        #[case] expected: AromaticSystemForm,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let edit = Edit::ModifyAromaticSystemField {
            id: AromaticSystemHandle::Id(AromaticSystemId(0)),
            change,
        };
        let undo = if journaled {
            molecule
                .apply_edit_with_undo(edit.clone(), &mut state)
                .unwrap()
        } else {
            molecule.apply_edit(edit.clone(), &mut state).unwrap();
            None
        };
        removal_entries.aromatic[0].1 = expected;
        let expected = Molecule::from_entries(removal_entries);
        assert_eq!(molecule, expected);

        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);

        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert!(molecule.normalized_eq(&initial));
        }
    }

    #[rstest]
    #[case::electrons(
        MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2, 2]) },
        MulticenterBondForm::from_electrons(vec![2, 2]),
    )]
    #[case::electron_count_length(
        MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 2]), new: ElectronCountsForm::Lit(vec![2, 1, 2]) },
        MulticenterBondForm::from_electrons(vec![2, 1, 2]),
    )]
    #[case::charge(
        MulticenterBondFieldChange::Charge { old: NumForm::Undetermined, new: NumForm::Lit(-1) },
        MulticenterBondForm { charge: NumForm::Lit(-1), ..MulticenterBondForm::from_electrons(vec![1, 2]) },
    )]
    #[case::unpaired_electrons(
        MulticenterBondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::default(), new: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } },
        MulticenterBondForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, ..MulticenterBondForm::from_electrons(vec![1, 2]) },
    )]
    fn test_molecule_apply_edit_modify_multicenter_bond_field(
        mut removal_entries: MoleculeEntries,
        #[case] change: MulticenterBondFieldChange,
        #[case] expected: MulticenterBondForm,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let edit = Edit::ModifyMulticenterBondField {
            id: MulticenterBondHandle::Id(MulticenterBondId(0)),
            change,
        };
        let undo = if journaled {
            molecule
                .apply_edit_with_undo(edit.clone(), &mut state)
                .unwrap()
        } else {
            molecule.apply_edit(edit.clone(), &mut state).unwrap();
            None
        };
        removal_entries.multicenter[0].1 = expected;
        let expected = Molecule::from_entries(removal_entries);
        assert_eq!(molecule, expected);

        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);

        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert!(molecule.normalized_eq(&initial));
        }
    }

    #[rstest]
    #[case::kind(
        NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic) },
        NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
    )]
    fn test_molecule_apply_edit_modify_noncovalent_bond_field(
        mut removal_entries: MoleculeEntries,
        #[case] change: NoncovalentBondFieldChange,
        #[case] expected: NoncovalentBondForm,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let edit = Edit::ModifyNoncovalentBondField {
            id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            change,
        };
        let undo = if journaled {
            molecule
                .apply_edit_with_undo(edit.clone(), &mut state)
                .unwrap()
        } else {
            molecule.apply_edit(edit.clone(), &mut state).unwrap();
            None
        };
        removal_entries.noncovalent[0].1 = expected;
        let expected = Molecule::from_entries(removal_entries);
        assert_eq!(molecule, expected);

        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);

        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert!(molecule.normalized_eq(&initial));
        }
    }

    #[rstest]
    #[case::configuration(
        StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::Undetermined, new: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)) },
        StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
    )]
    #[case::coset_range(
        StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::Undetermined, new: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(2)) },
        StereoAtomForm::new(StereoKind::Tetrahedral, 2_u32),
    )]
    fn test_molecule_apply_edit_modify_stereo_atom_field(
        mut removal_entries: MoleculeEntries,
        #[case] change: StereoAtomFieldChange,
        #[case] expected: StereoAtomForm,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let edit = Edit::ModifyStereoAtomField {
            id: StereoAtomHandle::Id(StereoAtomId(0)),
            change,
        };
        let undo = if journaled {
            molecule
                .apply_edit_with_undo(edit.clone(), &mut state)
                .unwrap()
        } else {
            molecule.apply_edit(edit.clone(), &mut state).unwrap();
            None
        };
        removal_entries.stereo_atoms[0].2 = expected;
        let expected = Molecule::from_entries(removal_entries);
        assert_eq!(molecule, expected);

        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);

        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert!(molecule.normalized_eq(&initial));
        }
    }

    #[rstest]
    #[case::configuration(
        StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)), new: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0)) },
        StereoBondForm::new(StereoKind::CisTrans, 0_u32),
    )]
    #[case::configuration_equivalent(
        StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::lit_set([1])), new: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0)) },
        StereoBondForm::new(StereoKind::CisTrans, 0_u32),
    )]
    #[case::coset_range(
        StereoBondFieldChange::Configuration { old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)), new: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(2)) },
        StereoBondForm::new(StereoKind::CisTrans, 2_u32),
    )]
    fn test_molecule_apply_edit_modify_stereo_bond_field(
        mut removal_entries: MoleculeEntries,
        #[case] change: StereoBondFieldChange,
        #[case] expected: StereoBondForm,
        #[values(false, true)] journaled: bool,
    ) {
        let initial = Molecule::from_entries(removal_entries.clone());
        let mut molecule = initial.clone();
        let mut state = ApplicationState::new(&molecule);
        let edit = Edit::ModifyStereoBondField {
            id: StereoBondHandle::Id(StereoBondId(0)),
            change,
        };
        let undo = if journaled {
            molecule
                .apply_edit_with_undo(edit.clone(), &mut state)
                .unwrap()
        } else {
            molecule.apply_edit(edit.clone(), &mut state).unwrap();
            None
        };
        removal_entries.stereo_bonds[0].2 = expected;
        let expected = Molecule::from_entries(removal_entries);
        assert_eq!(molecule, expected);

        let repeated = if journaled {
            molecule.apply_edit_with_undo(edit, &mut state).map(|_| ())
        } else {
            molecule.apply_edit(edit, &mut state)
        };
        assert_eq!(repeated, Err(TransactionError::OldStateMismatch));
        assert_eq!(molecule, expected);

        if journaled {
            molecule.apply_undo(undo.unwrap());
            assert!(molecule.normalized_eq(&initial));
        }
    }

    #[rstest]
    fn test_molecule_apply_undo_added_topology() {
        let initial = MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
                AtomForm::from_element(Element::O),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(2)),
            ],
            ..Default::default()
        };
        let mut extended = initial.clone();
        extended.atoms.push(AtomForm::from_element(Element::F));
        extended
            .bonds
            .push((AtomId(2), AtomId(3), BondForm::from_order(1)));
        let mut molecule = Molecule::from_entries(extended);

        molecule.apply_undo(Undo::RemoveAddedTopology {
            atoms: vec![AddedAtom {
                id: AtomId(3),
                attributes: AtomForm::from_element(Element::F),
            }],
            bonds: vec![AddedBond {
                id: BondId(2),
                endpoints: [AtomId(2), AtomId(3)],
                attributes: BondForm::from_order(1),
            }],
        });

        assert_eq!(molecule, Molecule::from_entries(initial));
    }

    #[rstest]
    #[case::topology(Undo::RemoveAddedTopology {
        atoms: vec![AddedAtom { id: AtomId(0), attributes: AtomForm::default() }],
        bonds: vec![AddedBond { id: BondId(0), endpoints: [AtomId(0), AtomId(1)], attributes: BondForm::default() }],
    })]
    #[case::dative_bond(Undo::RemoveAddedDativeBond(AddedDativeBond { id: DativeBondId(0), donors: vec![], acceptor: AtomId(0), attributes: DativeBondForm::default() }))]
    #[case::aromatic_system(Undo::RemoveAddedAromaticSystem(AddedAromaticSystem { id: AromaticSystemId(0), atoms: vec![], attributes: AromaticSystemForm::default() }))]
    #[case::multicenter_bond(Undo::RemoveAddedMulticenterBond(AddedMulticenterBond { id: MulticenterBondId(0), atoms: vec![], attributes: MulticenterBondForm::default() }))]
    #[case::noncovalent_bond(Undo::RemoveAddedNoncovalentBond(AddedNoncovalentBond { id: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(1)], attributes: NoncovalentBondForm::default() }))]
    #[case::stereo_atom(Undo::RemoveAddedStereoAtom(AddedStereoAtom { id: StereoAtomId(0), site: AtomId(0), ligands: vec![], attributes: StereoAtomForm::default() }))]
    #[case::stereo_bond(Undo::RemoveAddedStereoBond(AddedStereoBond { id: StereoBondId(0), site: BondId(0), ligands: vec![], attributes: StereoBondForm::default() }))]
    fn test_molecule_apply_undo_added_entry_manipulated(#[case] undo: Undo) {
        Molecule::default().apply_undo(undo);
    }

    #[rstest]
    #[case::atom(Undo::ModifyAtomField {
        id: AtomId(0),
        change: AtomFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(0) },
    })]
    #[case::bond(Undo::ModifyBondField {
        id: BondId(0),
        change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(1) },
    })]
    #[case::dative_bond(Undo::ModifyDativeBondField {
        id: DativeBondId(0),
        change: DativeBondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(1) },
    })]
    #[case::aromatic_system(Undo::ModifyAromaticSystemField {
        id: AromaticSystemId(0),
        change: AromaticSystemFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(0) },
    })]
    #[case::multicenter_bond(Undo::ModifyMulticenterBondField {
        id: MulticenterBondId(0),
        change: MulticenterBondFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(0) },
    })]
    #[case::noncovalent_bond(Undo::ModifyNoncovalentBondField {
        id: NoncovalentBondId(0),
        change: NoncovalentBondFieldChange::Kind {
            old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
        },
    })]
    #[case::stereo_atom(Undo::ModifyStereoAtomField {
        id: StereoAtomId(0),
        change: StereoAtomFieldChange::Configuration {
            old: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            new: StereoConfigurationForm::Undetermined,
        },
    })]
    #[case::stereo_bond(Undo::ModifyStereoBondField {
        id: StereoBondId(0),
        change: StereoBondFieldChange::Configuration {
            old: StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
            new: StereoConfigurationForm::Undetermined,
        },
    })]
    fn test_molecule_apply_undo_field_manipulated(#[case] undo: Undo) {
        Molecule::default().apply_undo(undo);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame(AtomId(0), vec![AtomId(1), AtomId(2)], true)]
    #[case::reordered_donors(AtomId(0), vec![AtomId(2), AtomId(1)], true)]
    #[case::different_acceptor(AtomId(1), vec![AtomId(1), AtomId(2)], false)]
    #[case::different_donors(AtomId(0), vec![AtomId(1), AtomId(3)], false)]
    fn test_molecule_dative_bond_equiv(
        #[case] acceptor: AtomId,
        #[case] donors: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let mut molecule = Molecule::default();
        for _ in 0..4 {
            molecule.add_atom(AtomForm::from_element(Element::C));
        }
        molecule.add_dative_bond(
            &[AtomId(1), AtomId(2)],
            AtomId(0),
            DativeBondForm::from_order(1),
        );
        let offered = DativeBondForm::from_order(1);

        assert_eq!(molecule.dative_bond_equiv(DativeBondId(0), acceptor, &donors, &offered), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame(vec![AtomId(0), AtomId(1), AtomId(2)], vec![10, 20, 30], true)]
    #[case::reordered_frame_carrying_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![30, 10, 20], true)]
    #[case::reordered_frame_keeping_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![10, 20, 30], false)]
    #[case::different_counts(vec![AtomId(0), AtomId(1), AtomId(2)], vec![10, 20, 99], false)]
    #[case::multiset_differs(vec![AtomId(0), AtomId(1), AtomId(3)], vec![10, 20, 30], false)]
    #[case::wrong_arity(vec![AtomId(0), AtomId(1)], vec![10, 20], false)]
    fn test_molecule_aromatic_system_equiv(
        #[case] atoms: Vec<AtomId>,
        #[case] electrons: Vec<i64>,
        #[case] expected: bool,
    ) {
        let mut molecule = Molecule::default();
        for _ in 0..4 {
            molecule.add_atom(AtomForm::from_element(Element::C));
        }
        molecule.add_aromatic_system(
            &[AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![10, 20, 30]),
        );
        let offered = AromaticSystemForm::from_electrons(electrons);

        assert_eq!(molecule.aromatic_system_equiv(AromaticSystemId(0), &atoms, &offered), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_atoms(vec![AtomId(0), AtomId(1), AtomId(2)], true)]
    #[case::reordered_atoms(vec![AtomId(2), AtomId(0), AtomId(1)], true)]
    #[case::different_atoms(vec![AtomId(0), AtomId(1), AtomId(3)], false)]
    #[case::wrong_arity(vec![AtomId(0), AtomId(1)], false)]
    fn test_molecule_aromatic_system_equiv_undetermined_electrons(
        #[case] atoms: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let mut molecule = Molecule::default();
        for _ in 0..4 {
            molecule.add_atom(AtomForm::from_element(Element::C));
        }
        molecule.add_aromatic_system(
            &[AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        );

        assert_eq!(
            molecule.aromatic_system_equiv(
                AromaticSystemId(0),
                &atoms,
                &AromaticSystemForm::default(),
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame(vec![AtomId(0), AtomId(1), AtomId(2)], vec![10, 20, 30], true)]
    #[case::reordered_frame_carrying_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![30, 10, 20], true)]
    #[case::reordered_frame_keeping_its_counts(vec![AtomId(2), AtomId(0), AtomId(1)], vec![10, 20, 30], false)]
    #[case::multiset_differs(vec![AtomId(0), AtomId(1), AtomId(3)], vec![10, 20, 30], false)]
    fn test_molecule_multicenter_bond_equiv(
        #[case] atoms: Vec<AtomId>,
        #[case] electrons: Vec<i64>,
        #[case] expected: bool,
    ) {
        let mut molecule = Molecule::default();
        for _ in 0..4 {
            molecule.add_atom(AtomForm::from_element(Element::C));
        }
        molecule.add_multicenter_bond(
            &[AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::from_electrons(vec![10, 20, 30]),
        );
        let offered = MulticenterBondForm::from_electrons(electrons);

        assert_eq!(molecule.multicenter_bond_equiv(MulticenterBondId(0), &atoms, &offered), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_atoms(vec![AtomId(0), AtomId(1), AtomId(2)], true)]
    #[case::different_atoms(vec![AtomId(0), AtomId(1), AtomId(3)], false)]
    fn test_molecule_multicenter_bond_equiv_undetermined_electrons(
        #[case] atoms: Vec<AtomId>,
        #[case] expected: bool,
    ) {
        let mut molecule = Molecule::default();
        for _ in 0..4 {
            molecule.add_atom(AtomForm::from_element(Element::C));
        }
        molecule.add_multicenter_bond(
            &[AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        );

        assert_eq!(
            molecule.multicenter_bond_equiv(
                MulticenterBondId(0),
                &atoms,
                &MulticenterBondForm::default(),
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame([AtomId(0), AtomId(1)], true)]
    #[case::reversed_frame([AtomId(1), AtomId(0)], true)]
    #[case::different_pair([AtomId(0), AtomId(2)], false)]
    fn test_molecule_noncovalent_bond_equiv(
        #[case] atoms: [AtomId; 2],
        #[case] expected: bool,
    ) {
        let mut molecule = Molecule::default();
        for _ in 0..3 {
            molecule.add_atom(AtomForm::from_element(Element::C));
        }
        molecule.add_noncovalent_bond(
            [AtomId(0), AtomId(1)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        let offered = NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond);

        assert_eq!(molecule.noncovalent_bond_equiv(NoncovalentBondId(0), atoms, &offered), expected);
    }

    #[fixture]
    fn stereo_molecule() -> Molecule {
        let mut b = Molecule::default();
        b.add_atom(AtomForm::from_element(Element::C));
        for element in [Element::F, Element::Cl, Element::Br, Element::I] {
            b.add_atom(AtomForm::from_element(element));
        }
        for ligand in 1..=4 {
            b.add_bond(AtomId(0), AtomId(ligand), BondForm::from_order(1));
        }
        b.add_stereo_atom(
            AtomId(0),
            &(1..=4)
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect::<Vec<_>>(),
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        );
        b
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored_frame([1, 2, 3, 4], 0, true)]
    #[case::stored_frame_other_coset([1, 2, 3, 4], 1, false)]
    #[case::transposed_frame_same_coset([2, 1, 3, 4], 0, false)]
    #[case::transposed_frame_other_coset([2, 1, 3, 4], 1, true)]
    #[case::multiset_differs([1, 2, 3, 5], 0, false)]
    fn test_molecule_stereo_atom_equiv(
        stereo_molecule: Molecule,
        #[case] ligands: [u32; 4],
        #[case] coset: u32,
        #[case] expected: bool,
    ) {
        let offered: Vec<StereoLigand> = ligands
            .into_iter()
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let attributes = StereoAtomForm::new(StereoKind::Tetrahedral, coset);

        assert_eq!(
            stereo_molecule.stereo_atom_equiv(StereoAtomId(0), AtomId(0), &offered, &attributes),
            expected,
        );
    }

    #[rstest]
    fn test_molecule_stereo_atom_equiv_reordered_frame(stereo_molecule: Molecule) {
        let stored: Vec<StereoLigand> = (1..=4)
            .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
            .collect();
        let configuration = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);

        assert!(
            stereo_molecule.stereo_atom_equiv(StereoAtomId(0), AtomId(0), &stored, &configuration),
            "the stored frame with its own configuration is equivalent to itself",
        );

        let transposed = Permutation::from_image(&[1, 0, 2, 3]);
        assert!(
            !stereo_molecule.stereo_atom_equiv(
                StereoAtomId(0),
                AtomId(0),
                &transposed.act(&stored),
                &configuration,
            ),
            "coset 0 against a transposed frame is the opposite arrangement, not the stored one",
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::stored(BondId(0), vec![2, 3, 4, 5], true)]
    #[case::within_endpoint(BondId(0), vec![3, 2, 4, 5], true)]
    #[case::endpoint_block_swap(BondId(0), vec![4, 5, 2, 3], true)]
    #[case::across_endpoints(BondId(0), vec![2, 4, 3, 5], false)]
    #[case::different_ligand(BondId(0), vec![2, 3, 4, 6], false)]
    #[case::different_site(BondId(1), vec![2, 3, 4, 5], false)]
    fn test_molecule_stereo_bond_equiv(
        #[case] site: BondId,
        #[case] ligand_ids: Vec<u32>,
        #[case] expected: bool,
    ) {
        let mut molecule = Molecule::default();
        for _ in 0..7 {
            molecule.add_atom(AtomForm::from_element(Element::C));
        }
        for (first, second) in [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)] {
            molecule.add_bond(AtomId(first), AtomId(second), BondForm::from_order(1));
        }
        molecule.add_stereo_bond(
            BondId(0),
            &[2, 3, 4, 5].map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)),
            StereoBondForm::default(),
        );
        let ligands = ligand_ids
            .into_iter()
            .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
            .collect::<Vec<_>>();

        assert_eq!(
            molecule.stereo_bond_equiv(
                StereoBondId(0),
                site,
                &ligands,
                &StereoBondForm::default(),
            ),
            expected,
        );
    }
}
