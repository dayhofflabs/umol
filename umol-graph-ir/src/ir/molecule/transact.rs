//! Transactional `Edit` application on `MoleculeEditor`.
//!
//! `transact(edits)` applies each `Edit` in order, records realized `Undo`
//! entries, and either returns a rollback-capable `Transaction` or reverse-
//! replays the journal before surfacing a `TransactionError`.
//!
//! `Id(n)` handles retain the identity of entity `n` in the initial host, while
//! `New(n)` handles retain the identity of the `n`th same-kind creation in the
//! edit sequence. Transaction application alone realizes those stable handles
//! against a host and tracks their liveness across compaction.

use std::collections::HashSet;
use std::hash::Hash;

use thiserror::Error;

use super::super::constraint::{Constraint, Constraints};
use super::super::edit::{
    AddBond, AddedAromaticSystem, AddedAtom, AddedBond, AddedDativeBond, AddedMulticenterBond,
    AddedNoncovalentBond, AddedStereoAtom, AddedStereoBond, AromaticSystemFieldChange,
    AromaticSystemHandle, AtomFieldChange, AtomHandle, BondFieldChange, BondHandle,
    CascadedConstraints, ConstraintEdit, DativeBondFieldChange, DativeBondHandle, Edit, Edits,
    MulticenterBondFieldChange, MulticenterBondHandle, NoncovalentBondFieldChange,
    NoncovalentBondHandle, RemovedAromaticSystem, RemovedAtom, RemovedBond, RemovedConstraint,
    RemovedDativeBond, RemovedMulticenterBond, RemovedNoncovalentBond, RemovedOverlays,
    RemovedStereoAtom, RemovedStereoBond, StereoAtomFieldChange, StereoAtomHandle,
    StereoBondFieldChange, StereoBondHandle, Undo,
};
use super::super::entity::EntityKind;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::ligand::StereoLigand;
use super::super::remap::{IdCompaction, UndoCompaction};
use super::super::traits::Canonicalize;
use super::MoleculeEditor;

#[derive(Debug, Error, PartialEq, Eq, Clone)]
pub enum TransactionError {
    #[error("{kind} handle {index} is out of range for {count} entries")]
    HandleOutOfRange {
        kind: EntityKind,
        index: usize,
        count: usize,
    },

    #[error("{kind} handle {index} refers to a removed entity")]
    HandleRemoved { kind: EntityKind, index: usize },

    #[error("duplicate {kind} in removal batch")]
    DuplicateRemoval { kind: EntityKind },

    /// `Set*Field` or `Set*Constraint`: current state does not match the
    /// edit's `old` payload.
    #[error("precondition failed: old state does not match current")]
    OldStateMismatch,

    /// `Remove*Constraint` with a value that's not present.
    #[error("missing constraint entry on remove")]
    MissingEntry,

    /// Edit shape is structurally invalid (e.g., `AddDativeBond` with no
    /// participants).
    #[error("malformed edit: {0}")]
    MalformedEdit(&'static str),

    #[error("rollback failed after apply error: apply={apply}; rollback={rollback}")]
    RollbackFailed {
        apply: Box<TransactionError>,
        rollback: Box<TransactionError>,
    },

    /// The rollback journal cannot be structurally applied to the supplied editor state.
    ///
    /// A transaction guarantees exact restoration only when rolled back against the exact
    /// post-transaction state (or the end of the consecutive chain represented by an appended
    /// journal). Other states are rejected when a required receiver or reconstruction slot is
    /// absent; structurally compatible but unrelated states are outside that guarantee.
    #[error("rollback journal does not match editor state")]
    RollbackStateMismatch,
}

/// Detached journal of the realized undos for one successfully applied edit batch.
///
/// Detachment permits journals for consecutive transactions to be appended and rolled back as a
/// unit. Exact restoration is guaranteed only for the exact post-transaction editor state, or the
/// end state of the consecutively appended transaction chain. A structurally incompatible state
/// returns [`TransactionError::RollbackStateMismatch`] without panicking; a structurally compatible
/// but unrelated state is outside the semantic guarantee.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Transaction {
    undo: Vec<Undo>,
}

impl Transaction {
    pub fn undos(&self) -> &[Undo] {
        &self.undo
    }

    /// Append the rollback journal for a transaction applied after this one.
    pub fn append(&mut self, later: Self) {
        self.undo.extend(later.undo);
    }

    /// Reverse the journal against its exact post-transaction editor state.
    ///
    /// This exactly restores the pre-transaction state when the editor is the state produced by
    /// this transaction (or by its appended consecutive chain). Structural incompatibility returns
    /// [`TransactionError::RollbackStateMismatch`] rather than panicking. No semantic result is
    /// guaranteed when a different editor happens to satisfy the journal's structural requirements.
    pub fn rollback(self, editor: &mut MoleculeEditor) -> Result<(), TransactionError> {
        rollback_journal(editor, self.undo)
    }
}

#[derive(Default)]
struct HandleTable<I> {
    initial: Vec<Option<I>>,
    created: Vec<Option<I>>,
}

impl<I: Copy> HandleTable<I> {
    fn new(initial: impl IntoIterator<Item = I>) -> Self {
        Self {
            initial: initial.into_iter().map(Some).collect(),
            created: Vec::new(),
        }
    }

    fn initial(&self, kind: EntityKind, index: usize) -> Result<I, TransactionError> {
        Self::resolve(&self.initial, kind, index)
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
        for id in self.initial.iter_mut().chain(&mut self.created) {
            if let Some(current) = *id {
                *id = compact(current);
            }
        }
    }
}

/// Host-specific realization of the stable handles in one [`Edits`] sequence.
///
/// Initial-host and newly created ordinals remain in separate per-kind tables. Removals compact
/// every surviving concrete id in both tables and leave tombstones for removed entities; the
/// transaction driver alone owns this realization and compaction state.
struct ApplicationState {
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
    fn new(editor: &MoleculeEditor) -> Self {
        Self {
            atoms: HandleTable::new((0..editor.atom_count()).map(AtomId::from)),
            bonds: HandleTable::new((0..editor.bond_count()).map(BondId::from)),
            dative_bonds: HandleTable::new((0..editor.dative_bond_count()).map(DativeBondId::from)),
            aromatic_systems: HandleTable::new(
                (0..editor.aromatic_system_count()).map(AromaticSystemId::from),
            ),
            multicenter_bonds: HandleTable::new(
                (0..editor.multicenter_bond_count()).map(MulticenterBondId::from),
            ),
            noncovalent_bonds: HandleTable::new(
                (0..editor.noncovalent_bond_count()).map(NoncovalentBondId::from),
            ),
            stereo_atoms: HandleTable::new((0..editor.stereo_atom_count()).map(StereoAtomId::from)),
            stereo_bonds: HandleTable::new((0..editor.stereo_bond_count()).map(StereoBondId::from)),
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
        ligands: Vec<(AtomHandle, super::super::ligand::StereoLigandKind)>,
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

    fn compact(&mut self, compaction: &IdCompaction) {
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

impl MoleculeEditor {
    /// Apply an ordered [`Edits`] batch atomically. On success, returns a rollback
    /// transaction. On any apply failure, reverse-replays the already-created
    /// undo journal.
    pub fn transact(&mut self, edits: Edits) -> Result<Transaction, TransactionError> {
        let mut journal: Vec<Undo> = Vec::with_capacity(edits.len());
        let mut state = ApplicationState::new(self);
        for edit in edits {
            match self.apply_edit(edit, &mut state) {
                Ok(undo) => journal.push(undo),
                Err(apply) => {
                    if let Err(rollback) = rollback_journal(self, journal) {
                        return Err(TransactionError::RollbackFailed {
                            apply: Box::new(apply),
                            rollback: Box::new(rollback),
                        });
                    }
                    return Err(apply);
                }
            }
        }
        Ok(Transaction { undo: journal })
    }

    pub fn transact_unchecked(&mut self, edits: Edits) {
        let mut state = ApplicationState::new(self);
        for edit in edits {
            if let Err(e) = self.apply_edit_unchecked(edit, &mut state) {
                panic!("invalid unchecked transaction edit: {e}");
            }
        }
    }

    fn apply_edit_unchecked(
        &mut self,
        edit: Edit,
        state: &mut ApplicationState,
    ) -> Result<(), TransactionError> {
        match edit {
            Edit::AddAtoms { atoms } => {
                for atom in atoms {
                    let id = self.add_atom(atom);
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
                             ast,
                         }| {
                            Ok(([state.atom(first)?, state.atom(second)?], ast))
                        },
                    )
                    .collect::<Result<_, TransactionError>>()?;
                for ([first, second], ast) in bonds {
                    let id = self.add_bond(first, second, ast);
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
                let compaction = self.remove(&atoms, &bonds);
                state.compact(&compaction);
                Ok(())
            }
            Edit::ModifyAtomField { id, change } => {
                let id = state.atom(id)?;
                self.apply_modify_atom_field(id, change)
            }
            Edit::ModifyBondField { id, change } => {
                let id = state.bond(id)?;
                self.apply_modify_bond_field(id, change)
            }
            Edit::AddDativeBond { atoms, ast } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let (acceptor, donors) =
                    resolved
                        .split_last()
                        .ok_or(TransactionError::MalformedEdit(
                            "AddDativeBond requires at least one participant atom",
                        ))?;
                let id = self.add_dative_bond(donors.to_vec(), *acceptor, ast);
                state.push_dative_bond(id);
                Ok(())
            }
            Edit::RemoveDativeBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, atoms, ast) in removes {
                    let id = state.dative_bond(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    let (acceptor, donors) =
                        saved_atoms
                            .split_last()
                            .ok_or(TransactionError::MalformedEdit(
                                "RemoveDativeBond requires at least one participant atom",
                            ))?;
                    if !self.dative_bond_equiv(id, *acceptor, donors, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::DativeBond)?;
                let forward = IdCompaction::relations(
                    ids.iter().map(|&id| id.into()).collect(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                );
                self.remove_dative_bonds(&ids);
                state.compact(&forward);
                Ok(())
            }
            Edit::ModifyDativeBondField { id, change } => {
                let id = state.dative_bond(id)?;
                self.apply_modify_dative_bond_field(id, change)
            }
            Edit::AddAromaticSystem { atoms, ast } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_aromatic_system(resolved, ast);
                state.push_aromatic_system(id);
                Ok(())
            }
            Edit::RemoveAromaticSystems { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, atoms, ast) in removes {
                    let id = state.aromatic_system(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    if !self.aromatic_system_equiv(id, &saved_atoms, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::AromaticSystem)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    ids.iter().map(|&id| id.into()).collect(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                );
                self.remove_aromatic_systems(&ids);
                state.compact(&forward);
                Ok(())
            }
            Edit::ModifyAromaticSystemField { id, change } => {
                let id = state.aromatic_system(id)?;
                self.apply_modify_aromatic_system_field(id, change)
            }
            Edit::AddMulticenterBond { atoms, ast } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_multicenter_bond(resolved, ast);
                state.push_multicenter_bond(id);
                Ok(())
            }
            Edit::RemoveMulticenterBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, atoms, ast) in removes {
                    let id = state.multicenter_bond(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    if !self.multicenter_bond_equiv(id, &saved_atoms, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::MulticenterBond)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    Vec::new(),
                    ids.iter().map(|&id| id.into()).collect(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                );
                self.remove_multicenter_bonds(&ids);
                state.compact(&forward);
                Ok(())
            }
            Edit::ModifyMulticenterBondField { id, change } => {
                let id = state.multicenter_bond(id)?;
                self.apply_modify_multicenter_bond_field(id, change)
            }
            Edit::AddNoncovalentBond { atoms, ast } => {
                let a = state.atom(atoms[0].clone())?;
                let b = state.atom(atoms[1].clone())?;
                let id = self.add_noncovalent_bond([a, b], ast);
                state.push_noncovalent_bond(id);
                Ok(())
            }
            Edit::RemoveNoncovalentBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, atoms, ast) in removes {
                    let id = state.noncovalent_bond(id)?;
                    let saved_atoms =
                        [state.atom(atoms[0].clone())?, state.atom(atoms[1].clone())?];
                    if !self.noncovalent_bond_equiv(id, saved_atoms, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::NoncovalentBond)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    ids.iter().map(|&id| id.into()).collect(),
                    Vec::new(),
                    Vec::new(),
                );
                self.remove_noncovalent_bonds(&ids);
                state.compact(&forward);
                Ok(())
            }
            Edit::ModifyNoncovalentBondField { id, change } => {
                let id = state.noncovalent_bond(id)?;
                self.apply_modify_noncovalent_bond_field(id, change)
            }
            Edit::AddStereoAtom { site, ligands, ast } => {
                let site = state.atom(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_atom(site, ligands, ast);
                state.push_stereo_atom(id);
                Ok(())
            }
            Edit::RemoveStereoAtoms { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, site, ligands, ast) in removes {
                    let id = state.stereo_atom(id)?;
                    let site = state.atom(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    if !self.stereo_atom_equiv(id, site, &ligands, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::StereoAtom)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    ids.iter().map(|&id| id.into()).collect(),
                    Vec::new(),
                );
                self.remove_stereo_atoms(&ids);
                state.compact(&forward);
                Ok(())
            }
            Edit::ModifyStereoAtomField { id, change } => {
                let id = state.stereo_atom(id)?;
                self.apply_modify_stereo_atom_field(id, change)
            }
            Edit::AddStereoBond { site, ligands, ast } => {
                let site = state.bond(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_bond(site, ligands, ast);
                state.push_stereo_bond(id);
                Ok(())
            }
            Edit::RemoveStereoBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                for (id, site, ligands, ast) in removes {
                    let id = state.stereo_bond(id)?;
                    let site = state.bond(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    if !self.stereo_bond_equiv(id, site, &ligands, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::StereoBond)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    ids.iter().map(|&id| id.into()).collect(),
                );
                self.remove_stereo_bonds(&ids);
                state.compact(&forward);
                Ok(())
            }
            Edit::ModifyStereoBondField { id, change } => {
                let id = state.stereo_bond(id)?;
                self.apply_modify_stereo_bond_field(id, change)
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

    fn apply_edit(
        &mut self,
        edit: Edit,
        state: &mut ApplicationState,
    ) -> Result<Undo, TransactionError> {
        match edit {
            Edit::AddAtoms { atoms } => {
                let mut added = Vec::with_capacity(atoms.len());
                for ast in atoms {
                    let id = self.add_atom(ast.clone());
                    state.push_atom(id);
                    added.push(AddedAtom { id, ast });
                }
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
                             ast,
                         }| {
                            Ok(([state.atom(first)?, state.atom(second)?], ast))
                        },
                    )
                    .collect::<Result<_, TransactionError>>()?;
                let mut added = Vec::with_capacity(bonds.len());
                for ([first, second], ast) in bonds {
                    let id = self.add_bond(first, second, ast.clone());
                    state.push_bond(id);
                    added.push(AddedBond {
                        id,
                        endpoints: [first, second],
                        ast,
                    });
                }
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
                let pre_constraints = self.constraints().clone();
                let compaction = if !atoms.is_empty() || !bonds.is_empty() {
                    self.remove(&atoms, &bonds)
                } else {
                    IdCompaction::empty()
                };
                state.compact(&compaction);
                let mut constraints = pre_constraints;
                let cascade = constraints.compact_with_update(&compaction);
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
                let undo = Undo::ModifyAtomField {
                    id,
                    change: change.clone().inverse(),
                };
                self.apply_modify_atom_field(id, change)?;
                Ok(undo)
            }
            Edit::ModifyBondField { id, change } => {
                let id = state.bond(id)?;
                let undo = Undo::ModifyBondField {
                    id,
                    change: change.clone().inverse(),
                };
                self.apply_modify_bond_field(id, change)?;
                Ok(undo)
            }
            Edit::AddDativeBond { atoms, ast } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let (acceptor, donors) =
                    resolved
                        .split_last()
                        .ok_or(TransactionError::MalformedEdit(
                            "AddDativeBond requires at least one participant atom",
                        ))?;
                let id = self.add_dative_bond(donors.to_vec(), *acceptor, ast);
                state.push_dative_bond(id);
                let view = self.dative_bond(id);
                Ok(Undo::RemoveAddedDativeBond(AddedDativeBond {
                    id,
                    atoms: view.atom_ids().collect(),
                    ast: view.ast.clone(),
                }))
            }
            Edit::RemoveDativeBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, atoms, ast) in removes {
                    let id = state.dative_bond(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    let (acceptor, donors) =
                        saved_atoms
                            .split_last()
                            .ok_or(TransactionError::MalformedEdit(
                                "RemoveDativeBond requires at least one participant atom",
                            ))?;
                    if !self.dative_bond_equiv(id, *acceptor, donors, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.dative_bond(id);
                    let current_atoms: Vec<AtomId> = view.atom_ids().collect();
                    removed.push(RemovedDativeBond {
                        id,
                        atoms: current_atoms,
                        ast: view.ast.clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::DativeBond)?;
                let forward = IdCompaction::relations(
                    ids.iter().map(|&i| i.into()).collect(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                );
                let mut pre_constraints = self.constraints().clone();
                self.remove_dative_bonds(&ids);
                state.compact(&forward);
                let cascade = pre_constraints.compact_with_update(&forward);
                Ok(Undo::RestoreRemovedDativeBonds {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
            }
            Edit::ModifyDativeBondField { id, change } => {
                let id = state.dative_bond(id)?;
                let undo = Undo::ModifyDativeBondField {
                    id,
                    change: change.clone().inverse(),
                };
                self.apply_modify_dative_bond_field(id, change)?;
                Ok(undo)
            }
            Edit::AddAromaticSystem { atoms, ast } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_aromatic_system(resolved, ast);
                state.push_aromatic_system(id);
                let view = self.aromatic_system(id);
                Ok(Undo::RemoveAddedAromaticSystem(AddedAromaticSystem {
                    id,
                    atoms: view.atom_ids().collect(),
                    ast: view.ast.clone(),
                }))
            }
            Edit::RemoveAromaticSystems { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, atoms, ast) in removes {
                    let id = state.aromatic_system(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    if !self.aromatic_system_equiv(id, &saved_atoms, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.aromatic_system(id);
                    let current_atoms: Vec<AtomId> = view.atom_ids().collect();
                    removed.push(RemovedAromaticSystem {
                        id,
                        atoms: current_atoms,
                        ast: view.ast.clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::AromaticSystem)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    ids.iter().map(|&i| i.into()).collect(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                );
                let mut pre_constraints = self.constraints().clone();
                self.remove_aromatic_systems(&ids);
                state.compact(&forward);
                let cascade = pre_constraints.compact_with_update(&forward);
                Ok(Undo::RestoreRemovedAromaticSystems {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
            }
            Edit::ModifyAromaticSystemField { id, change } => {
                let id = state.aromatic_system(id)?;
                let undo = Undo::ModifyAromaticSystemField {
                    id,
                    change: change.clone().inverse(),
                };
                self.apply_modify_aromatic_system_field(id, change)?;
                Ok(undo)
            }
            Edit::AddMulticenterBond { atoms, ast } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| state.atom(r))
                    .collect::<Result<_, _>>()?;
                let id = self.add_multicenter_bond(resolved, ast);
                state.push_multicenter_bond(id);
                let view = self.multicenter_bond(id);
                Ok(Undo::RemoveAddedMulticenterBond(AddedMulticenterBond {
                    id,
                    atoms: view.atom_ids().collect(),
                    ast: view.ast.clone(),
                }))
            }
            Edit::RemoveMulticenterBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, atoms, ast) in removes {
                    let id = state.multicenter_bond(id)?;
                    let saved_atoms: Vec<AtomId> = atoms
                        .iter()
                        .map(|r| state.atom(r.clone()))
                        .collect::<Result<_, _>>()?;
                    if !self.multicenter_bond_equiv(id, &saved_atoms, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.multicenter_bond(id);
                    let current_atoms: Vec<AtomId> = view.atom_ids().collect();
                    removed.push(RemovedMulticenterBond {
                        id,
                        atoms: current_atoms,
                        ast: view.ast.clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::MulticenterBond)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    Vec::new(),
                    ids.iter().map(|&i| i.into()).collect(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                );
                let mut pre_constraints = self.constraints().clone();
                self.remove_multicenter_bonds(&ids);
                state.compact(&forward);
                let cascade = pre_constraints.compact_with_update(&forward);
                Ok(Undo::RestoreRemovedMulticenterBonds {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
            }
            Edit::ModifyMulticenterBondField { id, change } => {
                let id = state.multicenter_bond(id)?;
                let undo = Undo::ModifyMulticenterBondField {
                    id,
                    change: change.clone().inverse(),
                };
                self.apply_modify_multicenter_bond_field(id, change)?;
                Ok(undo)
            }
            Edit::AddNoncovalentBond { atoms, ast } => {
                let a = state.atom(atoms[0].clone())?;
                let b = state.atom(atoms[1].clone())?;
                let id = self.add_noncovalent_bond([a, b], ast);
                state.push_noncovalent_bond(id);
                let view = self.noncovalent_bond(id);
                Ok(Undo::RemoveAddedNoncovalentBond(AddedNoncovalentBond {
                    id,
                    atoms: view.atoms,
                    ast: view.ast.clone(),
                }))
            }
            Edit::RemoveNoncovalentBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, atoms, ast) in removes {
                    let id = state.noncovalent_bond(id)?;
                    let saved_atoms =
                        [state.atom(atoms[0].clone())?, state.atom(atoms[1].clone())?];
                    if !self.noncovalent_bond_equiv(id, saved_atoms, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.noncovalent_bond(id);
                    removed.push(RemovedNoncovalentBond {
                        id,
                        atoms: view.atoms,
                        ast: view.ast.clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::NoncovalentBond)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    ids.iter().map(|&i| i.into()).collect(),
                    Vec::new(),
                    Vec::new(),
                );
                let mut pre_constraints = self.constraints().clone();
                self.remove_noncovalent_bonds(&ids);
                state.compact(&forward);
                let cascade = pre_constraints.compact_with_update(&forward);
                Ok(Undo::RestoreRemovedNoncovalentBonds {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
            }
            Edit::ModifyNoncovalentBondField { id, change } => {
                let id = state.noncovalent_bond(id)?;
                let undo = Undo::ModifyNoncovalentBondField {
                    id,
                    change: change.clone().inverse(),
                };
                self.apply_modify_noncovalent_bond_field(id, change)?;
                Ok(undo)
            }
            Edit::AddStereoAtom { site, ligands, ast } => {
                let site = state.atom(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_atom(site, ligands.clone(), ast.clone());
                state.push_stereo_atom(id);
                Ok(Undo::RemoveAddedStereoAtom(AddedStereoAtom {
                    id,
                    site,
                    ligands,
                    ast,
                }))
            }
            Edit::RemoveStereoAtoms { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, site, ligands, ast) in removes {
                    let id = state.stereo_atom(id)?;
                    let site = state.atom(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    if !self.stereo_atom_equiv(id, site, &ligands, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.stereo_atom(id);
                    removed.push(RemovedStereoAtom {
                        id,
                        site: view.site,
                        ligands: view.ligands.to_vec(),
                        ast: view.ast.clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::StereoAtom)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    ids.iter().map(|&i| i.into()).collect(),
                    Vec::new(),
                );
                let mut pre_constraints = self.constraints().clone();
                self.remove_stereo_atoms(&ids);
                state.compact(&forward);
                let cascade = pre_constraints.compact_with_update(&forward);
                Ok(Undo::RestoreRemovedStereoAtoms {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
            }
            Edit::ModifyStereoAtomField { id, change } => {
                let id = state.stereo_atom(id)?;
                let undo = Undo::ModifyStereoAtomField {
                    id,
                    change: change.clone().inverse(),
                };
                self.apply_modify_stereo_atom_field(id, change)?;
                Ok(undo)
            }
            Edit::AddStereoBond { site, ligands, ast } => {
                let site = state.bond(site)?;
                let ligands = state.stereo_ligands(ligands)?;
                let id = self.add_stereo_bond(site, ligands.clone(), ast.clone());
                state.push_stereo_bond(id);
                Ok(Undo::RemoveAddedStereoBond(AddedStereoBond {
                    id,
                    site,
                    ligands,
                    ast,
                }))
            }
            Edit::RemoveStereoBonds { removes } => {
                let mut ids = Vec::with_capacity(removes.len());
                let mut removed = Vec::with_capacity(removes.len());
                for (id, site, ligands, ast) in removes {
                    let id = state.stereo_bond(id)?;
                    let site = state.bond(site)?;
                    let ligands = state.stereo_ligands(ligands)?;
                    if !self.stereo_bond_equiv(id, site, &ligands, &ast) {
                        return Err(TransactionError::OldStateMismatch);
                    }
                    let view = self.stereo_bond(id);
                    removed.push(RemovedStereoBond {
                        id,
                        site: view.site,
                        ligands: view.ligands.to_vec(),
                        ast: view.ast.clone(),
                    });
                    ids.push(id);
                }
                ensure_unique(&ids, EntityKind::StereoBond)?;
                let forward = IdCompaction::relations(
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    ids.iter().map(|&i| i.into()).collect(),
                );
                let mut pre_constraints = self.constraints().clone();
                self.remove_stereo_bonds(&ids);
                state.compact(&forward);
                let cascade = pre_constraints.compact_with_update(&forward);
                Ok(Undo::RestoreRemovedStereoBonds {
                    removed,
                    undo_compaction: forward.undo_compaction(),
                    cascade,
                })
            }
            Edit::ModifyStereoBondField { id, change } => {
                let id = state.stereo_bond(id)?;
                let undo = Undo::ModifyStereoBondField {
                    id,
                    change: change.clone().inverse(),
                };
                self.apply_modify_stereo_bond_field(id, change)?;
                Ok(undo)
            }
            Edit::ModifyAtomConstraint { id, old, new } => {
                let id = state.atom(id.clone())?;
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
                list.remove_at(position);
                Ok(Undo::ApplyCascadedConstraints(CascadedConstraints {
                    removed: vec![RemovedConstraint {
                        position,
                        constraint,
                    }],
                    modified: Vec::new(),
                }))
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
                ast: self.atom(id).ast.clone(),
            })
            .collect();
        let removed_bonds = (0..self.bond_count())
            .map(BondId::from)
            .filter(|&id| {
                let view = self.bond(id);
                bond_set.contains(&id) || view.atoms.iter().any(|atom| atom_set.contains(atom))
            })
            .map(|id| {
                let view = self.bond(id);
                RemovedBond {
                    id,
                    endpoints: view.atoms,
                    ast: view.ast.clone(),
                }
            })
            .collect();

        let dative_bonds = (0..self.dative_bond_count())
            .map(DativeBondId::from)
            .filter_map(|id| {
                let view = self.dative_bond(id);
                let atoms: Vec<AtomId> = view.atom_ids().collect();
                atoms
                    .iter()
                    .any(|a| atom_set.contains(a))
                    .then(|| RemovedDativeBond {
                        id,
                        atoms,
                        ast: view.ast.clone(),
                    })
            })
            .collect();
        let aromatic_systems = (0..self.aromatic_system_count())
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
                        ast: view.ast.clone(),
                    })
            })
            .collect();
        let multicenter_bonds = (0..self.multicenter_bond_count())
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
                        ast: view.ast.clone(),
                    })
            })
            .collect();
        let noncovalent_bonds = (0..self.noncovalent_bond_count())
            .map(NoncovalentBondId::from)
            .filter_map(|id| {
                let view = self.noncovalent_bond(id);
                view.atoms
                    .iter()
                    .any(|a| atom_set.contains(a))
                    .then(|| RemovedNoncovalentBond {
                        id,
                        atoms: view.atoms,
                        ast: view.ast.clone(),
                    })
            })
            .collect();

        // A stereo atom drops when its site atom or any ligand atom is removed;
        // a stereo bond drops when its site bond (directly or via a removed
        // endpoint) or any ligand atom is removed. Mirrors `birelation_removed`.
        let stereo_atoms = (0..self.stereo_atom_count())
            .map(StereoAtomId::from)
            .filter_map(|id| {
                let view = self.stereo_atom(id);
                let dropped = atom_set.contains(&view.site)
                    || view.ligands.iter().any(|l| atom_set.contains(&l.atom_id));
                dropped.then(|| RemovedStereoAtom {
                    id,
                    site: view.site,
                    ligands: view.ligands.to_vec(),
                    ast: view.ast.clone(),
                })
            })
            .collect();
        let stereo_bonds = (0..self.stereo_bond_count())
            .map(StereoBondId::from)
            .filter_map(|id| {
                let view = self.stereo_bond(id);
                let site = view.site;
                let site_dropped = bond_set.contains(&site)
                    || self.bond(site).atoms.iter().any(|a| atom_set.contains(a));
                let ligand_dropped = view.ligands.iter().any(|l| atom_set.contains(&l.atom_id));
                (site_dropped || ligand_dropped).then(|| RemovedStereoBond {
                    id,
                    site,
                    ligands: view.ligands.to_vec(),
                    ast: view.ast.clone(),
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

    fn apply_modify_atom_field(
        &mut self,
        id: AtomId,
        change: AtomFieldChange,
    ) -> Result<(), TransactionError> {
        let atom = self.atom_mut(id);
        match change {
            AtomFieldChange::Element { old, new } => {
                if !atom.ast.element.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.element = new;
            }
            AtomFieldChange::IsotopeMass { old, new } => {
                if !atom.ast.isotope_mass.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.isotope_mass = new;
            }
            AtomFieldChange::Charge { old, new } => {
                if !atom.ast.charge.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.charge = new;
            }
            AtomFieldChange::ImplicitHydrogens { old, new } => {
                if !atom.ast.implicit_hydrogens.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.implicit_hydrogens = new;
            }
            AtomFieldChange::LonePairs { old, new } => {
                if !atom.ast.lone_pairs.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.lone_pairs = new;
            }
            AtomFieldChange::UnpairedElectrons { old, new } => {
                if !atom.ast.unpaired_electrons.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.unpaired_electrons = new;
            }
        }
        Ok(())
    }

    fn apply_modify_bond_field(
        &mut self,
        id: BondId,
        change: BondFieldChange,
    ) -> Result<(), TransactionError> {
        let bond = self.bond_mut(id);
        match change {
            BondFieldChange::Order { old, new } => {
                if !bond.ast.order.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                bond.ast.order = new;
            }
            BondFieldChange::Charge { old, new } => {
                if !bond.ast.charge.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                bond.ast.charge = new;
            }
            BondFieldChange::UnpairedElectrons { old, new } => {
                if !bond.ast.unpaired_electrons.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                bond.ast.unpaired_electrons = new;
            }
        }
        Ok(())
    }

    fn apply_modify_dative_bond_field(
        &mut self,
        id: DativeBondId,
        change: DativeBondFieldChange,
    ) -> Result<(), TransactionError> {
        let dat = self.dative_bond_mut(id);
        match change {
            DativeBondFieldChange::Order { old, new } => {
                if !dat.ast.order.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                dat.ast.order = new;
            }
        }
        Ok(())
    }

    fn apply_modify_aromatic_system_field(
        &mut self,
        id: AromaticSystemId,
        change: AromaticSystemFieldChange,
    ) -> Result<(), TransactionError> {
        let ar = self.aromatic_system_mut(id);
        match change {
            AromaticSystemFieldChange::Electrons { old, new } => {
                if !ar.ast.electrons.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                ar.ast.electrons = new;
            }
            AromaticSystemFieldChange::Charge { old, new } => {
                if !ar.ast.charge.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                ar.ast.charge = new;
            }
            AromaticSystemFieldChange::UnpairedElectrons { old, new } => {
                if !ar.ast.unpaired_electrons.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                ar.ast.unpaired_electrons = new;
            }
        }
        Ok(())
    }

    fn apply_modify_multicenter_bond_field(
        &mut self,
        id: MulticenterBondId,
        change: MulticenterBondFieldChange,
    ) -> Result<(), TransactionError> {
        let mc = self.multicenter_bond_mut(id);
        match change {
            MulticenterBondFieldChange::Electrons { old, new } => {
                if !mc.ast.electrons.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                mc.ast.electrons = new;
            }
            MulticenterBondFieldChange::Charge { old, new } => {
                if !mc.ast.charge.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                mc.ast.charge = new;
            }
            MulticenterBondFieldChange::UnpairedElectrons { old, new } => {
                if !mc.ast.unpaired_electrons.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                mc.ast.unpaired_electrons = new;
            }
        }
        Ok(())
    }

    fn apply_modify_noncovalent_bond_field(
        &mut self,
        id: NoncovalentBondId,
        change: NoncovalentBondFieldChange,
    ) -> Result<(), TransactionError> {
        let nc = self.noncovalent_bond_mut(id);
        match change {
            NoncovalentBondFieldChange::Kind { old, new } => {
                if !nc.ast.kind.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                nc.ast.kind = new;
            }
        }
        Ok(())
    }

    fn apply_modify_stereo_atom_field(
        &mut self,
        id: StereoAtomId,
        change: StereoAtomFieldChange,
    ) -> Result<(), TransactionError> {
        let sa = self.stereo_atom_mut(id);
        match change {
            StereoAtomFieldChange::Configuration { old, new } => {
                if !sa.ast.configuration.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                sa.ast.configuration = new;
            }
        }
        Ok(())
    }

    fn apply_modify_stereo_bond_field(
        &mut self,
        id: StereoBondId,
        change: StereoBondFieldChange,
    ) -> Result<(), TransactionError> {
        let sb = self.stereo_bond_mut(id);
        match change {
            StereoBondFieldChange::Configuration { old, new } => {
                if !sb.ast.configuration.canonical_eq(&old) {
                    return Err(TransactionError::OldStateMismatch);
                }
                sb.ast.configuration = new;
            }
        }
        Ok(())
    }

    fn apply_modify_atom_constraint(
        &mut self,
        id: AtomId,
        old: Option<super::super::constraint::AtomConstraintForm>,
        new: Option<super::super::constraint::AtomConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.atom_mut(id)
            .ast
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_bond_constraint(
        &mut self,
        id: BondId,
        old: Option<super::super::constraint::BondConstraintForm>,
        new: Option<super::super::constraint::BondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.bond_mut(id)
            .ast
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_dative_bond_constraint(
        &mut self,
        id: DativeBondId,
        old: Option<super::super::constraint::DativeBondConstraintForm>,
        new: Option<super::super::constraint::DativeBondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.dative_bond_mut(id)
            .ast
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_aromatic_system_constraint(
        &mut self,
        id: AromaticSystemId,
        old: Option<super::super::constraint::AromaticSystemConstraintForm>,
        new: Option<super::super::constraint::AromaticSystemConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.aromatic_system_mut(id)
            .ast
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_multicenter_bond_constraint(
        &mut self,
        id: MulticenterBondId,
        old: Option<super::super::constraint::MulticenterBondConstraintForm>,
        new: Option<super::super::constraint::MulticenterBondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.multicenter_bond_mut(id)
            .ast
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_noncovalent_bond_constraint(
        &mut self,
        id: NoncovalentBondId,
        old: Option<super::super::constraint::NoncovalentBondConstraintForm>,
        new: Option<super::super::constraint::NoncovalentBondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.noncovalent_bond_mut(id)
            .ast
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_stereo_atom_constraint(
        &mut self,
        id: StereoAtomId,
        old: Option<super::super::constraint::StereoAtomConstraintForm>,
        new: Option<super::super::constraint::StereoAtomConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.stereo_atom_mut(id)
            .ast
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }

    fn apply_modify_stereo_bond_constraint(
        &mut self,
        id: StereoBondId,
        old: Option<super::super::constraint::StereoBondConstraintForm>,
        new: Option<super::super::constraint::StereoBondConstraintForm>,
    ) -> Result<(), TransactionError> {
        // A key mismatch (old/new different kinds) and an old-value mismatch both surface as
        // `compare_and_set`'s `Contradiction` → `OldStateMismatch`.
        self.stereo_bond_mut(id)
            .ast
            .constraints
            .compare_and_set(old, new)
            .map_err(|_| TransactionError::OldStateMismatch)
    }
}

fn rollback_journal(
    editor: &mut MoleculeEditor,
    journal: Vec<Undo>,
) -> Result<(), TransactionError> {
    for undo in journal.into_iter().rev() {
        editor.apply_undo(undo)?;
    }
    Ok(())
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
    let Some(restored_count) = current_count.checked_add(removed.len()) else {
        return false;
    };
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
    occupied.into_iter().all(|slot| slot)
}

fn restored_constraints(
    update: &CascadedConstraints,
    current: &Constraints,
) -> Option<Constraints> {
    let restored_count = current.as_slice().len().checked_add(update.removed.len())?;
    let mut removed = vec![None; restored_count];
    let mut modified = vec![None; restored_count];

    for entry in &update.removed {
        let slot = removed.get_mut(entry.position)?;
        if slot.is_some() {
            return None;
        }
        *slot = Some(&entry.constraint);
    }
    for entry in &update.modified {
        if removed.get(entry.position)?.is_some() {
            return None;
        }
        let slot = modified.get_mut(entry.position)?;
        if slot.is_some() {
            return None;
        }
        *slot = Some(entry);
    }

    let mut current = current.as_slice().iter();
    let mut restored = Vec::with_capacity(restored_count);
    for position in 0..restored_count {
        if let Some(constraint) = removed[position] {
            restored.push(constraint.clone());
            continue;
        }
        let constraint = current.next()?;
        if let Some(change) = modified[position] {
            if constraint != &change.new {
                return None;
            }
            restored.push(change.old.clone());
        } else {
            restored.push(constraint.clone());
        }
    }
    if current.next().is_some() {
        return None;
    }
    Some(restored.into_iter().collect())
}

fn removed_dative_bonds_fit(removed: &[RemovedDativeBond], atom_count: usize) -> bool {
    removed.iter().all(|entry| {
        !entry.atoms.is_empty() && entry.atoms.iter().all(|id| id.index() < atom_count)
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
    editor: &MoleculeEditor,
    removed: &RemovedOverlays,
    undo_compaction: &UndoCompaction,
    atom_count: usize,
    bond_count: usize,
) -> bool {
    reconstruction_fits(
        editor.dative_bond_count(),
        removed.dative_bonds.iter().map(|entry| entry.id.index()),
        |id| {
            undo_compaction
                .uncompact_dative_bond(DativeBondId::from(id))
                .index()
        },
    ) && reconstruction_fits(
        editor.aromatic_system_count(),
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
        editor.multicenter_bond_count(),
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
        editor.noncovalent_bond_count(),
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
        editor.stereo_atom_count(),
        removed.stereo_atoms.iter().map(|entry| entry.id.index()),
        |id| {
            undo_compaction
                .uncompact_stereo_atom(StereoAtomId::from(id))
                .index()
        },
    ) && reconstruction_fits(
        editor.stereo_bond_count(),
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

impl MoleculeEditor {
    fn validate_undo(&self, undo: &Undo) -> Result<(), TransactionError> {
        let fits = match undo {
            Undo::RemoveAddedTopology { atoms, bonds } => {
                ids_fit(
                    atoms.iter().map(|entry| entry.id.index()),
                    self.atom_count(),
                ) && ids_fit(
                    bonds.iter().map(|entry| entry.id.index()),
                    self.bond_count(),
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
                let Some(atom_count) = self.atom_count().checked_add(atoms.len()) else {
                    return Err(rollback_mismatch());
                };
                let Some(bond_count) = self.bond_count().checked_add(bonds.len()) else {
                    return Err(rollback_mismatch());
                };
                undo_compaction.forward() == compaction
                    && reconstruction_fits(
                        self.atom_count(),
                        atoms.iter().map(|entry| entry.id.index()),
                        |id| undo_compaction.uncompact_atom(AtomId::from(id)).index(),
                    )
                    && reconstruction_fits(
                        self.bond_count(),
                        bonds.iter().map(|entry| entry.id.index()),
                        |id| undo_compaction.uncompact_bond(BondId::from(id)).index(),
                    )
                    && bonds
                        .iter()
                        .all(|entry| entry.endpoints.iter().all(|id| id.index() < atom_count))
                    && removed_overlays_fit(self, overlays, undo_compaction, atom_count, bond_count)
            }
            Undo::RemoveAddedDativeBond(entry) => entry.id.index() < self.dative_bond_count(),
            Undo::RestoreRemovedDativeBonds {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.dative_bond_count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_dative_bond(DativeBondId::from(id))
                            .index()
                    },
                ) && removed_dative_bonds_fit(removed, self.atom_count())
            }
            Undo::RemoveAddedAromaticSystem(entry) => {
                entry.id.index() < self.aromatic_system_count()
            }
            Undo::RestoreRemovedAromaticSystems {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.aromatic_system_count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_aromatic_system(AromaticSystemId::from(id))
                            .index()
                    },
                ) && removed_aromatic_systems_fit(removed, self.atom_count())
            }
            Undo::RemoveAddedMulticenterBond(entry) => {
                entry.id.index() < self.multicenter_bond_count()
            }
            Undo::RestoreRemovedMulticenterBonds {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.multicenter_bond_count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_multicenter_bond(MulticenterBondId::from(id))
                            .index()
                    },
                ) && removed_multicenter_bonds_fit(removed, self.atom_count())
            }
            Undo::RemoveAddedNoncovalentBond(entry) => {
                entry.id.index() < self.noncovalent_bond_count()
            }
            Undo::RestoreRemovedNoncovalentBonds {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.noncovalent_bond_count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_noncovalent_bond(NoncovalentBondId::from(id))
                            .index()
                    },
                ) && removed_noncovalent_bonds_fit(removed, self.atom_count())
            }
            Undo::RemoveAddedStereoAtom(entry) => entry.id.index() < self.stereo_atom_count(),
            Undo::RestoreRemovedStereoAtoms {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.stereo_atom_count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_stereo_atom(StereoAtomId::from(id))
                            .index()
                    },
                ) && removed_stereo_atoms_fit(removed, self.atom_count())
            }
            Undo::RemoveAddedStereoBond(entry) => entry.id.index() < self.stereo_bond_count(),
            Undo::RestoreRemovedStereoBonds {
                removed,
                undo_compaction,
                cascade: _,
            } => {
                reconstruction_fits(
                    self.stereo_bond_count(),
                    removed.iter().map(|entry| entry.id.index()),
                    |id| {
                        undo_compaction
                            .uncompact_stereo_bond(StereoBondId::from(id))
                            .index()
                    },
                ) && removed_stereo_bonds_fit(removed, self.atom_count(), self.bond_count())
            }
            Undo::ModifyAtomField { id, .. } => id.index() < self.atom_count(),
            Undo::ModifyBondField { id, .. } => id.index() < self.bond_count(),
            Undo::ModifyDativeBondField { id, .. } => id.index() < self.dative_bond_count(),
            Undo::ModifyAromaticSystemField { id, .. } => id.index() < self.aromatic_system_count(),
            Undo::ModifyMulticenterBondField { id, .. } => {
                id.index() < self.multicenter_bond_count()
            }
            Undo::ModifyNoncovalentBondField { id, .. } => {
                id.index() < self.noncovalent_bond_count()
            }
            Undo::ModifyStereoAtomField { id, .. } => id.index() < self.stereo_atom_count(),
            Undo::ModifyStereoBondField { id, .. } => id.index() < self.stereo_bond_count(),
            Undo::ApplyCascadedConstraints(_) => true,
            Undo::ApplyEdit(_) => true,
        };
        fits.then_some(()).ok_or_else(rollback_mismatch)
    }

    fn apply_undo(&mut self, undo: Undo) -> Result<(), TransactionError> {
        self.validate_undo(&undo)?;
        match undo {
            Undo::RemoveAddedTopology { atoms, bonds } => {
                self.remove_added_topology(&atoms, &bonds);
            }
            Undo::RestoreRemovedTopology {
                atoms,
                bonds,
                overlays,
                undo_compaction,
                cascade,
                ..
            } => {
                let constraints = restored_constraints(&cascade, self.constraints())
                    .ok_or_else(rollback_mismatch)?;
                self.restore_topology(atoms, bonds, overlays, &undo_compaction);
                *self.constraints_mut() = constraints;
            }
            Undo::RemoveAddedDativeBond(added) => self.remove_added_dative_bond(&added),
            Undo::RestoreRemovedDativeBonds {
                removed,
                undo_compaction,
                cascade,
            } => {
                let constraints = restored_constraints(&cascade, self.constraints())
                    .ok_or_else(rollback_mismatch)?;
                self.restore_dative_bonds(removed, &undo_compaction);
                *self.constraints_mut() = constraints;
            }
            Undo::RemoveAddedAromaticSystem(added) => self.remove_added_aromatic_system(&added),
            Undo::RestoreRemovedAromaticSystems {
                removed,
                undo_compaction,
                cascade,
            } => {
                let constraints = restored_constraints(&cascade, self.constraints())
                    .ok_or_else(rollback_mismatch)?;
                self.restore_aromatic_systems(removed, &undo_compaction);
                *self.constraints_mut() = constraints;
            }
            Undo::RemoveAddedMulticenterBond(added) => self.remove_added_multicenter_bond(&added),
            Undo::RestoreRemovedMulticenterBonds {
                removed,
                undo_compaction,
                cascade,
            } => {
                let constraints = restored_constraints(&cascade, self.constraints())
                    .ok_or_else(rollback_mismatch)?;
                self.restore_multicenter_bonds(removed, &undo_compaction);
                *self.constraints_mut() = constraints;
            }
            Undo::RemoveAddedNoncovalentBond(added) => self.remove_added_noncovalent_bond(&added),
            Undo::RestoreRemovedNoncovalentBonds {
                removed,
                undo_compaction,
                cascade,
            } => {
                let constraints = restored_constraints(&cascade, self.constraints())
                    .ok_or_else(rollback_mismatch)?;
                self.restore_noncovalent_bonds(removed, &undo_compaction);
                *self.constraints_mut() = constraints;
            }
            Undo::RemoveAddedStereoAtom(added) => self.remove_added_stereo_atom(&added),
            Undo::RestoreRemovedStereoAtoms {
                removed,
                undo_compaction,
                cascade,
            } => {
                let constraints = restored_constraints(&cascade, self.constraints())
                    .ok_or_else(rollback_mismatch)?;
                self.restore_stereo_atoms(removed, &undo_compaction);
                *self.constraints_mut() = constraints;
            }
            Undo::RemoveAddedStereoBond(added) => self.remove_added_stereo_bond(&added),
            Undo::RestoreRemovedStereoBonds {
                removed,
                undo_compaction,
                cascade,
            } => {
                let constraints = restored_constraints(&cascade, self.constraints())
                    .ok_or_else(rollback_mismatch)?;
                self.restore_stereo_bonds(removed, &undo_compaction);
                *self.constraints_mut() = constraints;
            }
            Undo::ModifyAtomField { id, change } => self
                .apply_modify_atom_field(id, change)
                .map_err(|_| rollback_mismatch())?,
            Undo::ModifyBondField { id, change } => self
                .apply_modify_bond_field(id, change)
                .map_err(|_| rollback_mismatch())?,
            Undo::ModifyDativeBondField { id, change } => self
                .apply_modify_dative_bond_field(id, change)
                .map_err(|_| rollback_mismatch())?,
            Undo::ModifyAromaticSystemField { id, change } => self
                .apply_modify_aromatic_system_field(id, change)
                .map_err(|_| rollback_mismatch())?,
            Undo::ModifyMulticenterBondField { id, change } => self
                .apply_modify_multicenter_bond_field(id, change)
                .map_err(|_| rollback_mismatch())?,
            Undo::ModifyNoncovalentBondField { id, change } => self
                .apply_modify_noncovalent_bond_field(id, change)
                .map_err(|_| rollback_mismatch())?,
            Undo::ModifyStereoAtomField { id, change } => self
                .apply_modify_stereo_atom_field(id, change)
                .map_err(|_| rollback_mismatch())?,
            Undo::ModifyStereoBondField { id, change } => self
                .apply_modify_stereo_bond_field(id, change)
                .map_err(|_| rollback_mismatch())?,
            Undo::ApplyCascadedConstraints(update) => {
                let constraints = restored_constraints(&update, self.constraints())
                    .ok_or_else(rollback_mismatch)?;
                *self.constraints_mut() = constraints;
            }
            Undo::ApplyEdit(edit) => {
                let mut state = ApplicationState::new(self);
                self.apply_edit(*edit, &mut state)
                    .map_err(|_| rollback_mismatch())?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use rstest::*;
    use umol_chem::element::Element;

    use super::super::super::aromatic::AromaticSystemForm;
    use super::super::super::atom::{AtomForm, ElementForm};
    use super::super::super::bond::BondForm;
    use super::super::super::constraint::{
        AromaticSystemConstraintForm, AtomConstraintForm, BondConstraintForm, Constraint,
        DativeBondConstraintForm, MoleculeConstraint, MulticenterBondConstraintForm,
        NoncovalentBondConstraintForm, RelationalConstraint, RingScope, StereoAtomConstraintForm,
        StereoBondConstraintForm, StereogenicityForm, SubPatternAnchor,
    };
    use super::super::super::dative::DativeBondForm;
    use super::super::super::edit::EntityHandle;
    use super::super::super::entity::Entity;
    use super::super::super::ligand::StereoLigandKind;
    use super::super::super::multicenter::MulticenterBondForm;
    use super::super::super::noncovalent::{
        NoncovalentBondForm, NoncovalentBondKind, NoncovalentBondKindForm,
    };
    use super::super::super::stereo::{
        CisTransStereoForm, StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoCoset,
        StereoKind,
    };
    use super::super::super::value::NumForm;
    use super::super::{MoleculeAst, MoleculeEntries};
    use super::*;
    use crate::ir::BooleanForm;

    #[fixture]
    fn empty() -> MoleculeEditor {
        MoleculeAst::default().edit()
    }

    #[fixture]
    fn one_atom() -> MoleculeEditor {
        let mut b = MoleculeAst::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b
    }

    #[fixture]
    fn diatomic() -> MoleculeEditor {
        let mut b = MoleculeAst::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b
    }

    #[rstest]
    fn test_molecule_editor_transact_add_atom(mut empty: MoleculeEditor) {
        let mut edits = Edits::new();
        edits.add_atom(AtomForm::from_element(Element::C));
        let tx = empty.transact(edits).unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedTopology { atoms, bonds }]
                if atoms.iter().map(|a| a.id).collect::<Vec<_>>() == vec![AtomId(0)]
                    && bonds.is_empty()
        ));
        let built = empty.build();
        assert_eq!(built.atoms().count(), 1);
        assert_eq!(
            built.atom(AtomId(0)).ast.element,
            ElementForm::Lit(Element::C)
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_add_atoms(mut empty: MoleculeEditor) {
        let mut edits = Edits::new();
        edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
        ]);
        let tx = empty.transact(edits).unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedTopology { atoms, bonds }]
                if atoms.iter().map(|a| a.id).collect::<Vec<_>>() == vec![AtomId(0), AtomId(1)]
                    && bonds.is_empty()
        ));
        let built = empty.build();
        assert_eq!(built.atoms().count(), 2);
        assert_eq!(
            built.atom(AtomId(0)).ast.element,
            ElementForm::Lit(Element::C)
        );
        assert_eq!(
            built.atom(AtomId(1)).ast.element,
            ElementForm::Lit(Element::N)
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_add_bond_via_handle(mut empty: MoleculeEditor) {
        let mut edits = Edits::new();
        let atoms = edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ]);
        edits.add_bond(atoms[0].clone(), atoms[1].clone(), BondForm::from_order(1));
        let tx = empty.transact(edits).unwrap();
        assert!(matches!(
            tx.undos(),
            [
                Undo::RemoveAddedTopology { atoms, bonds },
                Undo::RemoveAddedTopology { atoms: bond_atoms, bonds: added_bonds },
            ] if atoms.iter().map(|a| a.id).collect::<Vec<_>>() == vec![AtomId(0), AtomId(1)]
                && bonds.is_empty()
                && bond_atoms.is_empty()
                && added_bonds.iter().map(|b| b.id).collect::<Vec<_>>() == vec![BondId(0)]
        ));
    }

    #[rstest]
    fn test_molecule_editor_transact_rollback(mut one_atom: MoleculeEditor) {
        let before = one_atom.clone().build();
        // Mid-batch failure (out-of-range id on edit 2) rolls back the
        // already-applied AddAtom on edit 1.
        let mut edits = Edits::new();
        edits.add_atom(AtomForm::from_element(Element::N));
        edits.remove_atom(AtomHandle::Id(AtomId(99)));
        let err = one_atom.transact(edits).unwrap_err();
        assert_eq!(
            err,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 99,
                count: 1,
            }
        );
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_atom_field(mut one_atom: MoleculeEditor) {
        let tx = one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap();
        assert_eq!(
            tx.undos(),
            &[Undo::ModifyAtomField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::default(),
                },
            }],
        );
        assert_eq!(one_atom.build().atom(AtomId(0)).ast.charge, NumForm::Lit(1));
    }

    #[rstest]
    fn test_molecule_editor_transact_set_atom_field_error(mut one_atom: MoleculeEditor) {
        let err = one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(99),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    #[case::created_out_of_range(
        Edits::from_iter([Edit::AddBonds {
            bonds: vec![AddBond {
                endpoints: [AtomHandle::New(5), AtomHandle::New(6)],
                ast: BondForm::default(),
            }],
        }]),
        TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 5,
            count: 0,
        },
    )]
    #[case::initial_out_of_range(
        Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        }]),
        TransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 0,
            count: 0,
        },
    )]
    fn test_molecule_editor_transact_handle_error(
        mut empty: MoleculeEditor,
        #[case] edits: Edits,
        #[case] expected: TransactionError,
    ) {
        let before = empty.clone().build();
        assert_eq!(empty.transact(edits).unwrap_err(), expected);
        assert_eq!(empty.build(), before);
    }

    #[rstest]
    #[case::initial(
        1,
        Edits::from_iter([
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            },
        ]),
    )]
    #[case::created(
        0,
        Edits::from_iter([
            Edit::AddAtoms {
                atoms: vec![AtomForm::from_element(Element::C)],
            },
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::New(0)],
                bonds: Vec::new(),
            },
            Edit::ModifyAtomField {
                id: AtomHandle::New(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            },
        ]),
    )]
    fn test_molecule_editor_transact_handle_removed_error(
        #[case] initial_atom_count: usize,
        #[case] edits: Edits,
    ) {
        let mut editor = MoleculeAst::default().edit();
        for _ in 0..initial_atom_count {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        let before = editor.clone().build();

        assert_eq!(
            editor.transact(edits).unwrap_err(),
            TransactionError::HandleRemoved {
                kind: EntityKind::Atom,
                index: 0,
            }
        );
        assert_eq!(editor.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_handles_initial() {
        let mut editor = MoleculeAst::default().edit();
        editor.add_atom(AtomForm::from_element(Element::C));
        editor.add_atom(AtomForm::from_element(Element::N));
        editor.add_atom(AtomForm::from_element(Element::O));
        let edits = Edits::from_iter([
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(1)),
                change: AtomFieldChange::Element {
                    old: ElementForm::Lit(Element::N),
                    new: ElementForm::Lit(Element::F),
                },
            },
            Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(2)),
                change: AtomFieldChange::Element {
                    old: ElementForm::Lit(Element::O),
                    new: ElementForm::Lit(Element::Cl),
                },
            },
        ]);

        editor.transact(edits).unwrap();

        assert_eq!(
            (0..editor.atom_count())
                .map(|index| editor.atom(AtomId(index as u32)).ast.element.clone())
                .collect::<Vec<_>>(),
            vec![ElementForm::Lit(Element::F), ElementForm::Lit(Element::Cl)]
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_handles_created() {
        let mut editor = MoleculeAst::default().edit();
        let mut edits = Edits::new();
        let atoms = edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
        ]);
        edits.remove_atom(atoms[0].clone());
        edits.push(Edit::ModifyAtomField {
            id: atoms[1].clone(),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        });

        editor.transact(edits).unwrap();

        assert_eq!(editor.atom_count(), 1);
        assert_eq!(
            (
                editor.atom(AtomId(0)).ast.element.clone(),
                editor.atom(AtomId(0)).ast.charge.clone(),
            ),
            (ElementForm::Lit(Element::N), NumForm::Lit(1))
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_handles_reuse() {
        let mut editor = MoleculeAst::default().edit();
        let mut edits = Edits::new();
        let removed = edits.add_atom(AtomForm::from_element(Element::C));
        edits.remove_atom(removed);
        let surviving = edits.add_atom(AtomForm::from_element(Element::N));
        edits.push(Edit::ModifyAtomField {
            id: surviving,
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(-1),
            },
        });

        editor.transact(edits).unwrap();

        assert_eq!(editor.atom_count(), 1);
        assert_eq!(
            (
                editor.atom(AtomId(0)).ast.element.clone(),
                editor.atom(AtomId(0)).ast.charge.clone(),
            ),
            (ElementForm::Lit(Element::N), NumForm::Lit(-1))
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_handles_per_kind() {
        let mut editor = MoleculeAst::default().edit();
        let before = editor.clone().build();
        let mut edits = Edits::new();
        let atoms = edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::F),
        ]);
        let bonds = edits.add_bonds([
            AddBond {
                endpoints: [atoms[0].clone(), atoms[1].clone()],
                ast: BondForm::from_order(1),
            },
            AddBond {
                endpoints: [atoms[1].clone(), atoms[2].clone()],
                ast: BondForm::from_order(1),
            },
            AddBond {
                endpoints: [atoms[2].clone(), atoms[3].clone()],
                ast: BondForm::from_order(1),
            },
        ]);
        let dative = edits.add_dative_bond(
            vec![atoms[0].clone(), atoms[1].clone()],
            DativeBondForm::from_order(1),
        );
        let aromatic = edits.add_aromatic_system(
            vec![atoms[0].clone(), atoms[1].clone()],
            AromaticSystemForm::default(),
        );
        let multicenter = edits.add_multicenter_bond(
            vec![atoms[0].clone(), atoms[1].clone()],
            MulticenterBondForm::default(),
        );
        let noncovalent = edits.add_noncovalent_bond(
            [atoms[0].clone(), atoms[1].clone()],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        let stereo_atom = edits.add_stereo_atom(
            atoms[1].clone(),
            vec![
                (atoms[0].clone(), StereoLigandKind::Atom),
                (atoms[2].clone(), StereoLigandKind::Atom),
                (atoms[3].clone(), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let stereo_bond = edits.add_stereo_bond(
            bonds[1].clone(),
            vec![
                (atoms[0].clone(), StereoLigandKind::Atom),
                (atoms[3].clone(), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        edits.push(Edit::ModifyAtomField {
            id: atoms[0].clone(),
            change: AtomFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        });
        edits.push(Edit::ModifyBondField {
            id: bonds[0].clone(),
            change: BondFieldChange::Order {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        });
        edits.push(Edit::ModifyDativeBondField {
            id: dative,
            change: DativeBondFieldChange::Order {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        });
        edits.push(Edit::ModifyAromaticSystemField {
            id: aromatic,
            change: AromaticSystemFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(1),
            },
        });
        edits.push(Edit::ModifyMulticenterBondField {
            id: multicenter,
            change: MulticenterBondFieldChange::Charge {
                old: NumForm::default(),
                new: NumForm::Lit(-1),
            },
        });
        edits.push(Edit::ModifyNoncovalentBondField {
            id: noncovalent,
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
            },
        });
        edits.push(Edit::ModifyStereoAtomField {
            id: stereo_atom,
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            },
        });
        edits.push(Edit::ModifyStereoBondField {
            id: stereo_bond,
            change: StereoBondFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
            },
        });

        let transaction = editor.transact(edits).unwrap();

        assert_eq!(editor.atom(AtomId(0)).ast.charge, NumForm::Lit(1));
        assert_eq!(editor.bond(BondId(0)).ast.order, NumForm::Lit(2));
        assert_eq!(
            editor.dative_bond(DativeBondId(0)).ast.order,
            NumForm::Lit(2)
        );
        assert_eq!(
            editor.aromatic_system(AromaticSystemId(0)).ast.charge,
            NumForm::Lit(1)
        );
        assert_eq!(
            editor.multicenter_bond(MulticenterBondId(0)).ast.charge,
            NumForm::Lit(-1)
        );
        assert_eq!(
            editor.noncovalent_bond(NoncovalentBondId(0)).ast.kind,
            NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic)
        );
        assert_eq!(
            editor.stereo_atom(StereoAtomId(0)).ast.configuration,
            StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0))
        );
        assert_eq!(
            editor.stereo_bond(StereoBondId(0)).ast.configuration,
            StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0))
        );

        transaction.rollback(&mut editor).unwrap();
        assert_eq!(editor.build(), before);
    }

    #[rstest]
    #[case::first(0)]
    #[case::middle(1)]
    #[case::last(2)]
    fn test_molecule_editor_transact_add_bonds_error(
        mut diatomic: MoleculeEditor,
        #[case] invalid_position: usize,
    ) {
        let before = diatomic.clone().build();
        let mut bonds = vec![
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: BondForm::from_order(1),
            },
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: BondForm::from_order(2),
            },
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: BondForm::from_order(3),
            },
        ];
        bonds[invalid_position].endpoints[1] = AtomHandle::Id(AtomId(9));

        assert_eq!(
            diatomic
                .transact(Edits::from_iter([Edit::AddBonds { bonds }]))
                .unwrap_err(),
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 9,
                count: 2,
            }
        );
        assert_eq!(diatomic.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_topology(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::RemoveTopology {
                atoms: Vec::new(),
                bonds: vec![BondHandle::Id(BondId(0))],
            }]))
            .unwrap();

        assert_eq!(diatomic.bond_count(), 0);
        let [Undo::RestoreRemovedTopology {
            atoms,
            bonds,
            overlays,
            compaction,
            ..
        }] = tx.undos()
        else {
            panic!("RemoveTopology should produce one topology-restore undo")
        };
        assert!(atoms.is_empty());
        assert_eq!(
            bonds.iter().map(|b| b.id).collect::<Vec<_>>(),
            vec![BondId(0)]
        );
        assert!(overlays.dative_bonds.is_empty());
        assert_eq!(compaction.compact_bond(BondId(0)), None);
    }

    #[rstest]
    fn test_molecule_editor_transact_add_atom_constraint(mut one_atom: MoleculeEditor) {
        one_atom
            .transact(Edits::from_iter([
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old: None,
                    new: Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 1)),
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old: None,
                    new: Some(AtomConstraintForm::ring_membership(RingScope::Size(6), 1)),
                },
            ]))
            .unwrap();
        let next = one_atom.build();
        let cs: Vec<_> = next
            .atom(AtomId(0))
            .ast
            .constraints
            .iter()
            .cloned()
            .collect();
        assert_eq!(
            cs,
            vec![
                AtomConstraintForm::ring_membership(RingScope::Size(5), 1),
                AtomConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]
        );
    }

    #[rstest]
    #[case::singleton_set(NumForm::Lit(1), NumForm::lit_set([1]))]
    fn test_molecule_editor_transact_modify_atom_field_canonical(
        mut one_atom: MoleculeEditor,
        #[case] current: NumForm,
        #[case] old: NumForm,
    ) {
        // The modify's recorded `old` is canonically equal to — but structurally distinct from — the
        // stored charge, so the old-state check passes (structural `!=` would raise `OldStateMismatch`).
        one_atom.atom_mut(AtomId(0)).ast.charge = current;
        one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old,
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap();
        assert_eq!(one_atom.atom_mut(AtomId(0)).ast.charge, NumForm::Lit(2));
    }

    #[rstest]
    fn test_molecule_editor_transact_modify_atom_constraint_absent_error(
        mut one_atom: MoleculeEditor,
    ) {
        let err = one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old: Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 1)),
                new: None,
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    #[case::introduce(None, Some(AtomConstraintForm::valence(4)), Some(NumForm::Lit(4)))]
    #[case::replace(
        Some(AtomConstraintForm::valence(3)),
        Some(AtomConstraintForm::valence(4)),
        Some(NumForm::Lit(4))
    )]
    #[case::remove(Some(AtomConstraintForm::valence(3)), None, None)]
    fn test_molecule_editor_transact_set_atom_constraint(
        mut one_atom: MoleculeEditor,
        #[case] old: Option<AtomConstraintForm>,
        #[case] new: Option<AtomConstraintForm>,
        #[case] expected: Option<NumForm>,
    ) {
        if let Some(c) = old.clone() {
            one_atom.atom_mut(AtomId(0)).ast.constraints.set(c);
        }
        one_atom
            .transact(Edits::from_iter([Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old,
                new,
            }]))
            .unwrap();
        assert_eq!(
            one_atom.atom_mut(AtomId(0)).ast.constraints.valence(),
            expected.as_ref()
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_bond_constraint(mut diatomic: MoleculeEditor) {
        diatomic
            .transact(Edits::from_iter([Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(0)),
                old: None,
                new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            }]))
            .unwrap();
        assert!(diatomic
            .bond_mut(BondId(0))
            .ast
            .constraints
            .iter()
            .any(|c| *c == BondConstraintForm::Aromatic(BooleanForm::Lit(true))));
    }

    #[rstest]
    fn test_molecule_editor_transact_add_molecule_constraint(mut empty: MoleculeEditor) {
        let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        empty
            .transact(Edits::from_iter([Edit::AddMoleculeConstraint {
                constraint: c.clone().into(),
            }]))
            .unwrap();
        assert_eq!(empty.constraints_mut().as_slice(), &[c]);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_molecule_constraint(mut empty: MoleculeEditor) {
        let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        empty.push_constraint(c.clone());
        empty
            .transact(Edits::from_iter([Edit::RemoveMoleculeConstraint {
                constraint: c.clone().into(),
            }]))
            .unwrap();
        assert!(empty.constraints_mut().as_slice().is_empty());
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_molecule_constraint_absent_error(
        mut empty: MoleculeEditor,
    ) {
        let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        empty.push_constraint(c.clone());
        let err = empty
            .transact(Edits::from_iter([Edit::RemoveMoleculeConstraint {
                constraint: Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: None,
                    sum: NumForm::Lit(0),
                })
                .into(),
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::MissingEntry);
        assert_eq!(empty.constraints_mut().as_slice(), &[c]);
    }

    #[rstest]
    fn test_molecule_editor_transact_molecule_constraint_initial(
        mut batched_overlays: MoleculeEditor,
    ) {
        let constraint = Constraint::And(vec![
            Constraint::Atom(AtomId(1), AtomConstraintForm::valence(3_i64)),
            Constraint::Bond(BondId(1), BondConstraintForm::aromatic(true)),
            Constraint::DativeBond(DativeBondId(1), DativeBondConstraintForm::aromatic(true)),
            Constraint::AromaticSystem(
                AromaticSystemId(1),
                AromaticSystemConstraintForm::electron_count(6_i64),
            ),
            Constraint::MulticenterBond(
                MulticenterBondId(1),
                MulticenterBondConstraintForm::electron_count(2_i64),
            ),
            Constraint::NoncovalentBond(
                NoncovalentBondId(1),
                NoncovalentBondConstraintForm::intramolecular(true),
            ),
            Constraint::StereoAtom(
                StereoAtomId(1),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::StereoBond(
                StereoBondId(1),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
        ]);
        let before = batched_overlays.clone().build();
        let mut edits = Edits::new();
        edits.add_molecule_constraint(constraint.clone().into());

        let transaction = batched_overlays.transact(edits).unwrap();
        assert_eq!(batched_overlays.constraints().as_slice(), &[constraint]);

        transaction.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_molecule_constraint_created(mut empty: MoleculeEditor) {
        let before = empty.clone().build();
        let mut edits = Edits::new();
        let atoms = edits.add_atoms([
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
        ]);
        let bond = edits.add_bond(atoms[0].clone(), atoms[1].clone(), BondForm::from_order(1));
        let dative = edits.add_dative_bond(
            vec![atoms[0].clone(), atoms[1].clone()],
            DativeBondForm::from_order(1),
        );
        let aromatic = edits.add_aromatic_system(atoms.clone(), AromaticSystemForm::default());
        let multicenter = edits.add_multicenter_bond(atoms.clone(), MulticenterBondForm::default());
        let noncovalent = edits.add_noncovalent_bond(
            [atoms[0].clone(), atoms[1].clone()],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        let stereo_atom = edits.add_stereo_atom(
            atoms[0].clone(),
            vec![(atoms[1].clone(), StereoLigandKind::Atom)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let stereo_bond = edits.add_stereo_bond(
            bond.clone(),
            vec![
                (atoms[0].clone(), StereoLigandKind::Atom),
                (atoms[1].clone(), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let source = Constraint::And(vec![
            Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
            Constraint::Bond(BondId(7), BondConstraintForm::aromatic(true)),
            Constraint::DativeBond(DativeBondId(7), DativeBondConstraintForm::aromatic(true)),
            Constraint::AromaticSystem(
                AromaticSystemId(7),
                AromaticSystemConstraintForm::electron_count(6_i64),
            ),
            Constraint::MulticenterBond(
                MulticenterBondId(7),
                MulticenterBondConstraintForm::electron_count(2_i64),
            ),
            Constraint::NoncovalentBond(
                NoncovalentBondId(7),
                NoncovalentBondConstraintForm::intramolecular(true),
            ),
            Constraint::StereoAtom(
                StereoAtomId(7),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::StereoBond(
                StereoBondId(7),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::Relational(RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(7),
                parallel: BondId(7),
            }),
        ]);
        let mappings = HashMap::from([
            (
                Entity::Atom(AtomId(7)),
                EntityHandle::Atom(atoms[0].clone()),
            ),
            (Entity::Bond(BondId(7)), EntityHandle::Bond(bond)),
            (
                Entity::DativeBond(DativeBondId(7)),
                EntityHandle::DativeBond(dative),
            ),
            (
                Entity::AromaticSystem(AromaticSystemId(7)),
                EntityHandle::AromaticSystem(aromatic),
            ),
            (
                Entity::MulticenterBond(MulticenterBondId(7)),
                EntityHandle::MulticenterBond(multicenter),
            ),
            (
                Entity::NoncovalentBond(NoncovalentBondId(7)),
                EntityHandle::NoncovalentBond(noncovalent),
            ),
            (
                Entity::StereoAtom(StereoAtomId(7)),
                EntityHandle::StereoAtom(stereo_atom),
            ),
            (
                Entity::StereoBond(StereoBondId(7)),
                EntityHandle::StereoBond(stereo_bond),
            ),
        ]);
        edits.add_molecule_constraint(
            ConstraintEdit::new(source, |entity| mappings.get(&entity).cloned()).unwrap(),
        );
        let expected = Constraint::And(vec![
            Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3_i64)),
            Constraint::Bond(BondId(0), BondConstraintForm::aromatic(true)),
            Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(true)),
            Constraint::AromaticSystem(
                AromaticSystemId(0),
                AromaticSystemConstraintForm::electron_count(6_i64),
            ),
            Constraint::MulticenterBond(
                MulticenterBondId(0),
                MulticenterBondConstraintForm::electron_count(2_i64),
            ),
            Constraint::NoncovalentBond(
                NoncovalentBondId(0),
                NoncovalentBondConstraintForm::intramolecular(true),
            ),
            Constraint::StereoAtom(
                StereoAtomId(0),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::StereoBond(
                StereoBondId(0),
                StereoKind::CisTrans,
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
            ),
            Constraint::Relational(RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(0),
                parallel: BondId(0),
            }),
        ]);

        let transaction = empty.transact(edits).unwrap();
        assert_eq!(empty.constraints().as_slice(), &[expected]);

        transaction.rollback(&mut empty).unwrap();
        assert_eq!(empty.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_molecule_constraint_compaction(
        mut batched_overlays: MoleculeEditor,
    ) {
        let removed = Constraint::AromaticSystem(
            AromaticSystemId(1),
            AromaticSystemConstraintForm::electron_count(6_i64),
        );
        let added = Constraint::AromaticSystem(
            AromaticSystemId(1),
            AromaticSystemConstraintForm::electron_count(4_i64),
        );
        batched_overlays.push_constraint(removed.clone());
        let before = batched_overlays.clone().build();
        let mut edits = Edits::from_iter([Edit::RemoveAromaticSystems {
            removes: vec![(
                AromaticSystemHandle::Id(AromaticSystemId(0)),
                vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                AromaticSystemForm::default(),
            )],
        }]);
        edits.add_molecule_constraint(added.into());
        edits.remove_molecule_constraint(removed.into());

        let transaction = batched_overlays.transact(edits).unwrap();
        assert_eq!(
            batched_overlays.constraints().as_slice(),
            &[Constraint::AromaticSystem(
                AromaticSystemId(0),
                AromaticSystemConstraintForm::electron_count(4_i64),
            )],
        );

        transaction.rollback(&mut batched_overlays).unwrap();
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    #[case::forward(
        0,
        Edits::from_iter([Edit::AddMoleculeConstraint {
            constraint: ConstraintEdit::new(
                Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
                |_| Some(EntityHandle::Atom(AtomHandle::New(0))),
            ).unwrap(),
        }]),
        TransactionError::HandleOutOfRange { kind: EntityKind::Atom, index: 0, count: 0 },
    )]
    #[case::removed(
        1,
        Edits::from_iter([
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            Edit::AddMoleculeConstraint {
                constraint: ConstraintEdit::new(
                    Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
                    |_| Some(EntityHandle::Atom(AtomHandle::Id(AtomId(0)))),
                ).unwrap(),
            },
        ]),
        TransactionError::HandleRemoved { kind: EntityKind::Atom, index: 0 },
    )]
    fn test_molecule_editor_transact_molecule_constraint_error(
        #[case] initial_atom_count: usize,
        #[case] edits: Edits,
        #[case] expected: TransactionError,
    ) {
        let mut editor = MoleculeAst::default().edit();
        for _ in 0..initial_atom_count {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        let before = editor.clone().build();

        assert_eq!(editor.transact(edits), Err(expected));
        assert_eq!(editor.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_molecule_constraint_subpattern(mut one_atom: MoleculeEditor) {
        let pattern = MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 4],
            ..Default::default()
        });
        let mut source_anchor = SubPatternAnchor::new();
        source_anchor.push_atom(AtomId(7), AtomId(3));
        let edit = ConstraintEdit::new(
            Constraint::Molecule(MoleculeConstraint::SubPattern {
                anchor: source_anchor,
                pattern: Box::new(pattern.clone()),
            }),
            |_| Some(EntityHandle::Atom(AtomHandle::Id(AtomId(0)))),
        )
        .unwrap();
        let mut expected_anchor = SubPatternAnchor::new();
        expected_anchor.push_atom(AtomId(0), AtomId(3));
        let expected = Constraint::Molecule(MoleculeConstraint::SubPattern {
            anchor: expected_anchor,
            pattern: Box::new(pattern),
        });
        let before = one_atom.clone().build();

        let transaction = one_atom
            .transact(Edits::from_iter([Edit::AddMoleculeConstraint {
                constraint: edit,
            }]))
            .unwrap();
        assert_eq!(one_atom.constraints().as_slice(), &[expected]);

        transaction.rollback(&mut one_atom).unwrap();
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_topology_atom_error(mut one_atom: MoleculeEditor) {
        let before = one_atom.clone().build();
        let mut edits = Edits::new();
        edits.remove_atom(AtomHandle::Id(AtomId(9)));
        let err = one_atom.transact(edits).unwrap_err();
        assert_eq!(
            err,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 9,
                count: 1,
            }
        );
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_topology_bond_error(mut diatomic: MoleculeEditor) {
        let before = diatomic.clone().build();
        let mut edits = Edits::new();
        edits.remove_bond(BondHandle::Id(BondId(9)));
        let err = diatomic.transact(edits).unwrap_err();
        assert_eq!(
            err,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Bond,
                index: 9,
                count: 1,
            }
        );
        assert_eq!(diatomic.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_add_dative_bond_empty_atoms_error(
        mut one_atom: MoleculeEditor,
    ) {
        let err = one_atom
            .transact(Edits::from_iter([Edit::AddDativeBond {
                atoms: vec![],
                ast: DativeBondForm::from_order(1),
            }]))
            .unwrap_err();
        assert!(matches!(err, TransactionError::MalformedEdit(_)));
    }

    #[rstest]
    fn test_molecule_editor_transact_set_bond_field(mut diatomic: MoleculeEditor) {
        diatomic
            .transact(Edits::from_iter([Edit::ModifyBondField {
                id: BondHandle::Id(BondId(0)),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap();
        assert_eq!(diatomic.bond(BondId(0)).ast.order, NumForm::Lit(2));
    }

    #[rstest]
    fn test_molecule_editor_transact_set_bond_field_error(mut diatomic: MoleculeEditor) {
        let err = diatomic
            .transact(Edits::from_iter([Edit::ModifyBondField {
                id: BondHandle::Id(BondId(0)),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(99),
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    // Stereo elements — transactional add/remove + undo, and topology-cascade
    // restore. (D3i behavior, tested against D3j's view-based read path.)

    #[fixture]
    fn stereo_atom_skeleton() -> MoleculeEditor {
        let mut b = MoleculeAst::default().edit();
        for el in [Element::C, Element::F, Element::Cl, Element::Br, Element::I] {
            b.add_atom(AtomForm::from_element(el));
        }
        for t in 1u32..=4 {
            b.add_bond(AtomId(0), AtomId(t), BondForm::from_order(1));
        }
        b
    }

    fn tetrahedral_ligands() -> Vec<StereoLigand> {
        (1u32..=4)
            .map(|t| StereoLigand::new(AtomId(t), StereoLigandKind::Atom))
            .collect()
    }

    #[rstest]
    fn test_molecule_editor_transact_add_stereo_atom(mut stereo_atom_skeleton: MoleculeEditor) {
        let before = stereo_atom_skeleton.clone().build();
        let tx = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::AddStereoAtom {
                site: AtomHandle::Id(AtomId(0)),
                ligands: (1u32..=4)
                    .map(|t| (AtomHandle::Id(AtomId(t)), StereoLigandKind::Atom))
                    .collect(),
                ast: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            }]))
            .unwrap();
        assert_eq!(stereo_atom_skeleton.stereo_atom_count(), 1);
        tx.rollback(&mut stereo_atom_skeleton).unwrap();
        assert_eq!(stereo_atom_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_stereo_atom(mut stereo_atom_skeleton: MoleculeEditor) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let before = stereo_atom_skeleton.clone().build();
        let tx = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    (1u32..=4)
                        .map(|t| (AtomHandle::Id(AtomId(t)), StereoLigandKind::Atom))
                        .collect(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                )],
            }]))
            .unwrap();
        assert_eq!(stereo_atom_skeleton.stereo_atom_count(), 0);
        tx.rollback(&mut stereo_atom_skeleton).unwrap();
        assert_eq!(stereo_atom_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_stereo_atom_error(
        mut stereo_atom_skeleton: MoleculeEditor,
    ) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let err = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    (1u32..=4)
                        .map(|t| (AtomHandle::Id(AtomId(t)), StereoLigandKind::Atom))
                        .collect(),
                    // Wrong recorded coset (Th0 vs the stored Th1).
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                )],
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_topology_removal_restores_stereo_atom(
        mut stereo_atom_skeleton: MoleculeEditor,
    ) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let before = stereo_atom_skeleton.clone().build();
        // Removing a ligand atom cascades the stereo element away.
        let tx = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(1))],
                bonds: Vec::new(),
            }]))
            .unwrap();
        assert_eq!(stereo_atom_skeleton.stereo_atom_count(), 0);
        tx.rollback(&mut stereo_atom_skeleton).unwrap();
        assert_eq!(stereo_atom_skeleton.build(), before);
    }

    #[fixture]
    fn stereo_bond_skeleton() -> MoleculeEditor {
        let mut b = MoleculeAst::default().edit();
        for _ in 0..4 {
            b.add_atom(AtomForm::from_element(Element::C));
        }
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(2));
        b.add_bond(AtomId(2), AtomId(3), BondForm::from_order(1));
        b
    }

    #[rstest]
    fn test_molecule_editor_transact_add_stereo_bond(mut stereo_bond_skeleton: MoleculeEditor) {
        let before = stereo_bond_skeleton.clone().build();
        let tx = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::AddStereoBond {
                site: BondHandle::Id(BondId(1)),
                ligands: vec![
                    (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                    (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                ],
                ast: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            }]))
            .unwrap();
        assert_eq!(stereo_bond_skeleton.stereo_bond_count(), 1);
        tx.rollback(&mut stereo_bond_skeleton).unwrap();
        assert_eq!(stereo_bond_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_stereo_bond(mut stereo_bond_skeleton: MoleculeEditor) {
        stereo_bond_skeleton.add_stereo_bond(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let before = stereo_bond_skeleton.clone().build();
        let tx = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(1)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                )],
            }]))
            .unwrap();
        assert_eq!(stereo_bond_skeleton.stereo_bond_count(), 0);
        tx.rollback(&mut stereo_bond_skeleton).unwrap();
        assert_eq!(stereo_bond_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_stereo_bond_error(
        mut stereo_bond_skeleton: MoleculeEditor,
    ) {
        stereo_bond_skeleton.add_stereo_bond(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let error = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(1)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                    ],
                    // Wrong recorded coset (Ct0 vs the stored Ct1).
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                )],
            }]))
            .unwrap_err();

        assert_eq!(error, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_stereo_atom_field(
        mut stereo_atom_skeleton: MoleculeEditor,
    ) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let before = stereo_atom_skeleton.clone().build();
        let tx = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::ModifyStereoAtomField {
                id: StereoAtomHandle::Id(StereoAtomId(0)),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1),
                    ),
                    new: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0),
                    ),
                },
            }]))
            .unwrap();
        assert_eq!(
            stereo_atom_skeleton
                .stereo_atom(StereoAtomId(0))
                .ast
                .configuration,
            StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0),),
        );
        tx.rollback(&mut stereo_atom_skeleton).unwrap();
        assert_eq!(stereo_atom_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_stereo_atom_field_error(
        mut stereo_atom_skeleton: MoleculeEditor,
    ) {
        stereo_atom_skeleton.add_stereo_atom(
            AtomId(0),
            tetrahedral_ligands(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
        let err = stereo_atom_skeleton
            .transact(Edits::from_iter([Edit::ModifyStereoAtomField {
                id: StereoAtomHandle::Id(StereoAtomId(0)),
                change: StereoAtomFieldChange::Configuration {
                    // Wrong recorded coset (Th0 vs the stored Th1).
                    old: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0),
                    ),
                    new: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1),
                    ),
                },
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_stereo_bond_field(
        mut stereo_bond_skeleton: MoleculeEditor,
    ) {
        stereo_bond_skeleton.add_stereo_bond(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let before = stereo_bond_skeleton.clone().build();
        let tx = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::ModifyStereoBondField {
                id: StereoBondHandle::Id(StereoBondId(0)),
                change: StereoBondFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                    new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
                },
            }]))
            .unwrap();
        assert_eq!(
            stereo_bond_skeleton
                .stereo_bond(StereoBondId(0))
                .ast
                .configuration,
            StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
        );
        tx.rollback(&mut stereo_bond_skeleton).unwrap();
        assert_eq!(stereo_bond_skeleton.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_stereo_bond_field_error(
        mut stereo_bond_skeleton: MoleculeEditor,
    ) {
        stereo_bond_skeleton.add_stereo_bond(
            BondId(1),
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
        let err = stereo_bond_skeleton
            .transact(Edits::from_iter([Edit::ModifyStereoBondField {
                id: StereoBondHandle::Id(StereoBondId(0)),
                change: StereoBondFieldChange::Configuration {
                    // Wrong recorded coset (vs the stored 1).
                    old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                },
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[fixture]
    fn diatomic_with_overlays() -> MoleculeEditor {
        let mut b = MoleculeAst::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::N));
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b.add_dative_bond(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1));
        b.add_aromatic_system(vec![AtomId(0), AtomId(1)], AromaticSystemForm::default());
        b.add_multicenter_bond(vec![AtomId(0), AtomId(1)], MulticenterBondForm::default());
        b.add_noncovalent_bond(
            [AtomId(0), AtomId(1)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        b
    }

    #[fixture]
    fn batched_overlays() -> MoleculeEditor {
        let mut editor = MoleculeAst::default().edit();
        for _ in 0..6 {
            editor.add_atom(AtomForm::from_element(Element::C));
        }
        for index in 0..3_u32 {
            let first = AtomId(index * 2);
            let second = AtomId(index * 2 + 1);
            let bond = editor.add_bond(first, second, BondForm::from_order(1));
            editor.add_dative_bond(vec![first], second, DativeBondForm::from_order(1));
            editor.add_aromatic_system(vec![first, second], AromaticSystemForm::default());
            editor.add_multicenter_bond(vec![first, second], MulticenterBondForm::default());
            editor.add_noncovalent_bond(
                [first, second],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            );
            editor.add_stereo_atom(
                first,
                vec![StereoLigand::new(second, StereoLigandKind::Atom)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            );
            editor.add_stereo_bond(
                bond,
                vec![
                    StereoLigand::new(first, StereoLigandKind::Atom),
                    StereoLigand::new(second, StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            );
        }
        editor
    }

    #[rstest]
    #[case::dative_first(EntityKind::DativeBond, 0)]
    #[case::dative_middle(EntityKind::DativeBond, 1)]
    #[case::dative_last(EntityKind::DativeBond, 2)]
    #[case::aromatic_first(EntityKind::AromaticSystem, 0)]
    #[case::aromatic_middle(EntityKind::AromaticSystem, 1)]
    #[case::aromatic_last(EntityKind::AromaticSystem, 2)]
    #[case::multicenter_first(EntityKind::MulticenterBond, 0)]
    #[case::multicenter_middle(EntityKind::MulticenterBond, 1)]
    #[case::multicenter_last(EntityKind::MulticenterBond, 2)]
    #[case::noncovalent_first(EntityKind::NoncovalentBond, 0)]
    #[case::noncovalent_middle(EntityKind::NoncovalentBond, 1)]
    #[case::noncovalent_last(EntityKind::NoncovalentBond, 2)]
    #[case::stereo_atom_first(EntityKind::StereoAtom, 0)]
    #[case::stereo_atom_middle(EntityKind::StereoAtom, 1)]
    #[case::stereo_atom_last(EntityKind::StereoAtom, 2)]
    #[case::stereo_bond_first(EntityKind::StereoBond, 0)]
    #[case::stereo_bond_middle(EntityKind::StereoBond, 1)]
    #[case::stereo_bond_last(EntityKind::StereoBond, 2)]
    fn test_molecule_editor_transact_remove_overlays_error(
        mut batched_overlays: MoleculeEditor,
        #[case] kind: EntityKind,
        #[case] invalid_position: usize,
    ) {
        let before = batched_overlays.clone().build();
        let edit = match kind {
            EntityKind::DativeBond => Edit::RemoveDativeBonds {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            DativeBondHandle::Id(DativeBondId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            vec![
                                AtomHandle::Id(AtomId(index * 2)),
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                            ],
                            DativeBondForm::from_order(1),
                        )
                    })
                    .collect(),
            },
            EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            AromaticSystemHandle::Id(AromaticSystemId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            vec![
                                AtomHandle::Id(AtomId(index * 2)),
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                            ],
                            AromaticSystemForm::default(),
                        )
                    })
                    .collect(),
            },
            EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            MulticenterBondHandle::Id(MulticenterBondId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            vec![
                                AtomHandle::Id(AtomId(index * 2)),
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                            ],
                            MulticenterBondForm::default(),
                        )
                    })
                    .collect(),
            },
            EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            NoncovalentBondHandle::Id(NoncovalentBondId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            [
                                AtomHandle::Id(AtomId(index * 2)),
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                            ],
                            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                        )
                    })
                    .collect(),
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            StereoAtomHandle::Id(StereoAtomId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            AtomHandle::Id(AtomId(index * 2)),
                            vec![(
                                AtomHandle::Id(AtomId(index * 2 + 1)),
                                StereoLigandKind::Atom,
                            )],
                            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                        )
                    })
                    .collect(),
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: (0..3_u32)
                    .map(|index| {
                        (
                            StereoBondHandle::Id(StereoBondId(
                                if index as usize == invalid_position {
                                    9
                                } else {
                                    index
                                },
                            )),
                            BondHandle::Id(BondId(index)),
                            vec![
                                (AtomHandle::Id(AtomId(index * 2)), StereoLigandKind::Atom),
                                (
                                    AtomHandle::Id(AtomId(index * 2 + 1)),
                                    StereoLigandKind::Atom,
                                ),
                            ],
                            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                        )
                    })
                    .collect(),
            },
            EntityKind::Atom | EntityKind::Bond => unreachable!(),
        };

        assert_eq!(
            batched_overlays
                .transact(Edits::from_iter([edit]))
                .unwrap_err(),
            TransactionError::HandleOutOfRange {
                kind,
                index: 9,
                count: 3,
            }
        );
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    #[case::atom(EntityKind::Atom)]
    #[case::bond(EntityKind::Bond)]
    #[case::dative_bond(EntityKind::DativeBond)]
    #[case::aromatic_system(EntityKind::AromaticSystem)]
    #[case::multicenter_bond(EntityKind::MulticenterBond)]
    #[case::noncovalent_bond(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_molecule_editor_transact_duplicate_removal_error(
        mut batched_overlays: MoleculeEditor,
        #[case] kind: EntityKind,
    ) {
        let before = batched_overlays.clone().build();
        let edit = match kind {
            EntityKind::Atom => Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            EntityKind::Bond => Edit::RemoveTopology {
                atoms: Vec::new(),
                bonds: vec![BondHandle::Id(BondId(0)), BondHandle::Id(BondId(0))],
            },
            EntityKind::DativeBond => Edit::RemoveDativeBonds {
                removes: vec![
                    (
                        DativeBondHandle::Id(DativeBondId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        DativeBondForm::from_order(1),
                    ),
                    (
                        DativeBondHandle::Id(DativeBondId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        DativeBondForm::from_order(1),
                    ),
                ],
            },
            EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
                removes: vec![
                    (
                        AromaticSystemHandle::Id(AromaticSystemId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        AromaticSystemForm::default(),
                    ),
                    (
                        AromaticSystemHandle::Id(AromaticSystemId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        AromaticSystemForm::default(),
                    ),
                ],
            },
            EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
                removes: vec![
                    (
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        MulticenterBondForm::default(),
                    ),
                    (
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        MulticenterBondForm::default(),
                    ),
                ],
            },
            EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
                removes: vec![
                    (
                        NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                        [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    ),
                    (
                        NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                        [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    ),
                ],
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: vec![
                    (
                        StereoAtomHandle::Id(StereoAtomId(0)),
                        AtomHandle::Id(AtomId(0)),
                        vec![(AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom)],
                        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                    ),
                    (
                        StereoAtomHandle::Id(StereoAtomId(0)),
                        AtomHandle::Id(AtomId(0)),
                        vec![(AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom)],
                        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                    ),
                ],
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: vec![
                    (
                        StereoBondHandle::Id(StereoBondId(0)),
                        BondHandle::Id(BondId(0)),
                        vec![
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                            (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                        ],
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                    ),
                    (
                        StereoBondHandle::Id(StereoBondId(0)),
                        BondHandle::Id(BondId(0)),
                        vec![
                            (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                            (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                        ],
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                    ),
                ],
            },
        };

        assert_eq!(
            batched_overlays
                .transact(Edits::from_iter([edit]))
                .unwrap_err(),
            TransactionError::DuplicateRemoval { kind }
        );
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    #[case::dative_bond(EntityKind::DativeBond)]
    #[case::aromatic_system(EntityKind::AromaticSystem)]
    #[case::multicenter_bond(EntityKind::MulticenterBond)]
    #[case::noncovalent_bond(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_molecule_editor_transact_handle_removed_error_cascade(
        mut batched_overlays: MoleculeEditor,
        #[case] kind: EntityKind,
    ) {
        let before = batched_overlays.clone().build();
        let mut edits = Edits::from_iter([Edit::RemoveTopology {
            atoms: vec![AtomHandle::Id(AtomId(0))],
            bonds: Vec::new(),
        }]);
        edits.push(match kind {
            EntityKind::DativeBond => Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    DativeBondForm::from_order(1),
                )],
            },
            EntityKind::AromaticSystem => Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                )],
            },
            EntityKind::MulticenterBond => Edit::RemoveMulticenterBonds {
                removes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    MulticenterBondForm::default(),
                )],
            },
            EntityKind::NoncovalentBond => Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
            },
            EntityKind::StereoAtom => Edit::RemoveStereoAtoms {
                removes: vec![(
                    StereoAtomHandle::Id(StereoAtomId(0)),
                    AtomHandle::Id(AtomId(0)),
                    vec![(AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom)],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                )],
            },
            EntityKind::StereoBond => Edit::RemoveStereoBonds {
                removes: vec![(
                    StereoBondHandle::Id(StereoBondId(0)),
                    BondHandle::Id(BondId(0)),
                    vec![
                        (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                        (AtomHandle::Id(AtomId(1)), StereoLigandKind::Atom),
                    ],
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                )],
            },
            EntityKind::Atom | EntityKind::Bond => unreachable!(),
        });

        assert_eq!(
            batched_overlays.transact(edits).unwrap_err(),
            TransactionError::HandleRemoved { kind, index: 0 }
        );
        assert_eq!(batched_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_dative_bond_field(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyDativeBondField {
                id: DativeBondHandle::Id(DativeBondId(0)),
                change: DativeBondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .dative_bond(DativeBondId(0))
                .ast
                .order,
            NumForm::Lit(2),
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_aromatic_system_field(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyAromaticSystemField {
                id: AromaticSystemHandle::Id(AromaticSystemId(0)),
                change: AromaticSystemFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .aromatic_system(AromaticSystemId(0))
                .ast
                .charge,
            NumForm::Lit(1),
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_multicenter_bond_field(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                change: MulticenterBondFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(-1),
                },
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .multicenter_bond(MulticenterBondId(0))
                .ast
                .charge,
            NumForm::Lit(-1),
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_noncovalent_bond_field(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyNoncovalentBondField {
                id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                    new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                },
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .noncovalent_bond(NoncovalentBondId(0))
                .ast
                .kind,
            NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_add_dative_bond(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::AddDativeBond {
                atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: DativeBondForm::from_order(1),
            }]))
            .unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedDativeBond(added)] if added.id == DativeBondId(0)
        ));
        assert_eq!(diatomic.dative_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_add_aromatic_system(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::AddAromaticSystem {
                atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: AromaticSystemForm::default(),
            }]))
            .unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedAromaticSystem(added)] if added.id == AromaticSystemId(0)
        ));
        assert_eq!(diatomic.aromatic_system_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_add_multicenter_bond(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::AddMulticenterBond {
                atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: MulticenterBondForm::default(),
            }]))
            .unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedMulticenterBond(added)] if added.id == MulticenterBondId(0)
        ));
        assert_eq!(diatomic.multicenter_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_add_noncovalent_bond(mut diatomic: MoleculeEditor) {
        let tx = diatomic
            .transact(Edits::from_iter([Edit::AddNoncovalentBond {
                atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            }]))
            .unwrap();
        assert!(matches!(
            tx.undos(),
            [Undo::RemoveAddedNoncovalentBond(added)] if added.id == NoncovalentBondId(0)
        ));
        assert_eq!(diatomic.noncovalent_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_dative_bond(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let before = diatomic_with_overlays.clone().build();
        let transaction = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    DativeBondForm {
                        order: NumForm::Lit(1),
                        constraints: Default::default(),
                    },
                )],
            }]))
            .unwrap();
        assert_eq!(diatomic_with_overlays.dative_bond_count(), 0);
        transaction.rollback(&mut diatomic_with_overlays).unwrap();
        assert_eq!(diatomic_with_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_dative_bond_atoms_mismatch_error(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let err = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(0))], // wrong order
                    DativeBondForm {
                        order: NumForm::Lit(1),
                        constraints: Default::default(),
                    },
                )],
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
        assert_eq!(diatomic_with_overlays.dative_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_aromatic_system(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let before = diatomic_with_overlays.clone().build();
        let transaction = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                )],
            }]))
            .unwrap();
        assert_eq!(diatomic_with_overlays.aromatic_system_count(), 0);
        transaction.rollback(&mut diatomic_with_overlays).unwrap();
        assert_eq!(diatomic_with_overlays.build(), before);
    }

    // Batch removal of non-contiguous same-kind ids (0 and 2) in one edit: ids resolve against the
    // pre-removal state and compact once, so the survivor (former id 1) remaps to id 0. A single-id
    // sequence would stale id 2 after removing id 0.
    #[rstest]
    fn test_molecule_editor_transact_remove_aromatic_systems() {
        let mut b = MoleculeAst::default().edit();
        for _ in 0..6 {
            b.add_atom(AtomForm::from_element(Element::C));
        }
        b.add_aromatic_system(vec![AtomId(0), AtomId(1)], AromaticSystemForm::default());
        b.add_aromatic_system(vec![AtomId(2), AtomId(3)], AromaticSystemForm::default());
        b.add_aromatic_system(vec![AtomId(4), AtomId(5)], AromaticSystemForm::default());
        b.transact(Edits::from_iter([Edit::RemoveAromaticSystems {
            removes: vec![
                (
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                ),
                (
                    AromaticSystemHandle::Id(AromaticSystemId(2)),
                    vec![AtomHandle::Id(AtomId(4)), AtomHandle::Id(AtomId(5))],
                    AromaticSystemForm::default(),
                ),
            ],
        }]))
        .unwrap();
        assert_eq!(b.aromatic_system_count(), 1);
        assert_eq!(
            b.aromatic_system(AromaticSystemId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2), AtomId(3)],
        );
    }

    // Rolling back an aromatic-system removal restores a molecule constraint the removal dropped
    // (`dropped`) or remapped (`remapped`) — the overlay-remove undo captures the constraint cascade.
    #[rstest]
    #[case::dropped(AromaticSystemId(0), 0)]
    #[case::remapped(AromaticSystemId(1), 1)]
    fn test_molecule_editor_transact_remove_aromatic_system_rollback(
        #[case] constrained: AromaticSystemId,
        #[case] forward_constraint_count: usize,
    ) {
        let mut b = MoleculeAst::default().edit();
        for _ in 0..6 {
            b.add_atom(AtomForm::from_element(Element::C));
        }
        b.add_aromatic_system(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        );
        b.add_aromatic_system(
            vec![AtomId(3), AtomId(4), AtomId(5)],
            AromaticSystemForm::default(),
        );
        b.transact(Edits::from_iter([Edit::AddMoleculeConstraint {
            constraint: Constraint::AromaticSystem(
                constrained,
                AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6)),
            )
            .into(),
        }]))
        .unwrap();
        let before = b.clone().build();

        let tx = b
            .transact(Edits::from_iter([Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    vec![
                        AtomHandle::Id(AtomId(0)),
                        AtomHandle::Id(AtomId(1)),
                        AtomHandle::Id(AtomId(2)),
                    ],
                    AromaticSystemForm::default(),
                )],
            }]))
            .unwrap();
        assert_eq!(b.aromatic_system_count(), 1);
        assert_eq!(b.constraints().iter().count(), forward_constraint_count);

        tx.rollback(&mut b).unwrap();
        assert_eq!(b.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_multicenter_bond(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let before = diatomic_with_overlays.clone().build();
        let transaction = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveMulticenterBonds {
                removes: vec![(
                    MulticenterBondHandle::Id(MulticenterBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    MulticenterBondForm::default(),
                )],
            }]))
            .unwrap();
        assert_eq!(diatomic_with_overlays.multicenter_bond_count(), 0);
        transaction.rollback(&mut diatomic_with_overlays).unwrap();
        assert_eq!(diatomic_with_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_noncovalent_bond(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let before = diatomic_with_overlays.clone().build();
        let transaction = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
            }]))
            .unwrap();
        assert_eq!(diatomic_with_overlays.noncovalent_bond_count(), 0);
        transaction.rollback(&mut diatomic_with_overlays).unwrap();
        assert_eq!(diatomic_with_overlays.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_remove_noncovalent_bond_form_mismatch_error(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        let err = diatomic_with_overlays
            .transact(Edits::from_iter([Edit::RemoveNoncovalentBonds {
                removes: vec![(
                    NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic), // wrong
                )],
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_bond_constraint_value_bearing(
        mut diatomic: MoleculeEditor,
    ) {
        diatomic
            .transact(Edits::from_iter([Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(0)),
                old: None,
                new: Some(BondConstraintForm::cis_trans_stereo(
                    CisTransStereoForm::NotStereo,
                )),
            }]))
            .unwrap();
        assert_eq!(
            diatomic
                .bond(BondId(0))
                .ast
                .constraints
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            vec![BondConstraintForm::cis_trans_stereo(
                CisTransStereoForm::NotStereo
            )],
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_add_bond_constraint(mut diatomic: MoleculeEditor) {
        diatomic
            .transact(Edits::from_iter([Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(0)),
                old: None,
                new: Some(BondConstraintForm::ring_membership(RingScope::Size(5), 1)),
            }]))
            .unwrap();
        assert!(diatomic
            .bond(BondId(0))
            .ast
            .constraints
            .iter()
            .any(|c| *c == BondConstraintForm::ring_membership(RingScope::Size(5), 1)));
    }

    #[rstest]
    fn test_molecule_editor_transact_modify_bond_constraint_absent_error(
        mut diatomic: MoleculeEditor,
    ) {
        let err = diatomic
            .transact(Edits::from_iter([Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(0)),
                old: Some(BondConstraintForm::ring_membership(RingScope::Size(5), 1)),
                new: None,
            }]))
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_editor_transact_set_dative_bond_constraint(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyDativeBondConstraint {
                id: DativeBondHandle::Id(DativeBondId(0)),
                old: None,
                new: Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            }]))
            .unwrap();
        assert!(diatomic_with_overlays
            .dative_bond(DativeBondId(0))
            .ast
            .constraints
            .iter()
            .any(|c| *c == DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))));
    }

    #[rstest]
    fn test_molecule_editor_transact_set_aromatic_system_constraint(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyAromaticSystemConstraint {
                id: AromaticSystemHandle::Id(AromaticSystemId(0)),
                old: None,
                new: Some(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))),
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .aromatic_system(AromaticSystemId(0))
                .ast
                .constraints
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            vec![AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))],
        );
    }

    #[rstest]
    fn test_molecule_editor_transact_set_multicenter_bond_constraint(
        mut diatomic_with_overlays: MoleculeEditor,
    ) {
        diatomic_with_overlays
            .transact(Edits::from_iter([Edit::ModifyMulticenterBondConstraint {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                old: None,
                new: Some(MulticenterBondConstraintForm::ElectronCount(NumForm::Lit(
                    2,
                ))),
            }]))
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .multicenter_bond(MulticenterBondId(0))
                .ast
                .constraints
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            vec![MulticenterBondConstraintForm::ElectronCount(NumForm::Lit(
                2
            ))],
        );
    }

    #[fixture]
    fn triatomic_with_overlays() -> MoleculeEditor {
        let mut b = MoleculeAst::default().edit();
        b.add_atom(AtomForm::from_element(Element::C));
        b.add_atom(AtomForm::from_element(Element::N));
        b.add_atom(AtomForm::from_element(Element::O));
        b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
        b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(1));
        b.add_dative_bond(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1));
        b.add_aromatic_system(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        );
        b.add_multicenter_bond(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        );
        b.add_noncovalent_bond(
            [AtomId(0), AtomId(2)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        b
    }

    #[rstest]
    fn test_transaction_append(mut diatomic: MoleculeEditor) {
        let before = diatomic.clone().build();
        let first = diatomic
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap();
        let second = diatomic
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            }]))
            .unwrap();
        let expected_undos = [first.undos(), second.undos()].concat();

        let mut combined = Transaction::default();
        combined.append(first);
        combined.append(second);

        assert_eq!(combined.undos(), expected_undos);
        combined.rollback(&mut diatomic).unwrap();
        assert_eq!(diatomic.build(), before);
    }

    #[rstest]
    fn test_transaction_append_error(mut diatomic: MoleculeEditor) {
        let before = diatomic.clone().build();
        let first = diatomic
            .transact(Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]))
            .unwrap();
        let second = diatomic
            .transact(Edits::from_iter([Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old: None,
                new: Some(AtomConstraintForm::degree(1)),
            }]))
            .unwrap();
        let third = diatomic
            .transact(Edits::from_iter([Edit::AddDativeBond {
                atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: DativeBondForm::from_order(1),
            }]))
            .unwrap();
        let expected_undos = [first.undos(), second.undos(), third.undos()].concat();

        let mut combined = Transaction::default();
        combined.append(first);
        combined.append(second);
        combined.append(third);

        let materialized = diatomic.clone().build();
        let mut rejected = Edits::from_iter([Edit::ModifyAtomField {
            id: AtomHandle::Id(AtomId(0)),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        }]);
        rejected.remove_atom(AtomHandle::Id(AtomId(99)));
        let error = diatomic.transact(rejected).unwrap_err();
        assert_eq!(
            error,
            TransactionError::HandleOutOfRange {
                kind: EntityKind::Atom,
                index: 99,
                count: 2,
            }
        );
        assert_eq!(diatomic.clone().build(), materialized);
        assert_eq!(combined.undos(), expected_undos);

        combined.rollback(&mut diatomic).unwrap();
        assert_eq!(diatomic.build(), before);
    }

    enum RollbackCase {
        RemoveTopology,
        RemoveBond,
        AddTopology,
        Field,
        AddOverlay,
        RemoveOverlay,
        Constraint,
        CascadedConstraints,
    }

    #[rstest]
    #[case::remove_topology(RollbackCase::RemoveTopology)]
    #[case::remove_bond(RollbackCase::RemoveBond)]
    #[case::add_topology(RollbackCase::AddTopology)]
    #[case::field(RollbackCase::Field)]
    #[case::add_overlay(RollbackCase::AddOverlay)]
    #[case::remove_overlay(RollbackCase::RemoveOverlay)]
    #[case::constraint(RollbackCase::Constraint)]
    #[case::cascade(RollbackCase::CascadedConstraints)]
    fn test_transaction_rollback(#[case] case: RollbackCase) {
        let mut editor = match case {
            RollbackCase::AddTopology => MoleculeAst::default().edit(),
            RollbackCase::Field | RollbackCase::Constraint => {
                let mut b = MoleculeAst::default().edit();
                b.add_atom(AtomForm::from_element(Element::C));
                b
            }
            RollbackCase::AddOverlay => {
                let mut b = MoleculeAst::default().edit();
                b.add_atom(AtomForm::from_element(Element::C));
                b.add_atom(AtomForm::from_element(Element::C));
                b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
                b
            }
            RollbackCase::RemoveBond => {
                let mut b = MoleculeAst::default().edit();
                b.add_atom(AtomForm::from_element(Element::C));
                b.add_atom(AtomForm::from_element(Element::N));
                b.add_atom(AtomForm::from_element(Element::O));
                b.add_bond(AtomId(0), AtomId(1), BondForm::from_order(1));
                b.add_bond(AtomId(1), AtomId(2), BondForm::from_order(2));
                b
            }
            RollbackCase::CascadedConstraints => {
                let mut b = MoleculeAst::default().edit();
                b.add_atom(AtomForm::from_element(Element::C));
                b.add_atom(AtomForm::from_element(Element::N));
                b.push_constraint(Constraint::Atom(AtomId(1), AtomConstraintForm::degree(3)));
                b
            }
            RollbackCase::RemoveTopology | RollbackCase::RemoveOverlay => triatomic_with_overlays(),
        };
        let before = editor.clone().build();
        let edits = match case {
            RollbackCase::RemoveTopology => Edits::from_iter([Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(1))],
                bonds: vec![],
            }]),
            RollbackCase::RemoveBond => Edits::from_iter([Edit::RemoveTopology {
                atoms: Vec::new(),
                bonds: vec![BondHandle::Id(BondId(0))],
            }]),
            RollbackCase::AddTopology => Edits::from_iter([
                Edit::AddAtoms {
                    atoms: vec![
                        AtomForm::from_element(Element::C),
                        AtomForm::from_element(Element::O),
                    ],
                },
                Edit::AddBonds {
                    bonds: vec![AddBond {
                        endpoints: [AtomHandle::New(0), AtomHandle::New(1)],
                        ast: BondForm::from_order(2),
                    }],
                },
            ]),
            RollbackCase::Field => Edits::from_iter([Edit::ModifyAtomField {
                id: AtomHandle::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: NumForm::default(),
                    new: NumForm::Lit(1),
                },
            }]),
            RollbackCase::AddOverlay => Edits::from_iter([Edit::AddDativeBond {
                atoms: vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: DativeBondForm::from_order(1),
            }]),
            RollbackCase::RemoveOverlay => Edits::from_iter([Edit::RemoveDativeBonds {
                removes: vec![(
                    DativeBondHandle::Id(DativeBondId(0)),
                    vec![AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    DativeBondForm {
                        order: NumForm::Lit(1),
                        constraints: Default::default(),
                    },
                )],
            }]),
            RollbackCase::Constraint => Edits::from_iter([Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old: None,
                new: Some(AtomConstraintForm::ring_membership(RingScope::Size(5), 1)),
            }]),
            RollbackCase::CascadedConstraints => Edits::from_iter([Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(1))],
                bonds: Vec::new(),
            }]),
        };
        let tx = editor.transact(edits).unwrap();
        tx.rollback(&mut editor).unwrap();
        assert_eq!(editor.build(), before);
    }

    #[rstest]
    fn test_transaction_rollback_empty(mut one_atom: MoleculeEditor) {
        let before = one_atom.clone().build();
        Transaction::default().rollback(&mut one_atom).unwrap();
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    #[case::atom(EntityKind::Atom)]
    #[case::bond(EntityKind::Bond)]
    #[case::dative_bond(EntityKind::DativeBond)]
    #[case::aromatic_system(EntityKind::AromaticSystem)]
    #[case::multicenter_bond(EntityKind::MulticenterBond)]
    #[case::noncovalent_bond(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_transaction_rollback_field_receiver(
        mut empty: MoleculeEditor,
        #[case] kind: EntityKind,
    ) {
        let undo = match kind {
            EntityKind::Atom => Undo::ModifyAtomField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::default(),
                },
            },
            EntityKind::Bond => Undo::ModifyBondField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(2),
                    new: NumForm::Lit(1),
                },
            },
            EntityKind::DativeBond => Undo::ModifyDativeBondField {
                id: DativeBondId(0),
                change: DativeBondFieldChange::Order {
                    old: NumForm::Lit(2),
                    new: NumForm::Lit(1),
                },
            },
            EntityKind::AromaticSystem => Undo::ModifyAromaticSystemField {
                id: AromaticSystemId(0),
                change: AromaticSystemFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::default(),
                },
            },
            EntityKind::MulticenterBond => Undo::ModifyMulticenterBondField {
                id: MulticenterBondId(0),
                change: MulticenterBondFieldChange::Charge {
                    old: NumForm::Lit(1),
                    new: NumForm::default(),
                },
            },
            EntityKind::NoncovalentBond => Undo::ModifyNoncovalentBondField {
                id: NoncovalentBondId(0),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::Ionic),
                    new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                },
            },
            EntityKind::StereoAtom => Undo::ModifyStereoAtomField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0),
                    ),
                    new: StereoConfigurationForm::kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1),
                    ),
                },
            },
            EntityKind::StereoBond => Undo::ModifyStereoBondField {
                id: StereoBondId(0),
                change: StereoBondFieldChange::Configuration {
                    old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                },
            },
        };

        assert_eq!(
            (Transaction { undo: vec![undo] }).rollback(&mut empty),
            Err(TransactionError::RollbackStateMismatch),
        );
    }

    #[rstest]
    #[case::atom(EntityKind::Atom)]
    #[case::bond(EntityKind::Bond)]
    #[case::dative_bond(EntityKind::DativeBond)]
    #[case::aromatic_system(EntityKind::AromaticSystem)]
    #[case::multicenter_bond(EntityKind::MulticenterBond)]
    #[case::noncovalent_bond(EntityKind::NoncovalentBond)]
    #[case::stereo_atom(EntityKind::StereoAtom)]
    #[case::stereo_bond(EntityKind::StereoBond)]
    fn test_transaction_rollback_constraint_receiver(
        mut empty: MoleculeEditor,
        #[case] kind: EntityKind,
    ) {
        let edit = match kind {
            EntityKind::Atom => Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old: None,
                new: None,
            },
            EntityKind::Bond => Edit::ModifyBondConstraint {
                id: BondHandle::Id(BondId(0)),
                old: None,
                new: None,
            },
            EntityKind::DativeBond => Edit::ModifyDativeBondConstraint {
                id: DativeBondHandle::Id(DativeBondId(0)),
                old: None,
                new: None,
            },
            EntityKind::AromaticSystem => Edit::ModifyAromaticSystemConstraint {
                id: AromaticSystemHandle::Id(AromaticSystemId(0)),
                old: None,
                new: None,
            },
            EntityKind::MulticenterBond => Edit::ModifyMulticenterBondConstraint {
                id: MulticenterBondHandle::Id(MulticenterBondId(0)),
                old: None,
                new: None,
            },
            EntityKind::NoncovalentBond => Edit::ModifyNoncovalentBondConstraint {
                id: NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                old: None,
                new: None,
            },
            EntityKind::StereoAtom => Edit::ModifyStereoAtomConstraint {
                id: StereoAtomHandle::Id(StereoAtomId(0)),
                kind: None,
                old: None,
                new: None,
            },
            EntityKind::StereoBond => Edit::ModifyStereoBondConstraint {
                id: StereoBondHandle::Id(StereoBondId(0)),
                kind: None,
                old: None,
                new: None,
            },
        };

        assert_eq!(
            (Transaction {
                undo: vec![Undo::ApplyEdit(Box::new(edit))],
            })
            .rollback(&mut empty),
            Err(TransactionError::RollbackStateMismatch),
        );
    }

    #[rstest]
    fn test_transaction_rollback_added_topology_duplicate(mut one_atom: MoleculeEditor) {
        let before = one_atom.clone().build();
        let added = AddedAtom {
            id: AtomId(0),
            ast: AtomForm::from_element(Element::C),
        };
        let transaction = Transaction {
            undo: vec![Undo::RemoveAddedTopology {
                atoms: vec![added.clone(), added],
                bonds: Vec::new(),
            }],
        };

        assert_eq!(
            transaction.rollback(&mut one_atom),
            Err(TransactionError::RollbackStateMismatch),
        );
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    fn test_transaction_rollback_reconstruction_slot(mut empty: MoleculeEditor) {
        let transaction = Transaction {
            undo: vec![Undo::RestoreRemovedAromaticSystems {
                removed: vec![RemovedAromaticSystem {
                    id: AromaticSystemId(1),
                    atoms: Vec::new(),
                    ast: AromaticSystemForm::default(),
                }],
                undo_compaction: IdCompaction::empty().undo_compaction(),
                cascade: CascadedConstraints::default(),
            }],
        };

        assert_eq!(
            transaction.rollback(&mut empty),
            Err(TransactionError::RollbackStateMismatch),
        );
    }

    #[rstest]
    fn test_transaction_rollback_compaction_dimension(mut one_atom: MoleculeEditor) {
        let compaction = IdCompaction::empty();
        let transaction = Transaction {
            undo: vec![Undo::RestoreRemovedTopology {
                atoms: vec![RemovedAtom {
                    id: AtomId(0),
                    ast: AtomForm::from_element(Element::N),
                }],
                bonds: Vec::new(),
                overlays: RemovedOverlays::default(),
                undo_compaction: compaction.undo_compaction(),
                compaction,
                cascade: CascadedConstraints::default(),
            }],
        };

        assert_eq!(
            transaction.rollback(&mut one_atom),
            Err(TransactionError::RollbackStateMismatch),
        );
    }

    #[rstest]
    fn test_transaction_rollback_molecule_constraint_order(mut one_atom: MoleculeEditor) {
        let repeated = Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1));
        let middle = Constraint::Atom(AtomId(0), AtomConstraintForm::valence(4));
        one_atom.push_constraint(repeated.clone());
        one_atom.push_constraint(middle.clone());
        one_atom.push_constraint(repeated.clone());
        let before = one_atom.clone().build();

        let transaction = one_atom
            .transact(Edits::from_iter([Edit::RemoveMoleculeConstraint {
                constraint: repeated.clone().into(),
            }]))
            .unwrap();
        assert_eq!(one_atom.constraints().as_slice(), &[repeated, middle]);

        transaction.rollback(&mut one_atom).unwrap();
        assert_eq!(one_atom.build(), before);
    }

    #[rstest]
    fn test_molecule_editor_transact_unchecked(mut empty: MoleculeEditor) {
        let mut edits = Edits::new();
        edits.add_atom(AtomForm::from_element(Element::C));
        empty.transact_unchecked(edits);
        assert_eq!(empty.atom_count(), 1);
        assert_eq!(
            empty.atom(AtomId(0)).ast.element,
            ElementForm::Lit(Element::C)
        );
    }

    #[rstest]
    #[should_panic(expected = "invalid unchecked transaction edit")]
    fn test_molecule_editor_transact_unchecked_error(mut empty: MoleculeEditor) {
        let mut edits = Edits::new();
        edits.remove_atom(AtomHandle::Id(AtomId(0)));
        empty.transact_unchecked(edits);
    }
}
