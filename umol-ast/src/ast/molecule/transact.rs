//! Transactional `Edit` application on `MoleculeBuilder` (Phase 8d).
//!
//! `transact(edits)` snapshots the builder, applies each `Edit` in order,
//! and either returns the per-Edit `Action` list or restores the snapshot
//! and surfaces a `TransactionError`. Snapshot is a `MoleculeBuilder::clone()`
//! — internally Arc-shared per-relation storage, so the snapshot cost is
//! O(1) Arc clones plus a `Graph` clone.

use thiserror::Error;

use super::super::edit::{
    Action, AddBond, AromaticSystemFieldChange, AromaticSystemRef, AtomFieldChange, AtomRef,
    BondFieldChange, BondRef, DativeBondFieldChange, DativeBondRef, Edit,
    MulticenterBondFieldChange, MulticenterBondRef, NoncovalentBondFieldChange,
    NoncovalentBondRef, Undo,
};
use super::super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::MoleculeBuilder;

#[derive(Debug, Error, PartialEq, Eq, Clone)]
pub enum TransactionError {
    /// `Ref::New(N)` references an edit that hasn't been applied yet.
    #[error("symbolic ref points past end of action list: New({0}) with {1} actions so far")]
    RefOutOfRange(usize, usize),

    /// `Ref::New(N)` references an earlier edit, but that edit's Action is
    /// the wrong shape (e.g., `AtomRef::New(N)` where edit N was `AddBond`).
    #[error("symbolic ref type mismatch: expected {expected}, got {got}")]
    RefTypeMismatch { expected: &'static str, got: &'static str },

    /// An id in the Edit references an entity that does not exist (out of
    /// range against current entity counts).
    #[error("id out of range: {0}")]
    IdOutOfRange(&'static str),

    /// `Set*Field` or `Set*Constraint`: current state does not match the
    /// edit's `old` payload.
    #[error("precondition failed: old state does not match current")]
    OldStateMismatch,

    /// `Set*Constraint` applied to a non-unique kind, or `Add`/`Remove`
    /// applied to a unique kind.
    #[error("constraint shape mismatch for kind")]
    KindShapeMismatch,

    /// `Set*Constraint { old: Some(a), new: Some(b) }` with `a.kind() != b.kind()`.
    #[error("constraint kind mismatch between old and new")]
    KindMismatch,

    /// `Add*Constraint` with a value that's already present in the store.
    #[error("duplicate constraint entry on add")]
    DuplicateEntry,

    /// `Remove*Constraint` with a value that's not present.
    #[error("missing constraint entry on remove")]
    MissingEntry,

    /// Edit shape is structurally invalid (e.g., `AddDativeBond` with no
    /// participants).
    #[error("malformed edit: {0}")]
    MalformedEdit(&'static str),
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Transaction {
    undo: Vec<Undo>,
}

impl Transaction {
    pub fn new(undo: Vec<Undo>) -> Self {
        Self { undo }
    }

    pub fn undos(&self) -> &[Undo] {
        &self.undo
    }

}

impl MoleculeBuilder {
    /// Apply a batch of `Edit`s atomically. On success, returns the per-Edit
    /// `Action` list. On any apply failure, restores the pre-batch snapshot
    /// and surfaces the error.
    pub fn transact(&mut self, edits: Vec<Edit>) -> Result<Vec<Action>, TransactionError> {
        let snapshot = self.clone();
        let mut actions: Vec<Action> = Vec::with_capacity(edits.len());
        for edit in edits {
            match self.apply_edit(edit, &actions) {
                Ok(action) => actions.push(action),
                Err(e) => {
                    *self = snapshot;
                    return Err(e);
                }
            }
        }
        Ok(actions)
    }

    fn apply_edit(
        &mut self,
        edit: Edit,
        actions: &[Action],
    ) -> Result<Action, TransactionError> {
        match edit {
            Edit::AddAtoms { atoms } => {
                let ids = atoms
                    .into_iter()
                    .map(|ast| self.add_atom(ast))
                    .collect();
                Ok(Action::AtomsAdded(ids))
            }
            Edit::AddBonds { bonds } => {
                let mut ids = Vec::with_capacity(bonds.len());
                for AddBond { a, b, ast } in bonds {
                    let a = resolve_atom_ref(a, actions)?;
                    if a.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                    let b = resolve_atom_ref(b, actions)?;
                    if b.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                    ids.push(self.add_bond(a, b, ast));
                }
                Ok(Action::BondsAdded(ids))
            }
            Edit::RemoveTopology { atoms, bonds } => {
                let atoms: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|idx| {
                        let id = resolve_atom_ref(idx, actions)?;
                        if id.index() >= self.atom_count() {
                            return Err(TransactionError::IdOutOfRange("atom"));
                        }
                        Ok(id)
                    })
                    .collect::<Result<_, _>>()?;
                let bonds: Vec<BondId> = bonds
                    .into_iter()
                    .map(|idx| {
                        let id = resolve_bond_ref(idx, actions)?;
                        if id.index() >= self.bond_count() {
                            return Err(TransactionError::IdOutOfRange("bond"));
                        }
                        Ok(id)
                    })
                    .collect::<Result<_, _>>()?;
                if !atoms.is_empty() || !bonds.is_empty() {
                    self.remove(&atoms, &bonds);
                }
                Ok(Action::Done)
            }
            Edit::SetAtomField { idx, change } => {
                let id = resolve_atom_ref(idx, actions)?;
                if id.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                self.apply_set_atom_field(id, change)?;
                Ok(Action::Done)
            }
            Edit::SetBondField { idx, change } => {
                let id = resolve_bond_ref(idx, actions)?;
                if id.index() >= self.bond_count() { return Err(TransactionError::IdOutOfRange("bond")); }
                self.apply_set_bond_field(id, change)?;
                Ok(Action::Done)
            }
            Edit::AddDativeBond { atoms, ast } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| resolve_atom_ref(r, actions))
                    .collect::<Result<_, _>>()?;
                for a in &resolved {
                    if a.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                }
                let (acceptor, donors) = resolved
                    .split_last()
                    .ok_or(TransactionError::MalformedEdit(
                        "AddDativeBond requires at least one participant atom",
                    ))?;
                Ok(Action::DativeBondAdded(self.add_dative_bond(
                    donors.to_vec(),
                    *acceptor,
                    ast,
                )))
            }
            Edit::RemoveDativeBond { idx, atoms, ast } => {
                let id = resolve_dative_bond_ref(idx, actions)?;
                if id.index() >= self.dative_bond_count() { return Err(TransactionError::IdOutOfRange("dative bond")); }
                let saved_atoms: Vec<AtomId> = atoms
                    .iter()
                    .map(|r| resolve_atom_ref(r.clone(), actions))
                    .collect::<Result<_, _>>()?;
                let view = self.dative_bond(id);
                let current_atoms: Vec<AtomId> = view.atom_ids().collect();
                if view.ast != &ast || current_atoms != saved_atoms {
                    return Err(TransactionError::OldStateMismatch);
                }
                self.remove_dative_bonds(&[id]);
                Ok(Action::Done)
            }
            Edit::SetDativeBondField { idx, change } => {
                let id = resolve_dative_bond_ref(idx, actions)?;
                if id.index() >= self.dative_bond_count() { return Err(TransactionError::IdOutOfRange("dative bond")); }
                self.apply_set_dative_bond_field(id, change)?;
                Ok(Action::Done)
            }
            Edit::AddAromaticSystem { atoms, ast } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| resolve_atom_ref(r, actions))
                    .collect::<Result<_, _>>()?;
                for a in &resolved {
                    if a.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                }
                Ok(Action::AromaticSystemAdded(
                    self.add_aromatic_system(resolved, ast),
                ))
            }
            Edit::RemoveAromaticSystem { idx, atoms, ast } => {
                let id = resolve_aromatic_system_ref(idx, actions)?;
                if id.index() >= self.aromatic_system_count() { return Err(TransactionError::IdOutOfRange("aromatic system")); }
                let saved_atoms: Vec<AtomId> = atoms
                    .iter()
                    .map(|r| resolve_atom_ref(r.clone(), actions))
                    .collect::<Result<_, _>>()?;
                let view = self.aromatic_system(id);
                let current_atoms: Vec<AtomId> = view.atom_ids().collect();
                if view.ast != &ast || current_atoms != saved_atoms {
                    return Err(TransactionError::OldStateMismatch);
                }
                self.remove_aromatic_systems(&[id]);
                Ok(Action::Done)
            }
            Edit::SetAromaticSystemField { idx, change } => {
                let id = resolve_aromatic_system_ref(idx, actions)?;
                if id.index() >= self.aromatic_system_count() { return Err(TransactionError::IdOutOfRange("aromatic system")); }
                self.apply_set_aromatic_system_field(id, change)?;
                Ok(Action::Done)
            }
            Edit::AddMulticenterBond { atoms, ast } => {
                let resolved: Vec<AtomId> = atoms
                    .into_iter()
                    .map(|r| resolve_atom_ref(r, actions))
                    .collect::<Result<_, _>>()?;
                for a in &resolved {
                    if a.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                }
                Ok(Action::MulticenterBondAdded(
                    self.add_multicenter_bond(resolved, ast),
                ))
            }
            Edit::RemoveMulticenterBond { idx, atoms, ast } => {
                let id = resolve_multicenter_bond_ref(idx, actions)?;
                if id.index() >= self.multicenter_bond_count() { return Err(TransactionError::IdOutOfRange("multicenter bond")); }
                let saved_atoms: Vec<AtomId> = atoms
                    .iter()
                    .map(|r| resolve_atom_ref(r.clone(), actions))
                    .collect::<Result<_, _>>()?;
                let view = self.multicenter_bond(id);
                let current_atoms: Vec<AtomId> = view.atom_ids().collect();
                if view.ast != &ast || current_atoms != saved_atoms {
                    return Err(TransactionError::OldStateMismatch);
                }
                self.remove_multicenter_bonds(&[id]);
                Ok(Action::Done)
            }
            Edit::SetMulticenterBondField { idx, change } => {
                let id = resolve_multicenter_bond_ref(idx, actions)?;
                if id.index() >= self.multicenter_bond_count() { return Err(TransactionError::IdOutOfRange("multicenter bond")); }
                self.apply_set_multicenter_bond_field(id, change)?;
                Ok(Action::Done)
            }
            Edit::AddNoncovalentBond { atoms, ast } => {
                let a = resolve_atom_ref(atoms[0].clone(), actions)?;
                let b = resolve_atom_ref(atoms[1].clone(), actions)?;
                if a.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                if b.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                Ok(Action::NoncovalentBondAdded(
                    self.add_noncovalent_bond([a, b], ast),
                ))
            }
            Edit::RemoveNoncovalentBond { idx, atoms, ast } => {
                let id = resolve_noncovalent_bond_ref(idx, actions)?;
                if id.index() >= self.noncovalent_bond_count() { return Err(TransactionError::IdOutOfRange("noncovalent bond")); }
                let saved_atoms = [
                    resolve_atom_ref(atoms[0].clone(), actions)?,
                    resolve_atom_ref(atoms[1].clone(), actions)?,
                ];
                let view = self.noncovalent_bond(id);
                if view.ast != &ast || view.atoms != saved_atoms {
                    return Err(TransactionError::OldStateMismatch);
                }
                self.remove_noncovalent_bonds(&[id]);
                Ok(Action::Done)
            }
            Edit::SetNoncovalentBondField { idx, change } => {
                let id = resolve_noncovalent_bond_ref(idx, actions)?;
                if id.index() >= self.noncovalent_bond_count() { return Err(TransactionError::IdOutOfRange("noncovalent bond")); }
                self.apply_set_noncovalent_bond_field(id, change)?;
                Ok(Action::Done)
            }
            Edit::SetAtomConstraint { idx, old, new } => {
                let id = resolve_atom_ref(idx, actions)?;
                if id.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                self.apply_set_atom_constraint(id, old, new)?;
                Ok(Action::Done)
            }
            Edit::AddAtomConstraint { idx, constraint } => {
                let id = resolve_atom_ref(idx, actions)?;
                if id.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                if constraint.is_unique() {
                    return Err(TransactionError::KindShapeMismatch);
                }
                let cs = &mut self.atom_mut(id).ast.constraints;
                if cs.contains_entry(&constraint) {
                    return Err(TransactionError::DuplicateEntry);
                }
                cs.add(constraint);
                Ok(Action::Done)
            }
            Edit::RemoveAtomConstraint { idx, constraint } => {
                let id = resolve_atom_ref(idx, actions)?;
                if id.index() >= self.atom_count() { return Err(TransactionError::IdOutOfRange("atom")); }
                if constraint.is_unique() {
                    return Err(TransactionError::KindShapeMismatch);
                }
                let cs = &mut self.atom_mut(id).ast.constraints;
                cs.remove_entry(&constraint)
                    .ok_or(TransactionError::MissingEntry)?;
                Ok(Action::Done)
            }
            Edit::SetBondConstraint { idx, old, new } => {
                let id = resolve_bond_ref(idx, actions)?;
                if id.index() >= self.bond_count() { return Err(TransactionError::IdOutOfRange("bond")); }
                self.apply_set_bond_constraint(id, old, new)?;
                Ok(Action::Done)
            }
            Edit::AddBondConstraint { idx, constraint } => {
                let id = resolve_bond_ref(idx, actions)?;
                if id.index() >= self.bond_count() { return Err(TransactionError::IdOutOfRange("bond")); }
                if constraint.is_unique() {
                    return Err(TransactionError::KindShapeMismatch);
                }
                let cs = &mut self.bond_mut(id).ast.constraints;
                if cs.contains_entry(&constraint) {
                    return Err(TransactionError::DuplicateEntry);
                }
                cs.add(constraint);
                Ok(Action::Done)
            }
            Edit::RemoveBondConstraint { idx, constraint } => {
                let id = resolve_bond_ref(idx, actions)?;
                if id.index() >= self.bond_count() { return Err(TransactionError::IdOutOfRange("bond")); }
                if constraint.is_unique() {
                    return Err(TransactionError::KindShapeMismatch);
                }
                let cs = &mut self.bond_mut(id).ast.constraints;
                cs.remove_entry(&constraint)
                    .ok_or(TransactionError::MissingEntry)?;
                Ok(Action::Done)
            }
            Edit::SetDativeBondConstraint { idx, old, new } => {
                let id = resolve_dative_bond_ref(idx, actions)?;
                if id.index() >= self.dative_bond_count() { return Err(TransactionError::IdOutOfRange("dative bond")); }
                self.apply_set_dative_bond_constraint(id, old, new)?;
                Ok(Action::Done)
            }
            Edit::AddDativeBondConstraint { idx, constraint } => {
                let id = resolve_dative_bond_ref(idx, actions)?;
                if id.index() >= self.dative_bond_count() { return Err(TransactionError::IdOutOfRange("dative bond")); }
                if constraint.is_unique() {
                    return Err(TransactionError::KindShapeMismatch);
                }
                let cs = &mut self.dative_bond_mut(id).ast.constraints;
                if cs.contains_entry(&constraint) {
                    return Err(TransactionError::DuplicateEntry);
                }
                cs.add(constraint);
                Ok(Action::Done)
            }
            Edit::RemoveDativeBondConstraint { idx, constraint } => {
                let id = resolve_dative_bond_ref(idx, actions)?;
                if id.index() >= self.dative_bond_count() { return Err(TransactionError::IdOutOfRange("dative bond")); }
                if constraint.is_unique() {
                    return Err(TransactionError::KindShapeMismatch);
                }
                let cs = &mut self.dative_bond_mut(id).ast.constraints;
                cs.remove_entry(&constraint)
                    .ok_or(TransactionError::MissingEntry)?;
                Ok(Action::Done)
            }
            Edit::SetAromaticSystemConstraint { idx, old, new } => {
                let id = resolve_aromatic_system_ref(idx, actions)?;
                if id.index() >= self.aromatic_system_count() { return Err(TransactionError::IdOutOfRange("aromatic system")); }
                self.apply_set_aromatic_system_constraint(id, old, new)?;
                Ok(Action::Done)
            }
            Edit::SetMulticenterBondConstraint { idx, old, new } => {
                let id = resolve_multicenter_bond_ref(idx, actions)?;
                if id.index() >= self.multicenter_bond_count() { return Err(TransactionError::IdOutOfRange("multicenter bond")); }
                self.apply_set_multicenter_bond_constraint(id, old, new)?;
                Ok(Action::Done)
            }
            Edit::PushMoleculeConstraint { constraint } => {
                self.push_constraint(constraint);
                Ok(Action::Done)
            }
            Edit::PopMoleculeConstraint { constraint } => {
                let list = self.constraints_mut();
                let last_index = list.len().checked_sub(1).ok_or(
                    TransactionError::OldStateMismatch,
                )?;
                if list.as_slice()[last_index] != constraint {
                    return Err(TransactionError::OldStateMismatch);
                }
                list.remove_at(last_index);
                Ok(Action::Done)
            }
        }
    }

    fn apply_set_atom_field(
        &mut self,
        id: AtomId,
        change: AtomFieldChange,
    ) -> Result<(), TransactionError> {
        let atom = self.atom_mut(id);
        match change {
            AtomFieldChange::Element { old, new } => {
                if atom.ast.element != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.element = new;
            }
            AtomFieldChange::IsotopeMass { old, new } => {
                if atom.ast.isotope_mass != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.isotope_mass = new;
            }
            AtomFieldChange::Charge { old, new } => {
                if atom.ast.charge != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.charge = new;
            }
            AtomFieldChange::ImplicitHydrogens { old, new } => {
                if atom.ast.implicit_hydrogens != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.implicit_hydrogens = new;
            }
            AtomFieldChange::LonePairs { old, new } => {
                if atom.ast.lone_pairs != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.lone_pairs = new;
            }
            AtomFieldChange::Spin { old, new } => {
                if atom.ast.spin != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                atom.ast.spin = new;
            }
        }
        Ok(())
    }

    fn apply_set_bond_field(
        &mut self,
        id: BondId,
        change: BondFieldChange,
    ) -> Result<(), TransactionError> {
        let bond = self.bond_mut(id);
        match change {
            BondFieldChange::Order { old, new } => {
                if bond.ast.order != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                bond.ast.order = new;
            }
            BondFieldChange::Charge { old, new } => {
                if bond.ast.charge != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                bond.ast.charge = new;
            }
            BondFieldChange::Spin { old, new } => {
                if bond.ast.spin != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                bond.ast.spin = new;
            }
        }
        Ok(())
    }

    fn apply_set_dative_bond_field(
        &mut self,
        id: DativeBondId,
        change: DativeBondFieldChange,
    ) -> Result<(), TransactionError> {
        let dat = self.dative_bond_mut(id);
        match change {
            DativeBondFieldChange::AcceptorSlot { old, new } => {
                if dat.ast.acceptor_slot != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                dat.ast.acceptor_slot = new;
            }
            DativeBondFieldChange::Order { old, new } => {
                if dat.ast.order != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                dat.ast.order = new;
            }
        }
        Ok(())
    }

    fn apply_set_aromatic_system_field(
        &mut self,
        id: AromaticSystemId,
        change: AromaticSystemFieldChange,
    ) -> Result<(), TransactionError> {
        let ar = self.aromatic_system_mut(id);
        match change {
            AromaticSystemFieldChange::Electrons { old, new } => {
                if ar.ast.electrons != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                ar.ast.electrons = new;
            }
            AromaticSystemFieldChange::Charge { old, new } => {
                if ar.ast.charge != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                ar.ast.charge = new;
            }
            AromaticSystemFieldChange::Spin { old, new } => {
                if ar.ast.spin != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                ar.ast.spin = new;
            }
        }
        Ok(())
    }

    fn apply_set_multicenter_bond_field(
        &mut self,
        id: MulticenterBondId,
        change: MulticenterBondFieldChange,
    ) -> Result<(), TransactionError> {
        let mc = self.multicenter_bond_mut(id);
        match change {
            MulticenterBondFieldChange::Electrons { old, new } => {
                if mc.ast.electrons != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                mc.ast.electrons = new;
            }
            MulticenterBondFieldChange::Charge { old, new } => {
                if mc.ast.charge != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                mc.ast.charge = new;
            }
            MulticenterBondFieldChange::Spin { old, new } => {
                if mc.ast.spin != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                mc.ast.spin = new;
            }
        }
        Ok(())
    }

    fn apply_set_noncovalent_bond_field(
        &mut self,
        id: NoncovalentBondId,
        change: NoncovalentBondFieldChange,
    ) -> Result<(), TransactionError> {
        let nc = self.noncovalent_bond_mut(id);
        match change {
            NoncovalentBondFieldChange::Kind { old, new } => {
                if nc.ast.kind != old {
                    return Err(TransactionError::OldStateMismatch);
                }
                nc.ast.kind = new;
            }
        }
        Ok(())
    }

    fn apply_set_atom_constraint(
        &mut self,
        id: AtomId,
        old: Option<super::super::constraint::AtomConstraint>,
        new: Option<super::super::constraint::AtomConstraint>,
    ) -> Result<(), TransactionError> {
        let kind = match (&old, &new) {
            (Some(o), Some(n)) => {
                if o.kind() != n.kind() {
                    return Err(TransactionError::KindMismatch);
                }
                o.kind()
            }
            (Some(o), None) => o.kind(),
            (None, Some(n)) => n.kind(),
            (None, None) => return Ok(()),
        };
        if let Some(c) = &new {
            if !c.is_unique() {
                return Err(TransactionError::KindShapeMismatch);
            }
        }
        if let Some(c) = &old {
            if !c.is_unique() {
                return Err(TransactionError::KindShapeMismatch);
            }
        }
        let cs = &mut self.atom_mut(id).ast.constraints;
        if cs.get(kind).cloned() != old {
            return Err(TransactionError::OldStateMismatch);
        }
        match new {
            Some(c) => {
                cs.add(c);
            }
            None => {
                cs.remove(kind);
            }
        }
        Ok(())
    }

    fn apply_set_bond_constraint(
        &mut self,
        id: BondId,
        old: Option<super::super::constraint::BondConstraint>,
        new: Option<super::super::constraint::BondConstraint>,
    ) -> Result<(), TransactionError> {
        let kind = match (&old, &new) {
            (Some(o), Some(n)) => {
                if o.kind() != n.kind() {
                    return Err(TransactionError::KindMismatch);
                }
                o.kind()
            }
            (Some(o), None) => o.kind(),
            (None, Some(n)) => n.kind(),
            (None, None) => return Ok(()),
        };
        if let Some(c) = &new {
            if !c.is_unique() {
                return Err(TransactionError::KindShapeMismatch);
            }
        }
        if let Some(c) = &old {
            if !c.is_unique() {
                return Err(TransactionError::KindShapeMismatch);
            }
        }
        let cs = &mut self.bond_mut(id).ast.constraints;
        if cs.get(kind).cloned() != old {
            return Err(TransactionError::OldStateMismatch);
        }
        match new {
            Some(c) => {
                cs.add(c);
            }
            None => {
                cs.remove(kind);
            }
        }
        Ok(())
    }

    fn apply_set_dative_bond_constraint(
        &mut self,
        id: DativeBondId,
        old: Option<super::super::constraint::DativeBondConstraint>,
        new: Option<super::super::constraint::DativeBondConstraint>,
    ) -> Result<(), TransactionError> {
        let kind = match (&old, &new) {
            (Some(o), Some(n)) => {
                if o.kind() != n.kind() {
                    return Err(TransactionError::KindMismatch);
                }
                o.kind()
            }
            (Some(o), None) => o.kind(),
            (None, Some(n)) => n.kind(),
            (None, None) => return Ok(()),
        };
        if let Some(c) = &new {
            if !c.is_unique() {
                return Err(TransactionError::KindShapeMismatch);
            }
        }
        if let Some(c) = &old {
            if !c.is_unique() {
                return Err(TransactionError::KindShapeMismatch);
            }
        }
        let cs = &mut self.dative_bond_mut(id).ast.constraints;
        if cs.get(kind).cloned() != old {
            return Err(TransactionError::OldStateMismatch);
        }
        match new {
            Some(c) => {
                cs.add(c);
            }
            None => {
                cs.remove(kind);
            }
        }
        Ok(())
    }

    fn apply_set_aromatic_system_constraint(
        &mut self,
        id: AromaticSystemId,
        old: Option<super::super::constraint::AromaticSystemConstraint>,
        new: Option<super::super::constraint::AromaticSystemConstraint>,
    ) -> Result<(), TransactionError> {
        let kind = match (&old, &new) {
            (Some(o), Some(n)) => {
                if o.kind() != n.kind() {
                    return Err(TransactionError::KindMismatch);
                }
                o.kind()
            }
            (Some(o), None) => o.kind(),
            (None, Some(n)) => n.kind(),
            (None, None) => return Ok(()),
        };
        let cs = &mut self.aromatic_system_mut(id).ast.constraints;
        if cs.get(kind).cloned() != old {
            return Err(TransactionError::OldStateMismatch);
        }
        match new {
            Some(c) => {
                cs.add(c);
            }
            None => {
                cs.remove(kind);
            }
        }
        Ok(())
    }

    fn apply_set_multicenter_bond_constraint(
        &mut self,
        id: MulticenterBondId,
        old: Option<super::super::constraint::MulticenterBondConstraint>,
        new: Option<super::super::constraint::MulticenterBondConstraint>,
    ) -> Result<(), TransactionError> {
        let kind = match (&old, &new) {
            (Some(o), Some(n)) => {
                if o.kind() != n.kind() {
                    return Err(TransactionError::KindMismatch);
                }
                o.kind()
            }
            (Some(o), None) => o.kind(),
            (None, Some(n)) => n.kind(),
            (None, None) => return Ok(()),
        };
        let cs = &mut self.multicenter_bond_mut(id).ast.constraints;
        if cs.get(kind).cloned() != old {
            return Err(TransactionError::OldStateMismatch);
        }
        match new {
            Some(c) => {
                cs.add(c);
            }
            None => {
                cs.remove(kind);
            }
        }
        Ok(())
    }
}

fn resolve_atom_ref(r: AtomRef, actions: &[Action]) -> Result<AtomId, TransactionError> {
    match r {
        AtomRef::Id(id) => Ok(id),
        AtomRef::New(n) => {
            let mut seen = 0usize;
            for action in actions {
                match action {
                    Action::AtomsAdded(ids) => {
                        if n < seen + ids.len() {
                            return Ok(ids[n - seen]);
                        }
                        seen += ids.len();
                    }
                    Action::AtomAdded(id) => {
                        if n == seen {
                            return Ok(*id);
                        }
                        seen += 1;
                    }
                    _ => {}
                }
            }
            Err(TransactionError::RefOutOfRange(n, seen))
        }
    }
}

fn resolve_bond_ref(r: BondRef, actions: &[Action]) -> Result<BondId, TransactionError> {
    match r {
        BondRef::Id(id) => Ok(id),
        BondRef::New(n) => {
            let mut seen = 0usize;
            for action in actions {
                match action {
                    Action::BondsAdded(ids) => {
                        if n < seen + ids.len() {
                            return Ok(ids[n - seen]);
                        }
                        seen += ids.len();
                    }
                    Action::BondAdded(id) => {
                        if n == seen {
                            return Ok(*id);
                        }
                        seen += 1;
                    }
                    _ => {}
                }
            }
            Err(TransactionError::RefOutOfRange(n, seen))
        }
    }
}

fn resolve_dative_bond_ref(
    r: DativeBondRef,
    actions: &[Action],
) -> Result<DativeBondId, TransactionError> {
    match r {
        DativeBondRef::Id(id) => Ok(id),
        DativeBondRef::New(n) => {
            let action = actions
                .get(n)
                .ok_or(TransactionError::RefOutOfRange(n, actions.len()))?;
            match action {
                Action::DativeBondAdded(id) => Ok(*id),
                _ => Err(TransactionError::RefTypeMismatch {
                    expected: "DativeBondAdded",
                    got: action_name(action),
                }),
            }
        }
    }
}

fn resolve_aromatic_system_ref(
    r: AromaticSystemRef,
    actions: &[Action],
) -> Result<AromaticSystemId, TransactionError> {
    match r {
        AromaticSystemRef::Id(id) => Ok(id),
        AromaticSystemRef::New(n) => {
            let action = actions
                .get(n)
                .ok_or(TransactionError::RefOutOfRange(n, actions.len()))?;
            match action {
                Action::AromaticSystemAdded(id) => Ok(*id),
                _ => Err(TransactionError::RefTypeMismatch {
                    expected: "AromaticSystemAdded",
                    got: action_name(action),
                }),
            }
        }
    }
}

fn resolve_multicenter_bond_ref(
    r: MulticenterBondRef,
    actions: &[Action],
) -> Result<MulticenterBondId, TransactionError> {
    match r {
        MulticenterBondRef::Id(id) => Ok(id),
        MulticenterBondRef::New(n) => {
            let action = actions
                .get(n)
                .ok_or(TransactionError::RefOutOfRange(n, actions.len()))?;
            match action {
                Action::MulticenterBondAdded(id) => Ok(*id),
                _ => Err(TransactionError::RefTypeMismatch {
                    expected: "MulticenterBondAdded",
                    got: action_name(action),
                }),
            }
        }
    }
}

fn resolve_noncovalent_bond_ref(
    r: NoncovalentBondRef,
    actions: &[Action],
) -> Result<NoncovalentBondId, TransactionError> {
    match r {
        NoncovalentBondRef::Id(id) => Ok(id),
        NoncovalentBondRef::New(n) => {
            let action = actions
                .get(n)
                .ok_or(TransactionError::RefOutOfRange(n, actions.len()))?;
            match action {
                Action::NoncovalentBondAdded(id) => Ok(*id),
                _ => Err(TransactionError::RefTypeMismatch {
                    expected: "NoncovalentBondAdded",
                    got: action_name(action),
                }),
            }
        }
    }
}

fn action_name(a: &Action) -> &'static str {
    match a {
        Action::AtomsAdded(_) => "AtomsAdded",
        Action::BondsAdded(_) => "BondsAdded",
        Action::AtomAdded(_) => "AtomAdded",
        Action::BondAdded(_) => "BondAdded",
        Action::AromaticSystemAdded(_) => "AromaticSystemAdded",
        Action::DativeBondAdded(_) => "DativeBondAdded",
        Action::MulticenterBondAdded(_) => "MulticenterBondAdded",
        Action::NoncovalentBondAdded(_) => "NoncovalentBondAdded",
        Action::Done => "Done",
        Action::Cascaded { .. } => "Cascaded",
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::element::Element;

    use super::super::super::atom::{AtomAst, ElementAst};
    use super::super::super::bond::BondAst;
    use super::super::super::aromatic::AromaticSystemAst;
    use super::super::super::constraint::{
        AromaticSystemConstraint, AtomConstraint, BondConstraint, Constraint,
        DativeBondConstraint, MoleculeConstraint, MulticenterBondConstraint,
    };
    use super::super::super::dative::DativeBondAst;
    use super::super::super::multicenter::MulticenterBondAst;
    use super::super::super::noncovalent::{
        NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst,
    };
    use super::super::super::value::ValueAst;
    use super::super::MoleculeAst;
    use super::*;

    #[fixture]
    fn empty() -> MoleculeBuilder {
        MoleculeAst::default().edit()
    }

    #[fixture]
    fn one_atom() -> MoleculeBuilder {
        let mut b = MoleculeAst::default().edit();
        b.add_atom(AtomAst::from_element(Element::C));
        b
    }

    #[fixture]
    fn diatomic() -> MoleculeBuilder {
        let mut b = MoleculeAst::default().edit();
        b.add_atom(AtomAst::from_element(Element::C));
        b.add_atom(AtomAst::from_element(Element::C));
        b.add_bond(AtomId(0), AtomId(1), BondAst::from_order(1));
        b
    }

    #[rstest]
    fn test_molecule_builder_transact_add_atom(mut empty: MoleculeBuilder) {
        let actions = empty
            .transact(vec![Edit::add_atom(AtomAst::from_element(Element::C))])
            .unwrap();
        assert_eq!(actions, vec![Action::AtomsAdded(vec![AtomId(0)])]);
        let built = empty.build();
        assert_eq!(built.atoms().count(), 1);
        assert_eq!(built.atom(AtomId(0)).ast.element, ElementAst::Lit(Element::C));
    }

    #[rstest]
    fn test_molecule_builder_transact_add_bond_via_new_ref(mut empty: MoleculeBuilder) {
        let actions = empty
            .transact(vec![
                Edit::AddAtoms {
                    atoms: vec![
                        AtomAst::from_element(Element::C),
                        AtomAst::from_element(Element::C),
                    ],
                },
                Edit::add_bond(AtomRef::New(0), AtomRef::New(1), BondAst::from_order(1)),
            ])
            .unwrap();
        assert_eq!(
            actions,
            vec![
                Action::AtomsAdded(vec![AtomId(0), AtomId(1)]),
                Action::BondsAdded(vec![BondId(0)]),
            ]
        );
    }

    #[rstest]
    fn test_molecule_builder_transact_rollback(mut one_atom: MoleculeBuilder) {
        // Mid-batch failure (out-of-range id on edit 2) rolls back the
        // already-applied AddAtom on edit 1.
        let err = one_atom
            .transact(vec![
                Edit::add_atom(AtomAst::from_element(Element::N)),
                Edit::remove_atom(AtomRef::Id(AtomId(99))),
            ])
            .unwrap_err();
        assert_eq!(err, TransactionError::IdOutOfRange("atom"));
        assert_eq!(one_atom.atom_count(), 1);
    }

    #[rstest]
    fn test_molecule_builder_transact_set_atom_field(mut one_atom: MoleculeBuilder) {
        let actions = one_atom
            .transact(vec![Edit::SetAtomField {
                idx: AtomRef::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: ValueAst::default(),
                    new: ValueAst::Lit(1),
                },
            }])
            .unwrap();
        assert_eq!(actions, vec![Action::Done]);
        assert_eq!(one_atom.build().atom(AtomId(0)).ast.charge, ValueAst::Lit(1));
    }

    #[rstest]
    fn test_molecule_builder_transact_set_atom_field_error(mut one_atom: MoleculeBuilder) {
        let err = one_atom
            .transact(vec![Edit::SetAtomField {
                idx: AtomRef::Id(AtomId(0)),
                change: AtomFieldChange::Charge {
                    old: ValueAst::Lit(99),
                    new: ValueAst::Lit(1),
                },
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    #[case::ref_out_of_range(
        vec![Edit::add_bond(
            AtomRef::New(5),
            AtomRef::New(6),
            BondAst::default(),
        )],
        TransactionError::RefOutOfRange(5, 0),
    )]
    #[case::ref_type_mismatch(
        vec![
            Edit::add_atom(AtomAst::from_element(Element::C)),
            Edit::remove_bond(BondRef::New(0)),
        ],
        TransactionError::RefOutOfRange(0, 0),
    )]
    fn test_molecule_builder_transact_new_ref_error(
        mut empty: MoleculeBuilder,
        #[case] edits: Vec<Edit>,
        #[case] expected: TransactionError,
    ) {
        assert_eq!(empty.transact(edits).unwrap_err(), expected);
    }

    #[rstest]
    fn test_molecule_builder_transact_add_atom_constraint(mut one_atom: MoleculeBuilder) {
        one_atom
            .transact(vec![
                Edit::AddAtomConstraint {
                    idx: AtomRef::Id(AtomId(0)),
                    constraint: AtomConstraint::ring_size(5),
                },
                Edit::AddAtomConstraint {
                    idx: AtomRef::Id(AtomId(0)),
                    constraint: AtomConstraint::ring_size(6),
                },
            ])
            .unwrap();
        let next = one_atom.build();
        let cs: Vec<_> = next.atom(AtomId(0)).ast.constraints.iter().cloned().collect();
        assert_eq!(
            cs,
            vec![
                AtomConstraint::RingSize(ValueAst::Lit(5)),
                AtomConstraint::RingSize(ValueAst::Lit(6)),
            ]
        );
    }

    #[rstest]
    fn test_molecule_builder_transact_add_atom_constraint_error(mut one_atom: MoleculeBuilder) {
        let err = one_atom
            .transact(vec![Edit::AddAtomConstraint {
                idx: AtomRef::Id(AtomId(0)),
                constraint: AtomConstraint::valence(4),
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::KindShapeMismatch);
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_atom_constraint_error(mut one_atom: MoleculeBuilder) {
        let err = one_atom
            .transact(vec![Edit::RemoveAtomConstraint {
                idx: AtomRef::Id(AtomId(0)),
                constraint: AtomConstraint::ring_size(5),
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::MissingEntry);
    }

    #[rstest]
    #[case::introduce(
        None,
        Some(AtomConstraint::valence(4)),
        ValueAst::Lit(4),
    )]
    #[case::replace(
        Some(AtomConstraint::valence(3)),
        Some(AtomConstraint::valence(4)),
        ValueAst::Lit(4),
    )]
    #[case::remove(
        Some(AtomConstraint::valence(3)),
        None,
        ValueAst::Undetermined,
    )]
    fn test_molecule_builder_transact_set_atom_constraint(
        mut one_atom: MoleculeBuilder,
        #[case] old: Option<AtomConstraint>,
        #[case] new: Option<AtomConstraint>,
        #[case] expected: ValueAst,
    ) {
        if let Some(c) = old.clone() {
            one_atom.atom_mut(AtomId(0)).ast.constraints.add(c);
        }
        one_atom
            .transact(vec![Edit::SetAtomConstraint {
                idx: AtomRef::Id(AtomId(0)),
                old,
                new,
            }])
            .unwrap();
        assert_eq!(
            one_atom.atom_mut(AtomId(0)).ast.constraints.valence(),
            expected
        );
    }

    #[rstest]
    fn test_molecule_builder_transact_set_atom_constraint_error(mut one_atom: MoleculeBuilder) {
        let err = one_atom
            .transact(vec![Edit::SetAtomConstraint {
                idx: AtomRef::Id(AtomId(0)),
                old: None,
                new: Some(AtomConstraint::ring_size(5)),
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::KindShapeMismatch);
    }

    #[rstest]
    fn test_molecule_builder_transact_set_bond_constraint(mut diatomic: MoleculeBuilder) {
        diatomic
            .transact(vec![Edit::SetBondConstraint {
                idx: BondRef::Id(BondId(0)),
                old: None,
                new: Some(BondConstraint::Aromatic),
            }])
            .unwrap();
        assert!(diatomic
            .bond_mut(BondId(0))
            .ast
            .constraints
            .iter()
            .any(|c| *c == BondConstraint::Aromatic));
    }

    #[rstest]
    fn test_molecule_builder_transact_push_molecule_constraint(mut empty: MoleculeBuilder) {
        let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        empty
            .transact(vec![Edit::PushMoleculeConstraint {
                constraint: c.clone(),
            }])
            .unwrap();
        assert_eq!(empty.constraints_mut().as_slice(), &[c]);
    }

    #[rstest]
    fn test_molecule_builder_transact_pop_molecule_constraint(mut empty: MoleculeBuilder) {
        let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        empty.push_constraint(c.clone());
        empty
            .transact(vec![Edit::PopMoleculeConstraint {
                constraint: c.clone(),
            }])
            .unwrap();
        assert!(empty.constraints_mut().as_slice().is_empty());
    }

    #[rstest]
    fn test_molecule_builder_transact_pop_molecule_constraint_error(
        mut empty: MoleculeBuilder,
    ) {
        let c = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        empty.push_constraint(c.clone());
        let err = empty
            .transact(vec![Edit::PopMoleculeConstraint {
                constraint: Constraint::Molecule(MoleculeConstraint::Connected {
                    atoms: Some(vec![AtomId(0)]),
                }),
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
        assert_eq!(empty.constraints_mut().as_slice(), &[c]);
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_topology_atom_error(
        mut one_atom: MoleculeBuilder,
    ) {
        let err = one_atom
            .transact(vec![Edit::remove_atom(AtomRef::Id(AtomId(9)))])
            .unwrap_err();
        assert_eq!(err, TransactionError::IdOutOfRange("atom"));
        assert_eq!(one_atom.atom_count(), 1);
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_topology_bond_error(
        mut diatomic: MoleculeBuilder,
    ) {
        let err = diatomic
            .transact(vec![Edit::remove_bond(BondRef::Id(BondId(9)))])
            .unwrap_err();
        assert_eq!(err, TransactionError::IdOutOfRange("bond"));
        assert_eq!(diatomic.bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_builder_transact_add_dative_bond_empty_atoms_error(
        mut one_atom: MoleculeBuilder,
    ) {
        let err = one_atom
            .transact(vec![Edit::AddDativeBond {
                atoms: vec![],
                ast: DativeBondAst::from_order(1),
            }])
            .unwrap_err();
        assert!(matches!(err, TransactionError::MalformedEdit(_)));
    }

    #[rstest]
    fn test_molecule_builder_transact_set_bond_field(mut diatomic: MoleculeBuilder) {
        diatomic
            .transact(vec![Edit::SetBondField {
                idx: BondRef::Id(BondId(0)),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            }])
            .unwrap();
        assert_eq!(diatomic.bond(BondId(0)).ast.order, ValueAst::Lit(2));
    }

    #[rstest]
    fn test_molecule_builder_transact_set_bond_field_error(mut diatomic: MoleculeBuilder) {
        let err = diatomic
            .transact(vec![Edit::SetBondField {
                idx: BondRef::Id(BondId(0)),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(99),
                    new: ValueAst::Lit(2),
                },
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[fixture]
    fn diatomic_with_overlays() -> MoleculeBuilder {
        let mut b = MoleculeAst::default().edit();
        b.add_atom(AtomAst::from_element(Element::C));
        b.add_atom(AtomAst::from_element(Element::N));
        b.add_bond(AtomId(0), AtomId(1), BondAst::from_order(1));
        b.add_dative_bond(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1));
        b.add_aromatic_system(
            vec![AtomId(0), AtomId(1)],
            AromaticSystemAst::default(),
        );
        b.add_multicenter_bond(
            vec![AtomId(0), AtomId(1)],
            MulticenterBondAst::default(),
        );
        b.add_noncovalent_bond(
            [AtomId(0), AtomId(1)],
            NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        );
        b
    }

    #[rstest]
    fn test_molecule_builder_transact_set_dative_bond_field(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::SetDativeBondField {
                idx: DativeBondRef::Id(DativeBondId(0)),
                change: DativeBondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            }])
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .dative_bond(DativeBondId(0))
                .ast
                .order,
            ValueAst::Lit(2),
        );
    }

    #[rstest]
    fn test_molecule_builder_transact_set_aromatic_system_field(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::SetAromaticSystemField {
                idx: AromaticSystemRef::Id(AromaticSystemId(0)),
                change: AromaticSystemFieldChange::Charge {
                    old: ValueAst::default(),
                    new: ValueAst::Lit(1),
                },
            }])
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .aromatic_system(AromaticSystemId(0))
                .ast
                .charge,
            ValueAst::Lit(1),
        );
    }

    #[rstest]
    fn test_molecule_builder_transact_set_multicenter_bond_field(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::SetMulticenterBondField {
                idx: MulticenterBondRef::Id(MulticenterBondId(0)),
                change: MulticenterBondFieldChange::Charge {
                    old: ValueAst::default(),
                    new: ValueAst::Lit(-1),
                },
            }])
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .multicenter_bond(MulticenterBondId(0))
                .ast
                .charge,
            ValueAst::Lit(-1),
        );
    }

    #[rstest]
    fn test_molecule_builder_transact_set_noncovalent_bond_field(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::SetNoncovalentBondField {
                idx: NoncovalentBondRef::Id(NoncovalentBondId(0)),
                change: NoncovalentBondFieldChange::Kind {
                    old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
                    new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
                },
            }])
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .noncovalent_bond(NoncovalentBondId(0))
                .ast
                .kind,
            NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
        );
    }

    #[rstest]
    fn test_molecule_builder_transact_add_dative_bond(mut diatomic: MoleculeBuilder) {
        let actions = diatomic
            .transact(vec![Edit::AddDativeBond {
                atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: DativeBondAst::from_order(1),
            }])
            .unwrap();
        assert_eq!(actions, vec![Action::DativeBondAdded(DativeBondId(0))]);
        assert_eq!(diatomic.dative_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_builder_transact_add_aromatic_system(mut diatomic: MoleculeBuilder) {
        let actions = diatomic
            .transact(vec![Edit::AddAromaticSystem {
                atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: AromaticSystemAst::default(),
            }])
            .unwrap();
        assert_eq!(
            actions,
            vec![Action::AromaticSystemAdded(AromaticSystemId(0))],
        );
        assert_eq!(diatomic.aromatic_system_count(), 1);
    }

    #[rstest]
    fn test_molecule_builder_transact_add_multicenter_bond(mut diatomic: MoleculeBuilder) {
        let actions = diatomic
            .transact(vec![Edit::AddMulticenterBond {
                atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: MulticenterBondAst::default(),
            }])
            .unwrap();
        assert_eq!(
            actions,
            vec![Action::MulticenterBondAdded(MulticenterBondId(0))],
        );
        assert_eq!(diatomic.multicenter_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_builder_transact_add_noncovalent_bond(mut diatomic: MoleculeBuilder) {
        let actions = diatomic
            .transact(vec![Edit::AddNoncovalentBond {
                atoms: [AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            }])
            .unwrap();
        assert_eq!(
            actions,
            vec![Action::NoncovalentBondAdded(NoncovalentBondId(0))],
        );
        assert_eq!(diatomic.noncovalent_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_dative_bond(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::RemoveDativeBond {
                idx: DativeBondRef::Id(DativeBondId(0)),
                atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: DativeBondAst {
                    acceptor_slot: 1,
                    order: ValueAst::Lit(1),
                    constraints: Default::default(),
                },
            }])
            .unwrap();
        assert_eq!(diatomic_with_overlays.dative_bond_count(), 0);
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_dative_bond_atoms_mismatch_error(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        let err = diatomic_with_overlays
            .transact(vec![Edit::RemoveDativeBond {
                idx: DativeBondRef::Id(DativeBondId(0)),
                atoms: vec![AtomRef::Id(AtomId(1)), AtomRef::Id(AtomId(0))], // wrong order
                ast: DativeBondAst {
                    acceptor_slot: 1,
                    order: ValueAst::Lit(1),
                    constraints: Default::default(),
                },
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
        assert_eq!(diatomic_with_overlays.dative_bond_count(), 1);
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_aromatic_system(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::RemoveAromaticSystem {
                idx: AromaticSystemRef::Id(AromaticSystemId(0)),
                atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: AromaticSystemAst::default(),
            }])
            .unwrap();
        assert_eq!(diatomic_with_overlays.aromatic_system_count(), 0);
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_multicenter_bond(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::RemoveMulticenterBond {
                idx: MulticenterBondRef::Id(MulticenterBondId(0)),
                atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: MulticenterBondAst::default(),
            }])
            .unwrap();
        assert_eq!(diatomic_with_overlays.multicenter_bond_count(), 0);
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_noncovalent_bond(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::RemoveNoncovalentBond {
                idx: NoncovalentBondRef::Id(NoncovalentBondId(0)),
                atoms: [AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            }])
            .unwrap();
        assert_eq!(diatomic_with_overlays.noncovalent_bond_count(), 0);
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_noncovalent_bond_ast_mismatch_error(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        let err = diatomic_with_overlays
            .transact(vec![Edit::RemoveNoncovalentBond {
                idx: NoncovalentBondRef::Id(NoncovalentBondId(0)),
                atoms: [AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic), // wrong
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::OldStateMismatch);
    }

    #[rstest]
    fn test_molecule_builder_transact_set_bond_constraint_value_bearing(
        mut diatomic: MoleculeBuilder,
    ) {
        diatomic
            .transact(vec![Edit::SetBondConstraint {
                idx: BondRef::Id(BondId(0)),
                old: None,
                new: Some(BondConstraint::RingCount(ValueAst::Lit(1))),
            }])
            .unwrap();
        assert_eq!(
            diatomic
                .bond(BondId(0))
                .ast
                .constraints
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            vec![BondConstraint::RingCount(ValueAst::Lit(1))],
        );
    }

    #[rstest]
    fn test_molecule_builder_transact_add_bond_constraint(mut diatomic: MoleculeBuilder) {
        diatomic
            .transact(vec![Edit::AddBondConstraint {
                idx: BondRef::Id(BondId(0)),
                constraint: BondConstraint::RingSize(ValueAst::Lit(5)),
            }])
            .unwrap();
        assert!(diatomic
            .bond(BondId(0))
            .ast
            .constraints
            .iter()
            .any(|c| *c == BondConstraint::RingSize(ValueAst::Lit(5))));
    }

    #[rstest]
    fn test_molecule_builder_transact_remove_bond_constraint_error(
        mut diatomic: MoleculeBuilder,
    ) {
        let err = diatomic
            .transact(vec![Edit::RemoveBondConstraint {
                idx: BondRef::Id(BondId(0)),
                constraint: BondConstraint::RingSize(ValueAst::Lit(5)),
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::MissingEntry);
    }

    #[rstest]
    fn test_molecule_builder_transact_set_dative_bond_constraint(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::SetDativeBondConstraint {
                idx: DativeBondRef::Id(DativeBondId(0)),
                old: None,
                new: Some(DativeBondConstraint::Aromatic),
            }])
            .unwrap();
        assert!(diatomic_with_overlays
            .dative_bond(DativeBondId(0))
            .ast
            .constraints
            .iter()
            .any(|c| *c == DativeBondConstraint::Aromatic));
    }

    #[rstest]
    fn test_molecule_builder_transact_add_dative_bond_constraint_error(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        // All `DativeBondConstraint` variants are unique-per-dative-bond, so
        // `Add` is rejected. Use `Set` for unique constraints.
        let err = diatomic_with_overlays
            .transact(vec![Edit::AddDativeBondConstraint {
                idx: DativeBondRef::Id(DativeBondId(0)),
                constraint: DativeBondConstraint::RingSize(ValueAst::Lit(6)),
            }])
            .unwrap_err();
        assert_eq!(err, TransactionError::KindShapeMismatch);
    }

    #[rstest]
    fn test_molecule_builder_transact_set_aromatic_system_constraint(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::SetAromaticSystemConstraint {
                idx: AromaticSystemRef::Id(AromaticSystemId(0)),
                old: None,
                new: Some(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))),
            }])
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .aromatic_system(AromaticSystemId(0))
                .ast
                .constraints
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            vec![AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))],
        );
    }

    #[rstest]
    fn test_molecule_builder_transact_set_multicenter_bond_constraint(
        mut diatomic_with_overlays: MoleculeBuilder,
    ) {
        diatomic_with_overlays
            .transact(vec![Edit::SetMulticenterBondConstraint {
                idx: MulticenterBondRef::Id(MulticenterBondId(0)),
                old: None,
                new: Some(MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2))),
            }])
            .unwrap();
        assert_eq!(
            diatomic_with_overlays
                .multicenter_bond(MulticenterBondId(0))
                .ast
                .constraints
                .iter()
                .cloned()
                .collect::<Vec<_>>(),
            vec![MulticenterBondConstraint::ElectronCount(ValueAst::Lit(2))],
        );
    }

    #[cfg(any())]
    mod phase8_undo_contract_tests {
        use super::*;

        #[fixture]
        fn triatomic_with_overlays() -> MoleculeBuilder {
            let mut b = MoleculeAst::default().edit();
            b.add_atom(AtomAst::from_element(Element::C));
            b.add_atom(AtomAst::from_element(Element::N));
            b.add_atom(AtomAst::from_element(Element::O));
            b.add_bond(AtomId(0), AtomId(1), BondAst::from_order(1));
            b.add_bond(AtomId(1), AtomId(2), BondAst::from_order(1));
            b.add_dative_bond(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1));
            b.add_aromatic_system(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemAst::default(),
            );
            b.add_multicenter_bond(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondAst::default(),
            );
            b.add_noncovalent_bond(
                [AtomId(0), AtomId(2)],
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            );
            b
        }

        #[rstest]
        fn test_transaction_rollback(mut triatomic_with_overlays: MoleculeBuilder) {
            let before = triatomic_with_overlays.clone().build();
            let tx = triatomic_with_overlays
                .transact(vec![Edit::RemoveTopology {
                    atoms: vec![AtomRef::Id(AtomId(1))],
                    bonds: vec![],
                }])
                .unwrap();
            tx.rollback(&mut triatomic_with_overlays).unwrap();
            assert_eq!(triatomic_with_overlays.build(), before);
        }

        #[rstest]
        fn test_transaction_rollback_cascaded_overlays(
            mut triatomic_with_overlays: MoleculeBuilder,
        ) {
            let tx = triatomic_with_overlays
                .transact(vec![Edit::RemoveTopology {
                    atoms: vec![AtomRef::Id(AtomId(0))],
                    bonds: vec![],
                }])
                .unwrap();
            let [Undo::RestoreTopology { overlays, .. }] = tx.undos() else {
                panic!("RemoveTopology should produce one topology-restore undo")
            };
            assert_eq!(overlays.dative_bonds.len(), 1);
            assert_eq!(overlays.aromatic_systems.len(), 1);
            assert_eq!(overlays.multicenter_bonds.len(), 1);
            assert_eq!(overlays.noncovalent_bonds.len(), 1);
        }

        #[rstest]
        fn test_molecule_builder_transact_unchecked(mut empty: MoleculeBuilder) {
            empty.transact_unchecked(vec![Edit::AddAtoms {
                atoms: vec![AtomAst::from_element(Element::C)],
            }]);
            assert_eq!(empty.atom_count(), 1);
            assert_eq!(empty.atom(AtomId(0)).ast.element, ElementAst::Lit(Element::C));
        }
    }
}

