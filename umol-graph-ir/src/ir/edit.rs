//! Edit vocabulary for transactional molecule mutation.
//!
//! The `Edit` enum is the caller-facing data-form mutation vocabulary; realized
//! rollback data belongs to the `Undo` journal.
//!
//! Handles (`AtomHandle`, `BondHandle`, ...) are symbolic. `Id(n)` names entity
//! `n` in the transaction's initial host; `New(n)` names the `n`th same-kind
//! entity created in the same [`Edits`] sequence.

use std::collections::{BTreeSet, HashMap};
use std::hash::Hash;
use std::slice::Iter;
use std::vec::IntoIter;

use thiserror::Error;

use super::aromatic::{AromaticSystemForm, AromaticSystemUpdate};
use super::atom::{AtomForm, AtomUpdate, ElementForm, IsotopeMassForm};
use super::bond::{BondForm, BondUpdate};
use super::constraint::{
    AromaticSystemConstraintForm, AtomConstraintForm, BondConstraintForm, Constraint,
    DativeBondConstraintForm, MoleculeConstraint, MulticenterBondConstraintForm,
    NoncovalentBondConstraintForm, RelationalConstraint, StereoAtomConstraintForm,
    StereoBondConstraintForm,
};
use super::dative::{DativeBondForm, DativeBondUpdate};
use super::electrons::ElectronCountsForm;
use super::entity::{Entity, EntityKind};
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::{StereoLigand, StereoLigandKind};
use super::multicenter::{MulticenterBondForm, MulticenterBondUpdate};
use super::noncovalent::{NoncovalentBondForm, NoncovalentBondKindForm, NoncovalentBondUpdate};
use super::remap::{IdCompaction, IdRemapping, UndoCompaction};
use super::spin::UnpairedElectronsForm;
use super::stereo::{
    StereoAtomForm, StereoAtomUpdate, StereoBondForm, StereoBondUpdate, StereoConfigurationForm,
    StereoKind,
};
use super::traits::{Canonicalize, Lattice};
use super::value::NumForm;

/// One stereo-atom removal in a batched `RemoveStereoAtoms`: id, site, ligand frame, recorded attributes.
pub type StereoAtomRemoval = (
    StereoAtomHandle,
    AtomHandle,
    Vec<(AtomHandle, StereoLigandKind)>,
    StereoAtomForm,
);
/// One stereo-bond removal in a batched `RemoveStereoBonds`: id, site (a bond), ligand frame, attributes.
pub type StereoBondRemoval = (
    StereoBondHandle,
    BondHandle,
    Vec<(AtomHandle, StereoLigandKind)>,
    StereoBondForm,
);

/// Handle to an atom within an edit batch: either an initial-host `AtomId` or a
/// same-kind creation ordinal in the enclosing [`Edits`] sequence.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AtomHandle {
    Id(AtomId),
    New(usize),
}

/// Handle to a bond within an edit batch: either an initial-host `BondId` or a
/// same-kind creation ordinal in the enclosing [`Edits`] sequence.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum BondHandle {
    Id(BondId),
    New(usize),
}

/// Per-field old/new payload for an atom attribute mutation. Variant
/// discriminant identifies the field; `old` and `new` carry the typed values.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum AtomFieldChange {
    Element {
        old: ElementForm,
        new: ElementForm,
    },
    IsotopeMass {
        old: IsotopeMassForm,
        new: IsotopeMassForm,
    },
    Charge {
        old: NumForm,
        new: NumForm,
    },
    ImplicitHydrogens {
        old: NumForm,
        new: NumForm,
    },
    LonePairs {
        old: NumForm,
        new: NumForm,
    },
    UnpairedElectrons {
        old: UnpairedElectronsForm,
        new: UnpairedElectronsForm,
    },
}

impl AtomFieldChange {
    /// Swap `old` and `new` across all variants.
    pub fn inverse(self) -> Self {
        match self {
            Self::Element { old, new } => Self::Element { old: new, new: old },
            Self::IsotopeMass { old, new } => Self::IsotopeMass { old: new, new: old },
            Self::Charge { old, new } => Self::Charge { old: new, new: old },
            Self::ImplicitHydrogens { old, new } => Self::ImplicitHydrogens { old: new, new: old },
            Self::LonePairs { old, new } => Self::LonePairs { old: new, new: old },
            Self::UnpairedElectrons { old, new } => Self::UnpairedElectrons { old: new, new: old },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum BondFieldChange {
    Order {
        old: NumForm,
        new: NumForm,
    },
    Charge {
        old: NumForm,
        new: NumForm,
    },
    UnpairedElectrons {
        old: UnpairedElectronsForm,
        new: UnpairedElectronsForm,
    },
}

impl BondFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Order { old, new } => Self::Order { old: new, new: old },
            Self::Charge { old, new } => Self::Charge { old: new, new: old },
            Self::UnpairedElectrons { old, new } => Self::UnpairedElectrons { old: new, new: old },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum DativeBondFieldChange {
    Order { old: NumForm, new: NumForm },
}

impl DativeBondFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Order { old, new } => Self::Order { old: new, new: old },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum AromaticSystemFieldChange {
    Electrons {
        old: ElectronCountsForm,
        new: ElectronCountsForm,
    },
    Charge {
        old: NumForm,
        new: NumForm,
    },
    UnpairedElectrons {
        old: UnpairedElectronsForm,
        new: UnpairedElectronsForm,
    },
}

impl AromaticSystemFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Electrons { old, new } => Self::Electrons { old: new, new: old },
            Self::Charge { old, new } => Self::Charge { old: new, new: old },
            Self::UnpairedElectrons { old, new } => Self::UnpairedElectrons { old: new, new: old },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum MulticenterBondFieldChange {
    Electrons {
        old: ElectronCountsForm,
        new: ElectronCountsForm,
    },
    Charge {
        old: NumForm,
        new: NumForm,
    },
    UnpairedElectrons {
        old: UnpairedElectronsForm,
        new: UnpairedElectronsForm,
    },
}

impl MulticenterBondFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Electrons { old, new } => Self::Electrons { old: new, new: old },
            Self::Charge { old, new } => Self::Charge { old: new, new: old },
            Self::UnpairedElectrons { old, new } => Self::UnpairedElectrons { old: new, new: old },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum NoncovalentBondFieldChange {
    Kind {
        old: NoncovalentBondKindForm,
        new: NoncovalentBondKindForm,
    },
}

impl NoncovalentBondFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Kind { old, new } => Self::Kind { old: new, new: old },
        }
    }
}

/// Per-field old/new payload for an absolute stereo-atom configuration mutation.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum StereoAtomFieldChange {
    Configuration {
        old: StereoConfigurationForm,
        new: StereoConfigurationForm,
    },
}

impl StereoAtomFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Configuration { old, new } => Self::Configuration { old: new, new: old },
        }
    }
}

/// Per-field old/new payload for a stereo-bond mutation. Coset-only, for the
/// same reason as `StereoAtomFieldChange`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum StereoBondFieldChange {
    Configuration {
        old: StereoConfigurationForm,
        new: StereoConfigurationForm,
    },
}

impl StereoBondFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Configuration { old, new } => Self::Configuration { old: new, new: old },
        }
    }
}

/// Single bond addition inside an `Edit::AddBonds` batch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddBond {
    pub endpoints: [AtomHandle; 2],
    pub attributes: BondForm,
}

/// One modification of a molecule the caller holds.
///
/// An edit is an imperative instruction, not an algebraic value. It refers to entities by their
/// index in the transaction's initial host, and to entities created earlier in the same sequence
/// as `New(n)`, since the host's concrete numbering is not known when the sequence is written.
///
/// Removal uses cascade deletion: taking an atom out of a ring also removes its incident bonds and
/// any aromatic system it belonged to. This is an execution behavior, not a complete definition of
/// SqPO rewriting. Because the edit need not name what it discards, and because concrete ids and
/// compaction are only known during application, an edit cannot be inverted on its own. Checked
/// application records an [`Undo`] as these effects are realized. `RemoveTopology` removes atoms
/// and bonds together.
///
/// Edits have no canonical form. Sorting or deduplicating a sequence would invalidate the `New(n)`
/// references that depend on its order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Edit {
    // Atoms / bonds
    AddAtoms {
        atoms: Vec<AtomForm>,
    },
    AddBonds {
        bonds: Vec<AddBond>,
    },
    RemoveTopology {
        atoms: Vec<AtomHandle>,
        bonds: Vec<BondHandle>,
    },
    ModifyAtomField {
        id: AtomHandle,
        change: AtomFieldChange,
    },
    ModifyBondField {
        id: BondHandle,
        change: BondFieldChange,
    },

    // Dative bonds
    AddDativeBond {
        atoms: Vec<AtomHandle>,
        attributes: DativeBondForm,
    },
    RemoveDativeBonds {
        removes: Vec<(DativeBondHandle, Vec<AtomHandle>, DativeBondForm)>,
    },
    ModifyDativeBondField {
        id: DativeBondHandle,
        change: DativeBondFieldChange,
    },

    // Aromatic systems
    AddAromaticSystem {
        atoms: Vec<AtomHandle>,
        attributes: AromaticSystemForm,
    },
    RemoveAromaticSystems {
        removes: Vec<(AromaticSystemHandle, Vec<AtomHandle>, AromaticSystemForm)>,
    },
    ModifyAromaticSystemField {
        id: AromaticSystemHandle,
        change: AromaticSystemFieldChange,
    },

    // Multicenter bonds
    AddMulticenterBond {
        atoms: Vec<AtomHandle>,
        attributes: MulticenterBondForm,
    },
    RemoveMulticenterBonds {
        removes: Vec<(MulticenterBondHandle, Vec<AtomHandle>, MulticenterBondForm)>,
    },
    ModifyMulticenterBondField {
        id: MulticenterBondHandle,
        change: MulticenterBondFieldChange,
    },

    // Noncovalent bonds
    AddNoncovalentBond {
        atoms: [AtomHandle; 2],
        attributes: NoncovalentBondForm,
    },
    RemoveNoncovalentBonds {
        removes: Vec<(NoncovalentBondHandle, [AtomHandle; 2], NoncovalentBondForm)>,
    },
    ModifyNoncovalentBondField {
        id: NoncovalentBondHandle,
        change: NoncovalentBondFieldChange,
    },

    // Stereo elements. `ligands` carry their atom as an `AtomHandle` (Id or
    // same-batch New) plus the ligand kind; `site` is the atom/bond the
    // element is sited on.
    AddStereoAtom {
        site: AtomHandle,
        ligands: Vec<(AtomHandle, StereoLigandKind)>,
        attributes: StereoAtomForm,
    },
    RemoveStereoAtoms {
        removes: Vec<StereoAtomRemoval>,
    },
    ModifyStereoAtomField {
        id: StereoAtomHandle,
        change: StereoAtomFieldChange,
    },
    AddStereoBond {
        site: BondHandle,
        ligands: Vec<(AtomHandle, StereoLigandKind)>,
        attributes: StereoBondForm,
    },
    RemoveStereoBonds {
        removes: Vec<StereoBondRemoval>,
    },
    ModifyStereoBondField {
        id: StereoBondHandle,
        change: StereoBondFieldChange,
    },

    // Entity-inline constraints — keyed (one per `key()`), so a single modify
    // (old → new) covers add (old None), remove (new None), and replace.
    ModifyAtomConstraint {
        id: AtomHandle,
        old: Option<AtomConstraintForm>,
        new: Option<AtomConstraintForm>,
    },
    ModifyBondConstraint {
        id: BondHandle,
        old: Option<BondConstraintForm>,
        new: Option<BondConstraintForm>,
    },
    ModifyDativeBondConstraint {
        id: DativeBondHandle,
        old: Option<DativeBondConstraintForm>,
        new: Option<DativeBondConstraintForm>,
    },
    ModifyAromaticSystemConstraint {
        id: AromaticSystemHandle,
        old: Option<AromaticSystemConstraintForm>,
        new: Option<AromaticSystemConstraintForm>,
    },
    ModifyMulticenterBondConstraint {
        id: MulticenterBondHandle,
        old: Option<MulticenterBondConstraintForm>,
        new: Option<MulticenterBondConstraintForm>,
    },
    ModifyNoncovalentBondConstraint {
        id: NoncovalentBondHandle,
        old: Option<NoncovalentBondConstraintForm>,
        new: Option<NoncovalentBondConstraintForm>,
    },
    ModifyStereoAtomConstraint {
        id: StereoAtomHandle,
        /// Geometry context required to parse and render the constraint DSL.
        kind: Option<StereoKind>,
        old: Option<StereoAtomConstraintForm>,
        new: Option<StereoAtomConstraintForm>,
    },
    ModifyStereoBondConstraint {
        id: StereoBondHandle,
        /// Geometry context required to parse and render the constraint DSL.
        kind: Option<StereoKind>,
        old: Option<StereoBondConstraintForm>,
        new: Option<StereoBondConstraintForm>,
    },

    // Molecule-list constraints — a true multiset, so add/remove by value
    // (remove takes the last matching entry; its position is captured for undo).
    AddMoleculeConstraint {
        constraint: ConstraintEdit,
    },
    RemoveMoleculeConstraint {
        constraint: ConstraintEdit,
    },
}

/// An ordered batch of host-specific molecule edits.
///
/// Order is semantic: later entries may refer to entities created by earlier entries. For every
/// entity kind, `Id(n)` names entity `n` in the transaction's initial host and `New(n)` names the
/// `n`th same-kind creation in this sequence. Creation ordinals are independent between kinds,
/// never reused, and are not changed by removals. `Edits` issues these symbolic handles only;
/// transaction application owns their concrete ids, liveness, and compaction in a particular host.
///
/// The public mutation surface is append-only. Mutable iteration, insertion, removal, reordering,
/// and concatenation are deliberately absent because they could invalidate issued `New(n)` handles.
/// Checked application resolves the handles against one host and returns the realized undo journal
/// as a [`Transaction`](crate::ir::Transaction).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Edits {
    edits: Vec<Edit>,
    created_atoms: usize,
    created_bonds: usize,
    created_dative_bonds: usize,
    created_aromatic_systems: usize,
    created_multicenter_bonds: usize,
    created_noncovalent_bonds: usize,
    created_stereo_atoms: usize,
    created_stereo_bonds: usize,
}

impl Edits {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn as_slice(&self) -> &[Edit] {
        &self.edits
    }

    pub fn is_empty(&self) -> bool {
        self.edits.is_empty()
    }

    pub fn len(&self) -> usize {
        self.edits.len()
    }

    pub fn iter(&self) -> Iter<'_, Edit> {
        self.edits.iter()
    }

    /// Append one raw entry and account for every entity it creates.
    pub fn push(&mut self, edit: Edit) {
        match &edit {
            Edit::AddAtoms { atoms } => self.created_atoms += atoms.len(),
            Edit::AddBonds { bonds } => self.created_bonds += bonds.len(),
            Edit::AddDativeBond { .. } => self.created_dative_bonds += 1,
            Edit::AddAromaticSystem { .. } => self.created_aromatic_systems += 1,
            Edit::AddMulticenterBond { .. } => self.created_multicenter_bonds += 1,
            Edit::AddNoncovalentBond { .. } => self.created_noncovalent_bonds += 1,
            Edit::AddStereoAtom { .. } => self.created_stereo_atoms += 1,
            Edit::AddStereoBond { .. } => self.created_stereo_bonds += 1,
            _ => {}
        }
        self.edits.push(edit);
    }

    pub fn add_atom(&mut self, attributes: AtomForm) -> AtomHandle {
        let handle = AtomHandle::New(self.created_atoms);
        self.push(Edit::AddAtoms {
            atoms: vec![attributes],
        });
        handle
    }

    pub fn add_atoms(&mut self, atoms: impl IntoIterator<Item = AtomForm>) -> Vec<AtomHandle> {
        let atoms: Vec<_> = atoms.into_iter().collect();
        let handles = (self.created_atoms..self.created_atoms + atoms.len())
            .map(AtomHandle::New)
            .collect();
        self.push(Edit::AddAtoms { atoms });
        handles
    }

    pub fn add_bond(
        &mut self,
        first: AtomHandle,
        second: AtomHandle,
        attributes: BondForm,
    ) -> BondHandle {
        let handle = BondHandle::New(self.created_bonds);
        self.push(Edit::AddBonds {
            bonds: vec![AddBond {
                endpoints: [first, second],
                attributes,
            }],
        });
        handle
    }

    pub fn add_bonds(&mut self, bonds: impl IntoIterator<Item = AddBond>) -> Vec<BondHandle> {
        let bonds: Vec<_> = bonds.into_iter().collect();
        let handles = (self.created_bonds..self.created_bonds + bonds.len())
            .map(BondHandle::New)
            .collect();
        self.push(Edit::AddBonds { bonds });
        handles
    }

    pub fn add_dative_bond(
        &mut self,
        atoms: Vec<AtomHandle>,
        attributes: DativeBondForm,
    ) -> DativeBondHandle {
        let handle = DativeBondHandle::New(self.created_dative_bonds);
        self.push(Edit::AddDativeBond { atoms, attributes });
        handle
    }

    pub fn add_dative_bonds(
        &mut self,
        bonds: impl IntoIterator<Item = (Vec<AtomHandle>, DativeBondForm)>,
    ) -> Vec<DativeBondHandle> {
        bonds
            .into_iter()
            .map(|(atoms, attributes)| self.add_dative_bond(atoms, attributes))
            .collect()
    }

    pub fn add_aromatic_system(
        &mut self,
        atoms: Vec<AtomHandle>,
        attributes: AromaticSystemForm,
    ) -> AromaticSystemHandle {
        let handle = AromaticSystemHandle::New(self.created_aromatic_systems);
        self.push(Edit::AddAromaticSystem { atoms, attributes });
        handle
    }

    pub fn add_aromatic_systems(
        &mut self,
        systems: impl IntoIterator<Item = (Vec<AtomHandle>, AromaticSystemForm)>,
    ) -> Vec<AromaticSystemHandle> {
        systems
            .into_iter()
            .map(|(atoms, attributes)| self.add_aromatic_system(atoms, attributes))
            .collect()
    }

    pub fn add_multicenter_bond(
        &mut self,
        atoms: Vec<AtomHandle>,
        attributes: MulticenterBondForm,
    ) -> MulticenterBondHandle {
        let handle = MulticenterBondHandle::New(self.created_multicenter_bonds);
        self.push(Edit::AddMulticenterBond { atoms, attributes });
        handle
    }

    pub fn add_multicenter_bonds(
        &mut self,
        bonds: impl IntoIterator<Item = (Vec<AtomHandle>, MulticenterBondForm)>,
    ) -> Vec<MulticenterBondHandle> {
        bonds
            .into_iter()
            .map(|(atoms, attributes)| self.add_multicenter_bond(atoms, attributes))
            .collect()
    }

    pub fn add_noncovalent_bond(
        &mut self,
        atoms: [AtomHandle; 2],
        attributes: NoncovalentBondForm,
    ) -> NoncovalentBondHandle {
        let handle = NoncovalentBondHandle::New(self.created_noncovalent_bonds);
        self.push(Edit::AddNoncovalentBond { atoms, attributes });
        handle
    }

    pub fn add_noncovalent_bonds(
        &mut self,
        bonds: impl IntoIterator<Item = ([AtomHandle; 2], NoncovalentBondForm)>,
    ) -> Vec<NoncovalentBondHandle> {
        bonds
            .into_iter()
            .map(|(atoms, attributes)| self.add_noncovalent_bond(atoms, attributes))
            .collect()
    }

    pub fn add_stereo_atom(
        &mut self,
        site: AtomHandle,
        ligands: Vec<(AtomHandle, StereoLigandKind)>,
        attributes: StereoAtomForm,
    ) -> StereoAtomHandle {
        let handle = StereoAtomHandle::New(self.created_stereo_atoms);
        self.push(Edit::AddStereoAtom {
            site,
            ligands,
            attributes,
        });
        handle
    }

    pub fn add_stereo_atoms(
        &mut self,
        atoms: impl IntoIterator<
            Item = (
                AtomHandle,
                Vec<(AtomHandle, StereoLigandKind)>,
                StereoAtomForm,
            ),
        >,
    ) -> Vec<StereoAtomHandle> {
        atoms
            .into_iter()
            .map(|(site, ligands, attributes)| self.add_stereo_atom(site, ligands, attributes))
            .collect()
    }

    pub fn add_stereo_bond(
        &mut self,
        site: BondHandle,
        ligands: Vec<(AtomHandle, StereoLigandKind)>,
        attributes: StereoBondForm,
    ) -> StereoBondHandle {
        let handle = StereoBondHandle::New(self.created_stereo_bonds);
        self.push(Edit::AddStereoBond {
            site,
            ligands,
            attributes,
        });
        handle
    }

    pub fn add_stereo_bonds(
        &mut self,
        bonds: impl IntoIterator<
            Item = (
                BondHandle,
                Vec<(AtomHandle, StereoLigandKind)>,
                StereoBondForm,
            ),
        >,
    ) -> Vec<StereoBondHandle> {
        bonds
            .into_iter()
            .map(|(site, ligands, attributes)| self.add_stereo_bond(site, ligands, attributes))
            .collect()
    }

    pub fn remove_atom(&mut self, id: AtomHandle) {
        self.remove_topology(vec![id], Vec::new());
    }

    pub fn remove_bond(&mut self, id: BondHandle) {
        self.remove_topology(Vec::new(), vec![id]);
    }

    pub fn remove_topology(&mut self, atoms: Vec<AtomHandle>, bonds: Vec<BondHandle>) {
        self.push(Edit::RemoveTopology { atoms, bonds });
    }

    pub fn remove_dative_bonds(
        &mut self,
        removes: Vec<(DativeBondHandle, Vec<AtomHandle>, DativeBondForm)>,
    ) {
        self.push(Edit::RemoveDativeBonds { removes });
    }

    pub fn remove_aromatic_systems(
        &mut self,
        removes: Vec<(AromaticSystemHandle, Vec<AtomHandle>, AromaticSystemForm)>,
    ) {
        self.push(Edit::RemoveAromaticSystems { removes });
    }

    pub fn remove_multicenter_bonds(
        &mut self,
        removes: Vec<(MulticenterBondHandle, Vec<AtomHandle>, MulticenterBondForm)>,
    ) {
        self.push(Edit::RemoveMulticenterBonds { removes });
    }

    pub fn remove_noncovalent_bonds(
        &mut self,
        removes: Vec<(NoncovalentBondHandle, [AtomHandle; 2], NoncovalentBondForm)>,
    ) {
        self.push(Edit::RemoveNoncovalentBonds { removes });
    }

    pub fn remove_stereo_atoms(&mut self, removes: Vec<StereoAtomRemoval>) {
        self.push(Edit::RemoveStereoAtoms { removes });
    }

    pub fn remove_stereo_bonds(&mut self, removes: Vec<StereoBondRemoval>) {
        self.push(Edit::RemoveStereoBonds { removes });
    }

    pub fn add_molecule_constraint(&mut self, constraint: ConstraintEdit) {
        self.push(Edit::AddMoleculeConstraint { constraint });
    }

    pub fn remove_molecule_constraint(&mut self, constraint: ConstraintEdit) {
        self.push(Edit::RemoveMoleculeConstraint { constraint });
    }
}

impl FromIterator<Edit> for Edits {
    fn from_iter<I: IntoIterator<Item = Edit>>(iter: I) -> Self {
        let mut edits = Self::new();
        for edit in iter {
            edits.push(edit);
        }
        edits
    }
}

impl IntoIterator for Edits {
    type Item = Edit;
    type IntoIter = IntoIter<Edit>;

    fn into_iter(self) -> Self::IntoIter {
        self.edits.into_iter()
    }
}

impl Edits {
    /// Project an atom update into checked host-relative edits.
    pub fn update_atom(&mut self, id: AtomHandle, current: &AtomForm, update: &AtomUpdate) {
        if let Some(new) = &update.element {
            if !current.element.canonical_eq(new) {
                self.push(Edit::ModifyAtomField {
                    id: id.clone(),
                    change: AtomFieldChange::Element {
                        old: current.element.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.isotope_mass {
            if !current.isotope_mass.canonical_eq(new) {
                self.push(Edit::ModifyAtomField {
                    id: id.clone(),
                    change: AtomFieldChange::IsotopeMass {
                        old: current.isotope_mass.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.charge {
            if !current.charge.canonical_eq(new) {
                self.push(Edit::ModifyAtomField {
                    id: id.clone(),
                    change: AtomFieldChange::Charge {
                        old: current.charge.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.implicit_hydrogens {
            if !current.implicit_hydrogens.canonical_eq(new) {
                self.push(Edit::ModifyAtomField {
                    id: id.clone(),
                    change: AtomFieldChange::ImplicitHydrogens {
                        old: current.implicit_hydrogens.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.lone_pairs {
            if !current.lone_pairs.canonical_eq(new) {
                self.push(Edit::ModifyAtomField {
                    id: id.clone(),
                    change: AtomFieldChange::LonePairs {
                        old: current.lone_pairs.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        let new_unpaired_electrons = current
            .unpaired_electrons
            .update(&update.unpaired_electrons);
        if !current
            .unpaired_electrons
            .canonical_eq(&new_unpaired_electrons)
        {
            self.push(Edit::ModifyAtomField {
                id: id.clone(),
                change: AtomFieldChange::UnpairedElectrons {
                    old: current.unpaired_electrons.clone(),
                    new: new_unpaired_electrons,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            let unchanged = match (&old, &new) {
                (None, None) => true,
                (Some(old), Some(new)) => old.canonical_eq(new),
                _ => false,
            };
            if !unchanged {
                self.push(Edit::ModifyAtomConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
    }

    /// Project a localized-bond update into checked host-relative edits.
    pub fn update_bond(&mut self, id: BondHandle, current: &BondForm, update: &BondUpdate) {
        if let Some(new) = &update.order {
            if !current.order.canonical_eq(new) {
                self.push(Edit::ModifyBondField {
                    id: id.clone(),
                    change: BondFieldChange::Order {
                        old: current.order.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.charge {
            if !current.charge.canonical_eq(new) {
                self.push(Edit::ModifyBondField {
                    id: id.clone(),
                    change: BondFieldChange::Charge {
                        old: current.charge.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        let new_unpaired_electrons = current
            .unpaired_electrons
            .update(&update.unpaired_electrons);
        if !current
            .unpaired_electrons
            .canonical_eq(&new_unpaired_electrons)
        {
            self.push(Edit::ModifyBondField {
                id: id.clone(),
                change: BondFieldChange::UnpairedElectrons {
                    old: current.unpaired_electrons.clone(),
                    new: new_unpaired_electrons,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            let unchanged = match (&old, &new) {
                (None, None) => true,
                (Some(old), Some(new)) => old.canonical_eq(new),
                _ => false,
            };
            if !unchanged {
                self.push(Edit::ModifyBondConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
    }

    /// Project a dative-bond update into checked host-relative edits.
    pub fn update_dative_bond(
        &mut self,
        id: DativeBondHandle,
        current: &DativeBondForm,
        update: &DativeBondUpdate,
    ) {
        if let Some(new) = &update.order {
            if !current.order.canonical_eq(new) {
                self.push(Edit::ModifyDativeBondField {
                    id: id.clone(),
                    change: DativeBondFieldChange::Order {
                        old: current.order.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            let unchanged = match (&old, &new) {
                (None, None) => true,
                (Some(old), Some(new)) => old.canonical_eq(new),
                _ => false,
            };
            if !unchanged {
                self.push(Edit::ModifyDativeBondConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
    }

    /// Project an aromatic-system update into checked host-relative edits.
    pub fn update_aromatic_system(
        &mut self,
        id: AromaticSystemHandle,
        current: &AromaticSystemForm,
        update: &AromaticSystemUpdate,
    ) {
        if let Some(new) = &update.electrons {
            if !current.electrons.canonical_eq(new) {
                self.push(Edit::ModifyAromaticSystemField {
                    id: id.clone(),
                    change: AromaticSystemFieldChange::Electrons {
                        old: current.electrons.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.charge {
            if !current.charge.canonical_eq(new) {
                self.push(Edit::ModifyAromaticSystemField {
                    id: id.clone(),
                    change: AromaticSystemFieldChange::Charge {
                        old: current.charge.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        let new_unpaired_electrons = current
            .unpaired_electrons
            .update(&update.unpaired_electrons);
        if !current
            .unpaired_electrons
            .canonical_eq(&new_unpaired_electrons)
        {
            self.push(Edit::ModifyAromaticSystemField {
                id: id.clone(),
                change: AromaticSystemFieldChange::UnpairedElectrons {
                    old: current.unpaired_electrons.clone(),
                    new: new_unpaired_electrons,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            let unchanged = match (&old, &new) {
                (None, None) => true,
                (Some(old), Some(new)) => old.canonical_eq(new),
                _ => false,
            };
            if !unchanged {
                self.push(Edit::ModifyAromaticSystemConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
    }

    /// Project a multicenter-bond update into checked host-relative edits.
    pub fn update_multicenter_bond(
        &mut self,
        id: MulticenterBondHandle,
        current: &MulticenterBondForm,
        update: &MulticenterBondUpdate,
    ) {
        if let Some(new) = &update.electrons {
            if !current.electrons.canonical_eq(new) {
                self.push(Edit::ModifyMulticenterBondField {
                    id: id.clone(),
                    change: MulticenterBondFieldChange::Electrons {
                        old: current.electrons.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.charge {
            if !current.charge.canonical_eq(new) {
                self.push(Edit::ModifyMulticenterBondField {
                    id: id.clone(),
                    change: MulticenterBondFieldChange::Charge {
                        old: current.charge.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        let new_unpaired_electrons = current
            .unpaired_electrons
            .update(&update.unpaired_electrons);
        if !current
            .unpaired_electrons
            .canonical_eq(&new_unpaired_electrons)
        {
            self.push(Edit::ModifyMulticenterBondField {
                id: id.clone(),
                change: MulticenterBondFieldChange::UnpairedElectrons {
                    old: current.unpaired_electrons.clone(),
                    new: new_unpaired_electrons,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            let unchanged = match (&old, &new) {
                (None, None) => true,
                (Some(old), Some(new)) => old.canonical_eq(new),
                _ => false,
            };
            if !unchanged {
                self.push(Edit::ModifyMulticenterBondConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
    }

    /// Project a noncovalent-bond update into checked host-relative edits.
    pub fn update_noncovalent_bond(
        &mut self,
        id: NoncovalentBondHandle,
        current: &NoncovalentBondForm,
        update: &NoncovalentBondUpdate,
    ) {
        if let Some(new) = &update.kind {
            if !current.kind.canonical_eq(new) {
                self.push(Edit::ModifyNoncovalentBondField {
                    id: id.clone(),
                    change: NoncovalentBondFieldChange::Kind {
                        old: current.kind.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            let unchanged = match (&old, &new) {
                (None, None) => true,
                (Some(old), Some(new)) => old.canonical_eq(new),
                _ => false,
            };
            if !unchanged {
                self.push(Edit::ModifyNoncovalentBondConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
    }

    /// Project a stereo-atom update into checked host-relative edits.
    pub fn update_stereo_atom(
        &mut self,
        id: StereoAtomHandle,
        current: &StereoAtomForm,
        update: &StereoAtomUpdate,
    ) {
        let updated = current.update(update);
        let kind = update
            .configuration
            .kind()
            .or_else(|| current.configuration.kind())
            .or_else(|| updated.configuration.kind());
        if !current.configuration.canonical_eq(&updated.configuration) {
            self.push(Edit::ModifyStereoAtomField {
                id: id.clone(),
                change: StereoAtomFieldChange::Configuration {
                    old: current.configuration.clone(),
                    new: updated.configuration,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            let unchanged = match (&old, &new) {
                (None, None) => true,
                (Some(old), Some(new)) => old.canonical_eq(new),
                _ => false,
            };
            if !unchanged {
                self.push(Edit::ModifyStereoAtomConstraint {
                    id: id.clone(),
                    kind,
                    old,
                    new,
                });
            }
        }
    }

    /// Project a stereo-bond update into checked host-relative edits.
    pub fn update_stereo_bond(
        &mut self,
        id: StereoBondHandle,
        current: &StereoBondForm,
        update: &StereoBondUpdate,
    ) {
        let updated = current.update(update);
        let kind = update
            .configuration
            .kind()
            .or_else(|| current.configuration.kind())
            .or_else(|| updated.configuration.kind());
        if !current.configuration.canonical_eq(&updated.configuration) {
            self.push(Edit::ModifyStereoBondField {
                id: id.clone(),
                change: StereoBondFieldChange::Configuration {
                    old: current.configuration.clone(),
                    new: updated.configuration,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            let unchanged = match (&old, &new) {
                (None, None) => true,
                (Some(old), Some(new)) => old.canonical_eq(new),
                _ => false,
            };
            if !unchanged {
                self.push(Edit::ModifyStereoBondConstraint {
                    id: id.clone(),
                    kind,
                    old,
                    new,
                });
            }
        }
    }
}

// Handles for overlay relations: an initial-host id or a same-kind creation ordinal.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum DativeBondHandle {
    Id(DativeBondId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemHandle {
    Id(AromaticSystemId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MulticenterBondHandle {
    Id(MulticenterBondId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondHandle {
    Id(NoncovalentBondId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum StereoAtomHandle {
    Id(StereoAtomId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum StereoBondHandle {
    Id(StereoBondId),
    New(usize),
}

/// Typed handle to any molecule entity while constructing a [`ConstraintEdit`].
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum EntityHandle {
    Atom(AtomHandle),
    Bond(BondHandle),
    DativeBond(DativeBondHandle),
    AromaticSystem(AromaticSystemHandle),
    MulticenterBond(MulticenterBondHandle),
    NoncovalentBond(NoncovalentBondHandle),
    StereoAtom(StereoAtomHandle),
    StereoBond(StereoBondHandle),
}

impl EntityHandle {
    pub fn kind(&self) -> EntityKind {
        match self {
            Self::Atom(_) => EntityKind::Atom,
            Self::Bond(_) => EntityKind::Bond,
            Self::DativeBond(_) => EntityKind::DativeBond,
            Self::AromaticSystem(_) => EntityKind::AromaticSystem,
            Self::MulticenterBond(_) => EntityKind::MulticenterBond,
            Self::NoncovalentBond(_) => EntityKind::NoncovalentBond,
            Self::StereoAtom(_) => EntityKind::StereoAtom,
            Self::StereoBond(_) => EntityKind::StereoBond,
        }
    }
}

impl From<Entity> for EntityHandle {
    fn from(entity: Entity) -> Self {
        match entity {
            Entity::Atom(id) => Self::Atom(AtomHandle::Id(id)),
            Entity::Bond(id) => Self::Bond(BondHandle::Id(id)),
            Entity::DativeBond(id) => Self::DativeBond(DativeBondHandle::Id(id)),
            Entity::AromaticSystem(id) => Self::AromaticSystem(AromaticSystemHandle::Id(id)),
            Entity::MulticenterBond(id) => Self::MulticenterBond(MulticenterBondHandle::Id(id)),
            Entity::NoncovalentBond(id) => Self::NoncovalentBond(NoncovalentBondHandle::Id(id)),
            Entity::StereoAtom(id) => Self::StereoAtom(StereoAtomHandle::Id(id)),
            Entity::StereoBond(id) => Self::StereoBond(StereoBondHandle::Id(id)),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Error)]
pub enum ConstraintEditError {
    #[error("missing handle for {entity}")]
    MissingHandle { entity: Entity },

    #[error("handle kind mismatch for {entity}: found {actual}")]
    HandleKindMismatch { entity: Entity, actual: EntityKind },
}

/// Molecule-level constraint whose target-molecule references are stable edit handles.
///
/// The stored constraint uses normalized, per-kind slot ids. Each slot indexes the corresponding
/// private handle vector.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConstraintEdit {
    constraint: Constraint,
    atoms: Vec<AtomHandle>,
    bonds: Vec<BondHandle>,
    dative_bonds: Vec<DativeBondHandle>,
    aromatic_systems: Vec<AromaticSystemHandle>,
    multicenter_bonds: Vec<MulticenterBondHandle>,
    noncovalent_bonds: Vec<NoncovalentBondHandle>,
    stereo_atoms: Vec<StereoAtomHandle>,
    stereo_bonds: Vec<StereoBondHandle>,
}

impl ConstraintEdit {
    /// Normalize a concrete constraint through a complete mapping of its target-molecule entities.
    ///
    /// The mapping is requested once for each distinct referenced entity. Repeated mapped handles
    /// share one normalized slot within their entity kind.
    pub fn new(
        constraint: Constraint,
        mut handle_for: impl FnMut(Entity) -> Option<EntityHandle>,
    ) -> Result<Self, ConstraintEditError> {
        let mut atom_map = HashMap::new();
        let mut bond_map = HashMap::new();
        let mut dative_map = HashMap::new();
        let mut aromatic_map = HashMap::new();
        let mut multicenter_map = HashMap::new();
        let mut noncovalent_map = HashMap::new();
        let mut stereo_atom_map = HashMap::new();
        let mut stereo_bond_map = HashMap::new();
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut dative_bonds = Vec::new();
        let mut aromatic_systems = Vec::new();
        let mut multicenter_bonds = Vec::new();
        let mut noncovalent_bonds = Vec::new();
        let mut stereo_atoms = Vec::new();
        let mut stereo_bonds = Vec::new();

        let mut entities = BTreeSet::new();
        collect_constraint_entities(&constraint, &mut entities);
        for entity in entities {
            let handle = handle_for(entity).ok_or(ConstraintEditError::MissingHandle { entity })?;
            match (entity, handle) {
                (Entity::Atom(id), EntityHandle::Atom(handle)) => {
                    intern_handle(id, handle, &mut atoms, &mut atom_map);
                }
                (Entity::Bond(id), EntityHandle::Bond(handle)) => {
                    intern_handle(id, handle, &mut bonds, &mut bond_map);
                }
                (Entity::DativeBond(id), EntityHandle::DativeBond(handle)) => {
                    intern_handle(id, handle, &mut dative_bonds, &mut dative_map);
                }
                (Entity::AromaticSystem(id), EntityHandle::AromaticSystem(handle)) => {
                    intern_handle(id, handle, &mut aromatic_systems, &mut aromatic_map);
                }
                (Entity::MulticenterBond(id), EntityHandle::MulticenterBond(handle)) => {
                    intern_handle(id, handle, &mut multicenter_bonds, &mut multicenter_map);
                }
                (Entity::NoncovalentBond(id), EntityHandle::NoncovalentBond(handle)) => {
                    intern_handle(id, handle, &mut noncovalent_bonds, &mut noncovalent_map);
                }
                (Entity::StereoAtom(id), EntityHandle::StereoAtom(handle)) => {
                    intern_handle(id, handle, &mut stereo_atoms, &mut stereo_atom_map);
                }
                (Entity::StereoBond(id), EntityHandle::StereoBond(handle)) => {
                    intern_handle(id, handle, &mut stereo_bonds, &mut stereo_bond_map);
                }
                (entity, handle) => {
                    return Err(ConstraintEditError::HandleKindMismatch {
                        entity,
                        actual: handle.kind(),
                    });
                }
            }
        }

        let remapping = IdRemapping::new(
            atom_map,
            bond_map,
            dative_map,
            aromatic_map,
            multicenter_map,
            noncovalent_map,
            stereo_atom_map,
            stereo_bond_map,
        );
        Ok(Self {
            constraint: constraint.remap(&remapping),
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            stereo_atoms,
            stereo_bonds,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn resolve<E>(
        self,
        mut atom: impl FnMut(AtomHandle) -> Result<AtomId, E>,
        mut bond: impl FnMut(BondHandle) -> Result<BondId, E>,
        mut dative_bond: impl FnMut(DativeBondHandle) -> Result<DativeBondId, E>,
        mut aromatic_system: impl FnMut(AromaticSystemHandle) -> Result<AromaticSystemId, E>,
        mut multicenter_bond: impl FnMut(MulticenterBondHandle) -> Result<MulticenterBondId, E>,
        mut noncovalent_bond: impl FnMut(NoncovalentBondHandle) -> Result<NoncovalentBondId, E>,
        mut stereo_atom: impl FnMut(StereoAtomHandle) -> Result<StereoAtomId, E>,
        mut stereo_bond: impl FnMut(StereoBondHandle) -> Result<StereoBondId, E>,
    ) -> Result<Constraint, E> {
        let remapping = IdRemapping::new(
            self.atoms
                .into_iter()
                .enumerate()
                .map(|(slot, handle)| atom(handle).map(|id| (AtomId::from(slot), id)))
                .collect::<Result<_, _>>()?,
            self.bonds
                .into_iter()
                .enumerate()
                .map(|(slot, handle)| bond(handle).map(|id| (BondId::from(slot), id)))
                .collect::<Result<_, _>>()?,
            self.dative_bonds
                .into_iter()
                .enumerate()
                .map(|(slot, handle)| dative_bond(handle).map(|id| (DativeBondId::from(slot), id)))
                .collect::<Result<_, _>>()?,
            self.aromatic_systems
                .into_iter()
                .enumerate()
                .map(|(slot, handle)| {
                    aromatic_system(handle).map(|id| (AromaticSystemId::from(slot), id))
                })
                .collect::<Result<_, _>>()?,
            self.multicenter_bonds
                .into_iter()
                .enumerate()
                .map(|(slot, handle)| {
                    multicenter_bond(handle).map(|id| (MulticenterBondId::from(slot), id))
                })
                .collect::<Result<_, _>>()?,
            self.noncovalent_bonds
                .into_iter()
                .enumerate()
                .map(|(slot, handle)| {
                    noncovalent_bond(handle).map(|id| (NoncovalentBondId::from(slot), id))
                })
                .collect::<Result<_, _>>()?,
            self.stereo_atoms
                .into_iter()
                .enumerate()
                .map(|(slot, handle)| stereo_atom(handle).map(|id| (StereoAtomId::from(slot), id)))
                .collect::<Result<_, _>>()?,
            self.stereo_bonds
                .into_iter()
                .enumerate()
                .map(|(slot, handle)| stereo_bond(handle).map(|id| (StereoBondId::from(slot), id)))
                .collect::<Result<_, _>>()?,
        );
        Ok(self.constraint.remap(&remapping))
    }
}

impl From<Constraint> for ConstraintEdit {
    fn from(constraint: Constraint) -> Self {
        ConstraintEdit::new(constraint, |entity| Some(EntityHandle::from(entity)))
            .expect("an entity's identity handle has the same kind")
    }
}

fn intern_handle<I, H>(source: I, handle: H, handles: &mut Vec<H>, remapping: &mut HashMap<I, I>)
where
    I: Copy + Eq + Hash + From<usize>,
    H: Eq,
{
    let slot = handles
        .iter()
        .position(|candidate| candidate == &handle)
        .unwrap_or_else(|| {
            handles.push(handle);
            handles.len() - 1
        });
    remapping.insert(source, I::from(slot));
}

fn collect_constraint_entities(constraint: &Constraint, entities: &mut BTreeSet<Entity>) {
    match constraint {
        Constraint::Atom(id, _) => {
            entities.insert(Entity::Atom(*id));
        }
        Constraint::Bond(id, _) => {
            entities.insert(Entity::Bond(*id));
        }
        Constraint::DativeBond(id, _) => {
            entities.insert(Entity::DativeBond(*id));
        }
        Constraint::AromaticSystem(id, _) => {
            entities.insert(Entity::AromaticSystem(*id));
        }
        Constraint::MulticenterBond(id, _) => {
            entities.insert(Entity::MulticenterBond(*id));
        }
        Constraint::NoncovalentBond(id, _) => {
            entities.insert(Entity::NoncovalentBond(*id));
        }
        Constraint::StereoAtom(id, _, _) => {
            entities.insert(Entity::StereoAtom(*id));
        }
        Constraint::StereoBond(id, _, _) => {
            entities.insert(Entity::StereoBond(*id));
        }
        Constraint::Relational(constraint) => {
            collect_relational_constraint_entities(constraint, entities);
        }
        Constraint::Molecule(constraint) => {
            collect_molecule_constraint_entities(constraint, entities);
        }
        Constraint::And(constraints) | Constraint::Or(constraints) => {
            for constraint in constraints {
                collect_constraint_entities(constraint, entities);
            }
        }
        Constraint::Not(constraint) => collect_constraint_entities(constraint, entities),
    }
}

fn collect_relational_constraint_entities(
    constraint: &RelationalConstraint,
    entities: &mut BTreeSet<Entity>,
) {
    match constraint {
        RelationalConstraint::DativeBondDonors { bond, atoms }
        | RelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
            entities.insert(Entity::DativeBond(*bond));
            entities.extend(atoms.iter().copied().map(Entity::Atom));
        }
        RelationalConstraint::DativeBondDonor { bond, atom }
        | RelationalConstraint::DativeBondAcceptor { bond, atom } => {
            entities.insert(Entity::DativeBond(*bond));
            entities.insert(Entity::Atom(*atom));
        }
        RelationalConstraint::DativeBondAllDonors { bond, .. }
        | RelationalConstraint::DativeBondAnyDonor { bond, .. }
        | RelationalConstraint::DativeBondAcceptorSatisfies { bond, .. } => {
            entities.insert(Entity::DativeBond(*bond));
        }
        RelationalConstraint::DativeBondParallels { dative, parallel } => {
            entities.insert(Entity::DativeBond(*dative));
            entities.insert(Entity::Bond(*parallel));
        }
        RelationalConstraint::AromaticSystemAtoms { system, atoms }
        | RelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
            entities.insert(Entity::AromaticSystem(*system));
            entities.extend(atoms.iter().copied().map(Entity::Atom));
        }
        RelationalConstraint::AromaticSystemContains { system, atom } => {
            entities.insert(Entity::AromaticSystem(*system));
            entities.insert(Entity::Atom(*atom));
        }
        RelationalConstraint::AromaticSystemAllAtoms { system, .. }
        | RelationalConstraint::AromaticSystemAnyAtom { system, .. } => {
            entities.insert(Entity::AromaticSystem(*system));
        }
        RelationalConstraint::MulticenterBondAtoms { bond, atoms }
        | RelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
            entities.insert(Entity::MulticenterBond(*bond));
            entities.extend(atoms.iter().copied().map(Entity::Atom));
        }
        RelationalConstraint::MulticenterBondContains { bond, atom } => {
            entities.insert(Entity::MulticenterBond(*bond));
            entities.insert(Entity::Atom(*atom));
        }
        RelationalConstraint::MulticenterBondAllAtoms { bond, .. }
        | RelationalConstraint::MulticenterBondAnyAtom { bond, .. } => {
            entities.insert(Entity::MulticenterBond(*bond));
        }
        RelationalConstraint::NoncovalentBondEnds { bond, atoms } => {
            entities.insert(Entity::NoncovalentBond(*bond));
            entities.extend(atoms.iter().copied().map(Entity::Atom));
        }
        RelationalConstraint::NoncovalentBondContains { bond, atom } => {
            entities.insert(Entity::NoncovalentBond(*bond));
            entities.insert(Entity::Atom(*atom));
        }
        RelationalConstraint::NoncovalentBondEndsSatisfy { bond, .. } => {
            entities.insert(Entity::NoncovalentBond(*bond));
        }
        RelationalConstraint::StereoAtomSite { stereo_atom, atom }
        | RelationalConstraint::StereoAtomContains { stereo_atom, atom } => {
            entities.insert(Entity::StereoAtom(*stereo_atom));
            entities.insert(Entity::Atom(*atom));
        }
        RelationalConstraint::StereoAtomLigands { stereo_atom, atoms } => {
            entities.insert(Entity::StereoAtom(*stereo_atom));
            entities.extend(atoms.iter().copied().map(Entity::Atom));
        }
        RelationalConstraint::StereoAtomAllLigands { stereo_atom, .. }
        | RelationalConstraint::StereoAtomAnyLigand { stereo_atom, .. } => {
            entities.insert(Entity::StereoAtom(*stereo_atom));
        }
        RelationalConstraint::StereoBondSite { stereo_bond, bond } => {
            entities.insert(Entity::StereoBond(*stereo_bond));
            entities.insert(Entity::Bond(*bond));
        }
        RelationalConstraint::StereoBondContains { stereo_bond, atom } => {
            entities.insert(Entity::StereoBond(*stereo_bond));
            entities.insert(Entity::Atom(*atom));
        }
        RelationalConstraint::StereoBondLigands { stereo_bond, atoms } => {
            entities.insert(Entity::StereoBond(*stereo_bond));
            entities.extend(atoms.iter().copied().map(Entity::Atom));
        }
        RelationalConstraint::StereoBondAllLigands { stereo_bond, .. }
        | RelationalConstraint::StereoBondAnyLigand { stereo_bond, .. } => {
            entities.insert(Entity::StereoBond(*stereo_bond));
        }
    }
}

fn collect_molecule_constraint_entities(
    constraint: &MoleculeConstraint,
    entities: &mut BTreeSet<Entity>,
) {
    match constraint {
        MoleculeConstraint::ChargeSum { atoms, .. }
        | MoleculeConstraint::UnpairedElectronCoupling { atoms, .. }
        | MoleculeConstraint::Connected { atoms } => {
            if let Some(atoms) = atoms {
                entities.extend(atoms.iter().copied().map(Entity::Atom));
            }
        }
        MoleculeConstraint::BondOrderSum { bonds, .. } => {
            if let Some(bonds) = bonds {
                entities.extend(bonds.iter().copied().map(Entity::Bond));
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedAtom {
    pub id: AtomId,
    pub attributes: AtomForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedBond {
    pub id: BondId,
    pub endpoints: [AtomId; 2],
    pub attributes: BondForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedAtom {
    pub id: AtomId,
    pub attributes: AtomForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedBond {
    pub id: BondId,
    pub endpoints: [AtomId; 2],
    pub attributes: BondForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedDativeBond {
    pub id: DativeBondId,
    pub atoms: Vec<AtomId>,
    pub attributes: DativeBondForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedDativeBond {
    pub id: DativeBondId,
    pub atoms: Vec<AtomId>,
    pub attributes: DativeBondForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedAromaticSystem {
    pub id: AromaticSystemId,
    pub atoms: Vec<AtomId>,
    pub attributes: AromaticSystemForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedAromaticSystem {
    pub id: AromaticSystemId,
    pub atoms: Vec<AtomId>,
    pub attributes: AromaticSystemForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedMulticenterBond {
    pub id: MulticenterBondId,
    pub atoms: Vec<AtomId>,
    pub attributes: MulticenterBondForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedMulticenterBond {
    pub id: MulticenterBondId,
    pub atoms: Vec<AtomId>,
    pub attributes: MulticenterBondForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedNoncovalentBond {
    pub id: NoncovalentBondId,
    pub atoms: [AtomId; 2],
    pub attributes: NoncovalentBondForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedNoncovalentBond {
    pub id: NoncovalentBondId,
    pub atoms: [AtomId; 2],
    pub attributes: NoncovalentBondForm,
}

// Stereo elements carry both factors: the `site` (atom/bond) and the ordered
// `ligands`, unlike the single-atom-set overlays above.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedStereoAtom {
    pub id: StereoAtomId,
    pub site: AtomId,
    pub ligands: Vec<StereoLigand>,
    pub attributes: StereoAtomForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedStereoAtom {
    pub id: StereoAtomId,
    pub site: AtomId,
    pub ligands: Vec<StereoLigand>,
    pub attributes: StereoAtomForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedStereoBond {
    pub id: StereoBondId,
    pub site: BondId,
    pub ligands: Vec<StereoLigand>,
    pub attributes: StereoBondForm,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedStereoBond {
    pub id: StereoBondId,
    pub site: BondId,
    pub ligands: Vec<StereoLigand>,
    pub attributes: StereoBondForm,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemovedOverlays {
    pub dative_bonds: Vec<RemovedDativeBond>,
    pub aromatic_systems: Vec<RemovedAromaticSystem>,
    pub multicenter_bonds: Vec<RemovedMulticenterBond>,
    pub noncovalent_bonds: Vec<RemovedNoncovalentBond>,
    pub stereo_atoms: Vec<RemovedStereoAtom>,
    pub stereo_bonds: Vec<RemovedStereoBond>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedConstraint {
    pub position: usize,
    pub constraint: Constraint,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModifiedConstraint {
    pub position: usize,
    pub old: Constraint,
    pub new: Constraint,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CascadedConstraints {
    pub removed: Vec<RemovedConstraint>,
    pub modified: Vec<ModifiedConstraint>,
}

impl CascadedConstraints {
    pub fn is_empty(&self) -> bool {
        self.removed.is_empty() && self.modified.is_empty()
    }
}

/// The information needed to reverse one applied [`Edit`], captured as it was applied.
///
/// Checked application resolves symbolic handles and observes concrete allocations, compactions,
/// and cascaded removals that cannot be derived from the edit alone. Recorded undo entries are
/// replayed immediately if a later edit in the same batch fails. After successful application they
/// are returned in a [`Transaction`](crate::ir::Transaction) for explicit rollback against the
/// exact post-transaction state.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Undo {
    RemoveAddedTopology {
        atoms: Vec<AddedAtom>,
        bonds: Vec<AddedBond>,
    },
    RestoreRemovedTopology {
        atoms: Vec<RemovedAtom>,
        bonds: Vec<RemovedBond>,
        overlays: RemovedOverlays,
        compaction: IdCompaction,
        undo_compaction: UndoCompaction,
        cascade: CascadedConstraints,
    },
    RemoveAddedDativeBond(AddedDativeBond),
    RestoreRemovedDativeBonds {
        removed: Vec<RemovedDativeBond>,
        undo_compaction: UndoCompaction,
        cascade: CascadedConstraints,
    },
    RemoveAddedAromaticSystem(AddedAromaticSystem),
    RestoreRemovedAromaticSystems {
        removed: Vec<RemovedAromaticSystem>,
        undo_compaction: UndoCompaction,
        cascade: CascadedConstraints,
    },
    RemoveAddedMulticenterBond(AddedMulticenterBond),
    RestoreRemovedMulticenterBonds {
        removed: Vec<RemovedMulticenterBond>,
        undo_compaction: UndoCompaction,
        cascade: CascadedConstraints,
    },
    RemoveAddedNoncovalentBond(AddedNoncovalentBond),
    RestoreRemovedNoncovalentBonds {
        removed: Vec<RemovedNoncovalentBond>,
        undo_compaction: UndoCompaction,
        cascade: CascadedConstraints,
    },
    RemoveAddedStereoAtom(AddedStereoAtom),
    RestoreRemovedStereoAtoms {
        removed: Vec<RemovedStereoAtom>,
        undo_compaction: UndoCompaction,
        cascade: CascadedConstraints,
    },
    RemoveAddedStereoBond(AddedStereoBond),
    RestoreRemovedStereoBonds {
        removed: Vec<RemovedStereoBond>,
        undo_compaction: UndoCompaction,
        cascade: CascadedConstraints,
    },
    ModifyAtomField {
        id: AtomId,
        change: AtomFieldChange,
    },
    ModifyBondField {
        id: BondId,
        change: BondFieldChange,
    },
    ModifyDativeBondField {
        id: DativeBondId,
        change: DativeBondFieldChange,
    },
    ModifyAromaticSystemField {
        id: AromaticSystemId,
        change: AromaticSystemFieldChange,
    },
    ModifyMulticenterBondField {
        id: MulticenterBondId,
        change: MulticenterBondFieldChange,
    },
    ModifyNoncovalentBondField {
        id: NoncovalentBondId,
        change: NoncovalentBondFieldChange,
    },
    ModifyStereoAtomField {
        id: StereoAtomId,
        change: StereoAtomFieldChange,
    },
    ModifyStereoBondField {
        id: StereoBondId,
        change: StereoBondFieldChange,
    },
    ApplyCascadedConstraints(CascadedConstraints),
    ApplyEdit(Box<Edit>),
}

impl Undo {
    pub fn id_compaction(&self) -> Option<&IdCompaction> {
        match self {
            Self::RestoreRemovedTopology { compaction, .. } => Some(compaction),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use rstest::*;
    use umol_chem::element::Element;

    use super::super::boolean::BooleanForm;
    use super::super::constraint::{
        AromaticSystemConstraintsForm, AromaticValenceForm, AtomConstraintsForm,
        BondConstraintsForm, DativeBondConstraintsForm, MoleculeConstraint,
        MulticenterBondConstraintsForm, NoncovalentBondConstraintsForm, RelationalConstraint,
        RingScope, StereoAtomConstraintsForm, StereoBondConstraintsForm, StereogenicityForm,
    };
    use super::super::molecule::{Molecule, MoleculeEntries};
    use super::super::noncovalent::NoncovalentBondKind;
    use super::super::spin::UnpairedElectronsUpdate;
    use super::super::stereo::{
        StereoConfigurationForm, StereoConfigurationUpdate, StereoCoset, StereoKind, Stereogenicity,
    };
    use super::*;

    #[rstest]
    #[case::element(
        AtomFieldChange::Element {
            old: ElementForm::Lit(Element::C),
            new: ElementForm::Lit(Element::N),
        },
        AtomFieldChange::Element {
            old: ElementForm::Lit(Element::N),
            new: ElementForm::Lit(Element::C),
        },
    )]
    #[case::charge(
        AtomFieldChange::Charge {
            old: NumForm::Lit(0),
            new: NumForm::Lit(1),
        },
        AtomFieldChange::Charge {
            old: NumForm::Lit(1),
            new: NumForm::Lit(0),
        },
    )]
    fn test_atom_field_change_inverse(
        #[case] input: AtomFieldChange,
        #[case] expected: AtomFieldChange,
    ) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    #[case::configuration(
        StereoAtomFieldChange::Configuration {
            old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        },
        StereoAtomFieldChange::Configuration {
            old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        },
    )]
    fn test_stereo_atom_field_change_inverse(
        #[case] input: StereoAtomFieldChange,
        #[case] expected: StereoAtomFieldChange,
    ) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    #[case::configuration(
        StereoBondFieldChange::Configuration {
            old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
            new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
        },
        StereoBondFieldChange::Configuration {
            old: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
            new: StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
        },
    )]
    fn test_stereo_bond_field_change_inverse(
        #[case] input: StereoBondFieldChange,
        #[case] expected: StereoBondFieldChange,
    ) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    fn test_edits_new() {
        let edits = Edits::new();
        assert!(edits.is_empty());
        assert_eq!(
            edits,
            Edits {
                edits: Vec::new(),
                created_atoms: 0,
                created_bonds: 0,
                created_dative_bonds: 0,
                created_aromatic_systems: 0,
                created_multicenter_bonds: 0,
                created_noncovalent_bonds: 0,
                created_stereo_atoms: 0,
                created_stereo_bonds: 0,
            },
        );
    }

    #[rstest]
    fn test_edits_add_methods() {
        let atom = AtomForm::from_element(Element::C);
        let bond = BondForm::from_order(1);
        let dative = DativeBondForm::default();
        let aromatic = AromaticSystemForm::default();
        let multicenter = MulticenterBondForm::default();
        let noncovalent = NoncovalentBondForm::default();
        let stereo_atom = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0));
        let stereo_bond = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0));
        let mut edits = Edits::new();

        assert_eq!(edits.add_atom(atom.clone()), AtomHandle::New(0));
        assert_eq!(
            edits.add_dative_bond(vec![AtomHandle::Id(AtomId(0))], dative.clone()),
            DativeBondHandle::New(0),
        );
        assert_eq!(
            edits.add_bond(AtomHandle::Id(AtomId(0)), AtomHandle::New(0), bond.clone(),),
            BondHandle::New(0),
        );
        assert_eq!(
            edits.add_aromatic_system(vec![AtomHandle::New(0)], aromatic.clone()),
            AromaticSystemHandle::New(0),
        );
        assert_eq!(
            edits.add_multicenter_bond(vec![AtomHandle::New(0)], multicenter.clone()),
            MulticenterBondHandle::New(0),
        );
        assert_eq!(
            edits.add_noncovalent_bond(
                [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                noncovalent.clone(),
            ),
            NoncovalentBondHandle::New(0),
        );
        assert_eq!(
            edits.add_stereo_atom(
                AtomHandle::New(0),
                vec![(AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom)],
                stereo_atom.clone(),
            ),
            StereoAtomHandle::New(0),
        );
        assert_eq!(
            edits.add_stereo_bond(
                BondHandle::New(0),
                vec![(AtomHandle::New(0), StereoLigandKind::Atom)],
                stereo_bond.clone(),
            ),
            StereoBondHandle::New(0),
        );
        assert_eq!(
            edits.as_slice(),
            [
                Edit::AddAtoms { atoms: vec![atom] },
                Edit::AddDativeBond {
                    atoms: vec![AtomHandle::Id(AtomId(0))],
                    attributes: dative,
                },
                Edit::AddBonds {
                    bonds: vec![AddBond {
                        endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                        attributes: bond,
                    }],
                },
                Edit::AddAromaticSystem {
                    atoms: vec![AtomHandle::New(0)],
                    attributes: aromatic,
                },
                Edit::AddMulticenterBond {
                    atoms: vec![AtomHandle::New(0)],
                    attributes: multicenter,
                },
                Edit::AddNoncovalentBond {
                    atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                    attributes: noncovalent,
                },
                Edit::AddStereoAtom {
                    site: AtomHandle::New(0),
                    ligands: vec![(AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom)],
                    attributes: stereo_atom,
                },
                Edit::AddStereoBond {
                    site: BondHandle::New(0),
                    ligands: vec![(AtomHandle::New(0), StereoLigandKind::Atom)],
                    attributes: stereo_bond,
                },
            ],
        );
    }

    #[rstest]
    fn test_edits_add_atoms() {
        let carbon = AtomForm::from_element(Element::C);
        let nitrogen = AtomForm::from_element(Element::N);
        let oxygen = AtomForm::from_element(Element::O);
        let mut edits = Edits::new();

        assert_eq!(
            edits.add_atoms([carbon.clone(), nitrogen.clone()]),
            vec![AtomHandle::New(0), AtomHandle::New(1)],
        );
        assert_eq!(edits.add_atom(oxygen.clone()), AtomHandle::New(2));
        assert_eq!(
            edits.as_slice(),
            [
                Edit::AddAtoms {
                    atoms: vec![carbon, nitrogen],
                },
                Edit::AddAtoms {
                    atoms: vec![oxygen],
                },
            ],
        );
    }

    #[rstest]
    fn test_edits_add_bonds() {
        let single = BondForm::from_order(1);
        let double = BondForm::from_order(2);
        let triple = BondForm::from_order(3);
        let mut edits = Edits::new();
        let bonds = vec![
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                attributes: single.clone(),
            },
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
                attributes: double.clone(),
            },
        ];

        assert_eq!(
            edits.add_bonds(bonds.clone()),
            vec![BondHandle::New(0), BondHandle::New(1)],
        );
        assert_eq!(
            edits.add_bond(
                AtomHandle::Id(AtomId(2)),
                AtomHandle::Id(AtomId(3)),
                triple.clone(),
            ),
            BondHandle::New(2),
        );
        assert_eq!(
            edits.as_slice(),
            [
                Edit::AddBonds { bonds },
                Edit::AddBonds {
                    bonds: vec![AddBond {
                        endpoints: [AtomHandle::Id(AtomId(2)), AtomHandle::Id(AtomId(3)),],
                        attributes: triple,
                    }],
                },
            ],
        );
    }

    #[rstest]
    fn test_edits_add_overlay_batches() {
        let mut edits = Edits::new();

        assert_eq!(
            edits.add_dative_bonds([
                (vec![AtomHandle::Id(AtomId(0))], DativeBondForm::default()),
                (vec![AtomHandle::Id(AtomId(1))], DativeBondForm::default()),
            ]),
            vec![DativeBondHandle::New(0), DativeBondHandle::New(1)],
        );
        assert_eq!(
            edits.add_aromatic_systems([
                (
                    vec![AtomHandle::Id(AtomId(0))],
                    AromaticSystemForm::default(),
                ),
                (
                    vec![AtomHandle::Id(AtomId(1))],
                    AromaticSystemForm::default(),
                ),
            ]),
            vec![AromaticSystemHandle::New(0), AromaticSystemHandle::New(1)],
        );
        assert_eq!(
            edits.add_multicenter_bonds([
                (
                    vec![AtomHandle::Id(AtomId(0))],
                    MulticenterBondForm::default(),
                ),
                (
                    vec![AtomHandle::Id(AtomId(1))],
                    MulticenterBondForm::default(),
                ),
            ]),
            vec![MulticenterBondHandle::New(0), MulticenterBondHandle::New(1)],
        );
        assert_eq!(
            edits.add_noncovalent_bonds([
                (
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondForm::default(),
                ),
                (
                    [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
                    NoncovalentBondForm::default(),
                ),
            ]),
            vec![NoncovalentBondHandle::New(0), NoncovalentBondHandle::New(1),],
        );
        assert_eq!(
            edits.add_stereo_atoms([
                (
                    AtomHandle::Id(AtomId(0)),
                    Vec::new(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                ),
                (
                    AtomHandle::Id(AtomId(1)),
                    Vec::new(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                ),
            ]),
            vec![StereoAtomHandle::New(0), StereoAtomHandle::New(1)],
        );
        assert_eq!(
            edits.add_stereo_bonds([
                (
                    BondHandle::Id(BondId(0)),
                    Vec::new(),
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                ),
                (
                    BondHandle::Id(BondId(1)),
                    Vec::new(),
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                ),
            ]),
            vec![StereoBondHandle::New(0), StereoBondHandle::New(1)],
        );
        assert_eq!(edits.len(), 12);
        assert_eq!(
            edits.iter().cloned().collect::<Vec<_>>(),
            vec![
                Edit::AddDativeBond {
                    atoms: vec![AtomHandle::Id(AtomId(0))],
                    attributes: DativeBondForm::default(),
                },
                Edit::AddDativeBond {
                    atoms: vec![AtomHandle::Id(AtomId(1))],
                    attributes: DativeBondForm::default(),
                },
                Edit::AddAromaticSystem {
                    atoms: vec![AtomHandle::Id(AtomId(0))],
                    attributes: AromaticSystemForm::default(),
                },
                Edit::AddAromaticSystem {
                    atoms: vec![AtomHandle::Id(AtomId(1))],
                    attributes: AromaticSystemForm::default(),
                },
                Edit::AddMulticenterBond {
                    atoms: vec![AtomHandle::Id(AtomId(0))],
                    attributes: MulticenterBondForm::default(),
                },
                Edit::AddMulticenterBond {
                    atoms: vec![AtomHandle::Id(AtomId(1))],
                    attributes: MulticenterBondForm::default(),
                },
                Edit::AddNoncovalentBond {
                    atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    attributes: NoncovalentBondForm::default(),
                },
                Edit::AddNoncovalentBond {
                    atoms: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
                    attributes: NoncovalentBondForm::default(),
                },
                Edit::AddStereoAtom {
                    site: AtomHandle::Id(AtomId(0)),
                    ligands: Vec::new(),
                    attributes: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                },
                Edit::AddStereoAtom {
                    site: AtomHandle::Id(AtomId(1)),
                    ligands: Vec::new(),
                    attributes: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                },
                Edit::AddStereoBond {
                    site: BondHandle::Id(BondId(0)),
                    ligands: Vec::new(),
                    attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                },
                Edit::AddStereoBond {
                    site: BondHandle::Id(BondId(1)),
                    ligands: Vec::new(),
                    attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
                },
            ],
        );
    }

    #[rstest]
    fn test_edits_remove_topology() {
        let mut edits = Edits::new();
        edits.remove_atom(AtomHandle::Id(AtomId(0)));
        edits.remove_bond(BondHandle::New(0));
        edits.remove_topology(vec![AtomHandle::New(1)], vec![BondHandle::Id(BondId(2))]);

        assert_eq!(
            edits.as_slice(),
            [
                Edit::RemoveTopology {
                    atoms: vec![AtomHandle::Id(AtomId(0))],
                    bonds: Vec::new(),
                },
                Edit::RemoveTopology {
                    atoms: Vec::new(),
                    bonds: vec![BondHandle::New(0)],
                },
                Edit::RemoveTopology {
                    atoms: vec![AtomHandle::New(1)],
                    bonds: vec![BondHandle::Id(BondId(2))],
                },
            ],
        );
    }

    #[rstest]
    fn test_edits_remove_overlays() {
        let mut edits = Edits::new();
        edits.remove_dative_bonds(vec![(
            DativeBondHandle::Id(DativeBondId(0)),
            vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
            DativeBondForm::default(),
        )]);
        edits.remove_aromatic_systems(vec![(
            AromaticSystemHandle::New(0),
            vec![AtomHandle::Id(AtomId(0))],
            AromaticSystemForm::default(),
        )]);
        edits.remove_multicenter_bonds(vec![(
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            vec![AtomHandle::Id(AtomId(0))],
            MulticenterBondForm::default(),
        )]);
        edits.remove_noncovalent_bonds(vec![(
            NoncovalentBondHandle::New(0),
            [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
            NoncovalentBondForm::default(),
        )]);
        edits.remove_stereo_atoms(vec![(
            StereoAtomHandle::Id(StereoAtomId(0)),
            AtomHandle::Id(AtomId(0)),
            vec![(AtomHandle::New(0), StereoLigandKind::Atom)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )]);
        edits.remove_stereo_bonds(vec![(
            StereoBondHandle::New(0),
            BondHandle::Id(BondId(0)),
            vec![(AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom)],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )]);

        assert_eq!(
            edits.into_iter().collect::<Vec<_>>(),
            vec![
                Edit::RemoveDativeBonds {
                    removes: vec![(
                        DativeBondHandle::Id(DativeBondId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                        DativeBondForm::default(),
                    )],
                },
                Edit::RemoveAromaticSystems {
                    removes: vec![(
                        AromaticSystemHandle::New(0),
                        vec![AtomHandle::Id(AtomId(0))],
                        AromaticSystemForm::default(),
                    )],
                },
                Edit::RemoveMulticenterBonds {
                    removes: vec![(
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        vec![AtomHandle::Id(AtomId(0))],
                        MulticenterBondForm::default(),
                    )],
                },
                Edit::RemoveNoncovalentBonds {
                    removes: vec![(
                        NoncovalentBondHandle::New(0),
                        [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                        NoncovalentBondForm::default(),
                    )],
                },
                Edit::RemoveStereoAtoms {
                    removes: vec![(
                        StereoAtomHandle::Id(StereoAtomId(0)),
                        AtomHandle::Id(AtomId(0)),
                        vec![(AtomHandle::New(0), StereoLigandKind::Atom)],
                        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                    )],
                },
                Edit::RemoveStereoBonds {
                    removes: vec![(
                        StereoBondHandle::New(0),
                        BondHandle::Id(BondId(0)),
                        vec![(AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom)],
                        StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    )],
                },
            ],
        );
    }

    #[rstest]
    fn test_edits_molecule_constraint() {
        let constraint = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        let mut edits = Edits::new();
        edits.add_molecule_constraint(constraint.clone().into());
        edits.remove_molecule_constraint(constraint.clone().into());

        assert_eq!(
            edits.as_slice(),
            [
                Edit::AddMoleculeConstraint {
                    constraint: constraint.clone().into(),
                },
                Edit::RemoveMoleculeConstraint {
                    constraint: constraint.into(),
                },
            ],
        );
    }

    #[rstest]
    fn test_edits_push() {
        let entry = Edit::AddAtoms {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::N),
            ],
        };
        let mut edits = Edits::new();
        edits.push(entry.clone());

        assert_eq!(
            edits.add_atom(AtomForm::from_element(Element::O)),
            AtomHandle::New(2)
        );
        assert_eq!(edits.as_slice()[0], entry);
    }

    #[rstest]
    fn test_edits_from_iter() {
        let entries = vec![
            Edit::AddAtoms {
                atoms: vec![
                    AtomForm::from_element(Element::C),
                    AtomForm::from_element(Element::N),
                ],
            },
            Edit::AddBonds {
                bonds: vec![
                    AddBond {
                        endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        attributes: BondForm::from_order(1),
                    },
                    AddBond {
                        endpoints: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
                        attributes: BondForm::from_order(1),
                    },
                ],
            },
            Edit::AddDativeBond {
                atoms: Vec::new(),
                attributes: DativeBondForm::default(),
            },
            Edit::AddAromaticSystem {
                atoms: Vec::new(),
                attributes: AromaticSystemForm::default(),
            },
            Edit::AddMulticenterBond {
                atoms: Vec::new(),
                attributes: MulticenterBondForm::default(),
            },
            Edit::AddNoncovalentBond {
                atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                attributes: NoncovalentBondForm::default(),
            },
            Edit::AddStereoAtom {
                site: AtomHandle::Id(AtomId(0)),
                ligands: Vec::new(),
                attributes: StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            },
            Edit::AddStereoBond {
                site: BondHandle::Id(BondId(0)),
                ligands: Vec::new(),
                attributes: StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            },
        ];
        let mut edits: Edits = entries.clone().into_iter().collect();

        assert_eq!(edits.as_slice(), entries);
        assert_eq!(edits.add_atom(AtomForm::default()), AtomHandle::New(2));
        assert_eq!(
            edits.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
                BondForm::default(),
            ),
            BondHandle::New(2),
        );
        assert_eq!(
            edits.add_dative_bond(Vec::new(), DativeBondForm::default()),
            DativeBondHandle::New(1),
        );
        assert_eq!(
            edits.add_aromatic_system(Vec::new(), AromaticSystemForm::default()),
            AromaticSystemHandle::New(1),
        );
        assert_eq!(
            edits.add_multicenter_bond(Vec::new(), MulticenterBondForm::default()),
            MulticenterBondHandle::New(1),
        );
        assert_eq!(
            edits.add_noncovalent_bond(
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                NoncovalentBondForm::default(),
            ),
            NoncovalentBondHandle::New(1),
        );
        assert_eq!(
            edits.add_stereo_atom(
                AtomHandle::Id(AtomId(0)),
                Vec::new(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            StereoAtomHandle::New(1),
        );
        assert_eq!(
            edits.add_stereo_bond(
                BondHandle::Id(BondId(0)),
                Vec::new(),
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
            StereoBondHandle::New(1),
        );
    }

    #[rstest]
    fn test_edits_iter() {
        let entries = vec![
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(0))],
                bonds: Vec::new(),
            },
            Edit::AddMoleculeConstraint {
                constraint: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None })
                    .into(),
            },
        ];
        let edits: Edits = entries.clone().into_iter().collect();

        assert!(!edits.is_empty());
        assert_eq!(edits.len(), 2);
        assert_eq!(edits.as_slice(), entries);
        assert_eq!(edits.iter().cloned().collect::<Vec<_>>(), entries);
        assert_eq!(edits.into_iter().collect::<Vec<_>>(), entries);
    }

    #[rstest]
    fn test_edits_update_atom() {
        let current = AtomForm::from_element(Element::C)
            .with_isotope_mass(12_u32)
            .with_charge(0_i64)
            .with_implicit_hydrogens(4_i64)
            .with_lone_pairs(0_i64)
            .with_unpaired_electrons((2_u8, 3_u8))
            .with_constraint(AtomConstraintForm::valence(4_i64));
        let update = AtomUpdate {
            element: Some(ElementForm::Lit(Element::N)),
            isotope_mass: Some(IsotopeMassForm::Lit(13)),
            charge: Some(NumForm::Lit(1)),
            implicit_hydrogens: Some(NumForm::Lit(3)),
            lone_pairs: Some(NumForm::Lit(1)),
            unpaired_electrons: UnpairedElectronsUpdate {
                count: None,
                multiplicity: Some(NumForm::Lit(1)),
            },
            constraints: AtomConstraintsForm::from_iter([
                AtomConstraintForm::valence(NumForm::Undetermined),
                AtomConstraintForm::degree(2_i64),
            ]),
        };
        let mut edits = Edits::new();
        edits.update_atom(AtomHandle::Id(AtomId(7)), &current, &update);

        assert_eq!(
            edits.as_slice(),
            [
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::Element {
                        old: ElementForm::Lit(Element::C),
                        new: ElementForm::Lit(Element::N),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::IsotopeMass {
                        old: IsotopeMassForm::Lit(12),
                        new: IsotopeMassForm::Lit(13),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::Charge {
                        old: NumForm::Lit(0),
                        new: NumForm::Lit(1),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::ImplicitHydrogens {
                        old: NumForm::Lit(4),
                        new: NumForm::Lit(3),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::LonePairs {
                        old: NumForm::Lit(0),
                        new: NumForm::Lit(1),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::UnpairedElectrons {
                        old: UnpairedElectronsForm::from((2_u8, 3_u8)),
                        new: UnpairedElectronsForm::from((2_u8, 1_u8)),
                    },
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(7)),
                    old: Some(AtomConstraintForm::valence(4_i64)),
                    new: None,
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(7)),
                    old: None,
                    new: Some(AtomConstraintForm::degree(2_i64)),
                },
            ]
        );

        let expected = current.update(&update);
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![current.clone()],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_atom(AtomHandle::Id(AtomId(0)), &current, &update);
        editor
            .transact(applied_edits)
            .expect("atom update edits should apply");

        assert_eq!(editor.atom(AtomId(0)).attributes, &expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AtomForm::from_element(Element::C), AtomUpdate::default())]
    #[case::canonical_field(AtomForm::from_element(Element::C).with_charge(1_i64), AtomUpdate { charge: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(AtomForm::from_element(Element::C), AtomUpdate { constraints: AtomConstraintsForm::from(AtomConstraintForm::valence(NumForm::Undetermined)), ..Default::default() })]
    fn test_edits_update_atom_identity(#[case] current: AtomForm, #[case] update: AtomUpdate) {
        let mut edits = Edits::new();
        edits.update_atom(AtomHandle::Id(AtomId(0)), &current, &update);

        assert_eq!(edits, Edits::new());
    }

    #[rstest]
    fn test_edits_update_bond() {
        let current = BondForm::from_order(1)
            .with_charge(0_i64)
            .with_unpaired_electrons((2_u8, 3_u8))
            .with_constraint(BondConstraintForm::ring_membership(
                RingScope::Size(6),
                1_i64,
            ));
        let update = BondUpdate {
            order: Some(NumForm::Lit(2)),
            charge: Some(NumForm::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate {
                count: None,
                multiplicity: Some(NumForm::Lit(1)),
            },
            constraints: BondConstraintsForm::from_iter([
                BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined),
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            ]),
        };
        let mut edits = Edits::new();
        edits.update_bond(BondHandle::Id(BondId(7)), &current, &update);

        assert_eq!(
            edits.as_slice(),
            [
                Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(7)),
                    change: BondFieldChange::Order {
                        old: NumForm::Lit(1),
                        new: NumForm::Lit(2),
                    },
                },
                Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(7)),
                    change: BondFieldChange::Charge {
                        old: NumForm::Lit(0),
                        new: NumForm::Undetermined,
                    },
                },
                Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(7)),
                    change: BondFieldChange::UnpairedElectrons {
                        old: UnpairedElectronsForm::from((2_u8, 3_u8)),
                        new: UnpairedElectronsForm::from((2_u8, 1_u8)),
                    },
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(7)),
                    old: None,
                    new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(7)),
                    old: Some(BondConstraintForm::ring_membership(
                        RingScope::Size(6),
                        1_i64,
                    )),
                    new: None,
                },
            ]
        );

        let expected = current.update(&update);
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(), AtomForm::default()],
            bonds: vec![(AtomId(0), AtomId(1), current.clone())],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_bond(BondHandle::Id(BondId(0)), &current, &update);
        editor
            .transact(applied_edits)
            .expect("bond update edits should apply");

        assert_eq!(editor.bond(BondId(0)).attributes, &expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(BondForm::from_order(1), BondUpdate::default())]
    #[case::canonical_field(BondForm::from_order(1).with_charge(1_i64), BondUpdate { charge: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(BondForm::from_order(1), BondUpdate { constraints: BondConstraintsForm::from(BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() })]
    fn test_edits_update_bond_identity(#[case] current: BondForm, #[case] update: BondUpdate) {
        let mut edits = Edits::new();
        edits.update_bond(BondHandle::Id(BondId(0)), &current, &update);

        assert_eq!(edits, Edits::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)),
        DativeBondUpdate {
            order: Some(NumForm::Lit(2)),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined),
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            ]),
        },
        vec![
            Edit::ModifyDativeBondField {
                id: DativeBondHandle::Id(DativeBondId(7)),
                change: DativeBondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            },
            Edit::ModifyDativeBondConstraint {
                id: DativeBondHandle::Id(DativeBondId(7)),
                old: None,
                new: Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            },
            Edit::ModifyDativeBondConstraint {
                id: DativeBondHandle::Id(DativeBondId(7)),
                old: Some(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)),
                new: None,
            },
        ],
    )]
    fn test_edits_update_dative_bond(
        #[case] current: DativeBondForm,
        #[case] update: DativeBondUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        let mut edits = Edits::new();
        edits.update_dative_bond(
            DativeBondHandle::Id(DativeBondId(7)),
            &current,
            &update,
        );

        assert_eq!(edits.as_slice(), expected);

        let expected_attributes = current.update(&update);
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(), AtomForm::default()],
            dative: vec![(vec![AtomId(0)], AtomId(1), current.clone())],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_dative_bond(
            DativeBondHandle::Id(DativeBondId(0)),
            &current,
            &update,
        );
        editor
            .transact(applied_edits)
            .expect("dative-bond update edits should apply");

        assert_eq!(
            editor.dative_bond(DativeBondId(0)).attributes,
            &expected_attributes
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(DativeBondForm::from_order(1), DativeBondUpdate::default())]
    #[case::canonical_field(DativeBondForm::from_order(1), DativeBondUpdate { order: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(DativeBondForm::from_order(1), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() })]
    fn test_edits_update_dative_bond_identity(
        #[case] current: DativeBondForm,
        #[case] update: DativeBondUpdate,
    ) {
        let mut edits = Edits::new();
        edits.update_dative_bond(
            DativeBondHandle::Id(DativeBondId(0)),
            &current,
            &update,
        );

        assert_eq!(edits, Edits::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraint(
        AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(AromaticSystemConstraintForm::electron_count(6_i64)),
        AromaticSystemUpdate {
            electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])),
            charge: Some(NumForm::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) },
            constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)),
        },
        vec![
            Edit::ModifyAromaticSystemField {
                id: AromaticSystemHandle::Id(AromaticSystemId(7)),
                change: AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 1, 1]), new: ElectronCountsForm::Lit(vec![2, 2, 2]) },
            },
            Edit::ModifyAromaticSystemField {
                id: AromaticSystemHandle::Id(AromaticSystemId(7)),
                change: AromaticSystemFieldChange::Charge { old: NumForm::Lit(0), new: NumForm::Undetermined },
            },
            Edit::ModifyAromaticSystemField {
                id: AromaticSystemHandle::Id(AromaticSystemId(7)),
                change: AromaticSystemFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::from((2_u8, 3_u8)), new: UnpairedElectronsForm::from((2_u8, 1_u8)) },
            },
            Edit::ModifyAromaticSystemConstraint {
                id: AromaticSystemHandle::Id(AromaticSystemId(7)),
                old: Some(AromaticSystemConstraintForm::electron_count(6_i64)),
                new: None,
            },
        ],
    )]
    fn test_edits_update_aromatic_system(
        #[case] current: AromaticSystemForm,
        #[case] update: AromaticSystemUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        let mut edits = Edits::new();
        edits.update_aromatic_system(
            AromaticSystemHandle::Id(AromaticSystemId(7)),
            &current,
            &update,
        );

        assert_eq!(edits.as_slice(), expected);

        let expected_attributes = current.update(&update);
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::default(),
                AtomForm::default(),
                AtomForm::default(),
            ],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], current.clone())],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_aromatic_system(
            AromaticSystemHandle::Id(AromaticSystemId(0)),
            &current,
            &update,
        );
        editor
            .transact(applied_edits)
            .expect("aromatic-system update edits should apply");

        assert_eq!(
            editor.aromatic_system(AromaticSystemId(0)).attributes,
            &expected_attributes,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate::default())]
    #[case::canonical_field(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(1_i64), AromaticSystemUpdate { charge: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() })]
    fn test_edits_update_aromatic_system_identity(
        #[case] current: AromaticSystemForm,
        #[case] update: AromaticSystemUpdate,
    ) {
        let mut edits = Edits::new();
        edits.update_aromatic_system(
            AromaticSystemHandle::Id(AromaticSystemId(0)),
            &current,
            &update,
        );

        assert_eq!(edits, Edits::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraint(
        MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(MulticenterBondConstraintForm::electron_count(6_i64)),
        MulticenterBondUpdate {
            electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])),
            charge: Some(NumForm::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) },
            constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(NumForm::Undetermined)),
        },
        vec![
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(7)),
                change: MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 1, 1]), new: ElectronCountsForm::Lit(vec![2, 2, 2]) },
            },
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(7)),
                change: MulticenterBondFieldChange::Charge { old: NumForm::Lit(0), new: NumForm::Undetermined },
            },
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(7)),
                change: MulticenterBondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::from((2_u8, 3_u8)), new: UnpairedElectronsForm::from((2_u8, 1_u8)) },
            },
            Edit::ModifyMulticenterBondConstraint {
                id: MulticenterBondHandle::Id(MulticenterBondId(7)),
                old: Some(MulticenterBondConstraintForm::electron_count(6_i64)),
                new: None,
            },
        ],
    )]
    fn test_edits_update_multicenter_bond(
        #[case] current: MulticenterBondForm,
        #[case] update: MulticenterBondUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        let mut edits = Edits::new();
        edits.update_multicenter_bond(
            MulticenterBondHandle::Id(MulticenterBondId(7)),
            &current,
            &update,
        );

        assert_eq!(edits.as_slice(), expected);

        let expected_attributes = current.update(&update);
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::default(),
                AtomForm::default(),
                AtomForm::default(),
            ],
            multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], current.clone())],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_multicenter_bond(
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            &current,
            &update,
        );
        editor
            .transact(applied_edits)
            .expect("multicenter-bond update edits should apply");

        assert_eq!(
            editor.multicenter_bond(MulticenterBondId(0)).attributes,
            &expected_attributes,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate::default())]
    #[case::canonical_field(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(1_i64), MulticenterBondUpdate { charge: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() })]
    fn test_edits_update_multicenter_bond_identity(
        #[case] current: MulticenterBondForm,
        #[case] update: MulticenterBondUpdate,
    ) {
        let mut edits = Edits::new();
        edits.update_multicenter_bond(
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            &current,
            &update,
        );

        assert_eq!(edits, Edits::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::kind_and_constraint(
        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintForm::intramolecular(true)),
        NoncovalentBondUpdate {
            kind: Some(NoncovalentBondKindForm::Undetermined),
            constraints: NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(BooleanForm::Undetermined)),
        },
        vec![
            Edit::ModifyNoncovalentBondField {
                id: NoncovalentBondHandle::Id(NoncovalentBondId(7)),
                change: NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindForm::Undetermined },
            },
            Edit::ModifyNoncovalentBondConstraint {
                id: NoncovalentBondHandle::Id(NoncovalentBondId(7)),
                old: Some(NoncovalentBondConstraintForm::intramolecular(true)),
                new: None,
            },
        ],
    )]
    fn test_edits_update_noncovalent_bond(
        #[case] current: NoncovalentBondForm,
        #[case] update: NoncovalentBondUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        let mut edits = Edits::new();
        edits.update_noncovalent_bond(
            NoncovalentBondHandle::Id(NoncovalentBondId(7)),
            &current,
            &update,
        );

        assert_eq!(edits.as_slice(), expected);

        let expected_attributes = current.update(&update);
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(), AtomForm::default()],
            noncovalent: vec![(AtomId(0), AtomId(1), current.clone())],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_noncovalent_bond(
            NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            &current,
            &update,
        );
        editor
            .transact(applied_edits)
            .expect("noncovalent-bond update edits should apply");

        assert_eq!(
            editor.noncovalent_bond(NoncovalentBondId(0)).attributes,
            &expected_attributes,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate::default())]
    #[case::same_kind(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate { kind: Some(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)), ..Default::default() })]
    #[case::absent_constraint_removal(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(BooleanForm::Undetermined)), ..Default::default() })]
    fn test_edits_update_noncovalent_bond_identity(
        #[case] current: NoncovalentBondForm,
        #[case] update: NoncovalentBondUpdate,
    ) {
        let mut edits = Edits::new();
        edits.update_noncovalent_bond(
            NoncovalentBondHandle::Id(NoncovalentBondId(0)),
            &current,
            &update,
        );

        assert_eq!(edits, Edits::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::configuration_and_constraint(
        StereoAtomForm { configuration: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32), constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoAtomUpdate {
            configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) },
            constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        },
        vec![
            Edit::ModifyStereoAtomField {
                id: StereoAtomHandle::Id(StereoAtomId(7)),
                change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32) },
            },
            Edit::ModifyStereoAtomConstraint {
                id: StereoAtomHandle::Id(StereoAtomId(7)),
                kind: Some(StereoKind::Tetrahedral),
                old: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
                new: None,
            },
        ],
    )]
    fn test_edits_update_stereo_atom(
        #[case] current: StereoAtomForm,
        #[case] update: StereoAtomUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        let mut edits = Edits::new();
        edits.update_stereo_atom(
            StereoAtomHandle::Id(StereoAtomId(7)),
            &current,
            &update,
        );

        assert_eq!(edits.as_slice(), expected);

        let expected_attributes = current.update(&update);
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 5],
            stereo_atoms: vec![(
                AtomId(0),
                (1..=4)
                    .map(|index| StereoLigand::new(AtomId(index), StereoLigandKind::Atom))
                    .collect(),
                current.clone(),
            )],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_stereo_atom(
            StereoAtomHandle::Id(StereoAtomId(0)),
            &current,
            &update,
        );
        editor
            .transact(applied_edits)
            .expect("stereo-atom update edits should apply");

        assert_eq!(
            editor.stereo_atom(StereoAtomId(0)).attributes,
            &expected_attributes
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate::default())]
    #[case::relative(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() })]
    #[case::absent_constraint_removal(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() })]
    fn test_edits_update_stereo_atom_identity(
        #[case] current: StereoAtomForm,
        #[case] update: StereoAtomUpdate,
    ) {
        let mut edits = Edits::new();
        edits.update_stereo_atom(
            StereoAtomHandle::Id(StereoAtomId(0)),
            &current,
            &update,
        );

        assert_eq!(edits, Edits::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::configuration_and_constraint(
        StereoBondForm { configuration: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32), constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))) },
        StereoBondUpdate {
            configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) },
            constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        },
        vec![
            Edit::ModifyStereoBondField {
                id: StereoBondHandle::Id(StereoBondId(7)),
                change: StereoBondFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32) },
            },
            Edit::ModifyStereoBondConstraint {
                id: StereoBondHandle::Id(StereoBondId(7)),
                kind: Some(StereoKind::CisTrans),
                old: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
                new: None,
            },
        ],
    )]
    fn test_edits_update_stereo_bond(
        #[case] current: StereoBondForm,
        #[case] update: StereoBondUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        let mut edits = Edits::new();
        edits.update_stereo_bond(
            StereoBondHandle::Id(StereoBondId(7)),
            &current,
            &update,
        );

        assert_eq!(edits.as_slice(), expected);

        let expected_attributes = current.update(&update);
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::default(); 6],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            stereo_bonds: vec![(
                BondId(0),
                (2..=5)
                    .map(|index| StereoLigand::new(AtomId(index), StereoLigandKind::Atom))
                    .collect(),
                current.clone(),
            )],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_stereo_bond(
            StereoBondHandle::Id(StereoBondId(0)),
            &current,
            &update,
        );
        editor
            .transact(applied_edits)
            .expect("stereo-bond update edits should apply");

        assert_eq!(
            editor.stereo_bond(StereoBondId(0)).attributes,
            &expected_attributes
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoBondForm::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate::default())]
    #[case::relative(StereoBondForm::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() })]
    #[case::absent_constraint_removal(StereoBondForm::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() })]
    fn test_edits_update_stereo_bond_identity(
        #[case] current: StereoBondForm,
        #[case] update: StereoBondUpdate,
    ) {
        let mut edits = Edits::new();
        edits.update_stereo_bond(
            StereoBondHandle::Id(StereoBondId(0)),
            &current,
            &update,
        );

        assert_eq!(edits, Edits::new());
    }

    #[rstest]
    fn test_bond_field_change_inverse() {
        let change = BondFieldChange::Order {
            old: NumForm::Lit(1),
            new: NumForm::Lit(2),
        };
        assert_eq!(
            change.clone().inverse(),
            BondFieldChange::Order {
                old: NumForm::Lit(2),
                new: NumForm::Lit(1),
            },
        );
        assert_eq!(change.clone().inverse().inverse(), change);
    }

    #[rstest]
    #[case::atom(Entity::Atom(AtomId(1)), EntityHandle::Atom(AtomHandle::Id(AtomId(1))))]
    #[case::bond(Entity::Bond(BondId(2)), EntityHandle::Bond(BondHandle::Id(BondId(2))))]
    #[case::dative_bond(
        Entity::DativeBond(DativeBondId(3)),
        EntityHandle::DativeBond(DativeBondHandle::Id(DativeBondId(3)))
    )]
    #[case::aromatic_system(
        Entity::AromaticSystem(AromaticSystemId(4)),
        EntityHandle::AromaticSystem(AromaticSystemHandle::Id(AromaticSystemId(4)))
    )]
    #[case::multicenter_bond(
        Entity::MulticenterBond(MulticenterBondId(5)),
        EntityHandle::MulticenterBond(MulticenterBondHandle::Id(MulticenterBondId(5)))
    )]
    #[case::noncovalent_bond(
        Entity::NoncovalentBond(NoncovalentBondId(6)),
        EntityHandle::NoncovalentBond(NoncovalentBondHandle::Id(NoncovalentBondId(6)))
    )]
    #[case::stereo_atom(
        Entity::StereoAtom(StereoAtomId(7)),
        EntityHandle::StereoAtom(StereoAtomHandle::Id(StereoAtomId(7)))
    )]
    #[case::stereo_bond(
        Entity::StereoBond(StereoBondId(8)),
        EntityHandle::StereoBond(StereoBondHandle::Id(StereoBondId(8)))
    )]
    fn test_entity_handle_from(#[case] input: Entity, #[case] expected: EntityHandle) {
        assert_eq!(EntityHandle::from(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_leaf(
        Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
        vec![(Entity::Atom(AtomId(7)), EntityHandle::Atom(AtomHandle::New(2)))],
        ConstraintEdit {
            constraint: Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3_i64)),
            atoms: vec![AtomHandle::New(2)], bonds: vec![], dative_bonds: vec![], aromatic_systems: vec![],
            multicenter_bonds: vec![], noncovalent_bonds: vec![], stereo_atoms: vec![], stereo_bonds: vec![],
        },
    )]
    #[case::logical_shared_handle(
        Constraint::Or(vec![
            Constraint::And(vec![
                Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
                Constraint::Atom(AtomId(7), AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(2)))),
            ]),
            Constraint::And(vec![
                Constraint::Atom(AtomId(9), AtomConstraintForm::valence(2_i64)),
                Constraint::Not(Box::new(Constraint::Atom(AtomId(9), AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(1)))))),
            ]),
        ]),
        vec![
            (Entity::Atom(AtomId(7)), EntityHandle::Atom(AtomHandle::New(0))),
            (Entity::Atom(AtomId(9)), EntityHandle::Atom(AtomHandle::New(0))),
        ],
        ConstraintEdit {
            constraint: Constraint::Or(vec![
                Constraint::And(vec![
                    Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3_i64)),
                    Constraint::Atom(AtomId(0), AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(2)))),
                ]),
                Constraint::And(vec![
                    Constraint::Atom(AtomId(0), AtomConstraintForm::valence(2_i64)),
                    Constraint::Not(Box::new(Constraint::Atom(AtomId(0), AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(1)))))),
                ]),
            ]),
            atoms: vec![AtomHandle::New(0)], bonds: vec![], dative_bonds: vec![], aromatic_systems: vec![],
            multicenter_bonds: vec![], noncovalent_bonds: vec![], stereo_atoms: vec![], stereo_bonds: vec![],
        },
    )]
    #[case::relational_explicit(
        Constraint::Relational(RelationalConstraint::DativeBondParallels {
            dative: DativeBondId(5), parallel: BondId(8),
        }),
        vec![
            (Entity::DativeBond(DativeBondId(5)), EntityHandle::DativeBond(DativeBondHandle::New(1))),
            (Entity::Bond(BondId(8)), EntityHandle::Bond(BondHandle::Id(BondId(3)))),
        ],
        ConstraintEdit {
            constraint: Constraint::Relational(RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(0), parallel: BondId(0),
            }),
            atoms: vec![], bonds: vec![BondHandle::Id(BondId(3))], dative_bonds: vec![DativeBondHandle::New(1)], aromatic_systems: vec![],
            multicenter_bonds: vec![], noncovalent_bonds: vec![], stereo_atoms: vec![], stereo_bonds: vec![],
        },
    )]
    #[case::relational_quantified(
        Constraint::Relational(RelationalConstraint::DativeBondAllDonors {
            bond: DativeBondId(5), predicate: Box::new(AtomConstraintForm::valence(3_i64)),
        }),
        vec![(Entity::DativeBond(DativeBondId(5)), EntityHandle::DativeBond(DativeBondHandle::New(1)))],
        ConstraintEdit {
            constraint: Constraint::Relational(RelationalConstraint::DativeBondAllDonors {
                bond: DativeBondId(0), predicate: Box::new(AtomConstraintForm::valence(3_i64)),
            }),
            atoms: vec![], bonds: vec![], dative_bonds: vec![DativeBondHandle::New(1)], aromatic_systems: vec![],
            multicenter_bonds: vec![], noncovalent_bonds: vec![], stereo_atoms: vec![], stereo_bonds: vec![],
        },
    )]
    #[case::all_entity_leaves(
        Constraint::And(vec![
            Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
            Constraint::Bond(BondId(8), BondConstraintForm::aromatic(true)),
            Constraint::DativeBond(DativeBondId(9), DativeBondConstraintForm::aromatic(true)),
            Constraint::AromaticSystem(AromaticSystemId(10), AromaticSystemConstraintForm::electron_count(6_i64)),
            Constraint::MulticenterBond(MulticenterBondId(11), MulticenterBondConstraintForm::electron_count(2_i64)),
            Constraint::NoncovalentBond(NoncovalentBondId(12), NoncovalentBondConstraintForm::intramolecular(true)),
            Constraint::StereoAtom(StereoAtomId(13), StereoKind::Tetrahedral, StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
            Constraint::StereoBond(StereoBondId(14), StereoKind::CisTrans, StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        ]),
        vec![
            (Entity::Atom(AtomId(7)), EntityHandle::Atom(AtomHandle::New(0))),
            (Entity::Bond(BondId(8)), EntityHandle::Bond(BondHandle::New(1))),
            (Entity::DativeBond(DativeBondId(9)), EntityHandle::DativeBond(DativeBondHandle::New(2))),
            (Entity::AromaticSystem(AromaticSystemId(10)), EntityHandle::AromaticSystem(AromaticSystemHandle::New(3))),
            (Entity::MulticenterBond(MulticenterBondId(11)), EntityHandle::MulticenterBond(MulticenterBondHandle::New(4))),
            (Entity::NoncovalentBond(NoncovalentBondId(12)), EntityHandle::NoncovalentBond(NoncovalentBondHandle::New(5))),
            (Entity::StereoAtom(StereoAtomId(13)), EntityHandle::StereoAtom(StereoAtomHandle::New(6))),
            (Entity::StereoBond(StereoBondId(14)), EntityHandle::StereoBond(StereoBondHandle::New(7))),
        ],
        ConstraintEdit {
            constraint: Constraint::And(vec![
                Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3_i64)),
                Constraint::Bond(BondId(0), BondConstraintForm::aromatic(true)),
                Constraint::DativeBond(DativeBondId(0), DativeBondConstraintForm::aromatic(true)),
                Constraint::AromaticSystem(AromaticSystemId(0), AromaticSystemConstraintForm::electron_count(6_i64)),
                Constraint::MulticenterBond(MulticenterBondId(0), MulticenterBondConstraintForm::electron_count(2_i64)),
                Constraint::NoncovalentBond(NoncovalentBondId(0), NoncovalentBondConstraintForm::intramolecular(true)),
                Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
                Constraint::StereoBond(StereoBondId(0), StereoKind::CisTrans, StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
            ]),
            atoms: vec![AtomHandle::New(0)], bonds: vec![BondHandle::New(1)], dative_bonds: vec![DativeBondHandle::New(2)], aromatic_systems: vec![AromaticSystemHandle::New(3)],
            multicenter_bonds: vec![MulticenterBondHandle::New(4)], noncovalent_bonds: vec![NoncovalentBondHandle::New(5)], stereo_atoms: vec![StereoAtomHandle::New(6)], stereo_bonds: vec![StereoBondHandle::New(7)],
        },
    )]
    #[case::molecule_subsets(
        Constraint::And(vec![
            Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(7), AtomId(9)]), sum: NumForm::Lit(0) }),
            Constraint::Molecule(MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(4)]), sum: NumForm::Lit(2) }),
        ]),
        vec![
            (Entity::Atom(AtomId(7)), EntityHandle::Atom(AtomHandle::New(0))),
            (Entity::Atom(AtomId(9)), EntityHandle::Atom(AtomHandle::Id(AtomId(2)))),
            (Entity::Bond(BondId(4)), EntityHandle::Bond(BondHandle::New(1))),
        ],
        ConstraintEdit {
            constraint: Constraint::And(vec![
                Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: Some(vec![AtomId(0), AtomId(1)]), sum: NumForm::Lit(0) }),
                Constraint::Molecule(MoleculeConstraint::BondOrderSum { bonds: Some(vec![BondId(0)]), sum: NumForm::Lit(2) }),
            ]),
            atoms: vec![AtomHandle::New(0), AtomHandle::Id(AtomId(2))], bonds: vec![BondHandle::New(1)], dative_bonds: vec![], aromatic_systems: vec![],
            multicenter_bonds: vec![], noncovalent_bonds: vec![], stereo_atoms: vec![], stereo_bonds: vec![],
        },
    )]
    #[case::molecule_all_atoms(
        Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
        vec![],
        ConstraintEdit {
            constraint: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
            atoms: vec![], bonds: vec![], dative_bonds: vec![], aromatic_systems: vec![],
            multicenter_bonds: vec![], noncovalent_bonds: vec![], stereo_atoms: vec![], stereo_bonds: vec![],
        },
    )]
    fn test_constraint_edit_new(
        #[case] input: Constraint,
        #[case] mappings: Vec<(Entity, EntityHandle)>,
        #[case] expected: ConstraintEdit,
    ) {
        let mappings: HashMap<_, _> = mappings.into_iter().collect();

        assert_eq!(
            ConstraintEdit::new(input, |entity| Some(mappings[&entity].clone())),
            Ok(expected),
        );
    }

    #[rstest]
    #[case::atom_as_bond(
        Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
        Some(EntityHandle::Bond(BondHandle::New(0))),
        ConstraintEditError::HandleKindMismatch {
            entity: Entity::Atom(AtomId(7)),
            actual: EntityKind::Bond,
        },
    )]
    #[case::missing(
        Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
        None,
        ConstraintEditError::MissingHandle {
            entity: Entity::Atom(AtomId(7)),
        },
    )]
    fn test_constraint_edit_new_error(
        #[case] input: Constraint,
        #[case] handle: Option<EntityHandle>,
        #[case] expected: ConstraintEditError,
    ) {
        assert_eq!(
            ConstraintEdit::new(input, |_| handle.clone()),
            Err(expected)
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom_leaf(
        Constraint::Atom(AtomId(7), AtomConstraintForm::valence(3_i64)),
        ConstraintEdit {
            constraint: Constraint::Atom(AtomId(0), AtomConstraintForm::valence(3_i64)),
            atoms: vec![AtomHandle::Id(AtomId(7))], bonds: vec![], dative_bonds: vec![], aromatic_systems: vec![],
            multicenter_bonds: vec![], noncovalent_bonds: vec![], stereo_atoms: vec![], stereo_bonds: vec![],
        },
    )]
    #[case::relational(
        Constraint::Relational(RelationalConstraint::DativeBondParallels {
            dative: DativeBondId(5), parallel: BondId(8),
        }),
        ConstraintEdit {
            constraint: Constraint::Relational(RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(0), parallel: BondId(0),
            }),
            atoms: vec![], bonds: vec![BondHandle::Id(BondId(8))], dative_bonds: vec![DativeBondHandle::Id(DativeBondId(5))], aromatic_systems: vec![],
            multicenter_bonds: vec![], noncovalent_bonds: vec![], stereo_atoms: vec![], stereo_bonds: vec![],
        },
    )]
    fn test_constraint_edit_from(#[case] input: Constraint, #[case] expected: ConstraintEdit) {
        assert_eq!(ConstraintEdit::from(input), expected);
    }
}
