//! Edit vocabulary for transactional molecule mutation.
//!
//! The `Edit` enum is the caller-facing data-form mutation vocabulary; realized
//! rollback data belongs to the `Undo` journal.
//!
//! Handles (`AtomHandle`, `BondHandle`, ...) are symbolic. `Id(n)` names entity
//! `n` in the transaction's initial host; `New(n)` names the `n`th same-kind
//! entity created in the same [`Edits`] sequence.

use std::slice::Iter;
use std::vec::IntoIter;

use super::aromatic::{AromaticSystemAst, AromaticSystemUpdate};
use super::atom::{AtomAst, AtomUpdate, ElementAst, IsotopeMassAst};
use super::bond::{BondAst, BondUpdate};
use super::constraint::{
    AromaticSystemConstraintAst, AtomConstraintAst, BondConstraintAst, Constraint,
    DativeBondConstraintAst, MulticenterBondConstraintAst, NoncovalentBondConstraintAst,
    StereoAtomConstraintAst, StereoBondConstraintAst,
};
use super::dative::{DativeBondAst, DativeBondUpdate};
use super::electrons::ElectronCountsAst;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::{StereoLigand, StereoLigandKind};
use super::multicenter::{MulticenterBondAst, MulticenterBondUpdate};
use super::noncovalent::{NoncovalentBondAst, NoncovalentBondKindAst, NoncovalentBondUpdate};
use super::remap::{IdCompaction, UndoCompaction};
use super::spin::UnpairedElectronsAst;
use super::stereo::{
    StereoAtomAst, StereoAtomUpdate, StereoBondAst, StereoBondUpdate, StereoConfigurationAst,
};
use super::traits::{Canonicalize, Lattice};
use super::value::ValueAst;

/// One stereo-atom removal in a batched `RemoveStereoAtoms`: id, site, ligand frame, recorded ast.
pub type StereoAtomRemoval = (
    StereoAtomHandle,
    AtomHandle,
    Vec<(AtomHandle, StereoLigandKind)>,
    StereoAtomAst,
);
/// One stereo-bond removal in a batched `RemoveStereoBonds`: id, site (a bond), ligand frame, ast.
pub type StereoBondRemoval = (
    StereoBondHandle,
    BondHandle,
    Vec<(AtomHandle, StereoLigandKind)>,
    StereoBondAst,
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
        old: ElementAst,
        new: ElementAst,
    },
    IsotopeMass {
        old: IsotopeMassAst,
        new: IsotopeMassAst,
    },
    Charge {
        old: ValueAst,
        new: ValueAst,
    },
    ImplicitHydrogens {
        old: ValueAst,
        new: ValueAst,
    },
    LonePairs {
        old: ValueAst,
        new: ValueAst,
    },
    UnpairedElectrons {
        old: UnpairedElectronsAst,
        new: UnpairedElectronsAst,
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
        old: ValueAst,
        new: ValueAst,
    },
    Charge {
        old: ValueAst,
        new: ValueAst,
    },
    UnpairedElectrons {
        old: UnpairedElectronsAst,
        new: UnpairedElectronsAst,
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
    Order { old: ValueAst, new: ValueAst },
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
        old: ElectronCountsAst,
        new: ElectronCountsAst,
    },
    Charge {
        old: ValueAst,
        new: ValueAst,
    },
    UnpairedElectrons {
        old: UnpairedElectronsAst,
        new: UnpairedElectronsAst,
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
        old: ElectronCountsAst,
        new: ElectronCountsAst,
    },
    Charge {
        old: ValueAst,
        new: ValueAst,
    },
    UnpairedElectrons {
        old: UnpairedElectronsAst,
        new: UnpairedElectronsAst,
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
        old: NoncovalentBondKindAst,
        new: NoncovalentBondKindAst,
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
        old: StereoConfigurationAst,
        new: StereoConfigurationAst,
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
        old: StereoConfigurationAst,
        new: StereoConfigurationAst,
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
    pub ast: BondAst,
}

/// One raw mutation entry in an [`Edits`] transaction batch. Topology and
/// removal entries retain their semantic batching.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Edit {
    // Atoms / bonds
    AddAtoms {
        atoms: Vec<AtomAst>,
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
        ast: DativeBondAst,
    },
    RemoveDativeBonds {
        removes: Vec<(DativeBondHandle, Vec<AtomHandle>, DativeBondAst)>,
    },
    ModifyDativeBondField {
        id: DativeBondHandle,
        change: DativeBondFieldChange,
    },

    // Aromatic systems
    AddAromaticSystem {
        atoms: Vec<AtomHandle>,
        ast: AromaticSystemAst,
    },
    RemoveAromaticSystems {
        removes: Vec<(AromaticSystemHandle, Vec<AtomHandle>, AromaticSystemAst)>,
    },
    ModifyAromaticSystemField {
        id: AromaticSystemHandle,
        change: AromaticSystemFieldChange,
    },

    // Multicenter bonds
    AddMulticenterBond {
        atoms: Vec<AtomHandle>,
        ast: MulticenterBondAst,
    },
    RemoveMulticenterBonds {
        removes: Vec<(MulticenterBondHandle, Vec<AtomHandle>, MulticenterBondAst)>,
    },
    ModifyMulticenterBondField {
        id: MulticenterBondHandle,
        change: MulticenterBondFieldChange,
    },

    // Noncovalent bonds
    AddNoncovalentBond {
        atoms: [AtomHandle; 2],
        ast: NoncovalentBondAst,
    },
    RemoveNoncovalentBonds {
        removes: Vec<(NoncovalentBondHandle, [AtomHandle; 2], NoncovalentBondAst)>,
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
        ast: StereoAtomAst,
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
        ast: StereoBondAst,
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
        old: Option<AtomConstraintAst>,
        new: Option<AtomConstraintAst>,
    },
    ModifyBondConstraint {
        id: BondHandle,
        old: Option<BondConstraintAst>,
        new: Option<BondConstraintAst>,
    },
    ModifyDativeBondConstraint {
        id: DativeBondHandle,
        old: Option<DativeBondConstraintAst>,
        new: Option<DativeBondConstraintAst>,
    },
    ModifyAromaticSystemConstraint {
        id: AromaticSystemHandle,
        old: Option<AromaticSystemConstraintAst>,
        new: Option<AromaticSystemConstraintAst>,
    },
    ModifyMulticenterBondConstraint {
        id: MulticenterBondHandle,
        old: Option<MulticenterBondConstraintAst>,
        new: Option<MulticenterBondConstraintAst>,
    },
    ModifyNoncovalentBondConstraint {
        id: NoncovalentBondHandle,
        old: Option<NoncovalentBondConstraintAst>,
        new: Option<NoncovalentBondConstraintAst>,
    },
    ModifyStereoAtomConstraint {
        id: StereoAtomHandle,
        old: Option<StereoAtomConstraintAst>,
        new: Option<StereoAtomConstraintAst>,
    },
    ModifyStereoBondConstraint {
        id: StereoBondHandle,
        old: Option<StereoBondConstraintAst>,
        new: Option<StereoBondConstraintAst>,
    },

    // Molecule-list constraints — a true multiset, so add/remove by value
    // (remove takes the last matching entry; its position is captured for undo).
    AddMoleculeConstraint {
        constraint: Constraint,
    },
    RemoveMoleculeConstraint {
        constraint: Constraint,
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

    pub fn add_atom(&mut self, ast: AtomAst) -> AtomHandle {
        let handle = AtomHandle::New(self.created_atoms);
        self.push(Edit::AddAtoms { atoms: vec![ast] });
        handle
    }

    pub fn add_atoms(&mut self, atoms: impl IntoIterator<Item = AtomAst>) -> Vec<AtomHandle> {
        let atoms: Vec<_> = atoms.into_iter().collect();
        let handles = (self.created_atoms..self.created_atoms + atoms.len())
            .map(AtomHandle::New)
            .collect();
        self.push(Edit::AddAtoms { atoms });
        handles
    }

    pub fn add_bond(&mut self, first: AtomHandle, second: AtomHandle, ast: BondAst) -> BondHandle {
        let handle = BondHandle::New(self.created_bonds);
        self.push(Edit::AddBonds {
            bonds: vec![AddBond {
                endpoints: [first, second],
                ast,
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
        ast: DativeBondAst,
    ) -> DativeBondHandle {
        let handle = DativeBondHandle::New(self.created_dative_bonds);
        self.push(Edit::AddDativeBond { atoms, ast });
        handle
    }

    pub fn add_dative_bonds(
        &mut self,
        bonds: impl IntoIterator<Item = (Vec<AtomHandle>, DativeBondAst)>,
    ) -> Vec<DativeBondHandle> {
        bonds
            .into_iter()
            .map(|(atoms, ast)| self.add_dative_bond(atoms, ast))
            .collect()
    }

    pub fn add_aromatic_system(
        &mut self,
        atoms: Vec<AtomHandle>,
        ast: AromaticSystemAst,
    ) -> AromaticSystemHandle {
        let handle = AromaticSystemHandle::New(self.created_aromatic_systems);
        self.push(Edit::AddAromaticSystem { atoms, ast });
        handle
    }

    pub fn add_aromatic_systems(
        &mut self,
        systems: impl IntoIterator<Item = (Vec<AtomHandle>, AromaticSystemAst)>,
    ) -> Vec<AromaticSystemHandle> {
        systems
            .into_iter()
            .map(|(atoms, ast)| self.add_aromatic_system(atoms, ast))
            .collect()
    }

    pub fn add_multicenter_bond(
        &mut self,
        atoms: Vec<AtomHandle>,
        ast: MulticenterBondAst,
    ) -> MulticenterBondHandle {
        let handle = MulticenterBondHandle::New(self.created_multicenter_bonds);
        self.push(Edit::AddMulticenterBond { atoms, ast });
        handle
    }

    pub fn add_multicenter_bonds(
        &mut self,
        bonds: impl IntoIterator<Item = (Vec<AtomHandle>, MulticenterBondAst)>,
    ) -> Vec<MulticenterBondHandle> {
        bonds
            .into_iter()
            .map(|(atoms, ast)| self.add_multicenter_bond(atoms, ast))
            .collect()
    }

    pub fn add_noncovalent_bond(
        &mut self,
        atoms: [AtomHandle; 2],
        ast: NoncovalentBondAst,
    ) -> NoncovalentBondHandle {
        let handle = NoncovalentBondHandle::New(self.created_noncovalent_bonds);
        self.push(Edit::AddNoncovalentBond { atoms, ast });
        handle
    }

    pub fn add_noncovalent_bonds(
        &mut self,
        bonds: impl IntoIterator<Item = ([AtomHandle; 2], NoncovalentBondAst)>,
    ) -> Vec<NoncovalentBondHandle> {
        bonds
            .into_iter()
            .map(|(atoms, ast)| self.add_noncovalent_bond(atoms, ast))
            .collect()
    }

    pub fn add_stereo_atom(
        &mut self,
        site: AtomHandle,
        ligands: Vec<(AtomHandle, StereoLigandKind)>,
        ast: StereoAtomAst,
    ) -> StereoAtomHandle {
        let handle = StereoAtomHandle::New(self.created_stereo_atoms);
        self.push(Edit::AddStereoAtom { site, ligands, ast });
        handle
    }

    pub fn add_stereo_atoms(
        &mut self,
        atoms: impl IntoIterator<
            Item = (
                AtomHandle,
                Vec<(AtomHandle, StereoLigandKind)>,
                StereoAtomAst,
            ),
        >,
    ) -> Vec<StereoAtomHandle> {
        atoms
            .into_iter()
            .map(|(site, ligands, ast)| self.add_stereo_atom(site, ligands, ast))
            .collect()
    }

    pub fn add_stereo_bond(
        &mut self,
        site: BondHandle,
        ligands: Vec<(AtomHandle, StereoLigandKind)>,
        ast: StereoBondAst,
    ) -> StereoBondHandle {
        let handle = StereoBondHandle::New(self.created_stereo_bonds);
        self.push(Edit::AddStereoBond { site, ligands, ast });
        handle
    }

    pub fn add_stereo_bonds(
        &mut self,
        bonds: impl IntoIterator<
            Item = (
                BondHandle,
                Vec<(AtomHandle, StereoLigandKind)>,
                StereoBondAst,
            ),
        >,
    ) -> Vec<StereoBondHandle> {
        bonds
            .into_iter()
            .map(|(site, ligands, ast)| self.add_stereo_bond(site, ligands, ast))
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
        removes: Vec<(DativeBondHandle, Vec<AtomHandle>, DativeBondAst)>,
    ) {
        self.push(Edit::RemoveDativeBonds { removes });
    }

    pub fn remove_aromatic_systems(
        &mut self,
        removes: Vec<(AromaticSystemHandle, Vec<AtomHandle>, AromaticSystemAst)>,
    ) {
        self.push(Edit::RemoveAromaticSystems { removes });
    }

    pub fn remove_multicenter_bonds(
        &mut self,
        removes: Vec<(MulticenterBondHandle, Vec<AtomHandle>, MulticenterBondAst)>,
    ) {
        self.push(Edit::RemoveMulticenterBonds { removes });
    }

    pub fn remove_noncovalent_bonds(
        &mut self,
        removes: Vec<(NoncovalentBondHandle, [AtomHandle; 2], NoncovalentBondAst)>,
    ) {
        self.push(Edit::RemoveNoncovalentBonds { removes });
    }

    pub fn remove_stereo_atoms(&mut self, removes: Vec<StereoAtomRemoval>) {
        self.push(Edit::RemoveStereoAtoms { removes });
    }

    pub fn remove_stereo_bonds(&mut self, removes: Vec<StereoBondRemoval>) {
        self.push(Edit::RemoveStereoBonds { removes });
    }

    pub fn add_molecule_constraint(&mut self, constraint: Constraint) {
        self.push(Edit::AddMoleculeConstraint { constraint });
    }

    pub fn remove_molecule_constraint(&mut self, constraint: Constraint) {
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
    pub fn update_atom(&mut self, id: AtomHandle, current: &AtomAst, update: &AtomUpdate) {
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
    pub fn update_bond(&mut self, id: BondHandle, current: &BondAst, update: &BondUpdate) {
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
        current: &DativeBondAst,
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
        current: &AromaticSystemAst,
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
        current: &MulticenterBondAst,
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
        current: &NoncovalentBondAst,
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
        current: &StereoAtomAst,
        update: &StereoAtomUpdate,
    ) {
        let updated = current.update(update);
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
        current: &StereoBondAst,
        update: &StereoBondUpdate,
    ) {
        let updated = current.update(update);
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedAtom {
    pub id: AtomId,
    pub ast: AtomAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedBond {
    pub id: BondId,
    pub endpoints: [AtomId; 2],
    pub ast: BondAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedAtom {
    pub id: AtomId,
    pub ast: AtomAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedBond {
    pub id: BondId,
    pub endpoints: [AtomId; 2],
    pub ast: BondAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedDativeBond {
    pub id: DativeBondId,
    pub atoms: Vec<AtomId>,
    pub ast: DativeBondAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedDativeBond {
    pub id: DativeBondId,
    pub atoms: Vec<AtomId>,
    pub ast: DativeBondAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedAromaticSystem {
    pub id: AromaticSystemId,
    pub atoms: Vec<AtomId>,
    pub ast: AromaticSystemAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedAromaticSystem {
    pub id: AromaticSystemId,
    pub atoms: Vec<AtomId>,
    pub ast: AromaticSystemAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedMulticenterBond {
    pub id: MulticenterBondId,
    pub atoms: Vec<AtomId>,
    pub ast: MulticenterBondAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedMulticenterBond {
    pub id: MulticenterBondId,
    pub atoms: Vec<AtomId>,
    pub ast: MulticenterBondAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedNoncovalentBond {
    pub id: NoncovalentBondId,
    pub atoms: [AtomId; 2],
    pub ast: NoncovalentBondAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedNoncovalentBond {
    pub id: NoncovalentBondId,
    pub atoms: [AtomId; 2],
    pub ast: NoncovalentBondAst,
}

// Stereo elements carry both factors: the `site` (atom/bond) and the ordered
// `ligands`, unlike the single-atom-set overlays above.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedStereoAtom {
    pub id: StereoAtomId,
    pub site: AtomId,
    pub ligands: Vec<StereoLigand>,
    pub ast: StereoAtomAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedStereoAtom {
    pub id: StereoAtomId,
    pub site: AtomId,
    pub ligands: Vec<StereoLigand>,
    pub ast: StereoAtomAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddedStereoBond {
    pub id: StereoBondId,
    pub site: BondId,
    pub ligands: Vec<StereoLigand>,
    pub ast: StereoBondAst,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemovedStereoBond {
    pub id: StereoBondId,
    pub site: BondId,
    pub ligands: Vec<StereoLigand>,
    pub ast: StereoBondAst,
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

/// Realized rollback operation produced by the checked transaction path.
// TODO: Review
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
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::boolean::BooleanAst;
    use super::super::constraint::{
        AromaticSystemConstraintsAst, AtomConstraintsAst, BondConstraintsAst,
        DativeBondConstraintsAst, MoleculeConstraint, MulticenterBondConstraintsAst,
        NoncovalentBondConstraintsAst, RingScope, StereoAtomConstraintsAst,
        StereoBondConstraintsAst, StereogenicityAst,
    };
    use super::super::molecule::{MoleculeAst, MoleculeParts};
    use super::super::noncovalent::NoncovalentBondKind;
    use super::super::spin::UnpairedElectronsUpdate;
    use super::super::stereo::{
        StereoConfigurationAst, StereoConfigurationUpdate, StereoCoset, StereoKind, Stereogenicity,
    };
    use super::*;

    #[rstest]
    #[case::element(
        AtomFieldChange::Element {
            old: ElementAst::Lit(Element::C),
            new: ElementAst::Lit(Element::N),
        },
        AtomFieldChange::Element {
            old: ElementAst::Lit(Element::N),
            new: ElementAst::Lit(Element::C),
        },
    )]
    #[case::charge(
        AtomFieldChange::Charge {
            old: ValueAst::Lit(0),
            new: ValueAst::Lit(1),
        },
        AtomFieldChange::Charge {
            old: ValueAst::Lit(1),
            new: ValueAst::Lit(0),
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
            old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        },
        StereoAtomFieldChange::Configuration {
            old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
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
            old: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
            new: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
        },
        StereoBondFieldChange::Configuration {
            old: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
            new: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCoset::Lit(0)),
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
        let atom = AtomAst::from_element(Element::C);
        let bond = BondAst::from_order(1);
        let dative = DativeBondAst::default();
        let aromatic = AromaticSystemAst::default();
        let multicenter = MulticenterBondAst::default();
        let noncovalent = NoncovalentBondAst::default();
        let stereo_atom = StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0));
        let stereo_bond = StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0));
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
                    ast: dative,
                },
                Edit::AddBonds {
                    bonds: vec![AddBond {
                        endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                        ast: bond,
                    }],
                },
                Edit::AddAromaticSystem {
                    atoms: vec![AtomHandle::New(0)],
                    ast: aromatic,
                },
                Edit::AddMulticenterBond {
                    atoms: vec![AtomHandle::New(0)],
                    ast: multicenter,
                },
                Edit::AddNoncovalentBond {
                    atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                    ast: noncovalent,
                },
                Edit::AddStereoAtom {
                    site: AtomHandle::New(0),
                    ligands: vec![(AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom)],
                    ast: stereo_atom,
                },
                Edit::AddStereoBond {
                    site: BondHandle::New(0),
                    ligands: vec![(AtomHandle::New(0), StereoLigandKind::Atom)],
                    ast: stereo_bond,
                },
            ],
        );
    }

    #[rstest]
    fn test_edits_add_atoms() {
        let carbon = AtomAst::from_element(Element::C);
        let nitrogen = AtomAst::from_element(Element::N);
        let oxygen = AtomAst::from_element(Element::O);
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
        let single = BondAst::from_order(1);
        let double = BondAst::from_order(2);
        let triple = BondAst::from_order(3);
        let mut edits = Edits::new();
        let bonds = vec![
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: single.clone(),
            },
            AddBond {
                endpoints: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
                ast: double.clone(),
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
                        ast: triple,
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
                (vec![AtomHandle::Id(AtomId(0))], DativeBondAst::default()),
                (vec![AtomHandle::Id(AtomId(1))], DativeBondAst::default()),
            ]),
            vec![DativeBondHandle::New(0), DativeBondHandle::New(1)],
        );
        assert_eq!(
            edits.add_aromatic_systems([
                (
                    vec![AtomHandle::Id(AtomId(0))],
                    AromaticSystemAst::default(),
                ),
                (
                    vec![AtomHandle::Id(AtomId(1))],
                    AromaticSystemAst::default(),
                ),
            ]),
            vec![AromaticSystemHandle::New(0), AromaticSystemHandle::New(1)],
        );
        assert_eq!(
            edits.add_multicenter_bonds([
                (
                    vec![AtomHandle::Id(AtomId(0))],
                    MulticenterBondAst::default(),
                ),
                (
                    vec![AtomHandle::Id(AtomId(1))],
                    MulticenterBondAst::default(),
                ),
            ]),
            vec![MulticenterBondHandle::New(0), MulticenterBondHandle::New(1)],
        );
        assert_eq!(
            edits.add_noncovalent_bonds([
                (
                    [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    NoncovalentBondAst::default(),
                ),
                (
                    [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
                    NoncovalentBondAst::default(),
                ),
            ]),
            vec![NoncovalentBondHandle::New(0), NoncovalentBondHandle::New(1),],
        );
        assert_eq!(
            edits.add_stereo_atoms([
                (
                    AtomHandle::Id(AtomId(0)),
                    Vec::new(),
                    StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                ),
                (
                    AtomHandle::Id(AtomId(1)),
                    Vec::new(),
                    StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                ),
            ]),
            vec![StereoAtomHandle::New(0), StereoAtomHandle::New(1)],
        );
        assert_eq!(
            edits.add_stereo_bonds([
                (
                    BondHandle::Id(BondId(0)),
                    Vec::new(),
                    StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                ),
                (
                    BondHandle::Id(BondId(1)),
                    Vec::new(),
                    StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
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
                    ast: DativeBondAst::default(),
                },
                Edit::AddDativeBond {
                    atoms: vec![AtomHandle::Id(AtomId(1))],
                    ast: DativeBondAst::default(),
                },
                Edit::AddAromaticSystem {
                    atoms: vec![AtomHandle::Id(AtomId(0))],
                    ast: AromaticSystemAst::default(),
                },
                Edit::AddAromaticSystem {
                    atoms: vec![AtomHandle::Id(AtomId(1))],
                    ast: AromaticSystemAst::default(),
                },
                Edit::AddMulticenterBond {
                    atoms: vec![AtomHandle::Id(AtomId(0))],
                    ast: MulticenterBondAst::default(),
                },
                Edit::AddMulticenterBond {
                    atoms: vec![AtomHandle::Id(AtomId(1))],
                    ast: MulticenterBondAst::default(),
                },
                Edit::AddNoncovalentBond {
                    atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    ast: NoncovalentBondAst::default(),
                },
                Edit::AddNoncovalentBond {
                    atoms: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
                    ast: NoncovalentBondAst::default(),
                },
                Edit::AddStereoAtom {
                    site: AtomHandle::Id(AtomId(0)),
                    ligands: Vec::new(),
                    ast: StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                },
                Edit::AddStereoAtom {
                    site: AtomHandle::Id(AtomId(1)),
                    ligands: Vec::new(),
                    ast: StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                },
                Edit::AddStereoBond {
                    site: BondHandle::Id(BondId(0)),
                    ligands: Vec::new(),
                    ast: StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                },
                Edit::AddStereoBond {
                    site: BondHandle::Id(BondId(1)),
                    ligands: Vec::new(),
                    ast: StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
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
            DativeBondAst::default(),
        )]);
        edits.remove_aromatic_systems(vec![(
            AromaticSystemHandle::New(0),
            vec![AtomHandle::Id(AtomId(0))],
            AromaticSystemAst::default(),
        )]);
        edits.remove_multicenter_bonds(vec![(
            MulticenterBondHandle::Id(MulticenterBondId(0)),
            vec![AtomHandle::Id(AtomId(0))],
            MulticenterBondAst::default(),
        )]);
        edits.remove_noncovalent_bonds(vec![(
            NoncovalentBondHandle::New(0),
            [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
            NoncovalentBondAst::default(),
        )]);
        edits.remove_stereo_atoms(vec![(
            StereoAtomHandle::Id(StereoAtomId(0)),
            AtomHandle::Id(AtomId(0)),
            vec![(AtomHandle::New(0), StereoLigandKind::Atom)],
            StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )]);
        edits.remove_stereo_bonds(vec![(
            StereoBondHandle::New(0),
            BondHandle::Id(BondId(0)),
            vec![(AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom)],
            StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )]);

        assert_eq!(
            edits.into_iter().collect::<Vec<_>>(),
            vec![
                Edit::RemoveDativeBonds {
                    removes: vec![(
                        DativeBondHandle::Id(DativeBondId(0)),
                        vec![AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                        DativeBondAst::default(),
                    )],
                },
                Edit::RemoveAromaticSystems {
                    removes: vec![(
                        AromaticSystemHandle::New(0),
                        vec![AtomHandle::Id(AtomId(0))],
                        AromaticSystemAst::default(),
                    )],
                },
                Edit::RemoveMulticenterBonds {
                    removes: vec![(
                        MulticenterBondHandle::Id(MulticenterBondId(0)),
                        vec![AtomHandle::Id(AtomId(0))],
                        MulticenterBondAst::default(),
                    )],
                },
                Edit::RemoveNoncovalentBonds {
                    removes: vec![(
                        NoncovalentBondHandle::New(0),
                        [AtomHandle::Id(AtomId(0)), AtomHandle::New(0)],
                        NoncovalentBondAst::default(),
                    )],
                },
                Edit::RemoveStereoAtoms {
                    removes: vec![(
                        StereoAtomHandle::Id(StereoAtomId(0)),
                        AtomHandle::Id(AtomId(0)),
                        vec![(AtomHandle::New(0), StereoLigandKind::Atom)],
                        StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                    )],
                },
                Edit::RemoveStereoBonds {
                    removes: vec![(
                        StereoBondHandle::New(0),
                        BondHandle::Id(BondId(0)),
                        vec![(AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom)],
                        StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
                    )],
                },
            ],
        );
    }

    #[rstest]
    fn test_edits_molecule_constraint() {
        let constraint = Constraint::Molecule(MoleculeConstraint::Connected { atoms: None });
        let mut edits = Edits::new();
        edits.add_molecule_constraint(constraint.clone());
        edits.remove_molecule_constraint(constraint.clone());

        assert_eq!(
            edits.as_slice(),
            [
                Edit::AddMoleculeConstraint {
                    constraint: constraint.clone(),
                },
                Edit::RemoveMoleculeConstraint { constraint },
            ],
        );
    }

    #[rstest]
    fn test_edits_push() {
        let entry = Edit::AddAtoms {
            atoms: vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::N),
            ],
        };
        let mut edits = Edits::new();
        edits.push(entry.clone());

        assert_eq!(
            edits.add_atom(AtomAst::from_element(Element::O)),
            AtomHandle::New(2)
        );
        assert_eq!(edits.as_slice()[0], entry);
    }

    #[rstest]
    fn test_edits_from_iter() {
        let entries = vec![
            Edit::AddAtoms {
                atoms: vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                ],
            },
            Edit::AddBonds {
                bonds: vec![
                    AddBond {
                        endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                        ast: BondAst::from_order(1),
                    },
                    AddBond {
                        endpoints: [AtomHandle::Id(AtomId(1)), AtomHandle::Id(AtomId(2))],
                        ast: BondAst::from_order(1),
                    },
                ],
            },
            Edit::AddDativeBond {
                atoms: Vec::new(),
                ast: DativeBondAst::default(),
            },
            Edit::AddAromaticSystem {
                atoms: Vec::new(),
                ast: AromaticSystemAst::default(),
            },
            Edit::AddMulticenterBond {
                atoms: Vec::new(),
                ast: MulticenterBondAst::default(),
            },
            Edit::AddNoncovalentBond {
                atoms: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                ast: NoncovalentBondAst::default(),
            },
            Edit::AddStereoAtom {
                site: AtomHandle::Id(AtomId(0)),
                ligands: Vec::new(),
                ast: StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            },
            Edit::AddStereoBond {
                site: BondHandle::Id(BondId(0)),
                ligands: Vec::new(),
                ast: StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            },
        ];
        let mut edits: Edits = entries.clone().into_iter().collect();

        assert_eq!(edits.as_slice(), entries);
        assert_eq!(edits.add_atom(AtomAst::default()), AtomHandle::New(2));
        assert_eq!(
            edits.add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
                BondAst::default(),
            ),
            BondHandle::New(2),
        );
        assert_eq!(
            edits.add_dative_bond(Vec::new(), DativeBondAst::default()),
            DativeBondHandle::New(1),
        );
        assert_eq!(
            edits.add_aromatic_system(Vec::new(), AromaticSystemAst::default()),
            AromaticSystemHandle::New(1),
        );
        assert_eq!(
            edits.add_multicenter_bond(Vec::new(), MulticenterBondAst::default()),
            MulticenterBondHandle::New(1),
        );
        assert_eq!(
            edits.add_noncovalent_bond(
                [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                NoncovalentBondAst::default(),
            ),
            NoncovalentBondHandle::New(1),
        );
        assert_eq!(
            edits.add_stereo_atom(
                AtomHandle::Id(AtomId(0)),
                Vec::new(),
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            StereoAtomHandle::New(1),
        );
        assert_eq!(
            edits.add_stereo_bond(
                BondHandle::Id(BondId(0)),
                Vec::new(),
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
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
                constraint: Constraint::Molecule(MoleculeConstraint::Connected { atoms: None }),
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
        let current = AtomAst::from_element(Element::C)
            .with_isotope_mass(12_u32)
            .with_charge(0_i64)
            .with_implicit_hydrogens(4_i64)
            .with_lone_pairs(0_i64)
            .with_unpaired_electrons((2_u8, 3_u8))
            .with_constraint(AtomConstraintAst::valence(4_i64));
        let update = AtomUpdate {
            element: Some(ElementAst::Lit(Element::N)),
            isotope_mass: Some(IsotopeMassAst::Lit(13)),
            charge: Some(ValueAst::Lit(1)),
            implicit_hydrogens: Some(ValueAst::Lit(3)),
            lone_pairs: Some(ValueAst::Lit(1)),
            unpaired_electrons: UnpairedElectronsUpdate {
                count: None,
                multiplicity: Some(ValueAst::Lit(1)),
            },
            constraints: AtomConstraintsAst::from_iter([
                AtomConstraintAst::valence(ValueAst::Undetermined),
                AtomConstraintAst::degree(2_i64),
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
                        old: ElementAst::Lit(Element::C),
                        new: ElementAst::Lit(Element::N),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::IsotopeMass {
                        old: IsotopeMassAst::Lit(12),
                        new: IsotopeMassAst::Lit(13),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::Charge {
                        old: ValueAst::Lit(0),
                        new: ValueAst::Lit(1),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::ImplicitHydrogens {
                        old: ValueAst::Lit(4),
                        new: ValueAst::Lit(3),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::LonePairs {
                        old: ValueAst::Lit(0),
                        new: ValueAst::Lit(1),
                    },
                },
                Edit::ModifyAtomField {
                    id: AtomHandle::Id(AtomId(7)),
                    change: AtomFieldChange::UnpairedElectrons {
                        old: UnpairedElectronsAst::from((2_u8, 3_u8)),
                        new: UnpairedElectronsAst::from((2_u8, 1_u8)),
                    },
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(7)),
                    old: Some(AtomConstraintAst::valence(4_i64)),
                    new: None,
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(7)),
                    old: None,
                    new: Some(AtomConstraintAst::degree(2_i64)),
                },
            ]
        );

        let expected = current.update(&update);
        let molecule = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![current.clone()],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_atom(AtomHandle::Id(AtomId(0)), &current, &update);
        editor
            .transact(applied_edits)
            .expect("atom update edits should apply");

        assert_eq!(editor.atom(AtomId(0)).ast, &expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AtomAst::from_element(Element::C), AtomUpdate::default())]
    #[case::canonical_field(AtomAst::from_element(Element::C).with_charge(1_i64), AtomUpdate { charge: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(AtomAst::from_element(Element::C), AtomUpdate { constraints: AtomConstraintsAst::from(AtomConstraintAst::valence(ValueAst::Undetermined)), ..Default::default() })]
    fn test_edits_update_atom_identity(#[case] current: AtomAst, #[case] update: AtomUpdate) {
        let mut edits = Edits::new();
        edits.update_atom(AtomHandle::Id(AtomId(0)), &current, &update);

        assert_eq!(edits, Edits::new());
    }

    #[rstest]
    fn test_edits_update_bond() {
        let current = BondAst::from_order(1)
            .with_charge(0_i64)
            .with_unpaired_electrons((2_u8, 3_u8))
            .with_constraint(BondConstraintAst::ring_membership(
                RingScope::Size(6),
                1_i64,
            ));
        let update = BondUpdate {
            order: Some(ValueAst::Lit(2)),
            charge: Some(ValueAst::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate {
                count: None,
                multiplicity: Some(ValueAst::Lit(1)),
            },
            constraints: BondConstraintsAst::from_iter([
                BondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined),
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
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
                        old: ValueAst::Lit(1),
                        new: ValueAst::Lit(2),
                    },
                },
                Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(7)),
                    change: BondFieldChange::Charge {
                        old: ValueAst::Lit(0),
                        new: ValueAst::Undetermined,
                    },
                },
                Edit::ModifyBondField {
                    id: BondHandle::Id(BondId(7)),
                    change: BondFieldChange::UnpairedElectrons {
                        old: UnpairedElectronsAst::from((2_u8, 3_u8)),
                        new: UnpairedElectronsAst::from((2_u8, 1_u8)),
                    },
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(7)),
                    old: None,
                    new: Some(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(7)),
                    old: Some(BondConstraintAst::ring_membership(
                        RingScope::Size(6),
                        1_i64,
                    )),
                    new: None,
                },
            ]
        );

        let expected = current.update(&update);
        let molecule = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::default(), AtomAst::default()],
            bonds: vec![(AtomId(0), AtomId(1), current.clone())],
            ..Default::default()
        });
        let mut editor = molecule.edit();
        let mut applied_edits = Edits::new();
        applied_edits.update_bond(BondHandle::Id(BondId(0)), &current, &update);
        editor
            .transact(applied_edits)
            .expect("bond update edits should apply");

        assert_eq!(editor.bond(BondId(0)).ast, &expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(BondAst::from_order(1), BondUpdate::default())]
    #[case::canonical_field(BondAst::from_order(1).with_charge(1_i64), BondUpdate { charge: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(BondAst::from_order(1), BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)), ..Default::default() })]
    fn test_edits_update_bond_identity(#[case] current: BondAst, #[case] update: BondUpdate) {
        let mut edits = Edits::new();
        edits.update_bond(BondHandle::Id(BondId(0)), &current, &update);

        assert_eq!(edits, Edits::new());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)),
        DativeBondUpdate {
            order: Some(ValueAst::Lit(2)),
            constraints: DativeBondConstraintsAst::from_iter([
                DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined),
                DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            ]),
        },
        vec![
            Edit::ModifyDativeBondField {
                id: DativeBondHandle::Id(DativeBondId(7)),
                change: DativeBondFieldChange::Order { old: ValueAst::Lit(1), new: ValueAst::Lit(2) },
            },
            Edit::ModifyDativeBondConstraint {
                id: DativeBondHandle::Id(DativeBondId(7)),
                old: None,
                new: Some(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
            },
            Edit::ModifyDativeBondConstraint {
                id: DativeBondHandle::Id(DativeBondId(7)),
                old: Some(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)),
                new: None,
            },
        ],
    )]
    fn test_edits_update_dative_bond(
        #[case] current: DativeBondAst,
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

        let expected_ast = current.update(&update);
        let molecule = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::default(), AtomAst::default()],
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

        assert_eq!(editor.dative_bond(DativeBondId(0)).ast, &expected_ast);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(DativeBondAst::from_order(1), DativeBondUpdate::default())]
    #[case::canonical_field(DativeBondAst::from_order(1), DativeBondUpdate { order: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(DativeBondAst::from_order(1), DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)), ..Default::default() })]
    fn test_edits_update_dative_bond_identity(
        #[case] current: DativeBondAst,
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
        AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(AromaticSystemConstraintAst::electron_count(6_i64)),
        AromaticSystemUpdate {
            electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])),
            charge: Some(ValueAst::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Lit(1)) },
            constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)),
        },
        vec![
            Edit::ModifyAromaticSystemField {
                id: AromaticSystemHandle::Id(AromaticSystemId(7)),
                change: AromaticSystemFieldChange::Electrons { old: ElectronCountsAst::Lit(vec![1, 1, 1]), new: ElectronCountsAst::Lit(vec![2, 2, 2]) },
            },
            Edit::ModifyAromaticSystemField {
                id: AromaticSystemHandle::Id(AromaticSystemId(7)),
                change: AromaticSystemFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Undetermined },
            },
            Edit::ModifyAromaticSystemField {
                id: AromaticSystemHandle::Id(AromaticSystemId(7)),
                change: AromaticSystemFieldChange::UnpairedElectrons { old: UnpairedElectronsAst::from((2_u8, 3_u8)), new: UnpairedElectronsAst::from((2_u8, 1_u8)) },
            },
            Edit::ModifyAromaticSystemConstraint {
                id: AromaticSystemHandle::Id(AromaticSystemId(7)),
                old: Some(AromaticSystemConstraintAst::electron_count(6_i64)),
                new: None,
            },
        ],
    )]
    fn test_edits_update_aromatic_system(
        #[case] current: AromaticSystemAst,
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

        let expected_ast = current.update(&update);
        let molecule = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::default(),
                AtomAst::default(),
                AtomAst::default(),
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
            editor.aromatic_system(AromaticSystemId(0)).ast,
            &expected_ast,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AromaticSystemAst::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate::default())]
    #[case::canonical_field(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(1_i64), AromaticSystemUpdate { charge: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(AromaticSystemAst::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() })]
    fn test_edits_update_aromatic_system_identity(
        #[case] current: AromaticSystemAst,
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
        MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(MulticenterBondConstraintAst::electron_count(6_i64)),
        MulticenterBondUpdate {
            electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])),
            charge: Some(ValueAst::Undetermined),
            unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Lit(1)) },
            constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)),
        },
        vec![
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(7)),
                change: MulticenterBondFieldChange::Electrons { old: ElectronCountsAst::Lit(vec![1, 1, 1]), new: ElectronCountsAst::Lit(vec![2, 2, 2]) },
            },
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(7)),
                change: MulticenterBondFieldChange::Charge { old: ValueAst::Lit(0), new: ValueAst::Undetermined },
            },
            Edit::ModifyMulticenterBondField {
                id: MulticenterBondHandle::Id(MulticenterBondId(7)),
                change: MulticenterBondFieldChange::UnpairedElectrons { old: UnpairedElectronsAst::from((2_u8, 3_u8)), new: UnpairedElectronsAst::from((2_u8, 1_u8)) },
            },
            Edit::ModifyMulticenterBondConstraint {
                id: MulticenterBondHandle::Id(MulticenterBondId(7)),
                old: Some(MulticenterBondConstraintAst::electron_count(6_i64)),
                new: None,
            },
        ],
    )]
    fn test_edits_update_multicenter_bond(
        #[case] current: MulticenterBondAst,
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

        let expected_ast = current.update(&update);
        let molecule = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![
                AtomAst::default(),
                AtomAst::default(),
                AtomAst::default(),
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
            editor.multicenter_bond(MulticenterBondId(0)).ast,
            &expected_ast,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(MulticenterBondAst::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate::default())]
    #[case::canonical_field(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(1_i64), MulticenterBondUpdate { charge: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(MulticenterBondAst::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() })]
    fn test_edits_update_multicenter_bond_identity(
        #[case] current: MulticenterBondAst,
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
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(true)),
        NoncovalentBondUpdate {
            kind: Some(NoncovalentBondKindAst::Undetermined),
            constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(BooleanAst::Undetermined)),
        },
        vec![
            Edit::ModifyNoncovalentBondField {
                id: NoncovalentBondHandle::Id(NoncovalentBondId(7)),
                change: NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindAst::Undetermined },
            },
            Edit::ModifyNoncovalentBondConstraint {
                id: NoncovalentBondHandle::Id(NoncovalentBondId(7)),
                old: Some(NoncovalentBondConstraintAst::intramolecular(true)),
                new: None,
            },
        ],
    )]
    fn test_edits_update_noncovalent_bond(
        #[case] current: NoncovalentBondAst,
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

        let expected_ast = current.update(&update);
        let molecule = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::default(), AtomAst::default()],
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
            editor.noncovalent_bond(NoncovalentBondId(0)).ast,
            &expected_ast,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate::default())]
    #[case::same_kind(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)), ..Default::default() })]
    #[case::absent_constraint_removal(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(BooleanAst::Undetermined)), ..Default::default() })]
    fn test_edits_update_noncovalent_bond_identity(
        #[case] current: NoncovalentBondAst,
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
        StereoAtomAst { configuration: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 0_u32), constraints: StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))) },
        StereoAtomUpdate {
            configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) },
            constraints: StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)),
        },
        vec![
            Edit::ModifyStereoAtomField {
                id: StereoAtomHandle::Id(StereoAtomId(7)),
                change: StereoAtomFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 1_u32) },
            },
            Edit::ModifyStereoAtomConstraint {
                id: StereoAtomHandle::Id(StereoAtomId(7)),
                old: Some(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
                new: None,
            },
        ],
    )]
    fn test_edits_update_stereo_atom(
        #[case] current: StereoAtomAst,
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

        let expected_ast = current.update(&update);
        let molecule = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::default(); 5],
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

        assert_eq!(editor.stereo_atom(StereoAtomId(0)).ast, &expected_ast);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoAtomAst::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate::default())]
    #[case::relative(StereoAtomAst::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() })]
    #[case::absent_constraint_removal(StereoAtomAst::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { constraints: StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)), ..Default::default() })]
    fn test_edits_update_stereo_atom_identity(
        #[case] current: StereoAtomAst,
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
        StereoBondAst { configuration: StereoConfigurationAst::kinded(StereoKind::CisTrans, 0_u32), constraints: StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))) },
        StereoBondUpdate {
            configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) },
            constraints: StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)),
        },
        vec![
            Edit::ModifyStereoBondField {
                id: StereoBondHandle::Id(StereoBondId(7)),
                change: StereoBondFieldChange::Configuration { old: StereoConfigurationAst::kinded(StereoKind::CisTrans, 0_u32), new: StereoConfigurationAst::kinded(StereoKind::CisTrans, 1_u32) },
            },
            Edit::ModifyStereoBondConstraint {
                id: StereoBondHandle::Id(StereoBondId(7)),
                old: Some(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))),
                new: None,
            },
        ],
    )]
    fn test_edits_update_stereo_bond(
        #[case] current: StereoBondAst,
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

        let expected_ast = current.update(&update);
        let molecule = MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::default(); 6],
            bonds: vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
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

        assert_eq!(editor.stereo_bond(StereoBondId(0)).ast, &expected_ast);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoBondAst::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate::default())]
    #[case::relative(StereoBondAst::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() })]
    #[case::absent_constraint_removal(StereoBondAst::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { constraints: StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)), ..Default::default() })]
    fn test_edits_update_stereo_bond_identity(
        #[case] current: StereoBondAst,
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
            old: ValueAst::Lit(1),
            new: ValueAst::Lit(2),
        };
        assert_eq!(
            change.clone().inverse(),
            BondFieldChange::Order {
                old: ValueAst::Lit(2),
                new: ValueAst::Lit(1),
            },
        );
        assert_eq!(change.clone().inverse().inverse(), change);
    }
}
