//! Edit vocabulary for transactional molecule mutation.
//!
//! The `Edit` enum is the data-form vocabulary for `MoleculeEditor::transact`
//! `Edit` is caller-facing mutation data; realized rollback data belongs to
//! the `Undo` journal.
//!
//! Refs (`AtomHandle`, `BondHandle`, ...) are symbolic and appear only inside
//! `Edit`. `Id(_)` references an existing entity; `New(N)` references the
//! entity created by the Nth Edit earlier in the same batch.

use super::aromatic::{AromaticSystemAst, AromaticSystemUpdate};
use super::atom::{AtomAst, AtomUpdate, ElementAst, IsotopeMassAst};
use super::bond::{BondAst, BondUpdate};
use super::constraint::{
    AromaticSystemConstraintAst, AtomConstraintAst, BondConstraintAst, Constraint, Constraints,
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
use super::spin::SpinStateAst;
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

/// Handle to an atom within an edit batch: either an existing `AtomId` or the
/// Nth atom-creating Edit earlier in the same transaction batch.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AtomHandle {
    Id(AtomId),
    New(usize),
}

/// Handle to a bond within an edit batch (an existing `BondId` or the Nth
/// bond-creating Edit earlier in the batch).
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
    Spin {
        old: SpinStateAst,
        new: SpinStateAst,
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
            Self::Spin { old, new } => Self::Spin { old: new, new: old },
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
    Spin {
        old: SpinStateAst,
        new: SpinStateAst,
    },
}

impl BondFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Order { old, new } => Self::Order { old: new, new: old },
            Self::Charge { old, new } => Self::Charge { old: new, new: old },
            Self::Spin { old, new } => Self::Spin { old: new, new: old },
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
    Spin {
        old: SpinStateAst,
        new: SpinStateAst,
    },
}

impl AromaticSystemFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Electrons { old, new } => Self::Electrons { old: new, new: old },
            Self::Charge { old, new } => Self::Charge { old: new, new: old },
            Self::Spin { old, new } => Self::Spin { old: new, new: old },
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
    Spin {
        old: SpinStateAst,
        new: SpinStateAst,
    },
}

impl MulticenterBondFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Electrons { old, new } => Self::Electrons { old: new, new: old },
            Self::Charge { old, new } => Self::Charge { old: new, new: old },
            Self::Spin { old, new } => Self::Spin { old: new, new: old },
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

/// Single mutation operation. Compose `Vec<Edit>` into a transaction batch
/// Topology edits are bulk primitives; single-item helpers are constructors on
/// `Edit`.
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
    // (remove takes the first matching entry; its position is captured for undo).
    AddMoleculeConstraint {
        constraint: Constraint,
    },
    RemoveMoleculeConstraint {
        constraint: Constraint,
    },
}

impl Edit {
    pub fn add_atom(ast: AtomAst) -> Self {
        Self::AddAtoms { atoms: vec![ast] }
    }

    pub fn add_bond(first: AtomHandle, second: AtomHandle, ast: BondAst) -> Self {
        Self::AddBonds {
            bonds: vec![AddBond {
                endpoints: [first, second],
                ast,
            }],
        }
    }

    pub fn remove_atom(id: AtomHandle) -> Self {
        Self::RemoveTopology {
            atoms: vec![id],
            bonds: Vec::new(),
        }
    }

    pub fn remove_bond(id: BondHandle) -> Self {
        Self::RemoveTopology {
            atoms: Vec::new(),
            bonds: vec![id],
        }
    }

    /// Project an atom update into checked host-relative edits.
    pub fn for_atom_update(id: AtomHandle, current: &AtomAst, update: &AtomUpdate) -> Vec<Self> {
        let mut edits = Vec::new();
        if let Some(new) = &update.element {
            if !current.element.canonical_eq(new) {
                edits.push(Self::ModifyAtomField {
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
                edits.push(Self::ModifyAtomField {
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
                edits.push(Self::ModifyAtomField {
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
                edits.push(Self::ModifyAtomField {
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
                edits.push(Self::ModifyAtomField {
                    id: id.clone(),
                    change: AtomFieldChange::LonePairs {
                        old: current.lone_pairs.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        let new_spin = current.spin.update(&update.spin);
        if !current.spin.canonical_eq(&new_spin) {
            edits.push(Self::ModifyAtomField {
                id: id.clone(),
                change: AtomFieldChange::Spin {
                    old: current.spin.clone(),
                    new: new_spin,
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
                edits.push(Self::ModifyAtomConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
        edits
    }

    /// Project a localized-bond update into checked host-relative edits.
    pub fn for_bond_update(id: BondHandle, current: &BondAst, update: &BondUpdate) -> Vec<Self> {
        let mut edits = Vec::new();
        if let Some(new) = &update.order {
            if !current.order.canonical_eq(new) {
                edits.push(Self::ModifyBondField {
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
                edits.push(Self::ModifyBondField {
                    id: id.clone(),
                    change: BondFieldChange::Charge {
                        old: current.charge.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        let new_spin = current.spin.update(&update.spin);
        if !current.spin.canonical_eq(&new_spin) {
            edits.push(Self::ModifyBondField {
                id: id.clone(),
                change: BondFieldChange::Spin {
                    old: current.spin.clone(),
                    new: new_spin,
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
                edits.push(Self::ModifyBondConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
        edits
    }

    /// Project a dative-bond update into checked host-relative edits.
    pub fn for_dative_bond_update(
        id: DativeBondHandle,
        current: &DativeBondAst,
        update: &DativeBondUpdate,
    ) -> Vec<Self> {
        let mut edits = Vec::new();
        if let Some(new) = &update.order {
            if !current.order.canonical_eq(new) {
                edits.push(Self::ModifyDativeBondField {
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
                edits.push(Self::ModifyDativeBondConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
        edits
    }

    /// Project an aromatic-system update into checked host-relative edits.
    pub fn for_aromatic_system_update(
        id: AromaticSystemHandle,
        current: &AromaticSystemAst,
        update: &AromaticSystemUpdate,
    ) -> Vec<Self> {
        let mut edits = Vec::new();
        if let Some(new) = &update.electrons {
            if !current.electrons.canonical_eq(new) {
                edits.push(Self::ModifyAromaticSystemField {
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
                edits.push(Self::ModifyAromaticSystemField {
                    id: id.clone(),
                    change: AromaticSystemFieldChange::Charge {
                        old: current.charge.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        let new_spin = current.spin.update(&update.spin);
        if !current.spin.canonical_eq(&new_spin) {
            edits.push(Self::ModifyAromaticSystemField {
                id: id.clone(),
                change: AromaticSystemFieldChange::Spin {
                    old: current.spin.clone(),
                    new: new_spin,
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
                edits.push(Self::ModifyAromaticSystemConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
        edits
    }

    /// Project a multicenter-bond update into checked host-relative edits.
    pub fn for_multicenter_bond_update(
        id: MulticenterBondHandle,
        current: &MulticenterBondAst,
        update: &MulticenterBondUpdate,
    ) -> Vec<Self> {
        let mut edits = Vec::new();
        if let Some(new) = &update.electrons {
            if !current.electrons.canonical_eq(new) {
                edits.push(Self::ModifyMulticenterBondField {
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
                edits.push(Self::ModifyMulticenterBondField {
                    id: id.clone(),
                    change: MulticenterBondFieldChange::Charge {
                        old: current.charge.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        let new_spin = current.spin.update(&update.spin);
        if !current.spin.canonical_eq(&new_spin) {
            edits.push(Self::ModifyMulticenterBondField {
                id: id.clone(),
                change: MulticenterBondFieldChange::Spin {
                    old: current.spin.clone(),
                    new: new_spin,
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
                edits.push(Self::ModifyMulticenterBondConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
        edits
    }

    /// Project a noncovalent-bond update into checked host-relative edits.
    pub fn for_noncovalent_bond_update(
        id: NoncovalentBondHandle,
        current: &NoncovalentBondAst,
        update: &NoncovalentBondUpdate,
    ) -> Vec<Self> {
        let mut edits = Vec::new();
        if let Some(new) = &update.kind {
            if !current.kind.canonical_eq(new) {
                edits.push(Self::ModifyNoncovalentBondField {
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
                edits.push(Self::ModifyNoncovalentBondConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
        edits
    }

    /// Project a stereo-atom update into checked host-relative edits.
    pub fn for_stereo_atom_update(
        id: StereoAtomHandle,
        current: &StereoAtomAst,
        update: &StereoAtomUpdate,
    ) -> Vec<Self> {
        let mut edits = Vec::new();
        let updated = current.update(update);
        if !current.configuration.canonical_eq(&updated.configuration) {
            edits.push(Self::ModifyStereoAtomField {
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
                edits.push(Self::ModifyStereoAtomConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
        edits
    }

    /// Project a stereo-bond update into checked host-relative edits.
    pub fn for_stereo_bond_update(
        id: StereoBondHandle,
        current: &StereoBondAst,
        update: &StereoBondUpdate,
    ) -> Vec<Self> {
        let mut edits = Vec::new();
        let updated = current.update(update);
        if !current.configuration.canonical_eq(&updated.configuration) {
            edits.push(Self::ModifyStereoBondField {
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
                edits.push(Self::ModifyStereoBondConstraint {
                    id: id.clone(),
                    old,
                    new,
                });
            }
        }
        edits
    }
}

// Handles for overlay relations (an existing id or the Nth created earlier in the batch).
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

    pub fn rollback_into(self, constraints: &mut Constraints) {
        let mut items = constraints.take();
        for modified in self.modified {
            if let Some(pos) = items.iter().position(|c| *c == modified.new) {
                items[pos] = modified.old;
            }
        }
        for removed in self.removed {
            let position = removed.position.min(items.len());
            items.insert(position, removed.constraint);
        }
        *constraints = items.into_iter().collect();
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
        DativeBondConstraintsAst, MulticenterBondConstraintsAst, NoncovalentBondConstraintsAst,
        RingScope, StereoAtomConstraintsAst, StereoBondConstraintsAst, StereogenicityAst,
    };
    use super::super::noncovalent::NoncovalentBondKind;
    use super::super::spin::SpinStateUpdate;
    use super::super::stereo::{
        StereoConfigurationAst, StereoConfigurationUpdate, StereoCosetAst, StereoKind,
        Stereogenicity,
    };
    use super::*;

    #[fixture]
    fn carbon_atom() -> AtomAst {
        AtomAst::from_element(Element::C)
    }

    #[fixture]
    fn single_bond() -> BondAst {
        BondAst {
            order: ValueAst::Lit(1),
            ..BondAst::default()
        }
    }

    #[rstest]
    #[case::id(AtomHandle::Id(AtomId(3)))]
    #[case::new(AtomHandle::New(2))]
    fn test_atom_ref_variants(#[case] r: AtomHandle) {
        assert_eq!(r.clone(), r);
    }

    #[rstest]
    #[case::id(BondHandle::Id(BondId(5)))]
    #[case::new(BondHandle::New(0))]
    fn test_bond_ref_variants(#[case] r: BondHandle) {
        assert_eq!(r.clone(), r);
    }

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
            old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)),
            new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
        },
        StereoAtomFieldChange::Configuration {
            old: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
            new: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)),
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
            old: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCosetAst::Lit(0)),
            new: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
        },
        StereoBondFieldChange::Configuration {
            old: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
            new: StereoConfigurationAst::kinded(StereoKind::CisTrans, StereoCosetAst::Lit(0)),
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
    fn test_edit_add_atom(carbon_atom: AtomAst) {
        assert_eq!(
            Edit::add_atom(carbon_atom.clone()),
            Edit::AddAtoms {
                atoms: vec![carbon_atom],
            },
        );
    }

    #[rstest]
    fn test_edit_add_bond(single_bond: BondAst) {
        assert_eq!(
            Edit::add_bond(
                AtomHandle::Id(AtomId(0)),
                AtomHandle::Id(AtomId(1)),
                single_bond.clone()
            ),
            Edit::AddBonds {
                bonds: vec![AddBond {
                    endpoints: [AtomHandle::Id(AtomId(0)), AtomHandle::Id(AtomId(1))],
                    ast: single_bond,
                }],
            },
        );
    }

    #[rstest]
    fn test_edit_remove_atom() {
        assert_eq!(
            Edit::remove_atom(AtomHandle::Id(AtomId(2))),
            Edit::RemoveTopology {
                atoms: vec![AtomHandle::Id(AtomId(2))],
                bonds: Vec::new(),
            },
        );
    }

    #[rstest]
    fn test_edit_remove_bond() {
        assert_eq!(
            Edit::remove_bond(BondHandle::Id(BondId(4))),
            Edit::RemoveTopology {
                atoms: Vec::new(),
                bonds: vec![BondHandle::Id(BondId(4))],
            },
        );
    }

    #[rstest]
    fn test_edit_for_atom_update() {
        let current = AtomAst::from_element(Element::C)
            .with_isotope_mass(12_u32)
            .with_charge(0_i64)
            .with_implicit_hydrogens(4_i64)
            .with_lone_pairs(0_i64)
            .with_spin((2_u8, 3_u8))
            .with_constraint(AtomConstraintAst::valence(4_i64));
        let update = AtomUpdate {
            element: Some(ElementAst::Lit(Element::N)),
            isotope_mass: Some(IsotopeMassAst::Lit(13)),
            charge: Some(ValueAst::Lit(1)),
            implicit_hydrogens: Some(ValueAst::Lit(3)),
            lone_pairs: Some(ValueAst::Lit(1)),
            spin: SpinStateUpdate {
                unpaired: None,
                multiplicity: Some(ValueAst::Lit(1)),
            },
            constraints: AtomConstraintsAst::from_iter([
                AtomConstraintAst::valence(ValueAst::Undetermined),
                AtomConstraintAst::degree(2_i64),
            ]),
        };
        assert_eq!(
            Edit::for_atom_update(AtomHandle::Id(AtomId(7)), &current, &update),
            vec![
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
                    change: AtomFieldChange::Spin {
                        old: SpinStateAst::from((2_u8, 3_u8)),
                        new: SpinStateAst::from((2_u8, 1_u8)),
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
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AtomAst::from_element(Element::C), AtomUpdate::default())]
    #[case::canonical_field(AtomAst::from_element(Element::C).with_charge(1_i64), AtomUpdate { charge: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(AtomAst::from_element(Element::C), AtomUpdate { constraints: AtomConstraintsAst::from(AtomConstraintAst::valence(ValueAst::Undetermined)), ..Default::default() })]
    fn test_edit_for_atom_update_identity(#[case] current: AtomAst, #[case] update: AtomUpdate) {
        assert_eq!(
            Edit::for_atom_update(AtomHandle::Id(AtomId(0)), &current, &update),
            Vec::new()
        );
    }

    #[rstest]
    fn test_edit_for_bond_update() {
        let current = BondAst::from_order(1)
            .with_charge(0_i64)
            .with_spin((2_u8, 3_u8))
            .with_constraint(BondConstraintAst::ring_membership(
                RingScope::Size(6),
                1_i64,
            ));
        let update = BondUpdate {
            order: Some(ValueAst::Lit(2)),
            charge: Some(ValueAst::Undetermined),
            spin: SpinStateUpdate {
                unpaired: None,
                multiplicity: Some(ValueAst::Lit(1)),
            },
            constraints: BondConstraintsAst::from_iter([
                BondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined),
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            ]),
        };
        assert_eq!(
            Edit::for_bond_update(BondHandle::Id(BondId(7)), &current, &update),
            vec![
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
                    change: BondFieldChange::Spin {
                        old: SpinStateAst::from((2_u8, 3_u8)),
                        new: SpinStateAst::from((2_u8, 1_u8)),
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
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(BondAst::from_order(1), BondUpdate::default())]
    #[case::canonical_field(BondAst::from_order(1).with_charge(1_i64), BondUpdate { charge: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(BondAst::from_order(1), BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)), ..Default::default() })]
    fn test_edit_for_bond_update_identity(#[case] current: BondAst, #[case] update: BondUpdate) {
        assert_eq!(
            Edit::for_bond_update(BondHandle::Id(BondId(0)), &current, &update),
            Vec::new()
        );
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
    fn test_edit_for_dative_bond_update(
        #[case] current: DativeBondAst,
        #[case] update: DativeBondUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        assert_eq!(
            Edit::for_dative_bond_update(DativeBondHandle::Id(DativeBondId(7)), &current, &update),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(DativeBondAst::from_order(1), DativeBondUpdate::default())]
    #[case::canonical_field(DativeBondAst::from_order(1), DativeBondUpdate { order: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(DativeBondAst::from_order(1), DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)), ..Default::default() })]
    fn test_edit_for_dative_bond_update_identity(
        #[case] current: DativeBondAst,
        #[case] update: DativeBondUpdate,
    ) {
        assert_eq!(
            Edit::for_dative_bond_update(DativeBondHandle::Id(DativeBondId(0)), &current, &update),
            Vec::new(),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraint(
        AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_spin((2_u8, 3_u8)).with_constraint(AromaticSystemConstraintAst::electron_count(6_i64)),
        AromaticSystemUpdate {
            electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])),
            charge: Some(ValueAst::Undetermined),
            spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) },
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
                change: AromaticSystemFieldChange::Spin { old: SpinStateAst::from((2_u8, 3_u8)), new: SpinStateAst::from((2_u8, 1_u8)) },
            },
            Edit::ModifyAromaticSystemConstraint {
                id: AromaticSystemHandle::Id(AromaticSystemId(7)),
                old: Some(AromaticSystemConstraintAst::electron_count(6_i64)),
                new: None,
            },
        ],
    )]
    fn test_edit_for_aromatic_system_update(
        #[case] current: AromaticSystemAst,
        #[case] update: AromaticSystemUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        assert_eq!(
            Edit::for_aromatic_system_update(
                AromaticSystemHandle::Id(AromaticSystemId(7)),
                &current,
                &update,
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AromaticSystemAst::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate::default())]
    #[case::canonical_field(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(1_i64), AromaticSystemUpdate { charge: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(AromaticSystemAst::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() })]
    fn test_edit_for_aromatic_system_update_identity(
        #[case] current: AromaticSystemAst,
        #[case] update: AromaticSystemUpdate,
    ) {
        assert_eq!(
            Edit::for_aromatic_system_update(
                AromaticSystemHandle::Id(AromaticSystemId(0)),
                &current,
                &update,
            ),
            Vec::new(),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraint(
        MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_spin((2_u8, 3_u8)).with_constraint(MulticenterBondConstraintAst::electron_count(6_i64)),
        MulticenterBondUpdate {
            electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])),
            charge: Some(ValueAst::Undetermined),
            spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) },
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
                change: MulticenterBondFieldChange::Spin { old: SpinStateAst::from((2_u8, 3_u8)), new: SpinStateAst::from((2_u8, 1_u8)) },
            },
            Edit::ModifyMulticenterBondConstraint {
                id: MulticenterBondHandle::Id(MulticenterBondId(7)),
                old: Some(MulticenterBondConstraintAst::electron_count(6_i64)),
                new: None,
            },
        ],
    )]
    fn test_edit_for_multicenter_bond_update(
        #[case] current: MulticenterBondAst,
        #[case] update: MulticenterBondUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        assert_eq!(
            Edit::for_multicenter_bond_update(
                MulticenterBondHandle::Id(MulticenterBondId(7)),
                &current,
                &update,
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(MulticenterBondAst::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate::default())]
    #[case::canonical_field(MulticenterBondAst::from_electrons(vec![1, 1, 1]).with_charge(1_i64), MulticenterBondUpdate { charge: Some(ValueAst::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(MulticenterBondAst::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() })]
    fn test_edit_for_multicenter_bond_update_identity(
        #[case] current: MulticenterBondAst,
        #[case] update: MulticenterBondUpdate,
    ) {
        assert_eq!(
            Edit::for_multicenter_bond_update(
                MulticenterBondHandle::Id(MulticenterBondId(0)),
                &current,
                &update,
            ),
            Vec::new(),
        );
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
    fn test_edit_for_noncovalent_bond_update(
        #[case] current: NoncovalentBondAst,
        #[case] update: NoncovalentBondUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        assert_eq!(
            Edit::for_noncovalent_bond_update(
                NoncovalentBondHandle::Id(NoncovalentBondId(7)),
                &current,
                &update,
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate::default())]
    #[case::same_kind(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)), ..Default::default() })]
    #[case::absent_constraint_removal(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(BooleanAst::Undetermined)), ..Default::default() })]
    fn test_edit_for_noncovalent_bond_update_identity(
        #[case] current: NoncovalentBondAst,
        #[case] update: NoncovalentBondUpdate,
    ) {
        assert_eq!(
            Edit::for_noncovalent_bond_update(
                NoncovalentBondHandle::Id(NoncovalentBondId(0)),
                &current,
                &update,
            ),
            Vec::new(),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::configuration_and_constraint(
        StereoAtomAst { configuration: StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 0_u32), constraints: StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))) },
        StereoAtomUpdate {
            configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCosetAst::Lit(1)) },
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
    fn test_edit_for_stereo_atom_update(
        #[case] current: StereoAtomAst,
        #[case] update: StereoAtomUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        assert_eq!(
            Edit::for_stereo_atom_update(
                StereoAtomHandle::Id(StereoAtomId(7)),
                &current,
                &update,
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoAtomAst::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate::default())]
    #[case::relative(StereoAtomAst::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() })]
    #[case::absent_constraint_removal(StereoAtomAst::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { constraints: StereoAtomConstraintsAst::from(StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)), ..Default::default() })]
    fn test_edit_for_stereo_atom_update_identity(
        #[case] current: StereoAtomAst,
        #[case] update: StereoAtomUpdate,
    ) {
        assert_eq!(
            Edit::for_stereo_atom_update(
                StereoAtomHandle::Id(StereoAtomId(0)),
                &current,
                &update,
            ),
            Vec::new(),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::configuration_and_constraint(
        StereoBondAst { configuration: StereoConfigurationAst::kinded(StereoKind::CisTrans, 0_u32), constraints: StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic))) },
        StereoBondUpdate {
            configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCosetAst::Lit(1)) },
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
    fn test_edit_for_stereo_bond_update(
        #[case] current: StereoBondAst,
        #[case] update: StereoBondUpdate,
        #[case] expected: Vec<Edit>,
    ) {
        assert_eq!(
            Edit::for_stereo_bond_update(
                StereoBondHandle::Id(StereoBondId(7)),
                &current,
                &update,
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoBondAst::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate::default())]
    #[case::relative(StereoBondAst::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() })]
    #[case::absent_constraint_removal(StereoBondAst::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { constraints: StereoBondConstraintsAst::from(StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined)), ..Default::default() })]
    fn test_edit_for_stereo_bond_update_identity(
        #[case] current: StereoBondAst,
        #[case] update: StereoBondUpdate,
    ) {
        assert_eq!(
            Edit::for_stereo_bond_update(
                StereoBondHandle::Id(StereoBondId(0)),
                &current,
                &update,
            ),
            Vec::new(),
        );
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
