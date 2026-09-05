//! Resolved edit vocabulary: the `Delta` counterpart of the deferred `Edit`.
//!
//! A `Delta` is one resolved edit over a `Molecule`, referencing entities by stable
//! ids in the molecule's own id space (no positional `New`). The vocabulary is closed
//! under inversion — every delta's inverse is another delta.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::Hash;
use std::iter;
use std::mem::{discriminant, Discriminant};
use std::slice::{Iter, IterMut};
use std::vec::IntoIter;

use umol_perm::{DynPermutation, Permutation};

use super::aromatic::{AromaticSystemForm, AromaticSystemUpdate};
use super::atom::{AtomForm, AtomUpdate};
use super::bond::{BondForm, BondUpdate};
use super::constraint::{
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AtomConstraintForm,
    AtomConstraintKey, BondConstraintForm, BondConstraintKey, Constraint,
    ConstraintFrameActionDomain, ConstraintFrameActions, DativeBondConstraintForm,
    DativeBondConstraintKey, MulticenterBondConstraintForm, MulticenterBondConstraintKey,
    NoncovalentBondConstraintForm, NoncovalentBondConstraintKey, StereoAtomConstraintForm,
    StereoAtomConstraintKey, StereoBondConstraintForm, StereoBondConstraintKey,
};
use super::dative::{DativeBondForm, DativeBondUpdate};
use super::edit::{
    AromaticSystemFieldChange, AtomFieldChange, BondFieldChange, DativeBondFieldChange,
    MulticenterBondFieldChange, NoncovalentBondFieldChange, StereoAtomFieldChange,
    StereoBondFieldChange,
};
use super::entity::Entity;
use super::error::Contradiction;
use super::frame::OverlaysFrameAction;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::StereoLigand;
use super::multicenter::{MulticenterBondForm, MulticenterBondUpdate};
use super::noncovalent::{NoncovalentBondForm, NoncovalentBondUpdate};
#[cfg(test)]
use super::remap::IdRemapping;
use super::stereo::{
    StereoAtomForm, StereoAtomUpdate, StereoBondForm, StereoBondUpdate, StereoConfigurationForm,
    StereoKind,
};
use super::traits::{EntityPatch, FrameTransport, Lattice, Normalize};

/// A resolved edit to a single atom.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AtomDelta {
    Add {
        id: AtomId,
        attributes: AtomForm,
    },
    Remove {
        id: AtomId,
        attributes: AtomForm,
    },
    ModifyField {
        id: AtomId,
        change: AtomFieldChange,
    },
    ModifyConstraint {
        id: AtomId,
        old: Option<AtomConstraintForm>,
        new: Option<AtomConstraintForm>,
    },
}

impl AtomDelta {
    /// The inverse delta: `Add`↔`Remove`; `ModifyField` / `ModifyConstraint` swap old/new.
    pub fn inverse(self) -> Self {
        match self {
            Self::Add { id, attributes } => Self::Remove { id, attributes },
            Self::Remove { id, attributes } => Self::Add { id, attributes },
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

    /// Project an atom update into resolved deltas.
    pub fn for_update(id: AtomId, current: &AtomForm, update: &AtomUpdate) -> Vec<Self> {
        let mut deltas = Vec::new();
        if let Some(new) = &update.element {
            if !current.element.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
                    change: AtomFieldChange::Element {
                        old: current.element.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.isotope_mass {
            if !current.isotope_mass.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
                    change: AtomFieldChange::IsotopeMass {
                        old: current.isotope_mass.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.charge {
            if !current.charge.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
                    change: AtomFieldChange::Charge {
                        old: current.charge.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.implicit_hydrogens {
            if !current.implicit_hydrogens.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
                    change: AtomFieldChange::ImplicitHydrogens {
                        old: current.implicit_hydrogens.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.lone_pairs {
            if !current.lone_pairs.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
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
            .normalized_eq(&new_unpaired_electrons)
        {
            deltas.push(Self::ModifyField {
                id,
                change: AtomFieldChange::UnpairedElectrons {
                    old: current.unpaired_electrons.clone(),
                    new: new_unpaired_electrons,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            if !options_normalized_eq(&old, &new) {
                deltas.push(Self::ModifyConstraint { id, old, new });
            }
        }
        deltas
    }
}

/// A resolved edit to a single bond.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BondDelta {
    Add {
        id: BondId,
        atoms: [AtomId; 2],
        attributes: BondForm,
    },
    Remove {
        id: BondId,
        atoms: [AtomId; 2],
        attributes: BondForm,
    },
    ModifyField {
        id: BondId,
        change: BondFieldChange,
    },
    ModifyConstraint {
        id: BondId,
        old: Option<BondConstraintForm>,
        new: Option<BondConstraintForm>,
    },
}

impl BondDelta {
    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id,
                atoms,
                attributes,
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id,
                atoms,
                attributes,
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

    /// Project a localized-bond update into resolved deltas.
    pub fn for_update(id: BondId, current: &BondForm, update: &BondUpdate) -> Vec<Self> {
        let mut deltas = Vec::new();
        if let Some(new) = &update.order {
            if !current.order.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
                    change: BondFieldChange::Order {
                        old: current.order.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.charge {
            if !current.charge.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
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
            .normalized_eq(&new_unpaired_electrons)
        {
            deltas.push(Self::ModifyField {
                id,
                change: BondFieldChange::UnpairedElectrons {
                    old: current.unpaired_electrons.clone(),
                    new: new_unpaired_electrons,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            if !options_normalized_eq(&old, &new) {
                deltas.push(Self::ModifyConstraint { id, old, new });
            }
        }
        deltas
    }
}

/// A resolved edit to a single dative bond. `donors`/`acceptor` are the directed
/// participants (structural payload, like `BondDelta::atoms`); identity is the id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DativeBondDelta {
    Add {
        id: DativeBondId,
        donors: Vec<AtomId>,
        acceptor: AtomId,
        attributes: DativeBondForm,
    },
    Remove {
        id: DativeBondId,
        donors: Vec<AtomId>,
        acceptor: AtomId,
        attributes: DativeBondForm,
    },
    ModifyField {
        id: DativeBondId,
        change: DativeBondFieldChange,
    },
    ModifyConstraint {
        id: DativeBondId,
        old: Option<DativeBondConstraintForm>,
        new: Option<DativeBondConstraintForm>,
    },
}

impl DativeBondDelta {
    pub(crate) fn uses_participant_frame(&self) -> bool {
        match self {
            Self::Add { .. } | Self::Remove { .. } => true,
            Self::ModifyField { change, .. } => match change {
                DativeBondFieldChange::Order { .. } => false,
            },
            Self::ModifyConstraint { old, new, .. } => old
                .iter()
                .chain(new)
                .any(DativeBondConstraintForm::uses_participant_frame),
        }
    }

    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                donors,
                acceptor,
                attributes,
            } => Self::Remove {
                id,
                donors,
                acceptor,
                attributes,
            },
            Self::Remove {
                id,
                donors,
                acceptor,
                attributes,
            } => Self::Add {
                id,
                donors,
                acceptor,
                attributes,
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

    /// Project a dative-bond update into resolved deltas.
    pub fn for_update(
        id: DativeBondId,
        current: &DativeBondForm,
        update: &DativeBondUpdate,
    ) -> Vec<Self> {
        let mut deltas = Vec::new();
        if let Some(new) = &update.order {
            if !current.order.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
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
            if !options_normalized_eq(&old, &new) {
                deltas.push(Self::ModifyConstraint { id, old, new });
            }
        }
        deltas
    }
}

/// A resolved edit to a single aromatic system. `atoms` are the member atoms
/// (structural payload); identity is the id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticSystemDelta {
    Add {
        id: AromaticSystemId,
        atoms: Vec<AtomId>,
        attributes: AromaticSystemForm,
    },
    Remove {
        id: AromaticSystemId,
        atoms: Vec<AtomId>,
        attributes: AromaticSystemForm,
    },
    ModifyField {
        id: AromaticSystemId,
        change: AromaticSystemFieldChange,
    },
    ModifyConstraint {
        id: AromaticSystemId,
        old: Option<AromaticSystemConstraintForm>,
        new: Option<AromaticSystemConstraintForm>,
    },
}

impl AromaticSystemDelta {
    pub(crate) fn uses_participant_frame(&self) -> bool {
        match self {
            Self::Add { .. } | Self::Remove { .. } => true,
            Self::ModifyField { change, .. } => match change {
                AromaticSystemFieldChange::Electrons { .. } => true,
                AromaticSystemFieldChange::Charge { .. }
                | AromaticSystemFieldChange::UnpairedElectrons { .. } => false,
            },
            Self::ModifyConstraint { old, new, .. } => old
                .iter()
                .chain(new)
                .any(AromaticSystemConstraintForm::uses_participant_frame),
        }
    }

    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id,
                atoms,
                attributes,
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id,
                atoms,
                attributes,
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

    /// Project an aromatic-system update into resolved deltas.
    pub fn for_update(
        id: AromaticSystemId,
        current: &AromaticSystemForm,
        update: &AromaticSystemUpdate,
    ) -> Vec<Self> {
        let mut deltas = Vec::new();
        if let Some(new) = &update.electrons {
            if !current.electrons.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
                    change: AromaticSystemFieldChange::Electrons {
                        old: current.electrons.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.charge {
            if !current.charge.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
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
            .normalized_eq(&new_unpaired_electrons)
        {
            deltas.push(Self::ModifyField {
                id,
                change: AromaticSystemFieldChange::UnpairedElectrons {
                    old: current.unpaired_electrons.clone(),
                    new: new_unpaired_electrons,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            if !options_normalized_eq(&old, &new) {
                deltas.push(Self::ModifyConstraint { id, old, new });
            }
        }
        deltas
    }
}

/// A resolved edit to a single multicenter bond. `atoms` are the member atoms
/// (structural payload); identity is the id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MulticenterBondDelta {
    Add {
        id: MulticenterBondId,
        atoms: Vec<AtomId>,
        attributes: MulticenterBondForm,
    },
    Remove {
        id: MulticenterBondId,
        atoms: Vec<AtomId>,
        attributes: MulticenterBondForm,
    },
    ModifyField {
        id: MulticenterBondId,
        change: MulticenterBondFieldChange,
    },
    ModifyConstraint {
        id: MulticenterBondId,
        old: Option<MulticenterBondConstraintForm>,
        new: Option<MulticenterBondConstraintForm>,
    },
}

impl MulticenterBondDelta {
    pub(crate) fn uses_participant_frame(&self) -> bool {
        match self {
            Self::Add { .. } | Self::Remove { .. } => true,
            Self::ModifyField { change, .. } => match change {
                MulticenterBondFieldChange::Electrons { .. } => true,
                MulticenterBondFieldChange::Charge { .. }
                | MulticenterBondFieldChange::UnpairedElectrons { .. } => false,
            },
            Self::ModifyConstraint { old, new, .. } => old
                .iter()
                .chain(new)
                .any(MulticenterBondConstraintForm::uses_participant_frame),
        }
    }

    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id,
                atoms,
                attributes,
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id,
                atoms,
                attributes,
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

    /// Project a multicenter-bond update into resolved deltas.
    pub fn for_update(
        id: MulticenterBondId,
        current: &MulticenterBondForm,
        update: &MulticenterBondUpdate,
    ) -> Vec<Self> {
        let mut deltas = Vec::new();
        if let Some(new) = &update.electrons {
            if !current.electrons.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
                    change: MulticenterBondFieldChange::Electrons {
                        old: current.electrons.clone(),
                        new: new.clone(),
                    },
                });
            }
        }
        if let Some(new) = &update.charge {
            if !current.charge.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
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
            .normalized_eq(&new_unpaired_electrons)
        {
            deltas.push(Self::ModifyField {
                id,
                change: MulticenterBondFieldChange::UnpairedElectrons {
                    old: current.unpaired_electrons.clone(),
                    new: new_unpaired_electrons,
                },
            });
        }
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            if !options_normalized_eq(&old, &new) {
                deltas.push(Self::ModifyConstraint { id, old, new });
            }
        }
        deltas
    }
}

/// A resolved edit to a single noncovalent bond. `atoms` are its two participants
/// (structural payload); identity is the id.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NoncovalentBondDelta {
    Add {
        id: NoncovalentBondId,
        atoms: [AtomId; 2],
        attributes: NoncovalentBondForm,
    },
    Remove {
        id: NoncovalentBondId,
        atoms: [AtomId; 2],
        attributes: NoncovalentBondForm,
    },
    ModifyField {
        id: NoncovalentBondId,
        change: NoncovalentBondFieldChange,
    },
    ModifyConstraint {
        id: NoncovalentBondId,
        old: Option<NoncovalentBondConstraintForm>,
        new: Option<NoncovalentBondConstraintForm>,
    },
}

impl NoncovalentBondDelta {
    pub(crate) fn uses_participant_frame(&self) -> bool {
        match self {
            Self::Add { .. } | Self::Remove { .. } => true,
            Self::ModifyField { change, .. } => match change {
                NoncovalentBondFieldChange::Kind { .. } => false,
            },
            Self::ModifyConstraint { old, new, .. } => old
                .iter()
                .chain(new)
                .any(NoncovalentBondConstraintForm::uses_participant_frame),
        }
    }

    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id,
                atoms,
                attributes,
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id,
                atoms,
                attributes,
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

    /// Project a noncovalent-bond update into resolved deltas.
    pub fn for_update(
        id: NoncovalentBondId,
        current: &NoncovalentBondForm,
        update: &NoncovalentBondUpdate,
    ) -> Vec<Self> {
        let mut deltas = Vec::new();
        if let Some(new) = &update.kind {
            if !current.kind.normalized_eq(new) {
                deltas.push(Self::ModifyField {
                    id,
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
            if !options_normalized_eq(&old, &new) {
                deltas.push(Self::ModifyConstraint { id, old, new });
            }
        }
        deltas
    }
}

/// A resolved change to a single stereo atom.
///
/// `site` and `ligands` carry the structural incidence while `id` carries identity. Field and
/// constraint modifications state both the expected old value and the replacement value, matching
/// the absolute delta vocabulary used by the other entity kinds.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoAtomDelta {
    Add {
        id: StereoAtomId,
        site: AtomId,
        ligands: Vec<StereoLigand>,
        attributes: StereoAtomForm,
    },
    Remove {
        id: StereoAtomId,
        site: AtomId,
        ligands: Vec<StereoLigand>,
        attributes: StereoAtomForm,
    },
    ModifyField {
        id: StereoAtomId,
        change: StereoAtomFieldChange,
    },
    ModifyConstraint {
        id: StereoAtomId,
        /// Serialization context: the geometry kind the constraint renders/parses against (its
        /// permutation degree, `~` shortcut). `None` for a kind-free constraint on an
        /// `Undetermined`-geometry center. Not read by apply/normalize/diff.
        kind: Option<StereoKind>,
        old: Option<StereoAtomConstraintForm>,
        new: Option<StereoAtomConstraintForm>,
    },
}

impl StereoAtomDelta {
    pub fn id(&self) -> StereoAtomId {
        match self {
            Self::Add { id, .. }
            | Self::Remove { id, .. }
            | Self::ModifyField { id, .. }
            | Self::ModifyConstraint { id, .. } => *id,
        }
    }

    pub(crate) fn uses_participant_frame(&self) -> bool {
        match self {
            Self::Add { .. } | Self::Remove { .. } => true,
            Self::ModifyField { change, .. } => match change {
                StereoAtomFieldChange::Configuration { .. } => true,
            },
            Self::ModifyConstraint { old, new, .. } => old
                .iter()
                .chain(new)
                .any(StereoAtomConstraintForm::uses_participant_frame),
        }
    }

    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                site,
                ligands,
                attributes,
            } => Self::Remove {
                id,
                site,
                ligands,
                attributes,
            },
            Self::Remove {
                id,
                site,
                ligands,
                attributes,
            } => Self::Add {
                id,
                site,
                ligands,
                attributes,
            },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, kind, old, new } => Self::ModifyConstraint {
                id,
                kind,
                old: new,
                new: old,
            },
        }
    }

    /// Project a stereo-atom update into resolved deltas.
    pub fn for_update(
        id: StereoAtomId,
        current: &StereoAtomForm,
        update: &StereoAtomUpdate,
    ) -> Vec<Self> {
        let mut deltas = Vec::new();
        let updated = current.update(update);
        if !current.configuration.normalized_eq(&updated.configuration) {
            deltas.push(Self::ModifyField {
                id,
                change: StereoAtomFieldChange::Configuration {
                    old: current.configuration.clone(),
                    new: updated.configuration.clone(),
                },
            });
        }
        let kind = update
            .configuration
            .kind()
            .or_else(|| current.configuration.kind())
            .or_else(|| updated.configuration.kind());
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            if !options_normalized_eq(&old, &new) {
                deltas.push(Self::ModifyConstraint { id, kind, old, new });
            }
        }
        deltas
    }
}

/// A resolved change to a single stereo bond.
///
/// `site` and `ligands` carry the structural incidence while `id` carries identity. Field and
/// constraint modifications use the same absolute before/after vocabulary as `StereoAtomDelta`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoBondDelta {
    Add {
        id: StereoBondId,
        site: BondId,
        ligands: Vec<StereoLigand>,
        attributes: StereoBondForm,
    },
    Remove {
        id: StereoBondId,
        site: BondId,
        ligands: Vec<StereoLigand>,
        attributes: StereoBondForm,
    },
    ModifyField {
        id: StereoBondId,
        change: StereoBondFieldChange,
    },
    ModifyConstraint {
        id: StereoBondId,
        /// Serialization context — see `StereoAtomDelta::ModifyConstraint`.
        kind: Option<StereoKind>,
        old: Option<StereoBondConstraintForm>,
        new: Option<StereoBondConstraintForm>,
    },
}

impl StereoBondDelta {
    pub fn id(&self) -> StereoBondId {
        match self {
            Self::Add { id, .. }
            | Self::Remove { id, .. }
            | Self::ModifyField { id, .. }
            | Self::ModifyConstraint { id, .. } => *id,
        }
    }

    pub(crate) fn uses_participant_frame(&self) -> bool {
        match self {
            Self::Add { .. } | Self::Remove { .. } => true,
            Self::ModifyField { change, .. } => match change {
                StereoBondFieldChange::Configuration { .. } => true,
            },
            Self::ModifyConstraint { old, new, .. } => old
                .iter()
                .chain(new)
                .any(StereoBondConstraintForm::uses_participant_frame),
        }
    }

    pub fn inverse(self) -> Self {
        match self {
            Self::Add {
                id,
                site,
                ligands,
                attributes,
            } => Self::Remove {
                id,
                site,
                ligands,
                attributes,
            },
            Self::Remove {
                id,
                site,
                ligands,
                attributes,
            } => Self::Add {
                id,
                site,
                ligands,
                attributes,
            },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: change.inverse(),
            },
            Self::ModifyConstraint { id, kind, old, new } => Self::ModifyConstraint {
                id,
                kind,
                old: new,
                new: old,
            },
        }
    }

    /// Project a stereo-bond update into resolved deltas.
    pub fn for_update(
        id: StereoBondId,
        current: &StereoBondForm,
        update: &StereoBondUpdate,
    ) -> Vec<Self> {
        let mut deltas = Vec::new();
        let updated = current.update(update);
        if !current.configuration.normalized_eq(&updated.configuration) {
            deltas.push(Self::ModifyField {
                id,
                change: StereoBondFieldChange::Configuration {
                    old: current.configuration.clone(),
                    new: updated.configuration.clone(),
                },
            });
        }
        let kind = update
            .configuration
            .kind()
            .or_else(|| current.configuration.kind())
            .or_else(|| updated.configuration.kind());
        for constraint in update.constraints.iter() {
            let old = current.constraints.get(constraint.key()).cloned();
            let new = (!constraint.is_undetermined()).then(|| constraint.clone());
            if !options_normalized_eq(&old, &new) {
                deltas.push(Self::ModifyConstraint { id, kind, old, new });
            }
        }
        deltas
    }
}

fn transport_optional<T: FrameTransport>(
    value: Option<T>,
    action: &T::Action,
) -> Option<Option<T>> {
    match value {
        Some(value) => value.reframe_by(action).map(Some),
        None => Some(None),
    }
}

impl FrameTransport for DativeBondDelta {
    type Action = DynPermutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        Some(match self {
            Self::Add {
                id,
                donors,
                acceptor,
                attributes,
            } => Self::Add {
                id,
                donors: action.act(&donors)?,
                acceptor,
                attributes: attributes.reframe_by(action)?,
            },
            Self::Remove {
                id,
                donors,
                acceptor,
                attributes,
            } => Self::Remove {
                id,
                donors: action.act(&donors)?,
                acceptor,
                attributes: attributes.reframe_by(action)?,
            },
            Self::ModifyField { id, change } => Self::ModifyField { id, change },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: transport_optional(old, action)?,
                new: transport_optional(new, action)?,
            },
        })
    }
}

impl FrameTransport for AromaticSystemDelta {
    type Action = DynPermutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        Some(match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id,
                atoms: action.act(&atoms)?,
                attributes: attributes.reframe_by(action)?,
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id,
                atoms: action.act(&atoms)?,
                attributes: attributes.reframe_by(action)?,
            },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: match change {
                    AromaticSystemFieldChange::Electrons { old, new } => {
                        AromaticSystemFieldChange::Electrons {
                            old: old.reframe_by(action)?,
                            new: new.reframe_by(action)?,
                        }
                    }
                    AromaticSystemFieldChange::Charge { old, new } => {
                        AromaticSystemFieldChange::Charge { old, new }
                    }
                    AromaticSystemFieldChange::UnpairedElectrons { old, new } => {
                        AromaticSystemFieldChange::UnpairedElectrons { old, new }
                    }
                },
            },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: transport_optional(old, action)?,
                new: transport_optional(new, action)?,
            },
        })
    }
}

impl FrameTransport for MulticenterBondDelta {
    type Action = DynPermutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        Some(match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id,
                atoms: action.act(&atoms)?,
                attributes: attributes.reframe_by(action)?,
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id,
                atoms: action.act(&atoms)?,
                attributes: attributes.reframe_by(action)?,
            },
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: match change {
                    MulticenterBondFieldChange::Electrons { old, new } => {
                        MulticenterBondFieldChange::Electrons {
                            old: old.reframe_by(action)?,
                            new: new.reframe_by(action)?,
                        }
                    }
                    MulticenterBondFieldChange::Charge { old, new } => {
                        MulticenterBondFieldChange::Charge { old, new }
                    }
                    MulticenterBondFieldChange::UnpairedElectrons { old, new } => {
                        MulticenterBondFieldChange::UnpairedElectrons { old, new }
                    }
                },
            },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: transport_optional(old, action)?,
                new: transport_optional(new, action)?,
            },
        })
    }
}

impl FrameTransport for NoncovalentBondDelta {
    type Action = DynPermutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        if action.degree() != 2 {
            return None;
        }
        Some(match self {
            Self::Add {
                id,
                atoms,
                attributes,
            } => Self::Add {
                id,
                atoms: action.act(&atoms)?.try_into().ok()?,
                attributes: attributes.reframe_by(action)?,
            },
            Self::Remove {
                id,
                atoms,
                attributes,
            } => Self::Remove {
                id,
                atoms: action.act(&atoms)?.try_into().ok()?,
                attributes: attributes.reframe_by(action)?,
            },
            Self::ModifyField { id, change } => Self::ModifyField { id, change },
            Self::ModifyConstraint { id, old, new } => Self::ModifyConstraint {
                id,
                old: transport_optional(old, action)?,
                new: transport_optional(new, action)?,
            },
        })
    }
}

impl FrameTransport for StereoAtomDelta {
    type Action = Permutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        Some(match self {
            Self::Add {
                id,
                site,
                ligands,
                attributes,
            } => {
                if action.degree() != ligands.len() {
                    return None;
                }
                Self::Add {
                    id,
                    site,
                    ligands: action.act(&ligands),
                    attributes: attributes.reframe_by(action)?,
                }
            }
            Self::Remove {
                id,
                site,
                ligands,
                attributes,
            } => {
                if action.degree() != ligands.len() {
                    return None;
                }
                Self::Remove {
                    id,
                    site,
                    ligands: action.act(&ligands),
                    attributes: attributes.reframe_by(action)?,
                }
            }
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: match change {
                    StereoAtomFieldChange::Configuration { old, new } => {
                        StereoAtomFieldChange::Configuration {
                            old: old.reframe_by(action)?,
                            new: new.reframe_by(action)?,
                        }
                    }
                },
            },
            Self::ModifyConstraint { id, kind, old, new } => {
                if kind.is_some_and(|kind| kind.act(0, *action).is_none()) {
                    return None;
                }
                Self::ModifyConstraint {
                    id,
                    kind,
                    old: transport_optional(old, action)?,
                    new: transport_optional(new, action)?,
                }
            }
        })
    }
}

impl FrameTransport for StereoBondDelta {
    type Action = Permutation;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        StereoKind::CisTrans.act(0, *action)?;
        Some(match self {
            Self::Add {
                id,
                site,
                ligands,
                attributes,
            } => {
                if action.degree() != ligands.len() {
                    return None;
                }
                Self::Add {
                    id,
                    site,
                    ligands: action.act(&ligands),
                    attributes: attributes.reframe_by(action)?,
                }
            }
            Self::Remove {
                id,
                site,
                ligands,
                attributes,
            } => {
                if action.degree() != ligands.len() {
                    return None;
                }
                Self::Remove {
                    id,
                    site,
                    ligands: action.act(&ligands),
                    attributes: attributes.reframe_by(action)?,
                }
            }
            Self::ModifyField { id, change } => Self::ModifyField {
                id,
                change: match change {
                    StereoBondFieldChange::Configuration { old, new } => {
                        StereoBondFieldChange::Configuration {
                            old: old.reframe_by(action)?,
                            new: new.reframe_by(action)?,
                        }
                    }
                },
            },
            Self::ModifyConstraint { id, kind, old, new } => {
                if kind.is_some_and(|kind| kind.act(0, *action).is_none()) {
                    return None;
                }
                Self::ModifyConstraint {
                    id,
                    kind,
                    old: transport_optional(old, action)?,
                    new: transport_optional(new, action)?,
                }
            }
        })
    }
}

/// A resolved change to the molecule-level constraint set, as a set-diff.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

    pub(crate) fn collect_frame_action_domain(&self, domain: &mut ConstraintFrameActionDomain) {
        match self {
            Self::Add(constraint) | Self::Remove(constraint) => {
                constraint.collect_frame_action_domain(domain);
            }
        }
    }

    pub(crate) fn reframe_by_actions(
        self,
        actions: &impl ConstraintFrameActions,
    ) -> Result<Self, Entity> {
        Ok(match self {
            Self::Add(constraint) => Self::Add(constraint.reframe_by_actions(actions)?),
            Self::Remove(constraint) => Self::Remove(constraint.reframe_by_actions(actions)?),
        })
    }
}

impl FrameTransport for ConstraintDelta {
    type Action = OverlaysFrameAction;

    fn reframe_by(self, actions: &Self::Action) -> Option<Self> {
        self.reframe_by_actions(actions).ok()
    }
}

/// One rule-relative modification carried by a reaction, over any entity kind.
///
/// A delta is an algebraic value, not an instruction to a particular molecule. Existing entities
/// are identified in the reaction-owned id frame anchored by its left-hand side; additions extend
/// that frame with new ids. The frame remains meaningful before any host or match is selected.
///
/// Deltas are complete: one that removes an entity carries that entity, and one that changes a
/// field carries the old value as well as the new. The vocabulary is therefore closed under
/// inversion, [`Delta::inverse`] is total, and applying it twice returns the original delta.
/// Completeness is distinct from the DPO gluing condition: strict reaction application separately
/// rejects a match when host structure outside the explicit rule would be left dangling. An SqPO
/// application policy could instead cascade such host-relative removals, but the original delta
/// would not describe that additional realized effect.
///
/// Applying a reaction converts its deltas into [`Edit`](crate::ir::edit::Edit) values against the
/// matched host; the match supplies the translation between the two id spaces.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Delta {
    Atom(AtomDelta),
    Bond(BondDelta),
    DativeBond(DativeBondDelta),
    AromaticSystem(AromaticSystemDelta),
    MulticenterBond(MulticenterBondDelta),
    NoncovalentBond(NoncovalentBondDelta),
    StereoAtom(StereoAtomDelta),
    StereoBond(StereoBondDelta),
    Constraint(ConstraintDelta),
}

impl Delta {
    /// The inverse delta. Inversion is total and involutive.
    pub fn inverse(self) -> Self {
        match self {
            Self::Atom(delta) => Self::Atom(delta.inverse()),
            Self::Bond(delta) => Self::Bond(delta.inverse()),
            Self::DativeBond(delta) => Self::DativeBond(delta.inverse()),
            Self::AromaticSystem(delta) => Self::AromaticSystem(delta.inverse()),
            Self::MulticenterBond(delta) => Self::MulticenterBond(delta.inverse()),
            Self::NoncovalentBond(delta) => Self::NoncovalentBond(delta.inverse()),
            Self::StereoAtom(delta) => Self::StereoAtom(delta.inverse()),
            Self::StereoBond(delta) => Self::StereoBond(delta.inverse()),
            Self::Constraint(delta) => Self::Constraint(delta.inverse()),
        }
    }

    pub(crate) fn collect_frame_action_domain(&self, domain: &mut ConstraintFrameActionDomain) {
        match self {
            Self::Atom(_) | Self::Bond(_) => {}
            Self::DativeBond(delta) => {
                if delta.uses_participant_frame() {
                    let id = match delta {
                        DativeBondDelta::Add { id, .. }
                        | DativeBondDelta::Remove { id, .. }
                        | DativeBondDelta::ModifyField { id, .. }
                        | DativeBondDelta::ModifyConstraint { id, .. } => *id,
                    };
                    domain.insert_dative_bond(id);
                }
            }
            Self::AromaticSystem(delta) => {
                if delta.uses_participant_frame() {
                    let id = match delta {
                        AromaticSystemDelta::Add { id, .. }
                        | AromaticSystemDelta::Remove { id, .. }
                        | AromaticSystemDelta::ModifyField { id, .. }
                        | AromaticSystemDelta::ModifyConstraint { id, .. } => *id,
                    };
                    domain.insert_aromatic_system(id);
                }
            }
            Self::MulticenterBond(delta) => {
                if delta.uses_participant_frame() {
                    let id = match delta {
                        MulticenterBondDelta::Add { id, .. }
                        | MulticenterBondDelta::Remove { id, .. }
                        | MulticenterBondDelta::ModifyField { id, .. }
                        | MulticenterBondDelta::ModifyConstraint { id, .. } => *id,
                    };
                    domain.insert_multicenter_bond(id);
                }
            }
            Self::NoncovalentBond(delta) => {
                if delta.uses_participant_frame() {
                    let id = match delta {
                        NoncovalentBondDelta::Add { id, .. }
                        | NoncovalentBondDelta::Remove { id, .. }
                        | NoncovalentBondDelta::ModifyField { id, .. }
                        | NoncovalentBondDelta::ModifyConstraint { id, .. } => *id,
                    };
                    domain.insert_noncovalent_bond(id);
                }
            }
            Self::StereoAtom(delta) => {
                if delta.uses_participant_frame() {
                    domain.insert_stereo_atom(delta.id());
                }
            }
            Self::StereoBond(delta) => {
                if delta.uses_participant_frame() {
                    domain.insert_stereo_bond(delta.id());
                }
            }
            Self::Constraint(delta) => delta.collect_frame_action_domain(domain),
        }
    }

    pub(crate) fn reframe_by_actions(
        self,
        actions: &impl ConstraintFrameActions,
    ) -> Result<Self, Entity> {
        Ok(match self {
            Self::Atom(delta) => Self::Atom(delta),
            Self::Bond(delta) => Self::Bond(delta),
            Self::DativeBond(delta) => {
                if delta.uses_participant_frame() {
                    let id = match &delta {
                        DativeBondDelta::Add { id, .. }
                        | DativeBondDelta::Remove { id, .. }
                        | DativeBondDelta::ModifyField { id, .. }
                        | DativeBondDelta::ModifyConstraint { id, .. } => *id,
                    };
                    let entity = Entity::DativeBond(id);
                    let action = actions.dative_bond_action(id).ok_or(entity)?;
                    Self::DativeBond(delta.reframe_by(action).ok_or(entity)?)
                } else {
                    Self::DativeBond(delta)
                }
            }
            Self::AromaticSystem(delta) => {
                if delta.uses_participant_frame() {
                    let id = match &delta {
                        AromaticSystemDelta::Add { id, .. }
                        | AromaticSystemDelta::Remove { id, .. }
                        | AromaticSystemDelta::ModifyField { id, .. }
                        | AromaticSystemDelta::ModifyConstraint { id, .. } => *id,
                    };
                    let entity = Entity::AromaticSystem(id);
                    let action = actions.aromatic_system_action(id).ok_or(entity)?;
                    Self::AromaticSystem(delta.reframe_by(action).ok_or(entity)?)
                } else {
                    Self::AromaticSystem(delta)
                }
            }
            Self::MulticenterBond(delta) => {
                if delta.uses_participant_frame() {
                    let id = match &delta {
                        MulticenterBondDelta::Add { id, .. }
                        | MulticenterBondDelta::Remove { id, .. }
                        | MulticenterBondDelta::ModifyField { id, .. }
                        | MulticenterBondDelta::ModifyConstraint { id, .. } => *id,
                    };
                    let entity = Entity::MulticenterBond(id);
                    let action = actions.multicenter_bond_action(id).ok_or(entity)?;
                    Self::MulticenterBond(delta.reframe_by(action).ok_or(entity)?)
                } else {
                    Self::MulticenterBond(delta)
                }
            }
            Self::NoncovalentBond(delta) => {
                if delta.uses_participant_frame() {
                    let id = match &delta {
                        NoncovalentBondDelta::Add { id, .. }
                        | NoncovalentBondDelta::Remove { id, .. }
                        | NoncovalentBondDelta::ModifyField { id, .. }
                        | NoncovalentBondDelta::ModifyConstraint { id, .. } => *id,
                    };
                    let entity = Entity::NoncovalentBond(id);
                    let action = actions.noncovalent_bond_action(id).ok_or(entity)?;
                    Self::NoncovalentBond(delta.reframe_by(action).ok_or(entity)?)
                } else {
                    Self::NoncovalentBond(delta)
                }
            }
            Self::StereoAtom(delta) => {
                if delta.uses_participant_frame() {
                    let entity = Entity::StereoAtom(delta.id());
                    let action = actions.stereo_atom_action(delta.id()).ok_or(entity)?;
                    Self::StereoAtom(delta.reframe_by(action).ok_or(entity)?)
                } else {
                    Self::StereoAtom(delta)
                }
            }
            Self::StereoBond(delta) => {
                if delta.uses_participant_frame() {
                    let entity = Entity::StereoBond(delta.id());
                    let action = actions.stereo_bond_action(delta.id()).ok_or(entity)?;
                    Self::StereoBond(delta.reframe_by(action).ok_or(entity)?)
                } else {
                    Self::StereoBond(delta)
                }
            }
            Self::Constraint(delta) => Self::Constraint(delta.reframe_by_actions(actions)?),
        })
    }
}

impl FrameTransport for Delta {
    type Action = OverlaysFrameAction;

    fn reframe_by(self, actions: &Self::Action) -> Option<Self> {
        self.reframe_by_actions(actions).ok()
    }
}

/// Per-variant diff/apply ops for the `EntityPatch` impl, from the `(variant => attributes field)` map:
/// `apply_field`, `diff_field`, `diff_constraints`.
macro_rules! diff_field_ops {
    ($change:ident, $attributes:ident, $constraint:ident, { $($variant:ident => $field:ident),+ $(,)? }) => {
        fn apply_field(attributes: &mut $attributes, change: $change) -> Result<(), Contradiction> {
            match change {
                $(
                    $change::$variant { old, new } => {
                        if !attributes.$field.normalized_eq(&old) {
                            return Err(Contradiction);
                        }
                        attributes.$field = new;
                    }
                )+
            }
            Ok(())
        }

        fn diff_field(lhs: &$attributes, rhs: &$attributes) -> Vec<$change> {
            let mut out = Vec::new();
            $(
                if !lhs.$field.normalized_eq(&rhs.$field) {
                    out.push($change::$variant {
                        old: lhs.$field.clone(),
                        new: rhs.$field.clone(),
                    });
                }
            )+
            out
        }

        #[allow(clippy::type_complexity)]
        fn diff_constraints(
            lhs: &$attributes,
            rhs: &$attributes,
        ) -> Vec<(Option<$constraint>, Option<$constraint>)> {
            let mut lhs_by_key: HashMap<_, $constraint> = HashMap::new();
            for constraint in lhs.constraints.iter() {
                lhs_by_key.insert(constraint.key(), constraint.clone());
            }
            let mut rhs_by_key: HashMap<_, $constraint> = HashMap::new();
            for constraint in rhs.constraints.iter() {
                rhs_by_key.insert(constraint.key(), constraint.clone());
            }
            let mut keys: HashSet<_> = lhs_by_key.keys().cloned().collect();
            keys.extend(rhs_by_key.keys().cloned());
            let mut out = Vec::new();
            for key in keys {
                let l = lhs_by_key.get(&key).cloned();
                let r = rhs_by_key.get(&key).cloned();
                if !options_normalized_eq(&l, &r) {
                    out.push((l, r));
                }
            }
            out
        }
    };
}

/// Normalized equivalence over optional attributes: both absent is equal, both present compares by
/// `normalized_eq`, presence mismatch is unequal.
fn options_normalized_eq<T: Normalize>(l: &Option<T>, r: &Option<T>) -> bool {
    match (l, r) {
        (None, None) => true,
        (Some(a), Some(b)) => a.normalized_eq(b),
        _ => false,
    }
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
                    ) if prev_new.normalized_eq(&next_old) => Some($change::$variant { old, new }),
                )+
                #[allow(unreachable_patterns)]
                _ => None,
            }
        }

        fn field_is_identity(change: &$change) -> bool {
            match change {
                $( $change::$variant { old, new } => old.normalized_eq(new), )+
            }
        }
    };
}

/// One entity's span across a reaction — its slice of the superimposed `L`∪`K`∪`R`. A *state*, not
/// an operation (unlike `Edit` / `Delta`). `lhs()` / `rhs()` read the side values.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EntitySpan<T> {
    /// In the interface `K` — present and identical on both sides.
    Unchanged(T),
    /// In the interface `K` — present on both sides but relabeled (a dynamic entity).
    Modified { lhs: T, rhs: T },
    /// In `R` only — created.
    Added(T),
    /// In `L` only — deleted.
    Removed(T),
}

impl<T> EntitySpan<T> {
    /// The lhs (`L`) value, or `None` if the entity is created.
    pub fn lhs(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Removed(value) | Self::Modified { lhs: value, .. } => {
                Some(value)
            }
            Self::Added(_) => None,
        }
    }

    /// The rhs (`R`) value, or `None` if the entity is deleted.
    pub fn rhs(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Added(value) | Self::Modified { rhs: value, .. } => {
                Some(value)
            }
            Self::Removed(_) => None,
        }
    }

    /// Carry every present side through `f`, declining if any side declines.
    ///
    /// A `Modified` span applies the same `f` to both sides, which is what makes one selected frame
    /// action reach both: the span holds two values against a single participant list.
    pub fn try_map<U>(self, mut f: impl FnMut(T) -> Option<U>) -> Option<EntitySpan<U>> {
        Some(match self {
            Self::Unchanged(value) => EntitySpan::Unchanged(f(value)?),
            Self::Modified { lhs, rhs } => EntitySpan::Modified {
                lhs: f(lhs)?,
                rhs: f(rhs)?,
            },
            Self::Added(value) => EntitySpan::Added(f(value)?),
            Self::Removed(value) => EntitySpan::Removed(f(value)?),
        })
    }
}

impl<T: FrameTransport> FrameTransport for EntitySpan<T> {
    type Action = T::Action;

    fn reframe_by(self, action: &Self::Action) -> Option<Self> {
        self.try_map(|value| value.reframe_by(action))
    }
}

impl<T: Normalize> Normalize for EntitySpan<T> {
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Unchanged(value) => Self::Unchanged(value.normalize()?),
            Self::Modified { lhs, rhs } => {
                let lhs = lhs.normalize()?;
                let rhs = rhs.normalize()?;
                if lhs == rhs {
                    Self::Unchanged(lhs)
                } else {
                    Self::Modified { lhs, rhs }
                }
            }
            Self::Added(value) => Self::Added(value.normalize()?),
            Self::Removed(value) => Self::Removed(value.normalize()?),
        })
    }
}

impl<T: Normalize> EntitySpan<T> {
    /// Superimpose an entity's optional lhs and rhs values into a span — the per-entity kernel of
    /// `ReactionSpan::superimpose`: present-both maps to `Unchanged(lhs)` when the values are
    /// semantically equivalent and to `Modified` otherwise, lhs-only to `Removed`, rhs-only to
    /// `Added`, neither to `None`.
    pub fn superimpose(lhs: Option<T>, rhs: Option<T>) -> Option<Self> {
        match (lhs, rhs) {
            (Some(lhs), Some(rhs)) if lhs.normalized_eq(&rhs) => Some(Self::Unchanged(lhs)),
            (Some(lhs), Some(rhs)) => Some(Self::Modified { lhs, rhs }),
            (Some(lhs), None) => Some(Self::Removed(lhs)),
            (None, Some(rhs)) => Some(Self::Added(rhs)),
            (None, None) => None,
        }
    }
}

/// A molecule-level constraint's span across a reaction — its slice of the superimposed `L`∪`K`∪`R`.
/// A *state*, not an operation (unlike `ConstraintDelta`). `lhs()` / `rhs()` read the side values.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ConstraintSpan {
    /// In the interface `K` — present and identical on both sides.
    Unchanged(Constraint),
    /// In `R` only — created.
    Added(Constraint),
    /// In `L` only — deleted.
    Removed(Constraint),
}

impl ConstraintSpan {
    /// The lhs (`L`) value, or `None` if the constraint is created.
    pub fn lhs(&self) -> Option<&Constraint> {
        match self {
            Self::Unchanged(value) | Self::Removed(value) => Some(value),
            Self::Added(_) => None,
        }
    }

    /// The rhs (`R`) value, or `None` if the constraint is deleted.
    pub fn rhs(&self) -> Option<&Constraint> {
        match self {
            Self::Unchanged(value) | Self::Added(value) => Some(value),
            Self::Removed(_) => None,
        }
    }

    pub(crate) fn collect_frame_action_domain(&self, domain: &mut ConstraintFrameActionDomain) {
        match self {
            Self::Unchanged(constraint) | Self::Added(constraint) | Self::Removed(constraint) => {
                constraint.collect_frame_action_domain(domain);
            }
        }
    }

    pub(crate) fn reframe_by_actions(
        self,
        actions: &impl ConstraintFrameActions,
    ) -> Result<Self, Entity> {
        Ok(match self {
            Self::Unchanged(constraint) => Self::Unchanged(constraint.reframe_by_actions(actions)?),
            Self::Added(constraint) => Self::Added(constraint.reframe_by_actions(actions)?),
            Self::Removed(constraint) => Self::Removed(constraint.reframe_by_actions(actions)?),
        })
    }
}

impl Normalize for ConstraintSpan {
    fn normalize(self) -> Result<Self, Contradiction> {
        Ok(match self {
            Self::Unchanged(value) => Self::Unchanged(value.normalize()?),
            Self::Added(value) => Self::Added(value.normalize()?),
            Self::Removed(value) => Self::Removed(value.normalize()?),
        })
    }
}

impl FrameTransport for ConstraintSpan {
    type Action = OverlaysFrameAction;

    fn reframe_by(self, actions: &Self::Action) -> Option<Self> {
        self.reframe_by_actions(actions).ok()
    }
}

/// The per-entity op the fold operates on, abstracting `AtomDelta`/`BondDelta`. `atoms` carries
/// the entity's participant atoms in `Add`/`Remove` (`()` for an atom, its two ids for a bond).
pub(crate) enum EntityOp<F: EntityFold> {
    Add {
        atoms: F::Atoms,
        attributes: F::Attributes,
    },
    Remove {
        atoms: F::Atoms,
        attributes: F::Attributes,
    },
    ModifyField(F::FieldChange),
    ModifyConstraint {
        old: Option<F::Constraint>,
        new: Option<F::Constraint>,
    },
}

/// The normalize-fold extension of `EntityPatch` — the `EntityOp` deconstruction
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

    /// Recover this kind's deltas from its lhs/rhs state column: `Added`/`Removed` become
    /// structural `Add`/`Remove` (their `atoms` from `atoms(index)`), `Modified` becomes the
    /// field/constraint `diff`. `id(index)` selects the output id for each union-frame entry.
    fn append_deltas_from_states(
        states: &[EntitySpan<Self::Attributes>],
        id: impl Fn(usize) -> Self::Id,
        atoms: impl Fn(usize) -> Self::Atoms,
        deltas: &mut Deltas,
    ) {
        for (index, state) in states.iter().enumerate() {
            let id = id(index);
            match state {
                EntitySpan::Unchanged(_) => {}
                EntitySpan::Added(attributes) => deltas.push(
                    Self::rebuild(
                        id,
                        EntityOp::Add {
                            atoms: atoms(index),
                            attributes: attributes.clone(),
                        },
                    )
                    .into_delta(),
                ),
                EntitySpan::Removed(attributes) => deltas.push(
                    Self::rebuild(
                        id,
                        EntityOp::Remove {
                            atoms: atoms(index),
                            attributes: attributes.clone(),
                        },
                    )
                    .into_delta(),
                ),
                EntitySpan::Modified { lhs, rhs } => {
                    for delta in Self::diff(id, lhs, rhs) {
                        deltas.push(Self::into_delta(delta));
                    }
                }
            }
        }
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

/// Created entity: seed `attributes` from `Add`, absorb subsequent field/constraint changes; an
/// `Add`+`Remove` cancels. Yields one `Add` with the final attributes, or nothing.
fn fold_created<F: EntityFold>(ops: Vec<EntityOp<F>>) -> Result<Vec<EntityOp<F>>, Contradiction> {
    let mut state: Option<(F::Atoms, F::Attributes)> = None;
    let mut removed = false;
    for op in ops {
        if removed {
            return Err(Contradiction);
        }
        match op {
            EntityOp::Add { atoms, attributes } => {
                if state.is_some() {
                    return Err(Contradiction);
                }
                state = Some((atoms, attributes));
            }
            EntityOp::ModifyField(change) => {
                let (_, attributes) = state.as_mut().ok_or(Contradiction)?;
                F::apply_field(attributes, change)?;
            }
            EntityOp::ModifyConstraint { old, new } => {
                let (_, attributes) = state.as_mut().ok_or(Contradiction)?;
                F::apply_constraint(attributes, old, new)?;
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
        Some((atoms, attributes)) => vec![EntityOp::Add { atoms, attributes }],
        None => Vec::new(),
    })
}

/// Preserved entity: fuse `ModifyField` chains per field and `ModifyConstraint` chains per key. A
/// `Remove` subsumes the prior changes and carries the *original* value (the changes are
/// reverted on the removed attributes).
#[allow(clippy::type_complexity)]
fn fold_preserved<F: EntityFold>(ops: Vec<EntityOp<F>>) -> Result<Vec<EntityOp<F>>, Contradiction> {
    let mut fields: HashMap<Discriminant<F::FieldChange>, F::FieldChange> = HashMap::new();
    let mut constraints: HashMap<F::ConstraintKey, (Option<F::Constraint>, Option<F::Constraint>)> =
        HashMap::new();
    let mut removed: Option<(F::Atoms, F::Attributes)> = None;
    for op in ops {
        if removed.is_some() {
            return Err(Contradiction);
        }
        match op {
            EntityOp::Add { .. } => return Err(Contradiction),
            EntityOp::ModifyField(change) => {
                let field_key = discriminant(&change);
                let fused = match fields.remove(&field_key) {
                    Some(prev) => F::fuse_field(prev, change).ok_or(Contradiction)?,
                    None => change,
                };
                fields.insert(field_key, fused);
            }
            EntityOp::ModifyConstraint { old, new } => {
                let key = match old.as_ref().or(new.as_ref()) {
                    Some(constraint) => F::constraint_key(constraint),
                    None => continue,
                };
                match constraints.remove(&key) {
                    Some((first_old, prev_new)) => {
                        if !options_normalized_eq(&prev_new, &old) {
                            return Err(Contradiction);
                        }
                        constraints.insert(key, (first_old, new));
                    }
                    None => {
                        constraints.insert(key, (old, new));
                    }
                }
            }
            EntityOp::Remove { atoms, attributes } => {
                removed = Some((atoms, attributes));
            }
        }
    }
    if let Some((atoms, mut attributes)) = removed {
        for (_idx, change) in fields {
            F::apply_field(&mut attributes, F::field_inverse(change))?;
        }
        for (_key, (old, new)) in constraints {
            F::apply_constraint(&mut attributes, new, old)?;
        }
        return Ok(vec![EntityOp::Remove { atoms, attributes }]);
    }
    let mut out = Vec::new();
    for (_idx, change) in fields {
        if !F::field_is_identity(&change) {
            out.push(EntityOp::ModifyField(change));
        }
    }
    for (_key, (old, new)) in constraints {
        if !options_normalized_eq(&old, &new) {
            out.push(EntityOp::ModifyConstraint { old, new });
        }
    }
    Ok(out)
}

impl EntityPatch for AtomDelta {
    type Id = AtomId;
    type Attributes = AtomForm;
    type FieldChange = AtomFieldChange;
    type Constraint = AtomConstraintForm;

    fn modify_field(id: AtomId, change: AtomFieldChange) -> Self {
        AtomDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: AtomId,
        old: Option<AtomConstraintForm>,
        new: Option<AtomConstraintForm>,
    ) -> Self {
        AtomDelta::ModifyConstraint { id, old, new }
    }

    fn diff(id: AtomId, lhs: &AtomForm, rhs: &AtomForm) -> Vec<Self> {
        Self::for_update(id, lhs, &lhs.difference_to(rhs))
    }

    diff_field_ops!(AtomFieldChange, AtomForm, AtomConstraintForm, {
        Element => element,
        IsotopeMass => isotope_mass,
        Charge => charge,
        ImplicitHydrogens => implicit_hydrogens,
        LonePairs => lone_pairs,
        UnpairedElectrons => unpaired_electrons,
    });

    fn apply_constraint(
        attributes: &mut AtomForm,
        old: Option<AtomConstraintForm>,
        new: Option<AtomConstraintForm>,
    ) -> Result<(), Contradiction> {
        attributes.constraints.compare_and_set(old, new)
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
            AtomDelta::Add { attributes, .. } => EntityOp::Add {
                atoms: (),
                attributes,
            },
            AtomDelta::Remove { attributes, .. } => EntityOp::Remove {
                atoms: (),
                attributes,
            },
            AtomDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            AtomDelta::ModifyConstraint { old, new, .. } => EntityOp::ModifyConstraint { old, new },
        }
    }

    fn rebuild(id: AtomId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { attributes, .. } => AtomDelta::Add { id, attributes },
            EntityOp::Remove { attributes, .. } => AtomDelta::Remove { id, attributes },
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

    fn constraint_key(constraint: &AtomConstraintForm) -> AtomConstraintKey {
        constraint.key()
    }

    fold_field_ops!(AtomFieldChange, {
        Element,
        IsotopeMass,
        Charge,
        ImplicitHydrogens,
        LonePairs,
        UnpairedElectrons,
    });
}

impl EntityPatch for BondDelta {
    type Id = BondId;
    type Attributes = BondForm;
    type FieldChange = BondFieldChange;
    type Constraint = BondConstraintForm;

    fn modify_field(id: BondId, change: BondFieldChange) -> Self {
        BondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: BondId,
        old: Option<BondConstraintForm>,
        new: Option<BondConstraintForm>,
    ) -> Self {
        BondDelta::ModifyConstraint { id, old, new }
    }

    fn diff(id: BondId, lhs: &BondForm, rhs: &BondForm) -> Vec<Self> {
        Self::for_update(id, lhs, &lhs.difference_to(rhs))
    }

    diff_field_ops!(BondFieldChange, BondForm, BondConstraintForm, {
        Order => order,
        Charge => charge,
        UnpairedElectrons => unpaired_electrons,
    });

    fn apply_constraint(
        attributes: &mut BondForm,
        old: Option<BondConstraintForm>,
        new: Option<BondConstraintForm>,
    ) -> Result<(), Contradiction> {
        attributes.constraints.compare_and_set(old, new)
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
            BondDelta::Add {
                atoms, attributes, ..
            } => EntityOp::Add { atoms, attributes },
            BondDelta::Remove {
                atoms, attributes, ..
            } => EntityOp::Remove { atoms, attributes },
            BondDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            BondDelta::ModifyConstraint { old, new, .. } => EntityOp::ModifyConstraint { old, new },
        }
    }

    fn rebuild(id: BondId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { atoms, attributes } => BondDelta::Add {
                id,
                atoms,
                attributes,
            },
            EntityOp::Remove { atoms, attributes } => BondDelta::Remove {
                id,
                atoms,
                attributes,
            },
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

    fn constraint_key(constraint: &BondConstraintForm) -> BondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(BondFieldChange, {
        Order,
        Charge,
        UnpairedElectrons,
    });
}

impl EntityPatch for DativeBondDelta {
    type Id = DativeBondId;
    type Attributes = DativeBondForm;
    type FieldChange = DativeBondFieldChange;
    type Constraint = DativeBondConstraintForm;

    fn modify_field(id: DativeBondId, change: DativeBondFieldChange) -> Self {
        DativeBondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: DativeBondId,
        old: Option<DativeBondConstraintForm>,
        new: Option<DativeBondConstraintForm>,
    ) -> Self {
        DativeBondDelta::ModifyConstraint { id, old, new }
    }

    fn diff(id: DativeBondId, lhs: &DativeBondForm, rhs: &DativeBondForm) -> Vec<Self> {
        Self::for_update(id, lhs, &lhs.difference_to(rhs))
    }

    diff_field_ops!(DativeBondFieldChange, DativeBondForm, DativeBondConstraintForm, {
        Order => order,
    });

    fn apply_constraint(
        attributes: &mut DativeBondForm,
        old: Option<DativeBondConstraintForm>,
        new: Option<DativeBondConstraintForm>,
    ) -> Result<(), Contradiction> {
        attributes.constraints.compare_and_set(old, new)
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
                attributes,
                ..
            } => EntityOp::Add {
                atoms: (donors, acceptor),
                attributes,
            },
            DativeBondDelta::Remove {
                donors,
                acceptor,
                attributes,
                ..
            } => EntityOp::Remove {
                atoms: (donors, acceptor),
                attributes,
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
                attributes,
            } => DativeBondDelta::Add {
                id,
                donors,
                acceptor,
                attributes,
            },
            EntityOp::Remove {
                atoms: (donors, acceptor),
                attributes,
            } => DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                attributes,
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

    fn constraint_key(constraint: &DativeBondConstraintForm) -> DativeBondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(DativeBondFieldChange, { Order });
}

impl EntityPatch for AromaticSystemDelta {
    type Id = AromaticSystemId;
    type Attributes = AromaticSystemForm;
    type FieldChange = AromaticSystemFieldChange;
    type Constraint = AromaticSystemConstraintForm;

    fn modify_field(id: AromaticSystemId, change: AromaticSystemFieldChange) -> Self {
        AromaticSystemDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: AromaticSystemId,
        old: Option<AromaticSystemConstraintForm>,
        new: Option<AromaticSystemConstraintForm>,
    ) -> Self {
        AromaticSystemDelta::ModifyConstraint { id, old, new }
    }

    fn diff(id: AromaticSystemId, lhs: &AromaticSystemForm, rhs: &AromaticSystemForm) -> Vec<Self> {
        Self::for_update(id, lhs, &lhs.difference_to(rhs))
    }

    diff_field_ops!(
        AromaticSystemFieldChange,
        AromaticSystemForm,
        AromaticSystemConstraintForm,
        {
            Electrons => electrons,
            Charge => charge,
            UnpairedElectrons => unpaired_electrons,
        }
    );

    fn apply_constraint(
        attributes: &mut AromaticSystemForm,
        old: Option<AromaticSystemConstraintForm>,
        new: Option<AromaticSystemConstraintForm>,
    ) -> Result<(), Contradiction> {
        attributes.constraints.compare_and_set(old, new)
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
            AromaticSystemDelta::Add {
                atoms, attributes, ..
            } => EntityOp::Add { atoms, attributes },
            AromaticSystemDelta::Remove {
                atoms, attributes, ..
            } => EntityOp::Remove { atoms, attributes },
            AromaticSystemDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            AromaticSystemDelta::ModifyConstraint { old, new, .. } => {
                EntityOp::ModifyConstraint { old, new }
            }
        }
    }

    fn rebuild(id: AromaticSystemId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { atoms, attributes } => AromaticSystemDelta::Add {
                id,
                atoms,
                attributes,
            },
            EntityOp::Remove { atoms, attributes } => AromaticSystemDelta::Remove {
                id,
                atoms,
                attributes,
            },
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

    fn constraint_key(constraint: &AromaticSystemConstraintForm) -> AromaticSystemConstraintKey {
        constraint.key()
    }

    fold_field_ops!(AromaticSystemFieldChange, {
        Electrons,
        Charge,
        UnpairedElectrons
    });
}

impl EntityPatch for MulticenterBondDelta {
    type Id = MulticenterBondId;
    type Attributes = MulticenterBondForm;
    type FieldChange = MulticenterBondFieldChange;
    type Constraint = MulticenterBondConstraintForm;

    fn modify_field(id: MulticenterBondId, change: MulticenterBondFieldChange) -> Self {
        MulticenterBondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: MulticenterBondId,
        old: Option<MulticenterBondConstraintForm>,
        new: Option<MulticenterBondConstraintForm>,
    ) -> Self {
        MulticenterBondDelta::ModifyConstraint { id, old, new }
    }

    fn diff(
        id: MulticenterBondId,
        lhs: &MulticenterBondForm,
        rhs: &MulticenterBondForm,
    ) -> Vec<Self> {
        Self::for_update(id, lhs, &lhs.difference_to(rhs))
    }

    diff_field_ops!(
        MulticenterBondFieldChange,
        MulticenterBondForm,
        MulticenterBondConstraintForm,
        {
            Electrons => electrons,
            Charge => charge,
            UnpairedElectrons => unpaired_electrons,
        }
    );

    fn apply_constraint(
        attributes: &mut MulticenterBondForm,
        old: Option<MulticenterBondConstraintForm>,
        new: Option<MulticenterBondConstraintForm>,
    ) -> Result<(), Contradiction> {
        attributes.constraints.compare_and_set(old, new)
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
            MulticenterBondDelta::Add {
                atoms, attributes, ..
            } => EntityOp::Add { atoms, attributes },
            MulticenterBondDelta::Remove {
                atoms, attributes, ..
            } => EntityOp::Remove { atoms, attributes },
            MulticenterBondDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            MulticenterBondDelta::ModifyConstraint { old, new, .. } => {
                EntityOp::ModifyConstraint { old, new }
            }
        }
    }

    fn rebuild(id: MulticenterBondId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { atoms, attributes } => MulticenterBondDelta::Add {
                id,
                atoms,
                attributes,
            },
            EntityOp::Remove { atoms, attributes } => MulticenterBondDelta::Remove {
                id,
                atoms,
                attributes,
            },
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

    fn constraint_key(constraint: &MulticenterBondConstraintForm) -> MulticenterBondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(MulticenterBondFieldChange, {
        Electrons,
        Charge,
        UnpairedElectrons
    });
}

impl EntityPatch for NoncovalentBondDelta {
    type Id = NoncovalentBondId;
    type Attributes = NoncovalentBondForm;
    type FieldChange = NoncovalentBondFieldChange;
    type Constraint = NoncovalentBondConstraintForm;

    fn modify_field(id: NoncovalentBondId, change: NoncovalentBondFieldChange) -> Self {
        NoncovalentBondDelta::ModifyField { id, change }
    }

    fn modify_constraint(
        id: NoncovalentBondId,
        old: Option<NoncovalentBondConstraintForm>,
        new: Option<NoncovalentBondConstraintForm>,
    ) -> Self {
        NoncovalentBondDelta::ModifyConstraint { id, old, new }
    }

    fn diff(
        id: NoncovalentBondId,
        lhs: &NoncovalentBondForm,
        rhs: &NoncovalentBondForm,
    ) -> Vec<Self> {
        Self::for_update(id, lhs, &lhs.difference_to(rhs))
    }

    diff_field_ops!(
        NoncovalentBondFieldChange,
        NoncovalentBondForm,
        NoncovalentBondConstraintForm,
        {
            Kind => kind,
        }
    );

    fn apply_constraint(
        attributes: &mut NoncovalentBondForm,
        old: Option<NoncovalentBondConstraintForm>,
        new: Option<NoncovalentBondConstraintForm>,
    ) -> Result<(), Contradiction> {
        attributes.constraints.compare_and_set(old, new)
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
            NoncovalentBondDelta::Add {
                atoms, attributes, ..
            } => EntityOp::Add { atoms, attributes },
            NoncovalentBondDelta::Remove {
                atoms, attributes, ..
            } => EntityOp::Remove { atoms, attributes },
            NoncovalentBondDelta::ModifyField { change, .. } => EntityOp::ModifyField(change),
            NoncovalentBondDelta::ModifyConstraint { old, new, .. } => {
                EntityOp::ModifyConstraint { old, new }
            }
        }
    }

    fn rebuild(id: NoncovalentBondId, op: EntityOp<Self>) -> Self {
        match op {
            EntityOp::Add { atoms, attributes } => NoncovalentBondDelta::Add {
                id,
                atoms,
                attributes,
            },
            EntityOp::Remove { atoms, attributes } => NoncovalentBondDelta::Remove {
                id,
                atoms,
                attributes,
            },
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

    fn constraint_key(constraint: &NoncovalentBondConstraintForm) -> NoncovalentBondConstraintKey {
        constraint.key()
    }

    fold_field_ops!(NoncovalentBondFieldChange, { Kind });
}

impl EntityPatch for StereoAtomDelta {
    type Id = StereoAtomId;
    type Attributes = StereoAtomForm;
    type FieldChange = StereoAtomFieldChange;
    type Constraint = StereoAtomConstraintForm;

    fn modify_field(id: StereoAtomId, change: StereoAtomFieldChange) -> Self {
        StereoAtomDelta::ModifyField { id, change }
    }

    /// Kind-less fallback (the trait signature has no kind); the real producer is the overridden
    /// `diff`, which stamps `kind` from the entity's config. Stereo's flow never uses this arm.
    fn modify_constraint(
        id: StereoAtomId,
        old: Option<StereoAtomConstraintForm>,
        new: Option<StereoAtomConstraintForm>,
    ) -> Self {
        StereoAtomDelta::ModifyConstraint {
            id,
            kind: None,
            old,
            new,
        }
    }

    diff_field_ops!(StereoAtomFieldChange, StereoAtomForm, StereoAtomConstraintForm, {
        Configuration => configuration,
    });

    fn diff(id: StereoAtomId, lhs: &StereoAtomForm, rhs: &StereoAtomForm) -> Vec<Self> {
        Self::for_update(id, lhs, &lhs.difference_to(rhs))
    }

    fn apply_constraint(
        attributes: &mut StereoAtomForm,
        old: Option<StereoAtomConstraintForm>,
        new: Option<StereoAtomConstraintForm>,
    ) -> Result<(), Contradiction> {
        attributes.constraints.compare_and_set(old, new)
    }
}

impl EntityPatch for StereoBondDelta {
    type Id = StereoBondId;
    type Attributes = StereoBondForm;
    type FieldChange = StereoBondFieldChange;
    type Constraint = StereoBondConstraintForm;

    fn modify_field(id: StereoBondId, change: StereoBondFieldChange) -> Self {
        StereoBondDelta::ModifyField { id, change }
    }

    /// Kind-less fallback — see `StereoAtomDelta::modify_constraint`.
    fn modify_constraint(
        id: StereoBondId,
        old: Option<StereoBondConstraintForm>,
        new: Option<StereoBondConstraintForm>,
    ) -> Self {
        StereoBondDelta::ModifyConstraint {
            id,
            kind: None,
            old,
            new,
        }
    }

    diff_field_ops!(StereoBondFieldChange, StereoBondForm, StereoBondConstraintForm, {
        Configuration => configuration,
    });

    fn diff(id: StereoBondId, lhs: &StereoBondForm, rhs: &StereoBondForm) -> Vec<Self> {
        Self::for_update(id, lhs, &lhs.difference_to(rhs))
    }

    fn apply_constraint(
        attributes: &mut StereoBondForm,
        old: Option<StereoBondConstraintForm>,
        new: Option<StereoBondConstraintForm>,
    ) -> Result<(), Contradiction> {
        attributes.constraints.compare_and_set(old, new)
    }
}

/// Apply a resolved per-entity change to a form, reusing the `EntityPatch` apply that
/// `normalize` uses. `ModifyField` / `ModifyConstraint` mutate the attributes; `Add` / `Remove` are
/// no-ops (they carry a whole attributes, not a change). Materializes the rhs-hand value of a
/// preserved entity for a `ReactionSpan`.
pub(crate) fn apply_atom_change(
    attributes: &mut AtomForm,
    delta: &AtomDelta,
) -> Result<(), Contradiction> {
    match delta {
        AtomDelta::ModifyField { change, .. } => AtomDelta::apply_field(attributes, change.clone()),
        AtomDelta::ModifyConstraint { old, new, .. } => {
            AtomDelta::apply_constraint(attributes, old.clone(), new.clone())
        }
        AtomDelta::Add { .. } | AtomDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_bond_change(
    attributes: &mut BondForm,
    delta: &BondDelta,
) -> Result<(), Contradiction> {
    match delta {
        BondDelta::ModifyField { change, .. } => BondDelta::apply_field(attributes, change.clone()),
        BondDelta::ModifyConstraint { old, new, .. } => {
            BondDelta::apply_constraint(attributes, old.clone(), new.clone())
        }
        BondDelta::Add { .. } | BondDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_dative_change(
    attributes: &mut DativeBondForm,
    delta: &DativeBondDelta,
) -> Result<(), Contradiction> {
    match delta {
        DativeBondDelta::ModifyField { change, .. } => {
            DativeBondDelta::apply_field(attributes, change.clone())
        }
        DativeBondDelta::ModifyConstraint { old, new, .. } => {
            DativeBondDelta::apply_constraint(attributes, old.clone(), new.clone())
        }
        DativeBondDelta::Add { .. } | DativeBondDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_aromatic_change(
    attributes: &mut AromaticSystemForm,
    delta: &AromaticSystemDelta,
) -> Result<(), Contradiction> {
    match delta {
        AromaticSystemDelta::ModifyField { change, .. } => {
            AromaticSystemDelta::apply_field(attributes, change.clone())
        }
        AromaticSystemDelta::ModifyConstraint { old, new, .. } => {
            AromaticSystemDelta::apply_constraint(attributes, old.clone(), new.clone())
        }
        AromaticSystemDelta::Add { .. } | AromaticSystemDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_multicenter_change(
    attributes: &mut MulticenterBondForm,
    delta: &MulticenterBondDelta,
) -> Result<(), Contradiction> {
    match delta {
        MulticenterBondDelta::ModifyField { change, .. } => {
            MulticenterBondDelta::apply_field(attributes, change.clone())
        }
        MulticenterBondDelta::ModifyConstraint { old, new, .. } => {
            MulticenterBondDelta::apply_constraint(attributes, old.clone(), new.clone())
        }
        MulticenterBondDelta::Add { .. } | MulticenterBondDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_noncovalent_change(
    attributes: &mut NoncovalentBondForm,
    delta: &NoncovalentBondDelta,
) -> Result<(), Contradiction> {
    match delta {
        NoncovalentBondDelta::ModifyField { change, .. } => {
            NoncovalentBondDelta::apply_field(attributes, change.clone())
        }
        NoncovalentBondDelta::ModifyConstraint { old, new, .. } => {
            NoncovalentBondDelta::apply_constraint(attributes, old.clone(), new.clone())
        }
        NoncovalentBondDelta::Add { .. } | NoncovalentBondDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_stereo_atom_change(
    attributes: &mut StereoAtomForm,
    delta: &StereoAtomDelta,
) -> Result<(), Contradiction> {
    match delta {
        StereoAtomDelta::ModifyField { change, .. } => {
            StereoAtomDelta::apply_field(attributes, change.clone())
        }
        StereoAtomDelta::ModifyConstraint { old, new, .. } => {
            StereoAtomDelta::apply_constraint(attributes, old.clone(), new.clone())
        }
        StereoAtomDelta::Add { .. } | StereoAtomDelta::Remove { .. } => Ok(()),
    }
}

pub(crate) fn apply_stereo_bond_change(
    attributes: &mut StereoBondForm,
    delta: &StereoBondDelta,
) -> Result<(), Contradiction> {
    match delta {
        StereoBondDelta::ModifyField { change, .. } => {
            StereoBondDelta::apply_field(attributes, change.clone())
        }
        StereoBondDelta::ModifyConstraint { old, new, .. } => {
            StereoBondDelta::apply_constraint(attributes, old.clone(), new.clone())
        }
        StereoBondDelta::Add { .. } | StereoBondDelta::Remove { .. } => Ok(()),
    }
}

enum StereoConfigFold {
    Identity,
    Set {
        old: StereoConfigurationForm,
        new: StereoConfigurationForm,
    },
}

fn fold_stereo_config(
    changes: Vec<(StereoConfigurationForm, StereoConfigurationForm)>,
) -> Result<StereoConfigFold, Contradiction> {
    let mut changes = changes.into_iter();
    let Some((old, mut new)) = changes.next() else {
        return Ok(StereoConfigFold::Identity);
    };
    for (next_old, next_new) in changes {
        if new.clone().normalize()? != next_old.clone().normalize()? {
            return Err(Contradiction);
        }
        new = next_new;
    }
    if old.clone().normalize()? == new.clone().normalize()? {
        Ok(StereoConfigFold::Identity)
    } else {
        Ok(StereoConfigFold::Set { old, new })
    }
}

/// Fold one stereo atom's deltas to normal form in input order.
fn fold_stereo_atom_group(
    id: StereoAtomId,
    group: Vec<StereoAtomDelta>,
) -> Result<Vec<StereoAtomDelta>, Contradiction> {
    if group
        .iter()
        .any(|d| matches!(d, StereoAtomDelta::Add { .. }))
    {
        let mut state: Option<(AtomId, Vec<StereoLigand>, StereoAtomForm)> = None;
        let mut removed = false;
        for delta in group {
            if removed {
                return Err(Contradiction);
            }
            match delta {
                StereoAtomDelta::Add {
                    site,
                    ligands,
                    attributes,
                    ..
                } => {
                    if state.is_some() {
                        return Err(Contradiction);
                    }
                    state = Some((site, ligands, attributes));
                }
                StereoAtomDelta::Remove { .. } => {
                    if state.is_none() {
                        return Err(Contradiction);
                    }
                    state = None;
                    removed = true;
                }
                other => {
                    let (_, _, attributes) = state.as_mut().ok_or(Contradiction)?;
                    apply_stereo_atom_change(attributes, &other)?;
                }
            }
        }
        return Ok(match state {
            Some((site, ligands, attributes)) => vec![StereoAtomDelta::Add {
                id,
                site,
                ligands,
                attributes,
            }],
            None => Vec::new(),
        });
    }
    let mut kind: Option<StereoKind> = None;
    let mut config_changes = Vec::new();
    let mut constraints: HashMap<
        StereoAtomConstraintKey,
        (
            Option<StereoAtomConstraintForm>,
            Option<StereoAtomConstraintForm>,
        ),
    > = HashMap::new();
    let mut removed: Option<(AtomId, Vec<StereoLigand>, StereoAtomForm)> = None;
    for delta in group {
        if removed.is_some() {
            return Err(Contradiction);
        }
        match delta {
            StereoAtomDelta::Add { .. } => return Err(Contradiction),
            StereoAtomDelta::ModifyField {
                change: StereoAtomFieldChange::Configuration { old, new },
                ..
            } => {
                kind = kind.or_else(|| old.kind());
                config_changes.push((old, new));
            }
            StereoAtomDelta::ModifyConstraint {
                kind: constraint_kind,
                old,
                new,
                ..
            } => {
                kind = kind.or(constraint_kind);
                let key = match old.as_ref().or(new.as_ref()) {
                    Some(c) => c.key(),
                    None => continue,
                };
                match constraints.remove(&key) {
                    Some((first_old, prev_new)) => {
                        if !options_normalized_eq(&prev_new, &old) {
                            return Err(Contradiction);
                        }
                        constraints.insert(key, (first_old, new));
                    }
                    None => {
                        constraints.insert(key, (old, new));
                    }
                }
            }
            StereoAtomDelta::Remove {
                site,
                ligands,
                attributes,
                ..
            } => {
                removed = Some((site, ligands, attributes));
            }
        }
    }
    let config = fold_stereo_config(config_changes)?;
    if let Some((site, ligands, mut attributes)) = removed {
        match config {
            StereoConfigFold::Identity => {}
            StereoConfigFold::Set { old, new } => {
                if attributes.configuration.clone().normalize()? != new.clone().normalize()? {
                    return Err(Contradiction);
                }
                attributes.configuration = old;
            }
        }
        for (_key, (old, new)) in constraints {
            StereoAtomDelta::apply_constraint(&mut attributes, new, old)?;
        }
        return Ok(vec![StereoAtomDelta::Remove {
            id,
            site,
            ligands,
            attributes,
        }]);
    }
    let mut out = Vec::new();
    match config {
        StereoConfigFold::Identity => {}
        StereoConfigFold::Set { old, new } => out.push(StereoAtomDelta::ModifyField {
            id,
            change: StereoAtomFieldChange::Configuration { old, new },
        }),
    }
    for (_key, (old, new)) in constraints {
        if !options_normalized_eq(&old, &new) {
            out.push(StereoAtomDelta::ModifyConstraint { id, kind, old, new });
        }
    }
    Ok(out)
}

/// Fold one stereo bond's deltas to normal form — the `fold_stereo_atom_group` twin (bond ids/attributes).
fn fold_stereo_bond_group(
    id: StereoBondId,
    group: Vec<StereoBondDelta>,
) -> Result<Vec<StereoBondDelta>, Contradiction> {
    if group
        .iter()
        .any(|d| matches!(d, StereoBondDelta::Add { .. }))
    {
        let mut state: Option<(BondId, Vec<StereoLigand>, StereoBondForm)> = None;
        let mut removed = false;
        for delta in group {
            if removed {
                return Err(Contradiction);
            }
            match delta {
                StereoBondDelta::Add {
                    site,
                    ligands,
                    attributes,
                    ..
                } => {
                    if state.is_some() {
                        return Err(Contradiction);
                    }
                    state = Some((site, ligands, attributes));
                }
                StereoBondDelta::Remove { .. } => {
                    if state.is_none() {
                        return Err(Contradiction);
                    }
                    state = None;
                    removed = true;
                }
                other => {
                    let (_, _, attributes) = state.as_mut().ok_or(Contradiction)?;
                    apply_stereo_bond_change(attributes, &other)?;
                }
            }
        }
        return Ok(match state {
            Some((site, ligands, attributes)) => vec![StereoBondDelta::Add {
                id,
                site,
                ligands,
                attributes,
            }],
            None => Vec::new(),
        });
    }
    let mut kind: Option<StereoKind> = None;
    let mut config_changes = Vec::new();
    let mut constraints: HashMap<
        StereoBondConstraintKey,
        (
            Option<StereoBondConstraintForm>,
            Option<StereoBondConstraintForm>,
        ),
    > = HashMap::new();
    let mut removed: Option<(BondId, Vec<StereoLigand>, StereoBondForm)> = None;
    for delta in group {
        if removed.is_some() {
            return Err(Contradiction);
        }
        match delta {
            StereoBondDelta::Add { .. } => return Err(Contradiction),
            StereoBondDelta::ModifyField {
                change: StereoBondFieldChange::Configuration { old, new },
                ..
            } => {
                kind = kind.or_else(|| old.kind());
                config_changes.push((old, new));
            }
            StereoBondDelta::ModifyConstraint {
                kind: constraint_kind,
                old,
                new,
                ..
            } => {
                kind = kind.or(constraint_kind);
                let key = match old.as_ref().or(new.as_ref()) {
                    Some(c) => c.key(),
                    None => continue,
                };
                match constraints.remove(&key) {
                    Some((first_old, prev_new)) => {
                        if !options_normalized_eq(&prev_new, &old) {
                            return Err(Contradiction);
                        }
                        constraints.insert(key, (first_old, new));
                    }
                    None => {
                        constraints.insert(key, (old, new));
                    }
                }
            }
            StereoBondDelta::Remove {
                site,
                ligands,
                attributes,
                ..
            } => {
                removed = Some((site, ligands, attributes));
            }
        }
    }
    let config = fold_stereo_config(config_changes)?;
    if let Some((site, ligands, mut attributes)) = removed {
        match config {
            StereoConfigFold::Identity => {}
            StereoConfigFold::Set { old, new } => {
                if attributes.configuration.clone().normalize()? != new.clone().normalize()? {
                    return Err(Contradiction);
                }
                attributes.configuration = old;
            }
        }
        for (_key, (old, new)) in constraints {
            StereoBondDelta::apply_constraint(&mut attributes, new, old)?;
        }
        return Ok(vec![StereoBondDelta::Remove {
            id,
            site,
            ligands,
            attributes,
        }]);
    }
    let mut out = Vec::new();
    match config {
        StereoConfigFold::Identity => {}
        StereoConfigFold::Set { old, new } => out.push(StereoBondDelta::ModifyField {
            id,
            change: StereoBondFieldChange::Configuration { old, new },
        }),
    }
    for (_key, (old, new)) in constraints {
        if !options_normalized_eq(&old, &new) {
            out.push(StereoBondDelta::ModifyConstraint { id, kind, old, new });
        }
    }
    Ok(out)
}

/// Re-anchor a delta's ids and participant atoms through a total id relabeling. Used to move
/// deltas between id spaces (reverse re-anchoring, composition). The relabeling must cover every id
/// the delta references. Participant sequences and frame-relative payloads retain their supplied
/// frame; selecting another frame is a separate [`FrameTransport`] operation.
#[cfg(test)]
fn remap_delta(delta: Delta, map: &IdRemapping) -> Delta {
    match delta {
        Delta::Atom(a) => Delta::Atom(match a {
            AtomDelta::Add { id, attributes } => AtomDelta::Add {
                id: map.map_atom(id),
                attributes,
            },
            AtomDelta::Remove { id, attributes } => AtomDelta::Remove {
                id: map.map_atom(id),
                attributes,
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
            BondDelta::Add {
                id,
                atoms,
                attributes,
            } => BondDelta::Add {
                id: map.map_bond(id),
                atoms: [map.map_atom(atoms[0]), map.map_atom(atoms[1])],
                attributes,
            },
            BondDelta::Remove {
                id,
                atoms,
                attributes,
            } => BondDelta::Remove {
                id: map.map_bond(id),
                atoms: [map.map_atom(atoms[0]), map.map_atom(atoms[1])],
                attributes,
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
                attributes,
            } => DativeBondDelta::Add {
                id: map.map_dative(id),
                donors: donors.into_iter().map(|atom| map.map_atom(atom)).collect(),
                acceptor: map.map_atom(acceptor),
                attributes,
            },
            DativeBondDelta::Remove {
                id,
                donors,
                acceptor,
                attributes,
            } => DativeBondDelta::Remove {
                id: map.map_dative(id),
                donors: donors.into_iter().map(|atom| map.map_atom(atom)).collect(),
                acceptor: map.map_atom(acceptor),
                attributes,
            },
            DativeBondDelta::ModifyField { id, change } => DativeBondDelta::ModifyField {
                id: map.map_dative(id),
                change,
            },
            DativeBondDelta::ModifyConstraint { id, old, new } => {
                DativeBondDelta::ModifyConstraint {
                    id: map.map_dative(id),
                    old,
                    new,
                }
            }
        }),
        Delta::AromaticSystem(a) => Delta::AromaticSystem(match a {
            AromaticSystemDelta::Add {
                id,
                atoms,
                attributes,
            } => AromaticSystemDelta::Add {
                id: map.map_aromatic(id),
                atoms: atoms.into_iter().map(|atom| map.map_atom(atom)).collect(),
                attributes,
            },
            AromaticSystemDelta::Remove {
                id,
                atoms,
                attributes,
            } => AromaticSystemDelta::Remove {
                id: map.map_aromatic(id),
                atoms: atoms.into_iter().map(|atom| map.map_atom(atom)).collect(),
                attributes,
            },
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
            MulticenterBondDelta::Add {
                id,
                atoms,
                attributes,
            } => MulticenterBondDelta::Add {
                id: map.map_multicenter(id),
                atoms: atoms.into_iter().map(|atom| map.map_atom(atom)).collect(),
                attributes,
            },
            MulticenterBondDelta::Remove {
                id,
                atoms,
                attributes,
            } => MulticenterBondDelta::Remove {
                id: map.map_multicenter(id),
                atoms: atoms.into_iter().map(|atom| map.map_atom(atom)).collect(),
                attributes,
            },
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
            NoncovalentBondDelta::Add {
                id,
                atoms,
                attributes,
            } => NoncovalentBondDelta::Add {
                id: map.map_noncovalent(id),
                atoms: [map.map_atom(atoms[0]), map.map_atom(atoms[1])],
                attributes,
            },
            NoncovalentBondDelta::Remove {
                id,
                atoms,
                attributes,
            } => NoncovalentBondDelta::Remove {
                id: map.map_noncovalent(id),
                atoms: [map.map_atom(atoms[0]), map.map_atom(atoms[1])],
                attributes,
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
        // Stereo: ids relabel and the ligand-frame atom ids relabel in place; the frame is `Ordered`
        // (not re-sorted on remap), so the coset stays valid and constraints remain position-local.
        Delta::StereoAtom(s) => Delta::StereoAtom(match s {
            StereoAtomDelta::Add {
                id,
                site,
                ligands,
                attributes,
            } => StereoAtomDelta::Add {
                id: map.map_stereo_atom(id),
                site: map.map_atom(site),
                ligands: ligands
                    .into_iter()
                    .map(|l| StereoLigand::new(map.map_atom(l.atom_id), l.kind))
                    .collect(),
                attributes,
            },
            StereoAtomDelta::Remove {
                id,
                site,
                ligands,
                attributes,
            } => StereoAtomDelta::Remove {
                id: map.map_stereo_atom(id),
                site: map.map_atom(site),
                ligands: ligands
                    .into_iter()
                    .map(|l| StereoLigand::new(map.map_atom(l.atom_id), l.kind))
                    .collect(),
                attributes,
            },
            StereoAtomDelta::ModifyField { id, change } => StereoAtomDelta::ModifyField {
                id: map.map_stereo_atom(id),
                change,
            },
            StereoAtomDelta::ModifyConstraint { id, kind, old, new } => {
                StereoAtomDelta::ModifyConstraint {
                    id: map.map_stereo_atom(id),
                    kind,
                    old,
                    new,
                }
            }
        }),
        Delta::StereoBond(s) => Delta::StereoBond(match s {
            StereoBondDelta::Add {
                id,
                site,
                ligands,
                attributes,
            } => StereoBondDelta::Add {
                id: map.map_stereo_bond(id),
                site: map.map_bond(site),
                ligands: ligands
                    .into_iter()
                    .map(|l| StereoLigand::new(map.map_atom(l.atom_id), l.kind))
                    .collect(),
                attributes,
            },
            StereoBondDelta::Remove {
                id,
                site,
                ligands,
                attributes,
            } => StereoBondDelta::Remove {
                id: map.map_stereo_bond(id),
                site: map.map_bond(site),
                ligands: ligands
                    .into_iter()
                    .map(|l| StereoLigand::new(map.map_atom(l.atom_id), l.kind))
                    .collect(),
                attributes,
            },
            StereoBondDelta::ModifyField { id, change } => StereoBondDelta::ModifyField {
                id: map.map_stereo_bond(id),
                change,
            },
            StereoBondDelta::ModifyConstraint { id, kind, old, new } => {
                StereoBondDelta::ModifyConstraint {
                    id: map.map_stereo_bond(id),
                    kind,
                    old,
                    new,
                }
            }
        }),
        Delta::Constraint(c) => Delta::Constraint(match c {
            ConstraintDelta::Add(constraint) => ConstraintDelta::Add(constraint.remap(map)),
            ConstraintDelta::Remove(constraint) => ConstraintDelta::Remove(constraint.remap(map)),
        }),
    }
}

/// A collection of [`Delta`] values, as carried by a reaction.
///
/// Before normalization, input order matters within a chain of operations on the same entity:
/// an addition must precede modifications to that addition, and successive field changes must
/// connect through their `old` and `new` values. Cross-entity source order is not semantic because
/// deltas refer to entities directly through ids in the reaction-owned frame.
/// [`Normalize::normalize`] folds each entity's chain, rejects contradictions, and sorts the
/// normalized result. Unlike an edit sequence, the normal form retains no incidental source
/// ordering.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

impl IntoIterator for Deltas {
    type Item = Delta;
    type IntoIter = IntoIter<Delta>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl Normalize for Deltas {
    /// Per-entity fold to the normal form, then a stable sort. Different entities are
    /// independent and each entity's fold is deterministic over input order, so the result is
    /// a unique normal form; sequence order is not stored. `Err(Contradiction)` on an
    /// inconsistent set.
    #[allow(clippy::mutable_key_type)]
    fn normalize(self) -> Result<Self, Contradiction> {
        let mut atoms: HashMap<AtomId, Vec<AtomDelta>> = HashMap::new();
        let mut bonds: HashMap<BondId, Vec<BondDelta>> = HashMap::new();
        let mut dative: HashMap<DativeBondId, Vec<DativeBondDelta>> = HashMap::new();
        let mut aromatic: HashMap<AromaticSystemId, Vec<AromaticSystemDelta>> = HashMap::new();
        let mut multicenter: HashMap<MulticenterBondId, Vec<MulticenterBondDelta>> = HashMap::new();
        let mut noncovalent: HashMap<NoncovalentBondId, Vec<NoncovalentBondDelta>> = HashMap::new();
        let mut stereo_atom: HashMap<StereoAtomId, Vec<StereoAtomDelta>> = HashMap::new();
        let mut stereo_bond: HashMap<StereoBondId, Vec<StereoBondDelta>> = HashMap::new();
        let mut constraints: Vec<ConstraintDelta> = Vec::new();
        for delta in self.0 {
            match delta {
                Delta::Atom(d) => atoms.entry(d.id()).or_default().push(d),
                Delta::Bond(d) => bonds.entry(d.id()).or_default().push(d),
                Delta::DativeBond(d) => dative.entry(d.id()).or_default().push(d),
                Delta::AromaticSystem(d) => aromatic.entry(d.id()).or_default().push(d),
                Delta::MulticenterBond(d) => multicenter.entry(d.id()).or_default().push(d),
                Delta::NoncovalentBond(d) => noncovalent.entry(d.id()).or_default().push(d),
                Delta::StereoAtom(d) => stereo_atom.entry(d.id()).or_default().push(d),
                Delta::StereoBond(d) => stereo_bond.entry(d.id()).or_default().push(d),
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
        // Overlay kinds: same fold; a created overlay must not reference a net-removed atom.
        for (id, group) in dative {
            let folded = fold_group::<DativeBondDelta>(id, group)?;
            for delta in &folded {
                if let DativeBondDelta::Add {
                    donors, acceptor, ..
                } = delta
                {
                    if donors
                        .iter()
                        .chain(iter::once(acceptor))
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
        // Stereo: bespoke fold (coset ops), same created-references-removed-atom guard (site +
        // ligands for a stereo atom; ligands for a stereo bond, whose site is a bond).
        for (id, group) in stereo_atom {
            let folded = fold_stereo_atom_group(id, group)?;
            for delta in &folded {
                if let StereoAtomDelta::Add { site, ligands, .. } = delta {
                    if iter::once(site)
                        .chain(ligands.iter().map(|l| &l.atom_id))
                        .any(|atom| removed_atoms.contains(atom))
                    {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::StereoAtom));
        }
        for (id, group) in stereo_bond {
            let folded = fold_stereo_bond_group(id, group)?;
            for delta in &folded {
                if let StereoBondDelta::Add { ligands, .. } = delta {
                    if ligands.iter().any(|l| removed_atoms.contains(&l.atom_id)) {
                        return Err(Contradiction);
                    }
                }
            }
            out.extend(folded.into_iter().map(Delta::StereoBond));
        }
        // Molecule-level constraints are a set difference. Normalize each value, collapse repeated
        // operations, and cancel values present in both the add and remove sets.
        let mut net: BTreeMap<Constraint, (bool, bool)> = BTreeMap::new();
        for delta in constraints {
            match delta {
                ConstraintDelta::Add(constraint) => {
                    net.entry(constraint.normalize()?).or_default().0 = true;
                }
                ConstraintDelta::Remove(constraint) => {
                    net.entry(constraint.normalize()?).or_default().1 = true;
                }
            }
        }
        for (constraint, (add, remove)) in net {
            match (add, remove) {
                (true, false) => out.push(Delta::Constraint(ConstraintDelta::Add(constraint))),
                (false, true) => out.push(Delta::Constraint(ConstraintDelta::Remove(constraint))),
                _ => {}
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

    use super::super::constraint::{MoleculeConstraint, RelationalConstraint};
    use super::super::frame::{
        AromaticSystemsFrameAction, DativeBondsFrameAction, MulticenterBondsFrameAction,
        NoncovalentBondsFrameAction, StereoAtomsFrameAction, StereoBondsFrameAction,
    };
    use super::super::noncovalent::NoncovalentBondKind;
    use super::super::num::NumForm;
    use super::*;
    use crate::ir::{
        AromaticSystemConstraintsForm, AtomConstraintsForm, BondConstraintsForm, BooleanForm,
        DativeBondConstraintsForm, ElectronCountsForm, ElementForm, IsotopeMassForm,
        MulticenterBondConstraintsForm, NoncovalentBondConstraintsForm, NoncovalentBondKindForm,
        RingScope, StereoAtomConstraintsForm, StereoBondConstraintsForm, StereoConfigurationForm,
        StereoConfigurationUpdate, StereoCoset, StereoKind, StereoLigandKind, Stereogenicity,
        StereogenicityForm, UnpairedElectronsForm, UnpairedElectronsUpdate,
    };

    #[fixture]
    fn overlays_frame_action() -> OverlaysFrameAction {
        OverlaysFrameAction::new(
            DativeBondsFrameAction::from_vec(vec![]).expect("actions are admissible"),
            AromaticSystemsFrameAction::from_vec(vec![
                DynPermutation::try_from(vec![2, 0, 1]).expect("action is a permutation")
            ])
            .expect("action is admissible"),
            MulticenterBondsFrameAction::from_vec(vec![]).expect("actions are admissible"),
            NoncovalentBondsFrameAction::from_vec(vec![
                DynPermutation::try_from(vec![1, 0]).expect("action is a permutation")
            ])
            .expect("action is admissible"),
            StereoAtomsFrameAction::from_vec(vec![Permutation::from_image(&[1, 0, 2, 3])])
                .expect("action is admissible"),
            StereoBondsFrameAction::from_vec(vec![Permutation::from_image(&[1, 0, 2, 3])])
                .expect("action is admissible"),
        )
    }

    #[rstest]
    #[case::add_remove(
        AtomDelta::Add { id: AtomId(0), attributes: AtomForm::from_element(Element::C) },
        AtomDelta::Remove { id: AtomId(0), attributes: AtomForm::from_element(Element::C) }
    )]
    #[case::set_field(
        AtomDelta::ModifyField {
            id: AtomId(1),
            change: AtomFieldChange::Charge { old: NumForm::Lit(0), new: NumForm::Lit(1) },
        },
        AtomDelta::ModifyField {
            id: AtomId(1),
            change: AtomFieldChange::Charge { old: NumForm::Lit(1), new: NumForm::Lit(0) },
        }
    )]
    #[case::set_constraint(
        AtomDelta::ModifyConstraint {
            id: AtomId(2),
            old: Some(AtomConstraintForm::Valence(NumForm::Lit(4))),
            new: Some(AtomConstraintForm::Valence(NumForm::Lit(3))),
        },
        AtomDelta::ModifyConstraint {
            id: AtomId(2),
            old: Some(AtomConstraintForm::Valence(NumForm::Lit(3))),
            new: Some(AtomConstraintForm::Valence(NumForm::Lit(4))),
        }
    )]
    fn test_atom_delta_inverse(#[case] input: AtomDelta, #[case] expected: AtomDelta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    fn test_atom_delta_for_update() {
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
        assert_eq!(
            AtomDelta::for_update(AtomId(7), &current, &update),
            vec![
                AtomDelta::ModifyField {
                    id: AtomId(7),
                    change: AtomFieldChange::Element {
                        old: ElementForm::Lit(Element::C),
                        new: ElementForm::Lit(Element::N),
                    },
                },
                AtomDelta::ModifyField {
                    id: AtomId(7),
                    change: AtomFieldChange::IsotopeMass {
                        old: IsotopeMassForm::Lit(12),
                        new: IsotopeMassForm::Lit(13),
                    },
                },
                AtomDelta::ModifyField {
                    id: AtomId(7),
                    change: AtomFieldChange::Charge {
                        old: NumForm::Lit(0),
                        new: NumForm::Lit(1),
                    },
                },
                AtomDelta::ModifyField {
                    id: AtomId(7),
                    change: AtomFieldChange::ImplicitHydrogens {
                        old: NumForm::Lit(4),
                        new: NumForm::Lit(3),
                    },
                },
                AtomDelta::ModifyField {
                    id: AtomId(7),
                    change: AtomFieldChange::LonePairs {
                        old: NumForm::Lit(0),
                        new: NumForm::Lit(1),
                    },
                },
                AtomDelta::ModifyField {
                    id: AtomId(7),
                    change: AtomFieldChange::UnpairedElectrons {
                        old: UnpairedElectronsForm::from((2_u8, 3_u8)),
                        new: UnpairedElectronsForm::from((2_u8, 1_u8)),
                    },
                },
                AtomDelta::ModifyConstraint {
                    id: AtomId(7),
                    old: Some(AtomConstraintForm::valence(4_i64)),
                    new: None,
                },
                AtomDelta::ModifyConstraint {
                    id: AtomId(7),
                    old: None,
                    new: Some(AtomConstraintForm::degree(2_i64)),
                },
            ]
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AtomForm::from_element(Element::C), AtomUpdate::default())]
    #[case::normalized_field(AtomForm::from_element(Element::C).with_charge(1_i64), AtomUpdate { charge: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(AtomForm::from_element(Element::C), AtomUpdate { constraints: AtomConstraintsForm::from(AtomConstraintForm::valence(NumForm::Undetermined)), ..Default::default() })]
    fn test_atom_delta_for_update_identity(#[case] current: AtomForm, #[case] update: AtomUpdate) {
        assert_eq!(AtomDelta::for_update(AtomId(0), &current, &update), Vec::new());
    }

    #[rstest]
    #[case::singleton_set(NumForm::Lit(1), NumForm::lit_set([1]))]
    fn test_atom_delta_diff_canonical(#[case] lhs: NumForm, #[case] rhs: NumForm) {
        // Equivalent charges that are structurally distinct → `diff` emits nothing.
        let lhs = AtomForm::from_element(Element::C).with_charge(lhs);
        let rhs = AtomForm::from_element(Element::C).with_charge(rhs);
        assert_eq!(AtomDelta::diff(AtomId(0), &lhs, &rhs), Vec::new());
    }

    #[rstest]
    #[case::add_remove(
        BondDelta::Add {
            id: BondId(0),
            atoms: [AtomId(0), AtomId(1)],
            attributes: BondForm::default(),
        },
        BondDelta::Remove {
            id: BondId(0),
            atoms: [AtomId(0), AtomId(1)],
            attributes: BondForm::default(),
        }
    )]
    #[case::set_field(
        BondDelta::ModifyField {
            id: BondId(2),
            change: BondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
        },
        BondDelta::ModifyField {
            id: BondId(2),
            change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(1) },
        }
    )]
    #[case::set_constraint(
        BondDelta::ModifyConstraint {
            id: BondId(3),
            old: None,
            new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        },
        BondDelta::ModifyConstraint {
            id: BondId(3),
            old: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            new: None,
        }
    )]
    fn test_bond_delta_inverse(#[case] input: BondDelta, #[case] expected: BondDelta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    fn test_bond_delta_for_update() {
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
        assert_eq!(
            BondDelta::for_update(BondId(7), &current, &update),
            vec![
                BondDelta::ModifyField {
                    id: BondId(7),
                    change: BondFieldChange::Order {
                        old: NumForm::Lit(1),
                        new: NumForm::Lit(2),
                    },
                },
                BondDelta::ModifyField {
                    id: BondId(7),
                    change: BondFieldChange::Charge {
                        old: NumForm::Lit(0),
                        new: NumForm::Undetermined,
                    },
                },
                BondDelta::ModifyField {
                    id: BondId(7),
                    change: BondFieldChange::UnpairedElectrons {
                        old: UnpairedElectronsForm::from((2_u8, 3_u8)),
                        new: UnpairedElectronsForm::from((2_u8, 1_u8)),
                    },
                },
                BondDelta::ModifyConstraint {
                    id: BondId(7),
                    old: None,
                    new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                },
                BondDelta::ModifyConstraint {
                    id: BondId(7),
                    old: Some(BondConstraintForm::ring_membership(
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
    #[case::empty(BondForm::from_order(1), BondUpdate::default())]
    #[case::normalized_field(BondForm::from_order(1).with_charge(1_i64), BondUpdate { charge: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(BondForm::from_order(1), BondUpdate { constraints: BondConstraintsForm::from(BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() })]
    fn test_bond_delta_for_update_identity(#[case] current: BondForm, #[case] update: BondUpdate) {
        assert_eq!(BondDelta::for_update(BondId(0), &current, &update), Vec::new());
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
            DativeBondDelta::ModifyField {
                id: DativeBondId(7),
                change: DativeBondFieldChange::Order { old: NumForm::Lit(1), new: NumForm::Lit(2) },
            },
            DativeBondDelta::ModifyConstraint {
                id: DativeBondId(7),
                old: None,
                new: Some(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
            },
            DativeBondDelta::ModifyConstraint {
                id: DativeBondId(7),
                old: Some(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)),
                new: None,
            },
        ],
    )]
    fn test_dative_bond_delta_for_update(
        #[case] current: DativeBondForm,
        #[case] update: DativeBondUpdate,
        #[case] expected: Vec<DativeBondDelta>,
    ) {
        assert_eq!(DativeBondDelta::for_update(DativeBondId(7), &current, &update), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(DativeBondForm::from_order(1), DativeBondUpdate::default())]
    #[case::normalized_field(DativeBondForm::from_order(1), DativeBondUpdate { order: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(DativeBondForm::from_order(1), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() })]
    fn test_dative_bond_delta_for_update_identity(
        #[case] current: DativeBondForm,
        #[case] update: DativeBondUpdate,
    ) {
        assert_eq!(
            DativeBondDelta::for_update(DativeBondId(0), &current, &update),
            Vec::new(),
        );
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
            AromaticSystemDelta::ModifyField {
                id: AromaticSystemId(7),
                change: AromaticSystemFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 1, 1]), new: ElectronCountsForm::Lit(vec![2, 2, 2]) },
            },
            AromaticSystemDelta::ModifyField {
                id: AromaticSystemId(7),
                change: AromaticSystemFieldChange::Charge { old: NumForm::Lit(0), new: NumForm::Undetermined },
            },
            AromaticSystemDelta::ModifyField {
                id: AromaticSystemId(7),
                change: AromaticSystemFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::from((2_u8, 3_u8)), new: UnpairedElectronsForm::from((2_u8, 1_u8)) },
            },
            AromaticSystemDelta::ModifyConstraint {
                id: AromaticSystemId(7),
                old: Some(AromaticSystemConstraintForm::electron_count(6_i64)),
                new: None,
            },
        ],
    )]
    fn test_aromatic_system_delta_for_update(
        #[case] current: AromaticSystemForm,
        #[case] update: AromaticSystemUpdate,
        #[case] expected: Vec<AromaticSystemDelta>,
    ) {
        assert_eq!(AromaticSystemDelta::for_update(AromaticSystemId(7), &current, &update), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate::default())]
    #[case::normalized_field(AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(1_i64), AromaticSystemUpdate { charge: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(AromaticSystemForm::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() })]
    fn test_aromatic_system_delta_for_update_identity(
        #[case] current: AromaticSystemForm,
        #[case] update: AromaticSystemUpdate,
    ) {
        assert_eq!(
            AromaticSystemDelta::for_update(AromaticSystemId(0), &current, &update),
            Vec::new(),
        );
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
            MulticenterBondDelta::ModifyField {
                id: MulticenterBondId(7),
                change: MulticenterBondFieldChange::Electrons { old: ElectronCountsForm::Lit(vec![1, 1, 1]), new: ElectronCountsForm::Lit(vec![2, 2, 2]) },
            },
            MulticenterBondDelta::ModifyField {
                id: MulticenterBondId(7),
                change: MulticenterBondFieldChange::Charge { old: NumForm::Lit(0), new: NumForm::Undetermined },
            },
            MulticenterBondDelta::ModifyField {
                id: MulticenterBondId(7),
                change: MulticenterBondFieldChange::UnpairedElectrons { old: UnpairedElectronsForm::from((2_u8, 3_u8)), new: UnpairedElectronsForm::from((2_u8, 1_u8)) },
            },
            MulticenterBondDelta::ModifyConstraint {
                id: MulticenterBondId(7),
                old: Some(MulticenterBondConstraintForm::electron_count(6_i64)),
                new: None,
            },
        ],
    )]
    fn test_multicenter_bond_delta_for_update(
        #[case] current: MulticenterBondForm,
        #[case] update: MulticenterBondUpdate,
        #[case] expected: Vec<MulticenterBondDelta>,
    ) {
        assert_eq!(
            MulticenterBondDelta::for_update(MulticenterBondId(7), &current, &update),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate::default())]
    #[case::normalized_field(MulticenterBondForm::from_electrons(vec![1, 1, 1]).with_charge(1_i64), MulticenterBondUpdate { charge: Some(NumForm::lit_set([1])), ..Default::default() })]
    #[case::absent_constraint_removal(MulticenterBondForm::from_electrons(vec![1, 1, 1]), MulticenterBondUpdate { constraints: MulticenterBondConstraintsForm::from(MulticenterBondConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() })]
    fn test_multicenter_bond_delta_for_update_identity(
        #[case] current: MulticenterBondForm,
        #[case] update: MulticenterBondUpdate,
    ) {
        assert_eq!(
            MulticenterBondDelta::for_update(MulticenterBondId(0), &current, &update),
            Vec::new(),
        );
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
            NoncovalentBondDelta::ModifyField {
                id: NoncovalentBondId(7),
                change: NoncovalentBondFieldChange::Kind { old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond), new: NoncovalentBondKindForm::Undetermined },
            },
            NoncovalentBondDelta::ModifyConstraint {
                id: NoncovalentBondId(7),
                old: Some(NoncovalentBondConstraintForm::intramolecular(true)),
                new: None,
            },
        ],
    )]
    fn test_noncovalent_bond_delta_for_update(
        #[case] current: NoncovalentBondForm,
        #[case] update: NoncovalentBondUpdate,
        #[case] expected: Vec<NoncovalentBondDelta>,
    ) {
        assert_eq!(
            NoncovalentBondDelta::for_update(NoncovalentBondId(7), &current, &update),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate::default())]
    #[case::same_kind(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate { kind: Some(NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond)), ..Default::default() })]
    #[case::absent_constraint_removal(NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond), NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsForm::from(NoncovalentBondConstraintForm::intramolecular(BooleanForm::Undetermined)), ..Default::default() })]
    fn test_noncovalent_bond_delta_for_update_identity(
        #[case] current: NoncovalentBondForm,
        #[case] update: NoncovalentBondUpdate,
    ) {
        assert_eq!(
            NoncovalentBondDelta::for_update(NoncovalentBondId(0), &current, &update),
            Vec::new(),
        );
    }

    #[rstest]
    #[case::add_remove(
        StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        },
        StereoAtomDelta::Remove {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        }
    )]
    #[case::set_field(
        StereoAtomDelta::ModifyField {
            id: StereoAtomId(1),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                new: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            },
        },
        StereoAtomDelta::ModifyField {
            id: StereoAtomId(1),
            change: StereoAtomFieldChange::Configuration {
                old: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
                new: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            },
        }
    )]
    #[case::set_constraint(
        StereoAtomDelta::ModifyConstraint {
            id: StereoAtomId(2),
            kind: Some(StereoKind::Tetrahedral),
            old: None,
            new: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
        },
        StereoAtomDelta::ModifyConstraint {
            id: StereoAtomId(2),
            kind: Some(StereoKind::Tetrahedral),
            old: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)),
            new: None,
        }
    )]
    fn test_stereo_atom_delta_inverse(
        #[case] input: StereoAtomDelta,
        #[case] expected: StereoAtomDelta,
    ) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
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
            StereoAtomDelta::ModifyField {
                id: StereoAtomId(7),
                change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32) },
            },
            StereoAtomDelta::ModifyConstraint {
                id: StereoAtomId(7),
                kind: Some(StereoKind::Tetrahedral),
                old: Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
                new: None,
            },
        ],
    )]
    fn test_stereo_atom_delta_for_update(
        #[case] current: StereoAtomForm,
        #[case] update: StereoAtomUpdate,
        #[case] expected: Vec<StereoAtomDelta>,
    ) {
        assert_eq!(
            StereoAtomDelta::for_update(StereoAtomId(7), &current, &update),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate::default())]
    #[case::relative(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() })]
    #[case::absent_constraint_removal(StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32), StereoAtomUpdate { constraints: StereoAtomConstraintsForm::from(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() })]
    fn test_stereo_atom_delta_for_update_identity(
        #[case] current: StereoAtomForm,
        #[case] update: StereoAtomUpdate,
    ) {
        assert_eq!(
            StereoAtomDelta::for_update(StereoAtomId(0), &current, &update),
            Vec::new(),
        );
    }

    #[rstest]
    #[case::add_remove(
        StereoBondDelta::Add {
            id: StereoBondId(0),
            site: BondId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            attributes: StereoBondForm::new(StereoKind::CisTrans, 0u32),
        },
        StereoBondDelta::Remove {
            id: StereoBondId(0),
            site: BondId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            attributes: StereoBondForm::new(StereoKind::CisTrans, 0u32),
        }
    )]
    fn test_stereo_bond_delta_inverse(
        #[case] input: StereoBondDelta,
        #[case] expected: StereoBondDelta,
    ) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
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
            StereoBondDelta::ModifyField {
                id: StereoBondId(7),
                change: StereoBondFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32) },
            },
            StereoBondDelta::ModifyConstraint {
                id: StereoBondId(7),
                kind: Some(StereoKind::CisTrans),
                old: Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic))),
                new: None,
            },
        ],
    )]
    fn test_stereo_bond_delta_for_update(
        #[case] current: StereoBondForm,
        #[case] update: StereoBondUpdate,
        #[case] expected: Vec<StereoBondDelta>,
    ) {
        assert_eq!(
            StereoBondDelta::for_update(StereoBondId(7), &current, &update),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(StereoBondForm::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate::default())]
    #[case::relative(StereoBondForm::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() })]
    #[case::absent_constraint_removal(StereoBondForm::new(StereoKind::CisTrans, 1_u32), StereoBondUpdate { constraints: StereoBondConstraintsForm::from(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined)), ..Default::default() })]
    fn test_stereo_bond_delta_for_update_identity(
        #[case] current: StereoBondForm,
        #[case] update: StereoBondUpdate,
    ) {
        assert_eq!(
            StereoBondDelta::for_update(StereoBondId(0), &current, &update),
            Vec::new(),
        );
    }

    #[rstest]
    fn test_stereo_atom_delta_diff() {
        assert_eq!(
            StereoAtomDelta::diff(
                StereoAtomId(0),
                &StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
                &StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            ),
            vec![StereoAtomDelta::ModifyField {
                id: StereoAtomId(0),
                change: StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0)
                    ),
                    new: StereoConfigurationForm::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1)
                    ),
                },
            }],
        );
    }

    #[rstest]
    fn test_stereo_atom_delta_apply_field() {
        let mut attributes = StereoAtomForm::new(StereoKind::Tetrahedral, 0u32);
        StereoAtomDelta::apply_field(
            &mut attributes,
            StereoAtomFieldChange::Configuration {
                old: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                new: StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            },
        )
        .unwrap();
        assert_eq!(
            attributes,
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32)
        );
    }

    #[rstest]
    fn test_stereo_atom_delta_apply_field_error() {
        let mut attributes = StereoAtomForm::new(StereoKind::Tetrahedral, 1u32);
        assert_eq!(
            StereoAtomDelta::apply_field(
                &mut attributes,
                StereoAtomFieldChange::Configuration {
                    old: StereoConfigurationForm::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(0)
                    ),
                    new: StereoConfigurationForm::Kinded(
                        StereoKind::Tetrahedral,
                        StereoCoset::Lit(1)
                    ),
                },
            ),
            Err(Contradiction),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::configuration_chain(
        vec![
            StereoAtomDelta::ModifyField { id: StereoAtomId(0), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32) } },
            StereoAtomDelta::ModifyField { id: StereoAtomId(0), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32), new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32) } },
        ],
        vec![],
    )]
    #[case::addition(
        vec![
            StereoAtomDelta::Add {
                id: StereoAtomId(0),
                site: AtomId(0),
                ligands: vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            },
            StereoAtomDelta::ModifyField { id: StereoAtomId(0), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32) } },
        ],
        vec![StereoAtomDelta::Add {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        }],
    )]
    #[case::removal(
        vec![
            StereoAtomDelta::ModifyField { id: StereoAtomId(0), change: StereoAtomFieldChange::Configuration { old: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0_u32), new: StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1_u32) } },
            StereoAtomDelta::Remove {
                id: StereoAtomId(0),
                site: AtomId(0),
                ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
                attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            },
        ],
        vec![StereoAtomDelta::Remove {
            id: StereoAtomId(0),
            site: AtomId(0),
            ligands: vec![StereoLigand::new(AtomId(1), StereoLigandKind::Atom)],
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        }],
    )]
    fn test_deltas_normalize_stereo_atom(
        #[case] input: Vec<StereoAtomDelta>,
        #[case] expected: Vec<StereoAtomDelta>,
    ) {
        let canon = Deltas::from_iter(input.into_iter().map(Delta::StereoAtom))
            .normalize()
            .unwrap();
        assert_eq!(
            canon,
            Deltas::from_iter(expected.into_iter().map(Delta::StereoAtom)),
        );
    }

    #[rstest]
    fn test_constraint_delta_inverse() {
        let constraint = Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: NumForm::Lit(0),
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
        Delta::Atom(AtomDelta::Add { id: AtomId(0), attributes: AtomForm::from_element(Element::C) }),
        Delta::Atom(AtomDelta::Remove { id: AtomId(0), attributes: AtomForm::from_element(Element::C) })
    )]
    #[case::bond(
        Delta::Bond(BondDelta::Add {
            id: BondId(0),
            atoms: [AtomId(0), AtomId(1)],
            attributes: BondForm::default(),
        }),
        Delta::Bond(BondDelta::Remove {
            id: BondId(0),
            atoms: [AtomId(0), AtomId(1)],
            attributes: BondForm::default(),
        })
    )]
    fn test_delta_inverse(#[case] input: Delta, #[case] expected: Delta) {
        assert_eq!(input.clone().inverse(), expected);
        assert_eq!(input.clone().inverse().inverse(), input);
    }

    #[rstest]
    fn test_dative_bond_delta_reframe_by() {
        let action =
            DynPermutation::try_from([2, 0, 1].as_slice()).expect("test image is a permutation");
        let input = DativeBondDelta::Add {
            id: DativeBondId(4),
            donors: vec![AtomId(1), AtomId(2), AtomId(3)],
            acceptor: AtomId(9),
            attributes: DativeBondForm::from_order(2),
        };
        let expected = DativeBondDelta::Add {
            id: DativeBondId(4),
            donors: vec![AtomId(3), AtomId(1), AtomId(2)],
            acceptor: AtomId(9),
            attributes: DativeBondForm::from_order(2),
        };

        assert_eq!(input.reframe_by(&action), Some(expected));
    }

    #[rstest]
    fn test_aromatic_system_delta_reframe_by() {
        let action =
            DynPermutation::try_from([2, 0, 1].as_slice()).expect("test image is a permutation");
        let input = AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(3),
            change: AromaticSystemFieldChange::Electrons {
                old: ElectronCountsForm::Lit(vec![10, 20, 30]),
                new: ElectronCountsForm::Lit(vec![11, 21, 31]),
            },
        };
        let expected = AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(3),
            change: AromaticSystemFieldChange::Electrons {
                old: ElectronCountsForm::Lit(vec![30, 10, 20]),
                new: ElectronCountsForm::Lit(vec![31, 11, 21]),
            },
        };

        assert_eq!(input.reframe_by(&action), Some(expected));
    }

    #[rstest]
    fn test_multicenter_bond_delta_reframe_by() {
        let action =
            DynPermutation::try_from([2, 0, 1].as_slice()).expect("test image is a permutation");
        let input = MulticenterBondDelta::Remove {
            id: MulticenterBondId(2),
            atoms: vec![AtomId(1), AtomId(2), AtomId(3)],
            attributes: MulticenterBondForm::from_electrons(vec![2, 4, 6]),
        };
        let expected = MulticenterBondDelta::Remove {
            id: MulticenterBondId(2),
            atoms: vec![AtomId(3), AtomId(1), AtomId(2)],
            attributes: MulticenterBondForm::from_electrons(vec![6, 2, 4]),
        };

        assert_eq!(input.reframe_by(&action), Some(expected));
    }

    #[rstest]
    fn test_noncovalent_bond_delta_reframe_by() {
        let action =
            DynPermutation::try_from([1, 0].as_slice()).expect("test image is a permutation");
        let input = NoncovalentBondDelta::Add {
            id: NoncovalentBondId(5),
            atoms: [AtomId(1), AtomId(2)],
            attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        };
        let expected = NoncovalentBondDelta::Add {
            id: NoncovalentBondId(5),
            atoms: [AtomId(2), AtomId(1)],
            attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        };

        assert_eq!(input.reframe_by(&action), Some(expected));
    }

    #[rstest]
    fn test_noncovalent_bond_delta_reframe_by_error() {
        let action = DynPermutation::identity(3);
        let input = NoncovalentBondDelta::ModifyField {
            id: NoncovalentBondId(5),
            change: NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond),
                new: NoncovalentBondKindForm::Lit(NoncovalentBondKind::HalogenBond),
            },
        };

        assert_eq!(input.reframe_by(&action), None);
    }

    #[rstest]
    fn test_stereo_atom_delta_reframe_by() {
        let action = Permutation::from_image(&[1, 0, 2, 3]);
        let input = StereoAtomDelta::Add {
            id: StereoAtomId(6),
            site: AtomId(0),
            ligands: [1, 2, 3, 4]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        };
        let expected = StereoAtomDelta::Add {
            id: StereoAtomId(6),
            site: AtomId(0),
            ligands: [2, 1, 3, 4]
                .into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 1_u32),
        };

        assert_eq!(input.reframe_by(&action), Some(expected));
    }

    #[rstest]
    fn test_stereo_bond_delta_reframe_by() {
        let action = Permutation::from_image(&[1, 0, 2, 3]);
        let input = StereoBondDelta::ModifyField {
            id: StereoBondId(7),
            change: StereoBondFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32),
                new: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32),
            },
        };
        let expected = StereoBondDelta::ModifyField {
            id: StereoBondId(7),
            change: StereoBondFieldChange::Configuration {
                old: StereoConfigurationForm::kinded(StereoKind::CisTrans, 1_u32),
                new: StereoConfigurationForm::kinded(StereoKind::CisTrans, 0_u32),
            },
        };

        assert_eq!(input.reframe_by(&action), Some(expected));
    }

    #[rstest]
    fn test_stereo_bond_delta_reframe_by_error() {
        let action = Permutation::from_image(&[1, 2, 0, 3]);
        let input = StereoBondDelta::ModifyConstraint {
            id: StereoBondId(7),
            kind: None,
            old: None,
            new: None,
        };

        assert_eq!(input.reframe_by(&action), None);
    }

    #[rstest]
    fn test_entity_span_reframe_by() {
        let action =
            DynPermutation::try_from([2, 0, 1].as_slice()).expect("test image is a permutation");
        let input = EntitySpan::Modified {
            lhs: AromaticSystemForm::from_electrons(vec![10, 20, 30]),
            rhs: AromaticSystemForm::from_electrons(vec![11, 21, 31]),
        };
        let expected = EntitySpan::Modified {
            lhs: AromaticSystemForm::from_electrons(vec![30, 10, 20]),
            rhs: AromaticSystemForm::from_electrons(vec![31, 11, 21]),
        };

        assert_eq!(input.reframe_by(&action), Some(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unchanged(
        EntitySpan::Unchanged(NumForm::lit_set([2])),
        EntitySpan::Unchanged(NumForm::Lit(2)),
    )]
    #[case::modified_collapses(
        EntitySpan::Modified { lhs: NumForm::lit_set([2]), rhs: NumForm::Lit(2) },
        EntitySpan::Unchanged(NumForm::Lit(2)),
    )]
    #[case::modified_remains(
        EntitySpan::Modified { lhs: NumForm::lit_set([2]), rhs: NumForm::Lit(3) },
        EntitySpan::Modified { lhs: NumForm::Lit(2), rhs: NumForm::Lit(3) },
    )]
    fn test_entity_span_normalize(
        #[case] input: EntitySpan<NumForm>,
        #[case] expected: EntitySpan<NumForm>,
    ) {
        assert_eq!(input.normalize(), Ok(expected));
    }

    #[rstest]
    fn test_constraint_span_normalize() {
        let input = ConstraintSpan::Added(Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::Valence(NumForm::lit_set([4])),
        ));
        let expected = ConstraintSpan::Added(Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::Valence(NumForm::Lit(4)),
        ));

        assert_eq!(input.normalize(), Ok(expected));
    }

    #[rstest]
    fn test_constraint_delta_reframe_by(overlays_frame_action: OverlaysFrameAction) {
        let input = ConstraintDelta::Add(Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                ],
            },
        ));
        let expected = ConstraintDelta::Add(Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                ],
            },
        ));

        assert_eq!(input.reframe_by(&overlays_frame_action), Some(expected));
    }

    #[rstest]
    fn test_constraint_span_reframe_by(overlays_frame_action: OverlaysFrameAction) {
        let input = ConstraintSpan::Removed(Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                ],
            },
        ));
        let expected = ConstraintSpan::Removed(Constraint::Relational(
            RelationalConstraint::NoncovalentBondEndsSatisfy {
                bond: NoncovalentBondId(0),
                predicates: [
                    Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
                    Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                ],
            },
        ));

        assert_eq!(input.reframe_by(&overlays_frame_action), Some(expected));
    }

    #[rstest]
    fn test_delta_reframe_by(overlays_frame_action: OverlaysFrameAction) {
        let input = Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(0),
            change: AromaticSystemFieldChange::Electrons {
                old: ElectronCountsForm::Lit(vec![10, 20, 30]),
                new: ElectronCountsForm::Lit(vec![11, 21, 31]),
            },
        });
        let expected = Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(0),
            change: AromaticSystemFieldChange::Electrons {
                old: ElectronCountsForm::Lit(vec![30, 10, 20]),
                new: ElectronCountsForm::Lit(vec![31, 11, 21]),
            },
        });

        assert_eq!(input.reframe_by(&overlays_frame_action), Some(expected));
    }

    #[rstest]
    fn test_delta_reframe_by_error(overlays_frame_action: OverlaysFrameAction) {
        let input = Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(1),
            change: AromaticSystemFieldChange::Electrons {
                old: ElectronCountsForm::Lit(vec![10, 20, 30]),
                new: ElectronCountsForm::Lit(vec![11, 21, 31]),
            },
        });

        assert_eq!(input.reframe_by(&overlays_frame_action), None);
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(5), Some(&5))]
    #[case::modified(EntitySpan::Modified { lhs: 1, rhs: 2 }, Some(&1))]
    #[case::removed(EntitySpan::Removed(7), Some(&7))]
    #[case::added(EntitySpan::Added(9), None)]
    fn test_entity_span_lhs(#[case] state: EntitySpan<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(state.lhs(), expected);
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(5), Some(&5))]
    #[case::modified(EntitySpan::Modified { lhs: 1, rhs: 2 }, Some(&2))]
    #[case::added(EntitySpan::Added(9), Some(&9))]
    #[case::removed(EntitySpan::Removed(7), None)]
    fn test_entity_span_rhs(#[case] state: EntitySpan<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(state.rhs(), expected);
    }

    #[rstest]
    #[case::unchanged(EntitySpan::Unchanged(5), Some(EntitySpan::Unchanged(10)))]
    #[case::modified(
        EntitySpan::Modified { lhs: 1, rhs: 2 },
        Some(EntitySpan::Modified { lhs: 2, rhs: 4 }),
    )]
    #[case::added(EntitySpan::Added(9), Some(EntitySpan::Added(18)))]
    #[case::removed(EntitySpan::Removed(7), Some(EntitySpan::Removed(14)))]
    #[case::modified_lhs_declines(EntitySpan::Modified { lhs: 0, rhs: 2 }, None)]
    #[case::modified_rhs_declines(EntitySpan::Modified { lhs: 1, rhs: 0 }, None)]
    #[case::unchanged_declines(EntitySpan::Unchanged(0), None)]
    fn test_entity_span_try_map(
        #[case] span: EntitySpan<i32>,
        #[case] expected: Option<EntitySpan<i32>>,
    ) {
        assert_eq!(span.try_map(|v| (v != 0).then_some(v * 2)), expected);
    }

    #[rstest]
    #[case::singleton_set(NumForm::Lit(1), NumForm::lit_set([1]))]
    fn test_entity_span_superimpose_canonical(#[case] lhs: NumForm, #[case] rhs: NumForm) {
        // Equivalent sides that are structurally distinct superimpose to `Unchanged`,
        // not `Modified`.
        let lhs = AtomForm::from_element(Element::C).with_charge(lhs);
        let rhs = AtomForm::from_element(Element::C).with_charge(rhs);
        assert_eq!(
            EntitySpan::superimpose(Some(lhs.clone()), Some(rhs)),
            Some(EntitySpan::Unchanged(lhs))
        );
    }

    #[fixture]
    fn remapping() -> IdRemapping {
        IdRemapping::new(
            HashMap::from([
                (AtomId(0), AtomId(2)),
                (AtomId(1), AtomId(0)),
                (AtomId(2), AtomId(1)),
                (AtomId(3), AtomId(4)),
                (AtomId(4), AtomId(3)),
            ]),
            HashMap::from([(BondId(0), BondId(1)), (BondId(1), BondId(0))]),
            HashMap::from([(DativeBondId(0), DativeBondId(1))]),
            HashMap::from([(AromaticSystemId(0), AromaticSystemId(1))]),
            HashMap::from([(MulticenterBondId(0), MulticenterBondId(1))]),
            HashMap::from([(NoncovalentBondId(0), NoncovalentBondId(1))]),
            HashMap::from([(StereoAtomId(0), StereoAtomId(1))]),
            HashMap::from([(StereoBondId(0), StereoBondId(1))]),
        )
    }

    #[fixture]
    fn inverse_remapping() -> IdRemapping {
        IdRemapping::new(
            HashMap::from([
                (AtomId(0), AtomId(1)),
                (AtomId(1), AtomId(2)),
                (AtomId(2), AtomId(0)),
                (AtomId(3), AtomId(4)),
                (AtomId(4), AtomId(3)),
            ]),
            HashMap::from([(BondId(0), BondId(1)), (BondId(1), BondId(0))]),
            HashMap::from([(DativeBondId(1), DativeBondId(0))]),
            HashMap::from([(AromaticSystemId(1), AromaticSystemId(0))]),
            HashMap::from([(MulticenterBondId(1), MulticenterBondId(0))]),
            HashMap::from([(NoncovalentBondId(1), NoncovalentBondId(0))]),
            HashMap::from([(StereoAtomId(1), StereoAtomId(0))]),
            HashMap::from([(StereoBondId(1), StereoBondId(0))]),
        )
    }

    #[fixture]
    fn source_frame_action() -> OverlaysFrameAction {
        OverlaysFrameAction::new(
            DativeBondsFrameAction::from_vec(vec![
                DynPermutation::try_from(vec![2, 0, 1]).expect("action is a permutation")
            ])
            .expect("action is admissible"),
            AromaticSystemsFrameAction::from_vec(vec![
                DynPermutation::try_from(vec![2, 0, 1]).expect("action is a permutation")
            ])
            .expect("action is admissible"),
            MulticenterBondsFrameAction::from_vec(vec![
                DynPermutation::try_from(vec![2, 0, 1]).expect("action is a permutation")
            ])
            .expect("action is admissible"),
            NoncovalentBondsFrameAction::from_vec(vec![
                DynPermutation::try_from(vec![1, 0]).expect("action is a permutation")
            ])
            .expect("action is admissible"),
            StereoAtomsFrameAction::from_vec(vec![Permutation::from_image(&[1, 0, 2, 3])])
                .expect("action is admissible"),
            StereoBondsFrameAction::from_vec(vec![Permutation::from_image(&[1, 0, 2, 3])])
                .expect("action is admissible"),
        )
    }

    #[fixture]
    fn target_frame_action() -> OverlaysFrameAction {
        OverlaysFrameAction::new(
            DativeBondsFrameAction::from_vec(vec![
                DynPermutation::identity(3),
                DynPermutation::try_from(vec![2, 0, 1]).expect("action is a permutation"),
            ])
            .expect("actions are admissible"),
            AromaticSystemsFrameAction::from_vec(vec![
                DynPermutation::identity(3),
                DynPermutation::try_from(vec![2, 0, 1]).expect("action is a permutation"),
            ])
            .expect("actions are admissible"),
            MulticenterBondsFrameAction::from_vec(vec![
                DynPermutation::identity(3),
                DynPermutation::try_from(vec![2, 0, 1]).expect("action is a permutation"),
            ])
            .expect("actions are admissible"),
            NoncovalentBondsFrameAction::from_vec(vec![
                DynPermutation::identity(2),
                DynPermutation::try_from(vec![1, 0]).expect("action is a permutation"),
            ])
            .expect("actions are admissible"),
            StereoAtomsFrameAction::from_vec(vec![
                Permutation::identity(4),
                Permutation::from_image(&[1, 0, 2, 3]),
            ])
            .expect("actions are admissible"),
            StereoBondsFrameAction::from_vec(vec![
                Permutation::identity(4),
                Permutation::from_image(&[1, 0, 2, 3]),
            ])
            .expect("actions are admissible"),
        )
    }

    #[rstest]
    #[case::atom(
        Delta::Atom(AtomDelta::Add { id: AtomId(1), attributes: AtomForm::from_element(Element::C) }),
        Delta::Atom(AtomDelta::Add { id: AtomId(0), attributes: AtomForm::from_element(Element::C) })
    )]
    #[case::bond(
        Delta::Bond(BondDelta::Add {
            id: BondId(0),
            atoms: [AtomId(2), AtomId(1)],
            attributes: BondForm::default(),
        }),
        Delta::Bond(BondDelta::Add {
            id: BondId(1),
            atoms: [AtomId(1), AtomId(0)],
            attributes: BondForm::default(),
        })
    )]
    #[case::dative_bond(
        Delta::DativeBond(DativeBondDelta::Add {
            id: DativeBondId(0),
            donors: vec![AtomId(0), AtomId(2)],
            acceptor: AtomId(1),
            attributes: DativeBondForm::from_order(1),
        }),
        Delta::DativeBond(DativeBondDelta::Add {
            id: DativeBondId(1),
            donors: vec![AtomId(2), AtomId(1)],
            acceptor: AtomId(0),
            attributes: DativeBondForm::from_order(1),
        })
    )]
    #[case::aromatic_system_add(
        Delta::AromaticSystem(AromaticSystemDelta::Add {
            id: AromaticSystemId(0),
            atoms: vec![AtomId(0), AtomId(1)],
            attributes: AromaticSystemForm::from_electrons(vec![1, 2]),
        }),
        Delta::AromaticSystem(AromaticSystemDelta::Add {
            id: AromaticSystemId(1),
            atoms: vec![AtomId(2), AtomId(0)],
            attributes: AromaticSystemForm::from_electrons(vec![1, 2]),
        })
    )]
    #[case::aromatic_remove(
        Delta::AromaticSystem(AromaticSystemDelta::Remove {
            id: AromaticSystemId(0),
            atoms: vec![AtomId(0), AtomId(1)],
            attributes: AromaticSystemForm::from_electrons(vec![1, 2]),
        }),
        Delta::AromaticSystem(AromaticSystemDelta::Remove {
            id: AromaticSystemId(1),
            atoms: vec![AtomId(2), AtomId(0)],
            attributes: AromaticSystemForm::from_electrons(vec![1, 2]),
        })
    )]
    #[case::multicenter_bond(
        Delta::MulticenterBond(MulticenterBondDelta::Add {
            id: MulticenterBondId(0),
            atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
            attributes: MulticenterBondForm::from_electrons(vec![3, 5, 7]),
        }),
        Delta::MulticenterBond(MulticenterBondDelta::Add {
            id: MulticenterBondId(1),
            atoms: vec![AtomId(2), AtomId(0), AtomId(1)],
            attributes: MulticenterBondForm::from_electrons(vec![3, 5, 7]),
        })
    )]
    #[case::noncovalent_bond(
        Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(0),
            atoms: [AtomId(2), AtomId(1)],
            attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        }),
        Delta::NoncovalentBond(NoncovalentBondDelta::Add {
            id: NoncovalentBondId(1),
            atoms: [AtomId(1), AtomId(0)],
            attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        })
    )]
    #[case::overlay_modify_field(
        Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(0),
            change: AromaticSystemFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(1),
            },
        }),
        Delta::AromaticSystem(AromaticSystemDelta::ModifyField {
            id: AromaticSystemId(1),
            change: AromaticSystemFieldChange::Charge {
                old: NumForm::Lit(0),
                new: NumForm::Lit(1),
            },
        })
    )]
    #[case::constraint_add(
        Delta::Constraint(ConstraintDelta::Add(Constraint::Relational(
            RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(0),
                parallel: BondId(0),
            },
        ))),
        Delta::Constraint(ConstraintDelta::Add(Constraint::Relational(
            RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(1),
                parallel: BondId(1),
            },
        )))
    )]
    #[case::constraint_remove(
        Delta::Constraint(ConstraintDelta::Remove(Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::valence(3_i64),
        ))),
        Delta::Constraint(ConstraintDelta::Remove(Constraint::Atom(
            AtomId(2),
            AtomConstraintForm::valence(3_i64),
        )))
    )]
    fn test_remap_delta(remapping: IdRemapping, #[case] input: Delta, #[case] expected: Delta) {
        assert_eq!(remap_delta(input, &remapping), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_bond(Delta::DativeBond(DativeBondDelta::Add {
        id: DativeBondId(0),
        donors: vec![AtomId(0), AtomId(1), AtomId(2)],
        acceptor: AtomId(3),
        attributes: DativeBondForm::from_order(2),
    }))]
    #[case::aromatic_system(Delta::AromaticSystem(AromaticSystemDelta::Add {
        id: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
        attributes: AromaticSystemForm::from_electrons(vec![2, 4, 6]),
    }))]
    #[case::multicenter_bond(Delta::MulticenterBond(MulticenterBondDelta::Remove {
        id: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
        attributes: MulticenterBondForm::from_electrons(vec![1, 3, 5]),
    }))]
    #[case::noncovalent_bond(Delta::NoncovalentBond(NoncovalentBondDelta::Add {
        id: NoncovalentBondId(0),
        atoms: [AtomId(0), AtomId(1)],
        attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
    }))]
    #[case::stereo_atom(Delta::StereoAtom(StereoAtomDelta::Add {
        id: StereoAtomId(0),
        site: AtomId(0),
        ligands: [1, 2, 3, 4].into_iter().map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(),
        attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
    }))]
    #[case::stereo_bond(Delta::StereoBond(StereoBondDelta::Add {
        id: StereoBondId(0),
        site: BondId(0),
        ligands: [0, 1, 2, 3].into_iter().map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(),
        attributes: StereoBondForm::new(StereoKind::CisTrans, 0_u32),
    }))]
    #[case::constraint(Delta::Constraint(ConstraintDelta::Add(Constraint::Relational(
        RelationalConstraint::NoncovalentBondEndsSatisfy {
            bond: NoncovalentBondId(0),
            predicates: [
                Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
            ],
        },
    ))))]
    fn test_remap_delta_roundtrip(
        remapping: IdRemapping,
        inverse_remapping: IdRemapping,
        #[case] input: Delta,
    ) {
        let remapped = remap_delta(input.clone(), &remapping);

        assert_eq!(remap_delta(remapped, &inverse_remapping), input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_bond(Delta::DativeBond(DativeBondDelta::Add {
        id: DativeBondId(0),
        donors: vec![AtomId(0), AtomId(1), AtomId(2)],
        acceptor: AtomId(3),
        attributes: DativeBondForm::from_order(2),
    }))]
    #[case::aromatic_system(Delta::AromaticSystem(AromaticSystemDelta::Add {
        id: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
        attributes: AromaticSystemForm::from_electrons(vec![2, 4, 6]),
    }))]
    #[case::multicenter_bond(Delta::MulticenterBond(MulticenterBondDelta::Remove {
        id: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1), AtomId(2)],
        attributes: MulticenterBondForm::from_electrons(vec![1, 3, 5]),
    }))]
    #[case::noncovalent_bond(Delta::NoncovalentBond(NoncovalentBondDelta::Add {
        id: NoncovalentBondId(0),
        atoms: [AtomId(0), AtomId(1)],
        attributes: NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
    }))]
    #[case::stereo_atom(Delta::StereoAtom(StereoAtomDelta::Add {
        id: StereoAtomId(0),
        site: AtomId(0),
        ligands: [1, 2, 3, 4].into_iter().map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(),
        attributes: StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
    }))]
    #[case::stereo_bond(Delta::StereoBond(StereoBondDelta::Add {
        id: StereoBondId(0),
        site: BondId(0),
        ligands: [0, 1, 2, 3].into_iter().map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(),
        attributes: StereoBondForm::new(StereoKind::CisTrans, 0_u32),
    }))]
    #[case::constraint(Delta::Constraint(ConstraintDelta::Add(Constraint::Relational(
        RelationalConstraint::NoncovalentBondEndsSatisfy {
            bond: NoncovalentBondId(0),
            predicates: [
                Box::new(AtomConstraintForm::Valence(NumForm::Lit(4))),
                Box::new(AtomConstraintForm::Degree(NumForm::Lit(2))),
            ],
        },
    ))))]
    fn test_remap_delta_frame_transport(
        remapping: IdRemapping,
        source_frame_action: OverlaysFrameAction,
        target_frame_action: OverlaysFrameAction,
        #[case] input: Delta,
    ) {
        let reframed = input
            .clone()
            .reframe_by(&source_frame_action)
            .expect("source action covers the delta");
        let remapped = remap_delta(input, &remapping);

        assert_eq!(
            remap_delta(reframed, &remapping),
            remapped
                .reframe_by(&target_frame_action)
                .expect("target action covers the remapped delta"),
        );
    }

    fn charge_set(id: u32, old: i64, new: i64) -> Delta {
        Delta::Atom(AtomDelta::ModifyField {
            id: AtomId(id),
            change: AtomFieldChange::Charge {
                old: NumForm::Lit(old),
                new: NumForm::Lit(new),
            },
        })
    }

    #[rstest]
    fn test_deltas_into_iter() {
        let entries = vec![
            Delta::Atom(AtomDelta::ModifyField {
                id: AtomId(0),
                change: AtomFieldChange::Charge {
                    old: NumForm::Lit(0),
                    new: NumForm::Lit(1),
                },
            }),
            Delta::Bond(BondDelta::ModifyField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: NumForm::Lit(1),
                    new: NumForm::Lit(2),
                },
            }),
        ];
        let deltas = Deltas::from_iter(entries.clone());

        assert_eq!(deltas.into_iter().collect::<Vec<_>>(), entries);
    }

    #[rstest]
    fn test_deltas_normalize_field_fusion() {
        let deltas = Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 1, 2)]);
        assert_eq!(
            deltas.normalize().unwrap(),
            Deltas::from_iter([charge_set(0, 0, 2)]),
        );
    }

    #[rstest]
    fn test_deltas_normalize_field_noop_dropped() {
        let deltas = Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 1, 0)]);
        assert_eq!(deltas.normalize().unwrap(), Deltas::new());
    }

    #[rstest]
    fn test_deltas_normalize_created_absorbs_field() {
        let deltas = Deltas::from_iter([
            Delta::Atom(AtomDelta::Add {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C).with_charge(NumForm::Lit(0)),
            }),
            charge_set(0, 0, 1),
        ]);
        assert_eq!(
            deltas.normalize().unwrap(),
            Deltas::from_iter([Delta::Atom(AtomDelta::Add {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C).with_charge(NumForm::Lit(1)),
            })]),
        );
    }

    #[rstest]
    fn test_deltas_normalize_created_then_removed_cancels() {
        let deltas = Deltas::from_iter([
            Delta::Atom(AtomDelta::Add {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C),
            }),
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C),
            }),
        ]);
        assert_eq!(deltas.normalize().unwrap(), Deltas::new());
    }

    #[rstest]
    fn test_deltas_normalize_remove_subsumes_field() {
        // ModifyField then Remove must normalize to a Remove carrying the original value.
        let deltas = Deltas::from_iter([
            charge_set(0, 0, 1),
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C).with_charge(NumForm::Lit(1)),
            }),
        ]);
        assert_eq!(
            deltas.normalize().unwrap(),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C).with_charge(NumForm::Lit(0)),
            })]),
        );
    }

    #[rstest]
    fn test_deltas_normalize_constraint_chain() {
        let deltas = Deltas::from_iter([
            Delta::Atom(AtomDelta::ModifyConstraint {
                id: AtomId(0),
                old: None,
                new: Some(AtomConstraintForm::Valence(NumForm::Lit(4))),
            }),
            Delta::Atom(AtomDelta::ModifyConstraint {
                id: AtomId(0),
                old: Some(AtomConstraintForm::Valence(NumForm::Lit(4))),
                new: Some(AtomConstraintForm::Valence(NumForm::Lit(3))),
            }),
        ]);
        assert_eq!(
            deltas.normalize().unwrap(),
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyConstraint {
                id: AtomId(0),
                old: None,
                new: Some(AtomConstraintForm::Valence(NumForm::Lit(3))),
            })]),
        );
    }

    #[rstest]
    fn test_deltas_normalize_order_independent() {
        let order_set = Delta::Bond(BondDelta::ModifyField {
            id: BondId(0),
            change: BondFieldChange::Order {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        });
        let forward = Deltas::from_iter([charge_set(0, 0, 1), order_set.clone()]);
        let reverse = Deltas::from_iter([order_set, charge_set(0, 0, 1)]);
        assert_eq!(forward.normalize().unwrap(), reverse.normalize().unwrap());
    }

    #[rstest]
    fn test_deltas_normalize_idempotent() {
        let once = Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 1, 2)])
            .normalize()
            .unwrap();
        assert_eq!(once.clone().normalize().unwrap(), once);
    }

    #[rstest]
    fn test_deltas_normalize_dangling_bond_error() {
        let deltas = Deltas::from_iter([
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C),
            }),
            Delta::Bond(BondDelta::Add {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                attributes: BondForm::default(),
            }),
        ]);
        assert!(matches!(deltas.normalize(), Err(Contradiction)));
    }

    #[rstest]
    fn test_deltas_normalize_discontinuous_chain_error() {
        let deltas = Deltas::from_iter([charge_set(0, 0, 1), charge_set(0, 2, 3)]);
        assert!(matches!(deltas.normalize(), Err(Contradiction)));
    }

    #[rstest]
    #[case::add_remove(vec![
        Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(
            MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Lit(0) },
        ))),
        Delta::Constraint(ConstraintDelta::Remove(Constraint::Molecule(
            MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Lit(0) },
        ))),
    ])]
    #[case::duplicate_add_remove(vec![
        Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(
            MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Lit(0) },
        ))),
        Delta::Constraint(ConstraintDelta::Add(Constraint::Molecule(
            MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Lit(0) },
        ))),
        Delta::Constraint(ConstraintDelta::Remove(Constraint::Molecule(
            MoleculeConstraint::ChargeSum { atoms: None, sum: NumForm::Lit(0) },
        ))),
    ])]
    fn test_deltas_normalize_molecule_constraint_identity(#[case] values: Vec<Delta>) {
        assert_eq!(
            Deltas::from_iter(values).normalize().unwrap(),
            Deltas::new()
        );
    }

    #[rstest]
    #[case::duplicate_add(
        ConstraintDelta::Add(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: NumForm::Lit(0),
        })),
    )]
    #[case::duplicate_remove(
        ConstraintDelta::Remove(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: NumForm::Lit(0),
        })),
    )]
    fn test_deltas_normalize_molecule_constraint_set(#[case] value: ConstraintDelta) {
        let deltas = Deltas::from_iter([
            Delta::Constraint(value.clone()),
            Delta::Constraint(value.clone()),
        ]);
        assert_eq!(
            deltas.normalize().unwrap(),
            Deltas::from_iter([Delta::Constraint(value)]),
        );
    }
}
