//! Edit vocabulary for transactional molecule mutation.
//!
//! The `Edit` enum is the data-form vocabulary for `MoleculeBuilder::transact`
//! `Edit` is caller-facing mutation data; realized rollback data belongs to
//! the `Undo` journal.
//!
//! Refs (`AtomRef`, `BondRef`, ...) are symbolic and appear only inside
//! `Edit`. `Id(_)` references an existing entity; `New(N)` references the
//! entity created by the Nth Edit earlier in the same batch.

use super::aromatic::AromaticSystemAst;
use super::atom::{AtomAst, ElementAst, IsotopeMassAst};
use super::bond::BondAst;
use super::constraint::{
    AromaticSystemConstraint, AtomConstraint, BondConstraint, Constraint, Constraints,
    DativeBondConstraint, MulticenterBondConstraint,
};
use super::dative::DativeBondAst;
use super::electrons::ElectronCountsAst;
use super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::ligand::{StereoLigand, StereoLigandKind};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::{NoncovalentBondAst, NoncovalentBondKindAst};
use super::remap::{IdCompaction, UndoCompaction};
use super::spin::SpinStateAst;
use super::stereo::{StereoAtomAst, StereoBondAst, StereoConfigurationAst};
use super::value::ValueAst;

/// Symbolic reference to an atom: either an existing `AtomId` or the Nth
/// atom-creating Edit earlier in the same transaction batch.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AtomRef {
    Id(AtomId),
    New(usize),
}

/// Symbolic reference to a bond.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum BondRef {
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

/// Per-field old/new payload for a stereo-atom mutation. Only `coset` is
/// settable: `kind` fixes the coset's group, so changing it would desync the
/// configuration — kind changes go through remove + add.
#[derive(Clone, Debug, PartialEq, Eq)]
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
#[derive(Clone, Debug, PartialEq, Eq)]
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
    pub endpoints: [AtomRef; 2],
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
        atoms: Vec<AtomRef>,
        bonds: Vec<BondRef>,
    },
    ModifyAtomField {
        id: AtomRef,
        change: AtomFieldChange,
    },
    ModifyBondField {
        id: BondRef,
        change: BondFieldChange,
    },

    // Dative bonds
    AddDativeBond {
        atoms: Vec<AtomRef>,
        ast: DativeBondAst,
    },
    RemoveDativeBond {
        id: DativeBondRef,
        atoms: Vec<AtomRef>,
        ast: DativeBondAst,
    },
    ModifyDativeBondField {
        id: DativeBondRef,
        change: DativeBondFieldChange,
    },

    // Aromatic systems
    AddAromaticSystem {
        atoms: Vec<AtomRef>,
        ast: AromaticSystemAst,
    },
    RemoveAromaticSystem {
        id: AromaticSystemRef,
        atoms: Vec<AtomRef>,
        ast: AromaticSystemAst,
    },
    ModifyAromaticSystemField {
        id: AromaticSystemRef,
        change: AromaticSystemFieldChange,
    },

    // Multicenter bonds
    AddMulticenterBond {
        atoms: Vec<AtomRef>,
        ast: MulticenterBondAst,
    },
    RemoveMulticenterBond {
        id: MulticenterBondRef,
        atoms: Vec<AtomRef>,
        ast: MulticenterBondAst,
    },
    ModifyMulticenterBondField {
        id: MulticenterBondRef,
        change: MulticenterBondFieldChange,
    },

    // Noncovalent bonds
    AddNoncovalentBond {
        atoms: [AtomRef; 2],
        ast: NoncovalentBondAst,
    },
    RemoveNoncovalentBond {
        id: NoncovalentBondRef,
        atoms: [AtomRef; 2],
        ast: NoncovalentBondAst,
    },
    ModifyNoncovalentBondField {
        id: NoncovalentBondRef,
        change: NoncovalentBondFieldChange,
    },

    // Stereo elements. `ligands` carry their atom as an `AtomRef` (Id or
    // same-batch New) plus the ligand kind; `site` is the atom/bond the
    // element is sited on.
    AddStereoAtom {
        site: AtomRef,
        ligands: Vec<(AtomRef, StereoLigandKind)>,
        ast: StereoAtomAst,
    },
    RemoveStereoAtom {
        id: StereoAtomRef,
        site: AtomRef,
        ligands: Vec<(AtomRef, StereoLigandKind)>,
        ast: StereoAtomAst,
    },
    ModifyStereoAtomField {
        id: StereoAtomRef,
        change: StereoAtomFieldChange,
    },
    AddStereoBond {
        site: BondRef,
        ligands: Vec<(AtomRef, StereoLigandKind)>,
        ast: StereoBondAst,
    },
    RemoveStereoBond {
        id: StereoBondRef,
        site: BondRef,
        ligands: Vec<(AtomRef, StereoLigandKind)>,
        ast: StereoBondAst,
    },
    ModifyStereoBondField {
        id: StereoBondRef,
        change: StereoBondFieldChange,
    },

    // Entity-inline constraints — keyed (one per `key()`), so a single modify
    // (old → new) covers add (old None), remove (new None), and replace.
    ModifyAtomConstraint {
        id: AtomRef,
        old: Option<AtomConstraint>,
        new: Option<AtomConstraint>,
    },
    ModifyBondConstraint {
        id: BondRef,
        old: Option<BondConstraint>,
        new: Option<BondConstraint>,
    },
    ModifyDativeBondConstraint {
        id: DativeBondRef,
        old: Option<DativeBondConstraint>,
        new: Option<DativeBondConstraint>,
    },
    ModifyAromaticSystemConstraint {
        id: AromaticSystemRef,
        old: Option<AromaticSystemConstraint>,
        new: Option<AromaticSystemConstraint>,
    },
    ModifyMulticenterBondConstraint {
        id: MulticenterBondRef,
        old: Option<MulticenterBondConstraint>,
        new: Option<MulticenterBondConstraint>,
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

    pub fn add_bond(first: AtomRef, second: AtomRef, ast: BondAst) -> Self {
        Self::AddBonds {
            bonds: vec![AddBond {
                endpoints: [first, second],
                ast,
            }],
        }
    }

    pub fn remove_atom(id: AtomRef) -> Self {
        Self::RemoveTopology {
            atoms: vec![id],
            bonds: Vec::new(),
        }
    }

    pub fn remove_bond(id: BondRef) -> Self {
        Self::RemoveTopology {
            atoms: Vec::new(),
            bonds: vec![id],
        }
    }
}

// Symbolic refs for overlay relations.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum DativeBondRef {
    Id(DativeBondId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemRef {
    Id(AromaticSystemId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MulticenterBondRef {
    Id(MulticenterBondId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondRef {
    Id(NoncovalentBondId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum StereoAtomRef {
    Id(StereoAtomId),
    New(usize),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum StereoBondRef {
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
        remapping: IdCompaction,
        undo_compaction: UndoCompaction,
        cascade: CascadedConstraints,
    },
    RemoveAddedDativeBond(AddedDativeBond),
    RestoreRemovedDativeBond {
        removed: RemovedDativeBond,
        undo_compaction: UndoCompaction,
    },
    RemoveAddedAromaticSystem(AddedAromaticSystem),
    RestoreRemovedAromaticSystem {
        removed: RemovedAromaticSystem,
        undo_compaction: UndoCompaction,
    },
    RemoveAddedMulticenterBond(AddedMulticenterBond),
    RestoreRemovedMulticenterBond {
        removed: RemovedMulticenterBond,
        undo_compaction: UndoCompaction,
    },
    RemoveAddedNoncovalentBond(AddedNoncovalentBond),
    RestoreRemovedNoncovalentBond {
        removed: RemovedNoncovalentBond,
        undo_compaction: UndoCompaction,
    },
    RemoveAddedStereoAtom(AddedStereoAtom),
    RestoreRemovedStereoAtom {
        removed: RemovedStereoAtom,
        undo_compaction: UndoCompaction,
    },
    RemoveAddedStereoBond(AddedStereoBond),
    RestoreRemovedStereoBond {
        removed: RemovedStereoBond,
        undo_compaction: UndoCompaction,
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
    pub fn id_remapping(&self) -> Option<&IdCompaction> {
        match self {
            Self::RestoreRemovedTopology { remapping, .. } => Some(remapping),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::stereo::{StereoConfigurationAst, StereoCosetAst, StereoKind};
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
    #[case::id(AtomRef::Id(AtomId(3)))]
    #[case::new(AtomRef::New(2))]
    fn test_atom_ref_variants(#[case] r: AtomRef) {
        assert_eq!(r.clone(), r);
    }

    #[rstest]
    #[case::id(BondRef::Id(BondId(5)))]
    #[case::new(BondRef::New(0))]
    fn test_bond_ref_variants(#[case] r: BondRef) {
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
                AtomRef::Id(AtomId(0)),
                AtomRef::Id(AtomId(1)),
                single_bond.clone()
            ),
            Edit::AddBonds {
                bonds: vec![AddBond {
                    endpoints: [AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                    ast: single_bond,
                }],
            },
        );
    }

    #[rstest]
    fn test_edit_remove_atom() {
        assert_eq!(
            Edit::remove_atom(AtomRef::Id(AtomId(2))),
            Edit::RemoveTopology {
                atoms: vec![AtomRef::Id(AtomId(2))],
                bonds: Vec::new(),
            },
        );
    }

    #[rstest]
    fn test_edit_remove_bond() {
        assert_eq!(
            Edit::remove_bond(BondRef::Id(BondId(4))),
            Edit::RemoveTopology {
                atoms: Vec::new(),
                bonds: vec![BondRef::Id(BondId(4))],
            },
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
