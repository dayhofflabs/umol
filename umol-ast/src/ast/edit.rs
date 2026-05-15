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
use super::atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
use super::bond::BondAst;
use super::constraint::{
    AromaticSystemConstraint, AtomConstraint, BondConstraint, Constraint, Constraints,
    DativeBondConstraint, MulticenterBondConstraint,
};
use super::dative::DativeBondAst;
use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::{NoncovalentBondAst, NoncovalentBondKindAst};
use super::remap::{IdRemapping, UndoRemapping};
use super::spin::SpinStateAst;
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
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomFieldChange {
    Element {
        old: ElementAst,
        new: ElementAst,
    },
    IsotopeMass {
        old: IsotopeAst,
        new: IsotopeAst,
    },
    Charge {
        old: ValueAst,
        new: ValueAst,
    },
    ImplicitHydrogens {
        old: ImplicitHydrogensAst,
        new: ImplicitHydrogensAst,
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

#[derive(Clone, Debug, PartialEq, Eq)]
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DativeBondFieldChange {
    AcceptorSlot { old: u8, new: u8 },
    Order { old: ValueAst, new: ValueAst },
}

impl DativeBondFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::AcceptorSlot { old, new } => Self::AcceptorSlot { old: new, new: old },
            Self::Order { old, new } => Self::Order { old: new, new: old },
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AromaticSystemFieldChange {
    Electrons {
        old: Vec<ValueAst>,
        new: Vec<ValueAst>,
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MulticenterBondFieldChange {
    Electrons {
        old: Vec<ValueAst>,
        new: Vec<ValueAst>,
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

#[derive(Clone, Debug, PartialEq, Eq)]
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

/// Single bond addition inside an `Edit::AddBonds` batch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddBond {
    pub a: AtomRef,
    pub b: AtomRef,
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
    SetAtomField {
        idx: AtomRef,
        change: AtomFieldChange,
    },
    SetBondField {
        idx: BondRef,
        change: BondFieldChange,
    },

    // Dative bonds
    AddDativeBond {
        atoms: Vec<AtomRef>,
        ast: DativeBondAst,
    },
    RemoveDativeBond {
        idx: DativeBondRef,
        atoms: Vec<AtomRef>,
        ast: DativeBondAst,
    },
    SetDativeBondField {
        idx: DativeBondRef,
        change: DativeBondFieldChange,
    },

    // Aromatic systems
    AddAromaticSystem {
        atoms: Vec<AtomRef>,
        ast: AromaticSystemAst,
    },
    RemoveAromaticSystem {
        idx: AromaticSystemRef,
        atoms: Vec<AtomRef>,
        ast: AromaticSystemAst,
    },
    SetAromaticSystemField {
        idx: AromaticSystemRef,
        change: AromaticSystemFieldChange,
    },

    // Multicenter bonds
    AddMulticenterBond {
        atoms: Vec<AtomRef>,
        ast: MulticenterBondAst,
    },
    RemoveMulticenterBond {
        idx: MulticenterBondRef,
        atoms: Vec<AtomRef>,
        ast: MulticenterBondAst,
    },
    SetMulticenterBondField {
        idx: MulticenterBondRef,
        change: MulticenterBondFieldChange,
    },

    // Noncovalent bonds
    AddNoncovalentBond {
        atoms: [AtomRef; 2],
        ast: NoncovalentBondAst,
    },
    RemoveNoncovalentBond {
        idx: NoncovalentBondRef,
        atoms: [AtomRef; 2],
        ast: NoncovalentBondAst,
    },
    SetNoncovalentBondField {
        idx: NoncovalentBondRef,
        change: NoncovalentBondFieldChange,
    },

    // Entity-inline constraints — atom
    SetAtomConstraint {
        idx: AtomRef,
        old: Option<AtomConstraint>,
        new: Option<AtomConstraint>,
    },
    AddAtomConstraint {
        idx: AtomRef,
        constraint: AtomConstraint,
    },
    RemoveAtomConstraint {
        idx: AtomRef,
        constraint: AtomConstraint,
    },

    // Entity-inline constraints — bond
    SetBondConstraint {
        idx: BondRef,
        old: Option<BondConstraint>,
        new: Option<BondConstraint>,
    },
    AddBondConstraint {
        idx: BondRef,
        constraint: BondConstraint,
    },
    RemoveBondConstraint {
        idx: BondRef,
        constraint: BondConstraint,
    },

    // Entity-inline constraints — dative bond
    SetDativeBondConstraint {
        idx: DativeBondRef,
        old: Option<DativeBondConstraint>,
        new: Option<DativeBondConstraint>,
    },
    AddDativeBondConstraint {
        idx: DativeBondRef,
        constraint: DativeBondConstraint,
    },
    RemoveDativeBondConstraint {
        idx: DativeBondRef,
        constraint: DativeBondConstraint,
    },

    // Entity-inline constraints — aromatic system (no non-unique kinds)
    SetAromaticSystemConstraint {
        idx: AromaticSystemRef,
        old: Option<AromaticSystemConstraint>,
        new: Option<AromaticSystemConstraint>,
    },

    // Entity-inline constraints — multicenter bond (no non-unique kinds)
    SetMulticenterBondConstraint {
        idx: MulticenterBondRef,
        old: Option<MulticenterBondConstraint>,
        new: Option<MulticenterBondConstraint>,
    },

    // Molecule-list constraints (stack discipline; arbitrary-position removal
    // not part of the Edit grammar — use `constraints_mut().remove_at()` for
    // that, outside transact).
    PushMoleculeConstraint {
        constraint: Constraint,
    },
    PopMoleculeConstraint {
        constraint: Constraint,
    },
}

impl Edit {
    pub fn add_atom(ast: AtomAst) -> Self {
        Self::AddAtoms { atoms: vec![ast] }
    }

    pub fn add_bond(a: AtomRef, b: AtomRef, ast: BondAst) -> Self {
        Self::AddBonds {
            bonds: vec![AddBond { a, b, ast }],
        }
    }

    pub fn remove_atom(idx: AtomRef) -> Self {
        Self::RemoveTopology {
            atoms: vec![idx],
            bonds: Vec::new(),
        }
    }

    pub fn remove_bond(idx: BondRef) -> Self {
        Self::RemoveTopology {
            atoms: Vec::new(),
            bonds: vec![idx],
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

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RemovedOverlays {
    pub dative_bonds: Vec<RemovedDativeBond>,
    pub aromatic_systems: Vec<RemovedAromaticSystem>,
    pub multicenter_bonds: Vec<RemovedMulticenterBond>,
    pub noncovalent_bonds: Vec<RemovedNoncovalentBond>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DroppedConstraint {
    pub position: usize,
    pub constraint: Constraint,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RewrittenConstraint {
    pub position: usize,
    pub old: Constraint,
    pub new: Constraint,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ConstraintUpdate {
    pub dropped: Vec<DroppedConstraint>,
    pub rewritten: Vec<RewrittenConstraint>,
}

impl ConstraintUpdate {
    pub fn is_empty(&self) -> bool {
        self.dropped.is_empty() && self.rewritten.is_empty()
    }

    pub fn rollback_into(self, constraints: &mut Constraints) {
        let mut items = constraints.take();
        for rewritten in self.rewritten {
            if let Some(pos) = items.iter().position(|c| *c == rewritten.new) {
                items[pos] = rewritten.old;
            }
        }
        for dropped in self.dropped {
            let position = dropped.position.min(items.len());
            items.insert(position, dropped.constraint);
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
    RestoreTopology {
        atoms: Vec<RemovedAtom>,
        bonds: Vec<RemovedBond>,
        overlays: RemovedOverlays,
        remapping: IdRemapping,
        undo_remapping: UndoRemapping,
        constraint_update: ConstraintUpdate,
    },
    RemoveAddedDativeBond(AddedDativeBond),
    RestoreRemovedDativeBond {
        removed: RemovedDativeBond,
        undo_remapping: UndoRemapping,
    },
    RemoveAddedAromaticSystem(AddedAromaticSystem),
    RestoreRemovedAromaticSystem {
        removed: RemovedAromaticSystem,
        undo_remapping: UndoRemapping,
    },
    RemoveAddedMulticenterBond(AddedMulticenterBond),
    RestoreRemovedMulticenterBond {
        removed: RemovedMulticenterBond,
        undo_remapping: UndoRemapping,
    },
    RemoveAddedNoncovalentBond(AddedNoncovalentBond),
    RestoreRemovedNoncovalentBond {
        removed: RemovedNoncovalentBond,
        undo_remapping: UndoRemapping,
    },
    SetAtomField {
        id: AtomId,
        change: AtomFieldChange,
    },
    SetBondField {
        id: BondId,
        change: BondFieldChange,
    },
    SetDativeBondField {
        id: DativeBondId,
        change: DativeBondFieldChange,
    },
    SetAromaticSystemField {
        id: AromaticSystemId,
        change: AromaticSystemFieldChange,
    },
    SetMulticenterBondField {
        id: MulticenterBondId,
        change: MulticenterBondFieldChange,
    },
    SetNoncovalentBondField {
        id: NoncovalentBondId,
        change: NoncovalentBondFieldChange,
    },
    ApplyConstraintUpdate(ConstraintUpdate),
    ApplyEdit(Box<Edit>),
}

impl Undo {
    pub fn id_remapping(&self) -> Option<&IdRemapping> {
        match self {
            Self::RestoreTopology { remapping, .. } => Some(remapping),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::element::Element;

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
                    a: AtomRef::Id(AtomId(0)),
                    b: AtomRef::Id(AtomId(1)),
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
