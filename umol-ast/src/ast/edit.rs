//! Edit vocabulary for transactional molecule mutation (Phase 8).
//!
//! The `Edit` enum is the data-form vocabulary for `MoleculeBuilder::transact`
//! (Phase 8d). Each variant is self-inverting via [`Edit::inverse`]: `AddX`
//! ↔ `RemoveX`, and `Set*Field` swaps `old`/`new` on the inner `*FieldChange`.
//!
//! Refs (`AtomRef`, `BondRef`, ...) are symbolic and appear only inside
//! `Edit`. `Id(_)` references an existing entity; `New(N)` references the
//! entity created by the Nth Edit earlier in the same batch.

use super::aromatic::AromaticSystemAst;
use super::atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
use super::bond::BondAst;
use super::constraint::{
    AromaticSystemConstraint, AtomConstraint, BondConstraint, Constraint, DativeBondConstraint,
    MulticenterBondConstraint,
};
use super::dative::DativeBondAst;
use super::idx::{AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId};
use super::multicenter::MulticenterBondAst;
use super::noncovalent::{NoncovalentBondAst, NoncovalentBondKindAst};
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
    Element { old: ElementAst, new: ElementAst },
    IsotopeMass { old: IsotopeAst, new: IsotopeAst },
    Charge { old: ValueAst, new: ValueAst },
    ImplicitHydrogens { old: ImplicitHydrogensAst, new: ImplicitHydrogensAst },
    LonePairs { old: ValueAst, new: ValueAst },
    Spin { old: SpinStateAst, new: SpinStateAst },
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
    Order { old: ValueAst, new: ValueAst },
    Charge { old: ValueAst, new: ValueAst },
    Spin { old: SpinStateAst, new: SpinStateAst },
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
    Electrons { old: Vec<ValueAst>, new: Vec<ValueAst> },
    Charge { old: ValueAst, new: ValueAst },
    Spin { old: SpinStateAst, new: SpinStateAst },
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
    Electrons { old: Vec<ValueAst>, new: Vec<ValueAst> },
    Charge { old: ValueAst, new: ValueAst },
    Spin { old: SpinStateAst, new: SpinStateAst },
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
    Kind { old: NoncovalentBondKindAst, new: NoncovalentBondKindAst },
}

impl NoncovalentBondFieldChange {
    pub fn inverse(self) -> Self {
        match self {
            Self::Kind { old, new } => Self::Kind { old: new, new: old },
        }
    }
}

/// Single mutation operation. Compose `Vec<Edit>` into a transaction batch
/// (Phase 8d). Every variant carries the data needed to invert it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Edit {
    // Atoms / bonds
    AddAtom {
        ast: AtomAst,
    },
    RemoveAtom {
        idx: AtomRef,
        ast: AtomAst,
    },
    AddBond {
        a: AtomRef,
        b: AtomRef,
        ast: BondAst,
    },
    RemoveBond {
        idx: BondRef,
        endpoints: [AtomRef; 2],
        ast: BondAst,
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
    /// Edit that undoes `self`. Total: every variant has a defined inverse.
    /// `Remove*` returns `Add*` with the saved data; `Add*` returns
    /// `Remove*` with `*Ref::New(0)` as a sentinel (real id is unknown until
    /// apply time inside a transaction).
    pub fn inverse(self) -> Self {
        match self {
            Self::AddAtom { ast } => Self::RemoveAtom {
                idx: AtomRef::New(0),
                ast,
            },
            Self::RemoveAtom { idx: _, ast } => Self::AddAtom { ast },
            Self::AddBond { a, b, ast } => Self::RemoveBond {
                idx: BondRef::New(0),
                endpoints: [a, b],
                ast,
            },
            Self::RemoveBond {
                idx: _,
                endpoints: [a, b],
                ast,
            } => Self::AddBond { a, b, ast },
            Self::SetAtomField { idx, change } => Self::SetAtomField {
                idx,
                change: change.inverse(),
            },
            Self::SetBondField { idx, change } => Self::SetBondField {
                idx,
                change: change.inverse(),
            },

            Self::AddDativeBond { atoms, ast } => Self::RemoveDativeBond {
                idx: DativeBondRef::New(0),
                atoms,
                ast,
            },
            Self::RemoveDativeBond { idx: _, atoms, ast } => Self::AddDativeBond { atoms, ast },
            Self::SetDativeBondField { idx, change } => Self::SetDativeBondField {
                idx,
                change: change.inverse(),
            },

            Self::AddAromaticSystem { atoms, ast } => Self::RemoveAromaticSystem {
                idx: AromaticSystemRef::New(0),
                atoms,
                ast,
            },
            Self::RemoveAromaticSystem { idx: _, atoms, ast } => {
                Self::AddAromaticSystem { atoms, ast }
            }
            Self::SetAromaticSystemField { idx, change } => Self::SetAromaticSystemField {
                idx,
                change: change.inverse(),
            },

            Self::AddMulticenterBond { atoms, ast } => Self::RemoveMulticenterBond {
                idx: MulticenterBondRef::New(0),
                atoms,
                ast,
            },
            Self::RemoveMulticenterBond { idx: _, atoms, ast } => {
                Self::AddMulticenterBond { atoms, ast }
            }
            Self::SetMulticenterBondField { idx, change } => Self::SetMulticenterBondField {
                idx,
                change: change.inverse(),
            },

            Self::AddNoncovalentBond { atoms, ast } => Self::RemoveNoncovalentBond {
                idx: NoncovalentBondRef::New(0),
                atoms,
                ast,
            },
            Self::RemoveNoncovalentBond { idx: _, atoms, ast } => {
                Self::AddNoncovalentBond { atoms, ast }
            }
            Self::SetNoncovalentBondField { idx, change } => Self::SetNoncovalentBondField {
                idx,
                change: change.inverse(),
            },

            Self::SetAtomConstraint { idx, old, new } => Self::SetAtomConstraint {
                idx,
                old: new,
                new: old,
            },
            Self::AddAtomConstraint { idx, constraint } => Self::RemoveAtomConstraint {
                idx,
                constraint,
            },
            Self::RemoveAtomConstraint { idx, constraint } => Self::AddAtomConstraint {
                idx,
                constraint,
            },

            Self::SetBondConstraint { idx, old, new } => Self::SetBondConstraint {
                idx,
                old: new,
                new: old,
            },
            Self::AddBondConstraint { idx, constraint } => Self::RemoveBondConstraint {
                idx,
                constraint,
            },
            Self::RemoveBondConstraint { idx, constraint } => Self::AddBondConstraint {
                idx,
                constraint,
            },

            Self::SetDativeBondConstraint { idx, old, new } => Self::SetDativeBondConstraint {
                idx,
                old: new,
                new: old,
            },
            Self::AddDativeBondConstraint { idx, constraint } => Self::RemoveDativeBondConstraint {
                idx,
                constraint,
            },
            Self::RemoveDativeBondConstraint { idx, constraint } => {
                Self::AddDativeBondConstraint { idx, constraint }
            }

            Self::SetAromaticSystemConstraint { idx, old, new } => {
                Self::SetAromaticSystemConstraint {
                    idx,
                    old: new,
                    new: old,
                }
            }

            Self::SetMulticenterBondConstraint { idx, old, new } => {
                Self::SetMulticenterBondConstraint {
                    idx,
                    old: new,
                    new: old,
                }
            }

            Self::PushMoleculeConstraint { constraint } => {
                Self::PopMoleculeConstraint { constraint }
            }
            Self::PopMoleculeConstraint { constraint } => {
                Self::PushMoleculeConstraint { constraint }
            }
        }
    }
}

// Symbolic refs for overlay relations (Phase 8b).
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

/// Result of applying one `Edit` inside `transact`. `*Added` reports the new
/// id; `Done` is the no-new-id outcome of `Remove` / `Set*Field`; `Cascaded`
/// wraps the user's action plus the auto-generated R3-cascade Edits.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Action {
    AtomAdded(AtomId),
    BondAdded(BondId),
    AromaticSystemAdded(AromaticSystemId),
    DativeBondAdded(DativeBondId),
    MulticenterBondAdded(MulticenterBondId),
    NoncovalentBondAdded(NoncovalentBondId),
    Done,
    Cascaded {
        user: Box<Action>,
        cascade: Vec<Edit>,
    },
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use crate::ast::constraint::AtomConstraints;
    use crate::ast::value::Expr;
    use umol_shared::element::Element;

    #[fixture]
    fn carbon_atom() -> AtomAst {
        AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Undetermined,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ImplicitHydrogensAst::Lit(4),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::default(),
            constraints: AtomConstraints::new(),
        }
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
        let cloned = r.clone();
        assert_eq!(r, cloned);
    }

    #[rstest]
    #[case::id(BondRef::Id(BondId(5)))]
    #[case::new(BondRef::New(0))]
    fn test_bond_ref_variants(#[case] r: BondRef) {
        let cloned = r.clone();
        assert_eq!(r, cloned);
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
    #[case::isotope(
        AtomFieldChange::IsotopeMass {
            old: IsotopeAst::Undetermined,
            new: IsotopeAst::Lit(13),
        },
        AtomFieldChange::IsotopeMass {
            old: IsotopeAst::Lit(13),
            new: IsotopeAst::Undetermined,
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
    #[case::implicit_hydrogens(
        AtomFieldChange::ImplicitHydrogens {
            old: ImplicitHydrogensAst::Lit(3),
            new: ImplicitHydrogensAst::Lit(2),
        },
        AtomFieldChange::ImplicitHydrogens {
            old: ImplicitHydrogensAst::Lit(2),
            new: ImplicitHydrogensAst::Lit(3),
        },
    )]
    #[case::lone_pairs(
        AtomFieldChange::LonePairs {
            old: ValueAst::Lit(0),
            new: ValueAst::Lit(2),
        },
        AtomFieldChange::LonePairs {
            old: ValueAst::Lit(2),
            new: ValueAst::Lit(0),
        },
    )]
    #[case::spin(
        AtomFieldChange::Spin {
            old: SpinStateAst::default(),
            new: SpinStateAst {
                unpaired: ValueAst::Lit(1),
                multiplicity: ValueAst::Lit(2),
            },
        },
        AtomFieldChange::Spin {
            old: SpinStateAst {
                unpaired: ValueAst::Lit(1),
                multiplicity: ValueAst::Lit(2),
            },
            new: SpinStateAst::default(),
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
    fn test_edit_inverse_add_atom(carbon_atom: AtomAst) {
        let edit = Edit::AddAtom { ast: carbon_atom.clone() };
        assert_eq!(
            edit.inverse(),
            Edit::RemoveAtom {
                idx: AtomRef::New(0),
                ast: carbon_atom,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_remove_atom(carbon_atom: AtomAst) {
        let edit = Edit::RemoveAtom {
            idx: AtomRef::Id(AtomId(2)),
            ast: carbon_atom.clone(),
        };
        assert_eq!(edit.inverse(), Edit::AddAtom { ast: carbon_atom });
    }

    #[rstest]
    fn test_edit_inverse_add_bond(single_bond: BondAst) {
        let edit = Edit::AddBond {
            a: AtomRef::Id(AtomId(0)),
            b: AtomRef::Id(AtomId(1)),
            ast: single_bond.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::RemoveBond {
                idx: BondRef::New(0),
                endpoints: [AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast: single_bond,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_remove_bond(single_bond: BondAst) {
        let edit = Edit::RemoveBond {
            idx: BondRef::Id(BondId(4)),
            endpoints: [AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
            ast: single_bond.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::AddBond {
                a: AtomRef::Id(AtomId(0)),
                b: AtomRef::Id(AtomId(1)),
                ast: single_bond,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_set_atom_field() {
        let edit = Edit::SetAtomField {
            idx: AtomRef::Id(AtomId(7)),
            change: AtomFieldChange::Charge {
                old: ValueAst::Lit(0),
                new: ValueAst::Lit(-1),
            },
        };
        assert_eq!(
            edit.clone().inverse(),
            Edit::SetAtomField {
                idx: AtomRef::Id(AtomId(7)),
                change: AtomFieldChange::Charge {
                    old: ValueAst::Lit(-1),
                    new: ValueAst::Lit(0),
                },
            }
        );
        assert_eq!(edit.clone().inverse().inverse(), edit);
    }

    #[rstest]
    fn test_edit_inverse_set_atom_field_expr_payload() {
        let expr = ValueAst::Expr(Expr::Lit(2));
        let edit = Edit::SetAtomField {
            idx: AtomRef::New(1),
            change: AtomFieldChange::LonePairs {
                old: ValueAst::Lit(0),
                new: expr.clone(),
            },
        };
        assert_eq!(
            edit.inverse(),
            Edit::SetAtomField {
                idx: AtomRef::New(1),
                change: AtomFieldChange::LonePairs {
                    old: expr,
                    new: ValueAst::Lit(0),
                },
            }
        );
    }

    #[rstest]
    fn test_bond_field_change_inverse() {
        let c = BondFieldChange::Order {
            old: ValueAst::Lit(1),
            new: ValueAst::Lit(2),
        };
        assert_eq!(
            c.clone().inverse(),
            BondFieldChange::Order {
                old: ValueAst::Lit(2),
                new: ValueAst::Lit(1),
            }
        );
        assert_eq!(c.clone().inverse().inverse(), c);
    }

    #[rstest]
    fn test_dative_bond_field_change_inverse() {
        let c = DativeBondFieldChange::AcceptorSlot { old: 0, new: 1 };
        assert_eq!(
            c.clone().inverse(),
            DativeBondFieldChange::AcceptorSlot { old: 1, new: 0 }
        );
        assert_eq!(c.clone().inverse().inverse(), c);
    }

    #[rstest]
    fn test_aromatic_system_field_change_inverse() {
        let c = AromaticSystemFieldChange::Electrons {
            old: vec![ValueAst::Lit(1), ValueAst::Lit(1)],
            new: vec![ValueAst::Lit(2), ValueAst::Lit(0)],
        };
        assert_eq!(
            c.clone().inverse(),
            AromaticSystemFieldChange::Electrons {
                old: vec![ValueAst::Lit(2), ValueAst::Lit(0)],
                new: vec![ValueAst::Lit(1), ValueAst::Lit(1)],
            }
        );
        assert_eq!(c.clone().inverse().inverse(), c);
    }

    #[rstest]
    fn test_multicenter_bond_field_change_inverse() {
        let c = MulticenterBondFieldChange::Charge {
            old: ValueAst::Lit(0),
            new: ValueAst::Lit(-1),
        };
        assert_eq!(
            c.clone().inverse(),
            MulticenterBondFieldChange::Charge {
                old: ValueAst::Lit(-1),
                new: ValueAst::Lit(0),
            }
        );
        assert_eq!(c.clone().inverse().inverse(), c);
    }

    #[rstest]
    fn test_noncovalent_bond_field_change_inverse() {
        use crate::ast::noncovalent::NoncovalentBondKind;
        let c = NoncovalentBondFieldChange::Kind {
            old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
            new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::VanDerWaals),
        };
        assert_eq!(
            c.clone().inverse(),
            NoncovalentBondFieldChange::Kind {
                old: NoncovalentBondKindAst::Lit(NoncovalentBondKind::VanDerWaals),
                new: NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond),
            }
        );
        assert_eq!(c.clone().inverse().inverse(), c);
    }

    #[rstest]
    fn test_edit_inverse_set_bond_field() {
        let edit = Edit::SetBondField {
            idx: BondRef::Id(BondId(2)),
            change: BondFieldChange::Order {
                old: ValueAst::Lit(1),
                new: ValueAst::Lit(2),
            },
        };
        assert_eq!(
            edit.clone().inverse(),
            Edit::SetBondField {
                idx: BondRef::Id(BondId(2)),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(2),
                    new: ValueAst::Lit(1),
                },
            }
        );
        assert_eq!(edit.clone().inverse().inverse(), edit);
    }

    #[rstest]
    fn test_edit_inverse_add_dative_bond() {
        let ast = DativeBondAst::default();
        let edit = Edit::AddDativeBond {
            atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
            ast: ast.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::RemoveDativeBond {
                idx: DativeBondRef::New(0),
                atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_remove_dative_bond() {
        let ast = DativeBondAst::default();
        let edit = Edit::RemoveDativeBond {
            idx: DativeBondRef::Id(DativeBondId(3)),
            atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
            ast: ast.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::AddDativeBond {
                atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_add_aromatic_system() {
        let ast = AromaticSystemAst::default();
        let edit = Edit::AddAromaticSystem {
            atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1)), AtomRef::Id(AtomId(2))],
            ast: ast.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::RemoveAromaticSystem {
                idx: AromaticSystemRef::New(0),
                atoms: vec![
                    AtomRef::Id(AtomId(0)),
                    AtomRef::Id(AtomId(1)),
                    AtomRef::Id(AtomId(2)),
                ],
                ast,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_add_multicenter_bond() {
        let ast = MulticenterBondAst::default();
        let edit = Edit::AddMulticenterBond {
            atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
            ast: ast.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::RemoveMulticenterBond {
                idx: MulticenterBondRef::New(0),
                atoms: vec![AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(1))],
                ast,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_add_noncovalent_bond() {
        use crate::ast::noncovalent::NoncovalentBondKind;
        let ast = NoncovalentBondAst::new(NoncovalentBondKindAst::Lit(
            NoncovalentBondKind::HydrogenBond,
        ));
        let edit = Edit::AddNoncovalentBond {
            atoms: [AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(3))],
            ast: ast.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::RemoveNoncovalentBond {
                idx: NoncovalentBondRef::New(0),
                atoms: [AtomRef::Id(AtomId(0)), AtomRef::Id(AtomId(3))],
                ast,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_push_molecule_constraint() {
        use crate::ast::constraint::MoleculeConstraint;
        let c = Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
        });
        let edit = Edit::PushMoleculeConstraint { constraint: c.clone() };
        assert_eq!(
            edit.clone().inverse(),
            Edit::PopMoleculeConstraint { constraint: c.clone() }
        );
        assert_eq!(edit.clone().inverse().inverse(), edit);
    }

    #[rstest]
    fn test_edit_inverse_pop_molecule_constraint() {
        use crate::ast::constraint::MoleculeConstraint;
        let c = Constraint::Molecule(MoleculeConstraint::Connected {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
        });
        let edit = Edit::PopMoleculeConstraint { constraint: c.clone() };
        assert_eq!(
            edit.inverse(),
            Edit::PushMoleculeConstraint { constraint: c }
        );
    }

    #[rstest]
    fn test_edit_inverse_set_atom_constraint() {
        let edit = Edit::SetAtomConstraint {
            idx: AtomRef::Id(AtomId(2)),
            old: Some(AtomConstraint::valence(3)),
            new: Some(AtomConstraint::valence(4)),
        };
        assert_eq!(
            edit.clone().inverse(),
            Edit::SetAtomConstraint {
                idx: AtomRef::Id(AtomId(2)),
                old: Some(AtomConstraint::valence(4)),
                new: Some(AtomConstraint::valence(3)),
            }
        );
        assert_eq!(edit.clone().inverse().inverse(), edit);
    }

    #[rstest]
    fn test_edit_inverse_set_atom_constraint_introduce() {
        let edit = Edit::SetAtomConstraint {
            idx: AtomRef::Id(AtomId(0)),
            old: None,
            new: Some(AtomConstraint::valence(4)),
        };
        assert_eq!(
            edit.inverse(),
            Edit::SetAtomConstraint {
                idx: AtomRef::Id(AtomId(0)),
                old: Some(AtomConstraint::valence(4)),
                new: None,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_add_atom_constraint() {
        let c = AtomConstraint::ring_size(5);
        let edit = Edit::AddAtomConstraint {
            idx: AtomRef::Id(AtomId(1)),
            constraint: c.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::RemoveAtomConstraint {
                idx: AtomRef::Id(AtomId(1)),
                constraint: c,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_remove_atom_constraint() {
        let c = AtomConstraint::ring_size(6);
        let edit = Edit::RemoveAtomConstraint {
            idx: AtomRef::Id(AtomId(1)),
            constraint: c.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::AddAtomConstraint {
                idx: AtomRef::Id(AtomId(1)),
                constraint: c,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_set_bond_constraint_aromatic_marker() {
        let edit = Edit::SetBondConstraint {
            idx: BondRef::Id(BondId(0)),
            old: None,
            new: Some(BondConstraint::Aromatic),
        };
        assert_eq!(
            edit.inverse(),
            Edit::SetBondConstraint {
                idx: BondRef::Id(BondId(0)),
                old: Some(BondConstraint::Aromatic),
                new: None,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_add_bond_constraint() {
        let c = BondConstraint::RingSize(ValueAst::Lit(6));
        let edit = Edit::AddBondConstraint {
            idx: BondRef::Id(BondId(2)),
            constraint: c.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::RemoveBondConstraint {
                idx: BondRef::Id(BondId(2)),
                constraint: c,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_set_dative_bond_constraint() {
        let edit = Edit::SetDativeBondConstraint {
            idx: DativeBondRef::Id(DativeBondId(0)),
            old: None,
            new: Some(DativeBondConstraint::Aromatic),
        };
        assert_eq!(
            edit.inverse(),
            Edit::SetDativeBondConstraint {
                idx: DativeBondRef::Id(DativeBondId(0)),
                old: Some(DativeBondConstraint::Aromatic),
                new: None,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_add_dative_bond_constraint() {
        let c = DativeBondConstraint::RingSize(ValueAst::Lit(5));
        let edit = Edit::AddDativeBondConstraint {
            idx: DativeBondRef::Id(DativeBondId(0)),
            constraint: c.clone(),
        };
        assert_eq!(
            edit.inverse(),
            Edit::RemoveDativeBondConstraint {
                idx: DativeBondRef::Id(DativeBondId(0)),
                constraint: c,
            }
        );
    }

    #[rstest]
    fn test_edit_inverse_set_aromatic_system_constraint() {
        let edit = Edit::SetAromaticSystemConstraint {
            idx: AromaticSystemRef::Id(AromaticSystemId(0)),
            old: Some(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))),
            new: Some(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(10))),
        };
        assert_eq!(
            edit.clone().inverse(),
            Edit::SetAromaticSystemConstraint {
                idx: AromaticSystemRef::Id(AromaticSystemId(0)),
                old: Some(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(10))),
                new: Some(AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))),
            }
        );
        assert_eq!(edit.clone().inverse().inverse(), edit);
    }

    #[rstest]
    fn test_edit_inverse_set_multicenter_bond_constraint() {
        let edit = Edit::SetMulticenterBondConstraint {
            idx: MulticenterBondRef::Id(MulticenterBondId(0)),
            old: None,
            new: Some(MulticenterBondConstraint::ElectronCount(ValueAst::Lit(4))),
        };
        assert_eq!(
            edit.inverse(),
            Edit::SetMulticenterBondConstraint {
                idx: MulticenterBondRef::Id(MulticenterBondId(0)),
                old: Some(MulticenterBondConstraint::ElectronCount(ValueAst::Lit(4))),
                new: None,
            }
        );
    }

    #[rstest]
    fn test_action_cascaded_nested() {
        let inner = Action::Cascaded {
            user: Box::new(Action::AtomAdded(AtomId(0))),
            cascade: vec![Edit::RemoveAtom {
                idx: AtomRef::Id(AtomId(1)),
                ast: AtomAst::default(),
            }],
        };
        let outer = Action::Cascaded {
            user: Box::new(inner.clone()),
            cascade: vec![],
        };
        let Action::Cascaded { user, cascade } = outer else {
            unreachable!()
        };
        assert_eq!(*user, inner);
        assert!(cascade.is_empty());
    }
}
