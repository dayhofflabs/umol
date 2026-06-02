//! Cross-entity relational constraints.
//!
//! A `RelationalConstraint` relates one DAMN entity (dative bond, aromatic
//! system, multicenter bond, noncovalent bond) to one or more other entities
//! by reference. Every variant carries at least two indices (the outer entity
//! plus one or more inner atom/bond refs) or one index plus a delegated
//! predicate over a role slot.
//!
//! Relational constraints live **only** at molecule scope — as entries in
//! `MoleculeAst::constraints` (via `Constraint::Relational(...)`) or inside
//! `And`/`Or`/`Not` combinators. They cannot be inline on the entity AST:
//! the per-entity `XxxConstraints` containers are narrowed to value-only
//! variants so ref-bearing constraints are unrepresentable inline.
//!
//! Two sub-patterns share the enum:
//! - **Role identity**: `Donor`, `Acceptor`, `Parallels`, `Ends`, `Atoms`,
//!   `Contains`, `ContainsAll` — constrain an atom/bond identity to a role
//!   or set membership.
//! - **Role predicate**: `DonorSatisfies`, `AcceptorSatisfies`, `AllAtoms`,
//!   `AnyAtom`, `EndsSatisfy` — delegate an `AtomConstraint` to a role slot,
//!   quantified over the matching participants.

use super::super::ids::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::remap::IdRemapping;
use super::atom::AtomConstraint;

/// Cross-entity constraint relating one DAMN entity (dative bond, aromatic
/// system, multicenter bond, noncovalent bond) to atoms, bonds, or atom
/// predicates by reference. Lives only at molecule scope (in
/// [`Constraint::Relational`](super::Constraint::Relational) or inside
/// `And`/`Or`/`Not`); cannot appear inline on an entity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RelationalConstraint {
    /// The donor end of dative bond `bond` is `atom`.
    DativeBondDonor { bond: DativeBondId, atom: AtomId },
    /// The acceptor end of dative bond `bond` is `atom`.
    DativeBondAcceptor { bond: DativeBondId, atom: AtomId },
    /// Dative bond `dative` is parallel to a localized bond `parallel`
    /// (same atom pair).
    DativeBondParallels {
        dative: DativeBondId,
        parallel: BondId,
    },
    /// The donor end of dative bond `bond` satisfies `predicate`.
    DativeBondDonorSatisfies {
        bond: DativeBondId,
        predicate: Box<AtomConstraint>,
    },
    /// The acceptor end of dative bond `bond` satisfies `predicate`.
    DativeBondAcceptorSatisfies {
        bond: DativeBondId,
        predicate: Box<AtomConstraint>,
    },

    /// Aromatic system `system` consists of exactly `atoms` (as a set).
    AromaticSystemAtoms {
        system: AromaticSystemId,
        atoms: Vec<AtomId>,
    },
    /// `atom` is one of the atoms participating in aromatic system `system`.
    AromaticSystemContains {
        system: AromaticSystemId,
        atom: AtomId,
    },
    /// Every atom in `atoms` participates in aromatic system `system`
    /// (set inclusion, not equality).
    AromaticSystemContainsAll {
        system: AromaticSystemId,
        atoms: Vec<AtomId>,
    },
    /// Every atom of aromatic system `system` satisfies `predicate`.
    AromaticSystemAllAtoms {
        system: AromaticSystemId,
        predicate: Box<AtomConstraint>,
    },
    /// At least one atom of aromatic system `system` satisfies `predicate`.
    AromaticSystemAnyAtom {
        system: AromaticSystemId,
        predicate: Box<AtomConstraint>,
    },

    /// Multicenter bond `bond` consists of exactly `atoms` (as a set).
    MulticenterBondAtoms {
        bond: MulticenterBondId,
        atoms: Vec<AtomId>,
    },
    /// `atom` is one of the participants in multicenter bond `bond`.
    MulticenterBondContains {
        bond: MulticenterBondId,
        atom: AtomId,
    },
    /// Every atom in `atoms` participates in multicenter bond `bond`.
    MulticenterBondContainsAll {
        bond: MulticenterBondId,
        atoms: Vec<AtomId>,
    },
    /// Every participating atom of multicenter bond `bond` satisfies
    /// `predicate`.
    MulticenterBondAllAtoms {
        bond: MulticenterBondId,
        predicate: Box<AtomConstraint>,
    },
    /// At least one participating atom of multicenter bond `bond` satisfies
    /// `predicate`.
    MulticenterBondAnyAtom {
        bond: MulticenterBondId,
        predicate: Box<AtomConstraint>,
    },

    /// Noncovalent bond `bond` connects exactly the pair `atoms` (unordered).
    NoncovalentBondEnds {
        bond: NoncovalentBondId,
        atoms: [AtomId; 2],
    },
    /// `atom` is one of the two endpoints of noncovalent bond `bond`.
    NoncovalentBondContains {
        bond: NoncovalentBondId,
        atom: AtomId,
    },
    /// The two endpoints of noncovalent bond `bond` satisfy `predicates[0]`
    /// and `predicates[1]` respectively. Order is not symmetric: the bond
    /// stores its endpoints as an unordered pair, but each predicate is
    /// associated with one specific slot.
    NoncovalentBondEndsSatisfy {
        bond: NoncovalentBondId,
        predicates: [Box<AtomConstraint>; 2],
    },
}

impl RelationalConstraint {
    /// Simplify any inner `AtomConstraint` predicate's `ValueAst`. Refs
    /// (`bond`, `atom`, `system`, `parallel`) are unchanged.
    pub fn simplify(self) -> Self {
        match self {
            Self::DativeBondDonorSatisfies { bond, predicate } => Self::DativeBondDonorSatisfies {
                bond,
                predicate: Box::new((*predicate).simplify()),
            },
            Self::DativeBondAcceptorSatisfies { bond, predicate } => {
                Self::DativeBondAcceptorSatisfies {
                    bond,
                    predicate: Box::new((*predicate).simplify()),
                }
            }
            Self::AromaticSystemAllAtoms { system, predicate } => Self::AromaticSystemAllAtoms {
                system,
                predicate: Box::new((*predicate).simplify()),
            },
            Self::AromaticSystemAnyAtom { system, predicate } => Self::AromaticSystemAnyAtom {
                system,
                predicate: Box::new((*predicate).simplify()),
            },
            Self::MulticenterBondAllAtoms { bond, predicate } => Self::MulticenterBondAllAtoms {
                bond,
                predicate: Box::new((*predicate).simplify()),
            },
            Self::MulticenterBondAnyAtom { bond, predicate } => Self::MulticenterBondAnyAtom {
                bond,
                predicate: Box::new((*predicate).simplify()),
            },
            Self::NoncovalentBondEndsSatisfy { bond, predicates } => {
                let [a, b] = predicates;
                Self::NoncovalentBondEndsSatisfy {
                    bond,
                    predicates: [Box::new((*a).simplify()), Box::new((*b).simplify())],
                }
            }
            other @ (Self::DativeBondDonor { .. }
            | Self::DativeBondAcceptor { .. }
            | Self::DativeBondParallels { .. }
            | Self::AromaticSystemAtoms { .. }
            | Self::AromaticSystemContains { .. }
            | Self::AromaticSystemContainsAll { .. }
            | Self::MulticenterBondAtoms { .. }
            | Self::MulticenterBondContains { .. }
            | Self::MulticenterBondContainsAll { .. }
            | Self::NoncovalentBondEnds { .. }
            | Self::NoncovalentBondContains { .. }) => other,
        }
    }

    /// Remap all indices this constraint carries. Returns `None` if any
    /// referenced entity has been removed by the remapping.
    pub fn remap(self, remap: &IdRemapping) -> Option<Self> {
        Some(match self {
            Self::DativeBondDonor { bond, atom } => Self::DativeBondDonor {
                bond: remap.dative_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            Self::DativeBondAcceptor { bond, atom } => Self::DativeBondAcceptor {
                bond: remap.dative_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            Self::DativeBondParallels { dative, parallel } => Self::DativeBondParallels {
                dative: remap.dative_bond(dative)?,
                parallel: remap.bond(parallel)?,
            },
            Self::DativeBondDonorSatisfies { bond, predicate } => Self::DativeBondDonorSatisfies {
                bond: remap.dative_bond(bond)?,
                predicate,
            },
            Self::DativeBondAcceptorSatisfies { bond, predicate } => {
                Self::DativeBondAcceptorSatisfies {
                    bond: remap.dative_bond(bond)?,
                    predicate,
                }
            }
            Self::AromaticSystemAtoms { system, atoms } => {
                let system = remap.aromatic_system(system)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::AromaticSystemAtoms {
                    system,
                    atoms: atoms?,
                }
            }
            Self::AromaticSystemContains { system, atom } => Self::AromaticSystemContains {
                system: remap.aromatic_system(system)?,
                atom: remap.atom(atom)?,
            },
            Self::AromaticSystemContainsAll { system, atoms } => {
                let system = remap.aromatic_system(system)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::AromaticSystemContainsAll {
                    system,
                    atoms: atoms?,
                }
            }
            Self::AromaticSystemAllAtoms { system, predicate } => Self::AromaticSystemAllAtoms {
                system: remap.aromatic_system(system)?,
                predicate,
            },
            Self::AromaticSystemAnyAtom { system, predicate } => Self::AromaticSystemAnyAtom {
                system: remap.aromatic_system(system)?,
                predicate,
            },
            Self::MulticenterBondAtoms { bond, atoms } => {
                let bond = remap.multicenter_bond(bond)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::MulticenterBondAtoms {
                    bond,
                    atoms: atoms?,
                }
            }
            Self::MulticenterBondContains { bond, atom } => Self::MulticenterBondContains {
                bond: remap.multicenter_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            Self::MulticenterBondContainsAll { bond, atoms } => {
                let bond = remap.multicenter_bond(bond)?;
                let atoms: Option<Vec<_>> = atoms.into_iter().map(|a| remap.atom(a)).collect();
                Self::MulticenterBondContainsAll {
                    bond,
                    atoms: atoms?,
                }
            }
            Self::MulticenterBondAllAtoms { bond, predicate } => Self::MulticenterBondAllAtoms {
                bond: remap.multicenter_bond(bond)?,
                predicate,
            },
            Self::MulticenterBondAnyAtom { bond, predicate } => Self::MulticenterBondAnyAtom {
                bond: remap.multicenter_bond(bond)?,
                predicate,
            },
            Self::NoncovalentBondEnds { bond, atoms } => {
                let bond = remap.noncovalent_bond(bond)?;
                let [a, b] = atoms;
                Self::NoncovalentBondEnds {
                    bond,
                    atoms: [remap.atom(a)?, remap.atom(b)?],
                }
            }
            Self::NoncovalentBondContains { bond, atom } => Self::NoncovalentBondContains {
                bond: remap.noncovalent_bond(bond)?,
                atom: remap.atom(atom)?,
            },
            Self::NoncovalentBondEndsSatisfy { bond, predicates } => {
                Self::NoncovalentBondEndsSatisfy {
                    bond: remap.noncovalent_bond(bond)?,
                    predicates,
                }
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::Remapping;

    use super::*;
    use crate::ast::value::ValueAst;

    fn remapping(
        removed_nodes: Vec<u32>,
        removed_edges: Vec<u32>,
        removed_dative: Vec<u32>,
        removed_aromatic: Vec<u32>,
        removed_multicenter: Vec<u32>,
        removed_noncovalent: Vec<u32>,
    ) -> IdRemapping {
        IdRemapping::new(
            Remapping::new(removed_nodes, removed_edges),
            removed_dative,
            removed_aromatic,
            removed_multicenter,
            removed_noncovalent,
        )
    }

    fn val_pred() -> Box<AtomConstraint> {
        Box::new(AtomConstraint::Valence(ValueAst::Lit(4)))
    }

    /// Drop atom 1; drop dative 0; preserve other entities. Indices above
    /// the removed slot shift down by one.
    fn one_atom_one_dative() -> IdRemapping {
        remapping(vec![1], vec![], vec![0], vec![], vec![], vec![])
    }

    /// Drop bond 0; preserve other entities.
    fn drop_bond0() -> IdRemapping {
        remapping(vec![], vec![0], vec![], vec![], vec![], vec![])
    }

    /// Drop aromatic system 0.
    fn drop_aromatic0() -> IdRemapping {
        remapping(vec![], vec![], vec![], vec![0], vec![], vec![])
    }

    /// Drop multicenter bond 0.
    fn drop_multicenter0() -> IdRemapping {
        remapping(vec![], vec![], vec![], vec![], vec![0], vec![])
    }

    /// Drop noncovalent bond 0.
    fn drop_noncovalent0() -> IdRemapping {
        remapping(vec![], vec![], vec![], vec![], vec![], vec![0])
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::dative_donor_shifts(RelationalConstraint::DativeBondDonor { bond: DativeBondId(1), atom: AtomId(2) },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondDonor { bond: DativeBondId(0), atom: AtomId(1) }))]
    #[case::dative_donor_drops_when_atom_removed(RelationalConstraint::DativeBondDonor { bond: DativeBondId(1), atom: AtomId(1) },
        one_atom_one_dative(), None)]
    #[case::dative_donor_drops_when_bond_removed(RelationalConstraint::DativeBondDonor { bond: DativeBondId(0), atom: AtomId(2) },
        one_atom_one_dative(), None)]
    #[case::dative_acceptor_shifts(RelationalConstraint::DativeBondAcceptor { bond: DativeBondId(1), atom: AtomId(2) },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondAcceptor { bond: DativeBondId(0), atom: AtomId(1) }))]
    #[case::dative_parallels_shifts(RelationalConstraint::DativeBondParallels { dative: DativeBondId(1), parallel: BondId(2) },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondParallels { dative: DativeBondId(0), parallel: BondId(2) }))]
    #[case::dative_parallels_drops_when_bond_removed(RelationalConstraint::DativeBondParallels { dative: DativeBondId(0), parallel: BondId(0) },
        drop_bond0(), None)]
    #[case::dative_donor_satisfies_shifts(RelationalConstraint::DativeBondDonorSatisfies { bond: DativeBondId(1), predicate: val_pred() },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondDonorSatisfies { bond: DativeBondId(0), predicate: val_pred() }))]
    #[case::dative_acceptor_satisfies_shifts(RelationalConstraint::DativeBondAcceptorSatisfies { bond: DativeBondId(1), predicate: val_pred() },
        one_atom_one_dative(), Some(RelationalConstraint::DativeBondAcceptorSatisfies { bond: DativeBondId(0), predicate: val_pred() }))]
    #[case::aromatic_atoms_shifts(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(1), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![0], vec![], vec![]), Some(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::aromatic_atoms_drops_when_atom_removed(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(0), atoms: vec![AtomId(1), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::aromatic_atoms_drops_when_system_removed(RelationalConstraint::AromaticSystemAtoms { system: AromaticSystemId(0), atoms: vec![AtomId(0)] },
        drop_aromatic0(), None)]
    #[case::aromatic_contains_shifts(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(1), atom: AtomId(2) },
        remapping(vec![1], vec![], vec![], vec![0], vec![], vec![]), Some(RelationalConstraint::AromaticSystemContains { system: AromaticSystemId(0), atom: AtomId(1) }))]
    #[case::aromatic_contains_all_shifts(RelationalConstraint::AromaticSystemContainsAll { system: AromaticSystemId(0), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![]), Some(RelationalConstraint::AromaticSystemContainsAll { system: AromaticSystemId(0),
        atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::aromatic_contains_all_drops_when_atom_removed(RelationalConstraint::AromaticSystemContainsAll { system: AromaticSystemId(0), atoms: vec![AtomId(1)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::aromatic_all_atoms_shifts(RelationalConstraint::AromaticSystemAllAtoms { system: AromaticSystemId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![0], vec![], vec![]),
        Some(RelationalConstraint::AromaticSystemAllAtoms { system: AromaticSystemId(0), predicate: val_pred() }))]
    #[case::aromatic_any_atom_shifts(RelationalConstraint::AromaticSystemAnyAtom { system: AromaticSystemId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![0], vec![], vec![]), Some(RelationalConstraint::AromaticSystemAnyAtom { system: AromaticSystemId(0), predicate: val_pred() }))]
    #[case::multicenter_atoms_shifts(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(1), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![0], vec![]), Some(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::multicenter_atoms_drops_when_atom_removed(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(0), atoms: vec![AtomId(1)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::multicenter_atoms_drops_when_bond_removed(RelationalConstraint::MulticenterBondAtoms { bond: MulticenterBondId(0), atoms: vec![AtomId(0)] },
        drop_multicenter0(), None)]
    #[case::multicenter_contains_shifts(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondId(1), atom: AtomId(2) },
        remapping(vec![1], vec![], vec![], vec![], vec![0], vec![]), Some(RelationalConstraint::MulticenterBondContains { bond: MulticenterBondId(0), atom: AtomId(1) }))]
    #[case::multicenter_contains_all_shifts(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondId(0), atoms: vec![AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![]), Some(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(1)] }))]
    #[case::multicenter_contains_all_drops_when_atom_removed(RelationalConstraint::MulticenterBondContainsAll { bond: MulticenterBondId(0), atoms: vec![AtomId(1)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::multicenter_all_atoms_shifts(RelationalConstraint::MulticenterBondAllAtoms { bond: MulticenterBondId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![], vec![0], vec![]),
        Some(RelationalConstraint::MulticenterBondAllAtoms { bond: MulticenterBondId(0), predicate: val_pred() }))]
    #[case::multicenter_any_atom_shifts(RelationalConstraint::MulticenterBondAnyAtom { bond: MulticenterBondId(1), predicate: val_pred() },
        remapping(vec![], vec![], vec![], vec![], vec![0], vec![]),
        Some(RelationalConstraint::MulticenterBondAnyAtom { bond: MulticenterBondId(0), predicate: val_pred() }))]
    #[case::noncovalent_ends_shifts(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(1), atoms: [AtomId(0), AtomId(2)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![0]), Some(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(0),
        atoms: [AtomId(0), AtomId(1)] }))]
    #[case::noncovalent_ends_drops_when_atom_removed(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(1)] },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![]), None)]
    #[case::noncovalent_ends_drops_when_bond_removed(RelationalConstraint::NoncovalentBondEnds { bond: NoncovalentBondId(0), atoms: [AtomId(0), AtomId(1)] },
        drop_noncovalent0(), None)]
    #[case::noncovalent_contains_shifts(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(1), atom: AtomId(2) },
        remapping(vec![1], vec![], vec![], vec![], vec![], vec![0]),
        Some(RelationalConstraint::NoncovalentBondContains { bond: NoncovalentBondId(0), atom: AtomId(1) }))]
    #[case::noncovalent_ends_satisfy_shifts(RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondId(1), predicates: [val_pred(), val_pred()] },
        remapping(vec![], vec![], vec![], vec![], vec![], vec![0]), Some(RelationalConstraint::NoncovalentBondEndsSatisfy { bond: NoncovalentBondId(0),
        predicates: [val_pred(), val_pred()] }))]
    fn test_relational_constraint_remap(
        #[case] input: RelationalConstraint,
        #[case] remap: IdRemapping,
        #[case] expected: Option<RelationalConstraint>,
    ) {
        assert_eq!(input.remap(&remap), expected);
    }
}
