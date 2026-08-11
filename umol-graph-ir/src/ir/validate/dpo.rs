//! Reaction-application DPO check: a reaction's deletions must be dangling-free, that is
//! deleting an atom must also delete every bond and overlay incident to it.
//! Operates on the permissive `Reaction`. A reaction span establishes dangling-free projected
//! sides as a construction invariant.

use std::collections::HashSet;
use std::hash::Hash;

use thiserror::Error;

use super::super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta,
};
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::molecule::Molecule;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum DpoContradiction {
    #[error("deleted atom {atom:?} leaves dangling bond {bond:?}")]
    DanglingBond { atom: AtomId, bond: BondId },
    #[error("deleted atom {atom:?} leaves dangling dative bond {dative:?}")]
    DanglingDativeBond { atom: AtomId, dative: DativeBondId },
    #[error("deleted atom {atom:?} leaves dangling aromatic system {system:?}")]
    DanglingAromaticSystem {
        atom: AtomId,
        system: AromaticSystemId,
    },
    #[error("deleted atom {atom:?} leaves dangling multicenter bond {multicenter:?}")]
    DanglingMulticenterBond {
        atom: AtomId,
        multicenter: MulticenterBondId,
    },
    #[error("deleted atom {atom:?} leaves dangling noncovalent bond {noncovalent:?}")]
    DanglingNoncovalentBond {
        atom: AtomId,
        noncovalent: NoncovalentBondId,
    },
}

/// Check that every lhs bond or overlay incident to a deleted atom is also deleted.
pub fn check_reaction_dpo(lhs: &Molecule, deltas: &Deltas) -> Result<(), DpoContradiction> {
    let removed_atoms = removed(deltas, |d| match d {
        Delta::Atom(AtomDelta::Remove { id, .. }) => Some(*id),
        _ => None,
    });
    let removed_bonds = removed(deltas, |d| match d {
        Delta::Bond(BondDelta::Remove { id, .. }) => Some(*id),
        _ => None,
    });
    let removed_dative = removed(deltas, |d| match d {
        Delta::DativeBond(DativeBondDelta::Remove { id, .. }) => Some(*id),
        _ => None,
    });
    let removed_aromatic = removed(deltas, |d| match d {
        Delta::AromaticSystem(AromaticSystemDelta::Remove { id, .. }) => Some(*id),
        _ => None,
    });
    let removed_multicenter = removed(deltas, |d| match d {
        Delta::MulticenterBond(MulticenterBondDelta::Remove { id, .. }) => Some(*id),
        _ => None,
    });
    let removed_noncovalent = removed(deltas, |d| match d {
        Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, .. }) => Some(*id),
        _ => None,
    });

    for &atom in &removed_atoms {
        let view = lhs.atom(atom);
        for bond in view.bond_ids() {
            if !removed_bonds.contains(&bond) {
                return Err(DpoContradiction::DanglingBond { atom, bond });
            }
        }
        for dative in view.dative_bond_ids() {
            if !removed_dative.contains(&dative) {
                return Err(DpoContradiction::DanglingDativeBond { atom, dative });
            }
        }
        if let Some(system) = view.aromatic_system_id() {
            if !removed_aromatic.contains(&system) {
                return Err(DpoContradiction::DanglingAromaticSystem { atom, system });
            }
        }
        for multicenter in view.multicenter_bond_ids() {
            if !removed_multicenter.contains(&multicenter) {
                return Err(DpoContradiction::DanglingMulticenterBond { atom, multicenter });
            }
        }
        for noncovalent in view.noncovalent_bond_ids() {
            if !removed_noncovalent.contains(&noncovalent) {
                return Err(DpoContradiction::DanglingNoncovalentBond { atom, noncovalent });
            }
        }
    }
    Ok(())
}

fn removed<I: Eq + Hash>(deltas: &Deltas, extract: impl Fn(&Delta) -> Option<I>) -> HashSet<I> {
    deltas.iter().filter_map(extract).collect()
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use umol_chem::element::Element;

    use super::super::super::aromatic::AromaticSystemForm;
    use super::super::super::atom::AtomForm;
    use super::super::super::bond::BondForm;
    use super::super::super::constraint::Constraints;
    use super::super::super::dative::DativeBondForm;
    use super::super::super::molecule::{Molecule, MoleculeEntries};
    use super::super::super::multicenter::MulticenterBondForm;
    use super::super::super::noncovalent::{NoncovalentBondForm, NoncovalentBondKind};
    use super::super::super::reaction::Reaction;
    use super::*;

    #[rstest]
    #[case::co_deleted(Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C),
            }),
            Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                attributes: BondForm::from_order(1),
            }),
        ]),
    ))]
    #[case::no_deletion(Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        Deltas::new(),
    ))]
    #[case::isolated_atom(Reaction::new(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
            id: AtomId(0),
            attributes: AtomForm::from_element(Element::C),
        })]),
    ))]
    fn test_check_reaction_dpo(#[case] reaction: Reaction) {
        assert_eq!(check_reaction_dpo(&reaction.lhs, &reaction.deltas), Ok(()));
    }

    #[rstest]
    #[case::bond(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C),
            })]),
        ),
        DpoContradiction::DanglingBond { atom: AtomId(0), bond: BondId(0) }
    )]
    #[case::dative(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::N), AtomForm::from_element(Element::B)],
                dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::N),
            })]),
        ),
        DpoContradiction::DanglingDativeBond { atom: AtomId(0), dative: DativeBondId(0) }
    )]
    #[case::aromatic(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C)],
                aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemForm::from_electrons(vec![1, 2]))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::C),
            })]),
        ),
        DpoContradiction::DanglingAromaticSystem { atom: AtomId(0), system: AromaticSystemId(0) }
    )]
    #[case::multicenter(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::B),
                    AtomForm::from_element(Element::H),
                    AtomForm::from_element(Element::B),
                ],
                multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], MulticenterBondForm::from_electrons(vec![3, 5, 7]))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::B),
            })]),
        ),
        DpoContradiction::DanglingMulticenterBond { atom: AtomId(0), multicenter: MulticenterBondId(0) }
    )]
    #[case::noncovalent(
        Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
                noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                attributes: AtomForm::from_element(Element::O),
            })]),
        ),
        DpoContradiction::DanglingNoncovalentBond { atom: AtomId(0), noncovalent: NoncovalentBondId(0) }
    )]
    fn test_check_reaction_dpo_error(
        #[case] reaction: Reaction,
        #[case] expected: DpoContradiction,
    ) {
        assert_eq!(
            check_reaction_dpo(&reaction.lhs, &reaction.deltas),
            Err(expected)
        );
    }
}
