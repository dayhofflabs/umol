//! Tier-2 (invariant) DPO validator: a reaction's deletions must be dangling-free, that is
//! deleting an atom must also delete every bond and overlay incident to it.
//! Operates on the permissive reaction AST. A reaction span establishes dangling-free projected
//! sides as a construction invariant.

use std::collections::HashSet;
use std::hash::Hash;

use thiserror::Error;
use umol_utils::solution::Solution;

use super::super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, DativeBondDelta, Delta, Deltas,
    MulticenterBondDelta, NoncovalentBondDelta,
};
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::molecule::MoleculeAst;

/// Checks the DPO dangling invariant: a deleted atom leaves no incident bond or overlay.
#[derive(Clone, Copy, Debug, Default)]
pub struct DpoValidator;

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

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum DpoError {}

impl DpoValidator {
    /// Over a `ReactionAst`: every lhs bond/overlay incident to a deleted atom must also be deleted.
    /// Over a reaction's LHS and deltas: every LHS bond or overlay incident to a deleted atom must
    /// also be deleted.
    pub fn validate_reaction(
        &self,
        lhs: &MoleculeAst,
        deltas: &Deltas,
    ) -> Result<Solution<(), DpoContradiction>, DpoError> {
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
                    return contradiction(DpoContradiction::DanglingBond { atom, bond });
                }
            }
            for dative in view.dative_bond_ids() {
                if !removed_dative.contains(&dative) {
                    return contradiction(DpoContradiction::DanglingDativeBond { atom, dative });
                }
            }
            if let Some(system) = view.aromatic_system_id() {
                if !removed_aromatic.contains(&system) {
                    return contradiction(DpoContradiction::DanglingAromaticSystem {
                        atom,
                        system,
                    });
                }
            }
            for multicenter in view.multicenter_bond_ids() {
                if !removed_multicenter.contains(&multicenter) {
                    return contradiction(DpoContradiction::DanglingMulticenterBond {
                        atom,
                        multicenter,
                    });
                }
            }
            for noncovalent in view.noncovalent_bond_ids() {
                if !removed_noncovalent.contains(&noncovalent) {
                    return contradiction(DpoContradiction::DanglingNoncovalentBond {
                        atom,
                        noncovalent,
                    });
                }
            }
        }
        Ok(Solution::Determined(()))
    }
}

fn removed<I: Eq + Hash>(deltas: &Deltas, extract: impl Fn(&Delta) -> Option<I>) -> HashSet<I> {
    deltas.iter().filter_map(extract).collect()
}

fn contradiction(c: DpoContradiction) -> Result<Solution<(), DpoContradiction>, DpoError> {
    Ok(Solution::Contradictory(c))
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use umol_chem::element::Element;

    use super::super::super::aromatic::AromaticSystemAst;
    use super::super::super::atom::AtomForm;
    use super::super::super::bond::BondForm;
    use super::super::super::constraint::Constraints;
    use super::super::super::dative::DativeBondAst;
    use super::super::super::molecule::{MoleculeAst, MoleculeEntries};
    use super::super::super::multicenter::MulticenterBondAst;
    use super::super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use super::super::super::reaction::ReactionAst;
    use super::*;

    #[rstest]
    #[case::co_deleted(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomForm::from_element(Element::C) }),
            Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: BondForm::from_order(1),
            }),
        ]),
    ))]
    #[case::no_deletion(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        Deltas::new(),
    ))]
    #[case::isolated_atom(ReactionAst::new(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            bonds: vec![],
            ..Default::default()
        }),
        Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
            id: AtomId(0),
            ast: AtomForm::from_element(Element::C),
        })]),
    ))]
    fn test_dpo_validator_validate_reaction(#[case] reaction: ReactionAst) {
        assert_eq!(
            DpoValidator
                .validate_reaction(&reaction.lhs, &reaction.deltas)
                .unwrap(),
            Solution::Determined(())
        );
    }

    #[rstest]
    #[case::bond(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::O)],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomForm::from_element(Element::C),
            })]),
        ),
        DpoContradiction::DanglingBond { atom: AtomId(0), bond: BondId(0) }
    )]
    #[case::dative(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::N), AtomForm::from_element(Element::B)],
                dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomForm::from_element(Element::N),
            })]),
        ),
        DpoContradiction::DanglingDativeBond { atom: AtomId(0), dative: DativeBondId(0) }
    )]
    #[case::aromatic(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::C)],
                aromatic: vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomForm::from_element(Element::C),
            })]),
        ),
        DpoContradiction::DanglingAromaticSystem { atom: AtomId(0), system: AromaticSystemId(0) }
    )]
    #[case::multicenter(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![
                    AtomForm::from_element(Element::B),
                    AtomForm::from_element(Element::H),
                    AtomForm::from_element(Element::B),
                ],
                multicenter: vec![(vec![AtomId(0), AtomId(1), AtomId(2)], MulticenterBondAst::from_electrons(vec![3, 5, 7]))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomForm::from_element(Element::B),
            })]),
        ),
        DpoContradiction::DanglingMulticenterBond { atom: AtomId(0), multicenter: MulticenterBondId(0) }
    )]
    #[case::noncovalent(
        ReactionAst::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::O), AtomForm::from_element(Element::O)],
                noncovalent: vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
                constraints: Constraints::new(),
                ..Default::default()
            }),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomForm::from_element(Element::O),
            })]),
        ),
        DpoContradiction::DanglingNoncovalentBond { atom: AtomId(0), noncovalent: NoncovalentBondId(0) }
    )]
    fn test_dpo_validator_validate_reaction_error(
        #[case] reaction: ReactionAst,
        #[case] expected: DpoContradiction,
    ) {
        assert_eq!(
            DpoValidator
                .validate_reaction(&reaction.lhs, &reaction.deltas)
                .unwrap(),
            Solution::Contradictory(expected)
        );
    }
}
