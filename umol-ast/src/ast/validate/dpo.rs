//! Tier-2 (invariant) DPO validator: a reaction's deletions must be dangling-free, that is
//! deleting an atom must also delete every bond and overlay incident to it.
//! Operates on both the reaction AST and the reaction span AST.

use std::collections::HashSet;
use std::hash::Hash;
use std::iter;

use thiserror::Error;
use umol_graph_core::{EdgeId, RelationId};
use umol_utils::solution::Solution;

use super::super::delta::{
    AromaticSystemDelta, AtomDelta, BondDelta, DativeBondDelta, Delta, Deltas, EntitySpan,
    MulticenterBondDelta, NoncovalentBondDelta,
};
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::reaction::ReactionAst;
use super::super::reaction_span::ReactionSpanAst;

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
    pub fn validate_reaction(
        &self,
        reaction: &ReactionAst,
    ) -> Result<Solution<(), DpoContradiction>, DpoError> {
        let removed_atoms = removed(&reaction.deltas, |d| match d {
            Delta::Atom(AtomDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        });
        let removed_bonds = removed(&reaction.deltas, |d| match d {
            Delta::Bond(BondDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        });
        let removed_dative = removed(&reaction.deltas, |d| match d {
            Delta::DativeBond(DativeBondDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        });
        let removed_aromatic = removed(&reaction.deltas, |d| match d {
            Delta::AromaticSystem(AromaticSystemDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        });
        let removed_multicenter = removed(&reaction.deltas, |d| match d {
            Delta::MulticenterBond(MulticenterBondDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        });
        let removed_noncovalent = removed(&reaction.deltas, |d| match d {
            Delta::NoncovalentBond(NoncovalentBondDelta::Remove { id, .. }) => Some(*id),
            _ => None,
        });

        for &atom in &removed_atoms {
            let view = reaction.lhs.atom(atom);
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

    /// Over a `ReactionSpanAst`: a `Removed` atom must carry no surviving (non-`Removed`) incidence.
    pub fn validate_reaction_span(
        &self,
        span: &ReactionSpanAst,
    ) -> Result<Solution<(), DpoContradiction>, DpoError> {
        let removed: HashSet<usize> = span
            .atoms()
            .iter()
            .enumerate()
            .filter(|(_, state)| matches!(state, EntitySpan::Removed(_)))
            .map(|(index, _)| index)
            .collect();

        for (edge, state) in span.bonds().iter().enumerate() {
            if matches!(state, EntitySpan::Removed(_)) {
                continue;
            }
            for node in span.graph().edge_endpoints(EdgeId(edge as u32)) {
                if removed.contains(&node.index()) {
                    return contradiction(DpoContradiction::DanglingBond {
                        atom: AtomId::from(node),
                        bond: BondId(edge as u32),
                    });
                }
            }
        }
        for i in 0..span.dative_bonds().relation_count() {
            let rid = RelationId(i as u32);
            if matches!(span.dative_bonds().data(rid), EntitySpan::Removed(_)) {
                continue;
            }
            let acceptor = span.dative_bonds().participants_1(rid)[0];
            let donors = span.dative_bonds().participants_2(rid).iter().copied();
            for node in iter::once(acceptor).chain(donors) {
                if removed.contains(&node.index()) {
                    return contradiction(DpoContradiction::DanglingDativeBond {
                        atom: AtomId::from(node),
                        dative: DativeBondId(i as u32),
                    });
                }
            }
        }
        for i in 0..span.aromatic_systems().relation_count() {
            let rid = RelationId(i as u32);
            if matches!(span.aromatic_systems().data(rid), EntitySpan::Removed(_)) {
                continue;
            }
            for &node in span.aromatic_systems().participants(rid) {
                if removed.contains(&node.index()) {
                    return contradiction(DpoContradiction::DanglingAromaticSystem {
                        atom: AtomId::from(node),
                        system: AromaticSystemId(i as u32),
                    });
                }
            }
        }
        for i in 0..span.multicenter_bonds().relation_count() {
            let rid = RelationId(i as u32);
            if matches!(span.multicenter_bonds().data(rid), EntitySpan::Removed(_)) {
                continue;
            }
            for &node in span.multicenter_bonds().participants(rid) {
                if removed.contains(&node.index()) {
                    return contradiction(DpoContradiction::DanglingMulticenterBond {
                        atom: AtomId::from(node),
                        multicenter: MulticenterBondId(i as u32),
                    });
                }
            }
        }
        for i in 0..span.noncovalent_bonds().relation_count() {
            let rid = RelationId(i as u32);
            if matches!(span.noncovalent_bonds().data(rid), EntitySpan::Removed(_)) {
                continue;
            }
            for &node in span.noncovalent_bonds().participants(rid) {
                if removed.contains(&node.index()) {
                    return contradiction(DpoContradiction::DanglingNoncovalentBond {
                        atom: AtomId::from(node),
                        noncovalent: NoncovalentBondId(i as u32),
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
    use super::super::super::atom::AtomAst;
    use super::super::super::bond::BondAst;
    use super::super::super::constraint::Constraints;
    use super::super::super::dative::DativeBondAst;
    use super::super::super::molecule::MoleculeAst;
    use super::super::super::multicenter::MulticenterBondAst;
    use super::super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use super::*;

    #[rstest]
    #[case::co_deleted(ReactionAst::new(
        MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
        Deltas::from_iter([
            Delta::Atom(AtomDelta::Remove { id: AtomId(0), ast: AtomAst::from_element(Element::C) }),
            Delta::Bond(BondDelta::Remove {
                id: BondId(0),
                atoms: [AtomId(0), AtomId(1)],
                ast: BondAst::from_order(1),
            }),
        ]),
    ))]
    #[case::no_deletion(ReactionAst::new(
        MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
        Deltas::new(),
    ))]
    #[case::isolated_atom(ReactionAst::new(
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C)], vec![]),
        Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
            id: AtomId(0),
            ast: AtomAst::from_element(Element::C),
        })]),
    ))]
    fn test_dpo_validator_validate_reaction(#[case] reaction: ReactionAst) {
        assert_eq!(
            DpoValidator.validate_reaction(&reaction).unwrap(),
            Solution::Determined(())
        );
    }

    #[rstest]
    #[case::bond(
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::O)],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            })]),
        ),
        DpoContradiction::DanglingBond { atom: AtomId(0), bond: BondId(0) }
    )]
    #[case::dative(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::N), AtomAst::from_element(Element::B)],
                vec![], vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
                vec![], vec![], vec![], vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::N),
            })]),
        ),
        DpoContradiction::DanglingDativeBond { atom: AtomId(0), dative: DativeBondId(0) }
    )]
    #[case::aromatic(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::C), AtomAst::from_element(Element::C)],
                vec![],
                vec![],
                vec![(vec![AtomId(0), AtomId(1)], AromaticSystemAst::from_electrons(vec![1, 2]))],
                vec![], vec![], vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            })]),
        ),
        DpoContradiction::DanglingAromaticSystem { atom: AtomId(0), system: AromaticSystemId(0) }
    )]
    #[case::multicenter(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![
                    AtomAst::from_element(Element::B),
                    AtomAst::from_element(Element::H),
                    AtomAst::from_element(Element::B),
                ],
                vec![], vec![], vec![],
                vec![(vec![AtomId(0), AtomId(1), AtomId(2)], MulticenterBondAst::from_electrons(vec![3, 5, 7]))],
                vec![], vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::B),
            })]),
        ),
        DpoContradiction::DanglingMulticenterBond { atom: AtomId(0), multicenter: MulticenterBondId(0) }
    )]
    #[case::noncovalent(
        ReactionAst::new(
            MoleculeAst::from_parts(
                vec![AtomAst::from_element(Element::O), AtomAst::from_element(Element::O)],
                vec![], vec![], vec![], vec![],
                vec![(AtomId(0), AtomId(1), NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))],
                vec![], vec![],
                Constraints::new(),
            ),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::O),
            })]),
        ),
        DpoContradiction::DanglingNoncovalentBond { atom: AtomId(0), noncovalent: NoncovalentBondId(0) }
    )]
    fn test_dpo_validator_validate_reaction_error(
        #[case] reaction: ReactionAst,
        #[case] expected: DpoContradiction,
    ) {
        assert_eq!(
            DpoValidator.validate_reaction(&reaction).unwrap(),
            Solution::Contradictory(expected)
        );
    }

    #[rstest]
    fn test_dpo_validator_validate_reaction_span() {
        let span = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(0),
                    ast: AtomAst::from_element(Element::C),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    atoms: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        )
        .to_reaction_span()
        .unwrap();
        assert_eq!(
            DpoValidator.validate_reaction_span(&span).unwrap(),
            Solution::Determined(())
        );
    }

    #[rstest]
    fn test_dpo_validator_validate_reaction_span_error() {
        // Delete the C but keep the C-O bond: the span has a `Removed` atom carrying an unchanged
        // bond.
        let span = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Atom(AtomDelta::Remove {
                id: AtomId(0),
                ast: AtomAst::from_element(Element::C),
            })]),
        )
        .to_reaction_span()
        .unwrap();
        assert_eq!(
            DpoValidator.validate_reaction_span(&span).unwrap(),
            Solution::Contradictory(DpoContradiction::DanglingBond {
                atom: AtomId(0),
                bond: BondId(0),
            })
        );
    }
}
