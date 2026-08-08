//! Molecule-scope constraints relating stored entities to their participants.

use std::collections::BTreeSet;

use thiserror::Error;
use umol_graph_core::RelevantCycleEnumerationAlgorithm;
use umol_utils::solution::Solution;

use super::super::super::constraint::{AtomConstraintAst, RelationalConstraint};
use super::super::super::entity::Entity;
use super::super::super::id::AtomId;
use super::super::super::molecule::MoleculeAst;
use super::super::super::ring::{RingConfig, RingModel};
use super::super::super::traits::Lattice;
use super::super::super::view::RingViews;
use super::incidence::validate_atom_constraint;
use super::ConstraintError;

/// Evaluates one molecule-scope relational constraint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RelationalConstraintValidator;

#[derive(Clone, Debug, Error, PartialEq, Eq)]
#[error("relational constraint is not satisfied: {constraint:?}")]
pub struct RelationalConstraintContradiction {
    pub constraint: RelationalConstraint,
}

impl RelationalConstraintValidator {
    pub fn validate(
        &self,
        ast: &MoleculeAst,
        constraint: &RelationalConstraint,
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Result<Solution<(), RelationalConstraintContradiction>, ConstraintError> {
        let rings = uses_ring_predicate(constraint).then(|| {
            ast.rings(
                RingModel::default(),
                RingConfig {
                    relevant_cycle_algorithm,
                    ..RingConfig::default()
                },
            )
        });
        let rings = rings.as_ref();

        let truth = match constraint {
            RelationalConstraint::DativeBondDonors { bond, atoms } => {
                let view = ast
                    .dative_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::DativeBond(*bond)))?;
                require_atoms(ast, atoms.iter().copied())?;
                Truth::from_bool(atom_set(view.donor_ids()) == atom_set(atoms.iter().copied()))
            }
            RelationalConstraint::DativeBondDonor { bond, atom } => {
                let view = ast
                    .dative_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::DativeBond(*bond)))?;
                require_entity(ast, Entity::Atom(*atom))?;
                Truth::from_bool(view.donor_ids().any(|candidate| candidate == *atom))
            }
            RelationalConstraint::DativeBondContainsAllDonors { bond, atoms } => {
                let view = ast
                    .dative_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::DativeBond(*bond)))?;
                require_atoms(ast, atoms.iter().copied())?;
                Truth::from_bool(
                    atom_set(view.donor_ids()).is_superset(&atom_set(atoms.iter().copied())),
                )
            }
            RelationalConstraint::DativeBondAllDonors { bond, predicate } => {
                let view = ast
                    .dative_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::DativeBond(*bond)))?;
                all_atoms(ast, view.donor_ids(), predicate, rings)
            }
            RelationalConstraint::DativeBondAnyDonor { bond, predicate } => {
                let view = ast
                    .dative_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::DativeBond(*bond)))?;
                any_atom(ast, view.donor_ids(), predicate, rings)
            }
            RelationalConstraint::DativeBondAcceptor { bond, atom } => {
                let view = ast
                    .dative_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::DativeBond(*bond)))?;
                require_entity(ast, Entity::Atom(*atom))?;
                Truth::from_bool(view.acceptor_id() == *atom)
            }
            RelationalConstraint::DativeBondAcceptorSatisfies { bond, predicate } => {
                let view = ast
                    .dative_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::DativeBond(*bond)))?;
                evaluate_atom(ast, view.acceptor_id(), predicate, rings)
            }
            RelationalConstraint::DativeBondParallels { dative, parallel } => {
                let dative_view = ast
                    .dative_bonds()
                    .get(*dative)
                    .ok_or_else(|| invalid_reference(Entity::DativeBond(*dative)))?;
                let parallel_view = ast
                    .bonds()
                    .get(*parallel)
                    .ok_or_else(|| invalid_reference(Entity::Bond(*parallel)))?;
                let parallel_atoms = unordered_pair(parallel_view.atom_ids());
                Truth::from_bool(dative_view.donor_ids().any(|donor| {
                    unordered_pair([donor, dative_view.acceptor_id()]) == parallel_atoms
                }))
            }
            RelationalConstraint::AromaticSystemAtoms { system, atoms } => {
                let view = ast
                    .aromatic_systems()
                    .get(*system)
                    .ok_or_else(|| invalid_reference(Entity::AromaticSystem(*system)))?;
                require_atoms(ast, atoms.iter().copied())?;
                Truth::from_bool(atom_set(view.atom_ids()) == atom_set(atoms.iter().copied()))
            }
            RelationalConstraint::AromaticSystemContains { system, atom } => {
                let view = ast
                    .aromatic_systems()
                    .get(*system)
                    .ok_or_else(|| invalid_reference(Entity::AromaticSystem(*system)))?;
                require_entity(ast, Entity::Atom(*atom))?;
                Truth::from_bool(view.atom_ids().any(|candidate| candidate == *atom))
            }
            RelationalConstraint::AromaticSystemContainsAll { system, atoms } => {
                let view = ast
                    .aromatic_systems()
                    .get(*system)
                    .ok_or_else(|| invalid_reference(Entity::AromaticSystem(*system)))?;
                require_atoms(ast, atoms.iter().copied())?;
                Truth::from_bool(
                    atom_set(view.atom_ids()).is_superset(&atom_set(atoms.iter().copied())),
                )
            }
            RelationalConstraint::AromaticSystemAllAtoms { system, predicate } => {
                let view = ast
                    .aromatic_systems()
                    .get(*system)
                    .ok_or_else(|| invalid_reference(Entity::AromaticSystem(*system)))?;
                all_atoms(ast, view.atom_ids(), predicate, rings)
            }
            RelationalConstraint::AromaticSystemAnyAtom { system, predicate } => {
                let view = ast
                    .aromatic_systems()
                    .get(*system)
                    .ok_or_else(|| invalid_reference(Entity::AromaticSystem(*system)))?;
                any_atom(ast, view.atom_ids(), predicate, rings)
            }
            RelationalConstraint::MulticenterBondAtoms { bond, atoms } => {
                let view = ast
                    .multicenter_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::MulticenterBond(*bond)))?;
                require_atoms(ast, atoms.iter().copied())?;
                Truth::from_bool(atom_set(view.atom_ids()) == atom_set(atoms.iter().copied()))
            }
            RelationalConstraint::MulticenterBondContains { bond, atom } => {
                let view = ast
                    .multicenter_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::MulticenterBond(*bond)))?;
                require_entity(ast, Entity::Atom(*atom))?;
                Truth::from_bool(view.atom_ids().any(|candidate| candidate == *atom))
            }
            RelationalConstraint::MulticenterBondContainsAll { bond, atoms } => {
                let view = ast
                    .multicenter_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::MulticenterBond(*bond)))?;
                require_atoms(ast, atoms.iter().copied())?;
                Truth::from_bool(
                    atom_set(view.atom_ids()).is_superset(&atom_set(atoms.iter().copied())),
                )
            }
            RelationalConstraint::MulticenterBondAllAtoms { bond, predicate } => {
                let view = ast
                    .multicenter_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::MulticenterBond(*bond)))?;
                all_atoms(ast, view.atom_ids(), predicate, rings)
            }
            RelationalConstraint::MulticenterBondAnyAtom { bond, predicate } => {
                let view = ast
                    .multicenter_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::MulticenterBond(*bond)))?;
                any_atom(ast, view.atom_ids(), predicate, rings)
            }
            RelationalConstraint::NoncovalentBondEnds { bond, atoms } => {
                let view = ast
                    .noncovalent_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::NoncovalentBond(*bond)))?;
                require_atoms(ast, *atoms)?;
                Truth::from_bool(unordered_pair(view.atom_ids()) == unordered_pair(*atoms))
            }
            RelationalConstraint::NoncovalentBondContains { bond, atom } => {
                let view = ast
                    .noncovalent_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::NoncovalentBond(*bond)))?;
                require_entity(ast, Entity::Atom(*atom))?;
                Truth::from_bool(view.atom_ids().contains(atom))
            }
            RelationalConstraint::NoncovalentBondEndsSatisfy { bond, predicates } => {
                let view = ast
                    .noncovalent_bonds()
                    .get(*bond)
                    .ok_or_else(|| invalid_reference(Entity::NoncovalentBond(*bond)))?;
                let [first, second] = view.atom_ids();
                all_truth([
                    evaluate_atom(ast, first, &predicates[0], rings),
                    evaluate_atom(ast, second, &predicates[1], rings),
                ])
            }
            RelationalConstraint::StereoAtomSite { stereo_atom, atom } => {
                let view = ast
                    .stereo_atoms()
                    .get(*stereo_atom)
                    .ok_or_else(|| invalid_reference(Entity::StereoAtom(*stereo_atom)))?;
                require_entity(ast, Entity::Atom(*atom))?;
                Truth::from_bool(view.site_id() == *atom)
            }
            RelationalConstraint::StereoAtomContains { stereo_atom, atom } => {
                let view = ast
                    .stereo_atoms()
                    .get(*stereo_atom)
                    .ok_or_else(|| invalid_reference(Entity::StereoAtom(*stereo_atom)))?;
                require_entity(ast, Entity::Atom(*atom))?;
                Truth::from_bool(view.atom_ligand_ids().any(|candidate| candidate == *atom))
            }
            RelationalConstraint::StereoAtomLigands { stereo_atom, atoms } => {
                let view = ast
                    .stereo_atoms()
                    .get(*stereo_atom)
                    .ok_or_else(|| invalid_reference(Entity::StereoAtom(*stereo_atom)))?;
                require_atoms(ast, atoms.iter().copied())?;
                Truth::from_bool(
                    atom_set(view.atom_ligand_ids()) == atom_set(atoms.iter().copied()),
                )
            }
            RelationalConstraint::StereoAtomAllLigands {
                stereo_atom,
                predicate,
            } => {
                let view = ast
                    .stereo_atoms()
                    .get(*stereo_atom)
                    .ok_or_else(|| invalid_reference(Entity::StereoAtom(*stereo_atom)))?;
                all_atoms(ast, view.atom_ligand_ids(), predicate, rings)
            }
            RelationalConstraint::StereoAtomAnyLigand {
                stereo_atom,
                predicate,
            } => {
                let view = ast
                    .stereo_atoms()
                    .get(*stereo_atom)
                    .ok_or_else(|| invalid_reference(Entity::StereoAtom(*stereo_atom)))?;
                any_atom(ast, view.atom_ligand_ids(), predicate, rings)
            }
            RelationalConstraint::StereoBondSite { stereo_bond, bond } => {
                let view = ast
                    .stereo_bonds()
                    .get(*stereo_bond)
                    .ok_or_else(|| invalid_reference(Entity::StereoBond(*stereo_bond)))?;
                require_entity(ast, Entity::Bond(*bond))?;
                Truth::from_bool(view.site_id() == *bond)
            }
            RelationalConstraint::StereoBondContains { stereo_bond, atom } => {
                let view = ast
                    .stereo_bonds()
                    .get(*stereo_bond)
                    .ok_or_else(|| invalid_reference(Entity::StereoBond(*stereo_bond)))?;
                require_entity(ast, Entity::Atom(*atom))?;
                Truth::from_bool(view.atom_ligand_ids().any(|candidate| candidate == *atom))
            }
            RelationalConstraint::StereoBondLigands { stereo_bond, atoms } => {
                let view = ast
                    .stereo_bonds()
                    .get(*stereo_bond)
                    .ok_or_else(|| invalid_reference(Entity::StereoBond(*stereo_bond)))?;
                require_atoms(ast, atoms.iter().copied())?;
                Truth::from_bool(
                    atom_set(view.atom_ligand_ids()) == atom_set(atoms.iter().copied()),
                )
            }
            RelationalConstraint::StereoBondAllLigands {
                stereo_bond,
                predicate,
            } => {
                let view = ast
                    .stereo_bonds()
                    .get(*stereo_bond)
                    .ok_or_else(|| invalid_reference(Entity::StereoBond(*stereo_bond)))?;
                all_atoms(ast, view.atom_ligand_ids(), predicate, rings)
            }
            RelationalConstraint::StereoBondAnyLigand {
                stereo_bond,
                predicate,
            } => {
                let view = ast
                    .stereo_bonds()
                    .get(*stereo_bond)
                    .ok_or_else(|| invalid_reference(Entity::StereoBond(*stereo_bond)))?;
                any_atom(ast, view.atom_ligand_ids(), predicate, rings)
            }
        };

        Ok(match truth {
            Truth::True => Solution::Determined(()),
            Truth::Underdetermined => Solution::Underdetermined(()),
            Truth::False => Solution::Contradictory(RelationalConstraintContradiction {
                constraint: constraint.clone(),
            }),
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Truth {
    True,
    False,
    Underdetermined,
}

impl Truth {
    fn from_bool(value: bool) -> Self {
        if value {
            Self::True
        } else {
            Self::False
        }
    }
}

fn evaluate_atom(
    ast: &MoleculeAst,
    atom_id: AtomId,
    predicate: &AtomConstraintAst,
    rings: Option<&RingViews<'_>>,
) -> Truth {
    if is_ring_predicate(predicate) {
        let Some(rings) = rings else {
            return Truth::Underdetermined;
        };
        let ring_atom = rings.atom(atom_id);
        let (asserted, derived) = match predicate {
            AtomConstraintAst::RingDegree(asserted) => (asserted, ring_atom.ring_degree()),
            AtomConstraintAst::RingValence(asserted) => (asserted, ring_atom.ring_valence()),
            AtomConstraintAst::RingMembership(membership) => (
                &membership.count,
                ring_atom.ring_membership(membership.scope),
            ),
            _ => unreachable!("ring predicate classified above"),
        };
        if asserted.is_undetermined() {
            Truth::True
        } else if !derived.is_ground() {
            Truth::Underdetermined
        } else if asserted.matches(&derived) {
            Truth::True
        } else {
            Truth::False
        }
    } else {
        truth_from_solution(validate_atom_constraint(ast, atom_id, predicate))
    }
}

fn truth_from_solution<C>(outcome: Solution<(), C>) -> Truth {
    match outcome {
        Solution::Determined(()) => Truth::True,
        Solution::Underdetermined(()) => Truth::Underdetermined,
        Solution::Contradictory(_) => Truth::False,
    }
}

fn all_atoms<'a>(
    ast: &'a MoleculeAst,
    atoms: impl IntoIterator<Item = AtomId>,
    predicate: &AtomConstraintAst,
    rings: Option<&RingViews<'a>>,
) -> Truth {
    all_truth(
        atoms
            .into_iter()
            .map(|atom| evaluate_atom(ast, atom, predicate, rings)),
    )
}

fn any_atom<'a>(
    ast: &'a MoleculeAst,
    atoms: impl IntoIterator<Item = AtomId>,
    predicate: &AtomConstraintAst,
    rings: Option<&RingViews<'a>>,
) -> Truth {
    let mut any_underdetermined = false;
    for truth in atoms
        .into_iter()
        .map(|atom| evaluate_atom(ast, atom, predicate, rings))
    {
        match truth {
            Truth::True => return Truth::True,
            Truth::False => {}
            Truth::Underdetermined => any_underdetermined = true,
        }
    }
    if any_underdetermined {
        Truth::Underdetermined
    } else {
        Truth::False
    }
}

fn all_truth(truths: impl IntoIterator<Item = Truth>) -> Truth {
    let mut any_underdetermined = false;
    for truth in truths {
        match truth {
            Truth::True => {}
            Truth::False => return Truth::False,
            Truth::Underdetermined => any_underdetermined = true,
        }
    }
    if any_underdetermined {
        Truth::Underdetermined
    } else {
        Truth::True
    }
}

fn uses_ring_predicate(constraint: &RelationalConstraint) -> bool {
    match constraint {
        RelationalConstraint::DativeBondAllDonors { predicate, .. }
        | RelationalConstraint::DativeBondAnyDonor { predicate, .. }
        | RelationalConstraint::DativeBondAcceptorSatisfies { predicate, .. }
        | RelationalConstraint::AromaticSystemAllAtoms { predicate, .. }
        | RelationalConstraint::AromaticSystemAnyAtom { predicate, .. }
        | RelationalConstraint::MulticenterBondAllAtoms { predicate, .. }
        | RelationalConstraint::MulticenterBondAnyAtom { predicate, .. }
        | RelationalConstraint::StereoAtomAllLigands { predicate, .. }
        | RelationalConstraint::StereoAtomAnyLigand { predicate, .. }
        | RelationalConstraint::StereoBondAllLigands { predicate, .. }
        | RelationalConstraint::StereoBondAnyLigand { predicate, .. } => {
            is_ring_predicate(predicate)
        }
        RelationalConstraint::NoncovalentBondEndsSatisfy { predicates, .. } => predicates
            .iter()
            .any(|predicate| is_ring_predicate(predicate)),
        _ => false,
    }
}

fn is_ring_predicate(predicate: &AtomConstraintAst) -> bool {
    matches!(
        predicate,
        AtomConstraintAst::RingDegree(_)
            | AtomConstraintAst::RingValence(_)
            | AtomConstraintAst::RingMembership(_)
    )
}

fn require_entity(ast: &MoleculeAst, entity: Entity) -> Result<(), ConstraintError> {
    let present = match entity {
        Entity::Atom(id) => ast.atoms().contains(id),
        Entity::Bond(id) => ast.bonds().contains(id),
        Entity::DativeBond(id) => ast.dative_bonds().contains(id),
        Entity::AromaticSystem(id) => ast.aromatic_systems().contains(id),
        Entity::MulticenterBond(id) => ast.multicenter_bonds().contains(id),
        Entity::NoncovalentBond(id) => ast.noncovalent_bonds().contains(id),
        Entity::StereoAtom(id) => ast.stereo_atoms().contains(id),
        Entity::StereoBond(id) => ast.stereo_bonds().contains(id),
    };
    if present {
        Ok(())
    } else {
        Err(invalid_reference(entity))
    }
}

fn require_atoms(
    ast: &MoleculeAst,
    atoms: impl IntoIterator<Item = AtomId>,
) -> Result<(), ConstraintError> {
    for atom in atoms {
        require_entity(ast, Entity::Atom(atom))?;
    }
    Ok(())
}

fn invalid_reference(entity: Entity) -> ConstraintError {
    ConstraintError::InvalidReference { entity }
}

fn atom_set(atoms: impl IntoIterator<Item = AtomId>) -> BTreeSet<AtomId> {
    atoms.into_iter().collect()
}

fn unordered_pair(mut atoms: [AtomId; 2]) -> [AtomId; 2] {
    atoms.sort_unstable();
    atoms
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};

    use super::*;
    use crate::ir::constraint::RingScope;
    use crate::ir::id::{
        AromaticSystemId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId, StereoAtomId,
        StereoBondId,
    };
    use crate::ir::value::ValueAst;
    use crate::mol_dsl;

    #[fixture]
    fn relational_molecule() -> MoleculeAst {
        mol_dsl!(
            r#"{:atoms ["C" "N" "O" "F" "Cl" "Br" "I" "H"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"] [3 4 "2"] [6 7 "*"]]
                :dative-bonds [{:donors [1 2] :acceptor 0 :type "1"}]
                :aromatic-systems [{:atoms [0 1 2] :type "[1,1,1]"}]
                :multicenter-bonds [{:atoms [0 2 6] :type "[1,1,0]"}]
                :noncovalent-bonds [{:atoms [1 4] :type "Hbd"}]
                :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]
                :stereo-bonds [{:site 3 :ligands [0 5] :type "Ct1"}]}"#
        )
    }

    #[rstest]
    #[case::dative_donors(RelationalConstraint::DativeBondDonors {
        bond: DativeBondId(0),
        atoms: vec![AtomId(2), AtomId(1)],
    })]
    #[case::dative_donor(RelationalConstraint::DativeBondDonor {
        bond: DativeBondId(0),
        atom: AtomId(1),
    })]
    #[case::dative_contains_all(RelationalConstraint::DativeBondContainsAllDonors {
        bond: DativeBondId(0),
        atoms: vec![AtomId(2)],
    })]
    #[case::dative_all(RelationalConstraint::DativeBondAllDonors {
        bond: DativeBondId(0),
        predicate: Box::new(AtomConstraintAst::degree(2)),
    })]
    #[case::dative_any(RelationalConstraint::DativeBondAnyDonor {
        bond: DativeBondId(0),
        predicate: Box::new(AtomConstraintAst::degree(2)),
    })]
    #[case::dative_acceptor(RelationalConstraint::DativeBondAcceptor {
        bond: DativeBondId(0),
        atom: AtomId(0),
    })]
    #[case::dative_acceptor_satisfies(RelationalConstraint::DativeBondAcceptorSatisfies {
        bond: DativeBondId(0),
        predicate: Box::new(AtomConstraintAst::degree(2)),
    })]
    #[case::dative_parallels(RelationalConstraint::DativeBondParallels {
        dative: DativeBondId(0),
        parallel: BondId(0),
    })]
    #[case::aromatic_atoms(RelationalConstraint::AromaticSystemAtoms {
        system: AromaticSystemId(0),
        atoms: vec![AtomId(2), AtomId(0), AtomId(1)],
    })]
    #[case::aromatic_contains(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(2),
    })]
    #[case::aromatic_contains_all(RelationalConstraint::AromaticSystemContainsAll {
        system: AromaticSystemId(0),
        atoms: vec![AtomId(2), AtomId(0)],
    })]
    #[case::aromatic_all(RelationalConstraint::AromaticSystemAllAtoms {
        system: AromaticSystemId(0),
        predicate: Box::new(AtomConstraintAst::degree(2)),
    })]
    #[case::aromatic_any_ring(RelationalConstraint::AromaticSystemAnyAtom {
        system: AromaticSystemId(0),
        predicate: Box::new(AtomConstraintAst::ring_membership(RingScope::Size(3), 1)),
    })]
    #[case::multicenter_atoms(RelationalConstraint::MulticenterBondAtoms {
        bond: MulticenterBondId(0),
        atoms: vec![AtomId(6), AtomId(2), AtomId(0)],
    })]
    #[case::multicenter_contains(RelationalConstraint::MulticenterBondContains {
        bond: MulticenterBondId(0),
        atom: AtomId(6),
    })]
    #[case::multicenter_contains_all(RelationalConstraint::MulticenterBondContainsAll {
        bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(6)],
    })]
    #[case::multicenter_all(RelationalConstraint::MulticenterBondAllAtoms {
        bond: MulticenterBondId(0),
        predicate: Box::new(AtomConstraintAst::degree(ValueAst::RangeFrom(1))),
    })]
    #[case::multicenter_any(RelationalConstraint::MulticenterBondAnyAtom {
        bond: MulticenterBondId(0),
        predicate: Box::new(AtomConstraintAst::degree(1)),
    })]
    #[case::noncovalent_ends(RelationalConstraint::NoncovalentBondEnds {
        bond: NoncovalentBondId(0),
        atoms: [AtomId(4), AtomId(1)],
    })]
    #[case::noncovalent_contains(RelationalConstraint::NoncovalentBondContains {
        bond: NoncovalentBondId(0),
        atom: AtomId(4),
    })]
    #[case::noncovalent_satisfies(RelationalConstraint::NoncovalentBondEndsSatisfy {
        bond: NoncovalentBondId(0),
        predicates: [
            Box::new(AtomConstraintAst::degree(2)),
            Box::new(AtomConstraintAst::degree(1)),
        ],
    })]
    #[case::stereo_atom_site(RelationalConstraint::StereoAtomSite {
        stereo_atom: StereoAtomId(0),
        atom: AtomId(0),
    })]
    #[case::stereo_atom_contains(RelationalConstraint::StereoAtomContains {
        stereo_atom: StereoAtomId(0),
        atom: AtomId(3),
    })]
    #[case::stereo_atom_ligands(RelationalConstraint::StereoAtomLigands {
        stereo_atom: StereoAtomId(0),
        atoms: vec![AtomId(4), AtomId(3), AtomId(2), AtomId(1)],
    })]
    #[case::stereo_atom_all(RelationalConstraint::StereoAtomAllLigands {
        stereo_atom: StereoAtomId(0),
        predicate: Box::new(AtomConstraintAst::degree(ValueAst::RangeFrom(1))),
    })]
    #[case::stereo_atom_any(RelationalConstraint::StereoAtomAnyLigand {
        stereo_atom: StereoAtomId(0),
        predicate: Box::new(AtomConstraintAst::degree(1)),
    })]
    #[case::stereo_bond_site(RelationalConstraint::StereoBondSite {
        stereo_bond: StereoBondId(0),
        bond: BondId(3),
    })]
    #[case::stereo_bond_contains(RelationalConstraint::StereoBondContains {
        stereo_bond: StereoBondId(0),
        atom: AtomId(5),
    })]
    #[case::stereo_bond_ligands(RelationalConstraint::StereoBondLigands {
        stereo_bond: StereoBondId(0),
        atoms: vec![AtomId(5), AtomId(0)],
    })]
    #[case::stereo_bond_all(RelationalConstraint::StereoBondAllLigands {
        stereo_bond: StereoBondId(0),
        predicate: Box::new(AtomConstraintAst::degree(ValueAst::RangeFrom(0))),
    })]
    #[case::stereo_bond_any(RelationalConstraint::StereoBondAnyLigand {
        stereo_bond: StereoBondId(0),
        predicate: Box::new(AtomConstraintAst::degree(0)),
    })]
    #[case::vacuous_predicate(RelationalConstraint::DativeBondAllDonors {
        bond: DativeBondId(0),
        predicate: Box::new(AtomConstraintAst::valence(ValueAst::Undetermined)),
    })]
    fn test_relational_constraint_validator_validate(
        relational_molecule: MoleculeAst,
        #[case] constraint: RelationalConstraint,
    ) {
        assert_eq!(
            RelationalConstraintValidator.validate(
                &relational_molecule,
                &constraint,
                RelevantCycleEnumerationAlgorithm::Vismara,
            ),
            Ok(Solution::Determined(())),
        );
    }

    #[rstest]
    #[case::all(RelationalConstraint::MulticenterBondAllAtoms {
        bond: MulticenterBondId(0),
        predicate: Box::new(AtomConstraintAst::valence(2)),
    })]
    #[case::any(RelationalConstraint::MulticenterBondAnyAtom {
        bond: MulticenterBondId(0),
        predicate: Box::new(AtomConstraintAst::valence(3)),
    })]
    fn test_relational_constraint_validator_validate_partial(
        relational_molecule: MoleculeAst,
        #[case] constraint: RelationalConstraint,
    ) {
        assert_eq!(
            RelationalConstraintValidator.validate(
                &relational_molecule,
                &constraint,
                RelevantCycleEnumerationAlgorithm::Vismara,
            ),
            Ok(Solution::Underdetermined(())),
        );
    }

    #[rstest]
    #[case::exact_set(RelationalConstraint::DativeBondDonors {
        bond: DativeBondId(0),
        atoms: vec![AtomId(1)],
    })]
    #[case::contains(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(3),
    })]
    #[case::contains_all(RelationalConstraint::MulticenterBondContainsAll {
        bond: MulticenterBondId(0),
        atoms: vec![AtomId(0), AtomId(3)],
    })]
    #[case::all(RelationalConstraint::AromaticSystemAllAtoms {
        system: AromaticSystemId(0),
        predicate: Box::new(AtomConstraintAst::degree(1)),
    })]
    #[case::any(RelationalConstraint::DativeBondAnyDonor {
        bond: DativeBondId(0),
        predicate: Box::new(AtomConstraintAst::degree(1)),
    })]
    #[case::acceptor(RelationalConstraint::DativeBondAcceptor {
        bond: DativeBondId(0),
        atom: AtomId(1),
    })]
    #[case::acceptor_predicate(RelationalConstraint::DativeBondAcceptorSatisfies {
        bond: DativeBondId(0),
        predicate: Box::new(AtomConstraintAst::degree(1)),
    })]
    #[case::parallel(RelationalConstraint::DativeBondParallels {
        dative: DativeBondId(0),
        parallel: BondId(3),
    })]
    #[case::endpoints(RelationalConstraint::NoncovalentBondEnds {
        bond: NoncovalentBondId(0),
        atoms: [AtomId(0), AtomId(1)],
    })]
    #[case::ordered_endpoints(RelationalConstraint::NoncovalentBondEndsSatisfy {
        bond: NoncovalentBondId(0),
        predicates: [
            Box::new(AtomConstraintAst::degree(1)),
            Box::new(AtomConstraintAst::degree(2)),
        ],
    })]
    #[case::atom_site(RelationalConstraint::StereoAtomSite {
        stereo_atom: StereoAtomId(0),
        atom: AtomId(1),
    })]
    #[case::bond_site(RelationalConstraint::StereoBondSite {
        stereo_bond: StereoBondId(0),
        bond: BondId(0),
    })]
    #[case::atom_ligands(RelationalConstraint::StereoAtomLigands {
        stereo_atom: StereoAtomId(0),
        atoms: vec![AtomId(1), AtomId(2)],
    })]
    #[case::bond_ligands(RelationalConstraint::StereoBondLigands {
        stereo_bond: StereoBondId(0),
        atoms: vec![AtomId(0)],
    })]
    fn test_relational_constraint_validator_validate_contradiction(
        relational_molecule: MoleculeAst,
        #[case] constraint: RelationalConstraint,
    ) {
        assert_eq!(
            RelationalConstraintValidator.validate(
                &relational_molecule,
                &constraint,
                RelevantCycleEnumerationAlgorithm::Vismara,
            ),
            Ok(Solution::Contradictory(RelationalConstraintContradiction {
                constraint: constraint.clone(),
            },)),
        );
    }

    #[rstest]
    #[case::atom(
        RelationalConstraint::DativeBondDonor {
            bond: DativeBondId(0),
            atom: AtomId(99),
        },
        Entity::Atom(AtomId(99)),
    )]
    #[case::bond(
        RelationalConstraint::DativeBondParallels {
            dative: DativeBondId(0),
            parallel: BondId(99),
        },
        Entity::Bond(BondId(99)),
    )]
    #[case::dative(
        RelationalConstraint::DativeBondDonor {
            bond: DativeBondId(99),
            atom: AtomId(0),
        },
        Entity::DativeBond(DativeBondId(99)),
    )]
    #[case::aromatic(
        RelationalConstraint::AromaticSystemContains {
            system: AromaticSystemId(99),
            atom: AtomId(0),
        },
        Entity::AromaticSystem(AromaticSystemId(99)),
    )]
    #[case::multicenter(
        RelationalConstraint::MulticenterBondContains {
            bond: MulticenterBondId(99),
            atom: AtomId(0),
        },
        Entity::MulticenterBond(MulticenterBondId(99)),
    )]
    #[case::noncovalent(
        RelationalConstraint::NoncovalentBondContains {
            bond: NoncovalentBondId(99),
            atom: AtomId(0),
        },
        Entity::NoncovalentBond(NoncovalentBondId(99)),
    )]
    #[case::stereo_atom(
        RelationalConstraint::StereoAtomContains {
            stereo_atom: StereoAtomId(99),
            atom: AtomId(0),
        },
        Entity::StereoAtom(StereoAtomId(99)),
    )]
    #[case::stereo_bond(
        RelationalConstraint::StereoBondContains {
            stereo_bond: StereoBondId(99),
            atom: AtomId(0),
        },
        Entity::StereoBond(StereoBondId(99)),
    )]
    fn test_relational_constraint_validator_validate_error(
        relational_molecule: MoleculeAst,
        #[case] constraint: RelationalConstraint,
        #[case] entity: Entity,
    ) {
        assert_eq!(
            RelationalConstraintValidator.validate(
                &relational_molecule,
                &constraint,
                RelevantCycleEnumerationAlgorithm::Vismara,
            ),
            Err(ConstraintError::InvalidReference { entity }),
        );
    }
}
