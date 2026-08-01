//! Model-independent constraints derived from entity fields and directly incident entities.

use thiserror::Error;
use umol_graph_core::ConnectedComponentsAlgorithm;
use umol_utils::solution::Solution;

use super::super::super::constraint::{
    AromaticSystemConstraintAst, AromaticValenceAst, AtomConstraintAst, BondConstraintAst,
    DativeBondConstraintAst, MulticenterBondConstraintAst, MulticenterValenceAst,
    NoncovalentBondConstraintAst,
};
use super::super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::super::molecule::MoleculeAst;
use super::super::super::stereo::{CisTransStereoAst, StereoKind, TetrahedralStereoAst};
use super::super::super::traits::Lattice;
use super::super::super::view::{
    AromaticSystemView, AtomView, BondView, DativeBondView, MulticenterBondView,
    NoncovalentBondView,
};
use super::ConstraintError;

/// Evaluates model-independent incidence constraints; only noncovalent `#I` requires a graph
/// algorithm selector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IncidenceConstraintValidator;

impl IncidenceConstraintValidator {
    /// Validate every inline incidence constraint in entity order.
    pub fn validate(
        &self,
        ast: &MoleculeAst,
        connected_components_algorithm: ConnectedComponentsAlgorithm,
    ) -> Result<Solution<(), IncidenceConstraintContradiction>, ConstraintError> {
        let mut bond_components = None;
        let mut any_underdetermined = false;

        for atom in ast.atoms().iter() {
            if let Some(contradiction) = observe(
                validate_atom_constraints(self, atom),
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for bond in ast.bonds().iter() {
            if let Some(contradiction) = observe(validate_bond(bond), &mut any_underdetermined) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for bond in ast.dative_bonds().iter() {
            if let Some(contradiction) =
                observe(validate_dative_bond(bond), &mut any_underdetermined)
            {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for system in ast.aromatic_systems().iter() {
            if let Some(contradiction) =
                observe(validate_aromatic_system(system), &mut any_underdetermined)
            {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for bond in ast.multicenter_bonds().iter() {
            if let Some(contradiction) =
                observe(validate_multicenter_bond(bond), &mut any_underdetermined)
            {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for bond in ast.noncovalent_bonds().iter() {
            let [a, b] = bond.atom_ids();
            let intramolecular = if bond
                .constraints()
                .iter()
                .any(|constraint| !constraint.is_undetermined())
            {
                let components = bond_components.get_or_insert_with(|| {
                    bond_components_by_atom(ast, connected_components_algorithm)
                });
                components[a.index()] == components[b.index()]
            } else {
                false
            };
            if let Some(contradiction) = observe(
                validate_noncovalent_bond(bond, intramolecular),
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }

        Ok(if any_underdetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }

    /// Validate one non-ring atom constraint against its incidence-derived value.
    /// Ring constraints are outside this validator and are determined identities here.
    pub fn validate_atom(
        &self,
        atom: AtomView<'_>,
        constraint: &AtomConstraintAst,
    ) -> Solution<(), IncidenceConstraintContradiction> {
        match constraint {
            AtomConstraintAst::Valence(_) => evaluate(
                constraint,
                &AtomConstraintAst::valence(atom.valence()),
                atom_contradiction(atom, constraint),
            ),
            AtomConstraintAst::DonatedPairs(_) => {
                // Multi-donor dative incidence has no defined per-atom projection pending the
                // coordination/haptic entity split in discussion doc 117.
                let unsupported = atom.dative_bonds().any(|bond| {
                    bond.donor_count() != 1 && bond.donor_ids().any(|donor| donor == atom.id)
                });
                if !constraint.is_undetermined() && unsupported {
                    Solution::Underdetermined(())
                } else {
                    evaluate(
                        constraint,
                        &AtomConstraintAst::donated_pairs(atom.donated_pairs()),
                        atom_contradiction(atom, constraint),
                    )
                }
            }
            AtomConstraintAst::AcceptedPairs(_) => {
                // Multi-donor dative incidence has no defined per-atom projection pending the
                // coordination/haptic entity split in discussion doc 117.
                let unsupported = atom
                    .dative_bonds()
                    .any(|bond| bond.donor_count() != 1 && bond.acceptor_id() == atom.id);
                if !constraint.is_undetermined() && unsupported {
                    Solution::Underdetermined(())
                } else {
                    evaluate(
                        constraint,
                        &AtomConstraintAst::accepted_pairs(atom.accepted_pairs()),
                        atom_contradiction(atom, constraint),
                    )
                }
            }
            AtomConstraintAst::AromaticValence(_) => {
                let derived = if atom.is_in_aromatic_system() {
                    AromaticValenceAst::aromatic(atom.aromatic_valence())
                } else {
                    AromaticValenceAst::NotAromatic
                };
                evaluate(
                    constraint,
                    &AtomConstraintAst::aromatic_valence(derived),
                    atom_contradiction(atom, constraint),
                )
            }
            AtomConstraintAst::MulticenterValence(_) => {
                let derived = if atom.is_in_multicenter_bond() {
                    MulticenterValenceAst::multicenter(atom.multicenter_valence())
                } else {
                    MulticenterValenceAst::NotMulticenter
                };
                evaluate(
                    constraint,
                    &AtomConstraintAst::multicenter_valence(derived),
                    atom_contradiction(atom, constraint),
                )
            }
            AtomConstraintAst::TetrahedralStereo(_) => {
                let derived = match atom.stereo_atom() {
                    Some(stereo) => match (
                        stereo.ast.configuration.kind(),
                        stereo.ast.configuration.coset(),
                    ) {
                        (Some(StereoKind::Tetrahedral), Some(coset)) => {
                            TetrahedralStereoAst::stereo(coset.clone())
                        }
                        (Some(_), _) => TetrahedralStereoAst::NotStereo,
                        (None, _) => TetrahedralStereoAst::Undetermined,
                    },
                    None => TetrahedralStereoAst::NotStereo,
                };
                evaluate(
                    constraint,
                    &AtomConstraintAst::tetrahedral_stereo(derived),
                    atom_contradiction(atom, constraint),
                )
            }
            AtomConstraintAst::Degree(_) => evaluate(
                constraint,
                &AtomConstraintAst::degree(atom.degree()),
                atom_contradiction(atom, constraint),
            ),
            AtomConstraintAst::TotalDegree(_) => evaluate(
                constraint,
                &AtomConstraintAst::total_degree(atom.total_degree()),
                atom_contradiction(atom, constraint),
            ),
            AtomConstraintAst::TotalValence(_) => evaluate(
                constraint,
                &AtomConstraintAst::total_valence(atom.total_valence()),
                atom_contradiction(atom, constraint),
            ),
            AtomConstraintAst::TotalHydrogens(_) => evaluate(
                constraint,
                &AtomConstraintAst::total_hydrogens(atom.total_hydrogens()),
                atom_contradiction(atom, constraint),
            ),
            AtomConstraintAst::RingDegree(_)
            | AtomConstraintAst::RingValence(_)
            | AtomConstraintAst::RingMembership(_) => Solution::Determined(()),
        }
    }
}

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum IncidenceConstraintContradiction {
    #[error("atom {atom:?} does not satisfy incidence constraint {constraint:?}")]
    Atom {
        atom: AtomId,
        constraint: AtomConstraintAst,
    },
    #[error("bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    Bond {
        bond: BondId,
        constraint: BondConstraintAst,
    },
    #[error("dative bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    DativeBond {
        bond: DativeBondId,
        constraint: DativeBondConstraintAst,
    },
    #[error("aromatic system {system:?} does not satisfy incidence constraint {constraint:?}")]
    AromaticSystem {
        system: AromaticSystemId,
        constraint: AromaticSystemConstraintAst,
    },
    #[error("multicenter bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    MulticenterBond {
        bond: MulticenterBondId,
        constraint: MulticenterBondConstraintAst,
    },
    #[error("noncovalent bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    NoncovalentBond {
        bond: NoncovalentBondId,
        constraint: NoncovalentBondConstraintAst,
    },
}

fn validate_atom_constraints(
    validator: &IncidenceConstraintValidator,
    atom: AtomView<'_>,
) -> Solution<(), IncidenceConstraintContradiction> {
    let mut any_underdetermined = false;
    for constraint in atom.constraints().iter() {
        let outcome = validator.validate_atom(atom, constraint);
        if let Some(contradiction) = observe(outcome, &mut any_underdetermined) {
            return Solution::Contradictory(contradiction);
        }
    }
    finish(any_underdetermined)
}

fn validate_bond(bond: BondView<'_>) -> Solution<(), IncidenceConstraintContradiction> {
    let mut any_underdetermined = false;
    for constraint in bond.constraints().iter() {
        let outcome = match constraint {
            BondConstraintAst::Aromatic(_) => evaluate(
                constraint,
                &BondConstraintAst::aromatic(bond.is_in_aromatic_system()),
                bond_contradiction(bond, constraint),
            ),
            BondConstraintAst::CisTransStereo(_) => {
                let derived = match bond.stereo_bond() {
                    Some(stereo) => match (
                        stereo.ast.configuration.kind(),
                        stereo.ast.configuration.coset(),
                    ) {
                        (Some(StereoKind::CisTrans), Some(coset)) => {
                            CisTransStereoAst::stereo(coset.clone())
                        }
                        (Some(_), _) => CisTransStereoAst::NotStereo,
                        (None, _) => CisTransStereoAst::Undetermined,
                    },
                    None => CisTransStereoAst::NotStereo,
                };
                evaluate(
                    constraint,
                    &BondConstraintAst::cis_trans_stereo(derived),
                    bond_contradiction(bond, constraint),
                )
            }
            BondConstraintAst::RingMembership(_) => continue,
        };
        if let Some(contradiction) = observe(outcome, &mut any_underdetermined) {
            return Solution::Contradictory(contradiction);
        }
    }
    finish(any_underdetermined)
}

fn validate_dative_bond(
    bond: DativeBondView<'_>,
) -> Solution<(), IncidenceConstraintContradiction> {
    let mut any_underdetermined = false;
    for constraint in bond.constraints().iter() {
        let outcome = match constraint {
            DativeBondConstraintAst::Aromatic(_) if bond.donor_count() != 1 => {
                // Aromatic incidence is defined only for a binary dative bond pending the
                // coordination/haptic entity split in discussion doc 117.
                if constraint.is_undetermined() {
                    Solution::Determined(())
                } else {
                    Solution::Underdetermined(())
                }
            }
            DativeBondConstraintAst::Aromatic(_) => {
                let donor_system = bond
                    .donors()
                    .next()
                    .and_then(|donor| donor.aromatic_system_id());
                let derived =
                    donor_system.is_some() && donor_system == bond.acceptor().aromatic_system_id();
                evaluate(
                    constraint,
                    &DativeBondConstraintAst::aromatic(derived),
                    IncidenceConstraintContradiction::DativeBond {
                        bond: bond.id,
                        constraint: constraint.clone(),
                    },
                )
            }
            DativeBondConstraintAst::RingMembership(_) => continue,
        };
        if let Some(contradiction) = observe(outcome, &mut any_underdetermined) {
            return Solution::Contradictory(contradiction);
        }
    }
    finish(any_underdetermined)
}

fn validate_aromatic_system(
    system: AromaticSystemView<'_>,
) -> Solution<(), IncidenceConstraintContradiction> {
    let mut any_underdetermined = false;
    for constraint in system.constraints().iter() {
        let derived = AromaticSystemConstraintAst::electron_count(system.electron_count());
        let contradiction = IncidenceConstraintContradiction::AromaticSystem {
            system: system.id,
            constraint: constraint.clone(),
        };
        if let Some(contradiction) = observe(
            evaluate(constraint, &derived, contradiction),
            &mut any_underdetermined,
        ) {
            return Solution::Contradictory(contradiction);
        }
    }
    finish(any_underdetermined)
}

fn validate_multicenter_bond(
    bond: MulticenterBondView<'_>,
) -> Solution<(), IncidenceConstraintContradiction> {
    let mut any_underdetermined = false;
    for constraint in bond.constraints().iter() {
        let derived = MulticenterBondConstraintAst::electron_count(bond.electron_count());
        let contradiction = IncidenceConstraintContradiction::MulticenterBond {
            bond: bond.id,
            constraint: constraint.clone(),
        };
        if let Some(contradiction) = observe(
            evaluate(constraint, &derived, contradiction),
            &mut any_underdetermined,
        ) {
            return Solution::Contradictory(contradiction);
        }
    }
    finish(any_underdetermined)
}

fn validate_noncovalent_bond(
    bond: NoncovalentBondView<'_>,
    intramolecular: bool,
) -> Solution<(), IncidenceConstraintContradiction> {
    let mut any_underdetermined = false;
    for constraint in bond.constraints().iter() {
        let derived = NoncovalentBondConstraintAst::intramolecular(intramolecular);
        let contradiction = IncidenceConstraintContradiction::NoncovalentBond {
            bond: bond.id,
            constraint: constraint.clone(),
        };
        if let Some(contradiction) = observe(
            evaluate(constraint, &derived, contradiction),
            &mut any_underdetermined,
        ) {
            return Solution::Contradictory(contradiction);
        }
    }
    finish(any_underdetermined)
}

fn evaluate<C>(
    asserted: &C,
    derived: &C,
    contradiction: IncidenceConstraintContradiction,
) -> Solution<(), IncidenceConstraintContradiction>
where
    C: Lattice,
{
    if asserted.is_undetermined() {
        Solution::Determined(())
    } else if !derived.is_ground() {
        Solution::Underdetermined(())
    } else if asserted.matches(derived) {
        Solution::Determined(())
    } else {
        Solution::Contradictory(contradiction)
    }
}

fn atom_contradiction(
    atom: AtomView<'_>,
    constraint: &AtomConstraintAst,
) -> IncidenceConstraintContradiction {
    IncidenceConstraintContradiction::Atom {
        atom: atom.id,
        constraint: constraint.clone(),
    }
}

fn bond_contradiction(
    bond: BondView<'_>,
    constraint: &BondConstraintAst,
) -> IncidenceConstraintContradiction {
    IncidenceConstraintContradiction::Bond {
        bond: bond.id,
        constraint: constraint.clone(),
    }
}

fn observe<C>(outcome: Solution<(), C>, any_underdetermined: &mut bool) -> Option<C> {
    match outcome {
        Solution::Determined(()) => None,
        Solution::Underdetermined(()) => {
            *any_underdetermined = true;
            None
        }
        Solution::Contradictory(contradiction) => Some(contradiction),
    }
}

fn finish<C>(any_underdetermined: bool) -> Solution<(), C> {
    if any_underdetermined {
        Solution::Underdetermined(())
    } else {
        Solution::Determined(())
    }
}

fn bond_components_by_atom(
    ast: &MoleculeAst,
    algorithm: ConnectedComponentsAlgorithm,
) -> Vec<usize> {
    let atom_count = ast.atoms().count();
    let mut component_by_atom = vec![0; atom_count];
    for (component, atoms) in ast
        .graph()
        .connected_components(algorithm)
        .into_iter()
        .enumerate()
    {
        for atom in atoms {
            component_by_atom[atom.index()] = component;
        }
    }
    component_by_atom
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::mol_dsl;

    #[rstest]
    #[case::valence(r#"{:atoms ["C#v1" "C"] :bonds [[0 1 "1"]]}"#)]
    #[case::dative_pairs(
        r#"{:atoms ["N#d1" "B#t1"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1"}]}"#
    )]
    #[case::aromatic_valence(r#"{:atoms ["C#a1" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :type "[1,1]"}]}"#)]
    #[case::not_aromatic(r#"{:atoms ["C#a!"] :bonds []}"#)]
    #[case::multicenter_valence(r#"{:atoms ["C#m1" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "[1,1,0]"}]}"#)]
    #[case::not_multicenter(r#"{:atoms ["C#m!"] :bonds []}"#)]
    #[case::tetrahedral_stereo(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"#)]
    #[case::not_tetrahedral_stereo(r#"{:atoms ["C#T!"] :bonds []}"#)]
    #[case::degree_totals(
        r#"{:atoms ["C#h1#v3#D2#X3#V4#H2" "H" "C"] :bonds [[0 1 "1"] [0 2 "2"]]}"#
    )]
    #[case::bond_aromatic(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]] :aromatic-systems [{:atoms [0 1] :type "[1,1]"}]}"#)]
    #[case::bond_not_aromatic(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a!"]]}"#)]
    #[case::bond_cis_trans(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"#)]
    #[case::bond_not_cis_trans(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#C!"]]}"#)]
    #[case::dative_aromatic(r#"{:atoms ["N" "B"] :bonds [[0 1 "1"]] :dative-bonds [{:donors [0] :acceptor 1 :type "1#a"}] :aromatic-systems [{:atoms [0 1] :type "[1,1]"}]}"#)]
    #[case::dative_not_aromatic(
        r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#a!"}]}"#
    )]
    #[case::aromatic_electrons(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :type "[1,1]#e2"}]}"#)]
    #[case::multicenter_electrons(r#"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "[1,1,0]#e2"}]}"#)]
    #[case::noncovalent_intramolecular(r#"{:atoms ["N" "H"] :bonds [[0 1 "1"]] :noncovalent-bonds [{:atoms [0 1] :type "Hbd#I"}]}"#)]
    #[case::noncovalent_intermolecular(
        r#"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd#I!"}]}"#
    )]
    #[case::finite_set(r#"{:atoms ["C#v{1,2}" "C"] :bonds [[0 1 "1"]]}"#)]
    #[case::range(r#"{:atoms ["C#v(1..)" "C"] :bonds [[0 1 "1"]]}"#)]
    #[case::vacuous(r#"{:atoms ["C#v*"] :bonds []}"#)]
    fn test_incidence_constraint_validator_validate_determined(#[case] input: &str) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            IncidenceConstraintValidator.validate(&molecule, ConnectedComponentsAlgorithm::Bfs,),
            Ok(Solution::Determined(()))
        );
    }

    #[rstest]
    #[case::atom(
        r#"{:atoms ["C#v2" "C"] :bonds [[0 1 "1"]]}"#,
        IncidenceConstraintContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintAst::valence(2),
        }
    )]
    #[case::bond(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#,
        IncidenceConstraintContradiction::Bond {
            bond: BondId(0),
            constraint: BondConstraintAst::aromatic(true),
        }
    )]
    #[case::dative(
        r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :type "1#a"}]}"#,
        IncidenceConstraintContradiction::DativeBond {
            bond: DativeBondId(0),
            constraint: DativeBondConstraintAst::aromatic(true),
        }
    )]
    #[case::aromatic_system(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :type "[1,1]#e3"}]}"#,
        IncidenceConstraintContradiction::AromaticSystem {
            system: AromaticSystemId(0),
            constraint: AromaticSystemConstraintAst::electron_count(3),
        }
    )]
    #[case::multicenter_bond(
        r#"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "[1,1,0]#e3"}]}"#,
        IncidenceConstraintContradiction::MulticenterBond {
            bond: MulticenterBondId(0),
            constraint: MulticenterBondConstraintAst::electron_count(3),
        }
    )]
    #[case::noncovalent_bond(
        r#"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :type "Hbd#I"}]}"#,
        IncidenceConstraintContradiction::NoncovalentBond {
            bond: NoncovalentBondId(0),
            constraint: NoncovalentBondConstraintAst::intramolecular(true),
        }
    )]
    fn test_incidence_constraint_validator_validate_contradictory(
        #[case] input: &str,
        #[case] expected: IncidenceConstraintContradiction,
    ) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            IncidenceConstraintValidator.validate(&molecule, ConnectedComponentsAlgorithm::Bfs,),
            Ok(Solution::Contradictory(expected))
        );
    }

    #[rstest]
    #[case::bond_order(r#"{:atoms ["C#v1" "C"] :bonds [[0 1 "*"]]}"#)]
    #[case::aromatic_valence(
        r#"{:atoms ["C#a1" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :type "*"}]}"#
    )]
    #[case::multicenter_valence(
        r#"{:atoms ["C#m1" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "*"}]}"#
    )]
    #[case::tetrahedral_coset(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th*"}]}"#)]
    #[case::cis_trans_coset(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct*"}]}"#)]
    #[case::multi_donor_donated(r#"{:atoms ["N#d1" "N" "B"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :type "1"}]}"#)]
    #[case::multi_donor_accepted(r#"{:atoms ["N" "N" "B#t1"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :type "1"}]}"#)]
    #[case::multi_donor_aromatic(r#"{:atoms ["N" "N" "B"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :type "1#a"}]}"#)]
    #[case::aromatic_electrons(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :type "*#e2"}]}"#
    )]
    #[case::multicenter_electrons(
        r#"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :type "*#e2"}]}"#
    )]
    fn test_incidence_constraint_validator_validate_underdetermined(#[case] input: &str) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            IncidenceConstraintValidator.validate(&molecule, ConnectedComponentsAlgorithm::Bfs,),
            Ok(Solution::Underdetermined(()))
        );
    }
}
