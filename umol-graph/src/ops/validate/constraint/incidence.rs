//! Model-independent constraints derived from entity fields and directly incident entities.

use thiserror::Error;
use umol_graph_core::ConnectedComponentsAlgorithm;
use umol_graph_ir::ir::{
    AromaticSystemConstraintForm, AromaticSystemId, AromaticValenceForm, AtomConstraintForm,
    AtomConstraintKey, AtomId, BondConstraintForm, BondConstraintKey, BondId, CisTransStereoForm,
    DativeBondConstraintForm, DativeBondId, Entity, Lattice, Molecule,
    MulticenterBondConstraintForm, MulticenterBondId, MulticenterValenceForm,
    NoncovalentBondConstraintForm, NoncovalentBondId, StereoKind, TetrahedralStereoForm,
};
use umol_utils::solution::Solution;

use super::ConstraintInvariantsError;

/// Evaluates model-independent incidence constraints; only noncovalent `#I` requires a graph
/// algorithm selector.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IncidenceConstraintInvariantsValidator;

impl IncidenceConstraintInvariantsValidator {
    /// Validate every inline incidence constraint in entity order.
    pub fn validate(
        &self,
        molecule: &Molecule,
        connected_components_algorithm: ConnectedComponentsAlgorithm,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let mut bond_components = None;
        let mut any_underdetermined = false;

        for id in molecule.atoms().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_atom(molecule, id)?,
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for id in molecule.bonds().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_bond(molecule, id)?,
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for id in molecule.dative_bonds().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_dative_bond(molecule, id)?,
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for id in molecule.aromatic_systems().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_aromatic_system(molecule, id)?,
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for id in molecule.multicenter_bonds().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_multicenter_bond(molecule, id)?,
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for id in molecule.noncovalent_bonds().ids() {
            let bond = molecule.noncovalent_bond(id);
            let [a, b] = bond.atom_ids();
            let intramolecular = if bond
                .constraints()
                .iter()
                .any(|constraint| !constraint.is_undetermined())
            {
                let components = bond_components.get_or_insert_with(|| {
                    bond_components_by_atom(molecule, connected_components_algorithm)
                });
                components[a.index()] == components[b.index()]
            } else {
                false
            };
            if let Some(contradiction) = observe(
                validate_noncovalent_bond(molecule, id, intramolecular),
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

    /// Validate all inline incidence constraints on one molecule atom.
    pub fn validate_molecule_atom(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let atom =
            molecule
                .atoms()
                .get(atom_id)
                .ok_or(ConstraintInvariantsError::InvalidReference {
                    entity: Entity::Atom(atom_id),
                })?;
        Ok(conjunction(atom.constraints().iter().map(|constraint| {
            validate_atom_constraint(molecule, atom_id, constraint)
        })))
    }

    /// Validate one inline atom constraint selected by its container key.
    pub fn validate_molecule_atom_constraint(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
        key: AtomConstraintKey,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let atom =
            molecule
                .atoms()
                .get(atom_id)
                .ok_or(ConstraintInvariantsError::InvalidReference {
                    entity: Entity::Atom(atom_id),
                })?;
        Ok(atom
            .constraints()
            .get(key)
            .map_or(Solution::Determined(()), |constraint| {
                validate_atom_constraint(molecule, atom_id, constraint)
            }))
    }

    /// Validate all inline incidence constraints on one molecule bond.
    pub fn validate_molecule_bond(
        &self,
        molecule: &Molecule,
        bond_id: BondId,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let bond =
            molecule
                .bonds()
                .get(bond_id)
                .ok_or(ConstraintInvariantsError::InvalidReference {
                    entity: Entity::Bond(bond_id),
                })?;
        Ok(conjunction(bond.constraints().iter().filter_map(
            |constraint| validate_bond_constraint(molecule, bond_id, constraint),
        )))
    }

    /// Validate one inline localized-bond constraint selected by its container key.
    pub fn validate_molecule_bond_constraint(
        &self,
        molecule: &Molecule,
        bond_id: BondId,
        key: BondConstraintKey,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let bond =
            molecule
                .bonds()
                .get(bond_id)
                .ok_or(ConstraintInvariantsError::InvalidReference {
                    entity: Entity::Bond(bond_id),
                })?;
        Ok(bond
            .constraints()
            .get(key)
            .and_then(|constraint| validate_bond_constraint(molecule, bond_id, constraint))
            .unwrap_or(Solution::Determined(())))
    }

    /// Validate all inline incidence constraints on one molecule dative bond.
    pub fn validate_molecule_dative_bond(
        &self,
        molecule: &Molecule,
        bond_id: DativeBondId,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let bond = molecule.dative_bonds().get(bond_id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::DativeBond(bond_id),
            },
        )?;
        Ok(conjunction(bond.constraints().iter().filter_map(
            |constraint| validate_dative_bond_constraint(molecule, bond_id, constraint),
        )))
    }

    /// Validate all inline incidence constraints on one molecule aromatic system.
    pub fn validate_molecule_aromatic_system(
        &self,
        molecule: &Molecule,
        system_id: AromaticSystemId,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let system = molecule.aromatic_systems().get(system_id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::AromaticSystem(system_id),
            },
        )?;
        Ok(conjunction(system.constraints().iter().map(|constraint| {
            validate_aromatic_system_constraint(molecule, system_id, constraint)
        })))
    }

    /// Validate all inline incidence constraints on one molecule multicenter bond.
    pub fn validate_molecule_multicenter_bond(
        &self,
        molecule: &Molecule,
        bond_id: MulticenterBondId,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let bond = molecule.multicenter_bonds().get(bond_id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::MulticenterBond(bond_id),
            },
        )?;
        Ok(conjunction(bond.constraints().iter().map(|constraint| {
            validate_multicenter_bond_constraint(molecule, bond_id, constraint)
        })))
    }

    /// Validate all inline incidence constraints on one molecule noncovalent bond.
    pub fn validate_molecule_noncovalent_bond(
        &self,
        molecule: &Molecule,
        bond_id: NoncovalentBondId,
        connected_components_algorithm: ConnectedComponentsAlgorithm,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let bond = molecule.noncovalent_bonds().get(bond_id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::NoncovalentBond(bond_id),
            },
        )?;
        let intramolecular = if bond
            .constraints()
            .iter()
            .any(|constraint| !constraint.is_undetermined())
        {
            let components = bond_components_by_atom(molecule, connected_components_algorithm);
            let [a, b] = bond.atom_ids();
            components[a.index()] == components[b.index()]
        } else {
            false
        };
        Ok(validate_noncovalent_bond(molecule, bond_id, intramolecular))
    }
}

/// Validate one non-ring atom constraint against its incidence-derived value.
/// Ring constraints are determined identities here.
pub fn validate_atom_constraint(
    molecule: &Molecule,
    atom_id: AtomId,
    constraint: &AtomConstraintForm,
) -> Solution<(), IncidenceConstraintInvariantsContradiction> {
    let atom = molecule.atom(atom_id);
    match constraint {
        AtomConstraintForm::Valence(_) => evaluate(
            constraint,
            &AtomConstraintForm::valence(atom.valence()),
            atom_contradiction(atom_id, constraint),
        ),
        AtomConstraintForm::DonatedPairs(_) => {
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
                    &AtomConstraintForm::donated_pairs(atom.donated_pairs()),
                    atom_contradiction(atom_id, constraint),
                )
            }
        }
        AtomConstraintForm::AcceptedPairs(_) => {
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
                    &AtomConstraintForm::accepted_pairs(atom.accepted_pairs()),
                    atom_contradiction(atom_id, constraint),
                )
            }
        }
        AtomConstraintForm::AromaticValence(_) => {
            let derived = if atom.is_in_aromatic_system() {
                AromaticValenceForm::aromatic(atom.aromatic_valence())
            } else {
                AromaticValenceForm::NotAromatic
            };
            evaluate(
                constraint,
                &AtomConstraintForm::aromatic_valence(derived),
                atom_contradiction(atom_id, constraint),
            )
        }
        AtomConstraintForm::MulticenterValence(_) => {
            let derived = if atom.is_in_multicenter_bond() {
                MulticenterValenceForm::multicenter(atom.multicenter_valence())
            } else {
                MulticenterValenceForm::NotMulticenter
            };
            evaluate(
                constraint,
                &AtomConstraintForm::multicenter_valence(derived),
                atom_contradiction(atom_id, constraint),
            )
        }
        AtomConstraintForm::TetrahedralStereo(_) => {
            let derived = match atom.stereo_atom() {
                Some(stereo) => match (
                    stereo.attributes.configuration.kind(),
                    stereo.attributes.configuration.coset(),
                ) {
                    (Some(StereoKind::Tetrahedral), Some(coset)) => {
                        TetrahedralStereoForm::stereo(coset.clone())
                    }
                    (Some(_), _) => TetrahedralStereoForm::NotStereo,
                    (None, _) => TetrahedralStereoForm::Undetermined,
                },
                None => TetrahedralStereoForm::NotStereo,
            };
            evaluate(
                constraint,
                &AtomConstraintForm::tetrahedral_stereo(derived),
                atom_contradiction(atom_id, constraint),
            )
        }
        AtomConstraintForm::Degree(_) => evaluate(
            constraint,
            &AtomConstraintForm::degree(atom.degree()),
            atom_contradiction(atom_id, constraint),
        ),
        AtomConstraintForm::TotalDegree(_) => evaluate(
            constraint,
            &AtomConstraintForm::total_degree(atom.total_degree()),
            atom_contradiction(atom_id, constraint),
        ),
        AtomConstraintForm::TotalValence(_) => evaluate(
            constraint,
            &AtomConstraintForm::total_valence(atom.total_valence()),
            atom_contradiction(atom_id, constraint),
        ),
        AtomConstraintForm::TotalHydrogens(_) => evaluate(
            constraint,
            &AtomConstraintForm::total_hydrogens(atom.total_hydrogens()),
            atom_contradiction(atom_id, constraint),
        ),
        AtomConstraintForm::RingDegree(_)
        | AtomConstraintForm::RingValence(_)
        | AtomConstraintForm::RingMembership(_) => Solution::Determined(()),
    }
}

#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum IncidenceConstraintInvariantsContradiction {
    #[error("atom {atom:?} does not satisfy incidence constraint {constraint:?}")]
    Atom {
        atom: AtomId,
        constraint: AtomConstraintForm,
    },
    #[error("bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    Bond {
        bond: BondId,
        constraint: BondConstraintForm,
    },
    #[error("dative bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    DativeBond {
        bond: DativeBondId,
        constraint: DativeBondConstraintForm,
    },
    #[error("aromatic system {system:?} does not satisfy incidence constraint {constraint:?}")]
    AromaticSystem {
        system: AromaticSystemId,
        constraint: AromaticSystemConstraintForm,
    },
    #[error("multicenter bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    MulticenterBond {
        bond: MulticenterBondId,
        constraint: MulticenterBondConstraintForm,
    },
    #[error("noncovalent bond {bond:?} does not satisfy incidence constraint {constraint:?}")]
    NoncovalentBond {
        bond: NoncovalentBondId,
        constraint: NoncovalentBondConstraintForm,
    },
}

pub fn validate_bond_constraint(
    molecule: &Molecule,
    bond_id: BondId,
    constraint: &BondConstraintForm,
) -> Option<Solution<(), IncidenceConstraintInvariantsContradiction>> {
    let bond = molecule.bond(bond_id);
    Some(match constraint {
        BondConstraintForm::Aromatic(_) => evaluate(
            constraint,
            &BondConstraintForm::aromatic(bond.is_in_aromatic_system()),
            bond_contradiction(bond_id, constraint),
        ),
        BondConstraintForm::CisTransStereo(_) => {
            let derived = match bond.stereo_bond() {
                Some(stereo) => match (
                    stereo.attributes.configuration.kind(),
                    stereo.attributes.configuration.coset(),
                ) {
                    (Some(StereoKind::CisTrans), Some(coset)) => {
                        CisTransStereoForm::stereo(coset.clone())
                    }
                    (Some(_), _) => CisTransStereoForm::NotStereo,
                    (None, _) => CisTransStereoForm::Undetermined,
                },
                None => CisTransStereoForm::NotStereo,
            };
            evaluate(
                constraint,
                &BondConstraintForm::cis_trans_stereo(derived),
                bond_contradiction(bond_id, constraint),
            )
        }
        BondConstraintForm::RingMembership(_) => return None,
    })
}

pub fn validate_dative_bond_constraint(
    molecule: &Molecule,
    bond_id: DativeBondId,
    constraint: &DativeBondConstraintForm,
) -> Option<Solution<(), IncidenceConstraintInvariantsContradiction>> {
    let bond = molecule.dative_bond(bond_id);
    Some(match constraint {
        DativeBondConstraintForm::Aromatic(_) if bond.donor_count() != 1 => {
            // Aromatic incidence is defined only for a binary dative bond pending the
            // coordination/haptic entity split in discussion doc 117.
            if constraint.is_undetermined() {
                Solution::Determined(())
            } else {
                Solution::Underdetermined(())
            }
        }
        DativeBondConstraintForm::Aromatic(_) => {
            let donor_system = bond
                .donors()
                .next()
                .and_then(|donor| donor.aromatic_system_id());
            let derived =
                donor_system.is_some() && donor_system == bond.acceptor().aromatic_system_id();
            evaluate(
                constraint,
                &DativeBondConstraintForm::aromatic(derived),
                IncidenceConstraintInvariantsContradiction::DativeBond {
                    bond: bond.id,
                    constraint: constraint.clone(),
                },
            )
        }
        DativeBondConstraintForm::RingMembership(_) => return None,
    })
}

pub fn validate_aromatic_system_constraint(
    molecule: &Molecule,
    system_id: AromaticSystemId,
    constraint: &AromaticSystemConstraintForm,
) -> Solution<(), IncidenceConstraintInvariantsContradiction> {
    let derived = AromaticSystemConstraintForm::electron_count(
        molecule.aromatic_system(system_id).electron_count(),
    );
    evaluate(
        constraint,
        &derived,
        IncidenceConstraintInvariantsContradiction::AromaticSystem {
            system: system_id,
            constraint: constraint.clone(),
        },
    )
}

pub fn validate_multicenter_bond_constraint(
    molecule: &Molecule,
    bond_id: MulticenterBondId,
    constraint: &MulticenterBondConstraintForm,
) -> Solution<(), IncidenceConstraintInvariantsContradiction> {
    let derived = MulticenterBondConstraintForm::electron_count(
        molecule.multicenter_bond(bond_id).electron_count(),
    );
    evaluate(
        constraint,
        &derived,
        IncidenceConstraintInvariantsContradiction::MulticenterBond {
            bond: bond_id,
            constraint: constraint.clone(),
        },
    )
}

pub fn validate_noncovalent_bond_constraint(
    bond_id: NoncovalentBondId,
    constraint: &NoncovalentBondConstraintForm,
    intramolecular: bool,
) -> Solution<(), IncidenceConstraintInvariantsContradiction> {
    let derived = NoncovalentBondConstraintForm::intramolecular(intramolecular);
    evaluate(
        constraint,
        &derived,
        IncidenceConstraintInvariantsContradiction::NoncovalentBond {
            bond: bond_id,
            constraint: constraint.clone(),
        },
    )
}

fn validate_noncovalent_bond(
    molecule: &Molecule,
    bond_id: NoncovalentBondId,
    intramolecular: bool,
) -> Solution<(), IncidenceConstraintInvariantsContradiction> {
    let bond = molecule.noncovalent_bond(bond_id);
    let mut any_underdetermined = false;
    for constraint in bond.constraints().iter() {
        if let Some(contradiction) = observe(
            validate_noncovalent_bond_constraint(bond_id, constraint, intramolecular),
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
    contradiction: IncidenceConstraintInvariantsContradiction,
) -> Solution<(), IncidenceConstraintInvariantsContradiction>
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
    atom: AtomId,
    constraint: &AtomConstraintForm,
) -> IncidenceConstraintInvariantsContradiction {
    IncidenceConstraintInvariantsContradiction::Atom {
        atom,
        constraint: constraint.clone(),
    }
}

fn bond_contradiction(
    bond: BondId,
    constraint: &BondConstraintForm,
) -> IncidenceConstraintInvariantsContradiction {
    IncidenceConstraintInvariantsContradiction::Bond {
        bond,
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

fn conjunction<C>(outcomes: impl IntoIterator<Item = Solution<(), C>>) -> Solution<(), C> {
    let mut any_underdetermined = false;
    for outcome in outcomes {
        if let Some(contradiction) = observe(outcome, &mut any_underdetermined) {
            return Solution::Contradictory(contradiction);
        }
    }
    finish(any_underdetermined)
}

pub fn bond_components_by_atom(
    molecule: &Molecule,
    algorithm: ConnectedComponentsAlgorithm,
) -> Vec<usize> {
    let atom_count = molecule.atoms().count();
    let mut component_by_atom = vec![0; atom_count];
    for (component, atoms) in molecule
        .graph()
        .enumerate_connected_components(algorithm)
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
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    #[case::determined(
        mol_dsl!(r#"{:atoms ["C#v0"]}"#),
        AtomId(0),
        AtomConstraintKey::Valence,
        Ok(Solution::Determined(())),
    )]
    #[case::absent(
        mol_dsl!(r#"{:atoms ["C"]}"#),
        AtomId(0),
        AtomConstraintKey::Valence,
        Ok(Solution::Determined(())),
    )]
    #[case::contradictory(
        mol_dsl!(r#"{:atoms ["C#v1"]}"#),
        AtomId(0),
        AtomConstraintKey::Valence,
        Ok(Solution::Contradictory(IncidenceConstraintInvariantsContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintForm::valence(1),
        })),
    )]
    #[case::underdetermined(
        mol_dsl!(r#"{:atoms ["C#v1" "C"] :bonds [[0 1 "*"]]}"#),
        AtomId(0),
        AtomConstraintKey::Valence,
        Ok(Solution::Underdetermined(())),
    )]
    #[case::invalid_reference(
        mol_dsl!(r#"{:atoms ["C"]}"#),
        AtomId(1),
        AtomConstraintKey::Valence,
        Err(ConstraintInvariantsError::InvalidReference { entity: Entity::Atom(AtomId(1)) }),
    )]
    fn test_incidence_constraint_validator_validate_molecule_atom_constraint(
        #[case] molecule: Molecule,
        #[case] atom: AtomId,
        #[case] key: AtomConstraintKey,
        #[case] expected: Result<
            Solution<(), IncidenceConstraintInvariantsContradiction>,
            ConstraintInvariantsError,
        >,
    ) {
        assert_eq!(
            IncidenceConstraintInvariantsValidator
                .validate_molecule_atom_constraint(&molecule, atom, key),
            expected
        );
    }

    #[rstest]
    #[case::determined(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a!"]]}"#),
        BondId(0),
        BondConstraintKey::Aromatic,
        Ok(Solution::Determined(())),
    )]
    #[case::absent(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
        BondId(0),
        BondConstraintKey::Aromatic,
        Ok(Solution::Determined(())),
    )]
    #[case::contradictory(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#),
        BondId(0),
        BondConstraintKey::Aromatic,
        Ok(Solution::Contradictory(IncidenceConstraintInvariantsContradiction::Bond {
            bond: BondId(0),
            constraint: BondConstraintForm::aromatic(true),
        })),
    )]
    #[case::invalid_reference(
        mol_dsl!(r#"{:atoms ["C"]}"#),
        BondId(0),
        BondConstraintKey::Aromatic,
        Err(ConstraintInvariantsError::InvalidReference { entity: Entity::Bond(BondId(0)) }),
    )]
    fn test_incidence_constraint_validator_validate_molecule_bond_constraint(
        #[case] molecule: Molecule,
        #[case] bond: BondId,
        #[case] key: BondConstraintKey,
        #[case] expected: Result<
            Solution<(), IncidenceConstraintInvariantsContradiction>,
            ConstraintInvariantsError,
        >,
    ) {
        assert_eq!(
            IncidenceConstraintInvariantsValidator
                .validate_molecule_bond_constraint(&molecule, bond, key),
            expected
        );
    }

    #[rstest]
    #[case::valence(r#"{:atoms ["C#v1" "C"] :bonds [[0 1 "1"]]}"#)]
    #[case::dative_pairs(
        r#"{:atoms ["N#d1" "B#t1"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#
    )]
    #[case::aromatic_valence(r#"{:atoms ["C#a1" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]"}]}"#)]
    #[case::not_aromatic(r#"{:atoms ["C#a!"] :bonds []}"#)]
    #[case::multicenter_valence(r#"{:atoms ["C#m1" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,1,0]"}]}"#)]
    #[case::not_multicenter(r#"{:atoms ["C#m!"] :bonds []}"#)]
    #[case::tetrahedral_stereo(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#)]
    #[case::not_tetrahedral_stereo(r#"{:atoms ["C#T!"] :bonds []}"#)]
    #[case::degree_totals(
        r#"{:atoms ["C#h1#v3#D2#X3#V4#H2" "H" "C"] :bonds [[0 1 "1"] [0 2 "2"]]}"#
    )]
    #[case::bond_aromatic(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]"}]}"#)]
    #[case::bond_not_aromatic(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a!"]]}"#)]
    #[case::bond_cis_trans(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#)]
    #[case::bond_not_cis_trans(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#C!"]]}"#)]
    #[case::dative_aromatic(r#"{:atoms ["N" "B"] :bonds [[0 1 "1"]] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#a"}] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]"}]}"#)]
    #[case::dative_not_aromatic(
        r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#a!"}]}"#
    )]
    #[case::aromatic_electrons(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#e2"}]}"#)]
    #[case::multicenter_electrons(r#"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,1,0]#e2"}]}"#)]
    #[case::noncovalent_intramolecular(r#"{:atoms ["N" "H"] :bonds [[0 1 "1"]] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd#I"}]}"#)]
    #[case::noncovalent_intermolecular(
        r#"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd#I!"}]}"#
    )]
    #[case::finite_set(r#"{:atoms ["C#v{1,2}" "C"] :bonds [[0 1 "1"]]}"#)]
    #[case::range(r#"{:atoms ["C#v(1..)" "C"] :bonds [[0 1 "1"]]}"#)]
    #[case::vacuous(r#"{:atoms ["C#v*"] :bonds []}"#)]
    fn test_incidence_constraint_validator_validate_determined(#[case] input: &str) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            IncidenceConstraintInvariantsValidator
                .validate(&molecule, ConnectedComponentsAlgorithm::Bfs,),
            Ok(Solution::Determined(()))
        );
    }

    #[rstest]
    #[case::atom(
        r#"{:atoms ["C#v2" "C"] :bonds [[0 1 "1"]]}"#,
        IncidenceConstraintInvariantsContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintForm::valence(2),
        }
    )]
    #[case::bond(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "1#a"]]}"#,
        IncidenceConstraintInvariantsContradiction::Bond {
            bond: BondId(0),
            constraint: BondConstraintForm::aromatic(true),
        }
    )]
    #[case::dative(
        r#"{:atoms ["N" "B"] :bonds [] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#a"}]}"#,
        IncidenceConstraintInvariantsContradiction::DativeBond {
            bond: DativeBondId(0),
            constraint: DativeBondConstraintForm::aromatic(true),
        }
    )]
    #[case::aromatic_system(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#e3"}]}"#,
        IncidenceConstraintInvariantsContradiction::AromaticSystem {
            system: AromaticSystemId(0),
            constraint: AromaticSystemConstraintForm::electron_count(3),
        }
    )]
    #[case::multicenter_bond(
        r#"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,1,0]#e3"}]}"#,
        IncidenceConstraintInvariantsContradiction::MulticenterBond {
            bond: MulticenterBondId(0),
            constraint: MulticenterBondConstraintForm::electron_count(3),
        }
    )]
    #[case::noncovalent_bond(
        r#"{:atoms ["N" "H"] :bonds [] :noncovalent-bonds [{:atoms [0 1] :attrs "Hbd#I"}]}"#,
        IncidenceConstraintInvariantsContradiction::NoncovalentBond {
            bond: NoncovalentBondId(0),
            constraint: NoncovalentBondConstraintForm::intramolecular(true),
        }
    )]
    fn test_incidence_constraint_validator_validate_contradictory(
        #[case] input: &str,
        #[case] expected: IncidenceConstraintInvariantsContradiction,
    ) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            IncidenceConstraintInvariantsValidator
                .validate(&molecule, ConnectedComponentsAlgorithm::Bfs,),
            Ok(Solution::Contradictory(expected))
        );
    }

    #[rstest]
    #[case::bond_order(r#"{:atoms ["C#v1" "C"] :bonds [[0 1 "*"]]}"#)]
    #[case::aromatic_valence(
        r#"{:atoms ["C#a1" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :attrs "*"}]}"#
    )]
    #[case::multicenter_valence(
        r#"{:atoms ["C#m1" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "*"}]}"#
    )]
    #[case::tetrahedral_coset(r#"{:atoms ["C#T1" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th*"}]}"#)]
    #[case::cis_trans_coset(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct*"}]}"#)]
    #[case::multi_donor_donated(r#"{:atoms ["N#d1" "N" "B"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :attrs "1"}]}"#)]
    #[case::multi_donor_accepted(r#"{:atoms ["N" "N" "B#t1"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :attrs "1"}]}"#)]
    #[case::multi_donor_aromatic(r#"{:atoms ["N" "N" "B"] :bonds [] :dative-bonds [{:donors [0 1] :acceptor 2 :attrs "1#a"}]}"#)]
    #[case::aromatic_electrons(
        r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]] :aromatic-systems [{:atoms [0 1] :attrs "*#e2"}]}"#
    )]
    #[case::multicenter_electrons(
        r#"{:atoms ["C" "C" "C"] :bonds [] :multicenter-bonds [{:atoms [0 1 2] :attrs "*#e2"}]}"#
    )]
    fn test_incidence_constraint_validator_validate_underdetermined(#[case] input: &str) {
        let molecule = mol_dsl!(input);

        assert_eq!(
            IncidenceConstraintInvariantsValidator
                .validate(&molecule, ConnectedComponentsAlgorithm::Bfs,),
            Ok(Solution::Underdetermined(()))
        );
    }
}
