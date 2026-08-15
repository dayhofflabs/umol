//! Model-independent constraints derived from entity fields and directly incident entities.

use thiserror::Error;
use umol_graph_core::ConnectedComponentsAlgorithm;
use umol_graph_ir::ir::{
    AromaticSystemConstraintForm, AromaticSystemId, AtomConstraintForm, AtomConstraintKey, AtomId,
    BondConstraintForm, BondConstraintKey, BondId, DativeBondConstraintForm, DativeBondId, Entity,
    Lattice, Molecule, MulticenterBondConstraintForm, MulticenterBondId,
    NoncovalentBondConstraintForm, NoncovalentBondId,
};
use umol_utils::solution::Solution;

use super::{ConstraintInvariantsError, DerivedKind};

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
        reading: DerivedKind,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let mut bond_components = None;
        let mut any_underdetermined = false;

        for id in molecule.atoms().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_atom(molecule, id, reading)?,
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for id in molecule.bonds().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_bond(molecule, id, reading)?,
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for id in molecule.dative_bonds().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_dative_bond(molecule, id, reading)?,
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for id in molecule.aromatic_systems().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_aromatic_system(molecule, id, reading)?,
                &mut any_underdetermined,
            ) {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        for id in molecule.multicenter_bonds().ids() {
            if let Some(contradiction) = observe(
                self.validate_molecule_multicenter_bond(molecule, id, reading)?,
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
        reading: DerivedKind,
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
            validate_atom_constraint(molecule, atom_id, constraint, reading)
        })))
    }

    /// Validate one inline atom constraint selected by its container key.
    pub fn validate_molecule_atom_constraint(
        &self,
        molecule: &Molecule,
        atom_id: AtomId,
        key: AtomConstraintKey,
        reading: DerivedKind,
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
            .asserted(key)
            .map_or(Solution::Determined(()), |constraint| {
                validate_atom_constraint(molecule, atom_id, constraint, reading)
            }))
    }

    /// Validate all inline incidence constraints on one molecule bond.
    pub fn validate_molecule_bond(
        &self,
        molecule: &Molecule,
        bond_id: BondId,
        reading: DerivedKind,
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
            |constraint| validate_bond_constraint(molecule, bond_id, constraint, reading),
        )))
    }

    /// Validate one inline localized-bond constraint selected by its container key.
    pub fn validate_molecule_bond_constraint(
        &self,
        molecule: &Molecule,
        bond_id: BondId,
        key: BondConstraintKey,
        reading: DerivedKind,
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
            .asserted(key)
            .and_then(|constraint| validate_bond_constraint(molecule, bond_id, constraint, reading))
            .unwrap_or(Solution::Determined(())))
    }

    /// Validate all inline incidence constraints on one molecule dative bond.
    pub fn validate_molecule_dative_bond(
        &self,
        molecule: &Molecule,
        bond_id: DativeBondId,
        reading: DerivedKind,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let bond = molecule.dative_bonds().get(bond_id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::DativeBond(bond_id),
            },
        )?;
        Ok(conjunction(bond.constraints().iter().filter_map(
            |constraint| validate_dative_bond_constraint(molecule, bond_id, constraint, reading),
        )))
    }

    /// Validate all inline incidence constraints on one molecule aromatic system.
    pub fn validate_molecule_aromatic_system(
        &self,
        molecule: &Molecule,
        system_id: AromaticSystemId,
        reading: DerivedKind,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let system = molecule.aromatic_systems().get(system_id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::AromaticSystem(system_id),
            },
        )?;
        Ok(conjunction(system.constraints().iter().map(|constraint| {
            validate_aromatic_system_constraint(molecule, system_id, constraint, reading)
        })))
    }

    /// Validate all inline incidence constraints on one molecule multicenter bond.
    pub fn validate_molecule_multicenter_bond(
        &self,
        molecule: &Molecule,
        bond_id: MulticenterBondId,
        reading: DerivedKind,
    ) -> Result<Solution<(), IncidenceConstraintInvariantsContradiction>, ConstraintInvariantsError>
    {
        let bond = molecule.multicenter_bonds().get(bond_id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::MulticenterBond(bond_id),
            },
        )?;
        Ok(conjunction(bond.constraints().iter().map(|constraint| {
            validate_multicenter_bond_constraint(molecule, bond_id, constraint, reading)
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

/// Validate one non-ring atom constraint against the selected derived
/// reading: violation ⇔ the assertion and the reading fail to meet. Ring
/// constraints are determined identities here — the ring validator owns them.
pub fn validate_atom_constraint(
    molecule: &Molecule,
    atom_id: AtomId,
    constraint: &AtomConstraintForm,
    reading: DerivedKind,
) -> Solution<(), IncidenceConstraintInvariantsContradiction> {
    if matches!(
        constraint,
        AtomConstraintForm::RingDegree(_)
            | AtomConstraintForm::RingValence(_)
            | AtomConstraintForm::RingMembership(_)
    ) {
        return Solution::Determined(());
    }
    if constraint.is_undetermined() {
        return Solution::Determined(());
    }
    let atom = molecule.atom(atom_id);
    // Multi-donor dative incidence has no defined per-atom projection pending
    // the coordination/haptic entity split in discussion doc 117.
    let unsupported = match constraint {
        AtomConstraintForm::DonatedPairs(_) => atom
            .dative_bonds()
            .any(|bond| bond.donor_count() != 1 && bond.donor_ids().any(|donor| donor == atom.id)),
        AtomConstraintForm::AcceptedPairs(_) => atom
            .dative_bonds()
            .any(|bond| bond.donor_count() != 1 && bond.acceptor_id() == atom.id),
        _ => false,
    };
    if unsupported {
        return Solution::Underdetermined(());
    }
    match derived_atom(molecule, atom_id, constraint, reading) {
        Some(derived) => evaluate(
            constraint,
            &derived,
            atom_contradiction(atom_id, constraint),
        ),
        None => Solution::Underdetermined(()),
    }
}

fn derived_atom(
    molecule: &Molecule,
    atom_id: AtomId,
    constraint: &AtomConstraintForm,
    reading: DerivedKind,
) -> Option<AtomConstraintForm> {
    let constraints = molecule.atom(atom_id).constraints();
    match reading {
        DerivedKind::Derived => constraints.derived(constraint.key()),
        DerivedKind::DerivedComplete => constraints.derived_complete(constraint.key()),
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
    reading: DerivedKind,
) -> Option<Solution<(), IncidenceConstraintInvariantsContradiction>> {
    if matches!(constraint, BondConstraintForm::RingMembership(_)) {
        return None;
    }
    if constraint.is_undetermined() {
        return Some(Solution::Determined(()));
    }
    let constraints = molecule.bond(bond_id).constraints();
    let derived = match reading {
        DerivedKind::Derived => constraints.derived(constraint.key()),
        DerivedKind::DerivedComplete => constraints.derived_complete(constraint.key()),
    };
    Some(match derived {
        Some(derived) => evaluate(
            constraint,
            &derived,
            bond_contradiction(bond_id, constraint),
        ),
        None => Solution::Underdetermined(()),
    })
}

pub fn validate_dative_bond_constraint(
    molecule: &Molecule,
    bond_id: DativeBondId,
    constraint: &DativeBondConstraintForm,
    reading: DerivedKind,
) -> Option<Solution<(), IncidenceConstraintInvariantsContradiction>> {
    if matches!(constraint, DativeBondConstraintForm::RingMembership(_)) {
        return None;
    }
    if constraint.is_undetermined() {
        return Some(Solution::Determined(()));
    }
    // The views' derivation yields no value for multi-donor incidence,
    // pending the coordination/haptic entity split in discussion doc 117;
    // the `None` fallthrough reads it as undecided.
    let constraints = molecule.dative_bond(bond_id).constraints();
    let derived = match reading {
        DerivedKind::Derived => constraints.derived(constraint.key()),
        DerivedKind::DerivedComplete => constraints.derived_complete(constraint.key()),
    };
    Some(match derived {
        Some(derived) => evaluate(
            constraint,
            &derived,
            IncidenceConstraintInvariantsContradiction::DativeBond {
                bond: bond_id,
                constraint: constraint.clone(),
            },
        ),
        None => Solution::Underdetermined(()),
    })
}

pub fn validate_aromatic_system_constraint(
    molecule: &Molecule,
    system_id: AromaticSystemId,
    constraint: &AromaticSystemConstraintForm,
    reading: DerivedKind,
) -> Solution<(), IncidenceConstraintInvariantsContradiction> {
    if constraint.is_undetermined() {
        return Solution::Determined(());
    }
    let constraints = molecule.aromatic_system(system_id).constraints();
    let derived = match reading {
        DerivedKind::Derived => constraints.derived(constraint.key()),
        DerivedKind::DerivedComplete => constraints.derived_complete(constraint.key()),
    };
    match derived {
        Some(derived) => evaluate(
            constraint,
            &derived,
            IncidenceConstraintInvariantsContradiction::AromaticSystem {
                system: system_id,
                constraint: constraint.clone(),
            },
        ),
        None => Solution::Underdetermined(()),
    }
}

pub fn validate_multicenter_bond_constraint(
    molecule: &Molecule,
    bond_id: MulticenterBondId,
    constraint: &MulticenterBondConstraintForm,
    reading: DerivedKind,
) -> Solution<(), IncidenceConstraintInvariantsContradiction> {
    if constraint.is_undetermined() {
        return Solution::Determined(());
    }
    let constraints = molecule.multicenter_bond(bond_id).constraints();
    let derived = match reading {
        DerivedKind::Derived => constraints.derived(constraint.key()),
        DerivedKind::DerivedComplete => constraints.derived_complete(constraint.key()),
    };
    match derived {
        Some(derived) => evaluate(
            constraint,
            &derived,
            IncidenceConstraintInvariantsContradiction::MulticenterBond {
                bond: bond_id,
                constraint: constraint.clone(),
            },
        ),
        None => Solution::Underdetermined(()),
    }
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
    use umol_graph_ir::ir::{AromaticValenceForm, NumForm};
    use umol_graph_ir::mol_dsl;

    use super::*;

    #[rstest]
    #[case::complete_contradicts(
        DerivedKind::DerivedComplete,
        Solution::Contradictory(IncidenceConstraintInvariantsContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(
                NumForm::Undetermined,
            )),
        })
    )]
    #[case::derived_undecided(DerivedKind::Derived, Solution::Underdetermined(()))]
    fn test_validate_atom_constraint_reading(
        #[case] reading: DerivedKind,
        #[case] expected: Solution<(), IncidenceConstraintInvariantsContradiction>,
    ) {
        // An aromatic assertion on an atom in no stored system: the closure
        // reads the absence as `NotAromatic` and contradicts; the open
        // reading leaves the assertion undecided.
        let molecule = mol_dsl!(r#"{:atoms ["C#a+"] :bonds []}"#);
        assert_eq!(
            validate_atom_constraint(
                &molecule,
                AtomId(0),
                &AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(
                    NumForm::Undetermined,
                )),
                reading,
            ),
            expected
        );
    }

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
            IncidenceConstraintInvariantsValidator.validate_molecule_atom_constraint(
                &molecule,
                atom,
                key,
                DerivedKind::DerivedComplete
            ),
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
            IncidenceConstraintInvariantsValidator.validate_molecule_bond_constraint(
                &molecule,
                bond,
                key,
                DerivedKind::DerivedComplete
            ),
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
            IncidenceConstraintInvariantsValidator.validate(
                &molecule,
                ConnectedComponentsAlgorithm::Bfs,
                DerivedKind::DerivedComplete,
            ),
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
            IncidenceConstraintInvariantsValidator.validate(
                &molecule,
                ConnectedComponentsAlgorithm::Bfs,
                DerivedKind::DerivedComplete,
            ),
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
            IncidenceConstraintInvariantsValidator.validate(
                &molecule,
                ConnectedComponentsAlgorithm::Bfs,
                DerivedKind::DerivedComplete,
            ),
            Ok(Solution::Underdetermined(()))
        );
    }
}
