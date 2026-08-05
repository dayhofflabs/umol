//! Aromaticity perception primitive.
//!
//! [`AromaticityPerception`] dispatches to one of three algorithms (Hückel
//! rule, HMO, Clar) selected by [`AromaticityModel`] and runs perception
//! against an AST. It is the shared core used by three top-level entities:
//! the resolver (validates `#a` hints filled in by atom-typing), the
//! aromatizer (discovers aromatic systems from a Kekulé bond-order layout),
//! and the validator (verifies pre-existing aromatic systems against the
//! model). [`AromaticityPerception::derive`] is the standard AST-facing
//! operation; [`AromaticityPerception::find_systems`] remains available when
//! the caller supplies another per-atom electron source.
//!
//! System insertion and bond marking are exposed via
//! [`AromaticityPerception::add_systems`].

pub mod clar;
pub mod hmo;
pub mod hueckel_rule;

use std::collections::BTreeSet;

pub use clar::{ClarAromaticity, ClarError};
pub use hmo::{HmoAromaticity, HmoError, HmoOutput};
pub use hueckel_rule::HueckelRuleAromaticity;
use thiserror::Error;
use umol_ast::ast::{
    AromaticSystemAst, AromaticSystemId, AromaticValenceAst, AtomId, AtomView, BondConstraintAst,
    BondId, BooleanAst, ElectronCountsAst, MoleculeAst, RingConfig, RingModel, RingSetKind,
    TransactionError, ValueAst,
};
use umol_graph_core::{ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm};
use umol_utils::solution::Solution;

use crate::ops::model::AromaticityModel;

/// Chemistry-level rejection: the algorithm decided the input doesn't satisfy
/// the model. Carried inside `Solution::Contradictory`.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityContradiction {
    #[error("hmo: invalid input: {0}")]
    HmoInvalidInput(String),
    #[error("clar: non-benzenoid input: {0}")]
    ClarNonBenzenoid(String),
    #[error("aromaticity inconsistency: {0}")]
    Inconsistency(#[from] AromaticityInconsistency),
}

/// Setup-level failure returned in `Err`, never inside `Solution`.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    #[error("hmo: missing parameters: {0}")]
    HmoMissingParameters(String),
    #[error(transparent)]
    Transaction(#[from] TransactionError),
}

/// Algorithms used by aromaticity perception independently of the chemistry
/// semantics in [`AromaticityModel`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AromaticityConfig {
    /// Algorithms used to construct the relevant ring set.
    pub ring_config: RingConfig,
    /// Algorithm used to separate HMO candidate components.
    pub connected_components_algorithm: ConnectedComponentsAlgorithm,
    /// Algorithm used to select disjoint Clar sextets.
    pub maximum_independent_set_algorithm: MaximumIndependentSetAlgorithm,
}

impl Default for AromaticityConfig {
    fn default() -> Self {
        Self {
            ring_config: RingConfig::default(),
            connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
            maximum_independent_set_algorithm: MaximumIndependentSetAlgorithm::BranchAndBound,
        }
    }
}

/// Policy-free result of aromaticity perception and constraint/entity comparison.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticityDerivation {
    /// Aromatic systems accepted by the selected model.
    pub systems: Vec<(Vec<AtomId>, AromaticSystemAst)>,
    /// Constraint failures, entity failures, and independently valid mismatches.
    pub inconsistencies: Vec<AromaticityInconsistency>,
}

/// Policy-free classification of aromatic constraint and entity inconsistencies.
#[derive(Clone, Copy, Debug, Error, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticityInconsistency {
    #[error("aromatic valence at atom {atom:?} cannot produce a valid aromatic system")]
    AromaticValenceFailure { atom: AtomId },
    #[error("aromatic system {system:?} is not realizable under the selected model")]
    AromaticSystemFailure { system: AromaticSystemId },
    #[error("aromatic valence at atom {atom:?} disagrees with aromatic system {system:?}")]
    AromaticValenceMismatch {
        atom: AtomId,
        system: AromaticSystemId,
    },
    #[error("aromatic constraint at bond {bond:?} disagrees with aromatic system {system:?}")]
    AromaticBondConstraintMismatch {
        bond: BondId,
        system: AromaticSystemId,
    },
}

#[derive(Clone, Debug)]
pub enum AromaticityPerception {
    HueckelRule(HueckelRuleAromaticity),
    Hmo(HmoAromaticity),
    Clar(ClarAromaticity),
}

impl AromaticityPerception {
    pub fn new(model: &AromaticityModel) -> Self {
        match model {
            AromaticityModel::HueckelRule { scope, ring_limits } => Self::HueckelRule(
                HueckelRuleAromaticity::new(scope.clone(), ring_limits.clone()),
            ),
            AromaticityModel::Hmo {
                scope,
                stabilization_threshold,
            } => Self::Hmo(HmoAromaticity::new(scope.clone(), *stabilization_threshold)),
            AromaticityModel::Clar { .. } => Self::Clar(ClarAromaticity),
        }
    }

    /// Find candidate aromatic systems via the configured algorithm.
    /// The closure `electrons_at` returns each atom's π contribution if the
    /// atom is aromatic-eligible, else `None`.
    #[allow(clippy::complexity)]
    pub fn find_systems<F>(
        &self,
        ast: &MoleculeAst,
        config: AromaticityConfig,
        electrons_at: F,
    ) -> Result<
        Solution<Vec<(Vec<AtomId>, AromaticSystemAst)>, AromaticityContradiction>,
        AromaticityError,
    >
    where
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let model = self.ring_request();
        let rings = ast.rings(model, config.ring_config).into_ring_set();

        let systems = match self {
            Self::HueckelRule(m) => m.find_from_rings(ast, &rings, &electrons_at),
            Self::Hmo(m) => match m.find_from_rings(
                ast,
                &rings,
                config.connected_components_algorithm,
                &electrons_at,
            ) {
                Ok(systems) => systems,
                Err(HmoError::MissingParameters(s)) => {
                    return Err(AromaticityError::HmoMissingParameters(s));
                }
                Err(HmoError::InvalidInput(s)) => {
                    return Ok(Solution::Contradictory(
                        AromaticityContradiction::HmoInvalidInput(s),
                    ));
                }
                Err(HmoError::UndeterminedAtom(_)) => {
                    return Ok(Solution::Underdetermined(Vec::new()));
                }
            },
            Self::Clar(m) => match m.find_from_rings(
                ast,
                &rings,
                config.maximum_independent_set_algorithm,
                &electrons_at,
            ) {
                Ok(systems) => systems,
                Err(ClarError::NonBenzenoid(s)) => {
                    return Ok(Solution::Contradictory(
                        AromaticityContradiction::ClarNonBenzenoid(s),
                    ));
                }
            },
        };

        let mut sorted = systems;
        sorted.sort_by(|a, b| a.0.first().cmp(&b.0.first()));
        Ok(Solution::Determined(sorted))
    }

    /// Perceive aromatic systems from aromatic-valence constraints, independently assess stored
    /// systems from their electron contributions, and classify their relationship.
    pub fn derive(
        &self,
        ast: &MoleculeAst,
        config: AromaticityConfig,
    ) -> Result<Solution<AromaticityDerivation, AromaticityContradiction>, AromaticityError> {
        if ast.atoms().iter().any(|atom| {
            matches!(
                atom.ast.constraints.aromatic_valence(),
                Some(AromaticValenceAst::Aromatic(value))
                    if !matches!(value, ValueAst::Lit(_))
            )
        }) || ast
            .aromatic_systems()
            .iter()
            .any(|system| matches!(system.ast.electrons, ElectronCountsAst::Undetermined))
        {
            return Ok(Solution::Underdetermined(AromaticityDerivation::default()));
        }

        let systems = match self.find_systems(ast, config, |atom| {
            match atom.ast.constraints.aromatic_valence() {
                Some(AromaticValenceAst::Aromatic(ValueAst::Lit(value))) => {
                    u8::try_from(*value).ok()
                }
                Some(AromaticValenceAst::Aromatic(_)) | Some(AromaticValenceAst::NotAromatic) => {
                    None
                }
                Some(AromaticValenceAst::Undetermined) | None => match atom.aromatic_valence() {
                    ValueAst::Lit(value) => u8::try_from(value).ok(),
                    _ => None,
                },
            }
        })? {
            Solution::Determined(systems) => systems,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(AromaticityDerivation::default()));
            }
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };

        let system_members: Vec<BTreeSet<AtomId>> = systems
            .iter()
            .map(|(atoms, _)| atoms.iter().copied().collect())
            .collect();
        let accepted_atoms: BTreeSet<AtomId> = system_members.iter().flatten().copied().collect();
        let mut inconsistencies = BTreeSet::new();

        for atom in ast.atoms().iter() {
            if matches!(
                atom.ast.constraints.aromatic_valence(),
                Some(AromaticValenceAst::Aromatic(_))
            ) && !accepted_atoms.contains(&atom.id)
            {
                inconsistencies
                    .insert(AromaticityInconsistency::AromaticValenceFailure { atom: atom.id });
            }
        }

        let mut valid_existing = Vec::new();
        for existing in ast.aromatic_systems().iter() {
            let ElectronCountsAst::Lit(existing_electrons) = &existing.ast.electrons else {
                return Ok(Solution::Underdetermined(AromaticityDerivation::default()));
            };
            let existing_atoms: Vec<AtomId> = existing.atom_ids().collect();
            if existing_electrons.len() != existing_atoms.len() {
                inconsistencies.insert(AromaticityInconsistency::AromaticSystemFailure {
                    system: existing.id,
                });
                continue;
            }

            let existing_contributions: Vec<(AtomId, i64)> = existing_atoms
                .iter()
                .copied()
                .zip(existing_electrons.iter().copied())
                .collect();

            let perceived = match self.find_systems(ast, config, |atom| {
                existing_contributions
                    .iter()
                    .find_map(|&(candidate, electrons)| {
                        (candidate == atom.id).then(|| u8::try_from(electrons).ok())
                    })
                    .flatten()
            })? {
                Solution::Determined(perceived) => perceived,
                Solution::Underdetermined(_) => {
                    return Ok(Solution::Underdetermined(AromaticityDerivation::default()));
                }
                Solution::Contradictory(_) => {
                    inconsistencies.insert(AromaticityInconsistency::AromaticSystemFailure {
                        system: existing.id,
                    });
                    continue;
                }
            };

            let existing_members: BTreeSet<AtomId> = existing_atoms.iter().copied().collect();
            let valid = perceived.iter().any(|(perceived_atoms, perceived_system)| {
                let perceived_members: BTreeSet<AtomId> = perceived_atoms.iter().copied().collect();
                if perceived_members != existing_members {
                    return false;
                }
                let ElectronCountsAst::Lit(perceived_electrons) = &perceived_system.electrons
                else {
                    return false;
                };
                perceived_atoms
                    .iter()
                    .copied()
                    .zip(perceived_electrons.iter().copied())
                    .all(|(atom, electrons)| {
                        existing_contributions
                            .iter()
                            .find_map(|&(candidate, existing)| {
                                (candidate == atom).then_some(existing)
                            })
                            == Some(electrons)
                    })
            });
            if valid {
                valid_existing.push((existing.id, existing_contributions));
            } else {
                inconsistencies.insert(AromaticityInconsistency::AromaticSystemFailure {
                    system: existing.id,
                });
            }
        }

        for (system, contributions) in valid_existing {
            let members: BTreeSet<AtomId> = contributions.iter().map(|&(atom, _)| atom).collect();
            let has_matching_candidate =
                system_members.iter().any(|candidate| candidate == &members);
            for atom in ast.atoms().iter() {
                let Some(constraint) = atom.ast.constraints.aromatic_valence() else {
                    continue;
                };
                let mismatch = match constraint {
                    AromaticValenceAst::Aromatic(ValueAst::Lit(expected)) => {
                        has_matching_candidate
                            && contributions
                                .iter()
                                .find_map(|&(candidate, actual)| {
                                    (candidate == atom.id).then_some(actual)
                                })
                                .is_some_and(|actual| actual != *expected)
                    }
                    AromaticValenceAst::Aromatic(_)
                    | AromaticValenceAst::NotAromatic
                    | AromaticValenceAst::Undetermined => false,
                };
                if mismatch {
                    inconsistencies.insert(AromaticityInconsistency::AromaticValenceMismatch {
                        atom: atom.id,
                        system,
                    });
                }
            }

            for bond in ast.aromatic_system(system).bonds() {
                if matches!(bond.ast.constraints.aromatic(), BooleanAst::Lit(false)) {
                    inconsistencies.insert(
                        AromaticityInconsistency::AromaticBondConstraintMismatch {
                            bond: bond.id,
                            system,
                        },
                    );
                }
            }
        }

        Ok(Solution::Determined(AromaticityDerivation {
            systems,
            inconsistencies: inconsistencies.into_iter().collect(),
        }))
    }

    /// Add perceived systems to the AST.
    pub fn add_systems(
        &self,
        ast: &mut MoleculeAst,
        systems: Vec<(Vec<AtomId>, AromaticSystemAst)>,
    ) {
        if systems.is_empty() {
            return;
        }
        let mut builder = ast.edit();
        let new_indices: Vec<AromaticSystemId> = systems
            .into_iter()
            .map(|(atoms, system_ast)| builder.add_aromatic_system(atoms, system_ast))
            .collect();
        *ast = builder.build();

        let bond_ids: Vec<BondId> = new_indices
            .iter()
            .flat_map(|&idx| ast.aromatic_system(idx).bond_ids().collect::<Vec<_>>())
            .collect();
        for bond_id in bond_ids {
            let bond = ast.bond_mut(bond_id);
            bond.ast
                .constraints
                .set(BondConstraintAst::Aromatic(BooleanAst::Lit(true)));
        }
    }

    fn ring_request(&self) -> RingModel {
        let max_ring_size = match self {
            Self::HueckelRule(m) => m.ring_limits.max_ring_size,
            Self::Hmo(_) => 22,
            Self::Clar(_) => 6,
        };
        RingModel {
            kind: RingSetKind::Relevant,
            max_ring_size,
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AromaticSystemId, AromaticValenceAst, AtomAst, AtomConstraintAst, AtomConstraintKey,
        AtomId, BondAst, BondConstraintKey, ElectronCountsAst, MoleculeAst, MoleculeEntries,
        UnpairedElectronsAst, ValueAst,
    };
    use umol_ast::{mol_dsl, mol_dsl_ground};
    use umol_chem::element::Element;
    use umol_graph_core::{RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm};

    use super::*;
    use crate::ops::model::{ElementScope, RingLimits};

    #[rstest]
    fn test_aromaticity_config_default() {
        assert_eq!(
            AromaticityConfig::default(),
            AromaticityConfig {
                ring_config: RingConfig {
                    simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                },
                connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                maximum_independent_set_algorithm: MaximumIndependentSetAlgorithm::BranchAndBound,
            }
        );
    }

    fn any_hueckel() -> AromaticityPerception {
        AromaticityPerception::new(&AromaticityModel::HueckelRule {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        })
    }

    fn aromatic_valence_lit(ast: &MoleculeAst, idx: AtomId) -> Option<i64> {
        match ast
            .atom(idx)
            .ast
            .constraints
            .get(AtomConstraintKey::AromaticValence)?
        {
            AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(n))) => {
                Some(*n)
            }
            _ => None,
        }
    }

    fn aromatic(element: Element, pi: i64) -> AtomAst {
        let mut atom = AtomAst::from_element(element);
        atom.charge = ValueAst::Lit(0);
        atom.unpaired_electrons = UnpairedElectronsAst::closed_shell();
        atom.constraints.set(AtomConstraintAst::AromaticValence(
            AromaticValenceAst::Aromatic(ValueAst::Lit(pi)),
        ));
        atom
    }

    fn benzene() -> MoleculeAst {
        let atoms: Vec<AtomAst> = (0..6).map(|_| aromatic(Element::C, 1)).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn pyrrole() -> MoleculeAst {
        let atoms = vec![
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ];
        let bonds: Vec<_> = (0..5)
            .map(|i| (AtomId(i), AtomId((i + 1) % 5), BondAst::from_order(1)))
            .collect();
        MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn run_full(
        perception: &AromaticityPerception,
        ast: &mut MoleculeAst,
    ) -> Solution<(), AromaticityContradiction> {
        let outcome = perception
            .find_systems(ast, AromaticityConfig::default(), |v| {
                match v
                    .ast
                    .constraints
                    .aromatic_valence()
                    .unwrap_or(&AromaticValenceAst::Undetermined)
                {
                    AromaticValenceAst::Aromatic(ValueAst::Lit(n)) if *n >= 0 => Some(*n as u8),
                    _ => None,
                }
            })
            .unwrap();
        match outcome {
            Solution::Determined(systems) => {
                perception.add_systems(ast, systems);
                Solution::Determined(())
            }
            Solution::Underdetermined(_) => Solution::Underdetermined(()),
            Solution::Contradictory(c) => Solution::Contradictory(c),
        }
    }

    #[rstest]
    #[case::daylight_furan(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 0 "1#a"]]
        }"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![(
                vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
                AromaticSystemAst::from_electrons(vec![2, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsAst::closed_shell()),
            )],
            inconsistencies: vec![],
        }),
    )]
    #[case::mdl_furan(
        AromaticityModel::mdl(),
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 0 "1#a"]]
        }"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![],
            inconsistencies: vec![
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(1) },
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(2) },
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(3) },
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(4) },
            ],
        }),
    )]
    #[case::missing_and_extra_projections(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#h" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#a" "C"]
            :bonds [[0 1 "1#a!"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"] [6 7 "1#a"]]
            :aromatic-systems [{
                :atoms [0 1 2 3 4 5]
                :type "[1,1,1,1,1,1]"
            }]
        }"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![(
                vec![
                    AtomId(0),
                    AtomId(1),
                    AtomId(2),
                    AtomId(3),
                    AtomId(4),
                    AtomId(5),
                ],
                AromaticSystemAst::from_electrons(vec![1, 1, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsAst::closed_shell()),
            )],
            inconsistencies: vec![
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(6) },
                AromaticityInconsistency::AromaticBondConstraintMismatch {
                    bond: BondId(0),
                    system: AromaticSystemId(0),
                },
            ],
        }),
    )]
    #[case::vacuous_projections(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#a*" "C#a*" "C#a*" "C#a*" "C#a*" "C#a*"]
            :bonds [[0 1 "1#a*"] [1 2 "1#a*"] [2 3 "1#a*"] [3 4 "1#a*"]
                    [4 5 "1#a*"] [5 0 "1#a*"]]
            :aromatic-systems [{
                :atoms [0 1 2 3 4 5]
                :type "[1,1,1,1,1,1]"
            }]
        }"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![(
                vec![
                    AtomId(0),
                    AtomId(1),
                    AtomId(2),
                    AtomId(3),
                    AtomId(4),
                    AtomId(5),
                ],
                AromaticSystemAst::from_electrons(vec![1, 1, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsAst::closed_shell()),
            )],
            inconsistencies: vec![],
        }),
    )]
    #[case::electron_contribution(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
            :aromatic-systems [{
                :atoms [0 1 2 3 4 5]
                :type "[2,0,1,1,1,1]"
            }]
        }"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![(
                vec![
                    AtomId(0),
                    AtomId(1),
                    AtomId(2),
                    AtomId(3),
                    AtomId(4),
                    AtomId(5),
                ],
                AromaticSystemAst::from_electrons(vec![1, 1, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsAst::closed_shell()),
            )],
            inconsistencies: vec![
                AromaticityInconsistency::AromaticValenceMismatch {
                    atom: AtomId(0),
                    system: AromaticSystemId(0),
                },
                AromaticityInconsistency::AromaticValenceMismatch {
                    atom: AtomId(1),
                    system: AromaticSystemId(0),
                },
            ],
        }),
    )]
    #[case::conformant_system(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
            :aromatic-systems [{
                :atoms [0 1 2 3 4 5]
                :type "[1,1,1,1,1,1]"
            }]
        }"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![(
                vec![
                    AtomId(0),
                    AtomId(1),
                    AtomId(2),
                    AtomId(3),
                    AtomId(4),
                    AtomId(5),
                ],
                AromaticSystemAst::from_electrons(vec![1, 1, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsAst::closed_shell()),
            )],
            inconsistencies: vec![],
        }),
    )]
    #[case::existing_system_rejected(
        AromaticityModel::mdl(),
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 0 "1#a"]]
            :aromatic-systems [{
                :atoms [0 1 2 3 4]
                :type "[2,1,1,1,1]"
            }]
        }"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![],
            inconsistencies: vec![
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(1) },
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(2) },
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(3) },
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(4) },
                AromaticityInconsistency::AromaticSystemFailure {
                    system: AromaticSystemId(0),
                },
            ],
        }),
    )]
    #[case::non_ground_assertion(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#a+" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
        }"#),
        Solution::Underdetermined(AromaticityDerivation::default()),
    )]
    fn test_aromaticity_perception_derive(
        #[case] model: AromaticityModel,
        #[case] ast: MoleculeAst,
        #[case] expected: Solution<AromaticityDerivation, AromaticityContradiction>,
    ) {
        assert_eq!(
            AromaticityPerception::new(&model)
                .derive(&ast, AromaticityConfig::default())
                .unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_aromaticity_perception_hueckel_rule_benzene_writes_system() {
        let perception = AromaticityPerception::new(&AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        });
        let mut ast = benzene();
        let solution = run_full(&perception, &mut ast);
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(ast.aromatic_systems().count(), 1);
        let system = ast.aromatic_system(AromaticSystemId(0));
        let atoms: Vec<AtomId> = system.atom_ids().collect();
        assert_eq!(atoms.len(), 6);
        let aromatic_bond_count = ast
            .bonds()
            .iter()
            .filter(|view| view.ast.constraints.contains(BondConstraintKey::Aromatic))
            .count();
        assert_eq!(aromatic_bond_count, 6);
    }

    #[rstest]
    fn test_aromaticity_perception_clar_rejects_heterocycle() {
        let perception = AromaticityPerception::new(&AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        });
        let mut ast = pyrrole();
        let solution = run_full(&perception, &mut ast);
        assert!(matches!(
            solution,
            Solution::Contradictory(AromaticityContradiction::ClarNonBenzenoid(_))
        ));
    }

    #[rstest]
    #[case::cyclopropenium_cation(
        mol_dsl_ground!(r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#),
        0, vec![1, 1, 0], vec![0, 0, 1], vec![1, 1, 0],
    )]
    #[case::cot_dianion(
        mol_dsl_ground!(r#"{:atoms ["C #h #a" "C #c- #h #a2" "C #h #a" "C #h #a"
                                "C #h #a" "C #c- #h #a2" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"]
                               [4 5 "1"] [5 6 "1"] [6 7 "1"] [7 0 "1"]]}"#),
        0, vec![1, 2, 1, 1, 1, 2, 1, 1],
        vec![0, -1, 0, 0, 0, -1, 0, 0], vec![1, 2, 1, 1, 1, 2, 1, 1],
    )]
    #[case::s4_dication(
        mol_dsl_ground!(r#"{:atoms ["S #c+ #n1 #a" "S #n1 #a2" "S #c+ #n1 #a" "S #n1 #a2"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]]}"#),
        0, vec![1, 2, 1, 2], vec![1, 0, 1, 0], vec![1, 2, 1, 2],
    )]
    #[case::boratabenzene_anion(
        mol_dsl_ground!(r#"{:atoms ["B #c- #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        0, vec![1; 6], vec![-1, 0, 0, 0, 0, 0], vec![1; 6],
    )]
    #[case::borepin(
        mol_dsl_ground!(r#"{:atoms ["B #h #a0" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 0 "1"]]}"#),
        0, vec![0, 1, 1, 1, 1, 1, 1], vec![0; 7], vec![0, 1, 1, 1, 1, 1, 1],
    )]
    #[case::pyridinium(
        mol_dsl_ground!(r#"{:atoms ["N #c+ #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        0, vec![1; 6], vec![1, 0, 0, 0, 0, 0], vec![1; 6],
    )]
    #[case::pyrylium(
        mol_dsl_ground!(r#"{:atoms ["O #c+ #n1 #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        0, vec![1; 6], vec![1, 0, 0, 0, 0, 0], vec![1; 6],
    )]
    #[case::pyrrole(
        mol_dsl_ground!(r#"{:atoms ["N #h #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        0, vec![2, 1, 1, 1, 1], vec![0; 5], vec![2, 1, 1, 1, 1],
    )]
    #[case::furan(
        mol_dsl_ground!(r#"{:atoms ["O #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        0, vec![2, 1, 1, 1, 1], vec![0; 5], vec![2, 1, 1, 1, 1],
    )]
    #[case::thiophene(
        mol_dsl_ground!(r#"{:atoms ["S #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        0, vec![2, 1, 1, 1, 1], vec![0; 5], vec![2, 1, 1, 1, 1],
    )]
    fn test_aromaticity_perception_add_systems(
        #[case] mut ast: MoleculeAst,
        #[case] system_charge: i64,
        #[case] electrons: Vec<i64>,
        #[case] atom_charges: Vec<i64>,
        #[case] aromatic_valences: Vec<i64>,
    ) {
        let outcome = any_hueckel()
            .find_systems(&ast, AromaticityConfig::default(), |v| {
                match v
                    .ast
                    .constraints
                    .aromatic_valence()
                    .unwrap_or(&AromaticValenceAst::Undetermined)
                {
                    AromaticValenceAst::Aromatic(ValueAst::Lit(n)) if *n >= 0 => Some(*n as u8),
                    _ => None,
                }
            })
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut ast, systems);

        let system = ast.aromatic_system(AromaticSystemId(0));
        assert_eq!(system.ast.charge, ValueAst::Lit(system_charge));
        assert_eq!(system.ast.electrons, ElectronCountsAst::Lit(electrons));
        for (i, (q, k)) in atom_charges
            .iter()
            .zip(aromatic_valences.iter())
            .enumerate()
        {
            let idx = AtomId(i as u32);
            assert_eq!(ast.atom(idx).ast.charge, ValueAst::Lit(*q));
            assert_eq!(aromatic_valence_lit(&ast, idx), Some(*k));
        }
    }

    #[rstest]
    fn test_aromaticity_perception_hueckel_rule_no_aromatic_atom_returns_determined() {
        let perception = AromaticityPerception::new(&AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        });
        let atoms: Vec<AtomAst> = (0..6).map(|_| AtomAst::from_element(Element::C)).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        let mut ast = MoleculeAst::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        });
        let solution = run_full(&perception, &mut ast);
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(ast.aromatic_systems().count(), 0);
        let any_aromatic = ast
            .bonds()
            .iter()
            .any(|view| view.ast.constraints.contains(BondConstraintKey::Aromatic));
        assert!(!any_aromatic);
    }
}
