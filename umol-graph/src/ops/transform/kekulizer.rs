//! Kekulizer: aromatic-system form → Kekulé bond orders.
//!
//! For each aromatic system in the input, finds a perfect matching on the
//! induced subgraph (atoms in the system, bonds between them). Matched bonds
//! become order 2; unmatched bonds become order 1. The aromatic constraint is
//! removed from the system's bonds and from the system's atoms; the system
//! entry itself is removed at the end.
//!
//! The matching algorithm is configurable via [`KekulizationModel`]; the
//! atom processing order — which controls determinism — is fixed at
//! construction time and is the caller's responsibility (e.g., from a
//! nauty/Traces canonical labeling).
//!
//! TODO: Expand to charged systems.

use std::collections::HashSet;

use thiserror::Error;
use umol_ast::ast::{
    AromaticSystemId, AromaticSystemView, AtomConstraintKey, AtomId, BondConstraintKey, BondId,
    ElectronCountsAst, MoleculeAst, ValueAst,
};
use umol_graph_core::{NodeId, PerfectMatchingAlgorithm};

use crate::ops::transform::Transformer;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum KekulizerError {
    #[error("electron contributions are undetermined for aromatic system {0:?}")]
    UndeterminedElectrons(AromaticSystemId),
    #[error("charge is undetermined for aromatic system {0:?}")]
    UndeterminedCharge(AromaticSystemId),
    #[error("spin is undetermined for aromatic system {0:?}")]
    UndeterminedSpin(AromaticSystemId),
    #[error(
        "aromatic system {system:?} has {member_count} members but {electron_count} electron contributions"
    )]
    ElectronCountMismatch {
        system: AromaticSystemId,
        member_count: usize,
        electron_count: usize,
    },
    #[error(
        "atom {atom:?} in aromatic system {system:?} has unsupported electron contribution {contribution}"
    )]
    UnsupportedElectronContribution {
        system: AromaticSystemId,
        atom: AtomId,
        contribution: i64,
    },
    #[error("aromatic system {0:?} is open-shell")]
    OpenShell(AromaticSystemId),
    #[error("aromatic system {system:?} has {count} prescribed exposed atoms")]
    MultiplePrescribedExposures {
        system: AromaticSystemId,
        count: usize,
    },
    #[error("aromatic system {system:?} has unsupported charge {charge}")]
    UnsupportedCharge {
        system: AromaticSystemId,
        charge: i64,
    },
    #[error("aromatic system {0:?} mixes prescribed and mobile exposed atoms")]
    MixedExposureDemand(AromaticSystemId),
    #[error(
        "node order for aromatic system {system:?} has missing atoms {missing:?} and duplicate atoms {duplicates:?}"
    )]
    InvalidNodeOrder {
        system: AromaticSystemId,
        missing: Vec<AtomId>,
        duplicates: Vec<AtomId>,
    },
    #[error("no perfect matching exists for aromatic system {0:?}")]
    NoMatching(AromaticSystemId),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MatchingInputMode {
    Prescribed,
    OneMobileExposure,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct PrescribedExposure {
    atom: AtomId,
    electrons: i64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MatchingInput {
    required_covered: Vec<AtomId>,
    required_exposed: Vec<PrescribedExposure>,
    exposed_count: usize,
    mode: MatchingInputMode,
}

impl MatchingInput {
    fn from_system(view: AromaticSystemView<'_>) -> Result<Self, KekulizerError> {
        let system = view.id;
        let atoms: Vec<AtomId> = view.atom_ids().collect();
        let ElectronCountsAst::Lit(electrons) = view.electrons() else {
            return Err(KekulizerError::UndeterminedElectrons(system));
        };
        if electrons.len() != atoms.len() {
            return Err(KekulizerError::ElectronCountMismatch {
                system,
                member_count: atoms.len(),
                electron_count: electrons.len(),
            });
        }

        let ValueAst::Lit(charge) = view.charge() else {
            return Err(KekulizerError::UndeterminedCharge(system));
        };
        let (ValueAst::Lit(unpaired), ValueAst::Lit(multiplicity)) =
            (&view.spin().unpaired, &view.spin().multiplicity)
        else {
            return Err(KekulizerError::UndeterminedSpin(system));
        };
        if (*unpaired, *multiplicity) != (0, 1) {
            return Err(KekulizerError::OpenShell(system));
        }
        if charge.abs() > 1 {
            return Err(KekulizerError::UnsupportedCharge {
                system,
                charge: *charge,
            });
        }

        let mut required_covered = Vec::with_capacity(atoms.len());
        let mut required_exposed = Vec::new();
        for (&atom, &electrons) in atoms.iter().zip(electrons) {
            match electrons {
                1 => required_covered.push(atom),
                0 | 2 => required_exposed.push(PrescribedExposure { atom, electrons }),
                contribution => {
                    return Err(KekulizerError::UnsupportedElectronContribution {
                        system,
                        atom,
                        contribution,
                    });
                }
            }
        }

        if required_exposed.len() > 1 {
            return Err(KekulizerError::MultiplePrescribedExposures {
                system,
                count: required_exposed.len(),
            });
        }
        if !required_exposed.is_empty() && *charge != 0 {
            return Err(KekulizerError::MixedExposureDemand(system));
        }

        if required_exposed.is_empty() && *charge != 0 {
            Ok(Self {
                required_covered: Vec::new(),
                required_exposed,
                exposed_count: 1,
                mode: MatchingInputMode::OneMobileExposure,
            })
        } else {
            Ok(Self {
                exposed_count: required_exposed.len(),
                required_covered,
                required_exposed,
                mode: MatchingInputMode::Prescribed,
            })
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KekulizationModel {
    pub algorithm: PerfectMatchingAlgorithm,
}

impl KekulizationModel {
    pub fn new(algorithm: PerfectMatchingAlgorithm) -> Self {
        Self { algorithm }
    }
}

impl Default for KekulizationModel {
    fn default() -> Self {
        Self::new(PerfectMatchingAlgorithm::BacktrackingDfs)
    }
}

#[derive(Clone, Debug)]
pub struct Kekulizer {
    model: KekulizationModel,
    node_order: Vec<AtomId>,
}

impl Kekulizer {
    pub fn new(model: KekulizationModel, node_order: Vec<AtomId>) -> Self {
        Self { model, node_order }
    }
}

impl Transformer for Kekulizer {
    type Error = KekulizerError;

    fn transform_into(&self, ast: &mut MoleculeAst) -> Result<(), KekulizerError> {
        if ast.aromatic_systems().count() == 0 {
            return Ok(());
        }

        // Plan the per-system matching against an immutable AST snapshot, then
        // apply the bond-order writes and structural cleanup in passes that
        // require &mut.
        let plans = self.plan_systems(ast)?;

        // Pass 1: bond-order writes and Aromatic-constraint stripping.
        for plan in &plans {
            for &bid in &plan.matched_bonds {
                let bond = ast.bond_mut(bid).ast;
                bond.order = ValueAst::Lit(2);
                bond.constraints.remove(BondConstraintKey::Aromatic);
            }
            for &bid in &plan.unmatched_bonds {
                let bond = ast.bond_mut(bid).ast;
                bond.order = ValueAst::Lit(1);
                bond.constraints.remove(BondConstraintKey::Aromatic);
            }
            for &aidx in &plan.atoms {
                let atom = ast.atom_mut(aidx).ast;
                atom.constraints.remove(AtomConstraintKey::AromaticValence);
            }
        }

        // Pass 2: drop the aromatic system entries via the builder.
        let to_remove: Vec<AromaticSystemId> = plans.iter().map(|p| p.system_idx).collect();
        let mut builder = ast.edit();
        builder.remove_aromatic_systems(&to_remove);
        *ast = builder.build();

        Ok(())
    }

    fn generate_all<'a>(
        &'a self,
        ast: &'a MoleculeAst,
    ) -> Box<dyn Iterator<Item = MoleculeAst> + 'a> {
        Box::new(self.transform(ast).ok().into_iter())
    }
}

struct SystemPlan {
    system_idx: AromaticSystemId,
    atoms: Vec<AtomId>,
    matched_bonds: Vec<BondId>,
    unmatched_bonds: Vec<BondId>,
}

impl Kekulizer {
    /// Build the per-system matching plan against an immutable AST.
    fn plan_systems(&self, ast: &MoleculeAst) -> Result<Vec<SystemPlan>, KekulizerError> {
        let mut plans = Vec::with_capacity(ast.aromatic_systems().count());
        for view in ast.aromatic_systems().iter() {
            let system_idx = view.id;
            let _matching_input = MatchingInput::from_system(view)?;
            let system_atoms: Vec<AtomId> = view.atom_ids().collect();
            let bonds: Vec<BondId> = view.bond_ids().collect();
            let atom_set: HashSet<AtomId> = system_atoms.iter().copied().collect();
            let mut seen = HashSet::with_capacity(system_atoms.len());
            let mut ordered_host_atoms = Vec::with_capacity(system_atoms.len());
            let mut duplicates = Vec::new();
            for atom in self
                .node_order
                .iter()
                .copied()
                .filter(|atom| atom_set.contains(atom))
            {
                if seen.insert(atom) {
                    ordered_host_atoms.push(atom);
                } else if !duplicates.contains(&atom) {
                    duplicates.push(atom);
                }
            }
            let missing: Vec<AtomId> = system_atoms
                .iter()
                .copied()
                .filter(|atom| !seen.contains(atom))
                .collect();
            if !missing.is_empty() || !duplicates.is_empty() {
                return Err(KekulizerError::InvalidNodeOrder {
                    system: system_idx,
                    missing,
                    duplicates,
                });
            }

            let correspondence = view.induced_subgraph();
            let extracted = ast.extract(&correspondence);
            let sub_order: Vec<AtomId> = ordered_host_atoms
                .iter()
                .map(|&host| {
                    AtomId::from(
                        correspondence
                            .atoms()
                            .left_of(NodeId::from(host))
                            .expect("system atom maps to the extracted molecule"),
                    )
                })
                .collect();
            let matched: HashSet<BondId> = extracted
                .graph()
                .perfect_matching(&sub_order, self.model.algorithm)
                .ok_or(KekulizerError::NoMatching(system_idx))?
                .bonds()
                .map(|sub| {
                    correspondence
                        .bonds()
                        .right_of(sub)
                        .expect("matched bond maps to the host molecule")
                })
                .collect();
            let (matched_bonds, unmatched_bonds): (Vec<BondId>, Vec<BondId>) =
                bonds.iter().copied().partition(|b| matched.contains(b));

            plans.push(SystemPlan {
                system_idx,
                atoms: system_atoms,
                matched_bonds,
                unmatched_bonds,
            });
        }
        Ok(plans)
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{AromaticSystemId, AtomId, MoleculeAst};
    use umol_ast::{mol_dsl, mol_dsl_ground};

    use super::*;

    #[rstest]
    #[case::neutral_all_covered(
        mol_dsl_ground!(r#"{:atoms ["C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2 3] :type "[1,1,1,1]"}]}"#),
        MatchingInput {
            required_covered: vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
            required_exposed: vec![],
            exposed_count: 0,
            mode: MatchingInputMode::Prescribed,
        }
    )]
    #[case::prescribed_donor(
        mol_dsl_ground!(r#"{:atoms ["N" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :type "[2,1,1]"}]}"#),
        MatchingInput {
            required_covered: vec![AtomId(1), AtomId(2)],
            required_exposed: vec![PrescribedExposure { atom: AtomId(0), electrons: 2 }],
            exposed_count: 1,
            mode: MatchingInputMode::Prescribed,
        }
    )]
    #[case::prescribed_acceptor_at_positional_atom(
        mol_dsl_ground!(r#"{:atoms ["C" "B" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :type "[1,0,1]"}]}"#),
        MatchingInput {
            required_covered: vec![AtomId(0), AtomId(2)],
            required_exposed: vec![PrescribedExposure { atom: AtomId(1), electrons: 0 }],
            exposed_count: 1,
            mode: MatchingInputMode::Prescribed,
        }
    )]
    #[case::mobile_positive_charge(
        mol_dsl_ground!(r#"{:atoms ["C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :type "[1,1,1]#c+"}]}"#),
        MatchingInput {
            required_covered: vec![],
            required_exposed: vec![],
            exposed_count: 1,
            mode: MatchingInputMode::OneMobileExposure,
        }
    )]
    #[case::mobile_negative_charge(
        mol_dsl_ground!(r#"{:atoms ["C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :type "[1,1,1]#c-"}]}"#),
        MatchingInput {
            required_covered: vec![],
            required_exposed: vec![],
            exposed_count: 1,
            mode: MatchingInputMode::OneMobileExposure,
        }
    )]
    fn test_matching_input_from_system(
        #[case] input: MoleculeAst,
        #[case] expected: MatchingInput,
    ) {
        let system = input.aromatic_systems().iter().next().unwrap();
        assert_eq!(MatchingInput::from_system(system), Ok(expected));
    }

    #[rstest]
    #[case::undetermined_electrons(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "*#c0#u0#s1"}]}"#),
        KekulizerError::UndeterminedElectrons(AromaticSystemId(0))
    )]
    #[case::undetermined_charge(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "[1,1]#u0#s1"}]}"#),
        KekulizerError::UndeterminedCharge(AromaticSystemId(0))
    )]
    #[case::undetermined_spin(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "[1,1]#c0#u0"}]}"#),
        KekulizerError::UndeterminedSpin(AromaticSystemId(0))
    )]
    #[case::electron_count_mismatch(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "[1]#c0#u0#s1"}]}"#),
        KekulizerError::ElectronCountMismatch {
            system: AromaticSystemId(0),
            member_count: 2,
            electron_count: 1,
        }
    )]
    #[case::unsupported_contribution_at_positional_atom(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "[1,3]#c0#u0#s1"}]}"#),
        KekulizerError::UnsupportedElectronContribution {
            system: AromaticSystemId(0),
            atom: AtomId(1),
            contribution: 3,
        }
    )]
    #[case::open_shell(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "[1,1]#c0#u1#s2"}]}"#),
        KekulizerError::OpenShell(AromaticSystemId(0))
    )]
    #[case::multiple_prescribed_exposures(
        mol_dsl!(r#"{:atoms ["B" "C" "N"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :type "[0,1,2]#c0#u0#s1"}]}"#),
        KekulizerError::MultiplePrescribedExposures {
            system: AromaticSystemId(0),
            count: 2,
        }
    )]
    #[case::unsupported_charge(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "[1,1]#c+2#u0#s1"}]}"#),
        KekulizerError::UnsupportedCharge {
            system: AromaticSystemId(0),
            charge: 2,
        }
    )]
    #[case::mixed_prescribed_and_mobile_demand(
        mol_dsl!(r#"{:atoms ["B" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :type "[0,1]#c+#u0#s1"}]}"#),
        KekulizerError::MixedExposureDemand(AromaticSystemId(0))
    )]
    fn test_matching_input_from_system_error(
        #[case] input: MoleculeAst,
        #[case] expected: KekulizerError,
    ) {
        let system = input.aromatic_systems().iter().next().unwrap();
        assert_eq!(MatchingInput::from_system(system), Err(expected));
    }

    #[rstest]
    #[case::benzene(
        mol_dsl_ground!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]}"#),
        (0..6).map(AtomId).collect(),
        mol_dsl_ground!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#)
    )]
    #[case::embedded_benzene_nonidentity_correspondence(
        mol_dsl_ground!(r#"{:atoms ["O" "H" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :single] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 7 :aromatic] [2 7 :aromatic]] :aromatic-systems [{:atoms [2 3 4 5 6 7] :type "[1,1,1,1,1,1]"}]}"#),
        vec![AtomId(0), AtomId(2), AtomId(3), AtomId(4), AtomId(5), AtomId(6), AtomId(7), AtomId(1)],
        mol_dsl_ground!(r#"{:atoms ["O" "H" "C" "C" "C" "C" "C" "C"] :bonds [[0 1 :single] [2 3 :double] [3 4 :single] [4 5 :double] [5 6 :single] [6 7 :double] [2 7 :single]]}"#)
    )]
    #[case::multiple_disjoint_systems(
        mol_dsl_ground!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 7 :aromatic] [7 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3] :type "[1,1,1,1]"} {:atoms [4 5 6 7] :type "[1,1,1,1]"}]}"#),
        vec![AtomId(0), AtomId(4), AtomId(1), AtomId(5), AtomId(2), AtomId(6), AtomId(3), AtomId(7)],
        mol_dsl_ground!(r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 0 :single] [4 5 :double] [5 6 :single] [6 7 :double] [7 4 :single]]}"#)
    )]
    fn test_kekulizer_transform_into(
        #[case] input: MoleculeAst,
        #[case] node_order: Vec<AtomId>,
        #[case] expected: MoleculeAst,
    ) {
        let mut ast = input;
        Kekulizer::new(KekulizationModel::default(), node_order)
            .transform_into(&mut ast)
            .unwrap();
        assert_eq!(ast, expected);
    }

    #[rstest]
    #[case::kekule_benzene( mol_dsl_ground!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#))]
    fn test_kekulizer_transform_into_identity(#[case] input: MoleculeAst) {
        let mut ast = input.clone();
        Kekulizer::new(KekulizationModel::default(), (0..6).map(AtomId).collect())
            .transform_into(&mut ast)
            .unwrap();
        assert_eq!(ast, input);
    }

    #[rstest]
    #[case::no_matching(
        mol_dsl_ground!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [0 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :type "[1,1,1,1,1]"}]}"#),
        (0..5).map(AtomId).collect(),
        KekulizerError::NoMatching(AromaticSystemId(0))
    )]
    #[case::missing_system_atom(
        mol_dsl_ground!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3] :type "[1,1,1,1]"}]}"#),
        vec![AtomId(0), AtomId(1), AtomId(2)],
        KekulizerError::InvalidNodeOrder {
            system: AromaticSystemId(0),
            missing: vec![AtomId(3)],
            duplicates: vec![],
        }
    )]
    #[case::duplicate_system_atom(
        mol_dsl_ground!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3] :type "[1,1,1,1]"}]}"#),
        vec![AtomId(0), AtomId(1), AtomId(1), AtomId(2), AtomId(3)],
        KekulizerError::InvalidNodeOrder {
            system: AromaticSystemId(0),
            missing: vec![],
            duplicates: vec![AtomId(1)],
        }
    )]
    fn test_kekulizer_transform_into_error(
        #[case] input: MoleculeAst,
        #[case] node_order: Vec<AtomId>,
        #[case] expected: KekulizerError,
    ) {
        let mut ast = input;
        let result =
            Kekulizer::new(KekulizationModel::default(), node_order).transform_into(&mut ast);
        assert_eq!(result, Err(expected));
    }
}
