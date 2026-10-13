//! Kekulizer: aromatic-system form → Kekulé bond orders.
//!
//! For each aromatic system, excludes fixed zero- and two-electron contributors
//! and finds a perfect matching on the remaining one-electron contributors.
//! Each two-electron contributor regains one local lone pair. An otherwise
//! all-one-electron system with charge ±1 instead localizes that charge on the
//! single atom left exposed by a maximum matching. Matched bonds become order 2;
//! unmatched bonds become order 1. The aromatic constraint is
//! removed from the system's bonds and from the system's atoms; the system
//! entry itself is removed at the end.
//!
//! The matching algorithm is configurable via [`KekulizeConfig`]; the
//! atom processing order — which controls determinism — is fixed at
//! construction time and is the caller's responsibility (e.g., from a
//! nauty/Traces canonical labeling).

use std::collections::HashSet;
use std::{iter, mem};

use thiserror::Error;
use umol_graph_core::{
    BipartiteMaximumMatchingAlgorithm, GeneralMaximumMatchingAlgorithm, NonBipartiteGraphError,
};
use umol_graph_ir::ir::{
    AromaticSystemHandle, AromaticSystemId, AtomConstraintForm, AtomConstraintKey, AtomForm,
    AtomHandle, AtomId, AtomUpdate, BondConstraintKey, BondHandle, BondId, BondUpdate, Edit, Edits,
    ElectronCountsForm, Lattice, Molecule, MulticenterValenceForm, Normalize, NumForm,
};
use umol_utils::solution::Solution;

use crate::ops::invariant::{ValenceInvariantsValidator, ValenceMismatch};
use crate::ops::transform::Transformer;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum KekulizeError {
    #[error("electron contributions are undetermined for aromatic system {0:?}")]
    UndeterminedElectrons(AromaticSystemId),
    #[error("charge is undetermined for aromatic system {0:?}")]
    UndeterminedCharge(AromaticSystemId),
    #[error("unpaired electrons are undetermined for aromatic system {0:?}")]
    UndeterminedUnpairedElectrons(AromaticSystemId),
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
    #[error(
        "maximum matching for aromatic system {system:?} has deficiency {actual}; expected {expected}"
    )]
    MatchingDeficiency {
        system: AromaticSystemId,
        expected: usize,
        actual: usize,
    },
    #[error("charge is undetermined for exposed atom {atom:?} in aromatic system {system:?}")]
    UndeterminedExposedAtomCharge {
        system: AromaticSystemId,
        atom: AtomId,
    },
    #[error("lone pairs are undetermined for exposed atom {atom:?} in aromatic system {system:?}")]
    UndeterminedExposedAtomLonePairs {
        system: AromaticSystemId,
        atom: AtomId,
    },
    #[error(
        "atom {atom:?} contradicts {constraint:?} in planned aromatic system {aromatic_system:?}"
    )]
    ContradictedAssertion {
        aromatic_system: AromaticSystemId,
        atom: AtomId,
        constraint: AtomConstraintForm,
    },
    #[error("planned aromatic system {aromatic_system:?}: {mismatch}")]
    ValenceInvariant {
        aromatic_system: AromaticSystemId,
        mismatch: ValenceMismatch,
    },
    #[error("no perfect matching exists for aromatic system {0:?}")]
    NoMatching(AromaticSystemId),
    #[error("maximum matching requires a bipartite aromatic system {0:?}")]
    NonBipartiteMatching(AromaticSystemId),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MatchingInputMode {
    Prescribed,
    OneMobileExposure,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ExposedAtom {
    atom: AtomId,
    electrons: i64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MatchingInput {
    required_covered: Vec<AtomId>,
    required_exposed: Vec<ExposedAtom>,
    exposed_count: usize,
    mode: MatchingInputMode,
}

impl MatchingInput {
    fn from_system(molecule: &Molecule, system: AromaticSystemId) -> Result<Self, KekulizeError> {
        let view = molecule.aromatic_system(system);
        let atoms: Vec<AtomId> = view.atom_ids().collect();
        let ElectronCountsForm::Lit(electrons) = view.electrons() else {
            return Err(KekulizeError::UndeterminedElectrons(system));
        };
        if electrons.len() != atoms.len() {
            return Err(KekulizeError::ElectronCountMismatch {
                system,
                member_count: atoms.len(),
                electron_count: electrons.len(),
            });
        }

        let NumForm::Lit(charge) = view.charge() else {
            return Err(KekulizeError::UndeterminedCharge(system));
        };
        let (NumForm::Lit(unpaired_electrons), NumForm::Lit(multiplicity)) = (
            &view.unpaired_electrons().count,
            &view.unpaired_electrons().multiplicity,
        ) else {
            return Err(KekulizeError::UndeterminedUnpairedElectrons(system));
        };
        if (*unpaired_electrons, *multiplicity) != (0, 1) {
            return Err(KekulizeError::OpenShell(system));
        }
        if charge.abs() > 1 {
            return Err(KekulizeError::UnsupportedCharge {
                system,
                charge: *charge,
            });
        }

        let mut required_covered = Vec::with_capacity(atoms.len());
        let mut required_exposed = Vec::new();
        for (&atom, &electrons) in atoms.iter().zip(electrons) {
            match electrons {
                1 => required_covered.push(atom),
                0 | 2 => required_exposed.push(ExposedAtom { atom, electrons }),
                contribution => {
                    return Err(KekulizeError::UnsupportedElectronContribution {
                        system,
                        atom,
                        contribution,
                    });
                }
            }
        }

        if !required_exposed.is_empty() && *charge != 0 {
            return Err(KekulizeError::MixedExposureDemand(system));
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
pub enum MaximumMatchingAlgorithm {
    Edmonds,
    HopcroftKarp,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KekulizeConfig {
    pub algorithm: MaximumMatchingAlgorithm,
}

impl KekulizeConfig {
    pub fn new(algorithm: MaximumMatchingAlgorithm) -> Self {
        Self { algorithm }
    }
}

impl Default for KekulizeConfig {
    fn default() -> Self {
        Self::new(MaximumMatchingAlgorithm::Edmonds)
    }
}

/// Localizes aromatic systems into single/double bonds and atom-local electrons.
///
/// Fixed two-electron donors regain one lone pair; fixed zero-electron acceptors
/// retain their local state. The remaining one-electron contributors must admit
/// a perfect matching. A delocalized system charge of ±1 is supported separately
/// when all contributions are one, by localizing it on one exposed atom.
///
/// Transformations report [`KekulizeError::NoMatching`] for an impossible fixed
/// assignment and [`KekulizeError::UndeterminedExposedAtomLonePairs`] when restoring
/// a donor pair requires a non-literal local count. They do not select another
/// hydrogen placement or contribution assignment. Planned incidence must agree
/// with retained valence assertions, and planned atoms must not contradict
/// electron conservation. Unchanged spin data and unrelated atoms are preserved
/// without whole-molecule validation.
///
/// # Semantic properties
///
/// When the required readings are determined, every localized system atom
/// satisfies electron conservation and its retained valence assertions.
/// All errors occur before mutation: in-place failure leaves the input exactly
/// unchanged; consuming failure drops it. Unrelated undetermined state is retained.
#[derive(Clone, Debug)]
pub struct Kekulizer {
    config: KekulizeConfig,
    node_order: Vec<AtomId>,
}

impl Kekulizer {
    pub fn new(config: KekulizeConfig, node_order: Vec<AtomId>) -> Self {
        Self { config, node_order }
    }
}

impl Transformer for Kekulizer {
    type Error = KekulizeError;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, KekulizeError> {
        if molecule.aromatic_systems().count() == 0 {
            return Ok(molecule);
        }

        let edits = self.plan_transform(&molecule)?;
        let editor = molecule
            .edit()
            .apply(edits)
            .expect("kekulization plan applies to its input");
        Ok(editor
            .finish()
            .expect("kekulization plan preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), KekulizeError> {
        if molecule.aromatic_systems().count() == 0 {
            return Ok(());
        }

        let edits = self.plan_transform(molecule)?;
        *molecule = mem::take(molecule)
            .apply(edits)
            .expect("kekulization plan applies and preserves molecule integrity");
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        iter::once_with(move || self.transform(molecule.clone()).ok()).flatten()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct SystemPlan {
    system_idx: AromaticSystemId,
    atoms: Vec<AtomId>,
    // Follows atoms order; a mobile exposure has only one entry.
    exposed_atoms: Vec<ExposedAtom>,
    mobile_charge: Option<i64>,
    matched_bonds: Vec<BondId>,
    unmatched_bonds: Vec<BondId>,
}

impl Kekulizer {
    fn plan_transform(&self, molecule: &Molecule) -> Result<Edits, KekulizeError> {
        let plans = self.plan_systems(molecule)?;
        let mut edits = Edits::new();
        let mut orders = vec![None; molecule.bonds().count()];
        for plan in &plans {
            debug_assert!(plan.matched_bonds.iter().all(|&bond| {
                molecule.bond(bond).atom_ids().iter().all(|atom| {
                    !plan
                        .exposed_atoms
                        .iter()
                        .any(|exposed| exposed.atom == *atom)
                })
            }));
            for (id, order) in plan
                .matched_bonds
                .iter()
                .map(|&id| (id, 2))
                .chain(plan.unmatched_bonds.iter().map(|&id| (id, 1)))
            {
                orders[id.index()] = Some(order);
                let bond = molecule.bond(id).attributes();
                edits.update_bond(
                    BondHandle::Id(id),
                    bond,
                    &BondUpdate {
                        order: Some(NumForm::Lit(order)),
                        ..Default::default()
                    },
                );
                if let Some(old) = bond.constraints.get(BondConstraintKey::Aromatic) {
                    edits.push(Edit::ModifyBondConstraint {
                        id: BondHandle::Id(id),
                        old: Some(old.clone()),
                        new: None,
                    });
                }
            }
            let mut exposed_atoms = plan.exposed_atoms.iter().peekable();
            for &id in &plan.atoms {
                let atom = molecule.atom(id);
                let attributes = atom.attributes();
                let mut planned = AtomForm {
                    element: attributes.element.clone(),
                    isotope_mass: attributes.isotope_mass.clone(),
                    charge: attributes.charge.clone(),
                    implicit_hydrogens: attributes.implicit_hydrogens.clone(),
                    lone_pairs: attributes.lone_pairs.clone(),
                    unpaired_electrons: attributes.unpaired_electrons.clone(),
                    constraints: Default::default(),
                };
                if let Some(old) = attributes
                    .constraints
                    .get(AtomConstraintKey::AromaticValence)
                {
                    edits.push(Edit::ModifyAtomConstraint {
                        id: AtomHandle::Id(id),
                        old: Some(old.clone()),
                        new: None,
                    });
                }
                if let Some(exposed) = exposed_atoms.next_if(|exposed| exposed.atom == id) {
                    let mut update = AtomUpdate::default();
                    if let Some(system_charge) = plan.mobile_charge {
                        let NumForm::Lit(local_charge) = planned.charge else {
                            return Err(KekulizeError::UndeterminedExposedAtomCharge {
                                system: plan.system_idx,
                                atom: id,
                            });
                        };
                        planned.charge = NumForm::Lit(local_charge + system_charge);
                        update.charge = Some(planned.charge.clone());
                    }
                    if exposed.electrons == 2 {
                        let NumForm::Lit(lone_pairs) = planned.lone_pairs else {
                            return Err(KekulizeError::UndeterminedExposedAtomLonePairs {
                                system: plan.system_idx,
                                atom: id,
                            });
                        };
                        planned.lone_pairs = NumForm::Lit(lone_pairs + 1);
                        update.lone_pairs = Some(planned.lone_pairs.clone());
                    }
                    edits.update_atom(AtomHandle::Id(id), atom.attributes(), &update);
                }
                let mut valence = NumForm::Lit(0);
                for neighbor in atom.neighbors() {
                    let order = match orders[neighbor.bond_id().index()] {
                        Some(order) => NumForm::Lit(order),
                        None => neighbor.bond().order().clone().normalize().map_err(|_| {
                            KekulizeError::ValenceInvariant {
                                aromatic_system: plan.system_idx,
                                mismatch: ValenceMismatch::ContradictoryConstraint {
                                    atom_id: id,
                                    key: AtomConstraintKey::Valence,
                                },
                            }
                        })?,
                    };
                    valence = valence + order;
                }
                let mut determined_incidence = true;
                let has_dative = molecule.dative_bonds().has_incident(id);
                for key in [
                    AtomConstraintKey::Valence,
                    AtomConstraintKey::MulticenterValence,
                    AtomConstraintKey::DonatedPairs,
                    AtomConstraintKey::AcceptedPairs,
                ] {
                    let assertion = attributes.constraints.get(key);
                    let derived = match key {
                        AtomConstraintKey::Valence => {
                            Some(AtomConstraintForm::valence(valence.clone()))
                        }
                        AtomConstraintKey::MulticenterValence
                            if assertion.is_none() && !atom.is_in_multicenter_bond() =>
                        {
                            continue;
                        }
                        AtomConstraintKey::DonatedPairs | AtomConstraintKey::AcceptedPairs
                            if assertion.is_none() && !has_dative =>
                        {
                            continue;
                        }
                        _ => atom.constraints().derived_complete(key),
                    }
                    .ok_or(KekulizeError::ValenceInvariant {
                        aromatic_system: plan.system_idx,
                        mismatch: ValenceMismatch::ContradictoryConstraint { atom_id: id, key },
                    })?;
                    let value = match assertion {
                        None => derived,
                        Some(constraint) => attributes
                            .constraints
                            .asserted(key)
                            .and_then(|asserted| asserted.meet(&derived))
                            .ok_or_else(|| KekulizeError::ContradictedAssertion {
                                aromatic_system: plan.system_idx,
                                atom: id,
                                constraint: constraint.clone(),
                            })?,
                    };
                    determined_incidence &= !value.is_undetermined();
                    match value {
                        AtomConstraintForm::DonatedPairs(NumForm::Lit(0))
                        | AtomConstraintForm::AcceptedPairs(NumForm::Lit(0))
                        | AtomConstraintForm::MulticenterValence(
                            MulticenterValenceForm::NotMulticenter,
                        ) => {}
                        _ => planned.constraints.set(value),
                    }
                }
                // Bare-atom validation defaults absent incidence to zero.
                if determined_incidence {
                    if let Solution::Contradictory(mismatch) = ValenceInvariantsValidator
                        .validate_atom(&planned)
                        .unwrap_or_else(|error| match error {})
                    {
                        let mismatch = match mismatch {
                            ValenceMismatch::OrbitalCount {
                                orbital_count,
                                electron_count,
                                ..
                            } => ValenceMismatch::OrbitalCount {
                                atom_id: id,
                                orbital_count,
                                electron_count,
                            },
                            ValenceMismatch::ContradictoryConstraint { key, .. } => {
                                ValenceMismatch::ContradictoryConstraint { atom_id: id, key }
                            }
                        };
                        return Err(KekulizeError::ValenceInvariant {
                            aromatic_system: plan.system_idx,
                            mismatch,
                        });
                    }
                }
            }
            debug_assert!(exposed_atoms.next().is_none());
        }
        if !plans.is_empty() {
            edits.remove_aromatic_systems(
                plans
                    .into_iter()
                    .map(|plan| {
                        (
                            AromaticSystemHandle::Id(plan.system_idx),
                            plan.atoms.into_iter().map(AtomHandle::Id).collect(),
                            molecule
                                .aromatic_system(plan.system_idx)
                                .attributes()
                                .clone(),
                        )
                    })
                    .collect(),
            );
        }
        Ok(edits)
    }

    /// Build the per-system matching plan against an immutable molecule.
    fn plan_systems(&self, molecule: &Molecule) -> Result<Vec<SystemPlan>, KekulizeError> {
        let mut plans = Vec::with_capacity(molecule.aromatic_systems().count());
        for view in molecule.aromatic_systems().iter() {
            let system_idx = view.id();
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
                return Err(KekulizeError::InvalidNodeOrder {
                    system: system_idx,
                    missing,
                    duplicates,
                });
            }

            let matching_input = MatchingInput::from_system(molecule, view.id())?;
            let mobile_charge = match (matching_input.mode, view.charge()) {
                (MatchingInputMode::OneMobileExposure, NumForm::Lit(charge)) => Some(*charge),
                (MatchingInputMode::Prescribed, _) => None,
                (MatchingInputMode::OneMobileExposure, _) => {
                    unreachable!("matching input requires a literal system charge")
                }
            };
            let prescribed_exposed_set: HashSet<AtomId> = matching_input
                .required_exposed
                .iter()
                .map(|exposure| exposure.atom)
                .collect();
            let matching_host_atoms: Vec<AtomId> = match matching_input.mode {
                MatchingInputMode::Prescribed => ordered_host_atoms
                    .iter()
                    .copied()
                    .filter(|atom| !prescribed_exposed_set.contains(atom))
                    .collect(),
                MatchingInputMode::OneMobileExposure => ordered_host_atoms.clone(),
            };

            let correspondence_host_atoms: Vec<AtomId> = system_atoms
                .iter()
                .copied()
                .filter(|atom| !prescribed_exposed_set.contains(atom))
                .collect();
            let correspondence = molecule.induced_subgraph(&correspondence_host_atoms);
            let extracted = molecule.extract(&correspondence);
            let sub_order: Vec<AtomId> = matching_host_atoms
                .iter()
                .map(|&host| {
                    correspondence
                        .atoms()
                        .left_of(host)
                        .expect("system atom maps to the extracted molecule")
                })
                .collect();
            let matching = match self.config.algorithm {
                MaximumMatchingAlgorithm::Edmonds => extracted
                    .graph()
                    .general_maximum_matching(&sub_order, GeneralMaximumMatchingAlgorithm::Edmonds),
                MaximumMatchingAlgorithm::HopcroftKarp => extracted
                    .graph()
                    .bipartite_maximum_matching(
                        &sub_order,
                        BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
                    )
                    .map_err(|NonBipartiteGraphError| {
                        KekulizeError::NonBipartiteMatching(system_idx)
                    })?,
            };
            let deficiency = extracted.atoms().count() - 2 * matching.size();
            match matching_input.mode {
                MatchingInputMode::Prescribed if deficiency != 0 => {
                    return Err(KekulizeError::NoMatching(system_idx));
                }
                MatchingInputMode::OneMobileExposure
                    if deficiency != matching_input.exposed_count =>
                {
                    return Err(KekulizeError::MatchingDeficiency {
                        system: system_idx,
                        expected: matching_input.exposed_count,
                        actual: deficiency,
                    });
                }
                _ => {}
            }
            let exposed_atoms: Vec<ExposedAtom> = match matching_input.mode {
                MatchingInputMode::Prescribed => matching_input.required_exposed,
                MatchingInputMode::OneMobileExposure => extracted
                    .atoms()
                    .ids()
                    .filter(|&sub| !matching.is_matched(sub))
                    .map(|sub| ExposedAtom {
                        atom: correspondence
                            .atoms()
                            .right_of(sub)
                            .expect("extracted atom maps to the host molecule"),
                        electrons: if mobile_charge == Some(-1) { 2 } else { 0 },
                    })
                    .collect(),
            };
            let matched: HashSet<BondId> = matching
                .bonds()
                .map(|sub| {
                    correspondence
                        .bonds()
                        .right_of(sub)
                        .expect("matched bond maps to the host molecule")
                })
                .collect();
            let matched_atoms: HashSet<AtomId> = matched
                .iter()
                .flat_map(|&bond| molecule.bond(bond).atom_ids())
                .collect();
            let exposed_set: HashSet<AtomId> =
                exposed_atoms.iter().map(|exposed| exposed.atom).collect();
            let actual_exposed: HashSet<AtomId> = system_atoms
                .iter()
                .copied()
                .filter(|atom| !matched_atoms.contains(atom))
                .collect();
            let required_covered: HashSet<AtomId> = match matching_input.mode {
                MatchingInputMode::Prescribed => {
                    matching_input.required_covered.iter().copied().collect()
                }
                MatchingInputMode::OneMobileExposure => system_atoms
                    .iter()
                    .copied()
                    .filter(|atom| !exposed_set.contains(atom))
                    .collect(),
            };
            if exposed_atoms.len() != matching_input.exposed_count
                || matched_atoms != required_covered
                || actual_exposed != exposed_set
            {
                return Err(KekulizeError::NoMatching(system_idx));
            }
            let (matched_bonds, unmatched_bonds): (Vec<BondId>, Vec<BondId>) =
                bonds.iter().copied().partition(|b| matched.contains(b));

            plans.push(SystemPlan {
                system_idx,
                atoms: system_atoms,
                exposed_atoms,
                mobile_charge,
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
    use umol_graph_ir::ir::{AromaticSystemId, AtomId, Molecule};
    use umol_graph_ir::{mol_dsl, mol_dsl_concrete};

    use super::*;

    #[rstest]
    #[case::neutral_all_covered(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2 3] :attrs "[1,1,1,1]"}]}"#),
        MatchingInput {
            required_covered: vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
            required_exposed: vec![],
            exposed_count: 0,
            mode: MatchingInputMode::Prescribed,
        }
    )]
    #[case::prescribed_donor(
        mol_dsl_concrete!(r#"{:atoms ["N" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[2,1,1]"}]}"#),
        MatchingInput {
            required_covered: vec![AtomId(1), AtomId(2)],
            required_exposed: vec![ExposedAtom { atom: AtomId(0), electrons: 2 }],
            exposed_count: 1,
            mode: MatchingInputMode::Prescribed,
        }
    )]
    #[case::prescribed_acceptor_at_positional_atom(
        mol_dsl_concrete!(r#"{:atoms ["C" "B" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,0,1]"}]}"#),
        MatchingInput {
            required_covered: vec![AtomId(0), AtomId(2)],
            required_exposed: vec![ExposedAtom { atom: AtomId(1), electrons: 0 }],
            exposed_count: 1,
            mode: MatchingInputMode::Prescribed,
        }
    )]
    #[case::multiple_prescribed_exposures(
        mol_dsl_concrete!(r#"{:atoms ["B" "C" "N"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[0,1,2]"}]}"#),
        MatchingInput {
            required_covered: vec![AtomId(1)],
            required_exposed: vec![ExposedAtom { atom: AtomId(0), electrons: 0 }, ExposedAtom { atom: AtomId(2), electrons: 2 }],
            exposed_count: 2,
            mode: MatchingInputMode::Prescribed,
        }
    )]
    #[case::mobile_positive_charge(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]#c+"}]}"#),
        MatchingInput {
            required_covered: vec![],
            required_exposed: vec![],
            exposed_count: 1,
            mode: MatchingInputMode::OneMobileExposure,
        }
    )]
    #[case::mobile_negative_charge(
        mol_dsl_concrete!(r#"{:atoms ["C" "C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]#c-"}]}"#),
        MatchingInput {
            required_covered: vec![],
            required_exposed: vec![],
            exposed_count: 1,
            mode: MatchingInputMode::OneMobileExposure,
        }
    )]
    fn test_matching_input_from_system(#[case] input: Molecule, #[case] expected: MatchingInput) {
        let system = input.aromatic_systems().ids().next().unwrap();
        assert_eq!(MatchingInput::from_system(&input, system), Ok(expected));
    }

    #[rstest]
    #[case::short_counts(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,2]"}]}"#),
        KekulizeError::ElectronCountMismatch { system: AromaticSystemId(0), member_count: 3, electron_count: 2 }
    )]
    #[case::long_counts(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,2,0,5]"}]}"#),
        KekulizeError::ElectronCountMismatch { system: AromaticSystemId(0), member_count: 3, electron_count: 4 }
    )]
    #[case::undetermined_electrons(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "*#c0#u0#s1"}]}"#),
        KekulizeError::UndeterminedElectrons(AromaticSystemId(0))
    )]
    #[case::undetermined_charge(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#u0#s1"}]}"#),
        KekulizeError::UndeterminedCharge(AromaticSystemId(0))
    )]
    #[case::undetermined_unpaired_electrons(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#c0#u0"}]}"#),
        KekulizeError::UndeterminedUnpairedElectrons(AromaticSystemId(0))
    )]
    #[case::unsupported_contribution_at_positional_atom(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,3]#c0#u0#s1"}]}"#),
        KekulizeError::UnsupportedElectronContribution {
            system: AromaticSystemId(0),
            atom: AtomId(1),
            contribution: 3,
        }
    )]
    #[case::open_shell(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#c0#u1#s2"}]}"#),
        KekulizeError::OpenShell(AromaticSystemId(0))
    )]
    #[case::unsupported_charge(
        mol_dsl!(r#"{:atoms ["C" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[1,1]#c+2#u0#s1"}]}"#),
        KekulizeError::UnsupportedCharge {
            system: AromaticSystemId(0),
            charge: 2,
        }
    )]
    #[case::mixed_prescribed_and_mobile_demand(
        mol_dsl!(r#"{:atoms ["B" "C"] :bonds [] :aromatic-systems [{:atoms [0 1] :attrs "[0,1]#c+#u0#s1"}]}"#),
        KekulizeError::MixedExposureDemand(AromaticSystemId(0))
    )]
    fn test_matching_input_from_system_error(
        #[case] input: Molecule,
        #[case] expected: KekulizeError,
    ) {
        let system = input.aromatic_systems().ids().next().unwrap();
        assert_eq!(MatchingInput::from_system(&input, system), Err(expected));
    }

    #[rstest]
    #[case::edmonds(
        MaximumMatchingAlgorithm::Edmonds,
        KekulizeConfig {
            algorithm: MaximumMatchingAlgorithm::Edmonds,
        }
    )]
    #[case::hopcroft_karp(
        MaximumMatchingAlgorithm::HopcroftKarp,
        KekulizeConfig {
            algorithm: MaximumMatchingAlgorithm::HopcroftKarp,
        }
    )]
    fn test_kekulize_config_new(
        #[case] algorithm: MaximumMatchingAlgorithm,
        #[case] expected: KekulizeConfig,
    ) {
        assert_eq!(KekulizeConfig::new(algorithm), expected);
    }

    #[rstest]
    fn test_kekulize_config_default() {
        assert_eq!(
            KekulizeConfig::default(),
            KekulizeConfig {
                algorithm: MaximumMatchingAlgorithm::Edmonds,
            }
        );
    }

    #[rstest]
    #[case::benzene(
        mol_dsl_concrete!(r#"{:atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
        (0..6).map(AtomId).collect(),
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#)
    )]
    #[case::embedded_benzene_nonidentity_correspondence(
        mol_dsl_concrete!(r#"{:atoms ["O#h#n2" "H" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"] :bonds [[0 1 :single] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 7 :aromatic] [2 7 :aromatic]] :aromatic-systems [{:atoms [2 3 4 5 6 7] :attrs "[1,1,1,1,1,1]"}]}"#),
        vec![AtomId(0), AtomId(2), AtomId(3), AtomId(4), AtomId(5), AtomId(6), AtomId(7), AtomId(1)],
        mol_dsl_concrete!(r#"{:atoms ["O#h#n2" "H" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"] :bonds [[0 1 :single] [2 3 :double] [3 4 :single] [4 5 :double] [5 6 :single] [6 7 :double] [2 7 :single]]}"#)
    )]
    #[case::multiple_disjoint_systems(
        mol_dsl_concrete!(r#"{:atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 7 :aromatic] [7 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3] :attrs "[1,1,1,1]"} {:atoms [4 5 6 7] :attrs "[1,1,1,1]"}]}"#),
        vec![AtomId(0), AtomId(4), AtomId(1), AtomId(5), AtomId(2), AtomId(6), AtomId(3), AtomId(7)],
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 0 :single] [4 5 :double] [5 6 :single] [6 7 :double] [7 4 :single]]}"#)
    )]
    #[case::unchanged_spin(
        mol_dsl_concrete!(r#"{:atoms ["N#a#u2#s2" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
        (0..6).map(AtomId).collect(),
        mol_dsl_concrete!(r#"{:atoms ["N#u2#s2" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"]]}"#)
    )]
    #[case::unknown_external_order(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "*"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic] [0 6 "*"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
        (0..7).map(AtomId).collect(),
        mol_dsl_concrete!(r#"{:atoms ["C" "C#h" "C#h" "C#h" "C#h" "C#h" "*"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"] [0 6 "*"]]}"#)
    )]
    #[case::retained_incidence(
        mol_dsl_concrete!(r#"{:atoms ["N#a#v3#d1" "C#h#a" "C#h#a" "N#a#m2" "C#h#a" "C#h#a" "B#h3" "B#h3"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
            :dative-bonds [{:donors [0] :acceptor 6 :attrs "1"}]
            :multicenter-bonds [{:atoms [3 7] :attrs "[2,0]"}]}"#),
        (0..8).map(AtomId).collect(),
        mol_dsl_concrete!(r#"{:atoms ["N#v3#d1" "C#h" "C#h" "N#m2" "C#h" "C#h" "B#h3" "B#h3"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"]]
            :dative-bonds [{:donors [0] :acceptor 6 :attrs "1"}]
            :multicenter-bonds [{:atoms [3 7] :attrs "[2,0]"}]}"#)
    )]
    #[case::normalized_external_order(
        {
            let mut input = mol_dsl_concrete!(r#"{:atoms ["C#a#v4" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "F#n3"]
                :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic] [0 6 "1"]]
                :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#);
            input.bond_mut(BondId(6)).attributes_mut().order = NumForm::lit_set([1]);
            input
        },
        (0..7).map(AtomId).collect(),
        {
            let mut expected = mol_dsl_concrete!(r#"{:atoms ["C#v4" "C#h" "C#h" "C#h" "C#h" "C#h" "F#n3"]
                :bonds [[0 1 "2"] [1 2 "1"] [2 3 "2"] [3 4 "1"] [4 5 "2"] [0 5 "1"] [0 6 "1"]]}"#);
            expected.bond_mut(BondId(6)).attributes_mut().order = NumForm::lit_set([1]);
            expected
        }
    )]
    fn test_kekulizer_transform_into(
        #[case] input: Molecule,
        #[case] node_order: Vec<AtomId>,
        #[case] expected: Molecule,
    ) {
        let transformer = Kekulizer::new(KekulizeConfig::default(), node_order);
        assert_eq!(transformer.transform(input.clone()), Ok(expected.clone()));
        let mut molecule = input;
        transformer.transform_into(&mut molecule).unwrap();
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::kekule_benzene( mol_dsl_concrete!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#))]
    fn test_kekulizer_transform_into_identity(#[case] input: Molecule) {
        let transformer = Kekulizer::new(KekulizeConfig::default(), (0..6).map(AtomId).collect());
        assert_eq!(transformer.transform(input.clone()), Ok(input.clone()));
        let mut molecule = input.clone();
        transformer.transform_into(&mut molecule).unwrap();
        assert_eq!(molecule, input);
    }

    #[rstest]
    #[case::no_matching(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [0 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]"}]}"#),
        (0..5).map(AtomId).collect(),
        KekulizeError::NoMatching(AromaticSystemId(0))
    )]
    #[case::missing_system_atom(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3] :attrs "[1,1,1,1]"}]}"#),
        vec![AtomId(0), AtomId(1), AtomId(2)],
        KekulizeError::InvalidNodeOrder {
            system: AromaticSystemId(0),
            missing: vec![AtomId(3)],
            duplicates: vec![],
        }
    )]
    #[case::duplicate_system_atom(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3] :attrs "[1,1,1,1]"}]}"#),
        vec![AtomId(0), AtomId(1), AtomId(1), AtomId(2), AtomId(3)],
        KekulizeError::InvalidNodeOrder {
            system: AromaticSystemId(0),
            missing: vec![],
            duplicates: vec![AtomId(1)],
        }
    )]
    #[case::contradicted_valence(
        mol_dsl_concrete!(r#"{:atoms ["C#h#a" "C#h#a" "C#h#a#v2" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
        (0..6).map(AtomId).collect(),
        KekulizeError::ContradictedAssertion {
            aromatic_system: AromaticSystemId(0),
            atom: AtomId(2),
            constraint: AtomConstraintForm::valence(2),
        }
    )]
    #[case::later_system_mismatch(
        mol_dsl_concrete!(r#"{:atoms ["*#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"
                                    "C#h#a" "C#h#a" "C#h#a#n" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]
                    [6 7 :aromatic] [7 8 :aromatic] [8 9 :aromatic] [9 10 :aromatic] [10 11 :aromatic] [6 11 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}
                               {:atoms [6 7 8 9 10 11] :attrs "[1,1,1,1,1,1]"}]}"#),
        (0..12).map(AtomId).collect(),
        KekulizeError::ValenceInvariant {
            aromatic_system: AromaticSystemId(1),
            mismatch: ValenceMismatch::OrbitalCount {
                atom_id: AtomId(8), orbital_count: 10, electron_count: 8,
            },
        }
    )]
    #[case::contradicted_donation(
        mol_dsl_concrete!(r#"{:atoms ["C#h#a#d1" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
        (0..6).map(AtomId).collect(),
        KekulizeError::ContradictedAssertion {
            aromatic_system: AromaticSystemId(0),
            atom: AtomId(0),
            constraint: AtomConstraintForm::DonatedPairs(NumForm::Lit(1)),
        }
    )]
    #[case::contradictory_external_order(
        {
            let mut input = mol_dsl_concrete!(r#"{:atoms ["C#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "F#n3"]
                :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic] [0 6 "1"]]
                :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#);
            input.bond_mut(BondId(6)).attributes_mut().order = NumForm::lit_set([]);
            input
        },
        (0..7).map(AtomId).collect(),
        KekulizeError::ValenceInvariant {
            aromatic_system: AromaticSystemId(0),
            mismatch: ValenceMismatch::ContradictoryConstraint {
                atom_id: AtomId(0), key: AtomConstraintKey::Valence,
            },
        }
    )]
    fn test_kekulizer_transform_into_error(
        #[case] input: Molecule,
        #[case] node_order: Vec<AtomId>,
        #[case] expected: KekulizeError,
    ) {
        let transformer = Kekulizer::new(KekulizeConfig::default(), node_order);
        assert_eq!(transformer.transform(input.clone()), Err(expected.clone()));
        let mut molecule = input.clone();
        let result = transformer.transform_into(&mut molecule);
        assert_eq!(result, Err(expected));
        assert_eq!(molecule, input);
    }

    #[rstest]
    #[case::non_bipartite(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [0 2 :aromatic]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]}"#),
        (0..3).map(AtomId).collect(),
        KekulizeError::NonBipartiteMatching(AromaticSystemId(0)),
    )]
    #[case::mobile_exposure(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]}"#),
        (0..5).map(AtomId).collect(),
        KekulizeError::NonBipartiteMatching(AromaticSystemId(0)),
    )]
    fn test_kekulizer_transform_into_matching_error(
        #[case] input: Molecule,
        #[case] node_order: Vec<AtomId>,
        #[case] expected: KekulizeError,
    ) {
        let transformer = Kekulizer::new(
            KekulizeConfig::new(MaximumMatchingAlgorithm::HopcroftKarp),
            node_order,
        );
        assert_eq!(transformer.transform(input.clone()), Err(expected.clone()));
        let mut actual = input.clone();
        let result = transformer.transform_into(&mut actual);

        assert_eq!(result, Err(expected));
        assert_eq!(actual, input);
    }

    #[rstest]
    #[case::benzene(
        mol_dsl_concrete!(r#"{:atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [0 5 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"] :bonds [[0 1 :double] [1 2 :single] [2 3 :double] [3 4 :single] [4 5 :double] [0 5 :single]]}"#)
    )]
    fn test_kekulizer_transform_iter(#[case] input: Molecule, #[case] expected: Molecule) {
        let original = input.clone();
        let transformer = Kekulizer::new(KekulizeConfig::default(), input.atoms().ids().collect());
        let mut results = transformer.transform_iter(&input);
        assert_eq!(results.next(), Some(expected));
        assert_eq!(results.next(), None);
        assert_eq!(results.next(), None);
        assert_eq!(input, original);
    }

    #[rstest]
    #[case::no_matching(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [0 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]"}]}"#),
        KekulizeError::NoMatching(AromaticSystemId(0))
    )]
    fn test_kekulizer_transform_iter_error(
        #[case] input: Molecule,
        #[case] expected: KekulizeError,
    ) {
        let original = input.clone();
        let transformer = Kekulizer::new(KekulizeConfig::default(), input.atoms().ids().collect());
        {
            let mut results = transformer.transform_iter(&input);
            assert_eq!(results.next(), None);
            assert_eq!(results.next(), None);
        }
        assert_eq!(input, original);
        assert_eq!(transformer.transform(input), Err(expected));
    }

    #[rstest]
    #[case::five_membered_ring(
        mol_dsl_concrete!(r#"{:atoms ["N#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
        (0..5).map(AtomId).collect(),
        vec![SystemPlan {
            system_idx: AromaticSystemId(0),
            atoms: (0..5).map(AtomId).collect(),
            exposed_atoms: vec![ExposedAtom { atom: AtomId(0), electrons: 2 }],
            mobile_charge: None,
            matched_bonds: vec![BondId(1), BondId(3)],
            unmatched_bonds: vec![BondId(0), BondId(4), BondId(2)],
        }]
    )]
    #[case::seven_membered_ring(
        mol_dsl_concrete!(r#"{:atoms ["B#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[0,1,1,1,1,1,1]"}]}"#),
        (0..7).map(AtomId).collect(),
        vec![SystemPlan {
            system_idx: AromaticSystemId(0),
            atoms: (0..7).map(AtomId).collect(),
            exposed_atoms: vec![ExposedAtom { atom: AtomId(0), electrons: 0 }],
            mobile_charge: None,
            matched_bonds: vec![BondId(1), BondId(3), BondId(5)],
            unmatched_bonds: vec![BondId(0), BondId(6), BondId(2), BondId(4)],
        }]
    )]
    #[case::fused_heterocycle(
        mol_dsl_concrete!(r#"{:atoms ["N#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 7 :aromatic] [7 8 :aromatic] [8 3 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5 6 7 8] :attrs "[2,1,1,1,1,1,1,1,1]"}]}"#),
        (0..9).map(AtomId).collect(),
        vec![SystemPlan {
            system_idx: AromaticSystemId(0),
            atoms: (0..9).map(AtomId).collect(),
            exposed_atoms: vec![ExposedAtom { atom: AtomId(0), electrons: 2 }],
            mobile_charge: None,
            matched_bonds: vec![BondId(1), BondId(3), BondId(6), BondId(8)],
            unmatched_bonds: vec![BondId(0), BondId(4), BondId(2), BondId(9), BondId(5), BondId(7)],
        }]
    )]
    #[case::nontrivial_host_ids(
        mol_dsl_concrete!(r#"{:atoms ["O" "H" "N#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :single] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 2 :aromatic]] :aromatic-systems [{:atoms [2 3 4 5 6] :attrs "[2,1,1,1,1]"}]}"#),
        vec![AtomId(0), AtomId(2), AtomId(3), AtomId(4), AtomId(5), AtomId(6), AtomId(1)],
        vec![SystemPlan {
            system_idx: AromaticSystemId(0),
            atoms: (2..7).map(AtomId).collect(),
            exposed_atoms: vec![ExposedAtom { atom: AtomId(2), electrons: 2 }],
            mobile_charge: None,
            matched_bonds: vec![BondId(2), BondId(4)],
            unmatched_bonds: vec![BondId(1), BondId(5), BondId(3)],
        }]
    )]
    #[case::five_membered_mobile_hole_natural_order(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]}"#),
        (0..5).map(AtomId).collect(),
        vec![SystemPlan {
            system_idx: AromaticSystemId(0),
            atoms: (0..5).map(AtomId).collect(),
            exposed_atoms: vec![ExposedAtom { atom: AtomId(4), electrons: 2 }],
            mobile_charge: Some(-1),
            matched_bonds: vec![BondId(0), BondId(2)],
            unmatched_bonds: vec![BondId(4), BondId(1), BondId(3)],
        }]
    )]
    #[case::five_membered_mobile_hole_rotated_order(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]}"#),
        vec![AtomId(2), AtomId(3), AtomId(4), AtomId(0), AtomId(1)],
        vec![SystemPlan {
            system_idx: AromaticSystemId(0),
            atoms: (0..5).map(AtomId).collect(),
            exposed_atoms: vec![ExposedAtom { atom: AtomId(0), electrons: 2 }],
            mobile_charge: Some(-1),
            matched_bonds: vec![BondId(1), BondId(3)],
            unmatched_bonds: vec![BondId(0), BondId(4), BondId(2)],
        }]
    )]
    #[case::seven_membered_mobile_hole_natural_order(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[1,1,1,1,1,1,1]#c+"}]}"#),
        (0..7).map(AtomId).collect(),
        vec![SystemPlan {
            system_idx: AromaticSystemId(0),
            atoms: (0..7).map(AtomId).collect(),
            exposed_atoms: vec![ExposedAtom { atom: AtomId(6), electrons: 0 }],
            mobile_charge: Some(1),
            matched_bonds: vec![BondId(0), BondId(2), BondId(4)],
            unmatched_bonds: vec![BondId(6), BondId(1), BondId(3), BondId(5)],
        }]
    )]
    #[case::seven_membered_mobile_hole_rotated_order(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 4 :aromatic] [4 5 :aromatic] [5 6 :aromatic] [6 0 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[1,1,1,1,1,1,1]#c+"}]}"#),
        vec![AtomId(3), AtomId(4), AtomId(5), AtomId(6), AtomId(0), AtomId(1), AtomId(2)],
        vec![SystemPlan {
            system_idx: AromaticSystemId(0),
            atoms: (0..7).map(AtomId).collect(),
            exposed_atoms: vec![ExposedAtom { atom: AtomId(1), electrons: 0 }],
            mobile_charge: Some(1),
            matched_bonds: vec![BondId(6), BondId(2), BondId(4)],
            unmatched_bonds: vec![BondId(0), BondId(1), BondId(3), BondId(5)],
        }]
    )]
    fn test_kekulizer_plan_systems(
        #[case] input: Molecule,
        #[case] node_order: Vec<AtomId>,
        #[case] mut expected: Vec<SystemPlan>,
    ) {
        let kekulizer = Kekulizer::new(KekulizeConfig::default(), node_order);
        let mut actual = kekulizer.plan_systems(&input).unwrap();
        assert_eq!(kekulizer.plan_systems(&input), Ok(actual.clone()));
        for plan in actual.iter_mut().chain(&mut expected) {
            plan.matched_bonds.sort_unstable();
            plan.unmatched_bonds.sort_unstable();
        }
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::impossible_prescribed_hole(
        mol_dsl_concrete!(r#"{:atoms ["N#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic] [3 0 :aromatic] [0 4 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
        (0..5).map(AtomId).collect(),
        KekulizeError::NoMatching(AromaticSystemId(0))
    )]
    #[case::mobile_matching_deficiency(
        mol_dsl_concrete!(r#"{:atoms ["C#a" "C#a" "C#a" "C#a" "C#a"] :bonds [[0 1 :aromatic]] :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]#c-"}]}"#),
        (0..5).map(AtomId).collect(),
        KekulizeError::MatchingDeficiency {
            system: AromaticSystemId(0),
            expected: 1,
            actual: 3,
        }
    )]
    fn test_kekulizer_plan_systems_error(
        #[case] input: Molecule,
        #[case] node_order: Vec<AtomId>,
        #[case] expected: KekulizeError,
    ) {
        let kekulizer = Kekulizer::new(KekulizeConfig::default(), node_order);
        assert_eq!(kekulizer.plan_systems(&input), Err(expected));
    }
}
