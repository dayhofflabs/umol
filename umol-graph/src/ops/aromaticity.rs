//! Aromaticity perception primitive.
//!
//! [`AromaticityPerceiver`] dispatches to one of three algorithms (Hückel
//! rule, HMO, Clar) selected by [`AromaticityModel`] and runs perception
//! against an IR. It is the shared core used by three top-level entities:
//! the resolver (validates `#a` hints filled in by atom-typing), the
//! aromatizer (discovers aromatic systems from a Kekulé bond-order layout),
//! and the validator (verifies pre-existing aromatic systems against the
//! model). [`AromaticityPerceiver::derive`] is the standard IR-facing
//! operation; [`AromaticityPerceiver::find_systems`] remains available when
//! the caller supplies another per-atom electron source.
//!
//! System insertion and bond marking are exposed via
//! [`AromaticityPerceiver::add_systems`].

pub mod clar;
pub mod hmo;
pub mod hueckel;

use std::collections::BTreeSet;

pub use clar::{ClarAromaticity, ClarError};
pub use hmo::{HmoAromaticity, HmoError, HmoOutput};
pub use hueckel::HueckelAromaticity;
use thiserror::Error;
use umol_graph_core::{ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm};
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemId, AromaticValence, AromaticValenceForm, AsLit,
    AtomConstraintForm, AtomConstraintKey, AtomHandle, AtomId, BondConstraintForm, BondHandle,
    BondId, BondUpdate, BooleanForm, Edits, ElectronCountsForm, Molecule, MoleculeApplyError,
    RingConfig, RingId, RingModel, RingSet, RingSetKind,
};
use umol_utils::solution::Solution;

use crate::ops::model::{AromaticityModel, AromaticityRule};

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
    Apply(#[from] MoleculeApplyError),
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
    pub systems: Vec<(Vec<AtomId>, AromaticSystemForm)>,
    /// Constraint failures, entity failures, and independently valid mismatches.
    pub inconsistencies: Vec<AromaticityInconsistency>,
}

/// Policy-free classification of aromatic constraint and entity inconsistencies.
#[derive(Clone, Copy, Debug, Error, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AromaticityInconsistency {
    #[error("aromatic valence at atom {atom:?} cannot produce a valid aromatic system")]
    AromaticValenceFailure { atom: AtomId },
    #[error("aromatic constraint at bond {bond:?} cannot produce a valid aromatic system")]
    AromaticBondConstraintFailure { bond: BondId },
    #[error("aromatic system {system:?} is not realizable under the selected model")]
    AromaticSystemFailure { system: AromaticSystemId },
    #[error("aromatic valence at atom {atom:?} disagrees with aromatic system {system:?}")]
    AromaticValenceMismatch {
        atom: AtomId,
        system: AromaticSystemId,
    },
    #[error(
        "aromatic valence at atom {atom:?} disagrees with aromatic constraint at bond {bond:?}"
    )]
    AromaticValenceBondMismatch { atom: AtomId, bond: BondId },
    #[error("aromatic constraint at bond {bond:?} disagrees with aromatic system {system:?}")]
    AromaticBondConstraintMismatch {
        bond: BondId,
        system: AromaticSystemId,
    },
}

#[derive(Clone, Debug)]
pub enum AromaticityPerceiver {
    Hueckel(HueckelAromaticity),
    Hmo(HmoAromaticity),
    Clar(ClarAromaticity),
}

impl AromaticityPerceiver {
    pub fn new(model: &AromaticityModel) -> Self {
        match &model.rule {
            AromaticityRule::Hueckel { ring_limits } => Self::Hueckel(HueckelAromaticity::new(
                model.scope.clone(),
                ring_limits.clone(),
            )),
            AromaticityRule::Hmo {
                stabilization_threshold,
            } => Self::Hmo(HmoAromaticity::new(
                model.scope.clone(),
                *stabilization_threshold,
            )),
            AromaticityRule::Clar => Self::Clar(ClarAromaticity),
        }
    }

    /// The ring set `find_systems` perceives over — the rule's ring request
    /// under the configured algorithms.
    pub(crate) fn candidate_rings(
        &self,
        molecule: &Molecule,
        config: AromaticityConfig,
    ) -> RingSet {
        molecule
            .rings(self.ring_request(), config.ring_config)
            .into_ring_set()
    }

    /// Whether some total contribution reachable from the members' ranges is
    /// perceived as aromatic by the configured rule; rules without a usable
    /// bound accept every range.
    pub(crate) fn accepts_range(&self, members: &[(u32, u32)]) -> bool {
        match self {
            Self::Hueckel(m) => m.accepts_range(members),
            Self::Hmo(m) => m.accepts_range(members),
            Self::Clar(m) => m.accepts_range(members),
        }
    }

    /// The member sets a system could claim: every candidate ring, and —
    /// under the Hückel rule — every ring union within the ring limits.
    /// Member lists are sorted ascending.
    pub(crate) fn claim_candidates(&self, rings: &RingSet) -> Vec<Vec<AtomId>> {
        let mut candidates: Vec<Vec<AtomId>> = rings
            .iter()
            .map(|ring| {
                let mut atoms = ring.atoms().to_vec();
                atoms.sort_unstable();
                atoms
            })
            .collect();
        if let Self::Hueckel(m) = self {
            let eligible: Vec<RingId> = rings.ids().collect();
            for union in m.enumerate_unions(rings, &eligible) {
                let mut atoms: Vec<AtomId> = union.into_iter().collect();
                atoms.sort_unstable();
                candidates.push(atoms);
            }
        }
        candidates
    }

    /// Find candidate aromatic systems via the configured algorithm.
    /// The closure returns each atom's aromatic membership and π contribution.
    /// NotAromatic excludes an atom; Aromatic(0) remains eligible.
    /// No returned system contains both endpoints of a non-aromatic bond assertion.
    /// Conflicting combined candidates are rejected without an alternative partition search.
    #[allow(clippy::complexity)]
    pub fn find_systems<F>(
        &self,
        molecule: &Molecule,
        config: AromaticityConfig,
        electrons_at: F,
    ) -> Result<
        Solution<Vec<(Vec<AtomId>, AromaticSystemForm)>, AromaticityContradiction>,
        AromaticityError,
    >
    where
        F: Fn(AtomId) -> AromaticValence,
    {
        let rings = self.candidate_rings(molecule, config);
        self.find_systems_from_rings(molecule, &rings, config, electrons_at)
    }

    /// Find candidate aromatic systems from an existing candidate ring set.
    #[allow(clippy::complexity)]
    pub(crate) fn find_systems_from_rings<F>(
        &self,
        molecule: &Molecule,
        rings: &RingSet,
        config: AromaticityConfig,
        electrons_at: F,
    ) -> Result<
        Solution<Vec<(Vec<AtomId>, AromaticSystemForm)>, AromaticityContradiction>,
        AromaticityError,
    >
    where
        F: Fn(AtomId) -> AromaticValence,
    {
        self.find_systems_with_exclusions(
            molecule,
            rings,
            config,
            &bond_exclusions(molecule),
            electrons_at,
        )
    }

    #[allow(clippy::type_complexity)]
    fn find_systems_with_exclusions<F>(
        &self,
        molecule: &Molecule,
        rings: &RingSet,
        config: AromaticityConfig,
        excluded_bonds: &[[AtomId; 2]],
        electrons_at: F,
    ) -> Result<
        Solution<Vec<(Vec<AtomId>, AromaticSystemForm)>, AromaticityContradiction>,
        AromaticityError,
    >
    where
        F: Fn(AtomId) -> AromaticValence,
    {
        let systems = match self {
            Self::Hueckel(m) => {
                m.find_with_exclusions(molecule, rings, excluded_bonds, &electrons_at)
            }
            Self::Hmo(m) => match m.find_with_exclusions(
                molecule,
                rings,
                config.connected_components_algorithm,
                excluded_bonds,
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
            Self::Clar(m) => match m.find_with_exclusions(
                molecule,
                rings,
                config.maximum_independent_set_algorithm,
                excluded_bonds,
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

    /// Perceive aromatic systems from complete aromatic-valence evidence, independently assess
    /// stored systems from their electron contributions, and classify their relationship.
    /// Contradictory evidence is reported before perception; nonliteral aromatic counts leave
    /// perception underdetermined. An atom without aromatic evidence is excluded.
    pub fn derive(
        &self,
        molecule: &Molecule,
        config: AromaticityConfig,
    ) -> Result<Solution<AromaticityDerivation, AromaticityContradiction>, AromaticityError> {
        let mut contributions = Vec::with_capacity(molecule.atoms().count());
        let mut inconsistencies = BTreeSet::new();
        for system in molecule.aromatic_systems().iter() {
            if let ElectronCountsForm::Lit(electrons) = &system.attributes().electrons {
                if electrons.len() != system.atom_ids().len() {
                    inconsistencies.insert(AromaticityInconsistency::AromaticSystemFailure {
                        system: system.id(),
                    });
                }
            }
        }
        if !inconsistencies.is_empty() {
            return Ok(Solution::Determined(AromaticityDerivation {
                systems: Vec::new(),
                inconsistencies: inconsistencies.into_iter().collect(),
            }));
        }
        let mut underdetermined = false;
        for atom in molecule.atoms().ids() {
            match read_aromatic_valence(molecule, atom) {
                Ok(Some(value)) => contributions.push(value),
                Ok(None) => underdetermined = true,
                Err(inconsistency) => {
                    inconsistencies.insert(inconsistency);
                }
            }
        }
        if !inconsistencies.is_empty() {
            return Ok(Solution::Determined(AromaticityDerivation {
                systems: Vec::new(),
                inconsistencies: inconsistencies.into_iter().collect(),
            }));
        }
        if underdetermined
            || molecule.aromatic_systems().iter().any(|system| {
                matches!(
                    system.attributes().electrons,
                    ElectronCountsForm::Undetermined
                )
            })
        {
            return Ok(Solution::Underdetermined(AromaticityDerivation::default()));
        }

        let rings = self.candidate_rings(molecule, config);
        match self.check_stored_systems(molecule, &rings, config)? {
            Solution::Determined(existing) => inconsistencies.extend(existing),
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(AromaticityDerivation::default()));
            }
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction))
            }
        }
        if inconsistencies.iter().any(|inconsistency| {
            matches!(
                inconsistency,
                AromaticityInconsistency::AromaticBondConstraintMismatch { .. }
            )
        }) {
            return Ok(Solution::Determined(AromaticityDerivation {
                systems: Vec::new(),
                inconsistencies: inconsistencies.into_iter().collect(),
            }));
        }

        let systems = match self
            .find_systems_from_rings(molecule, &rings, config, |atom| contributions[atom.index()])?
        {
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
        for atom in molecule.atoms().iter() {
            if matches!(
                atom.attributes().constraints.aromatic_valence(),
                AromaticValenceForm::Aromatic(_)
            ) && !accepted_atoms.contains(&atom.id())
            {
                inconsistencies
                    .insert(AromaticityInconsistency::AromaticValenceFailure { atom: atom.id() });
            }
        }
        for bond in molecule.bonds().iter() {
            if matches!(
                bond.attributes().constraints.aromatic(),
                BooleanForm::Lit(true)
            ) && !system_members
                .iter()
                .any(|members| bond.atom_ids().iter().all(|atom| members.contains(atom)))
            {
                inconsistencies.insert(AromaticityInconsistency::AromaticBondConstraintFailure {
                    bond: bond.id(),
                });
            }
        }

        Ok(Solution::Determined(AromaticityDerivation {
            systems,
            inconsistencies: inconsistencies.into_iter().collect(),
        }))
    }

    #[allow(clippy::type_complexity)]
    pub(super) fn check_stored_systems(
        &self,
        molecule: &Molecule,
        rings: &RingSet,
        config: AromaticityConfig,
    ) -> Result<Solution<Vec<AromaticityInconsistency>, AromaticityContradiction>, AromaticityError>
    {
        let mut inconsistencies = BTreeSet::new();
        let mut valid_existing = Vec::new();
        for existing in molecule.aromatic_systems().iter() {
            let ElectronCountsForm::Lit(existing_electrons) = &existing.attributes().electrons
            else {
                return Ok(Solution::Underdetermined(Vec::new()));
            };
            let existing_atoms: Vec<AtomId> = existing.atom_ids().collect();
            if existing_electrons.len() != existing_atoms.len() {
                inconsistencies.insert(AromaticityInconsistency::AromaticSystemFailure {
                    system: existing.id(),
                });
                continue;
            }
            let existing_contributions: Vec<(AtomId, i64)> = existing_atoms
                .iter()
                .copied()
                .zip(existing_electrons.iter().copied())
                .collect();

            let perceived =
                match self.find_systems_with_exclusions(molecule, rings, config, &[], |atom| {
                    existing_contributions
                        .iter()
                        .find_map(|&(candidate, electrons)| {
                            (candidate == atom).then_some(AromaticValence::Aromatic(electrons))
                        })
                        .unwrap_or(AromaticValence::NotAromatic)
                })? {
                    Solution::Determined(perceived) => perceived,
                    Solution::Underdetermined(_) => {
                        return Ok(Solution::Underdetermined(Vec::new()));
                    }
                    Solution::Contradictory(_) => {
                        inconsistencies.insert(AromaticityInconsistency::AromaticSystemFailure {
                            system: existing.id(),
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
                let ElectronCountsForm::Lit(perceived_electrons) = &perceived_system.electrons
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
                valid_existing.push(existing.id());
            } else {
                inconsistencies.insert(AromaticityInconsistency::AromaticSystemFailure {
                    system: existing.id(),
                });
            }
        }

        for system in valid_existing {
            for bond in molecule.aromatic_system(system).bonds() {
                if matches!(
                    bond.attributes().constraints.aromatic(),
                    BooleanForm::Lit(false)
                ) {
                    inconsistencies.insert(
                        AromaticityInconsistency::AromaticBondConstraintMismatch {
                            bond: bond.id(),
                            system,
                        },
                    );
                }
            }
        }

        Ok(Solution::Determined(inconsistencies.into_iter().collect()))
    }

    /// Add perceived systems to the IR.
    pub fn add_systems(
        &self,
        molecule: &mut Molecule,
        systems: Vec<(Vec<AtomId>, AromaticSystemForm)>,
    ) {
        if systems.is_empty() {
            return;
        }
        let mut edits = Edits::new();
        for (atoms, attributes) in systems {
            let bond_ids = molecule.bonds().induced_ids(&atoms);
            edits.add_aromatic_system(atoms.into_iter().map(AtomHandle::Id).collect(), attributes);
            for id in bond_ids {
                edits.update_bond(
                    BondHandle::Id(id),
                    molecule.bond(id).attributes(),
                    &BondUpdate {
                        constraints: BondConstraintForm::aromatic(true).into(),
                        ..Default::default()
                    },
                );
            }
        }
        molecule
            .transact([edits])
            .expect("perceived aromatic systems preserve molecule integrity");
    }

    fn ring_request(&self) -> RingModel {
        let max_ring_size = match self {
            Self::Hueckel(m) => m.ring_limits.max_ring_size,
            Self::Hmo(_) => 22,
            Self::Clar(_) => 6,
        };
        RingModel {
            kind: RingSetKind::Relevant,
            max_ring_size,
        }
    }
}

fn bond_exclusions(molecule: &Molecule) -> Vec<[AtomId; 2]> {
    molecule
        .bonds()
        .iter()
        .filter(|bond| {
            matches!(
                bond.attributes().constraints.aromatic(),
                BooleanForm::Lit(false)
            )
        })
        .map(|bond| bond.atom_ids())
        .collect()
}

/// Read complete aromatic evidence; a nonliteral count is None and a contradiction is Err.
pub(super) fn read_aromatic_valence(
    molecule: &Molecule,
    atom: AtomId,
) -> Result<Option<AromaticValence>, AromaticityInconsistency> {
    let view = molecule.atom(atom);
    let constraints = view.constraints();
    let key = AtomConstraintKey::AromaticValence;
    match constraints.evidenced_complete(key) {
        Some(AtomConstraintForm::AromaticValence(value)) => Ok(value.as_lit()),
        None => {
            if constraints.asserted(key).is_none() {
                return Err(AromaticityInconsistency::AromaticValenceFailure { atom });
            }
            if let Some(system) = view.aromatic_system_id() {
                return Err(AromaticityInconsistency::AromaticValenceMismatch { atom, system });
            }
            if let Some(neighbor) = view.neighbors().find(|neighbor| {
                matches!(
                    neighbor.bond().constraints().aromatic(),
                    BooleanForm::Lit(true)
                )
            }) {
                return Err(AromaticityInconsistency::AromaticValenceBondMismatch {
                    atom,
                    bond: neighbor.bond_id(),
                });
            }
            Err(AromaticityInconsistency::AromaticValenceFailure { atom })
        }
        Some(_) => unreachable!("aromatic-valence reading preserves its key"),
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::{RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm};
    use umol_graph_ir::ir::{
        AromaticSystemId, AromaticValenceForm, AtomConstraintForm, AtomConstraintKey, AtomForm,
        AtomId, BondConstraintKey, BondForm, ElectronCountsForm, Molecule, MoleculeEntries,
        NumForm, UnpairedElectronsForm,
    };
    use umol_graph_ir::{mol_dsl, mol_dsl_concrete};

    use super::*;
    use crate::ops::model::{AromaticityRule, AromaticityTieBreak, ElementScope, RingLimits};

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

    #[rustfmt::skip]
    #[rstest]
    #[case::hueckel_rings_and_union(
        AromaticityRule::Hueckel { ring_limits: RingLimits::default() },
        vec![
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(5)],
            vec![AtomId(4), AtomId(5), AtomId(6), AtomId(7), AtomId(8), AtomId(9)],
            (0..10).map(AtomId).collect::<Vec<_>>(),
        ]
    )]
    #[case::hmo_rings_only(
        AromaticityRule::Hmo { stabilization_threshold: 0.0 },
        vec![
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(5)],
            vec![AtomId(4), AtomId(5), AtomId(6), AtomId(7), AtomId(8), AtomId(9)],
        ]
    )]
    #[case::clar_rings_only(
        AromaticityRule::Clar,
        vec![
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(5)],
            vec![AtomId(4), AtomId(5), AtomId(6), AtomId(7), AtomId(8), AtomId(9)],
        ]
    )]
    fn test_aromaticity_perceiver_claim_candidates(
        #[case] rule: AromaticityRule,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        let molecule = mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]}"#);
        let perception = AromaticityPerceiver::new(&AromaticityModel {
            scope: ElementScope::Any,
            rule,
            tie_break: AromaticityTieBreak::Strict,
        });
        let rings = perception.candidate_rings(&molecule, AromaticityConfig::default());
        assert_eq!(perception.claim_candidates(&rings), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hueckel_fused_merged(
        AromaticityRule::Hueckel { ring_limits: RingLimits::default() },
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]}"#),
        vec![(0..10).map(AtomId).collect::<Vec<_>>()]
    )]
    #[case::hueckel_coupled_ordered(
        AromaticityRule::Hueckel { ring_limits: RingLimits::default() },
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [5 6 "1"]
                    [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 6 "1"]]}"#),
        vec![
            (0..6).map(AtomId).collect::<Vec<_>>(),
            (6..12).map(AtomId).collect::<Vec<_>>(),
        ]
    )]
    #[case::hmo_ring(
        AromaticityRule::Hmo { stabilization_threshold: 0.0 },
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        vec![(0..6).map(AtomId).collect::<Vec<_>>()]
    )]
    #[case::clar_sextet(
        AromaticityRule::Clar,
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]}"#),
        vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(5)]]
    )]
    fn test_aromaticity_perceiver_find_systems_decomposition(
        #[case] rule: AromaticityRule,
        #[case] molecule: Molecule,
        #[case] expected: Vec<Vec<AtomId>>,
    ) {
        // The selection relies on the perception contract: one decomposition
        // per input — pairwise-disjoint systems with sorted member lists, in
        // a deterministic order.
        let perception = AromaticityPerceiver::new(&AromaticityModel {
            scope: ElementScope::Any,
            rule,
            tie_break: AromaticityTieBreak::Strict,
        });
        let electrons = |_: AtomId| AromaticValence::Aromatic(1);
        let config = AromaticityConfig::default();
        let rings = perception.candidate_rings(&molecule, config);
        let Solution::Determined(first) = perception
            .find_systems_from_rings(&molecule, &rings, config, electrons)
            .unwrap()
        else {
            panic!("perception did not determine");
        };
        let members: Vec<Vec<AtomId>> = first.iter().map(|(atoms, _)| atoms.clone()).collect();
        assert_eq!(members, expected);
        for (index, (atoms, _)) in first.iter().enumerate() {
            assert!(atoms.windows(2).all(|pair| pair[0] < pair[1]));
            for (other_index, (other, _)) in first.iter().enumerate() {
                if index != other_index {
                    assert!(atoms.iter().all(|atom| !other.contains(atom)));
                }
            }
        }
        let Solution::Determined(second) = perception
            .find_systems(&molecule, config, electrons)
            .unwrap()
        else {
            panic!("perception did not determine");
        };
        assert_eq!(first, second);
    }

    #[rstest]
    #[case::excluded_ring(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        vec![],
    )]
    #[case::peripheral_exclusion(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]}"#),
        vec![vec![AtomId(4), AtomId(5), AtomId(6), AtomId(7), AtomId(8), AtomId(9)]],
    )]
    #[case::separate_systems(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 6 "1"]
                    [0 6 "1#a!"]]}"#),
        vec![(0..6).map(AtomId).collect(), (6..12).map(AtomId).collect()],
    )]
    #[case::excluded_chord(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [0 3 "1#a!"]]}"#),
        vec![],
    )]
    #[case::coronene_hub_exclusion_with_separate_benzene(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C"
                    "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C"
                    "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 12 "1"]
                    [12 13 "1"] [13 14 "1"] [14 15 "1"] [15 16 "1"] [16 17 "1"] [17 18 "1"]
                    [18 19 "1"] [19 20 "1"] [20 21 "1"] [21 22 "1"] [22 23 "1"] [23 6 "1"]
                    [0 6 "1"] [1 9 "1"] [2 12 "1"] [3 15 "1"] [4 18 "1"] [5 21 "1"]
                    [24 25 "1"] [25 26 "1"] [26 27 "1"] [27 28 "1"] [28 29 "1"] [29 24 "1"]]}"#),
        vec![(24..30).map(AtomId).collect()],
    )]
    fn test_aromaticity_perceiver_find_systems_constraints(
        #[case] molecule: Molecule,
        #[case] members: Vec<Vec<AtomId>>,
        #[values(
            AromaticityRule::Hueckel { ring_limits: RingLimits::default() },
            AromaticityRule::Hmo { stabilization_threshold: 0.1 },
            AromaticityRule::Clar,
        )]
        rule: AromaticityRule,
    ) {
        let expected: Vec<_> = members
            .into_iter()
            .map(|atoms| {
                let form = AromaticSystemForm {
                    electrons: ElectronCountsForm::Lit(vec![1; atoms.len()]),
                    charge: NumForm::Lit(0),
                    unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                    constraints: Default::default(),
                };
                (atoms, form)
            })
            .collect();
        let perception = AromaticityPerceiver::new(&AromaticityModel {
            scope: ElementScope::Any,
            rule,
            tie_break: AromaticityTieBreak::Strict,
        });
        let config = AromaticityConfig::default();
        let electrons_at = |_| AromaticValence::Aromatic(1);
        assert_eq!(
            perception.find_systems(&molecule, config, electrons_at),
            Ok(Solution::Determined(expected.clone())),
        );
        let rings = perception.candidate_rings(&molecule, config);
        let mut direct = match perception {
            AromaticityPerceiver::Hueckel(rule) => {
                rule.find_from_rings(&molecule, &rings, &electrons_at)
            }
            AromaticityPerceiver::Hmo(rule) => rule
                .find_from_rings(
                    &molecule,
                    &rings,
                    config.connected_components_algorithm,
                    &electrons_at,
                )
                .unwrap(),
            AromaticityPerceiver::Clar(rule) => rule
                .find_from_rings(
                    &molecule,
                    &rings,
                    config.maximum_independent_set_algorithm,
                    &electrons_at,
                )
                .unwrap(),
        };
        direct.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(direct, expected);
    }

    fn any_hueckel() -> AromaticityPerceiver {
        AromaticityPerceiver::new(&AromaticityModel {
            scope: ElementScope::Any,
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
            tie_break: AromaticityTieBreak::Strict,
        })
    }

    fn aromatic_valence_lit(molecule: &Molecule, id: AtomId) -> Option<i64> {
        match molecule
            .atom(id)
            .attributes()
            .constraints
            .get(AtomConstraintKey::AromaticValence)?
        {
            AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(n))) => {
                Some(*n)
            }
            _ => None,
        }
    }

    fn aromatic(element: Element, pi: i64) -> AtomForm {
        let mut atom = AtomForm::from_element(element);
        atom.charge = NumForm::Lit(0);
        atom.unpaired_electrons = UnpairedElectronsForm::closed_shell();
        atom.constraints.set(AtomConstraintForm::AromaticValence(
            AromaticValenceForm::Aromatic(NumForm::Lit(pi)),
        ));
        atom
    }

    fn benzene() -> Molecule {
        let atoms: Vec<AtomForm> = (0..6).map(|_| aromatic(Element::C, 1)).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondForm::from_order(1)))
            .collect();
        Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn pyrrole() -> Molecule {
        let atoms = vec![
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ];
        let bonds: Vec<_> = (0..5)
            .map(|i| (AtomId(i), AtomId((i + 1) % 5), BondForm::from_order(1)))
            .collect();
        Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn run_full(
        perception: &AromaticityPerceiver,
        molecule: &mut Molecule,
    ) -> Solution<(), AromaticityContradiction> {
        let outcome = perception
            .find_systems(molecule, AromaticityConfig::default(), |v| {
                match molecule.atom(v).attributes().constraints.aromatic_valence() {
                    AromaticValenceForm::Aromatic(NumForm::Lit(n)) if *n >= 0 => {
                        AromaticValence::Aromatic(*n)
                    }
                    _ => AromaticValence::NotAromatic,
                }
            })
            .unwrap();
        match outcome {
            Solution::Determined(systems) => {
                perception.add_systems(molecule, systems);
                Solution::Determined(())
            }
            Solution::Underdetermined(_) => Solution::Underdetermined(()),
            Solution::Contradictory(c) => Solution::Contradictory(c),
        }
    }

    #[rstest]
    #[case::missing_evidence(
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#h2"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"]
                    [4 5 "1"] [5 6 "1"] [6 0 "1"]]
        }"#),
        AromaticityDerivation {
            systems: vec![],
            inconsistencies: (0..6).map(|id| {
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(id) }
            }).collect(),
        },
    )]
    #[case::zero_contribution(
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#c+#a0"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"]
                    [4 5 "1"] [5 6 "1"] [6 0 "1"]]
        }"#),
        AromaticityDerivation {
            systems: vec![(
                (0..7).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 1, 1, 0])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )],
            inconsistencies: vec![],
        },
    )]
    fn test_aromaticity_perceiver_derive_membership(
        #[case] molecule: Molecule,
        #[case] expected: AromaticityDerivation,
    ) {
        assert_eq!(
            AromaticityPerceiver::new(&AromaticityModel::daylight())
                .derive(&molecule, AromaticityConfig::default()),
            Ok(Solution::Determined(expected)),
        );
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
                AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
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
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(0) },
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(1) },
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(2) },
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(3) },
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(4) },
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
                :attrs "[1,1,1,1,1,1]"
            }]
        }"#),
        Solution::Underdetermined(AromaticityDerivation::default()),
    )]
    #[case::vacuous_projections(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#a*" "C#a*" "C#a*" "C#a*" "C#a*" "C#a*"]
            :bonds [[0 1 "1#a*"] [1 2 "1#a*"] [2 3 "1#a*"] [3 4 "1#a*"]
                    [4 5 "1#a*"] [5 0 "1#a*"]]
            :aromatic-systems [{
                :atoms [0 1 2 3 4 5]
                :attrs "[1,1,1,1,1,1]"
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
                AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
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
                :attrs "[2,0,1,1,1,1]"
            }]
        }"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![],
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
                :attrs "[1,1,1,1,1,1]"
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
                AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )],
            inconsistencies: vec![],
        }),
    )]
    #[case::short_counts(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,2]"}]}"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![],
            inconsistencies: vec![AromaticityInconsistency::AromaticSystemFailure { system: AromaticSystemId(0) }],
        }),
    )]
    #[case::long_counts(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,2,0,5]"}]}"#),
        Solution::Determined(AromaticityDerivation {
            systems: vec![],
            inconsistencies: vec![AromaticityInconsistency::AromaticSystemFailure { system: AromaticSystemId(0) }],
        }),
    )]
    #[case::existing_system_rejected(
        AromaticityModel::mdl(),
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 0 "1#a"]]
            :aromatic-systems [{
                :atoms [0 1 2 3 4]
                :attrs "[2,1,1,1,1]"
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
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(0) },
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(1) },
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(2) },
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(3) },
                AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(4) },
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
    fn test_aromaticity_perceiver_derive(
        #[case] model: AromaticityModel,
        #[case] molecule: Molecule,
        #[case] expected: Solution<AromaticityDerivation, AromaticityContradiction>,
    ) {
        assert_eq!(
            AromaticityPerceiver::new(&model)
                .derive(&molecule, AromaticityConfig::default())
                .unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_aromaticity_perceiver_hueckel_benzene() {
        let perception = AromaticityPerceiver::new(&AromaticityModel {
            scope: ElementScope::AllowList(vec![Element::C]),
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
            tie_break: AromaticityTieBreak::Strict,
        });
        let mut molecule = benzene();
        let solution = run_full(&perception, &mut molecule);
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(molecule.aromatic_systems().count(), 1);
        let system = molecule.aromatic_system(AromaticSystemId(0));
        let atoms: Vec<AtomId> = system.atom_ids().collect();
        assert_eq!(atoms.len(), 6);
        let aromatic_bond_count = molecule
            .bonds()
            .iter()
            .filter(|view| {
                view.attributes()
                    .constraints
                    .contains(BondConstraintKey::Aromatic)
            })
            .count();
        assert_eq!(aromatic_bond_count, 6);
    }

    #[rstest]
    fn test_aromaticity_perceiver_clar_heterocycle() {
        let perception = AromaticityPerceiver::new(&AromaticityModel {
            scope: ElementScope::Any,
            rule: AromaticityRule::Clar,
            tie_break: AromaticityTieBreak::Strict,
        });
        let mut molecule = pyrrole();
        let solution = run_full(&perception, &mut molecule);
        assert!(matches!(
            solution,
            Solution::Contradictory(AromaticityContradiction::ClarNonBenzenoid(_))
        ));
    }

    #[rstest]
    #[case::cyclopropenium_cation(
        mol_dsl_concrete!(r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#),
        0, vec![1, 1, 0], vec![0, 0, 1], vec![1, 1, 0],
    )]
    #[case::cot_dianion(
        mol_dsl_concrete!(r#"{:atoms ["C #h #a" "C #c- #h #a2" "C #h #a" "C #h #a"
                                "C #h #a" "C #c- #h #a2" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"]
                               [4 5 "1"] [5 6 "1"] [6 7 "1"] [7 0 "1"]]}"#),
        0, vec![1, 2, 1, 1, 1, 2, 1, 1],
        vec![0, -1, 0, 0, 0, -1, 0, 0], vec![1, 2, 1, 1, 1, 2, 1, 1],
    )]
    #[case::s4_dication(
        mol_dsl_concrete!(r#"{:atoms ["S #c+ #n1 #a" "S #n1 #a2" "S #c+ #n1 #a" "S #n1 #a2"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]]}"#),
        0, vec![1, 2, 1, 2], vec![1, 0, 1, 0], vec![1, 2, 1, 2],
    )]
    #[case::boratabenzene_anion(
        mol_dsl_concrete!(r#"{:atoms ["B #c- #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        0, vec![1; 6], vec![-1, 0, 0, 0, 0, 0], vec![1; 6],
    )]
    #[case::borepin(
        mol_dsl_concrete!(r#"{:atoms ["B #h #a0" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 0 "1"]]}"#),
        0, vec![0, 1, 1, 1, 1, 1, 1], vec![0; 7], vec![0, 1, 1, 1, 1, 1, 1],
    )]
    #[case::pyridinium(
        mol_dsl_concrete!(r#"{:atoms ["N #c+ #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        0, vec![1; 6], vec![1, 0, 0, 0, 0, 0], vec![1; 6],
    )]
    #[case::pyrylium(
        mol_dsl_concrete!(r#"{:atoms ["O #c+ #n1 #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        0, vec![1; 6], vec![1, 0, 0, 0, 0, 0], vec![1; 6],
    )]
    #[case::pyrrole(
        mol_dsl_concrete!(r#"{:atoms ["N #h #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        0, vec![2, 1, 1, 1, 1], vec![0; 5], vec![2, 1, 1, 1, 1],
    )]
    #[case::furan(
        mol_dsl_concrete!(r#"{:atoms ["O #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        0, vec![2, 1, 1, 1, 1], vec![0; 5], vec![2, 1, 1, 1, 1],
    )]
    #[case::thiophene(
        mol_dsl_concrete!(r#"{:atoms ["S #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        0, vec![2, 1, 1, 1, 1], vec![0; 5], vec![2, 1, 1, 1, 1],
    )]
    fn test_aromaticity_perceiver_add_systems(
        #[case] mut molecule: Molecule,
        #[case] system_charge: i64,
        #[case] electrons: Vec<i64>,
        #[case] atom_charges: Vec<i64>,
        #[case] aromatic_valences: Vec<i64>,
    ) {
        let outcome = any_hueckel()
            .find_systems(&molecule, AromaticityConfig::default(), |v| match molecule
                .atom(v)
                .attributes()
                .constraints
                .aromatic_valence()
            {
                AromaticValenceForm::Aromatic(NumForm::Lit(n)) if *n >= 0 => {
                    AromaticValence::Aromatic(*n)
                }
                _ => AromaticValence::NotAromatic,
            })
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut molecule, systems);

        let system = molecule.aromatic_system(AromaticSystemId(0));
        assert_eq!(system.attributes().charge, NumForm::Lit(system_charge));
        assert_eq!(
            system.attributes().electrons,
            ElectronCountsForm::Lit(electrons)
        );
        for (i, (q, k)) in atom_charges
            .iter()
            .zip(aromatic_valences.iter())
            .enumerate()
        {
            let id = AtomId(i as u32);
            assert_eq!(molecule.atom(id).attributes().charge, NumForm::Lit(*q));
            assert_eq!(aromatic_valence_lit(&molecule, id), Some(*k));
        }
    }

    #[rstest]
    #[case::disjoint_systems(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [4 5 "1"]]}"#),
        vec![
            (vec![AtomId(0), AtomId(1), AtomId(2)], AromaticSystemForm::from_electrons(vec![1, 1, 1])),
            (vec![AtomId(3), AtomId(4), AtomId(5)], AromaticSystemForm::from_electrons(vec![1, 1, 1])),
        ],
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "2#a"] [1 2 "1#a"] [2 3 "1"] [3 4 "2#a"] [4 5 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}
                               {:atoms [3 4 5] :attrs "[1,1,1]"}]}"#),
    )]
    fn test_aromaticity_perceiver_add_systems_entries(
        #[case] mut molecule: Molecule,
        #[case] systems: Vec<(Vec<AtomId>, AromaticSystemForm)>,
        #[case] expected: Molecule,
    ) {
        any_hueckel().add_systems(&mut molecule, systems);
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::empty(Molecule::new())]
    #[case::populated(mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#))]
    fn test_aromaticity_perceiver_add_systems_identity(#[case] molecule: Molecule) {
        let mut actual = molecule.clone();
        any_hueckel().add_systems(&mut actual, Vec::new());
        assert_eq!(actual, molecule);
    }

    #[rstest]
    fn test_aromaticity_perceiver_hueckel_no_aromatic_atom() {
        let perception = AromaticityPerceiver::new(&AromaticityModel {
            scope: ElementScope::AllowList(vec![Element::C]),
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
            tie_break: AromaticityTieBreak::Strict,
        });
        let atoms: Vec<AtomForm> = (0..6).map(|_| AtomForm::from_element(Element::C)).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondForm::from_order(1)))
            .collect();
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        });
        let solution = run_full(&perception, &mut molecule);
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(molecule.aromatic_systems().count(), 0);
        let any_aromatic = molecule.bonds().iter().any(|view| {
            view.attributes()
                .constraints
                .contains(BondConstraintKey::Aromatic)
        });
        assert!(!any_aromatic);
    }
}
