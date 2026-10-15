//! Aromaticity resolver. Perception reads aromatic valence constraints from the
//! materialized valence stage; planning emits complete aromatic systems and
//! discharges the assertions they realize without mutating the source molecule.

use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::{mem, slice};

use smallvec::{smallvec, SmallVec};
use thiserror::Error;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemHandle, AromaticSystemId, AromaticValence,
    AromaticValenceForm, AsLit, AtomConstraintForm, AtomForm, AtomHandle, AtomId, AtomUpdate,
    BondConstraintForm, BondConstraintKey, BondHandle, BondId, BondUpdate, BooleanForm, Edit,
    Edits, ElectronCountsForm, Lattice, Molecule, NumForm, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{
    read_aromatic_valence, AromaticityConfig, AromaticityContradiction, AromaticityError,
    AromaticityInconsistency, AromaticityPerceiver,
};
use crate::ops::model::{AromaticityModel, AromaticityTieBreak, ValenceTieBreak};
use crate::ops::resolve::ResolveState;
use crate::ops::valence::compare::compare_by_key;

/// Maximum number of search perception calls per rule-supplied component.
const MAX_PERCEPTION_CALLS: usize = 4096;

/// Search algorithm for joint atom completion and aromatic-system selection.
///
/// Both algorithms apply the same validity checks and tie-breaks. A search needing
/// more than 4096 perception calls in one component is underdetermined.
///
/// # Semantic properties
///
/// When both searches stay within the bound, resolution outcomes and reports agree.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AtomSelectionAlgorithm {
    /// Visit every completion combination, subject to the perception-call bound.
    FullDfs,
    /// Prune certainly invalid subtrees under the Error aromatic-valence policy.
    /// Under Keep, visit the same combinations as FullDfs.
    #[default]
    BacktrackingDfs,
}

/// How aromaticity resolution handles an independently invalid constraint or entity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityFailurePolicy {
    Error,
    Keep,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AromaticityResolveConfig {
    pub perception: AromaticityConfig,
    /// Completion search; defaults to BacktrackingDfs.
    pub atom_selection_algorithm: AtomSelectionAlgorithm,
    pub aromatic_valence_failure: AromaticityFailurePolicy,
    pub aromatic_system_failure: AromaticityFailurePolicy,
}

impl Default for AromaticityResolveConfig {
    fn default() -> Self {
        Self {
            perception: AromaticityConfig::default(),
            atom_selection_algorithm: AtomSelectionAlgorithm::default(),
            aromatic_valence_failure: AromaticityFailurePolicy::Error,
            aromatic_system_failure: AromaticityFailurePolicy::Error,
        }
    }
}

/// Resolves aromatic systems from the molecule's evidence.
///
/// # Semantic properties
///
/// Plans apply to their input and preserve molecule integrity. Consuming and in-place
/// operations agree on outcomes and determined molecules for the same input and configuration.
#[derive(Clone, Debug)]
pub struct AromaticityResolver {
    pub(super) perception: AromaticityPerceiver,
    tie_break: AromaticityTieBreak,
    config: AromaticityResolveConfig,
}

/// Failures to express aromatic systems as atom and bond assertions.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityProjectError {
    #[error("aromatic system {aromatic_system:?} has non-concrete fields")]
    NonConcreteAromaticSystem { aromatic_system: AromaticSystemId },
    #[error("aromatic system {aromatic_system:?} has nonzero charge")]
    ChargedAromaticSystem { aromatic_system: AromaticSystemId },
    #[error("aromatic system {aromatic_system:?} is not closed-shell singlet")]
    AromaticSystemSpin { aromatic_system: AromaticSystemId },
    #[error("atom {atom:?} has an assertion incompatible with its aromatic contribution")]
    AtomAssertion { atom: AtomId },
    #[error("bond {bond:?} has an assertion incompatible with aromaticity")]
    BondAssertion { bond: BondId },
}

impl AromaticityResolver {
    pub fn new(model: &AromaticityModel) -> Self {
        Self::with_config(model, AromaticityResolveConfig::default())
    }

    pub fn with_config(model: &AromaticityModel, config: AromaticityResolveConfig) -> Self {
        Self {
            perception: AromaticityPerceiver::new(model),
            tie_break: model.tie_break,
            config,
        }
    }

    /// Construct the complete aromaticity edit plan without mutating `molecule`.
    /// Newly added systems discharge their members' aromatic-valence assertions.
    /// Outside stored systems, confirmed bond-aromaticity assertions are removed.
    /// Assertions on stored systems and other constraint keys are preserved.
    /// Mismatched evidence is contradictory regardless of the failure policies.
    ///
    /// # Semantic properties
    ///
    /// Stored systems must be reproduced by perception of the whole molecule. Error
    /// rejects an unreproduced system during planning; Keep adds no system in its
    /// search component. Determined plans never introduce overlapping systems.
    pub fn plan_resolve(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<Edits, AromaticityContradiction>, AromaticityError> {
        let outcome = self.select(molecule, ResolveState::default(), ValenceTieBreak::Strict)?;

        match outcome {
            Solution::Determined(state) => {
                let existing: BTreeSet<Vec<AtomId>> = molecule
                    .aromatic_systems()
                    .iter()
                    .map(|aromatic_system| {
                        let mut atoms: Vec<AtomId> = aromatic_system.atom_ids().collect();
                        atoms.sort_unstable();
                        atoms
                    })
                    .collect();

                let mut edits = self.plan_bond_discharge(molecule, &state.aromatic_systems);
                for (atoms, aromatic_system) in state.aromatic_systems {
                    let mut key = atoms.clone();
                    key.sort_unstable();
                    if !existing.contains(&key) {
                        for edit in
                            self.perception
                                .plan_aromatic_system(molecule, atoms, aromatic_system)
                        {
                            edits.push(edit);
                        }
                    }
                }
                Ok(Solution::Determined(edits))
            }
            Solution::Underdetermined(_) => Ok(Solution::Underdetermined(Edits::new())),
            Solution::Contradictory(contradiction) => Ok(Solution::Contradictory(contradiction)),
        }
    }

    /// Resolves the molecule, consuming its storage without a recovery journal.
    ///
    /// Uses the same plan as [AromaticityResolver::resolve_into]. A determined result has checked
    /// molecule integrity. Rejection drops the input.
    ///
    /// # Errors
    ///
    /// Returns an error if the selected aromaticity model lacks required parameters.
    pub fn resolve(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, AromaticityContradiction, ()>, AromaticityError> {
        let edits = match self.plan_resolve(&molecule)? {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("aromaticity resolve plan applies to its input");
        let molecule = editor
            .finish()
            .expect("aromaticity resolve plan preserves molecule integrity");
        Ok(Solution::Determined(molecule))
    }

    /// Plans and applies aromaticity resolution without a recovery journal.
    ///
    /// # Semantic properties
    ///
    /// Planning finishes before mutation. Only Determined changes the molecule;
    /// every other outcome leaves it unchanged.
    ///
    /// # Errors
    ///
    /// Returns an error if the selected aromaticity model lacks required parameters.
    pub fn resolve_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), AromaticityContradiction>, AromaticityError> {
        let edits = match self.plan_resolve(molecule)? {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        if !edits.is_empty() {
            *molecule = mem::take(molecule)
                .apply(edits)
                .expect("aromaticity resolve plan applies and preserves molecule integrity");
        }
        Ok(Solution::Determined(()))
    }

    /// Projects the molecule, consuming its storage without a recovery journal.
    ///
    /// Uses the same plan as [AromaticityResolver::project_into]. A determined result has checked
    /// molecule integrity. Rejection drops the input.
    ///
    /// # Errors
    ///
    /// Returns the planning error described by [AromaticityResolver::plan_project].
    pub fn project(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, AromaticityContradiction, ()>, AromaticityProjectError> {
        let edits = self.plan_project(&molecule)?;
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("aromaticity project plan applies to its input");
        let molecule = editor
            .finish()
            .expect("aromaticity project plan preserves molecule integrity");
        Ok(Solution::Determined(molecule))
    }

    /// Replaces existing aromatic systems with atom contribution and bond assertions.
    ///
    /// Writes each member's electron contribution as `#a<n>` and each system bond as `#a`,
    /// then removes the systems. Atom electron fields and localized bonds are preserved.
    /// This operation does not perceive systems, admit valence candidates, or apply selection
    /// policies. Molecules without systems are unchanged.
    ///
    /// # Semantic properties
    ///
    /// Projection is idempotent and independent of system participant order when electron
    /// counts are transported with their participants. Successful projection publishes all
    /// edits together without a recovery journal; every failure leaves the input unchanged.
    /// System removal compacts constraints that refer to the removed systems.
    ///
    /// # Errors
    ///
    /// Rejects non-concrete system contributions, charge, or spin; nonzero system charge;
    /// non-singlet or nonzero-unpaired system spin; and incompatible atom or bond assertions.
    pub fn project_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), AromaticityContradiction>, AromaticityProjectError> {
        let edits = self.plan_project(molecule)?;
        if !edits.is_empty() {
            *molecule = mem::take(molecule)
                .apply(edits)
                .expect("aromaticity project plan applies and preserves molecule integrity");
        }
        Ok(Solution::Determined(()))
    }

    /// Plans aromatic atom and bond assertions and removal of the aromatic systems.
    ///
    /// Applying the returned edits to the same molecule produces the result of
    /// [AromaticityResolver::project]. The molecule remains unchanged while planning.
    ///
    /// # Errors
    ///
    /// Rejects nonconcrete system fields, nonzero charge, spin other than a closed-shell
    /// singlet, and conflicting assertions.
    pub fn plan_project(&self, molecule: &Molecule) -> Result<Edits, AromaticityProjectError> {
        if !molecule.has_aromatic_systems() {
            return Ok(Edits::new());
        }
        let mut edits = Edits::new();
        for aromatic_system in molecule.aromatic_systems().iter() {
            let AromaticSystemForm {
                electrons: ElectronCountsForm::Lit(electrons),
                charge: NumForm::Lit(charge),
                unpaired_electrons:
                    UnpairedElectronsForm {
                        count: NumForm::Lit(unpaired),
                        multiplicity: NumForm::Lit(multiplicity),
                    },
                constraints: _,
            } = aromatic_system.attributes()
            else {
                return Err(AromaticityProjectError::NonConcreteAromaticSystem {
                    aromatic_system: aromatic_system.id(),
                });
            };
            if *charge != 0 {
                return Err(AromaticityProjectError::ChargedAromaticSystem {
                    aromatic_system: aromatic_system.id(),
                });
            }
            if *unpaired != 0 || *multiplicity != 1 {
                return Err(AromaticityProjectError::AromaticSystemSpin {
                    aromatic_system: aromatic_system.id(),
                });
            }
            for (atom, &electrons) in aromatic_system.atom_ids().zip(electrons) {
                let attributes = molecule.atom(atom).attributes();
                let aromatic = AromaticValenceForm::aromatic(NumForm::Lit(electrons));
                let assertion = match attributes.constraints.aromatic_valence() {
                    AromaticValenceForm::Undetermined => aromatic,
                    existing => existing
                        .meet(&aromatic)
                        .ok_or(AromaticityProjectError::AtomAssertion { atom })?,
                };
                let mut update = AtomUpdate::default();
                update
                    .constraints
                    .set(AtomConstraintForm::aromatic_valence(assertion));
                edits.update_atom(AtomHandle::Id(atom), attributes, &update);
            }
            for bond in aromatic_system.bond_ids() {
                let attributes = molecule.bond(bond).attributes();
                let assertion = attributes
                    .constraints
                    .aromatic()
                    .meet(&BooleanForm::Lit(true))
                    .ok_or(AromaticityProjectError::BondAssertion { bond })?;
                let mut update = BondUpdate::default();
                update
                    .constraints
                    .set(BondConstraintForm::Aromatic(assertion));
                edits.update_bond(BondHandle::Id(bond), attributes, &update);
            }
        }
        edits.remove_aromatic_systems(
            molecule
                .aromatic_systems()
                .iter()
                .map(|aromatic_system| {
                    (
                        AromaticSystemHandle::Id(aromatic_system.id()),
                        aromatic_system.atom_ids().map(AtomHandle::Id).collect(),
                        aromatic_system.attributes().clone(),
                    )
                })
                .collect(),
        );
        Ok(edits)
    }

    /// Selects completions and aromatic systems for each candidate-ring component.
    ///
    /// Each combination chooses one completion per flexible atom and the systems
    /// perceived from those contributions. Atoms outside the candidate sets use
    /// their complete aromatic-valence evidence. Accepted systems remain
    /// disjoint: selection takes a whole combination's systems together.
    ///
    /// A combination must reproduce stored systems touching the component. Under
    /// the Error failure policy, it must also cover every atom whose completions
    /// all require aromaticity. MinElectronCount first prefers more covered atoms,
    /// then fewer electrons, then lexicographically ordered member sets. Surviving
    /// combinations with different completions for the same atoms compare by
    /// tie_break; equal completions use the lexicographically smallest partition.
    /// Ties and combinations restricting different atoms leave those atoms plural.
    ///
    /// More than 4096 search perception calls in a component, a
    /// nonliteral stored contribution outside the candidate sets, or undetermined
    /// perception yields Underdetermined. Contradictory evidence outside the candidate
    /// sets yields Contradictory. An uncovered atom requiring aromaticity
    /// yields Contradictory under the Error failure policy. If no positive aromatic
    /// assertion or stored aromatic system requires selection, returns state unchanged.
    /// Leaves state.tie_breaks unchanged; use [`Self::select_with_report`] to include
    /// this selection's tie-break uses.
    ///
    /// # Assumes
    ///
    /// State is contextually valid for the molecule: its candidate sets were admitted
    /// for that molecule, and its pending aromatic systems are compatible with it.
    ///
    /// # Errors
    ///
    /// Returns errors from aromaticity perception, including missing model parameters.
    ///
    /// # Semantic properties
    ///
    /// Within the enumeration bound, completions, accepted systems, and outcome equal
    /// exhaustive selection: pruning removes no valid combination. Cross-checked
    /// through Resolver::resolve_with_report under Hückel, HMO, and Clar.
    pub(crate) fn select(
        &self,
        molecule: &Molecule,
        state: ResolveState,
        tie_break: ValenceTieBreak,
    ) -> Result<Solution<ResolveState, AromaticityContradiction>, AromaticityError> {
        self.select_completions(molecule, state, tie_break, false)
    }

    /// Selects completions and aromatic systems, including tie-break uses in the state.
    ///
    /// Uses the same selection and failure policies as [`Self::select`]. Adds atoms
    /// selected by the valence key and members of systems selected by the structural
    /// tie-break to state.tie_breaks, sorting and deduplicating the result. If aromatic
    /// selection is unnecessary or underdetermined, keeps the existing tie-break list.
    ///
    /// # Assumes
    ///
    /// The same contextual validity requirements as [`Self::select`].
    ///
    /// # Errors
    ///
    /// Returns the same aromaticity perception errors as select.
    ///
    /// # Semantic properties
    ///
    /// Completions, accepted systems, contradictions, and errors agree with select.
    /// Within the enumeration bound, the full solution equals exhaustive selection,
    /// including tie-break uses, as checked by the property test target.
    pub(crate) fn select_with_report(
        &self,
        molecule: &Molecule,
        state: ResolveState,
        tie_break: ValenceTieBreak,
    ) -> Result<Solution<ResolveState, AromaticityContradiction>, AromaticityError> {
        self.select_completions(molecule, state, tie_break, true)
    }

    pub(crate) fn plan_bond_discharge(
        &self,
        molecule: &Molecule,
        aromatic_systems: &[(Vec<AtomId>, AromaticSystemForm)],
    ) -> Edits {
        let mut bonds = molecule
            .bonds()
            .iter()
            .filter_map(|bond| {
                let Some(BondConstraintForm::Aromatic(asserted)) = bond
                    .attributes()
                    .constraints
                    .get(BondConstraintKey::Aromatic)
                else {
                    return None;
                };
                (!bond.is_in_aromatic_system()).then_some((bond, asserted))
            })
            .peekable();
        let mut edits = Edits::new();
        if bonds.peek().is_none() {
            return edits;
        }
        let mut membership = vec![None; molecule.atoms().count()];
        for (aromatic_system, (atoms, _)) in aromatic_systems.iter().enumerate() {
            for atom in atoms {
                membership[atom.index()] = Some(aromatic_system);
            }
        }
        for (bond, asserted) in bonds {
            let [a, b] = bond.atom_ids();
            let aromatic =
                membership[a.index()].is_some() && membership[a.index()] == membership[b.index()];
            if asserted.meet(&BooleanForm::Lit(aromatic)).is_some() {
                edits.push(Edit::ModifyBondConstraint {
                    id: BondHandle::Id(bond.id()),
                    old: Some(BondConstraintForm::Aromatic(*asserted)),
                    new: None,
                });
            }
        }
        edits
    }

    fn select_completions(
        &self,
        molecule: &Molecule,
        state: ResolveState,
        tie_break: ValenceTieBreak,
        report: bool,
    ) -> Result<Solution<ResolveState, AromaticityContradiction>, AromaticityError> {
        let requires_aromaticity = molecule.aromatic_systems().count() != 0
            || molecule.atoms().iter().any(|atom| {
                matches!(
                    atom.attributes().constraints.aromatic_valence(),
                    AromaticValenceForm::Aromatic(_)
                )
            })
            || molecule.bonds().iter().any(|bond| {
                matches!(
                    bond.attributes().constraints.aromatic(),
                    BooleanForm::Lit(true)
                )
            })
            || molecule.dative_bonds().iter().any(|bond| {
                matches!(
                    bond.attributes().constraints.aromatic(),
                    BooleanForm::Lit(true)
                )
            });
        if !requires_aromaticity {
            return Ok(Solution::Determined(state));
        }

        let ResolveState {
            mut completions,
            mut aromatic_systems,
            mut tie_breaks,
        } = state;
        let mut stored = vec![AromaticValence::NotAromatic; molecule.atoms().count()];
        let mut underdetermined = false;
        for atom in molecule.atoms().ids() {
            if completions.get(atom).is_some() {
                continue;
            }
            match read_aromatic_valence(molecule, atom) {
                Ok(Some(value)) => stored[atom.index()] = value,
                Ok(None) => underdetermined = true,
                Err(AromaticityInconsistency::AromaticSystemFailure { .. })
                    if self.config.aromatic_system_failure == AromaticityFailurePolicy::Keep => {}
                Err(AromaticityInconsistency::AromaticValenceFailure { .. })
                    if self.config.aromatic_valence_failure == AromaticityFailurePolicy::Keep => {}
                Err(inconsistency) => {
                    return Ok(Solution::Contradictory(inconsistency.into()));
                }
            }
        }
        if underdetermined {
            return Ok(Solution::Underdetermined(ResolveState {
                completions,
                aromatic_systems,
                tie_breaks,
            }));
        }

        if let Some(inconsistency) = self
            .perception
            .check_aromatic_bonds(molecule)
            .into_iter()
            .next()
        {
            return Ok(Solution::Contradictory(inconsistency.into()));
        }
        let contribution = |form: &AtomForm| -> AromaticValence {
            form.constraints
                .aromatic_valence()
                .as_lit()
                .unwrap_or(AromaticValence::NotAromatic)
        };

        // Flexible: carrier atoms whose disjuncts differ in contribution.
        let flexible: Vec<(AtomId, Vec<AromaticValence>)> = completions
            .iter()
            .filter_map(|(atom, disjuncts)| {
                let contributions: Vec<AromaticValence> =
                    disjuncts.iter().map(contribution).collect();
                (contributions.iter().collect::<BTreeSet<_>>().len() > 1)
                    .then_some((atom, contributions))
            })
            .collect();

        let rings = self
            .perception
            .candidate_rings(molecule, self.config.perception);
        let required_atoms: Vec<_> = molecule
            .atoms()
            .iter()
            .filter_map(|atom| {
                completions
                    .get(atom.id())
                    .map_or_else(
                        || {
                            atom.aromatic_system_id().is_none()
                                && matches!(
                                    atom.attributes().constraints.aromatic_valence(),
                                    AromaticValenceForm::Aromatic(_)
                                )
                        },
                        |disjuncts| {
                            disjuncts.iter().all(|form| {
                                matches!(
                                    form.constraints.aromatic_valence(),
                                    AromaticValenceForm::Aromatic(_)
                                )
                            })
                        },
                    )
                    .then_some(atom.id())
            })
            .collect();
        if flexible.is_empty() {
            let mut derivation = match self.perception.derive_from_rings(
                molecule,
                &rings,
                self.config.perception,
                &required_atoms,
                |atom| {
                    completions
                        .get(atom)
                        .map_or(stored[atom.index()], |forms| contribution(&forms[0]))
                },
            )? {
                Solution::Determined(derivation) => derivation,
                Solution::Underdetermined(_) => {
                    return Ok(Solution::Underdetermined(ResolveState {
                        completions,
                        aromatic_systems,
                        tie_breaks,
                    }))
                }
                Solution::Contradictory(contradiction) => {
                    return Ok(Solution::Contradictory(contradiction))
                }
            };
            let mut failed_atoms = BTreeSet::new();
            for inconsistency in derivation.inconsistencies {
                match inconsistency {
                    AromaticityInconsistency::AromaticValenceFailure { .. }
                    | AromaticityInconsistency::AromaticBondConstraintFailure { .. }
                        if self.config.aromatic_valence_failure
                            == AromaticityFailurePolicy::Keep => {}
                    AromaticityInconsistency::AromaticSystemFailure { aromatic_system }
                        if self.config.aromatic_system_failure
                            == AromaticityFailurePolicy::Keep =>
                    {
                        failed_atoms.extend(molecule.aromatic_system(aromatic_system).atom_ids());
                    }
                    _ => return Ok(Solution::Contradictory(inconsistency.into())),
                }
            }
            if !failed_atoms.is_empty() {
                let components = self
                    .perception
                    .candidate_components(molecule, &rings, |atom| {
                        matches!(
                            completions
                                .get(atom)
                                .map_or(stored[atom.index()], |forms| contribution(&forms[0])),
                            AromaticValence::Aromatic(_)
                        )
                    });
                for component in components {
                    if !component.is_disjoint(&failed_atoms) {
                        failed_atoms.extend(component);
                    }
                }
                derivation
                    .aromatic_systems
                    .retain(|(atoms, _)| atoms.iter().all(|atom| !failed_atoms.contains(atom)));
            }
            aromatic_systems.extend(derivation.aromatic_systems);
            if report {
                tie_breaks.sort_unstable();
                tie_breaks.dedup();
            }
            return Ok(Solution::Determined(ResolveState {
                completions,
                aromatic_systems,
                tie_breaks,
            }));
        }
        let components = self
            .perception
            .candidate_components(molecule, &rings, |atom_id| match completions.get(atom_id) {
                Some(disjuncts) => disjuncts
                    .iter()
                    .any(|form| matches!(contribution(form), AromaticValence::Aromatic(_))),
                None => matches!(stored[atom_id.index()], AromaticValence::Aromatic(_)),
            });

        // A completion combination per index choice: the partition (the perceived
        // systems inside the component, sorted by member list) and the
        // restriction (the chosen forms of flexible member atoms, ascending).
        // Validity and selection act on whole completion combinations, never mixing
        // systems across them, so accepted systems are disjoint.
        type CompletionCombination = (
            Vec<(AtomId, AtomForm)>,
            Vec<(Vec<AtomId>, AromaticSystemForm)>,
        );
        let required_bonds: SmallVec<[(BondId, [AtomId; 2]); 16]> = molecule
            .bonds()
            .iter()
            .filter(|bond| {
                self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error
                    && matches!(
                        bond.attributes().constraints.aromatic(),
                        BooleanForm::Lit(true)
                    )
            })
            .map(|bond| (bond.id(), bond.atom_ids()))
            .collect();
        let mut accepted: Vec<(Vec<AtomId>, AromaticSystemForm)> = Vec::new();
        let mut tie_break_uses: BTreeSet<AtomId> = BTreeSet::new();
        let mut claimed: BTreeSet<AtomId> = BTreeSet::new();
        let mut claimed_bonds: SmallVec<[bool; 16]> = required_bonds
            .iter()
            .map(|(bond, _)| molecule.bond(*bond).is_in_aromatic_system())
            .collect();
        for component in &components {
            // Visit ring atoms first, followed by any other pi-component atoms.
            let mut ordered_atoms: Vec<AtomId> = Vec::new();
            for ring in rings.iter() {
                if !ring.atoms().iter().all(|atom| component.contains(atom)) {
                    continue;
                }
                let mut atoms = ring.atoms().to_vec();
                atoms.sort_unstable();
                for atom in atoms {
                    if flexible
                        .iter()
                        .any(|&(flexible_atom, _)| flexible_atom == atom)
                        && !ordered_atoms.contains(&atom)
                    {
                        ordered_atoms.push(atom);
                    }
                }
            }
            for &(atom, _) in &flexible {
                if component.contains(&atom) && !ordered_atoms.contains(&atom) {
                    ordered_atoms.push(atom);
                }
            }
            let component_flexible: Vec<(AtomId, Vec<AromaticValence>)> = ordered_atoms
                .iter()
                .map(|&atom| {
                    flexible
                        .iter()
                        .find(|&&(flexible_atom, _)| flexible_atom == atom)
                        .cloned()
                        .expect("ordered atom is flexible")
                })
                .collect();
            let flexible_positions: BTreeMap<AtomId, usize> = component_flexible
                .iter()
                .enumerate()
                .map(|(position, &(atom, _))| (atom, position))
                .collect();
            // Atoms whose every disjunct requires aromaticity: the carrier
            // elimination criterion, shared by the search prunes and the
            // validity filter below.
            let aromatic_only: Vec<AtomId> = required_atoms
                .iter()
                .copied()
                .filter(|atom| component.contains(atom))
                .collect();
            let fixed_contributions: BTreeMap<AtomId, AromaticValence> = component
                .iter()
                .filter(|atom| !flexible_positions.contains_key(atom))
                .map(|&atom| {
                    (
                        atom,
                        match completions.get(atom) {
                            Some(disjuncts) => contribution(&disjuncts[0]),
                            None => stored[atom.index()],
                        },
                    )
                })
                .collect();
            let remaining_contributions = |atom: AtomId, path: &[usize]| -> &[AromaticValence] {
                match flexible_positions.get(&atom) {
                    Some(&position) => {
                        let options = &component_flexible[position].1;
                        if position < path.len() {
                            slice::from_ref(&options[path[position]])
                        } else {
                            options
                        }
                    }
                    None => slice::from_ref(&fixed_contributions[&atom]),
                }
            };
            let can_contribute = |atom: AtomId, path: &[usize]| {
                if let Some(&position) = flexible_positions.get(&atom) {
                    if position < path.len() {
                        return matches!(
                            component_flexible[position].1[path[position]],
                            AromaticValence::Aromatic(_)
                        );
                    }
                }
                match completions.get(atom) {
                    Some(disjuncts) => disjuncts
                        .iter()
                        .any(|form| matches!(contribution(form), AromaticValence::Aromatic(_))),
                    None => matches!(stored[atom.index()], AromaticValence::Aromatic(_)),
                }
            };
            let require_aromaticity =
                self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error;
            let prune_enabled = require_aromaticity
                && self.config.atom_selection_algorithm == AtomSelectionAlgorithm::BacktrackingDfs;
            let candidate_aromatic_systems = |path: &[usize]| {
                let candidates: BTreeSet<_> = aromatic_only
                    .iter()
                    .flat_map(|&atom| {
                        self.perception
                            .candidate_aromatic_systems(molecule, &rings, atom, |atom| {
                                can_contribute(atom, path)
                            })
                    })
                    .collect();
                let mut candidates: Vec<_> = candidates.into_iter().collect();
                candidates.sort_by_key(BTreeSet::len);
                candidates
            };
            let candidates = if prune_enabled {
                candidate_aromatic_systems(&[])
            } else {
                Vec::new()
            };
            let certainly_invalid = |path: &[usize]| -> bool {
                let eligibility_changed = path.iter().enumerate().any(|(position, &choice)| {
                    component_flexible[position].1[choice] == AromaticValence::NotAromatic
                });
                let updated;
                let candidates = if eligibility_changed {
                    updated = candidate_aromatic_systems(path);
                    &updated
                } else {
                    &candidates
                };
                let mut unclaimed: SmallVec<[AtomId; 32]> = aromatic_only.iter().copied().collect();
                for members in candidates {
                    if !unclaimed.iter().any(|atom| members.contains(atom)) {
                        continue;
                    }
                    let contributions: SmallVec<[&[AromaticValence]; 16]> = members
                        .iter()
                        .map(|&atom| remaining_contributions(atom, path))
                        .collect();
                    if self.perception.accepts_contributions(&contributions) {
                        unclaimed.retain(|atom| !members.contains(atom));
                        if unclaimed.is_empty() {
                            return false;
                        }
                    }
                }
                !unclaimed.is_empty()
            };

            // Depth-first search over the flexible atoms; a subtree is cut
            // only when every completion below it is certainly invalid, so
            // the valid combinations equal the flat enumeration's.
            let mut combinations: BTreeMap<CompletionCombination, usize> = BTreeMap::new();
            let mut aromatic_system_failure = None;
            let mut path: Vec<usize> = Vec::new();
            let mut next = 0usize;
            let mut perception_calls = 0;
            loop {
                if path.len() == component_flexible.len() && next == 0 {
                    if perception_calls == MAX_PERCEPTION_CALLS {
                        return Ok(Solution::Underdetermined(ResolveState {
                            completions,
                            aromatic_systems,
                            tie_breaks,
                        }));
                    }
                    perception_calls += 1;
                    let choice: BTreeMap<AtomId, usize> = component_flexible
                        .iter()
                        .zip(&path)
                        .map(|(&(atom, _), &index)| (atom, index))
                        .collect();
                    let outcome = self.perception.derive_from_rings(
                        molecule,
                        &rings,
                        self.config.perception,
                        &required_atoms,
                        |atom| match choice.get(&atom) {
                            Some(&index) => {
                                contribution(&completions.get(atom).expect("flexible atom")[index])
                            }
                            None => match completions.get(atom) {
                                Some(disjuncts) => contribution(&disjuncts[0]),
                                None => stored[atom.index()],
                            },
                        },
                    )?;
                    let found = match outcome {
                        Solution::Determined(found) => found,
                        Solution::Underdetermined(_) => {
                            return Ok(Solution::Underdetermined(ResolveState {
                                completions,
                                aromatic_systems,
                                tie_breaks,
                            }));
                        }
                        Solution::Contradictory(contradiction) => {
                            return Ok(Solution::Contradictory(contradiction));
                        }
                    };
                    // Systems outside the component arise under this
                    // component's default indices; their own enumeration
                    // accumulates them.
                    let mut partition: Vec<(Vec<AtomId>, AromaticSystemForm)> = found
                        .aromatic_systems
                        .into_iter()
                        .filter(|(atoms, _)| atoms.iter().all(|atom| component.contains(atom)))
                        .collect();
                    partition.sort_by(|(a, _), (b, _)| a.cmp(b));
                    let restriction: Vec<(AtomId, AtomForm)> = choice
                        .iter()
                        .filter(|(atom, _)| partition.iter().any(|(atoms, _)| atoms.contains(atom)))
                        .map(|(&atom, &index)| {
                            (
                                atom,
                                completions.get(atom).expect("flexible atom")[index].clone(),
                            )
                        })
                        .collect();
                    let mut atom_failure = false;
                    let mut bond_failure = false;
                    let mut candidate_failure = None;
                    for inconsistency in found.inconsistencies {
                        match inconsistency {
                            AromaticityInconsistency::AromaticValenceFailure { atom } => {
                                atom_failure |= require_aromaticity && component.contains(&atom);
                            }
                            AromaticityInconsistency::AromaticBondConstraintFailure { bond } => {
                                bond_failure |= require_aromaticity
                                    && molecule
                                        .bond(bond)
                                        .atom_ids()
                                        .iter()
                                        .all(|atom| component.contains(atom));
                            }
                            AromaticityInconsistency::AromaticSystemFailure { aromatic_system } => {
                                if molecule
                                    .aromatic_system(aromatic_system)
                                    .atom_ids()
                                    .any(|atom| component.contains(&atom))
                                {
                                    candidate_failure.get_or_insert(inconsistency);
                                }
                            }
                            _ => return Ok(Solution::Contradictory(inconsistency.into())),
                        }
                    }
                    if !atom_failure {
                        if let Some(failure) = candidate_failure {
                            aromatic_system_failure = Some(
                                aromatic_system_failure
                                    .map_or(failure, |previous: AromaticityInconsistency| {
                                        previous.min(failure)
                                    }),
                            );
                        }
                    }
                    if !atom_failure && candidate_failure.is_none() {
                        for (atoms, _) in &partition {
                            claimed.extend(atoms.iter().copied());
                        }
                    }
                    let rejected = atom_failure || bond_failure || candidate_failure.is_some();
                    if !rejected {
                        let discovery_order = combinations.len();
                        combinations
                            .entry((restriction, partition))
                            .or_insert(discovery_order);
                    }
                    match path.pop() {
                        Some(index) => next = index + 1,
                        None => break,
                    }
                    continue;
                }
                if next >= component_flexible[path.len()].1.len() {
                    match path.pop() {
                        Some(index) => next = index + 1,
                        None => break,
                    }
                    continue;
                }
                path.push(next);
                if prune_enabled && certainly_invalid(&path) {
                    let index = path.pop().expect("pushed above");
                    next = index + 1;
                } else {
                    next = 0;
                }
            }

            // Validity under the `Error` failure policy first — the search
            // prunes by the same criterion: a combination may not leave an
            // aromatic-only component atom unclaimed.
            let mut ordered_combinations: Vec<_> = combinations.into_iter().collect();
            ordered_combinations.sort_unstable_by_key(|(_, discovery_order)| *discovery_order);
            let mut valid: Vec<_> = ordered_combinations
                .into_iter()
                .map(|(combination, _)| combination)
                .collect();
            if valid.is_empty() {
                if let Some(inconsistency) = aromatic_system_failure {
                    if self.config.aromatic_system_failure == AromaticityFailurePolicy::Error {
                        return Ok(Solution::Contradictory(inconsistency.into()));
                    }
                    claimed.extend(component.iter().copied());
                    for (index, &(_, endpoints)) in required_bonds.iter().enumerate() {
                        if endpoints.iter().any(|atom| component.contains(atom)) {
                            claimed_bonds[index] = true;
                        }
                    }
                }
                continue;
            }
            for (_, partition) in &valid {
                for (index, &(_, endpoints)) in required_bonds.iter().enumerate() {
                    if partition
                        .iter()
                        .any(|(atoms, _)| endpoints.iter().all(|atom| atoms.contains(atom)))
                    {
                        claimed_bonds[index] = true;
                    }
                }
            }

            // Structural order under `MinElectronCount`: most claimed atoms
            // first, then the smallest electron total, member-set lists
            // breaking ties lexicographically. The electron component acts
            // whenever the rule admits more than one total for the same
            // members; the coverage component only under a non-`Error`
            // failure policy.
            let mut structural_decided = false;
            if self.tie_break == AromaticityTieBreak::MinElectronCount {
                let structure = |partition: &Vec<(Vec<AtomId>, AromaticSystemForm)>| {
                    (
                        Reverse(
                            partition
                                .iter()
                                .map(|(atoms, _)| atoms.len())
                                .sum::<usize>(),
                        ),
                        partition
                            .iter()
                            .map(|(_, form)| match &form.electrons {
                                ElectronCountsForm::Lit(electrons) => electrons.iter().sum(),
                                _ => i64::MAX,
                            })
                            .sum::<i64>(),
                        partition
                            .iter()
                            .map(|(atoms, _)| atoms.clone())
                            .collect::<Vec<_>>(),
                    )
                };
                let best = valid
                    .iter()
                    .map(|(_, partition)| structure(partition))
                    .min()
                    .expect("non-empty survivors");
                structural_decided = report
                    && valid
                        .iter()
                        .any(|(_, partition)| structure(partition) != best);
                valid.retain(|(_, partition)| structure(partition) == best);
            }

            // Selection: a unique restriction is accepted; several compare by
            // the value key when they restrict the same atoms; otherwise the
            // members stay plural. Restriction-identical survivors take the
            // lexicographically smallest partition — representation
            // canonicalization, not policy.
            let mut restrictions: Vec<&Vec<(AtomId, AtomForm)>> =
                valid.iter().map(|(restriction, _)| restriction).collect();
            if restrictions.len() > 1 {
                let mut seen = HashSet::new();
                restrictions.retain(|restriction| seen.insert(*restriction));
            }
            let (winner_restriction, by_key) = if restrictions.len() == 1 {
                (restrictions[0], false)
            } else {
                if tie_break.key().is_empty() {
                    continue;
                }
                let domain: BTreeSet<AtomId> =
                    restrictions[0].iter().map(|(atom, _)| *atom).collect();
                if !restrictions.iter().all(|restriction| {
                    restriction
                        .iter()
                        .map(|(atom, _)| *atom)
                        .collect::<BTreeSet<_>>()
                        == domain
                }) {
                    continue;
                }
                let best = restrictions
                    .iter()
                    .copied()
                    .max_by(|a, b| compare_restrictions(a, b, tie_break))
                    .expect("non-empty restrictions");
                let unique = restrictions
                    .iter()
                    .filter(|restriction| {
                        compare_restrictions(restriction, best, tie_break).is_eq()
                    })
                    .count()
                    == 1;
                if !unique {
                    continue;
                }
                (best, true)
            };
            let (_, partition) = valid
                .iter()
                .filter(|(restriction, _)| restriction == winner_restriction)
                .min_by(|(_, a), (_, b)| {
                    a.iter()
                        .map(|(atoms, _)| atoms)
                        .cmp(b.iter().map(|(atoms, _)| atoms))
                })
                .expect("winner restriction present");
            accepted.extend(partition.iter().cloned());
            if structural_decided {
                for (atoms, _) in partition {
                    tie_break_uses.extend(atoms.iter().copied());
                }
            }
            for (atom, form) in winner_restriction {
                if report && by_key && completions.get(*atom).is_some_and(|entry| entry.len() > 1) {
                    tie_break_uses.insert(*atom);
                }
                completions.insert(*atom, smallvec![form.clone()]);
            }
        }

        let missing_atoms: Vec<_> = required_atoms
            .iter()
            .copied()
            .filter(|atom| !claimed.contains(atom))
            .collect();
        let outside_search = |aromatic_system| {
            molecule
                .aromatic_system(aromatic_system)
                .atom_ids()
                .all(|atom| {
                    components
                        .iter()
                        .all(|component| !component.contains(&atom))
                })
        };
        if (self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error
            && (!missing_atoms.is_empty() || claimed_bonds.iter().any(|claimed| !claimed)))
            || molecule.aromatic_systems().ids().any(outside_search)
        {
            let inconsistencies =
                self.perception
                    .check_aromatic_systems(molecule, &missing_atoms, &accepted);
            for inconsistency in inconsistencies {
                let reject = match inconsistency {
                    AromaticityInconsistency::AromaticValenceFailure { .. } => {
                        self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error
                    }
                    AromaticityInconsistency::AromaticBondConstraintFailure { bond } => {
                        required_bonds
                            .iter()
                            .zip(&claimed_bonds)
                            .any(|((required, _), claimed)| *required == bond && !claimed)
                    }
                    AromaticityInconsistency::AromaticSystemFailure { aromatic_system } => {
                        self.config.aromatic_system_failure == AromaticityFailurePolicy::Error
                            && outside_search(aromatic_system)
                    }
                    _ => true,
                };
                if reject {
                    return Ok(Solution::Contradictory(inconsistency.into()));
                }
            }
        }

        aromatic_systems.extend(accepted);
        if report {
            tie_breaks.extend(tie_break_uses);
            tie_breaks.sort_unstable();
            tie_breaks.dedup();
        }
        Ok(Solution::Determined(ResolveState {
            completions,
            aromatic_systems,
            tie_breaks,
        }))
    }
}

fn compare_restrictions(
    a: &[(AtomId, AtomForm)],
    b: &[(AtomId, AtomForm)],
    tie_break: ValenceTieBreak,
) -> Ordering {
    let key = tie_break.key();
    let b_forms: BTreeMap<AtomId, &AtomForm> = b.iter().map(|(atom, form)| (*atom, form)).collect();
    let mut a_sorted: Vec<&(AtomId, AtomForm)> = a.iter().collect();
    a_sorted.sort_unstable_by_key(|(atom, _)| *atom);
    for (atom, a_form) in a_sorted {
        let Some(b_form) = b_forms.get(atom) else {
            continue;
        };
        let ordering = compare_by_key(key, a_form, b_form);
        if ordering != Ordering::Equal {
            return ordering;
        }
    }
    Ordering::Equal
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_graph_core::{
        ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
        RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
    };
    use umol_graph_ir::ir::{
        AddAromaticSystem, AromaticSystemId, BondId, Edit, Edits, NumForm, RingConfig,
        UnpairedElectronsForm,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_concrete};

    use super::*;
    use crate::ops::model::{
        AromaticityRule, AromaticityTieBreak, ChemistryModel, ElementScope, RingLimits,
    };
    use crate::ops::resolve::{ResolveConfig, ResolveContradiction, Resolver};
    use crate::ops::valence::AtomCompletions;

    #[fixture]
    fn aromaticity_model() -> AromaticityModel {
        AromaticityModel {
            scope: ElementScope::Any,
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
            tie_break: AromaticityTieBreak::Strict,
        }
    }

    #[fixture]
    fn benzene() -> Molecule {
        mol_dsl!(
            r#"{
            :atoms ["C#i=#c0#h#n0#u0#s#v2#a" "C#i=#c0#h#n0#u0#s#v2#a"
                    "C#i=#c0#h#n0#u0#s#v2#a" "C#i=#c0#h#n0#u0#s#v2#a"
                    "C#i=#c0#h#n0#u0#s#v2#a" "C#i=#c0#h#n0#u0#s#v2#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#
        )
    }

    #[fixture]
    fn quinoline() -> Molecule {
        mol_dsl!(
            r#"{
            :atoms ["C#c0#a" "C#c0#a" "C#c0#a" "C#c0#a" "C#c0#a"
                    "C#c0#a" "C#c0#a" "C#c0#a" "N#c0" "C#c0#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]]}"#
        )
    }

    #[fixture]
    fn resolved_benzene() -> Molecule {
        mol_dsl_concrete!(
            r#"{
            :atoms ["C#h1" "C#h1" "C#h1" "C#h1" "C#h1" "C#h1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        )
    }

    #[fixture]
    fn aromatic_valence_mismatch() -> Molecule {
        mol_dsl!(
            r#"{
            :atoms ["C#a2" "C#a0" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        )
    }

    #[fixture]
    fn aromatic_bond_constraint_mismatch() -> Molecule {
        mol_dsl!(
            r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a!"] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        )
    }

    #[fixture]
    fn overlapping_aromatic_system() -> Molecule {
        mol_dsl!(
            r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        )
    }

    #[rstest]
    fn test_aromaticity_resolve_config_default() {
        assert_eq!(
            AromaticityResolveConfig::default(),
            AromaticityResolveConfig {
                perception: AromaticityConfig::default(),
                atom_selection_algorithm: AtomSelectionAlgorithm::BacktrackingDfs,
                aromatic_valence_failure: AromaticityFailurePolicy::Error,
                aromatic_system_failure: AromaticityFailurePolicy::Error,
            }
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_plan_resolve(
        aromaticity_model: AromaticityModel,
        benzene: Molecule,
    ) {
        assert_eq!(
            AromaticityResolver::with_config(
                &aromaticity_model,
                AromaticityResolveConfig {
                    perception: AromaticityConfig {
                        ring_config: RingConfig {
                            simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                        },
                        connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                        maximum_independent_set_algorithm:
                            MaximumIndependentSetAlgorithm::BranchAndBound,
                    },
                    ..AromaticityResolveConfig::default()
                },
            )
            .plan_resolve(&benzene),
            Ok(Solution::Determined(Edits::from_iter([
                Edit::AddAromaticSystem(Box::new(AddAromaticSystem {
                    atoms: (0..6).map(|id| AtomHandle::Id(AtomId(id))).collect(),
                    attributes: AromaticSystemForm::from_electrons(vec![1; 6])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                })),
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(0)),
                    old: Some(AtomConstraintForm::AromaticValence(
                        AromaticValenceForm::Aromatic(NumForm::Lit(1)),
                    )),
                    new: None,
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(1)),
                    old: Some(AtomConstraintForm::AromaticValence(
                        AromaticValenceForm::Aromatic(NumForm::Lit(1)),
                    )),
                    new: None,
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(2)),
                    old: Some(AtomConstraintForm::AromaticValence(
                        AromaticValenceForm::Aromatic(NumForm::Lit(1)),
                    )),
                    new: None,
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(3)),
                    old: Some(AtomConstraintForm::AromaticValence(
                        AromaticValenceForm::Aromatic(NumForm::Lit(1)),
                    )),
                    new: None,
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(4)),
                    old: Some(AtomConstraintForm::AromaticValence(
                        AromaticValenceForm::Aromatic(NumForm::Lit(1)),
                    )),
                    new: None,
                },
                Edit::ModifyAtomConstraint {
                    id: AtomHandle::Id(AtomId(5)),
                    old: Some(AtomConstraintForm::AromaticValence(
                        AromaticValenceForm::Aromatic(NumForm::Lit(1)),
                    )),
                    new: None,
                },
            ])))
        );
    }

    #[rstest]
    #[case::open_contribution(
        mol_dsl!(r#"{
            :atoms ["C#a+" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
        }"#)
    )]
    fn test_aromaticity_resolver_plan_resolve_partial(
        aromaticity_model: AromaticityModel,
        #[case] mut molecule: Molecule,
    ) {
        let original = molecule.clone();
        let resolver = AromaticityResolver::new(&aromaticity_model);
        assert_eq!(
            resolver.plan_resolve(&molecule),
            Ok(Solution::Underdetermined(Edits::new()))
        );
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_aromaticity_resolver_plan_resolve_aromatic_valence_mismatch(
        aromaticity_model: AromaticityModel,
        aromatic_valence_mismatch: Molecule,
        #[values(AromaticityFailurePolicy::Error, AromaticityFailurePolicy::Keep)]
        failure: AromaticityFailurePolicy,
    ) {
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_valence_failure: failure,
                aromatic_system_failure: failure,
                ..AromaticityResolveConfig::default()
            },
        );
        assert_eq!(
            resolver.plan_resolve(&aromatic_valence_mismatch),
            Ok(Solution::Contradictory(
                AromaticityInconsistency::AromaticValenceMismatch {
                    atom: AtomId(0),
                    aromatic_system: AromaticSystemId(0),
                }
                .into()
            )),
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_plan_resolve_aromatic_bond_constraint_mismatch(
        aromaticity_model: AromaticityModel,
        aromatic_bond_constraint_mismatch: Molecule,
        #[values(AromaticityFailurePolicy::Error, AromaticityFailurePolicy::Keep)]
        failure: AromaticityFailurePolicy,
    ) {
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_valence_failure: failure,
                aromatic_system_failure: failure,
                ..AromaticityResolveConfig::default()
            },
        );
        assert_eq!(
            resolver.plan_resolve(&aromatic_bond_constraint_mismatch),
            Ok(Solution::Contradictory(
                AromaticityInconsistency::AromaticBondConstraintMismatch {
                    bond: BondId(0),
                    aromatic_system: AromaticSystemId(0),
                }
                .into()
            )),
        );
    }

    #[rstest]
    #[case::missing_evidence(
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#h2"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"]
                    [4 5 "1"] [5 6 "1"] [6 0 "1"]]
        }"#),
        AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
    )]
    #[case::bond_evidence(
        mol_dsl!(r#"{:atoms ["C#a!" "C#a"] :bonds [[0 1 "1#a"]]}"#),
        AromaticityInconsistency::AromaticValenceBondMismatch { atom: AtomId(0), bond: BondId(0) },
    )]
    #[case::negative_aromatic_system_member(
        mol_dsl!(r#"{
            :atoms ["C#a!" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        AromaticityInconsistency::AromaticValenceMismatch { atom: AtomId(0), aromatic_system: AromaticSystemId(0) },
    )]
    fn test_aromaticity_resolver_plan_resolve_evidence(
        #[case] mut molecule: Molecule,
        #[case] expected: AromaticityInconsistency,
    ) {
        let original = molecule.clone();
        let resolver = AromaticityResolver::new(&AromaticityModel::daylight());
        assert_eq!(
            resolver.plan_resolve(&molecule),
            Ok(Solution::Contradictory(expected.into()))
        );
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Ok(Solution::Contradictory(expected.into()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::conformant(
        AromaticityModel::daylight(),
        AromaticityResolveConfig::default(),
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::rejected_projections(
        AromaticityModel::mdl(),
        AromaticityResolveConfig {
            aromatic_valence_failure: AromaticityFailurePolicy::Keep,
            ..AromaticityResolveConfig::default()
        },
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
        }"#)
    )]
    #[case::rejected_existing_aromatic_system(
        AromaticityModel::mdl(),
        AromaticityResolveConfig {
            aromatic_valence_failure: AromaticityFailurePolicy::Keep,
            aromatic_system_failure: AromaticityFailurePolicy::Keep,
            ..AromaticityResolveConfig::default()
        },
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]
        }"#)
    )]
    fn test_aromaticity_resolver_plan_resolve_identity(
        #[case] model: AromaticityModel,
        #[case] config: AromaticityResolveConfig,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            AromaticityResolver::with_config(&model, config).plan_resolve(&molecule),
            Ok(Solution::Determined(Edits::new()))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_resolve_overlap(
        overlapping_aromatic_system: Molecule,
        #[values(
            AtomSelectionAlgorithm::FullDfs,
            AtomSelectionAlgorithm::BacktrackingDfs
        )]
        atom_selection_algorithm: AtomSelectionAlgorithm,
    ) {
        let resolver = AromaticityResolver::with_config(
            &AromaticityModel::daylight(),
            AromaticityResolveConfig {
                atom_selection_algorithm,
                ..Default::default()
            },
        );
        assert_eq!(
            resolver.resolve(overlapping_aromatic_system.clone()),
            Ok(Solution::Contradictory(
                AromaticityInconsistency::AromaticSystemFailure {
                    aromatic_system: AromaticSystemId(0)
                }
                .into()
            ))
        );
        let mut output = overlapping_aromatic_system.clone();
        assert_eq!(
            resolver.resolve_into(&mut output),
            Ok(Solution::Contradictory(
                AromaticityInconsistency::AromaticSystemFailure {
                    aromatic_system: AromaticSystemId(0)
                }
                .into()
            ))
        );
        assert_eq!(output, overlapping_aromatic_system);
    }

    #[rstest]
    fn test_aromaticity_resolver_resolve_overlap_identity(
        overlapping_aromatic_system: Molecule,
        #[values(
            AtomSelectionAlgorithm::FullDfs,
            AtomSelectionAlgorithm::BacktrackingDfs
        )]
        atom_selection_algorithm: AtomSelectionAlgorithm,
    ) {
        let resolver = AromaticityResolver::with_config(
            &AromaticityModel::daylight(),
            AromaticityResolveConfig {
                aromatic_system_failure: AromaticityFailurePolicy::Keep,
                atom_selection_algorithm,
                ..Default::default()
            },
        );
        assert_eq!(
            resolver.resolve(overlapping_aromatic_system.clone()),
            Ok(Solution::Determined(overlapping_aromatic_system.clone()))
        );
        let mut output = overlapping_aromatic_system.clone();
        assert_eq!(
            resolver.resolve_into(&mut output),
            Ok(Solution::Determined(()))
        );
        assert_eq!(output, overlapping_aromatic_system);
        assert_eq!(output.clone().edit().finish(), Ok(output));
    }

    #[rstest]
    #[case::short(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1]"}]}"#))]
    #[case::long(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,0,0]"}]}"#))]
    fn test_aromaticity_resolver_resolve_electron_counts(#[case] molecule: Molecule) {
        let resolver = AromaticityResolver::new(&AromaticityModel::daylight());
        assert_eq!(
            resolver.resolve(molecule),
            Ok(Solution::Contradictory(
                AromaticityInconsistency::AromaticSystemFailure {
                    aromatic_system: AromaticSystemId(0)
                }
                .into()
            ))
        );
    }

    #[rstest]
    #[case::short(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1]"}]}"#))]
    #[case::long(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,0,0]"}]}"#))]
    fn test_aromaticity_resolver_resolve_electron_counts_identity(#[case] molecule: Molecule) {
        let resolver = AromaticityResolver::with_config(
            &AromaticityModel::daylight(),
            AromaticityResolveConfig {
                aromatic_system_failure: AromaticityFailurePolicy::Keep,
                ..Default::default()
            },
        );
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Determined(molecule))
        );
    }

    #[rstest]
    #[case::homogeneous_localized(
        AromaticityResolveConfig::default(),
        mol_dsl_concrete!(r#"{:atoms ["C#h#a" "C#h#a" "C#c+#h#a0"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#c+#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,0]"}]
        }"#)
    )]
    #[case::heterogeneous_localized(
        AromaticityResolveConfig::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["N#c+#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["N#c+#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::accepted_aromatic_system_with_rejected_projections(
        AromaticityResolveConfig {
            aromatic_valence_failure: AromaticityFailurePolicy::Keep,
            ..AromaticityResolveConfig::default()
        },
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h3#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h3#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    fn test_aromaticity_resolver_resolve(
        aromaticity_model: AromaticityModel,
        #[case] config: AromaticityResolveConfig,
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let resolver = AromaticityResolver::with_config(&aromaticity_model, config);
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    type SelectOutcome = Solution<ResolveState, AromaticityContradiction>;

    #[rstest]
    #[case::added_aromatic_system(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a#v2" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a#R+"] [1 2 "1#a"] [2 3 "1#a"]
                    [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#v2" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1#R+"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::biphenyl(
        mol_dsl_concrete!(r#"{
            :atoms ["C#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"
                    "C#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"]
                    [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]
                    [6 7 "1#a"] [7 8 "1#a"] [8 9 "1#a"]
                    [9 10 "1#a"] [10 11 "1#a"] [11 6 "1#a"] [0 6 "1#a!"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C#h" "C#h" "C#h" "C#h" "C#h"
                    "C" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [6 7 "1"] [7 8 "1"] [8 9 "1"]
                    [9 10 "1"] [10 11 "1"] [11 6 "1"] [0 6 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}
                               {:atoms [6 7 8 9 10 11] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    fn test_aromaticity_resolver_resolve_discharge(
        aromaticity_model: AromaticityModel,
        #[case] molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).resolve(molecule),
            Ok(Solution::Determined(expected))
        );
    }

    #[rstest]
    #[case::unclaimed_atom(
        mol_dsl_concrete!(r#"{:atoms ["C#h3#a"]}"#),
        mol_dsl_concrete!(r#"{:atoms ["C#h3#a"]}"#),
        AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
    )]
    #[case::unclaimed_bond(
        mol_dsl_concrete!(r#"{
            :atoms ["C#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"
                    "C#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 6 "1"]
                    [0 6 "1#a"]]}"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C#h" "C#h" "C#h" "C#h" "C#h"
                    "C" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 6 "1"]
                    [0 6 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}
                               {:atoms [6 7 8 9 10 11] :attrs "[1,1,1,1,1,1]"}]}"#),
        AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(12) },
    )]
    #[case::excluded_benzene(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
    )]
    #[case::naphthalene_peripheral_exclusion(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#a" "C#a"
                    "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]}"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C" "C"
                    "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]
            :aromatic-systems [{:atoms [4 5 6 7 8 9] :attrs "[1,1,1,1,1,1]"}]}"#),
        AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
    )]
    fn test_aromaticity_resolver_resolve_constraints(
        aromaticity_model: AromaticityModel,
        #[case] molecule: Molecule,
        #[case] kept: Molecule,
        #[case] inconsistency: AromaticityInconsistency,
        #[values(AromaticityFailurePolicy::Error, AromaticityFailurePolicy::Keep)]
        failure: AromaticityFailurePolicy,
    ) {
        let config = AromaticityResolveConfig {
            aromatic_valence_failure: failure,
            ..AromaticityResolveConfig::default()
        };
        let expected = match failure {
            AromaticityFailurePolicy::Error => Solution::Contradictory(inconsistency.into()),
            AromaticityFailurePolicy::Keep => Solution::Determined(kept.clone()),
        };
        let expected_into = match failure {
            AromaticityFailurePolicy::Error => Solution::Contradictory(inconsistency.into()),
            AromaticityFailurePolicy::Keep => Solution::Determined(()),
        };
        let expected_molecule = match failure {
            AromaticityFailurePolicy::Error => molecule.clone(),
            AromaticityFailurePolicy::Keep => kept,
        };
        let resolver = AromaticityResolver::with_config(&aromaticity_model, config);
        assert_eq!(resolver.resolve(molecule.clone()), Ok(expected.clone()));
        let mut standalone = molecule.clone();
        assert_eq!(
            resolver.resolve_into(&mut standalone),
            Ok(expected_into.clone())
        );
        assert_eq!(standalone, expected_molecule);

        let model = ChemistryModel {
            aromaticity: aromaticity_model,
            ..ChemistryModel::default()
        };
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                aromaticity: config,
                ..ResolveConfig::default()
            },
        );
        let expected = expected.map_contradiction(ResolveContradiction::Aromaticity);
        assert_eq!(resolver.resolve(molecule.clone()), Ok(expected.clone()));
        let mut composite = molecule;
        assert_eq!(
            resolver.resolve_into(&mut composite),
            Ok(expected_into.map_contradiction(ResolveContradiction::Aromaticity))
        );
        assert_eq!(composite, expected_molecule);
    }

    #[rstest]
    #[case::benzene(resolved_benzene(),
        mol_dsl_concrete!(r#"{:atoms ["C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]]}"#))]
    #[case::pyridine(mol_dsl_concrete!(r#"{
        :atoms ["N#n1" "C#h1" "C#h1" "C#h1" "C#h1" "C#h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#),
        mol_dsl_concrete!(r#"{
        :atoms ["N#n1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]]}"#))]
    #[case::pyrrole(mol_dsl_concrete!(r#"{
        :atoms ["N#h1" "C#h1" "C#h1" "C#h1" "C#h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
        mol_dsl_concrete!(r#"{
        :atoms ["N#h1#a2" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 0 "1#a"]]}"#))]
    #[case::furan(mol_dsl_concrete!(r#"{
        :atoms ["O#n1" "C#h1" "C#h1" "C#h1" "C#h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
        mol_dsl_concrete!(r#"{
        :atoms ["O#n1#a2" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 0 "1#a"]]}"#))]
    #[case::anion(mol_dsl_concrete!(r#"{
        :atoms ["C#c-#h1" "C#h1" "C#h1" "C#h1" "C#h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
        mol_dsl_concrete!(r#"{
        :atoms ["C#c-#h1#a2" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 0 "1#a"]]}"#))]
    #[case::cation(mol_dsl_concrete!(r#"{
        :atoms ["C#c+#h1" "C#h1" "C#h1" "C#h1" "C#h1" "C#h1" "C#h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[0,1,1,1,1,1,1]"}]}"#),
        mol_dsl_concrete!(r#"{
        :atoms ["C#c+#h1#a0" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 5 "1#a"] [5 6 "1#a"] [6 0 "1#a"]]}"#))]
    #[case::fused(mol_dsl_concrete!(r#"{
        :atoms ["C#h1" "C#h1" "C#h1" "C#h1" "C" "C" "C#h1" "C#h1" "C#h1" "C#h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5 6 7 8 9] :attrs "[1,1,1,1,1,1,1,1,1,1]"}]}"#),
        mol_dsl_concrete!(r#"{
        :atoms ["C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#a1" "C#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]
                [4 6 "1#a"] [6 7 "1#a"] [7 8 "1#a"] [8 9 "1#a"] [9 5 "1#a"]]}"#))]
    #[case::linked(mol_dsl_concrete!(r#"{
        :atoms ["C#h1" "C#h1" "C#h1" "C#h1" "C#h1" "C"
                "C" "C#h1" "C#h1" "C#h1" "C#h1" "C#h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"] [5 6 "1"]
                [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 6 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}
                          {:atoms [6 7 8 9 10 11] :attrs "[1,1,1,1,1,1]"}]}"#),
        mol_dsl_concrete!(r#"{
        :atoms ["C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#a1" "C#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
        :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"] [5 6 "1"]
                [6 7 "1#a"] [7 8 "1#a"] [8 9 "1#a"] [9 10 "1#a"] [10 11 "1#a"] [11 6 "1#a"]]}"#))]
    fn test_aromaticity_resolver_project(
        #[values(AromaticityModel::daylight(), AromaticityModel {
            scope: ElementScope::AllowList(vec![]), ..AromaticityModel::daylight()
        })]
        model: AromaticityModel,
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let resolver = AromaticityResolver::new(&model);
        assert_eq!(
            resolver.project(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            resolver.project_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
        assert_eq!(
            resolver.project_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::empty(Molecule::new())]
    #[case::ordinary(mol_dsl_concrete!(r#"{:atoms ["N#h3#n1" "C#h3#u1#s2" "O#n0#u2#s3"]}"#))]
    #[case::partial(mol_dsl!(r#"{:atoms ["C" "O"] :bonds [[0 1 "*"]]}"#))]
    fn test_aromaticity_resolver_project_identity(#[case] mut molecule: Molecule) {
        let model = AromaticityModel::daylight();
        let original = molecule.clone();
        assert_eq!(
            AromaticityResolver::new(&model).project(molecule.clone()),
            Ok(Solution::Determined(original.clone()))
        );
        assert_eq!(
            AromaticityResolver::new(&model).project_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::charge("[1,1,1,1,1,1]#c+#u0#s", AromaticityProjectError::ChargedAromaticSystem { aromatic_system: AromaticSystemId(0) })]
    #[case::radical("[1,1,1,1,1,1]#c0#u1#s2", AromaticityProjectError::AromaticSystemSpin { aromatic_system: AromaticSystemId(0) })]
    #[case::spin("[1,1,1,1,1,1]#c0#u0#s3", AromaticityProjectError::AromaticSystemSpin { aromatic_system: AromaticSystemId(0) })]
    #[case::charge_unknown("[1,1,1,1,1,1]#u0#s", AromaticityProjectError::NonConcreteAromaticSystem { aromatic_system: AromaticSystemId(0) })]
    #[case::electrons_unknown("*#c0#u0#s", AromaticityProjectError::NonConcreteAromaticSystem { aromatic_system: AromaticSystemId(0) })]
    #[case::unpaired_unknown("[1,1,1,1,1,1]#c0#s", AromaticityProjectError::NonConcreteAromaticSystem { aromatic_system: AromaticSystemId(0) })]
    #[case::multiplicity_unknown("[1,1,1,1,1,1]#c0#u0", AromaticityProjectError::NonConcreteAromaticSystem { aromatic_system: AromaticSystemId(0) })]
    fn test_aromaticity_resolver_project_error(
        resolved_benzene: Molecule,
        #[case] form: &str,
        #[case] expected: AromaticityProjectError,
    ) {
        let model = AromaticityModel::daylight();
        let mut editor = resolved_benzene.edit();
        *editor
            .aromatic_system_mut(AromaticSystemId(0))
            .attributes_mut() = form.parse().unwrap();
        let mut molecule = editor.finish().unwrap();
        let original = molecule.clone();
        assert_eq!(
            AromaticityResolver::new(&model).project(molecule.clone()),
            Err(expected.clone())
        );
        assert_eq!(
            AromaticityResolver::new(&model).project_into(&mut molecule),
            Err(expected)
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::atom(Some(AtomConstraintForm::aromatic_valence(AromaticValenceForm::NotAromatic)), None,
        AromaticityProjectError::AtomAssertion { atom: AtomId(5) })]
    #[case::contribution(Some(AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(2))), None,
        AromaticityProjectError::AtomAssertion { atom: AtomId(5) })]
    #[case::bond(None, Some(BondConstraintForm::Aromatic(BooleanForm::Lit(false))),
        AromaticityProjectError::BondAssertion { bond: BondId(5) })]
    fn test_aromaticity_resolver_project_constraints(
        resolved_benzene: Molecule,
        #[case] atom: Option<AtomConstraintForm>,
        #[case] bond: Option<BondConstraintForm>,
        #[case] expected: AromaticityProjectError,
    ) {
        let model = AromaticityModel::daylight();
        let mut editor = resolved_benzene.edit();
        if let Some(form) = atom {
            editor
                .atom_mut(AtomId(5))
                .attributes_mut()
                .constraints
                .set(form);
        }
        if let Some(form) = bond {
            editor
                .bond_mut(BondId(5))
                .attributes_mut()
                .constraints
                .set(form);
        }
        let mut molecule = editor.finish().unwrap();
        let original = molecule.clone();
        assert_eq!(
            AromaticityResolver::new(&model).project(molecule.clone()),
            Err(expected.clone())
        );
        assert_eq!(
            AromaticityResolver::new(&model).project_into(&mut molecule),
            Err(expected)
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::contributions(
        mol_dsl_concrete!(r#"{
            :atoms ["N#h1" "C#h1" "C#h1" "C#h1" "C#h1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#),
        mol_dsl_concrete!(r#"{
            :atoms ["N#h1#a2" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"] [4 0 "1#a"]]}"#),
    )]
    fn test_aromaticity_resolver_plan_project(
        #[case] molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let model = AromaticityModel::daylight();
        let edits = AromaticityResolver::new(&model)
            .plan_project(&molecule)
            .unwrap();
        let molecule = molecule.apply(edits).unwrap();
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_aromaticity_resolver_plan_project_error(mut resolved_benzene: Molecule) {
        resolved_benzene
            .atom_mut(AtomId(5))
            .attributes_mut()
            .constraints
            .set(AtomConstraintForm::aromatic_valence(
                AromaticValenceForm::NotAromatic,
            ));
        let model = AromaticityModel::daylight();
        assert_eq!(
            AromaticityResolver::new(&model).plan_project(&resolved_benzene),
            Err(AromaticityProjectError::AtomAssertion { atom: AtomId(5) }),
        );
    }

    #[rstest]
    #[case::aromatic_last(smallvec![
        atom_dsl!("N#c0#h2#n#a!"), atom_dsl!("N#c0#h2#n0#a2")
    ])]
    #[case::aromatic_first(smallvec![
        atom_dsl!("N#c0#h2#n0#a2"), atom_dsl!("N#c0#h2#n#a!")
    ])]
    fn test_aromaticity_resolver_select(#[case] options: SmallVec<[AtomForm; 1]>) {
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "N"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"] [0 6 "1"]]}"#
        );
        let model = AromaticityModel {
            scope: ElementScope::Any,
            rule: AromaticityRule::Hmo {
                stabilization_threshold: 0.0,
            },
            tie_break: AromaticityTieBreak::MinElectronCount,
        };
        let state = ResolveState {
            completions: AtomCompletions::from_iter([(AtomId(6), options)]),
            ..Default::default()
        };
        assert_eq!(
            AromaticityResolver::new(&model).select(&molecule, state, ValenceTieBreak::Strict),
            Ok(Solution::Determined(ResolveState {
                completions: AtomCompletions::from_iter([(
                    AtomId(6),
                    smallvec![atom_dsl!("N#c0#h2#n0#a2")]
                ),]),
                aromatic_systems: vec![(
                    (0..7).map(AtomId).collect(),
                    AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 1, 1, 2])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                )],
                tie_breaks: Vec::new(),
            }))
        );
    }

    #[rstest]
    #[case::unique_survivor(
        mol_dsl!(r#"{:atoms ["N#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                     :bonds [[0 1 "1#a+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![
                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
            ]),
            (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
        ]),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState { completions: AtomCompletions::from_iter([
                (AtomId(0), smallvec![atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2")]),
                (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            ]), aromatic_systems: vec![(
                (0..5).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )], tie_breaks: Vec::new() })
    )]
    #[case::quinoline(
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "N#c0" "C#c0"]
                     :bonds [[0 1 "1#a+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                             [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]]}"#),
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(6), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(7), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(8), smallvec![
                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
            ]),
            (AtomId(9), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
        ]),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState { completions: AtomCompletions::from_iter([
                (AtomId(0), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(6), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(7), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(8), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                (AtomId(9), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            ]), aromatic_systems: vec![(
                (0..10).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![1; 10])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )], tie_breaks: Vec::new() })
    )]
    #[case::tie_most_saturated(
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                     :bonds [[0 1 "1#a+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![
                atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
                atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
            ]),
            (AtomId(1), smallvec![
                atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
                atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
            ]),
            (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
        ]),
        ValenceTieBreak::MostSaturated,
        Solution::Determined(ResolveState { completions: AtomCompletions::from_iter([
                (AtomId(0), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0")]),
                (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2")]),
                (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            ]), aromatic_systems: vec![(
                (0..6).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![0, 2, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )], tie_breaks: vec![AtomId(0), AtomId(1)] })
    )]
    #[case::stored_assertions_empty_carrier(
        benzene(),
        AtomCompletions::new(),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState { completions: AtomCompletions::new(), aromatic_systems: vec![(
                (0..6).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![1; 6])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )], tie_breaks: Vec::new() })
    )]
    #[case::stored_aromatic_system_without_assertion(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        AtomCompletions::new(),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState {
            completions: AtomCompletions::new(),
            aromatic_systems: vec![(
                (0..6).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![1; 6])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )],
            tie_breaks: Vec::new(),
        })
    )]
    #[case::stored_aromatic_system_refines_assertion(
        mol_dsl!(r#"{
            :atoms ["C#a+" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        AtomCompletions::new(),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState {
            completions: AtomCompletions::new(),
            aromatic_systems: vec![(
                (0..6).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![1; 6])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )],
            tie_breaks: Vec::new(),
        })
    )]
    #[case::bond_evidence_mismatch(
        mol_dsl!(r#"{:atoms ["C#a!" "C#a"] :bonds [[0 1 "1#a"]]}"#),
        AtomCompletions::new(),
        ValenceTieBreak::Strict,
        Solution::Contradictory(AromaticityInconsistency::AromaticValenceBondMismatch {
            atom: AtomId(0), bond: BondId(0),
        }.into())
    )]
    #[case::aromatic_system_evidence_mismatch(
        mol_dsl!(r#"{
            :atoms ["C#a2" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        AtomCompletions::new(),
        ValenceTieBreak::Strict,
        Solution::Contradictory(AromaticityInconsistency::AromaticValenceMismatch {
            atom: AtomId(0), aromatic_system: AromaticSystemId(0),
        }.into())
    )]
    #[case::unclaimed_aromatic_contradiction(
        mol_dsl!(r#"{:atoms ["N#c0#a+"] :bonds []}"#),
        AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")])]),
        ValenceTieBreak::Strict,
        Solution::Contradictory(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }
        ))
    )]
    #[case::localized_bond_constraint(
        mol_dsl!(r#"{:atoms ["N#c0" "C#c0"] :bonds [[0 1 "1#a+"]]}"#),
        AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")])]),
        ValenceTieBreak::Strict,
        Solution::Underdetermined(ResolveState {
            completions: AtomCompletions::from_iter([(AtomId(0), smallvec![
                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
            ])]),
            aromatic_systems: vec![],
            tie_breaks: vec![],
        })
    )]
    #[case::dative_bond_constraint(
        mol_dsl!(r#"{
            :atoms ["N#c0" "C#c0"]
            :bonds []
            :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#a+"}]
        }"#),
        AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")])]),
        ValenceTieBreak::Strict,
        Solution::Contradictory(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }
        ))
    )]
    #[case::stored_undetermined_gate(
        mol_dsl!(r#"{:atoms ["C#a+" "C#a"] :bonds [[0 1 "1"]]}"#),
        AtomCompletions::new(),
        ValenceTieBreak::Strict,
        Solution::Underdetermined(ResolveState::default())
    )]
    fn test_aromaticity_resolver_select_with_report(
        aromaticity_model: AromaticityModel,
        #[case] molecule: Molecule,
        #[case] completions: AtomCompletions,
        #[case] tie_break: ValenceTieBreak,
        #[case] expected: SelectOutcome,
    ) {
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).select_with_report(
                &molecule,
                ResolveState {
                    completions,
                    ..ResolveState::default()
                },
                tie_break,
            ),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::absent_aromaticity_constraints(
        mol_dsl!(r#"{:atoms ["N#c0"] :bonds []}"#),
        AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")])]),
        ValenceTieBreak::Strict
    )]
    #[case::vacuous_aromaticity_constraints(
        mol_dsl!(r#"{
            :atoms ["N#c0#a*" "C#a*"]
            :bonds [[0 1 "1#a*"]]
            :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#a*"}]
        }"#),
        AtomCompletions::new(),
        ValenceTieBreak::Strict
    )]
    #[case::negative_aromaticity_constraints(
        mol_dsl!(r#"{
            :atoms ["N#c0#a!" "C#a!"]
            :bonds [[0 1 "1#a!"]]
            :dative-bonds [{:donors [0] :acceptor 1 :attrs "1#a!"}]
        }"#),
        AtomCompletions::new(),
        ValenceTieBreak::Strict
    )]
    #[case::tie_strict(
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                     :bonds [[0 1 "1#a+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![
                atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
                atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
            ]),
            (AtomId(1), smallvec![
                atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
                atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
            ]),
            (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
        ]),
        ValenceTieBreak::Strict
    )]
    fn test_aromaticity_resolver_select_with_report_identity(
        aromaticity_model: AromaticityModel,
        #[case] molecule: Molecule,
        #[case] completions: AtomCompletions,
        #[case] tie_break: ValenceTieBreak,
    ) {
        let state = ResolveState {
            completions,
            ..ResolveState::default()
        };
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).select_with_report(
                &molecule,
                state.clone(),
                tie_break
            ),
            Ok(Solution::Determined(state))
        );
    }

    #[rstest]
    #[case::shared_restriction(
        AtomCompletions::from_iter([
            (AtomId(8), smallvec![atom_dsl!("N#c0#a0"), atom_dsl!("N#c0#a2")]),
        ]),
        ValenceTieBreak::Strict,
        AtomCompletions::from_iter([
            (AtomId(8), smallvec![atom_dsl!("N#c0#a0"), atom_dsl!("N#c0#a2")]),
        ]),
        vec![((0..6).map(AtomId).collect(),
            AromaticSystemForm::from_electrons(vec![1; 6])
                .with_charge(0)
                .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()))],
        Vec::new()
    )]
    #[case::distinct_restrictions(
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![atom_dsl!("C#i=#c0#h1#n0#u0#s#v2#a0"), atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2")]),
            (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h1#n0#u0#s#v2#a0"), atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2")]),
            (AtomId(8), smallvec![atom_dsl!("N#c0#a0"), atom_dsl!("N#c0#a2")]),
        ]),
        ValenceTieBreak::MostSaturated,
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![atom_dsl!("C#i=#c0#h1#n0#u0#s#v2#a0")]),
            (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2")]),
            (AtomId(8), smallvec![atom_dsl!("N#c0#a0"), atom_dsl!("N#c0#a2")]),
        ]),
        vec![((0..6).map(AtomId).collect(),
            AromaticSystemForm::from_electrons(vec![0, 2, 1, 1, 1, 1])
                .with_charge(0)
                .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()))],
        (0..6).map(AtomId).collect()
    )]
    fn test_aromaticity_resolver_select_with_report_duplicates(
        quinoline: Molecule,
        #[case] completions: AtomCompletions,
        #[case] tie_break: ValenceTieBreak,
        #[case] expected_completions: AtomCompletions,
        #[case] aromatic_systems: Vec<(Vec<AtomId>, AromaticSystemForm)>,
        #[case] tie_breaks: Vec<AtomId>,
        #[values(
            AtomSelectionAlgorithm::FullDfs,
            AtomSelectionAlgorithm::BacktrackingDfs
        )]
        atom_selection_algorithm: AtomSelectionAlgorithm,
    ) {
        let model = AromaticityModel::daylight();
        let resolver = AromaticityResolver::with_config(
            &model,
            AromaticityResolveConfig {
                atom_selection_algorithm,
                aromatic_valence_failure: AromaticityFailurePolicy::Keep,
                ..Default::default()
            },
        );
        assert_eq!(
            resolver.select_with_report(
                &quinoline,
                ResolveState {
                    completions,
                    ..Default::default()
                },
                tie_break,
            ),
            Ok(Solution::Determined(ResolveState {
                completions: expected_completions,
                aromatic_systems,
                tie_breaks,
            }))
        );
    }

    #[rstest]
    #[case::strict(ValenceTieBreak::Strict)]
    fn test_aromaticity_resolver_select_with_report_duplicates_identity(
        quinoline: Molecule,
        #[case] tie_break: ValenceTieBreak,
        #[values(
            AtomSelectionAlgorithm::FullDfs,
            AtomSelectionAlgorithm::BacktrackingDfs
        )]
        atom_selection_algorithm: AtomSelectionAlgorithm,
    ) {
        let state = ResolveState {
            completions: AtomCompletions::from_iter([
                (
                    AtomId(0),
                    smallvec![
                        atom_dsl!("C#i=#c0#h1#n0#u0#s#v2#a0"),
                        atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2")
                    ],
                ),
                (
                    AtomId(1),
                    smallvec![
                        atom_dsl!("C#i=#c0#h1#n0#u0#s#v2#a0"),
                        atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2")
                    ],
                ),
                (
                    AtomId(8),
                    smallvec![atom_dsl!("N#c0#a0"), atom_dsl!("N#c0#a2")],
                ),
            ]),
            ..Default::default()
        };
        let model = AromaticityModel::daylight();
        let resolver = AromaticityResolver::with_config(
            &model,
            AromaticityResolveConfig {
                atom_selection_algorithm,
                aromatic_valence_failure: AromaticityFailurePolicy::Keep,
                ..Default::default()
            },
        );
        assert_eq!(
            resolver.select_with_report(&quinoline, state.clone(), tie_break),
            Ok(Solution::Determined(state))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_with_report_min_electron_count() {
        // Two totals pass the rule on the same members (all-pyridinic 6,
        // all-pyrrolic 10); the electron component picks the smaller and
        // records the members, with no value key consulted.
        let model = AromaticityModel {
            scope: ElementScope::Any,
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
            tie_break: AromaticityTieBreak::MinElectronCount,
        };
        let molecule = mol_dsl!(
            r#"{:atoms ["N#c0" "N#c0" "C#c0" "N#c0" "N#c0" "C#c0"]
                :bonds [[0 1 "1#a+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#
        );
        assert_eq!(
            AromaticityResolver::new(&model).select_with_report(
                &molecule,
                ResolveState {
                    completions: AtomCompletions::from_iter([
                        (
                            AtomId(0),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (
                            AtomId(1),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                        (
                            AtomId(3),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (
                            AtomId(4),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    ]),
                    ..ResolveState::default()
                },
                ValenceTieBreak::Strict,
            ),
            Ok(Solution::Determined(ResolveState {
                completions: AtomCompletions::from_iter([
                    (AtomId(0), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                    (AtomId(1), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                    (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    (AtomId(3), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                    (AtomId(4), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                    (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                ]),
                aromatic_systems: vec![(
                    (0..6).map(AtomId).collect(),
                    AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 1, 1])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                )],
                tie_breaks: (0..6).map(AtomId).collect(),
            }))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_with_report_tolerated_carrier() {
        // Keep policy admits the combination without a system alongside the full
        // ring; `MinElectronCount` realizes the full ring and records its members.
        let model = AromaticityModel {
            scope: ElementScope::Any,
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
            tie_break: AromaticityTieBreak::MinElectronCount,
        };
        let molecule = mol_dsl!(
            r#"{:atoms ["N#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                :bonds [[0 1 "1#a+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#
        );
        assert_eq!(
            AromaticityResolver::with_config(
                &model,
                AromaticityResolveConfig {
                    aromatic_valence_failure: AromaticityFailurePolicy::Keep,
                    ..AromaticityResolveConfig::default()
                },
            )
            .select_with_report(
                &molecule,
                ResolveState {
                    completions: AtomCompletions::from_iter([
                        (
                            AtomId(0),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                        (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                        (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                        (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    ]),
                    ..ResolveState::default()
                },
                ValenceTieBreak::Strict,
            ),
            Ok(Solution::Determined(ResolveState {
                completions: AtomCompletions::from_iter([
                    (AtomId(0), smallvec![atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2")]),
                    (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                ]),
                aromatic_systems: vec![(
                    (0..5).map(AtomId).collect(),
                    AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                )],
                tie_breaks: (0..5).map(AtomId).collect(),
            }))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_with_report_tolerated_carrier_identity(
        aromaticity_model: AromaticityModel,
    ) {
        // The same survivors under the model's `Strict`: structurally
        // distinct, so the members stay plural and the state passes through.
        let molecule = mol_dsl!(
            r#"{:atoms ["N#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                :bonds [[0 1 "1#a+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#
        );
        let state = ResolveState {
            completions: AtomCompletions::from_iter([
                (
                    AtomId(0),
                    smallvec![
                        atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                        atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                    ],
                ),
                (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            ]),
            ..ResolveState::default()
        };
        assert_eq!(
            AromaticityResolver::with_config(
                &aromaticity_model,
                AromaticityResolveConfig {
                    aromatic_valence_failure: AromaticityFailurePolicy::Keep,
                    ..AromaticityResolveConfig::default()
                },
            )
            .select_with_report(&molecule, state.clone(), ValenceTieBreak::Strict),
            Ok(Solution::Determined(state))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_with_report_tie_breaks(benzene: Molecule) {
        let resolver = AromaticityResolver::new(&AromaticityModel::daylight());
        let input = ResolveState {
            tie_breaks: vec![AtomId(2), AtomId(0), AtomId(2)],
            ..Default::default()
        };
        let mut expected = ResolveState {
            aromatic_systems: vec![(
                (0..6).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![1; 6])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )],
            ..input.clone()
        };
        assert_eq!(
            resolver.select(&benzene, input.clone(), ValenceTieBreak::Strict),
            Ok(Solution::Determined(expected.clone()))
        );
        expected.tie_breaks = vec![AtomId(0), AtomId(2)];
        assert_eq!(
            resolver.select_with_report(&benzene, input, ValenceTieBreak::Strict),
            Ok(Solution::Determined(expected))
        );
    }

    #[rstest]
    #[case::bond_between_aromatic_systems(
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 6 "1"] [0 6 "1#a"]]
        }"#),
        AtomId(0), smallvec![atom_dsl!("C#a"), atom_dsl!("C#a2")],
        AromaticityInconsistency::AromaticBondConstraintFailure { bond: BondId(12) }
    )]
    #[case::stored_aromatic_system_outside_search(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[3 4 "1"] [4 5 "1"] [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 3 "1"]]
            :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,0]"}]
        }"#),
        AtomId(3), smallvec![atom_dsl!("C#a"), atom_dsl!("C#a2")],
        AromaticityInconsistency::AromaticSystemFailure { aromatic_system: AromaticSystemId(0) }
    )]
    #[case::negative_bond_before_pruning(
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        AtomId(0), smallvec![atom_dsl!("C#a"), atom_dsl!("C#a2")],
        AromaticityInconsistency::AromaticBondConstraintMismatch { bond: BondId(0), aromatic_system: AromaticSystemId(0) }
    )]
    fn test_aromaticity_resolver_select_with_report_failure(
        #[case] molecule: Molecule,
        #[case] atom: AtomId,
        #[case] forms: SmallVec<[AtomForm; 1]>,
        #[case] expected: AromaticityInconsistency,
        #[values(
            AtomSelectionAlgorithm::FullDfs,
            AtomSelectionAlgorithm::BacktrackingDfs
        )]
        atom_selection_algorithm: AtomSelectionAlgorithm,
    ) {
        let resolver = AromaticityResolver::with_config(
            &AromaticityModel::daylight(),
            AromaticityResolveConfig {
                atom_selection_algorithm,
                ..Default::default()
            },
        );
        let state = ResolveState {
            completions: AtomCompletions::from_iter([(atom, forms)]),
            ..Default::default()
        };
        assert_eq!(
            resolver.select_with_report(&molecule, state, ValenceTieBreak::Strict),
            Ok(Solution::Contradictory(expected.into()))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_with_report_stored_conflict(
        aromaticity_model: AromaticityModel,
    ) {
        // No completion of the flexible atom reproduces the stored system:
        // `#a2` breaks the count, `#a!` removes the candidate ring.
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#c0" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        );
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).select_with_report(
                &molecule,
                ResolveState {
                    completions: AtomCompletions::from_iter([(
                        AtomId(0),
                        smallvec![
                            atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
                            atom_dsl!("C#i=#c0#h2#n0#u0#s#v2#a!"),
                        ]
                    )]),
                    ..ResolveState::default()
                },
                ValenceTieBreak::Strict,
            ),
            Ok(Solution::Contradictory(
                AromaticityContradiction::Inconsistency(
                    AromaticityInconsistency::AromaticSystemFailure {
                        aromatic_system: AromaticSystemId(0)
                    }
                )
            ))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_with_report_stored_conflict_identity(
        aromaticity_model: AromaticityModel,
    ) {
        // Same conflict under `Keep`: the component is inert and the state
        // passes through unchanged.
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#c0" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        );
        let state = ResolveState {
            completions: AtomCompletions::from_iter([(
                AtomId(0),
                smallvec![
                    atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
                    atom_dsl!("C#i=#c0#h2#n0#u0#s#v2#a!"),
                ],
            )]),
            ..ResolveState::default()
        };
        assert_eq!(
            AromaticityResolver::with_config(
                &aromaticity_model,
                AromaticityResolveConfig {
                    aromatic_system_failure: AromaticityFailurePolicy::Keep,
                    ..AromaticityResolveConfig::default()
                },
            )
            .select_with_report(&molecule, state.clone(), ValenceTieBreak::Strict),
            Ok(Solution::Determined(state))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_resolve_into_identity(
        aromaticity_model: AromaticityModel,
        mut benzene: Molecule,
    ) {
        let resolver = AromaticityResolver::new(&aromaticity_model);
        assert_eq!(
            resolver.resolve_into(&mut benzene),
            Ok(Solution::Determined(()))
        );
        let expected = benzene.clone();

        assert_eq!(
            resolver.resolve_into(&mut benzene),
            Ok(Solution::Determined(()))
        );
        assert_eq!(benzene, expected);
    }

    #[rstest]
    #[case::clar_heterocycle(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar, tie_break: AromaticityTieBreak::Strict },
        mol_dsl_concrete!(r#"{:atoms ["N #h #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        AromaticityContradiction::ClarNonBenzenoid(
            "Clar model requires benzenoid input but non-carbon aromatic atoms are present".to_string()
        )
    )]
    #[case::aromatic_valence_failure(
        AromaticityModel::mdl(),
        mol_dsl_concrete!(r#"{:atoms ["O #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }
        )
    )]
    #[case::aromatic_system_failure(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]"}]
        }"#),
        AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticSystemFailure {
                aromatic_system: AromaticSystemId(0)
            }
        )
    )]
    fn test_aromaticity_resolver_resolve_contradiction(
        #[case] model: AromaticityModel,
        #[case] mut molecule: Molecule,
        #[case] expected: AromaticityContradiction,
    ) {
        let original = molecule.clone();
        assert_eq!(
            AromaticityResolver::new(&model).resolve(molecule.clone()),
            Ok(Solution::Contradictory(expected.clone()))
        );
        assert_eq!(
            AromaticityResolver::new(&model).resolve_into(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
    }
}
