//! Aromaticity resolver. Perception reads aromatic valence constraints from the
//! materialized valence stage; planning emits complete aromatic systems and
//! discharges the assertions they realize without mutating the source molecule.

use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::slice;

use smallvec::{smallvec, SmallVec};
use thiserror::Error;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemHandle, AromaticSystemId, AromaticValence,
    AromaticValenceForm, AtomConstraintForm, AtomConstraintKey, AtomForm, AtomHandle, AtomId,
    BondConstraintForm, BondConstraintKey, BondHandle, BondId, BooleanForm, Edits,
    ElectronCountsForm, Lattice, Molecule, Normalize, NumForm, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{
    read_aromatic_valence, AromaticityConfig, AromaticityContradiction, AromaticityError,
    AromaticityInconsistency, AromaticityPerceiver,
};
use crate::ops::model::{AromaticityModel, AromaticityTieBreak, ElementScope, ValenceTieBreak};
use crate::ops::resolve::{ResolveState, WildcardPolicy};
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
    #[error(
        "aromatic system {aromatic_system:?} has {member_count} members but {electron_count} electron contributions"
    )]
    ElectronCountMismatch {
        aromatic_system: AromaticSystemId,
        member_count: usize,
        electron_count: usize,
    },
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
    /// Atom aromatic-valence and bond-aromaticity assertions are discharged when
    /// independently entailed by retained or newly added systems. Constraints that
    /// add information, including negative assertions outside systems, are preserved.
    /// A required Hückel contribution may be supplied by a wildcard atom in scope.
    /// Missing contributions and unresolved charge or spin on a retained system leave the
    /// plan underdetermined. Standalone resolution can retain an open electron vector.
    /// Member assertions that supply its contributions remain. Aggregate resolution
    /// requires [WildcardPolicy::Admit] to return such a system.
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
                let mut edits = self.plan_discharge(molecule, &state.aromatic_systems);
                if self.plan_systems(
                    molecule,
                    &state.aromatic_systems,
                    WildcardPolicy::Admit,
                    &mut edits,
                ) {
                    Ok(Solution::Determined(edits))
                } else {
                    Ok(Solution::Underdetermined(Edits::new()))
                }
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
            molecule.mutate(edits);
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
    /// non-singlet or nonzero-unpaired system spin; mismatched member/contribution counts;
    /// and incompatible atom or bond assertions.
    pub fn project_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), AromaticityContradiction>, AromaticityProjectError> {
        let edits = self.plan_project(molecule)?;
        if !edits.is_empty() {
            molecule.mutate(edits);
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
    /// singlet, mismatched member/contribution counts, and conflicting assertions.
    pub fn plan_project(&self, molecule: &Molecule) -> Result<Edits, AromaticityProjectError> {
        if !molecule.has_aromatic_systems() {
            return Ok(Edits::new());
        }
        let mut atom_constraints = Vec::new();
        let mut bond_constraints = Vec::new();
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
            let member_count = aromatic_system.atom_ids().len();
            if electrons.len() != member_count {
                return Err(AromaticityProjectError::ElectronCountMismatch {
                    aromatic_system: aromatic_system.id(),
                    member_count,
                    electron_count: electrons.len(),
                });
            }
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
                let new = AtomConstraintForm::aromatic_valence(assertion);
                let old = attributes
                    .constraints
                    .get(AtomConstraintKey::AromaticValence);
                if old.is_none_or(|old| !old.normalized_eq(&new)) {
                    atom_constraints.push((AtomHandle::Id(atom), old.cloned(), Some(new)));
                }
            }
            for bond in aromatic_system.bond_ids() {
                let attributes = molecule.bond(bond).attributes();
                let assertion = attributes
                    .constraints
                    .aromatic()
                    .meet(&BooleanForm::Lit(true))
                    .ok_or(AromaticityProjectError::BondAssertion { bond })?;
                let new = BondConstraintForm::Aromatic(assertion);
                let old = attributes.constraints.get(BondConstraintKey::Aromatic);
                if old.is_none_or(|old| !old.normalized_eq(&new)) {
                    bond_constraints.push((BondHandle::Id(bond), old.cloned(), Some(new)));
                }
            }
        }
        let mut edits = Edits::new();
        if !atom_constraints.is_empty() {
            edits.modify_atom_constraints(AtomConstraintKey::AromaticValence, atom_constraints);
        }
        if !bond_constraints.is_empty() {
            edits.modify_bond_constraints(BondConstraintKey::Aromatic, bond_constraints);
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
    /// A missing required contribution, undetermined perception, or more than 4096
    /// perception calls in a component yields
    /// Underdetermined. A retained system keeps its electron vector. Independent evidence
    /// can supply the count for an open contribution. An established failure to cover a
    /// required atom yields Contradictory under the Error failure policy. If no positive aromatic
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

    /// Adds missing systems and checks the attributes of each retained or new system.
    pub(super) fn plan_systems(
        &self,
        molecule: &Molecule,
        systems: &[(Vec<AtomId>, AromaticSystemForm)],
        wildcards: WildcardPolicy,
        edits: &mut Edits,
    ) -> bool {
        let accepts = |system: &AromaticSystemForm| {
            wildcards.admits(system.electrons.is_ground(), system.is_wildcard())
                && system.charge.is_ground()
                && system.unpaired_electrons.is_ground()
        };
        let mut determined = true;
        let mut existing = BTreeSet::new();
        for system in molecule.aromatic_systems().iter() {
            determined &= accepts(system.attributes());
            let mut atoms: Vec<_> = system.atom_ids().collect();
            atoms.sort_unstable();
            existing.insert(atoms);
        }
        for (atoms, system) in systems {
            let mut key = atoms.clone();
            key.sort_unstable();
            if existing.contains(&key) {
                continue;
            }
            determined &= accepts(system);
            for edit in
                self.perception
                    .plan_aromatic_system(molecule, atoms.clone(), system.clone())
            {
                edits.push(edit);
            }
        }
        determined
    }

    pub(crate) fn plan_discharge(
        &self,
        molecule: &Molecule,
        aromatic_systems: &[(Vec<AtomId>, AromaticSystemForm)],
    ) -> Edits {
        let mut edits = Edits::new();
        let mut atoms = Vec::new();
        for atom in molecule.atoms().iter() {
            let key = AtomConstraintKey::AromaticValence;
            let Some(asserted) = atom.attributes().constraints.get(key) else {
                continue;
            };
            if atom
                .constraints()
                .derived(key)
                .is_some_and(|derived| asserted.matches(&derived))
            {
                atoms.push((AtomHandle::Id(atom.id()), Some(asserted.clone()), None));
            }
        }
        if !atoms.is_empty() {
            edits.modify_atom_constraints(AtomConstraintKey::AromaticValence, atoms);
        }
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
                Some((bond, asserted))
            })
            .peekable();
        if bonds.peek().is_none() {
            return edits;
        }
        let mut membership = vec![None; molecule.atoms().count()];
        for aromatic_system in molecule.aromatic_systems().iter() {
            for atom in aromatic_system.atom_ids() {
                membership[atom.index()] = Some(aromatic_system.id().index());
            }
        }
        for (aromatic_system, (atoms, _)) in aromatic_systems.iter().enumerate() {
            for atom in atoms {
                if let Some(member) = membership.get_mut(atom.index()) {
                    *member = Some(molecule.aromatic_systems().count() + aromatic_system);
                }
            }
        }
        let mut changes = Vec::new();
        for (bond, asserted) in bonds {
            let [a, b] = bond.atom_ids();
            let aromatic =
                membership[a.index()].is_some() && membership[a.index()] == membership[b.index()];
            if aromatic && asserted.matches(&BooleanForm::Lit(true)) {
                changes.push((
                    BondHandle::Id(bond.id()),
                    Some(BondConstraintForm::Aromatic(*asserted)),
                    None,
                ));
            }
        }
        if !changes.is_empty() {
            edits.modify_bond_constraints(BondConstraintKey::Aromatic, changes);
        }
        edits
    }

    fn select_completions(
        &self,
        molecule: &Molecule,
        mut state: ResolveState,
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

        let completions = &mut state.completions;
        let mut fixed = vec![AromaticValence::NotAromatic; molecule.atoms().count()];
        let any_element_hueckel = match &self.perception {
            AromaticityPerceiver::Hueckel(model) if model.element_scope == ElementScope::Any => {
                Some(model)
            }
            _ => None,
        };
        let keep_failure = |failure| match failure {
            AromaticityInconsistency::AromaticSystemFailure { .. } => {
                self.config.aromatic_system_failure == AromaticityFailurePolicy::Keep
            }
            AromaticityInconsistency::AromaticValenceFailure { .. } => {
                self.config.aromatic_valence_failure == AromaticityFailurePolicy::Keep
            }
            _ => false,
        };
        let mut missing = BTreeSet::new();
        let mut flexible = Vec::new();
        let mut underdetermined = false;
        for atom in molecule.atoms().ids() {
            if let Some(forms) = completions.get(atom) {
                let mut options = Vec::with_capacity(forms.len());
                for form in forms {
                    let evidence = form.constraints.aromatic_valence();
                    let value = match evidence {
                        AromaticValenceForm::Aromatic(NumForm::Lit(value)) => {
                            AromaticValence::Aromatic(*value)
                        }
                        AromaticValenceForm::Aromatic(_) => {
                            match read_aromatic_valence(molecule, atom) {
                                Ok(Some(AromaticValence::Aromatic(value))) => {
                                    AromaticValence::Aromatic(value)
                                }
                                Ok(_) => {
                                    underdetermined = true;
                                    AromaticValence::NotAromatic
                                }
                                Err(failure) if keep_failure(failure) => {
                                    AromaticValence::NotAromatic
                                }
                                Err(failure) => return Ok(Solution::Contradictory(failure.into())),
                            }
                        }
                        _ => AromaticValence::NotAromatic,
                    };
                    options.push(value);
                }
                fixed[atom.index()] = options[0];
                if options.iter().any(|&value| value != options[0]) {
                    flexible.push((atom, options));
                }
                continue;
            }
            match read_aromatic_valence(molecule, atom) {
                Ok(Some(value)) => {
                    fixed[atom.index()] = value;
                    if value == AromaticValence::NotAromatic && any_element_hueckel.is_some() {
                        let form = molecule.atom(atom).attributes();
                        if form.is_wildcard()
                            && form.constraints.aromatic_valence().is_undetermined()
                        {
                            missing.insert(atom);
                        }
                    }
                }
                Ok(None) => match &self.perception {
                    AromaticityPerceiver::Hueckel(rule)
                        if molecule.atom(atom).attributes().is_wildcard() =>
                    {
                        if rule.element_scope == ElementScope::Any {
                            missing.insert(atom);
                        }
                    }
                    _ => underdetermined = true,
                },
                Err(failure) if keep_failure(failure) => {}
                Err(inconsistency) => {
                    return Ok(Solution::Contradictory(inconsistency.into()));
                }
            }
        }
        if let Some(inconsistency) = self
            .perception
            .check_aromatic_bonds(molecule)
            .into_iter()
            .next()
        {
            return Ok(Solution::Contradictory(inconsistency.into()));
        }
        if underdetermined {
            return Ok(Solution::Underdetermined(state));
        }
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
        // Only a request involving an eligible ring can depend on a bare wildcard.
        let incomplete_rings = any_element_hueckel.map_or_else(Vec::new, |model| {
            model.systems_missing_contributions(molecule, &rings, &missing, |atom| {
                completions.get(atom).map_or_else(
                    || matches!(fixed[atom.index()], AromaticValence::Aromatic(_)),
                    |forms| {
                        forms
                            .iter()
                            .any(|form| form.constraints.aromatic_valence().is_aromatic())
                    },
                )
            })
        });
        let needs_contribution = |failure: AromaticityInconsistency| {
            incomplete_rings.iter().any(|members| match failure {
                AromaticityInconsistency::AromaticValenceFailure { atom } => {
                    members.contains(&atom)
                }
                AromaticityInconsistency::AromaticBondConstraintFailure { bond } => molecule
                    .bond(bond)
                    .atom_ids()
                    .iter()
                    .all(|atom| members.contains(atom)),
                AromaticityInconsistency::AromaticSystemFailure { aromatic_system } => {
                    let system = molecule.aromatic_system(aromatic_system);
                    system.charge().is_compatible(&NumForm::Lit(0))
                        && system
                            .unpaired_electrons()
                            .is_compatible(&UnpairedElectronsForm::closed_shell())
                        && system.atom_ids().len() == members.len()
                        && system.atom_ids().all(|atom| members.contains(&atom))
                }
                _ => false,
            })
        };
        let finish = |mut state: ResolveState,
                      mut accepted: Vec<(Vec<AtomId>, AromaticSystemForm)>,
                      tie_break_uses: BTreeSet<AtomId>,
                      underdetermined: bool| {
            if underdetermined {
                return Solution::Underdetermined(state);
            }
            for system in molecule.aromatic_systems().iter() {
                if system.attributes().is_wildcard() {
                    if let Some((_, form)) = accepted.iter_mut().find(|(atoms, _)| {
                        atoms.len() == system.atom_ids().len()
                            && system
                                .atom_ids()
                                .all(|atom| atoms.binary_search(&atom).is_ok())
                    }) {
                        form.electrons = ElectronCountsForm::Undetermined;
                    }
                }
            }
            state.aromatic_systems.extend(accepted);
            if report {
                state.tie_breaks.extend(tie_break_uses);
                state.tie_breaks.sort_unstable();
                state.tie_breaks.dedup();
            }
            Solution::Determined(state)
        };
        if flexible.is_empty() {
            let mut derivation = match self.perception.derive_from_rings(
                molecule,
                &rings,
                self.config.perception,
                &required_atoms,
                |atom| fixed[atom.index()],
            )? {
                Solution::Determined(derivation) => derivation,
                Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(state)),
                Solution::Contradictory(contradiction) => {
                    return Ok(Solution::Contradictory(contradiction))
                }
            };
            let mut failed_atoms = BTreeSet::new();
            for inconsistency in derivation.inconsistencies {
                match inconsistency {
                    failure if needs_contribution(failure) => underdetermined = true,
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
                        matches!(fixed[atom.index()], AromaticValence::Aromatic(_))
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
            return Ok(finish(
                state,
                derivation.aromatic_systems,
                BTreeSet::new(),
                underdetermined,
            ));
        }
        let components = self
            .perception
            .candidate_components(molecule, &rings, |atom_id| match completions.get(atom_id) {
                Some(disjuncts) => disjuncts
                    .iter()
                    .any(|form| form.constraints.aromatic_valence().is_aromatic()),
                None => matches!(fixed[atom_id.index()], AromaticValence::Aromatic(_)),
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
                (self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error
                    || !incomplete_rings.is_empty())
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
                .filter(|atom| {
                    component.contains(atom)
                        && !incomplete_rings.iter().any(|ring| ring.contains(atom))
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
                    None => slice::from_ref(&fixed[atom.index()]),
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
                        .any(|form| form.constraints.aromatic_valence().is_aromatic()),
                    None => matches!(fixed[atom.index()], AromaticValence::Aromatic(_)),
                }
            };
            let require_aromaticity =
                self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error;
            let prune_enabled = require_aromaticity
                && !aromatic_only.is_empty()
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
                        return Ok(Solution::Underdetermined(state));
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
                        |atom| match flexible_positions.get(&atom) {
                            Some(&position) => component_flexible[position].1[path[position]],
                            None => fixed[atom.index()],
                        },
                    )?;
                    let found = match outcome {
                        Solution::Determined(found) => found,
                        Solution::Underdetermined(_) => {
                            return Ok(Solution::Underdetermined(state));
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
                            failure if needs_contribution(failure) => {}
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
        if ((self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error
            || !incomplete_rings.is_empty())
            && (!missing_atoms.is_empty() || claimed_bonds.iter().any(|claimed| !claimed)))
            || molecule.aromatic_systems().ids().any(outside_search)
            || (!incomplete_rings.is_empty() && molecule.has_aromatic_systems())
        {
            let inconsistencies =
                self.perception
                    .check_aromatic_systems(molecule, &missing_atoms, &accepted);
            for inconsistency in inconsistencies {
                if let AromaticityInconsistency::AromaticBondConstraintFailure { bond } =
                    inconsistency
                {
                    if !required_bonds
                        .iter()
                        .zip(&claimed_bonds)
                        .any(|((required, _), claimed)| *required == bond && !claimed)
                    {
                        continue;
                    }
                }
                let reject = match inconsistency {
                    failure if needs_contribution(failure) => {
                        underdetermined = true;
                        false
                    }
                    AromaticityInconsistency::AromaticValenceFailure { .. }
                    | AromaticityInconsistency::AromaticBondConstraintFailure { .. } => {
                        self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error
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

        Ok(finish(state, accepted, tie_break_uses, underdetermined))
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
    use umol_chem::element::Element;
    use umol_graph_core::{
        ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
        RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
    };
    use umol_graph_ir::ir::{
        AddAromaticSystem, AromaticSystemId, AtomConstraintKey, BondId, Edit, Edits, NumForm,
        RingConfig, RollbackMethod, UnpairedElectronsForm,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_concrete};

    use super::*;
    use crate::ops::model::{
        AromaticityRule, AromaticityTieBreak, ChemistryModel, ElementScope, RingLimits,
    };
    use crate::ops::resolve::{
        ProjectError, ProjectFlags, ResolveConfig, ResolveContradiction, Resolver,
    };
    use crate::ops::valence::{AtomCompletions, AtomTypingError, ValenceContradiction};

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
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}]
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
                Edit::ModifyAtomConstraints {
                    key: AtomConstraintKey::AromaticValence,
                    changes: (0..6)
                        .map(|id| (
                            AtomHandle::Id(AtomId(id)),
                            Some(AtomConstraintForm::AromaticValence(
                                AromaticValenceForm::Aromatic(NumForm::Lit(1))
                            )),
                            None,
                        ))
                        .collect(),
                },
            ])))
        );
    }

    #[rstest]
    #[case::bare(
        AromaticityModel::permissive(),
        "*",
        AromaticityFailurePolicy::Error,
        Solution::Underdetermined(Edits::new())
    )]
    #[case::bare_keep(
        AromaticityModel::permissive(),
        "*",
        AromaticityFailurePolicy::Keep,
        Solution::Underdetermined(Edits::new())
    )]
    #[case::open(
        AromaticityModel::permissive(),
        "*#a+",
        AromaticityFailurePolicy::Error,
        Solution::Underdetermined(Edits::new())
    )]
    #[case::open_keep(
        AromaticityModel::permissive(),
        "*#a+",
        AromaticityFailurePolicy::Keep,
        Solution::Underdetermined(Edits::new())
    )]
    #[case::daylight(AromaticityModel::daylight(), "*#a1", AromaticityFailurePolicy::Error,
        Solution::Contradictory(AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }.into()))]
    #[case::daylight_open(AromaticityModel::daylight(), "*#a+", AromaticityFailurePolicy::Error,
        Solution::Contradictory(AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }.into()))]
    #[case::daylight_keep(
        AromaticityModel::daylight(),
        "*#a1",
        AromaticityFailurePolicy::Keep,
        Solution::Determined(Edits::new())
    )]
    #[case::excluded(AromaticityModel::permissive(), "*#a!", AromaticityFailurePolicy::Error,
        Solution::Contradictory(AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(1) }.into()))]
    #[case::hmo(AromaticityModel { rule: AromaticityRule::Hmo { stabilization_threshold: 0.0 }, ..AromaticityModel::permissive() },
        "*#a1", AromaticityFailurePolicy::Error, Solution::Contradictory(AromaticityContradiction::HmoInvalidInput(
            "open-shell pi-system (odd electron count) not supported by HMO".to_string())))]
    #[case::clar(AromaticityModel { rule: AromaticityRule::Clar, ..AromaticityModel::permissive() },
        "*#a1", AromaticityFailurePolicy::Error, Solution::Contradictory(AromaticityContradiction::ClarNonBenzenoid(
            "Clar model requires benzenoid input but non-carbon aromatic atoms are present".to_string())))]
    fn test_aromaticity_resolver_plan_resolve_wildcard(
        mut benzene: Molecule,
        #[case] model: AromaticityModel,
        #[case] member: &str,
        #[case] failure: AromaticityFailurePolicy,
        #[case] expected: Solution<Edits, AromaticityContradiction>,
    ) {
        *benzene.atom_mut(AtomId(0)).attributes_mut() = atom_dsl!(member);
        let resolver = AromaticityResolver::with_config(
            &model,
            AromaticityResolveConfig {
                aromatic_valence_failure: failure,
                ..Default::default()
            },
        );
        assert_eq!(resolver.plan_resolve(&benzene), Ok(expected));
    }

    #[rstest]
    #[case::error(AromaticityFailurePolicy::Error, Solution::Contradictory(
        AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(6) }.into()))]
    #[case::keep(
        AromaticityFailurePolicy::Keep,
        Solution::Underdetermined(Edits::new())
    )]
    fn test_aromaticity_resolver_plan_resolve_partial_failure(
        #[case] failure: AromaticityFailurePolicy,
        #[case] expected: Solution<Edits, AromaticityContradiction>,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["*" "C#a1" "C#a1" "C#a1" "C#a1" "C#a1" "C#a1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#
        );
        let resolver = AromaticityResolver::with_config(
            &AromaticityModel::permissive(),
            AromaticityResolveConfig {
                aromatic_valence_failure: failure,
                ..Default::default()
            },
        );
        assert_eq!(resolver.plan_resolve(&molecule), Ok(expected));
    }

    #[rstest]
    #[case::open_contribution(
        mol_dsl!(r#"{
            :atoms ["C#a+" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
        }"#),
        Solution::Underdetermined(Edits::new()),
    )]
    #[case::bond_mismatch_unknown_first(
        mol_dsl!(r#"{
            :atoms ["C#a+" "C" "C" "C" "C" "C" "C"]
            :bonds [[1 2 "1#a!"] [2 3 "1"] [3 4 "1"]
                    [4 5 "1"] [5 6 "1"] [6 1 "1"]]
            :aromatic-systems [{:atoms [1 2 3 4 5 6] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        Solution::Contradictory(AromaticityInconsistency::AromaticBondConstraintMismatch {
            bond: BondId(0),
            aromatic_system: AromaticSystemId(0),
        }.into()),
    )]
    #[case::bond_mismatch_unknown_last(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C#a+"]
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        Solution::Contradictory(AromaticityInconsistency::AromaticBondConstraintMismatch {
            bond: BondId(0),
            aromatic_system: AromaticSystemId(0),
        }.into()),
    )]
    #[case::stored_open_charge(
        mol_dsl!(r#"{:atoms ["*#a1" "C#a1" "C#a1" "C#a1" "C#a1" "C#a1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*#u0#s"}]}"#),
        Solution::Underdetermined(Edits::new()),
    )]
    #[case::stored_open_spin(
        mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0"}]}"#),
        Solution::Underdetermined(Edits::new()),
    )]
    #[case::stored_missing_member(
        mol_dsl!(r#"{:atoms ["*" "C#a1" "C#a1" "C#a1" "C#a1" "C#a1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*#c0#u0#s"}]}"#),
        Solution::Underdetermined(Edits::new()),
    )]
    fn test_aromaticity_resolver_plan_resolve_partial(
        aromaticity_model: AromaticityModel,
        #[case] mut molecule: Molecule,
        #[case] expected: Solution<Edits, AromaticityContradiction>,
        #[values(AromaticityFailurePolicy::Error, AromaticityFailurePolicy::Keep)]
        failure: AromaticityFailurePolicy,
    ) {
        let original = molecule.clone();
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_valence_failure: failure,
                aromatic_system_failure: failure,
                ..Default::default()
            },
        );
        assert_eq!(resolver.plan_resolve(&molecule), Ok(expected.clone()));
        assert_eq!(molecule, original);
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Ok(expected.map(|_| ()))
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
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}]
        }"#),
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}]
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
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]#c0#u0#s"}]
        }"#),
        mol_dsl!(r#"{
            :atoms ["O#n1" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]#c0#u0#s"}]
        }"#)
    )]
    fn test_aromaticity_resolver_plan_resolve_constraints(
        #[case] model: AromaticityModel,
        #[case] config: AromaticityResolveConfig,
        #[case] molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        let resolver = AromaticityResolver::with_config(&model, config);
        let Solution::Determined(edits) = resolver.plan_resolve(&molecule).unwrap() else {
            panic!("aromaticity plan must be determined");
        };
        assert_eq!(molecule.apply(edits), Ok(expected));
    }

    #[rstest]
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
    #[case::negative_assertions(
        AromaticityModel::daylight(),
        AromaticityResolveConfig::default(),
        mol_dsl_concrete!(r#"{:atoms ["C#h3#a!" "C#h3#a!"] :bonds [[0 1 "1#a!"]]}"#)
    )]
    #[case::undetermined_contributions(
        AromaticityModel::daylight(),
        AromaticityResolveConfig::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a1" "C#h#a1" "C#h#a1" "C#h#a1" "C#h#a1" "C#h#a1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*"}]
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
    fn test_aromaticity_resolver_resolve_exocyclic(
        #[values(AromaticityModel::permissive(), AromaticityModel::daylight())]
        model: AromaticityModel,
    ) {
        let molecule = mol_dsl!(
            r#"{:atoms ["C#a1" "C#a1" "C#a1" "C#a1" "C#a1" "C#a1" "*"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"] [0 6 "1"]]}"#
        );
        let expected = mol_dsl!(
            r#"{:atoms ["C" "C" "C" "C" "C" "C" "*"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"] [0 6 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}]}"#
        );
        assert_eq!(
            AromaticityResolver::new(&model).resolve(molecule),
            Ok(Solution::Determined(expected))
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
    fn test_aromaticity_resolver_resolve_overlap_constraints(
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
        let expected = mol_dsl!(
            r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}]
        }"#
        );
        assert_eq!(
            resolver.resolve(overlapping_aromatic_system.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        let mut output = overlapping_aromatic_system.clone();
        assert_eq!(
            resolver.resolve_into(&mut output),
            Ok(Solution::Determined(()))
        );
        assert_eq!(output, expected);
        assert_eq!(output.clone().edit().finish(), Ok(output));
    }

    #[rstest]
    #[case::short(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1]#c0#u0#s"}]}"#))]
    #[case::long(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,0,0]#c0#u0#s"}]}"#))]
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
    #[case::short(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1]#c0#u0#s"}]}"#))]
    #[case::long(mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,0,0]#c0#u0#s"}]}"#))]
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
                    [9 10 "1"] [10 11 "1"] [11 6 "1"] [0 6 "1#a!"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}
                               {:atoms [6 7 8 9 10 11] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::stored_aromatic_system(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a{1,2}#D2" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a#R+"] [1 2 "1#a"] [2 3 "1#a"]
                    [3 4 "1#a"] [4 5 "1#a"] [5 0 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#D2" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1#R+"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::stored_participant_order(
        mol_dsl_concrete!(r#"{
            :atoms ["N#h#a{1,2}" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [2 1 0 4 3] :attrs "[1,1,2,1,1]"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["N#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [2 1 0 4 3] :attrs "[1,1,2,1,1]"}]
        }"#)
    )]
    #[case::undetermined_contributions(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*"}]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*"}]
        }"#)
    )]
    #[case::wildcard_member(
        mol_dsl!(r#"{:atoms ["*#a1" "C#a1" "C#a1" "C#a1" "C#a1" "C#a1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        mol_dsl!(r#"{:atoms ["*" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}]}"#),
    )]
    #[case::stored_wildcard_member(
        mol_dsl!(r#"{:atoms ["*#a1" "C#a1" "C#a1" "C#a1" "C#a1" "C#a1"]
            :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*#c0#u0#s"}]}"#),
        mol_dsl!(r#"{:atoms ["*#a1" "C#a1" "C#a1" "C#a1" "C#a1" "C#a1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "*#c0#u0#s"}]}"#),
    )]
    fn test_aromaticity_resolver_resolve_discharge(
        aromaticity_model: AromaticityModel,
        #[case] molecule: Molecule,
        #[case] expected: Molecule,
        #[values(AromaticityFailurePolicy::Error, AromaticityFailurePolicy::Keep)]
        failure: AromaticityFailurePolicy,
    ) {
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_valence_failure: failure,
                aromatic_system_failure: failure,
                ..Default::default()
            },
        );
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        let mut actual = molecule;
        assert_eq!(
            resolver.resolve_into(&mut actual),
            Ok(Solution::Determined(()))
        );
        assert_eq!(actual, expected);
        assert_eq!(resolver.resolve(actual), Ok(Solution::Determined(expected)));
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
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
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
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
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
            resolver.resolve_into(&mut composite, RollbackMethod::Journal),
            Ok(expected_into.map_contradiction(ResolveContradiction::Aromaticity))
        );
        assert_eq!(composite, expected_molecule);
    }

    #[rstest]
    #[case::atom(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a2" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        AromaticityInconsistency::AromaticValenceMismatch { atom: AtomId(0), aromatic_system: AromaticSystemId(0) },
        ResolveContradiction::Valence(ValenceContradiction::AtomTyping(AtomTypingError::NoMatch {
            atom_id: AtomId(0), element: Element::C, charge: Some(0),
        }))
    )]
    #[case::bond(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        AromaticityInconsistency::AromaticBondConstraintMismatch { bond: BondId(0), aromatic_system: AromaticSystemId(0) },
        ResolveContradiction::Aromaticity(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticBondConstraintMismatch { bond: BondId(0), aromatic_system: AromaticSystemId(0) }
        ))
    )]
    fn test_aromaticity_resolver_resolve_mismatch(
        aromaticity_model: AromaticityModel,
        #[case] molecule: Molecule,
        #[case] inconsistency: AromaticityInconsistency,
        #[case] aggregate_contradiction: ResolveContradiction,
        #[values(AromaticityFailurePolicy::Error, AromaticityFailurePolicy::Keep)]
        failure: AromaticityFailurePolicy,
        #[values(RollbackMethod::Clone, RollbackMethod::Journal)] rollback: RollbackMethod,
    ) {
        let config = AromaticityResolveConfig {
            aromatic_valence_failure: failure,
            aromatic_system_failure: failure,
            ..Default::default()
        };
        let resolver = AromaticityResolver::with_config(&aromaticity_model, config);
        let contradiction = AromaticityContradiction::from(inconsistency);
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Contradictory(contradiction.clone()))
        );
        let mut actual = molecule.clone();
        assert_eq!(
            resolver.resolve_into(&mut actual),
            Ok(Solution::Contradictory(contradiction.clone()))
        );
        assert_eq!(actual, molecule);
        let model = ChemistryModel {
            aromaticity: aromaticity_model,
            ..Default::default()
        };
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                aromaticity: config,
                ..Default::default()
            },
        );
        let contradiction = aggregate_contradiction;
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Contradictory(contradiction.clone()))
        );
        assert_eq!(
            resolver.resolve_with_report(molecule.clone()),
            Ok(Solution::Contradictory(contradiction.clone()))
        );
        assert_eq!(
            resolver.resolve_into(&mut actual, rollback),
            Ok(Solution::Contradictory(contradiction.clone()))
        );
        assert_eq!(actual, molecule);
        assert_eq!(
            resolver.resolve_into_with_report(&mut actual, rollback),
            Ok(Solution::Contradictory(contradiction))
        );
        assert_eq!(actual, molecule);
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
    #[case::short_counts("[1,1,1,1,1]#c0#u0#s", AromaticityProjectError::ElectronCountMismatch {
        aromatic_system: AromaticSystemId(0), member_count: 6, electron_count: 5,
    })]
    #[case::long_counts("[1,1,1,1,1,1,1]#c0#u0#s", AromaticityProjectError::ElectronCountMismatch {
        aromatic_system: AromaticSystemId(0), member_count: 6, electron_count: 7,
    })]
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
            AromaticityResolver::new(&model).plan_project(&molecule),
            Err(expected.clone())
        );
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
    #[case::short_counts(5)]
    #[case::long_counts(7)]
    fn test_aromaticity_resolver_project_electron_counts(
        resolved_benzene: Molecule,
        #[case] electron_count: usize,
        #[values(RollbackMethod::Clone, RollbackMethod::Journal)] rollback_method: RollbackMethod,
    ) {
        let mut molecule = Molecule::combine_all([&resolved_benzene, &resolved_benzene]);
        molecule
            .aromatic_system_mut(AromaticSystemId(1))
            .attributes_mut()
            .electrons = ElectronCountsForm::Lit(vec![1; electron_count]);
        let original = molecule.clone();
        let expected = AromaticityProjectError::ElectronCountMismatch {
            aromatic_system: AromaticSystemId(1),
            member_count: 6,
            electron_count,
        };
        let model = ChemistryModel::default();
        let aromaticity = AromaticityResolver::new(&model.aromaticity);
        assert_eq!(aromaticity.plan_project(&molecule), Err(expected.clone()));
        assert_eq!(aromaticity.project(molecule.clone()), Err(expected.clone()));
        assert_eq!(
            aromaticity.project_into(&mut molecule),
            Err(expected.clone())
        );
        assert_eq!(molecule, original);

        let resolver = Resolver::new(&model);
        assert_eq!(
            resolver.project(molecule.clone(), ProjectFlags::all()),
            Err(ProjectError::Aromaticity(expected.clone()))
        );
        assert_eq!(
            resolver.project_into(&mut molecule, ProjectFlags::all(), rollback_method),
            Err(ProjectError::Aromaticity(expected))
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
    fn test_aromaticity_resolver_plan_project_constraints(
        #[values(RollbackMethod::Clone, RollbackMethod::Journal)] rollback: RollbackMethod,
    ) {
        let mut molecule = mol_dsl_concrete!(
            r#"{
            :atoms ["C#c+#h#D2" "C#h#a+" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1#a"] [2 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2] :attrs "[0,1,1]"}]}"#
        );
        let raw = AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(
            NumForm::LitSet(Box::new(BTreeSet::from([0]))),
        ));
        molecule
            .atom_mut(AtomId(0))
            .attributes_mut()
            .constraints
            .set(raw.clone());
        molecule
            .bond_mut(BondId(0))
            .attributes_mut()
            .constraints
            .set(BondConstraintForm::Aromatic(BooleanForm::Undetermined));
        let model = AromaticityModel::daylight();
        let edits = AromaticityResolver::new(&model)
            .plan_project(&molecule)
            .unwrap();
        assert_eq!(
            edits,
            Edits::from_iter([
                Edit::ModifyAtomConstraints {
                    key: AtomConstraintKey::AromaticValence,
                    changes: vec![
                        (
                            AtomHandle::Id(AtomId(1)),
                            Some(AtomConstraintForm::AromaticValence(
                                AromaticValenceForm::Aromatic(NumForm::Undetermined)
                            )),
                            Some(AtomConstraintForm::AromaticValence(
                                AromaticValenceForm::Aromatic(NumForm::Lit(1))
                            ))
                        ),
                        (
                            AtomHandle::Id(AtomId(2)),
                            None,
                            Some(AtomConstraintForm::AromaticValence(
                                AromaticValenceForm::Aromatic(NumForm::Lit(1))
                            ))
                        ),
                    ],
                },
                Edit::ModifyBondConstraints {
                    key: BondConstraintKey::Aromatic,
                    changes: vec![
                        (
                            BondHandle::Id(BondId(0)),
                            Some(BondConstraintForm::Aromatic(BooleanForm::Undetermined)),
                            Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true)))
                        ),
                        (
                            BondHandle::Id(BondId(2)),
                            None,
                            Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true)))
                        ),
                    ],
                },
                Edit::RemoveAromaticSystems {
                    removes: vec![(
                        AromaticSystemHandle::Id(AromaticSystemId(0)),
                        vec![
                            AtomHandle::Id(AtomId(0)),
                            AtomHandle::Id(AtomId(1)),
                            AtomHandle::Id(AtomId(2))
                        ],
                        AromaticSystemForm::new(ElectronCountsForm::Lit(vec![0, 1, 1]))
                            .with_charge(0)
                            .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                    )],
                },
            ])
        );
        let mut expected = mol_dsl_concrete!(
            r#"{
            :atoms ["C#c+#h#D2" "C#h#a1" "C#h#a1"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 0 "1#a"]]}"#
        );
        expected
            .atom_mut(AtomId(0))
            .attributes_mut()
            .constraints
            .set(raw);
        assert_eq!(molecule.clone().apply(edits.clone()), Ok(expected.clone()));
        assert_eq!(molecule.transact([edits], rollback), Ok(()));
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

    #[fixture]
    fn ring_completions() -> AtomCompletions {
        AtomCompletions::from_iter([(AtomId(1), smallvec![atom_dsl!("C#a1"), atom_dsl!("C#a2")])])
    }

    #[rstest]
    #[case::literal("*#a1", AromaticityFailurePolicy::Error, Solution::Determined(ResolveState {
        completions: AtomCompletions::from_iter([(AtomId(1), smallvec![atom_dsl!("C#a1")])]),
        aromatic_systems: vec![((0..6).map(AtomId).collect(), AromaticSystemForm::from_electrons(vec![1; 6]).with_charge(0).with_unpaired_electrons(UnpairedElectronsForm::closed_shell()))],
        tie_breaks: Vec::new(),
    }))]
    #[case::bare("*", AromaticityFailurePolicy::Error, Solution::Underdetermined(ResolveState { completions: ring_completions(), ..Default::default() }))]
    #[case::bare_keep("*", AromaticityFailurePolicy::Keep, Solution::Underdetermined(ResolveState { completions: ring_completions(), ..Default::default() }))]
    #[case::open("*#a+", AromaticityFailurePolicy::Error, Solution::Underdetermined(ResolveState { completions: ring_completions(), ..Default::default() }))]
    #[case::open_keep("*#a+", AromaticityFailurePolicy::Keep, Solution::Underdetermined(ResolveState { completions: ring_completions(), ..Default::default() }))]
    fn test_aromaticity_resolver_select_wildcard(
        mut benzene: Molecule,
        #[from(ring_completions)] candidates: AtomCompletions,
        #[case] member: &str,
        #[case] failure: AromaticityFailurePolicy,
        #[case] expected: SelectOutcome,
        #[values(
            AtomSelectionAlgorithm::FullDfs,
            AtomSelectionAlgorithm::BacktrackingDfs
        )]
        atom_selection_algorithm: AtomSelectionAlgorithm,
    ) {
        *benzene.atom_mut(AtomId(0)).attributes_mut() = atom_dsl!(member);
        *benzene.atom_mut(AtomId(1)).attributes_mut() = atom_dsl!("C#a+");
        let state = ResolveState {
            completions: candidates,
            ..Default::default()
        };
        let resolver = AromaticityResolver::with_config(
            &AromaticityModel::permissive(),
            AromaticityResolveConfig {
                aromatic_valence_failure: failure,
                atom_selection_algorithm,
                ..Default::default()
            },
        );
        assert_eq!(
            resolver.select(&benzene, state.clone(), ValenceTieBreak::Strict),
            Ok(expected.clone())
        );
        assert_eq!(
            resolver.select_with_report(&benzene, state, ValenceTieBreak::Strict),
            Ok(expected)
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
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}]
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
    #[case::open_completion_with_member_evidence(
        benzene(),
        AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("C#a+")])]),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState {
            completions: AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("C#a+")])]),
            aromatic_systems: vec![((0..6).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![1; 6]).with_charge(0).with_unpaired_electrons(UnpairedElectronsForm::closed_shell()))],
            tie_breaks: Vec::new(),
        })
    )]
    #[case::open_completion_without_member_evidence(
        mol_dsl!(r#"{:atoms ["C#a+" "C#a1" "C#a1" "C#a1" "C#a1" "*#a1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("C#a+")])]),
        ValenceTieBreak::Strict,
        Solution::Underdetermined(ResolveState {
            completions: AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("C#a+")])]),
            ..Default::default()
        }),
    )]
    #[case::stored_open_electrons(
        mol_dsl!(r#"{:atoms ["*#a1" "C#a1" "C#a1" "C#a1" "C#a1" "C#a1"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [2 1 0 5 4 3] :attrs "*#c0#u0#s"}]}"#),
        AtomCompletions::new(),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState {
            aromatic_systems: vec![((0..6).map(AtomId).collect(),
                AromaticSystemForm::default().with_charge(0).with_unpaired_electrons(UnpairedElectronsForm::closed_shell()))],
            ..ResolveState::default()
        })
    )]
    #[case::stored_aromatic_system_refines_assertion(
        mol_dsl!(r#"{
            :atoms ["C#a+" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}]
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
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}]
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
