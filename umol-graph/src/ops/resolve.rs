//! Composite resolver: chains the per-entity resolvers (valence,
//! aromaticity, stereo, bonds, multicenter bonds) on a single `Molecule`.
//!
//! `Determined` requires every entity (atoms, bonds, dative bonds, aromatic
//! systems, multicenter bonds, noncovalent bonds) to be ground.

pub mod aromaticity;
pub mod bonds;
pub mod multicenter;
pub mod stereo;
pub mod valence;

use std::any::Any;
use std::collections::{BTreeMap, BTreeSet};

pub use aromaticity::{
    AromaticBondConstraintMismatchPolicy, AromaticityFailurePolicy, AromaticityMismatchPolicy,
    AromaticityResolveConfig, AromaticityResolver,
};
pub use bonds::{BondsContradiction, BondsError, BondsResolver};
pub use multicenter::{
    MulticenterBondsContradiction, MulticenterBondsError, MulticenterBondsResolver,
};
pub use stereo::{
    StereoContradiction, StereoError, StereoFailurePolicy, StereoMismatchPolicy,
    StereoResolveConfig, StereoResolver,
};
use thiserror::Error;
use umol_graph_ir::ir::{
    AromaticSystemConstraintForm, AromaticSystemConstraintKey, AromaticSystemForm,
    AromaticSystemHandle, AromaticSystemId, AromaticSystemUpdate, AromaticValenceForm,
    AtomConstraintForm, AtomConstraintKey, AtomHandle, AtomId, AtomUpdate, BondConstraintForm,
    BondConstraintKey, BondHandle, BondId, BondUpdate, BooleanForm, CisTransStereoForm, Constraint,
    ConstraintEdit, DativeBondConstraintForm, DativeBondConstraintKey, DativeBondHandle,
    DativeBondId, DativeBondUpdate, Edits, Lattice, Molecule, MulticenterBondConstraintForm,
    MulticenterBondConstraintKey, MulticenterBondHandle, MulticenterBondId, MulticenterBondUpdate,
    NoncovalentBondConstraintForm, NoncovalentBondConstraintKey, NoncovalentBondHandle,
    NoncovalentBondId, NoncovalentBondUpdate, Normalize, RingModel, RingSetKind,
    StereoAtomConstraintForm, StereoAtomConstraintKey, StereoAtomHandle, StereoAtomId,
    StereoAtomUpdate, StereoBondConstraintForm, StereoBondConstraintKey, StereoBondHandle,
    StereoBondId, StereoBondUpdate, StereoKind, TetrahedralStereoForm, Transaction,
    TransactionError,
};
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;
pub use valence::{ValenceContradiction, ValenceError, ValenceResolver};

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
use crate::ops::model::{ChemistryModel, ValenceTieBreak};
use crate::ops::valence::compare::compare_by_key;
use crate::ops::valence::{AtomCompletions, ResolveReport};
use crate::ops::validate::{
    ConstraintInvariantsContradiction, ConstraintInvariantsError, ConstraintInvariantsValidator,
    ConstraintValidateConfig, DerivedKind,
};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ResolveConfig {
    pub aromaticity: AromaticityResolveConfig,
    pub stereo: StereoResolveConfig,
}

/// Solver state threaded through the constitution round: the per-atom
/// candidate sets, the accepted aromatic systems pending materialization, and
/// the atoms selected by the tie-break key.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ResolveState {
    pub completions: AtomCompletions,
    pub systems: Vec<(Vec<AtomId>, AromaticSystemForm)>,
    pub tie_breaks: Vec<AtomId>,
}

impl ResolveState {
    /// The report projection: the plural survivors and the recorded
    /// tie-break uses.
    pub fn to_report(&self) -> ResolveReport {
        let mut unresolved = AtomCompletions::new();
        for (atom, disjuncts) in self.completions.iter() {
            if disjuncts.len() > 1 {
                unresolved.insert(atom, disjuncts.iter().cloned().collect());
            }
        }
        ResolveReport {
            unresolved,
            tie_breaks: self.tie_breaks.clone(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Resolver<'a> {
    pub valence: ValenceResolver<'a>,
    pub aromaticity: AromaticityResolver,
    pub stereo: StereoResolver,
    pub bonds: BondsResolver,
    pub multicenter_bonds: MulticenterBondsResolver,
    pub tie_break: ValenceTieBreak,
    pub config: ResolveConfig,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ResolveContradiction {
    #[error(transparent)]
    Placement(#[from] PlacementContradiction),
    #[error(transparent)]
    Valence(#[from] ValenceContradiction),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityContradiction),
    #[error(transparent)]
    Stereo(#[from] StereoContradiction),
    #[error(transparent)]
    Bonds(#[from] BondsContradiction),
    #[error(transparent)]
    MulticenterBonds(#[from] MulticenterBondsContradiction),
    #[error(transparent)]
    Discharge(#[from] DischargeContradiction),
}

/// Opening placement stage contradiction: an unsatisfiable molecule-scope
/// assertion, or colliding assertions whose meet is `⊥`.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum PlacementContradiction {
    #[error("placement: molecule constraint is unsatisfiable: {constraint:?}")]
    Normalize { constraint: Constraint },
    #[error("placement: colliding assertions meet to bottom: {constraint:?}")]
    Collision { constraint: Constraint },
}

/// Closing discharge pass contradiction: a stored assertion irreconcilable
/// with its derived value, or a molecule-scope constraint decided false.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum DischargeContradiction {
    #[error("discharge: assertion contradicts the derived value: {constraint:?}")]
    Assertion { constraint: Constraint },
    #[error("discharge: {0}")]
    Molecule(#[from] ConstraintInvariantsContradiction),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ResolveError {
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
    #[error(transparent)]
    Stereo(#[from] StereoError),
    #[error(transparent)]
    Bonds(#[from] BondsError),
    #[error(transparent)]
    MulticenterBonds(#[from] MulticenterBondsError),
    #[error("constitution commit failed: {0}")]
    Commit(TransactionError),
    #[error("placement commit failed: {0}")]
    Placement(TransactionError),
    #[error("discharge evaluation failed: {0}")]
    DischargeEvaluation(#[from] ConstraintInvariantsError),
    #[error("discharge commit failed: {0}")]
    Discharge(TransactionError),
    #[error("rollback failed after {cause}: {rollback}")]
    RollbackFailed {
        cause: ResolveRollbackCause,
        rollback: TransactionError,
    },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ResolveRollbackCause {
    #[error("resolve contradiction: {0}")]
    Contradiction(ResolveContradiction),
    #[error("resolve underdetermined")]
    Underdetermined,
    #[error("resolve error: {0}")]
    Error(Box<ResolveError>),
}

impl UmolError for ResolveError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl UmolError for ResolveContradiction {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// Resolution left the molecule underdetermined (no contradiction, but not
/// concrete). Surfaced as an error only at boundaries that require a
/// determined result; carries the report for inspection.
#[derive(Debug, Clone, Default, PartialEq, Eq, Error)]
#[error("resolution underdetermined")]
pub struct ResolveUnderdetermined {
    pub report: ResolveReport,
}

impl UmolError for ResolveUnderdetermined {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl<'a> Resolver<'a> {
    pub fn new(model: &'a ChemistryModel) -> Self {
        Self::with_config(model, ResolveConfig::default())
    }

    pub fn with_config(model: &'a ChemistryModel, config: ResolveConfig) -> Self {
        Self {
            valence: ValenceResolver::new(&model.valence),
            aromaticity: AromaticityResolver::with_config(&model.aromaticity, config.aromaticity),
            stereo: StereoResolver::with_config(&model.stereo, config.stereo),
            bonds: BondsResolver::new(),
            multicenter_bonds: MulticenterBondsResolver::new(),
            tie_break: model.valence.tie_break,
            config,
        }
    }

    pub fn resolve(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<ResolveReport, ResolveContradiction>, ResolveError> {
        // Opening placement stage: normalize the molecule-scope list and
        // inline bare entity leaves, collisions combining by meet.
        let placement = match plan_placement(molecule) {
            Ok(edits) => edits,
            Err(contradiction) => {
                return Ok(Solution::Contradictory(contradiction.into()));
            }
        };
        let mut editor = molecule.edit();
        let mut journal = Transaction::default();
        match editor.transact(placement) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => return Err(ResolveError::Placement(transaction)),
        }
        let placed = editor.build();
        let mut editor = placed.edit();

        let state = match self.valence.admit(&placed).map_err(ResolveError::Valence)? {
            Solution::Determined(state) => state,
            Solution::Underdetermined(_) => {
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Underdetermined,
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Ok(Solution::Underdetermined(ResolveReport::default()));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::from(contradiction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Contradiction(contradiction),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let outcome = match self.aromaticity.select(&placed, state, self.tie_break) {
            Ok(outcome) => outcome,
            Err(error) => {
                let error = ResolveError::Aromaticity(error);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Err(error);
            }
        };
        let mut state = match outcome {
            Solution::Determined(state) => state,
            Solution::Underdetermined(state) => {
                let report = state.to_report();
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Underdetermined,
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Ok(Solution::Underdetermined(report));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::from(contradiction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Contradiction(contradiction),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Ok(Solution::Contradictory(contradiction));
            }
        };

        // Finalization: the tie-break on plural atoms outside any candidate
        // system; a tie surviving the key stays plural.
        let key = self.tie_break.key();
        if !key.is_empty() {
            let plural: Vec<AtomId> = state
                .completions
                .iter()
                .filter_map(|(atom, disjuncts)| (disjuncts.len() > 1).then_some(atom))
                .collect();
            for atom in plural {
                let disjuncts = state.completions.get(atom).expect("plural atom").to_vec();
                let best = disjuncts
                    .iter()
                    .max_by(|a, b| compare_by_key(key, a, b))
                    .expect("non-empty entry")
                    .clone();
                let unique = disjuncts
                    .iter()
                    .filter(|form| compare_by_key(key, form, &best).is_eq())
                    .count()
                    == 1;
                if unique {
                    state.completions.insert(atom, smallvec::smallvec![best]);
                    state.tie_breaks.push(atom);
                }
            }
            state.tie_breaks.sort_unstable();
            state.tie_breaks.dedup();
        }

        let report = state.to_report();
        if !report.unresolved.is_empty() {
            if let Err(rollback) = journal.rollback(&mut editor) {
                return Err(ResolveError::RollbackFailed {
                    cause: ResolveRollbackCause::Underdetermined,
                    rollback,
                });
            }
            *molecule = editor.build();
            return Ok(Solution::Underdetermined(report));
        }

        // The single commit of the constitution round.
        let mut edits = Edits::new();
        for (atom, disjuncts) in state.completions.iter() {
            let current = placed.atom(atom).attributes;
            // The constraint channel holds assertions only: the commit
            // narrows fields; candidate constraints stay solver state.
            let mut selected = disjuncts[0].clone();
            selected.constraints = current.constraints.clone();
            let update = current.difference_to(&selected);
            edits.update_atom(AtomHandle::Id(atom), current, &update);
        }
        let existing: BTreeSet<Vec<AtomId>> = placed
            .aromatic_systems()
            .iter()
            .map(|system| {
                let mut atoms: Vec<AtomId> = system.atom_ids().collect();
                atoms.sort_unstable();
                atoms
            })
            .collect();
        for (atoms, system) in &state.systems {
            let mut key = atoms.clone();
            key.sort_unstable();
            if existing.contains(&key) {
                continue;
            }
            for edit in self
                .aromaticity
                .plan_system(&placed, atoms.clone(), system.clone())
            {
                edits.push(edit);
            }
        }
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error = ResolveError::Commit(transaction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Err(error);
            }
        }
        let working = editor.build();
        let mut editor = working.edit();

        let outcome = match self.stereo.plan(&working) {
            Ok(outcome) => outcome,
            Err(error) => {
                let error = ResolveError::Stereo(error);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Err(error);
            }
        };
        let edits = match outcome {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Underdetermined,
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Ok(Solution::Underdetermined(ResolveReport {
                    unresolved: AtomCompletions::new(),
                    tie_breaks: state.tie_breaks.clone(),
                }));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::Stereo(contradiction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Contradiction(contradiction),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error = ResolveError::Stereo(StereoError::Transaction(transaction));
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Err(error);
            }
        }
        let working = editor.build();
        let mut editor = working.edit();

        let edits = self.bonds.plan(&working);
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error = ResolveError::Bonds(BondsError::Transaction(transaction));
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Err(error);
            }
        }
        let working = editor.build();
        let mut editor = working.edit();

        let edits = match self.multicenter_bonds.plan(&working) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Underdetermined,
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Ok(Solution::Underdetermined(ResolveReport {
                    unresolved: AtomCompletions::new(),
                    tie_breaks: state.tie_breaks.clone(),
                }));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::MulticenterBonds(contradiction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Contradiction(contradiction),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error =
                    ResolveError::MulticenterBonds(MulticenterBondsError::Transaction(transaction));
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Err(error);
            }
        }

        // Closing discharge pass: remove determined-redundant assertions,
        // evaluate the remaining molecule-scope list.
        let working = editor.build();
        let mut editor = working.edit();
        let outcome = match self.plan_discharge(&working) {
            Ok(outcome) => outcome,
            Err(error) => {
                let error = ResolveError::DischargeEvaluation(error);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Err(error);
            }
        };
        let edits = match outcome {
            Ok(edits) => edits,
            Err(contradiction) => {
                let contradiction = ResolveContradiction::Discharge(contradiction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Contradiction(contradiction),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error = ResolveError::Discharge(transaction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolveError::RollbackFailed {
                        cause: ResolveRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *molecule = editor.build();
                return Err(error);
            }
        }

        let resolved = editor.build();
        if resolved.is_concrete() {
            *molecule = resolved;
            Ok(Solution::Determined(ResolveReport {
                unresolved: AtomCompletions::new(),
                tie_breaks: state.tie_breaks,
            }))
        } else {
            let mut editor = resolved.edit();
            if let Err(rollback) = journal.rollback(&mut editor) {
                return Err(ResolveError::RollbackFailed {
                    cause: ResolveRollbackCause::Underdetermined,
                    rollback,
                });
            }
            *molecule = editor.build();
            Ok(Solution::Underdetermined(ResolveReport {
                unresolved: AtomCompletions::new(),
                tie_breaks: state.tie_breaks,
            }))
        }
    }

    /// Plan the closing discharge pass: a stored assertion whose ground
    /// derived value refines it is redundant and removed; a meet to `⊥`
    /// contradicts unless the key's resolve failure policy tolerates the
    /// unrealized assertion; the remaining molecule-scope list is evaluated
    /// with the validator machinery (decided-true removed, decided-false
    /// contradictory, undecided kept).
    fn plan_discharge(
        &self,
        molecule: &Molecule,
    ) -> Result<Result<Edits, DischargeContradiction>, ConstraintInvariantsError> {
        let mut edits = Edits::new();

        let needs_rings = molecule.atoms().iter().any(|atom| {
            atom.attributes.constraints.iter().any(|c| {
                matches!(
                    c,
                    AtomConstraintForm::RingDegree(_)
                        | AtomConstraintForm::RingValence(_)
                        | AtomConstraintForm::RingMembership(_)
                )
            })
        }) || molecule.bonds().iter().any(|bond| {
            bond.attributes
                .constraints
                .iter()
                .any(|c| matches!(c, BondConstraintForm::RingMembership(_)))
        });
        let rings = needs_rings.then(|| {
            molecule
                .rings(
                    RingModel {
                        kind: RingSetKind::Relevant,
                        max_ring_size: 22,
                    },
                    self.config.aromaticity.perception.ring_config,
                )
                .into_ring_set()
        });

        for id in molecule.atoms().ids() {
            for asserted in molecule.atom(id).attributes.constraints.iter() {
                let mut reading = molecule.atom(id).constraints();
                if let Some(rings) = rings.as_ref() {
                    reading = reading.with_rings(rings);
                }
                let Some(derived) = reading.derived_complete(asserted.key()) else {
                    continue;
                };
                match asserted.meet(&derived) {
                    None => {
                        // The closure's negative reading marks the failure
                        // family (nothing was realized); a positive derived
                        // value marks the mismatch family (the realized
                        // entity disagrees).
                        let tolerated = match &derived {
                            AtomConstraintForm::AromaticValence(
                                AromaticValenceForm::NotAromatic,
                            ) => {
                                self.config.aromaticity.aromatic_valence_failure
                                    != AromaticityFailurePolicy::Error
                            }
                            AtomConstraintForm::AromaticValence(_) => {
                                self.config.aromaticity.aromatic_valence_mismatch
                                    != AromaticityMismatchPolicy::Error
                            }
                            AtomConstraintForm::TetrahedralStereo(
                                TetrahedralStereoForm::NotStereo,
                            ) => {
                                self.config.stereo.tetrahedral_stereo_failure
                                    != StereoFailurePolicy::Error
                            }
                            AtomConstraintForm::TetrahedralStereo(_) => {
                                self.config.stereo.tetrahedral_stereo_mismatch
                                    != StereoMismatchPolicy::Error
                            }
                            _ => false,
                        };
                        if !tolerated {
                            return Ok(Err(DischargeContradiction::Assertion {
                                constraint: Constraint::Atom(id, asserted.clone()),
                            }));
                        }
                    }
                    Some(_) if derived.is_ground() => {
                        let mut update = AtomUpdate::default();
                        update.constraints.set(asserted.as_undetermined());
                        edits.update_atom(
                            AtomHandle::Id(id),
                            molecule.atom(id).attributes,
                            &update,
                        );
                    }
                    Some(_) => {}
                }
            }
        }
        for id in molecule.bonds().ids() {
            for asserted in molecule.bond(id).attributes.constraints.iter() {
                let mut reading = molecule.bond(id).constraints();
                if let Some(rings) = rings.as_ref() {
                    reading = reading.with_rings(rings);
                }
                let Some(derived) = reading.derived_complete(asserted.key()) else {
                    continue;
                };
                match asserted.meet(&derived) {
                    None => {
                        let tolerated = match &derived {
                            BondConstraintForm::Aromatic(BooleanForm::Lit(false)) => {
                                self.config.aromaticity.aromatic_valence_failure
                                    != AromaticityFailurePolicy::Error
                            }
                            BondConstraintForm::Aromatic(_) => {
                                self.config.aromaticity.aromatic_bond_constraint_mismatch
                                    != AromaticBondConstraintMismatchPolicy::Error
                            }
                            BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo) => {
                                self.config.stereo.cis_trans_stereo_failure
                                    != StereoFailurePolicy::Error
                            }
                            BondConstraintForm::CisTransStereo(_) => {
                                self.config.stereo.cis_trans_stereo_mismatch
                                    != StereoMismatchPolicy::Error
                            }
                            _ => false,
                        };
                        if !tolerated {
                            return Ok(Err(DischargeContradiction::Assertion {
                                constraint: Constraint::Bond(id, asserted.clone()),
                            }));
                        }
                    }
                    Some(_) if derived.is_ground() => {
                        let mut update = BondUpdate::default();
                        update.constraints.set(asserted.as_undetermined());
                        edits.update_bond(
                            BondHandle::Id(id),
                            molecule.bond(id).attributes,
                            &update,
                        );
                    }
                    Some(_) => {}
                }
            }
        }
        for id in molecule.dative_bonds().ids() {
            for asserted in molecule.dative_bond(id).attributes.constraints.iter() {
                let Some(derived) = molecule
                    .dative_bond(id)
                    .constraints()
                    .derived_complete(asserted.key())
                else {
                    continue;
                };
                match asserted.meet(&derived) {
                    None => {
                        let tolerated = match &derived {
                            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)) => {
                                self.config.aromaticity.aromatic_valence_failure
                                    != AromaticityFailurePolicy::Error
                            }
                            DativeBondConstraintForm::Aromatic(_) => {
                                self.config.aromaticity.aromatic_bond_constraint_mismatch
                                    != AromaticBondConstraintMismatchPolicy::Error
                            }
                            _ => false,
                        };
                        if !tolerated {
                            return Ok(Err(DischargeContradiction::Assertion {
                                constraint: Constraint::DativeBond(id, asserted.clone()),
                            }));
                        }
                    }
                    Some(_) if derived.is_ground() => {
                        let mut update = DativeBondUpdate::default();
                        update.constraints.set(asserted.as_undetermined());
                        edits.update_dative_bond(
                            DativeBondHandle::Id(id),
                            molecule.dative_bond(id).attributes,
                            &update,
                        );
                    }
                    Some(_) => {}
                }
            }
        }
        for id in molecule.aromatic_systems().ids() {
            for asserted in molecule.aromatic_system(id).attributes.constraints.iter() {
                let Some(derived) = molecule
                    .aromatic_system(id)
                    .constraints()
                    .derived_complete(asserted.key())
                else {
                    continue;
                };
                match asserted.meet(&derived) {
                    None => {
                        return Ok(Err(DischargeContradiction::Assertion {
                            constraint: Constraint::AromaticSystem(id, asserted.clone()),
                        }));
                    }
                    Some(_) if derived.is_ground() => {
                        let mut update = AromaticSystemUpdate::default();
                        update.constraints.set(asserted.as_undetermined());
                        edits.update_aromatic_system(
                            AromaticSystemHandle::Id(id),
                            molecule.aromatic_system(id).attributes,
                            &update,
                        );
                    }
                    Some(_) => {}
                }
            }
        }
        for id in molecule.multicenter_bonds().ids() {
            for asserted in molecule.multicenter_bond(id).attributes.constraints.iter() {
                let Some(derived) = molecule
                    .multicenter_bond(id)
                    .constraints()
                    .derived_complete(asserted.key())
                else {
                    continue;
                };
                match asserted.meet(&derived) {
                    None => {
                        return Ok(Err(DischargeContradiction::Assertion {
                            constraint: Constraint::MulticenterBond(id, asserted.clone()),
                        }));
                    }
                    Some(_) if derived.is_ground() => {
                        let mut update = MulticenterBondUpdate::default();
                        update.constraints.set(asserted.as_undetermined());
                        edits.update_multicenter_bond(
                            MulticenterBondHandle::Id(id),
                            molecule.multicenter_bond(id).attributes,
                            &update,
                        );
                    }
                    Some(_) => {}
                }
            }
        }
        for id in molecule.noncovalent_bonds().ids() {
            for asserted in molecule.noncovalent_bond(id).attributes.constraints.iter() {
                let Some(derived) = molecule
                    .noncovalent_bond(id)
                    .constraints()
                    .derived_complete(asserted.key())
                else {
                    continue;
                };
                match asserted.meet(&derived) {
                    None => {
                        return Ok(Err(DischargeContradiction::Assertion {
                            constraint: Constraint::NoncovalentBond(id, asserted.clone()),
                        }));
                    }
                    Some(_) if derived.is_ground() => {
                        let mut update = NoncovalentBondUpdate::default();
                        update.constraints.set(asserted.as_undetermined());
                        edits.update_noncovalent_bond(
                            NoncovalentBondHandle::Id(id),
                            molecule.noncovalent_bond(id).attributes,
                            &update,
                        );
                    }
                    Some(_) => {}
                }
            }
        }
        // Stereo entity constraint keys derive vacuous, so their assertions
        // are always kept; the loops exist for the uniform surface.
        for id in molecule.stereo_atoms().ids() {
            for asserted in molecule.stereo_atom(id).attributes.constraints.iter() {
                if molecule
                    .stereo_atom(id)
                    .constraints()
                    .derived_complete(asserted.key())
                    .is_some()
                {
                    unreachable!("stereo atom constraint keys have no projection");
                }
            }
        }
        for id in molecule.stereo_bonds().ids() {
            for asserted in molecule.stereo_bond(id).attributes.constraints.iter() {
                if molecule
                    .stereo_bond(id)
                    .constraints()
                    .derived_complete(asserted.key())
                    .is_some()
                {
                    unreachable!("stereo bond constraint keys have no projection");
                }
            }
        }

        let validator = ConstraintInvariantsValidator::new(ConstraintValidateConfig {
            relevant_cycle_algorithm: self
                .config
                .aromaticity
                .perception
                .ring_config
                .relevant_cycle_algorithm,
            connected_components_algorithm: self
                .config
                .aromaticity
                .perception
                .connected_components_algorithm,
            derived_kind: DerivedKind::DerivedComplete,
        });
        for constraint in molecule.constraints().iter() {
            match validator.evaluate(molecule, constraint)? {
                Solution::Determined(()) => {
                    edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
                }
                Solution::Underdetermined(()) => {}
                Solution::Contradictory(contradiction) => {
                    return Ok(Err(DischargeContradiction::Molecule(contradiction)));
                }
            }
        }

        Ok(Ok(edits))
    }
}

/// Plan the opening placement stage: normalize every molecule-scope
/// constraint (trivial wrappers reduce to their element) and move bare
/// entity leaves into the targeted entity's store, collisions combining by
/// meet.
fn plan_placement(molecule: &Molecule) -> Result<Edits, PlacementContradiction> {
    let mut edits = Edits::new();
    let mut atoms: BTreeMap<(AtomId, AtomConstraintKey), AtomConstraintForm> = BTreeMap::new();
    let mut bonds: BTreeMap<(BondId, BondConstraintKey), BondConstraintForm> = BTreeMap::new();
    let mut dative: BTreeMap<(DativeBondId, DativeBondConstraintKey), DativeBondConstraintForm> =
        BTreeMap::new();
    let mut aromatic: BTreeMap<
        (AromaticSystemId, AromaticSystemConstraintKey),
        AromaticSystemConstraintForm,
    > = BTreeMap::new();
    let mut multicenter: BTreeMap<
        (MulticenterBondId, MulticenterBondConstraintKey),
        MulticenterBondConstraintForm,
    > = BTreeMap::new();
    let mut noncovalent: BTreeMap<
        (NoncovalentBondId, NoncovalentBondConstraintKey),
        NoncovalentBondConstraintForm,
    > = BTreeMap::new();
    let mut stereo_atoms: BTreeMap<
        (StereoAtomId, StereoKind, StereoAtomConstraintKey),
        StereoAtomConstraintForm,
    > = BTreeMap::new();
    let mut stereo_bonds: BTreeMap<
        (StereoBondId, StereoKind, StereoBondConstraintKey),
        StereoBondConstraintForm,
    > = BTreeMap::new();

    for constraint in molecule.constraints().iter() {
        let normalized =
            constraint
                .clone()
                .normalize()
                .map_err(|_| PlacementContradiction::Normalize {
                    constraint: constraint.clone(),
                })?;
        let collision = |constraint: Constraint| PlacementContradiction::Collision { constraint };
        match normalized {
            Constraint::Atom(id, inner) => {
                let stored = atoms.remove(&(id, inner.key())).or_else(|| {
                    molecule
                        .atom(id)
                        .attributes
                        .constraints
                        .get(inner.key())
                        .cloned()
                });
                let met = match stored {
                    Some(stored) => stored
                        .meet(&inner)
                        .ok_or_else(|| collision(Constraint::Atom(id, inner.clone())))?,
                    None => inner,
                };
                atoms.insert((id, met.key()), met);
                edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
            }
            Constraint::Bond(id, inner) => {
                let stored = bonds.remove(&(id, inner.key())).or_else(|| {
                    molecule
                        .bond(id)
                        .attributes
                        .constraints
                        .get(inner.key())
                        .cloned()
                });
                let met = match stored {
                    Some(stored) => stored
                        .meet(&inner)
                        .ok_or_else(|| collision(Constraint::Bond(id, inner.clone())))?,
                    None => inner,
                };
                bonds.insert((id, met.key()), met);
                edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
            }
            Constraint::DativeBond(id, inner) => {
                let stored = dative.remove(&(id, inner.key())).or_else(|| {
                    molecule
                        .dative_bond(id)
                        .attributes
                        .constraints
                        .get(inner.key())
                        .cloned()
                });
                let met = match stored {
                    Some(stored) => stored
                        .meet(&inner)
                        .ok_or_else(|| collision(Constraint::DativeBond(id, inner.clone())))?,
                    None => inner,
                };
                dative.insert((id, met.key()), met);
                edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
            }
            Constraint::AromaticSystem(id, inner) => {
                let stored = aromatic.remove(&(id, inner.key())).or_else(|| {
                    molecule
                        .aromatic_system(id)
                        .attributes
                        .constraints
                        .get(inner.key())
                        .cloned()
                });
                let met = match stored {
                    Some(stored) => stored
                        .meet(&inner)
                        .ok_or_else(|| collision(Constraint::AromaticSystem(id, inner.clone())))?,
                    None => inner,
                };
                aromatic.insert((id, met.key()), met);
                edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
            }
            Constraint::MulticenterBond(id, inner) => {
                let stored = multicenter.remove(&(id, inner.key())).or_else(|| {
                    molecule
                        .multicenter_bond(id)
                        .attributes
                        .constraints
                        .get(inner.key())
                        .cloned()
                });
                let met = match stored {
                    Some(stored) => stored
                        .meet(&inner)
                        .ok_or_else(|| collision(Constraint::MulticenterBond(id, inner.clone())))?,
                    None => inner,
                };
                multicenter.insert((id, met.key()), met);
                edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
            }
            Constraint::NoncovalentBond(id, inner) => {
                let stored = noncovalent.remove(&(id, inner.key())).or_else(|| {
                    molecule
                        .noncovalent_bond(id)
                        .attributes
                        .constraints
                        .get(inner.key())
                        .cloned()
                });
                let met = match stored {
                    Some(stored) => stored
                        .meet(&inner)
                        .ok_or_else(|| collision(Constraint::NoncovalentBond(id, inner.clone())))?,
                    None => inner,
                };
                noncovalent.insert((id, met.key()), met);
                edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
            }
            Constraint::StereoAtom(id, kind, inner) => {
                let stored = stereo_atoms.remove(&(id, kind, inner.key())).or_else(|| {
                    molecule
                        .stereo_atom(id)
                        .attributes
                        .constraints
                        .get(inner.key())
                        .cloned()
                });
                let met = match stored {
                    Some(stored) => stored.meet(&inner).ok_or_else(|| {
                        collision(Constraint::StereoAtom(id, kind, inner.clone()))
                    })?,
                    None => inner,
                };
                stereo_atoms.insert((id, kind, met.key()), met);
                edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
            }
            Constraint::StereoBond(id, kind, inner) => {
                let stored = stereo_bonds.remove(&(id, kind, inner.key())).or_else(|| {
                    molecule
                        .stereo_bond(id)
                        .attributes
                        .constraints
                        .get(inner.key())
                        .cloned()
                });
                let met = match stored {
                    Some(stored) => stored.meet(&inner).ok_or_else(|| {
                        collision(Constraint::StereoBond(id, kind, inner.clone()))
                    })?,
                    None => inner,
                };
                stereo_bonds.insert((id, kind, met.key()), met);
                edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
            }
            normalized => {
                if normalized != *constraint {
                    edits.remove_molecule_constraint(ConstraintEdit::from(constraint.clone()));
                    edits.add_molecule_constraint(ConstraintEdit::from(normalized));
                }
            }
        }
    }

    for ((id, _), form) in atoms {
        let mut update = AtomUpdate::default();
        update.constraints.set(form);
        edits.update_atom(AtomHandle::Id(id), molecule.atom(id).attributes, &update);
    }
    for ((id, _), form) in bonds {
        let mut update = BondUpdate::default();
        update.constraints.set(form);
        edits.update_bond(BondHandle::Id(id), molecule.bond(id).attributes, &update);
    }
    for ((id, _), form) in dative {
        let mut update = DativeBondUpdate::default();
        update.constraints.set(form);
        edits.update_dative_bond(
            DativeBondHandle::Id(id),
            molecule.dative_bond(id).attributes,
            &update,
        );
    }
    for ((id, _), form) in aromatic {
        let mut update = AromaticSystemUpdate::default();
        update.constraints.set(form);
        edits.update_aromatic_system(
            AromaticSystemHandle::Id(id),
            molecule.aromatic_system(id).attributes,
            &update,
        );
    }
    for ((id, _), form) in multicenter {
        let mut update = MulticenterBondUpdate::default();
        update.constraints.set(form);
        edits.update_multicenter_bond(
            MulticenterBondHandle::Id(id),
            molecule.multicenter_bond(id).attributes,
            &update,
        );
    }
    for ((id, _), form) in noncovalent {
        let mut update = NoncovalentBondUpdate::default();
        update.constraints.set(form);
        edits.update_noncovalent_bond(
            NoncovalentBondHandle::Id(id),
            molecule.noncovalent_bond(id).attributes,
            &update,
        );
    }
    for ((id, _, _), form) in stereo_atoms {
        let mut update = StereoAtomUpdate::default();
        update.constraints.set(form);
        edits.update_stereo_atom(
            StereoAtomHandle::Id(id),
            molecule.stereo_atom(id).attributes,
            &update,
        );
    }
    for ((id, _, _), form) in stereo_bonds {
        let mut update = StereoBondUpdate::default();
        update.constraints.set(form);
        edits.update_stereo_bond(
            StereoBondHandle::Id(id),
            molecule.stereo_bond(id).attributes,
            &update,
        );
    }
    Ok(edits)
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::{fixture, rstest};
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{
        AtomConstraintForm, AtomId, EntityKind, MoleculeConstraint, MulticenterValenceForm, NumForm,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_concrete};

    use super::*;
    use crate::ops::aromaticity::{AromaticityError, AromaticityInconsistency};
    use crate::ops::model::{
        AromaticityModel, AromaticityRule, AromaticityTieBreak, ChemistryModel, ElementScope,
        RingLimits, StereoModel, ValenceModel,
    };
    use crate::ops::stereo::StereoInconsistency;
    use crate::ops::valence::{AtomTypeRegistry, ValenceTable};
    use crate::ops::validate::{ConnectivityModel, IncidenceConstraintInvariantsContradiction};

    #[fixture]
    fn chemistry_model() -> ChemistryModel {
        ChemistryModel {
            connectivity: ConnectivityModel::default(),
            valence: ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
            aromaticity: AromaticityModel {
                scope: ElementScope::AllowList(vec![Element::C]),
                rule: AromaticityRule::Hueckel {
                    ring_limits: RingLimits::default(),
                },
                tie_break: AromaticityTieBreak::Strict,
            },
            stereo: StereoModel::default(),
        }
    }

    #[fixture]
    fn aromatic_molecule() -> Molecule {
        mol_dsl_concrete!(
            r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#
        )
    }

    #[fixture]
    fn stereo_molecule() -> Molecule {
        mol_dsl_concrete!(
            r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
                :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#
        )
    }

    #[rstest]
    fn test_resolve_config_default() {
        assert_eq!(
            ResolveConfig::default(),
            ResolveConfig {
                aromaticity: AromaticityResolveConfig {
                    reset_aromatic_valence: false,
                    ..AromaticityResolveConfig::default()
                },
                stereo: StereoResolveConfig {
                    reset_stereo_constraints: false,
                    ..StereoResolveConfig::default()
                },
            }
        );
    }

    #[rstest]
    #[case::contradiction(
        ResolveError::RollbackFailed {
            cause: ResolveRollbackCause::Contradiction(
                ResolveContradiction::Stereo(StereoContradiction::Inconsistency(
                    StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) },
                )),
            ),
            rollback: TransactionError::OldStateMismatch,
        },
        "rollback failed after resolve contradiction: stereo inconsistency: tetrahedral stereo constraint at atom AtomId(1) cannot be realized: precondition failed: old state does not match current"
    )]
    #[case::execution(
        ResolveError::RollbackFailed {
            cause: ResolveRollbackCause::Error(Box::new(ResolveError::Bonds(
                BondsError::Transaction(TransactionError::HandleOutOfRange {
                    kind: EntityKind::Bond,
                    index: 3,
                    count: 2,
                }),
            ))),
            rollback: TransactionError::OldStateMismatch,
        },
        "rollback failed after resolve error: bond handle 3 is out of range for 2 entries: precondition failed: old state does not match current"
    )]
    #[case::underdetermined(
        ResolveError::RollbackFailed {
            cause: ResolveRollbackCause::Underdetermined,
            rollback: TransactionError::OldStateMismatch,
        },
        "rollback failed after resolve underdetermined: precondition failed: old state does not match current"
    )]
    fn test_resolver_error(#[case] error: ResolveError, #[case] expected: &str) {
        assert_eq!(error.to_string(), expected);
    }

    #[rstest]
    fn test_resolver_new(
        chemistry_model: ChemistryModel,
        aromatic_molecule: Molecule,
        stereo_molecule: Molecule,
    ) {
        let resolver = Resolver::new(&chemistry_model);
        let explicit = Resolver::with_config(&chemistry_model, ResolveConfig::default());

        assert_eq!(
            resolver.aromaticity.plan(&aromatic_molecule),
            explicit.aromaticity.plan(&aromatic_molecule)
        );
        assert_eq!(
            resolver.stereo.plan(&stereo_molecule),
            explicit.stereo.plan(&stereo_molecule)
        );
    }

    #[rstest]
    #[case::reset_aromatic_valence(ResolveConfig {
        aromaticity: AromaticityResolveConfig {
            reset_aromatic_valence: true,
            ..AromaticityResolveConfig::default()
        },
        stereo: StereoResolveConfig::default(),
    })]
    #[case::reset_stereo_constraints(ResolveConfig {
        aromaticity: AromaticityResolveConfig::default(),
        stereo: StereoResolveConfig {
            reset_stereo_constraints: true,
            ..StereoResolveConfig::default()
        },
    })]
    fn test_resolver_with_config(
        chemistry_model: ChemistryModel,
        aromatic_molecule: Molecule,
        stereo_molecule: Molecule,
        #[case] config: ResolveConfig,
    ) {
        let resolver = Resolver::with_config(&chemistry_model, config);
        let expected_aromaticity =
            AromaticityResolver::with_config(&chemistry_model.aromaticity, config.aromaticity)
                .plan(&aromatic_molecule);
        let expected_stereo = StereoResolver::with_config(&chemistry_model.stereo, config.stereo)
            .plan(&stereo_molecule);

        assert_eq!(
            resolver.aromaticity.plan(&aromatic_molecule),
            expected_aromaticity
        );
        assert_eq!(resolver.stereo.plan(&stereo_molecule), expected_stereo);
        if config.aromaticity != AromaticityResolveConfig::default() {
            assert_ne!(
                expected_aromaticity,
                AromaticityResolver::new(&chemistry_model.aromaticity).plan(&aromatic_molecule)
            );
        }
        if config.stereo != StereoResolveConfig::default() {
            assert_ne!(
                expected_stereo,
                StereoResolver::new(&chemistry_model.stereo).plan(&stereo_molecule)
            );
        }
    }

    #[rstest]
    #[case::counts(ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())))]
    #[case::atom_typing(ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        )]))))]
    fn test_resolver_resolve(#[case] valence: ValenceModel) {
        let model = ChemistryModel {
            connectivity: ConnectivityModel::default(),
            valence,
            aromaticity: AromaticityModel {
                scope: ElementScope::AllowList(vec![Element::C]),
                rule: AromaticityRule::Hueckel {
                    ring_limits: RingLimits::default(),
                },
                tie_break: AromaticityTieBreak::Strict,
            },
            stereo: StereoModel::default(),
        };
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#v0#a!"]}"#);
        assert_eq!(
            Resolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#));
    }

    #[rstest]
    #[case::aromaticity(
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"]
            :bonds [[0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"]
                    [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"]
                    [4 5 "1#c0#u0#s"] [5 0 "1#c0#u0#s"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::aromaticity_bond_marks(
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"
                    "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"
                    "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"]
            :bonds [[0 1 "1#c0#u0#s#a"] [1 2 "1#c0#u0#s#a"]
                    [2 3 "1#c0#u0#s#a"] [3 4 "1#c0#u0#s#a"]
                    [4 5 "1#c0#u0#s#a"] [5 0 "1#c0#u0#s#a"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::aromaticity_both_marks(
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"]
            :bonds [[0 1 "1#c0#u0#s#a"] [1 2 "1#c0#u0#s#a"]
                    [2 3 "1#c0#u0#s#a"] [3 4 "1#c0#u0#s#a"]
                    [4 5 "1#c0#u0#s#a"] [5 0 "1#c0#u0#s#a"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::stereo(
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h*#n0#u0#s#T1" "F#i=#c0#h0#n0#u0#s"
                    "Cl#i=#c0#h0#n0#u0#s" "Br#i=#c0#h0#n0#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"]]
        }"#),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]
        }"#)
    )]
    fn test_resolver_resolve_stages(
        chemistry_model: ChemistryModel,
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
    ) {
        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::leaf(Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::Lit(0)),))]
    #[case::singleton_wrapper(Constraint::Or(vec![Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::Valence(NumForm::Lit(0)),
    )]))]
    fn test_resolver_resolve_placement(
        chemistry_model: ChemistryModel,
        #[case] constraint: Constraint,
    ) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#);
        molecule.constraints_mut().push(constraint);

        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#));
    }

    #[rstest]
    fn test_resolver_resolve_placement_collision(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#v0"]}"#);
        molecule.constraints_mut().push(Constraint::Atom(
            AtomId(0),
            AtomConstraintForm::Valence(NumForm::Lit(3)),
        ));
        let before = molecule.clone();

        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Contradictory(ResolveContradiction::Placement(
                PlacementContradiction::Collision {
                    constraint: Constraint::Atom(
                        AtomId(0),
                        AtomConstraintForm::Valence(NumForm::Lit(3)),
                    ),
                },
            )))
        );
        assert_eq!(molecule, before);
    }

    #[rstest]
    fn test_resolver_resolve_discharge_error(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#D5"]}"#);
        let before = molecule.clone();

        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Contradictory(ResolveContradiction::Discharge(
                DischargeContradiction::Assertion {
                    constraint: Constraint::Atom(
                        AtomId(0),
                        AtomConstraintForm::Degree(NumForm::Lit(5)),
                    ),
                },
            )))
        );
        assert_eq!(molecule, before);
    }

    #[rstest]
    #[case::decided_true(
        MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: NumForm::Lit(0),
        },
        Ok(Solution::Determined(ResolveReport::default()))
    )]
    fn test_resolver_resolve_discharge_molecule_scope(
        chemistry_model: ChemistryModel,
        #[case] constraint: MoleculeConstraint,
        #[case] expected: Result<Solution<ResolveReport, ResolveContradiction>, ResolveError>,
    ) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#);
        molecule
            .constraints_mut()
            .push(Constraint::Molecule(constraint));

        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            expected
        );
        assert!(molecule.constraints().is_empty());
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#));
    }

    #[rstest]
    fn test_resolver_resolve_discharge_molecule_scope_error(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#);
        molecule
            .constraints_mut()
            .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                atoms: None,
                sum: NumForm::Lit(5),
            }));
        let before = molecule.clone();

        let outcome = Resolver::new(&chemistry_model).resolve(&mut molecule);
        assert!(matches!(
            outcome,
            Ok(Solution::Contradictory(ResolveContradiction::Discharge(
                DischargeContradiction::Molecule(_)
            )))
        ));
        assert_eq!(molecule, before);
    }

    #[rstest]
    fn test_resolver_resolve_pyrrolyl() {
        let model = ChemistryModel {
            valence: ValenceModel {
                tie_break: ValenceTieBreak::MostSaturated,
                ..ValenceModel::default()
            },
            aromaticity: AromaticityModel::daylight(),
            ..ChemistryModel::default()
        };
        let mut molecule = mol_dsl!(
            r#"{:atoms ["C#i=#c0#n0#u0#s#a+" "C#i=#c0#n0#u0#s#a+"
                        "C#i=#c0#n0#u0#s#a+" "C#i=#c0#n0#u0#s#a+"
                        "N#i=#c0#h0#a+"]
                :bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"]
                        [1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"]
                        [3 4 "1#c0#u0#s#a"]]}"#
        );

        assert_eq!(
            Resolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(
            molecule,
            mol_dsl!(
                r#"{:aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,2]#c0#u0#s"}]
                    :atoms ["C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"
                            "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"
                            "N#i=#c0#h0#n0#u#s2"]
                    :bonds [[0 4 "1#c0#u0#s"] [0 1 "1#c0#u0#s"]
                            [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"]
                            [3 4 "1#c0#u0#s"]]}"#
            )
        );
    }

    #[rstest]
    fn test_resolver_resolve_underdetermined(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl!(
            r#"{
            :atoms ["C#i*#c0#h4#n0#u0#s#v0#a!" "C#i=#c0#h4#n0#u0#s"]
            :noncovalent-bonds [{:atoms [0 1] :attrs "*"}]
        }"#
        );
        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Underdetermined(ResolveReport::default()))
        );
        assert_eq!(
            molecule,
            mol_dsl!(
                r#"{
                :atoms ["C#i*#c0#h4#n0#u0#s#v0#a!" "C#i=#c0#h4#n0#u0#s"]
                :noncovalent-bonds [{:atoms [0 1] :attrs "*"}]
            }"#
            )
        );
    }

    #[rstest]
    #[case::counts(ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())))]
    #[case::atom_typing(ValenceModel::atom_typing(Cow::Borrowed(
        AtomTypeRegistry::default_registry()
    )))]
    fn test_resolver_resolve_partial(#[case] valence: ValenceModel) {
        let model = ChemistryModel {
            valence,
            ..ChemistryModel::default()
        };
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#);
        let original = molecule.clone();

        assert_eq!(
            Resolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Underdetermined(ResolveReport::default()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_resolver_resolve_later_underdetermined_rollback(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl!(
            r#"{
            :atoms ["C#i*#c0#h*#n0#u0#s#T+" "F#i=#c0#h0#n0#u0#s"
                    "Cl#i=#c0#h0#n0#u0#s" "Br#i=#c0#h0#n0#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"]]
        }"#
        );
        let original = molecule.clone();

        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Underdetermined(ResolveReport::default()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::aromaticity_setup(
        ChemistryModel {
            valence: ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
                    "C#i=#c0#h0#n0#u0#s#v2#a2"
                )]))),
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 }, tie_break: AromaticityTieBreak::Strict },
            ..ChemistryModel::default()
        },
        mol_dsl!(r#"{
            :atoms ["C#i*#v2#a2" "C#v2#a2" "C#v2#a2"
                    "C#v2#a2" "C#v2#a2" "C#v2#a2"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#),
        ResolveError::Aromaticity(AromaticityError::HmoMissingParameters(
            String::from("no Van-Catledge parameters for C with 2 pi-electrons"),
        )),
    )]
    fn test_resolver_resolve_error(
        #[case] model: ChemistryModel,
        #[case] mut molecule: Molecule,
        #[case] expected: ResolveError,
    ) {
        let original = molecule.clone();
        let resolver = Resolver::new(&model);

        assert_eq!(resolver.resolve(&mut molecule), Err(expected));
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::underdetermined(
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h0#n0#u0#s#v0#a!#m1"
                    "C#i=#c0#h0#n0#u0#s#v0#a!#m1"
                    "C#i=#c0#h0#n0#u0#s#v0#a!#m1"]
            :multicenter-bonds [{:atoms [0 1 2] :attrs "*"}]
        }"#),
        Solution::Underdetermined(ResolveReport::default())
    )]
    #[case::contradiction(
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h0#n0#u0#s#v0#a!#m1"]
        }"#),
        Solution::Contradictory(ResolveContradiction::MulticenterBonds(
            MulticenterBondsContradiction::Constraint(
                IncidenceConstraintInvariantsContradiction::Atom {
                    atom: AtomId(0),
                    constraint: AtomConstraintForm::multicenter_valence(
                        MulticenterValenceForm::multicenter(1),
                    ),
                },
            ),
        ))
    )]
    fn test_resolver_resolve_multicenter_constraint_precondition(
        #[case] mut molecule: Molecule,
        #[case] expected: Solution<ResolveReport, ResolveContradiction>,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([
                atom_dsl!("C#i=#c0#h0#n0#u0#s#v0#a!#m1"),
            ]))),
            ..ChemistryModel::default()
        };
        let original = molecule.clone();

        assert_eq!(Resolver::new(&model).resolve(&mut molecule), Ok(expected));
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::aromaticity(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar, tie_break: AromaticityTieBreak::Strict },
        mol_dsl!(r#"{
            :atoms ["N#i*#c0#h#n0#u0#s#a2" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a"]
            :bonds [[0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"]
                    [3 4 "1#c0#u0#s"] [4 0 "1#c0#u0#s"]]
        }"#),
        ResolveContradiction::Aromaticity(AromaticityContradiction::ClarNonBenzenoid(
            "Clar model requires benzenoid input but non-carbon aromatic atoms are present".to_string(),
        ))
    )]
    #[case::aromaticity_projection(
        AromaticityModel::mdl(),
        mol_dsl!(r#"{
            :atoms ["O#i*#c0#h0#n1#u0#s#a2" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a"]
            :bonds [[0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"]
                    [3 4 "1#c0#u0#s"] [4 0 "1#c0#u0#s"]]
        }"#),
        ResolveContradiction::Aromaticity(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }
        ))
    )]
    #[case::stereo(
        AromaticityModel { scope: ElementScope::AllowList(vec![Element::C]), rule: AromaticityRule::Hueckel { ring_limits: RingLimits::default() }, tie_break: AromaticityTieBreak::Strict },
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h#n0#u0#s#a#T1" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"]
            :bonds [[0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"]
                    [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"]
                    [4 5 "1#c0#u0#s"] [5 0 "1#c0#u0#s"]]
        }"#),
        ResolveContradiction::Stereo(StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(0) }
        ))
    )]
    fn test_resolver_resolve_contradiction(
        mut chemistry_model: ChemistryModel,
        #[case] aromaticity: AromaticityModel,
        #[case] mut molecule: Molecule,
        #[case] expected: ResolveContradiction,
    ) {
        chemistry_model.aromaticity = aromaticity;
        let original = molecule.clone();
        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_resolver_resolve_identity(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl_concrete!(
            r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                    [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        );
        let expected = molecule.clone();
        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, expected);
    }
}
