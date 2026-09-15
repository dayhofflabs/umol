//! Composite resolver: chains the per-entity resolvers (isotope, valence,
//! aromaticity, stereo, bonds, multicenter bonds) on a single `Molecule`.
//!
//! `Determined` requires every entity (atoms, bonds, dative bonds, aromatic
//! systems, multicenter bonds, noncovalent bonds) to be ground.

pub mod aromaticity;
pub mod bonds;
pub mod isotope;
pub mod multicenter;
pub mod stereo;
pub mod valence;

use std::any::Any;
use std::collections::{BTreeMap, BTreeSet};

use aromaticity::AromaticityProjectError;
pub use aromaticity::{
    AromaticBondConstraintMismatchPolicy, AromaticityFailurePolicy, AromaticityMismatchPolicy,
    AromaticityResolveConfig, AromaticityResolver,
};
use bitflags::bitflags;
pub use bonds::{BondsContradiction, BondsError, BondsResolver};
use isotope::IsotopeProjectError;
pub use isotope::{IsotopeContradiction, IsotopeError, IsotopePolicy, IsotopeResolver};
pub use multicenter::{
    MulticenterBondsContradiction, MulticenterBondsError, MulticenterBondsResolver,
};
use stereo::StereoProjectError;
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
    DativeBondId, DativeBondUpdate, Edits, Entity, Lattice, Molecule,
    MulticenterBondConstraintForm, MulticenterBondConstraintKey, MulticenterBondHandle,
    MulticenterBondId, MulticenterBondUpdate, NoncovalentBondConstraintForm,
    NoncovalentBondConstraintKey, NoncovalentBondHandle, NoncovalentBondId, NoncovalentBondUpdate,
    Normalize, NumForm, RingModel, RingSetKind, StereoAtomConstraintForm, StereoAtomConstraintKey,
    StereoAtomHandle, StereoAtomId, StereoAtomUpdate, StereoBondConstraintForm,
    StereoBondConstraintKey, StereoBondHandle, StereoBondId, StereoBondUpdate, StereoKind,
    TetrahedralStereoForm, TransactionError, UnpairedElectronsForm,
};
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;
use valence::ValenceProjectError;
pub use valence::{ValenceContradiction, ValenceError, ValenceResolver};

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
use crate::ops::model::{ChemistryModel, ValenceTieBreak};
use crate::ops::valence::compare::compare_by_key;
use crate::ops::valence::{AtomCompletions, ResolveReport};
use crate::ops::validate::{
    ConstraintInvariantsContradiction, ConstraintInvariantsError, ConstraintInvariantsValidator,
    ConstraintValidateConfig, DerivedKind,
};

/// Operational policies for resolution; unspecified isotope composition remains unresolved by default.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ResolveConfig {
    /// Completion of unspecified isotope composition; Strict by default.
    pub isotope: IsotopePolicy,
    pub aromaticity: AromaticityResolveConfig,
    pub stereo: StereoResolveConfig,
}

bitflags! {
    /// Stages selected for [`Resolver::project`], executed in declaration order.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub struct ProjectFlags: u8 {
        /// Replace stereo entities with fixed-frame assertions.
        const STEREO = 1 << 0;
        /// Replace aromatic systems with atom and bond assertions.
        const AROMATICITY = 1 << 1;
        /// Run valence projection, which preserves atom fields.
        const VALENCE = 1 << 2;
        /// Elide isotope defaults under the isotope policy.
        const ISOTOPE = 1 << 3;
    }
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
    pub isotope: IsotopeResolver,
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
    #[error("isotope commit failed: {0}")]
    Isotope(TransactionError),
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

/// Contradictions reported by the constituent projection phases.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ProjectContradiction {
    #[error(transparent)]
    Aromaticity(#[from] AromaticityContradiction),
    #[error(transparent)]
    Stereo(#[from] StereoContradiction),
    #[error(transparent)]
    Valence(#[from] ValenceContradiction),
}

/// Failure to project molecular information without implicit localization or loss.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ProjectError {
    #[error("{entity:?} requires concrete zero bond charge, found {charge:?}")]
    BondCharge { entity: Entity, charge: NumForm },
    #[error("{entity:?} requires concrete closed-shell singlet bond spin, found {spin:?}")]
    BondSpin {
        entity: Entity,
        spin: UnpairedElectronsForm,
    },
    #[error(transparent)]
    Stereo(#[from] StereoProjectError),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityProjectError),
    #[error(transparent)]
    Valence(#[from] ValenceProjectError),
    #[error(transparent)]
    Isotope(#[from] IsotopeProjectError),
}

impl UmolError for ProjectError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl UmolError for ProjectContradiction {
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
            isotope: IsotopeResolver::new(config.isotope),
            valence: ValenceResolver::new(&model.valence),
            aromaticity: AromaticityResolver::with_config(&model.aromaticity, config.aromaticity),
            stereo: StereoResolver::with_config(&model.stereo, config.stereo),
            bonds: BondsResolver::new(),
            multicenter_bonds: MulticenterBondsResolver::new(),
            tie_break: model.valence.tie_break,
            config,
        }
    }

    /// Resolves isotope composition before valence and aromaticity, then stereo and bonds.
    ///
    /// An unresolved isotope does not stop later phases. Only a completely determined
    /// result replaces the caller's molecule; every other outcome preserves it exactly.
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
        let editor = molecule
            .edit()
            .apply(placement)
            .map_err(ResolveError::Placement)?;
        // Placement changes assertions only, so the isotope fields still match this plan.
        let isotope_edits = match self.isotope.plan(molecule) {
            Solution::Determined(edits) | Solution::Underdetermined(edits) => edits,
            Solution::Contradictory(contradiction) => match contradiction {},
        };
        let editor = editor.apply(isotope_edits).map_err(ResolveError::Isotope)?;
        let placed = editor.build();
        let editor = placed.edit();

        let state = match self.valence.admit(&placed).map_err(ResolveError::Valence)? {
            Solution::Determined(state) => state,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(ResolveReport::default()));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::from(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let outcome = self
            .aromaticity
            .select(&placed, state, self.tie_break)
            .map_err(ResolveError::Aromaticity)?;
        let mut state = match outcome {
            Solution::Determined(state) => state,
            Solution::Underdetermined(state) => {
                let report = state.to_report();
                return Ok(Solution::Underdetermined(report));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::from(contradiction);
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
        let editor = editor.apply(edits).map_err(ResolveError::Commit)?;
        let working = editor.build();
        let editor = working.edit();

        let outcome = self.stereo.plan(&working).map_err(ResolveError::Stereo)?;
        let edits = match outcome {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(ResolveReport {
                    unresolved: AtomCompletions::new(),
                    tie_breaks: state.tie_breaks.clone(),
                }));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::Stereo(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        // These plans read independent domains of the same post-constitution snapshot. Apply
        // them in phase order without publishing between them.
        let bond_edits = self.bonds.plan(&working);
        let multicenter_outcome = self.multicenter_bonds.plan(&working);
        let editor = editor
            .apply(edits)
            .map_err(|error| ResolveError::Stereo(StereoError::Transaction(error)))?;
        let editor = editor
            .apply(bond_edits)
            .map_err(|error| ResolveError::Bonds(BondsError::Transaction(error)))?;

        let edits = match multicenter_outcome {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(ResolveReport {
                    unresolved: AtomCompletions::new(),
                    tie_breaks: state.tie_breaks.clone(),
                }));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::MulticenterBonds(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let editor = editor.apply(edits).map_err(|error| {
            ResolveError::MulticenterBonds(MulticenterBondsError::Transaction(error))
        })?;

        // Closing discharge pass: remove determined-redundant assertions,
        // evaluate the remaining molecule-scope list.
        let working = editor.build();
        let editor = working.edit();
        let outcome = self
            .plan_discharge(&working)
            .map_err(ResolveError::DischargeEvaluation)?;
        let edits = match outcome {
            Ok(edits) => edits,
            Err(contradiction) => {
                let contradiction = ResolveContradiction::Discharge(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let editor = editor.apply(edits).map_err(ResolveError::Discharge)?;

        let resolved = editor.build();
        if resolved.is_concrete() {
            *molecule = resolved;
            Ok(Solution::Determined(ResolveReport {
                unresolved: AtomCompletions::new(),
                tie_breaks: state.tie_breaks,
            }))
        } else {
            Ok(Solution::Underdetermined(ResolveReport {
                unresolved: AtomCompletions::new(),
                tie_breaks: state.tie_breaks,
            }))
        }
    }

    /// Projects selected stages within graph IR atomically.
    ///
    /// Recovers fixed-frame #T/#C assertions, writes member-aligned #a contributions and
    /// aromatic bond assertions, runs the no-op valence projection, then elides isotope
    /// defaults. Only stages selected by `flags` run, in that order. All H counts are
    /// preserved regardless of stage selection; use all flags for ordinary projection.
    /// Other fields and structures remain for the eventual format conversion. Successful
    /// projection does not establish format representability or require a concrete result:
    /// isotope default elision can leave isotope fields undetermined.
    ///
    /// # Semantic properties
    ///
    /// Only Determined publishes the candidate. Errors, contradictions, and underdetermination
    /// leave the caller's molecule unchanged. Projection does not invoke resolution or recovery
    /// comparisons. Implicit-H counts, atom electron fields, and localized charge are preserved.
    /// For inputs meeting the bond charge/spin requirements, the result agrees with executing
    /// the selected standalone projections in order. Empty flags leave such inputs unchanged.
    ///
    /// # Errors
    ///
    /// Localized and multicenter bonds require concrete zero charge and closed-shell singlet
    /// spin. Returns the owning phase's error for unprojectable stereo, aromatic-system fields,
    /// or isotopes. No charge/spin localization or unsupported-structure removal is implicit.
    pub fn project(
        &self,
        molecule: &mut Molecule,
        flags: ProjectFlags,
    ) -> Result<Solution<(), ProjectContradiction>, ProjectError> {
        let bonds = molecule
            .bonds()
            .iter()
            .map(|bond| {
                (
                    Entity::Bond(bond.id),
                    &bond.attributes.charge,
                    &bond.attributes.unpaired_electrons,
                )
            })
            .chain(molecule.multicenter_bonds().iter().map(|bond| {
                (
                    Entity::MulticenterBond(bond.id),
                    &bond.attributes.charge,
                    &bond.attributes.unpaired_electrons,
                )
            }));
        for (entity, charge, spin) in bonds {
            if !matches!(charge, NumForm::Lit(0)) {
                return Err(ProjectError::BondCharge {
                    entity,
                    charge: charge.clone(),
                });
            }
            if !matches!(
                spin,
                UnpairedElectronsForm {
                    count: NumForm::Lit(0),
                    multiplicity: NumForm::Lit(1)
                }
            ) {
                return Err(ProjectError::BondSpin {
                    entity,
                    spin: spin.clone(),
                });
            }
        }

        let mut candidate = molecule.clone();
        if flags.contains(ProjectFlags::STEREO) {
            match self.stereo.project(&mut candidate)? {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => return Ok(Solution::Underdetermined(())),
                Solution::Contradictory(contradiction) => {
                    return Ok(Solution::Contradictory(contradiction.into()))
                }
            }
        }
        if flags.contains(ProjectFlags::AROMATICITY) {
            match self.aromaticity.project(&mut candidate)? {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => return Ok(Solution::Underdetermined(())),
                Solution::Contradictory(contradiction) => {
                    return Ok(Solution::Contradictory(contradiction.into()))
                }
            }
        }
        if flags.contains(ProjectFlags::VALENCE) {
            match self.valence.project(&mut candidate, self.tie_break)? {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => return Ok(Solution::Underdetermined(())),
                Solution::Contradictory(contradiction) => {
                    return Ok(Solution::Contradictory(contradiction.into()))
                }
            }
        }
        if flags.contains(ProjectFlags::ISOTOPE) {
            match self.isotope.project(&mut candidate)? {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => return Ok(Solution::Underdetermined(())),
                Solution::Contradictory(contradiction) => match contradiction {},
            }
        }
        *molecule = candidate;
        Ok(Solution::Determined(()))
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
        AtomConstraintForm, AtomForm, AtomId, IsotopeMassForm, MoleculeConstraint, MoleculeEntries,
        MulticenterValenceForm, NumForm, StereoConfigurationForm, StereoCoset, StereoLigand,
        StereoLigandKind,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_concrete};

    use super::*;
    use crate::ingest::ingest_smiles;
    use crate::ops::aromaticity::{AromaticityError, AromaticityInconsistency};
    use crate::ops::model::{
        AromaticityModel, AromaticityRule, AromaticityTieBreak, ChemistryModel, ElementScope,
        RingLimits, StereoModel, ValenceModel,
    };
    use crate::ops::stereo::{StereoInconsistency, StereoPerception};
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
                isotope: IsotopePolicy::Strict,
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
        isotope: IsotopePolicy::Strict,
        aromaticity: AromaticityResolveConfig {
            reset_aromatic_valence: true,
            ..AromaticityResolveConfig::default()
        },
        stereo: StereoResolveConfig::default(),
    })]
    #[case::reset_stereo_constraints(ResolveConfig {
        isotope: IsotopePolicy::Strict,
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
            "C#c0#h4#n0#u0#s#v0#a!"
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
        let mut molecule = mol_dsl!(r#"{:atoms ["C#i=#c0#h4#v0#a!"]}"#);
        assert_eq!(
            Resolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#));
    }

    #[rstest]
    #[case::methane(atom_dsl!("C#i=#c0#h4#n0#u0#s"), atom_dsl!("C#i=#c0#h4"))]
    #[case::methyl(atom_dsl!("C#i=#c0#h3#n0#u1#s2"), atom_dsl!("C#i=#c0#h3"))]
    #[case::carbanion(atom_dsl!("C#i=#c-#h3#n1#u0#s"), atom_dsl!("C#i=#c-#h3"))]
    #[case::carbocation(atom_dsl!("C#i=#c+#h3#n0#u0#s"), atom_dsl!("C#i=#c+#h3"))]
    #[case::ammonia(atom_dsl!("N#i=#c0#h3#n1#u0#s"), atom_dsl!("N#i=#c0#h3"))]
    #[case::ammonium(atom_dsl!("N#i=#c+#h4#n0#u0#s"), atom_dsl!("N#i=#c+#h4"))]
    #[case::water(atom_dsl!("O#i=#c0#h2#n2#u0#s"), atom_dsl!("O#i=#c0#h2"))]
    #[case::hydroxide(atom_dsl!("O#i=#c-#h1#n3#u0#s"), atom_dsl!("O#i=#c-#h1"))]
    #[case::fluorane(atom_dsl!("F#i=#c0#h1#n3#u0#s"), atom_dsl!("F#i=#c0#h1"))]
    #[case::chloride(atom_dsl!("Cl#i=#c-#h0#n4#u0#s"), atom_dsl!("Cl#i=#c-#h0"))]
    #[case::phosphonium(atom_dsl!("P#i=#c+#h4#n0#u0#s"), atom_dsl!("P#i=#c+#h4"))]
    #[case::sulfane(atom_dsl!("S#i=#c0#h2#n2#u0#s"), atom_dsl!("S#i=#c0#h2"))]
    #[case::borane(atom_dsl!("B#i=#c0#h3#n0#u0#s"), atom_dsl!("B#i=#c0#h3"))]
    #[case::silane(atom_dsl!("Si#i=#c0#h4#n0#u0#s"), atom_dsl!("Si#i=#c0#h4"))]
    fn test_resolver_resolve_atoms(
        #[values(ValenceModel::smiles(), ValenceModel::default())] mut valence: ValenceModel,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)]
        tie_break: ValenceTieBreak,
        #[case] expected: AtomForm,
        #[case] input: AtomForm,
    ) {
        valence.tie_break = tie_break;
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![input],
            ..Default::default()
        });
        assert_eq!(
            Resolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(
            molecule,
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![expected],
                ..Default::default()
            })
        );
    }

    #[rstest]
    fn test_resolver_resolve_overlap(
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)]
        tie_break: ValenceTieBreak,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel {
                candidates: ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([
                    atom_dsl!("C#c0#h*#n0#u0#s"),
                    atom_dsl!("C#c0#h4#n0#u0#s"),
                ])))
                .candidates,
                tie_break,
            },
            ..Default::default()
        };
        let resolver = Resolver::new(&model);
        let mut molecule = mol_dsl!(r#"{:atoms ["C#i=#c0#h4"]}"#);
        assert_eq!(
            resolver.resolve(&mut molecule),
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
            Resolver::with_config(
                &chemistry_model,
                ResolveConfig {
                    isotope: IsotopePolicy::Natural,
                    ..Default::default()
                }
            )
            .resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::omitted(mol_dsl!(r#"{:atoms ["C#c0#h4"]}"#), mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#), true)]
    #[case::natural(mol_dsl!(r#"{:atoms ["C#i=#c0#h4"]}"#), mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#), false)]
    #[case::mass(mol_dsl!(r#"{:atoms ["C#i13#c0#h4"]}"#), mol_dsl!(r#"{:atoms ["C#i13#c0#h4#n0#u0#s"]}"#), false)]
    fn test_resolver_resolve_isotope(
        #[values(ValenceModel::smiles(), ValenceModel::default())] mut valence: ValenceModel,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)]
        tie_break: ValenceTieBreak,
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] isotope: IsotopePolicy,
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
        #[case] omitted: bool,
    ) {
        valence.tie_break = tie_break;
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                isotope,
                ..Default::default()
            },
        );
        let original = molecule.clone();
        let actual = resolver.resolve(&mut molecule);
        if omitted && isotope == IsotopePolicy::Strict {
            assert_eq!(
                actual,
                Ok(Solution::Underdetermined(ResolveReport::default()))
            );
            assert_eq!(molecule, original);
        } else {
            assert_eq!(actual, Ok(Solution::Determined(ResolveReport::default())));
            assert_eq!(molecule, expected);
        }
    }

    #[rstest]
    #[case::set(mol_dsl!(r#"{:atoms ["C#c0#h4" "C#i{12,13}#c0#h4"]}"#))]
    #[case::variable(mol_dsl!(r#"{:atoms ["C#c0#h4" "C#i?mass#c0#h4"]}"#))]
    fn test_resolver_resolve_isotope_partial(
        #[values(ValenceModel::smiles(), ValenceModel::default())] valence: ValenceModel,
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] isotope: IsotopePolicy,
        #[case] mut molecule: Molecule,
    ) {
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        let original = molecule.clone();
        assert_eq!(
            Resolver::with_config(
                &model,
                ResolveConfig {
                    isotope,
                    ..Default::default()
                }
            )
            .resolve(&mut molecule),
            Ok(Solution::Underdetermined(ResolveReport::default()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::undetermined(mol_dsl!(r#"{:atoms ["C#c0#h4#T1"]}"#))]
    #[case::variable(mol_dsl!(r#"{:atoms ["C#i?mass#c0#h4#T1"]}"#))]
    fn test_resolver_resolve_isotope_stereo(
        #[values(ValenceModel::smiles(), ValenceModel::default())] valence: ValenceModel,
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] isotope: IsotopePolicy,
        #[case] mut molecule: Molecule,
    ) {
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        let original = molecule.clone();
        assert_eq!(
            Resolver::with_config(
                &model,
                ResolveConfig {
                    isotope,
                    ..Default::default()
                }
            )
            .resolve(&mut molecule),
            Ok(Solution::Contradictory(ResolveContradiction::Stereo(
                StereoContradiction::Inconsistency(StereoInconsistency::TetrahedralStereoFailure {
                    atom: AtomId(0)
                })
            )))
        );
        assert_eq!(molecule, original);
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
        molecule
            .try_modify_constraints(|constraints| constraints.push(constraint))
            .expect("the test constraint references the molecule");

        assert_eq!(
            Resolver::with_config(
                &chemistry_model,
                ResolveConfig {
                    isotope: IsotopePolicy::Natural,
                    ..Default::default()
                }
            )
            .resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#));
    }

    #[rstest]
    fn test_resolver_resolve_placement_collision(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#v0"]}"#);
        molecule
            .try_modify_constraints(|constraints| {
                constraints.push(Constraint::Atom(
                    AtomId(0),
                    AtomConstraintForm::Valence(NumForm::Lit(3)),
                ));
            })
            .expect("the test constraint references the molecule");
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
            .try_modify_constraints(|constraints| {
                constraints.push(Constraint::Molecule(constraint));
            })
            .expect("the test constraint references the molecule");

        assert_eq!(
            Resolver::with_config(
                &chemistry_model,
                ResolveConfig {
                    isotope: IsotopePolicy::Natural,
                    ..Default::default()
                }
            )
            .resolve(&mut molecule),
            expected
        );
        assert!(molecule.constraints().is_empty());
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#));
    }

    #[rstest]
    fn test_resolver_resolve_discharge_molecule_scope_error(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s"]}"#);
        molecule
            .try_modify_constraints(|constraints| {
                constraints.push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
                    atoms: None,
                    sum: NumForm::Lit(5),
                }));
            })
            .expect("the test constraint references the molecule");
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
            Resolver::with_config(
                &chemistry_model,
                ResolveConfig {
                    isotope: IsotopePolicy::Natural,
                    ..Default::default()
                }
            )
            .resolve(&mut molecule),
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
    fn test_resolver_resolve_later_underdetermined(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl!(
            r#"{
            :atoms ["C#i*#c0#h*#n0#u0#s#T+" "F#i=#c0#h0#n0#u0#s"
                    "Cl#i=#c0#h0#n0#u0#s" "Br#i=#c0#h0#n0#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"]]
        }"#
        );
        let original = molecule.clone();

        assert_eq!(
            Resolver::with_config(
                &chemistry_model,
                ResolveConfig {
                    isotope: IsotopePolicy::Natural,
                    ..Default::default()
                }
            )
            .resolve(&mut molecule),
            Ok(Solution::Underdetermined(ResolveReport::default()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::aromaticity_setup(
        ChemistryModel {
            valence: ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
                    "C#c0#h0#n0#u0#s#v2#a2"
                )]))),
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 }, tie_break: AromaticityTieBreak::Strict },
            ..ChemistryModel::default()
        },
        mol_dsl!(r#"{
            :atoms ["C#i=#v2#a2" "C#i=#v2#a2" "C#i=#v2#a2"
                    "C#i=#v2#a2" "C#i=#v2#a2" "C#i=#v2#a2"]
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
            :atoms ["C#i=#c0#h0#n*#u0#s#v0#a!#m1"
                    "C#i=#c0#h0#n0#u0#s#v0#a!#m1"
                    "C#i=#c0#h0#n0#u0#s#v0#a!#m1"]
            :multicenter-bonds [{:atoms [0 1 2] :attrs "*"}]
        }"#),
        Solution::Underdetermined(ResolveReport::default())
    )]
    #[case::contradiction(
        mol_dsl!(r#"{
            :atoms ["C#i=#c0#h0#n*#u0#s#v0#a!#m1"]
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
                atom_dsl!("C#c0#h0#n0#u0#s#v0#a!#m1"),
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

    #[rstest]
    #[case::natural(
        mol_dsl_concrete!(r#"{:atoms ["C#h3#u1#s2" "C#i13#c-#h3#n1"]}"#),
        mol_dsl!(r#"{:atoms ["C#c0#h3#n0#u1#s2" "C#i13#c-#h3#n1#u0#s"]}"#))]
    fn test_resolver_project(#[case] mut molecule: Molecule, #[case] expected: Molecule) {
        let model = ChemistryModel::default();
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                isotope: IsotopePolicy::Natural,
                ..Default::default()
            },
        );
        assert_eq!(
            resolver.project(&mut molecule, ProjectFlags::all()),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::strict_natural(
        ValenceTieBreak::Strict,
        IsotopePolicy::Natural,
        r#"{:atoms ["C#c0#h4#n0#u0#s" "C#i13#c0#h4#n0#u0#s" "F#c0#h1#n3#u0#s"]}"#
    )]
    #[case::saturated_natural(
        ValenceTieBreak::MostSaturated,
        IsotopePolicy::Natural,
        r#"{:atoms ["C#c0#h4#n0#u0#s" "C#i13#c0#h4#n0#u0#s" "F#c0#h1#n3#u0#s"]}"#
    )]
    #[case::strict_isotope(
        ValenceTieBreak::Strict,
        IsotopePolicy::Strict,
        r#"{:atoms ["C#i=#c0#h4#n0#u0#s" "C#i13#c0#h4#n0#u0#s" "F#i=#c0#h1#n3#u0#s"]}"#
    )]
    #[case::saturated_isotope(
        ValenceTieBreak::MostSaturated,
        IsotopePolicy::Strict,
        r#"{:atoms ["C#i=#c0#h4#n0#u0#s" "C#i13#c0#h4#n0#u0#s" "F#i=#c0#h1#n3#u0#s"]}"#
    )]
    fn test_resolver_project_valence(
        #[case] policy: ValenceTieBreak,
        #[case] isotope: IsotopePolicy,
        #[case] expected: &str,
        #[values(false, true)] typing: bool,
    ) {
        let mut model = ChemistryModel {
            valence: if typing {
                ValenceModel::default()
            } else {
                ValenceModel::smiles()
            },
            ..Default::default()
        };
        model.valence.tie_break = policy;
        let mut molecule = mol_dsl_concrete!(r#"{:atoms ["C#h4" "C#i13#h4" "F#h1#n3"]}"#);
        assert_eq!(
            Resolver::with_config(
                &model,
                ResolveConfig {
                    isotope,
                    ..Default::default()
                }
            )
            .project(&mut molecule, ProjectFlags::all()),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, mol_dsl!(expected));
    }

    #[rstest]
    #[case::hydrogen("F[C@H](Cl)Br", 1, vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
    ])]
    #[case::atoms("F[C@](Cl)(Br)I", 0, vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
    ])]
    #[case::lone_pair("C[S@](=O)CC", 0, vec![
        StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
    ])]
    fn test_resolver_project_stereo(
        #[case] input: &str,
        #[case] hydrogens: i64,
        #[case] ligands: Vec<StereoLigand>,
        #[values(false, true)] typing: bool,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    ) {
        let mut model = ChemistryModel {
            valence: if typing {
                ValenceModel::default()
            } else {
                ValenceModel::smiles()
            },
            ..Default::default()
        };
        model.valence.tie_break = policy;
        let config = ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        };
        let mut molecule = ingest_smiles(input).unwrap();
        assert_eq!(
            Resolver::with_config(&model, config).project(&mut molecule, ProjectFlags::all()),
            Ok(Solution::Determined(()))
        );
        assert_eq!(
            molecule.atom(AtomId(1)).implicit_hydrogens(),
            &NumForm::Lit(hydrogens)
        );
        assert_eq!(
            molecule.atom(AtomId(1)).constraints().tetrahedral_stereo(),
            Some(&TetrahedralStereoForm::Stereo(StereoCoset::Lit(0)))
        );
        assert_eq!(
            StereoPerception::new(&model.stereo)
                .derive_stereo_atom(&molecule, AtomId(1), &StereoCoset::Lit(0))
                .map(|(ligands, _)| ligands),
            Some(ligands)
        );
    }

    #[rstest]
    #[case::empty(Molecule::new())]
    #[case::ordinary(mol_dsl_concrete!(r#"{:atoms ["N#h2#n1" "C#i13#c-#h3#n1" "C#h3#u1#s2"]}"#))]
    #[case::unresolved_order(mol_dsl!(r#"{:atoms ["C#i=" "C#i="] :bonds [[0 1 "*#c0#u0#s"]]}"#))]
    #[case::retained_relations(mol_dsl_concrete!(r#"{:atoms ["N#h3#n1" "B" "H" "B"]
        :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]
        :multicenter-bonds [{:atoms [1 2 3] :attrs "[1,0,1]#c0#u0#s"}]
        :noncovalent-bonds [{:atoms [0 3] :attrs "*"}]}"#))]
    fn test_resolver_project_identity(#[case] mut molecule: Molecule) {
        let model = ChemistryModel::default();
        let original = molecule.clone();
        assert_eq!(
            Resolver::new(&model).project(&mut molecule, ProjectFlags::all()),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::localized_charge(mol_dsl_concrete!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+"]]}"#),
        ProjectError::BondCharge { entity: Entity::Bond(BondId(0)), charge: NumForm::Lit(1) })]
    #[case::localized_negative_charge(mol_dsl_concrete!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c-"]]}"#),
        ProjectError::BondCharge { entity: Entity::Bond(BondId(0)), charge: NumForm::Lit(-1) })]
    #[case::localized_charge_unknown(mol_dsl!(r#"{:atoms ["C#i=" "C#i="] :bonds [[0 1 "1#u0#s"]]}"#),
        ProjectError::BondCharge { entity: Entity::Bond(BondId(0)), charge: NumForm::Undetermined })]
    #[case::localized_radical(mol_dsl_concrete!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#u1#s2"]]}"#),
        ProjectError::BondSpin { entity: Entity::Bond(BondId(0)), spin: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } })]
    #[case::localized_multiplicity(mol_dsl_concrete!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#u0#s3"]]}"#),
        ProjectError::BondSpin { entity: Entity::Bond(BondId(0)), spin: UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Lit(3) } })]
    #[case::localized_spin_unknown(mol_dsl!(r#"{:atoms ["C#i=" "C#i="] :bonds [[0 1 "1#c0"]]}"#),
        ProjectError::BondSpin { entity: Entity::Bond(BondId(0)), spin: UnpairedElectronsForm::default() })]
    #[case::multicenter_charge(mol_dsl_concrete!(r#"{:atoms ["B" "H" "B"] :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,0,1]#c-"}]}"#),
        ProjectError::BondCharge { entity: Entity::MulticenterBond(MulticenterBondId(0)), charge: NumForm::Lit(-1) })]
    #[case::multicenter_spin(mol_dsl_concrete!(r#"{:atoms ["B" "H" "B"] :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,0,1]#u1#s2"}]}"#),
        ProjectError::BondSpin { entity: Entity::MulticenterBond(MulticenterBondId(0)), spin: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) } })]
    fn test_resolver_project_bond_error(
        #[case] mut molecule: Molecule,
        #[case] expected: ProjectError,
    ) {
        let model = ChemistryModel::default();
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                isotope: IsotopePolicy::Natural,
                ..Default::default()
            },
        );
        let original = molecule.clone();
        assert_eq!(
            resolver.project(&mut molecule, ProjectFlags::all()),
            Err(expected)
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::stereo(StereoConfigurationForm::Undetermined, NumForm::Lit(0), UnpairedElectronsForm::closed_shell(), IsotopeMassForm::Natural,
        ProjectError::Stereo(StereoProjectError::UnsupportedStereoAtom { stereo_atom: StereoAtomId(0) }))]
    #[case::aromatic_charge(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0), NumForm::Lit(1), UnpairedElectronsForm::closed_shell(), IsotopeMassForm::Natural,
        ProjectError::Aromaticity(AromaticityProjectError::ChargedSystem { system: AromaticSystemId(0) }))]
    #[case::aromatic_spin(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0), NumForm::Lit(0), UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, IsotopeMassForm::Natural,
        ProjectError::Aromaticity(AromaticityProjectError::SystemSpin { system: AromaticSystemId(0) }))]
    #[case::isotope(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0), NumForm::Lit(0), UnpairedElectronsForm::closed_shell(), IsotopeMassForm::Undetermined,
        ProjectError::Isotope(IsotopeProjectError::NonGroundIsotope { atom: AtomId(0) }))]
    fn test_resolver_project_phase_error(
        #[case] stereo: StereoConfigurationForm,
        #[case] charge: NumForm,
        #[case] spin: UnpairedElectronsForm,
        #[case] isotope: IsotopeMassForm,
        #[case] expected: ProjectError,
        #[values(ProjectFlags::all(), ProjectFlags::all() - ProjectFlags::VALENCE)]
        flags: ProjectFlags,
    ) {
        let source = mol_dsl_concrete!(
            r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"] [6 7 "1"] [6 8 "1"] [6 9 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
            :stereo-atoms [{:site 6 :ligands [7 8 9 [:h 6]] :attrs "Th0"}]}"#
        );
        let mut editor = source.edit();
        editor
            .stereo_atom_mut(StereoAtomId(0))
            .attributes
            .configuration = stereo;
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .attributes
            .charge = charge;
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .attributes
            .unpaired_electrons = spin;
        editor.atom_mut(AtomId(0)).attributes.isotope_mass = isotope;
        let mut molecule = editor.build();
        let original = molecule.clone();
        let model = ChemistryModel::default();
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                isotope: IsotopePolicy::Natural,
                ..Default::default()
            },
        );
        assert_eq!(resolver.project(&mut molecule, flags), Err(expected));
        assert_eq!(molecule, original);
    }
}
