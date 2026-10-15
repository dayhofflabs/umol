//! Composite resolver: chains the per-entity resolvers (isotope, valence,
//! aromaticity, stereo, bonds, multicenter bonds) on a single `Molecule`.
//!
//! A determined result satisfies [`Molecule::is_concrete`]. Consuming operations
//! use an editor; borrowed operations use a transaction and retain changes only
//! on a determined result.

pub mod aromaticity;
pub mod bonds;
pub mod isotope;
pub mod multicenter;
pub mod stereo;
pub mod valence;

use std::any::Any;
use std::collections::BTreeSet;

use aromaticity::AromaticityProjectError;
pub use aromaticity::{
    AromaticityFailurePolicy, AromaticityResolveConfig, AromaticityResolver, AtomSelectionAlgorithm,
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
    StereoContradiction, StereoError, StereoFailurePolicy, StereoResolveConfig, StereoResolver,
};
use thiserror::Error;
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomConstraintKey, AtomHandle, AtomId, Edits, Entity, Molecule,
    MoleculeApplyError, Normalize, NumForm, RollbackMethod, Transaction, UnpairedElectronsForm,
};
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;
use valence::ValenceProjectError;
pub use valence::{ValenceContradiction, ValenceError, ValenceResolver};

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
use crate::ops::model::{ChemistryModel, ValenceTieBreak};
use crate::ops::valence::compare::compare_by_key;
use crate::ops::valence::{AtomCompletions, ResolveReport};

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
/// optional tie-break information.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ResolveState {
    pub completions: AtomCompletions,
    pub aromatic_systems: Vec<(Vec<AtomId>, AromaticSystemForm)>,
    pub tie_breaks: Vec<AtomId>,
}

impl ResolveState {
    /// Returns the plural surviving completions and tie-break uses.
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

/// Resolves or projects a molecule using the configured chemistry model and policies.
///
/// # Semantic properties
///
/// resolve and project consume the input and return a molecule only on Determined.
/// resolve_into and project_into retain changes only on Determined; other outcomes
/// restore the receiver under Molecule::normalized_eq, preserving participant order.
/// Every phase plan applies to its input and preserves molecule integrity. Both routes
/// check that invariant before publishing and agree on chemistry outcomes and determined
/// molecules for the same input and configuration.
/// Each phase discharges only the entity assertions its results realize.
/// Other assertions and the molecule-level constraint list are preserved,
/// except for removals requested by failure policies and their reference compaction.
/// Resolution does not normalize, inline, or evaluate the molecule-level list.
/// Constraint validation and Molecule::inline_constraints are explicit operations.
/// The with_report methods additionally collect unresolved candidates and tie-break
/// uses; omitting that report preserves the outcome, errors, and resulting molecule.
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
    Valence(#[from] ValenceContradiction),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityContradiction),
    #[error(transparent)]
    Stereo(#[from] StereoContradiction),
    #[error(transparent)]
    Bonds(#[from] BondsContradiction),
    #[error(transparent)]
    MulticenterBonds(#[from] MulticenterBondsContradiction),
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
/// determined result.
#[derive(Debug, Clone, Default, PartialEq, Eq, Error)]
#[error("resolution underdetermined")]
pub struct ResolveUnderdetermined;

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
    /// result returns the molecule. Every other outcome drops the consumed input.
    /// Uses one editor without copying the molecule or recording undo.
    ///
    /// # Errors
    ///
    /// Returns the owning phase's planning error.
    pub fn resolve(
        &self,
        molecule: Molecule,
    ) -> Result<Solution<Molecule, ResolveContradiction, ()>, ResolveError> {
        let isotope_edits = match self.isotope.plan_resolve(&molecule) {
            Solution::Determined(edits) | Solution::Underdetermined(edits) => edits,
            Solution::Contradictory(contradiction) => match contradiction {},
        };
        let mut editor = molecule.edit();
        editor
            .apply(isotope_edits)
            .expect("isotope resolve plan applies to its input");
        let isotopic = editor
            .probe()
            .expect("resolver phase plans preserve molecule integrity");

        let state = match self
            .valence
            .admit(isotopic)
            .map_err(ResolveError::Valence)?
        {
            Solution::Determined(state) => state,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(()));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::from(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let outcome = self
            .aromaticity
            .select(isotopic, state, self.tie_break)
            .map_err(ResolveError::Aromaticity)?;
        let mut state = match outcome {
            Solution::Determined(state) => state,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::from(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };

        self.select_atom_completions(&mut state, false);

        if state
            .completions
            .iter()
            .any(|(_, candidates)| candidates.len() > 1)
        {
            return Ok(Solution::Underdetermined(()));
        }

        let edits = self.plan_constitution(isotopic, &state);
        editor
            .apply(edits)
            .expect("constitution plan applies to its input");
        let working = editor
            .probe()
            .expect("resolver phase plans preserve molecule integrity");

        let outcome = self
            .stereo
            .plan_resolve(working)
            .map_err(ResolveError::Stereo)?;
        let edits = match outcome {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(()));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::Stereo(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        // Plan stereo and bond defaults from the same post-constitution molecule.
        let bond_edits = self.bonds.plan_resolve(working);
        let multicenter_outcome = self.multicenter_bonds.plan_resolve(working);
        editor
            .apply(edits)
            .expect("stereo resolve plan applies to its input");
        editor
            .apply(bond_edits)
            .expect("bond resolve plan applies to its input");

        let edits = match multicenter_outcome {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(()));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::MulticenterBonds(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        editor
            .apply(edits)
            .expect("multicenter-bond resolve plan applies to its input");

        let molecule = editor
            .finish()
            .expect("resolver phase plans preserve molecule integrity");
        if molecule.is_concrete() {
            Ok(Solution::Determined(molecule))
        } else {
            Ok(Solution::Underdetermined(()))
        }
    }

    /// Resolves a consumed molecule and reports unresolved candidates and tie-break uses.
    ///
    /// An unresolved isotope does not stop later phases. Only a completely determined
    /// result returns the molecule. Every other outcome drops the consumed input.
    /// Underdetermined returns the report alone. Uses one editor without copying
    /// the molecule or recording undo.
    ///
    /// # Errors
    ///
    /// Returns the owning phase's planning error.
    pub fn resolve_with_report(
        &self,
        molecule: Molecule,
    ) -> Result<
        Solution<(Molecule, ResolveReport), ResolveContradiction, ResolveReport>,
        ResolveError,
    > {
        let isotope_edits = match self.isotope.plan_resolve(&molecule) {
            Solution::Determined(edits) | Solution::Underdetermined(edits) => edits,
            Solution::Contradictory(contradiction) => match contradiction {},
        };
        let mut editor = molecule.edit();
        editor
            .apply(isotope_edits)
            .expect("isotope resolve plan applies to its input");
        let isotopic = editor
            .probe()
            .expect("resolver phase plans preserve molecule integrity");

        let state = match self
            .valence
            .admit(isotopic)
            .map_err(ResolveError::Valence)?
        {
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
            .select_with_report(isotopic, state, self.tie_break)
            .map_err(ResolveError::Aromaticity)?;
        let mut state = match outcome {
            Solution::Determined(state) => state,
            Solution::Underdetermined(state) => {
                return Ok(Solution::Underdetermined(state.to_report()));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::from(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };

        self.select_atom_completions(&mut state, true);

        if state
            .completions
            .iter()
            .any(|(_, candidates)| candidates.len() > 1)
        {
            return Ok(Solution::Underdetermined(state.to_report()));
        }

        let edits = self.plan_constitution(isotopic, &state);
        editor
            .apply(edits)
            .expect("constitution plan applies to its input");
        let working = editor
            .probe()
            .expect("resolver phase plans preserve molecule integrity");

        let outcome = self
            .stereo
            .plan_resolve(working)
            .map_err(ResolveError::Stereo)?;
        let edits = match outcome {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(ResolveReport {
                    unresolved: AtomCompletions::new(),
                    tie_breaks: state.tie_breaks,
                }));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::Stereo(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        // Plan stereo and bond defaults from the same post-constitution molecule.
        let bond_edits = self.bonds.plan_resolve(working);
        let multicenter_outcome = self.multicenter_bonds.plan_resolve(working);
        editor
            .apply(edits)
            .expect("stereo resolve plan applies to its input");
        editor
            .apply(bond_edits)
            .expect("bond resolve plan applies to its input");

        let edits = match multicenter_outcome {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(ResolveReport {
                    unresolved: AtomCompletions::new(),
                    tie_breaks: state.tie_breaks,
                }));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolveContradiction::MulticenterBonds(contradiction);
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        editor
            .apply(edits)
            .expect("multicenter-bond resolve plan applies to its input");

        let molecule = editor
            .finish()
            .expect("resolver phase plans preserve molecule integrity");
        let report = ResolveReport {
            unresolved: AtomCompletions::new(),
            tie_breaks: state.tie_breaks,
        };
        if molecule.is_concrete() {
            Ok(Solution::Determined((molecule, report)))
        } else {
            Ok(Solution::Underdetermined(report))
        }
    }

    /// Resolves isotope composition before valence and aromaticity, then stereo and bonds.
    ///
    /// An unresolved isotope does not stop later phases. Only a completely determined
    /// result retains the changes; every other outcome preserves the caller's molecule.
    ///
    /// # Errors
    ///
    /// Returns the owning phase's planning error.
    pub fn resolve_into(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), ResolveContradiction>, ResolveError> {
        let isotope_edits = match self.isotope.plan_resolve(molecule) {
            Solution::Determined(edits) | Solution::Underdetermined(edits) => edits,
            Solution::Contradictory(contradiction) => match contradiction {},
        };
        Transaction::run::<_, MoleculeApplyError>(
            molecule,
            RollbackMethod::Journal,
            |mut transaction| {
                Ok((|| {
                    transaction
                        .apply(isotope_edits)
                        .expect("isotope resolve plan applies to its input");
                    let isotopic = transaction
                        .probe()
                        .expect("resolver phase plans preserve molecule integrity");

                    let state = match self
                        .valence
                        .admit(isotopic)
                        .map_err(ResolveError::Valence)?
                    {
                        Solution::Determined(state) => state,
                        Solution::Underdetermined(_) => {
                            return Ok(Solution::Underdetermined(()));
                        }
                        Solution::Contradictory(contradiction) => {
                            let contradiction = ResolveContradiction::from(contradiction);
                            return Ok(Solution::Contradictory(contradiction));
                        }
                    };
                    let outcome = self
                        .aromaticity
                        .select(isotopic, state, self.tie_break)
                        .map_err(ResolveError::Aromaticity)?;
                    let mut state = match outcome {
                        Solution::Determined(state) => state,
                        Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
                        Solution::Contradictory(contradiction) => {
                            let contradiction = ResolveContradiction::from(contradiction);
                            return Ok(Solution::Contradictory(contradiction));
                        }
                    };

                    self.select_atom_completions(&mut state, false);

                    if state
                        .completions
                        .iter()
                        .any(|(_, candidates)| candidates.len() > 1)
                    {
                        return Ok(Solution::Underdetermined(()));
                    }

                    let edits = self.plan_constitution(isotopic, &state);
                    transaction
                        .apply(edits)
                        .expect("constitution plan applies to its input");
                    let working = transaction
                        .probe()
                        .expect("resolver phase plans preserve molecule integrity");

                    let outcome = self
                        .stereo
                        .plan_resolve(working)
                        .map_err(ResolveError::Stereo)?;
                    let edits = match outcome {
                        Solution::Determined(edits) => edits,
                        Solution::Underdetermined(_) => {
                            return Ok(Solution::Underdetermined(()));
                        }
                        Solution::Contradictory(contradiction) => {
                            let contradiction = ResolveContradiction::Stereo(contradiction);
                            return Ok(Solution::Contradictory(contradiction));
                        }
                    };
                    // Plan stereo and bond defaults from the same post-constitution molecule.
                    let bond_edits = self.bonds.plan_resolve(working);
                    let multicenter_outcome = self.multicenter_bonds.plan_resolve(working);
                    transaction
                        .apply(edits)
                        .expect("stereo resolve plan applies to its input");
                    transaction
                        .apply(bond_edits)
                        .expect("bond resolve plan applies to its input");

                    let edits = match multicenter_outcome {
                        Solution::Determined(edits) => edits,
                        Solution::Underdetermined(_) => {
                            return Ok(Solution::Underdetermined(()));
                        }
                        Solution::Contradictory(contradiction) => {
                            let contradiction =
                                ResolveContradiction::MulticenterBonds(contradiction);
                            return Ok(Solution::Contradictory(contradiction));
                        }
                    };
                    transaction
                        .apply(edits)
                        .expect("multicenter-bond resolve plan applies to its input");

                    if transaction
                        .probe()
                        .expect("resolver phase plans preserve molecule integrity")
                        .is_concrete()
                    {
                        transaction
                            .commit()
                            .expect("resolver phase plans preserve molecule integrity");
                        Ok(Solution::Determined(()))
                    } else {
                        Ok(Solution::Underdetermined(()))
                    }
                })())
            },
        )
        .expect("resolver applies and commits only its own valid plans")
    }

    /// Resolves a borrowed molecule and reports unresolved candidates and tie-break uses.
    ///
    /// An unresolved isotope does not stop later phases. Only a completely determined
    /// result retains the changes; every other outcome preserves the caller's molecule.
    ///
    /// # Errors
    ///
    /// Returns the owning phase's planning error.
    pub fn resolve_into_with_report(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<ResolveReport, ResolveContradiction>, ResolveError> {
        let isotope_edits = match self.isotope.plan_resolve(molecule) {
            Solution::Determined(edits) | Solution::Underdetermined(edits) => edits,
            Solution::Contradictory(contradiction) => match contradiction {},
        };
        Transaction::run::<_, MoleculeApplyError>(
            molecule,
            RollbackMethod::Journal,
            |mut transaction| {
                Ok((|| {
                    transaction
                        .apply(isotope_edits)
                        .expect("isotope resolve plan applies to its input");
                    let isotopic = transaction
                        .probe()
                        .expect("resolver phase plans preserve molecule integrity");

                    let state = match self
                        .valence
                        .admit(isotopic)
                        .map_err(ResolveError::Valence)?
                    {
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
                        .select_with_report(isotopic, state, self.tie_break)
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

                    self.select_atom_completions(&mut state, true);

                    if state
                        .completions
                        .iter()
                        .any(|(_, candidates)| candidates.len() > 1)
                    {
                        return Ok(Solution::Underdetermined(state.to_report()));
                    }

                    let edits = self.plan_constitution(isotopic, &state);
                    transaction
                        .apply(edits)
                        .expect("constitution plan applies to its input");
                    let working = transaction
                        .probe()
                        .expect("resolver phase plans preserve molecule integrity");

                    let outcome = self
                        .stereo
                        .plan_resolve(working)
                        .map_err(ResolveError::Stereo)?;
                    let edits = match outcome {
                        Solution::Determined(edits) => edits,
                        Solution::Underdetermined(_) => {
                            return Ok(Solution::Underdetermined(ResolveReport {
                                unresolved: AtomCompletions::new(),
                                tie_breaks: state.tie_breaks,
                            }));
                        }
                        Solution::Contradictory(contradiction) => {
                            let contradiction = ResolveContradiction::Stereo(contradiction);
                            return Ok(Solution::Contradictory(contradiction));
                        }
                    };
                    // Plan stereo and bond defaults from the same post-constitution molecule.
                    let bond_edits = self.bonds.plan_resolve(working);
                    let multicenter_outcome = self.multicenter_bonds.plan_resolve(working);
                    transaction
                        .apply(edits)
                        .expect("stereo resolve plan applies to its input");
                    transaction
                        .apply(bond_edits)
                        .expect("bond resolve plan applies to its input");

                    let edits = match multicenter_outcome {
                        Solution::Determined(edits) => edits,
                        Solution::Underdetermined(_) => {
                            return Ok(Solution::Underdetermined(ResolveReport {
                                unresolved: AtomCompletions::new(),
                                tie_breaks: state.tie_breaks,
                            }));
                        }
                        Solution::Contradictory(contradiction) => {
                            let contradiction =
                                ResolveContradiction::MulticenterBonds(contradiction);
                            return Ok(Solution::Contradictory(contradiction));
                        }
                    };
                    transaction
                        .apply(edits)
                        .expect("multicenter-bond resolve plan applies to its input");

                    if transaction
                        .probe()
                        .expect("resolver phase plans preserve molecule integrity")
                        .is_concrete()
                    {
                        transaction
                            .commit()
                            .expect("resolver phase plans preserve molecule integrity");
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
                })())
            },
        )
        .expect("resolver applies and commits only its own valid plans")
    }

    /// Projects selected stages within graph IR, consuming the molecule.
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
    /// Determined returns the molecule with checked integrity. Rejection drops the input.
    /// Execution uses one editor without copying the molecule or recording undo.
    /// Projection does not invoke resolution or recovery comparisons. Implicit-H counts,
    /// atom electron fields, and localized charge are preserved.
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
        molecule: Molecule,
        flags: ProjectFlags,
    ) -> Result<Solution<Molecule, ProjectContradiction, ()>, ProjectError> {
        let bonds = molecule
            .bonds()
            .iter()
            .map(|bond| {
                (
                    Entity::Bond(bond.id()),
                    &bond.attributes().charge,
                    &bond.attributes().unpaired_electrons,
                )
            })
            .chain(molecule.multicenter_bonds().iter().map(|bond| {
                (
                    Entity::MulticenterBond(bond.id()),
                    &bond.attributes().charge,
                    &bond.attributes().unpaired_electrons,
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

        let mut first_edits = if flags.contains(ProjectFlags::STEREO) {
            Some(self.stereo.plan_project(&molecule)?)
        } else if flags.contains(ProjectFlags::AROMATICITY) {
            Some(self.aromaticity.plan_project(&molecule)?)
        } else if flags.contains(ProjectFlags::ISOTOPE) {
            Some(self.isotope.plan_project(&molecule)?)
        } else {
            None
        };

        let mut editor = molecule.edit();
        if flags.contains(ProjectFlags::STEREO) {
            let edits = match first_edits.take() {
                Some(edits) => edits,
                None => self.stereo.plan_project(
                    editor
                        .probe()
                        .expect("resolver phase plans preserve molecule integrity"),
                )?,
            };
            editor
                .apply(edits)
                .expect("stereo project plan applies to its input");
        }
        if flags.contains(ProjectFlags::AROMATICITY) {
            let edits = match first_edits.take() {
                Some(edits) => edits,
                None => self.aromaticity.plan_project(
                    editor
                        .probe()
                        .expect("resolver phase plans preserve molecule integrity"),
                )?,
            };
            editor
                .apply(edits)
                .expect("aromaticity project plan applies to its input");
        }
        // Valence projection preserves all fields and requires no edits.
        if flags.contains(ProjectFlags::ISOTOPE) {
            let edits = match first_edits.take() {
                Some(edits) => edits,
                None => self.isotope.plan_project(
                    editor
                        .probe()
                        .expect("resolver phase plans preserve molecule integrity"),
                )?,
            };
            editor
                .apply(edits)
                .expect("isotope project plan applies to its input");
        }
        let molecule = editor
            .finish()
            .expect("resolver phase plans preserve molecule integrity");
        Ok(Solution::Determined(molecule))
    }

    /// Projects selected stages, retaining changes only on Determined.
    ///
    /// See [Resolver::project] for phase order, preserved fields, and preconditions.
    ///
    /// # Semantic properties
    ///
    /// Errors, contradictions, and underdetermination restore the receiver under
    /// Molecule::normalized_eq, preserving participant order. The result agrees with
    /// consuming projection of the same input under the same configuration and flags.
    ///
    /// # Errors
    ///
    /// Returns the same bond charge/spin and phase errors as [Resolver::project].
    pub fn project_into(
        &self,
        molecule: &mut Molecule,
        flags: ProjectFlags,
    ) -> Result<Solution<(), ProjectContradiction>, ProjectError> {
        let bonds = molecule
            .bonds()
            .iter()
            .map(|bond| {
                (
                    Entity::Bond(bond.id()),
                    &bond.attributes().charge,
                    &bond.attributes().unpaired_electrons,
                )
            })
            .chain(molecule.multicenter_bonds().iter().map(|bond| {
                (
                    Entity::MulticenterBond(bond.id()),
                    &bond.attributes().charge,
                    &bond.attributes().unpaired_electrons,
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

        let mut first_edits = if flags.contains(ProjectFlags::STEREO) {
            Some(self.stereo.plan_project(molecule)?)
        } else if flags.contains(ProjectFlags::AROMATICITY) {
            Some(self.aromaticity.plan_project(molecule)?)
        } else if flags.contains(ProjectFlags::ISOTOPE) {
            Some(self.isotope.plan_project(molecule)?)
        } else {
            None
        };

        Transaction::run::<_, MoleculeApplyError>(
            molecule,
            RollbackMethod::Journal,
            |mut transaction| {
                Ok((|| {
                    if flags.contains(ProjectFlags::STEREO) {
                        let edits = match first_edits.take() {
                            Some(edits) => edits,
                            None => self.stereo.plan_project(
                                transaction
                                    .probe()
                                    .expect("resolver phase plans preserve molecule integrity"),
                            )?,
                        };
                        transaction
                            .apply(edits)
                            .expect("stereo project plan applies to its input");
                    }
                    if flags.contains(ProjectFlags::AROMATICITY) {
                        let edits = match first_edits.take() {
                            Some(edits) => edits,
                            None => self.aromaticity.plan_project(
                                transaction
                                    .probe()
                                    .expect("resolver phase plans preserve molecule integrity"),
                            )?,
                        };
                        transaction
                            .apply(edits)
                            .expect("aromaticity project plan applies to its input");
                    }
                    // Valence projection preserves all fields and requires no edits.
                    if flags.contains(ProjectFlags::ISOTOPE) {
                        let edits = match first_edits.take() {
                            Some(edits) => edits,
                            None => self.isotope.plan_project(
                                transaction
                                    .probe()
                                    .expect("resolver phase plans preserve molecule integrity"),
                            )?,
                        };
                        transaction
                            .apply(edits)
                            .expect("isotope project plan applies to its input");
                    }
                    transaction
                        .commit()
                        .expect("resolver phase plans preserve molecule integrity");
                    Ok(Solution::Determined(()))
                })())
            },
        )
        .expect("resolver applies and commits only its own valid plans")
    }

    fn select_atom_completions(&self, state: &mut ResolveState, report: bool) {
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
                    if report {
                        state.tie_breaks.push(atom);
                    }
                }
            }
            if report {
                state.tie_breaks.sort_unstable();
                state.tie_breaks.dedup();
            }
        }
    }

    fn plan_constitution(&self, molecule: &Molecule, state: &ResolveState) -> Edits {
        let mut edits = self
            .aromaticity
            .plan_bond_discharge(molecule, &state.aromatic_systems);
        let mut elements = Vec::new();
        let mut isotope_masses = Vec::new();
        let mut charges = Vec::new();
        let mut implicit_hydrogens = Vec::new();
        let mut lone_pairs = Vec::new();
        let mut unpaired_electrons = Vec::new();
        let mut valences = Vec::new();
        for (atom, disjuncts) in state.completions.iter() {
            let current = molecule.atom(atom).attributes();
            let selected = &disjuncts[0];
            if !current.element.normalized_eq(&selected.element) {
                elements.push((
                    AtomHandle::Id(atom),
                    current.element.clone(),
                    selected.element.clone(),
                ));
            }
            if !current.isotope_mass.normalized_eq(&selected.isotope_mass) {
                isotope_masses.push((
                    AtomHandle::Id(atom),
                    current.isotope_mass.clone(),
                    selected.isotope_mass.clone(),
                ));
            }
            if !current.charge.normalized_eq(&selected.charge) {
                charges.push((
                    AtomHandle::Id(atom),
                    current.charge.clone(),
                    selected.charge.clone(),
                ));
            }
            if !current
                .implicit_hydrogens
                .normalized_eq(&selected.implicit_hydrogens)
            {
                implicit_hydrogens.push((
                    AtomHandle::Id(atom),
                    current.implicit_hydrogens.clone(),
                    selected.implicit_hydrogens.clone(),
                ));
            }
            if !current.lone_pairs.normalized_eq(&selected.lone_pairs) {
                lone_pairs.push((
                    AtomHandle::Id(atom),
                    current.lone_pairs.clone(),
                    selected.lone_pairs.clone(),
                ));
            }
            let selected_unpaired_electrons = current.unpaired_electrons.update(
                &current
                    .unpaired_electrons
                    .difference_to(&selected.unpaired_electrons),
            );
            if !current
                .unpaired_electrons
                .normalized_eq(&selected_unpaired_electrons)
            {
                unpaired_electrons.push((
                    AtomHandle::Id(atom),
                    current.unpaired_electrons.clone(),
                    selected_unpaired_electrons,
                ));
            }
            if let Some(asserted) = current.constraints.get(AtomConstraintKey::Valence) {
                valences.push((AtomHandle::Id(atom), Some(asserted.clone()), None));
            }
        }
        if !elements.is_empty() {
            edits.modify_atom_elements(elements);
        }
        if !isotope_masses.is_empty() {
            edits.modify_atom_isotope_masses(isotope_masses);
        }
        if !charges.is_empty() {
            edits.modify_atom_charges(charges);
        }
        if !implicit_hydrogens.is_empty() {
            edits.modify_atom_implicit_hydrogens(implicit_hydrogens);
        }
        if !lone_pairs.is_empty() {
            edits.modify_atom_lone_pairs(lone_pairs);
        }
        if !unpaired_electrons.is_empty() {
            edits.modify_atom_unpaired_electrons(unpaired_electrons);
        }
        if !valences.is_empty() {
            edits.modify_atom_constraints(AtomConstraintKey::Valence, valences);
        }
        let existing: BTreeSet<Vec<AtomId>> = molecule
            .aromatic_systems()
            .iter()
            .map(|aromatic_system| {
                let mut atoms: Vec<AtomId> = aromatic_system.atom_ids().collect();
                atoms.sort_unstable();
                atoms
            })
            .collect();
        for (atoms, aromatic_system) in &state.aromatic_systems {
            let mut key = atoms.clone();
            key.sort_unstable();
            if existing.contains(&key) {
                continue;
            }
            for edit in self.aromaticity.perception.plan_aromatic_system(
                molecule,
                atoms.clone(),
                aromatic_system.clone(),
            ) {
                edits.push(edit);
            }
        }
        edits
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::{fixture, rstest};
    use smallvec::smallvec;
    use umol_chem::element::Element;
    use umol_graph_core::Correspondence;
    use umol_graph_ir::ir::{
        AromaticSystemConstraintForm, AromaticSystemId, AtomConstraintForm, AtomForm, AtomId,
        BondId, Constraint, Edit, IsotopeMassForm, MoleculeConstraint, MoleculeCorrespondence,
        MoleculeEntries, MulticenterBondId, MulticenterValenceForm, NumForm,
        StereoAtomConstraintForm, StereoAtomId, StereoBondId, StereoConfigurationForm, StereoCoset,
        StereoKind, StereoLigand, StereoLigandKind, Stereogenicity, StereogenicityForm,
        TetrahedralStereoForm,
    };
    use umol_graph_ir::{atom_dsl, bond_dsl, mol_dsl, mol_dsl_concrete};

    use super::*;
    use crate::ingest::ingest_smiles;
    use crate::ops::aromaticity::{AromaticityError, AromaticityInconsistency};
    use crate::ops::model::{
        AromaticityModel, AromaticityRule, AromaticityTieBreak, ChemistryModel, ElementScope,
        RingLimits, StereoModel, ValenceModel,
    };
    use crate::ops::stereo::{StereoInconsistency, StereoPerception};
    use crate::ops::valence::{AtomTypeRegistry, CountsError, ValenceTable};
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

    #[fixture]
    fn indole_skeleton() -> Molecule {
        mol_dsl!(
            r#"{
            :atoms ["N#i=#c0#u0#s#a+" "C#i=#c0#u0#s#a" "C#i=#c0#u0#s#a"
                    "C#i=#c0#u0#s#a" "C#i=#c0#u0#s#a" "C#i=#c0#u0#s#a"
                    "C#i=#c0#u0#s#a" "C#i=#c0#u0#s#a" "C#i=#c0#u0#s#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]
                    [3 5 "1"] [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 4 "1"]]
        }"#
        )
    }

    #[rstest]
    fn test_resolve_config_default() {
        assert_eq!(
            ResolveConfig::default(),
            ResolveConfig {
                isotope: IsotopePolicy::Strict,
                aromaticity: AromaticityResolveConfig::default(),
                stereo: StereoResolveConfig::default(),
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
            resolver.aromaticity.plan_resolve(&aromatic_molecule),
            explicit.aromaticity.plan_resolve(&aromatic_molecule)
        );
        assert_eq!(
            resolver.stereo.plan_resolve(&stereo_molecule),
            explicit.stereo.plan_resolve(&stereo_molecule)
        );
    }

    #[rstest]
    #[case::keep_stereo_failure(
        ResolveConfig {
            isotope: IsotopePolicy::Strict,
            aromaticity: AromaticityResolveConfig::default(),
            stereo: StereoResolveConfig {
                tetrahedral_stereo_failure: StereoFailurePolicy::Keep,
                ..StereoResolveConfig::default()
            },
        },
        mol_dsl_concrete!(r#"{:atoms ["C#h4#T1"]}"#)
    )]
    fn test_resolver_with_config(
        chemistry_model: ChemistryModel,
        aromatic_molecule: Molecule,
        #[case] config: ResolveConfig,
        #[case] stereo_molecule: Molecule,
    ) {
        let resolver = Resolver::with_config(&chemistry_model, config);
        let expected_aromaticity =
            AromaticityResolver::with_config(&chemistry_model.aromaticity, config.aromaticity)
                .plan_resolve(&aromatic_molecule);
        let expected_stereo = StereoResolver::with_config(&chemistry_model.stereo, config.stereo)
            .plan_resolve(&stereo_molecule);

        assert_eq!(
            resolver.aromaticity.plan_resolve(&aromatic_molecule),
            expected_aromaticity
        );
        assert_eq!(
            resolver.stereo.plan_resolve(&stereo_molecule),
            expected_stereo
        );
        if config.stereo != StereoResolveConfig::default() {
            assert_ne!(
                expected_stereo,
                StereoResolver::new(&chemistry_model.stereo).plan_resolve(&stereo_molecule)
            );
        }
    }

    #[rstest]
    #[case::selected(
        mol_dsl!(r#"{:atoms ["C#i=#c0"]}"#),
        ValenceTieBreak::MostSaturated,
        Solution::Determined((
            mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#),
            ResolveReport { unresolved: AtomCompletions::new(), tie_breaks: vec![AtomId(0)] },
        ))
    )]
    #[case::plural(
        mol_dsl!(r#"{:atoms ["C#i=#c0"]}"#),
        ValenceTieBreak::Strict,
        Solution::Underdetermined(ResolveReport {
            unresolved: AtomCompletions::from_iter([(AtomId(0), vec![
                atom_dsl!("C#i=#c0#h4#n0#u0#s"),
                atom_dsl!("C#i=#c0#h3#n0#u1#s2"),
            ].into())]),
            tie_breaks: Vec::new(),
        })
    )]
    #[case::isotope(
        mol_dsl!(r#"{:atoms ["C#c0"]}"#),
        ValenceTieBreak::MostSaturated,
        Solution::Underdetermined(ResolveReport {
            unresolved: AtomCompletions::new(), tie_breaks: vec![AtomId(0)],
        })
    )]
    #[case::element(
        mol_dsl!(r#"{:atoms ["{C,N}#c0"]}"#),
        ValenceTieBreak::MostSaturated,
        Solution::Underdetermined(ResolveReport::default())
    )]
    fn test_resolver_resolve_with_report(
        #[case] molecule: Molecule,
        #[case] tie_break: ValenceTieBreak,
        #[case] expected: Solution<(Molecule, ResolveReport), ResolveContradiction, ResolveReport>,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel {
                tie_break,
                ..ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([
                    atom_dsl!("C#c0#h4#n0#u0#s"),
                    atom_dsl!("C#c0#h3#n0#u1#s2"),
                ])))
            },
            ..ChemistryModel::default()
        };
        let resolver = Resolver::new(&model);
        assert_eq!(
            resolver.resolve_with_report(molecule.clone()),
            Ok(expected.clone())
        );
        let (expected_molecule, expected_report, expected_owned) = match expected {
            Solution::Determined((resolved, report)) => (
                resolved.clone(),
                Solution::Determined(report),
                Solution::Determined(resolved),
            ),
            Solution::Underdetermined(report) => (
                molecule.clone(),
                Solution::Underdetermined(report),
                Solution::Underdetermined(()),
            ),
            Solution::Contradictory(contradiction) => (
                molecule.clone(),
                Solution::Contradictory(contradiction.clone()),
                Solution::Contradictory(contradiction),
            ),
        };
        assert_eq!(resolver.resolve(molecule.clone()), Ok(expected_owned));
        let mut with_report = molecule.clone();
        assert_eq!(
            resolver.resolve_into_with_report(&mut with_report),
            Ok(expected_report.clone())
        );
        assert_eq!(with_report, expected_molecule);
        let mut without_report = molecule;
        assert_eq!(
            resolver.resolve_into(&mut without_report),
            Ok(expected_report.map(|_| ()))
        );
        assert_eq!(without_report, expected_molecule);
    }

    #[rstest]
    #[case::higher_threshold(0.2)]
    #[case::lower_threshold(0.1)]
    fn test_resolver_resolve_with_report_hmo_error(
        indole_skeleton: Molecule,
        #[case] stabilization_threshold: f64,
        #[values(
            AtomSelectionAlgorithm::FullDfs,
            AtomSelectionAlgorithm::BacktrackingDfs
        )]
        atom_selection_algorithm: AtomSelectionAlgorithm,
    ) {
        let model = ChemistryModel {
            aromaticity: AromaticityModel {
                scope: ElementScope::Any,
                rule: AromaticityRule::Hmo {
                    stabilization_threshold,
                },
                tie_break: AromaticityTieBreak::MinElectronCount,
            },
            ..ChemistryModel::default()
        };
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                aromaticity: AromaticityResolveConfig {
                    atom_selection_algorithm,
                    ..AromaticityResolveConfig::default()
                },
                ..ResolveConfig::default()
            },
        );
        assert_eq!(
            resolver.resolve_with_report(indole_skeleton),
            Ok(Solution::Contradictory(ResolveContradiction::Aromaticity(
                AromaticityContradiction::HmoInvalidInput(
                    "open-shell pi-system (odd electron count) not supported by HMO".to_owned(),
                ),
            ))),
        );
    }

    #[rstest]
    #[case::full_at_bound(
        12,
        AtomSelectionAlgorithm::FullDfs,
        AromaticityFailurePolicy::Error,
        true
    )]
    #[case::full_above_bound(
        13,
        AtomSelectionAlgorithm::FullDfs,
        AromaticityFailurePolicy::Error,
        false
    )]
    #[case::pruned_above_product_bound(
        13,
        AtomSelectionAlgorithm::BacktrackingDfs,
        AromaticityFailurePolicy::Error,
        true
    )]
    #[case::keep_full(
        13,
        AtomSelectionAlgorithm::FullDfs,
        AromaticityFailurePolicy::Keep,
        false
    )]
    #[case::keep_backtracking(
        13,
        AtomSelectionAlgorithm::BacktrackingDfs,
        AromaticityFailurePolicy::Keep,
        false
    )]
    fn test_resolver_resolve_search_bound(
        #[case] nitrogen_count: usize,
        #[case] atom_selection_algorithm: AtomSelectionAlgorithm,
        #[case] aromatic_valence_failure: AromaticityFailurePolicy,
        #[case] determined: bool,
    ) {
        let bonds: Vec<_> = (0..14)
            .map(|atom| {
                (
                    AtomId(atom),
                    AtomId((atom + 1) % 14),
                    bond_dsl!("1#c0#u0#s"),
                )
            })
            .collect();
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: (0..14)
                .map(|atom| {
                    if atom < nitrogen_count {
                        atom_dsl!("N#i=#c0#u0#s#a+")
                    } else {
                        atom_dsl!("C#i=#c0#u0#s#a")
                    }
                })
                .collect(),
            bonds: bonds.clone(),
            ..MoleculeEntries::default()
        });
        let mut model = ChemistryModel::default();
        model.aromaticity.tie_break = AromaticityTieBreak::MinElectronCount;
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                aromaticity: AromaticityResolveConfig {
                    atom_selection_algorithm,
                    aromatic_valence_failure,
                    ..AromaticityResolveConfig::default()
                },
                ..ResolveConfig::default()
            },
        );
        let Solution::Determined(admitted) = resolver.valence.admit(&molecule).unwrap() else {
            panic!("the ring has complete valence candidates");
        };
        assert_eq!(
            admitted
                .completions
                .iter()
                .map(|(_, forms)| forms.len())
                .product::<usize>(),
            1 << nitrogen_count
        );
        let expected = if determined {
            Solution::Determined(Molecule::from_entries(MoleculeEntries {
                atoms: (0..14)
                    .map(|atom| {
                        if atom < nitrogen_count {
                            atom_dsl!("N#i=#c0#h0#n1#u0#s")
                        } else {
                            atom_dsl!("C#i=#c0#h1#n0#u0#s")
                        }
                    })
                    .collect(),
                bonds,
                aromatic: vec![(
                    (0..14).map(AtomId).collect(),
                    AromaticSystemForm::from_electrons(vec![1; 14])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                )],
                ..MoleculeEntries::default()
            }))
        } else {
            Solution::Underdetermined(())
        };
        assert_eq!(resolver.resolve(molecule), Ok(expected));
    }

    #[rstest]
    #[case::counts(ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())))]
    #[case::atom_typing(ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
            "C#c0#h4#n0#u0#s#v0#a!"
        )]))))]
    fn test_resolver_resolve_into(#[case] valence: ValenceModel) {
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
            Resolver::new(&model).resolve_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#a!"]}"#));
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
    fn test_resolver_resolve_into_atoms(
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
            Resolver::new(&model).resolve_into(&mut molecule),
            Ok(Solution::Determined(()))
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
    fn test_resolver_resolve_into_overlap(
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
            resolver.resolve_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#));
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_concrete!(r#"{
            :atoms ["C#T1" "F#n3" "Cl#n3" "Br#n3" "I#n3"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
        }"#),
        StereoInconsistency::TetrahedralStereoMismatch {
            atom: AtomId(0), stereo_atom: StereoAtomId(0),
        },
    )]
    #[case::cis_trans(
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "F#n3" "Cl#n3" "Br#n3" "I#n3"]
            :bonds [[0 1 "2#C1"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]
        }"#),
        StereoInconsistency::CisTransStereoMismatch {
            bond: BondId(0), stereo_bond: StereoBondId(0),
        },
    )]
    fn test_resolver_resolve_stereo_mismatch(
        chemistry_model: ChemistryModel,
        #[case] mut molecule: Molecule,
        #[case] inconsistency: StereoInconsistency,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        failure: StereoFailurePolicy,
    ) {
        let original = molecule.clone();
        let resolver = Resolver::with_config(
            &chemistry_model,
            ResolveConfig {
                stereo: StereoResolveConfig {
                    tetrahedral_stereo_failure: failure,
                    stereo_atom_failure: failure,
                    cis_trans_stereo_failure: failure,
                    stereo_bond_failure: failure,
                },
                ..ResolveConfig::default()
            },
        );
        let expected =
            ResolveContradiction::Stereo(StereoContradiction::Inconsistency(inconsistency));
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Contradictory(expected.clone()))
        );
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
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
        #[values(false, true)] consuming: bool,
        #[values(false, true)] shared: bool,
    ) {
        let retained = shared.then(|| molecule.clone());
        let resolver = Resolver::with_config(
            &chemistry_model,
            ResolveConfig {
                isotope: IsotopePolicy::Natural,
                ..Default::default()
            },
        );
        if consuming {
            assert_eq!(
                resolver.resolve(molecule),
                Ok(Solution::Determined(expected))
            );
        } else {
            assert_eq!(
                resolver.resolve_into(&mut molecule),
                Ok(Solution::Determined(()))
            );
            assert_eq!(molecule, expected);
        }
        drop(retained);
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
        let consumed = resolver.resolve(molecule.clone());
        let actual = resolver.resolve_into(&mut molecule);
        if omitted && isotope == IsotopePolicy::Strict {
            assert_eq!(consumed, Ok(Solution::Underdetermined(())));
            assert_eq!(actual, Ok(Solution::Underdetermined(())));
            assert_eq!(molecule, original);
        } else {
            assert_eq!(consumed, Ok(Solution::Determined(expected.clone())));
            assert_eq!(actual, Ok(Solution::Determined(())));
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
            .resolve(molecule.clone()),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(
            Resolver::with_config(
                &model,
                ResolveConfig {
                    isotope,
                    ..Default::default()
                }
            )
            .resolve_into(&mut molecule),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::undetermined(mol_dsl!(r#"{:atoms ["C#c0#h4#T1"]}"#))]
    #[case::variable(mol_dsl!(r#"{:atoms ["C#i?mass#c0#h4#T1"]}"#))]
    fn test_resolver_resolve_into_isotope_stereo(
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
            .resolve_into(&mut molecule),
            Ok(Solution::Contradictory(ResolveContradiction::Stereo(
                StereoContradiction::Inconsistency(StereoInconsistency::TetrahedralStereoFailure {
                    atom: AtomId(0)
                })
            )))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::atom(mol_dsl_concrete!(r#"{
        :atoms ["C#h" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1#o(0,1)="}]
    }"#))]
    #[case::bond(mol_dsl_concrete!(r#"{
        :atoms ["C" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1#o(0,1)="}]
    }"#))]
    fn test_resolver_resolve_stereo_assertions(
        chemistry_model: ChemistryModel,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            Resolver::new(&chemistry_model).resolve(molecule.clone()),
            Ok(Solution::Determined(molecule))
        );
    }

    #[rstest]
    #[case::matching_degree(mol_dsl_concrete!(r#"{:atoms ["C#h4#D0"]}"#))]
    #[case::conflicting_degree(mol_dsl_concrete!(r#"{:atoms ["C#h4#D1"]}"#))]
    #[case::uncompleted_valence(mol_dsl_concrete!(r#"{:atoms ["C#h4#v0"]}"#))]
    #[case::bond_ring(mol_dsl_concrete!(r#"{
        :atoms ["C#h3" "C#h3"] :bonds [[0 1 "1#R!"]]
    }"#))]
    #[case::stored_aromatic(mol_dsl_concrete!(r#"{
        :atoms ["C#h#a1" "C#h" "C#h" "C#h" "C#h" "C#h"]
        :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
    }"#))]
    #[case::stored_tetrahedral(mol_dsl_concrete!(r#"{
        :atoms ["C#h#T1" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th1"}]
    }"#))]
    #[case::stored_cis_trans(mol_dsl_concrete!(r#"{
        :atoms ["C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "2#C1"] [0 2 "1"] [1 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:h 1]] :attrs "Ct1"}]
    }"#))]
    fn test_resolver_resolve_constraints_identity(
        chemistry_model: ChemistryModel,
        #[case] molecule: Molecule,
    ) {
        let resolver = Resolver::new(&chemistry_model);
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Determined(molecule.clone()))
        );
        assert_eq!(
            resolver.resolve_with_report(molecule.clone()),
            Ok(Solution::Determined((
                molecule.clone(),
                ResolveReport::default()
            )))
        );
        let mut actual = molecule.clone();
        assert_eq!(
            resolver.resolve_into(&mut actual),
            Ok(Solution::Determined(()))
        );
        assert_eq!(actual, molecule);
        assert_eq!(
            resolver.resolve_into_with_report(&mut actual),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(actual, molecule);
    }

    #[rstest]
    #[case::true_sum(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: None, sum: NumForm::Lit(0),
    }))]
    #[case::false_sum(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: None, sum: NumForm::Lit(1),
    }))]
    #[case::open_sum(Constraint::Molecule(MoleculeConstraint::ChargeSum {
        atoms: None, sum: NumForm::LitSet(Box::new(BTreeSet::from([0, 1]))),
    }))]
    #[case::entity(Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::Lit(0))))]
    #[case::collision(Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::Lit(3))))]
    #[case::singleton(Constraint::Or(vec![Constraint::Atom(
        AtomId(0), AtomConstraintForm::Valence(NumForm::Lit(0)),
    )]))]
    #[case::conjunction(Constraint::And(vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::Lit(0))),
        Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::Lit(3))),
    ]))]
    #[case::disjunction(Constraint::Or(vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::Lit(0))),
        Constraint::Atom(AtomId(0), AtomConstraintForm::Valence(NumForm::Lit(3))),
    ]))]
    #[case::negation(Constraint::Not(Box::new(Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::Valence(NumForm::Lit(0)),
    ))))]
    fn test_resolver_resolve_molecule_constraints_identity(
        chemistry_model: ChemistryModel,
        #[case] constraint: Constraint,
    ) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom_dsl!("C#i=#c0#h4#n0#u0#s#v0")],
            constraints: vec![constraint].into(),
            ..MoleculeEntries::default()
        });
        let resolver = Resolver::new(&chemistry_model);
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Determined(molecule.clone()))
        );
        assert_eq!(
            resolver.resolve_with_report(molecule.clone()),
            Ok(Solution::Determined((
                molecule.clone(),
                ResolveReport::default()
            )))
        );
        let mut actual = molecule.clone();
        assert_eq!(
            resolver.resolve_into(&mut actual),
            Ok(Solution::Determined(()))
        );
        assert_eq!(actual, molecule);
        assert_eq!(
            resolver.resolve_into_with_report(&mut actual),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(actual, molecule);
    }

    #[rstest]
    fn test_resolver_resolve_into_pyrrolyl() {
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
            Resolver::new(&model).resolve_into(&mut molecule),
            Ok(Solution::Determined(()))
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
            .resolve(molecule.clone()),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(
            Resolver::with_config(
                &chemistry_model,
                ResolveConfig {
                    isotope: IsotopePolicy::Natural,
                    ..Default::default()
                }
            )
            .resolve_into(&mut molecule),
            Ok(Solution::Underdetermined(()))
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
            Resolver::new(&model).resolve(molecule.clone()),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(
            Resolver::new(&model).resolve_into(&mut molecule),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::isotope_and_hydrogen_completion(mol_dsl!(r#"{
            :atoms ["C#i*#c0#h*#n0#u0#s#T+" "F#i=#c0#h0#n0#u0#s"
                    "Cl#i=#c0#h0#n0#u0#s" "Br#i=#c0#h0#n0#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"]]
        }"#))]
    #[case::tetrahedral_open(mol_dsl_concrete!(r#"{
        :atoms ["C#h#T+" "F#n3" "Cl#n3" "Br#n3"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    }"#))]
    #[case::cis_trans_open(mol_dsl_concrete!(r#"{
        :atoms ["C#h" "C#h" "F#n3" "Cl#n3"]
        :bonds [[0 1 "2#C+"] [0 2 "1"] [1 3 "1"]]
    }"#))]
    #[case::stored_open(mol_dsl_concrete!(r#"{
        :atoms ["C#h" "C#h" "F#n3" "Cl#n3"]
        :bonds [[0 1 "2#C+"] [0 2 "1"] [1 3 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:h 1]] :attrs "Ct*"}]
    }"#))]
    fn test_resolver_resolve_later_underdetermined(
        chemistry_model: ChemistryModel,
        #[case] mut molecule: Molecule,
    ) {
        let original = molecule.clone();

        assert_eq!(
            Resolver::with_config(
                &chemistry_model,
                ResolveConfig {
                    isotope: IsotopePolicy::Natural,
                    ..Default::default()
                }
            )
            .resolve(molecule.clone()),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(
            Resolver::with_config(
                &chemistry_model,
                ResolveConfig {
                    isotope: IsotopePolicy::Natural,
                    ..Default::default()
                }
            )
            .resolve_into(&mut molecule),
            Ok(Solution::Underdetermined(()))
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

        assert_eq!(resolver.resolve(molecule.clone()), Err(expected.clone()));
        assert_eq!(
            resolver.resolve_with_report(molecule.clone()),
            Err(expected.clone())
        );
        assert_eq!(
            resolver.resolve_into_with_report(&mut molecule),
            Err(expected.clone())
        );
        assert_eq!(molecule, original);
        assert_eq!(resolver.resolve_into(&mut molecule), Err(expected));
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
        Solution::Underdetermined(())
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
        #[case] expected: Solution<(), ResolveContradiction>,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([
                atom_dsl!("C#c0#h0#n0#u0#s#v0#a!#m1"),
            ]))),
            ..ChemistryModel::default()
        };
        let original = molecule.clone();
        let consumed = Resolver::new(&model).resolve(molecule.clone());
        assert_eq!(
            consumed,
            Ok(match &expected {
                Solution::Determined(_) => Solution::Determined(molecule.clone()),
                Solution::Underdetermined(_) => Solution::Underdetermined(()),
                Solution::Contradictory(contradiction) =>
                    Solution::Contradictory(contradiction.clone()),
            })
        );

        assert_eq!(
            Resolver::new(&model).resolve_into(&mut molecule),
            Ok(expected)
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::atom_aromatic_system(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a2" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        ResolveContradiction::Valence(ValenceContradiction::Counts(CountsError::NoMatch)),
    )]
    #[case::atom_bond(
        mol_dsl_concrete!(r#"{:atoms ["C#h3#a!" "C#h3#a"] :bonds [[0 1 "1#a"]]}"#),
        ResolveContradiction::Valence(ValenceContradiction::Counts(CountsError::NoMatch)),
    )]
    #[case::bond_aromatic_system(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1#a!"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#),
        ResolveContradiction::Aromaticity(AromaticityInconsistency::AromaticBondConstraintMismatch {
            bond: BondId(0), aromatic_system: AromaticSystemId(0),
        }.into()),
    )]
    fn test_resolver_resolve_aromatic_mismatch(
        chemistry_model: ChemistryModel,
        #[case] mut molecule: Molecule,
        #[case] expected: ResolveContradiction,
        #[values(AromaticityFailurePolicy::Error, AromaticityFailurePolicy::Keep)]
        failure: AromaticityFailurePolicy,
    ) {
        let original = molecule.clone();
        let resolver = Resolver::with_config(
            &chemistry_model,
            ResolveConfig {
                aromaticity: AromaticityResolveConfig {
                    aromatic_valence_failure: failure,
                    aromatic_system_failure: failure,
                    ..AromaticityResolveConfig::default()
                },
                ..ResolveConfig::default()
            },
        );
        assert_eq!(
            resolver.resolve(molecule.clone()),
            Ok(Solution::Contradictory(expected.clone()))
        );
        assert_eq!(
            resolver.resolve_into(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
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
            Resolver::new(&chemistry_model).resolve(molecule.clone()),
            Ok(Solution::Contradictory(expected.clone()))
        );
        assert_eq!(
            Resolver::new(&chemistry_model).resolve_into(&mut molecule),
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
            Resolver::new(&chemistry_model).resolve(molecule.clone()),
            Ok(Solution::Determined(expected.clone()))
        );
        assert_eq!(
            Resolver::new(&chemistry_model).resolve_into(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::natural(
        mol_dsl_concrete!(r#"{:atoms ["C#h3#u1#s2" "C#i13#c-#h3#n1"]}"#),
        mol_dsl!(r#"{:atoms ["C#c0#h3#n0#u1#s2" "C#i13#c-#h3#n1#u0#s"]}"#))]
    #[case::all_stages(
        mol_dsl_concrete!(r#"{
            :atoms ["C#i13#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [6 7 "1"] [6 8 "1"] [6 9 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
            :stereo-atoms [{:site 6 :ligands [7 8 9 [:h 6]] :attrs "Th0"}]}"#),
        mol_dsl!(r#"{
            :atoms ["C#i13#c0#h#n0#u0#s#a1" "C#c0#h#n0#u0#s#a1"
                    "C#c0#h#n0#u0#s#a1" "C#c0#h#n0#u0#s#a1"
                    "C#c0#h#n0#u0#s#a1" "C#c0#h#n0#u0#s#a1"
                    "C#c0#h#n0#u0#s#T0" "F#c0#h0#n0#u0#s"
                    "Cl#c0#h0#n0#u0#s" "Br#c0#h0#n0#u0#s"]
            :bonds [[0 1 "1#c0#u0#s#a"] [1 2 "1#c0#u0#s#a"] [2 3 "1#c0#u0#s#a"]
                    [3 4 "1#c0#u0#s#a"] [4 5 "1#c0#u0#s#a"] [5 0 "1#c0#u0#s#a"]
                    [6 7 "1#c0#u0#s"] [6 8 "1#c0#u0#s"] [6 9 "1#c0#u0#s"]]}"#),
    )]
    fn test_resolver_project(
        #[case] mut molecule: Molecule,
        #[case] expected: Molecule,
        #[values(false, true)] consuming: bool,
        #[values(false, true)] shared: bool,
    ) {
        let original = Molecule::combine_all([&molecule]);
        let retained = shared.then(|| molecule.clone());
        let model = ChemistryModel::default();
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                isotope: IsotopePolicy::Natural,
                ..Default::default()
            },
        );
        if consuming {
            assert_eq!(
                resolver.project(molecule, ProjectFlags::all()),
                Ok(Solution::Determined(expected))
            );
        } else {
            assert_eq!(
                resolver.project_into(&mut molecule, ProjectFlags::all()),
                Ok(Solution::Determined(()))
            );
            assert_eq!(molecule, expected);
        }
        if let Some(retained) = retained {
            assert_eq!(retained, original);
        }
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
    fn test_resolver_project_into_valence(
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
            .project_into(&mut molecule, ProjectFlags::all()),
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
    fn test_resolver_project_into_stereo(
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
            Resolver::with_config(&model, config).project_into(&mut molecule, ProjectFlags::all()),
            Ok(Solution::Determined(()))
        );
        assert_eq!(
            molecule.atom(AtomId(1)).implicit_hydrogens(),
            &NumForm::Lit(hydrogens)
        );
        assert_eq!(
            molecule.atom(AtomId(1)).constraints().tetrahedral_stereo(),
            &TetrahedralStereoForm::Stereo(StereoCoset::Lit(0))
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
            Resolver::new(&model).project(molecule.clone(), ProjectFlags::all()),
            Ok(Solution::Determined(original.clone()))
        );
        assert_eq!(
            Resolver::new(&model).project_into(&mut molecule, ProjectFlags::all()),
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
        #[values(ProjectFlags::all(), ProjectFlags::empty())] flags: ProjectFlags,
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
            resolver.project(molecule.clone(), flags),
            Err(expected.clone())
        );
        assert_eq!(resolver.project_into(&mut molecule, flags), Err(expected));
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::stereo(StereoConfigurationForm::Undetermined, NumForm::Lit(0), UnpairedElectronsForm::closed_shell(), IsotopeMassForm::Natural,
        ProjectError::Stereo(StereoProjectError::UnsupportedStereoAtom { stereo_atom: StereoAtomId(0) }))]
    #[case::aromatic_charge(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0), NumForm::Lit(1), UnpairedElectronsForm::closed_shell(), IsotopeMassForm::Natural,
        ProjectError::Aromaticity(AromaticityProjectError::ChargedAromaticSystem { aromatic_system: AromaticSystemId(0) }))]
    #[case::aromatic_spin(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0), NumForm::Lit(0), UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Lit(2) }, IsotopeMassForm::Natural,
        ProjectError::Aromaticity(AromaticityProjectError::AromaticSystemSpin { aromatic_system: AromaticSystemId(0) }))]
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
        #[values((false, false), (false, true), (true, false), (true, true))] execution: (
            bool,
            bool,
        ),
    ) {
        let (consuming, shared) = execution;
        let source = mol_dsl_concrete!(
            r#"{:atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"] [6 7 "1"] [6 8 "1"] [6 9 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
            :stereo-atoms [{:site 6 :ligands [7 8 9 [:h 6]] :attrs "Th0"}]}"#
        );
        let mut editor = source.edit();
        editor
            .stereo_atom_mut(StereoAtomId(0))
            .attributes_mut()
            .configuration = stereo;
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .attributes_mut()
            .charge = charge;
        editor
            .aromatic_system_mut(AromaticSystemId(0))
            .attributes_mut()
            .unpaired_electrons = spin;
        editor.atom_mut(AtomId(0)).attributes_mut().isotope_mass = isotope;
        editor.constraints_mut().push(Constraint::And(vec![
            Constraint::AromaticSystem(
                AromaticSystemId(0),
                AromaticSystemConstraintForm::electron_count(6_i64),
            ),
            Constraint::StereoAtom(
                StereoAtomId(0),
                StereoKind::Tetrahedral,
                StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                    Stereogenicity::Stereogenic,
                )),
            ),
        ]));
        let mut molecule = editor.finish().unwrap();
        let original = Molecule::combine_all([&molecule]);
        let retained = shared.then(|| molecule.clone());
        let model = ChemistryModel::default();
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                isotope: IsotopePolicy::Natural,
                ..Default::default()
            },
        );
        if consuming {
            assert_eq!(resolver.project(molecule, flags), Err(expected));
        } else {
            assert_eq!(resolver.project_into(&mut molecule, flags), Err(expected));
            assert_eq!(molecule, original);
        }
        if let Some(retained) = retained {
            assert_eq!(retained, original);
        }
    }

    #[rstest]
    #[case::literal(Some(NumForm::Lit(0)))]
    #[case::set(Some(NumForm::LitSet(Box::new(BTreeSet::from([0, 1])))))]
    #[case::vacuous(Some(NumForm::Undetermined))]
    #[case::absent(None)]
    fn test_resolver_plan_constitution(
        chemistry_model: ChemistryModel,
        #[case] valence: Option<NumForm>,
    ) {
        let mut atom = atom_dsl!("C#i=#c0#h4#u0#s#a!#T*");
        if let Some(valence) = valence {
            atom.constraints.set(AtomConstraintForm::Valence(valence));
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom, atom_dsl!("C#i=#c0#h4#n0#u0#s#v0#a!")],
            ..MoleculeEntries::default()
        });
        let state = ResolveState {
            completions: [(AtomId(0), smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s#v0#a!")])]
                .into_iter()
                .collect(),
            ..ResolveState::default()
        };
        let edits = Resolver::new(&chemistry_model).plan_constitution(&molecule, &state);
        let mut editor = molecule.edit();
        editor.apply(edits).unwrap();
        assert_eq!(
            editor.finish().unwrap(),
            mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#a!#T*" "C#i=#c0#h4#n0#u0#s#v0#a!"]}"#)
        );
    }

    #[rstest]
    #[case::undetermined(NumForm::Undetermined, NumForm::Lit(0))]
    #[case::singleton(
        NumForm::LitSet(Box::new(BTreeSet::from([0]))),
        NumForm::LitSet(Box::new(BTreeSet::from([0]))),
    )]
    fn test_resolver_plan_constitution_fields(
        chemistry_model: ChemistryModel,
        #[case] count: NumForm,
        #[case] expected_count: NumForm,
        #[values(RollbackMethod::Clone, RollbackMethod::Journal)] rollback: RollbackMethod,
    ) {
        let mut carbon = atom_dsl!("C#i=#v0#a!#T*");
        carbon.unpaired_electrons.count = count.clone();
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![carbon, atom_dsl!("N#i=#v0#a!")],
            ..MoleculeEntries::default()
        });
        let state = ResolveState {
            completions: [
                (AtomId(0), smallvec![atom_dsl!("C#i=#c0#h4#n0#u0#s")]),
                (AtomId(1), smallvec![atom_dsl!("N#i=#c0#h3#n1#u0#s")]),
            ]
            .into_iter()
            .collect(),
            ..ResolveState::default()
        };
        let edits = Resolver::new(&chemistry_model).plan_constitution(&molecule, &state);
        assert_eq!(
            edits,
            Edits::from_iter([
                Edit::ModifyAtomCharges {
                    changes: vec![
                        (
                            AtomHandle::Id(AtomId(0)),
                            NumForm::Undetermined,
                            NumForm::Lit(0)
                        ),
                        (
                            AtomHandle::Id(AtomId(1)),
                            NumForm::Undetermined,
                            NumForm::Lit(0)
                        ),
                    ]
                },
                Edit::ModifyAtomImplicitHydrogens {
                    changes: vec![
                        (
                            AtomHandle::Id(AtomId(0)),
                            NumForm::Undetermined,
                            NumForm::Lit(4)
                        ),
                        (
                            AtomHandle::Id(AtomId(1)),
                            NumForm::Undetermined,
                            NumForm::Lit(3)
                        ),
                    ]
                },
                Edit::ModifyAtomLonePairs {
                    changes: vec![
                        (
                            AtomHandle::Id(AtomId(0)),
                            NumForm::Undetermined,
                            NumForm::Lit(0)
                        ),
                        (
                            AtomHandle::Id(AtomId(1)),
                            NumForm::Undetermined,
                            NumForm::Lit(1)
                        ),
                    ]
                },
                Edit::ModifyAtomUnpairedElectrons {
                    changes: vec![
                        (
                            AtomHandle::Id(AtomId(0)),
                            UnpairedElectronsForm {
                                count,
                                multiplicity: NumForm::Undetermined,
                            },
                            UnpairedElectronsForm {
                                count: expected_count.clone(),
                                multiplicity: NumForm::Lit(1),
                            }
                        ),
                        (
                            AtomHandle::Id(AtomId(1)),
                            UnpairedElectronsForm::default(),
                            UnpairedElectronsForm::closed_shell()
                        ),
                    ]
                },
                Edit::ModifyAtomConstraints {
                    key: AtomConstraintKey::Valence,
                    changes: vec![
                        (
                            AtomHandle::Id(AtomId(0)),
                            Some(AtomConstraintForm::Valence(NumForm::Lit(0))),
                            None
                        ),
                        (
                            AtomHandle::Id(AtomId(1)),
                            Some(AtomConstraintForm::Valence(NumForm::Lit(0))),
                            None
                        ),
                    ],
                },
            ])
        );
        let mut carbon = atom_dsl!("C#i=#c0#h4#n0#u0#s#a!#T*");
        carbon.unpaired_electrons.count = expected_count;
        let expected = Molecule::from_entries(MoleculeEntries {
            atoms: vec![carbon, atom_dsl!("N#i=#c0#h3#n1#u0#s#a!")],
            ..MoleculeEntries::default()
        });
        let correspondence = MoleculeCorrespondence::new(
            Correspondence::identity(2),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
            Correspondence::empty(),
        );
        assert_eq!(molecule.clone().apply(edits.clone()), Ok(expected.clone()));
        assert_eq!(
            molecule.clone().tracked_apply(edits.clone()),
            Ok((expected.clone(), correspondence.clone()))
        );
        assert_eq!(
            molecule.tracked_transact([edits], rollback),
            Ok(correspondence)
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::empty(Molecule::new())]
    #[case::complete(mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#a!"]}"#))]
    #[case::uncompleted(mol_dsl!(r#"{:atoms ["C#i=#v0#a!"]}"#))]
    fn test_resolver_plan_constitution_identity(
        chemistry_model: ChemistryModel,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            Resolver::new(&chemistry_model).plan_constitution(&molecule, &ResolveState::default()),
            Edits::new()
        );
    }
}
