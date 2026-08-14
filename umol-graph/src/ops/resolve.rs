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
use std::collections::BTreeSet;

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
    AromaticSystemForm, AtomHandle, AtomId, Edits, Molecule, Transaction, TransactionError,
};
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;
pub use valence::{ValenceContradiction, ValenceError, ValenceResolver};

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
use crate::ops::model::{ChemistryModel, ValenceTieBreak};
use crate::ops::valence::compare::compare_by_key;
use crate::ops::valence::{AtomCompletions, ResolveReport};

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
    #[error("constitution commit failed: {0}")]
    Commit(TransactionError),
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

/// Resolution left the molecule underdetermined (no contradiction, but not ground).
/// Surfaced as an error only at boundaries that require a determined result.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
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
            valence: ValenceResolver::new(&model.valence),
            aromaticity: AromaticityResolver::with_config(&model.aromaticity, config.aromaticity),
            stereo: StereoResolver::with_config(&model.stereo, config.stereo),
            bonds: BondsResolver::new(),
            multicenter_bonds: MulticenterBondsResolver::new(),
            tie_break: model.valence.tie_break,
        }
    }

    pub fn resolve(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<ResolveReport, ResolveContradiction>, ResolveError> {
        let state = match self
            .valence
            .admit(molecule)
            .map_err(ResolveError::Valence)?
        {
            Solution::Determined(state) => state,
            Solution::Underdetermined(_) => {
                return Ok(Solution::Underdetermined(ResolveReport::default()));
            }
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction.into()));
            }
        };
        let outcome = self
            .aromaticity
            .select(molecule, state, self.tie_break)
            .map_err(ResolveError::Aromaticity)?;
        let mut state = match outcome {
            Solution::Determined(state) => state,
            Solution::Underdetermined(state) => {
                return Ok(Solution::Underdetermined(state.to_report()));
            }
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction.into()));
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
        let mut editor = molecule.edit();
        let mut journal = Transaction::default();
        let mut edits = Edits::new();
        for (atom, disjuncts) in state.completions.iter() {
            let current = molecule.atom(atom).attributes;
            // The constraint channel holds assertions only: the commit
            // narrows fields; candidate constraints stay solver state.
            let mut selected = disjuncts[0].clone();
            selected.constraints = current.constraints.clone();
            let update = current.difference_to(&selected);
            edits.update_atom(AtomHandle::Id(atom), current, &update);
        }
        let existing: BTreeSet<Vec<AtomId>> = molecule
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
                .plan_system(molecule, atoms.clone(), system.clone())
            {
                edits.push(edit);
            }
        }
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => return Err(ResolveError::Commit(transaction)),
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

        let resolved = editor.build();
        if resolved.is_ground() {
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
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::{fixture, rstest};
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{AtomConstraintForm, AtomId, EntityKind, MulticenterValenceForm};
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_ground};

    use super::*;
    use crate::ops::aromaticity::{AromaticityError, AromaticityInconsistency};
    use crate::ops::model::{
        AromaticityModel, AromaticityRule, ChemistryModel, ElementScope, RingLimits, StereoModel,
        ValenceModel,
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
            },
            stereo: StereoModel::default(),
        }
    }

    #[fixture]
    fn aromatic_molecule() -> Molecule {
        mol_dsl_ground!(
            r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#
        )
    }

    #[fixture]
    fn stereo_molecule() -> Molecule {
        mol_dsl_ground!(
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
            },
            stereo: StereoModel::default(),
        };
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#v0#a!"]}"#);
        assert_eq!(
            Resolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Determined(ResolveReport::default()))
        );
        assert_eq!(
            molecule,
            mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#a!"]}"#)
        );
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
        mol_dsl_ground!(r#"{
            :atoms ["C#h#v2#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::stereo(
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h*#n0#u0#s#T1" "F#i=#c0#h0#n0#u0#s"
                    "Cl#i=#c0#h0#n0#u0#s" "Br#i=#c0#h0#n0#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"]]
        }"#),
        mol_dsl_ground!(r#"{
            :atoms ["C#h#v3#a!#T1" "F" "Cl" "Br"]
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
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 } },
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
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar { ring_limits: RingLimits::default() } },
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
        AromaticityModel { scope: ElementScope::AllowList(vec![Element::C]), rule: AromaticityRule::Hueckel { ring_limits: RingLimits::default() } },
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
        let mut molecule = mol_dsl_ground!(
            r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
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
