//! Composite resolver: chains the per-entity resolvers (valence,
//! aromaticity, stereo, bonds, multicenter bonds) on a single `MoleculeAst`.
//!
//! `Determined` requires every entity (atoms, bonds, dative bonds, aromatic
//! systems, multicenter bonds, noncovalent bonds) to be ground.

pub mod aromaticity;
pub mod bonds;
pub mod multicenter;
pub mod stereo;
pub mod valence;

use std::any::Any;

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
use umol_ast::ast::{MoleculeAst, Transaction, TransactionError};
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;
pub use valence::{ValenceContradiction, ValenceError, ValenceResolver};

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
use crate::ops::model::ChemistryModel;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ResolveConfig {
    pub aromaticity: AromaticityResolveConfig,
    pub stereo: StereoResolveConfig,
}

#[derive(Clone, Debug)]
pub struct Resolver<'a> {
    pub valence: ValenceResolver<'a>,
    pub aromaticity: AromaticityResolver,
    pub stereo: StereoResolver,
    pub bonds: BondsResolver,
    pub multicenter_bonds: MulticenterBondsResolver,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ResolverContradiction {
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
pub enum ResolverError {
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
    #[error("rollback failed after {cause}: {rollback}")]
    RollbackFailed {
        cause: ResolverRollbackCause,
        rollback: TransactionError,
    },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ResolverRollbackCause {
    #[error("resolver contradiction: {0}")]
    Contradiction(ResolverContradiction),
    #[error("resolver underdetermined")]
    Underdetermined,
    #[error("resolver error: {0}")]
    Error(Box<ResolverError>),
}

impl UmolError for ResolverError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl UmolError for ResolverContradiction {
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
        }
    }

    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), ResolverContradiction>, ResolverError> {
        let mut editor = ast.edit();
        let mut journal = Transaction::default();

        let edits = match self.valence.plan(ast) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolverContradiction::Valence(contradiction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Contradiction(contradiction),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error = ResolverError::Valence(ValenceError::Transaction(transaction));
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Err(error);
            }
        }
        let mut state = editor.build();
        editor = state.edit();

        let outcome = match self.aromaticity.plan(&state) {
            Ok(outcome) => outcome,
            Err(error) => {
                let error = ResolverError::Aromaticity(error);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Err(error);
            }
        };
        let edits = match outcome {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Underdetermined,
                        rollback,
                    });
                }
                *ast = editor.build();
                return Ok(Solution::Underdetermined(()));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolverContradiction::Aromaticity(contradiction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Contradiction(contradiction),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error = ResolverError::Aromaticity(AromaticityError::Transaction(transaction));
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Err(error);
            }
        }
        state = editor.build();
        editor = state.edit();

        let outcome = match self.stereo.plan(&state) {
            Ok(outcome) => outcome,
            Err(error) => {
                let error = ResolverError::Stereo(error);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Err(error);
            }
        };
        let edits = match outcome {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Underdetermined,
                        rollback,
                    });
                }
                *ast = editor.build();
                return Ok(Solution::Underdetermined(()));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolverContradiction::Stereo(contradiction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Contradiction(contradiction),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error = ResolverError::Stereo(StereoError::Transaction(transaction));
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Err(error);
            }
        }
        state = editor.build();
        editor = state.edit();

        let edits = self.bonds.plan(&state);
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error = ResolverError::Bonds(BondsError::Transaction(transaction));
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Err(error);
            }
        }
        state = editor.build();
        editor = state.edit();

        let edits = match self.multicenter_bonds.plan(&state) {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => {
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Underdetermined,
                        rollback,
                    });
                }
                *ast = editor.build();
                return Ok(Solution::Underdetermined(()));
            }
            Solution::Contradictory(contradiction) => {
                let contradiction = ResolverContradiction::MulticenterBonds(contradiction);
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Contradiction(contradiction),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        match editor.transact(edits) {
            Ok(transaction) => journal.append(transaction),
            Err(transaction) => {
                let error = ResolverError::MulticenterBonds(MulticenterBondsError::Transaction(
                    transaction,
                ));
                if let Err(rollback) = journal.rollback(&mut editor) {
                    return Err(ResolverError::RollbackFailed {
                        cause: ResolverRollbackCause::Error(Box::new(error)),
                        rollback,
                    });
                }
                *ast = editor.build();
                return Err(error);
            }
        }

        let resolved = editor.build();
        let solution = if resolved.is_ground() {
            Solution::Determined(())
        } else {
            Solution::Underdetermined(())
        };
        *ast = resolved;
        Ok(solution)
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::{fixture, rstest};
    use umol_ast::ast::{
        AtomConstraintAst, AtomId, IncidenceConstraintContradiction, MulticenterValenceAst,
    };
    use umol_ast::{atom_dsl, mol_dsl, mol_dsl_ground};
    use umol_chem::element::Element;

    use super::*;
    use crate::ops::aromaticity::AromaticityInconsistency;
    use crate::ops::model::{
        AromaticityModel, ChemistryModel, ElementScope, RingLimits, StereoModel, ValenceModel,
    };
    use crate::ops::stereo::StereoInconsistency;
    use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

    #[fixture]
    fn chemistry_model() -> ChemistryModel {
        ChemistryModel {
            valence: ValenceModel::Counts {
                table: Cow::Borrowed(ValenceTable::default_table()),
            },
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
            stereo: StereoModel::default(),
        }
    }

    #[fixture]
    fn aromatic_molecule() -> MoleculeAst {
        mol_dsl_ground!(
            r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#
        )
    }

    #[fixture]
    fn stereo_molecule() -> MoleculeAst {
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
        ResolverError::RollbackFailed {
            cause: ResolverRollbackCause::Contradiction(
                ResolverContradiction::Stereo(StereoContradiction::Inconsistency(
                    StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) },
                )),
            ),
            rollback: TransactionError::OldStateMismatch,
        },
        "rollback failed after resolver contradiction: stereo inconsistency: tetrahedral stereo constraint at atom AtomId(1) cannot be realized: precondition failed: old state does not match current"
    )]
    #[case::execution(
        ResolverError::RollbackFailed {
            cause: ResolverRollbackCause::Error(Box::new(ResolverError::Bonds(
                BondsError::Transaction(TransactionError::IdOutOfRange("bond")),
            ))),
            rollback: TransactionError::OldStateMismatch,
        },
        "rollback failed after resolver error: id out of range: bond: precondition failed: old state does not match current"
    )]
    #[case::underdetermined(
        ResolverError::RollbackFailed {
            cause: ResolverRollbackCause::Underdetermined,
            rollback: TransactionError::OldStateMismatch,
        },
        "rollback failed after resolver underdetermined: precondition failed: old state does not match current"
    )]
    fn test_resolver_error(#[case] error: ResolverError, #[case] expected: &str) {
        assert_eq!(error.to_string(), expected);
    }

    #[rstest]
    fn test_resolver_new(
        chemistry_model: ChemistryModel,
        aromatic_molecule: MoleculeAst,
        stereo_molecule: MoleculeAst,
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
        aromatic_molecule: MoleculeAst,
        stereo_molecule: MoleculeAst,
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
    #[case::counts(ValenceModel::Counts {
        table: Cow::Borrowed(ValenceTable::default_table()),
    })]
    #[case::atom_typing(ValenceModel::AtomTyping {
        registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        )])),
    })]
    fn test_resolver_resolve(#[case] valence: ValenceModel) {
        let model = ChemistryModel {
            valence,
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
            stereo: StereoModel::default(),
        };
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0#h4#n0#u0#s#v0#a!"]}"#);
        assert_eq!(
            Resolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Determined(()))
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
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
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
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :type "Th1"}]
        }"#)
    )]
    fn test_resolver_resolve_stages(
        chemistry_model: ChemistryModel,
        #[case] mut molecule: MoleculeAst,
        #[case] expected: MoleculeAst,
    ) {
        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    fn test_resolver_resolve_underdetermined(chemistry_model: ChemistryModel) {
        let mut molecule = mol_dsl!(
            r#"{
            :atoms ["C#i*#c0#h4#n0#u0#s#v0#a!" "C#i=#c0#h4#n0#u0#s"]
            :noncovalent-bonds [{:atoms [0 1] :type "*"}]
        }"#
        );
        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(
            molecule,
            mol_dsl!(
                r#"{
                :atoms ["C#i=#c0#h4#n0#u0#s#v0#a!" "C#i=#c0#h4#n0#u0#s"]
                :noncovalent-bonds [{:atoms [0 1] :type "*"}]
            }"#
            )
        );
    }

    #[rstest]
    #[case::counts(ValenceModel::Counts {
        table: Cow::Borrowed(ValenceTable::default_table()),
    })]
    #[case::atom_typing(ValenceModel::AtomTyping {
        registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
    })]
    fn test_resolver_resolve_partial(#[case] valence: ValenceModel) {
        let model = ChemistryModel {
            valence,
            ..ChemistryModel::default()
        };
        let mut molecule = mol_dsl!(r#"{:atoms ["C#c0" "{C,N}#c0"]}"#);
        let original = molecule.clone();

        assert_eq!(
            Resolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Underdetermined(()))
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
            Ok(Solution::Underdetermined(()))
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::underdetermined(
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h0#n0#u0#s#v0#a!#m1"
                    "C#i=#c0#h0#n0#u0#s#v0#a!#m1"
                    "C#i=#c0#h0#n0#u0#s#v0#a!#m1"]
            :multicenter-bonds [{:atoms [0 1 2] :type "*"}]
        }"#),
        Solution::Underdetermined(())
    )]
    #[case::contradiction(
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h0#n0#u0#s#v0#a!#m1"]
        }"#),
        Solution::Contradictory(ResolverContradiction::MulticenterBonds(
            MulticenterBondsContradiction::Constraint(
                IncidenceConstraintContradiction::Atom {
                    atom: AtomId(0),
                    constraint: AtomConstraintAst::multicenter_valence(
                        MulticenterValenceAst::multicenter(1),
                    ),
                },
            ),
        ))
    )]
    fn test_resolver_resolve_multicenter_constraint_precondition(
        #[case] mut molecule: MoleculeAst,
        #[case] expected: Solution<(), ResolverContradiction>,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel::AtomTyping {
                registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
                    "C#i=#c0#h0#n0#u0#s#v0#a!#m1"
                )])),
            },
            ..ChemistryModel::default()
        };
        let original = molecule.clone();

        assert_eq!(Resolver::new(&model).resolve(&mut molecule), Ok(expected));
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::aromaticity(
        AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        },
        mol_dsl!(r#"{
            :atoms ["N#i*#c0#h#n0#u0#s#a2" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a"]
            :bonds [[0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"]
                    [3 4 "1#c0#u0#s"] [4 0 "1#c0#u0#s"]]
        }"#),
        ResolverContradiction::Aromaticity(AromaticityContradiction::ClarNonBenzenoid(
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
        ResolverContradiction::Aromaticity(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }
        ))
    )]
    #[case::stereo(
        AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        },
        mol_dsl!(r#"{
            :atoms ["C#i*#c0#h#n0#u0#s#a#T1" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"
                    "C#i=#c0#h#n0#u0#s#a" "C#i=#c0#h#n0#u0#s#a"]
            :bonds [[0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"]
                    [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"]
                    [4 5 "1#c0#u0#s"] [5 0 "1#c0#u0#s"]]
        }"#),
        ResolverContradiction::Stereo(StereoContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(0) }
        ))
    )]
    fn test_resolver_resolve_contradiction(
        mut chemistry_model: ChemistryModel,
        #[case] aromaticity: AromaticityModel,
        #[case] mut molecule: MoleculeAst,
        #[case] expected: ResolverContradiction,
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
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#
        );
        let expected = molecule.clone();
        assert_eq!(
            Resolver::new(&chemistry_model).resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }
}
