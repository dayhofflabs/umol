//! Tier-1 constraint validator: cross-entity and molecule-scope constraint evaluation. Run at AST
//! construction/raise and available standalone; never consults a chemistry model.

use thiserror::Error;
use umol_graph_core::{
    ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_utils::solution::Solution;

pub mod incidence;
pub mod molecule;
pub mod relational;
pub mod ring;

use incidence::{
    bond_components_by_atom, validate_aromatic_system_constraint, validate_atom_constraint,
    validate_bond_constraint, validate_dative_bond_constraint,
    validate_multicenter_bond_constraint, validate_noncovalent_bond_constraint,
};
pub use incidence::{IncidenceConstraintContradiction, IncidenceConstraintValidator};
pub use molecule::{MoleculeConstraintContradiction, MoleculeConstraintValidator};
pub use relational::{RelationalConstraintContradiction, RelationalConstraintValidator};
pub use ring::{RingConstraintContradiction, RingConstraintValidator};

use super::super::constraint::{
    AtomConstraintForm, BondConstraintForm, Constraint, DativeBondConstraintForm,
    MoleculeConstraint, StereoAtomConstraintForm, StereoBondConstraintForm,
};
use super::super::entity::Entity;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};
use super::super::molecule::MoleculeAst;
use super::super::ring::{RingConfig, RingModel};
use super::super::stereo::StereoKind;
use super::super::substructure::SubstructureMatchAlgorithm;
use super::super::traits::Lattice;
use super::super::value::NumForm;
use super::super::view::RingViews;

/// Algorithm selectors used by complete model-independent constraint validation.
///
/// Focused validators take only the selectors they require. This bundle has no
/// default at the AST layer so every algorithm choice remains explicit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConstraintValidateConfig {
    pub relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    pub connected_components_algorithm: ConnectedComponentsAlgorithm,
    pub substructure_match_algorithm: SubstructureMatchAlgorithm,
    pub subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
}

/// Cross-check between inline and molecule-scope constraints and their
/// model-independent values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConstraintValidator {
    config: ConstraintValidateConfig,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConstraintContradiction {
    #[error(transparent)]
    Incidence(#[from] IncidenceConstraintContradiction),
    #[error(transparent)]
    Ring(#[from] RingConstraintContradiction),
    #[error(transparent)]
    Relational(#[from] RelationalConstraintContradiction),
    #[error(transparent)]
    Molecule(#[from] MoleculeConstraintContradiction),
    #[error("logical constraint is not satisfied: {constraint:?}")]
    Logical { constraint: Constraint },
    #[error("stereo atom {id:?} of kind {kind:?} does not satisfy constraint {constraint:?}")]
    StereoAtom {
        id: StereoAtomId,
        kind: StereoKind,
        constraint: StereoAtomConstraintForm,
    },
    #[error("stereo bond {id:?} of kind {kind:?} does not satisfy constraint {constraint:?}")]
    StereoBond {
        id: StereoBondId,
        kind: StereoKind,
        constraint: StereoBondConstraintForm,
    },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConstraintError {
    #[error("constraint references unavailable {entity}")]
    InvalidReference { entity: Entity },
    /// Dative-bond ring topology is deferred to the coordination/haptic entity split in doc 117.
    #[error("ring membership is not defined for dative bond {bond:?}")]
    DativeBondRingMembershipUnsupported { bond: DativeBondId },
}

impl ConstraintValidator {
    pub fn new(config: ConstraintValidateConfig) -> Self {
        Self { config }
    }

    pub fn validate(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        let mut evaluation = ConstraintEvaluation::new(ast, self.config);
        let mut any_underdetermined = false;

        for id in ast.atoms().ids() {
            for constraint in ast.atom(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_atom(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in ast.bonds().ids() {
            for constraint in ast.bond(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_bond(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in ast.dative_bonds().ids() {
            for constraint in ast.dative_bond(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_dative_bond(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in ast.aromatic_systems().ids() {
            for constraint in ast.aromatic_system(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_aromatic_system(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in ast.multicenter_bonds().ids() {
            for constraint in ast.multicenter_bond(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_multicenter_bond(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in ast.noncovalent_bonds().ids() {
            for constraint in ast.noncovalent_bond(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_noncovalent_bond(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for constraint in ast.constraints().iter() {
            if let Some(contradiction) =
                observe(evaluation.evaluate(constraint)?, &mut any_underdetermined)
            {
                return Ok(Solution::Contradictory(contradiction));
            }
        }

        Ok(if any_underdetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }
}

struct ConstraintEvaluation<'a> {
    ast: &'a MoleculeAst,
    config: ConstraintValidateConfig,
    rings: Option<RingViews<'a>>,
    component_by_atom: Option<Vec<usize>>,
}

impl<'a> ConstraintEvaluation<'a> {
    fn new(ast: &'a MoleculeAst, config: ConstraintValidateConfig) -> Self {
        Self {
            ast,
            config,
            rings: None,
            component_by_atom: None,
        }
    }

    fn evaluate(
        &mut self,
        constraint: &Constraint,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        match constraint {
            Constraint::Atom(id, constraint) => self.evaluate_atom(*id, constraint),
            Constraint::Bond(id, constraint) => self.evaluate_bond(*id, constraint),
            Constraint::DativeBond(id, constraint) => self.evaluate_dative_bond(*id, constraint),
            Constraint::AromaticSystem(id, constraint) => {
                self.evaluate_aromatic_system(*id, constraint)
            }
            Constraint::MulticenterBond(id, constraint) => {
                self.evaluate_multicenter_bond(*id, constraint)
            }
            Constraint::NoncovalentBond(id, constraint) => {
                self.evaluate_noncovalent_bond(*id, constraint)
            }
            Constraint::StereoAtom(id, kind, constraint) => {
                self.evaluate_stereo_atom(*id, *kind, constraint)
            }
            Constraint::StereoBond(id, kind, constraint) => {
                self.evaluate_stereo_bond(*id, *kind, constraint)
            }
            Constraint::Relational(constraint) => RelationalConstraintValidator
                .validate(self.ast, constraint, self.config.relevant_cycle_algorithm)
                .map(|outcome| outcome.map_contradiction(ConstraintContradiction::from)),
            Constraint::Molecule(constraint) => self.evaluate_molecule(constraint),
            Constraint::And(constraints) => self.evaluate_and(constraints),
            Constraint::Or(constraints) => self.evaluate_or(constraint, constraints),
            Constraint::Not(inner) => match self.evaluate(inner)? {
                Solution::Determined(()) => {
                    Ok(Solution::Contradictory(ConstraintContradiction::Logical {
                        constraint: constraint.clone(),
                    }))
                }
                Solution::Underdetermined(()) => Ok(Solution::Underdetermined(())),
                Solution::Contradictory(_) => Ok(Solution::Determined(())),
            },
        }
    }

    fn evaluate_and(
        &mut self,
        constraints: &[Constraint],
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        let mut any_underdetermined = false;
        for constraint in constraints {
            if let Some(contradiction) =
                observe(self.evaluate(constraint)?, &mut any_underdetermined)
            {
                return Ok(Solution::Contradictory(contradiction));
            }
        }
        Ok(if any_underdetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }

    fn evaluate_or(
        &mut self,
        constraint: &Constraint,
        alternatives: &[Constraint],
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        if alternatives.is_empty() {
            return Ok(Solution::Determined(()));
        }
        let mut any_underdetermined = false;
        for alternative in alternatives {
            match self.evaluate(alternative)? {
                Solution::Determined(()) => return Ok(Solution::Determined(())),
                Solution::Underdetermined(()) => any_underdetermined = true,
                Solution::Contradictory(_) => {}
            }
        }
        Ok(if any_underdetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Contradictory(ConstraintContradiction::Logical {
                constraint: constraint.clone(),
            })
        })
    }

    fn evaluate_atom(
        &mut self,
        id: AtomId,
        constraint: &AtomConstraintForm,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        self.require(Entity::Atom(id))?;
        if is_atom_ring_constraint(constraint) {
            if constraint.is_undetermined() {
                return Ok(Solution::Determined(()));
            }
            let ring_atom = self.rings().atom(id);
            let (asserted, derived) = match constraint {
                AtomConstraintForm::RingDegree(asserted) => (asserted, ring_atom.ring_degree()),
                AtomConstraintForm::RingValence(asserted) => (asserted, ring_atom.ring_valence()),
                AtomConstraintForm::RingMembership(membership) => (
                    &membership.count,
                    ring_atom.ring_membership(membership.scope),
                ),
                _ => unreachable!("ring constraint classified above"),
            };
            Ok(evaluate_value(asserted, &derived, || {
                RingConstraintContradiction::Atom {
                    atom: id,
                    constraint: constraint.clone(),
                }
                .into()
            }))
        } else {
            Ok(validate_atom_constraint(self.ast, id, constraint)
                .map_contradiction(ConstraintContradiction::from))
        }
    }

    fn evaluate_bond(
        &mut self,
        id: BondId,
        constraint: &BondConstraintForm,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        self.require(Entity::Bond(id))?;
        if let BondConstraintForm::RingMembership(membership) = constraint {
            if membership.is_undetermined() {
                return Ok(Solution::Determined(()));
            }
            let derived = self.rings().bond(id).ring_membership(membership.scope);
            Ok(evaluate_value(&membership.count, &derived, || {
                RingConstraintContradiction::Bond {
                    bond: id,
                    constraint: constraint.clone(),
                }
                .into()
            }))
        } else {
            Ok(validate_bond_constraint(self.ast, id, constraint)
                .expect("non-ring bond constraint")
                .map_contradiction(ConstraintContradiction::from))
        }
    }

    fn evaluate_dative_bond(
        &mut self,
        id: DativeBondId,
        constraint: &DativeBondConstraintForm,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        self.require(Entity::DativeBond(id))?;
        if let DativeBondConstraintForm::RingMembership(membership) = constraint {
            if membership.is_undetermined() {
                Ok(Solution::Determined(()))
            } else {
                Err(ConstraintError::DativeBondRingMembershipUnsupported { bond: id })
            }
        } else {
            Ok(validate_dative_bond_constraint(self.ast, id, constraint)
                .expect("non-ring dative-bond constraint")
                .map_contradiction(ConstraintContradiction::from))
        }
    }

    fn evaluate_aromatic_system(
        &self,
        id: AromaticSystemId,
        constraint: &super::super::constraint::AromaticSystemConstraintForm,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        self.require(Entity::AromaticSystem(id))?;
        Ok(
            validate_aromatic_system_constraint(self.ast, id, constraint)
                .map_contradiction(ConstraintContradiction::from),
        )
    }

    fn evaluate_multicenter_bond(
        &self,
        id: MulticenterBondId,
        constraint: &super::super::constraint::MulticenterBondConstraintForm,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        self.require(Entity::MulticenterBond(id))?;
        Ok(
            validate_multicenter_bond_constraint(self.ast, id, constraint)
                .map_contradiction(ConstraintContradiction::from),
        )
    }

    fn evaluate_noncovalent_bond(
        &mut self,
        id: NoncovalentBondId,
        constraint: &super::super::constraint::NoncovalentBondConstraintForm,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        self.require(Entity::NoncovalentBond(id))?;
        let intramolecular = if constraint.is_undetermined() {
            false
        } else {
            let [a, b] = self.ast.noncovalent_bond(id).atom_ids();
            self.components()[a.index()] == self.components()[b.index()]
        };
        Ok(
            validate_noncovalent_bond_constraint(id, constraint, intramolecular)
                .map_contradiction(ConstraintContradiction::from),
        )
    }

    fn evaluate_stereo_atom(
        &self,
        id: StereoAtomId,
        kind: StereoKind,
        constraint: &StereoAtomConstraintForm,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        let stereo = self
            .ast
            .stereo_atoms()
            .get(id)
            .ok_or(ConstraintError::InvalidReference {
                entity: Entity::StereoAtom(id),
            })?;
        Ok(evaluate_stereo_constraint(
            stereo.kind() == kind,
            constraint,
            stereo.constraints().get(constraint.key()),
            || ConstraintContradiction::StereoAtom {
                id,
                kind,
                constraint: constraint.clone(),
            },
        ))
    }

    fn evaluate_stereo_bond(
        &self,
        id: StereoBondId,
        kind: StereoKind,
        constraint: &StereoBondConstraintForm,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        let stereo = self
            .ast
            .stereo_bonds()
            .get(id)
            .ok_or(ConstraintError::InvalidReference {
                entity: Entity::StereoBond(id),
            })?;
        Ok(evaluate_stereo_constraint(
            stereo.kind() == kind,
            constraint,
            stereo.constraints().get(constraint.key()),
            || ConstraintContradiction::StereoBond {
                id,
                kind,
                constraint: constraint.clone(),
            },
        ))
    }

    fn evaluate_molecule(
        &mut self,
        constraint: &MoleculeConstraint,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        if let MoleculeConstraint::Connected { atoms } = constraint {
            let atoms: Vec<_> = match atoms {
                Some(atoms) => {
                    for &id in atoms {
                        self.require(Entity::Atom(id))?;
                    }
                    atoms.clone()
                }
                None => self.ast.atoms().ids().collect(),
            };
            let determined = atoms.len() < 2
                || atoms
                    .iter()
                    .all(|id| self.components()[id.index()] == self.components()[atoms[0].index()]);
            return Ok(if determined {
                Solution::Determined(())
            } else {
                Solution::Contradictory(
                    MoleculeConstraintContradiction {
                        constraint: constraint.clone(),
                    }
                    .into(),
                )
            });
        }
        MoleculeConstraintValidator
            .validate(self.ast, constraint, self.config)
            .map(|outcome| outcome.map_contradiction(ConstraintContradiction::from))
    }

    fn rings(&mut self) -> &RingViews<'a> {
        if self.rings.is_none() {
            self.rings = Some(self.ast.rings(
                RingModel::default(),
                RingConfig {
                    relevant_cycle_algorithm: self.config.relevant_cycle_algorithm,
                    ..RingConfig::default()
                },
            ));
        }
        self.rings.as_ref().expect("ring views initialized")
    }

    fn components(&mut self) -> &[usize] {
        if self.component_by_atom.is_none() {
            self.component_by_atom = Some(bond_components_by_atom(
                self.ast,
                self.config.connected_components_algorithm,
            ));
        }
        self.component_by_atom
            .as_deref()
            .expect("components initialized")
    }

    fn require(&self, entity: Entity) -> Result<(), ConstraintError> {
        let present = match entity {
            Entity::Atom(id) => self.ast.atoms().contains(id),
            Entity::Bond(id) => self.ast.bonds().contains(id),
            Entity::DativeBond(id) => self.ast.dative_bonds().contains(id),
            Entity::AromaticSystem(id) => self.ast.aromatic_systems().contains(id),
            Entity::MulticenterBond(id) => self.ast.multicenter_bonds().contains(id),
            Entity::NoncovalentBond(id) => self.ast.noncovalent_bonds().contains(id),
            Entity::StereoAtom(id) => self.ast.stereo_atoms().contains(id),
            Entity::StereoBond(id) => self.ast.stereo_bonds().contains(id),
        };
        if present {
            Ok(())
        } else {
            Err(ConstraintError::InvalidReference { entity })
        }
    }
}

fn is_atom_ring_constraint(constraint: &AtomConstraintForm) -> bool {
    matches!(
        constraint,
        AtomConstraintForm::RingDegree(_)
            | AtomConstraintForm::RingValence(_)
            | AtomConstraintForm::RingMembership(_)
    )
}

fn evaluate_value(
    asserted: &NumForm,
    derived: &NumForm,
    contradiction: impl FnOnce() -> ConstraintContradiction,
) -> Solution<(), ConstraintContradiction> {
    if asserted.is_undetermined() {
        Solution::Determined(())
    } else if !derived.is_ground() {
        Solution::Underdetermined(())
    } else if asserted.matches(derived) {
        Solution::Determined(())
    } else {
        Solution::Contradictory(contradiction())
    }
}

fn evaluate_stereo_constraint<C>(
    kind_matches: bool,
    asserted: &C,
    stored: Option<&C>,
    contradiction: impl FnOnce() -> ConstraintContradiction,
) -> Solution<(), ConstraintContradiction>
where
    C: Lattice,
{
    if asserted.is_undetermined() {
        Solution::Determined(())
    } else if !kind_matches {
        Solution::Contradictory(contradiction())
    } else if stored.is_none_or(|stored| !stored.is_ground()) {
        Solution::Underdetermined(())
    } else if asserted.matches(stored.expect("ground stored constraint")) {
        Solution::Determined(())
    } else {
        Solution::Contradictory(contradiction())
    }
}

fn observe(
    outcome: Solution<(), ConstraintContradiction>,
    any_underdetermined: &mut bool,
) -> Option<ConstraintContradiction> {
    match outcome {
        Solution::Determined(()) => None,
        Solution::Underdetermined(()) => {
            *any_underdetermined = true;
            None
        }
        Solution::Contradictory(contradiction) => Some(contradiction),
    }
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_graph_core::{
        ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm,
        SubgraphIsomorphismAlgorithm,
    };

    use super::*;
    use crate::ir::constraint::{
        AtomConstraintForm, BondConstraintForm, MoleculeConstraint, RelationalConstraint, RingScope,
    };
    use crate::ir::id::{AtomId, BondId, DativeBondId};
    use crate::ir::substructure::SubstructureMatchAlgorithm;
    use crate::mol_dsl;

    const CONFIG: ConstraintValidateConfig = ConstraintValidateConfig {
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
        substructure_match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
        subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
    };

    #[fixture]
    fn molecule() -> MoleculeAst {
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#)
    }

    #[rstest]
    #[case::incidence(
        IncidenceConstraintContradiction::Atom {
            atom: AtomId(2),
            constraint: AtomConstraintForm::valence(4),
        },
        ConstraintContradiction::Incidence(IncidenceConstraintContradiction::Atom {
            atom: AtomId(2),
            constraint: AtomConstraintForm::valence(4),
        })
    )]
    #[case::ring(
        RingConstraintContradiction::Bond {
            bond: BondId(3),
            constraint: BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        },
        ConstraintContradiction::Ring(RingConstraintContradiction::Bond {
            bond: BondId(3),
            constraint: BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        })
    )]
    #[case::relational(
        RelationalConstraintContradiction {
            constraint: RelationalConstraint::DativeBondDonor {
                bond: DativeBondId(1),
                atom: AtomId(2),
            },
        },
        ConstraintContradiction::Relational(RelationalConstraintContradiction {
            constraint: RelationalConstraint::DativeBondDonor {
                bond: DativeBondId(1),
                atom: AtomId(2),
            },
        })
    )]
    #[case::molecule(
        MoleculeConstraintContradiction {
            constraint: MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(1), AtomId(2)]),
            },
        },
        ConstraintContradiction::Molecule(MoleculeConstraintContradiction {
            constraint: MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(1), AtomId(2)]),
            },
        })
    )]
    fn test_constraint_contradiction_from<T>(
        #[case] input: T,
        #[case] expected: ConstraintContradiction,
    ) where
        T: Into<ConstraintContradiction>,
    {
        assert_eq!(input.into(), expected);
    }

    #[rstest]
    #[case::invalid_reference(
        ConstraintError::InvalidReference {
            entity: Entity::Atom(AtomId(4)),
        },
        "constraint references unavailable atom 4"
    )]
    #[case::dative_ring(
        ConstraintError::DativeBondRingMembershipUnsupported {
            bond: DativeBondId(5),
        },
        "ring membership is not defined for dative bond DativeBondId(5)"
    )]
    fn test_constraint_error_display(#[case] input: ConstraintError, #[case] expected: &str) {
        assert_eq!(input.to_string(), expected);
    }

    #[rstest]
    #[case::and_determined(Constraint::And(vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(0)),
        Constraint::Atom(AtomId(0), AtomConstraintForm::degree(0)),
    ]), Solution::Determined(()))]
    #[case::and_underdetermined(Constraint::And(vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(0)),
        Constraint::Atom(AtomId(0), AtomConstraintForm::total_hydrogens(1)),
    ]), Solution::Underdetermined(()))]
    #[case::or_determined(Constraint::Or(vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(1)),
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(0)),
    ]), Solution::Determined(()))]
    #[case::or_underdetermined(Constraint::Or(vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(1)),
        Constraint::Atom(AtomId(0), AtomConstraintForm::total_hydrogens(1)),
    ]), Solution::Underdetermined(()))]
    #[case::not_contradictory(Constraint::Not(Box::new(
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(1)),
    )), Solution::Determined(()))]
    fn test_constraint_validator_logical_outcomes(
        mut molecule: MoleculeAst,
        #[case] constraint: Constraint,
        #[case] expected: Solution<(), ConstraintContradiction>,
    ) {
        molecule.constraints_mut().push(constraint);

        assert_eq!(
            ConstraintValidator::new(CONFIG).validate(&molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::or(Constraint::Or(vec![
        Constraint::Atom(AtomId(0), AtomConstraintForm::valence(1)),
        Constraint::Atom(AtomId(0), AtomConstraintForm::degree(1)),
    ]))]
    #[case::not(Constraint::Not(Box::new(Constraint::Atom(
        AtomId(0),
        AtomConstraintForm::valence(0)
    ),)))]
    fn test_constraint_validator_logical_contradiction(
        mut molecule: MoleculeAst,
        #[case] constraint: Constraint,
    ) {
        molecule.constraints_mut().push(constraint.clone());

        assert_eq!(
            ConstraintValidator::new(CONFIG).validate(&molecule),
            Ok(Solution::Contradictory(ConstraintContradiction::Logical {
                constraint
            }))
        );
    }

    #[rstest]
    #[case::determined(AtomConstraintForm::valence(0))]
    #[case::underdetermined(AtomConstraintForm::total_hydrogens(1))]
    #[case::contradictory(AtomConstraintForm::valence(1))]
    fn test_constraint_validator_inline_top_level_agreement(
        molecule: MoleculeAst,
        #[case] constraint: AtomConstraintForm,
    ) {
        let mut inline = molecule.clone().edit();
        inline
            .atom_mut(AtomId(0))
            .ast
            .constraints
            .set(constraint.clone());
        let inline = inline.build();
        let mut top_level = molecule.edit();
        top_level
            .constraints_mut()
            .push(Constraint::Atom(AtomId(0), constraint));
        let top_level = top_level.build();
        let validator = ConstraintValidator::new(CONFIG);

        assert_eq!(validator.validate(&inline), validator.validate(&top_level));
    }
}
