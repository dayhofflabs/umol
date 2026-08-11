//! Tier-2 constraint invariants: cross-entity and molecule-scope constraint evaluation without a
//! chemistry model.

use thiserror::Error;
use umol_graph_core::{ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm};
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
pub use incidence::{
    IncidenceConstraintInvariantsContradiction, IncidenceConstraintInvariantsValidator,
};
pub use molecule::{
    MoleculeConstraintInvariantsContradiction, MoleculeConstraintInvariantsValidator,
};
pub use relational::{
    RelationalConstraintInvariantsContradiction, RelationalConstraintInvariantsValidator,
};
pub use ring::{RingConstraintInvariantsContradiction, RingConstraintInvariantsValidator};
use umol_graph_ir::ir::{
    AromaticSystemConstraintForm, AromaticSystemId, AtomConstraintForm, AtomId, BondConstraintForm,
    BondId, Constraint, DativeBondConstraintForm, DativeBondId, Entity, Lattice, Molecule,
    MoleculeConstraint, MulticenterBondConstraintForm, MulticenterBondId,
    NoncovalentBondConstraintForm, NoncovalentBondId, NumForm, RingConfig, RingModel, RingViews,
    StereoAtomConstraintForm, StereoAtomId, StereoBondConstraintForm, StereoBondId, StereoKind,
};

/// Algorithm selectors used by complete model-independent constraint validation.
///
/// Focused validators take only the selectors they require. This bundle has no
/// default at the graph-IR layer so every algorithm choice remains explicit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConstraintValidateConfig {
    pub relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    pub connected_components_algorithm: ConnectedComponentsAlgorithm,
}

/// Cross-check between inline and molecule-scope constraints and their
/// model-independent values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConstraintInvariantsValidator {
    config: ConstraintValidateConfig,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConstraintInvariantsContradiction {
    #[error(transparent)]
    Incidence(#[from] IncidenceConstraintInvariantsContradiction),
    #[error(transparent)]
    Ring(#[from] RingConstraintInvariantsContradiction),
    #[error(transparent)]
    Relational(#[from] RelationalConstraintInvariantsContradiction),
    #[error(transparent)]
    Molecule(#[from] MoleculeConstraintInvariantsContradiction),
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
pub enum ConstraintInvariantsError {
    #[error("constraint references unavailable {entity}")]
    InvalidReference { entity: Entity },
    /// Dative-bond ring topology is deferred to the coordination/haptic entity split in doc 117.
    #[error("ring membership is not defined for dative bond {bond:?}")]
    DativeBondRingMembershipUnsupported { bond: DativeBondId },
}

impl ConstraintInvariantsValidator {
    pub fn new(config: ConstraintValidateConfig) -> Self {
        Self { config }
    }

    pub fn validate(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
        let mut evaluation = ConstraintEvaluation::new(molecule, self.config);
        let mut any_underdetermined = false;

        for id in molecule.atoms().ids() {
            for constraint in molecule.atom(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_atom(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in molecule.bonds().ids() {
            for constraint in molecule.bond(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_bond(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in molecule.dative_bonds().ids() {
            for constraint in molecule.dative_bond(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_dative_bond(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in molecule.aromatic_systems().ids() {
            for constraint in molecule.aromatic_system(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_aromatic_system(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in molecule.multicenter_bonds().ids() {
            for constraint in molecule.multicenter_bond(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_multicenter_bond(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for id in molecule.noncovalent_bonds().ids() {
            for constraint in molecule.noncovalent_bond(id).constraints().iter() {
                if let Some(contradiction) = observe(
                    evaluation.evaluate_noncovalent_bond(id, constraint)?,
                    &mut any_underdetermined,
                ) {
                    return Ok(Solution::Contradictory(contradiction));
                }
            }
        }
        for constraint in molecule.constraints().iter() {
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
    molecule: &'a Molecule,
    config: ConstraintValidateConfig,
    rings: Option<RingViews<'a>>,
    component_by_atom: Option<Vec<usize>>,
}

impl<'a> ConstraintEvaluation<'a> {
    fn new(molecule: &'a Molecule, config: ConstraintValidateConfig) -> Self {
        Self {
            molecule,
            config,
            rings: None,
            component_by_atom: None,
        }
    }

    fn evaluate(
        &mut self,
        constraint: &Constraint,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
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
            Constraint::Relational(constraint) => RelationalConstraintInvariantsValidator
                .validate(
                    self.molecule,
                    constraint,
                    self.config.relevant_cycle_algorithm,
                )
                .map(|outcome| outcome.map_contradiction(ConstraintInvariantsContradiction::from)),
            Constraint::Molecule(constraint) => self.evaluate_molecule(constraint),
            Constraint::And(constraints) => self.evaluate_and(constraints),
            Constraint::Or(constraints) => self.evaluate_or(constraint, constraints),
            Constraint::Not(inner) => match self.evaluate(inner)? {
                Solution::Determined(()) => Ok(Solution::Contradictory(
                    ConstraintInvariantsContradiction::Logical {
                        constraint: constraint.clone(),
                    },
                )),
                Solution::Underdetermined(()) => Ok(Solution::Underdetermined(())),
                Solution::Contradictory(_) => Ok(Solution::Determined(())),
            },
        }
    }

    fn evaluate_and(
        &mut self,
        constraints: &[Constraint],
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
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
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
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
            Solution::Contradictory(ConstraintInvariantsContradiction::Logical {
                constraint: constraint.clone(),
            })
        })
    }

    fn evaluate_atom(
        &mut self,
        id: AtomId,
        constraint: &AtomConstraintForm,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
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
                RingConstraintInvariantsContradiction::Atom {
                    atom: id,
                    constraint: constraint.clone(),
                }
                .into()
            }))
        } else {
            Ok(validate_atom_constraint(self.molecule, id, constraint)
                .map_contradiction(ConstraintInvariantsContradiction::from))
        }
    }

    fn evaluate_bond(
        &mut self,
        id: BondId,
        constraint: &BondConstraintForm,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
        self.require(Entity::Bond(id))?;
        if let BondConstraintForm::RingMembership(membership) = constraint {
            if membership.is_undetermined() {
                return Ok(Solution::Determined(()));
            }
            let derived = self.rings().bond(id).ring_membership(membership.scope);
            Ok(evaluate_value(&membership.count, &derived, || {
                RingConstraintInvariantsContradiction::Bond {
                    bond: id,
                    constraint: constraint.clone(),
                }
                .into()
            }))
        } else {
            Ok(validate_bond_constraint(self.molecule, id, constraint)
                .expect("non-ring bond constraint")
                .map_contradiction(ConstraintInvariantsContradiction::from))
        }
    }

    fn evaluate_dative_bond(
        &mut self,
        id: DativeBondId,
        constraint: &DativeBondConstraintForm,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
        self.require(Entity::DativeBond(id))?;
        if let DativeBondConstraintForm::RingMembership(membership) = constraint {
            if membership.is_undetermined() {
                Ok(Solution::Determined(()))
            } else {
                Err(ConstraintInvariantsError::DativeBondRingMembershipUnsupported { bond: id })
            }
        } else {
            Ok(
                validate_dative_bond_constraint(self.molecule, id, constraint)
                    .expect("non-ring dative-bond constraint")
                    .map_contradiction(ConstraintInvariantsContradiction::from),
            )
        }
    }

    fn evaluate_aromatic_system(
        &self,
        id: AromaticSystemId,
        constraint: &AromaticSystemConstraintForm,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
        self.require(Entity::AromaticSystem(id))?;
        Ok(
            validate_aromatic_system_constraint(self.molecule, id, constraint)
                .map_contradiction(ConstraintInvariantsContradiction::from),
        )
    }

    fn evaluate_multicenter_bond(
        &self,
        id: MulticenterBondId,
        constraint: &MulticenterBondConstraintForm,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
        self.require(Entity::MulticenterBond(id))?;
        Ok(
            validate_multicenter_bond_constraint(self.molecule, id, constraint)
                .map_contradiction(ConstraintInvariantsContradiction::from),
        )
    }

    fn evaluate_noncovalent_bond(
        &mut self,
        id: NoncovalentBondId,
        constraint: &NoncovalentBondConstraintForm,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
        self.require(Entity::NoncovalentBond(id))?;
        let intramolecular = if constraint.is_undetermined() {
            false
        } else {
            let [a, b] = self.molecule.noncovalent_bond(id).atom_ids();
            self.components()[a.index()] == self.components()[b.index()]
        };
        Ok(
            validate_noncovalent_bond_constraint(id, constraint, intramolecular)
                .map_contradiction(ConstraintInvariantsContradiction::from),
        )
    }

    fn evaluate_stereo_atom(
        &self,
        id: StereoAtomId,
        kind: StereoKind,
        constraint: &StereoAtomConstraintForm,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
        let stereo = self.molecule.stereo_atoms().get(id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::StereoAtom(id),
            },
        )?;
        Ok(evaluate_stereo_constraint(
            stereo.kind() == kind,
            constraint,
            stereo.constraints().get(constraint.key()),
            || ConstraintInvariantsContradiction::StereoAtom {
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
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
        let stereo = self.molecule.stereo_bonds().get(id).ok_or(
            ConstraintInvariantsError::InvalidReference {
                entity: Entity::StereoBond(id),
            },
        )?;
        Ok(evaluate_stereo_constraint(
            stereo.kind() == kind,
            constraint,
            stereo.constraints().get(constraint.key()),
            || ConstraintInvariantsContradiction::StereoBond {
                id,
                kind,
                constraint: constraint.clone(),
            },
        ))
    }

    fn evaluate_molecule(
        &mut self,
        constraint: &MoleculeConstraint,
    ) -> Result<Solution<(), ConstraintInvariantsContradiction>, ConstraintInvariantsError> {
        if let MoleculeConstraint::Connected { atoms } = constraint {
            let atoms: Vec<_> = match atoms {
                Some(atoms) => {
                    for &id in atoms {
                        self.require(Entity::Atom(id))?;
                    }
                    atoms.clone()
                }
                None => self.molecule.atoms().ids().collect(),
            };
            let determined = atoms.len() < 2
                || atoms
                    .iter()
                    .all(|id| self.components()[id.index()] == self.components()[atoms[0].index()]);
            return Ok(if determined {
                Solution::Determined(())
            } else {
                Solution::Contradictory(
                    MoleculeConstraintInvariantsContradiction {
                        constraint: constraint.clone(),
                    }
                    .into(),
                )
            });
        }
        MoleculeConstraintInvariantsValidator
            .validate(self.molecule, constraint, self.config)
            .map(|outcome| outcome.map_contradiction(ConstraintInvariantsContradiction::from))
    }

    fn rings(&mut self) -> &RingViews<'a> {
        if self.rings.is_none() {
            self.rings = Some(self.molecule.rings(
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
                self.molecule,
                self.config.connected_components_algorithm,
            ));
        }
        self.component_by_atom
            .as_deref()
            .expect("components initialized")
    }

    fn require(&self, entity: Entity) -> Result<(), ConstraintInvariantsError> {
        let present = match entity {
            Entity::Atom(id) => self.molecule.atoms().contains(id),
            Entity::Bond(id) => self.molecule.bonds().contains(id),
            Entity::DativeBond(id) => self.molecule.dative_bonds().contains(id),
            Entity::AromaticSystem(id) => self.molecule.aromatic_systems().contains(id),
            Entity::MulticenterBond(id) => self.molecule.multicenter_bonds().contains(id),
            Entity::NoncovalentBond(id) => self.molecule.noncovalent_bonds().contains(id),
            Entity::StereoAtom(id) => self.molecule.stereo_atoms().contains(id),
            Entity::StereoBond(id) => self.molecule.stereo_bonds().contains(id),
        };
        if present {
            Ok(())
        } else {
            Err(ConstraintInvariantsError::InvalidReference { entity })
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
    contradiction: impl FnOnce() -> ConstraintInvariantsContradiction,
) -> Solution<(), ConstraintInvariantsContradiction> {
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
    contradiction: impl FnOnce() -> ConstraintInvariantsContradiction,
) -> Solution<(), ConstraintInvariantsContradiction>
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
    outcome: Solution<(), ConstraintInvariantsContradiction>,
    any_underdetermined: &mut bool,
) -> Option<ConstraintInvariantsContradiction> {
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
    use proptest::prelude::*;
    use rstest::{fixture, rstest};
    use umol_chem::element::Element;
    use umol_graph_core::{ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm};
    use umol_graph_ir::ir::{
        AtomConstraintForm, AtomForm, AtomId, BondConstraintForm, BondId, Constraints,
        DativeBondId, MoleculeConstraint, MoleculeEntries, RelationalConstraint, RingScope,
    };
    use umol_graph_ir::mol_dsl;

    use super::*;

    const CONFIG: ConstraintValidateConfig = ConstraintValidateConfig {
        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
        connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
    };

    #[fixture]
    fn molecule() -> Molecule {
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#)
    }

    #[rstest]
    #[case::incidence(
        IncidenceConstraintInvariantsContradiction::Atom {
            atom: AtomId(2),
            constraint: AtomConstraintForm::valence(4),
        },
        ConstraintInvariantsContradiction::Incidence(IncidenceConstraintInvariantsContradiction::Atom {
            atom: AtomId(2),
            constraint: AtomConstraintForm::valence(4),
        })
    )]
    #[case::ring(
        RingConstraintInvariantsContradiction::Bond {
            bond: BondId(3),
            constraint: BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        },
        ConstraintInvariantsContradiction::Ring(RingConstraintInvariantsContradiction::Bond {
            bond: BondId(3),
            constraint: BondConstraintForm::ring_membership(RingScope::Size(6), 1),
        })
    )]
    #[case::relational(
        RelationalConstraintInvariantsContradiction {
            constraint: RelationalConstraint::DativeBondDonor {
                bond: DativeBondId(1),
                atom: AtomId(2),
            },
        },
        ConstraintInvariantsContradiction::Relational(RelationalConstraintInvariantsContradiction {
            constraint: RelationalConstraint::DativeBondDonor {
                bond: DativeBondId(1),
                atom: AtomId(2),
            },
        })
    )]
    #[case::molecule(
        MoleculeConstraintInvariantsContradiction {
            constraint: MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(1), AtomId(2)]),
            },
        },
        ConstraintInvariantsContradiction::Molecule(MoleculeConstraintInvariantsContradiction {
            constraint: MoleculeConstraint::Connected {
                atoms: Some(vec![AtomId(1), AtomId(2)]),
            },
        })
    )]
    fn test_constraint_contradiction_from<T>(
        #[case] input: T,
        #[case] expected: ConstraintInvariantsContradiction,
    ) where
        T: Into<ConstraintInvariantsContradiction>,
    {
        assert_eq!(input.into(), expected);
    }

    #[rstest]
    #[case::invalid_reference(
        ConstraintInvariantsError::InvalidReference {
            entity: Entity::Atom(AtomId(4)),
        },
        "constraint references unavailable atom 4"
    )]
    #[case::dative_ring(
        ConstraintInvariantsError::DativeBondRingMembershipUnsupported {
            bond: DativeBondId(5),
        },
        "ring membership is not defined for dative bond DativeBondId(5)"
    )]
    fn test_constraint_error_display(
        #[case] input: ConstraintInvariantsError,
        #[case] expected: &str,
    ) {
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
        mut molecule: Molecule,
        #[case] constraint: Constraint,
        #[case] expected: Solution<(), ConstraintInvariantsContradiction>,
    ) {
        molecule.constraints_mut().push(constraint);

        assert_eq!(
            ConstraintInvariantsValidator::new(CONFIG).validate(&molecule),
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
        mut molecule: Molecule,
        #[case] constraint: Constraint,
    ) {
        molecule.constraints_mut().push(constraint.clone());

        assert_eq!(
            ConstraintInvariantsValidator::new(CONFIG).validate(&molecule),
            Ok(Solution::Contradictory(
                ConstraintInvariantsContradiction::Logical { constraint }
            ))
        );
    }

    #[rstest]
    #[case::determined(AtomConstraintForm::valence(0))]
    #[case::underdetermined(AtomConstraintForm::total_hydrogens(1))]
    #[case::contradictory(AtomConstraintForm::valence(1))]
    fn test_constraint_validator_inline_top_level_agreement(
        molecule: Molecule,
        #[case] constraint: AtomConstraintForm,
    ) {
        let mut inline = molecule.clone().edit();
        inline
            .atom_mut(AtomId(0))
            .attributes
            .constraints
            .set(constraint.clone());
        let inline = inline.build();
        let mut top_level = molecule.edit();
        top_level
            .constraints_mut()
            .push(Constraint::Atom(AtomId(0), constraint));
        let top_level = top_level.build();
        let validator = ConstraintInvariantsValidator::new(CONFIG);

        assert_eq!(validator.validate(&inline), validator.validate(&top_level));
    }

    fn molecule_with(constraint: Constraint) -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            constraints: Constraints::from(constraint),
            ..MoleculeEntries::default()
        })
    }

    proptest! {
        #[test]
        fn test_constraint_and_permutation_invariant(
            values in prop::collection::vec(0_i64..=3, 0..=8),
        ) {
            let constraints: Vec<_> = values
                .iter()
                .map(|&value| Constraint::Atom(AtomId(0), AtomConstraintForm::valence(value)))
                .collect();
            let mut reversed = constraints.clone();
            reversed.reverse();
            let validator = ConstraintInvariantsValidator::new(CONFIG);
            let forward = validator.validate(&molecule_with(Constraint::And(constraints))).unwrap();
            let reverse = validator.validate(&molecule_with(Constraint::And(reversed))).unwrap();

            prop_assert_eq!(forward.is_determined(), reverse.is_determined());
            prop_assert_eq!(forward.is_underdetermined(), reverse.is_underdetermined());
            prop_assert_eq!(forward.is_contradictory(), reverse.is_contradictory());
        }

        #[test]
        fn test_constraint_or_permutation_invariant(
            values in prop::collection::vec(0_i64..=3, 0..=8),
        ) {
            let constraints: Vec<_> = values
                .iter()
                .map(|&value| Constraint::Atom(AtomId(0), AtomConstraintForm::valence(value)))
                .collect();
            let mut reversed = constraints.clone();
            reversed.reverse();
            let validator = ConstraintInvariantsValidator::new(CONFIG);
            let forward = validator.validate(&molecule_with(Constraint::Or(constraints))).unwrap();
            let reverse = validator.validate(&molecule_with(Constraint::Or(reversed))).unwrap();

            prop_assert_eq!(forward.is_determined(), reverse.is_determined());
            prop_assert_eq!(forward.is_underdetermined(), reverse.is_underdetermined());
            prop_assert_eq!(forward.is_contradictory(), reverse.is_contradictory());
        }

        #[test]
        fn test_constraint_double_negation(value in 0_i64..=3) {
            let constraint = Constraint::Atom(AtomId(0), AtomConstraintForm::valence(value));
            let double_negation = Constraint::Not(Box::new(Constraint::Not(Box::new(
                constraint.clone(),
            ))));
            let validator = ConstraintInvariantsValidator::new(CONFIG);
            let direct = validator.validate(&molecule_with(constraint)).unwrap();
            let negated = validator.validate(&molecule_with(double_negation)).unwrap();

            prop_assert_eq!(direct.is_determined(), negated.is_determined());
            prop_assert_eq!(direct.is_underdetermined(), negated.is_underdetermined());
            prop_assert_eq!(direct.is_contradictory(), negated.is_contradictory());
        }

        #[test]
        fn test_constraint_inline_top_level_leaf_agreement(value in 0_i64..=3) {
            let constraint = AtomConstraintForm::valence(value);
            let inline = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C)
                    .with_constraint(constraint.clone())],
                ..MoleculeEntries::default()
            });
            let top_level = molecule_with(Constraint::Atom(AtomId(0), constraint));
            let validator = ConstraintInvariantsValidator::new(CONFIG);

            prop_assert_eq!(validator.validate(&inline), validator.validate(&top_level));
        }
    }
}
