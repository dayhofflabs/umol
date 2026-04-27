//! Aromaticity perception engine.
//!
//! `AromaticityResolver` dispatches to one of three perception algorithms
//! (Hückel rule, HMO, Clar) selected by the `AromaticityModel`. Each
//! algorithm enumerates aromatic systems on the AST's ring set and writes
//! results back: per-system atoms list + `AromaticSystemAst` (with per-atom
//! `electrons` filled in), plus a `BondConstraint::Aromatic` flag on each
//! induced bond.
//!
//! Each sub-algorithm declares one error enum covering its failure modes;
//! the dispatcher classifies variants into `Solution::Contradictory`,
//! `Solution::Underdetermined`, or `Err(AromaticityError)`.

pub mod clar;
pub mod hmo;
pub mod hueckel_rule;

use thiserror::Error;
use umol_ast::ast::{BondConstraint, BondIdx, MoleculeAst, RingFamily};

use crate::ops::config::AromaticityModel;
use crate::ops::solution::Solution;

pub use clar::{ClarAromaticity, ClarError};
pub use hmo::{HmoAromaticity, HmoError, HmoOutput};
pub use hueckel_rule::HueckelRuleAromaticity;

/// Chemistry-level rejection: the algorithm decided the input doesn't
/// satisfy the model. Carried inside `Solution::Contradictory`.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityContradiction {
    #[error("hmo: invalid input: {0}")]
    HmoInvalidInput(String),
    #[error("clar: non-benzenoid input: {0}")]
    ClarNonBenzenoid(String),
}

/// Setup-level failure: parameter table or configuration gap. Returned in
/// `Err`, never inside `Solution`.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    #[error("hmo: missing parameters: {0}")]
    HmoMissingParameters(String),
}

#[derive(Clone, Debug)]
pub enum AromaticityResolver {
    HueckelRule(HueckelRuleAromaticity),
    Hmo(HmoAromaticity),
    Clar(ClarAromaticity),
}

impl AromaticityResolver {
    pub fn new(model: &AromaticityModel) -> Self {
        match model {
            AromaticityModel::HueckelRule { scope, ring_limits } => Self::HueckelRule(
                HueckelRuleAromaticity::new(scope.clone(), ring_limits.clone()),
            ),
            AromaticityModel::Hmo {
                scope,
                stabilization_threshold,
            } => Self::Hmo(HmoAromaticity::new(scope.clone(), *stabilization_threshold)),
            AromaticityModel::Clar { .. } => Self::Clar(ClarAromaticity),
        }
    }

    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), AromaticityContradiction>, AromaticityError> {
        let (family, max_ring_size) = self.ring_request();
        let rings = ast.rings(family, max_ring_size).clone();

        let systems = match self {
            Self::HueckelRule(m) => m.find_from_rings(ast, &rings),
            Self::Hmo(m) => match m.find_from_rings(ast, &rings) {
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
                    return Ok(Solution::Underdetermined(()));
                }
            },
            Self::Clar(m) => match m.find_from_rings(ast, &rings) {
                Ok(systems) => systems,
                Err(ClarError::NonBenzenoid(s)) => {
                    return Ok(Solution::Contradictory(
                        AromaticityContradiction::ClarNonBenzenoid(s),
                    ));
                }
            },
        };

        if systems.is_empty() {
            return Ok(Solution::Determined(()));
        }

        let mut sorted = systems;
        sorted.sort_by(|a, b| a.0.first().cmp(&b.0.first()));

        let mut builder = ast.edit();
        for (atoms, system_ast) in sorted {
            let _ = builder.add_aromatic_system(atoms, system_ast);
        }
        *ast = builder.build();

        let bond_ids: Vec<BondIdx> = ast
            .aromatic_systems()
            .iter()
            .flat_map(|view| view.bonds().collect::<Vec<_>>())
            .collect();
        for bond_id in bond_ids {
            let bond = ast.bond_mut(bond_id);
            bond.data.constraints.add(BondConstraint::Aromatic);
        }

        Ok(Solution::Determined(()))
    }

    fn ring_request(&self) -> (RingFamily, usize) {
        match self {
            Self::HueckelRule(m) => (RingFamily::Simple, m.ring_limits.max_ring_size),
            Self::Hmo(_) => (RingFamily::Simple, 22),
            Self::Clar(_) => (RingFamily::Simple, 6),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AromaticValenceAst, AtomAst, AtomConstraint, AtomIdx, BondAst, Constraints, MoleculeAst,
        ValueAst,
    };
    use umol_shared::element::Element;

    use super::*;
    use crate::ops::config::{AromaticityModel, ElementScope, RingLimits};

    fn aromatic(element: Element, pi: i64) -> AtomAst {
        let mut atom = AtomAst::from_element(element);
        atom.constraints
            .add(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(
                ValueAst::Lit(pi),
            )));
        atom
    }

    fn benzene() -> MoleculeAst {
        let atoms: Vec<AtomAst> = (0..6).map(|_| aromatic(Element::C, 1)).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomIdx(i), AtomIdx((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    fn pyrrole() -> MoleculeAst {
        let atoms = vec![
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ];
        let bonds: Vec<_> = (0..5)
            .map(|i| (AtomIdx(i), AtomIdx((i + 1) % 5), BondAst::from_order(1)))
            .collect();
        MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_aromaticity_resolver_hueckel_rule_benzene_writes_system() {
        let resolver = AromaticityResolver::new(&AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        });
        let mut ast = benzene();
        let solution = resolver.resolve(&mut ast).unwrap();
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(ast.aromatic_systems().count(), 1);
        let system = ast.aromatic_system(umol_ast::ast::AromaticSystemIdx(0));
        let atoms: Vec<AtomIdx> = system.atoms().collect();
        assert_eq!(atoms.len(), 6);
        let aromatic_bond_count = ast
            .bonds()
            .iter()
            .filter(|view| {
                view.data
                    .constraints
                    .iter()
                    .any(|c| c.kind() == umol_ast::ast::BondConstraintKind::Aromatic)
            })
            .count();
        assert_eq!(aromatic_bond_count, 6);
    }

    #[rstest]
    fn test_aromaticity_resolver_clar_rejects_heterocycle() {
        let resolver = AromaticityResolver::new(&AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        });
        let mut ast = pyrrole();
        let solution = resolver.resolve(&mut ast).unwrap();
        assert!(matches!(
            solution,
            Solution::Contradictory(AromaticityContradiction::ClarNonBenzenoid(_))
        ));
    }

    #[rstest]
    fn test_aromaticity_resolver_hueckel_rule_no_aromatic_atom_returns_determined() {
        let resolver = AromaticityResolver::new(&AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        });
        let atoms: Vec<AtomAst> = (0..6).map(|_| AtomAst::from_element(Element::C)).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomIdx(i), AtomIdx((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        let mut ast = MoleculeAst::new(
            atoms,
            bonds,
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        let solution = resolver.resolve(&mut ast).unwrap();
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(ast.aromatic_systems().count(), 0);
        let any_aromatic = ast.bonds().iter().any(|view| {
            view.data
                .constraints
                .iter()
                .any(|c| c.kind() == umol_ast::ast::BondConstraintKind::Aromatic)
        });
        assert!(!any_aromatic);
    }
}
