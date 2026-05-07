//! Composite resolver: chains the per-entity resolvers (valence,
//! aromaticity, bonds, multicenter bonds) on a single `MoleculeAst`.
//! One-shot per pass — no fixpoint loop, no `ResolverCell`. Topology stays
//! invariant across narrowing, so intermediate views remain valid.
//!
//! `Determined` requires every entity (atoms, bonds, dative bonds, aromatic
//! systems, multicenter bonds, noncovalent bonds) to be ground. Sub-resolvers
//! report only completion or contradiction; the top-level `Resolver` decides
//! the global ground-status verdict.

pub mod aromaticity;
pub mod bonds;
pub mod multicenter;
pub mod valence;

pub use aromaticity::AromaticityResolver;
pub use bonds::{BondsContradiction, BondsError, BondsResolver};
pub use multicenter::{MulticenterBondsContradiction, MulticenterBondsError, MulticenterBondsResolver};
pub use valence::{ValenceContradiction, ValenceError, ValenceResolver};

use thiserror::Error;
use umol_ast::ast::MoleculeAst;

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
use crate::ops::config::ChemistryModel;
use crate::ops::solution::Solution;

#[derive(Clone, Debug)]
pub struct Resolver {
    pub valence: ValenceResolver,
    pub aromaticity: AromaticityResolver,
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
    Bonds(#[from] BondsError),
    #[error(transparent)]
    MulticenterBonds(#[from] MulticenterBondsError),
}

impl Resolver {
    pub fn new(model: &ChemistryModel) -> Self {
        Self {
            valence: ValenceResolver::new(&model.valence),
            aromaticity: AromaticityResolver::new(&model.aromaticity),
            bonds: BondsResolver::new(),
            multicenter_bonds: MulticenterBondsResolver::new(),
        }
    }

    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), ResolverContradiction>, ResolverError> {
        match self.valence.resolve(ast)? {
            Solution::Determined(()) | Solution::Underdetermined(()) => {}
            Solution::Contradictory(c) => {
                return Ok(Solution::Contradictory(ResolverContradiction::Valence(c)));
            }
        }
        match self.aromaticity.resolve(ast)? {
            Solution::Determined(()) | Solution::Underdetermined(()) => {}
            Solution::Contradictory(c) => {
                return Ok(Solution::Contradictory(
                    ResolverContradiction::Aromaticity(c),
                ));
            }
        }
        match self.bonds.resolve(ast)? {
            Solution::Determined(()) | Solution::Underdetermined(()) => {}
            Solution::Contradictory(c) => {
                return Ok(Solution::Contradictory(ResolverContradiction::Bonds(c)));
            }
        }
        match self.multicenter_bonds.resolve(ast)? {
            Solution::Determined(()) | Solution::Underdetermined(()) => {}
            Solution::Contradictory(c) => {
                return Ok(Solution::Contradictory(
                    ResolverContradiction::MulticenterBonds(c),
                ));
            }
        }

        Ok(if molecule_all_ground(ast) {
            Solution::Determined(())
        } else {
            Solution::Underdetermined(())
        })
    }
}

fn molecule_all_ground(ast: &MoleculeAst) -> bool {
    ast.atoms().iter().all(|v| v.data.is_ground())
        && ast.bonds().iter().all(|v| v.data.is_ground())
        && ast.dative_bonds().iter().all(|v| v.data.is_ground())
        && ast.aromatic_systems().iter().all(|v| v.data.is_ground())
        && ast.multicenter_bonds().iter().all(|v| v.data.is_ground())
        && ast.noncovalent_bonds().iter().all(|v| v.data.is_ground())
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::mol_zeroed;
    use umol_ast::ast::MoleculeAst;
    use umol_shared::element::Element;

    use super::*;
    use crate::ops::config::{AromaticityModel, ChemistryModel, ElementScope, RingLimits, ValenceModel};
    use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

    fn ground_methane() -> MoleculeAst {
        mol_zeroed!(r#"{:atoms ["C #h4"] :bonds []}"#)
    }

    fn benzene_pinned() -> MoleculeAst {
        mol_zeroed!(r#"{
            :atoms ["C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#)
    }

    fn counts_model() -> ChemistryModel {
        ChemistryModel {
            valence: ValenceModel::Counts {
                table: ValenceTable::default_table().clone(),
                allow_implicit_hydrogens: true,
            },
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
        }
    }

    fn typing_model() -> ChemistryModel {
        ChemistryModel {
            valence: ValenceModel::AtomTyping {
                registry: AtomTypeRegistry::default_registry().clone(),
            },
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
        }
    }

    #[rstest]
    fn test_resolver_resolve_ground_methane_determined() {
        let resolver = Resolver::new(&counts_model());
        let mut ast = ground_methane();
        let solution = resolver.resolve(&mut ast).unwrap();
        assert!(matches!(solution, Solution::Determined(())));
    }

    #[rstest]
    fn test_resolver_resolve_ground_methane_atom_typing_determined() {
        let resolver = Resolver::new(&typing_model());
        let mut ast = ground_methane();
        let solution = resolver.resolve(&mut ast).unwrap();
        assert!(matches!(solution, Solution::Determined(())));
    }

    #[rstest]
    fn test_resolver_resolve_benzene_writes_aromatic_system() {
        let resolver = Resolver::new(&counts_model());
        let mut ast = benzene_pinned();
        let solution = resolver.resolve(&mut ast).unwrap();
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(ast.aromatic_systems().count(), 1);
    }
}
