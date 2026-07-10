//! Composite resolver: chains the per-entity resolvers (valence,
//! aromaticity, bonds, multicenter bonds) on a single `MoleculeAst`.
//!
//! `Determined` requires every entity (atoms, bonds, dative bonds, aromatic
//! systems, multicenter bonds, noncovalent bonds) to be ground.

pub mod aromaticity;
pub mod bonds;
pub mod multicenter;
pub mod stereo;
pub mod valence;

use std::any::Any;

pub use aromaticity::AromaticityResolver;
pub use bonds::{BondsContradiction, BondsError, BondsResolver};
pub use multicenter::{
    MulticenterBondsContradiction, MulticenterBondsError, MulticenterBondsResolver,
};
pub use stereo::{StereoContradiction, StereoError, StereoResolver};
use thiserror::Error;
use umol_ast::ast::MoleculeAst;
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;
pub use valence::{ValenceContradiction, ValenceError, ValenceResolver};

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
use crate::ops::model::ChemistryModel;

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
        Self {
            valence: ValenceResolver::new(&model.valence),
            aromaticity: AromaticityResolver::new(&model.aromaticity),
            stereo: StereoResolver::new(&model.stereo),
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
                return Ok(Solution::Contradictory(ResolverContradiction::Aromaticity(
                    c,
                )));
            }
        }
        match self.stereo.resolve(ast)? {
            Solution::Determined(()) | Solution::Underdetermined(()) => {}
            Solution::Contradictory(c) => {
                return Ok(Solution::Contradictory(ResolverContradiction::Stereo(c)));
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

        Ok(if ast.is_ground() {
            Solution::Determined(())
        } else {
            Solution::Underdetermined(())
        })
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::*;
    use umol_ast::ast::MoleculeAst;
    use umol_ast::mol_dsl_ground;
    use umol_chem::element::Element;

    use super::*;
    use crate::ops::model::{
        AromaticityModel, AtomTypingModel, ChemistryModel, CountsModel, ElementScope, RingLimits,
        StereoModel, ValenceModel,
    };
    use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

    fn methane() -> MoleculeAst {
        mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#)
    }

    fn benzene() -> MoleculeAst {
        mol_dsl_ground!(
            r#"{
            :atoms ["C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#
        )
    }

    fn counts_model() -> ChemistryModel {
        ChemistryModel {
            valence: ValenceModel::Counts(CountsModel {
                table: Cow::Borrowed(ValenceTable::default_table()),
            }),
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
            stereo: StereoModel::default(),
        }
    }

    fn atom_typing_model() -> ChemistryModel {
        ChemistryModel {
            valence: ValenceModel::AtomTyping(AtomTypingModel {
                registry: Cow::Borrowed(AtomTypeRegistry::default_registry()),
            }),
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
            stereo: StereoModel::default(),
        }
    }

    #[rstest]
    fn test_resolver_resolve_ground_methane_determined() {
        let model = counts_model();
        let resolver = Resolver::new(&model);
        let mut ast = methane();
        let solution = resolver.resolve(&mut ast).unwrap();
        assert!(matches!(solution, Solution::Determined(())));
    }

    #[rstest]
    fn test_resolver_resolve_ground_methane_atom_typing_determined() {
        let model = atom_typing_model();
        let resolver = Resolver::new(&model);
        let mut ast = methane();
        let solution = resolver.resolve(&mut ast).unwrap();
        assert!(matches!(solution, Solution::Determined(())));
    }

    #[rstest]
    fn test_resolver_resolve_benzene_writes_aromatic_system() {
        let model = counts_model();
        let resolver = Resolver::new(&model);
        let mut ast = benzene();
        let solution = resolver.resolve(&mut ast).unwrap();
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(ast.aromatic_systems().count(), 1);
    }

    #[rstest]
    fn test_resolver_resolve_benzene_idempotent() {
        let model = counts_model();
        let resolver = Resolver::new(&model);
        let mut ast = benzene();
        resolver.resolve(&mut ast).unwrap();
        resolver.resolve(&mut ast).unwrap();
        assert_eq!(ast.aromatic_systems().count(), 1);
    }
}
