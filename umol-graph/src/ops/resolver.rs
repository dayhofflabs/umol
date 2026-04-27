//! Composite resolver: runs the valence resolver, then aromaticity
//! perception, on a single `MoleculeAst`. One-shot per pass — no fixpoint
//! loop, no `ResolverCell`. Topology stays invariant across narrowing, so
//! intermediate views remain valid.

use thiserror::Error;
use umol_ast::ast::MoleculeAst;

use crate::ops::aromaticity::{
    AromaticityContradiction, AromaticityError, AromaticityResolver,
};
use crate::ops::config::ChemistryModel;
use crate::ops::solution::Solution;
use crate::ops::valence::{ValenceContradiction, ValenceError, ValenceResolver};

#[derive(Clone, Debug)]
pub struct Resolver {
    pub valence: ValenceResolver,
    pub aromaticity: AromaticityResolver,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ResolverContradiction {
    #[error(transparent)]
    Valence(#[from] ValenceContradiction),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityContradiction),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ResolverError {
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
}

impl Resolver {
    pub fn new(model: &ChemistryModel) -> Self {
        Self {
            valence: ValenceResolver::new(&model.valence),
            aromaticity: AromaticityResolver::new(&model.aromaticity),
        }
    }

    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), ResolverContradiction>, ResolverError> {
        let mut underdetermined = false;

        match self.valence.resolve(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => underdetermined = true,
            Solution::Contradictory(c) => {
                return Ok(Solution::Contradictory(ResolverContradiction::Valence(c)));
            }
        }

        match self.aromaticity.resolve(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => underdetermined = true,
            Solution::Contradictory(c) => {
                return Ok(Solution::Contradictory(
                    ResolverContradiction::Aromaticity(c),
                ));
            }
        }

        Ok(if underdetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AromaticValenceAst, AtomAst, AtomConstraint, AtomIdx, BondAst, Constraints,
        ImplicitHydrogensAst, IsotopeAst, MoleculeAst, SpinStateAst, ValueAst,
    };
    use umol_shared::element::Element;

    use super::*;
    use crate::ops::config::{AromaticityModel, ChemistryModel, ElementScope, RingLimits, ValenceModel};
    use crate::ops::valence::{AtomTypeRegistry, ValenceTable};

    fn ground_methane() -> MoleculeAst {
        let mut c = AtomAst::from_element(Element::C);
        c.isotope_mass = IsotopeAst::Natural;
        c.charge = ValueAst::Lit(0);
        c.implicit_hydrogens = ImplicitHydrogensAst::Lit(4);
        c.lone_pairs = ValueAst::Lit(0);
        c.spin = SpinStateAst::new(0, 1);
        MoleculeAst::new(
            vec![c],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    fn benzene_pinned() -> MoleculeAst {
        let atoms: Vec<AtomAst> = (0..6)
            .map(|_| {
                let mut a = AtomAst::from_element(Element::C);
                a.isotope_mass = IsotopeAst::Natural;
                a.charge = ValueAst::Lit(0);
                a.implicit_hydrogens = ImplicitHydrogensAst::Lit(1);
                a.lone_pairs = ValueAst::Lit(0);
                a.spin = SpinStateAst::new(0, 1);
                a.constraints
                    .add(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(
                        ValueAst::Lit(1),
                    )));
                a
            })
            .collect();
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
