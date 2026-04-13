//! Molecule structural AST.

use umol_edn::{FromEdn, ToEdn};

use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::config::MoleculeAstConfig;
use crate::ast::constraint::MoleculeConstraint;
use crate::ast::error::GroundError;
use crate::ast::Ast;

#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct LocalizedBond {
    pub a: usize,
    pub b: usize,
    pub bond: BondAst,
}

#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct DativeBond {
    pub donor: usize,
    pub acceptor: usize,
    pub bond: BondAst,
}

#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct AromaticSystem {
    pub atoms: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct MulticenterBond {
    pub atoms: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct NoncovalentBond {
    pub a: usize,
    pub b: usize,
    pub bond: BondAst,
}

/// Molecule AST: structural representation of a molecule (ground or pattern).
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct MoleculeAst {
    pub atoms: Vec<AtomAst>,
    pub bonds: Vec<LocalizedBond>,
    pub dative_bonds: Vec<DativeBond>,
    pub aromatic_systems: Vec<AromaticSystem>,
    pub multicenter_bonds: Vec<MulticenterBond>,
    pub noncovalent_bonds: Vec<NoncovalentBond>,
    pub constraints: Vec<MoleculeConstraint>,
}

impl MoleculeAst {
    pub fn is_ground(&self) -> bool {
        self.atoms.iter().all(|a| a.is_ground())
            && self.bonds.iter().all(|b| b.bond.is_ground())
            && self.dative_bonds.iter().all(|b| b.bond.is_ground())
            && self.noncovalent_bonds.iter().all(|b| b.bond.is_ground())
            && self.constraints.iter().all(|c| c.is_ground_assertion())
    }
}

impl Ast for MoleculeAst {
    type Config = MoleculeAstConfig;
}

/// A `MoleculeAst` whose fields are all concrete and whose constraints are
/// all ground assertions. The invariant is checked once at construction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GroundMolecule(MoleculeAst);

impl GroundMolecule {
    pub fn new(ast: MoleculeAst) -> Result<Self, GroundError> {
        if ast.is_ground() {
            Ok(Self(ast))
        } else {
            Err(GroundError)
        }
    }

    pub fn as_ast(&self) -> &MoleculeAst {
        &self.0
    }

    pub fn into_ast(self) -> MoleculeAst {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::atom_ast::ElementAst;
    use umol_shared::spin_ast::SpinStateAst;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::constraint::{DerivedPred, RelationRefs};

    fn ground_ast() -> MoleculeAst {
        MoleculeAst {
            atoms: vec![AtomAst::from_element(umol_shared::element::Element::C)],
            bonds: vec![],
            dative_bonds: vec![],
            aromatic_systems: vec![],
            multicenter_bonds: vec![],
            noncovalent_bonds: vec![],
            constraints: vec![],
        }
    }

    #[rstest]
    #[case::empty(MoleculeAst::default(), true)]
    #[case::ground_atom(ground_ast(), true)]
    #[case::ground_with_lit_constraint(
        MoleculeAst {
            constraints: vec![MoleculeConstraint::Derived {
                predicate: DerivedPred::TotalCharge(ValueAst::Lit(-1)),
                refs: RelationRefs::default(),
            }],
            ..ground_ast()
        },
        true,
    )]
    #[case::wildcard_element(
        MoleculeAst {
            atoms: vec![AtomAst::new(ElementAst::Wildcard)],
            ..MoleculeAst::default()
        },
        false,
    )]
    #[case::wildcard_bond_order(
        MoleculeAst {
            atoms: vec![
                AtomAst::from_element(umol_shared::element::Element::C),
                AtomAst::from_element(umol_shared::element::Element::O),
            ],
            bonds: vec![LocalizedBond { a: 0, b: 1, bond: BondAst::new(ValueAst::Wildcard) }],
            ..MoleculeAst::default()
        },
        false,
    )]
    #[case::non_ground_constraint(
        MoleculeAst {
            constraints: vec![MoleculeConstraint::Derived {
                predicate: DerivedPred::TotalSpin(SpinStateAst::Wildcard),
                refs: RelationRefs::default(),
            }],
            ..ground_ast()
        },
        false,
    )]
    #[case::sub_pattern_constraint(
        MoleculeAst {
            constraints: vec![MoleculeConstraint::SubPattern {
                anchor: 0,
                pattern: Box::new(MoleculeAst::default()),
            }],
            ..ground_ast()
        },
        false,
    )]
    fn test_molecule_ast_is_ground(#[case] ast: MoleculeAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[test]
    fn test_ground_molecule_new() {
        let ast = ground_ast();
        let ground = GroundMolecule::new(ast.clone());
        assert!(ground.is_ok());
        assert_eq!(ground.unwrap().as_ast(), &ast);
    }

    #[test]
    fn test_ground_molecule_new_error() {
        let ast = MoleculeAst {
            atoms: vec![AtomAst::new(ElementAst::Wildcard)],
            ..MoleculeAst::default()
        };
        assert_eq!(GroundMolecule::new(ast), Err(GroundError));
    }

    #[test]
    fn test_ground_molecule_into_ast() {
        let ast = ground_ast();
        let ground = GroundMolecule::new(ast.clone()).unwrap();
        assert_eq!(ground.into_ast(), ast);
    }
}
