//! Molecule structural AST.

use index_vec::IndexVec;
use umol_edn::{FromEdn, ToEdn};
use umol_shared::value_ast::ValueAst;

use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::config::MoleculeAstConfig;
use crate::ast::constraint::MoleculeConstraint;
use crate::ast::error::GroundError;
use crate::ast::{Ast, AtomIdx};

/// Binary relation over two atoms with bond attributes.
///
/// For directed relations (dative, noncovalent), `source` is the donor/origin
/// and `target` is the acceptor/destination. For undirected relations
/// (localized bonds), the ordering is canonical (`source <= target`).
#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct BondTuple {
    pub source: AtomIdx,
    pub target: AtomIdx,
    pub bond: BondAst,
}

#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct AromaticSystem {
    pub atoms: Vec<AtomIdx>,
}

#[derive(Clone, Debug, PartialEq, Eq, FromEdn, ToEdn)]
pub struct MulticenterBond {
    pub atoms: Vec<AtomIdx>,
}

/// Molecule AST: structural representation of a molecule (ground or pattern).
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct MoleculeAst {
    pub atoms: IndexVec<AtomIdx, AtomAst>,
    pub bonds: Vec<BondTuple>,
    pub dative_bonds: Vec<BondTuple>,
    pub aromatic_systems: Vec<AromaticSystem>,
    pub multicenter_bonds: Vec<MulticenterBond>,
    pub noncovalent_bonds: Vec<BondTuple>,
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

    pub fn bond_order_sum(&self, atom: AtomIdx) -> Option<u8> {
        let mut sum: u8 = 0;
        for bond in &self.bonds {
            if bond.source == atom || bond.target == atom {
                match bond.bond.order {
                    ValueAst::Lit(n) => sum += n as u8,
                    _ => return None,
                }
            }
        }
        Some(sum)
    }

    pub fn dative_bond_order_sums(&self, atom: AtomIdx) -> (u8, u8) {
        let mut donated: u8 = 0;
        let mut accepted: u8 = 0;
        for bond in &self.dative_bonds {
            let order = match bond.bond.order {
                ValueAst::Lit(n) => n as u8,
                _ => continue,
            };
            if bond.source == atom {
                donated += order;
            } else if bond.target == atom {
                accepted += order;
            }
        }
        (donated, accepted)
    }

    pub fn is_in_aromatic_system(&self, atom: AtomIdx) -> bool {
        self.aromatic_systems
            .iter()
            .any(|sys| sys.atoms.contains(&atom))
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
    use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst, IsotopeAst};
    use umol_shared::spin::SpinState;
    use umol_shared::spin_ast::SpinStateAst;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::constraint::{DerivedPred, RelationRefs};

    fn ground_atom() -> AtomAst {
        AtomAst {
            element: ElementAst::Lit(umol_shared::element::Element::C),
            isotope_mass: IsotopeAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: HydrogenAst::Value(ValueAst::Lit(4)),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::Lit(SpinState::closed_shell()),
            valence: ValueAst::Lit(4),
            donated_pairs: ValueAst::Lit(0),
            accepted_pairs: ValueAst::Lit(0),
            aromatic_valence: AromaticValenceAst::NotAromatic,
            multicenter_valence: ValueAst::Lit(0),
        }
    }

    fn ground_ast() -> MoleculeAst {
        MoleculeAst {
            atoms: vec![ground_atom()].into(),
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
            atoms: vec![AtomAst::new(ElementAst::Undetermined)].into(),
            ..MoleculeAst::default()
        },
        false,
    )]
    #[case::wildcard_bond_order(
        MoleculeAst {
            atoms: vec![
                AtomAst::from_element(umol_shared::element::Element::C),
                AtomAst::from_element(umol_shared::element::Element::O),
            ].into(),
            bonds: vec![BondTuple { source: AtomIdx(0), target: AtomIdx(1), bond: BondAst::new(ValueAst::Undetermined) }],
            ..MoleculeAst::default()
        },
        false,
    )]
    #[case::non_ground_constraint(
        MoleculeAst {
            constraints: vec![MoleculeConstraint::Derived {
                predicate: DerivedPred::TotalSpin(SpinStateAst::default()),
                refs: RelationRefs::default(),
            }],
            ..ground_ast()
        },
        false,
    )]
    #[case::sub_pattern_constraint(
        MoleculeAst {
            constraints: vec![MoleculeConstraint::SubPattern {
                anchor: AtomIdx(0),
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
            atoms: vec![AtomAst::new(ElementAst::Undetermined)].into(),
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
