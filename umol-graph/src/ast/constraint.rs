//! Constraint AST: declarative facts over MoleculeAst consumed by the matcher and resolver.

use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::molecule::MoleculeAst;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraint {
    SubPattern {
        anchor: usize,
        pattern: Box<MoleculeAst>,
    },
    Derived {
        predicate: DerivedPred,
        refs: RelationRefs,
    },
    Matcher(MatcherFlag),
    And(Vec<MoleculeConstraint>),
    Or(Vec<MoleculeConstraint>),
    Not(Box<MoleculeConstraint>),
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RelationRefs {
    pub atoms: Vec<usize>,
    pub bonds: Vec<usize>,
    pub dative_bonds: Vec<usize>,
    pub aromatic_systems: Vec<usize>,
    pub multicenter_bonds: Vec<usize>,
    pub noncovalent_bonds: Vec<usize>,
}

impl RelationRefs {
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
            && self.bonds.is_empty()
            && self.dative_bonds.is_empty()
            && self.aromatic_systems.is_empty()
            && self.multicenter_bonds.is_empty()
            && self.noncovalent_bonds.is_empty()
    }

    pub fn atoms(atoms: Vec<usize>) -> Self {
        Self { atoms, ..Self::default() }
    }

    pub fn bonds(bonds: Vec<usize>) -> Self {
        Self { bonds, ..Self::default() }
    }

    pub fn dative_bonds(dative_bonds: Vec<usize>) -> Self {
        Self { dative_bonds, ..Self::default() }
    }

    pub fn aromatic_systems(aromatic_systems: Vec<usize>) -> Self {
        Self { aromatic_systems, ..Self::default() }
    }

    pub fn multicenter_bonds(multicenter_bonds: Vec<usize>) -> Self {
        Self { multicenter_bonds, ..Self::default() }
    }

    pub fn noncovalent_bonds(noncovalent_bonds: Vec<usize>) -> Self {
        Self { noncovalent_bonds, ..Self::default() }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DerivedPred {
    TotalCharge(ValueAst),
    TotalSpin(SpinStateAst),
    ValenceSum(ValueAst),
    AromaticElectronCount(ValueAst),
    RingSize(ValueAst),
    InRing,
    NotInRing,
    InRelation(RelationSym),
    NotInRelation(RelationSym),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MatcherFlag {
    Injective,
    NonInjective,
    Induced,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RelationSym {
    Atoms,
    Bonds,
    DativeBonds,
    AromaticSystems,
    MulticenterBonds,
    NoncovalentBonds,
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::value_ast::ValueAst;

    use super::*;

    fn charge(n: i64) -> MoleculeConstraint {
        MoleculeConstraint::Derived {
            predicate: DerivedPred::TotalCharge(ValueAst::Lit(n)),
            refs: RelationRefs::default(),
        }
    }

    fn in_ring_atom(atom: usize) -> MoleculeConstraint {
        MoleculeConstraint::Derived {
            predicate: DerivedPred::InRing,
            refs: RelationRefs::atoms(vec![atom]),
        }
    }

    #[rstest]
    #[case::and_pair(
        MoleculeConstraint::And(vec![charge(0), in_ring_atom(1)]),
        MoleculeConstraint::And(vec![charge(0), in_ring_atom(1)]),
        true,
    )]
    #[case::and_order_matters(
        MoleculeConstraint::And(vec![charge(0), in_ring_atom(1)]),
        MoleculeConstraint::And(vec![in_ring_atom(1), charge(0)]),
        false,
    )]
    #[case::or_distinct_payload(
        MoleculeConstraint::Or(vec![charge(0), charge(1)]),
        MoleculeConstraint::Or(vec![charge(0), charge(2)]),
        false,
    )]
    #[case::not_idempotent_eq(
        MoleculeConstraint::Not(Box::new(charge(-1))),
        MoleculeConstraint::Not(Box::new(charge(-1))),
        true,
    )]
    #[case::and_or_distinct(
        MoleculeConstraint::And(vec![charge(0)]),
        MoleculeConstraint::Or(vec![charge(0)]),
        false,
    )]
    fn test_molecule_constraint_combinators_eq(
        #[case] left: MoleculeConstraint,
        #[case] right: MoleculeConstraint,
        #[case] equal: bool,
    ) {
        assert_eq!(left == right, equal);
    }

    #[test]
    fn test_molecule_constraint_combinators_nested() {
        let inner = MoleculeConstraint::Or(vec![
            charge(-1),
            MoleculeConstraint::Not(Box::new(in_ring_atom(0))),
        ]);
        let outer = MoleculeConstraint::And(vec![charge(0), inner.clone()]);

        let MoleculeConstraint::And(children) = &outer else {
            panic!("expected And");
        };
        assert_eq!(children.len(), 2);
        assert_eq!(children[1], inner);
    }
}
