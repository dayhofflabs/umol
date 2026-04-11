//! Constraint AST: declarative facts over MoleculeAst consumed by the matcher and resolver.

use crate::ast::molecule::MoleculeAst;
use crate::ast::value::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MoleculeConstraint {
    SubPattern {
        anchor: usize,
        pattern: Box<MoleculeAst>,
    },
    Derived {
        predicate: DerivedPred,
        atoms: Vec<usize>,
    },
    Matcher(MatcherFlag),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DerivedPred {
    TotalCharge(ValueAst),
    TotalMultiplicity(ValueAst),
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
