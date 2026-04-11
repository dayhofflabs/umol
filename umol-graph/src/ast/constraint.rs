//! Constraint AST: declarative facts over MoleculeAst consumed by the matcher and resolver.

use umol_shared::{SpinStateAst, ValueAst};

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
