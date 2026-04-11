//! Molecule structural AST.

use umol_data::SpinState;
use umol_edn::{FromEdn, ToEdn};

use crate::ast::atom::AtomAst;
use crate::ast::bond::BondAst;
use crate::ast::config::MoleculeAstConfig;
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
    pub charge: Option<i64>,
    pub spin: Option<SpinState>,
}

impl Ast for MoleculeAst {
    type Config = MoleculeAstConfig;
}
