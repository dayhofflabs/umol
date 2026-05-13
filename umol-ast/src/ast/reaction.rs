//! Reaction rule AST: L ← K → R as two molecule ASTs plus an atom map.

use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::molecule::MoleculeAst;

/// Double-pushout reaction rule.
///
/// `lhs` and `rhs` are the left- and right-hand-side molecule patterns.
/// `atom_map` pairs define the interface K: atoms shared between L and R.
/// Everything in L not mapped is deleted; everything in R not mapped is
/// created. Bonds and overlay relations in K are inferred from the
/// topology of L and R over the mapped atoms.
#[derive(Clone, Debug)]
pub struct ReactionRuleAst {
    pub lhs: MoleculeAst,
    pub rhs: MoleculeAst,
    pub atom_map: Vec<(AtomId, AtomId)>,
}

/// Assignment mapping L-entity indices to target (G) entity indices.
/// Produced by substructure matching; consumed by `apply_rule`.
#[derive(Clone, Debug, Default)]
pub struct Assignment {
    pub atoms: Vec<(AtomId, AtomId)>,
    pub bonds: Vec<(BondId, BondId)>,
    pub dative_bonds: Vec<(DativeBondId, DativeBondId)>,
    pub aromatic_systems: Vec<(AromaticSystemId, AromaticSystemId)>,
    pub multicenter_bonds: Vec<(MulticenterBondId, MulticenterBondId)>,
    pub noncovalent_bonds: Vec<(NoncovalentBondId, NoncovalentBondId)>,
}
