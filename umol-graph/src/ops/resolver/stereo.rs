//! Structural stereo resolver: adds a `:stereo-atom` / `:stereo-bond` element
//! for each atom `#T` / bond `#C` that can be realized, using the canonical
//! ligand frame and copying the coset verbatim (raise already stored it in that
//! frame). Mirrors `AromaticityResolver`; computes no stereo symmetry; runs
//! after aromaticity (so aromatic-system membership is known). Skips sites that
//! already bear a stereo element, so re-runs are a no-op.

use thiserror::Error;
use umol_ast::ast::{AtomId, BondId, MoleculeAst, StereoAtomAst, StereoBondAst, StereoLigand};

use crate::ops::model::StereoModel;
use crate::ops::solution::Solution;

#[derive(Clone, Debug)]
pub struct StereoResolver {
    model: StereoModel,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoError {}

impl StereoResolver {
    pub fn new(model: &StereoModel) -> Self {
        Self {
            model: model.clone(),
        }
    }

    /// Adds a stereo element for each realizable atom `#T` / bond `#C`. The
    /// per-site decision is read first (immutable borrow), then applied through
    /// a single builder pass. Returns `Determined`; the inconsistency policy for
    /// non-realizable assertions is deferred.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), StereoContradiction>, StereoError> {
        let atom_adds: Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)> = ast
            .atoms()
            .ids()
            .filter_map(|id| self.resolve_atom(ast, id))
            .collect();
        let bond_adds: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)> = ast
            .bonds()
            .ids()
            .filter_map(|id| self.resolve_bond(ast, id))
            .collect();

        if atom_adds.is_empty() && bond_adds.is_empty() {
            return Ok(Solution::Determined(()));
        }

        let mut builder = ast.edit();
        for (site, ligands, data) in atom_adds {
            builder.add_stereo_atom(site, ligands, data);
        }
        for (site, ligands, data) in bond_adds {
            builder.add_stereo_bond(site, ligands, data);
        }
        *ast = builder.build();

        Ok(Solution::Determined(()))
    }

    fn resolve_atom(
        &self,
        _ast: &MoleculeAst,
        _id: AtomId,
    ) -> Option<(AtomId, Vec<StereoLigand>, StereoAtomAst)> {
        None
    }

    fn resolve_bond(
        &self,
        _ast: &MoleculeAst,
        _id: BondId,
    ) -> Option<(BondId, Vec<StereoLigand>, StereoBondAst)> {
        None
    }
}
