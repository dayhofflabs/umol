//! Structural stereo resolver: adds a `:stereo-atom` / `:stereo-bond` element
//! for each atom `#T` / bond `#C` that can be realized, using the canonical
//! ligand frame and copying the coset verbatim (raise already stored it in that
//! frame). Mirrors `AromaticityResolver`; computes no stereo symmetry; runs
//! after aromaticity (so aromatic-system membership is known). Skips sites that
//! already bear a stereo element, so re-runs are a no-op.

use thiserror::Error;
use umol_ast::ast::{
    AsLit, AtomId, BondId, MoleculeAst, StereoAtomAst, StereoBondAst, StereoConfigurationAst,
    StereoKind, StereoLigand, StereoLigandKind,
};

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
        ast: &MoleculeAst,
        id: AtomId,
    ) -> Option<(AtomId, Vec<StereoLigand>, StereoAtomAst)> {
        if ast.stereo_atoms().has_coincident(id) {
            return None;
        }
        let atom = ast.atom(id);
        if atom.is_in_aromatic_system() {
            return None;
        }

        let kind = StereoKind::Tetrahedral;
        let StereoConfigurationAst::Stereo(coset) = atom.ast.constraints.tetrahedral_stereo() else {
            return None;
        };
        let model = self.model.kind_model(kind)?;
        if !model.scope.contains(atom.element().as_lit()?) {
            return None;
        }

        let mut ligands: Vec<StereoLigand> = atom
            .neighbors()
            .map(|n| StereoLigand::new(n.atom_id(), StereoLigandKind::Atom))
            .collect();
        if ligands.len() + 1 == kind.degree() {
            let virtual_kind = if atom.implicit_hydrogens().as_lit()? >= 1 {
                StereoLigandKind::ImplicitHydrogen
            } else if atom.lone_pairs().as_lit()? >= 1 {
                StereoLigandKind::LonePair
            } else {
                return None;
            };
            ligands.push(StereoLigand::new(id, virtual_kind));
        }
        if ligands.len() != kind.degree() {
            return None;
        }

        Some((id, ligands, StereoAtomAst::new(kind, coset.simplify(kind))))
    }

    fn resolve_bond(
        &self,
        _ast: &MoleculeAst,
        _id: BondId,
    ) -> Option<(BondId, Vec<StereoLigand>, StereoBondAst)> {
        None
    }
}
