//! Structural stereo resolver: adds a `:stereo-atom` / `:stereo-bond` element
//! for each atom `#T` / bond `#C` that can be realized, using the canonical
//! ligand frame and copying the coset verbatim (raise already stored it in that
//! frame). Mirrors `AromaticityResolver`; computes no stereo symmetry; runs
//! after aromaticity (so aromatic-system membership is known). Skips sites that
//! already bear a stereo element, so re-runs are a no-op.

use thiserror::Error;
use umol_ast::ast::{
    AsLit, AtomId, BondId, CisTransStereoAst, MoleculeAst, StereoAtomAst, StereoBondAst,
    StereoKind, StereoLigand, StereoLigandKind, TetrahedralStereoAst,
};
use umol_utils::solution::Solution;

use crate::ops::model::StereoModel;

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
            .filter_map(|id| self.resolve_stereo_atom(ast, id))
            .collect();
        let bond_adds: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)> = ast
            .bonds()
            .ids()
            .filter_map(|id| self.resolve_stereo_bond(ast, id))
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

    fn resolve_stereo_atom(
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
        let TetrahedralStereoAst::Stereo(coset) = atom
            .ast
            .constraints
            .tetrahedral_stereo()
            .unwrap_or(&TetrahedralStereoAst::Undetermined)
        else {
            return None;
        };
        let coset = coset.clone();
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

        Some((id, ligands, StereoAtomAst::new(kind, coset)))
    }

    fn resolve_stereo_bond(
        &self,
        ast: &MoleculeAst,
        id: BondId,
    ) -> Option<(BondId, Vec<StereoLigand>, StereoBondAst)> {
        if ast.stereo_bonds().has_coincident(id) {
            return None;
        }
        let bond = ast.bond(id);

        let kind = StereoKind::CisTrans;
        let CisTransStereoAst::Stereo(coset) = bond
            .ast
            .constraints
            .cis_trans_stereo()
            .unwrap_or(&CisTransStereoAst::Undetermined)
        else {
            return None;
        };
        let coset = coset.clone();
        let model = self.model.kind_model(kind)?;
        // Endpoints are canonical (min, max) = raise's (start, end), so side_a/side_b
        // match the coset frame raise stored.
        let [a, b] = bond.atom_ids();
        if !model.scope.contains(ast.atom(a).element().as_lit()?)
            || !model.scope.contains(ast.atom(b).element().as_lit()?)
        {
            return None;
        }

        let side_a = self.bond_side_ligands(ast, a, b)?;
        let side_b = self.bond_side_ligands(ast, b, a)?;
        let ligands = vec![side_a[0], side_a[1], side_b[0], side_b[1]];

        Some((id, ligands, StereoBondAst::new(kind, coset)))
    }

    /// The two ligands of one double-bond end, in `cis_trans_side` order: the
    /// `atom`'s neighbors (ascending, excluding the `partner` across the bond),
    /// the first as `Atom`; the second as `Atom`, or a single virtual ligand
    /// (implicit hydrogen / lone pair) appended when the end has one substituent.
    fn bond_side_ligands(
        &self,
        ast: &MoleculeAst,
        atom: AtomId,
        partner: AtomId,
    ) -> Option<[StereoLigand; 2]> {
        let view = ast.atom(atom);
        let mut substituents = view
            .neighbors()
            .map(|n| n.atom_id())
            .filter(|&n| n != partner);
        let first = StereoLigand::new(substituents.next()?, StereoLigandKind::Atom);
        let second = match substituents.next() {
            Some(second) => StereoLigand::new(second, StereoLigandKind::Atom),
            None => {
                let virtual_kind = if view.implicit_hydrogens().as_lit()? >= 1 {
                    StereoLigandKind::ImplicitHydrogen
                } else if view.lone_pairs().as_lit()? >= 1 {
                    StereoLigandKind::LonePair
                } else {
                    return None;
                };
                StereoLigand::new(atom, virtual_kind)
            }
        };
        Some([first, second])
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AtomId, BondId, StereoAtomId, StereoCosetAst, StereoKind, StereoLigandKind,
    };
    use umol_ast::mol_dsl_ground;
    use umol_chem::element::Element;
    use umol_utils::solution::Solution;

    use super::StereoResolver;
    use crate::ops::model::{ElementScope, StereoKindModel, StereoModel};

    type StereoAtomData = (
        AtomId,
        StereoKind,
        StereoCosetAst,
        Vec<(AtomId, StereoLigandKind)>,
    );
    type StereoBondData = (
        BondId,
        StereoKind,
        StereoCosetAst,
        Vec<(AtomId, StereoLigandKind)>,
    );

    #[rstest]
    #[case::tetrahedral_atom(
        r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"] :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#,
        vec![(AtomId(1), StereoKind::Tetrahedral, StereoCosetAst::Lit(1),
        vec![(AtomId(0), StereoLigandKind::Atom), (AtomId(2), StereoLigandKind::Atom),
             (AtomId(3), StereoLigandKind::Atom), (AtomId(1), StereoLigandKind::ImplicitHydrogen)])], vec![])]
    #[case::cis_trans_bond(
        r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"] :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]}"#,
        vec![], vec![(BondId(1), StereoKind::CisTrans, StereoCosetAst::Lit(1),
        vec![(AtomId(0), StereoLigandKind::Atom), (AtomId(1), StereoLigandKind::ImplicitHydrogen),
             (AtomId(3), StereoLigandKind::Atom), (AtomId(2), StereoLigandKind::ImplicitHydrogen)])])]
    #[case::two_coordinate_skip(r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#, vec![], vec![])]
    #[case::aromatic_skip(
        r##"{:atoms ["C #h #T1" "C #h" "C #h" "C #h" "C #h" "C #h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "*#e6"}]}"##, vec![], vec![])]
    #[case::no_stereo( r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#, vec![], vec![])]
    fn test_stereo_resolver_resolve(
        #[case] input: &str,
        #[case] expected_atoms: Vec<StereoAtomData>,
        #[case] expected_bonds: Vec<StereoBondData>,
    ) {
        let mut ast = mol_dsl_ground!(input);
        let solution = StereoResolver::new(&StereoModel::default())
            .resolve(&mut ast)
            .unwrap();
        assert!(matches!(solution, Solution::Determined(())));

        let atoms: Vec<StereoAtomData> = ast
            .stereo_atoms()
            .iter()
            .map(|s| {
                (
                    s.site().id,
                    s.kind(),
                    s.coset().clone(),
                    s.ligands().map(|l| (l.atom_id(), l.kind())).collect(),
                )
            })
            .collect();
        let bonds: Vec<StereoBondData> = ast
            .stereo_bonds()
            .iter()
            .map(|s| {
                (
                    s.site().id,
                    s.kind(),
                    s.coset().clone(),
                    s.ligands().map(|l| (l.atom_id(), l.kind())).collect(),
                )
            })
            .collect();
        assert_eq!(atoms, expected_atoms);
        assert_eq!(bonds, expected_bonds);
    }

    #[rstest]
    fn test_stereo_resolver_resolve_out_of_scope() {
        let mut model = StereoModel::default();
        model.kind_models[StereoKind::Tetrahedral as usize] = Some(StereoKindModel {
            scope: ElementScope::AllowList(vec![Element::N]),
            fluxionality: false,
        });
        let mut ast = mol_dsl_ground!(
            r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"] :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#
        );
        StereoResolver::new(&model).resolve(&mut ast).unwrap();
        assert_eq!(ast.stereo_atoms().iter().count(), 0);
    }

    #[rstest]
    fn test_stereo_resolver_resolve_idempotent() {
        let resolver = StereoResolver::new(&StereoModel::default());
        let mut ast = mol_dsl_ground!(
            r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"] :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#
        );
        resolver.resolve(&mut ast).unwrap();
        resolver.resolve(&mut ast).unwrap();
        assert_eq!(ast.stereo_atoms().iter().count(), 1);
        let s = ast.stereo_atom(StereoAtomId(0));
        assert_eq!(s.kind(), StereoKind::Tetrahedral);
        assert_eq!(*s.coset(), StereoCosetAst::Lit(1));
    }
}
