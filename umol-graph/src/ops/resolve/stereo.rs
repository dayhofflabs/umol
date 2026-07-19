//! Structural stereo resolver. Planning reads `#T` / `#C` assertions from the
//! materialized aromaticity state and emits stereo-element additions plus
//! optional source-constraint removals without mutating the molecule.

use thiserror::Error;
use umol_ast::ast::{
    AsLit, AtomConstraintAst, AtomHandle, AtomId, AtomUpdate, BondConstraintAst, BondHandle,
    BondId, BondUpdate, CisTransStereoAst, Edit, MoleculeAst, StereoAtomAst, StereoBondAst,
    StereoCosetAst, StereoKind, StereoLigand, StereoLigandKind, TetrahedralStereoAst,
    TransactionError,
};
use umol_utils::solution::Solution;

use crate::ops::model::{InconsistencyPolicy, StereoModel};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StereoResolveConfig {
    pub reset_stereo_constraints: bool,
}

#[derive(Clone, Debug)]
pub struct StereoResolver {
    model: StereoModel,
    config: StereoResolveConfig,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoContradiction {
    #[error("tetrahedral stereo assertion at atom {0:?} cannot be realized")]
    UnrealizableAtom(AtomId),
    #[error("cis-trans stereo assertion at bond {0:?} cannot be realized")]
    UnrealizableBond(BondId),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoError {
    #[error(transparent)]
    Transaction(#[from] TransactionError),
}

impl StereoResolver {
    pub fn new(model: &StereoModel) -> Self {
        Self::with_config(model, StereoResolveConfig::default())
    }

    pub fn with_config(model: &StereoModel, config: StereoResolveConfig) -> Self {
        Self {
            model: model.clone(),
            config,
        }
    }

    /// Construct the complete stereo edit plan without mutating `ast`.
    pub fn plan(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<Vec<Edit>, StereoContradiction>, StereoError> {
        let mut edits = Vec::new();
        for id in ast.atoms().ids() {
            if ast.stereo_atoms().is_at(id) {
                continue;
            }
            let TetrahedralStereoAst::Stereo(coset) = ast
                .atom(id)
                .ast
                .constraints
                .tetrahedral_stereo()
                .unwrap_or(&TetrahedralStereoAst::Undetermined)
            else {
                continue;
            };
            let Some((ligands, stereo)) = self.derive_stereo_atom(ast, id, coset.clone()) else {
                match self.model.inconsistency {
                    InconsistencyPolicy::Keep => {}
                    InconsistencyPolicy::Strip => {
                        let mut update = AtomUpdate::default();
                        update.constraints.set(AtomConstraintAst::TetrahedralStereo(
                            TetrahedralStereoAst::Undetermined,
                        ));
                        edits.extend(Edit::for_atom_update(
                            AtomHandle::Id(id),
                            ast.atom(id).ast,
                            &update,
                        ));
                    }
                    InconsistencyPolicy::Error => {
                        return Ok(Solution::Contradictory(
                            StereoContradiction::UnrealizableAtom(id),
                        ));
                    }
                }
                continue;
            };
            edits.push(Edit::AddStereoAtom {
                site: AtomHandle::Id(id),
                ligands: ligands
                    .into_iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                ast: stereo,
            });
            if self.config.reset_stereo_constraints {
                let mut update = AtomUpdate::default();
                update.constraints.set(AtomConstraintAst::TetrahedralStereo(
                    TetrahedralStereoAst::Undetermined,
                ));
                edits.extend(Edit::for_atom_update(
                    AtomHandle::Id(id),
                    ast.atom(id).ast,
                    &update,
                ));
            }
        }
        for id in ast.bonds().ids() {
            if ast.stereo_bonds().is_at(id) {
                continue;
            }
            let CisTransStereoAst::Stereo(coset) = ast
                .bond(id)
                .ast
                .constraints
                .cis_trans_stereo()
                .unwrap_or(&CisTransStereoAst::Undetermined)
            else {
                continue;
            };
            let Some((ligands, stereo)) = self.derive_stereo_bond(ast, id, coset.clone()) else {
                match self.model.inconsistency {
                    InconsistencyPolicy::Keep => {}
                    InconsistencyPolicy::Strip => {
                        let mut update = BondUpdate::default();
                        update.constraints.set(BondConstraintAst::CisTransStereo(
                            CisTransStereoAst::Undetermined,
                        ));
                        edits.extend(Edit::for_bond_update(
                            BondHandle::Id(id),
                            ast.bond(id).ast,
                            &update,
                        ));
                    }
                    InconsistencyPolicy::Error => {
                        return Ok(Solution::Contradictory(
                            StereoContradiction::UnrealizableBond(id),
                        ));
                    }
                }
                continue;
            };
            edits.push(Edit::AddStereoBond {
                site: BondHandle::Id(id),
                ligands: ligands
                    .into_iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                ast: stereo,
            });
            if self.config.reset_stereo_constraints {
                let mut update = BondUpdate::default();
                update.constraints.set(BondConstraintAst::CisTransStereo(
                    CisTransStereoAst::Undetermined,
                ));
                edits.extend(Edit::for_bond_update(
                    BondHandle::Id(id),
                    ast.bond(id).ast,
                    &update,
                ));
            }
        }
        Ok(Solution::Determined(edits))
    }

    /// Plan and atomically apply structural stereo resolution.
    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), StereoContradiction>, StereoError> {
        let edits = match self.plan(ast)? {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let mut editor = ast.edit();
        editor.transact(edits)?;
        *ast = editor.build();
        Ok(Solution::Determined(()))
    }

    fn derive_stereo_atom(
        &self,
        ast: &MoleculeAst,
        id: AtomId,
        coset: StereoCosetAst,
    ) -> Option<(Vec<StereoLigand>, StereoAtomAst)> {
        let atom = ast.atom(id);
        if atom.is_in_aromatic_system() {
            return None;
        }

        let kind = StereoKind::Tetrahedral;
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

        Some((ligands, StereoAtomAst::new(kind, coset)))
    }

    fn derive_stereo_bond(
        &self,
        ast: &MoleculeAst,
        id: BondId,
        coset: StereoCosetAst,
    ) -> Option<(Vec<StereoLigand>, StereoBondAst)> {
        let bond = ast.bond(id);

        let kind = StereoKind::CisTrans;
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

        Some((ligands, StereoBondAst::new(kind, coset)))
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
    use rstest::{fixture, rstest};
    use umol_ast::ast::{AtomId, BondId};
    use umol_ast::mol_dsl_ground;

    use super::*;

    #[fixture]
    fn stereo_model() -> StereoModel {
        StereoModel::default()
    }

    #[rstest]
    fn test_stereo_resolve_config_default() {
        assert_eq!(
            StereoResolveConfig::default(),
            StereoResolveConfig {
                reset_stereo_constraints: false,
            }
        );
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
                             :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#),
        vec![Edit::AddStereoAtom {
            site: AtomHandle::Id(AtomId(1)),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(1)), StereoLigandKind::ImplicitHydrogen),
            ],
            ast: StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
        }]
    )]
    #[case::cis_trans(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]}"#),
        vec![Edit::AddStereoBond {
            site: BondHandle::Id(BondId(1)),
            ligands: vec![
                (AtomHandle::Id(AtomId(0)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(1)), StereoLigandKind::ImplicitHydrogen),
                (AtomHandle::Id(AtomId(3)), StereoLigandKind::Atom),
                (AtomHandle::Id(AtomId(2)), StereoLigandKind::ImplicitHydrogen),
            ],
            ast: StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
        }]
    )]
    fn test_stereo_resolver_plan(
        stereo_model: StereoModel,
        #[case] molecule: MoleculeAst,
        #[case] expected: Vec<Edit>,
    ) {
        assert_eq!(
            StereoResolver::new(&stereo_model).plan(&molecule),
            Ok(Solution::Determined(expected))
        );
    }

    #[rstest]
    #[case::no_assertion(mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h3"] :bonds [[0 1 "1"]]}"#))]
    #[case::existing_element(mol_dsl_ground!(r#"{
        :atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
        :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :type "Th1"}]
    }"#))]
    fn test_stereo_resolver_plan_identity(
        stereo_model: StereoModel,
        #[case] molecule: MoleculeAst,
    ) {
        assert_eq!(
            StereoResolver::new(&stereo_model).plan(&molecule),
            Ok(Solution::Determined(Vec::new()))
        );
    }

    #[rstest]
    #[case::keep(InconsistencyPolicy::Keep, Solution::Determined(Vec::new()))]
    #[case::strip(
        InconsistencyPolicy::Strip,
        Solution::Determined(vec![Edit::ModifyAtomConstraint {
            id: AtomHandle::Id(AtomId(1)),
            old: Some(AtomConstraintAst::TetrahedralStereo(
                TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1)),
            )),
            new: None,
        }])
    )]
    #[case::error(
        InconsistencyPolicy::Error,
        Solution::Contradictory(StereoContradiction::UnrealizableAtom(AtomId(1)))
    )]
    fn test_stereo_resolver_plan_inconsistency(
        mut stereo_model: StereoModel,
        #[case] policy: InconsistencyPolicy,
        #[case] expected: Solution<Vec<Edit>, StereoContradiction>,
    ) {
        stereo_model.inconsistency = policy;
        let molecule = mol_dsl_ground!(
            r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"] :bonds [[0 1 "1"] [1 2 "1"]]}"#
        );
        assert_eq!(
            StereoResolver::new(&stereo_model).plan(&molecule),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::tetrahedral(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1 #T1" "N #h2" "O #h1"]
                             :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]}"#),
        mol_dsl_ground!(r#"{
            :atoms ["C #h3" "C #h1" "N #h2" "O #h1"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :type "Th1"}]
        }"#)
    )]
    #[case::cis_trans(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]}"#),
        mol_dsl_ground!(r#"{
            :atoms ["C #h3" "C #h1" "C #h1" "C #h3"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :type "Ct1"}]
        }"#)
    )]
    fn test_stereo_resolver_resolve(
        stereo_model: StereoModel,
        #[case] mut molecule: MoleculeAst,
        #[case] expected: MoleculeAst,
    ) {
        let resolver = StereoResolver::with_config(
            &stereo_model,
            StereoResolveConfig {
                reset_stereo_constraints: true,
            },
        );
        assert_eq!(
            resolver.resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule, expected);
    }

    #[rstest]
    #[case::atom(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "S #h0 #T1" "C #h3"]
                             :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        StereoContradiction::UnrealizableAtom(AtomId(1))
    )]
    #[case::bond(
        mol_dsl_ground!(r#"{:atoms ["C #h3" "C #h2" "C #h1"]
                             :bonds [[0 1 "1"] [1 2 "2#C1"]]}"#),
        StereoContradiction::UnrealizableBond(BondId(1))
    )]
    fn test_stereo_resolver_resolve_error(
        stereo_model: StereoModel,
        #[case] mut molecule: MoleculeAst,
        #[case] expected: StereoContradiction,
    ) {
        let original = molecule.clone();
        assert_eq!(
            StereoResolver::new(&stereo_model).resolve(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
    }
}
