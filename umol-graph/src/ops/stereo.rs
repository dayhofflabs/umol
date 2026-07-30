//! Policy-free structural stereo perception.
//!
//! [`StereoPerception`] derives tetrahedral atom and cis-trans bond relations
//! from `#T` and `#C` projections and compares them with existing relations.
//! Resolver inconsistency policy and edit construction remain separate.

use std::collections::BTreeSet;

use umol_ast::ast::{
    AsLit, AtomId, BondId, CisTransStereoAst, Lattice, MoleculeAst, StereoAtomAst, StereoAtomId,
    StereoBondAst, StereoBondId, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    TetrahedralStereoAst,
};

use crate::ops::model::StereoModel;

/// Policy-free result of stereo perception and relation comparison.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StereoDerivation {
    /// Realizable tetrahedral stereo atoms.
    pub atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)>,
    /// Realizable cis-trans stereo bonds.
    pub bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)>,
    /// Unrealizable projections and existing relations that disagree with perception.
    pub mismatches: Vec<StereoMismatch>,
}

/// A stereo projection or relation that disagrees with perception.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoMismatch {
    UnrealizableAtom { atom: AtomId },
    UnrealizableBond { bond: BondId },
    AtomRelation { stereo_atom: StereoAtomId },
    BondRelation { stereo_bond: StereoBondId },
}

/// Structural stereo perception under a selected chemistry model.
#[derive(Clone, Debug)]
pub struct StereoPerception {
    model: StereoModel,
}

impl StereoPerception {
    pub fn new(model: &StereoModel) -> Self {
        Self {
            model: model.clone(),
        }
    }

    /// Derive stereo relations and compare them with every existing relation.
    pub fn derive(&self, ast: &MoleculeAst) -> StereoDerivation {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut mismatches = BTreeSet::new();

        for atom in ast.atoms().ids() {
            let relations: Vec<_> = ast
                .stereo_atoms()
                .iter()
                .filter(|relation| relation.site_id() == atom)
                .collect();
            let assertion = ast
                .atom(atom)
                .ast
                .constraints
                .tetrahedral_stereo()
                .unwrap_or(&TetrahedralStereoAst::Undetermined);

            let candidate = match assertion {
                TetrahedralStereoAst::Stereo(coset) => {
                    let candidate = self.derive_stereo_atom(ast, atom, coset);
                    if candidate.is_none() {
                        mismatches.insert(StereoMismatch::UnrealizableAtom { atom });
                    }
                    candidate
                }
                TetrahedralStereoAst::NotStereo => None,
                TetrahedralStereoAst::Undetermined => relations.first().and_then(|relation| {
                    if relation.ast.configuration.kind() != Some(StereoKind::Tetrahedral) {
                        return None;
                    }
                    let coset = relation.ast.configuration.coset()?;
                    let (ligands, _) = self.derive_stereo_atom(ast, atom, coset)?;
                    let coset = relation.coset_for(ligands.iter().copied())?;
                    self.derive_stereo_atom(ast, atom, &coset)
                }),
            };

            if let Some((ligands, stereo)) = &candidate {
                atoms.push((atom, ligands.clone(), stereo.clone()));
            }
            for relation in relations {
                let matches = candidate.as_ref().is_some_and(|(ligands, stereo)| {
                    if relation.ast.configuration.kind() != Some(StereoKind::Tetrahedral) {
                        return false;
                    }
                    let Some(expected) = stereo.configuration.coset() else {
                        return false;
                    };
                    relation
                        .coset_for(ligands.iter().copied())
                        .is_some_and(|actual| {
                            TetrahedralStereoAst::Stereo(expected.clone())
                                .matches(&TetrahedralStereoAst::Stereo(actual))
                        })
                });
                if !matches {
                    mismatches.insert(StereoMismatch::AtomRelation {
                        stereo_atom: relation.id,
                    });
                }
            }
        }

        for bond in ast.bonds().ids() {
            let relations: Vec<_> = ast
                .stereo_bonds()
                .iter()
                .filter(|relation| relation.site_id() == bond)
                .collect();
            let assertion = ast
                .bond(bond)
                .ast
                .constraints
                .cis_trans_stereo()
                .unwrap_or(&CisTransStereoAst::Undetermined);

            let candidate = match assertion {
                CisTransStereoAst::Stereo(coset) => {
                    let candidate = self.derive_stereo_bond(ast, bond, coset);
                    if candidate.is_none() {
                        mismatches.insert(StereoMismatch::UnrealizableBond { bond });
                    }
                    candidate
                }
                CisTransStereoAst::NotStereo => None,
                CisTransStereoAst::Undetermined => relations.first().and_then(|relation| {
                    if relation.ast.configuration.kind() != Some(StereoKind::CisTrans) {
                        return None;
                    }
                    let coset = relation.ast.configuration.coset()?;
                    let (ligands, _) = self.derive_stereo_bond(ast, bond, coset)?;
                    let coset = relation.coset_for(ligands.iter().copied())?;
                    self.derive_stereo_bond(ast, bond, &coset)
                }),
            };

            if let Some((ligands, stereo)) = &candidate {
                bonds.push((bond, ligands.clone(), stereo.clone()));
            }
            for relation in relations {
                let matches = candidate.as_ref().is_some_and(|(ligands, stereo)| {
                    if relation.ast.configuration.kind() != Some(StereoKind::CisTrans) {
                        return false;
                    }
                    let Some(expected) = stereo.configuration.coset() else {
                        return false;
                    };
                    relation
                        .coset_for(ligands.iter().copied())
                        .is_some_and(|actual| {
                            CisTransStereoAst::Stereo(expected.clone())
                                .matches(&CisTransStereoAst::Stereo(actual))
                        })
                });
                if !matches {
                    mismatches.insert(StereoMismatch::BondRelation {
                        stereo_bond: relation.id,
                    });
                }
            }
        }

        StereoDerivation {
            atoms,
            bonds,
            mismatches: mismatches.into_iter().collect(),
        }
    }

    /// Derive a tetrahedral stereo atom in the canonical neighbor frame.
    pub fn derive_stereo_atom(
        &self,
        ast: &MoleculeAst,
        atom: AtomId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoAtomAst)> {
        let view = ast.atom(atom);
        if view.is_in_aromatic_system() {
            return None;
        }

        let kind = StereoKind::Tetrahedral;
        let model = self.model.kind_model(kind)?;
        if !model.scope.contains(view.element().as_lit()?) {
            return None;
        }

        let mut ligands: Vec<StereoLigand> = view
            .neighbors()
            .map(|neighbor| StereoLigand::new(neighbor.atom_id(), StereoLigandKind::Atom))
            .collect();
        if ligands.len() + 1 == kind.degree() {
            let virtual_kind = if view.implicit_hydrogens().as_lit()? >= 1 {
                StereoLigandKind::ImplicitHydrogen
            } else if view.lone_pairs().as_lit()? >= 1 {
                StereoLigandKind::LonePair
            } else {
                return None;
            };
            ligands.push(StereoLigand::new(atom, virtual_kind));
        }
        if ligands.len() != kind.degree() {
            return None;
        }

        Some((ligands, StereoAtomAst::new(kind, coset.clone())))
    }

    /// Derive a cis-trans stereo bond in the canonical endpoint and side frames.
    pub fn derive_stereo_bond(
        &self,
        ast: &MoleculeAst,
        bond: BondId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoBondAst)> {
        let view = ast.bond(bond);
        let kind = StereoKind::CisTrans;
        let model = self.model.kind_model(kind)?;
        let [first, second] = view.atom_ids();
        if !model.scope.contains(ast.atom(first).element().as_lit()?)
            || !model.scope.contains(ast.atom(second).element().as_lit()?)
        {
            return None;
        }

        let first_side = self.bond_side_ligands(ast, first, second)?;
        let second_side = self.bond_side_ligands(ast, second, first)?;
        let ligands = vec![first_side[0], first_side[1], second_side[0], second_side[1]];
        Some((ligands, StereoBondAst::new(kind, coset.clone())))
    }

    fn bond_side_ligands(
        &self,
        ast: &MoleculeAst,
        atom: AtomId,
        partner: AtomId,
    ) -> Option<[StereoLigand; 2]> {
        let view = ast.atom(atom);
        let mut substituents = view
            .neighbors()
            .map(|neighbor| neighbor.atom_id())
            .filter(|&neighbor| neighbor != partner);
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
    use umol_ast::ast::{StereoAtomId, StereoBondId};
    use umol_ast::mol_dsl_ground;
    use umol_chem::element::Element;

    use super::*;
    use crate::ops::model::{ElementScope, StereoKindModel};

    #[rstest]
    #[case::materializations(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C1"] [6 7 "1"]]
        }"#),
        StereoDerivation {
            atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            mismatches: vec![],
        },
    )]
    #[case::existing_elements(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2"] [6 7 "1"]]
            :stereo-atoms [{
                :site 1
                :ligands [0 2 3 [:h 1]]
                :type "Th1"
            }]
            :stereo-bonds [{
                :site 4
                :ligands [4 [:h 5] 7 [:h 6]]
                :type "Ct1"
            }]
        }"#),
        StereoDerivation {
            atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            mismatches: vec![],
        },
    )]
    #[case::unrealizable_sites(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "S#h0#T1" "C#h3"
                    "C#h3" "C#h2" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [3 4 "1"] [4 5 "2#C1"]]
        }"#),
        StereoDerivation {
            atoms: vec![],
            bonds: vec![],
            mismatches: vec![
                StereoMismatch::UnrealizableAtom { atom: AtomId(1) },
                StereoMismatch::UnrealizableBond { bond: BondId(3) },
            ],
        },
    )]
    #[case::relation_mismatches(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C1"] [6 7 "1"]]
            :stereo-atoms [{
                :site 1
                :ligands [0 2 3 [:h 1]]
                :type "Th0"
            }]
            :stereo-bonds [{
                :site 4
                :ligands [4 [:h 5] 7 [:h 6]]
                :type "Ct0"
            }]
        }"#),
        StereoDerivation {
            atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            mismatches: vec![
                StereoMismatch::AtomRelation {
                    stereo_atom: StereoAtomId(0),
                },
                StereoMismatch::BondRelation {
                    stereo_bond: StereoBondId(0),
                },
            ],
        },
    )]
    fn test_stereo_perception_derive(#[case] ast: MoleculeAst, #[case] expected: StereoDerivation) {
        assert_eq!(
            StereoPerception::new(&StereoModel::default()).derive(&ast),
            expected
        );
    }

    #[rstest]
    #[case::tetrahedral(
        StereoModel::default(),
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        Some((
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        )),
    )]
    #[case::disabled_kind(
        {
            let mut model = StereoModel::default();
            model.kind_models[StereoKind::Tetrahedral as usize] = None;
            model
        },
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        None,
    )]
    #[case::element_scope(
        {
            let mut model = StereoModel::default();
            model.kind_models[StereoKind::Tetrahedral as usize] = Some(StereoKindModel {
                scope: ElementScope::AllowList(vec![Element::C]),
                fluxionality: false,
            });
            model
        },
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "S#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        None,
    )]
    #[case::ligand_arity(
        StereoModel::default(),
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h0#T1" "N#h2"]
            :bonds [[0 1 "1"] [1 2 "1"]]
        }"#),
        None,
    )]
    #[case::aromatic_exclusion(
        StereoModel::default(),
        mol_dsl_ground!(r#"{
            :atoms ["C" "C#T1" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [1 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2] :type "[1,1,1]"}]
        }"#),
        None,
    )]
    fn test_stereo_perception_derive_stereo_atom(
        #[case] model: StereoModel,
        #[case] ast: MoleculeAst,
        #[case] expected: Option<(Vec<StereoLigand>, StereoAtomAst)>,
    ) {
        assert_eq!(
            StereoPerception::new(&model)
                .derive_stereo_atom(&ast, AtomId(1), &StereoCoset::Lit(1),),
            expected
        );
    }

    #[rstest]
    #[case::cis_trans(
        StereoModel::default(),
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
        }"#),
        Some((
            vec![
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
            ],
            StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        )),
    )]
    #[case::disabled_kind(
        {
            let mut model = StereoModel::default();
            model.kind_models[StereoKind::CisTrans as usize] = None;
            model
        },
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
        }"#),
        None,
    )]
    #[case::element_scope(
        {
            let mut model = StereoModel::default();
            model.kind_models[StereoKind::CisTrans as usize] = Some(StereoKindModel {
                scope: ElementScope::AllowList(vec![Element::C]),
                fluxionality: false,
            });
            model
        },
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "N#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
        }"#),
        None,
    )]
    #[case::ligand_arity(
        StereoModel::default(),
        mol_dsl_ground!(r#"{
            :atoms ["C#h0#n0" "C#h" "C#h3"]
            :bonds [[0 1 "2#C1"] [1 2 "1"]]
        }"#),
        None,
    )]
    fn test_stereo_perception_derive_stereo_bond(
        #[case] model: StereoModel,
        #[case] ast: MoleculeAst,
        #[case] expected: Option<(Vec<StereoLigand>, StereoBondAst)>,
    ) {
        assert_eq!(
            StereoPerception::new(&model)
                .derive_stereo_bond(&ast, BondId(1), &StereoCoset::Lit(1),),
            expected
        );
    }
}
