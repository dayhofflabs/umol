//! Policy-free structural stereo perception.
//!
//! [`StereoPerception`] derives tetrahedral atom and cis-trans bond relations
//! from `#T` and `#C` projections and compares them with existing relations.
//! Resolver inconsistency policy and edit construction remain separate.

use std::collections::BTreeSet;

use thiserror::Error;
use umol_graph_ir::ir::{
    AsLit, AtomId, BondId, CisTransStereoForm, Lattice, Molecule, StereoAtomForm, StereoAtomId,
    StereoBondForm, StereoBondId, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    TetrahedralStereoForm,
};

use crate::ops::model::StereoModel;

/// Policy-free result of stereo perception and relation comparison.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StereoDerivation {
    /// Realizable tetrahedral stereo atoms.
    pub atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>,
    /// Realizable cis-trans stereo bonds.
    pub bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>,
    /// Constraint failures, entity failures, and independently valid mismatches.
    pub inconsistencies: Vec<StereoInconsistency>,
}

/// Policy-free classification of stereo constraint and entity inconsistencies.
#[derive(Clone, Copy, Debug, Error, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoInconsistency {
    #[error("tetrahedral stereo constraint at atom {atom:?} cannot be realized")]
    TetrahedralStereoFailure { atom: AtomId },
    #[error("stereo atom {stereo_atom:?} cannot be realized")]
    StereoAtomFailure { stereo_atom: StereoAtomId },
    #[error(
        "tetrahedral stereo constraint at atom {atom:?} disagrees with stereo atom {stereo_atom:?}"
    )]
    TetrahedralStereoMismatch {
        atom: AtomId,
        stereo_atom: StereoAtomId,
    },
    #[error("cis-trans stereo constraint at bond {bond:?} cannot be realized")]
    CisTransStereoFailure { bond: BondId },
    #[error("stereo bond {stereo_bond:?} cannot be realized")]
    StereoBondFailure { stereo_bond: StereoBondId },
    #[error(
        "cis-trans stereo constraint at bond {bond:?} disagrees with stereo bond {stereo_bond:?}"
    )]
    CisTransStereoMismatch {
        bond: BondId,
        stereo_bond: StereoBondId,
    },
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
    pub fn derive(&self, molecule: &Molecule) -> StereoDerivation {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut inconsistencies = BTreeSet::new();

        for atom in molecule.atoms().ids() {
            let relations: Vec<_> = molecule
                .stereo_atoms()
                .iter()
                .filter(|relation| relation.site_id() == atom)
                .collect();
            let assertion = molecule
                .atom(atom)
                .attributes
                .constraints
                .tetrahedral_stereo()
                .unwrap_or(&TetrahedralStereoForm::Undetermined);

            let constraint_candidate = match assertion {
                TetrahedralStereoForm::Stereo(coset) => {
                    let candidate = self.derive_stereo_atom(molecule, atom, coset);
                    if candidate.is_none() {
                        inconsistencies
                            .insert(StereoInconsistency::TetrahedralStereoFailure { atom });
                    }
                    candidate
                }
                TetrahedralStereoForm::NotStereo => None,
                TetrahedralStereoForm::Undetermined => None,
            };

            if let Some((ligands, stereo)) = &constraint_candidate {
                atoms.push((atom, ligands.clone(), stereo.clone()));
            }
            for relation in relations {
                let entity_candidate =
                    if relation.attributes.configuration.kind() == Some(StereoKind::Tetrahedral) {
                        self.derive_stereo_atom(molecule, atom, &StereoCoset::Undetermined)
                            .and_then(|(ligands, _)| {
                                let coset = relation.coset_for(ligands.iter().copied())?;
                                Some((ligands, StereoAtomForm::new(StereoKind::Tetrahedral, coset)))
                            })
                    } else {
                        None
                    };
                let Some((entity_ligands, entity_stereo)) = entity_candidate else {
                    inconsistencies.insert(StereoInconsistency::StereoAtomFailure {
                        stereo_atom: relation.id,
                    });
                    continue;
                };

                if constraint_candidate.is_none()
                    && matches!(assertion, TetrahedralStereoForm::Undetermined)
                {
                    atoms.push((atom, entity_ligands.clone(), entity_stereo.clone()));
                }
                let mismatch =
                    match assertion {
                        TetrahedralStereoForm::Stereo(_) => constraint_candidate
                            .as_ref()
                            .is_some_and(|(_, constraint_stereo)| {
                                let Some(expected) = constraint_stereo.configuration.coset() else {
                                    return false;
                                };
                                let Some(actual) = entity_stereo.configuration.coset() else {
                                    return false;
                                };
                                !TetrahedralStereoForm::Stereo(expected.clone())
                                    .matches(&TetrahedralStereoForm::Stereo(actual.clone()))
                            }),
                        TetrahedralStereoForm::NotStereo => true,
                        TetrahedralStereoForm::Undetermined => false,
                    };
                if mismatch {
                    inconsistencies.insert(StereoInconsistency::TetrahedralStereoMismatch {
                        atom,
                        stereo_atom: relation.id,
                    });
                }
            }
        }

        for bond in molecule.bonds().ids() {
            let relations: Vec<_> = molecule
                .stereo_bonds()
                .iter()
                .filter(|relation| relation.site_id() == bond)
                .collect();
            let assertion = molecule
                .bond(bond)
                .attributes
                .constraints
                .cis_trans_stereo()
                .unwrap_or(&CisTransStereoForm::Undetermined);

            let constraint_candidate = match assertion {
                CisTransStereoForm::Stereo(coset) => {
                    let candidate = self.derive_stereo_bond(molecule, bond, coset);
                    if candidate.is_none() {
                        inconsistencies.insert(StereoInconsistency::CisTransStereoFailure { bond });
                    }
                    candidate
                }
                CisTransStereoForm::NotStereo => None,
                CisTransStereoForm::Undetermined => None,
            };

            if let Some((ligands, stereo)) = &constraint_candidate {
                bonds.push((bond, ligands.clone(), stereo.clone()));
            }
            for relation in relations {
                let entity_candidate =
                    if relation.attributes.configuration.kind() == Some(StereoKind::CisTrans) {
                        self.derive_stereo_bond(molecule, bond, &StereoCoset::Undetermined)
                            .and_then(|(ligands, _)| {
                                let coset = relation.coset_for(ligands.iter().copied())?;
                                Some((ligands, StereoBondForm::new(StereoKind::CisTrans, coset)))
                            })
                    } else {
                        None
                    };
                let Some((entity_ligands, entity_stereo)) = entity_candidate else {
                    inconsistencies.insert(StereoInconsistency::StereoBondFailure {
                        stereo_bond: relation.id,
                    });
                    continue;
                };

                if constraint_candidate.is_none()
                    && matches!(assertion, CisTransStereoForm::Undetermined)
                {
                    bonds.push((bond, entity_ligands.clone(), entity_stereo.clone()));
                }
                let mismatch = match assertion {
                    CisTransStereoForm::Stereo(_) => {
                        constraint_candidate
                            .as_ref()
                            .is_some_and(|(_, constraint_stereo)| {
                                let Some(expected) = constraint_stereo.configuration.coset() else {
                                    return false;
                                };
                                let Some(actual) = entity_stereo.configuration.coset() else {
                                    return false;
                                };
                                !CisTransStereoForm::Stereo(expected.clone())
                                    .matches(&CisTransStereoForm::Stereo(actual.clone()))
                            })
                    }
                    CisTransStereoForm::NotStereo => true,
                    CisTransStereoForm::Undetermined => false,
                };
                if mismatch {
                    inconsistencies.insert(StereoInconsistency::CisTransStereoMismatch {
                        bond,
                        stereo_bond: relation.id,
                    });
                }
            }
        }

        StereoDerivation {
            atoms,
            bonds,
            inconsistencies: inconsistencies.into_iter().collect(),
        }
    }

    /// Derive a tetrahedral stereo atom in the canonical neighbor frame.
    pub fn derive_stereo_atom(
        &self,
        molecule: &Molecule,
        atom: AtomId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoAtomForm)> {
        let view = molecule.atom(atom);
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

        Some((ligands, StereoAtomForm::new(kind, coset.clone())))
    }

    /// Derive a cis-trans stereo bond in the canonical endpoint and side frames.
    pub fn derive_stereo_bond(
        &self,
        molecule: &Molecule,
        bond: BondId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoBondForm)> {
        let view = molecule.bond(bond);
        let kind = StereoKind::CisTrans;
        let model = self.model.kind_model(kind)?;
        let [first, second] = view.atom_ids();
        if !model
            .scope
            .contains(molecule.atom(first).element().as_lit()?)
            || !model
                .scope
                .contains(molecule.atom(second).element().as_lit()?)
        {
            return None;
        }

        let first_side = self.bond_side_ligands(molecule, first, second)?;
        let second_side = self.bond_side_ligands(molecule, second, first)?;
        let ligands = vec![first_side[0], first_side[1], second_side[0], second_side[1]];
        Some((ligands, StereoBondForm::new(kind, coset.clone())))
    }

    fn bond_side_ligands(
        &self,
        molecule: &Molecule,
        atom: AtomId,
        partner: AtomId,
    ) -> Option<[StereoLigand; 2]> {
        let view = molecule.atom(atom);
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
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{StereoAtomId, StereoBondId};
    use umol_graph_ir::mol_dsl_ground;

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
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            inconsistencies: vec![],
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
                :attrs "Th1"
            }]
            :stereo-bonds [{
                :site 4
                :ligands [4 [:h 5] 7 [:h 6]]
                :attrs "Ct1"
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
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            inconsistencies: vec![],
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
            inconsistencies: vec![
                StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) },
                StereoInconsistency::CisTransStereoFailure { bond: BondId(3) },
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
                :attrs "Th0"
            }]
            :stereo-bonds [{
                :site 4
                :ligands [4 [:h 5] 7 [:h 6]]
                :attrs "Ct0"
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
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            inconsistencies: vec![
                StereoInconsistency::TetrahedralStereoMismatch {
                    atom: AtomId(1),
                    stereo_atom: StereoAtomId(0),
                },
                StereoInconsistency::CisTransStereoMismatch {
                    bond: BondId(4),
                    stereo_bond: StereoBondId(0),
                },
            ],
        },
    )]
    #[case::entity_failures(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h" "N#h2" "O#h"
                    "C#h3" "C#h0" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2"] [6 7 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Sp1"}]
            :stereo-bonds [{:site 4 :ligands [4 [:h 5] 7 [:h 6]] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            atoms: vec![],
            bonds: vec![],
            inconsistencies: vec![
                StereoInconsistency::StereoAtomFailure {
                    stereo_atom: StereoAtomId(0),
                },
                StereoInconsistency::StereoBondFailure {
                    stereo_bond: StereoBondId(0),
                },
            ],
        },
    )]
    #[case::not_stereo_mismatches(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T!" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C!"] [6 7 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th1"}]
            :stereo-bonds [{:site 4 :ligands [4 [:h 5] 7 [:h 6]] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            atoms: vec![],
            bonds: vec![],
            inconsistencies: vec![
                StereoInconsistency::TetrahedralStereoMismatch {
                    atom: AtomId(1),
                    stereo_atom: StereoAtomId(0),
                },
                StereoInconsistency::CisTransStereoMismatch {
                    bond: BondId(4),
                    stereo_bond: StereoBondId(0),
                },
            ],
        },
    )]
    fn test_stereo_perception_derive(
        #[case] molecule: Molecule,
        #[case] expected: StereoDerivation,
    ) {
        assert_eq!(
            StereoPerception::new(&StereoModel::default()).derive(&molecule),
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
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
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
            :aromatic-systems [{:atoms [0 1 2] :attrs "[1,1,1]"}]
        }"#),
        None,
    )]
    fn test_stereo_perception_derive_stereo_atom(
        #[case] model: StereoModel,
        #[case] molecule: Molecule,
        #[case] expected: Option<(Vec<StereoLigand>, StereoAtomForm)>,
    ) {
        assert_eq!(
            StereoPerception::new(&model).derive_stereo_atom(
                &molecule,
                AtomId(1),
                &StereoCoset::Lit(1),
            ),
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
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
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
        #[case] molecule: Molecule,
        #[case] expected: Option<(Vec<StereoLigand>, StereoBondForm)>,
    ) {
        assert_eq!(
            StereoPerception::new(&model).derive_stereo_bond(
                &molecule,
                BondId(1),
                &StereoCoset::Lit(1),
            ),
            expected
        );
    }
}
