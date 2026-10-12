//! Policy-free structural stereo perception.
//!
//! [`StereoPerception`] derives tetrahedral atom and cis-trans bond relations
//! from complete `#T` and `#C` evidence and checks existing relations.
//! Resolver inconsistency policy and edit construction remain separate.

use std::collections::BTreeSet;

use thiserror::Error;
use umol_graph_ir::ir::{
    AsLit, AtomConstraintForm, AtomConstraintKey, AtomId, BondConstraintForm, BondConstraintKey,
    BondId, CisTransStereoForm, Molecule, RingConfig, RingModel, StereoAtomForm, StereoAtomId,
    StereoBondForm, StereoBondId, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    TetrahedralStereoForm,
};

use crate::ops::model::StereoModel;

/// Policy-free result of stereo perception and relation comparison.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StereoDerivation {
    /// Realizable tetrahedral stereo atoms.
    pub stereo_atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomForm)>,
    /// Realizable cis-trans stereo bonds.
    pub stereo_bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondForm)>,
    /// Sites of cis-trans assertions on double bonds whose smallest ring is below the model's
    /// `stereo_bond_minimum_ring_size`: no stereo bond is realized, and nothing is inconsistent.
    pub skipped_stereo_bonds: Vec<BondId>,
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

    // TODO: fold this back into `derive` once the resolver accepts undetermined inputs and the
    // partial-constraint gate in `StereoResolver::plan` is gone; it is separate only for that gate.
    /// Sites of cis-trans assertions the model does not realize: double bonds whose smallest
    /// ring is below `stereo_bond_minimum_ring_size`. Empty when the threshold is zero or no bond
    /// carries an assertion, without enumerating rings.
    pub fn skipped_stereo_bonds(&self, molecule: &Molecule) -> BTreeSet<BondId> {
        let minimum_ring_size = self.model.stereo_bond_minimum_ring_size as usize;
        let asserted: Vec<BondId> = molecule
            .bonds()
            .iter()
            .filter(|bond| {
                matches!(
                    bond.constraints().cis_trans_stereo(),
                    CisTransStereoForm::Stereo(_)
                )
            })
            .map(|bond| bond.id())
            .collect();
        if minimum_ring_size == 0 || asserted.is_empty() {
            return BTreeSet::new();
        }
        let rings = molecule.rings(RingModel::default(), RingConfig::default());
        asserted
            .into_iter()
            .filter(|&bond| {
                rings
                    .bond(bond)
                    .smallest_ring_size()
                    .is_some_and(|size| size < minimum_ring_size)
            })
            .collect()
    }

    /// Derive stereo relations from complete evidence and check existing relations.
    /// Positive evidence supplies the candidate configuration; negative evidence supplies none.
    /// Contradictory evidence between realizable sources is a mismatch. Independent constraint
    /// and entity realizability failures retain their separate diagnostics.
    pub fn derive(&self, molecule: &Molecule) -> StereoDerivation {
        let mut stereo_atoms = Vec::new();
        let mut stereo_bonds = Vec::new();
        let mut inconsistencies = BTreeSet::new();
        let skipped = self.skipped_stereo_bonds(molecule);

        for view in molecule.atoms().iter() {
            let atom = view.id();
            let relation = molecule.stereo_atoms().at(atom);
            let constraints = view.constraints();
            let asserted = constraints.tetrahedral_stereo();
            if relation.is_none() && !matches!(asserted, TetrahedralStereoForm::Stereo(_)) {
                continue;
            }
            let evidence = constraints.evidenced_complete(AtomConstraintKey::TetrahedralStereo);
            let coset = match &evidence {
                Some(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(
                    coset,
                ))) => Some(coset),
                _ => None,
            };
            let candidate = if coset.is_some() || relation.is_some() {
                self.derive_stereo_atom(molecule, atom, coset.unwrap_or(&StereoCoset::Undetermined))
            } else {
                None
            };
            if candidate.is_none()
                && (matches!(asserted, TetrahedralStereoForm::Stereo(_))
                    || evidence.is_none() && relation.is_none())
            {
                inconsistencies.insert(StereoInconsistency::TetrahedralStereoFailure { atom });
            }
            if let Some(relation) = relation {
                let realizable = relation.attributes().configuration.kind()
                    == Some(StereoKind::Tetrahedral)
                    && candidate.as_ref().is_some_and(|(ligands, _)| {
                        // Evidence already transports the configuration; tetrahedral permits
                        // every permutation of the same distinct ligands.
                        (evidence.is_some()
                            && relation.ligand_ids().len() == ligands.len()
                            && relation
                                .ligand_ids()
                                .iter()
                                .all(|ligand| ligands.contains(ligand)))
                            || relation.coset_for(ligands.iter().copied()).is_some()
                    });
                if !realizable {
                    inconsistencies.insert(StereoInconsistency::StereoAtomFailure {
                        stereo_atom: relation.id(),
                    });
                    continue;
                } else if evidence.is_none() {
                    inconsistencies.insert(StereoInconsistency::TetrahedralStereoMismatch {
                        atom,
                        stereo_atom: relation.id(),
                    });
                }
            }
            if coset.is_some() {
                if let Some((ligands, stereo)) = candidate {
                    stereo_atoms.push((atom, ligands, stereo));
                }
            }
        }

        for view in molecule.bonds().iter() {
            let bond = view.id();
            let relation = molecule.stereo_bonds().at(bond);
            let constraints = view.constraints();
            let asserted = constraints.cis_trans_stereo();
            if relation.is_none() && !matches!(asserted, CisTransStereoForm::Stereo(_)) {
                continue;
            }
            let evidence = constraints.evidenced_complete(BondConstraintKey::CisTransStereo);
            let coset = match &evidence {
                Some(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(coset))) => {
                    Some(coset)
                }
                _ => None,
            };
            let candidate = if coset.is_some() && !skipped.contains(&bond) || relation.is_some() {
                self.derive_stereo_bond(molecule, bond, coset.unwrap_or(&StereoCoset::Undetermined))
            } else {
                None
            };
            if candidate.is_none()
                && !skipped.contains(&bond)
                && (matches!(asserted, CisTransStereoForm::Stereo(_))
                    || evidence.is_none() && relation.is_none())
            {
                inconsistencies.insert(StereoInconsistency::CisTransStereoFailure { bond });
            }
            if let Some(relation) = relation {
                let realizable = relation.attributes().configuration.kind()
                    == Some(StereoKind::CisTrans)
                    && candidate.as_ref().is_some_and(|(ligands, _)| {
                        (evidence.is_some() && relation.ligand_ids() == ligands.as_slice())
                            || relation.coset_for(ligands.iter().copied()).is_some()
                    });
                if !realizable {
                    inconsistencies.insert(StereoInconsistency::StereoBondFailure {
                        stereo_bond: relation.id(),
                    });
                    continue;
                } else if evidence.is_none() {
                    inconsistencies.insert(StereoInconsistency::CisTransStereoMismatch {
                        bond,
                        stereo_bond: relation.id(),
                    });
                }
            }
            if coset.is_some() && !skipped.contains(&bond) {
                if let Some((ligands, stereo)) = candidate {
                    stereo_bonds.push((bond, ligands, stereo));
                }
            }
        }

        StereoDerivation {
            stereo_atoms,
            stereo_bonds,
            skipped_stereo_bonds: skipped.into_iter().collect(),
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
    use umol_graph_ir::mol_dsl_concrete;

    use super::*;
    use crate::ops::model::{ElementScope, StereoKindModel};

    #[rstest]
    #[case::materializations(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C1"] [6 7 "1"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            stereo_bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            skipped_stereo_bonds: vec![],
            inconsistencies: vec![],
        },
    )]
    #[case::existing_elements(
        mol_dsl_concrete!(r#"{
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
            stereo_atoms: vec![(
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            stereo_bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            skipped_stereo_bonds: vec![],
            inconsistencies: vec![],
        },
    )]
    #[case::unrealizable_sites(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "S#h0#T1" "C#h3"
                    "C#h3" "C#h2" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [3 4 "1"] [4 5 "2#C1"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],
            skipped_stereo_bonds: vec![],
            inconsistencies: vec![
                StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) },
                StereoInconsistency::CisTransStereoFailure { bond: BondId(3) },
            ],
        },
    )]
    #[case::relation_mismatches(
        mol_dsl_concrete!(r#"{
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
            stereo_atoms: vec![],
            stereo_bonds: vec![],
            skipped_stereo_bonds: vec![],
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
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "N#h2" "O#h"
                    "C#h3" "C#h0" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2"] [6 7 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Sp1"}]
            :stereo-bonds [{:site 4 :ligands [4 [:h 5] 7 [:h 6]] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],
            skipped_stereo_bonds: vec![],
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
    #[case::virtual_ligand_failures(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h" "F" "Cl" "Br" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [4 5 "2"] [4 6 "1"] [5 7 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:lp 0]] :attrs "Th1"}]
            :stereo-bonds [{:site 3 :ligands [6 [:lp 4] 7 [:h 5]] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            inconsistencies: vec![
                StereoInconsistency::StereoAtomFailure { stereo_atom: StereoAtomId(0) },
                StereoInconsistency::StereoBondFailure { stereo_bond: StereoBondId(0) },
            ],
            ..StereoDerivation::default()
        },
    )]
    #[case::negative_evidence(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#T!" "F" "Cl" "Br" "C#h" "C#h" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [4 5 "2#C!"] [4 6 "1"] [5 7 "1"]]
        }"#),
        StereoDerivation::default(),
    )]
    #[case::not_stereo_mismatches(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T!" "N#h2" "O#h"
                    "C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]
                    [4 5 "1"] [5 6 "2#C!"] [6 7 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 2 3 [:h 1]] :attrs "Th1"}]
            :stereo-bonds [{:site 4 :ligands [4 [:h 5] 7 [:h 6]] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],
            skipped_stereo_bonds: vec![],
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
    #[case::cyclohexene_asserted(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h1" "C#h1" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 5 "1"] [0 1 "2#C1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],
            skipped_stereo_bonds: vec![BondId(1)],
            inconsistencies: vec![],
        },
    )]
    #[case::cyclohexene_asserted_undetermined(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h1" "C#h1" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 5 "1"] [0 1 "2#C+"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![],
            stereo_bonds: vec![],
            skipped_stereo_bonds: vec![BondId(1)],
            inconsistencies: vec![],
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
    #[case::tetrahedral_refinement(
        mol_dsl_concrete!(r#"{
            :atoms ["C#T+" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]
        }"#),
        StereoDerivation {
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            ..StereoDerivation::default()
        },
    )]
    #[case::cis_trans_refinement(
        mol_dsl_concrete!(r#"{
            :atoms ["C" "C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "2#C+"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1"}]
        }"#),
        StereoDerivation {
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..StereoDerivation::default()
        },
    )]
    fn test_stereo_perception_derive_evidence(
        #[case] molecule: Molecule,
        #[case] expected: StereoDerivation,
    ) {
        assert_eq!(
            StereoPerception::new(&StereoModel::default()).derive(&molecule),
            expected,
        );
    }

    #[rstest]
    #[case::tetrahedral(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "S#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        None,
    )]
    #[case::ligand_arity(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h0#T1" "N#h2"]
            :bonds [[0 1 "1"] [1 2 "1"]]
        }"#),
        None,
    )]
    #[case::aromatic_exclusion(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "N#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
        }"#),
        None,
    )]
    #[case::ligand_arity(
        StereoModel::default(),
        mol_dsl_concrete!(r#"{
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
