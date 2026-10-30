use std::collections::{BTreeMap, BTreeSet};

use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomFieldKind, AtomForm, AtomId, BondConstraintForm, BondForm, BondId,
    CisTransStereoForm, ElementForm, Molecule, MoleculeEntries, NumForm, RingConfig,
    StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId, StereoConfigurationForm,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind, StereoTerm, Stereogenicity,
    TetrahedralStereoForm,
};
use umol_graph_ir::mol_dsl_concrete;

use super::StereoPerceptionError;
use crate::ops::model::{ElementScope, StereoKindModel, StereoModel};
use crate::ops::stereo::{
    StereoDerivation, StereoInconsistency, StereoPerception, StereoSelection, StereoSite,
};

#[fixture]
fn mixed_sites() -> Molecule {
    mol_dsl_concrete!(
        r#"{:atoms ["C" "H" "H" "F" "Cl" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2"] [5 7 "1"] [6 8 "1"]]}"#
    )
}

#[rstest]
#[case::tetrahedral(
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
    StereoSite::Atom(AtomId(0)),
    Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
    ]),
)]
#[case::implicit_hydrogen(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
    StereoSite::Atom(AtomId(0)),
    Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
    ]),
)]
#[case::phosphine(
    mol_dsl_concrete!(r#"{:atoms ["P#n" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
    StereoSite::Atom(AtomId(0)),
    Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
    ]),
)]
#[case::phosphonium(
    mol_dsl_concrete!(r#"{:atoms ["P#c+" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
    ]),
)]
#[case::ammonium(
    mol_dsl_concrete!(r#"{:atoms ["N#c+" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
    ]),
)]
#[case::arsine_hydrogen(
    mol_dsl_concrete!(r#"{:atoms ["As#h#n" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
    StereoSite::Atom(AtomId(0)),
    Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
    ]),
)]
#[case::sulfoxide(
    mol_dsl_concrete!(r#"{:atoms ["S#n" "O#n2" "C#h3" "C#h2" "C#h3"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [3 4 "1"]]}"#),
    StereoSite::Atom(AtomId(0)),
    Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
    ]),
)]
#[case::selenonium(
    mol_dsl_concrete!(r#"{:atoms ["Se#c+#n" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
    StereoSite::Atom(AtomId(0)),
    Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
    ]),
)]
#[case::amine(
    mol_dsl_concrete!(r#"{:atoms ["N#n" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), None,
)]
#[case::implicit_hydrogens(
    mol_dsl_concrete!(r#"{:atoms ["C#h2" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), None,
)]
#[case::lone_pairs(
    mol_dsl_concrete!(r#"{:atoms ["S#n2" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), None,
)]
#[case::carbonyl(
    mol_dsl_concrete!(r#"{:atoms ["C" "O#n2" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), None,
)]
#[case::oxime(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "N#n" "F" "O#h#n2"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]}"#),
    StereoSite::Bond(BondId(0)), Some(vec![
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
    ]),
)]
#[case::aziridine(
    mol_dsl_concrete!(r#"{:atoms ["N#n" "C#h2" "C#h2" "F"]
        :bonds [[0 1 "1"] [0 2 "1"] [1 2 "1"] [0 3 "1"]]}"#),
    StereoSite::Atom(AtomId(0)),
    Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
    ]),
)]
#[case::conjugated_nitrogen(
    mol_dsl_concrete!(r#"{:atoms ["N#n" "C#h2" "C#h2" "C#h" "O#n2"]
        :bonds [[0 1 "1"] [0 2 "1"] [1 2 "1"] [0 3 "1"] [3 4 "2"]]}"#),
    StereoSite::Atom(AtomId(0)), None,
)]
#[case::vinyl_nitrogen(
    mol_dsl_concrete!(r#"{:atoms ["N#n" "C#h2" "C#h2" "C#h" "C#h2"]
        :bonds [[0 1 "1"] [0 2 "1"] [1 2 "1"] [0 3 "1"] [3 4 "2"]]}"#),
    StereoSite::Atom(AtomId(0)), None,
)]
#[case::phosphoryl_nitrogen(
    mol_dsl_concrete!(r#"{:atoms ["N#n" "C#h2" "C#h2" "P" "O#n2" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [1 2 "1"] [0 3 "1"] [3 4 "2"] [3 5 "1"] [3 6 "1"]]}"#),
    StereoSite::Atom(AtomId(0)), Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
    ]),
)]
#[case::bridgehead_nitrogen(
    mol_dsl_concrete!(r#"{:atoms ["N#n" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h2" "C#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 7 "1"] [0 3 "1"] [3 4 "1"] [4 7 "1"]
                [0 5 "1"] [5 6 "1"] [6 7 "1"]]}"#),
    StereoSite::Atom(AtomId(0)),
    Some(vec![
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
    ]),
)]
#[case::cis_trans(
    mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#),
    StereoSite::Bond(BondId(0)),
    Some(vec![
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
    ]),
)]
#[case::imine(
    mol_dsl_concrete!(r#"{:atoms ["C#h" "N#n" "F" "C#h3"]
        :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]}"#),
    StereoSite::Bond(BondId(0)),
    Some(vec![
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
    ]),
)]
fn test_stereo_perception_recognize(
    #[case] molecule: Molecule,
    #[case] site: StereoSite,
    #[case] ligands: Option<Vec<StereoLigand>>,
) {
    let original = molecule.clone();
    let selection = StereoSelection {
        sites: Some(BTreeSet::from([site])),
        ..Default::default()
    };
    let mut expected = StereoDerivation::default();
    if let Some(ligands) = ligands {
        match site {
            StereoSite::Atom(atom) => expected.stereo_atoms.push((
                atom,
                ligands,
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            )),
            StereoSite::Bond(bond) => expected.stereo_bonds.push((
                bond,
                ligands,
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            )),
        }
    }
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).recognize(
            &molecule,
            &selection,
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(expected)
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::unrestricted(StereoSelection::default(), vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(4))])]
#[case::site(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(4))])), ..Default::default() }, vec![StereoSite::Bond(BondId(4))])]
#[case::tetrahedral(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])), ..Default::default() }, vec![StereoSite::Atom(AtomId(0))])]
#[case::cis_trans(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..Default::default() }, vec![StereoSite::Bond(BondId(4))])]
#[case::stereogenic(StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])), ..Default::default() }, vec![StereoSite::Bond(BondId(4))])]
#[case::prochiral(StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Prochiral])), ..Default::default() }, vec![StereoSite::Atom(AtomId(0))])]
#[case::symmetric(StereoSelection { stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])), ..Default::default() }, vec![])]
#[case::intersection(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])), kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..Default::default() }, vec![])]
#[case::empty_sites(StereoSelection { sites: Some(BTreeSet::new()), ..Default::default() }, vec![])]
#[case::empty_kinds(StereoSelection { kinds: Some(BTreeSet::new()), ..Default::default() }, vec![])]
#[case::empty_classes(StereoSelection { stereogenicities: Some(BTreeSet::new()), ..Default::default() }, vec![])]
fn test_stereo_perception_recognize_selection(
    mixed_sites: Molecule,
    #[case] selection: StereoSelection,
    #[case] sites: Vec<StereoSite>,
) {
    let original = mixed_sites.clone();
    let mut expected = StereoDerivation::default();
    for site in sites {
        match site {
            StereoSite::Atom(atom) => expected.stereo_atoms.push((
                atom,
                (1..=4)
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            )),
            StereoSite::Bond(bond) => expected.stereo_bonds.push((
                bond,
                vec![
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(8), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            )),
        }
    }
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).recognize(
            &mixed_sites,
            &selection,
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(expected)
    );
    assert_eq!(mixed_sites, original);
}

#[rstest]
#[case::tetrahedral(StereoSite::Atom(AtomId(0)))]
#[case::cis_trans(StereoSite::Bond(BondId(4)))]
fn test_stereo_perception_recognize_evidence(
    mut mixed_sites: Molecule,
    #[case] site: StereoSite,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0), StereoCoset::Lit(1), StereoCoset::lit_set([0, 1]))]
    coset: StereoCoset,
) {
    let mut expected = StereoDerivation::default();
    match site {
        StereoSite::Atom(atom) => {
            mixed_sites.atom_mut(atom).attributes_mut().constraints =
                AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::Stereo(
                    coset.clone(),
                ))
                .into();
            expected.stereo_atoms.push((
                atom,
                (1..=4)
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, coset),
            ));
        }
        StereoSite::Bond(bond) => {
            mixed_sites.bond_mut(bond).attributes_mut().constraints =
                BondConstraintForm::cis_trans_stereo(CisTransStereoForm::Stereo(coset.clone()))
                    .into();
            expected.stereo_bonds.push((
                bond,
                vec![
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(8), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, coset),
            ));
        }
    }
    let original = mixed_sites.clone();
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).recognize(
            &mixed_sites,
            &StereoSelection {
                sites: Some(BTreeSet::from([site])),
                ..Default::default()
            },
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(expected)
    );
    assert_eq!(mixed_sites, original);
}

#[rstest]
#[case::literal(StereoCoset::Lit(1))]
#[case::symbolic(StereoCoset::term(StereoTerm::var("configuration")))]
fn test_stereo_perception_recognize_assertion(#[case] coset: StereoCoset) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["N#n" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#
    );
    molecule.atom_mut(AtomId(0)).attributes_mut().constraints =
        AtomConstraintForm::tetrahedral_stereo(TetrahedralStereoForm::Stereo(coset.clone())).into();
    let original = molecule.clone();
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).recognize(
            &molecule,
            &StereoSelection::default(),
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(StereoDerivation {
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, coset)
            )],
            ..Default::default()
        })
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::negative_atom(mol_dsl_concrete!(r#"{:atoms ["C#T!" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
#[case::negative_bond(mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "2#C!"] [0 2 "1"] [1 3 "1"]]}"#))]
#[case::existing_atom(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [2 1 3 4] :attrs "Th1"}]}"#))]
#[case::existing_bond(mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
    :stereo-bonds [{:site 0 :ligands [[:h 0] 2 3 [:h 1]] :attrs "Ct1"}]}"#))]
fn test_stereo_perception_recognize_identity(#[case] molecule: Molecule) {
    let original = molecule.clone();
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).recognize(
            &molecule,
            &StereoSelection::default(),
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(StereoDerivation::default())
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::atom(mol_dsl_concrete!(r#"{:atoms ["C#T0" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#),
    StereoInconsistency::TetrahedralStereoMismatch { atom: AtomId(0), stereo_atom: StereoAtomId(0) })]
#[case::competing_kind(mol_dsl_concrete!(r#"{:atoms ["C#T0" "F" "Cl" "Br" "I"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]}"#),
    StereoInconsistency::TetrahedralStereoMismatch { atom: AtomId(0), stereo_atom: StereoAtomId(0) })]
#[case::bond(mol_dsl_concrete!(r#"{:atoms ["C#h" "C#h" "F" "Cl"]
    :bonds [[0 1 "2#C0"] [0 2 "1"] [1 3 "1"]]
    :stereo-bonds [{:site 0 :ligands [2 [:h 0] 3 [:h 1]] :attrs "Ct1"}]}"#),
    StereoInconsistency::CisTransStereoMismatch { bond: BondId(0), stereo_bond: StereoBondId(0) })]
#[case::atom_frame(mol_dsl_concrete!(r#"{:atoms ["C#h2#T1" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
    StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(0) })]
#[case::bond_frame(mol_dsl_concrete!(r#"{:atoms ["C#h2" "C#h" "Cl"] :bonds [[0 1 "2#C1"] [1 2 "1"]]}"#),
    StereoInconsistency::CisTransStereoFailure { bond: BondId(0) })]
fn test_stereo_perception_recognize_inconsistency(
    #[case] molecule: Molecule,
    #[case] inconsistency: StereoInconsistency,
) {
    let original = molecule.clone();
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).recognize(
            &molecule,
            &StereoSelection::default(),
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(StereoDerivation {
            inconsistencies: vec![inconsistency],
            ..Default::default()
        })
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::four(4, 8, false)]
#[case::six(6, 8, false)]
#[case::eight(8, 8, true)]
#[case::twelve(12, 8, true)]
#[case::unrestricted(6, 0, true)]
#[case::long_ring(24, 25, false)]
fn test_stereo_perception_recognize_ring(
    #[case] size: u32,
    #[case] minimum: u32,
    #[case] eligible: bool,
    #[values(false, true)] asserted: bool,
) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: (0..size)
            .map(|id| {
                AtomForm::from_element(Element::C)
                    .with_implicit_hydrogens(if id < 2 { 1 } else { 2 })
                    .into_concrete()
            })
            .collect(),
        bonds: (0..size)
            .map(|id| {
                (
                    AtomId(id),
                    AtomId((id + 1) % size),
                    if id == 0 {
                        let form = BondForm::from_order(2).into_concrete();
                        if asserted {
                            form.with_constraint(BondConstraintForm::cis_trans_stereo(
                                CisTransStereoForm::Stereo(StereoCoset::Undetermined),
                            ))
                        } else {
                            form
                        }
                    } else {
                        BondForm::from_order(1).into_concrete()
                    },
                )
            })
            .collect(),
        ..Default::default()
    });
    let model = StereoModel {
        stereo_bond_minimum_ring_size: minimum,
        ..Default::default()
    };
    let expected = StereoDerivation {
        stereo_bonds: if eligible {
            vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(size - 1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined),
            )]
        } else {
            vec![]
        },
        skipped_stereo_bonds: if asserted && !eligible {
            vec![BondId(0)]
        } else {
            vec![]
        },
        ..Default::default()
    };
    let original = molecule.clone();
    assert_eq!(
        StereoPerception::new(&model).recognize(
            &molecule,
            &StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Bond(BondId(0))])),
                ..Default::default()
            },
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(expected)
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::disabled(None)]
#[case::excluded(Some(StereoKindModel { scope: ElementScope::AllowList(vec![Element::N]), fluxionality: false }))]
fn test_stereo_perception_recognize_model(
    mixed_sites: Molecule,
    #[case] kind_model: Option<StereoKindModel>,
    #[values(StereoKind::Tetrahedral, StereoKind::CisTrans)] kind: StereoKind,
) {
    let mut model = StereoModel::default();
    model.kind_models[kind as usize] = kind_model;
    let selection = StereoSelection {
        kinds: Some(BTreeSet::from([kind])),
        ..Default::default()
    };
    assert_eq!(
        StereoPerception::new(&model).recognize(
            &mixed_sites,
            &selection,
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(StereoDerivation::default())
    );
}

#[rstest]
#[case::invalid_atom(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(9))])), ..Default::default() },
    StereoPerceptionError::InvalidSite { site: StereoSite::Atom(AtomId(9)) })]
#[case::invalid_bond(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Bond(BondId(7))])), ..Default::default() },
    StereoPerceptionError::InvalidSite { site: StereoSite::Bond(BondId(7)) })]
#[case::invalid_empty(StereoSelection { sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(9))])), kinds: Some(BTreeSet::new()), ..Default::default() },
    StereoPerceptionError::InvalidSite { site: StereoSite::Atom(AtomId(9)) })]
#[case::axial(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Axial])), ..Default::default() },
    StereoPerceptionError::UnsupportedKind { kind: StereoKind::Axial })]
#[case::square_planar(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::SquarePlanar])), ..Default::default() },
    StereoPerceptionError::UnsupportedKind { kind: StereoKind::SquarePlanar })]
#[case::trigonal_bipyramidal(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::TrigonalBipyramidal])), ..Default::default() },
    StereoPerceptionError::UnsupportedKind { kind: StereoKind::TrigonalBipyramidal })]
#[case::octahedral(StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Octahedral])), ..Default::default() },
    StereoPerceptionError::UnsupportedKind { kind: StereoKind::Octahedral })]
fn test_stereo_perception_recognize_error(
    mixed_sites: Molecule,
    #[case] selection: StereoSelection,
    #[case] error: StereoPerceptionError,
) {
    let original = mixed_sites.clone();
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).recognize(
            &mixed_sites,
            &selection,
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Err(error)
    );
    assert_eq!(mixed_sites, original);
}

#[rstest]
#[case::element(AtomFieldKind::Element)]
#[case::charge(AtomFieldKind::Charge)]
#[case::hydrogens(AtomFieldKind::ImplicitHydrogens)]
#[case::lone_pairs(AtomFieldKind::LonePairs)]
fn test_stereo_perception_recognize_partial(
    mut mixed_sites: Molecule,
    #[case] field: AtomFieldKind,
) {
    let mut atom = mixed_sites.atom_mut(AtomId(0));
    let attributes = atom.attributes_mut();
    match field {
        AtomFieldKind::Element => attributes.element = ElementForm::Undetermined,
        AtomFieldKind::Charge => attributes.charge = NumForm::Undetermined,
        AtomFieldKind::ImplicitHydrogens => attributes.implicit_hydrogens = NumForm::Undetermined,
        AtomFieldKind::LonePairs => attributes.lone_pairs = NumForm::Undetermined,
        _ => unreachable!(),
    }
    let perception = StereoPerception::new(&StereoModel::default());
    let original = mixed_sites.clone();
    assert_eq!(
        perception.recognize(
            &mixed_sites,
            &StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                ..Default::default()
            },
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Err(StereoPerceptionError::UnderdeterminedAtomField {
            atom: AtomId(0),
            field
        })
    );
    assert_eq!(
        perception.recognize(
            &mixed_sites,
            &StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Bond(BondId(4))])),
                ..Default::default()
            },
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(StereoDerivation {
            stereo_bonds: vec![(
                BondId(4),
                vec![
                    StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(8), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(6), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined)
            )],
            ..Default::default()
        })
    );
    assert_eq!(mixed_sites, original);
}

#[rstest]
#[case::atom(0, StereoSite::Atom(AtomId(0)))]
#[case::bond(4, StereoSite::Bond(BondId(4)))]
fn test_stereo_perception_recognize_bond_order(
    mut mixed_sites: Molecule,
    #[case] bond: u32,
    #[case] site: StereoSite,
) {
    mixed_sites.bond_mut(BondId(bond)).attributes_mut().order = NumForm::Undetermined;
    let original = mixed_sites.clone();
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).recognize(
            &mixed_sites,
            &StereoSelection {
                sites: Some(BTreeSet::from([site])),
                ..Default::default()
            },
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Err(StereoPerceptionError::UnderdeterminedBondOrder { bond: BondId(bond) })
    );
    assert_eq!(mixed_sites, original);
}

#[rstest]
#[case::tetrahedral(StereoKind::Tetrahedral)]
#[case::cis_trans(StereoKind::CisTrans)]
#[case::square_planar(StereoKind::SquarePlanar)]
fn test_stereo_perception_recognize_unsupported(mixed_sites: Molecule, #[case] kind: StereoKind) {
    let mut model = StereoModel::default();
    model.kind_models[kind as usize] = Some(StereoKindModel {
        scope: ElementScope::Any,
        fluxionality: true,
    });
    let expected = if kind == StereoKind::SquarePlanar {
        StereoPerceptionError::UnsupportedKind { kind }
    } else {
        StereoPerceptionError::UnsupportedFluxionality { kind }
    };
    assert_eq!(
        StereoPerception::new(&model).recognize(
            &mixed_sites,
            &StereoSelection::default(),
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Err(expected)
    );
}

#[rstest]
#[case::hydrogen_count(AtomFieldKind::ImplicitHydrogens, i64::MAX, true)]
#[case::charge(AtomFieldKind::Charge, i64::MIN, false)]
#[case::electron_count(AtomFieldKind::UnpairedElectrons, i64::MAX, true)]
fn test_stereo_perception_recognize_conjugation(
    #[case] field: AtomFieldKind,
    #[case] value: i64,
    #[case] eligible: bool,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{:atoms ["N#n" "C#h2" "C#h2" "C#h" "O#n2"]
        :bonds [[0 1 "1"] [0 2 "1"] [1 2 "1"] [0 3 "1"] [3 4 "2"]]}"#
    );
    let mut atom = molecule.atom_mut(AtomId(3));
    let attributes = atom.attributes_mut();
    match field {
        AtomFieldKind::ImplicitHydrogens => attributes.implicit_hydrogens = NumForm::Lit(value),
        AtomFieldKind::Charge => attributes.charge = NumForm::Lit(value),
        AtomFieldKind::UnpairedElectrons => {
            attributes.unpaired_electrons.count = NumForm::Lit(value)
        }
        _ => unreachable!(),
    }
    let original = molecule.clone();
    let expected = StereoDerivation {
        stereo_atoms: if eligible {
            vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined),
            )]
        } else {
            vec![]
        },
        ..Default::default()
    };
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).recognize(
            &molecule,
            &StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                ..Default::default()
            },
            AutomorphismAlgorithm::Nauty,
            RingConfig::default()
        ),
        Ok(expected)
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::atom(StereoSite::Atom(AtomId(0)), Stereogenicity::Stereogenic)]
#[case::bond(StereoSite::Bond(BondId(4)), Stereogenicity::Stereogenic)]
fn test_stereo_perception_classify(
    #[case] site: StereoSite,
    #[case] expected: Stereogenicity,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0), StereoCoset::Lit(1), StereoCoset::lit_set([0, 1]))]
    coset: StereoCoset,
) {
    let mut molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "F" "Cl" "Br" "I" "C#h" "C#h" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [5 6 "2"] [5 7 "1"] [6 8 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
        :stereo-bonds [{:site 4 :ligands [7 [:h 5] 8 [:h 6]] :attrs "Ct0"}]
    }"#
    );
    match site {
        StereoSite::Atom(_) => {
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset)
        }
        StereoSite::Bond(_) => {
            molecule
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .configuration = StereoConfigurationForm::Kinded(StereoKind::CisTrans, coset)
        }
    }
    let original = molecule.clone();
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).classify(
            &molecule,
            &BTreeSet::from([site]),
            AutomorphismAlgorithm::Nauty
        ),
        Ok(BTreeMap::from([(site, expected)]))
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::invalid_atom(Molecule::default(), StereoSite::Atom(AtomId(0)),
    StereoPerceptionError::InvalidSite { site: StereoSite::Atom(AtomId(0)) })]
#[case::invalid_bond(Molecule::default(), StereoSite::Bond(BondId(0)),
    StereoPerceptionError::InvalidSite { site: StereoSite::Bond(BondId(0)) })]
#[case::absent_entity(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#), StereoSite::Atom(AtomId(0)),
    StereoPerceptionError::MissingStereoEntity { site: StereoSite::Atom(AtomId(0)) })]
#[case::unknown_kind(
    mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "*"}]}"#),
    StereoSite::Atom(AtomId(0)),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) },
)]
fn test_stereo_perception_classify_error(
    #[case] molecule: Molecule,
    #[case] site: StereoSite,
    #[case] expected: StereoPerceptionError,
) {
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).classify(
            &molecule,
            &BTreeSet::from([site]),
            AutomorphismAlgorithm::Nauty
        ),
        Err(expected)
    );
}

#[rstest]
#[case::atom(mol_dsl_concrete!(r#"{:atoms ["C#h" "F" "Cl" "Br"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3] :attrs "Th0"}]}"#), StereoSite::Atom(AtomId(0)),
    StereoInconsistency::StereoAtomFailure { stereo_atom: StereoAtomId(0) })]
#[case::five_ligands(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4 5] :attrs "Th0"}]}"#), StereoSite::Atom(AtomId(0)),
    StereoInconsistency::StereoAtomFailure { stereo_atom: StereoAtomId(0) })]
fn test_stereo_perception_classify_frame(
    #[case] molecule: Molecule,
    #[case] site: StereoSite,
    #[case] inconsistency: StereoInconsistency,
) {
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).classify(
            &molecule,
            &BTreeSet::from([site]),
            AutomorphismAlgorithm::Nauty
        ),
        Err(StereoPerceptionError::Inconsistency(inconsistency))
    );
}

#[rstest]
#[case::distinct(mol_dsl_concrete!(r#"{:atoms ["N#n" "F" "Cl" "Br"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 [:lp 0]] :attrs "Th0"}]}"#), Stereogenicity::Stereogenic)]
#[case::equivalent_hydrogens(mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#), Stereogenicity::Prochiral)]
#[case::symmetric(mol_dsl_concrete!(r#"{:atoms ["C" "H" "H" "H" "H"]
    :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
    :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]}"#), Stereogenicity::Symmetric)]
fn test_stereo_perception_classify_model(
    #[case] mut molecule: Molecule,
    #[case] expected: Stereogenicity,
    #[values(StereoCoset::Undetermined, StereoCoset::Lit(0), StereoCoset::Lit(1), StereoCoset::lit_set([0, 1]))]
    coset: StereoCoset,
) {
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset);
    let original = molecule.clone();
    let model = StereoModel {
        kind_models: [None, None, None, None, None, None],
        ..Default::default()
    };
    assert_eq!(
        StereoPerception::new(&model).classify(
            &molecule,
            &BTreeSet::from([StereoSite::Atom(AtomId(0))]),
            AutomorphismAlgorithm::Nauty
        ),
        Ok(BTreeMap::from([(StereoSite::Atom(AtomId(0)), expected)]))
    );
    assert_eq!(molecule, original);
}

#[rstest]
#[case::unknown_kind(StereoConfigurationForm::Undetermined,
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::symbolic(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::var("a"))),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::empty(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::lit_set([])),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::invalid_coset(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(2)),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::invalid_set(StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::lit_set([0, 2])),
    StereoPerceptionError::UnsupportedConfiguration { site: StereoSite::Atom(AtomId(0)) })]
#[case::geometry(StereoConfigurationForm::Kinded(StereoKind::SquarePlanar, StereoCoset::Undetermined),
    StereoPerceptionError::UnsupportedKind { kind: StereoKind::SquarePlanar })]
fn test_stereo_perception_classify_configuration(
    mut mixed_sites: Molecule,
    #[case] configuration: StereoConfigurationForm,
    #[case] expected: StereoPerceptionError,
) {
    let mut editor = mixed_sites.edit();
    editor.add_stereo_atom(
        AtomId(0),
        &[1, 2, 3, 4].map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)),
        StereoAtomForm::default(),
    );
    mixed_sites = editor.finish().unwrap();
    mixed_sites
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = configuration;
    let original = mixed_sites.clone();
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).classify(
            &mixed_sites,
            &BTreeSet::from([StereoSite::Atom(AtomId(0))]),
            AutomorphismAlgorithm::Nauty
        ),
        Err(expected)
    );
    assert_eq!(mixed_sites, original);
}

#[rstest]
fn test_stereo_perception_classify_empty(mixed_sites: Molecule) {
    assert_eq!(
        StereoPerception::new(&StereoModel::default()).classify(
            &mixed_sites,
            &BTreeSet::new(),
            AutomorphismAlgorithm::Nauty
        ),
        Ok(BTreeMap::new())
    );
}
