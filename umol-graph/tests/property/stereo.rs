//! Recognition and classification on finite tetrahedral and cis/trans ligand domains.
//!
//! Four explicit terminal ligands are labeled H, F, Cl, Br, or I. The independent
//! reference uses label equality: four distinct Th ligands are stereogenic, one
//! equal pair is prochiral, and other partitions are symmetric. Ct is stereogenic
//! exactly when each endpoint's two labels differ. Recognition selects the exact
//! expected frame and open coset; classification is independent of the supplied
//! binary coset domain. Both queries preserve their input.

use std::collections::{BTreeMap, BTreeSet};
use std::iter;

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph::ops::model::StereoModel;
use umol_graph::ops::stereo::{StereoDerivation, StereoPerception, StereoSelection, StereoSite};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, BondId, Molecule, MoleculeEntries, RingConfig, StereoAtomForm,
    StereoBondForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind, Stereogenicity,
};

proptest! {
    #[test]
    fn test_stereo_perception_recognize_tetrahedral(
        labels in prop::array::uniform4(0usize..5), classes in 0u8..8,
        coset in prop_oneof![Just(StereoCoset::Undetermined), Just(StereoCoset::Lit(0)),
            Just(StereoCoset::Lit(1)), Just(StereoCoset::lit_set([0, 1]))],
    ) {
        let distinct = labels.into_iter().collect::<BTreeSet<_>>().len();
        let class = match distinct {
            4 => Stereogenicity::Stereogenic,
            3 => Stereogenicity::Prochiral,
            _ => Stereogenicity::Symmetric,
        };
        let selected: BTreeSet<_> = [Stereogenicity::Symmetric, Stereogenicity::Prochiral, Stereogenicity::Stereogenic]
            .into_iter().enumerate().filter_map(|(i, class)| (classes & (1 << i) != 0).then_some(class)).collect();
        let elements = [Element::H, Element::F, Element::Cl, Element::Br, Element::I];
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: iter::once(AtomForm::from_element(Element::C).into_concrete())
                .chain(labels.map(|label| AtomForm::from_element(elements[label]).into_concrete())).collect(),
            bonds: (1..=4).map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1).into_concrete())).collect(),
            ..Default::default()
        });
        let original = molecule.clone();
        let ligands: Vec<_> = (1..=4).map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect();
        let perception = StereoPerception::new(&StereoModel::default());
        let expected = StereoDerivation {
            stereo_atoms: if selected.contains(&class) {
                vec![(AtomId(0), ligands.clone(), StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined))]
            } else { vec![] }, ..Default::default()
        };
        prop_assert_eq!(perception.recognize(&molecule, &StereoSelection {
            sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
            kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])), stereogenicities: Some(selected),
        }, AutomorphismAlgorithm::Nauty, RingConfig::default()), Ok(expected));
        prop_assert_eq!(&molecule, &original);

        let mut editor = molecule.edit();
        editor.add_stereo_atom(AtomId(0), &ligands, StereoAtomForm::new(StereoKind::Tetrahedral, coset));
        let stored = editor.finish().unwrap();
        let original = stored.clone();
        prop_assert_eq!(perception.classify(&stored, &BTreeSet::from([StereoSite::Atom(AtomId(0))]),
            AutomorphismAlgorithm::Nauty), Ok(BTreeMap::from([(StereoSite::Atom(AtomId(0)), class)])));
        prop_assert_eq!(stored, original);
    }

    #[test]
    fn test_stereo_perception_recognize_cis_trans(
        labels in prop::array::uniform4(0usize..5), classes in 0u8..8,
        coset in prop_oneof![Just(StereoCoset::Undetermined), Just(StereoCoset::Lit(0)),
            Just(StereoCoset::Lit(1)), Just(StereoCoset::lit_set([0, 1]))],
    ) {
        let class = if labels[0] != labels[1] && labels[2] != labels[3] {
            Stereogenicity::Stereogenic
        } else { Stereogenicity::Symmetric };
        let selected: BTreeSet<_> = [Stereogenicity::Symmetric, Stereogenicity::Prochiral, Stereogenicity::Stereogenic]
            .into_iter().enumerate().filter_map(|(i, class)| (classes & (1 << i) != 0).then_some(class)).collect();
        let elements = [Element::H, Element::F, Element::Cl, Element::Br, Element::I];
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: [AtomForm::from_element(Element::C).into_concrete(), AtomForm::from_element(Element::C).into_concrete()]
                .into_iter().chain(labels.map(|label| AtomForm::from_element(elements[label]).into_concrete())).collect(),
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2).into_concrete()),
                (AtomId(0), AtomId(2), BondForm::from_order(1).into_concrete()),
                (AtomId(0), AtomId(3), BondForm::from_order(1).into_concrete()),
                (AtomId(1), AtomId(4), BondForm::from_order(1).into_concrete()),
                (AtomId(1), AtomId(5), BondForm::from_order(1).into_concrete())],
            ..Default::default()
        });
        let original = molecule.clone();
        let ligands: Vec<_> = (2..=5).map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect();
        let perception = StereoPerception::new(&StereoModel::default());
        let expected = StereoDerivation {
            stereo_bonds: if selected.contains(&class) {
                vec![(BondId(0), ligands.clone(), StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined))]
            } else { vec![] }, ..Default::default()
        };
        prop_assert_eq!(perception.recognize(&molecule, &StereoSelection {
            sites: Some(BTreeSet::from([StereoSite::Bond(BondId(0))])),
            kinds: Some(BTreeSet::from([StereoKind::CisTrans])), stereogenicities: Some(selected),
        }, AutomorphismAlgorithm::Nauty, RingConfig::default()), Ok(expected));
        prop_assert_eq!(&molecule, &original);

        let mut editor = molecule.edit();
        editor.add_stereo_bond(BondId(0), &ligands, StereoBondForm::new(StereoKind::CisTrans, coset));
        let stored = editor.finish().unwrap();
        let original = stored.clone();
        prop_assert_eq!(perception.classify(&stored, &BTreeSet::from([StereoSite::Bond(BondId(0))]),
            AutomorphismAlgorithm::Nauty), Ok(BTreeMap::from([(StereoSite::Bond(BondId(0)), class)])));
        prop_assert_eq!(stored, original);
    }
}
