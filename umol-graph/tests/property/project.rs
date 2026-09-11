//! Aromatic projection preserves member-aligned contributions and atom electron fields.
//! Independent five- and six-member ring expectations check exact projected values,
//! idempotence, and transport under rotated/reversed system frames. The same expectations
//! check projection of successfully ingested SMILES under both valence sources.
//! Stereo projection is checked against independent tetrahedral permutation parity and
//! cis-trans side-swap parity, preserving virtual-ligand kinds and open configurations.

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph::ingest::ingest_smiles_with;
use umol_graph::ops::model::{ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::{ResolveConfig, Resolver};
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticValenceForm, AtomConstraintForm, AtomForm, AtomId,
    BondConstraintForm, BondForm, BondId, BooleanForm, CisTransStereoForm, ElectronCountsForm,
    ElementForm, IsotopeMassForm, Molecule, MoleculeEntries, NumForm, StereoAtomForm,
    StereoBondForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind, TetrahedralStereoForm,
    UnpairedElectronsForm,
};
use umol_graph_ir::mol_dsl_concrete;
use umol_io::smiles::SmilesIoConfig;
use umol_utils::solution::Solution;

proptest! {
    #[test]
    fn test_aromaticity_resolver_project(
        rings in prop::collection::vec((0u8..5, 0usize..6, any::<bool>()), 1..5),
        typing in any::<bool>(),
    ) {
        let mut entries = MoleculeEntries::default();
        let mut expected = MoleculeEntries::default();
        let mut parsed_bonds = Vec::new();
        let mut smiles = Vec::new();
        for (kind, rotation, reverse) in rings {
            let size = if kind < 2 {6} else {5};
            smiles.push(match kind {
                0 => "c1ccccc1", 1 => "n1ccccc1", 2 => "[nH]1cccc1",
                3 => "o1cccc1", _ => "[cH-]1cccc1",
            });
            let start = entries.atoms.len();
            let mut contributions = Vec::new();
            for i in 0..size {
                let (element, charge, hydrogens, lone_pairs, contribution) = match (kind, i) {
                    (1, 0) => (Element::N, 0, 0, 1, 1),
                    (2, 0) => (Element::N, 0, 1, 0, 2),
                    (3, 0) => (Element::O, 0, 0, 1, 2),
                    (4, 0) => (Element::C, -1, 1, 0, 2),
                    _ => (Element::C, 0, 1, 0, 1),
                };
                entries.atoms.push(AtomForm {
                    element: ElementForm::Lit(element), isotope_mass: IsotopeMassForm::Natural,
                    charge: NumForm::Lit(charge), implicit_hydrogens: NumForm::Lit(hydrogens),
                    lone_pairs: NumForm::Lit(lone_pairs),
                    unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: Default::default(),
                });
                entries.bonds.push((AtomId((start + i) as u32), AtomId((start + (i + 1) % size) as u32), BondForm {
                    order: NumForm::Lit(1), charge: NumForm::Lit(0),
                    unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: Default::default(),
                }));
                expected.atoms.push(entries.atoms.last().unwrap().clone().with_constraint(
                    AtomConstraintForm::aromatic_valence(AromaticValenceForm::aromatic(contribution))));
                expected.bonds.push((AtomId((start + i) as u32), AtomId((start + (i + 1) % size) as u32), BondForm {
                    order: NumForm::Lit(1), charge: NumForm::Lit(0),
                    unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                    constraints: [BondConstraintForm::Aromatic(BooleanForm::Lit(true))].into_iter().collect(),
                }));
                contributions.push((AtomId((start + i) as u32), contribution));
            }
            let ring_bonds = &expected.bonds[expected.bonds.len() - size..];
            parsed_bonds.push(ring_bonds[size - 1].clone());
            parsed_bonds.extend_from_slice(&ring_bonds[..size - 1]);
            contributions.rotate_left(rotation % size);
            if reverse {contributions.reverse();}
            let (atoms, electrons) = contributions.into_iter().unzip();
            entries.aromatic.push((atoms, AromaticSystemForm {
                electrons: ElectronCountsForm::Lit(electrons), charge: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: Default::default(),
            }));
        }
        let model = ChemistryModel {valence: ValenceModel {
            tie_break: ValenceTieBreak::MostSaturated,
            ..if typing {ValenceModel::default()} else {ValenceModel::smiles()}
        }, ..Default::default()};
        let resolver = Resolver::new(&model);
        let mut molecule = Molecule::from_entries(entries);
        let expected_projection = Molecule::from_entries(expected.clone());
        prop_assert_eq!(resolver.aromaticity.project(&mut molecule), Ok(Solution::Determined(())));
        prop_assert_eq!(&molecule, &expected_projection);
        prop_assert_eq!(resolver.aromaticity.project(&mut molecule), Ok(Solution::Determined(())));
        prop_assert_eq!(molecule, expected_projection);

        let mut ingested = ingest_smiles_with(&smiles.join("."), &SmilesIoConfig::opensmiles(),
            &model, &ResolveConfig::default()).unwrap();
        prop_assert_eq!(resolver.aromaticity.project(&mut ingested), Ok(Solution::Determined(())));
        expected.bonds = parsed_bonds;
        prop_assert_eq!(ingested, Molecule::from_entries(expected));
    }

    #[test]
    fn test_stereo_resolver_project_tetrahedral(
        kind in 0u8..4,
        permutation in Just(vec![0usize, 1, 2, 3]).prop_shuffle(),
        coset in 0u32..2,
        open in any::<bool>(),
    ) {
        let base = match kind {
            0 => mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
            1 => mol_dsl_concrete!(r#"{:atoms ["C#h1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
            2 => mol_dsl_concrete!(r#"{:atoms ["N#n1" "F" "Cl" "Br"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#),
            _ => mol_dsl_concrete!(r#"{:atoms ["C" "F" "Cl" "Br" "H"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#),
        };
        let fixed = [
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            match kind {
                1 => StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                2 => StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                _ => StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            },
        ];
        let inversions = (0..4).flat_map(|i| (i + 1..4).map(move |j| (i, j)))
            .filter(|&(i, j)| permutation[i] > permutation[j]).count();
        let stored = if open { StereoCoset::Undetermined } else { StereoCoset::Lit(coset) };
        let projected = if open { StereoCoset::Undetermined } else { StereoCoset::Lit(coset ^ (inversions % 2) as u32) };
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: base.atoms().iter().map(|atom| atom.attributes.clone()).collect(),
            bonds: base.bonds().iter().map(|bond| {
                let [first, second] = bond.atom_ids();
                (first, second, bond.attributes.clone())
            }).collect(),
            stereo_atoms: vec![(AtomId(0), permutation.iter().map(|&i| fixed[i]).collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, stored))],
            ..Default::default()
        });
        let mut editor = base.edit();
        editor.atom_mut(AtomId(0)).attributes.constraints.set(
            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(projected)));
        let expected = editor.build();
        let model = ChemistryModel::default();
        let resolver = Resolver::new(&model);
        prop_assert_eq!(resolver.stereo.project(&mut molecule), Ok(Solution::Determined(())));
        prop_assert_eq!(&molecule, &expected);
        prop_assert_eq!(resolver.stereo.project(&mut molecule), Ok(Solution::Determined(())));
        prop_assert_eq!(molecule, expected);
    }

    #[test]
    fn test_stereo_resolver_project_cis_trans(
        kind in 0u8..3,
        first_swap in any::<bool>(), second_swap in any::<bool>(), endpoints_swap in any::<bool>(),
        coset in 0u32..2, open in any::<bool>(),
    ) {
        let (base, mut ligands) = match kind {
            0 => (mol_dsl_concrete!(r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"] :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]}"#),
                vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom)]),
            1 => (mol_dsl_concrete!(r#"{:atoms ["C#h1" "C#h1" "F" "Cl"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]}"#),
                vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen)]),
            _ => (mol_dsl_concrete!(r#"{:atoms ["N#n1" "C#h1" "C#h3" "F"] :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]}"#),
                vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen)]),
        };
        if first_swap { ligands.swap(0, 1); }
        if second_swap { ligands.swap(2, 3); }
        if endpoints_swap { ligands.rotate_left(2); }
        let stored = if open { StereoCoset::Undetermined } else { StereoCoset::Lit(coset) };
        let projected = if open { StereoCoset::Undetermined } else { StereoCoset::Lit(coset ^ u32::from(first_swap) ^ u32::from(second_swap)) };
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: base.atoms().iter().map(|atom| atom.attributes.clone()).collect(),
            bonds: base.bonds().iter().map(|bond| {
                let [first, second] = bond.atom_ids();
                (first, second, bond.attributes.clone())
            }).collect(),
            stereo_bonds: vec![(BondId(0), ligands, StereoBondForm::new(StereoKind::CisTrans, stored))],
            ..Default::default()
        });
        let mut editor = base.edit();
        editor.bond_mut(BondId(0)).attributes.constraints.set(
            BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(projected)));
        let expected = editor.build();
        let model = ChemistryModel::default();
        let resolver = Resolver::new(&model);
        prop_assert_eq!(resolver.stereo.project(&mut molecule), Ok(Solution::Determined(())));
        prop_assert_eq!(&molecule, &expected);
        prop_assert_eq!(resolver.stereo.project(&mut molecule), Ok(Solution::Determined(())));
        prop_assert_eq!(molecule, expected);
    }
}
