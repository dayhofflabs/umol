//! Aromatic projection preserves member-aligned contributions and atom electron fields.
//! Independent five- and six-member ring expectations check exact projected values,
//! idempotence, and transport under rotated/reversed system frames. The same expectations
//! check projection of successfully ingested SMILES under both valence sources.
//! Stereo projection is checked against independent tetrahedral permutation parity and
//! cis-trans side-swap parity, preserving virtual-ligand kinds and open configurations.
//! Composite projection of supported SMILES is compared with independent complete graph-IR
//! expectations, including aromatic contributions, both stereo assertions, and isotope elision.
//! Valence projection compares H elision with independent allowed-H sets for custom registries
//! and bounded counts tables, including isotope retention, field preservation, and idempotence.
//! Stage selection is compared with standalone composition for every flag combination on
//! ingested stereo/aromatic molecules; omitting valence must preserve every H count.

use std::borrow::Cow;

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph::export::Convey;
use umol_graph::ingest::{ingest_smiles, ingest_smiles_with};
use umol_graph::ops::model::{ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::valence::ValenceResolver;
use umol_graph::ops::resolve::{IsotopePolicy, ProjectFlags, ResolveConfig, Resolver};
use umol_graph::ops::valence::{AtomTypeRegistry, ValenceEntry, ValenceTable};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticValenceForm, AtomConstraintForm, AtomForm, AtomId,
    BondConstraintForm, BondForm, BondId, BooleanForm, Canonicalize, CanonicalizeContext,
    CisTransStereoForm, ElectronCountsForm, ElementForm, IsotopeMassForm, Molecule,
    MoleculeEntries, NumForm, StereoAtomForm, StereoBondForm, StereoCoset, StereoKind,
    StereoLigand, StereoLigandKind, TetrahedralStereoForm, UnpairedElectronsForm,
};
use umol_graph_ir::{atom_dsl, mol_dsl_concrete};
use umol_io::smiles::{Smiles, SmilesIoConfig};
use umol_utils::solution::Solution;

proptest! {
    #[test]
    fn test_smiles_convey_roundtrip(
        chain in 1usize..10, clockwise in any::<bool>(), trans in any::<bool>(),
        typing in any::<bool>(), explicit_h in any::<bool>(),
        ring in prop::sample::select(vec!["c1ccccc1", "c1cc[nH]c1"]),
    ) {
        let chirality = if clockwise { "@@" } else { "@" };
        let direction = if trans { "/" } else { "\\" };
        let hydrogen = if explicit_h { "[H]" } else { "" };
        let input = format!("{hydrogen}{}[C{chirality}H](F)/C=C{direction}{ring}", "C".repeat(chain));
        let mut model = ChemistryModel {
            valence: if typing { ValenceModel::default() } else { ValenceModel::smiles() },
            ..Default::default()
        };
        model.valence.tie_break = ValenceTieBreak::MostSaturated;
        let config = ResolveConfig { isotope: IsotopePolicy::Natural, ..Default::default() };
        let io = SmilesIoConfig::opensmiles();
        let source = ingest_smiles_with(&input, &io, &model, &config).unwrap();
        let original = source.clone();
        let resolver = Resolver::with_config(&model, config);
        let boundary = Smiles::convey(&source, &resolver, &io).unwrap();
        prop_assert_eq!(&source, &original);
        let text = boundary.render_with(&io).unwrap();
        let restored = ingest_smiles_with(&text, &io, &model, &config).unwrap();
        let context = CanonicalizeContext {
            para_stereo: false, automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        };
        prop_assert!(source.canonical_eq(&restored, &context));
        let repeated = Smiles::convey(&restored, &resolver, &io).unwrap().render_with(&io).unwrap();
        prop_assert_eq!(text, repeated);
    }

    #[test]
    fn test_valence_resolver_project_input(
        atoms in 2usize..25, mass in any::<bool>(), typing in any::<bool>(),
    ) {
        let input = format!("{}{}[CH3]", if mass { "[13CH3]" } else { "[CH3]" }, "[CH2]".repeat(atoms - 2));
        let model = ChemistryModel {
            valence: if typing { ValenceModel::default() } else { ValenceModel::smiles() },
            ..Default::default()
        };
        let original = ingest_smiles_with(&input, &SmilesIoConfig::opensmiles(), &model, &ResolveConfig::default()).unwrap();
        for policy in [ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated] {
            let mut expected = original.clone();
            if policy == ValenceTieBreak::MostSaturated {
                for index in usize::from(mass)..atoms {
                    expected.atom_mut(AtomId(index as u32)).attributes.implicit_hydrogens = NumForm::Undetermined;
                }
            }
            let mut molecule = original.clone();
            prop_assert_eq!(ValenceResolver::new(&model.valence).project(&mut molecule, policy), Ok(Solution::Determined(())));
            prop_assert_eq!(molecule, expected);
        }
    }

    #[test]
    fn test_valence_resolver_project_atom_typing(
        entries in prop::collection::vec((0i64..5, 0i64..5, 0i64..3), 0..24),
        valence in 0i64..5, aromatic in 0i64..3,
        stored_h in 0i64..5, stored_lp in 0i64..5,
        copies in 1usize..5, natural in any::<bool>(), is_aromatic in any::<bool>(),
    ) {
        let aromatic = if is_aromatic { aromatic } else { 0 };
        let rows = entries.iter().filter(|&&(_, _, a)| is_aromatic || a == 0).map(|&(h, v, a)| {
            let mut row = atom_dsl!("C#c0#u0#s");
            row.implicit_hydrogens = NumForm::Lit(h);
            row.constraints.set(AtomConstraintForm::valence(v));
            row.constraints.set(AtomConstraintForm::aromatic_valence(if is_aromatic { AromaticValenceForm::aromatic(a) } else { AromaticValenceForm::NotAromatic }));
            row
        }).collect::<Vec<_>>();
        let models = [
            ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms(rows.clone()))),
            ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms(
                (0..copies).flat_map(|_| rows.iter().rev().cloned()),
            ))),
        ];
        let mut allowed = entries.iter().filter(|(_, v, a)| *v == valence && *a == aromatic)
            .map(|(h, _, _)| *h).collect::<Vec<_>>();
        allowed.sort_unstable();
        allowed.dedup();
        for policy in [ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated] {
            let selected = match policy {
                ValenceTieBreak::Strict => match allowed.as_slice() { [h] => Some(*h), _ => None },
                ValenceTieBreak::MostSaturated => allowed.last().copied(),
            };
            for lp in [0, stored_lp] {
                let mut atom = atom_dsl!("C#c0#u0#s");
                atom.isotope_mass = if natural { IsotopeMassForm::Natural } else { IsotopeMassForm::Lit(13) };
                atom.implicit_hydrogens = NumForm::Lit(stored_h);
                atom.lone_pairs = NumForm::Lit(lp);
                atom.constraints.set(AtomConstraintForm::valence(valence));
                atom.constraints.set(AtomConstraintForm::aromatic_valence(if is_aromatic { AromaticValenceForm::aromatic(aromatic) } else { AromaticValenceForm::NotAromatic }));
                let original = Molecule::from_entries(MoleculeEntries { atoms: vec![atom.clone()], ..Default::default() });
                if natural && (!is_aromatic || stored_h == 0) && selected == Some(stored_h) {
                    atom.implicit_hydrogens = NumForm::Undetermined;
                }
                let expected = Molecule::from_entries(MoleculeEntries { atoms: vec![atom], ..Default::default() });
                for model in &models {
                    let resolver = ValenceResolver::new(model);
                    let mut molecule = original.clone();
                    prop_assert_eq!(resolver.project(&mut molecule, policy), Ok(Solution::Determined(())));
                    prop_assert_eq!(&molecule, &expected);
                    prop_assert_eq!(resolver.project(&mut molecule, policy), Ok(Solution::Determined(())));
                    prop_assert_eq!(&molecule, &expected);
                }
            }
        }
    }

    #[test]
    fn test_valence_resolver_project_counts(
        targets in prop::collection::vec(0u8..7, 0..5),
        valence in 0i64..5, aromatic in 0i64..3,
        stored_h in 0i64..5, stored_lp in 0i64..5,
        natural in any::<bool>(), is_aromatic in any::<bool>(),
    ) {
        let aromatic = if is_aromatic { aromatic } else { 0 };
        let mut table = ValenceTable::empty();
        table.insert(Element::C, ValenceEntry {
            target_covalences: targets.clone(),
            aromatic_valences: vec![1], fallback_aromatic_valences: vec![0],
        });
        let model = ValenceModel::counts(Cow::Owned(table));
        let target = targets.iter().map(|&v| i64::from(v)).filter(|&v| v >= valence)
            .min().unwrap_or(valence);
        // Enumerate H/LP pairs satisfying electron balance and the first-target bound.
        let mut allowed = Vec::new();
        for h in 0..=4 {
            for lp in 0..=4 {
                if valence + aromatic + h + 2 * lp == 4
                    && valence + h + i64::from(aromatic == 1) <= target
                {
                    allowed.push(h);
                }
            }
        }
        for policy in [ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated] {
            let selected = match policy {
                ValenceTieBreak::Strict => match allowed.as_slice() { [h] => Some(*h), _ => None },
                ValenceTieBreak::MostSaturated => allowed.last().copied(),
            };
            let mut atom = atom_dsl!("C#c0#u0#s");
            atom.isotope_mass = if natural { IsotopeMassForm::Natural } else { IsotopeMassForm::Lit(13) };
            atom.implicit_hydrogens = NumForm::Lit(stored_h);
            atom.lone_pairs = NumForm::Lit(stored_lp);
            atom.constraints.set(AtomConstraintForm::valence(valence));
            atom.constraints.set(AtomConstraintForm::aromatic_valence(if is_aromatic { AromaticValenceForm::aromatic(aromatic) } else { AromaticValenceForm::NotAromatic }));
            let mut molecule = Molecule::from_entries(MoleculeEntries { atoms: vec![atom.clone()], ..Default::default() });
            if natural && (!is_aromatic || stored_h == 0) && selected == Some(stored_h) {
                atom.implicit_hydrogens = NumForm::Undetermined;
            }
            let expected = Molecule::from_entries(MoleculeEntries { atoms: vec![atom], ..Default::default() });
            let resolver = ValenceResolver::new(&model);
            prop_assert_eq!(resolver.project(&mut molecule, policy), Ok(Solution::Determined(())));
            prop_assert_eq!(&molecule, &expected);
            prop_assert_eq!(resolver.project(&mut molecule, policy), Ok(Solution::Determined(())));
            prop_assert_eq!(&molecule, &expected);
        }
    }
}

proptest! {
    #[test]
    fn test_resolver_project_composition(
        chain in 1usize..12, clockwise in any::<bool>(), typing in any::<bool>(),
        saturated in any::<bool>(), natural in any::<bool>(),
    ) {
        let chirality = if clockwise { "@@" } else { "@" };
        let input = format!("[13CH3][C{chirality}H](F)/C=C/c1ccccc1.C[S{chirality}](=O){}", "C".repeat(chain));
        let mut model = ChemistryModel {
            valence: if typing { ValenceModel::default() } else { ValenceModel::smiles() },
            ..Default::default()
        };
        model.valence.tie_break = if saturated { ValenceTieBreak::MostSaturated } else { ValenceTieBreak::Strict };
        let config = ResolveConfig {
            isotope: if natural { IsotopePolicy::Natural } else { IsotopePolicy::Strict },
            ..Default::default()
        };
        let original = ingest_smiles(&input).unwrap();
        let hydrogens = original.atoms().iter().map(|atom| atom.implicit_hydrogens().clone()).collect::<Vec<_>>();
        let resolver = Resolver::with_config(&model, config);
        for bits in 0..=ProjectFlags::all().bits() {
            let flags = ProjectFlags::from_bits_retain(bits);
            let mut expected = original.clone();
            if flags.contains(ProjectFlags::STEREO) {
                prop_assert_eq!(resolver.stereo.project(&mut expected), Ok(Solution::Determined(())));
            }
            if flags.contains(ProjectFlags::AROMATICITY) {
                prop_assert_eq!(resolver.aromaticity.project(&mut expected), Ok(Solution::Determined(())));
            }
            if flags.contains(ProjectFlags::VALENCE) {
                prop_assert_eq!(resolver.valence.project(&mut expected, model.valence.tie_break), Ok(Solution::Determined(())));
            }
            if flags.contains(ProjectFlags::ISOTOPE) {
                prop_assert_eq!(resolver.isotope.project(&mut expected), Ok(Solution::Determined(())));
            }
            let mut projected = original.clone();
            prop_assert_eq!(resolver.project(&mut projected, flags), Ok(Solution::Determined(())));
            prop_assert_eq!(&projected, &expected);
            if !flags.contains(ProjectFlags::VALENCE) {
                let retained = projected.atoms().iter().map(|atom| atom.implicit_hydrogens().clone()).collect::<Vec<_>>();
                prop_assert_eq!(&retained, &hydrogens);
            }
        }
    }

    #[test]
    fn test_resolver_project_input(
        mass in prop::sample::select(vec![None, Some(13u32), Some(14u32)]),
        clockwise in any::<bool>(), trans in any::<bool>(),
        natural in any::<bool>(), typing in any::<bool>(),
        saturated in any::<bool>(),
    ) {
        let isotope = mass.map(|mass| mass.to_string()).unwrap_or_default();
        let chirality = if clockwise { "@@" } else { "@" };
        let direction = if trans { "/" } else { "\\" };
        let input = format!("[{isotope}CH3][C{chirality}H](F)/[CH]=[CH]{direction}[c]1[cH][cH][cH][cH][cH]1");
        let mut model = ChemistryModel {
            valence: if typing { ValenceModel::default() } else { ValenceModel::smiles() },
            ..Default::default()
        };
        model.valence.tie_break = if saturated { ValenceTieBreak::MostSaturated } else { ValenceTieBreak::Strict };
        let config = ResolveConfig {
            isotope: if natural { IsotopePolicy::Natural } else { IsotopePolicy::Strict },
            ..Default::default()
        };
        let mut molecule = ingest_smiles_with(&input, &SmilesIoConfig::opensmiles(), &model, &config).unwrap();
        let base = mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h1" "F#n3" "C#h1" "C#h1" "C#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1" "C#h1#a1"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"] [3 4 "2"] [4 5 "1"]
                [5 10 "1#a"] [5 6 "1#a"] [6 7 "1#a"] [7 8 "1#a"] [8 9 "1#a"] [9 10 "1#a"]]}"#);
        let mut editor = base.edit();
        for index in 2..6 {
            editor.atom_mut(AtomId(index)).attributes.implicit_hydrogens = NumForm::Undetermined;
        }
        if saturated && mass.is_none() {
            editor.atom_mut(AtomId(0)).attributes.implicit_hydrogens = NumForm::Undetermined;
        }
        if natural {
            for atom in base.atoms().ids() {
                editor.atom_mut(atom).attributes.isotope_mass = IsotopeMassForm::Undetermined;
            }
        }
        if let Some(mass) = mass {
            editor.atom_mut(AtomId(0)).attributes.isotope_mass = IsotopeMassForm::Lit(mass);
        }
        editor.atom_mut(AtomId(1)).attributes.constraints.set(AtomConstraintForm::TetrahedralStereo(
            TetrahedralStereoForm::Stereo(StereoCoset::Lit(u32::from(clockwise)))));
        editor.bond_mut(BondId(3)).attributes.constraints.set(BondConstraintForm::CisTransStereo(
            CisTransStereoForm::Stereo(StereoCoset::Lit(u32::from(trans)))));
        let expected = editor.build();
        prop_assert_eq!(Resolver::with_config(&model, config).project(&mut molecule, ProjectFlags::all()), Ok(Solution::Determined(())));
        prop_assert_eq!(molecule, expected);
    }

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
