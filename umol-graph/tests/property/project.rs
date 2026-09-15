//! Aromatic projection preserves member-aligned contributions and atom electron fields.
//! Independent five- and six-member ring expectations check exact projected values,
//! idempotence, and transport under rotated/reversed system frames. The same expectations
//! check projection of successfully ingested SMILES under both valence sources.
//! Stereo projection is checked against independent tetrahedral permutation parity and
//! cis-trans side-swap parity, preserving virtual-ligand kinds and open configurations.
//! Composite projection of supported SMILES is compared with independent complete graph-IR
//! expectations, including aromatic contributions, both stereo assertions, and isotope elision.
//! Convey compares H omission with independent allowed-H sets for custom registries and bounded
//! counts tables, including isotope retention, field preservation, and repeated conversion.
//! Stage selection is compared with standalone composition for every flag combination on
//! ingested stereo/aromatic molecules; every flag combination preserves H counts.
//! Reaction convey checks selective H omission, compacted atom-map indices, and reaction equivalence
//! for alcohol oxidation with partial correspondence and optional creation/deletion.
//! Export fixtures compare source orderings by canonical GraphIR equality across valence and
//! isotope policies; generated cases cover fused/linked aromatic systems and mapped reaction stereo.
//! MOL export compares coordinate-derived stereo with independently specified SMILES, including
//! absent coordinates and the resolver's undetermined Either outcome. Initial marker spelling
//! may normalize.

use std::borrow::Cow;
use std::collections::BTreeMap;

use proptest::prelude::*;
use rstest::rstest;
use umol_chem::element::Element;
use umol_chem::spin::SpinMultiplicity;
use umol_graph::export::{export_reaction_smiles_with, export_smiles_with, Convey};
use umol_graph::ingest::{ingest_reaction_smiles_with, ingest_smiles, ingest_smiles_with};
use umol_graph::ops::model::{ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::valence::ValenceResolver;
use umol_graph::ops::resolve::{IsotopePolicy, ProjectFlags, ResolveConfig, Resolver};
use umol_graph::ops::valence::{AtomTypeRegistry, ResolveReport, ValenceEntry, ValenceTable};
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticValenceForm, AtomConstraintForm, AtomForm, AtomId,
    BondConstraintForm, BondForm, BondId, BooleanForm, Canonicalize, CanonicalizeContext,
    CisTransStereoForm, ElectronCountsForm, ElementForm, IsotopeMassForm, Molecule,
    MoleculeEntries, NumForm, StereoAtomForm, StereoBondForm, StereoCoset, StereoKind,
    StereoLigand, StereoLigandKind, TetrahedralStereoForm, UnpairedElectronsForm,
};
use umol_graph_ir::{atom_dsl, mol_dsl_concrete};
use umol_io::ctfile::parser::parse_mol_to_ir;
use umol_io::smiles::{ReactionSmiles, Smiles, SmilesIoConfig};
use umol_io::table_ir::{Atom, Molecule as TableMolecule};
use umol_utils::solution::Solution;

proptest! {
    #[test]
    fn test_smiles_convey_roundtrip(
        chain in 1usize..10, clockwise in any::<bool>(), trans in any::<bool>(),
        typing in any::<bool>(), explicit_h in any::<bool>(),
        ring in prop::sample::select(vec!["c1ccccc1", "c1cc[nH]c1", "c1ccc2ccccc2c1", "c1ccccc1-c1ccccc1"]),
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
        let boundary = Smiles::convey(&source, &model, &config, &io).unwrap();
        prop_assert_eq!(&source, &original);
        let text = boundary.render_with(&io).unwrap();
        prop_assert_eq!(export_smiles_with(&source, &io, &model, &config).unwrap(), text.as_str());
        prop_assert_eq!(&source, &original);
        let restored = ingest_smiles_with(&text, &io, &model, &config).unwrap();
        let context = CanonicalizeContext {
            para_stereo: false, automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        };
        prop_assert!(source.canonical_eq(&restored, &context));
        let repeated = Smiles::convey(&restored, &model, &config, &io).unwrap().render_with(&io).unwrap();
        prop_assert_eq!(text, repeated);
    }

    #[test]
    fn test_reaction_smiles_convey_roundtrip(
        mut mapped in prop::collection::vec(any::<bool>(), 2..10),
        label_start in 1u32..1000, label_step in 1u32..10,
        components in any::<bool>(), typing in any::<bool>(), saturated in any::<bool>(),
    ) {
        let count = mapped.len();
        let mut left = String::new();
        let mut right = String::new();
        let mut left_h = Vec::new();
        let mut right_h = Vec::new();
        for (index, &paired) in mapped.iter().enumerate() {
            let label = if paired { format!(":{}", label_start + index as u32 * label_step) } else { String::new() };
            let before = if index == 0 { 3 } else { 2 };
            let after = if index == count - 1 { 1 } else { before };
            left.push_str(&format!("[CH{before}{label}]"));
            right.push_str(&format!("[CH{after}{label}]"));
            left_h.push(if paired || !saturated { Some(before) } else { None });
            right_h.push(if paired || (!saturated && after != 1) { Some(after) } else { None });
        }
        let oxygen_label = label_start + count as u32 * label_step;
        left.push_str(&format!("[OH:{oxygen_label}]"));
        right.push_str(&format!("=[O:{oxygen_label}]"));
        mapped.push(true);
        left_h.push(Some(1));
        right_h.push(Some(0));
        let input = if components { format!("[OH2].{left}>>{right}.[NH3]") } else { format!("{left}>>{right}") };
        let mut model = ChemistryModel {
            valence: if typing { ValenceModel::default() } else { ValenceModel::smiles() },
            ..Default::default()
        };
        model.valence.tie_break = if saturated { ValenceTieBreak::MostSaturated } else { ValenceTieBreak::Strict };
        let config = ResolveConfig { isotope: IsotopePolicy::Natural, ..Default::default() };
        let io = SmilesIoConfig::opensmiles();
        let source = ingest_reaction_smiles_with(&input, &io, &model, &config).unwrap();
        let original = source.clone();
        let boundary = ReactionSmiles::convey(&source, &model, &config, &io).unwrap();
        prop_assert_eq!(&source, &original);
        let mut expected_mapping = BTreeMap::new();
        for (rank, (index, _)) in mapped.iter().enumerate().filter(|(_, paired)| **paired).enumerate() {
            expected_mapping.insert(rank as u32 + 1, (vec![index as u32 + u32::from(components)], vec![rank as u32]));
        }
        prop_assert_eq!(&boundary.as_table_ir().atom_mapping, &expected_mapping);
        let mut expected_right_h = mapped.iter().zip(&right_h).filter(|(paired, _)| **paired).map(|(_, h)| *h)
            .chain(mapped.iter().zip(&right_h).filter(|(paired, _)| !**paired).map(|(_, h)| *h)).collect::<Vec<_>>();
        if components {
            left_h.insert(0, if saturated || typing { None } else { Some(2) });
            expected_right_h.push(if saturated { None } else { Some(3) });
        }
        prop_assert_eq!(boundary.as_table_ir().reactants.atoms.iter().map(|atom| atom.implicit_hydrogens).collect::<Vec<_>>(), left_h);
        prop_assert_eq!(boundary.as_table_ir().products.atoms.iter().map(|atom| atom.implicit_hydrogens).collect::<Vec<_>>(), expected_right_h);
        let text = boundary.render().unwrap();
        prop_assert_eq!(export_reaction_smiles_with(&source, &io, &model, &config).unwrap(), text.as_str());
        prop_assert_eq!(&source, &original);
        let restored = ingest_reaction_smiles_with(&text, &io, &model, &config).unwrap();
        let context = CanonicalizeContext { para_stereo: false, automorphism_algorithm: AutomorphismAlgorithm::Nauty };
        prop_assert!(source.canonical_eq(&restored, &context));
        let repeated = ReactionSmiles::convey(&restored, &model, &config, &io).unwrap().render().unwrap();
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
            let mut molecule = original.clone();
            prop_assert_eq!(ValenceResolver::new(&model.valence).project(&mut molecule, policy), Ok(Solution::Determined(())));
            prop_assert_eq!(&molecule, &original);
        }
    }

    #[test]
    fn test_smiles_convey_atom_typing(
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
                let expected = TableMolecule {
                    atoms: vec![Atom {
                        isotope_mass: if natural { None } else { Some(13) },
                        implicit_hydrogens: if natural && selected == Some(stored_h) { None } else { Some(stored_h as u8) },
                        charge: Some(0), lone_pairs: Some(lp as u8), unpaired_electrons: Some(0),
                        multiplicity: Some(SpinMultiplicity::SINGLET), valence: Some(valence as u8),
                        aromatic: Some(is_aromatic), ..Atom::from_element(Element::C)
                    }],
                    ..TableMolecule::empty()
                };
                for valence_model in &models {
                    let mut model = ChemistryModel { valence: valence_model.clone(), ..Default::default() };
                    model.valence.tie_break = policy;
                    let config = ResolveConfig { isotope: IsotopePolicy::Natural, ..Default::default() };
                    let molecule = original.clone();
                    let conveyed = Smiles::convey(&molecule, &model, &config, &SmilesIoConfig::opensmiles()).unwrap();
                    prop_assert_eq!(conveyed.as_table_ir(), &expected);
                    prop_assert_eq!(Smiles::convey(&molecule, &model, &config, &SmilesIoConfig::opensmiles()).unwrap(), conveyed);
                    prop_assert_eq!(&molecule, &original);
                }
            }
        }
    }

    #[test]
    fn test_smiles_convey_counts(
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
        let mut model = ChemistryModel {
            valence: ValenceModel::counts(Cow::Owned(table)), ..Default::default()
        };
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
            let original = Molecule::from_entries(MoleculeEntries { atoms: vec![atom], ..Default::default() });
            let expected = TableMolecule {
                atoms: vec![Atom {
                    isotope_mass: if natural { None } else { Some(13) },
                    implicit_hydrogens: if natural && selected == Some(stored_h) { None } else { Some(stored_h as u8) },
                    charge: Some(0), lone_pairs: Some(stored_lp as u8), unpaired_electrons: Some(0),
                    multiplicity: Some(SpinMultiplicity::SINGLET), valence: Some(valence as u8),
                    aromatic: Some(is_aromatic), ..Atom::from_element(Element::C)
                }],
                ..TableMolecule::empty()
            };
            model.valence.tie_break = policy;
            let config = ResolveConfig { isotope: IsotopePolicy::Natural, ..Default::default() };
            let molecule = original.clone();
            let conveyed = Smiles::convey(&molecule, &model, &config, &SmilesIoConfig::opensmiles()).unwrap();
            prop_assert_eq!(conveyed.as_table_ir(), &expected);
            prop_assert_eq!(Smiles::convey(&molecule, &model, &config, &SmilesIoConfig::opensmiles()).unwrap(), conveyed);
            prop_assert_eq!(&molecule, &original);
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
            let retained = projected.atoms().iter().map(|atom| atom.implicit_hydrogens().clone()).collect::<Vec<_>>();
            prop_assert_eq!(&retained, &hydrogens);
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

    #[test]
    fn test_export_smiles_with_mol(
        state in 0u8..4, typing in any::<bool>(),
        reverse_atoms in any::<bool>(), reverse_bonds in any::<bool>(),
        offset in -20i16..20, scale in 1u16..10,
    ) {
        let points = if state == 1 || state == 2 {
            [[0., 1.], [0., 0.], [2., 0.], [2., if state == 1 { 1. } else { -1. }]]
        } else { [[0., 0.]; 4] };
        let mut points = points.map(|point| point.map(|value| value * f64::from(scale) + f64::from(offset)));
        if reverse_atoms { points.reverse(); }
        let mut input = String::from("stereo\n  umol          2D\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n");
        for (index, [x, y]) in points.into_iter().enumerate() {
            input.push_str(&format!("{x:10.4}{y:10.4}{:10.4} {:3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
                0., if index == 0 || index == 3 { "F" } else { "C" }));
        }
        let mut bonds = [(1, 2, 1), (2, 3, 2), (3, 4, 1)];
        if reverse_bonds { bonds.reverse(); }
        for (a, b, order) in bonds {
            let (a, b) = if reverse_atoms { (5 - a, 5 - b) } else { (a, b) };
            let stereo = if order == 2 && state == 3 { 3 } else { 0 };
            input.push_str(&format!("{a:3}{b:3}{order:3}{stereo:3}  0  0  0\n"));
        }
        input.push_str("M  END\n");
        let model = ChemistryModel {
            valence: ValenceModel { tie_break: ValenceTieBreak::MostSaturated,
                ..if typing { ValenceModel::default() } else { ValenceModel::smiles() } },
            ..Default::default()
        };
        let config = ResolveConfig { isotope: IsotopePolicy::Natural, ..Default::default() };
        let io = SmilesIoConfig::opensmiles();
        let mut source = parse_mol_to_ir(&input).unwrap();
        let raised = source.clone();
        let result = Resolver::with_config(&model, config).resolve(&mut source).unwrap();
        if state == 3 {
            prop_assert_eq!(result, Solution::Underdetermined(ResolveReport::default()));
            prop_assert_eq!(source, raised);
        } else {
            prop_assert!(matches!(result, Solution::Determined(_)), "{result:?}");
            let original = source.clone();
            let text = export_smiles_with(&source, &io, &model, &config).unwrap();
            let expected = match state { 1 => "F/C=C\\F", 2 => "F/C=C/F", _ => "FC=CF" };
            let expected = ingest_smiles_with(expected, &io, &model, &config).unwrap();
            let restored = ingest_smiles_with(&text, &io, &model, &config).unwrap();
            let context = CanonicalizeContext { para_stereo: false, automorphism_algorithm: AutomorphismAlgorithm::Nauty };
            prop_assert!(source.canonical_eq(&expected, &context));
            prop_assert!(source.canonical_eq(&restored, &context));
            prop_assert_eq!(export_smiles_with(&source, &io, &model, &config).unwrap(), text);
            prop_assert_eq!(source, original);
        }
    }

    #[test]
    fn test_export_reaction_smiles_with_stereo(
        labels in Just(vec![0u32, 1, 2, 3]).prop_shuffle(), start in 1u32..1000,
        tetrahedral in any::<bool>(), flip in any::<bool>(), typing in any::<bool>(),
        saturated in any::<bool>(), natural in any::<bool>(),
    ) {
        let [a, b, c, d] = [0, 1, 2, 3].map(|index| labels[index] + start);
        let input = if tetrahedral {
            let winding = if flip { "@" } else { "@@" };
            format!("[F:{a}][C@H:{b}]([Cl:{c}])[Br:{d}]>>[Br:{d}][C{winding}H:{b}]([Cl:{c}])[F:{a}]")
        } else {
            let marker = if flip { "\\" } else { "/" };
            format!("[F:{a}]/[CH:{b}]=[CH:{c}]/[Cl:{d}]>>[Cl:{d}]/[CH:{c}]=[CH:{b}]{marker}[F:{a}]")
        };
        let model = ChemistryModel {
            valence: ValenceModel {
                tie_break: if saturated { ValenceTieBreak::MostSaturated } else { ValenceTieBreak::Strict },
                ..if typing { ValenceModel::default() } else { ValenceModel::smiles() }
            },
            ..Default::default()
        };
        let config = ResolveConfig { isotope: if natural { IsotopePolicy::Natural } else { IsotopePolicy::Strict }, ..Default::default() };
        let io = SmilesIoConfig::opensmiles();
        let source = ingest_reaction_smiles_with(&input, &io, &model, &config).unwrap();
        let original = source.clone();
        let text = export_reaction_smiles_with(&source, &io, &model, &config).unwrap();
        let restored = ingest_reaction_smiles_with(&text, &io, &model, &config).unwrap();
        let context = CanonicalizeContext { para_stereo: false, automorphism_algorithm: AutomorphismAlgorithm::Nauty };
        prop_assert!(source.canonical_eq(&restored, &context));
        prop_assert_eq!(export_reaction_smiles_with(&source, &io, &model, &config).unwrap(), text.as_str());
        prop_assert_eq!(export_reaction_smiles_with(&restored, &io, &model, &config).unwrap(), text);
        prop_assert_eq!(source, original);
    }
}

#[rstest]
#[case::chain("[CH3][CH2][OH]", "[OH][CH2][CH3]")]
#[case::components("[NH4+].[Cl-]", "[Cl-].[NH4+]")]
#[case::isotope("[13CH3][CH2][OH]", "[OH][CH2][13CH3]")]
#[case::radical("[CH3].[OH]", "[OH].[CH3]")]
#[case::implicit_h("[F][C@H]([Cl])[Br]", "[Br][C@@H]([Cl])[F]")]
#[case::actual_h("[H][C@]([F])([Cl])[Br]", "[Br][C@@]([F])([Cl])[H]")]
#[case::four_ligands("[F][C@]([Cl])([Br])[I]", "[I][C@@]([Cl])([Br])[F]")]
#[case::lone_pair("[CH3][S@](=[O])[CH2][CH3]", "[CH3][CH2][S@@](=[O])[CH3]")]
#[case::alkene("[F]/[CH]=[CH]/[Cl]", "[Cl]/[CH]=[CH]/[F]")]
#[case::four_substituents("[F]/[C]([Cl])=[C]([Br])/[I]", "[I]/[C]([Br])=[C]([Cl])/[F]")]
#[case::partial("[F]/[CH]=[CH][Cl]", "[F][CH]=[CH][Cl]")]
#[case::conjugated("[F]/[CH]=[CH]/[CH]=[CH]/[Cl]", "[Cl]\\[CH]=[CH]\\[CH]=[CH]\\[F]")]
#[case::ring_labels(
    "[CH2]%12[CH2][CH2][CH2][CH2][CH2]%12",
    "[CH2]1[CH2][CH2][CH2][CH2][CH2]1"
)]
#[case::benzene("[cH]1[cH][cH][cH][cH][cH]1", "[cH]1:[cH]:[cH]:[cH]:[cH]:[cH]:1")]
#[case::fused(
    "[cH]1[cH][cH][c]2[cH][cH][cH][cH][c]2[cH]1",
    "[c]12[cH][cH][cH][cH][c]1[cH][cH][cH][cH]2"
)]
#[case::linked(
    "[cH]1[cH][cH][cH][cH][c]1-[c]2[cH][cH][cH][cH][cH]2",
    "[cH]1[cH][cH][c](-[c]2[cH][cH][cH][cH][cH]2)[cH][cH]1"
)]
#[case::heteroaromatic("[nH]1[cH][cH][cH][cH]1", "[cH]1[cH][nH][cH][cH]1")]
fn test_export_smiles_with_order(
    #[case] first: &str,
    #[case] second: &str,
    #[values(false, true)] typing: bool,
    #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] isotope: IsotopePolicy,
) {
    let model = ChemistryModel {
        valence: ValenceModel {
            tie_break: policy,
            ..if typing {
                ValenceModel::default()
            } else {
                ValenceModel::smiles()
            }
        },
        ..Default::default()
    };
    let config = ResolveConfig {
        isotope,
        ..Default::default()
    };
    let io = SmilesIoConfig::opensmiles();
    let context = CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    let first = ingest_smiles_with(first, &io, &model, &config).unwrap();
    let second = ingest_smiles_with(second, &io, &model, &config).unwrap();
    assert!(first.canonical_eq(&second, &context));
    for source in [first, second] {
        let original = source.clone();
        let text = export_smiles_with(&source, &io, &model, &config).unwrap();
        let restored = ingest_smiles_with(&text, &io, &model, &config).unwrap();
        assert!(source.canonical_eq(&restored, &context), "{text}");
        assert_eq!(
            export_smiles_with(&source, &io, &model, &config),
            Ok(text.clone())
        );
        assert_eq!(
            export_smiles_with(&restored, &io, &model, &config),
            Ok(text)
        );
        assert_eq!(source, original);
    }
}

#[rstest]
#[case::methane("[CH4]", "[H][CH3]")]
#[case::tetrahedral("[F][C@H]([Cl])[Br]", "[F][C@]([H])([Cl])[Br]")]
fn test_export_smiles_with_hydrogens(
    #[case] implicit: &str,
    #[case] actual: &str,
    #[values(false, true)] typing: bool,
    #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
) {
    let model = ChemistryModel {
        valence: ValenceModel {
            tie_break: policy,
            ..if typing {
                ValenceModel::default()
            } else {
                ValenceModel::smiles()
            }
        },
        ..Default::default()
    };
    let config = ResolveConfig {
        isotope: IsotopePolicy::Natural,
        ..Default::default()
    };
    let io = SmilesIoConfig::opensmiles();
    let implicit = ingest_smiles_with(implicit, &io, &model, &config).unwrap();
    let actual = ingest_smiles_with(actual, &io, &model, &config).unwrap();
    let context = CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    assert!(!implicit.canonical_eq(&actual, &context));
    let implicit_text = export_smiles_with(&implicit, &io, &model, &config).unwrap();
    let actual_text = export_smiles_with(&actual, &io, &model, &config).unwrap();
    let implicit_restored = ingest_smiles_with(&implicit_text, &io, &model, &config).unwrap();
    let actual_restored = ingest_smiles_with(&actual_text, &io, &model, &config).unwrap();
    assert!(implicit.canonical_eq(&implicit_restored, &context));
    assert!(actual.canonical_eq(&actual_restored, &context));
    assert!(!implicit_restored.canonical_eq(&actual_restored, &context));
}
