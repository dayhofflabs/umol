//! Property-based coverage for SMILES parsing and rendering.

use std::array;
use std::collections::{BTreeMap, BTreeSet};
use std::panic::catch_unwind;

use proptest::collection::vec;
use proptest::prelude::*;
use proptest::sample::select;
use proptest::test_runner::{Config, FileFailurePersistence};
use rstest::rstest;
use umol_chem::element::Element;
use umol_chem::spin::SpinMultiplicity;
use umol_geometric_core::Point3D;
use umol_graph_ir::ir::{
    BondId, CisTransStereoForm, ElementForm, Molecule, StereoCoset, TryIntoIr,
};
use umol_io::smiles::config::SmilesIoConfig;
use umol_io::smiles::{
    parse_extended_smiles_bytes, ParseError, ReactionSmiles, Smiles, SmilesRenderError,
};
use umol_io::table_ir::{
    Atom, Bond, BondConfiguration, BondOrder, BondRelation, ExtendedMolecule,
    Molecule as TableMolecule, Span, StereoAtom, StereoBond, StereoLigand, Winding,
};

// Generate ASCII strings from a token-friendly alphabet to bias towards SMILES-like inputs.
// This is intentionally permissive; the property is "no panics".
fn smilesish() -> impl Strategy<Value = Vec<u8>> {
    // Common SMILES characters: letters, digits, bonds, ring, parens, brackets, slash/backslash, percent, dot
    const ALPHABET: &[u8] =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-=#$:/\\().%[]*";
    vec(select(ALPHABET.to_vec()), 0..256)
}

fn wildcard_chain() -> impl Strategy<Value = Vec<u8>> {
    (1_usize..=128).prop_flat_map(|atom_count| {
        (0..atom_count, vec(any::<bool>(), atom_count)).prop_map(
            |(wildcard_position, mut wildcards)| {
                wildcards[wildcard_position] = true;
                wildcards
                    .into_iter()
                    .map(|wildcard| if wildcard { b'*' } else { b'C' })
                    .collect()
            },
        )
    })
}

// TODO: Add more fine-grained property tests
// - Chains of arbitrary length are allowed (single atom, multiple atoms, different bond orders, brackets)
// - Branches of arbitrary length are allowed
// - Groups of arbitrary length are allowed
// - Rings of arbitrary length are allowed
// - Arbitrary number of components is allowed
// - Arbitrary number of branches is allowed
// - Arbitrary number of groups is allowed
// - Arbitrary number of rings is allowed
// - Arbitrary nesting of branches is allowed
// - Arbitrary nesting of rings is allowed
// - Arbitrary nesting of branches and rings is allowed
// - Arbitrary number of chiral centers is allowed
// - Arbitrary number of stereogenic bonds is allowed

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::WithSource("proptest-regressions"))),
        cases: 10000,
        ..Config::default()}
    )]

    // Crash-only: parser should never panic on arbitrary ASCII up to length 256
    #[test]
    fn never_panics_on_ascii(input in smilesish()) {
        let config = SmilesIoConfig::opensmiles();
        catch_unwind(|| {
            let _ = Smiles::parse_bytes_with(&input, &config);
        }).expect("SMILES parsing panicked");
    }

    // Error spans must point within the input bounds
    #[test]
    fn error_positions_within_bounds(input in smilesish()) {
        let res = Smiles::parse_bytes(&input);
        if let Err(err) = res {
            let len = input.len();
            let ok = match err {
                ParseError::DanglingBondDirection { .. }
                | ParseError::CisTransConflict { .. }
                | ParseError::UnsupportedStereoBond { .. }
                | ParseError::ConflictingBondConfiguration { .. }
                | ParseError::MissingPosition { .. }
                | ParseError::LeadingWhitespace => true,
                | ParseError::InvalidElement { pos } => pos < len,
                | ParseError::InvalidToken { pos } => pos < len,
                | ParseError::UnbalancedOpenParen { pos } => pos < len,
                | ParseError::UnbalancedCloseParen { pos } => pos < len,
                | ParseError::EmptyBranch { pos } => pos < len,
                | ParseError::EmptyGroup { pos } => pos < len,
                | ParseError::NonfinalGroup { pos } => pos < len,
                | ParseError::LeadingBond { pos } => pos < len,
                | ParseError::TrailingBond { pos } => pos < len,
                | ParseError::ConsecutiveBonds { pos } => pos < len,
                | ParseError::LeadingRing { pos } => pos < len,
                | ParseError::UnbalancedRingIndex { open_pos } => open_pos < len,
                | ParseError::InvalidRingIndex { pos } => pos < len,
                | ParseError::MismatchedRingBondOrders { pos, open_pos } => pos < len && open_pos < len,
                | ParseError::MismatchedRingBondDirections { pos, open_pos } => pos < len && open_pos < len,
                | ParseError::MismatchedRingBondDonations { pos, open_pos } => pos < len && open_pos < len,
                | ParseError::LeadingDot { pos } => pos < len,
                | ParseError::TrailingDot { pos } => pos < len,
                | ParseError::ConsecutiveDots { pos } => pos < len,
                | ParseError::DotBeforeRing { pos } => pos < len,
                | ParseError::EmptyBracket { pos } => pos < len,
                | ParseError::UnbalancedOpenBracket { pos } => pos < len,
                | ParseError::UnbalancedCloseBracket { pos } => pos < len,
                | ParseError::StrayBracketField { pos } => pos < len,
                | ParseError::DuplicateBracketField { pos } => pos < len,
                | ParseError::MissingClassIndex { pos } => pos < len,
                | ParseError::MissingChiralityIndex { pos } => pos < len,
                | ParseError::ChiralityOutOfRange { pos } => pos < len,
                | ParseError::BracketHwithHcount { pos } => pos < len,
                | ParseError::InvalidBracket { pos } => pos < len,
                | ParseError::InvalidCxTag { pos } => pos < len,
                | ParseError::AtomIndexOutOfBounds { .. } => true,
                | ParseError::BondIndexOutOfBounds { .. } => true,
                | ParseError::MismatchedAtomBondIndices { .. } => true,
                | ParseError::SgroupIndexOutOfBounds { .. } => true,
                | ParseError::MissingReactionArrow { pos } => pos <= len,
            };
            prop_assert!(ok, "error positions out of bounds: {:?}, len={}", err, len);
        }
    }

    // Bonds in successful parses must reference valid atom indices
    // Note: self-loop bonds (e.g., C11) are syntactically valid and checked during topology validation
    #[test]
    fn bonds_well_formed_on_success(input in smilesish()) {
        if let Ok(smiles) = Smiles::parse_bytes(&input) {
            let mol = smiles.as_table_ir();
            let n = mol.atoms.len() as u32;
            for b in &mol.bonds {
                let sa = b.start_atom();
                let ea = b.end_atom();
                // Bonds use 0-based atom indices
                prop_assert!(sa < n && ea < n, "bond endpoints out of bounds: {}-{} / n={}", sa, ea, n);
            }
        }
    }
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::WithSource(
            "proptest-regressions",
        ))),
        cases: 1_000,
        ..Config::default()
    })]

    #[test]
    fn test_smiles_wildcard_chain(input in wildcard_chain()) {
        let molecule = Smiles::parse_bytes(&input)
            .expect("generated wildcard chain should parse")
            .into_table_ir();

        prop_assert_eq!(molecule.atoms.len(), input.len());
        let expected_wildcards: Vec<_> = input
            .iter()
            .enumerate()
            .filter_map(|(position, symbol)| (*symbol == b'*').then_some(position))
            .collect();
        let actual_wildcards: Vec<_> = molecule
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(position, atom)| atom.element.is_none().then_some(position))
            .collect();
        prop_assert_eq!(actual_wildcards, expected_wildcards);

        let expected_atoms: Vec<_> = input
            .iter()
            .enumerate()
            .map(|(position, symbol)| {
                (
                    (*symbol == b'C').then_some(Element::C),
                    Some(Span::bytes(position as u32, position as u32 + 1)),
                )
            })
            .collect();
        let actual_atoms: Vec<_> = molecule
            .atoms
            .iter()
            .map(|atom| (atom.element, atom.span))
            .collect();
        prop_assert_eq!(actual_atoms, expected_atoms);

        let expected_bonds: Vec<_> = (1..input.len())
            .map(|position| {
                (
                    position as u32 - 1,
                    position as u32,
                    Some(Span::bytes(position as u32, position as u32 + 1)),
                )
            })
            .collect();
        let actual_bonds: Vec<_> = molecule
            .bonds
            .iter()
            .map(|bond| (bond.start_atom(), bond.end_atom(), bond.span))
            .collect();
        prop_assert_eq!(actual_bonds, expected_bonds);
    }

    #[test]
    fn test_smiles_wildcard_chain_differential(input in wildcard_chain()) {
        let basic = Smiles::parse_bytes(&input)
            .expect("generated wildcard chain should parse")
            .into_table_ir();
        let extended = parse_extended_smiles_bytes(&input)
            .expect("generated wildcard chain should parse in extended mode");
        prop_assert_eq!(ExtendedMolecule::from(basic), extended);
    }

    #[test]
    fn test_smiles_wildcard_chain_raise(input in wildcard_chain()) {
        let molecule = Smiles::parse_bytes(&input)
            .expect("generated wildcard chain should parse")
            .into_table_ir();
        let molecule = (&molecule)
            .try_into_ir(&())
            .expect("generated wildcard chain should raise");
        let expected: Vec<_> = input
            .iter()
            .map(|symbol| {
                if *symbol == b'C' {
                    ElementForm::Lit(Element::C)
                } else {
                    ElementForm::Undetermined
                }
            })
            .collect();
        let actual: Vec<_> = molecule
            .atoms()
            .iter()
            .map(|atom| atom.element().clone())
            .collect();
        prop_assert_eq!(actual, expected);
    }

    // Each alkene is definite exactly when both adjacent single bonds are marked.
    // Same glyphs along this linear traversal put the reference ligands opposite.
    #[test]
    fn test_smiles_parse_directional_chain(markers in vec(prop::option::of(any::<bool>()), 2..10)) {
        let mut input = String::from("F");
        for (index, marker) in markers.iter().enumerate() {
            match marker {
                Some(true) => input.push('/'),
                Some(false) => input.push('\\'),
                None => {},
            }
            input.push_str(if index + 1 == markers.len() { "F" } else { "C=C" });
        }
        let expected: Vec<_> = markers.windows(2).map(|pair| {
            match (pair[0], pair[1]) {
                (Some(left), Some(right)) => Some(CisTransStereoForm::stereo(
                    StereoCoset::Lit(u32::from(left == right)),
                )),
                _ => None,
            }
        }).collect();
        let table = Smiles::parse(&input).unwrap().into_table_ir();
        let molecule: Molecule = (&table).try_into_ir(&()).unwrap();
        let actual: Vec<_> = (0..markers.len() - 1).map(|index| {
            molecule.bond(BondId::from(2 * index + 1))
                .attributes.constraints.cis_trans_stereo().cloned()
        }).collect();
        prop_assert_eq!(actual, expected);
    }

    // Ring labels and atom-map numbers do not change the ligand encounter frame.
    #[test]
    fn test_smiles_parse_ring_frame(label in 10_u32..100, class in 0_u32..10000, clockwise in any::<bool>()) {
        let marker = if clockwise { "@@" } else { "@" };
        let input = format!("C[C{marker}H:{class}]%{label}CCCCO%{label}");
        let table = Smiles::parse(&input).unwrap().into_table_ir();
        prop_assert_eq!(table.stereo_atoms, vec![StereoAtom {
            atom: 1,
            ligands: vec![StereoLigand::Atom(0), StereoLigand::ImplicitHydrogen,
                          StereoLigand::Atom(6), StereoLigand::Atom(2)],
            winding: if clockwise { Winding::Clockwise } else { Winding::CounterClockwise },
        }]);
        prop_assert_eq!(table.atoms[1].class, Some(class));
        prop_assert_eq!(table.atoms[1].implicit_hydrogens, Some(1));
    }
}

proptest! {
    // Branch choices, redundant markers, and global reversal preserve the same stereo.
    #[test]
    fn test_smiles_parse_branched_stereo(
        marked in any::<[bool; 4]>(), same in any::<bool>(),
        swap_left in any::<bool>(), swap_right in any::<bool>(), reverse in any::<bool>(),
    ) {
        let glyph = |present: bool, forward: bool| {
            if !present { "" } else if forward ^ reverse { "/" } else { "\\" }
        };
        let (left, branch_left) = if swap_left { ("Cl", "F") } else { ("F", "Cl") };
        let (branch_right, right) = if swap_right { ("I", "Br") } else { ("Br", "I") };
        let orientation = same ^ swap_left ^ swap_right;
        let input = format!("{left}{}C({}{branch_left})=C({}{branch_right}){}{right}",
            glyph(marked[0], true), glyph(marked[1], true),
            glyph(marked[2], !orientation), glyph(marked[3], orientation));
        let table = Smiles::parse(&input).unwrap().into_table_ir();
        let expected = if (marked[0] || marked[1]) && (marked[2] || marked[3]) {
            vec![StereoBond { bond: 2, configuration: BondConfiguration::Framed {
                references: [0, 4],
                relation: if orientation { BondRelation::SameSide } else { BondRelation::OppositeSide },
            }}]
        } else { vec![] };
        prop_assert_eq!(table.stereo_bonds, expected);
    }

    #[test]
    fn test_smiles_parse_with_redundant_stereo(
        same in any::<bool>(), reverse in any::<bool>(), copies in 1_usize..32,
    ) {
        let first = if reverse { '\\' } else { '/' };
        let last = if same ^ reverse { '\\' } else { '/' };
        let input = format!("F{first}C=C{last}F");
        let annotation = if same { "c:1" } else { "t:1" };
        let extended = format!("{input} |{}|", vec![annotation; copies].join(","));
        let expected = Smiles::parse(&input).unwrap();
        let actual = Smiles::parse_with(&extended, &SmilesIoConfig::chemaxon()).unwrap();
        prop_assert_eq!(actual.as_table_ir().stereo_bonds.as_slice(), expected.as_table_ir().stereo_bonds.as_slice());
    }

    #[test]
    fn test_smiles_parse_with_coordinates(
        offset in any::<[i16; 3]>(), signs in any::<[bool; 3]>(),
        cycle in 0_usize..3, scale in 1_u16..1000, same in any::<bool>(),
    ) {
        let points = [[0.,1.,0.], [0.,0.,0.], [2.,0.,0.], [2.,if same {1.} else {-1.},0.]];
        let coordinates = points.map(|p| {
            let q: [f64; 3] = array::from_fn(|i| f64::from(offset[i])
                + f64::from(scale) * p[(i + cycle) % 3] * if signs[i] {1.} else {-1.});
            format!("{},{},{}", q[0], q[1], q[2])
        }).join(";");
        let input = format!("FC=CF |({coordinates})|");
        let table = Smiles::parse_with(&input, &SmilesIoConfig::chemaxon()).unwrap().into_table_ir();
        prop_assert_eq!(table.stereo_bonds, vec![StereoBond { bond: 1,
            configuration: BondConfiguration::Framed { references: [0, 3],
                relation: if same { BondRelation::SameSide } else { BondRelation::OppositeSide },
            },
        }]);
    }

    // Frame transport and basic/extended preservation; coordinates are not a second
    // authority for an already published bond frame.
    #[test]
    fn test_table_molecule_try_into_ir_bond_frames(
        first_swap in any::<bool>(), second_swap in any::<bool>(), same in any::<bool>(),
        positions in any::<[[i16; 3]; 6]>(),
    ) {
        let mut table = Smiles::parse("FC(Cl)=C(Br)I").unwrap().into_table_ir();
        table.stereo_bonds = vec![StereoBond {
            bond: 2,
            configuration: BondConfiguration::Framed {
                references: [if first_swap { 2 } else { 0 }, if second_swap { 5 } else { 4 }],
                relation: if same { BondRelation::SameSide } else { BondRelation::OppositeSide },
            },
        }];
        let expected = u32::from(!same ^ first_swap ^ second_swap);
        let raised: Molecule = (&table).try_into_ir(&()).unwrap();
        prop_assert_eq!(raised.bond(BondId(2)).attributes.constraints.cis_trans_stereo(),
            Some(&CisTransStereoForm::stereo(StereoCoset::Lit(expected))));
        table.positions = Some(positions.map(|[x,y,z]| Point3D::new(f64::from(x),f64::from(y),f64::from(z))).to_vec());
        let converted = TableMolecule::try_from(ExtendedMolecule::from(table.clone())).unwrap();
        prop_assert_eq!(&converted, &table);
        prop_assert_eq!((&converted).try_into_ir(&()), Ok(raised));
    }

    #[test]
    fn test_table_molecule_try_into_ir_open_frames(
        bond in 0_u32..8, references in any::<[u32; 2]>(), either in any::<bool>(), duplicate in any::<bool>(),
    ) {
        let mut table = Smiles::parse("FC=CF").unwrap().into_table_ir();
        let frame = StereoBond { bond, configuration: if either { BondConfiguration::Either } else {
            BondConfiguration::Framed { references, relation: BondRelation::SameSide }
        }};
        table.stereo_bonds.push(frame.clone());
        if duplicate { table.stereo_bonds.push(frame); }
        let result = catch_unwind(|| { let result: Result<Molecule, _> = (&table).try_into_ir(&()); result });
        prop_assert!(result.is_ok());
    }
}

fn atom_label(table: &TableMolecule, atom: u32) -> u32 {
    table.atoms[atom as usize].class.unwrap_or(atom + 1)
}

fn bonds(table: &TableMolecule) -> BTreeMap<(u32, u32), BondOrder> {
    table
        .bonds
        .iter()
        .map(|bond| {
            let mut endpoints = [
                atom_label(table, bond.atoms.first()),
                atom_label(table, bond.atoms.second()),
            ];
            endpoints.sort_unstable();
            ((endpoints[0], endpoints[1]), bond.order)
        })
        .collect()
}

fn bond_stereo(table: &TableMolecule) -> BTreeMap<(u32, u32), bool> {
    table
        .stereo_bonds
        .iter()
        .map(|frame| {
            let bond = &table.bonds[frame.bond as usize];
            let endpoints = [bond.atoms.first(), bond.atoms.second()];
            let BondConfiguration::Framed {
                references,
                relation,
            } = frame.configuration
            else {
                panic!("generated definite frame")
            };
            let mut opposite = relation == BondRelation::OppositeSide;
            for side in 0..2 {
                let lowest = table
                    .bonds
                    .iter()
                    .filter_map(|bond| bond.atoms.other(endpoints[side]))
                    .filter(|&atom| atom != endpoints[1 - side])
                    .map(|atom| atom_label(table, atom))
                    .min()
                    .unwrap();
                opposite ^= atom_label(table, references[side]) != lowest;
            }
            let mut labels = endpoints.map(|atom| atom_label(table, atom));
            labels.sort_unstable();
            ((labels[0], labels[1]), opposite)
        })
        .collect()
}

fn atom_stereo(table: &TableMolecule) -> BTreeMap<u32, bool> {
    table
        .stereo_atoms
        .iter()
        .map(|frame| {
            let labels: Vec<_> = frame
                .ligands
                .iter()
                .map(|ligand| match ligand {
                    StereoLigand::Atom(atom) => u64::from(atom_label(table, *atom)),
                    StereoLigand::ImplicitHydrogen => 0,
                    StereoLigand::LonePair => u64::MAX,
                })
                .collect();
            let inversions = (0..labels.len())
                .map(|i| {
                    labels[i + 1..]
                        .iter()
                        .filter(|&&label| label < labels[i])
                        .count()
                })
                .sum::<usize>();
            (
                atom_label(table, frame.atom),
                (frame.winding == Winding::Clockwise) ^ (inversions % 2 != 0),
            )
        })
        .collect()
}

proptest! {
    // Independent simple edge sets cover arbitrary branches, components, and cycles.
    #[test]
    fn test_smiles_render_connectivity(
        count in 0_usize..16, pairs in prop::collection::vec((0_usize..16, 0_usize..16), 0..45),
        reverse_bonds in any::<bool>(),
    ) {
        let pairs: BTreeSet<_> = pairs.into_iter().filter(|&(a,b)| a < count && b < count && a != b)
            .map(|(a,b)| (a.min(b), a.max(b))).collect();
        let mut table = TableMolecule {
            atoms: (0..count).map(|i| Atom {
                class: Some(i as u32 + 1), implicit_hydrogens: Some(0),
                ..Atom::aliphatic_atom(Element::C)
            }).collect(),
            bonds: pairs.iter().map(|&(a,b)| Bond::new(a as u32,b as u32,BondOrder::Single)).collect(),
            ..TableMolecule::empty()
        };
        if reverse_bonds { table.bonds.reverse(); }
        let smiles = Smiles::from_table_ir(table.clone());
        prop_assert_eq!(smiles.as_table_ir(), &table);
        let text = smiles.render().unwrap();
        prop_assert_eq!(smiles.render().unwrap(), text.clone());
        let parsed = Smiles::parse(&text).unwrap();
        prop_assert_eq!(bonds(parsed.as_table_ir()), bonds(&table));
        prop_assert_eq!(parsed.as_table_ir().bonds.len(), table.bonds.len());
        prop_assert_eq!(parsed.as_table_ir().atoms.len(), table.atoms.len());
        let labels: BTreeSet<_> = parsed.as_table_ir().atoms.iter().map(|atom| atom.class.unwrap()).collect();
        prop_assert_eq!(labels, (1..=count as u32).collect::<BTreeSet<_>>());
        let normalized = parsed.render().unwrap();
        let rerendered = Smiles::parse(&normalized).unwrap();
        prop_assert_eq!(bonds(rerendered.as_table_ir()), bonds(&table));
        prop_assert_eq!(rerendered.render().unwrap(), normalized);
        prop_assert_eq!(
            parsed.as_table_ir().atoms.iter().map(|atom| atom.class).collect::<Vec<_>>(),
            rerendered.as_table_ir().atoms.iter().map(|atom| atom.class).collect::<Vec<_>>(),
        );
    }

    #[test]
    fn test_smiles_render_atom_stereo(
        input in prop::sample::select(vec!["[C@](F)(Cl)(Br)I", "[C@H](F)(Cl)Br",
            "[H][C@@](F)(Cl)Br", "C[S@](=O)CC", "C[C@H]1CCCCO1",
            "O1CCC[C@](F)1Cl", "[C@]12(CCC1)CCC2", "O1CCC[C@]21CCCC2"]),
        flip in any::<bool>(), rotation in 0_usize..4,
    ) {
        let mut table = Smiles::parse(input).unwrap().into_table_ir();
        for (i, atom) in table.atoms.iter_mut().enumerate() {
            // A class requires brackets, so keep inferred-H atoms unlabelled; these fixtures
            // preserve atom discovery order. Label only atoms already carrying fixed H.
            if atom.implicit_hydrogens.is_some() { atom.class = Some(i as u32 + 1); }
        }
        for frame in &mut table.stereo_atoms {
            frame.ligands.rotate_left(rotation);
            if flip { frame.winding = Winding::Clockwise; }
        }
        let expected = atom_stereo(&table);
        let smiles = Smiles::from_table_ir(table.clone());
        let text = smiles.render().unwrap();
        let parsed = Smiles::parse(&text).unwrap();
        prop_assert_eq!(atom_stereo(parsed.as_table_ir()), expected);
        prop_assert_eq!(bonds(parsed.as_table_ir()), bonds(&table));
        prop_assert_eq!(parsed.render().unwrap(), text);
    }

    #[test]
    fn test_smiles_render_bond_chain(
        markers in prop::collection::vec(prop::option::of(any::<bool>()), 2..12),
    ) {
        let mut input = String::from("F");
        for (i, marker) in markers.iter().enumerate() {
            input.push_str(match marker { None => "", Some(true) => "/", Some(false) => "\\" });
            input.push_str(if i+1 == markers.len() { "F" } else { "C=C" });
        }
        let table = Smiles::parse(&input).unwrap().into_table_ir();
        let expected = bond_stereo(&table);
        let smiles = Smiles::from_table_ir(table);
        let text = smiles.render().unwrap();
        let parsed = Smiles::parse(&text).unwrap();
        prop_assert_eq!(bond_stereo(parsed.as_table_ir()), expected);
        prop_assert_eq!(parsed.render().unwrap(), text);
    }

    #[test]
    fn test_smiles_render_with(
        input in prop::sample::select(vec!["N->[Cu]", "C~C", "[te]", "[siH]", "[seH]", "CCO"]),
    ) {
        let config = SmilesIoConfig::lenient();
        let table = Smiles::parse_with(input, &config).unwrap().into_table_ir();
        let smiles = Smiles::from_table_ir(table.clone());
        let text = smiles.render_with(&config).unwrap();
        prop_assert_eq!(Smiles::parse_with(&text, &config).unwrap().render_with(&config).unwrap(), text);
        prop_assert_eq!(smiles.as_table_ir(), &table);
    }

    #[test]
    fn test_smiles_render_atom_fields(
        mass in prop::option::of(0_u32..1000), charge in -99_i8..100,
        hydrogens in 0_u8..10, class in any::<u32>(),
    ) {
        let atom = Atom { isotope_mass: mass, charge: Some(charge), implicit_hydrogens: Some(hydrogens),
            class: Some(class), ..Atom::aliphatic_atom(Element::C) };
        let smiles = Smiles::from_table_ir(TableMolecule { atoms: vec![atom.clone()], ..TableMolecule::empty() });
        let text = smiles.render().unwrap();
        let parsed = Smiles::parse(&text).unwrap();
        let actual = &parsed.as_table_ir().atoms[0];
        prop_assert_eq!((actual.element, actual.isotope_mass, actual.charge.unwrap_or(0), actual.implicit_hydrogens, actual.class),
            (atom.element, mass, charge, Some(hydrogens), Some(class)));
    }

    #[test]
    fn test_smiles_render_electrons(
        components in vec((select(vec![
            ("C", 0, 0, SpinMultiplicity::SINGLET),
            ("N", 1, 0, SpinMultiplicity::SINGLET),
            ("O", 2, 0, SpinMultiplicity::SINGLET),
            ("[CH4]", 0, 0, SpinMultiplicity::SINGLET),
            ("[13CH4]", 0, 0, SpinMultiplicity::SINGLET),
            ("[CH3]", 0, 1, SpinMultiplicity::DOUBLET),
            ("[CH2]", 0, 2, SpinMultiplicity::TRIPLET),
            ("[O]", 2, 2, SpinMultiplicity::TRIPLET),
            ("[NH4+]", 0, 0, SpinMultiplicity::SINGLET),
            ("[Cl-]", 4, 0, SpinMultiplicity::SINGLET),
            ("[He]", 1, 0, SpinMultiplicity::SINGLET),
        ]), any::<[bool; 3]>()), 1..32),
    ) {
        let input = components.iter().map(|((text, ..), _)| *text).collect::<Vec<_>>().join(".");
        let parsed = Smiles::parse(&input).unwrap();
        let mut table = parsed.clone().into_table_ir();
        for (atom, ((_, lone_pairs, unpaired, multiplicity), present)) in table.atoms.iter_mut().zip(components) {
            atom.lone_pairs = present[0].then_some(lone_pairs);
            atom.unpaired_electrons = present[1].then_some(unpaired);
            atom.multiplicity = present[2].then_some(multiplicity);
        }
        let smiles = Smiles::from_table_ir(table.clone());
        let text = smiles.render().unwrap();
        prop_assert_eq!(&text, &input);
        prop_assert_eq!(Smiles::parse(&text).unwrap(), parsed);
        prop_assert_eq!(smiles.into_table_ir(), table);
    }
}

// Every possible absent/slash/backslash assignment on the written single bonds is parsed.
// The oracle contains no production marker-selection or parity-constraint implementation.
#[rstest]
#[case::chain("C{}C=C{}C=C{}C=C{}C")]
#[case::branched("F{}C({}Cl)=C({}F){}C({}F)=C({}Cl){}Br")]
#[case::branched_system("C{}C=C({}C=C{}C){}C=C{}C")]
#[case::cycle("C1=C{}C=C{}C=C{}1")]
#[case::substituted_cycle("C1=C({}F){}C=C({}Cl){}C=C{}1")]
fn test_smiles_render_bond_stereo(#[case] template: &str) {
    let parts: Vec<_> = template.split("{}").collect();
    let candidates = parts.len() - 1;
    let mut accepted = BTreeSet::new();
    for assignment in 0..3_usize.pow(candidates as u32) {
        let mut assignment = assignment;
        let mut text = parts[0].to_owned();
        for part in &parts[1..] {
            text.push_str(["", "/", "\\"][assignment % 3]);
            assignment /= 3;
            text.push_str(part);
        }
        if let Ok(smiles) = Smiles::parse(&text) {
            accepted.insert(bond_stereo(smiles.as_table_ir()));
        }
    }
    let base = Smiles::parse(&parts.concat()).unwrap().into_table_ir();
    let sites: Vec<_> = base
        .bonds
        .iter()
        .enumerate()
        .filter(|(_, b)| b.order == BondOrder::Double)
        .map(|(i, b)| {
            let endpoints = [b.atoms.first(), b.atoms.second()];
            let references = [0, 1].map(|side| {
                base.bonds
                    .iter()
                    .filter_map(|b| b.atoms.other(endpoints[side]))
                    .filter(|&a| a != endpoints[1 - side])
                    .min()
                    .unwrap()
            });
            (i as u32, references)
        })
        .collect();
    for states in 0..3_usize.pow(sites.len() as u32) {
        let mut states = states;
        let mut table = base.clone();
        for &(bond, references) in &sites {
            let state = states % 3;
            states /= 3;
            if state != 0 {
                table.stereo_bonds.push(StereoBond {
                    bond,
                    configuration: BondConfiguration::Framed {
                        references,
                        relation: if state == 1 {
                            BondRelation::SameSide
                        } else {
                            BondRelation::OppositeSide
                        },
                    },
                });
            }
        }
        let expected = bond_stereo(&table);
        let result = Smiles::from_table_ir(table).render();
        assert_eq!(
            result.is_ok(),
            accepted.contains(&expected),
            "{template}: {expected:?}"
        );
        match result {
            Ok(text) => {
                let parsed = Smiles::parse(&text).unwrap();
                assert_eq!(bond_stereo(parsed.as_table_ir()), expected, "{text}");
            }
            Err(error) => assert_eq!(
                error,
                SmilesRenderError::NoMarkerAssignment {
                    bond: sites
                        .iter()
                        .find(|(bond, _)| {
                            let endpoints = base.bonds[*bond as usize].atoms;
                            expected.contains_key(&(endpoints.first() + 1, endpoints.second() + 1))
                        })
                        .unwrap()
                        .0,
                }
            ),
        }
    }
}

#[rstest]
#[case::ring_labels(
    "[C:1][C:8]12[C:6][C:11]2[C:9]1[C:12]",
    "[C:1][C:8]12[C:6][C:11]1[C:9]2[C:12]"
)]
fn test_smiles_render_idempotence(#[case] input: &str, #[case] expected: &str) {
    let parsed = Smiles::parse(input).unwrap();
    let output = parsed.render().unwrap();
    assert_eq!(output, expected);
    let reparsed = Smiles::parse(&output).unwrap();
    assert_eq!(reparsed.render().unwrap(), expected);
    assert_eq!(bonds(reparsed.as_table_ir()), bonds(parsed.as_table_ir()));
    for table in [parsed.as_table_ir(), reparsed.as_table_ir()] {
        assert_eq!(
            table
                .atoms
                .iter()
                .map(|atom| atom.class)
                .collect::<Vec<_>>(),
            vec![Some(1), Some(8), Some(6), Some(11), Some(9), Some(12)]
        );
    }
}

proptest! {
    #[test]
    fn test_reaction_smiles_render_roundtrip(
        sections in vec(vec(select(vec!["C", "N", "c1ccccc1", "F/C=C/F", "F/C=C\\F", "N[C@H](F)Cl", "[C:7]", "[O:19]", "[H]C"]), 0..5), 3),
    ) {
        let input = sections.iter().map(|section| section.join(".")).collect::<Vec<_>>().join(">");
        let parsed = ReactionSmiles::parse(&input).unwrap();
        let output = parsed.render().unwrap();
        prop_assert_eq!(&parsed.render().unwrap(), &output);
        let reparsed = ReactionSmiles::parse(&output).unwrap();
        prop_assert_eq!(reparsed.render().unwrap(), output);
        let original = parsed.as_table_ir();
        let result = reparsed.as_table_ir();
        prop_assert_eq!(&original.atom_mapping, &result.atom_mapping);
        for (before, after) in [(&original.reactants, &result.reactants), (&original.agents, &result.agents), (&original.products, &result.products)] {
            prop_assert_eq!(
                before.atoms.iter().map(|atom| (atom.element, atom.class, atom.implicit_hydrogens, atom.aromatic)).collect::<Vec<_>>(),
                after.atoms.iter().map(|atom| (atom.element, atom.class, atom.implicit_hydrogens, atom.aromatic)).collect::<Vec<_>>(),
            );
            prop_assert_eq!(bonds(before), bonds(after));
            prop_assert_eq!(atom_stereo(before), atom_stereo(after));
            prop_assert_eq!(bond_stereo(before), bond_stereo(after));
        }
    }
}
