//! Property-based coverage for the SMILES parser.

use std::panic::catch_unwind;

use proptest::collection::vec;
use proptest::prelude::*;
use proptest::sample::select;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_geometric_core::Point3D;
use umol_graph_ir::ir::{
    BondId, CisTransStereoForm, ElementForm, Molecule, StereoCoset, TryIntoIr,
};
use umol_io::smiles::config::SmilesIoConfig;
use umol_io::smiles::{parse_extended_smiles_bytes, ParseError, Smiles};
use umol_io::table_ir::{
    BondConfiguration, BondRelation, ExtendedMolecule, Molecule as TableMolecule, Span, StereoAtom,
    StereoBond, StereoLigand, Winding,
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
