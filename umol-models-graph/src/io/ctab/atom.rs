//! Atom block parser for CTab files.

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{alpha1, space0};
use nom::combinator::{all_consuming, complete, map, map_parser, map_res, opt, value};
use nom::error;
use nom::sequence::{delimited, preceded};
use nom::Parser;
use umol_data::{Element, NamedIsotope};

use crate::atom::{Atom, AtomLike, AtomList, AtomSymbol};
use crate::conformer::Point3D;

use super::convert::{
    convert_atom_charge_code, convert_atom_hydrogen_count_code, convert_atom_mass_diff_code,
    convert_atom_stereo_parity_code, convert_atom_valence_code,
};
use super::utils::{
    fixed_width_float, fixed_width_int, fixed_width_int_in_range, fixed_width_int_in_range_minus1,
};

#[derive(Debug, Clone)]
pub(crate) struct AtomInput {
    x: f64,
    y: f64,
    z: f64,
    symbol: AtomSymbol, // 'aaa' field: atom symbol (see AtomSymbol enum)
    mass_diff: i8, // 'dd' field: mass difference (-3, -2, -1, 0, 1, 2, 3, 4), 0 if value outside of this range
    charge_code: u8, // 'ccc' field: 0 = uncharged, 1 = +3, 2 = +2, 3 = +1, 4 = doublet radical, 5 = -1, 6 = -2, 7 = -3
    stereo_parity: u8, // 'sss' field: 0 = not stereo, 1 = odd, 2 = even, 3 = either or unmarked
    hydrogen_code: u8, // 'hhh' field: 1 = H0, 2 = H1, 3 = H2, 4 = H3, 5 = H4
    stereo_care: u8, // 'bbb' field: 0 = ignore stereo, 1 = stereo in query must match
    valence_code: u8, // 'vvv' field: 0 = default, 1-14 = explicit, 15 = explicit 0
    // Skipping obsolete or unused fields (HHH, rrr, iii)
    atom_mapping: u8, // 'mmm' field: 1..=number of atoms
    inversion: u8,    // 'nnn' field: 0 = property not applied, 2 = inverted, 3 = retained
    exact_change: u8, // 'eee' field: 0 = property not applied, 1 = charge in query must match
}

/// Parse atom symbol
/// Values: entry in periodic table or L for atom list, A, Q, * for unspecified atom, and LP for lone pair,
/// or R# for Rgroup label. Named isotopes (D, T) are supported as extension
fn atom_symbol<'a>() -> impl Parser<&'a [u8], Output = AtomSymbol, Error = error::Error<&'a [u8]>> {
    |input| {
        map_parser(
            take(3usize),
            all_consuming(alt((
                value(AtomSymbol::LonePair, delimited(space0, tag("LP"), space0)),
                value(
                    AtomSymbol::AtomList(AtomList { elements: vec![] }),
                    delimited(space0, tag("L"), space0),
                ),
                value(
                    AtomSymbol::Unspecified('A'),
                    delimited(space0, tag("A"), space0),
                ),
                value(
                    AtomSymbol::Unspecified('Q'),
                    delimited(space0, tag("Q"), space0),
                ),
                value(
                    AtomSymbol::Unspecified('*'),
                    delimited(space0, tag("*"), space0),
                ),
                map(
                    delimited(
                        space0,
                        (
                            tag("R"),
                            fixed_width_int_in_range_minus1::<usize, _>(2, 1..=31),
                        ),
                        space0,
                    ),
                    |(_, idx)| AtomSymbol::RGroup(idx),
                ),
                map(
                    map_res(delimited(space0, alpha1, space0), |s| {
                        Element::from_symbol_bytes(s)
                            .ok_or_else(|| error::Error::new(s, error::ErrorKind::MapRes))
                    }),
                    |element| AtomSymbol::Element(element),
                ),
                map(
                    map_res(delimited(space0, alpha1, space0), |s| {
                        NamedIsotope::from_symbol_bytes(s)
                            .ok_or_else(|| error::Error::new(s, error::ErrorKind::MapRes))
                    }),
                    |isotope| AtomSymbol::NamedIsotope(isotope),
                ),
            ))),
        )
        .parse(input)
    }
}

/// Parse atom line
/// xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee (69 characters wide)
///
/// *Values in the atom block*
/// -----------------------------------------------------------------------------------------
/// | Field | Meaning            | Values       | Notes                                     |
/// |-------|--------------------|--------------|-------------------------------------------|
/// | x,y,z | atom coordinates   |              | *[Generic]*, F10.4 format                 |
/// | aaa   | atom symbol        | s. above     | *[Generic, Query, 3D, RGroup]*            |
/// | dd    | mass difference    | -3..=4       | *[Generic]*, s. also M  ISO               |
/// | ccc   | charge code        | 0..=7        | *[Generic]*, s. also M  CHG/M  RAD        |
/// | sss   | stereo parity      | 0..=3        | *[Generic]*, ignored when read            |
/// | hhh   | hydrogen code      | 0..=5        | *[Query]*, H0 means no implicit Hs        |
/// |       |                    |              | Hn means >=n implicit Hs                  |
/// | bbb   | stereo care        | 0, 1         | *[Query]*, consider double bond stereo    |
/// |       |                    |              | when stereo care is 1 for both bond atoms |
/// | vvv   | valence code       | 0..=15       | *[Generic]*                               |
/// | mmm   | atom mapping       | 1..=#atoms   | *[Reaction]*                              |
/// | nnn   | inversion          | 0..=2        | *[Reaction]*                              |
/// | eee   | exact change       | 0, 1         | *[Reaction]*                              |
/// -----------------------------------------------------------------------------------------
///
pub(crate) fn atom_input<'a>(
) -> impl Parser<&'a [u8], Output = (Atom, Point3D), Error = error::Error<&'a [u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol();
    let mass_diff = map_res(
        fixed_width_int_in_range::<i8, _>(2, -3..=4),
        convert_atom_mass_diff_code,
    );
    let charge = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=7),
        convert_atom_charge_code,
    );
    let stereo_parity = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=3),
        convert_atom_stereo_parity_code,
    );
    let hydrogen_count = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=5),
        convert_atom_hydrogen_count_code,
    );
    let stereo_care = fixed_width_int_in_range_minus1::<u8, _>(3, 0..=1);
    let valence = map_res(
        fixed_width_int_in_range::<u8, _>(3, 0..=15),
        convert_atom_valence_code,
    );
    let atom_mapping = fixed_width_int::<u8>(3);
    let inversion = fixed_width_int_in_range::<u8, _>(3, 0..=2);
    let exact_change = fixed_width_int_in_range::<u8, _>(3, 0..=1);

    all_consuming(map(
        (
            x,
            y,
            z,
            take(1usize),
            symbol,
            opt(complete(mass_diff)),
            opt(complete(charge)),
            opt(complete(stereo_parity)),
            opt(complete(hydrogen_count)),
            opt(complete(stereo_care)),
            opt(complete(valence)),
            opt(complete(preceded(take(9usize), atom_mapping))),
            opt(complete(inversion)),
            opt(complete(exact_change)),
            space0,
        ),
        |(
            x,
            y,
            z,
            _,
            symbol,
            mass_diff,
            charge,
            stereo_parity,
            hydrogen_count,
            stereo_care,
            valence,
            atom_mapping,
            inversion,
            exact_change,
            _,
        )| {
            let atom = Atom {
                // element: match symbol {
                //     AtomSymbol::Element(e) => e,
                //     AtomSymbol::NamedIsotope(i) => i.element,
                //     AtomSymbol::AtomList(l) => l.elements[0],
                //     AtomSymbol::LonePair => Element::H,
                //     AtomSymbol::RGroup(i) => Element::H,
                // },
                // isotope_mass: mass_diff,
                // charge: charge.unwrap_or(0),
                // stereo_parity: stereo_parity.unwrap_or(AtomStereoParity::None),
                // hydrogen_count: hydrogen_count.unwrap_or(0),
                // valence: valence.unwrap_or(0),
                // atom_map_num: atom_mapping.unwrap_or(0),
                // radical: inversion.unwrap_or(0),
                // properties: std::collections::HashMap::new(),
            };
            (atom, Point3D::new(x, y, z))
        },
    ))
}

#[cfg(test)]
mod tests {
    use crate::AtomStereoParity;

    use super::*;
    use nom::{error::ErrorKind, Err};
    use rstest::rstest;

    #[rstest]
    #[case(b"A  ", AtomSymbol::Unspecified('A'))]
    #[case(b"Q  ", AtomSymbol::Unspecified('Q'))]
    #[case(b"*  ", AtomSymbol::Unspecified('*'))]
    #[case(b"L  ", AtomSymbol::AtomList(AtomList { elements: vec![] }))]
    #[case(b"LP ", AtomSymbol::LonePair)]
    #[case(b"R1 ", AtomSymbol::RGroup(0))]
    #[case(b"R3 ", AtomSymbol::RGroup(2))]
    #[case(b"H  ", AtomSymbol::Element(Element::H))]
    #[case(b"C  ", AtomSymbol::Element(Element::C))]
    #[case(b" C ", AtomSymbol::Element(Element::C))]
    #[case(b"  C", AtomSymbol::Element(Element::C))]
    #[case(b"Cu ", AtomSymbol::Element(Element::Cu))]
    #[case(b"cu ", AtomSymbol::Element(Element::Cu))]
    #[case(b"CU ", AtomSymbol::Element(Element::Cu))]
    #[case(b"D  ", AtomSymbol::NamedIsotope(NamedIsotope::D))]
    #[case(b"d  ", AtomSymbol::NamedIsotope(NamedIsotope::D))]
    #[case(b"T  ", AtomSymbol::NamedIsotope(NamedIsotope::T))]
    fn test_atom_symbol(#[case] input: &[u8], #[case] expected: AtomSymbol) {
        let (remaining, symbol) = atom_symbol().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(symbol, expected);
    }

    #[rstest]
    #[case(b"   ", "empty field", ErrorKind::Alpha)]
    #[case(b"H", "too short", ErrorKind::Eof)]
    #[case(b"R  ", "R group index missing", ErrorKind::MapRes)]
    #[case(b"R0 ", "R group index must be between 1 and 31", ErrorKind::MapRes)]
    #[case(b"R32", "R group index must be between 1 and 31", ErrorKind::MapRes)]
    #[case(b"Xx ", "Unknown atom symbol", ErrorKind::MapRes)]
    #[case(b"LQ ", "Unknown atom symbol", ErrorKind::Eof)]
    fn test_atom_symbol_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = atom_symbol().parse(input);
        assert!(result.is_err(), "{} should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    // From CTab spec (Figure 3)
    #[case(b"   -0.6622    0.5342    0.0000 C   0  0  2  0  0  0",
      AtomLine { x: -0.6622, y: 0.5342, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 2,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"    0.6622   -0.3000    0.0000 C   0  0  0  0  0  0",
      AtomLine { x: 0.6622, y: -0.3, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"   -0.7207    2.0817    0.0000 C   1  0  0  0  0  0",
      AtomLine { x: -0.7207, y: 2.0817, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 1, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"   -1.8622   -0.3695    0.0000 N   0  3  0  0  0  0",
      AtomLine { x: -1.8622, y: -0.3695, z: 0.0, symbol: AtomSymbol::Element(Element::N), mass_diff: 0, charge_code: 3, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"    0.6220   -1.8037    0.0000 O   0  0  0  0  0  0",
      AtomLine { x: 0.622, y: -1.8037, z: 0.0, symbol: AtomSymbol::Element(Element::O), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"    1.9464    0.4244    0.0000 O   0  5  0  0  0  0",
      AtomLine { x: 1.9464, y: 0.4244, z: 0.0, symbol: AtomSymbol::Element(Element::O), mass_diff: 0, charge_code: 5, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"    0.6622   -0.3000    0.0000   C 0  0  0  0  0  0  0  0  0  0  0  0",
      AtomLine { x: 0.6622, y: -0.3, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"    0.6622   -0.3000    0.0000   C 0  0  0  0  0",
      AtomLine { x: 0.6622, y: -0.3, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"    0.6622   -0.3000    0.0000   C",
      AtomLine { x: 0.6622, y: -0.3, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    // From RDKit test files
    #[case(b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
      AtomLine { x: 0.0, y: 0.0, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"   -3.0000   -7.8750    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0",
      AtomLine { x: -3.0, y: -7.875, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 1, inversion: 0, exact_change: 0 })]
    #[case(b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
      AtomLine { x: 0.0, y: 0.0, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"   -3.0000   -7.8750    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0",
      AtomLine { x: -3.0, y: -7.875, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 1, inversion: 0, exact_change: 0 })]
    #[case(b"   -1.5711   -7.8750    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0",
      AtomLine { x: -1.5711, y: -7.875, z: 0.0, symbol: AtomSymbol::Element(Element::Cl), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"   -1.2375    2.1509    0.0000 C   0  0  0  0  0  1  0  0  0  0  0  0",
      AtomLine { x: -1.2375, y: 2.1509, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 1, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"    3.6666   -2.5791    0.0000 C   0  0  3  0  0  0  0  0  0  0  0  0",
      AtomLine { x: 3.6666, y: -2.5791, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 3,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"    0.0895   -0.4313    0.0000 C   0  0  1  0  0  0  0  0  0  0  0  0",
      AtomLine { x: 0.0895, y: -0.4313, z: 0.0, symbol: AtomSymbol::Element(Element::C), mass_diff: 0, charge_code: 0, stereo_parity: 1,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    #[case(b"    0.7143   -0.2061    0.0000 D   0  0  0  0  0  0  0  0  0  0  0  0",
      AtomLine { x: 0.7143, y: -0.2061, z: 0.0, symbol: AtomSymbol::Element(Element::H), mass_diff: 1, charge_code: 0, stereo_parity: 0,
                 hydrogen_code: 0, stereo_care: 0, valence_code: 0, atom_mapping: 0, inversion: 0, exact_change: 0 })]
    fn test_atom_line(#[case] input: &[u8], #[case] expected_atom_line: AtomLine) {
        let (remaining, atom_line) = atom_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(atom_line, expected_atom_line);
    }

    #[rstest]
    #[case(b"   -0.6622    0.5342    0.0000 C   0  0  2  0  0  0",
           Some(Atom { element: Element::C, charge: 0, isotope_mass: None, stereo_parity: Some(AtomStereoParity::Even), hydrogen_count: Some(0), valence: Some(0), atom_map_num: Some(0), radical: None, properties: std::collections::HashMap::new() }),
           None,
           Point3D::new(-0.6622, 0.5342, 0.0))]
    #[case(b"    0.6622   -0.3000    0.0000 C   0  0  0  0  0  0",
           Some(Atom { element: Element::C, charge: 0, isotope_mass: None, stereo_parity: None, hydrogen_count: Some(0), valence: Some(0), atom_map_num: Some(0), radical: None, properties: std::collections::HashMap::new() }),
           None,
           Point3D::new(0.6622, -0.3, 0.0))]
    #[case(b"   -0.7207    2.0817    0.0000 C   1  0  0  0  0  0",
           Some(Atom { element: Element::C, charge: 0, isotope_mass: Some(13), stereo_parity: None, hydrogen_count: Some(0), valence: Some(0), atom_map_num: Some(0), radical: None, properties: std::collections::HashMap::new() }),
           None,
           Point3D::new(-0.7207, 2.0817, 0.0))]
    #[case(b"   -1.8622   -0.3695    0.0000 N   0  3  0  0  0  0",
           Some(Atom { element: Element::N, charge: 1, isotope_mass: None, stereo_parity: None, hydrogen_count: Some(0), valence: Some(0), atom_map_num: Some(0), radical: None, properties: std::collections::HashMap::new() }),
           None,
           Point3D::new(-1.8622, -0.3695, 0.0))]
    #[case(
        b"    0.0000    0.0000    0.0000 L   0  0  0  0  0  0",
        None,
        Some(AtomSymbol::AtomList(AtomList { elements: vec![] })),
        Point3D::new(0.0, 0.0, 0.0)
    )]
    #[case(
        b"    0.0000    0.0000    0.0000 A   0  0  0  0  0  0",
        None,
        Some(AtomSymbol::Unspecified('A')),
        Point3D::new(0.0, 0.0, 0.0)
    )]
    #[case(
        b"    0.0000    0.0000    0.0000 LP  0  0  0  0  0  0",
        None,
        Some(AtomSymbol::LonePair),
        Point3D::new(0.0, 0.0, 0.0)
    )]
    #[case(
        b"    0.0000    0.0000    0.0000 R1  0  0  0  0  0  0",
        None,
        Some(AtomSymbol::RGroup(0)),
        Point3D::new(0.0, 0.0, 0.0)
    )]
    #[case(b"    0.7143   -0.2061    0.0000 D   0  0  0  0  0  0",
           Some(Atom { element: Element::H, charge: 0, isotope_mass: Some(2), stereo_parity: None, hydrogen_count: Some(0), valence: Some(0), atom_map_num: Some(0), radical: None, properties: std::collections::HashMap::new() }),
           None,
           Point3D::new(0.7143, -0.2061, 0.0))]
    fn test_atom_line_parsed(
        #[case] input: &[u8],
        #[case] expected_atom: Option<Atom>,
        #[case] expected_atom_like: Option<AtomSymbol>,
        #[case] expected_position: Point3D,
    ) {
        let (remaining, (atom, atom_like, position)) = atom_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(atom, expected_atom);
        assert_eq!(atom_like, expected_atom_like);
        assert_eq!(position, expected_position);
    }
}
