//! Atom block parser for CTab files.

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{alpha1, space0};
use nom::combinator::{all_consuming, complete, map, map_parser, map_res, opt, value};
use nom::error;
use nom::sequence::{delimited, preceded};
use nom::Parser;
use umol_data::{Element, NamedIsotope};

use super::utils::{fixed_width_float, fixed_width_int, fixed_width_int_in_range_minus1};

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum AtomSymbol {
    Element(Element),
    AtomList,
    Unspecified(char),
    LonePair,
    RGroup(usize),
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AtomLine {
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

/// Parse atom symbol and named isotope
///
///
///
fn atom_symbol<'a>(
) -> impl Parser<&'a [u8], Output = (AtomSymbol, Option<i8>), Error = error::Error<&'a [u8]>> {
    |input| {
        map_parser(
            take(3usize),
            all_consuming(alt((
                value(
                    (AtomSymbol::LonePair, None),
                    delimited(space0, tag("LP"), space0),
                ),
                value(
                    (AtomSymbol::AtomList, None),
                    delimited(space0, tag("L"), space0),
                ),
                value(
                    (AtomSymbol::Unspecified('A'), None),
                    delimited(space0, tag("A"), space0),
                ),
                value(
                    (AtomSymbol::Unspecified('Q'), None),
                    delimited(space0, tag("Q"), space0),
                ),
                value(
                    (AtomSymbol::Unspecified('*'), None),
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
                    |(_, idx)| (AtomSymbol::RGroup(idx), None),
                ),
                map(
                    map_res(delimited(space0, alpha1, space0), |s| {
                        Element::from_symbol_bytes(s)
                            .ok_or_else(|| error::Error::new(s, error::ErrorKind::MapRes))
                    }),
                    |element| (AtomSymbol::Element(element), None),
                ),
                map(
                    map_res(delimited(space0, alpha1, space0), |s| {
                        NamedIsotope::from_symbol_bytes(s)
                            .map(|iso| (iso.element(), iso.mass_number()))
                            .ok_or_else(|| error::Error::new(s, error::ErrorKind::MapRes))
                    }),
                    |(element, mass_number)| {
                        (
                            AtomSymbol::Element(element),
                            Some(mass_number as i8 - element.reference_mass_number() as i8),
                        )
                    },
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
pub(crate) fn atom_line<'a>(
) -> impl Parser<&'a [u8], Output = AtomLine, Error = error::Error<&'a [u8]>> {
    let x = fixed_width_float::<f64>(10, 4);
    let y = fixed_width_float::<f64>(10, 4);
    let z = fixed_width_float::<f64>(10, 4);
    let symbol = atom_symbol();
    let mass_diff = fixed_width_int::<i8>(2);
    let charge_code = fixed_width_int::<u8>(3);
    let stereo_parity = fixed_width_int::<u8>(3);
    let hydrogen_code = fixed_width_int::<u8>(3);
    let stereo_care = fixed_width_int::<u8>(3);
    let valence_code = fixed_width_int::<u8>(3);
    let atom_mapping = fixed_width_int::<u8>(3);
    let inversion = fixed_width_int::<u8>(3);
    let exact_change = fixed_width_int::<u8>(3);
    all_consuming(map(
        (
            x,
            y,
            z,
            take(1usize),
            symbol,
            mass_diff,
            charge_code,
            stereo_parity,
            hydrogen_code,
            stereo_care,
            valence_code,
            opt(preceded(
                take(9usize),
                complete((
                    atom_mapping,
                    opt(complete((inversion, opt(complete(exact_change))))),
                )),
            )),
            space0,
        ),
        |(
            x,
            y,
            z,
            _,
            (symbol, mass_diff_named),
            mass_diff,
            charge_code,
            stereo_parity,
            hydrogen_code,
            stereo_care,
            valence_code,
            rest,
            _,
        )| {
            let (atom_mapping, inversion, exact_change) = match rest {
                Some((m, Some((n, Some(e))))) => (m, n, e),
                Some((m, Some((n, None)))) => (m, n, 0),
                Some((m, None)) => (m, 0, 0),
                None => (0, 0, 0),
            };
            AtomLine {
                x,
                y,
                z,
                symbol,
                mass_diff: mass_diff_named.unwrap_or(mass_diff),
                charge_code,
                stereo_parity,
                hydrogen_code,
                stereo_care,
                valence_code,
                atom_mapping,
                inversion,
                exact_change,
            }
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::error::ErrorKind;
    use nom::Err;
    use rstest::rstest;

    #[rstest]
    #[case(b"A  ", AtomSymbol::Unspecified('A'), None)]
    #[case(b"Q  ", AtomSymbol::Unspecified('Q'), None)]
    #[case(b"*  ", AtomSymbol::Unspecified('*'), None)]
    #[case(b"L  ", AtomSymbol::AtomList, None)]
    #[case(b"LP ", AtomSymbol::LonePair, None)]
    #[case(b"R1 ", AtomSymbol::RGroup(0), None)]
    #[case(b"R3 ", AtomSymbol::RGroup(2), None)]
    #[case(b"H  ", AtomSymbol::Element(Element::H), None)]
    #[case(b"C  ", AtomSymbol::Element(Element::C), None)]
    #[case(b" C ", AtomSymbol::Element(Element::C), None)]
    #[case(b"  C", AtomSymbol::Element(Element::C), None)]
    #[case(b"Cu ", AtomSymbol::Element(Element::Cu), None)]
    #[case(b"cu ", AtomSymbol::Element(Element::Cu), None)]
    #[case(b"CU ", AtomSymbol::Element(Element::Cu), None)]
    #[case(b"D  ", AtomSymbol::Element(Element::H), Some(1))]
    #[case(b"d  ", AtomSymbol::Element(Element::H), Some(1))]
    #[case(b"T  ", AtomSymbol::Element(Element::H), Some(2))]
    fn test_atom_symbol(
        #[case] input: &[u8],
        #[case] expected_symbol: AtomSymbol,
        #[case] expected_mass_diff: Option<i8>,
    ) {
        let (remaining, (symbol, mass_diff)) = atom_symbol().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(symbol, expected_symbol);
        assert_eq!(mass_diff, expected_mass_diff);
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
        let res = atom_symbol().parse(input);
        assert!(res.is_err(), "{}", desc);
        assert!(
            matches!(res.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            res.clone().unwrap_err().map(|e| e.code),
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
}
