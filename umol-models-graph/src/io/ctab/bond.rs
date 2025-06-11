//! Bond block parser for CTab files.

use super::utils::{fixed_width_int, fixed_width_int_minus1};
use nom::bytes::complete::take;
use nom::character::complete::space0;
use nom::combinator::{all_consuming, complete, map, opt};
use nom::error;
use nom::sequence::preceded;
use nom::Parser;

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct BondLine {
    first_atom: usize,
    second_atom: usize,
    bond_type: u8,
    bond_stereo: u8,
    bond_topology: u8,
    bond_reacting_center: u8,
}

/// Parse bond line
/// 111222tttsssxxxrrrccc (21 characters wide)
///
/// *Values in the bond block*
/// ------------------------------------------------------------------
/// | Field | Meaning              | Values     | Notes              |
/// |-------|----------------------|------------|--------------------|
/// | 111   | first atom           | 1..=aaa    | *[Generic]*        |
/// | 222   | second atom          | 1..=aaa    | *[Generic]*        |
/// | ttt   | bond type            | 1..=8      | *[Query]*          |
/// | sss   | bond stereo          | 0..=6      | *[Generic]*        |
/// | rrr   | bond topology        | 0..=2      | *[Query]*          |
/// | ccc   | bond reacting center | 0..=3      | *[Reaction,Query]* |
/// ------------------------------------------------------------------
///
fn bond_line<'a>() -> impl Parser<&'a [u8], Output = BondLine, Error = error::Error<&'a [u8]>> {
    let first_atom = fixed_width_int_minus1::<usize>(3);
    let second_atom = fixed_width_int_minus1::<usize>(3);
    let bond_type = fixed_width_int::<u8>(3);
    let bond_stereo = fixed_width_int::<u8>(3);
    all_consuming(map(
        (
            first_atom,
            second_atom,
            bond_type,
            bond_stereo,
            opt(preceded(
                take(3usize),
                (
                    complete(fixed_width_int::<u8>(3)),
                    opt(complete(fixed_width_int::<u8>(3))),
                ),
            )),
            space0,
        ),
        |(
            first_atom,  // '111' field: first atom
            second_atom, // '222' field: second atom
            bond_type,   // 'ttt' field: bond type
            bond_stereo, // 'sss' field: bond stereo
            rest,        // Optional fields:
            // 'rrr' field: bond topology
            // 'ccc' field: bond reacting center
            _,
        )| {
            let (bond_topology, bond_reacting_center) = match rest {
                Some((top, Some(rct))) => (top, rct),
                Some((top, None)) => (top, 0),
                None => (0, 0),
            };

            BondLine {
                first_atom,
                second_atom,
                bond_type,
                bond_stereo,
                bond_topology,
                bond_reacting_center,
            }
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::{error::ErrorKind, Err};
    use rstest::rstest;

    #[rstest]
    // From CTab spec (Figure 3)
    #[case(b"  1  2  1  0  0  0", BondLine { first_atom: 0, second_atom: 1, bond_type: 1, bond_stereo: 0, bond_topology: 0, bond_reacting_center: 0 })]
    #[case(b"  1  3  1  1  0  0", BondLine { first_atom: 0, second_atom: 2, bond_type: 1, bond_stereo: 1, bond_topology: 0, bond_reacting_center: 0 })]
    #[case(b"  1  4  1  0  0  0", BondLine { first_atom: 0, second_atom: 3, bond_type: 1, bond_stereo: 0, bond_topology: 0, bond_reacting_center: 0 })]
    #[case(b"  2  5  2  0  0  0", BondLine { first_atom: 1, second_atom: 4, bond_type: 2, bond_stereo: 0, bond_topology: 0, bond_reacting_center: 0 })]
    #[case(b"  2  6  1  0  0  0", BondLine { first_atom: 1, second_atom: 5, bond_type: 1, bond_stereo: 0, bond_topology: 0, bond_reacting_center: 0 })]
    // From RDKit test files
    #[case(b"  2  3  2  0  0  0  0", BondLine { first_atom: 1, second_atom: 2, bond_type: 2, bond_stereo: 0, bond_topology: 0, bond_reacting_center: 0 })]
    #[case(b"  3  5  1  0  0  2  0", BondLine { first_atom: 2, second_atom: 4, bond_type: 1, bond_stereo: 0, bond_topology: 2, bond_reacting_center: 0 })]
    #[case(b"  2  4  1  6  0  0  0", BondLine { first_atom: 1, second_atom: 3, bond_type: 1, bond_stereo: 6, bond_topology: 0, bond_reacting_center: 0 })]
    #[case(b"  2  3  1  0      ", BondLine { first_atom: 1, second_atom: 2, bond_type: 1, bond_stereo: 0, bond_topology: 0, bond_reacting_center: 0 })]
    fn test_bond_line(#[case] input: &[u8], #[case] expected: BondLine) {
        let (remaining, bond_line) = bond_line().parse(input).unwrap();
        assert_eq!(bond_line, expected);
        assert!(remaining.is_empty());
    }

    #[rstest]
    #[case(b"  1  2  1  ", "too few fields", ErrorKind::TakeWhileMN)]
    #[case(b"  1  2  1  0  0  0  0  0", "too many fields", ErrorKind::Eof)]
    fn test_bond_line_invalid(#[case] input: &[u8], #[case] desc: &str, #[case] expected_kind: ErrorKind) {
        let res = bond_line().parse(input);
        assert!(res.is_err(), "{} should have failed", desc);
        assert!(
            matches!(res.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            res.clone().unwrap_err().map(|e| e.code),
        );
    }
}
