//! Counts line parser for CTab files.

use bstr::ByteSlice;
use nom::branch::alt;
use nom::bytes::complete::take;
use nom::character::complete::space0;
use nom::combinator::{complete, map, verify};
use nom::sequence::{delimited, preceded, terminated};
use nom::{error, Parser};

use super::utils::{fixed_width_int, fixed_width_padding_n};
use crate::io::config::ParseFlags;

/// Parse counts line (39 characters wide)
/// aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
///
/// *Values in the counts block*
/// ---------------------------------------------------------------
/// | Field   | Meaning                    | Values     | Notes   |
/// |---------|----------------------------|------------|---------|
/// | aaa     | number of atoms            | >0         | Generic |
/// | bbb     | number of bonds            | >=0        | Generic |
/// | lll     | number of atom lists       | 0..=30     | Generic |
/// | ccc     | chiral flag                | 0, 1       | Generic |
/// | sss     | number of stext entries    | >=0        | Generic |
/// | mmm     | number of properties lines | >=0        | Generic |
/// | vvvvvvv | version stamp              | V2000      | Generic |
/// ---------------------------------------------------------------
///

pub fn counts_input<'a>(
    flags: ParseFlags,
) -> impl Parser<&'a [u8], Output = Counts, Error = error::Error<&'a [u8]>> {
    terminated(counts_input_inner(flags), space0)
}

/// Internal parser for counts_input
fn counts_input_inner<'a>(
    flags: ParseFlags,
) -> impl Parser<&'a [u8], Output = Counts, Error = error::Error<&'a [u8]>> + 'a {
    let allow_unicode = flags.contains(ParseFlags::UNICODE);
    let strict_padding = flags.contains(ParseFlags::STRICT_PADDING);
    map(
        (
            fixed_width_int::<i32>(3usize, allow_unicode),
            fixed_width_int::<i32>(3usize, allow_unicode),
            fixed_width_int::<i32>(3usize, allow_unicode),
            preceded(take(3usize), fixed_width_int::<i32>(3usize, allow_unicode)),
            fixed_width_int::<i32>(3usize, allow_unicode),
            alt((
                delimited(
                    fixed_width_padding_n(4, 3, allow_unicode, strict_padding),
                    fixed_width_int::<i32>(3usize, allow_unicode),
                    complete(verify(take(6usize), |s: &[u8]| s.find(b"V2000").is_some())),
                ),
                delimited(
                    take(42usize),
                    fixed_width_int::<i32>(3usize, allow_unicode),
                    complete(verify(take(6usize), |s: &[u8]| s.find(b"V2000").is_some())),
                ),
            )),
        ),
        |(atoms, bonds, atom_lists, chiral_flag, stext_entries, properties_lines)| Counts {
            atoms,
            bonds,
            atom_lists,
            chiral_flag,
            stext_entries,
            properties_lines,
        },
    )
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Counts {
    pub atoms: i32,      // 'aaa' - number of atoms (max 255)
    pub bonds: i32,      // 'bbb' - number of bonds (max 255)
    pub atom_lists: i32, // 'lll' - number of atom lists (max 30)
    // fff is obsolete, skipping
    pub chiral_flag: i32,   // 'ccc' - chiral flag (0=not chiral, 1=chiral)
    pub stext_entries: i32, // 'sss' - number of stext entries
    // xxx is obsolete, skipping
    // rrr is obsolete, skipping
    // ppp is obsolete, skipping
    // iii is obsolete, skipping
    pub properties_lines: i32, // 'mmm' - number of additional properties lines
                               // 'vvvvvv' - version stamp (V2000), fixed string, skipping
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::combinator::all_consuming;
    use nom::{error, Err};
    use rstest::*;

    #[rstest]
    #[case(b"  6  5  0  0  1                 3 V2000", "counts", false,
      Counts {atoms: 6, bonds: 5, atom_lists: 0, chiral_flag: 1, stext_entries: 0, properties_lines: 3})]
    #[case(b"  1  0  0  0  0  0  0  0  0  0999 V2000", "zeroes, 999 properties", false,
      Counts {atoms: 1, bonds: 0, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 999})]
    //       aaabbblllfffcccsss            mmmvvvvvv
    #[case("  6  5  0  0  1\u{00A0}               3 V2000".as_bytes(), "counts", true,
      Counts {atoms: 6, bonds: 5, atom_lists: 0, chiral_flag: 1, stext_entries: 0, properties_lines: 3})]
    #[case(b"  1  0  0  0  0  0  0  0  0  0000 V2000    ", "padded", false,
      Counts {atoms: 1, bonds: 0, atom_lists: 0, chiral_flag: 0, stext_entries: 0, properties_lines: 0})]
    //       aaabbblllfffcccsss                                          mmmvvvvvv
    #[case(b"  4  4  1  0  0  0                                          999 V2000", "extra spaces", false,
      Counts {atoms: 4, bonds: 4, atom_lists: 1, chiral_flag: 0, stext_entries: 0, properties_lines: 999})]
    fn test_counts_input(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] allow_unicode: bool,
        #[case] expected: Counts,
    ) {
        let flags = if allow_unicode { ParseFlags::UNICODE } else { ParseFlags::empty() };
        let res = all_consuming(counts_input(flags)).parse(input);
        assert!(res.is_ok(), "{} should have succeeded", desc);
        let (remaining, counts) = res.unwrap();
        assert!(
            remaining.is_empty(),
            "{} should have consumed all input",
            desc
        );
        assert_eq!(counts, expected, "{} should have parsed correctly", desc);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case(b"  4  2  0     0  0            999 V1000", "invalid version", error::ErrorKind::Eof)]
    #[case(b"  4  2  0     0  0            ", "too short", error::ErrorKind::Eof)]
    #[case(b" 1A  2  0     0  0            999 V2000", "non-numeric atom", error::ErrorKind::Eof)]
    #[case(b"  4 AA  0     0  0            999 V2000", "non-numeric bond", error::ErrorKind::Digit)]
    #[case("  6  5  0  0  1\u{00A0}               3 V2000".as_bytes(), "unicode whitespace", error::ErrorKind::Digit)]
    fn test_counts_input_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let res = counts_input(ParseFlags::empty()).parse(input);
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
