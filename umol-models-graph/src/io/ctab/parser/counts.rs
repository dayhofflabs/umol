//! Counts line parser for CTab files.

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::space0;
use nom::combinator::{map, opt, value};
use nom::error::Error as NomError;
use nom::sequence::{delimited, preceded, terminated};
use nom::Parser;

use super::utils::{fixed_width_int, fixed_width_padding_n};
use crate::io::ctab::config::CtabParseFlags;

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
pub fn counts_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = Counts, Error = NomError<&'inp [u8]>> + use<'inp> {
    let skip_padding = flags.contains(CtabParseFlags::SKIP_PADDING);
    let no_v2000_end_tags = flags.contains(CtabParseFlags::NO_V2000_END_TAGS);
    terminated(
        map(
            (
                fixed_width_int::<u32>(3),
                fixed_width_int::<u32>(3),
                fixed_width_int::<u32>(3),
                preceded(take(3usize), fixed_width_int::<u32>(3)),
                fixed_width_int::<u32>(3),
                delimited(
                    fixed_width_padding_n(4, 3, skip_padding),
                    fixed_width_int::<u32>(3),
                    version(no_v2000_end_tags),
                ),
            ),
            |(
                atom_count,
                bond_count,
                atom_list_count,
                chiral_flag,
                stext_entry_count,
                properties_lines,
            )| Counts {
                atom_count,
                bond_count,
                atom_list_count,
                chiral_flag,
                stext_entry_count,
                properties_lines,
            },
        ),
        space0,
    )
}

/// Parse version stamp
fn version<'inp>(
    no_v2000_end_tags: bool,
) -> impl Parser<&'inp [u8], Output = (), Error = NomError<&'inp [u8]>> {
    move |input: &'inp [u8]| {
        let v2000 = alt((tag(" V2000"), tag("V2000 ")));
        if no_v2000_end_tags {
            value((), opt(v2000)).parse(input)
        } else {
            value((), v2000).parse(input)
        }
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Counts {
    pub atom_count: u32,      // 'aaa' - number of atoms (max 255)
    pub bond_count: u32,      // 'bbb' - number of bonds (max 255)
    pub atom_list_count: u32, // 'lll' - number of atom lists (max 30)
    // fff is obsolete, skipping
    pub chiral_flag: u32,       // 'ccc' - chiral flag (0=not chiral, 1=chiral)
    pub stext_entry_count: u32, // 'sss' - number of stext entries
    // xxx is obsolete, skipping
    // rrr is obsolete, skipping
    // ppp is obsolete, skipping
    // iii is obsolete, skipping
    pub properties_lines: u32, // 'mmm' - number of additional properties lines
                               // 'vvvvvv' - version stamp (V2000), fixed string, skipping
}

#[cfg(test)]
mod tests {
    use nom::combinator::all_consuming;
    use nom::error::ErrorKind as NomErrorKind;
    use nom::Err;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case(b"  6  5  0  0  1                 3 V2000", "counts",
      Counts {atom_count: 6, bond_count: 5, atom_list_count: 0, chiral_flag: 1, stext_entry_count: 0, properties_lines: 3})]
    #[case(b"  1  0  0  0  0  0  0  0  0  0999 V2000", "zeroes, 999 properties",
      Counts {atom_count: 1, bond_count: 0, atom_list_count: 0, chiral_flag: 0, stext_entry_count: 0, properties_lines: 999})]
    #[case(b"  1  0  0  0  0  0  0  0  0  0000 V2000    ", "padded",
      Counts {atom_count: 1, bond_count: 0, atom_list_count: 0, chiral_flag: 0, stext_entry_count: 0, properties_lines: 0})]
    fn test_counts_input(#[case] input: &[u8], #[case] desc: &str, #[case] expected: Counts) {
        let res = all_consuming(counts_input(CtabParseFlags::STRICT)).parse(input);
        assert!(res.is_ok(), "{} should have succeeded", desc);
        let (remaining, counts) = res.unwrap();
        assert!(
            remaining.is_empty(),
            "{} should have consumed all input",
            desc
        );
        assert_eq!(counts, expected, "{} should have parsed correctly", desc);
    }

    #[rstest]
    #[case(b" 28 34  0  0  0  0  0  0  0  0  0", "no V2000 tag",
      Counts {atom_count: 28, bond_count: 34, atom_list_count: 0, chiral_flag: 0, stext_entry_count: 0, properties_lines: 0})]
    #[case(b"                                                                     ", "len 69, blank",
      Counts {atom_count: 0, bond_count: 0, atom_list_count: 0, chiral_flag: 0, stext_entry_count: 0, properties_lines: 0})]
    fn test_counts_input_no_v2000_tag(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected: Counts,
    ) {
        let res = all_consuming(counts_input(CtabParseFlags::NO_V2000_END_TAGS)).parse(input);
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
    #[case(b"  4  2  0     0  0            999 V1000", "invalid version", NomErrorKind::Tag)]
    #[case(b"  4  2  0     0  0            ", "too short", NomErrorKind::Tag)]
    #[case(b" 1A  2  0     0  0            999 V2000", "non-numeric atom", NomErrorKind::Eof)]
    #[case(b"  4 AA  0     0  0            999 V2000", "non-numeric bond", NomErrorKind::Digit)]
    fn test_counts_input_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: NomErrorKind,
    ) {
        let res = counts_input(CtabParseFlags::STRICT).parse(input);
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
