//! Counts line parser for CTab files.

use nom::character::complete::space0;
use nom::combinator::all_consuming;
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::sequence::terminated;
use nom::{Err, Parser};

use super::properties::{MoleculeChiralFlagEntry, PropertyEntries};
use super::utils::{parse_int_opt, validate_unused_n, LinesWithOffsetExt};
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;

/// Parse counts block
pub(super) fn counts_block<'inp>(
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Counts, Vec<PropertyEntries>, u32), Error = ParseError> + use<'inp>
{
    move |input: &'inp [u8]| {
        let (line, byte_len) = input.lines_with_offset().next().ok_or_else(|| {
            Err::Error(ParseError::UnexpectedEof {
                line: line_offset,
                block: "counts",
            })
        })?;

        let (_, (counts, properties)) = all_consuming(terminated(counts_input(flags), space0))
            .parse(line)
            .map_err(|e| Err::Error(ParseError::counts_from_nom(e, line_offset)))?;

        let remaining = &input[byte_len..];
        Ok((remaining, (counts, properties, line_offset + 1)))
    }
}

/// Parse counts line (39 characters wide)
/// aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
///
/// *Values in the counts block*
/// --------------------------------------------------------------------------
/// | Field   | Position | Meaning                    | Values     | Notes   |
/// |---------|----------|----------------------------|------------|---------|
/// | aaa     |  1- 3    | number of atoms            | >0         | Generic |
/// | bbb     |  4- 6    | number of bonds            | >=0        | Generic |
/// | lll     |  7- 9    | number of atom lists       | 0..=30     | Generic |
/// | fff     | 10-12    | obsolete                   | -          | Generic |
/// | ccc     | 13-15    | chiral flag                | 0, 1       | Generic |
/// | sss     | 16-18    | stext entries (unused)     | >=0        | Generic |
/// | xxx     | 19-21    | obsolete                   | -          | Generic |
/// | rrr     | 22-24    | obsolete                   | -          | Generic |
/// | ppp     | 25-27    | obsolete                   | -          | Generic |
/// | iii     | 28-30    | obsolete                   | -          | Generic |
/// | mmm     | 31-33    | properties lines (unused)  | >=0        | Generic |
/// | vvvvvv  | 34-39    | version stamp              | V2000      | Generic |
/// --------------------------------------------------------------------------
///
pub fn counts_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Counts, Vec<PropertyEntries>), Error = NomError<&'inp [u8]>>
       + use<'inp> {
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let no_v2000_end_tags = flags.contains(CtabParseFlags::NO_V2000_END_TAGS);

    move |input: &'inp [u8]| {
        let offset;

        if input.len() < 6 {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        // aaa: atom count (0-2)
        let atom_count: u32 = parse_int_opt::<u32>(input, &input[0..3])?.unwrap_or_default();

        // bbb: bond count (3-5)
        let bond_count: u32 = parse_int_opt::<u32>(input, &input[3..6])?.unwrap_or_default();

        // lll: atom list count (6-8)
        let atom_list_count: u32 = if input.len() >= 9 {
            parse_int_opt::<u32>(input, &input[6..9])?.unwrap_or_default()
        } else {
            0
        };

        // fff: obsolete field (9-11)
        if input.len() >= 12 {
            validate_unused_n(input, &input[9..12], 1, 3, skip_unused_fields)?;
        }

        // ccc: chiral flag (12-14)
        let chiral_flag = if input.len() >= 15 {
            let val = parse_int_opt::<u8>(input, &input[12..15])?.unwrap_or_default();
            if val > 1 {
                return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
            }
            val
        } else {
            0
        };

        // sss, xxx, rrr, ppp, iii: obsolete fields (15-29)
        let count = (input.len().saturating_sub(15) / 3).min(5);
        if count > 0 {
            validate_unused_n(
                input,
                &input[15..15 + count * 3],
                count,
                3,
                skip_unused_fields,
            )?;
        }

        // mmm: properties line count (30-32) - parsed as integer, value ignored
        if input.len() >= 33 {
            let _ = parse_int_opt::<u32>(input, &input[30..33])?;
        }

        // vvvvvv: version stamp (33-38)
        if no_v2000_end_tags {
            offset = input.len();
        } else {
            if input.len() < 39 {
                return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
            }
            offset = 39;
            if !(input[33..39] == *b" V2000" || input[33..39] == *b"V2000 ") {
                return Err(Err::Error(NomError::new(input, NomErrorKind::Tag)));
            }
        }

        let properties = if chiral_flag == 1 {
            vec![PropertyEntries::MoleculeChiralFlagEntry(
                MoleculeChiralFlagEntry { chiral_flag: true },
            )]
        } else {
            vec![]
        };

        Ok((
            &input[offset..],
            (
                Counts {
                    atom_count,
                    bond_count,
                    atom_list_count,
                },
                properties,
            ),
        ))
    }
}

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Counts {
    pub atom_count: u32,
    pub bond_count: u32,
    pub atom_list_count: u32,
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;
    use nom::error::ErrorKind as NomErrorKind;
    use nom::Err;
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::zeroes(b"  0  0  0  0  0  0  0  0  0  0  0 V2000", Counts {atom_count: 0, bond_count: 0, atom_list_count: 0}, vec![])]
    #[case::padded_newline(b"  1  0  0  0  0  0  0  0  0  0  0 V2000    \n", Counts {atom_count: 1, bond_count: 0, atom_list_count: 0}, vec![])]
    #[case::padded_crlf(b"  6  5  1     1                   V2000\r\n", Counts {atom_count: 6, bond_count: 5, atom_list_count: 1},
        vec![PropertyEntries::MoleculeChiralFlagEntry(MoleculeChiralFlagEntry { chiral_flag: true })])]
    #[case::properties_999(b"  6  5  1     1               999 V2000\n", Counts {atom_count: 6, bond_count: 5, atom_list_count: 1},
        vec![PropertyEntries::MoleculeChiralFlagEntry(MoleculeChiralFlagEntry { chiral_flag: true })])]
    #[case::no_terminator(b"  1  0  0  0  0  0  0  0  0  0  0 V2000", Counts {atom_count: 1, bond_count: 0, atom_list_count: 0}, vec![])]
    #[case::tag_only(b"                                  V2000\n", Counts {atom_count: 0, bond_count: 0, atom_list_count: 0}, vec![])]
    #[case::chiral_flag(b"  1  0  0  0  1  0  0  0  0  0  0 V2000", Counts {atom_count: 1, bond_count: 0, atom_list_count: 0},
        vec![PropertyEntries::MoleculeChiralFlagEntry(MoleculeChiralFlagEntry { chiral_flag: true })])]
    #[case::non_zero_unused(b"  4  2  0     0  1                V2000", Counts {atom_count: 4, bond_count: 2, atom_list_count: 0}, vec![])]
    fn test_counts_block(
        #[case] input: &[u8],
        #[case] expected_counts: Counts,
        #[case] expected_properties: Vec<PropertyEntries>,
    ) {
        let res = counts_block(0, CtabParseFlags::BASIC).parse(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_ok(), "{:?} should have succeeded", input_str);
        let (remaining, (counts, properties, line_offset)) = res.unwrap();
        assert!(
            remaining.is_empty(),
            "{:?} should have consumed all input",
            input_str
        );
        assert_eq!(
            counts, expected_counts,
            "{:?} should have parsed counts correctly",
            input_str
        );
        assert_eq!(
            properties, expected_properties,
            "{:?} should have parsed properties correctly",
            input_str
        );
        assert_eq!(
            line_offset, 1,
            "{:?} should have incremented line offset",
            input_str
        );
    }

    #[rstest]
    #[case::no_v2000_tag(b" 28 34                           ")]
    #[case::trailing_characters(b"  2  1  0  0  0  0  0  0  0  0  0    0")]
    #[case::rxn_header(b"$RXN\n")]
    fn test_counts_block_invalid(#[case] input: &[u8]) {
        let res = counts_block(0, CtabParseFlags::BASIC).parse(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_err(), "{:?} should have failed", input_str);
        assert!(matches!(
            res.unwrap_err(),
            Err::Error(ParseError::InvalidCountsLine { .. })
        ));
    }

    #[rstest]
    #[case::non_zero_unused(b"  4  2  0     0  1                V2000")]
    fn test_counts_block_strict_invalid(#[case] input: &[u8]) {
        let res = counts_block(0, CtabParseFlags::STRICT).parse(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_err(), "{:?} should have failed", input_str);
        assert!(matches!(
            res.unwrap_err(),
            Err::Error(ParseError::InvalidCountsLine { .. })
        ));
    }

    #[rstest]
    #[case::no_v2000_tag(b" 28 34                           ", Counts { atom_count: 28, bond_count: 34, atom_list_count: 0 }, vec![])]
    #[case::trailing_characters(b"  2  1  0  0  0  0  0  0  0  0  0    0", Counts { atom_count: 2, bond_count: 1, atom_list_count: 0 }, vec![])]
    fn test_counts_block_lenient(
        #[case] input: &[u8],
        #[case] expected_counts: Counts,
        #[case] expected_properties: Vec<PropertyEntries>,
    ) {
        let res = counts_block(0, CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT).parse(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_ok(), "{:?} should have succeeded", input_str);
        let (remaining, (counts, properties, line_offset)) = res.unwrap();
        assert!(
            remaining.is_empty(),
            "{:?} should have consumed all input",
            input_str
        );
        assert_eq!(
            counts, expected_counts,
            "{:?} should have parsed counts correctly",
            input_str
        );
        assert_eq!(
            properties, expected_properties,
            "{:?} should have parsed properties correctly",
            input_str
        );
        assert_eq!(
            line_offset, 1,
            "{:?} should have incremented line offset",
            input_str
        );
    }

    #[rstest]
    #[case::zeroes(b"  0  0  0  0  0  0  0  0  0  0  0 V2000", Counts {atom_count: 0, bond_count: 0, atom_list_count: 0}, vec![])]
    #[case::atom_count(b"  1  0  0  0  0  0  0  0  0  0  0 V2000", Counts {atom_count: 1, bond_count: 0, atom_list_count: 0}, vec![])]
    #[case::properties_999(b"  6  5  1     1               999 V2000", Counts {atom_count: 6, bond_count: 5, atom_list_count: 1},
        vec![PropertyEntries::MoleculeChiralFlagEntry(MoleculeChiralFlagEntry { chiral_flag: true })])]
    #[case::tag_only(b"                                  V2000", Counts {atom_count: 0, bond_count: 0, atom_list_count: 0}, vec![])]
    #[case::padded_blanks(b"  6  5  1                         V2000", Counts {atom_count: 6, bond_count: 5, atom_list_count: 1}, vec![])]
    #[case::chiral_flag(b"  1  0  0  0  1  0  0  0  0  0  0 V2000", Counts {atom_count: 1, bond_count: 0, atom_list_count: 0},
        vec![PropertyEntries::MoleculeChiralFlagEntry(MoleculeChiralFlagEntry { chiral_flag: true })])]
    #[case::invalid_unused(b"  4  2  0     0  1                V2000", Counts {atom_count: 4, bond_count: 2, atom_list_count: 0}, vec![])]
    fn test_counts_input(
        #[case] input: &[u8],
        #[case] expected_counts: Counts,
        #[case] expected_properties: Vec<PropertyEntries>,
    ) {
        let res = counts_input(CtabParseFlags::BASIC).parse(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_ok(), "{:?} should have succeeded", input_str);
        let (remaining, (counts, properties)) = res.unwrap();
        assert!(
            remaining.is_empty(),
            "{:?} should have consumed all input",
            input_str
        );
        assert_eq!(
            counts, expected_counts,
            "{:?} should have parsed counts correctly",
            input_str
        );
        assert_eq!(
            properties, expected_properties,
            "{:?} should have parsed properties correctly",
            input_str
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::blank(b"                                       ", NomErrorKind::Tag)]
    #[case::invalid_version(b"  4  2  0     0                   V1000", NomErrorKind::Tag)]
    #[case::len_9_too_short(b"  4  2  0", NomErrorKind::Eof)]
    #[case::len_32_too_short(b"  4  2  0     0                 ", NomErrorKind::Eof)]
    #[case::len_38_malformed(b" 6  6  0  0  0  0  0  0  0  0  1 V2000", NomErrorKind::Eof)]
    #[case::negative_atom_count(b" -1  2  0     0                   V2000", NomErrorKind::Digit)]
    #[case::non_numeric_atom_count(b"  a  2  0     0                   V2000", NomErrorKind::Digit)]
    #[case::trailing_chars_atom_count(b" 1a  2  0     0                   V2000", NomErrorKind::Eof)]
    #[case::negative_bond_count(b"  4 -2  0     0                   V2000", NomErrorKind::Digit)]
    #[case::non_numeric_bond_count(b"  4  a  0     0                   V2000", NomErrorKind::Digit)]
    #[case::non_numeric_atom_list_count(b"  4  2  a     0                   V2000", NomErrorKind::Digit)]
    #[case::chiral_flag_out_of_range(b"  1  0  0  0  2  0  0  0  0  0  0 V2000", NomErrorKind::Verify)]    
    #[case::non_numeric_chiral_flag(b"  4  2  0     a                   V2000", NomErrorKind::Digit)]
    #[case::non_numeric_properties(b"  4  2  0     0  1              a V2000", NomErrorKind::Digit)]
    fn test_counts_input_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
        let res = counts_input(CtabParseFlags::BASIC).parse(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_err(), "{:?} should have failed", input_str);
        assert!(
            matches!(res.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "{:?} should have failed with error kind {:?}, got {:?}",
            input_str,
            expected_kind,
            res.clone().unwrap_err().map(|e| e.code)
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::non_zero_unused(b"  4  2  0     0  1                V2000", NomErrorKind::Verify)]
    #[case::non_numeric_unused(b"  4  2  0     0  a                V2000", NomErrorKind::Verify)]
    fn test_counts_input_strict_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
        let res = counts_input(CtabParseFlags::STRICT).parse(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_err(), "{:?} should have failed", input_str);
        assert!(
            matches!(res.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "{:?} should have failed with error kind {:?}, got {:?}",
            input_str,
            NomErrorKind::Verify,
            res.clone().unwrap_err().map(|e| e.code)
        );
    }

    #[rstest]
    #[case::blank(b"                                 ", Counts { atom_count: 0, bond_count: 0, atom_list_count: 0 }, vec![])]
    #[case::padded_blanks(b" 28 34                           ", Counts { atom_count: 28, bond_count: 34, atom_list_count: 0 }, vec![])]
    #[case::padded_zeros(b" 28 34  0  0  0  0  0  0  0  0  0", Counts { atom_count: 28, bond_count: 34, atom_list_count: 0 }, vec![])]
    #[case::has_v2000_tag(b"  0  0  0  0  0  0  0  0  0  0  0 V2000", Counts { atom_count: 0, bond_count: 0, atom_list_count: 0 }, vec![])]
    fn test_counts_input_lenient(
        #[case] input: &[u8],
        #[case] expected_counts: Counts,
        #[case] expected_properties: Vec<PropertyEntries>,
    ) {
        let res = counts_input(CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT).parse(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_ok(), "{:?} should have succeeded", input_str);
        let (remaining, (counts, properties)) = res.unwrap();
        assert!(
            remaining.is_empty(),
            "{:?} should have consumed all input, remaining: {:?}",
            input_str,
            remaining.to_str_lossy()
        );
        assert_eq!(
            counts, expected_counts,
            "{:?} should have parsed counts correctly",
            input_str
        );
        assert_eq!(
            properties, expected_properties,
            "{:?} should have parsed properties correctly",
            input_str
        );
    }

    #[rstest]
    #[case::empty(b"", NomErrorKind::Eof)]
    fn test_counts_input_lenient_invalid(
        #[case] input: &[u8],
        #[case] expected_kind: NomErrorKind,
    ) {
        let res = counts_input(CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT).parse(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_err(), "{:?} should have failed", input_str);
        assert!(
            matches!(res.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "{:?} should have failed with error kind {:?}, got {:?}",
            input_str,
            expected_kind,
            res.clone().unwrap_err().map(|e| e.code)
        );
    }
}
