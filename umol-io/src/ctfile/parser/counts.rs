//! Counts line parser for CTab files.

use winnow::error::ErrMode;
use winnow::token::take;
use winnow::{ModalResult, Parser};

use super::properties::{MoleculeChiralFlagEntry, PropertyEntries};
use super::utils::{
    finish_line, input_error_column, next_line, parse_int_opt, validate_unused_n, Input, InputError,
};
use crate::ctfile::config::CtabParseFlags;
use crate::ctfile::error::ParseError;

/// Parse counts block
pub(super) fn counts_block(
    input: &mut &[u8],
    line_offset: u32,
    flags: CtabParseFlags,
) -> ModalResult<(Counts, Vec<PropertyEntries>, u32), ParseError> {
    let mut line = next_line(input).map_err(|_| {
        ErrMode::Cut(ParseError::UnexpectedEof {
            line: line_offset,
            block: "counts",
        })
    })?;
    let result = counts_input(flags).parse_next(&mut line).and_then(|value| {
        finish_line(&mut line)?;
        Ok(value)
    });
    let (counts, properties) = result.map_err(|error| {
        let col = input_error_column(error, &line);
        ErrMode::Cut(ParseError::InvalidCountsLine {
            line: line_offset,
            col,
        })
    })?;

    Ok((counts, properties, line_offset + 1))
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
fn counts_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<Input<'inp>, (Counts, Vec<PropertyEntries>), ErrMode<InputError>> + use<'inp> {
    let skip_unused_fields = flags.contains(CtabParseFlags::SKIP_UNUSED_FIELDS);
    let no_v2000_end_tags = flags.contains(CtabParseFlags::NO_V2000_END_TAGS);

    move |input: &mut Input<'inp>| {
        let bytes: &[u8] = input.as_ref();
        let offset;

        if bytes.len() < 6 {
            return Err(ErrMode::Backtrack(InputError {
                column: bytes.len() as u32,
            }));
        }

        // aaa: atom count (0-2)
        let atom_count: u32 = parse_int_opt::<u32>(&bytes[0..3], 0)?.unwrap_or_default();

        // bbb: bond count (3-5)
        let bond_count: u32 = parse_int_opt::<u32>(&bytes[3..6], 3)?.unwrap_or_default();

        // lll: atom list count (6-8)
        let atom_list_count: u32 = if bytes.len() >= 9 {
            parse_int_opt::<u32>(&bytes[6..9], 6)?.unwrap_or_default()
        } else {
            0
        };

        // fff: obsolete field (9-11)
        if bytes.len() >= 12 {
            validate_unused_n(&bytes[9..12], 1, 3, skip_unused_fields, 9)?;
        }

        // ccc: chiral flag (12-14)
        let chiral_flag = if bytes.len() >= 15 {
            let val = parse_int_opt::<u8>(&bytes[12..15], 12)?.unwrap_or_default();
            if val > 1 {
                return Err(ErrMode::Backtrack(InputError { column: 12 }));
            }
            val
        } else {
            0
        };

        // sss, xxx, rrr, ppp, iii: obsolete fields (15-29)
        let count = (bytes.len().saturating_sub(15) / 3).min(5);
        if count > 0 {
            validate_unused_n(&bytes[15..15 + count * 3], count, 3, skip_unused_fields, 15)?;
        }

        // mmm: properties line count (30-32) - parsed as integer, value ignored
        if bytes.len() >= 33 {
            let _ = parse_int_opt::<u32>(&bytes[30..33], 30)?;
        }

        // vvvvvv: version stamp (33-38)
        if no_v2000_end_tags {
            offset = bytes.len();
        } else {
            if bytes.len() < 39 {
                return Err(ErrMode::Backtrack(InputError {
                    column: bytes.len() as u32,
                }));
            }
            offset = 39;
            if !(bytes[33..39] == *b" V2000" || bytes[33..39] == *b"V2000 ") {
                return Err(ErrMode::Backtrack(InputError { column: 33 }));
            }
        }

        let properties = if chiral_flag == 1 {
            vec![PropertyEntries::MoleculeChiralFlagEntry(
                MoleculeChiralFlagEntry { chiral_flag: true },
            )]
        } else {
            vec![]
        };

        let _: &[u8] = take(offset).parse_next(input)?;
        Ok((
            Counts {
                atom_count,
                bond_count,
                atom_list_count,
            },
            properties,
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
    use pretty_assertions::assert_eq;
    use rstest::*;
    use winnow::error::ErrMode;
    use winnow::Parser;

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
        let mut remaining = input;
        assert_eq!(
            counts_block(&mut remaining, 0, CtabParseFlags::BASIC),
            Ok((expected_counts, expected_properties, 1))
        );
        assert!(remaining.is_empty());
    }

    #[rstest]
    #[case::no_v2000_tag(b" 28 34                           ", 33)]
    #[case::trailing_characters(b"  2  1  0  0  0  0  0  0  0  0  0    0", 38)]
    #[case::rxn_header(b"$RXN\n", 4)]
    fn test_counts_block_error(#[case] input: &[u8], #[case] col: u32) {
        let mut remaining = input;
        assert_eq!(
            counts_block(&mut remaining, 0, CtabParseFlags::BASIC),
            Err(ErrMode::Cut(ParseError::InvalidCountsLine { line: 0, col }))
        );
    }

    #[rstest]
    #[case::non_zero_unused(b"  4  2  0     0  1                V2000", 15)]
    fn test_counts_block_strict_error(#[case] input: &[u8], #[case] col: u32) {
        let mut remaining = input;
        assert_eq!(
            counts_block(&mut remaining, 0, CtabParseFlags::STRICT),
            Err(ErrMode::Cut(ParseError::InvalidCountsLine { line: 0, col }))
        );
    }

    #[rstest]
    #[case::empty(b"")]
    fn test_counts_block_eof_error(#[case] input: &[u8]) {
        let mut remaining = input;
        assert_eq!(
            counts_block(&mut remaining, 7, CtabParseFlags::BASIC),
            Err(ErrMode::Cut(ParseError::UnexpectedEof {
                line: 7,
                block: "counts",
            }))
        );
    }

    #[rstest]
    #[case::no_v2000_tag(b" 28 34                           ", Counts { atom_count: 28, bond_count: 34, atom_list_count: 0 }, vec![])]
    #[case::trailing_characters(b"  2  1  0  0  0  0  0  0  0  0  0    0", Counts { atom_count: 2, bond_count: 1, atom_list_count: 0 }, vec![])]
    fn test_counts_block_lenient(
        #[case] input: &[u8],
        #[case] expected_counts: Counts,
        #[case] expected_properties: Vec<PropertyEntries>,
    ) {
        let mut remaining = input;
        assert_eq!(
            counts_block(
                &mut remaining,
                0,
                CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT,
            ),
            Ok((expected_counts, expected_properties, 1))
        );
        assert!(remaining.is_empty());
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
        assert_eq!(
            counts_input(CtabParseFlags::BASIC).parse(Input::new(input)),
            Ok((expected_counts, expected_properties))
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::blank(b"                                       ", 33)]
    #[case::invalid_version(b"  4  2  0     0                   V1000", 33)]
    #[case::len_9_too_short(b"  4  2  0", 9)]
    #[case::len_32_too_short(b"  4  2  0     0                 ", 32)]
    #[case::len_38_malformed(b" 6  6  0  0  0  0  0  0  0  0  1 V2000", 38)]
    #[case::negative_atom_count(b" -1  2  0     0                   V2000", 0)]
    #[case::non_numeric_atom_count(b"  a  2  0     0                   V2000", 0)]
    #[case::trailing_chars_atom_count(b" 1a  2  0     0                   V2000", 0)]
    #[case::negative_bond_count(b"  4 -2  0     0                   V2000", 3)]
    #[case::non_numeric_bond_count(b"  4  a  0     0                   V2000", 3)]
    #[case::non_numeric_atom_list_count(b"  4  2  a     0                   V2000", 6)]
    #[case::chiral_flag_out_of_range(b"  1  0  0  0  2  0  0  0  0  0  0 V2000", 12)]
    #[case::non_numeric_chiral_flag(b"  4  2  0     a                   V2000", 12)]
    #[case::non_numeric_properties(b"  4  2  0     0  1              a V2000", 30)]
    fn test_counts_input_error(#[case] input: &[u8], #[case] column: u32) {
        let error = counts_input(CtabParseFlags::BASIC)
            .parse(Input::new(input))
            .unwrap_err()
            .into_inner();
        assert_eq!(error, InputError { column });
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::non_zero_unused(b"  4  2  0     0  1                V2000", 15)]
    #[case::non_numeric_unused(b"  4  2  0     0  a                V2000", 15)]
    fn test_counts_input_strict_error(#[case] input: &[u8], #[case] column: u32) {
        let error = counts_input(CtabParseFlags::STRICT)
            .parse(Input::new(input))
            .unwrap_err()
            .into_inner();
        assert_eq!(error, InputError { column });
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
        assert_eq!(
            counts_input(CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT)
                .parse(Input::new(input)),
            Ok((expected_counts, expected_properties))
        );
    }

    #[rstest]
    #[case::empty(b"", 0)]
    fn test_counts_input_lenient_error(#[case] input: &[u8], #[case] column: u32) {
        let error = counts_input(CtabParseFlags::BASIC_MAX & CtabParseFlags::LENIENT)
            .parse(Input::new(input))
            .unwrap_err()
            .into_inner();
        assert_eq!(error, InputError { column });
    }
}
