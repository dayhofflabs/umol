//! Parsers for legacy atom list entries for CTab files.

use nom::bytes::complete::{tag, take};
use nom::character::complete::space0;
use nom::combinator::{all_consuming, map, map_opt, map_res};
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::multi::length_count;
use nom::sequence::{delimited, preceded, terminated};
use nom::{Err, Parser};
use umol_shared::element::Element;

use super::utils::{
    fixed_width_int_in_range, fixed_width_int_minus1, fixed_width_int_partial, LinesWithOffsetExt,
};
use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::error::ParseError;
use crate::io::ctfile::parser::properties::{AtomListEntry, PropertyEntries};

// Parse legacy atom list block
pub(super) fn legacy_atom_list_block<'inp>(
    atom_list_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = (Vec<PropertyEntries>, u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let mut properties = Vec::with_capacity(atom_list_count as usize);
        let mut lines_iter = input.lines_with_offset();
        let mut byte_offset = 0;

        for line_index in 0..atom_list_count {
            let (line, byte_len) = lines_iter.next().ok_or_else(|| {
                Err::Error(ParseError::UnexpectedEof {
                    line: line_offset + line_index,
                    block: "legacy atom list",
                })
            })?;

            let (_, property) = all_consuming(terminated(legacy_atom_list_input(flags), space0))
                .parse(line)
                .map_err(|e| {
                    Err::Error(ParseError::legacy_atom_list_from_nom(
                        e,
                        line_offset + line_index,
                        line,
                    ))
                })?;
            properties.push(property);
            byte_offset += byte_len;
        }

        let remaining = &input[byte_offset..];
        Ok((remaining, (properties, line_offset + atom_list_count)))
    }
}

/// Parse a legacy atom list entry
/// aaa k    n 111 222 333 444 555
/// aaa: atom index
/// k: exclusion flag
/// n: count
/// 111 222 333 444 555: element symbols
pub fn legacy_atom_list_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<&'inp [u8], Output = PropertyEntries, Error = NomError<&'inp [u8]>> + use<'inp> {
    let legacy_atom_lists = flags.contains(CtabParseFlags::LEGACY_ATOM_LISTS);
    move |input: &'inp [u8]| {
        if !legacy_atom_lists {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
        }
        map(
            (
                fixed_width_int_minus1::<u32>(3),
                delimited(
                    tag(" "),
                    map_res(take(1usize), |b: &[u8]| match b {
                        b"T" => Ok(true),
                        b"F" | b" " => Ok(false),
                        _ => Err(NomError::new(b, NomErrorKind::Tag)),
                    }),
                    tag("    "),
                ),
                length_count(
                    fixed_width_int_in_range::<u8, _>(1, 1..=5),
                    preceded(
                        tag(" "),
                        map_opt(
                            fixed_width_int_partial::<u8>(3),
                            Element::from_atomic_number,
                        ),
                    ),
                ),
            ),
            |(atom_index, exclusion, elements)| {
                PropertyEntries::AtomListEntry(AtomListEntry {
                    atom_index,
                    exclusion,
                    elements,
                })
            },
        )
        .parse(input)
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;

    #[rstest]
    #[case::legacy_feature(b"  1 F    3   9   7   8")]
    fn test_legacy_atom_list_block_extended_invalid(#[case] input: &[u8]) {
        let result = legacy_atom_list_block(1, 0, CtabParseFlags::EXTENDED).parse(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_err(), "{:?} should have failed", input_str);
        assert!(
            matches!(
                result,
                Err(Err::Error(ParseError::InvalidLegacyAtomListLine { .. }))
            ),
            "Expected InvalidLegacyAtomListLine for {:?}",
            input_str
        );
    }

    #[rstest]
    #[case::newline(b"  1 F    3   9   7   8\n",
        PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::crlf(b"  1 F    3   9   7   8\r\n",
        PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::trailing_space_lf(b"  1 F    3   9   7   8  \n",
        PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::trailing_space_crlf(b"  1 F    3   9   7   8  \r\n",
        PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::no_terminator(b"  1 F    3   9   7   8",
        PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::no_terminator_trailing_space(b"  1 F    3   9   7   8  ",
        PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
    fn test_legacy_atom_list_block_lenient(
        #[case] input: &[u8],
        #[case] expected: PropertyEntries,
    ) {
        let result = legacy_atom_list_block(1, 0, CtabParseFlags::LENIENT).parse(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_ok(), "{:?} should have succeeded", input_str);

        let (remaining, (atom_list, _)) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");

        assert_eq!(atom_list.len(), 1, "Should have 1 atom list");
        assert_eq!(atom_list[0], expected);
    }

    #[rstest]
    #[case::count_is_zero(b"  1 F    0   9  ")]
    fn test_legacy_atom_list_block_lenient_invalid(#[case] input: &[u8]) {
        let result = legacy_atom_list_block(1, 0, CtabParseFlags::LENIENT).parse(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_err(), "{:?} should have failed", input_str);
        assert!(
            matches!(
                result,
                Err(Err::Error(ParseError::InvalidLegacyAtomListLine { .. }))
            ),
            "Expected InvalidLegacyAtomListLine for {:?}",
            input_str
        );
    }

    #[rstest]
    #[case::legacy_feature(b"  1 F    3   9   7   8", NomErrorKind::Verify)]
    fn test_legacy_atom_list_input_extended_invalid(
        #[case] input: &[u8],
        #[case] expected_kind: NomErrorKind,
    ) {
        let result = legacy_atom_list_input(CtabParseFlags::EXTENDED).parse(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_err(), "{:?} should have failed", input_str);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "{:?} should have matched error for {:?}, got {:?}",
            expected_kind,
            input_str,
            result.clone().unwrap_err().map(|e| e.code)
        );
    }

    #[rstest]
    #[case::exclusion_flag_false(b"  1 F    3   9   7   8",
    PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::exclusion_flag_true(b"  1 T    3   9   7   8",
       PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: true, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::empty_exclusion_flag(b"  1      3   9   7   8",
       PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::partial_field(b"  4 F    4   6   7   8  16",
       PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 3, exclusion: false, elements: vec![Element::C, Element::N, Element::O, Element::S] }))]
    fn test_legacy_atom_list_input_lenient(
        #[case] input: &[u8],
        #[case] expected: PropertyEntries,
    ) {
        let result = legacy_atom_list_input(CtabParseFlags::LENIENT).parse(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_ok(), "{:?} should have succeeded", input_str);
        let (remaining, result) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::count_is_zero(b"  1 F    0   9  ", NomErrorKind::Verify)]
    #[case::count_exceeds_5(b"  1 F    6    9  ", NomErrorKind::Verify)]
    #[case::invalid_exclusion_flag(b"  1 X    1   9  ", NomErrorKind::MapRes)]
    #[case::invalid_element_atomic_number(b"  1 F    1   0  ", NomErrorKind::MapOpt)]
    fn test_legacy_atom_list_input_lenient_invalid(
        #[case] input: &[u8],
        #[case] expected_kind: NomErrorKind,
    ) {
        let result = legacy_atom_list_input(CtabParseFlags::LENIENT).parse(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_err(), "{:?} should have failed", input_str);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "{:?} should have matched error for {:?}, got {:?}",
            expected_kind,
            input_str,
            result.clone().unwrap_err().map(|e| e.code)
        );
    }
}
