//! Parsers for legacy atom list entries

use bstr::ByteSlice;
use nom::bytes::complete::{tag, take};
use nom::character::complete::space0;
use nom::combinator::{map, map_opt, map_res};
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::multi::length_count;
use nom::sequence::{delimited, preceded, terminated};
use nom::{Err, Parser};
use umol_data::Element;

use super::utils::{fixed_width_int_in_range, fixed_width_int_minus1, fixed_width_int_partial};
use crate::io::ctab::config::CtabParseFlags;
use crate::io::ctab::parser::properties::{AtomListEntry, PropertyEntries};

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
        let input = input.trim_end_with(|c| c == '\r' || c == '\n');
        if !legacy_atom_lists {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
        }
        map(
            (
                fixed_width_int_minus1::<usize>(3),
                delimited(
                    tag(" "),
                    map_res(take(1usize), |b: &[u8]| match b {
                        b"T" => Ok(true),
                        b"F" | b" " => Ok(false),
                        _ => Err(NomError::new(b, NomErrorKind::Tag)),
                    }),
                    tag("    "),
                ),
                terminated(
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
                    space0,
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
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::Element;

    use super::*;

    #[rstest]
    #[case::exclusion_flag_false(b"  1 F    3   9   7   8  ",
    PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: false, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::exclusion_flag_true(b"  1 T    3   9   7   8  ",
       PropertyEntries::AtomListEntry(AtomListEntry { atom_index: 0, exclusion: true, elements: vec![Element::F, Element::N, Element::O] }))]
    #[case::empty_exclusion_flag(b"  1      3   9   7   8  ",
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

    #[rstest]
    #[case::exclusion_flag_false(b"  1 F    3   9   7   8  ", NomErrorKind::Verify)]
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
