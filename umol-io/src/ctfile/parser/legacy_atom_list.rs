//! Parsers for legacy atom list entries for CTab files.

use umol_chem::element::Element;
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::{ModalResult, Parser};

use super::utils::{finish_line, input_error_column, next_line, parse_int_opt, Input, InputError};
use crate::ctfile::config::CtabParseFlags;
use crate::ctfile::error::ParseError;
use crate::ctfile::parser::properties::{AtomListEntry, PropertyEntries};

// Parse legacy atom list block
pub(super) fn legacy_atom_list_block(
    input: &mut &[u8],
    atom_list_count: u32,
    line_offset: u32,
    flags: CtabParseFlags,
) -> ModalResult<(Vec<PropertyEntries>, u32), ParseError> {
    let mut properties = Vec::with_capacity(atom_list_count as usize);
    for line_index in 0..atom_list_count {
        let physical_line = line_offset + line_index;
        let mut line = next_line(input).map_err(|_| {
            ErrMode::Cut(ParseError::UnexpectedEof {
                line: physical_line,
                block: "legacy atom list",
            })
        })?;
        let result = legacy_atom_list_input(flags)
            .parse_next(&mut line)
            .and_then(|value| {
                finish_line(&mut line)?;
                Ok(value)
            });
        let property = result.map_err(|error| {
            ErrMode::Cut(ParseError::InvalidLegacyAtomListLine {
                line: physical_line,
                col: input_error_column(error, &line),
            })
        })?;
        properties.push(property);
    }
    Ok((properties, line_offset + atom_list_count))
}

/// Parse a legacy atom list entry
/// aaa k    n 111 222 333 444 555
/// aaa: atom index
/// k: exclusion flag
/// n: count
/// 111 222 333 444 555: element symbols
fn legacy_atom_list_input<'inp>(
    flags: CtabParseFlags,
) -> impl Parser<Input<'inp>, PropertyEntries, ErrMode<InputError>> + use<'inp> {
    let legacy_atom_lists = flags.contains(CtabParseFlags::LEGACY_ATOM_LISTS);
    move |input: &mut Input<'inp>| {
        if !legacy_atom_lists {
            return Err(ErrMode::Backtrack(InputError { column: 0 }));
        }
        let bytes: &[u8] = input.as_ref();
        if bytes.len() < 10 {
            return Err(ErrMode::Backtrack(InputError {
                column: bytes.len() as u32,
            }));
        }

        let atom_index = parse_int_opt::<u32>(&bytes[0..3], 0)?
            .and_then(|value| value.checked_sub(1))
            .ok_or(ErrMode::Backtrack(InputError { column: 0 }))?;
        if bytes[3] != b' ' {
            return Err(ErrMode::Backtrack(InputError { column: 3 }));
        }
        let exclusion = match bytes[4] {
            b'T' => true,
            b'F' | b' ' => false,
            _ => return Err(ErrMode::Backtrack(InputError { column: 4 })),
        };
        if bytes[5..9] != *b"    " {
            return Err(ErrMode::Backtrack(InputError { column: 5 }));
        }
        let element_count = parse_int_opt::<u8>(&bytes[9..10], 9)?
            .filter(|value| (1..=5).contains(value))
            .ok_or(ErrMode::Backtrack(InputError { column: 9 }))?;

        let mut offset = 10;
        let mut elements = Vec::with_capacity(element_count as usize);
        for _ in 0..element_count {
            if bytes.get(offset) != Some(&b' ') {
                return Err(ErrMode::Backtrack(InputError {
                    column: offset as u32,
                }));
            }
            offset += 1;
            let field_end = (offset + 3).min(bytes.len());
            let atomic_number = parse_int_opt::<u8>(&bytes[offset..field_end], offset as u32)?
                .ok_or(ErrMode::Backtrack(InputError {
                    column: offset as u32,
                }))?;
            let element = Element::from_atomic_number(atomic_number).ok_or(ErrMode::Backtrack(
                InputError {
                    column: offset as u32,
                },
            ))?;
            elements.push(element);
            offset = field_end;
        }

        let _: &[u8] = take(offset).parse_next(input)?;
        Ok(PropertyEntries::AtomListEntry(AtomListEntry {
            atom_index,
            exclusion,
            elements,
        }))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use winnow::error::ErrMode;
    use winnow::Parser;

    use super::*;

    #[rstest]
    #[case::legacy_feature(b"  1 F    3   9   7   8")]
    fn test_legacy_atom_list_block_extended_error(#[case] input: &[u8]) {
        let mut remaining = input;
        assert_eq!(
            legacy_atom_list_block(&mut remaining, 1, 0, CtabParseFlags::EXTENDED),
            Err(ErrMode::Cut(ParseError::InvalidLegacyAtomListLine {
                line: 0,
                col: 0,
            }))
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
        let mut remaining = input;
        let (atom_list, line_offset) =
            legacy_atom_list_block(&mut remaining, 1, 0, CtabParseFlags::LENIENT).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(line_offset, 1);
        assert_eq!(atom_list.len(), 1, "Should have 1 atom list");
        assert_eq!(atom_list[0], expected);
    }

    #[rstest]
    #[case::count_is_zero(b"  1 F    0   9  ")]
    fn test_legacy_atom_list_block_lenient_error(#[case] input: &[u8]) {
        let mut remaining = input;
        assert_eq!(
            legacy_atom_list_block(&mut remaining, 1, 0, CtabParseFlags::LENIENT),
            Err(ErrMode::Cut(ParseError::InvalidLegacyAtomListLine {
                line: 0,
                col: 9,
            }))
        );
    }

    #[rstest]
    #[case::empty(b"")]
    fn test_legacy_atom_list_block_eof_error(#[case] input: &[u8]) {
        let mut remaining = input;
        assert_eq!(
            legacy_atom_list_block(&mut remaining, 1, 3, CtabParseFlags::LENIENT),
            Err(ErrMode::Cut(ParseError::UnexpectedEof {
                line: 3,
                block: "legacy atom list",
            }))
        );
    }

    #[rstest]
    #[case::legacy_feature(b"  1 F    3   9   7   8", 0)]
    fn test_legacy_atom_list_input_extended_error(#[case] input: &[u8], #[case] column: u32) {
        let error = legacy_atom_list_input(CtabParseFlags::EXTENDED)
            .parse(Input::new(input))
            .unwrap_err()
            .into_inner();
        assert_eq!(error, InputError { column });
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
        let result = legacy_atom_list_input(CtabParseFlags::LENIENT)
            .parse(Input::new(input))
            .unwrap();
        assert_eq!(result, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::count_is_zero(b"  1 F    0   9  ", 9)]
    #[case::count_exceeds_5(b"  1 F    6    9  ", 9)]
    #[case::invalid_exclusion_flag(b"  1 X    1   9  ", 4)]
    #[case::invalid_element_atomic_number(b"  1 F    1   0  ", 11)]
    fn test_legacy_atom_list_input_lenient_error(
        #[case] input: &[u8],
        #[case] column: u32,
    ) {
        let error = legacy_atom_list_input(CtabParseFlags::LENIENT)
            .parse(Input::new(input))
            .unwrap_err()
            .into_inner();
        assert_eq!(error, InputError { column });
    }
}
