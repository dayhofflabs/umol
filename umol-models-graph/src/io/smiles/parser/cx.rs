//! CXSMILES annotation block parser
//!
//! Parses the `|...|` extension block in CXSMILES format.
//! Two parsers are provided:
//! - `parse_cx_annotations`: basic annotations only (for Molecule)
//! - `parse_extended_cx_annotations`: all annotations (for ExtendedMolecule)

#![allow(dead_code)]

use bstr::ByteSlice;
use nom::branch::alt;
use nom::bytes::complete::{tag, take_till1, take_while1};
use nom::character::complete::{char, u32 as nom_u32};
use nom::combinator::{opt, value};
use nom::error::{Error as NomError, ErrorKind};
use nom::multi::{many1, separated_list0};
use nom::number::complete::double;
use nom::sequence::{delimited, preceded, separated_pair};
use nom::{Err, IResult, Parser};
use umol_data::SpinMultiplicity;

use super::super::error::ParseError;
use super::utils::{split_escaped_semicolons, unescape_html_entities};
use crate::position::Point3D;
use crate::table_ir::{BondWedge, UnpairedElectrons};

/// Stereo group type for enhanced stereochemistry
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum StereoGroupType {
    /// Absolute stereochemistry (as drawn)
    Absolute,
    /// OR group - molecule is one of the stereoisomers (all centers flip together)
    Or(u32),
    /// AND group - mixture of stereoisomers (centers are independent)
    And(u32),
}

/// Enhanced stereo group
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoGroup {
    pub group_type: StereoGroupType,
    pub atoms: Vec<u32>,
}

/// A parsed CXSMILES annotation entry
#[derive(Clone, Debug, PartialEq)]
pub enum CxEntry {
    /// Atom coordinates: (x,y,z;...)
    Coordinates(Vec<Point3D>),
    /// Atom labels: $label;label;...$
    Labels(Vec<(u32, String)>),
    /// Atom values: $_AV:value;value;...$
    Values(Vec<(u32, String)>),
    /// Radical electrons: ^n:idx,idx,...
    Radicals(Vec<(u32, UnpairedElectrons)>),
    /// Wiggly bonds: w:, wU:, wD:
    WigglyBonds(Vec<(u32, BondWedge)>),
    /// Cis double bonds: c:
    CisBonds(Vec<u32>),
    /// Trans double bonds: t:
    TransBonds(Vec<u32>),
    /// Coordinate (dative) bonds: C:
    CoordinateBonds(Vec<(u32, u32)>),
    /// Hydrogen bonds: H:
    HydrogenBonds(Vec<(u32, u32)>),
    /// Fragment grouping: f: (extended only)
    FragmentGroups(Vec<Vec<u32>>),
    /// Enhanced stereo group: a:, o<n>:, &<n>: (extended only)
    StereoGroup(StereoGroup),
    /// Relative stereo marker: r (extended only)
    RelativeStereo,
    /// Atom properties: atomProp: (extended only)
    AtomProperties(Vec<(u32, String, String)>),
}

/// Parse basic CX annotations (for Molecule)
pub fn parse_cx_annotations(input: &[u8]) -> Result<Vec<CxEntry>, ParseError> {
    parse_cx_block(input, parse_basic_entry)
}

/// Parse extended CX annotations (for ExtendedMolecule)
pub fn parse_extended_cx_annotations(input: &[u8]) -> Result<Vec<CxEntry>, ParseError> {
    parse_cx_block(input, parse_extended_entry)
}

fn parse_cx_block<'inp>(
    input: &'inp [u8],
    entry_parser: impl Parser<&'inp [u8], Output = Option<CxEntry>, Error = NomError<&'inp [u8]>>,
) -> Result<Vec<CxEntry>, ParseError> {
    match delimited(
        opt(char('|')),
        separated_list0(char(','), entry_parser),
        opt(char('|')),
    )
    .parse(input)
    {
        Ok((remaining, options)) if remaining.is_empty() => {
            Ok(options.into_iter().flatten().collect())
        }
        Ok(_) => Err(ParseError::InvalidToken { pos: 0 }),
        Err(Err::Failure(e)) if e.code == ErrorKind::Verify => {
            Err(ParseError::InvalidCxProperty { pos: 0 })
        }
        Err(_) => Err(ParseError::InvalidToken { pos: 0 }),
    }
}

/// Parse coordinates (x,y) or (x,y,z) for a single atom, or empty for missing.
fn parse_atom_coordinates(input: &[u8]) -> IResult<&[u8], Point3D> {
    if input.is_empty() {
        return Ok((input, Point3D::new(f64::NAN, f64::NAN, f64::NAN)));
    }

    let (remaining, (x, y, z)) = (
        double,
        preceded(char(','), double),
        opt(preceded(char(','), double)),
    )
        .parse(input)?;

    Ok((remaining, Point3D::new(x, y, z.unwrap_or(0.0))))
}

/// Parse coordinates block: `(x,y,z;x,y,z;...)`
fn parse_coordinates(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, coords) = delimited(
        char('('),
        separated_list0(char(';'), parse_atom_coordinates),
        char(')'),
    )
    .parse(input)?;

    Ok((input, CxEntry::Coordinates(coords)))
}

/// Parse labels `$label;label;...$` or values `$_AV:value;value;...$`
fn parse_labels(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, inner) =
        delimited(char('$'), take_while1(|b| b != b'$'), char('$')).parse(input)?;

    let (is_values, data) = match inner.strip_prefix(b"_AV:") {
        Some(rest) => (true, rest),
        None => (false, inner),
    };

    let entries = split_escaped_semicolons(data);
    let result: Vec<_> = entries
        .into_iter()
        .enumerate()
        .filter(|(_, e)| !e.is_empty())
        .map(|(idx, e)| {
            (
                idx as u32,
                unescape_html_entities(e).to_str_lossy().into_owned(),
            )
        })
        .collect();

    if is_values {
        Ok((input, CxEntry::Values(result)))
    } else {
        Ok((input, CxEntry::Labels(result)))
    }
}

/// Convert CXSMILES radical code (1-7) to unpaired electrons.
fn convert_radical_code(code: u8) -> UnpairedElectrons {
    match code {
        1 => UnpairedElectrons::from_count(1),
        2 => UnpairedElectrons::from_count(2),
        3 => UnpairedElectrons::new(2, Some(SpinMultiplicity::Singlet)),
        4 => UnpairedElectrons::new(2, Some(SpinMultiplicity::Triplet)),
        5 => UnpairedElectrons::from_count(3),
        6 => UnpairedElectrons::new(3, Some(SpinMultiplicity::Doublet)),
        7 => UnpairedElectrons::new(3, Some(SpinMultiplicity::Quartet)),
        _ => UnpairedElectrons::from_count(1),
    }
}

/// Parse a single radical group: `^n:idx,idx,...`
fn parse_radical_group(input: &[u8]) -> IResult<&[u8], (u8, Vec<u32>)> {
    let (input, code) = delimited(
        char('^'),
        nom::character::complete::one_of("1234567"),
        char(':'),
    )
    .parse(input)?;
    let (input, indices) = separated_list0(comma_not_before_entry, nom_u32).parse(input)?;
    Ok((input, (code as u8 - b'0', indices)))
}

/// Parse radicals: `^n:idx,idx,...` (one or more groups).
fn parse_radicals(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, groups) = many1(parse_radical_group).parse(input)?;

    let result: Vec<_> = groups
        .into_iter()
        .flat_map(|(code, indices)| {
            let unpaired = convert_radical_code(code);
            indices.into_iter().map(move |idx| (idx, unpaired))
        })
        .collect();

    Ok((input, CxEntry::Radicals(result)))
}

/// Parse wiggly bonds: `w:`, `wU:`, `wD:` followed by atom.bond pairs.
fn parse_wiggly_bonds(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, wedge_type) = alt((
        value(BondWedge::EitherUp, tag("wU:")),
        value(BondWedge::EitherDown, tag("wD:")),
        value(BondWedge::Either, tag("w:")),
    ))
    .parse(input)?;

    let (input, pairs) = separated_list0(
        comma_not_before_entry,
        separated_pair(nom_u32, char('.'), nom_u32),
    )
    .parse(input)?;

    let result: Vec<_> = pairs.into_iter().map(|(a, _b)| (a, wedge_type)).collect();
    Ok((input, CxEntry::WigglyBonds(result)))
}

/// Parse cis/trans bond annotations: `c:`, `t:`, `ctu:`.
fn parse_cis_trans(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, is_cis) = alt((
        value(None, tag("ctu:")),
        value(Some(true), tag("c:")),
        value(Some(false), tag("t:")),
    ))
    .parse(input)?;

    let (input, indices) = separated_list0(comma_not_before_entry, nom_u32).parse(input)?;

    match is_cis {
        Some(true) => Ok((input, CxEntry::CisBonds(indices))),
        Some(false) => Ok((input, CxEntry::TransBonds(indices))),
        None => Ok((input, CxEntry::CisBonds(Vec::new()))),
    }
}

/// Parse coordinate (dative) bonds: `C:atom.bond,...`
fn parse_coordinate_bonds(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, pairs) = preceded(
        tag("C:"),
        separated_list0(
            comma_not_before_entry,
            separated_pair(nom_u32, char('.'), nom_u32),
        ),
    )
    .parse(input)?;

    Ok((input, CxEntry::CoordinateBonds(pairs)))
}

/// Parse hydrogen bonds: `H:atom.bond,...`
fn parse_hydrogen_bonds(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, pairs) = preceded(
        tag("H:"),
        separated_list0(
            comma_not_before_entry,
            separated_pair(nom_u32, char('.'), nom_u32),
        ),
    )
    .parse(input)?;

    Ok((input, CxEntry::HydrogenBonds(pairs)))
}

/// Parse fragment groups: `f:atom.atom.atom,...`
fn parse_fragment_groups(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let parse_group = separated_list0(char('.'), nom_u32);
    let (input, groups) = preceded(
        tag("f:"),
        separated_list0(comma_not_before_entry, parse_group),
    )
    .parse(input)?;

    let non_empty: Vec<_> = groups.into_iter().filter(|g| !g.is_empty()).collect();
    Ok((input, CxEntry::FragmentGroups(non_empty)))
}

/// Parse absolute stereo group: `a:idx,idx,...`
fn parse_stereo_absolute(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, atoms) =
        preceded(tag("a:"), separated_list0(comma_not_before_entry, nom_u32)).parse(input)?;

    Ok((
        input,
        CxEntry::StereoGroup(StereoGroup {
            group_type: StereoGroupType::Absolute,
            atoms,
        }),
    ))
}

/// Parse OR/AND stereo group: `o<n>:idx,idx,...` or `&<n>:idx,idx,...`
fn parse_stereo_or_and(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, (is_or, group_num, _, atoms)) = (
        alt((value(true, char('o')), value(false, char('&')))),
        nom_u32,
        char(':'),
        separated_list0(comma_not_before_entry, nom_u32),
    )
        .parse(input)?;

    let group_type = if is_or {
        StereoGroupType::Or(group_num)
    } else {
        StereoGroupType::And(group_num)
    };

    Ok((
        input,
        CxEntry::StereoGroup(StereoGroup { group_type, atoms }),
    ))
}

/// Parse a single atom property entry: `idx.key.value`
fn parse_atom_prop_entry(input: &[u8]) -> IResult<&[u8], (u32, String, String)> {
    let (input, (idx, _, key_bytes, _, value_bytes)) = (
        nom_u32,
        char('.'),
        take_while1(|b| b != b'.'),
        char('.'),
        take_while1(|b| b != b':' && b != b',' && b != b'|'),
    )
        .parse(input)?;

    let key = unescape_html_entities(key_bytes)
        .to_str_lossy()
        .into_owned();
    let value = unescape_html_entities(value_bytes)
        .to_str_lossy()
        .into_owned();

    Ok((input, (idx, key, value)))
}

/// Parse atom properties: `atomProp:idx.key.value:idx.key.value...`
fn parse_atom_properties(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, props) = preceded(
        tag("atomProp:"),
        separated_list0(char(':'), parse_atom_prop_entry),
    )
    .parse(input)?;

    Ok((input, CxEntry::AtomProperties(props)))
}

/// Parse relative stereo marker: `r` or `r:idx,idx,...`
fn parse_relative_stereo(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, _) = char('r').parse(input)?;

    if input.first() == Some(&b':') {
        let (input, _) = preceded(char(':'), separated_list0(char(','), nom_u32)).parse(input)?;
        Ok((input, CxEntry::RelativeStereo))
    } else {
        Ok((input, CxEntry::RelativeStereo))
    }
}

/// Check for extended-only feature markers.
fn parse_extended_marker(input: &[u8]) -> IResult<&[u8], ()> {
    value(
        (),
        alt((
            tag("f:"),
            tag("a:"),
            tag("atomProp:"),
            tag("Sg:"),
            tag("RG:"),
            tag("rb:"),
            tag("s:"),
            tag("u:"),
            tag("LN:"),
            tag("LO:"),
        )),
    )
    .parse(input)
}

/// Check for `r` as standalone extended marker.
fn parse_extended_r(input: &[u8]) -> IResult<&[u8], ()> {
    let (rest, _) = char('r').parse(input)?;
    if rest.is_empty() || matches!(rest[0], b',' | b'|' | b':') {
        Ok((rest, ()))
    } else {
        Err(Err::Error(NomError::new(input, ErrorKind::Tag)))
    }
}

/// Check for `o<n>:` or `&<n>:` extended markers.
fn parse_extended_o_and(input: &[u8]) -> IResult<&[u8], ()> {
    let (rest, _) = alt((char('o'), char('&'))).parse(input)?;
    let (rest, _) = nom_u32(rest)?;
    let (rest, _) = char(':').parse(rest)?;
    Ok((rest, ()))
}

/// Skip over an unknown/unrecognized entry until next delimiter.
fn skip_unknown_entry(input: &[u8]) -> IResult<&[u8], ()> {
    let (input, _) = take_till1(|b| b == b',' || b == b'|').parse(input)?;
    Ok((input, ()))
}

/// Parse comma only if not followed by an entry-start character.
fn comma_not_before_entry(input: &[u8]) -> IResult<&[u8], char> {
    let (rest, c) = char(',').parse(input)?;
    if is_entry_start(rest) {
        Err(Err::Error(NomError::new(input, ErrorKind::Char)))
    } else {
        Ok((rest, c))
    }
}

fn is_entry_start(input: &[u8]) -> bool {
    input
        .first()
        .map(|&b| is_entry_start_byte(b))
        .unwrap_or(false)
}

fn is_entry_start_byte(b: u8) -> bool {
    matches!(
        b,
        b'(' | b'$' | b'^' | b'w' | b'c' | b't' | b'C' | b'H' | b'f' | b'a' | b'o' | b'&' | b'r'
    )
}

/// Parse a basic CX entry (rejects extended features).
fn parse_basic_entry(input: &[u8]) -> IResult<&[u8], Option<CxEntry>> {
    if let Ok((rest, entry)) = alt((
        parse_coordinates,
        parse_labels,
        parse_radicals,
        parse_wiggly_bonds,
        parse_cis_trans,
        parse_coordinate_bonds,
        parse_hydrogen_bonds,
    ))
    .parse(input)
    {
        return Ok((rest, Some(entry)));
    }

    if parse_extended_marker(input).is_ok()
        || parse_extended_r(input).is_ok()
        || parse_extended_o_and(input).is_ok()
    {
        return Err(Err::Failure(NomError::new(input, ErrorKind::Verify)));
    }

    let (rest, _) = skip_unknown_entry(input)?;
    Ok((rest, None))
}

/// Parse an extended CX entry (all features allowed).
fn parse_extended_entry(input: &[u8]) -> IResult<&[u8], Option<CxEntry>> {
    if let Ok((rest, entry)) = alt((
        parse_coordinates,
        parse_labels,
        parse_radicals,
        parse_wiggly_bonds,
        parse_cis_trans,
        parse_coordinate_bonds,
        parse_hydrogen_bonds,
        parse_fragment_groups,
        parse_stereo_absolute,
        parse_stereo_or_and,
        parse_atom_properties,
        parse_relative_stereo,
    ))
    .parse(input)
    {
        return Ok((rest, Some(entry)));
    }

    let (rest, _) = skip_unknown_entry(input)?;
    Ok((rest, None))
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(b"||", vec![])]
    #[case::coordinate_bond(b"|C:0.1|", vec![CxEntry::CoordinateBonds(vec![(0, 1)])])]
    #[case::coordinate_bond_multiple(b"|C:0.1,2.3|", vec![CxEntry::CoordinateBonds(vec![(0, 1), (2, 3)])])]
    #[case::hydrogen_bond(b"|H:1.2|", vec![CxEntry::HydrogenBonds(vec![(1, 2)])])]
    #[case::radicals_multiple_atoms(b"|^1:0,1,2|", vec![CxEntry::Radicals(vec![(0, UnpairedElectrons { count: 1, multiplicity: None }),
        (1, UnpairedElectrons { count: 1, multiplicity: None }), (2, UnpairedElectrons { count: 1, multiplicity: None })])])]
    #[case::wiggly_bonds(b"|w:0.1,2.3|", vec![CxEntry::WigglyBonds(vec![(0, BondWedge::Either), (2, BondWedge::Either)])])]
    #[case::cis_trans(b"|c:0,1|", vec![CxEntry::CisBonds(vec![0, 1])])]
    #[case::trans_trans(b"|t:0,1|", vec![CxEntry::TransBonds(vec![0, 1])])]
    #[case::atom_labels(b"|$label1;label2;label3$|", vec![CxEntry::Labels(vec![(0, "label1".to_string()), (1, "label2".to_string()), (2, "label3".to_string())])])]
    #[case::atom_values(b"$_AV:value1;value2;value3$|", vec![CxEntry::Values(vec![(0, "value1".to_string()), (1, "value2".to_string()), (2, "value3".to_string())])])]
    #[case::coordinates_2d(b"|(1.5,2.5;3.5,4.5)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.5, 2.5, 0.0), Point3D::new(3.5, 4.5, 0.0)])])]
    #[case::coordinates_3d(b"|(1,2,3;4,5,6)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0), Point3D::new(4.0, 5.0, 6.0)])])]
    #[case::combined_entries(b"|^1:0,1,(1.0,2.0;3.0,4.0),C:2.3|", vec![CxEntry::Radicals(vec![(0, UnpairedElectrons { count: 1, multiplicity: None }),
        (1, UnpairedElectrons { count: 1, multiplicity: None })]), CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 0.0), Point3D::new(3.0, 4.0, 0.0)]), CxEntry::CoordinateBonds(vec![(2, 3)])])]
    #[case::unknown_tag_skipped(b"|xyz:123,C:0.1|", vec![CxEntry::CoordinateBonds(vec![(0, 1)])])]
    #[case::unknown_tag_only(b"|unknown|", vec![])]
    fn test_parse_cx_annotations(#[case] input: &[u8], #[case] expected: Vec<CxEntry>) {
        let result = parse_cx_annotations(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_ok(), "{:?} should have succeeded: {:?}", input_str, result);
        let entries = result.unwrap();
        assert_eq!(entries, expected, "{:?} should have parsed to {:?}", input_str, entries);
    }

    #[rstest]
    #[case::extended_feature_f(b"|f:0.1|", ParseError::InvalidCxProperty { pos: 0 })]
    #[case::extended_feature_a(b"|a:0,1|", ParseError::InvalidCxProperty { pos: 0 })]
    #[case::extended_feature_o(b"|o1:0,1|", ParseError::InvalidCxProperty { pos: 0 })]
    #[case::extended_feature_and(b"|&1:0,1|", ParseError::InvalidCxProperty { pos: 0 })]
    #[case::extended_feature_r(b"|r|", ParseError::InvalidCxProperty { pos: 0 })]
    #[case::extended_feature_atomprop(b"|atomProp:0.key.value|", ParseError::InvalidCxProperty { pos: 0 })]
    fn test_parse_cx_annotations_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let result = parse_cx_annotations(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_err(),
            "{:?} should have failed: {:?}",
            input_str,
            result
        );
        let error = result.unwrap_err();
        assert_eq!(
            error, expected,
            "{:?} should have returned an error: {:?}",
            input_str, expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(b"||", vec![])]
    #[case::coordinate_bond(b"|C:0.1|", vec![CxEntry::CoordinateBonds(vec![(0, 1)])])]
    #[case::hydrogen_bond(b"|H:1.2|", vec![CxEntry::HydrogenBonds(vec![(1, 2)])])]
    #[case::radicals(b"|^1:0|", vec![CxEntry::Radicals(vec![(0, UnpairedElectrons { count: 1, multiplicity: None })])])]
    #[case::wiggly_bonds(b"|w:0.1|", vec![CxEntry::WigglyBonds(vec![(0, BondWedge::Either)])])]
    #[case::cis_bonds(b"|c:0|", vec![CxEntry::CisBonds(vec![0])])]
    #[case::trans_bonds(b"|t:0|", vec![CxEntry::TransBonds(vec![0])])]
    #[case::atom_labels(b"|$label$|", vec![CxEntry::Labels(vec![(0, "label".to_string())])])]
    #[case::fragment_groups(b"|f:0.1.2,3.4|", vec![CxEntry::FragmentGroups(vec![vec![0, 1, 2], vec![3, 4]])])]
    #[case::stereo_absolute(b"|a:0,1,2|", vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Absolute, atoms: vec![0, 1, 2] })])]
    #[case::stereo_or(b"|o1:0,1|", vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Or(1), atoms: vec![0, 1] })])]
    #[case::stereo_and(b"|&1:0,1|", vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::And(1), atoms: vec![0, 1] })])]
    #[case::relative_stereo(b"|r|", vec![CxEntry::RelativeStereo])]
    #[case::atom_properties(b"|atomProp:0.key.value|", vec![CxEntry::AtomProperties(vec![(0, "key".to_string(), "value".to_string())])])]
    #[case::unknown_tag_skipped(b"|xyz:123,C:0.1|", vec![CxEntry::CoordinateBonds(vec![(0, 1)])])]
    #[case::coordinates_2d(b"|(1.5,2.5;3.5,4.5)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.5, 2.5, 0.0), Point3D::new(3.5, 4.5, 0.0)])])]
    #[case::coordinates_3d(b"|(1,2,3;4,5,6)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0), Point3D::new(4.0, 5.0, 6.0)])])]
    fn test_parse_extended_cx_annotations(#[case] input: &[u8], #[case] expected: Vec<CxEntry>) {
        let result = parse_extended_cx_annotations(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded: {:?}",
            input_str,
            result
        );
        let entries = result.unwrap();
        assert_eq!(
            entries, expected,
            "{:?} should have parsed to {:?}",
            input_str, entries
        );
    }
}
