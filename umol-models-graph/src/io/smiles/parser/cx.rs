//! CXSMILES annotation block parser
//!
//! Parses the `|...|` extension block in CXSMILES format.
//! Two parsers are provided:
//! - `parse_cx_annotations`: basic annotations only (for Molecule)
//! - `parse_extended_cx_annotations`: all annotations (for ExtendedMolecule)

use std::collections::HashMap;

use bstr::ByteSlice;
use nom::branch::alt;
use nom::bytes::complete::{tag, take_till1, take_while1};
use nom::character::complete::{char, satisfy, u32 as nom_u32};
use nom::combinator::{not, opt, value};
use nom::error::{Error as NomError, ErrorKind};
use nom::multi::{many1, separated_list0};
use nom::number::complete::double;
use nom::sequence::{delimited, preceded, separated_pair, terminated};
use nom::{Err, IResult, Parser};
use umol_data::SpinMultiplicity;

use super::super::error::ParseError;
use super::utils::{split_escaped_semicolons, unescape_html_entities};
use crate::position::Point3D;
use crate::table_ir::{
    Bond, BondDonation, BondNoncovalent, BondOrder, BondStereo, BondWedge, CxAnnotationData,
    ExtendedBond, ExtendedMolecule, Molecule, StereoMode, StereoSet, StereoSetMode,
    UnpairedElectrons,
};

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

/// Update Molecule with parsed CX entries
pub fn update_molecule(mol: &mut Molecule, entries: Vec<CxEntry>) {
    for entry in entries {
        match entry {
            CxEntry::Coordinates(coords) => {
                mol.positions = Some(coords);
            }
            CxEntry::Labels(labels) => {
                for (idx, label) in labels {
                    if let Some(atom) = mol.atoms.get_mut(idx as usize) {
                        atom.label = Some(label);
                    }
                }
            }
            CxEntry::Values(values) => {
                for (idx, value) in values {
                    if let Some(atom) = mol.atoms.get_mut(idx as usize) {
                        atom.value = Some(value);
                    }
                }
            }
            CxEntry::Radicals(radicals) => {
                for (idx, unpaired) in radicals {
                    if let Some(atom) = mol.atoms.get_mut(idx as usize) {
                        atom.unpaired_electrons = Some(unpaired);
                    }
                }
            }
            CxEntry::WigglyBonds(wiggly) => {
                for (idx, wedge) in wiggly {
                    if let Some(bond) = mol.bonds.get_mut(idx as usize) {
                        bond.wedge = Some(wedge);
                    }
                }
            }
            CxEntry::CisBonds(indices) => {
                for idx in indices {
                    if let Some(bond) = mol.bonds.get_mut(idx as usize) {
                        bond.stereo = Some(BondStereo::Cis);
                    }
                }
            }
            CxEntry::TransBonds(indices) => {
                for idx in indices {
                    if let Some(bond) = mol.bonds.get_mut(idx as usize) {
                        bond.stereo = Some(BondStereo::Trans);
                    }
                }
            }
            CxEntry::CoordinateBonds(pairs) => {
                for (from, to) in pairs {
                    let mut bond = Bond::new(from, to, BondOrder::Single);
                    bond.donation = Some(BondDonation::Donating);
                    mol.bonds.push(bond);
                }
            }
            CxEntry::HydrogenBonds(pairs) => {
                for (from, to) in pairs {
                    let mut bond = Bond::new(from, to, BondOrder::Zero);
                    bond.noncovalent = Some(BondNoncovalent::Hydrogen);
                    mol.bonds.push(bond);
                }
            }
            // Extended-only entries are ignored for basic Molecule
            CxEntry::FragmentGroups(_)
            | CxEntry::StereoGroup(_)
            | CxEntry::RelativeStereo
            | CxEntry::AtomProperties(_) => {}
        }
    }
}

/// Update ExtendedMolecule with parsed CX entries
pub fn update_extended_molecule(mol: &mut ExtendedMolecule, entries: Vec<CxEntry>) {
    let mut stereo_mode: Option<StereoMode> = None;
    let mut stereo_groups: HashMap<u32, StereoSet> = HashMap::new();
    let mut components: Option<Vec<Vec<u32>>> = None;

    for entry in entries {
        match entry {
            CxEntry::Coordinates(coords) => {
                mol.positions = Some(coords);
            }
            CxEntry::Labels(labels) => {
                for (idx, label) in labels {
                    if let Some(atom) = mol.atoms.get_mut(idx as usize) {
                        atom.label = Some(label);
                    }
                }
            }
            CxEntry::Values(values) => {
                for (idx, value) in values {
                    if let Some(atom) = mol.atoms.get_mut(idx as usize) {
                        atom.value = Some(value);
                    }
                }
            }
            CxEntry::Radicals(radicals) => {
                for (idx, unpaired) in radicals {
                    if let Some(atom) = mol.atoms.get_mut(idx as usize) {
                        atom.unpaired_electrons = Some(unpaired);
                    }
                }
            }
            CxEntry::WigglyBonds(wiggly) => {
                for (idx, wedge) in wiggly {
                    if let Some(bond) = mol.bonds.get_mut(idx as usize) {
                        bond.wedge = Some(wedge);
                    }
                }
            }
            CxEntry::CisBonds(indices) => {
                for idx in indices {
                    if let Some(bond) = mol.bonds.get_mut(idx as usize) {
                        bond.stereo = Some(BondStereo::Cis);
                    }
                }
            }
            CxEntry::TransBonds(indices) => {
                for idx in indices {
                    if let Some(bond) = mol.bonds.get_mut(idx as usize) {
                        bond.stereo = Some(BondStereo::Trans);
                    }
                }
            }
            CxEntry::CoordinateBonds(pairs) => {
                for (from, to) in pairs {
                    let mut bond = ExtendedBond::new(from, to, BondOrder::Single);
                    bond.donation = Some(BondDonation::Donating);
                    mol.bonds.push(bond);
                }
            }
            CxEntry::HydrogenBonds(pairs) => {
                for (from, to) in pairs {
                    let mut bond = ExtendedBond::new(from, to, BondOrder::Zero);
                    bond.noncovalent = Some(BondNoncovalent::Hydrogen);
                    mol.bonds.push(bond);
                }
            }
            CxEntry::FragmentGroups(groups) => {
                components = Some(groups);
            }
            CxEntry::StereoGroup(sg) => match sg.group_type {
                StereoGroupType::Absolute => {
                    // Absolute atoms don't need group storage; stereo_mode captures this
                    stereo_mode = Some(StereoMode::Absolute);
                }
                StereoGroupType::Or(n) => {
                    stereo_groups
                        .entry(n)
                        .and_modify(|s| s.atoms.extend(sg.atoms.iter().copied()))
                        .or_insert(StereoSet {
                            atoms: sg.atoms,
                            mode: StereoSetMode::Correlated,
                        });
                }
                StereoGroupType::And(n) => {
                    stereo_groups
                        .entry(n)
                        .and_modify(|s| s.atoms.extend(sg.atoms.iter().copied()))
                        .or_insert(StereoSet {
                            atoms: sg.atoms,
                            mode: StereoSetMode::Independent,
                        });
                }
            },
            CxEntry::RelativeStereo => {
                stereo_mode = Some(StereoMode::Relative);
            }
            CxEntry::AtomProperties(props) => {
                for (idx, key, value) in props {
                    if let Some(atom) = mol.atoms.get_mut(idx as usize) {
                        atom.properties.insert(key, value);
                    }
                }
            }
        }
    }

    // Store CX-specific data if any
    if stereo_mode.is_some() || !stereo_groups.is_empty() || components.is_some() {
        mol.cx_data = Some(CxAnnotationData {
            stereo_mode,
            stereo_groups,
            components,
        });
    }
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
        Ok(([], options)) => Ok(options.into_iter().flatten().collect()),
        Ok(_) => Err(ParseError::InvalidToken { pos: 0 }),
        Err(Err::Failure(e)) if e.code == ErrorKind::Verify => {
            Err(ParseError::InvalidCxProperty { pos: 0 })
        }
        Err(_) => Err(ParseError::InvalidToken { pos: 0 }),
    }
}

/// Parse coordinates (x,y) or (x,y,z) for a single atom.
/// Missing components default to 0.0.
fn parse_atom_coordinates(input: &[u8]) -> IResult<&[u8], Point3D> {
    let (input, coords) = separated_list0(char(','), opt(double)).parse(input)?;
    if coords.is_empty() {
        return Ok((input, Point3D::zero()));
    }
    if coords.len() > 3 {
        return Err(Err::Failure(NomError::new(input, ErrorKind::Tag)));
    }
    let x = coords.first().copied().flatten().unwrap_or(0.0);
    let y = coords.get(1).copied().flatten().unwrap_or(0.0);
    let z = coords.get(2).copied().flatten().unwrap_or(0.0);
    Ok((input, Point3D::new(x, y, z)))
}

/// Parse coordinates block: `(x,y,z;x,y,z;...)`
/// Empty parens `()` means no atoms have coordinates.
fn parse_coordinates(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, coords) = alt((
        value(vec![], tag("()")),
        delimited(
            char('('),
            separated_list0(char(';'), parse_atom_coordinates),
            char(')'),
        ),
    ))
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
    terminated(
        char(','),
        not(satisfy(|c| {
            matches!(
                c,
                '(' | '$' | '^' | 'w' | 'c' | 't' | 'C' | 'H' | 'f' | 'a' | 'o' | '&' | 'r'
            )
        })),
    )
    .parse(input)
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
    use umol_data::Element;

    use super::*;
    use crate::table_ir::{Atom, Chirality, ExtendedAtom};

    #[fixture]
    fn triatomic_molecule() -> Molecule {
        let mut mol = Molecule::empty();
        mol.atoms = vec![
            Atom::from_element(Element::C),
            Atom::from_element(Element::N),
            Atom::from_element(Element::O),
        ];
        mol.bonds = vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Double),
        ];
        mol
    }

    #[fixture]
    fn triatomic_extended_molecule() -> ExtendedMolecule {
        let mut mol = ExtendedMolecule::empty();
        let mut atom0 = ExtendedAtom::from_element(Element::C);
        atom0.chirality = Some(Chirality::Clockwise);
        let mut atom1 = ExtendedAtom::from_element(Element::N);
        atom1.chirality = Some(Chirality::CounterClockwise);
        mol.atoms = vec![atom0, atom1, ExtendedAtom::from_element(Element::O)];
        mol.bonds = vec![
            ExtendedBond::new(0, 1, BondOrder::Single),
            ExtendedBond::new(1, 2, BondOrder::Double),
        ];
        mol
    }

    #[rstest]
    #[case::blank(b"()", CxEntry::Coordinates(vec![]))]
    #[case::empty(b"(,,)", CxEntry::Coordinates(vec![Point3D::zero()]))]
    #[case::atom_2d(b"(1.0,2.0,)", CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 0.0)]))]
    #[case::atom_2d_nocomma(b"(1,2)", CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 0.0)]))]
    #[case::atom_3d(b"(1.0,2.0,3.0)", CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0)]))]
    #[case::atom_x(b"(1.0,,)", CxEntry::Coordinates(vec![Point3D::new(1.0, 0.0, 0.0)]))]
    #[case::atom_x_nocomma(b"(1)", CxEntry::Coordinates(vec![Point3D::new(1.0, 0.0, 0.0)]))]
    #[case::atom_y(b"(,2.0,)", CxEntry::Coordinates(vec![Point3D::new(0.0, 2.0, 0.0)]))]
    #[case::atom_y_nocomma(b"(,2)", CxEntry::Coordinates(vec![Point3D::new(0.0, 2.0, 0.0)]))]
    #[case::atom_z(b"(,,3.0)", CxEntry::Coordinates(vec![Point3D::new(0.0, 0.0, 3.0)]))]
    #[case::two_atoms_1(b"(;1,2)", CxEntry::Coordinates(vec![Point3D::new(0.0, 0.0, 0.0), Point3D::new(1.0, 2.0, 0.0)]))]
    #[case::two_atoms_2(b"(1,2;)", CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 0.0), Point3D::new(0.0, 0.0, 0.0)]))]
    fn test_parse_coordinates(#[case] input: &[u8], #[case] expected: CxEntry) {
        let result = parse_coordinates(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded: {:?}",
            input_str,
            result
        );
        let (_, entries) = result.unwrap();
        assert_eq!(
            entries, expected,
            "{:?} should have parsed to {:?}",
            input_str, entries
        );
    }

    #[rstest]
    #[case::atom_4d(b"(1.0,2.0,3.0,4.0)", ErrorKind::Tag)]
    fn test_parse_coordinates_invalid(#[case] input: &[u8], #[case] expected_kind: ErrorKind) {
        let result = parse_coordinates(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_err(),
            "{:?} should have failed: {:?}",
            input_str,
            result
        );
        assert!(
            matches!(result.clone(), Err(Err::Failure(e)) if e.code == expected_kind),
            "{:?} should have failed with error kind {:?}, got {:?}",
            input_str,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code)
        );
    }

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

    #[rustfmt::skip]
    #[rstest]
    #[case::coordinates(vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0)])], |mol: &Molecule| mol.positions == Some(vec![Point3D::new(1.0, 2.0, 3.0)]))]
    #[case::labels(vec![CxEntry::Labels(vec![(0, "C1".to_string()), (1, "N1".to_string())])], |mol: &Molecule| mol.atoms[0].label == Some("C1".to_string()) && mol.atoms[1].label == Some("N1".to_string()))]
    #[case::values(vec![CxEntry::Values(vec![(0, "val0".to_string())])], |mol: &Molecule| mol.atoms[0].value == Some("val0".to_string()))]
    #[case::radicals(vec![CxEntry::Radicals(vec![(0, UnpairedElectrons { count: 1, multiplicity: None })])],
        |mol: &Molecule| mol.atoms[0].unpaired_electrons == Some(UnpairedElectrons { count: 1, multiplicity: None }))]
    #[case::wiggly_bonds(vec![CxEntry::WigglyBonds(vec![(0, BondWedge::Either)])], |mol: &Molecule| mol.bonds[0].wedge == Some(BondWedge::Either))]
    #[case::cis_bonds(vec![CxEntry::CisBonds(vec![0])], |mol: &Molecule| mol.bonds[0].stereo == Some(BondStereo::Cis))]
    #[case::trans_bonds(vec![CxEntry::TransBonds(vec![1])], |mol: &Molecule| mol.bonds[1].stereo == Some(BondStereo::Trans))]
    #[case::coordinate_bonds(vec![CxEntry::CoordinateBonds(vec![(0, 2)])], |mol: &Molecule| mol.bonds.len() == 3 && mol.bonds[2].donation == Some(BondDonation::Donating))]
    #[case::hydrogen_bonds(vec![CxEntry::HydrogenBonds(vec![(0, 2)])], |mol: &Molecule| mol.bonds.len() == 3 && mol.bonds[2].noncovalent == Some(BondNoncovalent::Hydrogen))]
    #[case::extended_entries(vec![CxEntry::FragmentGroups(vec![vec![0, 1]]), CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Absolute, atoms: vec![0] }),
        CxEntry::RelativeStereo, CxEntry::AtomProperties(vec![(0, "k".to_string(), "v".to_string())])], |mol: &Molecule| mol.atoms[0].label.is_none())]
    fn test_update_molecule(
        triatomic_molecule: Molecule,
        #[case] entries: Vec<CxEntry>,
        #[case] check: fn(&Molecule) -> bool,
    ) {
        let mut mol = triatomic_molecule;
        update_molecule(&mut mol, entries);
        assert!(check(&mol));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::coordinates(vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0)])], |mol: &ExtendedMolecule| mol.positions == Some(vec![Point3D::new(1.0, 2.0, 3.0)]))]
    #[case::labels(vec![CxEntry::Labels(vec![(0, "C1".to_string())])], |mol: &ExtendedMolecule| mol.atoms[0].label == Some("C1".to_string()))]
    #[case::values(vec![CxEntry::Values(vec![(1, "val1".to_string())])], |mol: &ExtendedMolecule| mol.atoms[1].value == Some("val1".to_string()))]
    #[case::radicals(vec![CxEntry::Radicals(vec![(2, UnpairedElectrons { count: 2, multiplicity: None })])],
        |mol: &ExtendedMolecule| mol.atoms[2].unpaired_electrons == Some(UnpairedElectrons { count: 2, multiplicity: None }))]
    #[case::wiggly_bonds(vec![CxEntry::WigglyBonds(vec![(1, BondWedge::Either)])], |mol: &ExtendedMolecule| mol.bonds[1].wedge == Some(BondWedge::Either))]
    #[case::cis_bonds(vec![CxEntry::CisBonds(vec![1])], |mol: &ExtendedMolecule| mol.bonds[1].stereo == Some(BondStereo::Cis))]
    #[case::trans_bonds(vec![CxEntry::TransBonds(vec![0])], |mol: &ExtendedMolecule| mol.bonds[0].stereo == Some(BondStereo::Trans))]
    #[case::coordinate_bonds(vec![CxEntry::CoordinateBonds(vec![(1, 2)])], |mol: &ExtendedMolecule| mol.bonds.len() == 3 && mol.bonds[2].donation == Some(BondDonation::Donating))]
    #[case::hydrogen_bonds(vec![CxEntry::HydrogenBonds(vec![(0, 2)])], |mol: &ExtendedMolecule| mol.bonds.len() == 3 && mol.bonds[2].noncovalent == Some(BondNoncovalent::Hydrogen))]
    #[case::fragment_groups(vec![CxEntry::FragmentGroups(vec![vec![0, 1], vec![2]])], |mol: &ExtendedMolecule| mol.cx_data.as_ref().map(|d| d.components.as_ref()) == Some(Some(&vec![vec![0, 1], vec![2]])))]
    #[case::stereo_group_absolute(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Absolute, atoms: vec![0, 1] })],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_mode) == Some(StereoMode::Absolute))]
    #[case::stereo_group_or(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Or(1), atoms: vec![0] })],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_groups.get(&1)) == Some(&StereoSet { atoms: vec![0], mode: StereoSetMode::Correlated }))]
    #[case::stereo_group_and(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::And(2), atoms: vec![1] })],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_groups.get(&2)) == Some(&StereoSet { atoms: vec![1], mode: StereoSetMode::Independent }))]
    #[case::relative_stereo(vec![CxEntry::RelativeStereo],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_mode) == Some(StereoMode::Relative))]
    #[case::atom_properties(vec![CxEntry::AtomProperties(vec![(0, "key".to_string(), "value".to_string())])], |mol: &ExtendedMolecule| mol.atoms[0].properties.get("key") == Some(&"value".to_string()))]
    fn test_update_extended_molecule(
        triatomic_extended_molecule: ExtendedMolecule,
        #[case] entries: Vec<CxEntry>,
        #[case] check: fn(&ExtendedMolecule) -> bool,
    ) {
        let mut mol = triatomic_extended_molecule;
        update_extended_molecule(&mut mol, entries);
        assert!(check(&mol));
    }
}
