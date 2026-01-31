//! CXSMILES annotation block parser
//!
//! Parses the `|...|` extension block in CXSMILES format.
//! Two parsers are provided:
//! - `parse_cx_annotations`: basic annotations only (for Molecule)
//! - `parse_extended_cx_annotations`: all annotations (for ExtendedMolecule)

use std::collections::HashMap;

use bstr::ByteSlice;
use nom::branch::alt;
use nom::bytes::complete::{tag, take_while1};
use nom::character::complete::{char, one_of, satisfy, u32 as nom_u32};
use nom::combinator::{not, opt, value};
use nom::error::{Error as NomError, ErrorKind};
use nom::multi::{many1, separated_list0};
use nom::number::complete::double;
use nom::sequence::{delimited, preceded, separated_pair, terminated};
use nom::{Err, IResult, Parser};
use umol_data::SpinMultiplicity;

use super::super::config::SmilesParseFlags;
use super::super::error::ParseError;
use super::utils::{split_escaped_semicolons, unescape_html_entities};
use crate::position::Point3D;
use crate::table_ir::{
    BondDonation, BondNoncovalent, BondOrder, BondStereo, BondWedge, CxAnnotationData,
    ExtendedMolecule, Molecule, StereoInterpretation, StereoSet, StereoSetMode, UnpairedElectrons,
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
    /// Wiggly bonds: w:, wU:, wD: encoded as `<atom_idx>.<bond_idx>`
    WigglyBonds(Vec<(u32, u32, BondWedge)>),
    /// Cis double bonds: c:
    CisBonds(Vec<u32>),
    /// Trans double bonds: t:
    TransBonds(Vec<u32>),
    /// Unspecified (either) double bonds: ctu:
    UnspecBonds(Vec<u32>),
    /// Coordinate (dative) bonds: C: encoded as `<first_atom_idx>.<bond_idx>`
    CoordinateBonds(Vec<(u32, u32)>),
    /// Hydrogen bonds: H: encoded as `<first_atom_idx>.<bond_idx>`
    HydrogenBonds(Vec<(u32, u32)>),
    /// Fragment grouping: f: (extended only)
    FragmentGroups(Vec<Vec<u32>>),
    /// Enhanced stereo group: a:, o<n>:, &<n>: (extended only)
    StereoGroup(StereoGroup),
    /// Relative stereo tag: r (extended only)
    RelativeStereo,
    /// Atom properties: atomProp: (extended only)
    AtomProperties(Vec<(u32, String, String)>),
}

/// Parse basic CX annotations (for Molecule)
pub fn parse_cx_annotations(
    input: &[u8],
    flags: SmilesParseFlags,
) -> Result<Vec<CxEntry>, ParseError> {
    let skip_unknown_cx_tags = flags.contains(SmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS);
    parse_cx_block(input, |i| parse_basic_entry(i, skip_unknown_cx_tags))
}

/// Parse extended CX annotations (for ExtendedMolecule)
pub fn parse_extended_cx_annotations(
    input: &[u8],
    flags: SmilesParseFlags,
) -> Result<Vec<CxEntry>, ParseError> {
    let skip_unknown_cx_tags = flags.contains(SmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS);
    parse_cx_block(input, |i| parse_extended_entry(i, skip_unknown_cx_tags))
}

/// Update Molecule with parsed CX entries
pub fn update_molecule(mol: &mut Molecule, entries: Vec<CxEntry>) -> Result<(), ParseError> {
    for entry in entries {
        match entry {
            CxEntry::Coordinates(coords) => {
                if coords.len() > mol.atoms.len() {
                    return Err(ParseError::AtomIndexOutOfBounds {
                        atom_idx: mol.atoms.len() as u32,
                    });
                }
                mol.positions = Some(coords);
            }
            CxEntry::Labels(labels) => {
                for (idx, label) in labels {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.label = Some(label);
                }
            }
            CxEntry::Values(values) => {
                for (idx, value) in values {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.value = Some(value);
                }
            }
            CxEntry::Radicals(radicals) => {
                for (idx, unpaired) in radicals {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.unpaired_electrons = Some(unpaired);
                }
            }
            CxEntry::WigglyBonds(wiggly) => {
                for (atom_idx, bond_idx, wedge) in wiggly {
                    if atom_idx as usize >= mol.atoms.len() {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    }
                    let Some(bond) = mol.bonds.get_mut(bond_idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx });
                    };
                    let (a, b) = bond.atoms.as_tuple();
                    if atom_idx != a && atom_idx != b {
                        return Err(ParseError::MismatchedAtomBondIndices { atom_idx, bond_idx });
                    }
                    bond.wedge = Some(wedge);
                }
            }
            CxEntry::CisBonds(indices) => {
                for idx in indices {
                    let Some(bond) = mol.bonds.get_mut(idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: idx });
                    };
                    bond.stereo = Some(BondStereo::Cis);
                }
            }
            CxEntry::TransBonds(indices) => {
                for idx in indices {
                    let Some(bond) = mol.bonds.get_mut(idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: idx });
                    };
                    bond.stereo = Some(BondStereo::Trans);
                }
            }
            CxEntry::UnspecBonds(indices) => {
                for idx in indices {
                    let Some(bond) = mol.bonds.get_mut(idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: idx });
                    };
                    bond.stereo = Some(BondStereo::Either);
                }
            }
            CxEntry::CoordinateBonds(pairs) => {
                for (first_atom, bond_idx) in pairs {
                    if first_atom as usize >= mol.atoms.len() {
                        return Err(ParseError::AtomIndexOutOfBounds {
                            atom_idx: first_atom,
                        });
                    }
                    let Some(bond) = mol.bonds.get_mut(bond_idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx });
                    };
                    let (a, b) = bond.atoms.as_tuple();
                    if first_atom == a {
                        bond.donation = Some(BondDonation::Donating);
                    } else if first_atom == b {
                        bond.donation = Some(BondDonation::Accepting);
                    } else {
                        return Err(ParseError::MismatchedAtomBondIndices {
                            atom_idx: first_atom,
                            bond_idx,
                        });
                    }
                }
            }
            CxEntry::HydrogenBonds(pairs) => {
                for (first_atom, bond_idx) in pairs {
                    if first_atom as usize >= mol.atoms.len() {
                        return Err(ParseError::AtomIndexOutOfBounds {
                            atom_idx: first_atom,
                        });
                    }
                    let Some(bond) = mol.bonds.get_mut(bond_idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx });
                    };
                    let (a, b) = bond.atoms.as_tuple();
                    if first_atom != a && first_atom != b {
                        return Err(ParseError::MismatchedAtomBondIndices {
                            atom_idx: first_atom,
                            bond_idx,
                        });
                    }
                    bond.order = BondOrder::Zero;
                    bond.noncovalent = Some(BondNoncovalent::Hydrogen);
                }
            }
            _ => {}
        }
    }

    Ok(())
}

/// Update ExtendedMolecule with parsed CX entries
pub fn update_extended_molecule(
    mol: &mut ExtendedMolecule,
    entries: Vec<CxEntry>,
) -> Result<(), ParseError> {
    let mut stereo_interpretation: Option<StereoInterpretation> = None;
    let mut stereo_groups: HashMap<u32, StereoSet> = HashMap::new();
    let mut components: Option<Vec<Vec<u32>>> = None;

    for entry in entries {
        match entry {
            CxEntry::Coordinates(coords) => {
                if coords.len() > mol.atoms.len() {
                    return Err(ParseError::AtomIndexOutOfBounds {
                        atom_idx: mol.atoms.len() as u32,
                    });
                }
                mol.positions = Some(coords);
            }
            CxEntry::Labels(labels) => {
                for (idx, label) in labels {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.label = Some(label);
                }
            }
            CxEntry::Values(values) => {
                for (idx, value) in values {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.value = Some(value);
                }
            }
            CxEntry::Radicals(radicals) => {
                for (idx, unpaired) in radicals {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.unpaired_electrons = Some(unpaired);
                }
            }
            CxEntry::WigglyBonds(wiggly) => {
                for (atom_idx, bond_idx, wedge) in wiggly {
                    if atom_idx as usize >= mol.atoms.len() {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    }
                    let Some(bond) = mol.bonds.get_mut(bond_idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx });
                    };
                    let (a, b) = bond.atoms.as_tuple();
                    if atom_idx != a && atom_idx != b {
                        return Err(ParseError::MismatchedAtomBondIndices { atom_idx, bond_idx });
                    }
                    bond.wedge = Some(wedge);
                }
            }
            CxEntry::CisBonds(indices) => {
                for idx in indices {
                    let Some(bond) = mol.bonds.get_mut(idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: idx });
                    };
                    bond.stereo = Some(BondStereo::Cis);
                }
            }
            CxEntry::TransBonds(indices) => {
                for idx in indices {
                    let Some(bond) = mol.bonds.get_mut(idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: idx });
                    };
                    bond.stereo = Some(BondStereo::Trans);
                }
            }
            CxEntry::UnspecBonds(indices) => {
                for idx in indices {
                    let Some(bond) = mol.bonds.get_mut(idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: idx });
                    };
                    bond.stereo = Some(BondStereo::Either);
                }
            }
            CxEntry::CoordinateBonds(pairs) => {
                for (first_atom, bond_idx) in pairs {
                    if first_atom as usize >= mol.atoms.len() {
                        return Err(ParseError::AtomIndexOutOfBounds {
                            atom_idx: first_atom,
                        });
                    }
                    let Some(bond) = mol.bonds.get_mut(bond_idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx });
                    };
                    let (a, b) = bond.atoms.as_tuple();
                    if first_atom == a {
                        bond.donation = Some(BondDonation::Donating);
                    } else if first_atom == b {
                        bond.donation = Some(BondDonation::Accepting);
                    } else {
                        return Err(ParseError::MismatchedAtomBondIndices {
                            atom_idx: first_atom,
                            bond_idx,
                        });
                    }
                }
            }
            CxEntry::HydrogenBonds(pairs) => {
                for (first_atom, bond_idx) in pairs {
                    if first_atom as usize >= mol.atoms.len() {
                        return Err(ParseError::AtomIndexOutOfBounds {
                            atom_idx: first_atom,
                        });
                    }
                    let Some(bond) = mol.bonds.get_mut(bond_idx as usize) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx });
                    };
                    let (a, b) = bond.atoms.as_tuple();
                    if first_atom != a && first_atom != b {
                        return Err(ParseError::MismatchedAtomBondIndices {
                            atom_idx: first_atom,
                            bond_idx,
                        });
                    }
                    bond.order = BondOrder::Zero;
                    bond.noncovalent = Some(BondNoncovalent::Hydrogen);
                }
            }
            CxEntry::FragmentGroups(groups) => {
                for group in &groups {
                    for &atom_idx in group {
                        if atom_idx as usize >= mol.atoms.len() {
                            return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                        }
                    }
                }
                components = Some(groups);
            }
            CxEntry::StereoGroup(sg) => {
                for &atom_idx in &sg.atoms {
                    if atom_idx as usize >= mol.atoms.len() {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    }
                }

                match sg.group_type {
                    StereoGroupType::Absolute => {
                        // Absolute atoms don't need group storage; stereo_interpretation captures this
                        stereo_interpretation = Some(StereoInterpretation::Absolute);
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
                }
            }
            CxEntry::RelativeStereo => {
                stereo_interpretation = Some(StereoInterpretation::Relative);
            }
            CxEntry::AtomProperties(props) => {
                for (idx, key, value) in props {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.properties.insert(key, value);
                }
            }
        }
    }

    mol.stereo_interpretation = stereo_interpretation;

    // Store CX-specific data if any
    if !stereo_groups.is_empty() || components.is_some() {
        mol.cx_data = Some(CxAnnotationData {
            stereo_groups,
            components,
        });
    }

    Ok(())
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
            Err(ParseError::InvalidCxTag { pos: 0 })
        }
        Err(_) => Err(ParseError::InvalidToken { pos: 0 }),
    }
}

/// Parse a basic CX entry.
///
/// Supported tags:
/// - coordinates
/// - labels / values
/// - radicals
/// - wiggly bonds
/// - cis/trans/unspec (ctu)
/// - coordinate bonds
/// - hydrogen bonds
fn parse_basic_entry(input: &[u8], skip_unknown_cx_tags: bool) -> IResult<&[u8], Option<CxEntry>> {
    if input.is_empty() {
        return Err(Err::Error(NomError::new(input, ErrorKind::Tag)));
    }
    if input.first() == Some(&b'|') {
        // End of CX block.
        return Err(Err::Error(NomError::new(input, ErrorKind::Tag)));
    }

    match alt((
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
        Ok((rest, entry)) => Ok((rest, Some(entry))),
        Err(Err::Error(_)) => parse_unknown_entry(input, skip_unknown_cx_tags),
        Err(e) => Err(e),
    }
}

/// Parse an extended CX entry.
///
/// Supported tags:
/// - all basic tags, plus
/// - fragment groups
/// - enhanced stereo groups
/// - relative stereo tag
/// - atom properties
fn parse_extended_entry(
    input: &[u8],
    skip_unknown_cx_tags: bool,
) -> IResult<&[u8], Option<CxEntry>> {
    if input.is_empty() {
        return Err(Err::Error(NomError::new(input, ErrorKind::Tag)));
    }
    if input.first() == Some(&b'|') {
        // End of CX block.
        return Err(Err::Error(NomError::new(input, ErrorKind::Tag)));
    }

    match alt((
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
        Ok((rest, entry)) => Ok((rest, Some(entry))),
        Err(Err::Error(_)) => parse_unknown_entry(input, skip_unknown_cx_tags),
        Err(e) => Err(e),
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
    let (input, code) = delimited(char('^'), one_of("1234567"), char(':')).parse(input)?;
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

    let result: Vec<_> = pairs
        .into_iter()
        .map(|(atom_idx, bond_idx)| (atom_idx, bond_idx, wedge_type))
        .collect();
    Ok((input, CxEntry::WigglyBonds(result)))
}

/// Parse cis/trans bond annotations: `c:`, `t:`, `ctu:`.
fn parse_cis_trans(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, kind) = alt((
        value('u', tag("ctu:")),
        value('c', tag("c:")),
        value('t', tag("t:")),
    ))
    .parse(input)?;

    let (input, indices) = separated_list0(comma_not_before_entry, nom_u32).parse(input)?;

    match kind {
        'c' => Ok((input, CxEntry::CisBonds(indices))),
        't' => Ok((input, CxEntry::TransBonds(indices))),
        'u' => Ok((input, CxEntry::UnspecBonds(indices))),
        _ => unreachable!("unknown cis/trans/ctu tag"),
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

/// Parse relative stereo tag: `r`.
fn parse_relative_stereo(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (rest, _) = char('r').parse(input)?;

    // Only accept `r` as a standalone tag (`r` or `r:...`), not as the start of
    // another tag (e.g. `rb:`).
    if !rest.is_empty() && !matches!(rest[0], b':' | b',' | b'|') {
        return Err(Err::Error(NomError::new(input, ErrorKind::Tag)));
    }

    // `r:...` is used for reaction/multicomponent cases to list the fragment indices with relative
    // configuration. This isn't meaningful in our (molecule) TableIR parsing, so reject it.
    if rest.first() == Some(&b':') {
        return Err(Err::Failure(NomError::new(input, ErrorKind::Verify)));
    }

    Ok((rest, CxEntry::RelativeStereo))
}

/// Check whether a character can start a CX entry/tag.
fn is_cx_tag_start(c: char) -> bool {
    c.is_ascii_alphabetic() || matches!(c, '(' | '$' | '^' | '&')
}

/// Parse comma only if not followed by an entry-start character.
fn comma_not_before_entry(input: &[u8]) -> IResult<&[u8], char> {
    terminated(char(','), not(satisfy(is_cx_tag_start))).parse(input)
}

/// Skip over an unknown/unrecognized CX entry.
///
/// CXSMILES uses commas as both list separators *within* an entry and as entry separators.
/// We stop at:
/// - the closing `|`, or
/// - a comma that is followed by the start of another CX entry.
fn skip_unknown_entry(input: &[u8]) -> IResult<&[u8], ()> {
    if input.is_empty() {
        return Ok((input, ()));
    }

    let mut i = 0usize;
    while i < input.len() {
        match input[i] {
            b',' => {
                // A comma starts a new entry iff the next non-whitespace char looks like an entry start.
                let mut j = i + 1;
                while j < input.len() && input[j].is_ascii_whitespace() {
                    j += 1;
                }
                if j < input.len() {
                    let next = input[j] as char;
                    if is_cx_tag_start(next) {
                        break;
                    }
                }
            }
            _ => {}
        }
        i += 1;
    }

    if i == 0 {
        return Err(Err::Error(NomError::new(input, ErrorKind::TakeTill1)));
    }
    Ok((&input[i..], ()))
}

/// Parse a basic CX entry (rejects extended features).
fn parse_unknown_entry(
    input: &[u8],
    skip_unknown_cx_tags: bool,
) -> IResult<&[u8], Option<CxEntry>> {
    if skip_unknown_cx_tags {
        let (rest, _) = skip_unknown_entry(input)?;
        Ok((rest, None))
    } else {
        Err(Err::Failure(NomError::new(input, ErrorKind::Verify)))
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::Element;

    use super::*;
    use crate::table_ir::{Atom, Bond, Chirality, ExtendedAtom, ExtendedBond};

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
    #[case::wiggly_bonds(b"|w:0.1,2.3|", vec![CxEntry::WigglyBonds(vec![(0, 1, BondWedge::Either), (2, 3, BondWedge::Either)])])]
    #[case::cis_trans(b"|c:0,1|", vec![CxEntry::CisBonds(vec![0, 1])])]
    #[case::trans_trans(b"|t:0,1|", vec![CxEntry::TransBonds(vec![0, 1])])]
    #[case::atom_labels(b"|$label1;label2;label3$|", vec![CxEntry::Labels(vec![(0, "label1".to_string()), (1, "label2".to_string()), (2, "label3".to_string())])])]
    #[case::atom_values(b"$_AV:value1;value2;value3$|", vec![CxEntry::Values(vec![(0, "value1".to_string()), (1, "value2".to_string()), (2, "value3".to_string())])])]
    #[case::coordinates_2d(b"|(1.5,2.5;3.5,4.5)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.5, 2.5, 0.0), Point3D::new(3.5, 4.5, 0.0)])])]
    #[case::coordinates_3d(b"|(1,2,3;4,5,6)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0), Point3D::new(4.0, 5.0, 6.0)])])]
    #[case::combined_entries(b"|^1:0,1,(1.0,2.0;3.0,4.0),C:2.3|", vec![CxEntry::Radicals(vec![(0, UnpairedElectrons { count: 1, multiplicity: None }),
        (1, UnpairedElectrons { count: 1, multiplicity: None })]), CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 0.0), Point3D::new(3.0, 4.0, 0.0)]), CxEntry::CoordinateBonds(vec![(2, 3)])])]
    fn test_parse_cx_annotations(#[case] input: &[u8], #[case] expected: Vec<CxEntry>) {
        let result = parse_cx_annotations(input, SmilesParseFlags::default());
        let input_str = input.to_str_lossy();
        assert!(result.is_ok(), "{:?} should have succeeded: {:?}", input_str, result);
        let entries = result.unwrap();
        assert_eq!(entries, expected, "{:?} should have parsed to {:?}", input_str, entries);
    }

    #[rstest]
    #[case::unknown_and_known_tag(b"|xyz:123,C:0.1|", ParseError::InvalidCxTag { pos: 0 })]
    #[case::unknown_tag(b"|unknown|", ParseError::InvalidCxTag { pos: 0 })]
    #[case::extended_feature_f(b"|f:0.1|", ParseError::InvalidCxTag { pos: 0 })]
    #[case::extended_feature_a(b"|a:0,1|", ParseError::InvalidCxTag { pos: 0 })]
    #[case::extended_feature_o(b"|o1:0,1|", ParseError::InvalidCxTag { pos: 0 })]
    #[case::extended_feature_and(b"|&1:0,1|", ParseError::InvalidCxTag { pos: 0 })]
    #[case::extended_feature_r(b"|r|", ParseError::InvalidCxTag { pos: 0 })]
    #[case::extended_feature_atomprop(b"|atomProp:0.key.value|", ParseError::InvalidCxTag { pos: 0 })]
    fn test_parse_cx_annotations_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let result = parse_cx_annotations(input, SmilesParseFlags::default());
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

    #[rstest]
    #[case::unknown_and_known_tag(b"|xyz:123,C:0.1|", vec![CxEntry::CoordinateBonds(vec![(0, 1)])])]
    #[case::unknown_tag(b"|unknown|", vec![])]
    fn test_parse_cx_annotations_lenient(#[case] input: &[u8], #[case] expected: Vec<CxEntry>) {
        let flags = SmilesParseFlags::LENIENT;
        let result = parse_cx_annotations(input, flags);
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
    #[case::empty(b"||", vec![])]
    #[case::coordinate_bond(b"|C:0.1|", vec![CxEntry::CoordinateBonds(vec![(0, 1)])])]
    #[case::hydrogen_bond(b"|H:1.2|", vec![CxEntry::HydrogenBonds(vec![(1, 2)])])]
    #[case::radicals(b"|^1:0|", vec![CxEntry::Radicals(vec![(0, UnpairedElectrons { count: 1, multiplicity: None })])])]
    #[case::wiggly_bonds(b"|w:0.1|", vec![CxEntry::WigglyBonds(vec![(0, 1, BondWedge::Either)])])]
    #[case::cis_bonds(b"|c:0|", vec![CxEntry::CisBonds(vec![0])])]
    #[case::trans_bonds(b"|t:0|", vec![CxEntry::TransBonds(vec![0])])]
    #[case::atom_labels(b"|$label$|", vec![CxEntry::Labels(vec![(0, "label".to_string())])])]
    #[case::fragment_groups(b"|f:0.1.2,3.4|", vec![CxEntry::FragmentGroups(vec![vec![0, 1, 2], vec![3, 4]])])]
    #[case::stereo_absolute(b"|a:0,1,2|", vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Absolute, atoms: vec![0, 1, 2] })])]
    #[case::stereo_or(b"|o1:0,1|", vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Or(1), atoms: vec![0, 1] })])]
    #[case::stereo_and(b"|&1:0,1|", vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::And(1), atoms: vec![0, 1] })])]
    #[case::relative_stereo(b"|r|", vec![CxEntry::RelativeStereo])]
    #[case::atom_properties(b"|atomProp:0.key.value|", vec![CxEntry::AtomProperties(vec![(0, "key".to_string(), "value".to_string())])])]
    #[case::coordinates_2d(b"|(1.5,2.5;3.5,4.5)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.5, 2.5, 0.0), Point3D::new(3.5, 4.5, 0.0)])])]
    #[case::coordinates_3d(b"|(1,2,3;4,5,6)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0), Point3D::new(4.0, 5.0, 6.0)])])]
    fn test_parse_extended_cx_annotations(#[case] input: &[u8], #[case] expected: Vec<CxEntry>) {
        let result = parse_extended_cx_annotations(input, SmilesParseFlags::default());
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

    #[rstest]
    #[case::unknown_and_known_tag(b"|xyz:123,C:0.1|", ParseError::InvalidCxTag { pos: 0 })]
    #[case::unknown_tag(b"|unknown|", ParseError::InvalidCxTag { pos: 0 })]
    #[case::relative_stereo_with_fragment_list(b"|r:0|", ParseError::InvalidCxTag { pos: 0 })]
    fn test_parse_extended_cx_annotations_invalid(
        #[case] input: &[u8],
        #[case] expected: ParseError,
    ) {
        let result = parse_extended_cx_annotations(input, SmilesParseFlags::default());
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

    #[rstest]
    #[case::unknown_and_known_tag(b"|xyz:123,C:0.1|", vec![CxEntry::CoordinateBonds(vec![(0, 1)])])]
    #[case::unknown_tag(b"|unknown|", vec![])]
    fn test_parse_extended_cx_annotations_lenient(
        #[case] input: &[u8],
        #[case] expected: Vec<CxEntry>,
    ) {
        let flags = SmilesParseFlags::LENIENT;
        let result = parse_extended_cx_annotations(input, flags);
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

    #[rstest]
    #[case::relative_stereo_with_fragment_list(b"|r:0|", SmilesParseFlags::LENIENT)]
    fn test_parse_extended_cx_annotations_lenient_invalid(
        #[case] input: &[u8],
        #[case] flags: SmilesParseFlags,
    ) {
        let result = parse_extended_cx_annotations(input, flags);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_err(),
            "{:?} should have failed: {:?}",
            input_str,
            result
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::coordinates(vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0)])], |mol: &Molecule| mol.positions == Some(vec![Point3D::new(1.0, 2.0, 3.0)]))]
    #[case::labels(vec![CxEntry::Labels(vec![(0, "C1".to_string()), (1, "N1".to_string())])], |mol: &Molecule| mol.atoms[0].label == Some("C1".to_string()) && mol.atoms[1].label == Some("N1".to_string()))]
    #[case::values(vec![CxEntry::Values(vec![(0, "val0".to_string())])], |mol: &Molecule| mol.atoms[0].value == Some("val0".to_string()))]
    #[case::radicals(vec![CxEntry::Radicals(vec![(0, UnpairedElectrons { count: 1, multiplicity: None })])],
        |mol: &Molecule| mol.atoms[0].unpaired_electrons == Some(UnpairedElectrons { count: 1, multiplicity: None }))]
    #[case::wiggly_bonds(vec![CxEntry::WigglyBonds(vec![(0, 0, BondWedge::Either)])], |mol: &Molecule| mol.bonds[0].wedge == Some(BondWedge::Either))]
    #[case::cis_bonds(vec![CxEntry::CisBonds(vec![0])], |mol: &Molecule| mol.bonds[0].stereo == Some(BondStereo::Cis))]
    #[case::trans_bonds(vec![CxEntry::TransBonds(vec![1])], |mol: &Molecule| mol.bonds[1].stereo == Some(BondStereo::Trans))]
    #[case::coordinate_bonds(vec![CxEntry::CoordinateBonds(vec![(0, 0)])], |mol: &Molecule| mol.bonds[0].donation == Some(BondDonation::Donating))]
    #[case::hydrogen_bonds(vec![CxEntry::HydrogenBonds(vec![(0, 0)])], |mol: &Molecule| mol.bonds[0].noncovalent == Some(BondNoncovalent::Hydrogen) && mol.bonds[0].order == BondOrder::Zero)]
    #[case::extended_entries(vec![CxEntry::FragmentGroups(vec![vec![0, 1]]), CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Absolute, atoms: vec![0] }),
        CxEntry::RelativeStereo, CxEntry::AtomProperties(vec![(0, "k".to_string(), "v".to_string())])], |mol: &Molecule| mol.atoms[0].label.is_none())]
    fn test_update_molecule(
        triatomic_molecule: Molecule,
        #[case] entries: Vec<CxEntry>,
        #[case] check: fn(&Molecule) -> bool,
    ) {
        let mut mol = triatomic_molecule;
        update_molecule(&mut mol, entries).unwrap();
        assert!(check(&mol));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::coordinates(vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0)])], |mol: &ExtendedMolecule| mol.positions == Some(vec![Point3D::new(1.0, 2.0, 3.0)]))]
    #[case::labels(vec![CxEntry::Labels(vec![(0, "C1".to_string())])], |mol: &ExtendedMolecule| mol.atoms[0].label == Some("C1".to_string()))]
    #[case::values(vec![CxEntry::Values(vec![(1, "val1".to_string())])], |mol: &ExtendedMolecule| mol.atoms[1].value == Some("val1".to_string()))]
    #[case::radicals(vec![CxEntry::Radicals(vec![(2, UnpairedElectrons { count: 2, multiplicity: None })])],
        |mol: &ExtendedMolecule| mol.atoms[2].unpaired_electrons == Some(UnpairedElectrons { count: 2, multiplicity: None }))]
    #[case::wiggly_bonds(vec![CxEntry::WigglyBonds(vec![(1, 0, BondWedge::Either)])], |mol: &ExtendedMolecule| mol.bonds[0].wedge == Some(BondWedge::Either))]
    #[case::cis_bonds(vec![CxEntry::CisBonds(vec![1])], |mol: &ExtendedMolecule| mol.bonds[1].stereo == Some(BondStereo::Cis))]
    #[case::trans_bonds(vec![CxEntry::TransBonds(vec![0])], |mol: &ExtendedMolecule| mol.bonds[0].stereo == Some(BondStereo::Trans))]
    #[case::coordinate_bonds(vec![CxEntry::CoordinateBonds(vec![(1, 0)])], |mol: &ExtendedMolecule| mol.bonds[0].donation == Some(BondDonation::Accepting))]
    #[case::hydrogen_bonds(vec![CxEntry::HydrogenBonds(vec![(0, 0)])], |mol: &ExtendedMolecule| mol.bonds[0].noncovalent == Some(BondNoncovalent::Hydrogen) && mol.bonds[0].order == BondOrder::Zero)]
    #[case::fragment_groups(vec![CxEntry::FragmentGroups(vec![vec![0, 1], vec![2]])], |mol: &ExtendedMolecule| mol.cx_data.as_ref().map(|d| d.components.as_ref()) == Some(Some(&vec![vec![0, 1], vec![2]])))]
    #[case::stereo_group_absolute(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Absolute, atoms: vec![0, 1] })],
        |mol: &ExtendedMolecule| mol.stereo_interpretation == Some(StereoInterpretation::Absolute) && mol.cx_data.is_none())]
    #[case::stereo_group_or(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Or(1), atoms: vec![0] })],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_groups.get(&1)) == Some(&StereoSet { atoms: vec![0], mode: StereoSetMode::Correlated }))]
    #[case::stereo_group_and(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::And(2), atoms: vec![1] })],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_groups.get(&2)) == Some(&StereoSet { atoms: vec![1], mode: StereoSetMode::Independent }))]
    #[case::relative_stereo(vec![CxEntry::RelativeStereo],
        |mol: &ExtendedMolecule| mol.stereo_interpretation == Some(StereoInterpretation::Relative) && mol.cx_data.is_none())]
    #[case::atom_properties(vec![CxEntry::AtomProperties(vec![(0, "key".to_string(), "value".to_string())])], |mol: &ExtendedMolecule| mol.atoms[0].properties.get("key") == Some(&"value".to_string()))]
    fn test_update_extended_molecule(
        triatomic_extended_molecule: ExtendedMolecule,
        #[case] entries: Vec<CxEntry>,
        #[case] check: fn(&ExtendedMolecule) -> bool,
    ) {
        let mut mol = triatomic_extended_molecule;
        update_extended_molecule(&mut mol, entries).unwrap();
        assert!(check(&mol));
    }
}
