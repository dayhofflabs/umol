//! CXSMILES annotation block parser
//!
//! Parses the `|...|` extension block in CXSMILES format.
//! Two parsers are provided:
//! - `parse_cx_annotations`: basic annotations only (for Molecule)
//! - `parse_extended_cx_annotations`: all annotations (for ExtendedMolecule)

use std::collections::BTreeMap;

use bstr::ByteSlice;
use nom::branch::alt;
use nom::bytes::complete::{tag, take_while, take_while1};
use nom::character::complete::{char, one_of, satisfy, u32 as nom_u32};
use nom::combinator::{not, opt, peek, value};
use nom::error::{Error as NomError, ErrorKind};
use nom::multi::{count, many1, separated_list0};
use nom::number::complete::double;
use nom::sequence::{delimited, preceded, separated_pair, terminated};
use nom::{Err, IResult, Parser};
use umol_data::SpinMultiplicity;

use super::super::config::SmilesParseFlags;
use super::super::error::ParseError;
use super::utils::{split_escaped_semicolons, unescape_html_entities};
use crate::bond::BondNoncovalent;
use crate::position::Point3D;
use crate::table_ir::atom::{BicycloStereo, BicycloStereoData};
use crate::table_ir::{
    BondDonation, BondOrder, BondStereo, BondWedge, CxAnnotationData, ExtendedMolecule,
    ExtendedReaction, LinkAtom, Molecule, MulticenterBond, MulticenterSet, Reaction, RingBondCount,
    SGroup, SGroupBracketCoords, SGroupBracketOrientation, SGroupBracketStyle, SGroupConnectivity,
    SGroupData, SGroupDataType, SGroupSubtype, SGroupType, StereoInterpretation, StereoSet,
    StereoSetMode, SubstitutionCount, UnpairedElectrons, UnsaturatedAtom,
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
    /// Lone pairs: LP:idx,idx,... or lp:idx:count,...
    LonePairs(Vec<(u32, u8)>),
    /// Multicenter bonds: m:central:ligand.ligand,...
    MulticenterBonds(Vec<(u32, Vec<u32>)>),
    /// Ring bond count query: rb:idx:value,...
    RingBondCount(Vec<(u32, RingBondCount)>),
    /// Substitution count query: s:idx:value,...
    SubstitutionCount(Vec<(u32, SubstitutionCount)>),
    /// Unsaturated atom query: u:idx,idx,...
    Unsaturated(Vec<u32>),
    /// Ligand order: LO:centerIdx:idx1.idx2.idx3,...
    LigandOrder(Vec<(u32, Vec<u32>)>),
    /// Link nodes: LN:atom:min.max or LN:atom:min.max.outer1.outer2
    LinkNodes(Vec<(u32, LinkAtom)>),
    /// Polymer S-group: Sg:type:subtype:atoms:subscript:...
    Sgroup(SGroup),
    /// Data S-group: SgD:atomIndices:name:data:queryOp:unit:tag:coords
    SgroupData(SGroup),
    /// S-group hierarchy: SgH:parent:child.child,...
    SgroupHierarchy(Vec<(u32, Vec<u32>)>),
    /// Bicyclo stereo: THB:/TLB:/TEB:ligand:connection:lower:higher (one tag can have multiple entries)
    BicycloStereo(Vec<BicycloStereo>),
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
/// TODO: Add FragmentGroups, StereoGroups, RelativeStereo, LigandOrder to Molecule?
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
            CxEntry::LonePairs(pairs) => {
                for (idx, count) in pairs {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.lone_pairs = Some(count);
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
            CxEntry::MulticenterBonds(bonds) => {
                for (center, ligands) in bonds {
                    if center as usize >= mol.atoms.len() {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: center });
                    }
                    for &ligand in &ligands {
                        if ligand as usize >= mol.atoms.len() {
                            return Err(ParseError::AtomIndexOutOfBounds { atom_idx: ligand });
                        }
                    }
                    let bond = MulticenterBond::new(vec![
                        MulticenterSet::single(center),
                        MulticenterSet::new(ligands),
                    ]);
                    mol.multicenter_bonds.push(bond);
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
    let mut stereo_groups: BTreeMap<u32, StereoSet> = BTreeMap::new();
    let mut components: Option<Vec<Vec<u32>>> = None;
    let mut sgroups: BTreeMap<u32, SGroup> = BTreeMap::new();
    let mut sgroup_index: u32 = 0;
    let mut bicyclo_stereo: Vec<BicycloStereo> = vec![];

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
            CxEntry::LonePairs(pairs) => {
                for (idx, count) in pairs {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.lone_pairs = Some(count);
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
            CxEntry::MulticenterBonds(bonds) => {
                for (center, ligands) in bonds {
                    if center as usize >= mol.atoms.len() {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: center });
                    }
                    for &ligand in &ligands {
                        if ligand as usize >= mol.atoms.len() {
                            return Err(ParseError::AtomIndexOutOfBounds { atom_idx: ligand });
                        }
                    }
                    let bond = MulticenterBond::new(vec![
                        MulticenterSet::single(center),
                        MulticenterSet::new(ligands),
                    ]);
                    mol.multicenter_bonds.push(bond);
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
            CxEntry::RingBondCount(entries) => {
                for (idx, rbc) in entries {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.ring_bond_count = Some(rbc);
                }
            }
            CxEntry::SubstitutionCount(entries) => {
                for (idx, sc) in entries {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.substitution_count = Some(sc);
                }
            }
            CxEntry::Unsaturated(indices) => {
                for idx in indices {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.unsaturated = Some(UnsaturatedAtom);
                }
            }
            CxEntry::LinkNodes(entries) => {
                for (idx, link) in entries {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.link_atom = Some(link);
                }
            }
            CxEntry::LigandOrder(entries) => {
                let atom_count = mol.atoms.len();
                for (center_idx, neighbors) in entries {
                    for &neighbor_idx in &neighbors {
                        if neighbor_idx as usize >= atom_count {
                            return Err(ParseError::AtomIndexOutOfBounds {
                                atom_idx: neighbor_idx,
                            });
                        }
                    }
                    let Some(atom) = mol.atoms.get_mut(center_idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds {
                            atom_idx: center_idx,
                        });
                    };
                    let order: Vec<_> = neighbors
                        .into_iter()
                        .enumerate()
                        .map(|(i, n)| (n, (i + 1) as u8))
                        .collect();
                    atom.ligand_order = Some(order);
                }
            }
            CxEntry::Sgroup(sgroup) => {
                let atom_count = mol.atoms.len();
                for &idx in &sgroup.atom_indices {
                    if idx as usize >= atom_count {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    }
                }
                sgroups.insert(sgroup_index, sgroup);
                sgroup_index += 1;
            }
            CxEntry::SgroupData(sgroup) => {
                let atom_count = mol.atoms.len();
                for &idx in &sgroup.atom_indices {
                    if idx as usize >= atom_count {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    }
                }
                sgroups.insert(sgroup_index, sgroup);
                sgroup_index += 1;
            }
            CxEntry::SgroupHierarchy(pairs) => {
                for (parent_idx, children) in &pairs {
                    if !sgroups.contains_key(parent_idx) {
                        return Err(ParseError::SgroupIndexOutOfBounds {
                            sgroup_idx: *parent_idx,
                        });
                    }
                    for &child_idx in children {
                        if let Some(sg) = sgroups.get_mut(&child_idx) {
                            sg.hierarchy_parent = Some(*parent_idx);
                        } else {
                            return Err(ParseError::SgroupIndexOutOfBounds {
                                sgroup_idx: child_idx,
                            });
                        }
                    }
                }
            }
            CxEntry::BicycloStereo(entries) => {
                bicyclo_stereo.extend(entries);
            }
        }
    }

    mol.stereo_interpretation = stereo_interpretation;

    // Store CX-specific data if any
    if !stereo_groups.is_empty()
        || components.is_some()
        || !sgroups.is_empty()
        || !bicyclo_stereo.is_empty()
    {
        mol.cx_data = Some(CxAnnotationData {
            stereo_groups,
            components,
            sgroups,
            bicyclo_stereo: if bicyclo_stereo.is_empty() {
                None
            } else {
                Some(bicyclo_stereo)
            },
            ..Default::default()
        });
    }

    Ok(())
}

pub fn update_reaction(
    reaction: &mut Reaction,
    split: (Vec<CxEntry>, Vec<CxEntry>, Vec<CxEntry>),
) -> Result<(), ParseError> {
    let (reactant_entries, agent_entries, product_entries) = split;
    if !reactant_entries.is_empty() {
        update_molecule(&mut reaction.reactants, reactant_entries)?;
    }
    if !agent_entries.is_empty() {
        update_molecule(&mut reaction.agents, agent_entries)?;
    }
    if !product_entries.is_empty() {
        update_molecule(&mut reaction.products, product_entries)?;
    }
    Ok(())
}

pub fn update_extended_reaction(
    reaction: &mut ExtendedReaction,
    split: (Vec<CxEntry>, Vec<CxEntry>, Vec<CxEntry>),
) -> Result<(), ParseError> {
    let (reactant_entries, agent_entries, product_entries) = split;
    if !reactant_entries.is_empty() {
        update_extended_molecule(&mut reaction.reactants, reactant_entries)?;
    }
    if !agent_entries.is_empty() {
        update_extended_molecule(&mut reaction.agents, agent_entries)?;
    }
    if !product_entries.is_empty() {
        update_extended_molecule(&mut reaction.products, product_entries)?;
    }
    Ok(())
}

#[allow(clippy::type_complexity)]
pub fn split_reaction_cx_entries(
    entries: Vec<CxEntry>,
    reactant_atom_count: usize,
    reactant_bond_count: usize,
    agent_atom_count: usize,
    agent_bond_count: usize,
    product_atom_count: usize,
    product_bond_count: usize,
) -> Result<(Vec<CxEntry>, Vec<CxEntry>, Vec<CxEntry>), ParseError> {
    let atom_starts = [
        0u32,
        reactant_atom_count as u32,
        (reactant_atom_count + agent_atom_count) as u32,
    ];
    let atom_ends = [
        reactant_atom_count as u32,
        (reactant_atom_count + agent_atom_count) as u32,
        (reactant_atom_count + agent_atom_count + product_atom_count) as u32,
    ];
    let bond_starts = [
        0u32,
        reactant_bond_count as u32,
        (reactant_bond_count + agent_bond_count) as u32,
    ];
    let bond_ends = [
        reactant_bond_count as u32,
        (reactant_bond_count + agent_bond_count) as u32,
        (reactant_bond_count + agent_bond_count + product_bond_count) as u32,
    ];

    let atom_side = |idx: u32| -> Option<(usize, u32)> {
        for side in 0..3 {
            if idx >= atom_starts[side] && idx < atom_ends[side] {
                return Some((side, idx - atom_starts[side]));
            }
        }
        None
    };
    let bond_side = |idx: u32| -> Option<(usize, u32)> {
        for side in 0..3 {
            if idx >= bond_starts[side] && idx < bond_ends[side] {
                return Some((side, idx - bond_starts[side]));
            }
        }
        None
    };

    let mut out = (Vec::new(), Vec::new(), Vec::new());
    let mut push_for_side = |side: usize, entry: CxEntry| match side {
        0 => out.0.push(entry),
        1 => out.1.push(entry),
        _ => out.2.push(entry),
    };

    for entry in entries {
        match entry {
            CxEntry::Coordinates(coords) => {
                for side in 0..3 {
                    let start = atom_starts[side] as usize;
                    let end = atom_ends[side] as usize;
                    if start < coords.len() {
                        let side_coords = coords[start..coords.len().min(end)].to_vec();
                        if !side_coords.is_empty() {
                            push_for_side(side, CxEntry::Coordinates(side_coords));
                        }
                    }
                }
            }
            CxEntry::Labels(pairs) => {
                let mut side_pairs = [Vec::new(), Vec::new(), Vec::new()];
                for (idx, value) in pairs {
                    let Some((side, local_idx)) = atom_side(idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    side_pairs[side].push((local_idx, value));
                }
                for (side, pairs) in side_pairs.into_iter().enumerate() {
                    if !pairs.is_empty() {
                        push_for_side(side, CxEntry::Labels(pairs));
                    }
                }
            }
            CxEntry::Values(pairs) => {
                let mut side_pairs = [Vec::new(), Vec::new(), Vec::new()];
                for (idx, value) in pairs {
                    let Some((side, local_idx)) = atom_side(idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    side_pairs[side].push((local_idx, value));
                }
                for (side, pairs) in side_pairs.into_iter().enumerate() {
                    if !pairs.is_empty() {
                        push_for_side(side, CxEntry::Values(pairs));
                    }
                }
            }
            CxEntry::Radicals(pairs) => {
                let mut side_pairs = [Vec::new(), Vec::new(), Vec::new()];
                for (idx, value) in pairs {
                    let Some((side, local_idx)) = atom_side(idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    side_pairs[side].push((local_idx, value));
                }
                for (side, pairs) in side_pairs.into_iter().enumerate() {
                    if !pairs.is_empty() {
                        push_for_side(side, CxEntry::Radicals(pairs));
                    }
                }
            }
            CxEntry::WigglyBonds(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (atom_idx, bond_idx, wedge) in items {
                    let Some((atom_side_idx, local_atom)) = atom_side(atom_idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    };
                    let Some((bond_side_idx, local_bond)) = bond_side(bond_idx) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx });
                    };
                    if atom_side_idx != bond_side_idx {
                        return Err(ParseError::MismatchedAtomBondIndices { atom_idx, bond_idx });
                    }
                    side_items[atom_side_idx].push((local_atom, local_bond, wedge));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::WigglyBonds(items));
                    }
                }
            }
            CxEntry::CisBonds(indices) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for idx in indices {
                    let Some((side, local_idx)) = bond_side(idx) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: idx });
                    };
                    side_items[side].push(local_idx);
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::CisBonds(items));
                    }
                }
            }
            CxEntry::TransBonds(indices) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for idx in indices {
                    let Some((side, local_idx)) = bond_side(idx) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: idx });
                    };
                    side_items[side].push(local_idx);
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::TransBonds(items));
                    }
                }
            }
            CxEntry::UnspecBonds(indices) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for idx in indices {
                    let Some((side, local_idx)) = bond_side(idx) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: idx });
                    };
                    side_items[side].push(local_idx);
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::UnspecBonds(items));
                    }
                }
            }
            CxEntry::CoordinateBonds(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (atom_idx, bond_idx) in items {
                    let Some((atom_side_idx, local_atom)) = atom_side(atom_idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    };
                    let Some((bond_side_idx, local_bond)) = bond_side(bond_idx) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx });
                    };
                    if atom_side_idx != bond_side_idx {
                        return Err(ParseError::MismatchedAtomBondIndices { atom_idx, bond_idx });
                    }
                    side_items[atom_side_idx].push((local_atom, local_bond));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::CoordinateBonds(items));
                    }
                }
            }
            CxEntry::HydrogenBonds(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (atom_idx, bond_idx) in items {
                    let Some((atom_side_idx, local_atom)) = atom_side(atom_idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    };
                    let Some((bond_side_idx, local_bond)) = bond_side(bond_idx) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx });
                    };
                    if atom_side_idx != bond_side_idx {
                        return Err(ParseError::MismatchedAtomBondIndices { atom_idx, bond_idx });
                    }
                    side_items[atom_side_idx].push((local_atom, local_bond));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::HydrogenBonds(items));
                    }
                }
            }
            CxEntry::LonePairs(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (atom_idx, count) in items {
                    let Some((side, local_idx)) = atom_side(atom_idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    };
                    side_items[side].push((local_idx, count));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::LonePairs(items));
                    }
                }
            }
            CxEntry::MulticenterBonds(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (center, ligands) in items {
                    let Some((side, local_center)) = atom_side(center) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: center });
                    };
                    let mut local_ligands = Vec::with_capacity(ligands.len());
                    for ligand in ligands {
                        let Some((lig_side, local_lig)) = atom_side(ligand) else {
                            return Err(ParseError::AtomIndexOutOfBounds { atom_idx: ligand });
                        };
                        if lig_side != side {
                            return Err(ParseError::InvalidCxTag { pos: 0 });
                        }
                        local_ligands.push(local_lig);
                    }
                    side_items[side].push((local_center, local_ligands));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::MulticenterBonds(items));
                    }
                }
            }
            CxEntry::FragmentGroups(groups) => {
                let mut side_groups = [Vec::new(), Vec::new(), Vec::new()];
                for group in groups {
                    let mut local = Vec::with_capacity(group.len());
                    let mut side_opt = None;
                    for atom_idx in group {
                        let Some((side, local_idx)) = atom_side(atom_idx) else {
                            return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                        };
                        if let Some(existing) = side_opt {
                            if existing != side {
                                return Err(ParseError::InvalidCxTag { pos: 0 });
                            }
                        } else {
                            side_opt = Some(side);
                        }
                        local.push(local_idx);
                    }
                    if let Some(side) = side_opt {
                        side_groups[side].push(local);
                    }
                }
                for (side, groups) in side_groups.into_iter().enumerate() {
                    if !groups.is_empty() {
                        push_for_side(side, CxEntry::FragmentGroups(groups));
                    }
                }
            }
            CxEntry::StereoGroup(mut sg) => {
                let mut side_opt = None;
                for atom in &mut sg.atoms {
                    let Some((side, local_idx)) = atom_side(*atom) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: *atom });
                    };
                    if let Some(existing) = side_opt {
                        if existing != side {
                            return Err(ParseError::InvalidCxTag { pos: 0 });
                        }
                    } else {
                        side_opt = Some(side);
                    }
                    *atom = local_idx;
                }
                if let Some(side) = side_opt {
                    push_for_side(side, CxEntry::StereoGroup(sg));
                }
            }
            CxEntry::RelativeStereo => {
                for side in 0..3 {
                    if atom_ends[side] > atom_starts[side] {
                        push_for_side(side, CxEntry::RelativeStereo);
                    }
                }
            }
            CxEntry::AtomProperties(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (atom_idx, k, v) in items {
                    let Some((side, local_idx)) = atom_side(atom_idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    };
                    side_items[side].push((local_idx, k, v));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::AtomProperties(items));
                    }
                }
            }
            CxEntry::RingBondCount(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (atom_idx, value) in items {
                    let Some((side, local_idx)) = atom_side(atom_idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    };
                    side_items[side].push((local_idx, value));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::RingBondCount(items));
                    }
                }
            }
            CxEntry::SubstitutionCount(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (atom_idx, value) in items {
                    let Some((side, local_idx)) = atom_side(atom_idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    };
                    side_items[side].push((local_idx, value));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::SubstitutionCount(items));
                    }
                }
            }
            CxEntry::Unsaturated(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for atom_idx in items {
                    let Some((side, local_idx)) = atom_side(atom_idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    };
                    side_items[side].push(local_idx);
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::Unsaturated(items));
                    }
                }
            }
            CxEntry::LigandOrder(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (center, neigh) in items {
                    let Some((side, local_center)) = atom_side(center) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: center });
                    };
                    let mut local_neigh = Vec::with_capacity(neigh.len());
                    for atom_idx in neigh {
                        let Some((s, local_idx)) = atom_side(atom_idx) else {
                            return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                        };
                        if s != side {
                            return Err(ParseError::InvalidCxTag { pos: 0 });
                        }
                        local_neigh.push(local_idx);
                    }
                    side_items[side].push((local_center, local_neigh));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::LigandOrder(items));
                    }
                }
            }
            CxEntry::LinkNodes(items) => {
                let mut side_items = [Vec::new(), Vec::new(), Vec::new()];
                for (atom_idx, link) in items {
                    let Some((side, local_idx)) = atom_side(atom_idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx });
                    };
                    side_items[side].push((local_idx, link));
                }
                for (side, items) in side_items.into_iter().enumerate() {
                    if !items.is_empty() {
                        push_for_side(side, CxEntry::LinkNodes(items));
                    }
                }
            }
            CxEntry::Sgroup(mut sg) => {
                let mut side_opt = None;
                for idx in &mut sg.atom_indices {
                    let Some((side, local_idx)) = atom_side(*idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: *idx });
                    };
                    if let Some(existing) = side_opt {
                        if existing != side {
                            return Err(ParseError::InvalidCxTag { pos: 0 });
                        }
                    } else {
                        side_opt = Some(side);
                    }
                    *idx = local_idx;
                }
                for idx in &mut sg.bond_indices {
                    let Some((side, local_idx)) = bond_side(*idx) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: *idx });
                    };
                    if let Some(existing) = side_opt {
                        if existing != side {
                            return Err(ParseError::InvalidCxTag { pos: 0 });
                        }
                    } else {
                        side_opt = Some(side);
                    }
                    *idx = local_idx;
                }
                if let Some(indices) = &mut sg.parent_atom_indices {
                    for idx in indices {
                        let Some((side, local_idx)) = atom_side(*idx) else {
                            return Err(ParseError::AtomIndexOutOfBounds { atom_idx: *idx });
                        };
                        if let Some(existing) = side_opt {
                            if existing != side {
                                return Err(ParseError::InvalidCxTag { pos: 0 });
                            }
                        } else {
                            side_opt = Some(side);
                        }
                        *idx = local_idx;
                    }
                }
                if let Some(indices) = &mut sg.correspondence {
                    for idx in indices {
                        let Some((side, local_idx)) = bond_side(*idx) else {
                            return Err(ParseError::BondIndexOutOfBounds { bond_idx: *idx });
                        };
                        if let Some(existing) = side_opt {
                            if existing != side {
                                return Err(ParseError::InvalidCxTag { pos: 0 });
                            }
                        } else {
                            side_opt = Some(side);
                        }
                        *idx = local_idx;
                    }
                }
                if let Some(cb) = sg.connecting_bond.as_mut() {
                    let Some((side, local_idx)) = bond_side(cb.bond_index) else {
                        return Err(ParseError::BondIndexOutOfBounds {
                            bond_idx: cb.bond_index,
                        });
                    };
                    if let Some(existing) = side_opt {
                        if existing != side {
                            return Err(ParseError::InvalidCxTag { pos: 0 });
                        }
                    } else {
                        side_opt = Some(side);
                    }
                    cb.bond_index = local_idx;
                }
                if let Some(side) = side_opt {
                    push_for_side(side, CxEntry::Sgroup(sg));
                }
            }
            CxEntry::SgroupData(mut sg) => {
                let mut side_opt = None;
                for idx in &mut sg.atom_indices {
                    let Some((side, local_idx)) = atom_side(*idx) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: *idx });
                    };
                    if let Some(existing) = side_opt {
                        if existing != side {
                            return Err(ParseError::InvalidCxTag { pos: 0 });
                        }
                    } else {
                        side_opt = Some(side);
                    }
                    *idx = local_idx;
                }
                for idx in &mut sg.bond_indices {
                    let Some((side, local_idx)) = bond_side(*idx) else {
                        return Err(ParseError::BondIndexOutOfBounds { bond_idx: *idx });
                    };
                    if let Some(existing) = side_opt {
                        if existing != side {
                            return Err(ParseError::InvalidCxTag { pos: 0 });
                        }
                    } else {
                        side_opt = Some(side);
                    }
                    *idx = local_idx;
                }
                if let Some(indices) = &mut sg.parent_atom_indices {
                    for idx in indices {
                        let Some((side, local_idx)) = atom_side(*idx) else {
                            return Err(ParseError::AtomIndexOutOfBounds { atom_idx: *idx });
                        };
                        if let Some(existing) = side_opt {
                            if existing != side {
                                return Err(ParseError::InvalidCxTag { pos: 0 });
                            }
                        } else {
                            side_opt = Some(side);
                        }
                        *idx = local_idx;
                    }
                }
                if let Some(indices) = &mut sg.correspondence {
                    for idx in indices {
                        let Some((side, local_idx)) = bond_side(*idx) else {
                            return Err(ParseError::BondIndexOutOfBounds { bond_idx: *idx });
                        };
                        if let Some(existing) = side_opt {
                            if existing != side {
                                return Err(ParseError::InvalidCxTag { pos: 0 });
                            }
                        } else {
                            side_opt = Some(side);
                        }
                        *idx = local_idx;
                    }
                }
                if let Some(cb) = sg.connecting_bond.as_mut() {
                    let Some((side, local_idx)) = bond_side(cb.bond_index) else {
                        return Err(ParseError::BondIndexOutOfBounds {
                            bond_idx: cb.bond_index,
                        });
                    };
                    if let Some(existing) = side_opt {
                        if existing != side {
                            return Err(ParseError::InvalidCxTag { pos: 0 });
                        }
                    } else {
                        side_opt = Some(side);
                    }
                    cb.bond_index = local_idx;
                }
                if let Some(side) = side_opt {
                    push_for_side(side, CxEntry::SgroupData(sg));
                }
            }
            CxEntry::SgroupHierarchy(pairs) => {
                for side in 0..3 {
                    if atom_ends[side] > atom_starts[side] {
                        push_for_side(side, CxEntry::SgroupHierarchy(pairs.clone()));
                    }
                }
            }
            CxEntry::BicycloStereo(entries) => {
                let mut per_side = [Vec::new(), Vec::new(), Vec::new()];
                for entry in entries {
                    let (lig, con, lo, hi) = match &entry {
                        BicycloStereo::TowardsHigherBridge(d)
                        | BicycloStereo::TowardsLowerBridge(d)
                        | BicycloStereo::TowardsEitherBridge(d) => (
                            d.ligand_atom,
                            d.connection_atom,
                            d.lower_bridge_atoms.clone(),
                            d.higher_bridge_atoms.clone(),
                        ),
                    };
                    let Some((side, lig_local)) = atom_side(lig) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: lig });
                    };
                    let Some((side2, con_local)) = atom_side(con) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: con });
                    };
                    if side != side2 {
                        return Err(ParseError::InvalidCxTag { pos: 0 });
                    }
                    let map_vec = |v: Vec<u32>| -> Result<Vec<u32>, ParseError> {
                        let mut out = Vec::with_capacity(v.len());
                        for a in v {
                            let Some((s, local)) = atom_side(a) else {
                                return Err(ParseError::AtomIndexOutOfBounds { atom_idx: a });
                            };
                            if s != side {
                                return Err(ParseError::InvalidCxTag { pos: 0 });
                            }
                            out.push(local);
                        }
                        Ok(out)
                    };
                    let lower = map_vec(lo)?;
                    let higher = map_vec(hi)?;
                    let data = BicycloStereoData {
                        ligand_atom: lig_local,
                        connection_atom: con_local,
                        lower_bridge_atoms: lower,
                        higher_bridge_atoms: higher,
                    };
                    let remapped = match entry {
                        BicycloStereo::TowardsHigherBridge(_) => {
                            BicycloStereo::TowardsHigherBridge(data)
                        }
                        BicycloStereo::TowardsLowerBridge(_) => {
                            BicycloStereo::TowardsLowerBridge(data)
                        }
                        BicycloStereo::TowardsEitherBridge(_) => {
                            BicycloStereo::TowardsEitherBridge(data)
                        }
                    };
                    per_side[side].push(remapped);
                }
                for (side, entries) in per_side.into_iter().enumerate() {
                    if !entries.is_empty() {
                        push_for_side(side, CxEntry::BicycloStereo(entries));
                    }
                }
            }
        }
    }

    Ok(out)
}

/// Parse comma only when followed by an entry-start character (separator between CX entries).
fn comma_before_entry(input: &[u8]) -> IResult<&[u8], char> {
    terminated(char(','), peek(satisfy(is_cx_tag_start))).parse(input)
}

fn parse_cx_block<'inp>(
    input: &'inp [u8],
    entry_parser: impl Parser<&'inp [u8], Output = Option<CxEntry>, Error = NomError<&'inp [u8]>>,
) -> Result<Vec<CxEntry>, ParseError> {
    match delimited(
        opt(char('|')),
        separated_list0(comma_before_entry, entry_parser),
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
        parse_lone_pairs,
        parse_multicenter,
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
        alt((
            parse_wiggly_bonds,
            parse_cis_trans,
            parse_coordinate_bonds,
            parse_hydrogen_bonds,
        )),
        alt((parse_lone_pairs, parse_multicenter)),
        parse_fragment_groups,
        alt((
            parse_stereo_absolute,
            parse_stereo_or_and,
            parse_relative_stereo,
        )),
        parse_atom_properties,
        alt((
            parse_ring_bond_count,
            parse_substitution_count,
            parse_unsaturated,
        )),
        alt((parse_ligand_order, parse_link_nodes)),
        alt((parse_sgroup_data, parse_sgroup, parse_sgroup_hierarchy)),
        parse_bicyclo_stereo,
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

/// Parse lone pairs: `LP:idx,idx,...` (unspecified count) or `lp:idx:count,...` (explicit count).
fn parse_lone_pairs(input: &[u8]) -> IResult<&[u8], CxEntry> {
    // Try lp: format first (has explicit counts)
    if let Ok((input, _)) = tag::<_, _, NomError<&[u8]>>("lp:")(input) {
        let parse_entry = separated_pair(nom_u32, char(':'), nom_u32);
        let (input, entries) = separated_list0(comma_not_before_entry, parse_entry).parse(input)?;
        let result: Vec<_> = entries
            .into_iter()
            .map(|(idx, count)| (idx, count as u8))
            .collect();
        return Ok((input, CxEntry::LonePairs(result)));
    }

    // LP: format (implicit count, treated as 1 per atom)
    let (input, indices) =
        preceded(tag("LP:"), separated_list0(comma_not_before_entry, nom_u32)).parse(input)?;
    let result: Vec<_> = indices.into_iter().map(|idx| (idx, 1u8)).collect();
    Ok((input, CxEntry::LonePairs(result)))
}

/// Parse multicenter bonds: `m:central:ligand.ligand,...`
fn parse_multicenter(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, _) = tag("m:").parse(input)?;
    let parse_entry = |i| {
        (
            terminated(nom_u32, char(':')),
            separated_list0(char('.'), nom_u32),
        )
            .parse(i)
    };
    let (input, entries) = separated_list0(comma_not_before_entry, parse_entry).parse(input)?;

    Ok((input, CxEntry::MulticenterBonds(entries)))
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

/// Parse rb value: * (AsDrawn), -2, -1, 0, 2, 3, 4.
fn parse_rb_value(value: &[u8]) -> Result<Option<RingBondCount>, ()> {
    match value {
        b"*" | b"\\*" => Ok(Some(RingBondCount::AsDrawn)),
        b"-2" => Ok(Some(RingBondCount::AsDrawn)),
        b"-1" => Ok(Some(RingBondCount::NoRingBonds)),
        b"0" => Ok(None),
        b"2" => Ok(Some(RingBondCount::R2)),
        b"3" => Ok(Some(RingBondCount::R3)),
        b"4" => Ok(Some(RingBondCount::R4Plus)),
        _ => Err(()),
    }
}

/// Parse s value: * (AsDrawn), -2, -1, 0, 1-10.
fn parse_s_value(value: &[u8]) -> Result<Option<SubstitutionCount>, ()> {
    match value {
        b"*" | b"\\*" => Ok(Some(SubstitutionCount::AsDrawn)),
        b"-2" => Ok(Some(SubstitutionCount::AsDrawn)),
        b"-1" => Ok(Some(SubstitutionCount::NoSubstitution)),
        b"0" => Ok(None),
        b"1" => Ok(Some(SubstitutionCount::S1)),
        b"2" => Ok(Some(SubstitutionCount::S2)),
        b"3" => Ok(Some(SubstitutionCount::S3)),
        b"4" => Ok(Some(SubstitutionCount::S4)),
        b"5" => Ok(Some(SubstitutionCount::S5)),
        b"6" => Ok(Some(SubstitutionCount::S6Plus)),
        b"7" => Ok(Some(SubstitutionCount::S7)),
        b"8" => Ok(Some(SubstitutionCount::S8)),
        b"9" => Ok(Some(SubstitutionCount::S9)),
        b"10" => Ok(Some(SubstitutionCount::S10)),
        _ => Err(()),
    }
}

/// Parse ring bond count: `rb:idx:value,idx:value,...`
fn parse_ring_bond_count(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, pairs) = preceded(
        tag("rb:"),
        separated_list0(
            comma_not_before_entry,
            separated_pair(nom_u32, char(':'), take_while1(|b| b != b',' && b != b'|')),
        ),
    )
    .parse(input)?;

    let mut result = Vec::new();
    for (idx, val) in pairs {
        match parse_rb_value(val) {
            Ok(Some(rbc)) => result.push((idx, rbc)),
            Ok(None) => {}
            Err(()) => return Err(Err::Failure(NomError::new(input, ErrorKind::Verify))),
        }
    }

    Ok((input, CxEntry::RingBondCount(result)))
}

/// Parse substitution count: `s:idx:value,idx:value,...`
fn parse_substitution_count(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, pairs) = preceded(
        tag("s:"),
        separated_list0(
            comma_not_before_entry,
            separated_pair(nom_u32, char(':'), take_while1(|b| b != b',' && b != b'|')),
        ),
    )
    .parse(input)?;

    let mut result = Vec::new();
    for (idx, val) in pairs {
        match parse_s_value(val) {
            Ok(Some(sc)) => result.push((idx, sc)),
            Ok(None) => {}
            Err(()) => return Err(Err::Failure(NomError::new(input, ErrorKind::Verify))),
        }
    }

    Ok((input, CxEntry::SubstitutionCount(result)))
}

/// Parse unsaturated atoms: `u:idx,idx,...`
fn parse_unsaturated(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, indices) =
        preceded(tag("u:"), separated_list0(comma_not_before_entry, nom_u32)).parse(input)?;

    Ok((input, CxEntry::Unsaturated(indices)))
}

/// Parse ligand order: `LO:centerIdx:idx1.idx2.idx3,centerIdx2:idx1.idx2...`
fn parse_ligand_order(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let parse_entry = |i| {
        let (i, (center, neighbors)) = (
            terminated(nom_u32, char(':')),
            separated_list0(char('.'), nom_u32),
        )
            .parse(i)?;
        Ok((i, (center, neighbors)))
    };

    let (input, entries) = preceded(
        tag("LO:"),
        separated_list0(comma_not_before_entry, parse_entry),
    )
    .parse(input)?;

    Ok((input, CxEntry::LigandOrder(entries)))
}

/// Parse link nodes: `LN:atom:min.max` or `LN:atom:min.max.outer1.outer2`
fn parse_link_nodes(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let parse_entry = |i| {
        let (i, (atom_idx, values)) = (
            terminated(nom_u32, char(':')),
            separated_list0(char('.'), nom_u32),
        )
            .parse(i)?;
        if values.len() != 2 && values.len() != 4 {
            return Err(Err::Failure(NomError::new(i, ErrorKind::Verify)));
        }
        let min_repeat = values[0].min(255) as u8;
        let repeat_count = values[1].min(255) as u8;
        let (subs_index1, subs_index2) = if values.len() == 4 {
            (values[2], Some(values[3]))
        } else {
            (0u32, None)
        };
        let link = LinkAtom {
            min_repeat,
            repeat_count,
            subs_index1,
            subs_index2,
        };
        Ok((i, (atom_idx, link)))
    };

    let (input, entries) = preceded(
        tag("LN:"),
        separated_list0(comma_not_before_entry, parse_entry),
    )
    .parse(input)?;

    Ok((input, CxEntry::LinkNodes(entries)))
}

/// Parse data S-group: `SgD:atomIndices:name:data:queryOp:unit:tag:coords`
fn parse_sgroup_data(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, _) = tag("SgD:").parse(input)?;

    let parse_field = take_while(|b| b != b':');
    let (input, field_colons) = count((parse_field, char(':')), 7).parse(input)?;
    let (input, last_field) = take_until_entry_boundary(input)?;

    let mut segments: Vec<&[u8]> = field_colons.into_iter().map(|(s, _)| s).collect();
    segments.push(last_field);

    let atoms_input = segments[0];
    let atom_indices: Vec<u32> = atoms_input
        .split_str(",")
        .filter(|s| !s.is_empty())
        .filter_map(|s| {
            let s = s.to_str_lossy();
            s.trim().parse::<u32>().ok()
        })
        .collect();
    if atom_indices.is_empty() {
        return Err(Err::Failure(NomError::new(input, ErrorKind::Verify)));
    }

    let str_field = |i: usize| {
        unescape_html_entities(segments.get(i).copied().unwrap_or(&[]))
            .to_str_lossy()
            .into_owned()
    };

    let name = str_field(1);
    let data_content = str_field(2);
    let query_op = str_field(3);
    let unit = str_field(4);
    let tag_str = str_field(5);

    let data = SGroupData {
        field_type: SGroupDataType::Text,
        field_name: name,
        field_units: if unit.is_empty() { None } else { Some(unit) },
        query_identifier: if tag_str.is_empty() {
            None
        } else {
            Some(tag_str)
        },
        data_query_operator: if query_op.is_empty() {
            None
        } else {
            Some(query_op)
        },
        data_content: if data_content.is_empty() {
            None
        } else {
            Some(vec![data_content])
        },
    };

    let mut sgroup = SGroup::new(SGroupType::Data);
    sgroup.atom_indices = atom_indices;
    sgroup.data = Some(data);

    if let Some(coords) = segments.get(7) {
        if !coords.is_empty() && coords != b"(-1)" {
            let nums: Vec<f64> = coords
                .split_str(",")
                .filter_map(|p| p.to_str_lossy().trim().parse::<f64>().ok())
                .collect();
            if nums.len() == 4 {
                sgroup.bracket_coords = Some(SGroupBracketCoords {
                    bracket1: (nums[0], nums[1]),
                    bracket2: (nums[2], nums[3]),
                    bracket3: None,
                    bracket4: None,
                });
            } else if nums.len() == 8 {
                sgroup.bracket_coords = Some(SGroupBracketCoords {
                    bracket1: (nums[0], nums[1]),
                    bracket2: (nums[2], nums[3]),
                    bracket3: Some((nums[4], nums[5])),
                    bracket4: Some((nums[6], nums[7])),
                });
            }
        }
    }

    Ok((input, CxEntry::SgroupData(sgroup)))
}

fn parse_sgroup_type(s: &[u8]) -> Option<SGroupType> {
    match s {
        b"n" => Some(SGroupType::RepeatingUnit),
        b"mon" => Some(SGroupType::Monomer),
        b"mer" => Some(SGroupType::Mer),
        b"co" => Some(SGroupType::Copolymer),
        b"xl" => Some(SGroupType::Crosslink),
        b"mod" => Some(SGroupType::Modification),
        b"mix" => Some(SGroupType::Mixture),
        b"f" => Some(SGroupType::Formulation),
        b"any" => Some(SGroupType::AnyPolymer),
        b"gen" => Some(SGroupType::Generic),
        b"c" => Some(SGroupType::Component),
        b"grf" => Some(SGroupType::Graft),
        _ => None,
    }
}

fn parse_sgroup_subtype(s: &[u8]) -> Option<SGroupSubtype> {
    match s {
        b"alt" => Some(SGroupSubtype::Alternating),
        b"ran" => Some(SGroupSubtype::Random),
        b"blk" => Some(SGroupSubtype::Block),
        _ => None,
    }
}

fn parse_sgroup_connectivity(s: &[u8]) -> Option<SGroupConnectivity> {
    let base = s.split_str(",").next().unwrap_or(s);
    match base {
        b"hh" => Some(SGroupConnectivity::HeadToHead),
        b"ht" => Some(SGroupConnectivity::HeadToTail),
        b"eu" => Some(SGroupConnectivity::EitherUnknown),
        _ => None,
    }
}

fn parse_connectivity_flip(s: &[u8]) -> Option<bool> {
    let parts: Vec<&[u8]> = s.split_str(",").collect();
    if parts.len() < 2 {
        return None;
    }
    let flip_part = parts[1].to_str_lossy();
    match flip_part.trim() {
        "1" | "flip" | "true" => Some(true),
        "0" | "false" => Some(false),
        _ => None,
    }
}

fn parse_atom_list(s: &[u8]) -> Option<Vec<u32>> {
    let indices: Vec<u32> = s
        .split_str(",")
        .filter(|p| !p.is_empty())
        .filter_map(|p| p.to_str_lossy().trim().parse::<u32>().ok())
        .collect();
    if indices.is_empty() {
        None
    } else {
        Some(indices)
    }
}

fn parse_bond_list(s: &[u8]) -> Vec<u32> {
    s.split_str(".")
        .filter(|p| !p.is_empty())
        .filter_map(|p| p.to_str_lossy().trim().parse::<u32>().ok())
        .collect()
}

fn parse_bracket_orientation(s: &[u8]) -> Option<SGroupBracketOrientation> {
    match s {
        b"s" => Some(SGroupBracketOrientation::Straight),
        b"d" => Some(SGroupBracketOrientation::Down),
        _ => None,
    }
}

fn parse_bracket_style(s: &[u8]) -> Option<SGroupBracketStyle> {
    match s {
        b"b" => Some(SGroupBracketStyle::Default),
        b"c" => Some(SGroupBracketStyle::Curved),
        b"r" => Some(SGroupBracketStyle::TypeR),
        b"s" => Some(SGroupBracketStyle::TypeS),
        _ => None,
    }
}

fn parse_bracket_info(
    seg: &[u8],
) -> Option<(
    Option<SGroupBracketOrientation>,
    Option<SGroupBracketStyle>,
    Option<SGroupBracketCoords>,
)> {
    let parts: Vec<&[u8]> = seg.split(|&b| b == b',').collect();
    if parts.len() < 3 {
        return None;
    }
    let orientation = parse_bracket_orientation(parts[0]);
    let style = parse_bracket_style(parts[1]);
    let coords = if parts.len() >= 6 {
        let nums: Vec<f64> = parts[2..]
            .iter()
            .take(8)
            .filter_map(|p| p.to_str_lossy().trim().parse::<f64>().ok())
            .collect();
        if nums.len() == 8 {
            Some(SGroupBracketCoords {
                bracket1: (nums[0], nums[1]),
                bracket2: (nums[2], nums[3]),
                bracket3: Some((nums[4], nums[5])),
                bracket4: Some((nums[6], nums[7])),
            })
        } else if nums.len() == 4 {
            Some(SGroupBracketCoords {
                bracket1: (nums[0], nums[1]),
                bracket2: (nums[2], nums[3]),
                bracket3: None,
                bracket4: None,
            })
        } else {
            None
        }
    } else {
        None
    };
    if orientation.is_some() || style.is_some() || coords.is_some() {
        Some((orientation, style, coords))
    } else {
        None
    }
}

/// Parse polymer S-group: `Sg:type:subtype:atoms:subscript:connectivity:head:tail:bracket`
fn parse_sgroup(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let (input, _) = tag("Sg:").parse(input)?;
    let (input, content) = take_until_entry_boundary(input)?;

    let segments: Vec<&[u8]> = content.split(|&b| b == b':').collect();
    if segments.len() < 2 {
        return Err(Err::Failure(NomError::new(input, ErrorKind::Verify)));
    }

    let group_type = parse_sgroup_type(segments[0])
        .ok_or_else(|| Err::Failure(NomError::new(input, ErrorKind::Verify)))?;

    let mut atoms_idx = 1usize;
    let mut subtype = None;
    if segments.len() > 2 {
        if let Some(st) = parse_sgroup_subtype(segments[1]) {
            subtype = Some(st);
            atoms_idx = 2;
        } else if parse_atom_list(segments[1]).is_some() {
            atoms_idx = 1;
        }
    }

    let atom_indices = parse_atom_list(segments.get(atoms_idx).copied().unwrap_or(&[]))
        .ok_or_else(|| Err::Failure(NomError::new(input, ErrorKind::Verify)))?;

    let mut sgroup = SGroup::new(group_type);
    sgroup.group_subtype = subtype;
    sgroup.atom_indices = atom_indices;

    let mut idx = atoms_idx + 1;

    if idx < segments.len() {
        let seg = segments[idx];
        if let Some(conn) = parse_sgroup_connectivity(seg) {
            sgroup.connectivity = Some(conn);
            sgroup.connectivity_flip = parse_connectivity_flip(seg);
            idx += 1;
        } else {
            let s = unescape_html_entities(seg).to_str_lossy().into_owned();
            if !s.is_empty() {
                sgroup.subscript = Some(s);
            }
            idx += 1;
        }
    }

    if idx < segments.len() {
        let seg = segments[idx];
        if let Some(conn) = parse_sgroup_connectivity(seg) {
            sgroup.connectivity = Some(conn);
            if sgroup.connectivity_flip.is_none() {
                sgroup.connectivity_flip = parse_connectivity_flip(seg);
            }
            idx += 1;
        }
    }

    if idx < segments.len() {
        let head = parse_bond_list(segments[idx]);
        idx += 1;
        if idx < segments.len() {
            let tail = parse_bond_list(segments[idx]);
            sgroup.bond_indices = head.into_iter().chain(tail).collect();
            idx += 1;
        } else {
            sgroup.bond_indices = head;
        }
    }

    if idx < segments.len() {
        if let Some((orientation, style, coords)) = parse_bracket_info(segments[idx]) {
            sgroup.bracket_orientation = orientation;
            sgroup.bracket_style = style;
            sgroup.bracket_coords = coords;
        }
    }

    Ok((input, CxEntry::Sgroup(sgroup)))
}

/// Parse S-group hierarchy: `SgH:parentIdx1:child1.child2,parentIdx2:child1`
fn parse_sgroup_hierarchy(input: &[u8]) -> IResult<&[u8], CxEntry> {
    let parse_parent_children =
        separated_pair(nom_u32, char(':'), separated_list0(char('.'), nom_u32));
    let (input, pairs) = preceded(
        tag("SgH:"),
        separated_list0(comma_not_before_entry, parse_parent_children),
    )
    .parse(input)?;
    if pairs.is_empty() {
        return Err(Err::Failure(NomError::new(input, ErrorKind::Verify)));
    }
    Ok((input, CxEntry::SgroupHierarchy(pairs)))
}

/// Parse THB:/TLB:/TEB: tag. Format: THB:ligand:connection:lower:higher (comma-separated entries)
fn parse_bicyclo_stereo(input: &[u8]) -> IResult<&[u8], CxEntry> {
    fn parse_one(i: &[u8]) -> IResult<&[u8], BicycloStereoData> {
        let (i, ligand) = nom_u32.parse(i)?;
        let (i, _) = char(':').parse(i)?;
        let (i, connection) = nom_u32.parse(i)?;
        let (i, _) = char(':').parse(i)?;
        let (i, lower) = separated_list0(char('.'), nom_u32).parse(i)?;
        let (i, _) = char(':').parse(i)?;
        let (i, higher) = separated_list0(char('.'), nom_u32).parse(i)?;
        Ok((
            i,
            BicycloStereoData {
                ligand_atom: ligand,
                connection_atom: connection,
                lower_bridge_atoms: lower,
                higher_bridge_atoms: higher,
            },
        ))
    }

    let (input, (tag_bytes, entries)) = alt((
        (
            tag("THB:"),
            separated_list0(comma_not_before_entry, parse_one),
        ),
        (
            tag("TLB:"),
            separated_list0(comma_not_before_entry, parse_one),
        ),
        (
            tag("TEB:"),
            separated_list0(comma_not_before_entry, parse_one),
        ),
    ))
    .parse(input)?;

    let variant = match tag_bytes {
        b"THB:" => BicycloStereo::TowardsHigherBridge,
        b"TLB:" => BicycloStereo::TowardsLowerBridge,
        b"TEB:" => BicycloStereo::TowardsEitherBridge,
        _ => return Err(Err::Failure(NomError::new(input, ErrorKind::Tag))),
    };

    let entries: Vec<BicycloStereo> = entries.into_iter().map(variant).collect();
    if entries.is_empty() {
        return Err(Err::Failure(NomError::new(input, ErrorKind::Verify)));
    }
    Ok((input, CxEntry::BicycloStereo(entries)))
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

/// Take bytes until entry boundary: pipe, or comma followed by tag start and colon.
/// Requires ":" after the tag name so "s,b,1,2" (bracket coords) is not split at ",b".
fn take_until_entry_boundary(input: &[u8]) -> IResult<&[u8], &[u8]> {
    let mut i = 0;
    while i < input.len() {
        match input[i] {
            b'|' => return Ok((&input[i..], &input[..i])),
            b',' => {
                let mut j = i + 1;
                while j < input.len() && input[j].is_ascii_whitespace() {
                    j += 1;
                }
                if j < input.len() && is_cx_tag_start(input[j] as char) {
                    let mut k = j + 1;
                    while k < input.len() && (input[k].is_ascii_alphanumeric() || input[k] == b'_')
                    {
                        k += 1;
                    }
                    if k < input.len() && input[k] == b':' {
                        return Ok((&input[i..], &input[..i]));
                    }
                }
            }
            _ => {}
        }
        i += 1;
    }
    Ok((&input[i..], &input[..i]))
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
        if input[i] == b',' {
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
    use crate::table_ir::atom::{BicycloStereo, BicycloStereoData};
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
    #[case::cis_bonds(b"|c:0,1|", vec![CxEntry::CisBonds(vec![0, 1])])]
    #[case::trans_bonds(b"|t:0,1|", vec![CxEntry::TransBonds(vec![0, 1])])]
    #[case::unspec_bonds(b"|ctu:0,1|", vec![CxEntry::UnspecBonds(vec![0, 1])])]
    #[case::multicenter_bonds(b"|m:0:3.4,2:1.5|", vec![CxEntry::MulticenterBonds(vec![(0, vec![3, 4]), (2, vec![1, 5])])])]
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
    #[case::multicenter_bonds(b"|m:0:3.4,2:1.5|", vec![CxEntry::MulticenterBonds(vec![(0, vec![3, 4]), (2, vec![1, 5])])])]
    #[case::stereo_absolute(b"|a:0,1,2|", vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Absolute, atoms: vec![0, 1, 2] })])]
    #[case::stereo_or(b"|o1:0,1|", vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Or(1), atoms: vec![0, 1] })])]
    #[case::stereo_and(b"|&1:0,1|", vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::And(1), atoms: vec![0, 1] })])]
    #[case::relative_stereo(b"|r|", vec![CxEntry::RelativeStereo])]
    #[case::atom_properties(b"|atomProp:0.key.value|", vec![CxEntry::AtomProperties(vec![(0, "key".to_string(), "value".to_string())])])]
    #[case::coordinates_2d(b"|(1.5,2.5;3.5,4.5)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.5, 2.5, 0.0), Point3D::new(3.5, 4.5, 0.0)])])]
    #[case::coordinates_3d(b"|(1,2,3;4,5,6)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0), Point3D::new(4.0, 5.0, 6.0)])])]
    #[case::bicyclo_thb(b"|THB:12:11:2.4.3:7.10.8|", vec![CxEntry::BicycloStereo(vec![BicycloStereo::TowardsHigherBridge(BicycloStereoData {ligand_atom: 12, connection_atom: 11,
        lower_bridge_atoms: vec![2, 4, 3], higher_bridge_atoms: vec![7, 10, 8]})])])]
    #[case::bicyclo_tlb(b"|TLB:13:11:2.4.3:7.10.8|", vec![CxEntry::BicycloStereo(vec![BicycloStereo::TowardsLowerBridge(BicycloStereoData {ligand_atom: 13, connection_atom: 11,
        lower_bridge_atoms: vec![2, 4, 3], higher_bridge_atoms: vec![7, 10, 8]})])])]
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
    #[case::multicenter_bonds(vec![CxEntry::MulticenterBonds(vec![(0, vec![1, 2]), (2, vec![0, 1])])],
        |mol: &Molecule| mol.multicenter_bonds.len() == 2 &&
          mol.multicenter_bonds[0] == MulticenterBond::new(vec![MulticenterSet::new(vec![0]), MulticenterSet::new(vec![1, 2])]) &&
          mol.multicenter_bonds[1] == MulticenterBond::new(vec![MulticenterSet::new(vec![2]), MulticenterSet::new(vec![0, 1])]))]
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
    #[case::multicenter_bonds(vec![CxEntry::MulticenterBonds(vec![(0, vec![1, 2]), (2, vec![0, 1])])],
        |mol: &ExtendedMolecule| mol.multicenter_bonds.len() == 2 &&
          mol.multicenter_bonds[0] == MulticenterBond::new(vec![MulticenterSet::new(vec![0]), MulticenterSet::new(vec![1, 2])]) &&
          mol.multicenter_bonds[1] == MulticenterBond::new(vec![MulticenterSet::new(vec![2]), MulticenterSet::new(vec![0, 1])]))]
    #[case::stereo_group_absolute(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Absolute, atoms: vec![0, 1] })],
        |mol: &ExtendedMolecule| mol.stereo_interpretation == Some(StereoInterpretation::Absolute) && mol.cx_data.is_none())]
    #[case::stereo_group_or(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Or(1), atoms: vec![0] })],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_groups.get(&1)) == Some(&StereoSet { atoms: vec![0], mode: StereoSetMode::Correlated }))]
    #[case::stereo_group_and(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::And(2), atoms: vec![1] })],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_groups.get(&2)) == Some(&StereoSet { atoms: vec![1], mode: StereoSetMode::Independent }))]
    #[case::relative_stereo(vec![CxEntry::RelativeStereo],
        |mol: &ExtendedMolecule| mol.stereo_interpretation == Some(StereoInterpretation::Relative) && mol.cx_data.is_none())]
    #[case::atom_properties(vec![CxEntry::AtomProperties(vec![(0, "key".to_string(), "value".to_string())])], |mol: &ExtendedMolecule| mol.atoms[0].properties.get("key") == Some(&"value".to_string()))]
    #[case::bicyclo_stereo(vec![CxEntry::BicycloStereo(vec![BicycloStereo::TowardsHigherBridge(BicycloStereoData{ligand_atom: 12, connection_atom: 11,
        lower_bridge_atoms: vec![2, 4, 3], higher_bridge_atoms: vec![7, 10, 8], }), ])], |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.bicyclo_stereo.as_ref()).map(|v| v.len()) == Some(1))]
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
