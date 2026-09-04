//! CXSMILES annotation block parser
//!
//! Parses the `|...|` extension block in CXSMILES format.
//! Two parsers are provided:
//! - `parse_cx_annotations`: basic annotations only (for Molecule)
//! - `parse_extended_cx_annotations`: all annotations (for ExtendedMolecule)

use std::collections::BTreeMap;

use bstr::ByteSlice;
use umol_chem::spin::SpinMultiplicity;
use umol_geometric_core::Point3D;
use winnow::ascii::{dec_uint, float};
use winnow::combinator::{
    alt, delimited, not, opt, peek, preceded, repeat, separated, separated_pair, terminated,
};
use winnow::error::{ErrMode, ParserError};
use winnow::stream::Stream;
use winnow::token::{one_of, take_while};
use winnow::{ModalResult, Parser};

use super::super::config::SmilesSyntaxFlags;
use super::super::error::ParseError;
use super::utils::{split_escaped_semicolons, unescape_html_entities};
use crate::table_ir::bond::BondNoncovalent;
use crate::table_ir::{
    BicycloStereo, BicycloStereoData, BondDonation, BondOrder, BondStereo, BondWedge,
    ConfigurationScope, CxAnnotationData, ExtendedMolecule, ExtendedReaction, LinkAtom, Molecule,
    MulticenterBond, MulticenterSet, Reaction, RingBondCount, SGroup, SGroupBracketCoords,
    SGroupBracketOrientation, SGroupBracketStyle, SGroupConnectivity, SGroupData, SGroupDataType,
    SGroupSubtype, SGroupType, StereoSet, StereoSetRelation, SubstitutionCount, UnsaturatedAtom,
};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum CxParseError {
    Syntax,
    InvalidTag,
}

impl<I: Stream> ParserError<I> for CxParseError {
    type Inner = Self;

    fn from_input(_: &I) -> Self {
        Self::Syntax
    }

    fn into_inner(self) -> Result<Self::Inner, Self> {
        Ok(self)
    }
}

type PResult<T> = ModalResult<T, CxParseError>;

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
    Radicals(Vec<(u32, (u8, Option<SpinMultiplicity>))>),
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
pub fn parse_cx_annotations<'i>(
    input: &'i [u8],
    flags: SmilesSyntaxFlags,
) -> Result<Vec<CxEntry>, ParseError> {
    let skip_unknown_cx_tags = flags.contains(SmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS);
    parse_cx_block(input, |i: &mut &'i [u8]| {
        parse_basic_entry(i, skip_unknown_cx_tags)
    })
}

/// Parse extended CX annotations (for ExtendedMolecule)
pub fn parse_extended_cx_annotations<'i>(
    input: &'i [u8],
    flags: SmilesSyntaxFlags,
) -> Result<Vec<CxEntry>, ParseError> {
    let skip_unknown_cx_tags = flags.contains(SmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS);
    parse_cx_block(input, |i: &mut &'i [u8]| {
        parse_extended_entry(i, skip_unknown_cx_tags)
    })
}

/// Maps CXSMILES bond indices (close-order: ring-closure bonds counted at their
/// closing digit) onto our bond list (open-order: ring bonds at their opening
/// position). Built from the parser's per-ring `(close_rank, open_index)` record,
/// so it is O(ring count); empty (identity) for ring-free molecules.
pub struct BondIndexMap {
    /// `(close_rank, open_index)` per ring-closure bond, in closing order.
    ring_bonds: Vec<(usize, usize)>,
    /// Ring-closure bonds' open indices, ascending.
    ring_open_sorted: Vec<usize>,
    bond_count: usize,
}

impl BondIndexMap {
    pub fn new(ring_bonds: Vec<(usize, usize)>, bond_count: usize) -> Self {
        let mut ring_open_sorted: Vec<usize> = ring_bonds.iter().map(|&(_, open)| open).collect();
        ring_open_sorted.sort_unstable();
        Self {
            ring_bonds,
            ring_open_sorted,
            bond_count,
        }
    }

    /// No ring closures, so close-order equals open-order — translation is identity.
    pub fn is_identity(&self) -> bool {
        self.ring_bonds.is_empty()
    }

    /// Translate a CXSMILES (close-order) bond index to our (open-order) index,
    /// or `None` if out of range.
    pub fn translate(&self, cx_idx: u32) -> Option<u32> {
        let k = cx_idx as usize;
        if k >= self.bond_count {
            return None;
        }
        // A ring bond closing at rank k maps to its opening position.
        if let Some(&(_, open)) = self.ring_bonds.iter().find(|&&(rank, _)| rank == k) {
            return Some(open as u32);
        }
        // Otherwise it is a sequential bond: the s-th one in close order, whose
        // open index is the s-th position once the ring opening positions are skipped.
        let s = k - self
            .ring_bonds
            .iter()
            .filter(|&&(rank, _)| rank < k)
            .count();
        let mut open = s;
        for &r in &self.ring_open_sorted {
            if r <= open {
                open += 1;
            }
        }
        Some(open as u32)
    }
}

/// Rewrite the close-order bond indices in `entries` to our open-order indices.
/// A no-op for ring-free molecules.
pub fn remap_cx_bond_indices(
    entries: &mut [CxEntry],
    map: &BondIndexMap,
) -> Result<(), ParseError> {
    if map.is_identity() {
        return Ok(());
    }
    let translate = |bond_idx: &mut u32| -> Result<(), ParseError> {
        *bond_idx = map
            .translate(*bond_idx)
            .ok_or(ParseError::BondIndexOutOfBounds {
                bond_idx: *bond_idx,
            })?;
        Ok(())
    };
    for entry in entries {
        match entry {
            CxEntry::WigglyBonds(items) => {
                for (_, bond_idx, _) in items {
                    translate(bond_idx)?;
                }
            }
            CxEntry::CisBonds(indices)
            | CxEntry::TransBonds(indices)
            | CxEntry::UnspecBonds(indices) => {
                for bond_idx in indices {
                    translate(bond_idx)?;
                }
            }
            CxEntry::CoordinateBonds(pairs) | CxEntry::HydrogenBonds(pairs) => {
                for (_, bond_idx) in pairs {
                    translate(bond_idx)?;
                }
            }
            _ => {}
        }
    }
    Ok(())
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
                for (idx, (unpaired_electrons, multiplicity)) in radicals {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.unpaired_electrons = Some(unpaired_electrons);
                    atom.multiplicity = multiplicity;
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
    let mut configuration_scope: Option<ConfigurationScope> = None;
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
                for (idx, unpaired_electrons) in radicals {
                    let Some(atom) = mol.atoms.get_mut(idx as usize) else {
                        return Err(ParseError::AtomIndexOutOfBounds { atom_idx: idx });
                    };
                    atom.unpaired_electrons = Some(unpaired_electrons.0);
                    atom.multiplicity = unpaired_electrons.1;
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
                        // Absolute atoms don't need group storage; configuration_scope captures this
                        configuration_scope = Some(ConfigurationScope::Absolute);
                    }
                    StereoGroupType::Or(n) => {
                        stereo_groups
                            .entry(n)
                            .and_modify(|s| s.atoms.extend(sg.atoms.iter().copied()))
                            .or_insert(StereoSet {
                                atoms: sg.atoms,
                                relation: StereoSetRelation::Correlated,
                            });
                    }
                    StereoGroupType::And(n) => {
                        stereo_groups
                            .entry(n)
                            .and_modify(|s| s.atoms.extend(sg.atoms.iter().copied()))
                            .or_insert(StereoSet {
                                atoms: sg.atoms,
                                relation: StereoSetRelation::Independent,
                            });
                    }
                }
            }
            CxEntry::RelativeStereo => {
                configuration_scope = Some(ConfigurationScope::Relative);
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

    mol.configuration_scope = configuration_scope;

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
fn comma_before_entry(input: &mut &[u8]) -> PResult<u8> {
    terminated(b',', peek(one_of(is_cx_tag_start))).parse_next(input)
}

fn parse_u32(input: &mut &[u8]) -> PResult<u32> {
    dec_uint.parse_next(input)
}

fn parse_dot_indices(input: &mut &[u8]) -> PResult<Vec<u32>> {
    separated(0.., parse_u32, b'.').parse_next(input)
}

fn parse_cx_block<'inp>(
    input: &'inp [u8],
    entry_parser: impl Parser<&'inp [u8], Option<CxEntry>, ErrMode<CxParseError>>,
) -> Result<Vec<CxEntry>, ParseError> {
    let mut remaining = input;
    let result: PResult<Vec<Option<CxEntry>>> = delimited(
        opt(b'|'),
        separated(0.., entry_parser, comma_before_entry),
        opt(b'|'),
    )
    .parse_next(&mut remaining);
    match result {
        Ok(options) if remaining.is_empty() => Ok(options.into_iter().flatten().collect()),
        Ok(_) => Err(ParseError::InvalidToken { pos: 0 }),
        Err(ErrMode::Cut(CxParseError::InvalidTag)) => Err(ParseError::InvalidCxTag { pos: 0 }),
        Err(ErrMode::Cut(CxParseError::Syntax))
        | Err(ErrMode::Backtrack(_))
        | Err(ErrMode::Incomplete(_)) => Err(ParseError::InvalidToken { pos: 0 }),
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
fn parse_basic_entry(input: &mut &[u8], skip_unknown_cx_tags: bool) -> PResult<Option<CxEntry>> {
    if input.is_empty() {
        return Err(ErrMode::Backtrack(CxParseError::Syntax));
    }
    if input.first() == Some(&b'|') {
        // End of CX block.
        return Err(ErrMode::Backtrack(CxParseError::Syntax));
    }

    let start = *input;
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
    .parse_next(input)
    {
        Ok(entry) => Ok(Some(entry)),
        Err(ErrMode::Backtrack(_)) => {
            *input = start;
            parse_unknown_entry(input, skip_unknown_cx_tags)
        }
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
fn parse_extended_entry(input: &mut &[u8], skip_unknown_cx_tags: bool) -> PResult<Option<CxEntry>> {
    if input.is_empty() {
        return Err(ErrMode::Backtrack(CxParseError::Syntax));
    }
    if input.first() == Some(&b'|') {
        // End of CX block.
        return Err(ErrMode::Backtrack(CxParseError::Syntax));
    }

    let start = *input;
    match alt((
        alt((
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
        )),
        alt((parse_ligand_order, parse_link_nodes)),
        alt((parse_sgroup_data, parse_sgroup, parse_sgroup_hierarchy)),
        parse_bicyclo_stereo,
    ))
    .parse_next(input)
    {
        Ok(entry) => Ok(Some(entry)),
        Err(ErrMode::Backtrack(_)) => {
            *input = start;
            parse_unknown_entry(input, skip_unknown_cx_tags)
        }
        Err(e) => Err(e),
    }
}

/// Parse coordinates (x,y) or (x,y,z) for a single atom.
/// Missing components default to 0.0.
fn parse_atom_coordinates(input: &mut &[u8]) -> PResult<Point3D> {
    let coords: Vec<Option<f64>> = separated(0.., opt(float), b',').parse_next(input)?;
    if coords.is_empty() {
        return Ok(Point3D::zero());
    }
    if coords.len() > 3 {
        return Err(ErrMode::Cut(CxParseError::Syntax));
    }
    let x = coords.first().copied().flatten().unwrap_or(0.0);
    let y = coords.get(1).copied().flatten().unwrap_or(0.0);
    let z = coords.get(2).copied().flatten().unwrap_or(0.0);
    Ok(Point3D::new(x, y, z))
}

/// Parse coordinates block: `(x,y,z;x,y,z;...)`
/// Empty parens `()` means no atoms have coordinates.
fn parse_coordinates(input: &mut &[u8]) -> PResult<CxEntry> {
    let coords = alt((
        b"()".value(vec![]),
        delimited(b'(', separated(0.., parse_atom_coordinates, b';'), b')'),
    ))
    .parse_next(input)?;

    Ok(CxEntry::Coordinates(coords))
}

/// Parse labels `$label;label;...$` or values `$_AV:value;value;...$`
fn parse_labels(input: &mut &[u8]) -> PResult<CxEntry> {
    let inner = delimited(b'$', take_while(1.., |b| b != b'$'), b'$').parse_next(input)?;

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
        Ok(CxEntry::Values(result))
    } else {
        Ok(CxEntry::Labels(result))
    }
}

/// Convert CXSMILES radical code (1-7) to unpaired electrons.
fn convert_radical_code(code: u8) -> (u8, Option<SpinMultiplicity>) {
    match code {
        1 => (1, None),
        2 => (2, None),
        3 => (2, Some(SpinMultiplicity::SINGLET)),
        4 => (2, Some(SpinMultiplicity::TRIPLET)),
        5 => (3, None),
        6 => (3, Some(SpinMultiplicity::DOUBLET)),
        7 => (3, Some(SpinMultiplicity::QUARTET)),
        _ => (1, None),
    }
}

/// Parse a single radical group: `^n:idx,idx,...`
fn parse_radical_group(input: &mut &[u8]) -> PResult<(u8, Vec<u32>)> {
    let code = delimited(b'^', one_of(&b"1234567"[..]), b':').parse_next(input)?;
    let indices = separated(0.., parse_u32, comma_not_before_entry).parse_next(input)?;
    Ok((code - b'0', indices))
}

/// Parse radicals: `^n:idx,idx,...` (one or more groups).
fn parse_radicals(input: &mut &[u8]) -> PResult<CxEntry> {
    let groups: Vec<_> = repeat(1.., parse_radical_group).parse_next(input)?;

    let result: Vec<_> = groups
        .into_iter()
        .flat_map(|(code, indices)| {
            let unpaired_electrons = convert_radical_code(code);
            indices
                .into_iter()
                .map(move |idx| (idx, unpaired_electrons))
        })
        .collect();

    Ok(CxEntry::Radicals(result))
}

/// Parse wiggly bonds: `w:`, `wU:`, `wD:` followed by atom.bond pairs.
fn parse_wiggly_bonds(input: &mut &[u8]) -> PResult<CxEntry> {
    let wedge_type = alt((
        b"wU:".value(BondWedge::EitherUp),
        b"wD:".value(BondWedge::EitherDown),
        b"w:".value(BondWedge::Either),
    ))
    .parse_next(input)?;

    let pairs: Vec<(u32, u32)> = separated(
        0..,
        separated_pair(parse_u32, b'.', parse_u32),
        comma_not_before_entry,
    )
    .parse_next(input)?;

    let result: Vec<_> = pairs
        .into_iter()
        .map(|(atom_idx, bond_idx)| (atom_idx, bond_idx, wedge_type))
        .collect();
    Ok(CxEntry::WigglyBonds(result))
}

/// Parse cis/trans bond annotations: `c:`, `t:`, `ctu:`.
fn parse_cis_trans(input: &mut &[u8]) -> PResult<CxEntry> {
    let kind = alt((b"ctu:".value('u'), b"c:".value('c'), b"t:".value('t'))).parse_next(input)?;

    let indices = separated(0.., parse_u32, comma_not_before_entry).parse_next(input)?;

    match kind {
        'c' => Ok(CxEntry::CisBonds(indices)),
        't' => Ok(CxEntry::TransBonds(indices)),
        'u' => Ok(CxEntry::UnspecBonds(indices)),
        _ => unreachable!("unknown cis/trans/ctu tag"),
    }
}

/// Parse coordinate (dative) bonds: `C:atom.bond,...`
fn parse_coordinate_bonds(input: &mut &[u8]) -> PResult<CxEntry> {
    let pairs = preceded(
        b"C:",
        separated(
            0..,
            separated_pair(parse_u32, b'.', parse_u32),
            comma_not_before_entry,
        ),
    )
    .parse_next(input)?;

    Ok(CxEntry::CoordinateBonds(pairs))
}

/// Parse hydrogen bonds: `H:atom.bond,...`
fn parse_hydrogen_bonds(input: &mut &[u8]) -> PResult<CxEntry> {
    let pairs = preceded(
        b"H:",
        separated(
            0..,
            separated_pair(parse_u32, b'.', parse_u32),
            comma_not_before_entry,
        ),
    )
    .parse_next(input)?;

    Ok(CxEntry::HydrogenBonds(pairs))
}

/// Parse lone pairs: `LP:idx,idx,...` (unspecified count) or `lp:idx:count,...` (explicit count).
fn parse_lone_pairs(input: &mut &[u8]) -> PResult<CxEntry> {
    // Try lp: format first (has explicit counts)
    if input.starts_with(b"lp:") {
        *input = &input[3..];
        let entries: Vec<(u32, u32)> = separated(
            0..,
            separated_pair(parse_u32, b':', parse_u32),
            comma_not_before_entry,
        )
        .parse_next(input)?;
        let result: Vec<_> = entries
            .into_iter()
            .map(|(idx, count)| (idx, count as u8))
            .collect();
        return Ok(CxEntry::LonePairs(result));
    }

    // LP: format (implicit count, treated as 1 per atom)
    let indices: Vec<u32> =
        preceded(b"LP:", separated(0.., parse_u32, comma_not_before_entry)).parse_next(input)?;
    let result: Vec<_> = indices.into_iter().map(|idx| (idx, 1u8)).collect();
    Ok(CxEntry::LonePairs(result))
}

/// Parse multicenter bonds: `m:central:ligand.ligand,...`
fn parse_multicenter(input: &mut &[u8]) -> PResult<CxEntry> {
    b"m:".parse_next(input)?;
    let entries = separated(
        0..,
        (terminated(parse_u32, b':'), separated(0.., parse_u32, b'.')),
        comma_not_before_entry,
    )
    .parse_next(input)?;

    Ok(CxEntry::MulticenterBonds(entries))
}

/// Parse fragment groups: `f:atom.atom.atom,...`
fn parse_fragment_groups(input: &mut &[u8]) -> PResult<CxEntry> {
    let groups: Vec<Vec<u32>> = preceded(
        b"f:",
        separated(0.., parse_dot_indices, comma_not_before_entry),
    )
    .parse_next(input)?;

    let non_empty: Vec<_> = groups.into_iter().filter(|g| !g.is_empty()).collect();
    Ok(CxEntry::FragmentGroups(non_empty))
}

/// Parse absolute stereo group: `a:idx,idx,...`
fn parse_stereo_absolute(input: &mut &[u8]) -> PResult<CxEntry> {
    let atoms =
        preceded(b"a:", separated(0.., parse_u32, comma_not_before_entry)).parse_next(input)?;

    Ok(CxEntry::StereoGroup(StereoGroup {
        group_type: StereoGroupType::Absolute,
        atoms,
    }))
}

/// Parse OR/AND stereo group: `o<n>:idx,idx,...` or `&<n>:idx,idx,...`
fn parse_stereo_or_and(input: &mut &[u8]) -> PResult<CxEntry> {
    let (is_or, group_num, _, atoms) = (
        alt((b'o'.value(true), b'&'.value(false))),
        parse_u32,
        b':',
        separated(0.., parse_u32, comma_not_before_entry),
    )
        .parse_next(input)?;

    let group_type = if is_or {
        StereoGroupType::Or(group_num)
    } else {
        StereoGroupType::And(group_num)
    };

    Ok(CxEntry::StereoGroup(StereoGroup { group_type, atoms }))
}

/// Parse a single atom property entry: `idx.key.value`
fn parse_atom_prop_entry(input: &mut &[u8]) -> PResult<(u32, String, String)> {
    let (idx, _, key_bytes, _, value_bytes) = (
        parse_u32,
        b'.',
        take_while(1.., |b| b != b'.'),
        b'.',
        take_while(1.., |b| b != b':' && b != b',' && b != b'|'),
    )
        .parse_next(input)?;

    let key = unescape_html_entities(key_bytes)
        .to_str_lossy()
        .into_owned();
    let value = unescape_html_entities(value_bytes)
        .to_str_lossy()
        .into_owned();

    Ok((idx, key, value))
}

/// Parse atom properties: `atomProp:idx.key.value:idx.key.value...`
fn parse_atom_properties(input: &mut &[u8]) -> PResult<CxEntry> {
    let props =
        preceded(b"atomProp:", separated(0.., parse_atom_prop_entry, b':')).parse_next(input)?;

    Ok(CxEntry::AtomProperties(props))
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
fn parse_ring_bond_count(input: &mut &[u8]) -> PResult<CxEntry> {
    let pairs: Vec<(u32, &[u8])> = preceded(
        b"rb:",
        separated(
            0..,
            separated_pair(parse_u32, b':', take_while(1.., |b| b != b',' && b != b'|')),
            comma_not_before_entry,
        ),
    )
    .parse_next(input)?;

    let mut result = Vec::new();
    for (idx, val) in pairs {
        match parse_rb_value(val) {
            Ok(Some(rbc)) => result.push((idx, rbc)),
            Ok(None) => {}
            Err(()) => return Err(ErrMode::Cut(CxParseError::InvalidTag)),
        }
    }

    Ok(CxEntry::RingBondCount(result))
}

/// Parse substitution count: `s:idx:value,idx:value,...`
fn parse_substitution_count(input: &mut &[u8]) -> PResult<CxEntry> {
    let pairs: Vec<(u32, &[u8])> = preceded(
        b"s:",
        separated(
            0..,
            separated_pair(parse_u32, b':', take_while(1.., |b| b != b',' && b != b'|')),
            comma_not_before_entry,
        ),
    )
    .parse_next(input)?;

    let mut result = Vec::new();
    for (idx, val) in pairs {
        match parse_s_value(val) {
            Ok(Some(sc)) => result.push((idx, sc)),
            Ok(None) => {}
            Err(()) => return Err(ErrMode::Cut(CxParseError::InvalidTag)),
        }
    }

    Ok(CxEntry::SubstitutionCount(result))
}

/// Parse unsaturated atoms: `u:idx,idx,...`
fn parse_unsaturated(input: &mut &[u8]) -> PResult<CxEntry> {
    let indices =
        preceded(b"u:", separated(0.., parse_u32, comma_not_before_entry)).parse_next(input)?;

    Ok(CxEntry::Unsaturated(indices))
}

/// Parse ligand order: `LO:centerIdx:idx1.idx2.idx3,centerIdx2:idx1.idx2...`
fn parse_ligand_order(input: &mut &[u8]) -> PResult<CxEntry> {
    let entries = preceded(
        b"LO:",
        separated(
            0..,
            (terminated(parse_u32, b':'), separated(0.., parse_u32, b'.')),
            comma_not_before_entry,
        ),
    )
    .parse_next(input)?;

    Ok(CxEntry::LigandOrder(entries))
}

/// Parse link nodes: `LN:atom:min.max` or `LN:atom:min.max.outer1.outer2`
fn parse_link_nodes(input: &mut &[u8]) -> PResult<CxEntry> {
    let parse_entry = |i: &mut &[u8]| {
        let (atom_idx, values): (u32, Vec<u32>) =
            (terminated(parse_u32, b':'), separated(0.., parse_u32, b'.')).parse_next(i)?;
        if values.len() != 2 && values.len() != 4 {
            return Err(ErrMode::Cut(CxParseError::InvalidTag));
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
        Ok((atom_idx, link))
    };

    let entries =
        preceded(b"LN:", separated(0.., parse_entry, comma_not_before_entry)).parse_next(input)?;

    Ok(CxEntry::LinkNodes(entries))
}

/// Parse data S-group: `SgD:atomIndices:name:data:queryOp:unit:tag:coords`
fn parse_sgroup_data(input: &mut &[u8]) -> PResult<CxEntry> {
    b"SgD:".parse_next(input)?;

    let field_colons: Vec<(&[u8], u8)> =
        repeat(7, (take_while(0.., |b| b != b':'), b':')).parse_next(input)?;
    let last_field = take_until_entry_boundary(input)?;

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
        return Err(ErrMode::Cut(CxParseError::InvalidTag));
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

    Ok(CxEntry::SgroupData(sgroup))
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
fn parse_sgroup(input: &mut &[u8]) -> PResult<CxEntry> {
    b"Sg:".parse_next(input)?;
    let content = take_until_entry_boundary(input)?;

    let segments: Vec<&[u8]> = content.split(|&b| b == b':').collect();
    if segments.len() < 2 {
        return Err(ErrMode::Cut(CxParseError::InvalidTag));
    }

    let group_type =
        parse_sgroup_type(segments[0]).ok_or(ErrMode::Cut(CxParseError::InvalidTag))?;

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
        .ok_or(ErrMode::Cut(CxParseError::InvalidTag))?;

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

    Ok(CxEntry::Sgroup(sgroup))
}

/// Parse S-group hierarchy: `SgH:parentIdx1:child1.child2,parentIdx2:child1`
fn parse_sgroup_hierarchy(input: &mut &[u8]) -> PResult<CxEntry> {
    let parse_parent_children = separated_pair(parse_u32, b':', separated(0.., parse_u32, b'.'));
    let pairs: Vec<(u32, Vec<u32>)> = preceded(
        b"SgH:",
        separated(0.., parse_parent_children, comma_not_before_entry),
    )
    .parse_next(input)?;
    if pairs.is_empty() {
        return Err(ErrMode::Cut(CxParseError::InvalidTag));
    }
    Ok(CxEntry::SgroupHierarchy(pairs))
}

/// Parse THB:/TLB:/TEB: tag. Format: THB:ligand:connection:lower:higher (comma-separated entries)
fn parse_bicyclo_stereo(input: &mut &[u8]) -> PResult<CxEntry> {
    fn parse_one(input: &mut &[u8]) -> PResult<BicycloStereoData> {
        let ligand = parse_u32(input)?;
        b':'.parse_next(input)?;
        let connection = parse_u32(input)?;
        b':'.parse_next(input)?;
        let lower = separated(0.., parse_u32, b'.').parse_next(input)?;
        b':'.parse_next(input)?;
        let higher = separated(0.., parse_u32, b'.').parse_next(input)?;
        Ok(BicycloStereoData {
            ligand_atom: ligand,
            connection_atom: connection,
            lower_bridge_atoms: lower,
            higher_bridge_atoms: higher,
        })
    }

    let (tag_bytes, entries): (&[u8], Vec<BicycloStereoData>) = alt((
        (b"THB:", separated(0.., parse_one, comma_not_before_entry)),
        (b"TLB:", separated(0.., parse_one, comma_not_before_entry)),
        (b"TEB:", separated(0.., parse_one, comma_not_before_entry)),
    ))
    .parse_next(input)?;

    let variant = match tag_bytes {
        b"THB:" => BicycloStereo::TowardsHigherBridge,
        b"TLB:" => BicycloStereo::TowardsLowerBridge,
        b"TEB:" => BicycloStereo::TowardsEitherBridge,
        _ => return Err(ErrMode::Cut(CxParseError::Syntax)),
    };

    let entries: Vec<BicycloStereo> = entries.into_iter().map(variant).collect();
    if entries.is_empty() {
        return Err(ErrMode::Cut(CxParseError::InvalidTag));
    }
    Ok(CxEntry::BicycloStereo(entries))
}

/// Parse relative stereo tag: `r`.
fn parse_relative_stereo(input: &mut &[u8]) -> PResult<CxEntry> {
    let start = *input;
    b'r'.parse_next(input)?;

    // Only accept `r` as a standalone tag (`r` or `r:...`), not as the start of
    // another tag (e.g. `rb:`).
    if !input.is_empty() && !matches!(input[0], b':' | b',' | b'|') {
        *input = start;
        return Err(ErrMode::Backtrack(CxParseError::Syntax));
    }

    // `r:...` is used for reaction/multicomponent cases to list the fragment indices with relative
    // configuration. This isn't meaningful in our (molecule) TableIR parsing, so reject it.
    if input.first() == Some(&b':') {
        return Err(ErrMode::Cut(CxParseError::InvalidTag));
    }

    Ok(CxEntry::RelativeStereo)
}

/// Check whether a character can start a CX entry/tag.
fn is_cx_tag_start(c: u8) -> bool {
    c.is_ascii_alphabetic() || matches!(c, b'(' | b'$' | b'^' | b'&')
}

/// Take bytes until entry boundary: pipe, or comma followed by tag start and colon.
/// Requires ":" after the tag name so "s,b,1,2" (bracket coords) is not split at ",b".
fn take_until_entry_boundary<'i>(input: &mut &'i [u8]) -> PResult<&'i [u8]> {
    let start = *input;
    let mut i = 0;
    while i < start.len() {
        match start[i] {
            b'|' => {
                *input = &start[i..];
                return Ok(&start[..i]);
            }
            b',' => {
                let mut j = i + 1;
                while j < start.len() && start[j].is_ascii_whitespace() {
                    j += 1;
                }
                if j < start.len() && is_cx_tag_start(start[j]) {
                    let mut k = j + 1;
                    while k < start.len() && (start[k].is_ascii_alphanumeric() || start[k] == b'_')
                    {
                        k += 1;
                    }
                    if k < start.len() && start[k] == b':' {
                        *input = &start[i..];
                        return Ok(&start[..i]);
                    }
                }
            }
            _ => {}
        }
        i += 1;
    }
    *input = &start[i..];
    Ok(&start[..i])
}

/// Parse comma only if not followed by an entry-start character.
fn comma_not_before_entry(input: &mut &[u8]) -> PResult<u8> {
    terminated(b',', not(one_of(is_cx_tag_start))).parse_next(input)
}

/// Skip over an unknown/unrecognized CX entry.
///
/// CXSMILES uses commas as both list separators *within* an entry and as entry separators.
/// We stop at:
/// - the closing `|`, or
/// - a comma that is followed by the start of another CX entry.
fn skip_unknown_entry(input: &mut &[u8]) -> PResult<()> {
    if input.is_empty() {
        return Ok(());
    }

    let start = *input;
    let mut i = 0usize;
    while i < start.len() {
        if start[i] == b',' {
            // A comma starts a new entry iff the next non-whitespace char looks like an entry start.
            let mut j = i + 1;
            while j < start.len() && start[j].is_ascii_whitespace() {
                j += 1;
            }
            if j < start.len() && is_cx_tag_start(start[j]) {
                break;
            }
        }
        i += 1;
    }

    if i == 0 {
        return Err(ErrMode::Backtrack(CxParseError::Syntax));
    }
    *input = &start[i..];
    Ok(())
}

/// Parse a basic CX entry (rejects extended features).
fn parse_unknown_entry(input: &mut &[u8], skip_unknown_cx_tags: bool) -> PResult<Option<CxEntry>> {
    if skip_unknown_cx_tags {
        skip_unknown_entry(input)?;
        Ok(None)
    } else {
        Err(ErrMode::Cut(CxParseError::InvalidTag))
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;
    use crate::table_ir::atom::{BicycloStereo, BicycloStereoData, Chirality};
    use crate::table_ir::{Atom, Bond, ExtendedAtom, ExtendedBond};

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
        let mut remaining = input;
        let result = parse_coordinates(&mut remaining);
        assert_eq!(result, Ok(expected));
        assert_eq!(remaining, b"");
    }

    #[rstest]
    #[case::atom_4d(b"(1.0,2.0,3.0,4.0)")]
    fn test_parse_coordinates_invalid(#[case] input: &[u8]) {
        let mut remaining = input;
        let result = parse_coordinates(&mut remaining);
        assert_eq!(result, Err(ErrMode::Cut(CxParseError::Syntax)));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(b"||", vec![])]
    #[case::coordinate_bond(b"|C:0.1|", vec![CxEntry::CoordinateBonds(vec![(0, 1)])])]
    #[case::coordinate_bond_multiple(b"|C:0.1,2.3|", vec![CxEntry::CoordinateBonds(vec![(0, 1), (2, 3)])])]
    #[case::hydrogen_bond(b"|H:1.2|", vec![CxEntry::HydrogenBonds(vec![(1, 2)])])]
    #[case::radicals_multiple_atoms(b"|^1:0,1,2|", vec![CxEntry::Radicals(vec![(0, (1, None)),
        (1, (1, None)), (2, (1, None))])])]
    #[case::wiggly_bonds(b"|w:0.1,2.3|", vec![CxEntry::WigglyBonds(vec![(0, 1, BondWedge::Either), (2, 3, BondWedge::Either)])])]
    #[case::cis_bonds(b"|c:0,1|", vec![CxEntry::CisBonds(vec![0, 1])])]
    #[case::trans_bonds(b"|t:0,1|", vec![CxEntry::TransBonds(vec![0, 1])])]
    #[case::unspec_bonds(b"|ctu:0,1|", vec![CxEntry::UnspecBonds(vec![0, 1])])]
    #[case::multicenter_bonds(b"|m:0:3.4,2:1.5|", vec![CxEntry::MulticenterBonds(vec![(0, vec![3, 4]), (2, vec![1, 5])])])]
    #[case::atom_labels(b"|$label1;label2;label3$|", vec![CxEntry::Labels(vec![(0, "label1".to_string()), (1, "label2".to_string()), (2, "label3".to_string())])])]
    #[case::atom_values(b"$_AV:value1;value2;value3$|", vec![CxEntry::Values(vec![(0, "value1".to_string()), (1, "value2".to_string()), (2, "value3".to_string())])])]
    #[case::coordinates_2d(b"|(1.5,2.5;3.5,4.5)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.5, 2.5, 0.0), Point3D::new(3.5, 4.5, 0.0)])])]
    #[case::coordinates_3d(b"|(1,2,3;4,5,6)|", vec![CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 3.0), Point3D::new(4.0, 5.0, 6.0)])])]
    #[case::combined_entries(b"|^1:0,1,(1.0,2.0;3.0,4.0),C:2.3|", vec![CxEntry::Radicals(vec![(0, (1, None)),
        (1, (1, None))]), CxEntry::Coordinates(vec![Point3D::new(1.0, 2.0, 0.0), Point3D::new(3.0, 4.0, 0.0)]), CxEntry::CoordinateBonds(vec![(2, 3)])])]
    fn test_parse_cx_annotations(#[case] input: &[u8], #[case] expected: Vec<CxEntry>) {
        let result = parse_cx_annotations(input, SmilesSyntaxFlags::default());
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
        let result = parse_cx_annotations(input, SmilesSyntaxFlags::default());
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
    fn test_parse_cx_annotations_unknown_tags(
        #[case] input: &[u8],
        #[case] expected: Vec<CxEntry>,
    ) {
        let flags = SmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS;
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
    #[case::radicals(b"|^1:0|", vec![CxEntry::Radicals(vec![(0, (1, None))])])]
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
        let result = parse_extended_cx_annotations(input, SmilesSyntaxFlags::default());
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
        let result = parse_extended_cx_annotations(input, SmilesSyntaxFlags::default());
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
    fn test_parse_extended_cx_annotations_unknown_tags(
        #[case] input: &[u8],
        #[case] expected: Vec<CxEntry>,
    ) {
        let flags = SmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS;
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
    #[case::relative_stereo_with_fragment_list(
        b"|r:0|",
        SmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS
    )]
    fn test_parse_extended_cx_annotations_unknown_tags_error(
        #[case] input: &[u8],
        #[case] flags: SmilesSyntaxFlags,
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
    #[case::radicals(vec![CxEntry::Radicals(vec![(0, (1, None))])],
        |mol: &Molecule| mol.atoms[0].unpaired_electrons == Some(1) && mol.atoms[0].multiplicity.is_none())]
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
    #[case::radicals(vec![CxEntry::Radicals(vec![(2, (2, None))])],
        |mol: &ExtendedMolecule| mol.atoms[2].unpaired_electrons == Some(2) && mol.atoms[2].multiplicity.is_none())]
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
        |mol: &ExtendedMolecule| mol.configuration_scope == Some(ConfigurationScope::Absolute) && mol.cx_data.is_none())]
    #[case::stereo_group_or(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::Or(1), atoms: vec![0] })],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_groups.get(&1)) == Some(&StereoSet { atoms: vec![0], relation: StereoSetRelation::Correlated }))]
    #[case::stereo_group_and(vec![CxEntry::StereoGroup(StereoGroup { group_type: StereoGroupType::And(2), atoms: vec![1] })],
        |mol: &ExtendedMolecule| mol.cx_data.as_ref().and_then(|d| d.stereo_groups.get(&2)) == Some(&StereoSet { atoms: vec![1], relation: StereoSetRelation::Independent }))]
    #[case::relative_stereo(vec![CxEntry::RelativeStereo],
        |mol: &ExtendedMolecule| mol.configuration_scope == Some(ConfigurationScope::Relative) && mol.cx_data.is_none())]
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
