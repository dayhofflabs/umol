//! CXSMILES annotation block parser
//!
//! Parses the `|...|` extension block in CXSMILES format.
//! Two parsers are provided:
//! - `parse_cx_annotations`: basic annotations only (for Molecule)
//! - `parse_extended_cx_annotations`: all annotations (for ExtendedMolecule)

// TODO: Implement the parsers
#![allow(dead_code)]

use std::collections::{BTreeMap, HashMap};

use umol_data::SpinMultiplicity;

use super::super::error::ParseError;
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

/// Accumulator for CX annotation properties
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CxAccumulator {
    // Per-atom (indexed by atom position in SMILES)
    pub atom_labels: BTreeMap<u32, String>,
    pub atom_values: BTreeMap<u32, String>,
    pub atom_unpaired_electrons: BTreeMap<u32, UnpairedElectrons>,
    pub atom_properties: BTreeMap<u32, HashMap<String, String>>,

    // Per-bond (indexed by atom index for wiggly, or atom pairs for others)
    pub wiggly_bonds: BTreeMap<u32, BondWedge>,
    pub coordinate_bonds: Vec<(u32, u32)>,
    pub hydrogen_bonds: Vec<(u32, u32)>,
    pub cis_bonds: Vec<u32>,
    pub trans_bonds: Vec<u32>,

    // Molecule-level
    pub coordinates: Option<Vec<Point3D>>,
    pub stereo_groups: Vec<StereoGroup>,
    pub relative_stereo: bool,
    pub fragment_groups: Vec<Vec<u32>>,
}

impl CxAccumulator {
    pub fn new() -> Self {
        Self::default()
    }
}

/// Parse basic CX annotations (for Molecule)
///
/// Handles only basic annotations:
/// - Coordinates `(x,y,z;...)`
/// - Radicals `^1:` through `^7:`
/// - Atom labels `$...;...$`
/// - Atom values `$_AV:...$`
/// - Wiggly bonds `w:`, `wU:`, `wD:`
/// - CIS/TRANS `c:`, `t:`, `ctu:`
///
/// Returns error if extended-only annotations are encountered.
pub fn parse_cx_annotations(input: &[u8]) -> Result<CxAccumulator, ParseError> {
    let mut acc = CxAccumulator::new();

    // Strip leading/trailing pipes
    let inner = strip_pipes(input)?;
    if inner.is_empty() {
        return Ok(acc);
    }

    // Tokenize on ',' respecting nested structures
    let tokens = tokenize(inner)?;

    for (token, pos) in tokens {
        parse_basic_tag(&mut acc, token, pos)?;
    }

    Ok(acc)
}

/// Parse extended CX annotations (for ExtendedMolecule)
///
/// Handles all annotations including:
/// - All basic annotations
/// - Fragment grouping `f:`
/// - Enhanced stereo `a:`, `o<n>:`, `&<n>:`
/// - Atom properties `atomProp:`
/// - Coordinate bonds `C:`
/// - Hydrogen bonds `H:`
/// - Relative stereo `r`
pub fn parse_extended_cx_annotations(input: &[u8]) -> Result<CxAccumulator, ParseError> {
    let mut acc = CxAccumulator::new();

    // Strip leading/trailing pipes
    let inner = strip_pipes(input)?;
    if inner.is_empty() {
        return Ok(acc);
    }

    // Tokenize on ',' respecting nested structures
    let tokens = tokenize(inner)?;

    for (token, pos) in tokens {
        parse_extended_tag(&mut acc, token, pos)?;
    }

    Ok(acc)
}

/// Strip leading `|` and trailing `|` from input
fn strip_pipes(input: &[u8]) -> Result<&[u8], ParseError> {
    if input.is_empty() {
        return Ok(input);
    }

    let start = if input.first() == Some(&b'|') { 1 } else { 0 };
    let end = if input.last() == Some(&b'|') {
        input.len() - 1
    } else {
        input.len()
    };

    if start > end {
        return Ok(&[]);
    }

    Ok(&input[start..end])
}

/// Tokenize on ',' respecting nested `()`, `{}`, `$...$`
/// Returns (token, position) pairs
fn tokenize(input: &[u8]) -> Result<Vec<(&[u8], usize)>, ParseError> {
    let mut tokens = Vec::new();
    let mut start = 0;
    let mut i = 0;
    let mut paren_depth: u32 = 0;
    let mut brace_depth: u32 = 0;
    let mut in_dollar = false;

    while i < input.len() {
        let b = input[i];

        match b {
            b'(' if !in_dollar => paren_depth += 1,
            b')' if !in_dollar => paren_depth = paren_depth.saturating_sub(1),
            b'{' if !in_dollar => brace_depth += 1,
            b'}' if !in_dollar => brace_depth = brace_depth.saturating_sub(1),
            b'$' => in_dollar = !in_dollar,
            b',' if paren_depth == 0 && brace_depth == 0 && !in_dollar => {
                if i > start {
                    tokens.push((&input[start..i], start));
                }
                start = i + 1;
            }
            _ => {}
        }

        i += 1;
    }

    // Final token
    if start < input.len() {
        tokens.push((&input[start..], start));
    }

    Ok(tokens)
}

/// Unescape `&#code;` sequences in a string
fn unescape(input: &[u8]) -> Vec<u8> {
    let mut result = Vec::with_capacity(input.len());
    let mut i = 0;

    while i < input.len() {
        if i + 2 < input.len() && input[i] == b'&' && input[i + 1] == b'#' {
            // Find the semicolon
            let mut j = i + 2;
            while j < input.len() && input[j] != b';' && input[j].is_ascii_digit() {
                j += 1;
            }
            if j < input.len() && input[j] == b';' {
                // Parse the code
                if let Ok(s) = std::str::from_utf8(&input[i + 2..j]) {
                    if let Ok(code) = s.parse::<u8>() {
                        result.push(code);
                        i = j + 1;
                        continue;
                    }
                }
            }
        }
        result.push(input[i]);
        i += 1;
    }

    result
}

/// Parse a basic tag (errors on extended-only tags)
fn parse_basic_tag(acc: &mut CxAccumulator, token: &[u8], pos: usize) -> Result<(), ParseError> {
    // Coordinates: starts with '('
    if token.starts_with(b"(") {
        return parse_coordinates(acc, token, pos);
    }

    // Atom labels/values: starts with '$'
    if token.starts_with(b"$") {
        return parse_atom_labels_or_values(acc, token, pos);
    }

    // Radicals: ^1: through ^7:
    if token.starts_with(b"^") {
        return parse_radicals(acc, token, pos);
    }

    // Wiggly bonds: w:, wU:, wD:
    if token.starts_with(b"w:") || token.starts_with(b"wU:") || token.starts_with(b"wD:") {
        return parse_wiggly_bonds(acc, token, pos);
    }

    // CIS/TRANS: c:, t:, ctu:
    if token.starts_with(b"c:") || token.starts_with(b"t:") || token.starts_with(b"ctu:") {
        return parse_cis_trans(acc, token, pos);
    }

    // Extended-only tags - error in basic parser
    if token.starts_with(b"f:")
        || token.starts_with(b"a:")
        || token.starts_with(b"o")
        || token.starts_with(b"&")
        || token.starts_with(b"atomProp:")
        || token.starts_with(b"C:")
        || token.starts_with(b"H:")
        || token == b"r"
        || token.starts_with(b"r:")
        || token.starts_with(b"Sg:")
        || token.starts_with(b"RG:")
        || token.starts_with(b"rb:")
        || token.starts_with(b"s:")
        || token.starts_with(b"u:")
        || token.starts_with(b"LN:")
        || token.starts_with(b"LO:")
    {
        return Err(ParseError::InvalidCxProperty { pos });
    }

    // Unknown tag - ignore for forward compatibility
    Ok(())
}

/// Parse an extended tag (handles all tags)
fn parse_extended_tag(acc: &mut CxAccumulator, token: &[u8], pos: usize) -> Result<(), ParseError> {
    // Coordinates: starts with '('
    if token.starts_with(b"(") {
        return parse_coordinates(acc, token, pos);
    }

    // Atom labels/values: starts with '$'
    if token.starts_with(b"$") {
        return parse_atom_labels_or_values(acc, token, pos);
    }

    // Radicals: ^1: through ^7:
    if token.starts_with(b"^") {
        return parse_radicals(acc, token, pos);
    }

    // Wiggly bonds: w:, wU:, wD:
    if token.starts_with(b"w:") || token.starts_with(b"wU:") || token.starts_with(b"wD:") {
        return parse_wiggly_bonds(acc, token, pos);
    }

    // CIS/TRANS: c:, t:, ctu:
    if token.starts_with(b"c:") || token.starts_with(b"t:") || token.starts_with(b"ctu:") {
        return parse_cis_trans(acc, token, pos);
    }

    // Fragment grouping: f:
    if token.starts_with(b"f:") {
        return parse_fragment_groups(acc, token, pos);
    }

    // Enhanced stereo: a:, o<n>:, &<n>:
    if token.starts_with(b"a:") {
        return parse_stereo_absolute(acc, token, pos);
    }
    if token.starts_with(b"o") || token.starts_with(b"&") {
        return parse_stereo_group(acc, token, pos);
    }

    // Atom properties: atomProp:
    if token.starts_with(b"atomProp:") {
        return parse_atom_properties(acc, token, pos);
    }

    // Coordinate bonds: C:
    if token.starts_with(b"C:") {
        return parse_coordinate_bonds(acc, token, pos);
    }

    // Hydrogen bonds: H:
    if token.starts_with(b"H:") {
        return parse_hydrogen_bonds(acc, token, pos);
    }

    // Relative stereo: r or r:idx,...
    if token == b"r" || token.starts_with(b"r:") {
        return parse_relative_stereo(acc, token, pos);
    }

    // Deferred CTab legacy - ignore for now
    // Sg:, RG:, rb:, s:, u:, LN:, LO:, etc.

    // Unknown tag - ignore for forward compatibility
    Ok(())
}

/// Parse coordinates
fn parse_coordinates(acc: &mut CxAccumulator, token: &[u8], pos: usize) -> Result<(), ParseError> {
    // Format: (x,y,z;x2,y2,z2;...) - semicolon-separated triplets
    // Empty entries (missing coordinates) are just empty between semicolons

    // Strip outer parentheses
    if token.len() < 2 || token[0] != b'(' || token[token.len() - 1] != b')' {
        return Err(ParseError::InvalidToken { pos });
    }
    let inner = &token[1..token.len() - 1];

    if inner.is_empty() {
        acc.coordinates = Some(Vec::new());
        return Ok(());
    }

    let mut coords = Vec::new();

    // Split on semicolons
    for (entry_idx, entry) in inner.split(|&b| b == b';').enumerate() {
        if entry.is_empty() {
            // Empty entry - use NaN to indicate missing coordinates
            coords.push(Point3D::new(f64::NAN, f64::NAN, f64::NAN));
            continue;
        }

        // Split on commas to get x,y,z
        let parts: Vec<&[u8]> = entry.split(|&b| b == b',').collect();

        // Need exactly 3 parts for x,y,z (though z may be empty for 2D)
        if parts.is_empty() || parts.len() > 3 {
            return Err(ParseError::InvalidToken {
                pos: pos + 1 + entry_idx,
            });
        }

        let x = parts
            .first()
            .map_or(Ok(0.0), |bytes| parse_f64(bytes, pos))?;
        let y = parts
            .get(1)
            .map_or(Ok(0.0), |bytes| parse_f64(bytes, pos))?;
        let z = parts
            .get(2)
            .map_or(Ok(0.0), |bytes| parse_f64(bytes, pos))?;

        coords.push(Point3D::new(x, y, z));
    }

    acc.coordinates = Some(coords);
    Ok(())
}

/// Parse floating-point coordinate value
#[inline]
fn parse_f64(bytes: &[u8], pos: usize) -> Result<f64, ParseError> {
    // Default to 0 for missing component (e.g., z in 2D)
    if bytes.is_empty() {
        return Ok(0.0);
    }
    fast_float2::parse::<f64, _>(bytes).map_err(|_| ParseError::InvalidToken { pos })
}

/// Parse atom labels or values
///
/// Format: `$label1;label2;...$` for labels, or `$_AV:val1;val2;...$` for atom values
/// Empty entries (;;) mean no label/value for that atom position
///
/// Note: Escaped characters use `&#n;` format where n is ASCII code.
/// We must parse carefully since `;` is both a separator and part of escape sequences.
fn parse_atom_labels_or_values(
    acc: &mut CxAccumulator,
    token: &[u8],
    _pos: usize,
) -> Result<(), ParseError> {
    // Must start and end with $
    if token.len() < 2 || token[0] != b'$' || token[token.len() - 1] != b'$' {
        return Ok(()); // Not a valid label/value block, ignore
    }

    let inner = &token[1..token.len() - 1];
    if inner.is_empty() {
        return Ok(());
    }

    // Check if this is atom values ($_AV:...)
    let (is_values, data) = if inner.starts_with(b"_AV:") {
        (true, &inner[4..])
    } else {
        (false, inner)
    };

    // Split on semicolons, but respect escape sequences &#n;
    let entries = split_escaped_semicolons(data);

    for (idx, entry) in entries.into_iter().enumerate() {
        if !entry.is_empty() {
            let unescaped = unescape(&entry);
            if let Ok(s) = String::from_utf8(unescaped) {
                if is_values {
                    acc.atom_values.insert(idx as u32, s);
                } else {
                    acc.atom_labels.insert(idx as u32, s);
                }
            }
        }
    }

    Ok(())
}

/// Split on semicolons while respecting `&#n;` escape sequences
fn split_escaped_semicolons(input: &[u8]) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    let mut current = Vec::new();
    let mut i = 0;

    while i < input.len() {
        // Check for escape sequence &#...;
        if i + 2 < input.len() && input[i] == b'&' && input[i + 1] == b'#' {
            // Find the end of the escape sequence
            let mut j = i + 2;
            while j < input.len() && input[j].is_ascii_digit() {
                j += 1;
            }
            if j < input.len() && input[j] == b';' {
                // This semicolon is part of an escape sequence, include it
                current.extend_from_slice(&input[i..=j]);
                i = j + 1;
                continue;
            }
        }

        if input[i] == b';' {
            // Separator semicolon
            result.push(std::mem::take(&mut current));
            i += 1;
        } else {
            current.push(input[i]);
            i += 1;
        }
    }

    result.push(current);
    result
}

/// Convert CXSMILES radical code (1-7) to UnpairedElectrons
///
/// ChemAxon CXSMILES radical encoding:
/// - ^1: monovalent radical (1 electron)
/// - ^2: divalent radical (2 electrons)
/// - ^3: divalent singlet radical (2 electrons, singlet spin)
/// - ^4: divalent triplet radical (2 electrons, triplet spin)
/// - ^5: trivalent radical (3 electrons)
/// - ^6: trivalent doublet radical (3 electrons, doublet spin)
/// - ^7: trivalent quartet radical (3 electrons, quartet spin)
fn convert_radical_code(code: u8) -> UnpairedElectrons {
    match code {
        1 => UnpairedElectrons::from_count(1),
        2 => UnpairedElectrons::from_count(2),
        3 => UnpairedElectrons::new(2, Some(SpinMultiplicity::Singlet)),
        4 => UnpairedElectrons::new(2, Some(SpinMultiplicity::Triplet)),
        5 => UnpairedElectrons::from_count(3),
        6 => UnpairedElectrons::new(3, Some(SpinMultiplicity::Doublet)),
        7 => UnpairedElectrons::new(3, Some(SpinMultiplicity::Quartet)),
        _ => unreachable!("invalid radical code: {}", code),
    }
}

/// Parse radical annotations
fn parse_radicals(acc: &mut CxAccumulator, token: &[u8], pos: usize) -> Result<(), ParseError> {
    // Format: ^n:idx1,idx2,... where n is 1-7
    // Can have multiple radical specs: ^1:0,1^2:2,3
    let mut i = 0;

    while i < token.len() {
        // Expect '^'
        if token[i] != b'^' {
            return Err(ParseError::InvalidToken { pos: pos + i });
        }
        i += 1;

        // Expect digit 1-7
        if i >= token.len() {
            return Err(ParseError::InvalidToken { pos: pos + i });
        }
        let radical_type = token[i];
        if !(b'1'..=b'7').contains(&radical_type) {
            return Err(ParseError::InvalidToken { pos: pos + i });
        }
        let radical_code = radical_type - b'0';
        let unpaired = convert_radical_code(radical_code);
        i += 1;

        // Expect ':'
        if i >= token.len() || token[i] != b':' {
            return Err(ParseError::InvalidToken { pos: pos + i });
        }
        i += 1;

        // Parse comma-separated atom indices
        loop {
            let (idx, consumed) = parse_u32(&token[i..], pos + i)?;
            if consumed == 0 {
                // No digits found - might be end or next radical spec
                break;
            }
            i += consumed;

            acc.atom_unpaired_electrons.insert(idx, unpaired);

            // Check for comma (more indices) or end
            if i < token.len() && token[i] == b',' {
                i += 1;
                // If next char is '^', break to parse next radical spec
                if i < token.len() && token[i] == b'^' {
                    break;
                }
            } else {
                break;
            }
        }
    }

    Ok(())
}

/// Parse unsigned 32-bit integer
/// Returns (value, bytes_consumed). Returns (0, 0) if no digits found.
#[inline]
fn parse_u32(bytes: &[u8], pos: usize) -> Result<(u32, usize), ParseError> {
    let mut i = 0;
    let mut value: u32 = 0;

    while i < bytes.len() && bytes[i].is_ascii_digit() {
        let digit = (bytes[i] - b'0') as u32;
        value = value
            .checked_mul(10)
            .and_then(|v| v.checked_add(digit))
            .ok_or(ParseError::InvalidToken { pos })?;
        i += 1;
    }

    Ok((value, i))
}

/// Parse wiggly bond annotations
///
/// Format: `w:atomIdx.bondIdx,...` for wiggly (either), `wU:...` for up, `wD:...` for down
fn parse_wiggly_bonds(acc: &mut CxAccumulator, token: &[u8], pos: usize) -> Result<(), ParseError> {
    let (wedge_type, rest) = if token.starts_with(b"wU:") {
        (BondWedge::Up, &token[3..])
    } else if token.starts_with(b"wD:") {
        (BondWedge::Down, &token[3..])
    } else if token.starts_with(b"w:") {
        (BondWedge::Either, &token[2..])
    } else {
        return Ok(());
    };

    // Parse comma-separated atom.bond pairs
    for entry in rest.split(|&b| b == b',') {
        if entry.is_empty() {
            continue;
        }

        // Find the dot separator
        let dot_pos = entry.iter().position(|&b| b == b'.');
        let atom_bytes = dot_pos.map_or(entry, |p| &entry[..p]);

        let (atom_idx, consumed) = parse_u32(atom_bytes, pos)?;
        if consumed == 0 {
            continue;
        }

        acc.wiggly_bonds.insert(atom_idx, wedge_type);
    }

    Ok(())
}

/// Parse cis/trans bond annotations
///
/// Format: `c:bondIdx,...` for cis, `t:bondIdx,...` for trans, `ctu:...` for unspecified
fn parse_cis_trans(acc: &mut CxAccumulator, token: &[u8], pos: usize) -> Result<(), ParseError> {
    let (is_cis, rest) = if token.starts_with(b"ctu:") {
        // Unspecified - we ignore these for now
        return Ok(());
    } else if token.starts_with(b"c:") {
        (true, &token[2..])
    } else if token.starts_with(b"t:") {
        (false, &token[2..])
    } else {
        return Ok(());
    };

    // Parse comma-separated bond indices
    for entry in rest.split(|&b| b == b',') {
        if entry.is_empty() {
            continue;
        }

        let (bond_idx, consumed) = parse_u32(entry, pos)?;
        if consumed == 0 {
            continue;
        }

        if is_cis {
            acc.cis_bonds.push(bond_idx);
        } else {
            acc.trans_bonds.push(bond_idx);
        }
    }

    Ok(())
}

/// Parse fragment grouping
///
/// Format: `f:idx.idx,idx.idx,...` - dot-separated indices form a group, comma separates groups
fn parse_fragment_groups(
    acc: &mut CxAccumulator,
    token: &[u8],
    pos: usize,
) -> Result<(), ParseError> {
    if !token.starts_with(b"f:") {
        return Ok(());
    }
    let rest = &token[2..];

    // Parse comma-separated groups
    for group_bytes in rest.split(|&b| b == b',') {
        if group_bytes.is_empty() {
            continue;
        }

        let mut group = Vec::new();

        // Parse dot-separated fragment indices within each group
        for idx_bytes in group_bytes.split(|&b| b == b'.') {
            if idx_bytes.is_empty() {
                continue;
            }
            let (idx, consumed) = parse_u32(idx_bytes, pos)?;
            if consumed > 0 {
                group.push(idx);
            }
        }

        if !group.is_empty() {
            acc.fragment_groups.push(group);
        }
    }

    Ok(())
}

/// Parse absolute stereo group
///
/// Format: `a:idx,idx,...`
fn parse_stereo_absolute(
    acc: &mut CxAccumulator,
    token: &[u8],
    pos: usize,
) -> Result<(), ParseError> {
    if !token.starts_with(b"a:") {
        return Ok(());
    }
    let rest = &token[2..];

    let mut atoms = Vec::new();
    for idx_bytes in rest.split(|&b| b == b',') {
        if idx_bytes.is_empty() {
            continue;
        }
        let (idx, consumed) = parse_u32(idx_bytes, pos)?;
        if consumed > 0 {
            atoms.push(idx);
        }
    }

    if !atoms.is_empty() {
        acc.stereo_groups.push(StereoGroup {
            group_type: StereoGroupType::Absolute,
            atoms,
        });
    }

    Ok(())
}

/// Parse OR or AND stereo group
///
/// Format: `o<n>:idx,idx,...` for OR groups, `&<n>:idx,idx,...` for AND groups
fn parse_stereo_group(acc: &mut CxAccumulator, token: &[u8], pos: usize) -> Result<(), ParseError> {
    let is_or = token.starts_with(b"o");
    let is_and = token.starts_with(b"&");

    if !is_or && !is_and {
        return Ok(());
    }

    // Skip the initial character (o or &)
    let rest = &token[1..];

    // Parse group number until ':'
    let colon_pos = rest.iter().position(|&b| b == b':');
    let (group_bytes, atom_bytes) = match colon_pos {
        Some(p) => (&rest[..p], &rest[p + 1..]),
        None => return Ok(()), // No colon found, invalid
    };

    let (group_num, consumed) = parse_u32(group_bytes, pos)?;
    if consumed == 0 {
        return Ok(()); // No group number
    }

    // Parse atom indices
    let mut atoms = Vec::new();
    for idx_bytes in atom_bytes.split(|&b| b == b',') {
        if idx_bytes.is_empty() {
            continue;
        }
        let (idx, consumed) = parse_u32(idx_bytes, pos)?;
        if consumed > 0 {
            atoms.push(idx);
        }
    }

    if !atoms.is_empty() {
        let group_type = if is_or {
            StereoGroupType::Or(group_num)
        } else {
            StereoGroupType::And(group_num)
        };
        acc.stereo_groups.push(StereoGroup { group_type, atoms });
    }

    Ok(())
}

/// Parse atom properties
///
/// Format: `atomProp:idx.key.value:idx.key.value:...`
fn parse_atom_properties(
    acc: &mut CxAccumulator,
    token: &[u8],
    pos: usize,
) -> Result<(), ParseError> {
    if !token.starts_with(b"atomProp:") {
        return Ok(());
    }
    let rest = &token[9..];

    // Split on colons for each property entry
    for entry in rest.split(|&b| b == b':') {
        if entry.is_empty() {
            continue;
        }

        // Split on dots: idx.key.value
        let parts: Vec<&[u8]> = entry.splitn(3, |&b| b == b'.').collect();
        if parts.len() != 3 {
            continue; // Invalid entry, skip
        }

        let (idx, consumed) = parse_u32(parts[0], pos)?;
        if consumed == 0 {
            continue;
        }

        let key_unescaped = unescape(parts[1]);
        let value_unescaped = unescape(parts[2]);

        let key = match String::from_utf8(key_unescaped) {
            Ok(s) => s,
            Err(_) => continue,
        };
        let value = match String::from_utf8(value_unescaped) {
            Ok(s) => s,
            Err(_) => continue,
        };

        acc.atom_properties
            .entry(idx)
            .or_default()
            .insert(key, value);
    }

    Ok(())
}

/// Parse coordinate bonds
///
/// Format: `C:atomIdx.bondIdx,...`
fn parse_coordinate_bonds(
    acc: &mut CxAccumulator,
    token: &[u8],
    pos: usize,
) -> Result<(), ParseError> {
    if !token.starts_with(b"C:") {
        return Ok(());
    }
    let rest = &token[2..];

    for entry in rest.split(|&b| b == b',') {
        if entry.is_empty() {
            continue;
        }

        // Split on dot: atomIdx.bondIdx
        let dot_pos = entry.iter().position(|&b| b == b'.');
        let (atom_bytes, bond_bytes) = match dot_pos {
            Some(p) => (&entry[..p], &entry[p + 1..]),
            None => continue,
        };

        let (atom_idx, consumed1) = parse_u32(atom_bytes, pos)?;
        if consumed1 == 0 {
            continue;
        }

        let (bond_idx, consumed2) = parse_u32(bond_bytes, pos)?;
        if consumed2 == 0 {
            continue;
        }

        acc.coordinate_bonds.push((atom_idx, bond_idx));
    }

    Ok(())
}

/// Parse hydrogen bonds
///
/// Format: `H:atomIdx.bondIdx,...`
fn parse_hydrogen_bonds(
    acc: &mut CxAccumulator,
    token: &[u8],
    pos: usize,
) -> Result<(), ParseError> {
    if !token.starts_with(b"H:") {
        return Ok(());
    }
    let rest = &token[2..];

    for entry in rest.split(|&b| b == b',') {
        if entry.is_empty() {
            continue;
        }

        // Split on dot: atomIdx.bondIdx
        let dot_pos = entry.iter().position(|&b| b == b'.');
        let (atom_bytes, bond_bytes) = match dot_pos {
            Some(p) => (&entry[..p], &entry[p + 1..]),
            None => continue,
        };

        let (atom_idx, consumed1) = parse_u32(atom_bytes, pos)?;
        if consumed1 == 0 {
            continue;
        }

        let (bond_idx, consumed2) = parse_u32(bond_bytes, pos)?;
        if consumed2 == 0 {
            continue;
        }

        acc.hydrogen_bonds.push((atom_idx, bond_idx));
    }

    Ok(())
}

/// Parse relative stereo marker
///
/// Format: `r` (molecule-level) or `r:idx,idx,...` (fragment indices with relative stereo)
fn parse_relative_stereo(
    acc: &mut CxAccumulator,
    token: &[u8],
    _pos: usize,
) -> Result<(), ParseError> {
    if token == b"r" {
        // Molecule-level relative stereo
        acc.relative_stereo = true;
    }
    // Note: r:idx,idx,... format specifies which fragments have relative stereo
    // in a reaction context. For now we just set the flag.
    else if token.starts_with(b"r:") {
        acc.relative_stereo = true;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;
    use map_macro::*;
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::empty(b"||", CxAccumulator::new())]
    fn test_parse_cx_annotations(#[case] input: &[u8], #[case] expected: CxAccumulator) {
        let result = parse_cx_annotations(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_ok(), "{:?} should have succeeded", input_str);
        let acc = result.unwrap();
        assert_eq!(acc, expected);
    }

    #[rstest]
    #[case::invalid_coordinate(b"|(abc)|", ParseError::InvalidToken { pos: 0 })]
    #[case::too_many_coordinate_components(b"|(1,2,3,4)|", ParseError::InvalidToken { pos: 1 })]
    #[case::extended_feature(b"|f:0.1|", ParseError::InvalidCxProperty { pos: 0 })]
    fn test_parse_cx_annotations_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let result = parse_cx_annotations(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_err(), "{:?} should have failed", input_str);
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::empty(b"||", CxAccumulator::new())]
    fn test_parse_extended_cx_annotations(#[case] input: &[u8], #[case] expected: CxAccumulator) {
        let result = parse_extended_cx_annotations(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_ok(), "{:?} should have succeeded", input_str);
        let acc = result.unwrap();
        assert_eq!(acc, expected);
    }

    #[rstest]
    #[case::invalid_coordinate(b"|(abc)|", ParseError::InvalidToken { pos: 0 })]
    #[case::too_many_coordinate_components(b"|(1,2,3,4)|", ParseError::InvalidToken { pos: 1 })]
    fn test_parse_extended_cx_annotations_invalid(
        #[case] input: &[u8],
        #[case] expected: ParseError,
    ) {
        let result = parse_extended_cx_annotations(input);
        let input_str = input.to_str_lossy();
        assert!(result.is_err(), "{:?} should have failed", input_str);
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::both_sides(b"|foo|", b"foo")]
    #[case::no_pipes(b"foo", b"foo")]
    #[case::empty(b"||", b"")]
    #[case::blank(b"", b"")]
    fn test_strip_pipes(#[case] input: &[u8], #[case] expected: &[u8]) {
        let res = strip_pipes(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_ok(), "{:?} should have succeeded", input_str);
        let content = res.unwrap();
        assert_eq!(content, expected);
    }

    #[rstest]
    #[case::list(b"a,b,c", vec![(b"a".to_vec(), 0), (b"b".to_vec(), 2), (b"c".to_vec(), 4)])]
    #[case::nested(b"(1,2;3,4),foo", vec![(b"(1,2;3,4)".to_vec(), 0), (b"foo".to_vec(), 10)])]
    #[case::dollar(b"$a;b;c$,foo", vec![(b"$a;b;c$".to_vec(), 0), (b"foo".to_vec(), 8)])]
    fn test_tokenize(#[case] input: &[u8], #[case] expected_tokens: Vec<(Vec<u8>, usize)>) {
        let res = tokenize(input);
        let input_str = input.to_str_lossy();
        assert!(res.is_ok(), "{:?} should have succeeded", input_str);
        let tokens: Vec<(Vec<u8>, usize)> = res
            .unwrap()
            .into_iter()
            .map(|(b, pos)| (b.to_vec(), pos))
            .collect();
        assert_eq!(tokens, expected_tokens);
    }

    #[rstest]
    #[case::semicolon(b"foo&#59;bar", b"foo;bar")]
    #[case::comma(b"a&#44;b", b"a,b")]
    #[case::no_escape(b"plain", b"plain")]
    fn test_unescape(#[case] input: &[u8], #[case] expected: &[u8]) {
        let res = unescape(input);
        assert_eq!(res, expected);
    }

    #[rstest]
    #[case::empty(b"()", vec![])]
    #[case::atom_3d(b"(1.0,2.0,3.0)", vec![Point3D::new(1.0, 2.0, 3.0)])]
    #[case::diatomic_3d(b"(1,2,3;4,5,6)", vec![Point3D::new(1.0, 2.0, 3.0), Point3D::new(4.0, 5.0, 6.0)])]
    #[case::atom_3d_negative(b"(1.5,-2.5,3.5)", vec![Point3D::new(1.5, -2.5, 3.5)])]
    #[case::atom_2d(b"(1,2)", vec![Point3D::new(1.0, 2.0, 0.0)])]
    #[case::atom_1d(b"(1)", vec![Point3D::new(1.0, 0.0, 0.0)])]
    fn test_parse_coordinates(#[case] input: &[u8], #[case] expected: Vec<Point3D>) {
        let mut acc = CxAccumulator::new();
        let result = parse_coordinates(&mut acc, input, 0);
        let input_str = input.to_str_lossy();
        assert!(result.is_ok(), "{:?} should have succeeded", input_str);
        assert_eq!(acc.coordinates, Some(expected));
    }

    #[test]
    fn test_parse_coordinates_empty_entries() {
        let input = b"(;)";
        let mut acc = CxAccumulator::new();
        let result = parse_coordinates(&mut acc, input, 0);
        let input_str = input.to_str_lossy();
        assert!(result.is_ok(), "{:?} should have succeeded", input_str);
        let coords = acc.coordinates.unwrap();
        assert_eq!(coords.len(), 2);
        for coord in &coords {
            assert!(coord.x.is_nan());
        }
    }

    #[rstest]
    #[case::missing_parentheses(b"1,2,3", ParseError::InvalidToken { pos: 0 })]
    #[case::invalid_number(b"(abc)", ParseError::InvalidToken { pos: 0 })]
    #[case::too_many_components(b"(1,2,3,4)", ParseError::InvalidToken { pos: 1 })]
    fn test_parse_coordinates_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_coordinates(&mut acc, input, 0);
        let input_str = input.to_str_lossy();
        assert!(result.is_err(), "{:?} should have failed", input_str);
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::single_label(b"$foo$", btree_map! { 0 => "foo".to_string() })]
    #[case::multiple_labels(b"$foo;bar;baz$", btree_map! { 0 => "foo".to_string(), 1 => "bar".to_string(), 2 => "baz".to_string() })]
    #[case::empty_entries(b"$;foo;;bar$", btree_map! { 1 => "foo".to_string(), 3 => "bar".to_string() })]
    #[case::escaped(b"$a&#59;b$", btree_map! { 0 => "a;b".to_string() })]
    fn test_parse_atom_labels(#[case] input: &[u8], #[case] expected: BTreeMap<u32, String>) {
        let mut acc = CxAccumulator::new();
        let result = parse_atom_labels_or_values(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.atom_labels, expected);
    }

    // Note: parse_atom_labels_or_values is lenient and never returns errors

    #[rstest]
    #[case::single_value(b"$_AV:val1$", btree_map! { 0 => "val1".to_string() })]
    #[case::multiple_values(b"$_AV:v1;v2;v3$", btree_map! { 0 => "v1".to_string(), 1 => "v2".to_string(), 2 => "v3".to_string() })]
    fn test_parse_atom_values(#[case] input: &[u8], #[case] expected: BTreeMap<u32, String>) {
        let mut acc = CxAccumulator::new();
        let result = parse_atom_labels_or_values(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.atom_values, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::monovalent(b"^1:0", btree_map! { 0 => UnpairedElectrons::new(1, None) })]
    #[case::multiple_atoms(b"^1:0,1,2", btree_map! { 0 => UnpairedElectrons::new(1, None), 1 => UnpairedElectrons::new(1, None), 2 => UnpairedElectrons::new(1, None) })]
    #[case::divalent(b"^2:0", btree_map! { 0 => UnpairedElectrons::new(2, None) })]
    #[case::divalent_singlet(b"^3:0", btree_map! { 0 => UnpairedElectrons::new(2, Some(SpinMultiplicity::Singlet)) })]
    #[case::divalent_triplet(b"^4:0", btree_map! { 0 => UnpairedElectrons::new(2, Some(SpinMultiplicity::Triplet)) })]
    #[case::trivalent(b"^5:0", btree_map! { 0 => UnpairedElectrons::new(3, None) })]
    #[case::trivalent_doublet(b"^6:0", btree_map! { 0 => UnpairedElectrons::new(3, Some(SpinMultiplicity::Doublet)) })]
    #[case::trivalent_quartet(b"^7:0", btree_map! { 0 => UnpairedElectrons::new(3, Some(SpinMultiplicity::Quartet)) })]
    #[case::multiple_specs(b"^1:0^2:1", btree_map! { 0 => UnpairedElectrons::new(1, None), 1 => UnpairedElectrons::new(2, None) })]
    #[case::mixed(b"^1:0,1^3:2,3", btree_map! { 0 => UnpairedElectrons::new(1, None), 1 => UnpairedElectrons::new(1, None), 2 => UnpairedElectrons::new(2, Some(SpinMultiplicity::Singlet)), 3 => UnpairedElectrons::new(2, Some(SpinMultiplicity::Singlet)) })]
    #[case::trailing_comma(b"^1:0,1,", btree_map! { 0 => UnpairedElectrons::new(1, None), 1 => UnpairedElectrons::new(1, None) })]
    fn test_parse_radicals(
        #[case] input: &[u8],
        #[case] expected: BTreeMap<u32, UnpairedElectrons>,
    ) {
        let mut acc = CxAccumulator::new();
        let result = parse_radicals(&mut acc, input, 0);
        let input_str = input.to_str_lossy();
        assert!(result.is_ok(), "{:?} should have succeeded", input_str);
        assert_eq!(acc.atom_unpaired_electrons, expected);
    }

    #[rstest]
    #[case::missing_caret(b"1:0", ParseError::InvalidToken { pos: 0 })]
    #[case::invalid_type(b"^0:0", ParseError::InvalidToken { pos: 1 })]
    #[case::invalid_type_high(b"^8:0", ParseError::InvalidToken { pos: 1 })]
    #[case::missing_colon(b"^10", ParseError::InvalidToken { pos: 2 })]
    #[case::overflow(b"^1:99999999999", ParseError::InvalidToken { pos: 3 })]
    #[case::non_numeric_trailing(b"^1:0,abc", ParseError::InvalidToken { pos: 5 })]
    #[case::empty_entry(b"^1:0,,1", ParseError::InvalidToken { pos: 5 })]
    fn test_parse_radicals_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_radicals(&mut acc, input, 0);
        let input_str = input.to_str_lossy();
        assert!(result.is_err(), "{:?} should have failed", input_str);
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::wiggly(b"w:0.1,2.3", btree_map! { 0 => BondWedge::Either, 2 => BondWedge::Either })]
    #[case::up(b"wU:0.1", btree_map! { 0 => BondWedge::Up })]
    #[case::down(b"wD:1.2", btree_map! { 1 => BondWedge::Down })]
    #[case::non_numeric_skipped(b"w:abc.1,2.3", btree_map! { 2 => BondWedge::Either })]
    #[case::empty_entry(b"w:0.1,,2.3", btree_map! { 0 => BondWedge::Either, 2 => BondWedge::Either })]
    #[case::trailing_comma(b"w:0.1,", btree_map! { 0 => BondWedge::Either })]
    #[case::missing_bond(b"w:0.,2.3", btree_map! { 0 => BondWedge::Either, 2 => BondWedge::Either })]
    #[case::missing_atom(b"w:.1,2.3", btree_map! { 2 => BondWedge::Either })]
    fn test_parse_wiggly_bonds(#[case] input: &[u8], #[case] expected: BTreeMap<u32, BondWedge>) {
        let mut acc = CxAccumulator::new();
        let result = parse_wiggly_bonds(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.wiggly_bonds, expected);
    }

    #[rstest]
    #[case::overflow(b"w:99999999999.0", ParseError::InvalidToken { pos: 0 })]
    fn test_parse_wiggly_bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_wiggly_bonds(&mut acc, input, 0);
        assert!(
            result.is_err(),
            "{:?} should have failed",
            input.to_str_lossy()
        );
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::cis(b"c:0,1,2", vec![0, 1, 2], vec![])]
    #[case::trans(b"t:3,4", vec![], vec![3, 4])]
    #[case::non_numeric_skipped(b"c:0,abc,2", vec![0, 2], vec![])]
    #[case::empty_entry(b"c:0,,2", vec![0, 2], vec![])]
    #[case::trailing_comma(b"c:0,1,", vec![0, 1], vec![])]
    #[case::leading_comma(b"c:,0,1", vec![0, 1], vec![])]
    fn test_parse_cis_trans(
        #[case] input: &[u8],
        #[case] expected_cis: Vec<u32>,
        #[case] expected_trans: Vec<u32>,
    ) {
        let mut acc = CxAccumulator::new();
        let result = parse_cis_trans(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.cis_bonds, expected_cis);
        assert_eq!(acc.trans_bonds, expected_trans);
    }

    #[rstest]
    #[case::overflow(b"c:99999999999", ParseError::InvalidToken { pos: 0 })]
    fn test_parse_cis_trans_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_cis_trans(&mut acc, input, 0);
        assert!(
            result.is_err(),
            "{:?} should have failed",
            input.to_str_lossy()
        );
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::single_group(b"f:0.1", vec![vec![0, 1]])]
    #[case::multiple_groups(b"f:0.1,2.3.4", vec![vec![0, 1], vec![2, 3, 4]])]
    #[case::non_numeric(b"f:0.abc.2,3.4", vec![vec![0, 2], vec![3, 4]])]
    #[case::empty_entry(b"f:0.1,,2.3", vec![vec![0, 1], vec![2, 3]])]
    #[case::empty_dot(b"f:0..2,3.4", vec![vec![0, 2], vec![3, 4]])]
    #[case::trailing_dot(b"f:0.1.,2.3", vec![vec![0, 1], vec![2, 3]])]
    fn test_parse_fragment_groups(#[case] input: &[u8], #[case] expected: Vec<Vec<u32>>) {
        let mut acc = CxAccumulator::new();
        let result = parse_fragment_groups(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.fragment_groups, expected);
    }

    #[rstest]
    #[case::overflow(b"f:99999999999.0", ParseError::InvalidToken { pos: 0 })]
    fn test_parse_fragment_groups_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_fragment_groups(&mut acc, input, 0);
        assert!(
            result.is_err(),
            "{:?} should have failed",
            input.to_str_lossy()
        );
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::absolute(b"a:0,1,2", vec![0, 1, 2])]
    #[case::non_numeric(b"a:0,abc,2", vec![0, 2])]
    #[case::empty_entry(b"a:0,,2", vec![0, 2])]
    #[case::trailing_comma(b"a:0,1,", vec![0, 1])]
    #[case::leading_comma(b"a:,0,1", vec![0, 1])]
    fn test_parse_stereo_absolute(#[case] input: &[u8], #[case] expected_atoms: Vec<u32>) {
        let mut acc = CxAccumulator::new();
        let result = parse_stereo_absolute(&mut acc, input, 0);
        assert!(result.is_ok());
        if expected_atoms.is_empty() {
            assert!(acc.stereo_groups.is_empty());
        } else {
            assert_eq!(acc.stereo_groups.len(), 1);
            assert_eq!(acc.stereo_groups[0].group_type, StereoGroupType::Absolute);
            assert_eq!(acc.stereo_groups[0].atoms, expected_atoms);
        }
    }

    #[rstest]
    #[case::overflow(b"a:99999999999", ParseError::InvalidToken { pos: 0 })]
    fn test_parse_stereo_absolute_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_stereo_absolute(&mut acc, input, 0);
        assert!(
            result.is_err(),
            "{:?} should have failed",
            input.to_str_lossy()
        );
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::or_group(b"o1:0,1", StereoGroupType::Or(1), vec![0, 1])]
    #[case::and_group(b"&2:3,4,5", StereoGroupType::And(2), vec![3, 4, 5])]
    #[case::non_numeric(b"o1:0,abc,2", StereoGroupType::Or(1), vec![0, 2])]
    #[case::empty_entry(b"o1:0,,2", StereoGroupType::Or(1), vec![0, 2])]
    #[case::trailing_comma(b"&1:0,1,", StereoGroupType::And(1), vec![0, 1])]
    fn test_parse_stereo_group(
        #[case] input: &[u8],
        #[case] expected_type: StereoGroupType,
        #[case] expected_atoms: Vec<u32>,
    ) {
        let mut acc = CxAccumulator::new();
        let result = parse_stereo_group(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.stereo_groups.len(), 1);
        assert_eq!(acc.stereo_groups[0].group_type, expected_type);
        assert_eq!(acc.stereo_groups[0].atoms, expected_atoms);
    }

    #[rstest]
    #[case::group_overflow(b"o99999999999:0", ParseError::InvalidToken { pos: 0 })]
    #[case::atom_overflow(b"o1:99999999999", ParseError::InvalidToken { pos: 0 })]
    fn test_parse_stereo_group_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_stereo_group(&mut acc, input, 0);
        assert!(
            result.is_err(),
            "{:?} should have failed",
            input.to_str_lossy()
        );
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::atom_property(b"atomProp:0.key.val", btree_map! { 0 => hash_map! { "key".to_string() => "val".to_string() } })]
    #[case::multiple(b"atomProp:0.k1.v1:0.k2.v2:1.k.v", btree_map! { 0 => hash_map! { "k1".to_string() => "v1".to_string(), "k2".to_string() => "v2".to_string() }, 1 => hash_map! { "k".to_string() => "v".to_string() } })]
    #[case::non_numeric(b"atomProp:abc.key.val:1.k.v", btree_map! { 1 => hash_map! { "k".to_string() => "v".to_string() } })]
    #[case::empty_entry(b"atomProp:0.k.v::1.k.v", btree_map! { 0 => hash_map! { "k".to_string() => "v".to_string() }, 1 => hash_map! { "k".to_string() => "v".to_string() } })]
    #[case::missing_parts_skipped(b"atomProp:0.key:1.k.v", btree_map! { 1 => hash_map! { "k".to_string() => "v".to_string() } })]
    fn test_parse_atom_properties(
        #[case] input: &[u8],
        #[case] expected: BTreeMap<u32, HashMap<String, String>>,
    ) {
        let mut acc = CxAccumulator::new();
        let result = parse_atom_properties(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.atom_properties, expected);
    }

    #[rstest]
    #[case::overflow(b"atomProp:99999999999.key.val", ParseError::InvalidToken { pos: 0 })]
    fn test_parse_atom_properties_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_atom_properties(&mut acc, input, 0);
        assert!(
            result.is_err(),
            "{:?} should have failed",
            input.to_str_lossy()
        );
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::coordinate_bonds(b"C:0.1,2.3", vec![(0, 1), (2, 3)])]
    #[case::non_numeric_atom(b"C:abc.1,2.3", vec![(2, 3)])]
    #[case::non_numeric_bond(b"C:0.abc,2.3", vec![(2, 3)])]
    #[case::empty_entry(b"C:0.1,,2.3", vec![(0, 1), (2, 3)])]
    #[case::missing_bond(b"C:0.,2.3", vec![(2, 3)])]
    #[case::missing_atom(b"C:.1,2.3", vec![(2, 3)])]
    fn test_parse_coordinate_bonds(#[case] input: &[u8], #[case] expected: Vec<(u32, u32)>) {
        let mut acc = CxAccumulator::new();
        let result = parse_coordinate_bonds(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.coordinate_bonds, expected);
    }

    #[rstest]
    #[case::atom_overflow(b"C:99999999999.0", ParseError::InvalidToken { pos: 0 })]
    #[case::bond_overflow(b"C:0.99999999999", ParseError::InvalidToken { pos: 0 })]
    fn test_parse_coordinate_bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_coordinate_bonds(&mut acc, input, 0);
        assert!(
            result.is_err(),
            "{:?} should have failed",
            input.to_str_lossy()
        );
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::hydrogen_bonds(b"H:3.2,5.4", vec![(3, 2), (5, 4)])]
    #[case::non_numeric_atom(b"H:abc.2,5.4", vec![(5, 4)])]
    #[case::non_numeric_bond(b"H:3.abc,5.4", vec![(5, 4)])]
    #[case::empty_entry(b"H:3.2,,5.4", vec![(3, 2), (5, 4)])]
    #[case::missing_bond(b"H:3.,5.4", vec![(5, 4)])]
    #[case::missing_atom(b"H:.2,5.4", vec![(5, 4)])]
    fn test_parse_hydrogen_bonds(#[case] input: &[u8], #[case] expected: Vec<(u32, u32)>) {
        let mut acc = CxAccumulator::new();
        let result = parse_hydrogen_bonds(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.hydrogen_bonds, expected);
    }

    #[rstest]
    #[case::atom_overflow(b"H:99999999999.0", ParseError::InvalidToken { pos: 0 })]
    #[case::bond_overflow(b"H:0.99999999999", ParseError::InvalidToken { pos: 0 })]
    fn test_parse_hydrogen_bonds_invalid(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut acc = CxAccumulator::new();
        let result = parse_hydrogen_bonds(&mut acc, input, 0);
        assert!(
            result.is_err(),
            "{:?} should have failed",
            input.to_str_lossy()
        );
        let error = result.unwrap_err();
        assert_eq!(error, expected);
    }

    #[rstest]
    #[case::simple_r(b"r", true)]
    #[case::r_with_indices(b"r:1,2", true)]
    fn test_parse_relative_stereo(#[case] input: &[u8], #[case] expected: bool) {
        let mut acc = CxAccumulator::new();
        let result = parse_relative_stereo(&mut acc, input, 0);
        assert!(result.is_ok());
        assert_eq!(acc.relative_stereo, expected);
    }

}
