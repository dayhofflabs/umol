//! SMILES parser

use umol_data::Element;

mod builder;
mod utils;

use self::builder::{AtomData, ExtendedAtomData, ExtendedMoleculeBuilder, MoleculeBuilder};
use self::utils::{
    attach_atom, attach_atom_extended, invalid_ring_context, parse_bond, parse_bracket,
    parse_bracket_extended, parse_organic_aliphatic_element, parse_organic_aromatic_element,
    parse_ring_index, process_ring_closure, process_ring_closure_extended, Frame, OpenRing,
};
use super::config::{SmilesIoConfig, SmilesParseFlags};
use super::error::ParseError;
use crate::span::Span;
use crate::table_ir::{BondDirection, BondOrder, Molecule, WildcardAtom};

/// Parse SMILES string with strict OpenSMILES rules
pub fn parse_smiles(input: &str) -> Result<Molecule, ParseError> {
    parse_smiles_bytes(input.as_bytes())
}

/// Parse SMILES string with configuration
pub fn parse_smiles_with(input: &str, config: &SmilesIoConfig) -> Result<Molecule, ParseError> {
    parse_smiles_bytes_with(input.as_bytes(), config)
}

/// Parse SMILES bytes with strict OpenSMILES rules
pub fn parse_smiles_bytes(input: &[u8]) -> Result<Molecule, ParseError> {
    parse_smiles_bytes_with(input, &SmilesIoConfig::basic_opensmiles())
}

/// Parse SMILES bytes with configuration
pub fn parse_smiles_bytes_with(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<Molecule, ParseError> {
    let flags = config.parse_flags;

    // Strip leading/trailing whitespace, reject internal whitespace
    let mut start = 0usize;
    while start < input.len() && matches!(input[start], b' ' | b'\t' | b'\n' | b'\r') {
        start += 1;
    }
    if start == input.len() {
        return Ok(Molecule::empty());
    }
    if start > 0 {
        return Err(ParseError::InvalidWhitespace { pos: 0 });
    }
    let mut end = input.len();
    while end > 0 && matches!(input[end - 1], b' ' | b'\t' | b'\n' | b'\r') {
        end -= 1;
    }
    for (k, b) in input[start..end].iter().enumerate() {
        if matches!(*b, b' ' | b'\t' | b'\n' | b'\r') {
            return Err(ParseError::InvalidWhitespace { pos: start + k });
        }
    }

    parse_smiles_inner(&input[start..end], flags)
}

fn parse_smiles_inner(input: &[u8], _flags: SmilesParseFlags) -> Result<Molecule, ParseError> {
    let mut i = 0usize;
    let n = input.len();
    let mut builder = MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut branch_stack: Vec<Frame> = Vec::new();
    let mut ring_table: Vec<Option<OpenRing>> = Vec::new();
    let mut last_atom_idx: Option<u32> = None;
    let mut pending_bond: Option<(BondOrder, Option<BondDirection>, usize)> = None;
    let mut last_aromatic: bool = false;
    let mut just_closed_group: bool = false;

    while i < n {
        let b0 = input[i];

        // Reject whitespace (already stripped, but catch internal)
        if matches!(b0, b' ' | b'\t' | b'\n' | b'\r') {
            return Err(ParseError::InvalidWhitespace { pos: i });
        }

        if b0 != b'(' {
            just_closed_group = false;
        }
        if b0 == b'(' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos });
            }
            if just_closed_group {
                last_atom_idx = None;
                branch_stack.push(Frame::Group {
                    had_atom: false,
                    open_pos: i,
                });
                just_closed_group = false;
            } else {
                match last_atom_idx {
                    Some(idx) => branch_stack.push(Frame::Branch {
                        base: idx,
                        had_atom: false,
                        open_pos: i,
                    }),
                    None => branch_stack.push(Frame::Group {
                        had_atom: false,
                        open_pos: i,
                    }),
                }
            }
            i += 1;
            continue;
        }
        if b0 == b')' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos });
            }
            let Some(frame) = branch_stack.pop() else {
                return Err(ParseError::UnbalancedCloseParen { pos: i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(ParseError::EmptyBranch { pos: i });
                    }
                    last_atom_idx = Some(base);
                }
                Frame::Group {
                    had_atom, open_pos, ..
                } => {
                    if !had_atom {
                        if i + 1 != n {
                            return Err(ParseError::EmptyGroup { pos: i });
                        }
                        if i > 0 && input[i - 1] == b'.' {
                            return Err(ParseError::LeadingDot { pos: i - 1 });
                        }
                        if open_pos != 0 {
                            return Err(ParseError::EmptyGroup { pos: i });
                        }
                        last_atom_idx = None;
                        just_closed_group = false;
                    } else {
                        just_closed_group = true;
                        if branch_stack.is_empty() && i + 1 != n {
                            let next = input[i + 1];
                            if next != b'.' {
                                return Err(ParseError::NonfinalGroup { pos: i });
                            }
                        }
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'.' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos });
            }
            if i == 0 {
                return Err(ParseError::LeadingDot { pos: i });
            }
            if let Some(Frame::Group {
                had_atom: false, ..
            }) = branch_stack.last()
            {
                return Err(ParseError::LeadingDot { pos: i });
            }
            if i + 1 == n {
                return Err(ParseError::TrailingDot { pos: i });
            }
            if input[i + 1] == b'.' {
                return Err(ParseError::ConsecutiveDots { pos: i });
            }
            // Detect dot before ring (single digit ring index)
            if input[i + 1].is_ascii_digit() {
                return Err(ParseError::DotBeforeRing { pos: i });
            }
            // Detect dot before percent ring index
            if input[i + 1] == b'%' {
                return Err(ParseError::DotBeforeRing { pos: i });
            }
            last_atom_idx = None;
            last_aromatic = false;
            i += 1;
            continue;
        }
        match parse_ring_index(input, i) {
            Ok(Some((idx, next_i, _percent))) => {
                if last_atom_idx.is_none() {
                    return Err(ParseError::LeadingRing { pos: i });
                }
                if invalid_ring_context(&branch_stack) {
                    return Err(ParseError::LeadingRing { pos: 0 });
                }
                let bond = pending_bond.take();
                let (order_opt, dir_opt) = bond.map_or((None, None), |(o, d, _)| (Some(o), d));
                process_ring_closure(
                    &mut ring_table,
                    &mut builder,
                    last_aromatic,
                    last_atom_idx.unwrap(),
                    idx,
                    order_opt,
                    dir_opt,
                    i,
                    i + 1,
                )?;
                i = next_i;
                continue;
            }
            Err(e) => return Err(e),
            Ok(None) => {}
        }
        // percent branch is handled by parse_ring_index above
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\') {
            if pending_bond.is_some() {
                return Err(ParseError::ConsecutiveBonds { pos: i });
            }
            if last_atom_idx.is_none() {
                if let Some(Frame::Group {
                    had_atom: false, ..
                }) = branch_stack.last()
                {
                    return Err(ParseError::LeadingBond { pos: i });
                }
                return Err(ParseError::LeadingBond { pos: i });
            }
            let (order, bond_dir) = parse_bond(b0);
            pending_bond = Some((order, bond_dir, i));
            i += 1;
            continue;
        }
        if b0 == b'[' {
            let start = i + 1;
            let mut j = start;
            while j < n && input[j] != b']' {
                j += 1;
            }
            if j >= n {
                return Err(ParseError::UnbalancedOpenBracket { pos: i });
            }
            // Empty bracket []
            if j == start {
                return Err(ParseError::EmptyBracket { pos: i });
            }
            let inner = &input[start..j];
            let (elem_opt, aromatic, iso_opt, charge_opt, class_opt, h_opt, chir_opt) =
                parse_bracket(inner, i)?;
            let (element, aromatic) = match elem_opt {
                Some(e) => (e, aromatic),
                None => (Element::C, false),
            };
            let (s, e) = (Some(i as u32), Some((j + 1) as u32));
            let atom = AtomData {
                element,
                isotope: iso_opt,
                charge: charge_opt,
                hydrogen_count: h_opt,
                class: class_opt,
                aromatic,
                implicit_h: false,
                chirality: chir_opt,
                span: Span::from_bytes_opt(s, e),
            };
            let curr = builder.on_atom(atom);

            attach_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                aromatic,
                i as u32,
                (j + 1) as u32,
            );
            last_atom_idx = Some(curr);
            last_aromatic = aromatic;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i = j + 1;
            continue;
        }
        if b0 == b'C' {
            if i + 1 < n && input[i + 1] == b'l' {
                let (s, e) = (Some(i as u32), Some((i + 2) as u32));
                let curr = builder.on_atom_fast(Element::Cl, false, s, e);

                attach_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    false,
                    i as u32,
                    (i + 2) as u32,
                );
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += 2;
                continue;
            }
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            let curr = builder.on_atom_fast(Element::C, false, s, e);

            attach_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'B' {
            if i + 1 < n && input[i + 1] == b'r' {
                let (s, e) = (Some(i as u32), Some((i + 2) as u32));
                let curr = builder.on_atom_fast(Element::Br, false, s, e);

                attach_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    false,
                    i as u32,
                    (i + 2) as u32,
                );
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += 2;
                continue;
            }
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            let curr = builder.on_atom_fast(Element::B, false, s, e);

            attach_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }
        // Elements
        if b0.is_ascii_alphabetic() {
            if let Some((element, consumed)) = parse_organic_aliphatic_element(input, i) {
                let (s, e) = (Some(i as u32), Some((i + consumed) as u32));
                let curr = builder.on_atom_fast(element, false, s, e);

                attach_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    false,
                    i as u32,
                    (i + consumed) as u32,
                );
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += consumed;
                continue;
            }
            if let Some((element, consumed)) = parse_organic_aromatic_element(input, i) {
                let (s, e) = (Some(i as u32), Some((i + consumed) as u32));
                let curr = builder.on_atom_fast(element, true, s, e);

                attach_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    true,
                    i as u32,
                    (i + consumed) as u32,
                );
                last_atom_idx = Some(curr);
                last_aromatic = true;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += consumed;
                continue;
            }
            return Err(ParseError::InvalidElement { pos: i });
        }
        if b0 == b'*' {
            // Wildcards not supported in basic SMILES parser
            return Err(ParseError::InvalidElement { pos: i });
        }
        if b0 == b']' {
            return Err(ParseError::UnbalancedCloseBracket { pos: i });
        }
        // Bracket-only fields outside bracket
        if b0 == b'@' || b0 == b'+' {
            return Err(ParseError::StrayBracketField { pos: i });
        }
        return Err(ParseError::InvalidToken { pos: i });
    }

    if let Some((_, _, pos)) = pending_bond {
        return Err(ParseError::TrailingBond { pos });
    }
    if !branch_stack.is_empty() {
        let pos = match branch_stack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(ParseError::UnbalancedOpenParen { pos });
    }
    let mut last_open: Option<usize> = None;
    for entry in ring_table.iter().flatten() {
        match last_open {
            None => last_open = Some(entry.open_pos),
            Some(p) => {
                if entry.open_pos > p {
                    last_open = Some(entry.open_pos)
                }
            }
        }
    }
    if let Some(pos_open) = last_open {
        return Err(ParseError::UnbalancedRingIndex { open_pos: pos_open });
    }
    let mut mols = builder.finish();
    let mol = mols.pop().unwrap_or_else(Molecule::empty);
    Ok(mol)
}

use crate::table_ir::ExtendedMolecule;

/// Parse extended SMILES string with strict OpenSMILES rules
pub fn parse_extended_smiles(input: &str) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_smiles_bytes(input.as_bytes())
}

/// Parse extended SMILES string with configuration
pub fn parse_extended_smiles_with(
    input: &str,
    config: &SmilesIoConfig,
) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_smiles_bytes_with(input.as_bytes(), config)
}

/// Parse extended SMILES bytes with strict OpenSMILES rules
pub fn parse_extended_smiles_bytes(input: &[u8]) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_smiles_bytes_with(input, &SmilesIoConfig::basic_opensmiles())
}

/// Parse extended SMILES bytes with configuration
pub fn parse_extended_smiles_bytes_with(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<ExtendedMolecule, ParseError> {
    let flags = config.parse_flags;

    // Strip leading/trailing whitespace, reject internal whitespace
    let mut start = 0usize;
    while start < input.len() && matches!(input[start], b' ' | b'\t' | b'\n' | b'\r') {
        start += 1;
    }
    if start == input.len() {
        return Ok(ExtendedMolecule::empty());
    }
    if start > 0 {
        return Err(ParseError::InvalidWhitespace { pos: 0 });
    }
    let mut end = input.len();
    while end > 0 && matches!(input[end - 1], b' ' | b'\t' | b'\n' | b'\r') {
        end -= 1;
    }
    for (k, b) in input[start..end].iter().enumerate() {
        if matches!(*b, b' ' | b'\t' | b'\n' | b'\r') {
            return Err(ParseError::InvalidWhitespace { pos: start + k });
        }
    }

    parse_extended_smiles_inner(&input[start..end], flags)
}

fn parse_extended_smiles_inner(
    input: &[u8],
    _flags: SmilesParseFlags,
) -> Result<ExtendedMolecule, ParseError> {
    let mut i = 0usize;
    let n = input.len();
    let mut builder =
        ExtendedMoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1));
    let mut branch_stack: Vec<Frame> = Vec::new();
    let mut ring_table: Vec<Option<OpenRing>> = Vec::new();
    let mut last_atom_idx: Option<u32> = None;
    let mut pending_bond: Option<(BondOrder, Option<BondDirection>, usize)> = None;
    let mut last_aromatic: bool = false;
    let mut just_closed_group: bool = false;

    while i < n {
        let b0 = input[i];

        if matches!(b0, b' ' | b'\t' | b'\n' | b'\r') {
            return Err(ParseError::InvalidWhitespace { pos: i });
        }

        if b0 != b'(' {
            just_closed_group = false;
        }
        if b0 == b'(' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos });
            }
            if just_closed_group {
                last_atom_idx = None;
                branch_stack.push(Frame::Group {
                    had_atom: false,
                    open_pos: i,
                });
                just_closed_group = false;
            } else {
                match last_atom_idx {
                    Some(idx) => branch_stack.push(Frame::Branch {
                        base: idx,
                        had_atom: false,
                        open_pos: i,
                    }),
                    None => branch_stack.push(Frame::Group {
                        had_atom: false,
                        open_pos: i,
                    }),
                }
            }
            i += 1;
            continue;
        }
        if b0 == b')' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos });
            }
            let Some(frame) = branch_stack.pop() else {
                return Err(ParseError::UnbalancedCloseParen { pos: i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(ParseError::EmptyBranch { pos: i });
                    }
                    last_atom_idx = Some(base);
                }
                Frame::Group {
                    had_atom, open_pos, ..
                } => {
                    if !had_atom {
                        if i + 1 != n {
                            return Err(ParseError::EmptyGroup { pos: i });
                        }
                        if i > 0 && input[i - 1] == b'.' {
                            return Err(ParseError::LeadingDot { pos: i - 1 });
                        }
                        if open_pos != 0 {
                            return Err(ParseError::EmptyGroup { pos: i });
                        }
                        last_atom_idx = None;
                        just_closed_group = false;
                    } else {
                        just_closed_group = true;
                        if branch_stack.is_empty() && i + 1 != n {
                            let next = input[i + 1];
                            if next != b'.' {
                                return Err(ParseError::NonfinalGroup { pos: i });
                            }
                        }
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'.' {
            if let Some((_, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos });
            }
            if i == 0 {
                return Err(ParseError::LeadingDot { pos: i });
            }
            if let Some(Frame::Group {
                had_atom: false, ..
            }) = branch_stack.last()
            {
                return Err(ParseError::LeadingDot { pos: i });
            }
            if i + 1 == n {
                return Err(ParseError::TrailingDot { pos: i });
            }
            if input[i + 1] == b'.' {
                return Err(ParseError::ConsecutiveDots { pos: i });
            }
            if input[i + 1].is_ascii_digit() {
                return Err(ParseError::DotBeforeRing { pos: i });
            }
            if input[i + 1] == b'%' {
                return Err(ParseError::DotBeforeRing { pos: i });
            }
            last_atom_idx = None;
            last_aromatic = false;
            i += 1;
            continue;
        }
        match parse_ring_index(input, i) {
            Ok(Some((idx, next_i, _percent))) => {
                if last_atom_idx.is_none() {
                    return Err(ParseError::LeadingRing { pos: i });
                }
                if invalid_ring_context(&branch_stack) {
                    return Err(ParseError::LeadingRing { pos: 0 });
                }
                let bond = pending_bond.take();
                let (order_opt, dir_opt) = bond.map_or((None, None), |(o, d, _)| (Some(o), d));
                process_ring_closure_extended(
                    &mut ring_table,
                    &mut builder,
                    last_aromatic,
                    last_atom_idx.unwrap(),
                    idx,
                    order_opt,
                    dir_opt,
                    i,
                    i + 1,
                )?;
                i = next_i;
                continue;
            }
            Err(e) => return Err(e),
            Ok(None) => {}
        }
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\') {
            if pending_bond.is_some() {
                return Err(ParseError::ConsecutiveBonds { pos: i });
            }
            if last_atom_idx.is_none() {
                if let Some(Frame::Group {
                    had_atom: false, ..
                }) = branch_stack.last()
                {
                    return Err(ParseError::LeadingBond { pos: i });
                }
                return Err(ParseError::LeadingBond { pos: i });
            }
            let (order, bond_dir) = parse_bond(b0);
            pending_bond = Some((order, bond_dir, i));
            i += 1;
            continue;
        }
        if b0 == b'[' {
            let start = i + 1;
            let mut j = start;
            while j < n && input[j] != b']' {
                j += 1;
            }
            if j >= n {
                return Err(ParseError::UnbalancedOpenBracket { pos: i });
            }
            if j == start {
                return Err(ParseError::EmptyBracket { pos: i });
            }
            let inner = &input[start..j];
            let (symbol, aromatic, iso_opt, charge_opt, class_opt, h_opt, chir_opt) =
                parse_bracket_extended(inner, i)?;
            let (s, e) = (Some(i as u32), Some((j + 1) as u32));
            let atom = ExtendedAtomData {
                symbol,
                isotope: iso_opt,
                charge: charge_opt,
                hydrogen_count: h_opt,
                class: class_opt,
                aromatic,
                implicit_h: false,
                chirality: chir_opt,
                span: Span::from_bytes_opt(s, e),
            };
            let curr = builder.on_atom(atom);

            attach_atom_extended(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                aromatic,
                i as u32,
                (j + 1) as u32,
            );
            last_atom_idx = Some(curr);
            last_aromatic = aromatic;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i = j + 1;
            continue;
        }
        if b0 == b'C' {
            if i + 1 < n && input[i + 1] == b'l' {
                let (s, e) = (Some(i as u32), Some((i + 2) as u32));
                let curr = builder.on_atom_fast(Element::Cl, false, s, e);

                attach_atom_extended(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    false,
                    i as u32,
                    (i + 2) as u32,
                );
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += 2;
                continue;
            }
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            let curr = builder.on_atom_fast(Element::C, false, s, e);

            attach_atom_extended(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'B' {
            if i + 1 < n && input[i + 1] == b'r' {
                let (s, e) = (Some(i as u32), Some((i + 2) as u32));
                let curr = builder.on_atom_fast(Element::Br, false, s, e);

                attach_atom_extended(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    false,
                    i as u32,
                    (i + 2) as u32,
                );
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += 2;
                continue;
            }
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            let curr = builder.on_atom_fast(Element::B, false, s, e);

            attach_atom_extended(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0.is_ascii_alphabetic() {
            if let Some((element, consumed)) = parse_organic_aliphatic_element(input, i) {
                let (s, e) = (Some(i as u32), Some((i + consumed) as u32));
                let curr = builder.on_atom_fast(element, false, s, e);

                attach_atom_extended(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    false,
                    i as u32,
                    (i + consumed) as u32,
                );
                last_atom_idx = Some(curr);
                last_aromatic = false;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += consumed;
                continue;
            }
            if let Some((element, consumed)) = parse_organic_aromatic_element(input, i) {
                let (s, e) = (Some(i as u32), Some((i + consumed) as u32));
                let curr = builder.on_atom_fast(element, true, s, e);

                attach_atom_extended(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    last_aromatic,
                    true,
                    i as u32,
                    (i + consumed) as u32,
                );
                last_atom_idx = Some(curr);
                last_aromatic = true;
                if let Some(top) = branch_stack.last_mut() {
                    match top {
                        Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                            *had_atom = true
                        }
                    }
                }
                i += consumed;
                continue;
            }
            return Err(ParseError::InvalidElement { pos: i });
        }
        if b0 == b'*' {
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            let curr = builder.on_wildcard(WildcardAtom::Any, None, s, e);

            attach_atom_extended(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                last_aromatic,
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
            last_aromatic = false;
            if let Some(top) = branch_stack.last_mut() {
                match top {
                    Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                        *had_atom = true
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b']' {
            return Err(ParseError::UnbalancedCloseBracket { pos: i });
        }
        if b0 == b'@' || b0 == b'+' {
            return Err(ParseError::StrayBracketField { pos: i });
        }
        return Err(ParseError::InvalidToken { pos: i });
    }

    if let Some((_, _, pos)) = pending_bond {
        return Err(ParseError::TrailingBond { pos });
    }
    if !branch_stack.is_empty() {
        let pos = match branch_stack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(ParseError::UnbalancedOpenParen { pos });
    }
    let mut last_open: Option<usize> = None;
    for entry in ring_table.iter().flatten() {
        match last_open {
            None => last_open = Some(entry.open_pos),
            Some(p) => {
                if entry.open_pos > p {
                    last_open = Some(entry.open_pos)
                }
            }
        }
    }
    if let Some(pos_open) = last_open {
        return Err(ParseError::UnbalancedRingIndex { open_pos: pos_open });
    }
    let mut mols = builder.finish();
    let mol = mols.pop().unwrap_or_else(ExtendedMolecule::empty);
    Ok(mol)
}

#[cfg(test)]
mod tests;
