//! SMILES parser

use std::collections::BTreeMap;

use indexmap::IndexMap;
use umol_shared::element::Element;

mod builder;
mod cx;
mod utils;

use umol_ast::ast::{MoleculeAst, TryIntoAst};
use umol_shared::error::UmolError;

use self::builder::{AtomData, ExtendedAtomData, ExtendedMoleculeBuilder, MoleculeBuilder};
use self::cx::{
    parse_cx_annotations, parse_extended_cx_annotations, remap_cx_bond_indices,
    split_reaction_cx_entries, update_extended_molecule, update_extended_reaction, update_molecule,
    update_reaction, BondIndexMap,
};
use self::utils::{
    attach_atom, attach_extended_atom, invalid_ring_context, parse_bond, parse_bracket,
    parse_extended_bond, parse_extended_bracket, parse_organic_aliphatic_element,
    parse_organic_aromatic_element, parse_ring_index, Frame,
};
use super::config::{SmilesIoConfig, SmilesParseFlags};
use super::error::ParseError;
use crate::table_ir::{
    BondDonation, BondOrder, BondDirection, ChiralityFrame, ExtendedMolecule, ExtendedReaction,
    Molecule, Reaction, SourceFormat, Span, WildcardAtom,
};

/// Parse SMILES to [`MoleculeAst`] without running the solver.
pub fn parse_smiles_to_ast(input: &str) -> Result<MoleculeAst, Box<dyn UmolError>> {
    parse_smiles_bytes_to_ast(input.as_bytes())
}

/// Parse SMILES bytes to [`MoleculeAst`] without running the solver. Spans the parse (`ParseError`)
/// and raise (`RaiseError`) concerns, so per doc 065 it returns the boxed boundary error.
pub fn parse_smiles_bytes_to_ast(input: &[u8]) -> Result<MoleculeAst, Box<dyn UmolError>> {
    let table_mol = parse_smiles_bytes_to_table_ir(input)?;
    let ast: MoleculeAst = (&table_mol).try_into_ast(&())?;
    Ok(ast)
}

/// Parse SMILES string to `table_ir::Molecule` with basic OpenSMILES configuration.
pub fn parse_smiles_to_table_ir(input: &str) -> Result<Molecule, ParseError> {
    parse_smiles_bytes_to_table_ir(input.as_bytes())
}

/// Parse SMILES string to `table_ir::Molecule` with configuration.
pub fn parse_smiles_to_table_ir_with(
    input: &str,
    config: &SmilesIoConfig,
) -> Result<Molecule, ParseError> {
    parse_smiles_bytes_to_table_ir_with(input.as_bytes(), config)
}

/// Parse SMILES bytes to `table_ir::Molecule` with basic OpenSMILES rules.
pub fn parse_smiles_bytes_to_table_ir(input: &[u8]) -> Result<Molecule, ParseError> {
    parse_smiles_bytes_to_table_ir_with(input, &SmilesIoConfig::basic_opensmiles())
}

/// Parse SMILES bytes to `table_ir::Molecule` with configuration.
pub fn parse_smiles_bytes_to_table_ir_with(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<Molecule, ParseError> {
    let flags = config.parse_flags;
    debug_assert!(
        SmilesParseFlags::BASIC_MAX.contains(flags),
        "flags must be a subset of BASIC_MAX, got: {}",
        flags
    );

    if input.is_empty() {
        return Ok(Molecule::empty());
    }

    // Leading whitespace is not allowed (exception: whitespace-only input is allowed)
    if !input.is_empty() && input[0].is_ascii_whitespace() && !input.trim_ascii_start().is_empty() {
        return Err(ParseError::LeadingWhitespace);
    }

    // Check if the input contains a CX block, record ring bonds if it is present.
    let has_cx_annotations =
        flags.contains(SmilesParseFlags::CHEMAXON_EXTENSIONS) && input.contains(&b'|');
    let (remaining, (mut mol, ring_bonds, _)) =
        parse_smiles_inner(input, 0, false, has_cx_annotations, flags)?;

    // Inner parser stops at whitespace.
    let trimmed = remaining.trim_ascii_start();
    if trimmed.is_empty() {
        return Ok(mol);
    }

    // Chemaxon extensions
    if has_cx_annotations {
        let mut entries = parse_cx_annotations(trimmed, flags)?;
        let bond_map = BondIndexMap::new(ring_bonds, mol.bonds.len());
        remap_cx_bond_indices(&mut entries, &bond_map)?;
        update_molecule(&mut mol, entries)?;
    }

    Ok(mol)
}

/// Parse reaction SMILES with basic OpenSMILES configuration
pub fn parse_reaction_smiles(input: &str) -> Result<Reaction, ParseError> {
    parse_reaction_smiles_bytes(input.as_bytes())
}

/// Parse reaction SMILES string with configuration
pub fn parse_reaction_smiles_with(
    input: &str,
    config: &SmilesIoConfig,
) -> Result<Reaction, ParseError> {
    parse_reaction_smiles_bytes_with(input.as_bytes(), config)
}

/// Parse reaction SMILES bytes with basic OpenSMILES configuration
pub fn parse_reaction_smiles_bytes(input: &[u8]) -> Result<Reaction, ParseError> {
    parse_reaction_smiles_bytes_with(input, &SmilesIoConfig::basic_opensmiles())
}

/// Parse reaction SMILES bytes with configuration
pub fn parse_reaction_smiles_bytes_with(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<Reaction, ParseError> {
    let flags = config.parse_flags;
    debug_assert!(
        SmilesParseFlags::BASIC_MAX.contains(flags),
        "flags must be a subset of BASIC_MAX"
    );

    // Check if the input contains a CX block, record ring bonds if it is present.
    let has_cx_annotations =
        flags.contains(SmilesParseFlags::CHEMAXON_EXTENSIONS) && input.contains(&b'|');
    let mut remaining = input;
    let mut offset = 0usize;

    if remaining.starts_with(b".") {
        return Err(ParseError::LeadingDot { pos: 0 });
    }

    // Leading whitespace is not allowed
    if !input.is_empty() && input[0].is_ascii_whitespace() {
        return Err(ParseError::LeadingWhitespace);
    }

    // Reactants: parse one side-supermolecule until '>'.
    let (rest, (reactants, reactant_ring_bonds, new_offset)) =
        parse_smiles_inner(remaining, offset, true, has_cx_annotations, flags)?;
    offset = new_offset;
    remaining = rest;

    let agents;
    let agent_ring_bonds;
    if remaining.starts_with(b">>") {
        remaining = &remaining[2..];
        offset += 2;
        agents = Molecule::empty();
        agent_ring_bonds = Vec::new();
    } else if remaining.starts_with(b">") {
        remaining = &remaining[1..];
        offset += 1;

        // Agents: parse one side-supermolecule until '>'.
        let (rest, (agents_parsed, agents_ring_bonds, new_offset)) =
            parse_smiles_inner(remaining, offset, true, has_cx_annotations, flags)?;
        offset = new_offset;
        remaining = rest;
        agents = agents_parsed;
        agent_ring_bonds = agents_ring_bonds;

        if !remaining.starts_with(b">") {
            return Err(ParseError::MissingReactionArrow { pos: offset });
        }
        remaining = &remaining[1..];
        offset += 1;
    } else {
        return Err(ParseError::MissingReactionArrow { pos: offset });
    }

    // Products: parse one side-supermolecule until EOF/whitespace.
    let (rest, (products, product_ring_bonds, _new_offset)) =
        parse_smiles_inner(remaining, offset, true, has_cx_annotations, flags)?;

    let mut reaction = Reaction {
        reactants,
        products,
        agents,
        atom_mapping: BTreeMap::new(),
        comments: Vec::new(),
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    if has_cx_annotations {
        let entries = parse_cx_annotations(rest.trim_ascii_start(), flags)?;
        let mut split = split_reaction_cx_entries(
            entries,
            reaction.reactants.atom_count(),
            reaction.reactants.bond_count(),
            reaction.agents.atom_count(),
            reaction.agents.bond_count(),
            reaction.products.atom_count(),
            reaction.products.bond_count(),
        )?;
        remap_cx_bond_indices(
            &mut split.0,
            &BondIndexMap::new(reactant_ring_bonds, reaction.reactants.bond_count()),
        )?;
        remap_cx_bond_indices(
            &mut split.1,
            &BondIndexMap::new(agent_ring_bonds, reaction.agents.bond_count()),
        )?;
        remap_cx_bond_indices(
            &mut split.2,
            &BondIndexMap::new(product_ring_bonds, reaction.products.bond_count()),
        )?;
        update_reaction(&mut reaction, split)?;
    }
    collect_atom_mapping(&mut reaction);
    Ok(reaction)
}

// The third tuple slot carries the ring-closure records for CX bond-index remapping
// (a tuple for now; a named parse-result struct is the eventual home).
#[allow(clippy::type_complexity)]
fn parse_smiles_inner(
    input: &[u8],
    offset: usize,
    as_reaction: bool,
    store_rings: bool,
    flags: SmilesParseFlags,
) -> Result<(&[u8], (Molecule, Vec<(usize, usize)>, usize)), ParseError> {
    let extended_bonds = flags.contains(SmilesParseFlags::EXTENDED_BONDS);
    let mut i = 0usize;
    let n = input.len();
    let mut builder =
        MoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1), store_rings);
    let mut branch_stack: Vec<Frame> = Vec::new();
    let mut last_atom_idx: Option<usize> = None;
    let mut pending_bond: Option<(BondOrder, Option<BondDirection>, Option<BondDonation>, usize)> =
        None;
    let mut after_closed_group: bool = false;

    while i < n {
        let b0 = input[i];

        // Stop at whitespace - return remaining input
        if b0.is_ascii_whitespace() {
            break;
        }
        if as_reaction && b0 == b'>' {
            break;
        }

        if b0 != b'(' {
            after_closed_group = false;
        }
        if b0 == b'(' {
            if let Some((_, _, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos: offset + pos });
            }
            if after_closed_group {
                last_atom_idx = None;
                branch_stack.push(Frame::Group {
                    had_atom: false,
                    open_pos: i,
                });
                after_closed_group = false;
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
            if let Some((_, _, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos: offset + pos });
            }
            let Some(frame) = branch_stack.pop() else {
                return Err(ParseError::UnbalancedCloseParen { pos: offset + i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(ParseError::EmptyBranch { pos: offset + i });
                    }
                    last_atom_idx = Some(base);
                }
                Frame::Group { had_atom, .. } => {
                    if !had_atom {
                        return Err(ParseError::EmptyGroup { pos: offset + i });
                    }
                    after_closed_group = true;
                    if let Some(parent) = branch_stack.last_mut() {
                        match parent {
                            Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                                *had_atom = true
                            }
                        }
                    }
                    if branch_stack.is_empty() && i + 1 != n {
                        let next = input[i + 1];
                        if next != b'.' {
                            return Err(ParseError::NonfinalGroup { pos: offset + i });
                        }
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'.' {
            if let Some((_, _, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos: offset + pos });
            }
            if i == 0 {
                return Err(ParseError::LeadingDot { pos: offset + i });
            }
            if let Some(Frame::Group {
                had_atom: false, ..
            }) = branch_stack.last()
            {
                return Err(ParseError::LeadingDot { pos: offset + i });
            }
            if i + 1 == n {
                return Err(ParseError::TrailingDot { pos: offset + i });
            }
            if as_reaction && input[i + 1] == b'>' {
                return Err(ParseError::TrailingDot { pos: offset + i });
            }
            if input[i + 1] == b'.' {
                return Err(ParseError::ConsecutiveDots { pos: offset + i });
            }
            // Detect dot before ring (single digit ring index)
            if input[i + 1].is_ascii_digit() {
                return Err(ParseError::DotBeforeRing { pos: offset + i });
            }
            // Detect dot before percent ring index
            if input[i + 1] == b'%' {
                return Err(ParseError::DotBeforeRing { pos: offset + i });
            }
            last_atom_idx = None;
            i += 1;
            continue;
        }
        match parse_ring_index(input, i, offset) {
            Ok(Some((idx, next_i, _percent))) => {
                if last_atom_idx.is_none() {
                    return Err(ParseError::LeadingRing { pos: offset + i });
                }
                if invalid_ring_context(&branch_stack) {
                    return Err(ParseError::LeadingRing { pos: offset });
                }
                let bond = pending_bond.take();
                let (order_opt, direction_opt, donation_opt) =
                    bond.map_or((None, None, None), |(o, d, don, _)| (Some(o), d, don));
                builder.on_ring_bond(
                    last_atom_idx.unwrap(),
                    idx,
                    order_opt,
                    direction_opt,
                    donation_opt,
                    i,
                    i + 1,
                    offset,
                )?;
                i = next_i;
                continue;
            }
            Err(e) => return Err(e),
            Ok(None) => {}
        }
        // percent branch is handled by parse_ring_index above
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\')
            || (extended_bonds && matches!(b0, b'~' | b'<'))
        {
            if pending_bond.is_some() {
                return Err(ParseError::ConsecutiveBonds { pos: offset + i });
            }
            if last_atom_idx.is_none() {
                if let Some(Frame::Group {
                    had_atom: false, ..
                }) = branch_stack.last()
                {
                    return Err(ParseError::LeadingBond { pos: offset + i });
                }
                return Err(ParseError::LeadingBond { pos: offset + i });
            }
            // Use extended bond parsing for ->, <-, ~ when EXTENDED_BONDS is set
            if extended_bonds {
                let (order, direction, donation, consumed) = parse_extended_bond(input, i);
                pending_bond = Some((order, direction, donation, i));
                i += consumed;
            } else {
                let (order, bond_direction) = parse_bond(b0);
                pending_bond = Some((order, bond_direction, None, i));
                i += 1;
            }
            continue;
        }
        if b0 == b'[' {
            let start = i + 1;
            let mut j = start;
            while j < n && input[j] != b']' {
                j += 1;
            }
            if j >= n {
                return Err(ParseError::UnbalancedOpenBracket { pos: offset + i });
            }
            // Empty bracket []
            if j == start {
                return Err(ParseError::EmptyBracket { pos: offset + i });
            }
            let inner = &input[start..j];
            let (element, aromatic, iso_opt, charge_opt, class_opt, h_opt, chir_opt) =
                parse_bracket(inner, offset + i, flags)?;
            let (s, e) = (Some(i as u32), Some((j + 1) as u32));
            let atom = AtomData {
                element,
                isotope: iso_opt,
                charge: charge_opt,
                implicit_hydrogens: Some(h_opt.unwrap_or(0)),
                class: class_opt,
                aromatic,
                chirality: chir_opt,
                span: Span::from_bytes_opt(s, e),
            };
            let curr = builder.on_atom(atom);

            attach_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                aromatic,
                i as u32,
                (j + 1) as u32,
            );
            last_atom_idx = Some(curr);
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
                    false,
                    i as u32,
                    (i + 2) as u32,
                );
                last_atom_idx = Some(curr);
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
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
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
                    false,
                    i as u32,
                    (i + 2) as u32,
                );
                last_atom_idx = Some(curr);
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
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
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
                    false,
                    i as u32,
                    (i + consumed) as u32,
                );
                last_atom_idx = Some(curr);
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
                    true,
                    i as u32,
                    (i + consumed) as u32,
                );
                last_atom_idx = Some(curr);
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
            return Err(ParseError::InvalidElement { pos: offset + i });
        }
        if b0 == b'*' {
            // Wildcards not supported in basic SMILES parser
            return Err(ParseError::InvalidElement { pos: offset + i });
        }
        if b0 == b']' {
            return Err(ParseError::UnbalancedCloseBracket { pos: offset + i });
        }
        // Bracket-only fields outside bracket
        if b0 == b'@' || b0 == b'+' {
            return Err(ParseError::StrayBracketField { pos: offset + i });
        }
        return Err(ParseError::InvalidToken { pos: offset + i });
    }

    if let Some((_, _, _, pos)) = pending_bond {
        return Err(ParseError::TrailingBond { pos: offset + pos });
    }
    if !branch_stack.is_empty() {
        let pos = match branch_stack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(ParseError::UnbalancedOpenParen { pos: offset + pos });
    }
    if let Some(pos_open) = builder.unclosed_ring_pos() {
        return Err(ParseError::UnbalancedRingIndex {
            open_pos: offset + pos_open,
        });
    }
    let mut mols = builder.finish();
    let mol = mols
        .pop()
        .map(|mut mol| {
            mol.source_format = SourceFormat::SMILES;
            mol.chirality_frame = Some(ChiralityFrame::FirstNeighborToward);
            mol
        })
        .unwrap_or_else(Molecule::empty);
    let ring_bonds = builder.take_ring_bonds();
    let new_offset = offset + i;
    Ok((&input[i..], (mol, ring_bonds, new_offset)))
}

/// Parse extended SMILES string with basic OpenSMILES configuration
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

/// Parse extended SMILES bytes with basic OpenSMILES configuration
pub fn parse_extended_smiles_bytes(input: &[u8]) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_smiles_bytes_with(input, &SmilesIoConfig::basic_opensmiles())
}

/// Parse extended SMILES bytes with configuration
pub fn parse_extended_smiles_bytes_with(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<ExtendedMolecule, ParseError> {
    let flags = config.parse_flags;

    if input.is_empty() {
        return Ok(ExtendedMolecule::empty());
    }

    // Check if the input contains a CX block, record ring bonds if it is present.
    let has_cx_annotations =
        flags.contains(SmilesParseFlags::CHEMAXON_EXTENSIONS) && input.contains(&b'|');
    let (remaining, (mut mol, ring_bonds, _)) =
        parse_extended_smiles_inner(input, 0, false, has_cx_annotations, flags)?;

    // Inner parser stops at whitespace. Leading whitespace is not allowed
    // (exception: whitespace-only input is allowed)
    if mol.atoms.is_empty() && !remaining.trim_ascii_start().is_empty() {
        return Err(ParseError::LeadingWhitespace);
    }

    if remaining.is_empty() {
        return Ok(mol);
    }

    // Chemaxon annotations
    if has_cx_annotations {
        let mut entries = parse_extended_cx_annotations(remaining.trim_ascii_start(), flags)?;
        let bond_map = BondIndexMap::new(ring_bonds, mol.bonds.len());
        remap_cx_bond_indices(&mut entries, &bond_map)?;
        update_extended_molecule(&mut mol, entries)?;
    }

    Ok(mol)
}

/// Parse extended reaction SMILES string with basic OpenSMILES configuration
pub fn parse_extended_reaction_smiles(input: &str) -> Result<ExtendedReaction, ParseError> {
    parse_extended_reaction_smiles_bytes(input.as_bytes())
}

/// Parse extended reaction SMILES string with configuration
pub fn parse_extended_reaction_smiles_with(
    input: &str,
    config: &SmilesIoConfig,
) -> Result<ExtendedReaction, ParseError> {
    parse_extended_reaction_smiles_bytes_with(input.as_bytes(), config)
}

/// Parse extended reaction SMILES bytes with basic OpenSMILES configuration
pub fn parse_extended_reaction_smiles_bytes(input: &[u8]) -> Result<ExtendedReaction, ParseError> {
    parse_extended_reaction_smiles_bytes_with(input, &SmilesIoConfig::basic_opensmiles())
}

/// Parse extended reaction SMILES bytes with configuration
pub fn parse_extended_reaction_smiles_bytes_with(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<ExtendedReaction, ParseError> {
    let flags = config.parse_flags;

    // Check if the input contains a CX block, record ring bonds if it is present.
    let has_cx_annotations =
        flags.contains(SmilesParseFlags::CHEMAXON_EXTENSIONS) && input.contains(&b'|');
    let mut remaining = input;
    let mut offset = 0usize;

    if remaining.starts_with(b".") {
        return Err(ParseError::LeadingDot { pos: 0 });
    }

    if !input.is_empty() && input[0].is_ascii_whitespace() {
        return Err(ParseError::LeadingWhitespace);
    }

    // Reactants: parse one side-supermolecule until '>'.
    let (rest, (reactants, reactant_ring_bonds, new_offset)) =
        parse_extended_smiles_inner(remaining, offset, true, has_cx_annotations, flags)?;
    offset = new_offset;
    remaining = rest;

    let agents;
    let agent_ring_bonds;
    if remaining.starts_with(b">>") {
        remaining = &remaining[2..];
        offset += 2;
        agents = ExtendedMolecule::empty();
        agent_ring_bonds = Vec::new();
    } else if remaining.starts_with(b">") {
        remaining = &remaining[1..];
        offset += 1;

        let (rest, (agents_parsed, agents_ring_bonds, new_offset)) =
            parse_extended_smiles_inner(remaining, offset, true, has_cx_annotations, flags)?;
        offset = new_offset;
        remaining = rest;
        agents = agents_parsed;
        agent_ring_bonds = agents_ring_bonds;

        if !remaining.starts_with(b">") {
            return Err(ParseError::MissingReactionArrow { pos: offset });
        }
        remaining = &remaining[1..];
        offset += 1;
    } else {
        return Err(ParseError::MissingReactionArrow { pos: offset });
    }

    // Products: parse one side-supermolecule until EOF/whitespace.
    let (rest, (products, product_ring_bonds, _new_offset)) =
        parse_extended_smiles_inner(remaining, offset, true, has_cx_annotations, flags)?;

    let mut reaction = ExtendedReaction {
        reactants,
        products,
        agents,
        atom_mapping: BTreeMap::new(),
        comments: Vec::new(),
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    if has_cx_annotations {
        let entries = parse_extended_cx_annotations(rest.trim_ascii_start(), flags)?;
        let mut split = split_reaction_cx_entries(
            entries,
            reaction.reactants.atom_count(),
            reaction.reactants.bond_count(),
            reaction.agents.atom_count(),
            reaction.agents.bond_count(),
            reaction.products.atom_count(),
            reaction.products.bond_count(),
        )?;
        remap_cx_bond_indices(
            &mut split.0,
            &BondIndexMap::new(reactant_ring_bonds, reaction.reactants.bond_count()),
        )?;
        remap_cx_bond_indices(
            &mut split.1,
            &BondIndexMap::new(agent_ring_bonds, reaction.agents.bond_count()),
        )?;
        remap_cx_bond_indices(
            &mut split.2,
            &BondIndexMap::new(product_ring_bonds, reaction.products.bond_count()),
        )?;
        update_extended_reaction(&mut reaction, split)?;
    }
    collect_extended_atom_mapping(&mut reaction);
    Ok(reaction)
}

#[allow(clippy::type_complexity)]
fn parse_extended_smiles_inner(
    input: &[u8],
    offset: usize,
    as_reaction: bool,
    store_rings: bool,
    flags: SmilesParseFlags,
) -> Result<(&[u8], (ExtendedMolecule, Vec<(usize, usize)>, usize)), ParseError> {
    let extended_bonds = flags.contains(SmilesParseFlags::EXTENDED_BONDS);
    let mut i = 0usize;
    let n = input.len();
    let mut builder =
        ExtendedMoleculeBuilder::with_capacity(n.max(1), n.max(1).saturating_sub(1), store_rings);
    let mut branch_stack: Vec<Frame> = Vec::new();
    let mut last_atom_idx: Option<usize> = None;
    let mut pending_bond: Option<(BondOrder, Option<BondDirection>, Option<BondDonation>, usize)> =
        None;
    let mut just_closed_group: bool = false;

    while i < n {
        let b0 = input[i];

        // Stop at whitespace - return remaining input
        if b0.is_ascii_whitespace() {
            break;
        }
        if as_reaction && b0 == b'>' {
            break;
        }

        if b0 != b'(' {
            just_closed_group = false;
        }
        if b0 == b'(' {
            if let Some((_, _, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos: offset + pos });
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
            if let Some((_, _, _, pos)) = pending_bond {
                return Err(ParseError::TrailingBond { pos: offset + pos });
            }
            let Some(frame) = branch_stack.pop() else {
                return Err(ParseError::UnbalancedCloseParen { pos: offset + i });
            };
            match frame {
                Frame::Branch { base, had_atom, .. } => {
                    if !had_atom {
                        return Err(ParseError::EmptyBranch { pos: offset + i });
                    }
                    last_atom_idx = Some(base);
                }
                Frame::Group { had_atom, .. } => {
                    if !had_atom {
                        return Err(ParseError::EmptyGroup { pos: offset + i });
                    }
                    just_closed_group = true;
                    if let Some(parent) = branch_stack.last_mut() {
                        match parent {
                            Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. } => {
                                *had_atom = true
                            }
                        }
                    }
                    if branch_stack.is_empty() && i + 1 != n {
                        let next = input[i + 1];
                        if next != b'.' {
                            return Err(ParseError::NonfinalGroup { pos: offset + i });
                        }
                    }
                }
            }
            i += 1;
            continue;
        }
        if b0 == b'.' {
            if let Some((_, _, _, pos)) = pending_bond {
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
            if as_reaction && input[i + 1] == b'>' {
                return Err(ParseError::TrailingDot { pos: offset + i });
            }
            if input[i + 1] == b'.' {
                return Err(ParseError::ConsecutiveDots { pos: offset + i });
            }
            if input[i + 1].is_ascii_digit() {
                return Err(ParseError::DotBeforeRing { pos: offset + i });
            }
            if input[i + 1] == b'%' {
                return Err(ParseError::DotBeforeRing { pos: offset + i });
            }
            last_atom_idx = None;
            i += 1;
            continue;
        }
        match parse_ring_index(input, i, offset) {
            Ok(Some((idx, next_i, _percent))) => {
                if last_atom_idx.is_none() {
                    return Err(ParseError::LeadingRing { pos: offset + i });
                }
                if invalid_ring_context(&branch_stack) {
                    return Err(ParseError::LeadingRing { pos: offset });
                }
                let bond = pending_bond.take();
                let (order_opt, direction_opt, donation_opt) =
                    bond.map_or((None, None, None), |(o, d, don, _)| (Some(o), d, don));
                builder.on_ring_bond(
                    last_atom_idx.unwrap(),
                    idx,
                    order_opt,
                    direction_opt,
                    donation_opt,
                    i,
                    i + 1,
                    offset,
                )?;
                i = next_i;
                continue;
            }
            Err(e) => return Err(e),
            Ok(None) => {}
        }
        // Extended bonds: also match ~ and < (for <- dative)
        if matches!(b0, b'-' | b'=' | b'#' | b'$' | b':' | b'/' | b'\\')
            || (extended_bonds && matches!(b0, b'~' | b'<'))
        {
            if pending_bond.is_some() {
                return Err(ParseError::ConsecutiveBonds { pos: offset + i });
            }
            if last_atom_idx.is_none() {
                if let Some(Frame::Group {
                    had_atom: false, ..
                }) = branch_stack.last()
                {
                    return Err(ParseError::LeadingBond { pos: offset + i });
                }
                return Err(ParseError::LeadingBond { pos: offset + i });
            }
            // Use extended bond parsing for ->, <-, ~ when EXTENDED_BONDS is set
            if extended_bonds {
                let (order, direction, donation, consumed) = parse_extended_bond(input, i);
                pending_bond = Some((order, direction, donation, i));
                i += consumed;
            } else {
                let (order, bond_direction) = parse_bond(b0);
                pending_bond = Some((order, bond_direction, None, i));
                i += 1;
            }
            continue;
        }
        if b0 == b'[' {
            let start = i + 1;
            let mut j = start;
            while j < n && input[j] != b']' {
                j += 1;
            }
            if j >= n {
                return Err(ParseError::UnbalancedOpenBracket { pos: offset + i });
            }
            if j == start {
                return Err(ParseError::EmptyBracket { pos: offset + i });
            }
            let inner = &input[start..j];
            let (symbol, aromatic, iso_opt, charge_opt, class_opt, h_opt, chir_opt) =
                parse_extended_bracket(inner, offset + i, flags)?;
            let (s, e) = (Some(i as u32), Some((j + 1) as u32));
            let atom = ExtendedAtomData {
                symbol,
                isotope: iso_opt,
                charge: charge_opt,
                implicit_hydrogens: Some(h_opt.unwrap_or(0)),
                class: class_opt,
                aromatic,
                chirality: chir_opt,
                span: Span::from_bytes_opt(s, e),
            };
            let curr = builder.on_atom(atom);

            attach_extended_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                aromatic,
                i as u32,
                (j + 1) as u32,
            );
            last_atom_idx = Some(curr);
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

                attach_extended_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    false,
                    i as u32,
                    (i + 2) as u32,
                );
                last_atom_idx = Some(curr);
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

            attach_extended_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
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

                attach_extended_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    false,
                    i as u32,
                    (i + 2) as u32,
                );
                last_atom_idx = Some(curr);
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

            attach_extended_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
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

                attach_extended_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    false,
                    i as u32,
                    (i + consumed) as u32,
                );
                last_atom_idx = Some(curr);
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

                attach_extended_atom(
                    &mut builder,
                    last_atom_idx,
                    curr,
                    &mut pending_bond,
                    true,
                    i as u32,
                    (i + consumed) as u32,
                );
                last_atom_idx = Some(curr);
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
            return Err(ParseError::InvalidElement { pos: offset + i });
        }
        if b0 == b'*' {
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            let curr = builder.on_wildcard(WildcardAtom::Any, None, s, e);

            attach_extended_atom(
                &mut builder,
                last_atom_idx,
                curr,
                &mut pending_bond,
                false,
                i as u32,
                (i + 1) as u32,
            );
            last_atom_idx = Some(curr);
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
            return Err(ParseError::UnbalancedCloseBracket { pos: offset + i });
        }
        if b0 == b'@' || b0 == b'+' {
            return Err(ParseError::StrayBracketField { pos: offset + i });
        }
        return Err(ParseError::InvalidToken { pos: offset + i });
    }

    if let Some((_, _, _, pos)) = pending_bond {
        return Err(ParseError::TrailingBond { pos: offset + pos });
    }
    if !branch_stack.is_empty() {
        let pos = match branch_stack.last().unwrap() {
            Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. } => *open_pos,
        };
        return Err(ParseError::UnbalancedOpenParen { pos: offset + pos });
    }
    if let Some(pos_open) = builder.unclosed_ring_pos() {
        return Err(ParseError::UnbalancedRingIndex {
            open_pos: offset + pos_open,
        });
    }
    let mut mols = builder.finish();
    let mol = mols
        .pop()
        .map(|mut mol| {
            mol.source_format = SourceFormat::SMILES;
            mol.chirality_frame = Some(ChiralityFrame::FirstNeighborToward);
            mol
        })
        .unwrap_or_else(ExtendedMolecule::empty);
    let ring_bonds = builder.take_ring_bonds();
    let new_offset = offset + i;
    Ok((&input[i..], (mol, ring_bonds, new_offset)))
}

fn collect_atom_mapping(reaction: &mut Reaction) {
    for (at_idx, atom) in reaction.reactants.atoms.iter().enumerate() {
        if let Some(class) = atom.class {
            reaction
                .atom_mapping
                .entry(class)
                .or_default()
                .0
                .push(at_idx as u32);
        }
    }

    for (at_idx, atom) in reaction.products.atoms.iter().enumerate() {
        if let Some(class) = atom.class {
            reaction
                .atom_mapping
                .entry(class)
                .or_default()
                .1
                .push(at_idx as u32);
        }
    }
}

fn collect_extended_atom_mapping(reaction: &mut ExtendedReaction) {
    for (at_idx, atom) in reaction.reactants.atoms.iter().enumerate() {
        if let Some(class) = atom.class {
            reaction
                .atom_mapping
                .entry(class)
                .or_default()
                .0
                .push(at_idx as u32);
        }
    }

    for (at_idx, atom) in reaction.products.atoms.iter().enumerate() {
        if let Some(class) = atom.class {
            reaction
                .atom_mapping
                .entry(class)
                .or_default()
                .1
                .push(at_idx as u32);
        }
    }
}

#[cfg(test)]
mod tests;
