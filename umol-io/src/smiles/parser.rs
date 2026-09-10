//! SMILES parser

use std::collections::BTreeMap;

use indexmap::IndexMap;
use umol_chem::element::Element;

mod builder;
mod cx;
mod utils;

use self::builder::{
    AtomData, AtomMapping, ExtendedAtomData, ExtendedMoleculeBuilder, MoleculeEditor,
};
use self::cx::{
    parse_cx_annotations, parse_extended_cx_annotations, remap_cx_bond_indices,
    split_reaction_cx_entries, update_extended_molecule, update_extended_reaction, update_molecule,
    update_reaction, BondIndexMap,
};
use self::utils::{
    parse_bond, parse_bracket, parse_extended_bond, parse_extended_bracket,
    parse_organic_aliphatic_element, parse_organic_aromatic_element, parse_ring_index,
};
use super::config::{SmilesIoConfig, SmilesSyntaxFlags};
use super::error::ParseError;
use crate::table_ir::{ExtendedMolecule, ExtendedReaction, Molecule, Reaction, SourceFormat, Span};

/// Parse a molecular SMILES byte slice into `table_ir::Molecule`.
pub(crate) fn parse_molecule(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<Molecule, ParseError> {
    let flags = config.syntax_flags;
    if input.is_empty() {
        return Ok(Molecule::empty());
    }

    // Leading whitespace is not allowed (exception: whitespace-only input is allowed)
    if !input.is_empty() && input[0].is_ascii_whitespace() && !input.trim_ascii_start().is_empty() {
        return Err(ParseError::LeadingWhitespace);
    }

    // Check if the input contains a CX block, record ring bonds if it is present.
    let has_cx_annotations =
        flags.contains(SmilesSyntaxFlags::CHEMAXON_EXTENSIONS) && input.contains(&b'|');
    let (remaining, (mut mol, ring_bonds, _)) =
        parse_smiles_inner(input, 0, false, has_cx_annotations, flags, None)?;

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

/// Parse a reaction SMILES byte slice into `table_ir::Reaction`.
pub(crate) fn parse_reaction(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<Reaction, ParseError> {
    let flags = config.syntax_flags;
    // Check if the input contains a CX block, record ring bonds if it is present.
    let has_cx_annotations =
        flags.contains(SmilesSyntaxFlags::CHEMAXON_EXTENSIONS) && input.contains(&b'|');
    let mut atom_mapping = BTreeMap::new();
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
    let (rest, (reactants, reactant_ring_bonds, new_offset)) = parse_smiles_inner(
        remaining,
        offset,
        true,
        has_cx_annotations,
        flags,
        Some((&mut atom_mapping, false)),
    )?;
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
            parse_smiles_inner(remaining, offset, true, has_cx_annotations, flags, None)?;
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
    let (rest, (products, product_ring_bonds, _new_offset)) = parse_smiles_inner(
        remaining,
        offset,
        true,
        has_cx_annotations,
        flags,
        Some((&mut atom_mapping, true)),
    )?;

    let mut reaction = Reaction {
        reactants,
        products,
        agents,
        atom_mapping,
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
    Ok(reaction)
}

#[allow(clippy::type_complexity)]
fn parse_smiles_inner<'a>(
    input: &'a [u8],
    offset: usize,
    as_reaction: bool,
    store_rings: bool,
    flags: SmilesSyntaxFlags,
    mapping: Option<(&mut AtomMapping, bool)>,
) -> Result<(&'a [u8], (Molecule, Vec<(usize, usize)>, usize)), ParseError> {
    let extended_bonds = flags.contains(SmilesSyntaxFlags::EXTENDED_BONDS);
    let mut i = 0usize;
    let n = input.len();
    let mut builder =
        MoleculeEditor::with_capacity(n.max(1), n.max(1).saturating_sub(1), store_rings, mapping);
    while i < n {
        let b0 = input[i];

        // Stop at whitespace - return remaining input
        if b0.is_ascii_whitespace() {
            break;
        }
        if as_reaction && b0 == b'>' {
            break;
        }

        builder.token(b0);
        if b0 == b'(' {
            builder.open_branch(i, offset)?;
            i += 1;
            continue;
        }
        if b0 == b')' {
            builder.close_branch(input, i, offset)?;
            i += 1;
            continue;
        }
        if b0 == b'.' {
            builder.dot(input, i, offset, as_reaction)?;
            i += 1;
            continue;
        }
        match parse_ring_index(input, i, offset) {
            Ok(Some((idx, next_i, _percent))) => {
                builder.on_ring_bond(idx, i, i + 1, offset)?;
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
            // Use extended bond parsing for ->, <-, ~ when EXTENDED_BONDS is set
            if extended_bonds {
                let (order, direction, donation, consumed) = parse_extended_bond(input, i);
                builder.on_bond(order, direction, donation, i, offset)?;
                i += consumed;
            } else {
                let (order, bond_direction) = parse_bond(b0);
                builder.on_bond(order, bond_direction, None, i, offset)?;
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
            builder.on_atom(atom);
            i = j + 1;
            continue;
        }
        if b0 == b'C' {
            if i + 1 < n && input[i + 1] == b'l' {
                let (s, e) = (Some(i as u32), Some((i + 2) as u32));
                builder.on_atom_fast(Element::Cl, false, s, e);
                i += 2;
                continue;
            }
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            builder.on_atom_fast(Element::C, false, s, e);
            i += 1;
            continue;
        }
        if b0 == b'B' {
            if i + 1 < n && input[i + 1] == b'r' {
                let (s, e) = (Some(i as u32), Some((i + 2) as u32));
                builder.on_atom_fast(Element::Br, false, s, e);
                i += 2;
                continue;
            }
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            builder.on_atom_fast(Element::B, false, s, e);
            i += 1;
            continue;
        }
        // Elements
        if b0.is_ascii_alphabetic() {
            if let Some((element, consumed)) = parse_organic_aliphatic_element(input, i) {
                let (s, e) = (Some(i as u32), Some((i + consumed) as u32));
                builder.on_atom_fast(element, false, s, e);
                i += consumed;
                continue;
            }
            if let Some((element, consumed)) = parse_organic_aromatic_element(input, i) {
                let (s, e) = (Some(i as u32), Some((i + consumed) as u32));
                builder.on_atom_fast(element, true, s, e);
                i += consumed;
                continue;
            }
            return Err(ParseError::InvalidElement { pos: offset + i });
        }
        if b0 == b'*' {
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            builder.on_wildcard(s, e);
            i += 1;
            continue;
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

    let (mol, ring_bonds) = builder.finish(offset)?;
    let new_offset = offset + i;
    Ok((&input[i..], (mol, ring_bonds, new_offset)))
}

/// Parse extended SMILES text with the OpenSMILES configuration.
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

/// Parse extended SMILES bytes with the OpenSMILES configuration.
pub fn parse_extended_smiles_bytes(input: &[u8]) -> Result<ExtendedMolecule, ParseError> {
    parse_extended_smiles_bytes_with(input, &SmilesIoConfig::opensmiles())
}

/// Parse extended SMILES bytes with configuration
pub fn parse_extended_smiles_bytes_with(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<ExtendedMolecule, ParseError> {
    let flags = config.syntax_flags;

    if input.is_empty() {
        return Ok(ExtendedMolecule::empty());
    }

    // Check if the input contains a CX block, record ring bonds if it is present.
    let has_cx_annotations =
        flags.contains(SmilesSyntaxFlags::CHEMAXON_EXTENSIONS) && input.contains(&b'|');
    let (remaining, (mut mol, ring_bonds, _)) =
        parse_extended_smiles_inner(input, 0, false, has_cx_annotations, flags, None)?;

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

/// Parse extended reaction SMILES text with the OpenSMILES configuration.
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

/// Parse extended reaction SMILES bytes with the OpenSMILES configuration.
pub fn parse_extended_reaction_smiles_bytes(input: &[u8]) -> Result<ExtendedReaction, ParseError> {
    parse_extended_reaction_smiles_bytes_with(input, &SmilesIoConfig::opensmiles())
}

/// Parse extended reaction SMILES bytes with configuration
pub fn parse_extended_reaction_smiles_bytes_with(
    input: &[u8],
    config: &SmilesIoConfig,
) -> Result<ExtendedReaction, ParseError> {
    let flags = config.syntax_flags;

    // Check if the input contains a CX block, record ring bonds if it is present.
    let has_cx_annotations =
        flags.contains(SmilesSyntaxFlags::CHEMAXON_EXTENSIONS) && input.contains(&b'|');
    let mut atom_mapping = BTreeMap::new();
    let mut remaining = input;
    let mut offset = 0usize;

    if remaining.starts_with(b".") {
        return Err(ParseError::LeadingDot { pos: 0 });
    }

    if !input.is_empty() && input[0].is_ascii_whitespace() {
        return Err(ParseError::LeadingWhitespace);
    }

    // Reactants: parse one side-supermolecule until '>'.
    let (rest, (reactants, reactant_ring_bonds, new_offset)) = parse_extended_smiles_inner(
        remaining,
        offset,
        true,
        has_cx_annotations,
        flags,
        Some((&mut atom_mapping, false)),
    )?;
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
            parse_extended_smiles_inner(remaining, offset, true, has_cx_annotations, flags, None)?;
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
    let (rest, (products, product_ring_bonds, _new_offset)) = parse_extended_smiles_inner(
        remaining,
        offset,
        true,
        has_cx_annotations,
        flags,
        Some((&mut atom_mapping, true)),
    )?;

    let mut reaction = ExtendedReaction {
        reactants,
        products,
        agents,
        atom_mapping,
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
    Ok(reaction)
}

#[allow(clippy::type_complexity)]
fn parse_extended_smiles_inner<'a>(
    input: &'a [u8],
    offset: usize,
    as_reaction: bool,
    store_rings: bool,
    flags: SmilesSyntaxFlags,
    mapping: Option<(&mut AtomMapping, bool)>,
) -> Result<(&'a [u8], (ExtendedMolecule, Vec<(usize, usize)>, usize)), ParseError> {
    let extended_bonds = flags.contains(SmilesSyntaxFlags::EXTENDED_BONDS);
    let mut i = 0usize;
    let n = input.len();
    let mut builder = ExtendedMoleculeBuilder::with_capacity(
        n.max(1),
        n.max(1).saturating_sub(1),
        store_rings,
        mapping,
    );
    while i < n {
        let b0 = input[i];

        // Stop at whitespace - return remaining input
        if b0.is_ascii_whitespace() {
            break;
        }
        if as_reaction && b0 == b'>' {
            break;
        }

        builder.token(b0);
        if b0 == b'(' {
            builder.open_branch(i, offset)?;
            i += 1;
            continue;
        }
        if b0 == b')' {
            builder.close_branch(input, i, offset)?;
            i += 1;
            continue;
        }
        if b0 == b'.' {
            builder.dot(input, i, offset, as_reaction)?;
            i += 1;
            continue;
        }
        match parse_ring_index(input, i, offset) {
            Ok(Some((idx, next_i, _percent))) => {
                builder.on_ring_bond(idx, i, i + 1, offset)?;
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
            // Use extended bond parsing for ->, <-, ~ when EXTENDED_BONDS is set
            if extended_bonds {
                let (order, direction, donation, consumed) = parse_extended_bond(input, i);
                builder.on_bond(order, direction, donation, i, offset)?;
                i += consumed;
            } else {
                let (order, bond_direction) = parse_bond(b0);
                builder.on_bond(order, bond_direction, None, i, offset)?;
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
            builder.on_atom(atom);
            i = j + 1;
            continue;
        }
        if b0 == b'C' {
            if i + 1 < n && input[i + 1] == b'l' {
                let (s, e) = (Some(i as u32), Some((i + 2) as u32));
                builder.on_atom_fast(Element::Cl, false, s, e);
                i += 2;
                continue;
            }
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            builder.on_atom_fast(Element::C, false, s, e);
            i += 1;
            continue;
        }
        if b0 == b'B' {
            if i + 1 < n && input[i + 1] == b'r' {
                let (s, e) = (Some(i as u32), Some((i + 2) as u32));
                builder.on_atom_fast(Element::Br, false, s, e);
                i += 2;
                continue;
            }
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            builder.on_atom_fast(Element::B, false, s, e);
            i += 1;
            continue;
        }
        if b0.is_ascii_alphabetic() {
            if let Some((element, consumed)) = parse_organic_aliphatic_element(input, i) {
                let (s, e) = (Some(i as u32), Some((i + consumed) as u32));
                builder.on_atom_fast(element, false, s, e);
                i += consumed;
                continue;
            }
            if let Some((element, consumed)) = parse_organic_aromatic_element(input, i) {
                let (s, e) = (Some(i as u32), Some((i + consumed) as u32));
                builder.on_atom_fast(element, true, s, e);
                i += consumed;
                continue;
            }
            return Err(ParseError::InvalidElement { pos: offset + i });
        }
        if b0 == b'*' {
            let (s, e) = (Some(i as u32), Some((i + 1) as u32));
            builder.on_wildcard(s, e);
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

    let (mol, ring_bonds) = builder.finish(offset)?;
    let new_offset = offset + i;
    Ok((&input[i..], (mol, ring_bonds, new_offset)))
}

#[cfg(test)]
mod tests;
