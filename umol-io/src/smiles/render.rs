//! SMILES output traversal and formatting.

use std::collections::BTreeMap;
use std::fmt::Write;

use smallvec::SmallVec;
use umol_chem::element::Element;
use umol_chem::spin::SpinMultiplicity;
use umol_perm::Permutation;

use self::stereo::{assign_markers, MarkerAssignmentError, MarkerComponentError};
use self::traversal::Traversal;
use super::config::{SmilesIoConfig, SmilesSyntaxFlags};
use crate::table_ir::{
    Atom, BondConfiguration, BondDirection, BondDonation, BondOrder, Chirality, Molecule,
    StereoAtom, StereoLigand, Winding,
};

mod stereo;
mod traversal;

#[derive(Debug, PartialEq, Eq)]
pub(super) enum RenderError {
    AtomIndexOutOfBounds { atom: u32 },
    BondIndexOutOfBounds { bond: u32 },
    SelfBond { bond: u32 },
    DuplicateBond { bond: u32 },
    UnsupportedMolecule { field: &'static str },
    UnsupportedAtom { atom: u32, field: &'static str },
    UnsupportedBond { bond: u32, field: &'static str },
    InferredHydrogens { atom: u32 },
    DuplicateStereoAtom { atom: u32 },
    InvalidStereoAtom { atom: u32 },
    RingLabel { label: usize },
    DuplicateStereoBond { bond: u32 },
    UnsupportedStereoBond { bond: u32 },
    MissingMarkerCandidate { bond: u32, atom: u32 },
    InvalidStereoBondReference { bond: u32, atom: u32 },
    NoMarkerAssignment { bond: u32 },
}

impl From<MarkerAssignmentError> for RenderError {
    fn from(error: MarkerAssignmentError) -> Self {
        match error {
            MarkerAssignmentError::InvalidReference { bond, atom } => {
                Self::InvalidStereoBondReference { bond, atom }
            }
            MarkerAssignmentError::NoAssignment { bond } => Self::NoMarkerAssignment { bond },
            MarkerAssignmentError::Component(error) => match error {
                MarkerComponentError::BondIndexOutOfBounds { bond } => {
                    Self::BondIndexOutOfBounds { bond }
                }
                MarkerComponentError::AtomIndexOutOfBounds { atom } => {
                    Self::AtomIndexOutOfBounds { atom }
                }
                MarkerComponentError::DuplicateStereoBond { bond } => {
                    Self::DuplicateStereoBond { bond }
                }
                MarkerComponentError::UnsupportedSite { bond } => {
                    Self::UnsupportedStereoBond { bond }
                }
                MarkerComponentError::MissingCandidate { bond, atom } => {
                    Self::MissingMarkerCandidate { bond, atom }
                }
            },
        }
    }
}

/// Formats the supported table fields without chemical resolution or source-spelling recovery.
pub(super) fn render(molecule: &Molecule, config: &SmilesIoConfig) -> Result<String, RenderError> {
    for (present, field) in [
        (molecule.positions.is_some(), "positions"),
        (!molecule.multicenter_bonds.is_empty(), "multicenter_bonds"),
        (
            molecule.configuration_scope.is_some(),
            "configuration_scope",
        ),
        (!molecule.comments.is_empty(), "comments"),
        (!molecule.properties.is_empty(), "properties"),
    ] {
        if present {
            return Err(RenderError::UnsupportedMolecule { field });
        }
    }
    for (index, bond) in molecule.bonds.iter().enumerate() {
        let bond_index = index as u32;
        for atom in [bond.atoms.first(), bond.atoms.second()] {
            if atom as usize >= molecule.atoms.len() {
                return Err(RenderError::AtomIndexOutOfBounds { atom });
            }
        }
        if bond.atoms.first() == bond.atoms.second() {
            return Err(RenderError::SelfBond { bond: bond_index });
        }
        for (present, field) in [
            (bond.charge.is_some_and(|charge| charge != 0), "charge"),
            (
                bond.unpaired_electrons.is_some_and(|count| count != 0),
                "unpaired_electrons",
            ),
            (
                bond.multiplicity
                    .is_some_and(|spin| spin != SpinMultiplicity::SINGLET),
                "multiplicity",
            ),
            (bond.noncovalent.is_some(), "noncovalent"),
            (bond.wedge.is_some(), "wedge"),
        ] {
            if present {
                return Err(RenderError::UnsupportedBond {
                    bond: bond_index,
                    field,
                });
            }
        }
    }
    for frame in &molecule.stereo_bonds {
        if frame.configuration == BondConfiguration::Either {
            return Err(RenderError::UnsupportedBond {
                bond: frame.bond,
                field: "Either",
            });
        }
    }
    let mut frames = BTreeMap::new();
    for frame in &molecule.stereo_atoms {
        if frame.atom as usize >= molecule.atoms.len() {
            return Err(RenderError::AtomIndexOutOfBounds { atom: frame.atom });
        }
        if frames.insert(frame.atom, frame).is_some() {
            return Err(RenderError::DuplicateStereoAtom { atom: frame.atom });
        }
    }
    let traversal = Traversal::new(molecule);
    let markers = assign_markers(molecule, &traversal)?;
    let mut output = String::with_capacity(molecule.atoms.len());
    let mut ends = Vec::new();
    for (position, visit) in traversal.atoms.iter().enumerate() {
        while ends.last() == Some(&position) {
            ends.pop();
            if !ends.is_empty() {
                output.push(')');
            }
        }
        if visit.parent.is_none() {
            if position != 0 {
                output.push('.');
            }
            ends.push(visit.subtree_end);
        } else if ends.last().is_some_and(|&end| visit.subtree_end < end) {
            output.push('(');
            ends.push(visit.subtree_end);
        }
        if let Some(parent) = visit.parent {
            append_bond(
                &mut output,
                molecule,
                parent.bond,
                parent.atom,
                &markers,
                config,
            )?;
        }
        let mut neighbors = SmallVec::<[u32; 4]>::new();
        for neighbor in traversal.neighbors(position) {
            if neighbors.contains(&neighbor.atom) {
                return Err(RenderError::DuplicateBond {
                    bond: neighbor.bond,
                });
            }
            neighbors.push(neighbor.atom);
        }
        let winding = frames
            .get(&visit.atom)
            .map(|frame| {
                atom_winding(
                    frame,
                    &molecule.atoms[visit.atom as usize],
                    &neighbors,
                    visit.parent.is_some(),
                )
            })
            .transpose()?;
        append_atom(
            &mut output,
            visit.atom,
            &molecule.atoms[visit.atom as usize],
            winding,
            config,
        )?;
        for ring in &traversal.rings[visit.rings.clone()] {
            if ring.opening {
                append_bond(
                    &mut output,
                    molecule,
                    ring.neighbor.bond,
                    visit.atom,
                    &markers,
                    config,
                )?;
            }
            if ring.label > 99 {
                return Err(RenderError::RingLabel { label: ring.label });
            }
            if ring.label >= 10 {
                output.push('%');
            }
            write!(output, "{}", ring.label).expect("String write");
        }
    }
    while ends.pop().is_some() {
        if !ends.is_empty() {
            output.push(')');
        }
    }
    Ok(output)
}

fn atom_winding(
    frame: &StereoAtom,
    atom: &Atom,
    neighbors: &[u32],
    has_parent: bool,
) -> Result<Winding, RenderError> {
    let error = || RenderError::InvalidStereoAtom { atom: frame.atom };
    if frame.ligands.len() != 4 {
        return Err(error());
    }
    let mut ligands: SmallVec<[StereoLigand; 4]> =
        neighbors.iter().copied().map(StereoLigand::Atom).collect();
    match (neighbors.len(), atom.implicit_hydrogens) {
        (4, Some(0)) => {}
        (3, Some(1)) => ligands.insert(usize::from(has_parent), StereoLigand::ImplicitHydrogen),
        (3, Some(0)) => ligands.insert(usize::from(has_parent), StereoLigand::LonePair),
        _ => return Err(error()),
    }
    let permutation = Permutation::between(&frame.ligands, &ligands).ok_or_else(error)?;
    Ok(match (frame.winding, permutation.sign() < 0) {
        (Winding::Clockwise, false) | (Winding::CounterClockwise, true) => Winding::Clockwise,
        _ => Winding::CounterClockwise,
    })
}

fn append_atom(
    output: &mut String,
    index: u32,
    atom: &Atom,
    winding: Option<Winding>,
    config: &SmilesIoConfig,
) -> Result<(), RenderError> {
    let unsupported = |field| RenderError::UnsupportedAtom { atom: index, field };
    for (present, field) in [
        (atom.valence.is_some(), "valence"),
        (atom.lone_pairs.is_some(), "lone_pairs"),
        (atom.unpaired_electrons.is_some(), "unpaired_electrons"),
        (atom.multiplicity.is_some(), "multiplicity"),
        (atom.label.is_some(), "label"),
        (atom.value.is_some(), "value"),
    ] {
        if present {
            return Err(unsupported(field));
        }
    }
    match atom.chirality {
        None => {}
        Some(
            Chirality::Clockwise
            | Chirality::CounterClockwise
            | Chirality::Tetrahedral { arr: 1 | 2 },
        ) if winding.is_some() => {}
        _ => return Err(unsupported("chirality")),
    }
    if atom
        .implicit_hydrogens
        .is_some_and(|h| h > 9 || (h > 0 && atom.element == Some(Element::H)))
    {
        return Err(unsupported("implicit_hydrogens"));
    }
    if atom.charge.is_some_and(|charge| charge.unsigned_abs() > 99) {
        return Err(unsupported("charge"));
    }
    let aromatic = atom.aromatic == Some(true);
    let symbol = if aromatic {
        match atom.element {
            Some(Element::B) => "b",
            Some(Element::C) => "c",
            Some(Element::N) => "n",
            Some(Element::O) => "o",
            Some(Element::P) => "p",
            Some(Element::S) => "s",
            Some(Element::As) => "as",
            Some(Element::Se) => "se",
            Some(Element::Si)
                if config
                    .syntax_flags
                    .contains(SmilesSyntaxFlags::EXTENDED_AROMATICS) =>
            {
                "si"
            }
            Some(Element::Te)
                if config
                    .syntax_flags
                    .contains(SmilesSyntaxFlags::EXTENDED_AROMATICS) =>
            {
                "te"
            }
            _ => return Err(unsupported("aromatic")),
        }
    } else {
        atom.element.map_or("*", |element| element.symbol())
    };
    let organic = if aromatic {
        matches!(
            atom.element,
            Some(Element::B | Element::C | Element::N | Element::O | Element::P | Element::S)
        )
    } else {
        matches!(
            atom.element,
            None | Some(
                Element::B
                    | Element::C
                    | Element::N
                    | Element::O
                    | Element::P
                    | Element::S
                    | Element::F
                    | Element::Cl
                    | Element::Br
                    | Element::I
            )
        )
    };
    let bracket = atom.implicit_hydrogens.is_some()
        || atom.isotope_mass.is_some()
        || atom.charge.is_some_and(|charge| charge != 0)
        || atom.class.is_some()
        || winding.is_some()
        || !organic;
    if !bracket {
        output.push_str(symbol);
        return Ok(());
    }
    let hydrogens = atom
        .implicit_hydrogens
        .ok_or(RenderError::InferredHydrogens { atom: index })?;
    output.push('[');
    if let Some(mass) = atom.isotope_mass {
        write!(output, "{mass}").expect("String write");
    }
    output.push_str(symbol);
    if let Some(winding) = winding {
        output.push_str(if winding == Winding::Clockwise {
            "@@"
        } else {
            "@"
        });
    }
    if hydrogens > 0 {
        output.push('H');
        if hydrogens > 1 {
            write!(output, "{hydrogens}").expect("String write");
        }
    }
    if let Some(charge) = atom.charge.filter(|&charge| charge != 0) {
        output.push(if charge > 0 { '+' } else { '-' });
        if charge.unsigned_abs() > 1 {
            write!(output, "{}", charge.unsigned_abs()).expect("String write");
        }
    }
    if let Some(class) = atom.class {
        write!(output, ":{class}").expect("String write");
    }
    output.push(']');
    Ok(())
}

fn append_bond(
    output: &mut String,
    molecule: &Molecule,
    index: u32,
    from: u32,
    markers: &[(u32, BondDirection)],
    config: &SmilesIoConfig,
) -> Result<(), RenderError> {
    let bond = &molecule.bonds[index as usize];
    let unsupported = |field| RenderError::UnsupportedBond { bond: index, field };
    if let Some(donation) = bond.donation {
        if !config
            .syntax_flags
            .contains(SmilesSyntaxFlags::EXTENDED_BONDS)
            || bond.order != BondOrder::Single
        {
            return Err(unsupported("donation"));
        }
        let donation = if bond.atoms.first() == from {
            donation
        } else {
            donation.flip()
        };
        output.push_str(match donation {
            BondDonation::Donating => "->",
            BondDonation::Accepting => "<-",
            BondDonation::Shared => return Err(unsupported("donation")),
        });
        return Ok(());
    }
    if let Ok(position) = markers.binary_search_by_key(&index, |&(bond, _)| bond) {
        let direction = if bond.atoms.first() == from {
            markers[position].1
        } else {
            markers[position].1.flip()
        };
        output.push(if direction == BondDirection::Rising {
            '/'
        } else {
            '\\'
        });
        return Ok(());
    }
    let both_aromatic = [bond.atoms.first(), bond.atoms.second()]
        .iter()
        .all(|&atom| molecule.atoms[atom as usize].aromatic == Some(true));
    output.push_str(match bond.order {
        BondOrder::Single if both_aromatic => "-",
        BondOrder::Single => "",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Quadruple => "$",
        BondOrder::Aromatic if both_aromatic => "",
        BondOrder::Aromatic => ":",
        BondOrder::Any
            if config
                .syntax_flags
                .contains(SmilesSyntaxFlags::EXTENDED_BONDS) =>
        {
            "~"
        }
        _ => return Err(unsupported("order")),
    });
    Ok(())
}

#[cfg(test)]
mod tests;
