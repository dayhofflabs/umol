//! Molecule assembly from CTab atom and bond data

use crate::atom::{Atom, AtomLike, AtomSymbol};
use crate::bond::Bond;
use crate::conformer::{Conformer, Point3D};
use crate::molecule::Molecule;
use nom::error;
use nom::multi::count;
use nom::sequence::terminated;
use nom::Parser;
use umol::error::{DataError, Result};

use super::atom::{atom_input, atom_like_input};
use super::bond::{bond_line, BondLine};

/// Result of molecule construction - either strict Molecule or general MoleculeLike
#[derive(Debug)]
pub enum MoleculeConstructionResult {
    Molecule(Molecule),
    MoleculeLike(MoleculeLike),
}

/// General molecule type that can handle all atom types
#[derive(Debug, Clone)]
pub struct MoleculeLike {
    pub atoms: Vec<AtomLike>,
    pub bonds: Vec<Bond>,
    pub conformers: Vec<Conformer>,
    pub properties: std::collections::HashMap<String, String>,
}

impl MoleculeLike {
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            conformers: Vec::new(),
            properties: std::collections::HashMap::new(),
        }
    }

    /// Try to convert to strict Molecule if all atoms are standard
    pub fn try_into_molecule(self) -> Option<Molecule> {
        let mut atoms = Vec::new();
        
        for atom_like in self.atoms {
            match atom_like.symbol {
                AtomSymbol::Element(element) => {
                    atoms.push(Atom {
                        element,
                        charge: atom_like.charge,
                        isotope_mass: atom_like.isotope_mass,
                        stereo_parity: atom_like.stereo_parity,
                        hydrogen_count: atom_like.hydrogen_count,
                        valence: atom_like.valence,
                        atom_map_num: atom_like.atom_map_num,
                        radical: None,
                        properties: std::collections::HashMap::new(),
                    });
                }
                AtomSymbol::NamedIsotope(isotope) => {
                    atoms.push(Atom {
                        element: isotope.element(),
                        charge: atom_like.charge,
                        isotope_mass: Some(isotope.mass_number()),
                        stereo_parity: atom_like.stereo_parity,
                        hydrogen_count: atom_like.hydrogen_count,
                        valence: atom_like.valence,
                        atom_map_num: atom_like.atom_map_num,
                        radical: None,
                        properties: std::collections::HashMap::new(),
                    });
                }
                _ => return None, // Non-standard atom found
            }
        }

        // Create molecule with converted atoms
        let mut molecule = Molecule::new();
        
        // Add atoms and collect positions for conformer
        let mut positions: Vec<Point3D> = Vec::new();
        for atom in atoms {
            molecule.add_atom(atom);
        }
        
        // Add bonds
        for bond in self.bonds {
            // Note: This assumes bond indices are correct
            // In a real implementation, we'd need to handle index mapping
        }
        
        // Add conformers
        for conformer in self.conformers {
            let _ = molecule.add_conformer(conformer);
        }
        
        // Add properties
        for (key, value) in self.properties {
            molecule.set_property(key, value);
        }

        Some(molecule)
    }
}

/// Parse molecule using strict parsing (fails on non-standard atoms)
pub(crate) fn molecule_input<'a>(
    atom_count: usize,
    bond_count: usize,
) -> impl Parser<&'a [u8], Output = Result<Molecule>, Error = error::Error<&'a [u8]>> {
    move |input| {
        // Parse atoms
        let (input, atom_position_pairs): (_, Vec<(Atom, Point3D)>) = 
            count(terminated(atom_input(), nom::character::complete::line_ending), atom_count)
                .parse(input)?;
        
        // Separate atoms and positions for conformer
        let mut atoms = Vec::new();
        let mut positions: Vec<Point3D> = Vec::new();
        
        for (atom, position) in atom_position_pairs {
            atoms.push(atom);
            positions.push(position);
        }
        
        // Parse bonds
        let (input, bond_lines): (_, Vec<BondLine>) = 
            count(terminated(bond_line(), nom::character::complete::line_ending), bond_count)
                .parse(input)?;
        
        // Convert BondLine to Bond (simplified conversion)
        let bonds: Vec<Bond> = bond_lines.into_iter().map(|_bond_line| {
            // For now, create a default bond - in practice would convert properly
            Bond::new(crate::bond::BondType::Single)
        }).collect();
        
        // Construct molecule
        let mut molecule = Molecule::new();
        
        // Add atoms
        for atom in atoms {
            molecule.add_atom(atom);
        }
        
        // Add bonds (simplified - assumes correct indexing)
        for bond in bonds {
            // In real implementation, would need proper index handling
        }
        
        // Create conformer from positions
        let conformer = Conformer { positions };
        let _ = molecule.add_conformer(conformer);
        
        Ok((input, Ok(molecule)))
    }
}

/// Parse molecule using permissive parsing (handles all atom types)
pub(crate) fn molecule_like_input<'a>(
    atom_count: usize,
    bond_count: usize,
) -> impl Parser<&'a [u8], Output = MoleculeLike, Error = error::Error<&'a [u8]>> {
    move |input| {
        // Parse atoms with permissive parser
        let (input, atom_results): (_, Vec<(AtomLike, Point3D)>) = 
            count(terminated(atom_like_input(), nom::character::complete::line_ending), atom_count)
                .parse(input)?;
        
        let (atoms, positions): (Vec<_>, Vec<_>) = atom_results.into_iter().unzip();
        
        // Parse bonds
        let (input, bond_lines): (_, Vec<BondLine>) = 
            count(terminated(bond_line(), nom::character::complete::line_ending), bond_count)
                .parse(input)?;
        
        // Convert BondLine to Bond (simplified conversion)
        let bonds: Vec<Bond> = bond_lines.into_iter().map(|_bond_line| {
            // For now, create a default bond - in practice would convert properly
            Bond::new(crate::bond::BondType::Single)
        }).collect();
        
        // Construct molecule-like
        let mut molecule_like = MoleculeLike::new();
        molecule_like.atoms = atoms;
        molecule_like.bonds = bonds;
        
        // Create conformer from positions
        let conformer = Conformer { positions };
        molecule_like.conformers.push(conformer);
        
        Ok((input, molecule_like))
    }
}

/// High-level API: Parse molecule with automatic fallback
pub fn parse_molecule_smart(input: &[u8], atom_count: usize, bond_count: usize) -> MoleculeConstructionResult {
    // Try strict parsing first
    if let Ok((_, Ok(molecule))) = molecule_input(atom_count, bond_count).parse(input) {
        return MoleculeConstructionResult::Molecule(molecule);
    }
    
    // Fall back to permissive parsing
    if let Ok((_, molecule_like)) = molecule_like_input(atom_count, bond_count).parse(input) {
        return MoleculeConstructionResult::MoleculeLike(molecule_like);
    }
    
    // If both fail, return empty MoleculeLike (could be improved with better error handling)
    MoleculeConstructionResult::MoleculeLike(MoleculeLike::new())
} 