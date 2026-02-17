//! Atom type for GraphIR.

use umol_data::{Element, SpinMultiplicity};

/// Basic atom IR
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    element: Element,
    isotope_mass: Option<u32>,
    charge: i8,
    hydrogens: u8,
    lone_pairs: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    unpaired_electrons: u8,
    multiplicity: SpinMultiplicity,
    valence: u8,
    aromatic_valence: u8,
    multicenter_valence: u8,
}

impl Atom {
    pub fn element(&self) -> Element {
        self.element
    }

    pub fn isotope_mass(&self) -> Option<u32> {
        self.isotope_mass
    }
    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn hydrogens(&self) -> u8 {
        self.hydrogens
    }

    pub fn lone_pairs(&self) -> u8 {
        self.lone_pairs
    }

    pub fn donated_pairs(&self) -> u8 {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> u8 {
        self.accepted_pairs
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.multiplicity
    }
    pub fn valence(&self) -> u8 {
        self.valence
    }

    pub fn aromatic_valence(&self) -> u8 {
        self.aromatic_valence
    }

    pub fn multicenter_valence(&self) -> u8 {
        self.multicenter_valence
    }
}
