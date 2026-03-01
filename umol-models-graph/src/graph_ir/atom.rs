//! Atom type and builder for GraphIR.

use smallvec::SmallVec;
use umol_data::{Element, SpinMultiplicity};

use super::error::ResolutionError;
use super::valence::AtomTypeSpec;
use crate::table_ir::atom::{Atom as TableAtom, Chirality};

/// Resolved atom in GraphIR. All fields are definite.
/// Created via `AtomBuilder::build()` after resolution phases complete.
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    element: Element,
    isotope_mass: Option<u32>,
    charge: i8,
    hydrogens: u8,
    lone_pairs: u8,
    unpaired_electrons: u8,
    multiplicity: SpinMultiplicity,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
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

    pub fn unpaired_electrons(&self) -> u8 {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.multiplicity
    }

    pub fn valence(&self) -> u8 {
        self.valence
    }

    pub fn donated_pairs(&self) -> u8 {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> u8 {
        self.accepted_pairs
    }

    pub fn aromatic_valence(&self) -> u8 {
        self.aromatic_valence
    }

    pub fn multicenter_valence(&self) -> u8 {
        self.multicenter_valence
    }

    pub fn to_builder(&self) -> AtomBuilder {
        let spec = AtomTypeSpec::new(
            self.element,
            self.charge,
            self.hydrogens,
            self.lone_pairs,
            self.unpaired_electrons,
            self.multiplicity,
            self.valence,
            self.donated_pairs,
            self.accepted_pairs,
            self.aromatic_valence,
            self.multicenter_valence,
        )
        .expect("resolved Atom fields always form a valid AtomTypeSpec");
        AtomBuilder {
            element: self.element,
            isotope_mass: self.isotope_mass,
            charge: Some(self.charge),
            hydrogens: Some(self.hydrogens),
            lone_pairs: Some(self.lone_pairs),
            unpaired_electrons: Some(self.unpaired_electrons),
            multiplicity: Some(self.multiplicity),
            aromatic_hint: None,
            chirality_hint: None,
            candidates: SmallVec::from_elem(spec, 1),
        }
    }
}

/// Builder for constructing `Atom` values. Used as graph node weights during
/// resolution phases and for manual molecule construction. Fields progress
/// from `None` (unknown) to `Some` as resolution phases fill them in.
#[derive(Debug, Clone)]
pub struct AtomBuilder {
    element: Element,
    isotope_mass: Option<u32>,
    charge: Option<i8>,
    hydrogens: Option<u8>,
    lone_pairs: Option<u8>,
    unpaired_electrons: Option<u8>,
    multiplicity: Option<SpinMultiplicity>,
    aromatic_hint: Option<bool>,
    chirality_hint: Option<Chirality>,
    candidates: SmallVec<[AtomTypeSpec; 4]>,
}

impl AtomBuilder {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            hydrogens: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic_hint: None,
            chirality_hint: None,
            candidates: SmallVec::new(),
        }
    }

    pub fn from_table_atom(atom: &TableAtom) -> Self {
        Self {
            element: atom.element,
            isotope_mass: atom.isotope_mass,
            charge: atom.charge,
            hydrogens: atom.hydrogens,
            lone_pairs: atom.lone_pairs,
            unpaired_electrons: atom.unpaired_electrons.map(|u| u.count),
            multiplicity: atom.unpaired_electrons.and_then(|u| u.multiplicity),
            aromatic_hint: atom.aromatic,
            chirality_hint: atom.chirality,
            candidates: SmallVec::new(),
        }
    }

    pub fn element(&self) -> Element {
        self.element
    }

    pub fn isotope_mass(&self) -> Option<u32> {
        self.isotope_mass
    }

    pub fn charge(&self) -> Option<i8> {
        self.charge
    }

    pub fn hydrogens(&self) -> Option<u8> {
        self.hydrogens
    }

    pub fn lone_pairs(&self) -> Option<u8> {
        self.lone_pairs
    }

    pub fn unpaired_electrons(&self) -> Option<u8> {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> Option<SpinMultiplicity> {
        self.multiplicity
    }

    pub fn aromatic_hint(&self) -> Option<bool> {
        self.aromatic_hint
    }

    pub fn chirality_hint(&self) -> Option<&Chirality> {
        self.chirality_hint.as_ref()
    }

    pub fn candidates(&self) -> &[AtomTypeSpec] {
        &self.candidates
    }

    pub fn set_isotope_mass(&mut self, isotope_mass: u32) -> &mut Self {
        self.isotope_mass = Some(isotope_mass);
        self
    }

    pub fn set_charge(&mut self, charge: i8) -> &mut Self {
        self.charge = Some(charge);
        self
    }

    pub fn set_hydrogens(&mut self, hydrogens: u8) -> &mut Self {
        self.hydrogens = Some(hydrogens);
        self
    }

    pub fn set_lone_pairs(&mut self, lone_pairs: u8) -> &mut Self {
        self.lone_pairs = Some(lone_pairs);
        self
    }

    // TODO: Check consistency with multiplicity
    pub fn set_unpaired_electrons(&mut self, unpaired_electrons: u8) -> &mut Self {
        self.unpaired_electrons = Some(unpaired_electrons);
        self
    }

    // TODO: Check consistency with unpaired electrons
    pub fn set_multiplicity(&mut self, multiplicity: SpinMultiplicity) -> &mut Self {
        self.multiplicity = Some(multiplicity);
        self
    }

    pub fn set_aromatic_hint(&mut self, aromatic: bool) -> &mut Self {
        self.aromatic_hint = Some(aromatic);
        self
    }

    pub fn set_chirality_hint(&mut self, chirality: Chirality) -> &mut Self {
        self.chirality_hint = Some(chirality);
        self
    }

    pub fn set_candidates(&mut self, candidates: SmallVec<[AtomTypeSpec; 4]>) -> &mut Self {
        self.candidates = candidates;
        self
    }

    pub fn add_candidate(&mut self, candidate: AtomTypeSpec) -> &mut Self {
        if !self.candidates.contains(&candidate) {
            self.candidates.push(candidate);
        }
        self
    }

    /// Build the final `Atom` from the builder state.
    ///
    /// Requires exactly one valence candidate remaining (resolution phases
    /// must have narrowed the set). All fields on that candidate become the
    /// definite atom properties.
    pub fn build(&self) -> Result<Atom, ResolutionError> {
        let candidate = match self.candidates.len() {
            0 => {
                return Err(ResolutionError::ValenceNoMatch(format!(
                    "no valence match for {:?}",
                    self.element
                )))
            }
            1 => &self.candidates[0],
            n => {
                return Err(ResolutionError::ValenceAmbiguous(format!(
                    "{} valence matches for {:?}",
                    n, self.element
                )))
            }
        };

        Ok(Atom {
            element: self.element,
            isotope_mass: self.isotope_mass,
            charge: candidate.charge(),
            hydrogens: candidate.hydrogens(),
            lone_pairs: candidate.lone_pairs(),
            unpaired_electrons: candidate.unpaired_electrons(),
            multiplicity: candidate.multiplicity(),
            valence: candidate.valence(),
            donated_pairs: candidate.donated_pairs(),
            accepted_pairs: candidate.accepted_pairs(),
            aromatic_valence: candidate.aromatic_valence(),
            multicenter_valence: candidate.multicenter_valence(),
        })
    }
}
