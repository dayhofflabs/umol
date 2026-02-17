//! Atom types for GraphIR.

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
    deloc_electrons: u8,
    multiplicity: SpinMultiplicity,
    class: Option<u32>,
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

    pub fn deloc_electrons(&self) -> u8 {
        self.deloc_electrons
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.multiplicity
    }

    pub fn class(&self) -> Option<u32> {
        self.class
    }
}

impl From<Atom> for AtomBuilder {
    fn from(atom: Atom) -> Self {
        Self {
            element: atom.element,
            isotope_mass: atom.isotope_mass,
            charge: Some(atom.charge),
            hydrogens: Some(atom.hydrogens),
            lone_pairs: Some(atom.lone_pairs),
            donated_pairs: Some(atom.donated_pairs),
            accepted_pairs: Some(atom.accepted_pairs),
            unpaired_electrons: Some(atom.unpaired_electrons),
            deloc_electrons: Some(atom.deloc_electrons),
            multiplicity: Some(atom.multiplicity),
            class: atom.class,
        }
    }
}

/// Builder type for creating and mutating `Atom` types including strict typing.
#[derive(Debug)]
pub struct AtomBuilder {
    element: Element,
    isotope_mass: Option<u32>,
    charge: Option<i8>,
    hydrogens: Option<u8>,
    lone_pairs: Option<u8>,
    donated_pairs: Option<u8>,
    accepted_pairs: Option<u8>,
    unpaired_electrons: Option<u8>,
    deloc_electrons: Option<u8>,
    multiplicity: Option<SpinMultiplicity>,
    class: Option<u32>,
}

impl AtomBuilder {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            hydrogens: None,
            lone_pairs: None,
            donated_pairs: None,
            accepted_pairs: None,
            unpaired_electrons: None,
            deloc_electrons: None,
            multiplicity: None,
            class: None,
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

    pub fn donated_pairs(&self) -> Option<u8> {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> Option<u8> {
        self.accepted_pairs
    }

    pub fn unpaired_electrons(&self) -> Option<u8> {
        self.unpaired_electrons
    }

    pub fn deloc_electrons(&self) -> Option<u8> {
        self.deloc_electrons
    }

    pub fn multiplicity(&self) -> Option<SpinMultiplicity> {
        self.multiplicity
    }

    pub fn class(&self) -> Option<u32> {
        self.class
    }

    pub fn set_element(&mut self, element: Element) -> &mut Self {
        self.element = element;
        self
    }

    pub fn set_isotope_mass(&mut self, isotope_mass: u32) -> &mut Self {
        self.isotope_mass = Some(isotope_mass);
        self
    }

    pub fn set_charge(&mut self, charge: i8) -> &mut Self {
        self.charge = Some(charge);
        self
    }

    pub fn set_lone_pairs(&mut self, count: u8) -> &mut Self {
        self.lone_pairs = Some(count);
        self
    }

    pub fn set_donated_pairs(&mut self, count: u8) -> &mut Self {
        self.donated_pairs = Some(count);
        self
    }

    pub fn set_accepted_pairs(&mut self, count: u8) -> &mut Self {
        self.accepted_pairs = Some(count);
        self
    }

    pub fn set_unpaired_electrons(&mut self, count: u8) -> &mut Self {
        self.unpaired_electrons = Some(count);
        self
    }

    pub fn set_deloc_electrons(&mut self, count: u8) -> &mut Self {
        self.deloc_electrons = Some(count);
        self
    }

    pub fn set_multiplicity(&mut self, multiplicity: SpinMultiplicity) -> &mut Self {
        self.multiplicity = Some(multiplicity);
        self
    }

    pub fn set_class(&mut self, class: u32) -> &mut Self {
        self.class = Some(class);
        self
    }

    pub fn update_element(&mut self, f: impl FnOnce(Element) -> Element) -> &mut Self {
        self.element = f(self.element);
        self
    }

    pub fn update_isotope_mass(&mut self, f: impl FnOnce(u32) -> u32) -> &mut Self {
        self.isotope_mass = self.isotope_mass.map(f);
        self
    }

    pub fn update_charge(&mut self, f: impl FnOnce(i8) -> i8) -> &mut Self {
        self.charge = Some(f(self.charge.unwrap_or(0)));
        self
    }

    pub fn update_hydrogens(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.hydrogens = Some(f(self.hydrogens.unwrap_or(0)));
        self
    }

    pub fn update_lone_pairs(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.lone_pairs = Some(f(self.lone_pairs.unwrap_or(0)));
        self
    }

    pub fn update_donated_pairs(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.donated_pairs = Some(f(self.donated_pairs.unwrap_or(0)));
        self
    }

    pub fn update_accepted_pairs(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.accepted_pairs = Some(f(self.accepted_pairs.unwrap_or(0)));
        self
    }

    pub fn update_unpaired_electrons(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.unpaired_electrons = Some(f(self.unpaired_electrons.unwrap_or(0)));
        self
    }

    pub fn update_deloc_electrons(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.deloc_electrons = Some(f(self.deloc_electrons.unwrap_or(0)));
        self
    }

    pub fn update_multiplicity(
        &mut self,
        f: impl FnOnce(SpinMultiplicity) -> SpinMultiplicity,
    ) -> &mut Self {
        self.multiplicity = self.multiplicity.map(f);
        self
    }

    pub fn update_class(&mut self, f: impl FnOnce(Option<u32>) -> Option<u32>) -> &mut Self {
        self.class = f(self.class);
        self
    }

        // pub fn build(self) -> Result<Atom, ResolutionError> {
        //     self.build_with(&DEFAULT_ATOM_VALIDATOR, &DEFAULT_ATOM_MATCHER)
        // }

    //     pub fn build_with(
    //         self,
    //         validator: &AtomValidator,
    //         matcher: &AtomMatcher,
    //     ) -> Result<Atom, ResolutionError> {
    //         let atom_specs = matcher.find(&self)?;
    //         if atom_specs.is_empty() {
    //             return Err(ResolutionError::InvalidAtomSpec(format!("{:?}", self)));
    //         } else if atom_specs.len() > 1 {
    //             return Err(ResolutionError::InvalidAtomSpec(format!(
    //                 "{:?}: {}",
    //                 self,
    //                 atom_specs
    //                     .iter()
    //                     .map(|s| s.to_string())
    //                     .collect::<Vec<String>>()
    //                     .join(", ")
    //             )));
    //         }
    //         let atom_spec = atom_specs.first().unwrap();
    //         let atom = Atom {
    //             element: atom_spec.element(),
    //             charge: self.charge.unwrap_or(atom_spec.charge()),
    //             valence: self.valence.unwrap_or_else(|| atom_spec.valence()),
    //             multiplicity: self
    //                 .multiplicity
    //                 .unwrap_or_else(|| atom_spec.multiplicity()),
    //             lone_pairs: self.lone_pairs.unwrap_or_else(|| atom_spec.lone_pairs()),
    //             donated_pairs: self
    //                 .donated_pairs
    //                 .unwrap_or_else(|| atom_spec.donated_pairs()),
    //             accepted_pairs: self
    //                 .accepted_pairs
    //                 .unwrap_or_else(|| atom_spec.accepted_pairs()),
    //             unpaired_electrons: self
    //                 .unpaired_electrons
    //                 .unwrap_or_else(|| atom_spec.unpaired_electrons()),
    //             implicit_hydrogens: self
    //                 .implicit_hydrogens
    //                 .unwrap_or_else(|| atom_spec.implicit_hydrogens()),
    //             position: self.position,
    //             isotope_mass: self.isotope_mass,
    //             aromatic: self.aromatic,
    //             chirality: self.chirality,
    //             class: self.class,
    //             span: self.span,
    //         };
    //         validator.validate(&atom)?;
    //         Ok(atom)
    //     }
}
