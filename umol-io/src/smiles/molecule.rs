use std::str::FromStr;

use super::config::SmilesIoConfig;
use super::error::{ParseError, SmilesRenderError};
use super::parser::parse_molecule;
use super::render;
use crate::table_ir::{Molecule, SourceFormat};

/// Semantic value of a molecular SMILES representation.
///
/// Owns a TableIR for parsing and rendering. Wrapping a table preserves its contents;
/// rendering checks the properties needed to produce SMILES. Original spelling is not retained.
#[derive(Clone, Debug, PartialEq)]
pub struct Smiles {
    table_ir: Molecule,
}

impl Smiles {
    /// Parse SMILES text with the OpenSMILES configuration.
    ///
    /// Directional evidence is normalized into explicit bond frames.
    ///
    /// # Errors
    ///
    /// Rejects invalid syntax, conflicting direction markers, and markers outside the
    /// supported local double-bond domain. Consistent incomplete evidence adds no frame.
    pub fn parse(input: &str) -> Result<Self, ParseError> {
        Self::parse_bytes(input.as_bytes())
    }

    /// Parse SMILES bytes with the OpenSMILES configuration.
    pub fn parse_bytes(input: &[u8]) -> Result<Self, ParseError> {
        Self::parse_bytes_with(input, &SmilesIoConfig::opensmiles())
    }

    /// Parse SMILES text with an explicit IO configuration.
    pub fn parse_with(input: &str, config: &SmilesIoConfig) -> Result<Self, ParseError> {
        Self::parse_bytes_with(input.as_bytes(), config)
    }

    /// Parse SMILES bytes with an explicit IO configuration.
    pub fn parse_bytes_with(input: &[u8], config: &SmilesIoConfig) -> Result<Self, ParseError> {
        let mut table_ir = parse_molecule(input, config)?;
        table_ir.source_format = SourceFormat::SMILES;
        Ok(Self::from_table_ir(table_ir))
    }

    /// Render with the OpenSMILES configuration.
    ///
    /// # Errors
    ///
    /// Fails if the value contains fields or stereo that OpenSMILES output cannot express.
    pub fn render(&self) -> Result<String, SmilesRenderError> {
        self.render_with(&SmilesIoConfig::opensmiles())
    }

    /// Render under an explicit IO configuration, recomputing traversal and stereo markers.
    ///
    /// # Errors
    ///
    /// Fails for unsupported fields or unrepresentable stereo. Parsed CX annotations can
    /// fail even under the configuration that accepted them; CX output is not implemented.
    ///
    /// # Semantic properties
    ///
    /// Output is deterministic for the table and configuration. It preserves supported
    /// molecular semantics, including fixed versus inferred H counts and stereo assertions,
    /// but need not preserve source spelling, ring labels, or redundant direction markers.
    pub fn render_with(&self, config: &SmilesIoConfig) -> Result<String, SmilesRenderError> {
        render::render(&self.table_ir, config)
    }

    /// Take ownership of a TableIR for SMILES rendering.
    ///
    /// The table is retained unchanged. Rendering checks the properties it requires.
    pub fn from_table_ir(table_ir: Molecule) -> Self {
        Self { table_ir }
    }

    /// Consume the SMILES value and return its neutral TableIR boundary value.
    pub fn into_table_ir(self) -> Molecule {
        self.table_ir
    }

    /// Borrow the neutral TableIR boundary value.
    pub fn as_table_ir(&self) -> &Molecule {
        &self.table_ir
    }
}

impl FromStr for Smiles {
    type Err = ParseError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        Self::parse(input)
    }
}

#[cfg(test)]
mod tests {
    use indexmap::IndexMap;
    use rstest::rstest;
    use umol_chem::element::Element;

    use super::*;
    use crate::table_ir::{
        Atom, Bond, BondConfiguration, BondOrder, BondRelation, Span, StereoBond, StereoLigand,
    };

    #[rstest]
    #[case::empty(Molecule::empty(), "")]
    #[case::atom(Molecule { atoms: vec![Atom::aliphatic_atom(Element::C)], ..Molecule::empty() }, "C")]
    #[case::bond(Molecule { atoms: vec![Atom::aliphatic_atom(Element::C); 2],
        bonds: vec![Bond::new(0,1,BondOrder::Double)], ..Molecule::empty() }, "C=C")]
    fn test_smiles_from_table_ir(#[case] table: Molecule, #[case] expected: &str) {
        let config = SmilesIoConfig::opensmiles();
        let smiles = Smiles::from_table_ir(table.clone());
        assert_eq!(smiles.as_table_ir(), &table);
        assert_eq!(smiles.render_with(&config), Ok(expected.to_owned()));
        assert_eq!(smiles.into_table_ir(), table);
    }

    #[rstest]
    #[case::carbon(
        "C",
        Molecule {
            atoms: vec![Atom::aliphatic_atom_with_span(Element::C, Span::bytes(0, 1))],
            bonds: Vec::new(),
            positions: None,
            multicenter_bonds: Vec::new(),
            configuration_scope: None,
            stereo_atoms: Vec::new(),
            stereo_bonds: Vec::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::SMILES,
        }
    )]
    #[case::wildcard(
        "*",
        Molecule {
            atoms: vec![Atom::wildcard_with_span(Span::bytes(0, 1))],
            bonds: Vec::new(),
            positions: None,
            multicenter_bonds: Vec::new(),
            configuration_scope: None,
            stereo_atoms: Vec::new(),
            stereo_bonds: Vec::new(),
            comments: Vec::new(),
            properties: IndexMap::new(),
            source_format: SourceFormat::SMILES,
        }
    )]
    #[case::empty(
        "",
        Molecule {
            source_format: SourceFormat::SMILES,
            ..Molecule::empty()
        }
    )]
    fn test_smiles_parse(#[case] input: &str, #[case] expected: Molecule) {
        let smiles = Smiles::parse(input).unwrap();
        assert_eq!(smiles.as_table_ir(), &expected);
    }

    #[rstest]
    #[case::leading_whitespace(" C", ParseError::LeadingWhitespace)]
    #[case::invalid_element("Q", ParseError::InvalidElement { pos: 0 })]
    fn test_smiles_parse_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(Smiles::parse(input), Err(expected));
    }

    #[rstest]
    #[case::biphenyl("c1ccccc1-c2ccccc2", "c1ccccc1-c1ccccc1")]
    #[case::fixed_h("[CH3][CH2][OH]", "[CH3][CH2][OH]")]
    #[case::inferred_h("CCO", "CCO")]
    #[case::explicit_h("[H]CO", "[H]CO")]
    #[case::partial("F/C=CF", "FC=CF")]
    #[case::redundant("F/C(/Cl)=C(/Br)I", "FC(/Cl)=C(Br)\\I")]
    fn test_smiles_render(#[case] input: &str, #[case] expected: &str) {
        let smiles = Smiles::parse(input).unwrap();
        let original = smiles.clone();
        assert_eq!(smiles.render(), Ok(expected.to_owned()));
        assert_eq!(smiles, original);
    }

    #[rstest]
    #[case::endpoint(vec![Bond::new(0,2,BondOrder::Single)], vec![], SmilesRenderError::AtomIndexOutOfBounds { atom: 2 })]
    #[case::frame_bond(vec![Bond::new(0,1,BondOrder::Double)], vec![StereoBond { bond: 9, configuration: BondConfiguration::Framed {
        references: [0,1], relation: BondRelation::SameSide } }], SmilesRenderError::BondIndexOutOfBounds { bond: 9 })]
    fn test_smiles_render_table_error(
        #[case] bonds: Vec<Bond>,
        #[case] stereo_bonds: Vec<StereoBond>,
        #[case] expected: SmilesRenderError,
    ) {
        let table = Molecule {
            atoms: vec![Atom::aliphatic_atom(Element::C); 2],
            bonds,
            stereo_bonds,
            ..Molecule::empty()
        };
        assert_eq!(Smiles::from_table_ir(table).render(), Err(expected));
    }

    #[rstest]
    #[case::missing(5, SmilesRenderError::AtomIndexOutOfBounds { atom: 5 })]
    #[case::duplicate(1, SmilesRenderError::DuplicateStereoAtom { atom: 1 })]
    fn test_smiles_render_stereo_site_error(
        #[case] atom: u32,
        #[case] expected: SmilesRenderError,
    ) {
        let mut molecule = Smiles::parse("F[C@](Cl)(Br)I").unwrap().into_table_ir();
        let mut frame = molecule.stereo_atoms[0].clone();
        frame.atom = atom;
        molecule.stereo_atoms.push(frame);
        assert_eq!(Smiles::from_table_ir(molecule).render(), Err(expected));
    }

    #[rstest]
    #[case::short(vec![StereoLigand::Atom(0), StereoLigand::Atom(2), StereoLigand::Atom(3)], SmilesRenderError::InvalidStereoAtom { atom: 1 })]
    #[case::duplicate(vec![StereoLigand::Atom(0), StereoLigand::Atom(2), StereoLigand::Atom(3), StereoLigand::Atom(3)], SmilesRenderError::InvalidStereoAtom { atom: 1 })]
    #[case::absent(vec![StereoLigand::Atom(0), StereoLigand::Atom(2), StereoLigand::Atom(3), StereoLigand::Atom(9)], SmilesRenderError::InvalidStereoAtom { atom: 1 })]
    fn test_smiles_render_stereo_frame_error(
        #[case] ligands: Vec<StereoLigand>,
        #[case] expected: SmilesRenderError,
    ) {
        let mut table = Smiles::parse("F[C@](Cl)(Br)I").unwrap().into_table_ir();
        table.stereo_atoms[0].ligands = ligands;
        assert_eq!(Smiles::from_table_ir(table).render(), Err(expected));
    }

    #[rstest]
    #[case::virtual_ligand(StereoLigand::ImplicitHydrogen)]
    fn test_smiles_render_stereo_atom_error(#[case] ligand: StereoLigand) {
        let mut table = Smiles::parse("F[C@](Cl)(Br)I").unwrap().into_table_ir();
        table.stereo_atoms[0].ligands[3] = ligand;
        let smiles = Smiles::from_table_ir(table);
        assert_eq!(
            smiles.render(),
            Err(SmilesRenderError::InvalidStereoAtom { atom: 1 })
        );
    }

    #[rstest]
    #[case::limit(102, 100)]
    fn test_smiles_render_ring_error(#[case] count: usize, #[case] label: usize) {
        let table = Molecule {
            atoms: vec![Atom::aliphatic_atom(Element::C); count],
            bonds: (1..count)
                .map(|atom| Bond::new(atom as u32 - 1, atom as u32, BondOrder::Single))
                .chain((2..count).map(|atom| Bond::new(0, atom as u32, BondOrder::Single)))
                .collect(),
            ..Molecule::empty()
        };
        assert_eq!(
            Smiles::from_table_ir(table).render(),
            Err(SmilesRenderError::RingLabel { label })
        );
    }

    #[rstest]
    #[case::dative("N->[Cu]", SmilesRenderError::UnsupportedBond { bond: 0, field: "donation" })]
    #[case::any("C~C", SmilesRenderError::UnsupportedBond { bond: 0, field: "order" })]
    #[case::aromatic("[te]", SmilesRenderError::UnsupportedAtom { atom: 0, field: "aromatic" })]
    fn test_smiles_render_with(#[case] input: &str, #[case] narrower: SmilesRenderError) {
        let config = SmilesIoConfig::lenient();
        let table = Smiles::parse_with(input, &config).unwrap().into_table_ir();
        let smiles = Smiles::from_table_ir(table.clone());
        assert_eq!(smiles.render_with(&config), Ok(input.to_owned()));
        assert_eq!(smiles.render(), Err(narrower.clone()));
        assert_eq!(Smiles::from_table_ir(table).render(), Err(narrower));
    }

    #[rstest]
    #[case::label("C |$name$|", SmilesRenderError::UnsupportedAtom { atom: 0, field: "label" })]
    #[case::either("FC=CF |ctu:1|", SmilesRenderError::UnsupportedBond { bond: 1, field: "Either" })]
    fn test_smiles_render_with_error(#[case] input: &str, #[case] reason: SmilesRenderError) {
        let config = SmilesIoConfig::chemaxon();
        let smiles = Smiles::parse_with(input, &config).unwrap();
        assert_eq!(smiles.render_with(&config), Err(reason.clone()));
        assert_eq!(
            Smiles::from_table_ir(smiles.into_table_ir()).render_with(&config),
            Err(reason)
        );
    }
}
