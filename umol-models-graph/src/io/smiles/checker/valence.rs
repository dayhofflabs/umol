use std::collections::HashMap;

use crate::io::smiles::diagnostics::{Category, Code, Diagnostic, DiagnosticList, Severity, Span};
use crate::io::ir::{AtomSymbol, BondOrder, BondSymbol, Molecule};
use umol_data::{e, Element};

#[derive(Clone, Copy)]
pub enum ValencePolicy {
    Off,
    Warn,
    Error,
}

#[derive(Clone, Copy)]
pub struct ValenceConfig {
    pub enabled: bool,
    pub overflow_policy: ValencePolicy,
    pub check_bracket: bool,
    pub infer_bracket_implicit: bool,
    pub aromatic_as_one: bool,
    pub patterns_enabled: bool,
    pub no_match_policy: ValencePolicy,
    pub ambiguous_match_policy: ValencePolicy,
}

impl Default for ValenceConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            overflow_policy: ValencePolicy::Warn,
            check_bracket: true,
            infer_bracket_implicit: false,
            aromatic_as_one: true,
            patterns_enabled: true,
            no_match_policy: ValencePolicy::Off,
            ambiguous_match_policy: ValencePolicy::Off,
        }
    }
}

#[derive(Default)]
pub struct ValenceModel {
    states: HashMap<Element, Vec<u8>>,
    pub patterns: ValencePatternTable,
}

impl ValenceModel {
    pub fn simple_organic() -> Self {
        Self { states: HashMap::new(), patterns: ValencePatternTable::organic_strict() }
    }
    pub fn states_for(&self, e: Element) -> Option<&[u8]> { self.states.get(&e).map(|v| &v[..]) }
    pub fn set_states(&mut self, e: Element, states: Vec<u8>) { self.states.insert(e, states); }
}

#[derive(Clone, Copy, Debug)]
pub struct ValencePattern {
    pub element: Option<Element>,
    pub bond_sum: Option<u8>,
    pub charge: Option<i8>,
    pub implicit_h: Option<u8>,
    pub unpaired: Option<u8>,
}

#[derive(Default)]
pub struct ValencePatternTable { pub patterns: Vec<ValencePattern> }

impl ValencePatternTable {
    #[rustfmt::skip]
    pub fn organic_strict() -> Self {
        Self { patterns: vec![
            // H
            ValencePattern { element: Some(e!(H)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // HR
            ValencePattern { element: Some(e!(H)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // H^
            ValencePattern { element: Some(e!(H)), bond_sum: Some(0), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // H+
            ValencePattern { element: Some(e!(H)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // H-
        
            // B
            ValencePattern { element: Some(e!(B)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // BR3, B(OH)3, B2O3
            ValencePattern { element: Some(e!(B)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HBR2
            ValencePattern { element: Some(e!(B)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(2), unpaired: Some(0) }, // H2BR
            ValencePattern { element: Some(e!(B)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(3), unpaired: Some(0) }, // H3B
            ValencePattern { element: Some(e!(B)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // BR2^
            ValencePattern { element: Some(e!(B)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(1), unpaired: Some(1) }, // HBR^
            ValencePattern { element: Some(e!(B)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(2), unpaired: Some(1) }, // H2B^
            ValencePattern { element: Some(e!(B)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(3) }, // B^^^
            ValencePattern { element: Some(e!(B)), bond_sum: Some(4), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // BR4-
            ValencePattern { element: Some(e!(B)), bond_sum: Some(3), charge: Some(-1), implicit_h: Some(1), unpaired: Some(0) }, // HBR3-
            ValencePattern { element: Some(e!(B)), bond_sum: Some(2), charge: Some(-1), implicit_h: Some(2), unpaired: Some(0) }, // H2BR2-
            ValencePattern { element: Some(e!(B)), bond_sum: Some(1), charge: Some(-1), implicit_h: Some(3), unpaired: Some(0) }, // H3BR-
            ValencePattern { element: Some(e!(B)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(4), unpaired: Some(0) }, // H4B-
        
            // C
            ValencePattern { element: Some(e!(C)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // CR4, CO2, HOC(=O)R
            ValencePattern { element: Some(e!(C)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HCR3, HC(=O)R, HC(=O)OH
            ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(2), unpaired: Some(0) }, // H2CR2, H2C=O
            ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(3), unpaired: Some(0) }, // H3CR
            ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(4), unpaired: Some(0) }, // H4C
            ValencePattern { element: Some(e!(C)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // CR3^
            ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(1), unpaired: Some(1) }, // HCR2^
            ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(2), unpaired: Some(1) }, // H2CR^
            ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(3), unpaired: Some(1) }, // H3C^
            ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(2) }, // CR2^^
            ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(1), unpaired: Some(2) }, // HCR^^
            ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(2), unpaired: Some(2) }, // H2C^^
            ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(2) }, // C^^, CO
            ValencePattern { element: Some(e!(C)), bond_sum: Some(3), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // CR3+
            ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(1), implicit_h: Some(1), unpaired: Some(0) }, // HCR2+
            ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(1), implicit_h: Some(2), unpaired: Some(0) }, // H2CR+
            ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(1), implicit_h: Some(3), unpaired: Some(0) }, // H3C+
            ValencePattern { element: Some(e!(C)), bond_sum: Some(3), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // CR3-
            ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(-1), implicit_h: Some(1), unpaired: Some(0) }, // HCR2-
            ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(-1), implicit_h: Some(2), unpaired: Some(0) }, // H2CR-
            ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(3), unpaired: Some(0) }, // H3C-
        
            // N
            ValencePattern { element: Some(e!(N)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // NR3, HONO
            ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HNR2
            ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(2), unpaired: Some(0) }, // H2NR
            ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(3), unpaired: Some(0) }, // H3N
            ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // NR2^, NO
            ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(1), unpaired: Some(1) }, // HNR^
            ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(2), unpaired: Some(1) }, // H2N^
            ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(1), unpaired: Some(2) }, // NR^^
            ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(2) }, // HN^^
            ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(3) }, // N^^^
            ValencePattern { element: Some(e!(N)), bond_sum: Some(4), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // NR4+
            ValencePattern { element: Some(e!(N)), bond_sum: Some(3), charge: Some(1), implicit_h: Some(1), unpaired: Some(0) }, // HNR3+
            ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(1), implicit_h: Some(2), unpaired: Some(0) }, // H2NR2+
            ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(1), implicit_h: Some(3), unpaired: Some(0) }, // H3NR+
            ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(1), implicit_h: Some(4), unpaired: Some(0) }, // H4N+
            ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // NO+
            ValencePattern { element: Some(e!(N)), bond_sum: Some(2), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // NR2-
            ValencePattern { element: Some(e!(N)), bond_sum: Some(1), charge: Some(-1), implicit_h: Some(1), unpaired: Some(0) }, // HNR-
            ValencePattern { element: Some(e!(N)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(2), unpaired: Some(0) }, // H2N-
            ValencePattern { element: Some(e!(N)), bond_sum: Some(5), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // NRO2, HONO2
            ValencePattern { element: Some(e!(N)), bond_sum: Some(4), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // NO2+
            ValencePattern { element: Some(e!(N)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // NO2
        
            // O
            ValencePattern { element: Some(e!(O)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // OR2
            ValencePattern { element: Some(e!(O)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HOR
            ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(2), unpaired: Some(0) }, // H2O
            ValencePattern { element: Some(e!(O)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // OR^
            ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // HO^
            ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(2) }, // O^^
            ValencePattern { element: Some(e!(O)), bond_sum: Some(3), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // OR3+
            ValencePattern { element: Some(e!(O)), bond_sum: Some(2), charge: Some(1), implicit_h: Some(1), unpaired: Some(0) }, // HOR2+
            ValencePattern { element: Some(e!(O)), bond_sum: Some(1), charge: Some(1), implicit_h: Some(2), unpaired: Some(0) }, // H2OR+
            ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(1), implicit_h: Some(3), unpaired: Some(0) }, // H3O+
            ValencePattern { element: Some(e!(O)), bond_sum: Some(1), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // OR-
            ValencePattern { element: Some(e!(O)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(1), unpaired: Some(0) }, // HO-
        
            // F
            ValencePattern { element: Some(e!(F)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // FR
            ValencePattern { element: Some(e!(F)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HF
            ValencePattern { element: Some(e!(F)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // F^
            ValencePattern { element: Some(e!(F)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // F-
        
            // P
            ValencePattern { element: Some(e!(P)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // PR3, P(OH)3, P4O6
            ValencePattern { element: Some(e!(P)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HPR2
            ValencePattern { element: Some(e!(P)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(2), unpaired: Some(0) }, // H2PR
            ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(3), unpaired: Some(0) }, // H3P
            ValencePattern { element: Some(e!(P)), bond_sum: Some(4), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // PR4+
            ValencePattern { element: Some(e!(P)), bond_sum: Some(3), charge: Some(1), implicit_h: Some(1), unpaired: Some(0) }, // HPR3+
            ValencePattern { element: Some(e!(P)), bond_sum: Some(2), charge: Some(1), implicit_h: Some(2), unpaired: Some(0) }, // H2PR2+
            ValencePattern { element: Some(e!(P)), bond_sum: Some(1), charge: Some(1), implicit_h: Some(3), unpaired: Some(0) }, // H3PR+
            ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(1), implicit_h: Some(4), unpaired: Some(0) }, // H4P+
            ValencePattern { element: Some(e!(P)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // PR2^
            ValencePattern { element: Some(e!(P)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(1), unpaired: Some(1) }, // HPR^
            ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(2), unpaired: Some(1) }, // H2P^
            ValencePattern { element: Some(e!(P)), bond_sum: Some(2), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // PR2-
            ValencePattern { element: Some(e!(P)), bond_sum: Some(1), charge: Some(-1), implicit_h: Some(1), unpaired: Some(0) }, // HPR-
            ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(2), unpaired: Some(0) }, // H2P-
            ValencePattern { element: Some(e!(P)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // P^^^
            ValencePattern { element: Some(e!(P)), bond_sum: Some(5), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // PR3(=O), PO(OH)3, P4O10
        
            // S
            ValencePattern { element: Some(e!(S)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // SR2
            ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HSR
            ValencePattern { element: Some(e!(S)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(2), unpaired: Some(0) }, // H2S
            ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // SR^
            ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(1), unpaired: Some(1) }, // SH^
            ValencePattern { element: Some(e!(S)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(2) }, // S^^
            ValencePattern { element: Some(e!(S)), bond_sum: Some(3), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // SR3+
            ValencePattern { element: Some(e!(S)), bond_sum: Some(2), charge: Some(1), implicit_h: Some(1), unpaired: Some(0) }, // HSR2+
            ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(1), implicit_h: Some(2), unpaired: Some(0) }, // H2SR+
            ValencePattern { element: Some(e!(S)), bond_sum: Some(0), charge: Some(1), implicit_h: Some(3), unpaired: Some(0) }, // H3S+
            ValencePattern { element: Some(e!(S)), bond_sum: Some(1), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // SR-
            ValencePattern { element: Some(e!(S)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(1), unpaired: Some(0) }, // HS-
            ValencePattern { element: Some(e!(S)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // S(=O)R2, SO2, SO(OH)2
            ValencePattern { element: Some(e!(S)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HS(=O)R
            ValencePattern { element: Some(e!(S)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(2), unpaired: Some(0) }, // H2S=O
            ValencePattern { element: Some(e!(S)), bond_sum: Some(6), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // S(=O)(=O)R2, SO3, SO(OH)2
            ValencePattern { element: Some(e!(S)), bond_sum: Some(5), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HS(=O)(=O)R
            ValencePattern { element: Some(e!(S)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(2), unpaired: Some(0) }, // H2S(=O)=O
        
            // Cl
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // ClR
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HCl, HOCl, Cl2O
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // Cl^
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // Cl-
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // ClO
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(2), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // ClO+
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(2), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // ClO-
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // Cl(=O)OH, RClO
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HClO
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // ClO2
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(4), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // ClO2+
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(4), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // ClO2-
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(5), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // Cl(=O)(=O)OH, RClO2
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HCl(=O)(=O)
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(6), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // ClO3
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(6), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // ClO3+
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(6), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // ClO3-
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(7), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // Cl(=O)(=O)(=O)OH, RClO3, ClO4-
            ValencePattern { element: Some(e!(Cl)), bond_sum: Some(6), charge: Some(0), implicit_h: Some(1), unpaired: Some(1) }, // HClO3
        
            // Br
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // BrR
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HBr, HOBr
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // Br^
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // Br-
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // BrO
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(2), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // BrO+
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(2), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // BrO-
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // Br(=O)OH, RBrO
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HBrO
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // BrO2
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(4), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // BrO2+
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(4), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // BrO2-
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(5), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // Br(=O)(=O)OH, RBrO2
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HBrO2
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(6), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // BrO3
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(6), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // BrO3+
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(6), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // BrO3-
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(7), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // Br(=O)(=O)(=O)OH, RBrO3, BrO4-
            ValencePattern { element: Some(e!(Br)), bond_sum: Some(6), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HBrO3
        
            // I
            ValencePattern { element: Some(e!(I)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // IR
            ValencePattern { element: Some(e!(I)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HI, HOI
            ValencePattern { element: Some(e!(I)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // I^
            ValencePattern { element: Some(e!(I)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // I-
            ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // IO
            ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // IO+
            ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // IO-
            ValencePattern { element: Some(e!(I)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // I(=O)OH, RIO
            ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HIO
            ValencePattern { element: Some(e!(I)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // IO2
            ValencePattern { element: Some(e!(I)), bond_sum: Some(4), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // IO2+
            ValencePattern { element: Some(e!(I)), bond_sum: Some(4), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // IO2-
            ValencePattern { element: Some(e!(I)), bond_sum: Some(5), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // I(=O)(=O)OH, RIO2
            ValencePattern { element: Some(e!(I)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HIO
            ValencePattern { element: Some(e!(I)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HIO2
            ValencePattern { element: Some(e!(I)), bond_sum: Some(6), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) }, // IO3
            ValencePattern { element: Some(e!(I)), bond_sum: Some(6), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) }, // IO3+
            ValencePattern { element: Some(e!(I)), bond_sum: Some(6), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) }, // IO3-
            ValencePattern { element: Some(e!(I)), bond_sum: Some(7), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) }, // I(=O)(=O)(=O)OH, RIO3, IO4-
            ValencePattern { element: Some(e!(I)), bond_sum: Some(6), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) }, // HIO3
        ]}
    }
}

fn match_score(p: &ValencePattern, element: Element, bond_sum: u8, charge: i8, unpaired: Option<u8>) -> Option<usize> {
    if let Some(e) = p.element { if e != element { return None; } }
    if let Some(bs) = p.bond_sum { if bs != bond_sum { return None; } }
    if let Some(ch) = p.charge { if ch != charge { return None; } }
    match (p.unpaired, unpaired) { (Some(pu), Some(u)) if pu != u => return None, _ => {} }
    let mut score = 0usize;
    if p.element.is_some() { score += 1; }
    if p.bond_sum.is_some() { score += 1; }
    if p.charge.is_some() { score += 1; }
    if p.unpaired.is_some() && unpaired.is_some() { score += 1; }
    if p.implicit_h.is_some() { score += 1; }
    Some(score)
}

#[derive(Default)]
struct PatternDecision { implicit_h: Option<u8>, matches: usize }

fn resolve_valence_pattern(
    element: Element,
    bond_sum: u8,
    charge: i8,
    unpaired: Option<u8>,
    tbl: &ValencePatternTable,
) -> PatternDecision {
    let mut best_score = 0usize; let mut best_idx: Option<usize> = None; let mut match_count = 0usize;
    for (idx, p) in tbl.patterns.iter().enumerate() {
        if let Some(score) = match_score(p, element, bond_sum, charge, unpaired) {
            match_count += 1; if score > best_score { best_score = score; best_idx = Some(idx); }
        }
    }
    if let Some(i) = best_idx { let p = tbl.patterns[i]; PatternDecision { implicit_h: p.implicit_h, matches: match_count } } else { PatternDecision::default() }
}

pub struct ValenceArtifacts {
    pub atoms_checked: usize,
    pub overflow_count: usize,
    pub bracket_mismatch_count: usize,
}

pub fn check_valence(
    mol: &Molecule,
    report: &mut DiagnosticList,
    model: &ValenceModel,
    cfg: &ValenceConfig,
) -> ValenceArtifacts {
    if !cfg.enabled {
        return ValenceArtifacts { atoms_checked: 0, overflow_count: 0, bracket_mismatch_count: 0 };
    }
    let atom_len = mol.atoms.len();
    let mut order_sum: Vec<u8> = vec![0; atom_len];
    let order_weight = |o: BondOrder| -> u8 { match o {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => if cfg.aromatic_as_one { 1 } else { 1 },
        BondOrder::Unknown => 0,
    }};
    for b in &mol.bonds {
        let (a, c) = (b.start_atom, b.end_atom);
        if let BondSymbol::Bond(ord) = b.symbol {
            let w = order_weight(ord);
            // Bond atoms are 1-based, convert to 0-based for array indexing
            let a_idx = (a as usize).saturating_sub(1);
            let c_idx = (c as usize).saturating_sub(1);
            if a_idx < atom_len { order_sum[a_idx] = order_sum[a_idx].saturating_add(w); }
            if c_idx < atom_len { order_sum[c_idx] = order_sum[c_idx].saturating_add(w); }
        }
    }
    let mut artifacts = ValenceArtifacts { atoms_checked: 0, overflow_count: 0, bracket_mismatch_count: 0 };
    let atom_element = |idx: usize| -> Option<Element> { match &mol.atoms[idx].symbol { AtomSymbol::Element(e) => Some(*e), _ => None } };
    for (i, atom) in mol.atoms.iter().enumerate() {
        let Some(elem) = atom_element(i) else { continue; };
        let charge_i32 = atom.charge.unwrap_or(0);
        let charge_i8 = charge_i32.clamp(i8::MIN as i32, i8::MAX as i32) as i8;
        let sum_orders_u8 = order_sum[i];
        let sum_orders = sum_orders_u8 as i32;
        let total_valence = sum_orders - charge_i32;
        artifacts.atoms_checked += 1;
        let mut implicit_h_opt: Option<i32> = None;
        
        // Get atom span for diagnostics: prefer IR start; fallback to whole input
        let atom_span = mol
            .atoms
            .get(i)
            .and_then(|a| a.span_start.map(|s| Span::new(s as usize, s as usize + 1)))
            .unwrap_or_else(|| Span::new(0, 0));
        let atom_id = (i + 1) as u32;
        
        if cfg.patterns_enabled {
            let d = resolve_valence_pattern(elem, sum_orders_u8, charge_i8, None, &model.patterns);
            if d.matches == 0 {
                match cfg.no_match_policy {
                    ValencePolicy::Off => {}
                    ValencePolicy::Warn => report.push(Diagnostic {
                        code: Code::NoMatch,
                        category: Category::Valence,
                        severity: Severity::Warning,
                        span: atom_span,
                        message: "No valence pattern matched",
                        details: Some(format!("atom_id={}", atom_id)),
                    }),
                    ValencePolicy::Error => report.push(Diagnostic {
                        code: Code::NoMatch,
                        category: Category::Valence,
                        severity: Severity::Error,
                        span: atom_span,
                        message: "No valence pattern matched",
                        details: Some(format!("atom_id={}", atom_id)),
                    }),
                }
            } else {
                if d.matches > 1 {
                    match cfg.ambiguous_match_policy {
                        ValencePolicy::Off => {}
                        ValencePolicy::Warn => report.push(Diagnostic {
                            code: Code::AmbiguousMatch,
                            category: Category::Valence,
                            severity: Severity::Warning,
                            span: atom_span,
                            message: "Multiple valence patterns matched; selected most specific",
                            details: Some(format!("atom_id={}", atom_id)),
                        }),
                        ValencePolicy::Error => report.push(Diagnostic {
                            code: Code::AmbiguousMatch,
                            category: Category::Valence,
                            severity: Severity::Error,
                            span: atom_span,
                            message: "Multiple valence patterns matched; selected most specific",
                            details: Some(format!("atom_id={}", atom_id)),
                        }),
                    }
                }
                if let Some(h) = d.implicit_h { implicit_h_opt = Some(h as i32); }
            }
        }
        
        if implicit_h_opt.is_none() {
            let states_opt = model.states_for(elem);
            if states_opt.is_none() { continue; }
            let states = states_opt.unwrap();
            if states.is_empty() { continue; }
            let mut chosen: Option<u8> = None;
            for &s in states { if (s as i32) >= total_valence { chosen = Some(s); break; } }
            if chosen.is_none() {
                artifacts.overflow_count += 1;
                match cfg.overflow_policy {
                    ValencePolicy::Off => {}
                    ValencePolicy::Warn => report.push(Diagnostic {
                        code: Code::ValenceOutOfElementRange,
                        category: Category::Valence,
                        severity: Severity::Warning,
                        span: atom_span,
                        message: "Valence exceeds maximum allowed state",
                        details: Some(format!("atom_id={}, total_valence={}", atom_id, total_valence)),
                    }),
                    ValencePolicy::Error => report.push(Diagnostic {
                        code: Code::ValenceOutOfElementRange,
                        category: Category::Valence,
                        severity: Severity::Error,
                        span: atom_span,
                        message: "Valence exceeds maximum allowed state",
                        details: Some(format!("atom_id={}, total_valence={}", atom_id, total_valence)),
                    }),
                }
                continue;
            }
            let chosen_state = chosen.unwrap() as i32;
            implicit_h_opt = Some(chosen_state - total_valence);
        }
        
        let implicit_h = implicit_h_opt.unwrap_or(0);
        if cfg.check_bracket {
            if let Some(h_explicit) = atom.hydrogen_count {
                let implied = implicit_h.max(0) as u32;
                if h_explicit != implied {
                    artifacts.bracket_mismatch_count += 1;
                    report.push(Diagnostic {
                        code: Code::HcountMismatch,
                        category: Category::Valence,
                        severity: Severity::Error,
                        span: atom_span,
                        message: "Bracket H count mismatches valence-based implicit H",
                        details: Some(format!(
                            "atom_id={}, expected_H={}, found_H={}",
                            atom_id, implied, h_explicit
                        )),
                    });
                }
            }
        }
    }
    artifacts
}


