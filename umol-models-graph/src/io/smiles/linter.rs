//! SMILES linter

pub mod aromaticity; // now located at io/smiles/linter/aromaticity.rs
pub mod stereo_double; // now located at io/smiles/linter/stereo_double.rs
pub mod stereo_chiral; // now located at io/smiles/linter/stereo_chiral.rs
pub mod topology; // now located at io/smiles/linter/topology.rs
pub mod valence; // now located at io/smiles/linter/valence.rs

pub use valence::lint_valence;

use crate::io::ir::Molecule;
use crate::io::smiles::config::{SmilesCheckFlags, SmilesIoConfig, SmilesLintConfig};
use crate::io::smiles::diagnostics::DiagnosticList;
use crate::io::smiles::parser::parse_smiles_to_ir;

use self::aromaticity::{lint_aromaticity, AromaticityConfig, AromaticityModel};
use self::stereo_chiral::lint_stereo_chiral;
use self::stereo_double::lint_stereo_double;
use self::topology::lint_topology;
use self::valence::{ValenceConfig, ValenceModel};

pub struct SmilesModels {
    pub valence: ValenceModel,
    pub aromaticity: AromaticityModel,
}

impl Default for SmilesModels {
    fn default() -> Self {
        Self {
            valence: ValenceModel::simple_organic(),
            aromaticity: AromaticityModel::default(),
        }
    }
}

#[derive(Default, Clone)]
pub struct TopologyAnnotations {
    pub self_loops: Vec<u32>,                       // bond_ids
    pub parallel_groups: Vec<(u32, u32, Vec<u32>)>, // (atom_a, atom_b, bond_ids)
    pub component_id: Vec<u32>,                     // per-atom component index (1-based group id)
    pub components: Vec<Vec<u32>>,                  // list of components as atom_ids
    pub bridges: Vec<u32>,                          // bond_ids
    pub articulation_points: Vec<u32>,              // atom_ids
    pub is_in_cycle_atom: Vec<bool>,                // length = num_atoms
    pub is_in_cycle_bond: Vec<bool>,                // length = num_bonds
}
#[derive(Default, Clone)]
pub struct ValenceAnnotations {
    pub total_valence_observed: Vec<i32>,      // per-atom
    pub implicit_h_suggested: Vec<Option<u8>>, // per-atom
    pub chosen_valence_state: Vec<Option<u8>>, // per-atom
    pub unknown_valence_atoms: Vec<u32>,       // atom_ids with unverified/unknown valence state
    pub bracket_mismatch_atoms: Vec<u32>,      // atom_ids
    // New neutral per-atom annotations
    pub valence_status: Vec<ValenceStatus>, // per-atom final status
    pub source: Vec<ValenceSource>,         // how the state was chosen
    pub match_count: Vec<u16>,              // number of matching patterns
    pub overflow_amount: Vec<Option<i32>>,  // (total_valence - max_allowed) when out of range
    pub required_bracket_h: Vec<Option<u8>>, // implied H if bracket omitted and >0
    pub has_unknown_bond_order: Vec<bool>,  // incident unknown bond order
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValenceStatus {
    Valid,
    NoMatch,
    Ambiguous,
    OutOfRange,
    MissingStates,
    BracketMismatch,
    UnknownBondOrder,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValenceSource {
    Pattern,
    Numeric,
}
pub struct AromaticityAnnotations;
#[derive(Default, Clone)]
pub struct StereoAnnotations {
    pub checked_double_bonds: Vec<u32>,
    pub candidates: Vec<u32>,
    pub insufficient: Vec<u32>,
    pub conflicts: Vec<u32>,
}

#[derive(Default)]
pub struct Annotations {
    pub topology: Option<TopologyAnnotations>,
    pub valence: Option<ValenceAnnotations>,
    pub arom: Option<AromaticityAnnotations>,
    pub stereo: Option<StereoAnnotations>,
}

#[derive(Default)]
pub struct CheckOutput {
    pub diagnostics: DiagnosticList,
    pub annotations: Annotations,
}

pub fn lint_ir(
    mol: &Molecule,
    check_flags: &SmilesCheckFlags,
    _lint_config: &SmilesLintConfig,
    models: &SmilesModels,
) -> CheckOutput {
    let mut diagnostics = DiagnosticList::new();
    let mut annotations = Annotations::default();
    
    // Run topology checks
    lint_topology(mol, check_flags, &mut annotations, &mut diagnostics);
    
    // Run valence checks
    let valence_cfg = ValenceConfig::default();
    lint_valence(
        mol,
        check_flags,
        &mut annotations,
        &mut diagnostics,
        &models.valence,
        &valence_cfg,
    );
    
    // Run stereo checks
    lint_stereo_double(mol, check_flags, &mut annotations, &mut diagnostics);
    lint_stereo_chiral(mol, check_flags, &mut annotations, &mut diagnostics);
    // Run aromaticity checks
    let arom_cfg = AromaticityConfig::default();
    lint_aromaticity(
        mol,
        check_flags,
        &mut annotations,
        &mut diagnostics,
        &models.aromaticity,
        &arom_cfg,
    );
    
    // STRICT mode: upgrade warnings to errors
    if check_flags.contains(SmilesCheckFlags::STRICT) {
        for d in diagnostics.diagnostics.iter_mut() {
            if d.severity == crate::io::smiles::diagnostics::Severity::Warning {
                d.severity = crate::io::smiles::diagnostics::Severity::Error;
            }
        }
    }
    
    CheckOutput {
        diagnostics,
        annotations,
    }
}

pub fn lint_smiles(input: &[u8]) -> DiagnosticList {
    let io_config = SmilesIoConfig::default();
    let models = SmilesModels::default();
    lint_smiles_with(input, &io_config, &models)
}

pub fn lint_smiles_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    models: &SmilesModels,
) -> DiagnosticList {
    let parse_output = parse_smiles_to_ir(input, &io_config.parse_flags);
    match parse_output {
        Ok(parse_output) => {
            let mol = parse_output.ir;
            let _report = DiagnosticList::new();
            let _annotations = Annotations::default();
            let check_output =
                lint_ir(&mol, &io_config.check_flags, &io_config.lint_config, models);
            check_output.diagnostics
        }
        Err(e) => DiagnosticList::from(e.as_diagnostic("")),
    }
}

#[cfg(test)]
mod tests;
