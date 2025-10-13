//! Post-parse SMILES checkers: re-exports of category modules

// pub mod aromaticity;
// pub mod linalg;
// pub mod stereo_double;
pub mod topology;
pub mod valence;

// pub use aromaticity::{
//     check_aromaticity, AromaticityArtifacts, AromaticityConfig, AromaticityMethod, AromaticityModel,
// };
// pub use stereo_double::{check_stereo_double, StereoArtifacts};
// pub use topology::{check_topology, TopologyArtifacts};
// pub use valence::{
//     check_valence, ValencePolicy, ValenceArtifacts, ValenceConfig, ValenceModel, ValencePattern,
//     ValencePatternTable,
// };
pub use valence::check_valence;

use crate::io::ir::Molecule;
use crate::io::smiles::checker::topology::check_topology;
use crate::io::smiles::checker::valence::{ValenceConfig, ValenceModel};
use crate::io::smiles::config::{SmilesCheckFlags, SmilesIoConfig, SmilesLintConfig};
use crate::io::smiles::diagnostics::DiagnosticList;
use crate::io::smiles::parser::parse_smiles_inner;

pub struct SmilesModels {
    pub valence: ValenceModel,
}

impl Default for SmilesModels {
    fn default() -> Self {
        Self {
            valence: ValenceModel::simple_organic(),
        }
    }
}

pub struct TopologyAnnotations;
pub struct ValenceAnnotations;
pub struct AromaticityAnnotations;
pub struct StereoAnnotations;

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

pub fn check_parsed(
    mol: &Molecule,
    check_flags: &SmilesCheckFlags,
    _lint_config: &SmilesLintConfig,
    models: &SmilesModels,
) -> CheckOutput {
    let mut report = DiagnosticList::new();
    let annotations = Annotations::default();
    
    // Run topology checks
    report.extend(check_topology(mol, check_flags));
    
    // Run valence checks
    let valence_cfg = ValenceConfig::default();
    check_valence(mol, &mut report, &models.valence, &valence_cfg);
    
    // TODO: Run stereo checks
    // check_stereo_double(mol, metadata, &mut report);
    
    CheckOutput {
        diagnostics: report,
        annotations,
    }
}

pub fn check_smiles(input: &[u8]) -> DiagnosticList {
    let io_config = SmilesIoConfig::default();
    let models = SmilesModels::default();
    check_smiles_with(input, &io_config, &models)
}

pub fn check_smiles_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    models: &SmilesModels,
) -> DiagnosticList {
    let parse_output = parse_smiles_inner(input, &io_config.parse_flags);
    match parse_output {
        Ok(parse_output) => {
            let mol = parse_output.sir;
            let _report = DiagnosticList::new();
            let _annotations = Annotations::default();
            let check_output = check_parsed(
                &mol,
                &io_config.check_flags,
                &io_config.lint_config,
                models,
            );
            check_output.diagnostics
        }
        Err(e) => {
            let mut report = DiagnosticList::new();
            report.push(e.as_diagnostic(""));
            report
        }
    }
}

#[cfg(test)]
mod tests;
