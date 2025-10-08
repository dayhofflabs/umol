//! Post-parse SMILES checkers: re-exports of category modules

// pub mod aromaticity;
// pub mod linalg;
// pub mod stereo_double;
// pub mod topology;
// pub mod valence;

// pub use aromaticity::{
//     check_aromaticity, AromaticityArtifacts, AromaticityConfig, AromaticityMethod, AromaticityModel,
// };
// pub use stereo_double::{check_stereo_double, StereoArtifacts};
// pub use topology::{check_topology, TopologyArtifacts};
// pub use valence::{
//     check_valence, ValencePolicy, ValenceArtifacts, ValenceConfig, ValenceModel, ValencePattern,
//     ValencePatternTable,
// };

use super::diagnostics::DiagnosticsReport;
use super::parser::ParseMetadata;
use crate::io::ir::Molecule;
use crate::io::smiles::config::{SmilesParseFlags, SmilesCheckFlags, SmilesLintConfig, SmilesIoConfig};
use crate::io::smiles::parser::parse_smiles_inner;

#[derive(Default)]
pub struct SmilesModels;

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
    pub diagnostics: DiagnosticsReport,
    pub annotations: Annotations,
}

pub fn check_parsed(
    mol: &Molecule,
    meta: &ParseMetadata,
    check_flags: &SmilesCheckFlags,
    lint_config: &SmilesLintConfig,
    models: &SmilesModels,
) -> CheckOutput {
    let report = DiagnosticsReport::new();
    let annotations = Annotations::default();
    // check_topology(mol, side, &mut report, input_len);
    // check_stereo_double(mol, side, &mut report, input_len);
    CheckOutput {
        diagnostics: report,
        annotations: annotations,
    }
}

pub fn check_smiles(input: &[u8]) -> DiagnosticsReport {
    let io_config = SmilesIoConfig::default();
    let models = SmilesModels::default();
    check_smiles_with(input, &io_config, &models)
}

pub fn check_smiles_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    models: &SmilesModels,
) -> DiagnosticsReport {
    let parse_output = parse_smiles_inner(input, &io_config.parse_flags);
    match parse_output {
        Ok(parse_output) => {
            let mol = parse_output.sir;
            let meta = parse_output.meta.unwrap_or_default();
            let _report = DiagnosticsReport::new();
            let _annotations = Annotations::default();
            let check_output = check_parsed(&mol, &meta, &io_config.check_flags, &io_config.lint_config, models);
            check_output.diagnostics
        }
        Err(e) => {
            let mut report = DiagnosticsReport::new();
            report.push(e.as_diagnostic(""));
            report
        }
    }
}

#[cfg(test)]
mod tests;
