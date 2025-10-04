//! Post-parse SMILES checkers: re-exports of category modules

pub mod aromaticity;
pub mod linalg;
pub mod stereo_double;
pub mod topology;
pub mod valence;

pub use aromaticity::{
    check_aromaticity, AromaticityArtifacts, AromaticityConfig, AromaticityMethod, AromaticityModel,
};
pub use stereo_double::{check_stereo_double, StereoArtifacts};
pub use topology::{check_topology, TopologyArtifacts};
pub use valence::{
    check_valence, ValencePolicy, ValenceArtifacts, ValenceConfig, ValenceModel, ValencePattern,
    ValencePatternTable,
};

// Placeholder for future side-channel (ring events, bond/atom spans, etc.)
pub struct SideChannel;

#[derive(Default)]
pub struct LintConfig {
    pub enabled_codes: Vec<&'static str>,
    pub disabled_codes: Vec<&'static str>,
}

pub struct ModelProfile;

pub struct CheckOptions<'a> {
    pub profile: &'a ModelProfile,
    pub lint: &'a LintConfig,
}

use super::diagnostics::DiagnosticsReport;
use crate::io::ir::Molecule;

pub fn check_smiles(
    mol: &Molecule,
    side: Option<&SideChannel>,
    input_len: usize,
    _opts: &CheckOptions,
) -> DiagnosticsReport {
    let mut report = DiagnosticsReport::new();
    let _topo = topology::check_topology(mol, side, &mut report, input_len);
    let _st = stereo_double::check_stereo_double(mol, side, &mut report, input_len);
    // valence/aromaticity are invoked from the linter today
    report
}

#[cfg(test)]
mod tests;
