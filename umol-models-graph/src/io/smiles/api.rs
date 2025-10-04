//! SMILES API scaffolding: stage outputs, annotations, and helpers

use crate::io::ir::Molecule;
use crate::io::smiles::diagnostics::{DiagnosticCategory, Diagnostic, DiagnosticCode, DiagnosticsReport, Span};
// use super::checker::{
//     check_aromaticity, check_stereo_double, check_topology, check_valence, AromaticityConfig,
//     AromaticityModel, SideChannel, ValenceConfig, ValenceModel,
// };
use crate::io::smiles::parse_smiles;

/// Lightweight parse metadata (explicit side-channel replacement)
#[derive(Debug, Default, Clone)]
pub struct ParseMeta {
    pub token_spans: Vec<(usize, usize)>,
    pub ring_events: Vec<u32>,
}

/// Parse stage output (scaffold; parser returns Molecule today)
#[derive(Debug, Default, Clone)]
pub struct ParseOutput {
    pub sir: Option<Molecule>,
    pub meta: ParseMeta,
    pub diagnostics: DiagnosticsReport,
}

/// Typed annotations produced by checker passes (scaffold)
#[derive(Debug, Default, Clone)]
pub struct Annotations {
    pub has_topology: bool,
    pub has_valence: bool,
    pub has_aromaticity: bool,
    pub has_stereo: bool,
}

/// Check stage output (scaffold)
#[derive(Debug, Default, Clone)]
pub struct CheckOutput {
    pub annotations: Annotations,
    pub diagnostics: DiagnosticsReport,
}

/// Graph stage output (scaffold)
#[derive(Debug, Default, Clone)]
pub struct GraphOutput {}

/// Convenience: wrap existing parse into scaffolding (no behavior change)
pub fn parse_smiles_meta(input: &str) -> ParseOutput {
    let mut out = ParseOutput::default();
    match parse_smiles(input.as_bytes()) {
        Ok(m) => out.sir = Some(m),
        Err(e) => {
            // Map will be centralized later; store a terse detail for now
            // let span = Span::new(0, input.len());
            // out.diagnostics.push(Diagnostic::error(
                
            //     Category::Syn,
            //     span,
            //     "SMILES parse error",
            // ));
        }
    }
    out
}

// /// Run standard checks and produce typed annotations (scaffold wrapper)
// pub fn check_and_annotate(mol: &Molecule, input_len: usize) -> CheckOutput {
//     let mut report = DiagnosticsReport::new();
//     let side: Option<&SideChannel> = None;

//     // Topology
//     let _topo = check_topology(mol, side, &mut report, input_len);

//     // Valence
//     let v_cfg = ValenceConfig::default();
//     let v_model = ValenceModel::simple_organic();
//     let _val = check_valence(mol, side, &mut report, input_len, &v_model, &v_cfg);

//     // Stereo
//     let _st = check_stereo_double(mol, side, &mut report, input_len);

//     // Aromaticity
//     let a_cfg = AromaticityConfig::default();
//     let a_model = AromaticityModel::default();
//     let _ar = check_aromaticity(mol, side, &mut report, input_len, &a_model, &a_cfg);

//     CheckOutput {
//         annotations: Annotations {
//             has_topology: true,
//             has_valence: true,
//             has_aromaticity: true,
//             has_stereo: true,
//         },
//         diagnostics: report,
//     }
// }


