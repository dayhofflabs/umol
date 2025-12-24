//! Linting operations on GraphIR molecules.

mod gir;

use serde::{Deserialize, Serialize};

use self::gir::lint_gir;
use crate::diagnostics::{Diagnostic, DiagnosticKind, DiagnosticList};
use crate::edits::EditList;
use crate::graph_ir::convert::sir_to_gir;
use crate::io::smiles::config::SmilesIoConfig;
use crate::io::smiles::parser::parse_smiles_to_sir;
use crate::valence::ValenceModel;

#[derive(Default, Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LintOutput {
    pub diagnostics: DiagnosticList,
    pub edits: EditList,
}

pub fn lint_smiles(input: &[u8]) -> LintOutput {
    let io_config = SmilesIoConfig::default();
    let valence_model = ValenceModel::simple_organic();
    lint_smiles_with(input, &io_config, &valence_model)
}

pub fn lint_smiles_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    valence_model: &ValenceModel,
) -> LintOutput {
    let parse_output = parse_smiles_to_sir(input, &io_config.parse_flags);
    match parse_output {
        Ok(parse_output) => {
            let mol = parse_output.ir;
            match sir_to_gir(&mol) {
                Ok(gir) => lint_gir(
                    &gir,
                    &io_config.lint_flags,
                    &io_config.lint_config,
                    valence_model,
                ),
                Err(_e) => LintOutput {
                    diagnostics: DiagnosticList::from(Diagnostic::from_kind(
                        DiagnosticKind::GraphConversionUnknown,
                    )),
                    edits: EditList::new(),
                },
            }
        }
        Err(e) => LintOutput {
            diagnostics: DiagnosticList::from(e.into()),
            edits: EditList::new(),
        },
    }
}
