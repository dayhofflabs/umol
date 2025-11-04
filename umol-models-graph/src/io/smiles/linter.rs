//! Linting operations on GraphIR molecules.

mod gir;

use self::gir::lint_gir;
use crate::graph_ir::convert::sir_to_gir;
use crate::io::smiles::config::{SmilesIoConfig, SmilesLintFlags};
use crate::io::smiles::diagnostics::{
    Category, Code, Diagnostic, DiagnosticList, EditList, Severity,
};
use crate::io::smiles::parser::parse_smiles_to_sir;
use crate::span::Span;
use crate::valence::ValenceModel;

#[derive(Default)]
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
                Ok(gir) => {
                    let mut output = lint_gir(
                        &gir,
                        &io_config.lint_flags,
                        &io_config.lint_config,
                        valence_model,
                    );
                    if io_config.lint_flags.contains(SmilesLintFlags::STRICT) {
                        output.diagnostics.upgrade_warnings();
                    }
                    output
                }
                Err(err) => {
                    let mut output = LintOutput {
                        diagnostics: DiagnosticList::from(Diagnostic {
                            code: Code::InternalError,
                            category: Category::Internal,
                            severity: Severity::Error,
                            span: Span::bytes(0, 0),
                            message: "GraphIR conversion failed",
                            details: Some(err.to_string()),
                        }),
                        edits: EditList::new(),
                    };
                    if io_config.lint_flags.contains(SmilesLintFlags::STRICT) {
                        output.diagnostics.upgrade_warnings();
                    }
                    output
                }
            }
        }
        Err(e) => {
            let mut output = LintOutput {
                diagnostics: DiagnosticList::from(e.as_diagnostic("")),
                edits: EditList::new(),
            };
            if io_config.lint_flags.contains(SmilesLintFlags::STRICT) {
                output.diagnostics.upgrade_warnings();
            }
            output
        }
    }
}
