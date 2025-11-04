//! SMILES linter

mod gir;

use umol::error::DataError;
use umol::Error as UmolError;

use self::gir::{lint_gir, lint_topology_from_sir};
use crate::graph_ir::convert::sir_to_gir;
use crate::io::smiles::config::{SmilesIoConfig, SmilesLintFlags};
use crate::io::smiles::diagnostics::{
    Category, Code, Diagnostic, DiagnosticList, EditList, Severity,
};
use crate::io::smiles::parser::parse_smiles_to_sir;
use crate::simple_ir;
use crate::span::Span;
use crate::valence::ValenceModel;

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

#[derive(Default)]
pub struct LintOutput {
    pub diagnostics: DiagnosticList,
    pub edits: EditList,
}

pub fn lint_smiles(input: &[u8]) -> LintOutput {
    let io_config = SmilesIoConfig::default();
    let models = SmilesModels::default();
    lint_smiles_with(input, &io_config, &models)
}

pub fn lint_smiles_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    models: &SmilesModels,
) -> LintOutput {
    let parse_output = parse_smiles_to_sir(input, &io_config.parse_flags);
    match parse_output {
        Ok(parse_output) => {
            let mol = parse_output.ir;
            match sir_to_gir(&mol) {
                Ok(gir) => {
                    let mut output =
                        lint_gir(&gir, &io_config.lint_flags, &io_config.lint_config, models);
                    if io_config.lint_flags.contains(SmilesLintFlags::STRICT) {
                        output.diagnostics.upgrade_warnings();
                    }
                    output
                }
                Err(err) => {
                    let mut output = LintOutput {
                        diagnostics: conversion_error_diagnostics(err, &mol),
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

fn conversion_error_diagnostics(err: UmolError, sir: &simple_ir::Molecule) -> DiagnosticList {
    match err {
        UmolError::Data(DataError::LoopBond(atom_idx)) => {
            let mut diagnostics = DiagnosticList::new();
            lint_topology_from_sir(sir, &mut diagnostics);
            diagnostics
        }
        UmolError::Data(DataError::DuplicateBondIndex(a, b)) => {
            let mut diagnostics = DiagnosticList::new();
            lint_topology_from_sir(sir, &mut diagnostics);
            diagnostics
        }
        other => DiagnosticList::from(Diagnostic {
            code: Code::InternalError,
            category: Category::Internal,
            severity: Severity::Error,
            span: Span::bytes(0, 0),
            message: "GraphIR conversion failed",
            details: Some(other.to_string()),
        }),
    }
}
