use super::super::diagnostics::{Diagnostic, DiagnosticList};
use super::SideChannel;
use crate::io::ir::{BondDir, BondOrder, BondSymbol, Molecule};

pub struct StereoArtifacts {
    pub checked_double_bonds: usize,
    pub insufficient_count: usize,
    pub conflict_count: usize,
}

pub fn check_stereo_double(
    mol: &crate::io::ir::Molecule,
    _side: Option<&()>,
    report: &mut DiagnosticList,
    _input_len: usize,
) {
    let _ = (mol, report);
}
