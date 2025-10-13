use std::collections::HashMap;

use super::super::diagnostics::{Diagnostic, DiagnosticList};
use super::linalg::hmo_density_from_adjacency;
use super::SideChannel;
use crate::io::ir::{BondOrder, BondSymbol, Molecule};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AromaticityMethod {
    Hmo,
    Clar,
    Combined,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum HmoBackend {
    Nalgebra,
    Faer,
}

pub struct AromaticityConfig {
    pub enabled: bool,
    pub method: AromaticityMethod,
    pub hmo_backend: HmoBackend,
    pub fractional_bond_threshold: f32,
    pub annotate_symmetry: bool,
}

impl Default for AromaticityConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            method: AromaticityMethod::Hmo,
            hmo_backend: HmoBackend::Nalgebra,
            fractional_bond_threshold: 0.20,
            annotate_symmetry: false,
        }
    }
}

pub struct AromaticityModel {}

impl Default for AromaticityModel {
    fn default() -> Self {
        Self {}
    }
}

pub struct AromaticityArtifacts {
    pub checked_atoms: usize,
    pub pi_centers: usize,
    pub pi_edges: usize,
    pub method_used: AromaticityMethod,
}

struct PiGraph {
    size: usize,
    atom_to_pi: Vec<Option<usize>>, // map molecule atom index -> π index
    edges: Vec<(usize, usize)>,
}

fn build_pi_graph(mol: &Molecule) -> PiGraph {
    let _ = mol;
    PiGraph { size: 0, atom_to_pi: Vec::new(), edges: Vec::new() }
}

pub fn check_aromaticity(
    mol: &Molecule,
    _side: Option<&SideChannel>,
    report: &mut DiagnosticList,
    _input_len: usize,
    _model: &AromaticityModel,
    _cfg: &AromaticityConfig,
) -> AromaticityArtifacts {
    let _ = (mol, report);
    AromaticityArtifacts { checked_atoms: 0, pi_centers: 0, pi_edges: 0, method_used: AromaticityMethod::Hmo }
}
