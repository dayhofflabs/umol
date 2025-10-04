use std::collections::HashMap;

use super::super::diagnostics::{Category, Diagnostic, DiagnosticsReport, Severity, Span};
use super::linalg::hmo_density_from_adjacency;
use super::SideChannel;
use crate::io::ir::{BondOrder, BondSymbol, Molecule};
use crate::io::smiles::diagnostics::DiagnosticCode;

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
    let n_atoms = mol.atoms.len();
    let mut is_pi_center = vec![false; n_atoms];
    for b in &mol.bonds {
        let (Some(a), Some(c)) = (b.start_atom, b.end_atom) else {
            continue;
        };
        if let BondSymbol::Bond(ord) = b.symbol {
            match ord {
                BondOrder::Double | BondOrder::Aromatic => {
                    let ai = a as usize;
                    let ci = c as usize;
                    if ai < n_atoms {
                        is_pi_center[ai] = true;
                    }
                    if ci < n_atoms {
                        is_pi_center[ci] = true;
                    }
                }
                _ => {}
            }
        }
    }
    let mut atom_to_pi = vec![None; n_atoms];
    let mut idx = 0usize;
    for (i, flag) in is_pi_center.iter().copied().enumerate() {
        if flag {
            atom_to_pi[i] = Some(idx);
            idx += 1;
        }
    }
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for b in &mol.bonds {
        let (Some(a), Some(c)) = (b.start_atom, b.end_atom) else {
            continue;
        };
        let ai = a as usize;
        let ci = c as usize;
        if ai >= n_atoms || ci >= n_atoms {
            continue;
        }
        let (Some(pa), Some(pc)) = (atom_to_pi[ai], atom_to_pi[ci]) else {
            continue;
        };
        edges.push(if pa < pc { (pa, pc) } else { (pc, pa) });
    }
    PiGraph {
        size: idx,
        atom_to_pi,
        edges,
    }
}

pub fn check_aromaticity(
    mol: &Molecule,
    _side: Option<&SideChannel>,
    report: &mut DiagnosticsReport,
    input_len: usize,
    _amodel: &AromaticityModel,
    cfg: &AromaticityConfig,
) -> AromaticityArtifacts {
    if !cfg.enabled {
        return AromaticityArtifacts {
            checked_atoms: 0,
            pi_centers: 0,
            pi_edges: 0,
            method_used: cfg.method,
        };
    }
    let n_atoms = mol.atoms.len();
    let pg = build_pi_graph(mol);

    let has_aromatic_atoms = mol.atoms.iter().any(|a| a.aromatic == Some(true));
    let has_aromatic_bonds = mol
        .bonds
        .iter()
        .any(|b| matches!(b.symbol, BondSymbol::Bond(BondOrder::Aromatic)));
    let has_aromatic_tokens = has_aromatic_atoms || has_aromatic_bonds;

    let mut has_significant = false;
    let mut significant_count = 0usize;
    if matches!(
        cfg.method,
        AromaticityMethod::Hmo | AromaticityMethod::Combined
    ) {
        let p = hmo_density_from_adjacency(&pg.edges, pg.size);
        // Create bond order map for π-edge list
        let mut bmap: HashMap<(usize, usize), f32> = HashMap::new();
        for &(u, v) in &pg.edges {
            if u < pg.size && v < pg.size && u != v {
                let b = 2.0 * p[(u, v)] as f32;
                let key = if u < v { (u, v) } else { (v, u) };
                bmap.insert(key, b);
            }
        }
        let mut pi_to_atom: Vec<usize> = vec![0; pg.size];
        for (ai, opt_pi) in pg.atom_to_pi.iter().enumerate() {
            if let Some(pi) = opt_pi {
                if *pi < pi_to_atom.len() {
                    pi_to_atom[*pi] = ai;
                }
            }
        }
        let nearest_int = |x: f32| -> f32 {
            let r = x.round();
            if r < 0.0 {
                0.0
            } else if r > 2.0 {
                2.0
            } else {
                r
            }
        };
        let mut bond_lookup: HashMap<(usize, usize), (usize, BondOrder)> = HashMap::new();
        for (bi, b) in mol.bonds.iter().enumerate() {
            if let (Some(a), Some(c)) = (b.start_atom, b.end_atom) {
                let (ai, ci) = (a as usize, c as usize);
                let key = if ai < ci { (ai, ci) } else { (ci, ai) };
                if let BondSymbol::Bond(ord) = b.symbol {
                    bond_lookup.insert(key, (bi, ord));
                }
            }
        }
        for &(u, v) in &pg.edges {
            let key = if u < v { (u, v) } else { (v, u) };
            if let Some(&b) = bmap.get(&key) {
                let ai = pi_to_atom[u];
                let ci = pi_to_atom[v];
                let mkey = if ai < ci { (ai, ci) } else { (ci, ai) };
                if let Some((_bi, _ord)) = bond_lookup.get(&mkey) {
                    let diff = (b - nearest_int(b)).abs();
                    if diff >= cfg.fractional_bond_threshold {
                        significant_count += 1;
                    }
                }
            }
        }
        has_significant = significant_count > 0;
    }

    // if has_aromatic_tokens && (!has_significant || pg.size == 0) {
    //     report.push(Diagnostic {
    //         code: DiagnosticCode::,
    //         category: Category::Arom,
    //         severity: Severity::Warning,
    //         span: Span::new(0, input_len),
    //         message: "Aromatic annotations inconsistent with HMO delocalization",
    //         details: Some(format!(
    //             "pi_centers={}, significant_edges={}",
    //             pg.size, significant_count
    //         )),
    //     });
    // }
    // if !has_aromatic_tokens && has_significant {
    //     report.push(Diagnostic {
    //         code: DiagnosticCode::StylePreferAromaticForm,
    //         category: Category::Style,
    //         severity: Severity::Warning,
    //         span: Span::new(0, input_len),
    //         message: "Consider aromatic form for delocalized system",
    //         details: Some(format!("significant_edges={}", significant_count)),
    //     });
    // }

    AromaticityArtifacts {
        checked_atoms: n_atoms,
        pi_centers: pg.size,
        pi_edges: pg.edges.len(),
        method_used: cfg.method,
    }
}
