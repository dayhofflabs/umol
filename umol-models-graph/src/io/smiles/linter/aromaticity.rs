use super::super::diagnostics::DiagnosticList;
use crate::io::ir::Molecule;
use crate::io::smiles::linter::Annotations;
use crate::io::smiles::config::SmilesCheckFlags;
use super::super::diagnostics::{Category, Code, Diagnostic, Severity, Span};
use crate::io::ir::AtomSymbol;
use umol_data::Element;
use std::fs::File;
use std::io::Read;
use std::path::Path;
use umol::{Error as UmolError, Result as UmolResult};
use umol::error::ParseError as UmolParseError;
use nalgebra::{DMatrix, SymmetricEigen};

pub fn lint_aromaticity(
    mol: &Molecule,
    check_flags: &SmilesCheckFlags,
    _annotations: &mut Annotations,
    diagnostics: &mut DiagnosticList,
    _model: &AromaticityModel,
    _cfg: &AromaticityConfig,
) {
    let _ = mol;
    if !check_flags.contains(SmilesCheckFlags::AROMATICITY) { return; }
    let pi = build_pi_graph(mol);
    let (pi_comp_id, comps) = pi_components(&pi);

    // Method/backend gating (stubbed config for now; TODO: wire from models/config)
    let cfg = AromaticityConfig::default();
    match cfg.method {
        AromaticityMethod::Hmo => {
            match cfg.hmo_backend {
                HmoBackend::Nalgebra => { let _ = hmo_eigenvalues(&pi); }
                HmoBackend::Faer => { /* TODO: implement faer backend */ }
            }
        }
        AromaticityMethod::Clar => { clar_stub(&pi); }
        AromaticityMethod::Combined => {
            match cfg.hmo_backend {
                HmoBackend::Nalgebra => { let _ = hmo_eigenvalues(&pi); }
                HmoBackend::Faer => { /* TODO: implement faer backend */ }
            }
            clar_stub(&pi);
        }
    }

    // Emit basic ring-consistency diagnostics for existing aromatic annotations
    // Atoms annotated aromatic must belong to a cyclic pi component
    for (atom_idx, atom) in mol.atoms.iter().enumerate() {
        if atom.aromatic == Some(true) {
            let span = match (atom.span_start, atom.span_end) {
                (Some(s), Some(e)) if e >= s => Span::new(s as usize, e as usize),
                (Some(s), _) => Span::new(s as usize, s as usize + 1),
                _ => Span::new(0, 0),
            };
            // Element eligibility check
            if let AtomSymbol::Element(e) = atom.symbol {
                if !supports_aromatic_element(e) {
                    diagnostics.push(Diagnostic {
                        code: Code::InvalidAromaticAtom,
                        category: Category::Aromaticity,
                        severity: Severity::Error,
                        span,
                        message: "Element is not eligible for aromaticity",
                        details: Some(format!("atom_id={}, element={:?}", (atom_idx as u32) + 1, e)),
                    });
                }
            }
            let comp_ok = if let Some(Some(pi_idx)) = pi.atom_to_pi.get(atom_idx) {
                let cid = pi_comp_id[*pi_idx];
                let (nodes, edges) = comps[cid];
                edges >= nodes // has at least one cycle
            } else { false };
            if !comp_ok {
                diagnostics.push(Diagnostic {
                    code: Code::AromaticAtomNotInRing,
                    category: Category::Aromaticity,
                    severity: Severity::Error,
                    span,
                    message: "Aromatic atom not part of any aromatic ring",
                    details: Some(format!("atom_id={}", (atom_idx as u32) + 1)),
                });
            }
        }
    }
    // Bonds annotated aromatic must lie in a cyclic pi component
    for (bond_id, b) in (1u32..).zip(mol.bonds.iter()) {
        use crate::io::ir::{BondOrder, BondSymbol};
        if !matches!(b.symbol, BondSymbol::Bond(BondOrder::Aromatic)) { continue; }
        let span = match (b.span_start, b.span_end) {
            (Some(s), Some(e)) if e >= s => Span::new(s as usize, e as usize),
            (Some(s), _) => Span::new(s as usize, s as usize + 1),
            _ => Span::new(0, 0),
        };
        let n = mol.atoms.len();
        let a = (b.start_atom as usize).saturating_sub(1);
        let c = (b.end_atom as usize).saturating_sub(1);
        let mut in_cyclic_pi = false;
        if a < n && c < n {
            // Endpoint eligibility
            let mut invalid_endpoint = false;
            if let Some(at) = mol.atoms.get(a) {
                if let AtomSymbol::Element(e) = at.symbol { if !supports_aromatic_element(e) { invalid_endpoint = true; } }
            }
            if let Some(ct) = mol.atoms.get(c) {
                if let AtomSymbol::Element(e) = ct.symbol { if !supports_aromatic_element(e) { invalid_endpoint = true; } }
            }
            if invalid_endpoint {
                diagnostics.push(Diagnostic {
                    code: Code::InvalidAromaticBondAtom,
                    category: Category::Aromaticity,
                    severity: Severity::Error,
                    span,
                    message: "Aromatic bond has non-eligible endpoint atom",
                    details: Some(format!("bond_id={}", bond_id)),
                });
            }
            if let (Some(Some(i)), Some(Some(j))) = (pi.atom_to_pi.get(a), pi.atom_to_pi.get(c)) {
                let cid_i = pi_comp_id[*i];
                let cid_j = pi_comp_id[*j];
                if cid_i == cid_j {
                    let (nodes, edges) = comps[cid_i];
                    in_cyclic_pi = edges >= nodes;
                }
            }
        }
        if !in_cyclic_pi {
            diagnostics.push(Diagnostic {
                code: Code::AromaticBondNotInRing,
                category: Category::Aromaticity,
                severity: Severity::Error,
                span,
                message: "Aromatic bond not part of any aromatic ring",
                details: Some(format!("bond_id={}", bond_id)),
            });
        }
    }

    // Component-level style/consistency checks and Hückel criterion (4n+2 vs 4n)
    let mut comp_nodes: Vec<Vec<usize>> = vec![Vec::new(); comps.len()];
    for (atom_idx, map) in pi.atom_to_pi.iter().enumerate() {
        if let Some(pi_idx) = map {
            let cid = pi_comp_id[*pi_idx];
            if cid < comp_nodes.len() { comp_nodes[cid].push(atom_idx); }
        }
    }
    for cid in 0..comps.len() {
        let (nodes, edges) = comps[cid];
        if nodes == 0 { continue; }
        let mut aromatic_atoms = 0usize;
        for &ai in &comp_nodes[cid] {
            if mol.atoms[ai].aromatic == Some(true) { aromatic_atoms += 1; }
        }
        let mut aromatic_bonds = 0usize;
        // count aromatic bonds within comp by scanning edges mapped from bonds
        for b in &mol.bonds {
            use crate::io::ir::{BondOrder, BondSymbol};
            if !matches!(b.symbol, BondSymbol::Bond(BondOrder::Aromatic)) { continue; }
            let a = (b.start_atom as usize).saturating_sub(1);
            let c = (b.end_atom as usize).saturating_sub(1);
            if let (Some(Some(i)), Some(Some(j))) = (pi.atom_to_pi.get(a), pi.atom_to_pi.get(c)) {
                if pi_comp_id[*i] == cid && pi_comp_id[*j] == cid { aromatic_bonds += 1; }
            }
        }
        if aromatic_bonds > 0 && aromatic_atoms < comp_nodes[cid].len() {
            // Mixed/inconsistent aromaticity style warning
            let span = comp_nodes[cid]
                .first()
                .and_then(|&ai| mol.atoms[ai].span_start.map(|s| Span::new(s as usize, s as usize + 1)))
                .unwrap_or_else(|| Span::new(0, 0));
            diagnostics.push(Diagnostic {
                code: Code::AvoidInconsistentAromaticity,
                category: Category::Aromaticity,
                severity: Severity::Warning,
                span,
                message: "Inconsistent aromaticity within π-component",
                details: Some(format!("component_id={}, aromatic_atoms={}/{}", cid, aromatic_atoms, comp_nodes[cid].len())),
            });
        }
        // Huckel criterion: simplistic per-component check
        if edges >= nodes && aromatic_bonds > 0 {
            let ok = (nodes % 4) == 2;
            if !ok {
                let span = comp_nodes[cid]
                    .first()
                    .and_then(|&ai| mol.atoms[ai].span_start.map(|s| Span::new(s as usize, s as usize + 1)))
                    .unwrap_or_else(|| Span::new(0, 0));
                diagnostics.push(Diagnostic {
                    code: Code::HuckelFail,
                    category: Category::Aromaticity,
                    severity: Severity::Error,
                    span,
                    message: "Hückel (4n+2) rule not satisfied for aromatic component",
                    details: Some(format!("component_id={}, nodes={}, edges={}", cid, nodes, edges)),
                });
            }
        }
    }
}

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

#[derive(Default, serde::Deserialize, serde::Serialize)]
pub struct AromaticityPatternTable {
    pub eligible_elements: Vec<Element>,
}

pub struct AromaticityModel {
    pub patterns: AromaticityPatternTable,
}

impl Default for AromaticityModel {
    fn default() -> Self {
        Self { patterns: AromaticityPatternTable { eligible_elements: vec![
            Element::C, Element::N, Element::O, Element::S, Element::P, Element::Si, Element::Ge, Element::Sn, Element::B,
        ] } }
    }
}

impl AromaticityModel {
    pub fn with_patterns(patterns: AromaticityPatternTable) -> Self { Self { patterns } }
    pub fn set_patterns(&mut self, patterns: AromaticityPatternTable) { self.patterns = patterns; }
    pub fn from_patterns_reader<R: Read>(mut r: R) -> UmolResult<Self> {
        let mut buf = String::new();
        r.read_to_string(&mut buf)
            .map_err(|e| UmolError::from(UmolParseError::IoError(e)))?;
        let table: AromaticityPatternTable = toml::from_str(&buf)
            .map_err(|e| UmolError::from(UmolParseError::Invalid(format!("TOML parse error: {}", e))))?;
        Ok(Self::with_patterns(table))
    }
    pub fn from_patterns_file(path: &Path) -> UmolResult<Self> {
        let f = File::open(path)
            .map_err(|e| UmolError::from(UmolParseError::IoError(e)))?;
        Self::from_patterns_reader(f)
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
    use crate::io::ir::{BondOrder, BondSymbol};
    let n = mol.atoms.len();
    let mut marked: Vec<bool> = vec![false; n];
    for b in &mol.bonds {
        let is_pi = matches!(b.symbol, BondSymbol::Bond(BondOrder::Double) | BondSymbol::Bond(BondOrder::Aromatic));
        if !is_pi { continue; }
        let a = (b.start_atom as usize).saturating_sub(1);
        let c = (b.end_atom as usize).saturating_sub(1);
        if a < n { marked[a] = true; }
        if c < n { marked[c] = true; }
    }
    let mut atom_to_pi = vec![None; n];
    let mut idx = 0usize;
    for i in 0..n { if marked[i] { atom_to_pi[i] = Some(idx); idx += 1; } }
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for b in &mol.bonds {
        let is_pi = matches!(b.symbol, BondSymbol::Bond(BondOrder::Double) | BondSymbol::Bond(BondOrder::Aromatic));
        if !is_pi { continue; }
        let a = (b.start_atom as usize).saturating_sub(1);
        let c = (b.end_atom as usize).saturating_sub(1);
        if a < n && c < n {
            if let (Some(i), Some(j)) = (atom_to_pi[a], atom_to_pi[c]) {
                if i != j { edges.push((i, j)); }
            }
        }
    }
    PiGraph { size: idx, atom_to_pi, edges }
}

fn pi_components(pi: &PiGraph) -> (Vec<usize>, Vec<(usize, usize)>) {
    let n = pi.size;
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(u, v) in &pi.edges {
        if u < n && v < n {
            adj[u].push(v);
            adj[v].push(u);
        }
    }
    let mut vis = vec![false; n];
    let mut comps = Vec::new();
    let mut comp_id = vec![0usize; n];
    let mut cid = 0usize;
    for s in 0..n {
        if vis[s] { continue; }
        let this_cid = cid; cid += 1;
        let mut stack = vec![s];
        let mut nodes = 0usize;
        let mut edges = 0usize;
        vis[s] = true;
        while let Some(u) = stack.pop() {
            nodes += 1;
            edges += adj[u].len();
            comp_id[u] = this_cid;
            for &v in &adj[u] {
                if !vis[v] { vis[v] = true; stack.push(v); }
            }
        }
        // Each edge is counted twice in undirected adjacency
        comps.push((nodes, edges / 2));
    }
    (comp_id, comps)
}

fn hmo_eigenvalues(pi: &PiGraph) -> Vec<f64> {
    let n = pi.size;
    if n == 0 { return Vec::new(); }
    let mut a = DMatrix::<f64>::zeros(n, n);
    for &(u, v) in &pi.edges {
        if u < n && v < n {
            a[(u, v)] = 1.0;
            a[(v, u)] = 1.0;
        }
    }
    let eig = SymmetricEigen::new(a);
    let mut vals = eig.eigenvalues.data.as_vec().clone();
    vals.sort_by(|x, y| x.partial_cmp(y).unwrap());
    vals
}

fn supports_aromatic_element(e: Element) -> bool {
    match e {
        Element::C | Element::N | Element::O | Element::S | Element::P | Element::Si | Element::Ge | Element::Sn | Element::B => true,
        _ => false,
    }
}

fn clar_stub(pi: &PiGraph) {
    let _ = pi;
    // Placeholder for Clar sextet analysis; to be implemented in a future todo
}


