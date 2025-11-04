// GIR linting logic operating on GraphIR molecules.

use std::collections::HashMap;

use super::{LintOutput, SmilesModels};
use crate::graph_ir::molecule::{AtomIndex, BondIndex};
use crate::graph_ir::Molecule as GraphMolecule;
use crate::io::smiles::config::{SmilesLintConfig, SmilesLintFlags};
use crate::io::smiles::diagnostics::{
    Category, Code, Diagnostic, DiagnosticList, EditList, Severity,
};
use crate::simple_ir;
use crate::span::Span;

pub(crate) fn lint_gir(
    gir: &GraphMolecule,
    lint_flags: &SmilesLintFlags,
    _lint_config: &SmilesLintConfig,
    models: &SmilesModels,
) -> LintOutput {
    let mut diagnostics = DiagnosticList::new();
    let edits = EditList::new();

    if lint_flags.contains(SmilesLintFlags::TOPOLOGY) {
        lint_gir_topology(gir, &mut diagnostics);
    }

    if lint_flags.contains(SmilesLintFlags::VALENCE) {
        lint_gir_valence(gir, models, &mut diagnostics);
    }

    if lint_flags.contains(SmilesLintFlags::AROMATICITY) {
        lint_gir_aromaticity(gir, &mut diagnostics);
    }

    if lint_flags.contains(SmilesLintFlags::STEREO) {
        lint_gir_stereo(gir, &mut diagnostics);
    }

    LintOutput { diagnostics, edits }
}

fn lint_gir_topology(gir: &GraphMolecule, diagnostics: &mut DiagnosticList) {
    for bond_idx in gir.bond_indices() {
        let Some((a_idx, b_idx)) = gir.bond_atom_indices(bond_idx) else {
            continue;
        };

        if a_idx == b_idx {
            diagnostics.push(Diagnostic {
                code: Code::SelfLoopRing,
                category: Category::Topology,
                severity: Severity::Error,
                span: bond_span_gir(gir, bond_idx),
                message: "Self-loop bond",
                details: Some(format!("atom_id={}", display_atom(a_idx))),
            });
            continue;
        }

        let mut bonds_between: Vec<_> = gir.bonds_between(a_idx, b_idx).collect();
        if bonds_between.len() < 2 {
            continue;
        }

        bonds_between.sort_by_key(|idx| idx.index());
        if bond_idx != bonds_between[0] {
            continue;
        }

        let bonds_display: Vec<usize> =
            bonds_between.iter().map(|idx| display_bond(*idx)).collect();

        diagnostics.push(Diagnostic {
            code: Code::ParallelEdges,
            category: Category::Topology,
            severity: Severity::Error,
            span: bond_span_gir(gir, bond_idx),
            message: "Multiple bonds between the same atom pair",
            details: Some(format!(
                "atom1={} atom2={} bonds={:?}",
                display_atom(a_idx),
                display_atom(b_idx),
                bonds_display
            )),
        });
    }
}

fn lint_gir_valence(gir: &GraphMolecule, models: &SmilesModels, diagnostics: &mut DiagnosticList) {
    for atom_idx in gir.atom_indices() {
        let atom = match gir.atom(atom_idx) {
            Some(atom) => atom,
            None => continue,
        };

        let bond_sum: u32 = gir
            .atom_bond_indices(atom_idx)
            .filter_map(|bond_idx| gir.bond(bond_idx))
            .map(|bond| u32::from(bond.order().value()))
            .sum();

        let implicit_h = atom.implicit_h();
        let effective_valence = bond_sum + implicit_h;

        if let Some(states) = models.valence.states_for(atom.element()) {
            if !states
                .iter()
                .any(|&state| u32::from(state) == effective_valence)
            {
                diagnostics.push(Diagnostic {
                    code: Code::ValenceOutOfElementRange,
                    category: Category::Valence,
                    severity: Severity::Error,
                    span: atom_span_gir(gir, atom_idx),
                    message: "Observed valence is not permitted for this element",
                    details: Some(format!(
                        "atom_id={} element={:?} valence={} allowed={:?}",
                        display_atom(atom_idx),
                        atom.element(),
                        effective_valence,
                        states
                    )),
                });
            }
        }
    }
}

fn lint_gir_aromaticity(gir: &GraphMolecule, diagnostics: &mut DiagnosticList) {
    for atom_idx in gir.atom_indices() {
        let atom = match gir.atom(atom_idx) {
            Some(atom) => atom,
            None => continue,
        };

        if atom.aromatic() == Some(true) {
            let degree = gir.atom_neighbor_indices(atom_idx).count();
            if degree < 2 {
                diagnostics.push(Diagnostic {
                    code: Code::AromaticAtomNotInRing,
                    category: Category::Aromaticity,
                    severity: Severity::Warning,
                    span: atom_span_gir(gir, atom_idx),
                    message: "Aromatic atom does not participate in a ring",
                    details: Some(format!(
                        "atom_id={} degree={}",
                        display_atom(atom_idx),
                        degree
                    )),
                });
            }
        }
    }
}

fn lint_gir_stereo(_gir: &GraphMolecule, _diagnostics: &mut DiagnosticList) {
    // Stereo validation for GIR will be implemented in a dedicated pass.
}

fn atom_span_gir(gir: &GraphMolecule, idx: AtomIndex) -> Span {
    gir.atom(idx)
        .and_then(|atom| atom.span())
        .unwrap_or_else(|| Span::bytes(0, 0))
}

fn bond_span_gir(gir: &GraphMolecule, idx: BondIndex) -> Span {
    gir.bond(idx)
        .and_then(|bond| bond.span())
        .unwrap_or_else(|| Span::bytes(0, 0))
}

fn display_atom(idx: AtomIndex) -> usize {
    idx.index() + 1
}

fn display_bond(idx: BondIndex) -> usize {
    idx.index() + 1
}

fn atom_span(sir: &simple_ir::Molecule, idx: usize) -> Span {
    sir.atoms
        .get(idx)
        .and_then(|atom| atom.span)
        .unwrap_or_else(|| Span::bytes(0, 0))
}

fn bond_span(sir: &simple_ir::Molecule, idx: usize) -> Span {
    sir.bonds
        .get(idx)
        .and_then(|bond| bond.span)
        .unwrap_or_else(|| Span::bytes(0, 0))
}

fn canonical_pair(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

pub(crate) fn lint_topology_from_sir(sir: &simple_ir::Molecule, diagnostics: &mut DiagnosticList) {
    let mut parallel_edges: HashMap<(usize, usize), Vec<usize>> = HashMap::new();

    for (bond_idx, bond) in sir.bonds.iter().enumerate() {
        let Some(start) = usize::try_from(bond.start_atom).ok() else {
            continue;
        };
        let Some(end) = usize::try_from(bond.end_atom).ok() else {
            continue;
        };

        if start == end {
            diagnostics.push(Diagnostic {
                code: Code::SelfLoopRing,
                category: Category::Topology,
                severity: Severity::Error,
                span: bond_span(sir, bond_idx),
                message: "Self-loop bond",
                details: Some(format!("atom_id={}", start + 1)),
            });
            continue;
        }

        parallel_edges
            .entry(canonical_pair(start, end))
            .or_default()
            .push(bond_idx);
    }

    for ((a, b), bond_indices) in parallel_edges.into_iter() {
        if bond_indices.len() < 2 {
            continue;
        }

        let span = bond_span(sir, bond_indices[0]);
        let bonds_display: Vec<usize> = bond_indices.iter().map(|idx| idx + 1).collect();

        diagnostics.push(Diagnostic {
            code: Code::ParallelEdges,
            category: Category::Topology,
            severity: Severity::Error,
            span,
            message: "Multiple bonds between the same atom pair",
            details: Some(format!(
                "atom1={} atom2={} bonds={:?}",
                a + 1,
                b + 1,
                bonds_display
            )),
        });
    }
}
