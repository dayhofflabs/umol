// Linting operations on GraphIR molecules.

use super::LintOutput;
use crate::graph_ir::{AtomIndex, BondIndex, Molecule};
use crate::io::smiles::config::{SmilesLintConfig, SmilesLintFlags};
use crate::io::smiles::diagnostics::{
    Category, Code, Diagnostic, DiagnosticList, EditList, Severity,
};
use crate::span::Span;
use crate::valence::ValenceModel;

pub(crate) fn lint_gir(
    gir: &Molecule,
    lint_flags: &SmilesLintFlags,
    _lint_config: &SmilesLintConfig,
    valence_model: &ValenceModel,
) -> LintOutput {
    let mut diagnostics = DiagnosticList::new();
    let edits = EditList::new();

    if lint_flags.contains(SmilesLintFlags::TOPOLOGY) {
        lint_gir_topology(gir, &mut diagnostics);
    }

    if lint_flags.contains(SmilesLintFlags::VALENCE) {
        lint_gir_valence(gir, valence_model, &mut diagnostics);
    }

    if lint_flags.contains(SmilesLintFlags::AROMATICITY) {
        lint_gir_aromaticity(gir, &mut diagnostics);
    }

    if lint_flags.contains(SmilesLintFlags::STEREO) {
        lint_gir_stereo(gir, &mut diagnostics);
    }

    LintOutput { diagnostics, edits }
}

fn lint_gir_topology(gir: &Molecule, diagnostics: &mut DiagnosticList) {
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
                details: Some(format!("atom={}", a_idx.index())),
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

        let bonds_display: Vec<usize> = bonds_between.iter().map(|idx| idx.index()).collect();

        diagnostics.push(Diagnostic {
            code: Code::ParallelEdges,
            category: Category::Topology,
            severity: Severity::Error,
            span: bond_span_gir(gir, bond_idx),
            message: "Multiple bonds between the same atom pair",
            details: Some(format!(
                "atom1={} atom2={} bonds={:?}",
                a_idx.index(),
                b_idx.index(),
                bonds_display
            )),
        });
    }
}

fn lint_gir_valence(
    gir: &Molecule,
    valence_model: &ValenceModel,
    diagnostics: &mut DiagnosticList,
) {
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

        if let Some(states) = valence_model.states_for(atom.element()) {
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
                        "atom={} element={:?} valence={} allowed={:?}",
                        atom_idx.index(),
                        atom.element(),
                        effective_valence,
                        states
                    )),
                });
            }
        }
    }
}

fn lint_gir_aromaticity(gir: &Molecule, diagnostics: &mut DiagnosticList) {
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
                    details: Some(format!("atom={} degree={}", atom_idx.index(), degree)),
                });
            }
        }
    }
}

fn lint_gir_stereo(_gir: &Molecule, _diagnostics: &mut DiagnosticList) {
    // Stereo validation for GIR will be implemented in a dedicated pass.
}

fn atom_span_gir(gir: &Molecule, idx: AtomIndex) -> Span {
    gir.atom(idx)
        .and_then(|atom| atom.span())
        .unwrap_or_else(|| Span::bytes(0, 0))
}

fn bond_span_gir(gir: &Molecule, idx: BondIndex) -> Span {
    gir.bond(idx)
        .and_then(|bond| bond.span())
        .unwrap_or_else(|| Span::bytes(0, 0))
}
