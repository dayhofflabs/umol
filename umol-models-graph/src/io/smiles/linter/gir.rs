// Linting operations on GraphIR molecules.

use super::LintOutput;
use crate::diagnostics::{Diagnostic, DiagnosticKind, DiagnosticList, Severity};
use crate::edits::EditList;
use crate::graph_ir::{AtomIndex, BondIndex, Molecule};
use crate::io::smiles::config::{SmilesLintConfig, SmilesLintFlags};
use crate::span::Span;
use crate::valence::ValenceModel;

pub(super) fn lint_gir(
    molecule: &Molecule,
    lint_flags: &SmilesLintFlags,
    _lint_config: &SmilesLintConfig,
    valence_model: &ValenceModel,
) -> LintOutput {
    let mut diagnostics = DiagnosticList::new();
    let edits = EditList::new();

    if lint_flags.contains(SmilesLintFlags::TOPOLOGY) {
        lint_gir_topology(molecule, &mut diagnostics);
    }

    if lint_flags.contains(SmilesLintFlags::VALENCE) {
        lint_gir_valence(molecule, valence_model, &mut diagnostics);
    }

    if lint_flags.contains(SmilesLintFlags::AROMATICITY) {
        lint_gir_aromaticity(molecule, &mut diagnostics);
    }

    if lint_flags.contains(SmilesLintFlags::STEREO) {
        lint_gir_stereo(molecule, &mut diagnostics);
    }

    LintOutput { diagnostics, edits }
}

fn lint_gir_topology(molecule: &Molecule, diagnostics: &mut DiagnosticList) {
    for bond_idx in molecule.bond_indices() {
        let Some((a_idx, b_idx)) = molecule.bond_atom_indices(bond_idx) else {
            continue;
        };

        if a_idx == b_idx {
            diagnostics.push(Diagnostic::new(
                DiagnosticKind::GraphTopologySelfLoopRing,
                Severity::Error,
                Some(bond_span_by_idx(molecule, bond_idx)),
                Some(format!("atom={}", a_idx.index())),
            ));
            continue;
        }

        let mut bonds_between: Vec<_> = molecule.bonds_between(a_idx, b_idx).collect();
        if bonds_between.len() < 2 {
            continue;
        }

        bonds_between.sort_by_key(|idx| idx.index());
        if bond_idx != bonds_between[0] {
            continue;
        }

        let bonds_display: Vec<usize> = bonds_between.iter().map(|idx| idx.index()).collect();

        diagnostics.push(Diagnostic::new(
            DiagnosticKind::GraphTopologyParallelEdges,
            Severity::Error,
            Some(bond_span_by_idx(molecule, bond_idx)),
            Some(format!(
                "atom1={} atom2={} bonds={:?}",
                a_idx.index(),
                b_idx.index(),
                bonds_display
            )),
        ));
    }
}

fn lint_gir_valence(
    molecule: &Molecule,
    valence_model: &ValenceModel,
    diagnostics: &mut DiagnosticList,
) {
    for atom_idx in molecule.atom_indices() {
        let atom = match molecule.atom(atom_idx) {
            Some(atom) => atom,
            None => continue,
        };

        let bond_sum: u32 = molecule
            .atom_bond_indices(atom_idx)
            .filter_map(|bond_idx| molecule.bond(bond_idx))
            .map(|bond| u32::from(bond.order().value()))
            .sum();

        let implicit_hydrogens = atom.implicit_hydrogens();
        let effective_valence = bond_sum + implicit_hydrogens;

        if let Some(states) = valence_model.states_for(atom.element()) {
            if !states
                .iter()
                .any(|&state| u32::from(state) == effective_valence)
            {
                diagnostics.push(Diagnostic::new(
                    DiagnosticKind::GraphValenceOutOfElementRange,
                    Severity::Error,
                    Some(atom_span_by_idx(molecule, atom_idx)),
                    Some(format!(
                        "atom={} element={:?} valence={} allowed={:?}",
                        atom_idx.index(),
                        atom.element(),
                        effective_valence,
                        states
                    )),
                ));
            }
        }
    }
}

fn lint_gir_aromaticity(molecule: &Molecule, diagnostics: &mut DiagnosticList) {
    for atom_idx in molecule.atom_indices() {
        let atom = match molecule.atom(atom_idx) {
            Some(atom) => atom,
            None => continue,
        };

        if atom.aromatic() == Some(true) {
            let degree = molecule.atom_neighbor_indices(atom_idx).count();
            if degree < 2 {
                diagnostics.push(Diagnostic::new(
                    DiagnosticKind::GraphAromaticityAromaticAtomNotInRing,
                    Severity::Warning,
                    Some(atom_span_by_idx(molecule, atom_idx)),
                    Some(format!("atom={} degree={}", atom_idx.index(), degree)),
                ));
            }
        }
    }
}

fn lint_gir_stereo(_molecule: &Molecule, _diagnostics: &mut DiagnosticList) {
    // Stereo validation for GIR will be implemented in a dedicated pass.
}

fn atom_span_by_idx(molecule: &Molecule, idx: AtomIndex) -> Span {
    molecule
        .atom(idx)
        .and_then(|atom| atom.span())
        .unwrap_or_else(|| Span::bytes(0, 0))
}

fn bond_span_by_idx(molecule: &Molecule, idx: BondIndex) -> Span {
    molecule
        .bond(idx)
        .and_then(|bond| bond.span())
        .unwrap_or_else(|| Span::bytes(0, 0))
}
