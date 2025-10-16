use std::collections::HashMap;

use crate::io::ir::Molecule;
use crate::io::smiles::checker::Annotations;
use crate::io::smiles::config::SmilesCheckFlags;
use crate::io::smiles::diagnostics::{Category, Code, Diagnostic, DiagnosticList, Severity, Span};

pub fn check_topology(
    mol: &Molecule,
    check_flags: &SmilesCheckFlags,
    _annotations: &mut Annotations,
    diagnostics: &mut DiagnosticList,
) {
    if !check_flags.contains(SmilesCheckFlags::TOPOLOGY) {
        return;
    }
    let mut edge_mult: HashMap<(u32, u32), Vec<u32>> = HashMap::new();

    // Self-loops
    for (bond_id, bond) in (1u32..).zip(mol.bonds.iter()) {
        let (atom1, atom2) = (bond.start_atom, bond.end_atom);
        if atom1 == atom2 {
            let span = match (bond.span_start, bond.span_end) {
                (Some(s), Some(e)) if e >= s => Span::new(s as usize, e as usize),
                (Some(s), _) => Span::new(s as usize, s as usize + 1),
                _ => Span::new(0, 0),
            };
            diagnostics.push(Diagnostic {
                code: Code::SelfLoopRing,
                category: Category::Topology,
                severity: Severity::Error,
                span,
                message: "Self-loop bond",
                details: Some(format!("bond_id={}", bond_id)),
            });
            continue;
        }
        let key = if atom1 < atom2 {
            (atom1, atom2)
        } else {
            (atom2, atom1)
        };
        edge_mult
            .entry(key)
            .and_modify(|ids| ids.push(bond_id))
            .or_insert_with(|| vec![bond_id]);
    }

    // Parallel edges
    for ((atom1, atom2), bond_ids) in edge_mult.into_iter() {
        if bond_ids.len() >= 2 {
            let span = if let Some(&first_id) = bond_ids.first() {
                let idx = (first_id - 1) as usize;
                if let Some(b) = mol.bonds.get(idx) {
                    match (b.span_start, b.span_end) {
                        (Some(s), Some(e)) if e >= s => Span::new(s as usize, e as usize),
                        (Some(s), _) => Span::new(s as usize, s as usize + 1),
                        _ => Span::new(0, 0),
                    }
                } else {
                    Span::new(0, 0)
                }
            } else {
                Span::new(0, 0)
            };
            diagnostics.push(Diagnostic {
                code: Code::ParallelEdges,
                category: Category::Topology,
                severity: Severity::Error,
                span,
                message: "Multiple bonds between the same atom pair",
                details: Some(format!(
                    "atom1={}, atom2={}, bond_ids={:?}",
                    atom1, atom2, bond_ids
                )),
            });
        }
    }
}
