use std::collections::HashMap;

use crate::io::ir::Molecule;
use crate::io::smiles::config::SmilesCheckFlags;
use crate::io::smiles::diagnostics::{
    Category, Code, Diagnostic, Severity, Span,
};

pub fn check_topology(mol: &Molecule, flags: SmilesCheckFlags) -> impl Iterator<Item = Diagnostic> {
    let mut diags: Vec<Diagnostic> = Vec::new();
    if !flags.contains(SmilesCheckFlags::TOPOLOGY) {
        return diags.into_iter();
    }
    let mut edge_mult: HashMap<(u32, u32), Vec<u32>> = HashMap::new();

    for (bond_id, bond) in (1u32..).zip(mol.bonds.iter()) {
        let (atom1, atom2) = (bond.start_atom, bond.end_atom);
        if atom1 == atom2 {
            let span = bond
                .span_start
                .map(|s| Span::new(s as usize, s as usize + 1))
                .unwrap_or_else(|| Span::new(0, 0));
            diags.push(Diagnostic {
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

    for ((atom1, atom2), bond_ids) in edge_mult.into_iter() {
        if bond_ids.len() >= 2 {
            // Use the span of the first bond for diagnostic location; prefer IR start
            let span = if let Some(&first_id) = bond_ids.first() {
                let idx = (first_id - 1) as usize;
                mol.bonds
                    .get(idx)
                    .and_then(|b| b.span_start.map(|s| Span::new(s as usize, s as usize + 1)))
                    .unwrap_or_else(|| Span::new(0, 0))
            } else {
                Span::new(0, 0)
            };
            diags.push(Diagnostic {
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
    diags.into_iter()
}
