use std::collections::HashMap;

use crate::io::ir::Molecule;
use crate::io::smiles::diagnostics::{
    Category, Code, Diagnostic, DiagnosticsReport, Severity, Span,
};
use crate::io::smiles::parser::ParseMetadata;

pub struct TopologyArtifacts {
    pub self_loops: usize,
    pub parallel_pairs: usize,
}

pub fn check_topology(
    mol: &Molecule,
    _metadata: &ParseMetadata,
    report: &mut DiagnosticsReport,
    input_len: usize,
) -> TopologyArtifacts {
    let mut artifacts = TopologyArtifacts {
        self_loops: 0,
        parallel_pairs: 0,
    };
    // Track multiplicity of edges between unordered atom pairs
    let mut edge_mult: HashMap<(u32, u32), usize> = HashMap::new();

    for (index, bond) in mol.bonds.iter().enumerate() {
        let (atom1, atom2) = (bond.start_atom, bond.end_atom);
        if atom1 == atom2 {
            artifacts.self_loops += 1;
            report.push(Diagnostic {
                code: Code::SelfLoopRing,
                category: Category::Topology,
                severity: Severity::Error,
                span: Span::new(0, input_len),
                message: "Self-loop bond",
                details: Some(format!("bond_index={}", index)),
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
            .and_modify(|count| *count += 1)
            .or_insert(1);
    }

    for ((atom1, atom2), count) in edge_mult.into_iter() {
        if count >= 2 {
            artifacts.parallel_pairs += 1;
            report.push(Diagnostic {
                code: Code::ParallelEdges,
                category: Category::Topology,
                severity: Severity::Error,
                span: Span::new(0, input_len),
                message: "Multiple bonds between the same atom pair",
                details: Some(format!("atom1={atom1}, atom2={atom2}, count={count}")),
            });
        }
    }

    artifacts
}
