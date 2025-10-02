use std::collections::hash_map::Entry;
use std::collections::HashMap;

use crate::diagnostics::{Category, Code, DiagnosticsReport, Severity, Span};
use crate::io::ir::Molecule;

use super::SideChannel;

pub struct TopologyArtifacts {
    pub self_loops: usize,
    pub parallel_pairs: usize,
}

pub fn check_topology(
    mol: &Molecule,
    _side: Option<&SideChannel>,
    report: &mut DiagnosticsReport,
    input_len: usize,
) -> TopologyArtifacts {
    let mut artifacts = TopologyArtifacts {
        self_loops: 0,
        parallel_pairs: 0,
    };
    // Track multiplicity of edges between unordered atom pairs
    let mut edge_mult: HashMap<(u32, u32), (usize, (usize, usize))> = HashMap::new();

    for (bi, b) in mol.bonds.iter().enumerate() {
        let (Some(a), Some(c)) = (b.start_atom, b.end_atom) else {
            continue;
        };
        if a == c {
            artifacts.self_loops += 1;
            report.push(crate::diagnostics::Diagnostic {
                code: Code("TOPO_SELF_LOOP"),
                category: Category::Topology,
                severity: Severity::Error,
                span: Span::new(0, input_len),
                message: "Self-loop bond",
                details: Some(format!("bond_index={}", bi)),
            });
            continue;
        }
        let key = if a < c { (a, c) } else { (c, a) };
        match edge_mult.entry(key) {
            Entry::Vacant(v) => {
                v.insert((1, (a as usize, c as usize)));
            }
            Entry::Occupied(mut o) => {
                let (ref mut cnt, _ab) = o.get_mut();
                *cnt += 1;
            }
        }
    }

    for ((_a, _b), (cnt, _ab)) in edge_mult.into_iter() {
        if cnt >= 2 {
            artifacts.parallel_pairs += 1;
            report.push(crate::diagnostics::Diagnostic {
                code: Code("TOPO_PARALLEL_EDGES"),
                category: Category::Topology,
                severity: Severity::Error,
                span: Span::new(0, input_len),
                message: "Multiple bonds between the same atom pair",
                details: None,
            });
        }
    }

    artifacts
}


