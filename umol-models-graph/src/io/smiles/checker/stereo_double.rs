use super::super::diagnostics::{Category, Diagnostic, DiagnosticsReport, Severity, Span};
use super::SideChannel;
use crate::io::ir::{BondDir, BondOrder, BondSymbol, Molecule};
use crate::io::smiles::diagnostics::DiagnosticCode;

pub struct StereoArtifacts {
    pub checked_double_bonds: usize,
    pub insufficient_count: usize,
    pub conflict_count: usize,
}

pub fn check_stereo_double(
    mol: &Molecule,
    _side: Option<&SideChannel>,
    report: &mut DiagnosticsReport,
    input_len: usize,
) -> StereoArtifacts {
    let mut artifacts = StereoArtifacts {
        checked_double_bonds: 0,
        insufficient_count: 0,
        conflict_count: 0,
    };

    // Helper: count oriented single-bond directions attached to `end_atom`, excluding the partner `other`
    let collect_dirs = |end_atom: u32, other: u32| -> (usize, usize) {
        let mut up = 0usize;
        let mut down = 0usize;
        for b in &mol.bonds {
            let (Some(s), Some(e)) = (b.start_atom, b.end_atom) else {
                continue;
            };
            if (s == end_atom && e == other) || (s == other && e == end_atom) {
                continue;
            }
            if s == end_atom || e == end_atom {
                if let BondSymbol::Bond(BondOrder::Single) = b.symbol {
                    if let Some(dir) = b.direction {
                        match dir {
                            BondDir::Up => up += 1,
                            BondDir::Down => down += 1,
                            _ => {}
                        }
                    }
                }
            }
        }
        (up, down)
    };

    for b in &mol.bonds {
        let (Some(a), Some(c)) = (b.start_atom, b.end_atom) else {
            continue;
        };
        // if let BondSymbol::Bond(BondOrder::Double) = b.symbol {
        //     artifacts.checked_double_bonds += 1;
        //     let (up_a, down_a) = collect_dirs(a, c);
        //     let (up_c, down_c) = collect_dirs(c, a);

        //     // Insufficient: need at least one oriented substituent on each end
        //     if (up_a + down_a) == 0 || (up_c + down_c) == 0 {
        //         artifacts.insufficient_count += 1;
        //         report.push(Diagnostic {
        //             code: DiagnosticCode::ParseError,
        //             category: Category::Stereo,
        //             severity: Severity::Error,
        //             span: Span::new(0, input_len),
        //             message: "Insufficient stereo information around double bond",
        //             details: None,
        //         });
        //         continue;
        //     }

        //     // Conflict: two or more oriented substituents with the same direction on the same end
        //     if up_a >= 2 || down_a >= 2 || up_c >= 2 || down_c >= 2 {
        //         artifacts.conflict_count += 1;
        //         report.push(Diagnostic {
        //             code: DiagnosticCode::ParseError,
        //             category: Category::Stereo,
        //             severity: Severity::Error,
        //             span: Span::new(0, input_len),
        //             message: "Conflicting same-direction stereo marks on one double-bond end",
        //             details: None,
        //         });
        //         continue;
        //     }
        // }
    }

    artifacts
}
