use super::super::diagnostics::{Category, Code, Diagnostic, DiagnosticList, Severity, Span};
use crate::io::ir::{AtomSymbol, Molecule};
use crate::io::smiles::checker::{Annotations, StereoAnnotations};
use crate::io::smiles::config::SmilesCheckFlags;
use umol_data::Element;

pub fn check_stereo_chiral(
    mol: &Molecule,
    check_flags: &SmilesCheckFlags,
    annotations: &mut Annotations,
    diagnostics: &mut DiagnosticList,
) {
    if !check_flags.contains(SmilesCheckFlags::STEREO) { return; }
    if annotations.stereo.is_none() { annotations.stereo = Some(StereoAnnotations::default()); }
    let _st = annotations.stereo.as_mut().unwrap();

    for (atom_id, atom) in (1u32..).zip(mol.atoms.iter()) {
        let span = match (atom.span_start, atom.span_end) {
            (Some(s), Some(e)) if e >= s => Span::new(s as usize, e as usize),
            (Some(s), _) => Span::new(s as usize, s as usize + 1),
            _ => Span::new(0, 0),
        };
        // Only consider element atoms
        let element = match atom.symbol { AtomSymbol::Element(e) => e, _ => continue };
        let deg = degree(mol, atom_id);
        // Unsupported central chirality (very rough heuristic; expand later)
        if !supports_central_chirality(element) && atom.chirality.is_some() {
            diagnostics.push(Diagnostic {
                code: Code::UnsupportedCentralChiralityElement,
                category: Category::Stereo,
                severity: Severity::Warning,
                span,
                message: "Element typically does not support central chirality",
                details: Some(format!("atom_id={}, element={:?}", atom_id, element)),
            });
        }
        // Substituent count vs chirality kind (rough)
        if let Some(_chi) = atom.chirality {
            if deg < 3 {
                diagnostics.push(Diagnostic {
                    code: Code::ChiralitySubstituentMismatch,
                    category: Category::Stereo,
                    severity: Severity::Warning,
                    span,
                    message: "Chirality annotation incompatible with substituent count",
                    details: Some(format!("atom_id={}, degree={}", atom_id, deg)),
                });
            }
        }
        // Non-chiral despite annotation (placeholder: degenerate neighbors)
        if atom.chirality.is_some() && !is_potentially_stereogenic_atom(mol, atom_id) {
            diagnostics.push(Diagnostic {
                code: Code::NonChiralAnnotated,
                category: Category::Stereo,
                severity: Severity::Warning,
                span,
                message: "Atom annotated as chiral but appears non-stereogenic",
                details: Some(format!("atom_id={}", atom_id)),
            });
        }
        // Stereogenic annotation (neutral): use stereo annotations to flag candidates
        if is_potentially_stereogenic_atom(mol, atom_id) {
            // Could add a dedicated list in StereoAnnotations later
            // For now, no-op; this is a stub per request
        }
    }
}

fn degree(mol: &Molecule, atom_id: u32) -> usize {
    let idx = (atom_id as usize).saturating_sub(1);
    let mut d = 0usize;
    for b in &mol.bonds {
        if b.start_atom == atom_id || b.end_atom == atom_id { d += 1; }
    }
    d
}

fn supports_central_chirality(element: Element) -> bool {
    match element {
        Element::C | Element::N | Element::P | Element::S | Element::Si | Element::Ge | Element::Sn => true,
        _ => false,
    }
}

fn is_potentially_stereogenic_atom(mol: &Molecule, atom_id: u32) -> bool {
    // Cheap proxy: degree >= 4 (or 3 for planar N), and at least two distinct neighbor ids
    let idx = (atom_id as usize).saturating_sub(1);
    if idx >= mol.atoms.len() { return false; }
    let mut neigh: Vec<u32> = Vec::new();
    for b in &mol.bonds {
        if b.start_atom == atom_id { neigh.push(b.end_atom); }
        else if b.end_atom == atom_id { neigh.push(b.start_atom); }
    }
    if neigh.len() < 3 { return false; }
    neigh.sort_unstable(); neigh.dedup();
    neigh.len() >= 3
}
