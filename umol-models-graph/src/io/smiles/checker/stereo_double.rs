use super::super::diagnostics::{Category, Code, Diagnostic, DiagnosticList, Severity, Span};
use crate::io::ir::Molecule;
use crate::io::smiles::checker::{Annotations, StereoAnnotations};
use crate::io::smiles::config::SmilesCheckFlags;
use umol_data::Element;

pub struct StereoArtifacts {
    pub checked_double_bonds: usize,
    pub insufficient_count: usize,
    pub conflict_count: usize,
}

pub fn check_stereo_double(
    mol: &Molecule,
    check_flags: &SmilesCheckFlags,
    annotations: &mut Annotations,
    diagnostics: &mut DiagnosticList,
) {
    if !check_flags.contains(SmilesCheckFlags::STEREO) { return; }
    if annotations.stereo.is_none() { annotations.stereo = Some(StereoAnnotations::default()); }
    let stereo = annotations.stereo.as_mut().unwrap();
    for (bond_id, b) in (1u32..).zip(mol.bonds.iter()) {
        let is_double = matches!(b.symbol, crate::io::ir::BondSymbol::Bond(crate::io::ir::BondOrder::Double));
        if !is_double { continue; }
        if is_stereogenic_candidate(mol, (bond_id - 1) as usize) {
            stereo.candidates.push(bond_id);
        }
        stereo.checked_double_bonds.push(bond_id);
        let _span = match (b.span_start, b.span_end) {
            (Some(s), Some(e)) if e >= s => Span::new(s as usize, e as usize),
            (Some(s), _) => Span::new(s as usize, s as usize + 1),
            _ => Span::new(0, 0),
        };
        // Placeholder: real insufficient/conflict detection to be implemented
        let _ = diagnostics; // silence for now
    }
}

// Minimal stereogenic candidate heuristic:
// A C=C bond is a candidate if each sp2 center has at least two distinct substituents
// (by immediate neighbor identity), ignoring hydrogens and the other alkene carbon.
fn is_stereogenic_candidate(mol: &Molecule, bond_idx: usize) -> bool {
    if let Some(b) = mol.bonds.get(bond_idx) {
        use crate::io::ir::{BondOrder, BondSymbol};
        if !matches!(b.symbol, BondSymbol::Bond(BondOrder::Double)) { return false; }
        let a = (b.start_atom as usize).saturating_sub(1);
        let c = (b.end_atom as usize).saturating_sub(1);
        let n_atoms = mol.atoms.len();
        if a >= n_atoms || c >= n_atoms { return false; }
        let mut neigh_a: Vec<u32> = Vec::new();
        let mut neigh_c: Vec<u32> = Vec::new();
        for bb in &mol.bonds {
            if bb.start_atom == b.start_atom && bb.end_atom != b.end_atom { neigh_a.push(bb.end_atom); }
            if bb.end_atom == b.start_atom && bb.start_atom != b.end_atom { neigh_a.push(bb.start_atom); }
            if bb.start_atom == b.end_atom && bb.end_atom != b.start_atom { neigh_c.push(bb.end_atom); }
            if bb.end_atom == b.end_atom && bb.start_atom != b.start_atom { neigh_c.push(bb.start_atom); }
        }
        // Remove hydrogens from consideration (by element)
        let is_h = |id: u32| -> bool {
            let idx = (id as usize).saturating_sub(1);
            if let Some(atom) = mol.atoms.get(idx) {
                if let crate::io::ir::AtomSymbol::Element(el) = atom.symbol { return el == Element::H; }
            }
            false
        };
        neigh_a.retain(|&id| id != b.end_atom && !is_h(id));
        neigh_c.retain(|&id| id != b.start_atom && !is_h(id));
        // Distinctness by atom id (cheap placeholder; full CIP not implemented here)
        let distinct_a = neigh_a.len() >= 2 && (neigh_a[0] != *neigh_a.get(1).unwrap_or(&neigh_a[0]));
        let distinct_c = neigh_c.len() >= 2 && (neigh_c[0] != *neigh_c.get(1).unwrap_or(&neigh_c[0]));
        return distinct_a && distinct_c;
    }
    false
}
