//! Temporary experimental doc-99 counts resolver, run in parallel with
//! [`CountsValenceResolver`] to compare invariants-based resolution against
//! the current `compute_implicit_hydrogens` approach on the conformance
//! corpus. The `covalence_set` bound is computed under a selectable
//! [`ValenceScheme`] so the charge/unpaired adjustment ordering can be tested
//! against examples. Delete once the scheme is settled.
//!
//! In this branch `#a! #m! #d0 #t0` hold, so `total_valence` collapses to
//! `v + h` (topology bond valence + implicit hydrogens). Charge and unpaired
//! signs are pinned to old counts' isoelectronic behavior: `O⁻` (charge −1,
//! 0 unpaired) maps to `F`'s `covalence_set = [1]`, i.e. a cation raises the
//! allowed valence (`+charge`) and an unpaired electron lowers it
//! (`−unpaired`).

use umol_ast::ast::{
    AromaticValenceAst, AsLit, AtomAst, AtomConstraint, AtomId, AtomView, Lattice, MoleculeAst,
    ValueAst,
};
use umol_shared::element::Element;

use crate::ops::valence::{ValenceEntry, ValenceInvariants, ValenceTable};

/// How the per-element `covalence_set` bound on `v + h` is adjusted for charge
/// and unpaired electrons.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValenceScheme {
    /// Scheme 1: adjust the looked-up value arithmetically —
    /// `v + h` is allowed iff `(v + h) − charge + unpaired ∈ covalence_set(element)`.
    ArithmeticAdjust,
    /// Scheme 2: shift the element, then look up —
    /// `v + h` is allowed iff `(v + h) ∈ covalence_set(element.shift(−charge − unpaired))`.
    ElementShift,
}

#[derive(Clone, Debug)]
pub struct CountsNewResolver {
    table: ValenceTable,
    scheme: ValenceScheme,
}

impl CountsNewResolver {
    pub fn new(table: ValenceTable, scheme: ValenceScheme) -> Self {
        Self { table, scheme }
    }

    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), String> {
        for i in 0..ast.atoms().count() as u32 {
            self.resolve_atom(ast, AtomId(i))?;
        }
        Ok(())
    }

    fn resolve_atom(&self, ast: &mut MoleculeAst, atom_id: AtomId) -> Result<(), String> {
        let view = ast.atom(atom_id);
        if view.ast.is_ground() {
            return Ok(());
        }
        let Some(element) = view.element().as_lit() else {
            return Ok(());
        };
        let entry = self.table.entry(element);
        let aromatic_trials = aromatic_trials(&view, entry);
        let topology_v = view.valence();
        let saved = view.ast.clone();

        let mut candidates: Vec<AtomAst> = Vec::new();
        for av in &aromatic_trials {
            {
                let m = ast.atom_mut(atom_id).ast;
                *m = saved.clone();
                m.constraints
                    .add(AtomConstraint::AromaticValence(av.clone()));
                if let ValueAst::Lit(v) = topology_v {
                    m.constraints.add(AtomConstraint::Valence(ValueAst::Lit(v)));
                }
            }
            candidates.extend(ValenceInvariants::solve_atom(ast, atom_id));
        }
        {
            let m = ast.atom_mut(atom_id).ast;
            *m = saved;
        }

        if let Some(entry) = entry {
            if !entry.covalence_set.is_empty() {
                let topology_v = topology_v.as_lit().unwrap_or(0);
                candidates.retain(|c| {
                    matches!(
                        c.constraints.aromatic_valence(),
                        AromaticValenceAst::Aromatic(_)
                    ) || self.valence_allowed(element, &entry.covalence_set, topology_v, c)
                });
            }
        }

        if let Some(min_u) = candidates.iter().filter_map(|c| c.spin.unpaired.as_lit()).min() {
            candidates.retain(|c| c.spin.unpaired.as_lit() == Some(min_u));
        }
        if let Some(max_n) = candidates.iter().filter_map(|c| c.lone_pairs.as_lit()).max() {
            candidates.retain(|c| c.lone_pairs.as_lit() == Some(max_n));
        }

        match candidates.len() {
            0 => Err(format!("atom {atom_id:?}: no candidate ({element})")),
            1 => {
                let cand = candidates.into_iter().next().unwrap();
                ast.atom_mut(atom_id).ast.narrow_from(&cand);
                Ok(())
            }
            n => Err(format!("atom {atom_id:?}: {n} ambiguous candidates ({element})")),
        }
    }

    fn valence_allowed(
        &self,
        element: Element,
        covalence_set: &[u8],
        topology_v: i64,
        candidate: &AtomAst,
    ) -> bool {
        let Some(h) = candidate.implicit_hydrogens.as_lit() else {
            return false;
        };
        let charge = candidate.charge.as_lit().unwrap_or(0);
        let unpaired = candidate.spin.unpaired.as_lit().unwrap_or(0);
        let total_valence = topology_v + h;
        match self.scheme {
            ValenceScheme::ArithmeticAdjust => u8::try_from(total_valence - charge + unpaired)
                .is_ok_and(|t| covalence_set.contains(&t)),
            ValenceScheme::ElementShift => {
                let Ok(delta) = i8::try_from(-charge - unpaired) else {
                    return false;
                };
                let Some(shifted) = element.shift(delta) else {
                    return false;
                };
                let Some(shifted_entry) = self.table.entry(shifted) else {
                    return false;
                };
                u8::try_from(total_valence).is_ok_and(|t| shifted_entry.covalence_set.contains(&t))
            }
        }
    }
}

fn aromatic_trials(view: &AtomView<'_>, entry: Option<&ValenceEntry>) -> Vec<AromaticValenceAst> {
    let aromatic_constraint = view.constraints().aromatic_valence();
    let is_aromatic = view.is_in_aromatic_system()
        || AromaticValenceAst::aromatic(ValueAst::Undetermined).matches(&aromatic_constraint);
    if !is_aromatic {
        return vec![AromaticValenceAst::NotAromatic];
    }
    let set = entry
        .map(|e| e.aromatic_valence_set.as_slice())
        .unwrap_or(&[]);
    if set.is_empty() {
        return vec![
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
            AromaticValenceAst::Aromatic(ValueAst::Lit(2)),
        ];
    }
    set.iter()
        .map(|&v| AromaticValenceAst::Aromatic(ValueAst::Lit(v as i64)))
        .collect()
}
