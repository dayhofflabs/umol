//! Counts-based valence resolver: per-element valence table with optional
//! implicit-hydrogen inference (RDKit-style).

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AtomAst, AtomConstraint, AtomId, AtomView, ElementAst,
    ImplicitHydrogensAst, MoleculeAst, ValueAst,
};
use umol_shared::element::Element;

use crate::ops::valence::shared::{
    aromatic_pi_pinned, atom_is_aromatic, charge_or_zero, infer_normal_aromatic_implicit_hydrogens,
    lift_constraints, narrow_atom, try_build_candidate, AtomCandidate,
};
use crate::ops::valence::table::ValenceTable;

#[derive(Clone, Debug)]
pub struct CountsValenceResolver {
    pub table: ValenceTable,
    pub allow_implicit_hydrogens: bool,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum CountsError {
    #[error("no valid valence state for {atom:?} (element {element}, charge {charge}, valence {valence})")]
    NoValidValenceState {
        atom: AtomId,
        element: Element,
        charge: i8,
        valence: u8,
    },
}

impl CountsValenceResolver {
    pub fn new(table: ValenceTable, allow_implicit_hydrogens: bool) -> Self {
        Self {
            table,
            allow_implicit_hydrogens,
        }
    }

    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), CountsError> {
        for i in 0..ast.atoms().count() as u32 {
            let idx = AtomId(i);
            let view = ast.atom(idx);
            if view.ast.is_ground() {
                continue;
            }
            let ElementAst::Lit(element) = view.ast.element else {
                continue;
            };
            let Some(valence) = view.valence().literal().and_then(|n| u8::try_from(n).ok()) else {
                continue;
            };

            let candidates = self.candidates_for(&view, element, valence);
            match candidates.len() {
                0 => {
                    return Err(CountsError::NoValidValenceState {
                        atom: idx,
                        element,
                        charge: charge_or_zero(view.ast),
                        valence,
                    });
                }
                1 => {
                    let cand = candidates.into_iter().next().unwrap();
                    let atom_mut = ast.atom_mut(idx).ast;
                    narrow_atom(atom_mut, &cand.ast);
                    lift_constraints(atom_mut, &cand.lifted);
                }
                _ => {}
            }
        }
        Ok(())
    }

    fn candidates_for(
        &self,
        view: &AtomView<'_>,
        element: Element,
        valence: u8,
    ) -> Vec<AtomCandidate> {
        let atom = view.ast;
        let charge = charge_or_zero(atom);
        let entry = match self.table.entry(element) {
            Some(e) => e,
            None => return Vec::new(),
        };

        if atom_is_aromatic(view) {
            let aromatic_valences = if charge != 0 {
                element
                    .shift(-charge)
                    .and_then(|e| self.table.entry(e))
                    .map(|e| e.allowed_aromatic_valences.as_slice())
                    .unwrap_or(entry.allowed_aromatic_valences.as_slice())
            } else {
                entry.allowed_aromatic_valences.as_slice()
            };
            return build_aromatic_candidates(
                aromatic_valences,
                atom,
                element,
                charge,
                valence,
                self.allow_implicit_hydrogens,
                aromatic_pi_pinned(atom),
            );
        }

        let implicit_hydrogens = match &atom.implicit_hydrogens {
            ImplicitHydrogensAst::Lit(n) => *n as u8,
            _ if self.allow_implicit_hydrogens => {
                match self
                    .table
                    .compute_implicit_hydrogens(element, charge, valence)
                {
                    Some(h) => h,
                    None => return Vec::new(),
                }
            }
            _ => 0,
        };

        try_build_candidate(element, charge, implicit_hydrogens, valence, 0, atom)
            .into_iter()
            .map(|ast| AtomCandidate {
                ast,
                lifted: vec![AtomConstraint::Valence(ValueAst::Lit(valence as i64))],
            })
            .collect()
    }
}

fn build_aromatic_candidates(
    allowed_aromatic_valences: &[u8],
    atom: &AtomAst,
    element: Element,
    charge: i8,
    valence: u8,
    allow_implicit_hydrogens: bool,
    aromatic_pi_pinned: Option<u8>,
) -> Vec<AtomCandidate> {
    if allowed_aromatic_valences.is_empty() {
        return Vec::new();
    }

    let effective_electrons = (element.valence_electrons() as i16) - (charge as i16);
    let mut candidates = Vec::new();

    for &a in allowed_aromatic_valences {
        if let Some(pinned) = aromatic_pi_pinned {
            if a != pinned {
                continue;
            }
        }
        let sigma_budget = effective_electrons - (a as i16);
        if sigma_budget < valence as i16 {
            continue;
        }
        let implicit_hydrogens = match &atom.implicit_hydrogens {
            ImplicitHydrogensAst::Lit(n) => *n as u8,
            ImplicitHydrogensAst::Normal => {
                let Some(h) = infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
                else {
                    continue;
                };
                h
            }
            ImplicitHydrogensAst::Undetermined => {
                if allow_implicit_hydrogens {
                    (sigma_budget - valence as i16) as u8
                } else {
                    0
                }
            }
            _ => continue,
        };
        if implicit_hydrogens > 1 {
            continue;
        }
        let total_sigma = valence + implicit_hydrogens;
        let remaining = effective_electrons - total_sigma as i16 - a as i16;
        if remaining < 0 || remaining % 2 != 0 {
            continue;
        }
        if let Some(candidate) =
            try_build_candidate(element, charge, implicit_hydrogens, valence, a, atom)
        {
            candidates.push(AtomCandidate {
                ast: candidate,
                lifted: vec![
                    AtomConstraint::Valence(ValueAst::Lit(valence as i64)),
                    AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(
                        a as i64,
                    ))),
                ],
            });
        }
    }

    candidates
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{AtomAst, AtomId, BondAst, Constraints, ImplicitHydrogensAst, MoleculeAst};
    use umol_ast::{mol, mol_zeroed};
    use umol_shared::element::Element;

    use super::*;
    use crate::valence_table;

    fn carbon_methane_with_undetermined() -> MoleculeAst {
        // Carbon with implicit_hydrogens left undetermined; valence = 0 (no
        // bonds). Counts resolver should infer 4 implicit Hs.
        mol!(r#"{:atoms ["C"] :bonds []}"#)
    }

    fn ethane() -> MoleculeAst {
        let mut a = AtomAst::from_element(Element::C);
        let mut b = AtomAst::from_element(Element::C);
        a.implicit_hydrogens = ImplicitHydrogensAst::Undetermined;
        b.implicit_hydrogens = ImplicitHydrogensAst::Undetermined;
        MoleculeAst::from_parts(
            vec![a, b],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_counts_valence_resolver_resolve_ground_passthrough() {
        let table = ValenceTable::default_table().clone();
        let resolver = CountsValenceResolver::new(table, true);
        let mut ast = mol_zeroed!(r#"{:atoms ["C #h4"] :bonds []}"#);
        resolver.resolve(&mut ast).unwrap();
    }

    #[rstest]
    fn test_counts_valence_resolver_resolve_methane_implicit_h() {
        let table = ValenceTable::default_table().clone();
        let resolver = CountsValenceResolver::new(table, true);
        let mut ast = carbon_methane_with_undetermined();
        resolver.resolve(&mut ast).unwrap();
        let atom = ast.atom(AtomId(0)).ast;
        assert!(matches!(
            atom.implicit_hydrogens,
            ImplicitHydrogensAst::Lit(4)
        ));
    }

    #[rstest]
    fn test_counts_valence_resolver_resolve_ethane_implicit_h() {
        let table = ValenceTable::default_table().clone();
        let resolver = CountsValenceResolver::new(table, true);
        let mut ast = ethane();
        resolver.resolve(&mut ast).unwrap();
        for i in 0..2 {
            let atom = ast.atom(AtomId(i)).ast;
            assert!(matches!(
                atom.implicit_hydrogens,
                ImplicitHydrogensAst::Lit(3)
            ));
        }
    }

    #[rstest]
    fn test_counts_valence_resolver_resolve_unknown_element() {
        // Custom table without Si entry: a Si atom with valence 0 yields no
        // candidates → contradiction.
        let table = valence_table! { C => [4] };
        let resolver = CountsValenceResolver::new(table, true);
        let si = AtomAst::from_element(Element::Si);
        let mut ast = MoleculeAst::from_atoms_and_bonds(vec![si], vec![]);
        let err = resolver.resolve(&mut ast).unwrap_err();
        assert!(matches!(
            err,
            CountsError::NoValidValenceState {
                element: Element::Si,
                ..
            }
        ));
    }
}
