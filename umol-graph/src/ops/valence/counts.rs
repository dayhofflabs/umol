//! Counts-based valence resolver: per-element valence table with optional
//! implicit-hydrogen inference (RDKit-style).

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AsLit, AtomAst, AtomConstraint, AtomId, AtomView, ElementAst,
    ImplicitHydrogensAst, IsotopeAst, Lattice, MoleculeAst, SpinStateAst, ValueAst,
};
use umol_shared::element::Element;
use umol_shared::spin::{SpinMultiplicity, SpinState, MAX_UNPAIRED_ELECTRONS};

use crate::ops::valence::normal_valence::NormalValenceTable;
use crate::ops::valence::table::ValenceTable;

#[derive(Clone, Debug)]
pub struct CountsValenceResolver {
    pub table: ValenceTable,
    pub normal_valence: NormalValenceTable,
    pub allow_implicit_hydrogens: bool,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum CountsError {
    #[error("no valid valence state for {atom:?} (element {element}, charge {charge:?}, valence {valence:?})")]
    NoValidValenceState {
        atom: AtomId,
        element: Element,
        charge: ValueAst,
        valence: ValueAst,
    },
}

impl CountsValenceResolver {
    pub fn new(
        table: ValenceTable,
        normal_valence: NormalValenceTable,
        allow_implicit_hydrogens: bool,
    ) -> Self {
        Self {
            table,
            normal_valence,
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
            let valence_ast = view.valence();
            let Some(valence) = valence_ast.as_lit().and_then(|n| u8::try_from(n).ok()) else {
                continue;
            };

            let candidates = self.candidates_for(&view, element, valence);
            match candidates.len() {
                0 => {
                    return Err(CountsError::NoValidValenceState {
                        atom: idx,
                        element,
                        charge: view.ast.charge.clone(),
                        valence: valence_ast,
                    });
                }
                1 => {
                    let cand = candidates.into_iter().next().unwrap();
                    let atom_mut = ast.atom_mut(idx).ast;
                    atom_mut.narrow_from(&cand);
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
    ) -> Vec<AtomAst> {
        let atom = view.ast;
        let charge = atom.charge.as_lit_or(0) as i8;
        let entry = match self.table.entry(element) {
            Some(e) => e,
            None => return Vec::new(),
        };

        let aromatic_constraint = atom.constraints.aromatic_valence();
        let is_aromatic = view.is_in_aromatic_system()
            || AromaticValenceAst::aromatic(ValueAst::Undetermined).matches(&aromatic_constraint);
        if is_aromatic {
            let aromatic_valences = if charge != 0 {
                element
                    .shift(-charge)
                    .and_then(|e| self.table.entry(e))
                    .map(|e| e.allowed_aromatic_valences.as_slice())
                    .unwrap_or(entry.allowed_aromatic_valences.as_slice())
            } else {
                entry.allowed_aromatic_valences.as_slice()
            };
            return self.build_aromatic_candidates(
                aromatic_valences,
                atom,
                element,
                charge,
                valence,
                &aromatic_constraint,
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
            .map(|mut ast| {
                ast.constraints
                    .add(AtomConstraint::Valence(ValueAst::Lit(valence as i64)));
                ast
            })
            .collect()
    }

    fn build_aromatic_candidates(
        &self,
        allowed_aromatic_valences: &[u8],
        atom: &AtomAst,
        element: Element,
        charge: i8,
        valence: u8,
        aromatic_constraint: &AromaticValenceAst,
    ) -> Vec<AtomAst> {
        if allowed_aromatic_valences.is_empty() {
            return Vec::new();
        }

        let effective_electrons = (element.valence_electrons() as i16) - (charge as i16);
        let mut candidates = Vec::new();

        for &a in allowed_aromatic_valences {
            let candidate_aromatic = AromaticValenceAst::Aromatic(ValueAst::Lit(a as i64));
            if !aromatic_constraint.matches(&candidate_aromatic) {
                continue;
            }
            let sigma_budget = effective_electrons - (a as i16);
            if sigma_budget < valence as i16 {
                continue;
            }
            let implicit_hydrogens = match &atom.implicit_hydrogens {
                ImplicitHydrogensAst::Lit(n) => *n as u8,
                ImplicitHydrogensAst::Normal => {
                    let Some(h) =
                        self.normal_valence
                            .implicit_hydrogens_for(element, charge, valence, true)
                    else {
                        continue;
                    };
                    h
                }
                ImplicitHydrogensAst::Undetermined => {
                    if self.allow_implicit_hydrogens {
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
            if let Some(mut candidate) =
                try_build_candidate(element, charge, implicit_hydrogens, valence, a, atom)
            {
                candidate
                    .constraints
                    .add(AtomConstraint::Valence(ValueAst::Lit(valence as i64)));
                candidate.constraints.add(AtomConstraint::AromaticValence(
                    AromaticValenceAst::Aromatic(ValueAst::Lit(a as i64)),
                ));
                candidates.push(candidate);
            }
        }

        candidates
    }
}

fn try_build_candidate(
    element: Element,
    charge: i8,
    implicit_hydrogens: u8,
    valence: u8,
    aromatic_pi: u8,
    atom_ast: &AtomAst,
) -> Option<AtomAst> {
    let total_valence = valence + implicit_hydrogens;
    let num_electrons = (element.valence_electrons() as i16) - (charge as i16);
    let unassigned = num_electrons - (total_valence as i16) - (aromatic_pi as i16);
    if unassigned < 0 {
        return None;
    }

    let (unpaired, lone_pairs) = resolve_unpaired_lone_pairs(atom_ast, unassigned)?;
    if unpaired > MAX_UNPAIRED_ELECTRONS {
        return None;
    }

    let spin = if let Some(g) = atom_ast.spin.as_lit() {
        if g.unpaired() != unpaired {
            return None;
        }
        g
    } else if let ValueAst::Lit(m) = &atom_ast.spin.multiplicity {
        let mult = SpinMultiplicity::from_repr(*m as u8)?;
        SpinState::try_new(unpaired, mult).ok()?
    } else {
        SpinState::max_multiplicity(unpaired)?
    };

    Some(AtomAst {
        element: ElementAst::Lit(element),
        isotope_mass: match &atom_ast.isotope_mass {
            IsotopeAst::Undetermined => IsotopeAst::Natural,
            other => other.clone(),
        },
        charge: ValueAst::Lit(charge as i64),
        implicit_hydrogens: ImplicitHydrogensAst::Lit(implicit_hydrogens as i64),
        lone_pairs: ValueAst::Lit(lone_pairs as i64),
        spin: SpinStateAst::from(spin),
        constraints: atom_ast.constraints.clone(),
    })
}

fn resolve_unpaired_lone_pairs(atom_ast: &AtomAst, unassigned: i16) -> Option<(u8, u8)> {
    let fixed_unpaired = match (atom_ast.spin.as_lit(), &atom_ast.spin.unpaired) {
        (Some(s), _) => Some(s.unpaired()),
        (None, ValueAst::Lit(u)) => Some(*u as u8),
        _ => None,
    };

    let fixed_lone_pairs = match &atom_ast.lone_pairs {
        ValueAst::Lit(lp) => Some(*lp as u8),
        _ => None,
    };

    match (fixed_unpaired, fixed_lone_pairs) {
        (None, None) => Some(((unassigned % 2) as u8, (unassigned / 2) as u8)),
        (Some(unpaired), None) => {
            let remaining = unassigned - (unpaired as i16);
            if remaining < 0 || remaining % 2 != 0 {
                return None;
            }
            Some((unpaired, (remaining / 2) as u8))
        }
        (None, Some(lone_pairs)) => {
            let remaining = unassigned - (2 * lone_pairs as i16);
            if remaining < 0 {
                return None;
            }
            Some((remaining as u8, lone_pairs))
        }
        (Some(unpaired), Some(lone_pairs)) => {
            if (unpaired as i16) + (2 * lone_pairs as i16) != unassigned {
                return None;
            }
            Some((unpaired, lone_pairs))
        }
    }
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
        let resolver =
            CountsValenceResolver::new(table, NormalValenceTable::default_table().clone(), true);
        let mut ast = mol_zeroed!(r#"{:atoms ["C #h4"] :bonds []}"#);
        resolver.resolve(&mut ast).unwrap();
    }

    #[rstest]
    fn test_counts_valence_resolver_resolve_methane_implicit_h() {
        let table = ValenceTable::default_table().clone();
        let resolver =
            CountsValenceResolver::new(table, NormalValenceTable::default_table().clone(), true);
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
        let resolver =
            CountsValenceResolver::new(table, NormalValenceTable::default_table().clone(), true);
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
        let table = valence_table! { C => [4] };
        let resolver =
            CountsValenceResolver::new(table, NormalValenceTable::default_table().clone(), true);
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
