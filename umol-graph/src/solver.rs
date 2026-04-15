//! Constraint solver: resolution, validation, and matching post-filter.

use smallvec::SmallVec;
use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::spin::{SpinMultiplicity, SpinState, MAX_UNPAIRED_ELECTRONS};
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::AtomIdx;
use crate::ast::atom::AtomAst;
use crate::ast::matcher::Assignment;
use crate::ast::molecule::MoleculeAst;
use crate::atom::{AromaticValence, IsotopeMass};
use crate::graph_ir::atom::Atom;
use crate::graph_ir::config_data::{AtomTypeRegistry, NormalValenceTable, ValenceTable};

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Solution<T> {
    Determined(T),
    Underdetermined(T),
    Contradictory,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Progress {
    Advanced,
    Fixpoint,
    Contradictory,
}

#[derive(Clone, Debug)]
pub enum ValenceStrategy {
    AtomTyping {
        registry: AtomTypeRegistry,
    },
    Counts {
        table: ValenceTable,
        allow_implicit_hydrogens: bool,
    },
}

#[derive(Clone, Debug)]
pub struct Solver {
    pub valence: ValenceStrategy,
}

impl Default for Solver {
    fn default() -> Self {
        Self {
            valence: ValenceStrategy::AtomTyping {
                registry: AtomTypeRegistry::default_registry().clone(),
            },
        }
    }
}

impl Solver {
    pub fn resolve(&self, ast: &mut MoleculeAst) -> Solution<()> {
        loop {
            match self.valence.refine(ast) {
                Progress::Advanced => continue,
                Progress::Fixpoint => break,
                Progress::Contradictory => return Solution::Contradictory,
            }
        }
        if ast.atoms().all(|(_, a)| a.is_ground()) {
            Solution::Determined(())
        } else {
            Solution::Underdetermined(())
        }
    }

    pub fn filter(
        &self,
        _query: &MoleculeAst,
        target: &MoleculeAst,
        assignments: Vec<Assignment>,
    ) -> Vec<Assignment> {
        assignments
            .into_iter()
            .filter(|a| {
                a.0.iter()
                    .all(|&t_idx| self.valence.validate(target, t_idx))
            })
            .collect()
    }

    pub fn validate(&self, ast: &MoleculeAst) -> Solution<()> {
        for i in 0..ast.atom_count() {
            if !self.valence.validate(ast, i) {
                return Solution::Contradictory;
            }
        }
        if ast.atoms().all(|(_, a)| a.is_ground()) {
            Solution::Determined(())
        } else {
            Solution::Underdetermined(())
        }
    }
}

impl ValenceStrategy {
    pub fn candidates_for(&self, ast: &MoleculeAst, idx: AtomIdx) -> SmallVec<[Atom; 4]> {
        let atom = ast.atom(idx);
        let element = match atom.element {
            ElementAst::Lit(e) => e,
            _ => return SmallVec::new(),
        };
        let Some(valence) = ast.bond_order_sum(idx) else {
            return SmallVec::new();
        };
        let charge = atom.charge_or_zero();
        let (donated_pairs, accepted_pairs) = ast.dative_bond_order_sums(idx);
        let is_aromatic = ast.is_in_aromatic_system(idx);

        match self {
            Self::AtomTyping { registry } => atom_typing_candidates(
                registry,
                atom,
                element,
                charge,
                valence,
                is_aromatic,
            ),
            Self::Counts {
                table,
                allow_implicit_hydrogens,
            } => counts_candidates(
                table,
                *allow_implicit_hydrogens,
                atom,
                element,
                charge,
                valence,
                donated_pairs,
                accepted_pairs,
                is_aromatic,
            ),
        }
    }

    pub fn refine(&self, ast: &mut MoleculeAst) -> Progress {
        let mut advanced = false;
        for i in 0..ast.atom_count() as u32 {
            let idx = AtomIdx(i);
            if ast.atom(idx).is_ground() {
                continue;
            }
            if !matches!(ast.atom(idx).element, ElementAst::Lit(_)) {
                continue;
            }
            if ast.bond_order_sum(idx).is_none() {
                continue;
            }
            let candidates = self.candidates_for(ast, idx);
            match candidates.len() {
                0 => return Progress::Contradictory,
                1 => {
                    advanced |= narrow_atom(ast.atom_mut(idx), &candidates[0]);
                }
                _ => {}
            }
        }
        if advanced {
            Progress::Advanced
        } else {
            Progress::Fixpoint
        }
    }

    pub fn validate(&self, ast: &MoleculeAst, atom_index: usize) -> bool {
        let idx = AtomIdx(atom_index as u32);
        if !matches!(ast.atom(idx).element, ElementAst::Lit(_)) {
            return true;
        }
        if ast.bond_order_sum(idx).is_none() {
            return true;
        }
        !self.candidates_for(ast, idx).is_empty()
    }
}

fn atom_typing_candidates(
    registry: &AtomTypeRegistry,
    atom_ast: &AtomAst,
    element: Element,
    charge: i8,
    valence: u8,
    is_aromatic: bool,
) -> SmallVec<[Atom; 4]> {
    let implicit_h_constraint = match &atom_ast.implicit_hydrogens {
        HydrogenAst::Value(ValueAst::Lit(n)) => Some(*n as u8),
        HydrogenAst::Normal => {
            let Some(h) = infer_normal_implicit_hydrogens(element, charge, valence, is_aromatic)
            else {
                return SmallVec::new();
            };
            Some(h)
        }
        HydrogenAst::Undetermined => {
            infer_normal_implicit_hydrogens(element, charge, valence, is_aromatic)
        }
        _ => None,
    };

    let charge_key = match &atom_ast.charge {
        ValueAst::Lit(n) => Some(*n as i8),
        _ => None,
    };

    registry
        .lookup(element, charge_key)
        .iter()
        .filter(|candidate| {
            (match implicit_h_constraint {
                Some(h) => candidate.implicit_hydrogens() == h,
                None => true,
            }) && atom_matches_constraints(candidate, atom_ast)
        })
        .cloned()
        .collect()
}

fn infer_normal_implicit_hydrogens(
    element: Element,
    charge: i8,
    explicit_valence: u8,
    is_aromatic: bool,
) -> Option<u8> {
    if is_aromatic {
        if charge != 0 {
            return None;
        }
        return if element == Element::C {
            Some(3_u8.saturating_sub(explicit_valence))
        } else if matches!(
            element,
            Element::B
                | Element::N
                | Element::O
                | Element::P
                | Element::S
                | Element::Se
                | Element::As
        ) {
            Some(0)
        } else {
            None
        };
    }

    let normal_valence =
        NormalValenceTable::default_table().normal_valence_for(element, charge)?;
    Some(normal_valence.saturating_sub(explicit_valence))
}

#[allow(clippy::too_many_arguments)]
fn counts_candidates(
    table: &ValenceTable,
    allow_implicit_hydrogens: bool,
    atom_ast: &AtomAst,
    element: Element,
    charge: i8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    is_aromatic: bool,
) -> SmallVec<[Atom; 4]> {
    let entry = match table.entry(element) {
        Some(e) => e,
        None => return SmallVec::new(),
    };

    if is_aromatic {
        let aromatic_valences = if charge != 0 {
            element
                .shift(-charge)
                .and_then(|e| table.entry(e))
                .map(|e| e.allowed_aromatic_valences.as_slice())
                .unwrap_or(entry.allowed_aromatic_valences.as_slice())
        } else {
            entry.allowed_aromatic_valences.as_slice()
        };
        return build_aromatic_candidates(
            aromatic_valences,
            atom_ast,
            element,
            charge,
            valence,
            donated_pairs,
            accepted_pairs,
            allow_implicit_hydrogens,
        );
    }

    let implicit_hydrogens = match &atom_ast.implicit_hydrogens {
        HydrogenAst::Value(ValueAst::Lit(n)) => *n as u8,
        _ if allow_implicit_hydrogens => {
            match table.compute_implicit_hydrogens(element, charge, valence) {
                Some(h) => h,
                None => return SmallVec::new(),
            }
        }
        _ => 0,
    };

    try_build_candidate(
        element,
        charge,
        implicit_hydrogens,
        valence,
        donated_pairs,
        accepted_pairs,
        AromaticValence::NotAromatic,
        atom_ast,
    )
    .map(|a| SmallVec::from_elem(a, 1))
    .unwrap_or_default()
}

#[allow(clippy::too_many_arguments)]
fn build_aromatic_candidates(
    allowed_aromatic_valences: &[u8],
    atom_ast: &AtomAst,
    element: Element,
    charge: i8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    allow_implicit_hydrogens: bool,
) -> SmallVec<[Atom; 4]> {
    if allowed_aromatic_valences.is_empty() {
        return SmallVec::new();
    }

    let effective_electrons = (element.valence_electrons() as i16) - (charge as i16);
    let mut candidates = SmallVec::new();

    for &a in allowed_aromatic_valences {
        let sigma_budget = effective_electrons - (a as i16);
        if sigma_budget < valence as i16 {
            continue;
        }
        let implicit_hydrogens = match &atom_ast.implicit_hydrogens {
            HydrogenAst::Value(ValueAst::Lit(n)) => *n as u8,
            HydrogenAst::Normal => {
                let Some(h) =
                    infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
                else {
                    continue;
                };
                h
            }
            HydrogenAst::Undetermined => {
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
        if let Some(atom) = try_build_candidate(
            element,
            charge,
            implicit_hydrogens,
            valence,
            donated_pairs,
            accepted_pairs,
            AromaticValence::Valence(a),
            atom_ast,
        ) {
            candidates.push(atom);
        }
    }

    candidates
}

fn infer_normal_aromatic_implicit_hydrogens(
    element: Element,
    charge: i8,
    valence: u8,
) -> Option<u8> {
    if charge != 0 {
        return None;
    }
    if element == Element::C {
        Some(3 - valence)
    } else if matches!(
        element,
        Element::B | Element::N | Element::O | Element::P | Element::S | Element::Se | Element::As
    ) {
        Some(0)
    } else {
        None
    }
}

#[allow(clippy::too_many_arguments)]
fn try_build_candidate(
    element: Element,
    charge: i8,
    implicit_hydrogens: u8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    aromatic_valence: AromaticValence,
    atom_ast: &AtomAst,
) -> Option<Atom> {
    let total_valence = valence + implicit_hydrogens;
    let num_electrons = (element.valence_electrons() as i16) - (charge as i16);
    let unassigned =
        num_electrons - (total_valence as i16) - (aromatic_valence.valence() as i16);
    if unassigned < 0 {
        return None;
    }

    let (unpaired, lone_pairs) = resolve_unpaired_lone_pairs(atom_ast, unassigned)?;
    if unpaired > MAX_UNPAIRED_ELECTRONS {
        return None;
    }

    let spin = match &atom_ast.spin {
        SpinStateAst::Lit(s) => {
            if s.unpaired_electrons() != unpaired {
                return None;
            }
            *s
        }
        SpinStateAst::Pair {
            multiplicity: ValueAst::Lit(m),
            ..
        } => {
            let mult = SpinMultiplicity::from_multiplicity(*m as u8)?;
            SpinState::try_new(unpaired, mult).ok()?
        }
        _ => SpinState::max_multiplicity(unpaired)?,
    };

    Atom::try_new(
        element,
        None,
        charge,
        implicit_hydrogens,
        lone_pairs,
        unpaired,
        spin.multiplicity(),
        valence,
        donated_pairs,
        accepted_pairs,
        aromatic_valence,
        0,
    )
    .ok()
}

fn resolve_unpaired_lone_pairs(atom_ast: &AtomAst, unassigned: i16) -> Option<(u8, u8)> {
    let fixed_unpaired = match &atom_ast.spin {
        SpinStateAst::Lit(s) => Some(s.unpaired_electrons()),
        SpinStateAst::Pair {
            unpaired: ValueAst::Lit(u),
            ..
        } => Some(*u as u8),
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

fn atom_matches_constraints(candidate: &Atom, atom_ast: &AtomAst) -> bool {
    match_value_constraint(&atom_ast.charge, candidate.charge() as i64)
        && match_value_constraint(&atom_ast.lone_pairs, candidate.lone_pairs() as i64)
        && match_value_constraint(&atom_ast.valence, candidate.valence() as i64)
        && match_value_constraint(
            &atom_ast.donated_pairs,
            candidate.donated_pairs() as i64,
        )
        && match_value_constraint(
            &atom_ast.accepted_pairs,
            candidate.accepted_pairs() as i64,
        )
        && match_value_constraint(
            &atom_ast.multicenter_valence,
            candidate.multicenter_valence() as i64,
        )
        && match_spin_constraint(&atom_ast.spin, candidate.spin())
        && match_aromatic_constraint(&atom_ast.aromatic_valence, candidate.aromatic_valence())
}

fn match_value_constraint(constraint: &ValueAst, value: i64) -> bool {
    constraint.matches(value)
}

fn match_spin_constraint(constraint: &SpinStateAst, value: SpinState) -> bool {
    constraint.matches(value)
}

fn match_aromatic_constraint(constraint: &AromaticValenceAst, value: AromaticValence) -> bool {
    match constraint {
        AromaticValenceAst::Undetermined => true,
        AromaticValenceAst::NotAromatic => value == AromaticValence::NotAromatic,
        AromaticValenceAst::Value(v) => v.matches(value.valence() as i64),
    }
}


fn narrow_atom(atom_ast: &mut AtomAst, candidate: &Atom) -> bool {
    let mut changed = false;
    if matches!(atom_ast.isotope_mass, IsotopeAst::Undetermined) {
        atom_ast.isotope_mass = match candidate.isotope_mass() {
            IsotopeMass::Natural => IsotopeAst::Natural,
            IsotopeMass::MassNumber(n) => IsotopeAst::Lit(n),
        };
        changed = true;
    }
    if matches!(atom_ast.charge, ValueAst::Undetermined) {
        atom_ast.charge = ValueAst::Lit(candidate.charge() as i64);
        changed = true;
    }
    if matches!(atom_ast.implicit_hydrogens, HydrogenAst::Undetermined) {
        atom_ast.implicit_hydrogens = HydrogenAst::Value(ValueAst::Lit(
            candidate.implicit_hydrogens() as i64,
        ));
        changed = true;
    }
    if matches!(atom_ast.lone_pairs, ValueAst::Undetermined) {
        atom_ast.lone_pairs = ValueAst::Lit(candidate.lone_pairs() as i64);
        changed = true;
    }
    if matches!(atom_ast.spin, SpinStateAst::Pair { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Undetermined }) {
        atom_ast.spin = SpinStateAst::Lit(candidate.spin());
        changed = true;
    }
    if matches!(atom_ast.valence, ValueAst::Undetermined) {
        atom_ast.valence = ValueAst::Lit(candidate.valence() as i64);
        changed = true;
    }
    if matches!(atom_ast.donated_pairs, ValueAst::Undetermined) {
        atom_ast.donated_pairs = ValueAst::Lit(candidate.donated_pairs() as i64);
        changed = true;
    }
    if matches!(atom_ast.accepted_pairs, ValueAst::Undetermined) {
        atom_ast.accepted_pairs = ValueAst::Lit(candidate.accepted_pairs() as i64);
        changed = true;
    }
    if matches!(atom_ast.aromatic_valence, AromaticValenceAst::Undetermined) {
        atom_ast.aromatic_valence = match candidate.aromatic_valence() {
            AromaticValence::NotAromatic => AromaticValenceAst::NotAromatic,
            AromaticValence::Valence(v) => AromaticValenceAst::Value(ValueAst::Lit(v as i64)),
        };
        changed = true;
    }
    if matches!(atom_ast.multicenter_valence, ValueAst::Undetermined) {
        atom_ast.multicenter_valence =
            ValueAst::Lit(candidate.multicenter_valence() as i64);
        changed = true;
    }
    changed
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::bond::BondAst;
    use crate::registry;

    use crate::ast::matcher::Assignment;

    fn h2() -> MoleculeAst {
        MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::H),
                AtomAst::from_element(Element::H),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        )
    }

    #[test]
    fn test_solver_resolve_h2_atom_typing() {
        let solver = Solver {
            valence: ValenceStrategy::AtomTyping {
                registry: registry!["H #v"],
            },
        };
        let mut ast = h2();
        let result = solver.resolve(&mut ast);
        assert_eq!(result, Solution::Determined(()));
        assert!(ast.atoms().all(|(_, a)| a.is_ground()));
        assert_eq!(
            ast.atom(AtomIdx(0)).implicit_hydrogens,
            HydrogenAst::Value(ValueAst::Lit(0))
        );
    }

    #[test]
    fn test_solver_resolve_contradictory_empty_registry() {
        let solver = Solver {
            valence: ValenceStrategy::AtomTyping {
                registry: AtomTypeRegistry::new(),
            },
        };
        let mut ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let result = solver.resolve(&mut ast);
        assert_eq!(result, Solution::Contradictory);
    }

    #[test]
    fn test_solver_resolve_already_ground() {
        let solver = Solver::default();
        let mut ast = MoleculeAst::new(vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
        let result = solver.resolve(&mut ast);
        assert_eq!(result, Solution::Determined(()));
    }

    #[test]
    fn test_solver_resolve_wildcard_element_underdetermined() {
        let solver = Solver::default();
        let mut ast = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let result = solver.resolve(&mut ast);
        assert_eq!(result, Solution::Underdetermined(()));
    }

    #[test]
    fn test_solver_resolve_h2_counts() {
        let solver = Solver {
            valence: ValenceStrategy::Counts {
                table: ValenceTable::default_table().clone(),
                allow_implicit_hydrogens: true,
            },
        };
        let mut ast = h2();
        let result = solver.resolve(&mut ast);
        assert_eq!(result, Solution::Determined(()));
        assert!(ast.atoms().all(|(_, a)| a.is_ground()));
        assert_eq!(
            ast.atom(AtomIdx(0)).implicit_hydrogens,
            HydrogenAst::Value(ValueAst::Lit(0))
        );
    }

    #[test]
    fn test_solver_resolve_bare_carbon_counts() {
        let solver = Solver {
            valence: ValenceStrategy::Counts {
                table: ValenceTable::default_table().clone(),
                allow_implicit_hydrogens: true,
            },
        };
        let mut ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let result = solver.resolve(&mut ast);
        assert_eq!(result, Solution::Determined(()));
        assert_eq!(
            ast.atom(AtomIdx(0)).implicit_hydrogens,
            HydrogenAst::Value(ValueAst::Lit(4))
        );
    }

    #[test]
    fn test_solver_validate_ground_h2() {
        let solver = Solver {
            valence: ValenceStrategy::AtomTyping {
                registry: registry!["H #v"],
            },
        };
        let mut ast = h2();
        solver.resolve(&mut ast);
        let result = solver.validate(&ast);
        assert_eq!(result, Solution::Determined(()));
    }

    #[test]
    fn test_solver_validate_non_ground() {
        let solver = Solver::default();
        let ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let result = solver.validate(&ast);
        assert!(matches!(result, Solution::Underdetermined(())));
    }

    #[test]
    fn test_bond_order_sum_no_bonds() {
        let ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(ast.bond_order_sum(AtomIdx(0)), Some(0));
    }

    #[test]
    fn test_bond_order_sum_single() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(ast.bond_order_sum(AtomIdx(0)), Some(1));
    }

    #[test]
    fn test_bond_order_sum_double() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(ast.bond_order_sum(AtomIdx(0)), Some(2));
    }

    #[test]
    fn test_bond_order_sum_wildcard() {
        let ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::new(ValueAst::Undetermined))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(ast.bond_order_sum(AtomIdx(0)), None);
    }

    #[test]
    fn test_solver_filter_all_valid() {
        let solver = Solver {
            valence: ValenceStrategy::AtomTyping {
                registry: registry!["H #v"],
            },
        };
        let mut target = h2();
        solver.resolve(&mut target);
        let assignments = vec![Assignment(vec![0, 1])];
        let result = solver.filter(&MoleculeAst::default(), &target, assignments.clone());
        assert_eq!(result, assignments);
    }

    #[test]
    fn test_solver_filter_empty() {
        let solver = Solver::default();
        let target = MoleculeAst::default();
        let result = solver.filter(&MoleculeAst::default(), &target, vec![]);
        assert_eq!(result, vec![]);
    }
}
