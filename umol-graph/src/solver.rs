//! Constraint solver: resolution, validation, and matching post-filter.

use smallvec::SmallVec;
use umol_shared::atom_ast::{AromaticValenceAst, ElementAst, HydrogenAst};
use umol_shared::element::Element;
use umol_shared::spin::{SpinMultiplicity, SpinState, MAX_UNPAIRED_ELECTRONS};
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::atom::AtomAst;
use crate::ast::molecule::MoleculeAst;
use crate::atom::AromaticValence;
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
        if ast.atoms.iter().all(|a| !needs_narrowing(a)) {
            Solution::Determined(())
        } else {
            Solution::Underdetermined(())
        }
    }

    pub fn validate(&self, ast: &MoleculeAst) -> Solution<()> {
        for (i, atom) in ast.atoms.iter().enumerate() {
            if !self.valence.validate(atom, ast, i) {
                return Solution::Contradictory;
            }
        }
        if ast.atoms.iter().all(|a| !needs_narrowing(a)) {
            Solution::Determined(())
        } else {
            Solution::Underdetermined(())
        }
    }
}

impl ValenceStrategy {
    pub fn refine(&self, ast: &mut MoleculeAst) -> Progress {
        let mut advanced = false;
        for i in 0..ast.atoms.len() {
            if !needs_narrowing(&ast.atoms[i]) {
                continue;
            }
            let element = match ast.atoms[i].element {
                ElementAst::Lit(e) => e,
                _ => continue,
            };
            let Some(valence) = bond_order_sum(ast, i) else {
                continue;
            };
            let charge = extract_charge(&ast.atoms[i]);
            let (donated_pairs, accepted_pairs) = dative_bond_order_sums(ast, i);
            let is_aromatic = is_in_aromatic_system(ast, i);

            let candidates = match self {
                Self::AtomTyping { registry } => atom_typing_candidates(
                    registry,
                    &ast.atoms[i],
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
                    &ast.atoms[i],
                    element,
                    charge,
                    valence,
                    donated_pairs,
                    accepted_pairs,
                    is_aromatic,
                ),
            };

            match candidates.len() {
                0 => return Progress::Contradictory,
                1 => {
                    narrow_atom(&mut ast.atoms[i], &candidates[0]);
                    advanced = true;
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

    pub fn validate(&self, atom: &AtomAst, ast: &MoleculeAst, atom_index: usize) -> bool {
        let element = match atom.element {
            ElementAst::Lit(e) => e,
            _ => return true, // non-ground element: can't validate, pass through
        };
        let Some(valence) = bond_order_sum(ast, atom_index) else {
            return true;
        };
        let charge = extract_charge(atom);
        let (donated_pairs, accepted_pairs) = dative_bond_order_sums(ast, atom_index);
        let is_aromatic = is_in_aromatic_system(ast, atom_index);

        let candidates = match self {
            Self::AtomTyping { registry } => {
                atom_typing_candidates(registry, atom, element, charge, valence, is_aromatic)
            }
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
        };

        !candidates.is_empty()
    }
}

fn bond_order_sum(ast: &MoleculeAst, atom_index: usize) -> Option<u8> {
    let mut sum: u8 = 0;
    for bond in &ast.bonds {
        if bond.source == atom_index || bond.target == atom_index {
            match bond.bond.order {
                ValueAst::Lit(n) => sum += n as u8,
                _ => return None,
            }
        }
    }
    Some(sum)
}

fn dative_bond_order_sums(ast: &MoleculeAst, atom_index: usize) -> (u8, u8) {
    let mut donated: u8 = 0;
    let mut accepted: u8 = 0;
    for bond in &ast.dative_bonds {
        let order = match bond.bond.order {
            ValueAst::Lit(n) => n as u8,
            _ => continue,
        };
        if bond.source == atom_index {
            donated += order;
        } else if bond.target == atom_index {
            accepted += order;
        }
    }
    (donated, accepted)
}

fn is_in_aromatic_system(ast: &MoleculeAst, atom_index: usize) -> bool {
    ast.aromatic_systems
        .iter()
        .any(|sys| sys.atoms.contains(&atom_index))
}

fn extract_charge(atom: &AtomAst) -> i8 {
    match &atom.charge {
        Some(ValueAst::Lit(n)) => *n as i8,
        _ => 0,
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
    // Determine implicit hydrogen constraint: Some(n) = must match n, None = unconstrained
    let implicit_h_constraint = match &atom_ast.implicit_hydrogens {
        Some(HydrogenAst::Value(ValueAst::Lit(n))) => Some(*n as u8),
        Some(HydrogenAst::Normal) => {
            let Some(h) = infer_normal_implicit_hydrogens(element, charge, valence, is_aromatic)
            else {
                return SmallVec::new();
            };
            Some(h)
        }
        // None = field not specified → unconstrained, accept any candidate
        None => infer_normal_implicit_hydrogens(element, charge, valence, is_aromatic),
        _ => None,
    };

    let charge_key = match &atom_ast.charge {
        Some(ValueAst::Lit(n)) => Some(*n as i8),
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
        Some(HydrogenAst::Value(ValueAst::Lit(n))) => *n as u8,
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
    let has_normal_h = matches!(&atom_ast.implicit_hydrogens, Some(HydrogenAst::Normal));
    let mut candidates = SmallVec::new();

    for &a in allowed_aromatic_valences {
        let sigma_budget = effective_electrons - (a as i16);
        if sigma_budget < valence as i16 {
            continue;
        }
        let implicit_hydrogens = match &atom_ast.implicit_hydrogens {
            Some(HydrogenAst::Value(ValueAst::Lit(n))) => *n as u8,
            Some(HydrogenAst::Normal) => {
                let Some(h) =
                    infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
                else {
                    continue;
                };
                h
            }
            None if has_normal_h => {
                let Some(h) =
                    infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
                else {
                    continue;
                };
                h
            }
            None => {
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
        Some(SpinStateAst::Lit(s)) => {
            if s.unpaired_electrons() != unpaired {
                return None;
            }
            *s
        }
        Some(SpinStateAst::Pair {
            multiplicity: ValueAst::Lit(m),
            ..
        }) => {
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
        Some(SpinStateAst::Lit(s)) => Some(s.unpaired_electrons()),
        Some(SpinStateAst::Pair {
            unpaired: ValueAst::Lit(u),
            ..
        }) => Some(*u as u8),
        _ => None,
    };

    let fixed_lone_pairs = match &atom_ast.lone_pairs {
        Some(ValueAst::Lit(lp)) => Some(*lp as u8),
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

fn match_value_constraint(constraint: &Option<ValueAst>, value: i64) -> bool {
    match constraint {
        None => true,
        Some(v) => v.matches(value),
    }
}

fn match_spin_constraint(constraint: &Option<SpinStateAst>, value: SpinState) -> bool {
    match constraint {
        None => true,
        Some(s) => s.matches(value),
    }
}

fn match_aromatic_constraint(
    constraint: &Option<AromaticValenceAst>,
    value: AromaticValence,
) -> bool {
    match constraint {
        None => true,
        Some(AromaticValenceAst::Undetermined) => true,
        Some(AromaticValenceAst::NotAromatic) => value == AromaticValence::NotAromatic,
        Some(AromaticValenceAst::Value(v)) => v.matches(value.valence() as i64),
    }
}

/// An atom needs narrowing if any valence-relevant field is absent.
/// This differs from `!is_ground()` because `is_ground` treats None as
/// vacuously ground (correct for pattern matching, not for resolution).
fn needs_narrowing(atom: &AtomAst) -> bool {
    atom.charge.is_none()
        || atom.implicit_hydrogens.is_none()
        || atom.lone_pairs.is_none()
        || atom.spin.is_none()
        || atom.valence.is_none()
        || atom.donated_pairs.is_none()
        || atom.accepted_pairs.is_none()
        || atom.aromatic_valence.is_none()
        || atom.multicenter_valence.is_none()
}

fn narrow_atom(atom_ast: &mut AtomAst, candidate: &Atom) {
    if atom_ast.charge.is_none() {
        atom_ast.charge = Some(ValueAst::Lit(candidate.charge() as i64));
    }
    if atom_ast.implicit_hydrogens.is_none() {
        atom_ast.implicit_hydrogens = Some(HydrogenAst::Value(ValueAst::Lit(
            candidate.implicit_hydrogens() as i64,
        )));
    }
    if atom_ast.lone_pairs.is_none() {
        atom_ast.lone_pairs = Some(ValueAst::Lit(candidate.lone_pairs() as i64));
    }
    if atom_ast.spin.is_none() {
        atom_ast.spin = Some(SpinStateAst::Lit(candidate.spin()));
    }
    if atom_ast.valence.is_none() {
        atom_ast.valence = Some(ValueAst::Lit(candidate.valence() as i64));
    }
    if atom_ast.donated_pairs.is_none() {
        atom_ast.donated_pairs = Some(ValueAst::Lit(candidate.donated_pairs() as i64));
    }
    if atom_ast.accepted_pairs.is_none() {
        atom_ast.accepted_pairs = Some(ValueAst::Lit(candidate.accepted_pairs() as i64));
    }
    if atom_ast.aromatic_valence.is_none() {
        atom_ast.aromatic_valence = Some(match candidate.aromatic_valence() {
            AromaticValence::NotAromatic => AromaticValenceAst::NotAromatic,
            AromaticValence::Valence(v) => AromaticValenceAst::Value(ValueAst::Lit(v as i64)),
        });
    }
    if atom_ast.multicenter_valence.is_none() {
        atom_ast.multicenter_valence =
            Some(ValueAst::Lit(candidate.multicenter_valence() as i64));
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::bond::BondAst;
    use crate::ast::molecule::BondTuple;
    use crate::registry;

    fn atom(e: Element) -> AtomAst {
        AtomAst::from_element(e)
    }

    fn bond(source: usize, target: usize, order: u8) -> BondTuple {
        BondTuple {
            source,
            target,
            bond: BondAst::from_order(order),
        }
    }

    fn h2() -> MoleculeAst {
        MoleculeAst {
            atoms: vec![atom(Element::H), atom(Element::H)],
            bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        }
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
        assert!(ast.is_ground());
        assert_eq!(
            ast.atoms[0].implicit_hydrogens,
            Some(HydrogenAst::Value(ValueAst::Lit(0)))
        );
    }

    #[test]
    fn test_solver_resolve_contradictory_empty_registry() {
        let solver = Solver {
            valence: ValenceStrategy::AtomTyping {
                registry: AtomTypeRegistry::new(),
            },
        };
        let mut ast = MoleculeAst {
            atoms: vec![atom(Element::C)],
            ..Default::default()
        };
        let result = solver.resolve(&mut ast);
        assert_eq!(result, Solution::Contradictory);
    }

    #[test]
    fn test_solver_resolve_already_ground() {
        let solver = Solver::default();
        let mut ast = MoleculeAst::default();
        let result = solver.resolve(&mut ast);
        assert_eq!(result, Solution::Determined(()));
    }

    #[test]
    fn test_solver_resolve_wildcard_element_underdetermined() {
        let solver = Solver::default();
        let mut ast = MoleculeAst {
            atoms: vec![AtomAst::new(ElementAst::Undetermined)],
            ..Default::default()
        };
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
        assert!(ast.is_ground());
        assert_eq!(
            ast.atoms[0].implicit_hydrogens,
            Some(HydrogenAst::Value(ValueAst::Lit(0)))
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
        let mut ast = MoleculeAst {
            atoms: vec![atom(Element::C)],
            ..Default::default()
        };
        let result = solver.resolve(&mut ast);
        assert_eq!(result, Solution::Determined(()));
        // C with no bonds → 4 implicit H
        assert_eq!(
            ast.atoms[0].implicit_hydrogens,
            Some(HydrogenAst::Value(ValueAst::Lit(4)))
        );
    }

    #[test]
    fn test_solver_validate_ground_h2() {
        let solver = Solver {
            valence: ValenceStrategy::AtomTyping {
                registry: registry!["H #v"],
            },
        };
        // Fully resolve first, then validate
        let mut ast = h2();
        solver.resolve(&mut ast);
        let result = solver.validate(&ast);
        assert_eq!(result, Solution::Determined(()));
    }

    #[test]
    fn test_solver_validate_non_ground() {
        let solver = Solver::default();
        let ast = MoleculeAst {
            atoms: vec![atom(Element::C)],
            ..Default::default()
        };
        let result = solver.validate(&ast);
        // Non-ground but consistent → Underdetermined
        assert!(matches!(result, Solution::Underdetermined(())));
    }

    #[rstest]
    #[case::no_bonds(
        MoleculeAst { atoms: vec![atom(Element::C)], ..Default::default() },
        0,
        Some(0),
    )]
    #[case::single_bond(
        MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::C)],
            bonds: vec![bond(0, 1, 1)],
            ..Default::default()
        },
        0,
        Some(1),
    )]
    #[case::double_bond_sum(
        MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::O)],
            bonds: vec![bond(0, 1, 2)],
            ..Default::default()
        },
        0,
        Some(2),
    )]
    #[case::non_ground_bond(
        MoleculeAst {
            atoms: vec![atom(Element::C), atom(Element::O)],
            bonds: vec![BondTuple { source: 0, target: 1, bond: BondAst::new(ValueAst::Undetermined) }],
            ..Default::default()
        },
        0,
        None,
    )]
    fn test_bond_order_sum(
        #[case] ast: MoleculeAst,
        #[case] atom_index: usize,
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(bond_order_sum(&ast, atom_index), expected);
    }
}
