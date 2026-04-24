//! Valence theory and candidate enumeration shared by resolver, validator, and matcher.

use umol_ast::ast::atom::{ElementAst, ImplicitHydrogensAst};
use umol_shared::element::Element;
use umol_shared::spin::{SpinMultiplicity, SpinState, MAX_UNPAIRED_ELECTRONS};
use umol_ast::ast::spin::SpinStateAst;
use umol_ast::ast::value::ValueAst;

use crate::ast::AtomIdx;
use crate::ast::atom::AtomAst;
use crate::ast::constraint::{AromaticValenceConstraint, AtomConstraint, MoleculeConstraint};
use crate::ast::molecule::MoleculeAst;
use crate::ops::resolve::Progress;
use crate::ops::valence::{AtomTypeRegistry, NormalValenceTable, ValenceTable};

#[derive(Clone, Debug)]
pub enum ValenceTheory {
    AtomTyping {
        registry: AtomTypeRegistry,
    },
    Counts {
        table: ValenceTable,
        allow_implicit_hydrogens: bool,
    },
}

/// One concrete result for narrowing: the resolved base atom and the per-atom
/// constraints (e.g. `AromaticValence(Lit(n))`) that should be lifted into the
/// molecule's constraint vec to pin the chosen interpretation.
#[derive(Clone, Debug)]
pub struct AtomCandidate {
    pub ast: AtomAst,
    pub lifted: Vec<AtomConstraint>,
}

impl ValenceTheory {
    pub fn candidates_for(&self, ast: &MoleculeAst, idx: AtomIdx) -> Vec<AtomCandidate> {
        let atom = ast.atom(idx).data;
        let element = match atom.element {
            ElementAst::Lit(e) => e,
            _ => return Vec::new(),
        };
        let Some(valence) = ast.valence(idx) else {
            return Vec::new();
        };
        let charge = atom.charge_or_zero();
        let (donated_pairs, accepted_pairs) = ast.dative_bond_order_sums(idx);
        let is_aromatic = ast.atom_is_aromatic(idx);
        let aromatic_pi_pinned = ast.atom_aromatic_valence(idx);

        match self {
            Self::AtomTyping { registry } => atom_typing_candidates(
                registry,
                ast,
                idx,
                atom,
                element,
                charge,
                valence,
                donated_pairs,
                accepted_pairs,
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
                is_aromatic,
                aromatic_pi_pinned,
            ),
        }
    }

    pub fn refine(&self, ast: &mut MoleculeAst) -> Progress {
        let mut advanced = false;
        for i in 0..ast.atoms().count() as u32 {
            let idx = AtomIdx(i);
            if ast.atom(idx).data.is_ground() {
                continue;
            }
            if !matches!(ast.atom(idx).data.element, ElementAst::Lit(_)) {
                continue;
            }
            if ast.valence(idx).is_none() {
                continue;
            }
            let candidates = self.candidates_for(ast, idx);
            match candidates.len() {
                0 => return Progress::Contradictory,
                1 => {
                    let cand = &candidates[0];
                    advanced |= narrow_atom(ast.atom_mut(idx).data, &cand.ast);
                    advanced |= lift_constraints(ast, idx, &cand.lifted);
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

    /// Theory-driven per-atom feasibility check: `true` iff at least one
    /// candidate exists under this theory's matching rules (atom-typing
    /// registry or counts table) for the atom in its current environment.
    /// Not universal; depends on the chosen theory.
    pub fn validate(&self, ast: &MoleculeAst, atom_index: usize) -> bool {
        let idx = AtomIdx(atom_index as u32);
        if !matches!(ast.atom(idx).data.element, ElementAst::Lit(_)) {
            return true;
        }
        if ast.valence(idx).is_none() {
            return true;
        }
        !self.candidates_for(ast, idx).is_empty()
    }
}

/// Per-atom electron-conservation invariant. Theory-independent: orbital-side
/// occupancy must equal source-side electron count. Propagator invoked
/// implicitly on every atom by `Resolver` and `Validator`; not a
/// `MoleculeConstraint` variant. Returns true for non-ground atoms (not yet
/// evaluable).
#[derive(Clone, Copy, Debug, Default)]
pub struct ElectronInvariant;

impl ElectronInvariant {
    pub fn validate(&self, ast: &MoleculeAst, atom_index: usize) -> bool {
        let idx = AtomIdx(atom_index as u32);
        let atom = ast.atom(idx).data;
        let ElementAst::Lit(element) = atom.element else {
            return true;
        };
        let ValueAst::Lit(charge) = atom.charge else {
            return true;
        };
        let ImplicitHydrogensAst::Value(ValueAst::Lit(implicit_h)) = atom.implicit_hydrogens else {
            return true;
        };
        let ValueAst::Lit(lone_pairs) = atom.lone_pairs else {
            return true;
        };
        let SpinStateAst::from_state(spin) = atom.spin else {
            return true;
        };
        let Some(valence) = ast.valence(idx) else {
            return true;
        };
        let (donated_pairs, accepted_pairs) = ast.dative_bond_order_sums(idx);
        let aromatic_valence = ast.atom_aromatic_valence(idx).unwrap_or(0) as i32;
        let aromatic_increment: i32 = if aromatic_valence == 1 { 1 } else { 0 };
        let multicenter_valence = ast.atom_multicenter_valence(idx).unwrap_or(0) as i32;
        let unpaired = spin.unpaired_electrons() as i32;

        let total_orbital = unpaired
            + 2 * lone_pairs as i32
            + 2 * donated_pairs as i32
            + 2 * accepted_pairs as i32
            + 2 * implicit_h as i32
            + 2 * valence as i32
            + aromatic_valence
            + aromatic_increment
            + multicenter_valence;

        let total_source = (element.valence_electrons() as i32) - (charge as i32)
            + implicit_h as i32
            + valence as i32
            + aromatic_increment
            + 2 * accepted_pairs as i32;

        total_orbital == total_source
    }
}

/// Per-entity spin-coupling invariant. Tier-2 physics: a literal
/// `(unpaired, multiplicity)` pair must satisfy `multiplicity = unpaired −
/// 2k + 1` for some `k ∈ 0..=unpaired/2`. Runs on any entity carrying a
/// `SpinStateAst` (atom, aromatic system, multicenter bond). Returns true
/// for non-ground pairs (not yet evaluable).
///
/// Stub: evaluator not implemented. Until wired in, the parser admits
/// physically incompatible literal pairs and the solver will need to reject
/// them here. Migration item tracked in doc 86 §"Invariants" and the
/// Implementation-status list.
#[derive(Clone, Copy, Debug, Default)]
pub struct SpinCouplingInvariant;

impl SpinCouplingInvariant {
    pub fn validate(&self, _spin: &SpinStateAst) -> bool {
        // TODO(doc 86): lift the parity rule from SpinState::are_compatible
        // and wire into Validator::validate + matcher post-filter, matching
        // the ElectronInvariant layout above.
        true
    }
}

#[allow(clippy::too_many_arguments)]
fn atom_typing_candidates(
    registry: &AtomTypeRegistry,
    ast: &MoleculeAst,
    idx: AtomIdx,
    atom_ast: &AtomAst,
    element: Element,
    charge: i8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    is_aromatic: bool,
) -> Vec<AtomCandidate> {
    let implicit_h_constraint = match &atom_ast.implicit_hydrogens {
        ImplicitHydrogensAst::Value(ValueAst::Lit(n)) => Some(*n as u8),
        ImplicitHydrogensAst::Normal => {
            let Some(h) = infer_normal_implicit_hydrogens(element, charge, valence, is_aromatic)
            else {
                return Vec::new();
            };
            Some(h)
        }
        ImplicitHydrogensAst::Undetermined => {
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
        .filter(|pattern| {
            (match implicit_h_constraint {
                Some(h) => match &pattern.ast.implicit_hydrogens {
                    ImplicitHydrogensAst::Value(ValueAst::Lit(n)) => *n as u8 == h,
                    _ => false,
                },
                None => true,
            }) && pattern_constraints_compatible(
                ast,
                idx,
                &pattern.constraints,
                valence,
                donated_pairs,
                accepted_pairs,
            ) && base_atom_compatible(atom_ast, &pattern.ast)
        })
        .map(|pattern| AtomCandidate {
            ast: pattern.ast.clone(),
            lifted: pattern.constraints.clone(),
        })
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
    is_aromatic: bool,
    aromatic_pi_pinned: Option<u8>,
) -> Vec<AtomCandidate> {
    let entry = match table.entry(element) {
        Some(e) => e,
        None => return Vec::new(),
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
            allow_implicit_hydrogens,
            aromatic_pi_pinned,
        );
    }

    let implicit_hydrogens = match &atom_ast.implicit_hydrogens {
        ImplicitHydrogensAst::Value(ValueAst::Lit(n)) => *n as u8,
        _ if allow_implicit_hydrogens => {
            match table.compute_implicit_hydrogens(element, charge, valence) {
                Some(h) => h,
                None => return Vec::new(),
            }
        }
        _ => 0,
    };

    try_build_candidate(element, charge, implicit_hydrogens, valence, 0, atom_ast)
        .into_iter()
        .map(|ast| AtomCandidate {
            ast,
            lifted: vec![AtomConstraint::Valence(ValueAst::Lit(valence as i64))],
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn build_aromatic_candidates(
    allowed_aromatic_valences: &[u8],
    atom_ast: &AtomAst,
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
        let implicit_hydrogens = match &atom_ast.implicit_hydrogens {
            ImplicitHydrogensAst::Value(ValueAst::Lit(n)) => *n as u8,
            ImplicitHydrogensAst::Normal => {
                let Some(h) =
                    infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
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
            try_build_candidate(element, charge, implicit_hydrogens, valence, a, atom_ast)
        {
            candidates.push(AtomCandidate {
                ast: candidate,
                lifted: vec![
                    AtomConstraint::Valence(ValueAst::Lit(valence as i64)),
                    AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(
                        ValueAst::Lit(a as i64),
                    )),
                ],
            });
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

fn base_atom_compatible(query: &AtomAst, candidate: &AtomAst) -> bool {
    value_matches(&query.charge, &candidate.charge)
        && value_matches(&query.lone_pairs, &candidate.lone_pairs)
        && spin_matches(&query.spin, &candidate.spin)
}

fn pattern_constraints_compatible(
    ast: &MoleculeAst,
    idx: AtomIdx,
    constraints: &[AtomConstraint],
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
) -> bool {
    constraints
        .iter()
        .all(|c| atom_constraint_holds(ast, idx, c, valence, donated_pairs, accepted_pairs))
}

fn atom_constraint_holds(
    ast: &MoleculeAst,
    idx: AtomIdx,
    constraint: &AtomConstraint,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
) -> bool {
    match constraint {
        AtomConstraint::Valence(query) => match query {
            ValueAst::Undetermined => true,
            ValueAst::Lit(q) => *q as u8 == valence,
            _ => false,
        },
        AtomConstraint::DonatedPairs(query) => match query {
            ValueAst::Undetermined => true,
            ValueAst::Lit(q) => *q as u8 == donated_pairs,
            _ => false,
        },
        AtomConstraint::AcceptedPairs(query) => match query {
            ValueAst::Undetermined => true,
            ValueAst::Lit(q) => *q as u8 == accepted_pairs,
            _ => false,
        },
        AtomConstraint::AromaticValence(query) => {
            let actual_pi = ast.atom_aromatic_valence(idx);
            let actual_is_aromatic = ast.atom_is_aromatic(idx);
            match query {
                AromaticValenceConstraint::NotAromatic => !actual_is_aromatic,
                AromaticValenceConstraint::Value(ValueAst::Undetermined) => actual_is_aromatic,
                AromaticValenceConstraint::Value(ValueAst::Lit(q)) => match actual_pi {
                    Some(actual) => actual == *q as u8,
                    None => actual_is_aromatic,
                },
                _ => false,
            }
        }
        _ => true,
    }
}

fn value_matches(query: &ValueAst, candidate: &ValueAst) -> bool {
    match (query, candidate) {
        (ValueAst::Undetermined, _) => true,
        (ValueAst::Lit(q), ValueAst::Lit(c)) => q == c,
        _ => false,
    }
}

fn spin_matches(query: &SpinStateAst, candidate: &SpinStateAst) -> bool {
    match (query, candidate) {
        (SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Undetermined }, _) => true,
        (SpinStateAst::from_state(q), SpinStateAst::from_state(c)) => q == c,
        _ => true,
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

    let spin = match &atom_ast.spin {
        SpinStateAst::from_state(s) => {
            if s.unpaired_electrons() != unpaired {
                return None;
            }
            *s
        }
        SpinStateAst {
            multiplicity: ValueAst::Lit(m),
            ..
        } => {
            let mult = SpinMultiplicity::from_multiplicity(*m as u8)?;
            SpinState::try_new(unpaired, mult).ok()?
        }
        _ => SpinState::max_multiplicity(unpaired)?,
    };

    Some(AtomAst {
        element: ElementAst::Lit(element),
        isotope_mass: atom_ast.isotope_mass.clone(),
        charge: ValueAst::Lit(charge as i64),
        implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Lit(implicit_hydrogens as i64)),
        lone_pairs: ValueAst::Lit(lone_pairs as i64),
        spin: SpinStateAst::from_state(spin),
    })
}

fn resolve_unpaired_lone_pairs(atom_ast: &AtomAst, unassigned: i16) -> Option<(u8, u8)> {
    let fixed_unpaired = match &atom_ast.spin {
        SpinStateAst::from_state(s) => Some(s.unpaired_electrons()),
        SpinStateAst {
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

fn lift_constraints(ast: &mut MoleculeAst, idx: AtomIdx, lifted: &[AtomConstraint]) -> bool {
    if lifted.is_empty() {
        return false;
    }
    let mut changed = false;
    for c in lifted {
        if narrow_atom_constraint(ast, idx, c) {
            changed = true;
        }
    }
    changed
}

fn narrow_atom_constraint(
    ast: &mut MoleculeAst,
    idx: AtomIdx,
    new_constraint: &AtomConstraint,
) -> bool {
    let kind = new_constraint.kind();
    let existing = ast.constraints().atoms().get(&idx).and_then(|s| s.get(kind));
    match existing {
        Some(e) if !narrowable(e, new_constraint) => false,
        _ => {
            ast.constraints_mut().insert(MoleculeConstraint::AtomPred(
                idx,
                new_constraint.clone(),
            ));
            true
        }
    }
}

fn narrowable(existing: &AtomConstraint, new_c: &AtomConstraint) -> bool {
    use AromaticValenceConstraint as A;
    use AtomConstraint as C;
    matches!(
        (existing, new_c),
        (
            C::AromaticValence(A::Value(ValueAst::Undetermined)),
            C::AromaticValence(A::Value(ValueAst::Lit(_)))
        ) | (
            C::Valence(ValueAst::Undetermined),
            C::Valence(ValueAst::Lit(_))
        ) | (
            C::DonatedPairs(ValueAst::Undetermined),
            C::DonatedPairs(ValueAst::Lit(_))
        ) | (
            C::AcceptedPairs(ValueAst::Undetermined),
            C::AcceptedPairs(ValueAst::Lit(_))
        ) | (
            C::MulticenterValence(ValueAst::Undetermined),
            C::MulticenterValence(ValueAst::Lit(_))
        )
    )
}

fn narrow_atom(atom_ast: &mut AtomAst, candidate: &AtomAst) -> bool {
    let mut changed = false;
    changed |= narrow_value(&mut atom_ast.charge, &candidate.charge);
    if matches!(
        atom_ast.implicit_hydrogens,
        ImplicitHydrogensAst::Undetermined | ImplicitHydrogensAst::Normal
    ) && atom_ast.implicit_hydrogens != candidate.implicit_hydrogens
    {
        atom_ast.implicit_hydrogens = candidate.implicit_hydrogens.clone();
        changed = true;
    }
    changed |= narrow_value(&mut atom_ast.lone_pairs, &candidate.lone_pairs);
    if !matches!(candidate.spin, SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Undetermined })
        && matches!(atom_ast.spin, SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Undetermined })
    {
        atom_ast.spin = candidate.spin.clone();
        changed = true;
    }
    changed
}

fn narrow_value(target: &mut ValueAst, source: &ValueAst) -> bool {
    if matches!(target, ValueAst::Undetermined) && matches!(source, ValueAst::Lit(_)) {
        *target = source.clone();
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use umol_shared::element::Element;

    use super::*;
    use umol_ast::ast::atom::IsotopeAst;

    use crate::ast::bond::BondAst;
    use crate::ast::config::AtomAstConfig;
    use crate::registry;
    use crate::ops::aromaticity::AromaticityTheory;
    use crate::ops::chemistry::Chemistry;
    use crate::ops::resolve::Resolver;
    use crate::ops::solution::Solution;
    use crate::ops::validate::Validator;

    fn coerce_zeroed(ast: &mut MoleculeAst) {
        let cfg = AtomAstConfig::zeroed();
        for i in 0..ast.atoms().count() as u32 {
            ast.atom_mut(AtomIdx(i)).data.coerce(&cfg);
        }
    }

    fn ground_bond(order: i64) -> BondAst {
        BondAst {
            order: ValueAst::Lit(order),
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from_state(SpinState::closed_shell()),
        }
    }

    fn h2() -> MoleculeAst {
        let mut ast = MoleculeAst::new(
            vec![
                AtomAst::from_element(Element::H),
                AtomAst::from_element(Element::H),
            ],
            vec![(AtomIdx(0), AtomIdx(1), ground_bond(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        coerce_zeroed(&mut ast);
        ast
    }

    fn determined_ast(
        solver: &Chemistry,
        ast: MoleculeAst,
    ) -> MoleculeAst {
        match Resolver::new(solver).resolve(ast).unwrap() {
            Solution::Determined(a) => a,
            other => panic!("expected Determined, got {:?}", other),
        }
    }

    #[test]
    fn test_solver_resolve_h2_atom_typing() {
        let solver = Chemistry {
            valence: ValenceTheory::AtomTyping {
                registry: registry!["H #v"],
            },
            aromaticity: AromaticityTheory::daylight(),
        };
        let ast = determined_ast(&solver, h2());
        assert_eq!(
            ast.atom(AtomIdx(0)).data.implicit_hydrogens,
            ImplicitHydrogensAst::Value(ValueAst::Lit(0))
        );
    }

    #[test]
    fn test_solver_resolve_contradictory_empty_registry() {
        let solver = Chemistry {
            valence: ValenceTheory::AtomTyping {
                registry: AtomTypeRegistry::new(),
            },
            aromaticity: AromaticityTheory::daylight(),
        };
        let ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let outcome = Resolver::new(&solver).resolve(ast).unwrap();
        assert!(
            matches!(outcome, Solution::Contradictory | Solution::Underdetermined(_)),
            "got {:?}",
            outcome
        );
    }

    #[test]
    fn test_solver_resolve_already_ground() {
        let solver = Chemistry::default();
        let ast = MoleculeAst::new(vec![], vec![], vec![], vec![], vec![], vec![], vec![]);
        assert!(matches!(
            Resolver::new(&solver).resolve(ast).unwrap(),
            Solution::Determined(_)
        ));
    }

    #[test]
    fn test_solver_resolve_wildcard_element_underdetermined() {
        let solver = Chemistry::default();
        let ast = MoleculeAst::new(
            vec![AtomAst::new(ElementAst::Undetermined)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert!(matches!(
            Resolver::new(&solver).resolve(ast).unwrap(),
            Solution::Underdetermined(_)
        ));
    }

    #[test]
    fn test_solver_resolve_h2_counts() {
        let solver = Chemistry {
            valence: ValenceTheory::Counts {
                table: ValenceTable::default_table().clone(),
                allow_implicit_hydrogens: true,
            },
            aromaticity: AromaticityTheory::daylight(),
        };
        let ast = determined_ast(&solver, h2());
        assert_eq!(
            ast.atom(AtomIdx(0)).data.implicit_hydrogens,
            ImplicitHydrogensAst::Value(ValueAst::Lit(0))
        );
    }

    #[test]
    fn test_solver_resolve_bare_carbon_counts() {
        let solver = Chemistry {
            valence: ValenceTheory::Counts {
                table: ValenceTable::default_table().clone(),
                allow_implicit_hydrogens: true,
            },
            aromaticity: AromaticityTheory::daylight(),
        };
        let ast = MoleculeAst::new(
            vec![AtomAst {
                element: ElementAst::Lit(Element::C),
                isotope_mass: IsotopeAst::Natural,
                ..Default::default()
            }],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let ast = determined_ast(&solver, ast);
        assert_eq!(
            ast.atom(AtomIdx(0)).data.implicit_hydrogens,
            ImplicitHydrogensAst::Value(ValueAst::Lit(4))
        );
    }

    #[test]
    fn test_solver_validate_ground_h2() {
        let solver = Chemistry {
            valence: ValenceTheory::AtomTyping {
                registry: registry!["H #v"],
            },
            aromaticity: AromaticityTheory::daylight(),
        };
        let ast = determined_ast(&solver, h2());
        assert!(Validator::new(&solver).validate(&ast).is_determined());
    }

    #[test]
    fn test_solver_validate_non_ground() {
        let solver = Chemistry::default();
        let ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert!(matches!(
            Validator::new(&solver).validate(&ast),
            Solution::Underdetermined(())
        ));
    }

    #[test]
    fn test_valence_no_bonds() {
        let ast = MoleculeAst::new(
            vec![AtomAst::from_element(Element::C)],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert_eq!(ast.valence(AtomIdx(0)), Some(0));
    }

    #[test]
    fn test_valence_single() {
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
        assert_eq!(ast.valence(AtomIdx(0)), Some(1));
    }

    #[test]
    fn test_valence_double() {
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
        assert_eq!(ast.valence(AtomIdx(0)), Some(2));
    }

    #[test]
    fn test_valence_wildcard() {
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
        assert_eq!(ast.valence(AtomIdx(0)), None);
    }

    #[test]
    fn test_electron_invariant_multicenter_term() {
        let atom = AtomAst {
            element: ElementAst::Lit(Element::H),
            isotope_mass: IsotopeAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Lit(0)),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::from_state(SpinState::closed_shell()),
        };
        let mut ast = MoleculeAst::new(
            vec![atom],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        assert!(!ElectronInvariant.validate(&ast, 0));
        ast.constraints_mut().insert(MoleculeConstraint::AtomPred(
            AtomIdx(0),
            AtomConstraint::MulticenterValence(ValueAst::Lit(1)),
        ));
        assert!(ElectronInvariant.validate(&ast, 0));
    }
}
