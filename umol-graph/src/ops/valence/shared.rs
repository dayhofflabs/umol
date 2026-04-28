//! Shared helpers for the AtomTyping and Counts valence resolvers.

use umol_ast::ast::{
    AromaticValenceAst, AtomAst, AtomConstraint, AtomConstraintKind, AtomIdx, AtomView, BondAst,
    ElementAst, ImplicitHydrogensAst, MoleculeAst, SpinStateAst, ValueAst,
};
use umol_shared::element::Element;
use umol_shared::spin::{SpinMultiplicity, SpinState, MAX_UNPAIRED_ELECTRONS};

use crate::ops::valence::table::ValenceTable;

/// One concrete narrowing target: the resolved base atom and the per-atom
/// constraints to lift onto `AtomAst.constraints` to pin the chosen
/// interpretation.
#[derive(Clone, Debug)]
pub struct AtomCandidate {
    pub ast: AtomAst,
    pub lifted: Vec<AtomConstraint>,
}

pub fn charge_or_zero(atom: &AtomAst) -> i8 {
    match atom.charge {
        ValueAst::Lit(c) => c as i8,
        _ => 0,
    }
}

/// Extracts the per-atom aromatic π-electron count from the inline
/// `AromaticValence` constraint, if pinned to a literal.
pub fn aromatic_pi_pinned(atom: &AtomAst) -> Option<u8> {
    match atom.constraints.get(AtomConstraintKind::AromaticValence)? {
        AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(n)))
            if *n >= 0 =>
        {
            Some(*n as u8)
        }
        _ => None,
    }
}

/// True iff the atom is in any aromatic system (relation membership) or
/// carries an `Aromatic(_)` constraint declaration.
pub fn atom_is_aromatic(view: &AtomView<'_>) -> bool {
    if view.is_in_aromatic_system() {
        return true;
    }
    matches!(
        view.data.constraints.get(AtomConstraintKind::AromaticValence),
        Some(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(_)))
    )
}

/// Ground SpinState if both fields are concrete, else `None`.
pub fn ground_spin_state(spin: &SpinStateAst) -> Option<SpinState> {
    let (ValueAst::Lit(u), ValueAst::Lit(m)) = (&spin.unpaired, &spin.multiplicity) else {
        return None;
    };
    let mult = SpinMultiplicity::from_multiplicity(*m as u8)?;
    SpinState::try_new(*u as u8, mult).ok()
}

pub fn spin_state_undetermined(spin: &SpinStateAst) -> bool {
    matches!(spin.unpaired, ValueAst::Undetermined)
        && matches!(spin.multiplicity, ValueAst::Undetermined)
}

pub fn infer_normal_implicit_hydrogens(
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

    let normal_valence = ValenceTable::default_table().normal_valence_for(element, charge)?;
    Some(normal_valence.saturating_sub(explicit_valence))
}

pub fn infer_normal_aromatic_implicit_hydrogens(
    element: Element,
    charge: i8,
    valence: u8,
) -> Option<u8> {
    if charge != 0 {
        return None;
    }
    if element == Element::C {
        Some(3_u8.saturating_sub(valence))
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
    }
}

pub fn base_atom_compatible(query: &AtomAst, candidate: &AtomAst) -> bool {
    value_matches(&query.charge, &candidate.charge)
        && value_matches(&query.lone_pairs, &candidate.lone_pairs)
        && spin_matches(&query.spin, &candidate.spin)
}

pub fn value_matches(query: &ValueAst, candidate: &ValueAst) -> bool {
    match (query, candidate) {
        (ValueAst::Undetermined, _) => true,
        (ValueAst::Lit(q), ValueAst::Lit(c)) => q == c,
        _ => false,
    }
}

pub fn spin_matches(query: &SpinStateAst, candidate: &SpinStateAst) -> bool {
    if spin_state_undetermined(query) {
        return true;
    }
    match (ground_spin_state(query), ground_spin_state(candidate)) {
        (Some(q), Some(c)) => q == c,
        _ => true,
    }
}

pub fn pattern_constraints_compatible(
    view: &AtomView<'_>,
    constraints: &[AtomConstraint],
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
) -> bool {
    constraints
        .iter()
        .all(|c| atom_constraint_holds(view, c, valence, donated_pairs, accepted_pairs))
}

fn atom_constraint_holds(
    view: &AtomView<'_>,
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
            let actual_pi = aromatic_pi_pinned(view.data);
            let actual_is_aromatic = atom_is_aromatic(view);
            match query {
                AromaticValenceAst::NotAromatic => !actual_is_aromatic,
                AromaticValenceAst::Aromatic(ValueAst::Undetermined) => actual_is_aromatic,
                AromaticValenceAst::Aromatic(ValueAst::Lit(q)) => match actual_pi {
                    Some(actual) => actual == *q as u8,
                    None => actual_is_aromatic,
                },
                _ => false,
            }
        }
        _ => true,
    }
}

pub fn try_build_candidate(
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

    let spin = if let Some(g) = ground_spin_state(&atom_ast.spin) {
        if g.unpaired() != unpaired {
            return None;
        }
        g
    } else if let ValueAst::Lit(m) = &atom_ast.spin.multiplicity {
        let mult = SpinMultiplicity::from_multiplicity(*m as u8)?;
        SpinState::try_new(unpaired, mult).ok()?
    } else {
        SpinState::max_multiplicity(unpaired)?
    };

    Some(AtomAst {
        element: ElementAst::Lit(element),
        isotope_mass: atom_ast.isotope_mass.clone(),
        charge: ValueAst::Lit(charge as i64),
        implicit_hydrogens: ImplicitHydrogensAst::Lit(implicit_hydrogens as i64),
        lone_pairs: ValueAst::Lit(lone_pairs as i64),
        spin: SpinStateAst::from_state(spin),
        constraints: atom_ast.constraints.clone(),
    })
}

fn resolve_unpaired_lone_pairs(atom_ast: &AtomAst, unassigned: i16) -> Option<(u8, u8)> {
    let fixed_unpaired = match (
        ground_spin_state(&atom_ast.spin),
        &atom_ast.spin.unpaired,
    ) {
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

/// Apply candidate fields to the atom in place. Only narrows `Undetermined`
/// fields onto literal candidates; never overwrites existing literals.
pub fn narrow_atom(atom: &mut AtomAst, candidate: &AtomAst) -> bool {
    let mut changed = false;
    changed |= narrow_value(&mut atom.charge, &candidate.charge);
    if matches!(
        atom.implicit_hydrogens,
        ImplicitHydrogensAst::Undetermined | ImplicitHydrogensAst::Normal
    ) && atom.implicit_hydrogens != candidate.implicit_hydrogens
    {
        atom.implicit_hydrogens = candidate.implicit_hydrogens.clone();
        changed = true;
    }
    changed |= narrow_value(&mut atom.lone_pairs, &candidate.lone_pairs);
    if !spin_state_undetermined(&candidate.spin) && spin_state_undetermined(&atom.spin) {
        atom.spin = candidate.spin.clone();
        changed = true;
    }
    changed
}

pub fn narrow_value(target: &mut ValueAst, source: &ValueAst) -> bool {
    if matches!(target, ValueAst::Undetermined) && matches!(source, ValueAst::Lit(_)) {
        *target = source.clone();
        true
    } else {
        false
    }
}

/// Adds each lifted constraint onto the atom, narrowing where possible.
/// Constraints are stored inline on `AtomAst.constraints`; "lifted" is a
/// legacy term — there is no separate molecule-level container.
pub fn lift_constraints(atom: &mut AtomAst, lifted: &[AtomConstraint]) -> bool {
    let mut changed = false;
    for c in lifted {
        if narrow_atom_constraint(atom, c) {
            changed = true;
        }
    }
    changed
}

fn narrow_atom_constraint(atom: &mut AtomAst, new_c: &AtomConstraint) -> bool {
    let kind = new_c.kind();
    let existing = atom.constraints.get(kind);
    match existing {
        Some(e) if !narrowable(e, new_c) => false,
        _ => {
            atom.constraints.add(new_c.clone());
            true
        }
    }
}

fn narrowable(existing: &AtomConstraint, new_c: &AtomConstraint) -> bool {
    use AromaticValenceAst as A;
    use AtomConstraint as C;
    matches!(
        (existing, new_c),
        (
            C::AromaticValence(A::Aromatic(ValueAst::Undetermined)),
            C::AromaticValence(A::Aromatic(ValueAst::Lit(_)))
        ) | (
            C::Valence(ValueAst::Undetermined),
            C::Valence(ValueAst::Lit(_))
        ) | (
            C::DonatedPairs(ValueAst::Undetermined),
            C::DonatedPairs(ValueAst::Lit(_))
        ) | (
            C::AcceptedPairs(ValueAst::Undetermined),
            C::AcceptedPairs(ValueAst::Lit(_))
        )
    )
}

/// σ-bond order sum, restricted to the value the legacy resolver used — the
/// AtomView aggregate over neighbors. Returns `None` for non-ground bond
/// orders.
pub fn atom_sigma_valence(view: &AtomView<'_>) -> Option<u8> {
    let v = view.bond_order_sum()?;
    u8::try_from(v).ok()
}

/// Donor/acceptor pair counts on incident dative bonds. Either component is
/// `None` if any contributing dative's `order` is non-ground (or if the donor
/// side aggregates over multi-donor datives, which have no per-atom share).
pub fn atom_dative_counts(view: &AtomView<'_>) -> (Option<u8>, Option<u8>) {
    (
        view.donated_pairs().and_then(|v| u8::try_from(v).ok()),
        view.accepted_pairs().and_then(|v| u8::try_from(v).ok()),
    )
}

/// Bond order field on `BondAst` (for callers that already hold the bond).
#[allow(dead_code)]
pub fn bond_order_lit(bond: &BondAst) -> Option<u8> {
    match &bond.order {
        ValueAst::Lit(n) if *n >= 0 => u8::try_from(*n).ok(),
        _ => None,
    }
}

/// Convenience: fetch `(idx, &AtomAst)` views with the molecule back-ref.
#[allow(dead_code)]
pub fn atom_view<'a>(ast: &'a MoleculeAst, idx: AtomIdx) -> AtomView<'a> {
    ast.atom(idx)
}
