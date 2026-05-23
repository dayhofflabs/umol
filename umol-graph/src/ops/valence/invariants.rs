//! Per-atom electron-conservation equation: `check_atom` evaluates on a
//! fully-constrained atom; `solve_atom` returns the ground `AtomAst`
//! candidates satisfying the equation when fields are `Undetermined`.

use std::ops::RangeInclusive;

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AsLit, AtomAst, AtomId, ElementAst, IsotopeAst, MoleculeAst,
    MulticenterValenceAst, SpinStateAst, ValueAst,
};
use umol_shared::element::Element;
use umol_shared::spin::{SpinState, MAX_UNPAIRED_ELECTRONS};

pub struct ValenceInvariants;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum Mismatch {
    #[error("atom {atom:?}: orbital count {orbital_count} != electron count {electron_count}")]
    OrbitalCount {
        atom: AtomId,
        orbital_count: i64,
        electron_count: i64,
    },
}

impl ValenceInvariants {
    /// Returns `Ok(())` when every relevant field is ground and the
    /// orbital == electron equation holds, or when any field is non-`Lit`
    /// (the check cannot fire on a non-ground atom).
    pub fn check_atom(ast: &MoleculeAst, atom: AtomId) -> Result<(), Mismatch> {
        let view = ast.atom(atom);
        let Some(element) = view.element().as_lit() else {
            return Ok(());
        };
        let Some(charge) = view.charge().as_lit() else {
            return Ok(());
        };
        let Some(implicit_h) = view.implicit_hydrogens().as_lit() else {
            return Ok(());
        };
        let Some(lone_pairs) = view.lone_pairs().as_lit() else {
            return Ok(());
        };
        let Some(unpaired) = view.spin().unpaired.as_lit() else {
            return Ok(());
        };
        let valence = match resolve_topology_or_constraint(
            &view.constraints().valence(),
            &view.valence(),
        ) {
            Some(v) => v,
            None => return Ok(()),
        };
        let donated = match resolve_topology_or_constraint(
            &view.constraints().donated_pairs(),
            &view.donated_pairs(),
        ) {
            Some(v) => v,
            None => return Ok(()),
        };
        let accepted = match resolve_topology_or_constraint(
            &view.constraints().accepted_pairs(),
            &view.accepted_pairs(),
        ) {
            Some(v) => v,
            None => return Ok(()),
        };
        let aromatic_valence = match resolve_aromatic(
            &view.constraints().aromatic_valence(),
            &view.aromatic_valence(),
        ) {
            Some(v) => v,
            None => return Ok(()),
        };
        let multicenter_valence = match resolve_multicenter(
            &view.constraints().multicenter_valence(),
            &view.multicenter_valence(),
        ) {
            Some(v) => v,
            None => return Ok(()),
        };

        let orbital = orbital_count(
            unpaired,
            lone_pairs,
            donated,
            accepted,
            implicit_h,
            valence,
            aromatic_valence,
            multicenter_valence,
        );
        let electron = electron_count(element, charge, implicit_h, valence, aromatic_valence, accepted);
        if orbital != electron {
            return Err(Mismatch::OrbitalCount {
                atom,
                orbital_count: orbital,
                electron_count: electron,
            });
        }
        Ok(())
    }

    /// Run `check_atom` for every atom in the molecule. Stops at the first
    /// mismatch.
    pub fn check(ast: &MoleculeAst) -> Result<(), Mismatch> {
        for i in 0..ast.atoms().count() as u32 {
            Self::check_atom(ast, AtomId(i))?;
        }
        Ok(())
    }

    /// Enumerate all ground `AtomAst` candidates for `atom` consistent with
    /// the per-atom electron-conservation equation. Returns an empty `Vec`
    /// if no candidate exists (including the cases where `element` is non-
    /// `Lit` or any non-enumerated field cannot be resolved from constraint
    /// or topology).
    ///
    /// Free-variable enumeration ranges come from element-level bounds:
    /// `charge_bounds`, `max_implicit_hydrogens`, `max_unpaired_electrons`,
    /// and `valence_capacity / 2` for `lone_pairs`. `Lit` fields use their
    /// value directly. `valence`, `donated_pairs`, `accepted_pairs`,
    /// `aromatic_valence`, `multicenter_valence` are read from the view
    /// (constraint or topology); if non-`Lit`, the enumeration yields no
    /// candidates.
    ///
    /// Callers narrow the molecule themselves: this method does not mutate.
    pub fn solve_atom(ast: &MoleculeAst, atom: AtomId) -> Vec<AtomAst> {
        let view = ast.atom(atom);
        let Some(element) = view.element().as_lit() else {
            return Vec::new();
        };

        let Some(valence) = resolve_topology_or_constraint(
            &view.constraints().valence(),
            &view.valence(),
        ) else {
            return Vec::new();
        };
        let Some(donated) = resolve_topology_or_constraint(
            &view.constraints().donated_pairs(),
            &view.donated_pairs(),
        ) else {
            return Vec::new();
        };
        let Some(accepted) = resolve_topology_or_constraint(
            &view.constraints().accepted_pairs(),
            &view.accepted_pairs(),
        ) else {
            return Vec::new();
        };
        let Some(aromatic_valence) = resolve_aromatic(
            &view.constraints().aromatic_valence(),
            &view.aromatic_valence(),
        ) else {
            return Vec::new();
        };
        let Some(multicenter_valence) = resolve_multicenter(
            &view.constraints().multicenter_valence(),
            &view.multicenter_valence(),
        ) else {
            return Vec::new();
        };

        let charges = enumeration_range_i64(view.charge(), element_charge_range(element));
        let implicit_h_range = enumeration_range_i64(
            view.implicit_hydrogens(),
            0..=(element.max_implicit_hydrogens() as i64),
        );
        let lone_pair_range = enumeration_range_i64(
            view.lone_pairs(),
            0..=((element.valence_capacity() / 2) as i64),
        );
        let unpaired_range = enumeration_range_unpaired(&view.spin().unpaired, element);

        let isotope = view.isotope_mass().clone();
        let constraints = view.constraints().clone();
        let original_spin = view.spin().clone();

        let mut candidates: Vec<AtomAst> = Vec::new();
        for charge in charges.iter().copied() {
            for implicit_h in implicit_h_range.iter().copied() {
                for lone_pairs in lone_pair_range.iter().copied() {
                    for unpaired in unpaired_range.iter().copied() {
                        let orbital = orbital_count(
                            unpaired as i64,
                            lone_pairs,
                            donated,
                            accepted,
                            implicit_h,
                            valence,
                            aromatic_valence,
                            multicenter_valence,
                        );
                        let electron = electron_count(
                            element,
                            charge,
                            implicit_h,
                            valence,
                            aromatic_valence,
                            accepted,
                        );
                        if orbital != electron {
                            continue;
                        }
                        let Some(spin) = build_spin(&original_spin, unpaired) else {
                            continue;
                        };
                        candidates.push(AtomAst {
                            element: ElementAst::Lit(element),
                            isotope_mass: match &isotope {
                                IsotopeAst::Undetermined => IsotopeAst::Natural,
                                other => other.clone(),
                            },
                            charge: ValueAst::Lit(charge),
                            implicit_hydrogens: ValueAst::Lit(implicit_h),
                            lone_pairs: ValueAst::Lit(lone_pairs),
                            spin: SpinStateAst::from(spin),
                            constraints: constraints.clone(),
                        });
                    }
                }
            }
        }

        candidates
    }
}

fn orbital_count(
    unpaired: i64,
    lone_pairs: i64,
    donated_pairs: i64,
    accepted_pairs: i64,
    implicit_h: i64,
    valence: i64,
    aromatic_valence: i64,
    multicenter_valence: i64,
) -> i64 {
    let aromatic_increment = if aromatic_valence == 1 { 1 } else { 0 };
    unpaired
        + 2 * lone_pairs
        + 2 * donated_pairs
        + 2 * accepted_pairs
        + 2 * implicit_h
        + 2 * valence
        + aromatic_valence
        + aromatic_increment
        + multicenter_valence
}

fn electron_count(
    element: Element,
    charge: i64,
    implicit_h: i64,
    valence: i64,
    aromatic_valence: i64,
    accepted_pairs: i64,
) -> i64 {
    let aromatic_increment = if aromatic_valence == 1 { 1 } else { 0 };
    (element.valence_electrons() as i64) - charge
        + implicit_h
        + valence
        + aromatic_increment
        + 2 * accepted_pairs
}

/// Resolve a value to a non-negative integer, preferring the constraint
/// store's `Lit` over the topology-derived `Lit`. Returns `None` when
/// neither path produces a non-negative literal.
fn resolve_topology_or_constraint(constraint: &ValueAst, topology: &ValueAst) -> Option<i64> {
    match (constraint, topology) {
        (ValueAst::Lit(v), _) if *v >= 0 => Some(*v),
        (ValueAst::Undetermined, ValueAst::Lit(t)) if *t >= 0 => Some(*t),
        _ => None,
    }
}

/// Resolve an aromatic-valence integer from the constraint container
/// (`Aromatic(Lit)`, `NotAromatic`, or `Undetermined`) falling back to the
/// topology-derived aromatic-system count. `NotAromatic` and absent both
/// resolve to `0`.
fn resolve_aromatic(constraint: &AromaticValenceAst, topology: &ValueAst) -> Option<i64> {
    match (constraint, topology) {
        (AromaticValenceAst::Aromatic(ValueAst::Lit(v)), _) if *v >= 0 => Some(*v),
        (AromaticValenceAst::NotAromatic, _) => Some(0),
        (AromaticValenceAst::Undetermined, ValueAst::Lit(t)) if *t >= 0 => Some(*t),
        _ => None,
    }
}

/// Resolve a multicenter-valence integer from the constraint container
/// (`Multicenter(Lit)`, `NotMulticenter`, or `Undetermined`) falling back
/// to the topology-derived sum.
fn resolve_multicenter(
    constraint: &MulticenterValenceAst,
    topology: &ValueAst,
) -> Option<i64> {
    match (constraint, topology) {
        (MulticenterValenceAst::Multicenter(ValueAst::Lit(v)), _) if *v >= 0 => Some(*v),
        (MulticenterValenceAst::NotMulticenter, _) => Some(0),
        (MulticenterValenceAst::Undetermined, ValueAst::Lit(t)) if *t >= 0 => Some(*t),
        _ => None,
    }
}

fn element_charge_range(element: Element) -> RangeInclusive<i64> {
    let (lo, hi) = element.charge_bounds();
    (lo as i64)..=(hi as i64)
}

fn enumeration_range_i64(
    field: &ValueAst,
    bound: RangeInclusive<i64>,
) -> Vec<i64> {
    match field {
        ValueAst::Lit(n) => vec![*n],
        ValueAst::Undetermined => bound.collect(),
        _ => Vec::new(),
    }
}

fn enumeration_range_unpaired(field: &ValueAst, element: Element) -> Vec<u8> {
    match field {
        ValueAst::Lit(n) => {
            if let Ok(u) = u8::try_from(*n) {
                if u <= MAX_UNPAIRED_ELECTRONS {
                    vec![u]
                } else {
                    Vec::new()
                }
            } else {
                Vec::new()
            }
        }
        ValueAst::Undetermined => (0..=element.max_unpaired_electrons()).collect(),
        _ => Vec::new(),
    }
}

/// Build a `SpinState` consistent with `unpaired` and any multiplicity
/// constraint on the original atom. Returns `None` if no compatible
/// multiplicity exists.
fn build_spin(original: &SpinStateAst, unpaired: u8) -> Option<SpinState> {
    if let Some(g) = original.as_lit() {
        if g.unpaired() == unpaired {
            return Some(g);
        }
        return None;
    }
    if let ValueAst::Lit(m) = &original.multiplicity {
        let mult = umol_shared::spin::SpinMultiplicity::from_repr(*m as u8)?;
        return SpinState::try_new(unpaired, mult).ok();
    }
    SpinState::max_multiplicity(unpaired)
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AtomAst, AtomId, ElementAst, IsotopeAst, MoleculeAst, SpinStateAst, ValueAst,
    };
    use umol_shared::element::Element;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::ground_methane(
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], vec![]),
        AtomId(0),
    )]
    #[case::undetermined_h(
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], vec![]),
        AtomId(0),
    )]
    fn test_valence_invariants_check_atom(#[case] ast: MoleculeAst, #[case] atom: AtomId) {
        assert_eq!(ValenceInvariants::check_atom(&ast, atom), Ok(()));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::orbital_count_mismatch(
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(99),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], vec![]),
        AtomId(0),
        Mismatch::OrbitalCount { atom: AtomId(0), orbital_count: 198, electron_count: 103 },
    )]
    fn test_valence_invariants_check_atom_error(
        #[case] ast: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Mismatch,
    ) {
        assert_eq!(ValenceInvariants::check_atom(&ast, atom), Err(expected));
    }

    #[rstest]
    #[case::ground_methane(
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], vec![]),
    )]
    fn test_valence_invariants_check(#[case] ast: MoleculeAst) {
        assert_eq!(ValenceInvariants::check(&ast), Ok(()));
    }

    #[rstest]
    #[case::methane(
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], vec![]),
        AtomId(0),
        vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            ..Default::default()
        }],
    )]
    #[case::infeasible_h(
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(99),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], vec![]),
        AtomId(0),
        vec![],
    )]
    fn test_valence_invariants_solve_atom(
        #[case] ast: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<AtomAst>,
    ) {
        assert_eq!(ValenceInvariants::solve_atom(&ast, atom), expected);
    }
}
