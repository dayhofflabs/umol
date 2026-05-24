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
use umol_shared::spin::{SpinMultiplicity, SpinState};

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
        let valence = match view.constraints().valence() {
            ValueAst::Lit(v) if v >= 0 => v,
            ValueAst::Undetermined => match view.valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Ok(()),
            },
            _ => return Ok(()),
        };
        let donated = match view.constraints().donated_pairs() {
            ValueAst::Lit(v) if v >= 0 => v,
            ValueAst::Undetermined => match view.donated_pairs() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Ok(()),
            },
            _ => return Ok(()),
        };
        let accepted = match view.constraints().accepted_pairs() {
            ValueAst::Lit(v) if v >= 0 => v,
            ValueAst::Undetermined => match view.accepted_pairs() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Ok(()),
            },
            _ => return Ok(()),
        };
        let aromatic_valence = match view.constraints().aromatic_valence() {
            AromaticValenceAst::Aromatic(ValueAst::Lit(v)) if v >= 0 => v,
            AromaticValenceAst::NotAromatic => 0,
            AromaticValenceAst::Undetermined => match view.aromatic_valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Ok(()),
            },
            _ => return Ok(()),
        };
        let multicenter_valence = match view.constraints().multicenter_valence() {
            MulticenterValenceAst::Multicenter(ValueAst::Lit(v)) if v >= 0 => v,
            MulticenterValenceAst::NotMulticenter => 0,
            MulticenterValenceAst::Undetermined => match view.multicenter_valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Ok(()),
            },
            _ => return Ok(()),
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

    /// Enumerate all ground `AtomAst` candidates for `atom` consistent with the per-atom
    /// electron conservation. Returns an empty `Vec` if no candidate exists.
    ///
    /// Per-field bounds:
    /// - `min_charge <= charge <= max_charge`
    /// - `0 <= implicit_hydrogens <= max_implicit_hydrogens`
    /// - `0 <= unpaired_electrons <= max_unpaired_electrons`
    /// - `0 <= lone_pairs <= valence_capacity / 2`
    ///
    /// Topology-derived constraints from atom view (constraint or topology):
    /// `valence`, `donated_pairs`, `accepted_pairs`, `aromatic_valence`, `multicenter_valence`.
    ///
    pub fn solve_atom(ast: &MoleculeAst, atom: AtomId) -> Vec<AtomAst> {
        let view = ast.atom(atom);
        let Some(element) = view.element().as_lit() else {
            return Vec::new();
        };

        let valence = match view.constraints().valence() {
            ValueAst::Lit(v) if v >= 0 => v,
            ValueAst::Undetermined => match view.valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };
        let donated = match view.constraints().donated_pairs() {
            ValueAst::Lit(v) if v >= 0 => v,
            ValueAst::Undetermined => match view.donated_pairs() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };
        let accepted = match view.constraints().accepted_pairs() {
            ValueAst::Lit(v) if v >= 0 => v,
            ValueAst::Undetermined => match view.accepted_pairs() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };
        let aromatic_valence = match view.constraints().aromatic_valence() {
            AromaticValenceAst::Aromatic(ValueAst::Lit(v)) if v >= 0 => v,
            AromaticValenceAst::NotAromatic => 0,
            AromaticValenceAst::Undetermined => match view.aromatic_valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };
        let multicenter_valence = match view.constraints().multicenter_valence() {
            MulticenterValenceAst::Multicenter(ValueAst::Lit(v)) if v >= 0 => v,
            MulticenterValenceAst::NotMulticenter => 0,
            MulticenterValenceAst::Undetermined => match view.multicenter_valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };

        let charge_range = enumeration_range(view.charge(), element_charge_range(element));
        let implicit_h_range = enumeration_range(
            view.implicit_hydrogens(),
            0..=(element.max_implicit_hydrogens() as i64),
        );
        let lone_pair_range = enumeration_range(
            view.lone_pairs(),
            0..=((element.valence_capacity() / 2) as i64),
        );
        let unpaired_range = enumeration_range(
            &view.spin().unpaired,
            0..=(element.max_unpaired_electrons() as i64),
        );

        let isotope = view.isotope_mass().clone();
        let constraints = view.constraints().clone();
        let original_spin = view.spin().clone();

        let mut candidates: Vec<AtomAst> = Vec::new();
        for charge in charge_range {
            for implicit_h in implicit_h_range.clone() {
                for lone_pairs in lone_pair_range.clone() {
                    for unpaired in unpaired_range.clone() {
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
                        let Ok(unpaired) = u8::try_from(unpaired) else {
                            continue;
                        };
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

/// Compute number of atomic spin orbitals involved in topology (shared and unshared)
#[allow(clippy::complexity)]
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

/// Compute number of electrons assigned to atom
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

fn element_charge_range(element: Element) -> RangeInclusive<i64> {
    let (lo, hi) = element.charge_bounds();
    (lo as i64)..=(hi as i64)
}

/// `Lit(n)` enumerates the single value `n`; `Undetermined` enumerates the
/// element-derived `bound`; any other form enumerates nothing.
fn enumeration_range(field: &ValueAst, bound: RangeInclusive<i64>) -> RangeInclusive<i64> {
    match field {
        ValueAst::Lit(n) => *n..=*n,
        ValueAst::Undetermined => bound,
        _ => RangeInclusive::new(1, 0),
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
        let mult = SpinMultiplicity::from_repr(*m as u8)?;
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
