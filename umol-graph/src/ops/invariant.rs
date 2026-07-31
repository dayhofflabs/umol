//! Per-atom electron-conservation equation. `check` reports the molecule-wide
//! verdict; `check_atom` evaluates a standalone `AtomAst` (topology-derived
//! valences default to zero); `enumerate_atom` returns the ground `AtomAst`
//! candidates satisfying the equation when fields are `Undetermined`. Shared
//! physics, not tied to a specific valence model.

use std::ops::RangeInclusive;

use thiserror::Error;
use umol_ast::ast::{
    aromatic_covalence, AromaticValenceAst, AsLit, AtomAst, AtomConstraintAst, AtomConstraintsAst,
    AtomId, ElementAst, Lattice, MoleculeAst, MulticenterValenceAst, UnpairedElectronsAst,
    ValueAst,
};
use umol_chem::element::Element;
use umol_chem::spin::{SpinState, UnpairedElectrons};
use umol_utils::solution::Solution;

pub struct ValenceInvariants;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceMismatch {
    #[error("atom {atom_id:?}: orbital count {orbital_count} != electron count {electron_count}")]
    OrbitalCount {
        atom_id: AtomId,
        orbital_count: i64,
        electron_count: i64,
    },
}

impl ValenceInvariants {
    /// Molecule-wide verdict: `Underdetermined` if any atom has a non-`Lit`
    /// field the check can't fire on, `Contradictory` on the first orbital !=
    /// electron mismatch, else `Determined`.
    pub fn check(ast: &MoleculeAst) -> Solution<(), ValenceMismatch> {
        for id in ast.atoms().ids() {
            match Self::check_molecule_atom(ast, id) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => return Solution::Underdetermined(()),
                Solution::Contradictory(mismatch) => return Solution::Contradictory(mismatch),
            }
        }
        Solution::Determined(())
    }

    /// Verdict for a standalone atom. Topology-derived valences default to zero;
    /// only a non-negative literal constraint raises them. `Underdetermined`
    /// when element / charge / implicit-H / lone-pairs / unpaired electrons are not all
    /// `Lit`.
    pub fn check_atom(atom: &AtomAst) -> Solution<(), ValenceMismatch> {
        let (
            Some(element),
            Some(charge),
            Some(implicit_h),
            Some(lone_pairs),
            Some(unpaired_electrons),
        ) = (
            atom.element.as_lit(),
            atom.charge.as_lit(),
            atom.implicit_hydrogens.as_lit(),
            atom.lone_pairs.as_lit(),
            atom.unpaired_electrons.count.as_lit(),
        )
        else {
            return Solution::Underdetermined(());
        };
        let valence = match atom
            .constraints
            .valence()
            .unwrap_or(&ValueAst::Undetermined)
        {
            ValueAst::Lit(v) if *v >= 0 => *v,
            ValueAst::Undetermined => 0,
            _ => return Solution::Underdetermined(()),
        };
        let donated = match atom
            .constraints
            .donated_pairs()
            .unwrap_or(&ValueAst::Undetermined)
        {
            ValueAst::Lit(v) if *v >= 0 => *v,
            ValueAst::Undetermined => 0,
            _ => return Solution::Underdetermined(()),
        };
        let accepted = match atom
            .constraints
            .accepted_pairs()
            .unwrap_or(&ValueAst::Undetermined)
        {
            ValueAst::Lit(v) if *v >= 0 => *v,
            ValueAst::Undetermined => 0,
            _ => return Solution::Underdetermined(()),
        };
        let aromatic_constraint = atom
            .constraints
            .aromatic_valence()
            .unwrap_or(&AromaticValenceAst::Undetermined);
        let aromatic_valence = match aromatic_constraint
            .as_lit()
            .map(|valence| valence.valence_count())
        {
            Some(valence) if valence >= 0 => valence,
            None if aromatic_constraint.is_undetermined() => 0,
            _ => return Solution::Underdetermined(()),
        };
        let multicenter_constraint = atom
            .constraints
            .multicenter_valence()
            .unwrap_or(&MulticenterValenceAst::Undetermined);
        let multicenter_valence = match multicenter_constraint
            .as_lit()
            .map(|valence| valence.valence_count())
        {
            Some(valence) if valence >= 0 => valence,
            None if multicenter_constraint.is_undetermined() => 0,
            _ => return Solution::Underdetermined(()),
        };
        let orbital = orbital_count(
            unpaired_electrons,
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
        if orbital == electron {
            Solution::Determined(())
        } else {
            Solution::Contradictory(ValenceMismatch::OrbitalCount {
                atom_id: AtomId(0),
                orbital_count: orbital,
                electron_count: electron,
            })
        }
    }

    /// Per-atom verdict reading the atom in its molecule context: each valence
    /// is taken from a literal constraint, else the topology-derived value.
    /// `Underdetermined` when any required field is non-`Lit`.
    fn check_molecule_atom(ast: &MoleculeAst, atom_id: AtomId) -> Solution<(), ValenceMismatch> {
        let atom = ast.atom(atom_id);
        let Some(element) = atom.element().as_lit() else {
            return Solution::Underdetermined(());
        };
        let Some(charge) = atom.charge().as_lit() else {
            return Solution::Underdetermined(());
        };
        let Some(implicit_h) = atom.implicit_hydrogens().as_lit() else {
            return Solution::Underdetermined(());
        };
        let Some(lone_pairs) = atom.lone_pairs().as_lit() else {
            return Solution::Underdetermined(());
        };
        let Some(unpaired_electrons) = atom.unpaired_electrons().count.as_lit() else {
            return Solution::Underdetermined(());
        };
        let valence = match atom
            .constraints()
            .valence()
            .unwrap_or(&ValueAst::Undetermined)
        {
            ValueAst::Lit(v) if *v >= 0 => *v,
            ValueAst::Undetermined => match atom.valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Solution::Underdetermined(()),
            },
            _ => return Solution::Underdetermined(()),
        };
        let donated = match atom
            .constraints()
            .donated_pairs()
            .unwrap_or(&ValueAst::Undetermined)
        {
            ValueAst::Lit(v) if *v >= 0 => *v,
            ValueAst::Undetermined => match atom.donated_pairs() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Solution::Underdetermined(()),
            },
            _ => return Solution::Underdetermined(()),
        };
        let accepted = match atom
            .constraints()
            .accepted_pairs()
            .unwrap_or(&ValueAst::Undetermined)
        {
            ValueAst::Lit(v) if *v >= 0 => *v,
            ValueAst::Undetermined => match atom.accepted_pairs() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Solution::Underdetermined(()),
            },
            _ => return Solution::Underdetermined(()),
        };
        let aromatic_constraint = atom
            .constraints()
            .aromatic_valence()
            .unwrap_or(&AromaticValenceAst::Undetermined);
        let aromatic_valence = match aromatic_constraint
            .as_lit()
            .map(|valence| valence.valence_count())
        {
            Some(valence) if valence >= 0 => valence,
            None if aromatic_constraint.is_undetermined() => match atom.aromatic_valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Solution::Underdetermined(()),
            },
            _ => return Solution::Underdetermined(()),
        };
        let multicenter_constraint = atom
            .constraints()
            .multicenter_valence()
            .unwrap_or(&MulticenterValenceAst::Undetermined);
        let multicenter_valence = match multicenter_constraint
            .as_lit()
            .map(|valence| valence.valence_count())
        {
            Some(valence) if valence >= 0 => valence,
            None if multicenter_constraint.is_undetermined() => match atom.multicenter_valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Solution::Underdetermined(()),
            },
            _ => return Solution::Underdetermined(()),
        };
        let orbital = orbital_count(
            unpaired_electrons,
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
        if orbital == electron {
            Solution::Determined(())
        } else {
            Solution::Contradictory(ValenceMismatch::OrbitalCount {
                atom_id,
                orbital_count: orbital,
                electron_count: electron,
            })
        }
    }

    /// Enumerate all ground `AtomAst` candidates for `atom_id` consistent with the per-atom
    /// electron conservation. Returns an empty `Vec` if no candidate exists.
    ///
    /// Per-field bounds:
    /// - `min_charge <= charge <= max_charge`
    /// - `0 <= implicit_hydrogens <= max_implicit_hydrogens`
    /// - `0 <= unpaired_electrons <= max_unpaired_electrons`
    /// - `0 <= lone_pairs <= valence_capacity / 2`
    ///
    /// Structural inputs read from the atom (constraint, else topology):
    /// `valence`, `donated_pairs`, `accepted_pairs`, `aromatic_valence`, `multicenter_valence`.
    ///
    pub fn enumerate_atom(ast: &MoleculeAst, atom_id: AtomId) -> Vec<AtomAst> {
        let atom = ast.atom(atom_id);
        let Some(element) = atom.element().as_lit() else {
            return Vec::new();
        };

        let valence = match atom
            .constraints()
            .valence()
            .unwrap_or(&ValueAst::Undetermined)
        {
            ValueAst::Lit(v) if *v >= 0 => *v,
            ValueAst::Undetermined => match atom.valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };
        let donated = match atom
            .constraints()
            .donated_pairs()
            .unwrap_or(&ValueAst::Undetermined)
        {
            ValueAst::Lit(v) if *v >= 0 => *v,
            ValueAst::Undetermined => match atom.donated_pairs() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };
        let accepted = match atom
            .constraints()
            .accepted_pairs()
            .unwrap_or(&ValueAst::Undetermined)
        {
            ValueAst::Lit(v) if *v >= 0 => *v,
            ValueAst::Undetermined => match atom.accepted_pairs() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };
        let aromatic_constraint = atom
            .constraints()
            .aromatic_valence()
            .unwrap_or(&AromaticValenceAst::Undetermined);
        let aromatic_valence = match aromatic_constraint
            .as_lit()
            .map(|valence| valence.valence_count())
        {
            Some(valence) if valence >= 0 => valence,
            None if aromatic_constraint.is_undetermined() => match atom.aromatic_valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };
        let multicenter_constraint = atom
            .constraints()
            .multicenter_valence()
            .unwrap_or(&MulticenterValenceAst::Undetermined);
        let multicenter_valence = match multicenter_constraint
            .as_lit()
            .map(|valence| valence.valence_count())
        {
            Some(valence) if valence >= 0 => valence,
            None if multicenter_constraint.is_undetermined() => match atom.multicenter_valence() {
                ValueAst::Lit(t) if t >= 0 => t,
                _ => return Vec::new(),
            },
            _ => return Vec::new(),
        };

        let charge_range = enumeration_values(atom.charge(), element_charge_range(element));
        let implicit_h_range = enumeration_values(
            atom.implicit_hydrogens(),
            0..=(element.max_implicit_hydrogens() as i64),
        );
        let lone_pair_range = enumeration_values(
            atom.lone_pairs(),
            0..=((element.valence_capacity() / 2) as i64),
        );
        let unpaired_electron_range = enumeration_values(
            &atom.unpaired_electrons().count,
            0..=(element.max_unpaired_electrons() as i64),
        );

        let mut candidates: Vec<AtomAst> = Vec::new();
        for charge in charge_range {
            for &implicit_h in &implicit_h_range {
                for &lone_pairs in &lone_pair_range {
                    for &unpaired_electrons in &unpaired_electron_range {
                        let orbital = orbital_count(
                            unpaired_electrons,
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
                        let Some(multiplicity) =
                            enumerate_multiplicity(atom.unpaired_electrons(), unpaired_electrons)
                        else {
                            continue;
                        };
                        let assignment = AtomAst {
                            element: ElementAst::Lit(element),
                            charge: ValueAst::Lit(charge),
                            implicit_hydrogens: ValueAst::Lit(implicit_h),
                            lone_pairs: ValueAst::Lit(lone_pairs),
                            unpaired_electrons: UnpairedElectronsAst {
                                count: ValueAst::Lit(unpaired_electrons),
                                multiplicity: ValueAst::Lit(multiplicity),
                            },
                            constraints: AtomConstraintsAst::from(AtomConstraintAst::Valence(
                                ValueAst::Lit(valence),
                            )),
                            ..Default::default()
                        };
                        let Some(candidate) = atom.ast.meet(&assignment) else {
                            continue;
                        };
                        candidates.push(candidate.into_ground());
                    }
                }
            }
        }

        candidates
    }
}

fn enumerate_multiplicity(unpaired_electrons: &UnpairedElectronsAst, count: i64) -> Option<i64> {
    let multiplicity = match unpaired_electrons.multiplicity {
        ValueAst::Lit(m) => m,
        ValueAst::Undetermined => count + 1,
        _ => return None,
    };
    SpinState::try_from(UnpairedElectrons {
        count,
        multiplicity,
    })
    .ok()
    .map(|_| multiplicity)
}

/// Compute number of atomic spin orbitals involved in topology (shared and unshared)
#[allow(clippy::complexity)]
fn orbital_count(
    unpaired_electrons: i64,
    lone_pairs: i64,
    donated_pairs: i64,
    accepted_pairs: i64,
    implicit_h: i64,
    valence: i64,
    aromatic_valence: i64,
    multicenter_valence: i64,
) -> i64 {
    let aromatic_covalence = aromatic_covalence(aromatic_valence);
    unpaired_electrons
        + 2 * lone_pairs
        + 2 * donated_pairs
        + 2 * accepted_pairs
        + 2 * implicit_h
        + 2 * valence
        + aromatic_valence
        + aromatic_covalence
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
    let aromatic_covalence = aromatic_covalence(aromatic_valence);
    (element.valence_electrons() as i64) - charge
        + implicit_h
        + valence
        + aromatic_covalence
        + 2 * accepted_pairs
}

fn element_charge_range(element: Element) -> RangeInclusive<i64> {
    let (lo, hi) = element.charge_bounds();
    (lo as i64)..=(hi as i64)
}

fn enumeration_values(field: &ValueAst, bound: RangeInclusive<i64>) -> Vec<i64> {
    match field {
        ValueAst::Lit(n) => vec![*n],
        ValueAst::Undetermined => bound.collect(),
        _ => bound
            .filter(|value| field.matches(&ValueAst::Lit(*value)))
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AtomAst, AtomConstraintAst, AtomConstraintsAst, AtomId, ElementAst, IsotopeMassAst,
        MoleculeAst, MoleculeParts, UnpairedElectronsAst, ValueAst,
    };
    use umol_chem::element::Element;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::ground_methane(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        Solution::Determined(()),
    )]
    #[case::undetermined_h(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        Solution::Underdetermined(()),
    )]
    #[case::orbital_count_mismatch(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(99),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        Solution::Contradictory(ValenceMismatch::OrbitalCount { atom_id: AtomId(0), orbital_count: 198, electron_count: 103 }),
    )]
    fn test_valence_invariants_check_molecule_atom(
        #[case] ast: MoleculeAst,
        #[case] atom_id: AtomId,
        #[case] expected: Solution<(), ValenceMismatch>,
    ) {
        assert_eq!(ValenceInvariants::check_molecule_atom(&ast, atom_id), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ground_methane(
        AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::not_aromatic(
        AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: AtomConstraintsAst::from(AtomConstraintAst::aromatic_valence(
                AromaticValenceAst::NotAromatic,
            )),
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::aromatic_zero(
        AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: AtomConstraintsAst::from(AtomConstraintAst::aromatic_valence(
                AromaticValenceAst::Aromatic(ValueAst::Lit(0)),
            )),
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::undetermined_charge(
        AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Undetermined,
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        },
        Solution::Underdetermined(()),
    )]
    #[case::orbital_count_mismatch(
        AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(99),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        },
        Solution::Contradictory(ValenceMismatch::OrbitalCount { atom_id: AtomId(0), orbital_count: 198, electron_count: 103 }),
    )]
    fn test_valence_invariants_check_atom(
        #[case] atom: AtomAst,
        #[case] expected: Solution<(), ValenceMismatch>,
    ) {
        assert_eq!(ValenceInvariants::check_atom(&atom), expected);
    }

    #[rstest]
    #[case::ground_methane(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
    )]
    fn test_valence_invariants_check(#[case] ast: MoleculeAst) {
        assert_eq!(ValenceInvariants::check(&ast), Solution::Determined(()));
    }

    #[rstest]
    #[case::methane(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeMassAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            lone_pairs: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: AtomConstraintsAst::from(AtomConstraintAst::Valence(ValueAst::Lit(0))),
        }],
    )]
    #[case::infeasible_h(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(99),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![],
    )]
    // Every enumerated field already ground: the single combination is returned
    // verbatim (oxygen atom, triplet, two lone pairs).
    #[case::all_fields_ground(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::O),
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(2),
            unpaired_electrons: UnpairedElectronsAst::from((2_u8, 3_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomAst {
            element: ElementAst::Lit(Element::O),
            isotope_mass: IsotopeMassAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(2),
            unpaired_electrons: UnpairedElectronsAst::from((2_u8, 3_u8)),
            constraints: AtomConstraintsAst::from(AtomConstraintAst::Valence(ValueAst::Lit(0))),
        }],
    )]
    // Unpaired electrons fully ground to a non-maximal but valid coupling (3 electrons as a
    // doublet); lone pairs open and fixed by conservation to 1.
    #[case::ground_unpaired_electrons(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::N),
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((3_u8, 2_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomAst {
            element: ElementAst::Lit(Element::N),
            isotope_mass: IsotopeMassAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(1),
            unpaired_electrons: UnpairedElectronsAst::from((3_u8, 2_u8)),
            constraints: AtomConstraintsAst::from(AtomConstraintAst::Valence(ValueAst::Lit(0))),
        }],
    )]
    // Multiplicity is ground (singlet), unpaired-electron count open: conservation fixes the
    // count to 2, and the meet keeps the pinned multiplicity (open-shell singlet).
    #[case::ground_multiplicity(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(2),
            lone_pairs: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst {
                count: ValueAst::Undetermined,
                multiplicity: ValueAst::Lit(1),
            },
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeMassAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(2),
            lone_pairs: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((2_u8, 1_u8)),
            constraints: AtomConstraintsAst::from(AtomConstraintAst::Valence(ValueAst::Lit(0))),
        }],
    )]
    // Unpaired electrons pinned to a physically impossible pair (count 2, multiplicity 2):
    // the only conservation-valid count is incompatible, so no candidate.
    #[case::inconsistent_unpaired_electrons(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(2),
            lone_pairs: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![],
    )]
    // Both components open: unpaired-electron count fixed by conservation to 2, multiplicity
    // defaulted by `into_ground` to the maximum (triplet carbene).
    #[case::open_unpaired_electrons(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(2),
            lone_pairs: ValueAst::Lit(0),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeMassAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(2),
            lone_pairs: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((2_u8, 3_u8)),
            constraints: AtomConstraintsAst::from(AtomConstraintAst::Valence(ValueAst::Lit(0))),
        }],
    )]
    // Nonzero (given) charge: oxide anion, 7 electrons, resolves to a doublet.
    #[case::nonzero_charge(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::O),
            charge: ValueAst::Lit(-1),
            implicit_hydrogens: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(3),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomAst {
            element: ElementAst::Lit(Element::O),
            isotope_mass: IsotopeMassAst::Natural,
            charge: ValueAst::Lit(-1),
            implicit_hydrogens: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(3),
            unpaired_electrons: UnpairedElectronsAst::from((1_u8, 2_u8)),
            constraints: AtomConstraintsAst::from(AtomConstraintAst::Valence(ValueAst::Lit(0))),
        }],
    )]
    // Specified isotope survives the meet (Natural can't, Lit(13) does) and is
    // preserved through `into_ground`.
    #[case::specified_isotope(
        MoleculeAst::from_parts(MoleculeParts { atoms: vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeMassAst::Lit(13),
            charge: ValueAst::Lit(0),
            lone_pairs: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeMassAst::Lit(13),
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ValueAst::Lit(4),
            lone_pairs: ValueAst::Lit(0),
            unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: AtomConstraintsAst::from(AtomConstraintAst::Valence(ValueAst::Lit(0))),
        }],
    )]
    fn test_valence_invariants_enumerate_atom(
        #[case] ast: MoleculeAst,
        #[case] atom_id: AtomId,
        #[case] expected: Vec<AtomAst>,
    ) {
        assert_eq!(ValenceInvariants::enumerate_atom(&ast, atom_id), expected);
    }
}
