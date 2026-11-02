//! Per-atom electron-conservation equation. `check` reports the molecule-wide
//! solution; `check_atom` evaluates a bare `AtomForm` (topology-derived
//! valences default to zero); `enumerate_atom` returns the ground `AtomForm`
//! candidates satisfying the equation when fields are `Undetermined`. Shared
//! physics, not tied to a specific valence model.

use std::ops::RangeInclusive;

use thiserror::Error;
use umol_chem::element::Element;
use umol_chem::spin::{SpinState, UnpairedElectrons};
use umol_graph_ir::ir::{
    aromatic_covalence, AsLit, AtomConstraintForm, AtomConstraintKey, AtomConstraintsForm,
    AtomForm, AtomId, ElementForm, Lattice, Molecule, NumForm, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

/// Low-level electron-conservation checks and candidate enumeration.
///
/// Validation callers should normally use [`ValenceInvariantsValidator`]. Candidate enumeration
/// is exposed separately because it returns possible ground atom forms rather than a validation
/// solution.
pub struct ValenceInvariants;

/// Validates the model-independent electron-conservation invariant.
#[derive(Clone, Copy, Debug, Default)]
pub struct ValenceInvariantsValidator;

/// Setup failures from valence-invariant validation.
///
/// The current implementation has no setup failure, but the distinct error channel preserves the
/// common validator result contract.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceInvariantsError {}

/// A model-independent electron-conservation contradiction.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValenceMismatch {
    #[error("atom {atom_id:?}: contradictory {key:?} reading")]
    ContradictoryConstraint {
        atom_id: AtomId,
        key: AtomConstraintKey,
    },
    #[error("atom {atom_id:?}: orbital count {orbital_count} != electron count {electron_count}")]
    OrbitalCount {
        atom_id: AtomId,
        orbital_count: i64,
        electron_count: i64,
    },
}

impl ValenceInvariants {
    /// Number of electrons assigned to an atom, including shared electrons.
    ///
    /// Uses the actual element's valence electrons and the supplied contributions.
    /// This is a signed calculation, not a validity or capacity check.
    #[inline]
    pub fn electron_count(
        element: Element,
        charge: i64,
        implicit_hydrogens: i64,
        valence: i64,
        aromatic_valence: i64,
        accepted_pairs: i64,
    ) -> i64 {
        i64::from(element.valence_electrons()) - charge
            + implicit_hydrogens
            + valence
            + aromatic_covalence(aromatic_valence)
            + 2 * accepted_pairs
    }

    /// Electrons available for lone pairs and unpaired electrons (`2*n + u`).
    ///
    /// Subtracts bonding occupation from the electron count. Accepted pairs
    /// cancel from that balance. A negative result is retained, indicating
    /// that the supplied bonding contributions exceed the available electrons.
    #[inline]
    pub fn nonbonding_electrons(
        element: Element,
        charge: i64,
        implicit_hydrogens: i64,
        valence: i64,
        aromatic_valence: i64,
        donated_pairs: i64,
        multicenter_valence: i64,
    ) -> i64 {
        Self::electron_count(
            element,
            charge,
            implicit_hydrogens,
            valence,
            aromatic_valence,
            0,
        ) - orbital_count(
            0,
            0,
            donated_pairs,
            0,
            implicit_hydrogens,
            valence,
            aromatic_valence,
            multicenter_valence,
        )
    }

    /// Molecule-wide solution: `Underdetermined` if any atom has a non-`Lit`
    /// field the check can't fire on, `Contradictory` on the first orbital !=
    /// electron mismatch, else `Determined`.
    fn check(molecule: &Molecule) -> Solution<(), ValenceMismatch> {
        for id in molecule.atoms().ids() {
            match Self::check_molecule_atom(molecule, id) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => return Solution::Underdetermined(()),
                Solution::Contradictory(mismatch) => return Solution::Contradictory(mismatch),
            }
        }
        Solution::Determined(())
    }

    /// Solution for a bare `AtomForm`. Topology-derived valences default to zero;
    /// only a non-negative literal constraint raises them. `Underdetermined`
    /// when element / charge / implicit-H / lone-pairs / unpaired electrons are not all
    /// `Lit`.
    fn check_atom(atom: &AtomForm) -> Solution<(), ValenceMismatch> {
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
        let valence = match atom.constraints.valence() {
            NumForm::Lit(v) if *v >= 0 => *v,
            NumForm::Undetermined => 0,
            _ => return Solution::Underdetermined(()),
        };
        let donated = match atom.constraints.donated_pairs() {
            NumForm::Lit(v) if *v >= 0 => *v,
            NumForm::Undetermined => 0,
            _ => return Solution::Underdetermined(()),
        };
        let accepted = match atom.constraints.accepted_pairs() {
            NumForm::Lit(v) if *v >= 0 => *v,
            NumForm::Undetermined => 0,
            _ => return Solution::Underdetermined(()),
        };
        let aromatic_constraint = atom.constraints.aromatic_valence();
        let aromatic_valence = match aromatic_constraint
            .as_lit()
            .map(|valence| valence.valence_count())
        {
            Some(valence) if valence >= 0 => valence,
            None if aromatic_constraint.is_undetermined() => 0,
            _ => return Solution::Underdetermined(()),
        };
        let multicenter_constraint = atom.constraints.multicenter_valence();
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
        let electron = Self::electron_count(
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

    /// Per-atom solution using complete evidence for each contribution.
    /// `Underdetermined` when any required field is non-`Lit`.
    fn check_molecule_atom(molecule: &Molecule, atom_id: AtomId) -> Solution<(), ValenceMismatch> {
        let atom = molecule.atom(atom_id);
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
            .evidenced_complete(AtomConstraintKey::Valence)
        {
            Some(AtomConstraintForm::Valence(form)) => form.as_lit(),
            None => {
                return Solution::Contradictory(ValenceMismatch::ContradictoryConstraint {
                    atom_id,
                    key: AtomConstraintKey::Valence,
                })
            }
            _ => unreachable!(),
        };
        let Some(valence) = valence.filter(|value| *value >= 0) else {
            return Solution::Underdetermined(());
        };
        let donated = match atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::DonatedPairs)
        {
            Some(AtomConstraintForm::DonatedPairs(form)) => form.as_lit(),
            None => {
                return Solution::Contradictory(ValenceMismatch::ContradictoryConstraint {
                    atom_id,
                    key: AtomConstraintKey::DonatedPairs,
                })
            }
            _ => unreachable!(),
        };
        let Some(donated) = donated.filter(|value| *value >= 0) else {
            return Solution::Underdetermined(());
        };
        let accepted = match atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::AcceptedPairs)
        {
            Some(AtomConstraintForm::AcceptedPairs(form)) => form.as_lit(),
            None => {
                return Solution::Contradictory(ValenceMismatch::ContradictoryConstraint {
                    atom_id,
                    key: AtomConstraintKey::AcceptedPairs,
                })
            }
            _ => unreachable!(),
        };
        let Some(accepted) = accepted.filter(|value| *value >= 0) else {
            return Solution::Underdetermined(());
        };
        let aromatic_valence = match atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::AromaticValence)
        {
            Some(AtomConstraintForm::AromaticValence(form)) => {
                form.as_lit().map(|value| value.valence_count())
            }
            None => {
                return Solution::Contradictory(ValenceMismatch::ContradictoryConstraint {
                    atom_id,
                    key: AtomConstraintKey::AromaticValence,
                })
            }
            _ => unreachable!(),
        };
        let Some(aromatic_valence) = aromatic_valence.filter(|value| *value >= 0) else {
            return Solution::Underdetermined(());
        };
        let multicenter_valence = match atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::MulticenterValence)
        {
            Some(AtomConstraintForm::MulticenterValence(form)) => {
                form.as_lit().map(|value| value.valence_count())
            }
            None => {
                return Solution::Contradictory(ValenceMismatch::ContradictoryConstraint {
                    atom_id,
                    key: AtomConstraintKey::MulticenterValence,
                })
            }
            _ => unreachable!(),
        };
        let Some(multicenter_valence) = multicenter_valence.filter(|value| *value >= 0) else {
            return Solution::Underdetermined(());
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
        let electron = Self::electron_count(
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

    /// Enumerate all ground `AtomForm` candidates for `atom_id` consistent with the per-atom
    /// electron conservation. Returns an empty `Vec` if no candidate exists.
    ///
    /// Per-field bounds:
    /// - `min_charge <= charge <= max_charge`
    /// - `0 <= implicit_hydrogens <= max_implicit_hydrogens`
    /// - `0 <= unpaired_electrons <= max_unpaired_electrons`
    /// - `0 <= lone_pairs <= valence_capacity / 2`
    ///
    /// Structural inputs use complete evidence, meeting assertions with relations:
    /// `valence`, `donated_pairs`, `accepted_pairs`, `aromatic_valence`, `multicenter_valence`.
    ///
    pub fn enumerate_atom(molecule: &Molecule, atom_id: AtomId) -> Vec<AtomForm> {
        let atom = molecule.atom(atom_id);
        let Some(element) = atom.element().as_lit() else {
            return Vec::new();
        };

        let valence = match atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::Valence)
        {
            Some(AtomConstraintForm::Valence(form)) => form.as_lit(),
            None => return Vec::new(),
            _ => unreachable!(),
        };
        let Some(valence) = valence.filter(|value| *value >= 0) else {
            return Vec::new();
        };
        let donated = match atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::DonatedPairs)
        {
            Some(AtomConstraintForm::DonatedPairs(form)) => form.as_lit(),
            None => return Vec::new(),
            _ => unreachable!(),
        };
        let Some(donated) = donated.filter(|value| *value >= 0) else {
            return Vec::new();
        };
        let accepted = match atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::AcceptedPairs)
        {
            Some(AtomConstraintForm::AcceptedPairs(form)) => form.as_lit(),
            None => return Vec::new(),
            _ => unreachable!(),
        };
        let Some(accepted) = accepted.filter(|value| *value >= 0) else {
            return Vec::new();
        };
        let aromatic_valence = match atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::AromaticValence)
        {
            Some(AtomConstraintForm::AromaticValence(form)) => {
                form.as_lit().map(|value| value.valence_count())
            }
            None => return Vec::new(),
            _ => unreachable!(),
        };
        let Some(aromatic_valence) = aromatic_valence.filter(|value| *value >= 0) else {
            return Vec::new();
        };
        let multicenter_valence = match atom
            .constraints()
            .evidenced_complete(AtomConstraintKey::MulticenterValence)
        {
            Some(AtomConstraintForm::MulticenterValence(form)) => {
                form.as_lit().map(|value| value.valence_count())
            }
            None => return Vec::new(),
            _ => unreachable!(),
        };
        let Some(multicenter_valence) = multicenter_valence.filter(|value| *value >= 0) else {
            return Vec::new();
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

        let mut candidates: Vec<AtomForm> = Vec::new();
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
                        let electron = Self::electron_count(
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
                        let assignment = AtomForm {
                            element: ElementForm::Lit(element),
                            charge: NumForm::Lit(charge),
                            implicit_hydrogens: NumForm::Lit(implicit_h),
                            lone_pairs: NumForm::Lit(lone_pairs),
                            unpaired_electrons: UnpairedElectronsForm {
                                count: NumForm::Lit(unpaired_electrons),
                                multiplicity: NumForm::Lit(multiplicity),
                            },
                            constraints: AtomConstraintsForm::from(AtomConstraintForm::Valence(
                                NumForm::Lit(valence),
                            )),
                            ..Default::default()
                        };
                        let Some(candidate) = atom.attributes().meet(&assignment) else {
                            continue;
                        };
                        candidates.push(candidate.into_concrete());
                    }
                }
            }
        }

        candidates
    }
}

impl ValenceInvariantsValidator {
    /// Validate electron conservation for every molecule atom.
    ///
    /// Non-literal required values produce [`Solution::Underdetermined`]; a literal orbital/electron
    /// mismatch produces [`Solution::Contradictory`].
    pub fn validate(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<(), ValenceMismatch>, ValenceInvariantsError> {
        Ok(ValenceInvariants::check(molecule))
    }

    /// Validate electron conservation for a free-standing atom with zero topology contribution.
    pub fn validate_atom(
        &self,
        atom: &AtomForm,
    ) -> Result<Solution<(), ValenceMismatch>, ValenceInvariantsError> {
        Ok(ValenceInvariants::check_atom(atom))
    }
}

fn enumerate_multiplicity(unpaired_electrons: &UnpairedElectronsForm, count: i64) -> Option<i64> {
    let multiplicity = match unpaired_electrons.multiplicity {
        NumForm::Lit(m) => m,
        NumForm::Undetermined => count + 1,
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

fn element_charge_range(element: Element) -> RangeInclusive<i64> {
    let (lo, hi) = element.charge_bounds();
    (lo as i64)..=(hi as i64)
}

fn enumeration_values(field: &NumForm, bound: RangeInclusive<i64>) -> Vec<i64> {
    match field {
        NumForm::Lit(n) => vec![*n],
        NumForm::Undetermined => bound.collect(),
        _ => bound
            .filter(|value| field.matches(&NumForm::Lit(*value)))
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{
        AromaticValenceForm, AtomConstraintForm, AtomConstraintsForm, AtomForm, AtomId, BondForm,
        ElementForm, IsotopeMassForm, Molecule, MoleculeEntries, NumForm, UnpairedElectronsForm,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl};

    use super::*;

    #[rstest]
    #[case::methane(Element::C, 0, 4, 0, 0, 0, 8)]
    #[case::aromatic_carbon(Element::C, 0, 1, 2, 1, 0, 8)]
    #[case::aromatic_oxygen(Element::O, 0, 0, 2, 2, 0, 8)]
    #[case::borane_acceptor(Element::B, 0, 3, 0, 0, 1, 8)]
    #[case::hydride(Element::H, -1, 0, 0, 0, 0, 2)]
    fn test_valence_invariants_electron_count(
        #[case] element: Element,
        #[case] charge: i64,
        #[case] implicit_hydrogens: i64,
        #[case] valence: i64,
        #[case] aromatic_valence: i64,
        #[case] accepted_pairs: i64,
        #[case] expected: i64,
    ) {
        assert_eq!(
            ValenceInvariants::electron_count(
                element,
                charge,
                implicit_hydrogens,
                valence,
                aromatic_valence,
                accepted_pairs,
            ),
            expected,
        );
    }

    #[rstest]
    #[case::methane((Element::C, 0), 4, 0, 0, 0, 0, 0)]
    #[case::ammonia((Element::N, 0), 3, 0, 0, 0, 0, 2)]
    #[case::ammonia_donor((Element::N, 0), 3, 0, 0, 1, 0, 0)]
    #[case::bridging_hydride((Element::H, -1), 0, 0, 0, 0, 2, 0)]
    #[case::xenon_difluoride((Element::Xe, 0), 0, 0, 0, 0, 2, 6)]
    #[case::aromatic_oxygen((Element::O, 0), 0, 2, 2, 0, 0, 2)]
    #[case::combined_contributions((Element::O, 0), 0, 1, 1, 1, 1, 1)]
    #[case::negative_budget((Element::H, 0), 0, 1, 0, 0, 1, -1)]
    fn test_valence_invariants_nonbonding_electrons(
        #[case] atom: (Element, i64),
        #[case] implicit_hydrogens: i64,
        #[case] valence: i64,
        #[case] aromatic_valence: i64,
        #[case] donated_pairs: i64,
        #[case] multicenter_valence: i64,
        #[case] expected: i64,
    ) {
        let (element, charge) = atom;
        assert_eq!(
            ValenceInvariants::nonbonding_electrons(
                element,
                charge,
                implicit_hydrogens,
                valence,
                aromatic_valence,
                donated_pairs,
                multicenter_valence,
            ),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ground_methane(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(4),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        Solution::Determined(()),
    )]
    #[case::undetermined_h(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        Solution::Underdetermined(()),
    )]
    #[case::orbital_count_mismatch(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(99),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        Solution::Contradictory(ValenceMismatch::OrbitalCount { atom_id: AtomId(0), orbital_count: 198, electron_count: 103 }),
    )]
    #[case::contradictory_valence(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm {
                element: ElementForm::Lit(Element::C), charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(3), lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
                ..Default::default()
            }; 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::new(NumForm::lit_set([])))],
            ..Default::default()
        }),
        AtomId(0),
        Solution::Contradictory(ValenceMismatch::ContradictoryConstraint { atom_id: AtomId(0), key: AtomConstraintKey::Valence }),
    )]
    fn test_valence_invariants_check_molecule_atom(
        #[case] molecule: Molecule,
        #[case] atom_id: AtomId,
        #[case] expected: Solution<(), ValenceMismatch>,
    ) {
        assert_eq!(ValenceInvariants::check_molecule_atom(&molecule, atom_id), expected);
    }

    #[rstest]
    #[case::donor(
        r#"{:atoms ["N#c0#h3#n0#u0#s" "B"]
            :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#,
        Solution::Determined(())
    )]
    #[case::unknown_donor(
        r#"{:atoms ["N#c0#h3#n0#u0#s" "B"]
            :dative-bonds [{:donors [0] :acceptor 1 :attrs "*"}]}"#,
        Solution::Underdetermined(())
    )]
    #[case::asserted_donor(
        r#"{:atoms ["N#c0#h3#n0#u0#s#d" "B"]
            :dative-bonds [{:donors [0] :acceptor 1 :attrs "*"}]}"#,
        Solution::Determined(())
    )]
    #[case::conflicting_donor(
        r#"{:atoms ["N#c0#h3#n1#u0#s#d0" "B"]
            :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#,
        Solution::Contradictory(ValenceMismatch::ContradictoryConstraint {
            atom_id: AtomId(0), key: AtomConstraintKey::DonatedPairs,
        })
    )]
    #[case::multicenter(
        r#"{:atoms ["H#c-#h0#n0#u0#s" "Cr" "Cr"]
            :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
        Solution::Determined(())
    )]
    #[case::unknown_multicenter(
        r#"{:atoms ["H#c-#h0#n0#u0#s" "Cr" "Cr"]
            :multicenter-bonds [{:atoms [0 1 2] :attrs "*"}]}"#,
        Solution::Underdetermined(())
    )]
    #[case::asserted_multicenter(
        r#"{:atoms ["H#c-#h0#n0#u0#s#m2" "Cr" "Cr"]
            :multicenter-bonds [{:atoms [0 1 2] :attrs "*"}]}"#,
        Solution::Determined(())
    )]
    #[case::conflicting_multicenter(
        r#"{:atoms ["H#c-#h0#n1#u0#s#m0" "Cr" "Cr"]
            :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
        Solution::Contradictory(ValenceMismatch::ContradictoryConstraint {
            atom_id: AtomId(0), key: AtomConstraintKey::MulticenterValence,
        })
    )]
    #[case::conflicting_aromatic(
        r#"{:atoms ["C#c0#h1#n0#u0#s#a1" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"]]
            :aromatic-systems [{:atoms [0 1 2] :attrs "[0,1,1]"}]}"#,
        Solution::Contradictory(ValenceMismatch::ContradictoryConstraint {
            atom_id: AtomId(0), key: AtomConstraintKey::AromaticValence,
        })
    )]
    fn test_valence_invariants_check_molecule_atom_contributions(
        #[case] input: &str,
        #[case] expected: Solution<(), ValenceMismatch>,
    ) {
        let molecule = mol_dsl!(input);
        assert_eq!(
            ValenceInvariants::check_molecule_atom(&molecule, AtomId(0)),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ground_methane(
        AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(4),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::not_aromatic(
        AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(4),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: AtomConstraintsForm::from(AtomConstraintForm::aromatic_valence(
                AromaticValenceForm::NotAromatic,
            )),
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::aromatic_zero(
        AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(4),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: AtomConstraintsForm::from(AtomConstraintForm::aromatic_valence(
                AromaticValenceForm::Aromatic(NumForm::Lit(0)),
            )),
            ..Default::default()
        },
        Solution::Determined(()),
    )]
    #[case::undetermined_charge(
        AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Undetermined,
            lone_pairs: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(4),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        },
        Solution::Underdetermined(()),
    )]
    #[case::orbital_count_mismatch(
        AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(99),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        },
        Solution::Contradictory(ValenceMismatch::OrbitalCount { atom_id: AtomId(0), orbital_count: 198, electron_count: 103 }),
    )]
    fn test_valence_invariants_check_atom(
        #[case] atom: AtomForm,
        #[case] expected: Solution<(), ValenceMismatch>,
    ) {
        assert_eq!(ValenceInvariants::check_atom(&atom), expected);
    }

    #[rstest]
    #[case::ground_methane(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(4),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
    )]
    fn test_valence_invariants_check(#[case] molecule: Molecule) {
        assert_eq!(
            ValenceInvariants::check(&molecule),
            Solution::Determined(())
        );
    }

    #[rstest]
    #[case::methane(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            isotope_mass: IsotopeMassForm::Natural,
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(4),
            lone_pairs: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: AtomConstraintsForm::from(AtomConstraintForm::Valence(NumForm::Lit(0))),
        }],
    )]
    #[case::infeasible_h(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(99),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![],
    )]
    // Every enumerated field already ground: the single combination is returned
    // verbatim (oxygen atom, triplet, two lone pairs).
    #[case::all_fields_ground(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::O),
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(2),
            unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomForm {
            element: ElementForm::Lit(Element::O),
            isotope_mass: IsotopeMassForm::Natural,
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(2),
            unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
            constraints: AtomConstraintsForm::from(AtomConstraintForm::Valence(NumForm::Lit(0))),
        }],
    )]
    // Unpaired electrons fully ground to a non-maximal but valid coupling (3 electrons as a
    // doublet); lone pairs open and fixed by conservation to 1.
    #[case::ground_unpaired_electrons(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::N),
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((3_u8, 2_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomForm {
            element: ElementForm::Lit(Element::N),
            isotope_mass: IsotopeMassForm::Natural,
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(1),
            unpaired_electrons: UnpairedElectronsForm::from((3_u8, 2_u8)),
            constraints: AtomConstraintsForm::from(AtomConstraintForm::Valence(NumForm::Lit(0))),
        }],
    )]
    // Multiplicity is ground (singlet), unpaired-electron count open: conservation fixes the
    // count to 2, and the meet keeps the pinned multiplicity (open-shell singlet).
    #[case::ground_multiplicity(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(2),
            lone_pairs: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm {
                count: NumForm::Undetermined,
                multiplicity: NumForm::Lit(1),
            },
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            isotope_mass: IsotopeMassForm::Natural,
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(2),
            lone_pairs: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((2_u8, 1_u8)),
            constraints: AtomConstraintsForm::from(AtomConstraintForm::Valence(NumForm::Lit(0))),
        }],
    )]
    // Unpaired electrons pinned to a physically impossible pair (count 2, multiplicity 2):
    // the only conservation-valid count is incompatible, so no candidate.
    #[case::inconsistent_unpaired_electrons(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(2),
            lone_pairs: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![],
    )]
    // Both components open: unpaired-electron count fixed by conservation to 2, multiplicity
    // defaulted by `into_concrete` to the maximum (triplet carbene).
    #[case::open_unpaired_electrons(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(2),
            lone_pairs: NumForm::Lit(0),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            isotope_mass: IsotopeMassForm::Natural,
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(2),
            lone_pairs: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
            constraints: AtomConstraintsForm::from(AtomConstraintForm::Valence(NumForm::Lit(0))),
        }],
    )]
    // Nonzero (given) charge: oxide anion, 7 electrons, resolves to a doublet.
    #[case::nonzero_charge(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::O),
            charge: NumForm::Lit(-1),
            implicit_hydrogens: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(3),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomForm {
            element: ElementForm::Lit(Element::O),
            isotope_mass: IsotopeMassForm::Natural,
            charge: NumForm::Lit(-1),
            implicit_hydrogens: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(3),
            unpaired_electrons: UnpairedElectronsForm::from((1_u8, 2_u8)),
            constraints: AtomConstraintsForm::from(AtomConstraintForm::Valence(NumForm::Lit(0))),
        }],
    )]
    // Specified isotope survives the meet (Natural can't, Lit(13) does) and is
    // preserved through `into_concrete`.
    #[case::specified_isotope(
        Molecule::from_entries(MoleculeEntries { atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            isotope_mass: IsotopeMassForm::Lit(13),
            charge: NumForm::Lit(0),
            lone_pairs: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            ..Default::default()
        }], bonds: vec![], ..Default::default() }),
        AtomId(0),
        vec![AtomForm {
            element: ElementForm::Lit(Element::C),
            isotope_mass: IsotopeMassForm::Lit(13),
            charge: NumForm::Lit(0),
            implicit_hydrogens: NumForm::Lit(4),
            lone_pairs: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: AtomConstraintsForm::from(AtomConstraintForm::Valence(NumForm::Lit(0))),
        }],
    )]
    fn test_valence_invariants_enumerate_atom(
        #[case] molecule: Molecule,
        #[case] atom_id: AtomId,
        #[case] expected: Vec<AtomForm>,
    ) {
        assert_eq!(
            ValenceInvariants::enumerate_atom(&molecule, atom_id),
            expected
        );
    }

    #[rstest]
    #[case::multicenter(
        r#"{:atoms ["H#c-#h0#u0#s" "Cr" "Cr"]
            :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
        vec!["H#i=#c-#h0#n0#u0#s#v0"]
    )]
    #[case::unknown_multicenter(
        r#"{:atoms ["H#c-#h0#u0#s" "Cr" "Cr"]
            :multicenter-bonds [{:atoms [0 1 2] :attrs "*"}]}"#, vec![]
    )]
    #[case::conflicting_multicenter(
        r#"{:atoms ["H#c-#h0#u0#s#m0" "Cr" "Cr"]
            :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#, vec![]
    )]
    #[case::donor(
        r#"{:atoms ["N#c0#h3#u0#s" "B"]
            :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#,
        vec!["N#i=#c0#h3#n0#u0#s#v0"]
    )]
    #[case::conflicting_donor(
        r#"{:atoms ["N#c0#h3#u0#s#d0" "B"]
            :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#, vec![]
    )]
    fn test_valence_invariants_enumerate_atom_contributions(
        #[case] input: &str,
        #[case] expected: Vec<&str>,
    ) {
        let molecule = mol_dsl!(input);
        assert_eq!(
            ValenceInvariants::enumerate_atom(&molecule, AtomId(0)),
            expected
                .iter()
                .map(|atom| atom_dsl!(atom))
                .collect::<Vec<_>>(),
        );
    }
}
