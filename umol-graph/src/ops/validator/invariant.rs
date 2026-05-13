//! Electron-invariant validator. For every atom, the orbital-occupancy count
//! and the source-electron count (`Z − q + neighbor contributions`) must be
//! equal.

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AtomAst, AtomConstraint, AtomConstraintKind, AtomId, MoleculeAst,
    MulticenterValenceAst, ValueAst,
};

use crate::ops::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct ElectronInvariantValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ElectronInvariantContradiction {
    #[error("atom {atom:?}: orbital count {orbital_count} != electron count {electron_count}")]
    AtomInvariantMismatch {
        atom: AtomId,
        orbital_count: i64,
        electron_count: i64,
    },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ElectronInvariantError {}

impl ElectronInvariantValidator {
    pub fn validate(
        &self,
        ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), ElectronInvariantContradiction>, ElectronInvariantError> {
        let ast = ast.as_ref();
        let mut any_undetermined = false;
        for view in ast.atoms().iter() {
            let atom = view.ast;
            let Some(element) = atom.element.literal() else {
                any_undetermined = true;
                continue;
            };
            let Some(charge) = atom.charge.literal() else {
                any_undetermined = true;
                continue;
            };
            let Some(implicit_h) = atom.implicit_hydrogens.literal() else {
                any_undetermined = true;
                continue;
            };
            let Some(lone_pairs) = atom.lone_pairs.literal() else {
                any_undetermined = true;
                continue;
            };
            let Some(unpaired) = atom.spin.unpaired.literal() else {
                any_undetermined = true;
                continue;
            };
            let valence_electrons = element.valence_electrons() as i64;

            let valence: i64 = match (view.valence_constraint(), view.valence()) {
                (Some(ValueAst::Lit(v)), _) if *v >= 0 => *v,
                (None | Some(ValueAst::Undetermined), Some(t)) => t as i64,
                _ => {
                    any_undetermined = true;
                    continue;
                }
            };
            let donated_pairs: i64 =
                match (view.donated_pairs_constraint(), view.donated_pairs()) {
                    (Some(ValueAst::Lit(v)), _) if *v >= 0 => *v,
                    (None | Some(ValueAst::Undetermined), Some(t)) => t as i64,
                    _ => {
                        any_undetermined = true;
                        continue;
                    }
                };
            let accepted_pairs: i64 =
                match (view.accepted_pairs_constraint(), view.accepted_pairs()) {
                    (Some(ValueAst::Lit(v)), _) if *v >= 0 => *v,
                    (None | Some(ValueAst::Undetermined), Some(t)) => t as i64,
                    _ => {
                        any_undetermined = true;
                        continue;
                    }
                };
            let aromatic_valence: i64 = match (
                view.aromatic_valence_constraint(),
                view.aromatic_valence(),
            ) {
                (Some(AromaticValenceAst::Aromatic(ValueAst::Lit(v))), _) if *v >= 0 => *v,
                (Some(AromaticValenceAst::NotAromatic), _) => 0,
                (None | Some(AromaticValenceAst::Undetermined), Some(t)) => t as i64,
                _ => {
                    any_undetermined = true;
                    continue;
                }
            };
            let multicenter_valence: i64 = match (
                view.multicenter_valence_constraint(),
                view.multicenter_valence(),
            ) {
                (Some(MulticenterValenceAst::Multicenter(ValueAst::Lit(v))), _) if *v >= 0 => *v,
                (Some(MulticenterValenceAst::NotMulticenter), _) => 0,
                (None | Some(MulticenterValenceAst::Undetermined), Some(t)) => t as i64,
                _ => {
                    any_undetermined = true;
                    continue;
                }
            };

            let aromatic_increment = if aromatic_valence == 1 { 1 } else { 0 };
            let orbital_count = unpaired
                + 2 * lone_pairs
                + 2 * donated_pairs
                + 2 * accepted_pairs
                + 2 * implicit_h
                + 2 * valence
                + aromatic_valence
                + aromatic_increment
                + multicenter_valence;
            let electron_count = valence_electrons - charge
                + implicit_h
                + valence
                + aromatic_increment
                + 2 * accepted_pairs;

            if orbital_count != electron_count {
                return Ok(Solution::Contradictory(
                    ElectronInvariantContradiction::AtomInvariantMismatch {
                        atom: view.id,
                        orbital_count,
                        electron_count,
                    },
                ));
            }
        }
        Ok(if any_undetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }

    pub fn validate_atom(
        &self,
        atom: &AtomAst,
    ) -> Result<Solution<(), ElectronInvariantContradiction>, ElectronInvariantError> {
        let Some(element) = atom.element.literal() else {
            return Ok(Solution::Underdetermined(()));
        };
        let Some(charge) = atom.charge.literal() else {
            return Ok(Solution::Underdetermined(()));
        };
        let Some(implicit_h) = atom.implicit_hydrogens.literal() else {
            return Ok(Solution::Underdetermined(()));
        };
        let Some(lone_pairs) = atom.lone_pairs.literal() else {
            return Ok(Solution::Underdetermined(()));
        };
        let Some(unpaired) = atom.spin.unpaired.literal() else {
            return Ok(Solution::Underdetermined(()));
        };
        let valence_electrons = element.valence_electrons() as i64;

        // Atom-only mode: no topology. Each valence defaults to 0 unless a
        // literal constraint pins it.
        let valence: i64 = match atom.constraints.get(AtomConstraintKind::Valence) {
            Some(AtomConstraint::Valence(ValueAst::Lit(v))) if *v >= 0 => *v,
            None | Some(AtomConstraint::Valence(ValueAst::Undetermined)) => 0,
            _ => return Ok(Solution::Underdetermined(())),
        };
        let donated_pairs: i64 = match atom.constraints.get(AtomConstraintKind::DonatedPairs) {
            Some(AtomConstraint::DonatedPairs(ValueAst::Lit(v))) if *v >= 0 => *v,
            None | Some(AtomConstraint::DonatedPairs(ValueAst::Undetermined)) => 0,
            _ => return Ok(Solution::Underdetermined(())),
        };
        let accepted_pairs: i64 = match atom.constraints.get(AtomConstraintKind::AcceptedPairs) {
            Some(AtomConstraint::AcceptedPairs(ValueAst::Lit(v))) if *v >= 0 => *v,
            None | Some(AtomConstraint::AcceptedPairs(ValueAst::Undetermined)) => 0,
            _ => return Ok(Solution::Underdetermined(())),
        };
        let aromatic_valence: i64 = match atom.constraints.get(AtomConstraintKind::AromaticValence)
        {
            Some(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(
                v,
            )))) if *v >= 0 => *v,
            Some(AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic)) => 0,
            None | Some(AtomConstraint::AromaticValence(AromaticValenceAst::Undetermined)) => 0,
            _ => return Ok(Solution::Underdetermined(())),
        };
        let multicenter_valence: i64 =
            match atom.constraints.get(AtomConstraintKind::MulticenterValence) {
                Some(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(
                    ValueAst::Lit(v),
                ))) if *v >= 0 => *v,
                Some(AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter)) => {
                    0
                }
                None
                | Some(AtomConstraint::MulticenterValence(MulticenterValenceAst::Undetermined)) => {
                    0
                }
                _ => return Ok(Solution::Underdetermined(())),
            };

        let aromatic_increment = if aromatic_valence == 1 { 1 } else { 0 };
        let orbital_count = unpaired
            + 2 * lone_pairs
            + 2 * donated_pairs
            + 2 * accepted_pairs
            + 2 * implicit_h
            + 2 * valence
            + aromatic_valence
            + aromatic_increment
            + multicenter_valence;
        let electron_count = valence_electrons - charge
            + implicit_h
            + valence
            + aromatic_increment
            + 2 * accepted_pairs;

        if orbital_count == electron_count {
            Ok(Solution::Determined(()))
        } else {
            Ok(Solution::Contradictory(
                ElectronInvariantContradiction::AtomInvariantMismatch {
                    atom: AtomId(0),
                    orbital_count,
                    electron_count,
                },
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AtomAst, AtomId, BondAst, Constraints, ImplicitHydrogensAst, MoleculeAst, SpinStateAst,
    };
    use umol_shared::element::Element;

    use super::*;

    fn ground_methane_atom() -> AtomAst {
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = ValueAst::Lit(0);
        atom.lone_pairs = ValueAst::Lit(0);
        atom.implicit_hydrogens = ImplicitHydrogensAst::Lit(4);
        atom.spin = SpinStateAst::from((0_u8, 1_u8));
        atom
    }

    fn ethane() -> MoleculeAst {
        let mut ch3_a = AtomAst::from_element(Element::C);
        ch3_a.charge = ValueAst::Lit(0);
        ch3_a.lone_pairs = ValueAst::Lit(0);
        ch3_a.implicit_hydrogens = ImplicitHydrogensAst::Lit(3);
        ch3_a.spin = SpinStateAst::from((0_u8, 1_u8));
        let ch3_b = ch3_a.clone();
        MoleculeAst::from_parts(
            vec![ch3_a, ch3_b],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_electron_invariant_validator_validate_atom_determined() {
        let v = ElectronInvariantValidator;
        let atom = ground_methane_atom();
        let result = v.validate_atom(&atom).unwrap();
        assert!(matches!(result, Solution::Determined(())));
    }

    #[rstest]
    fn test_electron_invariant_validator_validate_atom_underdetermined() {
        let v = ElectronInvariantValidator;
        let mut atom = ground_methane_atom();
        atom.charge = ValueAst::Undetermined;
        let result = v.validate_atom(&atom).unwrap();
        assert!(matches!(result, Solution::Underdetermined(())));
    }

    #[rstest]
    fn test_electron_invariant_validator_validate_atom_contradictory() {
        let v = ElectronInvariantValidator;
        let mut atom = ground_methane_atom();
        atom.implicit_hydrogens = ImplicitHydrogensAst::Lit(99);
        let result = v.validate_atom(&atom).unwrap();
        assert!(matches!(
            result,
            Solution::Contradictory(ElectronInvariantContradiction::AtomInvariantMismatch { .. })
        ));
    }

    #[rstest]
    fn test_electron_invariant_validator_validate_determined() {
        let v = ElectronInvariantValidator;
        let result = v.validate(ethane()).unwrap();
        assert!(matches!(result, Solution::Determined(())));
    }

    #[rstest]
    fn test_electron_invariant_validator_validate_contradictory() {
        let v = ElectronInvariantValidator;
        let mut ast = ethane();
        ast.atoms_mut().next().unwrap().implicit_hydrogens = ImplicitHydrogensAst::Lit(99);
        let result = v.validate(ast).unwrap();
        assert!(matches!(
            result,
            Solution::Contradictory(ElectronInvariantContradiction::AtomInvariantMismatch { .. })
        ));
    }
}
