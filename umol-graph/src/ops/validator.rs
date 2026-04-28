//! Tier-2 validators: physics invariants (electron count, spin coupling),
//! constraint cross-checks, and entity-structure shape checks.
//!
//! Each validator borrows a `MoleculeAst` (or `AtomAst`) and returns
//! `Result<Solution<(), C>, E>` per doc 92. Determined and Underdetermined
//! are both successful outcomes; only `Contradictory(C)` is a failure on the
//! `Solution` side. Setup-level failures (parameter-table gaps, etc.) live in
//! `Err(E)`; tier-2 validators have no setup so their `Error` types are
//! uninhabited.
//!
//! The composite [`Validator`] runs the four sub-validators in order and
//! lifts their per-engine `Contradiction` and `Error` types into unions via
//! `From` impls. `validate_atom` runs only those sub-validators that make
//! sense without a surrounding molecule (atom-typing registry use).

use thiserror::Error;
use umol_ast::ast::{
    AromaticSystemView, AromaticValenceAst, AtomAst, AtomConstraintKind, AtomIdx, AtomView,
    ImplicitHydrogensAst, MoleculeAst, MulticenterBondView, MulticenterValenceAst, ValueAst,
};
use umol_ast::ast::ElementAst;

use crate::ops::solution::Solution;

// region: ElectronInvariantValidator

#[derive(Clone, Copy, Debug, Default)]
pub struct ElectronInvariantValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ElectronInvariantContradiction {
    #[error("atom {atom:?}: orbital count {orbital_count} != electron count {electron_count}")]
    AtomInvariantMismatch {
        atom: AtomIdx,
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
            match check_electron_invariant_in_molecule(&view) {
                AtomInvariantCheck::Balanced => {}
                AtomInvariantCheck::Underdetermined => any_undetermined = true,
                AtomInvariantCheck::Mismatch {
                    orbital_count,
                    electron_count,
                } => {
                    return Ok(Solution::Contradictory(
                        ElectronInvariantContradiction::AtomInvariantMismatch {
                            atom: view.idx,
                            orbital_count,
                            electron_count,
                        },
                    ));
                }
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
        Ok(match check_electron_invariant_standalone(atom) {
            AtomInvariantCheck::Balanced => Solution::Determined(()),
            AtomInvariantCheck::Underdetermined => Solution::Underdetermined(()),
            AtomInvariantCheck::Mismatch {
                orbital_count,
                electron_count,
            } => Solution::Contradictory(ElectronInvariantContradiction::AtomInvariantMismatch {
                atom: AtomIdx(0),
                orbital_count,
                electron_count,
            }),
        })
    }
}

enum AtomInvariantCheck {
    Balanced,
    Underdetermined,
    Mismatch {
        orbital_count: i64,
        electron_count: i64,
    },
}

fn check_electron_invariant_in_molecule(view: &AtomView<'_>) -> AtomInvariantCheck {
    let atom = view.data;
    let Some(intrinsic) = read_atom_intrinsics(atom) else {
        return AtomInvariantCheck::Underdetermined;
    };

    let valence = match resolve_value(view.valence_constraint(), view.bond_order_sum()) {
        Some(v) => v,
        None => return AtomInvariantCheck::Underdetermined,
    };
    let donated_pairs =
        match resolve_value(view.donated_pairs_constraint(), view.donated_pairs()) {
            Some(v) => v,
            None => return AtomInvariantCheck::Underdetermined,
        };
    let accepted_pairs =
        match resolve_value(view.accepted_pairs_constraint(), view.accepted_pairs()) {
            Some(v) => v,
            None => return AtomInvariantCheck::Underdetermined,
        };
    let aromatic_valence = match resolve_aromatic_valence(
        view.aromatic_valence_constraint(),
        view.aromatic_contribution(),
    ) {
        Some(v) => v,
        None => return AtomInvariantCheck::Underdetermined,
    };
    let multicenter_valence = match resolve_multicenter_valence(
        view.multicenter_valence_constraint(),
        view.multicenter_contribution(),
    ) {
        Some(v) => v,
        None => return AtomInvariantCheck::Underdetermined,
    };

    evaluate_invariant(
        intrinsic,
        valence,
        donated_pairs,
        accepted_pairs,
        aromatic_valence,
        multicenter_valence,
    )
}

fn check_electron_invariant_standalone(atom: &AtomAst) -> AtomInvariantCheck {
    let Some(intrinsic) = read_atom_intrinsics(atom) else {
        return AtomInvariantCheck::Underdetermined;
    };

    // Atom-only mode: no topology. Valences default to 0 unless asserted
    // via constraints.
    let Some(valence) = resolve_value(constraint_value(atom, AtomConstraintKind::Valence), Some(0))
    else {
        return AtomInvariantCheck::Underdetermined;
    };
    let Some(donated_pairs) = resolve_value(
        constraint_value(atom, AtomConstraintKind::DonatedPairs),
        Some(0),
    ) else {
        return AtomInvariantCheck::Underdetermined;
    };
    let Some(accepted_pairs) = resolve_value(
        constraint_value(atom, AtomConstraintKind::AcceptedPairs),
        Some(0),
    ) else {
        return AtomInvariantCheck::Underdetermined;
    };
    let aromatic_valence = match resolve_aromatic_valence(
        atom_aromatic_valence_constraint(atom),
        Some(0),
    ) {
        Some(v) => v,
        None => return AtomInvariantCheck::Underdetermined,
    };
    let multicenter_valence = match resolve_multicenter_valence(
        atom_multicenter_valence_constraint(atom),
        Some(0),
    ) {
        Some(v) => v,
        None => return AtomInvariantCheck::Underdetermined,
    };

    evaluate_invariant(
        intrinsic,
        valence,
        donated_pairs,
        accepted_pairs,
        aromatic_valence,
        multicenter_valence,
    )
}

struct AtomIntrinsics {
    valence_electrons: i64,
    charge: i64,
    implicit_h: i64,
    lone_pairs: i64,
    unpaired: i64,
}

fn read_atom_intrinsics(atom: &AtomAst) -> Option<AtomIntrinsics> {
    let ElementAst::Lit(element) = atom.element else {
        return None;
    };
    let ValueAst::Lit(charge) = atom.charge else {
        return None;
    };
    let ImplicitHydrogensAst::Lit(implicit_h) = atom.implicit_hydrogens else {
        return None;
    };
    let ValueAst::Lit(lone_pairs) = atom.lone_pairs else {
        return None;
    };
    let ValueAst::Lit(unpaired) = atom.spin.unpaired else {
        return None;
    };
    Some(AtomIntrinsics {
        valence_electrons: element.valence_electrons() as i64,
        charge,
        implicit_h,
        lone_pairs,
        unpaired,
    })
}

fn evaluate_invariant(
    intrinsic: AtomIntrinsics,
    valence: u32,
    donated_pairs: u32,
    accepted_pairs: u32,
    aromatic_valence: u32,
    multicenter_valence: u32,
) -> AtomInvariantCheck {
    let valence = valence as i64;
    let donated_pairs = donated_pairs as i64;
    let accepted_pairs = accepted_pairs as i64;
    let aromatic_valence = aromatic_valence as i64;
    let multicenter_valence = multicenter_valence as i64;
    let aromatic_increment = if aromatic_valence == 1 { 1 } else { 0 };

    // Two independent counts of total valence electrons at this atom:
    //   `orbital_count`  — electrons by orbital occupancy (lone pairs,
    //                      bond pairs, π contributions, unpaired)
    //   `electron_count` — electrons by source (atom's own Z−q + electrons
    //                      contributed by each neighbor)
    // Equal for every chemically valid atom; mismatch is the tier-2 violation.
    let orbital_count = intrinsic.unpaired
        + 2 * intrinsic.lone_pairs
        + 2 * donated_pairs
        + 2 * accepted_pairs
        + 2 * intrinsic.implicit_h
        + 2 * valence
        + aromatic_valence
        + aromatic_increment
        + multicenter_valence;

    let electron_count = intrinsic.valence_electrons - intrinsic.charge
        + intrinsic.implicit_h
        + valence
        + aromatic_increment
        + 2 * accepted_pairs;

    if orbital_count == electron_count {
        AtomInvariantCheck::Balanced
    } else {
        AtomInvariantCheck::Mismatch {
            orbital_count,
            electron_count,
        }
    }
}

fn resolve_value(constraint: Option<&ValueAst>, topology: Option<u32>) -> Option<u32> {
    match constraint {
        Some(ValueAst::Lit(v)) if *v >= 0 => Some(*v as u32),
        Some(ValueAst::Lit(_)) => None,
        Some(ValueAst::Undetermined) | None => topology,
        Some(_) => None,
    }
}

fn resolve_aromatic_valence(
    constraint: Option<&AromaticValenceAst>,
    topology: Option<u32>,
) -> Option<u32> {
    match constraint {
        Some(AromaticValenceAst::Aromatic(ValueAst::Lit(v))) if *v >= 0 => Some(*v as u32),
        Some(AromaticValenceAst::NotAromatic) => Some(0),
        Some(AromaticValenceAst::Undetermined) | None => topology,
        Some(_) => None,
    }
}

fn resolve_multicenter_valence(
    constraint: Option<&MulticenterValenceAst>,
    topology: Option<u32>,
) -> Option<u32> {
    match constraint {
        Some(MulticenterValenceAst::Multicenter(ValueAst::Lit(v))) if *v >= 0 => Some(*v as u32),
        Some(MulticenterValenceAst::NotMulticenter) => Some(0),
        Some(MulticenterValenceAst::Undetermined) | None => topology,
        Some(_) => None,
    }
}

fn constraint_value(atom: &AtomAst, kind: AtomConstraintKind) -> Option<&ValueAst> {
    match atom.constraints.get(kind)? {
        umol_ast::ast::AtomConstraint::Valence(v)
        | umol_ast::ast::AtomConstraint::DonatedPairs(v)
        | umol_ast::ast::AtomConstraint::AcceptedPairs(v)
        | umol_ast::ast::AtomConstraint::Degree(v)
        | umol_ast::ast::AtomConstraint::Connectivity(v)
        | umol_ast::ast::AtomConstraint::RingConnectivity(v)
        | umol_ast::ast::AtomConstraint::TotalHydrogens(v)
        | umol_ast::ast::AtomConstraint::RingCount(v)
        | umol_ast::ast::AtomConstraint::RingSize(v) => Some(v),
        _ => None,
    }
}

fn atom_aromatic_valence_constraint(atom: &AtomAst) -> Option<&AromaticValenceAst> {
    match atom.constraints.get(AtomConstraintKind::AromaticValence)? {
        umol_ast::ast::AtomConstraint::AromaticValence(v) => Some(v),
        _ => None,
    }
}

fn atom_multicenter_valence_constraint(atom: &AtomAst) -> Option<&MulticenterValenceAst> {
    match atom.constraints.get(AtomConstraintKind::MulticenterValence)? {
        umol_ast::ast::AtomConstraint::MulticenterValence(v) => Some(v),
        _ => None,
    }
}

// endregion: ElectronInvariantValidator

// region: SpinCouplingValidator

/// Per-entity spin-coupling parity check: a literal `(unpaired,
/// multiplicity)` pair must satisfy `multiplicity = unpaired - 2k + 1` for
/// some `k ∈ 0..=unpaired/2`. Runs on any entity carrying a `SpinStateAst`
/// (atom, aromatic system, multicenter bond).
///
/// Stub: always returns `Determined`. Implementation pending; the parity
/// rule is in `umol_shared::spin::SpinState::are_compatible`.
#[derive(Clone, Copy, Debug, Default)]
pub struct SpinCouplingValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinCouplingContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinCouplingError {}

impl SpinCouplingValidator {
    pub fn validate(
        &self,
        _ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), SpinCouplingContradiction>, SpinCouplingError> {
        Ok(Solution::Determined(()))
    }

    pub fn validate_atom(
        &self,
        _atom: &AtomAst,
    ) -> Result<Solution<(), SpinCouplingContradiction>, SpinCouplingError> {
        Ok(Solution::Determined(()))
    }
}

// endregion: SpinCouplingValidator

// region: ConstraintValidator

/// Cross-check between local atom constraints and topology-derived values
/// across all entity types, plus molecule-scope constraint evaluation
/// (`:connected`, `:total-charge`, etc.).
///
/// Stub: always returns `Determined`. Filled in once the per-relation
/// constraint evaluators land.
#[derive(Clone, Copy, Debug, Default)]
pub struct ConstraintValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConstraintContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConstraintError {}

impl ConstraintValidator {
    pub fn validate(
        &self,
        _ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        Ok(Solution::Determined(()))
    }
}

// endregion: ConstraintValidator

// region: EntityStructureValidator

/// Structural shape checks on per-relation entities: the
/// `electrons: Vec<ValueAst>` field on each aromatic system and multicenter
/// bond must match the participants list in length.
#[derive(Clone, Copy, Debug, Default)]
pub struct EntityStructureValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum EntityStructureContradiction {
    #[error(
        "aromatic system: electrons.len() = {electrons_len} but atoms.len() = {atoms_len}"
    )]
    AromaticSystemElectronsLengthMismatch {
        electrons_len: usize,
        atoms_len: usize,
    },
    #[error(
        "multicenter bond: electrons.len() = {electrons_len} but atoms.len() = {atoms_len}"
    )]
    MulticenterElectronsLengthMismatch {
        electrons_len: usize,
        atoms_len: usize,
    },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum EntityStructureError {}

impl EntityStructureValidator {
    pub fn validate(
        &self,
        ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), EntityStructureContradiction>, EntityStructureError> {
        let ast = ast.as_ref();
        for view in ast.aromatic_systems().iter() {
            if let Some(c) = aromatic_system_length_check(&view) {
                return Ok(Solution::Contradictory(c));
            }
        }
        for view in ast.multicenter_bonds().iter() {
            if let Some(c) = multicenter_length_check(&view) {
                return Ok(Solution::Contradictory(c));
            }
        }
        Ok(Solution::Determined(()))
    }
}

fn aromatic_system_length_check(
    view: &AromaticSystemView<'_>,
) -> Option<EntityStructureContradiction> {
    let atoms_len = view.atoms().count();
    let electrons_len = view.data.electrons.len();
    if electrons_len != 0 && electrons_len != atoms_len {
        Some(
            EntityStructureContradiction::AromaticSystemElectronsLengthMismatch {
                electrons_len,
                atoms_len,
            },
        )
    } else {
        None
    }
}

fn multicenter_length_check(
    view: &MulticenterBondView<'_>,
) -> Option<EntityStructureContradiction> {
    let atoms_len = view.atoms().count();
    let electrons_len = view.data.electrons.len();
    if electrons_len != 0 && electrons_len != atoms_len {
        Some(
            EntityStructureContradiction::MulticenterElectronsLengthMismatch {
                electrons_len,
                atoms_len,
            },
        )
    } else {
        None
    }
}

// endregion: EntityStructureValidator

// region: composite Validator

#[derive(Clone, Copy, Debug, Default)]
pub struct Validator {
    pub electron_invariant: ElectronInvariantValidator,
    pub spin_coupling: SpinCouplingValidator,
    pub constraint: ConstraintValidator,
    pub entity_structure: EntityStructureValidator,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValidatorContradiction {
    #[error(transparent)]
    ElectronInvariant(#[from] ElectronInvariantContradiction),
    #[error(transparent)]
    SpinCoupling(#[from] SpinCouplingContradiction),
    #[error(transparent)]
    Constraint(#[from] ConstraintContradiction),
    #[error(transparent)]
    EntityStructure(#[from] EntityStructureContradiction),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValidatorError {
    #[error(transparent)]
    ElectronInvariant(#[from] ElectronInvariantError),
    #[error(transparent)]
    SpinCoupling(#[from] SpinCouplingError),
    #[error(transparent)]
    Constraint(#[from] ConstraintError),
    #[error(transparent)]
    EntityStructure(#[from] EntityStructureError),
}

impl Validator {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn validate(
        &self,
        ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), ValidatorContradiction>, ValidatorError> {
        let ast = ast.as_ref();
        let mut any_undetermined = false;

        // Run validators in order. First contradiction wins.
        match self.entity_structure.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.electron_invariant.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.spin_coupling.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.constraint.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
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
    ) -> Result<Solution<(), ValidatorContradiction>, ValidatorError> {
        let mut any_undetermined = false;

        match self.electron_invariant.validate_atom(atom)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.spin_coupling.validate_atom(atom)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }

        Ok(if any_undetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }
}

// endregion: composite Validator

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemAst, AtomAst, AtomConstraint, AtomIdx, BondAst, Constraints, ImplicitHydrogensAst,
        MoleculeAst, MulticenterBondAst, SpinStateAst, ValueAst,
    };
    use umol_shared::element::Element;

    use super::*;

    fn ground_methane_atom() -> AtomAst {
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = ValueAst::Lit(0);
        atom.lone_pairs = ValueAst::Lit(0);
        atom.implicit_hydrogens = ImplicitHydrogensAst::Lit(4);
        atom.spin = SpinStateAst::new(0, 1);
        atom
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

    fn ethane() -> MoleculeAst {
        let mut ch3_a = AtomAst::from_element(Element::C);
        ch3_a.charge = ValueAst::Lit(0);
        ch3_a.lone_pairs = ValueAst::Lit(0);
        ch3_a.implicit_hydrogens = ImplicitHydrogensAst::Lit(3);
        ch3_a.spin = SpinStateAst::new(0, 1);
        let ch3_b = ch3_a.clone();
        MoleculeAst::new(
            vec![ch3_a, ch3_b],
            vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
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

    #[rstest]
    fn test_entity_structure_validator_aromatic_length_mismatch() {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ];
        let aromatic = vec![(
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
            AromaticSystemAst::new(
                vec![ValueAst::Lit(1), ValueAst::Lit(1)],
                ValueAst::Lit(0),
                SpinStateAst::default(),
            ),
        )];
        let ast = MoleculeAst::new(
            atoms,
            vec![],
            vec![],
            aromatic,
            vec![],
            vec![],
            Constraints::default(),
        );
        let v = EntityStructureValidator;
        let result = v.validate(ast).unwrap();
        assert!(matches!(
            result,
            Solution::Contradictory(
                EntityStructureContradiction::AromaticSystemElectronsLengthMismatch {
                    electrons_len: 2,
                    atoms_len: 3,
                }
            )
        ));
    }

    #[rstest]
    fn test_entity_structure_validator_multicenter_length_mismatch() {
        let atoms = vec![
            AtomAst::from_element(Element::B),
            AtomAst::from_element(Element::B),
            AtomAst::from_element(Element::H),
        ];
        let multicenter = vec![(
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
            MulticenterBondAst::new(
                vec![ValueAst::Lit(1)],
                ValueAst::Lit(0),
                SpinStateAst::default(),
            ),
        )];
        let ast = MoleculeAst::new(
            atoms,
            vec![],
            vec![],
            vec![],
            multicenter,
            vec![],
            Constraints::default(),
        );
        let v = EntityStructureValidator;
        let result = v.validate(ast).unwrap();
        assert!(matches!(
            result,
            Solution::Contradictory(
                EntityStructureContradiction::MulticenterElectronsLengthMismatch {
                    electrons_len: 1,
                    atoms_len: 3,
                }
            )
        ));
    }

    #[rstest]
    fn test_entity_structure_validator_empty_electrons_passes() {
        let atoms = vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ];
        let aromatic = vec![(
            vec![AtomIdx(0), AtomIdx(1)],
            AromaticSystemAst::default(),
        )];
        let ast = MoleculeAst::new(
            atoms,
            vec![],
            vec![],
            aromatic,
            vec![],
            vec![],
            Constraints::default(),
        );
        let v = EntityStructureValidator;
        assert!(matches!(v.validate(ast).unwrap(), Solution::Determined(())));
    }

    #[rstest]
    fn test_validator_composite_validate_determined() {
        let v = Validator::new();
        let result = v.validate(ethane()).unwrap();
        assert!(matches!(result, Solution::Determined(())));
    }

    #[rstest]
    fn test_validator_composite_validate_atom_determined() {
        let v = Validator::new();
        let atom = ground_methane_atom();
        let result = v.validate_atom(&atom).unwrap();
        assert!(matches!(result, Solution::Determined(())));
    }

    #[rstest]
    fn test_validator_composite_validate_atom_with_constraint_only() {
        let v = Validator::new();
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = ValueAst::Lit(0);
        atom.lone_pairs = ValueAst::Lit(0);
        atom.implicit_hydrogens = ImplicitHydrogensAst::Lit(3);
        atom.spin = SpinStateAst::new(0, 1);
        atom.constraints
            .add(AtomConstraint::Valence(ValueAst::Lit(1)));
        let result = v.validate_atom(&atom).unwrap();
        assert!(matches!(result, Solution::Determined(())));
    }
}
