//! Aromaticity perception primitive.
//!
//! [`AromaticityPerception`] dispatches to one of three algorithms (Hückel
//! rule, HMO, Clar) selected by [`AromaticityModel`] and runs perception
//! against an AST. It is the shared core used by three top-level entities:
//! the resolver (validates `#a` hints filled in by atom-typing), the
//! aromatizer (discovers aromatic systems from a Kekulé bond-order layout),
//! and the validator (verifies pre-existing aromatic systems against the
//! model). Each entity composes [`AromaticityPerception::find_systems`] with
//! its own per-atom electron-counting closure.
//!
//! Mutation — pattern-match charge equalization, system insertion, bond
//! marking — is exposed via [`AromaticityPerception::add_systems`].

pub mod clar;
pub mod hmo;
pub mod hueckel_rule;

pub use clar::{ClarAromaticity, ClarError};
pub use hmo::{HmoAromaticity, HmoError, HmoOutput};
pub use hueckel_rule::HueckelRuleAromaticity;
use thiserror::Error;
use umol_ast::ast::{
    AromaticSystemAst, AromaticSystemIdx, AromaticValenceAst, AtomConstraint, AtomConstraintKind,
    AtomIdx, AtomView, BondConstraint, BondIdx, MoleculeAst, RingFamily, ValueAst,
};

use crate::ops::config::AromaticityModel;
use crate::ops::solution::Solution;

/// Chemistry-level rejection: the algorithm decided the input doesn't satisfy
/// the model. Carried inside `Solution::Contradictory`.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityContradiction {
    #[error("hmo: invalid input: {0}")]
    HmoInvalidInput(String),
    #[error("clar: non-benzenoid input: {0}")]
    ClarNonBenzenoid(String),
}

/// Setup-level failure: parameter table or configuration gap. Returned in
/// `Err`, never inside `Solution`.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    #[error("hmo: missing parameters: {0}")]
    HmoMissingParameters(String),
}

#[derive(Clone, Debug)]
pub enum AromaticityPerception {
    HueckelRule(HueckelRuleAromaticity),
    Hmo(HmoAromaticity),
    Clar(ClarAromaticity),
}

impl AromaticityPerception {
    pub fn new(model: &AromaticityModel) -> Self {
        match model {
            AromaticityModel::HueckelRule { scope, ring_limits } => Self::HueckelRule(
                HueckelRuleAromaticity::new(scope.clone(), ring_limits.clone()),
            ),
            AromaticityModel::Hmo {
                scope,
                stabilization_threshold,
            } => Self::Hmo(HmoAromaticity::new(scope.clone(), *stabilization_threshold)),
            AromaticityModel::Clar { .. } => Self::Clar(ClarAromaticity),
        }
    }

    /// Find candidate aromatic systems via the configured algorithm. Takes
    /// `&mut MoleculeAst` only so that the AST's ring cache can populate; no
    /// chemistry mutation happens here. The closure `electrons_at` returns
    /// each atom's π contribution if the atom is aromatic-eligible, else
    /// `None`. Each caller passes a different closure: see
    /// [`electrons_from_aromatic_constraint`] for the resolver / validator
    /// case (reads `#a` from the atom constraints); the aromatizer derives π
    /// from bond orders.
    pub fn find_systems<F>(
        &self,
        ast: &mut MoleculeAst,
        electrons_at: F,
    ) -> Result<
        Solution<Vec<(Vec<AtomIdx>, AromaticSystemAst)>, AromaticityContradiction>,
        AromaticityError,
    >
    where
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let (family, max_ring_size) = self.ring_request();
        let rings = ast.rings(family, max_ring_size).clone();

        let systems = match self {
            Self::HueckelRule(m) => m.find_from_rings(ast, &rings, &electrons_at),
            Self::Hmo(m) => match m.find_from_rings(ast, &rings, &electrons_at) {
                Ok(systems) => systems,
                Err(HmoError::MissingParameters(s)) => {
                    return Err(AromaticityError::HmoMissingParameters(s));
                }
                Err(HmoError::InvalidInput(s)) => {
                    return Ok(Solution::Contradictory(
                        AromaticityContradiction::HmoInvalidInput(s),
                    ));
                }
                Err(HmoError::UndeterminedAtom(_)) => {
                    return Ok(Solution::Underdetermined(Vec::new()));
                }
            },
            Self::Clar(m) => match m.find_from_rings(ast, &rings, &electrons_at) {
                Ok(systems) => systems,
                Err(ClarError::NonBenzenoid(s)) => {
                    return Ok(Solution::Contradictory(
                        AromaticityContradiction::ClarNonBenzenoid(s),
                    ));
                }
            },
        };

        let mut sorted = systems;
        sorted.sort_by(|a, b| a.0.first().cmp(&b.0.first()));
        Ok(Solution::Determined(sorted))
    }

    /// Add perceived systems to the AST. Insert system entries, apply
    /// pattern-match charge equalization, mark induced bonds aromatic.
    /// Mutates `ast`.
    pub fn add_systems(
        &self,
        ast: &mut MoleculeAst,
        systems: Vec<(Vec<AtomIdx>, AromaticSystemAst)>,
    ) {
        if systems.is_empty() {
            return;
        }
        let mut builder = ast.edit();
        let new_indices: Vec<AromaticSystemIdx> = systems
            .into_iter()
            .map(|(atoms, system_ast)| builder.add_aromatic_system(atoms, system_ast))
            .collect();
        *ast = builder.build();

        for &idx in &new_indices {
            equalize_charges(ast, idx);
        }

        let bond_ids: Vec<BondIdx> = new_indices
            .iter()
            .flat_map(|&idx| ast.aromatic_system(idx).bonds().collect::<Vec<_>>())
            .collect();
        for bond_id in bond_ids {
            let bond = ast.bond_mut(bond_id);
            bond.data.constraints.add(BondConstraint::Aromatic);
        }
    }

    fn ring_request(&self) -> (RingFamily, usize) {
        match self {
            Self::HueckelRule(m) => (RingFamily::Simple, m.ring_limits.max_ring_size),
            Self::Hmo(_) => (RingFamily::Simple, 22),
            Self::Clar(_) => (RingFamily::Simple, 6),
        }
    }
}

/// Per-atom electron-counting closure for the resolver and validator: reads
/// the `AromaticValence::Aromatic(Lit(n))` constraint. Returns `None` if
/// the constraint is missing or non-numeric.
pub(crate) fn electrons_from_aromatic_constraint(view: &AtomView<'_>) -> Option<u8> {
    match view.data.constraints.get(AtomConstraintKind::AromaticValence)? {
        AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(n)))
            if *n >= 0 =>
        {
            Some(*n as u8)
        }
        _ => None,
    }
}

/// Pattern-match charge equalization on a single aromatic system already
/// inserted into `ast`. Examines `(atom.charge, system.electrons[i])` for
/// each member atom; the two patterns
///
/// - `(+1, 0)` (e.g., tropylium C⁺)
/// - `(-1, 2)` (e.g., Cp⁻ C⁻)
///
/// rewrite to `(0, 1)` on the atom side (`charge`, `AromaticValence`
/// constraint) and on the system side (`electrons[i]`, `system.charge`).
/// Other atoms — including pyridinium-style `(+1, 1)` and any non-`Lit`
/// values — are left untouched. Spin is not modified.
fn equalize_charges(ast: &mut MoleculeAst, system_idx: AromaticSystemIdx) {
    let atoms: Vec<AtomIdx> = ast.aromatic_system(system_idx).atoms().collect();
    let mut accumulated = match ast.aromatic_system(system_idx).data.charge {
        ValueAst::Lit(c) => c,
        _ => 0,
    };
    for (i, atom_idx) in atoms.iter().copied().enumerate() {
        let c = match ast.atom(atom_idx).data.charge {
            ValueAst::Lit(c) => c,
            _ => continue,
        };
        let e = match ast.aromatic_system(system_idx).data.electrons[i] {
            ValueAst::Lit(e) => e,
            _ => continue,
        };
        if !matches!((c, e), (1, 0) | (-1, 2)) {
            continue;
        }
        ast.aromatic_system_mut(system_idx).electrons[i] = ValueAst::Lit(1);
        accumulated += c;
        let atom_mut = ast.atom_mut(atom_idx).data;
        atom_mut.charge = ValueAst::Lit(0);
        atom_mut.constraints.add(AtomConstraint::AromaticValence(
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
        ));
    }
    ast.aromatic_system_mut(system_idx).charge = ValueAst::Lit(accumulated);
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::mol;
    use umol_ast::ast::{
        AromaticSystemIdx, AromaticValenceAst, AtomAst, AtomConstraint, AtomConstraintKind,
        AtomIdx, BondAst, BondConstraintKind, MoleculeAst, SpinStateAst, ValueAst,
    };
    use umol_shared::element::Element;

    use super::*;
    use crate::ops::config::{ElementScope, RingLimits};

    fn any_hueckel() -> AromaticityPerception {
        AromaticityPerception::new(&AromaticityModel::HueckelRule {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        })
    }

    fn aromatic_valence_lit(ast: &MoleculeAst, idx: AtomIdx) -> Option<i64> {
        match ast.atom(idx).data.constraints.get(AtomConstraintKind::AromaticValence)? {
            AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(n))) => {
                Some(*n)
            }
            _ => None,
        }
    }

    fn aromatic(element: Element, pi: i64) -> AtomAst {
        let mut atom = AtomAst::from_element(element);
        atom.charge = ValueAst::Lit(0);
        atom.spin = SpinStateAst::closed_shell();
        atom.constraints.add(AtomConstraint::AromaticValence(
            AromaticValenceAst::Aromatic(ValueAst::Lit(pi)),
        ));
        atom
    }

    fn benzene() -> MoleculeAst {
        let atoms: Vec<AtomAst> = (0..6).map(|_| aromatic(Element::C, 1)).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomIdx(i), AtomIdx((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        MoleculeAst::from_atoms_and_bonds(
            atoms,
            bonds,
        )
    }

    fn pyrrole() -> MoleculeAst {
        let atoms = vec![
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ];
        let bonds: Vec<_> = (0..5)
            .map(|i| (AtomIdx(i), AtomIdx((i + 1) % 5), BondAst::from_order(1)))
            .collect();
        MoleculeAst::from_atoms_and_bonds(
            atoms,
            bonds,
        )
    }

    fn run_full(
        perception: &AromaticityPerception,
        ast: &mut MoleculeAst,
    ) -> Solution<(), AromaticityContradiction> {
        let outcome = perception
            .find_systems(ast, electrons_from_aromatic_constraint)
            .unwrap();
        match outcome {
            Solution::Determined(systems) => {
                perception.add_systems(ast, systems);
                Solution::Determined(())
            }
            Solution::Underdetermined(_) => Solution::Underdetermined(()),
            Solution::Contradictory(c) => Solution::Contradictory(c),
        }
    }

    #[rstest]
    fn test_aromaticity_perception_hueckel_rule_benzene_writes_system() {
        let perception = AromaticityPerception::new(&AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        });
        let mut ast = benzene();
        let solution = run_full(&perception, &mut ast);
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(ast.aromatic_systems().count(), 1);
        let system = ast.aromatic_system(AromaticSystemIdx(0));
        let atoms: Vec<AtomIdx> = system.atoms().collect();
        assert_eq!(atoms.len(), 6);
        let aromatic_bond_count = ast
            .bonds()
            .iter()
            .filter(|view| {
                view.data
                    .constraints
                    .iter()
                    .any(|c| c.kind() == BondConstraintKind::Aromatic)
            })
            .count();
        assert_eq!(aromatic_bond_count, 6);
    }

    #[rstest]
    fn test_aromaticity_perception_clar_rejects_heterocycle() {
        let perception = AromaticityPerception::new(&AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        });
        let mut ast = pyrrole();
        let solution = run_full(&perception, &mut ast);
        assert!(matches!(
            solution,
            Solution::Contradictory(AromaticityContradiction::ClarNonBenzenoid(_))
        ));
    }

    // region: equalize_charges

    #[rstest]
    fn test_equalize_charges_cyclopropenium_cation() {
        // [C₃H₃]⁺: one atom at (+1, 0), two at (0, 1). After equalization,
        // all three atoms are (0, 1) and the system carries charge +1.
        let mut ast = mol!(r#"{
            :atoms ["C #c0 #a" "C #c0 #a" "C #c+ #a0"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]
        }"#);
        let outcome = any_hueckel()
            .find_systems(&mut ast, electrons_from_aromatic_constraint)
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut ast, systems);

        let system = ast.aromatic_system(AromaticSystemIdx(0));
        assert_eq!(system.data.charge, ValueAst::Lit(1));
        assert_eq!(system.data.electrons, vec![ValueAst::Lit(1); 3]);
        for i in 0..3 {
            assert_eq!(ast.atom(AtomIdx(i)).data.charge, ValueAst::Lit(0));
            assert_eq!(aromatic_valence_lit(&ast, AtomIdx(i)), Some(1));
        }
    }

    #[rstest]
    fn test_equalize_charges_cot_dianion_accumulates_to_minus_two() {
        // [C₈H₈]²⁻: six atoms at (0, 1), two at (-1, 2). Equalization fires
        // on each (-1, 2) atom; the per-atom -1 charges accumulate into
        // system.charge = -2. No "n=2" branch is required.
        let mut ast = mol!(r#"{
            :atoms ["C #c0 #a" "C #c- #a2" "C #c0 #a" "C #c0 #a"
                    "C #c0 #a" "C #c- #a2" "C #c0 #a" "C #c0 #a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"]
                    [4 5 "1"] [5 6 "1"] [6 7 "1"] [7 0 "1"]]
        }"#);
        let outcome = any_hueckel()
            .find_systems(&mut ast, electrons_from_aromatic_constraint)
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut ast, systems);

        let system = ast.aromatic_system(AromaticSystemIdx(0));
        assert_eq!(system.data.charge, ValueAst::Lit(-2));
        assert_eq!(system.data.electrons, vec![ValueAst::Lit(1); 8]);
        for i in 0..8 {
            assert_eq!(ast.atom(AtomIdx(i)).data.charge, ValueAst::Lit(0));
            assert_eq!(aromatic_valence_lit(&ast, AtomIdx(i)), Some(1));
        }
    }

    #[rstest]
    fn test_equalize_charges_pyridinium_leaves_n_charged() {
        // Pyridinium [C₅H₅NH]⁺: N at (+1, 1) — σ-side charge, π already
        // canonical. Rule skips (q+e=2 ≠ 1). N keeps its +1 charge.
        let mut ast = mol!(r#"{
            :atoms ["N #c+ #a" "C #c0 #a" "C #c0 #a" "C #c0 #a" "C #c0 #a" "C #c0 #a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#);
        let outcome = any_hueckel()
            .find_systems(&mut ast, electrons_from_aromatic_constraint)
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut ast, systems);

        let system = ast.aromatic_system(AromaticSystemIdx(0));
        assert_eq!(system.data.charge, ValueAst::Lit(0));
        assert_eq!(ast.atom(AtomIdx(0)).data.charge, ValueAst::Lit(1));
        assert_eq!(aromatic_valence_lit(&ast, AtomIdx(0)), Some(1));
        for i in 1..6 {
            assert_eq!(ast.atom(AtomIdx(i)).data.charge, ValueAst::Lit(0));
        }
    }

    #[rstest]
    fn test_equalize_charges_pyrylium_leaves_o_charged() {
        // Pyrylium [C₅H₅O]⁺: O at (+1, 1). Same shape as pyridinium.
        let mut ast = mol!(r#"{
            :atoms ["O #c+ #a" "C #c0 #a" "C #c0 #a" "C #c0 #a" "C #c0 #a" "C #c0 #a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#);
        let outcome = any_hueckel()
            .find_systems(&mut ast, electrons_from_aromatic_constraint)
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut ast, systems);

        let system = ast.aromatic_system(AromaticSystemIdx(0));
        assert_eq!(system.data.charge, ValueAst::Lit(0));
        assert_eq!(ast.atom(AtomIdx(0)).data.charge, ValueAst::Lit(1));
        assert_eq!(aromatic_valence_lit(&ast, AtomIdx(0)), Some(1));
    }

    #[rstest]
    fn test_equalize_charges_pyrrole_leaves_n_at_pi_two() {
        // Pyrrole [C₄H₄NH]: neutral N at (0, 2) — k=2 atom in canonical
        // form, no offset. Rule skips. N keeps its π=2 contribution.
        let mut ast = mol!(r#"{
            :atoms ["N #c0 #a2" "C #c0 #a" "C #c0 #a" "C #c0 #a" "C #c0 #a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        }"#);
        let outcome = any_hueckel()
            .find_systems(&mut ast, electrons_from_aromatic_constraint)
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut ast, systems);

        let system = ast.aromatic_system(AromaticSystemIdx(0));
        assert_eq!(system.data.charge, ValueAst::Lit(0));
        assert_eq!(system.data.electrons[0], ValueAst::Lit(2));
        assert_eq!(ast.atom(AtomIdx(0)).data.charge, ValueAst::Lit(0));
        assert_eq!(aromatic_valence_lit(&ast, AtomIdx(0)), Some(2));
    }

    #[rstest]
    #[case::furan(Element::O)]
    #[case::thiophene(Element::S)]
    fn test_equalize_charges_furan_thiophene_leave_heteroatom_alone(
        #[case] heteroatom: Element,
    ) {
        // Furan / thiophene: heteroatom at (0, 2). Same shape as pyrrole.
        let symbol = heteroatom.symbol();
        let source = format!(
            r#"{{
                :atoms ["{symbol} #c0 #a2" "C #c0 #a" "C #c0 #a" "C #c0 #a" "C #c0 #a"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            }}"#
        );
        let mut ast: MoleculeAst = source.parse().unwrap();
        let outcome = any_hueckel()
            .find_systems(&mut ast, electrons_from_aromatic_constraint)
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut ast, systems);

        let system = ast.aromatic_system(AromaticSystemIdx(0));
        assert_eq!(system.data.charge, ValueAst::Lit(0));
        assert_eq!(system.data.electrons[0], ValueAst::Lit(2));
        assert_eq!(ast.atom(AtomIdx(0)).data.charge, ValueAst::Lit(0));
        assert_eq!(aromatic_valence_lit(&ast, AtomIdx(0)), Some(2));
    }

    #[rstest]
    fn test_equalize_charges_s4_dication_known_gap() {
        // S₄²⁺ (square-planar): two S at (+1, 1), two S at (0, 2). Total π
        // = 6 (4n+2, n=1), Hueckel-aromatic. Per-atom q+e = 2 for every
        // atom; the k=1 rule never fires. The molecule passes through with
        // its mixed (+1,1)/(0,2) distribution intact and system.charge = 0.
        // Reaching the (0,2)-uniform equalized form would require k=2
        // equalization, which is not implemented.
        let mut ast = mol!(r#"{
            :atoms ["S #c+ #a" "S #c0 #a2" "S #c+ #a" "S #c0 #a2"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]]
        }"#);
        let outcome = any_hueckel()
            .find_systems(&mut ast, electrons_from_aromatic_constraint)
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut ast, systems);

        let system = ast.aromatic_system(AromaticSystemIdx(0));
        assert_eq!(system.data.charge, ValueAst::Lit(0));
        assert_eq!(
            system.data.electrons,
            vec![
                ValueAst::Lit(1),
                ValueAst::Lit(2),
                ValueAst::Lit(1),
                ValueAst::Lit(2),
            ],
        );
        assert_eq!(ast.atom(AtomIdx(0)).data.charge, ValueAst::Lit(1));
        assert_eq!(ast.atom(AtomIdx(1)).data.charge, ValueAst::Lit(0));
        assert_eq!(ast.atom(AtomIdx(2)).data.charge, ValueAst::Lit(1));
        assert_eq!(ast.atom(AtomIdx(3)).data.charge, ValueAst::Lit(0));
    }

    // endregion: equalize_charges

    #[rstest]
    fn test_aromaticity_perception_hueckel_rule_no_aromatic_atom_returns_determined() {
        let perception = AromaticityPerception::new(&AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        });
        let atoms: Vec<AtomAst> = (0..6).map(|_| AtomAst::from_element(Element::C)).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomIdx(i), AtomIdx((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        let mut ast = MoleculeAst::from_atoms_and_bonds(
            atoms,
            bonds,
        );
        let solution = run_full(&perception, &mut ast);
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(ast.aromatic_systems().count(), 0);
        let any_aromatic = ast.bonds().iter().any(|view| {
            view.data
                .constraints
                .iter()
                .any(|c| c.kind() == BondConstraintKind::Aromatic)
        });
        assert!(!any_aromatic);
    }
}
