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
    AromaticSystemAst, AromaticSystemId, AromaticValenceAst, AtomConstraint, AtomConstraintKind,
    AtomId, AtomView, BondConstraint, BondId, ElementAst, ImplicitHydrogensAst, MoleculeAst,
    RingFamily, ValueAst,
};
use umol_shared::element::Element;

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
        Solution<Vec<(Vec<AtomId>, AromaticSystemAst)>, AromaticityContradiction>,
        AromaticityError,
    >
    where
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let (family, max_ring_size) = self.ring_request();
        let rings = ast.rings_with(family, max_ring_size, |_| true);

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
        systems: Vec<(Vec<AtomId>, AromaticSystemAst)>,
    ) {
        if systems.is_empty() {
            return;
        }
        let mut builder = ast.edit();
        let new_indices: Vec<AromaticSystemId> = systems
            .into_iter()
            .map(|(atoms, system_ast)| builder.add_aromatic_system(atoms, system_ast))
            .collect();
        *ast = builder.build();

        for &idx in &new_indices {
            equalize_charges(ast, idx);
        }

        let bond_ids: Vec<BondId> = new_indices
            .iter()
            .flat_map(|&idx| ast.aromatic_system(idx).bond_ids().collect::<Vec<_>>())
            .collect();
        for bond_id in bond_ids {
            let bond = ast.bond_mut(bond_id);
            bond.ast.constraints.add(BondConstraint::Aromatic);
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
    match view.ast.constraints.get(AtomConstraintKind::AromaticValence)? {
        AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(n)))
            if *n >= 0 =>
        {
            Some(*n as u8)
        }
        _ => None,
    }
}

/// Charge-delocalization equalization on a single aromatic system already
/// inserted into `ast`.
///
/// Equalization captures **delocalization**: when every atom in the system
/// is the same element, the atoms are π-MO-equivalent and any single-atom
/// localization of the system's charge is an arbitrary symmetry-breaking
/// choice. The rule rewrites each atom to its canonical neutral π state
/// `(q=0, π=K)` where `K = V(element) − σ_bonds − 2·σ_lone_pairs` (closed-
/// shell electron accounting), accumulating the per-atom charge into
/// `system.charge`.
///
/// In a heterogeneous system (the aromatic ring contains more than one
/// element), the heteroatom is the natural locus of the molecule's charge
/// — pyridinium's `+1` lives on N, boratabenzene anion's `−1` lives on B,
/// pyrylium's `+1` lives on O. Equalization would erase that chemistry.
/// The rule skips heterogeneous systems entirely.
///
/// Concrete behavior:
///
/// - Carbocyclic ions (Cp⁻, tropylium, cyclopropenium, COT²⁻, etc.):
///   monoelement, all C with K=1 → equalize.
/// - `[S₄]²⁺` (square planar, all S with K=2): monoelement → equalize.
///   The two `(+1, 1)` atoms become `(0, 2)`; system.charge = +2.
/// - Boratabenzene anion `[C₅H₅BR]⁻`: heterogeneous (B + C) → skip.
///   B keeps `(−1, 1)`; system.charge = 0.
/// - Pyridinium, pyrylium, pyrrole, furan, thiophene, borepin, borazine,
///   1,2-azaborine: heterogeneous → skip. Charges (when present) stay on
///   the heteroatom.
///
/// Spin is not modified.
fn equalize_charges(ast: &mut MoleculeAst, system_idx: AromaticSystemId) {
    let atoms: Vec<AtomId> = ast.aromatic_system(system_idx).atom_ids().collect();
    let Some(element) = monoelement(ast, &atoms) else {
        return;
    };
    let v = element.valence_electrons() as i64;

    let mut accumulated = match ast.aromatic_system(system_idx).ast.charge {
        ValueAst::Lit(c) => c,
        _ => 0,
    };
    for (i, atom_idx) in atoms.iter().copied().enumerate() {
        let view = ast.atom(atom_idx);
        let atom = view.ast;
        let lp = match atom.lone_pairs {
            ValueAst::Lit(n) => n,
            _ => continue,
        };
        // TODO: replace with `view.sigma_bond_count()` once the AtomView
        // API is fleshed out — the inline `degree + implicit_h` reach-around
        // is a stopgap.
        let implicit_h = match atom.implicit_hydrogens {
            ImplicitHydrogensAst::Lit(n) => n,
            _ => 0,
        };
        let sigma_bonds = view.neighbors().count() as i64 + implicit_h;
        let k = v - sigma_bonds - 2 * lp;
        let c = match atom.charge {
            ValueAst::Lit(c) => c,
            _ => continue,
        };
        let e = match ast.aromatic_system(system_idx).ast.electrons[i] {
            ValueAst::Lit(e) => e,
            _ => continue,
        };
        if e == k {
            continue;
        }
        ast.aromatic_system_mut(system_idx).electrons[i] = ValueAst::Lit(k);
        accumulated += c;
        let atom_mut = ast.atom_mut(atom_idx).ast;
        atom_mut.charge = ValueAst::Lit(0);
        atom_mut.constraints.add(AtomConstraint::AromaticValence(
            AromaticValenceAst::Aromatic(ValueAst::Lit(k)),
        ));
    }
    ast.aromatic_system_mut(system_idx).charge = ValueAst::Lit(accumulated);
}

/// Returns the shared element if every atom in `atoms` has a literal,
/// matching element. `None` if any atom's element is undetermined or the
/// system is heterogeneous.
fn monoelement(ast: &MoleculeAst, atoms: &[AtomId]) -> Option<Element> {
    let mut iter = atoms.iter();
    let first = match ast.atom(*iter.next()?).ast.element {
        ElementAst::Lit(el) => el,
        _ => return None,
    };
    for &idx in iter {
        match ast.atom(idx).ast.element {
            ElementAst::Lit(el) if el == first => {}
            _ => return None,
        }
    }
    Some(first)
}


#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::mol_zeroed;
    use umol_ast::ast::{
        AromaticSystemId, AromaticValenceAst, AtomAst, AtomConstraint, AtomConstraintKind,
        AtomId, BondAst, BondConstraintKind, MoleculeAst, SpinStateAst, ValueAst,
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

    fn aromatic_valence_lit(ast: &MoleculeAst, idx: AtomId) -> Option<i64> {
        match ast.atom(idx).ast.constraints.get(AtomConstraintKind::AromaticValence)? {
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
            .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondAst::from_order(1)))
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
            .map(|i| (AtomId(i), AtomId((i + 1) % 5), BondAst::from_order(1)))
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
        let system = ast.aromatic_system(AromaticSystemId(0));
        let atoms: Vec<AtomId> = system.atom_ids().collect();
        assert_eq!(atoms.len(), 6);
        let aromatic_bond_count = ast
            .bonds()
            .iter()
            .filter(|view| {
                view.ast
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

    #[rstest]
    #[case::cyclopropenium_cation(
        mol_zeroed!(r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#),
        1, vec![1, 1, 1], vec![0, 0, 0], vec![1, 1, 1],
    )]
    #[case::cot_dianion(
        mol_zeroed!(r#"{:atoms ["C #h #a" "C #c- #h #a2" "C #h #a" "C #h #a"
                                "C #h #a" "C #c- #h #a2" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"]
                               [4 5 "1"] [5 6 "1"] [6 7 "1"] [7 0 "1"]]}"#),
        -2, vec![1; 8], vec![0; 8], vec![1; 8],
    )]
    #[case::s4_dication(
        mol_zeroed!(r#"{:atoms ["S #c+ #n1 #a" "S #n1 #a2" "S #c+ #n1 #a" "S #n1 #a2"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"]]}"#),
        2, vec![2; 4], vec![0; 4], vec![2; 4],
    )]
    #[case::boratabenzene_anion(
        mol_zeroed!(r#"{:atoms ["B #c- #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        0, vec![1; 6], vec![-1, 0, 0, 0, 0, 0], vec![1; 6],
    )]
    #[case::borepin(
        mol_zeroed!(r#"{:atoms ["B #h #a0" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"] [6 0 "1"]]}"#),
        0, vec![0, 1, 1, 1, 1, 1, 1], vec![0; 7], vec![0, 1, 1, 1, 1, 1, 1],
    )]
    #[case::pyridinium(
        mol_zeroed!(r#"{:atoms ["N #c+ #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        0, vec![1; 6], vec![1, 0, 0, 0, 0, 0], vec![1; 6],
    )]
    #[case::pyrylium(
        mol_zeroed!(r#"{:atoms ["O #c+ #n1 #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        0, vec![1; 6], vec![1, 0, 0, 0, 0, 0], vec![1; 6],
    )]
    #[case::pyrrole(
        mol_zeroed!(r#"{:atoms ["N #h #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        0, vec![2, 1, 1, 1, 1], vec![0; 5], vec![2, 1, 1, 1, 1],
    )]
    #[case::furan(
        mol_zeroed!(r#"{:atoms ["O #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        0, vec![2, 1, 1, 1, 1], vec![0; 5], vec![2, 1, 1, 1, 1],
    )]
    #[case::thiophene(
        mol_zeroed!(r#"{:atoms ["S #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                       :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        0, vec![2, 1, 1, 1, 1], vec![0; 5], vec![2, 1, 1, 1, 1],
    )]
    fn test_equalize_charges(
        #[case] mut ast: MoleculeAst,
        #[case] system_charge: i64,
        #[case] electrons: Vec<i64>,
        #[case] atom_charges: Vec<i64>,
        #[case] aromatic_valences: Vec<i64>,
    ) {
        let outcome = any_hueckel()
            .find_systems(&mut ast, electrons_from_aromatic_constraint)
            .unwrap();
        let Solution::Determined(systems) = outcome else {
            panic!("expected Determined, got {outcome:?}");
        };
        any_hueckel().add_systems(&mut ast, systems);

        let system = ast.aromatic_system(AromaticSystemId(0));
        assert_eq!(system.ast.charge, ValueAst::Lit(system_charge));
        assert_eq!(
            system.ast.electrons,
            electrons.into_iter().map(ValueAst::Lit).collect::<Vec<_>>(),
        );
        for (i, (q, k)) in atom_charges.iter().zip(aromatic_valences.iter()).enumerate() {
            let idx = AtomId(i as u32);
            assert_eq!(ast.atom(idx).ast.charge, ValueAst::Lit(*q));
            assert_eq!(aromatic_valence_lit(&ast, idx), Some(*k));
        }
    }

    #[rstest]
    fn test_aromaticity_perception_hueckel_rule_no_aromatic_atom_returns_determined() {
        let perception = AromaticityPerception::new(&AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        });
        let atoms: Vec<AtomAst> = (0..6).map(|_| AtomAst::from_element(Element::C)).collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        let mut ast = MoleculeAst::from_atoms_and_bonds(
            atoms,
            bonds,
        );
        let solution = run_full(&perception, &mut ast);
        assert!(matches!(solution, Solution::Determined(())));
        assert_eq!(ast.aromatic_systems().count(), 0);
        let any_aromatic = ast.bonds().iter().any(|view| {
            view.ast
                .constraints
                .iter()
                .any(|c| c.kind() == BondConstraintKind::Aromatic)
        });
        assert!(!any_aromatic);
    }
}
