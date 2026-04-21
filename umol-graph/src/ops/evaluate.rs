//! Constraint evaluation against a ground `Molecule` target.
//!
//! Each variant of [`AtomConstraint`], [`BondConstraint`], and [`MoleculeConstraint`]
//! has an `evaluate` method that returns `true` iff the constraint holds on the
//! supplied target. Non-ground values (e.g. an unresolved bond order) are
//! treated as failures — the evaluator rejects when required data is absent
//! rather than returning a three-valued outcome.

use umol_graph_core::NodeId;
use umol_ast::ast::atom::{ElementAst, ImplicitHydrogensAst};
use umol_shared::element::Element;
use umol_shared::spin::{SpinMultiplicity, SpinState, MAX_UNPAIRED_ELECTRONS};
use umol_ast::ast::spin::SpinStateAst;
use umol_ast::ast::value::ValueAst;

use crate::api::molecule::Molecule;
use crate::ast::constraint::{
    AromaticValenceConstraint, AtomConstraint, BondConstraint, MoleculeConstraint,
};
use crate::ast::{AromaticSystemIdx, AtomIdx, BondIdx, MulticenterBondIdx};

impl AtomConstraint {
    pub fn evaluate(&self, target: &Molecule, atom: AtomIdx) -> bool {
        match self {
            Self::Valence(v) => match target.atom_valence(atom) {
                Some(n) => v.matches_value(n as i64),
                None => false,
            },
            Self::AromaticValence(c) => match c {
                AromaticValenceConstraint::NotAromatic => !target.atom_in_aromatic_system(atom),
                AromaticValenceConstraint::Value(v) => match target.ast().atom_aromatic_valence(atom)
                {
                    Some(n) => v.matches_value(n as i64),
                    None => false,
                },
            },
            Self::MulticenterValence(v) => match target.ast().atom_multicenter_valence(atom) {
                Some(n) => v.matches_value(n as i64),
                None => false,
            },
            Self::DonatedPairs(v) => {
                let (donated, _) = target.atom_dative_bond_order_sums(atom);
                v.matches_value(donated as i64)
            }
            Self::AcceptedPairs(v) => {
                let (_, accepted) = target.atom_dative_bond_order_sums(atom);
                v.matches_value(accepted as i64)
            }
            Self::Degree(v) => v.matches_value(target.graph().degree(NodeId(atom.0)) as i64),
            Self::Connectivity(v) => match implicit_h(target, atom) {
                Some(h) => {
                    let d = target.graph().degree(NodeId(atom.0));
                    v.matches_value((d + h as usize) as i64)
                }
                None => false,
            },
            Self::TotalHCount(v) => match total_h(target, atom) {
                Some(h) => v.matches_value(h as i64),
                None => false,
            },
            Self::InRing => target.rings().contains_atom(atom),
            Self::RingCount(v) => {
                let n = target
                    .rings()
                    .iter()
                    .filter(|r| r.atoms().contains(&atom))
                    .count();
                v.matches_value(n as i64)
            }
            Self::RingSize(v) => match target.rings().atom_smallest_ring_size(atom) {
                Some(size) => v.matches_value(size as i64),
                None => false,
            },
        }
    }
}

impl BondConstraint {
    pub fn evaluate(&self, target: &Molecule, bond: BondIdx) -> bool {
        match self {
            Self::RingBond => target.rings().contains_bond(bond),
            Self::Aromatic => {
                let view = target.bond(bond);
                let (src, tgt) = (view.src, view.tgt);
                target.aromatic_systems().iter().any(|sys| {
                    let atoms: Vec<_> = sys.atoms().collect();
                    atoms.contains(&src) && atoms.contains(&tgt)
                })
            }
        }
    }
}

impl MoleculeConstraint {
    pub fn evaluate(&self, target: &Molecule) -> bool {
        match self {
            Self::AtomPred(atom, c) => c.evaluate(target, *atom),
            Self::BondPred(bond, c) => c.evaluate(target, *bond),
            Self::TotalCharge(v) => match total_charge(target) {
                Some(q) => v.matches_value(q),
                None => false,
            },
            Self::TotalSpin(expected) => match total_spin(target) {
                Some(state) => expected.matches_state(state),
                None => false,
            },
            Self::AromaticElectronCount(idx, v) => match aromatic_electron_count(target, *idx) {
                Some(n) => v.matches_value(n as i64),
                None => false,
            },
            Self::MulticenterElectronCount(idx, v) => {
                match multicenter_electron_count(target, *idx) {
                    Some(n) => v.matches_value(n as i64),
                    None => false,
                }
            }
            Self::BondOrderSum(bonds, v) => match bond_order_sum(target, bonds) {
                Some(sum) => v.matches_value(sum),
                None => false,
            },
            Self::Connected(atoms) => all_in_same_component(target, atoms),
            Self::SubPattern {
                target_anchor,
                pattern_anchor,
                pattern,
            } => {
                if pattern.atoms().count() == 0 {
                    return false;
                }
                let query = crate::api::pattern::MoleculePattern::new((**pattern).clone());
                !crate::ops::matcher::Matcher::new()
                    .find_at(&query, target, *pattern_anchor, *target_anchor)
                    .is_empty()
            }
            Self::And(xs) => xs.iter().all(|c| c.evaluate(target)),
            Self::Or(xs) => xs.iter().any(|c| c.evaluate(target)),
            Self::Not(inner) => !inner.evaluate(target),
        }
    }
}

fn implicit_h(target: &Molecule, atom: AtomIdx) -> Option<u8> {
    match &target.ast()[atom].implicit_hydrogens {
        ImplicitHydrogensAst::Value(ValueAst::Lit(n)) => Some(*n as u8),
        _ => None,
    }
}

fn total_h(target: &Molecule, atom: AtomIdx) -> Option<u8> {
    let implicit = implicit_h(target, atom)?;
    let explicit = target
        .atom_neighbors(atom)
        .filter(|n| {
            matches!(
                target.ast()[n.atom].element,
                ElementAst::Lit(Element::H)
            )
        })
        .count() as u8;
    Some(implicit + explicit)
}

fn total_charge(target: &Molecule) -> Option<i64> {
    let mut sum = 0i64;
    for atom in target.atoms().iter() {
        match atom.data.charge {
            ValueAst::Lit(n) => sum += n,
            _ => return None,
        }
    }
    Some(sum)
}

fn total_spin(target: &Molecule) -> Option<SpinState> {
    let mut unpaired = 0u32;
    for atom in target.atoms().iter() {
        match &atom.data.spin {
            SpinStateAst::from_state(s) => unpaired += s.unpaired_electrons() as u32,
            _ => return None,
        }
    }
    if unpaired > MAX_UNPAIRED_ELECTRONS as u32 {
        return None;
    }
    let u = unpaired as u8;
    let multiplicity = SpinMultiplicity::from_multiplicity(u + 1)?;
    Some(SpinState::new(u, multiplicity))
}

fn aromatic_electron_count(target: &Molecule, idx: AromaticSystemIdx) -> Option<u8> {
    let sys = target.aromatic_systems().iter().find(|s| s.idx == idx)?;
    let mut sum = 0u32;
    for atom in sys.atoms() {
        let n = target.ast().atom_aromatic_valence(atom)?;
        sum += n as u32;
    }
    if sum > u8::MAX as u32 { None } else { Some(sum as u8) }
}

fn multicenter_electron_count(target: &Molecule, idx: MulticenterBondIdx) -> Option<u8> {
    let mc = target.multicenter_bonds().iter().find(|m| m.idx == idx)?;
    let mut sum = 0u32;
    for atom in mc.atoms() {
        let n = target.ast().atom_multicenter_valence(atom)?;
        sum += n as u32;
    }
    if sum > u8::MAX as u32 { None } else { Some(sum as u8) }
}

fn bond_order_sum(target: &Molecule, bonds: &[BondIdx]) -> Option<i64> {
    let mut sum = 0i64;
    for &b in bonds {
        match target.bond(b).data.order {
            ValueAst::Lit(n) => sum += n,
            _ => return None,
        }
    }
    Some(sum)
}

fn all_in_same_component(target: &Molecule, atoms: &[AtomIdx]) -> bool {
    if atoms.len() <= 1 {
        return true;
    }
    let components = target
        .graph()
        .connected_components(umol_graph_core::ConnectedComponentsAlgorithm::Bfs);
    let mut owner = vec![usize::MAX; target.graph().node_count()];
    for (i, comp) in components.iter().enumerate() {
        for n in comp {
            owner[n.index()] = i;
        }
    }
    let first = owner[atoms[0].index()];
    atoms.iter().all(|a| owner[a.index()] == first)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_ast::ast::atom::{ElementAst, ImplicitHydrogensAst, IsotopeAst};
    use umol_shared::element::Element;
    use umol_shared::spin::SpinState;
    use umol_ast::ast::spin::SpinStateAst;
    use umol_ast::ast::value::ValueAst;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::molecule::MoleculeAst;

    #[fixture]
    fn methane() -> Molecule {
        let atom = AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Lit(4)),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::from_state(SpinState::closed_shell()),
        };
        let ast = MoleculeAst::new(
            vec![atom],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        Molecule::new(ast).unwrap()
    }

    #[fixture]
    fn cyclopropane() -> Molecule {
        let c = AtomAst {
            element: ElementAst::Lit(Element::C),
            isotope_mass: IsotopeAst::Natural,
            charge: ValueAst::Lit(0),
            implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Lit(2)),
            lone_pairs: ValueAst::Lit(0),
            spin: SpinStateAst::from_state(SpinState::closed_shell()),
        };
        let b = BondAst {
            order: ValueAst::Lit(1),
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from_state(SpinState::closed_shell()),
        };
        let ast = MoleculeAst::new(
            vec![c.clone(), c.clone(), c.clone()],
            vec![
                (AtomIdx(0), AtomIdx(1), b.clone()),
                (AtomIdx(1), AtomIdx(2), b.clone()),
                (AtomIdx(2), AtomIdx(0), b.clone()),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        Molecule::new(ast).unwrap()
    }

    #[rstest]
    #[case::valence_zero(AtomConstraint::Valence(ValueAst::Lit(0)), true)]
    #[case::valence_mismatch(AtomConstraint::Valence(ValueAst::Lit(1)), false)]
    #[case::degree_zero(AtomConstraint::Degree(ValueAst::Lit(0)), true)]
    #[case::connectivity_four(AtomConstraint::Connectivity(ValueAst::Lit(4)), true)]
    #[case::total_h_four(AtomConstraint::TotalHCount(ValueAst::Lit(4)), true)]
    #[case::total_h_mismatch(AtomConstraint::TotalHCount(ValueAst::Lit(3)), false)]
    #[case::in_ring_false(AtomConstraint::InRing, false)]
    #[case::ring_count_zero(AtomConstraint::RingCount(ValueAst::Lit(0)), true)]
    #[case::not_aromatic(
        AtomConstraint::AromaticValence(AromaticValenceConstraint::NotAromatic),
        true,
    )]
    #[case::donated_zero(AtomConstraint::DonatedPairs(ValueAst::Lit(0)), true)]
    #[case::accepted_zero(AtomConstraint::AcceptedPairs(ValueAst::Lit(0)), true)]
    fn test_atom_constraint_evaluate_methane(
        methane: Molecule,
        #[case] c: AtomConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.evaluate(&methane, AtomIdx(0)), expected);
    }

    #[rstest]
    #[case::in_ring(AtomConstraint::InRing, true)]
    #[case::ring_count(AtomConstraint::RingCount(ValueAst::Lit(1)), true)]
    #[case::ring_size_three(AtomConstraint::RingSize(ValueAst::Lit(3)), true)]
    #[case::ring_size_four(AtomConstraint::RingSize(ValueAst::Lit(4)), false)]
    #[case::degree_two(AtomConstraint::Degree(ValueAst::Lit(2)), true)]
    #[case::connectivity_four(AtomConstraint::Connectivity(ValueAst::Lit(4)), true)]
    #[case::valence_two(AtomConstraint::Valence(ValueAst::Lit(2)), true)]
    fn test_atom_constraint_evaluate_cyclopropane(
        cyclopropane: Molecule,
        #[case] c: AtomConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.evaluate(&cyclopropane, AtomIdx(0)), expected);
    }

    #[rstest]
    #[case::ring_bond_true(BondConstraint::RingBond, true)]
    #[case::aromatic_false(BondConstraint::Aromatic, false)]
    fn test_bond_constraint_evaluate(
        cyclopropane: Molecule,
        #[case] c: BondConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.evaluate(&cyclopropane, BondIdx(0)), expected);
    }

    #[rstest]
    #[case::total_charge_zero(MoleculeConstraint::TotalCharge(ValueAst::Lit(0)), true)]
    #[case::total_charge_one(MoleculeConstraint::TotalCharge(ValueAst::Lit(1)), false)]
    #[case::total_spin_singlet(
        MoleculeConstraint::TotalSpin(SpinStateAst::from_state(SpinState::closed_shell())),
        true,
    )]
    #[case::connected_single(MoleculeConstraint::Connected(vec![AtomIdx(0)]), true)]
    fn test_molecule_constraint_evaluate_methane(
        methane: Molecule,
        #[case] c: MoleculeConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.evaluate(&methane), expected);
    }

    #[rstest]
    #[case::connected_all(
        MoleculeConstraint::Connected(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]),
        true,
    )]
    #[case::bond_order_sum_three(
        MoleculeConstraint::BondOrderSum(
            vec![BondIdx(0), BondIdx(1), BondIdx(2)],
            ValueAst::Lit(3),
        ),
        true,
    )]
    #[case::bond_order_sum_mismatch(
        MoleculeConstraint::BondOrderSum(
            vec![BondIdx(0), BondIdx(1), BondIdx(2)],
            ValueAst::Lit(6),
        ),
        false,
    )]
    fn test_molecule_constraint_evaluate_cyclopropane(
        cyclopropane: Molecule,
        #[case] c: MoleculeConstraint,
        #[case] expected: bool,
    ) {
        assert_eq!(c.evaluate(&cyclopropane), expected);
    }

    #[rstest]
    fn test_molecule_constraint_evaluate_combinators(methane: Molecule) {
        let ring = MoleculeConstraint::AtomPred(AtomIdx(0), AtomConstraint::InRing);
        let charge = MoleculeConstraint::TotalCharge(ValueAst::Lit(0));
        assert!(!MoleculeConstraint::And(vec![ring.clone(), charge.clone()]).evaluate(&methane));
        assert!(MoleculeConstraint::Or(vec![ring.clone(), charge.clone()]).evaluate(&methane));
        assert!(MoleculeConstraint::Not(Box::new(ring)).evaluate(&methane));
    }

    #[rstest]
    fn test_molecule_constraint_evaluate_sub_pattern_empty_pattern(methane: Molecule) {
        let c = MoleculeConstraint::SubPattern {
            target_anchor: AtomIdx(0),
            pattern_anchor: AtomIdx(0),
            pattern: Box::new(MoleculeAst::default()),
        };
        assert!(!c.evaluate(&methane));
    }

    #[rstest]
    fn test_molecule_constraint_evaluate_sub_pattern_single_atom_matches(methane: Molecule) {
        let pattern = MoleculeAst::new(
            vec![AtomAst {
                element: ElementAst::Lit(Element::C),
                isotope_mass: IsotopeAst::Natural,
                charge: ValueAst::Undetermined,
                implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Undetermined),
                lone_pairs: ValueAst::Undetermined,
                spin: SpinStateAst::default(),
            }],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let c = MoleculeConstraint::SubPattern {
            target_anchor: AtomIdx(0),
            pattern_anchor: AtomIdx(0),
            pattern: Box::new(pattern),
        };
        assert!(c.evaluate(&methane));
    }

    #[rstest]
    fn test_molecule_constraint_evaluate_sub_pattern_element_mismatch(methane: Molecule) {
        let pattern = MoleculeAst::new(
            vec![AtomAst {
                element: ElementAst::Lit(Element::N),
                isotope_mass: IsotopeAst::Natural,
                charge: ValueAst::Undetermined,
                implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Undetermined),
                lone_pairs: ValueAst::Undetermined,
                spin: SpinStateAst::default(),
            }],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        );
        let c = MoleculeConstraint::SubPattern {
            target_anchor: AtomIdx(0),
            pattern_anchor: AtomIdx(0),
            pattern: Box::new(pattern),
        };
        assert!(!c.evaluate(&methane));
    }
}

