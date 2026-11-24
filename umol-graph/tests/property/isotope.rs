//! IsotopeResolver's exact inverse and idempotence laws for each isotope and wildcard policy pair.
//! Generated chains carry independently varied elements, isotopes and charges, with
//! other valence fields unresolved. The inverse property checks the projected
//! molecule as well as the roundtrip, so paired no-op implementations cannot pass.
//! Raw isotope forms additionally exercise planning agreement and rollback;
//! these include unnormalized singleton sets and unconstrained variables.

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph::ops::model::{ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::{
    IsotopePolicy, IsotopeResolver, ResolveConfig, Resolver, WildcardPolicy,
};
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, ElementForm, IsotopeMassForm, Molecule, MoleculeEntries, NumForm,
    RollbackMethod, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

proptest! {
    // Fixed-H carbon chains separate isotope completion from wildcard admission.
    #[test]
    fn test_resolver_resolve_into_isotope(
        isotopes in prop::collection::vec(prop_oneof![
            Just(IsotopeMassForm::Undetermined), Just(IsotopeMassForm::Natural),
            (12u32..15).prop_map(IsotopeMassForm::Lit),
        ], 1..16),
        symbolic in prop_oneof![Just(IsotopeMassForm::lit_set([12, 13])),
            Just(IsotopeMassForm::Var(Box::new(("mass".into(), None))))],
    ) {
        let size = isotopes.len();
        let omitted = isotopes.contains(&IsotopeMassForm::Undetermined);
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: isotopes.into_iter().enumerate().map(|(i, isotope_mass)| AtomForm {
                element: ElementForm::Lit(Element::C), isotope_mass, charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(4 - i64::from(i > 0) - i64::from(i + 1 < size)),
                ..Default::default()
            }).collect(),
            bonds: (1..size).map(|i| (AtomId((i - 1) as u32), AtomId(i as u32), BondForm {
                order: NumForm::Lit(1), charge: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::closed_shell(), ..Default::default()
            })).collect(),
            ..Default::default()
        });
        let mut symbolic_source = source.clone();
        symbolic_source.atom_mut(AtomId(0)).attributes_mut().isotope_mass = symbolic;
        let mut expected = source.clone();
        for index in 0..expected.atoms().count() {
            let mut view = expected.atom_mut(AtomId(index as u32));
            let atom = view.attributes_mut();
            if atom.isotope_mass == IsotopeMassForm::Undetermined {
                atom.isotope_mass = IsotopeMassForm::Natural;
            }
            atom.lone_pairs = NumForm::Lit(0);
            atom.unpaired_electrons = UnpairedElectronsForm::closed_shell();
        }
        for valence in [ValenceModel::smiles(), ValenceModel::default()] {
            for tie_break in [ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated] {
                let model = ChemistryModel {
                    valence: ValenceModel { tie_break, ..valence.clone() }, ..Default::default()
                };
                for isotope in [IsotopePolicy::Strict, IsotopePolicy::Natural] {
                    for wildcards in [WildcardPolicy::Strict, WildcardPolicy::Admit] {
                        let resolver = Resolver::with_config(&model, ResolveConfig { isotope, wildcards, ..Default::default() });
                        let mut molecule = source.clone();
                        let result = resolver.resolve_into(&mut molecule, RollbackMethod::Journal).unwrap();
                        if isotope == IsotopePolicy::Strict && omitted {
                            prop_assert_eq!(&result, &Solution::Underdetermined(()));
                            prop_assert_eq!(&molecule, &source);
                        } else {
                            prop_assert_eq!(&result, &Solution::Determined(()));
                            prop_assert_eq!(&molecule, &expected);
                        }
                        let once = molecule.clone();
                        prop_assert_eq!(resolver.resolve_into(&mut molecule, RollbackMethod::Journal), Ok(result));
                        prop_assert_eq!(molecule, once);
                        let mut symbolic = symbolic_source.clone();
                        prop_assert_eq!(resolver.resolve_into(&mut symbolic, RollbackMethod::Journal), Ok(Solution::Underdetermined(())));
                        prop_assert_eq!(&symbolic, &symbolic_source);
                    }
                }
            }
        }
    }

    #[test]
    fn test_isotope_resolver_project_into_roundtrip(
        fields in prop::collection::vec((
            prop_oneof![Just(ElementForm::Lit(Element::C)), Just(ElementForm::Undetermined)],
            prop_oneof![Just(IsotopeMassForm::Natural), any::<u32>().prop_map(IsotopeMassForm::Lit)],
            -3i64..4,
        ), 0..24),
    ) {
        let size = fields.len();
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: fields.into_iter().map(|(element, isotope_mass, charge)| AtomForm {
                element, isotope_mass, charge: NumForm::Lit(charge), ..Default::default()
            }).collect(),
            bonds: (1..size).map(|i| (AtomId((i - 1) as u32), AtomId(i as u32), BondForm {
                order: NumForm::Lit(1), ..Default::default()
            })).collect(),
            ..Default::default()
        });
        for policy in [IsotopePolicy::Strict, IsotopePolicy::Natural] {
            for wildcards in [WildcardPolicy::Strict, WildcardPolicy::Admit] {
                let resolver = IsotopeResolver::with_config(policy, wildcards);
                let mut projected = source.clone();
                let mut expected = source.clone();
                if policy == IsotopePolicy::Natural {
                    for index in 0..expected.atoms().count() {
                        let mut view = expected.atom_mut(AtomId(index as u32));
                        let atom = view.attributes_mut();
                        if atom.isotope_mass == IsotopeMassForm::Natural
                            && (wildcards == WildcardPolicy::Strict || atom.element != ElementForm::Undetermined)
                        {
                            atom.isotope_mass = IsotopeMassForm::Undetermined;
                        }
                    }
                }
                prop_assert_eq!(resolver.project_into(&mut projected), Ok(Solution::Determined(())));
                prop_assert_eq!(&projected, &expected);
                prop_assert_eq!(resolver.resolve_into(&mut projected), Ok(Solution::Determined(())));
                prop_assert_eq!(&projected, &source);
            }
        }
    }

    #[test]
    fn test_isotope_resolver_resolve_into_idempotence(
        fields in prop::collection::vec((
            prop_oneof![Just(ElementForm::Lit(Element::C)), Just(ElementForm::Undetermined)],
            prop_oneof![
                Just(IsotopeMassForm::Undetermined),
                Just(IsotopeMassForm::Natural),
                any::<u32>().prop_map(IsotopeMassForm::Lit),
                prop::collection::btree_set(any::<u32>(), 0..4)
                    .prop_map(|values| IsotopeMassForm::LitSet(Box::new(values))),
                Just(IsotopeMassForm::Var(Box::new(("isotope".into(), None)))),
            ],
        ), 0..24),
    ) {
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: fields.into_iter().map(|(element, isotope_mass)| AtomForm {
                element, isotope_mass, ..Default::default()
            }).collect(),
            ..Default::default()
        });
        for policy in [IsotopePolicy::Strict, IsotopePolicy::Natural] {
            for wildcards in [WildcardPolicy::Strict, WildcardPolicy::Admit] {
                let resolver = IsotopeResolver::with_config(policy, wildcards);
                let mut resolved = source.clone();
                let result = resolver.resolve_into(&mut resolved).unwrap();
                match resolver.plan_resolve(&source) {
                    Solution::Determined(edits) => {
                        let expected = source.clone().apply(edits).unwrap();
                        prop_assert_eq!(&resolved, &expected);
                        prop_assert_eq!(&result, &Solution::Determined(()));
                    }
                    Solution::Underdetermined(_) => {
                        prop_assert_eq!(&resolved, &source);
                        prop_assert_eq!(&result, &Solution::Underdetermined(()));
                    }
                    Solution::Contradictory(contradiction) => match contradiction {},
                }
                let once = resolved.clone();
                prop_assert_eq!(resolver.resolve_into(&mut resolved), Ok(result));
                prop_assert_eq!(resolved, once);
            }
        }
    }
}
