//! IsotopeResolver's exact inverse and idempotence laws under both policies.
//! Generated chains carry independently varied isotope and charge fields, with
//! other valence fields unresolved. The inverse property checks the projected
//! molecule as well as the roundtrip, so paired no-op implementations cannot pass.
//! Raw isotope forms additionally exercise planning agreement and nonpublication;
//! these include unnormalized singleton sets and unconstrained variables.

use proptest::prelude::*;
use umol_graph::ops::resolve::{IsotopePolicy, IsotopeResolver};
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, IsotopeMassForm, Molecule, MoleculeEntries, NumForm,
};
use umol_utils::solution::Solution;

proptest! {
    #[test]
    fn test_isotope_resolver_project_roundtrip(
        fields in prop::collection::vec((
            prop_oneof![Just(IsotopeMassForm::Natural), any::<u32>().prop_map(IsotopeMassForm::Lit)],
            -3i64..4,
        ), 0..24),
    ) {
        let size = fields.len();
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: fields.into_iter().map(|(isotope_mass, charge)| AtomForm {
                isotope_mass, charge: NumForm::Lit(charge), ..Default::default()
            }).collect(),
            bonds: (1..size).map(|i| (AtomId((i - 1) as u32), AtomId(i as u32), BondForm {
                order: NumForm::Lit(1), ..Default::default()
            })).collect(),
            ..Default::default()
        });
        for policy in [IsotopePolicy::Strict, IsotopePolicy::Natural] {
            let resolver = IsotopeResolver::new(policy);
            let mut projected = source.clone();
            let mut expected = source.clone();
            if policy == IsotopePolicy::Natural {
                expected.modify_atoms(|atom| AtomForm {
                    isotope_mass: match atom.isotope_mass {
                        IsotopeMassForm::Natural => IsotopeMassForm::Undetermined,
                        isotope => isotope,
                    },
                    ..atom
                });
            }
            prop_assert_eq!(resolver.project(&mut projected), Ok(Solution::Determined(())));
            prop_assert_eq!(&projected, &expected);
            prop_assert_eq!(resolver.resolve(&mut projected), Ok(Solution::Determined(())));
            prop_assert_eq!(&projected, &source);
        }
    }

    #[test]
    fn test_isotope_resolver_resolve_idempotence(
        isotopes in prop::collection::vec(prop_oneof![
            Just(IsotopeMassForm::Undetermined),
            Just(IsotopeMassForm::Natural),
            any::<u32>().prop_map(IsotopeMassForm::Lit),
            prop::collection::btree_set(any::<u32>(), 0..4)
                .prop_map(|values| IsotopeMassForm::LitSet(Box::new(values))),
            Just(IsotopeMassForm::Var(Box::new(("isotope".into(), None)))),
        ], 0..24),
    ) {
        let source = Molecule::from_entries(MoleculeEntries {
            atoms: isotopes.into_iter().map(|isotope_mass| AtomForm {
                isotope_mass, ..Default::default()
            }).collect(),
            ..Default::default()
        });
        for policy in [IsotopePolicy::Strict, IsotopePolicy::Natural] {
            let resolver = IsotopeResolver::new(policy);
            let mut resolved = source.clone();
            let result = resolver.resolve(&mut resolved).unwrap();
            match resolver.plan(&source) {
                Solution::Determined(edits) => {
                    let mut editor = source.edit();
                    editor.transact(edits).unwrap();
                    prop_assert_eq!(&resolved, &editor.build());
                    prop_assert_eq!(&result, &Solution::Determined(()));
                }
                Solution::Underdetermined(_) => {
                    prop_assert_eq!(&resolved, &source);
                    prop_assert_eq!(&result, &Solution::Underdetermined(()));
                }
                Solution::Contradictory(contradiction) => match contradiction {},
            }
            let once = resolved.clone();
            prop_assert_eq!(resolver.resolve(&mut resolved), Ok(result));
            prop_assert_eq!(resolved, once);
        }
    }
}
