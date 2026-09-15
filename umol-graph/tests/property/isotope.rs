//! IsotopeResolver's exact inverse and idempotence laws under both policies.
//! Generated chains carry independently varied isotope and charge fields, with
//! other valence fields unresolved. The inverse property checks the projected
//! molecule as well as the roundtrip, so paired no-op implementations cannot pass.
//! Raw isotope forms additionally exercise planning agreement and nonpublication;
//! these include unnormalized singleton sets and unconstrained variables.

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph::ops::model::{ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::{IsotopePolicy, IsotopeResolver, ResolveConfig, Resolver};
use umol_graph::ops::valence::ResolveReport;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, ElementForm, IsotopeMassForm, Molecule, MoleculeEntries, NumForm,
    UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

proptest! {
    // Fixed-H carbon chains separate isotope policy from valence selection. Exact
    // publication/nonpublication and idempotence are checked for every policy/model pair.
    #[test]
    fn test_resolver_resolve_isotope(
        isotopes in prop::collection::vec(prop_oneof![
            Just(IsotopeMassForm::Undetermined), Just(IsotopeMassForm::Natural),
            (12u32..15).prop_map(IsotopeMassForm::Lit),
        ], 1..16),
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
        let mut expected = source.clone();
        expected.modify_atoms(|atom| AtomForm {
            isotope_mass: match atom.isotope_mass {
                IsotopeMassForm::Undetermined => IsotopeMassForm::Natural,
                value => value,
            },
            lone_pairs: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::closed_shell(),
            ..atom
        });
        for valence in [ValenceModel::smiles(), ValenceModel::default()] {
            for tie_break in [ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated] {
                let model = ChemistryModel {
                    valence: ValenceModel { tie_break, ..valence.clone() }, ..Default::default()
                };
                for isotope in [IsotopePolicy::Strict, IsotopePolicy::Natural] {
                    let resolver = Resolver::with_config(&model, ResolveConfig { isotope, ..Default::default() });
                    let mut molecule = source.clone();
                    let result = resolver.resolve(&mut molecule).unwrap();
                    if isotope == IsotopePolicy::Strict && omitted {
                        prop_assert_eq!(&result, &Solution::Underdetermined(ResolveReport::default()));
                        prop_assert_eq!(&molecule, &source);
                    } else {
                        prop_assert_eq!(&result, &Solution::Determined(ResolveReport::default()));
                        prop_assert_eq!(&molecule, &expected);
                    }
                    let once = molecule.clone();
                    prop_assert_eq!(resolver.resolve(&mut molecule), Ok(result));
                    prop_assert_eq!(molecule, once);
                }
            }
        }
    }

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
