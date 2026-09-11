//! Ordinary valence projection preserves retained fields and recovers the source atom state.
//! Carbon chains with closed shells or single radicals exercise both candidate sources and
//! both tie-breaks with natural composition and carbon-13. Expected states are built
//! from carbon's electron budget, independently of
//! admission. A separate counts-only pairing test exercises exact rejection and atomic failure.

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph::ops::model::{ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::valence::ValenceProjectError;
use umol_graph::ops::resolve::Resolver;
use umol_graph::ops::valence::ResolveReport;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondForm, ElementForm, IsotopeMassForm, Molecule, MoleculeEntries, NumForm,
    UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

proptest! {
    #[test]
    fn test_valence_resolver_project_roundtrip(
        radicals in prop::collection::vec(any::<bool>(), 1..24),
        isotope in prop_oneof![Just(IsotopeMassForm::Natural), Just(IsotopeMassForm::Lit(13))],
        typing in any::<bool>(),
        most_saturated in any::<bool>(),
    ) {
        let size = radicals.len();
        let atoms: Vec<_> = radicals.iter().enumerate().map(|(i, &radical)| {
            let valence = i64::from(i > 0) + i64::from(i + 1 < size);
            AtomForm {
                element: ElementForm::Lit(Element::C),
                isotope_mass: isotope.clone(),
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(4 - valence - i64::from(radical)),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm {
                    count: NumForm::Lit(i64::from(radical)),
                    multiplicity: NumForm::Lit(1 + i64::from(radical)),
                },
                constraints: Default::default(),
            }
        }).collect();
        let bonds = (1..size).map(|i| (AtomId((i-1) as u32), AtomId(i as u32), BondForm {
            order: NumForm::Lit(1), charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: Default::default(),
        })).collect();
        let original = Molecule::from_entries(MoleculeEntries {atoms, bonds, ..Default::default()});
        let mut expected = original.clone();
        expected.modify_atoms(|atom| AtomForm {
            lone_pairs: NumForm::Undetermined,
            unpaired_electrons: UnpairedElectronsForm::default(),
            ..atom
        });
        let model = ChemistryModel {valence: ValenceModel {
            tie_break: if most_saturated {ValenceTieBreak::MostSaturated} else {ValenceTieBreak::Strict},
            ..if typing {ValenceModel::default()} else {ValenceModel::smiles()}
        }, ..Default::default()};
        let resolver = Resolver::new(&model);
        let mut projected = original.clone();
        prop_assert_eq!(resolver.valence.project(&mut projected, resolver.tie_break), Ok(Solution::Determined(ResolveReport::default())));
        prop_assert_eq!(&projected, &expected);
        prop_assert_eq!(resolver.resolve(&mut projected), Ok(Solution::Determined(ResolveReport::default())));
        prop_assert_eq!(projected, original);
    }

    #[test]
    fn test_valence_resolver_project_pairing(hydrogens in 0i64..3, altered in any::<bool>()) {
        let nonbonding = 4 - hydrogens;
        let unpaired = nonbonding % 2 + if altered {2} else {0};
        let source = Molecule::from_entries(MoleculeEntries {atoms: vec![AtomForm {
            element: ElementForm::Lit(Element::C), isotope_mass: IsotopeMassForm::Natural,
            charge: NumForm::Lit(0), implicit_hydrogens: NumForm::Lit(hydrogens),
            lone_pairs: NumForm::Lit((nonbonding - unpaired) / 2),
            unpaired_electrons: UnpairedElectronsForm {count: NumForm::Lit(unpaired), multiplicity: NumForm::Lit(unpaired+1)},
            constraints: Default::default(),
        }], ..Default::default()});
        let model = ChemistryModel {valence: ValenceModel::smiles(), ..Default::default()};
        let resolver = Resolver::new(&model);
        let mut molecule = source.clone();
        let result = resolver.valence.project(&mut molecule, resolver.tie_break);
        if altered {
            prop_assert_eq!(result, Err(ValenceProjectError::AtomMismatch {atom: AtomId(0)}));
        } else {
            prop_assert_eq!(result, Ok(Solution::Determined(ResolveReport::default())));
            prop_assert_eq!(resolver.resolve(&mut molecule), Ok(Solution::Determined(ResolveReport::default())));
        }
        prop_assert_eq!(molecule, source);
    }
}
