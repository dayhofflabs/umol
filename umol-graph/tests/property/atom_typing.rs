//! Registry lookup agrees with filtering the supplied rows by element and valence.
//! Repeated entries preserve admission results and order.

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph::ops::valence::{AtomTypeRegistry, AtomTypingValence};
use umol_graph_ir::ir::{AtomConstraintForm, AtomForm};
use umol_graph_ir::{atom_dsl, mol_dsl};

proptest! {
    // Arbitrary insertion order, mixed charges, gaps and duplicates use a flat-list reference.
    #[test]
    fn test_atom_type_registry_lookup(
        entries in prop::collection::vec((0usize..3, 0u8..10, -3i8..4, 0u8..5), 0..60),
    ) {
        let elements = [Element::C, Element::N, Element::O];
        let rows = entries.iter().map(|&(element, valence, charge, hydrogens)| {
            AtomForm::from_element(elements[element])
                .with_charge(i64::from(charge))
                .with_implicit_hydrogens(i64::from(hydrogens))
                .with_constraints([AtomConstraintForm::valence(i64::from(valence))])
        }).collect::<Vec<_>>();
        let bulk = AtomTypeRegistry::from_atoms(rows.clone());
        let mut incremental = AtomTypeRegistry::new();
        for row in &rows {
            incremental.add(row.clone());
        }
        prop_assert_eq!(&bulk, &incremental);
        prop_assert_eq!(bulk.content_hash(), incremental.content_hash());
        for (index, &element) in elements.iter().enumerate() {
            let mut expected_all = Vec::new();
            for valence in 0..=10 {
                let expected = entries.iter().zip(&rows)
                    .filter(|((e, v, _, _), _)| *e == index && *v == valence)
                    .map(|(_, row)| row.clone()).collect::<Vec<_>>();
                prop_assert_eq!(bulk.lookup(element, valence), expected.as_slice());
                expected_all.extend(expected);
            }
            prop_assert_eq!(bulk.lookup_by_element(element), expected_all.as_slice());
        }
    }

    // Repeating any registry row preserves the admitted set and its first-occurrence order.
    #[test]
    fn test_atom_typing_valence_admit_duplicates(copies in 1usize..12) {
        let rows = [atom_dsl!("C#i*#c0#h2#n1#u0#s#v0#a!"), atom_dsl!("C#i*#c0#h2#n0#u2#s3#v0#a!")];
        let registry = AtomTypeRegistry::from_atoms(rows.clone());
        let repeated = AtomTypeRegistry::from_atoms((0..copies).flat_map(|_| rows.clone()));
        let molecule = mol_dsl!(r#"{:atoms ["C#i=#c0#h2"]}"#);
        prop_assert_eq!(AtomTypingValence::new(&repeated).admit(&molecule), AtomTypingValence::new(&registry).admit(&molecule));
    }
}
