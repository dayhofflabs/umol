//! Repeated registry entries preserve admission results and order.

use proptest::prelude::*;
use umol_graph::ops::valence::{AtomTypeRegistry, AtomTypingValence};
use umol_graph_ir::{atom_dsl, mol_dsl};

proptest! {
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
