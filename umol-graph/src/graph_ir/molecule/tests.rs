use pretty_assertions::assert_eq;
use rstest::*;
use smallvec::SmallVec;
use umol_data::Element;

use super::*;
use crate::graph_ir::atom::Atom;
use crate::graph_ir::atom_pattern::AtomPattern;
use crate::graph_ir::bond_pattern::BondPattern;
use crate::graph_ir::config::ResolveConfig;
use crate::graph_ir::molecule::Molecule;
use crate::graph_ir::molecule_builder::MoleculeBuilder;

#[fixture]
fn naphthalene_molecule() -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let atoms: Vec<AtomIndex> = (0..10)
        .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
        .collect();
    let ring1_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    for (a, b) in ring1_edges {
        builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
    }
    let ring2_edges = [(3, 6), (6, 7), (7, 8), (8, 9), (9, 4)];
    for (a, b) in ring2_edges {
        builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
    }
    let carbon: Atom = "C#v4".parse().unwrap();
    for atom in builder.atom_indices().collect::<Vec<_>>() {
        builder
            .set_atom_candidates(atom, SmallVec::from_elem(carbon, 1))
            .expect("atom should exist");
    }
    builder
        .build(&ResolveConfig::default())
        .expect("test molecule should build")
}

#[rstest]
#[case::naphthalene(naphthalene_molecule(), vec![10])]
fn test_biconnected_components(#[case] molecule: Molecule, #[case] expected_sizes: Vec<usize>) {
    let mut actual_sizes: Vec<usize> = molecule
        .biconnected_components()
        .iter()
        .map(|c| c.len())
        .collect();
    actual_sizes.sort_unstable();
    assert_eq!(actual_sizes, expected_sizes);
}

#[test]
fn test_atom_aromatic_valence_resolved_semantics() {
    let mut aromatic_builder = MoleculeBuilder::new();
    let aromatic_atom = aromatic_builder.add_atom(AtomPattern::new(Element::C));
    aromatic_builder
        .set_atom_candidates(
            aromatic_atom,
            SmallVec::from_elem("C#h#v2#a".parse::<Atom>().unwrap(), 1),
        )
        .expect("atom should exist");
    let aromatic = aromatic_builder
        .build(&ResolveConfig::default())
        .expect("aromatic molecule should build");
    assert_eq!(aromatic.atom_aromatic_valence(aromatic_atom), 1);

    let mut non_aromatic_builder = MoleculeBuilder::new();
    let non_aromatic_atom = non_aromatic_builder.add_atom(AtomPattern::new(Element::C));
    non_aromatic_builder
        .set_atom_candidates(
            non_aromatic_atom,
            SmallVec::from_elem("C#v4".parse::<Atom>().unwrap(), 1),
        )
        .expect("atom should exist");
    let non_aromatic = non_aromatic_builder
        .build(&ResolveConfig::default())
        .expect("non-aromatic molecule should build");
    assert_eq!(non_aromatic.atom_aromatic_valence(non_aromatic_atom), 0);
    assert_eq!(non_aromatic.atom_aromatic_valence(AtomIndex::new(999)), 0);
}
