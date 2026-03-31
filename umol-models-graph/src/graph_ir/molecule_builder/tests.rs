use rstest::*;
use smallvec::SmallVec;
use umol_data::Element;

use super::super::atom_pattern::AtomPattern;
use super::super::bond_pattern::BondPattern;
use super::super::config::ResolveConfig;
use super::super::molecule::Molecule;
use super::*;
use crate::graph_ir::atom::Atom;

#[fixture]
fn empty_builder() -> MoleculeBuilder {
    MoleculeBuilder::new()
}

#[fixture]
fn single_atom_builder() -> MoleculeBuilder {
    let mut builder = MoleculeBuilder::new();
    builder.add_atom(AtomPattern::new(Element::C));
    builder
}

#[fixture]
fn ring_builder(#[default(6)] n: usize) -> MoleculeBuilder {
    let mut builder = MoleculeBuilder::new();
    let atoms: Vec<AtomIndex> = (0..n)
        .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
        .collect();
    for i in 0..n {
        builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % n], BondPattern::new(1));
    }
    builder
}

#[fixture]
fn chain_builder(#[default(5)] n: usize) -> MoleculeBuilder {
    let mut builder = MoleculeBuilder::new();
    let atoms: Vec<AtomIndex> = (0..n)
        .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
        .collect();
    for i in 0..n - 1 {
        builder.add_bond_unchecked(atoms[i], atoms[i + 1], BondPattern::new(1));
    }
    builder
}

#[fixture]
fn naphthalene_builder() -> MoleculeBuilder {
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
    builder
}

#[rustfmt::skip]
    #[fixture]
    fn cubane_builder() -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..8)
            .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
            .collect();
        let edges = [
            (0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6),
            (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7),
        ];
        for (a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
        }
        builder
    }

#[fixture]
fn spiro_builder() -> MoleculeBuilder {
    let mut builder = MoleculeBuilder::new();
    let atoms: Vec<AtomIndex> = (0..5)
        .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
        .collect();
    let edges = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)];
    for (a, b) in edges {
        builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
    }
    builder
}

#[fixture]
fn bridged_builder() -> MoleculeBuilder {
    let mut builder = MoleculeBuilder::new();
    let atoms: Vec<AtomIndex> = (0..6)
        .map(|_| builder.add_atom(AtomPattern::new(Element::C)))
        .collect();
    let ring1_edges = [(0, 2), (2, 1), (1, 3), (3, 0)];
    for (a, b) in ring1_edges {
        builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
    }
    let ring2_edges = [(0, 4), (4, 1), (1, 5), (5, 0)];
    for (a, b) in ring2_edges {
        builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
    }
    builder
}

#[fixture]
fn naphthalene_molecule(mut naphthalene_builder: MoleculeBuilder) -> Molecule {
    let carbon: Atom = "C#v4".parse().unwrap();
    for atom in naphthalene_builder.atom_indices().collect::<Vec<_>>() {
        naphthalene_builder
            .set_atom_candidates(atom, SmallVec::from_elem(carbon, 1))
            .expect("atom should exist");
    }
    naphthalene_builder
        .build(&ResolveConfig::default())
        .expect("test molecule should build")
}

#[rstest]
#[case::empty(empty_builder(), vec![])]
#[case::single_atom(single_atom_builder(), vec![])]
#[case::chain(chain_builder(5), vec![])]
#[case::single_ring(ring_builder(6), vec![6])]
#[case::naphthalene(naphthalene_builder(), vec![10])]
#[case::spiro(spiro_builder(), vec![3, 3])]
#[case::cubane(cubane_builder(), vec![8])]
fn test_biconnected_components(
    #[case] builder: MoleculeBuilder,
    #[case] expected_sizes: Vec<usize>,
) {
    let mut sizes: Vec<usize> = builder
        .biconnected_components()
        .iter()
        .map(|c| c.len())
        .collect();
    sizes.sort_unstable();
    assert_eq!(sizes, expected_sizes);
}

#[test]
fn test_atom_aromatic_valence_finds_aromatic_candidate_not_just_first() {
    let mut builder = MoleculeBuilder::new();
    let atom = builder.add_atom(AtomPattern::new(Element::C));
    builder
        .set_atom_candidates(
            atom,
            SmallVec::from_vec(vec![
                "C#v4".parse::<Atom>().unwrap(),
                "C#h#v2#a".parse::<Atom>().unwrap(),
            ]),
        )
        .expect("atom should exist");
    assert_eq!(builder.atom_aromatic_valence(atom), 1);
}

#[test]
fn test_atom_aromatic_valence_zero_for_non_aromatic_or_missing() {
    let mut builder = MoleculeBuilder::new();
    let atom = builder.add_atom(AtomPattern::new(Element::C));
    builder
        .set_atom_candidates(
            atom,
            SmallVec::from_elem("C#v4".parse::<Atom>().unwrap(), 1),
        )
        .expect("atom should exist");
    assert_eq!(builder.atom_aromatic_valence(atom), 0);
    assert_eq!(builder.atom_aromatic_valence(AtomIndex::new(999)), 0);
}
