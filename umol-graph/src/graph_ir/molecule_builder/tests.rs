use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use rstest::*;
use smallvec::{smallvec, SmallVec};
use umol_shared::Element;

use super::super::atom_pattern::AtomPattern;
use super::super::bond_pattern::BondPattern;
use super::super::config::ResolveConfig;
use super::super::molecule::{AtomIndex, BondIndex, Molecule};
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

#[rstest]
#[case::non_aromatic(vec!["C#v4"], 0)]
#[case::aromatic_second(vec!["C#v4", "C#h#v2#a"], 1)]
#[case::missing(vec![], 0)]
fn test_molecule_builder_atom_aromatic_valence(
    #[case] candidates: Vec<&str>,
    #[case] expected: u8,
) {
    let mut builder = MoleculeBuilder::new();
    let atom = builder.add_atom(AtomPattern::new(Element::C));
    if !candidates.is_empty() {
        builder
            .set_atom_candidates(
                atom,
                candidates
                    .into_iter()
                    .map(|s| s.parse::<Atom>().unwrap())
                    .collect(),
            )
            .expect("atom should exist");
    }
    assert_eq!(builder.atom_aromatic_valence(atom), expected);
}

#[rstest]
#[case::empty(ResolutionContext::default())]
#[case::nonaromatic(ResolutionContext { atom_candidates: HashMap::from([(AtomIndex::new(0), smallvec!["C#h4".parse::<Atom>().unwrap()])]),
    atom_aromatic_hints: HashMap::from([(AtomIndex::new(0), false)]), bond_aromatic_hints: HashMap::new(), atom_normal_implicit_hydrogens: HashSet::from([AtomIndex::new(0)]) })]
fn test_resolution_context_display(#[case] ctx: ResolutionContext) {
    let roundtripped = ResolutionContext::from_str(&ctx.to_string()).unwrap();
    assert_eq!(roundtripped.atom_candidates, ctx.atom_candidates);
    assert_eq!(roundtripped.atom_aromatic_hints, ctx.atom_aromatic_hints);
    assert_eq!(roundtripped.bond_aromatic_hints, ctx.bond_aromatic_hints);
    assert_eq!(roundtripped.atom_normal_implicit_hydrogens, ctx.atom_normal_implicit_hydrogens);
}

#[rstest]
#[case::nonaromatic("{:atom-candidates {0 [\"C#h4\"]} :atom-aromatic-hints {0 false} :bond-aromatic-hints {} :atom-normal-implicit-hydrogens #{0}}",
    ResolutionContext { atom_candidates: HashMap::from([(AtomIndex::new(0), smallvec!["C#h4".parse::<Atom>().unwrap()])]), atom_aromatic_hints: HashMap::from([(AtomIndex::new(0), false)]),
    bond_aromatic_hints: HashMap::new(), atom_normal_implicit_hydrogens: HashSet::from([AtomIndex::new(0)]) })]
#[case::aromatic("{:atom-candidates {0 [\"C#h#v2#a\"] 1 [\"C#h#v2#a\"] 2 [\"C#h#v2#a\"] 3 [\"C#h#v2#a\"] 4 [\"C#h#v2#a\"] 5 [\"C#h#v2#a\"]} :atom-aromatic-hints {0 true 1 true 2 true 3 true 4 true 5 true}
                   :bond-aromatic-hints {0 true 1 true 2 true 3 true 4 true 5 true} :atom-normal-implicit-hydrogens #{0 1 2 3 4 5}}",
    ResolutionContext { atom_candidates: HashMap::from([ (AtomIndex::new(0), smallvec!["C#h#v2#a".parse::<Atom>().unwrap()]), (AtomIndex::new(1), smallvec!["C#h#v2#a".parse::<Atom>().unwrap()]),
                                                         (AtomIndex::new(2), smallvec!["C#h#v2#a".parse::<Atom>().unwrap()]), (AtomIndex::new(3), smallvec!["C#h#v2#a".parse::<Atom>().unwrap()]),
                                                         (AtomIndex::new(4), smallvec!["C#h#v2#a".parse::<Atom>().unwrap()]), (AtomIndex::new(5), smallvec!["C#h#v2#a".parse::<Atom>().unwrap()])]),
                       atom_aromatic_hints: HashMap::from([ (AtomIndex::new(0), true), (AtomIndex::new(1), true), (AtomIndex::new(2), true), (AtomIndex::new(3), true), (AtomIndex::new(4), true), (AtomIndex::new(5), true) ]),
                       bond_aromatic_hints: HashMap::from([ (BondIndex::new(0), true), (BondIndex::new(1), true), (BondIndex::new(2), true), (BondIndex::new(3), true), (BondIndex::new(4), true), (BondIndex::new(5), true) ]),
                       atom_normal_implicit_hydrogens: HashSet::from([ AtomIndex::new(0), AtomIndex::new(1), AtomIndex::new(2), AtomIndex::new(3), AtomIndex::new(4), AtomIndex::new(5) ]) })]
fn test_resolution_context_from_str(#[case] edn: &str, #[case] expected: ResolutionContext) {
    let ctx = ResolutionContext::from_str(edn).unwrap();
    assert_eq!(ctx.atom_candidates, expected.atom_candidates);
    assert_eq!(ctx.atom_aromatic_hints, expected.atom_aromatic_hints);
    assert_eq!(ctx.bond_aromatic_hints, expected.bond_aromatic_hints);
    assert_eq!(
        ctx.atom_normal_implicit_hydrogens,
        expected.atom_normal_implicit_hydrogens
    );
}
