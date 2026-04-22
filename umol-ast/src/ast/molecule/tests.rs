use pretty_assertions::assert_eq;
use rstest::*;
use umol_graph_core::{
    BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm, CycleEnumerationAlgorithm,
    MaxIndependentSetAlgorithm, MaxMatchingAlgorithm, MatchingEnumerationAlgorithm,
    ShortestCycleAlgorithm,
};
use umol_shared::element::Element;

use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::atom::{AtomAst, ElementAst};
use crate::ast::bond::BondAst;
use crate::ast::constraint::{Constraint, Constraints, MoleculeConstraint};
use crate::ast::dative::DativeBondAst;
use crate::ast::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentKind};
use crate::ast::rings::RingFamily;
use crate::ast::value::ValueAst;

use super::MoleculeAst;

fn ground_atom() -> AtomAst {
    let mut a = AtomAst::from_element(Element::C);
    a.isotope_mass = crate::ast::atom::IsotopeAst::Natural;
    a.charge = ValueAst::Lit(0);
    a.implicit_hydrogens = crate::ast::atom::ImplicitHydrogensAst::Value(ValueAst::Lit(4));
    a.lone_pairs = ValueAst::Lit(0);
    a.spin = crate::ast::spin::SpinStateAst::new(0, 1);
    a
}

#[test]
fn test_molecule_ast_is_ground_empty() {
    assert!(MoleculeAst::default().is_ground());
}

#[test]
fn test_molecule_ast_is_ground_atom() {
    let ast = MoleculeAst::new(
        vec![ground_atom()],
        vec![],
        vec![],
        vec![],
        vec![],
        vec![],
        Constraints::default(),
    );
    assert!(ast.is_ground());
}

#[test]
fn test_molecule_ast_is_ground_wildcard_element() {
    let ast = MoleculeAst::new(
        vec![AtomAst::new(ElementAst::Undetermined)],
        vec![],
        vec![],
        vec![],
        vec![],
        vec![],
        Constraints::default(),
    );
    assert!(!ast.is_ground());
}

#[test]
fn test_molecule_ast_is_ground_wildcard_bond() {
    let ast = MoleculeAst::new(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
        ],
        vec![(AtomIdx(0), AtomIdx(1), BondAst::new(ValueAst::Undetermined))],
        vec![],
        vec![],
        vec![],
        vec![],
        Constraints::default(),
    );
    assert!(!ast.is_ground());
}

#[test]
fn test_molecule_ast_is_ground_ignores_constraints() {
    let mut ast = MoleculeAst::new(
        vec![ground_atom()],
        vec![],
        vec![],
        vec![],
        vec![],
        vec![],
        Constraints::default(),
    );
    ast.constraints_mut()
        .push_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: vec![],
            sum: ValueAst::Undetermined,
        }));
    assert!(ast.is_ground());
}

#[test]
fn test_molecule_ast_neighbors() {
    let ast = MoleculeAst::new(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::N),
        ],
        vec![
            (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
            (AtomIdx(0), AtomIdx(2), BondAst::from_order(2)),
        ],
        vec![],
        vec![],
        vec![],
        vec![],
        Constraints::default(),
    );
    assert_eq!(ast.neighbors(AtomIdx(0)).count(), 2);
    assert_eq!(ast.neighbors(AtomIdx(1)).count(), 1);
    assert_eq!(ast.neighbors(AtomIdx(2)).count(), 1);
}

#[test]
fn test_molecule_ast_edit_add_aromatic_system() {
    let ast = MoleculeAst::new(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
        vec![],
        vec![],
        vec![],
        vec![],
        Constraints::default(),
    );
    let mut b = ast.edit();
    let id = b.add_aromatic_system(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default());
    let new_ast = b.build();
    assert_eq!(id, AromaticSystemIdx(0));
    assert_eq!(new_ast.aromatic_systems().count(), 1);
    assert_eq!(ast.aromatic_systems().count(), 0);
}

#[test]
fn test_molecule_ast_counts() {
    let ast = MoleculeAst::new(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
        ],
        vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(2))],
        vec![],
        vec![(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default())],
        vec![],
        vec![],
        Constraints::default(),
    );
    assert_eq!(ast.atom_count(), 2);
    assert_eq!(ast.bond_count(), 1);
    assert_eq!(ast.aromatic_system_count(), 1);
    assert_eq!(ast.dative_bond_count(), 0);
    assert_eq!(ast.multicenter_bond_count(), 0);
    assert_eq!(ast.noncovalent_bond_count(), 0);
}

fn rich_molecule() -> MoleculeAst {
    MoleculeAst::new(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::N),
            AtomAst::from_element(Element::O),
        ],
        vec![
            (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
            (AtomIdx(1), AtomIdx(2), BondAst::from_order(2)),
            (AtomIdx(2), AtomIdx(3), BondAst::from_order(1)),
        ],
        vec![(AtomIdx(2), AtomIdx(3), DativeBondAst::new())],
        vec![(
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
            AromaticSystemAst::default(),
        )],
        vec![(
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
            MulticenterBondAst::default(),
        )],
        vec![(
            AtomIdx(0),
            AtomIdx(3),
            NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond),
        )],
        Constraints::default(),
    )
}

#[test]
fn test_molecule_ast_bond_view() {
    let ast = rich_molecule();
    let bv = ast.bond(BondIdx(0));
    assert_eq!(bv.idx, BondIdx(0));
    assert_eq!(bv.src, AtomIdx(0));
    assert_eq!(bv.tgt, AtomIdx(1));
    assert_eq!(bv.data.order, ValueAst::Lit(1));

    let bv2 = ast.bond(BondIdx(2));
    assert_eq!(bv2.src, AtomIdx(2));
    assert_eq!(bv2.tgt, AtomIdx(3));
}

#[test]
fn test_molecule_ast_bond_views_iter() {
    let ast = rich_molecule();
    let views: Vec<_> = ast.bonds().iter().collect();
    assert_eq!(views.len(), 3);
    assert_eq!(views[0].src, AtomIdx(0));
    assert_eq!(views[1].src, AtomIdx(1));
    assert_eq!(views[2].src, AtomIdx(2));
}

#[test]
fn test_molecule_ast_dative_bond_view() {
    let ast = rich_molecule();
    let dv = ast.dative_bond(DativeBondIdx(0));
    assert_eq!(dv.idx, DativeBondIdx(0));
    assert_eq!(dv.donor, AtomIdx(2));
    assert_eq!(dv.acceptor, AtomIdx(3));
}

#[test]
fn test_molecule_ast_dative_bond_views_iter() {
    let ast = rich_molecule();
    let views: Vec<_> = ast.dative_bonds().iter().collect();
    assert_eq!(views.len(), 1);
    assert_eq!(views[0].donor, AtomIdx(2));
    assert_eq!(views[0].acceptor, AtomIdx(3));
}

#[test]
fn test_molecule_ast_aromatic_system_view() {
    let ast = rich_molecule();
    let av = ast.aromatic_system(AromaticSystemIdx(0));
    assert_eq!(av.idx, AromaticSystemIdx(0));
    let mut atoms: Vec<_> = av.atoms().collect();
    atoms.sort_unstable();
    assert_eq!(atoms, vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
    let mut bonds: Vec<_> = av.bonds().collect();
    bonds.sort_unstable();
    assert_eq!(bonds, vec![BondIdx(0), BondIdx(1)]);
}

#[test]
fn test_molecule_ast_aromatic_system_views_iter() {
    let ast = rich_molecule();
    let views: Vec<_> = ast.aromatic_systems().iter().collect();
    assert_eq!(views.len(), 1);
    assert_eq!(views[0].atoms().count(), 3);
    assert_eq!(views[0].bonds().count(), 2);
}

#[test]
fn test_molecule_ast_multicenter_bond_view() {
    let ast = rich_molecule();
    let mv = ast.multicenter_bond(MulticenterBondIdx(0));
    assert_eq!(mv.idx, MulticenterBondIdx(0));
    let mut atoms: Vec<_> = mv.atoms().collect();
    atoms.sort_unstable();
    assert_eq!(atoms, vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
}

#[test]
fn test_molecule_ast_multicenter_bond_views_iter() {
    let ast = rich_molecule();
    let views: Vec<_> = ast.multicenter_bonds().iter().collect();
    assert_eq!(views.len(), 1);
    assert_eq!(views[0].atoms().count(), 3);
}

#[test]
fn test_molecule_ast_noncovalent_bond_view() {
    let ast = rich_molecule();
    let nv = ast.noncovalent_bond(NoncovalentBondIdx(0));
    assert_eq!(nv.idx, NoncovalentBondIdx(0));
    let mut atoms = nv.atoms;
    atoms.sort_unstable();
    assert_eq!(atoms, [AtomIdx(0), AtomIdx(3)]);
}

#[test]
fn test_molecule_ast_noncovalent_bond_views_iter() {
    let ast = rich_molecule();
    let views: Vec<_> = ast.noncovalent_bonds().iter().collect();
    assert_eq!(views.len(), 1);
}

#[test]
fn test_molecule_ast_connecting_bond() {
    let ast = rich_molecule();
    assert_eq!(ast.connecting_bond(AtomIdx(0), AtomIdx(1)), Some(BondIdx(0)));
    assert_eq!(ast.connecting_bond(AtomIdx(1), AtomIdx(0)), Some(BondIdx(0)));
    assert_eq!(ast.connecting_bond(AtomIdx(0), AtomIdx(3)), None);
}

#[test]
fn test_molecule_ast_dative_bonds_incident() {
    let ast = rich_molecule();
    let inc: Vec<_> = ast.dative_bonds_incident(AtomIdx(2)).collect();
    assert_eq!(inc, vec![DativeBondIdx(0)]);
    let inc: Vec<_> = ast.dative_bonds_incident(AtomIdx(3)).collect();
    assert_eq!(inc, vec![DativeBondIdx(0)]);
    let inc: Vec<_> = ast.dative_bonds_incident(AtomIdx(0)).collect();
    assert!(inc.is_empty());
}

#[test]
fn test_molecule_ast_aromatic_systems_incident() {
    let ast = rich_molecule();
    let inc: Vec<_> = ast.aromatic_systems_incident(AtomIdx(1)).collect();
    assert_eq!(inc, vec![AromaticSystemIdx(0)]);
    let inc: Vec<_> = ast.aromatic_systems_incident(AtomIdx(3)).collect();
    assert!(inc.is_empty());
}

#[test]
fn test_molecule_ast_multicenter_bonds_incident() {
    let ast = rich_molecule();
    let inc: Vec<_> = ast.multicenter_bonds_incident(AtomIdx(0)).collect();
    assert_eq!(inc, vec![MulticenterBondIdx(0)]);
    let inc: Vec<_> = ast.multicenter_bonds_incident(AtomIdx(3)).collect();
    assert!(inc.is_empty());
}

#[test]
fn test_molecule_ast_noncovalent_bonds_incident() {
    let ast = rich_molecule();
    let inc: Vec<_> = ast.noncovalent_bonds_incident(AtomIdx(0)).collect();
    assert_eq!(inc, vec![NoncovalentBondIdx(0)]);
    let inc: Vec<_> = ast.noncovalent_bonds_incident(AtomIdx(3)).collect();
    assert_eq!(inc, vec![NoncovalentBondIdx(0)]);
    let inc: Vec<_> = ast.noncovalent_bonds_incident(AtomIdx(1)).collect();
    assert!(inc.is_empty());
}

#[test]
fn test_molecule_ast_induced_dative_bonds() {
    let ast = rich_molecule();
    assert_eq!(
        ast.induced_dative_bonds(&[AtomIdx(2), AtomIdx(3)]),
        vec![DativeBondIdx(0)]
    );
    assert!(ast.induced_dative_bonds(&[AtomIdx(0), AtomIdx(2)]).is_empty());
}

#[test]
fn test_molecule_ast_induced_aromatic_systems() {
    let ast = rich_molecule();
    assert_eq!(
        ast.induced_aromatic_systems(&[AtomIdx(0), AtomIdx(1), AtomIdx(2)]),
        vec![AromaticSystemIdx(0)]
    );
    assert!(ast.induced_aromatic_systems(&[AtomIdx(0), AtomIdx(1)]).is_empty());
}

#[test]
fn test_molecule_ast_induced_multicenter_bonds() {
    let ast = rich_molecule();
    assert_eq!(
        ast.induced_multicenter_bonds(&[AtomIdx(0), AtomIdx(1), AtomIdx(2)]),
        vec![MulticenterBondIdx(0)]
    );
    assert!(ast.induced_multicenter_bonds(&[AtomIdx(0), AtomIdx(1)]).is_empty());
}

#[test]
fn test_molecule_ast_induced_noncovalent_bonds() {
    let ast = rich_molecule();
    assert_eq!(
        ast.induced_noncovalent_bonds(&[AtomIdx(0), AtomIdx(3)]),
        vec![NoncovalentBondIdx(0)]
    );
    assert!(ast.induced_noncovalent_bonds(&[AtomIdx(0), AtomIdx(1)]).is_empty());
}

#[test]
fn test_molecule_ast_neighbor_view() {
    let ast = rich_molecule();
    let nbrs: Vec<_> = ast.neighbors(AtomIdx(1)).collect();
    assert_eq!(nbrs.len(), 2);
    assert!(nbrs.iter().any(|n| n.atom == AtomIdx(0) && n.bond == BondIdx(0)));
    assert!(nbrs.iter().any(|n| n.atom == AtomIdx(2) && n.bond == BondIdx(1)));
}

#[test]
fn test_molecule_ast_atom_view() {
    let ast = rich_molecule();
    let av = ast.atom(AtomIdx(2));
    assert_eq!(av.idx, AtomIdx(2));
    assert_eq!(av.data.element, ElementAst::Lit(Element::N));
}

#[test]
fn test_molecule_ast_atom_views_iter() {
    let ast = rich_molecule();
    let views: Vec<_> = ast.atoms().iter().collect();
    assert_eq!(views.len(), 4);
    assert_eq!(views[0].idx, AtomIdx(0));
    assert_eq!(views[3].idx, AtomIdx(3));
}

#[test]
fn test_molecule_ast_induced_bonds() {
    let ast = MoleculeAst::new(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        vec![
            (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
            (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
            (AtomIdx(0), AtomIdx(2), BondAst::from_order(1)),
        ],
        vec![],
        vec![],
        vec![],
        vec![],
        Constraints::default(),
    );
    let bonds = ast.induced_bonds(&[AtomIdx(0), AtomIdx(1)]);
    assert_eq!(bonds, vec![BondIdx(0)]);

    let mut all = ast.induced_bonds(&[AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
    all.sort_unstable();
    assert_eq!(all, vec![BondIdx(0), BondIdx(1), BondIdx(2)]);
}

fn chain(n: usize) -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); n];
    let bonds: Vec<_> = (0..n.saturating_sub(1))
        .map(|i| {
            (
                AtomIdx(i as u32),
                AtomIdx((i + 1) as u32),
                BondAst::from_order(1),
            )
        })
        .collect();
    MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default())
}

fn ring(n: usize) -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); n];
    let bonds: Vec<_> = (0..n)
        .map(|i| {
            (
                AtomIdx(i as u32),
                AtomIdx(((i + 1) % n) as u32),
                BondAst::from_order(1),
            )
        })
        .collect();
    MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default())
}

fn two_components() -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); 4];
    let bonds = vec![
        (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
        (AtomIdx(2), AtomIdx(3), BondAst::from_order(1)),
    ];
    MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default())
}

#[rstest]
#[case::isolated(chain(1), AtomIdx(0), 0)]
#[case::chain_end(chain(3), AtomIdx(0), 1)]
#[case::chain_mid(chain(3), AtomIdx(1), 2)]
#[case::ring_vertex(ring(6), AtomIdx(0), 2)]
fn test_molecule_ast_degree(
    #[case] ast: MoleculeAst,
    #[case] atom: AtomIdx,
    #[case] expected: usize,
) {
    assert_eq!(ast.degree(atom), expected);
}

#[rstest]
#[case::single(chain(3), 1)]
#[case::two(two_components(), 2)]
#[case::empty(MoleculeAst::default(), 0)]
fn test_molecule_ast_connected_components(#[case] ast: MoleculeAst, #[case] expected: usize) {
    let cc = ast.connected_components(ConnectedComponentsAlgorithm::Bfs);
    assert_eq!(cc.len(), expected);
}

#[rstest]
#[case::ring_6(ring(6), 1)]
#[case::chain(chain(5), 0)]
fn test_molecule_ast_biconnected_components(
    #[case] ast: MoleculeAst,
    #[case] expected: usize,
) {
    let bcc = ast.biconnected_components(BiconnectedComponentsAlgorithm::Tarjan);
    assert_eq!(bcc.len(), expected);
}

#[rstest]
#[case::ring_bond(ring(6), BondIdx(0), Some(6))]
#[case::chain_bond(chain(3), BondIdx(0), None)]
fn test_molecule_ast_shortest_cycle_through_bond(
    #[case] ast: MoleculeAst,
    #[case] bond: BondIdx,
    #[case] expected: Option<usize>,
) {
    assert_eq!(
        ast.shortest_cycle_through_bond(bond, ShortestCycleAlgorithm::Bfs),
        expected
    );
}

#[rstest]
#[case::ring_atom(ring(6), AtomIdx(0), Some(6))]
#[case::chain_atom(chain(3), AtomIdx(1), None)]
fn test_molecule_ast_shortest_cycle_through_atom(
    #[case] ast: MoleculeAst,
    #[case] atom: AtomIdx,
    #[case] expected: Option<usize>,
) {
    assert_eq!(
        ast.shortest_cycle_through_atom(atom, ShortestCycleAlgorithm::Bfs),
        expected
    );
}

#[rstest]
#[case::hexagon(ring(6), 6, 1)]
#[case::hexagon_cutoff(ring(6), 5, 0)]
#[case::chain(chain(5), 10, 0)]
#[case::empty(MoleculeAst::default(), 10, 0)]
fn test_molecule_ast_enumerate_cycles(
    #[case] ast: MoleculeAst,
    #[case] max_size: usize,
    #[case] expected: usize,
) {
    let cycles = ast.enumerate_cycles(max_size, CycleEnumerationAlgorithm::Vismara);
    assert_eq!(cycles.len(), expected);
}

#[rstest]
#[case::triangle(ring(3), 1)]
#[case::chain_3(chain(3), 2)]
fn test_molecule_ast_maximum_independent_set(
    #[case] ast: MoleculeAst,
    #[case] expected: usize,
) {
    let mis = ast.maximum_independent_set(MaxIndependentSetAlgorithm::BranchAndBound);
    assert_eq!(mis.len(), expected);
}

#[rstest]
#[case::chain_4(chain(4), 2)]
#[case::ring_6(ring(6), 3)]
#[case::single(chain(1), 0)]
fn test_molecule_ast_maximum_matching(
    #[case] ast: MoleculeAst,
    #[case] expected_size: usize,
) {
    let m = ast.maximum_matching(MaxMatchingAlgorithm::Edmonds);
    assert_eq!(m.size(), expected_size);
}

#[test]
fn test_bond_matching_mate() {
    let ast = chain(4);
    let m = ast.maximum_matching(MaxMatchingAlgorithm::Edmonds);
    assert!(m.is_matched(AtomIdx(0)));
    let mate = m.mate(AtomIdx(0));
    assert!(mate.is_some());
}

#[rstest]
#[case::ring_6(ring(6), 2)]
fn test_molecule_ast_enumerate_perfect_matchings(
    #[case] ast: MoleculeAst,
    #[case] expected: usize,
) {
    let ms = ast.enumerate_perfect_matchings(MatchingEnumerationAlgorithm::BranchAndBound);
    assert_eq!(ms.len(), expected);
    for m in &ms {
        assert!(m.is_perfect(ast.atom_count()));
    }
}

#[rstest]
#[case::ring_6(ring(6), 1)]
#[case::chain_3(chain(3), 2)]
fn test_molecule_ast_automorphisms(
    #[case] ast: MoleculeAst,
    #[case] expected_orbits: usize,
) {
    let auto = ast.automorphisms(|_| 0u8, umol_graph_core::AutomorphismAlgorithm::Nauty);
    assert_eq!(auto.num_orbits(), expected_orbits);
    assert_eq!(auto.atom_count(), ast.atom_count());
}

#[test]
fn test_atom_automorphism_same_orbit() {
    let ast = ring(6);
    let auto = ast.automorphisms(|_| 0u8, umol_graph_core::AutomorphismAlgorithm::Nauty);
    assert!(auto.same_orbit(AtomIdx(0), AtomIdx(3)));
}

#[test]
fn test_molecule_ast_subgraph_isomorphisms() {
    let target = ring(6);
    let query = chain(2);
    let matches = target.subgraph_isomorphisms(
        &query,
        &mut |_, _| true,
        &mut |_, _| true,
        umol_graph_core::SubgraphIsomorphismAlgorithm::Vf2,
    );
    assert_eq!(matches.len(), 12);
}

#[test]
fn test_molecule_ast_subgraph_isomorphisms_at() {
    let target = ring(6);
    let query = chain(2);
    let matches = target.subgraph_isomorphisms_at(
        &query,
        (AtomIdx(0), AtomIdx(0)),
        &mut |_, _| true,
        &mut |_, _| true,
        umol_graph_core::SubgraphIsomorphismAlgorithm::Vf2,
    );
    assert_eq!(matches.len(), 2);
}

#[test]
fn test_molecule_ast_induced_subgraph() {
    let ast = rich_molecule();
    let sub = ast.induced_subgraph(&[AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
    assert_eq!(sub.ast.atom_count(), 3);
    assert_eq!(sub.ast.bond_count(), 2);
    assert_eq!(sub.atom_map, vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
    assert_eq!(sub.bond_map, vec![BondIdx(0), BondIdx(1)]);
    assert_eq!(sub.ast.aromatic_system_count(), 1);
    assert_eq!(sub.aromatic_system_map, vec![AromaticSystemIdx(0)]);
    assert_eq!(sub.ast.multicenter_bond_count(), 1);
    assert_eq!(sub.multicenter_bond_map, vec![MulticenterBondIdx(0)]);
    assert_eq!(sub.ast.dative_bond_count(), 0);
    assert!(sub.dative_bond_map.is_empty());
    assert_eq!(sub.ast.noncovalent_bond_count(), 0);
    assert!(sub.noncovalent_bond_map.is_empty());
}

#[test]
fn test_molecule_ast_induced_subgraph_preserves_dative() {
    let ast = rich_molecule();
    let sub = ast.induced_subgraph(&[AtomIdx(2), AtomIdx(3)]);
    assert_eq!(sub.ast.atom_count(), 2);
    assert_eq!(sub.ast.dative_bond_count(), 1);
    assert_eq!(sub.dative_bond_map, vec![DativeBondIdx(0)]);
}

#[test]
fn test_builder_remove_aromatic_systems() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.remove_aromatic_systems(&[AromaticSystemIdx(0)]);
    let result = b.build();
    assert_eq!(result.aromatic_system_count(), 0);
    assert_eq!(result.atom_count(), 4);
    assert_eq!(result.bond_count(), 3);
}

#[test]
fn test_builder_remove_dative_bonds() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.remove_dative_bonds(&[DativeBondIdx(0)]);
    let result = b.build();
    assert_eq!(result.dative_bond_count(), 0);
    assert_eq!(result.atom_count(), 4);
}

#[test]
fn test_builder_remove_multicenter_bonds() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.remove_multicenter_bonds(&[MulticenterBondIdx(0)]);
    let result = b.build();
    assert_eq!(result.multicenter_bond_count(), 0);
}

#[test]
fn test_builder_remove_noncovalent_bonds() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.remove_noncovalent_bonds(&[NoncovalentBondIdx(0)]);
    let result = b.build();
    assert_eq!(result.noncovalent_bond_count(), 0);
}

#[test]
fn test_builder_atom_mut() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.atom_mut(AtomIdx(0)).element = ElementAst::Lit(Element::N);
    let result = b.build();
    assert_eq!(result[AtomIdx(0)].element, ElementAst::Lit(Element::N));
    assert_eq!(ast[AtomIdx(0)].element, ElementAst::Lit(Element::C));
}

#[test]
fn test_builder_bond_mut() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.bond_mut(BondIdx(0)).order = ValueAst::Lit(3);
    let result = b.build();
    assert_eq!(result[BondIdx(0)].order, ValueAst::Lit(3));
    assert_eq!(ast[BondIdx(0)].order, ValueAst::Lit(1));
}

#[test]
fn test_builder_constraints_mut() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.constraints_mut()
        .push_atom(AtomIdx(0), crate::ast::constraint::AtomConstraint::Degree(ValueAst::Lit(2)));
    let result = b.build();
    assert_eq!(result.constraints().atom(AtomIdx(0)).len(), 1);
    assert!(ast.constraints().atom(AtomIdx(0)).is_empty());
}

#[rstest]
#[case::hexagon(ring(6), 6, 1)]
#[case::hexagon_cutoff(ring(6), 5, 0)]
#[case::chain(chain(5), 10, 0)]
#[case::empty(MoleculeAst::default(), 10, 0)]
fn test_molecule_ast_rings(
    #[case] ast: MoleculeAst,
    #[case] max_ring_size: usize,
    #[case] expected: usize,
) {
    let rs = ast.rings(RingFamily::Simple, max_ring_size, |_| true);
    assert_eq!(rs.count(), expected);
}

#[test]
fn test_molecule_ast_rings_atom_filter() {
    let ast = ring(6);
    let rs = ast.rings(RingFamily::Simple, 10, |a| a.0 < 3);
    assert_eq!(rs.count(), 0);
}

#[test]
fn test_molecule_ast_rings_induced() {
    let atoms = vec![AtomAst::from_element(Element::C); 4];
    let bonds = vec![
        (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
        (AtomIdx(0), AtomIdx(2), BondAst::from_order(1)),
        (AtomIdx(0), AtomIdx(3), BondAst::from_order(1)),
        (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
        (AtomIdx(1), AtomIdx(3), BondAst::from_order(1)),
        (AtomIdx(2), AtomIdx(3), BondAst::from_order(1)),
    ];
    let ast = MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default());
    let simple = ast.rings(RingFamily::Simple, 4, |_| true);
    let induced = ast.rings(RingFamily::Induced, 4, |_| true);
    assert_eq!(simple.count(), 4);
    assert_eq!(induced.count(), 4);
}

#[test]
fn test_molecule_ast_rings_induced_naphthalene() {
    let atoms = vec![AtomAst::from_element(Element::C); 10];
    #[rustfmt::skip]
    let bonds = vec![
        (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
        (AtomIdx(1), AtomIdx(2), BondAst::from_order(1)),
        (AtomIdx(2), AtomIdx(3), BondAst::from_order(1)),
        (AtomIdx(3), AtomIdx(4), BondAst::from_order(1)),
        (AtomIdx(4), AtomIdx(5), BondAst::from_order(1)),
        (AtomIdx(5), AtomIdx(0), BondAst::from_order(1)),
        (AtomIdx(3), AtomIdx(6), BondAst::from_order(1)),
        (AtomIdx(6), AtomIdx(7), BondAst::from_order(1)),
        (AtomIdx(7), AtomIdx(8), BondAst::from_order(1)),
        (AtomIdx(8), AtomIdx(9), BondAst::from_order(1)),
        (AtomIdx(9), AtomIdx(4), BondAst::from_order(1)),
    ];
    let ast = MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], Constraints::default());
    let simple = ast.rings(RingFamily::Simple, 10, |_| true);
    assert_eq!(simple.count(), 2);
    let induced = ast.rings(RingFamily::Induced, 10, |_| true);
    assert_eq!(induced.count(), 2);
}

#[test]
fn test_rings_membership() {
    let ast = ring(6);
    let rs = ast.rings(RingFamily::Simple, 6, |_| true);
    assert!(rs.contains_atom(AtomIdx(0)));
    assert!(rs.contains_bond(BondIdx(0)));
    assert_eq!(rs.atom_smallest_ring_size(AtomIdx(0)), Some(6));
}

#[test]
fn test_dpo_add_then_remove() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    let new_a = b.add_atom(AtomAst::from_element(Element::Br));
    b.add_bond(AtomIdx(0), new_a, BondAst::from_order(1));
    b.remove_aromatic_systems(&[AromaticSystemIdx(0)]);
    let _remap = b.remove(&[AtomIdx(3)], &[BondIdx(2)]);
    let result = b.build();
    assert_eq!(result.atom_count(), 4);
    assert_eq!(result.bond_count(), 3);
    assert_eq!(result.aromatic_system_count(), 0);
    assert_eq!(result.dative_bond_count(), 0);
    assert_eq!(result.noncovalent_bond_count(), 0);
}
