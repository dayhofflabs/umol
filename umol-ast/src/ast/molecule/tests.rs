use std::collections::HashSet;

use pretty_assertions::assert_eq;
use rstest::*;
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, MatchingEnumerationAlgorithm, MaxIndependentSetAlgorithm,
    MaxMatchingAlgorithm, NodeId, ShortestCycleAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_shared::element::Element;

use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
use crate::ast::bond::BondAst;
use crate::ast::constraint::{
    AtomConstraint, Constraint, Constraints, DativeBondConstraint, MoleculeConstraint,
};
use crate::ast::dative::{DativeBondAst, DativeDirection};
use crate::ast::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentKind, NoncovalentKindAst};
use crate::ast::rings::RingFamily;
use crate::ast::spin::SpinStateAst;
use crate::ast::value::ValueAst;

use super::MoleculeAst;

fn ground_atom() -> AtomAst {
    let mut a = AtomAst::from_element(Element::C);
    a.isotope_mass = IsotopeAst::Natural;
    a.charge = ValueAst::Lit(0);
    a.implicit_hydrogens = ImplicitHydrogensAst::Lit(4);
    a.lone_pairs = ValueAst::Lit(0);
    a.spin = SpinStateAst::new(0, 1);
    a
}

fn constraints_with_molecule(c: Constraint) -> Constraints {
    let mut out = Constraints::new();
    out.push_molecule(c);
    out
}

#[rstest]
#[case::empty(MoleculeAst::default(), true)]
#[case::ground_atom(
    MoleculeAst::new(
        vec![ground_atom()],
        vec![], vec![], vec![], vec![], vec![],
        Constraints::default(),
    ),
    true,
)]
#[case::wildcard_element(
    MoleculeAst::new(
        vec![AtomAst::new(ElementAst::Undetermined)],
        vec![], vec![], vec![], vec![], vec![],
        Constraints::default(),
    ),
    false,
)]
#[case::wildcard_bond(
    MoleculeAst::new(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
        ],
        vec![(AtomIdx(0), AtomIdx(1), BondAst::new(ValueAst::Undetermined))],
        vec![], vec![], vec![], vec![],
        Constraints::default(),
    ),
    false,
)]
#[case::ground_atom_with_undetermined_constraint(
    MoleculeAst::new(
        vec![ground_atom()],
        vec![], vec![], vec![], vec![], vec![],
        constraints_with_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: vec![],
            sum: ValueAst::Undetermined,
        })),
    ),
    true,
)]
fn test_molecule_ast_is_ground(#[case] ast: MoleculeAst, #[case] expected: bool) {
    assert_eq!(ast.is_ground(), expected);
}

#[rstest]
#[case::hub(AtomIdx(0), vec![(AtomIdx(1), BondIdx(0)), (AtomIdx(2), BondIdx(1))])]
#[case::leaf_o(AtomIdx(1), vec![(AtomIdx(0), BondIdx(0))])]
#[case::leaf_n(AtomIdx(2), vec![(AtomIdx(0), BondIdx(1))])]
fn test_molecule_ast_neighbors(
    #[case] atom: AtomIdx,
    #[case] expected: Vec<(AtomIdx, BondIdx)>,
) {
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
    let nbrs: Vec<(AtomIdx, BondIdx)> =
        ast.neighbors(atom).map(|n| (n.atom, n.bond)).collect();
    assert_eq!(nbrs, expected);
}

#[rstest]
fn test_molecule_builder_add_aromatic_system() {
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
    let new_atoms: Vec<AtomIdx> = new_ast
        .aromatic_system(AromaticSystemIdx(0))
        .atoms()
        .collect();
    assert_eq!(new_atoms, vec![AtomIdx(0), AtomIdx(1)]);
    assert_eq!(
        new_ast
            .aromatic_systems()
            .ids()
            .collect::<Vec<_>>(),
        vec![AromaticSystemIdx(0)]
    );
    assert_eq!(
        ast.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemIdx>::new()
    );
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

#[rstest]
#[case::c_c(BondIdx(0), AtomIdx(0), AtomIdx(1), ValueAst::Lit(1))]
#[case::c_n(BondIdx(1), AtomIdx(1), AtomIdx(2), ValueAst::Lit(2))]
#[case::n_o(BondIdx(2), AtomIdx(2), AtomIdx(3), ValueAst::Lit(1))]
fn test_molecule_ast_bond(
    #[case] idx: BondIdx,
    #[case] src: AtomIdx,
    #[case] tgt: AtomIdx,
    #[case] order: ValueAst,
) {
    let ast = rich_molecule();
    let bv = ast.bond(idx);
    assert_eq!(bv.idx, idx);
    assert_eq!(bv.src, src);
    assert_eq!(bv.tgt, tgt);
    assert_eq!(bv.data.order, order);
}

#[rstest]
fn test_molecule_ast_bonds() {
    let ast = rich_molecule();
    let projected: Vec<(BondIdx, AtomIdx, AtomIdx, ValueAst)> = ast
        .bonds()
        .iter()
        .map(|v| (v.idx, v.src, v.tgt, v.data.order.clone()))
        .collect();
    assert_eq!(
        projected,
        vec![
            (BondIdx(0), AtomIdx(0), AtomIdx(1), ValueAst::Lit(1)),
            (BondIdx(1), AtomIdx(1), AtomIdx(2), ValueAst::Lit(2)),
            (BondIdx(2), AtomIdx(2), AtomIdx(3), ValueAst::Lit(1)),
        ]
    );
}

#[rstest]
fn test_molecule_ast_dative_bond() {
    let ast = rich_molecule();
    let dv = ast.dative_bond(DativeBondIdx(0));
    assert_eq!(dv.idx, DativeBondIdx(0));
    assert_eq!(dv.donor, AtomIdx(2));
    assert_eq!(dv.acceptor, AtomIdx(3));
    assert_eq!(dv.data.direction, DativeDirection::Forward);
}

#[rstest]
fn test_molecule_ast_dative_bonds() {
    let ast = rich_molecule();
    let projected: Vec<(DativeBondIdx, AtomIdx, AtomIdx, DativeDirection)> = ast
        .dative_bonds()
        .iter()
        .map(|v| (v.idx, v.donor, v.acceptor, v.data.direction))
        .collect();
    assert_eq!(
        projected,
        vec![(
            DativeBondIdx(0),
            AtomIdx(2),
            AtomIdx(3),
            DativeDirection::Forward,
        )]
    );
}

#[rstest]
fn test_molecule_ast_aromatic_system() {
    let ast = rich_molecule();
    let av = ast.aromatic_system(AromaticSystemIdx(0));
    assert_eq!(av.idx, AromaticSystemIdx(0));
    assert_eq!(
        av.atoms().collect::<Vec<_>>(),
        vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]
    );
    assert_eq!(
        av.bonds().collect::<Vec<_>>(),
        vec![BondIdx(0), BondIdx(1)]
    );
}

#[rstest]
fn test_molecule_ast_aromatic_systems() {
    let ast = rich_molecule();
    let projected: Vec<(AromaticSystemIdx, Vec<AtomIdx>, Vec<BondIdx>)> = ast
        .aromatic_systems()
        .iter()
        .map(|v| (v.idx, v.atoms().collect(), v.bonds().collect()))
        .collect();
    assert_eq!(
        projected,
        vec![(
            AromaticSystemIdx(0),
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
            vec![BondIdx(0), BondIdx(1)],
        )]
    );
}

#[rstest]
fn test_molecule_ast_multicenter_bond() {
    let ast = rich_molecule();
    let mv = ast.multicenter_bond(MulticenterBondIdx(0));
    assert_eq!(mv.idx, MulticenterBondIdx(0));
    assert_eq!(
        mv.atoms().collect::<Vec<_>>(),
        vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]
    );
}

#[rstest]
fn test_molecule_ast_multicenter_bonds() {
    let ast = rich_molecule();
    let projected: Vec<(MulticenterBondIdx, Vec<AtomIdx>)> = ast
        .multicenter_bonds()
        .iter()
        .map(|v| (v.idx, v.atoms().collect()))
        .collect();
    assert_eq!(
        projected,
        vec![(
            MulticenterBondIdx(0),
            vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)],
        )]
    );
}

#[rstest]
fn test_molecule_ast_noncovalent_bond() {
    let ast = rich_molecule();
    let nv = ast.noncovalent_bond(NoncovalentBondIdx(0));
    assert_eq!(nv.idx, NoncovalentBondIdx(0));
    assert_eq!(nv.atoms, [AtomIdx(0), AtomIdx(3)]);
}

#[rstest]
fn test_molecule_ast_noncovalent_bonds() {
    let ast = rich_molecule();
    let projected: Vec<(NoncovalentBondIdx, [AtomIdx; 2])> = ast
        .noncovalent_bonds()
        .iter()
        .map(|v| (v.idx, v.atoms))
        .collect();
    assert_eq!(
        projected,
        vec![(NoncovalentBondIdx(0), [AtomIdx(0), AtomIdx(3)])]
    );
}

#[rstest]
#[case::forward(AtomIdx(0), AtomIdx(1), Some(BondIdx(0)))]
#[case::reverse(AtomIdx(1), AtomIdx(0), Some(BondIdx(0)))]
#[case::non_adjacent(AtomIdx(0), AtomIdx(3), None)]
fn test_molecule_ast_connecting_bond(
    #[case] a: AtomIdx,
    #[case] b: AtomIdx,
    #[case] expected: Option<BondIdx>,
) {
    let ast = rich_molecule();
    assert_eq!(ast.connecting_bond(a, b), expected);
}

#[rstest]
#[case::donor(AtomIdx(2), AtomIdx(3), Some(DativeBondIdx(0)))]
#[case::reverse_rejects(AtomIdx(3), AtomIdx(2), None)]
#[case::unrelated(AtomIdx(0), AtomIdx(3), None)]
fn test_molecule_ast_connecting_dative_bond(
    #[case] donor: AtomIdx,
    #[case] acceptor: AtomIdx,
    #[case] expected: Option<DativeBondIdx>,
) {
    let ast = rich_molecule();
    assert_eq!(ast.connecting_dative_bond(donor, acceptor), expected);
}

#[rstest]
#[case::forward(AtomIdx(0), AtomIdx(3), Some(NoncovalentBondIdx(0)))]
#[case::reverse(AtomIdx(3), AtomIdx(0), Some(NoncovalentBondIdx(0)))]
#[case::unrelated(AtomIdx(0), AtomIdx(1), None)]
fn test_molecule_ast_connecting_noncovalent_bond(
    #[case] a: AtomIdx,
    #[case] b: AtomIdx,
    #[case] expected: Option<NoncovalentBondIdx>,
) {
    let ast = rich_molecule();
    assert_eq!(ast.connecting_noncovalent_bond(a, b), expected);
}

#[rstest]
#[case::donor(AtomIdx(2), vec![DativeBondIdx(0)])]
#[case::acceptor(AtomIdx(3), vec![DativeBondIdx(0)])]
#[case::outside(AtomIdx(0), vec![])]
fn test_molecule_ast_dative_bonds_incident(
    #[case] atom: AtomIdx,
    #[case] expected: Vec<DativeBondIdx>,
) {
    let ast = rich_molecule();
    let inc: Vec<_> = ast.dative_bonds_incident(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomIdx(1), vec![AromaticSystemIdx(0)])]
#[case::outside(AtomIdx(3), vec![])]
fn test_molecule_ast_aromatic_systems_incident(
    #[case] atom: AtomIdx,
    #[case] expected: Vec<AromaticSystemIdx>,
) {
    let ast = rich_molecule();
    let inc: Vec<_> = ast.aromatic_systems_incident(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomIdx(0), vec![MulticenterBondIdx(0)])]
#[case::outside(AtomIdx(3), vec![])]
fn test_molecule_ast_multicenter_bonds_incident(
    #[case] atom: AtomIdx,
    #[case] expected: Vec<MulticenterBondIdx>,
) {
    let ast = rich_molecule();
    let inc: Vec<_> = ast.multicenter_bonds_incident(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::first(AtomIdx(0), vec![NoncovalentBondIdx(0)])]
#[case::second(AtomIdx(3), vec![NoncovalentBondIdx(0)])]
#[case::outside(AtomIdx(1), vec![])]
fn test_molecule_ast_noncovalent_bonds_incident(
    #[case] atom: AtomIdx,
    #[case] expected: Vec<NoncovalentBondIdx>,
) {
    let ast = rich_molecule();
    let inc: Vec<_> = ast.noncovalent_bonds_incident(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::full(vec![AtomIdx(2), AtomIdx(3)], vec![DativeBondIdx(0)])]
#[case::partial_only(vec![AtomIdx(0), AtomIdx(2)], vec![])]
#[case::disjoint(vec![AtomIdx(0), AtomIdx(1)], vec![])]
fn test_molecule_ast_induced_dative_bonds(
    #[case] atoms: Vec<AtomIdx>,
    #[case] expected: Vec<DativeBondIdx>,
) {
    let ast = rich_molecule();
    assert_eq!(ast.induced_dative_bonds(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)], vec![AromaticSystemIdx(0)])]
#[case::partial(vec![AtomIdx(0), AtomIdx(1)], vec![])]
#[case::disjoint(vec![AtomIdx(3)], vec![])]
fn test_molecule_ast_induced_aromatic_systems(
    #[case] atoms: Vec<AtomIdx>,
    #[case] expected: Vec<AromaticSystemIdx>,
) {
    let ast = rich_molecule();
    assert_eq!(ast.induced_aromatic_systems(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)], vec![MulticenterBondIdx(0)])]
#[case::partial(vec![AtomIdx(0), AtomIdx(1)], vec![])]
#[case::disjoint(vec![AtomIdx(3)], vec![])]
fn test_molecule_ast_induced_multicenter_bonds(
    #[case] atoms: Vec<AtomIdx>,
    #[case] expected: Vec<MulticenterBondIdx>,
) {
    let ast = rich_molecule();
    assert_eq!(ast.induced_multicenter_bonds(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomIdx(0), AtomIdx(3)], vec![NoncovalentBondIdx(0)])]
#[case::partial(vec![AtomIdx(0), AtomIdx(1)], vec![])]
#[case::disjoint(vec![AtomIdx(1), AtomIdx(2)], vec![])]
fn test_molecule_ast_induced_noncovalent_bonds(
    #[case] atoms: Vec<AtomIdx>,
    #[case] expected: Vec<NoncovalentBondIdx>,
) {
    let ast = rich_molecule();
    assert_eq!(ast.induced_noncovalent_bonds(&atoms), expected);
}

#[rstest]
#[case::atom_0(AtomIdx(0), Element::C)]
#[case::atom_1(AtomIdx(1), Element::C)]
#[case::atom_2(AtomIdx(2), Element::N)]
#[case::atom_3(AtomIdx(3), Element::O)]
fn test_molecule_ast_atom(#[case] idx: AtomIdx, #[case] element: Element) {
    let ast = rich_molecule();
    let av = ast.atom(idx);
    assert_eq!(av.idx, idx);
    assert_eq!(av.data.element, ElementAst::Lit(element));
}

#[rstest]
fn test_molecule_ast_atoms() {
    let ast = rich_molecule();
    let projected: Vec<(AtomIdx, ElementAst)> = ast
        .atoms()
        .iter()
        .map(|v| (v.idx, v.data.element.clone()))
        .collect();
    assert_eq!(
        projected,
        vec![
            (AtomIdx(0), ElementAst::Lit(Element::C)),
            (AtomIdx(1), ElementAst::Lit(Element::C)),
            (AtomIdx(2), ElementAst::Lit(Element::N)),
            (AtomIdx(3), ElementAst::Lit(Element::O)),
        ]
    );
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
    let auto = ast.automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty);
    assert_eq!(auto.num_orbits(), expected_orbits);
    assert_eq!(auto.atom_count(), ast.atom_count());
}

#[test]
fn test_atom_automorphism_same_orbit() {
    let ast = ring(6);
    let auto = ast.automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty);
    assert!(auto.same_orbit(AtomIdx(0), AtomIdx(3)));
}

#[rstest]
fn test_molecule_ast_subgraph_isomorphisms() {
    let target = ring(6);
    let query = chain(2);
    let mut matches = target.subgraph_isomorphisms(
        &query,
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
    );
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![
            vec![AtomIdx(0), AtomIdx(1)],
            vec![AtomIdx(0), AtomIdx(5)],
            vec![AtomIdx(1), AtomIdx(0)],
            vec![AtomIdx(1), AtomIdx(2)],
            vec![AtomIdx(2), AtomIdx(1)],
            vec![AtomIdx(2), AtomIdx(3)],
            vec![AtomIdx(3), AtomIdx(2)],
            vec![AtomIdx(3), AtomIdx(4)],
            vec![AtomIdx(4), AtomIdx(3)],
            vec![AtomIdx(4), AtomIdx(5)],
            vec![AtomIdx(5), AtomIdx(0)],
            vec![AtomIdx(5), AtomIdx(4)],
        ]
    );
}

#[rstest]
fn test_molecule_ast_subgraph_isomorphisms_at() {
    let target = ring(6);
    let query = chain(2);
    let mut matches = target.subgraph_isomorphisms_at(
        &query,
        (AtomIdx(0), AtomIdx(0)),
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
    );
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![
            vec![AtomIdx(0), AtomIdx(1)],
            vec![AtomIdx(0), AtomIdx(5)],
        ]
    );
}

#[rstest]
fn test_molecule_ast_induced_subgraph() {
    let ast = rich_molecule();
    let sub = ast.induced_subgraph(&[AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
    let atom_elements: Vec<_> = sub
        .ast
        .atoms()
        .iter()
        .map(|v| v.data.element.clone())
        .collect();
    assert_eq!(
        atom_elements,
        vec![
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::N),
        ]
    );
    let bonds: Vec<(AtomIdx, AtomIdx, ValueAst)> = sub
        .ast
        .bonds()
        .iter()
        .map(|v| (v.src, v.tgt, v.data.order.clone()))
        .collect();
    assert_eq!(
        bonds,
        vec![
            (AtomIdx(0), AtomIdx(1), ValueAst::Lit(1)),
            (AtomIdx(1), AtomIdx(2), ValueAst::Lit(2)),
        ]
    );
    assert_eq!(sub.atom_map, vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]);
    assert_eq!(sub.bond_map, vec![BondIdx(0), BondIdx(1)]);
    assert_eq!(sub.aromatic_system_map, vec![AromaticSystemIdx(0)]);
    assert_eq!(sub.multicenter_bond_map, vec![MulticenterBondIdx(0)]);
    assert_eq!(sub.dative_bond_map, Vec::<DativeBondIdx>::new());
    assert_eq!(sub.noncovalent_bond_map, Vec::<NoncovalentBondIdx>::new());
}

#[rstest]
fn test_molecule_ast_induced_subgraph_preserves_dative() {
    let ast = rich_molecule();
    let sub = ast.induced_subgraph(&[AtomIdx(2), AtomIdx(3)]);
    assert_eq!(sub.atom_map, vec![AtomIdx(2), AtomIdx(3)]);
    assert_eq!(sub.dative_bond_map, vec![DativeBondIdx(0)]);
    let dv = sub.ast.dative_bond(DativeBondIdx(0));
    assert_eq!(dv.donor, AtomIdx(0));
    assert_eq!(dv.acceptor, AtomIdx(1));
    assert_eq!(dv.data.direction, DativeDirection::Forward);
}

#[rstest]
fn test_molecule_builder_remove_aromatic_systems() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.remove_aromatic_systems(&[AromaticSystemIdx(0)]);
    let result = b.build();
    assert_eq!(
        result
            .aromatic_systems()
            .ids()
            .collect::<Vec<_>>(),
        Vec::<AromaticSystemIdx>::new()
    );
    assert_eq!(
        result.atoms().iter().map(|v| v.idx).collect::<Vec<_>>(),
        vec![AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)]
    );
    assert_eq!(
        result.bonds().iter().map(|v| v.idx).collect::<Vec<_>>(),
        vec![BondIdx(0), BondIdx(1), BondIdx(2)]
    );
}

#[rstest]
fn test_molecule_builder_remove_dative_bonds() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.remove_dative_bonds(&[DativeBondIdx(0)]);
    let result = b.build();
    assert_eq!(
        result.dative_bonds().ids().collect::<Vec<_>>(),
        Vec::<DativeBondIdx>::new()
    );
}

#[rstest]
fn test_molecule_builder_remove_multicenter_bonds() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.remove_multicenter_bonds(&[MulticenterBondIdx(0)]);
    let result = b.build();
    assert_eq!(
        result.multicenter_bonds().ids().collect::<Vec<_>>(),
        Vec::<MulticenterBondIdx>::new()
    );
}

#[rstest]
fn test_molecule_builder_remove_noncovalent_bonds() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.remove_noncovalent_bonds(&[NoncovalentBondIdx(0)]);
    let result = b.build();
    assert_eq!(
        result.noncovalent_bonds().ids().collect::<Vec<_>>(),
        Vec::<NoncovalentBondIdx>::new()
    );
}

#[rstest]
fn test_molecule_builder_atom_mut() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.atom_mut(AtomIdx(0)).element = ElementAst::Lit(Element::N);
    let result = b.build();
    assert_eq!(result[AtomIdx(0)].element, ElementAst::Lit(Element::N));
    assert_eq!(ast[AtomIdx(0)].element, ElementAst::Lit(Element::C));
}

#[rstest]
fn test_molecule_builder_bond_mut() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.bond_mut(BondIdx(0)).order = ValueAst::Lit(3);
    let result = b.build();
    assert_eq!(result[BondIdx(0)].order, ValueAst::Lit(3));
    assert_eq!(ast[BondIdx(0)].order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_builder_constraints_mut() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    b.constraints_mut()
        .push_atom(AtomIdx(0), AtomConstraint::Degree(ValueAst::Lit(2)));
    let result = b.build();
    assert_eq!(
        result.constraints().atom(AtomIdx(0)),
        &[AtomConstraint::Degree(ValueAst::Lit(2))]
    );
    assert_eq!(
        ast.constraints().atom(AtomIdx(0)),
        &[] as &[AtomConstraint]
    );
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

#[rstest]
fn test_molecule_builder_add_and_remove() {
    let ast = rich_molecule();
    let mut b = ast.edit();
    let new_a = b.add_atom(AtomAst::from_element(Element::Br));
    b.add_bond(AtomIdx(0), new_a, BondAst::from_order(1));
    b.remove_aromatic_systems(&[AromaticSystemIdx(0)]);
    let _remap = b.remove(&[AtomIdx(3)], &[BondIdx(2)]);
    let result = b.build();
    let atoms: Vec<Element> = result
        .atoms()
        .iter()
        .map(|v| match v.data.element {
            ElementAst::Lit(e) => e,
            _ => panic!("non-ground element in builder result"),
        })
        .collect();
    assert_eq!(atoms, vec![Element::C, Element::C, Element::N, Element::Br]);
    let bonds: Vec<(AtomIdx, AtomIdx, ValueAst)> = result
        .bonds()
        .iter()
        .map(|v| (v.src, v.tgt, v.data.order.clone()))
        .collect();
    assert_eq!(
        bonds,
        vec![
            (AtomIdx(0), AtomIdx(1), ValueAst::Lit(1)),
            (AtomIdx(1), AtomIdx(2), ValueAst::Lit(2)),
            (AtomIdx(0), AtomIdx(3), ValueAst::Lit(1)),
        ]
    );
    assert_eq!(
        result.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemIdx>::new()
    );
    assert_eq!(
        result.dative_bonds().ids().collect::<Vec<_>>(),
        Vec::<DativeBondIdx>::new()
    );
    assert_eq!(
        result.noncovalent_bonds().ids().collect::<Vec<_>>(),
        Vec::<NoncovalentBondIdx>::new()
    );
}

#[rstest]
#[case::forward_order(AtomIdx(0), AtomIdx(1), DativeDirection::Forward)]
#[case::reverse_order(AtomIdx(1), AtomIdx(0), DativeDirection::Reverse)]
fn test_molecule_ast_dative_direction_normalization(
    #[case] donor: AtomIdx,
    #[case] acceptor: AtomIdx,
    #[case] expected_direction: DativeDirection,
) {
    let atoms = vec![ground_atom(), ground_atom()];
    let ast = MoleculeAst::new(
        atoms,
        Vec::new(),
        vec![(donor, acceptor, DativeBondAst::new())],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    let view = ast.dative_bond(DativeBondIdx(0));
    assert_eq!(view.donor, donor);
    assert_eq!(view.acceptor, acceptor);
    assert_eq!(view.data.direction, expected_direction);
}

#[rstest]
fn test_molecule_ast_eq_canonical_across_bond_order() {
    let atoms_a = vec![ground_atom(), ground_atom()];
    let atoms_b = vec![ground_atom(), ground_atom()];
    let bond = BondAst {
        order: ValueAst::Lit(1),
        charge: ValueAst::Lit(0),
        spin: SpinStateAst::closed_shell(),
        constraints: Vec::new(),
    };
    let forward = MoleculeAst::new(
        atoms_a,
        vec![(AtomIdx(0), AtomIdx(1), bond.clone())],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    let reverse = MoleculeAst::new(
        atoms_b,
        vec![(AtomIdx(1), AtomIdx(0), bond)],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    assert_eq!(forward, reverse);
}

#[rstest]
fn test_molecule_ast_eq_canonical_across_dative_order() {
    let atoms_a = vec![ground_atom(), ground_atom()];
    let atoms_b = vec![ground_atom(), ground_atom()];
    let forward = MoleculeAst::new(
        atoms_a,
        Vec::new(),
        vec![(AtomIdx(0), AtomIdx(1), DativeBondAst::new())],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    let reverse = MoleculeAst::new(
        atoms_b,
        Vec::new(),
        vec![(AtomIdx(1), AtomIdx(0), DativeBondAst::new())],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    assert_ne!(
        forward, reverse,
        "dative direction is part of identity; reversed donor/acceptor should differ"
    );
}

#[rstest]
fn test_molecule_ast_graph() {
    let ast = rich_molecule();
    let g = ast.graph();
    assert_eq!(g.node_count(), 4);
    assert_eq!(g.edge_count(), 3);
    assert_eq!(
        g.edge_endpoints(EdgeId(0)),
        [NodeId(0), NodeId(1)]
    );
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomIdx(0), AtomIdx(1), AtomIdx(2)]),
    Some(AromaticSystemIdx(0)),
)]
#[case::subset(
    HashSet::from([AtomIdx(0), AtomIdx(1)]),
    None,
)]
#[case::disjoint(
    HashSet::from([AtomIdx(3)]),
    None,
)]
fn test_molecule_ast_connecting_aromatic_system(
    #[case] atoms: HashSet<AtomIdx>,
    #[case] expected: Option<AromaticSystemIdx>,
) {
    let ast = rich_molecule();
    assert_eq!(ast.connecting_aromatic_system(&atoms), expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomIdx(0), AtomIdx(1), AtomIdx(2)]),
    Some(MulticenterBondIdx(0)),
)]
#[case::subset(
    HashSet::from([AtomIdx(0), AtomIdx(1)]),
    None,
)]
fn test_molecule_ast_connecting_multicenter_bond(
    #[case] atoms: HashSet<AtomIdx>,
    #[case] expected: Option<MulticenterBondIdx>,
) {
    let ast = rich_molecule();
    assert_eq!(ast.connecting_multicenter_bond(&atoms), expected);
}

#[rstest]
fn test_molecule_ast_enumerate_maximum_matchings() {
    let ast = ring(4);
    let mut ms: Vec<Vec<(AtomIdx, AtomIdx)>> = ast
        .enumerate_maximum_matchings(MatchingEnumerationAlgorithm::BranchAndBound)
        .into_iter()
        .map(|m| {
            let mut pairs: Vec<_> = (0..ast.atom_count())
                .map(AtomIdx::from)
                .filter_map(|a| m.mate(a).filter(|b| a < *b).map(|b| (a, b)))
                .collect();
            pairs.sort_unstable();
            pairs
        })
        .collect();
    ms.sort_unstable();
    assert_eq!(
        ms,
        vec![
            vec![(AtomIdx(0), AtomIdx(1)), (AtomIdx(2), AtomIdx(3))],
            vec![(AtomIdx(0), AtomIdx(3)), (AtomIdx(1), AtomIdx(2))],
        ]
    );
}

#[rstest]
fn test_molecule_ast_index_atom() {
    let ast = rich_molecule();
    assert_eq!(ast[AtomIdx(2)].element, ElementAst::Lit(Element::N));
}

#[rstest]
fn test_molecule_ast_index_bond() {
    let ast = rich_molecule();
    assert_eq!(ast[BondIdx(1)].order, ValueAst::Lit(2));
}

#[rstest]
fn test_molecule_ast_index_dative_bond() {
    let ast = rich_molecule();
    assert_eq!(ast[DativeBondIdx(0)].direction, DativeDirection::Forward);
}

#[rstest]
fn test_molecule_ast_index_aromatic_system() {
    let ast = rich_molecule();
    assert_eq!(ast[AromaticSystemIdx(0)].electrons, ValueAst::Undetermined);
}

#[rstest]
fn test_molecule_ast_index_multicenter_bond() {
    let ast = rich_molecule();
    assert_eq!(ast[MulticenterBondIdx(0)].electrons, ValueAst::Undetermined);
}

#[rstest]
fn test_molecule_ast_index_noncovalent_bond() {
    let ast = rich_molecule();
    assert_eq!(
        ast[NoncovalentBondIdx(0)].kind,
        NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond)
    );
}

#[rstest]
fn test_molecule_ast_atoms_mut() {
    let mut ast = rich_molecule();
    for a in ast.atoms_mut() {
        a.charge = ValueAst::Lit(1);
    }
    let charges: Vec<ValueAst> = ast.atoms().iter().map(|v| v.data.charge.clone()).collect();
    assert_eq!(
        charges,
        vec![
            ValueAst::Lit(1),
            ValueAst::Lit(1),
            ValueAst::Lit(1),
            ValueAst::Lit(1),
        ]
    );
}

#[rstest]
fn test_molecule_ast_bonds_mut() {
    let mut ast = rich_molecule();
    for b in ast.bonds_mut() {
        b.order = ValueAst::Lit(1);
    }
    let orders: Vec<ValueAst> = ast.bonds().iter().map(|v| v.data.order.clone()).collect();
    assert_eq!(
        orders,
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(1)]
    );
}

#[rstest]
fn test_molecule_ast_dative_bond_mut() {
    let mut ast = rich_molecule();
    ast.dative_bond_mut(DativeBondIdx(0))
        .constraints
        .push(DativeBondConstraint::RingSize(
            ValueAst::Lit(6),
        ));
    assert_eq!(
        ast[DativeBondIdx(0)].constraints,
        vec![DativeBondConstraint::RingSize(
            ValueAst::Lit(6)
        )]
    );
}

#[rstest]
fn test_molecule_ast_aromatic_system_mut() {
    let mut ast = rich_molecule();
    ast.aromatic_system_mut(AromaticSystemIdx(0)).electrons = ValueAst::Lit(6);
    assert_eq!(ast[AromaticSystemIdx(0)].electrons, ValueAst::Lit(6));
}

#[rstest]
fn test_molecule_ast_aromatic_systems_mut() {
    let mut ast = rich_molecule();
    for a in ast.aromatic_systems_mut() {
        a.electrons = ValueAst::Lit(6);
    }
    let electrons: Vec<ValueAst> = ast
        .aromatic_systems()
        .iter()
        .map(|v| v.data.electrons.clone())
        .collect();
    assert_eq!(electrons, vec![ValueAst::Lit(6)]);
}

#[rstest]
fn test_molecule_ast_multicenter_bond_mut() {
    let mut ast = rich_molecule();
    ast.multicenter_bond_mut(MulticenterBondIdx(0)).electrons = ValueAst::Lit(2);
    assert_eq!(ast[MulticenterBondIdx(0)].electrons, ValueAst::Lit(2));
}

#[rstest]
fn test_molecule_ast_multicenter_bonds_mut() {
    let mut ast = rich_molecule();
    for m in ast.multicenter_bonds_mut() {
        m.electrons = ValueAst::Lit(2);
    }
    let electrons: Vec<ValueAst> = ast
        .multicenter_bonds()
        .iter()
        .map(|v| v.data.electrons.clone())
        .collect();
    assert_eq!(electrons, vec![ValueAst::Lit(2)]);
}

#[rstest]
fn test_molecule_ast_noncovalent_bond_mut() {
    let mut ast = rich_molecule();
    ast.noncovalent_bond_mut(NoncovalentBondIdx(0)).kind =
        NoncovalentKindAst::Lit(NoncovalentKind::Ionic);
    assert_eq!(
        ast[NoncovalentBondIdx(0)].kind,
        NoncovalentKindAst::Lit(NoncovalentKind::Ionic)
    );
}
