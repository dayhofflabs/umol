use std::collections::HashSet;

use pretty_assertions::assert_eq;
use rstest::*;
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, MatchingEnumerationAlgorithm, MaxIndependentSetAlgorithm,
    MaxMatchingAlgorithm, NodeId, ShortestCycleAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_shared::element::Element;

use super::super::aromatic::AromaticSystemAst;
use super::super::atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
use super::super::bond::BondAst;
use super::super::constraint::{
    AtomConstraint, AtomConstraints, BondConstraint, BondConstraints, Constraint, Constraints,
    DativeBondConstraint, DativeBondConstraints, MoleculeConstraint, RelationalConstraint,
    SubPatternAnchor,
};
use super::super::dative::DativeBondAst;
use super::super::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};
use super::super::multicenter::MulticenterBondAst;
use super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst};
use super::super::rings::RingFamily;
use super::super::spin::SpinStateAst;
use super::super::value::{Expr, ValueAst};
use super::MoleculeAst;
use crate::{mol, mol_zeroed};

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
    out.push(c);
    out
}

#[rstest]
fn test_molecule_ast_new() {
    let m = MoleculeAst::new();
    assert_eq!(m.atom_count(), 0);
    assert_eq!(m.bond_count(), 0);
    assert_eq!(m.dative_bond_count(), 0);
    assert_eq!(m.aromatic_system_count(), 0);
    assert_eq!(m.multicenter_bond_count(), 0);
    assert_eq!(m.noncovalent_bond_count(), 0);
    assert_eq!(m.constraints().len(), 0);
}

#[rstest]
fn test_molecule_ast_default_equals_new() {
    assert_eq!(MoleculeAst::default(), MoleculeAst::new());
}

#[rstest]
fn test_molecule_ast_from_atoms_and_bonds() {
    let atoms = vec![
        AtomAst::from_element(Element::C),
        AtomAst::from_element(Element::O),
    ];
    let bonds = vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))];
    let m = MoleculeAst::from_atoms_and_bonds(atoms, bonds);
    assert_eq!(m.atom_count(), 2);
    assert_eq!(m.bond_count(), 1);
    assert_eq!(m.dative_bond_count(), 0);
    assert_eq!(m.aromatic_system_count(), 0);
    assert_eq!(m.multicenter_bond_count(), 0);
    assert_eq!(m.noncovalent_bond_count(), 0);
    assert_eq!(m.atom(AtomIdx(0)).data.element, ElementAst::Lit(Element::C));
    assert_eq!(m.atom(AtomIdx(1)).data.element, ElementAst::Lit(Element::O));
    assert_eq!(m.bond(BondIdx(0)).data.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_builder() {
    assert_eq!(MoleculeAst::builder().build(), MoleculeAst::new());
}

#[rstest]
#[case::empty(MoleculeAst::default(), true)]
#[case::ground_atom(
    mol_zeroed!(r#"{:atoms ["C #h4"] :bonds []}"#),
    true,
)]
#[case::wildcard_element(
    mol!(r#"{:atoms ["*"] :bonds []}"#),
    false,
)]
#[case::wildcard_bond(
    mol!(r#"{:atoms ["C" "O"] :bonds [[0 1 "*"]]}"#),
    false,
)]
#[case::ground_atom_with_undetermined_constraint(
    MoleculeAst::from_parts(
        vec![ground_atom()],
        vec![], vec![], vec![], vec![], vec![],
        constraints_with_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![]),
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
fn test_molecule_ast_neighbors(#[case] atom: AtomIdx, #[case] expected: Vec<(AtomIdx, BondIdx)>) {
    let ast = MoleculeAst::from_parts(
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
    let nbrs: Vec<(AtomIdx, BondIdx)> = ast.neighbors(atom).map(|n| (n.atom, n.bond)).collect();
    assert_eq!(nbrs, expected);
}

#[rstest]
fn test_molecule_builder_add_aromatic_system() {
    let ast = MoleculeAst::from_parts(
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
        new_ast.aromatic_systems().ids().collect::<Vec<_>>(),
        vec![AromaticSystemIdx(0)]
    );
    assert_eq!(
        ast.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemIdx>::new()
    );
}

#[fixture]
fn rich_molecule() -> MoleculeAst {
    MoleculeAst::from_parts(
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
        vec![(vec![AtomIdx(2)], AtomIdx(3), DativeBondAst::from_order(1))],
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
            NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        Constraints::default(),
    )
}

#[rstest]
#[case::c_c(BondIdx(0), AtomIdx(0), AtomIdx(1), ValueAst::Lit(1))]
#[case::c_n(BondIdx(1), AtomIdx(1), AtomIdx(2), ValueAst::Lit(2))]
#[case::n_o(BondIdx(2), AtomIdx(2), AtomIdx(3), ValueAst::Lit(1))]
fn test_molecule_ast_bond(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] idx: BondIdx,
    #[case] src: AtomIdx,
    #[case] tgt: AtomIdx,
    #[case] order: ValueAst,
) {
    let bv = ast.bond(idx);
    assert_eq!(bv.idx, idx);
    assert_eq!(bv.src, src);
    assert_eq!(bv.tgt, tgt);
    assert_eq!(bv.data.order, order);
}

#[rstest]
fn test_molecule_ast_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
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
fn test_molecule_ast_dative_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let dv = ast.dative_bond(DativeBondIdx(0));
    assert_eq!(dv.idx, DativeBondIdx(0));
    assert_eq!(dv.acceptor, AtomIdx(3));
    assert_eq!(dv.donors().collect::<Vec<_>>(), vec![AtomIdx(2)]);
    assert_eq!(dv.atoms().collect::<Vec<_>>(), vec![AtomIdx(2), AtomIdx(3)]);
    assert_eq!(dv.data.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_dative_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(DativeBondIdx, Vec<AtomIdx>, AtomIdx)> = ast
        .dative_bonds()
        .iter()
        .map(|v| (v.idx, v.donors().collect(), v.acceptor))
        .collect();
    assert_eq!(
        projected,
        vec![(DativeBondIdx(0), vec![AtomIdx(2)], AtomIdx(3))]
    );
}

#[rstest]
fn test_molecule_ast_aromatic_system(#[from(rich_molecule)] ast: MoleculeAst) {
    let av = ast.aromatic_system(AromaticSystemIdx(0));
    assert_eq!(av.idx, AromaticSystemIdx(0));
    assert_eq!(
        av.atoms().collect::<Vec<_>>(),
        vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]
    );
    assert_eq!(av.bonds().collect::<Vec<_>>(), vec![BondIdx(0), BondIdx(1)]);
}

#[rstest]
fn test_molecule_ast_aromatic_systems(#[from(rich_molecule)] ast: MoleculeAst) {
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
fn test_molecule_ast_multicenter_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mv = ast.multicenter_bond(MulticenterBondIdx(0));
    assert_eq!(mv.idx, MulticenterBondIdx(0));
    assert_eq!(
        mv.atoms().collect::<Vec<_>>(),
        vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)]
    );
}

#[rstest]
fn test_molecule_ast_multicenter_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
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
fn test_molecule_ast_noncovalent_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let nv = ast.noncovalent_bond(NoncovalentBondIdx(0));
    assert_eq!(nv.idx, NoncovalentBondIdx(0));
    assert_eq!(nv.atoms, [AtomIdx(0), AtomIdx(3)]);
}

#[rstest]
fn test_molecule_ast_noncovalent_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
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
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomIdx,
    #[case] b: AtomIdx,
    #[case] expected: Option<BondIdx>,
) {
    assert_eq!(ast.connecting_bond(a, b), expected);
}

#[rstest]
#[case::full_match(HashSet::from([AtomIdx(2), AtomIdx(3)]), Some(DativeBondIdx(0)))]
#[case::subset(HashSet::from([AtomIdx(2)]), None)]
#[case::disjoint(HashSet::from([AtomIdx(0), AtomIdx(1)]), None)]
fn test_molecule_ast_connecting_dative_bond(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomIdx>,
    #[case] expected: Option<DativeBondIdx>,
) {
    assert_eq!(ast.connecting_dative_bond(&atoms), expected);
}

#[rstest]
#[case::forward(AtomIdx(0), AtomIdx(3), Some(NoncovalentBondIdx(0)))]
#[case::reverse(AtomIdx(3), AtomIdx(0), Some(NoncovalentBondIdx(0)))]
#[case::unrelated(AtomIdx(0), AtomIdx(1), None)]
fn test_molecule_ast_connecting_noncovalent_bond(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomIdx,
    #[case] b: AtomIdx,
    #[case] expected: Option<NoncovalentBondIdx>,
) {
    assert_eq!(ast.connecting_noncovalent_bond(a, b), expected);
}

#[rstest]
#[case::donor(AtomIdx(2), vec![DativeBondIdx(0)])]
#[case::acceptor(AtomIdx(3), vec![DativeBondIdx(0)])]
#[case::outside(AtomIdx(0), vec![])]
fn test_molecule_ast_dative_bonds_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomIdx,
    #[case] expected: Vec<DativeBondIdx>,
) {
    let inc: Vec<_> = ast.dative_bonds_incident(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomIdx(1), vec![AromaticSystemIdx(0)])]
#[case::outside(AtomIdx(3), vec![])]
fn test_molecule_ast_aromatic_systems_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomIdx,
    #[case] expected: Vec<AromaticSystemIdx>,
) {
    let inc: Vec<_> = ast.aromatic_systems_incident(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomIdx(0), vec![MulticenterBondIdx(0)])]
#[case::outside(AtomIdx(3), vec![])]
fn test_molecule_ast_multicenter_bonds_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomIdx,
    #[case] expected: Vec<MulticenterBondIdx>,
) {
    let inc: Vec<_> = ast.multicenter_bonds_incident(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::first(AtomIdx(0), vec![NoncovalentBondIdx(0)])]
#[case::second(AtomIdx(3), vec![NoncovalentBondIdx(0)])]
#[case::outside(AtomIdx(1), vec![])]
fn test_molecule_ast_noncovalent_bonds_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomIdx,
    #[case] expected: Vec<NoncovalentBondIdx>,
) {
    let inc: Vec<_> = ast.noncovalent_bonds_incident(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::full(vec![AtomIdx(2), AtomIdx(3)], vec![DativeBondIdx(0)])]
#[case::partial_only(vec![AtomIdx(0), AtomIdx(2)], vec![])]
#[case::disjoint(vec![AtomIdx(0), AtomIdx(1)], vec![])]
fn test_molecule_ast_induced_dative_bonds(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomIdx>,
    #[case] expected: Vec<DativeBondIdx>,
) {
    assert_eq!(ast.induced_dative_bonds(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)], vec![AromaticSystemIdx(0)])]
#[case::partial(vec![AtomIdx(0), AtomIdx(1)], vec![])]
#[case::disjoint(vec![AtomIdx(3)], vec![])]
fn test_molecule_ast_induced_aromatic_systems(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomIdx>,
    #[case] expected: Vec<AromaticSystemIdx>,
) {
    assert_eq!(ast.induced_aromatic_systems(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomIdx(0), AtomIdx(1), AtomIdx(2)], vec![MulticenterBondIdx(0)])]
#[case::partial(vec![AtomIdx(0), AtomIdx(1)], vec![])]
#[case::disjoint(vec![AtomIdx(3)], vec![])]
fn test_molecule_ast_induced_multicenter_bonds(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomIdx>,
    #[case] expected: Vec<MulticenterBondIdx>,
) {
    assert_eq!(ast.induced_multicenter_bonds(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomIdx(0), AtomIdx(3)], vec![NoncovalentBondIdx(0)])]
#[case::partial(vec![AtomIdx(0), AtomIdx(1)], vec![])]
#[case::disjoint(vec![AtomIdx(1), AtomIdx(2)], vec![])]
fn test_molecule_ast_induced_noncovalent_bonds(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomIdx>,
    #[case] expected: Vec<NoncovalentBondIdx>,
) {
    assert_eq!(ast.induced_noncovalent_bonds(&atoms), expected);
}

#[rstest]
#[case::atom_0(AtomIdx(0), Element::C)]
#[case::atom_1(AtomIdx(1), Element::C)]
#[case::atom_2(AtomIdx(2), Element::N)]
#[case::atom_3(AtomIdx(3), Element::O)]
fn test_molecule_ast_atom(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] idx: AtomIdx,
    #[case] element: Element,
) {
    let av = ast.atom(idx);
    assert_eq!(av.idx, idx);
    assert_eq!(av.data.element, ElementAst::Lit(element));
}

#[rstest]
fn test_molecule_ast_atoms(#[from(rich_molecule)] ast: MoleculeAst) {
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
    let ast = MoleculeAst::from_parts(
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
    MoleculeAst::from_atoms_and_bonds(
        atoms,
        bonds,
    )
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
    MoleculeAst::from_atoms_and_bonds(
        atoms,
        bonds,
    )
}

fn two_components() -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); 4];
    let bonds = vec![
        (AtomIdx(0), AtomIdx(1), BondAst::from_order(1)),
        (AtomIdx(2), AtomIdx(3), BondAst::from_order(1)),
    ];
    MoleculeAst::from_parts(
        atoms,
        bonds,
        vec![],
        vec![],
        vec![],
        vec![],
        Constraints::default(),
    )
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
fn test_molecule_ast_biconnected_components(#[case] ast: MoleculeAst, #[case] expected: usize) {
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
fn test_molecule_ast_maximum_independent_set(#[case] ast: MoleculeAst, #[case] expected: usize) {
    let mis = ast.maximum_independent_set(MaxIndependentSetAlgorithm::BranchAndBound);
    assert_eq!(mis.len(), expected);
}

#[rstest]
#[case::chain_4(chain(4), 2)]
#[case::ring_6(ring(6), 3)]
#[case::single(chain(1), 0)]
fn test_molecule_ast_maximum_matching(#[case] ast: MoleculeAst, #[case] expected_size: usize) {
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
fn test_molecule_ast_automorphisms(#[case] ast: MoleculeAst, #[case] expected_orbits: usize) {
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
        vec![vec![AtomIdx(0), AtomIdx(1)], vec![AtomIdx(0), AtomIdx(5)],]
    );
}

#[rstest]
fn test_molecule_ast_induced_subgraph(#[from(rich_molecule)] ast: MoleculeAst) {
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
fn test_molecule_ast_induced_subgraph_preserves_dative(#[from(rich_molecule)] ast: MoleculeAst) {
    let sub = ast.induced_subgraph(&[AtomIdx(2), AtomIdx(3)]);
    assert_eq!(sub.atom_map, vec![AtomIdx(2), AtomIdx(3)]);
    assert_eq!(sub.dative_bond_map, vec![DativeBondIdx(0)]);
    let dv = sub.ast.dative_bond(DativeBondIdx(0));
    assert_eq!(dv.acceptor, AtomIdx(1));
    assert_eq!(dv.donors().collect::<Vec<_>>(), vec![AtomIdx(0)]);
    assert_eq!(dv.data.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_builder_remove_aromatic_systems(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_aromatic_systems(&[AromaticSystemIdx(0)]);
    let result = b.build();
    assert_eq!(
        result.aromatic_systems().ids().collect::<Vec<_>>(),
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
fn test_molecule_builder_remove_dative_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_dative_bonds(&[DativeBondIdx(0)]);
    let result = b.build();
    assert_eq!(
        result.dative_bonds().ids().collect::<Vec<_>>(),
        Vec::<DativeBondIdx>::new()
    );
}

#[rstest]
fn test_molecule_builder_remove_multicenter_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_multicenter_bonds(&[MulticenterBondIdx(0)]);
    let result = b.build();
    assert_eq!(
        result.multicenter_bonds().ids().collect::<Vec<_>>(),
        Vec::<MulticenterBondIdx>::new()
    );
}

#[rstest]
fn test_molecule_builder_remove_noncovalent_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_noncovalent_bonds(&[NoncovalentBondIdx(0)]);
    let result = b.build();
    assert_eq!(
        result.noncovalent_bonds().ids().collect::<Vec<_>>(),
        Vec::<NoncovalentBondIdx>::new()
    );
}

#[rstest]
fn test_molecule_builder_atom_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.atom_mut(AtomIdx(0)).element = ElementAst::Lit(Element::N);
    let result = b.build();
    assert_eq!(result[AtomIdx(0)].element, ElementAst::Lit(Element::N));
    assert_eq!(ast[AtomIdx(0)].element, ElementAst::Lit(Element::C));
}

#[rstest]
fn test_molecule_builder_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.bond_mut(BondIdx(0)).order = ValueAst::Lit(3);
    let result = b.build();
    assert_eq!(result[BondIdx(0)].order, ValueAst::Lit(3));
    assert_eq!(ast[BondIdx(0)].order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_builder_atom_constraint_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.atom_mut(AtomIdx(0))
        .constraints
        .add(AtomConstraint::Degree(ValueAst::Lit(2)));
    let result = b.build();
    assert_eq!(
        result[AtomIdx(0)].constraints,
        AtomConstraints::from_iter([AtomConstraint::Degree(ValueAst::Lit(2))])
    );
    assert!(ast[AtomIdx(0)].constraints.is_empty());
}

#[rstest]
fn test_molecule_builder_add_dative_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let id = b.add_dative_bond(vec![AtomIdx(1)], AtomIdx(0), DativeBondAst::from_order(1));
    let result = b.build();
    assert_eq!(id, DativeBondIdx(1));
    let view = result.dative_bond(id);
    assert_eq!(view.acceptor, AtomIdx(0));
    assert_eq!(view.donors().collect::<Vec<_>>(), vec![AtomIdx(1)]);
    // Participants are sorted by NodeId; acceptor=0 lands at slot 0.
    assert_eq!(view.data.acceptor_slot, 0);
}

#[rstest]
fn test_molecule_builder_add_multicenter_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let id = b.add_multicenter_bond(
        vec![AtomIdx(1), AtomIdx(2), AtomIdx(3)],
        MulticenterBondAst::default(),
    );
    let result = b.build();
    assert_eq!(id, MulticenterBondIdx(1));
    let atoms: Vec<AtomIdx> = result.multicenter_bond(id).atoms().collect();
    assert_eq!(atoms, vec![AtomIdx(1), AtomIdx(2), AtomIdx(3)]);
}

#[rstest]
fn test_molecule_builder_add_noncovalent_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let id = b.add_noncovalent_bond(
        [AtomIdx(1), AtomIdx(2)],
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
    );
    let result = b.build();
    assert_eq!(id, NoncovalentBondIdx(1));
    let view = result.noncovalent_bond(id);
    assert_eq!(view.atoms, [AtomIdx(1), AtomIdx(2)]);
}

#[rstest]
fn test_molecule_builder_push_constraint_and_constraints_mut(
    #[from(rich_molecule)] ast: MoleculeAst,
) {
    let mut b = ast.edit();
    b.push_constraint(Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomIdx(0), AtomIdx(1)]),
    }));
    b.constraints_mut()
        .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomIdx(0)]),
            sum: ValueAst::Lit(0),
        }));
    let result = b.build();
    assert_eq!(result.constraints().len(), 2);
}

#[rstest]
fn test_molecule_builder_dative_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.dative_bond_mut(DativeBondIdx(0))
        .constraints
        .add(DativeBondConstraint::RingSize(ValueAst::Lit(5)));
    let result = b.build();
    assert!(!result[DativeBondIdx(0)].constraints.is_empty());
    assert!(ast[DativeBondIdx(0)].constraints.is_empty());
}

#[rstest]
fn test_molecule_builder_aromatic_system_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.aromatic_system_mut(AromaticSystemIdx(0)).charge = ValueAst::Lit(0);
    let result = b.build();
    assert_eq!(result[AromaticSystemIdx(0)].charge, ValueAst::Lit(0));
}

#[rstest]
fn test_molecule_builder_multicenter_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.multicenter_bond_mut(MulticenterBondIdx(0)).electrons =
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(0)];
    let result = b.build();
    assert_eq!(
        result[MulticenterBondIdx(0)].electrons,
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(0)],
    );
}

#[rstest]
fn test_molecule_builder_noncovalent_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.noncovalent_bond_mut(NoncovalentBondIdx(0)).kind =
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic);
    let result = b.build();
    assert_eq!(
        result[NoncovalentBondIdx(0)].kind,
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic),
    );
}

#[rstest]
fn test_molecule_builder_remove_empty_is_noop(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_dative_bonds(&[]);
    b.remove_aromatic_systems(&[]);
    b.remove_multicenter_bonds(&[]);
    b.remove_noncovalent_bonds(&[]);
    let result = b.build();
    assert_eq!(result.dative_bonds().count(), 1);
    assert_eq!(result.aromatic_systems().count(), 1);
    assert_eq!(result.multicenter_bonds().count(), 1);
    assert_eq!(result.noncovalent_bonds().count(), 1);
}

#[rstest]
#[case::hexagon(ring(6), 6, 1)]
#[case::hexagon_cutoff(ring(6), 5, 0)]
#[case::chain(chain(5), 10, 0)]
#[case::empty(MoleculeAst::default(), 10, 0)]
fn test_molecule_ast_rings(
    #[case] mut ast: MoleculeAst,
    #[case] max_ring_size: usize,
    #[case] expected: usize,
) {
    let rs = ast.rings(RingFamily::Simple, max_ring_size);
    assert_eq!(rs.count(), expected);
}

#[test]
fn test_molecule_ast_enumerate_rings_atom_filter() {
    let ast = ring(6);
    let rs = ast.enumerate_rings(RingFamily::Simple, 10, |a| a.0 < 3);
    assert_eq!(rs.count(), 0);
}

#[test]
fn test_molecule_ast_rings_induced() {
    let mut ast = mol!(r#"{
        :atoms ["C" "C" "C" "C"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 2 "1"] [1 3 "1"] [2 3 "1"]]
    }"#);
    let simple_count = ast.rings(RingFamily::Simple, 4).count();
    let induced_count = ast.rings(RingFamily::Induced, 4).count();
    assert_eq!(simple_count, 4);
    assert_eq!(induced_count, 4);
}

#[test]
fn test_molecule_ast_rings_induced_naphthalene() {
    let mut ast = mol!(r#"{
        :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
        :bonds [
            [0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
            [3 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]
        ]
    }"#);
    let simple_count = ast.rings(RingFamily::Simple, 10).count();
    assert_eq!(simple_count, 2);
    let induced_count = ast.rings(RingFamily::Induced, 10).count();
    assert_eq!(induced_count, 2);
}

#[test]
fn test_rings_membership() {
    let mut ast = ring(6);
    let rs = ast.rings(RingFamily::Simple, 6);
    assert!(rs.contains_atom(AtomIdx(0)));
    assert!(rs.contains_bond(BondIdx(0)));
    assert_eq!(rs.atom_smallest_ring_size(AtomIdx(0)), Some(6));
}

#[rstest]
fn test_molecule_builder_add_and_remove(#[from(rich_molecule)] ast: MoleculeAst) {
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
#[case::donor_below_acceptor(AtomIdx(0), AtomIdx(1), 1)]
#[case::donor_above_acceptor(AtomIdx(1), AtomIdx(0), 0)]
fn test_molecule_ast_dative_acceptor_slot(
    #[case] donor: AtomIdx,
    #[case] acceptor: AtomIdx,
    #[case] expected_slot: u8,
) {
    let atoms = vec![ground_atom(), ground_atom()];
    let ast = MoleculeAst::from_parts(
        atoms,
        Vec::new(),
        vec![(vec![donor], acceptor, DativeBondAst::from_order(1))],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    let view = ast.dative_bond(DativeBondIdx(0));
    assert_eq!(view.acceptor, acceptor);
    assert_eq!(view.donors().collect::<Vec<_>>(), vec![donor]);
    assert_eq!(view.data.acceptor_slot, expected_slot);
}

#[rstest]
fn test_molecule_ast_eq_canonical_across_bond_order() {
    let atoms_a = vec![ground_atom(), ground_atom()];
    let atoms_b = vec![ground_atom(), ground_atom()];
    let bond = BondAst {
        order: ValueAst::Lit(1),
        charge: ValueAst::Lit(0),
        spin: SpinStateAst::closed_shell(),
        constraints: BondConstraints::new(),
    };
    let forward = MoleculeAst::from_parts(
        atoms_a,
        vec![(AtomIdx(0), AtomIdx(1), bond.clone())],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    let reverse = MoleculeAst::from_parts(
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
    let forward = MoleculeAst::from_parts(
        atoms_a,
        Vec::new(),
        vec![(vec![AtomIdx(0)], AtomIdx(1), DativeBondAst::from_order(1))],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    let reverse = MoleculeAst::from_parts(
        atoms_b,
        Vec::new(),
        vec![(vec![AtomIdx(1)], AtomIdx(0), DativeBondAst::from_order(1))],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    assert_ne!(
        forward, reverse,
        "acceptor identity is part of dative bond; swapping donor/acceptor should differ"
    );
}

#[rstest]
fn test_molecule_ast_graph(#[from(rich_molecule)] ast: MoleculeAst) {
    let g = ast.graph();
    assert_eq!(g.node_count(), 4);
    assert_eq!(g.edge_count(), 3);
    assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
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
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomIdx>,
    #[case] expected: Option<AromaticSystemIdx>,
) {
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
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomIdx>,
    #[case] expected: Option<MulticenterBondIdx>,
) {
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
fn test_molecule_ast_index_atom(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(ast[AtomIdx(2)].element, ElementAst::Lit(Element::N));
}

#[rstest]
fn test_molecule_ast_index_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(ast[BondIdx(1)].order, ValueAst::Lit(2));
}

#[rstest]
fn test_molecule_ast_index_dative_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(ast[DativeBondIdx(0)].order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_index_aromatic_system(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast[AromaticSystemIdx(0)].electrons.is_empty());
}

#[rstest]
fn test_molecule_ast_index_multicenter_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast[MulticenterBondIdx(0)].electrons.is_empty());
}

#[rstest]
fn test_molecule_ast_index_noncovalent_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(
        ast[NoncovalentBondIdx(0)].kind,
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)
    );
}

#[rstest]
fn test_molecule_ast_atoms_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
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
fn test_molecule_ast_bonds_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
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
fn test_molecule_ast_dative_bond_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.dative_bond_mut(DativeBondIdx(0))
        .constraints
        .add(DativeBondConstraint::RingSize(ValueAst::Lit(6)));
    assert_eq!(
        ast[DativeBondIdx(0)].constraints,
        DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(6))])
    );
}

#[rstest]
fn test_molecule_ast_aromatic_system_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.aromatic_system_mut(AromaticSystemIdx(0)).electrons = vec![ValueAst::Lit(1); 3];
    assert_eq!(
        ast[AromaticSystemIdx(0)].electrons,
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(1)],
    );
}

#[rstest]
fn test_molecule_ast_aromatic_systems_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    for a in ast.aromatic_systems_mut() {
        a.electrons = vec![ValueAst::Lit(1); 3];
    }
    let electrons: Vec<Vec<ValueAst>> = ast
        .aromatic_systems()
        .iter()
        .map(|v| v.data.electrons.clone())
        .collect();
    assert_eq!(electrons, vec![vec![ValueAst::Lit(1); 3]]);
}

#[rstest]
fn test_molecule_ast_multicenter_bond_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.multicenter_bond_mut(MulticenterBondIdx(0)).electrons =
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(0)];
    assert_eq!(
        ast[MulticenterBondIdx(0)].electrons,
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(0)],
    );
}

#[rstest]
fn test_molecule_ast_multicenter_bonds_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    for m in ast.multicenter_bonds_mut() {
        m.electrons = vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(0)];
    }
    let electrons: Vec<Vec<ValueAst>> = ast
        .multicenter_bonds()
        .iter()
        .map(|v| v.data.electrons.clone())
        .collect();
    assert_eq!(
        electrons,
        vec![vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(0)]],
    );
}

#[rstest]
fn test_molecule_ast_noncovalent_bond_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.noncovalent_bond_mut(NoncovalentBondIdx(0)).kind =
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic);
    assert_eq!(
        ast[NoncovalentBondIdx(0)].kind,
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic)
    );
}

// -- lift_constraints / inline_constraints ---------------------

/// Set-equality assertion for constraint vecs: order is unspecified, so the
/// test compares as multisets via sort+eq.
fn assert_same_constraints(a: &Constraints, b: &Constraints) {
    let mut x: Vec<&Constraint> = a.iter().collect();
    let mut y: Vec<&Constraint> = b.iter().collect();
    x.sort_by_key(|c| format!("{c:?}"));
    y.sort_by_key(|c| format!("{c:?}"));
    assert_eq!(x, y);
}

#[rstest]
fn test_molecule_ast_lift_constraints_empty() {
    let mut ast = MoleculeAst::default();
    ast.lift_constraints();
    assert!(ast.constraints().is_empty());
}

#[rstest]
fn test_molecule_ast_lift_constraints_drains_inline_stores(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.atom_mut(AtomIdx(0))
        .data
        .constraints
        .add(AtomConstraint::Valence(ValueAst::Lit(4)));
    ast.atom_mut(AtomIdx(2))
        .data
        .constraints
        .add(AtomConstraint::Degree(ValueAst::Lit(3)));
    ast.bond_mut(BondIdx(0))
        .data
        .constraints
        .add(BondConstraint::Aromatic);
    ast.dative_bond_mut(DativeBondIdx(0))
        .constraints
        .add(DativeBondConstraint::RingCount(ValueAst::Lit(1)));

    ast.lift_constraints();

    assert!(ast[AtomIdx(0)].constraints.is_empty());
    assert!(ast[AtomIdx(2)].constraints.is_empty());
    assert!(ast[BondIdx(0)].constraints.is_empty());
    assert!(ast[DativeBondIdx(0)].constraints.is_empty());

    let mut expected = Constraints::new();
    expected.push(Constraint::Atom(
        AtomIdx(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));
    expected.push(Constraint::Atom(
        AtomIdx(2),
        AtomConstraint::Degree(ValueAst::Lit(3)),
    ));
    expected.push(Constraint::Bond(BondIdx(0), BondConstraint::Aromatic));
    expected.push(Constraint::DativeBond(
        DativeBondIdx(0),
        DativeBondConstraint::RingCount(ValueAst::Lit(1)),
    ));
    assert_same_constraints(ast.constraints(), &expected);
}

#[rstest]
fn test_molecule_ast_lift_constraints_appends_to_existing(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    let prior = Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemIdx(0),
        atom: AtomIdx(0),
    });
    ast.constraints_mut().push(prior.clone());
    ast.atom_mut(AtomIdx(0))
        .data
        .constraints
        .add(AtomConstraint::Valence(ValueAst::Lit(4)));

    ast.lift_constraints();

    let mut expected = Constraints::new();
    expected.push(prior);
    expected.push(Constraint::Atom(
        AtomIdx(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));
    assert_same_constraints(ast.constraints(), &expected);
}

#[rstest]
fn test_molecule_ast_inline_constraints_drains_top_level_leaves(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.constraints_mut().push(Constraint::Atom(
        AtomIdx(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));
    ast.constraints_mut()
        .push(Constraint::Bond(BondIdx(0), BondConstraint::Aromatic));
    ast.constraints_mut().push(Constraint::DativeBond(
        DativeBondIdx(0),
        DativeBondConstraint::RingSize(ValueAst::Lit(5)),
    ));

    ast.inline_constraints();

    assert!(ast.constraints().is_empty());
    assert_eq!(
        ast[AtomIdx(0)].constraints,
        AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Lit(4))])
    );
    assert_eq!(
        ast[BondIdx(0)].constraints,
        BondConstraints::from_iter([BondConstraint::Aromatic])
    );
    assert_eq!(
        ast[DativeBondIdx(0)].constraints,
        DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(5))])
    );
}

#[rstest]
fn test_molecule_ast_inline_constraints_last_wins_on_collision(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.constraints_mut().push(Constraint::Atom(
        AtomIdx(0),
        AtomConstraint::Valence(ValueAst::Lit(3)),
    ));
    ast.constraints_mut().push(Constraint::Atom(
        AtomIdx(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));

    ast.inline_constraints();

    // Only one Valence survives; with two competing inserts of the same kind,
    // exactly one wins (which one is unspecified). Verify count and kind.
    assert_eq!(ast[AtomIdx(0)].constraints.len(), 1);
    let v = ast[AtomIdx(0)].constraints.iter().next().unwrap().clone();
    assert!(matches!(v, AtomConstraint::Valence(_)));
}

#[rstest]
fn test_molecule_ast_inline_constraints_skips_combinator_nested(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    let leaf = Constraint::Atom(AtomIdx(0), AtomConstraint::Valence(ValueAst::Lit(4)));
    let nested = Constraint::And(vec![
        leaf.clone(),
        Constraint::Bond(BondIdx(0), BondConstraint::Aromatic),
    ]);
    ast.constraints_mut().push(nested.clone());

    ast.inline_constraints();

    let mut expected = Constraints::new();
    expected.push(nested);
    assert_same_constraints(ast.constraints(), &expected);
    assert!(ast[AtomIdx(0)].constraints.is_empty());
    assert!(ast[BondIdx(0)].constraints.is_empty());
}

#[rstest]
fn test_molecule_ast_inline_constraints_skips_relational_and_molecule(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    let rel = Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemIdx(0),
        atom: AtomIdx(0),
    });
    let mol = Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomIdx(0), AtomIdx(1)]),
    });
    ast.constraints_mut().push(rel.clone());
    ast.constraints_mut().push(mol.clone());
    ast.constraints_mut().push(Constraint::Atom(
        AtomIdx(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));

    ast.inline_constraints();

    let mut expected = Constraints::new();
    expected.push(rel);
    expected.push(mol);
    assert_same_constraints(ast.constraints(), &expected);
    assert_eq!(
        ast[AtomIdx(0)].constraints,
        AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Lit(4))])
    );
}

#[rstest]
fn test_molecule_ast_lift_then_inline_roundtrips_inline_state(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.atom_mut(AtomIdx(0))
        .data
        .constraints
        .add(AtomConstraint::Valence(ValueAst::Lit(4)));
    ast.atom_mut(AtomIdx(0))
        .data
        .constraints
        .add(AtomConstraint::Degree(ValueAst::Lit(3)));
    ast.bond_mut(BondIdx(0))
        .data
        .constraints
        .add(BondConstraint::Aromatic);
    ast.dative_bond_mut(DativeBondIdx(0))
        .constraints
        .add(DativeBondConstraint::RingCount(ValueAst::Lit(1)));

    let original = ast.clone();

    ast.lift_constraints();
    assert!(ast[AtomIdx(0)].constraints.is_empty());
    ast.inline_constraints();

    assert_eq!(ast, original);
}

// region: simplify_values

/// Walks every value-bearing slot the simplifier touches: atom fields
/// (charge, isotope, implicit-h, lone-pairs, spin), inline atom constraints,
/// bond order/charge/spin, dative ring constraint, aromatic-system
/// charge/spin/electrons, multicenter charge/spin/electrons, molecule-scope
/// `ChargeSum::sum`, an `And` combinator with non-canonical inner shapes,
/// a `Relational` predicate, and a `SubPattern` whose pattern atom carries
/// non-canonical values too. Each non-canonical shape simplifies to its
/// canonical form, exercising the full recursion.
#[rstest]
fn test_molecule_ast_simplify_values_reduces_throughout() {
    let mut ast = MoleculeAst::from_parts(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::N),
        ],
        vec![(AtomIdx(0), AtomIdx(1), BondAst::from_order(1))],
        vec![(vec![AtomIdx(0)], AtomIdx(1), DativeBondAst::from_order(1))],
        vec![(vec![AtomIdx(0), AtomIdx(1)], AromaticSystemAst::default())],
        vec![(vec![AtomIdx(0), AtomIdx(1)], MulticenterBondAst::default())],
        vec![(
            AtomIdx(0),
            AtomIdx(1),
            NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        Constraints::default(),
    );

    // -- Non-canonical shapes seeded across the structure --------------
    // Atom 0: charge as Expr::Lit (lifts to ValueAst::Lit), isotope_mass as
    // Expr::Neg(Lit) (lifts to IsotopeAst::Lit), implicit_hydrogens as
    // Expr::Lit, lone_pairs as Expr::Neg(Neg(_)) (folds), spin both fields
    // wrapped in Expr.
    {
        let atom = ast.atom_mut(AtomIdx(0)).data;
        atom.charge = ValueAst::Expr(Expr::Lit(2));
        atom.isotope_mass = IsotopeAst::Expr(Expr::Neg(Box::new(Expr::Lit(13))));
        atom.implicit_hydrogens = ImplicitHydrogensAst::Expr(Expr::Lit(3));
        atom.lone_pairs = ValueAst::Expr(Expr::Neg(Box::new(Expr::Neg(Box::new(Expr::Lit(1))))));
        atom.spin =
            SpinStateAst::from_values(ValueAst::Expr(Expr::Lit(0)), ValueAst::Expr(Expr::Lit(1)));
        // And an inline atom constraint with a non-canonical Expr.
        atom.constraints
            .add(AtomConstraint::Valence(ValueAst::Expr(Expr::Lit(4))));
    }

    // Bond 0: order/charge/spin all wrapped, plus an inline bond ring-count
    // with a non-canonical Expr.
    {
        let bond = ast.bond_mut(BondIdx(0)).data;
        bond.order = ValueAst::Expr(Expr::Lit(1));
        bond.charge = ValueAst::Expr(Expr::Neg(Box::new(Expr::Lit(0))));
        bond.spin =
            SpinStateAst::from_values(ValueAst::Expr(Expr::Lit(0)), ValueAst::Expr(Expr::Lit(1)));
        bond.constraints
            .add(BondConstraint::RingCount(ValueAst::Expr(Expr::Lit(1))));
    }

    // Dative bond inline ring-size with non-canonical Expr.
    ast.dative_bond_mut(DativeBondIdx(0))
        .constraints
        .add(DativeBondConstraint::RingSize(ValueAst::Expr(Expr::Lit(5))));

    // Aromatic system 0: charge/electrons/spin wrapped. Three member atoms,
    // so electrons has three entries.
    {
        let ar = ast.aromatic_system_mut(AromaticSystemIdx(0));
        ar.charge = ValueAst::Expr(Expr::Lit(0));
        ar.electrons = vec![
            ValueAst::Expr(Expr::Lit(1)),
            ValueAst::Expr(Expr::Lit(1)),
            ValueAst::Expr(Expr::Lit(1)),
        ];
        ar.spin =
            SpinStateAst::from_values(ValueAst::Expr(Expr::Lit(0)), ValueAst::Expr(Expr::Lit(1)));
    }

    // Multicenter bond 0: same pattern, three member atoms.
    {
        let mc = ast.multicenter_bond_mut(MulticenterBondIdx(0));
        mc.charge = ValueAst::Expr(Expr::Lit(0));
        mc.electrons = vec![
            ValueAst::Expr(Expr::Lit(1)),
            ValueAst::Expr(Expr::Lit(1)),
            ValueAst::Expr(Expr::Lit(0)),
        ];
        mc.spin =
            SpinStateAst::from_values(ValueAst::Expr(Expr::Lit(0)), ValueAst::Expr(Expr::Lit(1)));
    }

    // Molecule-scope constraints: a ChargeSum, a Relational predicate,
    // an And combinator wrapping non-canonical leaves, and a SubPattern
    // whose pattern atom 0 carries a non-canonical charge.
    let mut pattern = MoleculeAst::from_parts(
        vec![AtomAst::from_element(Element::C)],
        vec![],
        vec![],
        vec![],
        vec![],
        vec![],
        Constraints::default(),
    );
    pattern.atom_mut(AtomIdx(0)).data.charge = ValueAst::Expr(Expr::Lit(-3));

    ast.constraints_mut()
        .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomIdx(0), AtomIdx(1)]),
            sum: ValueAst::Expr(Expr::Lit(0)),
        }));
    ast.constraints_mut().push(Constraint::Relational(
        RelationalConstraint::AromaticSystemAllAtoms {
            system: AromaticSystemIdx(0),
            predicate: Box::new(AtomConstraint::Valence(ValueAst::Expr(Expr::Lit(4)))),
        },
    ));
    ast.constraints_mut()
        .push(Constraint::And(vec![Constraint::Atom(
            AtomIdx(1),
            AtomConstraint::Degree(ValueAst::Expr(Expr::Lit(3))),
        )]));
    ast.constraints_mut()
        .push(Constraint::Molecule(MoleculeConstraint::SubPattern {
            anchor: {
                let mut a = SubPatternAnchor::new();
                a.push_atom(AtomIdx(0), AtomIdx(0));
                a
            },
            pattern: Box::new(pattern),
        }));

    ast.simplify_values();

    // -- Atom 0 ---------------------------------------------------------
    let atom = ast.atom(AtomIdx(0)).data;
    assert_eq!(atom.charge, ValueAst::Lit(2));
    assert_eq!(atom.isotope_mass, IsotopeAst::Lit(-13));
    assert_eq!(atom.implicit_hydrogens, ImplicitHydrogensAst::Lit(3));
    assert_eq!(atom.lone_pairs, ValueAst::Lit(1));
    assert_eq!(atom.spin, SpinStateAst::new(0, 1));
    assert_eq!(
        atom.constraints,
        AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Lit(4))]),
    );

    // -- Bond 0 ---------------------------------------------------------
    let bond = ast.bond(BondIdx(0)).data;
    assert_eq!(bond.order, ValueAst::Lit(1));
    // Neg(Lit(0)) is preserved by Expr::simplify but the Expr is not at the
    // ValueAst-Expr top, so ValueAst::simplify lifts via Lit(-0) = Lit(0).
    assert_eq!(bond.charge, ValueAst::Lit(0));
    assert_eq!(bond.spin, SpinStateAst::new(0, 1));
    assert_eq!(
        bond.constraints,
        BondConstraints::from_iter([BondConstraint::RingCount(ValueAst::Lit(1))]),
    );

    // -- Dative bond 0 --------------------------------------------------
    assert_eq!(
        ast[DativeBondIdx(0)].constraints,
        DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(5))]),
    );

    // -- Aromatic system 0 ---------------------------------------------
    let ar = &ast[AromaticSystemIdx(0)];
    assert_eq!(ar.charge, ValueAst::Lit(0));
    assert_eq!(
        ar.electrons,
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(1)],
    );
    assert_eq!(ar.spin, SpinStateAst::new(0, 1));

    // -- Multicenter bond 0 ---------------------------------------------
    let mc = &ast[MulticenterBondIdx(0)];
    assert_eq!(mc.charge, ValueAst::Lit(0));
    assert_eq!(
        mc.electrons,
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(0)],
    );
    assert_eq!(mc.spin, SpinStateAst::new(0, 1));

    // -- Molecule-scope constraints ------------------------------------
    let cs: Vec<&Constraint> = ast.constraints().iter().collect();
    match cs[0] {
        Constraint::Molecule(MoleculeConstraint::ChargeSum { atoms: _, sum }) => {
            assert_eq!(sum, &ValueAst::Lit(0));
        }
        c => panic!("expected ChargeSum, got {c:?}"),
    }
    match cs[1] {
        Constraint::Relational(RelationalConstraint::AromaticSystemAllAtoms {
            predicate, ..
        }) => {
            assert_eq!(**predicate, AtomConstraint::Valence(ValueAst::Lit(4)));
        }
        c => panic!("expected Relational AromaticSystemAllAtoms, got {c:?}"),
    }
    match cs[2] {
        Constraint::And(xs) => match &xs[0] {
            Constraint::Atom(idx, AtomConstraint::Degree(v)) => {
                assert_eq!(*idx, AtomIdx(1));
                assert_eq!(v, &ValueAst::Lit(3));
            }
            c => panic!("expected Atom Degree leaf, got {c:?}"),
        },
        c => panic!("expected And, got {c:?}"),
    }
    match cs[3] {
        Constraint::Molecule(MoleculeConstraint::SubPattern { pattern, .. }) => {
            // The SubPattern's pattern molecule was simplified recursively.
            assert_eq!(pattern.atom(AtomIdx(0)).data.charge, ValueAst::Lit(-3));
        }
        c => panic!("expected SubPattern, got {c:?}"),
    }
}

// endregion: simplify_values
