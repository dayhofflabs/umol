use std::collections::HashSet;

use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;
use umol_graph_core::{
    AutomorphismAlgorithm, BiconnectedComponentsAlgorithm, ConnectedComponentsAlgorithm,
    CycleEnumerationAlgorithm, EdgeId, MatchingEnumerationAlgorithm, MaxIndependentSetAlgorithm,
    MaxMatchingAlgorithm, NodeId, ShortestCycleAlgorithm, SubgraphIsomorphismAlgorithm,
};

use super::super::aromatic::AromaticSystemAst;
use super::super::atom::{AtomAst, ElementAst, IsotopeMassAst};
use super::super::bond::BondAst;
use super::super::boolean::BooleanAst;
use super::super::constraint::{
    AtomConstraint, AtomConstraints, BondConstraint, BondConstraints, Constraint, Constraints,
    DativeBondConstraint, DativeBondConstraints, MoleculeConstraint, RelationalConstraint,
    RingScope,
};
use super::super::dative::DativeBondAst;
use super::super::electrons::ElectronCountsAst;
use super::super::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::multicenter::MulticenterBondAst;
use super::super::noncovalent::{NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst};
use super::super::ring::{RingFamily, RingSet};
use super::super::spin::SpinStateAst;
use super::super::stereo::{StereoAtomAst, StereoCosetAst, StereoKind};
use super::super::value::ValueAst;
use super::MoleculeAst;
use crate::{mol, mol_ground};

fn ground_atom() -> AtomAst {
    let mut a = AtomAst::from_element(Element::C);
    a.isotope_mass = IsotopeMassAst::Natural;
    a.charge = ValueAst::Lit(0);
    a.implicit_hydrogens = ValueAst::Lit(4);
    a.lone_pairs = ValueAst::Lit(0);
    a.spin = SpinStateAst::from((0_u8, 1_u8));
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
    assert_eq!(m.atoms().count(), 0);
    assert_eq!(m.bonds().count(), 0);
    assert_eq!(m.dative_bonds().count(), 0);
    assert_eq!(m.aromatic_systems().count(), 0);
    assert_eq!(m.multicenter_bonds().count(), 0);
    assert_eq!(m.noncovalent_bonds().count(), 0);
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
    let bonds = vec![(AtomId(0), AtomId(1), BondAst::from_order(1))];
    let m = MoleculeAst::from_atoms_and_bonds(atoms, bonds);
    assert_eq!(m.atoms().count(), 2);
    assert_eq!(m.bonds().count(), 1);
    assert_eq!(m.dative_bonds().count(), 0);
    assert_eq!(m.aromatic_systems().count(), 0);
    assert_eq!(m.multicenter_bonds().count(), 0);
    assert_eq!(m.noncovalent_bonds().count(), 0);
    assert_eq!(m.atom(AtomId(0)).ast.element, ElementAst::Lit(Element::C));
    assert_eq!(m.atom(AtomId(1)).ast.element, ElementAst::Lit(Element::O));
    assert_eq!(m.bond(BondId(0)).ast.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_builder() {
    assert_eq!(MoleculeAst::builder().build(), MoleculeAst::new());
}

#[rstest]
#[case::empty(MoleculeAst::default(), true)]
#[case::ground_atom(
    mol_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
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
        Vec::new(),
        Vec::new(),
        constraints_with_molecule(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![]),
            sum: ValueAst::Undetermined,
        })),
    ),
    true,
)]
#[case::stereo_atom_ground_coset(
    MoleculeAst::from_parts(
        vec![ground_atom()],
        vec![], vec![], vec![], vec![], vec![],
        vec![(AtomId(0), vec![], StereoAtomAst::new(StereoKind::Tetrahedral, 1u32))],
        Vec::new(),
        Constraints::new(),
    ),
    true,
)]
#[case::stereo_atom_undetermined_coset(
    MoleculeAst::from_parts(
        vec![ground_atom()],
        vec![], vec![], vec![], vec![], vec![],
        vec![(AtomId(0), vec![], StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined))],
        Vec::new(),
        Constraints::new(),
    ),
    false,
)]
fn test_molecule_ast_is_ground(#[case] ast: MoleculeAst, #[case] expected: bool) {
    assert_eq!(ast.is_ground(), expected);
}

#[rstest]
#[case::hub(AtomId(0), vec![(AtomId(1), BondId(0)), (AtomId(2), BondId(1))])]
#[case::leaf_o(AtomId(1), vec![(AtomId(0), BondId(0))])]
#[case::leaf_n(AtomId(2), vec![(AtomId(0), BondId(1))])]
fn test_molecule_ast_neighbors(#[case] atom: AtomId, #[case] expected: Vec<(AtomId, BondId)>) {
    let ast = MoleculeAst::from_parts(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::O),
            AtomAst::from_element(Element::N),
        ],
        vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(0), AtomId(2), BondAst::from_order(2)),
        ],
        vec![],
        vec![],
        vec![],
        vec![],
        Vec::new(),
        Vec::new(),
        Constraints::default(),
    );
    let nbrs: Vec<(AtomId, BondId)> = ast
        .neighbors(atom)
        .map(|n| (n.atom_id(), n.bond_id()))
        .collect();
    assert_eq!(nbrs, expected);
}

#[rstest]
fn test_molecule_builder_add_aromatic_system() {
    let ast = MoleculeAst::from_parts(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        vec![],
        vec![],
        vec![],
        vec![],
        Vec::new(),
        Vec::new(),
        Constraints::default(),
    );
    let mut b = ast.edit();
    let id = b.add_aromatic_system(vec![AtomId(0), AtomId(1)], AromaticSystemAst::default());
    let new_ast = b.build();
    assert_eq!(id, AromaticSystemId(0));
    let new_atoms: Vec<AtomId> = new_ast
        .aromatic_system(AromaticSystemId(0))
        .atom_ids()
        .collect();
    assert_eq!(new_atoms, vec![AtomId(0), AtomId(1)]);
    assert_eq!(
        new_ast.aromatic_systems().ids().collect::<Vec<_>>(),
        vec![AromaticSystemId(0)]
    );
    assert_eq!(
        ast.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemId>::new()
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
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(1), AtomId(2), BondAst::from_order(2)),
            (AtomId(2), AtomId(3), BondAst::from_order(1)),
        ],
        vec![(vec![AtomId(2)], AtomId(3), DativeBondAst::from_order(1))],
        vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemAst::default(),
        )],
        vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondAst::default(),
        )],
        vec![(
            AtomId(0),
            AtomId(3),
            NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        Vec::new(),
        Vec::new(),
        Constraints::default(),
    )
}

#[rstest]
#[case::c_c(BondId(0), AtomId(0), AtomId(1), ValueAst::Lit(1))]
#[case::c_n(BondId(1), AtomId(1), AtomId(2), ValueAst::Lit(2))]
#[case::n_o(BondId(2), AtomId(2), AtomId(3), ValueAst::Lit(1))]
fn test_molecule_ast_bond(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] id: BondId,
    #[case] first: AtomId,
    #[case] second: AtomId,
    #[case] order: ValueAst,
) {
    let bv = ast.bond(id);
    assert_eq!(bv.id, id);
    assert_eq!(bv.atom_ids()[0], first);
    assert_eq!(bv.atom_ids()[1], second);
    assert_eq!(bv.ast.order, order);
}

#[rstest]
fn test_molecule_ast_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(BondId, AtomId, AtomId, ValueAst)> = ast
        .bonds()
        .iter()
        .map(|v| (v.id, v.atom_ids()[0], v.atom_ids()[1], v.ast.order.clone()))
        .collect();
    assert_eq!(
        projected,
        vec![
            (BondId(0), AtomId(0), AtomId(1), ValueAst::Lit(1)),
            (BondId(1), AtomId(1), AtomId(2), ValueAst::Lit(2)),
            (BondId(2), AtomId(2), AtomId(3), ValueAst::Lit(1)),
        ]
    );
}

#[rstest]
fn test_molecule_ast_dative_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let dv = ast.dative_bond(DativeBondId(0));
    assert_eq!(dv.id, DativeBondId(0));
    assert_eq!(dv.acceptor_id(), AtomId(3));
    assert_eq!(dv.donor_ids().collect::<Vec<_>>(), vec![AtomId(2)]);
    assert_eq!(
        dv.atom_ids().collect::<Vec<_>>(),
        vec![AtomId(2), AtomId(3)]
    );
    assert_eq!(dv.ast.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_dative_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(DativeBondId, Vec<AtomId>, AtomId)> = ast
        .dative_bonds()
        .iter()
        .map(|v| (v.id, v.donor_ids().collect(), v.acceptor_id()))
        .collect();
    assert_eq!(
        projected,
        vec![(DativeBondId(0), vec![AtomId(2)], AtomId(3))]
    );
}

#[rstest]
fn test_molecule_ast_aromatic_system(#[from(rich_molecule)] ast: MoleculeAst) {
    let av = ast.aromatic_system(AromaticSystemId(0));
    assert_eq!(av.id, AromaticSystemId(0));
    assert_eq!(
        av.atom_ids().collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2)]
    );
    assert_eq!(
        av.bond_ids().collect::<Vec<_>>(),
        vec![BondId(0), BondId(1)]
    );
}

#[rstest]
fn test_molecule_ast_aromatic_systems(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(AromaticSystemId, Vec<AtomId>, Vec<BondId>)> = ast
        .aromatic_systems()
        .iter()
        .map(|v| (v.id, v.atom_ids().collect(), v.bond_ids().collect()))
        .collect();
    assert_eq!(
        projected,
        vec![(
            AromaticSystemId(0),
            vec![AtomId(0), AtomId(1), AtomId(2)],
            vec![BondId(0), BondId(1)],
        )]
    );
}

#[rstest]
fn test_molecule_ast_multicenter_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mv = ast.multicenter_bond(MulticenterBondId(0));
    assert_eq!(mv.id, MulticenterBondId(0));
    assert_eq!(
        mv.atom_ids().collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2)]
    );
}

#[rstest]
fn test_molecule_ast_multicenter_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(MulticenterBondId, Vec<AtomId>)> = ast
        .multicenter_bonds()
        .iter()
        .map(|v| (v.id, v.atom_ids().collect()))
        .collect();
    assert_eq!(
        projected,
        vec![(MulticenterBondId(0), vec![AtomId(0), AtomId(1), AtomId(2)],)]
    );
}

#[rstest]
fn test_molecule_ast_noncovalent_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let nv = ast.noncovalent_bond(NoncovalentBondId(0));
    assert_eq!(nv.id, NoncovalentBondId(0));
    assert_eq!(nv.atom_ids(), [AtomId(0), AtomId(3)]);
}

#[rstest]
fn test_molecule_ast_noncovalent_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(NoncovalentBondId, [AtomId; 2])> = ast
        .noncovalent_bonds()
        .iter()
        .map(|v| (v.id, v.atom_ids()))
        .collect();
    assert_eq!(
        projected,
        vec![(NoncovalentBondId(0), [AtomId(0), AtomId(3)])]
    );
}

#[rstest]
#[case::forward(AtomId(0), AtomId(1), Some(BondId(0)))]
#[case::reverse(AtomId(1), AtomId(0), Some(BondId(0)))]
#[case::non_adjacent(AtomId(0), AtomId(3), None)]
fn test_bond_views_connecting_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<BondId>,
) {
    assert_eq!(ast.bonds().connecting_id(a, b), expected);
}

#[rstest]
#[case::matched(AtomId(3), vec![AtomId(2)], Some(DativeBondId(0)))]
#[case::role_swap(AtomId(2), vec![AtomId(3)], None)]
#[case::wrong_donor(AtomId(3), vec![AtomId(1)], None)]
fn test_dative_bond_views_connecting_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] acceptor: AtomId,
    #[case] donors: Vec<AtomId>,
    #[case] expected: Option<DativeBondId>,
) {
    assert_eq!(
        ast.dative_bonds().connecting_id(acceptor, &donors),
        expected
    );
}

#[rstest]
#[case::forward(AtomId(0), AtomId(3), Some(NoncovalentBondId(0)))]
#[case::reverse(AtomId(3), AtomId(0), Some(NoncovalentBondId(0)))]
#[case::unrelated(AtomId(0), AtomId(1), None)]
fn test_noncovalent_bond_views_connecting_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<NoncovalentBondId>,
) {
    assert_eq!(ast.noncovalent_bonds().connecting_id(a, b), expected);
}

#[rstest]
#[case::donor(AtomId(2), vec![DativeBondId(0)])]
#[case::acceptor(AtomId(3), vec![DativeBondId(0)])]
#[case::outside(AtomId(0), vec![])]
fn test_dative_bond_views_incident_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<DativeBondId>,
) {
    let inc: Vec<_> = ast.dative_bonds().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomId(1), vec![AromaticSystemId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_aromatic_system_views_incident_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<AromaticSystemId>,
) {
    let inc: Vec<_> = ast.aromatic_systems().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::member(AtomId(0), vec![MulticenterBondId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_multicenter_bond_views_incident_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<MulticenterBondId>,
) {
    let inc: Vec<_> = ast.multicenter_bonds().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::first(AtomId(0), vec![NoncovalentBondId(0)])]
#[case::second(AtomId(3), vec![NoncovalentBondId(0)])]
#[case::outside(AtomId(1), vec![])]
fn test_noncovalent_bond_views_incident_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    let inc: Vec<_> = ast.noncovalent_bonds().incident_ids(atom).collect();
    assert_eq!(inc, expected);
}

#[rstest]
#[case::full(vec![AtomId(2), AtomId(3)], vec![DativeBondId(0)])]
#[case::partial_only(vec![AtomId(0), AtomId(2)], vec![])]
#[case::disjoint(vec![AtomId(0), AtomId(1)], vec![])]
fn test_dative_bond_views_induced_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<DativeBondId>,
) {
    assert_eq!(ast.dative_bonds().induced_ids(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AromaticSystemId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
#[case::disjoint(vec![AtomId(3)], vec![])]
fn test_aromatic_system_views_induced_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<AromaticSystemId>,
) {
    assert_eq!(ast.aromatic_systems().induced_ids(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![MulticenterBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
#[case::disjoint(vec![AtomId(3)], vec![])]
fn test_multicenter_bond_views_induced_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<MulticenterBondId>,
) {
    assert_eq!(ast.multicenter_bonds().induced_ids(&atoms), expected);
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(3)], vec![NoncovalentBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
#[case::disjoint(vec![AtomId(1), AtomId(2)], vec![])]
fn test_noncovalent_bond_views_induced_ids(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    assert_eq!(ast.noncovalent_bonds().induced_ids(&atoms), expected);
}

#[rstest]
#[case::forward(AtomId(0), AtomId(1), Some(BondId(0)))]
#[case::reverse(AtomId(1), AtomId(0), Some(BondId(0)))]
#[case::non_adjacent(AtomId(0), AtomId(3), None)]
fn test_bond_views_connecting(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<BondId>,
) {
    assert_eq!(ast.bonds().connecting(a, b).map(|v| v.id), expected);
}

#[rstest]
#[case::pair(vec![AtomId(0), AtomId(1)], vec![BondId(0)])]
#[case::triangle(vec![AtomId(0), AtomId(1), AtomId(2)], vec![BondId(0), BondId(1)])]
#[case::singleton(vec![AtomId(0)], vec![])]
fn test_bond_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<BondId>,
) {
    let mut got: Vec<BondId> = ast
        .bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    got.sort_unstable();
    assert_eq!(got, expected);
}

#[rstest]
#[case::donor(AtomId(2), vec![DativeBondId(0)])]
#[case::acceptor(AtomId(3), vec![DativeBondId(0)])]
#[case::outside(AtomId(0), vec![])]
fn test_dative_bond_views_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<DativeBondId>,
) {
    let got: Vec<DativeBondId> = ast.dative_bonds().incident(atom).map(|v| v.id).collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::matched(AtomId(3), vec![AtomId(2)], Some(DativeBondId(0)))]
#[case::role_swap(AtomId(2), vec![AtomId(3)], None)]
fn test_dative_bond_views_connecting(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] acceptor: AtomId,
    #[case] donors: Vec<AtomId>,
    #[case] expected: Option<DativeBondId>,
) {
    assert_eq!(
        ast.dative_bonds()
            .connecting(acceptor, &donors)
            .map(|v| v.id),
        expected
    );
}

#[rstest]
#[case::full(vec![AtomId(2), AtomId(3)], vec![DativeBondId(0)])]
#[case::partial_only(vec![AtomId(0), AtomId(2)], vec![])]
fn test_dative_bond_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<DativeBondId>,
) {
    let got: Vec<DativeBondId> = ast
        .dative_bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::member(AtomId(1), vec![AromaticSystemId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_aromatic_system_views_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<AromaticSystemId>,
) {
    let got: Vec<AromaticSystemId> = ast
        .aromatic_systems()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(AromaticSystemId(0)),
)]
#[case::subset(HashSet::from([AtomId(0), AtomId(1)]), None)]
fn test_aromatic_system_views_connecting(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<AromaticSystemId>,
) {
    assert_eq!(
        ast.aromatic_systems().connecting(atoms).map(|v| v.id),
        expected
    );
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AromaticSystemId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
fn test_aromatic_system_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<AromaticSystemId>,
) {
    let got: Vec<AromaticSystemId> = ast
        .aromatic_systems()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::member(AtomId(0), vec![MulticenterBondId(0)])]
#[case::outside(AtomId(3), vec![])]
fn test_multicenter_bond_views_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<MulticenterBondId>,
) {
    let got: Vec<MulticenterBondId> = ast
        .multicenter_bonds()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(MulticenterBondId(0)),
)]
#[case::subset(HashSet::from([AtomId(0), AtomId(1)]), None)]
fn test_multicenter_bond_views_connecting(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<MulticenterBondId>,
) {
    assert_eq!(
        ast.multicenter_bonds().connecting(atoms).map(|v| v.id),
        expected,
    );
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(1), AtomId(2)], vec![MulticenterBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
fn test_multicenter_bond_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<MulticenterBondId>,
) {
    let got: Vec<MulticenterBondId> = ast
        .multicenter_bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::first(AtomId(0), vec![NoncovalentBondId(0)])]
#[case::second(AtomId(3), vec![NoncovalentBondId(0)])]
#[case::outside(AtomId(1), vec![])]
fn test_noncovalent_bond_views_incident(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    let got: Vec<NoncovalentBondId> = ast
        .noncovalent_bonds()
        .incident(atom)
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::forward(AtomId(0), AtomId(3), Some(NoncovalentBondId(0)))]
#[case::reverse(AtomId(3), AtomId(0), Some(NoncovalentBondId(0)))]
#[case::unrelated(AtomId(0), AtomId(1), None)]
fn test_noncovalent_bond_views_connecting(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] a: AtomId,
    #[case] b: AtomId,
    #[case] expected: Option<NoncovalentBondId>,
) {
    assert_eq!(
        ast.noncovalent_bonds().connecting(a, b).map(|v| v.id),
        expected,
    );
}

#[rstest]
#[case::full(vec![AtomId(0), AtomId(3)], vec![NoncovalentBondId(0)])]
#[case::partial(vec![AtomId(0), AtomId(1)], vec![])]
fn test_noncovalent_bond_views_induced(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: Vec<AtomId>,
    #[case] expected: Vec<NoncovalentBondId>,
) {
    let got: Vec<NoncovalentBondId> = ast
        .noncovalent_bonds()
        .induced(&atoms)
        .into_iter()
        .map(|v| v.id)
        .collect();
    assert_eq!(got, expected);
}

#[rstest]
#[case::atom_0(AtomId(0), Element::C)]
#[case::atom_1(AtomId(1), Element::C)]
#[case::atom_2(AtomId(2), Element::N)]
#[case::atom_3(AtomId(3), Element::O)]
fn test_molecule_ast_atom(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] id: AtomId,
    #[case] element: Element,
) {
    let av = ast.atom(id);
    assert_eq!(av.id, id);
    assert_eq!(av.ast.element, ElementAst::Lit(element));
}

#[rstest]
fn test_molecule_ast_is_empty() {
    assert!(MoleculeAst::default().is_empty());
}

#[rstest]
fn test_molecule_ast_is_empty_rich(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(!ast.is_empty());
}

#[rstest]
fn test_molecule_ast_has_constraints_empty() {
    assert!(!MoleculeAst::default().has_constraints());
}

#[rstest]
fn test_molecule_ast_has_constraints_rich(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(!ast.has_constraints());
}

#[rstest]
fn test_molecule_ast_has_dative_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast.has_dative_bonds());
}

#[rstest]
fn test_molecule_ast_has_aromatic_systems(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast.has_aromatic_systems());
}

#[rstest]
fn test_molecule_ast_has_multicenter_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast.has_multicenter_bonds());
}

#[rstest]
fn test_molecule_ast_has_noncovalent_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast.has_noncovalent_bonds());
}

#[rstest]
fn test_molecule_ast_has_overlays(#[from(rich_molecule)] ast: MoleculeAst) {
    assert!(ast.has_overlays());
}

#[rstest]
fn test_molecule_ast_has_overlays_empty() {
    assert!(!MoleculeAst::default().has_overlays());
}

#[rstest]
fn test_molecule_ast_has_stereo_atoms() {
    let ast = mol!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]] :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"#
    );
    assert!(ast.has_stereo_atoms());
    assert!(!ast.has_stereo_bonds());
    assert!(ast.has_overlays());
}

#[rstest]
fn test_molecule_ast_has_stereo_bonds() {
    let ast = mol!(
        r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]] :stereo-bonds [{:site 1 :ligands [0 3] :type "Ct1"}]}"#
    );
    assert!(ast.has_stereo_bonds());
    assert!(!ast.has_stereo_atoms());
    assert!(ast.has_overlays());
}

#[rstest]
fn test_molecule_ast_atoms(#[from(rich_molecule)] ast: MoleculeAst) {
    let projected: Vec<(AtomId, ElementAst)> = ast
        .atoms()
        .iter()
        .map(|v| (v.id, v.ast.element.clone()))
        .collect();
    assert_eq!(
        projected,
        vec![
            (AtomId(0), ElementAst::Lit(Element::C)),
            (AtomId(1), ElementAst::Lit(Element::C)),
            (AtomId(2), ElementAst::Lit(Element::N)),
            (AtomId(3), ElementAst::Lit(Element::O)),
        ]
    );
}

#[test]
fn test_bond_views_induced_ids() {
    let ast = MoleculeAst::from_parts(
        vec![
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
            AtomAst::from_element(Element::C),
        ],
        vec![
            (AtomId(0), AtomId(1), BondAst::from_order(1)),
            (AtomId(1), AtomId(2), BondAst::from_order(1)),
            (AtomId(0), AtomId(2), BondAst::from_order(1)),
        ],
        vec![],
        vec![],
        vec![],
        vec![],
        Vec::new(),
        Vec::new(),
        Constraints::default(),
    );
    let bonds = ast.bonds().induced_ids(&[AtomId(0), AtomId(1)]);
    assert_eq!(bonds, vec![BondId(0)]);

    let mut all = ast.bonds().induced_ids(&[AtomId(0), AtomId(1), AtomId(2)]);
    all.sort_unstable();
    assert_eq!(all, vec![BondId(0), BondId(1), BondId(2)]);
}

fn chain(n: usize) -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); n];
    let bonds: Vec<_> = (0..n.saturating_sub(1))
        .map(|i| {
            (
                AtomId(i as u32),
                AtomId((i + 1) as u32),
                BondAst::from_order(1),
            )
        })
        .collect();
    MoleculeAst::from_atoms_and_bonds(atoms, bonds)
}

fn ring(n: usize) -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); n];
    let bonds: Vec<_> = (0..n)
        .map(|i| {
            (
                AtomId(i as u32),
                AtomId(((i + 1) % n) as u32),
                BondAst::from_order(1),
            )
        })
        .collect();
    MoleculeAst::from_atoms_and_bonds(atoms, bonds)
}

fn two_components() -> MoleculeAst {
    let atoms = vec![AtomAst::from_element(Element::C); 4];
    let bonds = vec![
        (AtomId(0), AtomId(1), BondAst::from_order(1)),
        (AtomId(2), AtomId(3), BondAst::from_order(1)),
    ];
    MoleculeAst::from_parts(
        atoms,
        bonds,
        vec![],
        vec![],
        vec![],
        vec![],
        Vec::new(),
        Vec::new(),
        Constraints::default(),
    )
}

#[rstest]
#[case::isolated(chain(1), AtomId(0), 0)]
#[case::chain_end(chain(3), AtomId(0), 1)]
#[case::chain_mid(chain(3), AtomId(1), 2)]
#[case::ring_vertex(ring(6), AtomId(0), 2)]
fn test_molecule_ast_degree(
    #[case] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: usize,
) {
    assert_eq!(ast.graph().degree(atom), expected);
}

#[rstest]
#[case::single(chain(3), 1)]
#[case::two(two_components(), 2)]
#[case::empty(MoleculeAst::default(), 0)]
fn test_molecule_ast_connected_components(#[case] ast: MoleculeAst, #[case] expected: usize) {
    let cc = ast
        .graph()
        .connected_components(ConnectedComponentsAlgorithm::Bfs);
    assert_eq!(cc.len(), expected);
}

#[rstest]
#[case::ring_6(ring(6), 1)]
#[case::chain(chain(5), 0)]
fn test_molecule_ast_biconnected_components(#[case] ast: MoleculeAst, #[case] expected: usize) {
    let bcc = ast
        .graph()
        .biconnected_components(BiconnectedComponentsAlgorithm::Tarjan);
    assert_eq!(bcc.len(), expected);
}

#[rstest]
#[case::ring_bond(ring(6), BondId(0), Some(6))]
#[case::chain_bond(chain(3), BondId(0), None)]
fn test_molecule_ast_shortest_cycle_through_bond(
    #[case] ast: MoleculeAst,
    #[case] bond: BondId,
    #[case] expected: Option<usize>,
) {
    assert_eq!(
        ast.graph()
            .shortest_cycle_through_bond(bond, ShortestCycleAlgorithm::Bfs),
        expected
    );
}

#[rstest]
#[case::ring_atom(ring(6), AtomId(0), Some(6))]
#[case::chain_atom(chain(3), AtomId(1), None)]
fn test_molecule_ast_shortest_cycle_through_atom(
    #[case] ast: MoleculeAst,
    #[case] atom: AtomId,
    #[case] expected: Option<usize>,
) {
    assert_eq!(
        ast.graph()
            .shortest_cycle_through_atom(atom, ShortestCycleAlgorithm::Bfs),
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
    let cycles = ast
        .graph()
        .enumerate_cycles(max_size, CycleEnumerationAlgorithm::Vismara);
    assert_eq!(cycles.len(), expected);
}

#[rstest]
#[case::triangle(ring(3), 1)]
#[case::chain_3(chain(3), 2)]
fn test_molecule_ast_maximum_independent_set(#[case] ast: MoleculeAst, #[case] expected: usize) {
    let mis = ast
        .graph()
        .maximum_independent_set(MaxIndependentSetAlgorithm::BranchAndBound);
    assert_eq!(mis.len(), expected);
}

#[rstest]
#[case::chain_4(chain(4), 2)]
#[case::ring_6(ring(6), 3)]
#[case::single(chain(1), 0)]
fn test_molecule_ast_maximum_matching(#[case] ast: MoleculeAst, #[case] expected_size: usize) {
    let m = ast.graph().maximum_matching(MaxMatchingAlgorithm::Edmonds);
    assert_eq!(m.size(), expected_size);
}

#[test]
fn test_bond_matching_mate() {
    let ast = chain(4);
    let m = ast.graph().maximum_matching(MaxMatchingAlgorithm::Edmonds);
    assert!(m.is_matched(AtomId(0)));
    let mate = m.mate(AtomId(0));
    assert!(mate.is_some());
}

#[rstest]
#[case::ring_6(ring(6), 2)]
fn test_molecule_ast_enumerate_perfect_matchings(
    #[case] ast: MoleculeAst,
    #[case] expected: usize,
) {
    let ms = ast
        .graph()
        .enumerate_perfect_matchings(MatchingEnumerationAlgorithm::BranchAndBound);
    assert_eq!(ms.len(), expected);
    for m in &ms {
        assert!(m.is_perfect(ast.atoms().count()));
    }
}

#[rstest]
#[case::ring_6(ring(6), 1)]
#[case::chain_3(chain(3), 2)]
fn test_molecule_ast_automorphisms(#[case] ast: MoleculeAst, #[case] expected_orbits: usize) {
    let auto = ast
        .graph()
        .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty);
    assert_eq!(auto.orbit_count(), expected_orbits);
    assert_eq!(auto.atom_count(), ast.atoms().count());
}

#[test]
fn test_atom_automorphism_same_orbit() {
    let ast = ring(6);
    let auto = ast
        .graph()
        .automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty);
    assert!(auto.same_orbit(AtomId(0), AtomId(3)));
}

#[rstest]
fn test_molecule_ast_subgraph_isomorphisms() {
    let target = ring(6);
    let query = chain(2);
    let mut matches = target.graph().subgraph_isomorphisms(
        &query.graph(),
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
    );
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![
            vec![AtomId(0), AtomId(1)],
            vec![AtomId(0), AtomId(5)],
            vec![AtomId(1), AtomId(0)],
            vec![AtomId(1), AtomId(2)],
            vec![AtomId(2), AtomId(1)],
            vec![AtomId(2), AtomId(3)],
            vec![AtomId(3), AtomId(2)],
            vec![AtomId(3), AtomId(4)],
            vec![AtomId(4), AtomId(3)],
            vec![AtomId(4), AtomId(5)],
            vec![AtomId(5), AtomId(0)],
            vec![AtomId(5), AtomId(4)],
        ]
    );
}

#[rstest]
fn test_molecule_ast_subgraph_isomorphisms_at() {
    let target = ring(6);
    let query = chain(2);
    let mut matches = target.graph().subgraph_isomorphisms_at(
        &query.graph(),
        (AtomId(0), AtomId(0)),
        &mut |_, _| true,
        &mut |_, _| true,
        SubgraphIsomorphismAlgorithm::Vf2,
    );
    matches.sort_unstable();
    assert_eq!(
        matches,
        vec![vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(5)],]
    );
}

#[rstest]
fn test_molecule_ast_induced_subgraph(#[from(rich_molecule)] ast: MoleculeAst) {
    let sub = ast.induced_subgraph(&[AtomId(0), AtomId(1), AtomId(2)]);
    let extracted = ast.extract(&sub);
    let atom_elements: Vec<_> = extracted
        .atoms()
        .iter()
        .map(|v| v.ast.element.clone())
        .collect();
    assert_eq!(
        atom_elements,
        vec![
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::N),
        ]
    );
    let bonds: Vec<(AtomId, AtomId, ValueAst)> = extracted
        .bonds()
        .iter()
        .map(|v| (v.atom_ids()[0], v.atom_ids()[1], v.ast.order.clone()))
        .collect();
    assert_eq!(
        bonds,
        vec![
            (AtomId(0), AtomId(1), ValueAst::Lit(1)),
            (AtomId(1), AtomId(2), ValueAst::Lit(2)),
        ]
    );
    assert_eq!(
        sub.atoms().mates(),
        &[
            (NodeId(0), NodeId(0)),
            (NodeId(1), NodeId(1)),
            (NodeId(2), NodeId(2))
        ]
    );
    assert_eq!(
        sub.bonds().mates(),
        &[(BondId(0), BondId(0)), (BondId(1), BondId(1))]
    );
    assert_eq!(
        sub.aromatic_systems().mates(),
        &[(AromaticSystemId(0), AromaticSystemId(0))]
    );
    assert_eq!(
        sub.multicenter_bonds().mates(),
        &[(MulticenterBondId(0), MulticenterBondId(0))]
    );
    assert_eq!(
        sub.dative_bonds().mates(),
        &[] as &[(DativeBondId, DativeBondId)]
    );
    assert_eq!(
        sub.noncovalent_bonds().mates(),
        &[] as &[(NoncovalentBondId, NoncovalentBondId)]
    );
}

#[rstest]
fn test_molecule_ast_induced_subgraph_preserves_dative(#[from(rich_molecule)] ast: MoleculeAst) {
    let sub = ast.induced_subgraph(&[AtomId(2), AtomId(3)]);
    assert_eq!(
        sub.atoms().mates(),
        &[(NodeId(0), NodeId(2)), (NodeId(1), NodeId(3))]
    );
    assert_eq!(
        sub.dative_bonds().mates(),
        &[(DativeBondId(0), DativeBondId(0))]
    );
    let extracted = ast.extract(&sub);
    let dv = extracted.dative_bond(DativeBondId(0));
    assert_eq!(dv.acceptor_id(), AtomId(1));
    assert_eq!(dv.donor_ids().collect::<Vec<_>>(), vec![AtomId(0)]);
    assert_eq!(dv.ast.order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_edits(#[from(rich_molecule)] ast: MoleculeAst) {
    use super::super::edit::{AtomRef, BondRef, Edit};
    let sub = ast.induced_subgraph(&[AtomId(0), AtomId(1), AtomId(2)]);
    assert_eq!(
        ast.edits(&sub),
        vec![Edit::RemoveTopology {
            atoms: vec![AtomRef::Id(AtomId(3))],
            bonds: vec![BondRef::Id(BondId(2))],
        }]
    );
}

#[rstest]
fn test_molecule_ast_edits_identity(#[from(rich_molecule)] ast: MoleculeAst) {
    let atom_ids: Vec<AtomId> = ast.atoms().iter().map(|v| v.id).collect();
    let sub = ast.induced_subgraph(&atom_ids);
    assert_eq!(ast.edits(&sub), Vec::new());
}

#[rstest]
fn test_molecule_ast_extract(#[from(rich_molecule)] ast: MoleculeAst) {
    let sub = ast.induced_subgraph(&[AtomId(0), AtomId(1)]);
    let extracted = ast.extract(&sub);
    assert_eq!(extracted.atoms().count(), 2);
}

#[rstest]
fn test_molecule_builder_remove_aromatic_systems(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_aromatic_systems(&[AromaticSystemId(0)]);
    let result = b.build();
    assert_eq!(
        result.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemId>::new()
    );
    assert_eq!(
        result.atoms().iter().map(|v| v.id).collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]
    );
    assert_eq!(
        result.bonds().iter().map(|v| v.id).collect::<Vec<_>>(),
        vec![BondId(0), BondId(1), BondId(2)]
    );
}

#[rstest]
fn test_molecule_builder_remove_dative_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_dative_bonds(&[DativeBondId(0)]);
    let result = b.build();
    assert_eq!(
        result.dative_bonds().ids().collect::<Vec<_>>(),
        Vec::<DativeBondId>::new()
    );
}

#[rstest]
fn test_molecule_builder_remove_multicenter_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_multicenter_bonds(&[MulticenterBondId(0)]);
    let result = b.build();
    assert_eq!(
        result.multicenter_bonds().ids().collect::<Vec<_>>(),
        Vec::<MulticenterBondId>::new()
    );
}

#[rstest]
fn test_molecule_builder_remove_noncovalent_bonds(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.remove_noncovalent_bonds(&[NoncovalentBondId(0)]);
    let result = b.build();
    assert_eq!(
        result.noncovalent_bonds().ids().collect::<Vec<_>>(),
        Vec::<NoncovalentBondId>::new()
    );
}

#[rstest]
fn test_molecule_builder_atom_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.atom_mut(AtomId(0)).ast.element = ElementAst::Lit(Element::N);
    let result = b.build();
    assert_eq!(result[AtomId(0)].element, ElementAst::Lit(Element::N));
    assert_eq!(ast[AtomId(0)].element, ElementAst::Lit(Element::C));
}

#[rstest]
fn test_molecule_builder_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.bond_mut(BondId(0)).ast.order = ValueAst::Lit(3);
    let result = b.build();
    assert_eq!(result[BondId(0)].order, ValueAst::Lit(3));
    assert_eq!(ast[BondId(0)].order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_builder_atom_constraint_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.atom_mut(AtomId(0))
        .ast
        .constraints
        .add(AtomConstraint::Degree(ValueAst::Lit(2)));
    let result = b.build();
    assert_eq!(
        result[AtomId(0)].constraints,
        AtomConstraints::from_iter([AtomConstraint::Degree(ValueAst::Lit(2))])
    );
    assert!(ast[AtomId(0)].constraints.is_empty());
}

#[rstest]
fn test_molecule_builder_add_dative_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let id = b.add_dative_bond(vec![AtomId(1)], AtomId(0), DativeBondAst::from_order(1));
    let result = b.build();
    assert_eq!(id, DativeBondId(1));
    let view = result.dative_bond(id);
    assert_eq!(view.acceptor_id(), AtomId(0));
    assert_eq!(view.donor_ids().collect::<Vec<_>>(), vec![AtomId(1)]);
}

#[rstest]
fn test_molecule_builder_add_multicenter_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let id = b.add_multicenter_bond(
        vec![AtomId(1), AtomId(2), AtomId(3)],
        MulticenterBondAst::default(),
    );
    let result = b.build();
    assert_eq!(id, MulticenterBondId(1));
    let atoms: Vec<AtomId> = result.multicenter_bond(id).atom_ids().collect();
    assert_eq!(atoms, vec![AtomId(1), AtomId(2), AtomId(3)]);
}

#[rstest]
fn test_molecule_builder_add_noncovalent_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let id = b.add_noncovalent_bond(
        [AtomId(1), AtomId(2)],
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
    );
    let result = b.build();
    assert_eq!(id, NoncovalentBondId(1));
    let view = result.noncovalent_bond(id);
    assert_eq!(view.atom_ids(), [AtomId(1), AtomId(2)]);
}

#[rstest]
fn test_molecule_builder_push_constraint_and_constraints_mut(
    #[from(rich_molecule)] ast: MoleculeAst,
) {
    let mut b = ast.edit();
    b.push_constraint(Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
    }));
    b.constraints_mut()
        .push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0)]),
            sum: ValueAst::Lit(0),
        }));
    let result = b.build();
    assert_eq!(result.constraints().len(), 2);
}

#[rstest]
fn test_molecule_builder_dative_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.dative_bond_mut(DativeBondId(0))
        .ast
        .constraints
        .add(DativeBondConstraint::ring_membership(RingScope::Size(5), 1));
    let result = b.build();
    assert!(!result[DativeBondId(0)].constraints.is_empty());
    assert!(ast[DativeBondId(0)].constraints.is_empty());
}

#[rstest]
fn test_molecule_builder_aromatic_system_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.aromatic_system_mut(AromaticSystemId(0)).ast.charge = ValueAst::Lit(0);
    let result = b.build();
    assert_eq!(result[AromaticSystemId(0)].charge, ValueAst::Lit(0));
}

#[rstest]
fn test_molecule_builder_multicenter_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.multicenter_bond_mut(MulticenterBondId(0)).ast.electrons =
        ElectronCountsAst::Lit(vec![1, 1, 0]);
    let result = b.build();
    assert_eq!(
        result[MulticenterBondId(0)].electrons,
        ElectronCountsAst::Lit(vec![1, 1, 0]),
    );
}

#[rstest]
fn test_molecule_builder_noncovalent_bond_mut(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    b.noncovalent_bond_mut(NoncovalentBondId(0)).ast.kind =
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic);
    let result = b.build();
    assert_eq!(
        result[NoncovalentBondId(0)].kind,
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
    #[case] ast: MoleculeAst,
    #[case] max_ring_size: usize,
    #[case] expected: usize,
) {
    let rs = ast.rings_with(RingFamily::Simple, max_ring_size, |_| true);
    assert_eq!(rs.count(), expected);
}

#[test]
fn test_molecule_ast_rings_with_atom_filter() {
    let ast = ring(6);
    let rs = ast.rings_with(RingFamily::Simple, 10, |a| a.0 < 3);
    assert_eq!(rs.count(), 0);
}

#[test]
fn test_molecule_ast_rings_returns_same_slot() {
    let ast = ring(6);
    let first: *const RingSet = ast.rings();
    let second: *const RingSet = ast.rings();
    assert_eq!(first, second);
}

#[test]
fn test_molecule_ast_rings_cache_survives_attribute_mutation() {
    let mut ast = ring(6);
    let first: *const RingSet = ast.rings();
    ast.atom_mut(AtomId(0)).ast.charge = ValueAst::Lit(1);
    let second: *const RingSet = ast.rings();
    assert_eq!(first, second);
}

#[test]
fn test_molecule_ast_rings_cache_reset_after_build() {
    let ast = ring(6);
    let count_before = ast.rings().count();
    let mut b = ast.edit();
    b.atom_mut(AtomId(0)).ast.element = ElementAst::Lit(Element::N);
    let next = b.build();
    assert_eq!(next.rings().count(), count_before);
}

#[test]
fn test_molecule_ast_rings_induced() {
    let ast = mol!(
        r#"{
        :atoms ["C" "C" "C" "C"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 2 "1"] [1 3 "1"] [2 3 "1"]]
    }"#
    );
    let simple_count = ast.rings_with(RingFamily::Simple, 4, |_| true).count();
    let induced_count = ast.rings_with(RingFamily::Relevant, 4, |_| true).count();
    assert_eq!(simple_count, 4);
    assert_eq!(induced_count, 4);
}

#[test]
fn test_molecule_ast_rings_induced_naphthalene() {
    let ast = mol!(
        r#"{
        :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
        :bonds [
            [0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
            [3 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]
        ]
    }"#
    );
    let simple_count = ast.rings_with(RingFamily::Simple, 10, |_| true).count();
    assert_eq!(simple_count, 2);
    let induced_count = ast.rings_with(RingFamily::Relevant, 10, |_| true).count();
    assert_eq!(induced_count, 2);
}

#[test]
fn test_rings_membership() {
    let ast = ring(6);
    let rs = ast.rings_with(RingFamily::Simple, 6, |_| true);
    assert!(rs.contains_atom(AtomId(0)));
    assert!(rs.contains_bond(BondId(0)));
    assert_eq!(rs.atom_smallest_ring_size(AtomId(0)), Some(6));
}

#[rstest]
fn test_molecule_builder_add_and_remove(#[from(rich_molecule)] ast: MoleculeAst) {
    let mut b = ast.edit();
    let new_a = b.add_atom(AtomAst::from_element(Element::Br));
    b.add_bond(AtomId(0), new_a, BondAst::from_order(1));
    b.remove_aromatic_systems(&[AromaticSystemId(0)]);
    let _compaction = b.remove(&[AtomId(3)], &[BondId(2)]);
    let result = b.build();
    let atoms: Vec<Element> = result
        .atoms()
        .iter()
        .map(|v| match v.ast.element {
            ElementAst::Lit(e) => e,
            _ => panic!("non-ground element in builder result"),
        })
        .collect();
    assert_eq!(atoms, vec![Element::C, Element::C, Element::N, Element::Br]);
    let bonds: Vec<(AtomId, AtomId, ValueAst)> = result
        .bonds()
        .iter()
        .map(|v| (v.atom_ids()[0], v.atom_ids()[1], v.ast.order.clone()))
        .collect();
    assert_eq!(
        bonds,
        vec![
            (AtomId(0), AtomId(1), ValueAst::Lit(1)),
            (AtomId(1), AtomId(2), ValueAst::Lit(2)),
            (AtomId(0), AtomId(3), ValueAst::Lit(1)),
        ]
    );
    assert_eq!(
        result.aromatic_systems().ids().collect::<Vec<_>>(),
        Vec::<AromaticSystemId>::new()
    );
    assert_eq!(
        result.dative_bonds().ids().collect::<Vec<_>>(),
        Vec::<DativeBondId>::new()
    );
    assert_eq!(
        result.noncovalent_bonds().ids().collect::<Vec<_>>(),
        Vec::<NoncovalentBondId>::new()
    );
}

#[rstest]
#[case::donor_below_acceptor(AtomId(0), AtomId(1))]
#[case::donor_above_acceptor(AtomId(1), AtomId(0))]
fn test_molecule_ast_dative_acceptor_donor(#[case] donor: AtomId, #[case] acceptor: AtomId) {
    let atoms = vec![ground_atom(), ground_atom()];
    let ast = MoleculeAst::from_parts(
        atoms,
        Vec::new(),
        vec![(vec![donor], acceptor, DativeBondAst::from_order(1))],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    let view = ast.dative_bond(DativeBondId(0));
    assert_eq!(view.acceptor_id(), acceptor);
    assert_eq!(view.donor_ids().collect::<Vec<_>>(), vec![donor]);
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
        vec![(AtomId(0), AtomId(1), bond.clone())],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    let reverse = MoleculeAst::from_parts(
        atoms_b,
        vec![(AtomId(1), AtomId(0), bond)],
        Vec::new(),
        Vec::new(),
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
        vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Constraints::new(),
    );
    let reverse = MoleculeAst::from_parts(
        atoms_b,
        Vec::new(),
        vec![(vec![AtomId(1)], AtomId(0), DativeBondAst::from_order(1))],
        Vec::new(),
        Vec::new(),
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
fn test_molecule_ast_raw_graph(#[from(rich_molecule)] ast: MoleculeAst) {
    let g = ast.raw_graph();
    assert_eq!(g.node_count(), 4);
    assert_eq!(g.edge_count(), 3);
    assert_eq!(g.edge_endpoints(EdgeId(0)), [NodeId(0), NodeId(1)]);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(AromaticSystemId(0)),
)]
#[case::subset(
    HashSet::from([AtomId(0), AtomId(1)]),
    None,
)]
#[case::disjoint(
    HashSet::from([AtomId(3)]),
    None,
)]
fn test_aromatic_system_views_connecting_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<AromaticSystemId>,
) {
    assert_eq!(ast.aromatic_systems().connecting_id(atoms), expected);
}

#[rstest]
#[case::full_match(
    HashSet::from([AtomId(0), AtomId(1), AtomId(2)]),
    Some(MulticenterBondId(0)),
)]
#[case::subset(
    HashSet::from([AtomId(0), AtomId(1)]),
    None,
)]
fn test_multicenter_bond_views_connecting_id(
    #[from(rich_molecule)] ast: MoleculeAst,
    #[case] atoms: HashSet<AtomId>,
    #[case] expected: Option<MulticenterBondId>,
) {
    assert_eq!(ast.multicenter_bonds().connecting_id(atoms), expected);
}

#[rstest]
fn test_molecule_ast_enumerate_maximum_matchings() {
    let ast = ring(4);
    let mut ms: Vec<Vec<(AtomId, AtomId)>> = ast
        .graph()
        .enumerate_maximum_matchings(MatchingEnumerationAlgorithm::BranchAndBound)
        .into_iter()
        .map(|m| {
            let mut pairs: Vec<_> = (0..ast.atoms().count())
                .map(AtomId::from)
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
            vec![(AtomId(0), AtomId(1)), (AtomId(2), AtomId(3))],
            vec![(AtomId(0), AtomId(3)), (AtomId(1), AtomId(2))],
        ]
    );
}

#[rstest]
fn test_molecule_ast_index_atom(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(ast[AtomId(2)].element, ElementAst::Lit(Element::N));
}

#[rstest]
fn test_molecule_ast_index_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(ast[BondId(1)].order, ValueAst::Lit(2));
}

#[rstest]
fn test_molecule_ast_index_dative_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(ast[DativeBondId(0)].order, ValueAst::Lit(1));
}

#[rstest]
fn test_molecule_ast_index_aromatic_system(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(
        ast[AromaticSystemId(0)].electrons,
        ElectronCountsAst::Undetermined
    );
}

#[rstest]
fn test_molecule_ast_index_multicenter_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(
        ast[MulticenterBondId(0)].electrons,
        ElectronCountsAst::Undetermined
    );
}

#[rstest]
fn test_molecule_ast_index_noncovalent_bond(#[from(rich_molecule)] ast: MoleculeAst) {
    assert_eq!(
        ast[NoncovalentBondId(0)].kind,
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)
    );
}

#[rstest]
fn test_molecule_ast_atoms_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    for a in ast.atoms_mut() {
        a.charge = ValueAst::Lit(1);
    }
    let charges: Vec<ValueAst> = ast.atoms().iter().map(|v| v.ast.charge.clone()).collect();
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
    let orders: Vec<ValueAst> = ast.bonds().iter().map(|v| v.ast.order.clone()).collect();
    assert_eq!(
        orders,
        vec![ValueAst::Lit(1), ValueAst::Lit(1), ValueAst::Lit(1)]
    );
}

#[rstest]
fn test_molecule_ast_dative_bond_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.dative_bond_mut(DativeBondId(0))
        .constraints
        .add(DativeBondConstraint::ring_membership(RingScope::Size(6), 1));
    assert_eq!(
        ast[DativeBondId(0)].constraints,
        DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(
            RingScope::Size(6),
            1
        )])
    );
}

#[rstest]
fn test_molecule_ast_aromatic_system_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.aromatic_system_mut(AromaticSystemId(0)).electrons = ElectronCountsAst::Lit(vec![1; 3]);
    assert_eq!(
        ast[AromaticSystemId(0)].electrons,
        ElectronCountsAst::Lit(vec![1, 1, 1]),
    );
}

#[rstest]
fn test_molecule_ast_aromatic_systems_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    for a in ast.aromatic_systems_mut() {
        a.electrons = ElectronCountsAst::Lit(vec![1; 3]);
    }
    let electrons: Vec<ElectronCountsAst> = ast
        .aromatic_systems()
        .iter()
        .map(|v| v.ast.electrons.clone())
        .collect();
    assert_eq!(electrons, vec![ElectronCountsAst::Lit(vec![1; 3])]);
}

#[rstest]
fn test_molecule_ast_multicenter_bond_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.multicenter_bond_mut(MulticenterBondId(0)).electrons =
        ElectronCountsAst::Lit(vec![1, 1, 0]);
    assert_eq!(
        ast[MulticenterBondId(0)].electrons,
        ElectronCountsAst::Lit(vec![1, 1, 0]),
    );
}

#[rstest]
fn test_molecule_ast_multicenter_bonds_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    for m in ast.multicenter_bonds_mut() {
        m.electrons = ElectronCountsAst::Lit(vec![1, 1, 0]);
    }
    let electrons: Vec<ElectronCountsAst> = ast
        .multicenter_bonds()
        .iter()
        .map(|v| v.ast.electrons.clone())
        .collect();
    assert_eq!(electrons, vec![ElectronCountsAst::Lit(vec![1, 1, 0])],);
}

#[rstest]
fn test_molecule_ast_noncovalent_bond_mut(#[from(rich_molecule)] mut ast: MoleculeAst) {
    ast.noncovalent_bond_mut(NoncovalentBondId(0)).kind =
        NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic);
    assert_eq!(
        ast[NoncovalentBondId(0)].kind,
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
    ast.atom_mut(AtomId(0))
        .ast
        .constraints
        .add(AtomConstraint::Valence(ValueAst::Lit(4)));
    ast.atom_mut(AtomId(2))
        .ast
        .constraints
        .add(AtomConstraint::Degree(ValueAst::Lit(3)));
    ast.bond_mut(BondId(0))
        .ast
        .constraints
        .add(BondConstraint::Aromatic(BooleanAst::Lit(true)));
    ast.dative_bond_mut(DativeBondId(0))
        .constraints
        .add(DativeBondConstraint::ring_membership(
            RingScope::All,
            ValueAst::Lit(1),
        ));

    ast.lift_constraints();

    assert!(ast[AtomId(0)].constraints.is_empty());
    assert!(ast[AtomId(2)].constraints.is_empty());
    assert!(ast[BondId(0)].constraints.is_empty());
    assert!(ast[DativeBondId(0)].constraints.is_empty());

    let mut expected = Constraints::new();
    expected.push(Constraint::Atom(
        AtomId(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));
    expected.push(Constraint::Atom(
        AtomId(2),
        AtomConstraint::Degree(ValueAst::Lit(3)),
    ));
    expected.push(Constraint::Bond(
        BondId(0),
        BondConstraint::Aromatic(BooleanAst::Lit(true)),
    ));
    expected.push(Constraint::DativeBond(
        DativeBondId(0),
        DativeBondConstraint::ring_membership(RingScope::All, ValueAst::Lit(1)),
    ));
    assert_same_constraints(ast.constraints(), &expected);
}

#[rstest]
fn test_molecule_ast_lift_constraints_appends_to_existing(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    let prior = Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(0),
    });
    ast.constraints_mut().push(prior.clone());
    ast.atom_mut(AtomId(0))
        .ast
        .constraints
        .add(AtomConstraint::Valence(ValueAst::Lit(4)));

    ast.lift_constraints();

    let mut expected = Constraints::new();
    expected.push(prior);
    expected.push(Constraint::Atom(
        AtomId(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));
    assert_same_constraints(ast.constraints(), &expected);
}

#[rstest]
fn test_molecule_ast_inline_constraints_drains_top_level_leaves(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));
    ast.constraints_mut().push(Constraint::Bond(
        BondId(0),
        BondConstraint::Aromatic(BooleanAst::Lit(true)),
    ));
    ast.constraints_mut().push(Constraint::DativeBond(
        DativeBondId(0),
        DativeBondConstraint::ring_membership(RingScope::Size(5), 1),
    ));

    ast.inline_constraints();

    assert!(ast.constraints().is_empty());
    assert_eq!(
        ast[AtomId(0)].constraints,
        AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Lit(4))])
    );
    assert_eq!(
        ast[BondId(0)].constraints,
        BondConstraints::from_iter([BondConstraint::Aromatic(BooleanAst::Lit(true))])
    );
    assert_eq!(
        ast[DativeBondId(0)].constraints,
        DativeBondConstraints::from_iter([DativeBondConstraint::ring_membership(
            RingScope::Size(5),
            1
        )])
    );
}

#[rstest]
fn test_molecule_ast_inline_constraints_last_wins_on_collision(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraint::Valence(ValueAst::Lit(3)),
    ));
    ast.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));

    ast.inline_constraints();

    // Only one Valence survives; with two competing inserts of the same kind,
    // exactly one wins (which one is unspecified). Verify count and kind.
    assert_eq!(ast[AtomId(0)].constraints.len(), 1);
    let v = ast[AtomId(0)].constraints.iter().next().unwrap().clone();
    assert!(matches!(v, AtomConstraint::Valence(_)));
}

#[rstest]
fn test_molecule_ast_inline_constraints_skips_combinator_nested(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    let leaf = Constraint::Atom(AtomId(0), AtomConstraint::Valence(ValueAst::Lit(4)));
    let nested = Constraint::And(vec![
        leaf.clone(),
        Constraint::Bond(BondId(0), BondConstraint::Aromatic(BooleanAst::Lit(true))),
    ]);
    ast.constraints_mut().push(nested.clone());

    ast.inline_constraints();

    let mut expected = Constraints::new();
    expected.push(nested);
    assert_same_constraints(ast.constraints(), &expected);
    assert!(ast[AtomId(0)].constraints.is_empty());
    assert!(ast[BondId(0)].constraints.is_empty());
}

#[rstest]
fn test_molecule_ast_inline_constraints_skips_relational_and_molecule(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    let rel = Constraint::Relational(RelationalConstraint::AromaticSystemContains {
        system: AromaticSystemId(0),
        atom: AtomId(0),
    });
    let mol = Constraint::Molecule(MoleculeConstraint::Connected {
        atoms: Some(vec![AtomId(0), AtomId(1)]),
    });
    ast.constraints_mut().push(rel.clone());
    ast.constraints_mut().push(mol.clone());
    ast.constraints_mut().push(Constraint::Atom(
        AtomId(0),
        AtomConstraint::Valence(ValueAst::Lit(4)),
    ));

    ast.inline_constraints();

    let mut expected = Constraints::new();
    expected.push(rel);
    expected.push(mol);
    assert_same_constraints(ast.constraints(), &expected);
    assert_eq!(
        ast[AtomId(0)].constraints,
        AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Lit(4))])
    );
}

#[rstest]
fn test_molecule_ast_lift_then_inline_roundtrips_inline_state(
    #[from(rich_molecule)] mut ast: MoleculeAst,
) {
    ast.atom_mut(AtomId(0))
        .ast
        .constraints
        .add(AtomConstraint::Valence(ValueAst::Lit(4)));
    ast.atom_mut(AtomId(0))
        .ast
        .constraints
        .add(AtomConstraint::Degree(ValueAst::Lit(3)));
    ast.bond_mut(BondId(0))
        .ast
        .constraints
        .add(BondConstraint::Aromatic(BooleanAst::Lit(true)));
    ast.dative_bond_mut(DativeBondId(0))
        .constraints
        .add(DativeBondConstraint::ring_membership(
            RingScope::All,
            ValueAst::Lit(1),
        ));

    let original = ast.clone();

    ast.lift_constraints();
    assert!(ast[AtomId(0)].constraints.is_empty());
    ast.inline_constraints();

    assert_eq!(ast, original);
}
