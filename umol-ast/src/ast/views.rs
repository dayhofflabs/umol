//! Read-only views over `MoleculeAst` topology and relations.
//!
//! View records bundle an index with the underlying data so consumers
//! never assemble (id, data, participants) tuples by hand. Namespace
//! types group per-relation accessors (`count`, `ids`, `iter`, `get`,
//! and `Index`) without burying them on `MoleculeAst` itself.
//!
//! Production code split per entity into submodules; tests aggregated
//! here against a shared `molecule()` fixture.

mod aromatic_system;
mod atom;
mod bond;
mod dative_bond;
mod graph;
mod multicenter_bond;
mod neighbor;
mod noncovalent_bond;

pub use aromatic_system::{
    AromaticSystemBuilderView, AromaticSystemBuilderViewMut, AromaticSystemView,
    AromaticSystemViews,
};
pub use atom::{AtomBuilderView, AtomBuilderViewMut, AtomView, AtomViewMut, AtomViews};
pub use bond::{BondBuilderView, BondBuilderViewMut, BondView, BondViewMut, BondViews};
pub use dative_bond::{
    DativeBondBuilderView, DativeBondBuilderViewMut, DativeBondView, DativeBondViews,
};
pub use graph::GraphView;
pub use multicenter_bond::{
    MulticenterBondBuilderView, MulticenterBondBuilderViewMut, MulticenterBondView,
    MulticenterBondViews,
};
pub use neighbor::NeighborView;
pub use noncovalent_bond::{
    NoncovalentBondBuilderView, NoncovalentBondBuilderViewMut, NoncovalentBondView,
    NoncovalentBondViews,
};

#[cfg(test)]
use super::atom::AtomAst;
#[cfg(test)]
use super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::aromatic::AromaticSystemAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::{
        AromaticValenceAst, AtomConstraint, Constraints, MulticenterValenceAst,
    };
    use crate::ast::dative::DativeBondAst;
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::multicenter::MulticenterBondAst;
    use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind};
    use crate::ast::rings::RingFamily;
    use crate::ast::value::ValueAst;
    use crate::mol;

    /// 4-atom molecule with one of every relation kind:
    /// atoms C C N O; bonds 0-1 single, 1-2 double, 2-3 single;
    /// dative donor=2 → acceptor=3; aromatic system [0,1,2];
    /// multicenter bond [0,1,2]; noncovalent H-bond 0-3.
    #[fixture]
    fn molecule() -> MoleculeAst {
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
            Constraints::default(),
        )
    }

    // --- AtomViews ---

    #[rstest]
    fn test_atom_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.atoms().count(), 4);
    }

    #[rstest]
    fn test_atom_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.atoms().ids().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
        );
    }

    #[rstest]
    fn test_atom_views_iter(molecule: MoleculeAst) {
        let views = molecule.atoms();
        let collected: Vec<(AtomId, AtomAst)> =
            views.iter().map(|v| (v.id, v.ast.clone())).collect();
        assert_eq!(
            collected,
            vec![
                (AtomId(0), AtomAst::from_element(Element::C)),
                (AtomId(1), AtomAst::from_element(Element::C)),
                (AtomId(2), AtomAst::from_element(Element::N)),
                (AtomId(3), AtomAst::from_element(Element::O)),
            ],
        );
    }

    #[rstest]
    fn test_atom_views_get(molecule: MoleculeAst) {
        let view = molecule.atoms().get(AtomId(2));
        assert_eq!(view.id, AtomId(2));
        assert_eq!(*view.ast, AtomAst::from_element(Element::N));
    }

    #[rstest]
    fn test_atom_views_index(molecule: MoleculeAst) {
        let atom: &AtomAst = &molecule.atoms()[AtomId(2)];
        assert_eq!(*atom, AtomAst::from_element(Element::N));
    }

    // --- AtomView ---

    #[rstest]
    fn test_atom_view_neighbors(molecule: MoleculeAst) {
        let view = molecule.atom(AtomId(1));
        let collected: Vec<(BondId, AtomId, BondAst)> = view
            .neighbors()
            .map(|n| (n.bond, n.atom, n.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(0), AtomId(0), BondAst::from_order(1)),
                (BondId(1), AtomId(2), BondAst::from_order(2)),
            ],
        );
    }

    #[rstest]
    #[case::no_incident(
        mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(3),
        ValueAst::Lit(0),
    )]
    #[case::single(
        mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(0),
        ValueAst::Lit(1),
    )]
    #[case::three_around_center(
        mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(1),
        ValueAst::Lit(3),
    )]
    #[case::double(
        mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "2"]]}"#),
        AtomId(2),
        ValueAst::Lit(2),
    )]
    #[case::undetermined_bond(
        mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "*"]]}"#),
        AtomId(0),
        ValueAst::Undetermined,
    )]
    fn test_atom_view_valence(
        #[case] molecule: MoleculeAst,
        #[case] center: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(molecule.atom(center).valence(), expected);
    }

    #[rstest]
    #[case::with_constraint(Some(AtomConstraint::valence(4)), ValueAst::Lit(4))]
    #[case::absent(None, ValueAst::Undetermined)]
    fn test_atom_view_valence_constraint(
        #[case] constraint: Option<AtomConstraint>,
        #[case] expected: ValueAst,
    ) {
        let mut atom = AtomAst::from_element(Element::C);
        if let Some(c) = constraint {
            atom.constraints.add(c);
        }
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(molecule.atom(AtomId(0)).constraints().valence(), expected);
    }

    #[rstest]
    #[case::donor(AtomId(0), ValueAst::Lit(1))]
    #[case::acceptor(AtomId(1), ValueAst::Lit(0))]
    fn test_atom_view_donated_pairs(#[case] atom: AtomId, #[case] expected: ValueAst) {
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(molecule.atom(atom).donated_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_donated_pairs_constraint() {
        let mut atom = AtomAst::from_element(Element::N);
        atom.constraints.add(AtomConstraint::donated_pairs(1));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().donated_pairs(),
            ValueAst::Lit(1),
        );
    }

    #[rstest]
    #[case::donor(AtomId(0), ValueAst::Lit(0))]
    #[case::acceptor(AtomId(1), ValueAst::Lit(1))]
    fn test_atom_view_accepted_pairs(#[case] atom: AtomId, #[case] expected: ValueAst) {
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::N),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![(vec![AtomId(0)], AtomId(1), DativeBondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        );
        assert_eq!(molecule.atom(atom).accepted_pairs(), expected);
    }

    #[rstest]
    fn test_atom_view_accepted_pairs_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::accepted_pairs(2));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().accepted_pairs(),
            ValueAst::Lit(2),
        );
    }

    #[rstest]
    fn test_atom_view_aromatic_valence_not_in_system() {
        let molecule = mol!(r#"{:atoms ["C"] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).aromatic_valence(),
            ValueAst::Lit(0)
        );
    }

    #[rstest]
    #[case::in_system(AtomId(0), true)]
    #[case::not_in_system(AtomId(3), false)]
    fn test_atom_view_is_in_aromatic_system(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.atom(atom).is_in_aromatic_system(), expected);
    }

    #[rstest]
    #[case::participant(AtomId(0), Some(AromaticSystemId(0)))]
    #[case::not_participant(AtomId(3), None)]
    fn test_atom_view_aromatic_system(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Option<AromaticSystemId>,
    ) {
        let id = molecule.atom(atom).aromatic_system().map(|v| v.id);
        assert_eq!(id, expected);
    }

    #[rstest]
    #[case::donor(AtomId(2), vec![DativeBondId(0)])]
    #[case::acceptor(AtomId(3), vec![DativeBondId(0)])]
    #[case::uninvolved(AtomId(0), vec![])]
    fn test_atom_view_dative_bonds(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<DativeBondId>,
    ) {
        let ids: Vec<DativeBondId> = molecule.atom(atom).dative_bonds().map(|v| v.id).collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case::participant(AtomId(0), vec![MulticenterBondId(0)])]
    #[case::uninvolved(AtomId(3), vec![])]
    fn test_atom_view_multicenter_bonds(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<MulticenterBondId>,
    ) {
        let ids: Vec<MulticenterBondId> = molecule
            .atom(atom)
            .multicenter_bonds()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case::endpoint_0(AtomId(0), vec![NoncovalentBondId(0)])]
    #[case::endpoint_3(AtomId(3), vec![NoncovalentBondId(0)])]
    #[case::uninvolved(AtomId(1), vec![])]
    fn test_atom_view_noncovalent_bonds(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<NoncovalentBondId>,
    ) {
        let ids: Vec<NoncovalentBondId> = molecule
            .atom(atom)
            .noncovalent_bonds()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    /// Cyclohexane with one chain atom: 0-1-2-3-4-5-0 closing the ring, plus 0-6 dangling.
    #[fixture]
    fn ring_with_chain() -> MoleculeAst {
        MoleculeAst::from_atoms_and_bonds(
            vec![AtomAst::from_element(Element::C); 7],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
                (AtomId(3), AtomId(4), BondAst::from_order(1)),
                (AtomId(4), AtomId(5), BondAst::from_order(1)),
                (AtomId(5), AtomId(0), BondAst::from_order(1)),
                (AtomId(0), AtomId(6), BondAst::from_order(1)),
            ],
        )
    }

    #[rstest]
    #[case::ring_atom_0(AtomId(0), true)]
    #[case::ring_atom_3(AtomId(3), true)]
    #[case::ring_atom_5(AtomId(5), true)]
    #[case::chain_atom_6(AtomId(6), false)]
    fn test_atom_view_is_in_ring(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(ring_with_chain.atom(atom).is_in_ring(), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), true)]
    #[case::chain_atom(AtomId(6), false)]
    fn test_atom_view_is_in_ring_from(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        assert_eq!(ring_with_chain.atom(atom).is_in_ring_from(&rings), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), 1)]
    #[case::chain_atom(AtomId(6), 0)]
    fn test_atom_view_rings_from(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected_count: usize,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        let count = ring_with_chain.atom(atom).rings_from(&rings).count();
        assert_eq!(count, expected_count);
    }

    #[rstest]
    #[case::aromatic_and_multicenter(molecule(), AtomId(0), true)]
    #[case::aromatic_only_in_rich(molecule(), AtomId(1), true)]
    #[case::dative_donor(molecule(), AtomId(2), true)]
    #[case::dative_acceptor(molecule(), AtomId(3), true)]
    #[case::bare_atom_0(
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
        AtomId(0),
        false,
    )]
    #[case::bare_atom_1(
        MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
        ),
        AtomId(1),
        false,
    )]
    fn test_atom_view_is_in_overlays(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(mol.atom(atom).is_in_overlays(), expected);
    }

    #[rstest]
    fn test_atom_view_aromatic_valence_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::aromatic_valence(
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
        ));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().aromatic_valence(),
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
        );
    }

    #[rstest]
    #[case::single_bond(
        vec![(vec![AtomId(0), AtomId(1)], vec![ValueAst::Lit(2), ValueAst::Lit(2)])],
        ValueAst::Lit(2),
    )]
    #[case::two_bonds(
        vec![
            (vec![AtomId(0), AtomId(1)], vec![ValueAst::Lit(2), ValueAst::Lit(2)]),
            (vec![AtomId(0), AtomId(2)], vec![ValueAst::Lit(1), ValueAst::Lit(1)]),
        ],
        ValueAst::Lit(3),
    )]
    #[case::undetermined_aborts(
        vec![(vec![AtomId(0), AtomId(1)], vec![ValueAst::Undetermined, ValueAst::Lit(2)])],
        ValueAst::Undetermined,
    )]
    fn test_atom_view_multicenter_valence(
        #[case] bonds: Vec<(Vec<AtomId>, Vec<ValueAst>)>,
        #[case] expected: ValueAst,
    ) {
        let multicenter: Vec<_> = bonds
            .into_iter()
            .map(|(parts, electrons)| (parts, MulticenterBondAst::new(electrons)))
            .collect();
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![],
            vec![],
            multicenter,
            vec![],
            Constraints::default(),
        );
        assert_eq!(molecule.atom(AtomId(0)).multicenter_valence(), expected);
    }

    #[rstest]
    fn test_atom_view_multicenter_valence_constraint() {
        let mut atom = AtomAst::from_element(Element::C);
        atom.constraints.add(AtomConstraint::multicenter_valence(
            MulticenterValenceAst::Multicenter(ValueAst::Lit(2)),
        ));
        let molecule = MoleculeAst::from_atoms_and_bonds(vec![atom], vec![]);
        assert_eq!(
            molecule.atom(AtomId(0)).constraints().multicenter_valence(),
            MulticenterValenceAst::Multicenter(ValueAst::Lit(2)),
        );
    }

    #[rstest]
    #[case::ethane_carbon(mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#), AtomId(0), ValueAst::Lit(1))]
    #[case::ethene_carbon(mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#), AtomId(0), ValueAst::Lit(1))]
    #[case::three_bonds(mol!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]}"#), AtomId(0), ValueAst::Lit(3))]
    fn test_atom_view_degree(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(mol.atom(atom).degree(), expected);
    }

    #[rstest]
    fn test_atom_view_total_degree() {
        // Methane: 0 incident bonds in graph + implicit_h=4 + no multicenter.
        let molecule = mol!(r#"{:atoms ["C#h4"] :bonds []}"#);
        assert_eq!(molecule.atom(AtomId(0)).total_degree(), ValueAst::Lit(4),);
    }

    #[rstest]
    fn test_atom_view_total_degree_undetermined() {
        // implicit_hydrogens = Normal (placeholder) collapses to Undetermined.
        let molecule = mol!(r#"{:atoms ["C#h="] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).total_degree(),
            ValueAst::Undetermined,
        );
    }

    #[rstest]
    #[case::all_heavy(
        mol!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        AtomId(0),
        ValueAst::Lit(2),
    )]
    #[case::one_h_neighbor(
        mol!(r#"{:atoms ["C" "C" "H"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        AtomId(0),
        ValueAst::Lit(1),
    )]
    fn test_atom_view_heavy_atom_degree(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(mol.atom(atom).heavy_atom_degree(), expected);
    }

    #[rstest]
    #[case::all_heavy(
        mol!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [0 2 "2"]]}"#),
        AtomId(0),
        ValueAst::Lit(3),
    )]
    #[case::skips_h(
        mol!(r#"{:atoms ["C" "C" "H"] :bonds [[0 1 "2"] [0 2 "1"]]}"#),
        AtomId(0),
        ValueAst::Lit(2),
    )]
    fn test_atom_view_heavy_atom_valence(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(mol.atom(atom).heavy_atom_valence(), expected);
    }

    #[rstest]
    #[case::implicit_only(
        mol!(r#"{:atoms ["C#h4"] :bonds []}"#),
        AtomId(0),
        ValueAst::Lit(4),
    )]
    #[case::implicit_and_explicit(
        mol!(r#"{:atoms ["C#h2" "H" "H"] :bonds [[0 1 "1"] [0 2 "1"]]}"#),
        AtomId(0),
        ValueAst::Lit(4),
    )]
    #[case::implicit_normal_collapses(
        mol!(r#"{:atoms ["C#h="] :bonds []}"#),
        AtomId(0),
        ValueAst::Undetermined,
    )]
    fn test_atom_view_total_hydrogens(
        #[case] mol: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(mol.atom(atom).total_hydrogens(), expected);
    }

    #[rstest]
    fn test_atom_view_total_valence_sum_of_terms() {
        // Methane with implicit_h=4: valence=0, implicit=4, aromatic=0,
        // multicenter=0 → total=4.
        let molecule = mol!(r#"{:atoms ["C#h4"] :bonds []}"#);
        assert_eq!(molecule.atom(AtomId(0)).total_valence(), ValueAst::Lit(4),);
    }

    #[rstest]
    fn test_atom_view_total_valence_implicit_normal_collapses() {
        let molecule = mol!(r#"{:atoms ["C#h="] :bonds []}"#);
        assert_eq!(
            molecule.atom(AtomId(0)).total_valence(),
            ValueAst::Undetermined,
        );
    }

    #[rstest]
    fn test_atom_view_multicenter_degree() {
        // 3-atom multicenter bond: atom 0's multicenter_degree = co-participant
        // count = 2.
        let molecule = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![],
            vec![],
            vec![],
            vec![(
                vec![AtomId(0), AtomId(1), AtomId(2)],
                MulticenterBondAst::new(vec![ValueAst::Lit(2), ValueAst::Lit(2), ValueAst::Lit(2)]),
            )],
            vec![],
            Constraints::default(),
        );
        assert_eq!(
            molecule.atom(AtomId(0)).multicenter_degree(),
            ValueAst::Lit(2),
        );
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), ValueAst::Lit(1))]
    #[case::ring_atom_alt(AtomId(3), ValueAst::Lit(1))]
    #[case::chain_atom(AtomId(6), ValueAst::Lit(0))]
    fn test_atom_view_ring_count(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.atom(atom).ring_count(), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), vec![6])]
    #[case::chain_atom(AtomId(6), vec![])]
    fn test_atom_view_ring_size(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<usize>,
    ) {
        let sizes: Vec<_> = ring_with_chain.atom(atom).ring_size().collect();
        assert_eq!(sizes, expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), Some(6))]
    #[case::chain_atom(AtomId(6), None)]
    fn test_atom_view_smallest_ring_size(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Option<usize>,
    ) {
        assert_eq!(ring_with_chain.atom(atom).smallest_ring_size(), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), ValueAst::Lit(2))]
    #[case::chain_atom(AtomId(6), ValueAst::Lit(0))]
    fn test_atom_view_ring_degree(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.atom(atom).ring_degree(), expected);
    }

    #[rstest]
    #[case::ring_atom(AtomId(0), ValueAst::Lit(2))]
    #[case::chain_atom(AtomId(6), ValueAst::Lit(0))]
    fn test_atom_view_ring_valence(
        ring_with_chain: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.atom(atom).ring_valence(), expected);
    }

    // --- BondViews ---

    #[rstest]
    fn test_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.bonds().count(), 3);
    }

    #[rstest]
    fn test_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.bonds().ids().collect::<Vec<_>>(),
            vec![BondId(0), BondId(1), BondId(2)],
        );
    }

    #[rstest]
    fn test_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(BondId, [AtomId; 2], BondAst)> = molecule
            .bonds()
            .iter()
            .map(|v| (v.id, v.atom_ids(), v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(0), [AtomId(0), AtomId(1)], BondAst::from_order(1)),
                (BondId(1), [AtomId(1), AtomId(2)], BondAst::from_order(2)),
                (BondId(2), [AtomId(2), AtomId(3)], BondAst::from_order(1)),
            ],
        );
    }

    #[rstest]
    fn test_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.bonds().get(BondId(1));
        assert_eq!(view.id, BondId(1));
        assert_eq!(view.atom_ids(), [AtomId(1), AtomId(2)]);
        assert_eq!(*view.ast, BondAst::from_order(2));
    }

    #[rstest]
    fn test_bond_views_index(molecule: MoleculeAst) {
        let bond: &BondAst = &molecule.bonds()[BondId(1)];
        assert_eq!(*bond, BondAst::from_order(2));
    }

    // --- BondView ---

    #[rstest]
    fn test_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(molecule.bond(BondId(1)).atom_ids(), [AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_bond_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule.bond(BondId(1)).atoms().map(|a| a.id).collect();
        assert_eq!(ids, vec![AtomId(1), AtomId(2)]);
    }

    #[rstest]
    #[case::both_endpoints_aromatic(BondId(0), Some(AromaticSystemId(0)))]
    #[case::both_endpoints_aromatic_alt(BondId(1), Some(AromaticSystemId(0)))]
    #[case::one_endpoint_outside(BondId(2), None)]
    fn test_bond_view_aromatic_system(
        molecule: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: Option<AromaticSystemId>,
    ) {
        let id = molecule.bond(bond).aromatic_system().map(|v| v.id);
        assert_eq!(id, expected);
    }

    #[rstest]
    #[case::both_endpoints_aromatic(BondId(0), true)]
    #[case::both_endpoints_aromatic_alt(BondId(1), true)]
    #[case::one_endpoint_outside(BondId(2), false)]
    fn test_bond_view_is_in_aromatic_system(
        molecule: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.bond(bond).is_in_aromatic_system(), expected);
    }

    #[rstest]
    #[case::ring_bond_0_1(BondId(0), true)]
    #[case::ring_bond_5_0(BondId(5), true)]
    #[case::chain_bond_0_6(BondId(6), false)]
    fn test_bond_view_is_in_ring(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        assert_eq!(ring_with_chain.bond(bond).is_in_ring(), expected);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), true)]
    #[case::chain_bond(BondId(6), false)]
    fn test_bond_view_is_in_ring_from(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: bool,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        assert_eq!(ring_with_chain.bond(bond).is_in_ring_from(&rings), expected);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), 1)]
    #[case::chain_bond(BondId(6), 0)]
    fn test_bond_view_rings_from(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected_count: usize,
    ) {
        let rings = ring_with_chain.rings_with(RingFamily::Relevant, 22, |_| true);
        let count = ring_with_chain.bond(bond).rings_from(&rings).count();
        assert_eq!(count, expected_count);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), ValueAst::Lit(1))]
    #[case::chain_bond(BondId(6), ValueAst::Lit(0))]
    fn test_bond_view_ring_count(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: ValueAst,
    ) {
        assert_eq!(ring_with_chain.bond(bond).ring_count(), expected);
    }

    #[rstest]
    #[case::ring_bond(BondId(0), vec![6])]
    #[case::chain_bond(BondId(6), vec![])]
    fn test_bond_view_ring_size(
        ring_with_chain: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: Vec<usize>,
    ) {
        let sizes: Vec<_> = ring_with_chain.bond(bond).ring_size().collect();
        assert_eq!(sizes, expected);
    }

    // --- DativeBondViews ---

    #[rstest]
    fn test_dative_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bonds().count(), 1);
    }

    #[rstest]
    fn test_dative_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.dative_bonds().ids().collect::<Vec<_>>(),
            vec![DativeBondId(0)],
        );
    }

    #[rstest]
    fn test_dative_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(DativeBondId, AtomId, DativeBondAst)> = molecule
            .dative_bonds()
            .iter()
            .map(|v| (v.id, v.acceptor_id, v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                DativeBondId(0),
                AtomId(3),
                DativeBondAst::from_order(1).with_acceptor_slot(1),
            )],
        );
    }

    #[rstest]
    fn test_dative_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.dative_bonds().get(DativeBondId(0));
        assert_eq!(view.id, DativeBondId(0));
        assert_eq!(view.acceptor_id, AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_views_index(molecule: MoleculeAst) {
        let dative: &DativeBondAst = &molecule.dative_bonds()[DativeBondId(0)];
        assert_eq!(dative.order, ValueAst::Lit(1));
    }

    // --- DativeBondView ---

    #[rstest]
    fn test_dative_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2), AtomId(3)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_donor_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .dative_bond(DativeBondId(0))
                .donor_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2)],
        );
    }

    #[rstest]
    fn test_dative_bond_view_acceptor_id(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bond(DativeBondId(0)).acceptor_id, AtomId(3));
    }

    #[rstest]
    fn test_dative_bond_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .dative_bond(DativeBondId(0))
            .atoms()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(2), AtomId(3)]);
    }

    #[rstest]
    fn test_dative_bond_view_donors(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .dative_bond(DativeBondId(0))
            .donors()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(2)]);
    }

    #[rstest]
    fn test_dative_bond_view_acceptor(molecule: MoleculeAst) {
        assert_eq!(
            molecule.dative_bond(DativeBondId(0)).acceptor().id,
            AtomId(3),
        );
    }

    #[rstest]
    fn test_dative_bond_view_atom_count(molecule: MoleculeAst) {
        assert_eq!(molecule.dative_bond(DativeBondId(0)).atom_count(), 2);
    }

    // --- AromaticSystemViews ---

    #[rstest]
    fn test_aromatic_system_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.aromatic_systems().count(), 1);
    }

    #[rstest]
    fn test_aromatic_system_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.aromatic_systems().ids().collect::<Vec<_>>(),
            vec![AromaticSystemId(0)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(AromaticSystemId, Vec<AtomId>)> = molecule
            .aromatic_systems()
            .iter()
            .map(|v| (v.id, v.atom_ids().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(AromaticSystemId(0), vec![AtomId(0), AtomId(1), AtomId(2)])],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_get(molecule: MoleculeAst) {
        let view = molecule.aromatic_systems().get(AromaticSystemId(0));
        assert_eq!(view.id, AromaticSystemId(0));
        assert_eq!(
            view.atom_ids().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_aromatic_system_views_index(molecule: MoleculeAst) {
        let _: &AromaticSystemAst = &molecule.aromatic_systems()[AromaticSystemId(0)];
    }

    // --- AromaticSystemView ---

    #[rstest]
    fn test_aromatic_system_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bond_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .bond_ids()
                .collect::<Vec<_>>(),
            vec![BondId(0), BondId(1)],
        );
    }

    #[rstest]
    fn test_aromatic_system_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .atoms()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(0), AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_aromatic_system_view_bonds(molecule: MoleculeAst) {
        let ids: Vec<BondId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .bonds()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![BondId(0), BondId(1)]);
    }

    #[rstest]
    fn test_aromatic_system_view_electron_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .electron_count(),
            ValueAst::Lit(0),
        );
    }

    #[rstest]
    fn test_aromatic_system_view_atom_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule.aromatic_system(AromaticSystemId(0)).atom_count(),
            3
        );
    }

    #[rstest]
    fn test_aromatic_system_view_bond_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule.aromatic_system(AromaticSystemId(0)).bond_count(),
            2
        );
    }

    #[rstest]
    #[case::two_in(vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(1)])]
    #[case::all_in(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AtomId(0), AtomId(1), AtomId(2)])]
    #[case::disjoint(vec![AtomId(3)], vec![])]
    fn test_aromatic_system_view_overlapping_atoms(
        molecule: MoleculeAst,
        #[case] subset: Vec<AtomId>,
        #[case] expected: Vec<AtomId>,
    ) {
        let ids: Vec<AtomId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .overlapping_atoms(&subset)
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    #[case::one(vec![BondId(0)], vec![BondId(0)])]
    #[case::both(vec![BondId(0), BondId(1)], vec![BondId(0), BondId(1)])]
    #[case::other(vec![BondId(2)], vec![])]
    fn test_aromatic_system_view_overlapping_bonds(
        molecule: MoleculeAst,
        #[case] subset: Vec<BondId>,
        #[case] expected: Vec<BondId>,
    ) {
        let ids: Vec<BondId> = molecule
            .aromatic_system(AromaticSystemId(0))
            .overlapping_bonds(&subset)
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    #[rstest]
    fn test_aromatic_system_view_overlapping_rings(molecule: MoleculeAst) {
        let ids: Vec<usize> = molecule
            .aromatic_system(AromaticSystemId(0))
            .overlapping_rings()
            .map(|r| r.len())
            .collect();
        assert_eq!(ids, Vec::<usize>::new());
    }

    // --- MulticenterBondViews ---

    #[rstest]
    fn test_multicenter_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.multicenter_bonds().count(), 1);
    }

    #[rstest]
    fn test_multicenter_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.multicenter_bonds().ids().collect::<Vec<_>>(),
            vec![MulticenterBondId(0)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(MulticenterBondId, Vec<AtomId>)> = molecule
            .multicenter_bonds()
            .iter()
            .map(|v| (v.id, v.atom_ids().collect()))
            .collect();
        assert_eq!(
            collected,
            vec![(MulticenterBondId(0), vec![AtomId(0), AtomId(1), AtomId(2)],)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.multicenter_bonds().get(MulticenterBondId(0));
        assert_eq!(view.id, MulticenterBondId(0));
        assert_eq!(
            view.atom_ids().collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_views_index(molecule: MoleculeAst) {
        let _: &MulticenterBondAst = &molecule.multicenter_bonds()[MulticenterBondId(0)];
    }

    // --- MulticenterBondView ---

    #[rstest]
    fn test_multicenter_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2)],
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_atoms(molecule: MoleculeAst) {
        let ids: Vec<AtomId> = molecule
            .multicenter_bond(MulticenterBondId(0))
            .atoms()
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, vec![AtomId(0), AtomId(1), AtomId(2)]);
    }

    #[rstest]
    fn test_multicenter_bond_view_electron_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .multicenter_bond(MulticenterBondId(0))
                .electron_count(),
            ValueAst::Lit(0),
        );
    }

    #[rstest]
    fn test_multicenter_bond_view_atom_count(molecule: MoleculeAst) {
        assert_eq!(
            molecule.multicenter_bond(MulticenterBondId(0)).atom_count(),
            3,
        );
    }

    #[rstest]
    #[case::two_in(vec![AtomId(0), AtomId(1)], vec![AtomId(0), AtomId(1)])]
    #[case::all_in(vec![AtomId(0), AtomId(1), AtomId(2)], vec![AtomId(0), AtomId(1), AtomId(2)])]
    #[case::disjoint(vec![AtomId(3)], vec![])]
    fn test_multicenter_bond_view_overlapping_atoms(
        molecule: MoleculeAst,
        #[case] subset: Vec<AtomId>,
        #[case] expected: Vec<AtomId>,
    ) {
        let ids: Vec<AtomId> = molecule
            .multicenter_bond(MulticenterBondId(0))
            .overlapping_atoms(&subset)
            .map(|v| v.id)
            .collect();
        assert_eq!(ids, expected);
    }

    // --- NoncovalentBondViews ---

    #[rstest]
    fn test_noncovalent_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.noncovalent_bonds().count(), 1);
    }

    #[rstest]
    fn test_noncovalent_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.noncovalent_bonds().ids().collect::<Vec<_>>(),
            vec![NoncovalentBondId(0)],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_views_iter(molecule: MoleculeAst) {
        let collected: Vec<(NoncovalentBondId, [AtomId; 2], NoncovalentBondAst)> = molecule
            .noncovalent_bonds()
            .iter()
            .map(|v| (v.id, v.atom_ids(), v.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![(
                NoncovalentBondId(0),
                [AtomId(0), AtomId(3)],
                NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
            )],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_views_get(molecule: MoleculeAst) {
        let view = molecule.noncovalent_bonds().get(NoncovalentBondId(0));
        assert_eq!(view.id, NoncovalentBondId(0));
        assert_eq!(view.atom_ids(), [AtomId(0), AtomId(3)]);
    }

    #[rstest]
    fn test_noncovalent_bond_views_index(molecule: MoleculeAst) {
        let _: &NoncovalentBondAst = &molecule.noncovalent_bonds()[NoncovalentBondId(0)];
    }

    // --- NoncovalentBondView ---

    #[rstest]
    fn test_noncovalent_bond_view_atom_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.noncovalent_bond(NoncovalentBondId(0)).atom_ids(),
            [AtomId(0), AtomId(3)],
        );
    }

    #[rstest]
    fn test_noncovalent_bond_view_atoms(molecule: MoleculeAst) {
        let ids = molecule
            .noncovalent_bond(NoncovalentBondId(0))
            .atoms()
            .map(|v| v.id);
        assert_eq!(ids, [AtomId(0), AtomId(3)]);
    }

    // --- NeighborView ---

    #[rstest]
    fn test_neighbor_view_fields(molecule: MoleculeAst) {
        let collected: Vec<(BondId, AtomId, BondAst)> = molecule
            .neighbors(AtomId(2))
            .map(|n| (n.bond, n.atom, n.ast.clone()))
            .collect();
        assert_eq!(
            collected,
            vec![
                (BondId(1), AtomId(1), BondAst::from_order(2)),
                (BondId(2), AtomId(3), BondAst::from_order(1)),
            ],
        );
    }
}
