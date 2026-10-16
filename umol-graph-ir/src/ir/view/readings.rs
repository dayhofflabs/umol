use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_chem::element::Element;

use crate::ir::{
    AromaticValenceForm, AtomConstraintForm, AtomConstraintKey, AtomForm, AtomId,
    BondConstraintForm, BondForm, BondId, ElementForm, Molecule, MoleculeEntries, Normalize,
    NumForm, RingConfig, RingModel,
};

#[rstest]
#[case::absent(None, Some(AtomConstraintForm::Valence(NumForm::Undetermined)))]
#[case::top(
    Some(NumForm::Undetermined),
    Some(AtomConstraintForm::Valence(NumForm::Undetermined))
)]
#[case::literal(
    Some(NumForm::Lit(2)),
    Some(AtomConstraintForm::Valence(NumForm::Lit(2)))
)]
#[case::singleton(Some(NumForm::lit_set([2])), Some(AtomConstraintForm::Valence(NumForm::Lit(2))))]
#[case::bottom(Some(NumForm::lit_set([])), None)]
fn test_atom_constraints_view_asserted_normalization(
    #[case] assertion: Option<NumForm>,
    #[case] expected: Option<AtomConstraintForm>,
) {
    let mut atom = AtomForm::from_element(Element::C);
    if let Some(value) = assertion {
        atom.constraints.set(AtomConstraintForm::Valence(value));
    }
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom],
        ..Default::default()
    });
    assert_eq!(
        molecule
            .atom(AtomId(0))
            .constraints()
            .asserted(AtomConstraintKey::Valence),
        expected
    );
}

#[rstest]
#[case::singleton(vec![NumForm::lit_set([2])], Some(NumForm::Lit(2)))]
#[case::empty(vec![], Some(NumForm::Lit(0)))]
#[case::unknown(vec![NumForm::Undetermined], Some(NumForm::Undetermined))]
#[case::bottom(vec![NumForm::lit_set([])], None)]
#[case::unknown_then_bottom(vec![NumForm::Undetermined, NumForm::lit_set([])], None)]
#[case::bottom_then_unknown(vec![NumForm::lit_set([]), NumForm::Undetermined], None)]
fn test_atom_view_valence_normalization(
    #[case] orders: Vec<NumForm>,
    #[case] expected: Option<NumForm>,
) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); orders.len() + 1],
        bonds: orders
            .into_iter()
            .enumerate()
            .map(|(i, order)| {
                (
                    AtomId(0),
                    AtomId(i as u32 + 1),
                    BondForm {
                        order,
                        ..Default::default()
                    },
                )
            })
            .collect(),
        ..Default::default()
    });
    let atom = molecule.atom(AtomId(0));
    assert_eq!(atom.valence(), expected);
    assert_eq!(
        atom.constraints().derived(AtomConstraintKey::Valence),
        expected.map(AtomConstraintForm::Valence)
    );
}

#[rstest]
#[case::absent(
    None,
    false,
    Some(AromaticValenceForm::NotAromatic),
    Some(AromaticValenceForm::NotAromatic)
)]
#[case::vacuous(
    Some(AromaticValenceForm::Undetermined),
    false,
    Some(AromaticValenceForm::NotAromatic),
    Some(AromaticValenceForm::NotAromatic)
)]
#[case::bond_only(
    None,
    true,
    Some(AromaticValenceForm::NotAromatic),
    Some(AromaticValenceForm::Aromatic(NumForm::Undetermined))
)]
#[case::own_positive(
    Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
    true,
    Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
    Some(AromaticValenceForm::Aromatic(NumForm::Lit(1)))
)]
#[case::own_singleton(
    Some(AromaticValenceForm::Aromatic(NumForm::lit_set([1]))),
    true,
    Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
    Some(AromaticValenceForm::Aromatic(NumForm::Lit(1)))
)]
#[case::own_bottom(Some(AromaticValenceForm::Aromatic(NumForm::lit_set([]))), true, None, None)]
#[case::conflict(
    Some(AromaticValenceForm::NotAromatic),
    true,
    Some(AromaticValenceForm::NotAromatic),
    None
)]
fn test_atom_constraints_view_evidenced_complete_sources(
    #[case] own: Option<AromaticValenceForm>,
    #[case] bond_marked: bool,
    #[case] expected_asserted: Option<AromaticValenceForm>,
    #[case] expected: Option<AromaticValenceForm>,
) {
    let mut atom = AtomForm::from_element(Element::C);
    if let Some(value) = own {
        atom.constraints
            .set(AtomConstraintForm::AromaticValence(value));
    }
    let mut bond = BondForm::from_order(1);
    if bond_marked {
        bond.constraints.set(BondConstraintForm::aromatic(true));
    }
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom, AtomForm::from_element(Element::C)],
        bonds: vec![(AtomId(0), AtomId(1), bond)],
        ..Default::default()
    });
    let constraints = molecule.atom(AtomId(0)).constraints();
    assert_eq!(
        constraints.evidenced_complete(AtomConstraintKey::AromaticValence),
        expected.map(AtomConstraintForm::AromaticValence)
    );
    assert_eq!(
        constraints.asserted_complete(AtomConstraintKey::AromaticValence),
        expected_asserted.map(AtomConstraintForm::AromaticValence)
    );
    assert_eq!(
        constraints.derived(AtomConstraintKey::DonatedPairs),
        Some(AtomConstraintForm::DonatedPairs(NumForm::Undetermined))
    );
    assert_eq!(
        constraints.derived_complete(AtomConstraintKey::DonatedPairs),
        Some(AtomConstraintForm::DonatedPairs(NumForm::Lit(0)))
    );
}

#[rstest]
#[case::hydrogen(ElementForm::lit_set([Element::H]), Some(NumForm::Lit(0)), Some(NumForm::Lit(1)))]
#[case::unknown(
    ElementForm::Undetermined,
    Some(NumForm::Lit(1)),
    Some(NumForm::Undetermined)
)]
#[case::bottom(ElementForm::lit_set([]), None, None)]
fn test_atom_view_heavy_atom_degree_normalization(
    #[case] neighbor: ElementForm,
    #[case] heavy: Option<NumForm>,
    #[case] hydrogens: Option<NumForm>,
) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_implicit_hydrogens(0),
            AtomForm {
                element: neighbor,
                ..Default::default()
            },
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    let atom = molecule.atom(AtomId(0));
    assert_eq!(atom.heavy_atom_degree(), heavy);
    assert_eq!(atom.heavy_atom_valence(), heavy);
    assert_eq!(atom.total_hydrogens(), hydrogens);
}

#[rstest]
fn test_atom_constraints_view_derived_mutation() {
    let mut molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    molecule.bond_mut(BondId(0)).attributes_mut().order = NumForm::lit_set([]);
    assert_eq!(molecule.atom(AtomId(0)).valence(), None);
    assert_eq!(molecule.bond(BondId(0)).order(), &NumForm::lit_set([]));
    molecule.bond_mut(BondId(0)).attributes_mut().order = NumForm::lit_set([2]);
    assert_eq!(molecule.atom(AtomId(0)).valence(), Some(NumForm::Lit(2)));
    assert_eq!(molecule.bond(BondId(0)).order(), &NumForm::lit_set([2]));
    molecule.atom_mut(AtomId(1)).attributes_mut().charge = NumForm::lit_set([]);
    assert_eq!(molecule.atom(AtomId(0)).valence(), Some(NumForm::Lit(2)));
}

#[rstest]
#[case::literal(NumForm::Lit(2))]
#[case::singleton(NumForm::lit_set([2]))]
#[case::set(NumForm::lit_set([1, 2]))]
#[case::unknown(NumForm::Undetermined)]
fn test_atom_constraints_view_normalization(#[case] value: NumForm) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C).with_implicit_hydrogens(value.clone()); 3],
        bonds: vec![
            (
                AtomId(0),
                AtomId(1),
                BondForm {
                    order: value.clone(),
                    ..Default::default()
                },
            ),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(0), BondForm::from_order(1)),
        ],
        ..Default::default()
    });
    let normal = molecule.clone().normalize().unwrap();
    let rings = molecule
        .rings(RingModel::default(), RingConfig::default())
        .into_ring_set();
    let before = molecule.atom(AtomId(0)).constraints().with_rings(&rings);
    let after = normal.atom(AtomId(0)).constraints().with_rings(&rings);
    for key in [
        AtomConstraintKey::Valence,
        AtomConstraintKey::TotalDegree,
        AtomConstraintKey::TotalHydrogens,
        AtomConstraintKey::TotalValence,
        AtomConstraintKey::RingValence,
    ] {
        assert_eq!(before.asserted(key), after.asserted(key));
        assert_eq!(before.asserted_complete(key), after.asserted_complete(key));
        assert_eq!(before.derived(key), after.derived(key));
        assert_eq!(before.derived_complete(key), after.derived_complete(key));
        assert_eq!(before.evidenced(key), after.evidenced(key));
        assert_eq!(
            before.evidenced_complete(key),
            after.evidenced_complete(key)
        );
    }
}
