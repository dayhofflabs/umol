use rstest::rstest;
use umol_graph_ir::ir::{
    AromaticSystemId, AtomForm, AtomId, BondForm, BondId, Molecule, StereoAtomId, StereoBondId,
};
use umol_graph_ir::mol;

#[test]
fn test_mol_builds_molecule() {
    let molecule: Molecule = mol! {
        (c1: C) - (c2: C) = (o: O),
        (c1) - (n: N),
    };

    assert_eq!(molecule.atoms().count(), 4);
    assert_eq!(molecule.bonds().count(), 3);
    // bonds in declaration order: c1-c2 single, c2=o double, c1-n single
    assert_eq!(
        molecule.bond(BondId(0)).attributes,
        &BondForm::from_order(1)
    );
    assert_eq!(
        molecule.bond(BondId(1)).attributes,
        &BondForm::from_order(2)
    );
    assert_eq!(
        molecule.bond(BondId(2)).attributes,
        &BondForm::from_order(1)
    );
    assert_eq!(molecule.bond(BondId(1)).atom_ids(), [AtomId(1), AtomId(2)]);
}

#[test]
fn test_mol_quoted_spec() {
    // a rich DSL spec rides a string literal
    let molecule: Molecule = mol! {
        (c: "C#h3") - (n: N),
    };

    assert_eq!(molecule.atoms().count(), 2);
    assert_eq!(molecule.bonds().count(), 1);
}

#[test]
fn test_mol_bond_spec() {
    // a rich bond via the DSL spec: order 1 + aromatic flag
    let molecule: Molecule = mol! {
        (c1: C) -[ "1#a" ]- (c2: C),
    };

    assert_eq!(molecule.bonds().count(), 1);
    assert_eq!(
        molecule.bond(BondId(0)).attributes,
        &"1#a".parse::<BondForm>().unwrap()
    );
}

#[rstest]
fn test_mol_named_bond() {
    // a `-[name: "spec"]-` bond binds a label (inert until referenced) and carries the spec
    let molecule: Molecule = mol!((c: C) -[ b: "2" ]- (o: O));

    assert_eq!(molecule.bonds().count(), 1);
    assert_eq!(
        molecule.bond(BondId(0)).attributes,
        &"2".parse::<BondForm>().unwrap()
    );
}

#[rstest]
fn test_mol_anonymous_atoms() {
    // bare element idents introduce anonymous, unreferenceable atoms
    let molecule: Molecule = mol!(C - O);

    assert_eq!(molecule.atoms().count(), 2);
    assert_eq!(molecule.bonds().count(), 1);
    assert_eq!(
        molecule.atom(AtomId(0)).attributes,
        &"C".parse::<AtomForm>().unwrap()
    );
    assert_eq!(
        molecule.atom(AtomId(1)).attributes,
        &"O".parse::<AtomForm>().unwrap()
    );
    assert_eq!(molecule.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
}

#[rstest]
fn test_mol_anonymous_spec() {
    // the intended use: anonymous atoms carrying rich DSL specs
    let molecule: Molecule = mol!("C#h3" - "O#h");

    assert_eq!(
        molecule.atom(AtomId(0)).attributes,
        &"C#h3".parse::<AtomForm>().unwrap()
    );
    assert_eq!(
        molecule.atom(AtomId(1)).attributes,
        &"O#h".parse::<AtomForm>().unwrap()
    );
    assert_eq!(
        molecule.bond(BondId(0)).attributes,
        &BondForm::from_order(1)
    );
}

#[rstest]
fn test_mol_anonymous_mixed() {
    // a named atom referenced across paths, wired to anonymous terminals by position
    let molecule: Molecule = mol! {
        (c: C) = O,
        (c) - N,
    };

    assert_eq!(molecule.atoms().count(), 3);
    assert_eq!(molecule.bonds().count(), 2);
    // c=O double: position 0 (c) to position 1 (O)
    assert_eq!(
        molecule.bond(BondId(0)).attributes,
        &BondForm::from_order(2)
    );
    assert_eq!(molecule.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    // c-N single: position 0 (c) to position 2 (N)
    assert_eq!(molecule.bond(BondId(1)).atom_ids(), [AtomId(0), AtomId(2)]);
}

#[rstest]
fn test_mol_aromatic() {
    let molecule: Molecule = mol! {
        (c1: C) - (c2: C),
        aromatic [(c1) (c2)],
    };

    assert_eq!(molecule.aromatic_systems().count(), 1);
    assert_eq!(
        molecule
            .aromatic_system(AromaticSystemId(0))
            .atoms()
            .map(|view| view.id)
            .collect::<Vec<_>>(),
        vec![AtomId(0), AtomId(1)]
    );
}

#[rstest]
fn test_mol_dative() {
    let molecule: Molecule = mol! {
        (n: N), (b: B),
        dative [(n)] (b),
    };

    assert_eq!(molecule.dative_bonds().count(), 1);
}

#[rstest]
fn test_mol_multicenter() {
    let molecule: Molecule = mol! {
        (b1: B), (b2: B), (h: H),
        multicenter [(b1) (b2) (h)],
    };

    assert_eq!(molecule.multicenter_bonds().count(), 1);
}

#[rstest]
fn test_mol_noncovalent() {
    let molecule: Molecule = mol! {
        (o: O), (h: H),
        noncovalent [(o) (h)],
    };

    assert_eq!(molecule.noncovalent_bonds().count(), 1);
}

#[rstest]
fn test_mol_stereo_atom() {
    let molecule: Molecule = mol! {
        (c: C), (f: F), (cl: Cl), (br: Br), (i: I),
        stereo atom (c) [(f) (cl) (br) (i)] : "Th0",
    };

    assert_eq!(molecule.stereo_atoms().count(), 1);
    assert_eq!(molecule.stereo_atom(StereoAtomId(0)).site_id(), AtomId(0));
}

#[rstest]
fn test_mol_stereo_bond() {
    let molecule: Molecule = mol! {
        (c1: C) -[db: "2"]- (c2: C), (f: F), (h: H),
        stereo bond (db) [(f) (h) (f) (h)] : "Ct1",
    };

    assert_eq!(molecule.stereo_bonds().count(), 1);
    assert_eq!(molecule.stereo_bond(StereoBondId(0)).site_id(), BondId(0));
}
