use rstest::rstest;

use umol_ast::ast::{AtomAst, AtomId, BondAst, BondId, MoleculeAst};
use umol_ast::mol;

#[test]
fn test_mol_builds_molecule() {
    let molecule: MoleculeAst = mol! {
        (c1: C) - (c2: C) = (o: O),
        (c1) - (n: N),
    };

    assert_eq!(molecule.atoms().count(), 4);
    assert_eq!(molecule.bonds().count(), 3);
    // bonds in declaration order: c1-c2 single, c2=o double, c1-n single
    assert_eq!(molecule.bond(BondId(0)).ast, &BondAst::from_order(1));
    assert_eq!(molecule.bond(BondId(1)).ast, &BondAst::from_order(2));
    assert_eq!(molecule.bond(BondId(2)).ast, &BondAst::from_order(1));
    assert_eq!(molecule.bond(BondId(1)).atom_ids(), [AtomId(1), AtomId(2)]);
}

#[test]
fn test_mol_quoted_spec() {
    // a rich DSL spec rides a string literal
    let molecule: MoleculeAst = mol! {
        (c: "C#h3") - (n: N),
    };

    assert_eq!(molecule.atoms().count(), 2);
    assert_eq!(molecule.bonds().count(), 1);
}

#[test]
fn test_mol_bond_spec() {
    // a rich bond via the DSL spec: order 1 + aromatic flag
    let molecule: MoleculeAst = mol! {
        (c1: C) -[ "1#a" ]- (c2: C),
    };

    assert_eq!(molecule.bonds().count(), 1);
    assert_eq!(
        molecule.bond(BondId(0)).ast,
        &"1#a".parse::<BondAst>().unwrap()
    );
}

#[rstest]
fn test_mol_anonymous_atoms() {
    // bare element idents introduce anonymous, unreferenceable atoms
    let molecule: MoleculeAst = mol!(C - O);

    assert_eq!(molecule.atoms().count(), 2);
    assert_eq!(molecule.bonds().count(), 1);
    assert_eq!(molecule.atom(AtomId(0)).ast, &"C".parse::<AtomAst>().unwrap());
    assert_eq!(molecule.atom(AtomId(1)).ast, &"O".parse::<AtomAst>().unwrap());
    assert_eq!(molecule.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
}

#[rstest]
fn test_mol_anonymous_spec() {
    // the intended use: anonymous atoms carrying rich DSL specs
    let molecule: MoleculeAst = mol!("C#h3" - "O#h");

    assert_eq!(
        molecule.atom(AtomId(0)).ast,
        &"C#h3".parse::<AtomAst>().unwrap()
    );
    assert_eq!(
        molecule.atom(AtomId(1)).ast,
        &"O#h".parse::<AtomAst>().unwrap()
    );
    assert_eq!(molecule.bond(BondId(0)).ast, &BondAst::from_order(1));
}

#[rstest]
fn test_mol_anonymous_mixed() {
    // a named atom referenced across paths, wired to anonymous terminals by position
    let molecule: MoleculeAst = mol! {
        (c: C) = O,
        (c) - N,
    };

    assert_eq!(molecule.atoms().count(), 3);
    assert_eq!(molecule.bonds().count(), 2);
    // c=O double: position 0 (c) to position 1 (O)
    assert_eq!(molecule.bond(BondId(0)).ast, &BondAst::from_order(2));
    assert_eq!(molecule.bond(BondId(0)).atom_ids(), [AtomId(0), AtomId(1)]);
    // c-N single: position 0 (c) to position 2 (N)
    assert_eq!(molecule.bond(BondId(1)).atom_ids(), [AtomId(0), AtomId(2)]);
}
