use pretty_assertions::assert_eq;
use rstest::*;

use super::super::*;

#[rstest]
#[case::wildcard_bare(b"*", 1, 0)]
#[case::wildcard_in_chain(b"C*C", 3, 2)]
#[case::wildcard_branch(b"C(*)C", 3, 2)]
#[case::wildcard_with_class(b"[*:1]", 1, 0)]
#[case::wildcard_bonded(b"C-*", 2, 1)]
#[case::multiple_wildcards(b"*.*", 2, 0)]
fn extended_wildcard(#[case] input: &[u8], #[case] atoms: usize, #[case] bonds: usize) {
    let res = parse_extended_smiles_bytes(input);
    assert!(res.is_ok(), "{:?} should have succeeded: {:?}", input, res);
    let mol = res.unwrap();
    assert_eq!(
        mol.atom_count(),
        atoms,
        "atom count mismatch for {:?}",
        input
    );
    assert_eq!(
        mol.bond_count(),
        bonds,
        "bond count mismatch for {:?}",
        input
    );
}
