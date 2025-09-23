use umol_models_graph::io::smiles::{parse_smiles_m0, M0Error};

#[test]
fn m0_accepts_chain_c() {
    let mol = parse_smiles_m0(b"CCCCC").expect("m0 should parse chain");
    assert_eq!(mol.atoms.len(), 5);
    assert_eq!(mol.bonds.len(), 4);
}

#[test]
fn m0_accepts_halogen_pairs() {
    let mol = parse_smiles_m0(b"CClCBrC").expect("m0 should parse halogens");
    assert_eq!(mol.atoms.len(), 5);
    assert_eq!(mol.bonds.len(), 4);
}

#[test]
fn m0_rejects_brackets() {
    let err = parse_smiles_m0(b"[C]").unwrap_err();
    match err { M0Error::UnsupportedToken { .. } => {}, _ => panic!("unexpected error") }
}

#[test]
fn m0_rejects_rings() {
    let err = parse_smiles_m0(b"C1CC1").unwrap_err();
    match err { M0Error::UnsupportedToken { .. } => {}, _ => panic!("unexpected error") }
}


