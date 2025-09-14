//! Valid tests for OpenSMILES (UMOL)

use rstest::rstest;

use super::fixtures::{parse_and_assert_invariants, rng};

fn atom_pool() -> &'static [&'static str] {
    &[
        // conservative aliphatic organic subset
        "C", "N", "O", "S", "P", "F", "Cl", "Br", "I",
    ]
}

fn bond_pool() -> &'static [&'static str] {
    &["-", "="]
}

fn gen_atom(rng: &mut fastrand::Rng) -> String {
    let a = atom_pool()[rng.usize(..atom_pool().len())];
    a.to_string()
}

fn maybe_bond(rng: &mut fastrand::Rng) -> Option<&'static str> {
    if rng.u8(..10) < 6 {
        None
    } else {
        Some(bond_pool()[rng.usize(..bond_pool().len())])
    }
}

#[rstest]
fn valid_acceptance_set() {
    // c1ccccc1-c2ccccc2 (aliphatic single bond between aromatic rings)
    if let Err(e) = parse_and_assert_invariants("c1ccccc1-c2ccccc2") { panic!("{}", e); }

    // [*] fields: charge, chirality, H count, class
    for s in ["[*-]", "[*@H:0]", "[*H2]", "[*+:1]"] {
        if let Err(e) = parse_and_assert_invariants(s) { panic!("{}", e); }
    }

    // [H+] valid; [H][H] valid; bridging H examples
    for s in ["[H+]", "[H][H]", "[BH2]1[H][BH2][H]1"] {
        if let Err(e) = parse_and_assert_invariants(s) { panic!("{}", e); }
    }

    // Named and numeric isotopes
    for s in ["[2H]", "[3H]", "[13C]"] {
        if let Err(e) = parse_and_assert_invariants(s) { panic!("{}", e); }
    }

    // [H][CH3] valid
    if let Err(e) = parse_and_assert_invariants("[H][CH3]") { panic!("{}", e); }
}

fn maybe_branch_tokens(rng: &mut fastrand::Rng, depth: usize) -> Option<Vec<String>> {
    if depth == 0 || rng.u8(..10) < 7 {
        return None;
    }
    let inner = gen_chain_tokens(rng, depth - 1, 1);
    let mut v = Vec::with_capacity(inner.len() + 2);
    v.push("(".to_string());
    v.extend(inner);
    v.push(")".to_string());
    Some(v)
}

fn gen_chain_tokens(rng: &mut fastrand::Rng, depth: usize, min_len: usize) -> Vec<String> {
    let len = min_len + rng.usize(..=2); // 1..=3 nodes
    let mut tokens: Vec<String> = Vec::new();

    for i in 0..len {
        if i > 0 {
            if let Some(b) = maybe_bond(rng) {
                tokens.push(b.to_string());
            }
        }
        let a = gen_atom(rng);
        tokens.push(a);
        if let Some(branch) = maybe_branch_tokens(rng, depth) {
            tokens.extend(branch);
        }
    }

    tokens
}

fn gen_simple_ring(rng: &mut fastrand::Rng) -> String {
    // Generate a simple aliphatic ring: C1 C{n} 1 with n in [2..5] (size >= 3)
    let n = 2 + rng.usize(..=3);
    let mut s = String::from("C1");
    for _ in 0..n { s.push('C'); }
    s.push('1');
    s
}

fn gen_component(rng: &mut fastrand::Rng, depth: usize) -> String {
    // Occasionally produce a known-safe ring structure
    if rng.u8(..10) < 2 {
        return gen_simple_ring(rng);
    }
    gen_chain_tokens(rng, depth, 1).join("")
}

fn gen_molecule(rng: &mut fastrand::Rng, depth: usize, width: usize) -> String {
    let comps = 1 + rng.usize(..=width);
    let mut parts: Vec<String> = Vec::new();
    for _ in 0..comps {
        parts.push(gen_component(rng, depth));
    }
    parts.join(".")
}

#[rstest]
fn valid_unknown_atom_forms() {
    for s in ["*", "C*C", "*.*", "[*]", "[*H]", "[*H2]", "[*-]", "[*+:1]"] {
        if let Err(e) = parse_and_assert_invariants(s) {
            panic!("{}", e);
        }
    }
}

#[rstest]
fn valid_star_in_aromatic_contexts() {
    for s in [
        // star in aromatic ring (semantic aromaticity not enforced here)
        "c1*cccc1",
        // star adjacent to aromatic atoms with explicit ':' bonds
        "c:*:c",
        // bracketed star with fields inside aromatic ring
        "c1[*H+:2]cccc1",
    ] {
        if let Err(e) = parse_and_assert_invariants(s) {
            panic!("{}", e);
        }
    }
}

#[rstest]
fn valid_aromatic_ring_minimal() {
    // R25: basic aromatic ring with lowercase tokens should be accepted
    if let Err(e) = parse_and_assert_invariants("c1ccc1") {
        panic!("{}", e);
    }
}

#[rstest]
fn valid_specific_cases() {
    // H0 allowed in brackets
    // Using literal strings to avoid temporary String lifetimes
    if let Err(e) = parse_and_assert_invariants("[CH0]") { panic!("{}", e); }
    if let Err(e) = parse_and_assert_invariants("[NH0]") { panic!("{}", e); }
    // ring index 0 valid
    for s in ["C0CCCCC0"] {
        if let Err(e) = parse_and_assert_invariants(s) { panic!("{}", e); }
    }
    // Percent ring indices close with matching %NN occurrences
    if let Err(e) = parse_and_assert_invariants("C%12CCCC%12") { panic!("{}", e); }
    // C1CCC%01 valid (1 == %01)
    for s in ["C1CCC%01"] {
        if let Err(e) = parse_and_assert_invariants(s) { panic!("{}", e); }
    }
    // Multiple ring closures per atom (various forms)
    for s in [
        "C1CC12CC2",
        "C1CC1%10CC%10",
        "C%10CC%101CC1",
        "C%10CC%10%11CC%11",
    ] {
        if let Err(e) = parse_and_assert_invariants(s) { panic!("{}", e); }
    }
}

#[rstest]
fn valid_generated_accepts_and_invariants(mut rng: fastrand::Rng) {
    for _ in 0..5000 {
        let s = gen_molecule(&mut rng, 2, 1);
        if let Err(e) = parse_and_assert_invariants(&s) {
            panic!("{}", e);
        }
    }
}
