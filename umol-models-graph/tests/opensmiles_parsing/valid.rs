use rstest::rstest;

use super::fixtures::{rng, accepts};

fn atom_pool() -> &'static [&'static str] {
    &[
        // conservative organic subset only
        "C", "N", "O", "S", "P", "F", "Cl", "Br", "I",
    ]
}

fn bond_pool() -> &'static [&'static str] { &["-", "="] }

fn ring_indices() -> &'static [&'static str] { &["1", "2", "3", "%12"] }

fn gen_atom(rng: &mut fastrand::Rng) -> String {
    let a = atom_pool()[rng.usize(..atom_pool().len())];
    a.to_string()
}

fn maybe_bond(rng: &mut fastrand::Rng) -> Option<&'static str> {
    if rng.bool() { None } else { Some(bond_pool()[rng.usize(..bond_pool().len())]) }
}

fn maybe_branch_tokens(rng: &mut fastrand::Rng, depth: usize) -> Option<Vec<String>> {
    if depth == 0 || rng.u8(..10) < 7 { return None; }
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
    let mut atom_pos: Vec<usize> = Vec::new();

    for i in 0..len {
        if i > 0 {
            if let Some(b) = maybe_bond(rng) { tokens.push(b.to_string()); }
        }
        let a = gen_atom(rng);
        tokens.push(a);
        atom_pos.push(tokens.len() - 1);
        if let Some(branch) = maybe_branch_tokens(rng, depth) { tokens.extend(branch); }
    }

    // Initially skip rings to keep generator strictly valid; add later once stable

    tokens
}

fn gen_component(rng: &mut fastrand::Rng, depth: usize) -> String {
    gen_chain_tokens(rng, depth, 1).join("")
}

fn gen_molecule(rng: &mut fastrand::Rng, depth: usize, width: usize) -> String {
    let comps = 1 + rng.usize(..=width);
    let mut parts: Vec<String> = Vec::new();
    for _ in 0..comps { parts.push(gen_component(rng, depth)); }
    parts.join(".")
}

#[rstest]
fn valid_generated_accepts(mut rng: fastrand::Rng) {
    let mut accepted = 0usize;
    for _ in 0..5000 {
        let s = gen_molecule(&mut rng, 2, 1);
        if accepts(&s) { accepted += 1; } else {
            panic!("Generated valid SMILES was rejected: {}", s);
        }
    }
    assert!(accepted > 0);
}
