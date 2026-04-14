//! VF2 substructure matching benchmark against 9k SMILES corpus.
//!
//! Three query patterns of increasing complexity:
//! 1. Branched: C(C)C(C)N  — 5 atoms, 4 single bonds
//! 2. Phenol:   c1ccccc1O  — 7 atoms, 7 any-bonds (6-ring + O)
//! 3. Bicyclic: C1~C(~C~C~C2)~C2~C~C~C1 — 9 atoms, 10 any-bonds (fused 5-6)

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};

use umol_graph::ast::atom::AtomAst;
use umol_graph::ast::bond::BondAst;
use umol_graph::ast::config::MoleculeAstConfig;
use umol_graph::ast::matcher::{find_matches, MatchQuery, MatchTarget};
use umol_graph::ast::molecule::{BondTuple, MoleculeAst};
use umol_graph::ast::ToAst;
use umol_graph::graph_ir::molecule_builder::MoleculeBuilder;
use umol_graph::io::smiles::parse_smiles;
use umol_shared::element::Element;
use umol_shared::value_ast::ValueAst;

fn load_smiles() -> Vec<String> {
    let data_dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/smiles_parsing/data/basic_opensmiles"
    );
    let mut smiles_list = Vec::new();
    for entry in walkdir::WalkDir::new(data_dir)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| {
            e.path()
                .extension()
                .map_or(false, |ext| ext == "smiles")
        })
    {
        let content = std::fs::read_to_string(entry.path()).unwrap();
        if let Some(s) = content.lines().nth(1) {
            if !s.is_empty() {
                smiles_list.push(s.to_string());
            }
        }
    }
    smiles_list
}

fn load_corpus(smiles_list: &[String]) -> Vec<MoleculeAst> {
    let cfg = MoleculeAstConfig::zeroed();
    let mut asts = Vec::new();
    for s in smiles_list {
        if let Ok(table_mol) = parse_smiles(s) {
            let builder = MoleculeBuilder::from_table_molecule(&table_mol);
            asts.push(builder.to_ast(&cfg));
        }
    }
    asts
}

fn wb() -> BondAst {
    BondAst::new(ValueAst::Undetermined)
}

fn sb() -> BondAst {
    BondAst::from_order(1)
}

fn a(e: Element) -> AtomAst {
    AtomAst::from_element(e)
}

fn bt(source: usize, target: usize, bond: BondAst) -> BondTuple {
    BondTuple {
        source,
        target,
        bond,
    }
}

// C(C)C(C)N — 5 atoms, 4 single bonds
fn pattern_branched() -> MoleculeAst {
    MoleculeAst {
        atoms: vec![a(Element::C), a(Element::C), a(Element::C), a(Element::C), a(Element::N)],
        bonds: vec![bt(0, 1, wb()), bt(0, 2, wb()), bt(2, 3, wb()), bt(2, 4, wb())],
        ..Default::default()
    }
}

// 6-membered C ring + O, any bonds
fn pattern_phenol() -> MoleculeAst {
    MoleculeAst {
        atoms: vec![
            a(Element::C), a(Element::C), a(Element::C),
            a(Element::C), a(Element::C), a(Element::C),
            a(Element::O),
        ],
        bonds: vec![
            bt(0, 1, wb()), bt(1, 2, wb()), bt(2, 3, wb()),
            bt(3, 4, wb()), bt(4, 5, wb()), bt(5, 0, wb()),
            bt(5, 6, wb()),
        ],
        ..Default::default()
    }
}

// Fused 5-6 bicyclic, all C, any bonds
// C1~C(~C~C~C2)~C2~C~C~C1
// Ring 1: 0-1-5-6-7-8 (6-membered)
// Ring 2: 1-2-3-4-5   (5-membered)
// Fused edge: 1-5
fn pattern_bicyclic() -> MoleculeAst {
    MoleculeAst {
        atoms: (0..9).map(|_| a(Element::C)).collect(),
        bonds: vec![
            bt(0, 1, wb()), bt(1, 2, wb()), bt(2, 3, wb()),
            bt(3, 4, wb()), bt(4, 5, wb()), bt(1, 5, wb()),
            bt(5, 6, wb()), bt(6, 7, wb()), bt(7, 8, wb()),
            bt(8, 0, wb()),
        ],
        ..Default::default()
    }
}

fn substructure_benchmark(c: &mut Criterion) {
    let smiles_list = load_smiles();
    let asts = load_corpus(&smiles_list);
    let corpus_size = asts.len();

    let targets: Vec<MatchTarget> = asts.iter().map(MatchTarget::new).collect();

    let patterns: Vec<(&str, MoleculeAst)> = vec![
        ("branched", pattern_branched()),
        ("phenol", pattern_phenol()),
        ("bicyclic", pattern_bicyclic()),
    ];

    let mut group = c.benchmark_group("substructure_vf2");

    for (name, pattern) in &patterns {
        let query = MatchQuery::new(pattern);

        let match_count = targets
            .iter()
            .filter(|t| !find_matches(&query, t).is_empty())
            .count();

        group.bench_function(
            BenchmarkId::new(*name, format!("{corpus_size}_hits_{match_count}")),
            |b| {
                b.iter(|| {
                    for target in &targets {
                        black_box(find_matches(&query, target));
                    }
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, substructure_benchmark);
criterion_main!(benches);
