//! End-to-end SMILES ingest benchmarks — parse, raise, and resolve — over
//! representative single molecules: an unbranched alkane (localized valence
//! only), benzene and pyridine (joint aromatic selection), bare methane
//! (plural admission collapsed by the tie-break), and the bare-nitrogen
//! fused heteroaromatics quinoline and purine (assignment search over a
//! flexible fused component).
//!
//! The `smiles_roundtrip/resolve` group isolates resolution with the SMILES counts
//! model. Parsing, raising, resolver construction, and input cloning are outside timing.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, Criterion};
use umol_graph::ingest::ingest_smiles;
use umol_graph::ops::model::{ChemistryModel, ValenceModel};
use umol_graph::ops::resolve::Resolver;
use umol_graph_ir::ir::{Molecule, TryIntoIr};
use umol_io::smiles::Smiles;
use umol_utils::solution::Solution;

fn bench_ingest_smiles(c: &mut Criterion) {
    let mut g = c.benchmark_group("resolve_ingest_smiles");

    for (name, smiles) in [
        ("methane", "C"),
        ("octane", "CCCCCCCC"),
        ("benzene", "c1ccccc1"),
        ("pyridine", "c1ccncc1"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("quinoline", "c1ccc2ccccc2n1"),
        ("purine", "c1ncc2ncnc2n1"),
    ] {
        g.bench_function(name, |b| {
            b.iter(|| ingest_smiles(black_box(smiles)).unwrap())
        });
    }

    g.finish();
}

fn bench_resolve(c: &mut Criterion) {
    let model = ChemistryModel {
        valence: ValenceModel::smiles(),
        ..ChemistryModel::default()
    };
    let resolver = Resolver::new(&model);
    let mut group = c.benchmark_group("smiles_roundtrip/resolve");
    for (name, input) in [
        (
            "chain_64",
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
        ),
        ("branched", "CC(C)(F)C(Cl)CC"),
        ("components", "CC.O.[Na+]"),
        ("aromatic_fused", "c1ccc2ccccc2c1"),
        ("aromatic_lone_pair", "[nH]1cccc1"),
        ("charged", "[NH4+]"),
        ("radical", "[CH3]"),
        ("tetra_four", "[C@](F)(Cl)(Br)I"),
        ("tetra_ring", "C[C@H]1CCCCO1"),
        ("tetra_explicit_h", "C[C@]1([H])CCCCO1"),
        ("tetra_lone_pair", "C[S@](=O)CC"),
        ("alkene_four", "F/C(Cl)=C(/Br)I"),
        ("shared_chain", "C/C=C/C=C/C"),
        ("partial_triene", "C/C=C/C=CC=C/C"),
        ("shared_branch", "F/C=C(/C=C/F)C=C"),
        ("shared_cycle", "C1/C=C/C=C/CCC1"),
    ] {
        let table = Smiles::parse(input).unwrap().into_table_ir();
        let molecule: Molecule = (&table).try_into_ir(&()).unwrap();
        let mut checked = molecule.clone();
        assert!(
            matches!(
                resolver.resolve(&mut checked).unwrap(),
                Solution::Determined(_)
            ),
            "{name}"
        );
        group.bench_function(name, |b| {
            b.iter_batched_ref(
                || molecule.clone(),
                |molecule| resolver.resolve(black_box(molecule)).unwrap(),
                BatchSize::SmallInput,
            );
        });
    }
    group.finish();
}

criterion_group!(resolve, bench_ingest_smiles, bench_resolve);
criterion_main!(resolve);
