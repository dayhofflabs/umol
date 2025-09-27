//! Benchmarks for SMILES parsing

use std::hint::black_box;

use bstr::ByteVec;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rand_chacha::ChaCha12Rng;
use umol_models_graph::io::smiles::{
    parse_smiles_m0, parse_smiles_m1, parse_smiles_m2, parse_smiles_m3, parse_smiles_m4,
};

fn opensmiles_parsing(c: &mut Criterion) {
    // Mixed elements chain (organic-only mix omitting bare H)
    let mut rng = ChaCha12Rng::seed_from_u64(20250922);
    let mut mix_20 = Vec::from_slice(b"CNOFPSICNOFPSICNOFPS");
    mix_20.shuffle(&mut rng);
    let mut mix_50 = Vec::from_slice(b"CNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSIC");
    mix_50.shuffle(&mut rng);
    let mut mix_100 = Vec::from_slice(b"CNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICNOFPSICN");
    mix_100.shuffle(&mut rng);

    // Chain-only corpus, bare atoms (organic-only mix omitting bare H)
    let chain_inputs = [
        ("empty", &b""[..]),
        ("c_1", &b"C"[..]),
        ("c_5", &b"CCCCC"[..]),
        ("c_10", &b"CCCCCCCCCC"[..]),
        ("c_50", &b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_100", &b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_1000", &b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("mix_20", &mix_20[..]),
        ("mix_50", &mix_50[..]),
        ("mix_100", &mix_100[..])
    ];

    // Tree corpus, bare atoms
    let tree_inputs = [
        ("empty", &b"()"[..]),
        ("c_1", &b"(C)"[..]),
        ("c_5_1", &b"(CCCCC)"[..]),
        ("c_5_2", &b"C(C(C(C(C))))"[..]),
        ("c_5_3", &b"((((C)C)C)C)C"[..]),
        ("c_5_4", &b"C(CCC)C"[..]),
        ("c_5_5", &b"CC(C)CC"[..]),
        ("c_5_6", &b"CC(C)(C)C"[..]),
        ("c_10_1", &b"(CCCCCCCCCC)"[..]),
        ("c_10_2", &b"C(C(C(C(C(C(C(C(C(C)))))))))"[..]),
        ("c_10_3", &b"(((((((((C)C)C)C)C)C)C)C)C)C"[..]),
        ("c_10_4", &b"C(CCCCCCCC)C"[..]),
        ("c_10_5", &b"CCC(C)C(CC)CCC"[..]),
        ("c_10_6", &b"CC(C)(C)C(C)C(C)C"[..]),
        ("c_10_7", &b"CCC(C(C)C)CC(C)C"[..]),
        ("c_10_8", &b"C(C)C(C)C(C)C(C)C(C)CCCCC"[..]),
        ("c_50_r1", &b"C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_50_r5", &b"C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_50_r10", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_50_r25", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_100_r1", &b"C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_100_r5", &b"C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_100_r10", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_100_r25", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_100_r50", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_1000_r1", &b"C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_1000_r5", &b"C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_1000_r10", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_1000_r25", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_1000_r50", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_1000_r100", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
        ("c_1000_r250", &b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"[..]),
    ];

    // Chain bonds corpus, bare atoms
    let chain_bonds_inputs = [
        ("c_2", &b"C-C"[..]),
        ("c_2_stereo", &b"C/C"[..]),
        ("c_2_double", &b"C=C"[..]),
        ("c_2_triple", &b"C#C"[..]),
        ("c_2_quadruple", &b"C$C"[..]),
        ("c_5", &b"C=C=C=C=C"[..]),
        ("c_5_stereo", &b"C/C=C\\C=C/C"[..]),
        ("c_10", &b"C=C=C=C=C=C=C=C=C=C"[..]),
        ("c_10_aromatic", &b"C:C:C:C:C:C:C:C:C:C"[..]),
        ("c_10_stereo", &b"C/C=C\\C=C/C=C\\C=C/C"[..]),
        ("c_10_mixed", &b"C-C=C#C$C:C/C\\C-C=C"[..]),
        ("mixed_10", &b"C-N=C-Br-N=F-S-O-Cl=N"[..]),
        ("c_50", &b"C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C"[..]),
        ("c_50_mixed", &b"C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C"[..]),
        ("c_100", &b"C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C"[..]),
        ("c_100_mixed", &b"C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C"[..]),
        ("c_1000", &b"C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C"[..]),
    ];

    // Tree bonds corpus, bare atoms
    let tree_bonds_inputs = [
        ("c_5_1", &b"(C=C=C=C=C)"[..]),
        ("c_5_2", &b"C(-C(-C(-C(-C))))"[..]),
        ("c_5_3", &b"((((C)-C)-C)-C)-C"[..]),
        ("c_10_1", &b"(C/C=C\\C=C/C=C\\C=C/C)"[..]),
        ("c_10_2", &b"C(-C(#C(-C(#C(-C(#C(-C(#C(-C)))))))))"[..]),
        ("c_10_3", &b"(((((((((C)=C)=C)=C)=C)=C)=C)=C)=C)=C"[..]),
        ("c_10_4", &b"C(-C=C-C=C-C=C-C=C)=C"[..]),
        ("c_10_5", &b"C=CC(C)=C(C=C)C=CC"[..]),
        ("c_10_6", &b"CC(C)(C)C(=C)C(=C)C"[..]),
        ("c_10_7", &b"CCC(-C(-C)-C)C=C(-C)C"[..]),
        ("c_10_8", &b"C(C)=C(C)C(C)=C(C)C(C)=CC=CC=C"[..]),
        ("c_10_9", &b"(C=C=C=C=C=C=C=C=C=C)"[..]),
        ("c_50", &b"C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)"[..]),
        ("c_100", &b"C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)"[..]),
    ];

    // Ring corpus, simple cycles
    let ring_cycles_inputs = [
        ("c3", &b"C1CC1"[..]),
        ("c4", &b"C1CCC1"[..]),
        ("c5", &b"C1CCCC1"[..]),
        ("c6", &b"C1CCCCC1"[..]),
        ("c8", &b"C1CCCCCCC1"[..]),
        ("pct12", &b"C%12CCCC%12"[..]),
    ];

    // Ring corpus, fused and spiro examples
    let ring_fused_spiro_inputs = [
        ("fused_1", &b"C1CC2CCCCC2CC1"[..]),
        ("fused_2", &b"C1CCC2CCCC2CC1"[..]),
        ("spiro_1", &b"C1CCC2(CC1)CCC2"[..]),
        ("spiro_2", &b"C1CC2(C1)CC2C"[..]),
    ];

    // Ring corpus, directed closures and percent indices
    let ring_stereo_inputs = [
        ("dir_up_open", &b"C/1CC1"[..]),
        ("dir_up_close", &b"C1CC/1"[..]),
        ("dir_up_both", &b"C/1CC/1"[..]),
        ("dir_down_both", &b"C\\1CC\\1"[..]),
        ("pct_dir_up_open", &b"C/%12CC%12"[..]),
        ("pct_dir_up_close", &b"C%12CC/%12"[..]),
        ("pct_dir_down_both", &b"C\\%12CC\\%12"[..]),
    ];

    // // Lex-only baseline
    // let mut group_lex = c.benchmark_group("opensmiles_parsing/lex_only");
    // for (name, s) in inputs.iter() {
    //     group_lex.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
    //         b.iter(|| {
    //             let input = black_box(input);
    //             let mut n = 0usize;
    //             for tok in Lexer::new(input) {
    //                 let _ = tok;
    //                 n += 1;
    //             }
    //             std::hint::black_box(n);
    //         })
    //     });
    // }
    // group_lex.finish();

    // // Parse only
    // let mut group_parse = c.benchmark_group("opensmiles_parsing/parse_only");
    // for (name, s) in inputs.iter() {
    //     group_parse.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
    //         b.iter(|| {
    //             let input = black_box(input);
    //             let mut state = ParseState::default();
    //             let parser = MoleculeParser::new();
    //             let lexer = Lexer::new(input);
    //             let _ = parser.parse(&mut state, lexer);
    //         })
    //     });
    // }
    // group_parse.finish();

    // // Parser minimal (no IR, no diags, increment counter in every action)
    // let mut group_min = c.benchmark_group("opensmiles_parsing/parse_minimal");
    // for (name, s) in inputs.iter() {
    //     group_min.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
    //         b.iter(|| {
    //             let input = black_box(input);
    //             let mut state = ParseState::with_mode(ParserMode::Minimal);
    //             let parser = MoleculeParser::new();
    //             let lexer = Lexer::new(input);
    //             let _ = parser.parse(&mut state, lexer);
    //             std::hint::black_box(state.action_count);
    //         })
    //     });
    // }
    // group_min.finish();

    // FSM M0 chain-only
    let mut group_m0_chain = c.benchmark_group("opensmiles_parsing/parse_m0_chain");
    for (name, s) in chain_inputs.iter() {
        group_m0_chain.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m0(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m0_chain.finish();

    // FSM M1 chain
    let mut group_m1_chain = c.benchmark_group("opensmiles_parsing/parse_m1_chain");
    for (name, s) in chain_inputs.iter() {
        group_m1_chain.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m1(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m1_chain.finish();

    // FSM M1 branch
    let mut group_m1_tree = c.benchmark_group("opensmiles_parsing/parse_m1_tree");
    for (name, s) in tree_inputs.iter() {
        group_m1_tree.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m1(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m1_tree.finish();

    // FSM M2 chain
    let mut group_m2_chain = c.benchmark_group("opensmiles_parsing/parse_m2_chain");
    for (name, s) in chain_inputs.iter() {
        group_m2_chain.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m2(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m2_chain.finish();

    // FSM M2 tree
    let mut group_m2_tree = c.benchmark_group("opensmiles_parsing/parse_m2_tree");
    for (name, s) in tree_inputs.iter() {
        group_m2_tree.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m2(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m2_tree.finish();

    // FSM M2 chain bonds
    let mut group_m2_chain_bonds = c.benchmark_group("opensmiles_parsing/parse_m2_chain_bonds");
    for (name, s) in chain_bonds_inputs.iter() {
        group_m2_chain_bonds.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m2(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m2_chain_bonds.finish();

    // FSM M2 tree bonds
    let mut group_m2_tree_bonds = c.benchmark_group("opensmiles_parsing/parse_m2_tree_bonds");
    for (name, s) in tree_bonds_inputs.iter() {
        group_m2_tree_bonds.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m2(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m2_tree_bonds.finish();

    let mut group_m3_cycles = c.benchmark_group("opensmiles_parsing/parse_m3_cycles");
    for (name, s) in ring_cycles_inputs.iter() {
        group_m3_cycles.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m3(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m3_cycles.finish();

    let mut group_m3_fused = c.benchmark_group("opensmiles_parsing/parse_m3_fused_spiro");
    for (name, s) in ring_fused_spiro_inputs.iter() {
        group_m3_fused.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m3(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m3_fused.finish();

    let mut group_m3_stereo = c.benchmark_group("opensmiles_parsing/parse_m3_ring_stereo");
    for (name, s) in ring_stereo_inputs.iter() {
        group_m3_stereo.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m3(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m3_stereo.finish();

    // M4: components
    let component_inputs = [
        ("two", &b"CC.CC"[..]),
        ("three", &b"C.C.C"[..]),
        ("grouped", &b"C(C).C(C)"[..]),
        ("rings_across_digit", &b"C1.CC1"[..]),
        ("rings_across_pct", &b"C%12.CC%12"[..]),
        ("stereo_up", &b"C/1.CC/1"[..]),
        ("stereo_down", &b"C\\1.CC\\1"[..]),
        ("stereo_up_pct", &b"C/%12.CC/%12"[..]),
        ("stereo_down_pct", &b"C\\%12.CC\\%12"[..]),
    ];
    let mut group_m4_components = c.benchmark_group("opensmiles_parsing/parse_m4_components");
    for (name, s) in component_inputs.iter() {
        group_m4_components.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_m4(input);
                assert!(result.is_ok());
            })
        });
    }
    group_m4_components.finish();
}

criterion_group!(benches, opensmiles_parsing);
criterion_main!(benches);
