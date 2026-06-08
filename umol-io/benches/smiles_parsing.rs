//! Benchmarks for SMILES parsing

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use umol_io::smiles::parse_extended_smiles_bytes;
use umol_io::smiles::parser::parse_smiles_bytes_to_table_ir;

// Chain-only corpus, bare atoms (organic-only mix omitting bare H)
fn chain_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("empty", b""),
        ("c_1", b"C"),
        ("c_5", b"CCCCC"),
        ("c_10", b"CCCCCCCCCC"),
        // ("c_50", b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_100", b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_1000", b"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("mix_20", b"COSPSNPOSNCCIOFIPFNF"),
        // ("mix_50", b"PSNPPSICSFFFIPISONNNPIFCOOINCFSOOFICCCCSOCNONSFPIP"),
        // ("mix_100", b"NSOOPNFNSIINOPCNNONCFOPPPNFISSFSPCFSIPNOCSISIFCFPOPCOCSCSICFFICFNPPOSNONSSFOCPCOSNOCCPNPFIIOIFFNCIII"),
    ]
}

// Tree corpus, bare atoms
fn tree_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("c_1", b"(C)"),
        ("c_5_1", b"(CCCCC)"),
        ("c_5_2", b"C(C(C(C(C))))"),
        ("c_5_3", b"(((((C)C)C)C)C)"),
        ("c_5_4", b"C(CCC)C"),
        ("c_5_5", b"CC(C)CC"),
        ("c_5_6", b"CC(C)(C)C"),
        ("c_10_1", b"(CCCCCCCCCC)"),
        ("c_10_2", b"C(C(C(C(C(C(C(C(C(C)))))))))"),
        ("c_10_3", b"((((((((((C)C)C)C)C)C)C)C)C)C)"),
        ("c_10_4", b"C(CCCCCCCC)C"),
        ("c_10_5", b"CCC(C)C(CC)CCC"),
        ("c_10_6", b"CC(C)(C)C(C)C(C)C"),
        ("c_10_7", b"CCC(C(C)C)CC(C)C"),
        ("c_10_8", b"C(C)C(C)C(C)C(C)C(C)CCCCC"),
        // ("c_50_r1", b"C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_50_r5", b"C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_50_r10", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_50_r25", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_100_r1", b"C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_100_r5", b"C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_100_r10", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_100_r25", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_100_r50", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_1000_r1", b"C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_1000_r5", b"C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_1000_r10", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_1000_r25", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_1000_r50", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_1000_r100", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
        // ("c_1000_r250", b"C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"),
    ]
}

// Chain bonds corpus, bare atoms
fn chain_bonds_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("c_2", b"C-C"),
        ("c_2_stereo", b"C/C"),
        ("c_2_double", b"C=C"),
        ("c_2_triple", b"C#C"),
        ("c_2_quadruple", b"C$C"),
        ("c_5", b"C=C=C=C=C"),
        ("c_5_stereo", b"C/C=C\\C=C/C"),
        ("c_10", b"C=C=C=C=C=C=C=C=C=C"),
        ("c_10_aromatic", b"C:C:C:C:C:C:C:C:C:C"),
        ("c_10_stereo", b"C/C=C\\C=C/C=C\\C=C/C"),
        ("c_10_mixed", b"C-C=C#C$C:C/C\\C-C=C"),
        ("mixed_10", b"C-N=C-Br-N=F-S-O-Cl=N"),
        // ("c_50", b"C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C=C"),
        // ("c_50_mixed", b"C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C"),
        // ("c_100", b"C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C-C#C"),
        // ("c_100_mixed", b"C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C-C=C#C$C:C/C\\C-C=C"),
        // ("c_1000", b"C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C"),
    ]
}

// Tree bonds corpus, bare atoms
fn tree_bonds_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("c_5_1", b"(C=C=C=C=C)"),
        ("c_5_2", b"C(-C(-C(-C(-C))))"),
        ("c_5_3", b"(((((C)-C)-C)-C)-C)"),
        ("c_10_1", b"(C/C=C\\C=C/C=C\\C=C/C)"),
        ("c_10_2", b"C(-C(#C(-C(#C(-C(#C(-C(#C(-C)))))))))"),
        ("c_10_3", b"((((((((((C)=C)=C)=C)=C)=C)=C)=C)=C)=C)"),
        ("c_10_4", b"C(-C=C-C=C-C=C-C=C)=C"),
        ("c_10_5", b"C=CC(C)=C(C=C)C=CC"),
        ("c_10_6", b"CC(C)(C)C(=C)C(=C)C"),
        ("c_10_7", b"CCC(-C(-C)-C)C=C(-C)C"),
        ("c_10_8", b"C(C)=C(C)C(C)=C(C)C(C)=CC=CC=C"),
        ("c_10_9", b"(C=C=C=C=C=C=C=C=C=C)"),
        // ("c_50", b"C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)"),
        // ("c_100", b"C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)C(=C)"),
    ]
}

// Ring corpus, simple cycles
fn ring_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("c3", b"C1CC1"),
        ("c4", b"C1CCC1"),
        ("c5", b"C1CCCC1"),
        ("c6", b"C1CCCCC1"),
        ("c8", b"C1CCCCCCC1"),
        ("pct12", b"C%12CCCC%12"),
    ]
}

// Ring corpus, fused and spiro examples
fn ring_fused_spiro_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("fused_1", b"C1CC2CCCCC2CC1"),
        ("fused_2", b"C1CCC2CCCC2CC1"),
        ("spiro_1", b"C1CCC2(CC1)CCC2"),
        ("spiro_2", b"C1CC2(C1)CC2C"),
    ]
}

// Ring corpus, directed closures and percent indices
fn ring_stereo_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("dir_up_open", b"C/1CC1"),
        // ("dir_up_close", b"C1CC/1"),
        // ("dir_up_both", b"C/1CC/1"),
        // ("dir_down_both", b"C\\1CC\\1"),
        // ("pct_dir_up_open", b"C/%12CC%12"),
        // ("pct_dir_up_close", b"C%12CC/%12"),
        // ("pct_dir_down_both", b"C\\%12CC\\%12"),
    ]
}
// Ring corpus, many rings
fn ring_complex_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("adamantane", b"C1C2CC3CC1CC(C2)C3"),
        // ("dodecahedrane", b"C12C3C4C5C1C1C6C2C2C3C3C4C4C5C1C1C6C2C3C41"),
        // ("closo-dodecaborane", b"[H]B1234B567([H])B189([H])B21%10([H])B32%11([H])B453([H])B645([H])B786([H])B478([H])B12([H])([B-]9%1067[H])[B-]3%1158[H]"),
        // ("c60_fullerene", b"C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23"),
        // ("vitamin_b12", b"[H][C@]12[C@H](CC(N)=O)[C@@]3(C)CCC(=O)NC[C@@H](C)OP(=O)([O-])O[C@H]4[C@@H](O)[C@H](O[C@@H]4COP(=O)(O)O)n4c[n+](c5cc(C)c(C)cc54)[Co-3]456([CH2][C@H]7O[C@@H](n8cnc9c(N)ncnc98)[C@H](O)[C@@H]7O)[N]1C3=C(C)C1=[N+]4C(=CC3=[N+]5C(=C(C)C4=[N+]6[C@]2(C)[C@@](C)(CC(N)=O)[C@@H]4CCC(N)=O)[C@@](C)(CC(N)=O)[C@@H]3CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O"),
    ]
}

// Components corpus
fn component_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("two", b"CC.CC"),
        // ("three", b"C.C.C"),
        // ("grouped", b"C(C).C(C)"),
        // ("rings_across_digit", b"C1.CC1"),
        // ("rings_across_pct", b"C%12.CC%12"),
        // ("stereo_up", b"C/1.CC/1"),
        // ("stereo_down", b"C\\1.CC\\1"),
        // ("stereo_up_pct", b"C/%12.CC/%12"),
        // ("stereo_down_pct", b"C\\%12.CC\\%12"),
    ]
}

// Bracket corpus
fn bracket_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("brkt_C_50", b"[C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C][C]"),
        // ("brkt_c_50", b"[c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c][c]"),
        // ("brkt_chiral_50", b"[C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H][C@H]"),
        // ("brkt_charge_50", b"[C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1][C+1]"),
        // ("brkt_class_50", b"[C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1][C:1]"),
        // ("brkt_hcount_50", b"[CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3][CH3]"),
        // ("brkt_mixed_100", b"[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C[C]C"),
    ]
}

// Whitespace corpus
fn whitespace_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        ("empty_space", b" "),
        // ("empty_tabs", b"\t\t"),
        // ("empty_newlines", b"\n\n"),
        // ("empty_crlf", b"\r\n"),
        ("trail_space_short", b"CC "),
        ("trail_space_long", b"CCCCCCCCCC     "),
        // ("trail_tab", b"CC\t"),
        // ("trail_tabs", b"CC\t\t\t"),
        // ("trail_newline", b"CC\n"),
        // ("trail_newlines", b"CC\n\n"),
        // ("trail_cr", b"CC\r"),
        // ("trail_crlf", b"CC\r\n"),
        // ("trail_mixed", b"CC \t\r\n\t "),
    ]
}
// Wildcard corpus
fn wildcard_inputs() -> Vec<(&'static str, &'static [u8])> {
    vec![
        (
            "star_50",
            b"**************************************************",
        ),
        // ("star_200", b"********************************************************************************************************************************************************************************************************"),
        // ("brkt_star_50", b"[*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*][*]"),
        // ("star_in_chain_100", b"C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*C*"),
        // ("star_in_chain_bonds_100", b"C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*C-*"),
    ]
}

fn smiles_parsing(c: &mut Criterion) {
    // Chains
    let mut group_chain = c.benchmark_group("smiles_parsing/chain");
    for (name, s) in chain_inputs().iter() {
        group_chain.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_smiles_bytes_to_table_ir(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_chain.finish();

    // Branches
    let mut group_tree = c.benchmark_group("smiles_parsing/tree");
    for (name, s) in tree_inputs().iter() {
        group_tree.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_smiles_bytes_to_table_ir(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_tree.finish();

    // Chain bonds
    let mut group_chain_bonds = c.benchmark_group("smiles_parsing/chain_bonds");
    for (name, s) in chain_bonds_inputs().iter() {
        group_chain_bonds.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_smiles_bytes_to_table_ir(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_chain_bonds.finish();

    // Tree bonds
    let mut group_tree_bonds = c.benchmark_group("smiles_parsing/tree_bonds");
    for (name, s) in tree_bonds_inputs().iter() {
        group_tree_bonds.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_smiles_bytes_to_table_ir(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_tree_bonds.finish();

    // Rings
    let mut group_rings = c.benchmark_group("smiles_parsing/rings");
    for (name, s) in ring_inputs().iter() {
        group_rings.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_smiles_bytes_to_table_ir(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_rings.finish();

    // Fused & spiro rings
    let mut group_fused = c.benchmark_group("smiles_parsing/fused_spiro");
    for (name, s) in ring_fused_spiro_inputs().iter() {
        group_fused.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_bytes_to_table_ir(input);
                assert!(result.is_ok());
            })
        });
    }
    group_fused.finish();

    // Ring stereo
    let mut group_ring_stereo = c.benchmark_group("smiles_parsing/ring_stereo");
    for (name, s) in ring_stereo_inputs().iter() {
        group_ring_stereo.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_bytes_to_table_ir(input);
                assert!(result.is_ok());
            })
        });
    }
    group_ring_stereo.finish();

    // Complex rings
    let mut group_complex_rings = c.benchmark_group("smiles_parsing/complex_rings");
    for (name, s) in ring_complex_inputs().iter() {
        group_complex_rings.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_bytes_to_table_ir(input);
                assert!(result.is_ok());
            })
        });
    }
    group_complex_rings.finish();

    // Components
    let mut group_components = c.benchmark_group("smiles_parsing/components");
    for (name, s) in component_inputs().iter() {
        group_components.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_bytes_to_table_ir(input);
                assert!(result.is_ok());
            })
        });
    }
    group_components.finish();

    // Brackets
    let mut group_brackets = c.benchmark_group("smiles_parsing/brackets");
    for (name, bytes) in bracket_inputs().iter() {
        group_brackets.bench_with_input(BenchmarkId::from_parameter(name), bytes, |b, input| {
            b.iter(|| {
                let input = black_box(input);
                let result = parse_smiles_bytes_to_table_ir(input);
                assert!(result.is_ok());
            })
        });
    }
    group_brackets.finish();

    // Whitespace
    let mut group_whitespace = c.benchmark_group("smiles_parsing/whitespace");
    for (name, s) in whitespace_inputs().iter() {
        group_whitespace.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_smiles_bytes_to_table_ir(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_whitespace.finish();
}

fn extended_smiles_parsing(c: &mut Criterion) {
    // Chains
    let mut group_chain = c.benchmark_group("extended_smiles_parsing/chain");
    for (name, s) in chain_inputs().iter() {
        group_chain.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_chain.finish();

    // Branches
    let mut group_tree = c.benchmark_group("extended_smiles_parsing/tree");
    for (name, s) in tree_inputs().iter() {
        group_tree.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_tree.finish();

    // Chain bonds
    let mut group_chain_bonds = c.benchmark_group("extended_smiles_parsing/chain_bonds");
    for (name, s) in chain_bonds_inputs().iter() {
        group_chain_bonds.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_chain_bonds.finish();

    // Tree bonds
    let mut group_tree_bonds = c.benchmark_group("extended_smiles_parsing/tree_bonds");
    for (name, s) in tree_bonds_inputs().iter() {
        group_tree_bonds.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_tree_bonds.finish();

    // Rings
    let mut group_rings = c.benchmark_group("extended_smiles_parsing/rings");
    for (name, s) in ring_inputs().iter() {
        group_rings.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_rings.finish();

    // Fused & spiro rings
    let mut group_fused = c.benchmark_group("extended_smiles_parsing/fused_spiro");
    for (name, s) in ring_fused_spiro_inputs().iter() {
        group_fused.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_fused.finish();

    // Ring stereo
    let mut group_ring_stereo = c.benchmark_group("extended_smiles_parsing/ring_stereo");
    for (name, s) in ring_stereo_inputs().iter() {
        group_ring_stereo.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_ring_stereo.finish();

    // Complex rings
    let mut group_complex_rings = c.benchmark_group("extended_smiles_parsing/complex_rings");
    for (name, s) in ring_complex_inputs().iter() {
        group_complex_rings.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_complex_rings.finish();

    // Components
    let mut group_components = c.benchmark_group("extended_smiles_parsing/components");
    for (name, s) in component_inputs().iter() {
        group_components.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_components.finish();

    // Brackets
    let mut group_brackets = c.benchmark_group("extended_smiles_parsing/brackets");
    for (name, bytes) in bracket_inputs().iter() {
        group_brackets.bench_with_input(BenchmarkId::from_parameter(name), bytes, |b, input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_brackets.finish();

    // Whitespace
    let mut group_whitespace = c.benchmark_group("extended_smiles_parsing/whitespace");
    for (name, s) in whitespace_inputs().iter() {
        group_whitespace.bench_with_input(BenchmarkId::from_parameter(name), s, |b, &input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_whitespace.finish();

    // Wildcards
    let mut group_wild = c.benchmark_group("extended_smiles_parsing/wildcards");
    for (name, bytes) in wildcard_inputs().iter() {
        group_wild.bench_with_input(BenchmarkId::from_parameter(name), bytes, |b, input| {
            b.iter(|| {
                let result = parse_extended_smiles_bytes(black_box(input));
                assert!(result.is_ok());
            })
        });
    }
    group_wild.finish();
}

criterion_group!(benches, smiles_parsing, extended_smiles_parsing);
criterion_main!(benches);
