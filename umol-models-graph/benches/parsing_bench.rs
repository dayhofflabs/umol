//! Benchmark parsing of MOL atom.inputs

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use nom::Parser;
use umol_models_graph::io::ctab::{
    atom::atom_input,
    bond::bond_input,
    counts::counts_input,
    properties::property_input,
};

fn parsing_benchmarks(c: &mut Criterion) {
    {
        let mut group = c.benchmark_group("mol_parsing/counts");
        let test_cases = [
            ("valid", &b"  6  5  0  0  1  0  0  0  0  0999 V2000"[..]),
            (
                "invalid",
                &b"  4  2  0     0  0            999 V1000"[..],
            ),
        ];
        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| counts_input().parse(std::hint::black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/atom");

        let test_cases = [
            (
                "len69",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0  0  0  0  0  0  0"[..],
            ),
            (
                "len69_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0  0  0  X  0  0  0"[..],
            ),
            (
                "len51",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0"[..],
            ),
            (
                "len51_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  X"[..],
            ),
            (
                "len42",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0"[..],
            ),
            (
                "len42_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  X"[..],
            ),
            (
                "len39",
                &b"   -0.1234    0.4560    0.7890 C   0  0"[..],
            ),
            (
                "len39_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  X"[..],
            ),
            (
                "len36",
                &b"   -0.1234    0.4560    0.7890 C   0"[..],
            ),
            (
                "len36_invalid",
                &b"   -0.1234    0.4560    0.7890 C   X"[..],
            ),
            (
                "len34",
                &b"   -0.1234    0.4560    0.7890 C  "[..],
            ),
            (
                "len34_invalid",
                &b"   -0.1234    0.4560    0.7890 X  "[..],
            ),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| atom_input().parse(std::hint::black_box(input)))
        });
        }

        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/bond");
        let test_cases = [
            ("len21", &b"  1  2  1  1  0  0  0"[..]),
            ("len21_invalid", &b"  1  2  1  A  0  0  0"[..]),
            ("len12", &b"  1  3  1  1"[..]),
            ("len12_invalid", &b"  1  2  1  A"[..]),
            ("len9", &b"  1  2  1"[..]),
            ("len9_invalid", &b"  1  2  A"[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| bond_input().parse(std::hint::black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/properties");
        let test_cases = [
            ("chg1", &b"M  CHG  1   1  -1"[..]),
            (
                "chg6",
                &b"M  CHG  6   1  -1   2   1   3  -1   4   1   5  -1   6   1"[..],
            ),
            ("chg_invalid", &b"M  CHG  1   1  -1  a"[..]),
            ("rad1", &b"M  RAD  1   1   2"[..]),
            (
                "rad6",
                &b"M  RAD  6   1   1   2   2   3   3   4   1   5   2   6   3"[..],
            ),
            ("rad_invalid", &b"M  RAD  1   1   4"[..]),
            ("iso1", &b"M  ISO  1   1  13"[..]),
            (
                "iso6",
                &b"M  ISO  6   1  13   2  14   3  12   4  13   5  14   6  12"[..],
            ),
            ("iso_invalid", &b"M  ISO  1   1 130"[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| property_input(std::hint::black_box(input)))
            });
        }
        group.finish();
    }
}

criterion_group!(benches, parsing_benchmarks);
criterion_main!(benches);
