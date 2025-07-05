//! Benchmark parsing of MOL atom.inputs

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use nom::Parser;
use umol_models_graph::io::ctab::parser::{
    atom::{atom_input, atom_input_standard},
    bond::{bond_input, bond_input_standard},
    counts::counts_input,
    properties::{property_input, property_input_standard},
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
        let mut group = c.benchmark_group("mol_parsing/atom_standard");

        let test_cases = [
            (
                "len69_standard",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0  0  0  0  0  0  0"[..],
            ),
            (
                "len69_standard_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0  0  0  0  0  0  X"[..],
            ),
            (
                "len51_standard",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0"[..],
            ),
            (
                "len51_standard_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  X"[..],
            ),
            (
                "len42_standard",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0"[..],
            ),
            (
                "len42_standard_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  X"[..],
            ),
            (
                "len39_standard",
                &b"   -0.1234    0.4560    0.7890 C   0  0"[..],
            ),
            (
                "len39_standard_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  X"[..],
            ),
            (
                "len36_standard",
                &b"   -0.1234    0.4560    0.7890 C   0"[..],
            ),
            (
                "len36_standard_invalid",
                &b"   -0.1234    0.4560    0.7890 C   X"[..],
            ),
            (
                "len34_standard",
                &b"   -0.1234    0.4560    0.7890 C  "[..],
            ),
            (
                "len34_standard_invalid",
                &b"   -0.1234    0.4560    0.7890 X  "[..],
            ),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| atom_input_standard().parse(std::hint::black_box(input)))
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
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0  0  0  0  0  0  X"[..],
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
        let mut group = c.benchmark_group("mol_parsing/bond_standard");
        let test_cases = [
            ("len21_standard", &b"  1  2  1  1  0  0  0"[..]),
            ("len21_standard_invalid", &b"  1  2  1  A  0  0  0"[..]),
            ("len12_standard", &b"  1  3  1  1"[..]),
            ("len12_standard_invalid", &b"  1  2  1  A"[..]),
            ("len9_standard", &b"  1  2  1"[..]),
            ("len9_standard_invalid", &b"  1  2  A"[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| bond_input_standard().parse(std::hint::black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/bond");
        let test_cases = [
            ("len21", &b"  1  2  1  0     2  1"[..]),
            ("len21_invalid", &b"  1  2  1  0     4  1"[..]),
            ("len18_ring", &b"  1  2  8  0     1"[..]),
            ("len18_invalid", &b"  1  2  8  0     X"[..]),
            ("len12", &b"  1  3  2  3"[..]),
            ("len12_invalid", &b"  1  2  1  A"[..]),
            ("len9", &b"  1  2  1"[..]),
            ("len9_invalid", &b"  1  2  9"[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| bond_input().parse(std::hint::black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/properties_standard");
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
            ("sty1", &b"M  STY  1   1 SUP"[..]),
            ("sty2", &b"M  STY  2   1 SUP   2 DAT"[..]),
            ("slb1", &b"M  SLB  1   1   1"[..]),
            ("slb2", &b"M  SLB  2   1  14   2  15"[..]),
            ("sal_simple", &b"M  SAL  1  1   5"[..]),
            ("sal_multi", &b"M  SAL  1  3   1   2   3"[..]),
            ("sbl_simple", &b"M  SBL  1  1   3"[..]),
            ("sbl_multi", &b"M  SBL  2  2   1   2"[..]),
            ("alias_simple", &b"A    1 CF3"[..]),
            ("alias_long", &b"A   15 Ph"[..]),
            ("value_simple", &b"V    1 *"[..]),
            ("value_long", &b"V   15 query"[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| property_input_standard(std::hint::black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/properties");
        let test_cases = [
            ("chg1", &b"M  CHG  1   1  -1"[..]),
            (
                "chg8",
                &b"M  CHG  8   1  -1   2   1   3  -1   4   1   5  -1   6   1   7  -1   8   1"[..],
            ),
            ("rad1", &b"M  RAD  1   1   2"[..]),
            ("iso1", &b"M  ISO  1   1  13"[..]),
            ("sty1", &b"M  STY  1   1 SUP"[..]),
            ("sty8", &b"M  STY  8   1 SUP   2 DAT   3 MUL   4 SRU   5 GEN   6 SUP   7 DAT   8 MUL"[..]),
            ("slb1", &b"M  SLB  1   1 Et "[..]),
            ("sal_complex", &b"M  SAL  3  5   1   2   3   4   5"[..]),
            ("sbl_complex", &b"M  SBL  2  4   1   2   3   4"[..]),
            ("alias", &b"A    1 CF3"[..]),
            ("value", &b"V    1 *"[..]),
            ("als1", &b"M  ALS  1  3FC   N   O   "[..]),
            ("als_complex", &b"M  ALS  5  4TF   C   N   O   S"[..]),
            ("apo1", &b"M  APO  1   1   1"[..]),
            ("apo4", &b"M  APO  4   1   1   2   2   3   3   4   1"[..]),
            ("aal1", &b"M  AAL  1 1   2   1"[..]),
            ("aal2", &b"M  AAL  3 2   4   1   5   2"[..]),
            ("rbc1", &b"M  RBC  1   1   2"[..]),
            ("rbc8", &b"M  RBC  8   1  -2   2  -1   3   0   4   2   5   3   6   4   7   5   8   6"[..]),
            ("sub1", &b"M  SUB  1   1   3"[..]),
            ("sub8", &b"M  SUB  8   1  -2   2  -1   3   0   4   1   5   2   6   3   7   4   8   5"[..]),
            ("uns1", &b"M  UNS  1   1   1"[..]),
            ("uns4", &b"M  UNS  4   1   1   2   1   3   0   4   1"[..]),
            ("lin1", &b"M  LIN  1   1   2   5   7"[..]),
            ("lin2", &b"M  LIN  2   1   2   5   7   3   3   0   0"[..]),
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
