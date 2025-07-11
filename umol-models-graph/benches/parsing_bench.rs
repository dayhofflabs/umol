//! Benchmark for of MOL format parsing

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
        let mut group = c.benchmark_group("mol_parsing/atom_general");

        let test_cases = [
            (
                "len69_general",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0  0  0  0  0  0  0"[..],
            ),
            (
                "len69_general_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0  0  0  0  0  0  X"[..],
            ),
            (
                "len51_general",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0"[..],
            ),
            (
                "len51_general_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  X"[..],
            ),
            (
                "len42_general",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0"[..],
            ),
            (
                "len42_general_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  X"[..],
            ),
            (
                "len39_general",
                &b"   -0.1234    0.4560    0.7890 C   0  0"[..],
            ),
            (
                "len39_general_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  X"[..],
            ),
            (
                "len36_general",
                &b"   -0.1234    0.4560    0.7890 C   0"[..],
            ),
            (
                "len36_general_invalid",
                &b"   -0.1234    0.4560    0.7890 C   X"[..],
            ),
            (
                "len34_general",
                &b"   -0.1234    0.4560    0.7890 C  "[..],
            ),
            (
                "len34_general_invalid",
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
        let mut group = c.benchmark_group("mol_parsing/bond_general");
        let test_cases = [
            ("len21_general", &b"  1  2  1  0     2  1"[..]),
            ("len21_general_invalid", &b"  1  2  1  0     4  1"[..]),
            ("len18_ring_general", &b"  1  2  8  0     1"[..]),
            ("len18_general_invalid", &b"  1  2  8  0     X"[..]),
            ("len12_general", &b"  1  3  2  3"[..]),
            ("len12_general_invalid", &b"  1  2  1  A"[..]),
            ("len9_general", &b"  1  2  1"[..]),
            ("len9_general_invalid", &b"  1  2  9"[..]),
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
            ("chg1_standard", &b"M  CHG  1   1  -1"[..]),
            (
                "chg6_standard",
                &b"M  CHG  6   1  -1   2   1   3  -1   4   1   5  -1   6   1"[..],
            ),
            ("chg_standard_invalid", &b"M  CHG  1   1  -1  a"[..]),
            ("rad1_standard", &b"M  RAD  1   1   2"[..]),
            (
                "rad6_standard",
                &b"M  RAD  6   1   1   2   2   3   3   4   1   5   2   6   3"[..],
            ),
            ("rad_standard_invalid", &b"M  RAD  1   1   4"[..]),
            ("iso1_standard", &b"M  ISO  1   1  13"[..]),
            ("iso8_standard", &b"M  ISO  8   1  13   2  14   3  12   4  13   5  14   6  12   7  13   8  14"[..]),
            ("iso_standard_invalid", &b"M  ISO  1   1  40"[..]),
            ("sty1_standard", &b"M  STY  1   1 SUP"[..]),
            ("sty2_standard", &b"M  STY  2   1 SUP   2 DAT"[..]),
            ("slb1_standard", &b"M  SLB  1   1   1"[..]),
            ("slb2_standard", &b"M  SLB  2   1  14   2  15"[..]),
            ("sal_simple_standard", &b"M  SAL  1  1   5"[..]),
            ("sal_multi_standard", &b"M  SAL  1  3   1   2   3"[..]),
            ("sbl_simple_standard", &b"M  SBL  1  1   3"[..]),
            ("sbl_multi_standard", &b"M  SBL  2  2   1   2"[..]),
            ("alias_simple_standard", &b"A    1 CF3"[..]),
            ("alias_long_standard", &b"A   15 Ph"[..]),
            ("value_simple_standard", &b"V    1 *"[..]),
            ("value_long_standard", &b"V   15 query"[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| property_input_standard(std::hint::black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/properties_general");
        let test_cases = [
            ("chg1_general", &b"M  CHG  1   1  -1"[..]),
            (
                "chg8_general",
                &b"M  CHG  8   1  -1   2   1   3  -1   4   1   5  -1   6   1   7  -1   8   1"[..],
            ),
            ("rad1_general", &b"M  RAD  1   1   2"[..]),
            ("sty1_general", &b"M  STY  1   1 SUP"[..]),
            ("sty8_general", &b"M  STY  8   1 SUP   2 DAT   3 MUL   4 SRU   5 GEN   6 SUP   7 DAT   8 MUL"[..]),
            ("slb1_general", &b"M  SLB  1   1 Et "[..]),
            ("sal_complex_general", &b"M  SAL  3  5   1   2   3   4   5"[..]),
            ("sbl_complex_general", &b"M  SBL  2  4   1   2   3   4"[..]),
            ("alias_general", &b"A    1 CF3"[..]),
            ("value_general", &b"V    1 *"[..]),
            ("als1_general", &b"M  ALS  1  3FC   N   O   "[..]),
            ("als_complex_general", &b"M  ALS  5  4TF   C   N   O   S"[..]),
            ("apo1_general", &b"M  APO  1   1   1"[..]),
            ("apo4_general", &b"M  APO  4   1   1   2   2   3   3   4   1"[..]),
            ("aal1_general", &b"M  AAL  1 1   2   1"[..]),
            ("aal2_general", &b"M  AAL  3 2   4   1   5   2"[..]),
            ("rbc1_general", &b"M  RBC  1   1   2"[..]),
            ("rbc8_general", &b"M  RBC  8   1  -2   2  -1   3   0   4   2   5   3   6   4   7   5   8   6"[..]),
            ("sub1_general", &b"M  SUB  1   1   3"[..]),
            ("sub8_general", &b"M  SUB  8   1  -2   2  -1   3   0   4   1   5   2   6   3   7   4   8   5"[..]),
            ("uns1_general", &b"M  UNS  1   1   1"[..]),
            ("uns4_general", &b"M  UNS  4   1   1   2   1   3   0   4   1"[..]),
            ("lin1_general", &b"M  LIN  1   1   2   5   7"[..]),
            ("lin2_general", &b"M  LIN  2   1   2   5   7   3   3   0   0"[..]),
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
