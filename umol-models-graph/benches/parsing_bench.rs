//! Benchmark for of MOL format parsing

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use nom::Parser;
use umol_models_graph::io::config::ParseFlags;
use umol_models_graph::io::ctab::parser::{
    atom_input, atomlike_input, basic_property_input, bond_input, bondlike_input, counts_input,
    legacy_atom_list_input, property_input,
};

fn parsing_benchmarks(c: &mut Criterion) {
    {
        let mut group = c.benchmark_group("mol_parsing/counts");
        let test_cases = [
            ("valid", &b"  6  5  0  0  1  0  0  0  0  0999 V2000"[..]),
            ("invalid", &b"  4  2  0     0  0            999 V1000"[..]),
        ];
        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| counts_input(ParseFlags::BASIC).parse(black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/atom");

        let test_cases = [
            (
                "len69_basic",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0  0  0  0  0  0  0"[..],
            ),
            (
                "len69_basic_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0  0  0  0  0  0  X"[..],
            ),
            (
                "len51_basic",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  0"[..],
            ),
            (
                "len51_basic_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0  0  0  X"[..],
            ),
            (
                "len42_basic",
                &b"   -0.1234    0.4560    0.7890 C   0  0  0"[..],
            ),
            (
                "len42_basic_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  X"[..],
            ),
            (
                "len39_basic",
                &b"   -0.1234    0.4560    0.7890 C   0  0"[..],
            ),
            (
                "len39_basic_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  X"[..],
            ),
            ("len36_basic", &b"   -0.1234    0.4560    0.7890 C   0"[..]),
            (
                "len36_basic_invalid",
                &b"   -0.1234    0.4560    0.7890 C   X"[..],
            ),
            ("len34_basic", &b"   -0.1234    0.4560    0.7890 C  "[..]),
            (
                "len34_basic_invalid",
                &b"   -0.1234    0.4560    0.7890 X  "[..],
            ),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| atom_input(ParseFlags::BASIC).parse(black_box(input)))
            });
        }

        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/atomlike");

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
            ("len42", &b"   -0.1234    0.4560    0.7890 C   0  0  0"[..]),
            (
                "len42_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  0  X"[..],
            ),
            ("len39", &b"   -0.1234    0.4560    0.7890 C   0  0"[..]),
            (
                "len39_invalid",
                &b"   -0.1234    0.4560    0.7890 C   0  X"[..],
            ),
            ("len36", &b"   -0.1234    0.4560    0.7890 C   0"[..]),
            (
                "len36_invalid",
                &b"   -0.1234    0.4560    0.7890 C   X"[..],
            ),
            ("len34", &b"   -0.1234    0.4560    0.7890 C  "[..]),
            ("len34_invalid", &b"   -0.1234    0.4560    0.7890 X  "[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| atomlike_input(ParseFlags::LENIENT).parse(std::hint::black_box(input)))
            });
        }

        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/bond");
        let test_cases = [
            ("len21_basic", &b"  1  2  1  1  0  0  0"[..]),
            ("len21_basic_invalid", &b"  1  2  1  A  0  0  0"[..]),
            ("len12_basic", &b"  1  3  1  1"[..]),
            ("len12_basic_invalid", &b"  1  2  1  A"[..]),
            ("len9_basic", &b"  1  2  1"[..]),
            ("len9_basic_invalid", &b"  1  2  A"[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| bond_input(ParseFlags::BASIC).parse(black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/bondlike");
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
                b.iter(|| bondlike_input(ParseFlags::LENIENT).parse(black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/legacy_atom_list");
        let test_cases = [
            ("no_exclusion", &b"  1 F    3   9   7   8  "[..]),
            ("exclusion", &b"  1 T    3   9   7   8  "[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| legacy_atom_list_input().parse(black_box(input)))
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("mol_parsing/basic_properties");
        let test_cases = [
            ("chg_basic", &b"M  CHG  1   1  -1"[..]),
            (
                "chg6_basic",
                &b"M  CHG  6   1  -1   2   1   3  -1   4   1   5  -1   6   1"[..],
            ),
            ("chg_basic_invalid", &b"M  CHG  1   1  -1  a"[..]),
            ("rad1_basic", &b"M  RAD  1   1   2"[..]),
            (
                "rad6_basic",
                &b"M  RAD  6   1   1   2   2   3   3   4   1   5   2   6   3"[..],
            ),
            ("rad_basic_invalid", &b"M  RAD  1   1   4"[..]),
            ("iso1_basic", &b"M  ISO  1   1  13"[..]),
            (
                "iso8_basic",
                &b"M  ISO  8   1  13   2  14   3  12   4  13   5  14   6  12   7  13   8  14"[..],
            ),
            ("iso_basic_invalid", &b"M  ISO  1   1  40"[..]),
            ("sty1_basic", &b"M  STY  1   1 SUP"[..]),
            ("sty2_basic", &b"M  STY  2   1 SUP   2 DAT"[..]),
            ("slb1_basic", &b"M  SLB  1   1   1"[..]),
            ("slb2_basic", &b"M  SLB  2   1  14   2  15"[..]),
            ("sal_simple_basic", &b"M  SAL  1  1   5"[..]),
            ("sal_multi_basic", &b"M  SAL  1  3   1   2   3"[..]),
            ("sbl_simple_basic", &b"M  SBL  1  1   3"[..]),
            ("sbl_multi_basic", &b"M  SBL  2  2   1   2"[..]),
            ("alias_simple_basic", &b"A    1 CF3"[..]),
            ("alias_long_basic", &b"A   15 Ph"[..]),
            ("value_simple_basic", &b"V    1 *"[..]),
            ("value_long_basic", &b"V   15 query"[..]),
            ("sst1_basic", &b"M  SST  1   1 ALT"[..]),
            ("sst2_basic", &b"M  SST  2   1 RAN   2 BLO"[..]),
            ("smt_simple_basic", &b"M  SMT   1 n"[..]),
            ("smt_long_basic", &b"M  SMT   2 CH2CH2"[..]),
            ("zbo1_basic", &b"M  ZBO  1   1   0"[..]),
            (
                "zbo4_basic",
                &b"M  ZBO  4   1   0   2   0   3   0   4   0"[..],
            ),
            ("zch1_basic", &b"M  ZCH  1   1  -1"[..]),
            (
                "zch4_basic",
                &b"M  ZCH  4   1  -1   2   1   3  -1   4   1"[..],
            ),
            ("hyd1_basic", &b"M  HYD  1   1   3"[..]),
            (
                "hyd6_basic",
                &b"M  HYD  6   1   3   2   2   3   1   4   0   5   3   6   2"[..],
            ),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| basic_property_input(ParseFlags::BASIC).parse(black_box(input)))
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
            (
                "iso8",
                &b"M  ISO  8   1  13   2  14   3  12   4  13   5  14   6  12   7  13   8  14"[..],
            ),
            ("sty1", &b"M  STY  1   1 SUP"[..]),
            (
                "sty8",
                &b"M  STY  8   1 SUP   2 DAT   3 MUL   4 SRU   5 GEN   6 SUP   7 DAT   8 MUL"[..],
            ),
            ("sst1", &b"M  SST  1   1 ALT"[..]),
            ("smt_multiplier", &b"M  SMT   1 n"[..]),
            ("smt_subscript", &b"M  SMT   2 CH2CH2"[..]),
            ("slb1", &b"M  SLB  1   1 Et "[..]),
            ("sal5", &b"M  SAL  3  5   1   2   3   4   5"[..]),
            ("sbl4", &b"M  SBL  2  4   1   2   3   4"[..]),
            ("alias", &b"A    1 CF3"[..]),
            ("value", &b"V    1 *"[..]),
            ("als1", &b"M  ALS  1  3FC   N   O   "[..]),
            ("als4", &b"M  ALS  5  4TF   C   N   O   S"[..]),
            ("apo1", &b"M  APO  1   1   1"[..]),
            ("aal1", &b"M  AAL  1 1   2   1"[..]),
            ("aal2", &b"M  AAL  3 2   4   1   5   2"[..]),
            ("rbc1", &b"M  RBC  1   1   2"[..]),
            ("sub1", &b"M  SUB  1   1   3"[..]),
            ("uns1", &b"M  UNS  1   1   1"[..]),
            ("lin1", &b"M  LIN  1   1   2   5   7"[..]),
            ("zbo1", &b"M  ZBO  1   1   0"[..]),
            ("zbo4", &b"M  ZBO  4   1   0   2   0   3   0   4   0"[..]),
            ("zch1", &b"M  ZCH  1   1  -1"[..]),
            ("zch4", &b"M  ZCH  4   1  -1   2   1   3  -1   4   1"[..]),
            ("hyd1", &b"M  HYD  1   1   3"[..]),
            (
                "hyd6",
                &b"M  HYD  6   1   3   2   2   3   1   4   0   5   3   6   2"[..],
            ),
            ("scn1", &b"M  SCN  1   1 HH"[..]),
            ("scn2", &b"M  SCN  2   1 HT   2 HH"[..]),
            ("sds1", &b"M  SDS EXP  1   1"[..]),
            (
                "spa12",
                &b"M  SPA  1 12   3   4   5   6   9  10  11  12  13  14  15  16"[..],
            ),
            ("crs", &b"M  CRS   1  3  10   9   4"[..]),
            ("sdi1", &b"M  SDI  1   1    1.2    2.3"[..]),
            (
                "sdi2",
                &b"M  SDI  2   1    1.2    2.3   2    3.4    4.5"[..],
            ),
            ("sbv", &b"M  SBV  1   1    1.0    2.0"[..]),
            ("sdt", &b"  SDT   1 pH   "[..]),
            (
                "sdd",
                &b"M  SDD   1     0.0000    0.0000    DR    ALL  1       6"[..],
            ),
            ("scd", &b"M  SCD   1 4.6"[..]),
            ("sed", &b"M  SED   2 E/Z unknown"[..]),
            ("spl1", &b"M  SPL  1   1   2"[..]),
            ("snc1", &b"M  SNC  1   1   5"[..]),
            ("rgp1", &b"M  RGP  1   1   1"[..]),
            ("log1", &b"M  LOG  1   1   0   0  >2"[..]),
        ];

        for (id, data) in test_cases.iter() {
            group.bench_with_input(BenchmarkId::from_parameter(id), data, |b, &input| {
                b.iter(|| property_input(ParseFlags::LENIENT).parse(black_box(input)))
            });
        }
        group.finish();
    }
}

criterion_group!(benches, parsing_benchmarks);
criterion_main!(benches);
