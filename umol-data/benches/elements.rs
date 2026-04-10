//! Benchmark for data structures

use criterion::{criterion_group, criterion_main, Criterion};
use umol_data::{e, Element, Isotope};

fn elements(c: &mut Criterion) {
    {
        let mut group = c.benchmark_group("element/parse_symbol");
        for symbol in &["C".to_string(), "Cu".to_string(), "Zn".to_string()] {
            group.bench_with_input(symbol, symbol, |b, symbol| {
                b.iter(|| {
                    let _ = Element::from_symbol(symbol);
                });
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("element/is_element");
        for symbol in &["C".to_string(), "Cu".to_string(), "Zn".to_string()] {
            group.bench_with_input(symbol, symbol, |b, symbol| {
                b.iter(|| {
                    let _ = Element::is_element(symbol);
                });
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("element/to_string");
        for element in &[e!(C), e!(Cu), e!(Zn)] {
            group.bench_with_input(element.symbol(), element, |b, element| {
                b.iter(|| {
                    let _ = element.to_string();
                });
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("element/next_element");
        for element in &[e!(C), e!(Ne), e!(Og)] {
            group.bench_with_input(element.symbol(), element, |b, element| {
                b.iter(|| {
                    let _ = element.next();
                });
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("isotope/from_symbol");
        for symbol in &["12C".to_string(), "63Cu".to_string(), "209Bi".to_string()] {
            group.bench_with_input(symbol, symbol, |b, symbol| {
                b.iter(|| {
                    let _ = Isotope::from_symbol(symbol);
                });
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("isotope/to_string");
        for (symbol, isotope) in &[
            ("12C".to_string(), Isotope::checked_new(e!(C), 12)),
            ("63Cu".to_string(), Isotope::checked_new(e!(Cu), 63)),
            ("209Bi".to_string(), Isotope::checked_new(e!(Bi), 209)),
        ] {
            group.bench_with_input(symbol, isotope, |b, isotope| {
                b.iter(|| {
                    let _ = isotope.map(|i| i.to_string());
                });
            });
        }
        group.finish();
    }

    {
        let mut group = c.benchmark_group("isotope/is_catalogued");
        for (symbol, inputs) in &[
            ("12C".to_string(), (e!(C), 12)),
            ("40C".to_string(), (e!(C), 40)),
            ("31P".to_string(), (e!(P), 31)),
            ("63Cu".to_string(), (e!(Cu), 63)),
            ("209Bi".to_string(), (e!(Bi), 209)),
        ] {
            group.bench_with_input(symbol, inputs, |b, &(element, mass_number)| {
                b.iter(|| {
                    let _ = Isotope::is_catalogued(element, mass_number);
                });
            });
        }
        group.finish();
    }
}

criterion_group!(benches, elements);
criterion_main!(benches);
