//! SVG rendering benchmarks over already constructed depictions.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use umol_chem::element::Element;
use umol_geometric_core::Point2D;
use umol_graph_ir::ir::{AtomForm, AtomId, BondForm, Molecule, MoleculeEntries};
use umol_io::depiction::molecule::depict;
use umol_io::depiction::Depiction;
use umol_io::layout::MoleculeLayout;
use umol_io::svg::render;

fn chain_depiction(atom_count: usize) -> Depiction {
    let atoms = vec![AtomForm::from_element(Element::C); atom_count];
    let bonds = (1..atom_count)
        .map(|index| {
            (
                AtomId::from(index - 1),
                AtomId::from(index),
                BondForm::from_order(if index % 3 == 0 { 2 } else { 1 }),
            )
        })
        .collect();
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    });
    let positions = (0..atom_count)
        .map(|index| Point2D::new(index as f64, f64::from((index % 2) as u8) * 0.5))
        .collect();
    let layout = MoleculeLayout::try_new(positions).unwrap();
    depict(&molecule, &layout).unwrap()
}

fn bench_svg_render(c: &mut Criterion) {
    let mut group = c.benchmark_group("svg/render");
    for atom_count in [8, 128] {
        let depiction = chain_depiction(atom_count);
        group.throughput(Throughput::Elements(depiction.items().len() as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(atom_count),
            &depiction,
            |b, depiction| b.iter(|| black_box(render(black_box(depiction)))),
        );
    }
    group.finish();
}

criterion_group!(svg, bench_svg_render);
criterion_main!(svg);
