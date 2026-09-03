//! SVG rendering benchmarks over already constructed depictions.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use umol_chem::element::Element;
use umol_geometric_core::Point2D;
use umol_graph_core::Correspondence;
use umol_graph_ir::ir::{AtomForm, AtomId, BondForm, Molecule, MoleculeEntries};
use umol_graph_ir::mol_dsl;
use umol_io::depiction::molecule::depict;
use umol_io::depiction::reaction::depict_from_sides;
use umol_io::depiction::Depiction;
use umol_io::layout::MoleculeLayout;
use umol_io::svg::render;

struct SvgCase {
    name: &'static str,
    depiction: Depiction,
}

impl SvgCase {
    fn new(name: &'static str, depiction: Depiction) -> Self {
        Self { name, depiction }
    }
}

fn molecule_depiction(molecule: &Molecule, positions: &[[f64; 2]]) -> Depiction {
    let layout =
        MoleculeLayout::try_new(positions.iter().map(|&[x, y]| Point2D::new(x, y)).collect())
            .expect("benchmark coordinates must be finite");
    depict(molecule, &layout).expect("benchmark layout must match the molecule frame")
}

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

fn svg_cases() -> Vec<SvgCase> {
    let labeled = mol_dsl!(r#"{:atoms ["C#i13#c+#h2" "O#c-"] :bonds [[0 1 "2"]]}"#);
    let tetrahedral = mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    );
    let fused_aromatic = mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]
                    [5 0 "1"] [2 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"]
                    [9 3 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5 6 7 8 9] :attrs "*"}]}"#
    );

    let lhs = mol_dsl!(
        r#"{:atoms ["C" "O" "N" "Cl"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]}"#
    );
    let rhs = mol_dsl!(
        r#"{:atoms ["C" "O" "N" "Cl" "F"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [0 4 "1"]]}"#
    );
    let lhs_layout = MoleculeLayout::try_new(vec![
        Point2D::new(-1.5, 0.0),
        Point2D::new(-0.5, 0.4),
        Point2D::new(0.5, -0.4),
        Point2D::new(1.5, 0.0),
    ])
    .expect("benchmark coordinates must be finite");
    let rhs_layout = MoleculeLayout::try_new(vec![
        Point2D::new(-1.0, 0.0),
        Point2D::new(0.0, 0.4),
        Point2D::new(1.0, 0.0),
        Point2D::new(2.0, 0.4),
        Point2D::new(-1.0, -1.0),
    ])
    .expect("benchmark coordinates must be finite");
    let correspondence = Correspondence::new(
        (0..4)
            .map(|index| (AtomId::from(index), AtomId::from(index)))
            .collect(),
        4,
        5,
    )
    .expect("benchmark atom correspondence must be a partial bijection");
    let reaction = depict_from_sides(&lhs, &lhs_layout, &rhs, &rhs_layout, &correspondence)
        .expect("benchmark frames must agree");

    vec![
        SvgCase::new(
            "labeled_atoms",
            molecule_depiction(&labeled, &[[0.0, 0.0], [1.0, 0.0]]),
        ),
        SvgCase::new(
            "tetrahedral_stereo",
            molecule_depiction(
                &tetrahedral,
                &[
                    [0.0, 0.0],
                    [1.0, 0.0],
                    [-0.5, 0.866],
                    [-0.5, -0.866],
                    [0.0, -1.25],
                ],
            ),
        ),
        SvgCase::new(
            "fused_aromatic",
            molecule_depiction(
                &fused_aromatic,
                &[
                    [-1.732, 0.5],
                    [-0.866, 1.0],
                    [0.0, 0.5],
                    [0.0, -0.5],
                    [-0.866, -1.0],
                    [-1.732, -0.5],
                    [0.866, 1.0],
                    [1.732, 0.5],
                    [1.732, -0.5],
                    [0.866, -1.0],
                ],
            ),
        ),
        SvgCase::new("mapped_reaction", reaction),
    ]
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

    let mut group = c.benchmark_group("svg/render/representative");
    for case in svg_cases() {
        group.throughput(Throughput::Elements(case.depiction.items().len() as u64));
        group.bench_with_input(case.name, &case.depiction, |b, depiction| {
            b.iter(|| black_box(render(black_box(depiction))))
        });
    }
    group.finish();
}

criterion_group!(svg, bench_svg_render);
criterion_main!(svg);
