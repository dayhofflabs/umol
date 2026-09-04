//! SVG rendering benchmarks over already constructed opaque depictions.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AtomForm, AtomId, BondDelta, BondFieldChange, BondForm, BondId, Delta, Deltas, Molecule,
    MoleculeEntries, NumForm, Reaction,
};
use umol_graph_ir::mol_dsl;
use umol_io::depict::{Depict, Depiction};

struct SvgCase {
    name: &'static str,
    depiction: Depiction,
}

impl SvgCase {
    fn new(name: &'static str, depiction: Depiction) -> Self {
        Self { name, depiction }
    }
}

fn molecule_depiction(molecule: &Molecule) -> Depiction {
    molecule
        .depict()
        .expect("benchmark molecule must be depictable")
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
    molecule_depiction(&molecule)
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
    let reaction = Reaction::new(
        mol_dsl!(
            r#"{:atoms ["C" "O" "N" "Cl"]
                :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]}"#
        ),
        Deltas::from_iter([Delta::Bond(BondDelta::ModifyField {
            id: BondId(0),
            change: BondFieldChange::Order {
                old: NumForm::Lit(1),
                new: NumForm::Lit(2),
            },
        })]),
    );

    vec![
        SvgCase::new("labeled_atoms", molecule_depiction(&labeled)),
        SvgCase::new("tetrahedral_stereo", molecule_depiction(&tetrahedral)),
        SvgCase::new("fused_aromatic", molecule_depiction(&fused_aromatic)),
        SvgCase::new(
            "mapped_reaction",
            reaction
                .depict()
                .expect("benchmark reaction must be depictable"),
        ),
    ]
}

fn bench_svg_render(c: &mut Criterion) {
    let mut group = c.benchmark_group("svg/render");
    for atom_count in [8, 128] {
        let depiction = chain_depiction(atom_count);
        group.throughput(Throughput::Elements(atom_count as u64));
        group.bench_with_input(
            BenchmarkId::from_parameter(atom_count),
            &depiction,
            |b, depiction| {
                b.iter(|| black_box(black_box(depiction).render_svg()));
            },
        );
    }
    group.finish();

    let mut group = c.benchmark_group("svg/render/representative");
    for case in svg_cases() {
        group.bench_with_input(case.name, &case.depiction, |b, depiction| {
            b.iter(|| black_box(black_box(depiction).render_svg()));
        });
    }
    group.finish();
}

criterion_group!(svg, bench_svg_render);
criterion_main!(svg);
