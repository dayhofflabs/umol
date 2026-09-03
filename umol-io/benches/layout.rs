//! CoordGen molecule-layout benchmarks over fixed graph-IR fixtures.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use umol_chem::element::Element;
use umol_coordgen_sys::COORDGEN_VERSION;
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomForm, AtomId, BondForm, ElementForm, Molecule, MoleculeEntries,
};
use umol_graph_ir::mol_dsl;
use umol_io::layout::{layout_molecule, MoleculeLayoutAlgorithm};

const ALGORITHM: MoleculeLayoutAlgorithm = MoleculeLayoutAlgorithm::CoordGen;

struct LayoutCase {
    category: &'static str,
    name: &'static str,
    molecule: Molecule,
}

impl LayoutCase {
    fn new(category: &'static str, name: &'static str, molecule: Molecule) -> Self {
        let layout = layout_molecule(&molecule, ALGORITHM).expect("benchmark fixture must lay out");
        layout
            .check_frame(&molecule)
            .expect("benchmark fixture must preserve its atom frame");
        Self {
            category,
            name,
            molecule,
        }
    }

    fn benchmark_id(&self) -> BenchmarkId {
        BenchmarkId::new(self.category, self.name)
    }
}

fn atom(element: Element) -> AtomForm {
    AtomForm::from_element(element)
}

fn molecule(atoms: Vec<AtomForm>, edges: &[(usize, usize, u8)]) -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds: edges
            .iter()
            .map(|&(atom_0, atom_1, order)| {
                (
                    AtomId::from(atom_0),
                    AtomId::from(atom_1),
                    BondForm::from_order(order),
                )
            })
            .collect(),
        ..Default::default()
    })
}

fn cases() -> Vec<LayoutCase> {
    let acyclic = molecule(
        [
            Element::C,
            Element::N,
            Element::O,
            Element::C,
            Element::N,
            Element::O,
            Element::C,
            Element::N,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        &[
            (0, 1, 1),
            (1, 2, 1),
            (1, 3, 2),
            (3, 4, 1),
            (4, 5, 2),
            (4, 6, 1),
            (6, 7, 2),
        ],
    );

    let cyclic = molecule(
        vec![atom(Element::C); 8],
        &[
            (0, 1, 1),
            (1, 2, 1),
            (2, 3, 1),
            (3, 4, 1),
            (4, 5, 1),
            (5, 6, 1),
            (6, 7, 1),
            (7, 0, 1),
        ],
    );

    let mut aromatic_entries = MoleculeEntries {
        atoms: vec![atom(Element::C); 6],
        bonds: (0..6)
            .map(|atom_0| {
                (
                    AtomId::from(atom_0),
                    AtomId::from((atom_0 + 1) % 6),
                    BondForm::from_order(1),
                )
            })
            .collect(),
        ..Default::default()
    };
    aromatic_entries.aromatic.push((
        (0..6).map(AtomId::from).collect(),
        AromaticSystemForm::from_electrons(vec![1; 6]),
    ));
    let aromatic = Molecule::from_entries(aromatic_entries);

    let disconnected = molecule(
        [
            Element::C,
            Element::C,
            Element::C,
            Element::O,
            Element::N,
            Element::O,
            Element::H,
            Element::H,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        &[(0, 1, 1), (1, 2, 1), (3, 4, 2), (4, 5, 2), (6, 7, 1)],
    );

    let underdetermined = molecule(
        vec![AtomForm::new(ElementForm::Undetermined); 8],
        &[
            (0, 1, 1),
            (1, 2, 2),
            (2, 3, 1),
            (3, 4, 3),
            (4, 5, 1),
            (5, 6, 2),
            (6, 7, 2),
        ],
    );

    let cis_trans_z = mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]}"#
    );
    let cis_trans_e = mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
    );

    // These shapes come from the atom-mapping benchmark. Their mapping hardness is caused by
    // symmetry and output multiplicity; the layout benchmark records how CoordGen handles the same
    // molecular inputs without treating mapping cost as a layout property.
    let complete_edges = (0..7)
        .flat_map(|atom_0| (atom_0 + 1..7).map(move |atom_1| (atom_0, atom_1, 1)))
        .collect::<Vec<_>>();
    let high_symmetry = molecule(vec![atom(Element::C); 7], &complete_edges);
    let repeated_components = molecule(
        vec![atom(Element::C); 6],
        &[(0, 1, 1), (2, 3, 1), (4, 5, 1)],
    );

    vec![
        LayoutCase::new("acyclic", "asymmetric_tree_8", acyclic),
        LayoutCase::new("cyclic", "cyclooctane", cyclic),
        LayoutCase::new("aromatic", "benzene", aromatic),
        LayoutCase::new("disconnected", "mixed_components_8", disconnected),
        LayoutCase::new("underdetermined", "wildcard_path_8", underdetermined),
        LayoutCase::new("cis_trans", "z_but_2_ene", cis_trans_z),
        LayoutCase::new("cis_trans", "e_but_2_ene", cis_trans_e),
        LayoutCase::new(
            "mapping_hard_tail",
            "high_symmetry_complete_7",
            high_symmetry,
        ),
        LayoutCase::new(
            "mapping_hard_tail",
            "repeated_components_3x2",
            repeated_components,
        ),
    ]
}

fn bench_layout(c: &mut Criterion) {
    let cases = cases();
    let mut group = c.benchmark_group(format!("molecule_layout/coordgen-{COORDGEN_VERSION}"));

    for case in &cases {
        group.throughput(Throughput::Elements(case.molecule.atoms().count() as u64));
        group.bench_with_input(case.benchmark_id(), &case.molecule, |b, molecule| {
            b.iter(|| {
                black_box(
                    layout_molecule(black_box(molecule), ALGORITHM)
                        .expect("validated benchmark fixture must lay out"),
                )
            })
        });
    }
    group.finish();
}

criterion_group!(layout, bench_layout);
criterion_main!(layout);
