//! Separate common and rare reading costs; no assumed workload frequencies.

use std::alloc::{GlobalAlloc, Layout, System};
use std::env;
use std::hint::black_box;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use criterion::Criterion;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    ArithExpr, AtomConstraintForm, AtomConstraintKey, AtomConstraintsForm, AtomForm, AtomId,
    BondForm, BondId, GraphSymmetryConfig, Lattice, LigandPermutation, Molecule, MoleculeEntries,
    Normalize, NumForm, OrientedLigandPermutation, StereoAtomConstraintForm,
    StereoAtomConstraintKey, StereoAtomForm, StereoAtomId, StereoBondConstraintForm,
    StereoBondConstraintKey, StereoBondForm, StereoBondId, StereoConfigurationForm, StereoCoset,
    StereoKind, StereoLigand, StereoLigandKind, StereoLigandPosition, Stereogenicity,
    StereogenicityForm,
};
use umol_graph_ir::mol_dsl_concrete;
use umol_perm::{Orientation, Permutation};

struct CountingAllocator;
static COUNTING: AtomicBool = AtomicBool::new(false);
static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);
static BYTES: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if COUNTING.load(Ordering::Relaxed) {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
            BYTES.fetch_add(layout.size(), Ordering::Relaxed);
        }
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        if COUNTING.load(Ordering::Relaxed) {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
            BYTES.fetch_add(size, Ordering::Relaxed);
        }
        unsafe { System.realloc(ptr, layout, size) }
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

fn measure<T>(
    criterion: &mut Criterion,
    allocations: bool,
    id: &str,
    mut operation: impl FnMut() -> T,
) {
    if allocations {
        black_box(operation());
        ALLOCATIONS.store(0, Ordering::Relaxed);
        BYTES.store(0, Ordering::Relaxed);
        COUNTING.store(true, Ordering::Relaxed);
        for _ in 0..1000 {
            black_box(operation());
        }
        COUNTING.store(false, Ordering::Relaxed);
        println!(
            "{id}\t{}\t{}",
            ALLOCATIONS.load(Ordering::Relaxed) / 1000,
            BYTES.load(Ordering::Relaxed) / 1000
        );
    } else {
        criterion.bench_function(id, |b| b.iter(&mut operation));
    }
}

macro_rules! stereo_readings {
    ($criterion:expr, $allocations:expr, $name:literal, $kind:expr, $form:ident,
     $constraint:ident, $key:ident, $entity:ident, $id:expr, $entries:ident, $site:expr,
     $bonds:expr, $ligands:expr) => {
        for (group, name, configuration, assertion) in [
            (
                "common",
                "literal_empty",
                StereoConfigurationForm::kinded($kind, 0_u32),
                None,
            ),
            (
                "common",
                "unknown_coset",
                StereoConfigurationForm::kinded($kind, StereoCoset::Undetermined),
                None,
            ),
            (
                "common",
                "literal_assertion",
                StereoConfigurationForm::kinded($kind, 0_u32),
                Some(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
            ),
            (
                "common",
                "kindless_empty",
                StereoConfigurationForm::Undetermined,
                None,
            ),
            (
                "semantic_change",
                "kindless_assertion",
                StereoConfigurationForm::Undetermined,
                Some(StereogenicityForm::Lit(Stereogenicity::Stereogenic)),
            ),
            (
                "rare",
                "singleton_assertion",
                StereoConfigurationForm::kinded($kind, 0_u32),
                Some(StereogenicityForm::lit_set([Stereogenicity::Stereogenic])),
            ),
            (
                "semantic_change",
                "kindless_bottom",
                StereoConfigurationForm::Undetermined,
                Some(StereogenicityForm::lit_set([])),
            ),
        ] {
            let mut form = $form {
                configuration,
                ..Default::default()
            };
            if let Some(value) = assertion {
                form.constraints.set($constraint::Stereogenicity(value));
            }
            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 6],
                bonds: $bonds
                    .into_iter()
                    .map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1)))
                    .collect(),
                $entries: vec![(
                    $site,
                    $ligands
                        .into_iter()
                        .map(|i| StereoLigand::new(AtomId(i), StereoLigandKind::Atom))
                        .collect(),
                    form.clone(),
                )],
                ..Default::default()
            });
            let key = $key::Stereogenicity;
            let top = Some($constraint::Stereogenicity(
                StereogenicityForm::Undetermined,
            ));
            let expected = match name {
                "literal_assertion" | "singleton_assertion" => Some($constraint::Stereogenicity(
                    StereogenicityForm::Lit(Stereogenicity::Stereogenic),
                )),
                _ => top.clone(),
            };
            let reading = molecule.$entity($id).constraints();
            assert_eq!(reading.asserted(key), expected);
            assert_eq!(reading.asserted_complete(key), expected);
            assert_eq!(reading.evidenced(key), expected);
            assert_eq!(reading.evidenced_complete(key), expected);
            assert_eq!(reading.derived(key), top);
            assert_eq!(reading.derived_complete(key), top);
            let prefix = format!("stereo_readings/{}/{group}/{name}", $name);
            macro_rules! reading {
                ($method:ident) => {
                    measure(
                        $criterion,
                        $allocations,
                        &format!("{prefix}/{}", stringify!($method)),
                        || {
                            black_box(&molecule)
                                .$entity($id)
                                .constraints()
                                .$method(black_box(key))
                        },
                    );
                };
            }
            reading!(asserted);
            reading!(asserted_complete);
            reading!(derived);
            reading!(derived_complete);
            reading!(evidenced);
            reading!(evidenced_complete);
            if name != "kindless_bottom" {
                let other = form.clone();
                let normalized = form.clone().normalize().unwrap();
                assert!(form.normalized_eq(&form));
                assert_eq!(form.meet(&form), Some(normalized.clone()));
                assert_eq!(form.join(&form).unwrap(), normalized);
                assert!(form.matches(&form));
                assert!(form.is_compatible(&form));
                macro_rules! operation {
                    ($method:ident) => {
                        measure(
                            $criterion,
                            $allocations,
                            &format!("{prefix}/{}", stringify!($method)),
                            || black_box(&form).$method(black_box(&other)),
                        );
                    };
                }
                operation!(normalized_eq);
                operation!(meet);
                operation!(join);
                operation!(matches);
                operation!(is_compatible);
            }
        }
    };
}

fn derived_stereo_readings(criterion: &mut Criterion, allocations: bool) {
    let cases = [
        (
            "chlorofluoropropane",
            mol_dsl_concrete!(
                r#"{
        :atoms ["C" "C" "C#h3" "F" "Cl" "H" "H" "H"]
        :bonds [[0 1 "1"] [1 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 3 4 5] :attrs "Th0"}
                       {:site 1 :ligands [0 2 6 7] :attrs "Th0"}]
    }"#
            ),
            StereoAtomId(1),
            [2, 3],
        ),
        (
            "meso",
            mol_dsl_concrete!(
                r#"{
        :atoms ["C" "C" "H" "F" "Cl" "H" "F" "Cl"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [1 5 "1"] [1 6 "1"] [1 7 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                       {:site 1 :ligands [0 5 6 7] :attrs "Th1"}]
    }"#
            ),
            StereoAtomId(0),
            [0, 1],
        ),
        (
            "square_planar",
            mol_dsl_concrete!(
                r#"{
           :atoms ["Pt#c-2" "N#c+#h3" "N#c+#h3" "Cl" "Cl"]
           :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
           :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Sp0"}]
       }"#
            ),
            StereoAtomId(0),
            [0, 1],
        ),
        (
            "trigonal_bipyramidal",
            mol_dsl_concrete!(
                r#"{
           :atoms ["P" "Cl" "F" "F" "F" "F"]
           :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [0 5 "1"]]
           :stereo-atoms [{:site 0 :ligands [1 2 3 4 5] :attrs "Tb0"}]
       }"#
            ),
            StereoAtomId(0),
            [1, 2],
        ),
    ];
    for (name, molecule, id, pair) in &cases {
        let symmetry = molecule
            .graph_symmetry(&GraphSymmetryConfig::default())
            .into_graph_symmetry();
        let stereo = molecule.stereo_atom(*id).with_symmetry(&symmetry);
        let permutation = OrientedLigandPermutation {
            permutation: LigandPermutation(Permutation::identity(stereo.ligand_count())),
            orientation: Orientation::Proper,
        };
        measure(
            criterion,
            allocations,
            &format!("stereo_derived_readings/{name}/stereogenicity"),
            || black_box(black_box(stereo).stereogenicity()),
        );
        measure(
            criterion,
            allocations,
            &format!("stereo_derived_readings/{name}/topicity"),
            || {
                black_box(
                    black_box(stereo)
                        .topicity(StereoLigandPosition(pair[0]), StereoLigandPosition(pair[1])),
                )
            },
        );
        measure(
            criterion,
            allocations,
            &format!("stereo_derived_readings/{name}/ligand_symmetry"),
            || black_box(black_box(stereo).ligand_symmetry(black_box(permutation))),
        );
    }
    let molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "C" "F" "Cl" "F" "Cl"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]
    }"#
    );
    let symmetry = molecule
        .graph_symmetry(&GraphSymmetryConfig::default())
        .into_graph_symmetry();
    let stereo = molecule
        .stereo_bond(StereoBondId(0))
        .with_symmetry(&symmetry);
    let permutation = OrientedLigandPermutation {
        permutation: LigandPermutation(Permutation::from_image(&[2, 3, 0, 1])),
        orientation: Orientation::Proper,
    };
    measure(
        criterion,
        allocations,
        "stereo_derived_readings/cis_trans/stereogenicity",
        || black_box(black_box(stereo).stereogenicity()),
    );
    measure(
        criterion,
        allocations,
        "stereo_derived_readings/cis_trans/topicity",
        || black_box(black_box(stereo).topicity(StereoLigandPosition(0), StereoLigandPosition(2))),
    );
    measure(
        criterion,
        allocations,
        "stereo_derived_readings/cis_trans/ligand_symmetry",
        || black_box(black_box(stereo).ligand_symmetry(black_box(permutation))),
    );
}

fn stereogenic_site_kinds(criterion: &mut Criterion, allocations: bool) {
    for (name, molecule) in [
        (
            "inositol",
            mol_dsl_concrete!(
                r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h" "C#h" "O#h" "O#h" "O#h" "O#h" "O#h" "O#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [0 6 "1"] [1 7 "1"] [2 8 "1"] [3 9 "1"] [4 10 "1"] [5 11 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 6 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 2 7 [:h 1]] :attrs "Th*"}
                       {:site 2 :ligands [1 3 8 [:h 2]] :attrs "Th*"}
                       {:site 3 :ligands [2 4 9 [:h 3]] :attrs "Th*"}
                       {:site 4 :ligands [3 5 10 [:h 4]] :attrs "Th*"}
                       {:site 5 :ligands [0 4 11 [:h 5]] :attrs "Th*"}]
    }"#
            ),
        ),
        (
            "tris_hydroxyethyl_methane",
            mol_dsl_concrete!(
                r#"{
        :atoms ["C#h" "C#h" "C#h" "C#h" "C#h3" "C#h3" "C#h3" "O#h" "O#h" "O#h"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [2 5 "1"] [3 6 "1"]
                [1 7 "1"] [2 8 "1"] [3 9 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 [:h 0]] :attrs "Th*"}
                       {:site 1 :ligands [0 4 7 [:h 1]] :attrs "Th*"}
                       {:site 2 :ligands [0 5 8 [:h 2]] :attrs "Th*"}
                       {:site 3 :ligands [0 6 9 [:h 3]] :attrs "Th*"}]
    }"#
            ),
        ),
        (
            "bis_methylcyclohexyl_methanol",
            mol_dsl_concrete!(
                r#"{
        :atoms ["C#h" "C#h2" "C#h2" "C#h" "C#h2" "C#h2" "C#h" "C#h2" "C#h2" "C#h" "C#h2" "C#h2"
                "C#h" "C#h3" "C#h3" "O#h"]
        :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 6 "1"]
                [0 12 "1"] [6 12 "1"] [3 13 "1"] [9 14 "1"] [12 15 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 5 12 [:h 0]] :attrs "Th*"}
                       {:site 3 :ligands [2 4 13 [:h 3]] :attrs "Th*"}
                       {:site 6 :ligands [7 11 12 [:h 6]] :attrs "Th*"}
                       {:site 9 :ligands [8 10 14 [:h 9]] :attrs "Th*"}
                       {:site 12 :ligands [0 6 15 [:h 12]] :attrs "Th*"}]
    }"#
            ),
        ),
        (
            "bis_fluoroethyl_difluorohexane",
            mol_dsl_concrete!(
                r#"{
        :atoms ["C#h3" "C#h3" "C#h3" "C#h3" "C#h" "C#h" "C#h" "C#h" "C#h" "C#h"
                "F" "F" "F" "F"]
        :bonds [[8 9 "1"] [8 4 "1"] [8 5 "1"] [9 6 "1"] [9 7 "1"] [4 0 "1"] [5 1 "1"]
                [6 2 "1"] [7 3 "1"] [4 10 "1"] [5 11 "1"] [6 12 "1"] [7 13 "1"]]
        :stereo-atoms [{:site 4 :ligands [0 8 10 [:h 4]] :attrs "Th*"}
                       {:site 5 :ligands [1 8 11 [:h 5]] :attrs "Th*"}
                       {:site 6 :ligands [2 9 12 [:h 6]] :attrs "Th*"}
                       {:site 7 :ligands [3 9 13 [:h 7]] :attrs "Th*"}
                       {:site 8 :ligands [4 5 9 [:h 8]] :attrs "Th*"}
                       {:site 9 :ligands [6 7 8 [:h 9]] :attrs "Th*"}]
    }"#
            ),
        ),
    ] {
        let symmetry = molecule
            .graph_symmetry(&GraphSymmetryConfig::default())
            .into_graph_symmetry();
        measure(
            criterion,
            allocations,
            &format!("stereogenic_site_kinds/{name}"),
            || molecule.stereogenic_site_kinds(black_box(&symmetry)),
        );
    }
}

fn chfclbr_readings(criterion: &mut Criterion, allocations: bool) {
    let molecule = mol_dsl_concrete!(
        r#"{
        :atoms ["C" "H" "F" "Cl" "Br"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}]
    }"#
    );
    let config = GraphSymmetryConfig::default();
    let symmetry = molecule.graph_symmetry(&config).into_graph_symmetry();
    let stereo = molecule
        .stereo_atom(StereoAtomId(0))
        .with_symmetry(&symmetry);
    measure(criterion, allocations, "chfclbr/stereogenicity", || {
        black_box(stereo).stereogenicity()
    });
    measure(
        criterion,
        allocations,
        "chfclbr/stereogenic_site_kinds",
        || molecule.stereogenic_site_kinds(black_box(&symmetry)),
    );
    measure(criterion, allocations, "chfclbr/prepare", || {
        black_box(&molecule)
            .graph_symmetry(black_box(&config))
            .into_graph_symmetry()
    });
    measure(
        criterion,
        allocations,
        "chfclbr/prepare_and_stereogenic_site_kinds",
        || {
            black_box(&molecule)
                .graph_symmetry(black_box(&config))
                .stereogenic_site_kinds()
        },
    );
}

fn main() {
    let allocations = env::args().any(|arg| arg == "--allocations");
    let mut criterion = if allocations {
        Criterion::default()
    } else {
        Criterion::default().configure_from_args()
    };
    for (group, name, assertion, order) in [
        ("common", "absent", None, NumForm::Lit(1)),
        (
            "common",
            "top",
            Some(NumForm::Undetermined),
            NumForm::Undetermined,
        ),
        ("common", "literal", Some(NumForm::Lit(1)), NumForm::Lit(1)),
        (
            "rare",
            "singleton",
            Some(NumForm::lit_set([1])),
            NumForm::lit_set([1]),
        ),
        (
            "rare",
            "set",
            Some(NumForm::lit_set([1, 2])),
            NumForm::lit_set([1, 2]),
        ),
        (
            "rare",
            "expression",
            Some(NumForm::arith_expr(ArithExpr::Sum(vec![
                ArithExpr::Lit(0),
                ArithExpr::Lit(1),
            ]))),
            NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(0), ArithExpr::Lit(1)])),
        ),
    ] {
        let mut atom = AtomForm::from_element(Element::C);
        if let Some(assertion) = assertion {
            atom.constraints.set(AtomConstraintForm::valence(assertion));
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom, AtomForm::from_element(Element::C)],
            bonds: vec![(
                AtomId(0),
                AtomId(1),
                BondForm {
                    order,
                    ..Default::default()
                },
            )],
            ..Default::default()
        });
        let pattern = AtomConstraintsForm::from(AtomConstraintForm::valence(1));
        for operation in [
            "asserted",
            "derived",
            "derived_complete",
            "valence",
            "satisfies",
        ] {
            let read = || {
                let atom = black_box(&molecule).atom(AtomId(0));
                let constraints = atom.constraints();
                let key = black_box(AtomConstraintKey::Valence);
                match operation {
                    "asserted" => {
                        black_box(constraints.asserted(key));
                    }
                    "derived" => {
                        black_box(constraints.derived(key));
                    }
                    "derived_complete" => {
                        black_box(constraints.derived_complete(key));
                    }
                    "valence" => {
                        black_box(atom.valence());
                    }
                    "satisfies" => {
                        black_box(constraints.satisfies(black_box(&pattern)));
                    }
                    _ => unreachable!(),
                }
            };
            let id = format!("constraint_readings/{group}/{name}/{operation}");
            if allocations {
                read();
                ALLOCATIONS.store(0, Ordering::Relaxed);
                BYTES.store(0, Ordering::Relaxed);
                COUNTING.store(true, Ordering::Relaxed);
                for _ in 0..1000 {
                    read();
                }
                COUNTING.store(false, Ordering::Relaxed);
                println!(
                    "{id}\t{}\t{}",
                    ALLOCATIONS.load(Ordering::Relaxed) / 1000,
                    BYTES.load(Ordering::Relaxed) / 1000
                );
            } else {
                criterion.bench_function(&id, |b| b.iter(read));
            }
        }
    }
    stereo_readings!(
        &mut criterion,
        allocations,
        "atom",
        StereoKind::Tetrahedral,
        StereoAtomForm,
        StereoAtomConstraintForm,
        StereoAtomConstraintKey,
        stereo_atom,
        StereoAtomId(0),
        stereo_atoms,
        AtomId(0),
        [(0, 1), (0, 2), (0, 3), (0, 4)],
        [1, 2, 3, 4]
    );
    stereo_readings!(
        &mut criterion,
        allocations,
        "bond",
        StereoKind::CisTrans,
        StereoBondForm,
        StereoBondConstraintForm,
        StereoBondConstraintKey,
        stereo_bond,
        StereoBondId(0),
        stereo_bonds,
        BondId(0),
        [(0, 1), (0, 2), (0, 3), (1, 4), (1, 5)],
        [2, 3, 4, 5]
    );
    derived_stereo_readings(&mut criterion, allocations);
    stereogenic_site_kinds(&mut criterion, allocations);
    chfclbr_readings(&mut criterion, allocations);
    if !allocations {
        criterion.final_summary();
    }
}
