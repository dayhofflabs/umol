//! Separate common and rare reading costs; no assumed workload frequencies.

use std::alloc::{GlobalAlloc, Layout, System};
use std::env;
use std::hint::black_box;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use criterion::Criterion;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    ArithExpr, AtomConstraintForm, AtomConstraintKey, AtomConstraintsForm, AtomForm, AtomId,
    BondForm, BondId, Lattice, Molecule, MoleculeEntries, Normalize, NumForm,
    StereoAtomConstraintForm, StereoAtomConstraintKey, StereoAtomForm, StereoAtomId,
    StereoBondConstraintForm, StereoBondConstraintKey, StereoBondForm, StereoBondId,
    StereoConfigurationForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    Stereogenicity, StereogenicityForm,
};

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
    if !allocations {
        criterion.final_summary();
    }
}
