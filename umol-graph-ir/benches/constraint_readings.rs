//! Separate common and rare reading costs; no assumed workload frequencies.

use std::alloc::{GlobalAlloc, Layout, System};
use std::env;
use std::hint::black_box;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use criterion::Criterion;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    ArithExpr, AtomConstraintForm, AtomConstraintKey, AtomConstraintsForm, AtomForm, AtomId,
    BondForm, Molecule, MoleculeEntries, NumForm,
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
    if !allocations {
        criterion.final_summary();
    }
}
