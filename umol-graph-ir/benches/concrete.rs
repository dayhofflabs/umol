//! `Molecule::is_concrete()` benchmarks.
//!
//! Exercises the concreteness check — a hot path during Molecule construction
//! and in matcher inner loops — across three representative shapes:
//!
//! - `indole_concrete`: realistic, fully-lowered concrete molecule (nearly
//!   every field is `Lit`). Represents the common case.
//! - `indole_bool_expr`: indole with a few fields replaced by boolean-domain
//!   `NumForm::PredExpr` patterns (`Rel`, `Mem`).
//! - `arith_expr_heavy`: every numeric field on every atom carries an
//!   arithmetic `NumForm::ArithExpr` of depth 3 (non-ground symbolic values).

use std::collections::BTreeSet;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_chem::element::Element;
use umol_edn::FromEdn;
use umol_graph_ir::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_graph_ir::ir::{
    ArithExpr, AtomForm, AtomId, BondForm, ElementForm, IntoIr, IsotopeMassForm, MemOp, Molecule,
    MoleculeEntries, NumForm, PredExpr, RelOp, UnpairedElectronsForm,
};

#[path = "fixtures.rs"]
#[allow(
    dead_code,
    reason = "this bench includes a shared fixture module but uses only the indole fixture"
)]
mod fixtures;
use fixtures::MOL_INDOLE;

fn indole_concrete() -> Molecule {
    let dsl = MoleculeDsl::from_edn_str(MOL_INDOLE).unwrap();
    let cfg = MoleculeDefaults::concrete();
    dsl.into_ir(&cfg)
}

fn indole_with_bool_expr_fields() -> Molecule {
    // Realistic pattern: take the concrete indole, then splat a few atom
    // fields with boolean-domain ValueExpr patterns (which short-circuit via
    // is_arithmetic()).
    let mut molecule = indole_concrete();
    let mut b = molecule.edit();
    b.atom_mut(AtomId(0)).attributes.charge = NumForm::pred_expr(PredExpr::Rel(
        ArithExpr::Var("c".into()),
        RelOp::Eq,
        ArithExpr::Lit(0),
    ));
    b.atom_mut(AtomId(2)).attributes.lone_pairs = NumForm::pred_expr(PredExpr::Mem(
        ArithExpr::Var("n".into()),
        MemOp::In,
        BTreeSet::from([0, 1, 2]),
    ));
    molecule = b.build();
    molecule
}

fn arith_expr_heavy() -> Molecule {
    // Every numeric field is an arithmetic `ArithExpr` of depth 3 — a
    // non-literal symbolic value, so `is_concrete` (literal-only) returns false.
    let arith = || {
        NumForm::arith_expr(ArithExpr::Product(vec![
            ArithExpr::Sum(vec![ArithExpr::Lit(2), ArithExpr::Lit(3)]),
            ArithExpr::Neg(Box::new(ArithExpr::Lit(1))),
        ]))
    };
    let make_atom = |el: Element| AtomForm {
        element: ElementForm::Lit(el),
        isotope_mass: IsotopeMassForm::Lit(12),
        charge: arith(),
        implicit_hydrogens: NumForm::arith_expr(ArithExpr::Neg(Box::new(ArithExpr::Lit(1)))),
        lone_pairs: arith(),
        unpaired_electrons: UnpairedElectronsForm {
            count: arith(),
            multiplicity: arith(),
        },
        constraints: Default::default(),
    };
    let atoms: Vec<AtomForm> = (0..20)
        .map(|i| {
            make_atom(match i % 4 {
                0 => Element::C,
                1 => Element::N,
                2 => Element::O,
                _ => Element::H,
            })
        })
        .collect();
    let bonds: Vec<(AtomId, AtomId, BondForm)> = (0..19)
        .map(|i| {
            let bond = BondForm {
                order: arith(),
                charge: arith(),
                unpaired_electrons: UnpairedElectronsForm::default(),
                constraints: Default::default(),
            };
            (AtomId(i as u32), AtomId(i as u32 + 1), bond)
        })
        .collect();
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

fn bench_is_concrete(c: &mut Criterion) {
    let mut g = c.benchmark_group("molecule_is_concrete");

    let indole = indole_concrete();
    g.bench_function("indole_concrete", |b| {
        b.iter(|| black_box(&indole).is_concrete())
    });

    let indole_expr = indole_with_bool_expr_fields();
    g.bench_function("indole_bool_expr", |b| {
        b.iter(|| black_box(&indole_expr).is_concrete())
    });

    let heavy = arith_expr_heavy();
    g.bench_function("arith_expr_heavy", |b| {
        b.iter(|| black_box(&heavy).is_concrete())
    });

    g.finish();
}

criterion_group!(is_concrete, bench_is_concrete);
criterion_main!(is_concrete);
