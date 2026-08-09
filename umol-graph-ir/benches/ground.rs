//! `MoleculeAst::is_ground()` benchmarks.
//!
//! Exercises the groundness check — a hot path during Molecule construction
//! and in matcher inner loops — across three representative shapes:
//!
//! - `indole_ground`: realistic, fully-lowered ground molecule (nearly every
//!   field is `Lit`). Represents the common case.
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
    ArithExpr, AtomAst, AtomId, BondAst, ElementAst, IntoIr, IsotopeMassAst, MemOp, MoleculeAst,
    MoleculeEntries, NumForm, PredExpr, RelOp, UnpairedElectronsAst,
};

#[path = "fixtures.rs"]
#[allow(dead_code)]
mod fixtures;
use fixtures::MOL_INDOLE;

fn indole_ground() -> MoleculeAst {
    let dsl = MoleculeDsl::from_edn_str(MOL_INDOLE).unwrap();
    let cfg = MoleculeDefaults::zeroed();
    dsl.into_ir(&cfg)
}

fn indole_with_bool_expr_fields() -> MoleculeAst {
    // Realistic pattern: take the ground indole, then splat a few atom
    // fields with boolean-domain ValueExpr patterns (which short-circuit via
    // is_arithmetic()).
    let mut ast = indole_ground();
    let mut b = ast.edit();
    b.atom_mut(AtomId(0)).ast.charge = NumForm::pred_expr(PredExpr::Rel(
        ArithExpr::Var("c".into()),
        RelOp::Eq,
        ArithExpr::Lit(0),
    ));
    b.atom_mut(AtomId(2)).ast.lone_pairs = NumForm::pred_expr(PredExpr::Mem(
        ArithExpr::Var("n".into()),
        MemOp::In,
        BTreeSet::from([0, 1, 2]),
    ));
    ast = b.build();
    ast
}

fn arith_expr_heavy() -> MoleculeAst {
    // Every numeric field is an arithmetic `ArithExpr` of depth 3 — a
    // non-ground symbolic value, so `is_ground` (literal-only) returns false.
    let arith = || {
        NumForm::arith_expr(ArithExpr::Product(vec![
            ArithExpr::Sum(vec![ArithExpr::Lit(2), ArithExpr::Lit(3)]),
            ArithExpr::Neg(Box::new(ArithExpr::Lit(1))),
        ]))
    };
    let make_atom = |el: Element| AtomAst {
        element: ElementAst::Lit(el),
        isotope_mass: IsotopeMassAst::Lit(12),
        charge: arith(),
        implicit_hydrogens: NumForm::arith_expr(ArithExpr::Neg(Box::new(ArithExpr::Lit(1)))),
        lone_pairs: arith(),
        unpaired_electrons: UnpairedElectronsAst {
            count: arith(),
            multiplicity: arith(),
        },
        constraints: Default::default(),
    };
    let atoms: Vec<AtomAst> = (0..20)
        .map(|i| {
            make_atom(match i % 4 {
                0 => Element::C,
                1 => Element::N,
                2 => Element::O,
                _ => Element::H,
            })
        })
        .collect();
    let bonds: Vec<(AtomId, AtomId, BondAst)> = (0..19)
        .map(|i| {
            let bond = BondAst {
                order: arith(),
                charge: arith(),
                unpaired_electrons: UnpairedElectronsAst::default(),
                constraints: Default::default(),
            };
            (AtomId(i as u32), AtomId(i as u32 + 1), bond)
        })
        .collect();
    MoleculeAst::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

fn bench_is_ground(c: &mut Criterion) {
    let mut g = c.benchmark_group("molecule_ast_is_ground");

    let indole = indole_ground();
    g.bench_function("indole_ground", |b| {
        b.iter(|| black_box(&indole).is_ground())
    });

    let indole_expr = indole_with_bool_expr_fields();
    g.bench_function("indole_bool_expr", |b| {
        b.iter(|| black_box(&indole_expr).is_ground())
    });

    let heavy = arith_expr_heavy();
    g.bench_function("arith_expr_heavy", |b| {
        b.iter(|| black_box(&heavy).is_ground())
    });

    g.finish();
}

criterion_group!(is_ground, bench_is_ground);
criterion_main!(is_ground);
