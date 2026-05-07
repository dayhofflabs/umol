//! `MoleculeAst::is_ground()` benchmarks.
//!
//! Exercises the groundness check — a hot path during Molecule construction
//! and in matcher inner loops — across three representative shapes:
//!
//! - `indole_ground`: realistic, fully-lowered ground molecule (nearly every
//!   field is `Lit`). Represents the common case.
//! - `indole_bool_expr`: indole with a few fields replaced by boolean-domain
//!   `ValueAst::Expr` patterns (`Rel`, `Mem`). Simulates matcher workloads
//!   where `is_arithmetic()` short-circuits.
//! - `arith_expr_heavy`: pathological upper bound — every numeric field on
//!   every atom carries an arithmetic `Expr` of depth 3. Exercises the
//!   evaluator walk.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use umol_ast::ast::{
    ArithOp, AtomAst, AtomConstraint, AtomIdx, BondAst, ElementAst, Expr, ImplicitHydrogensAst,
    IntoAst, IsotopeAst, MoleculeAst, RelOp, SpinStateAst, ValueAst,
};
use umol_ast::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_edn::FromEdn;
use umol_shared::element::Element;

#[path = "fixtures.rs"]
#[allow(dead_code)]
mod fixtures;
use fixtures::MOL_INDOLE;

fn indole_ground() -> MoleculeAst {
    let dsl = MoleculeDsl::from_edn_str(MOL_INDOLE).unwrap();
    let cfg = MoleculeDefaults::zeroed();
    dsl.into_ast(&cfg)
}

fn indole_with_bool_expr_fields() -> MoleculeAst {
    // Realistic pattern: take the ground indole, then splat a few atom
    // fields with boolean-domain Expr patterns (which short-circuit via
    // is_arithmetic()).
    let mut ast = indole_ground();
    let mut b = ast.edit();
    b.atom_mut(AtomIdx(0)).charge = ValueAst::Expr(Expr::Rel(
        Box::new(Expr::Var("c".into())),
        RelOp::Eq,
        Box::new(Expr::Lit(0)),
    ));
    b.atom_mut(AtomIdx(2)).lone_pairs =
        ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("n".into())), vec![0, 1, 2]));
    ast = b.build();
    ast
}

fn arith_expr_heavy() -> MoleculeAst {
    // Upper bound: every numeric field is an arithmetic Expr of depth 3.
    // The tree is constant-valued (evaluator can fold it), so semantic
    // is_ground will walk the full tree; syntactic is_ground fails on the
    // first `Expr` encountered.
    let arith = || {
        ValueAst::Expr(Expr::BinOp(
            Box::new(Expr::BinOp(
                Box::new(Expr::Lit(2)),
                ArithOp::Add,
                Box::new(Expr::Lit(3)),
            )),
            ArithOp::Mul,
            Box::new(Expr::Neg(Box::new(Expr::Lit(1)))),
        ))
    };
    let make_atom = |el: Element| AtomAst {
        element: ElementAst::Lit(el),
        isotope_mass: IsotopeAst::Expr(Expr::BinOp(
            Box::new(Expr::Lit(12)),
            ArithOp::Add,
            Box::new(Expr::Lit(0)),
        )),
        charge: arith(),
        implicit_hydrogens: ImplicitHydrogensAst::Expr(Expr::Neg(Box::new(Expr::Lit(1)))),
        lone_pairs: arith(),
        spin: SpinStateAst::from_values(arith(), arith()),
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
    let bonds: Vec<(AtomIdx, AtomIdx, BondAst)> = (0..19)
        .map(|i| {
            let bond = BondAst {
                order: arith(),
                charge: arith(),
                spin: SpinStateAst::default(),
                constraints: Default::default(),
            };
            (AtomIdx(i as u32), AtomIdx(i as u32 + 1), bond)
        })
        .collect();
    MoleculeAst::from_atoms_and_bonds(
        atoms,
        bonds,
    )
}

// Silence `unused` warnings on pre-switch build: `AtomConstraint` is kept
// here for future bench cases involving atom-constraint patterns.
#[allow(dead_code)]
fn _keep_import(c: AtomConstraint) -> AtomConstraint {
    c
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
