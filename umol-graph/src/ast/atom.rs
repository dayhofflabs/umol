//! Atom-level mechanical lifting from table IR.

use umol_ast::ast::atom::{AtomAst, ImplicitHydrogensAst, IsotopeAst};
use umol_ast::ast::spin::SpinStateAst;
use umol_ast::ast::value::ValueAst;

pub use umol_ast::ast::atom::AtomAst as _AtomAstReexport;

use crate::table_ir::atom::{Atom as TableAtom, ImplicitHydrogens};

/// Lift the base fields of a `table_ir::Atom` to an `AtomAst`. Topology- and
/// chemistry-derived fields (`valence`, `aromatic`) are lifted at the molecule
/// level so they land in the molecule constraint vec.
pub fn atom_from_table(atom: &TableAtom) -> AtomAst {
    let mut ast = AtomAst::from_element(atom.element);
    if let Some(mass) = atom.isotope_mass {
        ast.isotope_mass = IsotopeAst::Lit(mass);
    }
    if let Some(charge) = atom.charge {
        ast.charge = ValueAst::Lit(charge as i64);
    }
    match atom.implicit_hydrogens {
        Some(ImplicitHydrogens::Hydrogens(h)) => {
            ast.implicit_hydrogens = ImplicitHydrogensAst::Value(ValueAst::Lit(h as i64));
        }
        Some(ImplicitHydrogens::Normal) => {
            ast.implicit_hydrogens = ImplicitHydrogensAst::Normal;
        }
        None => {}
    }
    if let Some(lp) = atom.lone_pairs {
        ast.lone_pairs = ValueAst::Lit(lp as i64);
    }
    let u = atom
        .unpaired_electrons
        .map_or(ValueAst::Undetermined, |u| ValueAst::Lit(u as i64));
    let m = atom
        .multiplicity
        .map_or(ValueAst::Undetermined, |m| ValueAst::Lit(m.multiplicity() as i64));
    if !matches!(u, ValueAst::Undetermined) || !matches!(m, ValueAst::Undetermined) {
        ast.spin = SpinStateAst::from_values(u, m);
    }
    ast
}

// TODO: Verify Aromatic valence mode: Aromatic
pub fn coerce_atom(ast: &mut AtomAst, cfg: &AtomAstConfig) {
    if matches!(ast.isotope_mass, IsotopeAst::Undetermined) {
        ast.isotope_mass = match cfg.isotope_mode {
            IsotopeMode::Natural => IsotopeAst::Natural,
            IsotopeMode::Required => IsotopeAst::Undetermined,
        };
    }
    if matches!(ast.charge, ValueAst::Undetermined) {
        ast.charge = match cfg.charge_mode {
            NumericMode::Zero => ValueAst::Lit(0),
            NumericMode::Required => ValueAst::Undetermined,
        };
    }
    if matches!(ast.implicit_hydrogens, ImplicitHydrogensAst::Undetermined) {
        ast.implicit_hydrogens = match cfg.implicit_h_mode {
            ImplicitHydrogenMode::Normal => ImplicitHydrogensAst::Normal,
            ImplicitHydrogenMode::Zero => ImplicitHydrogensAst::Value(ValueAst::Lit(0)),
            ImplicitHydrogenMode::Required => ImplicitHydrogensAst::Undetermined,
        };
    }
    let SpinStateAst {
        unpaired,
        multiplicity,
    } = mem::take(&mut ast.spin);
    let resolved_u = if matches!(unpaired, ValueAst::Undetermined) {
        match cfg.unpaired_electrons_mode {
            UnpairedElectronsMode::Zero => ValueAst::Lit(0),
            UnpairedElectronsMode::Required => ValueAst::Undetermined,
            UnpairedElectronsMode::Derived => match &multiplicity {
                ValueAst::Lit(m) => ValueAst::Lit(m - 1),
                _ => ValueAst::Undetermined,
            },
        }
    } else {
        unpaired
    };
    let resolved_m = if matches!(multiplicity, ValueAst::Undetermined) {
        match cfg.multiplicity_mode {
            MultiplicityMode::Required => ValueAst::Undetermined,
            MultiplicityMode::Derived => match &resolved_u {
                ValueAst::Lit(u) => ValueAst::Lit(u + 1),
                _ => ValueAst::Undetermined,
            },
        }
    } else {
        multiplicity
    };
    ast.spin = SpinStateAst::from_values(resolved_u, resolved_m);
    let coerce_numeric = |v: &mut ValueAst, mode: &NumericMode| {
        if matches!(*v, ValueAst::Undetermined) {
            *v = match mode {
                NumericMode::Zero => ValueAst::Lit(0),
                NumericMode::Required => ValueAst::Undetermined,
            };
        }
    };
    coerce_numeric(&mut ast.lone_pairs, &cfg.lone_pairs_mode);
}

/// Collapse fields back to `Undetermined` where the current value is what
/// `coerce_atom` would have produced. Call after solving to restore roundtrip
/// fidelity with the DSL.
pub fn release_atom(ast: &mut AtomAst, cfg: &AtomAstConfig) {
    if matches!(
        (&cfg.isotope_mode, &ast.isotope_mass),
        (IsotopeMode::Natural | IsotopeMode::Required, IsotopeAst::Natural)
    ) {
        ast.isotope_mass = IsotopeAst::Undetermined;
    }
    if matches!(
        (&cfg.charge_mode, &ast.charge),
        (NumericMode::Zero, ValueAst::Lit(0))
    ) {
        ast.charge = ValueAst::Undetermined;
    }
    match (&cfg.implicit_h_mode, &ast.implicit_hydrogens) {
        (ImplicitHydrogenMode::Normal, ImplicitHydrogensAst::Normal) => {
            ast.implicit_hydrogens = ImplicitHydrogensAst::Undetermined;
        }
        (ImplicitHydrogenMode::Zero, ImplicitHydrogensAst::Value(ValueAst::Lit(0))) => {
            ast.implicit_hydrogens = ImplicitHydrogensAst::Undetermined;
        }
        _ => {}
    }
    match mem::take(&mut ast.spin) {
        SpinStateAst::from_state(state) => {
            let u_value = state.unpaired_electrons();
            let m_value = state.multiplicity();
            let derived = m_value.multiplicity() == u_value + 1;
            let u_ast = match cfg.unpaired_electrons_mode {
                UnpairedElectronsMode::Zero if u_value == 0 => ValueAst::Undetermined,
                UnpairedElectronsMode::Derived if derived => match cfg.multiplicity_mode {
                    MultiplicityMode::Required => ValueAst::Undetermined,
                    MultiplicityMode::Derived if u_value == 0 => ValueAst::Undetermined,
                    MultiplicityMode::Derived => ValueAst::Lit(u_value as i64),
                },
                _ => ValueAst::Lit(u_value as i64),
            };
            let m_ast = match cfg.multiplicity_mode {
                MultiplicityMode::Derived if derived => ValueAst::Undetermined,
                _ => ValueAst::Lit(m_value.multiplicity() as i64),
            };
            ast.spin = SpinStateAst::from_values(u_ast, m_ast);
        }
        pair => ast.spin = pair,
    }
    let release_numeric = |v: &mut ValueAst, mode: &NumericMode| {
        if matches!((mode, &*v), (NumericMode::Zero, ValueAst::Lit(0))) {
            *v = ValueAst::Undetermined;
        }
    };
    release_numeric(&mut ast.lone_pairs, &cfg.lone_pairs_mode);
}
