use super::*;
use crate::diagnostics::Category;
use crate::io::ir::BondOrder;
use crate::io::ir::builder::{AtomData, BondData, MoleculeBuilder};
use umol_data::Element;

fn build_linear(m: &mut MoleculeBuilder, elements: &[Element], orders: &[BondOrder]) -> Molecule {
    assert!(elements.len() >= 1);
    assert_eq!(orders.len() + 1, elements.len());
    let mut ids = Vec::with_capacity(elements.len());
    for &e in elements {
        let id = m.on_atom(AtomData {
            element: e,
            isotope: None,
            charge: None,
            hydrogen_count: None,
            class: None,
            aromatic: false,
            implicit_h: false,
            chirality: None,
            unknown_symbol: false,
        });
        ids.push(id);
    }
    for (i, &ord) in orders.iter().enumerate() {
        m.on_bond(ids[i], ids[i + 1], BondData { order: ord, dir: None });
    }
    let mut mols = m.finish();
    mols.pop().unwrap()
}

#[test]
fn topo_self_loop_and_parallel_edges() {
    let mut mb = MoleculeBuilder::with_capacity(4, 4);
    let c1 = mb.on_atom_fast(Element::C, true, false);
    let c2 = mb.on_atom_fast(Element::C, true, false);
    // self-loop on c1
    mb.on_bond(c1, c1, BondData { order: BondOrder::Single, dir: None });
    // two edges between c1-c2
    mb.on_bond(c1, c2, BondData { order: BondOrder::Single, dir: None });
    mb.on_bond(c1, c2, BondData { order: BondOrder::Double, dir: None });
    let mut mols = mb.finish();
    let mol = mols.pop().unwrap();

    let mut report = DiagnosticsReport::new();
    let _ = check_topology(&mol, None, &mut report, 0);
    let codes: Vec<&str> = report.diagnostics.iter().map(|d| d.code.0).collect();
    assert!(codes.contains(&"TOPO_SELF_LOOP"));
    assert!(codes.contains(&"TOPO_PARALLEL_EDGES"));
}

#[test]
fn valence_pattern_match() {
    // Ethane: C-C, each C total_valence = 1; states (C: [2,4]) → pick 2 → implicit H = 1
    let mut mb = MoleculeBuilder::with_capacity(2, 1);
    let mol = build_linear(&mut mb, &[Element::C, Element::C], &[BondOrder::Single]);

    let model = ValenceModel::simple_organic();
    let cfg = ValenceConfig { enabled: true, overflow_policy: OverflowPolicy::Error, check_bracket: true, infer_bracket_implicit: false, aromatic_as_one: true, patterns_enabled: true, no_match_policy: OverflowPolicy::Off, ambiguous_match_policy: OverflowPolicy::Off };
    let mut report = DiagnosticsReport::new();
    check_valence(&mol, None, &mut report, 0, &model, &cfg);
    assert!(report.diagnostics.iter().all(|d| d.category != Category::Valence));
}

#[test]
fn valence_pattern_mismatch() {
    // Carbon with triple and double bond (spurious): total_valence for C0 = 5
    let mut mb = MoleculeBuilder::with_capacity(3, 2);
    let c1 = mb.on_atom_fast(Element::C, true, false);
    let c2 = mb.on_atom_fast(Element::C, true, false);
    let c3 = mb.on_atom_fast(Element::C, true, false);
    mb.on_bond(c1, c2, BondData { order: BondOrder::Triple, dir: None });
    mb.on_bond(c1, c3, BondData { order: BondOrder::Double, dir: None });
    let mut mols = mb.finish();
    let mol = mols.pop().unwrap();

    let mut model = ValenceModel::simple_organic();
    // Use numeric fallback only: set states for C and disable patterns so we get overflow
    model.set_states(Element::C, vec![4]);
    let cfg = ValenceConfig { enabled: true, overflow_policy: OverflowPolicy::Warn, check_bracket: true, infer_bracket_implicit: false, aromatic_as_one: true, patterns_enabled: false, no_match_policy: OverflowPolicy::Off, ambiguous_match_policy: OverflowPolicy::Off };
    let mut report = DiagnosticsReport::new();
    check_valence(&mol, None, &mut report, 0, &model, &cfg);
    assert!(report.diagnostics.iter().any(|d| d.code.0 == "VALENCE_EXCEEDS_MAX"));
}

#[test]
fn valence_bracket_h_match_and_mismatch() {
    // Build bracket carbon with explicit H where single bond sum=1; pattern table suggests implicit 3 for neutral C.
    // Match: H3; Mismatch: H1
    // We'll create two molecules: one with H1 (ok) and one with H3 (mismatch).

    // Helper to create a single bracket carbon with explicit Hn followed by a plain carbon
    let build_with_h = |h: u8| {
        let mut mb = MoleculeBuilder::with_capacity(2, 1);
        // Bracket carbon with explicit H
        let c0 = mb.on_atom(AtomData {
            element: Element::C,
            isotope: None,
            charge: None,
            hydrogen_count: Some(h),
            class: None,
            aromatic: false,
            implicit_h: false,
            chirality: None,
            unknown_symbol: false,
        });
        // Neighbor carbon
        let c1 = mb.on_atom(AtomData {
            element: Element::C,
            isotope: None,
            charge: None,
            hydrogen_count: None,
            class: None,
            aromatic: false,
            implicit_h: false,
            chirality: None,
            unknown_symbol: false,
        });
        mb.on_bond(c0, c1, BondData { order: BondOrder::Single, dir: None });
        mb.finish().pop().unwrap()
    };

    let model = ValenceModel::simple_organic();
    let cfg = ValenceConfig { enabled: true, overflow_policy: OverflowPolicy::Error, check_bracket: true, infer_bracket_implicit: false, aromatic_as_one: true, patterns_enabled: true, no_match_policy: OverflowPolicy::Off, ambiguous_match_policy: OverflowPolicy::Off };

    // Match: H3 OK
    let mol_ok = build_with_h(3);
    let mut report_ok = DiagnosticsReport::new();
    check_valence(&mol_ok, None, &mut report_ok, 0, &model, &cfg);
    assert!(report_ok.diagnostics.iter().all(|d| d.code.0 != "VALENCE_BRACKET_H_MISMATCH"));

    // Mismatch: H1 vs implied 3
    let mol_bad = build_with_h(1);
    let mut report_bad = DiagnosticsReport::new();
    let _ = check_valence(&mol_bad, None, &mut report_bad, 0, &model, &cfg);
    assert!(report_bad.diagnostics.iter().any(|d| d.code.0 == "VALENCE_BRACKET_H_MISMATCH"));
}

#[test]
fn arom_inconsistent_lowercase_when_no_fractional_bonding() {
    // Cyclohexane marked aromatic should be inconsistent under HMO (no π edges here)
    let mut mb = MoleculeBuilder::with_capacity(6, 6);
    let a = (0..6)
        .map(|_| mb.on_atom_fast(Element::C, true, true)) // aromatic=true
        .collect::<Vec<_>>();
    for i in 0..6 {
        let j = (i + 1) % 6;
        mb.on_bond(a[i], a[j], BondData { order: BondOrder::Single, dir: None });
    }
    let mut mols = mb.finish();
    let mol = mols.pop().unwrap();

    let a_cfg = AromaticityConfig { enabled: true, ..Default::default() };
    let a_model = AromaticityModel::default();
    let mut report = DiagnosticsReport::new();
    let _ = check_aromaticity(&mol, None, &mut report, 0, &a_model, &a_cfg);
    assert!(report.diagnostics.iter().any(|d| d.code.0 == "AROM_INCONSISTENT_LOWERCASE"));
}

#[test]
fn style_prefer_aromatic_form_when_fractional_bonding() {
    // Benzene in Kekulé: alternating double/single; expect preference for aromatic form
    let mut mb = MoleculeBuilder::with_capacity(6, 6);
    let a = (0..6)
        .map(|_| mb.on_atom_fast(Element::C, true, false)) // aromatic=false
        .collect::<Vec<_>>();
    for i in 0..6 {
        let j = (i + 1) % 6;
        let ord = if i % 2 == 0 { BondOrder::Double } else { BondOrder::Single };
        mb.on_bond(a[i], a[j], BondData { order: ord, dir: None });
    }
    let mut mols = mb.finish();
    let mol = mols.pop().unwrap();

    let a_cfg = AromaticityConfig { enabled: true, ..Default::default() };
    let a_model = AromaticityModel::default();
    let mut report = DiagnosticsReport::new();
    let _ = check_aromaticity(&mol, None, &mut report, 0, &a_model, &a_cfg);
    assert!(report.diagnostics.iter().any(|d| d.code.0 == "STYLE_PREFER_AROMATIC_FORM"));
}

