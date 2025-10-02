use std::collections::HashMap;

use crate::diagnostics::{Category, Code, DiagnosticsReport, Severity, Span};
use crate::io::ir::{AtomSymbol, BondOrder, BondSymbol, Molecule};
use umol_data::Element;

use super::SideChannel;

#[derive(Clone, Copy)]
pub enum OverflowPolicy {
    Off,
    Warn,
    Error,
}

#[derive(Clone, Copy)]
pub struct ValenceConfig {
    pub enabled: bool,
    pub overflow_policy: OverflowPolicy,
    pub check_bracket: bool,
    pub infer_bracket_implicit: bool,
    pub aromatic_as_one: bool,
    pub patterns_enabled: bool,
    pub no_match_policy: OverflowPolicy,
    pub ambiguous_match_policy: OverflowPolicy,
}

impl Default for ValenceConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            overflow_policy: OverflowPolicy::Warn,
            check_bracket: true,
            infer_bracket_implicit: false,
            aromatic_as_one: true,
            patterns_enabled: true,
            no_match_policy: OverflowPolicy::Off,
            ambiguous_match_policy: OverflowPolicy::Off,
        }
    }
}

#[derive(Default)]
pub struct ValenceModel {
    states: HashMap<Element, Vec<u8>>, // numeric fallback (empty initially)
    pub patterns: ValencePatternTable, // explicit pattern table
}

impl ValenceModel {
    pub fn simple_organic() -> Self {
        Self { states: HashMap::new(), patterns: ValencePatternTable::organic_strict() }
    }
    pub fn states_for(&self, e: Element) -> Option<&[u8]> { self.states.get(&e).map(|v| &v[..]) }
    pub fn set_states(&mut self, e: Element, states: Vec<u8>) { self.states.insert(e, states); }
}

#[derive(Clone, Copy, Debug)]
pub struct ValencePattern {
    pub element: Option<Element>,
    pub bond_sum: Option<u8>,
    pub charge: Option<i8>,
    pub implicit_h: Option<u8>,
    pub unpaired: Option<u8>,
}

#[derive(Default)]
pub struct ValencePatternTable { pub patterns: Vec<ValencePattern> }

impl ValencePatternTable {
    #[rustfmt::skip]
    pub fn organic_strict() -> Self {
        use umol_data::e;
        Self { patterns: vec![
            // H
            ValencePattern { element: Some(e!(H)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) },
            ValencePattern { element: Some(e!(H)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(0), unpaired: Some(1) },
            ValencePattern { element: Some(e!(H)), bond_sum: Some(0), charge: Some(1), implicit_h: Some(0), unpaired: Some(0) },
            ValencePattern { element: Some(e!(H)), bond_sum: Some(0), charge: Some(-1), implicit_h: Some(0), unpaired: Some(0) },
            // C (subset)
            ValencePattern { element: Some(e!(C)), bond_sum: Some(4), charge: Some(0), implicit_h: Some(0), unpaired: Some(0) },
            ValencePattern { element: Some(e!(C)), bond_sum: Some(3), charge: Some(0), implicit_h: Some(1), unpaired: Some(0) },
            ValencePattern { element: Some(e!(C)), bond_sum: Some(2), charge: Some(0), implicit_h: Some(2), unpaired: Some(0) },
            ValencePattern { element: Some(e!(C)), bond_sum: Some(1), charge: Some(0), implicit_h: Some(3), unpaired: Some(0) },
            ValencePattern { element: Some(e!(C)), bond_sum: Some(0), charge: Some(0), implicit_h: Some(4), unpaired: Some(0) },
        ]}
    }
}

fn match_score(p: &ValencePattern, element: Element, bond_sum: u8, charge: i8, unpaired: Option<u8>) -> Option<usize> {
    if let Some(e) = p.element { if e != element { return None; } }
    if let Some(bs) = p.bond_sum { if bs != bond_sum { return None; } }
    if let Some(ch) = p.charge { if ch != charge { return None; } }
    match (p.unpaired, unpaired) { (Some(pu), Some(u)) if pu != u => return None, _ => {} }
    let mut score = 0usize;
    if p.element.is_some() { score += 1; }
    if p.bond_sum.is_some() { score += 1; }
    if p.charge.is_some() { score += 1; }
    if p.unpaired.is_some() && unpaired.is_some() { score += 1; }
    if p.implicit_h.is_some() { score += 1; }
    Some(score)
}

#[derive(Default)]
struct PatternDecision { implicit_h: Option<u8>, matches: usize }

fn resolve_valence_pattern(
    element: Element,
    bond_sum: u8,
    charge: i8,
    unpaired: Option<u8>,
    tbl: &ValencePatternTable,
) -> PatternDecision {
    let mut best_score = 0usize; let mut best_idx: Option<usize> = None; let mut match_count = 0usize;
    for (idx, p) in tbl.patterns.iter().enumerate() {
        if let Some(score) = match_score(p, element, bond_sum, charge, unpaired) {
            match_count += 1; if score > best_score { best_score = score; best_idx = Some(idx); }
        }
    }
    if let Some(i) = best_idx { let p = tbl.patterns[i]; PatternDecision { implicit_h: p.implicit_h, matches: match_count } } else { PatternDecision::default() }
}

pub struct ValenceArtifacts {
    pub atoms_checked: usize,
    pub overflow_count: usize,
    pub bracket_mismatch_count: usize,
}

pub fn check_valence(
    mol: &Molecule,
    _side: Option<&SideChannel>,
    report: &mut DiagnosticsReport,
    input_len: usize,
    model: &ValenceModel,
    cfg: &ValenceConfig,
) -> ValenceArtifacts {
    if !cfg.enabled {
        return ValenceArtifacts { atoms_checked: 0, overflow_count: 0, bracket_mismatch_count: 0 };
    }
    let atom_len = mol.atoms.len();
    let mut order_sum: Vec<u8> = vec![0; atom_len];
    let order_weight = |o: BondOrder| -> u8 { match o {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => if cfg.aromatic_as_one { 1 } else { 1 },
        BondOrder::Unknown => 0,
    }};
    for b in &mol.bonds {
        let (Some(a), Some(c)) = (b.start_atom, b.end_atom) else { continue; };
        if let BondSymbol::Bond(ord) = b.symbol {
            let w = order_weight(ord);
            if (a as usize) < atom_len { order_sum[a as usize] = order_sum[a as usize].saturating_add(w); }
            if (c as usize) < atom_len { order_sum[c as usize] = order_sum[c as usize].saturating_add(w); }
        }
    }
    let mut artifacts = ValenceArtifacts { atoms_checked: 0, overflow_count: 0, bracket_mismatch_count: 0 };
    let atom_element = |idx: usize| -> Option<Element> { match &mol.atoms[idx].symbol { AtomSymbol::Element(e) => Some(*e), _ => None } };
    for (i, atom) in mol.atoms.iter().enumerate() {
        let Some(elem) = atom_element(i) else { continue; };
        let charge_i32 = atom.charge.unwrap_or(0);
        let charge_i8 = charge_i32.clamp(i8::MIN as i32, i8::MAX as i32) as i8;
        let sum_orders_u8 = order_sum[i];
        let sum_orders = sum_orders_u8 as i32;
        let total_valence = sum_orders - charge_i32;
        artifacts.atoms_checked += 1;
        let mut implicit_h_opt: Option<i32> = None;
        if cfg.patterns_enabled {
            let d = resolve_valence_pattern(elem, sum_orders_u8, charge_i8, None, &model.patterns);
            if d.matches == 0 {
                match cfg.no_match_policy {
                    OverflowPolicy::Off => {}
                    OverflowPolicy::Warn => report.push(crate::diagnostics::Diagnostic { code: Code("VALENCE_NO_MATCH"), category: Category::Valence, severity: Severity::Warning, span: Span::new(0, input_len), message: "No valence pattern matched", details: Some(format!("atom_index={}", i)) }),
                    OverflowPolicy::Error => report.push(crate::diagnostics::Diagnostic { code: Code("VALENCE_NO_MATCH"), category: Category::Valence, severity: Severity::Error, span: Span::new(0, input_len), message: "No valence pattern matched", details: Some(format!("atom_index={}", i)) }),
                }
            } else {
                if d.matches > 1 {
                    match cfg.ambiguous_match_policy {
                        OverflowPolicy::Off => {}
                        OverflowPolicy::Warn => report.push(crate::diagnostics::Diagnostic { code: Code("VALENCE_AMBIGUOUS_MATCH"), category: Category::Valence, severity: Severity::Warning, span: Span::new(0, input_len), message: "Multiple valence patterns matched; selected most specific", details: Some(format!("atom_index={}", i)) }),
                        OverflowPolicy::Error => report.push(crate::diagnostics::Diagnostic { code: Code("VALENCE_AMBIGUOUS_MATCH"), category: Category::Valence, severity: Severity::Error, span: Span::new(0, input_len), message: "Multiple valence patterns matched; selected most specific", details: Some(format!("atom_index={}", i)) }),
                    }
                }
                if let Some(h) = d.implicit_h { implicit_h_opt = Some(h as i32); }
            }
        }
        if implicit_h_opt.is_none() {
            let states_opt = model.states_for(elem);
            if states_opt.is_none() { continue; }
            let states = states_opt.unwrap();
            if states.is_empty() { continue; }
            let mut chosen: Option<u8> = None;
            for &s in states { if (s as i32) >= total_valence { chosen = Some(s); break; } }
            if chosen.is_none() {
                artifacts.overflow_count += 1;
                match cfg.overflow_policy {
                    OverflowPolicy::Off => {}
                    OverflowPolicy::Warn => report.push(crate::diagnostics::Diagnostic { code: Code("VALENCE_EXCEEDS_MAX"), category: Category::Valence, severity: Severity::Warning, span: Span::new(0, input_len), message: "Valence exceeds maximum allowed state", details: Some(format!("atom_index={}", i)) }),
                    OverflowPolicy::Error => report.push(crate::diagnostics::Diagnostic { code: Code("VALENCE_EXCEEDS_MAX"), category: Category::Valence, severity: Severity::Error, span: Span::new(0, input_len), message: "Valence exceeds maximum allowed state", details: Some(format!("atom_index={}", i)) }),
                }
                continue;
            }
            let chosen_state = chosen.unwrap() as i32;
            implicit_h_opt = Some(chosen_state - total_valence);
        }
        let implicit_h = implicit_h_opt.unwrap_or(0);
        if cfg.check_bracket {
            if let Some(h_explicit) = atom.hydrogen_count {
                let implied = implicit_h.max(0) as u32;
                if h_explicit != implied {
                    artifacts.bracket_mismatch_count += 1;
                    report.push(crate::diagnostics::Diagnostic { code: Code("VALENCE_BRACKET_H_MISMATCH"), category: Category::Valence, severity: Severity::Error, span: Span::new(0, input_len), message: "Bracket H count mismatches valence-based implicit H", details: Some(format!("atom_index={}, expected_H={}, found_H={}", i, implied, h_explicit)) });
                }
            }
        }
    }
    artifacts
}


