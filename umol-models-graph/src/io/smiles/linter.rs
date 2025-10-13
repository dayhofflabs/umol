//! Linting for SMILES: collect diagnostics from lexing/parsing.

mod context;
mod emitter;
pub use context::LintContext;
pub use emitter::{Emitter, LintConfig};

use super::checker::{
    check_aromaticity, check_stereo_double, check_topology, check_valence, AromaticityConfig,
    AromaticityModel, ValenceConfig, ValenceModel,
};
use super::diagnostics::{Category, Diagnostic, DiagnosticList, Severity, Span};
use super::parser::parse_smiles;

// SMILES linter, runs post-parse
pub fn lint_smiles(input: &str) -> DiagnosticList {
    let ctx = LintContext::new(input);
    let mut report = DiagnosticList::new();
    let mut emitter = Emitter::new(&mut report);

    // Map parser errors into diagnostics
    let parse_res = parse_smiles(input.as_bytes());
    if let Err(err) = &parse_res {
        // Central mapping for parser errors (full coverage)
        // TODO: Generate the diagnostic here
        emitter.emit(d);
    }
    if let Err(_err) = &parse_res {
        // Legacy block retained intentionally empty after centralization
        match *_err {
            // Ring diagnostics
            _ => {}
        }
    }

    // Style and numeric advisories
    run_style_and_numeric_checks(&ctx, &mut emitter, parse_res.is_ok());

    // Post-parse: topology + stereo checks when parsing succeeds (after emitter releases &mut report)
    if let Ok(ref mol) = parse_res {
        let input_len = input.len();
        report.extend(check_topology(mol, super::config::SmilesCheckFlags::ALL));
        let v_cfg = ValenceConfig::default();
        let v_model = ValenceModel::simple_organic();
        check_valence(mol, &mut report, &v_model, &v_cfg);
        check_stereo_double(mol, None, &mut report, input_len);
        // Aromaticity verification scaffold (HMO/Clar config only)
        let a_cfg = AromaticityConfig::default();
        let a_model = AromaticityModel::default();
        let _ = check_aromaticity(mol, None, &mut report, input_len, &a_model, &a_cfg);
    }

    report
}

fn run_style_and_numeric_checks(ctx: &LintContext, emit: &mut Emitter, only_when_parse_ok: bool) {
    let input = ctx.input;
    let bytes = input.as_bytes();

    // Percent ring index style: %01..%09
    let mut i = 0usize;
    while i + 2 < bytes.len() {
        if bytes[i] == b'%'
            && bytes[i + 1] == b'0'
            && (bytes[i + 2] >= b'1' && bytes[i + 2] <= b'9')
        {
            // emit.emit(Diagnostic {
            //     code: Code("STYLE_UNNECESSARY_PERCENT_RING_INDEX"),
            //     category: Category::Style,
            //     severity: Severity::Warning,
            //     span: Span::new(i, i + 3),
            //     message: "Prefer single-digit ring index for 1..9",
            //     details: None,
            // });
        }
        i += 1;
    }

    // Ring numbering style checks only when parsing succeeded to avoid noise
    if only_when_parse_ok {
        let seq = ring_indices_sequence(bytes);
        // if let Some((first_num, s, e)) = seq.first().copied() {
        //     if first_num != 1 {
        //         emit.emit(Diagnostic {
        //             code: Code("STYLE_FIRST_RING_NOT_ONE"),
        //             category: Category::Style,
        //             severity: Severity::Warning,
        //             span: Span::new(s, e),
        //             message: "Prefer starting ring numbering at 1",
        //             details: None,
        //         });
        //     }
        // }
        // // Non-consecutive jumps
        // let mut last: Option<u32> = None;
        // for (num, s, e) in seq.into_iter() {
        //     if let Some(p) = last {
        //         if num > p + 1 {
        //             emit.emit(Diagnostic {
        //                 code: Code("STYLE_NONCONSECUTIVE_RING_NUMBERING"),
        //                 category: Category::Style,
        //                 severity: Severity::Warning,
        //                 span: Span::new(s, e),
        //                 message: "Non-consecutive ring numbering",
        //                 details: None,
        //             });
        //             break;
        //         }
        //     }
        //     last = Some(num);
        // }
    }
}

fn ring_indices_sequence(bytes: &[u8]) -> Vec<(u32, usize, usize)> {
    let mut res = Vec::new();
    let mut i = 0usize;
    let n = bytes.len();
    let mut in_brkt = false;
    while i < n {
        let b = bytes[i];
        if b == b'[' {
            in_brkt = true;
            i += 1;
            continue;
        }
        if b == b']' {
            in_brkt = false;
            i += 1;
            continue;
        }
        if in_brkt {
            i += 1;
            continue;
        }
        if b == b'%' {
            if i + 2 < n && bytes[i + 1].is_ascii_digit() && bytes[i + 2].is_ascii_digit() {
                let num = ((bytes[i + 1] - b'0') as u32) * 10 + (bytes[i + 2] - b'0') as u32;
                res.push((num, i, i + 3));
                i += 3;
                continue;
            }
        } else if b.is_ascii_digit() {
            let num = (b - b'0') as u32;
            res.push((num, i, i + 1));
            i += 1;
            continue;
        }
        i += 1;
    }
    res
}
