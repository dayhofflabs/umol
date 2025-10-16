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
use super::parser::parse_smiles_core;

// SMILES linter, runs post-parse
pub fn lint_smiles(input: &str) -> DiagnosticList {
    let ctx = LintContext::new(input);
    let mut report = DiagnosticList::new();
    let mut emitter = Emitter::new(&mut report);

    // Map parser errors into diagnostics
    let parse_res = parse_smiles_core(input.as_bytes(), &super::config::SmilesParseFlags::default());
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
    if let Ok(ref parse_ok) = parse_res {
        let mol = &parse_ok.sir;
        let _input_len = input.len();
        let mut annotations = super::checker::Annotations::default();
        // Topology
        super::checker::topology::check_topology(
            mol,
            &super::config::SmilesCheckFlags::ALL,
            &mut annotations,
            &mut report,
        );
        // Valence
        let v_cfg = ValenceConfig::default();
        let v_model = ValenceModel::simple_organic();
        check_valence(
            mol,
            &super::config::SmilesCheckFlags::ALL,
            &mut annotations,
            &mut report,
            &v_model,
            &v_cfg,
        );
        // Stereo
        super::checker::stereo_double::check_stereo_double(
            mol,
            &super::config::SmilesCheckFlags::ALL,
            &mut annotations,
            &mut report,
        );
        // Aromaticity verification scaffold (HMO/Clar config only)
        // let a_cfg = AromaticityConfig::default();
        // let a_model = AromaticityModel::default();
        // super::checker::aromaticity::check_aromaticity(mol, &super::config::SmilesCheckFlags::ALL, &mut annotations, &mut report);
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

// -----------------------------
// Linter API (policy layer)
// -----------------------------

use super::checker::{self, SmilesModels};
use super::config::{SmilesCheckFlags, SmilesParseFlags, SmilesLintConfig};
use super::checker::valence::ValenceModel;
use super::checker::aromaticity::AromaticityModel;
use super::ir::Molecule;
use super::diagnostics::{Code};

#[derive(Clone, Default)]
pub struct LinterConfig {
    pub parse_flags: SmilesParseFlags,
    pub check_flags: SmilesCheckFlags,
    pub profile: &'static str,
    pub enabled: Vec<Code>,
    pub disabled: Vec<Code>,
    pub severity_overrides: Vec<(Code, Severity)>,
    pub suppress_codes: Vec<Code>,
    pub suppress_categories: Vec<Category>,
    pub suppress_spans: Vec<(usize, usize)>,
}

#[derive(Clone)]
pub struct LintModels {
    pub valence: ValenceModel,
    pub aromaticity: AromaticityModel,
}

#[derive(Clone, Default)]
pub struct LintStats {
    pub num_errors: usize,
    pub num_warnings: usize,
    pub num_by_category: Vec<(Category, usize)>,
}

#[derive(Default)]
pub struct LintOutput {
    pub diagnostics: DiagnosticList,
    pub annotations: checker::Annotations,
    pub stats: LintStats,
}

pub struct Linter {
    pub config: LinterConfig,
    pub models: LintModels,
}

impl Linter {
    pub fn new(config: LinterConfig, models: LintModels) -> Self { Self { config, models } }

    pub fn run(&self, input: &[u8]) -> LintOutput {
        let parse_res = parse_smiles_core(input, &self.config.parse_flags);
        match parse_res {
            Ok(parsed) => self.run_parsed(&parsed.sir),
            Err(e) => LintOutput { diagnostics: DiagnosticList::from(e.as_diagnostic("")), annotations: checker::Annotations::default(), stats: LintStats::default() },
        }
    }

    pub fn run_parsed(&self, mol: &Molecule) -> LintOutput {
        // Run checkers
        let mut diagnostics = DiagnosticList::new();
        let mut annotations = checker::Annotations::default();

        checker::topology::check_topology(mol, &self.config.check_flags, &mut annotations, &mut diagnostics);
        let v_cfg = checker::valence::ValenceConfig::default();
        checker::valence::check_valence(
            mol,
            &self.config.check_flags,
            &mut annotations,
            &mut diagnostics,
            &self.models.valence,
            &v_cfg,
        );
        checker::stereo_double::check_stereo_double(mol, &self.config.check_flags, &mut annotations, &mut diagnostics);
        checker::stereo_chiral::check_stereo_chiral(mol, &self.config.check_flags, &mut annotations, &mut diagnostics);
        let a_cfg = checker::aromaticity::AromaticityConfig::default();
        checker::aromaticity::check_aromaticity(mol, &self.config.check_flags, &mut annotations, &mut diagnostics, &self.models.aromaticity, &a_cfg);

        // Apply policy
        self.apply_policy(&mut diagnostics);

        // Stats
        let stats = self.compute_stats(&diagnostics);
        LintOutput { diagnostics, annotations, stats }
    }

    fn apply_policy(&self, diagnostics: &mut DiagnosticList) {
        // Filter by suppressions
        diagnostics.diagnostics.retain(|d| {
            if self.config.suppress_codes.iter().any(|c| *c == d.code) { return false; }
            if self.config.suppress_categories.iter().any(|c| *c == d.category) { return false; }
            if self.config.suppress_spans.iter().any(|(s, e)| d.span.start >= *s && d.span.end <= *e) { return false; }
            if self.config.disabled.iter().any(|c| *c == d.code) { return false; }
            true
        });
        // Severity overrides
        for d in diagnostics.diagnostics.iter_mut() {
            if let Some((_, sev)) = self.config.severity_overrides.iter().find(|(c, _)| *c == d.code) {
                d.severity = *sev;
            }
        }
        // STRICT profile or flag upgrades warnings to errors
        let strict_active = self.config.check_flags.contains(SmilesCheckFlags::STRICT) || self.config.profile.eq_ignore_ascii_case("strict");
        if strict_active {
            for d in diagnostics.diagnostics.iter_mut() {
                if d.severity == Severity::Warning { d.severity = Severity::Error; }
            }
        }
        // Sorting (stable): by span.start, then category, then code
        diagnostics.diagnostics.sort_by(|a, b| {
            a.span.start.cmp(&b.span.start)
                .then(a.category.as_str().cmp(b.category.as_str()))
                .then(a.code.as_str().cmp(b.code.as_str()))
        });
    }

    fn compute_stats(&self, diagnostics: &DiagnosticList) -> LintStats {
        let mut stats = LintStats::default();
        use std::collections::HashMap;
        let mut by_cat: HashMap<Category, usize> = HashMap::new();
        for d in diagnostics.iter() {
            match d.severity {
                Severity::Error => stats.num_errors += 1,
                Severity::Warning => stats.num_warnings += 1,
            }
            *by_cat.entry(d.category).or_insert(0) += 1;
        }
        stats.num_by_category = by_cat.into_iter().collect();
        stats
    }
}
