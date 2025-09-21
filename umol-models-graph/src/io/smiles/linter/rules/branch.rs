//! Branch structure rules (pre-parse) using Segments.

use super::{Rule, RuleMeta};
use crate::diagnostics::{Category, Code, Severity, Span};
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::LintContext;
use crate::io::smiles::iterators::{BranchEventKind, Branches, Segment};

pub struct BranchRule;
static META_BRANCH: RuleMeta = RuleMeta {
    id: "BRANCH_RULE",
    category: Category::Branch,
    default_severity: Severity::Error,
};

impl Rule for BranchRule {
    fn meta(&self) -> &'static RuleMeta { &META_BRANCH }

    fn check(&self, ctx: &LintContext, emit: &mut Emitter) {
        let segs = ctx.segments();
        let mut depth: usize = 0;
        let mut saw_stack: Vec<bool> = Vec::new();
        let mut open_stack: Vec<Span> = Vec::new();

        // First pass: structure events using Branches
        for ev in Branches::new(&segs).into_iter() {
            match ev.kind {
                BranchEventKind::Open => {
                    depth += 1;
                    saw_stack.push(false);
                    open_stack.push(ev.span);
                }
                BranchEventKind::Close => {
                    if depth == 0 {
                        emit.candidate(DiagnosticCandidate {
                            code: Code("BRCH_UNEXPECTED_CLOSE"),
                            category: Category::Branch,
                            severity: Severity::Error,
                            span: ev.span,
                            message: "Unmatched ')' with no open branch",
                            scope: Scope::Global,
                        });
                    } else {
                        let saw = saw_stack.pop().unwrap_or(false);
                        let _open = open_stack.pop().unwrap_or(ev.span);
                        if !saw {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("BRCH_EMPTY_BRANCH"),
                                category: Category::Branch,
                                severity: Severity::Error,
                                span: ev.span,
                                message: "Empty branch",
                                scope: Scope::Global,
                            });
                        }
                        depth -= 1;
                    }
                }
                BranchEventKind::NewComponent => {
                    // Empty component since last sep/open
                    if !saw_stack.last().copied().unwrap_or(true) {
                        emit.candidate(DiagnosticCandidate {
                            code: Code("BRCH_EMPTY_BRANCH"),
                            category: Category::Branch,
                            severity: Severity::Error,
                            span: ev.span,
                            message: "Empty branch",
                            scope: Scope::Global,
                        });
                    }
                    // Branch cannot cross components
                    emit.candidate(DiagnosticCandidate {
                        code: Code("BRCH_UNCLOSED"),
                        category: Category::Branch,
                        severity: Severity::Error,
                        span: ev.span,
                        message: "Open branch not closed before component separator",
                        scope: Scope::Global,
                    });
                    // Reset component-local flag
                    if let Some(top) = saw_stack.last_mut() { *top = false; }
                }
            }
        }

        // Second pass: detect dangling bonds inside branches with minimal scanning
        let segs_len = segs.len();
        let mut i = 0usize;
        while i < segs_len {
            match segs[i] {
                Segment::Bond { span, .. } => {
                    // Find next non-whitespace significant terminator within branch scope
                    let mut j = i + 1;
                    while j < segs_len {
                        match segs[j] {
                            Segment::WhitespaceBlock { .. } => { j += 1; }
                            Segment::BranchClose { span: nsp } | Segment::NewComponent { span: nsp } => {
                                emit.candidate(DiagnosticCandidate {
                                    code: Code("BRCH_DANGLING_BOND"),
                                    category: Category::Branch,
                                    severity: Severity::Error,
                                    span: Span::new(span.start, nsp.end),
                                    message: "Dangling bond inside branch",
                                    scope: Scope::Global,
                                });
                                break;
                            }
                            _ => { break; }
                        }
                    }
                }
                Segment::AtomSimple { .. } | Segment::AtomBracket { .. } | Segment::RingClosure { .. } => {
                    if let Some(top) = saw_stack.last_mut() { *top = true; }
                }
                _ => {}
            }
            i += 1;
        }

        // EOF: any remaining open branches
        if depth > 0 {
            if let Some(open) = open_stack.first().copied() {
                emit.candidate(DiagnosticCandidate {
                    code: Code("BRCH_UNCLOSED"),
                    category: Category::Branch,
                    severity: Severity::Error,
                    span: Span::new(open.start, ctx.input.len()),
                    message: "Unclosed branch at end of input",
                    scope: Scope::Global,
                });
            } else {
                emit.candidate(DiagnosticCandidate {
                    code: Code("BRCH_UNCLOSED"),
                    category: Category::Branch,
                    severity: Severity::Error,
                    span: Span::new(ctx.input.len(), ctx.input.len()),
                    message: "Unclosed branch at end of input",
                    scope: Scope::Global,
                });
            }
        }
    }
}

pub static BRANCH_RULE: BranchRule = BranchRule;


