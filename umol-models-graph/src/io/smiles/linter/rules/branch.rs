//! Branch structure rules (pre-parse) using Segments.

use super::{Rule, RuleMeta};
use crate::diagnostics::{Category, Code, Severity, Span};
use crate::io::smiles::linter::emitter::{DiagnosticCandidate, Emitter, Scope};
use crate::io::smiles::linter::LintContext;
use crate::io::smiles::iterators::Segment;

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

        let mark_content = |saw_stack: &mut Vec<bool>| {
            if let Some(top) = saw_stack.last_mut() { *top = true; }
        };

        let next_non_ws = |segs: &Vec<Segment>, mut i: usize| -> Option<(usize, Span)> {
            while i < segs.len() {
                match &segs[i] {
                    Segment::WhitespaceBlock { .. } => i += 1,
                    Segment::BranchClose { span }
                    | Segment::ComponentSeparator { span } => return Some((i, *span)),
                    _ => return None,
                }
            }
            None
        };

        let mut i = 0usize;
        while i < segs.len() {
            match segs[i] {
                Segment::WhitespaceBlock { .. } => {}
                Segment::BranchOpen { span } => {
                    depth += 1;
                    saw_stack.push(false);
                    open_stack.push(span);
                }
                Segment::BranchClose { span } => {
                    if depth == 0 {
                        emit.candidate(DiagnosticCandidate {
                            code: Code("BRCH_UNEXPECTED_CLOSE"),
                            category: Category::Branch,
                            severity: Severity::Error,
                            span,
                            message: "Unmatched ')' with no open branch",
                            scope: Scope::Global,
                        });
                    } else {
                        let saw = saw_stack.pop().unwrap_or(false);
                        let _open = open_stack.pop().unwrap_or(span);
                        if !saw {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("BRCH_EMPTY_BRANCH"),
                                category: Category::Branch,
                                severity: Severity::Error,
                                span,
                                message: "Empty branch",
                                scope: Scope::Global,
                            });
                        }
                        depth -= 1;
                    }
                }
                Segment::ComponentSeparator { span } => {
                    if depth > 0 {
                        // Empty component since last sep/open
                        if !saw_stack.last().copied().unwrap_or(true) {
                            emit.candidate(DiagnosticCandidate {
                                code: Code("BRCH_EMPTY_BRANCH"),
                                category: Category::Branch,
                                severity: Severity::Error,
                                span,
                                message: "Empty branch",
                                scope: Scope::Global,
                            });
                        }
                        // Branch cannot cross components
                        emit.candidate(DiagnosticCandidate {
                            code: Code("BRCH_UNCLOSED"),
                            category: Category::Branch,
                            severity: Severity::Error,
                            span,
                            message: "Open branch not closed before component separator",
                            scope: Scope::Global,
                        });
                        // Reset component-local flag
                        if let Some(top) = saw_stack.last_mut() { *top = false; }
                    }
                }
                Segment::Bond { span, .. } => {
                    if depth > 0 {
                        if let Some((j, nspan)) = next_non_ws(&segs, i + 1) {
                            match segs[j] {
                                Segment::BranchClose { .. } | Segment::ComponentSeparator { .. } => {
                                    emit.candidate(DiagnosticCandidate {
                                        code: Code("BRCH_DANGLING_BOND"),
                                        category: Category::Branch,
                                        severity: Severity::Error,
                                        span: Span::new(span.start, nspan.end),
                                        message: "Dangling bond inside branch",
                                        scope: Scope::Global,
                                    });
                                }
                                _ => {}
                            }
                        }
                    }
                }
                Segment::AtomSimple { .. } | Segment::AtomBracket { .. } | Segment::RingClosure { .. } => {
                    // Content inside current subcomponent
                    if depth > 0 { mark_content(&mut saw_stack); }
                }
                _ => {}
            }
            i += 1;
        }

        // EOF: any remaining open branches
        if depth > 0 {
            // Use the earliest unclosed open to the end
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


