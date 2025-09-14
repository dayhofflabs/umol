//! Diagnostics types for SMILES linting (errors and warnings).

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    Error,
    Warning,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Category {
    Lex,
    Syn,
    Ring,
    Num,
    Brkt,
    Stereo,
    Style,
    Internal,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Code(pub &'static str);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
}

impl Span {
    pub fn new(start: usize, end: usize) -> Self { Self { start, end } }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    pub code: Code,
    pub severity: Severity,
    pub category: Category,
    pub span: Span,
    pub message: &'static str,
    pub details: Option<String>,
}

impl Diagnostic {
    pub fn error(code: Code, category: Category, span: Span, message: &'static str) -> Self {
        Self { code, severity: Severity::Error, category, span, message, details: None }
    }
    pub fn warning(code: Code, category: Category, span: Span, message: &'static str) -> Self {
        Self { code, severity: Severity::Warning, category, span, message, details: None }
    }
    pub fn with_details(mut self, details: impl Into<String>) -> Self {
        self.details = Some(details.into());
        self
    }
}

#[derive(Debug, Default, Clone)]
pub struct DiagnosticsReport {
    pub diagnostics: Vec<Diagnostic>,
}

impl DiagnosticsReport {
    pub fn new() -> Self { Self { diagnostics: Vec::new() } }
    pub fn push(&mut self, d: Diagnostic) { self.diagnostics.push(d); }
    pub fn has_errors(&self) -> bool { self.diagnostics.iter().any(|d| d.severity == Severity::Error) }
    pub fn errors(&self) -> impl Iterator<Item = &Diagnostic> { self.diagnostics.iter().filter(|d| d.severity == Severity::Error) }
    pub fn warnings(&self) -> impl Iterator<Item = &Diagnostic> { self.diagnostics.iter().filter(|d| d.severity == Severity::Warning) }
}

impl fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} [{}:{:?}] @{}..{}: {}", self.message, self.code.0, self.category, self.span.start, self.span.end, self.details.as_deref().unwrap_or(""))
    }
}


