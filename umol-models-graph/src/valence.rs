//! Valence models for graph-based molecular models.

use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::path::Path;

use thiserror::Error;
use toml::de::Error as TomlError;
use umol_data::Element;

use crate::diagnostics::{Diagnostic, DiagnosticKind, Severity};

#[derive(Debug, Error)]
pub enum ValenceError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),
    #[error("TOML parse error: {0}")]
    Toml(#[from] TomlError),
    #[error("No matching valence pattern for {0}")]
    NoMatch(String),
    #[error("Ambiguous valence pattern match for {0}")]
    AmbiguousMatch(String),
}

impl From<ValenceError> for Diagnostic {
    fn from(error: ValenceError) -> Self {
        use DiagnosticKind::*;
        let (kind, details) = match error {
            ValenceError::Io(ref e) => (phValenceError, Some(e.to_string())),
            ValenceError::Toml(ref e) => (GraphValenceError, Some(e.to_string())),
            ValenceError::NoMatch(ref s) => (GraphNoMatch, Some(s.clone())),
            ValenceError::AmbiguousMatch(ref s) => (GraphAmbiguousMatch, Some(s.clone())),
        };
        Diagnostic {
            kind,
            category: kind.category(),
            severity: Severity::Error,
            span: None,
            details,
        }
    }
}

impl From<ValenceError> for umol::error::ParseError {
    fn from(error: ValenceError) -> Self {
        umol::error::ParseError::Format(Box::new(error))
    }
}

impl From<ValenceError> for umol::Error {
    fn from(error: ValenceError) -> Self {
        umol::Error::Parse(error.into())
    }
}

type Result<T> = std::result::Result<T, ValenceError>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ValencePolicy {
    Ignore,
    Warn,
    Error,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValenceConfig {
    pub enabled: bool,
    pub check_brackets: bool,
    pub infer_bracket_implicit: bool,
    pub match_patterns: bool,
    pub no_match_policy: ValencePolicy,
    pub out_of_range_policy: ValencePolicy,
    pub ambiguous_match_policy: ValencePolicy,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ValenceModel {
    states: HashMap<Element, Vec<u8>>,
    pub patterns: ValencePatternTable,
}

impl ValenceModel {
    pub fn simple_organic() -> Self {
        Self {
            states: HashMap::new(),
            patterns: ValencePatternTable::new(),
        }
    }
    pub fn states_for(&self, e: Element) -> Option<&[u8]> {
        self.states.get(&e).map(|v| &v[..])
    }
    pub fn set_states(&mut self, e: Element, states: Vec<u8>) {
        self.states.insert(e, states);
    }
    pub fn set_patterns(&mut self, patterns: ValencePatternTable) {
        self.patterns = patterns;
    }
    pub fn load_patterns_from_vec(&mut self, patterns: Vec<ValencePattern>) {
        self.patterns = ValencePatternTable { patterns };
    }
    pub fn with_patterns(table: ValencePatternTable) -> Self {
        Self {
            states: HashMap::new(),
            patterns: table,
        }
    }
    pub fn from_patterns_reader<R: io::Read>(mut r: R) -> Result<Self> {
        let mut buf = String::new();
        r.read_to_string(&mut buf)?;
        let table: ValencePatternTable = toml::from_str(&buf)?;
        Ok(Self::with_patterns(table))
    }
    pub fn from_patterns_file(path: &Path) -> Result<Self> {
        let f = File::open(path)?;
        Self::from_patterns_reader(f)
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ValencePattern {
    pub element: Option<Element>,
    pub bond_sum: Option<u8>,
    pub charge: Option<i8>,
    pub implicit_hydrogens: Option<u8>,
    pub unpaired: Option<u8>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ValencePatternTable {
    pub patterns: Vec<ValencePattern>,
}

impl ValencePatternTable {
    pub fn new() -> Self {
        Self::default()
    }
}
