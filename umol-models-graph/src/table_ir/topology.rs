//! Topology types for TableIR.

use super::atom::Atom;
use super::bond::Bond;
use crate::span::Span;

/// Ring type for TableIR
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ring {
    pub ring_idx: u32,
    pub start_atom: Option<u32>,
    pub end_atom: Option<u32>,
    pub open_span: Option<Span>,
    pub close_span: Option<Span>,
}

/// Fragment type for TableIR
#[derive(Clone, Debug, PartialEq)]
pub struct Fragment {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

/// Link type for TableIR
#[derive(Clone, Debug, PartialEq)]
pub struct Link {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}
