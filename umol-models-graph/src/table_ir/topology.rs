//! Topology types for TableIR.

use crate::span::Span;

/// Ring type for TableIR
#[derive(Clone, Debug, PartialEq)]
pub struct Ring {
    pub ring_idx: u32,
    pub start_atom: Option<u32>,
    pub end_atom: Option<u32>,
    pub open_span: Option<Span>,
    pub close_span: Option<Span>,
}

impl Ring {
    pub fn is_closed(&self) -> bool {
        self.start_atom.is_some() && self.end_atom.is_some()
    }

    pub fn update_atoms(&self, a: u32, b: u32) -> Self {
        let mut updated = self.clone();
        updated.start_atom = self.start_atom.map(|_| a);
        updated.end_atom = self.end_atom.map(|_| b);
        updated
    }
}
