//! Context for SMILES linting.

use std::cell::{Ref, RefCell};

use crate::io::smiles::lexer_old::Lexer;
use crate::io::smiles::iterators::{Segment, Segments};

pub struct LintContext<'a> {
    pub input: &'a str,
    pub lexer: Lexer<'a>,
    // Lazily available resources as needed later
    segments: RefCell<Option<Vec<Segment<'a>>>>,
}

impl<'a> LintContext<'a> {
    pub fn new(input: &'a str) -> Self {
        Self {
            input,
            lexer: Lexer::new(input),
            segments: RefCell::new(None),
        }
    }

    pub fn segments(&self) -> Ref<'_, Vec<Segment<'a>>> {
        if self.segments.borrow().is_none() {
            let v = Segments::new(self.input).collect::<Vec<_>>();
            *self.segments.borrow_mut() = Some(v);
        }
        Ref::map(self.segments.borrow(), |opt: &Option<Vec<Segment<'a>>>| opt.as_ref().unwrap())
    }
}
