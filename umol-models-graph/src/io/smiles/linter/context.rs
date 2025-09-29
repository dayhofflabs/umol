//! Context for SMILES linting.

pub struct LintContext<'a> { pub input: &'a str }

impl<'a> LintContext<'a> { pub fn new(input: &'a str) -> Self { Self { input } } }
