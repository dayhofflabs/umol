//! Generic property container for TableIR.

/// Generic property container
#[derive(Clone, Debug, PartialEq)]
pub struct Property {
    pub name: String,
    pub value: String,
}
