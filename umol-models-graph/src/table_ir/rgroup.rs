//! RGroup for TableIR.
//! RGroup is a temporary type inherited from CTFile.
// TODO: Replace by semantically defined structures that RGroup type combines.

#[derive(Debug, Clone, PartialEq)]
pub struct RGroup {
    pub label: Option<u32>,
    pub dependent_label: Option<u32>,
    pub rgroup_or_h: bool,
    pub occurrence: Vec<RGroupOccurrence>,
}

impl RGroup {
    pub fn new(label: Option<u32>) -> Self {
        Self {
            label,
            dependent_label: None,
            rgroup_or_h: false,
            occurrence: vec![RGroupOccurrence::GreaterThan(0)],
        }
    }
}

/// RGroup occurrence type
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RGroupOccurrence {
    Exactly(u8),
    Range(u8, u8),
    GreaterThan(u8),
    FewerThan(u8),
}

impl RGroupOccurrence {
    pub fn contains(&self, count: u8) -> bool {
        match self {
            RGroupOccurrence::Exactly(n) => *n == count,
            RGroupOccurrence::Range(n, m) => count >= *n && count <= *m,
            RGroupOccurrence::GreaterThan(n) => count > *n,
            RGroupOccurrence::FewerThan(n) => count < *n,
        }
    }
}
