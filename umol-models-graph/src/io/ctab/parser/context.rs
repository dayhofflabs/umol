//! Parser context

use std::collections::HashMap;

use umol::error::DataError;
use umol::Result;

use crate::io::ctab::molecule::Molecule;

#[derive(Debug, Clone, PartialEq)]
pub struct Context {
    // SGroup data, keyed by SGroup index and data field name
    sgroup_data: HashMap<(usize, String), Vec<String>>,

    // Accumulate warnings
    // warnings: Vec<String>,
}

impl Context {
    pub fn new() -> Self {
        Self {
            sgroup_data: HashMap::new(),
        }
    }

    pub fn finalize(self, _molecule: &mut Molecule) -> Result<()> {
        // Validate that no incomplete SGroups are left
        // TODO: Return list of incomplete SGroups
        if !self.sgroup_data.is_empty() {
            return Err(DataError::InvalidFeature(
                "SGroup data is incomplete".to_string(),
            )
            .into());
        }
        Ok(())
    }
}
