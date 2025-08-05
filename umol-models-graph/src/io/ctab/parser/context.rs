//! Parser context

use crate::io::ctab::molecule::Molecule;
use std::collections::HashMap;
use umol::error::DataError;
use umol::Result;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Context {
    // SGroup data, keyed by SGroup index and field name
    sgroup_data: HashMap<(usize, String), String>,
    pub current_data_field: Option<String>,
}

impl Context {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a piece of SGroup data
    pub fn add_sgroup_data(
        &mut self,
        sgroup_index: usize,
        data: String,
        is_end: bool,
    ) -> Result<()> {
        if let Some(field_name) = &self.current_data_field {
            let entry = self
                .sgroup_data
                .entry((sgroup_index, field_name.clone()))
                .or_default();
            entry.push_str(&data);
            if is_end {
                self.current_data_field = None;
            }
            Ok(())
        } else {
            Err(DataError::InvalidFragment(format!(
                "SGroup data found for SGroup {}, but no data description (SDT) was provided.",
                sgroup_index
            ))
            .into())
        }
    }

    pub fn finalize(&self, molecule: &mut Molecule) -> Result<()> {
        for ((sgroup_index, field_name), data_content) in &self.sgroup_data {
            if let Some(sgroup) = molecule
                .sgroups
                .values_mut()
                .find(|s| s.label == Some(*sgroup_index as u32))
            {
                if let Some(data) = sgroup.data.get_mut(field_name) {
                    if let Some(dc) = &mut data.data_content {
                        dc.push(data_content.clone());
                    } else {
                        data.data_content = Some(vec![data_content.clone()]);
                    }
                } else {
                    return Err(DataError::InvalidFeature(format!(
                        "SGroup data found for SGroup {} with field name '{}', but no data description (SDT) was provided.",
                        sgroup_index, field_name
                    ))
                    .into());
                }
            } else {
                return Err(DataError::InvalidFeature(format!(
                    "SGroup data found for non-existent SGroup index {}",
                    sgroup_index
                ))
                .into());
            }
        }
        Ok(())
    }
}
