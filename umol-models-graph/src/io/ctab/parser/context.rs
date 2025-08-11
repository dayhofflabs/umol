//! Parser context

use crate::io::ctab::sgroup::SGroupType;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Context {
    pub current_sgroup_index: Option<usize>,
    pub current_data_field: Option<String>,
    pub current_data_content: Option<Vec<String>>,
    pub current_sgroup_type: Option<SGroupType>,
}

impl Context {
    pub fn new() -> Self {
        Self::default()
    }
}
