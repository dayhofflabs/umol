//! Parser context

use std::collections::HashMap;

use crate::table_ir::SGroupType;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Context {
    pub current_sgroup_index: Option<usize>,
    pub current_data_field: Option<String>,
    pub current_data_content: Option<Vec<String>>,
    pub sgroup_types: HashMap<usize, SGroupType>,
}

impl Context {
    pub fn new() -> Self {
        Self::default()
    }
}
