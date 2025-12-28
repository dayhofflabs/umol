//! Edit definitions for GraphIR.

use std::fmt;

use strum::{Display, EnumDiscriminants, EnumIter};

#[derive(Clone, Debug, Display, PartialEq, EnumDiscriminants, EnumIter)]
#[strum(serialize_all = "SCREAMING_SNAKE_CASE")]
pub enum Edit {
    #[strum(message = "No operation")]
    NoOp,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct EditList {
    pub edits: Vec<Edit>,
}

impl EditList {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, edit: Edit) {
        self.edits.push(edit);
    }

    pub fn extend<I: IntoIterator<Item = Edit>>(&mut self, edits: I) {
        self.edits.extend(edits);
    }

    pub fn append_list(&mut self, other: &mut EditList) {
        self.edits.append(&mut other.edits);
    }

    pub fn iter(&self) -> impl Iterator<Item = &Edit> {
        self.edits.iter()
    }

    pub fn into_vec(self) -> Vec<Edit> {
        self.edits
    }
}

impl IntoIterator for EditList {
    type Item = Edit;
    type IntoIter = std::vec::IntoIter<Edit>;

    fn into_iter(self) -> Self::IntoIter {
        self.edits.into_iter()
    }
}

impl fmt::Display for EditList {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for edit in &self.edits {
            writeln!(f, "- {}", edit)?;
        }
        Ok(())
    }
}
