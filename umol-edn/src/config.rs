//! Parser configuration.

use crate::edn::Edn;
use crate::error::EdnError;

/// Behavior on duplicate map keys.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DuplicateKeyPolicy {
    #[default]
    Error,
    LastWins,
}

/// A tag reader transforms the parsed value after a tag into an `Edn` value.
pub type TagFn = fn(Edn) -> Result<Edn, EdnError>;

/// Registry of tag readers.
///
/// When the parser encounters `#tag value`, it looks up `tag` in this registry.
/// If found, the handler is called with the parsed value. If not found, the
/// value is wrapped as `Edn::Tagged(tag, value)`.
#[derive(Clone, Debug)]
pub struct TagReaders {
    readers: Vec<(Box<str>, TagFn)>,
}

impl TagReaders {
    /// Create a registry with no readers.
    pub fn empty() -> Self {
        Self {
            readers: Vec::new(),
        }
    }

    /// Register a tag reader. The tag should not include `#`.
    pub fn insert(&mut self, tag: impl Into<String>, f: TagFn) {
        let tag: String = tag.into();
        if let Some(entry) = self.readers.iter_mut().find(|(k, _)| **k == *tag) {
            entry.1 = f;
        } else {
            self.readers.push((tag.into_boxed_str(), f));
        }
    }

    /// Look up a reader for the given tag.
    pub fn get(&self, tag: &str) -> Option<&TagFn> {
        self.readers.iter().find(|(k, _)| &**k == tag).map(|(_, f)| f)
    }
}

impl Default for TagReaders {
    #[allow(unused_mut)]
    fn default() -> Self {
        let mut r = Self::empty();
        #[cfg(feature = "chrono")]
        r.insert("inst", crate::tags::read_inst);
        #[cfg(feature = "uuid")]
        r.insert("uuid", crate::tags::read_uuid);
        r
    }
}

/// Parser configuration.
#[derive(Clone, Debug, Default)]
pub struct ParseConfig {
    pub duplicate_keys: DuplicateKeyPolicy,
    pub tag_readers: TagReaders,
}
