//! Parser configuration.

use std::collections::HashMap;

use crate::edn::Edn;
use crate::error::EdnError;

/// Which dialect to parse.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Dialect {
    /// Strict EDN spec: no `\b`, `\f`, octal escapes, `\formfeed`, `\backspace`.
    Edn,
    /// Clojure-compatible extensions: `\b`, `\f`, octal string escapes,
    /// `\formfeed` and `\backspace` character literals.
    #[default]
    Clojure,
}

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
    readers: HashMap<String, TagFn>,
}

impl TagReaders {
    /// Create a registry with no readers.
    pub fn empty() -> Self {
        Self {
            readers: HashMap::new(),
        }
    }

    /// Register a tag reader. The tag should not include `#`.
    pub fn insert(&mut self, tag: impl Into<String>, f: TagFn) {
        self.readers.insert(tag.into(), f);
    }

    /// Look up a reader for the given tag.
    pub fn get(&self, tag: &str) -> Option<&TagFn> {
        self.readers.get(tag)
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

/// Namespace resolution for `::keyword` and `::alias/name` syntax (Clojure only).
///
/// In Clojure, `::foo` resolves to `:current-ns/foo` and `::str/foo` resolves
/// to `:fully.qualified.ns/foo` where `str` is looked up in `aliases`.
#[derive(Clone, Debug, Default)]
pub struct AutoResolve {
    pub current_ns: String,
    pub aliases: HashMap<String, String>,
}

/// Parser configuration.
#[derive(Clone, Debug)]
pub struct ParseConfig {
    pub dialect: Dialect,
    pub duplicate_keys: DuplicateKeyPolicy,
    pub tag_readers: TagReaders,
    pub auto_resolve: Option<AutoResolve>,
}

impl Default for ParseConfig {
    fn default() -> Self {
        ParseConfig {
            dialect: Dialect::default(),
            duplicate_keys: DuplicateKeyPolicy::default(),
            tag_readers: TagReaders::default(),
            auto_resolve: None,
        }
    }
}
