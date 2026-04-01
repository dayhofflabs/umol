//! Parser configuration.

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

/// Parser configuration.
#[derive(Clone, Debug)]
pub struct ParseConfig {
    pub dialect: Dialect,
    pub duplicate_keys: DuplicateKeyPolicy,
}

impl Default for ParseConfig {
    fn default() -> Self {
        ParseConfig {
            dialect: Dialect::default(),
            duplicate_keys: DuplicateKeyPolicy::default(),
        }
    }
}
