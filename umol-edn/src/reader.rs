//! Public read API: read_string, read_all, Reader iterator.

use std::sync::OnceLock;

use crate::config::ParseConfig;
use crate::edn::Edn;
use crate::error::EdnError;
use crate::parser::{parse_all, parse_value, parse_value_strict};

pub(crate) fn default_config() -> &'static ParseConfig {
    static CONFIG: OnceLock<ParseConfig> = OnceLock::new();
    CONFIG.get_or_init(ParseConfig::default)
}

/// Parse a single EDN value from a string, rejecting trailing content.
pub fn read_string<'a>(input: &'a str) -> Result<Edn<'a>, EdnError> {
    parse_value_strict(input, default_config()).map_err(Into::into)
}

/// Parse a single EDN value with custom config, rejecting trailing content.
pub fn read_string_with<'a>(input: &'a str, config: &ParseConfig) -> Result<Edn<'a>, EdnError> {
    parse_value_strict(input, config).map_err(Into::into)
}

/// Parse all EDN values from a string.
pub fn read_all<'a>(input: &'a str) -> Result<Vec<Edn<'a>>, EdnError> {
    parse_all(input, default_config()).map_err(Into::into)
}

/// Parse all EDN values with custom config.
pub fn read_all_with<'a>(input: &'a str, config: &ParseConfig) -> Result<Vec<Edn<'a>>, EdnError> {
    parse_all(input, config).map_err(Into::into)
}

/// Streaming iterator over EDN values in a string.
pub struct Reader<'a> {
    remaining: &'a str,
    config: ParseConfig,
}

impl<'a> Reader<'a> {
    pub fn new(input: &'a str) -> Self {
        Reader {
            remaining: input,
            config: ParseConfig::default(),
        }
    }

    pub fn with_config(input: &'a str, config: ParseConfig) -> Self {
        Reader {
            remaining: input,
            config,
        }
    }
}

impl<'a> Iterator for Reader<'a> {
    type Item = Result<Edn<'a>, EdnError>;

    fn next(&mut self) -> Option<Self::Item> {
        // Skip whitespace/comments to check if we're at the end
        let trimmed = self
            .remaining
            .trim_start_matches([' ', '\t', '\n', '\r', ',']);
        if trimmed.is_empty() {
            return None;
        }
        self.remaining = trimmed;

        match parse_value(self.remaining, &self.config) {
            Ok((val, rest)) => {
                self.remaining = rest;
                Some(Ok(val))
            }
            Err(e) => {
                self.remaining = "";
                Some(Err(e.into()))
            }
        }
    }
}
