//! Parsing configuration

#[derive(Debug, Clone)]
pub struct MolParsingConfig {
    pub allow_unicode: bool,
    pub allow_truncated_lines: bool,
    pub allow_queries: bool,
    pub allow_superatoms: bool,
    pub allow_named_isotopes: bool,
}

impl Default for MolParsingConfig {
    fn default() -> Self {
        Self {
            allow_unicode: false,
            allow_truncated_lines: true,
            allow_queries: true,
            allow_superatoms: true,
            allow_named_isotopes: true,
        }
    }
}