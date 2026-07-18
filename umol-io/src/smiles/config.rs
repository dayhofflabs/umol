//! SMILES parsing configuration

use std::fmt;

use bitflags::bitflags;

bitflags! {
    /// Flags for parsing SMILES strings
    #[derive(Debug, Clone, Copy, PartialEq, Default)]
    pub struct SmilesParseFlags: u32 {
        // Format extensions
        const EXTENDED_AROMATICS = 1 << 1; // se, te, as, si
        const EXTENDED_BONDS = 1 << 2;     // ->, <-, ~
        const CHEMAXON_EXTENSIONS = 1 << 3; // |...| CXSMILES extension block

        // Input validation strictness
        const SKIP_UNKNOWN_CHEMAXON_TAGS = 1 << 10; // Skip unknown ChemAxon tags

        // Presets
        const OPENSMILES = 0;

        // Lenient parser
        const LENIENT = Self::EXTENDED_AROMATICS.bits() | Self::EXTENDED_BONDS.bits() |
            Self::CHEMAXON_EXTENSIONS.bits() | Self::SKIP_UNKNOWN_CHEMAXON_TAGS.bits();

        // ChemAxon-compatible parser
        const CHEMAXON = Self::CHEMAXON_EXTENSIONS.bits();
    }
}

impl fmt::Display for SmilesParseFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts = Vec::new();
        if *self == SmilesParseFlags::OPENSMILES {
            parts.push("OPENSMILES");
        } else if *self == SmilesParseFlags::LENIENT {
            parts.push("LENIENT");
        } else if *self == SmilesParseFlags::CHEMAXON {
            parts.push("CHEMAXON");
        } else {
            if self.contains(SmilesParseFlags::EXTENDED_AROMATICS) {
                parts.push("EXTENDED_AROMATICS");
            }
            if self.contains(SmilesParseFlags::EXTENDED_BONDS) {
                parts.push("EXTENDED_BONDS");
            }
            if self.contains(SmilesParseFlags::CHEMAXON_EXTENSIONS) {
                parts.push("CHEMAXON_EXTENSIONS");
            }
            if self.contains(SmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS) {
                parts.push("SKIP_UNKNOWN_CHEMAXON_TAGS");
            }
        }
        write!(f, "{}", parts.join(" | "))
    }
}

bitflags! {
    /// Configuration for SMILES checking
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct SmilesLintFlags: u32 {
        const TOPOLOGY = 1;
        const VALENCE = 2;
        const AROMATICITY = 4;
        const STEREO = 8;
        const STRICT = 65536;

        const NONE = 0;
        const ALL = Self::TOPOLOGY.bits() | Self::VALENCE.bits() |
                    Self::AROMATICITY.bits() | Self::STEREO.bits();
    }
}

impl fmt::Display for SmilesLintFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts = Vec::new();

        if *self == SmilesLintFlags::NONE {
            parts.push("NONE");
        } else if *self == SmilesLintFlags::ALL {
            parts.push("ALL");
        } else {
            if self.contains(SmilesLintFlags::TOPOLOGY) {
                parts.push("TOPOLOGY");
            }
            if self.contains(SmilesLintFlags::VALENCE) {
                parts.push("VALENCE");
            }
            if self.contains(SmilesLintFlags::AROMATICITY) {
                parts.push("AROMATICITY");
            }
            if self.contains(SmilesLintFlags::STEREO) {
                parts.push("STEREO");
            }
        }

        write!(f, "{}", parts.join(" | "))
    }
}

impl Default for SmilesLintFlags {
    fn default() -> Self {
        Self::ALL
    }
}

/// Configuration for SMILES linting
#[derive(Debug, Clone, Default)]
pub struct SmilesLintConfig {
    pub enabled: Vec<&'static str>,
    pub disabled: Vec<&'static str>,
    pub enable_gir: bool,
}

impl fmt::Display for SmilesLintConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SmilesLintConfig(enabled: {:?}, disabled: {:?}, enable_gir: {})",
            self.enabled, self.disabled, self.enable_gir
        )
    }
}

/// Configuration for SMILES parsing/writing
#[derive(Debug, Clone)]
pub struct SmilesIoConfig {
    pub parse_flags: SmilesParseFlags,
    pub lint_flags: SmilesLintFlags,
    pub lint_config: SmilesLintConfig,
}

impl SmilesIoConfig {
    pub fn with_parse_flags(flags: SmilesParseFlags) -> Self {
        Self {
            parse_flags: flags,
            lint_flags: SmilesLintFlags::default(),
            lint_config: SmilesLintConfig::default(),
        }
    }
    pub fn opensmiles() -> Self {
        Self::with_parse_flags(SmilesParseFlags::OPENSMILES)
    }
    pub fn lenient() -> Self {
        Self::with_parse_flags(SmilesParseFlags::LENIENT)
    }
    pub fn chemaxon() -> Self {
        Self::with_parse_flags(SmilesParseFlags::CHEMAXON)
    }
    pub fn with_lint_flags(flags: SmilesLintFlags) -> Self {
        Self {
            parse_flags: SmilesParseFlags::default(),
            lint_flags: flags,
            lint_config: SmilesLintConfig::default(),
        }
    }
    pub fn with_lint_config(config: SmilesLintConfig) -> Self {
        Self {
            parse_flags: SmilesParseFlags::default(),
            lint_flags: SmilesLintFlags::default(),
            lint_config: config,
        }
    }
}

impl Default for SmilesIoConfig {
    fn default() -> Self {
        Self::opensmiles()
    }
}

impl fmt::Display for SmilesIoConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SmilesIoConfig(parse: {}, lint: {}, lint_config: {})",
            self.parse_flags, self.lint_flags, self.lint_config
        )
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::retired_wildcard(1, None)]
    #[case::extended_aromatics(1 << 1, Some(SmilesParseFlags::EXTENDED_AROMATICS))]
    #[case::extended_bonds(1 << 2, Some(SmilesParseFlags::EXTENDED_BONDS))]
    #[case::chemaxon_extensions(1 << 3, Some(SmilesParseFlags::CHEMAXON_EXTENSIONS))]
    #[case::extended_syntax(
        (1 << 1) | (1 << 2),
        Some(SmilesParseFlags::EXTENDED_AROMATICS | SmilesParseFlags::EXTENDED_BONDS)
    )]
    fn test_smiles_parse_flags_from_bits(
        #[case] bits: u32,
        #[case] expected: Option<SmilesParseFlags>,
    ) {
        assert_eq!(SmilesParseFlags::from_bits(bits), expected);
    }

    #[rstest]
    #[case::opensmiles(SmilesParseFlags::OPENSMILES, "OPENSMILES")]
    #[case::extended_aromatics(SmilesParseFlags::EXTENDED_AROMATICS, "EXTENDED_AROMATICS")]
    #[case::extended_bonds(SmilesParseFlags::EXTENDED_BONDS, "EXTENDED_BONDS")]
    #[case::chemaxon(SmilesParseFlags::CHEMAXON, "CHEMAXON")]
    #[case::skip_unknown_chemaxon_tags(
        SmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS,
        "SKIP_UNKNOWN_CHEMAXON_TAGS"
    )]
    #[case::extended_syntax(
        SmilesParseFlags::EXTENDED_AROMATICS | SmilesParseFlags::EXTENDED_BONDS,
        "EXTENDED_AROMATICS | EXTENDED_BONDS"
    )]
    #[case::chemaxon_skip_unknown(
        SmilesParseFlags::CHEMAXON | SmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS,
        "CHEMAXON_EXTENSIONS | SKIP_UNKNOWN_CHEMAXON_TAGS"
    )]
    #[case::lenient(SmilesParseFlags::LENIENT, "LENIENT")]
    fn test_smiles_parse_flags_display(#[case] flags: SmilesParseFlags, #[case] expected: &str) {
        assert_eq!(flags.to_string(), expected);
    }

    #[rstest]
    #[case::opensmiles(SmilesIoConfig::opensmiles(), SmilesParseFlags::OPENSMILES)]
    #[case::lenient(SmilesIoConfig::lenient(), SmilesParseFlags::LENIENT)]
    #[case::chemaxon(SmilesIoConfig::chemaxon(), SmilesParseFlags::CHEMAXON)]
    fn test_smiles_io_config(#[case] config: SmilesIoConfig, #[case] expected: SmilesParseFlags) {
        assert_eq!(config.parse_flags, expected);
        assert_eq!(config.lint_flags, SmilesLintFlags::ALL);
        assert!(config.lint_config.enabled.is_empty());
        assert!(config.lint_config.disabled.is_empty());
        assert!(!config.lint_config.enable_gir);
    }

    #[rstest]
    fn test_smiles_io_config_default() {
        let config = SmilesIoConfig::default();

        assert_eq!(config.parse_flags, SmilesParseFlags::OPENSMILES);
        assert_eq!(config.lint_flags, SmilesLintFlags::ALL);
        assert!(config.lint_config.enabled.is_empty());
        assert!(config.lint_config.disabled.is_empty());
        assert!(!config.lint_config.enable_gir);
    }
}
