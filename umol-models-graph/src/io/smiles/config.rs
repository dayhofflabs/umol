//! SMILES parsing configuration

use std::fmt;

use bitflags::bitflags;

bitflags! {
    /// Flags for parsing SMILES strings
    #[derive(Debug, Clone, Copy, PartialEq, Default)]
    pub struct SmilesParseFlags: u32 {
        const STRICT_OPENSMILES = 0;
    }
}

impl fmt::Display for SmilesParseFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "STRICT_OPENSMILES")
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
#[derive(Debug, Clone, Default)]
pub struct SmilesIoConfig {
    pub parse_flags: SmilesParseFlags,
    pub lint_flags: SmilesLintFlags,
    pub lint_config: SmilesLintConfig,
}

impl SmilesIoConfig {
    pub fn strict_opensmiles() -> Self {
        Self::default()
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
