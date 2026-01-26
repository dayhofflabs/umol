//! SMILES parsing configuration

use std::fmt;

use bitflags::bitflags;

bitflags! {
    /// Flags for parsing SMILES strings
    #[derive(Debug, Clone, Copy, PartialEq, Default)]
    pub struct SmilesParseFlags: u32 {
        // Parser capabilities
        const WILDCARDS = 1; // *

        // Presets
        const BASIC_OPENSMILES = 0;

        // Maximum capabilities for basic parser
        const BASIC_MAX = Self::BASIC_OPENSMILES.bits();

        // Maximum capabilities for extended parser (everything)
        const EXTENDED_MAX = Self::BASIC_MAX.bits() | Self::WILDCARDS.bits();

        // Strict parser: only additional capabilities of extended parser over basic parser according to OpenSMILES spec
        // Basic strict: BASIC & STRICT = OPENSMILES_STRICT
        const STRICT = Self::WILDCARDS.bits();

        // Default for basic parser
        const BASIC = Self::BASIC_OPENSMILES.bits();

        // Default for extended parser
        const EXTENDED = Self::BASIC.bits() | Self::STRICT.bits();

        // Lenient parser: Currently same as extended parser
        const LENIENT = Self::EXTENDED.bits();

        // OpenSMILES parser: BASIC_OPENSMILES | WILDCARDS
        const OPENSMILES = Self::BASIC_OPENSMILES.bits() | Self::WILDCARDS.bits();
    }
}

impl fmt::Display for SmilesParseFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts = Vec::new();
        if *self == SmilesParseFlags::BASIC_OPENSMILES {
            parts.push("BASIC_OPENSMILES");
        } else if *self == SmilesParseFlags::OPENSMILES {
            parts.push("OPENSMILES");
        } else if *self == SmilesParseFlags::BASIC {
            parts.push("BASIC");
        } else if *self == SmilesParseFlags::BASIC_MAX {
            parts.push("BASIC_MAX");
        } else if *self == SmilesParseFlags::STRICT {
            parts.push("STRICT");
        } else if *self == SmilesParseFlags::EXTENDED {
            parts.push("EXTENDED");
        } else if *self == SmilesParseFlags::EXTENDED_MAX {
            parts.push("EXTENDED_MAX");
        } else if *self == SmilesParseFlags::LENIENT {
            parts.push("LENIENT");
        } else {
            if self.contains(SmilesParseFlags::WILDCARDS) {
                parts.push("WILDCARDS");
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
#[derive(Debug, Clone, Default)]
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
    pub fn basic_opensmiles() -> Self {
        Self::with_parse_flags(SmilesParseFlags::BASIC_OPENSMILES)
    }
    pub fn opensmiles() -> Self {
        Self::with_parse_flags(SmilesParseFlags::OPENSMILES)
    }
    pub fn basic() -> Self {
        Self::with_parse_flags(SmilesParseFlags::BASIC)
    }
    pub fn basic_max() -> Self {
        Self::with_parse_flags(SmilesParseFlags::BASIC_MAX)
    }
    pub fn strict() -> Self {
        Self::with_parse_flags(SmilesParseFlags::STRICT)
    }
    pub fn extended() -> Self {
        Self::with_parse_flags(SmilesParseFlags::EXTENDED)
    }
    pub fn extended_max() -> Self {
        Self::with_parse_flags(SmilesParseFlags::EXTENDED_MAX)
    }
    pub fn lenient() -> Self {
        Self::with_parse_flags(SmilesParseFlags::LENIENT)
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

impl fmt::Display for SmilesIoConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SmilesIoConfig(parse: {}, lint: {}, lint_config: {})",
            self.parse_flags, self.lint_flags, self.lint_config
        )
    }
}
