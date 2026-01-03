//! SMILES parsing configuration

use std::fmt;

use bitflags::bitflags;

bitflags! {
    /// Flags for parsing SMILES strings
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct SmilesParseFlags: u32 {
        // Extensions (lex/syntax only)
        const EXTENDED_WS = 1;       // allow ASCII inter-token whitespace
        const ALLOWS_COMMENTS = 1 << 1;   // // line and /* block */ comments
        const EXPLICIT_EOI = 1 << 2;      // explicit end-of-input marker token

        // Reserved for extensions (bits 3-15)
        // const RESERVED_1 = 1 << 3;
        // const RESERVED_2 = 1 << 4;
        // const RESERVED_3 = 1 << 5;
        // const RESERVED_4 = 1 << 6;
        // const RESERVED_5 = 1 << 7;
        // const RESERVED_6 = 1 << 8;
        // const RESERVED_7 = 1 << 9;
        // const RESERVED_8 = 1 << 10;
        // const RESERVED_9 = 1 << 11;
        // const RESERVED_10 = 1 << 12;
        // const RESERVED_11 = 1 << 13;
        // const RESERVED_12 = 1 << 14;
        // const RESERVED_13 = 1 << 15;

        // Metadata generation
        const NO_METADATA = 1 << 16;     // no metadata emitted

        // Dialects
        // Core OpenSMILES behavior: terminator-only WS, no comments
        const STRICT_OPENSMILES = 0;

        // Umol dialect: allow inter-token whitespace and comments
        const UMOL_DIALECT = Self::EXTENDED_WS.bits() | Self::ALLOWS_COMMENTS.bits();

        // Presets
        const LENIENT = Self::UMOL_DIALECT.bits();
    }
}

impl Default for SmilesParseFlags {
    fn default() -> Self {
        Self::STRICT_OPENSMILES
    }
}

impl fmt::Display for SmilesParseFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts = Vec::new();

        if *self == SmilesParseFlags::STRICT_OPENSMILES {
            parts.push("STRICT_OPENSMILES");
        } else if *self == SmilesParseFlags::UMOL_DIALECT {
            parts.push("UMOL_DIALECT");
        } else if *self == SmilesParseFlags::LENIENT {
            parts.push("LENIENT");
        } else {
            if self.contains(SmilesParseFlags::EXTENDED_WS) {
                parts.push("EXTENDED_WS");
            }
            if self.contains(SmilesParseFlags::ALLOWS_COMMENTS) {
                parts.push("COMMENTS");
            }
            if self.contains(SmilesParseFlags::EXPLICIT_EOI) {
                parts.push("EXPLICIT_EOI");
            }
            if self.contains(SmilesParseFlags::NO_METADATA) {
                parts.push("NO_METADATA");
            }
        }

        write!(f, "{}", parts.join(" | "))
    }
}

bitflags! {
    /// Configuration for SMILES checking
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct SmilesLintFlags: u32 {

        // Lint categories
        const TOPOLOGY = 1;
        const VALENCE = 2;
        const AROMATICITY = 4;
        const STEREO = 8;

        // Strict mode: convert warnings to errors
        const STRICT = 65536;

        // Presets
        const NONE = 0;
        const ALL = (Self::TOPOLOGY.bits() | Self::VALENCE.bits() | Self::AROMATICITY.bits() |
                    Self::STEREO.bits());
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

    pub fn strict_opensmiles() -> Self {
        Self::with_parse_flags(SmilesParseFlags::STRICT_OPENSMILES)
    }

    pub fn umol_dialect() -> Self {
        Self::with_parse_flags(SmilesParseFlags::UMOL_DIALECT)
    }
}

impl fmt::Display for SmilesIoConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SmilesIoConfig(SmilesParseFlags: {}, SmilesLintFlags: {}, SmilesLintConfig: {})",
            self.parse_flags, self.lint_flags, self.lint_config
        )
    }
}
