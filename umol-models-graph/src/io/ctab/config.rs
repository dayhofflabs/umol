//! Parsing configuration for CTab-based formats

use std::fmt;

use bitflags::bitflags;

bitflags! {
    /// Flags for parsing CTab-based formats
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct CtabParseFlags: u32 {
        // Core chemical features (bits 0-7)
        const NAMED_ISOTOPES = 1;         // D, T recognition
        const PSEUDOATOMS = 2;            // General pseudoatoms (Ala)
        const QUERIES = 4;                // RBC, SUB, UNS, LIN, ALS
        const EXTENDED_QUERIES = 8;       // AH, QH, XH, MH (CXSMILES)
        const ELECTRONS = 16;             // LP (lone pairs)
        const RGROUPS = 32;               // APO, AAL, RGP, LOG
        const SGROUPS = 64;               // STY, SST, SLB, SAL, SBL, SMT, SCN
        const ADVANCED_SGROUPS = 128;     // SDS, SPA, CRS, SDI, SBV, SDT, SDD, SCD, SED, SPL, SNC
        const EXTENDED_RANGE = 256;       // Extended range of values (bond orders etc.)
        const EXTENDED_ISOTOPES = 512;    // Extended isotopes (no catalog check)
        const CLARK_EXTENSIONS = 1024;    // ZBO, ZCH, HYD
        const LEGACY_FEATURES = 2048;     // Legacy atom list, missing V2000 tag

        // Reserved for extensions (bits 10-15)
        // const RESERVED_1 = 4096;
        // const RESERVED_2 = 8192;
        // const RESERVED_3 = 16384;
        // const RESERVED_4 = 32768;

        // Ergonomic features (bits 16-23)
        const UNICODE = 65536;            // Unicode whitespace handling
        const STRICT_PADDING = 131072;    // Extra field validation
        const DEBUG = 262144;             // Debug output during parsing

        // Presets
        const MINIMAL = 0;
        const BASIC = Self::MINIMAL .bits() | Self::NAMED_ISOTOPES.bits() | Self::SGROUPS.bits() | Self::CLARK_EXTENSIONS.bits();
        const EXTENDED = (Self::BASIC.bits() | Self::QUERIES.bits() | Self::RGROUPS.bits() | Self::ELECTRONS.bits() |
                         Self::PSEUDOATOMS.bits() | Self::EXTENDED_RANGE.bits() | Self::EXTENDED_ISOTOPES.bits());
        const FULL = Self::EXTENDED.bits() | Self::ADVANCED_SGROUPS.bits() | Self::EXTENDED_QUERIES.bits();
        const STRICT = Self::MINIMAL.bits() | Self::STRICT_PADDING.bits();
        const LENIENT = Self::FULL.bits() | Self::UNICODE.bits();
    }
}

impl fmt::Display for CtabParseFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut parts = Vec::new();

        if *self == CtabParseFlags::MINIMAL {
            parts.push("MINIMAL");
        } else if *self == CtabParseFlags::BASIC {
            parts.push("BASIC");
        } else if *self == CtabParseFlags::EXTENDED {
            parts.push("EXTENDED");
        } else if *self == CtabParseFlags::FULL {
            parts.push("FULL");
        } else {
            // Show individual flags
            if self.contains(CtabParseFlags::NAMED_ISOTOPES) {
                parts.push("NAMED_ISOTOPES");
            }
            if self.contains(CtabParseFlags::RGROUPS) {
                parts.push("RGROUPS");
            }
            if self.contains(CtabParseFlags::QUERIES) {
                parts.push("QUERIES");
            }
            if self.contains(CtabParseFlags::EXTENDED_QUERIES) {
                parts.push("EXTENDED_QUERIES");
            }
            if self.contains(CtabParseFlags::ELECTRONS) {
                parts.push("ELECTRONS");
            }
            if self.contains(CtabParseFlags::PSEUDOATOMS) {
                parts.push("PSEUDOATOMS");
            }
            if self.contains(CtabParseFlags::UNICODE) {
                parts.push("UNICODE");
            }
            if self.contains(CtabParseFlags::STRICT_PADDING) {
                parts.push("STRICT_PADDING");
            }
            if self.contains(CtabParseFlags::LEGACY_FEATURES) {
                parts.push("LEGACY_FEATURES");
            }
        }

        write!(f, "{}", parts.join(" | "))
    }
}

impl Default for CtabParseFlags {
    fn default() -> Self {
        Self::MINIMAL
    }
}

/// Configuration for MOL file parsing/writing
#[derive(Debug, Clone, Default)]
pub struct MolIoConfig {
    pub parse_flags: CtabParseFlags,
}

impl MolIoConfig {
    pub fn with_parse_flags(flags: CtabParseFlags) -> Self {
        Self { parse_flags: flags }
    }

    pub fn minimal() -> Self {
        Self::with_parse_flags(CtabParseFlags::MINIMAL)
    }

    pub fn basic() -> Self {
        Self::with_parse_flags(CtabParseFlags::BASIC)
    }

    pub fn extended() -> Self {
        Self::with_parse_flags(CtabParseFlags::EXTENDED)
    }

    pub fn full() -> Self {
        Self::with_parse_flags(CtabParseFlags::FULL)
    }

    pub fn strict() -> Self {
        Self::with_parse_flags(CtabParseFlags::STRICT)
    }

    pub fn lenient() -> Self {
        Self::with_parse_flags(CtabParseFlags::LENIENT)
    }
}

impl fmt::Display for MolIoConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MolIoConfig(CtabParseFlags: {})", self.parse_flags)
    }
}
