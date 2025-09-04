//! Parsing configuration

use bitflags::bitflags;
use std::fmt;

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct ParseFlags: u32 {
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

#[derive(Debug, Clone)]
pub struct ParseConfig {
    pub parse_flags: ParseFlags,
}

impl ParseConfig {
    pub fn to_flags(&self) -> ParseFlags {
        self.parse_flags
    }

    pub fn with_flags(flags: ParseFlags) -> Self {
        Self { parse_flags: flags }
    }

    pub fn minimal() -> Self {
        Self::with_flags(ParseFlags::MINIMAL)
    }

    pub fn basic() -> Self {
        Self::with_flags(ParseFlags::BASIC)
    }

    pub fn extended() -> Self {
        Self::with_flags(ParseFlags::EXTENDED)
    }

    pub fn full() -> Self {
        Self::with_flags(ParseFlags::FULL)
    }

    pub fn strict() -> Self {
        Self::with_flags(ParseFlags::STRICT)
    }

    pub fn lenient() -> Self {
        Self::with_flags(ParseFlags::LENIENT)
    }
}

impl Default for ParseConfig {
    fn default() -> Self {
        Self::basic()
    }
}

impl fmt::Display for ParseConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ParseConfig({})", self.parse_flags)
    }
}

impl fmt::Display for ParseFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            return write!(f, "EMPTY");
        }

        let mut parts = Vec::new();

        if *self == ParseFlags::MINIMAL {
            parts.push("MINIMAL");
        } else if *self == ParseFlags::BASIC {
            parts.push("BASIC");
        } else if *self == ParseFlags::EXTENDED {
            parts.push("EXTENDED");
        } else if *self == ParseFlags::FULL {
            parts.push("FULL");
        } else {
            // Show individual flags
            if self.contains(ParseFlags::NAMED_ISOTOPES) {
                parts.push("NAMED_ISOTOPES");
            }
            if self.contains(ParseFlags::RGROUPS) {
                parts.push("RGROUPS");
            }
            if self.contains(ParseFlags::QUERIES) {
                parts.push("QUERIES");
            }
            if self.contains(ParseFlags::EXTENDED_QUERIES) {
                parts.push("EXTENDED_QUERIES");
            }
            if self.contains(ParseFlags::ELECTRONS) {
                parts.push("ELECTRONS");
            }
            if self.contains(ParseFlags::PSEUDOATOMS) {
                parts.push("PSEUDOATOMS");
            }
            if self.contains(ParseFlags::UNICODE) {
                parts.push("UNICODE");
            }
            if self.contains(ParseFlags::STRICT_PADDING) {
                parts.push("STRICT_PADDING");
            }
            if self.contains(ParseFlags::LEGACY_FEATURES) {
                parts.push("LEGACY_FEATURES");
            }
        }

        write!(f, "{}", parts.join(" | "))
    }
}
