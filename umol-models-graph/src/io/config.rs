//! Parsing configuration

use bitflags::bitflags;

bitflags! {
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct ParseFlags: u32 {
        // Core chemical features (bits 0-7)
        const NAMED_ISOTOPES = 1;         // D, T recognition
        const QUERIES = 2;                // RBC, SUB, UNS, LIN, ALS
        const EXTENDED_QUERIES = 4;       // AH, QH, XH, MH (CXSMILES)
        const ELECTRONS = 8;              // LP (lone pairs)
        const RGROUPS = 16;               // APO, AAL, RGP, LOG
        const SGROUPS = 32;               // STY, SST, SLB, SAL, SBL, SMT, SCN
        const ADVANCED_SGROUPS = 64;      // SDS, SPA, CRS, SDI, SBV, SDT, SDD, SCD, SED, SPL, SNC
        const CLARK_EXTENSIONS = 128;     // ZBO, ZCH, HYD
        const LEGACY_FEATURES = 256;      // Legacy atom list

        // Reserved for extensions (bits 9-15)
        // const RESERVED_1 = 512;
        // const RESERVED_2 = 1024;
        // const RESERVED_3 = 2048;
        // const RESERVED_4 = 4096;
        // const RESERVED_5 = 8192;
        // const RESERVED_6 = 16384;
        // const RESERVED_7 = 32768;

        // Non-standard ergonomic features (bits 16-23)
        const UNICODE = 65536;            // Unicode whitespace handling
        const STRICT_PADDING = 131072;    // Extra field validation

        // Presets
        const MINIMAL = Self::NAMED_ISOTOPES.bits();
        const BASIC = Self::MINIMAL .bits() | Self::SGROUPS.bits() | Self::CLARK_EXTENSIONS.bits();
        const EXTENDED = Self::BASIC.bits() | Self::QUERIES.bits() | Self::RGROUPS.bits() | Self::ELECTRONS.bits() | Self::LEGACY_FEATURES.bits();
        const FULL = Self::EXTENDED.bits() | Self::ADVANCED_SGROUPS.bits();
        const ALL = Self::FULL.bits() | Self::EXTENDED_QUERIES.bits();
        const STRICT = Self::MINIMAL.bits() | Self::STRICT_PADDING.bits();
        const LENIENT = Self::ALL.bits() | Self::UNICODE.bits();
    }
}

#[derive(Debug, Clone)]
pub struct ParsingConfig {
    pub parse_flags: ParseFlags,
}

impl ParsingConfig {
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
}

impl Default for ParsingConfig {
    fn default() -> Self {
        Self::basic()
    }
}
