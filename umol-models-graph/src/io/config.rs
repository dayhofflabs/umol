//! Parsing configuration

use std::fmt;

use bitflags::bitflags;

bitflags! {
    /// Flags for parsing MOL files
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct MolParseFlags: u32 {
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

impl fmt::Display for MolParseFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            return write!(f, "EMPTY");
        }

        let mut parts = Vec::new();

        if *self == MolParseFlags::MINIMAL {
            parts.push("MINIMAL");
        } else if *self == MolParseFlags::BASIC {
            parts.push("BASIC");
        } else if *self == MolParseFlags::EXTENDED {
            parts.push("EXTENDED");
        } else if *self == MolParseFlags::FULL {
            parts.push("FULL");
        } else {
            // Show individual flags
            if self.contains(MolParseFlags::NAMED_ISOTOPES) {
                parts.push("NAMED_ISOTOPES");
            }
            if self.contains(MolParseFlags::RGROUPS) {
                parts.push("RGROUPS");
            }
            if self.contains(MolParseFlags::QUERIES) {
                parts.push("QUERIES");
            }
            if self.contains(MolParseFlags::EXTENDED_QUERIES) {
                parts.push("EXTENDED_QUERIES");
            }
            if self.contains(MolParseFlags::ELECTRONS) {
                parts.push("ELECTRONS");
            }
            if self.contains(MolParseFlags::PSEUDOATOMS) {
                parts.push("PSEUDOATOMS");
            }
            if self.contains(MolParseFlags::UNICODE) {
                parts.push("UNICODE");
            }
            if self.contains(MolParseFlags::STRICT_PADDING) {
                parts.push("STRICT_PADDING");
            }
            if self.contains(MolParseFlags::LEGACY_FEATURES) {
                parts.push("LEGACY_FEATURES");
            }
        }

        write!(f, "{}", parts.join(" | "))
    }
}

/// Configuration for MOL file parsing/writing
#[derive(Debug, Clone)]
pub struct MolIoConfig {
    pub parse_flags: MolParseFlags,
}

impl MolIoConfig {
    pub fn to_flags(&self) -> MolParseFlags {
        self.parse_flags
    }

    pub fn with_flags(flags: MolParseFlags) -> Self {
        Self { parse_flags: flags }
    }

    pub fn minimal() -> Self {
        Self::with_flags(MolParseFlags::MINIMAL)
    }

    pub fn basic() -> Self {
        Self::with_flags(MolParseFlags::BASIC)
    }

    pub fn extended() -> Self {
        Self::with_flags(MolParseFlags::EXTENDED)
    }

    pub fn full() -> Self {
        Self::with_flags(MolParseFlags::FULL)
    }

    pub fn strict() -> Self {
        Self::with_flags(MolParseFlags::STRICT)
    }

    pub fn lenient() -> Self {
        Self::with_flags(MolParseFlags::LENIENT)
    }
}

impl Default for MolIoConfig {
    fn default() -> Self {
        Self::basic()
    }
}

impl fmt::Display for MolIoConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "MolIoConfig({})", self.parse_flags)
    }
}

bitflags! {
    /// Flags for parsing SMILES strings
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct SmilesParseFlags: u32 {
        // Core OpenSMILES behavior
        const STRICT_OPENSMILES = 0;   // terminator-only WS, no comments

        // Extensions (lex/syntax only)
        const INTERTOKEN_WS = 1;       // allow ASCII inter-token whitespace
        const COMMENTS = 2;            // // line and /* block */ comments
        const EXPLICIT_EOI = 4;        // explicit end-of-input marker token
        const NO_LINTS = 64;           // no lints emitted

        // Presets
        const UMOL_DIALECT = Self::INTERTOKEN_WS.bits() | Self::COMMENTS.bits();
    }
}

impl fmt::Display for SmilesParseFlags {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            return write!(f, "EMPTY");
        }

        let mut parts = Vec::new();

        if self.contains(SmilesParseFlags::STRICT_OPENSMILES) {
            parts.push("STRICT_OPENSMILES");
        }
        if self.contains(SmilesParseFlags::INTERTOKEN_WS) {
            parts.push("INTERTOKEN_WS");
        }
        if self.contains(SmilesParseFlags::COMMENTS) {
            parts.push("COMMENTS");
        }
        if self.contains(SmilesParseFlags::EXPLICIT_EOI) {
            parts.push("EXPLICIT_EOI");
        }
        if self.contains(SmilesParseFlags::NO_LINTS) {
            parts.push("NO_LINTS");
        }

        write!(f, "{}", parts.join(" | "))
    }
}

pub struct SmilesIoConfig {
    pub parse_flags: SmilesParseFlags,
    pub lint: Option<SmilesLintFlags>,
    pub check: SmilesCheckFlags,
}

impl SmilesIoConfig {
    pub fn to_flags(&self) -> SmilesParseFlags {
        self.parse_flags
    }

    pub fn with_flags(flags: SmilesParseFlags) -> Self {
        Self {
            parse_flags: flags,
            lint: None,
            check: SmilesCheckFlags::default(),
        }
    }

    pub fn strict_opensmiles() -> Self {
        Self::with_flags(SmilesParseFlags::STRICT_OPENSMILES)
    }
}

impl Default for SmilesIoConfig {
    fn default() -> Self {
        Self::strict_opensmiles()
    }
}

impl fmt::Display for SmilesIoConfig {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SmilesIoConfig({})", self.parse_flags)
    }
}

#[derive(Debug, Clone, Default)]
pub struct SmilesLintFlags {
    pub enabled: Vec<&'static str>,
    pub disabled: Vec<&'static str>,
}

#[derive(Debug, Clone, Default)]
pub struct SmilesCheckFlags {
    pub enable_topology: bool,
    pub enable_valence: bool,
    pub enable_aromaticity: bool,
    pub enable_stereo: bool,
}
