//! Python bindings for configured SMILES ingestion.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_io::smiles::config::{
    SmilesIoConfig as IoSmilesIoConfig, SmilesSyntaxFlags as IoSmilesSyntaxFlags,
};

/// Parser capabilities and named SMILES acceptance-policy presets.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmilesParseFlags {
    flags: IoSmilesSyntaxFlags,
}

#[pymethods]
impl SmilesParseFlags {
    #[new]
    fn new(bits: u32) -> PyResult<Self> {
        IoSmilesSyntaxFlags::from_bits(bits)
            .map(Self::from_rust)
            .ok_or_else(|| PyValueError::new_err(format!("unknown SMILES parse flag bits: {bits}")))
    }

    #[classattr]
    #[pyo3(name = "EXTENDED_AROMATICS")]
    fn extended_aromatics() -> Self {
        Self::from_rust(IoSmilesSyntaxFlags::EXTENDED_AROMATICS)
    }

    #[classattr]
    #[pyo3(name = "EXTENDED_BONDS")]
    fn extended_bonds() -> Self {
        Self::from_rust(IoSmilesSyntaxFlags::EXTENDED_BONDS)
    }

    #[classattr]
    #[pyo3(name = "CHEMAXON_EXTENSIONS")]
    fn chemaxon_extensions() -> Self {
        Self::from_rust(IoSmilesSyntaxFlags::CHEMAXON_EXTENSIONS)
    }

    #[classattr]
    #[pyo3(name = "SKIP_UNKNOWN_CHEMAXON_TAGS")]
    fn skip_unknown_chemaxon_tags() -> Self {
        Self::from_rust(IoSmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS)
    }

    #[classattr]
    #[pyo3(name = "OPENSMILES")]
    fn opensmiles() -> Self {
        Self::from_rust(IoSmilesSyntaxFlags::OPENSMILES)
    }

    #[classattr]
    #[pyo3(name = "LENIENT")]
    fn lenient() -> Self {
        Self::from_rust(IoSmilesSyntaxFlags::LENIENT)
    }

    #[classattr]
    #[pyo3(name = "CHEMAXON")]
    fn chemaxon() -> Self {
        Self::from_rust(IoSmilesSyntaxFlags::CHEMAXON)
    }

    #[getter]
    fn bits(&self) -> u32 {
        self.flags.bits()
    }

    fn __or__(&self, other: &Self) -> Self {
        Self::from_rust(self.flags | other.flags)
    }

    fn __repr__(&self) -> String {
        let name = self.flags.to_string();
        if self.flags == IoSmilesSyntaxFlags::LENIENT || self.flags.bits().count_ones() <= 1 {
            format!("SmilesParseFlags.{name}")
        } else {
            format!("SmilesParseFlags({name})")
        }
    }
}

impl SmilesParseFlags {
    pub(crate) fn from_rust(flags: IoSmilesSyntaxFlags) -> Self {
        Self { flags }
    }

    pub(crate) fn to_rust(self) -> IoSmilesSyntaxFlags {
        self.flags
    }
}

/// Owned configuration for SMILES parser acceptance policy.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmilesIoConfig {
    parse_flags: SmilesParseFlags,
}

#[pymethods]
impl SmilesIoConfig {
    #[staticmethod]
    fn opensmiles() -> Self {
        Self::from_rust(&IoSmilesIoConfig::opensmiles())
    }

    #[staticmethod]
    fn lenient() -> Self {
        Self::from_rust(&IoSmilesIoConfig::lenient())
    }

    #[staticmethod]
    fn chemaxon() -> Self {
        Self::from_rust(&IoSmilesIoConfig::chemaxon())
    }

    #[staticmethod]
    fn with_parse_flags(parse_flags: SmilesParseFlags) -> Self {
        Self { parse_flags }
    }

    #[getter]
    fn parse_flags(&self) -> SmilesParseFlags {
        self.parse_flags
    }

    fn __repr__(&self) -> String {
        if self.parse_flags == SmilesParseFlags::opensmiles() {
            "SmilesIoConfig.opensmiles()".to_owned()
        } else if self.parse_flags == SmilesParseFlags::lenient() {
            "SmilesIoConfig.lenient()".to_owned()
        } else if self.parse_flags == SmilesParseFlags::chemaxon() {
            "SmilesIoConfig.chemaxon()".to_owned()
        } else {
            format!(
                "SmilesIoConfig.with_parse_flags({})",
                self.parse_flags.__repr__()
            )
        }
    }
}

impl SmilesIoConfig {
    pub(crate) fn from_rust(config: &IoSmilesIoConfig) -> Self {
        Self {
            parse_flags: SmilesParseFlags::from_rust(config.syntax_flags),
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for configured SMILES ingestion"
    )]
    pub(crate) fn to_rust(self) -> IoSmilesIoConfig {
        IoSmilesIoConfig::with_syntax_flags(self.parse_flags.to_rust())
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_io::smiles::config::SmilesLintFlags as IoSmilesLintFlags;

    use super::*;

    #[rstest]
    #[case::opensmiles(0, IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::extended_aromatics(1 << 1, IoSmilesSyntaxFlags::EXTENDED_AROMATICS)]
    #[case::extended_bonds(1 << 2, IoSmilesSyntaxFlags::EXTENDED_BONDS)]
    #[case::chemaxon_extensions(1 << 3, IoSmilesSyntaxFlags::CHEMAXON_EXTENSIONS)]
    #[case::skip_unknown_chemaxon_tags(
        1 << 10,
        IoSmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS
    )]
    #[case::lenient(6, IoSmilesSyntaxFlags::LENIENT)]
    fn test_smiles_parse_flags_new(#[case] bits: u32, #[case] expected: IoSmilesSyntaxFlags) {
        assert_eq!(SmilesParseFlags::new(bits).unwrap().to_rust(), expected);
    }

    #[rstest]
    #[case::retired_wildcard(1, "unknown SMILES parse flag bits: 1")]
    #[case::unknown_high_bit(1 << 31, "unknown SMILES parse flag bits: 2147483648")]
    fn test_smiles_parse_flags_new_error(#[case] bits: u32, #[case] expected: &str) {
        Python::attach(|py| {
            let error = SmilesParseFlags::new(bits).unwrap_err();
            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::extended_aromatics(
        SmilesParseFlags::extended_aromatics(),
        IoSmilesSyntaxFlags::EXTENDED_AROMATICS
    )]
    #[case::extended_bonds(
        SmilesParseFlags::extended_bonds(),
        IoSmilesSyntaxFlags::EXTENDED_BONDS
    )]
    #[case::chemaxon_extensions(
        SmilesParseFlags::chemaxon_extensions(),
        IoSmilesSyntaxFlags::CHEMAXON_EXTENSIONS
    )]
    #[case::skip_unknown_chemaxon_tags(
        SmilesParseFlags::skip_unknown_chemaxon_tags(),
        IoSmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS
    )]
    #[case::opensmiles(SmilesParseFlags::opensmiles(), IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::lenient(SmilesParseFlags::lenient(), IoSmilesSyntaxFlags::LENIENT)]
    #[case::chemaxon(SmilesParseFlags::chemaxon(), IoSmilesSyntaxFlags::CHEMAXON)]
    fn test_smiles_parse_flags_classattrs(
        #[case] flags: SmilesParseFlags,
        #[case] expected: IoSmilesSyntaxFlags,
    ) {
        assert_eq!(flags.to_rust(), expected);
    }

    #[rstest]
    #[case::opensmiles(SmilesParseFlags::opensmiles(), 0)]
    #[case::extended_bonds(SmilesParseFlags::extended_bonds(), 1 << 2)]
    #[case::lenient(SmilesParseFlags::lenient(), 6)]
    fn test_smiles_parse_flags_bits(#[case] flags: SmilesParseFlags, #[case] expected: u32) {
        assert_eq!(flags.bits(), expected);
    }

    #[rstest]
    #[case::extended_syntax(
        SmilesParseFlags::extended_aromatics(),
        SmilesParseFlags::extended_bonds(),
        IoSmilesSyntaxFlags::EXTENDED_AROMATICS | IoSmilesSyntaxFlags::EXTENDED_BONDS
    )]
    #[case::chemaxon_lenient_tags(
        SmilesParseFlags::chemaxon(),
        SmilesParseFlags::skip_unknown_chemaxon_tags(),
        IoSmilesSyntaxFlags::CHEMAXON | IoSmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS
    )]
    fn test_smiles_parse_flags_or(
        #[case] left: SmilesParseFlags,
        #[case] right: SmilesParseFlags,
        #[case] expected: IoSmilesSyntaxFlags,
    ) {
        assert_eq!(left.__or__(&right).to_rust(), expected);
    }

    #[rstest]
    #[case::opensmiles(SmilesParseFlags::opensmiles(), "SmilesParseFlags.OPENSMILES")]
    #[case::extended_aromatics(
        SmilesParseFlags::extended_aromatics(),
        "SmilesParseFlags.EXTENDED_AROMATICS"
    )]
    #[case::chemaxon_alias(SmilesParseFlags::chemaxon_extensions(), "SmilesParseFlags.CHEMAXON")]
    #[case::lenient(SmilesParseFlags::lenient(), "SmilesParseFlags.LENIENT")]
    fn test_smiles_parse_flags_repr(#[case] flags: SmilesParseFlags, #[case] expected: &str) {
        assert_eq!(flags.__repr__(), expected);
    }

    #[rstest]
    #[case::opensmiles(IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::extended_aromatics(IoSmilesSyntaxFlags::EXTENDED_AROMATICS)]
    #[case::extended_bonds(IoSmilesSyntaxFlags::EXTENDED_BONDS)]
    #[case::chemaxon_extensions(IoSmilesSyntaxFlags::CHEMAXON_EXTENSIONS)]
    #[case::skip_unknown_chemaxon_tags(IoSmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS)]
    #[case::lenient(IoSmilesSyntaxFlags::LENIENT)]
    fn test_smiles_parse_flags_from_rust(#[case] flags: IoSmilesSyntaxFlags) {
        assert_eq!(SmilesParseFlags::from_rust(flags).flags, flags);
    }

    #[rstest]
    #[case::opensmiles(IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::extended_aromatics(IoSmilesSyntaxFlags::EXTENDED_AROMATICS)]
    #[case::extended_bonds(IoSmilesSyntaxFlags::EXTENDED_BONDS)]
    #[case::chemaxon_extensions(IoSmilesSyntaxFlags::CHEMAXON_EXTENSIONS)]
    #[case::skip_unknown_chemaxon_tags(IoSmilesSyntaxFlags::SKIP_UNKNOWN_CHEMAXON_TAGS)]
    #[case::lenient(IoSmilesSyntaxFlags::LENIENT)]
    fn test_smiles_parse_flags_to_rust(#[case] expected: IoSmilesSyntaxFlags) {
        assert_eq!(SmilesParseFlags::from_rust(expected).to_rust(), expected);
    }

    #[rstest]
    fn test_smiles_io_config_opensmiles() {
        assert_eq!(
            SmilesIoConfig::opensmiles(),
            SmilesIoConfig {
                parse_flags: SmilesParseFlags::from_rust(IoSmilesSyntaxFlags::OPENSMILES),
            }
        );
    }

    #[rstest]
    fn test_smiles_io_config_lenient() {
        assert_eq!(
            SmilesIoConfig::lenient(),
            SmilesIoConfig {
                parse_flags: SmilesParseFlags::from_rust(IoSmilesSyntaxFlags::LENIENT),
            }
        );
    }

    #[rstest]
    fn test_smiles_io_config_chemaxon() {
        assert_eq!(
            SmilesIoConfig::chemaxon(),
            SmilesIoConfig {
                parse_flags: SmilesParseFlags::from_rust(IoSmilesSyntaxFlags::CHEMAXON),
            }
        );
    }

    #[rstest]
    #[case::opensmiles(SmilesParseFlags::opensmiles())]
    #[case::extended_syntax(
        SmilesParseFlags::extended_aromatics().__or__(&SmilesParseFlags::extended_bonds())
    )]
    fn test_smiles_io_config_with_parse_flags(#[case] parse_flags: SmilesParseFlags) {
        assert_eq!(
            SmilesIoConfig::with_parse_flags(parse_flags),
            SmilesIoConfig { parse_flags }
        );
    }

    #[rstest]
    #[case::opensmiles(SmilesIoConfig::opensmiles(), SmilesParseFlags::opensmiles())]
    #[case::lenient(SmilesIoConfig::lenient(), SmilesParseFlags::lenient())]
    #[case::chemaxon(SmilesIoConfig::chemaxon(), SmilesParseFlags::chemaxon())]
    fn test_smiles_io_config_parse_flags(
        #[case] config: SmilesIoConfig,
        #[case] expected: SmilesParseFlags,
    ) {
        assert_eq!(config.parse_flags(), expected);
    }

    #[rstest]
    #[case::opensmiles(SmilesIoConfig::opensmiles(), "SmilesIoConfig.opensmiles()")]
    #[case::lenient(SmilesIoConfig::lenient(), "SmilesIoConfig.lenient()")]
    #[case::chemaxon(SmilesIoConfig::chemaxon(), "SmilesIoConfig.chemaxon()")]
    fn test_smiles_io_config_repr(#[case] config: SmilesIoConfig, #[case] expected: &str) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::opensmiles(IoSmilesIoConfig::opensmiles(), IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::lenient(IoSmilesIoConfig::lenient(), IoSmilesSyntaxFlags::LENIENT)]
    #[case::chemaxon(IoSmilesIoConfig::chemaxon(), IoSmilesSyntaxFlags::CHEMAXON)]
    #[case::extended_syntax(
        IoSmilesIoConfig::with_syntax_flags(
            IoSmilesSyntaxFlags::EXTENDED_AROMATICS | IoSmilesSyntaxFlags::EXTENDED_BONDS
        ),
        IoSmilesSyntaxFlags::EXTENDED_AROMATICS | IoSmilesSyntaxFlags::EXTENDED_BONDS
    )]
    fn test_smiles_io_config_from_rust(
        #[case] config: IoSmilesIoConfig,
        #[case] expected: IoSmilesSyntaxFlags,
    ) {
        assert_eq!(
            SmilesIoConfig::from_rust(&config).parse_flags.to_rust(),
            expected
        );
    }

    #[rstest]
    #[case::opensmiles(SmilesIoConfig::opensmiles(), IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::lenient(SmilesIoConfig::lenient(), IoSmilesSyntaxFlags::LENIENT)]
    #[case::chemaxon(SmilesIoConfig::chemaxon(), IoSmilesSyntaxFlags::CHEMAXON)]
    #[case::extended_syntax(
        SmilesIoConfig::with_parse_flags(
            SmilesParseFlags::extended_aromatics().__or__(&SmilesParseFlags::extended_bonds())
        ),
        IoSmilesSyntaxFlags::EXTENDED_AROMATICS | IoSmilesSyntaxFlags::EXTENDED_BONDS
    )]
    fn test_smiles_io_config_to_rust(
        #[case] config: SmilesIoConfig,
        #[case] expected: IoSmilesSyntaxFlags,
    ) {
        let config = config.to_rust();
        assert_eq!(config.syntax_flags, expected);
        assert_eq!(config.lint_flags, IoSmilesLintFlags::ALL);
        assert!(config.lint_config.enabled.is_empty());
        assert!(config.lint_config.disabled.is_empty());
        assert!(!config.lint_config.enable_gir);
    }
}
