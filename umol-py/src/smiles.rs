//! Python bindings for configured SMILES ingestion.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_io::smiles::config::{
    SmilesIoConfig as IoSmilesIoConfig, SmilesSyntaxFlags as IoSmilesSyntaxFlags,
};

/// Ordinary SMILES syntax capabilities and named presets.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmilesSyntaxFlags {
    flags: IoSmilesSyntaxFlags,
}

#[pymethods]
impl SmilesSyntaxFlags {
    #[new]
    fn new(bits: u32) -> PyResult<Self> {
        IoSmilesSyntaxFlags::from_bits(bits)
            .filter(|flags| flags.bits() & !IoSmilesSyntaxFlags::LENIENT.bits() == 0)
            .map(Self::from_rust)
            .ok_or_else(|| {
                PyValueError::new_err(format!("unknown SMILES syntax flag bits: {bits}"))
            })
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
    #[pyo3(name = "OPENSMILES")]
    fn opensmiles() -> Self {
        Self::from_rust(IoSmilesSyntaxFlags::OPENSMILES)
    }

    #[classattr]
    #[pyo3(name = "LENIENT")]
    fn lenient() -> Self {
        Self::from_rust(IoSmilesSyntaxFlags::LENIENT)
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
            format!("SmilesSyntaxFlags.{name}")
        } else {
            format!("SmilesSyntaxFlags({name})")
        }
    }
}

impl SmilesSyntaxFlags {
    pub(crate) fn from_rust(flags: IoSmilesSyntaxFlags) -> Self {
        Self { flags }
    }

    pub(crate) fn to_rust(self) -> IoSmilesSyntaxFlags {
        self.flags
    }
}

/// Owned configuration for ordinary SMILES syntax.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmilesIoConfig {
    syntax_flags: SmilesSyntaxFlags,
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
    #[pyo3(signature = (*, syntax_flags))]
    fn with_syntax_flags(syntax_flags: SmilesSyntaxFlags) -> Self {
        Self { syntax_flags }
    }

    #[getter]
    fn syntax_flags(&self) -> SmilesSyntaxFlags {
        self.syntax_flags
    }

    fn __repr__(&self) -> String {
        if self.syntax_flags == SmilesSyntaxFlags::opensmiles() {
            "SmilesIoConfig.opensmiles()".to_owned()
        } else if self.syntax_flags == SmilesSyntaxFlags::lenient() {
            "SmilesIoConfig.lenient()".to_owned()
        } else {
            format!(
                "SmilesIoConfig.with_syntax_flags(syntax_flags={})",
                self.syntax_flags.__repr__()
            )
        }
    }
}

impl SmilesIoConfig {
    pub(crate) fn from_rust(config: &IoSmilesIoConfig) -> Self {
        Self {
            syntax_flags: SmilesSyntaxFlags::from_rust(config.syntax_flags),
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for configured SMILES ingestion"
    )]
    pub(crate) fn to_rust(self) -> IoSmilesIoConfig {
        IoSmilesIoConfig::with_syntax_flags(self.syntax_flags.to_rust())
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
    #[case::lenient(6, IoSmilesSyntaxFlags::LENIENT)]
    fn test_smiles_syntax_flags_new(#[case] bits: u32, #[case] expected: IoSmilesSyntaxFlags) {
        assert_eq!(SmilesSyntaxFlags::new(bits).unwrap().to_rust(), expected);
    }

    #[rstest]
    #[case::retired_wildcard(1, "unknown SMILES syntax flag bits: 1")]
    #[case::chemaxon_extensions(1 << 3, "unknown SMILES syntax flag bits: 8")]
    #[case::skip_unknown_chemaxon_tags(1 << 10, "unknown SMILES syntax flag bits: 1024")]
    #[case::unknown_high_bit(1 << 31, "unknown SMILES syntax flag bits: 2147483648")]
    fn test_smiles_syntax_flags_new_error(#[case] bits: u32, #[case] expected: &str) {
        Python::attach(|py| {
            let error = SmilesSyntaxFlags::new(bits).unwrap_err();
            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::extended_aromatics(
        SmilesSyntaxFlags::extended_aromatics(),
        IoSmilesSyntaxFlags::EXTENDED_AROMATICS
    )]
    #[case::extended_bonds(
        SmilesSyntaxFlags::extended_bonds(),
        IoSmilesSyntaxFlags::EXTENDED_BONDS
    )]
    #[case::opensmiles(SmilesSyntaxFlags::opensmiles(), IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::lenient(SmilesSyntaxFlags::lenient(), IoSmilesSyntaxFlags::LENIENT)]
    fn test_smiles_syntax_flags_classattrs(
        #[case] flags: SmilesSyntaxFlags,
        #[case] expected: IoSmilesSyntaxFlags,
    ) {
        assert_eq!(flags.to_rust(), expected);
    }

    #[rstest]
    #[case::opensmiles(SmilesSyntaxFlags::opensmiles(), 0)]
    #[case::extended_aromatics(SmilesSyntaxFlags::extended_aromatics(), 1 << 1)]
    #[case::extended_bonds(SmilesSyntaxFlags::extended_bonds(), 1 << 2)]
    #[case::lenient(SmilesSyntaxFlags::lenient(), 6)]
    fn test_smiles_syntax_flags_bits(#[case] flags: SmilesSyntaxFlags, #[case] expected: u32) {
        assert_eq!(flags.bits(), expected);
    }

    #[rstest]
    #[case::extended_syntax(
        SmilesSyntaxFlags::extended_aromatics(),
        SmilesSyntaxFlags::extended_bonds(),
        IoSmilesSyntaxFlags::EXTENDED_AROMATICS | IoSmilesSyntaxFlags::EXTENDED_BONDS
    )]
    fn test_smiles_syntax_flags_or(
        #[case] left: SmilesSyntaxFlags,
        #[case] right: SmilesSyntaxFlags,
        #[case] expected: IoSmilesSyntaxFlags,
    ) {
        assert_eq!(left.__or__(&right).to_rust(), expected);
    }

    #[rstest]
    #[case::opensmiles(SmilesSyntaxFlags::opensmiles(), "SmilesSyntaxFlags.OPENSMILES")]
    #[case::extended_aromatics(
        SmilesSyntaxFlags::extended_aromatics(),
        "SmilesSyntaxFlags.EXTENDED_AROMATICS"
    )]
    #[case::extended_bonds(
        SmilesSyntaxFlags::extended_bonds(),
        "SmilesSyntaxFlags.EXTENDED_BONDS"
    )]
    #[case::lenient(SmilesSyntaxFlags::lenient(), "SmilesSyntaxFlags.LENIENT")]
    fn test_smiles_syntax_flags_repr(#[case] flags: SmilesSyntaxFlags, #[case] expected: &str) {
        assert_eq!(flags.__repr__(), expected);
    }

    #[rstest]
    #[case::opensmiles(IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::extended_aromatics(IoSmilesSyntaxFlags::EXTENDED_AROMATICS)]
    #[case::extended_bonds(IoSmilesSyntaxFlags::EXTENDED_BONDS)]
    #[case::lenient(IoSmilesSyntaxFlags::LENIENT)]
    fn test_smiles_syntax_flags_from_rust(#[case] flags: IoSmilesSyntaxFlags) {
        assert_eq!(SmilesSyntaxFlags::from_rust(flags).flags, flags);
    }

    #[rstest]
    #[case::opensmiles(IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::extended_aromatics(IoSmilesSyntaxFlags::EXTENDED_AROMATICS)]
    #[case::extended_bonds(IoSmilesSyntaxFlags::EXTENDED_BONDS)]
    #[case::lenient(IoSmilesSyntaxFlags::LENIENT)]
    fn test_smiles_syntax_flags_to_rust(#[case] expected: IoSmilesSyntaxFlags) {
        assert_eq!(SmilesSyntaxFlags::from_rust(expected).to_rust(), expected);
    }

    #[rstest]
    fn test_smiles_io_config_opensmiles() {
        assert_eq!(
            SmilesIoConfig::opensmiles(),
            SmilesIoConfig {
                syntax_flags: SmilesSyntaxFlags::from_rust(IoSmilesSyntaxFlags::OPENSMILES),
            }
        );
    }

    #[rstest]
    fn test_smiles_io_config_lenient() {
        assert_eq!(
            SmilesIoConfig::lenient(),
            SmilesIoConfig {
                syntax_flags: SmilesSyntaxFlags::from_rust(IoSmilesSyntaxFlags::LENIENT),
            }
        );
    }

    #[rstest]
    #[case::opensmiles(SmilesSyntaxFlags::opensmiles())]
    #[case::composition(
        SmilesSyntaxFlags::extended_aromatics().__or__(&SmilesSyntaxFlags::extended_bonds())
    )]
    fn test_smiles_io_config_with_syntax_flags(#[case] syntax_flags: SmilesSyntaxFlags) {
        assert_eq!(
            SmilesIoConfig::with_syntax_flags(syntax_flags),
            SmilesIoConfig { syntax_flags }
        );
    }

    #[rstest]
    #[case::opensmiles(SmilesIoConfig::opensmiles(), SmilesSyntaxFlags::opensmiles())]
    #[case::lenient(SmilesIoConfig::lenient(), SmilesSyntaxFlags::lenient())]
    #[case::extended_aromatics(
        SmilesIoConfig::with_syntax_flags(SmilesSyntaxFlags::extended_aromatics()),
        SmilesSyntaxFlags::extended_aromatics()
    )]
    fn test_smiles_io_config_syntax_flags(
        #[case] config: SmilesIoConfig,
        #[case] expected: SmilesSyntaxFlags,
    ) {
        assert_eq!(config.syntax_flags(), expected);
    }

    #[rstest]
    #[case::opensmiles(SmilesIoConfig::opensmiles(), "SmilesIoConfig.opensmiles()")]
    #[case::lenient(SmilesIoConfig::lenient(), "SmilesIoConfig.lenient()")]
    #[case::extended_aromatics(
        SmilesIoConfig::with_syntax_flags(SmilesSyntaxFlags::extended_aromatics()),
        "SmilesIoConfig.with_syntax_flags(syntax_flags=SmilesSyntaxFlags.EXTENDED_AROMATICS)"
    )]
    fn test_smiles_io_config_repr(#[case] config: SmilesIoConfig, #[case] expected: &str) {
        assert_eq!(config.__repr__(), expected);
    }

    #[rstest]
    #[case::opensmiles(IoSmilesIoConfig::opensmiles(), IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::lenient(IoSmilesIoConfig::lenient(), IoSmilesSyntaxFlags::LENIENT)]
    #[case::extended_aromatics(
        IoSmilesIoConfig::with_syntax_flags(IoSmilesSyntaxFlags::EXTENDED_AROMATICS),
        IoSmilesSyntaxFlags::EXTENDED_AROMATICS
    )]
    #[case::extended_bonds(
        IoSmilesIoConfig::with_syntax_flags(IoSmilesSyntaxFlags::EXTENDED_BONDS),
        IoSmilesSyntaxFlags::EXTENDED_BONDS
    )]
    #[case::composition(
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
            SmilesIoConfig::from_rust(&config),
            SmilesIoConfig {
                syntax_flags: SmilesSyntaxFlags::from_rust(expected),
            },
        );
    }

    #[rstest]
    #[case::opensmiles(SmilesIoConfig::opensmiles(), IoSmilesSyntaxFlags::OPENSMILES)]
    #[case::lenient(SmilesIoConfig::lenient(), IoSmilesSyntaxFlags::LENIENT)]
    #[case::extended_aromatics(
        SmilesIoConfig::with_syntax_flags(SmilesSyntaxFlags::extended_aromatics()),
        IoSmilesSyntaxFlags::EXTENDED_AROMATICS
    )]
    #[case::extended_bonds(
        SmilesIoConfig::with_syntax_flags(SmilesSyntaxFlags::extended_bonds()),
        IoSmilesSyntaxFlags::EXTENDED_BONDS
    )]
    #[case::composition(
        SmilesIoConfig::with_syntax_flags(
            SmilesSyntaxFlags::extended_aromatics().__or__(&SmilesSyntaxFlags::extended_bonds())
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
