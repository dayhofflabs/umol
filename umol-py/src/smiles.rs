//! Python bindings for configured SMILES ingestion.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_io::smiles::config::SmilesParseFlags as IoSmilesParseFlags;

/// Parser capabilities and named SMILES acceptance-policy presets.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmilesParseFlags {
    flags: IoSmilesParseFlags,
}

#[pymethods]
impl SmilesParseFlags {
    #[new]
    fn new(bits: u32) -> PyResult<Self> {
        IoSmilesParseFlags::from_bits(bits)
            .map(Self::from_rust)
            .ok_or_else(|| PyValueError::new_err(format!("unknown SMILES parse flag bits: {bits}")))
    }

    #[classattr]
    #[pyo3(name = "EXTENDED_AROMATICS")]
    fn extended_aromatics() -> Self {
        Self::from_rust(IoSmilesParseFlags::EXTENDED_AROMATICS)
    }

    #[classattr]
    #[pyo3(name = "EXTENDED_BONDS")]
    fn extended_bonds() -> Self {
        Self::from_rust(IoSmilesParseFlags::EXTENDED_BONDS)
    }

    #[classattr]
    #[pyo3(name = "CHEMAXON_EXTENSIONS")]
    fn chemaxon_extensions() -> Self {
        Self::from_rust(IoSmilesParseFlags::CHEMAXON_EXTENSIONS)
    }

    #[classattr]
    #[pyo3(name = "SKIP_UNKNOWN_CHEMAXON_TAGS")]
    fn skip_unknown_chemaxon_tags() -> Self {
        Self::from_rust(IoSmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS)
    }

    #[classattr]
    #[pyo3(name = "OPENSMILES")]
    fn opensmiles() -> Self {
        Self::from_rust(IoSmilesParseFlags::OPENSMILES)
    }

    #[classattr]
    #[pyo3(name = "LENIENT")]
    fn lenient() -> Self {
        Self::from_rust(IoSmilesParseFlags::LENIENT)
    }

    #[classattr]
    #[pyo3(name = "CHEMAXON")]
    fn chemaxon() -> Self {
        Self::from_rust(IoSmilesParseFlags::CHEMAXON)
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
        if self.flags == IoSmilesParseFlags::LENIENT || self.flags.bits().count_ones() <= 1 {
            format!("SmilesParseFlags.{name}")
        } else {
            format!("SmilesParseFlags({name})")
        }
    }
}

impl SmilesParseFlags {
    pub(crate) fn from_rust(flags: IoSmilesParseFlags) -> Self {
        Self { flags }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for SMILES configuration values"
    )]
    pub(crate) fn to_rust(self) -> IoSmilesParseFlags {
        self.flags
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::opensmiles(0, IoSmilesParseFlags::OPENSMILES)]
    #[case::extended_aromatics(1 << 1, IoSmilesParseFlags::EXTENDED_AROMATICS)]
    #[case::extended_bonds(1 << 2, IoSmilesParseFlags::EXTENDED_BONDS)]
    #[case::chemaxon_extensions(1 << 3, IoSmilesParseFlags::CHEMAXON_EXTENSIONS)]
    #[case::skip_unknown_chemaxon_tags(
        1 << 10,
        IoSmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS
    )]
    #[case::combination(
        (1 << 1) | (1 << 2),
        IoSmilesParseFlags::EXTENDED_AROMATICS | IoSmilesParseFlags::EXTENDED_BONDS
    )]
    #[case::lenient(1038, IoSmilesParseFlags::LENIENT)]
    fn test_smiles_parse_flags_new(#[case] bits: u32, #[case] expected: IoSmilesParseFlags) {
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
        IoSmilesParseFlags::EXTENDED_AROMATICS
    )]
    #[case::extended_bonds(SmilesParseFlags::extended_bonds(), IoSmilesParseFlags::EXTENDED_BONDS)]
    #[case::chemaxon_extensions(
        SmilesParseFlags::chemaxon_extensions(),
        IoSmilesParseFlags::CHEMAXON_EXTENSIONS
    )]
    #[case::skip_unknown_chemaxon_tags(
        SmilesParseFlags::skip_unknown_chemaxon_tags(),
        IoSmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS
    )]
    #[case::opensmiles(SmilesParseFlags::opensmiles(), IoSmilesParseFlags::OPENSMILES)]
    #[case::lenient(SmilesParseFlags::lenient(), IoSmilesParseFlags::LENIENT)]
    #[case::chemaxon(SmilesParseFlags::chemaxon(), IoSmilesParseFlags::CHEMAXON)]
    fn test_smiles_parse_flags_classattrs(
        #[case] flags: SmilesParseFlags,
        #[case] expected: IoSmilesParseFlags,
    ) {
        assert_eq!(flags.to_rust(), expected);
    }

    #[rstest]
    #[case::opensmiles(SmilesParseFlags::opensmiles(), 0)]
    #[case::extended_bonds(SmilesParseFlags::extended_bonds(), 1 << 2)]
    #[case::lenient(SmilesParseFlags::lenient(), 1038)]
    fn test_smiles_parse_flags_bits(#[case] flags: SmilesParseFlags, #[case] expected: u32) {
        assert_eq!(flags.bits(), expected);
    }

    #[rstest]
    #[case::extended_syntax(
        SmilesParseFlags::extended_aromatics(),
        SmilesParseFlags::extended_bonds(),
        IoSmilesParseFlags::EXTENDED_AROMATICS | IoSmilesParseFlags::EXTENDED_BONDS
    )]
    #[case::chemaxon_lenient_tags(
        SmilesParseFlags::chemaxon(),
        SmilesParseFlags::skip_unknown_chemaxon_tags(),
        IoSmilesParseFlags::CHEMAXON | IoSmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS
    )]
    fn test_smiles_parse_flags_or(
        #[case] left: SmilesParseFlags,
        #[case] right: SmilesParseFlags,
        #[case] expected: IoSmilesParseFlags,
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
    #[case::combination(
        SmilesParseFlags::extended_aromatics().__or__(&SmilesParseFlags::extended_bonds()),
        "SmilesParseFlags(EXTENDED_AROMATICS | EXTENDED_BONDS)"
    )]
    #[case::lenient(SmilesParseFlags::lenient(), "SmilesParseFlags.LENIENT")]
    fn test_smiles_parse_flags_repr(#[case] flags: SmilesParseFlags, #[case] expected: &str) {
        assert_eq!(flags.__repr__(), expected);
    }

    #[rstest]
    #[case::opensmiles(IoSmilesParseFlags::OPENSMILES)]
    #[case::extended_aromatics(IoSmilesParseFlags::EXTENDED_AROMATICS)]
    #[case::extended_bonds(IoSmilesParseFlags::EXTENDED_BONDS)]
    #[case::chemaxon_extensions(IoSmilesParseFlags::CHEMAXON_EXTENSIONS)]
    #[case::skip_unknown_chemaxon_tags(IoSmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS)]
    #[case::lenient(IoSmilesParseFlags::LENIENT)]
    fn test_smiles_parse_flags_from_rust(#[case] flags: IoSmilesParseFlags) {
        assert_eq!(SmilesParseFlags::from_rust(flags).flags, flags);
    }

    #[rstest]
    #[case::opensmiles(IoSmilesParseFlags::OPENSMILES)]
    #[case::extended_aromatics(IoSmilesParseFlags::EXTENDED_AROMATICS)]
    #[case::extended_bonds(IoSmilesParseFlags::EXTENDED_BONDS)]
    #[case::chemaxon_extensions(IoSmilesParseFlags::CHEMAXON_EXTENSIONS)]
    #[case::skip_unknown_chemaxon_tags(IoSmilesParseFlags::SKIP_UNKNOWN_CHEMAXON_TAGS)]
    #[case::lenient(IoSmilesParseFlags::LENIENT)]
    fn test_smiles_parse_flags_to_rust(#[case] expected: IoSmilesParseFlags) {
        assert_eq!(SmilesParseFlags::from_rust(expected).to_rust(), expected);
    }
}
