//! Surface encoding for handles in standalone edit documents.

use umol_edn::{DeError, Edn, EdnMapHelper, FromEdn, ToEdn};

use super::edn_utils::single_key_map;
use crate::ast::edit::{
    AromaticSystemHandle, AtomHandle, BondHandle, DativeBondHandle, MulticenterBondHandle,
    NoncovalentBondHandle, StereoAtomHandle, StereoBondHandle,
};
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

/// Surface form shared by every typed handle in a standalone edit document.
///
/// A bare integer identifies an entity in the initial host. `{:new n}` identifies the `n`th
/// same-kind entity created earlier in the edit sequence.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub(crate) enum EditHandleDsl {
    Id(u32),
    New(usize),
}

impl<'de> FromEdn<'de> for EditHandleDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Int(_) => u32::from_edn(edn).map(Self::Id),
            Edn::Map(map) => {
                let mut helper = EdnMapHelper::new(map);
                let index = helper.required("new")?;
                helper.finalize()?;
                if map.len() != 1 {
                    return Err(DeError::Custom(
                        "edit handle map keys must be keywords".to_string(),
                    ));
                }
                Ok(Self::New(index))
            }
            other => Err(DeError::TypeMismatch {
                expected: "edit handle (non-negative integer or {:new n} map)",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for EditHandleDsl {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::Id(index) => index.to_edn(),
            Self::New(index) => single_key_map("new", index.to_edn()),
        }
    }
}

macro_rules! impl_typed_handle_conversion {
    ($handle:ident, $id:ident) => {
        impl From<$handle> for EditHandleDsl {
            fn from(handle: $handle) -> Self {
                match handle {
                    $handle::Id(id) => Self::Id(id.0),
                    $handle::New(index) => Self::New(index),
                }
            }
        }

        impl From<EditHandleDsl> for $handle {
            fn from(handle: EditHandleDsl) -> Self {
                match handle {
                    EditHandleDsl::Id(index) => Self::Id($id(index)),
                    EditHandleDsl::New(index) => Self::New(index),
                }
            }
        }

        impl<'de> FromEdn<'de> for $handle {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                EditHandleDsl::from_edn(edn).map(Self::from)
            }
        }

        impl ToEdn for $handle {
            fn to_edn(&self) -> Edn<'static> {
                EditHandleDsl::from(self.clone()).to_edn()
            }
        }
    };
}

impl_typed_handle_conversion!(AtomHandle, AtomId);
impl_typed_handle_conversion!(BondHandle, BondId);
impl_typed_handle_conversion!(DativeBondHandle, DativeBondId);
impl_typed_handle_conversion!(AromaticSystemHandle, AromaticSystemId);
impl_typed_handle_conversion!(MulticenterBondHandle, MulticenterBondId);
impl_typed_handle_conversion!(NoncovalentBondHandle, NoncovalentBondId);
impl_typed_handle_conversion!(StereoAtomHandle, StereoAtomId);
impl_typed_handle_conversion!(StereoBondHandle, StereoBondId);

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::{read_string, EdnError};

    use super::*;

    #[rstest]
    #[case::id("7", EditHandleDsl::Id(7))]
    #[case::new("{:new 7}", EditHandleDsl::New(7))]
    fn test_edit_handle_dsl_edn_roundtrip(#[case] input: &str, #[case] expected: EditHandleDsl) {
        let parsed = EditHandleDsl::from_edn_str(input).unwrap();

        assert_eq!(parsed, expected);
        assert_eq!(parsed.to_edn(), read_string(input).unwrap());
    }

    #[rstest]
    #[case::negative_id(
        "-1",
        EdnError::De(DeError::OutOfRange {
            value: "-1".to_string(),
            target: "u32",
            path: Vec::new(),
        }),
    )]
    #[case::keyword(
        ":carbon",
        EdnError::De(DeError::TypeMismatch {
            expected: "edit handle (non-negative integer or {:new n} map)",
            got: "keyword",
            path: Vec::new(),
        }),
    )]
    #[case::structural_ref(
        "{:atoms [0 1]}",
        EdnError::De(DeError::MissingField {
            key: "new".to_string(),
            path: Vec::new(),
        }),
    )]
    #[case::empty_map(
        "{}",
        EdnError::De(DeError::MissingField {
            key: "new".to_string(),
            path: Vec::new(),
        }),
    )]
    #[case::negative_new(
        "{:new -1}",
        EdnError::De(DeError::OutOfRange {
            value: "-1".to_string(),
            target: "usize",
            path: Vec::new(),
        }),
    )]
    #[case::non_integer_new(
        "{:new :first}",
        EdnError::De(DeError::TypeMismatch {
            expected: "int",
            got: "keyword",
            path: Vec::new(),
        }),
    )]
    #[case::extra_keyword(
        "{:new 0 :atoms [0 1]}",
        EdnError::De(DeError::UnknownField {
            key: "atoms".to_string(),
            path: Vec::new(),
        }),
    )]
    #[case::extra_non_keyword(
        "{:new 0 \"atoms\" [0 1]}",
        EdnError::De(DeError::Custom(
            "edit handle map keys must be keywords".to_string(),
        )),
    )]
    fn test_edit_handle_dsl_from_edn_error(#[case] input: &str, #[case] expected: EdnError) {
        assert_eq!(EditHandleDsl::from_edn_str(input), Err(expected));
    }

    #[rstest]
    #[case::id(EditHandleDsl::Id(7))]
    #[case::new(EditHandleDsl::New(7))]
    fn test_typed_edit_handle_edn_roundtrip(#[case] handle: EditHandleDsl) {
        match handle {
            EditHandleDsl::Id(index) => {
                assert_eq!(
                    AtomHandle::from_edn_str("7"),
                    Ok(AtomHandle::Id(AtomId(index)))
                );
                assert_eq!(
                    AtomHandle::Id(AtomId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    BondHandle::from_edn_str("7"),
                    Ok(BondHandle::Id(BondId(index)))
                );
                assert_eq!(
                    BondHandle::Id(BondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    DativeBondHandle::from_edn_str("7"),
                    Ok(DativeBondHandle::Id(DativeBondId(index)))
                );
                assert_eq!(
                    DativeBondHandle::Id(DativeBondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    AromaticSystemHandle::from_edn_str("7"),
                    Ok(AromaticSystemHandle::Id(AromaticSystemId(index)))
                );
                assert_eq!(
                    AromaticSystemHandle::Id(AromaticSystemId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    MulticenterBondHandle::from_edn_str("7"),
                    Ok(MulticenterBondHandle::Id(MulticenterBondId(index)))
                );
                assert_eq!(
                    MulticenterBondHandle::Id(MulticenterBondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    NoncovalentBondHandle::from_edn_str("7"),
                    Ok(NoncovalentBondHandle::Id(NoncovalentBondId(index)))
                );
                assert_eq!(
                    NoncovalentBondHandle::Id(NoncovalentBondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    StereoAtomHandle::from_edn_str("7"),
                    Ok(StereoAtomHandle::Id(StereoAtomId(index)))
                );
                assert_eq!(
                    StereoAtomHandle::Id(StereoAtomId(index)).to_edn(),
                    Edn::Int(index.into())
                );
                assert_eq!(
                    StereoBondHandle::from_edn_str("7"),
                    Ok(StereoBondHandle::Id(StereoBondId(index)))
                );
                assert_eq!(
                    StereoBondHandle::Id(StereoBondId(index)).to_edn(),
                    Edn::Int(index.into())
                );
            }
            EditHandleDsl::New(index) => {
                assert_eq!(
                    AtomHandle::from_edn_str("{:new 7}"),
                    Ok(AtomHandle::New(index))
                );
                assert_eq!(
                    AtomHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    BondHandle::from_edn_str("{:new 7}"),
                    Ok(BondHandle::New(index))
                );
                assert_eq!(
                    BondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    DativeBondHandle::from_edn_str("{:new 7}"),
                    Ok(DativeBondHandle::New(index))
                );
                assert_eq!(
                    DativeBondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    AromaticSystemHandle::from_edn_str("{:new 7}"),
                    Ok(AromaticSystemHandle::New(index))
                );
                assert_eq!(
                    AromaticSystemHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    MulticenterBondHandle::from_edn_str("{:new 7}"),
                    Ok(MulticenterBondHandle::New(index))
                );
                assert_eq!(
                    MulticenterBondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    NoncovalentBondHandle::from_edn_str("{:new 7}"),
                    Ok(NoncovalentBondHandle::New(index))
                );
                assert_eq!(
                    NoncovalentBondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    StereoAtomHandle::from_edn_str("{:new 7}"),
                    Ok(StereoAtomHandle::New(index))
                );
                assert_eq!(
                    StereoAtomHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
                assert_eq!(
                    StereoBondHandle::from_edn_str("{:new 7}"),
                    Ok(StereoBondHandle::New(index))
                );
                assert_eq!(
                    StereoBondHandle::New(index).to_edn(),
                    read_string("{:new 7}").unwrap()
                );
            }
        }
    }
}
