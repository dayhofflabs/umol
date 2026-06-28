//! Surface-level entity references. Each ref is a positional index (`Edn::Int`)
//! or a symbolic id keyword (`Edn::Keyword`); resolution to / from the AST id
//! consults the surrounding `MoleculeMetadata`.

use indexmap::IndexMap;
use umol_edn::{DeError, Edn, EdnError, EdnStreamDeserializer, FromEdn, ToEdn};

use super::edn_utils::eof_err;
use super::error::ParseError;
use super::molecule::MoleculeMetadata;
use crate::ast::id::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
    StereoAtomId, StereoBondId,
};

macro_rules! define_ref {
    ($name:ident, $id:ident, $accessor:ident, $kind:literal, $reader:ident) => {
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub enum $name {
            Index(usize),
            Id(String),
        }

        impl $name {
            /// Build a ref from an AST index, preferring an id from `metadata`
            /// if one is recorded for this index.
            pub fn from_ast(id: $id, metadata: &MoleculeMetadata) -> Self {
                if let Some(name) = metadata.$accessor(id) {
                    Self::Id(name.to_string())
                } else {
                    Self::Index(id.index())
                }
            }

            /// Resolve this ref to an AST index against `metadata`. Fails on
            /// unknown id or out-of-range numeric index.
            pub fn into_ast(
                self,
                count: usize,
                metadata: &MoleculeMetadata,
            ) -> Result<$id, ParseError> {
                match self {
                    Self::Index(i) => {
                        if i < count {
                            Ok($id::from(i))
                        } else {
                            Err(ParseError::InvalidRef {
                                kind: $kind,
                                value: i.to_string(),
                            })
                        }
                    }
                    Self::Id(name) => {
                        for i in 0..count {
                            let id = $id::from(i);
                            if metadata.$accessor(id) == Some(name.as_str()) {
                                return Ok(id);
                            }
                        }
                        Err(ParseError::InvalidRef {
                            kind: $kind,
                            value: name,
                        })
                    }
                }
            }

            /// Resolve this ref against a pre-built id → index map. O(1) id
            /// lookup; intended for entity-loop resolution where cloning the
            /// full `MoleculeMetadata` per call is wasteful.
            pub fn resolve(
                self,
                count: usize,
                id_to_idx: &IndexMap<String, $id>,
            ) -> Result<$id, ParseError> {
                match self {
                    Self::Index(i) => {
                        if i < count {
                            Ok($id::from(i))
                        } else {
                            Err(ParseError::InvalidRef {
                                kind: $kind,
                                value: i.to_string(),
                            })
                        }
                    }
                    Self::Id(name) => id_to_idx.get(&name).copied().ok_or(ParseError::InvalidRef {
                        kind: $kind,
                        value: name,
                    }),
                }
            }
        }

        impl<'de> FromEdn<'de> for $name {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                match edn {
                    Edn::Int(n) => {
                        let i = usize::try_from(*n).map_err(|_| DeError::OutOfRange {
                            value: n.to_string(),
                            target: "usize",
                            path: Vec::new(),
                        })?;
                        Ok(Self::Index(i))
                    }
                    Edn::Keyword(k) => Ok(Self::Id(k.name().to_string())),
                    other => Err(DeError::TypeMismatch {
                        expected: concat!($kind, " ref (int or keyword)"),
                        got: other.kind(),
                        path: Vec::new(),
                    }),
                }
            }
        }

        impl ToEdn for $name {
            fn to_edn(&self) -> Edn<'static> {
                match self {
                    Self::Index(i) => Edn::Int(*i as i64),
                    Self::Id(name) => Edn::Keyword(umol_edn::EdnKeyword::owned(name.clone())),
                }
            }
        }

        pub(super) fn $reader(de: &mut EdnStreamDeserializer<'_>) -> Result<$name, EdnError> {
            match de.peek_byte()?.ok_or_else(eof_err)? {
                b':' => Ok($name::Id(de.read_keyword_name()?.into_owned())),
                _ => {
                    let n = de.read_i64()?;
                    let i = usize::try_from(n).map_err(|_| DeError::OutOfRange {
                        value: n.to_string(),
                        target: "usize",
                        path: Vec::new(),
                    })?;
                    Ok($name::Index(i))
                }
            }
        }
    };
}

define_ref!(AtomRef, AtomId, atom_id, "atom", read_atom_ref);
define_ref!(BondRef, BondId, bond_id, "bond", read_bond_ref);
define_ref!(
    DativeBondRef,
    DativeBondId,
    dative_bond_id,
    "dative-bond",
    read_dative_bond_ref
);
define_ref!(
    AromaticSystemRef,
    AromaticSystemId,
    aromatic_system_id,
    "aromatic-system",
    read_aromatic_system_ref
);
define_ref!(
    MulticenterBondRef,
    MulticenterBondId,
    multicenter_bond_id,
    "multicenter-bond",
    read_multicenter_bond_ref
);
define_ref!(
    NoncovalentBondRef,
    NoncovalentBondId,
    noncovalent_bond_id,
    "noncovalent-bond",
    read_noncovalent_bond_ref
);
define_ref!(
    StereoAtomRef,
    StereoAtomId,
    stereo_atom_id,
    "stereo-atom",
    read_stereo_atom_ref
);
define_ref!(
    StereoBondRef,
    StereoBondId,
    stereo_bond_id,
    "stereo-bond",
    read_stereo_bond_ref
);

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::{read_string, EdnKeyword};

    use super::*;

    #[fixture]
    fn meta_with_atom_id() -> MoleculeMetadata {
        MoleculeMetadata::new().with_atom_id(AtomId(2), "c1")
    }

    #[rstest]
    #[case::int(Edn::Int(3), AtomRef::Index(3))]
    #[case::keyword(Edn::Keyword(EdnKeyword::owned("c1".into())), AtomRef::Id("c1".into()))]
    fn test_atom_ref_from_edn(#[case] input: Edn<'static>, #[case] expected: AtomRef) {
        assert_eq!(AtomRef::from_edn(&input).unwrap(), expected);
    }

    #[rstest]
    fn test_atom_ref_from_edn_rejects_other_kinds() {
        let err = AtomRef::from_edn(&Edn::Str("x".into())).unwrap_err();
        assert!(matches!(
            err,
            DeError::TypeMismatch {
                expected: "atom ref (int or keyword)",
                ..
            }
        ));
    }

    #[rstest]
    #[case::index(AtomRef::Index(5), Edn::Int(5))]
    #[case::id(AtomRef::Id("c1".into()), Edn::Keyword(EdnKeyword::owned("c1".into())))]
    fn test_atom_ref_to_edn(#[case] input: AtomRef, #[case] expected: Edn<'static>) {
        assert_eq!(input.to_edn(), expected);
    }

    #[rstest]
    #[case::int("3", AtomRef::Index(3))]
    #[case::keyword(":c1", AtomRef::Id("c1".into()))]
    fn test_atom_ref_roundtrip_edn_string(#[case] input: &str, #[case] expected: AtomRef) {
        let tree = read_string(input).unwrap();
        let parsed = AtomRef::from_edn(&tree).unwrap();
        assert_eq!(parsed, expected);
        let rendered = parsed.to_edn();
        let reparsed = AtomRef::from_edn(&rendered).unwrap();
        assert_eq!(reparsed, expected);
    }

    #[rstest]
    fn test_atom_ref_from_ast_uses_id_when_present(meta_with_atom_id: MoleculeMetadata) {
        let r = AtomRef::from_ast(AtomId(2), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Id("c1".into()));
    }

    #[rstest]
    fn test_atom_ref_from_ast_falls_back_to_index_without_id(meta_with_atom_id: MoleculeMetadata) {
        let r = AtomRef::from_ast(AtomId(4), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Index(4));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_id(meta_with_atom_id: MoleculeMetadata) {
        let id = AtomRef::Id("c1".into())
            .into_ast(5, &meta_with_atom_id)
            .unwrap();
        assert_eq!(id, AtomId(2));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_index(meta_with_atom_id: MoleculeMetadata) {
        let id = AtomRef::Index(3).into_ast(5, &meta_with_atom_id).unwrap();
        assert_eq!(id, AtomId(3));
    }

    #[rstest]
    fn test_atom_ref_into_ast_out_of_range_index(meta_with_atom_id: MoleculeMetadata) {
        let err = AtomRef::Index(9)
            .into_ast(5, &meta_with_atom_id)
            .unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "9".into(),
            }
        );
    }

    #[rstest]
    fn test_atom_ref_into_ast_unknown_id(meta_with_atom_id: MoleculeMetadata) {
        let err = AtomRef::Id("nope".into())
            .into_ast(5, &meta_with_atom_id)
            .unwrap_err();
        assert_eq!(
            err,
            ParseError::InvalidRef {
                kind: "atom",
                value: "nope".into(),
            }
        );
    }
}
