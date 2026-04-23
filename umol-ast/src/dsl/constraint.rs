//! Tree-shaped constraint DSLs.
//!
//! Boundary types between the AST `Constraint` tree and its EDN form. Refs in
//! the tree carry either an integer index or a symbolic id; resolution to /
//! from the `AtomIdx` / `BondIdx` / ... on the AST is a separate fallible
//! step that consults the surrounding `Metadata`.

use umol_edn::{DeError, Edn, FromEdn, ToEdn};

use super::error::ParseError;
use super::molecule::Metadata;
use crate::ast::idx::{
    AromaticSystemIdx, AtomIdx, BondIdx, DativeBondIdx, MulticenterBondIdx, NoncovalentBondIdx,
};

macro_rules! define_ref {
    ($name:ident, $idx:ident, $field:ident, $kind:literal) => {
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub enum $name {
            Index(usize),
            Id(String),
        }

        impl $name {
            /// Build a ref from an AST index, preferring an id from `metadata`
            /// if one is recorded for this index.
            pub fn from_ast(idx: $idx, metadata: &Metadata) -> Self {
                if let Some(name) = metadata.$field.get(&idx) {
                    Self::Id(name.clone())
                } else {
                    Self::Index(idx.index())
                }
            }

            /// Resolve this ref to an AST index against `metadata`. Fails on
            /// unknown id or out-of-range numeric index.
            pub fn into_ast(self, count: usize, metadata: &Metadata) -> Result<$idx, ParseError> {
                match self {
                    Self::Index(i) => {
                        if i < count {
                            Ok($idx::from(i))
                        } else {
                            Err(ParseError::InvalidRef {
                                kind: $kind,
                                value: i.to_string(),
                            })
                        }
                    }
                    Self::Id(name) => metadata
                        .$field
                        .iter()
                        .find(|(_, n)| n.as_str() == name)
                        .map(|(idx, _)| *idx)
                        .ok_or(ParseError::InvalidRef {
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
                    Self::Id(name) => {
                        Edn::Keyword(umol_edn::EdnKeyword::owned(name.clone()))
                    }
                }
            }
        }
    };
}

define_ref!(AtomRef, AtomIdx, atom_ids, "atom");
define_ref!(BondRef, BondIdx, bond_ids, "bond");
define_ref!(DativeBondRef, DativeBondIdx, dative_bond_ids, "dative-bond");
define_ref!(
    AromaticSystemRef,
    AromaticSystemIdx,
    aromatic_system_ids,
    "aromatic-system"
);
define_ref!(
    MulticenterBondRef,
    MulticenterBondIdx,
    multicenter_bond_ids,
    "multicenter-bond"
);
define_ref!(
    NoncovalentBondRef,
    NoncovalentBondIdx,
    noncovalent_bond_ids,
    "noncovalent-bond"
);

#[cfg(test)]
mod tests {
    use bimap::BiMap;
    use indexmap::IndexMap;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::{read_string, EdnKeyword};

    use super::*;

    #[fixture]
    fn meta_with_atom_id() -> Metadata {
        let mut atom_ids = IndexMap::new();
        atom_ids.insert(AtomIdx(2), "c1".to_string());
        Metadata {
            atom_ids,
            atom_aliases: BiMap::new(),
            bond_ids: IndexMap::new(),
            dative_bond_ids: IndexMap::new(),
            aromatic_system_ids: IndexMap::new(),
            multicenter_bond_ids: IndexMap::new(),
            noncovalent_bond_ids: IndexMap::new(),
        }
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
        assert!(matches!(err, DeError::TypeMismatch { expected: "atom ref (int or keyword)", .. }));
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
    fn test_atom_ref_from_ast_uses_id_when_present(meta_with_atom_id: Metadata) {
        let r = AtomRef::from_ast(AtomIdx(2), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Id("c1".into()));
    }

    #[rstest]
    fn test_atom_ref_from_ast_falls_back_to_index_without_id(meta_with_atom_id: Metadata) {
        let r = AtomRef::from_ast(AtomIdx(4), &meta_with_atom_id);
        assert_eq!(r, AtomRef::Index(4));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_id(meta_with_atom_id: Metadata) {
        let idx = AtomRef::Id("c1".into())
            .into_ast(5, &meta_with_atom_id)
            .unwrap();
        assert_eq!(idx, AtomIdx(2));
    }

    #[rstest]
    fn test_atom_ref_into_ast_resolves_index(meta_with_atom_id: Metadata) {
        let idx = AtomRef::Index(3)
            .into_ast(5, &meta_with_atom_id)
            .unwrap();
        assert_eq!(idx, AtomIdx(3));
    }

    #[rstest]
    fn test_atom_ref_into_ast_out_of_range_index(meta_with_atom_id: Metadata) {
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
    fn test_atom_ref_into_ast_unknown_id(meta_with_atom_id: Metadata) {
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
