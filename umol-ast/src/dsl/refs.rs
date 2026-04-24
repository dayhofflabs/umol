//! Surface-level references to atoms, bonds, aromatic systems, and
//! multicenter bonds. Each ref is either a positional index (`Edn::Int`) or a
//! symbolic id (`Edn::Symbol`).

use umol_edn::{DeError, Edn, EdnSymbol, FromEdn, ToEdn};

macro_rules! define_ref_dsl {
    ($name:ident, $kind:literal) => {
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub enum $name {
            Index(usize),
            Id(String),
        }

        impl<'de> FromEdn<'de> for $name {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                match edn {
                    Edn::Int(n) => {
                        usize::try_from(*n)
                            .map(Self::Index)
                            .map_err(|_| DeError::OutOfRange {
                                value: n.to_string(),
                                target: "usize",
                                path: Vec::new(),
                            })
                    }
                    Edn::Symbol(s) => Ok(Self::Id(s.as_str().to_string())),
                    other => Err(DeError::TypeMismatch {
                        expected: concat!($kind, "-ref (int or symbol)"),
                        got: other.kind(),
                        path: Vec::new(),
                    }),
                }
            }
        }

        impl ToEdn for $name {
            fn to_edn(&self) -> Edn<'static> {
                match self {
                    Self::Index(n) => Edn::Int(*n as i64),
                    Self::Id(s) => Edn::Symbol(EdnSymbol::owned(s.clone())),
                }
            }
        }
    };
}

define_ref_dsl!(AtomRefDsl, "atom");
define_ref_dsl!(BondRefDsl, "bond");
define_ref_dsl!(AromaticSystemRefDsl, "aromatic-system");
define_ref_dsl!(MulticenterBondRefDsl, "multicenter-bond");
define_ref_dsl!(DativeBondRefDsl, "dative-bond");
define_ref_dsl!(NoncovalentBondRefDsl, "noncovalent-bond");

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::{read_string, EdnSymbol};

    use super::*;

    #[rstest]
    #[case::index("3", AtomRefDsl::Index(3))]
    #[case::id("foo", AtomRefDsl::Id("foo".to_string()))]
    #[case::id_compound("a1", AtomRefDsl::Id("a1".to_string()))]
    fn test_atom_ref_dsl_from_edn(#[case] input: &str, #[case] expected: AtomRefDsl) {
        let edn = read_string(input).unwrap();
        let parsed = AtomRefDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed, expected);
    }

    #[rstest]
    #[case::index(AtomRefDsl::Index(3), Edn::Int(3))]
    #[case::id(AtomRefDsl::Id("foo".to_string()), Edn::Symbol(EdnSymbol::new("foo")))]
    fn test_atom_ref_dsl_to_edn(#[case] input: AtomRefDsl, #[case] expected: Edn<'static>) {
        assert_eq!(input.to_edn(), expected);
    }

    #[rstest]
    #[case::index(AtomRefDsl::Index(0))]
    #[case::index_large(AtomRefDsl::Index(1_000_000))]
    #[case::id(AtomRefDsl::Id("a".to_string()))]
    #[case::id_complex(AtomRefDsl::Id("atom-1".to_string()))]
    fn test_atom_ref_dsl_roundtrip(#[case] input: AtomRefDsl) {
        let edn = input.to_edn();
        let parsed = AtomRefDsl::from_edn(&edn).unwrap();
        assert_eq!(input, parsed);
    }

    #[rstest]
    #[case::nil("nil")]
    #[case::string("\"foo\"")]
    #[case::keyword(":foo")]
    #[case::negative_int("-1")]
    fn test_atom_ref_dsl_from_edn_error(#[case] input: &str) {
        let edn = read_string(input).unwrap();
        assert!(AtomRefDsl::from_edn(&edn).is_err());
    }

    #[rstest]
    fn test_bond_ref_dsl_roundtrip() {
        let r = BondRefDsl::Id("b1".to_string());
        let edn = r.to_edn();
        let parsed = BondRefDsl::from_edn(&edn).unwrap();
        assert_eq!(r, parsed);
    }

    #[rstest]
    fn test_aromatic_system_ref_dsl_roundtrip() {
        let r = AromaticSystemRefDsl::Index(2);
        let edn = r.to_edn();
        let parsed = AromaticSystemRefDsl::from_edn(&edn).unwrap();
        assert_eq!(r, parsed);
    }

    #[rstest]
    fn test_multicenter_bond_ref_dsl_roundtrip() {
        let r = MulticenterBondRefDsl::Id("mc".to_string());
        let edn = r.to_edn();
        let parsed = MulticenterBondRefDsl::from_edn(&edn).unwrap();
        assert_eq!(r, parsed);
    }

    #[rstest]
    fn test_dative_bond_ref_dsl_roundtrip() {
        let r = DativeBondRefDsl::Index(5);
        let edn = r.to_edn();
        let parsed = DativeBondRefDsl::from_edn(&edn).unwrap();
        assert_eq!(r, parsed);
    }

    #[rstest]
    fn test_noncovalent_bond_ref_dsl_roundtrip() {
        let r = NoncovalentBondRefDsl::Id("nc".to_string());
        let edn = r.to_edn();
        let parsed = NoncovalentBondRefDsl::from_edn(&edn).unwrap();
        assert_eq!(r, parsed);
    }
}
