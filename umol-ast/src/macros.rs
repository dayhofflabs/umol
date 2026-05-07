//! Construction macros for the AST and DSL surface.
//!
//! All four are runtime parsers that panic on malformed input — appropriate
//! for tests, examples, and inline literal data. They wrap the corresponding
//! `FromStr` impl with `.unwrap()`.

/// Parse a molecule-EDN string into a `MoleculeAst`. Metadata (atom IDs,
/// aliases, etc.) in the input is dropped. Use [`dsl!`] to keep metadata.
///
/// ```ignore
/// let m = mol!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
/// ```
#[macro_export]
macro_rules! mol {
    ($s:expr $(,)?) => {{
        <$crate::ast::MoleculeAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a molecule-EDN string into a `MoleculeDsl`, preserving metadata
/// (atom IDs, aliases, etc.). Use [`mol!`] when you only need the AST.
///
/// ```ignore
/// let d = dsl!(r#"{:atom-aliases [:c "C"] :atoms [:c :c] :bonds [[0 1 "1"]]}"#);
/// ```
#[macro_export]
macro_rules! dsl {
    ($s:expr $(,)?) => {{
        <$crate::dsl::MoleculeDsl as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a compact atom-string into an `AtomAst`.
///
/// ```ignore
/// let a = atom!("C#h=#a+");
/// ```
#[macro_export]
macro_rules! atom {
    ($s:expr $(,)?) => {{
        <$crate::ast::AtomAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a compact bond-string into a `BondAst`.
///
/// ```ignore
/// let b = bond!("1#a");
/// ```
#[macro_export]
macro_rules! bond {
    ($s:expr $(,)?) => {{
        <$crate::ast::BondAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::element::Element;

    use crate::ast::{AtomIdx, BondIdx, ElementAst, MoleculeAst, ValueAst};

    #[rstest]
    fn test_mol_macro_parses_to_molecule_ast() {
        let m: MoleculeAst = mol!(r#"{:atoms ["C" "O"] :bonds [[0 1 "2"]]}"#);
        assert_eq!(m.atom_count(), 2);
        assert_eq!(m.bond_count(), 1);
        assert_eq!(m.atom(AtomIdx(0)).data.element, ElementAst::Lit(Element::C));
        assert_eq!(m.atom(AtomIdx(1)).data.element, ElementAst::Lit(Element::O));
        assert_eq!(m.bond(BondIdx(0)).data.order, ValueAst::Lit(2));
    }

    #[rstest]
    fn test_dsl_macro_preserves_metadata() {
        let d = dsl!(r#"{:atom-aliases [:c "C"] :atoms [:c :c] :bonds [[0 1 "1"]]}"#);
        assert_eq!(d.ast().atom_count(), 2);
        assert!(d.metadata().has_atom_alias("c"));
    }

    #[rstest]
    fn test_atom_macro_parses_to_atom_ast() {
        let a = atom!("C#c+");
        assert_eq!(a.element, ElementAst::Lit(Element::C));
        assert_eq!(a.charge, ValueAst::Lit(1));
    }

    #[rstest]
    fn test_bond_macro_parses_to_bond_ast() {
        let b = bond!("2");
        assert_eq!(b.order, ValueAst::Lit(2));
    }

    #[rstest]
    #[should_panic]
    fn test_mol_macro_panics_on_malformed() {
        let _: MoleculeAst = mol!("not valid edn");
    }

    #[rstest]
    #[should_panic]
    fn test_atom_macro_panics_on_malformed() {
        let _ = atom!("definitely not an atom");
    }
}
