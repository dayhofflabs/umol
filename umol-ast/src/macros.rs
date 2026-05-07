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

/// Parse a molecule-EDN string into a `MoleculeAst` with `MoleculeDefaults::zeroed()`
/// applied — every undetermined per-atom and per-bond field is filled with
/// the zero-policy default (charge 0, lone pairs 0, normal H, closed-shell
/// spin, etc.). Use this for ground-state fixture inputs; use [`mol!`] when
/// you want the input to pass through verbatim.
///
/// ```ignore
/// let methane = mol_zeroed!(r#"{:atoms ["C #h4"] :bonds []}"#);
/// ```
#[macro_export]
macro_rules! mol_zeroed {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::MoleculeDsl =
            <$crate::dsl::MoleculeDsl as ::core::str::FromStr>::from_str($s).unwrap();
        let (ast, _meta) = dsl.into_parts();
        <$crate::dsl::MoleculeDsl as $crate::ast::IntoAst<$crate::ast::MoleculeAst>>::into_ast(
            $crate::dsl::MoleculeDsl::from_parts(ast, $crate::dsl::Metadata::default()),
            &$crate::dsl::MoleculeDefaults::zeroed(),
        )
    }};
}

/// Parse a compact atom-string into an `AtomAst` with `AtomDefaults::zeroed()`
/// applied.
#[macro_export]
macro_rules! atom_zeroed {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::AtomDsl =
            <$crate::dsl::AtomDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::AtomDsl as $crate::ast::IntoAst<$crate::ast::AtomAst>>::into_ast(
            dsl,
            &$crate::dsl::AtomDefaults::zeroed(),
        )
    }};
}

/// Parse a compact bond-string into a `BondAst` with `BondDefaults::zeroed()`
/// applied.
#[macro_export]
macro_rules! bond_zeroed {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::BondDsl =
            <$crate::dsl::BondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::BondDsl as $crate::ast::IntoAst<$crate::ast::BondAst>>::into_ast(
            dsl,
            &$crate::dsl::BondDefaults::zeroed(),
        )
    }};
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::element::Element;

    use crate::ast::{
        AtomIdx, BondIdx, ElementAst, ImplicitHydrogensAst, IsotopeAst, MoleculeAst, ValueAst,
    };

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

    #[rstest]
    fn test_mol_zeroed_macro_fills_ground_state_defaults() {
        // Methane via mol_zeroed: only #h4 supplied; other fields filled
        // by zeroed defaults (isotope=Natural, charge=0, lone_pairs=0,
        // closed-shell spin).
        let m = mol_zeroed!(r#"{:atoms ["C #h4"] :bonds []}"#);
        let atom = m.atom(AtomIdx(0)).data;
        assert_eq!(atom.element, ElementAst::Lit(Element::C));
        assert_eq!(atom.isotope_mass, IsotopeAst::Natural);
        assert_eq!(atom.charge, ValueAst::Lit(0));
        assert_eq!(atom.implicit_hydrogens, ImplicitHydrogensAst::Lit(4));
        assert_eq!(atom.lone_pairs, ValueAst::Lit(0));
    }

    #[rstest]
    fn test_atom_zeroed_macro_fills_ground_state_defaults() {
        let a = atom_zeroed!("C #h4");
        assert_eq!(a.element, ElementAst::Lit(Element::C));
        assert_eq!(a.isotope_mass, IsotopeAst::Natural);
        assert_eq!(a.charge, ValueAst::Lit(0));
        assert_eq!(a.implicit_hydrogens, ImplicitHydrogensAst::Lit(4));
        assert_eq!(a.lone_pairs, ValueAst::Lit(0));
    }

    #[rstest]
    fn test_bond_zeroed_macro_fills_ground_state_defaults() {
        let b = bond_zeroed!("1");
        assert_eq!(b.order, ValueAst::Lit(1));
        assert_eq!(b.charge, ValueAst::Lit(0));
    }
}
