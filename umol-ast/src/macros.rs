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
/// let a = atom!("C#h*#a+");
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

/// Parse molecule-EDN string into a `MoleculeAst` with `MoleculeDefaults::ground()` applied.
/// Mirrors `AtomAst::into_ground()` at the molecule scope.
#[macro_export]
macro_rules! mol_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::MoleculeDsl =
            <$crate::dsl::MoleculeDsl as ::core::str::FromStr>::from_str($s).unwrap();
        let (ast, _meta) = dsl.into_parts();
        <$crate::dsl::MoleculeDsl as $crate::ast::IntoAst<$crate::ast::MoleculeAst>>::into_ast(
            $crate::dsl::MoleculeDsl::from_parts(ast, $crate::dsl::Metadata::default()),
            &$crate::dsl::MoleculeDefaults::ground(),
        )
    }};
}

/// Parse atom DSL into an `AtomAst` with `AtomDefaults::ground()` applied.
#[macro_export]
macro_rules! atom_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::AtomDsl =
            <$crate::dsl::AtomDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::AtomDsl as $crate::ast::IntoAst<$crate::ast::AtomAst>>::into_ast(
            dsl,
            &$crate::dsl::AtomDefaults::ground(),
        )
    }};
}

/// Parse a bond DSL into a `BondAst` with `BondDefaults::ground()` applied.
#[macro_export]
macro_rules! bond_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::BondDsl =
            <$crate::dsl::BondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::BondDsl as $crate::ast::IntoAst<$crate::ast::BondAst>>::into_ast(
            dsl,
            &$crate::dsl::BondDefaults::ground(),
        )
    }};
}

/// Parse a compact dative-bond-string into a `DativeBondAst`.
#[macro_export]
macro_rules! dative {
    ($s:expr $(,)?) => {{
        <$crate::ast::DativeBondAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a dative-bond DSL string into a `DativeBondAst` with `DativeBondDefaults::ground()` applied.
#[macro_export]
macro_rules! dative_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::DativeBondDsl =
            <$crate::dsl::DativeBondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::DativeBondDsl as $crate::ast::IntoAst<$crate::ast::DativeBondAst>>::into_ast(
            dsl,
            &$crate::dsl::DativeBondDefaults::ground(),
        )
    }};
}

/// Parse a compact aromatic-system-string into an `AromaticSystemAst`.
#[macro_export]
macro_rules! aromatic {
    ($s:expr $(,)?) => {{
        <$crate::ast::AromaticSystemAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse an aromatic-system DSL string into an `AromaticSystemAst` with
/// `AromaticSystemDefaults::ground()` applied.
#[macro_export]
macro_rules! aromatic_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::AromaticSystemDsl =
            <$crate::dsl::AromaticSystemDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::AromaticSystemDsl as $crate::ast::IntoAst<$crate::ast::AromaticSystemAst>>::into_ast(
            dsl,
            &$crate::dsl::AromaticSystemDefaults::ground(),
        )
    }};
}

/// Parse a compact multicenter-bond-string into a `MulticenterBondAst`.
#[macro_export]
macro_rules! multicenter {
    ($s:expr $(,)?) => {{
        <$crate::ast::MulticenterBondAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a multicenter-bond DSL string into a `MulticenterBondAst` with
/// `MulticenterBondDefaults::ground()` applied.
#[macro_export]
macro_rules! multicenter_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::MulticenterBondDsl =
            <$crate::dsl::MulticenterBondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::MulticenterBondDsl as $crate::ast::IntoAst<$crate::ast::MulticenterBondAst>>::into_ast(
            dsl,
            &$crate::dsl::MulticenterBondDefaults::ground(),
        )
    }};
}

/// Parse a compact noncovalent-bond-string into a `NoncovalentBondAst`.
#[macro_export]
macro_rules! noncovalent {
    ($s:expr $(,)?) => {{
        <$crate::ast::NoncovalentBondAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a noncovalent-bond DSL string into a `NoncovalentBondAst` with
/// `NoncovalentBondDefaults::ground()` applied.
#[macro_export]
macro_rules! noncovalent_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::NoncovalentBondDsl =
            <$crate::dsl::NoncovalentBondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::NoncovalentBondDsl as $crate::ast::IntoAst<$crate::ast::NoncovalentBondAst>>::into_ast(
            dsl,
            &$crate::dsl::NoncovalentBondDefaults::ground(),
        )
    }};
}

/// Parse a compact stereo-atom-string into a `StereoAtomAst`.
#[macro_export]
macro_rules! stereo_atom {
    ($s:expr $(,)?) => {{
        <$crate::ast::StereoAtomAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a stereo-atom DSL string into a `StereoAtomAst` with `StereoAtomDefaults::ground()` applied.
#[macro_export]
macro_rules! stereo_atom_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::StereoAtomDsl =
            <$crate::dsl::StereoAtomDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::StereoAtomDsl as $crate::ast::IntoAst<$crate::ast::StereoAtomAst>>::into_ast(
            dsl,
            &$crate::dsl::StereoAtomDefaults::ground(),
        )
    }};
}

/// Parse a compact stereo-bond-string into a `StereoBondAst`.
#[macro_export]
macro_rules! stereo_bond {
    ($s:expr $(,)?) => {{
        <$crate::ast::StereoBondAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a stereo-bond DSL string into a `StereoBondAst` with `StereoBondDefaults::ground()` applied.
#[macro_export]
macro_rules! stereo_bond_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::StereoBondDsl =
            <$crate::dsl::StereoBondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::StereoBondDsl as $crate::ast::IntoAst<$crate::ast::StereoBondAst>>::into_ast(
            dsl,
            &$crate::dsl::StereoBondDefaults::ground(),
        )
    }};
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;

    use crate::ast::constraint::RingScope;
    use crate::ast::{
        AromaticSystemAst, AromaticSystemConstraint, AtomAst, AtomConstraint, AtomId, BondAst,
        BondConstraint, Constraints, DativeBondAst, DativeBondConstraint, MoleculeAst,
        MulticenterBondAst, NoncovalentBondAst, NoncovalentBondKind, StereoAtomAst, StereoBondAst,
        StereoCosetAst, StereoKind,
    };
    use crate::dsl::molecule::Metadata;
    use crate::dsl::{AtomDsl, MoleculeDsl};

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("{}", MoleculeAst::default())]
    #[case::carbon_oxygen(r#"{:atoms ["C #h2" "O"] :bonds [[0 1 "2"]]}"#,
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C).with_implicit_hydrogens(2_i64), AtomAst::from_element(Element::O)],
        vec![(AtomId(0), AtomId(1), BondAst::from_order(2))]))]
    #[case::aromatic_system(r##"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :electrons [1 1 1] :type "#e3"}]}"##,
        MoleculeAst::from_parts(vec![AtomAst::from_element(Element::C); 3],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1)), (AtomId(1), AtomId(2), BondAst::from_order(1)), (AtomId(2), AtomId(0), BondAst::from_order(1))],
            vec![], vec![(vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemAst::from_counts(vec![1; 3]).with_constraint(AromaticSystemConstraint::electron_count(3)))],
            vec![], vec![],
            Vec::new(), Vec::new(), Constraints::default()))]
    fn test_mol_macro(#[case] input: &str, #[case] expected: MoleculeAst) {
        assert_eq!(mol!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_mol_macro_error() {
        let _: MoleculeAst = mol!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("{}", MoleculeDsl::default())]
    #[case::with_alias(
        r#"{:atom-aliases [:c "C"] :atoms [:c :c] :bonds [[0 1 "1"]]}"#,
        MoleculeDsl::from_parts(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C); 2],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Metadata::new().with_atom_alias("c", "C".parse::<AtomDsl>().unwrap()),
        ),
    )]
    #[case::with_atom_ids(
        r#"{:atoms [[:a "C"] [:b "C"]] :bonds []}"#,
        MoleculeDsl::from_parts(
            MoleculeAst::from_atoms_and_bonds(
                vec![AtomAst::from_element(Element::C); 2],
                vec![],
            ),
            Metadata::new().with_atom_id(AtomId(0), "a").with_atom_id(AtomId(1), "b"),
        ),
    )]
    fn test_dsl_macro(#[case] input: &str, #[case] expected: MoleculeDsl) {
        assert_eq!(dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_dsl_macro_error() {
        let _ = dsl!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon_charge("C#c+", AtomAst::from_element(Element::C).with_charge(1_i64))]
    fn test_atom_macro(#[case] input: &str, #[case] expected: AtomAst) {
        assert_eq!(atom!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_atom_macro_error() {
        let _ = atom!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::double("2", BondAst::from_order(2))]
    #[case::aromatic("1#a", BondAst::from_order(1).with_constraint(BondConstraint::Aromatic))]
    fn test_bond_macro(#[case] input: &str, #[case] expected: BondAst) {
        assert_eq!(bond!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_bond_macro_error() {
        let _ = bond!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::methane(r#"{:atoms ["C #h4"] :bonds []}"#,
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C).with_implicit_hydrogens(4_i64).into_ground()], vec![]))]
    #[case::carbon_charged(r#"{:atoms ["C #c+"] :bonds []}"#,
        MoleculeAst::from_atoms_and_bonds(vec![AtomAst::from_element(Element::C).with_charge(1_i64).into_ground()], vec![]))]
    fn test_mol_ground_macro(#[case] input: &str, #[case] expected: MoleculeAst) {
        assert_eq!(mol_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon_h4("C #h4", AtomAst::from_element(Element::C).with_implicit_hydrogens(4_i64).into_ground())]
    #[case::carbon("C", AtomAst::from_element(Element::C).into_ground())]
    #[case::carbon_v4("C #v4", AtomAst::from_element(Element::C).with_constraint(AtomConstraint::valence(4_i64)).into_ground())]
    fn test_atom_ground_macro(#[case] input: &str, #[case] expected: AtomAst) {
        assert_eq!(atom_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::double("2", BondAst::from_order(2).into_ground())]
    fn test_bond_ground_macro(#[case] input: &str, #[case] expected: BondAst) {
        assert_eq!(bond_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", DativeBondAst::from_order(1))]
    #[case::with_ring_size( "2 #R(6)", DativeBondAst::from_order(2).with_constraint(DativeBondConstraint::ring_membership(RingScope::Size(6), 1)),)]
    fn test_dative_macro(#[case] input: &str, #[case] expected: DativeBondAst) {
        assert_eq!(dative!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_dative_macro_error() {
        let _ = dative!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::double("2", DativeBondAst::from_order(2).into_ground())]
    fn test_dative_ground_macro(#[case] input: &str, #[case] expected: DativeBondAst) {
        assert_eq!(dative_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge("#c+", AromaticSystemAst::default().with_charge(1_i64))]
    #[case::electrons("#e6", AromaticSystemAst::default().with_constraint(AromaticSystemConstraint::electron_count(6)))]
    fn test_aromatic_macro(#[case] input: &str, #[case] expected: AromaticSystemAst) {
        assert_eq!(aromatic!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_aromatic_macro_error() {
        let _ = aromatic!("not_a_predicate");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons("#e6", AromaticSystemAst::default().with_constraint(AromaticSystemConstraint::electron_count(6)).into_ground())]
    fn test_aromatic_ground_macro(#[case] input: &str, #[case] expected: AromaticSystemAst) {
        assert_eq!(aromatic_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge("#c-", MulticenterBondAst::default().with_charge(-1_i64))]
    fn test_multicenter_macro(#[case] input: &str, #[case] expected: MulticenterBondAst) {
        assert_eq!(multicenter!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_multicenter_macro_error() {
        let _ = multicenter!("not_a_predicate");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charged("#c-", MulticenterBondAst::default().with_charge(-1_i64).into_ground())]
    fn test_multicenter_ground_macro(#[case] input: &str, #[case] expected: MulticenterBondAst) {
        assert_eq!(multicenter_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hbond("Hbd", NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    #[case::ionic("Ion", NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic))]
    fn test_noncovalent_macro(#[case] input: &str, #[case] expected: NoncovalentBondAst) {
        assert_eq!(noncovalent!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_noncovalent_macro_error() {
        let _ = noncovalent!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hbond("Hbd", NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).into_ground())]
    fn test_noncovalent_ground_macro(#[case] input: &str, #[case] expected: NoncovalentBondAst) {
        assert_eq!(noncovalent_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ccw("Th0", StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)))]
    #[case::undetermined("Th*", StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined))]
    #[case::square_planar("Sp2", StereoAtomAst::new(StereoKind::SquarePlanar, StereoCosetAst::Lit(2)))]
    fn test_stereo_atom_macro(#[case] input: &str, #[case] expected: StereoAtomAst) {
        assert_eq!(stereo_atom!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_stereo_atom_macro_error() {
        let _ = stereo_atom!("Th!");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ccw("Th0", StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)))]
    fn test_stereo_atom_ground_macro(#[case] input: &str, #[case] expected: StereoAtomAst) {
        assert_eq!(stereo_atom_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::z("Ct1", StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)))]
    #[case::undetermined("Ct*", StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Undetermined))]
    fn test_stereo_bond_macro(#[case] input: &str, #[case] expected: StereoBondAst) {
        assert_eq!(stereo_bond!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_stereo_bond_macro_error() {
        let _ = stereo_bond!("Ct!");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::z("Ct1", StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)))]
    fn test_stereo_bond_ground_macro(#[case] input: &str, #[case] expected: StereoBondAst) {
        assert_eq!(stereo_bond_ground!(input), expected);
    }
}
