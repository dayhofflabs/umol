//! Construction macros for the AST and DSL objects.

/// Parse a molecule-EDN string into a `MoleculeAst`. MoleculeMetadata (atom IDs,
/// aliases, etc.) in the input is dropped. Use `MoleculeDsl::from_str` to keep metadata.
///
/// ```ignore
/// let m = mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#);
/// ```
#[macro_export]
macro_rules! mol_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::MoleculeAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse molecule-EDN string into a `MoleculeAst` with `MoleculeDefaults::ground()` applied.
/// Mirrors `AtomForm::into_ground()` at the molecule scope.
#[macro_export]
macro_rules! mol_dsl_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::MoleculeDsl =
            <$crate::dsl::MoleculeDsl as ::core::str::FromStr>::from_str($s).unwrap();
        let (ast, _meta) = dsl.into_parts();
        <$crate::dsl::MoleculeDsl as $crate::ir::IntoIr<$crate::ir::MoleculeAst>>::into_ir(
            $crate::dsl::MoleculeDsl::new(ast, $crate::dsl::MoleculeMetadata::default())
                .expect("empty metadata is coherent"),
            &$crate::dsl::MoleculeDefaults::ground(),
        )
    }};
}

/// Parse a compact atom-string into an `AtomForm`.
///
/// ```ignore
/// let a = atom_dsl!("C#h*#a+");
/// ```
#[macro_export]
macro_rules! atom_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::AtomForm as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse atom DSL into an `AtomForm` with `AtomDefaults::ground()` applied.
#[macro_export]
macro_rules! atom_dsl_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::AtomDsl =
            <$crate::dsl::AtomDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::AtomDsl as $crate::ir::IntoIr<$crate::ir::AtomForm>>::into_ir(
            dsl,
            &$crate::dsl::AtomDefaults::ground(),
        )
    }};
}

/// Parse a compact atom-update string into `AtomUpdate`.
#[macro_export]
macro_rules! atom_update_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::AtomUpdate as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a compact bond-string into a `BondForm`.
///
/// ```ignore
/// let b = bond_dsl!("1#a");
/// ```
#[macro_export]
macro_rules! bond_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::BondForm as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a bond DSL into a `BondForm` with `BondDefaults::ground()` applied.
#[macro_export]
macro_rules! bond_dsl_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::BondDsl =
            <$crate::dsl::BondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::BondDsl as $crate::ir::IntoIr<$crate::ir::BondForm>>::into_ir(
            dsl,
            &$crate::dsl::BondDefaults::ground(),
        )
    }};
}

/// Parse a compact bond-update string into `BondUpdate`.
#[macro_export]
macro_rules! bond_update_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::BondUpdate as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a compact dative-bond-string into a `DativeBondForm`.
#[macro_export]
macro_rules! dative_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::DativeBondForm as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a dative-bond DSL string into a `DativeBondForm` with `DativeBondDefaults::ground()` applied.
#[macro_export]
macro_rules! dative_dsl_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::DativeBondDsl =
            <$crate::dsl::DativeBondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::DativeBondDsl as $crate::ir::IntoIr<$crate::ir::DativeBondForm>>::into_ir(
            dsl,
            &$crate::dsl::DativeBondDefaults::ground(),
        )
    }};
}

/// Parse a compact aromatic-system-string into an `AromaticSystemForm`.
#[macro_export]
macro_rules! aromatic_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::AromaticSystemForm as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse an aromatic-system DSL string into an `AromaticSystemForm` with
/// `AromaticSystemDefaults::ground()` applied.
#[macro_export]
macro_rules! aromatic_dsl_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::AromaticSystemDsl =
            <$crate::dsl::AromaticSystemDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::AromaticSystemDsl as $crate::ir::IntoIr<$crate::ir::AromaticSystemForm>>::into_ir(
            dsl,
            &$crate::dsl::AromaticSystemDefaults::ground(),
        )
    }};
}

/// Parse a compact multicenter-bond-string into a `MulticenterBondForm`.
#[macro_export]
macro_rules! multicenter_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::MulticenterBondForm as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a multicenter-bond DSL string into a `MulticenterBondForm` with
/// `MulticenterBondDefaults::ground()` applied.
#[macro_export]
macro_rules! multicenter_dsl_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::MulticenterBondDsl =
            <$crate::dsl::MulticenterBondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::MulticenterBondDsl as $crate::ir::IntoIr<$crate::ir::MulticenterBondForm>>::into_ir(
            dsl,
            &$crate::dsl::MulticenterBondDefaults::ground(),
        )
    }};
}

/// Parse a compact noncovalent-bond-string into a `NoncovalentBondForm`.
#[macro_export]
macro_rules! noncovalent_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::NoncovalentBondForm as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a noncovalent-bond DSL string into a `NoncovalentBondForm` with
/// `NoncovalentBondDefaults::ground()` applied.
#[macro_export]
macro_rules! noncovalent_dsl_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::NoncovalentBondDsl =
            <$crate::dsl::NoncovalentBondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::NoncovalentBondDsl as $crate::ir::IntoIr<$crate::ir::NoncovalentBondForm>>::into_ir(
            dsl,
            &$crate::dsl::NoncovalentBondDefaults::ground(),
        )
    }};
}

/// Parse a compact stereo-atom-string into a `StereoAtomAst`.
#[macro_export]
macro_rules! stereo_atom_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::StereoAtomAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a stereo-atom DSL string into a `StereoAtomAst` with `StereoAtomDefaults::ground()` applied.
#[macro_export]
macro_rules! stereo_atom_dsl_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::StereoAtomDsl =
            <$crate::dsl::StereoAtomDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::StereoAtomDsl as $crate::ir::IntoIr<$crate::ir::StereoAtomAst>>::into_ir(
            dsl,
            &$crate::dsl::StereoAtomDefaults::ground(),
        )
    }};
}

/// Parse a compact stereo-bond-string into a `StereoBondAst`.
#[macro_export]
macro_rules! stereo_bond_dsl {
    ($s:expr $(,)?) => {{
        <$crate::ir::StereoBondAst as ::core::str::FromStr>::from_str($s).unwrap()
    }};
}

/// Parse a stereo-bond DSL string into a `StereoBondAst` with `StereoBondDefaults::ground()` applied.
#[macro_export]
macro_rules! stereo_bond_dsl_ground {
    ($s:expr $(,)?) => {{
        let dsl: $crate::dsl::StereoBondDsl =
            <$crate::dsl::StereoBondDsl as ::core::str::FromStr>::from_str($s).unwrap();
        <$crate::dsl::StereoBondDsl as $crate::ir::IntoIr<$crate::ir::StereoBondAst>>::into_ir(
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

    use crate::dsl::{AtomDsl, MoleculeDsl, MoleculeMetadata};
    use crate::ir::constraint::RingScope;
    use crate::ir::{
        AromaticSystemConstraintAst, AromaticSystemForm, AtomConstraintAst, AtomConstraintsAst,
        AtomForm, AtomId, AtomUpdate, BondConstraintAst, BondConstraintsAst, BondForm, BondUpdate,
        BooleanForm, DativeBondConstraintAst, DativeBondForm, ElementForm, Entity, MoleculeAst,
        MoleculeEntries, MulticenterBondForm, NoncovalentBondForm, NoncovalentBondKind, NumForm,
        StereoAtomAst, StereoBondAst, StereoCoset, StereoKind,
    };

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("{}", MoleculeDsl::default())]
    #[case::with_alias(
        r#"{:atom-aliases [:c "C"] :atoms [:c :c] :bonds [[0 1 "1"]]}"#,
        MoleculeDsl::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 2],
                bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
                ..Default::default()
            }),
            {
                let mut metadata = MoleculeMetadata::new();
                metadata
                    .add_atom_alias("c", "C".parse::<AtomDsl>().unwrap())
                    .unwrap();
                metadata
            },
        ).unwrap(),
    )]
    #[case::with_atom_ids(
        r#"{:atoms [[:a "C"] [:b "C"]] :bonds []}"#,
        MoleculeDsl::new(
            MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); 2],
                bonds: vec![],
                ..Default::default()
            }),
            {
                let mut metadata = MoleculeMetadata::new();
                metadata.set_keyword(Entity::Atom(AtomId(0)), "a").unwrap();
                metadata.set_keyword(Entity::Atom(AtomId(1)), "b").unwrap();
                metadata
            },
        ).unwrap(),
    )]
    fn test_dsl_macro(#[case] input: &str, #[case] expected: MoleculeDsl) {
        assert_eq!(input.parse::<MoleculeDsl>().unwrap(), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_dsl_macro_error() {
        let _ = "invalid".parse::<MoleculeDsl>().unwrap();
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("{}", MoleculeAst::default())]
    #[case::carbon_oxygen(r#"{:atoms ["C #h2" "O"] :bonds [[0 1 "2"]]}"#,
        MoleculeAst::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C).with_implicit_hydrogens(2_i64), AtomForm::from_element(Element::O)],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(2))], ..Default::default() }))]
    #[case::aromatic_system(r##"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]] :aromatic-systems [{:atoms [0 1 2] :type "[1,1,1]#e3"}]}"##,
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 3],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1)), (AtomId(1), AtomId(2), BondForm::from_order(1)), (AtomId(2), AtomId(0), BondForm::from_order(1))],
            aromatic: vec![(vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![1; 3]).with_constraint(AromaticSystemConstraintAst::electron_count(3)))],
            ..Default::default()
        }))]
    fn test_mol_macro(#[case] input: &str, #[case] expected: MoleculeAst) {
        assert_eq!(mol_dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_mol_macro_error() {
        let _: MoleculeAst = mol_dsl!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::methane(r#"{:atoms ["C #h4"] :bonds []}"#,
        MoleculeAst::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C).with_implicit_hydrogens(4_i64).into_ground()], bonds: vec![], ..Default::default() }))]
    #[case::carbon_charged(r#"{:atoms ["C #c+"] :bonds []}"#,
        MoleculeAst::from_entries(MoleculeEntries { atoms: vec![AtomForm::from_element(Element::C).with_charge(1_i64).into_ground()], bonds: vec![], ..Default::default() }))]
    fn test_mol_ground_macro(#[case] input: &str, #[case] expected: MoleculeAst) {
        assert_eq!(mol_dsl_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon_charge("C#c+", AtomForm::from_element(Element::C).with_charge(1_i64))]
    fn test_atom_macro(#[case] input: &str, #[case] expected: AtomForm) {
        assert_eq!(atom_dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_atom_macro_error() {
        let _ = atom_dsl!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon_h4("C #h4", AtomForm::from_element(Element::C).with_implicit_hydrogens(4_i64).into_ground())]
    #[case::carbon("C", AtomForm::from_element(Element::C).into_ground())]
    #[case::carbon_v4("C #v4", AtomForm::from_element(Element::C).with_constraint(AtomConstraintAst::valence(4_i64)).into_ground())]
    fn test_atom_ground_macro(#[case] input: &str, #[case] expected: AtomForm) {
        assert_eq!(atom_dsl_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", AtomUpdate::default())]
    #[case::element("C#h4", AtomUpdate { element: Some(ElementForm::Lit(Element::C)), implicit_hydrogens: Some(NumForm::Lit(4)), ..Default::default() })]
    #[case::field_only("#h4", AtomUpdate { implicit_hydrogens: Some(NumForm::Lit(4)), ..Default::default() })]
    #[case::constraint_only("#v4", AtomUpdate { constraints: AtomConstraintsAst::from(AtomConstraintAst::valence(4_i64)), ..Default::default() })]
    fn test_atom_update_macro(#[case] input: &str, #[case] expected: AtomUpdate) {
        assert_eq!(atom_update_dsl!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::double("2", BondForm::from_order(2))]
    #[case::aromatic("1#a", BondForm::from_order(1).with_constraint(BondConstraintAst::Aromatic(BooleanForm::Lit(true))))]
    fn test_bond_macro(#[case] input: &str, #[case] expected: BondForm) {
        assert_eq!(bond_dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_bond_macro_error() {
        let _ = bond_dsl!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::double("2", BondForm::from_order(2).into_ground())]
    fn test_bond_ground_macro(#[case] input: &str, #[case] expected: BondForm) {
        assert_eq!(bond_dsl_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", BondUpdate::default())]
    #[case::order("1", BondUpdate { order: Some(NumForm::Lit(1)), ..Default::default() })]
    #[case::field_only("#c+", BondUpdate { charge: Some(NumForm::Lit(1)), ..Default::default() })]
    #[case::constraint_only("#a", BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanForm::Lit(true))), ..Default::default() })]
    fn test_bond_update_macro(#[case] input: &str, #[case] expected: BondUpdate) {
        assert_eq!(bond_update_dsl!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", DativeBondForm::from_order(1))]
    #[case::with_ring_size( "2 #R(6)", DativeBondForm::from_order(2).with_constraint(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)),)]
    fn test_dative_macro(#[case] input: &str, #[case] expected: DativeBondForm) {
        assert_eq!(dative_dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_dative_macro_error() {
        let _ = dative_dsl!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::double("2", DativeBondForm::from_order(2).into_ground())]
    fn test_dative_ground_macro(#[case] input: &str, #[case] expected: DativeBondForm) {
        assert_eq!(dative_dsl_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge("[1,1,1]#c+", AromaticSystemForm::from_electrons(vec![1, 1, 1]).with_charge(1_i64))]
    #[case::electron_count("*#e6", AromaticSystemForm::default().with_constraint(AromaticSystemConstraintAst::electron_count(6)))]
    fn test_aromatic_macro(#[case] input: &str, #[case] expected: AromaticSystemForm) {
        assert_eq!(aromatic_dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_aromatic_macro_error() {
        let _ = aromatic_dsl!("not_a_predicate");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count("[1,1,1,1,1,1]#e6", AromaticSystemForm::from_electrons(vec![1; 6]).with_constraint(AromaticSystemConstraintAst::electron_count(6)).into_ground())]
    fn test_aromatic_ground_macro(#[case] input: &str, #[case] expected: AromaticSystemForm) {
        assert_eq!(aromatic_dsl_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge("[1,1,1,1,1]#c-", MulticenterBondForm::from_electrons(vec![1; 5]).with_charge(-1_i64))]
    fn test_multicenter_macro(#[case] input: &str, #[case] expected: MulticenterBondForm) {
        assert_eq!(multicenter_dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_multicenter_macro_error() {
        let _ = multicenter_dsl!("not_a_predicate");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charged("[1,1,1,1,1]#c-", MulticenterBondForm::from_electrons(vec![1; 5]).with_charge(-1_i64).into_ground())]
    fn test_multicenter_ground_macro(#[case] input: &str, #[case] expected: MulticenterBondForm) {
        assert_eq!(multicenter_dsl_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hbond("Hbd", NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond))]
    #[case::ionic("Ion", NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic))]
    fn test_noncovalent_macro(#[case] input: &str, #[case] expected: NoncovalentBondForm) {
        assert_eq!(noncovalent_dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_noncovalent_macro_error() {
        let _ = noncovalent_dsl!("invalid");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hbond("Hbd", NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond).into_ground())]
    fn test_noncovalent_ground_macro(#[case] input: &str, #[case] expected: NoncovalentBondForm) {
        assert_eq!(noncovalent_dsl_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ccw("Th0", StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)))]
    #[case::undetermined("Th*", StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Undetermined))]
    #[case::square_planar("Sp2", StereoAtomAst::new(StereoKind::SquarePlanar, StereoCoset::Lit(2)))]
    fn test_stereo_atom_macro(#[case] input: &str, #[case] expected: StereoAtomAst) {
        assert_eq!(stereo_atom_dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_stereo_atom_macro_error() {
        let _ = stereo_atom_dsl!("Th!");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ccw("Th0", StereoAtomAst::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)))]
    fn test_stereo_atom_ground_macro(#[case] input: &str, #[case] expected: StereoAtomAst) {
        assert_eq!(stereo_atom_dsl_ground!(input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::z("Ct1", StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)))]
    #[case::undetermined("Ct*", StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Undetermined))]
    fn test_stereo_bond_macro(#[case] input: &str, #[case] expected: StereoBondAst) {
        assert_eq!(stereo_bond_dsl!(input), expected);
    }

    #[rstest]
    #[should_panic]
    fn test_stereo_bond_macro_error() {
        let _ = stereo_bond_dsl!("Ct!");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::z("Ct1", StereoBondAst::new(StereoKind::CisTrans, StereoCoset::Lit(1)))]
    fn test_stereo_bond_ground_macro(#[case] input: &str, #[case] expected: StereoBondAst) {
        assert_eq!(stereo_bond_dsl_ground!(input), expected);
    }
}
