//! Stereo DSL serialization properties.
//!
//! Lossless generated forms retain exact roundtrips. Raw normalizable entities use
//! normalized_eq, while standalone kindless updates omit constraints explicitly.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]
    #[test]
    fn test_stereo_atom_dsl_display_from_str_roundtrip(
        stereo in stereo_atom_form_strategy(),
    ) {
        let dsl = StereoAtomDsl(stereo);
        let rendered = dsl.to_string();
        let parsed: StereoAtomDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_stereo_atom_dsl_roundtrip_raw(form in raw_stereo_atom_form_strategy()) {
        if form.clone().normalize().is_ok() {
            let dsl = StereoAtomDsl(form.clone());
            let text: StereoAtomDsl = dsl.to_string().parse().unwrap();
            let tree = StereoAtomDsl::from_edn(&dsl.to_edn()).unwrap();
            let stream = StereoAtomDsl::from_edn_str(&dsl.to_edn().to_string()).unwrap();
            prop_assert!(text.0.normalized_eq(&form));
            prop_assert!(tree.0.normalized_eq(&form));
            prop_assert!(stream.0.normalized_eq(&form));
        }
    }

    #[test]
    fn test_stereo_atom_update_dsl_roundtrip_kindless(
        payload in raw_stereo_atom_form_strategy(),
        clear in any::<bool>(),
    ) {
        let configuration = if clear { StereoConfigurationUpdate::Undetermined } else { StereoConfigurationUpdate::Unchanged };
        let dsl = StereoAtomUpdateDsl(StereoAtomUpdate { configuration: configuration.clone(), constraints: payload.constraints });
        let expected = StereoAtomUpdateDsl(StereoAtomUpdate { configuration, ..Default::default() });
        prop_assert_eq!(dsl.to_string(), if clear { "*" } else { "" });
        prop_assert_eq!(dsl.to_string().parse::<StereoAtomUpdateDsl>().unwrap(), expected.clone());
        prop_assert_eq!(StereoAtomUpdateDsl::from_edn(&dsl.to_edn()).unwrap(), expected);
    }

    #[test]
    fn test_stereo_atom_update_dsl_display_from_str_roundtrip(
        update in stereo_atom_update_strategy(),
    ) {
        let dsl = StereoAtomUpdateDsl(update);
        let rendered = dsl.to_string();
        let parsed: StereoAtomUpdateDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_stereo_atom_update_display_from_str_roundtrip(
        update in stereo_atom_update_strategy(),
    ) {
        let rendered = update.to_string();
        let parsed: StereoAtomUpdate = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(update, parsed);
    }

    #[test]
    fn test_stereo_atom_update_dsl_to_edn_from_edn_roundtrip(
        update in stereo_atom_update_strategy(),
    ) {
        let dsl = StereoAtomUpdateDsl(update);
        let edn = dsl.to_edn();
        let parsed = StereoAtomUpdateDsl::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("parse failed: {e}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_stereo_bond_dsl_display_from_str_roundtrip(
        stereo in stereo_bond_form_strategy(),
    ) {
        let dsl = StereoBondDsl(stereo);
        let rendered = dsl.to_string();
        let parsed: StereoBondDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_stereo_bond_dsl_roundtrip_raw(form in raw_stereo_bond_form_strategy()) {
        if form.clone().normalize().is_ok() {
            let dsl = StereoBondDsl(form.clone());
            let text: StereoBondDsl = dsl.to_string().parse().unwrap();
            let tree = StereoBondDsl::from_edn(&dsl.to_edn()).unwrap();
            let stream = StereoBondDsl::from_edn_str(&dsl.to_edn().to_string()).unwrap();
            prop_assert!(text.0.normalized_eq(&form));
            prop_assert!(tree.0.normalized_eq(&form));
            prop_assert!(stream.0.normalized_eq(&form));
        }
    }

    #[test]
    fn test_stereo_bond_update_dsl_roundtrip_kindless(
        payload in raw_stereo_bond_form_strategy(),
        clear in any::<bool>(),
    ) {
        let configuration = if clear { StereoConfigurationUpdate::Undetermined } else { StereoConfigurationUpdate::Unchanged };
        let dsl = StereoBondUpdateDsl(StereoBondUpdate { configuration: configuration.clone(), constraints: payload.constraints });
        let expected = StereoBondUpdateDsl(StereoBondUpdate { configuration, ..Default::default() });
        prop_assert_eq!(dsl.to_string(), if clear { "*" } else { "" });
        prop_assert_eq!(dsl.to_string().parse::<StereoBondUpdateDsl>().unwrap(), expected.clone());
        prop_assert_eq!(StereoBondUpdateDsl::from_edn(&dsl.to_edn()).unwrap(), expected);
    }

    #[test]
    fn test_stereo_bond_update_dsl_display_from_str_roundtrip(
        update in stereo_bond_update_strategy(),
    ) {
        let dsl = StereoBondUpdateDsl(update);
        let rendered = dsl.to_string();
        let parsed: StereoBondUpdateDsl = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_stereo_bond_update_display_from_str_roundtrip(
        update in stereo_bond_update_strategy(),
    ) {
        let rendered = update.to_string();
        let parsed: StereoBondUpdate = rendered.parse().map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}\nrendered: {rendered}"))
        })?;
        prop_assert_eq!(update, parsed);
    }

    #[test]
    fn test_stereo_bond_update_dsl_to_edn_from_edn_roundtrip(
        update in stereo_bond_update_strategy(),
    ) {
        let dsl = StereoBondUpdateDsl(update);
        let edn = dsl.to_edn();
        let parsed = StereoBondUpdateDsl::from_edn(&edn)
            .map_err(|e| TestCaseError::fail(format!("parse failed: {e}")))?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Canonical `Th0`/`Th1` render to the `:ccw`/`:cw` EDN keyword shorthand
    /// and parse back to the same form.
    #[test]
    fn test_stereo_atom_dsl_keyword_to_edn_from_edn_roundtrip(
        coset in prop_oneof![Just(0u32), Just(1u32)],
    ) {
        let dsl = StereoAtomDsl(StereoAtomForm::new(
            StereoKind::Tetrahedral,
            StereoCoset::Lit(coset),
        ));
        let edn = dsl.to_edn();
        prop_assert!(
            matches!(&edn, Edn::Keyword(_)),
            "expected keyword render for canonical stereo atom, got {edn:?}",
        );
        let parsed = StereoAtomDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Canonical `Ct0`/`Ct1` render to the `:z`/`:e` EDN keyword shorthand.
    #[test]
    fn test_stereo_bond_dsl_keyword_to_edn_from_edn_roundtrip(
        coset in prop_oneof![Just(0u32), Just(1u32)],
    ) {
        let dsl = StereoBondDsl(StereoBondForm::new(
            StereoKind::CisTrans,
            StereoCoset::Lit(coset),
        ));
        let edn = dsl.to_edn();
        prop_assert!(
            matches!(&edn, Edn::Keyword(_)),
            "expected keyword render for canonical stereo bond, got {edn:?}",
        );
        let parsed = StereoBondDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("parse failed: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    /// Molecule-scope (EDN-shaped) stereo atom constraint: `to_edn` → `from_edn`
    /// roundtrips for every kind and constraint variant.
    #[test]
    fn test_stereo_atom_constraint_dsl_to_edn_from_edn_roundtrip(
        args in stereo_atom_kind_strategy().prop_flat_map(|kind| {
            stereo_atom_constraint_strategy(kind).prop_map(move |c| (kind, c))
        }),
    ) {
        let (kind, constraint) = args;
        let dsl = StereoAtomConstraintDsl(kind, constraint);
        let edn = dsl.to_edn();
        let parsed = StereoAtomConstraintDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("from_edn failed for {edn:?}: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }

    #[test]
    fn test_stereo_bond_constraint_dsl_to_edn_from_edn_roundtrip(
        constraint in stereo_bond_constraint_strategy(StereoKind::CisTrans),
    ) {
        let dsl = StereoBondConstraintDsl(StereoKind::CisTrans, constraint);
        let edn = dsl.to_edn();
        let parsed = StereoBondConstraintDsl::from_edn(&edn).map_err(|e| {
            TestCaseError::fail(format!("from_edn failed for {edn:?}: {e}"))
        })?;
        prop_assert_eq!(dsl, parsed);
    }
}
