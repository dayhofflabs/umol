//! MOL parsing into a resolved [`Molecule`] — parse + raise + resolve.
//!
//! SMILES ingestion lives in [`crate::ingest`].

use umol_graph_ir::ir::{Molecule, TryIntoIr};
use umol_io::ctfile::config::CtfileIoConfig;
use umol_io::ctfile::parser::parse_mol_bytes_to_table_ir_with;
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;

use crate::ops::model::{ChemistryModel, ValenceModel};
use crate::ops::resolve::{ResolveConfig, ResolveUnderdetermined, Resolver};

/// Parse MOL to a resolved [`Molecule`] using default IO config and model.
pub fn parse_mol(input: &str) -> Result<Molecule, Box<dyn UmolError>> {
    parse_mol_bytes(input.as_bytes())
}

/// Parse MOL bytes to a resolved [`Molecule`] using default IO config and model.
pub fn parse_mol_bytes(input: &[u8]) -> Result<Molecule, Box<dyn UmolError>> {
    parse_mol_bytes_with(
        input,
        &CtfileIoConfig::basic(),
        &ChemistryModel {
            valence: ValenceModel::mdl(),
            ..ChemistryModel::default()
        },
        &ResolveConfig::default(),
    )
}

/// Parse MOL to a resolved [`Molecule`] with explicit IO, model, and resolve config.
pub fn parse_mol_with(
    input: &str,
    io_config: &CtfileIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<Molecule, Box<dyn UmolError>> {
    parse_mol_bytes_with(input.as_bytes(), io_config, model, resolve_config)
}

/// Parse MOL bytes to a resolved [`Molecule`] with explicit IO, model, and resolve config.
pub fn parse_mol_bytes_with(
    input: &[u8],
    io_config: &CtfileIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<Molecule, Box<dyn UmolError>> {
    let table_mol = parse_mol_bytes_to_table_ir_with(input, io_config)?;
    let molecule: Molecule = (&table_mol).try_into_ir(&())?;
    match Resolver::with_config(model, *resolve_config).resolve(molecule)? {
        Solution::Determined(molecule) => Ok(molecule),
        Solution::Underdetermined(()) => Err(Box::new(ResolveUnderdetermined)),
        Solution::Contradictory(c) => Err(Box::new(c)),
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{
        AromaticValenceForm, AtomId, BondId, RollbackMethod, StereoBondForm, StereoBondId,
        StereoConfigurationForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    };
    use umol_graph_ir::mol_dsl_concrete;
    use umol_io::ctfile::config::CtfileIoConfig;
    use umol_io::ctfile::parse_mol_to_ir;

    use super::{parse_mol_bytes, parse_mol_bytes_with, parse_mol_with};
    use crate::ops::model::{
        AromaticityModel, AromaticityRule, AromaticityTieBreak, ChemistryModel, ElementScope,
        RingLimits, StereoModel, ValenceModel, ValenceTieBreak,
    };
    use crate::ops::resolve::{
        ResolveConfig, ResolveUnderdetermined, Resolver, StereoFailurePolicy, StereoResolveConfig,
    };
    use crate::ops::valence::ValenceTable;
    use crate::ops::validate::ConnectivityModel;

    const METHANE_MOL: &str = "Methane\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    const BENZENE_AROMATIC_MOL: &str = "benzene\n\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  4  0  0  0  0\n  2  3  4  0  0  0  0\n  3  4  4  0  0  0  0\n  4  5  4  0  0  0  0\n  5  6  4  0  0  0  0\n  6  1  4  0  0  0  0\nM  END\n";

    #[fixture]
    fn valence_table() -> &'static ValenceTable {
        ValenceTable::default_table()
    }

    #[rstest]
    #[case::methane(METHANE_MOL, 1, "C#i=#c0#h4#n0#u0#s")]
    #[case::benzene(BENZENE_AROMATIC_MOL, 6, "C#i=#c0#h#n0#u0#s")]
    fn test_parse_mol_to_ir_counts_resolve(
        valence_table: &'static ValenceTable,
        #[case] input: &str,
        #[case] atom_count: u32,
        #[case] expected_atom: &str,
    ) {
        let mut molecule = parse_mol_to_ir(input).unwrap();
        let model = ChemistryModel {
            valence: ValenceModel {
                tie_break: ValenceTieBreak::MostSaturated,
                ..ValenceModel::counts(Cow::Borrowed(valence_table))
            },
            ..ChemistryModel::default()
        };
        Resolver::new(&model)
            .resolve_into(&mut molecule, RollbackMethod::Journal)
            .unwrap();
        assert_eq!(molecule.atoms().count(), atom_count as usize);
        for i in 0..atom_count {
            assert_eq!(
                molecule.atom(AtomId(i)).attributes().to_string(),
                expected_atom
            );
        }
    }

    #[rstest]
    #[case::methane(METHANE_MOL)]
    fn test_parse_mol_bytes(#[case] input: &str) {
        assert_eq!(
            parse_mol_bytes(input.as_bytes()).unwrap(),
            parse_mol_bytes_with(
                input.as_bytes(),
                &CtfileIoConfig::basic(),
                &ChemistryModel {
                    valence: ValenceModel::mdl(),
                    ..ChemistryModel::default()
                },
                &ResolveConfig::default(),
            )
            .unwrap()
        );
    }

    #[rstest]
    #[case::z(1.0, StereoCoset::Lit(0))]
    #[case::e(-1.0, StereoCoset::Lit(1))]
    fn test_parse_mol_bytes_with_stereo(
        #[case] ligand_y: f64,
        #[case] coset: StereoCoset,
        #[values(0, 8, 64)] stereo_bond_min_ring_size: u32,
        #[values(
            StereoFailurePolicy::Error,
            StereoFailurePolicy::Keep,
            StereoFailurePolicy::Remove
        )]
        failure: StereoFailurePolicy,
    ) {
        let input = format!(
            r#"cyclohexene
  umol           2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    2.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    2.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000{ligand_y:10.4}    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  6  1  1  0  0  0  0
M  END
"#
        );
        let original = input.clone();
        let raised = parse_mol_to_ir(&input).unwrap();
        let frame = vec![
            StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Virtual),
        ];
        assert_eq!(
            raised
                .stereo_bonds()
                .iter()
                .map(|view| (
                    view.site_id(),
                    view.ligand_ids().to_vec(),
                    view.attributes().clone(),
                ))
                .collect::<Vec<_>>(),
            vec![(
                BondId(0),
                frame,
                StereoBondForm::new(StereoKind::CisTrans, coset.clone())
            )]
        );
        let mut expected = mol_dsl_concrete!(
            r#"{
            :atoms ["C#i=#h" "C#i=#h" "C#i=#h2" "C#i=#h2" "C#i=#h2" "C#i=#h2"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :stereo-bonds [{:site 0 :ligands [5 [:h 0] 2 [:h 1]] :attrs "Ct0"}]
        }"#
        );
        expected
            .stereo_bond_mut(StereoBondId(0))
            .attributes_mut()
            .configuration = StereoConfigurationForm::kinded(StereoKind::CisTrans, coset);
        let model = ChemistryModel {
            valence: ValenceModel::mdl(),
            stereo: StereoModel {
                stereo_bond_min_ring_size,
                stereo_bond_e_min_ring_size: 64,
                ..Default::default()
            },
            ..Default::default()
        };
        let config = ResolveConfig {
            stereo: StereoResolveConfig {
                stereo_atom_failure: failure,
                stereo_bond_failure: failure,
                ..Default::default()
            },
            ..Default::default()
        };
        assert_eq!(
            parse_mol_bytes_with(input.as_bytes(), &CtfileIoConfig::basic(), &model, &config)
                .unwrap(),
            expected
        );
        assert_eq!(
            parse_mol_with(&input, &CtfileIoConfig::basic(), &model, &config).unwrap(),
            expected
        );
        assert_eq!(input, original);
    }

    #[rstest]
    #[case::counts(
        CtfileIoConfig::basic(),
        ChemistryModel {
            connectivity: ConnectivityModel::default(),
            valence: ValenceModel {
                tie_break: ValenceTieBreak::MostSaturated,
                ..ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table()))
            },
            aromaticity: AromaticityModel { scope: ElementScope::AllowList(vec![Element::C]), rule: AromaticityRule::Hueckel { ring_limits: RingLimits::default() }, tie_break: AromaticityTieBreak::Strict },
            stereo: StereoModel::default(),
        },
        ResolveConfig::default(),
        "C#i=#c0#h4#n0#u0#s"
    )]
    fn test_parse_mol_bytes_with(
        #[case] io_config: CtfileIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: &str,
    ) {
        let molecule =
            parse_mol_bytes_with(METHANE_MOL.as_bytes(), &io_config, &model, &resolve_config)
                .unwrap();
        assert_eq!(molecule.atom(AtomId(0)).attributes().to_string(), expected);
        assert_eq!(molecule.clone().edit().finish(), Ok(molecule));
    }

    #[rstest]
    #[case::counts(METHANE_MOL)]
    fn test_parse_mol_bytes_with_underdetermined(#[case] input: &str) {
        let model = ChemistryModel {
            valence: ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table())),
            ..ChemistryModel::default()
        };
        let error = parse_mol_bytes_with(
            input.as_bytes(),
            &CtfileIoConfig::basic(),
            &model,
            &ResolveConfig::default(),
        )
        .unwrap_err();
        assert_eq!(
            error.as_any().downcast_ref::<ResolveUnderdetermined>(),
            Some(&ResolveUnderdetermined)
        );
    }

    #[rstest]
    #[case::default(
        CtfileIoConfig::basic(),
        ChemistryModel::default(),
        ResolveConfig::default(),
        vec![AromaticValenceForm::Undetermined; 6]
    )]
    fn test_parse_mol_bytes_with_aromatic_valence(
        #[case] io_config: CtfileIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: Vec<AromaticValenceForm>,
    ) {
        let molecule = parse_mol_bytes_with(
            BENZENE_AROMATIC_MOL.as_bytes(),
            &io_config,
            &model,
            &resolve_config,
        )
        .unwrap();

        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.attributes().constraints.aromatic_valence().clone())
                .collect::<Vec<_>>(),
            expected
        );
        assert_eq!(molecule.clone().edit().finish(), Ok(molecule));
    }
}
