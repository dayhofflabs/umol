//! MOL parsing into a resolved [`Molecule`] — parse + raise + resolve.
//!
//! SMILES ingestion lives in [`crate::ingest`].

use umol_graph_ir::ir::{Molecule, TryIntoIr};
use umol_io::ctfile::config::CtfileIoConfig;
use umol_io::ctfile::parser::parse_mol_bytes_to_table_ir_with;
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;

use crate::ops::model::ChemistryModel;
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
        &ChemistryModel::default(),
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
    let mut molecule: Molecule = (&table_mol).try_into_ir(&())?;
    match Resolver::with_config(model, *resolve_config).resolve(&mut molecule)? {
        Solution::Determined(()) => Ok(molecule),
        Solution::Underdetermined(()) => Err(Box::new(ResolveUnderdetermined)),
        Solution::Contradictory(c) => Err(Box::new(c)),
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{AromaticValenceForm, AtomId, NumForm};
    use umol_io::ctfile::config::CtfileIoConfig;
    use umol_io::ctfile::parse_mol_to_ast;

    use super::{parse_mol_bytes, parse_mol_bytes_with};
    use crate::ops::model::{
        AromaticityModel, ChemistryModel, ElementScope, RingLimits, StereoModel, ValenceModel,
    };
    use crate::ops::resolve::{AromaticityResolveConfig, ResolveConfig, StereoResolveConfig};
    use crate::ops::valence::{CountsValence, ValenceTable};

    const METHANE_MOL: &str = "Methane\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    const BENZENE_AROMATIC_MOL: &str = "benzene\n\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  4  0  0  0  0\n  2  3  4  0  0  0  0\n  3  4  4  0  0  0  0\n  4  5  4  0  0  0  0\n  5  6  4  0  0  0  0\n  6  1  4  0  0  0  0\nM  END\n";

    #[fixture]
    fn valence_table() -> &'static ValenceTable {
        ValenceTable::default_table()
    }

    #[rstest]
    #[case::methane(METHANE_MOL, 1, "C#i=#c0#h4#n0#u0#s#v0#a!")]
    #[case::benzene(BENZENE_AROMATIC_MOL, 6, "C#i=#c0#h#n0#u0#s#v2#a")]
    fn test_parse_mol_to_ast_counts_resolve(
        valence_table: &'static ValenceTable,
        #[case] input: &str,
        #[case] atom_count: u32,
        #[case] expected_atom: &str,
    ) {
        let mut molecule = parse_mol_to_ast(input).unwrap();
        CountsValence::new(valence_table)
            .resolve(&mut molecule)
            .unwrap();
        assert_eq!(molecule.atoms().count(), atom_count as usize);
        for i in 0..atom_count {
            assert_eq!(
                molecule.atom(AtomId(i)).attributes.to_string(),
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
                &ChemistryModel::default(),
                &ResolveConfig::default(),
            )
            .unwrap()
        );
    }

    #[rstest]
    #[case::counts(
        CtfileIoConfig::basic(),
        ChemistryModel {
            valence: ValenceModel::Counts {
                table: Cow::Borrowed(ValenceTable::default_table()),
            },
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
            stereo: StereoModel::default(),
        },
        ResolveConfig::default(),
        "C#i=#c0#h4#n0#u0#s#v0#a!"
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
        assert_eq!(molecule.atom(AtomId(0)).attributes.to_string(), expected);
    }

    #[rstest]
    #[case::retained(
        CtfileIoConfig::basic(),
        ChemistryModel::default(),
        ResolveConfig::default(),
        vec![Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))); 6]
    )]
    #[case::reset(
        CtfileIoConfig::basic(),
        ChemistryModel::default(),
        ResolveConfig {
            aromaticity: AromaticityResolveConfig {
                reset_aromatic_valence: true,
                ..AromaticityResolveConfig::default()
            },
            stereo: StereoResolveConfig::default(),
        },
        vec![None; 6]
    )]
    fn test_parse_mol_bytes_with_aromatic_valence(
        #[case] io_config: CtfileIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: Vec<Option<AromaticValenceForm>>,
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
                .map(|atom| atom.attributes.constraints.aromatic_valence().cloned())
                .collect::<Vec<_>>(),
            expected
        );
    }
}
