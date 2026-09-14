//! Conversion of graph models into external-format boundary values.

use std::any::Any;

use thiserror::Error;
use umol_chem::spin::SpinMultiplicity;
use umol_graph_ir::ir::{
    AromaticValenceForm, AsLit, AtomConstraintForm, AtomId, BondConstraintForm, BondId,
    BooleanForm, CisTransStereoForm, Constraint, ElementForm, Entity, IsotopeMassForm, Lattice,
    Molecule, NoncovalentBondKind, NoncovalentBondKindForm, NumForm, StereoCoset,
    TetrahedralStereoForm,
};
use umol_io::smiles::{Smiles, SmilesIoConfig};
use umol_io::table_ir::{
    Atom, Bond, BondConfiguration, BondDonation, BondNoncovalent, BondOrder, BondRelation,
    Molecule as TableMolecule, StereoAtom, StereoBond, StereoLigand, Winding,
};
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;

use crate::ops::resolve::{ProjectContradiction, ProjectError, ProjectFlags, Resolver};

/// Convert a graph model into an external-format boundary value.
pub trait Convey: Sized {
    type Input;
    type Config;
    type Error;

    /// Project a private copy and convert its fields into the boundary representation.
    fn convey(
        input: &Self::Input,
        resolver: &Resolver<'_>,
        config: &Self::Config,
    ) -> Result<Self, Self::Error>;
}

/// Failure to project a molecule or express its projected values in TableIR.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ConveyError {
    #[error(transparent)]
    Projection(#[from] ProjectError),
    #[error(transparent)]
    Contradiction(#[from] ProjectContradiction),
    #[error("projection is underdetermined")]
    Underdetermined,
    #[error("{entity:?} {field} cannot be represented in TableIR: {value}")]
    Value {
        entity: Entity,
        field: &'static str,
        value: String,
    },
    #[error("constraint cannot be represented in TableIR: {0:?}")]
    Constraint(Constraint),
    #[error("{entity:?} cannot be represented in TableIR")]
    Entity { entity: Entity },
    #[error("tetrahedral frame cannot be constructed at atom {atom:?}")]
    StereoAtom { atom: AtomId },
    #[error("stereo references cannot be constructed at bond {bond:?}")]
    StereoBond { bond: BondId },
}

impl UmolError for ConveyError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl Convey for Smiles {
    type Input = Molecule;
    type Config = SmilesIoConfig;
    type Error = ConveyError;

    /// Convert projected molecular fields and stereo assertions into a SMILES boundary value.
    ///
    /// H counts are copied after projection: a literal becomes Some, including zero, and
    /// Undetermined becomes None. Atom electron fields remain in the table. Format syntax
    /// options are consumed by `render_with`; convey performs no rendering or H inference.
    ///
    /// # Semantic properties
    ///
    /// Success and failure leave the input unchanged. Atom and localized-bond order are
    /// preserved. Actual H atoms, implicit H, and lone-pair stereo ligands remain distinct.
    /// Supported ingested structures survive convey, render, parse, and interpretation under
    /// the same resolver policy, up to molecular equivalence and permitted notation changes.
    ///
    /// # Errors
    ///
    /// Reports projection failures and projected values or constraints that TableIR cannot
    /// encode. A successfully constructed boundary may still fail its format's rendering rules.
    fn convey(
        input: &Molecule,
        resolver: &Resolver<'_>,
        _config: &SmilesIoConfig,
    ) -> Result<Self, ConveyError> {
        convey_molecule(input, resolver, ProjectFlags::all()).map(Self::from_table_ir)
    }
}

fn convey_molecule(
    input: &Molecule,
    resolver: &Resolver<'_>,
    flags: ProjectFlags,
) -> Result<TableMolecule, ConveyError> {
    let mut projected = input.clone();
    match resolver.project(&mut projected, flags)? {
        Solution::Determined(()) => {}
        Solution::Underdetermined(()) => return Err(ConveyError::Underdetermined),
        Solution::Contradictory(error) => return Err(error.into()),
    }
    let molecule = &projected;
    let mut table = TableMolecule::empty();
    table.atoms.reserve(molecule.atoms().count());
    table.bonds.reserve(molecule.bonds().count());
    for atom in molecule.atoms().iter() {
        let entity = Entity::Atom(atom.id);
        let form = atom.attributes;
        let mut lowered = match &form.element {
            ElementForm::Undetermined => Atom::wildcard(),
            ElementForm::Lit(element) => Atom::from_element(*element),
            value => {
                return Err(ConveyError::Value {
                    entity,
                    field: "element",
                    value: format!("{value:?}"),
                })
            }
        };
        lowered.isotope_mass = match &form.isotope_mass {
            IsotopeMassForm::Undetermined | IsotopeMassForm::Natural => None,
            IsotopeMassForm::Lit(mass) => Some(*mass),
            value => {
                return Err(ConveyError::Value {
                    entity,
                    field: "isotope",
                    value: format!("{value:?}"),
                })
            }
        };
        lowered.charge = lower_number(&form.charge, entity, "charge")?;
        lowered.implicit_hydrogens =
            lower_number(&form.implicit_hydrogens, entity, "implicit hydrogens")?;
        lowered.lone_pairs = lower_number(&form.lone_pairs, entity, "lone pairs")?;
        lowered.unpaired_electrons =
            lower_number(&form.unpaired_electrons.count, entity, "unpaired electrons")?;
        lowered.multiplicity = lower_multiplicity(&form.unpaired_electrons.multiplicity, entity)?;
        for constraint in form.constraints.iter() {
            match constraint {
                AtomConstraintForm::Valence(value) => {
                    lowered.valence = lower_number(value, entity, "valence")?
                }
                AtomConstraintForm::AromaticValence(value) => {
                    lowered.aromatic = match value {
                        AromaticValenceForm::Undetermined => None,
                        AromaticValenceForm::NotAromatic => Some(false),
                        AromaticValenceForm::Aromatic(_) => Some(true),
                    }
                }
                AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(coset)) => {
                    table
                        .stereo_atoms
                        .push(lower_stereo_atom(molecule, atom.id, coset)?);
                }
                value if value.is_undetermined() => {}
                value => {
                    return Err(ConveyError::Constraint(Constraint::Atom(
                        atom.id,
                        value.clone(),
                    )))
                }
            }
        }
        table.atoms.push(lowered);
    }
    for bond in molecule.bonds().iter() {
        let entity = Entity::Bond(bond.id);
        let [first, second] = bond.atom_ids();
        let mut lowered = Bond::new(first.0, second.0, lower_order(bond.order(), entity)?);
        lowered.charge = lower_number(&bond.attributes.charge, entity, "charge")?;
        lowered.unpaired_electrons = lower_number(
            &bond.attributes.unpaired_electrons.count,
            entity,
            "unpaired electrons",
        )?;
        lowered.multiplicity =
            lower_multiplicity(&bond.attributes.unpaired_electrons.multiplicity, entity)?;
        for constraint in bond.attributes.constraints.iter() {
            match constraint {
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)) => {
                    lowered.order = BondOrder::Aromatic
                }
                BondConstraintForm::Aromatic(BooleanForm::Lit(false)) => {}
                BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(coset)) => {
                    table
                        .stereo_bonds
                        .push(lower_stereo_bond(molecule, bond.id, coset)?);
                }
                value if value.is_undetermined() => {}
                value => {
                    return Err(ConveyError::Constraint(Constraint::Bond(
                        bond.id,
                        value.clone(),
                    )))
                }
            }
        }
        table.bonds.push(lowered);
    }
    for bond in molecule.dative_bonds().iter() {
        let entity = Entity::DativeBond(bond.id);
        let mut donors = bond.donor_ids();
        let Some(donor) = donors.next() else {
            return Err(ConveyError::Entity { entity });
        };
        if donors.next().is_some() {
            return Err(ConveyError::Entity { entity });
        }
        for constraint in bond.attributes.constraints.iter() {
            if !constraint.is_undetermined() {
                return Err(ConveyError::Constraint(Constraint::DativeBond(
                    bond.id,
                    constraint.clone(),
                )));
            }
        }
        table.bonds.push(Bond::new_dative(
            donor.0,
            bond.acceptor_id().0,
            lower_order(&bond.attributes.order, entity)?,
            BondDonation::Donating,
        ));
    }
    if let Some(bond) = molecule.multicenter_bonds().iter().next() {
        return Err(ConveyError::Entity {
            entity: Entity::MulticenterBond(bond.id),
        });
    }
    for bond in molecule.noncovalent_bonds().iter() {
        let entity = Entity::NoncovalentBond(bond.id);
        if bond.kind() != &NoncovalentBondKindForm::Lit(NoncovalentBondKind::HydrogenBond) {
            return Err(ConveyError::Entity { entity });
        }
        for constraint in bond.attributes.constraints.iter() {
            if !constraint.is_undetermined() {
                return Err(ConveyError::Constraint(Constraint::NoncovalentBond(
                    bond.id,
                    constraint.clone(),
                )));
            }
        }
        let [first, second] = bond.atom_ids();
        table.bonds.push(Bond::new_noncovalent(
            first.0,
            second.0,
            BondNoncovalent::Hydrogen,
        ));
    }
    if let Some(constraint) = molecule.constraints().iter().next() {
        return Err(ConveyError::Constraint(constraint.clone()));
    }
    Ok(table)
}

fn lower_number<T: TryFrom<i64>>(
    form: &NumForm,
    entity: Entity,
    field: &'static str,
) -> Result<Option<T>, ConveyError> {
    match form {
        NumForm::Undetermined => Ok(None),
        _ => form
            .as_lit()
            .and_then(|value| T::try_from(value).ok())
            .map(Some)
            .ok_or_else(|| ConveyError::Value {
                entity,
                field,
                value: format!("{form:?}"),
            }),
    }
}

fn lower_multiplicity(
    form: &NumForm,
    entity: Entity,
) -> Result<Option<SpinMultiplicity>, ConveyError> {
    lower_number::<u8>(form, entity, "multiplicity")?
        .map(|value| {
            SpinMultiplicity::new(value).ok_or_else(|| ConveyError::Value {
                entity,
                field: "multiplicity",
                value: format!("{form:?}"),
            })
        })
        .transpose()
}

fn lower_order(form: &NumForm, entity: Entity) -> Result<BondOrder, ConveyError> {
    if matches!(form, NumForm::Undetermined) {
        return Ok(BondOrder::Any);
    }
    if let NumForm::LitSet(values) = form {
        if values.iter().copied().eq([1, 2]) {
            return Ok(BondOrder::SingleOrDouble);
        }
    }
    match form.as_lit() {
        Some(0) => Ok(BondOrder::Zero),
        Some(1) => Ok(BondOrder::Single),
        Some(2) => Ok(BondOrder::Double),
        Some(3) => Ok(BondOrder::Triple),
        Some(4) => Ok(BondOrder::Quadruple),
        Some(5) => Ok(BondOrder::Quintuple),
        Some(6) => Ok(BondOrder::Sextuple),
        _ => Err(ConveyError::Value {
            entity,
            field: "order",
            value: format!("{form:?}"),
        }),
    }
}

fn lower_stereo_atom(
    molecule: &Molecule,
    atom: AtomId,
    coset: &StereoCoset,
) -> Result<StereoAtom, ConveyError> {
    let winding = match coset.as_lit() {
        Some(0) => Winding::CounterClockwise,
        Some(1) => Winding::Clockwise,
        _ => return Err(ConveyError::StereoAtom { atom }),
    };
    let view = molecule.atom(atom);
    let mut ligands: Vec<_> = view
        .neighbors()
        .map(|neighbor| StereoLigand::Atom(neighbor.atom_id().0))
        .collect();
    if ligands.len() == 3 {
        ligands.push(match view.implicit_hydrogens().as_lit() {
            Some(1) => StereoLigand::ImplicitHydrogen,
            Some(0) if view.lone_pairs().as_lit().is_some_and(|count| count > 0) => {
                StereoLigand::LonePair
            }
            _ => return Err(ConveyError::StereoAtom { atom }),
        });
    }
    if ligands.len() != 4 {
        return Err(ConveyError::StereoAtom { atom });
    }
    Ok(StereoAtom {
        atom: atom.0,
        ligands,
        winding,
    })
}

fn lower_stereo_bond(
    molecule: &Molecule,
    bond: BondId,
    coset: &StereoCoset,
) -> Result<StereoBond, ConveyError> {
    let configuration = if matches!(coset, StereoCoset::Undetermined) {
        BondConfiguration::Either
    } else {
        let relation = match coset.as_lit() {
            Some(0) => BondRelation::SameSide,
            Some(1) => BondRelation::OppositeSide,
            _ => return Err(ConveyError::StereoBond { bond }),
        };
        let [first, second] = molecule.bond(bond).atom_ids();
        let mut references = [0; 2];
        for (side, (atom, partner)) in [(first, second), (second, first)].into_iter().enumerate() {
            references[side] = molecule
                .atom(atom)
                .neighbors()
                .map(|neighbor| neighbor.atom_id())
                .find(|&neighbor| neighbor != partner)
                .ok_or(ConveyError::StereoBond { bond })?
                .0;
        }
        BondConfiguration::Framed {
            references,
            relation,
        }
    };
    Ok(StereoBond {
        bond: bond.0,
        configuration,
    })
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element;
    use umol_graph_core::AutomorphismAlgorithm;
    use umol_graph_ir::ir::{Canonicalize, CanonicalizeContext};
    use umol_graph_ir::mol_dsl_concrete;

    use super::*;
    use crate::ingest::{ingest_smiles, ingest_smiles_with};
    use crate::ops::model::{ChemistryModel, ValenceModel};
    use crate::ops::resolve::{IsotopePolicy, ResolveConfig};

    #[rstest]
    #[case::ordinary("C", None, None, 0, SpinMultiplicity::SINGLET)]
    #[case::radical("[CH3]", None, Some(3), 1, SpinMultiplicity::DOUBLET)]
    #[case::isotope("[13CH4]", Some(13), Some(4), 0, SpinMultiplicity::SINGLET)]
    fn test_smiles_convey(
        #[case] input: &str,
        #[case] isotope_mass: Option<u32>,
        #[case] implicit_hydrogens: Option<u8>,
        #[case] unpaired_electrons: u8,
        #[case] multiplicity: SpinMultiplicity,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..Default::default()
        };
        let resolver = Resolver::with_config(
            &model,
            ResolveConfig {
                isotope: IsotopePolicy::Natural,
                ..Default::default()
            },
        );
        let source = ingest_smiles(input).unwrap();
        let mut expected = TableMolecule::empty();
        expected.atoms.push(Atom {
            isotope_mass,
            implicit_hydrogens,
            charge: Some(0),
            lone_pairs: Some(0),
            unpaired_electrons: Some(unpaired_electrons),
            multiplicity: Some(multiplicity),
            ..Atom::from_element(Element::C)
        });
        assert_eq!(
            Smiles::convey(&source, &resolver, &SmilesIoConfig::opensmiles())
                .unwrap()
                .into_table_ir(),
            expected
        );
    }

    #[rstest]
    #[case::empty("", "")]
    #[case::chain("CCO", "CCO")]
    #[case::branched("CC(C)O", "CC(C)O")]
    #[case::aromatic("c1ccccc1", "[cH]1[cH][cH][cH][cH][cH]1")]
    #[case::heteroaromatic("n1ccccc1", "n1[cH][cH][cH][cH][cH]1")]
    #[case::pyrrole("[nH]1cccc1", "[nH]1[cH][cH][cH][cH]1")]
    #[case::biphenyl("c1ccccc1-c2ccccc2", "[cH]1[cH][cH][cH][cH]c1-c1[cH][cH][cH][cH][cH]1")]
    #[case::radical("[CH3]", "[CH3]")]
    #[case::charge("[NH4+]", "[NH4+]")]
    #[case::isotope("[13CH4]", "[13CH4]")]
    #[case::tetrahedral_h("F[C@H](Cl)Br", "F[C@H](Cl)Br")]
    #[case::tetrahedral_atoms("F[C@](Cl)(Br)I", "F[C@](Cl)(Br)I")]
    #[case::lone_pair("C[S@](=O)CC", "C[S@](=O)CC")]
    #[case::explicit_h("[H][C@](F)(Cl)Br", "[H][C@](F)(Cl)Br")]
    #[case::trans("F/C=C/Cl", "F/C=C/Cl")]
    #[case::cis("F/C=C\\Cl", "F/C=C\\Cl")]
    #[case::branched_stereo("F/C(Cl)=C(/Br)I", "FC(/Cl)=C(Br)\\I")]
    #[case::conjugated("C/C=C/C=C/C", "C/C=C/C=C/C")]
    fn test_smiles_convey_roundtrip(#[case] input: &str, #[case] expected: &str) {
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..Default::default()
        };
        let config = ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        };
        let io = SmilesIoConfig::opensmiles();
        let resolver = Resolver::with_config(&model, config);
        let original = ingest_smiles(input).unwrap();
        let source = original.clone();
        let smiles = Smiles::convey(&source, &resolver, &io).unwrap();
        assert_eq!(source, original);
        let text = smiles.render().unwrap();
        assert_eq!(text, expected);
        let restored = ingest_smiles_with(&text, &io, &model, &config).unwrap();
        assert!(original.canonical_eq(
            &restored,
            &CanonicalizeContext {
                para_stereo: false,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        ));
    }

    #[rstest]
    #[case::either("2#C+", BondConfiguration::Either)]
    #[case::same("2#C0", BondConfiguration::Framed { references: [0, 3], relation: BondRelation::SameSide })]
    #[case::opposite("2#C1", BondConfiguration::Framed { references: [0, 3], relation: BondRelation::OppositeSide })]
    fn test_smiles_convey_stereo_bond(
        #[case] bond: &str,
        #[case] configuration: BondConfiguration,
    ) {
        let source = mol_dsl_concrete!(&format!(
            r#"{{:atoms ["F#n3" "C#h" "C#h" "F#n3"] :bonds [[0 1 "1"] [1 2 "{bond}"] [2 3 "1"]]}}"#
        ));
        let model = ChemistryModel::default();
        let table = Smiles::convey(
            &source,
            &Resolver::new(&model),
            &SmilesIoConfig::opensmiles(),
        )
        .unwrap()
        .into_table_ir();
        assert_eq!(
            table.stereo_bonds,
            vec![StereoBond {
                bond: 1,
                configuration
            }]
        );
    }

    #[rstest]
    #[case::hydrogens("C#i13#h256", "implicit hydrogens", "Lit(256)")]
    #[case::lone_pairs("C#n256", "lone pairs", "Lit(256)")]
    #[case::charge("C#c128", "charge", "Lit(128)")]
    #[case::unpaired("C#u256#s", "unpaired electrons", "Lit(256)")]
    #[case::multiplicity("C#s256", "multiplicity", "Lit(256)")]
    #[case::zero_multiplicity("C#s0", "multiplicity", "Lit(0)")]
    fn test_smiles_convey_error(
        #[case] atom: &str,
        #[case] field: &'static str,
        #[case] value: &str,
    ) {
        let source = mol_dsl_concrete!(&format!("{{:atoms [\"{atom}\"]}}"));
        let original = source.clone();
        let model = ChemistryModel::default();
        assert_eq!(
            Smiles::convey(
                &source,
                &Resolver::new(&model),
                &SmilesIoConfig::opensmiles()
            ),
            Err(ConveyError::Value {
                entity: Entity::Atom(AtomId(0)),
                field,
                value: value.to_owned()
            })
        );
        assert_eq!(source, original);
    }

    #[rstest]
    #[case::bond_charge(
        mol_dsl_concrete!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+"]]}"#),
        ConveyError::Projection(ProjectError::BondCharge { entity: Entity::Bond(BondId(0)), charge: NumForm::Lit(1) })
    )]
    #[case::atom_constraint(
        mol_dsl_concrete!(r#"{:atoms ["C#d3"]}"#),
        ConveyError::Constraint(Constraint::Atom(AtomId(0), AtomConstraintForm::DonatedPairs(NumForm::Lit(3))))
    )]
    fn test_smiles_convey_projection_error(
        #[case] source: Molecule,
        #[case] expected: ConveyError,
    ) {
        let original = source.clone();
        let model = ChemistryModel::default();
        assert_eq!(
            Smiles::convey(
                &source,
                &Resolver::new(&model),
                &SmilesIoConfig::opensmiles()
            ),
            Err(expected)
        );
        assert_eq!(source, original);
    }
}
