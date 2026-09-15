//! Conversion of graph models into external-format boundary values and text.
//!
//! [`Convey`] constructs a resolver from the supplied chemistry model and resolve configuration,
//! projects a private GraphIR copy, and converts it to the output boundary. The boundary's renderer
//! then formats its table. [`export_smiles`] and [`export_reaction_smiles`] compose those operations
//! with the same OpenSMILES, SMILES-valence, and Natural-isotope defaults as ingestion. Their
//! `_with` variants take `(input, io_config, model, resolve_config)`, matching ingestion.
//!
//! Supported output includes ordinary atoms, isotopes and charges, radicals encoded through
//! bracket H counts, neutral closed-shell aromatic systems, tetrahedral stereo, and local
//! double-bond stereo with a consistent slash assignment. Projection preserves H counts; Convey
//! can omit a table H count when notation permits and inference under the supplied model and
//! policy reproduces it. Rendering performs no chemical resolution.
//!
//! Projection rejects nonzero bond/aromatic-system charge and non-singlet or unpaired spin there.
//! Conversion rejects values outside TableIR capacity and unsupported entities or constraints.
//! Rendering owns narrower syntax limits, including Either, cumulene stereo, unsupported bond
//! orders, and configurations with no slash assignment. Atom lone-pair counts remain in TableIR
//! but have no independent ordinary-SMILES token. No coordinates are generated.
//!
//! # Examples
//!
//! Default text ingestion and export:
//!
//! ```
//! use umol_graph::export::{export_reaction_smiles, export_smiles};
//! use umol_graph::ingest::{ingest_reaction_smiles, ingest_smiles};
//!
//! let molecule = ingest_smiles("F/C=C/F")?;
//! assert_eq!(export_smiles(&molecule)?, "F/C=C/F");
//! let reaction = ingest_reaction_smiles("[CH4:9]>>[CH4:9]")?;
//! assert_eq!(export_reaction_smiles(&reaction)?, "[CH4:1]>>[CH4:1]");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Explicit parsing, interpretation, conversion, and rendering use the same configurations:
//!
//! ```
//! use umol_graph::export::{export_smiles_with, Convey};
//! use umol_graph::ingest::{ingest_smiles_with, Interpret};
//! use umol_graph::ops::model::{ChemistryModel, ValenceModel};
//! use umol_graph::ops::resolve::{IsotopePolicy, ResolveConfig};
//! use umol_io::smiles::{Smiles, SmilesIoConfig};
//!
//! let io_config = SmilesIoConfig::opensmiles();
//! let model = ChemistryModel { valence: ValenceModel::smiles(), ..Default::default() };
//! let resolve_config = ResolveConfig { isotope: IsotopePolicy::Natural, ..Default::default() };
//! let parsed = Smiles::parse_with("F[C@H](Cl)Br", &io_config)?;
//! let molecule = parsed.interpret(&model, &resolve_config)?;
//! let output = Smiles::convey(&molecule, &model, &resolve_config, &io_config)?;
//! let text = output.render_with(&io_config)?;
//! assert_eq!(text, "F[C@H](Cl)Br");
//! assert_eq!(
//!     ingest_smiles_with("F[C@H](Cl)Br", &io_config, &model, &resolve_config)?,
//!     molecule,
//! );
//! assert_eq!(export_smiles_with(&molecule, &io_config, &model, &resolve_config)?, text);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use std::any::Any;
use std::iter;

use thiserror::Error;
use umol_chem::element::Element;
use umol_chem::spin::SpinMultiplicity;
use umol_graph_ir::ir::{
    AromaticValenceForm, AsLit, AtomConstraintForm, AtomId, BondConstraintForm, BondId,
    BooleanForm, CisTransStereoForm, Constraint, Contradiction, ElementForm, Entity,
    IsotopeMassForm, Lattice, Molecule, NoncovalentBondKind, NoncovalentBondKindForm, NumForm,
    Reaction, StereoCoset, TetrahedralStereoForm,
};
use umol_io::smiles::{
    ReactionSmiles, ReactionSmilesRenderError, Smiles, SmilesIoConfig, SmilesRenderError,
};
use umol_io::table_ir::{
    Atom, Bond, BondConfiguration, BondDonation, BondNoncovalent, BondOrder, BondRelation,
    Molecule as TableMolecule, Reaction as TableReaction, StereoAtom, StereoBond, StereoLigand,
    Winding,
};
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;

use crate::ops::model::{ChemistryModel, ValenceModel};
use crate::ops::resolve::{
    IsotopePolicy, ProjectContradiction, ProjectError, ProjectFlags, ResolveConfig, Resolver,
};

/// Convert a graph model into an external-format boundary value.
pub trait Convey: Sized {
    type Input;
    type Config;
    type Error;

    /// Project a private copy and convert its fields into the boundary representation.
    fn convey(
        input: &Self::Input,
        model: &ChemistryModel,
        resolve_config: &ResolveConfig,
        io_config: &Self::Config,
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

/// Failure to materialize a reaction or convey either molecular side.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionConveyError {
    #[error("reaction cannot be materialized: {0}")]
    Materialization(#[from] Contradiction),
    #[error("reactants: {0}")]
    Reactants(#[source] ConveyError),
    #[error("products: {0}")]
    Products(#[source] ConveyError),
    #[error("atom correspondence pair {index} has no representable one-based map label")]
    AtomMapLabel { index: usize },
}

impl UmolError for ReactionConveyError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// Failure to convey or render a molecule as SMILES text.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum SmilesOutputError {
    #[error("{0}")]
    Convey(#[from] ConveyError),
    #[error("{0}")]
    Render(#[from] SmilesRenderError),
}

impl UmolError for SmilesOutputError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// Failure to convey or render a reaction as SMILES text, retaining side context.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionSmilesOutputError {
    #[error("{0}")]
    Convey(#[from] ReactionConveyError),
    #[error("{0}")]
    Render(#[from] ReactionSmilesRenderError),
}

impl UmolError for ReactionSmilesOutputError {
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
    /// Constructs a resolver from the model and resolve configuration. Projection preserves H;
    /// conversion omits a count only when brackets are unnecessary and inference under the
    /// selected valence policy reproduces it. Otherwise it retains the count, including zero.
    /// Undetermined H remains absent. Atom electron fields remain in the table. Format syntax
    /// options are consumed by `render_with`.
    ///
    /// # Semantic properties
    ///
    /// Success and failure leave the input unchanged. Atom and localized-bond order are
    /// preserved. Actual H atoms, implicit H, and lone-pair stereo ligands remain distinct.
    /// For the supported ingested domain, convey, render, parse, and interpretation under
    /// the same model and policy recover a molecule equal under Molecule::canonical_eq with
    /// para_stereo=false. Source spelling, ring labels, and redundant slash markers may change.
    /// H omission uses joint valence/aromatic evidence: Strict requires one allowed count,
    /// MostSaturated selects the greatest. Missing or ambiguous evidence retains the count.
    /// Stored lone pairs and atom-typing row order or repetition do not change that decision.
    ///
    /// # Errors
    ///
    /// Reports projection failures and projected values or constraints that TableIR cannot
    /// encode. A successfully constructed boundary may still fail its format's rendering rules.
    fn convey(
        input: &Molecule,
        model: &ChemistryModel,
        resolve_config: &ResolveConfig,
        _io_config: &SmilesIoConfig,
    ) -> Result<Self, ConveyError> {
        let resolver = Resolver::with_config(model, *resolve_config);
        convey_molecule(input, &resolver, iter::empty()).map(Self::from_table_ir)
    }
}

impl Convey for ReactionSmiles {
    type Input = Reaction;
    type Config = SmilesIoConfig;
    type Error = ReactionConveyError;

    /// Materialize and project both reaction sides, retaining H counts for atom-map labels.
    ///
    /// Constructs one resolver for both sides and runs all projection stages on private copies.
    /// Surviving atom pairs receive one-based labels in the materialized correspondence's order.
    /// Unmatched atoms remain unlabeled. The atom_mapping index is populated alongside the
    /// Atom.class labels.
    /// Agents and source-only labels or metadata are not available from graph IR.
    ///
    /// # Semantic properties
    ///
    /// Success and failure leave the source unchanged. The boundary preserves supported side
    /// semantics and atom correspondence through side compaction, including creation and deletion.
    /// Mapped atoms retain their H counts. Unmatched atoms use the same omission decision as
    /// molecular conversion. Repeated convey produces the same boundary. For the supported
    /// ingested domain, convey/render followed by parse/interpret under the same model and policy
    /// recovers a reaction equal under Reaction::canonical_eq with para_stereo=false, including
    /// correspondence. Rendering still owns format-specific support failures.
    ///
    /// # Errors
    ///
    /// Reports materialization failures, side-specific projection or conversion failures, and
    /// map labels that exceed TableIR's capacity.
    fn convey(
        input: &Reaction,
        model: &ChemistryModel,
        resolve_config: &ResolveConfig,
        _io_config: &SmilesIoConfig,
    ) -> Result<Self, ReactionConveyError> {
        let span = input.to_reaction_span()?;
        let correspondence = span.correspondence();
        let pairs = correspondence.atoms().matched_pairs();
        let count = u32::try_from(pairs.len()).map_err(|_| ReactionConveyError::AtomMapLabel {
            index: u32::MAX as usize,
        })?;
        let labels = 1..=count;
        let resolver = Resolver::with_config(model, *resolve_config);
        let reactants = convey_molecule(
            &span.lhs(),
            &resolver,
            pairs
                .iter()
                .zip(labels.clone())
                .map(|(&(left, _), label)| (left, label)),
        )
        .map_err(ReactionConveyError::Reactants)?;
        let products = convey_molecule(
            &span.rhs(),
            &resolver,
            pairs
                .iter()
                .zip(labels.clone())
                .map(|(&(_, right), label)| (right, label)),
        )
        .map_err(ReactionConveyError::Products)?;
        let mut table = TableReaction::from_molecules(reactants, products, TableMolecule::empty());
        for (&(left, right), class) in pairs.iter().zip(labels) {
            table
                .atom_mapping
                .insert(class, (vec![left.0], vec![right.0]));
        }
        Ok(Self::from_table_ir(table))
    }
}

/// Export a molecule with the OpenSMILES configuration, SMILES valence preset,
/// and Natural isotope policy.
///
/// The input is unchanged on success and failure.
///
/// # Errors
///
/// Preserves the error from molecular conversion or rendering.
pub fn export_smiles(input: &Molecule) -> Result<String, SmilesOutputError> {
    export_smiles_with(
        input,
        &SmilesIoConfig::opensmiles(),
        &ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        },
        &ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        },
    )
}

/// Export a molecule as SMILES text with explicit IO, chemistry, and resolve configuration.
///
/// Composes Smiles::convey and Smiles::render_with using the same IO configuration.
///
/// # Semantic properties
///
/// The input is unchanged on success and failure. Repeated output is text-identical for the
/// same input and configurations. Supported ingested molecules survive export and ingestion
/// under the same configurations with canonical GraphIR equality (para_stereo=false).
///
/// # Errors
///
/// Preserves the error from molecular conversion or rendering.
pub fn export_smiles_with(
    input: &Molecule,
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<String, SmilesOutputError> {
    Ok(Smiles::convey(input, model, resolve_config, io_config)?.render_with(io_config)?)
}

/// Export a reaction with the OpenSMILES configuration, SMILES valence preset,
/// and Natural isotope policy.
///
/// The input is unchanged on success and failure.
///
/// # Errors
///
/// Preserves materialization, side-specific conversion, and rendering errors.
pub fn export_reaction_smiles(input: &Reaction) -> Result<String, ReactionSmilesOutputError> {
    export_reaction_smiles_with(
        input,
        &SmilesIoConfig::opensmiles(),
        &ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        },
        &ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        },
    )
}

/// Export a reaction as SMILES text with explicit IO, chemistry, and resolve configuration.
///
/// Composes ReactionSmiles::convey and ReactionSmiles::render_with using the same IO configuration.
///
/// # Semantic properties
///
/// The input is unchanged on success and failure. Repeated output is text-identical for the
/// same input and configurations. Supported ingested reactions survive export and ingestion
/// under the same configurations with canonical GraphIR equality (para_stereo=false), including
/// atom correspondence. Original numeric map labels are not retained.
///
/// # Errors
///
/// Preserves materialization, side-specific conversion, and rendering errors.
pub fn export_reaction_smiles_with(
    input: &Reaction,
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<String, ReactionSmilesOutputError> {
    Ok(ReactionSmiles::convey(input, model, resolve_config, io_config)?.render_with(io_config)?)
}

fn convey_molecule(
    input: &Molecule,
    resolver: &Resolver<'_>,
    labels: impl Iterator<Item = (AtomId, u32)>,
) -> Result<TableMolecule, ConveyError> {
    let mut projected = input.clone();
    match resolver.project(&mut projected, ProjectFlags::all())? {
        Solution::Determined(()) => {}
        Solution::Underdetermined(()) => return Err(ConveyError::Underdetermined),
        Solution::Contradictory(error) => return Err(error.into()),
    }
    let molecule = &projected;
    let mut table = TableMolecule::empty();
    table.atoms.reserve(molecule.atoms().count());
    table.bonds.reserve(molecule.bonds().count());
    let mut labels = labels.peekable();
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
        if labels.peek().is_some_and(|&(id, _)| id == atom.id) {
            lowered.class = labels.next().map(|(_, class)| class);
        }
        elide_implicit_hydrogens(molecule, atom.id, resolver, &mut lowered);
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

fn elide_implicit_hydrogens(
    molecule: &Molecule,
    atom_id: AtomId,
    resolver: &Resolver<'_>,
    lowered: &mut Atom,
) {
    let Some(hydrogens) = lowered.implicit_hydrogens else {
        return;
    };
    let atom = molecule.atom(atom_id);
    if !matches!(atom.attributes.isotope_mass, IsotopeMassForm::Undetermined)
        || lowered.charge != Some(0)
        || lowered.unpaired_electrons != Some(0)
        || (lowered.aromatic == Some(true) && lowered.element != Some(Element::C))
        || lowered.class.is_some()
        || matches!(
            atom.constraints().tetrahedral_stereo(),
            Some(TetrahedralStereoForm::Stereo(_))
        )
        || !matches!(
            lowered.element,
            Some(
                Element::B
                    | Element::C
                    | Element::N
                    | Element::O
                    | Element::P
                    | Element::S
                    | Element::F
                    | Element::Cl
                    | Element::Br
                    | Element::I
            )
        )
    {
        return;
    }
    if resolver
        .valence
        .infer_implicit_hydrogens(molecule, atom_id, resolver.tie_break)
        == Some(i64::from(hydrogens))
    {
        lowered.implicit_hydrogens = None;
    }
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
    use std::collections::BTreeMap;
    use std::error::Error as _;

    use rstest::rstest;
    use umol_chem::element::Element;
    use umol_graph_core::AutomorphismAlgorithm;
    use umol_graph_ir::ir::{
        AtomDelta, AtomFieldChange, BondDelta, BondFieldChange, Canonicalize, CanonicalizeContext,
        Delta, Deltas,
    };
    use umol_graph_ir::mol_dsl_concrete;

    use super::*;
    use crate::ingest::{ingest_reaction_smiles, ingest_smiles, ingest_smiles_with};
    use crate::ops::model::{ChemistryModel, ValenceModel, ValenceTieBreak};
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
        let config = ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        };
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
            Smiles::convey(&source, &model, &config, &SmilesIoConfig::opensmiles())
                .unwrap()
                .into_table_ir(),
            expected
        );
    }

    #[rstest]
    #[case::first("[H]C", 0)]
    #[case::last("C[H]", 1)]
    fn test_smiles_convey_atom(
        #[case] input: &str,
        #[case] index: usize,
        #[values(false, true)] typing: bool,
    ) {
        let mut model = ChemistryModel {
            valence: if typing {
                ValenceModel::default()
            } else {
                ValenceModel::smiles()
            },
            ..Default::default()
        };
        model.valence.tie_break = ValenceTieBreak::MostSaturated;
        let config = ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        };
        let source = ingest_smiles(input).unwrap();
        let smiles =
            Smiles::convey(&source, &model, &config, &SmilesIoConfig::opensmiles()).unwrap();
        assert_eq!(
            smiles.as_table_ir().atoms[index],
            Atom {
                implicit_hydrogens: Some(0),
                charge: Some(0),
                lone_pairs: Some(0),
                unpaired_electrons: Some(0),
                multiplicity: Some(SpinMultiplicity::SINGLET),
                ..Atom::from_element(Element::H)
            }
        );
    }

    #[rstest]
    #[case::empty("", "")]
    #[case::chain("CCO", "CCO")]
    #[case::branched("CC(C)O", "CC(C)O")]
    #[case::aromatic("c1ccccc1", "c1ccccc1")]
    #[case::heteroaromatic("n1ccccc1", "[n]1ccccc1")]
    #[case::pyrrole("[nH]1cccc1", "[nH]1cccc1")]
    #[case::biphenyl("c1ccccc1-c2ccccc2", "c1ccccc1-c1ccccc1")]
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
        let original = ingest_smiles(input).unwrap();
        let source = original.clone();
        let smiles = Smiles::convey(&source, &model, &config, &io).unwrap();
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
    #[case::methane("C", "[CH4]", "C")]
    #[case::chain("CC", "[CH3][CH3]", "CC")]
    #[case::carbon_lone_pair("[CH2]", "[CH2]", "[CH2]")]
    #[case::substituted_lone_pair("[C](F)Cl", "[C](F)Cl", "[C](F)Cl")]
    #[case::nitrogen_lone_pairs("[NH]", "[NH]", "[NH]")]
    #[case::ammonia("N", "[NH3]", "N")]
    #[case::aromatic("c1ccccc1", "c1ccccc1", "c1ccccc1")]
    #[case::pyridine("n1ccccc1", "[n]1ccccc1", "[n]1ccccc1")]
    #[case::pyrrole("[nH]1cccc1", "[nH]1cccc1", "[nH]1cccc1")]
    #[case::radical("[CH3]", "[CH3]", "[CH3]")]
    #[case::charge("[NH4+]", "[NH4+]", "[NH4+]")]
    #[case::isotope("[13CH4]", "[13CH4]", "[13CH4]")]
    #[case::hydrogen("[H]C", "[H][CH3]", "[H]C")]
    #[case::bracket_element("[SiH4]", "[SiH4]", "[SiH4]")]
    #[case::tetrahedral_h("F[C@H](Cl)Br", "F[C@H](Cl)Br", "F[C@H](Cl)Br")]
    #[case::tetrahedral_atoms("F[C@](Cl)(Br)I", "F[C@](Cl)(Br)I", "F[C@](Cl)(Br)I")]
    #[case::tetrahedral_lone_pair("C[S@](=O)CC", "[CH3][S@](=O)[CH2][CH3]", "C[S@](=O)CC")]
    fn test_smiles_convey_hydrogens(
        #[case] input: &str,
        #[case] strict: &str,
        #[case] saturated: &str,
        #[values(false, true)] typing: bool,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)] policy: ValenceTieBreak,
    ) {
        let mut model = ChemistryModel {
            valence: if typing {
                ValenceModel::default()
            } else {
                ValenceModel::smiles()
            },
            ..Default::default()
        };
        model.valence.tie_break = policy;
        let config = ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        };
        let source = ingest_smiles(input).unwrap();
        let original = source.clone();
        let smiles =
            Smiles::convey(&source, &model, &config, &SmilesIoConfig::opensmiles()).unwrap();
        assert_eq!(
            smiles.render(),
            Ok(match policy {
                ValenceTieBreak::Strict => strict,
                ValenceTieBreak::MostSaturated => saturated,
            }
            .to_owned())
        );
        assert_eq!(source, original);
    }

    #[rstest]
    #[case::boron("[BH3]", "[BH3]", "B")]
    #[case::water("[OH2]", "[OH2]", "O")]
    fn test_smiles_convey_model(
        #[case] input: &str,
        #[case] counts: &str,
        #[case] atom_typing: &str,
        #[values(false, true)] typing: bool,
    ) {
        let mut model = ChemistryModel {
            valence: if typing {
                ValenceModel::default()
            } else {
                ValenceModel::smiles()
            },
            ..Default::default()
        };
        model.valence.tie_break = ValenceTieBreak::Strict;
        let config = ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        };
        let source = ingest_smiles(input).unwrap();
        let smiles =
            Smiles::convey(&source, &model, &config, &SmilesIoConfig::opensmiles()).unwrap();
        assert_eq!(
            smiles.render(),
            Ok(if typing { atom_typing } else { counts }.to_owned())
        );
    }

    #[rstest]
    #[case::natural(IsotopePolicy::Natural, None)]
    #[case::strict(IsotopePolicy::Strict, Some(4))]
    fn test_smiles_convey_isotope(
        #[case] isotope: IsotopePolicy,
        #[case] hydrogens: Option<u8>,
        #[values(false, true)] typing: bool,
    ) {
        let mut model = ChemistryModel {
            valence: if typing {
                ValenceModel::default()
            } else {
                ValenceModel::smiles()
            },
            ..Default::default()
        };
        model.valence.tie_break = ValenceTieBreak::MostSaturated;
        let source = ingest_smiles("C.[13CH4]").unwrap();
        let config = ResolveConfig {
            isotope,
            ..Default::default()
        };
        let table = Smiles::convey(&source, &model, &config, &SmilesIoConfig::opensmiles())
            .unwrap()
            .into_table_ir();
        let expected = TableMolecule {
            atoms: vec![
                Atom {
                    implicit_hydrogens: hydrogens,
                    charge: Some(0),
                    lone_pairs: Some(0),
                    unpaired_electrons: Some(0),
                    multiplicity: Some(SpinMultiplicity::SINGLET),
                    ..Atom::from_element(Element::C)
                },
                Atom {
                    isotope_mass: Some(13),
                    implicit_hydrogens: Some(4),
                    charge: Some(0),
                    lone_pairs: Some(0),
                    unpaired_electrons: Some(0),
                    multiplicity: Some(SpinMultiplicity::SINGLET),
                    ..Atom::from_element(Element::C)
                },
            ],
            ..TableMolecule::empty()
        };
        assert_eq!(table, expected);
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
            &model,
            &ResolveConfig::default(),
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
                &model,
                &ResolveConfig::default(),
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
                &model,
                &ResolveConfig::default(),
                &SmilesIoConfig::opensmiles()
            ),
            Err(expected)
        );
        assert_eq!(source, original);
    }

    #[rstest]
    #[case::empty(">>", ">>")]
    #[case::creation(">>C", ">>C")]
    #[case::deletion("C>>", "C>>")]
    #[case::unmapped("C>>O", "C>>O")]
    #[case::labels("[CH4:19]>>[CH4:19]", "[CH4:1]>>[CH4:1]")]
    #[case::unpaired_labels("[CH4:19]>>[OH2:7]", "C>>O")]
    #[case::substitution("[CH3:7][Cl:3]>>[CH3:7][OH:9]", "[CH3:1]Cl>>[CH3:1]O")]
    #[case::crossing("[CH3:7][OH:2]>>[OH:2][CH3:7]", "[CH3:1][OH:2]>>[CH3:1][OH:2]")]
    #[case::bond_order("[CH3:9][CH3:4]>>[CH2:9]=[CH2:4]", "[CH3:1][CH3:2]>>[CH2:1]=[CH2:2]")]
    #[case::compaction(
        "O.[CH3:19][CH3:2]>>[CH2:19]=[CH2:2].N",
        "O.[CH3:1][CH3:2]>>[CH2:1]=[CH2:2].N"
    )]
    #[case::disconnected("[Na+:18].[Cl-]>>[Na+:18].[Br-]", "[Na+:1].[Cl-]>>[Na+:1].[Br-]")]
    #[case::radical("[CH3:8]>>[CH3:8]", "[CH3:1]>>[CH3:1]")]
    #[case::isotope("[13CH4:8]>>[13CH4:8]", "[13CH4:1]>>[13CH4:1]")]
    #[case::tetrahedral(
        "[F:9][C@H:2]([Cl:7])[Br:3]>>[F:9][C@H:2]([Cl:7])[Br:3]",
        "[F:1][C@H:2]([Cl:3])[Br:4]>>[F:1][C@H:2]([Cl:3])[Br:4]"
    )]
    #[case::explicit_h(
        "[H:9][C@:5]([F:3])([Cl:7])[Br:8]>>[H:9][C@:5]([F:3])([Cl:7])[Br:8]",
        "[H:1][C@:2]([F:3])([Cl:4])[Br:5]>>[H:1][C@:2]([F:3])([Cl:4])[Br:5]"
    )]
    #[case::alkene(
        "[F:9]/[CH:2]=[CH:7]/[Cl:3]>>[F:9]/[CH:2]=[CH:7]\\[Cl:3]",
        "[F:1]/[CH:2]=[CH:3]/[Cl:4]>>[F:1]/[CH:2]=[CH:3]\\[Cl:4]"
    )]
    #[case::aromatic(
        "[cH:19]1[cH][cH][cH][cH][cH]1>>[cH:19]1[cH][cH][cH][cH][cH]1",
        "[cH:1]1ccccc1>>[cH:1]1ccccc1"
    )]
    fn test_reaction_smiles_convey_roundtrip(#[case] input: &str, #[case] expected: &str) {
        let source = ingest_reaction_smiles(input).unwrap();
        let original = source.clone();
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..Default::default()
        };
        let config = ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        };
        let boundary =
            ReactionSmiles::convey(&source, &model, &config, &SmilesIoConfig::opensmiles())
                .unwrap();
        assert_eq!(source, original);
        assert_eq!(boundary.as_table_ir().agents, TableMolecule::empty());
        assert_eq!(boundary.as_table_ir().comments, Vec::<String>::new());
        assert_eq!(
            boundary.as_table_ir().properties,
            TableReaction::empty().properties
        );
        let text = boundary.render().unwrap();
        assert_eq!(text, expected);
        let restored = ingest_reaction_smiles(&text).unwrap();
        assert!(source.canonical_eq(
            &restored,
            &CanonicalizeContext {
                para_stereo: false,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            }
        ));
        assert_eq!(
            ReactionSmiles::convey(&restored, &model, &config, &SmilesIoConfig::opensmiles())
                .unwrap(),
            boundary
        );
    }

    #[rstest]
    #[case::compaction(
        "O.[CH3:19][CH3:2]>>[CH2:2]=[CH2:19].N",
        vec![None, Some(1), Some(2)], vec![Some(1), Some(2), None],
        [(1, (vec![1], vec![0])), (2, (vec![2], vec![1]))].into()
    )]
    #[case::partial(
        "[CH3:7][Cl:3]>>[CH3:7][OH:9]",
        vec![Some(1), None], vec![Some(1), None],
        [(1, (vec![0], vec![0]))].into()
    )]
    fn test_reaction_smiles_convey_mapping(
        #[case] input: &str,
        #[case] reactants: Vec<Option<u32>>,
        #[case] products: Vec<Option<u32>>,
        #[case] mapping: BTreeMap<u32, (Vec<u32>, Vec<u32>)>,
    ) {
        let source = ingest_reaction_smiles(input).unwrap();
        let model = ChemistryModel::default();
        let table = ReactionSmiles::convey(
            &source,
            &model,
            &ResolveConfig::default(),
            &SmilesIoConfig::opensmiles(),
        )
        .unwrap()
        .into_table_ir();
        assert_eq!(
            table
                .reactants
                .atoms
                .iter()
                .map(|atom| atom.class)
                .collect::<Vec<_>>(),
            reactants
        );
        assert_eq!(
            table
                .products
                .atoms
                .iter()
                .map(|atom| atom.class)
                .collect::<Vec<_>>(),
            products
        );
        assert_eq!(table.atom_mapping, mapping);
    }

    #[rstest]
    #[case::materialization(
        Reaction::new(
            mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h3"] :bonds [[0 1 "1"]]}"#),
            Deltas::from_iter([Delta::Bond(BondDelta::ModifyField { id: BondId(0), change: BondFieldChange::Order { old: NumForm::Lit(2), new: NumForm::Lit(3) } })]),
        ),
        ReactionConveyError::Materialization(Contradiction)
    )]
    #[case::reactants(
        Reaction::new(mol_dsl_concrete!(r#"{:atoms ["C#h256"]}"#), Deltas::default()),
        ReactionConveyError::Reactants(ConveyError::Value { entity: Entity::Atom(AtomId(0)), field: "implicit hydrogens", value: "Lit(256)".to_owned() })
    )]
    #[case::products(
        Reaction::new(
            mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#),
            Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::ImplicitHydrogens { old: NumForm::Lit(4), new: NumForm::Lit(256) } })]),
        ),
        ReactionConveyError::Products(ConveyError::Value { entity: Entity::Atom(AtomId(0)), field: "implicit hydrogens", value: "Lit(256)".to_owned() })
    )]
    #[case::projection(
        Reaction::new(mol_dsl_concrete!(r#"{:atoms ["C#h3" "C#h3"] :bonds [[0 1 "1#c+"]]}"#), Deltas::default()),
        ReactionConveyError::Reactants(ConveyError::Projection(ProjectError::BondCharge { entity: Entity::Bond(BondId(0)), charge: NumForm::Lit(1) }))
    )]
    fn test_reaction_smiles_convey_error(
        #[case] source: Reaction,
        #[case] expected: ReactionConveyError,
    ) {
        let original = source.clone();
        let model = ChemistryModel::default();
        assert_eq!(
            ReactionSmiles::convey(
                &source,
                &model,
                &ResolveConfig::default(),
                &SmilesIoConfig::opensmiles()
            ),
            Err(expected)
        );
        assert_eq!(source, original);
    }

    #[rstest]
    #[case::ordinary("CCO", "CCO")]
    #[case::radical("[CH3]", "[CH3]")]
    #[case::isotope("[13CH4]", "[13CH4]")]
    #[case::aromatic("c1ccccc1", "c1ccccc1")]
    #[case::tetrahedral("F[C@H](Cl)Br", "F[C@H](Cl)Br")]
    #[case::explicit_h("[H][C@](F)(Cl)Br", "[H][C@](F)(Cl)Br")]
    #[case::coupled("F/C=C/C=C/F", "F/C=C/C=C/F")]
    fn test_export_smiles(#[case] input: &str, #[case] expected: &str) {
        let source = ingest_smiles(input).unwrap();
        let original = source.clone();
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..Default::default()
        };
        let resolve_config = ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        };
        let output = export_smiles(&source).unwrap();
        assert_eq!(output, expected);
        assert_eq!(
            export_smiles_with(
                &source,
                &SmilesIoConfig::opensmiles(),
                &model,
                &resolve_config
            ),
            Ok(output.clone())
        );
        assert_eq!(export_smiles(&source), Ok(output.clone()));
        let restored = ingest_smiles(&output).unwrap();
        let context = CanonicalizeContext {
            para_stereo: false,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        };
        assert!(source.canonical_eq(&restored, &context));
        assert_eq!(source, original);
    }

    #[rstest]
    #[case::conversion(
        mol_dsl_concrete!(r#"{:atoms ["C#h256"]}"#),
        SmilesOutputError::Convey(ConveyError::Value { entity: Entity::Atom(AtomId(0)), field: "implicit hydrogens", value: "Lit(256)".to_owned() })
    )]
    #[case::projection(
        mol_dsl_concrete!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1#c+"]]}"#),
        SmilesOutputError::Convey(ConveyError::Projection(ProjectError::BondCharge { entity: Entity::Bond(BondId(0)), charge: NumForm::Lit(1) }))
    )]
    #[case::render(
        mol_dsl_concrete!(r#"{:atoms ["C#h10"]}"#),
        SmilesOutputError::Render(SmilesRenderError::UnsupportedAtom { atom: 0, field: "implicit_hydrogens" })
    )]
    fn test_export_smiles_error(#[case] source: Molecule, #[case] expected: SmilesOutputError) {
        let original = source.clone();
        let model = ChemistryModel::default();
        let resolve_config = ResolveConfig::default();
        let error = export_smiles(&source).unwrap_err();
        assert_eq!(error, expected);
        assert_eq!(
            export_smiles_with(
                &source,
                &SmilesIoConfig::opensmiles(),
                &model,
                &resolve_config
            ),
            Err(expected.clone())
        );
        match &expected {
            SmilesOutputError::Convey(expected) => assert_eq!(
                error.source().unwrap().downcast_ref::<ConveyError>(),
                Some(expected)
            ),
            SmilesOutputError::Render(expected) => assert_eq!(
                error.source().unwrap().downcast_ref::<SmilesRenderError>(),
                Some(expected)
            ),
        }
        assert_eq!(source, original);
    }

    #[rstest]
    #[case::opensmiles(SmilesIoConfig::opensmiles(), Err(SmilesOutputError::Render(SmilesRenderError::UnsupportedBond { bond: 0, field: "donation" })))]
    #[case::lenient(SmilesIoConfig::lenient(), Ok("[NH3]->[BH3]".to_owned()))]
    fn test_export_smiles_with(
        #[case] io_config: SmilesIoConfig,
        #[case] expected: Result<String, SmilesOutputError>,
    ) {
        let source = mol_dsl_concrete!(
            r#"{:atoms ["N#h3#n" "B#h3"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#
        );
        let original = source.clone();
        let model = ChemistryModel::default();
        let resolve_config = ResolveConfig::default();
        let output = export_smiles_with(&source, &io_config, &model, &resolve_config);
        assert_eq!(output, expected);
        let explicit = Smiles::convey(&source, &model, &resolve_config, &io_config)
            .unwrap()
            .render_with(&io_config)
            .map_err(SmilesOutputError::from);
        assert_eq!(output, explicit);
        assert_eq!(source, original);
    }

    #[rstest]
    #[case::empty(">>", ">>")]
    #[case::creation(">>C", ">>C")]
    #[case::substitution("[CH3:7][Cl:3]>>[CH3:7][OH:9]", "[CH3:1]Cl>>[CH3:1]O")]
    #[case::aromatic("c1ccccc1>>c1ccccc1", "c1ccccc1>>c1ccccc1")]
    #[case::tetrahedral("F[C@H](Cl)Br>>F[C@@H](Cl)Br", "F[C@H](Cl)Br>>F[C@@H](Cl)Br")]
    #[case::coupled("F/C=C/C=C/F>>F/C=C/C=C\\F", "F/C=C/C=C/F>>F/C=C/C=C\\F")]
    fn test_export_reaction_smiles(#[case] input: &str, #[case] expected: &str) {
        let source = ingest_reaction_smiles(input).unwrap();
        let original = source.clone();
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..Default::default()
        };
        let resolve_config = ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        };
        let output = export_reaction_smiles(&source).unwrap();
        assert_eq!(output, expected);
        assert_eq!(
            export_reaction_smiles_with(
                &source,
                &SmilesIoConfig::opensmiles(),
                &model,
                &resolve_config
            ),
            Ok(output.clone())
        );
        assert_eq!(export_reaction_smiles(&source), Ok(output.clone()));
        let restored = ingest_reaction_smiles(&output).unwrap();
        let context = CanonicalizeContext {
            para_stereo: false,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        };
        assert!(source.canonical_eq(&restored, &context));
        assert_eq!(source, original);
    }

    #[rstest]
    #[case::materialization(
        Reaction::new(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#), Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::ImplicitHydrogens { old: NumForm::Lit(3), new: NumForm::Lit(2) } })])),
        ReactionSmilesOutputError::Convey(ReactionConveyError::Materialization(Contradiction))
    )]
    #[case::reactants_conversion(
        Reaction::new(mol_dsl_concrete!(r#"{:atoms ["C#h256"]}"#), Deltas::default()),
        ReactionSmilesOutputError::Convey(ReactionConveyError::Reactants(ConveyError::Value { entity: Entity::Atom(AtomId(0)), field: "implicit hydrogens", value: "Lit(256)".to_owned() }))
    )]
    #[case::products_conversion(
        Reaction::new(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#), Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::ImplicitHydrogens { old: NumForm::Lit(4), new: NumForm::Lit(256) } })])),
        ReactionSmilesOutputError::Convey(ReactionConveyError::Products(ConveyError::Value { entity: Entity::Atom(AtomId(0)), field: "implicit hydrogens", value: "Lit(256)".to_owned() }))
    )]
    #[case::reactants_render(
        Reaction::new(mol_dsl_concrete!(r#"{:atoms ["C#h10"]}"#), Deltas::default()),
        ReactionSmilesOutputError::Render(ReactionSmilesRenderError::Reactants(SmilesRenderError::UnsupportedAtom { atom: 0, field: "implicit_hydrogens" }))
    )]
    #[case::products_render(
        Reaction::new(mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#), Deltas::from_iter([Delta::Atom(AtomDelta::ModifyField { id: AtomId(0), change: AtomFieldChange::ImplicitHydrogens { old: NumForm::Lit(4), new: NumForm::Lit(10) } })])),
        ReactionSmilesOutputError::Render(ReactionSmilesRenderError::Products(SmilesRenderError::UnsupportedAtom { atom: 0, field: "implicit_hydrogens" }))
    )]
    fn test_export_reaction_smiles_error(
        #[case] source: Reaction,
        #[case] expected: ReactionSmilesOutputError,
    ) {
        let original = source.clone();
        let model = ChemistryModel::default();
        let resolve_config = ResolveConfig::default();
        let error = export_reaction_smiles(&source).unwrap_err();
        assert_eq!(error, expected);
        assert_eq!(
            export_reaction_smiles_with(
                &source,
                &SmilesIoConfig::opensmiles(),
                &model,
                &resolve_config
            ),
            Err(expected.clone())
        );
        match &expected {
            ReactionSmilesOutputError::Convey(expected) => assert_eq!(
                error
                    .source()
                    .unwrap()
                    .downcast_ref::<ReactionConveyError>(),
                Some(expected)
            ),
            ReactionSmilesOutputError::Render(expected) => assert_eq!(
                error
                    .source()
                    .unwrap()
                    .downcast_ref::<ReactionSmilesRenderError>(),
                Some(expected)
            ),
        }
        assert_eq!(
            error.source().unwrap().source().unwrap().to_string(),
            expected.source().unwrap().source().unwrap().to_string()
        );
        assert_eq!(source, original);
    }

    #[rstest]
    #[case::opensmiles(SmilesIoConfig::opensmiles(), Err(ReactionSmilesOutputError::Render(ReactionSmilesRenderError::Reactants(SmilesRenderError::UnsupportedBond { bond: 0, field: "donation" }))))]
    #[case::lenient(SmilesIoConfig::lenient(), Ok("[NH3:1]->[BH3:2]>>[NH3:1]->[BH3:2]".to_owned()))]
    fn test_export_reaction_smiles_with(
        #[case] io_config: SmilesIoConfig,
        #[case] expected: Result<String, ReactionSmilesOutputError>,
    ) {
        let source = Reaction::new(
            mol_dsl_concrete!(
                r#"{:atoms ["N#h3#n" "B#h3"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#
            ),
            Deltas::default(),
        );
        let original = source.clone();
        let model = ChemistryModel::default();
        let resolve_config = ResolveConfig::default();
        let output = export_reaction_smiles_with(&source, &io_config, &model, &resolve_config);
        assert_eq!(output, expected);
        let explicit = ReactionSmiles::convey(&source, &model, &resolve_config, &io_config)
            .unwrap()
            .render_with(&io_config)
            .map_err(ReactionSmilesOutputError::from);
        assert_eq!(output, explicit);
        assert_eq!(source, original);
    }
}
