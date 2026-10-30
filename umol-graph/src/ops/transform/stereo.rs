//! Stereo transformations and their configuration.

mod enumeration;

use std::collections::{BTreeMap, BTreeSet};
use std::iter;

use thiserror::Error;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AtomConstraintKey, AtomHandle, BondConstraintKey, BondHandle, ConstitutionColoring, Edits,
    GraphSymmetryConfig, Molecule, Normalize, RingConfig, StereoAtomHandle, StereoBondHandle,
    StereoConfigurationForm, StereoCoset, StereoKind, Stereogenicity,
};

use self::enumeration::{AbsoluteStereoAssignments, StereoAssignmentChecks};
use crate::ops::model::StereoModel;
use crate::ops::stereo::{
    StereoInconsistency, StereoPerception, StereoPerceptionError, StereoSelection, StereoSite,
};
use crate::ops::transform::Transformer;
use crate::ops::validate::constraint::ConstraintInvariantsError;

/// Site selection and algorithms for stereo addition.
///
/// Default selects tetrahedral and cis/trans kinds, stereogenic sites only,
/// and all topology sites. It uses Nauty and RingConfig::default().
///
/// Construction preserves the supplied selection without checking ids against
/// a molecule. The executing operation checks ids before indexing or mutation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AddStereoConfig {
    pub selection: StereoSelection,
    pub automorphism_algorithm: AutomorphismAlgorithm,
    pub ring_config: RingConfig,
}

impl Default for AddStereoConfig {
    fn default() -> Self {
        Self {
            selection: StereoSelection {
                sites: None,
                kinds: Some(BTreeSet::from([
                    StereoKind::Tetrahedral,
                    StereoKind::CisTrans,
                ])),
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            },
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            ring_config: RingConfig::default(),
        }
    }
}

/// Stereo addition could not interpret its input or realize selected evidence.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum AddStereoError {
    #[error(transparent)]
    Perception(#[from] StereoPerceptionError),
    #[error(transparent)]
    Inconsistency(#[from] StereoInconsistency),
}

/// Add selected missing tetrahedral atoms and cis/trans bonds.
///
/// Uses StereoPerception's recognition and classification with the supplied
/// model and configuration. New cosets remain undetermined unless source
/// evidence supplies them. Existing entities, their frames and configurations,
/// source assertions, and molecular topology are preserved.
///
/// All rejection occurs during planning. Undetermined fields outside the
/// selected processing scope are retained. Recognition failures and incompatible
/// selected evidence return AddStereoError. The iterator lazily produces one
/// independent result on success and no result on rejection.
///
/// # Semantic properties
///
/// Addition is idempotent for the same model and configuration. Consuming and
/// in-place results agree under structural equality. An in-place rejection
/// leaves the molecule unchanged, and every result preserves its integrity.
#[derive(Clone, Debug)]
pub struct StereoAdder {
    perception: StereoPerception,
    config: AddStereoConfig,
}

impl StereoAdder {
    /// Construct an adder with the default selection and algorithms.
    pub fn new(model: &StereoModel) -> Self {
        Self::with_config(model, AddStereoConfig::default())
    }

    /// Store a model snapshot and explicit configuration for stereo addition.
    ///
    /// Site ids are checked against each input molecule during transformation.
    pub fn with_config(model: &StereoModel, config: AddStereoConfig) -> Self {
        Self {
            perception: StereoPerception::new(model),
            config,
        }
    }

    fn plan_transform(&self, molecule: &Molecule) -> Result<Edits, AddStereoError> {
        let derivation = self.perception.recognize(
            molecule,
            &self.config.selection,
            self.config.automorphism_algorithm,
            self.config.ring_config,
        )?;
        if let Some(&inconsistency) = derivation.inconsistencies.first() {
            return Err(inconsistency.into());
        }

        let mut edits = Edits::new();
        for (site, ligands, form) in derivation.stereo_atoms {
            edits.add_stereo_atom(
                AtomHandle::Id(site),
                ligands
                    .into_iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                form,
            );
        }
        for (site, ligands, form) in derivation.stereo_bonds {
            edits.add_stereo_bond(
                BondHandle::Id(site),
                ligands
                    .into_iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                form,
            );
        }
        Ok(edits)
    }
}

impl Transformer for StereoAdder {
    type Error = AddStereoError;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, AddStereoError> {
        let edits = self.plan_transform(&molecule)?;
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("stereo addition plan applies to its input");
        Ok(editor
            .finish()
            .expect("stereo addition plan preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), AddStereoError> {
        let edits = self.plan_transform(molecule)?;
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        iter::once_with(move || self.transform(molecule.clone()).ok()).flatten()
    }
}

/// Which configuration states stereo removal accepts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RemoveStereoMode {
    /// Accept every configuration, including an undetermined kind.
    All,
    /// Accept only a known kind with an explicitly undetermined coset.
    Undetermined,
}

impl RemoveStereoMode {
    /// Whether the configuration satisfies this removal mode.
    ///
    /// Undetermined tests the coset state, independently of the known kind.
    /// An unknown kind has no coset to test; finite sets and symbolic cosets
    /// are not undetermined cosets. No classification or resolution is performed.
    pub fn matches(self, configuration: &StereoConfigurationForm) -> bool {
        match self {
            Self::All => true,
            Self::Undetermined => matches!(configuration.coset(), Some(StereoCoset::Undetermined)),
        }
    }
}

/// Site selection, removal mode, and algorithm for stereo removal.
///
/// Default accepts all configurations and leaves sites, kinds, and
/// stereogenicities unrestricted. It uses Nauty. One mode applies uniformly
/// to every selected site and kind.
///
/// Construction preserves the supplied selection without checking ids against
/// a molecule or calculating symmetry. The executing operation checks ids
/// before indexing or mutation and classifies sites only when selection needs it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RemoveStereoConfig {
    pub selection: StereoSelection,
    pub mode: RemoveStereoMode,
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

impl Default for RemoveStereoConfig {
    fn default() -> Self {
        Self {
            selection: StereoSelection::default(),
            mode: RemoveStereoMode::All,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        }
    }
}

/// Stereo removal could not interpret its selection or classify selected entities.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum RemoveStereoError {
    #[error(transparent)]
    Perception(#[from] StereoPerceptionError),
}

/// Remove selected stereo entities and their corresponding #T/#C assertions.
///
/// The mode and site/kind selectors apply before stereogenicity is calculated.
/// Classification uses StereoPerception with the unchanged molecular context
/// and is omitted when the class selector is unrestricted. Such removal accepts
/// kindless, symbolic, and other-geometry records without interpreting them.
/// Undetermined tests only an explicitly undetermined coset of a known kind.
///
/// Unrelated assertions, topology, and retained entity data are preserved;
/// graph IR compacts entity ids and their references after removal. Invalid
/// selected topology ids or unsupported required classification produce
/// RemoveStereoError before any mutation. The iterator lazily yields one
/// independent result on success and none on rejection.
///
/// # Semantic properties
///
/// Consuming and in-place results agree under structural equality. In-place
/// rejection preserves the input under structural equality. Successful results
/// preserve molecule integrity. Without a class filter, removal is idempotent
/// for the same config.
#[derive(Clone, Debug)]
pub struct StereoRemover {
    perception: StereoPerception,
    config: RemoveStereoConfig,
}

impl StereoRemover {
    /// Construct a remover with unrestricted All selection and default algorithms.
    pub fn new(model: &StereoModel) -> Self {
        Self::with_config(model, RemoveStereoConfig::default())
    }

    /// Store a model snapshot and explicit configuration for stereo removal.
    ///
    /// Site ids are checked against each input molecule during transformation.
    pub fn with_config(model: &StereoModel, config: RemoveStereoConfig) -> Self {
        Self {
            perception: StereoPerception::new(model),
            config,
        }
    }

    fn plan_transform(&self, molecule: &Molecule) -> Result<Edits, RemoveStereoError> {
        let selection = &self.config.selection;
        if let Some(sites) = &selection.sites {
            for &site in sites {
                let available = match site {
                    StereoSite::Atom(atom) => molecule.atoms().get(atom).is_some(),
                    StereoSite::Bond(bond) => molecule.bonds().get(bond).is_some(),
                };
                if !available {
                    return Err(StereoPerceptionError::InvalidSite { site }.into());
                }
            }
        }
        if selection.sites.as_ref().is_some_and(BTreeSet::is_empty)
            || selection.kinds.as_ref().is_some_and(BTreeSet::is_empty)
            || selection
                .stereogenicities
                .as_ref()
                .is_some_and(BTreeSet::is_empty)
        {
            return Ok(Edits::new());
        }

        let mut sites: BTreeSet<_> = molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| {
                (
                    StereoSite::Atom(stereo.site_id()),
                    &stereo.attributes().configuration,
                )
            })
            .chain(molecule.stereo_bonds().iter().map(|stereo| {
                (
                    StereoSite::Bond(stereo.site_id()),
                    &stereo.attributes().configuration,
                )
            }))
            .filter(|(site, configuration)| {
                self.config.mode.matches(configuration)
                    && selection
                        .sites
                        .as_ref()
                        .is_none_or(|sites| sites.contains(site))
                    && selection.kinds.as_ref().is_none_or(|kinds| {
                        configuration
                            .kind()
                            .is_some_and(|kind| kinds.contains(&kind))
                    })
            })
            .map(|(site, _)| site)
            .collect();
        if let Some(classes) = &selection.stereogenicities {
            let classifications =
                self.perception
                    .classify(molecule, &sites, self.config.automorphism_algorithm)?;
            sites.retain(|site| classes.contains(&classifications[site]));
        }

        let mut edits = Edits::new();
        let mut atoms = Vec::new();
        let mut atom_constraints = Vec::new();
        for stereo in molecule.stereo_atoms().iter() {
            let atom = stereo.site_id();
            if !sites.contains(&StereoSite::Atom(atom)) {
                continue;
            }
            atoms.push((
                StereoAtomHandle::Id(stereo.id()),
                AtomHandle::Id(atom),
                stereo
                    .ligand_ids()
                    .iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                stereo.attributes().clone(),
            ));
            if let Some(asserted) = molecule
                .atom(atom)
                .attributes()
                .constraints
                .get(AtomConstraintKey::TetrahedralStereo)
            {
                atom_constraints.push((AtomHandle::Id(atom), Some(asserted.clone()), None));
            }
        }
        let mut bonds = Vec::new();
        let mut bond_constraints = Vec::new();
        for stereo in molecule.stereo_bonds().iter() {
            let bond = stereo.site_id();
            if !sites.contains(&StereoSite::Bond(bond)) {
                continue;
            }
            bonds.push((
                StereoBondHandle::Id(stereo.id()),
                BondHandle::Id(bond),
                stereo
                    .ligand_ids()
                    .iter()
                    .map(|ligand| (AtomHandle::Id(ligand.atom_id), ligand.kind))
                    .collect(),
                stereo.attributes().clone(),
            ));
            if let Some(asserted) = molecule
                .bond(bond)
                .attributes()
                .constraints
                .get(BondConstraintKey::CisTransStereo)
            {
                bond_constraints.push((BondHandle::Id(bond), Some(asserted.clone()), None));
            }
        }
        if !atoms.is_empty() {
            edits.remove_stereo_atoms(atoms);
        }
        if !bonds.is_empty() {
            edits.remove_stereo_bonds(bonds);
        }
        if !atom_constraints.is_empty() {
            edits.modify_atom_constraints(AtomConstraintKey::TetrahedralStereo, atom_constraints);
        }
        if !bond_constraints.is_empty() {
            edits.modify_bond_constraints(BondConstraintKey::CisTransStereo, bond_constraints);
        }
        Ok(edits)
    }
}

impl Transformer for StereoRemover {
    type Error = RemoveStereoError;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, RemoveStereoError> {
        let edits = self.plan_transform(&molecule)?;
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("stereo removal plan applies to its input");
        Ok(editor
            .finish()
            .expect("stereo removal plan preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), RemoveStereoError> {
        let edits = self.plan_transform(molecule)?;
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        iter::once_with(move || self.transform(molecule.clone()).ok()).flatten()
    }
}

/// Site selection and algorithm for stereo enumeration.
///
/// Default selects tetrahedral and cis/trans kinds, all topology sites, and
/// all stereogenicity values. It uses Nauty. Processing order is supplied
/// separately from this configuration.
///
/// Construction preserves the supplied selection without checking ids against
/// a molecule. The executing operation checks ids before indexing or mutation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EnumerateStereoConfig {
    pub selection: StereoSelection,
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

impl Default for EnumerateStereoConfig {
    fn default() -> Self {
        Self {
            selection: StereoSelection {
                sites: None,
                kinds: Some(BTreeSet::from([
                    StereoKind::Tetrahedral,
                    StereoKind::CisTrans,
                ])),
                stereogenicities: None,
            },
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        }
    }
}

/// Stereo enumeration could not prepare or satisfy the selected input.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum EnumerateStereoError {
    #[error(transparent)]
    Perception(#[from] StereoPerceptionError),
    #[error(transparent)]
    Constraint(#[from] ConstraintInvariantsError),
    #[error("stereo site order has missing sites {missing:?} and duplicate sites {duplicates:?}")]
    InvalidSiteOrder {
        missing: Vec<StereoSite>,
        duplicates: Vec<StereoSite>,
    },
    #[error("no consistent absolute stereo assignment")]
    NoAssignment,
}

/// Determine selected absolute tetrahedral and cis/trans configurations.
///
/// Visits existing selected entities in the supplied site order, trying admitted
/// cosets in ascending order in their stored frames. Fixed cosets are preserved;
/// undetermined cosets and finite literal sets supply their finite domains.
/// Constant expressions are normalized under their kind. Variable expressions
/// and other geometries are unsupported when selected for processing.
///
/// Site, kind, and stereogenicity selectors intersect. Classification uses
/// StereoPerception in the whole molecular context when requested. Automatic
/// recognition rules do not reject existing arrangements. Selected sites must
/// occur exactly once in the order; entries outside that set are ignored.
/// Source assertions, topology, retained configurations, and ligand frames are
/// preserved. Consistency checks allow unrelated undetermined fields.
///
/// # Semantic properties
///
/// Consuming and in-place methods select the same first consistent assignment
/// under structural equality. In-place rejection leaves the input unchanged.
/// The iterator starts with that same result and emits every consistent
/// assignment without symmetry deduplication. Empty work yields one unchanged
/// result when its checks are consistent; no solution rejects singular calls
/// and yields no iterator results. Every result preserves molecule integrity.
#[derive(Clone, Debug)]
pub struct StereoEnumerator {
    model: StereoModel,
    config: EnumerateStereoConfig,
    site_order: Vec<StereoSite>,
}

impl StereoEnumerator {
    /// Store a model snapshot and processing order with the default config.
    ///
    /// Contextual site and order checks occur during execution.
    pub fn new(model: &StereoModel, site_order: Vec<StereoSite>) -> Self {
        Self::with_config(model, EnumerateStereoConfig::default(), site_order)
    }

    /// Store a model snapshot, owned configuration, and processing order.
    ///
    /// Construction does not check the order or bind ids to a molecule.
    pub fn with_config(
        model: &StereoModel,
        config: EnumerateStereoConfig,
        site_order: Vec<StereoSite>,
    ) -> Self {
        Self {
            model: model.clone(),
            config,
            site_order,
        }
    }

    /// Issue an owned iterator over the input's consistent absolute assignments.
    ///
    /// Captures the source snapshot and prepares selection and domains now;
    /// constructs candidate molecules lazily on next. The issued run can outlive
    /// both engine and source and is unaffected by their subsequent mutation.
    /// Preparation rejection or no consistent assignment produces an empty
    /// iterator; singular Transformer methods report the rejection cause.
    pub fn transform_iter(&self, molecule: &Molecule) -> EnumerateStereoIter {
        EnumerateStereoIter {
            source: molecule.clone(),
            enumeration: self.prepare(molecule).ok(),
        }
    }

    fn prepare(&self, molecule: &Molecule) -> Result<StereoEnumeration, EnumerateStereoError> {
        let selection = &self.config.selection;
        if let Some(sites) = &selection.sites {
            for &site in sites {
                let available = match site {
                    StereoSite::Atom(atom) => molecule.atoms().get(atom).is_some(),
                    StereoSite::Bond(bond) => molecule.bonds().get(bond).is_some(),
                };
                if !available {
                    return Err(StereoPerceptionError::InvalidSite { site }.into());
                }
            }
        }
        let mut configurations: BTreeMap<_, _> = molecule
            .stereo_atoms()
            .iter()
            .map(|stereo| {
                (
                    StereoSite::Atom(stereo.site_id()),
                    stereo.attributes().configuration.clone(),
                )
            })
            .chain(molecule.stereo_bonds().iter().map(|stereo| {
                (
                    StereoSite::Bond(stereo.site_id()),
                    stereo.attributes().configuration.clone(),
                )
            }))
            .filter(|(site, configuration)| {
                selection
                    .sites
                    .as_ref()
                    .is_none_or(|sites| sites.contains(site))
                    && selection.kinds.as_ref().is_none_or(|kinds| {
                        configuration
                            .kind()
                            .is_some_and(|kind| kinds.contains(&kind))
                    })
            })
            .collect();
        if let Some(classes) = &selection.stereogenicities {
            if classes.is_empty() {
                configurations.clear();
            } else {
                let mut normalized = molecule.clone();
                for (&site, configuration) in &configurations {
                    if matches!(
                        (site, configuration.kind()),
                        (StereoSite::Atom(_), Some(StereoKind::CisTrans))
                            | (StereoSite::Bond(_), Some(StereoKind::Tetrahedral))
                    ) {
                        return Err(StereoPerceptionError::UnsupportedConfiguration { site }.into());
                    }
                    if let Ok(configuration) = configuration.clone().normalize() {
                        match site {
                            StereoSite::Atom(atom) => {
                                let id = normalized
                                    .stereo_atoms()
                                    .at(atom)
                                    .expect("selected stereo atom")
                                    .id();
                                normalized
                                    .stereo_atom_mut(id)
                                    .attributes_mut()
                                    .configuration = configuration;
                            }
                            StereoSite::Bond(bond) => {
                                let id = normalized
                                    .stereo_bonds()
                                    .at(bond)
                                    .expect("selected stereo bond")
                                    .id();
                                normalized
                                    .stereo_bond_mut(id)
                                    .attributes_mut()
                                    .configuration = configuration;
                            }
                        }
                    }
                }
                let sites = configurations.keys().copied().collect();
                let classifications = StereoPerception::new(&self.model).classify(
                    &normalized,
                    &sites,
                    self.config.automorphism_algorithm,
                )?;
                configurations.retain(|site, _| classes.contains(&classifications[site]));
            }
        }
        let sites: BTreeSet<_> = configurations.keys().copied().collect();
        let mut seen = BTreeSet::new();
        let mut duplicates = Vec::new();
        let mut ordered = Vec::new();
        for &site in self.site_order.iter().filter(|site| sites.contains(site)) {
            if seen.insert(site) {
                ordered.push(site);
            } else if !duplicates.contains(&site) {
                duplicates.push(site);
            }
        }
        let missing: Vec<_> = sites.difference(&seen).copied().collect();
        if !missing.is_empty() || !duplicates.is_empty() {
            return Err(EnumerateStereoError::InvalidSiteOrder {
                missing,
                duplicates,
            });
        }
        for configuration in configurations.values() {
            if let Some(kind) = configuration.kind() {
                if self
                    .model
                    .kind_model(kind)
                    .is_some_and(|model| model.fluxionality)
                {
                    return Err(StereoPerceptionError::UnsupportedFluxionality { kind }.into());
                }
            }
        }
        let assignments = AbsoluteStereoAssignments::new(ordered.into_iter().map(|site| {
            (
                site,
                configurations
                    .remove(&site)
                    .expect("ordered selected configuration"),
            )
        }))?;
        let checks = StereoAssignmentChecks::new(molecule, &sites)?;
        Ok(StereoEnumeration {
            assignments,
            checks,
            symmetry_config: GraphSymmetryConfig {
                coloring: ConstitutionColoring::full(),
                iterate_to_fixpoint: self.model.para_stereo,
                max_iterations: 16,
                automorphism_algorithm: self.config.automorphism_algorithm,
            },
        })
    }

    fn plan_transform(&self, molecule: &Molecule) -> Result<Edits, EnumerateStereoError> {
        let (_, edits) = self
            .prepare(molecule)?
            .next(molecule)?
            .ok_or(EnumerateStereoError::NoAssignment)?;
        Ok(edits)
    }
}

impl Transformer for StereoEnumerator {
    type Error = EnumerateStereoError;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, Self::Error> {
        let edits = self.plan_transform(&molecule)?;
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("stereo enumeration plan applies to its input");
        Ok(editor
            .finish()
            .expect("stereo enumeration plan preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error> {
        let edits = self.plan_transform(molecule)?;
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        StereoEnumerator::transform_iter(self, molecule)
    }
}

/// Owned consistent assignments issued by StereoEnumerator::transform_iter.
///
/// Has no independent public constructor. Owns a source snapshot, prepared
/// domains, order, and cursor. Candidates are constructed lazily and yielded
/// results are independently mutable. Rejection and exhaustion terminate the
/// iterator permanently. Symmetry-equivalent assignments remain distinct items.
pub struct EnumerateStereoIter {
    source: Molecule,
    enumeration: Option<StereoEnumeration>,
}

impl Iterator for EnumerateStereoIter {
    type Item = Molecule;

    fn next(&mut self) -> Option<Molecule> {
        match self.enumeration.as_mut()?.next(&self.source) {
            Ok(Some((molecule, _))) => Some(molecule),
            Ok(None) | Err(_) => {
                self.enumeration = None;
                None
            }
        }
    }
}

impl iter::FusedIterator for EnumerateStereoIter {}

struct StereoEnumeration {
    assignments: AbsoluteStereoAssignments,
    checks: StereoAssignmentChecks,
    symmetry_config: GraphSymmetryConfig<ConstitutionColoring>,
}

impl StereoEnumeration {
    fn next(
        &mut self,
        source: &Molecule,
    ) -> Result<Option<(Molecule, Edits)>, EnumerateStereoError> {
        for assignment in self.assignments.by_ref() {
            let mut candidate = source.clone();
            for (site, configuration) in &assignment {
                match site {
                    StereoSite::Atom(atom) => {
                        let id = candidate
                            .stereo_atoms()
                            .at(*atom)
                            .expect("prepared stereo atom")
                            .id();
                        if candidate.stereo_atom(id).attributes().configuration != *configuration {
                            candidate.stereo_atom_mut(id).attributes_mut().configuration =
                                configuration.clone();
                        }
                    }
                    StereoSite::Bond(bond) => {
                        let id = candidate
                            .stereo_bonds()
                            .at(*bond)
                            .expect("prepared stereo bond")
                            .id();
                        if candidate.stereo_bond(id).attributes().configuration != *configuration {
                            candidate.stereo_bond_mut(id).attributes_mut().configuration =
                                configuration.clone();
                        }
                    }
                }
            }
            if !self
                .checks
                .is_consistent(&candidate, &self.symmetry_config)?
            {
                continue;
            }
            let mut atoms = Vec::new();
            let mut bonds = Vec::new();
            for (site, configuration) in assignment {
                match site {
                    StereoSite::Atom(atom) => {
                        let stereo = source
                            .stereo_atoms()
                            .at(atom)
                            .expect("prepared stereo atom");
                        let old = &stereo.attributes().configuration;
                        if *old != configuration {
                            atoms.push((
                                StereoAtomHandle::Id(stereo.id()),
                                old.clone(),
                                configuration,
                            ));
                        }
                    }
                    StereoSite::Bond(bond) => {
                        let stereo = source
                            .stereo_bonds()
                            .at(bond)
                            .expect("prepared stereo bond");
                        let old = &stereo.attributes().configuration;
                        if *old != configuration {
                            bonds.push((
                                StereoBondHandle::Id(stereo.id()),
                                old.clone(),
                                configuration,
                            ));
                        }
                    }
                }
            }
            let mut edits = Edits::new();
            if !atoms.is_empty() {
                edits.modify_stereo_atom_configurations(atoms);
            }
            if !bonds.is_empty() {
                edits.modify_stereo_bond_configurations(bonds);
            }
            return Ok(Some((candidate, edits)));
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    mod adder;
    mod enumerator;
    mod remover;

    use rstest::rstest;
    use umol_graph_ir::ir::{AtomId, BondId, StereoTerm};

    use super::*;
    use crate::ops::stereo::StereoSite;

    #[rstest]
    fn test_add_stereo_config_default() {
        assert_eq!(
            AddStereoConfig::default(),
            AddStereoConfig {
                selection: StereoSelection {
                    sites: None,
                    kinds: Some(BTreeSet::from([
                        StereoKind::Tetrahedral,
                        StereoKind::CisTrans,
                    ])),
                    stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
                },
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                ring_config: RingConfig::default(),
            },
        );
    }

    #[rstest]
    #[case::default(
        AddStereoConfig::default(),
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0))],
    )]
    #[case::atom_site(
        AddStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0))])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0))],
    )]
    #[case::bond_site(
        AddStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Bond(BondId(0))])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::tetrahedral(
        AddStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0))],
    )]
    #[case::cis_trans(
        AddStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::intersection(
        AddStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([
                    StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)),
                    StereoSite::Bond(BondId(0)),
                ])),
                kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::prochiral(
        AddStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Prochiral])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(1))],
    )]
    #[case::symmetric(
        AddStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(2)), StereoSite::Bond(BondId(1))],
    )]
    #[case::all_stereogenicities(
        AddStereoConfig {
            selection: StereoSelection {
                stereogenicities: None,
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)),
            StereoSite::Atom(AtomId(2)), StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
    )]
    #[case::empty_kinds(
        AddStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::new()),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![],
    )]
    #[case::unbound_sites(
        AddStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([
                    StereoSite::Atom(AtomId(u32::MAX)), StereoSite::Bond(BondId(u32::MAX)),
                ])),
                ..AddStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![],
    )]
    fn test_add_stereo_config_selection(
        #[case] config: AddStereoConfig,
        #[case] expected: Vec<StereoSite>,
    ) {
        let candidates = [
            (
                StereoSite::Atom(AtomId(0)),
                StereoKind::Tetrahedral,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Atom(AtomId(1)),
                StereoKind::Tetrahedral,
                Stereogenicity::Prochiral,
            ),
            (
                StereoSite::Atom(AtomId(2)),
                StereoKind::Tetrahedral,
                Stereogenicity::Symmetric,
            ),
            (
                StereoSite::Atom(AtomId(3)),
                StereoKind::SquarePlanar,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Bond(BondId(0)),
                StereoKind::CisTrans,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Bond(BondId(1)),
                StereoKind::CisTrans,
                Stereogenicity::Symmetric,
            ),
        ];
        let selected: Vec<_> = candidates
            .into_iter()
            .filter(|&(site, kind, stereogenicity)| {
                config
                    .selection
                    .matches(site, Some(kind), Some(stereogenicity))
            })
            .map(|(site, _, _)| site)
            .collect();
        assert_eq!(selected, expected);
    }

    #[rstest]
    #[case::undetermined(StereoCoset::Undetermined, true)]
    #[case::literal_zero(StereoCoset::Lit(0), false)]
    #[case::literal_one(StereoCoset::Lit(1), false)]
    #[case::finite_set(StereoCoset::lit_set([0, 1]), false)]
    #[case::singleton_set(StereoCoset::lit_set([1]), false)]
    #[case::empty_set(StereoCoset::lit_set([]), false)]
    #[case::symbolic(StereoCoset::term(StereoTerm::var("configuration")), false)]
    fn test_remove_stereo_mode_matches(
        #[case] coset: StereoCoset,
        #[case] undetermined: bool,
        #[values(
            StereoKind::Tetrahedral,
            StereoKind::CisTrans,
            StereoKind::SquarePlanar
        )]
        kind: StereoKind,
        #[values(RemoveStereoMode::All, RemoveStereoMode::Undetermined)] mode: RemoveStereoMode,
    ) {
        let configuration = StereoConfigurationForm::Kinded(kind, coset);
        let expected = match mode {
            RemoveStereoMode::All => true,
            RemoveStereoMode::Undetermined => undetermined,
        };
        assert_eq!(mode.matches(&configuration), expected);
    }

    #[rstest]
    #[case::all(RemoveStereoMode::All, true)]
    #[case::undetermined(RemoveStereoMode::Undetermined, false)]
    fn test_remove_stereo_mode_matches_kindless(
        #[case] mode: RemoveStereoMode,
        #[case] expected: bool,
    ) {
        assert_eq!(
            mode.matches(&StereoConfigurationForm::Undetermined),
            expected
        );
    }

    #[rstest]
    fn test_remove_stereo_config_default() {
        assert_eq!(
            RemoveStereoConfig::default(),
            RemoveStereoConfig {
                selection: StereoSelection {
                    sites: None,
                    kinds: None,
                    stereogenicities: None,
                },
                mode: RemoveStereoMode::All,
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            },
        );
    }

    #[rstest]
    #[case::unrestricted(
        StereoSelection::default(),
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)),
            StereoSite::Atom(AtomId(2)), StereoSite::Atom(AtomId(3)),
            StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1)),
            StereoSite::Atom(AtomId(4)), StereoSite::Bond(BondId(2))],
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(3)), StereoSite::Bond(BondId(0))],
    )]
    #[case::tetrahedral(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])), ..Default::default() },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2))],
        vec![StereoSite::Atom(AtomId(0))],
    )]
    #[case::cis_trans(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::CisTrans])), ..Default::default() },
        vec![StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::both_kinds(
        StereoSelection { kinds: Some(BTreeSet::from([StereoKind::Tetrahedral, StereoKind::CisTrans])),
            ..Default::default() },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2)),
            StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0))],
    )]
    #[case::site_kind_intersection(
        StereoSelection {
            sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0)),
                StereoSite::Bond(BondId(1))])),
            kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
            stereogenicities: None,
        },
        vec![StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::stereogenic(
        StereoSelection {
            kinds: Some(BTreeSet::from([StereoKind::Tetrahedral, StereoKind::CisTrans])),
            stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0))],
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Bond(BondId(0))],
    )]
    #[case::non_stereogenic(
        StereoSelection {
            kinds: Some(BTreeSet::from([StereoKind::Tetrahedral, StereoKind::CisTrans])),
            stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric, Stereogenicity::Prochiral])),
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2)), StereoSite::Bond(BondId(1))],
        vec![],
    )]
    #[case::empty_kinds(
        StereoSelection { kinds: Some(BTreeSet::new()), ..Default::default() },
        vec![], vec![],
    )]
    #[case::empty_sites(
        StereoSelection { sites: Some(BTreeSet::new()), ..Default::default() },
        vec![], vec![],
    )]
    fn test_remove_stereo_config_selection(
        #[case] selection: StereoSelection,
        #[case] all: Vec<StereoSite>,
        #[case] undetermined: Vec<StereoSite>,
        #[values(RemoveStereoMode::All, RemoveStereoMode::Undetermined)] mode: RemoveStereoMode,
    ) {
        let config = RemoveStereoConfig {
            selection,
            mode,
            ..Default::default()
        };
        let candidates = [
            (
                StereoSite::Atom(AtomId(0)),
                StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined),
                Some(Stereogenicity::Stereogenic),
            ),
            (
                StereoSite::Atom(AtomId(1)),
                StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                Some(Stereogenicity::Prochiral),
            ),
            (
                StereoSite::Atom(AtomId(2)),
                StereoConfigurationForm::Kinded(
                    StereoKind::Tetrahedral,
                    StereoCoset::lit_set([0, 1]),
                ),
                Some(Stereogenicity::Symmetric),
            ),
            (
                StereoSite::Atom(AtomId(3)),
                StereoConfigurationForm::Kinded(
                    StereoKind::SquarePlanar,
                    StereoCoset::Undetermined,
                ),
                None,
            ),
            (
                StereoSite::Bond(BondId(0)),
                StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Undetermined),
                Some(Stereogenicity::Stereogenic),
            ),
            (
                StereoSite::Bond(BondId(1)),
                StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)),
                Some(Stereogenicity::Symmetric),
            ),
            (
                StereoSite::Atom(AtomId(4)),
                StereoConfigurationForm::Undetermined,
                None,
            ),
            (
                StereoSite::Bond(BondId(2)),
                StereoConfigurationForm::Undetermined,
                None,
            ),
        ];
        let selected: Vec<_> = candidates
            .iter()
            .filter(|(site, configuration, stereogenicity)| {
                config.mode.matches(configuration)
                    && config
                        .selection
                        .matches(*site, configuration.kind(), *stereogenicity)
            })
            .map(|(site, _, _)| *site)
            .collect();
        let expected = match mode {
            RemoveStereoMode::All => all,
            RemoveStereoMode::Undetermined => undetermined,
        };
        assert_eq!(selected, expected);
    }

    #[rstest]
    fn test_enumerate_stereo_config_default() {
        assert_eq!(
            EnumerateStereoConfig::default(),
            EnumerateStereoConfig {
                selection: StereoSelection {
                    sites: None,
                    kinds: Some(BTreeSet::from([
                        StereoKind::Tetrahedral,
                        StereoKind::CisTrans,
                    ])),
                    stereogenicities: None,
                },
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            },
        );
    }

    #[rstest]
    #[case::default(
        EnumerateStereoConfig::default(),
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2)),
            StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
    )]
    #[case::tetrahedral(
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::Tetrahedral])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(0)), StereoSite::Atom(AtomId(1)), StereoSite::Atom(AtomId(2))],
    )]
    #[case::cis_trans(
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))],
    )]
    #[case::sites(
        EnumerateStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(1)), StereoSite::Bond(BondId(1))])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(1)), StereoSite::Bond(BondId(1))],
    )]
    #[case::intersection(
        EnumerateStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(0)),
                    StereoSite::Bond(BondId(0)), StereoSite::Bond(BondId(1))])),
                kinds: Some(BTreeSet::from([StereoKind::CisTrans])),
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Stereogenic])),
            },
            ..Default::default()
        },
        vec![StereoSite::Bond(BondId(0))],
    )]
    #[case::prochiral(
        EnumerateStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Prochiral])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(1))],
    )]
    #[case::symmetric(
        EnumerateStereoConfig {
            selection: StereoSelection {
                stereogenicities: Some(BTreeSet::from([Stereogenicity::Symmetric])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![StereoSite::Atom(AtomId(2)), StereoSite::Bond(BondId(1))],
    )]
    #[case::empty_kinds(
        EnumerateStereoConfig {
            selection: StereoSelection {
                kinds: Some(BTreeSet::new()),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![],
    )]
    #[case::unbound_sites(
        EnumerateStereoConfig {
            selection: StereoSelection {
                sites: Some(BTreeSet::from([StereoSite::Atom(AtomId(u32::MAX)),
                    StereoSite::Bond(BondId(u32::MAX))])),
                ..EnumerateStereoConfig::default().selection
            },
            ..Default::default()
        },
        vec![],
    )]
    fn test_enumerate_stereo_config_selection(
        #[case] config: EnumerateStereoConfig,
        #[case] expected: Vec<StereoSite>,
    ) {
        let candidates = [
            (
                StereoSite::Atom(AtomId(0)),
                StereoKind::Tetrahedral,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Atom(AtomId(1)),
                StereoKind::Tetrahedral,
                Stereogenicity::Prochiral,
            ),
            (
                StereoSite::Atom(AtomId(2)),
                StereoKind::Tetrahedral,
                Stereogenicity::Symmetric,
            ),
            (
                StereoSite::Atom(AtomId(3)),
                StereoKind::SquarePlanar,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Bond(BondId(0)),
                StereoKind::CisTrans,
                Stereogenicity::Stereogenic,
            ),
            (
                StereoSite::Bond(BondId(1)),
                StereoKind::CisTrans,
                Stereogenicity::Symmetric,
            ),
        ];
        let selected: Vec<_> = candidates
            .into_iter()
            .filter(|&(site, kind, stereogenicity)| {
                config
                    .selection
                    .matches(site, Some(kind), Some(stereogenicity))
            })
            .map(|(site, _, _)| site)
            .collect();
        assert_eq!(selected, expected);
    }
}
