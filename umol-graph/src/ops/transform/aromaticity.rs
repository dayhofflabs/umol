//! Transformations between localized and aromatic-system representations.

use std::collections::HashSet;
use std::convert::Infallible;
use std::iter;

use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_core::{
    BipartiteMaximumMatchingAlgorithm, GeneralMaximumMatchingAlgorithm, NonBipartiteGraphError,
};
use umol_graph_ir::ir::{
    AromaticSystemHandle, AromaticSystemId, AromaticSystemUpdate, AromaticValence,
    AromaticValenceForm, AsLit, AtomConstraintForm, AtomConstraintKey, AtomForm, AtomHandle,
    AtomId, AtomUpdate, BondConstraintKey, BondHandle, BondId, BondUpdate, Edit, Edits,
    ElectronCountsForm, ElementForm, Lattice, Molecule, MulticenterValenceForm, Normalize, NumForm,
    RingSet,
};
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{
    aromatic_valence_inconsistency, AromaticityConfig, AromaticityContradiction, AromaticityError,
    AromaticityInconsistency, AromaticityPerceiver,
};
use crate::ops::invariant::{ValenceInvariantsValidator, ValenceMismatch};
use crate::ops::model::AromaticityModel;
use crate::ops::transform::Transformer;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromatizeError {
    #[error("aromaticity setup: {0}")]
    Setup(#[from] AromaticityError),
    #[error("aromaticity contradiction: {0}")]
    Contradiction(#[from] AromaticityContradiction),
    #[error("aromatization input is underdetermined")]
    Underdetermined,
    #[error(
        "atom {atom:?} contradicts {constraint:?} in planned aromatic system {aromatic_system:?}"
    )]
    ContradictedAssertion {
        aromatic_system: AromaticSystemId,
        atom: AtomId,
        constraint: AtomConstraintForm,
    },
    #[error("planned aromatic system {aromatic_system:?}: {mismatch}")]
    ValenceInvariant {
        aromatic_system: AromaticSystemId,
        mismatch: ValenceMismatch,
    },
}

/// Transfers localized π electrons into perceived aromatic systems.
///
/// Charge, hydrogens, spin, external bonds, and unrelated state are preserved.
/// A molecule already containing an aromatic system is unchanged. Fresh perception
/// need not recover the same system membership as an earlier Kekulé localization.
///
/// Required nonliteral contribution or transfer data produces
/// [`AromatizeError::Underdetermined`]; the molecule need not be wholly ground.
/// Planned valence contradictions are reported before mutation, so an in-place
/// failure leaves the input unchanged.
///
/// # Semantic properties
///
/// On inputs without stored systems, every added system agrees with its members' open evidence
/// and contains no bond asserting non-aromaticity. Non-aromatic atom evidence
/// excludes that atom from perception. An incompatible selected contribution or
/// an unclaimed positive atom or bond assertion is an aromaticity contradiction.
#[derive(Clone, Debug)]
pub struct Aromatizer {
    perception: AromaticityPerceiver,
    config: AromaticityConfig,
}

impl Aromatizer {
    pub fn new(model: &AromaticityModel) -> Self {
        Self::with_config(model, AromaticityConfig::default())
    }

    pub fn with_config(model: &AromaticityModel, config: AromaticityConfig) -> Self {
        Self {
            perception: AromaticityPerceiver::new(model),
            config,
        }
    }

    fn plan_transform(&self, molecule: &Molecule) -> Result<Edits, AromatizeError> {
        let mut edits = Edits::new();
        if molecule.aromatic_systems().count() > 0 {
            return Ok(edits);
        }
        let rings = self.perception.candidate_rings(molecule, self.config);
        let mut contributions = Vec::with_capacity(if rings.count() == 0 {
            0
        } else {
            molecule.atoms().count()
        });
        let mut required_atoms = Vec::new();
        let mut mismatches = Vec::new();
        for atom in molecule.atoms().iter() {
            let evidence = atom
                .constraints()
                .evidenced(AtomConstraintKey::AromaticValence)
                .ok_or_else(|| {
                    AromaticityContradiction::from(aromatic_valence_inconsistency(
                        molecule,
                        atom.id(),
                    ))
                })?;
            if matches!(
                atom.attributes().constraints.aromatic_valence(),
                AromaticValenceForm::Aromatic(_)
            ) {
                required_atoms.push(atom.id());
            }
            if rings.count() == 0 {
                continue;
            }
            let contribution = if matches!(
                evidence,
                AtomConstraintForm::AromaticValence(AromaticValenceForm::NotAromatic)
            ) {
                AromaticValence::NotAromatic
            } else {
                match electrons_from_kekule(molecule, &rings, atom.id()) {
                    Err(AromatizeError::Underdetermined) if !rings.contains_atom(atom.id()) => {
                        AromaticValence::NotAromatic
                    }
                    result => result?,
                }
            };
            if let AromaticValence::Aromatic(electrons) = contribution {
                if !evidence.is_undetermined()
                    && evidence
                        .meet(&AtomConstraintForm::aromatic_valence(
                            AromaticValenceForm::aromatic(electrons),
                        ))
                        .is_none()
                {
                    mismatches.push(atom.id());
                }
            }
            contributions.push(contribution);
        }
        if rings.count() == 0 {
            if let Some(inconsistency) = self
                .perception
                .check_aromatic_systems(molecule, &required_atoms, &[])
                .into_iter()
                .next()
            {
                return Err(AromaticityContradiction::from(inconsistency).into());
            }
            return Ok(edits);
        }
        let derivation = self
            .perception
            .derive_from_rings(molecule, &rings, self.config, &required_atoms, |atom| {
                contributions[atom.index()]
            })?
            .into_decisive(AromatizeError::Underdetermined)?;
        if let Some(inconsistency) = derivation.inconsistencies.into_iter().next() {
            return Err(AromaticityContradiction::from(inconsistency).into());
        }
        for (index, (atoms, attributes)) in derivation.aromatic_systems.into_iter().enumerate() {
            let aromatic_system = AromaticSystemId(index as u32);
            if !mismatches.is_empty() {
                if let Some(&atom) = atoms
                    .iter()
                    .find(|atom| mismatches.binary_search(atom).is_ok())
                {
                    return Err(AromaticityContradiction::from(
                        AromaticityInconsistency::AromaticValenceMismatch {
                            atom,
                            aromatic_system,
                        },
                    )
                    .into());
                }
            }
            let ElectronCountsForm::Lit(electrons) = &attributes.electrons else {
                unreachable!("perception supplies literal aromatic contributions");
            };
            for (&id, &electrons) in atoms.iter().zip(electrons) {
                let atom = molecule.atom(id);
                let mut planned = atom.attributes().clone();
                planned
                    .constraints
                    .remove(AtomConstraintKey::AromaticValence);
                let mut valence = 0;
                for neighbor in atom.neighbors() {
                    let bond = neighbor.bond();
                    let mut order = bond.order().as_lit().expect(
                        "Kekulé contribution recovery requires literal incident bond orders",
                    );
                    if order == 2 && atoms.binary_search(&neighbor.atom_id()).is_ok() {
                        order = 1;
                        if id < neighbor.atom_id() {
                            edits.update_bond(
                                BondHandle::Id(bond.id()),
                                bond.attributes(),
                                &BondUpdate {
                                    order: Some(NumForm::Lit(1)),
                                    ..Default::default()
                                },
                            );
                        }
                    }
                    valence += order;
                }
                if electrons == 2 {
                    let pairs = atom
                        .lone_pairs()
                        .as_lit()
                        .ok_or(AromatizeError::Underdetermined)?;
                    planned.lone_pairs = NumForm::Lit(pairs - 1);
                    edits.update_atom(
                        AtomHandle::Id(id),
                        atom.attributes(),
                        &AtomUpdate {
                            lone_pairs: Some(planned.lone_pairs.clone()),
                            ..Default::default()
                        },
                    );
                }
                for key in [
                    AtomConstraintKey::Valence,
                    AtomConstraintKey::AromaticValence,
                    AtomConstraintKey::MulticenterValence,
                    AtomConstraintKey::DonatedPairs,
                    AtomConstraintKey::AcceptedPairs,
                ] {
                    let derived = match key {
                        AtomConstraintKey::Valence => Some(AtomConstraintForm::valence(valence)),
                        AtomConstraintKey::AromaticValence => {
                            Some(AtomConstraintForm::aromatic_valence(
                                AromaticValenceForm::aromatic(electrons),
                            ))
                        }
                        _ => atom.constraints().derived_complete(key),
                    }
                    .ok_or(AromatizeError::ValenceInvariant {
                        aromatic_system,
                        mismatch: ValenceMismatch::ContradictoryConstraint { atom_id: id, key },
                    })?;
                    let value = planned
                        .constraints
                        .asserted(key)
                        .and_then(|asserted| asserted.meet(&derived))
                        .ok_or_else(|| AromatizeError::ContradictedAssertion {
                            aromatic_system,
                            atom: id,
                            constraint: planned
                                .constraints
                                .get(key)
                                .cloned()
                                .unwrap_or_else(|| derived.clone()),
                        })?;
                    match value {
                        AtomConstraintForm::DonatedPairs(NumForm::Lit(0))
                        | AtomConstraintForm::AcceptedPairs(NumForm::Lit(0))
                        | AtomConstraintForm::MulticenterValence(
                            MulticenterValenceForm::NotMulticenter,
                        ) => {
                            planned.constraints.remove(key);
                        }
                        _ => planned.constraints.set(value),
                    }
                }
                if let Solution::Contradictory(mismatch) = ValenceInvariantsValidator
                    .validate_atom(&planned)
                    .unwrap_or_else(|error| match error {})
                {
                    let mismatch = match mismatch {
                        ValenceMismatch::OrbitalCount {
                            orbital_count,
                            electron_count,
                            ..
                        } => ValenceMismatch::OrbitalCount {
                            atom_id: id,
                            orbital_count,
                            electron_count,
                        },
                        ValenceMismatch::ContradictoryConstraint { key, .. } => {
                            ValenceMismatch::ContradictoryConstraint { atom_id: id, key }
                        }
                    };
                    return Err(AromatizeError::ValenceInvariant {
                        aromatic_system,
                        mismatch,
                    });
                }
            }
            for edit in self
                .perception
                .plan_aromatic_system(molecule, atoms, attributes)
            {
                edits.push(edit);
            }
        }
        Ok(edits)
    }
}

impl Transformer for Aromatizer {
    type Error = AromatizeError;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, AromatizeError> {
        let edits = self.plan_transform(&molecule)?;
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("aromatization plan applies to its input");
        Ok(editor
            .finish()
            .expect("aromatization plan preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), AromatizeError> {
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

/// Read localized contributions using the perception's candidate ring set.
/// Ring double bonds supply one electron; an external double bond remains
/// localized and leaves a carbon site with zero aromatic electrons.
/// Donors require an available pair; saturated and cumulated sites are excluded.
pub(crate) fn electrons_from_kekule(
    molecule: &Molecule,
    rings: &RingSet,
    atom: AtomId,
) -> Result<AromaticValence, AromatizeError> {
    let view = molecule.atom(atom);
    let ElementForm::Lit(element) = view.attributes().element else {
        return Err(AromatizeError::Underdetermined);
    };
    let mut double = None;
    for neighbor in view.neighbors() {
        match neighbor
            .bond()
            .order()
            .as_lit()
            .ok_or(AromatizeError::Underdetermined)?
        {
            1 => {}
            2 if double.is_none() => double = Some(neighbor.bond_id()),
            _ => return Ok(AromaticValence::NotAromatic),
        }
    }
    if let Some(double) = double {
        return Ok(
            if !rings.contains_atom(atom) || rings.contains_bond(double) {
                AromaticValence::Aromatic(1)
            } else if element == Element::C {
                AromaticValence::Aromatic(0)
            } else {
                AromaticValence::NotAromatic
            },
        );
    }
    let donor = match element {
        Element::C => match view
            .charge()
            .as_lit()
            .ok_or(AromatizeError::Underdetermined)?
        {
            1 => return Ok(AromaticValence::Aromatic(0)),
            -1 => true,
            _ => false,
        },
        Element::N | Element::O | Element::S | Element::Se | Element::P | Element::As => true,
        _ => false,
    };
    if donor {
        match view
            .lone_pairs()
            .as_lit()
            .ok_or(AromatizeError::Underdetermined)?
        {
            pairs if pairs > 0 => return Ok(AromaticValence::Aromatic(2)),
            _ => {}
        }
    }
    Ok(AromaticValence::NotAromatic)
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum KekulizeError {
    #[error("electron contributions are undetermined for aromatic system {0:?}")]
    UndeterminedElectrons(AromaticSystemId),
    #[error("charge is undetermined for aromatic system {0:?}")]
    UndeterminedCharge(AromaticSystemId),
    #[error("unpaired electrons are undetermined for aromatic system {0:?}")]
    UndeterminedUnpairedElectrons(AromaticSystemId),
    #[error(
        "aromatic system {aromatic_system:?} has {member_count} members but {electron_count} electron contributions"
    )]
    ElectronCountMismatch {
        aromatic_system: AromaticSystemId,
        member_count: usize,
        electron_count: usize,
    },
    #[error(
        "atom {atom:?} in aromatic system {aromatic_system:?} has unsupported electron contribution {contribution}"
    )]
    UnsupportedElectronContribution {
        aromatic_system: AromaticSystemId,
        atom: AtomId,
        contribution: i64,
    },
    #[error("aromatic system {0:?} is open-shell")]
    OpenShell(AromaticSystemId),
    #[error("aromatic system {aromatic_system:?} has unsupported charge {charge}")]
    UnsupportedCharge {
        aromatic_system: AromaticSystemId,
        charge: i64,
    },
    #[error("aromatic system {0:?} mixes prescribed and mobile exposed atoms")]
    MixedExposureDemand(AromaticSystemId),
    #[error(
        "node order for aromatic system {aromatic_system:?} has missing atoms {missing:?} and duplicate atoms {duplicates:?}"
    )]
    InvalidNodeOrder {
        aromatic_system: AromaticSystemId,
        missing: Vec<AtomId>,
        duplicates: Vec<AtomId>,
    },
    #[error(
        "maximum matching for aromatic system {aromatic_system:?} has deficiency {actual}; expected {expected}"
    )]
    MatchingDeficiency {
        aromatic_system: AromaticSystemId,
        expected: usize,
        actual: usize,
    },
    #[error(
        "charge is undetermined for exposed atom {atom:?} in aromatic system {aromatic_system:?}"
    )]
    UndeterminedExposedAtomCharge {
        aromatic_system: AromaticSystemId,
        atom: AtomId,
    },
    #[error("lone pairs are undetermined for exposed atom {atom:?} in aromatic system {aromatic_system:?}")]
    UndeterminedExposedAtomLonePairs {
        aromatic_system: AromaticSystemId,
        atom: AtomId,
    },
    #[error(
        "atom {atom:?} contradicts {constraint:?} in planned aromatic system {aromatic_system:?}"
    )]
    ContradictedAssertion {
        aromatic_system: AromaticSystemId,
        atom: AtomId,
        constraint: AtomConstraintForm,
    },
    #[error("planned aromatic system {aromatic_system:?}: {mismatch}")]
    ValenceInvariant {
        aromatic_system: AromaticSystemId,
        mismatch: ValenceMismatch,
    },
    #[error("no perfect matching exists for aromatic system {0:?}")]
    NoMatching(AromaticSystemId),
    #[error("maximum matching requires a bipartite aromatic system {0:?}")]
    NonBipartiteMatching(AromaticSystemId),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MatchingInputMode {
    Prescribed,
    OneMobileExposure,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct ExposedAtom {
    atom: AtomId,
    electrons: i64,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct MatchingInput {
    required_covered: Vec<AtomId>,
    required_exposed: Vec<ExposedAtom>,
    exposed_count: usize,
    mode: MatchingInputMode,
}

impl MatchingInput {
    fn from_aromatic_system(
        molecule: &Molecule,
        aromatic_system: AromaticSystemId,
    ) -> Result<Self, KekulizeError> {
        let view = molecule.aromatic_system(aromatic_system);
        let atoms: Vec<AtomId> = view.atom_ids().collect();
        let ElectronCountsForm::Lit(electrons) = view.electrons() else {
            return Err(KekulizeError::UndeterminedElectrons(aromatic_system));
        };
        if electrons.len() != atoms.len() {
            return Err(KekulizeError::ElectronCountMismatch {
                aromatic_system,
                member_count: atoms.len(),
                electron_count: electrons.len(),
            });
        }

        let NumForm::Lit(charge) = view.charge() else {
            return Err(KekulizeError::UndeterminedCharge(aromatic_system));
        };
        let (NumForm::Lit(unpaired_electrons), NumForm::Lit(multiplicity)) = (
            &view.unpaired_electrons().count,
            &view.unpaired_electrons().multiplicity,
        ) else {
            return Err(KekulizeError::UndeterminedUnpairedElectrons(
                aromatic_system,
            ));
        };
        if (*unpaired_electrons, *multiplicity) != (0, 1) {
            return Err(KekulizeError::OpenShell(aromatic_system));
        }
        if charge.abs() > 1 {
            return Err(KekulizeError::UnsupportedCharge {
                aromatic_system,
                charge: *charge,
            });
        }

        let mut required_covered = Vec::with_capacity(atoms.len());
        let mut required_exposed = Vec::new();
        for (&atom, &electrons) in atoms.iter().zip(electrons) {
            match electrons {
                1 => required_covered.push(atom),
                0 | 2 => required_exposed.push(ExposedAtom { atom, electrons }),
                contribution => {
                    return Err(KekulizeError::UnsupportedElectronContribution {
                        aromatic_system,
                        atom,
                        contribution,
                    });
                }
            }
        }

        if !required_exposed.is_empty() && *charge != 0 {
            return Err(KekulizeError::MixedExposureDemand(aromatic_system));
        }

        if required_exposed.is_empty() && *charge != 0 {
            Ok(Self {
                required_covered: Vec::new(),
                required_exposed,
                exposed_count: 1,
                mode: MatchingInputMode::OneMobileExposure,
            })
        } else {
            Ok(Self {
                exposed_count: required_exposed.len(),
                required_covered,
                required_exposed,
                mode: MatchingInputMode::Prescribed,
            })
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MaximumMatchingAlgorithm {
    Edmonds,
    HopcroftKarp,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KekulizeConfig {
    pub algorithm: MaximumMatchingAlgorithm,
}

impl KekulizeConfig {
    pub fn new(algorithm: MaximumMatchingAlgorithm) -> Self {
        Self { algorithm }
    }
}

impl Default for KekulizeConfig {
    fn default() -> Self {
        Self::new(MaximumMatchingAlgorithm::Edmonds)
    }
}

/// Localizes aromatic systems into single/double bonds and atom-local electrons.
///
/// Fixed two-electron donors regain one lone pair; fixed zero-electron acceptors
/// retain their local state. The remaining one-electron contributors must admit
/// a perfect matching. A delocalized system charge of ±1 is supported separately
/// when all contributions are one, by localizing it on one exposed atom.
///
/// Transformations report [`KekulizeError::NoMatching`] for an impossible fixed
/// assignment and [`KekulizeError::UndeterminedExposedAtomLonePairs`] when restoring
/// a donor pair requires a non-literal local count. They do not select another
/// hydrogen placement or contribution assignment. Planned incidence must agree
/// with retained valence assertions, and planned atoms must not contradict
/// electron conservation. Unchanged spin data and unrelated atoms are preserved
/// without whole-molecule validation.
///
/// The system entries and their atoms' and bonds' aromatic assertions are removed.
/// The matching algorithm is selected by [`KekulizeConfig`]. The atom processing
/// order is supplied at construction and determines which localization is chosen.
///
/// # Semantic properties
///
/// When the required readings are determined, every localized system atom
/// satisfies electron conservation and its retained valence assertions.
/// All errors occur before mutation: in-place failure leaves the input exactly
/// unchanged; consuming failure drops it. Unrelated undetermined state is retained.
#[derive(Clone, Debug)]
pub struct Kekulizer {
    config: KekulizeConfig,
    node_order: Vec<AtomId>,
}

impl Kekulizer {
    pub fn new(config: KekulizeConfig, node_order: Vec<AtomId>) -> Self {
        Self { config, node_order }
    }
}

impl Transformer for Kekulizer {
    type Error = KekulizeError;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, KekulizeError> {
        if molecule.aromatic_systems().count() == 0 {
            return Ok(molecule);
        }

        let edits = self.plan_transform(&molecule)?;
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("kekulization plan applies to its input");
        Ok(editor
            .finish()
            .expect("kekulization plan preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), KekulizeError> {
        if molecule.aromatic_systems().count() == 0 {
            return Ok(());
        }

        let edits = self.plan_transform(molecule)?;
        molecule.mutate(edits);
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        iter::once_with(move || self.transform(molecule.clone()).ok()).flatten()
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct AromaticSystemPlan {
    aromatic_system_idx: AromaticSystemId,
    atoms: Vec<AtomId>,
    // Follows atoms order; a mobile exposure has only one entry.
    exposed_atoms: Vec<ExposedAtom>,
    mobile_charge: Option<i64>,
    matched_bonds: Vec<BondId>,
    unmatched_bonds: Vec<BondId>,
}

impl Kekulizer {
    fn plan_transform(&self, molecule: &Molecule) -> Result<Edits, KekulizeError> {
        let plans = self.plan_aromatic_systems(molecule)?;
        let mut edits = Edits::new();
        let mut orders = vec![None; molecule.bonds().count()];
        for plan in &plans {
            debug_assert!(plan.matched_bonds.iter().all(|&bond| {
                molecule.bond(bond).atom_ids().iter().all(|atom| {
                    !plan
                        .exposed_atoms
                        .iter()
                        .any(|exposed| exposed.atom == *atom)
                })
            }));
            for (id, order) in plan
                .matched_bonds
                .iter()
                .map(|&id| (id, 2))
                .chain(plan.unmatched_bonds.iter().map(|&id| (id, 1)))
            {
                orders[id.index()] = Some(order);
                let bond = molecule.bond(id).attributes();
                edits.update_bond(
                    BondHandle::Id(id),
                    bond,
                    &BondUpdate {
                        order: Some(NumForm::Lit(order)),
                        ..Default::default()
                    },
                );
                if let Some(old) = bond.constraints.get(BondConstraintKey::Aromatic) {
                    edits.push(Edit::ModifyBondConstraint {
                        id: BondHandle::Id(id),
                        old: Some(old.clone()),
                        new: None,
                    });
                }
            }
            let mut exposed_atoms = plan.exposed_atoms.iter().peekable();
            for &id in &plan.atoms {
                let atom = molecule.atom(id);
                let attributes = atom.attributes();
                let mut planned = AtomForm {
                    element: attributes.element.clone(),
                    isotope_mass: attributes.isotope_mass.clone(),
                    charge: attributes.charge.clone(),
                    implicit_hydrogens: attributes.implicit_hydrogens.clone(),
                    lone_pairs: attributes.lone_pairs.clone(),
                    unpaired_electrons: attributes.unpaired_electrons.clone(),
                    constraints: Default::default(),
                };
                if let Some(old) = attributes
                    .constraints
                    .get(AtomConstraintKey::AromaticValence)
                {
                    edits.push(Edit::ModifyAtomConstraint {
                        id: AtomHandle::Id(id),
                        old: Some(old.clone()),
                        new: None,
                    });
                }
                if let Some(exposed) = exposed_atoms.next_if(|exposed| exposed.atom == id) {
                    let mut update = AtomUpdate::default();
                    if let Some(aromatic_system_charge) = plan.mobile_charge {
                        let NumForm::Lit(local_charge) = planned.charge else {
                            return Err(KekulizeError::UndeterminedExposedAtomCharge {
                                aromatic_system: plan.aromatic_system_idx,
                                atom: id,
                            });
                        };
                        planned.charge = NumForm::Lit(local_charge + aromatic_system_charge);
                        update.charge = Some(planned.charge.clone());
                    }
                    if exposed.electrons == 2 {
                        let NumForm::Lit(lone_pairs) = planned.lone_pairs else {
                            return Err(KekulizeError::UndeterminedExposedAtomLonePairs {
                                aromatic_system: plan.aromatic_system_idx,
                                atom: id,
                            });
                        };
                        planned.lone_pairs = NumForm::Lit(lone_pairs + 1);
                        update.lone_pairs = Some(planned.lone_pairs.clone());
                    }
                    edits.update_atom(AtomHandle::Id(id), atom.attributes(), &update);
                }
                let mut valence = NumForm::Lit(0);
                for neighbor in atom.neighbors() {
                    let order = match orders[neighbor.bond_id().index()] {
                        Some(order) => NumForm::Lit(order),
                        None => neighbor.bond().order().clone().normalize().map_err(|_| {
                            KekulizeError::ValenceInvariant {
                                aromatic_system: plan.aromatic_system_idx,
                                mismatch: ValenceMismatch::ContradictoryConstraint {
                                    atom_id: id,
                                    key: AtomConstraintKey::Valence,
                                },
                            }
                        })?,
                    };
                    valence = valence + order;
                }
                let mut determined_incidence = true;
                let has_dative = molecule.dative_bonds().has_incident(id);
                for key in [
                    AtomConstraintKey::Valence,
                    AtomConstraintKey::MulticenterValence,
                    AtomConstraintKey::DonatedPairs,
                    AtomConstraintKey::AcceptedPairs,
                ] {
                    let assertion = attributes.constraints.get(key);
                    let derived = match key {
                        AtomConstraintKey::Valence => {
                            Some(AtomConstraintForm::valence(valence.clone()))
                        }
                        AtomConstraintKey::MulticenterValence
                            if assertion.is_none() && !atom.is_in_multicenter_bond() =>
                        {
                            continue;
                        }
                        AtomConstraintKey::DonatedPairs | AtomConstraintKey::AcceptedPairs
                            if assertion.is_none() && !has_dative =>
                        {
                            continue;
                        }
                        _ => atom.constraints().derived_complete(key),
                    }
                    .ok_or(KekulizeError::ValenceInvariant {
                        aromatic_system: plan.aromatic_system_idx,
                        mismatch: ValenceMismatch::ContradictoryConstraint { atom_id: id, key },
                    })?;
                    let value = match assertion {
                        None => derived,
                        Some(constraint) => attributes
                            .constraints
                            .asserted(key)
                            .and_then(|asserted| asserted.meet(&derived))
                            .ok_or_else(|| KekulizeError::ContradictedAssertion {
                                aromatic_system: plan.aromatic_system_idx,
                                atom: id,
                                constraint: constraint.clone(),
                            })?,
                    };
                    determined_incidence &= !value.is_undetermined();
                    match value {
                        AtomConstraintForm::DonatedPairs(NumForm::Lit(0))
                        | AtomConstraintForm::AcceptedPairs(NumForm::Lit(0))
                        | AtomConstraintForm::MulticenterValence(
                            MulticenterValenceForm::NotMulticenter,
                        ) => {}
                        _ => planned.constraints.set(value),
                    }
                }
                // Bare-atom validation defaults absent incidence to zero.
                if determined_incidence {
                    if let Solution::Contradictory(mismatch) = ValenceInvariantsValidator
                        .validate_atom(&planned)
                        .unwrap_or_else(|error| match error {})
                    {
                        let mismatch = match mismatch {
                            ValenceMismatch::OrbitalCount {
                                orbital_count,
                                electron_count,
                                ..
                            } => ValenceMismatch::OrbitalCount {
                                atom_id: id,
                                orbital_count,
                                electron_count,
                            },
                            ValenceMismatch::ContradictoryConstraint { key, .. } => {
                                ValenceMismatch::ContradictoryConstraint { atom_id: id, key }
                            }
                        };
                        return Err(KekulizeError::ValenceInvariant {
                            aromatic_system: plan.aromatic_system_idx,
                            mismatch,
                        });
                    }
                }
            }
            debug_assert!(exposed_atoms.next().is_none());
        }
        if !plans.is_empty() {
            edits.remove_aromatic_systems(
                plans
                    .into_iter()
                    .map(|plan| {
                        (
                            AromaticSystemHandle::Id(plan.aromatic_system_idx),
                            plan.atoms.into_iter().map(AtomHandle::Id).collect(),
                            molecule
                                .aromatic_system(plan.aromatic_system_idx)
                                .attributes()
                                .clone(),
                        )
                    })
                    .collect(),
            );
        }
        Ok(edits)
    }

    /// Build the per-system matching plan against an immutable molecule.
    fn plan_aromatic_systems(
        &self,
        molecule: &Molecule,
    ) -> Result<Vec<AromaticSystemPlan>, KekulizeError> {
        let mut plans = Vec::with_capacity(molecule.aromatic_systems().count());
        for view in molecule.aromatic_systems().iter() {
            let aromatic_system_idx = view.id();
            let aromatic_system_atoms: Vec<AtomId> = view.atom_ids().collect();
            let bonds: Vec<BondId> = view.bond_ids().collect();
            let atom_set: HashSet<AtomId> = aromatic_system_atoms.iter().copied().collect();
            let mut seen = HashSet::with_capacity(aromatic_system_atoms.len());
            let mut ordered_host_atoms = Vec::with_capacity(aromatic_system_atoms.len());
            let mut duplicates = Vec::new();
            for atom in self
                .node_order
                .iter()
                .copied()
                .filter(|atom| atom_set.contains(atom))
            {
                if seen.insert(atom) {
                    ordered_host_atoms.push(atom);
                } else if !duplicates.contains(&atom) {
                    duplicates.push(atom);
                }
            }
            let missing: Vec<AtomId> = aromatic_system_atoms
                .iter()
                .copied()
                .filter(|atom| !seen.contains(atom))
                .collect();
            if !missing.is_empty() || !duplicates.is_empty() {
                return Err(KekulizeError::InvalidNodeOrder {
                    aromatic_system: aromatic_system_idx,
                    missing,
                    duplicates,
                });
            }

            let matching_input = MatchingInput::from_aromatic_system(molecule, view.id())?;
            let mobile_charge = match (matching_input.mode, view.charge()) {
                (MatchingInputMode::OneMobileExposure, NumForm::Lit(charge)) => Some(*charge),
                (MatchingInputMode::Prescribed, _) => None,
                (MatchingInputMode::OneMobileExposure, _) => {
                    unreachable!("matching input requires a literal system charge")
                }
            };
            let prescribed_exposed_set: HashSet<AtomId> = matching_input
                .required_exposed
                .iter()
                .map(|exposure| exposure.atom)
                .collect();
            let matching_host_atoms: Vec<AtomId> = match matching_input.mode {
                MatchingInputMode::Prescribed => ordered_host_atoms
                    .iter()
                    .copied()
                    .filter(|atom| !prescribed_exposed_set.contains(atom))
                    .collect(),
                MatchingInputMode::OneMobileExposure => ordered_host_atoms.clone(),
            };

            let correspondence_host_atoms: Vec<AtomId> = aromatic_system_atoms
                .iter()
                .copied()
                .filter(|atom| !prescribed_exposed_set.contains(atom))
                .collect();
            let correspondence = molecule.induced_subgraph(&correspondence_host_atoms);
            let extracted = molecule.extract(&correspondence);
            let sub_order: Vec<AtomId> = matching_host_atoms
                .iter()
                .map(|&host| {
                    correspondence
                        .atoms()
                        .left_of(host)
                        .expect("system atom maps to the extracted molecule")
                })
                .collect();
            let matching = match self.config.algorithm {
                MaximumMatchingAlgorithm::Edmonds => extracted
                    .graph()
                    .general_maximum_matching(&sub_order, GeneralMaximumMatchingAlgorithm::Edmonds),
                MaximumMatchingAlgorithm::HopcroftKarp => extracted
                    .graph()
                    .bipartite_maximum_matching(
                        &sub_order,
                        BipartiteMaximumMatchingAlgorithm::HopcroftKarp,
                    )
                    .map_err(|NonBipartiteGraphError| {
                        KekulizeError::NonBipartiteMatching(aromatic_system_idx)
                    })?,
            };
            let deficiency = extracted.atoms().count() - 2 * matching.size();
            match matching_input.mode {
                MatchingInputMode::Prescribed if deficiency != 0 => {
                    return Err(KekulizeError::NoMatching(aromatic_system_idx));
                }
                MatchingInputMode::OneMobileExposure
                    if deficiency != matching_input.exposed_count =>
                {
                    return Err(KekulizeError::MatchingDeficiency {
                        aromatic_system: aromatic_system_idx,
                        expected: matching_input.exposed_count,
                        actual: deficiency,
                    });
                }
                _ => {}
            }
            let exposed_atoms: Vec<ExposedAtom> = match matching_input.mode {
                MatchingInputMode::Prescribed => matching_input.required_exposed,
                MatchingInputMode::OneMobileExposure => extracted
                    .atoms()
                    .ids()
                    .filter(|&sub| !matching.is_matched(sub))
                    .map(|sub| ExposedAtom {
                        atom: correspondence
                            .atoms()
                            .right_of(sub)
                            .expect("extracted atom maps to the host molecule"),
                        electrons: if mobile_charge == Some(-1) { 2 } else { 0 },
                    })
                    .collect(),
            };
            let matched: HashSet<BondId> = matching
                .bonds()
                .map(|sub| {
                    correspondence
                        .bonds()
                        .right_of(sub)
                        .expect("matched bond maps to the host molecule")
                })
                .collect();
            let matched_atoms: HashSet<AtomId> = matched
                .iter()
                .flat_map(|&bond| molecule.bond(bond).atom_ids())
                .collect();
            let exposed_set: HashSet<AtomId> =
                exposed_atoms.iter().map(|exposed| exposed.atom).collect();
            let actual_exposed: HashSet<AtomId> = aromatic_system_atoms
                .iter()
                .copied()
                .filter(|atom| !matched_atoms.contains(atom))
                .collect();
            let required_covered: HashSet<AtomId> = match matching_input.mode {
                MatchingInputMode::Prescribed => {
                    matching_input.required_covered.iter().copied().collect()
                }
                MatchingInputMode::OneMobileExposure => aromatic_system_atoms
                    .iter()
                    .copied()
                    .filter(|atom| !exposed_set.contains(atom))
                    .collect(),
            };
            if exposed_atoms.len() != matching_input.exposed_count
                || matched_atoms != required_covered
                || actual_exposed != exposed_set
            {
                return Err(KekulizeError::NoMatching(aromatic_system_idx));
            }
            let (matched_bonds, unmatched_bonds): (Vec<BondId>, Vec<BondId>) =
                bonds.iter().copied().partition(|b| matched.contains(b));

            plans.push(AromaticSystemPlan {
                aromatic_system_idx,
                atoms: aromatic_system_atoms,
                exposed_atoms,
                mobile_charge,
                matched_bonds,
                unmatched_bonds,
            });
        }
        Ok(plans)
    }
}

/// Delocalizes atom-localized charge over homogeneous aromatic systems.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct DelocalizeCharge;

#[derive(Clone, Debug, PartialEq, Eq)]
struct DelocalizationPlan {
    aromatic_system: AromaticSystemId,
    charge: i64,
    electrons: Vec<i64>,
    atoms: Vec<(AtomId, i64)>,
}

impl DelocalizationPlan {
    fn derive(molecule: &Molecule, aromatic_system: AromaticSystemId) -> Option<Self> {
        let view = molecule.aromatic_system(aromatic_system);
        let atom_ids: Vec<AtomId> = view.atom_ids().collect();
        let (&first, rest) = atom_ids.split_first()?;
        let ElementForm::Lit(element) = molecule.atom(first).attributes().element else {
            return None;
        };
        if rest
            .iter()
            .any(|&atom| molecule.atom(atom).attributes().element != ElementForm::Lit(element))
        {
            return None;
        }

        let ElectronCountsForm::Lit(old_electrons) = &view.attributes().electrons else {
            return None;
        };
        if old_electrons.len() != atom_ids.len() {
            return None;
        }
        let NumForm::Lit(mut charge) = view.attributes().charge else {
            return None;
        };

        let mut electrons = Vec::with_capacity(atom_ids.len());
        let mut atoms = Vec::with_capacity(atom_ids.len());
        for atom_id in atom_ids {
            let atom = molecule.atom(atom_id);
            let NumForm::Lit(atom_charge) = atom.charge() else {
                return None;
            };
            let NumForm::Lit(degree) = atom.degree() else {
                return None;
            };
            let NumForm::Lit(implicit_hydrogens) = atom.implicit_hydrogens() else {
                return None;
            };
            let NumForm::Lit(lone_pairs) = atom.lone_pairs() else {
                return None;
            };
            let contribution = i64::from(element.valence_electrons())
                - degree
                - implicit_hydrogens
                - 2 * lone_pairs;
            charge += atom_charge;
            electrons.push(contribution);
            atoms.push((atom_id, contribution));
        }

        Some(Self {
            aromatic_system,
            charge,
            electrons,
            atoms,
        })
    }
}

impl DelocalizeCharge {
    fn plan_transform(&self, molecule: &Molecule) -> Edits {
        let mut edits = Edits::new();
        for aromatic_system in molecule.aromatic_systems().ids() {
            let Some(plan) = DelocalizationPlan::derive(molecule, aromatic_system) else {
                continue;
            };
            for (id, contribution) in plan.atoms {
                edits.update_atom(
                    AtomHandle::Id(id),
                    molecule.atom(id).attributes(),
                    &AtomUpdate {
                        charge: Some(NumForm::Lit(0)),
                        constraints: AtomConstraintForm::AromaticValence(
                            AromaticValenceForm::Aromatic(NumForm::Lit(contribution)),
                        )
                        .into(),
                        ..Default::default()
                    },
                );
            }
            edits.update_aromatic_system(
                AromaticSystemHandle::Id(plan.aromatic_system),
                molecule.aromatic_system(plan.aromatic_system).attributes(),
                &AromaticSystemUpdate {
                    charge: Some(NumForm::Lit(plan.charge)),
                    electrons: Some(ElectronCountsForm::Lit(plan.electrons)),
                    ..Default::default()
                },
            );
        }
        edits
    }
}

impl Transformer for DelocalizeCharge {
    type Error = Infallible;

    fn transform(&self, molecule: Molecule) -> Result<Molecule, Self::Error> {
        let edits = self.plan_transform(&molecule);
        let mut editor = molecule.edit();
        editor
            .apply(edits)
            .expect("charge delocalization plan applies to its input");
        Ok(editor
            .finish()
            .expect("charge delocalization plan preserves molecule integrity"))
    }

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error> {
        let edits = self.plan_transform(molecule);
        if !edits.is_empty() {
            molecule.mutate(edits);
        }
        Ok(())
    }

    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a {
        iter::once_with(move || match self.transform(molecule.clone()) {
            Ok(transformed) => transformed,
            Err(never) => match never {},
        })
    }
}

#[cfg(test)]
mod tests;
