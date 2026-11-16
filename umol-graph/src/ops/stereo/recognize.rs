//! Recognition and classification of stereo sites.

use std::collections::BTreeSet;

use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AsLit, AtomConstraintForm, AtomConstraintKey, AtomForm, AtomId, BondConstraintForm,
    BondConstraintKey, BondId, CisTransStereoForm, Entity, GraphSymmetry, GraphSymmetryConfig,
    Molecule, Normalize, RingConfig, RingModel, RingSet, StereoCoset, StereoKind,
    StereogenicSiteKind, TetrahedralStereoForm,
};

use super::{StereoDerivation, StereoInconsistency, StereoPerceiver, StereoSelection, StereoSite};

/// A stereo query cannot interpret the selected input or requested processing.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum StereoPerceptionError {
    #[error("stereo selection references unavailable site {site:?}")]
    InvalidSite { site: StereoSite },
    #[error("site {site:?} has no stereo entity to classify")]
    MissingStereoEntity { site: StereoSite },
    #[error("stereo processing does not support kind {kind:?}")]
    UnsupportedKind { kind: StereoKind },
    #[error("stereo configuration at {site:?} is not supported")]
    UnsupportedConfiguration { site: StereoSite },
    #[error("stereo recognition requires literal fields at atom {atom:?}: {attributes}")]
    UnderdeterminedAtomField {
        atom: AtomId,
        attributes: Box<AtomForm>,
    },
    #[error("stereo recognition requires a literal order at bond {bond:?}")]
    UnderdeterminedBondOrder { bond: BondId },
    #[error(transparent)]
    Inconsistency(#[from] StereoInconsistency),
}

impl StereoPerceiver {
    /// Recognize selected missing tetrahedral atoms and cis/trans bonds.
    ///
    /// Positive source evidence requests the kind; new cosets are undetermined.
    /// Negative evidence excludes a site. Automatic local chemistry rules apply
    /// to unasserted sites; asserted sites use the structural rules of
    /// derive_stereo_atom and derive_stereo_bond.
    /// Existing entities are retained and checked for conflicting site evidence.
    /// Missing stereo bonds on double bonds follow the model's minimum ring size.
    /// Existing entities and source assertions are preserved.
    ///
    /// Stereogenic-site-kind filtering classifies all eligible prospective entities
    /// together with existing records in a temporary molecule, then applies the output
    /// selectors. Required recognition fields outside the selected output sites must
    /// therefore also be literal. With no site-kind filter, recognition processes only
    /// the selected sites and geometries. Empty sets select nothing. Supplied site ids
    /// are checked before processing.
    ///
    /// # Errors
    ///
    /// Reports unavailable selected sites, unsupported requested kinds,
    /// required non-literal input, and unsupported configurations
    /// when classification is requested. Evidence failures use the returned
    /// StereoDerivation's existing StereoInconsistency diagnostics.
    ///
    /// # Semantic properties
    ///
    /// The input is unchanged. Returned additions refer to its topology ids and
    /// have undetermined cosets in the recognized ligand order. Existing stereo
    /// entities are never returned as additions. Site, geometry and site-kind filters
    /// intersect, while classification retains the rest of the molecular context.
    pub fn recognize(
        &self,
        molecule: &Molecule,
        selection: &StereoSelection,
        automorphism_algorithm: AutomorphismAlgorithm,
        ring_config: RingConfig,
        max_classification_trials: u32,
    ) -> Result<StereoDerivation, StereoPerceptionError> {
        if let Some(sites) = &selection.sites {
            for &site in sites {
                check_site(molecule, site)?;
            }
        }
        if selection.sites.as_ref().is_some_and(BTreeSet::is_empty)
            || selection.kinds.as_ref().is_some_and(BTreeSet::is_empty)
            || selection
                .stereogenic_site_kinds
                .as_ref()
                .is_some_and(BTreeSet::is_empty)
        {
            return Ok(StereoDerivation::default());
        }
        for kind in [
            StereoKind::Tetrahedral,
            StereoKind::CisTrans,
            StereoKind::Axial,
            StereoKind::SquarePlanar,
            StereoKind::TrigonalBipyramidal,
            StereoKind::Octahedral,
        ] {
            let requested = match &selection.kinds {
                Some(kinds) => kinds.contains(&kind),
                None => self.model.kind_model(kind).is_some(),
            };
            if requested {
                check_kind(kind)?;
            }
        }

        let mut result = StereoDerivation::default();
        let mut atom_rings = None;
        let mut bond_rings = None;
        if self.model.kind_model(StereoKind::Tetrahedral).is_some()
            && (selection.stereogenic_site_kinds.is_some()
                || selection
                    .kinds
                    .as_ref()
                    .is_none_or(|kinds| kinds.contains(&StereoKind::Tetrahedral)))
        {
            for atom in molecule.atoms().iter() {
                let site = StereoSite::Atom(atom.id());
                if selection.stereogenic_site_kinds.is_none()
                    && selection
                        .sites
                        .as_ref()
                        .is_some_and(|sites| !sites.contains(&site))
                {
                    continue;
                }
                let evidence = atom
                    .constraints()
                    .evidenced(AtomConstraintKey::TetrahedralStereo);
                if let Some(existing) = molecule.stereo_atoms().at(atom.id()) {
                    if evidence.is_none() {
                        result.inconsistencies.push(
                            StereoInconsistency::TetrahedralStereoMismatch {
                                atom: atom.id(),
                                stereo_atom: existing.id(),
                            },
                        );
                    }
                    continue;
                }
                let asserted = match evidence {
                    Some(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo)) => {
                        true
                    }
                    Some(AtomConstraintForm::TetrahedralStereo(
                        TetrahedralStereoForm::NotStereo,
                    )) => continue,
                    Some(_) => false,
                    None => {
                        result.inconsistencies.push(
                            StereoInconsistency::TetrahedralStereoFailure { atom: atom.id() },
                        );
                        continue;
                    }
                };
                let eligible = asserted
                    || self.tetrahedral_eligible(
                        molecule,
                        atom.id(),
                        ring_config,
                        &mut atom_rings,
                    )?;
                if eligible {
                    require_frame_fields(molecule, atom.id())?;
                    let candidate = self
                        .derive_stereo_atom(molecule, atom.id(), &StereoCoset::Undetermined)
                        .filter(|(_, form)| form.configuration.clone().normalize().is_ok());
                    if let Some((ligands, form)) = candidate {
                        result.stereo_atoms.push((atom.id(), ligands, form));
                    } else if asserted {
                        result.inconsistencies.push(
                            StereoInconsistency::TetrahedralStereoFailure { atom: atom.id() },
                        );
                    }
                }
            }
        }
        if self.model.kind_model(StereoKind::CisTrans).is_some()
            && (selection.stereogenic_site_kinds.is_some()
                || selection
                    .kinds
                    .as_ref()
                    .is_none_or(|kinds| kinds.contains(&StereoKind::CisTrans)))
        {
            for bond in molecule.bonds().iter() {
                let site = StereoSite::Bond(bond.id());
                if selection.stereogenic_site_kinds.is_none()
                    && selection
                        .sites
                        .as_ref()
                        .is_some_and(|sites| !sites.contains(&site))
                {
                    continue;
                }
                let evidence = bond
                    .constraints()
                    .evidenced(BondConstraintKey::CisTransStereo);
                if let Some(existing) = molecule.stereo_bonds().at(bond.id()) {
                    if evidence.is_none() {
                        result
                            .inconsistencies
                            .push(StereoInconsistency::CisTransStereoMismatch {
                                bond: bond.id(),
                                stereo_bond: existing.id(),
                            });
                    }
                    continue;
                }
                let asserted = match evidence {
                    Some(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo)) => true,
                    Some(BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo)) => {
                        continue
                    }
                    Some(_) => false,
                    None => {
                        result
                            .inconsistencies
                            .push(StereoInconsistency::CisTransStereoFailure { bond: bond.id() });
                        continue;
                    }
                };
                let double = match bond.order().as_lit() {
                    Some(order) => order == 2,
                    None if asserted => false,
                    None => {
                        return Err(StereoPerceptionError::UnderdeterminedBondOrder {
                            bond: bond.id(),
                        })
                    }
                };
                if !asserted && !double {
                    continue;
                }
                if !self.stereo_bond_ring_eligible(
                    molecule,
                    bond.id(),
                    ring_config,
                    &mut bond_rings,
                ) {
                    continue;
                }
                let [first, second] = bond.atom_ids();
                let model = self
                    .model
                    .kind_model(StereoKind::CisTrans)
                    .expect("enabled cis/trans recognition");
                let mut in_scope = true;
                for atom in [first, second] {
                    let element = molecule
                        .atom(atom)
                        .element()
                        .as_lit()
                        .ok_or_else(|| underdetermined_atom_field(molecule, atom))?;
                    in_scope &= model.scope.contains(element);
                }
                if !in_scope {
                    if asserted {
                        result
                            .inconsistencies
                            .push(StereoInconsistency::CisTransStereoFailure { bond: bond.id() });
                    }
                    continue;
                }
                require_frame_fields(molecule, first)?;
                require_frame_fields(molecule, second)?;
                let candidate = self
                    .derive_stereo_bond(molecule, bond.id(), &StereoCoset::Undetermined)
                    .filter(|(_, form)| form.configuration.clone().normalize().is_ok());
                if let Some((ligands, form)) = candidate {
                    result.stereo_bonds.push((bond.id(), ligands, form));
                } else if asserted {
                    result
                        .inconsistencies
                        .push(StereoInconsistency::CisTransStereoFailure { bond: bond.id() });
                }
            }
        }
        if let Some(classes) = &selection.stereogenic_site_kinds {
            let sites: Vec<_> = result
                .stereo_atoms
                .iter()
                .map(|(site, _, _)| StereoSite::Atom(*site))
                .chain(
                    result
                        .stereo_bonds
                        .iter()
                        .map(|(site, _, _)| StereoSite::Bond(*site)),
                )
                .collect();
            if !sites.is_empty() {
                let mut editor = molecule.clone().edit();
                for (site, ligands, form) in &result.stereo_atoms {
                    editor.add_stereo_atom(*site, ligands, form.clone());
                }
                for (site, ligands, form) in &result.stereo_bonds {
                    editor.add_stereo_bond(*site, ligands, form.clone());
                }
                let prospective = editor
                    .finish()
                    .expect("recognized stereo additions preserve molecule integrity");
                let symmetry = prospective
                    .graph_symmetry(&GraphSymmetryConfig {
                        automorphism_algorithm,
                        max_classification_trials,
                        ..Default::default()
                    })
                    .into_graph_symmetry();
                let classifications = self.classify(&prospective, &sites, &symmetry)?;
                result.stereo_atoms.retain(|(site, _, _)| {
                    let index = classifications
                        .binary_search_by_key(&StereoSite::Atom(*site), |&(site, _)| site)
                        .expect("recognized site was classified");
                    classes.contains(&classifications[index].1)
                });
                result.stereo_bonds.retain(|(site, _, _)| {
                    let index = classifications
                        .binary_search_by_key(&StereoSite::Bond(*site), |&(site, _)| site)
                        .expect("recognized site was classified");
                    classes.contains(&classifications[index].1)
                });
            }
        }
        result.stereo_atoms.retain(|(site, _, form)| {
            selection
                .sites
                .as_ref()
                .is_none_or(|sites| sites.contains(&StereoSite::Atom(*site)))
                && selection.kinds.as_ref().is_none_or(|kinds| {
                    form.configuration
                        .kind()
                        .is_some_and(|kind| kinds.contains(&kind))
                })
        });
        result.stereo_bonds.retain(|(site, _, form)| {
            selection
                .sites
                .as_ref()
                .is_none_or(|sites| sites.contains(&StereoSite::Bond(*site)))
                && selection.kinds.as_ref().is_none_or(|kinds| {
                    form.configuration
                        .kind()
                        .is_some_and(|kind| kinds.contains(&kind))
                })
        });
        Ok(result)
    }

    /// Classify requested topology sites using the supplied graph-symmetry data.
    ///
    /// Input order is unrestricted and duplicate requests are coalesced. The result is sorted
    /// by atom site IDs followed by bond site IDs. Classification uses one complete kind batch,
    /// independently of the model's recognition scopes and stored coset values.
    /// An empty request returns an empty vector. The supplied data must describe this molecule.
    ///
    /// # Errors
    ///
    /// Reports unavailable topology sites, missing stereo records, unknown or misplaced kinds,
    /// and ligand frames whose degree differs from their kind. Stored cosets are not interpreted.
    ///
    /// # Semantic properties
    ///
    /// There is one entry per distinct requested site. The molecule and its exact representation
    /// are preserved. Valid stored coset assignments do not change the result.
    /// Relabeling topology sites transports the result keys without changing their kinds.
    /// Compatible ligand-frame actions preserve the kinds at the same sites.
    pub fn classify(
        &self,
        molecule: &Molecule,
        sites: &[StereoSite],
        symmetry: &GraphSymmetry,
    ) -> Result<Vec<(StereoSite, StereogenicSiteKind)>, StereoPerceptionError> {
        let mut sites = sites.to_vec();
        sites.sort_unstable();
        sites.dedup();
        let mut records = Vec::with_capacity(sites.len());
        for &site in &sites {
            check_site(molecule, site)?;
            let (entity, configuration, degree, failure) = match site {
                StereoSite::Atom(atom) => {
                    let stereo = molecule
                        .stereo_atoms()
                        .at(atom)
                        .ok_or(StereoPerceptionError::MissingStereoEntity { site })?;
                    (
                        Entity::StereoAtom(stereo.id()),
                        &stereo.attributes().configuration,
                        stereo.ligand_count(),
                        StereoInconsistency::StereoAtomFailure {
                            stereo_atom: stereo.id(),
                        },
                    )
                }
                StereoSite::Bond(bond) => {
                    let stereo = molecule
                        .stereo_bonds()
                        .at(bond)
                        .ok_or(StereoPerceptionError::MissingStereoEntity { site })?;
                    (
                        Entity::StereoBond(stereo.id()),
                        &stereo.attributes().configuration,
                        stereo.ligand_count(),
                        StereoInconsistency::StereoBondFailure {
                            stereo_bond: stereo.id(),
                        },
                    )
                }
            };
            let Some(kind) = configuration.kind() else {
                return Err(StereoPerceptionError::UnsupportedConfiguration { site });
            };
            if match site {
                StereoSite::Atom(_) => kind == StereoKind::CisTrans,
                StereoSite::Bond(_) => kind != StereoKind::CisTrans,
            } {
                return Err(StereoPerceptionError::UnsupportedConfiguration { site });
            }
            if degree != kind.degree() {
                return Err(failure.into());
            }
            records.push(entity);
        }
        if sites.is_empty() {
            return Ok(Vec::new());
        }
        let kinds = molecule.stereogenic_site_kinds(symmetry);
        Ok(sites
            .into_iter()
            .zip(records)
            .map(|(site, entity)| {
                let index = kinds
                    .binary_search_by_key(&entity, |&(entity, _)| entity)
                    .expect("the kinded requested record is included in the batch");
                (site, kinds[index].1)
            })
            .collect())
    }

    fn tetrahedral_eligible(
        &self,
        molecule: &Molecule,
        atom: AtomId,
        ring_config: RingConfig,
        rings: &mut Option<RingSet>,
    ) -> Result<bool, StereoPerceptionError> {
        let view = molecule.atom(atom);
        let neighbors = view.neighbors().len();
        if !(2..=4).contains(&neighbors) || view.is_in_aromatic_system() {
            return Ok(false);
        }
        let element = view
            .element()
            .as_lit()
            .ok_or_else(|| underdetermined_atom_field(molecule, atom))?;
        if !self
            .model
            .kind_model(StereoKind::Tetrahedral)
            .expect("enabled tetrahedral recognition")
            .scope
            .contains(element)
        {
            return Ok(false);
        }
        require_frame_fields(molecule, atom)?;
        let hydrogens = view
            .implicit_hydrogens()
            .as_lit()
            .expect("literal frame fields");
        let lone_pairs = view.lone_pairs().as_lit().expect("literal frame fields");
        if !(0..=1).contains(&hydrogens)
            || !(0..=1).contains(&lone_pairs)
            || neighbors as i64 + hydrogens + lone_pairs != 4
        {
            return Ok(false);
        }
        let mut valence = hydrogens;
        for neighbor in view.neighbors() {
            let order = neighbor.bond().order().as_lit().ok_or(
                StereoPerceptionError::UnderdeterminedBondOrder {
                    bond: neighbor.bond_id(),
                },
            )?;
            if !(1..=2).contains(&order) {
                return Ok(false);
            }
            valence += order;
        }
        let charge = view
            .charge()
            .as_lit()
            .ok_or_else(|| underdetermined_atom_field(molecule, atom))?;
        let bonding = neighbors as i64 + hydrogens;
        let eligible = match element {
            Element::B => charge == -1 && valence == 4 && bonding == 4 && lone_pairs == 0,
            Element::C | Element::Si | Element::Ge | Element::Sn => {
                charge == 0 && valence == 4 && bonding == 4 && lone_pairs == 0
            }
            Element::P | Element::As => {
                (charge == 0 && valence == 3 && bonding == 3 && lone_pairs == 1)
                    || (bonding == 4
                        && lone_pairs == 0
                        && (charge == 1 && valence == 4 || charge == 0 && valence == 5))
            }
            Element::S | Element::Se => {
                hydrogens == 0
                    && (bonding == 3
                        && lone_pairs == 1
                        && (charge == 0 && valence == 4 || charge == 1 && valence == 3)
                        || bonding == 4
                            && lone_pairs == 0
                            && (charge == 0 && valence == 6 || charge == 1 && valence == 5))
            }
            Element::N if bonding == 3 && lone_pairs == 1 && charge == 0 && valence == 3 => {
                let ids: Vec<_> = view
                    .neighbors()
                    .map(|neighbor| neighbor.atom_id())
                    .collect();
                let triangle = ids.iter().enumerate().any(|(i, &first)| {
                    ids[i + 1..].iter().any(|&second| {
                        molecule
                            .neighbors(first)
                            .any(|neighbor| neighbor.atom_id() == second)
                    })
                });
                let rigid = if triangle {
                    true
                } else if ids.len() < 3 || ids.iter().any(|&id| molecule.neighbors(id).len() < 2) {
                    false
                } else {
                    let ring_set = rings.get_or_insert_with(|| {
                        molecule
                            .rings(RingModel::default(), ring_config)
                            .into_ring_set()
                    });
                    let incident: Vec<_> = ring_set
                        .iter()
                        .filter(|ring| ring.atoms().contains(&atom))
                        .collect();
                    incident.iter().enumerate().any(|(i, first)| {
                        incident[i + 1..].iter().any(|second| {
                            ring_set.are_bridged(first.id, second.id)
                                && view
                                    .neighbors()
                                    .filter(|neighbor| {
                                        first.bonds().contains(&neighbor.bond_id())
                                            && second.bonds().contains(&neighbor.bond_id())
                                    })
                                    .count()
                                    == 1
                        })
                    })
                };
                rigid && !nitrogen_conjugated(molecule, atom)?
            }
            Element::N => {
                bonding == 4
                    && lone_pairs == 0
                    && (charge == 1 && valence == 4 || charge == 0 && valence == 5)
            }
            _ => bonding == 4 && lone_pairs == 0 && valence == 4,
        };
        Ok(eligible)
    }
}

fn check_site(molecule: &Molecule, site: StereoSite) -> Result<(), StereoPerceptionError> {
    let available = match site {
        StereoSite::Atom(atom) => molecule.atoms().get(atom).is_some(),
        StereoSite::Bond(bond) => molecule.bonds().get(bond).is_some(),
    };
    if available {
        Ok(())
    } else {
        Err(StereoPerceptionError::InvalidSite { site })
    }
}

fn check_kind(kind: StereoKind) -> Result<(), StereoPerceptionError> {
    if !matches!(kind, StereoKind::Tetrahedral | StereoKind::CisTrans) {
        return Err(StereoPerceptionError::UnsupportedKind { kind });
    }

    Ok(())
}

fn require_frame_fields(molecule: &Molecule, atom: AtomId) -> Result<(), StereoPerceptionError> {
    let view = molecule.atom(atom);
    if view.element().as_lit().is_none()
        || view.implicit_hydrogens().as_lit().is_none()
        || view.lone_pairs().as_lit().is_none()
    {
        return Err(underdetermined_atom_field(molecule, atom));
    }
    Ok(())
}

fn underdetermined_atom_field(molecule: &Molecule, atom: AtomId) -> StereoPerceptionError {
    StereoPerceptionError::UnderdeterminedAtomField {
        atom,
        attributes: Box::new(molecule.atom(atom).attributes().clone()),
    }
}

fn nitrogen_conjugated(molecule: &Molecule, atom: AtomId) -> Result<bool, StereoPerceptionError> {
    if !conjugation_eligible(molecule, atom)? {
        return Ok(false);
    }
    for neighbor in molecule.neighbors(atom) {
        let id = neighbor.atom_id();
        let view = molecule.atom(id);
        if view.is_in_aromatic_system() {
            return Ok(true);
        }
        if !(2..=3).contains(&view.neighbors().len()) || !conjugation_eligible(molecule, id)? {
            continue;
        }
        for other in view.neighbors() {
            if other.atom_id() == atom {
                continue;
            }
            let order = other.bond().order().as_lit().ok_or(
                StereoPerceptionError::UnderdeterminedBondOrder {
                    bond: other.bond_id(),
                },
            )?;
            if order > 1 && conjugation_eligible(molecule, other.atom_id())? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn conjugation_eligible(molecule: &Molecule, atom: AtomId) -> Result<bool, StereoPerceptionError> {
    let view = molecule.atom(atom);
    if view.neighbors().len() > 3 {
        return Ok(false);
    }
    let element = view
        .element()
        .as_lit()
        .ok_or_else(|| underdetermined_atom_field(molecule, atom))?;
    if !matches!(element.group_base18(), 2 | 13..=16) {
        return Ok(false);
    }
    let hydrogens = view
        .implicit_hydrogens()
        .as_lit()
        .ok_or_else(|| underdetermined_atom_field(molecule, atom))?;
    let count = view.neighbors().len() as i128 + i128::from(hydrogens);
    if !(0..=3).contains(&hydrogens)
        || count > 3
        || (element.period() > 2
            && matches!(element.group_base18(), 15 | 16)
            && !(element.group_base18() == 16 && count < 2))
    {
        return Ok(false);
    }
    let charge = view
        .charge()
        .as_lit()
        .ok_or_else(|| underdetermined_atom_field(molecule, atom))?;
    let unpaired = view
        .unpaired_electrons()
        .count
        .as_lit()
        .ok_or_else(|| underdetermined_atom_field(molecule, atom))?;
    let outer = i128::from(element.valence_electrons());
    let minimum = outer.min(8 - outer);
    let mut valence = i128::from(hydrogens);
    for neighbor in view.neighbors() {
        valence += i128::from(neighbor.bond().order().as_lit().ok_or(
            StereoPerceptionError::UnderdeterminedBondOrder {
                bond: neighbor.bond_id(),
            },
        )?);
    }
    Ok((charge != 0 || valence <= minimum)
        && minimum - count + (outer - minimum - i128::from(charge)).max(0) - i128::from(unpaired)
            > 0)
}

#[cfg(test)]
mod tests;
