//! Recognition and classification of stereo sites.

use std::collections::{BTreeMap, BTreeSet};

use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AsLit, AtomConstraintForm, AtomConstraintKey, AtomFieldKind, AtomId, BondConstraintForm,
    BondConstraintKey, BondId, CisTransStereoForm, ConstitutionColoring, GraphSymmetryConfig,
    Molecule, Normalize, RingConfig, RingModel, RingSet, StereoConfigurationForm, StereoCoset,
    StereoKind, Stereogenicity, TetrahedralStereoForm,
};

use super::{StereoDerivation, StereoInconsistency, StereoPerception, StereoSelection, StereoSite};

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
    #[error("fluxional stereo processing is not implemented for {kind:?}")]
    UnsupportedFluxionality { kind: StereoKind },
    #[error("stereo recognition requires a literal {field:?} at atom {atom:?}")]
    UnderdeterminedAtomField { atom: AtomId, field: AtomFieldKind },
    #[error("stereo recognition requires a literal order at bond {bond:?}")]
    UnderdeterminedBondOrder { bond: BondId },
    #[error(transparent)]
    Inconsistency(#[from] StereoInconsistency),
}

impl StereoPerception {
    /// Recognize selected missing tetrahedral atoms and cis/trans bonds.
    ///
    /// Positive source evidence supplies the configuration; otherwise the coset
    /// remains undetermined. Negative evidence excludes a site. Automatic local
    /// chemistry rules apply to unasserted sites; explicit arrangements retain
    /// the structural rules used by derive_stereo_atom and derive_stereo_bond.
    /// Existing entities are retained and checked for conflicting site evidence.
    ///
    /// Stereogenicity filtering uses molecular symmetry with prospective entities
    /// present in a temporary molecule. It is omitted when that filter is None.
    /// Empty sets select nothing. Supplied site ids are checked before processing.
    ///
    /// # Errors
    ///
    /// Reports unavailable selected sites, unsupported requested kinds or
    /// fluxionality, required non-literal input, and unsupported configurations
    /// when classification is requested. Evidence failures use the returned
    /// StereoDerivation's existing StereoInconsistency diagnostics.
    ///
    /// # Semantic properties
    ///
    /// The input is unchanged. Returned additions refer to its topology ids and
    /// preserve source cosets in the recognized ligand order. Existing stereo
    /// entities are never returned as additions. Site/kind/stereogenicity filters
    /// intersect, while classification retains the rest of the molecular context.
    pub fn recognize(
        &self,
        molecule: &Molecule,
        selection: &StereoSelection,
        automorphism_algorithm: AutomorphismAlgorithm,
        ring_config: RingConfig,
    ) -> Result<StereoDerivation, StereoPerceptionError> {
        if let Some(sites) = &selection.sites {
            for &site in sites {
                check_site(molecule, site)?;
            }
        }
        if selection.sites.as_ref().is_some_and(BTreeSet::is_empty)
            || selection.kinds.as_ref().is_some_and(BTreeSet::is_empty)
            || selection
                .stereogenicities
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
                check_kind(self, kind)?;
            }
        }

        let mut result = StereoDerivation::default();
        let mut rings = None;
        if self.model.kind_model(StereoKind::Tetrahedral).is_some()
            && selection
                .kinds
                .as_ref()
                .is_none_or(|kinds| kinds.contains(&StereoKind::Tetrahedral))
        {
            for atom in molecule.atoms().iter() {
                let site = StereoSite::Atom(atom.id());
                if selection
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
                let (coset, asserted) = match evidence {
                    Some(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(
                        coset,
                    ))) => (coset, true),
                    Some(AtomConstraintForm::TetrahedralStereo(
                        TetrahedralStereoForm::NotStereo,
                    )) => continue,
                    Some(_) => (StereoCoset::Undetermined, false),
                    None => {
                        result.inconsistencies.push(
                            StereoInconsistency::TetrahedralStereoFailure { atom: atom.id() },
                        );
                        continue;
                    }
                };
                let eligible = asserted
                    || self.tetrahedral_eligible(molecule, atom.id(), ring_config, &mut rings)?;
                if eligible {
                    require_frame_fields(molecule, atom.id())?;
                    let candidate = self
                        .derive_stereo_atom(molecule, atom.id(), &coset)
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
            && selection
                .kinds
                .as_ref()
                .is_none_or(|kinds| kinds.contains(&StereoKind::CisTrans))
        {
            for bond in molecule.bonds().iter() {
                let site = StereoSite::Bond(bond.id());
                if selection
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
                let (coset, asserted) = match evidence {
                    Some(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(coset))) => {
                        (coset, true)
                    }
                    Some(BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo)) => {
                        continue
                    }
                    Some(_) => (StereoCoset::Undetermined, false),
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
                if double && self.model.stereo_bond_minimum_ring_size > 0 {
                    let ring_set = self.recognition_rings(molecule, ring_config, &mut rings);
                    if ring_set
                        .bond_smallest_ring_size(bond.id())
                        .is_some_and(|size| {
                            size < self.model.stereo_bond_minimum_ring_size as usize
                        })
                    {
                        if asserted {
                            result.skipped_stereo_bonds.push(bond.id());
                        }
                        continue;
                    }
                }
                let [first, second] = bond.atom_ids();
                let model = self
                    .model
                    .kind_model(StereoKind::CisTrans)
                    .expect("enabled cis/trans recognition");
                let mut in_scope = true;
                for atom in [first, second] {
                    let element = molecule.atom(atom).element().as_lit().ok_or(
                        StereoPerceptionError::UnderdeterminedAtomField {
                            atom,
                            field: AtomFieldKind::Element,
                        },
                    )?;
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
                    .derive_stereo_bond(molecule, bond.id(), &coset)
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
        if let Some(classes) = &selection.stereogenicities {
            let sites: BTreeSet<_> = result
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
                let classifications =
                    self.classify(&prospective, &sites, automorphism_algorithm)?;
                result.stereo_atoms.retain(|(site, _, _)| {
                    classes.contains(&classifications[&StereoSite::Atom(*site)])
                });
                result.stereo_bonds.retain(|(site, _, _)| {
                    classes.contains(&classifications[&StereoSite::Bond(*site)])
                });
            }
        }
        Ok(result)
    }

    /// Classify existing stereo entities at the supplied topology sites.
    ///
    /// Uses the whole molecule's symmetry and the model's para-stereo setting.
    /// Automatic recognition rules and element scopes do not reject existing
    /// arrangements. Empty input returns an empty map without a symmetry query.
    ///
    /// # Errors
    ///
    /// Reports unavailable topology ids, missing stereo entities, unsupported
    /// kinds/configurations or fluxionality, and incompatible ligand arity.
    /// Classification supports Th/Ct with literal, undetermined, or
    /// nonempty finite absolute coset domains. Symbolic cosets are unsupported.
    ///
    /// # Semantic properties
    ///
    /// The returned keys equal the supplied sites. Classification does not
    /// change configurations or any other input data. It uses the projected
    /// ligand symmetry without choosing a coset.
    pub fn classify(
        &self,
        molecule: &Molecule,
        sites: &BTreeSet<StereoSite>,
        automorphism_algorithm: AutomorphismAlgorithm,
    ) -> Result<BTreeMap<StereoSite, Stereogenicity>, StereoPerceptionError> {
        for &site in sites {
            check_site(molecule, site)?;
            let (configuration, degree, failure) = match site {
                StereoSite::Atom(atom) => {
                    let record = molecule
                        .stereo_atoms()
                        .at(atom)
                        .ok_or(StereoPerceptionError::MissingStereoEntity { site })?;
                    (
                        &record.attributes().configuration,
                        record.ligand_ids().len(),
                        StereoInconsistency::StereoAtomFailure {
                            stereo_atom: record.id(),
                        },
                    )
                }
                StereoSite::Bond(bond) => {
                    let record = molecule
                        .stereo_bonds()
                        .at(bond)
                        .ok_or(StereoPerceptionError::MissingStereoEntity { site })?;
                    (
                        &record.attributes().configuration,
                        record.ligand_ids().len(),
                        StereoInconsistency::StereoBondFailure {
                            stereo_bond: record.id(),
                        },
                    )
                }
            };
            let StereoConfigurationForm::Kinded(kind, coset) = configuration else {
                return Err(StereoPerceptionError::UnsupportedConfiguration { site });
            };
            check_kind(self, *kind)?;
            if degree != kind.degree() {
                return Err(failure.into());
            }
            let supported = match coset {
                StereoCoset::Undetermined => true,
                StereoCoset::Lit(value) => *value < 2,
                StereoCoset::LitSet(values) => {
                    !values.is_empty() && values.iter().all(|&value| value < 2)
                }
                StereoCoset::Term(_) => false,
            };
            if !supported {
                return Err(StereoPerceptionError::UnsupportedConfiguration { site });
            }
        }
        if sites.is_empty() {
            return Ok(BTreeMap::new());
        }
        let symmetry = molecule.graph_symmetry(&GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: self.model.para_stereo,
            max_iterations: 16,
            automorphism_algorithm,
        });
        Ok(sites
            .iter()
            .map(|&site| {
                let class = match site {
                    StereoSite::Atom(atom) => molecule
                        .stereo_atom_symmetry(
                            &symmetry,
                            molecule
                                .stereo_atoms()
                                .at(atom)
                                .expect("validated stereo atom")
                                .id(),
                        )
                        .stereogenicity(),
                    StereoSite::Bond(bond) => molecule
                        .stereo_bond_symmetry(
                            &symmetry,
                            molecule
                                .stereo_bonds()
                                .at(bond)
                                .expect("validated stereo bond")
                                .id(),
                        )
                        .stereogenicity(),
                };
                (site, class)
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
        let element =
            view.element()
                .as_lit()
                .ok_or(StereoPerceptionError::UnderdeterminedAtomField {
                    atom,
                    field: AtomFieldKind::Element,
                })?;
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
        let charge =
            view.charge()
                .as_lit()
                .ok_or(StereoPerceptionError::UnderdeterminedAtomField {
                    atom,
                    field: AtomFieldKind::Charge,
                })?;
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
                    let ring_set = self.recognition_rings(molecule, ring_config, rings);
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

    fn recognition_rings<'a>(
        &self,
        molecule: &Molecule,
        config: RingConfig,
        rings: &'a mut Option<RingSet>,
    ) -> &'a RingSet {
        rings.get_or_insert_with(|| {
            let mut model = RingModel::default();
            if self.model.stereo_bond_minimum_ring_size > 0 {
                model.max_ring_size = model
                    .max_ring_size
                    .max(self.model.stereo_bond_minimum_ring_size as usize - 1);
            }
            molecule.rings(model, config).into_ring_set()
        })
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

fn check_kind(
    perception: &StereoPerception,
    kind: StereoKind,
) -> Result<(), StereoPerceptionError> {
    if !matches!(kind, StereoKind::Tetrahedral | StereoKind::CisTrans) {
        return Err(StereoPerceptionError::UnsupportedKind { kind });
    }
    if perception
        .model
        .kind_model(kind)
        .is_some_and(|model| model.fluxionality)
    {
        return Err(StereoPerceptionError::UnsupportedFluxionality { kind });
    }
    Ok(())
}

fn require_frame_fields(molecule: &Molecule, atom: AtomId) -> Result<(), StereoPerceptionError> {
    let view = molecule.atom(atom);
    for (field, literal) in [
        (AtomFieldKind::Element, view.element().as_lit().is_some()),
        (
            AtomFieldKind::ImplicitHydrogens,
            view.implicit_hydrogens().as_lit().is_some(),
        ),
        (
            AtomFieldKind::LonePairs,
            view.lone_pairs().as_lit().is_some(),
        ),
    ] {
        if !literal {
            return Err(StereoPerceptionError::UnderdeterminedAtomField { atom, field });
        }
    }
    Ok(())
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
    let element =
        view.element()
            .as_lit()
            .ok_or(StereoPerceptionError::UnderdeterminedAtomField {
                atom,
                field: AtomFieldKind::Element,
            })?;
    if !matches!(element.group_base18(), 2 | 13..=16) {
        return Ok(false);
    }
    let hydrogens = view.implicit_hydrogens().as_lit().ok_or(
        StereoPerceptionError::UnderdeterminedAtomField {
            atom,
            field: AtomFieldKind::ImplicitHydrogens,
        },
    )?;
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
        .ok_or(StereoPerceptionError::UnderdeterminedAtomField {
            atom,
            field: AtomFieldKind::Charge,
        })?;
    let unpaired = view.unpaired_electrons().count.as_lit().ok_or(
        StereoPerceptionError::UnderdeterminedAtomField {
            atom,
            field: AtomFieldKind::UnpairedElectrons,
        },
    )?;
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
