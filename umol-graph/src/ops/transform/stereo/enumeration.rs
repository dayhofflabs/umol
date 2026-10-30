//! Ordered absolute domains and stereo assignment checks.

use std::collections::{BTreeMap, BTreeSet};

use umol_graph_core::{ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm};
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomConstraintKey, BondConstraintForm, BondConstraintKey, BooleanForm,
    CisTransStereoForm, ConstitutionColoring, Constraint, GraphSymmetry, GraphSymmetryConfig,
    Lattice, Molecule, Normalize, StereoAtomConstraintForm, StereoBondConstraintForm,
    StereoConfigurationForm, StereoCoset, StereoKind, StereoSymmetry, StereoTerm,
    StereogenicityForm, TetrahedralStereoForm, TopicityRelationForm,
};
use umol_perm::OrientedPermutation;
use umol_utils::solution::Solution;

use crate::ops::stereo::{StereoInconsistency, StereoPerceptionError, StereoSite};
use crate::ops::validate::constraint::{
    ConstraintInvariantsError, ConstraintInvariantsValidator, ConstraintValidateConfig, DerivedKind,
};

/// The caller supplies the selected configurations in processing order.
/// Domains and cursor storage grow with the number of sites, not their product.
pub(super) struct AbsoluteStereoAssignments {
    domains: Vec<(StereoSite, StereoKind, Vec<u32>)>,
    positions: Vec<usize>,
    exhausted: bool,
}

impl AbsoluteStereoAssignments {
    pub(super) fn new(
        configurations: impl IntoIterator<Item = (StereoSite, StereoConfigurationForm)>,
    ) -> Result<Self, StereoPerceptionError> {
        let mut domains = Vec::new();
        for (site, configuration) in configurations {
            let (kind, values) = absolute_domain(site, configuration)?;
            domains.push((site, kind, values));
        }
        let exhausted = domains.iter().any(|(_, _, domain)| domain.is_empty());
        Ok(Self {
            positions: vec![0; domains.len()],
            domains,
            exhausted,
        })
    }
}

fn absolute_domain(
    site: StereoSite,
    configuration: StereoConfigurationForm,
) -> Result<(StereoKind, Vec<u32>), StereoPerceptionError> {
    let StereoConfigurationForm::Kinded(kind, coset) = &configuration else {
        return Err(StereoPerceptionError::UnsupportedConfiguration { site });
    };
    let kind = *kind;
    if !matches!(kind, StereoKind::Tetrahedral | StereoKind::CisTrans) {
        return Err(StereoPerceptionError::UnsupportedKind { kind });
    }
    if !matches!(
        (site, kind),
        (StereoSite::Atom(_), StereoKind::Tetrahedral)
            | (StereoSite::Bond(_), StereoKind::CisTrans)
    ) {
        return Err(StereoPerceptionError::UnsupportedConfiguration { site });
    }
    let valid = match coset {
        StereoCoset::Undetermined => true,
        StereoCoset::Lit(value) => *value < 2,
        StereoCoset::LitSet(values) => values.iter().all(|&value| value < 2),
        StereoCoset::Term(term) => {
            let mut base = term.as_ref();
            loop {
                match base {
                    StereoTerm::Var(_) => break false,
                    StereoTerm::Lit(value) => break *value < 2,
                    StereoTerm::LitSet(values) => break values.iter().all(|&value| value < 2),
                    StereoTerm::Swap(inner)
                    | StereoTerm::Mirror(inner)
                    | StereoTerm::Apply(inner, _) => base = inner,
                }
            }
        }
    };
    if !valid {
        return Err(StereoPerceptionError::UnsupportedConfiguration { site });
    }
    let values = match configuration.normalize() {
        Err(_) => Vec::new(),
        Ok(StereoConfigurationForm::Kinded(_, StereoCoset::Undetermined)) => vec![0, 1],
        Ok(StereoConfigurationForm::Kinded(_, StereoCoset::Lit(value))) => vec![value],
        Ok(StereoConfigurationForm::Kinded(_, StereoCoset::LitSet(values))) => {
            values.into_iter().collect()
        }
        _ => return Err(StereoPerceptionError::UnsupportedConfiguration { site }),
    };
    Ok((kind, values))
}

impl Iterator for AbsoluteStereoAssignments {
    type Item = Vec<(StereoSite, StereoConfigurationForm)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }
        let assignment = self
            .domains
            .iter()
            .zip(&self.positions)
            .map(|((site, kind, values), &position)| {
                (
                    *site,
                    StereoConfigurationForm::Kinded(*kind, StereoCoset::Lit(values[position])),
                )
            })
            .collect();
        self.exhausted = true;
        for (position, (_, _, values)) in self.positions.iter_mut().zip(&self.domains).rev() {
            *position += 1;
            if *position < values.len() {
                self.exhausted = false;
                break;
            }
            *position = 0;
        }
        Some(assignment)
    }
}

/// Checks source assertions at processed sites and symmetry assertions that
/// assignments can affect, including assertions on retained stereo entities.
/// Candidates retain the source topology, entity ids, frames, and assertions;
/// only the chosen configuration values change.
pub(super) struct StereoAssignmentChecks {
    constraints: Vec<Constraint>,
}

impl StereoAssignmentChecks {
    pub(super) fn new(
        molecule: &Molecule,
        sites: &BTreeSet<StereoSite>,
    ) -> Result<Self, StereoPerceptionError> {
        let mut constraints = Vec::new();
        for &site in sites {
            check_frame(molecule, site)?;
            match site {
                StereoSite::Atom(atom) => {
                    if !molecule.stereo_atoms().is_at(atom) {
                        return Err(StereoPerceptionError::MissingStereoEntity { site });
                    }
                    if let Some(constraint) = molecule
                        .atom(atom)
                        .attributes()
                        .constraints
                        .get(AtomConstraintKey::TetrahedralStereo)
                    {
                        constraints.push(Constraint::Atom(atom, constraint.clone()));
                    }
                }
                StereoSite::Bond(bond) => {
                    if !molecule.stereo_bonds().is_at(bond) {
                        return Err(StereoPerceptionError::MissingStereoEntity { site });
                    }
                    if let Some(constraint) = molecule
                        .bond(bond)
                        .attributes()
                        .constraints
                        .get(BondConstraintKey::CisTransStereo)
                    {
                        constraints.push(Constraint::Bond(bond, constraint.clone()));
                    }
                }
            }
        }
        for stereo in molecule.stereo_atoms().iter() {
            if let Some(kind @ (StereoKind::Tetrahedral | StereoKind::CisTrans)) =
                stereo.attributes().configuration.kind()
            {
                constraints.extend(
                    stereo
                        .attributes()
                        .constraints
                        .iter()
                        .cloned()
                        .map(|constraint| Constraint::StereoAtom(stereo.id(), kind, constraint)),
                );
            }
        }
        for stereo in molecule.stereo_bonds().iter() {
            if let Some(kind @ (StereoKind::Tetrahedral | StereoKind::CisTrans)) =
                stereo.attributes().configuration.kind()
            {
                constraints.extend(
                    stereo
                        .attributes()
                        .constraints
                        .iter()
                        .cloned()
                        .map(|constraint| Constraint::StereoBond(stereo.id(), kind, constraint)),
                );
            }
        }
        constraints.extend(
            molecule
                .constraints()
                .iter()
                .filter(|constraint| affected_constraint(constraint, sites))
                .cloned(),
        );
        for constraint in &constraints {
            check_configuration_assertion(molecule, constraint)?;
        }
        Ok(Self { constraints })
    }

    pub(super) fn is_consistent(
        &self,
        molecule: &Molecule,
        symmetry_config: &GraphSymmetryConfig<ConstitutionColoring>,
    ) -> Result<bool, ConstraintInvariantsError> {
        let mut evaluation = StereoAssignmentEvaluation {
            molecule,
            symmetry_config,
            symmetry: None,
            local_symmetries: BTreeMap::new(),
            invariants: ConstraintInvariantsValidator::new(ConstraintValidateConfig {
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                derived_kind: DerivedKind::DerivedComplete,
            }),
        };
        let mut consistent = true;
        for constraint in &self.constraints {
            consistent &= evaluation.evaluate(constraint)? != Some(false);
        }
        Ok(consistent)
    }
}

fn check_frame(molecule: &Molecule, site: StereoSite) -> Result<(), StereoPerceptionError> {
    let (configuration, degree, failure) = match site {
        StereoSite::Atom(atom) => {
            if molecule.atoms().get(atom).is_none() {
                return Err(StereoPerceptionError::InvalidSite { site });
            }
            let Some(stereo) = molecule.stereo_atoms().at(atom) else {
                return Ok(());
            };
            (
                &stereo.attributes().configuration,
                stereo.ligand_ids().len(),
                StereoInconsistency::StereoAtomFailure {
                    stereo_atom: stereo.id(),
                },
            )
        }
        StereoSite::Bond(bond) => {
            if molecule.bonds().get(bond).is_none() {
                return Err(StereoPerceptionError::InvalidSite { site });
            }
            let Some(stereo) = molecule.stereo_bonds().at(bond) else {
                return Ok(());
            };
            (
                &stereo.attributes().configuration,
                stereo.ligand_ids().len(),
                StereoInconsistency::StereoBondFailure {
                    stereo_bond: stereo.id(),
                },
            )
        }
    };
    if configuration.kind().is_some_and(|kind| {
        matches!(kind, StereoKind::Tetrahedral | StereoKind::CisTrans) && degree != kind.degree()
    }) {
        return Err(failure.into());
    }
    Ok(())
}

fn check_configuration_assertion(
    molecule: &Molecule,
    constraint: &Constraint,
) -> Result<(), StereoPerceptionError> {
    match constraint {
        Constraint::Atom(atom, AtomConstraintForm::TetrahedralStereo(asserted)) => {
            check_frame(molecule, StereoSite::Atom(*atom))?;
            if let Some(stereo) = molecule.stereo_atoms().at(*atom) {
                let configuration = &stereo.attributes().configuration;
                if configuration.kind() == Some(StereoKind::Tetrahedral) {
                    absolute_domain(StereoSite::Atom(*atom), configuration.clone())?;
                }
            }
            if let TetrahedralStereoForm::Stereo(coset) = asserted {
                absolute_domain(
                    StereoSite::Atom(*atom),
                    StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset.clone()),
                )?;
            }
            if let Some(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(
                coset,
            ))) = molecule
                .atom(*atom)
                .attributes()
                .constraints
                .get(AtomConstraintKey::TetrahedralStereo)
            {
                absolute_domain(
                    StereoSite::Atom(*atom),
                    StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, coset.clone()),
                )?;
            }
        }
        Constraint::Bond(bond, BondConstraintForm::CisTransStereo(asserted)) => {
            check_frame(molecule, StereoSite::Bond(*bond))?;
            if let Some(stereo) = molecule.stereo_bonds().at(*bond) {
                let configuration = &stereo.attributes().configuration;
                if configuration.kind() == Some(StereoKind::CisTrans) {
                    absolute_domain(StereoSite::Bond(*bond), configuration.clone())?;
                }
            }
            if let CisTransStereoForm::Stereo(coset) = asserted {
                absolute_domain(
                    StereoSite::Bond(*bond),
                    StereoConfigurationForm::Kinded(StereoKind::CisTrans, coset.clone()),
                )?;
            }
            if let Some(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(coset))) =
                molecule
                    .bond(*bond)
                    .attributes()
                    .constraints
                    .get(BondConstraintKey::CisTransStereo)
            {
                absolute_domain(
                    StereoSite::Bond(*bond),
                    StereoConfigurationForm::Kinded(StereoKind::CisTrans, coset.clone()),
                )?;
            }
        }
        Constraint::And(children) | Constraint::Or(children) => {
            for child in children {
                check_configuration_assertion(molecule, child)?;
            }
        }
        Constraint::StereoAtom(id, _, _) => check_frame(
            molecule,
            StereoSite::Atom(molecule.stereo_atom(*id).site_id()),
        )?,
        Constraint::StereoBond(id, _, _) => check_frame(
            molecule,
            StereoSite::Bond(molecule.stereo_bond(*id).site_id()),
        )?,
        Constraint::Not(child) => check_configuration_assertion(molecule, child)?,
        _ => {}
    }
    Ok(())
}

fn affected_constraint(constraint: &Constraint, sites: &BTreeSet<StereoSite>) -> bool {
    match constraint {
        Constraint::Atom(atom, AtomConstraintForm::TetrahedralStereo(_)) => {
            sites.contains(&StereoSite::Atom(*atom))
        }
        Constraint::Bond(bond, BondConstraintForm::CisTransStereo(_)) => {
            sites.contains(&StereoSite::Bond(*bond))
        }
        Constraint::StereoAtom(_, StereoKind::Tetrahedral | StereoKind::CisTrans, _)
        | Constraint::StereoBond(_, StereoKind::Tetrahedral | StereoKind::CisTrans, _) => true,
        Constraint::And(children) | Constraint::Or(children) => children
            .iter()
            .any(|child| affected_constraint(child, sites)),
        Constraint::Not(child) => affected_constraint(child, sites),
        _ => false,
    }
}

struct StereoAssignmentEvaluation<'a> {
    molecule: &'a Molecule,
    symmetry_config: &'a GraphSymmetryConfig<ConstitutionColoring>,
    symmetry: Option<GraphSymmetry>,
    local_symmetries: BTreeMap<StereoSite, StereoSymmetry>,
    invariants: ConstraintInvariantsValidator,
}

impl StereoAssignmentEvaluation<'_> {
    fn evaluate(
        &mut self,
        constraint: &Constraint,
    ) -> Result<Option<bool>, ConstraintInvariantsError> {
        if constraint.is_vacuous() {
            return Ok(Some(true));
        }
        Ok(match constraint {
            Constraint::Atom(atom, asserted @ AtomConstraintForm::TetrahedralStereo(_)) => {
                evaluate_reading(
                    asserted,
                    self.molecule
                        .atom(*atom)
                        .constraints()
                        .evidenced_complete(AtomConstraintKey::TetrahedralStereo),
                )
            }
            Constraint::Bond(bond, asserted @ BondConstraintForm::CisTransStereo(_)) => {
                evaluate_reading(
                    asserted,
                    self.molecule
                        .bond(*bond)
                        .constraints()
                        .evidenced_complete(BondConstraintKey::CisTransStereo),
                )
            }
            Constraint::StereoAtom(id, kind, constraint) => {
                let stereo = self.molecule.stereo_atom(*id);
                if stereo.attributes().configuration.kind() != Some(*kind) {
                    return Ok(Some(false));
                }
                let symmetry = self.local_symmetry(StereoSite::Atom(stereo.site_id()));
                match constraint {
                    StereoAtomConstraintForm::Stereogenicity(asserted) => {
                        symmetry.map(|symmetry| {
                            asserted.matches(&StereogenicityForm::Lit(symmetry.stereogenicity()))
                        })
                    }
                    StereoAtomConstraintForm::Topicity(asserted) => symmetry.map(|symmetry| {
                        let pair = asserted.pair;
                        pair.first().index() < symmetry.group().degree()
                            && pair.second().index() < symmetry.group().degree()
                            && asserted.relation.matches(&TopicityRelationForm::Lit(
                                symmetry.topicity(pair.first(), pair.second()),
                            ))
                    }),
                    StereoAtomConstraintForm::LigandSymmetry(asserted) => {
                        symmetry.map(|symmetry| {
                            let op = OrientedPermutation::new(
                                asserted.permutation.permutation.0,
                                asserted.permutation.orientation,
                            );
                            op.degree() == symmetry.group().degree()
                                && asserted
                                    .invariant
                                    .matches(&BooleanForm::Lit(symmetry.group().contains(op)))
                        })
                    }
                    StereoAtomConstraintForm::Fluxionality(_) => None,
                }
            }
            Constraint::StereoBond(id, kind, constraint) => {
                let stereo = self.molecule.stereo_bond(*id);
                if stereo.attributes().configuration.kind() != Some(*kind) {
                    return Ok(Some(false));
                }
                let symmetry = self.local_symmetry(StereoSite::Bond(stereo.site_id()));
                match constraint {
                    StereoBondConstraintForm::Stereogenicity(asserted) => {
                        symmetry.map(|symmetry| {
                            asserted.matches(&StereogenicityForm::Lit(symmetry.stereogenicity()))
                        })
                    }
                    StereoBondConstraintForm::Topicity(asserted) => symmetry.map(|symmetry| {
                        let pair = asserted.pair;
                        pair.first().index() < symmetry.group().degree()
                            && pair.second().index() < symmetry.group().degree()
                            && asserted.relation.matches(&TopicityRelationForm::Lit(
                                symmetry.topicity(pair.first(), pair.second()),
                            ))
                    }),
                    StereoBondConstraintForm::LigandSymmetry(asserted) => {
                        symmetry.map(|symmetry| {
                            let op = OrientedPermutation::new(
                                asserted.permutation.permutation.0,
                                asserted.permutation.orientation,
                            );
                            op.degree() == symmetry.group().degree()
                                && asserted
                                    .invariant
                                    .matches(&BooleanForm::Lit(symmetry.group().contains(op)))
                        })
                    }
                    StereoBondConstraintForm::Fluxionality(_) => None,
                }
            }
            Constraint::And(children) | Constraint::Or(children) => {
                let is_and = matches!(constraint, Constraint::And(_));
                let mut decisive = false;
                let mut undetermined = false;
                for child in children {
                    let value = self.evaluate(child)?;
                    decisive |= value == Some(!is_and);
                    undetermined |= value.is_none();
                }
                if decisive {
                    Some(!is_and)
                } else if undetermined {
                    None
                } else {
                    Some(is_and || children.is_empty())
                }
            }
            Constraint::Not(child) => self.evaluate(child)?.map(|value| !value),
            _ => match self.invariants.evaluate(self.molecule, constraint)? {
                Solution::Determined(()) => Some(true),
                Solution::Underdetermined(()) => None,
                Solution::Contradictory(_) => Some(false),
            },
        })
    }

    fn local_symmetry(&mut self, site: StereoSite) -> Option<&StereoSymmetry> {
        if !self.local_symmetries.contains_key(&site) {
            let (kind, degree) = match site {
                StereoSite::Atom(atom) => {
                    let stereo = self.molecule.stereo_atoms().at(atom)?;
                    (
                        stereo.attributes().configuration.kind()?,
                        stereo.ligand_ids().len(),
                    )
                }
                StereoSite::Bond(bond) => {
                    let stereo = self.molecule.stereo_bonds().at(bond)?;
                    (
                        stereo.attributes().configuration.kind()?,
                        stereo.ligand_ids().len(),
                    )
                }
            };
            if !matches!(kind, StereoKind::Tetrahedral | StereoKind::CisTrans)
                || degree != kind.degree()
            {
                return None;
            }
            let symmetry = self
                .symmetry
                .get_or_insert_with(|| self.molecule.graph_symmetry(self.symmetry_config));
            let local = match site {
                StereoSite::Atom(atom) => self
                    .molecule
                    .stereo_atom_symmetry(symmetry, self.molecule.stereo_atoms().at(atom)?.id()),
                StereoSite::Bond(bond) => self
                    .molecule
                    .stereo_bond_symmetry(symmetry, self.molecule.stereo_bonds().at(bond)?.id()),
            };
            self.local_symmetries.insert(site, local);
        }
        self.local_symmetries.get(&site)
    }
}

fn evaluate_reading<C: Lattice>(asserted: &C, evidenced: Option<C>) -> Option<bool> {
    let Some(evidenced) = evidenced else {
        return Some(false);
    };
    let Ok(asserted) = asserted.clone().normalize() else {
        return Some(false);
    };
    if evidenced.satisfies(&asserted) {
        Some(true)
    } else if evidenced.is_compatible(&asserted) {
        None
    } else {
        Some(false)
    }
}

#[cfg(test)]
mod tests;
