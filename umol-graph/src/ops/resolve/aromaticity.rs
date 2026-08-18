//! Aromaticity resolver. Perception reads aromatic valence constraints from the
//! materialized valence stage; planning emits complete aromatic systems and
//! localized bond aromatic constraints without mutating the source molecule.

use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet};

use smallvec::smallvec;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemHandle, AromaticSystemId, AromaticValenceForm, AsLit,
    AtomConstraintForm, AtomForm, AtomHandle, AtomId, AtomUpdate, BondConstraintForm, BondHandle,
    BondUpdate, BooleanForm, DativeBondConstraintForm, Edits, ElectronCountsForm, Molecule,
    NumForm, RingSet,
};
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{
    AromaticityConfig, AromaticityContradiction, AromaticityError, AromaticityInconsistency,
    AromaticityPerceiver,
};
use crate::ops::model::{AromaticityModel, AromaticityTieBreak, ValenceTieBreak};
use crate::ops::resolve::ResolveState;
use crate::ops::valence::compare::compare_by_key;

/// Per-component enumeration bound for assignments over aromatic-flexible
/// atoms; an exceeding component leaves the molecule underdetermined rather
/// than being sampled.
const MAX_ASSIGNMENTS: usize = 4096;

/// How aromaticity resolution handles an independently invalid constraint or entity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityFailurePolicy {
    Error,
    Keep,
}

/// How aromaticity resolution handles a valid aromatic-valence constraint that disagrees with a
/// valid aromatic system.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticityMismatchPolicy {
    Error,
    Keep,
    RemoveConstraint,
    ReplaceEntity,
}

/// How aromaticity resolution handles a valid localized-bond aromatic constraint that disagrees
/// with a valid aromatic system.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AromaticBondConstraintMismatchPolicy {
    Error,
    Keep,
    RemoveConstraint,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AromaticityResolveConfig {
    pub perception: AromaticityConfig,
    pub aromatic_valence_failure: AromaticityFailurePolicy,
    pub aromatic_system_failure: AromaticityFailurePolicy,
    pub aromatic_valence_mismatch: AromaticityMismatchPolicy,
    pub aromatic_bond_constraint_mismatch: AromaticBondConstraintMismatchPolicy,
    pub reset_aromatic_valence: bool,
}

impl Default for AromaticityResolveConfig {
    fn default() -> Self {
        Self {
            perception: AromaticityConfig::default(),
            aromatic_valence_failure: AromaticityFailurePolicy::Error,
            aromatic_system_failure: AromaticityFailurePolicy::Error,
            aromatic_valence_mismatch: AromaticityMismatchPolicy::Error,
            aromatic_bond_constraint_mismatch: AromaticBondConstraintMismatchPolicy::Error,
            reset_aromatic_valence: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct AromaticityResolver {
    perception: AromaticityPerceiver,
    tie_break: AromaticityTieBreak,
    config: AromaticityResolveConfig,
}

impl AromaticityResolver {
    pub fn new(model: &AromaticityModel) -> Self {
        Self::with_config(model, AromaticityResolveConfig::default())
    }

    pub fn with_config(model: &AromaticityModel, config: AromaticityResolveConfig) -> Self {
        Self {
            perception: AromaticityPerceiver::new(model),
            tie_break: model.tie_break,
            config,
        }
    }

    /// Construct the complete aromaticity edit plan without mutating `molecule`.
    pub fn plan(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<Edits, AromaticityContradiction>, AromaticityError> {
        let outcome = self.perception.derive(molecule, self.config.perception)?;

        match outcome {
            Solution::Determined(derivation) => {
                for &inconsistency in &derivation.inconsistencies {
                    let error = match inconsistency {
                        AromaticityInconsistency::AromaticValenceFailure { .. } => {
                            self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error
                        }
                        AromaticityInconsistency::AromaticSystemFailure { .. } => {
                            self.config.aromatic_system_failure == AromaticityFailurePolicy::Error
                        }
                        AromaticityInconsistency::AromaticValenceMismatch { .. } => {
                            self.config.aromatic_valence_mismatch
                                == AromaticityMismatchPolicy::Error
                        }
                        AromaticityInconsistency::AromaticBondConstraintMismatch { .. } => {
                            self.config.aromatic_bond_constraint_mismatch
                                == AromaticBondConstraintMismatchPolicy::Error
                        }
                    };
                    if error {
                        return Ok(Solution::Contradictory(inconsistency.into()));
                    }
                }

                let existing: BTreeSet<Vec<AtomId>> = molecule
                    .aromatic_systems()
                    .iter()
                    .map(|system| {
                        let mut atoms: Vec<AtomId> = system.atom_ids().collect();
                        atoms.sort_unstable();
                        atoms
                    })
                    .collect();

                let mut edits = Edits::new();
                let mut remove_constraints = BTreeSet::new();
                let mut remove_bond_constraints = BTreeSet::new();
                let mut replacements = BTreeSet::new();
                let mut suppressed = BTreeSet::new();

                for inconsistency in derivation.inconsistencies {
                    match inconsistency {
                        AromaticityInconsistency::AromaticValenceMismatch { atom, system } => {
                            let existing_members: BTreeSet<AtomId> =
                                molecule.aromatic_system(system).atom_ids().collect();
                            let candidate = derivation.systems.iter().position(|(atoms, _)| {
                                atoms.iter().copied().collect::<BTreeSet<_>>() == existing_members
                            });
                            match self.config.aromatic_valence_mismatch {
                                AromaticityMismatchPolicy::Error => unreachable!(),
                                AromaticityMismatchPolicy::Keep => {
                                    if let Some(candidate) = candidate {
                                        suppressed.insert(candidate);
                                    }
                                }
                                AromaticityMismatchPolicy::RemoveConstraint => {
                                    remove_constraints.insert(atom);
                                    if let Some(candidate) = candidate {
                                        suppressed.insert(candidate);
                                    }
                                }
                                AromaticityMismatchPolicy::ReplaceEntity => {
                                    if let Some(candidate) = candidate {
                                        replacements.insert((system, candidate));
                                    }
                                }
                            }
                        }
                        AromaticityInconsistency::AromaticBondConstraintMismatch {
                            bond, ..
                        } => match self.config.aromatic_bond_constraint_mismatch {
                            AromaticBondConstraintMismatchPolicy::Error => unreachable!(),
                            AromaticBondConstraintMismatchPolicy::Keep => {}
                            AromaticBondConstraintMismatchPolicy::RemoveConstraint => {
                                remove_bond_constraints.insert(bond);
                            }
                        },
                        AromaticityInconsistency::AromaticValenceFailure { .. }
                        | AromaticityInconsistency::AromaticSystemFailure { .. } => {}
                    }
                }

                if !replacements.is_empty() {
                    let removes = replacements
                        .iter()
                        .map(|&(system, _)| {
                            let view = molecule.aromatic_system(system);
                            (
                                AromaticSystemHandle::Id(system),
                                view.atom_ids().map(AtomHandle::Id).collect(),
                                view.attributes.clone(),
                            )
                        })
                        .collect();
                    edits.remove_aromatic_systems(removes);
                }

                for id in remove_constraints {
                    let mut update = AtomUpdate::default();
                    update.constraints.set(AtomConstraintForm::AromaticValence(
                        AromaticValenceForm::Undetermined,
                    ));
                    edits.update_atom(AtomHandle::Id(id), molecule.atom(id).attributes, &update);
                }
                for bond in remove_bond_constraints {
                    let mut update = BondUpdate::default();
                    update
                        .constraints
                        .set(BondConstraintForm::Aromatic(BooleanForm::Undetermined));
                    edits.update_bond(
                        BondHandle::Id(bond),
                        molecule.bond(bond).attributes,
                        &update,
                    );
                }

                let replaced_candidates: BTreeSet<usize> = replacements
                    .iter()
                    .map(|&(_, candidate)| candidate)
                    .collect();
                let replaced_entities: BTreeSet<AromaticSystemId> =
                    replacements.iter().map(|&(system, _)| system).collect();
                let retained_existing: BTreeSet<Vec<AtomId>> = molecule
                    .aromatic_systems()
                    .iter()
                    .filter(|system| !replaced_entities.contains(&system.id))
                    .map(|system| {
                        let mut atoms: Vec<AtomId> = system.atom_ids().collect();
                        atoms.sort_unstable();
                        atoms
                    })
                    .collect();

                for (candidate, (atoms, system)) in derivation.systems.into_iter().enumerate() {
                    let mut key = atoms.clone();
                    key.sort_unstable();
                    if replaced_candidates.contains(&candidate)
                        || (!suppressed.contains(&candidate)
                            && !existing.contains(&key)
                            && !retained_existing.contains(&key))
                    {
                        for edit in self.plan_system(molecule, atoms, system) {
                            edits.push(edit);
                        }
                    }
                }
                Ok(Solution::Determined(edits))
            }
            Solution::Underdetermined(_) => Ok(Solution::Underdetermined(Edits::new())),
            Solution::Contradictory(contradiction) => Ok(Solution::Contradictory(contradiction)),
        }
    }

    /// Plan and atomically apply aromaticity resolution.
    pub fn resolve(
        &self,
        molecule: &mut Molecule,
    ) -> Result<Solution<(), AromaticityContradiction>, AromaticityError> {
        let edits = match self.plan(molecule)? {
            Solution::Determined(edits) => edits,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => {
                return Ok(Solution::Contradictory(contradiction));
            }
        };
        let mut editor = molecule.edit();
        editor.transact(edits)?;
        *molecule = editor.build();
        Ok(Solution::Determined(()))
    }

    /// Selection among assignments per candidate-ring component, mutating
    /// nothing: an assignment is one completion choice per flexible atom of
    /// the component together with the systems perception finds under that
    /// narrowing; its restriction covers system members only. Contribution
    /// sourcing is uniform — an atom's carrier entry if present, else its
    /// stored input assertion (with the overlay-derived fallback). An
    /// assignment is valid iff every stored aromatic system touching the
    /// component reappears with the same member set, and — under the `Error`
    /// failure policy — no component atom whose every disjunct requires
    /// aromaticity is left unclaimed. Under the model's `MinElectronCount`
    /// tie-break, the structural order runs on the valid assignments first —
    /// claimed-atom count descending, then electron total ascending, then
    /// member-set lists lexicographic — and the chosen systems' members are
    /// recorded as tie-break uses when this order decided; under the model's
    /// `Strict` it never runs. A unique
    /// surviving restriction is accepted; several fall to `tie_break`
    /// member-wise. When the key leaves a tie, or when survivors restrict
    /// different atoms, the members stay plural and nothing is accepted.
    /// Survivors identical in restriction take the lexicographically smallest
    /// partition. The
    /// winner's systems are accepted as a whole — never mixed across
    /// assignments, so accepted systems are disjoint. Returns the narrowed
    /// carrier, the accepted systems, and the atoms selected by the key.
    ///
    /// More than `MAX_ASSIGNMENTS` assignments in one candidate-ring
    /// component, a non-literal
    /// stored `#a` outside the carrier, or an undetermined perception yields
    /// `Underdetermined` with the carrier unchanged. A carrier atom whose
    /// every disjunct requires aromaticity but which no accepted or tied
    /// system claims is `Contradictory`.
    ///
    /// # Semantic properties
    ///
    /// The outcome is search-independent: for carriers whose components stay
    /// within `MAX_ASSIGNMENTS`, the result (compared by `==` on the
    /// returned solution) equals that of a selection enumerating every
    /// assignment of every component exhaustively — the pruned search never
    /// removes a valid assignment. Cross-checked in the `property` test
    /// target against a definition-level flat enumeration over generated
    /// one- and two-ring Hückel scenarios under every policy combination.
    pub fn select(
        &self,
        molecule: &Molecule,
        state: ResolveState,
        tie_break: ValenceTieBreak,
    ) -> Result<Solution<ResolveState, AromaticityContradiction>, AromaticityError> {
        let ResolveState {
            mut completions,
            mut systems,
            mut tie_breaks,
        } = state;
        let carrier_atoms: BTreeSet<AtomId> = completions.iter().map(|(atom, _)| atom).collect();
        let stored_gate = molecule.atoms().iter().any(|atom| {
            !carrier_atoms.contains(&atom.id)
                && matches!(
                    atom.attributes.constraints.aromatic_valence(),
                    Some(AromaticValenceForm::Aromatic(valence)) if valence.as_lit().is_none()
                )
        });
        if stored_gate {
            return Ok(Solution::Underdetermined(ResolveState {
                completions,
                systems,
                tie_breaks,
            }));
        }

        let contribution = |form: &AtomForm| -> Option<u8> {
            match form.constraints.aromatic_valence() {
                Some(AromaticValenceForm::Aromatic(NumForm::Lit(valence))) => {
                    u8::try_from(*valence).ok()
                }
                _ => None,
            }
        };

        // Flexible: carrier atoms whose disjuncts differ in contribution.
        let flexible: Vec<(AtomId, Vec<Option<u8>>)> = completions
            .iter()
            .filter_map(|(atom, disjuncts)| {
                let contributions: Vec<Option<u8>> = disjuncts.iter().map(contribution).collect();
                (contributions.iter().collect::<BTreeSet<_>>().len() > 1)
                    .then_some((atom, contributions))
            })
            .collect();

        // The rule's acceptance couples atoms only within a candidate-ring
        // component: enumeration and the assignment bound are per component.
        // A flexible atom outside every component has no candidate ring and
        // falls through to the finalization tie-break.
        let rings = self
            .perception
            .candidate_rings(molecule, self.config.perception);
        let components = candidate_components(&rings, |atom_id| match completions.get(atom_id) {
            Some(disjuncts) => disjuncts.iter().any(|form| contribution(form).is_some()),
            None => stored_contribution(molecule, atom_id).is_some(),
        });
        for component in &components {
            let component_flexible: Vec<&(AtomId, Vec<Option<u8>>)> = flexible
                .iter()
                .filter(|(atom, _)| component.contains(atom))
                .collect();
            let assignment_count: usize = component_flexible
                .iter()
                .map(|(_, contributions)| contributions.len())
                .product();
            if assignment_count > MAX_ASSIGNMENTS {
                return Ok(Solution::Underdetermined(ResolveState {
                    completions,
                    systems,
                    tie_breaks,
                }));
            }
        }

        // An assignment per index choice: the partition (the perceived
        // systems inside the component, sorted by member list) and the
        // restriction (the chosen forms of flexible member atoms, ascending).
        // Validity and selection act on whole assignments, never mixing
        // systems across them, so accepted systems are disjoint.
        type Assignment = (
            Vec<(AtomId, AtomForm)>,
            Vec<(Vec<AtomId>, AromaticSystemForm)>,
        );
        let stored_systems: Vec<(AromaticSystemId, Vec<AtomId>)> = molecule
            .aromatic_systems()
            .iter()
            .map(|system| {
                let mut atoms: Vec<AtomId> = system.atom_ids().collect();
                atoms.sort_unstable();
                (system.id, atoms)
            })
            .collect();

        let claim_candidates = self.perception.claim_candidates(&rings);
        let mut accepted: Vec<(Vec<AtomId>, AromaticSystemForm)> = Vec::new();
        let mut tie_break_uses: BTreeSet<AtomId> = BTreeSet::new();
        let mut claimed: BTreeSet<AtomId> = BTreeSet::new();
        for component in &components {
            // Flexible atoms ordered ring-by-ring for the search; every
            // component atom lies in a candidate ring, so the order is total.
            let mut ordered_atoms: Vec<AtomId> = Vec::new();
            for ring in rings.iter() {
                if !ring.atoms().iter().all(|atom| component.contains(atom)) {
                    continue;
                }
                let mut atoms = ring.atoms().to_vec();
                atoms.sort_unstable();
                for atom in atoms {
                    if flexible
                        .iter()
                        .any(|&(flexible_atom, _)| flexible_atom == atom)
                        && !ordered_atoms.contains(&atom)
                    {
                        ordered_atoms.push(atom);
                    }
                }
            }
            let component_flexible: Vec<(AtomId, Vec<Option<u8>>)> = ordered_atoms
                .iter()
                .map(|&atom| {
                    flexible
                        .iter()
                        .find(|&&(flexible_atom, _)| flexible_atom == atom)
                        .cloned()
                        .expect("ordered atom is flexible")
                })
                .collect();
            let flexible_positions: BTreeMap<AtomId, usize> = component_flexible
                .iter()
                .enumerate()
                .map(|(position, &(atom, _))| (atom, position))
                .collect();
            // Atoms whose every disjunct requires aromaticity: the carrier
            // elimination criterion, shared by the search prunes and the
            // validity filter below.
            let aromatic_only: Vec<AtomId> = component
                .iter()
                .copied()
                .filter(|&atom| {
                    completions.get(atom).is_some_and(|disjuncts| {
                        disjuncts.iter().all(|form| {
                            matches!(
                                form.constraints.aromatic_valence(),
                                Some(AromaticValenceForm::Aromatic(_))
                            )
                        })
                    })
                })
                .collect();
            let fixed_contributions: BTreeMap<AtomId, Option<u8>> = component
                .iter()
                .filter(|atom| !flexible_positions.contains_key(atom))
                .map(|&atom| {
                    (
                        atom,
                        match completions.get(atom) {
                            Some(disjuncts) => contribution(&disjuncts[0]),
                            None => stored_contribution(molecule, atom),
                        },
                    )
                })
                .collect();
            let component_candidates: Vec<&Vec<AtomId>> = claim_candidates
                .iter()
                .filter(|members| members.iter().all(|atom| component.contains(atom)))
                .collect();

            // A candidate is settled-rejected under a partial assignment when
            // no completion can make the rule accept it: a member fixed to a
            // non-contribution, or a reachable total range the rule refuses.
            let settled_rejected = |members: &[AtomId], path: &[usize]| -> bool {
                let mut ranges: Vec<(u32, u32)> = Vec::with_capacity(members.len());
                for atom in members {
                    let range = match flexible_positions.get(atom) {
                        Some(&position) => {
                            let options = &component_flexible[position].1;
                            let assigned = (position < path.len())
                                .then(|| options[path[position]])
                                .map(|value| value.map(|value| (value, value)));
                            match assigned {
                                Some(value) => value,
                                None => {
                                    let values: Vec<u8> =
                                        options.iter().filter_map(|&value| value).collect();
                                    values
                                        .iter()
                                        .copied()
                                        .min()
                                        .zip(values.iter().copied().max())
                                }
                            }
                        }
                        None => fixed_contributions[atom].map(|value| (value, value)),
                    };
                    match range {
                        Some((lower, upper)) => {
                            ranges.push((u32::from(lower), u32::from(upper)));
                        }
                        None => return true,
                    }
                }
                !self.perception.accepts_range(&ranges)
            };
            // A partial assignment is certainly invalid when some
            // aromatic-only atom has every claim candidate settled-rejected.
            let certainly_invalid = |path: &[usize]| -> bool {
                aromatic_only.iter().any(|&atom| {
                    component_candidates
                        .iter()
                        .filter(|members| members.contains(&atom))
                        .all(|members| settled_rejected(members, path))
                })
            };
            let prune_enabled =
                self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error;

            // Depth-first search over the flexible atoms; a subtree is cut
            // only when every completion below it is certainly invalid, so
            // the valid-assignment set equals the flat enumeration's.
            let mut assignments: Vec<Assignment> = Vec::new();
            let mut path: Vec<usize> = Vec::new();
            let mut next = 0usize;
            loop {
                if path.len() == component_flexible.len() && next == 0 {
                    let choice: BTreeMap<AtomId, usize> = component_flexible
                        .iter()
                        .zip(&path)
                        .map(|(&(atom, _), &index)| (atom, index))
                        .collect();
                    let outcome = self.perception.find_systems_from_rings(
                        molecule,
                        &rings,
                        self.config.perception,
                        |atom| match choice.get(&atom) {
                            Some(&index) => {
                                contribution(&completions.get(atom).expect("flexible atom")[index])
                            }
                            None => match completions.get(atom) {
                                Some(disjuncts) => contribution(&disjuncts[0]),
                                None => stored_contribution(molecule, atom),
                            },
                        },
                    )?;
                    let found = match outcome {
                        Solution::Determined(found) => found,
                        Solution::Underdetermined(_) => {
                            return Ok(Solution::Underdetermined(ResolveState {
                                completions,
                                systems,
                                tie_breaks,
                            }));
                        }
                        Solution::Contradictory(contradiction) => {
                            return Ok(Solution::Contradictory(contradiction));
                        }
                    };
                    // Systems outside the component arise under this
                    // component's default indices; their own enumeration
                    // accumulates them.
                    let mut partition: Vec<(Vec<AtomId>, AromaticSystemForm)> = found
                        .into_iter()
                        .filter(|(atoms, _)| atoms.iter().all(|atom| component.contains(atom)))
                        .collect();
                    partition.sort_by(|(a, _), (b, _)| a.cmp(b));
                    let restriction: Vec<(AtomId, AtomForm)> = choice
                        .iter()
                        .filter(|(atom, _)| partition.iter().any(|(atoms, _)| atoms.contains(atom)))
                        .map(|(&atom, &index)| {
                            (
                                atom,
                                completions.get(atom).expect("flexible atom")[index].clone(),
                            )
                        })
                        .collect();
                    let assignment = (restriction, partition);
                    if !assignments.contains(&assignment) {
                        assignments.push(assignment);
                    }
                    match path.pop() {
                        Some(index) => next = index + 1,
                        None => break,
                    }
                    continue;
                }
                if next >= component_flexible[path.len()].1.len() {
                    match path.pop() {
                        Some(index) => next = index + 1,
                        None => break,
                    }
                    continue;
                }
                path.push(next);
                if prune_enabled && certainly_invalid(&path) {
                    let index = path.pop().expect("pushed above");
                    next = index + 1;
                } else {
                    next = 0;
                }
            }

            // Validity under the `Error` failure policy first — the search
            // prunes by the same criterion: an assignment may not leave an
            // aromatic-only component atom unclaimed.
            let mut valid = assignments;
            if prune_enabled {
                valid.retain(|(_, partition)| {
                    aromatic_only
                        .iter()
                        .all(|atom| partition.iter().any(|(atoms, _)| atoms.contains(atom)))
                });
            }
            if valid.is_empty() {
                continue;
            }
            // Validity: every stored system touching the component must
            // reappear with the same member set; when no carrier-valid
            // assignment reproduces one, the failure policy decides between
            // contradiction and an inert component.
            for (system, members) in stored_systems
                .iter()
                .filter(|(_, members)| members.iter().any(|atom| component.contains(atom)))
            {
                valid.retain(|(_, partition)| partition.iter().any(|(atoms, _)| atoms == members));
                if valid.is_empty() {
                    if self.config.aromatic_system_failure == AromaticityFailurePolicy::Error {
                        return Ok(Solution::Contradictory(
                            AromaticityInconsistency::AromaticSystemFailure { system: *system }
                                .into(),
                        ));
                    }
                    break;
                }
            }
            if valid.is_empty() {
                continue;
            }
            for (_, partition) in &valid {
                for (atoms, _) in partition {
                    claimed.extend(atoms.iter().copied());
                }
            }

            // Structural order under `MinElectronCount`: most claimed atoms
            // first, then the smallest electron total, member-set lists
            // breaking ties lexicographically. The electron component acts
            // whenever the rule admits more than one total for the same
            // members; the coverage component only under a non-`Error`
            // failure policy.
            let mut structural_decided = false;
            if self.tie_break == AromaticityTieBreak::MinElectronCount {
                let structure = |partition: &Vec<(Vec<AtomId>, AromaticSystemForm)>| {
                    (
                        Reverse(
                            partition
                                .iter()
                                .map(|(atoms, _)| atoms.len())
                                .sum::<usize>(),
                        ),
                        partition
                            .iter()
                            .map(|(_, form)| match &form.electrons {
                                ElectronCountsForm::Lit(electrons) => electrons.iter().sum(),
                                _ => i64::MAX,
                            })
                            .sum::<i64>(),
                        partition
                            .iter()
                            .map(|(atoms, _)| atoms.clone())
                            .collect::<Vec<_>>(),
                    )
                };
                let best = valid
                    .iter()
                    .map(|(_, partition)| structure(partition))
                    .min()
                    .expect("non-empty survivors");
                structural_decided = valid
                    .iter()
                    .any(|(_, partition)| structure(partition) != best);
                valid.retain(|(_, partition)| structure(partition) == best);
            }

            // Selection: a unique restriction is accepted; several compare by
            // the value key when they restrict the same atoms; otherwise the
            // members stay plural. Restriction-identical survivors take the
            // lexicographically smallest partition — representation
            // canonicalization, not policy.
            let mut restrictions: Vec<&Vec<(AtomId, AtomForm)>> = Vec::new();
            for (restriction, _) in &valid {
                if !restrictions.contains(&restriction) {
                    restrictions.push(restriction);
                }
            }
            let (winner_restriction, by_key) = if restrictions.len() == 1 {
                (restrictions[0], false)
            } else {
                if tie_break.key().is_empty() {
                    continue;
                }
                let domain: BTreeSet<AtomId> =
                    restrictions[0].iter().map(|(atom, _)| *atom).collect();
                if !restrictions.iter().all(|restriction| {
                    restriction
                        .iter()
                        .map(|(atom, _)| *atom)
                        .collect::<BTreeSet<_>>()
                        == domain
                }) {
                    continue;
                }
                let best = restrictions
                    .iter()
                    .copied()
                    .max_by(|a, b| compare_restrictions(a, b, tie_break))
                    .expect("non-empty restrictions");
                let unique = restrictions
                    .iter()
                    .filter(|restriction| {
                        compare_restrictions(restriction, best, tie_break).is_eq()
                    })
                    .count()
                    == 1;
                if !unique {
                    continue;
                }
                (best, true)
            };
            let (_, partition) = valid
                .iter()
                .filter(|(restriction, _)| restriction == winner_restriction)
                .min_by(|(_, a), (_, b)| {
                    a.iter()
                        .map(|(atoms, _)| atoms)
                        .cmp(b.iter().map(|(atoms, _)| atoms))
                })
                .expect("winner restriction present");
            accepted.extend(partition.iter().cloned());
            if structural_decided {
                for (atoms, _) in partition {
                    tie_break_uses.extend(atoms.iter().copied());
                }
            }
            for (atom, form) in winner_restriction {
                if by_key && completions.get(*atom).is_some_and(|entry| entry.len() > 1) {
                    tie_break_uses.insert(*atom);
                }
                completions.insert(*atom, smallvec![form.clone()]);
            }
        }

        // A carrier atom whose every disjunct requires aromaticity but which
        // no accepted or tied system claims cannot be completed; the failure
        // policy decides between contradiction and keeping the assertion.
        if self.config.aromatic_valence_failure == AromaticityFailurePolicy::Error {
            for (atom, disjuncts) in completions.iter() {
                if claimed.contains(&atom) {
                    continue;
                }
                if disjuncts.iter().all(|form| {
                    matches!(
                        form.constraints.aromatic_valence(),
                        Some(AromaticValenceForm::Aromatic(_))
                    )
                }) {
                    return Ok(Solution::Contradictory(
                        AromaticityInconsistency::AromaticValenceFailure { atom }.into(),
                    ));
                }
            }
        }

        systems.extend(accepted);
        tie_breaks.extend(tie_break_uses);
        tie_breaks.sort_unstable();
        tie_breaks.dedup();
        Ok(Solution::Determined(ResolveState {
            completions,
            systems,
            tie_breaks,
        }))
    }

    pub(crate) fn plan_system(
        &self,
        molecule: &Molecule,
        atoms: Vec<AtomId>,
        system: AromaticSystemForm,
    ) -> Edits {
        let mut atom_updates = Vec::new();
        if self.config.reset_aromatic_valence {
            for &atom_id in &atoms {
                let mut update = AtomUpdate::default();
                update.constraints.set(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Undetermined,
                ));
                atom_updates.push((atom_id, update));
            }
        }

        let mut edits = Edits::new();
        edits.add_aromatic_system(atoms.iter().copied().map(AtomHandle::Id).collect(), system);
        for (atom_id, update) in atom_updates {
            edits.update_atom(
                AtomHandle::Id(atom_id),
                molecule.atom(atom_id).attributes,
                &update,
            );
        }

        let members: BTreeSet<AtomId> = atoms.iter().copied().collect();
        let mut bond_ids = BTreeSet::new();
        for &atom_id in &atoms {
            for neighbor in molecule.atom(atom_id).neighbors() {
                if members.contains(&neighbor.atom_id()) {
                    bond_ids.insert(neighbor.bond_id());
                }
            }
        }
        for bond_id in bond_ids {
            if matches!(
                molecule.bond(bond_id).attributes.constraints.aromatic(),
                BooleanForm::Lit(_)
            ) {
                continue;
            }
            let mut update = BondUpdate::default();
            update
                .constraints
                .set(BondConstraintForm::Aromatic(BooleanForm::Lit(true)));
            edits.update_bond(
                BondHandle::Id(bond_id),
                molecule.bond(bond_id).attributes,
                &update,
            );
        }
        edits
    }
}

/// Member-wise lexicographic comparison of two assignment restrictions for
/// the same system, in ascending member order, each member compared by the
/// tie-break key over its candidate forms.
/// The contribution of an atom outside the carrier: a literal stored
/// assertion, else — with no assertion opinion — the stored aromatic
/// system's literal electron count.
fn stored_contribution(molecule: &Molecule, atom: AtomId) -> Option<u8> {
    let view = molecule.atom(atom);
    match view.attributes.constraints.aromatic_valence() {
        Some(AromaticValenceForm::Aromatic(NumForm::Lit(valence))) => u8::try_from(*valence).ok(),
        Some(AromaticValenceForm::Aromatic(_) | AromaticValenceForm::NotAromatic) => None,
        Some(AromaticValenceForm::Undetermined) | None => match view.aromatic_valence() {
            NumForm::Lit(valence) => u8::try_from(valence).ok(),
            _ => None,
        },
    }
}

/// Connected components of the aromatic-candidate graph: the perception's
/// rings whose members are all aromatic-capable, connected over shared
/// atoms. The aromaticity rule's acceptance couples atoms only within a
/// component, so enumeration, validity, selection, and the assignment bound
/// are all per component.
fn candidate_components<F>(rings: &RingSet, capable: F) -> Vec<BTreeSet<AtomId>>
where
    F: Fn(AtomId) -> bool,
{
    let mut components: Vec<BTreeSet<AtomId>> = Vec::new();
    for ring in rings.iter() {
        if !ring.atoms().iter().all(|&atom| capable(atom)) {
            continue;
        }
        let ring: BTreeSet<AtomId> = ring.atoms().iter().copied().collect();
        let (connected, disjoint): (Vec<_>, Vec<_>) = components
            .into_iter()
            .partition(|component| !component.is_disjoint(&ring));
        let mut merged = ring;
        for component in connected {
            merged.extend(component);
        }
        components = disjoint;
        components.push(merged);
    }
    components.sort();
    components
}

fn compare_restrictions(
    a: &[(AtomId, AtomForm)],
    b: &[(AtomId, AtomForm)],
    tie_break: ValenceTieBreak,
) -> Ordering {
    let key = tie_break.key();
    let b_forms: BTreeMap<AtomId, &AtomForm> = b.iter().map(|(atom, form)| (*atom, form)).collect();
    let mut a_sorted: Vec<&(AtomId, AtomForm)> = a.iter().collect();
    a_sorted.sort_unstable_by_key(|(atom, _)| *atom);
    for (atom, a_form) in a_sorted {
        let Some(b_form) = b_forms.get(atom) else {
            continue;
        };
        let ordering = compare_by_key(key, a_form, b_form);
        if ordering != Ordering::Equal {
            return ordering;
        }
    }
    Ordering::Equal
}

#[cfg(test)]
mod tests {
    use rstest::{fixture, rstest};
    use umol_graph_core::{
        ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
        RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
    };
    use umol_graph_ir::ir::{
        AromaticSystemId, BondConstraintKey, BondId, Edit, Edits, NumForm, RingConfig, RingModel,
        RingSetKind, UnpairedElectronsForm,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_concrete};

    use super::*;
    use crate::ops::model::{AromaticityRule, AromaticityTieBreak, ElementScope, RingLimits};
    use crate::ops::valence::AtomCompletions;

    #[rustfmt::skip]
    #[rstest]
    #[case::fused_pair(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]}"#),
        vec![(0..10).map(AtomId).collect::<BTreeSet<_>>()]
    )]
    #[case::coupled_rings(
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [5 6 "1"]
                    [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 10 "1"] [10 11 "1"] [11 6 "1"]]}"#),
        vec![
            (0..6).map(AtomId).collect::<BTreeSet<_>>(),
            (6..12).map(AtomId).collect::<BTreeSet<_>>(),
        ]
    )]
    #[case::chain(
        mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
        vec![]
    )]
    fn test_candidate_components(
        #[case] molecule: Molecule,
        #[case] expected: Vec<BTreeSet<AtomId>>,
    ) {
        let rings = molecule
            .rings(
                RingModel {
                    kind: RingSetKind::Relevant,
                    max_ring_size: 22,
                },
                RingConfig::default(),
            )
            .into_ring_set();
        assert_eq!(candidate_components(&rings, |_| true), expected);
    }

    #[rstest]
    fn test_candidate_components_capability() {
        // One incapable atom removes its rings; the remaining candidate ring
        // is its own component.
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [4 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 5 "1"]]}"#
        );
        let rings = molecule
            .rings(
                RingModel {
                    kind: RingSetKind::Relevant,
                    max_ring_size: 22,
                },
                RingConfig::default(),
            )
            .into_ring_set();
        assert_eq!(
            candidate_components(&rings, |atom| atom != AtomId(0)),
            vec![[4, 5, 6, 7, 8, 9]
                .map(AtomId)
                .into_iter()
                .collect::<BTreeSet<_>>()]
        );
    }

    #[fixture]
    fn aromaticity_model() -> AromaticityModel {
        AromaticityModel {
            scope: ElementScope::Any,
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
            tie_break: AromaticityTieBreak::Strict,
        }
    }

    #[fixture]
    fn benzene() -> Molecule {
        mol_dsl!(
            r#"{
            :atoms ["C#i=#c0#h#n0#u0#s#v2#a" "C#i=#c0#h#n0#u0#s#v2#a"
                    "C#i=#c0#h#n0#u0#s#v2#a" "C#i=#c0#h#n0#u0#s#v2#a"
                    "C#i=#c0#h#n0#u0#s#v2#a" "C#i=#c0#h#n0#u0#s#v2#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#
        )
    }

    #[fixture]
    fn aromatic_valence_mismatch() -> Molecule {
        mol_dsl!(
            r#"{
            :atoms ["C#a2" "C#a0" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        )
    }

    #[fixture]
    fn aromatic_bond_constraint_mismatch() -> Molecule {
        mol_dsl!(
            r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a!"] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        )
    }

    #[rstest]
    fn test_aromaticity_resolve_config_default() {
        assert_eq!(
            AromaticityResolveConfig::default(),
            AromaticityResolveConfig {
                perception: AromaticityConfig::default(),
                aromatic_valence_failure: AromaticityFailurePolicy::Error,
                aromatic_system_failure: AromaticityFailurePolicy::Error,
                aromatic_valence_mismatch: AromaticityMismatchPolicy::Error,
                aromatic_bond_constraint_mismatch: AromaticBondConstraintMismatchPolicy::Error,
                reset_aromatic_valence: false,
            }
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_plan(aromaticity_model: AromaticityModel, benzene: Molecule) {
        assert_eq!(
            AromaticityResolver::with_config(
                &aromaticity_model,
                AromaticityResolveConfig {
                    perception: AromaticityConfig {
                        ring_config: RingConfig {
                            simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                            relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                        },
                        connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                        maximum_independent_set_algorithm:
                            MaximumIndependentSetAlgorithm::BranchAndBound,
                    },
                    reset_aromatic_valence: false,
                    ..AromaticityResolveConfig::default()
                },
            )
            .plan(&benzene),
            Ok(Solution::Determined(Edits::from_iter([
                Edit::AddAromaticSystem {
                    atoms: (0..6).map(|id| AtomHandle::Id(AtomId(id))).collect(),
                    attributes: AromaticSystemForm::from_electrons(vec![1; 6])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(0)),
                    old: None,
                    new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(1)),
                    old: None,
                    new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(2)),
                    old: None,
                    new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(3)),
                    old: None,
                    new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(4)),
                    old: None,
                    new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                },
                Edit::ModifyBondConstraint {
                    id: BondHandle::Id(BondId(5)),
                    old: None,
                    new: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
                },
            ])))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_plan_partial(aromaticity_model: AromaticityModel) {
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#a+" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
        }"#
        );

        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).plan(&molecule),
            Ok(Solution::Underdetermined(Edits::new()))
        );
    }

    #[rstest]
    #[case::error(
        AromaticityMismatchPolicy::Error,
        Solution::Contradictory(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceMismatch {
                atom: AtomId(0),
                system: AromaticSystemId(0),
            }
        ))
    )]
    #[case::keep(AromaticityMismatchPolicy::Keep, Solution::Determined(Edits::new()))]
    #[case::remove_constraint(
        AromaticityMismatchPolicy::RemoveConstraint,
        Solution::Determined(Edits::from_iter([
            Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(0)),
                old: Some(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(NumForm::Lit(2)),
                )),
                new: None,
            },
            Edit::ModifyAtomConstraint {
                id: AtomHandle::Id(AtomId(1)),
                old: Some(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(NumForm::Lit(0)),
                )),
                new: None,
            },
        ]))
    )]
    #[case::replace_entity(
        AromaticityMismatchPolicy::ReplaceEntity,
        Solution::Determined(Edits::from_iter([
            Edit::RemoveAromaticSystems {
                removes: vec![(
                    AromaticSystemHandle::Id(AromaticSystemId(0)),
                    (0..6).map(|id| AtomHandle::Id(AtomId(id))).collect(),
                    AromaticSystemForm::from_electrons(vec![1; 6]),
                )],
            },
            Edit::AddAromaticSystem {
                atoms: (0..6).map(|id| AtomHandle::Id(AtomId(id))).collect(),
                attributes: AromaticSystemForm::from_electrons(vec![2, 0, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            },
        ]))
    )]
    fn test_aromaticity_resolver_plan_aromatic_valence_mismatch(
        aromaticity_model: AromaticityModel,
        aromatic_valence_mismatch: Molecule,
        #[case] policy: AromaticityMismatchPolicy,
        #[case] expected: Solution<Edits, AromaticityContradiction>,
    ) {
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_valence_mismatch: policy,
                ..AromaticityResolveConfig::default()
            },
        );

        assert_eq!(resolver.plan(&aromatic_valence_mismatch), Ok(expected));
    }

    #[rstest]
    fn test_aromaticity_resolver_resolve_aromatic_valence_mismatch_reset(
        aromaticity_model: AromaticityModel,
        mut aromatic_valence_mismatch: Molecule,
    ) {
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_valence_mismatch: AromaticityMismatchPolicy::ReplaceEntity,
                reset_aromatic_valence: true,
                ..AromaticityResolveConfig::default()
            },
        );
        let expected = mol_dsl!(
            r#"{
            :atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{
                :atoms [0 1 2 3 4 5]
                :attrs "[2,0,1,1,1,1]#c0#u0#s"
            }]
        }"#
        );

        assert_eq!(
            resolver.resolve(&mut aromatic_valence_mismatch),
            Ok(Solution::Determined(()))
        );
        assert_eq!(aromatic_valence_mismatch, expected);
    }

    #[rstest]
    #[case::error(
        AromaticBondConstraintMismatchPolicy::Error,
        Solution::Contradictory(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticBondConstraintMismatch {
                bond: BondId(0),
                system: AromaticSystemId(0),
            }
        ))
    )]
    #[case::keep(
        AromaticBondConstraintMismatchPolicy::Keep,
        Solution::Determined(Edits::new())
    )]
    #[case::remove_constraint(
        AromaticBondConstraintMismatchPolicy::RemoveConstraint,
        Solution::Determined(Edits::from_iter([Edit::ModifyBondConstraint {
            id: BondHandle::Id(BondId(0)),
            old: Some(BondConstraintForm::Aromatic(BooleanForm::Lit(false))),
            new: None,
        }]))
    )]
    fn test_aromaticity_resolver_plan_aromatic_bond_constraint_mismatch(
        aromaticity_model: AromaticityModel,
        aromatic_bond_constraint_mismatch: Molecule,
        #[case] policy: AromaticBondConstraintMismatchPolicy,
        #[case] expected: Solution<Edits, AromaticityContradiction>,
    ) {
        let resolver = AromaticityResolver::with_config(
            &aromaticity_model,
            AromaticityResolveConfig {
                aromatic_bond_constraint_mismatch: policy,
                ..AromaticityResolveConfig::default()
            },
        );

        assert_eq!(
            resolver.plan(&aromatic_bond_constraint_mismatch),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::conformant(
        AromaticityModel::daylight(),
        AromaticityResolveConfig::default(),
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 5 :aromatic] [5 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#)
    )]
    #[case::rejected_projections(
        AromaticityModel::mdl(),
        AromaticityResolveConfig {
            aromatic_valence_failure: AromaticityFailurePolicy::Keep,
            ..AromaticityResolveConfig::default()
        },
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
        }"#)
    )]
    #[case::rejected_existing_system(
        AromaticityModel::mdl(),
        AromaticityResolveConfig {
            aromatic_valence_failure: AromaticityFailurePolicy::Keep,
            aromatic_system_failure: AromaticityFailurePolicy::Keep,
            ..AromaticityResolveConfig::default()
        },
        mol_dsl!(r#"{
            :atoms ["O#n1#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 :aromatic] [1 2 :aromatic] [2 3 :aromatic]
                    [3 4 :aromatic] [4 0 :aromatic]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]
        }"#)
    )]
    fn test_aromaticity_resolver_plan_identity(
        #[case] model: AromaticityModel,
        #[case] config: AromaticityResolveConfig,
        #[case] molecule: Molecule,
    ) {
        assert_eq!(
            AromaticityResolver::with_config(&model, config).plan(&molecule),
            Ok(Solution::Determined(Edits::new()))
        );
    }

    #[rstest]
    #[case::homogeneous_localized(
        AromaticityResolveConfig::default(),
        mol_dsl_concrete!(r#"{:atoms ["C #h #a" "C #h #a" "C #c+ #h #a0"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"]]}"#),
        NumForm::Lit(0),
        vec![NumForm::Lit(0), NumForm::Lit(0), NumForm::Lit(1)],
        vec![
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(0))),
        ]
    )]
    #[case::heterogeneous_localized(
        AromaticityResolveConfig::default(),
        mol_dsl_concrete!(r#"{:atoms ["N #c+ #h #a" "C #h #a" "C #h #a"
                                      "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]
                                      [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        NumForm::Lit(0),
        vec![
            NumForm::Lit(1), NumForm::Lit(0), NumForm::Lit(0),
            NumForm::Lit(0), NumForm::Lit(0), NumForm::Lit(0),
        ],
        vec![Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))); 6]
    )]
    #[case::accepted_system_with_rejected_projections(
        AromaticityResolveConfig {
            aromatic_valence_failure: AromaticityFailurePolicy::Keep,
            ..AromaticityResolveConfig::default()
        },
        mol_dsl_concrete!(r#"{
            :atoms ["C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a" "C#h#a"
                    "C#h3#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#),
        NumForm::Lit(0),
        vec![NumForm::Lit(0); 7],
        vec![
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
            Some(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
        ]
    )]
    #[case::reset_source_constraints(
        AromaticityResolveConfig {
            reset_aromatic_valence: true,
            ..AromaticityResolveConfig::default()
        },
        benzene(),
        NumForm::Lit(0),
        vec![NumForm::Lit(0); 6],
        vec![None; 6]
    )]
    fn test_aromaticity_resolver_resolve(
        aromaticity_model: AromaticityModel,
        #[case] config: AromaticityResolveConfig,
        #[case] mut molecule: Molecule,
        #[case] expected_system_charge: NumForm,
        #[case] expected_atom_charges: Vec<NumForm>,
        #[case] expected_aromatic_valences: Vec<Option<AromaticValenceForm>>,
    ) {
        assert_eq!(
            AromaticityResolver::with_config(&aromaticity_model, config).resolve(&mut molecule),
            Ok(Solution::Determined(()))
        );
        assert_eq!(molecule.aromatic_systems().count(), 1);
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .attributes
                .charge,
            expected_system_charge
        );
        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.attributes.charge.clone())
                .collect::<Vec<_>>(),
            expected_atom_charges
        );
        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.attributes.constraints.aromatic_valence().cloned())
                .collect::<Vec<_>>(),
            expected_aromatic_valences
        );
        assert!(molecule.bonds().iter().all(|bond| matches!(
            bond.attributes.constraints.get(BondConstraintKey::Aromatic),
            Some(BondConstraintForm::Aromatic(BooleanForm::Lit(true)))
        )));
    }

    type SelectOutcome = Solution<ResolveState, AromaticityContradiction>;

    #[rstest]
    #[case::unique_survivor(
        mol_dsl!(r#"{:atoms ["N#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                     :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![
                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
            ]),
            (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
        ]),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState { completions: AtomCompletions::from_iter([
                (AtomId(0), smallvec![atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2")]),
                (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            ]), systems: vec![(
                (0..5).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )], tie_breaks: Vec::new() })
    )]
    #[case::quinoline(
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "N#c0" "C#c0"]
                     :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                             [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]]}"#),
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(6), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(7), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(8), smallvec![
                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
            ]),
            (AtomId(9), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
        ]),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState { completions: AtomCompletions::from_iter([
                (AtomId(0), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(6), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(7), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(8), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                (AtomId(9), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            ]), systems: vec![(
                (0..10).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![1; 10])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )], tie_breaks: Vec::new() })
    )]
    #[case::tie_most_saturated(
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                     :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![
                atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
                atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
            ]),
            (AtomId(1), smallvec![
                atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
                atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
            ]),
            (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
        ]),
        ValenceTieBreak::MostSaturated,
        Solution::Determined(ResolveState { completions: AtomCompletions::from_iter([
                (AtomId(0), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0")]),
                (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2")]),
                (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            ]), systems: vec![(
                (0..6).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![0, 2, 1, 1, 1, 1])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )], tie_breaks: vec![AtomId(0), AtomId(1)] })
    )]
    #[case::stored_assertions_empty_carrier(
        benzene(),
        AtomCompletions::new(),
        ValenceTieBreak::Strict,
        Solution::Determined(ResolveState { completions: AtomCompletions::new(), systems: vec![(
                (0..6).map(AtomId).collect(),
                AromaticSystemForm::from_electrons(vec![1; 6])
                    .with_charge(0)
                    .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
            )], tie_breaks: Vec::new() })
    )]
    #[case::unclaimed_aromatic_contradiction(
        mol_dsl!(r#"{:atoms ["N#c0"] :bonds []}"#),
        AtomCompletions::from_iter([(AtomId(0), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")])]),
        ValenceTieBreak::Strict,
        Solution::Contradictory(AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }
        ))
    )]
    #[case::stored_undetermined_gate(
        mol_dsl!(r#"{:atoms ["C#a+" "C#a"] :bonds [[0 1 "1"]]}"#),
        AtomCompletions::new(),
        ValenceTieBreak::Strict,
        Solution::Underdetermined(ResolveState::default())
    )]
    fn test_aromaticity_resolver_select(
        aromaticity_model: AromaticityModel,
        #[case] molecule: Molecule,
        #[case] completions: AtomCompletions,
        #[case] tie_break: ValenceTieBreak,
        #[case] expected: SelectOutcome,
    ) {
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).select(
                &molecule,
                ResolveState {
                    completions,
                    ..ResolveState::default()
                },
                tie_break,
            ),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::tie_strict(
        mol_dsl!(r#"{:atoms ["C#c0" "C#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                     :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
        AtomCompletions::from_iter([
            (AtomId(0), smallvec![
                atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
                atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
            ]),
            (AtomId(1), smallvec![
                atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a0"),
                atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
            ]),
            (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
        ]),
        ValenceTieBreak::Strict
    )]
    fn test_aromaticity_resolver_select_identity(
        aromaticity_model: AromaticityModel,
        #[case] molecule: Molecule,
        #[case] completions: AtomCompletions,
        #[case] tie_break: ValenceTieBreak,
    ) {
        let state = ResolveState {
            completions,
            ..ResolveState::default()
        };
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).select(
                &molecule,
                state.clone(),
                tie_break
            ),
            Ok(Solution::Determined(state))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_min_electron_count() {
        // Two totals pass the rule on the same members (all-pyridinic 6,
        // all-pyrrolic 10); the electron component picks the smaller and
        // records the members, with no value key consulted.
        let model = AromaticityModel {
            scope: ElementScope::Any,
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
            tie_break: AromaticityTieBreak::MinElectronCount,
        };
        let molecule = mol_dsl!(
            r#"{:atoms ["N#c0" "N#c0" "C#c0" "N#c0" "N#c0" "C#c0"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#
        );
        assert_eq!(
            AromaticityResolver::new(&model).select(
                &molecule,
                ResolveState {
                    completions: AtomCompletions::from_iter([
                        (
                            AtomId(0),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (
                            AtomId(1),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                        (
                            AtomId(3),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (
                            AtomId(4),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    ]),
                    ..ResolveState::default()
                },
                ValenceTieBreak::Strict,
            ),
            Ok(Solution::Determined(ResolveState {
                completions: AtomCompletions::from_iter([
                    (AtomId(0), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                    (AtomId(1), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                    (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    (AtomId(3), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                    (AtomId(4), smallvec![atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a")]),
                    (AtomId(5), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                ]),
                systems: vec![(
                    (0..6).map(AtomId).collect(),
                    AromaticSystemForm::from_electrons(vec![1, 1, 1, 1, 1, 1])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                )],
                tie_breaks: (0..6).map(AtomId).collect(),
            }))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_tolerated_carrier() {
        // Keep policy admits the no-system assignment alongside the full
        // ring; `MinElectronCount` realizes the full ring and records its members.
        let model = AromaticityModel {
            scope: ElementScope::Any,
            rule: AromaticityRule::Hueckel {
                ring_limits: RingLimits::default(),
            },
            tie_break: AromaticityTieBreak::MinElectronCount,
        };
        let molecule = mol_dsl!(
            r#"{:atoms ["N#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#
        );
        assert_eq!(
            AromaticityResolver::with_config(
                &model,
                AromaticityResolveConfig {
                    aromatic_valence_failure: AromaticityFailurePolicy::Keep,
                    ..AromaticityResolveConfig::default()
                },
            )
            .select(
                &molecule,
                ResolveState {
                    completions: AtomCompletions::from_iter([
                        (
                            AtomId(0),
                            smallvec![
                                atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                                atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                            ]
                        ),
                        (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                        (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                        (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                        (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    ]),
                    ..ResolveState::default()
                },
                ValenceTieBreak::Strict,
            ),
            Ok(Solution::Determined(ResolveState {
                completions: AtomCompletions::from_iter([
                    (AtomId(0), smallvec![atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2")]),
                    (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                    (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                ]),
                systems: vec![(
                    (0..5).map(AtomId).collect(),
                    AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1])
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                )],
                tie_breaks: (0..5).map(AtomId).collect(),
            }))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_tolerated_carrier_identity(
        aromaticity_model: AromaticityModel,
    ) {
        // The same survivors under the model's `Strict`: structurally
        // distinct, so the members stay plural and the state passes through.
        let molecule = mol_dsl!(
            r#"{:atoms ["N#c0" "C#c0" "C#c0" "C#c0" "C#c0"]
                :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#
        );
        let state = ResolveState {
            completions: AtomCompletions::from_iter([
                (
                    AtomId(0),
                    smallvec![
                        atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                        atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                    ],
                ),
                (AtomId(1), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(2), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(3), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
                (AtomId(4), smallvec![atom_dsl!("C#i=#c0#h#n0#u0#s#v2#a")]),
            ]),
            ..ResolveState::default()
        };
        assert_eq!(
            AromaticityResolver::with_config(
                &aromaticity_model,
                AromaticityResolveConfig {
                    aromatic_valence_failure: AromaticityFailurePolicy::Keep,
                    ..AromaticityResolveConfig::default()
                },
            )
            .select(&molecule, state.clone(), ValenceTieBreak::Strict),
            Ok(Solution::Determined(state))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_stored_conflict(aromaticity_model: AromaticityModel) {
        // No completion of the flexible atom reproduces the stored system:
        // `#a2` breaks the count, `#a!` removes the candidate ring.
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#c0" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        );
        assert_eq!(
            AromaticityResolver::new(&aromaticity_model).select(
                &molecule,
                ResolveState {
                    completions: AtomCompletions::from_iter([(
                        AtomId(0),
                        smallvec![
                            atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
                            atom_dsl!("C#i=#c0#h2#n0#u0#s#v2#a!"),
                        ]
                    )]),
                    ..ResolveState::default()
                },
                ValenceTieBreak::Strict,
            ),
            Ok(Solution::Contradictory(
                AromaticityContradiction::Inconsistency(
                    AromaticityInconsistency::AromaticSystemFailure {
                        system: AromaticSystemId(0)
                    }
                )
            ))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_select_stored_conflict_identity(
        aromaticity_model: AromaticityModel,
    ) {
        // Same conflict under `Keep`: the component is inert and the state
        // passes through unchanged.
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#c0" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]
        }"#
        );
        let state = ResolveState {
            completions: AtomCompletions::from_iter([(
                AtomId(0),
                smallvec![
                    atom_dsl!("C#i=#c0#h0#n0#u0#s#v2#a2"),
                    atom_dsl!("C#i=#c0#h2#n0#u0#s#v2#a!"),
                ],
            )]),
            ..ResolveState::default()
        };
        assert_eq!(
            AromaticityResolver::with_config(
                &aromaticity_model,
                AromaticityResolveConfig {
                    aromatic_system_failure: AromaticityFailurePolicy::Keep,
                    ..AromaticityResolveConfig::default()
                },
            )
            .select(&molecule, state.clone(), ValenceTieBreak::Strict),
            Ok(Solution::Determined(state))
        );
    }

    #[rstest]
    fn test_aromaticity_resolver_resolve_identity(
        aromaticity_model: AromaticityModel,
        mut benzene: Molecule,
    ) {
        let resolver = AromaticityResolver::new(&aromaticity_model);
        assert_eq!(resolver.resolve(&mut benzene), Ok(Solution::Determined(())));
        let expected = benzene.clone();

        assert_eq!(resolver.resolve(&mut benzene), Ok(Solution::Determined(())));
        assert_eq!(benzene, expected);
    }

    #[rstest]
    #[case::clar_heterocycle(
        AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar, tie_break: AromaticityTieBreak::Strict },
        mol_dsl_concrete!(r#"{:atoms ["N #h #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        AromaticityContradiction::ClarNonBenzenoid(
            "Clar model requires benzenoid input but non-carbon aromatic atoms are present".to_string()
        )
    )]
    #[case::aromatic_valence_failure(
        AromaticityModel::mdl(),
        mol_dsl_concrete!(r#"{:atoms ["O #n1 #a2" "C #h #a" "C #h #a" "C #h #a" "C #h #a"]
                              :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]}"#),
        AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) }
        )
    )]
    #[case::aromatic_system_failure(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]"}]
        }"#),
        AromaticityContradiction::Inconsistency(
            AromaticityInconsistency::AromaticSystemFailure {
                system: AromaticSystemId(0)
            }
        )
    )]
    fn test_aromaticity_resolver_resolve_contradiction(
        #[case] model: AromaticityModel,
        #[case] mut molecule: Molecule,
        #[case] expected: AromaticityContradiction,
    ) {
        let original = molecule.clone();
        assert_eq!(
            AromaticityResolver::new(&model).resolve(&mut molecule),
            Ok(Solution::Contradictory(expected))
        );
        assert_eq!(molecule, original);
    }
}
