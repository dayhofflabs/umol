//! Forward atom completion and search-independence of aromatic completion selection.
//!
//! Completion checks resolve against independent atom states for carbon/radical/isotope
//! chains and mixed-element chains with charges, lone pairs, and localized bonds. The
//! input leaves electron fields undetermined; both valence sources and tie-breaks must
//! produce the expected complete molecule.
//!
//! Aromatic discharge is checked against independently assembled results for
//! separate or linked benzene rings with mixed stored and newly added systems.
//! Stored-system assertions and unrelated keys survive; newly realized assertions
//! and negative linking-bond assertions are removed. Resolution is idempotent.
//!
//! Added aromatic systems agree with the input atom evidence and contain no bond
//! asserting non-aromaticity. The skeleton domain retains open, contradictory,
//! and stored-system cases. A second rendering supplies complete carbon atoms
//! with one-electron contributions while preserving topology and bond assertions,
//! so successful additions are exercised independently of atom completion.
//!
//! Stereo discharge is checked against independent results for mixtures of stored
//! and newly added tetrahedral and cis-trans entities with literal or open configurations.
//! Existing assertions and unrelated keys are preserved, and resolution is idempotent.
//!
//! Skeleton properties compare each determined standalone or composite result's exact
//! remaining assertions with those realized by its additions and completions. SMILES and
//! MOL renderings that resolve under Error policies leave no aromatic or stereo assertion.
//! Projection followed by resolution reproduces the systems modulo participant frames
//! and removes the assertions written by projection under the same chemistry model.
//!
//! Semantic property: for carriers whose components stay within the
//! enumeration bound, both aromatic selectors return the outcome of a
//! selection that enumerates every completion combination of every component
//! exhaustively — the pruned search never removes a valid completion combination.
//!
//! Operational domain: `strategies::select_scenario` — one- and two-ring
//! Hückel skeletons, all ring atoms in the carrier with literal
//! contributions, no stored systems, one positive localized-bond aromaticity
//! constraint, both failure policies and both tie-breaks on each axis.
//!
//! Validation method: comparison with the definition-level selection below,
//! which enumerates completion combinations flat with no pruning. Perception
//! (`find_systems`) and the value-key comparison (`compare_by_key`) are
//! shared with production — the property targets the search and selection
//! stages, not the perception rule. The explicit report method also matches the
//! reference tie-break list; ordinary selection leaves that list empty.

use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet};

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph::ops::aromaticity::{
    AromaticityConfig, AromaticityContradiction, AromaticityInconsistency, AromaticityPerceiver,
};
use umol_graph::ops::model::{
    AromaticityModel, AromaticityRule, AromaticityTieBreak, ChemistryModel, ElementScope,
    RingLimits, StereoModel, ValenceModel, ValenceTieBreak,
};
use umol_graph::ops::resolve::{
    AromaticityFailurePolicy, AromaticityResolveConfig, AromaticityResolver, IsotopePolicy,
    ProjectFlags, ResolveConfig, ResolveState, Resolver, StereoResolver,
};
use umol_graph::ops::valence::compare::compare_by_key;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticValence, AromaticValenceForm, AtomConstraintForm,
    AtomConstraintKey, AtomForm, AtomId, BondConstraintForm, BondConstraintKey, BondForm, BondId,
    BooleanForm, CisTransStereoForm, Constraint, ElectronCountsForm, ElementForm, IsotopeMassForm,
    Lattice, Molecule, MoleculeConstraint, MoleculeEntries, NumForm, Reframe, RingConfig,
    RingModel, RingSetKind, StereoAtomForm, StereoBondForm, StereoCoset, StereoKind, StereoLigand,
    StereoLigandKind, TetrahedralStereoForm, TryIntoIr, UnpairedElectronsForm,
};
use umol_graph_ir::{atom_dsl, bond_dsl};
use umol_io::ctfile::parser::parse_mol_to_ir;
use umol_io::smiles::Smiles;
use umol_utils::solution::Solution;

use crate::strategies::{aromatic_skeleton, select_scenario, SelectScenario};

fn contribution(form: &AtomForm) -> Option<u8> {
    match form.constraints.aromatic_valence() {
        AromaticValenceForm::Aromatic(NumForm::Lit(valence)) => u8::try_from(*valence).ok(),
        _ => None,
    }
}

fn compare_restrictions(
    a: &[(AtomId, AtomForm)],
    b: &[(AtomId, AtomForm)],
    tie_break: ValenceTieBreak,
) -> Ordering {
    let key = tie_break.key();
    let b_forms: BTreeMap<AtomId, &AtomForm> = b.iter().map(|(atom, form)| (*atom, form)).collect();
    for (atom, a_form) in a {
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

/// Definition-level selection: enumerates completion combinations per component
/// without pruning or an enumeration bound, then applies the documented validity
/// and selection stages. Scoped to the generated domain: no stored systems,
/// no atom-level assertions, every ring atom in the carrier.
fn exhaustive_select(
    scenario: &SelectScenario,
) -> Solution<ResolveState, AromaticityContradiction> {
    let SelectScenario {
        molecule,
        completions,
        model,
        config,
        tie_break,
    } = scenario;
    let perception = AromaticityPerceiver::new(model);
    let mut completions = completions.clone();
    let mut systems: Vec<(Vec<AtomId>, AromaticSystemForm)> = Vec::new();
    let mut tie_breaks: Vec<AtomId> = Vec::new();

    let flexible: Vec<(AtomId, Vec<Option<u8>>)> = completions
        .iter()
        .filter_map(|(atom, disjuncts)| {
            let contributions: Vec<Option<u8>> = disjuncts.iter().map(contribution).collect();
            (contributions.iter().collect::<BTreeSet<_>>().len() > 1)
                .then_some((atom, contributions))
        })
        .collect();

    let rings = molecule
        .rings(
            RingModel {
                kind: RingSetKind::Relevant,
                max_ring_size: 22,
            },
            RingConfig::default(),
        )
        .into_ring_set();
    let capable = |atom: AtomId| {
        completions
            .get(atom)
            .is_some_and(|disjuncts| disjuncts.iter().any(|form| contribution(form).is_some()))
    };
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

    type CompletionCombination = (
        Vec<(AtomId, AtomForm)>,
        Vec<(Vec<AtomId>, AromaticSystemForm)>,
    );
    let mut accepted: Vec<(Vec<AtomId>, AromaticSystemForm)> = Vec::new();
    let mut tie_break_uses: BTreeSet<AtomId> = BTreeSet::new();
    let mut claimed: BTreeSet<AtomId> = BTreeSet::new();
    let required_bonds: Vec<_> = molecule
        .bonds()
        .iter()
        .filter(|bond| {
            matches!(
                bond.attributes().constraints.aromatic(),
                BooleanForm::Lit(true)
            )
        })
        .map(|bond| (bond.id(), bond.atom_ids()))
        .collect();
    let mut claimed_bonds = BTreeSet::new();
    for component in &components {
        let component_flexible: Vec<(AtomId, Vec<Option<u8>>)> = flexible
            .iter()
            .filter(|(atom, _)| component.contains(atom))
            .cloned()
            .collect();
        let aromatic_only: Vec<AtomId> = component
            .iter()
            .copied()
            .filter(|&atom| {
                completions.get(atom).is_some_and(|disjuncts| {
                    disjuncts.iter().all(|form| {
                        matches!(
                            form.constraints.aromatic_valence(),
                            AromaticValenceForm::Aromatic(_)
                        )
                    })
                })
            })
            .collect();

        let mut combinations: Vec<CompletionCombination> = Vec::new();
        let mut completion_indices = vec![0usize; component_flexible.len()];
        loop {
            let choice: BTreeMap<AtomId, usize> = component_flexible
                .iter()
                .zip(&completion_indices)
                .map(|(&(atom, _), &index)| (atom, index))
                .collect();
            let outcome = perception
                .find_systems(molecule, AromaticityConfig::default(), |atom| {
                    match choice.get(&atom) {
                        Some(&index) => {
                            contribution(&completions.get(atom).expect("flexible atom")[index])
                        }
                        None => completions
                            .get(atom)
                            .and_then(|disjuncts| contribution(&disjuncts[0])),
                    }
                    .map(|value| AromaticValence::Aromatic(i64::from(value)))
                    .unwrap_or(AromaticValence::NotAromatic)
                })
                .expect("Hückel perception is infallible");
            let found = match outcome {
                Solution::Determined(found) => found,
                Solution::Underdetermined(_) => {
                    return Solution::Underdetermined(ResolveState {
                        completions,
                        systems,
                        tie_breaks,
                    });
                }
                Solution::Contradictory(contradiction) => {
                    return Solution::Contradictory(contradiction);
                }
            };
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
            let combination = (restriction, partition);
            if !combinations.contains(&combination) {
                combinations.push(combination);
            }

            let mut position = component_flexible.len();
            loop {
                if position == 0 {
                    break;
                }
                position -= 1;
                completion_indices[position] += 1;
                if completion_indices[position] < component_flexible[position].1.len() {
                    break;
                }
                completion_indices[position] = 0;
            }
            if component_flexible.is_empty() || completion_indices.iter().all(|&index| index == 0) {
                break;
            }
        }

        let mut valid = combinations;
        if config.aromatic_valence_failure == AromaticityFailurePolicy::Error {
            valid.retain(|(_, partition)| {
                aromatic_only
                    .iter()
                    .all(|atom| partition.iter().any(|(atoms, _)| atoms.contains(atom)))
            });
        }
        if valid.is_empty() {
            continue;
        }
        for (_, partition) in &valid {
            for (atoms, _) in partition {
                claimed.extend(atoms.iter().copied());
            }
        }
        if config.aromatic_valence_failure == AromaticityFailurePolicy::Error {
            valid.retain(|(_, partition)| {
                required_bonds.iter().all(|(_, endpoints)| {
                    !endpoints.iter().all(|atom| component.contains(atom))
                        || partition
                            .iter()
                            .any(|(members, _)| endpoints.iter().all(|atom| members.contains(atom)))
                })
            });
            if valid.is_empty() {
                continue;
            }
            for &(bond, endpoints) in &required_bonds {
                if valid.iter().any(|(_, partition)| {
                    partition
                        .iter()
                        .any(|(members, _)| endpoints.iter().all(|atom| members.contains(atom)))
                }) {
                    claimed_bonds.insert(bond);
                }
            }
        }

        let mut structural_decided = false;
        if model.tie_break == AromaticityTieBreak::MinElectronCount {
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
            let domain: BTreeSet<AtomId> = restrictions[0].iter().map(|(atom, _)| *atom).collect();
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
                .max_by(|a, b| compare_restrictions(a, b, *tie_break))
                .expect("non-empty restrictions");
            let unique = restrictions
                .iter()
                .filter(|restriction| compare_restrictions(restriction, best, *tie_break).is_eq())
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
        let narrowed: Vec<(AtomId, AtomForm)> = winner_restriction.clone();
        for (atom, form) in &narrowed {
            if by_key && completions.get(*atom).is_some_and(|entry| entry.len() > 1) {
                tie_break_uses.insert(*atom);
            }
            completions.insert(*atom, smallvec::smallvec![form.clone()]);
        }
    }

    if config.aromatic_valence_failure == AromaticityFailurePolicy::Error {
        for (atom, disjuncts) in completions.iter() {
            if claimed.contains(&atom) {
                continue;
            }
            if disjuncts.iter().all(|form| {
                matches!(
                    form.constraints.aromatic_valence(),
                    AromaticValenceForm::Aromatic(_)
                )
            }) {
                return Solution::Contradictory(
                    AromaticityInconsistency::AromaticValenceFailure { atom }.into(),
                );
            }
        }
        for (bond, _) in required_bonds {
            if !claimed_bonds.contains(&bond) {
                return Solution::Contradictory(
                    AromaticityInconsistency::AromaticBondConstraintFailure { bond }.into(),
                );
            }
        }
    }

    systems.extend(accepted);
    tie_breaks.extend(tie_break_uses);
    tie_breaks.sort_unstable();
    tie_breaks.dedup();
    Solution::Determined(ResolveState {
        completions,
        systems,
        tie_breaks,
    })
}

proptest! {
    #[test]
    fn test_stereo_resolver_resolve_discharge(
        sites in prop::collection::vec((any::<bool>(), any::<bool>(), 0u32..3), 1..5),
    ) {
        let mut source = MoleculeEntries::default();
        let mut expected = MoleculeEntries::default();
        let mut added_atoms = Vec::new();
        let mut added_bonds = Vec::new();
        for (cis_trans, stored, coset) in sites {
            let coset = if coset == 2 { StereoCoset::Undetermined } else { StereoCoset::Lit(coset) };
            let offset = source.atoms.len();
            let a = AtomId::from(offset);
            let b = AtomId::from(offset + 1);
            let c = AtomId::from(offset + 2);
            let d = AtomId::from(offset + 3);
            let mut atoms = if cis_trans {
                vec![atom_dsl!("C#h"), atom_dsl!("C#h"), atom_dsl!("F"), atom_dsl!("Cl")]
            } else {
                vec![atom_dsl!("C#h"), atom_dsl!("F"), atom_dsl!("Cl"), atom_dsl!("Br")]
            };
            atoms[0].constraints.set(AtomConstraintForm::Valence(NumForm::Lit(3)));
            let mut bonds = if cis_trans {
                vec![(a, b, bond_dsl!("2#R!")), (a, c, bond_dsl!("1")), (b, d, bond_dsl!("1"))]
            } else {
                vec![(a, b, bond_dsl!("1#R!")), (a, c, bond_dsl!("1")), (a, d, bond_dsl!("1"))]
            };
            let mut expected_atoms = atoms.clone();
            let mut expected_bonds = bonds.clone();
            if cis_trans {
                let assertion = BondConstraintForm::CisTransStereo(
                    CisTransStereoForm::Stereo(coset.clone()));
                bonds[0].2.constraints.set(assertion.clone());
                if stored { expected_bonds[0].2.constraints.set(assertion); }
                let entity = (BondId::from(source.bonds.len()), vec![
                    StereoLigand::new(c, StereoLigandKind::Atom),
                    StereoLigand::new(a, StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(d, StereoLigandKind::Atom),
                    StereoLigand::new(b, StereoLigandKind::ImplicitHydrogen),
                ], StereoBondForm::new(StereoKind::CisTrans, coset.clone()));
                if stored {
                    source.stereo_bonds.push(entity.clone());
                    expected.stereo_bonds.push(entity);
                } else { added_bonds.push(entity); }
            } else {
                let assertion = AtomConstraintForm::TetrahedralStereo(
                    TetrahedralStereoForm::Stereo(coset.clone()));
                atoms[0].constraints.set(assertion.clone());
                if stored { expected_atoms[0].constraints.set(assertion); }
                let entity = (a, vec![
                    StereoLigand::new(b, StereoLigandKind::Atom),
                    StereoLigand::new(c, StereoLigandKind::Atom),
                    StereoLigand::new(d, StereoLigandKind::Atom),
                    StereoLigand::new(a, StereoLigandKind::ImplicitHydrogen),
                ], StereoAtomForm::new(StereoKind::Tetrahedral, coset.clone()));
                if stored {
                    source.stereo_atoms.push(entity.clone());
                    expected.stereo_atoms.push(entity);
                } else { added_atoms.push(entity); }
            }
            source.atoms.extend(atoms);
            source.bonds.extend(bonds);
            expected.atoms.extend(expected_atoms);
            expected.bonds.extend(expected_bonds);
        }
        expected.stereo_atoms.extend(added_atoms);
        expected.stereo_bonds.extend(added_bonds);
        let mut molecule = Molecule::from_entries(source);
        let expected = Molecule::from_entries(expected);
        let resolver = StereoResolver::new(&StereoModel::default());
        prop_assert_eq!(resolver.resolve(molecule.clone()), Ok(Solution::Determined(expected.clone())));
        prop_assert_eq!(resolver.resolve_into(&mut molecule), Ok(Solution::Determined(())));
        prop_assert_eq!(&molecule, &expected);
        prop_assert_eq!(resolver.resolve(molecule), Ok(Solution::Determined(expected)));
    }

    #[test]
    fn test_aromaticity_resolver_resolve_discharge(
        rings in prop::collection::vec((any::<bool>(), 0u8..3), 1..5),
        linked in any::<bool>(),
    ) {
        let mut source = MoleculeEntries::default();
        let mut expected = MoleculeEntries::default();
        let mut additions = Vec::new();
        for (ring, &(stored, bond_assertions)) in rings.iter().enumerate() {
            let offset = ring * 6;
            for position in 0..6 {
                let degree = 2 + i64::from(linked
                    && ((position == 0 && ring > 0)
                        || (position == 3 && ring + 1 < rings.len())));
                let mut atom = atom_dsl!("C#i=#c0#n0#u0#s");
                atom.implicit_hydrogens = NumForm::Lit(3 - degree);
                atom.constraints.set(AtomConstraintForm::Valence(NumForm::Lit(degree)));
                let mut asserted = atom.clone();
                asserted.constraints.set(AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(NumForm::Lit(1))));
                expected.atoms.push(if stored { asserted.clone() } else { atom });
                source.atoms.push(asserted);

                let a = AtomId::from(offset + position);
                let b = AtomId::from(offset + (position + 1) % 6);
                let bond = bond_dsl!("1#c0#u0#s#R+");
                let mut asserted = bond.clone();
                if bond_assertions != 0 {
                    asserted.constraints.set(BondConstraintForm::Aromatic(
                        if bond_assertions == 1 { BooleanForm::Undetermined }
                        else { BooleanForm::Lit(true) }
                    ));
                }
                expected.bonds.push((a, b, if stored { asserted.clone() } else { bond }));
                source.bonds.push((a, b, asserted));
            }
            let system = ((offset..offset + 6).map(AtomId::from).collect(), AromaticSystemForm {
                electrons: ElectronCountsForm::Lit(vec![1; 6]),
                charge: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                constraints: Default::default(),
            });
            if stored {
                source.aromatic.push(system.clone());
                expected.aromatic.push(system);
            } else {
                additions.push(system);
            }
            if linked && ring > 0 {
                let a = AtomId::from(offset - 3);
                let b = AtomId::from(offset);
                let bond = bond_dsl!("1#c0#u0#s");
                let mut asserted = bond.clone();
                asserted.constraints.set(BondConstraintForm::Aromatic(BooleanForm::Lit(false)));
                source.bonds.push((a, b, asserted));
                expected.bonds.push((a, b, bond));
            }
        }
        expected.aromatic.extend(additions);
        let constraint = Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None, sum: NumForm::Lit(1),
        });
        source.constraints.push(constraint.clone());
        expected.constraints.push(constraint);
        let mut molecule = Molecule::from_entries(source);
        let expected = Molecule::from_entries(expected);
        let resolver = AromaticityResolver::new(&AromaticityModel::daylight());
        let model = ChemistryModel::default();
        let composite = Resolver::new(&model);
        prop_assert_eq!(composite.resolve(molecule.clone()), Ok(Solution::Determined(expected.clone())));
        let mut borrowed = molecule.clone();
        prop_assert_eq!(composite.resolve_into(&mut borrowed), Ok(Solution::Determined(())));
        prop_assert_eq!(borrowed, expected.clone());
        prop_assert_eq!(resolver.resolve(molecule.clone()), Ok(Solution::Determined(expected.clone())));
        prop_assert_eq!(resolver.resolve_into(&mut molecule), Ok(Solution::Determined(())));
        prop_assert_eq!(&molecule, &expected);
        prop_assert_eq!(resolver.resolve(molecule), Ok(Solution::Determined(expected)));
    }

    #[test]
    fn test_resolver_resolve_discharge(
        skeleton in aromatic_skeleton(),
        typing in any::<bool>(),
    ) {
        let raw = skeleton.molecule();
        let mut editor = raw.clone().edit();
        for atom in raw.atoms().iter() {
            let degree = atom.neighbors().count() as i64;
            let mut target = editor.atom_mut(atom.id());
            let constraints = &mut target.attributes_mut().constraints;
            constraints.set(AtomConstraintForm::Valence(NumForm::Lit(degree)));
            constraints.set(AtomConstraintForm::Degree(NumForm::Lit(degree + 1)));
        }
        editor.constraints_mut().push(Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None, sum: NumForm::Lit(1),
        }));
        let input = editor.finish().unwrap();
        let mut model = ChemistryModel {
            valence: if typing { ValenceModel::default() } else { ValenceModel::smiles() },
            ..ChemistryModel::default()
        };
        model.valence.tie_break = ValenceTieBreak::MostSaturated;
        let resolver = Resolver::with_config(&model, ResolveConfig {
            isotope: IsotopePolicy::Natural, ..ResolveConfig::default()
        });
        let mut results = Vec::new();
        if let Ok(Solution::Determined(output)) = resolver.aromaticity.resolve(input.clone()) {
            results.push((true, false, false, output));
        }
        if let Ok(Solution::Determined(output)) = resolver.stereo.resolve(input.clone()) {
            results.push((false, true, false, output));
        }
        if let Ok(Solution::Determined(output)) = resolver.resolve(input.clone()) {
            results.push((true, true, true, output));
        }
        for (aromatic, stereo, completes_atoms, output) in results {
            prop_assert_eq!(output.constraints(), input.constraints());
            let added_members: BTreeSet<_> = output.aromatic_systems().iter()
                .skip(input.aromatic_systems().count()).flat_map(|system| system.atom_ids()).collect();
            let added_atoms: BTreeSet<_> = output.stereo_atoms().iter()
                .skip(input.stereo_atoms().count()).map(|site| site.site_id()).collect();
            let added_bonds: BTreeSet<_> = output.stereo_bonds().iter()
                .skip(input.stereo_bonds().count()).map(|site| site.site_id()).collect();
            for atom in input.atoms().iter() {
                let mut expected = atom.attributes().constraints.clone();
                if aromatic && added_members.contains(&atom.id()) {
                    expected.remove(AtomConstraintKey::AromaticValence);
                }
                if stereo && added_atoms.contains(&atom.id()) {
                    expected.remove(AtomConstraintKey::TetrahedralStereo);
                }
                if completes_atoms {
                    expected.remove(AtomConstraintKey::Valence);
                }
                prop_assert_eq!(&output.atom(atom.id()).attributes().constraints, &expected);
            }
            for bond in input.bonds().iter() {
                let mut expected = bond.attributes().constraints.clone();
                if aromatic && !bond.is_in_aromatic_system() {
                    if let Some(BondConstraintForm::Aromatic(asserted)) = expected.get(BondConstraintKey::Aromatic) {
                        let realized = BooleanForm::Lit(output.bond(bond.id()).is_in_aromatic_system());
                        if asserted.meet(&realized).is_some() {
                            expected.remove(BondConstraintKey::Aromatic);
                        }
                    }
                }
                if stereo && added_bonds.contains(&bond.id()) {
                    expected.remove(BondConstraintKey::CisTransStereo);
                }
                prop_assert_eq!(&output.bond(bond.id()).attributes().constraints, &expected);
            }
        }
    }

    #[test]
    fn test_resolver_resolve_renderings(
        skeleton in aromatic_skeleton(),
        typing in any::<bool>(),
    ) {
        let mut model = ChemistryModel {
            valence: if typing { ValenceModel::default() } else { ValenceModel::smiles() },
            ..ChemistryModel::default()
        };
        model.valence.tie_break = ValenceTieBreak::MostSaturated;
        let resolver = Resolver::with_config(&model, ResolveConfig {
            isotope: IsotopePolicy::Natural, ..ResolveConfig::default()
        });
        let smiles: Molecule = Smiles::parse(&skeleton.smiles()).unwrap()
            .as_table_ir().try_into_ir(&()).unwrap();
        let molfile = parse_mol_to_ir(&skeleton.molfile()).unwrap();
        for input in [smiles, molfile] {
            let Ok(Solution::Determined(output)) = resolver.resolve(input) else { continue; };
            for atom in output.atoms().iter() {
                prop_assert_eq!(atom.attributes().constraints.get(AtomConstraintKey::AromaticValence), None);
                prop_assert_eq!(atom.attributes().constraints.get(AtomConstraintKey::TetrahedralStereo), None);
            }
            for bond in output.bonds().iter() {
                prop_assert_eq!(bond.attributes().constraints.get(BondConstraintKey::Aromatic), None);
                prop_assert_eq!(bond.attributes().constraints.get(BondConstraintKey::CisTransStereo), None);
            }
        }
    }

    #[test]
    fn test_resolver_project_roundtrip(
        skeleton in aromatic_skeleton(),
        typing in any::<bool>(),
    ) {
        let mut model = ChemistryModel {
            valence: if typing { ValenceModel::default() } else { ValenceModel::smiles() },
            ..ChemistryModel::default()
        };
        model.valence.tie_break = ValenceTieBreak::MostSaturated;
        let resolver = Resolver::with_config(&model, ResolveConfig {
            isotope: IsotopePolicy::Natural, ..ResolveConfig::default()
        });
        let smiles: Molecule = Smiles::parse(&skeleton.smiles()).unwrap()
            .as_table_ir().try_into_ir(&()).unwrap();
        let molfile = parse_mol_to_ir(&skeleton.molfile()).unwrap();
        for input in [skeleton.molecule(), smiles, molfile] {
            let Ok(Solution::Determined(source)) = resolver.resolve(input) else { continue; };
            let Solution::Determined(projected) = resolver.project(source.clone(),
                ProjectFlags::AROMATICITY | ProjectFlags::STEREO).unwrap() else {
                return Err(TestCaseError::fail("projection of a resolved skeleton was not determined"));
            };
            let mut expected = source.clone().edit();
            for system in source.aromatic_systems().iter() {
                for atom in system.atom_ids() {
                    expected.atom_mut(atom).attributes_mut().constraints.remove(AtomConstraintKey::AromaticValence);
                }
            }
            for bond in source.bonds().iter().filter(|bond| bond.is_in_aromatic_system()) {
                expected.bond_mut(bond.id()).attributes_mut().constraints.remove(BondConstraintKey::Aromatic);
            }
            for site in source.stereo_atoms().iter() {
                expected.atom_mut(site.site_id()).attributes_mut().constraints.remove(AtomConstraintKey::TetrahedralStereo);
            }
            for site in source.stereo_bonds().iter() {
                expected.bond_mut(site.site_id()).attributes_mut().constraints.remove(BondConstraintKey::CisTransStereo);
            }
            let expected = expected.finish().unwrap();
            let Solution::Determined(restored) = resolver.resolve(projected).unwrap() else {
                return Err(TestCaseError::fail("projected skeleton did not resolve under the same model"));
            };
            prop_assert!(restored.framed_eq(&expected));
        }
    }

    #[test]
    fn test_resolver_resolve_aromatic_constraints(
        skeleton in aromatic_skeleton(),
        rule in 0usize..3,
        keep in any::<bool>(),
    ) {
        let source = skeleton.molecule();
        let concrete = Molecule::from_entries(MoleculeEntries {
            atoms: source.atoms().iter().map(|atom| AtomForm {
                element: ElementForm::Lit(Element::C),
                isotope_mass: IsotopeMassForm::Natural,
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(3 - atom.neighbors().count() as i64),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                constraints: AtomConstraintForm::AromaticValence(
                    AromaticValenceForm::Aromatic(NumForm::Lit(1)),
                ).into(),
            }).collect(),
            bonds: source.bonds().iter().map(|bond| {
                let [a, b] = bond.atom_ids();
                (a, b, BondForm {
                    order: NumForm::Lit(1), charge: NumForm::Lit(0),
                    unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                    constraints: bond.attributes().constraints.clone(),
                })
            }).collect(),
            ..MoleculeEntries::default()
        });
        let model = ChemistryModel {
            aromaticity: AromaticityModel {
                scope: ElementScope::Any,
                rule: match rule {
                    0 => AromaticityRule::Hueckel { ring_limits: RingLimits::default() },
                    1 => AromaticityRule::Hmo { stabilization_threshold: 0.1 },
                    _ => AromaticityRule::Clar,
                },
                tie_break: AromaticityTieBreak::Strict,
            },
            ..ChemistryModel::default()
        };
        let config = AromaticityResolveConfig {
            aromatic_valence_failure: if keep { AromaticityFailurePolicy::Keep } else { AromaticityFailurePolicy::Error },
            aromatic_system_failure: if keep { AromaticityFailurePolicy::Keep } else { AromaticityFailurePolicy::Error },
            ..AromaticityResolveConfig::default()
        };
        let standalone = AromaticityResolver::with_config(&model.aromaticity, config);
        let composite = Resolver::with_config(&model, ResolveConfig {
            aromaticity: config, ..ResolveConfig::default()
        });
        for input in [source, concrete] {
            let mut outputs = Vec::new();
            if let Ok(Solution::Determined(output)) = standalone.resolve(input.clone()) {
                outputs.push(output);
            }
            if let Ok(Solution::Determined(output)) = composite.resolve(input.clone()) {
                outputs.push(output);
            }
            for output in outputs {
                for system in output.aromatic_systems().iter().skip(input.aromatic_systems().count()) {
                    let ElectronCountsForm::Lit(electrons) = &system.attributes().electrons else {
                        return Err(TestCaseError::fail("new system has nonliteral electron counts"));
                    };
                    let members: Vec<_> = system.atom_ids().collect();
                    prop_assert_eq!(electrons.len(), members.len());
                    for (&atom, &electrons) in members.iter().zip(electrons) {
                        let literal = AtomConstraintForm::AromaticValence(
                            AromaticValenceForm::Aromatic(NumForm::Lit(electrons)),
                        );
                        let evidence = input.atom(atom).constraints().evidenced(AtomConstraintKey::AromaticValence);
                        prop_assert_eq!(evidence.and_then(|value| value.meet(&literal)), Some(literal));
                    }
                    for bond in input.bonds().iter() {
                        if matches!(bond.attributes().constraints.aromatic(), BooleanForm::Lit(false)) {
                            prop_assert!(!bond.atom_ids().iter().all(|atom| members.contains(atom)));
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_aromaticity_resolver_select_with_report_search_independence(
        scenario in select_scenario()
    ) {
        let resolver = AromaticityResolver::with_config(&scenario.model, scenario.config);
        let actual = resolver
            .select_with_report(
                &scenario.molecule,
                ResolveState {
                    completions: scenario.completions.clone(),
                    ..ResolveState::default()
                },
                scenario.tie_break,
            )
            .expect("Hückel perception is infallible");
        let expected = exhaustive_select(&scenario);
        prop_assert_eq!(actual, expected.clone());
        let without_report = resolver.select(
            &scenario.molecule,
            ResolveState {
                completions: scenario.completions.clone(),
                ..ResolveState::default()
            },
            scenario.tie_break,
        ).expect("Hückel perception is infallible");
        prop_assert_eq!(without_report, expected.map(|mut state| {
            state.tie_breaks.clear();
            state
        }));
    }
}

proptest! {
    #[test]
    fn test_resolver_resolve_carbon(
        radicals in prop::collection::vec(any::<bool>(), 1..24),
        isotope in prop_oneof![Just(IsotopeMassForm::Natural), Just(IsotopeMassForm::Lit(13))],
        typing in any::<bool>(),
        most_saturated in any::<bool>(),
    ) {
        let size = radicals.len();
        let atoms: Vec<_> = radicals.iter().enumerate().map(|(i, &radical)| {
            let valence = i64::from(i > 0) + i64::from(i + 1 < size);
            AtomForm {
                element: ElementForm::Lit(Element::C),
                isotope_mass: isotope.clone(),
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(4 - valence - i64::from(radical)),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm {
                    count: NumForm::Lit(i64::from(radical)),
                    multiplicity: NumForm::Lit(1 + i64::from(radical)),
                },
                constraints: Default::default(),
            }
        }).collect();
        let bonds = (1..size).map(|i| (AtomId((i-1) as u32), AtomId(i as u32), BondForm {
            order: NumForm::Lit(1), charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: Default::default(),
        })).collect();
        let original = Molecule::from_entries(MoleculeEntries {atoms, bonds, ..Default::default()});
        let mut input = original.clone();
        for index in 0..input.atoms().count() {
            let mut view = input.atom_mut(AtomId(index as u32));
            let atom = view.attributes_mut();
            atom.lone_pairs = NumForm::Undetermined;
            atom.unpaired_electrons = UnpairedElectronsForm::default();
        }
        let model = ChemistryModel {valence: ValenceModel {
            tie_break: if most_saturated {ValenceTieBreak::MostSaturated} else {ValenceTieBreak::Strict},
            ..if typing {ValenceModel::default()} else {ValenceModel::smiles()}
        }, ..Default::default()};
        let resolver = Resolver::new(&model);
        prop_assert_eq!(resolver.resolve(input.clone()), Ok(Solution::Determined(original.clone())));
        prop_assert_eq!(resolver.resolve_into(&mut input), Ok(Solution::Determined(())));
        prop_assert_eq!(input, original);
    }

    #[test]
    fn test_resolver_resolve_into_heteroatoms(
        sites in prop::collection::vec((0usize..8, 0i64..3), 1..24),
    ) {
        let states = [
            (Element::C, 0, 4, 0),
            (Element::C, 1, 3, 0),
            (Element::C, -1, 3, 1),
            (Element::N, 0, 3, 1),
            (Element::N, 1, 4, 0),
            (Element::O, 0, 2, 2),
            (Element::O, -1, 1, 3),
            (Element::S, 0, 2, 2),
        ];
        let atoms = sites.iter().enumerate().map(|(i, &(choice, order))| {
            let valence = if i > 0 {order} else {0}
                + sites.get(i + 1).map_or(0, |&(_, next_order)| next_order);
            let eligible: Vec<_> = states.iter().filter(|&&(_, _, capacity, _)| capacity >= valence).collect();
            let &(element, charge, capacity, lone_pairs) = eligible[choice % eligible.len()];
            AtomForm {
                element: ElementForm::Lit(element),
                isotope_mass: IsotopeMassForm::Natural,
                charge: NumForm::Lit(charge),
                implicit_hydrogens: NumForm::Lit(capacity - valence),
                lone_pairs: NumForm::Lit(lone_pairs),
                unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                constraints: Default::default(),
            }
        }).collect();
        let bonds = sites.iter().enumerate().skip(1).filter(|&(_, &(_, order))| order > 0)
            .map(|(i, &(_, order))| (AtomId((i - 1) as u32), AtomId(i as u32), BondForm {
                order: NumForm::Lit(order), charge: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: Default::default(),
            })).collect();
        let source = Molecule::from_entries(MoleculeEntries {atoms, bonds, ..Default::default()});
        let mut input = source.clone();
        for index in 0..input.atoms().count() {
            let mut view = input.atom_mut(AtomId(index as u32));
            let atom = view.attributes_mut();
            atom.lone_pairs = NumForm::Undetermined;
            atom.unpaired_electrons = UnpairedElectronsForm::default();
        }
        for valence in [ValenceModel::smiles(), ValenceModel::default()] {
            for tie_break in [ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated] {
                let model = ChemistryModel {
                    valence: ValenceModel {tie_break, ..valence.clone()},
                    ..Default::default()
                };
                let resolver = Resolver::new(&model);
                let mut molecule = input.clone();
                prop_assert_eq!(resolver.resolve_into(&mut molecule), Ok(Solution::Determined(())));
                prop_assert_eq!(&molecule, &source);
            }
        }
    }

}
