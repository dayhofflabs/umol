//! Resolution of molecular structure.

use std::collections::HashMap;

use smallvec::SmallVec;
use umol_shared::element::Element;
use umol_shared::spin::{SpinState, MAX_UNPAIRED_ELECTRONS};

use super::aromaticity::AromaticityModel;
use super::atom::Atom;
use super::atom_pattern::{AtomPattern, HydrogenPattern, Pattern};
use super::config::{
    AromaticityHintPolicy, AromaticityStrategy, ResolveConfig, TopologyResolveFlags,
    ValenceMatchPolicy, ValenceStrategy,
};
use super::config_data::{NormalValenceTable, ValenceTable};
use super::error::{GraphIrError, ResolutionError};
use super::molecule::{AtomIndex, BondIndex};
use super::molecule_builder::MoleculeBuilder;
use super::rings::{RingEnumerator, RingFamily};
use crate::atom::AromaticValence;

/// Resolve a molecular structure to a ground `Molecule` using default configuration.
pub fn resolve_molecule(builder: &mut MoleculeBuilder) -> Result<(), GraphIrError> {
    resolve_molecule_with(builder, &ResolveConfig::default())
}

/// Run all resolution phases on a `MoleculeBuilder` in place.
///
/// Modifies the builder by populating atom candidates and aromatic systems.
/// Call `builder.build(config)` afterwards to ground to a `Molecule`.
///
/// Phases:
/// 1. **Topology** — validate self-loops, parallel edges, connectivity
/// 2. **Valence** — match atoms against valid valence states
/// 3. **Aromaticity** — select from valence candidates by ring membership
/// 4. **Stereochemistry** — validate chiral centers and bond stereo
pub fn resolve_molecule_with(
    builder: &mut MoleculeBuilder,
    config: &ResolveConfig,
) -> Result<(), GraphIrError> {
    resolve_topology_with(builder, config)?;
    resolve_valence_with(builder, config)?;
    resolve_aromaticity_with(builder, config)?;
    resolve_stereo_with(builder, config)?;
    Ok(())
}

fn resolve_topology_with(
    builder: &mut MoleculeBuilder,
    config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    if !config
        .topology
        .flags
        .contains(TopologyResolveFlags::DISCONNECTED_MOLECULES)
        && builder.component_count() > 1
    {
        return Err(ResolutionError::TopologyDisconnected);
    }

    for bond_index in builder.bond_indices() {
        let (a, b) = builder
            .bond_atom_indices(bond_index)
            .expect("bond_index from bond_indices must be valid");
        if a == b {
            return Err(ResolutionError::TopologySelfLoop(bond_index.index() as u32));
        }
    }

    let mut bond_pairs: HashMap<(AtomIndex, AtomIndex), Vec<BondIndex>> = HashMap::new();
    for bond_index in builder.bond_indices() {
        let (a, b) = builder
            .bond_atom_indices(bond_index)
            .expect("bond_index from bond_indices must be valid");
        let key = if a <= b { (a, b) } else { (b, a) };
        bond_pairs.entry(key).or_default().push(bond_index);
    }
    for indices in bond_pairs.values() {
        if indices.len() >= 2 {
            return Err(ResolutionError::TopologyParallelEdges(
                indices[0].index() as u32,
                indices[1].index() as u32,
            ));
        }
    }

    Ok(())
}

fn resolve_valence_with(
    builder: &mut MoleculeBuilder,
    config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    let atom_indices: Vec<AtomIndex> = builder.atom_indices().collect();

    for &atom_index in &atom_indices {
        let candidates = valence_candidates(&config.valence.strategy, builder, atom_index);

        if candidates.is_empty() {
            if config.valence.no_match_policy == ValenceMatchPolicy::Ignore {
                continue;
            }
            let element = builder
                .atom(atom_index)
                .expect("atom_index must be valid")
                .element();
            return Err(ResolutionError::ValenceNoMatch(format!(
                "atom {:?} at index {}",
                element,
                atom_index.index()
            )));
        }

        if candidates.len() > 1 && config.valence.ambiguous_policy != ValenceMatchPolicy::Ignore {
            let element = builder
                .atom(atom_index)
                .expect("atom_index must be valid")
                .element();
            let specs: Vec<String> = candidates.iter().map(ToString::to_string).collect();
            return Err(ResolutionError::ValenceAmbiguous(format!(
                "atom {:?} at index {} has {} valence matches: {}",
                element,
                atom_index.index(),
                candidates.len(),
                specs.join(", ")
            )));
        }

        builder
            .set_atom_candidates(atom_index, candidates)
            .expect("atom_index from atom_indices must be valid");
    }

    Ok(())
}

fn resolve_aromaticity_with(
    builder: &mut MoleculeBuilder,
    config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    let has_hints = builder
        .atom_indices()
        .any(|i| builder.atom_aromatic_hint(i));
    if !has_hints {
        return Ok(());
    }

    let model = AromaticityModel::new(&config.aromaticity.aromaticity_strategy);
    let ring_family = match config.aromaticity.aromaticity_strategy {
        AromaticityStrategy::Clar => RingFamily::InducedBenzenoid,
        AromaticityStrategy::HueckelRule { .. } | AromaticityStrategy::Hmo { .. } => {
            RingFamily::Simple
        }
    };
    let enumerator = RingEnumerator::new(ring_family, &config.aromaticity.enumeration_strategy);
    let rings = enumerator.enumerate_builder(builder);
    for system in model
        .aromatic_systems(builder, &rings)
        .map_err(|e| ResolutionError::AromaticityInconsistent(e.to_string()))?
    {
        builder.add_aromatic_system(system);
    }

    validate_aromatic_hints(builder, config.aromaticity.hint_policy)?;

    Ok(())
}

fn validate_aromatic_hints(
    builder: &MoleculeBuilder,
    policy: AromaticityHintPolicy,
) -> Result<(), ResolutionError> {
    if policy == AromaticityHintPolicy::Ignore {
        return Ok(());
    }

    for atom_index in builder.atom_indices() {
        if builder.atom_aromatic_hint(atom_index) && !builder.atom_has_aromatic_systems(atom_index)
        {
            let element = builder
                .atom(atom_index)
                .expect("atom_index must be valid")
                .element();
            return Err(ResolutionError::AromaticityInconsistent(format!(
                "atom {:?} at index {} has aromatic hint but is not in any detected aromatic system",
                element,
                atom_index.index()
            )));
        }
    }

    for bond_index in builder.bond_indices() {
        if builder.bond_aromatic_hint(bond_index) != Some(true) {
            continue;
        }
        let (a, b) = builder
            .bond_atom_indices(bond_index)
            .expect("bond_index must be valid");
        if !builder.atom_has_aromatic_systems(a) || !builder.atom_has_aromatic_systems(b) {
            return Err(ResolutionError::AromaticityInconsistent(format!(
                "bond at index {} has aromatic hint but its endpoints are not both in an aromatic system",
                bond_index.index()
            )));
        }
    }

    Ok(())
}

fn resolve_stereo_with(
    _builder: &mut MoleculeBuilder,
    _config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    Ok(())
}

fn valence_candidates(
    strategy: &ValenceStrategy,
    builder: &MoleculeBuilder,
    atom_index: AtomIndex,
) -> SmallVec<[Atom; 4]> {
    match strategy {
        ValenceStrategy::AtomTyping { registry } => {
            let mut pattern = AtomPattern::from_builder_atom(builder, atom_index);
            if pattern.implicit_hydrogens == HydrogenPattern::Normal {
                let Some(hydrogens) = infer_normal_implicit_hydrogens(builder, atom_index) else {
                    return SmallVec::new();
                };
                pattern.implicit_hydrogens = HydrogenPattern::Is(hydrogens);
            }
            todo!("old pipeline: registry.candidates_for removed, use solver.rs")
        }
        ValenceStrategy::Counts {
            table,
            allow_implicit_hydrogens,
        } => counts_candidates(table, *allow_implicit_hydrogens, builder, atom_index),
    }
}

fn infer_normal_implicit_hydrogens(
    builder: &MoleculeBuilder,
    atom_index: AtomIndex,
) -> Option<u8> {
    let atom = builder.atom(atom_index).expect("atom_index must be valid");
    let element = atom.element();
    let charge = match atom.charge {
        Pattern::Is(c) => c,
        Pattern::Any => 0,
    };
    let explicit_valence = builder.atom_bond_order_sum(atom_index);

    if builder.atom_aromatic_hint(atom_index) {
        if charge != 0 {
            return None;
        }
        return if element == Element::C {
            Some(3_u8.saturating_sub(explicit_valence))
        } else if matches!(
            element,
            Element::B
                | Element::N
                | Element::O
                | Element::P
                | Element::S
                | Element::Se
                | Element::As
        ) {
            Some(0)
        } else {
            None
        };
    }

    let normal_valence =
        NormalValenceTable::default_table().normal_valence_for(element, charge)?;
    Some(normal_valence.saturating_sub(explicit_valence))
}

fn counts_candidates(
    table: &ValenceTable,
    allow_implicit_hydrogens: bool,
    builder: &MoleculeBuilder,
    atom_index: AtomIndex,
) -> SmallVec<[Atom; 4]> {
    let atom = builder.atom(atom_index).expect("atom_index must be valid");
    let element = atom.element();
    let charge = match atom.charge {
        Pattern::Is(c) => c,
        Pattern::Any => 0,
    };
    let valence = builder.atom_bond_order_sum(atom_index);
    let (donated_pairs, accepted_pairs) = builder.atom_dative_bond_order_sums(atom_index);

    let entry = match table.entry(element) {
        Some(e) => e,
        None => return SmallVec::new(),
    };

    if builder.atom_aromatic_hint(atom_index) {
        let aromatic_valences = if charge != 0 {
            element
                .shift(-charge)
                .and_then(|e| table.entry(e))
                .map(|e| e.allowed_aromatic_valences.as_slice())
                .unwrap_or(entry.allowed_aromatic_valences.as_slice())
        } else {
            entry.allowed_aromatic_valences.as_slice()
        };
        return build_aromatic_candidates(
            aromatic_valences,
            element,
            charge,
            valence,
            donated_pairs,
            accepted_pairs,
            allow_implicit_hydrogens,
            builder.atom_has_normal_implicit_hydrogens(atom_index),
            atom,
        );
    }

    let implicit_hydrogens = if let HydrogenPattern::Is(h) = &atom.implicit_hydrogens {
        *h
    } else if allow_implicit_hydrogens {
        match table.compute_implicit_hydrogens(element, charge, valence) {
            Some(h) => h,
            None => return SmallVec::new(),
        }
    } else {
        0
    };
    try_build_atom(
        element,
        charge,
        implicit_hydrogens,
        valence,
        donated_pairs,
        accepted_pairs,
        AromaticValence::NotAromatic,
        atom,
    )
    .into_iter()
    .collect()
}

#[allow(clippy::too_many_arguments)]
fn build_aromatic_candidates(
    allowed_aromatic_valences: &[u8],
    element: Element,
    charge: i8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    allow_implicit_hydrogens: bool,
    normal_implicit_hydrogens: bool,
    atom: &AtomPattern,
) -> SmallVec<[Atom; 4]> {
    if allowed_aromatic_valences.is_empty() {
        return SmallVec::new();
    }

    let effective_electrons = (element.valence_electrons() as i16) - (charge as i16);
    let mut candidates = SmallVec::new();

    for &a in allowed_aromatic_valences {
        let sigma_budget = effective_electrons - (a as i16);
        if sigma_budget < valence as i16 {
            continue;
        }
        let implicit_hydrogens = match &atom.implicit_hydrogens {
            HydrogenPattern::Is(h) => *h,
            HydrogenPattern::Normal if normal_implicit_hydrogens => {
                let Some(h) =
                    infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
                else {
                    continue;
                };
                h
            }
            HydrogenPattern::Normal => continue,
            HydrogenPattern::Any if normal_implicit_hydrogens => {
                let Some(h) =
                    infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
                else {
                    continue;
                };
                h
            }
            HydrogenPattern::Any => {
                if allow_implicit_hydrogens {
                    (sigma_budget - valence as i16) as u8
                } else {
                    0
                }
            }
        };
        if implicit_hydrogens > 1 {
            continue;
        }
        let total_sigma = valence + implicit_hydrogens;
        let remaining = effective_electrons - total_sigma as i16 - a as i16;
        if remaining < 0 || remaining % 2 != 0 {
            continue;
        }
        if let Some(atom_out) = try_build_atom(
            element,
            charge,
            implicit_hydrogens,
            valence,
            donated_pairs,
            accepted_pairs,
            AromaticValence::Valence(a),
            atom,
        ) {
            candidates.push(atom_out);
        }
    }

    candidates
}

fn infer_normal_aromatic_implicit_hydrogens(
    element: Element,
    charge: i8,
    valence: u8,
) -> Option<u8> {
    if charge != 0 {
        return None;
    }
    if element == Element::C {
        Some(3 - valence)
    } else if matches!(
        element,
        Element::B | Element::N | Element::O | Element::P | Element::S | Element::Se | Element::As
    ) {
        Some(0)
    } else {
        None
    }
}

#[allow(clippy::too_many_arguments)]
fn try_build_atom(
    element: Element,
    charge: i8,
    implicit_hydrogens: u8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    aromatic_valence: AromaticValence,
    atom: &AtomPattern,
) -> Option<Atom> {
    let total_valence = valence + implicit_hydrogens;
    let num_electrons = (element.valence_electrons() as i16) - (charge as i16);
    let unassigned = num_electrons - (total_valence as i16) - (aromatic_valence.valence() as i16);
    if unassigned < 0 {
        return None;
    }

    let (unpaired, lone_pairs) = match (atom.unpaired_electrons, atom.lone_pairs) {
        (Pattern::Any, Pattern::Any) => ((unassigned % 2) as u8, (unassigned / 2) as u8),
        (Pattern::Is(unpaired), Pattern::Any) => {
            let remaining = unassigned - (unpaired as i16);
            if remaining < 0 || remaining % 2 != 0 {
                return None;
            }
            (unpaired, (remaining / 2) as u8)
        }
        (Pattern::Any, Pattern::Is(lone_pairs)) => {
            let remaining = unassigned - (2 * lone_pairs as i16);
            if remaining < 0 {
                return None;
            }
            (remaining as u8, lone_pairs)
        }
        (Pattern::Is(unpaired), Pattern::Is(lone_pairs)) => {
            if (unpaired as i16) + (2 * lone_pairs as i16) != unassigned {
                return None;
            }
            (unpaired, lone_pairs)
        }
    };
    if unpaired > MAX_UNPAIRED_ELECTRONS {
        return None;
    }

    let spin = match atom.multiplicity {
        Pattern::Is(m) => SpinState::try_new(unpaired, m).ok()?,
        Pattern::Any => SpinState::max_multiplicity(unpaired)?,
    };
    Atom::try_new(
        element,
        None,
        charge,
        implicit_hydrogens,
        lone_pairs,
        unpaired,
        spin.multiplicity(),
        valence,
        donated_pairs,
        accepted_pairs,
        aromatic_valence,
        0,
    )
    .ok()
}

// TODO: Remove TableMolecule from the tests, just work on builder directly
#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::element::Element;

    use super::super::config::ValenceStrategy;
    use super::super::config_data::{AtomTypeRegistry, ValenceTable};
    use super::super::molecule_builder::MoleculeBuilder;
    use super::*;
    use crate::registry;
    use crate::table_ir::atom::Atom as TableAtom;
    use crate::table_ir::bond::{Bond as TableBond, BondOrder};
    use crate::table_ir::Molecule as TableMolecule;

    #[fixture]
    fn h_atom() -> TableAtom {
        TableAtom::from_element(Element::H)
    }

    #[fixture]
    fn c_atom() -> TableAtom {
        TableAtom::from_element(Element::C)
    }

    #[fixture]
    fn h2_molecule(h_atom: TableAtom) -> TableMolecule {
        let mut table = TableMolecule::empty();
        table.atoms.push(h_atom.clone());
        table.atoms.push(h_atom.clone());
        table.bonds.push(TableBond::new(0, 1, BondOrder::Single));
        table
    }

    #[fixture]
    fn c_molecule(c_atom: TableAtom) -> TableMolecule {
        let mut table = TableMolecule::empty();
        table.atoms.push(c_atom);
        table
    }

    #[fixture]
    fn ch3_molecule(c_atom: TableAtom, h_atom: TableAtom) -> TableMolecule {
        let mut table = TableMolecule::empty();
        table.atoms.push(c_atom.clone());
        table.atoms.push(h_atom.clone());
        table.atoms.push(h_atom.clone());
        table.atoms.push(h_atom.clone());
        table.bonds.push(TableBond::new(0, 1, BondOrder::Single));
        table.bonds.push(TableBond::new(0, 2, BondOrder::Single));
        table.bonds.push(TableBond::new(0, 3, BondOrder::Single));
        table
    }

    #[fixture]
    fn config_with_empty_registry() -> ResolveConfig {
        let mut config = ResolveConfig::default();
        config.valence.strategy = ValenceStrategy::AtomTyping {
            registry: AtomTypeRegistry::new(),
        };
        config
    }

    #[fixture]
    fn config_with_h_registry() -> ResolveConfig {
        let mut config = ResolveConfig::default();
        config.valence.strategy = ValenceStrategy::AtomTyping {
            registry: registry!["H #v"],
        };
        config
    }

    #[fixture]
    fn config_with_counts_strategy() -> ResolveConfig {
        let mut config = ResolveConfig::default();
        config.valence.strategy = ValenceStrategy::Counts {
            table: ValenceTable::default_table().clone(),
            allow_implicit_hydrogens: true,
        };
        config
    }

    #[fixture]
    fn config_with_ch_registry() -> ResolveConfig {
        let mut config = ResolveConfig::default();
        config.valence.strategy = ValenceStrategy::AtomTyping {
            registry: registry!["H #v", "C #u #v3", "C #c+ #v3"],
        };
        config
    }

    #[rstest]
    fn resolve_molecule(h2_molecule: TableMolecule, config_with_h_registry: ResolveConfig) {
        let mut builder = MoleculeBuilder::from_table_molecule(&h2_molecule);
        resolve_molecule_with(&mut builder, &config_with_h_registry).unwrap();
        let mol = builder.build(&config_with_h_registry).unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(
            mol.atom(mol.atom_indices().next().unwrap())
                .unwrap()
                .implicit_hydrogens(),
            0
        );
    }

    #[rstest]
    fn resolve_molecule_counts_strategy(
        h2_molecule: TableMolecule,
        config_with_counts_strategy: ResolveConfig,
    ) {
        let mut builder = MoleculeBuilder::from_table_molecule(&h2_molecule);
        resolve_molecule_with(&mut builder, &config_with_counts_strategy).unwrap();
        let mol = builder.build(&config_with_counts_strategy).unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(
            mol.atom(mol.atom_indices().next().unwrap())
                .unwrap()
                .implicit_hydrogens(),
            0
        );
    }

    // TODO: Add tests for aromaticity phase

    #[rstest]
    fn resolve_molecule_no_match(
        c_molecule: TableMolecule,
        config_with_empty_registry: ResolveConfig,
    ) {
        let mut builder = MoleculeBuilder::from_table_molecule(&c_molecule);
        let result = resolve_molecule_with(&mut builder, &config_with_empty_registry);
        assert!(matches!(
            result,
            Err(GraphIrError::Resolution(ResolutionError::ValenceNoMatch(_)))
        ));
    }

    #[rstest]
    fn resolve_molecule_ambiguous(
        ch3_molecule: TableMolecule,
        config_with_ch_registry: ResolveConfig,
    ) {
        let mut builder = MoleculeBuilder::from_table_molecule(&ch3_molecule);
        let result = resolve_molecule_with(&mut builder, &config_with_ch_registry);
        assert!(matches!(
            result,
            Err(GraphIrError::Resolution(ResolutionError::ValenceAmbiguous(
                _
            )))
        ));
    }
}
