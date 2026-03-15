//! Resolution of TableIR molecules to GraphIR molecules.

pub mod aromaticity;
pub mod valence;

use std::collections::HashMap;

use self::aromaticity::AromaticityModel;
use self::valence::ValenceValidator;
use super::atom::AtomBuilder;
use super::bond::BondBuilder;
use super::config::{
    AromaticityHintPolicy, ResolveConfig, RingMethod, TopologyResolveFlags, ValenceMatchPolicy,
    ValenceStrategy,
};
use super::dative::DativeBond;
use super::error::ResolutionError;
use super::molecule::builder::MoleculeBuilder;
use super::molecule::{AtomIndex, Molecule};
use super::multicenter::{MulticenterBond, MulticenterContribution, MulticenterSet};
use super::noncovalent::NoncovalentBond;
use super::rings::MoleculeRings;
use crate::table_ir::{BondDonation, Molecule as TableMolecule};

/// Resolve a TableIR molecule to a GraphIR molecule using default configuration.
pub fn resolve_molecule(molecule: &TableMolecule) -> Result<Molecule, ResolutionError> {
    resolve_molecule_with(molecule, &ResolveConfig::default())
}

/// Resolve a TableIR molecule to a GraphIR molecule with the given configuration.
///
/// Populates a `MoleculeBuilder` from the TableIR data, runs the resolution
/// phases in order, then builds the final `Molecule`.
///
/// 1. **Topology** — build graph from TableIR, validate indices, self-loops,
///    parallel edges, connectivity.
/// 2. **Valence** — match atoms against valid valence states
/// 3. **Aromaticity** — select from valence candidates by ring membership
/// 4. **Stereochemistry** — validate chiral centers and bond stereo
pub fn resolve_molecule_with(
    molecule: &TableMolecule,
    config: &ResolveConfig,
) -> Result<Molecule, ResolutionError> {
    let mut builder = resolve_topology_with(molecule, config)?;
    resolve_valence_with(&mut builder, config)?;
    resolve_aromaticity_with(&mut builder, config)?;
    resolve_stereo_with(&mut builder, config)?;
    builder.build(config)
}

fn resolve_topology_with(
    molecule: &TableMolecule,
    config: &ResolveConfig,
) -> Result<MoleculeBuilder, ResolutionError> {
    if config.topology.enabled {
        if !config
            .topology
            .flags
            .contains(TopologyResolveFlags::DISCONNECTED_MOLECULES)
            && molecule.component_count() > 1
        {
            return Err(ResolutionError::TopologyDisconnected);
        }

        for (i, bond) in molecule.bonds.iter().enumerate() {
            if bond.atoms.first() == bond.atoms.second() {
                return Err(ResolutionError::TopologySelfLoop(i as u32));
            }
        }

        let mut atom_bonds: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
        for (i, bond) in molecule.bonds.iter().enumerate() {
            let key = (bond.atoms.first(), bond.atoms.second());
            atom_bonds.entry(key).or_default().push(i as u32);
        }
        for (_pair, indices) in atom_bonds {
            if indices.len() >= 2 {
                return Err(ResolutionError::TopologyParallelEdges(
                    indices[0], indices[1],
                ));
            }
        }
    }

    let n = molecule.atom_count();
    let m = molecule.bond_count();
    let mut builder = MoleculeBuilder::with_capacity(n, m);
    let mut node_indices: Vec<AtomIndex> = Vec::with_capacity(n);
    for atom in &molecule.atoms {
        node_indices.push(builder.add_atom(AtomBuilder::from_table_atom(atom)));
    }
    for bond in &molecule.bonds {
        let a = node_indices[bond.atoms.first() as usize];
        let b = node_indices[bond.atoms.second() as usize];
        if bond.noncovalent.is_some() {
            builder.add_noncovalent_bond(NoncovalentBond::from_table_bond(bond, &node_indices));
        } else if matches!(
            bond.donation,
            Some(BondDonation::Donating | BondDonation::Accepting)
        ) {
            builder.add_dative_bond(DativeBond::from_table_bond(bond, &node_indices));
        } else {
            builder.add_bond_unchecked(a, b, BondBuilder::from_table_bond(bond));
        }
    }

    for mc in &molecule.multicenter_bonds {
        let sets: Vec<MulticenterSet> = mc
            .contributions()
            .iter()
            .map(|contrib| {
                let contributions: Vec<MulticenterContribution> = contrib
                    .atoms()
                    .iter()
                    .map(|&idx| {
                        if (idx as usize) >= n {
                            Err(ResolutionError::AtomIndexOutOfRange(idx as u32))
                        } else {
                            Ok(MulticenterContribution::topology_only(
                                node_indices[idx as usize],
                            ))
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(MulticenterSet::topology_only(
                    contributions.iter().map(|c| c.atom()),
                ))
            })
            .collect::<Result<Vec<_>, _>>()?;
        builder.add_multicenter_bond(MulticenterBond::new(sets));
    }

    Ok(builder)
}

fn resolve_valence_with(
    builder: &mut MoleculeBuilder,
    config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    let validator = match config.valence.strategy {
        ValenceStrategy::AtomTyping => {
            ValenceValidator::AtomTyping(config.valence.atom_type_registry.clone())
        }
        ValenceStrategy::Counts => ValenceValidator::Counts {
            table: config.valence.valence_table.clone(),
            enable_implicit_hydrogens: config.valence.enable_implicit_hydrogens,
        },
    };

    let atom_indices: Vec<AtomIndex> = builder.atom_indices().collect();
    for atom_index in atom_indices {
        let candidates = validator.candidates_for(builder, atom_index);

        if candidates.is_empty() {
            if config.valence.no_match_policy == ValenceMatchPolicy::Ignore {
                continue;
            }
            let element = builder
                .atom(atom_index)
                .expect("atom_index must be valid")
                .element();
            return Err(ResolutionError::ValenceNoMatch(format!(
                "atom {:?} at index {} has no valence match",
                element,
                atom_index.index()
            )));
        }

        if candidates.len() > 1 && config.valence.ambiguous_policy != ValenceMatchPolicy::Ignore {
            let element = builder
                .atom(atom_index)
                .expect("atom_index must be valid")
                .element();
            let specs: Vec<String> = candidates.iter().map(|s| s.to_string()).collect();
            return Err(ResolutionError::ValenceAmbiguous(format!(
                "atom {:?} at index {} has {} valence matches: {}",
                element,
                atom_index.index(),
                candidates.len(),
                specs.join(", ")
            )));
        }

        builder
            .atom_mut(atom_index)
            .expect("atom_index from atom_indices must be valid")
            .set_candidates(candidates);
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

    let model = AromaticityModel::from_strategy(&config.aromaticity.perception_strategy);
    let ring_cfg = &config.aromaticity.ring_strategy;
    let rings = match ring_cfg.method {
        RingMethod::GlobalAllCycles => MoleculeRings::from_builder(
            builder,
            ring_cfg.max_ring_size,
            ring_cfg.max_rings_per_component,
        ),
        RingMethod::PiSubgraph => MoleculeRings::from_pi_subgraph(
            builder,
            ring_cfg.max_ring_size,
            ring_cfg.max_rings_per_component,
        ),
    };
    for system in model.aromatic_systems(builder, &rings)? {
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
        let bond = builder.bond(bond_index).expect("bond_index must be valid");
        if bond.aromatic_hint() != Some(true) {
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
    config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    if !config.stereo.enabled {
        return Ok(());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_data::Element;

    use super::super::config_data::AtomTypeRegistry;
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
        config.valence.atom_type_registry = AtomTypeRegistry::new();
        config
    }

    #[fixture]
    fn config_with_h_registry() -> ResolveConfig {
        let mut config = ResolveConfig::default();
        config.valence.atom_type_registry = registry!["{H+0v1}"];
        config
    }

    #[fixture]
    fn config_with_counts_strategy() -> ResolveConfig {
        let mut config = ResolveConfig::default();
        config.valence.strategy = ValenceStrategy::Counts;
        config
    }

    #[fixture]
    fn config_with_ch_registry() -> ResolveConfig {
        let mut config = ResolveConfig::default();
        config.valence.atom_type_registry = registry!["{H+0v1}", "{C+0v3}", "{C+1v3}"];
        config
    }

    #[rstest]
    fn resolve_molecule(h2_molecule: TableMolecule, config_with_h_registry: ResolveConfig) {
        let resolved = resolve_molecule_with(&h2_molecule, &config_with_h_registry);
        assert!(resolved.is_ok());
        let mol = resolved.unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(
            mol.atom(mol.atom_indices().next().unwrap())
                .unwrap()
                .hydrogens(),
            0
        );
    }

    #[rstest]
    fn resolve_molecule_counts_strategy(
        h2_molecule: TableMolecule,
        config_with_counts_strategy: ResolveConfig,
    ) {
        let resolved = resolve_molecule_with(&h2_molecule, &config_with_counts_strategy);
        assert!(resolved.is_ok());
        let mol = resolved.unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(
            mol.atom(mol.atom_indices().next().unwrap())
                .unwrap()
                .hydrogens(),
            0
        );
    }

    // TODO: Add tests for aromaticity phase

    #[rstest]
    fn resolve_molecule_no_match(
        c_molecule: TableMolecule,
        config_with_empty_registry: ResolveConfig,
    ) {
        let resolved = resolve_molecule_with(&c_molecule, &config_with_empty_registry);
        assert!(matches!(resolved, Err(ResolutionError::ValenceNoMatch(_))));
    }

    #[rstest]
    fn resolve_molecule_ambiguous(
        ch3_molecule: TableMolecule,
        config_with_ch_registry: ResolveConfig,
    ) {
        let resolved = resolve_molecule_with(&ch3_molecule, &config_with_ch_registry);
        assert!(matches!(
            resolved,
            Err(ResolutionError::ValenceAmbiguous(_))
        ));
    }
}
