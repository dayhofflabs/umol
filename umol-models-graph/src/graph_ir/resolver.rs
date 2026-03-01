//! Resolution of TableIR molecules to GraphIR molecules.

use std::collections::HashMap;

use super::config::{ResolveConfig, TopologyResolveFlags, ValenceMatchPolicy, ValenceStrategyKind};
use super::dative::DativeBond;
use super::error::ResolutionError;
use super::molecule::{AtomIndex, MoleculeBuilder};
use super::multicenter::{MulticenterBond, MulticenterContribution, MulticenterSet};
use super::noncovalent::NoncovalentBond;
use super::valence::ValenceValidator;
use super::{AtomBuilder, Bond, Molecule};
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
/// 2. **Valence** — match atoms against valid valence states; for aromatic atoms,
///    enumerate candidate states with σ-bond sum and π-contribution.
/// 3. **Aromaticity** — select from valence candidates by ring membership and
///    Kekulé/Hückel feasibility.
/// 4. **Stereochemistry** — validate chiral centers and bond stereo against the
///    resolved topology and bond orders.
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
            builder.add_bond_unchecked(a, b, Bond::from_table_bond(bond)?);
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
    if !config.valence.enabled {
        return Ok(());
    }

    let validator = match config.valence.strategy {
        ValenceStrategyKind::AtomTyping => {
            ValenceValidator::AtomTyping(config.valence.registry.clone())
        }
        ValenceStrategyKind::Counts => ValenceValidator::Counts,
    };

    let atom_indices: Vec<AtomIndex> = builder.atom_indices().collect();
    for atom_index in atom_indices {
        let candidates = validator.candidates_for(builder, atom_index);

        if candidates.is_empty() {
            if config.valence.no_match_policy == ValenceMatchPolicy::Ignore {
                continue;
            }
            let element = builder.atom(atom_index).expect("atom_index must be valid").element();
            return Err(ResolutionError::ValenceNoMatch(format!(
                "atom {:?} at index {} has no valence match",
                element,
                atom_index.index()
            )));
        }

        if candidates.len() > 1 && config.valence.ambiguous_policy != ValenceMatchPolicy::Ignore {
            let element = builder.atom(atom_index).expect("atom_index must be valid").element();
            return Err(ResolutionError::ValenceAmbiguous(format!(
                "atom {:?} at index {} has {} valence matches",
                element,
                atom_index.index(),
                candidates.len()
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
    _builder: &mut MoleculeBuilder,
    config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    if !config.aromaticity.enabled {
        return Ok(());
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
    use umol_data::Element;

    use super::*;
    use crate::registry;
    use super::super::valence::AtomTypeRegistry;
    use crate::table_ir::atom::Atom as TableAtom;
    use crate::table_ir::bond::{Bond as TableBond, BondOrder};
    use crate::table_ir::Molecule as TableMolecule;

    fn empty_atom(element: Element) -> TableAtom {
        TableAtom::from_element(element)
    }

    fn config_with(registry: AtomTypeRegistry) -> ResolveConfig {
        let mut config = ResolveConfig::default();
        config.valence.registry = registry;
        config
    }

    #[test]
    fn resolve_molecule_h2_succeeds() {
        let mut table = TableMolecule::empty();
        table.atoms.push(empty_atom(Element::H));
        table.atoms.push(empty_atom(Element::H));
        table.bonds.push(TableBond::new(0, 1, BondOrder::Single));

        let resolved = resolve_molecule_with(&table, &config_with(registry!["[H+0v1]"]));
        assert!(resolved.is_ok());
    }

    #[test]
    #[should_panic(expected = "counts-based valence is not implemented yet")]
    fn resolve_molecule_counts_strategy_panics() {
        let mut table = TableMolecule::empty();
        table.atoms.push(empty_atom(Element::H));
        table.atoms.push(empty_atom(Element::H));
        table.bonds.push(TableBond::new(0, 1, BondOrder::Single));

        let mut config = config_with(registry!["[H+0v1]"]);
        config.valence.strategy = ValenceStrategyKind::Counts;

        let _ = resolve_molecule_with(&table, &config);
    }

    #[test]
    fn resolve_molecule_no_match_errors() {
        let mut table = TableMolecule::empty();
        table.atoms.push(empty_atom(Element::C));

        let resolved = resolve_molecule_with(&table, &config_with(AtomTypeRegistry::new()));
        assert!(matches!(resolved, Err(ResolutionError::ValenceNoMatch(_))));
    }

    #[test]
    fn resolve_molecule_ambiguous_errors() {
        let mut table = TableMolecule::empty();
        table.atoms.push(empty_atom(Element::C));
        table.atoms.push(empty_atom(Element::H));
        table.atoms.push(empty_atom(Element::H));
        table.atoms.push(empty_atom(Element::H));
        table.bonds.push(TableBond::new(0, 1, BondOrder::Single));
        table.bonds.push(TableBond::new(0, 2, BondOrder::Single));
        table.bonds.push(TableBond::new(0, 3, BondOrder::Single));

        let resolved = resolve_molecule_with(
            &table,
            &config_with(registry!["[H+0v1]", "[C+0v3]", "[C+1v3]"]),
        );
        assert!(matches!(
            resolved,
            Err(ResolutionError::ValenceAmbiguous(_))
        ));
    }
}
