//! Resolution of molecular structure.

use std::collections::HashMap;

use super::aromaticity::AromaticityModel;
use super::config::{
    AromaticityHintPolicy, AromaticityStrategy, ResolveConfig, TopologyResolveFlags,
    ValenceMatchPolicy,
};
use super::error::{GraphIrError, ResolutionError};
use super::molecule::{AtomIndex, BondIndex};
use super::molecule_builder::MoleculeBuilder;
use super::rings::{RingEnumerator, RingFamily};
use super::valence::ValenceMatcher;

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
    let matcher = ValenceMatcher::new(&config.valence.strategy);
    let atom_indices: Vec<AtomIndex> = builder.atom_indices().collect();
    for atom_index in atom_indices {
        let candidates = matcher.candidates_for(builder, atom_index);

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

// TODO: Remove TableMolecule from the tests, just work on builder directly
#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::Element;

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
