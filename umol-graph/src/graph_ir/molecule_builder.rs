//! GraphIR molecule builder.

use std::collections::{HashMap, HashSet};
use std::fmt;
use std::str::FromStr;

use petgraph::prelude::*;
use petgraph::stable_graph::StableGraph;
use petgraph::visit::{EdgeRef, NodeIndexable};
use smallvec::SmallVec;
use umol_data::{SpinMultiplicity, SpinState};
use umol_edn::{DeError, Edn, EdnMap, EdnMapHelper, EdnSet, FormatConfig, FromEdn, ToEdn};

use super::aromaticity::{AromaticContribution, AromaticSystem};
use super::atom::Atom;
use super::atom_pattern::{AtomPattern, HydrogenPattern};
use super::bond_pattern::BondPattern;
use super::config::ResolveConfig;
use super::dative::DativeBond;
use super::error::{GraphIrError, ResolutionError};
use super::molecule::{
    AromaticSystemIndex, AtomIndex, BondIndex, DativeBondIndex, Molecule, MulticenterBondIndex,
    NoncovalentBondIndex,
};
use super::multicenter::{MulticenterBond, MulticenterContribution, MulticenterSet};
use super::noncovalent::NoncovalentBond;
use crate::algorithms::biconnected_components;
use crate::atom::AromaticValence;
use crate::dsl::ast::{FromAst, ToAst};
use crate::dsl::bond::BondAst;
use crate::dsl::config::MoleculeDslConfig;
use crate::dsl::error::{LoweringError, ParseError};
use crate::dsl::molecule::{
    AromaticSystem as AromaticSystemAst, DativeBond as DativeBondAst,
    LocalizedBond as LocalizedBondAst, MoleculeAst, MoleculeAstWrapper,
    MulticenterBond as MulticenterBondAst, NoncovalentBond as NoncovalentBondAst,
};
use crate::table_ir::atom::ImplicitHydrogens;
use crate::table_ir::bond::BondOrder;
use crate::table_ir::{BondDonation, Molecule as TableMolecule};

/// Transient resolution state carried by `MoleculeBuilder` during the resolution
/// pipeline. Not part of the final `Molecule` or the `MoleculeAst`.
#[derive(Debug, Clone, Default)]
pub struct ResolutionContext {
    pub atom_candidates: HashMap<AtomIndex, SmallVec<[Atom; 4]>>,
    pub atom_aromatic_hints: HashMap<AtomIndex, bool>,
    pub bond_aromatic_hints: HashMap<BondIndex, bool>,
    pub atom_normal_implicit_hydrogens: HashSet<AtomIndex>,
}

impl<'de> FromEdn<'de> for ResolutionContext {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let m = match edn {
            Edn::Map(m) => m,
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "ResolutionContext map",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        let mut h = EdnMapHelper::new(m);

        let atom_candidates = match h.optional::<Edn>("atom-candidates")? {
            Some(Edn::Map(m)) => parse_atom_candidates(&m)?,
            Some(other) => {
                return Err(DeError::TypeMismatch {
                    expected: "map",
                    got: other.kind(),
                    path: vec!["atom-candidates".into()],
                });
            }
            None => HashMap::new(),
        };

        let atom_aromatic_hints = match h.optional::<Edn>("atom-aromatic-hints")? {
            Some(Edn::Map(m)) => parse_node_map(&m, edn_to_bool)?,
            Some(other) => {
                return Err(DeError::TypeMismatch {
                    expected: "map",
                    got: other.kind(),
                    path: vec!["atom-aromatic-hints".into()],
                });
            }
            None => HashMap::new(),
        };

        let bond_aromatic_hints = match h.optional::<Edn>("bond-aromatic-hints")? {
            Some(Edn::Map(m)) => parse_edge_map(&m, edn_to_bool)?,
            Some(other) => {
                return Err(DeError::TypeMismatch {
                    expected: "map",
                    got: other.kind(),
                    path: vec!["bond-aromatic-hints".into()],
                });
            }
            None => HashMap::new(),
        };

        // Accept both EDN sets and vectors (serde legacy format)
        let atom_normal_implicit_hydrogens = match h.optional::<Edn>("atom-normal-implicit-hydrogens")? {
            Some(edn) => parse_index_set(&edn)?,
            None => HashSet::new(),
        };

        Ok(Self {
            atom_candidates,
            atom_aromatic_hints,
            bond_aromatic_hints,
            atom_normal_implicit_hydrogens,
        })
    }
}

impl ToEdn for ResolutionContext {
    fn to_edn(&self) -> Edn<'static> {
        let mut m = EdnMap::with_capacity(4);

        let mut ac = EdnMap::with_capacity(self.atom_candidates.len());
        for (k, v) in &self.atom_candidates {
            let atoms: Vec<_> = v.iter().map(|a| a.to_edn()).collect();
            ac.insert(Edn::Int(k.index() as i64), Edn::Vector(atoms.into()));
        }
        m.insert(Edn::keyword("atom-candidates"), Edn::Map(ac));

        m.insert(
            Edn::keyword("atom-aromatic-hints"),
            Edn::Map(index_map_to_edn(&self.atom_aromatic_hints)),
        );
        m.insert(
            Edn::keyword("bond-aromatic-hints"),
            Edn::Map(edge_map_to_edn(&self.bond_aromatic_hints)),
        );

        let mut set = EdnSet::new();
        for idx in &self.atom_normal_implicit_hydrogens {
            set.insert(Edn::Int(idx.index() as i64));
        }
        m.insert(
            Edn::keyword("atom-normal-implicit-hydrogens"),
            Edn::Set(set),
        );

        Edn::Map(m)
    }
}

fn edn_to_usize(edn: &Edn<'_>) -> Result<usize, DeError> {
    match edn {
        Edn::Int(n) => Ok(*n as usize),
        other => Err(DeError::TypeMismatch {
            expected: "integer",
            got: other.kind(),
            path: Vec::new(),
        }),
    }
}

fn edn_to_bool(edn: &Edn<'_>) -> Result<bool, DeError> {
    match edn {
        Edn::Bool(b) => Ok(*b),
        other => Err(DeError::TypeMismatch {
            expected: "boolean",
            got: other.kind(),
            path: Vec::new(),
        }),
    }
}

fn parse_node_map<V>(
    map: &EdnMap<'_>,
    parse_value: fn(&Edn<'_>) -> Result<V, DeError>,
) -> Result<HashMap<AtomIndex, V>, DeError> {
    let mut result = HashMap::with_capacity(map.len());
    for (k, v) in map.iter() {
        result.insert(NodeIndex::new(edn_to_usize(k)?), parse_value(v)?);
    }
    Ok(result)
}

fn parse_edge_map<V>(
    map: &EdnMap<'_>,
    parse_value: fn(&Edn<'_>) -> Result<V, DeError>,
) -> Result<HashMap<BondIndex, V>, DeError> {
    let mut result = HashMap::with_capacity(map.len());
    for (k, v) in map.iter() {
        result.insert(EdgeIndex::new(edn_to_usize(k)?), parse_value(v)?);
    }
    Ok(result)
}

fn parse_atom_candidates(map: &EdnMap<'_>) -> Result<HashMap<AtomIndex, SmallVec<[Atom; 4]>>, DeError> {
    let mut result = HashMap::with_capacity(map.len());
    for (k, v) in map.iter() {
        let atoms: Vec<Atom> = Vec::from_edn(v)?;
        result.insert(NodeIndex::new(edn_to_usize(k)?), SmallVec::from_vec(atoms));
    }
    Ok(result)
}

/// Parse EDN set or vector of integers into a HashSet of AtomIndex.
fn parse_index_set(edn: &Edn<'_>) -> Result<HashSet<AtomIndex>, DeError> {
    let mut result = HashSet::new();
    match edn {
        Edn::Set(s) => {
            for e in s.iter() {
                result.insert(NodeIndex::new(edn_to_usize(e)?));
            }
        }
        Edn::Vector(v) => {
            for e in v.iter() {
                result.insert(NodeIndex::new(edn_to_usize(e)?));
            }
        }
        other => {
            return Err(DeError::TypeMismatch {
                expected: "set or vector",
                got: other.kind(),
                path: Vec::new(),
            });
        }
    }
    Ok(result)
}

fn index_map_to_edn(map: &HashMap<AtomIndex, bool>) -> EdnMap<'static> {
    let mut m = EdnMap::with_capacity(map.len());
    for (k, v) in map {
        m.insert(Edn::Int(k.index() as i64), Edn::Bool(*v));
    }
    m
}

fn edge_map_to_edn(map: &HashMap<BondIndex, bool>) -> EdnMap<'static> {
    let mut m = EdnMap::with_capacity(map.len());
    for (k, v) in map {
        m.insert(Edn::Int(k.index() as i64), Edn::Bool(*v));
    }
    m
}

impl fmt::Display for ResolutionContext {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let edn = self.to_edn();
        edn.to_string_with(&FormatConfig::default()).fmt(f)
    }
}

impl FromStr for ResolutionContext {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let tree = umol_edn::read_string(s)?;
        Self::from_edn(&tree).map_err(|e| ParseError::EdnParse(e.to_string()))
    }
}

/// Builder for constructing a `Molecule`. Carries `AtomPattern` nodes during
/// resolution phases; `build()` finalizes each atom and produces a `Molecule`.
///
/// Used both by the resolution pipeline (from TableIR) and for manual
/// molecule construction.
#[derive(Debug, Clone)]
pub struct MoleculeBuilder {
    graph: StableGraph<AtomPattern, BondPattern, Undirected, u32>,
    resolution: ResolutionContext,
    dative_bonds: Vec<DativeBond>,
    aromatic_systems: Vec<AromaticSystem>,
    multicenter_bonds: Vec<MulticenterBond>,
    noncovalent_bonds: Vec<NoncovalentBond>,
    charge: Option<i8>,
    spin: Option<SpinState>,
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self {
            graph: StableGraph::default(),
            resolution: ResolutionContext::default(),
            dative_bonds: Vec::new(),
            aromatic_systems: Vec::new(),
            multicenter_bonds: Vec::new(),
            noncovalent_bonds: Vec::new(),
            charge: None,
            spin: None,
        }
    }

    pub fn with_capacity(atom_capacity: usize, bond_capacity: usize) -> Self {
        Self {
            graph: StableGraph::with_capacity(atom_capacity, bond_capacity),
            resolution: ResolutionContext {
                atom_candidates: HashMap::with_capacity(atom_capacity),
                atom_aromatic_hints: HashMap::with_capacity(atom_capacity),
                bond_aromatic_hints: HashMap::with_capacity(bond_capacity),
                atom_normal_implicit_hydrogens: HashSet::with_capacity(atom_capacity),
            },
            dative_bonds: Vec::new(),
            aromatic_systems: Vec::new(),
            multicenter_bonds: Vec::new(),
            noncovalent_bonds: Vec::new(),
            charge: None,
            spin: None,
        }
    }

    pub fn from_molecule(molecule: &Molecule) -> Self {
        let mut builder = Self::with_capacity(molecule.atom_count(), molecule.bond_count());

        let mut atom_map: HashMap<AtomIndex, AtomIndex> = HashMap::new();
        for atom_idx in molecule.atom_indices() {
            let atom = molecule.atom(atom_idx).expect("atom index must be valid");
            let new_idx = builder.add_atom(AtomPattern::from_atom(atom));
            builder
                .set_atom_candidates(new_idx, SmallVec::from_elem(*atom, 1))
                .expect("newly added atom index must be valid");
            atom_map.insert(atom_idx, new_idx);
        }

        for bond_idx in molecule.bond_indices() {
            let bond = molecule.bond(bond_idx).expect("bond index must be valid");
            let (a, b) = molecule
                .bond_atom_indices(bond_idx)
                .expect("bond index must be valid");
            let new_a = *atom_map.get(&a).expect("source atom must be mapped");
            let new_b = *atom_map.get(&b).expect("source atom must be mapped");
            builder.add_bond_unchecked(new_a, new_b, BondPattern::from_bond(bond));
        }

        builder.dative_bonds = molecule.dative_bonds().cloned().collect();
        builder.aromatic_systems = molecule.aromatic_systems().cloned().collect();
        builder.multicenter_bonds = molecule.multicenter_bonds().cloned().collect();
        builder.noncovalent_bonds = molecule.noncovalent_bonds().cloned().collect();
        builder.charge = Some(molecule.charge());
        builder.spin = Some(molecule.spin());
        builder
    }

    /// Build a `MoleculeBuilder` from a `TableMolecule` without running any topology checks.
    ///
    /// Topology validation (self-loops, parallel edges, connectivity) is a separate pass
    /// that operates on the returned builder; see `resolve_topology_with` in `resolution.rs`.
    pub fn from_table_molecule(molecule: &TableMolecule) -> Self {
        let n = molecule.atom_count();
        let m = molecule.bond_count();
        let mut builder = Self::with_capacity(n, m);
        let mut node_indices: Vec<AtomIndex> = Vec::with_capacity(n);

        for atom in &molecule.atoms {
            let idx = builder.add_atom(AtomPattern::from_table_atom(atom));
            if let Some(aromatic) = atom.aromatic {
                builder
                    .set_atom_aromatic_hint(idx, aromatic)
                    .expect("newly added atom index must be valid");
            }
            if atom.implicit_hydrogens == Some(ImplicitHydrogens::Normal) {
                builder
                    .set_atom_normal_implicit_hydrogens(idx)
                    .expect("newly added atom index must be valid");
            }
            node_indices.push(idx);
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
                let bond_idx = builder.add_bond_unchecked(a, b, BondPattern::from_table_bond(bond));
                if bond.order == BondOrder::Aromatic {
                    builder.set_bond_aromatic_hint(bond_idx, true);
                }
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
                            assert!(
                                (idx as usize) < node_indices.len(),
                                "multicenter bond references atom index {} which is out of range \
                                 (molecule has {} atoms)",
                                idx,
                                node_indices.len()
                            );
                            MulticenterContribution::topology_only(node_indices[idx as usize])
                        })
                        .collect();
                    MulticenterSet::topology_only(contributions.iter().map(|c| c.atom()))
                })
                .collect();
            builder.add_multicenter_bond(MulticenterBond::new(sets));
        }

        builder
    }

    pub fn biconnected_components(&self) -> Vec<Vec<AtomIndex>> {
        let mut atoms: Vec<AtomIndex> = self.atom_indices().collect();
        atoms.sort_unstable();
        if atoms.is_empty() {
            return Vec::new();
        }

        let atom_to_id: HashMap<AtomIndex, usize> = atoms
            .iter()
            .copied()
            .enumerate()
            .map(|(i, a)| (a, i))
            .collect();
        let adj = self.adjacency_list();
        let mut adj_int: Vec<Vec<usize>> = vec![Vec::new(); atoms.len()];
        for &atom in &atoms {
            let mut neighbors = adj
                .get(&atom)
                .map(|ns| {
                    ns.iter()
                        .filter_map(|&n| atom_to_id.get(&n).copied())
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            neighbors.sort_unstable();
            neighbors.dedup();
            let u = atom_to_id[&atom];
            adj_int[u] = neighbors;
        }

        biconnected_components(atoms.len(), &adj_int)
            .into_iter()
            .map(|component| component.into_iter().map(|i| atoms[i]).collect())
            .collect()
    }

    // Atoms
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Number of connected components in the covalent bond graph.
    pub fn component_count(&self) -> usize {
        let mut visited: HashSet<AtomIndex> = HashSet::new();
        let mut count = 0;
        for start in self.graph.node_indices() {
            if visited.contains(&start) {
                continue;
            }
            count += 1;
            let mut stack = vec![start];
            while let Some(node) = stack.pop() {
                if !visited.insert(node) {
                    continue;
                }
                for edge in self.graph.edges(node) {
                    let neighbor = if edge.source() == node {
                        edge.target()
                    } else {
                        edge.source()
                    };
                    if !visited.contains(&neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }
        count
    }

    pub fn atom_indices(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.node_indices()
    }

    pub fn atoms(&self) -> impl Iterator<Item = &AtomPattern> + '_ {
        self.graph.node_weights()
    }

    pub fn atom(&self, index: AtomIndex) -> Option<&AtomPattern> {
        self.graph.node_weight(index)
    }

    pub fn atom_mut(&mut self, index: AtomIndex) -> Option<&mut AtomPattern> {
        self.graph.node_weight_mut(index)
    }

    pub fn add_atom(&mut self, atom: impl Into<AtomPattern>) -> AtomIndex {
        self.graph.add_node(atom.into())
    }

    /// Add a fully-resolved atom as both pattern and sole candidate.
    pub fn add_resolved_atom(&mut self, atom: Atom) -> AtomIndex {
        let idx = self.add_atom(AtomPattern::from_atom(&atom));
        self.resolution
            .atom_candidates
            .insert(idx, SmallVec::from_elem(atom, 1));
        idx
    }

    pub fn remove_atom(&mut self, index: AtomIndex) -> Option<AtomPattern> {
        self.resolution.atom_candidates.remove(&index);
        self.resolution.atom_aromatic_hints.remove(&index);
        self.resolution
            .atom_normal_implicit_hydrogens
            .remove(&index);
        self.graph.remove_node(index)
    }

    pub fn replace_atom(
        &mut self,
        index: AtomIndex,
        atom: impl Into<AtomPattern>,
    ) -> Option<AtomPattern> {
        self.graph
            .node_weight_mut(index)
            .map(|old| std::mem::replace(old, atom.into()))
    }

    pub fn set_atom_candidates(
        &mut self,
        index: AtomIndex,
        candidates: SmallVec<[Atom; 4]>,
    ) -> Option<()> {
        if !self.graph.contains_node(index) {
            return None;
        }
        self.resolution.atom_candidates.insert(index, candidates);
        Some(())
    }

    pub fn atom_candidates(&self, index: AtomIndex) -> &[Atom] {
        self.resolution
            .atom_candidates
            .get(&index)
            .map_or(&[], |v| v.as_slice())
    }

    pub fn atom_candidates_mut(&mut self, index: AtomIndex) -> Option<&mut SmallVec<[Atom; 4]>> {
        self.resolution.atom_candidates.get_mut(&index)
    }

    pub fn add_atom_candidate(&mut self, index: AtomIndex, candidate: Atom) -> Option<()> {
        if !self.graph.contains_node(index) {
            return None;
        }
        let entry = self.resolution.atom_candidates.entry(index).or_default();
        if !entry.contains(&candidate) {
            entry.push(candidate);
        }
        Some(())
    }

    pub fn set_atom_aromatic_hint(&mut self, index: AtomIndex, hint: bool) -> Option<()> {
        if !self.graph.contains_node(index) {
            return None;
        }
        self.resolution.atom_aromatic_hints.insert(index, hint);
        Some(())
    }

    pub fn atom_explicit_aromatic_hint(&self, index: AtomIndex) -> Option<bool> {
        self.resolution.atom_aromatic_hints.get(&index).copied()
    }

    pub fn set_bond_aromatic_hint(&mut self, index: BondIndex, hint: bool) {
        self.resolution.bond_aromatic_hints.insert(index, hint);
    }

    pub fn bond_aromatic_hint(&self, index: BondIndex) -> Option<bool> {
        self.resolution.bond_aromatic_hints.get(&index).copied()
    }

    pub fn set_atom_normal_implicit_hydrogens(&mut self, index: AtomIndex) -> Option<()> {
        if !self.graph.contains_node(index) {
            return None;
        }
        self.resolution.atom_normal_implicit_hydrogens.insert(index);
        Some(())
    }

    pub fn clear_atom_normal_implicit_hydrogens(&mut self, index: AtomIndex) {
        self.resolution
            .atom_normal_implicit_hydrogens
            .remove(&index);
    }

    pub fn atom_has_normal_implicit_hydrogens(&self, index: AtomIndex) -> bool {
        self.resolution
            .atom_normal_implicit_hydrogens
            .contains(&index)
    }

    // Atom properties
    fn atom_candidate_property<T>(
        &self,
        index: AtomIndex,
        getter: impl Fn(&Atom) -> T,
    ) -> Option<T> {
        self.resolution
            .atom_candidates
            .get(&index)
            .and_then(|c| c.iter().next().map(getter))
    }

    pub fn atom_valence(&self, index: AtomIndex) -> u8 {
        self.atom_candidate_property(index, |c| c.valence())
            .unwrap_or(0)
    }

    pub fn atom_aromatic_valence(&self, index: AtomIndex) -> u8 {
        self.resolution
            .atom_candidates
            .get(&index)
            .and_then(|candidates| {
                candidates.iter().find_map(|c| match c.aromatic_valence() {
                    AromaticValence::Valence(n) => Some(n),
                    AromaticValence::NotAromatic => None,
                })
            })
            .unwrap_or(0)
    }

    // Atom aromatic hints
    /// Returns true if this atom should be treated as aromatic based on
    /// its own aromatic_hint or any incident bond's aromatic_hint.
    pub fn atom_aromatic_hint(&self, index: AtomIndex) -> bool {
        if self.atom_explicit_aromatic_hint(index) == Some(true) {
            return true;
        }
        self.graph
            .edges(index)
            .any(|e| self.resolution.bond_aromatic_hints.get(&e.id()) == Some(&true))
    }

    /// Atoms that have at least one candidate with a non-None aromatic valence.
    pub fn aromatic_candidate_atoms(&self) -> impl Iterator<Item = AtomIndex> + '_ {
        self.atom_indices()
            .filter(|&atom| self.atom_has_aromatic_candidate(atom))
    }

    pub fn atom_has_aromatic_candidate(&self, index: AtomIndex) -> bool {
        self.resolution
            .atom_candidates
            .get(&index)
            .map(|candidates| {
                candidates
                    .iter()
                    .any(|c| c.aromatic_valence().is_aromatic())
            })
            .unwrap_or(false)
    }

    // Bonds
    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    pub fn bond_indices(&self) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edge_indices()
    }

    pub fn bond(&self, index: BondIndex) -> Option<&BondPattern> {
        self.graph.edge_weight(index)
    }

    pub fn bond_mut(&mut self, index: BondIndex) -> Option<&mut BondPattern> {
        self.graph.edge_weight_mut(index)
    }

    pub fn add_bond(&mut self, a: AtomIndex, b: AtomIndex, bond: BondPattern) -> Option<BondIndex> {
        if !self.graph.contains_node(a) || !self.graph.contains_node(b) {
            return None;
        }
        Some(self.graph.add_edge(a, b, bond))
    }

    pub fn add_bond_unchecked(
        &mut self,
        a: AtomIndex,
        b: AtomIndex,
        bond: BondPattern,
    ) -> BondIndex {
        debug_assert!(
            self.graph.contains_node(a),
            "atom index {:?} not in builder",
            a
        );
        debug_assert!(
            self.graph.contains_node(b),
            "atom index {:?} not in builder",
            b
        );
        self.graph.add_edge(a, b, bond)
    }

    pub fn remove_bond(&mut self, index: BondIndex) -> Option<BondPattern> {
        self.graph.remove_edge(index)
    }

    pub fn replace_bond(&mut self, index: BondIndex, bond: BondPattern) -> Option<BondPattern> {
        self.graph
            .edge_weight_mut(index)
            .map(|old| std::mem::replace(old, bond))
    }

    // Dative bonds
    pub fn dative_bond_count(&self) -> usize {
        self.dative_bonds.len()
    }

    pub fn dative_bond_indices(&self) -> impl Iterator<Item = DativeBondIndex> + '_ {
        (0..self.dative_bond_count()).map(|i| DativeBondIndex(i as u32))
    }

    pub fn dative_bonds(&self) -> impl Iterator<Item = &DativeBond> + '_ {
        self.dative_bonds.iter()
    }

    pub fn dative_bond(&self, index: DativeBondIndex) -> Option<&DativeBond> {
        self.dative_bonds.get(index.index())
    }

    pub fn dative_bond_mut(&mut self, index: DativeBondIndex) -> Option<&mut DativeBond> {
        self.dative_bonds.get_mut(index.index())
    }

    pub fn add_dative_bond(&mut self, bond: DativeBond) {
        self.dative_bonds.push(bond);
    }

    pub fn remove_dative_bond(&mut self, index: DativeBondIndex) -> Option<DativeBond> {
        let i = index.index();
        if i >= self.dative_bonds.len() {
            return None;
        }
        Some(self.dative_bonds.remove(i))
    }

    pub fn replace_dative_bond(
        &mut self,
        index: DativeBondIndex,
        bond: DativeBond,
    ) -> Option<DativeBond> {
        self.dative_bonds
            .get_mut(index.index())
            .map(|b| std::mem::replace(b, bond))
    }

    // Aromatic systems
    pub fn aromatic_system_count(&self) -> usize {
        self.aromatic_systems.len()
    }

    pub fn aromatic_system_indices(&self) -> impl Iterator<Item = AromaticSystemIndex> + '_ {
        (0..self.aromatic_system_count()).map(|i| AromaticSystemIndex(i as u32))
    }

    pub fn aromatic_systems(&self) -> impl Iterator<Item = &AromaticSystem> + '_ {
        self.aromatic_systems.iter()
    }

    pub fn aromatic_system(&self, index: AromaticSystemIndex) -> Option<&AromaticSystem> {
        self.aromatic_systems.get(index.index())
    }

    pub fn aromatic_system_mut(
        &mut self,
        index: AromaticSystemIndex,
    ) -> Option<&mut AromaticSystem> {
        self.aromatic_systems.get_mut(index.index())
    }

    pub fn add_aromatic_system(&mut self, system: AromaticSystem) {
        self.aromatic_systems.push(system);
    }

    pub fn clear_aromatic_systems(&mut self) {
        self.aromatic_systems.clear();
    }

    pub fn remove_aromatic_system(&mut self, index: AromaticSystemIndex) -> Option<AromaticSystem> {
        let i = index.index();
        if i >= self.aromatic_systems.len() {
            return None;
        }
        Some(self.aromatic_systems.remove(i))
    }

    pub fn replace_aromatic_system(
        &mut self,
        index: AromaticSystemIndex,
        system: AromaticSystem,
    ) -> Option<AromaticSystem> {
        self.aromatic_systems
            .get_mut(index.index())
            .map(|s| std::mem::replace(s, system))
    }

    // Multicenter bonds
    pub fn multicenter_bond_count(&self) -> usize {
        self.multicenter_bonds.len()
    }

    pub fn multicenter_bond_indices(&self) -> impl Iterator<Item = MulticenterBondIndex> + '_ {
        (0..self.multicenter_bond_count()).map(|i| MulticenterBondIndex(i as u32))
    }

    pub fn multicenter_bonds(&self) -> impl Iterator<Item = &MulticenterBond> + '_ {
        self.multicenter_bonds.iter()
    }

    pub fn multicenter_bond(&self, index: MulticenterBondIndex) -> Option<&MulticenterBond> {
        self.multicenter_bonds.get(index.index())
    }

    pub fn multicenter_bond_mut(
        &mut self,
        index: MulticenterBondIndex,
    ) -> Option<&mut MulticenterBond> {
        self.multicenter_bonds.get_mut(index.index())
    }

    pub fn add_multicenter_bond(&mut self, bond: MulticenterBond) {
        self.multicenter_bonds.push(bond);
    }

    pub fn remove_multicenter_bond(
        &mut self,
        index: MulticenterBondIndex,
    ) -> Option<MulticenterBond> {
        let i = index.index();
        if i >= self.multicenter_bonds.len() {
            return None;
        }
        Some(self.multicenter_bonds.remove(i))
    }

    pub fn replace_multicenter_bond(
        &mut self,
        index: MulticenterBondIndex,
        bond: MulticenterBond,
    ) -> Option<MulticenterBond> {
        self.multicenter_bonds
            .get_mut(index.index())
            .map(|b| std::mem::replace(b, bond))
    }

    // Non-covalent bonds
    pub fn noncovalent_bond_count(&self) -> usize {
        self.noncovalent_bonds.len()
    }

    pub fn noncovalent_bond_indices(&self) -> impl Iterator<Item = NoncovalentBondIndex> + '_ {
        (0..self.noncovalent_bond_count()).map(|i| NoncovalentBondIndex(i as u32))
    }

    pub fn noncovalent_bonds(&self) -> impl Iterator<Item = &NoncovalentBond> + '_ {
        self.noncovalent_bonds.iter()
    }

    pub fn noncovalent_bond(&self, index: NoncovalentBondIndex) -> Option<&NoncovalentBond> {
        self.noncovalent_bonds.get(index.index())
    }

    pub fn noncovalent_bond_mut(
        &mut self,
        index: NoncovalentBondIndex,
    ) -> Option<&mut NoncovalentBond> {
        self.noncovalent_bonds.get_mut(index.index())
    }

    pub fn add_noncovalent_bond(&mut self, bond: NoncovalentBond) {
        self.noncovalent_bonds.push(bond);
    }

    pub fn remove_noncovalent_bond(
        &mut self,
        index: NoncovalentBondIndex,
    ) -> Option<NoncovalentBond> {
        let i = index.index();
        if i >= self.noncovalent_bonds.len() {
            return None;
        }
        Some(self.noncovalent_bonds.remove(i))
    }

    pub fn replace_noncovalent_bond(
        &mut self,
        index: NoncovalentBondIndex,
        bond: NoncovalentBond,
    ) -> Option<NoncovalentBond> {
        self.noncovalent_bonds
            .get_mut(index.index())
            .map(|b| std::mem::replace(b, bond))
    }

    // Molecular charge and spin
    pub fn charge(&self) -> Option<i8> {
        self.charge
    }

    pub fn spin(&self) -> Option<SpinState> {
        self.spin
    }

    pub fn set_charge(&mut self, charge: i8) {
        self.charge = Some(charge);
    }

    pub fn clear_charge(&mut self) {
        self.charge = None;
    }

    pub fn set_spin(&mut self, spin: SpinState) {
        self.spin = Some(spin);
    }

    pub fn clear_spin(&mut self) {
        self.spin = None;
    }

    pub fn resolution_context(&self) -> &ResolutionContext {
        &self.resolution
    }

    pub fn set_resolution_context(&mut self, ctx: ResolutionContext) {
        self.resolution = ctx;
    }

    // Atom-atom relationships
    pub fn adjacency_list(&self) -> HashMap<AtomIndex, Vec<AtomIndex>> {
        let mut adj = HashMap::with_capacity(self.graph.node_count());
        for atom in self.graph.node_indices() {
            adj.insert(atom, Vec::new());
        }
        for bond in self.graph.edge_indices() {
            let (a, b) = self.graph.edge_endpoints(bond).unwrap();
            adj.get_mut(&a).unwrap().push(b);
            adj.get_mut(&b).unwrap().push(a);
        }
        adj
    }

    pub fn atom_neighbor_indices(&self, index: AtomIndex) -> impl Iterator<Item = AtomIndex> + '_ {
        self.graph.neighbors(index)
    }

    pub fn atom_neighbors(&self, index: AtomIndex) -> impl Iterator<Item = &AtomPattern> + '_ {
        self.graph
            .neighbors(index)
            .map(|n| self.graph.node_weight(n).unwrap())
    }

    // TODO: Add dative and noncovalent neighbors (+indices)
    // TODO: Add aromatic system and multicenter system partners (+indices)

    // Atom-bond relationships
    pub fn atom_bond_count(&self, index: AtomIndex) -> usize {
        self.graph.edges(index).count()
    }

    pub fn atom_bond_indices(&self, index: AtomIndex) -> impl Iterator<Item = BondIndex> + '_ {
        self.graph.edges(index).map(|e| e.id())
    }

    pub fn atom_bonds(&self, index: AtomIndex) -> impl Iterator<Item = &BondPattern> + '_ {
        self.graph.edges(index).map(|e| e.weight())
    }

    pub fn atom_bond_order_sum(&self, index: AtomIndex) -> u8 {
        self.graph.edges(index).map(|e| e.weight().order()).sum()
    }

    pub fn connecting_bond_index(&self, a: AtomIndex, b: AtomIndex) -> Option<BondIndex> {
        self.graph.edges_connecting(a, b).next().map(|e| e.id())
    }

    pub fn connecting_bond(&self, a: AtomIndex, b: AtomIndex) -> Option<&BondPattern> {
        self.graph.edges_connecting(a, b).next().map(|e| e.weight())
    }

    pub fn bond_atom_indices(&self, index: BondIndex) -> Option<(AtomIndex, AtomIndex)> {
        self.graph.edge_endpoints(index)
    }

    pub fn bond_atoms(&self, index: BondIndex) -> Option<(&AtomPattern, &AtomPattern)> {
        self.graph.edge_endpoints(index).map(|(a, b)| {
            (
                self.graph.node_weight(a).unwrap(),
                self.graph.node_weight(b).unwrap(),
            )
        })
    }

    // Atom-dative bond relationships
    pub fn atom_has_dative_bonds(&self, index: AtomIndex) -> bool {
        self.dative_bonds.iter().any(|b| b.contains_atom(index))
    }

    pub fn atom_dative_bond_counts(&self, index: AtomIndex) -> (usize, usize) {
        // TODO: Be forgiving about atoms not in the graph
        debug_assert!(
            self.graph.contains_node(index),
            "atom index {:?} not in builder",
            index
        );
        let mut donated = 0;
        let mut accepted = 0;
        for db in &self.dative_bonds {
            if db.donor() == index {
                donated += 1;
            } else if db.acceptor() == index {
                accepted += 1;
            }
        }
        (donated, accepted)
    }

    pub fn atom_dative_bond_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = DativeBondIndex> + '_ {
        self.dative_bond_indices()
            .filter(move |&i| self.dative_bond(i).unwrap().contains_atom(index))
    }

    pub fn atom_dative_bonds(&self, index: AtomIndex) -> impl Iterator<Item = &DativeBond> + '_ {
        self.dative_bonds().filter(move |b| b.contains_atom(index))
    }

    pub fn atom_dative_bond_order_sums(&self, index: AtomIndex) -> (u8, u8) {
        // TODO: Be forgiving about atoms not in the graph
        debug_assert!(
            self.graph.contains_node(index),
            "atom index {:?} not in builder",
            index
        );

        let mut donated = 0;
        let mut accepted = 0;
        for db in &self.dative_bonds {
            if db.donor() == index {
                donated += db.order();
            } else if db.acceptor() == index {
                accepted += db.order();
            }
        }
        (donated, accepted)
    }

    // Atom-aromatic system relationships
    pub fn atom_has_aromatic_systems(&self, index: AtomIndex) -> bool {
        // TODO: Be forgiving about atoms not in the graph
        debug_assert!(
            self.graph.contains_node(index),
            "atom index {:?} not in builder",
            index
        );
        self.atom_aromatic_systems_indices(index).next().is_some()
    }

    pub fn atom_aromatic_systems_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AromaticSystemIndex> + '_ {
        self.aromatic_system_indices()
            .filter(move |&i| self.aromatic_system(i).unwrap().contains_atom(index))
    }

    pub fn atom_aromatic_systems(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AromaticSystem> + '_ {
        self.aromatic_systems()
            .filter(move |s| s.contains_atom(index))
            .cloned()
    }

    // Atom-multicenter bond relationships
    pub fn atom_has_multicenter_bonds(&self, index: AtomIndex) -> bool {
        // TODO: Be forgiving about atoms not in the graph
        debug_assert!(
            self.graph.contains_node(index),
            "atom index {:?} not in builder",
            index
        );
        self.atom_multicenter_bonds_indices(index).next().is_some()
    }

    pub fn atom_multicenter_bonds_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = MulticenterBondIndex> + '_ {
        self.multicenter_bond_indices()
            .filter(move |&i| self.multicenter_bond(i).unwrap().contains_atom(index))
    }

    pub fn atom_multicenter_bonds(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = MulticenterBond> + '_ {
        self.multicenter_bonds()
            .filter(move |b| b.contains_atom(index))
            .cloned()
    }

    // Atom-noncovalent bond relationships
    pub fn atom_has_noncovalent_bonds(&self, index: AtomIndex) -> bool {
        self.noncovalent_bonds
            .iter()
            .any(|b| b.contains_atom(index))
    }

    pub fn atom_noncovalent_bond_indices(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = NoncovalentBondIndex> + '_ {
        self.noncovalent_bond_indices()
            .filter(move |&i| self.noncovalent_bond(i).unwrap().contains_atom(index))
    }

    pub fn atom_noncovalent_bonds(
        &self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &NoncovalentBond> + '_ {
        self.noncovalent_bonds()
            .filter(move |b| b.contains_atom(index))
    }

    /// Build the final `Molecule` by finalizing each `AtomPattern` into an `Atom`.
    ///
    /// Requires all atoms to have exactly one valence candidate
    /// remaining (i.e., resolution phases must have been run).
    pub fn build(self, config: &ResolveConfig) -> Result<Molecule, GraphIrError> {
        self.build_inner(config).map_err(GraphIrError::from)
    }

    fn build_inner(self, _config: &ResolveConfig) -> Result<Molecule, ResolutionError> {
        let mut graph =
            StableGraph::with_capacity(self.graph.node_count(), self.graph.edge_count());

        let mut index_map = Vec::with_capacity(self.graph.node_bound());
        index_map.resize(self.graph.node_bound(), None);

        for old_idx in self.graph.node_indices() {
            let pattern = self.graph.node_weight(old_idx).unwrap();
            let candidates = self
                .resolution
                .atom_candidates
                .get(&old_idx)
                .ok_or_else(|| {
                    ResolutionError::ValenceNoMatch(format!(
                        "no valence match for {:?}",
                        pattern.element()
                    ))
                })?;
            let candidate = match candidates.as_slice() {
                [] => {
                    return Err(ResolutionError::ValenceNoMatch(format!(
                        "no valence match for {:?}",
                        pattern.element()
                    )))
                }
                [single] => single,
                many => {
                    let specs: Vec<String> = many.iter().map(ToString::to_string).collect();
                    return Err(ResolutionError::ValenceAmbiguous(format!(
                        "{} valence matches for {:?}: {}",
                        many.len(),
                        pattern.element(),
                        specs.join(", ")
                    )));
                }
            };

            // Resolve `HydrogenPattern::Normal` before matching: by this point the registry
            // has selected the candidate, so we substitute Normal with the candidate's
            // concrete hydrogen count.  Normal is a deferred constraint that cannot be
            // evaluated without registry context; leaving it unresolved would silently
            // accept any hydrogen count.
            let resolved_pattern = AtomPattern {
                implicit_hydrogens: match &pattern.implicit_hydrogens {
                    HydrogenPattern::Normal => HydrogenPattern::Is(candidate.implicit_hydrogens()),
                    h => h.clone(),
                },
                ..pattern.clone()
            };
            if !resolved_pattern.matches_atom(candidate) {
                return Err(ResolutionError::ValenceViolation(
                    pattern.element(),
                    format!("atom candidate mismatch for {}", candidate),
                ));
            }

            if let Err(error) = candidate.check_invariants() {
                return Err(ResolutionError::ValenceViolation(
                    pattern.element(),
                    format!(
                        "atom invariant verification failed for {}: {}",
                        candidate, error
                    ),
                ));
            }

            let atom = *candidate;
            let new_idx = graph.add_node(atom);
            index_map[old_idx.index()] = Some(new_idx);
        }

        for old_edge in self.graph.edge_indices() {
            let (a, b) = self.graph.edge_endpoints(old_edge).unwrap();
            let bond_builder = self.graph.edge_weight(old_edge).unwrap();
            let new_a = index_map[a.index()].unwrap();
            let new_b = index_map[b.index()].unwrap();
            graph.add_edge(
                new_a,
                new_b,
                bond_builder
                    .to_bond()
                    .map_err(|e| ResolutionError::InvalidBond(e.to_string()))?,
            );
        }

        let atom_charge: i8 = graph.node_weights().map(|a| a.charge()).sum();
        let bond_charge: i8 = graph.edge_weights().map(|b| b.charge()).sum();
        let aromatic_charge: i8 = self
            .aromatic_systems
            .iter()
            .map(|system| system.charge())
            .sum();
        let multicenter_charge: i8 = self
            .multicenter_bonds
            .iter()
            .flat_map(|bond| bond.sets().iter())
            .map(|set| set.charge())
            .sum();
        let charge = atom_charge + bond_charge + aromatic_charge + multicenter_charge;

        if let Some(explicit) = self.charge {
            if explicit != charge {
                return Err(ResolutionError::MolecularChargeMismatch {
                    explicit,
                    atom_sum: charge,
                });
            }
        }

        let mut feature_spins: Vec<SpinState> = graph.node_weights().map(|a| a.spin()).collect();
        feature_spins.extend(self.aromatic_systems.iter().map(|s| s.spin()));

        let spin = match self.spin {
            Some(explicit) => {
                if !explicit.is_constructible_from(&feature_spins) {
                    let atom_unpaired_sum: u16 = feature_spins
                        .iter()
                        .map(|s| s.unpaired_electrons() as u16)
                        .sum();
                    return Err(ResolutionError::MolecularSpinIncompatible {
                        explicit_unpaired: explicit.unpaired_electrons(),
                        explicit_multiplicity: explicit.multiplicity().multiplicity(),
                        atom_unpaired_sum,
                    });
                }
                explicit
            }
            None => {
                let atom_unpaired_sum: u16 = feature_spins
                    .iter()
                    .map(|s| s.unpaired_electrons() as u16)
                    .sum();
                let compatible =
                    compatible_molecular_multiplicities(&feature_spins).ok_or_else(|| {
                        let element = graph
                            .node_weights()
                            .next()
                            .map(|a| a.element())
                            .unwrap_or(umol_data::Element::C);
                        ResolutionError::ValenceViolation(
                            element,
                            "molecular spin exceeds maximum representable".to_string(),
                        )
                    })?;
                match compatible.as_slice() {
                    [single] => {
                        let multiplicity = SpinMultiplicity::from_multiplicity(*single)
                            .expect("compatible multiplicity is always in 1..=10");
                        SpinState::new(atom_unpaired_sum as u8, multiplicity)
                    }
                    [] => {
                        return Err(ResolutionError::ValenceViolation(
                            graph.node_weights().next().unwrap().element(),
                            "no compatible molecular spin for atom-level spins".to_string(),
                        ));
                    }
                    _ => {
                        return Err(ResolutionError::MolecularSpinIncomplete {
                            atom_unpaired_sum,
                            compatible_multiplicities: compatible,
                        });
                    }
                }
            }
        };

        Ok(Molecule::from_parts(
            graph,
            self.dative_bonds,
            self.aromatic_systems,
            self.multicenter_bonds,
            self.noncovalent_bonds,
            charge,
            spin,
        ))
    }
}

impl Default for MoleculeBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl FromAst<MoleculeAst> for MoleculeBuilder {
    fn from_ast(ast: &MoleculeAst, cfg: &MoleculeDslConfig) -> Result<Self, LoweringError> {
        let mut builder = Self::new();
        let mut indices: Vec<AtomIndex> = Vec::with_capacity(ast.atoms.len());

        for atom_ast in &ast.atoms {
            let pattern = AtomPattern::from_ast(atom_ast, &cfg.atom)?;
            indices.push(builder.add_atom(pattern));
        }

        let resolve = |i: usize| -> Result<AtomIndex, LoweringError> {
            indices
                .get(i)
                .copied()
                .ok_or_else(|| LoweringError::UnknownLabel(i.to_string()))
        };

        for bond in &ast.bonds {
            let a = resolve(bond.a)?;
            let b = resolve(bond.b)?;
            let pattern = BondPattern::from_ast(&bond.bond, &cfg.bond)?;
            builder.add_bond_unchecked(a, b, pattern);
        }

        for db in &ast.dative_bonds {
            let donor = resolve(db.donor)?;
            let acceptor = resolve(db.acceptor)?;
            let order = BondPattern::from_ast(&db.bond, &cfg.bond)?.order();
            builder.add_dative_bond(DativeBond::new(donor, acceptor, order));
        }

        for sys in &ast.aromatic_systems {
            let atoms: Vec<AtomIndex> = sys
                .atoms
                .iter()
                .map(|i| resolve(*i))
                .collect::<Result<_, _>>()?;
            let contributions: Vec<AromaticContribution> = atoms
                .into_iter()
                .map(|a| AromaticContribution::new(a, 0))
                .collect();
            builder.add_aromatic_system(AromaticSystem::new(contributions));
        }

        for mc in &ast.multicenter_bonds {
            let atoms: Vec<AtomIndex> = mc
                .atoms
                .iter()
                .map(|i| resolve(*i))
                .collect::<Result<_, _>>()?;
            let set = MulticenterSet::topology_only(atoms);
            builder.add_multicenter_bond(MulticenterBond::new(std::iter::once(set)));
        }

        for nc in &ast.noncovalent_bonds {
            let a = resolve(nc.a)?;
            let b = resolve(nc.b)?;
            builder.add_noncovalent_bond(NoncovalentBond::new(
                a,
                b,
                crate::bond::BondNoncovalent::Hydrogen,
            ));
        }

        if let Some(charge) = ast.charge {
            let c = i8::try_from(charge).map_err(|_| LoweringError::OutOfRange {
                field: "charge",
                value: charge,
            })?;
            builder.set_charge(c);
        }

        if let Some(spin) = ast.spin {
            builder.set_spin(spin);
        }

        Ok(builder)
    }
}

impl ToAst<MoleculeAst> for MoleculeBuilder {
    fn to_ast(&self, cfg: &MoleculeDslConfig) -> MoleculeAst {
        let atom_indices: Vec<AtomIndex> = self.atom_indices().collect();
        let position_of: HashMap<AtomIndex, usize> = atom_indices
            .iter()
            .enumerate()
            .map(|(i, &idx)| (idx, i))
            .collect();

        let atoms: Vec<_> = atom_indices
            .iter()
            .map(|&idx| self.atom(idx).unwrap().to_ast(&cfg.atom))
            .collect();

        let pos = |idx: AtomIndex| -> usize { *position_of.get(&idx).unwrap() };

        let bonds: Vec<LocalizedBondAst> = self
            .bond_indices()
            .map(|bi| {
                let (a, b) = self.bond_atom_indices(bi).unwrap();
                LocalizedBondAst {
                    a: pos(a),
                    b: pos(b),
                    bond: self.bond(bi).unwrap().to_ast(&cfg.bond),
                }
            })
            .collect();

        let dative_bonds: Vec<DativeBondAst> = self
            .dative_bonds()
            .map(|db| DativeBondAst {
                donor: pos(db.donor()),
                acceptor: pos(db.acceptor()),
                bond: BondAst::from_order(db.order()),
            })
            .collect();

        let aromatic_systems: Vec<AromaticSystemAst> = self
            .aromatic_systems()
            .map(|sys| AromaticSystemAst {
                atoms: sys.atoms().map(pos).collect(),
            })
            .collect();

        let multicenter_bonds: Vec<MulticenterBondAst> = self
            .multicenter_bonds()
            .map(|mc| MulticenterBondAst {
                atoms: mc.all_atoms().into_iter().map(pos).collect(),
            })
            .collect();

        let noncovalent_bonds: Vec<NoncovalentBondAst> = self
            .noncovalent_bonds()
            .map(|nc| NoncovalentBondAst {
                a: pos(nc.a()),
                b: pos(nc.b()),
                bond: BondAst::from_order(1),
            })
            .collect();

        MoleculeAst {
            atoms,
            bonds,
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            noncovalent_bonds,
            charge: self.charge().map(|c| c as i64),
            spin: self.spin(),
        }
    }
}

impl fmt::Display for MoleculeBuilder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        MoleculeAstWrapper::from_ast(self.to_ast(&MoleculeDslConfig::zeroed())).fmt(f)
    }
}

impl FromStr for MoleculeBuilder {
    type Err = LoweringError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let dsl = MoleculeAstWrapper::from_str(s).map_err(|e| LoweringError::Molecule(e.to_string()))?;
        Self::from_ast(dsl.ast(), &MoleculeDslConfig::zeroed())
    }
}

impl<'de> FromEdn<'de> for MoleculeBuilder {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let s = match edn {
            Edn::Str(s) => s.as_ref(),
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "string",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        let dsl = MoleculeAstWrapper::from_str(s).map_err(|e| DeError::subgrammar("molecule", e))?;
        Self::from_ast(dsl.ast(), &MoleculeDslConfig::zeroed())
            .map_err(|e| DeError::subgrammar("molecule", e))
    }
}

impl ToEdn for MoleculeBuilder {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(std::borrow::Cow::Owned(self.to_string()))
    }
}

fn compatible_molecular_multiplicities(states: &[SpinState]) -> Option<Vec<u8>> {
    let unpaired_total: u32 = states.iter().map(|s| s.unpaired_electrons() as u32).sum();
    if unpaired_total > u8::MAX as u32 {
        return None;
    }
    let total_u8 = unpaired_total as u8;
    let mut compatible = Vec::new();
    for m in 1..=10 {
        let Some(mult) = SpinMultiplicity::from_multiplicity(m) else {
            continue;
        };
        let Ok(candidate) = SpinState::try_new(total_u8, mult) else {
            continue;
        };
        if candidate.is_constructible_from(states) {
            compatible.push(m);
        }
    }
    Some(compatible)
}

#[cfg(test)]
mod tests;
