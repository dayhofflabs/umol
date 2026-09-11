//! Candidate-bond components for local double-bond direction assignment.

use std::collections::BTreeMap;
use std::mem;
use std::ops::ControlFlow;

use smallvec::SmallVec;
use umol_graph_core::{BreadthFirstEvent, Graph, NodeId};

use crate::table_ir::{AtomNeighbors, Bond, BondConfiguration, BondOrder, Molecule};

/// One connected set of candidate bonds, containing at least one definite site.
#[derive(Debug, PartialEq, Eq)]
pub(super) struct MarkerComponent {
    pub(super) candidates: Vec<u32>,
    pub(super) sites: Vec<MarkerSite>,
}

/// Endpoint blocks use the bond row's first/second endpoint order.
#[derive(Debug, PartialEq, Eq)]
pub(super) struct MarkerSite {
    pub(super) bond: u32,
    pub(super) candidates: [SmallVec<[u32; 2]>; 2],
    pub(super) configuration: Option<BondConfiguration>,
}

#[derive(Debug, PartialEq, Eq)]
pub(super) enum MarkerComponentError {
    BondIndexOutOfBounds { bond: u32 },
    AtomIndexOutOfBounds { atom: u32 },
    DuplicateStereoBond { bond: u32 },
    UnsupportedSite { bond: u32 },
    MissingCandidate { bond: u32, atom: u32 },
}

/// Derives components before marker selection, preserving both candidate endpoint blocks.
///
/// Ordinary double bonds with one or two actual substituents per endpoint group their
/// ordinary single side bonds. Either and unasserted sites couple choices just like definite
/// sites, but only components containing definite assertions are returned. Sites without
/// candidates at both ends do not couple groups; definite ones report MissingCandidate.
/// No element, chemical-equivalence, ring-size, or direction-preference rule is applied.
///
/// Candidate, site, and component order follow bond-table indices. Frames are retained as
/// supplied; reference validation and parity transport belong to assignment. No markers,
/// molecular paths, ring openings, or preferred reference substituents are selected here.
/// This result belongs to the unchanged input table; it does not certify renderability.
pub(super) fn derive_marker_components(
    molecule: &Molecule,
) -> Result<Vec<MarkerComponent>, MarkerComponentError> {
    if !molecule
        .stereo_bonds
        .iter()
        .any(|frame| matches!(frame.configuration, BondConfiguration::Framed { .. }))
    {
        return Ok(Vec::new());
    }
    let mut configurations = BTreeMap::new();
    for frame in &molecule.stereo_bonds {
        let bond = molecule
            .bonds
            .get(frame.bond as usize)
            .ok_or(MarkerComponentError::BondIndexOutOfBounds { bond: frame.bond })?;
        if !ordinary(bond) || bond.order != BondOrder::Double {
            return Err(MarkerComponentError::UnsupportedSite { bond: frame.bond });
        }
        if configurations
            .insert(frame.bond, frame.configuration)
            .is_some()
        {
            return Err(MarkerComponentError::DuplicateStereoBond { bond: frame.bond });
        }
    }
    let neighbors = molecule.atom_neighbors();
    let mut sites = Vec::new();
    for (index, bond) in molecule.bonds.iter().enumerate() {
        if !ordinary(bond) || bond.order != BondOrder::Double {
            continue;
        }
        let index = index as u32;
        let endpoints = [bond.atoms.first(), bond.atoms.second()];
        for atom in endpoints {
            if atom as usize >= molecule.atoms.len() {
                return Err(MarkerComponentError::AtomIndexOutOfBounds { atom });
            }
        }
        let configuration = configurations.get(&index).copied();
        if local_substituents(molecule, &neighbors, index).is_none() {
            if matches!(configuration, Some(BondConfiguration::Framed { .. })) {
                return Err(MarkerComponentError::UnsupportedSite { bond: index });
            }
            continue;
        }
        let candidates: [SmallVec<[u32; 2]>; 2] = endpoints.map(|atom| {
            neighbors
                .neighbors(atom)
                .iter()
                .filter_map(|neighbor| {
                    let side_bond = &molecule.bonds[neighbor.bond as usize];
                    (ordinary(side_bond) && side_bond.order == BondOrder::Single)
                        .then_some(neighbor.bond)
                })
                .collect()
        });
        for (side, block) in candidates.iter().enumerate() {
            if block.is_empty() && matches!(configuration, Some(BondConfiguration::Framed { .. })) {
                return Err(MarkerComponentError::MissingCandidate {
                    bond: index,
                    atom: endpoints[side],
                });
            }
        }
        if candidates.iter().all(|side| !side.is_empty()) {
            sites.push(MarkerSite {
                bond: index,
                candidates,
                configuration,
            });
        }
    }

    let mut candidates: Vec<_> = sites
        .iter()
        .flat_map(|site| site.candidates.iter().flatten().copied())
        .collect();
    candidates.sort_unstable();
    candidates.dedup();
    let candidate_node = |bond| {
        candidates
            .binary_search(&bond)
            .expect("collected candidate") as u32
    };
    let mut edges = Vec::new();
    for site in &sites {
        let mut group = site
            .candidates
            .iter()
            .flatten()
            .copied()
            .map(candidate_node);
        if let Some(first) = group.next() {
            edges.extend(group.map(|other| [first, other]));
        }
    }
    let graph = Graph::new(candidates.len(), &edges);
    let roots = sites
        .iter()
        .filter(|site| matches!(site.configuration, Some(BondConfiguration::Framed { .. })))
        .map(|site| NodeId(candidate_node(site.candidates[0][0])));
    let mut membership = vec![None; candidates.len()];
    let mut components = Vec::new();
    let mut members = Vec::new();
    let _: ControlFlow<()> = graph.visit_breadth_first(roots, None, |event| {
        match event {
            BreadthFirstEvent::Discover { node, .. } => {
                membership[node.index()] = Some(components.len());
                members.push(candidates[node.index()]);
            }
            BreadthFirstEvent::FinishTree { .. } => {
                members.sort_unstable();
                components.push(MarkerComponent {
                    candidates: mem::take(&mut members),
                    sites: Vec::new(),
                });
            }
            BreadthFirstEvent::Finish { .. } => {}
        }
        ControlFlow::Continue(())
    });
    for site in sites {
        if let Some(component) = membership[candidate_node(site.candidates[0][0]) as usize] {
            components[component].sites.push(site);
        }
    }
    components.sort_unstable_by_key(|component| component.candidates[0]);
    Ok(components)
}

fn ordinary(bond: &Bond) -> bool {
    bond.donation.is_none() && bond.noncovalent.is_none()
}

fn local_substituents(
    molecule: &Molecule,
    neighbors: &AtomNeighbors,
    site: u32,
) -> Option<[SmallVec<[u32; 2]>; 2]> {
    let pair = molecule.bonds[site as usize].atoms;
    let endpoints = [pair.first(), pair.second()];
    if endpoints[0] == endpoints[1] {
        return None;
    }
    let mut substituents = [SmallVec::<[u32; 2]>::new(), SmallVec::new()];
    for (side, &atom) in endpoints.iter().enumerate() {
        for neighbor in neighbors.neighbors(atom) {
            let bond = &molecule.bonds[neighbor.bond as usize];
            if !ordinary(bond) || neighbor.bond == site {
                continue;
            }
            if bond.order == BondOrder::Double || neighbor.atom == atom {
                return None;
            }
            if neighbor.atom == endpoints[1 - side] || substituents[side].contains(&neighbor.atom) {
                continue;
            }
            if substituents[side].len() == 2 {
                return None;
            }
            substituents[side].push(neighbor.atom);
        }
        if substituents[side].is_empty() {
            return None;
        }
    }
    if substituents[0]
        .iter()
        .any(|atom| substituents[1].contains(atom))
    {
        return None;
    }
    Some(substituents)
}

#[cfg(test)]
mod tests;
