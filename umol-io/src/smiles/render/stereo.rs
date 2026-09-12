//! Candidate selection and direction assignment for local double-bond stereo.

use std::cell::OnceCell;
use std::collections::BTreeMap;
use std::mem;
use std::ops::ControlFlow;

use smallvec::SmallVec;
use umol_graph_core::{BreadthFirstEvent, Graph, NodeId};

use super::traversal::Traversal;
use crate::table_ir::{
    AtomNeighbors, Bond, BondConfiguration, BondDirection, BondOrder, BondRelation, Molecule,
};

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
/// Component candidates, sites, and components follow bond-table indices.
/// Per-side candidates follow neighbor atom order. Frames are retained as
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

#[derive(Debug, PartialEq, Eq)]
pub(super) enum MarkerAssignmentError {
    Component(MarkerComponentError),
    InvalidReference {
        bond: u32,
        atom: u32,
    },
    /// Lowest-index definite site in the unsatisfiable component.
    NoAssignment {
        bond: u32,
    },
}

/// Returns selected bonds in table order, with directions from each bond's first endpoint.
/// The traversal belongs to this unchanged table and supplies each group's initial `/`.
/// Either sites constrain coverage here; notation support for Either belongs to rendering.
pub(super) fn assign_markers(
    molecule: &Molecule,
    traversal: &Traversal,
) -> Result<Vec<(u32, BondDirection)>, MarkerAssignmentError> {
    let components =
        derive_marker_components(molecule).map_err(MarkerAssignmentError::Component)?;
    if components.is_empty() {
        return Ok(Vec::new());
    }
    let mut reversed: BTreeMap<_, _> = components
        .iter()
        .flat_map(|component| component.candidates.iter().map(|&bond| (bond, false)))
        .collect();
    for (bond, endpoints) in traversal.bonds() {
        if let Some(reverse) = reversed.get_mut(&bond) {
            *reverse = molecule.bonds[bond as usize].atoms.first() != endpoints[0];
        }
    }
    let neighbors = OnceCell::new();
    let mut output = Vec::new();
    for component in components {
        let mut sites = Vec::with_capacity(component.sites.len());
        let mut constraints =
            vec![SmallVec::<[(usize, bool); 4]>::new(); component.candidates.len()];
        for site in &component.sites {
            let pair = molecule.bonds[site.bond as usize].atoms;
            let endpoints = [pair.first(), pair.second()];
            let candidates = site.candidates.each_ref().map(|side| {
                side.iter()
                    .map(|bond| {
                        component
                            .candidates
                            .binary_search(bond)
                            .expect("component candidate")
                    })
                    .collect::<SmallVec<[usize; 2]>>()
            });
            let mut connect = |first: usize, second: usize, opposite: bool| {
                constraints[first].push((second, opposite));
                constraints[second].push((first, opposite));
            };
            let outgoing_reversed = |candidate: usize, endpoint| {
                molecule.bonds[component.candidates[candidate] as usize]
                    .atoms
                    .first()
                    != endpoint
            };
            for (side, block) in candidates.iter().enumerate() {
                if let &[first, second] = block.as_slice() {
                    connect(
                        first,
                        second,
                        true ^ outgoing_reversed(first, endpoints[side])
                            ^ outgoing_reversed(second, endpoints[side]),
                    );
                }
            }
            let definite = matches!(site.configuration, Some(BondConfiguration::Framed { .. }));
            if let Some(BondConfiguration::Framed {
                references,
                relation,
            }) = site.configuration
            {
                for side in 0..2 {
                    let reference = references[side];
                    let is_candidate = site.candidates[side].iter().any(|&bond| {
                        let pair = molecule.bonds[bond as usize].atoms;
                        pair.first() == reference || pair.second() == reference
                    });
                    if reference == endpoints[side]
                        || reference == endpoints[1 - side]
                        || (!is_candidate
                            && !neighbors
                                .get_or_init(|| molecule.atom_neighbors())
                                .neighbors(endpoints[side])
                                .iter()
                                .any(|neighbor| {
                                    neighbor.atom == reference
                                        && ordinary(&molecule.bonds[neighbor.bond as usize])
                                }))
                    {
                        return Err(MarkerAssignmentError::InvalidReference {
                            bond: site.bond,
                            atom: reference,
                        });
                    }
                }
                let reference_reversed = |candidate: usize, side: usize| {
                    let pair = molecule.bonds[component.candidates[candidate] as usize].atoms;
                    let other = if pair.first() == endpoints[side] {
                        pair.second()
                    } else {
                        pair.first()
                    };
                    outgoing_reversed(candidate, endpoints[side]) ^ (other != references[side])
                };
                for &first in &candidates[0] {
                    for &second in &candidates[1] {
                        connect(
                            first,
                            second,
                            (relation == BondRelation::OppositeSide)
                                ^ reference_reversed(first, 0)
                                ^ reference_reversed(second, 1),
                        );
                    }
                }
            }
            sites.push(SelectionSite {
                candidates,
                definite,
            });
        }
        let seeds: Vec<_> = component
            .candidates
            .iter()
            .map(|bond| reversed[bond])
            .collect();
        let Some(directions) = select_markers(&sites, &constraints, &seeds) else {
            return Err(MarkerAssignmentError::NoAssignment {
                bond: component
                    .sites
                    .iter()
                    .find(|site| {
                        matches!(site.configuration, Some(BondConfiguration::Framed { .. }))
                    })
                    .expect("definite component")
                    .bond,
            });
        };
        output.extend(component.candidates.into_iter().zip(directions).filter_map(
            |(bond, direction)| {
                direction.map(|falling| {
                    (
                        bond,
                        if falling {
                            BondDirection::Falling
                        } else {
                            BondDirection::Rising
                        },
                    )
                })
            },
        ));
    }
    output.sort_unstable_by_key(|&(bond, _)| bond);
    Ok(output)
}

struct SelectionSite {
    candidates: [SmallVec<[usize; 2]>; 2],
    definite: bool,
}

fn select_markers(
    sites: &[SelectionSite],
    constraints: &[SmallVec<[(usize, bool); 4]>],
    seeds: &[bool],
) -> Option<Vec<Option<bool>>> {
    let mut selected = vec![None; seeds.len()];
    let mut trail = Vec::with_capacity(seeds.len());
    let mut decisions = Vec::new();
    let mut directions = vec![None; seeds.len()];
    let mut pending = Vec::new();
    loop {
        if propagate_selection(sites, &mut selected, &mut trail) {
            if let Some(candidate) = selected.iter().position(Option::is_none) {
                decisions.push((candidate, trail.len(), false));
                selected[candidate] = Some(false);
                trail.push(candidate);
                continue;
            }
            if assign_directions(&selected, constraints, seeds, &mut directions, &mut pending) {
                return Some(directions);
            }
        }
        loop {
            let (candidate, checkpoint, tried_marked) = decisions.pop()?;
            for index in trail.drain(checkpoint..) {
                selected[index] = None;
            }
            if !tried_marked {
                decisions.push((candidate, checkpoint, true));
                selected[candidate] = Some(true);
                trail.push(candidate);
                break;
            }
        }
    }
}

fn propagate_selection(
    sites: &[SelectionSite],
    selected: &mut [Option<bool>],
    trail: &mut Vec<usize>,
) -> bool {
    loop {
        let start = trail.len();
        for site in sites {
            let covered = site.candidates.each_ref().map(|block| {
                block
                    .iter()
                    .any(|&candidate| selected[candidate] == Some(true))
            });
            for side in 0..2 {
                if site.definite && !covered[side] {
                    let mut remaining = site.candidates[side]
                        .iter()
                        .copied()
                        .filter(|&candidate| selected[candidate].is_none());
                    let Some(candidate) = remaining.next() else {
                        return false;
                    };
                    if remaining.next().is_none() {
                        selected[candidate] = Some(true);
                        trail.push(candidate);
                    }
                } else if !site.definite && covered[1 - side] {
                    for &candidate in &site.candidates[side] {
                        match selected[candidate] {
                            Some(true) => return false,
                            Some(false) => {}
                            None => {
                                selected[candidate] = Some(false);
                                trail.push(candidate);
                            }
                        }
                    }
                }
            }
        }
        if trail.len() == start {
            return true;
        }
    }
}

fn assign_directions(
    selected: &[Option<bool>],
    constraints: &[SmallVec<[(usize, bool); 4]>],
    seeds: &[bool],
    directions: &mut [Option<bool>],
    pending: &mut Vec<usize>,
) -> bool {
    directions.fill(None);
    pending.clear();
    for root in 0..selected.len() {
        if selected[root] != Some(true) || directions[root].is_some() {
            continue;
        }
        directions[root] = Some(seeds[root]);
        pending.push(root);
        while let Some(first) = pending.pop() {
            let direction = directions[first].expect("assigned direction");
            for &(second, opposite) in &constraints[first] {
                if selected[second] != Some(true) {
                    continue;
                }
                let expected = direction ^ opposite;
                match directions[second] {
                    Some(actual) if actual != expected => return false,
                    Some(_) => {}
                    None => {
                        directions[second] = Some(expected);
                        pending.push(second);
                    }
                }
            }
        }
    }
    true
}

#[cfg(test)]
mod tests;
