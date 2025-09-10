//! SMILES parse state.

use std::collections::HashMap;

use smallvec::SmallVec;
use umol::error::DataError;
use umol::Result;

use crate::io::ir::{BondDir, BondOrder, Chirality};

/// Default depth of branching in SMILES.
/// Can be exceeded (incurs extra memory allocation) if needed.
const DEFAULT_BRANCH_DEPTH: usize = 4;

/// Default number of neighbors in a stereocenter.
/// Can be exceeded (incurs extra memory allocation) if needed.
const DEFAULT_NEIGHBOR_COUNT: usize = 4;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ParseState {
    // Lexer tracking
    pub byte_pos: usize,
    pub token_pos: usize,

    // Atom/bond tracking
    pub current_atom: usize,
    pub next_atom: usize,
    pub next_bond: usize,

    // Structure tracking
    pub rings: HashMap<u32, RingState>,
    pub branches: SmallVec<[BranchState; DEFAULT_BRANCH_DEPTH]>,
    pub stereocenters: HashMap<usize, StereocenterState>,

    // Bond context
    pub pending_bond: Option<Bond>,
}

impl ParseState {
    pub fn open_ring(&mut self, ring_index: u32, bond: Option<Bond>) -> Result<()> {
        if self.rings.contains_key(&ring_index) {
            return Err(DataError::InvalidRing(format!(
                "Ring index {} already exists",
                ring_index
            ))
            .into());
        }
        self.rings.insert(
            ring_index,
            RingState {
                start_atom: self.current_atom,
                pending_bond: bond,
            },
        );
        Ok(())
    }

    pub fn close_ring(&mut self, ring_index: u32) -> Result<(usize, Option<Bond>)> {
        self.rings
            .remove(&ring_index)
            .map(|state| (state.start_atom, state.pending_bond))
            .ok_or(
                DataError::InvalidRing(format!("Ring index {} does not exist", ring_index)).into(),
            )
    }

    pub fn open_branch(&mut self, bond: Option<Bond>) -> Result<()> {
        self.branches.push(BranchState {
            parent_atom: self.current_atom,
            return_bond: self.pending_bond.clone(),
        });
        self.pending_bond = bond;
        Ok(())
    }

    pub fn close_branch(&mut self) -> Result<()> {
        let state = self.branches.pop().ok_or_else(|| {
            DataError::InvalidFeature("Unmatched ')' with no open branch".to_string())
        })?;
        self.current_atom = state.parent_atom;
        self.pending_bond = state.return_bond;
        Ok(())
    }

    pub fn open_stereocenter(
        &mut self,
        atom: usize,
        chirality: Option<Chirality>,
        expected_neighbors: u8,
    ) -> Result<()> {
        if self.stereocenters.contains_key(&atom) {
            return Err(DataError::InvalidFeature(format!(
                "Stereocenter at atom {} already open",
                atom
            ))
            .into());
        }
        self.stereocenters.insert(
            atom,
            StereocenterState {
                chirality,
                neighbors: SmallVec::new(),
                expected_neighbors,
            },
        );
        Ok(())
    }

    pub fn add_substituent(&mut self, atom: usize, neighbor: usize) -> Result<()> {
        let sc = self.stereocenters.get_mut(&atom).ok_or_else(|| {
            DataError::InvalidFeature(format!(
                "No open stereocenter at atom {} to add substituent",
                atom
            ))
        })?;
        if sc.neighbors.contains(&neighbor) {
            return Err(DataError::InvalidFeature(format!(
                "Duplicate substituent {} at stereocenter {}",
                neighbor, atom
            ))
            .into());
        }
        if (sc.neighbors.len() as u8) >= sc.expected_neighbors {
            return Err(DataError::InvalidFeature(format!(
                "Too many substituents for stereocenter {}: expected {}",
                atom, sc.expected_neighbors
            ))
            .into());
        }
        sc.neighbors.push(neighbor);
        Ok(())
    }

    pub fn close_stereocenter(&mut self, atom: usize) -> Result<StereocenterState> {
        let sc = self.stereocenters.remove(&atom).ok_or_else(|| {
            DataError::InvalidFeature(format!(
                "Stereocenter at atom {} was not open",
                atom
            ))
        })?;
        if sc.expected_neighbors > 0 && sc.neighbors.len() as u8 != sc.expected_neighbors {
            return Err(DataError::InvalidFeature(format!(
                "Incomplete stereocenter {}: got {} of {} substituents",
                atom,
                sc.neighbors.len(),
                sc.expected_neighbors
            ))
            .into());
        }
        Ok(sc)
    }

    pub fn separate_component(&mut self) -> Result<()> {
        self.pending_bond = None;
        Ok(())
    }

    pub fn merge_bond_dir(&mut self, dir: BondDir) -> Result<()> {
        match &mut self.pending_bond {
            Some(b) => match b.dir {
                None => b.dir = Some(dir),
                Some(existing) if existing == dir => {},
                Some(_) => {
                    return Err(DataError::InvalidBondSpec(
                        "Conflicting bond directions on the same bond".to_string(),
                    )
                    .into())
                }
            },
            None => {
                self.pending_bond = Some(Bond { order: BondOrder::Single, dir: Some(dir) });
            }
        }
        Ok(())
    }

    pub fn next_bond_index(&mut self) -> usize {
        let idx = self.next_bond;
        self.next_bond += 1;
        idx
    }

    pub fn expect_neighbors(&mut self, atom: usize, n: u8) -> Result<()> {
        let sc = self.stereocenters.get_mut(&atom).ok_or_else(|| {
            DataError::InvalidFeature(format!(
                "No open stereocenter at atom {} to set expectation",
                atom
            ))
        })?;
        if sc.neighbors.len() as u8 > n {
            return Err(DataError::InvalidFeature(format!(
                "Too many neighbors already recorded for stereocenter {}: {} > {}",
                atom,
                sc.neighbors.len(),
                n
            ))
            .into());
        }
        sc.expected_neighbors = n;
        Ok(())
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct RingState {
    pub start_atom: usize,
    pub pending_bond: Option<Bond>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BranchState {
    pub parent_atom: usize,
    pub return_bond: Option<Bond>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StereocenterState {
    pub chirality: Option<Chirality>,
    pub neighbors: SmallVec<[usize; DEFAULT_NEIGHBOR_COUNT]>,
    pub expected_neighbors: u8,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Bond {
    pub order: BondOrder,
    pub dir: Option<BondDir>,
}


