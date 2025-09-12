//! SMILES parse state.

use std::collections::HashMap;

use smallvec::SmallVec;
use slog::Logger;
use umol::error::DataError;
use umol::Result;

use crate::io::ir::{Atom as IRAtom, Bond as IRBond, BondDir, BondOrder, Chirality, Molecule as IRMolecule, SourceFormat};

/// Default depth of branching in SMILES.
/// Can be exceeded (incurs extra memory allocation) if needed.
const DEFAULT_BRANCH_DEPTH: usize = 4;

/// Default nesting depth for component stacks (dot-separated fragments inside branches).
const DEFAULT_COMPONENT_DEPTH: usize = 4;

/// Default number of neighbors in a stereocenter.
/// Can be exceeded (incurs extra memory allocation) if needed.
const DEFAULT_NEIGHBOR_COUNT: usize = 4;

#[derive(Clone, Debug, Default)]
pub struct ParseState {
    // Lexer tracking
    pub byte_pos: usize,
    pub tok_pos: usize,

    // Atom/bond tracking
    pub last_atom_idx: usize,
    pub next_atom_idx: usize,
    pub next_bond_idx: usize,

    // Structure tracking
    pub rings: HashMap<u32, RingState>,
    pub branches: SmallVec<[BranchState; DEFAULT_BRANCH_DEPTH]>,
    pub stereocenters: HashMap<usize, StereoState>,

    // Bond context
    pub staged_bond: Option<BondSpec>,

    // Logger
    pub log: Option<Logger>,

    // Deferred parse error (chain mode)
    pub first_err: Option<String>,

    // IR building (parser mode)
    pub buf_atoms: Vec<IRAtom>,
    pub buf_bonds: Vec<IRBond>,
    pub molecules: Vec<IRMolecule>,

    // Nested component molecules collected inside branches; appended after outer completes
    pub pending_molecules: Vec<IRMolecule>,

    // Component stack to support '.' inside branches
    pub comp_stack: SmallVec<[(Vec<IRAtom>, Vec<IRBond>); DEFAULT_COMPONENT_DEPTH]>,
}

impl ParseState {
    pub fn open_ring(&mut self, ring_index: u32, bond: Option<BondSpec>) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "open_ring"; "ring" => ring_index); }
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
                start_atom: self.last_atom_idx,
                pending_bond: bond,
            },
        );
        Ok(())
    }

    pub fn close_ring(&mut self, ring_index: u32) -> Result<(usize, Option<BondSpec>)> {
        if let Some(log) = &self.log { slog::debug!(log, "close_ring"; "ring" => ring_index); }
        self.rings
            .remove(&ring_index)
            .map(|state| (state.start_atom, state.pending_bond))
            .ok_or(
                DataError::InvalidRing(format!("Ring index {} does not exist", ring_index)).into(),
            )
    }

    pub fn open_branch(&mut self, bond: Option<BondSpec>) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "open_branch"); }
        self.branches.push(BranchState {
            parent_atom: self.last_atom_idx,
            return_bond: self.staged_bond.clone(),
        });
        self.staged_bond = bond;
        Ok(())
    }

    pub fn close_branch(&mut self) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "close_branch"); }
        let state = self.branches.pop().ok_or_else(|| {
            DataError::InvalidFeature("Unmatched ')' with no open branch".to_string())
        })?;
        self.last_atom_idx = state.parent_atom;
        self.staged_bond = state.return_bond;
        Ok(())
    }

    pub fn stereo_open(
        &mut self,
        atom: usize,
        chirality: Option<Chirality>,
        expected_neighbors: u8,
    ) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "open_stereocenter"; "atom" => atom as i64); }
        if self.stereocenters.contains_key(&atom) {
            return Err(DataError::InvalidFeature(format!(
                "Stereocenter at atom {} already open",
                atom
            ))
            .into());
        }
        self.stereocenters.insert(
            atom,
            StereoState {
                chirality,
                neighbors: SmallVec::new(),
                expected_neighbors,
            },
        );
        Ok(())
    }

    pub fn stereo_add(&mut self, atom: usize, neighbor: usize) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "add_substituent"; "atom" => atom as i64, "neighbor" => neighbor as i64); }
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

    pub fn stereo_close(&mut self, atom: usize) -> Result<StereoState> {
        if let Some(log) = &self.log { slog::debug!(log, "close_stereocenter"; "atom" => atom as i64); }
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

    pub fn split_component(&mut self) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "split_component"); }
        self.staged_bond = None;
        Ok(())
    }

    /// Begin parsing a nested component (used for '.' inside branches)
    pub fn enter_comp(&mut self) {
        if let Some(log) = &self.log { slog::debug!(log, "enter_comp"); }
        let saved = (std::mem::take(&mut self.buf_atoms), std::mem::take(&mut self.buf_bonds));
        self.comp_stack.push(saved);
    }

    /// Finalize the current component into a molecule and restore the previous context
    pub fn exit_comp(&mut self) {
        if let Some(log) = &self.log { slog::debug!(log, "exit_comp"); }
        // finalize the subcomponent we're in now
        self.finish_molecule();
        if let Some(m) = self.molecules.pop() {
            self.pending_molecules.push(m);
        }
        if let Some((atoms, bonds)) = self.comp_stack.pop() {
            self.buf_atoms = atoms;
            self.buf_bonds = bonds;
        }
    }

    pub fn stage_bond_dir(&mut self, dir: BondDir) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "stage_bond_dir"); }
        match &mut self.staged_bond {
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
                self.staged_bond = Some(BondSpec { order: BondOrder::Single, dir: Some(dir) });
            }
        }
        Ok(())
    }

    pub fn stage_bond_order(&mut self, order: BondOrder) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "stage_bond_order"; "order" => format!("{:?}", order)); }
        if self.last_atom_idx == 0 && self.next_atom_idx == 0 {
            // Defer error reporting in chain mode; record first error
            if self.first_err.is_none() { self.first_err = Some("Bond cannot start a chain".to_string()); }
            return Err(DataError::InvalidBondSpec("Bond cannot start a chain".to_string()).into());
        }
        if self.staged_bond.is_some() {
            if self.first_err.is_none() { self.first_err = Some("Consecutive bond symbols are not allowed".to_string()); }
            return Err(DataError::InvalidBondSpec("Consecutive bond symbols are not allowed".to_string()).into());
        }
        self.staged_bond = Some(BondSpec { order, dir: None });
        Ok(())
    }

    pub fn bump_bond_idx(&mut self) -> usize {
        if let Some(log) = &self.log { slog::debug!(log, "bump_bond_idx"; "before" => self.next_bond_idx as i64); }
        let idx = self.next_bond_idx;
        self.next_bond_idx += 1;
        idx
    }

    pub fn bump_atom_idx(&mut self) -> usize {
        if let Some(log) = &self.log { slog::debug!(log, "bump_atom_idx"; "before" => self.next_atom_idx as i64); }
        let idx = self.next_atom_idx;
        self.next_atom_idx += 1;
        idx
    }

    pub fn stereo_expect(&mut self, atom: usize, n: u8) -> Result<()> {
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

    pub fn link_to(&mut self, new_atom: usize) -> Result<Option<ResolvedBond>> {
        if let Some(log) = &self.log { slog::debug!(log, "link_to"; "from" => self.last_atom_idx as i64, "to" => new_atom as i64); }
        if new_atom == 0 { return Ok(None); }
        let start_atom = self.last_atom_idx;
        // Resolve bond to use (explicit staged or implicit single)
        let bond = self
            .staged_bond
            .take()
            .unwrap_or(BondSpec { order: BondOrder::Single, dir: None });
        // Reserve bond index
        let bond_index = self.bump_bond_idx();
        // Stereo neighbor tracking (best-effort; ignore errors here)
        if self.stereocenters.contains_key(&start_atom) {
            let _ = self.stereo_add(start_atom, new_atom);
        }
        if self.stereocenters.contains_key(&new_atom) {
            let _ = self.stereo_add(new_atom, start_atom);
        }
        // Trailing bond guard note remains
        Ok(Some(ResolvedBond { start_atom, end_atom: new_atom, bond_index, bond }))
    }

    pub fn finish_chain(&mut self) -> Result<()> {
        if let Some(ref msg) = self.first_err {
            return Err(DataError::InvalidBondSpec(msg.clone()).into());
        }
        if self.staged_bond.is_some() {
            return Err(DataError::InvalidBondSpec("Trailing bond symbol".to_string()).into());
        }
        Ok(())
    }

    pub fn push_atom(&mut self, atom: IRAtom) {
        self.buf_atoms.push(atom);
    }

    pub fn push_resolved_bond(&mut self, rb: ResolvedBond) {
        let mut bond = IRBond::from_order(rb.bond.order);
        bond.start_atom = Some(rb.start_atom as u32);
        bond.end_atom = Some(rb.end_atom as u32);
        bond.direction = rb.bond.dir;
        self.buf_bonds.push(bond);
    }

    pub fn finish_molecule(&mut self) {
        if self.buf_atoms.is_empty() {
            return;
        }
        let mut mol = IRMolecule::default();
        mol.source_format = SourceFormat::SMILES;
        mol.atoms = std::mem::take(&mut self.buf_atoms);
        mol.bonds = std::mem::take(&mut self.buf_bonds);
        self.molecules.push(mol);
    }

    pub fn drain_molecules(&mut self) -> Vec<IRMolecule> {
        std::mem::take(&mut self.molecules)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBond {
    pub start_atom: usize,
    pub end_atom: usize,
    pub bond_index: usize,
    pub bond: BondSpec,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct RingState {
    pub start_atom: usize,
    pub pending_bond: Option<BondSpec>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BranchState {
    pub parent_atom: usize,
    pub return_bond: Option<BondSpec>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StereoState {
    pub chirality: Option<Chirality>,
    pub neighbors: SmallVec<[usize; DEFAULT_NEIGHBOR_COUNT]>,
    pub expected_neighbors: u8,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BondSpec {
    pub order: BondOrder,
    pub dir: Option<BondDir>,
}


