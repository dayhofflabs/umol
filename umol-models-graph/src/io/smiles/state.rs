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
    pub staged_bond: Option<BondInfo>,

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
    pub fn open_ring(&mut self, ring_index: u32, bond: Option<BondInfo>) -> Result<()> {
        if self.rings.contains_key(&ring_index) {
            return Err(DataError::InvalidRing(format!(
                "Ring index {} already exists",
                ring_index
            ))
            .into());
        }
        let start_atom = self.last_atom_idx;
        self.rings.insert(
            ring_index,
            RingState {
                start_atom,
                pending_bond: bond,
            },
        );
        if let Some(log) = &self.log {
            let (ord, dir) = match &self.rings.get(&ring_index).unwrap().pending_bond {
                Some(b) => (format!("{:?}", b.order), format!("{:?}", b.dir)),
                None => ("<implicit>".to_string(), "None".to_string()),
            };
            slog::debug!(log, "ring_opened";
                "ring" => ring_index as i64,
                "start_atom" => start_atom as i64,
                "pending_order" => ord,
                "pending_dir" => dir,
            );
        }
        Ok(())
    }

    pub fn close_ring(&mut self, ring_index: u32) -> Result<(usize, Option<BondInfo>)> {
        let removed = self.rings
            .remove(&ring_index)
            .map(|state| (state.start_atom, state.pending_bond));
        match removed {
            Some((start_atom, pending)) => {
                if let Some(log) = &self.log {
                    let (ord, dir) = match &pending {
                        Some(b) => (format!("{:?}", b.order), format!("{:?}", b.dir)),
                        None => ("<implicit>".to_string(), "None".to_string()),
                    };
                    slog::debug!(log, "ring_closed";
                        "ring" => ring_index as i64,
                        "start_atom" => start_atom as i64,
                        "end_atom" => self.last_atom_idx as i64,
                        "pending_order" => ord,
                        "pending_dir" => dir,
                    );
                }
                Ok((start_atom, pending))
            }
            None => Err(
                DataError::InvalidRing(format!("Ring index {} does not exist", ring_index)).into(),
            ),
        }
    }

    pub fn open_branch(&mut self, bond: Option<BondInfo>) -> Result<()> {
        let parent_atom = self.last_atom_idx;
        let return_bond = self.staged_bond.clone();
        self.branches.push(BranchState { parent_atom, return_bond });
        self.staged_bond = bond.clone();
        if let Some(log) = &self.log {
            let rb = match &self.branches.last().unwrap().return_bond {
                Some(b) => format!("{:?}/{:?}", b.order, b.dir),
                None => "None".to_string(),
            };
            let nb = match &bond {
                Some(b) => format!("{:?}/{:?}", b.order, b.dir),
                None => "None".to_string(),
            };
            slog::debug!(log, "branch_opened";
                "parent_atom" => parent_atom as i64,
                "stack_depth" => self.branches.len() as i64,
                "return_bond" => rb,
                "branch_bond" => nb,
            );
        }
        Ok(())
    }

    pub fn close_branch(&mut self) -> Result<()> {
        let state = self.branches.pop().ok_or_else(|| {
            DataError::InvalidFeature("Unmatched ')' with no open branch".to_string())
        })?;
        self.last_atom_idx = state.parent_atom;
        self.staged_bond = state.return_bond;
        if let Some(log) = &self.log {
            let sb = match &self.staged_bond {
                Some(b) => format!("{:?}/{:?}", b.order, b.dir),
                None => "None".to_string(),
            };
            slog::debug!(log, "branch_closed";
                "restored_parent_atom" => self.last_atom_idx as i64,
                "stack_depth" => self.branches.len() as i64,
                "restored_bond" => sb,
            );
        }
        Ok(())
    }

    pub fn stereo_open(
        &mut self,
        atom: usize,
        chirality: Option<Chirality>,
        expected_neighbors: u8,
    ) -> Result<()> {
        if let Some(log) = &self.log {
            slog::debug!(log, "stereo_open";
                "atom" => atom as i64,
                "chirality" => format!("{:?}", chirality),
                "expected_neighbors" => expected_neighbors as i64,
            );
        }
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
        if let Some(log) = &self.log { slog::debug!(log, "stereo_add_neighbor"; "atom" => atom as i64, "neighbor" => neighbor as i64); }
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
        if let Some(log) = &self.log { slog::debug!(log, "stereo_close"; "atom" => atom as i64); }
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
        if let Some(log) = &self.log { slog::debug!(log, "component_split_by_dot"; "atom_at_split" => self.last_atom_idx as i64); }
        self.staged_bond = None;
        Ok(())
    }

    /// Begin parsing a nested component (used for '.' inside branches)
    pub fn enter_comp(&mut self) {
        if let Some(log) = &self.log { slog::debug!(log, "enter_component_context"; "depth_before" => self.comp_stack.len() as i64, "atoms" => self.buf_atoms.len() as i64, "bonds" => self.buf_bonds.len() as i64); }
        let saved = (std::mem::take(&mut self.buf_atoms), std::mem::take(&mut self.buf_bonds));
        self.comp_stack.push(saved);
    }

    /// Finalize the current component into a molecule and restore the previous context
    pub fn exit_comp(&mut self) {
        // finalize the subcomponent we're in now
        self.finish_molecule();
        if let Some(m) = self.molecules.pop() {
            self.pending_molecules.push(m);
        }
        if let Some((atoms, bonds)) = self.comp_stack.pop() {
            self.buf_atoms = atoms;
            self.buf_bonds = bonds;
        }
        if let Some(log) = &self.log { slog::debug!(log, "exit_component_context"; "depth_after" => self.comp_stack.len() as i64, "pending_molecules" => self.pending_molecules.len() as i64); }
    }

    pub fn stage_bond_dir(&mut self, dir: BondDir) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "stage_bond_dir"; "new_dir" => format!("{:?}", dir)); }
        match &mut self.staged_bond {
            Some(b) => match b.dir {
                None => b.dir = Some(dir),
                Some(existing) if existing == dir => {},
                Some(_) => {
                    return Err(DataError::InvalidBond(
                        "Conflicting bond directions on the same bond".to_string(),
                    )
                    .into())
                }
            },
            None => {
                self.staged_bond = Some(BondInfo { order: BondOrder::Single, dir: Some(dir) });
            }
        }
        Ok(())
    }

    pub fn stage_bond_order(&mut self, order: BondOrder) -> Result<()> {
        if let Some(log) = &self.log { slog::debug!(log, "stage_bond_order"; "order" => format!("{:?}", order), "has_staged" => self.staged_bond.is_some()); }
        if self.last_atom_idx == 0 && self.next_atom_idx == 0 {
            // Defer error reporting in chain mode; record first error
            if self.first_err.is_none() { self.first_err = Some("Bond cannot start a chain".to_string()); }
            return Err(DataError::InvalidBond("Bond cannot start a chain".to_string()).into());
        }
        if self.staged_bond.is_some() {
            if self.first_err.is_none() { self.first_err = Some("Consecutive bond symbols are not allowed".to_string()); }
            return Err(DataError::InvalidBond("Consecutive bond symbols are not allowed".to_string()).into());
        }
        self.staged_bond = Some(BondInfo { order, dir: None });
        Ok(())
    }

    pub fn bump_bond_idx(&mut self) -> usize {
        let idx = self.next_bond_idx;
        self.next_bond_idx += 1;
        if let Some(log) = &self.log { slog::debug!(log, "alloc_bond_index"; "allocated" => idx as i64, "next" => self.next_bond_idx as i64); }
        idx
    }

    pub fn bump_atom_idx(&mut self) -> usize {
        let idx = self.next_atom_idx;
        self.next_atom_idx += 1;
        if let Some(log) = &self.log { slog::debug!(log, "alloc_atom_index"; "allocated" => idx as i64, "next" => self.next_atom_idx as i64, "link_from" => self.last_atom_idx as i64); }
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
        if new_atom == 0 { return Ok(None); }
        let start_atom = self.last_atom_idx;
        // Resolve bond to use (explicit staged or implicit single)
        let bond = self
            .staged_bond
            .take()
            .unwrap_or(BondInfo { order: BondOrder::Single, dir: None });
        // Reserve bond index
        let bond_index = self.bump_bond_idx();
        if let Some(log) = &self.log {
            slog::debug!(log, "create_bond";
                "start_atom" => start_atom as i64,
                "end_atom" => new_atom as i64,
                "order" => format!("{:?}", bond.order),
                "dir" => format!("{:?}", bond.dir),
                "bond_index" => bond_index as i64,
            );
        }
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
            if let Some(log) = &self.log { slog::debug!(log, "finish_chain_error"; "error" => msg.clone()); }
            return Err(DataError::InvalidBond(msg.clone()).into());
        }
        if self.staged_bond.is_some() {
            if let Some(log) = &self.log { slog::debug!(log, "finish_chain_error"; "error" => "Trailing bond symbol"); }
            return Err(DataError::InvalidBond("Trailing bond symbol".to_string()).into());
        }
        if let Some(log) = &self.log { slog::debug!(log, "finish_chain_ok"; "atoms" => self.next_atom_idx as i64, "bonds" => self.next_bond_idx as i64); }
        Ok(())
    }

    pub fn push_atom(&mut self, atom: IRAtom) {
        let idx = atom.index.unwrap_or_default();
        if let Some(log) = &self.log { slog::debug!(log, "push_atom"; "atom_index" => idx as i64); }
        self.buf_atoms.push(atom);
    }

    pub fn push_resolved_bond(&mut self, rb: ResolvedBond) {
        if let Some(log) = &self.log {
            slog::debug!(log, "push_bond";
                "start_atom" => rb.start_atom as i64,
                "end_atom" => rb.end_atom as i64,
                "order" => format!("{:?}", rb.bond.order),
                "dir" => format!("{:?}", rb.bond.dir),
                "bond_index" => rb.bond_index as i64,
            );
        }
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
        if let Some(log) = &self.log {
            slog::debug!(log, "finish_molecule";
                "atoms" => mol.atoms.len() as i64,
                "bonds" => mol.bonds.len() as i64,
            );
        }
        self.molecules.push(mol);
    }

    pub fn drain_molecules(&mut self) -> Vec<IRMolecule> {
        let out = std::mem::take(&mut self.molecules);
        if let Some(log) = &self.log { slog::debug!(log, "drain_molecules"; "count" => out.len() as i64); }
        out
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBond {
    pub start_atom: usize,
    pub end_atom: usize,
    pub bond_index: usize,
    pub bond: BondInfo,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct RingState {
    pub start_atom: usize,
    pub pending_bond: Option<BondInfo>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BranchState {
    pub parent_atom: usize,
    pub return_bond: Option<BondInfo>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StereoState {
    pub chirality: Option<Chirality>,
    pub neighbors: SmallVec<[usize; DEFAULT_NEIGHBOR_COUNT]>,
    pub expected_neighbors: u8,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BondInfo {
    pub order: BondOrder,
    pub dir: Option<BondDir>,
}


