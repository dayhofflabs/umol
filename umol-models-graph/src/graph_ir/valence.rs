//! Valence resolution strategies and atom typing support for GraphIR.

use smallvec::SmallVec;
use umol_data::{SpinState, MAX_UNPAIRED_ELECTRONS};

use super::atom_type::{AtomTypeQuery, AtomTypeSpec};
use super::config_data::{AtomTypeRegistry, ValenceTable};
use super::molecule::{AtomIndex, MoleculeBuilder};

/// Atom-level valence validator.
#[derive(Debug, Clone)]
pub enum ValenceValidator {
    AtomTyping(AtomTypeRegistry),
    Counts(ValenceTable),
}

impl ValenceValidator {
    pub fn candidates_for(
        &self,
        builder: &MoleculeBuilder,
        atom_index: AtomIndex,
    ) -> SmallVec<[AtomTypeSpec; 4]> {
        match self {
            ValenceValidator::AtomTyping(registry) => {
                registry.candidates_for(&AtomTypeQuery::from_builder_atom(builder, atom_index))
            }
            ValenceValidator::Counts(table) => Self::counts_candidates(table, builder, atom_index),
        }
    }

    fn counts_candidates(
        table: &ValenceTable,
        builder: &MoleculeBuilder,
        atom_index: AtomIndex,
    ) -> SmallVec<[AtomTypeSpec; 4]> {
        let atom = builder.atom(atom_index).expect("atom_index must be valid");
        let element = atom.element();
        let charge = atom.charge().unwrap_or(0);
        let explicit_valence = builder.atom_bond_order_sum(atom_index);
        let (donated, accepted) = builder.atom_dative_bond_order_sums(atom_index);

        let entry = match table.entry(element) {
            Some(e) => e,
            None => return SmallVec::new(),
        };
        let implicit_hydrogens =
            match table.compute_implicit_hydrogens(element, charge, explicit_valence) {
                Some(h) => h,
                None => return SmallVec::new(),
            };
        let total_valence = explicit_valence + implicit_hydrogens;
        let num_electrons = (entry.outer_electrons as i16) - (charge as i16);
        let unassigned_electrons = num_electrons - (total_valence as i16);
        if unassigned_electrons < 0 {
            return SmallVec::new();
        }
        let unpaired_unassigned = (unassigned_electrons % 2) as u8;
        let lone_pairs_unassigned = (unassigned_electrons / 2) as u8;
        let unpaired = atom.unpaired_electrons().unwrap_or(unpaired_unassigned);
        if unpaired > MAX_UNPAIRED_ELECTRONS {
            return SmallVec::new();
        }
        let lone_pairs = atom.lone_pairs().unwrap_or(lone_pairs_unassigned);
        let spin = match atom.multiplicity() {
            Some(m) => match SpinState::try_new(unpaired, m) {
                Some(s) => s,
                None => return SmallVec::new(),
            },
            None => match SpinState::max_multiplicity(unpaired) {
                Some(s) => s,
                None => return SmallVec::new(),
            },
        };
        let multiplicity = spin.multiplicity();

        match AtomTypeSpec::new(
            element,
            charge,
            implicit_hydrogens,
            lone_pairs,
            unpaired,
            multiplicity,
            total_valence,
            donated,
            accepted,
            0,
            0,
        ) {
            Ok(spec) => SmallVec::from_elem(spec, 1),
            Err(_) => SmallVec::new(),
        }
    }
}
