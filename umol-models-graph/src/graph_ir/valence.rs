//! Valence resolution strategies and atom typing support for GraphIR.

use smallvec::SmallVec;
use umol_data::{Element, SpinState, MAX_UNPAIRED_ELECTRONS};

use crate::graph_ir::atom::AtomBuilder;
use crate::graph_ir::atom_type::{AromaticValence, AtomTypeQuery, AtomTypeSpec};
use crate::graph_ir::config::ValenceStrategy;
use crate::graph_ir::config_data::{AtomTypeRegistry, ValenceEntry, ValenceTable};
use crate::graph_ir::molecule::{AtomIndex, MoleculeBuilder};

pub enum ValenceMatcher {
    AtomTyping {
        registry: AtomTypeRegistry,
    },
    Counts {
        table: ValenceTable,
        allow_implicit_hydrogens: bool,
    },
}

impl ValenceMatcher {
    pub fn new(strategy: &ValenceStrategy) -> Self {
        match strategy {
            ValenceStrategy::AtomTyping { registry } => Self::AtomTyping {
                registry: registry.clone(),
            },
            ValenceStrategy::Counts {
                table,
                allow_implicit_hydrogens,
            } => Self::Counts {
                table: table.clone(),
                allow_implicit_hydrogens: *allow_implicit_hydrogens,
            },
        }
    }

    pub fn candidates_for(
        &self,
        builder: &MoleculeBuilder,
        atom_index: AtomIndex,
    ) -> SmallVec<[AtomTypeSpec; 4]> {
        match self {
            Self::AtomTyping { registry } => {
                registry.candidates_for(&AtomTypeQuery::from_builder_atom(builder, atom_index))
            }
            Self::Counts {
                table,
                allow_implicit_hydrogens,
            } => counts_candidates(table, *allow_implicit_hydrogens, builder, atom_index),
        }
    }
}

fn counts_candidates(
    table: &ValenceTable,
    allow_implicit_hydrogens: bool,
    builder: &MoleculeBuilder,
    atom_index: AtomIndex,
) -> SmallVec<[AtomTypeSpec; 4]> {
    let atom = builder.atom(atom_index).expect("atom_index must be valid");
    let element = atom.element();
    let charge = atom.charge().unwrap_or(0);
    let explicit_valence = builder.atom_bond_order_sum(atom_index);
    let (donated_pairs, accepted_pairs) = builder.atom_dative_bond_order_sums(atom_index);

    let entry = match table.entry(element) {
        Some(e) => e,
        None => return SmallVec::new(),
    };

    if builder.atom_aromatic_hint(atom_index) {
        return build_aromatic_spec(
            entry,
            element,
            charge,
            explicit_valence,
            donated_pairs,
            accepted_pairs,
            &atom,
        );
    }

    let implicit_hydrogens = if let Some(h) = atom.hydrogens() {
        h
    } else if allow_implicit_hydrogens {
        match table.compute_implicit_hydrogens(element, charge, explicit_valence) {
            Some(h) => h,
            None => return SmallVec::new(),
        }
    } else {
        0
    };
    build_spec(
        element,
        charge,
        implicit_hydrogens,
        explicit_valence,
        donated_pairs,
        accepted_pairs,
        AromaticValence::None,
        &atom,
        entry,
    )
}

fn build_aromatic_spec(
    entry: &ValenceEntry,
    element: Element,
    charge: i8,
    explicit_valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    atom: &AtomBuilder,
) -> SmallVec<[AtomTypeSpec; 4]> {
    if entry.allowed_aromatic_valences.is_empty() {
        return SmallVec::new();
    }

    let effective_electrons = (entry.outer_electrons as i16) - (charge as i16);
    let mut candidates = SmallVec::new();

    for &a in &entry.allowed_aromatic_valences {
        let sigma_budget = effective_electrons - (a as i16);
        if sigma_budget < explicit_valence as i16 {
            continue;
        }
        let implicit_h = if let Some(h) = atom.hydrogens() {
            h
        } else {
            (sigma_budget - explicit_valence as i16) as u8
        };
        let total_sigma = explicit_valence + implicit_h;
        let remaining = effective_electrons - total_sigma as i16 - a as i16;
        if remaining < 0 || remaining % 2 != 0 {
            continue;
        }
        if let Some(spec) = try_build_spec(
            element,
            charge,
            implicit_h,
            explicit_valence,
            donated_pairs,
            accepted_pairs,
            AromaticValence::Valence(a),
            atom,
            entry,
        ) {
            candidates.push(spec);
        }
    }

    candidates
}

fn build_spec(
    element: Element,
    charge: i8,
    implicit_hydrogens: u8,
    explicit_valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    aromatic_valence: AromaticValence,
    atom: &AtomBuilder,
    entry: &ValenceEntry,
) -> SmallVec<[AtomTypeSpec; 4]> {
    match try_build_spec(
        element,
        charge,
        implicit_hydrogens,
        explicit_valence,
        donated_pairs,
        accepted_pairs,
        aromatic_valence,
        atom,
        entry,
    ) {
        Some(spec) => SmallVec::from_elem(spec, 1),
        None => SmallVec::new(),
    }
}

fn try_build_spec(
    element: Element,
    charge: i8,
    implicit_hydrogens: u8,
    explicit_valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    aromatic_valence: AromaticValence,
    atom: &AtomBuilder,
    entry: &ValenceEntry,
) -> Option<AtomTypeSpec> {
    let total_valence = explicit_valence + implicit_hydrogens;
    let num_electrons = (entry.outer_electrons as i16) - (charge as i16);
    let unassigned = num_electrons - (total_valence as i16) - (aromatic_valence.valence() as i16);
    if unassigned < 0 {
        return None;
    }
    let unpaired_unassigned = (unassigned % 2) as u8;
    let lone_pairs_unassigned = (unassigned / 2) as u8;
    let unpaired = atom.unpaired_electrons().unwrap_or(unpaired_unassigned);
    if unpaired > MAX_UNPAIRED_ELECTRONS {
        return None;
    }
    let lone_pairs = atom.lone_pairs().unwrap_or(lone_pairs_unassigned);
    let spin = match atom.multiplicity() {
        Some(m) => SpinState::try_new(unpaired, m)?,
        None => SpinState::max_multiplicity(unpaired)?,
    };
    AtomTypeSpec::new(
        element,
        charge,
        implicit_hydrogens,
        lone_pairs,
        unpaired,
        spin.multiplicity(),
        explicit_valence,
        donated_pairs,
        accepted_pairs,
        aromatic_valence,
        0,
    )
    .ok()
}
