//! Valence resolution strategies and atom typing support for GraphIR.

use smallvec::SmallVec;
use umol_data::{Element, SpinState, MAX_UNPAIRED_ELECTRONS};

use crate::graph_ir::atom::AtomBuilder;
use crate::graph_ir::atom_type::{AtomTypeQuery, AtomTypeSpec, HydrogenConstraint};
use crate::atom::{AromaticValence, ImplicitHydrogens};
use crate::graph_ir::config::ValenceStrategy;
use crate::graph_ir::config_data::{AtomTypeRegistry, NormalValenceTable, ValenceEntry, ValenceTable};
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
            Self::AtomTyping { registry } => atom_typing_candidates(registry, builder, atom_index),
            Self::Counts {
                table,
                allow_implicit_hydrogens,
            } => counts_candidates(table, *allow_implicit_hydrogens, builder, atom_index),
        }
    }
}

fn atom_typing_candidates(
    registry: &AtomTypeRegistry,
    builder: &MoleculeBuilder,
    atom_index: AtomIndex,
) -> SmallVec<[AtomTypeSpec; 4]> {
    let mut query = AtomTypeQuery::from_builder_atom(builder, atom_index);
    if query.implicit_hydrogens == Some(HydrogenConstraint::Normal) {
        let inferred = infer_normal_implicit_hydrogens(builder, atom_index);
        let Some(hydrogens) = inferred else {
            return SmallVec::new();
        };
        query.implicit_hydrogens = Some(HydrogenConstraint::Hydrogens(hydrogens));
    }
    registry.candidates_for(&query)
}

fn infer_normal_implicit_hydrogens(builder: &MoleculeBuilder, atom_index: AtomIndex) -> Option<u8> {
    let atom = builder.atom(atom_index).expect("atom_index must be valid");
    let element = atom.element();
    let charge = atom.charge().unwrap_or(0);
    let explicit_valence = builder.atom_bond_order_sum(atom_index);

    if builder.atom_aromatic_hint(atom_index) {
        if charge != 0 {
            return None;
        }
        return if element == Element::C {
            Some(3_u8.saturating_sub(explicit_valence))
        } else if matches!(
            element,
            Element::B | Element::N | Element::O | Element::P | Element::S | Element::Se | Element::As
        ) {
            Some(0)
        } else {
            None
        };
    }

    let normal_valence = NormalValenceTable::default_table().normal_valence_for(element, charge)?;
    Some(normal_valence.saturating_sub(explicit_valence))
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
    let valence = builder.atom_bond_order_sum(atom_index);
    let (donated_pairs, accepted_pairs) = builder.atom_dative_bond_order_sums(atom_index);

    let entry = match table.entry(element) {
        Some(e) => e,
        None => return SmallVec::new(),
    };

    if builder.atom_aromatic_hint(atom_index) {
        let aromatic_valences = if charge != 0 {
            element
                .shift(-charge)
                .and_then(|e| table.entry(e))
                .map(|e| e.allowed_aromatic_valences.as_slice())
                .unwrap_or(entry.allowed_aromatic_valences.as_slice())
        } else {
            entry.allowed_aromatic_valences.as_slice()
        };
        return build_aromatic_spec(
            aromatic_valences,
            entry,
            element,
            charge,
            valence,
            donated_pairs,
            accepted_pairs,
            allow_implicit_hydrogens,
            &atom,
        );
    }

    let implicit_hydrogens = if let Some(h) = atom.hydrogen_count() {
        h
    } else if allow_implicit_hydrogens {
        match table.compute_implicit_hydrogens(element, charge, valence) {
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
        valence,
        donated_pairs,
        accepted_pairs,
        AromaticValence::None,
        &atom,
        entry,
    )
}

fn build_aromatic_spec(
    allowed_aromatic_valences: &[u8],
    entry: &ValenceEntry,
    element: Element,
    charge: i8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    allow_implicit_hydrogens: bool,
    atom: &AtomBuilder,
) -> SmallVec<[AtomTypeSpec; 4]> {
    if allowed_aromatic_valences.is_empty() {
        return SmallVec::new();
    }

    let effective_electrons = (entry.outer_electrons as i16) - (charge as i16);
    let mut candidates = SmallVec::new();

    for &a in allowed_aromatic_valences {
        let sigma_budget = effective_electrons - (a as i16);
        if sigma_budget < valence as i16 {
            continue;
        }
        let implicit_hydrogens = match atom.implicit_hydrogens() {
            Some(ImplicitHydrogens::Hydrogens(h)) => h,
            Some(ImplicitHydrogens::Normal) => {
                let Some(h) =
                    infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
                else {
                    continue;
                };
                h
            }
            None => {
                if allow_implicit_hydrogens {
                    (sigma_budget - valence as i16) as u8
                } else {
                    0
                }
            }
        };
        if implicit_hydrogens > 1 {
            continue;
        }
        let total_sigma = valence + implicit_hydrogens;
        let remaining = effective_electrons - total_sigma as i16 - a as i16;
        if remaining < 0 || remaining % 2 != 0 {
            continue;
        }
        if let Some(spec) = try_build_spec(
            element,
            charge,
            implicit_hydrogens,
            valence,
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

fn infer_normal_aromatic_implicit_hydrogens(
    element: Element,
    charge: i8,
    valence: u8,
) -> Option<u8> {
    if charge != 0 {
        return None;
    }

    if element == Element::C {
        Some(3 - valence)
    } else if matches!(
        element,
        Element::B | Element::N | Element::O | Element::P | Element::S | Element::Se | Element::As
    ) {
        Some(0)
    } else {
        None
    }
}

fn build_spec(
    element: Element,
    charge: i8,
    implicit_hydrogens: u8,
    valence: u8,
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
        valence,
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
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    aromatic_valence: AromaticValence,
    atom: &AtomBuilder,
    entry: &ValenceEntry,
) -> Option<AtomTypeSpec> {
    let total_valence = valence + implicit_hydrogens;
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
        valence,
        donated_pairs,
        accepted_pairs,
        aromatic_valence,
        0,
    )
    .ok()
}
