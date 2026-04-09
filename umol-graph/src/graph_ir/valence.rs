//! Valence resolution strategies and atom typing support for GraphIR.

use smallvec::SmallVec;
use umol_data::{Element, SpinState, MAX_UNPAIRED_ELECTRONS};

use super::atom::Atom;
use super::atom_pattern::{AtomPattern, HydrogenPattern, Pattern};
use super::config::ValenceStrategy;
use super::config_data::{AtomTypeRegistry, NormalValenceTable, ValenceTable};
use super::molecule::AtomIndex;
use super::molecule_builder::MoleculeBuilder;
use crate::atom::AromaticValence;

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
    ) -> SmallVec<[Atom; 4]> {
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
) -> SmallVec<[Atom; 4]> {
    let mut pattern = AtomPattern::from_builder_atom(builder, atom_index);
    if pattern.implicit_hydrogens == HydrogenPattern::Normal {
        let inferred = infer_normal_implicit_hydrogens(builder, atom_index);
        let Some(hydrogens) = inferred else {
            return SmallVec::new();
        };
        pattern.implicit_hydrogens = HydrogenPattern::Is(hydrogens);
    }
    registry.candidates_for(&pattern)
}

fn infer_normal_implicit_hydrogens(builder: &MoleculeBuilder, atom_index: AtomIndex) -> Option<u8> {
    let atom = builder.atom(atom_index).expect("atom_index must be valid");
    let element = atom.element();
    let charge = match atom.charge {
        Pattern::Is(c) => c,
        Pattern::Any => 0,
    };
    let explicit_valence = builder.atom_bond_order_sum(atom_index);

    if builder.atom_aromatic_hint(atom_index) {
        if charge != 0 {
            return None;
        }
        return if element == Element::C {
            Some(3_u8.saturating_sub(explicit_valence))
        } else if matches!(
            element,
            Element::B
                | Element::N
                | Element::O
                | Element::P
                | Element::S
                | Element::Se
                | Element::As
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
) -> SmallVec<[Atom; 4]> {
    let atom = builder.atom(atom_index).expect("atom_index must be valid");
    let element = atom.element();
    let charge = match atom.charge {
        Pattern::Is(c) => c,
        Pattern::Any => 0,
    };
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
            element,
            charge,
            valence,
            donated_pairs,
            accepted_pairs,
            allow_implicit_hydrogens,
            builder.atom_has_normal_implicit_hydrogens(atom_index),
            atom,
        );
    }

    let implicit_hydrogens = if let HydrogenPattern::Is(h) = &atom.implicit_hydrogens {
        *h
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
        AromaticValence::NotAromatic,
        atom,
    )
}

#[allow(clippy::too_many_arguments)]
fn build_aromatic_spec(
    allowed_aromatic_valences: &[u8],
    element: Element,
    charge: i8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    allow_implicit_hydrogens: bool,
    normal_implicit_hydrogens: bool,
    atom: &AtomPattern,
) -> SmallVec<[Atom; 4]> {
    if allowed_aromatic_valences.is_empty() {
        return SmallVec::new();
    }

    // Element metadata is the canonical source of valence electron counts.
    let effective_electrons = (element.valence_electrons() as i16) - (charge as i16);
    let mut candidates = SmallVec::new();

    for &a in allowed_aromatic_valences {
        let sigma_budget = effective_electrons - (a as i16);
        if sigma_budget < valence as i16 {
            continue;
        }
        let implicit_hydrogens = match &atom.implicit_hydrogens {
            HydrogenPattern::Is(h) => *h,
            HydrogenPattern::Normal if normal_implicit_hydrogens => {
                let Some(h) = infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
                else {
                    continue;
                };
                h
            }
            HydrogenPattern::Normal => continue,
            HydrogenPattern::Any if normal_implicit_hydrogens => {
                let Some(h) = infer_normal_aromatic_implicit_hydrogens(element, charge, valence)
                else {
                    continue;
                };
                h
            }
            HydrogenPattern::Any => {
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
        if let Some(atom_out) = try_build_atom(
            element,
            charge,
            implicit_hydrogens,
            valence,
            donated_pairs,
            accepted_pairs,
            AromaticValence::Valence(a),
            atom,
        ) {
            candidates.push(atom_out);
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

#[allow(clippy::too_many_arguments)]
fn build_spec(
    element: Element,
    charge: i8,
    implicit_hydrogens: u8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    aromatic_valence: AromaticValence,
    atom: &AtomPattern,
) -> SmallVec<[Atom; 4]> {
    match try_build_atom(
        element,
        charge,
        implicit_hydrogens,
        valence,
        donated_pairs,
        accepted_pairs,
        aromatic_valence,
        atom,
    ) {
        Some(atom) => SmallVec::from_elem(atom, 1),
        None => SmallVec::new(),
    }
}

#[allow(clippy::too_many_arguments)]
fn try_build_atom(
    element: Element,
    charge: i8,
    implicit_hydrogens: u8,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    aromatic_valence: AromaticValence,
    atom: &AtomPattern,
) -> Option<Atom> {
    let total_valence = valence + implicit_hydrogens;
    // Element metadata is the canonical source of valence electron counts.
    let num_electrons = (element.valence_electrons() as i16) - (charge as i16);
    let unassigned = num_electrons - (total_valence as i16) - (aromatic_valence.valence() as i16);
    if unassigned < 0 {
        return None;
    }

    // Resolve (unpaired, lone_pairs) from one shared electron budget.
    // If the input fixes either value, infer the other consistently.
    let (unpaired, lone_pairs) = match (atom.unpaired_electrons, atom.lone_pairs) {
        (Pattern::Any, Pattern::Any) => ((unassigned % 2) as u8, (unassigned / 2) as u8),
        (Pattern::Is(unpaired), Pattern::Any) => {
            let remaining = unassigned - (unpaired as i16);
            if remaining < 0 || remaining % 2 != 0 {
                return None;
            }
            (unpaired, (remaining / 2) as u8)
        }
        (Pattern::Any, Pattern::Is(lone_pairs)) => {
            let remaining = unassigned - (2 * lone_pairs as i16);
            if remaining < 0 {
                return None;
            }
            (remaining as u8, lone_pairs)
        }
        (Pattern::Is(unpaired), Pattern::Is(lone_pairs)) => {
            if (unpaired as i16) + (2 * lone_pairs as i16) != unassigned {
                return None;
            }
            (unpaired, lone_pairs)
        }
    };
    if unpaired > MAX_UNPAIRED_ELECTRONS {
        return None;
    }

    let spin = match atom.multiplicity {
        Pattern::Is(m) => SpinState::try_new(unpaired, m).ok()?,
        Pattern::Any => SpinState::max_multiplicity(unpaired)?,
    };
    Atom::try_new(
        element,
        None,
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
