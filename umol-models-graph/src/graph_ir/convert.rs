//! Conversions between SimpleIR and GraphIR.

use umol::error::DataError;
use umol::Result;
use umol_data::Element;

use super::atom_matcher::STRICT_ATOM_MATCHER;
use super::atom_validator::STRICT_ATOM_VALIDATOR;
use super::bond_matcher::STRICT_BOND_MATCHER;
use super::{BondBuilder, BondDonation, BondOrder, Molecule, MoleculeBuilder};
use crate::simple_ir::{self as sir, AtomSymbol, BondSymbol};

/// Convert a SimpleIR molecule into a strictly validated GraphIR molecule.
pub fn sir_to_gir(src: &sir::Molecule) -> Result<Molecule> {
    let mut builder = MoleculeBuilder::with_capacity(src.atoms.len(), src.bonds.len());
    let mut atom_indices = Vec::with_capacity(src.atoms.len());

    for (idx, atom) in src.atoms.iter().enumerate() {
        let (element, isotope_hint) = match atom.symbol {
            AtomSymbol::Element(el) => (el, None),
            AtomSymbol::NamedIsotope(named) => (Element::from(named), Some(named.mass_number())),
            _ => {
                return Err(DataError::InvalidAtom(format!(
                    "SimpleIR atom {} is not supported in GraphIR",
                    idx
                ))
                .into())
            }
        };

        let (builder_idx, atom_builder) = builder.create_atom(element);

        if let Some(charge) = atom.charge {
            atom_builder.set_charge(charge);
        }

        atom_builder.set_span_opt(atom.span);

        if let Some(isotope) = atom.isotope.or(isotope_hint) {
            atom_builder.set_isotope(isotope);
        }

        if let Some(h) = atom.hydrogens {
            atom_builder.set_implicit_h(h);
        }

        if let Some(unpaired) = atom.unpaired_e {
            atom_builder.set_unpaired_e(unpaired);
        }

        if let Some(aromatic) = atom.aromatic {
            atom_builder.set_aromatic(aromatic);
        }

        if let Some(chirality) = atom.chirality {
            atom_builder.set_chirality(chirality);
        }

        if let Some(class_num) = atom.class {
            atom_builder.set_class_num(class_num);
        }

        if let Some(position) = atom.position.clone() {
            atom_builder.set_position(position);
        }

        atom_indices.push(builder_idx);
    }

    for (idx, bond) in src.bonds.iter().enumerate() {
        let a = usize::try_from(bond.start_atom).map_err(|_| {
            DataError::InvalidBond(format!(
                "bond {} references invalid start atom index {}",
                idx, bond.start_atom
            ))
        })?;
        let b = usize::try_from(bond.end_atom).map_err(|_| {
            DataError::InvalidBond(format!(
                "bond {} references invalid end atom index {}",
                idx, bond.end_atom
            ))
        })?;

        let &a_idx = atom_indices
            .get(a)
            .ok_or_else(|| DataError::MissingAtomIndex(a))?;
        let &b_idx = atom_indices
            .get(b)
            .ok_or_else(|| DataError::MissingAtomIndex(b))?;

        let order = convert_bond_order(bond.symbol, idx)?;

        let mut bond_builder = BondBuilder::new(order);
        bond_builder.set_donation(BondDonation::Shared);
        bond_builder.set_symbol(bond.symbol);
        bond_builder.set_direction(bond.direction);
        bond_builder.set_stereo(bond.stereo);
        bond_builder.set_ring(bond.ring);
        bond_builder.set_span_opt(bond.span);
        builder.add_bond(a_idx, b_idx, bond_builder)?;
    }

    builder.build_with(
        &STRICT_ATOM_VALIDATOR,
        &STRICT_ATOM_MATCHER,
        &STRICT_BOND_MATCHER,
    )
}

fn convert_bond_order(symbol: BondSymbol, bond_idx: usize) -> Result<BondOrder> {
    let order = match symbol {
        BondSymbol::Bond(order) => order,
        _ => {
            return Err(DataError::InvalidBond(format!(
                "bond {} uses unsupported symbol {:?}",
                bond_idx, symbol
            ))
            .into())
        }
    };

    let bond_order = match order {
        sir::BondOrder::Zero => BondOrder::Zero,
        sir::BondOrder::Single => BondOrder::Single,
        sir::BondOrder::Double => BondOrder::Double,
        sir::BondOrder::Triple => BondOrder::Triple,
        sir::BondOrder::Quadruple => BondOrder::Quadruple,
        sir::BondOrder::Aromatic | sir::BondOrder::Unknown => {
            return Err(DataError::InvalidBondOrder(format!(
                "unsupported bond order {:?} at bond {}",
                order, bond_idx
            ))
            .into())
        }
    };

    Ok(bond_order)
}

#[cfg(test)]
mod tests {
    use umol_data::Element;

    use super::*;

    #[test]
    fn convert_water() {
        let mut sir = sir::Molecule::default();

        sir.atoms.push(sir::Atom {
            symbol: AtomSymbol::Element(Element::O),
            ..Default::default()
        });
        sir.atoms.push(sir::Atom {
            symbol: AtomSymbol::Element(Element::H),
            ..Default::default()
        });
        sir.atoms.push(sir::Atom {
            symbol: AtomSymbol::Element(Element::H),
            ..Default::default()
        });

        sir.bonds.push(sir::Bond {
            start_atom: 0,
            end_atom: 1,
            symbol: BondSymbol::Bond(sir::BondOrder::Single),
            ..Default::default()
        });
        sir.bonds.push(sir::Bond {
            start_atom: 0,
            end_atom: 2,
            symbol: BondSymbol::Bond(sir::BondOrder::Single),
            ..Default::default()
        });

        let gir = sir_to_gir(&sir).expect("conversion must succeed");
        assert_eq!(gir.atom_count(), 3);
        assert_eq!(gir.bond_count(), 2);

        let oxygen = gir
            .atoms()
            .find(|a| a.element() == Element::O)
            .expect("oxygen must exist");
        assert_eq!(oxygen.valence(), 2);
        assert_eq!(oxygen.charge(), 0);
    }

    #[test]
    fn convert_rejects_non_element_atoms() {
        let mut sir = sir::Molecule::default();
        sir.atoms.push(sir::Atom::default());

        let result = sir_to_gir(&sir);
        assert!(result.is_err());
    }

    #[test]
    fn convert_rejects_aromatic_bonds() {
        let mut sir = sir::Molecule::default();
        sir.atoms.push(sir::Atom {
            symbol: AtomSymbol::Element(Element::C),
            ..Default::default()
        });
        sir.atoms.push(sir::Atom {
            symbol: AtomSymbol::Element(Element::C),
            ..Default::default()
        });
        sir.bonds.push(sir::Bond {
            start_atom: 0,
            end_atom: 1,
            symbol: BondSymbol::Bond(sir::BondOrder::Aromatic),
            ..Default::default()
        });

        let result = sir_to_gir(&sir);
        assert!(result.is_err());
    }
}
