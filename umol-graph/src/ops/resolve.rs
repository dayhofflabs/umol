//! Resolution engine: refines a non-ground `MoleculeAst` to a ground `Molecule`.

use umol_shared::atom_ast::IsotopeAst;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::molecule::MoleculeAst;
use crate::ast::rings::{RingCache, RingFamily};
use crate::api::Molecule;
use crate::ops::aromaticity::{AromaticityStrategy, AromaticityTheory};
use crate::ops::chemistry::Chemistry;
use crate::ops::error::ResolutionError;
use crate::ops::propagate::ValenceTheory;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Progress {
    Advanced,
    Fixpoint,
    Contradictory,
}

pub struct Resolver<'s> {
    chemistry: &'s Chemistry,
}

impl<'s> Resolver<'s> {
    pub fn new(chemistry: &'s Chemistry) -> Self {
        Self { chemistry }
    }

    pub fn resolve(self, ast: MoleculeAst) -> Result<Molecule, ResolutionError> {
        let mut cell = ResolverCell::new(ast);
        cell.refine_valence(&self.chemistry.valence)?;
        cell.refine_aromaticity(&self.chemistry.aromaticity)?;
        cell.finalize()
    }
}

/// Transient container during resolution.
///
/// Owns the `MoleculeAst` while the resolver mutates it and stores topology
/// caches populated as a side effect of perception. On finalize, transfers
/// the populated cache slots into the resulting `Molecule`.
struct ResolverCell {
    ast: MoleculeAst,
    rings: RingCache,
}

impl ResolverCell {
    fn new(ast: MoleculeAst) -> Self {
        Self {
            ast,
            rings: RingCache::new(),
        }
    }

    fn refine_valence(&mut self, valence: &ValenceTheory) -> Result<(), ResolutionError> {
        match valence.refine(&mut self.ast) {
            Progress::Contradictory => Err(ResolutionError::Contradictory),
            Progress::Advanced | Progress::Fixpoint => Ok(()),
        }
    }

    fn refine_aromaticity(
        &mut self,
        aromaticity: &AromaticityTheory,
    ) -> Result<(), ResolutionError> {
        let family = match &aromaticity.strategy {
            AromaticityStrategy::Clar(_) => RingFamily::InducedBenzenoid,
            AromaticityStrategy::HueckelRule(_) | AromaticityStrategy::Hmo(_) => {
                RingFamily::Simple
            }
        };
        let rings = self.rings.get(&self.ast, family, &aromaticity.ring_enumeration);
        match aromaticity.refine(&mut self.ast, &rings)? {
            Progress::Contradictory => Err(ResolutionError::Contradictory),
            Progress::Advanced | Progress::Fixpoint => Ok(()),
        }
    }

    fn finalize(mut self) -> Result<Molecule, ResolutionError> {
        // TODO: resolve correctly. Stopgap: ground the fields the resolver
        // currently leaves wildcard. Atoms get isotope=Natural, bonds get
        // charge=0 and spin=closed-shell.
        for atom in self.ast.atoms_mut() {
            if matches!(atom.isotope_mass, IsotopeAst::Undetermined) {
                atom.isotope_mass = IsotopeAst::Natural;
            }
        }
        for bond in self.ast.bonds_mut() {
            if matches!(bond.charge, ValueAst::Undetermined) {
                bond.charge = ValueAst::Lit(0);
            }
            if let SpinStateAst::Pair { unpaired, multiplicity } = &mut bond.spin {
                match (&*unpaired, &*multiplicity) {
                    (ValueAst::Undetermined, ValueAst::Undetermined) => {
                        *unpaired = ValueAst::Lit(0);
                        *multiplicity = ValueAst::Lit(1);
                    }
                    (ValueAst::Undetermined, ValueAst::Lit(m)) => {
                        *unpaired = ValueAst::Lit((m - 1).max(0));
                    }
                    (ValueAst::Lit(u), ValueAst::Undetermined) => {
                        *multiplicity = ValueAst::Lit(u + 1);
                    }
                    _ => {}
                }
            }
            if let Ok(Some(state)) = bond.spin.try_into_ground() {
                bond.spin = SpinStateAst::Lit(state);
            }
        }
        for sys in self.ast.aromatic_systems_mut() {
            if matches!(sys.charge, ValueAst::Undetermined) {
                sys.charge = ValueAst::Lit(0);
            }
            if let SpinStateAst::Pair { unpaired, multiplicity } = &mut sys.spin {
                match (&*unpaired, &*multiplicity) {
                    (ValueAst::Undetermined, ValueAst::Undetermined) => {
                        *unpaired = ValueAst::Lit(0);
                        *multiplicity = ValueAst::Lit(1);
                    }
                    (ValueAst::Undetermined, ValueAst::Lit(m)) => {
                        *unpaired = ValueAst::Lit((m - 1).max(0));
                    }
                    (ValueAst::Lit(u), ValueAst::Undetermined) => {
                        *multiplicity = ValueAst::Lit(u + 1);
                    }
                    _ => {}
                }
            }
            if let Ok(Some(state)) = sys.spin.try_into_ground() {
                sys.spin = SpinStateAst::Lit(state);
            }
        }
        Molecule::from_parts(self.ast, self.rings).map_err(|_| ResolutionError::Underdetermined)
    }
}
