//! Resolution engine: refines a non-ground `MoleculeAst` toward a ground term.

use umol_ast::ast::atom::IsotopeAst;
use umol_ast::ast::spin::SpinStateAst;
use umol_ast::ast::value::ValueAst;

use crate::ast::molecule::MoleculeAst;
use crate::ast::rings::{RingCache, RingFamily};
use crate::ops::aromaticity::{AromaticityStrategy, AromaticityTheory};
use crate::ops::chemistry::Chemistry;
use crate::ops::error::ResolutionError;
use crate::ops::propagate::ValenceTheory;
use crate::ops::solution::Solution;

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

    pub fn resolve(
        self,
        ast: MoleculeAst,
    ) -> Result<Solution<MoleculeAst>, ResolutionError> {
        let mut cell = ResolverCell::new(ast);
        match cell.refine_valence(&self.chemistry.valence) {
            Progress::Contradictory => return Ok(Solution::Contradictory),
            Progress::Advanced | Progress::Fixpoint => {}
        }
        match cell.refine_aromaticity(&self.chemistry.aromaticity)? {
            Progress::Contradictory => return Ok(Solution::Contradictory),
            Progress::Advanced | Progress::Fixpoint => {}
        }
        Ok(cell.finalize())
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

    fn refine_valence(&mut self, valence: &ValenceTheory) -> Progress {
        valence.refine(&mut self.ast)
    }

    fn refine_aromaticity(
        &mut self,
        aromaticity: &AromaticityTheory,
    ) -> Result<Progress, ResolutionError> {
        let family = match &aromaticity.strategy {
            AromaticityStrategy::Clar(_) => RingFamily::InducedBenzenoid,
            AromaticityStrategy::HueckelRule(_) | AromaticityStrategy::Hmo(_) => {
                RingFamily::Simple
            }
        };
        let rings = self.rings.get(&self.ast, family, &aromaticity.ring_enumeration);
        Ok(aromaticity.refine(&mut self.ast, &rings)?)
    }

    fn finalize(mut self) -> Solution<MoleculeAst> {
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
            if let SpinStateAst { unpaired, multiplicity } = &mut bond.spin {
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
                bond.spin = SpinStateAst::from_state(state);
            }
        }
        for sys in self.ast.aromatic_systems_mut() {
            if matches!(sys.charge, ValueAst::Undetermined) {
                sys.charge = ValueAst::Lit(0);
            }
            if let SpinStateAst { unpaired, multiplicity } = &mut sys.spin {
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
                sys.spin = SpinStateAst::from_state(state);
            }
        }
        if self.ast.atoms().iter().all(|v| v.data.is_ground()) {
            Solution::Determined(self.ast)
        } else {
            Solution::Underdetermined(self.ast)
        }
    }
}
