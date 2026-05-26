//! Atom-typing valence resolver: narrows each atom against a registry of
//! `AtomAst` patterns keyed by element and (optionally) charge.

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AsLit, AtomAst, AtomConstraint, AtomId, AtomView, Lattice,
    MoleculeAst, MulticenterValenceAst, ValueAst,
};
use umol_shared::element::Element;

use crate::ops::valence::registry::AtomTypeRegistry;

#[derive(Clone, Debug)]
pub struct AtomTypingValence {
    pub registry: AtomTypeRegistry,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AtomTypingError {
    #[error("no atom-typing match for {atom:?} (element {element}, charge {charge:?})")]
    NoMatch {
        atom: AtomId,
        element: Element,
        charge: Option<i8>,
    },
}

impl AtomTypingValence {
    pub fn new(registry: AtomTypeRegistry) -> Self {
        Self { registry }
    }

    /// Iterates atoms, narrowing each non-ground atom against the registry.
    /// Returns `Err` on the first atom that has zero matching patterns;
    /// returns `Ok` if every atom either narrowed or stayed underdetermined
    /// without contradiction.
    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), AtomTypingError> {
        for i in 0..ast.atoms().count() as u32 {
            let idx = AtomId(i);
            let atom = ast.atom(idx);
            if atom.ast.is_ground() {
                continue;
            }

            let Some(element) = atom.element().as_lit() else {
                continue;
            };

            // Topological derived predicates should return Lit.
            let valence = atom.valence().as_lit_expect("valence should be Lit");
            let donated = atom.donated_pairs().as_lit_expect("donated pairs should be Lit");
            let accepted = atom.accepted_pairs().as_lit_expect("accepted pairs should be Lit");

            let match_input = self.build_match_input(&atom, valence, donated, accepted);
            let charge_key = match_input.charge.as_lit().map(|n| n as i8);
            let candidates: Vec<&AtomAst> = self
                .registry
                .lookup(element, charge_key)
                .iter()
                .filter(|pat| match_input.matches(pat))
                .collect();

            match candidates.len() {
                0 => {
                    return Err(AtomTypingError::NoMatch {
                        atom: idx,
                        element,
                        charge: charge_key,
                    });
                }
                1 => {
                    let cand = candidates[0];
                    let atom_mut = ast.atom_mut(idx).ast;
                    atom_mut.narrow_from(cand);
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// Project an atom into the lattice shape registry patterns are written
    /// against: pre-narrow `implicit_hydrogens` via the normal-valence table,
    /// and synthesize ground constraints for the topology-derived counts
    /// (localized valence, dative donated/accepted pairs) plus the
    /// membership-derived aromatic valence. `pattern.meet(&match_input)` then filters by those
    /// constraints directly via `AtomConstraints::meet`.
    fn build_match_input(
        &self,
        atom: &AtomView<'_>,
        valence: i64,
        donated: i64,
        accepted: i64,
    ) -> AtomAst {
        let mut match_input = atom.ast.clone();
        // First arm: idempotency — aromatic-system membership from a prior sweep.
        // Second arm: declared Aromatic(_) from the parser, before aromaticity perception runs.
        let is_aromatic = atom.is_in_aromatic_system()
            || AromaticValenceAst::aromatic(ValueAst::Undetermined)
                .matches(&match_input.constraints.aromatic_valence());

        let aromatic_constraint = if atom.is_in_aromatic_system() {
            atom.aromatic_valence().as_lit().map(|pi| {
                AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(pi)))
            })
        } else if !is_aromatic {
            Some(AtomConstraint::AromaticValence(
                AromaticValenceAst::NotAromatic,
            ))
        } else {
            None
        };

        let multicenter_constraint = atom.multicenter_valence().as_lit().map(|mc| {
            let mc_constraint = if mc == 0 {
                MulticenterValenceAst::NotMulticenter
            } else {
                MulticenterValenceAst::Multicenter(ValueAst::Lit(mc))
            };
            AtomConstraint::MulticenterValence(mc_constraint)
        });

        for constraint in [
            Some(AtomConstraint::Valence(ValueAst::Lit(valence))),
            Some(AtomConstraint::DonatedPairs(ValueAst::Lit(donated))),
            Some(AtomConstraint::AcceptedPairs(ValueAst::Lit(accepted))),
            aromatic_constraint,
            multicenter_constraint,
        ]
        .into_iter()
        .flatten()
        {
            match_input.constraints.add(constraint);
        }
        match_input
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::MoleculeAst;
    use umol_ast::{mol, mol_ground};

    use super::*;
    use crate::registry;

    fn methane() -> MoleculeAst {
        mol_ground!(r#"{:atoms ["C #h4"] :bonds []}"#)
    }

    fn methyl_chloride_partial() -> MoleculeAst {
        mol!(r#"{:atoms ["C #h3" "Cl"] :bonds [[0 1 "1"]]}"#)
    }

    #[rstest]
    fn test_atom_typing_valence_resolver_resolve_ground_methane_passthrough() {
        let reg = AtomTypeRegistry::default_registry().clone();
        let resolver =
            AtomTypingValence::new(reg);
        let mut ast = methane();
        resolver.resolve(&mut ast).unwrap();
    }

    #[rstest]
    fn test_atom_typing_valence_resolver_resolve_no_match() {
        let reg = registry!["C#c0#h4#n0#u0"];
        let resolver =
            AtomTypingValence::new(reg);
        let mut ast = methyl_chloride_partial();
        let err = resolver.resolve(&mut ast).unwrap_err();
        assert!(matches!(err, AtomTypingError::NoMatch { .. }));
    }

    #[rstest]
    fn test_atom_typing_valence_resolver_empty_registry_passes_through_ground_atoms() {
        let reg = AtomTypeRegistry::new();
        let resolver =
            AtomTypingValence::new(reg);
        let mut ast = methane();
        resolver.resolve(&mut ast).unwrap();
    }
}
