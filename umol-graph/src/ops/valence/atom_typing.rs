//! Atom-typing valence resolver: narrows each atom against a registry of
//! `AtomAst` patterns keyed by element and (optionally) charge.

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AsLit, AtomAst, AtomConstraint, AtomId, AtomView, ElementAst, Lattice,
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
    NoMatchingPattern {
        atom: AtomId,
        element: Element,
        charge: ValueAst,
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
            let view = ast.atom(idx);
            if view.ast.is_ground() {
                continue;
            }
            let ElementAst::Lit(element) = view.ast.element else {
                continue;
            };
            let Some(valence) = view.valence().as_lit().and_then(|n| u8::try_from(n).ok()) else {
                continue;
            };
            let Some(donated) = view
                .donated_pairs()
                .as_lit()
                .and_then(|n| u8::try_from(n).ok())
            else {
                continue;
            };
            let Some(accepted) = view
                .accepted_pairs()
                .as_lit()
                .and_then(|n| u8::try_from(n).ok())
            else {
                continue;
            };

            let prepared = self.prepare_atom(&view, element, valence, donated, accepted);
            let charge_key = prepared.charge.as_lit().map(|n| n as i8);
            let compatibles: Vec<&AtomAst> = self
                .registry
                .lookup(element, charge_key)
                .iter()
                .filter(|pat| prepared.matches(pat))
                .collect();

            match compatibles.len() {
                0 => {
                    return Err(AtomTypingError::NoMatchingPattern {
                        atom: idx,
                        element,
                        charge: view.ast.charge.clone(),
                    });
                }
                1 => {
                    let cand = compatibles[0].clone();
                    let atom_mut = ast.atom_mut(idx).ast;
                    atom_mut.narrow_from(&cand);
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
    /// membership-derived aromatic valence. `pattern.meet(&prepared)` then filters by those
    /// constraints directly via `AtomConstraints::meet`.
    fn prepare_atom(
        &self,
        view: &AtomView<'_>,
        element: Element,
        valence: u8,
        donated: u8,
        accepted: u8,
    ) -> AtomAst {
        let mut prepared = view.ast.clone();
        // First arm: idempotency — aromatic-system membership from a prior sweep.
        // Second arm: declared Aromatic(_) from the parser, before aromaticity perception runs.
        // Neither arm requires aromaticity to run ahead of valence.
        let is_aromatic = view.is_in_aromatic_system()
            || AromaticValenceAst::aromatic(ValueAst::Undetermined)
                .matches(&prepared.constraints.aromatic_valence());
        prepared
            .constraints
            .add(AtomConstraint::Valence(ValueAst::Lit(valence as i64)));
        prepared
            .constraints
            .add(AtomConstraint::DonatedPairs(ValueAst::Lit(donated as i64)));
        prepared
            .constraints
            .add(AtomConstraint::AcceptedPairs(ValueAst::Lit(
                accepted as i64,
            )));
        if view.is_in_aromatic_system() {
            if let Some(pi) = view
                .aromatic_valence()
                .as_lit()
                .and_then(|n| u8::try_from(n).ok())
            {
                prepared.constraints.add(AtomConstraint::AromaticValence(
                    AromaticValenceAst::Aromatic(ValueAst::Lit(pi as i64)),
                ));
            }
        } else if !is_aromatic {
            prepared.constraints.add(AtomConstraint::AromaticValence(
                AromaticValenceAst::NotAromatic,
            ));
        }
        if let Some(mc) = view
            .multicenter_valence()
            .as_lit()
            .and_then(|n| u8::try_from(n).ok())
        {
            let mc_constraint = if mc == 0 {
                MulticenterValenceAst::NotMulticenter
            } else {
                MulticenterValenceAst::Multicenter(ValueAst::Lit(mc as i64))
            };
            prepared
                .constraints
                .add(AtomConstraint::MulticenterValence(mc_constraint));
        }
        prepared
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
        assert!(matches!(err, AtomTypingError::NoMatchingPattern { .. }));
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
