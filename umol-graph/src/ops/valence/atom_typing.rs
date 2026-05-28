//! Atom-typing valence resolver: narrows each atom against a registry of
//! `AtomAst` patterns keyed by element and (optionally) charge.

use thiserror::Error;
use umol_ast::ast::{AsLit, AtomAst, AtomId, Lattice, MoleculeAst};
use umol_shared::element::Element;

use super::compare::compare_valence_preference;
use super::registry::AtomTypeRegistry;

#[derive(Clone, Debug)]
pub struct AtomTypingValence {
    pub registry: AtomTypeRegistry,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AtomTypingError {
    #[error("no atom-typing match for {atom_id:?} (element {element}, charge {charge:?})")]
    NoMatch {
        atom_id: AtomId,
        element: Element,
        charge: Option<i8>,
    },
}

impl AtomTypingValence {
    pub fn new(registry: AtomTypeRegistry) -> Self {
        Self { registry }
    }

    /// Iterates atoms, narrowing each non-ground atom against the registry.
    /// Returns `Err` on the first atom that has zero matching patterns.
    /// Multiple matches are resolved via [`compare_valence_preference`].
    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), AtomTypingError> {
        for id in ast.atoms().ids() {
            let atom = ast.atom(id);
            if atom.is_ground() {
                continue;
            }

            let Some(element) = atom.element().as_lit() else {
                continue;
            };

            // Topological derived predicates should return Lit.
            let constraints = atom.derive_constraints();
            let pattern = atom.ast.clone().with_constraints(constraints);
            let charge = pattern.charge.as_lit().map(|n| n as i8);
            let candidates: Vec<&AtomAst> = self
                .registry
                .lookup(element, charge)
                .iter()
                .filter(|entry| pattern.matches(entry))
                .collect();

            match candidates.len() {
                0 => {
                    return Err(AtomTypingError::NoMatch {
                        atom_id: id,
                        element,
                        charge,
                    });
                }
                1 => {
                    let cand = candidates[0];
                    let atom_mut = ast.atom_mut(id).ast;
                    atom_mut.narrow_from(cand);
                }
                _ => {
                    let best = candidates
                        .into_iter()
                        .max_by(|a, b| compare_valence_preference(a, b))
                        .unwrap();
                    let atom_mut = ast.atom_mut(id).ast;
                    atom_mut.narrow_from(best);
                }
            }
        }
        Ok(())
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
        let resolver = AtomTypingValence::new(reg);
        let mut ast = methane();
        resolver.resolve(&mut ast).unwrap();
    }

    #[rstest]
    fn test_atom_typing_valence_resolver_resolve_no_match() {
        let reg = registry!["C#c0#h4#n0#u0"];
        let resolver = AtomTypingValence::new(reg);
        let mut ast = methyl_chloride_partial();
        let err = resolver.resolve(&mut ast).unwrap_err();
        assert!(matches!(err, AtomTypingError::NoMatch { .. }));
    }

    #[rstest]
    fn test_atom_typing_valence_resolver_empty_registry_passes_through_ground_atoms() {
        let reg = AtomTypeRegistry::new();
        let resolver = AtomTypingValence::new(reg);
        let mut ast = methane();
        resolver.resolve(&mut ast).unwrap();
    }
}
