//! Atom-typing valence resolver: narrows each atom against a registry of
//! `AtomAst` patterns keyed by element and (optionally) charge.

use thiserror::Error;
use umol_ast::ast::{
    AromaticValenceAst, AsLit, AtomAst, AtomConstraint, AtomId, AtomView, Lattice, MoleculeAst,
    MulticenterValenceAst, ValueAst,
};
use umol_shared::element::Element;

use super::compare::compare_valence_preference;
use super::registry::AtomTypeRegistry;

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
    /// Returns `Err` on the first atom that has zero matching patterns.
    /// Multiple matches are resolved via [`compare_valence_preference`].
    pub fn resolve(&self, ast: &mut MoleculeAst) -> Result<(), AtomTypingError> {
        for i in 0..ast.atoms().count() as u32 {
            let idx = AtomId(i);
            let atom = ast.atom(idx);
            if atom.is_ground() {
                continue;
            }

            let Some(element) = atom.element().as_lit() else {
                continue;
            };

            // Topological derived predicates should return Lit.
            let pattern = add_constraints(&atom);
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
                        atom: idx,
                        element,
                        charge,
                    });
                }
                1 => {
                    let cand = candidates[0];
                    let atom_mut = ast.atom_mut(idx).ast;
                    atom_mut.narrow_from(cand);
                }
                _ => {
                    let best = candidates
                        .into_iter()
                        .max_by(|a, b| compare_valence_preference(a, b))
                        .unwrap();
                    let atom_mut = ast.atom_mut(idx).ast;
                    atom_mut.narrow_from(best);
                }
            }
        }
        Ok(())
    }
}

/// Add constraints to atom AST for registry pattern matching.
fn add_constraints(atom: &AtomView<'_>) -> AtomAst {
    let mut updated = atom.ast.clone();

    updated
        .constraints
        .add(AtomConstraint::valence(ValueAst::Lit(
            atom.valence().as_lit_expect("valence should be Lit"),
        )));
    updated
        .constraints
        .add(AtomConstraint::donated_pairs(ValueAst::Lit(
            atom.donated_pairs()
                .as_lit_expect("donated pairs should be Lit"),
        )));
    updated
        .constraints
        .add(AtomConstraint::accepted_pairs(ValueAst::Lit(
            atom.accepted_pairs()
                .as_lit_expect("accepted pairs should be Lit"),
        )));

    if atom.is_in_aromatic_system() {
        let aromatic = atom.aromatic_valence().as_lit_or(0);
        updated.constraints.add(AtomConstraint::aromatic_valence(
            AromaticValenceAst::aromatic(aromatic),
        ));
    } else if atom.neighbors().any(|n| n.bond().constraints().aromatic()) {
        updated.constraints.add(AtomConstraint::aromatic_valence(
            AromaticValenceAst::aromatic(ValueAst::Undetermined),
        ));
    };

    if atom.is_in_multicenter_bond() {
        let multicenter = atom.multicenter_valence().as_lit_or(0);
        updated.constraints.add(AtomConstraint::multicenter_valence(
            MulticenterValenceAst::multicenter(multicenter),
        ));
    };

    updated
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
