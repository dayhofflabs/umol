//! Atom-typing valence resolver: narrows each atom against a registry of
//! `AtomAst` patterns keyed by element and (optionally) charge.

use thiserror::Error;
use umol_ast::ast::{
    AtomAst, AtomConstraint, AtomId, AtomView, ElementAst, ImplicitHydrogensAst, MoleculeAst,
    ValueAst,
};
use umol_shared::element::Element;

use crate::ops::valence::registry::AtomTypeRegistry;
use crate::ops::valence::shared::{
    atom_dative_counts, atom_is_aromatic, base_atom_compatible, charge_or_zero,
    infer_normal_implicit_hydrogens, lift_constraints, narrow_atom, pattern_constraints_compatible,
    AtomCandidate,
};

#[derive(Clone, Debug)]
pub struct AtomTypingValenceResolver {
    pub registry: AtomTypeRegistry,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AtomTypingError {
    #[error("no atom-typing match for {atom:?} (element {element}, charge {charge})")]
    NoMatchingPattern {
        atom: AtomId,
        element: Element,
        charge: i8,
    },
}

impl AtomTypingValenceResolver {
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
            if view
                .valence()
                .literal()
                .and_then(|n| u8::try_from(n).ok())
                .is_none()
            {
                continue;
            }

            let candidates = self.candidates_for(&view, element);
            match candidates.len() {
                0 => {
                    return Err(AtomTypingError::NoMatchingPattern {
                        atom: idx,
                        element,
                        charge: charge_or_zero(view.ast),
                    });
                }
                1 => {
                    let cand = candidates.into_iter().next().unwrap();
                    let atom_mut = ast.atom_mut(idx).ast;
                    narrow_atom(atom_mut, &cand.ast);
                    lift_constraints(atom_mut, &cand.lifted);
                }
                _ => {}
            }
        }
        Ok(())
    }

    fn candidates_for(&self, view: &AtomView<'_>, element: Element) -> Vec<AtomCandidate> {
        let atom = view.ast;
        let valence = match view.valence().literal().and_then(|n| u8::try_from(n).ok()) {
            Some(v) => v,
            None => return Vec::new(),
        };
        let charge = charge_or_zero(atom);
        let (donated_pairs, accepted_pairs) = match atom_dative_counts(view) {
            (Some(d), Some(a)) => (d, a),
            _ => return Vec::new(),
        };
        let is_aromatic = atom_is_aromatic(view);

        let implicit_h_constraint = match &atom.implicit_hydrogens {
            ImplicitHydrogensAst::Lit(n) => Some(*n as u8),
            ImplicitHydrogensAst::Normal => {
                let Some(h) =
                    infer_normal_implicit_hydrogens(element, charge, valence, is_aromatic)
                else {
                    return Vec::new();
                };
                Some(h)
            }
            ImplicitHydrogensAst::Undetermined => {
                infer_normal_implicit_hydrogens(element, charge, valence, is_aromatic)
            }
            _ => None,
        };

        let charge_key = match &atom.charge {
            ValueAst::Lit(n) => Some(*n as i8),
            _ => None,
        };

        self.registry
            .lookup(element, charge_key)
            .iter()
            .filter(|pattern| {
                pattern_implicit_h_compatible(pattern, implicit_h_constraint)
                    && pattern_constraints_compatible(
                        view,
                        &collect_pattern_constraints(pattern),
                        valence,
                        donated_pairs,
                        accepted_pairs,
                    )
                    && base_atom_compatible(atom, pattern)
            })
            .map(|pattern| AtomCandidate {
                ast: pattern.clone(),
                lifted: collect_pattern_constraints(pattern),
            })
            .collect()
    }
}

fn pattern_implicit_h_compatible(pattern: &AtomAst, implicit_h: Option<u8>) -> bool {
    match implicit_h {
        Some(h) => match &pattern.implicit_hydrogens {
            ImplicitHydrogensAst::Lit(n) => *n as u8 == h,
            _ => false,
        },
        None => true,
    }
}

fn collect_pattern_constraints(pattern: &AtomAst) -> Vec<AtomConstraint> {
    pattern.constraints.iter().cloned().collect()
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::MoleculeAst;
    use umol_ast::{mol, mol_zeroed};

    use super::*;
    use crate::registry;

    fn methane() -> MoleculeAst {
        mol_zeroed!(r#"{:atoms ["C #h4"] :bonds []}"#)
    }

    fn methyl_chloride_partial() -> MoleculeAst {
        mol!(r#"{:atoms ["C #h3" "Cl"] :bonds [[0 1 "1"]]}"#)
    }

    #[rstest]
    fn test_atom_typing_valence_resolver_resolve_ground_methane_passthrough() {
        let reg = AtomTypeRegistry::default_registry().clone();
        let resolver = AtomTypingValenceResolver::new(reg);
        let mut ast = methane();
        resolver.resolve(&mut ast).unwrap();
    }

    #[rstest]
    fn test_atom_typing_valence_resolver_resolve_no_match() {
        let reg = registry!["C#c0#h4#n0#u0"];
        let resolver = AtomTypingValenceResolver::new(reg);
        let mut ast = methyl_chloride_partial();
        let err = resolver.resolve(&mut ast).unwrap_err();
        assert!(matches!(err, AtomTypingError::NoMatchingPattern { .. }));
    }

    #[rstest]
    fn test_atom_typing_valence_resolver_empty_registry_passes_through_ground_atoms() {
        let reg = AtomTypeRegistry::new();
        let resolver = AtomTypingValenceResolver::new(reg);
        let mut ast = methane();
        resolver.resolve(&mut ast).unwrap();
    }
}
