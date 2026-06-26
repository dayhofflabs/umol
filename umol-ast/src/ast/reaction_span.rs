//! Reaction span AST: the superimposed `L ∪_K R` graph encoding a reaction's DPO rule span.
//!
//! Distinct from `ReactionAst` (operational — `lhs` + `deltas`): this is a *materialized*
//! superimposed graph carrying, per atom/bond, both its before and after state plus a
//! membership tag. The DPO span `L ←K─ R` is read off the tags — `K = Unchanged ∪ Modified`,
//! `L = K ∪ Removed`, `R = K ∪ Added` — and `right()` / `left()` project the two sides back to
//! a `MoleculeAst`. `Modified` (a preserved element relabeled across the reaction) is the
//! relabeling-DPO reading: the element persists in `K`, its label resolved per side.
//!
//! Localized topology only. Molecule-level constraints and overlays stay on the operational
//! `ReactionAst`; they are not represented here yet.

use std::collections::{BTreeMap, HashMap};

use umol_graph_core::{EdgeId, Graph};

use super::atom::AtomAst;
use super::bond::BondAst;
use super::delta::{apply_atom_change, apply_bond_change, AtomDelta, BondDelta, Delta};
use super::error::Contradiction;
use super::id::{AtomId, BondId};
use super::molecule::MoleculeAst;
use super::reaction::ReactionAst;
use super::traits::Canonicalize;

/// One element's superimposed state: its DPO membership plus the value(s) it carries.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Change<T> {
    /// In the interface `K` — present and identical on both sides.
    Unchanged(T),
    /// In the interface `K` — present on both sides but relabeled (a dynamic element).
    Modified { left: T, right: T },
    /// In `R` only — created.
    Added(T),
    /// In `L` only — deleted.
    Removed(T),
}

impl<T> Change<T> {
    /// The left-side (`L`) value, or `None` if the element is created.
    pub fn left(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Removed(value) | Self::Modified { left: value, .. } => {
                Some(value)
            }
            Self::Added(_) => None,
        }
    }

    /// The right-side (`R`) value, or `None` if the element is deleted.
    pub fn right(&self) -> Option<&T> {
        match self {
            Self::Unchanged(value) | Self::Added(value) | Self::Modified { right: value, .. } => {
                Some(value)
            }
            Self::Removed(_) => None,
        }
    }
}

/// The superimposed reaction graph — the reaction's DPO rule span, materialized. The union
/// topology is the `lhs` frame (deleted elements kept as nodes/edges) with created elements
/// appended; `atoms` / `bonds` are indexed parallel to the graph's nodes / edges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionSpanAst {
    graph: Graph,
    atoms: Vec<Change<AtomAst>>,
    bonds: Vec<Change<BondAst>>,
}

impl ReactionSpanAst {
    pub fn graph(&self) -> &Graph {
        &self.graph
    }

    pub fn atoms(&self) -> &[Change<AtomAst>] {
        &self.atoms
    }

    pub fn bonds(&self) -> &[Change<BondAst>] {
        &self.bonds
    }

    /// The left-hand (reactant) molecule: every element present on the left, in a compacted
    /// frame (created elements dropped).
    pub fn left(&self) -> MoleculeAst {
        self.project(|atom| atom.left(), |bond| bond.left())
    }

    /// The right-hand (product) molecule: every element present on the right, in a compacted
    /// frame (deleted elements dropped).
    pub fn right(&self) -> MoleculeAst {
        self.project(|atom| atom.right(), |bond| bond.right())
    }

    /// Project one side to a `MoleculeAst`. `atom_side` / `bond_side` pick the left or right
    /// value of each element; absent elements are dropped and the survivors are renumbered.
    fn project(
        &self,
        atom_side: impl Fn(&Change<AtomAst>) -> Option<&AtomAst>,
        bond_side: impl Fn(&Change<BondAst>) -> Option<&BondAst>,
    ) -> MoleculeAst {
        let mut compacted: Vec<Option<AtomId>> = vec![None; self.atoms.len()];
        let mut atoms: Vec<AtomAst> = Vec::new();
        for (node, change) in self.atoms.iter().enumerate() {
            if let Some(ast) = atom_side(change) {
                compacted[node] = Some(AtomId(atoms.len() as u32));
                atoms.push(ast.clone());
            }
        }
        let mut bonds: Vec<(AtomId, AtomId, BondAst)> = Vec::new();
        for (edge, change) in self.bonds.iter().enumerate() {
            if let Some(ast) = bond_side(change) {
                let [a, b] = self.graph.edge_endpoints(EdgeId(edge as u32));
                if let (Some(a), Some(b)) = (compacted[a.index()], compacted[b.index()]) {
                    bonds.push((a, b, ast.clone()));
                }
            }
        }
        MoleculeAst::from_atoms_and_bonds(atoms, bonds)
    }
}

impl ReactionAst {
    /// Materialize the superimposed reaction span. Canonicalizes the deltas, then annotates
    /// each `lhs` element (in its own frame) with its before/after state — `Removed` /
    /// `Added` / `Modified` / `Unchanged` — appending created elements. A `Modified` element's
    /// right value is its left value with the entity's field/constraint changes applied.
    /// `Err(Contradiction)` if the deltas are inconsistent (or inconsistent with `lhs`).
    pub fn to_reaction_span(&self) -> Result<ReactionSpanAst, Contradiction> {
        let deltas = self.deltas.clone().canonicalize()?;
        let lhs = &self.lhs;
        let atom_count = lhs.atoms().count();
        let bond_count = lhs.bonds().count();

        let mut removed_atoms: HashMap<AtomId, AtomAst> = HashMap::new();
        let mut added_atoms: BTreeMap<AtomId, AtomAst> = BTreeMap::new();
        let mut atom_changes: HashMap<AtomId, Vec<AtomDelta>> = HashMap::new();
        let mut removed_bonds: HashMap<BondId, BondAst> = HashMap::new();
        let mut added_bonds: BTreeMap<BondId, ([AtomId; 2], BondAst)> = BTreeMap::new();
        let mut bond_changes: HashMap<BondId, Vec<BondDelta>> = HashMap::new();

        for delta in deltas.iter() {
            match delta {
                Delta::Atom(atom) => match atom {
                    AtomDelta::Remove { id, ast } => {
                        removed_atoms.insert(*id, ast.clone());
                    }
                    AtomDelta::Add { id, ast } => {
                        added_atoms.insert(*id, ast.clone());
                    }
                    AtomDelta::SetField { id, .. } | AtomDelta::SetConstraint { id, .. } => {
                        atom_changes.entry(*id).or_default().push(atom.clone());
                    }
                },
                Delta::Bond(bond) => match bond {
                    BondDelta::Remove { id, ast, .. } => {
                        removed_bonds.insert(*id, ast.clone());
                    }
                    BondDelta::Add { id, endpoints, ast } => {
                        added_bonds.insert(*id, (*endpoints, ast.clone()));
                    }
                    BondDelta::SetField { id, .. } | BondDelta::SetConstraint { id, .. } => {
                        bond_changes.entry(*id).or_default().push(bond.clone());
                    }
                },
                // Molecule-level constraints are carried by the operational form, not the span.
                Delta::Constraint(_) => {}
            }
        }

        // Union node index per frame atom id: lhs atoms keep their id, created atoms append.
        let mut atom_index: HashMap<AtomId, usize> =
            HashMap::with_capacity(atom_count + added_atoms.len());
        for node in 0..atom_count {
            atom_index.insert(AtomId(node as u32), node);
        }
        for (offset, &id) in added_atoms.keys().enumerate() {
            atom_index.insert(id, atom_count + offset);
        }

        let mut atoms: Vec<Change<AtomAst>> = Vec::with_capacity(atom_count + added_atoms.len());
        for node in 0..atom_count {
            let id = AtomId(node as u32);
            if let Some(ast) = removed_atoms.get(&id) {
                atoms.push(Change::Removed(ast.clone()));
            } else if let Some(changes) = atom_changes.get(&id) {
                let left = lhs.atom(id).ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_atom_change(&mut right, change)?;
                }
                atoms.push(Change::Modified { left, right });
            } else {
                atoms.push(Change::Unchanged(lhs.atom(id).ast.clone()));
            }
        }
        for ast in added_atoms.into_values() {
            atoms.push(Change::Added(ast));
        }

        let mut bonds: Vec<Change<BondAst>> = Vec::with_capacity(bond_count + added_bonds.len());
        let mut edges: Vec<[u32; 2]> = Vec::with_capacity(bond_count + added_bonds.len());
        for edge in 0..bond_count {
            let id = BondId(edge as u32);
            let [a, b] = lhs.raw_graph().edge_endpoints(EdgeId(edge as u32));
            edges.push([a.0, b.0]);
            if let Some(ast) = removed_bonds.get(&id) {
                bonds.push(Change::Removed(ast.clone()));
            } else if let Some(changes) = bond_changes.get(&id) {
                let left = lhs.bond(id).ast.clone();
                let mut right = left.clone();
                for change in changes {
                    apply_bond_change(&mut right, change)?;
                }
                bonds.push(Change::Modified { left, right });
            } else {
                bonds.push(Change::Unchanged(lhs.bond(id).ast.clone()));
            }
        }
        for (endpoints, ast) in added_bonds.into_values() {
            edges.push([
                atom_index[&endpoints[0]] as u32,
                atom_index[&endpoints[1]] as u32,
            ]);
            bonds.push(Change::Added(ast));
        }

        let graph = Graph::new(atoms.len(), &edges);
        Ok(ReactionSpanAst {
            graph,
            atoms,
            bonds,
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::super::delta::Deltas;
    use super::super::edit::BondFieldChange;
    use super::super::value::ValueAst;
    use super::*;

    #[rstest]
    #[case::unchanged(Change::Unchanged(5), Some(&5))]
    #[case::modified(Change::Modified { left: 1, right: 2 }, Some(&1))]
    #[case::removed(Change::Removed(7), Some(&7))]
    #[case::added(Change::Added(9), None)]
    fn test_change_left(#[case] change: Change<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(change.left(), expected);
    }

    #[rstest]
    #[case::unchanged(Change::Unchanged(5), Some(&5))]
    #[case::modified(Change::Modified { left: 1, right: 2 }, Some(&2))]
    #[case::added(Change::Added(9), Some(&9))]
    #[case::removed(Change::Removed(7), None)]
    fn test_change_right(#[case] change: Change<i32>, #[case] expected: Option<&i32>) {
        assert_eq!(change.right(), expected);
    }

    #[rstest]
    fn test_reaction_ast_to_reaction_span() {
        let reaction = ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::C),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([Delta::Bond(BondDelta::SetField {
                id: BondId(0),
                change: BondFieldChange::Order {
                    old: ValueAst::Lit(1),
                    new: ValueAst::Lit(2),
                },
            })]),
        );
        let span = reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.atoms(),
            [
                Change::Unchanged(AtomAst::from_element(Element::C)),
                Change::Unchanged(AtomAst::from_element(Element::C)),
            ],
        );
        assert_eq!(
            span.bonds(),
            [Change::Modified {
                left: BondAst::from_order(1),
                right: BondAst::from_order(2),
            }],
        );
    }

    #[rstest]
    fn test_reaction_span_ast_right(substitution_reaction: ReactionAst) {
        let span = substitution_reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.atoms(),
            [
                Change::Unchanged(AtomAst::from_element(Element::C)),
                Change::Removed(AtomAst::from_element(Element::O)),
                Change::Added(AtomAst::from_element(Element::N)),
            ],
        );
        assert_eq!(
            span.bonds(),
            [
                Change::Removed(BondAst::from_order(1)),
                Change::Added(BondAst::from_order(1)),
            ],
        );
        assert_eq!(
            span.right(),
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::N),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
        );
    }

    #[rstest]
    fn test_reaction_span_ast_left(substitution_reaction: ReactionAst) {
        let span = substitution_reaction.to_reaction_span().unwrap();
        assert_eq!(
            span.left(),
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
        );
    }

    // C-O with atom 1 (O) and its bond removed, replaced by a new N (atom 2) bonded to C.
    #[fixture]
    fn substitution_reaction() -> ReactionAst {
        ReactionAst::new(
            MoleculeAst::from_atoms_and_bonds(
                vec![
                    AtomAst::from_element(Element::C),
                    AtomAst::from_element(Element::O),
                ],
                vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            ),
            Deltas::from_iter([
                Delta::Atom(AtomDelta::Remove {
                    id: AtomId(1),
                    ast: AtomAst::from_element(Element::O),
                }),
                Delta::Bond(BondDelta::Remove {
                    id: BondId(0),
                    endpoints: [AtomId(0), AtomId(1)],
                    ast: BondAst::from_order(1),
                }),
                Delta::Atom(AtomDelta::Add {
                    id: AtomId(2),
                    ast: AtomAst::from_element(Element::N),
                }),
                Delta::Bond(BondDelta::Add {
                    id: BondId(1),
                    endpoints: [AtomId(0), AtomId(2)],
                    ast: BondAst::from_order(1),
                }),
            ]),
        )
    }
}
