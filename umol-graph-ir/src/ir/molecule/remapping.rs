use std::sync::Arc;

use umol_graph_core::{EdgeId, FixedRelationSet, FixedVarBirelationSet, Graph, NodeId, Remapping};

use super::super::aromatic::AromaticSystems;
use super::super::dative::DativeBonds;
use super::super::multicenter::MulticenterBonds;
use super::Molecule;
use crate::ir::{Constraints, MoleculeCorrespondence};

impl Molecule {
    /// Relabel this molecule into the dense target id spaces described by `correspondence`.
    ///
    /// The correspondence must describe every entity in this molecule and be total on both sides
    /// for all eight entity families. The operation transports topology, relation participants,
    /// position-sensitive relation data, stereo frames, entity forms, and constraint references.
    /// It does not validate chemistry, normalize attributes, repair references, compact tables, or
    /// remove entities.
    ///
    /// # Panics
    ///
    /// Panics when this molecule fails its representation-integrity contract or `correspondence`
    /// does not describe a complete dense renumbering of it. Use [`Self::try_remap`] for
    /// independently supplied values.
    ///
    /// # Semantic properties
    ///
    /// The result is equivalent to `self` under `correspondence`. Identity remapping is exact,
    /// inverse remapping recovers the original molecule, and sequential remapping agrees with
    /// correspondence composition.
    pub fn remap(&self, correspondence: &MoleculeCorrespondence) -> Self {
        self.try_remap(correspondence).expect(
            "molecule remapping requires an integrity-valid source and a complete dense correspondence",
        )
    }

    /// Checked form of [`Self::remap`].
    ///
    /// Returns `None` when this molecule fails its representation-integrity contract, when the
    /// correspondence's source counts differ from the molecule's entity counts, or when any entity
    /// family is not a bijection onto a dense target id space.
    pub fn try_remap(&self, correspondence: &MoleculeCorrespondence) -> Option<Self> {
        self.check_integrity().ok()?;

        let counts_match = [
            (correspondence.atoms().left_count(), self.atoms.len()),
            (correspondence.bonds().left_count(), self.bonds.len()),
            (
                correspondence.dative_bonds().left_count(),
                self.dative_bonds.count(),
            ),
            (
                correspondence.aromatic_systems().left_count(),
                self.aromatic_systems.count(),
            ),
            (
                correspondence.multicenter_bonds().left_count(),
                self.multicenter_bonds.count(),
            ),
            (
                correspondence.noncovalent_bonds().left_count(),
                self.noncovalent_bonds.count(),
            ),
            (
                correspondence.stereo_atoms().left_count(),
                self.stereo_atoms.count(),
            ),
            (
                correspondence.stereo_bonds().left_count(),
                self.stereo_bonds.count(),
            ),
        ]
        .into_iter()
        .all(|(mapped, actual)| mapped == actual);
        if !counts_match || !correspondence.is_total() {
            return None;
        }

        let graph_remapping = Remapping::new(
            correspondence
                .atoms()
                .matched_pairs()
                .iter()
                .map(|&(_, right)| NodeId::from(right))
                .collect(),
            correspondence
                .bonds()
                .matched_pairs()
                .iter()
                .map(|&(_, right)| EdgeId::from(right))
                .collect(),
        );
        let id_remapping = correspondence.to_remapping()?;

        let atoms = reorder(
            self.atoms.as_ref().clone(),
            correspondence.atoms().right_count(),
            correspondence
                .atoms()
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index())),
        )?;
        let bonds = reorder(
            self.bonds.as_ref().clone(),
            correspondence.bonds().right_count(),
            correspondence
                .bonds()
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index())),
        )?;
        let edges = reorder(
            self.graph
                .edge_ids()
                .map(|edge| {
                    let [first, second] = self.graph.edge_endpoints(edge);
                    [
                        graph_remapping.map_node(first).0,
                        graph_remapping.map_node(second).0,
                    ]
                })
                .collect(),
            correspondence.bonds().right_count(),
            correspondence
                .bonds()
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index())),
        )?;
        let graph = Graph::new(atoms.len(), &edges);

        let dative_bonds = DativeBonds::new(reorder(
            self.dative_bonds.remap(&graph_remapping).into_entries(),
            correspondence.dative_bonds().right_count(),
            correspondence
                .dative_bonds()
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index())),
        )?);
        let aromatic_systems = AromaticSystems::new(reorder(
            self.aromatic_systems.remap(&graph_remapping).into_entries(),
            correspondence.aromatic_systems().right_count(),
            correspondence
                .aromatic_systems()
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index())),
        )?);
        let multicenter_bonds = MulticenterBonds::new(reorder(
            self.multicenter_bonds
                .remap(&graph_remapping)
                .into_entries(),
            correspondence.multicenter_bonds().right_count(),
            correspondence
                .multicenter_bonds()
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index())),
        )?);
        let noncovalent_bonds = FixedRelationSet::new(reorder(
            self.noncovalent_bonds
                .remap(&graph_remapping)
                .into_entries(),
            correspondence.noncovalent_bonds().right_count(),
            correspondence
                .noncovalent_bonds()
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index())),
        )?);
        let stereo_atoms = FixedVarBirelationSet::new(reorder(
            self.stereo_atoms.remap(&graph_remapping).into_entries(),
            correspondence.stereo_atoms().right_count(),
            correspondence
                .stereo_atoms()
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index())),
        )?);
        let stereo_bonds = FixedVarBirelationSet::new(reorder(
            self.stereo_bonds.remap(&graph_remapping).into_entries(),
            correspondence.stereo_bonds().right_count(),
            correspondence
                .stereo_bonds()
                .matched_pairs()
                .iter()
                .map(|&(left, right)| (left.index(), right.index())),
        )?);
        let constraints: Constraints = self
            .constraints
            .clone()
            .into_iter()
            .map(|constraint| constraint.remap(&id_remapping))
            .collect();

        Some(Self::from_arcs(
            graph,
            Arc::new(atoms),
            Arc::new(bonds),
            dative_bonds,
            aromatic_systems,
            multicenter_bonds,
            Arc::new(noncovalent_bonds).into(),
            Arc::new(stereo_atoms).into(),
            Arc::new(stereo_bonds).into(),
            constraints,
        ))
    }
}

fn reorder<T>(
    values: Vec<T>,
    target_count: usize,
    pairs: impl IntoIterator<Item = (usize, usize)>,
) -> Option<Vec<T>> {
    let mut source = values.into_iter().map(Some).collect::<Vec<_>>();
    let mut target = (0..target_count).map(|_| None).collect::<Vec<_>>();
    for (left, right) in pairs {
        let value = source.get_mut(left)?.take()?;
        let slot = target.get_mut(right)?;
        if slot.replace(value).is_some() {
            return None;
        }
    }
    if source.iter().any(Option::is_some) {
        return None;
    }
    target.into_iter().collect()
}
