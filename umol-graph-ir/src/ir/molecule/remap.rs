//! Molecule remapping operations.

use std::sync::Arc;

use umol_graph_core::Graph;

use super::super::aromatic::AromaticSystems;
use super::super::dative::DativeBonds;
use super::super::multicenter::MulticenterBonds;
use super::super::noncovalent::NoncovalentBonds;
use super::super::stereo::{StereoAtoms, StereoBonds};
use super::Molecule;
use crate::ir::{Constraints, MoleculeCorrespondence, MoleculeRemapping};

impl Molecule {
    /// Renumber every entity table and reference using the supplied permutations.
    ///
    /// Transports topology, relation participants, entity forms, and constraint references.
    /// Participant sequences and their positional payloads are preserved. This does not validate
    /// chemistry, normalize attributes, repair references, or add or remove entities.
    ///
    /// # Panics
    ///
    /// Panics when any component length differs from this molecule's corresponding entity count.
    /// Use [`Self::try_remap`] for independently supplied permutations.
    ///
    /// # Semantic properties
    ///
    /// Identity is exact; applying the inverse permutations recovers the original molecule.
    /// Sequential renumbering agrees with the composition of the component permutations.
    /// The result satisfies [`Self::framed_eq_under`] under the same remapping.
    pub fn remap(&self, remapping: &MoleculeRemapping) -> Self {
        self.try_remap(remapping)
            .expect("molecule remapping requires matching entity counts")
    }

    /// Checked form of [`Self::remap`].
    ///
    /// Returns `None` when any component length differs from the molecule's entity count.
    pub fn try_remap(&self, remapping: &MoleculeRemapping) -> Option<Self> {
        let counts_match = [
            (remapping.graph().nodes().len(), self.atoms.len()),
            (remapping.graph().edges().len(), self.bonds.len()),
            (remapping.dative_bonds().len(), self.dative_bonds.count()),
            (
                remapping.aromatic_systems().len(),
                self.aromatic_systems.count(),
            ),
            (
                remapping.multicenter_bonds().len(),
                self.multicenter_bonds.count(),
            ),
            (
                remapping.noncovalent_bonds().len(),
                self.noncovalent_bonds.count(),
            ),
            (remapping.stereo_atoms().len(), self.stereo_atoms.count()),
            (remapping.stereo_bonds().len(), self.stereo_bonds.count()),
        ]
        .into_iter()
        .all(|(mapped, actual)| mapped == actual);
        if !counts_match {
            return None;
        }

        let graph_remapping = remapping.graph();

        let atoms = remapping
            .graph()
            .nodes()
            .remap_vec(self.atoms.as_ref().clone());
        let bonds = remapping
            .graph()
            .edges()
            .remap_vec(self.bonds.as_ref().clone());
        let edges = remapping.graph().edges().remap_vec(
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
        );
        let graph = Graph::new(atoms.len(), &edges);

        let dative_bonds = DativeBonds::new(
            remapping
                .dative_bonds()
                .remap_vec(self.dative_bonds.remap(graph_remapping).into_entries()),
        );
        let aromatic_systems = AromaticSystems::new(
            remapping
                .aromatic_systems()
                .remap_vec(self.aromatic_systems.remap(graph_remapping).into_entries()),
        );
        let multicenter_bonds = MulticenterBonds::new(
            remapping
                .multicenter_bonds()
                .remap_vec(self.multicenter_bonds.remap(graph_remapping).into_entries()),
        );
        let noncovalent_bonds = NoncovalentBonds::new(
            remapping
                .noncovalent_bonds()
                .remap_vec(self.noncovalent_bonds.remap(graph_remapping).into_entries()),
        );
        let stereo_atoms = StereoAtoms::new(
            remapping
                .stereo_atoms()
                .remap_vec(self.stereo_atoms.remap(graph_remapping).into_entries()),
        );
        let stereo_bonds = StereoBonds::new(
            remapping
                .stereo_bonds()
                .remap_vec(self.stereo_bonds.remap(graph_remapping).into_entries()),
        );
        let correspondence = MoleculeCorrespondence::from(remapping);
        let constraints: Constraints = self
            .constraints
            .clone()
            .into_iter()
            .map(|constraint| constraint.map(&correspondence))
            .collect();

        Some(
            Self::try_from_arcs(
                graph,
                Arc::new(atoms),
                Arc::new(bonds),
                dative_bonds,
                aromatic_systems,
                multicenter_bonds,
                noncovalent_bonds,
                stereo_atoms,
                stereo_bonds,
                constraints,
            )
            .expect("dense molecule remapping preserves representation integrity"),
        )
    }
}
