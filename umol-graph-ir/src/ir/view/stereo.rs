//! Stereo atom and stereo bond views.

use std::collections::HashSet;
use std::iter;

use umol_perm::{OrientedPermutationGroup, Permutation};

use super::super::id::{AtomId, BondId, StereoAtomId, StereoBondId, StereoLigandPosition};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::molecule::Molecule;
use super::super::stereo::{
    StereoAtomForm, StereoAtoms, StereoBondForm, StereoBonds, StereoKind, Stereogenicity, Topicity,
};
use super::super::symmetry::StereoSymmetry;
use super::super::traits::{FrameTransport, Lattice};
use super::atom::AtomView;
use super::bond::BondView;
use super::constraints::{StereoAtomConstraintsView, StereoBondConstraintsView};
use super::ligand::StereoLigandView;
use crate::ir::{
    StereoAtomConstraintForm, StereoAtomConstraintKey, StereoAtomConstraintsForm,
    StereoBondConstraintForm, StereoBondConstraintKey, StereoBondConstraintsForm, StereoCoset,
};

/// Namespace accessor for stereo-atom views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct StereoAtomViews<'a> {
    molecule: &'a Molecule,
}

impl<'a> StereoAtomViews<'a> {
    /// Ids of stereo atoms incident on `atom` (site or ligand).
    pub fn incident_ids(&self, atom: AtomId) -> impl ExactSizeIterator<Item = StereoAtomId> + 'a {
        self.molecule.raw_stereo_atoms().incident_ids(atom)
    }

    /// Id of the stereo atom on `site` with exactly this ligand set, if any. The frame order is not
    /// matched.
    pub fn of_id(&self, site: AtomId, ligands: &[StereoLigand]) -> Option<StereoAtomId> {
        self.molecule
            .raw_stereo_atoms()
            .coincident_id(site, ligands)
    }

    /// Any stereo atom is incident on `atom` (site or ligand).
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.molecule.raw_stereo_atoms().has_incident(atom)
    }

    /// Views of stereo atoms incident on `atom` (site or ligand).
    pub fn incident(&self, atom: AtomId) -> impl ExactSizeIterator<Item = StereoAtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.incident_ids(atom)
            .map(move |id| StereoAtomView { molecule, id })
    }

    // Ids of stereo atoms incident, in which `atom` is ligand.
    pub fn incident_as_ligand_ids(&self, atom: AtomId) -> impl Iterator<Item = StereoAtomId> + 'a {
        let set = self.molecule.raw_stereo_atoms();
        let ligand = StereoLigand {
            atom_id: atom,
            kind: StereoLigandKind::Atom,
        };
        set.incident_ids(atom)
            .filter(move |&id| set.ligands(id).contains(&ligand))
    }

    /// Any stereo atom is incident, in which `atom` is ligand.
    pub fn has_incident_as_ligand(&self, atom: AtomId) -> bool {
        let set = self.molecule.raw_stereo_atoms();
        let ligand = StereoLigand {
            atom_id: atom,
            kind: StereoLigandKind::Atom,
        };
        set.incident_ids(atom)
            .any(|id| set.ligands(id).contains(&ligand))
    }

    /// Views of stereo atoms incident, in which `atom` is ligand.
    pub fn incident_as_ligand(
        &self,
        atom: AtomId,
    ) -> impl Iterator<Item = StereoAtomView<'a>> + 'a {
        let molecule = self.molecule;
        self.incident_as_ligand_ids(atom)
            .map(move |id| StereoAtomView { molecule, id })
    }

    /// Id of the stereo atom sited on `atom`, if any.
    pub fn at_id(&self, atom: AtomId) -> Option<StereoAtomId> {
        let set = self.molecule.raw_stereo_atoms();
        set.incident_ids(atom).find(move |&id| set.site(id) == atom)
    }

    /// Whether a stereo atom is sited on `atom`.
    pub fn is_at(&self, atom: AtomId) -> bool {
        let set = self.molecule.raw_stereo_atoms();
        set.incident_ids(atom).any(move |id| set.site(id) == atom)
    }

    /// View of the stereo atom sited on `atom`, if any.
    pub fn at(&self, atom: AtomId) -> Option<StereoAtomView<'a>> {
        let molecule = self.molecule;
        self.at_id(atom)
            .map(move |id| StereoAtomView { molecule, id })
    }
}

/// Borrowed view of a stereo atom: the site atom, its ordered ligands, and data.
#[derive(Clone, Copy, Debug)]
pub struct StereoAtomView<'a> {
    molecule: &'a Molecule,
    id: StereoAtomId,
}

impl<'a> StereoAtomView<'a> {
    /// Constraint reading of this stereo atom: the container's read API
    /// (asserted side, meanings intact) plus the keyed accessors. Mutation
    /// stays on the stored container.
    #[inline]
    pub fn constraints(&self) -> StereoAtomConstraintsView<'a> {
        StereoAtomConstraintsView::new(self.molecule, self.id)
    }

    /// View of the stereo site atom.
    pub fn site(&self) -> AtomView<'a> {
        self.molecule.atom(self.site_id())
    }

    /// Site atom followed by the distinct ligand atoms — the relation's atom
    /// incidence. Deduped: a virtual ligand's bearing atom is the site, so it is
    /// not repeated.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let site = self.site_id();
        let ligands = self.molecule.raw_stereo_atoms().ligands(self.id);
        let mut seen = HashSet::new();
        iter::once(site)
            .chain(ligands.iter().map(|l| l.atom_id))
            .filter(move |id| seen.insert(*id))
    }
}

/// Read-only editor access to a stereo atom.
pub struct StereoAtomEditorView<'a> {
    stereo_atoms: &'a StereoAtoms,
    id: StereoAtomId,
}

/// Mutable attribute access to a stereo atom.
#[derive(Debug)]
pub struct StereoAtomViewMut<'a> {
    stereo_atoms: &'a mut StereoAtoms,
    id: StereoAtomId,
}

/// Mutable editor access to a stereo atom.
///
/// Structural mutations update incidence and preserve attributes and constraints.
/// Ligand mutations preserve the site. Publication checks structural integrity.
#[derive(Debug)]
pub struct StereoAtomEditorViewMut<'a> {
    stereo_atoms: &'a mut StereoAtoms,
    id: StereoAtomId,
}

/// Namespace accessor for stereo-bond views on a `Molecule`.
#[derive(Clone, Copy)]
pub struct StereoBondViews<'a> {
    molecule: &'a Molecule,
}

impl<'a> StereoBondViews<'a> {
    /// Id of the stereo bond on `site` with exactly this ligand set, if any. The frame order is not
    /// matched.
    pub fn of_id(&self, site: BondId, ligands: &[StereoLigand]) -> Option<StereoBondId> {
        self.molecule
            .raw_stereo_bonds()
            .coincident_id(site, ligands)
    }

    /// Ids of stereo bonds incident on `atom` (site endpoint or ligand). The
    /// site is an edge, so node incidence covers only ligands; site-endpoint
    /// membership is unioned in (and deduped) explicitly.
    pub fn incident_to_atom_ids(&self, atom: AtomId) -> impl Iterator<Item = StereoBondId> + 'a {
        let ligand_ids = self.molecule.raw_stereo_bonds().incident_to_atom_ids(atom);
        let mut seen = HashSet::new();
        self.incident_as_site_ids(atom)
            .chain(ligand_ids)
            .filter(move |id| seen.insert(*id))
    }

    /// Any stereo bond is incident on `atom` (site endpoint or ligand).
    pub fn has_incident_to_atom(&self, atom: AtomId) -> bool {
        self.molecule.raw_stereo_bonds().has_incident_to_atom(atom)
            || self.has_incident_as_site(atom)
    }

    /// Views of stereo bonds incident on `atom` (site endpoint or ligand anchor).
    pub fn incident_to_atom(&self, atom: AtomId) -> impl Iterator<Item = StereoBondView<'a>> + 'a {
        let molecule = self.molecule;
        self.incident_to_atom_ids(atom)
            .map(move |id| StereoBondView { molecule, id })
    }

    /// Ids of stereo bonds incident, in which `atom` is ligand.
    pub fn incident_as_ligand_ids(&self, atom: AtomId) -> impl Iterator<Item = StereoBondId> + 'a {
        let set = self.molecule.raw_stereo_bonds();
        let ligand = StereoLigand {
            atom_id: atom,
            kind: StereoLigandKind::Atom,
        };
        set.incident_to_atom_ids(atom)
            .filter(move |&id| set.ligands(id).contains(&ligand))
    }

    /// Any stereo bond is incident, in which `atom` is ligand.
    pub fn has_incident_as_ligand(&self, atom: AtomId) -> bool {
        let set = self.molecule.raw_stereo_bonds();
        let ligand = StereoLigand {
            atom_id: atom,
            kind: StereoLigandKind::Atom,
        };
        set.incident_to_atom_ids(atom)
            .any(|id| set.ligands(id).contains(&ligand))
    }

    /// Views of stereo bonds incident, in which `atom` is ligand.
    pub fn incident_as_ligand(
        &self,
        atom: AtomId,
    ) -> impl Iterator<Item = StereoBondView<'a>> + 'a {
        let molecule = self.molecule;
        self.incident_as_ligand_ids(atom)
            .map(move |id| StereoBondView { molecule, id })
    }

    /// Ids of stereo bonds, in which `atom` is a site endpoint.
    pub fn incident_as_site_ids(&self, atom: AtomId) -> impl Iterator<Item = StereoBondId> + 'a {
        let set = self.molecule.raw_stereo_bonds();
        self.molecule
            .neighbors(atom)
            .flat_map(move |n| set.incident_to_bond_ids(n.bond_id()))
    }

    /// Any stereo bond, in which `atom` is a site endpoint.
    pub fn has_incident_as_site(&self, atom: AtomId) -> bool {
        let set = self.molecule.raw_stereo_bonds();
        self.molecule
            .neighbors(atom)
            .any(move |n| set.has_incident_to_bond(n.bond_id()))
    }

    /// Views of stereo bonds, in which `atom` is a site endpoint.
    pub fn incident_as_site(&self, atom: AtomId) -> impl Iterator<Item = StereoBondView<'a>> + 'a {
        let molecule = self.molecule;
        self.incident_as_site_ids(atom)
            .map(move |id| StereoBondView { molecule, id })
    }

    /// Id of the stereo bond sited on `bond`, if any.
    pub fn at_id(&self, bond: BondId) -> Option<StereoBondId> {
        self.molecule
            .raw_stereo_bonds()
            .incident_to_bond_ids(bond)
            .next()
    }

    /// Whether a stereo bond is sited on `bond`.
    pub fn is_at(&self, bond: BondId) -> bool {
        self.molecule.raw_stereo_bonds().has_incident_to_bond(bond)
    }

    /// View of the stereo bond sited on `bond`, if any.
    pub fn at(&self, bond: BondId) -> Option<StereoBondView<'a>> {
        let molecule = self.molecule;
        self.at_id(bond)
            .map(move |id| StereoBondView { molecule, id })
    }
}

/// Borrowed view of a stereo bond: the site bond, its ordered ligands, and data.
#[derive(Clone, Copy, Debug)]
pub struct StereoBondView<'a> {
    molecule: &'a Molecule,
    id: StereoBondId,
}

impl<'a> StereoBondView<'a> {
    /// Constraint reading of this stereo bond: the container's read API
    /// (asserted side, meanings intact) plus the keyed accessors. Mutation
    /// stays on the stored container.
    #[inline]
    pub fn constraints(&self) -> StereoBondConstraintsView<'a> {
        StereoBondConstraintsView::new(self.molecule, self.id)
    }

    /// View of the stereo site bond.
    pub fn site(&self) -> BondView<'a> {
        self.molecule.bond(self.site_id())
    }

    /// The site bond's two atoms followed by the distinct ligand atoms — the
    /// relation's atom incidence. Deduped: a virtual ligand's bearing atom is a
    /// site endpoint, so it is not repeated.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let [a, b] = self.site().atom_ids();
        let ligands = self.molecule.raw_stereo_bonds().ligands(self.id);
        let mut seen = HashSet::new();
        [a, b]
            .into_iter()
            .chain(ligands.iter().map(|l| l.atom_id))
            .filter(move |id| seen.insert(*id))
    }
}

/// Read-only editor access to a stereo bond.
pub struct StereoBondEditorView<'a> {
    stereo_bonds: &'a StereoBonds,
    id: StereoBondId,
}

/// Mutable attribute access to a stereo bond.
#[derive(Debug)]
pub struct StereoBondViewMut<'a> {
    stereo_bonds: &'a mut StereoBonds,
    id: StereoBondId,
}

/// Mutable editor access to a stereo bond.
///
/// Structural mutations update incidence and preserve attributes and constraints.
/// Ligand mutations preserve the site. Publication checks structural integrity.
#[derive(Debug)]
pub struct StereoBondEditorViewMut<'a> {
    stereo_bonds: &'a mut StereoBonds,
    id: StereoBondId,
}

macro_rules! stereo_views {
    ($views:ident, $view:ident, $id:ty, $raw:ident) => {
        impl<'a> $views<'a> {
            pub(crate) fn new(molecule: &'a Molecule) -> Self {
                Self { molecule }
            }

            pub fn count(&self) -> usize {
                self.molecule.$raw().count()
            }

            pub fn ids(&self) -> impl ExactSizeIterator<Item = $id> {
                self.molecule.$raw().ids()
            }

            pub fn iter(&self) -> impl ExactSizeIterator<Item = $view<'a>> {
                let molecule = self.molecule;
                molecule.$raw().ids().map(move |id| $view { molecule, id })
            }

            pub fn contains(&self, id: $id) -> bool {
                self.molecule.$raw().contains(id)
            }

            pub fn get(&self, id: $id) -> Option<$view<'a>> {
                self.contains(id).then_some($view {
                    molecule: self.molecule,
                    id,
                })
            }
        }
    };
}

stereo_views!(
    StereoAtomViews,
    StereoAtomView,
    StereoAtomId,
    raw_stereo_atoms
);
stereo_views!(
    StereoBondViews,
    StereoBondView,
    StereoBondId,
    raw_stereo_bonds
);

macro_rules! stereo_editor_view {
    ($view:ident, $set:ident, $sets:ty, $id:ty, $site:ty, $form:ty) => {
        impl<'a> $view<'a> {
            pub(crate) fn new($set: &'a $sets, id: $id) -> Self {
                Self { $set, id }
            }

            #[inline]
            pub fn id(&self) -> $id {
                self.id
            }

            #[inline]
            pub fn attributes(&self) -> &'a $form {
                self.$set.attributes(self.id)
            }

            #[inline]
            pub fn site_id(&self) -> $site {
                self.$set.site(self.id)
            }

            /// The ordered ligand identifiers.
            #[inline]
            pub fn ligand_ids(&self) -> &'a [StereoLigand] {
                self.$set.ligands(self.id)
            }
        }
    };
}

stereo_editor_view!(
    StereoAtomEditorView,
    stereo_atoms,
    StereoAtoms,
    StereoAtomId,
    AtomId,
    StereoAtomForm
);
stereo_editor_view!(
    StereoBondEditorView,
    stereo_bonds,
    StereoBonds,
    StereoBondId,
    BondId,
    StereoBondForm
);

macro_rules! stereo_view_mut {
    ($view:ident, $set:ident, $sets:ty, $id:ty, $site:ty, $form:ty, $constraints:ty) => {
        impl<'a> $view<'a> {
            pub(crate) fn new($set: &'a mut $sets, id: $id) -> Self {
                Self { $set, id }
            }

            #[inline]
            pub fn id(&self) -> $id {
                self.id
            }

            #[inline]
            pub fn attributes(&self) -> &$form {
                self.$set.attributes(self.id)
            }

            #[inline]
            pub fn attributes_mut(&mut self) -> &mut $form {
                self.$set.attributes_mut(self.id)
            }

            #[inline]
            pub fn site_id(&self) -> $site {
                self.$set.site(self.id)
            }

            /// The ordered ligand identifiers.
            #[inline]
            pub fn ligand_ids(&self) -> &[StereoLigand] {
                self.$set.ligands(self.id)
            }

            #[inline]
            pub fn constraints(&self) -> &$constraints {
                &self.attributes().constraints
            }
        }
    };
}

stereo_view_mut!(
    StereoAtomViewMut,
    stereo_atoms,
    StereoAtoms,
    StereoAtomId,
    AtomId,
    StereoAtomForm,
    StereoAtomConstraintsForm
);
stereo_view_mut!(
    StereoAtomEditorViewMut,
    stereo_atoms,
    StereoAtoms,
    StereoAtomId,
    AtomId,
    StereoAtomForm,
    StereoAtomConstraintsForm
);
stereo_view_mut!(
    StereoBondViewMut,
    stereo_bonds,
    StereoBonds,
    StereoBondId,
    BondId,
    StereoBondForm,
    StereoBondConstraintsForm
);
stereo_view_mut!(
    StereoBondEditorViewMut,
    stereo_bonds,
    StereoBonds,
    StereoBondId,
    BondId,
    StereoBondForm,
    StereoBondConstraintsForm
);

/// Methods shared by the readonly molecule stereo views.
macro_rules! stereo_view_queries {
    ($view:ident, $id:ty, $form:ty, $site:ty, $raw:ident, $ligand_error:literal) => {
        impl<'a> $view<'a> {
            #[inline]
            pub fn id(&self) -> $id {
                self.id
            }

            #[inline]
            pub fn attributes(&self) -> &'a $form {
                self.molecule.$raw().attributes(self.id)
            }

            /// Id of the stereo site.
            #[inline]
            pub fn site_id(&self) -> $site {
                self.molecule.$raw().site(self.id)
            }

            /// The ordered ligands defining the stereo configuration.
            pub fn ligands(&self) -> impl ExactSizeIterator<Item = StereoLigandView<'a>> + 'a {
                let molecule = self.molecule;
                let ligands = self.molecule.$raw().ligands(self.id);
                ligands
                    .iter()
                    .map(move |ligand| StereoLigandView::new(*ligand, molecule))
            }

            /// View of the ligand at the given coordination position. Panics if it is
            /// not a position in this stereo frame.
            pub fn ligand(&self, ligand_id: StereoLigandPosition) -> StereoLigandView<'a> {
                let ligand = *self
                    .molecule
                    .$raw()
                    .ligands(self.id)
                    .get(ligand_id.index())
                    .expect($ligand_error);
                StereoLigandView::new(ligand, self.molecule)
            }

            pub fn atom_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
                self.ligands()
                    .filter(|ligand| ligand.kind() == StereoLigandKind::Atom)
            }

            pub fn implicit_hydrogen_ligands(
                &self,
            ) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
                self.ligands()
                    .filter(|ligand| ligand.kind() == StereoLigandKind::ImplicitHydrogen)
            }

            pub fn lone_pair_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
                self.ligands()
                    .filter(|ligand| ligand.kind() == StereoLigandKind::LonePair)
            }

            /// The stored configuration read against `ligands`, which must hold the same distinct
            /// participants as the stored frame.
            pub fn coset_for(
                &self,
                ligands: impl IntoIterator<Item = StereoLigand>,
            ) -> Option<StereoCoset> {
                let requested: Vec<StereoLigand> = ligands.into_iter().collect();
                let action =
                    Permutation::between(self.molecule.$raw().ligands(self.id), &requested)?;
                self.attributes()
                    .configuration
                    .clone()
                    .reframe_by(&action)?
                    .coset()
                    .cloned()
            }

            pub fn is_ground(&self) -> bool {
                self.attributes().is_ground()
            }

            /// The local oriented ligand-position symmetry group.
            pub fn ligand_symmetry(&self, symmetry: &StereoSymmetry) -> OrientedPermutationGroup {
                symmetry.group().clone()
            }

            /// Stereogenicity classification of this carrier.
            pub fn stereogenicity(&self, symmetry: &StereoSymmetry) -> Stereogenicity {
                symmetry.stereogenicity()
            }

            /// Whether this carrier is a genuine stereocenter.
            pub fn is_stereogenic(&self, symmetry: &StereoSymmetry) -> bool {
                symmetry.is_stereogenic()
            }

            /// Whether this carrier is prochiral (some enantiotopic ligand pair).
            pub fn is_prochiral(&self, symmetry: &StereoSymmetry) -> bool {
                symmetry.stereogenicity() == Stereogenicity::Prochiral
            }

            /// Kind-level chirality — whether the geometry can encode handedness.
            /// No symmetry computation.
            pub fn is_chiral(&self) -> bool {
                self.kind().is_chiral_class()
            }

            /// Topicity of two ligand positions.
            pub fn topicity(
                &self,
                a: StereoLigandPosition,
                b: StereoLigandPosition,
                symmetry: &StereoSymmetry,
            ) -> Topicity {
                symmetry.topicity(a, b)
            }

            pub fn is_homotopic(
                &self,
                a: StereoLigandPosition,
                b: StereoLigandPosition,
                symmetry: &StereoSymmetry,
            ) -> bool {
                symmetry.topicity(a, b) == Topicity::Homotopic
            }

            pub fn is_enantiotopic(
                &self,
                a: StereoLigandPosition,
                b: StereoLigandPosition,
                symmetry: &StereoSymmetry,
            ) -> bool {
                symmetry.topicity(a, b) == Topicity::Enantiotopic
            }

            pub fn is_diastereotopic(
                &self,
                a: StereoLigandPosition,
                b: StereoLigandPosition,
                symmetry: &StereoSymmetry,
            ) -> bool {
                symmetry.topicity(a, b) == Topicity::Diastereotopic
            }

            /// The ordered ligand identifiers.
            #[inline]
            pub fn ligand_ids(&self) -> &'a [StereoLigand] {
                self.molecule.$raw().ligands(self.id)
            }
        }
    };
}

stereo_view_queries!(
    StereoAtomView,
    StereoAtomId,
    StereoAtomForm,
    AtomId,
    raw_stereo_atoms,
    "ligand id must refer to a ligand of this stereo atom"
);
stereo_view_queries!(
    StereoBondView,
    StereoBondId,
    StereoBondForm,
    BondId,
    raw_stereo_bonds,
    "ligand id must refer to a ligand of this stereo bond"
);

macro_rules! stereo_local_queries {
    ($view:ident, $borrow:lifetime) => {
        impl<'a> $view<'a> {
            /// The coordination-geometry kind.
            ///
            /// # Panics
            ///
            /// Panics when the configuration is undetermined.
            #[inline]
            pub fn kind(&self) -> StereoKind {
                self.attributes()
                    .configuration
                    .kind()
                    .expect("stereo view has a concrete kind")
            }

            /// The stereo coset.
            ///
            /// # Panics
            ///
            /// Panics when the configuration is undetermined.
            #[inline]
            pub fn coset(&self) -> &$borrow StereoCoset {
                self.attributes()
                    .configuration
                    .coset()
                    .expect("stereo view has a concrete coset")
            }

            pub fn ligand_count(&self) -> usize {
                self.ligand_ids().len()
            }

            /// The frame position of an actual-atom ligand with this id, if present.
            pub fn ligand_position(&self, id: AtomId) -> Option<StereoLigandPosition> {
                self.ligand_ids()
                    .iter()
                    .position(|ligand| {
                        ligand.kind == StereoLigandKind::Atom && ligand.atom_id == id
                    })
                    .map(|index| StereoLigandPosition(index as u32))
            }

            pub fn atom_ligand_ids(&self) -> impl Iterator<Item = AtomId> + $borrow {
                self.ligand_ids()
                    .iter()
                    .filter(|ligand| ligand.kind == StereoLigandKind::Atom)
                    .map(|ligand| ligand.atom_id)
            }

            pub fn implicit_hydrogen_atom_ids(&self) -> impl Iterator<Item = AtomId> + $borrow {
                self.ligand_ids()
                    .iter()
                    .filter(|ligand| ligand.kind == StereoLigandKind::ImplicitHydrogen)
                    .map(|ligand| ligand.atom_id)
            }

            pub fn lone_pair_atom_ids(&self) -> impl Iterator<Item = AtomId> + $borrow {
                self.ligand_ids()
                    .iter()
                    .filter(|ligand| ligand.kind == StereoLigandKind::LonePair)
                    .map(|ligand| ligand.atom_id)
            }

            pub fn atom_ligand_count(&self) -> usize {
                self.atom_ligand_ids().count()
            }

            pub fn implicit_hydrogen_count(&self) -> usize {
                self.implicit_hydrogen_atom_ids().count()
            }

            pub fn lone_pair_count(&self) -> usize {
                self.lone_pair_atom_ids().count()
            }
        }
    };
}

stereo_local_queries!(StereoAtomView, 'a);
stereo_local_queries!(StereoAtomEditorView, 'a);
stereo_local_queries!(StereoAtomViewMut, '_);
stereo_local_queries!(StereoAtomEditorViewMut, '_);
stereo_local_queries!(StereoBondView, 'a);
stereo_local_queries!(StereoBondEditorView, 'a);
stereo_local_queries!(StereoBondViewMut, '_);
stereo_local_queries!(StereoBondEditorViewMut, '_);

macro_rules! stereo_editor_mutation {
    ($view:ident, $set:ident, $site:ty) => {
        impl<'a> $view<'a> {
            /// Replace the site, preserving ligands, attributes, and constraints.
            pub fn replace_site(&mut self, site: $site) {
                self.$set.replace_site(self.id, site);
            }

            /// Replace the ligand list, preserving the supplied order.
            pub fn replace_ligands(&mut self, ligands: &[StereoLigand]) {
                self.$set.replace_ligands(self.id, ligands);
            }

            /// Replace the ligand at `position` without changing the ligand count.
            ///
            /// # Panics
            ///
            /// Panics when `position` is outside the ligand list.
            pub fn replace_ligand(&mut self, position: StereoLigandPosition, ligand: StereoLigand) {
                self.$set.replace_ligand(self.id, position.index(), ligand);
            }

            /// Insert a ligand at `position`, preserving the order of the existing ligands.
            ///
            /// # Panics
            ///
            /// Panics when `position` exceeds the ligand count. Insertion at the end is allowed.
            pub fn insert_ligand(&mut self, position: StereoLigandPosition, ligand: StereoLigand) {
                self.$set.insert_ligand(self.id, position.index(), ligand);
            }

            /// Remove the ligand at `position`, preserving the order of the remaining ligands.
            ///
            /// # Panics
            ///
            /// Panics when `position` is outside the ligand list.
            pub fn remove_ligand(&mut self, position: StereoLigandPosition) {
                self.$set.remove_ligand(self.id, position.index());
            }
        }
    };
}

stereo_editor_mutation!(StereoAtomEditorViewMut, stereo_atoms, AtomId);
stereo_editor_mutation!(StereoBondEditorViewMut, stereo_bonds, BondId);

/// Stored constraint container of the stereo atom `id`.
pub(crate) fn stereo_atom_asserted_constraints(
    molecule: &Molecule,
    id: StereoAtomId,
) -> &StereoAtomConstraintsForm {
    &molecule.stereo_atom(id).attributes().constraints
}

/// Derived side of one stereo-atom constraint key: always vacuous.
pub(crate) fn stereo_atom_derived_constraint(
    _molecule: &Molecule,
    _id: StereoAtomId,
    _key: StereoAtomConstraintKey,
    _complete: bool,
) -> Option<StereoAtomConstraintForm> {
    None
}

/// Stored constraint container of the stereo bond `id`.
pub(crate) fn stereo_bond_asserted_constraints(
    molecule: &Molecule,
    id: StereoBondId,
) -> &StereoBondConstraintsForm {
    &molecule.stereo_bond(id).attributes().constraints
}

/// Derived side of one stereo-bond constraint key: always vacuous.
pub(crate) fn stereo_bond_derived_constraint(
    _molecule: &Molecule,
    _id: StereoBondId,
    _key: StereoBondConstraintKey,
    _complete: bool,
) -> Option<StereoBondConstraintForm> {
    None
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_core::AutomorphismAlgorithm;

    use super::super::assert_exact_size_by;
    use crate::ir::atom::AtomForm;
    use crate::ir::bond::BondForm;
    use crate::ir::coloring::ConstitutionColoring;
    use crate::ir::constraint::{
        StereoAtomConstraintForm, StereoBondConstraintForm, StereogenicityForm,
    };
    use crate::ir::entity::Entity;
    use crate::ir::id::{AtomId, BondId, StereoAtomId, StereoBondId, StereoLigandPosition};
    use crate::ir::ligand::{StereoLigand, StereoLigandKind};
    use crate::ir::molecule::{Molecule, MoleculeEntries, MoleculeIntegrityError};
    use crate::ir::stereo::{
        StereoAtomForm, StereoBondForm, StereoCoset, StereoKind, Stereogenicity, Topicity,
    };
    use crate::ir::symmetry::GraphSymmetryConfig;

    #[fixture]
    fn molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            // Atom 6 is an unbonded spare: a node that is neither stereo-bond site nor ligand.
            atoms: vec![AtomForm::from_element(Element::C); 7],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(2)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
                (AtomId(0), AtomId(3), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(2), AtomId(4), BondForm::from_order(1)),
                (AtomId(2), AtomId(5), BondForm::from_order(1)),
                (AtomId(3), AtomId(1), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        })
    }

    #[fixture]
    fn virtual_ligand_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(2)),
                (AtomId(4), AtomId(5), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(2), AtomId(4), BondForm::from_order(1)),
                (AtomId(3), AtomId(5), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            stereo_bonds: vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        })
    }

    // A 4-membered ring (atoms 0-1-2-3) with two pendant atoms (4 on 0, 5 on 1);
    // a stereo atom on ring atom 0 and a stereo bond on ring bond 0-1.
    #[fixture]
    fn ring_molecule() -> Molecule {
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(1), AtomId(2), BondForm::from_order(1)),
                (AtomId(2), AtomId(3), BondForm::from_order(1)),
                (AtomId(3), AtomId(0), BondForm::from_order(1)),
                (AtomId(0), AtomId(4), BondForm::from_order(1)),
                (AtomId(1), AtomId(5), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            )],
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            )],
            ..Default::default()
        })
    }

    #[fixture]
    fn stereo_entries() -> MoleculeEntries {
        let mut bonds = vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(5), AtomId(6), BondForm::from_order(2)),
            (AtomId(0), AtomId(5), BondForm::from_order(1)),
        ];
        for site in [AtomId(0), AtomId(1), AtomId(5), AtomId(6)] {
            for ligand in [AtomId(2), AtomId(3), AtomId(4), AtomId(7)] {
                bonds.push((site, ligand, BondForm::from_order(1)));
            }
        }
        let ligands = vec![
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
        ];
        MoleculeEntries {
            atoms: vec![AtomForm::default(); 8],
            bonds,
            stereo_atoms: vec![
                (
                    AtomId(0),
                    ligands.clone(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))
                        .with_constraint(StereoAtomConstraintForm::Stereogenicity(
                            StereogenicityForm::Lit(Stereogenicity::Prochiral),
                        )),
                ),
                (AtomId(6), vec![], StereoAtomForm::default()),
            ],
            stereo_bonds: vec![
                (
                    BondId(0),
                    ligands.clone(),
                    StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)).with_constraint(
                        StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(
                            Stereogenicity::Prochiral,
                        )),
                    ),
                ),
                (BondId(2), ligands, StereoBondForm::default()),
            ],
            ..Default::default()
        }
    }

    #[rstest]
    fn test_stereo_atom_views_count(molecule: Molecule) {
        assert_eq!(molecule.stereo_atoms().count(), 1);
    }

    #[rstest]
    fn test_stereo_atom_views_ids(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().stereo_atoms().ids(), vec![], |id| id);
        assert_exact_size_by(molecule.stereo_atoms().ids(), vec![StereoAtomId(0)], |id| {
            id
        });
    }

    #[rstest]
    fn test_stereo_atom_views_iter(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().stereo_atoms().iter(), vec![], |view| {
            (view.id(), view.site_id())
        });
        assert_exact_size_by(
            molecule.stereo_atoms().iter(),
            vec![(StereoAtomId(0), AtomId(0))],
            |view| (view.id(), view.site_id()),
        );
    }

    #[rstest]
    #[case::present(StereoAtomId(0), true)]
    #[case::absent(StereoAtomId(99), false)]
    fn test_stereo_atom_views_contains(
        molecule: Molecule,
        #[case] id: StereoAtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.stereo_atoms().contains(id), expected);
    }

    #[rstest]
    fn test_stereo_atom_views_get(molecule: Molecule) {
        let res = molecule.stereo_atoms().get(StereoAtomId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id(), StereoAtomId(0));
        assert_eq!(view.site_id(), AtomId(0));
        assert_eq!(view.kind(), StereoKind::Tetrahedral);
        assert_eq!(
            view.ligands()
                .map(|ligand| (ligand.kind(), ligand.atom_id()))
                .collect::<Vec<_>>(),
            vec![
                (StereoLigandKind::Atom, AtomId(1)),
                (StereoLigandKind::Atom, AtomId(2)),
                (StereoLigandKind::Atom, AtomId(3)),
                (StereoLigandKind::Atom, AtomId(4)),
            ],
        );
        assert_eq!(
            view.attributes(),
            &StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
        );
    }

    #[rstest]
    fn test_stereo_atom_views_get_none(molecule: Molecule) {
        let res = molecule.stereo_atoms().get(StereoAtomId(99));
        assert!(res.is_none());
    }

    #[rstest]
    #[case::site(AtomId(0), vec![StereoAtomId(0)])]
    #[case::ligand(AtomId(2), vec![StereoAtomId(0)])]
    #[case::unrelated(AtomId(5), vec![])]
    fn test_stereo_atom_views_incident(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoAtomId>,
    ) {
        assert_exact_size_by(
            molecule.stereo_atoms().incident_ids(atom),
            expected.clone(),
            |id| id,
        );
        assert_exact_size_by(molecule.stereo_atoms().incident(atom), expected, |view| {
            view.id()
        });
    }

    #[rstest]
    #[case::ligand(AtomId(2), vec![StereoAtomId(0)])]
    #[case::site_not_ligand(AtomId(0), vec![])]
    #[case::unrelated(AtomId(5), vec![])]
    fn test_stereo_atom_views_incident_as_ligand(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoAtomId>,
    ) {
        assert_eq!(
            molecule
                .stereo_atoms()
                .incident_as_ligand_ids(atom)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::site(AtomId(0), Some(StereoAtomId(0)))]
    #[case::ligand_not_site(AtomId(2), None)]
    #[case::unrelated(AtomId(5), None)]
    fn test_stereo_atom_views_at(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Option<StereoAtomId>,
    ) {
        assert_eq!(molecule.stereo_atoms().at_id(atom), expected);
    }

    #[rstest]
    fn test_stereo_atom_view_site_id(molecule: Molecule) {
        assert_eq!(molecule.stereo_atom(StereoAtomId(0)).site_id(), AtomId(0));
    }

    #[rstest]
    fn test_stereo_atom_view_site(molecule: Molecule) {
        let view = molecule.stereo_atom(StereoAtomId(0)).site();
        assert_eq!(view.id(), AtomId(0));
        assert_eq!(view.attributes(), &AtomForm::from_element(Element::C));
    }

    #[rstest]
    fn test_stereo_atom_view_coset(molecule: Molecule) {
        assert_eq!(
            molecule.stereo_atom(StereoAtomId(0)).coset(),
            &StereoCoset::Lit(1),
        );
    }

    #[rstest]
    fn test_stereo_atom_view_ligand_count(molecule: Molecule) {
        assert_eq!(molecule.stereo_atom(StereoAtomId(0)).ligand_count(), 4);
    }

    #[rstest]
    fn test_stereo_atom_view_ligands(molecule: Molecule) {
        assert_exact_size_by(
            molecule.stereo_atom(StereoAtomId(0)).ligands(),
            vec![
                (StereoLigandKind::Atom, AtomId(1)),
                (StereoLigandKind::Atom, AtomId(2)),
                (StereoLigandKind::Atom, AtomId(3)),
                (StereoLigandKind::Atom, AtomId(4)),
            ],
            |ligand| (ligand.kind(), ligand.atom_id()),
        );
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), StereoLigandKind::Atom, AtomId(1))]
    #[case::last(StereoLigandPosition(3), StereoLigandKind::Atom, AtomId(4))]
    fn test_stereo_atom_view_ligand(
        molecule: Molecule,
        #[case] ligand_id: StereoLigandPosition,
        #[case] kind: StereoLigandKind,
        #[case] atom: AtomId,
    ) {
        let ligand = molecule.stereo_atom(StereoAtomId(0)).ligand(ligand_id);
        assert_eq!(ligand.kind(), kind);
        assert_eq!(ligand.atom_id(), atom);
    }

    #[rstest]
    fn test_stereo_atom_view_stereo_queries() {
        // A clean stereocenter: C bonded to four distinct halogens.
        let mol = Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
                AtomForm::from_element(Element::Br),
                AtomForm::from_element(Element::I),
            ],
            bonds: (1..=4)
                .map(|i| (AtomId(0), AtomId(i), BondForm::from_order(1)))
                .collect(),
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        });
        let gs = mol.graph_symmetry(&GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: true,
            max_iterations: 16,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        });
        let symmetry = mol.stereo_atom_symmetry(&gs, StereoAtomId(0));
        let view = mol.stereo_atom(StereoAtomId(0));

        assert!(view.is_chiral()); // tetrahedral kind, no symmetry computation
        assert!(view.is_stereogenic(&symmetry));
        assert_eq!(view.stereogenicity(&symmetry), Stereogenicity::Stereogenic);
        assert_eq!(
            view.topicity(StereoLigandPosition(0), StereoLigandPosition(1), &symmetry),
            Topicity::Diastereotopic,
        );
        assert!(view.is_diastereotopic(
            StereoLigandPosition(0),
            StereoLigandPosition(1),
            &symmetry
        ));
        assert_eq!(view.ligand_symmetry(&symmetry).order(), 1);
        assert_eq!(
            view.ligand_position(AtomId(2)),
            Some(StereoLigandPosition(1))
        );
        assert_eq!(view.ligand_position(AtomId(99)), None);
        assert_eq!(view.ligand_ids().len(), 4);
    }

    #[rstest]
    #[should_panic]
    fn test_stereo_atom_view_ligand_error(molecule: Molecule) {
        molecule
            .stereo_atom(StereoAtomId(0))
            .ligand(StereoLigandPosition(4));
    }

    #[rstest]
    fn test_stereo_ligand_view_atom(molecule: Molecule) {
        let ligand = molecule
            .stereo_atom(StereoAtomId(0))
            .ligands()
            .next()
            .unwrap();
        let atom = ligand.atom();
        assert_eq!(atom.id(), AtomId(1));
        assert_eq!(atom.attributes(), &AtomForm::from_element(Element::C));
    }

    #[rstest]
    fn test_stereo_atom_view_atom_ligands(virtual_ligand_molecule: Molecule) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .atom_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .atom_ligand_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(1), AtomId(4)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .atom_ligand_count(),
            2,
        );
    }

    #[rstest]
    fn test_stereo_atom_view_implicit_hydrogen_ligands(virtual_ligand_molecule: Molecule) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .implicit_hydrogen_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(0)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .implicit_hydrogen_atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .implicit_hydrogen_count(),
            1,
        );
    }

    #[rstest]
    fn test_stereo_atom_view_lone_pair_ligands(virtual_ligand_molecule: Molecule) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .lone_pair_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(0)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .lone_pair_atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .lone_pair_count(),
            1,
        );
    }

    #[rstest]
    fn test_stereo_atom_view_coset_for_none(molecule: Molecule) {
        let view = molecule.stereo_atom(StereoAtomId(0));
        let ligands = [
            StereoLigand {
                atom_id: AtomId(1),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(2),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(3),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(view.coset_for(ligands[..3].iter().copied()), None);
        assert_eq!(
            view.coset_for([ligands[0], ligands[0], ligands[2], ligands[3]]),
            None,
        );
        assert_eq!(
            view.coset_for([
                ligands[0],
                ligands[1],
                ligands[2],
                StereoLigand {
                    atom_id: AtomId(99),
                    kind: StereoLigandKind::Atom,
                },
            ]),
            None,
        );
    }

    #[rstest]
    fn test_stereo_atom_view_coset_for(molecule: Molecule) {
        let view = molecule.stereo_atom(StereoAtomId(0));
        let ligands = vec![
            StereoLigand {
                atom_id: AtomId(1),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(2),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(3),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(view.coset_for(ligands.clone()), Some(StereoCoset::Lit(1)),);

        let reordered = vec![ligands[1], ligands[0], ligands[2], ligands[3]];
        assert_eq!(view.coset_for(reordered), Some(StereoCoset::Lit(0)));
    }

    #[rstest]
    fn test_stereo_atom_view_atom_ids(molecule: Molecule, virtual_ligand_molecule: Molecule) {
        assert_eq!(
            molecule
                .stereo_atom(StereoAtomId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
        );
        // Virtual ligands carry the site atom, so it is not repeated.
        assert_eq!(
            virtual_ligand_molecule
                .stereo_atom(StereoAtomId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(0), AtomId(1), AtomId(4)],
        );
    }

    #[rstest]
    fn test_stereo_atom_view_attributes(molecule: Molecule) {
        let (attributes, frame) = {
            let view = molecule.stereo_atom(StereoAtomId(0));
            assert_eq!(view.id(), StereoAtomId(0));
            (view.attributes(), view.ligand_ids())
        };
        assert_eq!(
            attributes,
            &StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))
        );
        assert_eq!(
            frame,
            &[
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom)
            ]
        );
    }

    #[rstest]
    fn test_stereo_atom_editor_view_attributes(molecule: Molecule) {
        let editor = molecule.edit();
        let (attributes, frame) = {
            let view = editor.stereo_atom(StereoAtomId(0));
            assert_eq!(view.id(), StereoAtomId(0));
            (view.attributes(), view.ligand_ids())
        };
        assert_eq!(
            attributes,
            &StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))
        );
        assert_eq!(
            frame,
            &[
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom)
            ]
        );
    }

    #[rstest]
    fn test_stereo_atom_editor_view_local_queries(virtual_ligand_molecule: Molecule) {
        let editor = virtual_ligand_molecule.edit();
        let (coset, atoms, hydrogens, lone_pairs) = {
            let view = editor.stereo_atom(StereoAtomId(0));
            assert_eq!(view.kind(), StereoKind::Tetrahedral);
            assert_eq!(view.ligand_count(), 4);
            assert_eq!(view.atom_ligand_count(), 2);
            assert_eq!(view.implicit_hydrogen_count(), 1);
            assert_eq!(view.lone_pair_count(), 1);
            assert_eq!(
                view.ligand_position(AtomId(1)),
                Some(StereoLigandPosition(0))
            );
            assert_eq!(
                view.ligand_position(AtomId(4)),
                Some(StereoLigandPosition(3))
            );
            assert_eq!(view.ligand_position(AtomId(0)), None);
            assert_eq!(view.ligand_position(AtomId(99)), None);
            (
                view.coset(),
                view.atom_ligand_ids(),
                view.implicit_hydrogen_atom_ids(),
                view.lone_pair_atom_ids(),
            )
        };
        assert_eq!(coset, &StereoCoset::Lit(1));
        assert_eq!(atoms.collect::<Vec<_>>(), vec![AtomId(1), AtomId(4)]);
        assert_eq!(hydrogens.collect::<Vec<_>>(), vec![AtomId(0)]);
        assert_eq!(lone_pairs.collect::<Vec<_>>(), vec![AtomId(0)]);
    }

    #[rstest]
    #[case::kind(false)]
    #[case::coset(true)]
    #[should_panic(expected = "stereo view has a concrete")]
    fn test_stereo_atom_editor_view_configuration_error(
        mut stereo_entries: MoleculeEntries,
        #[case] coset: bool,
    ) {
        stereo_entries.stereo_atoms[0].2 = StereoAtomForm::default();
        let virtual_ligand_molecule = Molecule::from_entries(stereo_entries);
        let editor = virtual_ligand_molecule.edit();
        let view = editor.stereo_atom(StereoAtomId(0));
        if coset {
            let _ = view.coset();
        } else {
            let _ = view.kind();
        }
    }

    #[rstest]
    #[should_panic(expected = "invalid stereo atom id")]
    fn test_stereo_atom_editor_view_id_error(molecule: Molecule) {
        molecule.edit().stereo_atom(StereoAtomId(1));
    }

    #[rstest]
    fn test_stereo_atom_view_mut_local_queries(virtual_ligand_molecule: Molecule) {
        let mut molecule = virtual_ligand_molecule;
        let view = molecule.stereo_atom_mut(StereoAtomId(0));
        assert_eq!(view.kind(), StereoKind::Tetrahedral);
        assert_eq!(view.ligand_count(), 4);
        assert_eq!(view.atom_ligand_count(), 2);
        assert_eq!(view.implicit_hydrogen_count(), 1);
        assert_eq!(view.lone_pair_count(), 1);
        assert_eq!(
            view.ligand_position(AtomId(1)),
            Some(StereoLigandPosition(0))
        );
        assert_eq!(
            view.ligand_position(AtomId(4)),
            Some(StereoLigandPosition(3))
        );
        assert_eq!(view.ligand_position(AtomId(0)), None);
        assert_eq!(view.ligand_position(AtomId(99)), None);
        let (coset, atoms, hydrogens, lone_pairs) = (
            view.coset(),
            view.atom_ligand_ids(),
            view.implicit_hydrogen_atom_ids(),
            view.lone_pair_atom_ids(),
        );
        assert_eq!(coset, &StereoCoset::Lit(1));
        assert_eq!(atoms.collect::<Vec<_>>(), vec![AtomId(1), AtomId(4)]);
        assert_eq!(hydrogens.collect::<Vec<_>>(), vec![AtomId(0)]);
        assert_eq!(lone_pairs.collect::<Vec<_>>(), vec![AtomId(0)]);
    }

    #[rstest]
    #[case::kind(false)]
    #[case::coset(true)]
    #[should_panic(expected = "stereo view has a concrete")]
    fn test_stereo_atom_view_mut_configuration_error(
        mut stereo_entries: MoleculeEntries,
        #[case] coset: bool,
    ) {
        stereo_entries.stereo_atoms[0].2 = StereoAtomForm::default();
        let virtual_ligand_molecule = Molecule::from_entries(stereo_entries);
        let mut molecule = virtual_ligand_molecule;
        let view = molecule.stereo_atom_mut(StereoAtomId(0));
        if coset {
            let _ = view.coset();
        } else {
            let _ = view.kind();
        }
    }

    #[rstest]
    fn test_stereo_atom_view_mut_attributes_mut(mut molecule: Molecule) {
        let expected = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0));
        {
            let mut view = molecule.stereo_atom_mut(StereoAtomId(0));
            assert_eq!(view.id(), StereoAtomId(0));
            view.attributes_mut().configuration = expected.configuration.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.kind(), StereoKind::Tetrahedral);
            assert_eq!(view.coset(), &StereoCoset::Lit(0));
            assert_eq!(view.constraints(), &expected.constraints);
            assert_eq!(view.site_id(), AtomId(0));
            assert_eq!(
                view.ligand_ids(),
                &[
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom)
                ]
            );
        }
        assert_eq!(
            molecule.stereo_atom(StereoAtomId(0)).attributes(),
            &expected
        );
        assert_eq!(molecule.stereo_atom(StereoAtomId(0)).site_id(), AtomId(0));
    }

    #[rstest]
    fn test_stereo_atom_editor_view_mut_local_queries(virtual_ligand_molecule: Molecule) {
        let mut editor = virtual_ligand_molecule.edit();
        let view = editor.stereo_atom_mut(StereoAtomId(0));
        assert_eq!(view.kind(), StereoKind::Tetrahedral);
        assert_eq!(view.ligand_count(), 4);
        assert_eq!(view.atom_ligand_count(), 2);
        assert_eq!(view.implicit_hydrogen_count(), 1);
        assert_eq!(view.lone_pair_count(), 1);
        assert_eq!(
            view.ligand_position(AtomId(1)),
            Some(StereoLigandPosition(0))
        );
        assert_eq!(
            view.ligand_position(AtomId(4)),
            Some(StereoLigandPosition(3))
        );
        assert_eq!(view.ligand_position(AtomId(0)), None);
        assert_eq!(view.ligand_position(AtomId(99)), None);
        let (coset, atoms, hydrogens, lone_pairs) = (
            view.coset(),
            view.atom_ligand_ids(),
            view.implicit_hydrogen_atom_ids(),
            view.lone_pair_atom_ids(),
        );
        assert_eq!(coset, &StereoCoset::Lit(1));
        assert_eq!(atoms.collect::<Vec<_>>(), vec![AtomId(1), AtomId(4)]);
        assert_eq!(hydrogens.collect::<Vec<_>>(), vec![AtomId(0)]);
        assert_eq!(lone_pairs.collect::<Vec<_>>(), vec![AtomId(0)]);
    }

    #[rstest]
    #[case::kind(false)]
    #[case::coset(true)]
    #[should_panic(expected = "stereo view has a concrete")]
    fn test_stereo_atom_editor_view_mut_configuration_error(
        mut stereo_entries: MoleculeEntries,
        #[case] coset: bool,
    ) {
        stereo_entries.stereo_atoms[0].2 = StereoAtomForm::default();
        let virtual_ligand_molecule = Molecule::from_entries(stereo_entries);
        let mut editor = virtual_ligand_molecule.edit();
        let view = editor.stereo_atom_mut(StereoAtomId(0));
        if coset {
            let _ = view.coset();
        } else {
            let _ = view.kind();
        }
    }

    #[rstest]
    fn test_stereo_atom_editor_view_mut_attributes_mut(molecule: Molecule) {
        let mut editor = molecule.edit();
        let expected = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0));
        {
            let mut view = editor.stereo_atom_mut(StereoAtomId(0));
            assert_eq!(view.id(), StereoAtomId(0));
            view.attributes_mut().configuration = expected.configuration.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.kind(), StereoKind::Tetrahedral);
            assert_eq!(view.coset(), &StereoCoset::Lit(0));
            assert_eq!(view.constraints(), &expected.constraints);
            assert_eq!(view.site_id(), AtomId(0));
            assert_eq!(
                view.ligand_ids(),
                &[
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom)
                ]
            );
        }
        assert_eq!(editor.stereo_atom(StereoAtomId(0)).attributes(), &expected);
        assert_eq!(editor.stereo_atom(StereoAtomId(0)).site_id(), AtomId(0));
    }

    #[rstest]
    fn test_stereo_atom_editor_view_mut_replace_site(mut stereo_entries: MoleculeEntries) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        {
            let mut view = editor.stereo_atom_mut(StereoAtomId(0));
            view.replace_site(AtomId(5));
            assert_eq!(view.site_id(), AtomId(5));
            assert_eq!(view.ligand_ids(), stereo_entries.stereo_atoms[0].1);
            assert_eq!(view.attributes(), &stereo_entries.stereo_atoms[0].2);
        }
        stereo_entries.stereo_atoms[0].0 = AtomId(5);
        let molecule = editor.try_build().unwrap();
        assert_eq!(
            molecule.stereo_atoms().at_id(AtomId(5)),
            Some(StereoAtomId(0))
        );
        assert_eq!(molecule.stereo_atoms().at_id(AtomId(0)), None);
        assert_eq!(molecule, Molecule::from_entries(stereo_entries));
    }

    #[rstest]
    #[case::missing(AtomId(99), MoleculeIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(99)) })]
    #[case::duplicate(AtomId(6), MoleculeIntegrityError::DuplicateStereoAtomSites { atom: AtomId(6) })]
    fn test_stereo_atom_editor_view_mut_replace_site_publication(
        stereo_entries: MoleculeEntries,
        #[case] site: AtomId,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        editor.stereo_atom_mut(StereoAtomId(0)).replace_site(site);
        assert_eq!(editor.try_build(), Err(expected));
    }

    #[rstest]
    #[case::reordered(vec![StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)])]
    #[case::virtual_ligands(vec![StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::LonePair)])]
    #[case::empty(vec![])]
    #[case::single(vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom)])]
    fn test_stereo_atom_editor_view_mut_replace_ligands(
        mut stereo_entries: MoleculeEntries,
        #[case] ligands: Vec<StereoLigand>,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        {
            let mut view = editor.stereo_atom_mut(StereoAtomId(0));
            view.replace_ligands(&ligands);
            assert_eq!(view.ligand_ids(), ligands);
            assert_eq!(view.ligand_count(), ligands.len());
            assert_eq!(view.site_id(), stereo_entries.stereo_atoms[0].0);
            assert_eq!(view.attributes(), &stereo_entries.stereo_atoms[0].2);
        }
        stereo_entries.stereo_atoms[0].1 = ligands;
        assert_eq!(
            editor.try_build(),
            Ok(Molecule::from_entries(stereo_entries))
        );
    }

    #[rstest]
    #[case::missing_atom(vec![StereoLigand::new(AtomId(99), StereoLigandKind::Atom)], MoleculeIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(99)) })]
    #[case::duplicate(vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom); 2], MoleculeIntegrityError::DuplicateStereoLigand { entity: Entity::StereoAtom(StereoAtomId(0)), ligand: StereoLigand::new(AtomId(2), StereoLigandKind::Atom) })]
    #[case::wrong_anchor(vec![StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)], MoleculeIntegrityError::StereoLigandIncidenceMismatch { entity: Entity::StereoAtom(StereoAtomId(0)) })]
    fn test_stereo_atom_editor_view_mut_replace_ligands_publication(
        stereo_entries: MoleculeEntries,
        #[case] ligands: Vec<StereoLigand>,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        editor
            .stereo_atom_mut(StereoAtomId(0))
            .replace_ligands(&ligands);
        assert_eq!(editor.try_build(), Err(expected));
    }

    #[rstest]
    #[case::attributes_first(true)]
    #[case::ligands_first(false)]
    fn test_stereo_atom_editor_view_mut_replace_ligands_attributes(
        mut stereo_entries: MoleculeEntries,
        #[case] attributes_first: bool,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        let mut attributes = stereo_entries.stereo_atoms[0].2.clone();
        attributes.configuration =
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)).configuration;
        let ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
        ];
        {
            let mut view = editor.stereo_atom_mut(StereoAtomId(0));
            if attributes_first {
                *view.attributes_mut() = attributes.clone();
                view.replace_ligands(&ligands);
            } else {
                view.replace_ligands(&ligands);
                *view.attributes_mut() = attributes.clone();
            }
            assert_eq!(view.coset(), &StereoCoset::Lit(0));
        }
        stereo_entries.stereo_atoms[0].1 = ligands;
        stereo_entries.stereo_atoms[0].2 = attributes;
        assert_eq!(
            editor.try_build(),
            Ok(Molecule::from_entries(stereo_entries))
        );
    }

    #[rstest]
    fn test_stereo_atom_editor_view_mut_replace_ligands_incidence(
        mut stereo_entries: MoleculeEntries,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        {
            let mut view = editor.stereo_atom_mut(StereoAtomId(0));
            view.replace_ligand(
                StereoLigandPosition(0),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            );
            view.insert_ligand(
                StereoLigandPosition(1),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            );
            view.remove_ligand(StereoLigandPosition(2));
        }
        stereo_entries.stereo_atoms[0].1 = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
        ];
        let molecule = editor.try_build().unwrap();
        assert_eq!(
            molecule
                .stereo_atoms()
                .incident_as_ligand_ids(AtomId(3))
                .collect::<Vec<_>>(),
            vec![]
        );
        assert_eq!(
            molecule
                .stereo_atoms()
                .incident_as_ligand_ids(AtomId(2))
                .collect::<Vec<_>>(),
            vec![StereoAtomId(0)]
        );
        assert_eq!(molecule, Molecule::from_entries(stereo_entries));
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), vec![StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)])]
    #[case::last(StereoLigandPosition(3), vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen)])]
    fn test_stereo_atom_editor_view_mut_replace_ligand(
        stereo_entries: MoleculeEntries,
        #[case] position: StereoLigandPosition,
        #[case] expected: Vec<StereoLigand>,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        let mut view = editor.stereo_atom_mut(StereoAtomId(0));
        view.replace_ligand(
            position,
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        );
        assert_eq!(view.ligand_ids(), expected);
        assert_eq!(view.ligand_count(), expected.len());
        assert_eq!(view.site_id(), stereo_entries.stereo_atoms[0].0);
        assert_eq!(view.attributes(), &stereo_entries.stereo_atoms[0].2);
    }

    #[rstest]
    #[case::empty(0, StereoLigandPosition(0))]
    #[case::end(4, StereoLigandPosition(4))]
    #[should_panic]
    fn test_stereo_atom_editor_view_mut_replace_ligand_error(
        stereo_entries: MoleculeEntries,
        #[case] length: usize,
        #[case] position: StereoLigandPosition,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        let mut view = editor.stereo_atom_mut(StereoAtomId(0));
        let ligands = view.ligand_ids()[..length].to_vec();
        view.replace_ligands(&ligands);
        view.replace_ligand(
            position,
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        );
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), vec![StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)])]
    #[case::middle(StereoLigandPosition(2), vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)])]
    #[case::end(StereoLigandPosition(4), vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen)])]
    fn test_stereo_atom_editor_view_mut_insert_ligand(
        stereo_entries: MoleculeEntries,
        #[case] position: StereoLigandPosition,
        #[case] expected: Vec<StereoLigand>,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        let mut view = editor.stereo_atom_mut(StereoAtomId(0));
        view.insert_ligand(
            position,
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        );
        assert_eq!(view.ligand_ids(), expected);
        assert_eq!(view.ligand_count(), expected.len());
        assert_eq!(view.site_id(), stereo_entries.stereo_atoms[0].0);
        assert_eq!(view.attributes(), &stereo_entries.stereo_atoms[0].2);
    }

    #[rstest]
    fn test_stereo_atom_editor_view_mut_insert_ligand_empty(stereo_entries: MoleculeEntries) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        let mut view = editor.stereo_atom_mut(StereoAtomId(0));
        view.replace_ligands(&[]);
        view.insert_ligand(
            StereoLigandPosition(0),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        );
        assert_eq!(
            view.ligand_ids(),
            &[StereoLigand::new(AtomId(2), StereoLigandKind::Atom)]
        );
        view.remove_ligand(StereoLigandPosition(0));
        assert_eq!(view.ligand_ids(), &[]);
        assert_eq!(view.ligand_count(), 0);
    }

    #[rstest]
    #[case::empty(0, StereoLigandPosition(1))]
    #[case::end(4, StereoLigandPosition(5))]
    #[should_panic]
    fn test_stereo_atom_editor_view_mut_insert_ligand_error(
        stereo_entries: MoleculeEntries,
        #[case] length: usize,
        #[case] position: StereoLigandPosition,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        let mut view = editor.stereo_atom_mut(StereoAtomId(0));
        let ligands = view.ligand_ids()[..length].to_vec();
        view.replace_ligands(&ligands);
        view.insert_ligand(
            position,
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        );
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), vec![StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)])]
    #[case::last(StereoLigandPosition(3), vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)])]
    fn test_stereo_atom_editor_view_mut_remove_ligand(
        stereo_entries: MoleculeEntries,
        #[case] position: StereoLigandPosition,
        #[case] expected: Vec<StereoLigand>,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        let mut view = editor.stereo_atom_mut(StereoAtomId(0));
        view.remove_ligand(position);
        assert_eq!(view.ligand_ids(), expected);
        assert_eq!(view.ligand_count(), expected.len());
        assert_eq!(view.site_id(), stereo_entries.stereo_atoms[0].0);
        assert_eq!(view.attributes(), &stereo_entries.stereo_atoms[0].2);
    }

    #[rstest]
    #[case::empty(0, StereoLigandPosition(0))]
    #[case::end(4, StereoLigandPosition(4))]
    #[should_panic]
    fn test_stereo_atom_editor_view_mut_remove_ligand_error(
        stereo_entries: MoleculeEntries,
        #[case] length: usize,
        #[case] position: StereoLigandPosition,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        let mut view = editor.stereo_atom_mut(StereoAtomId(0));
        let ligands = view.ligand_ids()[..length].to_vec();
        view.replace_ligands(&ligands);
        view.remove_ligand(position);
    }

    #[rstest]
    fn test_stereo_bond_views_count(molecule: Molecule) {
        assert_eq!(molecule.stereo_bonds().count(), 1);
    }

    #[rstest]
    fn test_stereo_bond_views_ids(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().stereo_bonds().ids(), vec![], |id| id);
        assert_exact_size_by(molecule.stereo_bonds().ids(), vec![StereoBondId(0)], |id| {
            id
        });
    }

    #[rstest]
    fn test_stereo_bond_views_iter(molecule: Molecule) {
        assert_exact_size_by(Molecule::default().stereo_bonds().iter(), vec![], |view| {
            (view.id(), view.site_id())
        });
        assert_exact_size_by(
            molecule.stereo_bonds().iter(),
            vec![(StereoBondId(0), BondId(1))],
            |view| (view.id(), view.site_id()),
        );
    }

    #[rstest]
    #[case::present(StereoBondId(0), true)]
    #[case::absent(StereoBondId(99), false)]
    fn test_stereo_bond_views_contains(
        molecule: Molecule,
        #[case] id: StereoBondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.stereo_bonds().contains(id), expected);
    }

    #[rstest]
    fn test_stereo_bond_views_get(molecule: Molecule) {
        let res = molecule.stereo_bonds().get(StereoBondId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id(), StereoBondId(0));
        assert_eq!(view.site_id(), BondId(1));
        assert_eq!(view.kind(), StereoKind::CisTrans);
        assert_eq!(
            view.attributes(),
            &StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
        );
    }

    #[rstest]
    fn test_stereo_bond_views_get_none(molecule: Molecule) {
        let res = molecule.stereo_bonds().get(StereoBondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    #[case::site_endpoint(AtomId(2), vec![StereoBondId(0)])]
    #[case::ligand(AtomId(4), vec![StereoBondId(0)])]
    #[case::unrelated(AtomId(6), vec![])]
    fn test_stereo_bond_views_incident_to_atom_ids(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoBondId>,
    ) {
        assert_eq!(
            molecule
                .stereo_bonds()
                .incident_to_atom_ids(atom)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::site_endpoint(AtomId(2), vec![StereoBondId(0)])]
    #[case::ligand_not_site(AtomId(4), vec![])]
    #[case::unrelated(AtomId(6), vec![])]
    fn test_stereo_bond_views_incident_as_site(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoBondId>,
    ) {
        assert_eq!(
            molecule
                .stereo_bonds()
                .incident_as_site_ids(atom)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::ligand(AtomId(4), vec![StereoBondId(0)])]
    #[case::site_not_ligand(AtomId(2), vec![])]
    #[case::unrelated(AtomId(6), vec![])]
    fn test_stereo_bond_views_incident_as_ligand(
        molecule: Molecule,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoBondId>,
    ) {
        assert_eq!(
            molecule
                .stereo_bonds()
                .incident_as_ligand_ids(atom)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::site(BondId(1), Some(StereoBondId(0)))]
    #[case::non_site(BondId(0), None)]
    fn test_stereo_bond_views_at(
        molecule: Molecule,
        #[case] bond: BondId,
        #[case] expected: Option<StereoBondId>,
    ) {
        assert_eq!(molecule.stereo_bonds().at_id(bond), expected);
    }

    #[rstest]
    fn test_stereo_bond_view_site_id(molecule: Molecule) {
        assert_eq!(molecule.stereo_bond(StereoBondId(0)).site_id(), BondId(1));
    }

    #[rstest]
    fn test_stereo_bond_view_site(molecule: Molecule) {
        let view = molecule.stereo_bond(StereoBondId(0)).site();
        assert_eq!(view.id(), BondId(1));
        assert_eq!(view.atom_ids(), [AtomId(2), AtomId(3)]);
    }

    #[rstest]
    fn test_stereo_bond_view_coset(molecule: Molecule) {
        assert_eq!(
            molecule.stereo_bond(StereoBondId(0)).coset(),
            &StereoCoset::Lit(1),
        );
    }

    #[rstest]
    fn test_stereo_bond_view_ligand_count(molecule: Molecule) {
        assert_eq!(molecule.stereo_bond(StereoBondId(0)).ligand_count(), 4);
    }

    #[rstest]
    fn test_stereo_bond_view_ligands(molecule: Molecule) {
        assert_exact_size_by(
            molecule.stereo_bond(StereoBondId(0)).ligands(),
            vec![
                (StereoLigandKind::Atom, AtomId(4)),
                (StereoLigandKind::Atom, AtomId(5)),
                (StereoLigandKind::Atom, AtomId(0)),
                (StereoLigandKind::Atom, AtomId(1)),
            ],
            |ligand| (ligand.kind(), ligand.atom_id()),
        );
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), StereoLigandKind::Atom, AtomId(4))]
    #[case::second(StereoLigandPosition(1), StereoLigandKind::Atom, AtomId(5))]
    #[case::third(StereoLigandPosition(2), StereoLigandKind::Atom, AtomId(0))]
    #[case::fourth(StereoLigandPosition(3), StereoLigandKind::Atom, AtomId(1))]
    fn test_stereo_bond_view_ligand(
        molecule: Molecule,
        #[case] ligand_id: StereoLigandPosition,
        #[case] kind: StereoLigandKind,
        #[case] atom: AtomId,
    ) {
        let ligand = molecule.stereo_bond(StereoBondId(0)).ligand(ligand_id);
        assert_eq!(ligand.kind(), kind);
        assert_eq!(ligand.atom_id(), atom);
    }

    #[rstest]
    #[should_panic]
    fn test_stereo_bond_view_ligand_error(molecule: Molecule) {
        molecule
            .stereo_bond(StereoBondId(0))
            .ligand(StereoLigandPosition(4));
    }

    #[rstest]
    fn test_stereo_bond_view_atom_ligands(virtual_ligand_molecule: Molecule) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .atom_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(4), AtomId(5)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .atom_ligand_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(4), AtomId(5)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .atom_ligand_count(),
            2,
        );
    }

    #[rstest]
    fn test_stereo_bond_view_implicit_hydrogen_ligands(virtual_ligand_molecule: Molecule) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .implicit_hydrogen_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(2)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .implicit_hydrogen_atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .implicit_hydrogen_count(),
            1,
        );
    }

    #[rstest]
    fn test_stereo_bond_view_lone_pair_ligands(virtual_ligand_molecule: Molecule) {
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .lone_pair_ligands()
                .map(|ligand| ligand.atom_id())
                .collect::<Vec<_>>(),
            vec![AtomId(3)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .lone_pair_atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(3)],
        );
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .lone_pair_count(),
            1,
        );
    }

    #[rstest]
    fn test_stereo_bond_view_coset_for_none(molecule: Molecule) {
        let view = molecule.stereo_bond(StereoBondId(0));
        let ligands = [
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(5),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(0),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(1),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(view.coset_for(ligands[..1].iter().copied()), None);
        assert_eq!(
            view.coset_for([ligands[0], ligands[0], ligands[2], ligands[3]]),
            None,
        );
        assert_eq!(
            view.coset_for([
                ligands[0],
                ligands[1],
                ligands[2],
                StereoLigand {
                    atom_id: AtomId(99),
                    kind: StereoLigandKind::Atom,
                },
            ]),
            None,
        );
    }

    #[rstest]
    fn test_stereo_bond_view_coset_for(molecule: Molecule) {
        let view = molecule.stereo_bond(StereoBondId(0));
        let ligands = vec![
            StereoLigand {
                atom_id: AtomId(4),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(5),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(0),
                kind: StereoLigandKind::Atom,
            },
            StereoLigand {
                atom_id: AtomId(1),
                kind: StereoLigandKind::Atom,
            },
        ];
        assert_eq!(view.coset_for(ligands.clone()), Some(StereoCoset::Lit(1)),);

        let reordered = vec![ligands[1], ligands[0], ligands[2], ligands[3]];
        assert_eq!(view.coset_for(reordered), Some(StereoCoset::Lit(0)));
    }

    #[rstest]
    fn test_stereo_bond_view_atom_ids(molecule: Molecule, virtual_ligand_molecule: Molecule) {
        assert_eq!(
            molecule
                .stereo_bond(StereoBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![
                AtomId(2),
                AtomId(3),
                AtomId(4),
                AtomId(5),
                AtomId(0),
                AtomId(1)
            ],
        );
        // Virtual ligands sit on the site-bond endpoints, so they are not repeated.
        assert_eq!(
            virtual_ligand_molecule
                .stereo_bond(StereoBondId(0))
                .atom_ids()
                .collect::<Vec<_>>(),
            vec![AtomId(2), AtomId(3), AtomId(4), AtomId(5)],
        );
    }

    #[rstest]
    fn test_stereo_bond_view_attributes(molecule: Molecule) {
        let (attributes, frame) = {
            let view = molecule.stereo_bond(StereoBondId(0));
            assert_eq!(view.id(), StereoBondId(0));
            (view.attributes(), view.ligand_ids())
        };
        assert_eq!(
            attributes,
            &StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))
        );
        assert_eq!(
            frame,
            &[
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom)
            ]
        );
    }

    #[rstest]
    fn test_stereo_bond_editor_view_attributes(molecule: Molecule) {
        let editor = molecule.edit();
        let (attributes, frame) = {
            let view = editor.stereo_bond(StereoBondId(0));
            assert_eq!(view.id(), StereoBondId(0));
            (view.attributes(), view.ligand_ids())
        };
        assert_eq!(
            attributes,
            &StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))
        );
        assert_eq!(
            frame,
            &[
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom)
            ]
        );
    }

    #[rstest]
    fn test_stereo_bond_editor_view_local_queries(virtual_ligand_molecule: Molecule) {
        let editor = virtual_ligand_molecule.edit();
        let (coset, atoms, hydrogens, lone_pairs) = {
            let view = editor.stereo_bond(StereoBondId(0));
            assert_eq!(view.kind(), StereoKind::CisTrans);
            assert_eq!(view.ligand_count(), 4);
            assert_eq!(view.atom_ligand_count(), 2);
            assert_eq!(view.implicit_hydrogen_count(), 1);
            assert_eq!(view.lone_pair_count(), 1);
            assert_eq!(
                view.ligand_position(AtomId(4)),
                Some(StereoLigandPosition(0))
            );
            assert_eq!(
                view.ligand_position(AtomId(5)),
                Some(StereoLigandPosition(2))
            );
            assert_eq!(view.ligand_position(AtomId(2)), None);
            assert_eq!(view.ligand_position(AtomId(3)), None);
            assert_eq!(view.ligand_position(AtomId(99)), None);
            (
                view.coset(),
                view.atom_ligand_ids(),
                view.implicit_hydrogen_atom_ids(),
                view.lone_pair_atom_ids(),
            )
        };
        assert_eq!(coset, &StereoCoset::Lit(1));
        assert_eq!(atoms.collect::<Vec<_>>(), vec![AtomId(4), AtomId(5)]);
        assert_eq!(hydrogens.collect::<Vec<_>>(), vec![AtomId(2)]);
        assert_eq!(lone_pairs.collect::<Vec<_>>(), vec![AtomId(3)]);
    }

    #[rstest]
    #[case::kind(false)]
    #[case::coset(true)]
    #[should_panic(expected = "stereo view has a concrete")]
    fn test_stereo_bond_editor_view_configuration_error(
        mut stereo_entries: MoleculeEntries,
        #[case] coset: bool,
    ) {
        stereo_entries.stereo_bonds[0].2 = StereoBondForm::default();
        let virtual_ligand_molecule = Molecule::from_entries(stereo_entries);
        let editor = virtual_ligand_molecule.edit();
        let view = editor.stereo_bond(StereoBondId(0));
        if coset {
            let _ = view.coset();
        } else {
            let _ = view.kind();
        }
    }

    #[rstest]
    #[should_panic(expected = "invalid stereo bond id")]
    fn test_stereo_bond_editor_view_id_error(molecule: Molecule) {
        molecule.edit().stereo_bond(StereoBondId(1));
    }

    #[rstest]
    fn test_stereo_bond_view_mut_local_queries(virtual_ligand_molecule: Molecule) {
        let mut molecule = virtual_ligand_molecule;
        let view = molecule.stereo_bond_mut(StereoBondId(0));
        assert_eq!(view.kind(), StereoKind::CisTrans);
        assert_eq!(view.ligand_count(), 4);
        assert_eq!(view.atom_ligand_count(), 2);
        assert_eq!(view.implicit_hydrogen_count(), 1);
        assert_eq!(view.lone_pair_count(), 1);
        assert_eq!(
            view.ligand_position(AtomId(4)),
            Some(StereoLigandPosition(0))
        );
        assert_eq!(
            view.ligand_position(AtomId(5)),
            Some(StereoLigandPosition(2))
        );
        assert_eq!(view.ligand_position(AtomId(2)), None);
        assert_eq!(view.ligand_position(AtomId(3)), None);
        assert_eq!(view.ligand_position(AtomId(99)), None);
        let (coset, atoms, hydrogens, lone_pairs) = (
            view.coset(),
            view.atom_ligand_ids(),
            view.implicit_hydrogen_atom_ids(),
            view.lone_pair_atom_ids(),
        );
        assert_eq!(coset, &StereoCoset::Lit(1));
        assert_eq!(atoms.collect::<Vec<_>>(), vec![AtomId(4), AtomId(5)]);
        assert_eq!(hydrogens.collect::<Vec<_>>(), vec![AtomId(2)]);
        assert_eq!(lone_pairs.collect::<Vec<_>>(), vec![AtomId(3)]);
    }

    #[rstest]
    #[case::kind(false)]
    #[case::coset(true)]
    #[should_panic(expected = "stereo view has a concrete")]
    fn test_stereo_bond_view_mut_configuration_error(
        mut stereo_entries: MoleculeEntries,
        #[case] coset: bool,
    ) {
        stereo_entries.stereo_bonds[0].2 = StereoBondForm::default();
        let virtual_ligand_molecule = Molecule::from_entries(stereo_entries);
        let mut molecule = virtual_ligand_molecule;
        let view = molecule.stereo_bond_mut(StereoBondId(0));
        if coset {
            let _ = view.coset();
        } else {
            let _ = view.kind();
        }
    }

    #[rstest]
    fn test_stereo_bond_view_mut_attributes_mut(mut molecule: Molecule) {
        let expected = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0));
        {
            let mut view = molecule.stereo_bond_mut(StereoBondId(0));
            assert_eq!(view.id(), StereoBondId(0));
            view.attributes_mut().configuration = expected.configuration.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.kind(), StereoKind::CisTrans);
            assert_eq!(view.coset(), &StereoCoset::Lit(0));
            assert_eq!(view.constraints(), &expected.constraints);
            assert_eq!(view.site_id(), BondId(1));
            assert_eq!(
                view.ligand_ids(),
                &[
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom)
                ]
            );
        }
        assert_eq!(
            molecule.stereo_bond(StereoBondId(0)).attributes(),
            &expected
        );
        assert_eq!(molecule.stereo_bond(StereoBondId(0)).site_id(), BondId(1));
    }

    #[rstest]
    fn test_stereo_bond_editor_view_mut_local_queries(virtual_ligand_molecule: Molecule) {
        let mut editor = virtual_ligand_molecule.edit();
        let view = editor.stereo_bond_mut(StereoBondId(0));
        assert_eq!(view.kind(), StereoKind::CisTrans);
        assert_eq!(view.ligand_count(), 4);
        assert_eq!(view.atom_ligand_count(), 2);
        assert_eq!(view.implicit_hydrogen_count(), 1);
        assert_eq!(view.lone_pair_count(), 1);
        assert_eq!(
            view.ligand_position(AtomId(4)),
            Some(StereoLigandPosition(0))
        );
        assert_eq!(
            view.ligand_position(AtomId(5)),
            Some(StereoLigandPosition(2))
        );
        assert_eq!(view.ligand_position(AtomId(2)), None);
        assert_eq!(view.ligand_position(AtomId(3)), None);
        assert_eq!(view.ligand_position(AtomId(99)), None);
        let (coset, atoms, hydrogens, lone_pairs) = (
            view.coset(),
            view.atom_ligand_ids(),
            view.implicit_hydrogen_atom_ids(),
            view.lone_pair_atom_ids(),
        );
        assert_eq!(coset, &StereoCoset::Lit(1));
        assert_eq!(atoms.collect::<Vec<_>>(), vec![AtomId(4), AtomId(5)]);
        assert_eq!(hydrogens.collect::<Vec<_>>(), vec![AtomId(2)]);
        assert_eq!(lone_pairs.collect::<Vec<_>>(), vec![AtomId(3)]);
    }

    #[rstest]
    #[case::kind(false)]
    #[case::coset(true)]
    #[should_panic(expected = "stereo view has a concrete")]
    fn test_stereo_bond_editor_view_mut_configuration_error(
        mut stereo_entries: MoleculeEntries,
        #[case] coset: bool,
    ) {
        stereo_entries.stereo_bonds[0].2 = StereoBondForm::default();
        let virtual_ligand_molecule = Molecule::from_entries(stereo_entries);
        let mut editor = virtual_ligand_molecule.edit();
        let view = editor.stereo_bond_mut(StereoBondId(0));
        if coset {
            let _ = view.coset();
        } else {
            let _ = view.kind();
        }
    }

    #[rstest]
    fn test_stereo_bond_editor_view_mut_attributes_mut(molecule: Molecule) {
        let mut editor = molecule.edit();
        let expected = StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0));
        {
            let mut view = editor.stereo_bond_mut(StereoBondId(0));
            assert_eq!(view.id(), StereoBondId(0));
            view.attributes_mut().configuration = expected.configuration.clone();
            assert_eq!(view.attributes(), &expected);
            assert_eq!(view.kind(), StereoKind::CisTrans);
            assert_eq!(view.coset(), &StereoCoset::Lit(0));
            assert_eq!(view.constraints(), &expected.constraints);
            assert_eq!(view.site_id(), BondId(1));
            assert_eq!(
                view.ligand_ids(),
                &[
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom)
                ]
            );
        }
        assert_eq!(editor.stereo_bond(StereoBondId(0)).attributes(), &expected);
        assert_eq!(editor.stereo_bond(StereoBondId(0)).site_id(), BondId(1));
    }

    #[rstest]
    fn test_stereo_bond_editor_view_mut_replace_site(mut stereo_entries: MoleculeEntries) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        {
            let mut view = editor.stereo_bond_mut(StereoBondId(0));
            view.replace_site(BondId(1));
            assert_eq!(view.site_id(), BondId(1));
            assert_eq!(view.ligand_ids(), stereo_entries.stereo_bonds[0].1);
            assert_eq!(view.attributes(), &stereo_entries.stereo_bonds[0].2);
        }
        stereo_entries.stereo_bonds[0].0 = BondId(1);
        let molecule = editor.try_build().unwrap();
        assert_eq!(
            molecule.stereo_bonds().at_id(BondId(1)),
            Some(StereoBondId(0))
        );
        assert_eq!(molecule.stereo_bonds().at_id(BondId(0)), None);
        assert_eq!(molecule, Molecule::from_entries(stereo_entries));
    }

    #[rstest]
    #[case::missing(BondId(99), MoleculeIntegrityError::InvalidReference { entity: Entity::Bond(BondId(99)) })]
    #[case::duplicate(BondId(2), MoleculeIntegrityError::DuplicateStereoBondSites { bond: BondId(2) })]
    fn test_stereo_bond_editor_view_mut_replace_site_publication(
        stereo_entries: MoleculeEntries,
        #[case] site: BondId,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        editor.stereo_bond_mut(StereoBondId(0)).replace_site(site);
        assert_eq!(editor.try_build(), Err(expected));
    }

    #[rstest]
    #[case::reordered(vec![StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)])]
    #[case::virtual_ligands(vec![StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(1), StereoLigandKind::LonePair)])]
    fn test_stereo_bond_editor_view_mut_replace_ligands(
        mut stereo_entries: MoleculeEntries,
        #[case] ligands: Vec<StereoLigand>,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        {
            let mut view = editor.stereo_bond_mut(StereoBondId(0));
            view.replace_ligands(&ligands);
            assert_eq!(view.ligand_ids(), ligands);
            assert_eq!(view.ligand_count(), ligands.len());
            assert_eq!(view.site_id(), stereo_entries.stereo_bonds[0].0);
            assert_eq!(view.attributes(), &stereo_entries.stereo_bonds[0].2);
        }
        stereo_entries.stereo_bonds[0].1 = ligands;
        assert_eq!(
            editor.try_build(),
            Ok(Molecule::from_entries(stereo_entries))
        );
    }

    #[rstest]
    #[case::missing_atom(vec![StereoLigand::new(AtomId(99), StereoLigandKind::Atom)], MoleculeIntegrityError::InvalidReference { entity: Entity::Atom(AtomId(99)) })]
    #[case::duplicate(vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom); 2], MoleculeIntegrityError::DuplicateStereoLigand { entity: Entity::StereoBond(StereoBondId(0)), ligand: StereoLigand::new(AtomId(2), StereoLigandKind::Atom) })]
    #[case::wrong_anchor(vec![StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)], MoleculeIntegrityError::StereoLigandIncidenceMismatch { entity: Entity::StereoBond(StereoBondId(0)) })]
    fn test_stereo_bond_editor_view_mut_replace_ligands_publication(
        stereo_entries: MoleculeEntries,
        #[case] ligands: Vec<StereoLigand>,
        #[case] expected: MoleculeIntegrityError,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        editor
            .stereo_bond_mut(StereoBondId(0))
            .replace_ligands(&ligands);
        assert_eq!(editor.try_build(), Err(expected));
    }

    #[rstest]
    #[case::attributes_first(true)]
    #[case::ligands_first(false)]
    fn test_stereo_bond_editor_view_mut_replace_ligands_attributes(
        mut stereo_entries: MoleculeEntries,
        #[case] attributes_first: bool,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        let mut attributes = stereo_entries.stereo_bonds[0].2.clone();
        attributes.configuration =
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)).configuration;
        let ligands = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
        ];
        {
            let mut view = editor.stereo_bond_mut(StereoBondId(0));
            if attributes_first {
                *view.attributes_mut() = attributes.clone();
                view.replace_ligands(&ligands);
            } else {
                view.replace_ligands(&ligands);
                *view.attributes_mut() = attributes.clone();
            }
            assert_eq!(view.coset(), &StereoCoset::Lit(0));
        }
        stereo_entries.stereo_bonds[0].1 = ligands;
        stereo_entries.stereo_bonds[0].2 = attributes;
        assert_eq!(
            editor.try_build(),
            Ok(Molecule::from_entries(stereo_entries))
        );
    }

    #[rstest]
    fn test_stereo_bond_editor_view_mut_replace_ligands_incidence(
        mut stereo_entries: MoleculeEntries,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        {
            let mut view = editor.stereo_bond_mut(StereoBondId(0));
            view.replace_ligand(
                StereoLigandPosition(0),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            );
            view.insert_ligand(
                StereoLigandPosition(1),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            );
            view.remove_ligand(StereoLigandPosition(2));
        }
        stereo_entries.stereo_bonds[0].1 = vec![
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(7), StereoLigandKind::Atom),
        ];
        let molecule = editor.try_build().unwrap();
        assert_eq!(
            molecule
                .stereo_bonds()
                .incident_as_ligand_ids(AtomId(3))
                .collect::<Vec<_>>(),
            vec![StereoBondId(1)]
        );
        assert_eq!(
            molecule
                .stereo_bonds()
                .incident_as_ligand_ids(AtomId(2))
                .collect::<Vec<_>>(),
            vec![StereoBondId(0), StereoBondId(1)]
        );
        assert_eq!(molecule, Molecule::from_entries(stereo_entries));
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), vec![StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)])]
    #[case::last(StereoLigandPosition(3), vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen)])]
    fn test_stereo_bond_editor_view_mut_replace_ligand(
        stereo_entries: MoleculeEntries,
        #[case] position: StereoLigandPosition,
        #[case] expected: Vec<StereoLigand>,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        let mut view = editor.stereo_bond_mut(StereoBondId(0));
        view.replace_ligand(
            position,
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        );
        assert_eq!(view.ligand_ids(), expected);
        assert_eq!(view.ligand_count(), expected.len());
        assert_eq!(view.site_id(), stereo_entries.stereo_bonds[0].0);
        assert_eq!(view.attributes(), &stereo_entries.stereo_bonds[0].2);
    }

    #[rstest]
    #[case::empty(0, StereoLigandPosition(0))]
    #[case::end(4, StereoLigandPosition(4))]
    #[should_panic]
    fn test_stereo_bond_editor_view_mut_replace_ligand_error(
        stereo_entries: MoleculeEntries,
        #[case] length: usize,
        #[case] position: StereoLigandPosition,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        let mut view = editor.stereo_bond_mut(StereoBondId(0));
        let ligands = view.ligand_ids()[..length].to_vec();
        view.replace_ligands(&ligands);
        view.replace_ligand(
            position,
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        );
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), vec![StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)])]
    #[case::middle(StereoLigandPosition(2), vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)])]
    #[case::end(StereoLigandPosition(4), vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom), StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen)])]
    fn test_stereo_bond_editor_view_mut_insert_ligand(
        stereo_entries: MoleculeEntries,
        #[case] position: StereoLigandPosition,
        #[case] expected: Vec<StereoLigand>,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        let mut view = editor.stereo_bond_mut(StereoBondId(0));
        view.insert_ligand(
            position,
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        );
        assert_eq!(view.ligand_ids(), expected);
        assert_eq!(view.ligand_count(), expected.len());
        assert_eq!(view.site_id(), stereo_entries.stereo_bonds[0].0);
        assert_eq!(view.attributes(), &stereo_entries.stereo_bonds[0].2);
    }

    #[rstest]
    fn test_stereo_bond_editor_view_mut_insert_ligand_empty(stereo_entries: MoleculeEntries) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        let mut view = editor.stereo_bond_mut(StereoBondId(0));
        view.replace_ligands(&[]);
        view.insert_ligand(
            StereoLigandPosition(0),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        );
        assert_eq!(
            view.ligand_ids(),
            &[StereoLigand::new(AtomId(2), StereoLigandKind::Atom)]
        );
        view.remove_ligand(StereoLigandPosition(0));
        assert_eq!(view.ligand_ids(), &[]);
        assert_eq!(view.ligand_count(), 0);
    }

    #[rstest]
    #[case::empty(0, StereoLigandPosition(1))]
    #[case::end(4, StereoLigandPosition(5))]
    #[should_panic]
    fn test_stereo_bond_editor_view_mut_insert_ligand_error(
        stereo_entries: MoleculeEntries,
        #[case] length: usize,
        #[case] position: StereoLigandPosition,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        let mut view = editor.stereo_bond_mut(StereoBondId(0));
        let ligands = view.ligand_ids()[..length].to_vec();
        view.replace_ligands(&ligands);
        view.insert_ligand(
            position,
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
        );
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), vec![StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom), StereoLigand::new(AtomId(7), StereoLigandKind::Atom)])]
    #[case::last(StereoLigandPosition(3), vec![StereoLigand::new(AtomId(2), StereoLigandKind::Atom), StereoLigand::new(AtomId(3), StereoLigandKind::Atom), StereoLigand::new(AtomId(4), StereoLigandKind::Atom)])]
    fn test_stereo_bond_editor_view_mut_remove_ligand(
        stereo_entries: MoleculeEntries,
        #[case] position: StereoLigandPosition,
        #[case] expected: Vec<StereoLigand>,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries.clone()).edit();
        let mut view = editor.stereo_bond_mut(StereoBondId(0));
        view.remove_ligand(position);
        assert_eq!(view.ligand_ids(), expected);
        assert_eq!(view.ligand_count(), expected.len());
        assert_eq!(view.site_id(), stereo_entries.stereo_bonds[0].0);
        assert_eq!(view.attributes(), &stereo_entries.stereo_bonds[0].2);
    }

    #[rstest]
    #[case::empty(0, StereoLigandPosition(0))]
    #[case::end(4, StereoLigandPosition(4))]
    #[should_panic]
    fn test_stereo_bond_editor_view_mut_remove_ligand_error(
        stereo_entries: MoleculeEntries,
        #[case] length: usize,
        #[case] position: StereoLigandPosition,
    ) {
        let mut editor = Molecule::from_entries(stereo_entries).edit();
        let mut view = editor.stereo_bond_mut(StereoBondId(0));
        let ligands = view.ligand_ids()[..length].to_vec();
        view.replace_ligands(&ligands);
        view.remove_ligand(position);
    }

    #[rstest]
    #[case::exact(
        AtomId(0),
        vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ],
        Some(StereoAtomId(0))
    )]
    #[case::reordered(
        AtomId(0),
        vec![
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        ],
        Some(StereoAtomId(0))
    )]
    #[case::wrong_site(
        AtomId(6),
        vec![
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
        ],
        None
    )]
    fn test_stereo_atom_views_of_id(
        molecule: Molecule,
        #[case] site: AtomId,
        #[case] ligands: Vec<StereoLigand>,
        #[case] expected: Option<StereoAtomId>,
    ) {
        assert_eq!(molecule.stereo_atoms().of_id(site, &ligands), expected);
    }

    #[rstest]
    #[case::exact(
        BondId(1),
        vec![
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        ],
        Some(StereoBondId(0))
    )]
    #[case::wrong_site(
        BondId(0),
        vec![
            StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        ],
        None
    )]
    fn test_stereo_bond_views_of_id(
        molecule: Molecule,
        #[case] site: BondId,
        #[case] ligands: Vec<StereoLigand>,
        #[case] expected: Option<StereoBondId>,
    ) {
        assert_eq!(molecule.stereo_bonds().of_id(site, &ligands), expected);
    }
}
