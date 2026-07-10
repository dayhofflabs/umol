//! Stereo atom and stereo bond views.

use std::collections::HashSet;
use std::iter;
use std::ops::Index;

use umol_graph_core::{EdgeId, FixedVarBirelationSet, NodeId, Ordered, RelationId};
use umol_perm::{OrientedPermutationGroup, Permutation};

use super::super::id::{AtomId, BondId, StereoAtomId, StereoBondId, StereoLigandPosition};
use super::super::ligand::{StereoLigand, StereoLigandKind};
use super::super::molecule::MoleculeAst;
use super::super::ring::RingView;
use super::super::stereo::{
    coset_apply_permutation, StereoAtomAst, StereoBondAst, StereoKind, Stereogenicity, Topicity,
};
use super::super::symmetry::StereoSymmetry;
use super::super::traits::Lattice;
use super::atom::AtomView;
use super::bond::BondView;
use super::ligand::StereoLigandView;
use crate::ast::{StereoAtomConstraintsAst, StereoBondConstraintsAst, StereoCosetAst};

type StereoAtomSet =
    FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, StereoAtomAst>;
type StereoBondSet =
    FixedVarBirelationSet<EdgeId, Ordered, 1, StereoLigand, Ordered, StereoBondAst>;

/// Namespace accessor for stereo-atom views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct StereoAtomViews<'a> {
    molecule: &'a MoleculeAst,
    stereo_atoms: &'a StereoAtomSet,
}

impl<'a> StereoAtomViews<'a> {
    pub(crate) fn new(molecule: &'a MoleculeAst, stereo_atoms: &'a StereoAtomSet) -> Self {
        Self {
            molecule,
            stereo_atoms,
        }
    }

    pub fn count(&self) -> usize {
        self.stereo_atoms.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = StereoAtomId> {
        self.stereo_atoms.relation_ids().map(StereoAtomId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = StereoAtomView<'a>> {
        let molecule = self.molecule;
        let set = self.stereo_atoms;
        set.relation_ids().map(move |rid| StereoAtomView {
            id: StereoAtomId::from(rid),
            site: set.participants_1(rid)[0],
            ligands: set.participants_2(rid),
            ast: set.data(rid),
            molecule,
        })
    }

    /// Whether two stereo atoms share a site — the "≤1 stereo center per site" structural conflict.
    /// The per-entity emit-compliance primitive: the entity-structure validator and the generating ops
    /// (`apply_at` / `meet_pushout`) both consult it, and no generating op may emit a molecule for
    /// which it holds.
    pub fn has_conflict(&self) -> bool {
        let mut sites: HashSet<AtomId> = HashSet::new();
        self.iter().any(|view| !sites.insert(view.site_id()))
    }

    pub fn contains(&self, id: StereoAtomId) -> bool {
        self.stereo_atoms.contains(RelationId::from(id))
    }

    pub fn get(&self, id: StereoAtomId) -> Option<StereoAtomView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let rid = RelationId::from(id);
        Some(StereoAtomView {
            id,
            site: self.stereo_atoms.participants_1(rid)[0],
            ligands: self.stereo_atoms.participants_2(rid),
            ast: self.stereo_atoms.data(rid),
            molecule: self.molecule,
        })
    }

    /// Ids of stereo atoms incident on `atom` (site or ligand).
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = StereoAtomId> + 'a {
        self.stereo_atoms
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| StereoAtomId::from(rid))
    }

    /// Id of the stereo atom on `site` with exactly this ligand multiset, if any. Ligands are a
    /// multiset (virtual ligands repeat); their frame order is not matched.
    pub fn connecting_id(&self, site: AtomId, ligands: &[StereoLigand]) -> Option<StereoAtomId> {
        self.stereo_atoms
            .find_by_participants(&[NodeId::from(site)], ligands)
            .map(StereoAtomId::from)
    }

    /// Any stereo atom is incident on `atom` (site or ligand).
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.stereo_atoms.has_incident(NodeId::from(atom))
    }

    /// Views of stereo atoms incident on `atom` (site or ligand).
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = StereoAtomView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.stereo_atoms;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            StereoAtomView {
                id,
                site: set.participants_1(rid)[0],
                ligands: set.participants_2(rid),
                ast: set.data(rid),
                molecule,
            }
        })
    }

    // Ids of stereo atoms incident, in which `atom` is ligand.
    pub fn incident_at_ligand_ids(&self, atom: AtomId) -> impl Iterator<Item = StereoAtomId> + 'a {
        let set = self.stereo_atoms;
        let ligand = StereoLigand {
            atom_id: atom,
            kind: StereoLigandKind::Atom,
        };
        set.incident(NodeId::from(atom))
            .iter()
            .filter(move |&&rid| set.participants_2(rid).contains(&ligand))
            .map(|&rid| StereoAtomId::from(rid))
    }

    /// Any stereo atom is incident, in which `atom` is ligand.
    pub fn has_incident_at_ligand(&self, atom: AtomId) -> bool {
        let set = self.stereo_atoms;
        let ligand = StereoLigand {
            atom_id: atom,
            kind: StereoLigandKind::Atom,
        };
        set.incident(NodeId::from(atom))
            .iter()
            .any(|&rid| set.participants_2(rid).contains(&ligand))
    }

    /// Views of stereo atoms incident, in which `atom` is ligand.
    pub fn incident_at_ligand(
        &self,
        atom: AtomId,
    ) -> impl Iterator<Item = StereoAtomView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.stereo_atoms;
        self.incident_at_ligand_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            StereoAtomView {
                id,
                site: set.participants_1(rid)[0],
                ligands: set.participants_2(rid),
                ast: set.data(rid),
                molecule,
            }
        })
    }

    /// Id of the stereo atom coincident with `atom` (its site).
    pub fn coincident_id(&self, atom: AtomId) -> Option<StereoAtomId> {
        let set = self.stereo_atoms;
        let site = NodeId::from(atom);
        set.incident(site)
            .iter()
            .find(move |&&rid| set.participants_1(rid).contains(&site))
            .map(|&rid| StereoAtomId::from(rid))
    }

    /// Whether any stereo atom is coincident with `atom` (its site).
    pub fn has_coincident(&self, atom: AtomId) -> bool {
        let set = self.stereo_atoms;
        let site = NodeId::from(atom);
        set.incident(site)
            .iter()
            .any(move |&rid| set.participants_1(rid).contains(&site))
    }

    /// View of the stereo atom coincident with `atom` (its site).
    pub fn coincident(&self, atom: AtomId) -> Option<StereoAtomView<'a>> {
        let molecule = self.molecule;
        let set = self.stereo_atoms;
        self.coincident_id(atom).map(move |id| {
            let rid = RelationId::from(id);
            StereoAtomView {
                id,
                site: set.participants_1(rid)[0],
                ligands: set.participants_2(rid),
                ast: set.data(rid),
                molecule,
            }
        })
    }
}

impl<'a> Index<StereoAtomId> for StereoAtomViews<'a> {
    type Output = StereoAtomAst;
    fn index(&self, id: StereoAtomId) -> &StereoAtomAst {
        self.stereo_atoms.data(RelationId::from(id))
    }
}

/// Borrowed view of a stereo atom: the site atom, its ordered ligands, and data.
#[derive(Clone, Copy, Debug)]
pub struct StereoAtomView<'a> {
    pub id: StereoAtomId,
    site: NodeId,
    ligands: &'a [StereoLigand],
    pub ast: &'a StereoAtomAst,
    molecule: &'a MoleculeAst,
}

impl<'a> StereoAtomView<'a> {
    #[inline]
    /// The coordination-geometry kind.
    pub fn kind(&self) -> StereoKind {
        self.ast
            .configuration
            .kind()
            .expect("stereo view has a concrete kind")
    }

    #[inline]
    /// The stereo coset.
    pub fn coset(&self) -> &'a StereoCosetAst {
        self.ast
            .configuration
            .coset()
            .expect("stereo view has a concrete coset")
    }

    #[inline]
    /// The stereo atom constraints.
    pub fn constraints(&self) -> &'a StereoAtomConstraintsAst {
        &self.ast.constraints
    }

    /// ID of the stereo site atom.
    pub fn site_id(&self) -> AtomId {
        AtomId::from(self.site)
    }

    /// View of the stereo site atom.
    pub fn site(&self) -> AtomView<'a> {
        self.molecule.atom(self.site_id())
    }

    pub fn ligand_count(&self) -> usize {
        self.ligands.len()
    }

    /// The ordered ligands occupying the site's coordination positions.
    pub fn ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        let molecule = self.molecule;
        let ligands = self.ligands;
        ligands
            .iter()
            .map(move |ligand| StereoLigandView::new(*ligand, molecule))
    }

    /// View of the ligand at the given coordination position. Panics if it is
    /// not a coordination position of this stereo atom.
    pub fn ligand(&self, ligand_id: StereoLigandPosition) -> StereoLigandView<'a> {
        let ligand = *self
            .ligands
            .get(ligand_id.index())
            .expect("ligand id must refer to a ligand of this stereo atom");
        StereoLigandView::new(ligand, self.molecule)
    }

    pub fn atom_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::Atom)
    }

    pub fn atom_ligand_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atom_ligands().map(|ligand| ligand.atom_id())
    }

    pub fn atom_ligand_count(&self) -> usize {
        self.atom_ligands().count()
    }

    pub fn implicit_hydrogen_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::ImplicitHydrogen)
    }

    pub fn implicit_hydrogen_atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.implicit_hydrogen_ligands()
            .map(|ligand| ligand.atom_id())
    }

    pub fn implicit_hydrogen_count(&self) -> usize {
        self.implicit_hydrogen_ligands().count()
    }

    pub fn lone_pair_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::LonePair)
    }

    pub fn lone_pair_atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.lone_pair_ligands().map(|ligand| ligand.atom_id())
    }

    pub fn lone_pair_count(&self) -> usize {
        self.lone_pair_ligands().count()
    }

    pub fn permutation_for(
        &self,
        ligands: impl IntoIterator<Item = StereoLigand>,
    ) -> Option<Permutation> {
        permutation_for_ligands(self.ligands, ligands)
    }

    pub fn coset_for(
        &self,
        ligands: impl IntoIterator<Item = StereoLigand>,
    ) -> Option<StereoCosetAst> {
        let permutation = self.permutation_for(ligands)?;
        coset_apply_permutation(self.coset(), permutation, self.kind())
    }

    /// Site atom followed by the distinct ligand atoms — the relation's atom
    /// incidence. Deduped: a virtual ligand's bearing atom is the site, so it is
    /// not repeated.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let site = self.site_id();
        let ligands = self.ligands;
        let mut seen = HashSet::new();
        iter::once(site)
            .chain(ligands.iter().map(|l| l.atom_id))
            .filter(move |id| seen.insert(*id))
    }

    /// Rings from the molecule's canonical `RingSet` sharing at least one atom
    /// with this stereo atom (site or ligand).
    pub fn overlapping_rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let atoms: Vec<AtomId> = self.atom_ids().collect();
        self.molecule
            .rings()
            .iter()
            .filter(move |r| r.atoms().iter().any(|a| atoms.contains(a)))
    }

    pub fn is_ground(&self) -> bool {
        self.ast.is_ground()
    }
}

/// Namespace accessor for stereo-bond views on a `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct StereoBondViews<'a> {
    molecule: &'a MoleculeAst,
    stereo_bonds: &'a StereoBondSet,
}

impl<'a> StereoBondViews<'a> {
    pub(crate) fn new(molecule: &'a MoleculeAst, stereo_bonds: &'a StereoBondSet) -> Self {
        Self {
            molecule,
            stereo_bonds,
        }
    }

    pub fn count(&self) -> usize {
        self.stereo_bonds.relation_count()
    }

    pub fn ids(&self) -> impl Iterator<Item = StereoBondId> {
        self.stereo_bonds.relation_ids().map(StereoBondId::from)
    }

    /// Id of the stereo bond on `site` with exactly this ligand multiset, if any. Ligands are a
    /// multiset (virtual ligands repeat); their frame order is not matched.
    pub fn connecting_id(&self, site: BondId, ligands: &[StereoLigand]) -> Option<StereoBondId> {
        self.stereo_bonds
            .find_by_participants(&[EdgeId::from(site)], ligands)
            .map(StereoBondId::from)
    }

    pub fn iter(&self) -> impl Iterator<Item = StereoBondView<'a>> {
        let molecule = self.molecule;
        let set = self.stereo_bonds;
        set.relation_ids().map(move |rid| StereoBondView {
            id: StereoBondId::from(rid),
            site: set.participants_1(rid)[0],
            ligands: set.participants_2(rid),
            ast: set.data(rid),
            molecule,
        })
    }

    /// Whether two stereo bonds share a site — the "≤1 stereo center per site" structural conflict; the
    /// stereo-bond twin of [`StereoAtomViews::has_conflict`].
    pub fn has_conflict(&self) -> bool {
        let mut sites: HashSet<BondId> = HashSet::new();
        self.iter().any(|view| !sites.insert(view.site_id()))
    }

    pub fn contains(&self, id: StereoBondId) -> bool {
        self.stereo_bonds.contains(RelationId::from(id))
    }

    pub fn get(&self, id: StereoBondId) -> Option<StereoBondView<'a>> {
        if !self.contains(id) {
            return None;
        }
        let rid = RelationId::from(id);
        Some(StereoBondView {
            id,
            site: self.stereo_bonds.participants_1(rid)[0],
            ligands: self.stereo_bonds.participants_2(rid),
            ast: self.stereo_bonds.data(rid),
            molecule: self.molecule,
        })
    }

    /// Ids of stereo bonds incident on `atom` (site endpoint or ligand). The
    /// site is an edge, so node incidence covers only ligands; site-endpoint
    /// membership is unioned in (and deduped) explicitly.
    pub fn incident_ids(&self, atom: AtomId) -> impl Iterator<Item = StereoBondId> + 'a {
        let ligand_ids = self
            .stereo_bonds
            .incident(NodeId::from(atom))
            .iter()
            .map(|&rid| StereoBondId::from(rid));
        let mut seen = HashSet::new();
        self.incident_at_site_ids(atom)
            .chain(ligand_ids)
            .filter(move |id| seen.insert(*id))
    }

    /// Any stereo bond is incident on `atom` (site endpoint or ligand).
    pub fn has_incident(&self, atom: AtomId) -> bool {
        self.stereo_bonds.has_incident(NodeId::from(atom)) || self.has_incident_at_site(atom)
    }

    /// Views of stereo bonds incident on `atom`.
    pub fn incident(&self, atom: AtomId) -> impl Iterator<Item = StereoBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.stereo_bonds;
        self.incident_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            StereoBondView {
                id,
                site: set.participants_1(rid)[0],
                ligands: set.participants_2(rid),

                ast: set.data(rid),
                molecule,
            }
        })
    }

    /// Ids of stereo bonds incident, in which `atom` is ligand.
    pub fn incident_at_ligand_ids(&self, atom: AtomId) -> impl Iterator<Item = StereoBondId> + 'a {
        let set = self.stereo_bonds;
        let ligand = StereoLigand {
            atom_id: atom,
            kind: StereoLigandKind::Atom,
        };
        set.incident(NodeId::from(atom))
            .iter()
            .filter(move |&&rid| set.participants_2(rid).contains(&ligand))
            .map(|&rid| StereoBondId::from(rid))
    }

    /// Any stereo bond is incident, in which `atom` is ligand.
    pub fn has_incident_at_ligand(&self, atom: AtomId) -> bool {
        let set = self.stereo_bonds;
        let ligand = StereoLigand {
            atom_id: atom,
            kind: StereoLigandKind::Atom,
        };
        set.incident(NodeId::from(atom))
            .iter()
            .any(|&rid| set.participants_2(rid).contains(&ligand))
    }

    /// Views of stereo bonds incident, in which `atom` is ligand.
    pub fn incident_at_ligand(
        &self,
        atom: AtomId,
    ) -> impl Iterator<Item = StereoBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.stereo_bonds;
        self.incident_at_ligand_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            StereoBondView {
                id,
                site: set.participants_1(rid)[0],
                ligands: set.participants_2(rid),
                ast: set.data(rid),
                molecule,
            }
        })
    }

    /// Ids of stereo bonds, in which `atom` is a site endpoint.
    pub fn incident_at_site_ids(&self, atom: AtomId) -> impl Iterator<Item = StereoBondId> + 'a {
        let set = self.stereo_bonds;
        self.molecule
            .neighbors(atom)
            .flat_map(move |n| set.incident_edge(EdgeId::from(n.bond_id())).iter().copied())
            .map(StereoBondId::from)
    }

    /// Any stereo bond, in which `atom` is a site endpoint.
    pub fn has_incident_at_site(&self, atom: AtomId) -> bool {
        let set = self.stereo_bonds;
        self.molecule
            .neighbors(atom)
            .any(move |n| set.has_incident_edge(EdgeId::from(n.bond_id())))
    }

    /// Views of stereo bonds, in which `atom` is a site endpoint.
    pub fn incident_at_site(&self, atom: AtomId) -> impl Iterator<Item = StereoBondView<'a>> + 'a {
        let molecule = self.molecule;
        let set = self.stereo_bonds;
        self.incident_at_site_ids(atom).map(move |id| {
            let rid = RelationId::from(id);
            StereoBondView {
                id,
                site: set.participants_1(rid)[0],
                ligands: set.participants_2(rid),
                ast: set.data(rid),
                molecule,
            }
        })
    }

    /// Id of the stereo bond coincident with `bond` (its site).
    pub fn coincident_id(&self, bond: BondId) -> Option<StereoBondId> {
        self.stereo_bonds
            .incident_edge(EdgeId::from(bond))
            .first()
            .map(|&rid| StereoBondId::from(rid))
    }

    /// Whether any stereo bond is coincident with `bond` (its site).
    pub fn has_coincident(&self, bond: BondId) -> bool {
        self.stereo_bonds.has_incident_edge(EdgeId::from(bond))
    }

    /// View of the stereo bond coincident with `bond` (its site).
    pub fn coincident(&self, bond: BondId) -> Option<StereoBondView<'a>> {
        let molecule = self.molecule;
        let set = self.stereo_bonds;
        self.coincident_id(bond).map(move |id| {
            let rid = RelationId::from(id);
            StereoBondView {
                id,
                site: set.participants_1(rid)[0],
                ligands: set.participants_2(rid),
                ast: set.data(rid),
                molecule,
            }
        })
    }
}

impl<'a> Index<StereoBondId> for StereoBondViews<'a> {
    type Output = StereoBondAst;
    fn index(&self, id: StereoBondId) -> &StereoBondAst {
        self.stereo_bonds.data(RelationId::from(id))
    }
}

/// Borrowed view of a stereo bond: the site bond, its ordered ligands, and data.
#[derive(Clone, Copy, Debug)]
pub struct StereoBondView<'a> {
    pub id: StereoBondId,
    site: EdgeId,
    ligands: &'a [StereoLigand],
    pub ast: &'a StereoBondAst,
    molecule: &'a MoleculeAst,
}

impl<'a> StereoBondView<'a> {
    #[inline]
    /// The coordination-geometry kind.
    pub fn kind(&self) -> StereoKind {
        self.ast
            .configuration
            .kind()
            .expect("stereo view has a concrete kind")
    }

    #[inline]
    /// The stereo coset.
    pub fn coset(&self) -> &'a StereoCosetAst {
        self.ast
            .configuration
            .coset()
            .expect("stereo view has a concrete coset")
    }

    #[inline]
    /// The stereo bond constraints.
    pub fn constraints(&self) -> &'a StereoBondConstraintsAst {
        &self.ast.constraints
    }

    /// ID of the stereo site bond.
    pub fn site_id(&self) -> BondId {
        BondId::from(self.site)
    }

    /// View of the stereo site bond.
    pub fn site(&self) -> BondView<'a> {
        self.molecule.bond(self.site_id())
    }

    pub fn ligand_count(&self) -> usize {
        self.ligands.len()
    }

    /// The ordered ligands defining the bond's configuration.
    pub fn ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        let molecule = self.molecule;
        let ligands = self.ligands;
        ligands
            .iter()
            .map(move |ligand| StereoLigandView::new(*ligand, molecule))
    }

    /// View of the ligand at the given coordination position. Panics if it is
    /// not a coordination position of this stereo bond.
    pub fn ligand(&self, ligand_id: StereoLigandPosition) -> StereoLigandView<'a> {
        let ligand = *self
            .ligands
            .get(ligand_id.index())
            .expect("ligand id must refer to a ligand of this stereo bond");
        StereoLigandView::new(ligand, self.molecule)
    }

    pub fn atom_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::Atom)
    }

    pub fn atom_ligand_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.atom_ligands().map(|ligand| ligand.atom_id())
    }

    pub fn atom_ligand_count(&self) -> usize {
        self.atom_ligands().count()
    }

    pub fn implicit_hydrogen_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::ImplicitHydrogen)
    }

    pub fn implicit_hydrogen_atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.implicit_hydrogen_ligands()
            .map(|ligand| ligand.atom_id())
    }

    pub fn implicit_hydrogen_count(&self) -> usize {
        self.implicit_hydrogen_ligands().count()
    }

    pub fn lone_pair_ligands(&self) -> impl Iterator<Item = StereoLigandView<'a>> + 'a {
        self.ligands()
            .filter(|ligand| ligand.kind() == StereoLigandKind::LonePair)
    }

    pub fn lone_pair_atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        self.lone_pair_ligands().map(|ligand| ligand.atom_id())
    }

    pub fn lone_pair_count(&self) -> usize {
        self.lone_pair_ligands().count()
    }

    pub fn permutation_for(
        &self,
        ligands: impl IntoIterator<Item = StereoLigand>,
    ) -> Option<Permutation> {
        permutation_for_ligands(self.ligands, ligands)
    }

    pub fn coset_for(
        &self,
        ligands: impl IntoIterator<Item = StereoLigand>,
    ) -> Option<StereoCosetAst> {
        let permutation = self.permutation_for(ligands)?;
        coset_apply_permutation(self.coset(), permutation, self.kind())
    }

    /// The site bond's two atoms followed by the distinct ligand atoms — the
    /// relation's atom incidence. Deduped: a virtual ligand's bearing atom is a
    /// site endpoint, so it is not repeated.
    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + 'a {
        let [a, b] = self.site().atom_ids();
        let ligands = self.ligands;
        let mut seen = HashSet::new();
        [a, b]
            .into_iter()
            .chain(ligands.iter().map(|l| l.atom_id))
            .filter(move |id| seen.insert(*id))
    }

    /// Rings from the molecule's canonical `RingSet` sharing at least one atom
    /// with this stereo bond (site endpoints or ligands).
    pub fn overlapping_rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let atoms: Vec<AtomId> = self.atom_ids().collect();
        self.molecule
            .rings()
            .iter()
            .filter(move |r| r.atoms().iter().any(|a| atoms.contains(a)))
    }

    pub fn is_ground(&self) -> bool {
        self.ast.is_ground()
    }
}

/// Stereo query methods shared by `StereoAtomView` and `StereoBondView`. The
/// orbit/stereogenicity queries are pure reads over a per-carrier `StereoSymmetry`
/// the caller has already computed (`MoleculeAst::stereo_atom_symmetry` /
/// `stereo_bond_symmetry`); ops compute it once and never share it across ops.
macro_rules! stereo_view_queries {
    ($view:ident) => {
        impl $view<'_> {
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

            /// The frame position of an atom ligand, if `atom` is one.
            pub fn ligand_position(&self, atom: AtomId) -> Option<StereoLigandPosition> {
                self.ligands()
                    .position(|l| l.kind() == StereoLigandKind::Atom && l.atom_id() == atom)
                    .map(|i| StereoLigandPosition(i as u32))
            }

            /// The ordered ligand frame (atom ligands + virtual implicit-Hs / lone-pairs).
            pub fn ligand_frame(&self) -> Vec<StereoLigand> {
                self.ligands()
                    .map(|l| StereoLigand::new(l.atom_id(), l.kind()))
                    .collect()
            }
        }
    };
}

stereo_view_queries!(StereoAtomView);
stereo_view_queries!(StereoBondView);

fn has_unique_ligands(ligands: &[StereoLigand]) -> bool {
    ligands.iter().copied().collect::<HashSet<_>>().len() == ligands.len()
}

fn permutation_for_ligands(
    current: &[StereoLigand],
    ligands: impl IntoIterator<Item = StereoLigand>,
) -> Option<Permutation> {
    let current: Vec<StereoLigand> = current.to_vec();
    let requested: Vec<StereoLigand> = ligands.into_iter().collect();
    if current.len() != requested.len()
        || !has_unique_ligands(&current)
        || !has_unique_ligands(&requested)
    {
        return None;
    }
    let current_set: HashSet<StereoLigand> = current.iter().copied().collect();
    let requested_set: HashSet<StereoLigand> = requested.iter().copied().collect();
    (current_set == requested_set).then(|| Permutation::between(&current, &requested))
}

// Builder-scope view bundles for stereo elements. `ligands` is a borrow into
// builder storage so old-state checks compare without cloning; callers clone
// only what they keep (the `ast`).

pub struct StereoAtomEditorView<'a> {
    pub id: StereoAtomId,
    pub ast: &'a StereoAtomAst,
    pub site: AtomId,
    pub ligands: &'a [StereoLigand],
}

pub struct StereoBondEditorView<'a> {
    pub id: StereoBondId,
    pub ast: &'a StereoBondAst,
    pub site: BondId,
    pub ligands: &'a [StereoLigand],
}

pub struct StereoAtomEditorViewMut<'a> {
    pub id: StereoAtomId,
    pub ast: &'a mut StereoAtomAst,
    pub site: AtomId,
    pub ligands: &'a [StereoLigand],
}

pub struct StereoBondEditorViewMut<'a> {
    pub id: StereoBondId,
    pub ast: &'a mut StereoBondAst,
    pub site: BondId,
    pub ligands: &'a [StereoLigand],
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_perm::Permutation;

    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::coloring::ConstitutionColoring;
    use crate::ast::constraint::Constraints;
    use crate::ast::id::{AtomId, BondId, StereoAtomId, StereoBondId, StereoLigandPosition};
    use crate::ast::ligand::{StereoLigand, StereoLigandKind};
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::stereo::{
        StereoAtomAst, StereoBondAst, StereoCosetAst, StereoKind, Stereogenicity, Topicity,
    };
    use crate::ast::symmetry::GraphSymmetryConfig;

    #[fixture]
    fn molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            // Atom 6 is an unbonded spare: a node that is neither stereo-bond site nor ligand.
            vec![AtomAst::from_element(Element::C); 7],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(2)),
                (AtomId(4), AtomId(5), BondAst::from_order(1)),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
            )],
            vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
            )],
            Constraints::default(),
        )
    }

    #[fixture]
    fn virtual_ligand_molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 6],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(2)),
                (AtomId(4), AtomId(5), BondAst::from_order(1)),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
            )],
            vec![(
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
            )],
            Constraints::default(),
        )
    }

    // A 4-membered ring (atoms 0-1-2-3) with two pendant atoms (4 on 0, 5 on 1);
    // a stereo atom on ring atom 0 and a stereo bond on ring bond 0-1.
    #[fixture]
    fn ring_molecule() -> MoleculeAst {
        MoleculeAst::from_parts(
            vec![AtomAst::from_element(Element::C); 6],
            vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
                (AtomId(3), AtomId(0), BondAst::from_order(1)),
                (AtomId(0), AtomId(4), BondAst::from_order(1)),
                (AtomId(1), AtomId(5), BondAst::from_order(1)),
            ],
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
            )],
            vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
            )],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_stereo_atom_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_atoms().count(), 1);
    }

    #[rstest]
    fn test_stereo_atom_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.stereo_atoms().ids().collect::<Vec<_>>(),
            vec![StereoAtomId(0)],
        );
    }

    #[rstest]
    #[case::present(StereoAtomId(0), true)]
    #[case::absent(StereoAtomId(99), false)]
    fn test_stereo_atom_views_contains(
        molecule: MoleculeAst,
        #[case] id: StereoAtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.stereo_atoms().contains(id), expected);
    }

    #[rstest]
    fn test_stereo_atom_views_get(molecule: MoleculeAst) {
        let res = molecule.stereo_atoms().get(StereoAtomId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id, StereoAtomId(0));
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
            view.ast,
            &StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
        );
    }

    #[rstest]
    fn test_stereo_atom_views_get_none(molecule: MoleculeAst) {
        let res = molecule.stereo_atoms().get(StereoAtomId(99));
        assert!(res.is_none());
    }

    #[rstest]
    #[case::site(AtomId(0), vec![StereoAtomId(0)])]
    #[case::ligand(AtomId(2), vec![StereoAtomId(0)])]
    #[case::unrelated(AtomId(5), vec![])]
    fn test_stereo_atom_views_incident(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoAtomId>,
    ) {
        assert_eq!(
            molecule
                .stereo_atoms()
                .incident_ids(atom)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::ligand(AtomId(2), vec![StereoAtomId(0)])]
    #[case::site_not_ligand(AtomId(0), vec![])]
    #[case::unrelated(AtomId(5), vec![])]
    fn test_stereo_atom_views_incident_at_ligand(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoAtomId>,
    ) {
        assert_eq!(
            molecule
                .stereo_atoms()
                .incident_at_ligand_ids(atom)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::site(AtomId(0), Some(StereoAtomId(0)))]
    #[case::ligand_not_site(AtomId(2), None)]
    #[case::unrelated(AtomId(5), None)]
    fn test_stereo_atom_views_coincident(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Option<StereoAtomId>,
    ) {
        assert_eq!(molecule.stereo_atoms().coincident_id(atom), expected);
    }

    #[rstest]
    fn test_stereo_atom_view_site_id(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_atom(StereoAtomId(0)).site_id(), AtomId(0));
    }

    #[rstest]
    fn test_stereo_atom_view_site(molecule: MoleculeAst) {
        let view = molecule.stereo_atom(StereoAtomId(0)).site();
        assert_eq!(view.id, AtomId(0));
        assert_eq!(view.ast, &AtomAst::from_element(Element::C));
    }

    #[rstest]
    fn test_stereo_atom_view_coset(molecule: MoleculeAst) {
        assert_eq!(
            molecule.stereo_atom(StereoAtomId(0)).coset(),
            &StereoCosetAst::Lit(1),
        );
    }

    #[rstest]
    fn test_stereo_atom_view_ligand_count(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_atom(StereoAtomId(0)).ligand_count(), 4);
    }

    #[rstest]
    fn test_stereo_atom_view_ligands(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .stereo_atom(StereoAtomId(0))
                .ligands()
                .map(|ligand| (ligand.kind(), ligand.atom_id()))
                .collect::<Vec<_>>(),
            vec![
                (StereoLigandKind::Atom, AtomId(1)),
                (StereoLigandKind::Atom, AtomId(2)),
                (StereoLigandKind::Atom, AtomId(3)),
                (StereoLigandKind::Atom, AtomId(4)),
            ],
        );
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), StereoLigandKind::Atom, AtomId(1))]
    #[case::last(StereoLigandPosition(3), StereoLigandKind::Atom, AtomId(4))]
    fn test_stereo_atom_view_ligand(
        molecule: MoleculeAst,
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
        let mol = MoleculeAst::from_parts(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::F),
                AtomAst::from_element(Element::Cl),
                AtomAst::from_element(Element::Br),
                AtomAst::from_element(Element::I),
            ],
            (1..=4)
                .map(|i| (AtomId(0), AtomId(i), BondAst::from_order(1)))
                .collect(),
            vec![],
            vec![],
            vec![],
            vec![],
            vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                ],
                StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(0)),
            )],
            vec![],
            Constraints::default(),
        );
        let gs = mol.graph_symmetry(&GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: true,
            max_iterations: 16,
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
        assert_eq!(view.ligand_frame().len(), 4);
    }

    #[rstest]
    #[should_panic]
    fn test_stereo_atom_view_ligand_error(molecule: MoleculeAst) {
        molecule
            .stereo_atom(StereoAtomId(0))
            .ligand(StereoLigandPosition(4));
    }

    #[rstest]
    fn test_stereo_ligand_view_atom(molecule: MoleculeAst) {
        let ligand = molecule
            .stereo_atom(StereoAtomId(0))
            .ligands()
            .next()
            .unwrap();
        let atom = ligand.atom();
        assert_eq!(atom.id, AtomId(1));
        assert_eq!(atom.ast, &AtomAst::from_element(Element::C));
    }

    #[rstest]
    fn test_stereo_atom_view_atom_ligands(virtual_ligand_molecule: MoleculeAst) {
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
    fn test_stereo_atom_view_implicit_hydrogen_ligands(virtual_ligand_molecule: MoleculeAst) {
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
    fn test_stereo_atom_view_lone_pair_ligands(virtual_ligand_molecule: MoleculeAst) {
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
    fn test_stereo_atom_view_permutation_for(molecule: MoleculeAst) {
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
        assert_eq!(
            view.permutation_for(ligands.clone()),
            Some(Permutation::identity(4)),
        );

        let reordered = vec![ligands[1], ligands[0], ligands[2], ligands[3]];
        assert_eq!(
            view.permutation_for(reordered),
            Some(Permutation::from_image(4, &[1, 0, 2, 3])),
        );
    }

    #[rstest]
    fn test_stereo_atom_view_permutation_for_none(molecule: MoleculeAst) {
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
        assert_eq!(view.permutation_for(ligands[..3].iter().copied()), None);
        assert_eq!(
            view.permutation_for([ligands[0], ligands[0], ligands[2], ligands[3]]),
            None,
        );
        assert_eq!(
            view.permutation_for([
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
    fn test_stereo_atom_view_coset_for(molecule: MoleculeAst) {
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
        assert_eq!(
            view.coset_for(ligands.clone()),
            Some(StereoCosetAst::Lit(1)),
        );

        let reordered = vec![ligands[1], ligands[0], ligands[2], ligands[3]];
        assert_eq!(view.coset_for(reordered), Some(StereoCosetAst::Lit(0)));
    }

    #[rstest]
    fn test_stereo_atom_view_atom_ids(molecule: MoleculeAst, virtual_ligand_molecule: MoleculeAst) {
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
    fn test_stereo_atom_view_overlapping_rings(ring_molecule: MoleculeAst) {
        let rings: Vec<Vec<AtomId>> = ring_molecule
            .stereo_atom(StereoAtomId(0))
            .overlapping_rings()
            .map(|r| {
                let mut atoms = r.atoms().to_vec();
                atoms.sort_by_key(|a| a.0);
                atoms
            })
            .collect();
        assert_eq!(
            rings,
            vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]],
        );
    }

    #[rstest]
    fn test_stereo_atom_views_index(molecule: MoleculeAst) {
        assert_eq!(
            &molecule.stereo_atoms()[StereoAtomId(0)],
            &StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1)),
        );
    }

    #[rstest]
    fn test_stereo_bond_views_count(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_bonds().count(), 1);
    }

    #[rstest]
    fn test_stereo_bond_views_ids(molecule: MoleculeAst) {
        assert_eq!(
            molecule.stereo_bonds().ids().collect::<Vec<_>>(),
            vec![StereoBondId(0)],
        );
    }

    #[rstest]
    #[case::present(StereoBondId(0), true)]
    #[case::absent(StereoBondId(99), false)]
    fn test_stereo_bond_views_contains(
        molecule: MoleculeAst,
        #[case] id: StereoBondId,
        #[case] expected: bool,
    ) {
        assert_eq!(molecule.stereo_bonds().contains(id), expected);
    }

    #[rstest]
    fn test_stereo_bond_views_get(molecule: MoleculeAst) {
        let res = molecule.stereo_bonds().get(StereoBondId(0));
        assert!(res.is_some());
        let view = res.unwrap();
        assert_eq!(view.id, StereoBondId(0));
        assert_eq!(view.site_id(), BondId(1));
        assert_eq!(view.kind(), StereoKind::CisTrans);
        assert_eq!(
            view.ast,
            &StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
        );
    }

    #[rstest]
    fn test_stereo_bond_views_get_none(molecule: MoleculeAst) {
        let res = molecule.stereo_bonds().get(StereoBondId(99));
        assert!(res.is_none());
    }

    #[rstest]
    #[case::site_endpoint(AtomId(2), vec![StereoBondId(0)])]
    #[case::ligand(AtomId(4), vec![StereoBondId(0)])]
    #[case::unrelated(AtomId(6), vec![])]
    fn test_stereo_bond_views_incident(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoBondId>,
    ) {
        assert_eq!(
            molecule
                .stereo_bonds()
                .incident_ids(atom)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::site_endpoint(AtomId(2), vec![StereoBondId(0)])]
    #[case::ligand_not_site(AtomId(4), vec![])]
    #[case::unrelated(AtomId(6), vec![])]
    fn test_stereo_bond_views_incident_at_site(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoBondId>,
    ) {
        assert_eq!(
            molecule
                .stereo_bonds()
                .incident_at_site_ids(atom)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::ligand(AtomId(4), vec![StereoBondId(0)])]
    #[case::site_not_ligand(AtomId(2), vec![])]
    #[case::unrelated(AtomId(6), vec![])]
    fn test_stereo_bond_views_incident_at_ligand(
        molecule: MoleculeAst,
        #[case] atom: AtomId,
        #[case] expected: Vec<StereoBondId>,
    ) {
        assert_eq!(
            molecule
                .stereo_bonds()
                .incident_at_ligand_ids(atom)
                .collect::<Vec<_>>(),
            expected,
        );
    }

    #[rstest]
    #[case::site(BondId(1), Some(StereoBondId(0)))]
    #[case::non_site(BondId(0), None)]
    fn test_stereo_bond_views_coincident(
        molecule: MoleculeAst,
        #[case] bond: BondId,
        #[case] expected: Option<StereoBondId>,
    ) {
        assert_eq!(molecule.stereo_bonds().coincident_id(bond), expected);
    }

    #[rstest]
    fn test_stereo_bond_view_site_id(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_bond(StereoBondId(0)).site_id(), BondId(1));
    }

    #[rstest]
    fn test_stereo_bond_view_site(molecule: MoleculeAst) {
        let view = molecule.stereo_bond(StereoBondId(0)).site();
        assert_eq!(view.id, BondId(1));
        assert_eq!(view.atom_ids(), [AtomId(2), AtomId(3)]);
    }

    #[rstest]
    fn test_stereo_bond_view_coset(molecule: MoleculeAst) {
        assert_eq!(
            molecule.stereo_bond(StereoBondId(0)).coset(),
            &StereoCosetAst::Lit(1),
        );
    }

    #[rstest]
    fn test_stereo_bond_view_ligand_count(molecule: MoleculeAst) {
        assert_eq!(molecule.stereo_bond(StereoBondId(0)).ligand_count(), 4);
    }

    #[rstest]
    fn test_stereo_bond_view_ligands(molecule: MoleculeAst) {
        assert_eq!(
            molecule
                .stereo_bond(StereoBondId(0))
                .ligands()
                .map(|ligand| (ligand.kind(), ligand.atom_id()))
                .collect::<Vec<_>>(),
            vec![
                (StereoLigandKind::Atom, AtomId(4)),
                (StereoLigandKind::Atom, AtomId(5)),
                (StereoLigandKind::Atom, AtomId(0)),
                (StereoLigandKind::Atom, AtomId(1)),
            ],
        );
    }

    #[rstest]
    #[case::first(StereoLigandPosition(0), StereoLigandKind::Atom, AtomId(4))]
    #[case::second(StereoLigandPosition(1), StereoLigandKind::Atom, AtomId(5))]
    #[case::third(StereoLigandPosition(2), StereoLigandKind::Atom, AtomId(0))]
    #[case::fourth(StereoLigandPosition(3), StereoLigandKind::Atom, AtomId(1))]
    fn test_stereo_bond_view_ligand(
        molecule: MoleculeAst,
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
    fn test_stereo_bond_view_ligand_error(molecule: MoleculeAst) {
        molecule
            .stereo_bond(StereoBondId(0))
            .ligand(StereoLigandPosition(4));
    }

    #[rstest]
    fn test_stereo_bond_view_atom_ligands(virtual_ligand_molecule: MoleculeAst) {
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
    fn test_stereo_bond_view_implicit_hydrogen_ligands(virtual_ligand_molecule: MoleculeAst) {
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
    fn test_stereo_bond_view_lone_pair_ligands(virtual_ligand_molecule: MoleculeAst) {
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
    fn test_stereo_bond_view_permutation_for(molecule: MoleculeAst) {
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
        assert_eq!(
            view.permutation_for(ligands.clone()),
            Some(Permutation::identity(4)),
        );

        let reordered = vec![ligands[1], ligands[0], ligands[2], ligands[3]];
        assert_eq!(
            view.permutation_for(reordered),
            Some(Permutation::from_image(4, &[1, 0, 2, 3])),
        );
    }

    #[rstest]
    fn test_stereo_bond_view_permutation_for_none(molecule: MoleculeAst) {
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
        assert_eq!(view.permutation_for(ligands[..1].iter().copied()), None);
        assert_eq!(
            view.permutation_for([ligands[0], ligands[0], ligands[2], ligands[3]]),
            None,
        );
        assert_eq!(
            view.permutation_for([
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
    fn test_stereo_bond_view_coset_for(molecule: MoleculeAst) {
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
        assert_eq!(
            view.coset_for(ligands.clone()),
            Some(StereoCosetAst::Lit(1)),
        );

        let reordered = vec![ligands[1], ligands[0], ligands[2], ligands[3]];
        assert_eq!(view.coset_for(reordered), Some(StereoCosetAst::Lit(0)));
    }

    #[rstest]
    fn test_stereo_bond_view_atom_ids(molecule: MoleculeAst, virtual_ligand_molecule: MoleculeAst) {
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
    fn test_stereo_bond_view_overlapping_rings(ring_molecule: MoleculeAst) {
        let rings: Vec<Vec<AtomId>> = ring_molecule
            .stereo_bond(StereoBondId(0))
            .overlapping_rings()
            .map(|r| {
                let mut atoms = r.atoms().to_vec();
                atoms.sort_by_key(|a| a.0);
                atoms
            })
            .collect();
        assert_eq!(
            rings,
            vec![vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]],
        );
    }

    #[rstest]
    fn test_stereo_bond_views_index(molecule: MoleculeAst) {
        assert_eq!(
            &molecule.stereo_bonds()[StereoBondId(0)],
            &StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1)),
        );
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
    fn test_stereo_atom_views_connecting_id(
        molecule: MoleculeAst,
        #[case] site: AtomId,
        #[case] ligands: Vec<StereoLigand>,
        #[case] expected: Option<StereoAtomId>,
    ) {
        assert_eq!(
            molecule.stereo_atoms().connecting_id(site, &ligands),
            expected
        );
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
    fn test_stereo_bond_views_connecting_id(
        molecule: MoleculeAst,
        #[case] site: BondId,
        #[case] ligands: Vec<StereoLigand>,
        #[case] expected: Option<StereoBondId>,
    ) {
        assert_eq!(
            molecule.stereo_bonds().connecting_id(site, &ligands),
            expected
        );
    }
}
