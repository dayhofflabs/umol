//! Atom views: `AtomViews` namespace, `AtomView` / `AtomViewMut` AST bundles,
//! `AtomBuilderView` / `AtomBuilderViewMut` builder bundles.

use std::ops::Index;

use umol_shared::element::Element;

use super::super::atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
use super::super::constraint::AtomConstraints;
use super::super::idx::{
    AromaticSystemId, AtomId, BondId, DativeBondId, MulticenterBondId, NoncovalentBondId,
};
use super::super::molecule::MoleculeAst;
use super::super::rings::{RingSet, RingView};
use super::super::spin::SpinStateAst;
use super::super::value::ValueAst;
use super::aromatic_system::AromaticSystemView;
use super::dative_bond::DativeBondView;
use super::multicenter_bond::MulticenterBondView;
use super::neighbor::NeighborView;
use super::noncovalent_bond::NoncovalentBondView;

/// Namespace accessor for atom views on a `MoleculeAst`. Provides `count`,
/// `ids`, `iter`, `get`, and `Index` without burying them on `MoleculeAst`.
#[derive(Clone, Copy)]
pub struct AtomViews<'a> {
    molecule: &'a MoleculeAst,
    atoms: &'a [AtomAst],
}

impl<'a> AtomViews<'a> {
    pub(in crate::ast) fn new(molecule: &'a MoleculeAst, atoms: &'a [AtomAst]) -> Self {
        Self { molecule, atoms }
    }

    pub fn count(&self) -> usize {
        self.atoms.len()
    }

    pub fn ids(&self) -> impl Iterator<Item = AtomId> {
        (0..self.atoms.len() as u32).map(AtomId)
    }

    pub fn iter(&self) -> impl Iterator<Item = AtomView<'a>> {
        let molecule = self.molecule;
        self.atoms.iter().enumerate().map(move |(i, ast)| AtomView {
            id: AtomId(i as u32),
            ast,
            molecule,
        })
    }

    pub fn get(&self, id: AtomId) -> AtomView<'a> {
        AtomView {
            id,
            ast: &self.atoms[id.index()],
            molecule: self.molecule,
        }
    }
}

impl<'a> Index<AtomId> for AtomViews<'a> {
    type Output = AtomAst;
    fn index(&self, id: AtomId) -> &AtomAst {
        &self.atoms[id.index()]
    }
}

/// Borrowed view of an atom: index, underlying `AtomAst`, and the parent
/// `MoleculeAst` for cross-relation chemistry methods.
///
/// Chemistry methods come in pairs: the topology-derived value (summed from
/// incident bonds / dative bonds / aromatic system / multicenter bonds) and
/// the matching local-constraint value carried in `data.constraints`. The
/// validator cross-checks the two when both are ground.
#[derive(Clone, Copy, Debug)]
pub struct AtomView<'a> {
    pub id: AtomId,
    pub ast: &'a AtomAst,
    molecule: &'a MoleculeAst,
}

impl<'a> AtomView<'a> {
    #[inline]
    pub fn element(&self) -> &'a ElementAst {
        &self.ast.element
    }

    #[inline]
    pub fn isotope_mass(&self) -> &'a IsotopeAst {
        &self.ast.isotope_mass
    }

    #[inline]
    pub fn charge(&self) -> &'a ValueAst {
        &self.ast.charge
    }

    #[inline]
    pub fn implicit_hydrogens(&self) -> &'a ImplicitHydrogensAst {
        &self.ast.implicit_hydrogens
    }

    #[inline]
    pub fn lone_pairs(&self) -> &'a ValueAst {
        &self.ast.lone_pairs
    }

    #[inline]
    pub fn spin(&self) -> &'a SpinStateAst {
        &self.ast.spin
    }

    #[inline]
    pub fn constraints(&self) -> &'a AtomConstraints {
        &self.ast.constraints
    }

    /// Iterator over incident bonds and their neighbor atoms. Equivalent to
    /// `self.molecule.neighbors(self.id)` but exposed on the view so closures
    /// that take `&AtomView` (e.g. perception electron-counting) can inspect
    /// bonds without reaching back to the molecule.
    pub fn neighbors(&self) -> impl Iterator<Item = NeighborView<'a>> {
        self.molecule.neighbors(self.id)
    }

    /// IDs of incident localized bonds, in iteration order of `neighbors`.
    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + 'a {
        self.molecule.neighbors(self.id).map(|n| n.bond)
    }

    /// Localized valence: sum of incident `Bond.order` values. Returns
    /// `ValueAst::Lit(n)` when every incident bond order is `Lit`; collapses
    /// to `Undetermined` if any bond order is non-`Lit`.
    pub fn valence(&self) -> ValueAst {
        self.neighbors()
            .map(|n| n.ast.order.clone())
            .fold(ValueAst::Lit(0), |acc, order| acc + order)
    }

    /// Sum of `order` over incident dative bonds where this atom is the sole
    /// donor (multi-donor datives contribute nothing per individual donor —
    /// the donated pair is collective and has no well-defined per-atom
    /// share). Returns `ValueAst::Lit(0)` when this atom donates to no
    /// single-donor dative bonds; collapses to `Undetermined` if any
    /// contributing dative's `order` is non-`Lit`.
    pub fn donated_pairs(&self) -> ValueAst {
        let mut sum = ValueAst::Lit(0);
        for view in self.dative_bonds() {
            let donor_ids: Vec<AtomId> = view.donor_ids().collect();
            if donor_ids.len() != 1 || donor_ids[0] != self.id {
                continue;
            }
            sum = sum + view.ast.order.clone();
        }
        sum
    }

    /// Sum of `order` over incident dative bonds where this atom is the
    /// acceptor. Returns `ValueAst::Lit(0)` when this atom is not an
    /// acceptor; collapses to `Undetermined` if any contributing dative's
    /// `order` is non-`Lit`.
    pub fn accepted_pairs(&self) -> ValueAst {
        let mut sum = ValueAst::Lit(0);
        for view in self.dative_bonds() {
            if view.acceptor_id != self.id {
                continue;
            }
            sum = sum + view.ast.order.clone();
        }
        sum
    }

    /// Electron contribution from the aromatic system this atom belongs to.
    /// `ValueAst::Lit(0)` if the atom is not in any aromatic system;
    /// `Undetermined` if the system's per-atom electron count is non-`Lit`.
    pub fn aromatic_valence(&self) -> ValueAst {
        let Some(sys) = self.aromatic_system() else {
            return ValueAst::Lit(0);
        };
        let Some(pos) = sys.atom_ids().position(|a| a == self.id) else {
            return ValueAst::Undetermined;
        };
        sys.ast
            .electrons
            .get(pos)
            .cloned()
            .unwrap_or(ValueAst::Undetermined)
    }

    pub fn is_in_aromatic_system(&self) -> bool {
        self.aromatic_system().is_some()
    }

    /// The aromatic system containing this atom, if any. Per-perception
    /// design an atom belongs to at most one aromatic system; this
    /// returns the first incident system.
    pub fn aromatic_system(&self) -> Option<AromaticSystemView<'a>> {
        self.aromatic_system_id()
            .map(|id| self.molecule.aromatic_system(id))
    }

    pub fn aromatic_system_id(&self) -> Option<AromaticSystemId> {
        self.molecule
            .aromatic_systems()
            .incident_ids(self.id)
            .next()
    }

    pub fn dative_bonds(&self) -> impl Iterator<Item = DativeBondView<'a>> + 'a {
        self.molecule.dative_bonds().incident(self.id)
    }

    pub fn dative_bond_ids(&self) -> impl Iterator<Item = DativeBondId> + 'a {
        self.molecule.dative_bonds().incident_ids(self.id)
    }

    pub fn multicenter_bonds(&self) -> impl Iterator<Item = MulticenterBondView<'a>> + 'a {
        self.molecule.multicenter_bonds().incident(self.id)
    }

    pub fn multicenter_bond_ids(&self) -> impl Iterator<Item = MulticenterBondId> + 'a {
        self.molecule.multicenter_bonds().incident_ids(self.id)
    }

    pub fn noncovalent_bonds(&self) -> impl Iterator<Item = NoncovalentBondView<'a>> + 'a {
        self.molecule.noncovalent_bonds().incident(self.id)
    }

    pub fn noncovalent_bond_ids(&self) -> impl Iterator<Item = NoncovalentBondId> + 'a {
        self.molecule.noncovalent_bonds().incident_ids(self.id)
    }

    /// True if this atom participates in any of the four overlay relations
    /// (aromatic system, dative bond, multicenter bond, noncovalent bond).
    /// Mirror of `MoleculeAst::has_overlays` scoped to a single atom; useful
    /// as a pre-mutation predicate before structural removal.
    pub fn is_in_overlays(&self) -> bool {
        self.aromatic_system().is_some()
            || self.dative_bonds().next().is_some()
            || self.multicenter_bonds().next().is_some()
            || self.noncovalent_bonds().next().is_some()
    }

    /// True if this atom belongs to any ring in the molecule's canonical
    /// ring set (Vismara relevant cycles, max ring size 22). Uses the
    /// molecule's cached canonical `RingSet`.
    pub fn is_in_ring(&self) -> bool {
        self.molecule.rings().contains_atom(self.id)
    }

    /// True if this atom appears in any ring of the supplied set.
    pub fn is_in_ring_from(&self, rings: &RingSet) -> bool {
        rings.contains_atom(self.id)
    }

    /// Rings containing this atom drawn from the molecule's canonical
    /// `RingSet` (Vismara relevant cycles, max ring size 22).
    pub fn rings(&self) -> impl Iterator<Item = RingView<'a>> + 'a {
        let id = self.id;
        self.molecule
            .rings()
            .iter()
            .filter(move |v| v.atoms().contains(&id))
    }

    /// Rings from the supplied set that contain this atom.
    pub fn rings_from<'r>(&self, rings: &'r RingSet) -> impl Iterator<Item = RingView<'r>> + 'r {
        let id = self.id;
        rings.iter().filter(move |v| v.atoms().contains(&id))
    }

    /// Sum of per-atom contributions across incident multicenter bonds.
    /// `ValueAst::Lit(0)` when not in any multicenter bond; collapses to
    /// `Undetermined` if any contribution is non-`Lit`.
    pub fn multicenter_valence(&self) -> ValueAst {
        let mut sum = ValueAst::Lit(0);
        for view in self.multicenter_bonds() {
            let Some(pos) = view.atom_ids().position(|a| a == self.id) else {
                return ValueAst::Undetermined;
            };
            let term = view
                .ast
                .electrons
                .get(pos)
                .cloned()
                .unwrap_or(ValueAst::Undetermined);
            sum = sum + term;
        }
        sum
    }

    /// Count of incident localized bonds, each weighted 1. Always `Lit`.
    pub fn degree(&self) -> ValueAst {
        ValueAst::Lit(self.neighbors().count() as i64)
    }

    /// `degree` + `implicit_hydrogens` + `multicenter_degree`. Collapses to
    /// `Undetermined` if any term is non-`Lit`.
    pub fn total_degree(&self) -> ValueAst {
        self.degree()
            + ValueAst::from(self.ast.implicit_hydrogens.clone())
            + self.multicenter_degree()
    }

    /// Count of incident localized bonds whose neighbor is not a literal
    /// hydrogen atom (Element::H). Always `Lit`; non-`Lit` neighbor
    /// elements count as heavy (i.e., not filtered out).
    pub fn heavy_atom_degree(&self) -> ValueAst {
        let count = self
            .neighbors()
            .filter(|n| {
                !matches!(
                    self.molecule.atom(n.atom).ast.element,
                    ElementAst::Lit(Element::H),
                )
            })
            .count();
        ValueAst::Lit(count as i64)
    }

    /// `valence` over incident bonds whose neighbor is not a literal
    /// hydrogen. Collapses to `Undetermined` if any contributing bond order
    /// is non-`Lit`.
    pub fn heavy_atom_valence(&self) -> ValueAst {
        self.neighbors()
            .filter(|n| {
                !matches!(
                    self.molecule.atom(n.atom).ast.element,
                    ElementAst::Lit(Element::H),
                )
            })
            .map(|n| n.ast.order.clone())
            .fold(ValueAst::Lit(0), |acc, order| acc + order)
    }

    /// Explicit hydrogens (incident neighbors with `Element::H`) plus
    /// `implicit_hydrogens`. Collapses to `Undetermined` if `implicit_hydrogens`
    /// is non-`Lit` (including `Normal`).
    pub fn total_hydrogens(&self) -> ValueAst {
        let explicit = self
            .neighbors()
            .filter(|n| {
                matches!(
                    self.molecule.atom(n.atom).ast.element,
                    ElementAst::Lit(Element::H),
                )
            })
            .count() as i64;
        ValueAst::Lit(explicit) + ValueAst::from(self.ast.implicit_hydrogens.clone())
    }

    /// Full electron-sharing sum at this atom:
    /// `valence + implicit_hydrogens + aromatic_valence + multicenter_valence`.
    /// Diverges from SMARTS `v<n>` for aromatic lone-pair donors (pyrrole N,
    /// furan O) which contribute the donated pair via `aromatic_valence`.
    pub fn total_valence(&self) -> ValueAst {
        self.valence()
            + ValueAst::from(self.ast.implicit_hydrogens.clone())
            + self.aromatic_valence()
            + self.multicenter_valence()
    }

    /// Count of multicenter co-participants across all incident multicenter
    /// bonds. Per the no-overlap structural rule these are not localized-
    /// bond neighbors. Always `Lit`.
    pub fn multicenter_degree(&self) -> ValueAst {
        let count: usize = self
            .multicenter_bonds()
            .map(|mc| mc.ast.electrons.len().saturating_sub(1))
            .sum();
        ValueAst::Lit(count as i64)
    }

    /// Count of canonical rings (Vismara / max ring size 22) containing
    /// this atom. Always `Lit`.
    pub fn ring_count(&self) -> ValueAst {
        ValueAst::Lit(self.rings().count() as i64)
    }

    /// Sizes of canonical rings containing this atom, in iteration order.
    /// Multi-valued: an atom in fused rings yields one size per ring.
    pub fn ring_size(&self) -> impl Iterator<Item = usize> + 'a {
        self.rings().map(|r| r.len())
    }

    /// Smallest containing canonical ring size, or `None` if not in any
    /// ring. Chemistry-classification helper; not the constraint
    /// counterpart of `RingSize` (which uses interpretation B — "in some
    /// ring of size n").
    pub fn smallest_ring_size(&self) -> Option<usize> {
        self.ring_size().min()
    }

    /// Count of incident bonds that participate in any canonical ring.
    /// Always `Lit`.
    pub fn ring_degree(&self) -> ValueAst {
        let count = self
            .neighbors()
            .filter(|n| self.molecule.bond(n.bond).is_in_ring())
            .count();
        ValueAst::Lit(count as i64)
    }

    /// Sum of bond orders of incident bonds that participate in any
    /// canonical ring. Collapses to `Undetermined` if any contributing
    /// bond's `order` is non-`Lit`.
    pub fn ring_valence(&self) -> ValueAst {
        self.neighbors()
            .filter(|n| self.molecule.bond(n.bond).is_in_ring())
            .map(|n| n.ast.order.clone())
            .fold(ValueAst::Lit(0), |acc, order| acc + order)
    }
}

/// Mutable borrowed view of an atom.
#[derive(Debug)]
pub struct AtomViewMut<'a> {
    pub id: AtomId,
    pub ast: &'a mut AtomAst,
}

// Builder-scope view bundles for atoms.

pub struct AtomBuilderView<'a> {
    pub id: AtomId,
    pub ast: &'a AtomAst,
}

pub struct AtomBuilderViewMut<'a> {
    pub id: AtomId,
    pub ast: &'a mut AtomAst,
}
